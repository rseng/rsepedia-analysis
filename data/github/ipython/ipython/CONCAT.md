# Security Policy

## Reporting a Vulnerability

All IPython and Jupyter security are handled via security@ipython.org. 
You can find more informations on the Jupyter website. https://jupyter.org/security
## Triaging Issues

On the IPython repository,  we strive to trust users and give them responsibility.
By using one of our bots, any user can close issues or add/remove
labels by mentioning the bot and asking it to do things on your behalf.

To close an issue (or PR), even if you did not create it, use the following:

> @meeseeksdev close

This command can be in the middle of another comment, but must start on its
own line. 

To add labels to an issue, ask the bot to `tag` with a comma-separated list of
tags to add:

> @meeseeksdev tag windows, documentation

Only already pre-created tags can be added.  So far, the list is limited to:
`async/await`, `backported`, `help wanted`, `documentation`, `notebook`,
`tab-completion`, `windows`

To remove a label, use the `untag` command:

> @meeseeksdev untag windows, documentation

We'll be adding additional capabilities for the bot and will share them here
when they are ready to be used.

## Opening an Issue

When opening a new Issue, please take the following steps:

1. Search GitHub and/or Google for your issue to avoid duplicate reports.
   Keyword searches for your error messages are most helpful.
2. If possible, try updating to master and reproducing your issue,
   because we may have already fixed it.
3. Try to include a minimal reproducible test case.
4. Include relevant system information.  Start with the output of:

        python -c "import IPython; print(IPython.sys_info())"

   And include any relevant package versions, depending on the issue, such as
   matplotlib, numpy, Qt, Qt bindings (PyQt/PySide), tornado, web browser, etc.

## Pull Requests

Some guidelines on contributing to IPython:

* All work is submitted via Pull Requests.
* Pull Requests can be submitted as soon as there is code worth discussing.
  Pull Requests track the branch, so you can continue to work after the PR is submitted.
  Review and discussion can begin well before the work is complete,
  and the more discussion the better.
  The worst case is that the PR is closed.
* Pull Requests should generally be made against master
* Pull Requests should be tested, if feasible:
    - bugfixes should include regression tests.
    - new behavior should at least get minimal exercise.
* New features and backwards-incompatible changes should be documented by adding
  a new file to the [pr](docs/source/whatsnew/pr) directory, see [the README.md
  there](docs/source/whatsnew/pr/README.md) for details.
* Don't make 'cleanup' pull requests just to change code style.
  We don't follow any style guide strictly, and we consider formatting changes
  unnecessary noise.
  If you're making functional changes, you can clean up the specific pieces of
  code you're working on.

[Travis](http://travis-ci.org/#!/ipython/ipython) does a pretty good job testing
IPython and Pull Requests, but it may make sense to manually perform tests,
particularly for PRs that affect `IPython.parallel` or Windows.

For more detailed information, see our [GitHub Workflow](https://github.com/ipython/ipython/wiki/Dev:-GitHub-workflow).

## Running Tests

All the tests can be run by using
```shell
pytest
```

All the tests for a single module (for example **test_alias**) can be run by using the fully qualified path to the module.
```shell
pytest IPython/core/tests/test_alias.py
```

Only a single test (for example **test_alias_lifecycle**) within a single file can be run by adding the specific test after a `::` at the end:
```shell
pytest IPython/core/tests/test_alias.py::test_alias_lifecycle
```
---
name: Bug report / Question / Feature
about: Anything related to IPython itsel
title: ''
labels: ''
assignees: ''

---

<!-- This is the repository for IPython command line, if you can try to make sure this question/bug/feature belong here and not on one of the Jupyter repositories. 

If it's a generic Python/Jupyter question, try other forums or discourse.jupyter.org.

If you are unsure, it's ok to post here, though, there are few maintainer so you might not get a fast response. 

-->
Documenting What's New
----------------------

When making a new pull request that either adds a new feature, or makes a
backwards-incompatible change to IPython, please add a new `.rst` file in this
directory documenting this change as a part of your Pull Request.

This will allow multiple Pull Requests to do the same without conflicting with
one another. Periodically, IPython developers with commit rights will run a
script and populate [development.rst](../development.rst)
with the contents of this directory, and clean it up.

Files which describe new features can have any name, such as
`antigravity-feature.rst`, whereas backwards incompatible changes **must have**
have a filename starting with `incompat-`, such as
`incompat-switching-to-perl.rst`.  Our "What's new" files always have two
sections,  and this prefix scheme will make sure that the backwards incompatible
changes get routed to their proper section.

To merge these files into :file:`whatsnew/development.rst`, run the script :file:`tools/update_whatsnew.py`.
=============================
 The IPython licensing terms
=============================

IPython is licensed under the terms of the Modified BSD License (also known as
New or Revised or 3-Clause BSD). See the LICENSE file.


About the IPython Development Team
----------------------------------

Fernando Perez began IPython in 2001 based on code from Janko Hauser
<jhauser@zscout.de> and Nathaniel Gray <n8gray@caltech.edu>.  Fernando is still
the project lead.

The IPython Development Team is the set of all contributors to the IPython
project.  This includes all of the IPython subprojects. 

The core team that coordinates development on GitHub can be found here:
https://github.com/ipython/.

Our Copyright Policy
--------------------

IPython uses a shared copyright model. Each contributor maintains copyright
over their contributions to IPython. But, it is important to note that these
contributions are typically only changes to the repositories. Thus, the IPython
source code, in its entirety is not the copyright of any single person or
institution.  Instead, it is the collective copyright of the entire IPython
Development Team.  If individual contributors want to maintain a record of what
changes/contributions they have specific copyright on, they should indicate
their copyright in the commit message of the change, when they commit the
change to one of the IPython repositories.

With this in mind, the following banner should be used in any source code file 
to indicate the copyright and license terms:

::

    # Copyright (c) IPython Development Team.
    # Distributed under the terms of the Modified BSD License.
IPython provides a rich toolkit to help you make the most out of using Python
interactively.  Its main components are:

 * A powerful interactive Python shell
 * A `Jupyter <https://jupyter.org/>`_ kernel to work with Python code in Jupyter
   notebooks and other interactive frontends.

The enhanced interactive Python shells have the following main features:

 * Comprehensive object introspection.

 * Input history, persistent across sessions.

 * Caching of output results during a session with automatically generated
   references.

 * Extensible tab completion, with support by default for completion of python
   variables and keywords, filenames and function keywords.

 * Extensible system of 'magic' commands for controlling the environment and
   performing many tasks related either to IPython or the operating system.

 * A rich configuration system with easy switching between different setups
   (simpler than changing $PYTHONSTARTUP environment variables every time).

 * Session logging and reloading.

 * Extensible syntax processing for special purpose situations.

 * Access to the system shell with user-extensible alias system.

 * Easily embeddable in other Python programs and GUIs.

 * Integrated access to the pdb debugger and the Python profiler.

The latest development version is always available from IPython's `GitHub
site <http://github.com/ipython>`_.
.. image:: https://codecov.io/github/ipython/ipython/coverage.svg?branch=master
    :target: https://codecov.io/github/ipython/ipython?branch=master

.. image:: https://img.shields.io/pypi/v/IPython.svg
    :target: https://pypi.python.org/pypi/ipython

.. image:: https://github.com/ipython/ipython/actions/workflows/test.yml/badge.svg
    :target: https://github.com/ipython/ipython/actions/workflows/test.yml)

.. image:: https://www.codetriage.com/ipython/ipython/badges/users.svg
    :target: https://www.codetriage.com/ipython/ipython/

.. image:: https://raster.shields.io/badge/Follows-NEP29-brightgreen.png
    :target: https://numpy.org/neps/nep-0029-deprecation_policy.html


===========================================
 IPython: Productive Interactive Computing
===========================================

Overview
========

Welcome to IPython.  Our full documentation is available on `ipython.readthedocs.io
<https://ipython.readthedocs.io/en/stable/>`_ and contains information on how to install, use, and
contribute to the project.
IPython (Interactive Python) is a command shell for interactive computing in multiple programming languages, originally developed for the Python programming language, that offers introspection, rich media, shell syntax, tab completion, and history.

**IPython versions and Python Support**

Starting with IPython 7.10, IPython follows `NEP 29 <https://numpy.org/neps/nep-0029-deprecation_policy.html>`_

**IPython 7.17+** requires Python version 3.7 and above.

**IPython 7.10+** requires Python version 3.6 and above.

**IPython 7.0** requires Python version 3.5 and above.

**IPython 6.x** requires Python version 3.3 and above.

**IPython 5.x LTS** is the compatible release for Python 2.7.
If you require Python 2 support, you **must** use IPython 5.x LTS. Please
update your project configurations and requirements as necessary.


The Notebook, Qt console and a number of other pieces are now parts of *Jupyter*.
See the `Jupyter installation docs <https://jupyter.readthedocs.io/en/latest/install.html>`__
if you want to use these.

Main features of IPython
========================
Comprehensive object introspection.

Input history, persistent across sessions.

Caching of output results during a session with automatically generated references.

Extensible tab completion, with support by default for completion of python variables and keywords, filenames and function keywords.

Extensible system of ‘magic’ commands for controlling the environment and performing many tasks related to IPython or the operating system.

A rich configuration system with easy switching between different setups (simpler than changing $PYTHONSTARTUP environment variables every time).

Session logging and reloading.

Extensible syntax processing for special purpose situations.

Access to the system shell with user-extensible alias system.

Easily embeddable in other Python programs and GUIs.

Integrated access to the pdb debugger and the Python profiler.


Development and Instant running
===============================

You can find the latest version of the development documentation on `readthedocs
<https://ipython.readthedocs.io/en/latest/>`_.

You can run IPython from this directory without even installing it system-wide
by typing at the terminal::

   $ python -m IPython

Or see the `development installation docs
<https://ipython.readthedocs.io/en/latest/install/install.html#installing-the-development-version>`_
for the latest revision on read the docs.

Documentation and installation instructions for older version of IPython can be
found on the `IPython website <https://ipython.org/documentation.html>`_



IPython requires Python version 3 or above
==========================================

Starting with version 6.0, IPython does not support Python 2.7, 3.0, 3.1, or
3.2.

For a version compatible with Python 2.7, please install the 5.x LTS Long Term
Support version.

If you are encountering this error message you are likely trying to install or
use IPython from source. You need to checkout the remote 5.x branch. If you are
using git the following should work::

  $ git fetch origin
  $ git checkout 5.x

If you encounter this error message with a regular install of IPython, then you
likely need to update your package manager, for example if you are using `pip`
check the version of pip with::

  $ pip --version

You will need to update pip to the version 9.0.1 or greater. If you are not using
pip, please inquiry with the maintainers of the package for your package
manager.

For more information see one of our blog posts:

    https://blog.jupyter.org/release-of-ipython-5-0-8ce60b8d2e8e

As well as the following Pull-Request for discussion:

    https://github.com/ipython/ipython/pull/9900

This error does also occur if you are invoking ``setup.py`` directly – which you
should not – or are using ``easy_install`` If this is the case, use ``pip
install .`` instead of ``setup.py install`` , and ``pip install -e .`` instead
of ``setup.py develop`` If you are depending on IPython as a dependency you may
also want to have a conditional dependency on IPython depending on the Python
version::

    install_req = ['ipython']
    if sys.version_info[0] < 3 and 'bdist_wheel' not in sys.argv:
        install_req.remove('ipython')
        install_req.append('ipython<6')

    setup(
        ...
        install_requires=install_req
    )

Alternatives to IPython
=======================

IPython may not be to your taste; if that's the case there might be similar
project that you might want to use:

- The classic Python REPL.
- `bpython <https://bpython-interpreter.org/>`_
- `mypython <https://www.asmeurer.com/mypython/>`_
- `ptpython and ptipython <https://pypi.org/project/ptpython/>`_
- `Xonsh <https://xon.sh/>`_

Ignoring commits with git blame.ignoreRevsFile
==============================================

As of git 2.23, it is possible to make formatting changes without breaking
``git blame``. See the `git documentation
<https://git-scm.com/docs/git-config#Documentation/git-config.txt-blameignoreRevsFile>`_
for more details.

To use this feature you must:

- Install git >= 2.23
- Configure your local git repo by running:
   - POSIX: ``tools\configure-git-blame-ignore-revs.sh``
   - Windows:  ``tools\configure-git-blame-ignore-revs.bat``
IPython Documentation
---------------------

This directory contains the majority of the documentation for IPython.


Deploy docs
-----------

Documentation is automatically deployed on ReadTheDocs on every push or merged
Pull requests.


Requirements
------------

The documentation must be built using Python 3.

In addition to :ref:`devinstall`,
the following tools are needed to build the documentation:

 - sphinx
 - sphinx_rtd_theme
 - docrepr

In a conda environment, or a Python 3 ``venv``, you should be able to run::

  cd ipython
  pip install -U -r docs/requirements.txt


Build Commands
--------------

The documentation gets built using ``make``, and comes in several flavors.

``make html`` - build the API and narrative documentation web pages, this is
the default ``make`` target, so running just ``make`` is equivalent to ``make
html``.

``make html_noapi`` - same as above, but without running the auto-generated API
docs. When you are working on the narrative documentation, the most time
consuming portion  of the build process is the processing and rendering of the
API documentation. This build target skips that.

``make pdf`` will compile a pdf from the documentation.

You can run ``make help`` to see information on all possible make targets.

To save time,
the make targets above only process the files that have been changed since the
previous docs build.
To remove the previous docs build you can use ``make clean``.
You can also combine ``clean`` with other `make` commands;
for example,
``make clean html`` will do a complete rebuild of the docs or ``make clean pdf`` will do a complete build of the pdf.


Continuous Integration
----------------------

Documentation builds are included in the Travis-CI continuous integration process,
so you can see the results of the docs build for any pull request at
https://travis-ci.org/ipython/ipython/pull_requests.

.. _ipython_directive:

========================
IPython Sphinx Directive
========================

.. note::

   The IPython Sphinx Directive is in 'beta' and currently under
   active development. Improvements to the code or documentation are welcome!

.. |rst| replace:: reStructured text

The :rst:dir:`ipython` directive is a stateful shell that can be used
in |rst| files.

It knows about standard ipython prompts, and extracts the input and output
lines.  These prompts will be renumbered starting at ``1``.  The inputs will be
fed to an embedded ipython interpreter and the outputs from that interpreter
will be inserted as well.  For example, code blocks like the following::

  .. ipython::

     In [136]: x = 2

     In [137]: x**3
     Out[137]: 8

will be rendered as

.. ipython::

   In [136]: x = 2

   In [137]: x**3
   Out[137]: 8

.. note::

   This tutorial should be read side-by-side with the Sphinx source
   for this document because otherwise you will see only the rendered
   output and not the code that generated it.  Excepting the example
   above, we will not in general be showing the literal ReST in this
   document that generates the rendered output.


Directive and options
=====================

The IPython directive takes a number of options detailed here.

.. rst:directive:: ipython

   Create an IPython directive.

   .. rst:directive:option:: doctest

      Run a doctest on IPython code blocks in rst.

   .. rst:directive:option:: python

      Used to indicate that the relevant code block does not have IPython prompts.

   .. rst:directive:option:: okexcept

      Allow the code block to raise an exception.

   .. rst:directive:option:: okwarning

      Allow the code block to emit an warning.

   .. rst:directive:option:: suppress

      Silence any warnings or expected errors.

   .. rst:directive:option:: verbatim

      A noop that allows for any text to be syntax highlighted as valid IPython code.

   .. rst:directive:option:: savefig: OUTFILE [IMAGE_OPTIONS]

      Save output from matplotlib to *outfile*.

It's important to note that all of these options can be used for the entire
directive block or they can decorate individual lines of code as explained
in :ref:`pseudo-decorators`.


Persisting the Python session across IPython directive blocks
=============================================================

The state from previous sessions is stored, and standard error is
trapped. At doc build time, ipython's output and std err will be
inserted, and prompts will be renumbered. So the prompt below should
be renumbered in the rendered docs, and pick up where the block above
left off.

.. ipython::
  :verbatim:

  In [138]: z = x*3   # x is recalled from previous block

  In [139]: z
  Out[139]: 6

  In [142]: print(z)
  6

  In [141]: q = z[)   # this is a syntax error -- we trap ipy exceptions
  ------------------------------------------------------------
     File "<ipython console>", line 1
       q = z[)   # this is a syntax error -- we trap ipy exceptions
	     ^
  SyntaxError: invalid syntax


Adding documentation tests to your IPython directive
====================================================

The embedded interpreter supports some limited markup.  For example,
you can put comments in your ipython sessions, which are reported
verbatim.  There are some handy "pseudo-decorators" that let you
doctest the output.  The inputs are fed to an embedded ipython
session and the outputs from the ipython session are inserted into
your doc.  If the output in your doc and in the ipython session don't
match on a doctest assertion, an error will occur.


.. ipython::

   In [1]: x = 'hello world'

   # this will raise an error if the ipython output is different
   @doctest
   In [2]: x.upper()
   Out[2]: 'HELLO WORLD'

   # some readline features cannot be supported, so we allow
   # "verbatim" blocks, which are dumped in verbatim except prompts
   # are continuously numbered
   @verbatim
   In [3]: x.st<TAB>
   x.startswith  x.strip

For more information on @doctest decorator, please refer to the end of this page in Pseudo-Decorators section.

Multi-line input
================

Multi-line input is supported.

.. ipython::
   :verbatim:

   In [130]: url = 'http://ichart.finance.yahoo.com/table.csv?s=CROX\
      .....: &d=9&e=22&f=2009&g=d&a=1&br=8&c=2006&ignore=.csv'

   In [131]: print(url.split('&'))
   ['http://ichart.finance.yahoo.com/table.csv?s=CROX', 'd=9', 'e=22',

Testing directive outputs
=========================

The IPython Sphinx Directive makes it possible to test the outputs that you provide with your code. To do this,
decorate the contents in your directive block with one of the options listed
above.

If an IPython doctest decorator is found, it will take these steps when your documentation is built:

1. Run the *input* lines in your IPython directive block against the current Python kernel (remember that the session
persists across IPython directive blocks);

2. Compare the *output* of this with the output text that you've put in the IPython directive block (what comes
after `Out[NN]`);

3. If there is a difference, the directive will raise an error and your documentation build will fail.

You can do doctesting on multi-line output as well.  Just be careful
when using non-deterministic inputs like random numbers in the ipython
directive, because your inputs are run through a live interpreter, so
if you are doctesting random output you will get an error.  Here we
"seed" the random number generator for deterministic output, and we
suppress the seed line so it doesn't show up in the rendered output

.. ipython::

   In [133]: import numpy.random

   @suppress
   In [134]: numpy.random.seed(2358)

   @doctest
   In [135]: numpy.random.rand(10,2)
   Out[135]:
   array([[0.64524308, 0.59943846],
          [0.47102322, 0.8715456 ],
          [0.29370834, 0.74776844],
          [0.99539577, 0.1313423 ],
          [0.16250302, 0.21103583],
          [0.81626524, 0.1312433 ],
          [0.67338089, 0.72302393],
          [0.7566368 , 0.07033696],
          [0.22591016, 0.77731835],
          [0.0072729 , 0.34273127]])

For more information on @supress and @doctest decorators, please refer to the end of this file in
Pseudo-Decorators section.

Another demonstration of multi-line input and output

.. ipython::
   :verbatim:

   In [106]: print(x)
   jdh

   In [109]: for i in range(10):
      .....:     print(i)
      .....:
      .....:
   0
   1
   2
   3
   4
   5
   6
   7
   8
   9


Most of the "pseudo-decorators" can be used an options to ipython
mode.  For example, to setup matplotlib pylab but suppress the output,
you can do.  When using the matplotlib ``use`` directive, it should
occur before any import of pylab.  This will not show up in the
rendered docs, but the commands will be executed in the embedded
interpreter and subsequent line numbers will be incremented to reflect
the inputs::


  .. ipython::
     :suppress:

     In [144]: from matplotlib.pylab import *

     In [145]: ion()

.. ipython::
   :suppress:

   In [144]: from matplotlib.pylab import *

   In [145]: ion()

Likewise, you can set ``:doctest:`` or ``:verbatim:`` to apply these
settings to the entire block.  For example,

.. ipython::
   :verbatim:

   In [9]: cd mpl/examples/
   /home/jdhunter/mpl/examples

   In [10]: pwd
   Out[10]: '/home/jdhunter/mpl/examples'


   In [14]: cd mpl/examples/<TAB>
   mpl/examples/animation/        mpl/examples/misc/
   mpl/examples/api/              mpl/examples/mplot3d/
   mpl/examples/axes_grid/        mpl/examples/pylab_examples/
   mpl/examples/event_handling/   mpl/examples/widgets

   In [14]: cd mpl/examples/widgets/
   /home/msierig/mpl/examples/widgets

   In [15]: !wc *
       2    12    77 README.txt
      40    97   884 buttons.py
      26    90   712 check_buttons.py
      19    52   416 cursor.py
     180   404  4882 menu.py
      16    45   337 multicursor.py
      36   106   916 radio_buttons.py
      48   226  2082 rectangle_selector.py
      43   118  1063 slider_demo.py
      40   124  1088 span_selector.py
     450  1274 12457 total

You can create one or more pyplot plots and insert them with the
``@savefig`` decorator.

For more information on @savefig decorator, please refer to the end of this page in Pseudo-Decorators section.

.. ipython::

   @savefig plot_simple.png width=4in
   In [151]: plot([1,2,3]);

   # use a semicolon to suppress the output
   @savefig hist_simple.png width=4in
   In [151]: hist(np.random.randn(10000), 100);

In a subsequent session, we can update the current figure with some
text, and then resave

.. ipython::


   In [151]: ylabel('number')

   In [152]: title('normal distribution')

   @savefig hist_with_text.png width=4in
   In [153]: grid(True)

You can also have function definitions included in the source.

.. ipython::

   In [3]: def square(x):
      ...:     """
      ...:     An overcomplicated square function as an example.
      ...:     """
      ...:     if x < 0:
      ...:         x = abs(x)
      ...:     y = x * x
      ...:     return y
      ...:

Then call it from a subsequent section.

.. ipython::

   In [4]: square(3)
   Out [4]: 9

   In [5]: square(-2)
   Out [5]: 4


Writing Pure Python Code
------------------------

Pure python code is supported by the optional argument `python`. In this pure
python syntax you do not include the output from the python interpreter. The
following markup::

   .. ipython:: python

      foo = 'bar'
      print(foo)
      foo = 2
      foo**2

Renders as

.. ipython:: python

   foo = 'bar'
   print(foo)
   foo = 2
   foo**2

We can even plot from python, using the savefig decorator, as well as, suppress
output with a semicolon

.. ipython:: python

   @savefig plot_simple_python.png width=4in
   plot([1,2,3]);

For more information on @savefig decorator, please refer to the end of this page in Pseudo-Decorators section.

Similarly, std err is inserted

.. ipython:: python
   :okexcept:

   foo = 'bar'
   foo[)

Handling Comments
==================

Comments are handled and state is preserved

.. ipython:: python

   # comments are handled
   print(foo)

If you don't see the next code block then the options work.

.. ipython:: python
   :suppress:

   ioff()
   ion()

Splitting Python statements across lines
========================================

Multi-line input is handled.

.. ipython:: python

   line = 'Multi\
           line &\
           support &\
           works'
   print(line.split('&'))

Functions definitions are correctly parsed

.. ipython:: python

   def square(x):
       """
       An overcomplicated square function as an example.
       """
       if x < 0:
           x = abs(x)
       y = x * x
       return y

And persist across sessions

.. ipython:: python

   print(square(3))
   print(square(-2))

Pretty much anything you can do with the ipython code, you can do with
a simple python script. Obviously, though it doesn't make sense
to use the doctest option.

.. _pseudo-decorators:

Pseudo-Decorators
=================

Here are the supported decorators, and any optional arguments they
take.  Some of the decorators can be used as options to the entire
block (eg ``verbatim`` and ``suppress``), and some only apply to the
line just below them (eg ``savefig``).

@suppress

    execute the ipython input block, but suppress the input and output
    block from the rendered output.  Also, can be applied to the entire
    ``.. ipython`` block as a directive option with ``:suppress:``.

@verbatim

    insert the input and output block in verbatim, but auto-increment
    the line numbers. Internally, the interpreter will be fed an empty
    string, so it is a no-op that keeps line numbering consistent.
    Also, can be applied to the entire ``.. ipython`` block as a
    directive option with ``:verbatim:``.

@savefig OUTFILE [IMAGE_OPTIONS]

    save the figure to the static directory and insert it into the
    document, possibly binding it into a minipage and/or putting
    code/figure label/references to associate the code and the
    figure. Takes args to pass to the image directive (*scale*,
    *width*, etc can be kwargs); see `image options
    <http://docutils.sourceforge.net/docs/ref/rst/directives.html#image>`_
    for details.

@doctest

    Compare the pasted in output in the ipython block with the output
    generated at doc build time, and raise errors if they don't
    match. Also, can be applied to the entire ``.. ipython`` block as a
    directive option with ``:doctest:``.

Configuration Options
=====================

ipython_savefig_dir

    The directory in which to save the figures. This is relative to the
    Sphinx source directory. The default is `html_static_path`.

ipython_rgxin

    The compiled regular expression to denote the start of IPython input
    lines. The default is `re.compile('In \[(\d+)\]:\s?(.*)\s*')`. You
    shouldn't need to change this.

ipython_rgxout

    The compiled regular expression to denote the start of IPython output
    lines. The default is `re.compile('Out\[(\d+)\]:\s?(.*)\s*')`. You
    shouldn't need to change this.


ipython_promptin

    The string to represent the IPython input prompt in the generated ReST.
    The default is `'In [%d]:'`. This expects that the line numbers are used
    in the prompt.

ipython_promptout

    The string to represent the IPython prompt in the generated ReST. The
    default is `'Out [%d]:'`. This expects that the line numbers are used
    in the prompt.


Automatically generated documentation
=====================================

.. automodule:: IPython.sphinxext.ipython_directive

.. _overview:

========
Overview
========

One of Python's most useful features is its interactive interpreter.
It allows for very fast testing of ideas without the overhead of
creating test files as is typical in most programming languages.
However, the interpreter supplied with the standard Python distribution
is somewhat limited for extended interactive use.

The goal of IPython is to create a comprehensive environment for
interactive and exploratory computing.  To support this goal, IPython
has three main components:

* An enhanced interactive Python shell.

* A decoupled :ref:`two-process communication model <ipythonzmq>`, which
  allows for multiple clients to connect to a computation kernel, most notably
  the web-based notebook provided with `Jupyter <https://jupyter.org>`_.

* An architecture for interactive parallel computing now part of the
  `ipyparallel` package.

All of IPython is open source (released under the revised BSD license).

Enhanced interactive Python shell
=================================

IPython's interactive shell (:command:`ipython`), has the following goals,
amongst others:

1. Provide an interactive shell superior to Python's default. IPython
   has many features for tab-completion, object introspection, system shell
   access, command history retrieval across sessions, and its own special
   command system for adding functionality when working interactively. It
   tries to be a very efficient environment both for Python code development
   and for exploration of problems using Python objects (in situations like
   data analysis).
  
2. Serve as an embeddable, ready to use interpreter for your own
   programs. An interactive IPython shell can be started with a single call
   from inside another program, providing access to the current namespace.
   This can be very useful both for debugging purposes and for situations
   where a blend of batch-processing and interactive exploration are needed.
  
3. Offer a flexible framework which can be used as the base
   environment for working with other systems, with Python as the underlying
   bridge language. Specifically scientific environments like Mathematica,
   IDL and Matlab inspired its design, but similar ideas can be
   useful in many fields.
  
4. Allow interactive testing of threaded graphical toolkits. IPython
   has support for interactive, non-blocking control of GTK, Qt, WX, GLUT, and
   OS X applications via special threading flags. The normal Python
   shell can only do this for Tkinter applications.

Main features of the interactive shell
--------------------------------------

* Dynamic object introspection. One can access docstrings, function
  definition prototypes, source code, source files and other details
  of any object accessible to the interpreter with a single
  keystroke (:samp:`?`, and using :samp:`??` provides additional detail).
  
* Searching through modules and namespaces with :samp:`*` wildcards, both
  when using the :samp:`?` system and via the :samp:`%psearch` command.

* Completion in the local namespace, by typing :kbd:`TAB` at the prompt.
  This works for keywords, modules, methods, variables and files in the
  current directory. This is supported via the ``prompt_toolkit`` library.
  Custom completers can be implemented easily for different purposes
  (system commands, magic arguments etc.)

* Numbered input/output prompts with command history (persistent
  across sessions and tied to each profile), full searching in this
  history and caching of all input and output.

* User-extensible 'magic' commands. A set of commands prefixed with
  :samp:`%`  or :samp:`%%` is available for controlling IPython itself and provides
  directory control, namespace information and many aliases to
  common system shell commands.

* Alias facility for defining your own system aliases.

* Complete system shell access. Lines starting with :samp:`!` are passed
  directly to the system shell, and using :samp:`!!` or :samp:`var = !cmd` 
  captures shell output into python variables for further use.

* The ability to expand python variables when calling the system shell. In a
  shell command, any python variable prefixed with :samp:`$` is expanded. A
  double :samp:`$$` allows passing a literal :samp:`$` to the shell (for access
  to shell and environment variables like :envvar:`PATH`).

* Filesystem navigation, via a magic :samp:`%cd` command, along with a
  persistent bookmark system (using :samp:`%bookmark`) for fast access to
  frequently visited directories.

* A lightweight persistence framework via the :samp:`%store` command, which
  allows you to save arbitrary Python variables. These get restored
  when you run the :samp:`%store -r` command.

* Automatic indentation and highlighting of code as you type (through the
  `prompt_toolkit` library).

* Macro system for quickly re-executing multiple lines of previous
  input with a single name via the :samp:`%macro` command. Macros can be
  stored persistently via :samp:`%store` and edited via :samp:`%edit`.

* Session logging (you can then later use these logs as code in your
  programs). Logs can optionally timestamp all input, and also store
  session output (marked as comments, so the log remains valid
  Python source code).

* Session restoring: logs can be replayed to restore a previous
  session to the state where you left it.

* Verbose and colored exception traceback printouts. Easier to parse
  visually, and in verbose mode they produce a lot of useful
  debugging information (basically a terminal version of the cgitb
  module).

* Auto-parentheses via the :samp:`%autocall` command: callable objects can be
  executed without parentheses: :samp:`sin 3` is automatically converted to
  :samp:`sin(3)`

* Auto-quoting: using :samp:`,`, or :samp:`;` as the first character forces
  auto-quoting of the rest of the line: :samp:`,my_function a b` becomes
  automatically :samp:`my_function("a","b")`, while :samp:`;my_function a b`
  becomes :samp:`my_function("a b")`.

* Extensible input syntax. You can define filters that pre-process
  user input to simplify input in special situations. This allows
  for example pasting multi-line code fragments which start with
  :samp:`>>>` or :samp:`...` such as those from other python sessions or the
  standard Python documentation.

* Flexible :ref:`configuration system <config_overview>`. It uses a
  configuration file which allows permanent setting of all command-line
  options, module loading, code and file execution. The system allows
  recursive file inclusion, so you can have a base file with defaults and
  layers which load other customizations for particular projects.

* Embeddable. You can call IPython as a python shell inside your own
  python programs. This can be used both for debugging code or for
  providing interactive abilities to your programs with knowledge
  about the local namespaces (very useful in debugging and data
  analysis situations).

* Easy debugger access. You can set IPython to call up an enhanced version of
  the Python debugger (pdb) every time there is an uncaught exception. This
  drops you inside the code which triggered the exception with all the data
  live and it is possible to navigate the stack to rapidly isolate the source
  of a bug. The :samp:`%run` magic command (with the :samp:`-d` option) can run
  any script under pdb's control, automatically setting initial breakpoints for
  you.  This version of pdb has IPython-specific improvements, including
  tab-completion and traceback coloring support. For even easier debugger
  access, try :samp:`%debug` after seeing an exception.

* Profiler support. You can run single statements (similar to
  :samp:`profile.run()`) or complete programs under the profiler's control.
  While this is possible with standard cProfile or profile modules,
  IPython wraps this functionality with magic commands (see :samp:`%prun`
  and :samp:`%run -p`) convenient for rapid interactive work.

* Simple timing information. You can use the :samp:`%timeit` command to get
  the execution time of a Python statement or expression. This machinery is
  intelligent enough to do more repetitions for commands that finish very
  quickly in order to get a better estimate of their running time. 

.. sourcecode:: ipython

    In [1]: %timeit 1+1
    10000000 loops, best of 3: 25.5 ns per loop

    In [2]: %timeit [math.sin(x) for x in range(5000)]
    1000 loops, best of 3: 719 µs per loop

.. 

  To get the timing information for more than one expression, use the
  :samp:`%%timeit` cell magic command.
  

* Doctest support. The special :samp:`%doctest_mode` command toggles a mode
  to use doctest-compatible prompts, so you can use IPython sessions as
  doctest code. By default, IPython also allows you to paste existing
  doctests, and strips out the leading :samp:`>>>` and :samp:`...` prompts in
  them.

.. _ipythonzmq:

Decoupled two-process model
==============================

IPython has abstracted and extended the notion of a traditional
*Read-Evaluate-Print Loop* (REPL) environment by decoupling the *evaluation*
into its own process. We call this process a **kernel**: it receives execution
instructions from clients and communicates the results back to them.

This decoupling allows us to have several clients connected to the same
kernel, and even allows clients and kernels to live on different machines.
With the exclusion of the traditional single process terminal-based IPython
(what you start if you run ``ipython`` without any subcommands), all
other IPython machinery uses this two-process model. Most of this is now part
of the `Jupyter` project, which includes ``jupyter console``,  ``jupyter
qtconsole``, and ``jupyter notebook``.

As an example, this means that when you start ``jupyter qtconsole``, you're
really starting two processes, a kernel and a Qt-based client which can send
commands to and receive results from that kernel. If there is already a kernel
running that you want to connect to, you can pass the  ``--existing`` flag
which will skip initiating a new kernel and connect to the most recent kernel,
instead. To connect to a specific kernel once you have several kernels
running, use the ``%connect_info`` magic to get the unique connection file,
which will be something like ``--existing kernel-19732.json`` but with
different numbers which correspond to the Process ID of the kernel.

You can read more about using `jupyter qtconsole
<https://jupyter.org/qtconsole/>`_, and
`jupyter notebook <https://jupyter-notebook.readthedocs.io/en/latest/>`_. There
is also a :ref:`message spec <messaging>` which documents the protocol for
communication between kernels
and clients.

.. seealso::

    `Frontend/Kernel Model`_ example notebook


Interactive parallel computing
==============================


This functionality is optional and now part of the `ipyparallel
<https://ipyparallel.readthedocs.io/>`_ project.

Portability and Python requirements
-----------------------------------

Version 7.0+ supports Python 3.4 and higher.
Versions 6.x support Python 3.3 and higher.
Versions 2.0 to 5.x work with Python 2.7.x releases and Python 3.3 and higher.
Version 1.0 additionally worked with Python 2.6 and 3.2.
Version 0.12 was the first version to fully support Python 3.

IPython is known to work on the following operating systems:

	* Linux
	* Most other Unix-like OSs (AIX, Solaris, BSD, etc.)
	* Mac OS X
	* Windows (CygWin, XP, Vista, etc.)

See :ref:`here <install_index>` for instructions on how to install IPython.

.. include:: links.txt
.. _introduction:

=====================
IPython Documentation
=====================

.. only:: html

   :Release: |release|
   :Date: |today|

Welcome to the official IPython documentation.

IPython provides a rich toolkit to help you make the most of using Python
interactively.  Its main components are:

* A powerful interactive Python shell.


.. image:: ./_images/ipython-6-screenshot.png
    :alt: Screenshot of IPython 6.0
    :align: center


* A `Jupyter <https://jupyter.org/>`_ kernel to work with Python code in Jupyter
  notebooks and other interactive frontends.

The enhanced interactive Python shells and kernel have the following main
features:

* Comprehensive object introspection.

* Input history, persistent across sessions.

* Caching of output results during a session with automatically generated
  references.

* Extensible tab completion, with support by default for completion of python
  variables and keywords, filenames and function keywords.

* Extensible system of 'magic' commands for controlling the environment and
  performing many tasks related to IPython or the operating system.

* A rich configuration system with easy switching between different setups
  (simpler than changing ``$PYTHONSTARTUP`` environment variables every time).

* Session logging and reloading.

* Extensible syntax processing for special purpose situations.

* Access to the system shell with user-extensible alias system.

* Easily embeddable in other Python programs and GUIs.

* Integrated access to the pdb debugger and the Python profiler.


The Command line interface inherits the above functionality and adds 
 
* real multi-line editing thanks to `prompt_toolkit <https://python-prompt-toolkit.readthedocs.io/en/stable/>`_.
 
* syntax highlighting as you type.

* integration with command line editor for a better workflow.

The kernel also has its share of features. When used with a compatible frontend,
it allows:

* the object to create a rich display of Html, Images, Latex, Sound and
  Video.

* interactive widgets with the use of the `ipywidgets <https://ipywidgets.readthedocs.io/en/stable/>`_ package.


This documentation will walk you through most of the features of the IPython
command line and kernel, as well as describe the internal mechanisms in order
to improve your Python workflow.

You can find the table of content for this documentation in the left
sidebar, allowing you to come back to previous sections or skip ahead, if needed. 


The latest development version is always available from IPython's `GitHub
repository <http://github.com/ipython/ipython>`_.


.. toctree::
   :maxdepth: 1
   :hidden:

   self
   overview
   whatsnew/index
   install/index
   interactive/index
   config/index
   development/index
   coredev/index
   api/index
   sphinxext
   about/index

.. seealso::

   `Jupyter documentation <https://jupyter.readthedocs.io/en/latest/>`__
     The Jupyter documentation provides information about the Notebook code and other Jupyter sub-projects.
   `ipyparallel documentation <https://ipyparallel.readthedocs.io/en/latest/>`__
     Formerly ``IPython.parallel``.


.. only:: html

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`

:orphan:

.. _parallel_index:

====================================
Using IPython for parallel computing
====================================

IPython.parallel has moved to `ipyparallel <https://github.com/ipython/ipyparallel>`_.

.. _shell_mimerenderer:


Mime Renderer Extensions
========================

Like it's cousins, Jupyter Notebooks and JupyterLab, Terminal IPython can be
thought to render a number of mimetypes in the shell. This can be used to either
display inline images if your terminal emulator supports it; or open some
display results with external file viewers.

Registering new mimetype handlers can so far only be done my extensions and
requires 4 steps:

   - Define a callable that takes 2 parameters:``data`` and ``metadata``; return
     value of the callable is so far ignored. This callable is responsible for
     "displaying" the given mimetype. Which can be sending the right escape
     sequences and bytes to the current terminal; or open an external program. -
   - Appending the right mimetype to ``ipython.display_formatter.active_types``
     for IPython to know it should not ignore those mimetypes.
   - Enabling the given mimetype: ``ipython.display_formatter.formatters[mime].enabled = True``
   - Registering above callable with mimetype handler:
     ``ipython.mime_renderers[mime] = handler``


Here is a complete IPython extension to display images inline and convert math
to png, before displaying it inline for iterm2 on macOS ::


    from base64 import encodebytes
    from IPython.lib.latextools import latex_to_png


    def mathcat(data, meta):
        png = latex_to_png(f'$${data}$$'.replace('\displaystyle', '').replace('$$$', '$$'))
        imcat(png, meta)

    IMAGE_CODE = '\033]1337;File=name=name;inline=true;:{}\a'

    def imcat(image_data, metadata):
        try:
            print(IMAGE_CODE.format(encodebytes(image_data).decode()))
        # bug workaround
        except:
            print(IMAGE_CODE.format(image_data))

    def register_mimerenderer(ipython, mime, handler):
        ipython.display_formatter.active_types.append(mime)
        ipython.display_formatter.formatters[mime].enabled = True
        ipython.mime_renderers[mime] = handler

    def load_ipython_extension(ipython):
        register_mimerenderer(ipython, 'image/png', imcat)
        register_mimerenderer(ipython, 'image/jpeg', imcat)
        register_mimerenderer(ipython, 'text/latex', mathcat)

This example only work for iterm2 on macOS and skip error handling for brevity.
One could also invoke an external viewer with ``subprocess.run()`` and a
temporary file, which is left as an exercise.
================================
Integrating with GUI event loops
================================

When the user types ``%gui qt``, IPython integrates itself with the Qt event
loop, so you can use both a GUI and an interactive prompt together. IPython
supports a number of common GUI toolkits, but from IPython 3.0, it is possible
to integrate other event loops without modifying IPython itself.

Supported event loops include ``qt4``, ``qt5``, ``gtk2``, ``gtk3``, ``gtk4``,
``wx``, ``osx`` and ``tk``. Make sure the event loop you specify matches the
GUI toolkit used by your own code.

To make IPython GUI event loop integration occur automatically at every
startup, set the ``c.InteractiveShellApp.gui`` configuration key in your
IPython profile (see :ref:`setting_config`).

If the event loop you use is supported by IPython, turning on event loop
integration follows the steps just described whether you use Terminal IPython
or an IPython kernel.

However, the way Terminal IPython handles event loops is very different from
the way IPython kernel does, so if you need to integrate with a new kind of
event loop, different steps are needed to integrate with each.

Integrating with a new event loop in the terminal
-------------------------------------------------

.. versionchanged:: 5.0

   There is a new API for event loop integration using prompt_toolkit.

In the terminal, IPython uses prompt_toolkit to prompt the user for input.
prompt_toolkit provides hooks to integrate with an external event loop.

To integrate an event loop, define a function which runs the GUI event loop
until there is input waiting for prompt_toolkit to process. There are two ways
to detect this condition::

    # Polling for input.
    def inputhook(context):
        while not context.input_is_ready():
            # Replace this with the appropriate call for the event loop:
            iterate_loop_once()

    # Using a file descriptor to notify the event loop to stop.
    def inputhook2(context):
        fd = context.fileno()
        # Replace the functions below with those for the event loop.
        add_file_reader(fd, callback=stop_the_loop)
        run_the_loop()

Once you have defined this function, register it with IPython:

.. currentmodule:: IPython.terminal.pt_inputhooks

.. function:: register(name, inputhook)

   Register the function *inputhook* as the event loop integration for the
   GUI *name*. If ``name='foo'``, then the user can enable this integration
   by running ``%gui foo``.


Integrating with a new event loop in the kernel
-----------------------------------------------

The kernel runs its own event loop, so it's simpler to integrate with others.
IPython allows the other event loop to take control, but it must call
:meth:`IPython.kernel.zmq.kernelbase.Kernel.do_one_iteration` periodically.

To integrate with this, write a function that takes a single argument,
the IPython kernel instance, arranges for your event loop to call
``kernel.do_one_iteration()`` at least every ``kernel._poll_interval`` seconds,
and starts the event loop.

Decorate this function with :func:`IPython.kernel.zmq.eventloops.register_integration`,
passing in the names you wish to register it for. Here is a slightly simplified
version of the Tkinter integration already included in IPython::

    @register_integration('tk')
    def loop_tk(kernel):
        """Start a kernel with the Tk event loop."""
        from tkinter import Tk

        # Tk uses milliseconds
        poll_interval = int(1000*kernel._poll_interval)
        # For Tkinter, we create a Tk object and call its withdraw method.
        class Timer(object):
            def __init__(self, func):
                self.app = Tk()
                self.app.withdraw()
                self.func = func

            def on_timer(self):
                self.func()
                self.app.after(poll_interval, self.on_timer)

            def start(self):
                self.on_timer()  # Call it once to get things going.
                self.app.mainloop()

        kernel.timer = Timer(kernel.do_one_iteration)
        kernel.timer.start()

Some event loops can go one better, and integrate checking for messages on the
kernel's ZMQ sockets, making the kernel more responsive than plain polling. How
to do this is outside the scope of this document; if you are interested, look at
the integration with Qt in :mod:`IPython.kernel.zmq.eventloops`.
.. _events:
.. _callbacks:

==============
IPython Events
==============

Extension code can register callbacks functions which will be called on specific
events within the IPython code. You can see the current list of available
callbacks, and the parameters that will be passed with each, in the callback
prototype functions defined in :mod:`IPython.core.events`.

To register callbacks, use :meth:`IPython.core.events.EventManager.register`.
For example::

    class VarWatcher(object):
        def __init__(self, ip):
            self.shell = ip
            self.last_x = None
        
        def pre_execute(self):
            self.last_x = self.shell.user_ns.get('x', None)
        
        def pre_run_cell(self, info):
            print('Cell code: "%s"' % info.raw_cell)
        
        def post_execute(self):
            if self.shell.user_ns.get('x', None) != self.last_x:
                print("x changed!")
        
        def post_run_cell(self, result):
            print('Cell code: "%s"' % result.info.raw_cell)
            if result.error_before_exec:
                print('Error before execution: %s' % result.error_before_exec)
        
    def load_ipython_extension(ip):
        vw = VarWatcher(ip)
        ip.events.register('pre_execute', vw.pre_execute)
        ip.events.register('pre_run_cell', vw.pre_run_cell)
        ip.events.register('post_execute', vw.post_execute)
        ip.events.register('post_run_cell', vw.post_run_cell)


Events
======

These are the events IPython will emit. Callbacks will be passed no arguments, unless otherwise specified.

shell_initialized
-----------------

.. code-block:: python

    def shell_initialized(ipython):
        ...

This event is triggered only once, at the end of setting up IPython.
Extensions registered to load by default as part of configuration can use this to execute code to finalize setup.
Callbacks will be passed the InteractiveShell instance.

pre_run_cell
------------

``pre_run_cell`` fires prior to interactive execution (e.g. a cell in a notebook).
It can be used to note the state prior to execution, and keep track of changes.
An object containing information used for the code execution is provided as an argument.

pre_execute
-----------

``pre_execute`` is like ``pre_run_cell``, but is triggered prior to *any* execution.
Sometimes code can be executed by libraries, etc. which
skipping the history/display mechanisms, in which cases ``pre_run_cell`` will not fire.

post_run_cell
-------------

``post_run_cell`` runs after interactive execution (e.g. a cell in a notebook).
It can be used to cleanup or notify or perform operations on any side effects produced during execution.
For instance, the inline matplotlib backend uses this event to display any figures created but not explicitly displayed during the course of the cell.
The object which will be returned as the execution result is provided as an
argument.

post_execute
------------

The same as ``pre_execute``, ``post_execute`` is like ``post_run_cell``,
but fires for *all* executions, not just interactive ones.


.. seealso::

   Module :mod:`IPython.core.hooks`
     The older 'hooks' system allows end users to customise some parts of
     IPython's behaviour.
   
   :doc:`inputtransforms`
     By registering input transformers that don't change code, you can monitor
     what is being executed.
.. _integrating:

=====================================
Integrating your objects with IPython
=====================================

Tab completion
==============

To change the attributes displayed by tab-completing your object, define a
``__dir__(self)`` method for it. For more details, see the documentation of the
built-in `dir() function <http://docs.python.org/library/functions.html#dir>`_.

You can also customise key completions for your objects, e.g. pressing tab after
``obj["a``. To do so, define a method ``_ipython_key_completions_()``, which
returns a list of objects which are possible keys in a subscript expression
``obj[key]``.

.. versionadded:: 5.0
   Custom key completions

.. _integrating_rich_display:

Rich display
============

Custom methods
----------------------
IPython can display richer representations of objects.
To do this, you can define ``_ipython_display_()``, or any of a number of 
``_repr_*_()`` methods. 
Note that these are surrounded by single, not double underscores.

.. list-table:: Supported ``_repr_*_`` methods
   :widths: 20 15 15 15
   :header-rows: 1
   
   * - Format
     - REPL
     - Notebook
     - Qt Console
   * - ``_repr_pretty_``
     - yes
     - yes
     - yes
   * - ``_repr_svg_``
     - no
     - yes
     - yes
   * - ``_repr_png_``
     - no
     - yes
     - yes
   * - ``_repr_jpeg_``
     - no
     - yes
     - yes
   * - ``_repr_html_``
     - no
     - yes
     - no
   * - ``_repr_javascript_``
     - no
     - yes
     - no
   * - ``_repr_markdown_``
     - no
     - yes
     - no
   * - ``_repr_latex_``
     - no
     - yes
     - no
   * - ``_repr_mimebundle_``
     - no
     - ?
     - ?

If the methods don't exist, or return ``None``, the standard ``repr()`` is used.

For example::

    class Shout(object):
        def __init__(self, text):
            self.text = text
        
        def _repr_html_(self):
            return "<h1>" + self.text + "</h1>"


Special methods
^^^^^^^^^^^^^^^

Pretty printing
"""""""""""""""

To customize how your object is pretty-printed, add a ``_repr_pretty_`` method
to the class. 
The method should accept a pretty printer, and a boolean that indicates whether
the printer detected a cycle.
The method should act on the printer to produce your customized pretty output.
Here is an example::

    class MyObject(object):

        def _repr_pretty_(self, p, cycle):
            if cycle:
                p.text('MyObject(...)')
            else:
                p.text('MyObject[...]')

For details on how to use the pretty printer, see :py:mod:`IPython.lib.pretty`.

More powerful methods
"""""""""""""""""""""

.. class:: MyObject

   .. method:: _repr_mimebundle_(include=None, exclude=None)

     Should return a dictionary of multiple formats, keyed by mimetype, or a tuple
     of two dictionaries: *data, metadata* (see :ref:`Metadata`).
     If this returns something, other ``_repr_*_`` methods are ignored.
     The method should take keyword arguments ``include`` and ``exclude``, though 
     it is not required to respect them.

   .. method:: _ipython_display_()

      Displays the object as a side effect; the return value is ignored. If this
      is defined, all other display methods are ignored.
      This method is ignored in the REPL.


Metadata
^^^^^^^^

We often want to provide frontends with guidance on how to display the data. To
support this, ``_repr_*_()`` methods (except `_repr_pretty_``?) can also return a ``(data, metadata)``
tuple where ``metadata`` is a dictionary containing arbitrary key-value pairs for
the frontend to interpret. An example use case is ``_repr_jpeg_()``, which can
be set to return a jpeg image and a ``{'height': 400, 'width': 600}`` dictionary
to inform the frontend how to size the image.



Formatters for third-party types
--------------------------------

The user can also register formatters for types without modifying the class::

    from bar.baz import Foo

    def foo_html(obj):
        return '<marquee>Foo object %s</marquee>' % obj.name

    html_formatter = get_ipython().display_formatter.formatters['text/html']
    html_formatter.for_type(Foo, foo_html)

    # Or register a type without importing it - this does the same as above:
    html_formatter.for_type_by_name('bar.baz', 'Foo', foo_html)

Custom exception tracebacks
===========================

Rarely, you might want to display a custom traceback when reporting an
exception. To do this, define the custom traceback using
`_render_traceback_(self)` method which returns a list of strings, one string
for each line of the traceback. For example, the `ipyparallel
<https://ipyparallel.readthedocs.io/>`__ a parallel computing framework for
IPython, does this to display errors from multiple engines.

Please be conservative in using this feature; by replacing the default traceback
you may hide important information from the user.

===========================
Custom input transformation
===========================

IPython extends Python syntax to allow things like magic commands, and help with
the ``?`` syntax. There are several ways to customise how the user's input is
processed into Python code to be executed.

These hooks are mainly for other projects using IPython as the core of their
interactive interface. Using them carelessly can easily break IPython!

String based transformations
============================

.. currentmodule:: IPython.core.inputtransforms

When the user enters code, it is first processed as a string. By the
end of this stage, it must be valid Python syntax.

.. versionchanged:: 7.0

   The API for string and token-based transformations has been completely
   redesigned. Any third party code extending input transformation will need to
   be rewritten. The new API is, hopefully, simpler.

String based transformations are functions which accept a list of strings:
each string is a single line of the input cell, including its line ending.
The transformation function should return output in the same structure.

These transformations are in two groups, accessible as attributes of
the :class:`~IPython.core.interactiveshell.InteractiveShell` instance.
Each group is a list of transformation functions.

* ``input_transformers_cleanup`` run first on input, to do things like stripping
  prompts and leading indents from copied code. It may not be possible at this
  stage to parse the input as valid Python code.
* Then IPython runs its own transformations to handle its special syntax, like
  ``%magics`` and ``!system`` commands. This part does not expose extension
  points.
* ``input_transformers_post`` run as the last step, to do things like converting
  float literals into decimal objects. These may attempt to parse the input as
  Python code.  

These transformers may raise :exc:`SyntaxError` if the input code is invalid, but
in most cases it is clearer to pass unrecognised code through unmodified and let
Python's own parser decide whether it is valid.

For example, imagine we want to obfuscate our code by reversing each line, so
we'd write ``)5(f =+ a`` instead of ``a += f(5)``. Here's how we could swap it
back the right way before IPython tries to run it::

    def reverse_line_chars(lines):
        new_lines = []
        for line in lines:
            chars = line[:-1]  # the newline needs to stay at the end
            new_lines.append(chars[::-1] + '\n')
        return new_lines

To start using this::

    ip = get_ipython()
    ip.input_transformers_cleanup.append(reverse_line_chars)

.. versionadded:: 7.17

    input_transformers can now have an attribute ``has_side_effects`` set to
    `True`, which will prevent the transformers from being ran when IPython is
    trying to guess whether the user input is complete. 



AST transformations
===================

After the code has been parsed as Python syntax, you can use Python's powerful
*Abstract Syntax Tree* tools to modify it. Subclass :class:`ast.NodeTransformer`,
and add an instance to ``shell.ast_transformers``.

This example wraps integer literals in an ``Integer`` class, which is useful for
mathematical frameworks that want to handle e.g. ``1/3`` as a precise fraction::


    class IntegerWrapper(ast.NodeTransformer):
        """Wraps all integers in a call to Integer()"""
        def visit_Num(self, node):
            if isinstance(node.n, int):
                return ast.Call(func=ast.Name(id='Integer', ctx=ast.Load()),
                                args=[node], keywords=[])
            return node
.. _config_index:

===============================
Configuration and customization
===============================

Configuring IPython
-------------------

.. toctree::
   :maxdepth: 2

   intro
   options/index
   shortcuts/index
   details

.. seealso::

   :doc:`/development/config`
      Technical details of the config system.

Extending and integrating with IPython
--------------------------------------

.. toctree::
   :maxdepth: 2

   extensions/index
   integrating
   custommagics
   shell_mimerenderer
   inputtransforms
   callbacks
   eventloops
.. _defining_magics:

Defining custom magics
======================

There are two main ways to define your own magic functions: from standalone
functions and by inheriting from a base class provided by IPython:
:class:`IPython.core.magic.Magics`. Below we show code you can place in a file
that you load from your configuration, such as any file in the ``startup``
subdirectory of your default IPython profile.

First, let us see the simplest case. The following shows how to create a line
magic, a cell one and one that works in both modes, using just plain functions:

.. sourcecode:: python

    from IPython.core.magic import (register_line_magic, register_cell_magic,
                                    register_line_cell_magic)

    @register_line_magic
    def lmagic(line):
        "my line magic"
        return line

    @register_cell_magic
    def cmagic(line, cell):
        "my cell magic"
        return line, cell

    @register_line_cell_magic
    def lcmagic(line, cell=None):
        "Magic that works both as %lcmagic and as %%lcmagic"
        if cell is None:
            print("Called as line magic")
            return line
        else:
            print("Called as cell magic")
            return line, cell

    # In an interactive session, we need to delete these to avoid
    # name conflicts for automagic to work on line magics.
    del lmagic, lcmagic


You can also create magics of all three kinds by inheriting from the
:class:`IPython.core.magic.Magics` class.  This lets you create magics that can
potentially hold state in between calls, and that have full access to the main
IPython object:
    
.. sourcecode:: python

    # This code can be put in any Python module, it does not require IPython
    # itself to be running already.  It only creates the magics subclass but
    # doesn't instantiate it yet.
    from __future__ import print_function
    from IPython.core.magic import (Magics, magics_class, line_magic,
                                    cell_magic, line_cell_magic)

    # The class MUST call this class decorator at creation time
    @magics_class
    class MyMagics(Magics):

        @line_magic
        def lmagic(self, line):
            "my line magic"
            print("Full access to the main IPython object:", self.shell)
            print("Variables in the user namespace:", list(self.shell.user_ns.keys()))
            return line

        @cell_magic
        def cmagic(self, line, cell):
            "my cell magic"
            return line, cell

        @line_cell_magic
        def lcmagic(self, line, cell=None):
            "Magic that works both as %lcmagic and as %%lcmagic"
            if cell is None:
                print("Called as line magic")
                return line
            else:
                print("Called as cell magic")
                return line, cell


    # In order to actually use these magics, you must register them with a
    # running IPython.

    def load_ipython_extension(ipython):
        """
        Any module file that define a function named `load_ipython_extension`
        can be loaded via `%load_ext module.path` or be configured to be
        autoloaded by IPython at startup time.
        """
        # You can register the class itself without instantiating it.  IPython will
        # call the default constructor on it.
        ipython.register_magics(MyMagics)

If you want to create a class with a different constructor that holds
additional state, then you should always call the parent constructor and
instantiate the class yourself before registration:

.. sourcecode:: python

    @magics_class
    class StatefulMagics(Magics):
        "Magics that hold additional state"

        def __init__(self, shell, data):
            # You must call the parent constructor
            super(StatefulMagics, self).__init__(shell)
            self.data = data
        
        # etc...

    def load_ipython_extension(ipython):
        """
        Any module file that define a function named `load_ipython_extension`
        can be loaded via `%load_ext module.path` or be configured to be
        autoloaded by IPython at startup time.
        """
        # This class must then be registered with a manually created instance,
        # since its constructor has different arguments from the default:
        magics = StatefulMagics(ipython, some_data)
        ipython.register_magics(magics)


.. note::

   In early IPython versions 0.12 and before the line magics were
   created using a :func:`define_magic` API function.  This API has been
   replaced with the above in IPython 0.13 and then completely removed
   in IPython 5.  Maintainers of IPython extensions that still use the
   :func:`define_magic` function are advised to adjust their code
   for the current API.


Accessing user namespace and local scope
========================================

When creating line magics, you may need to access surrounding scope  to get user
variables (e.g when called inside functions). IPython provide the
``@needs_local_scope`` decorator that can be imported from
``IPython.core.magics``. When decorated with ``@needs_local_scope`` a magic will
be passed ``local_ns`` as an argument. As a convenience ``@needs_local_scope``
can also be applied to cell magics even if cell magics cannot appear at local
scope context.

Complete Example
================

Here is a full example of a magic package. You can distribute magics using
setuptools, distutils, or any other distribution tools like `flit
<https://flit.readthedocs.io>`_ for pure Python packages.

When distributing magics as part of a package, recommended best practice is to
execute the registration inside the `load_ipython_extension` as demonstrated in
the example below, instead of directly in the module (as in the initial example
with the ``@register_*`` decorators). This means a user will need to explicitly
choose to load your magic with ``%load_ext``. instead implicitly getting it when
importing the module. This is particularly relevant if loading your magic has 
side effects, if it is slow to load, or if it might override another magic with
the same name. 

.. sourcecode:: bash

   .
   ├── example_magic
   │   ├── __init__.py
   │   └── abracadabra.py
   └── setup.py

.. sourcecode:: bash

   $ cat example_magic/__init__.py
   """An example magic"""
   __version__ = '0.0.1'
   
   from .abracadabra import Abracadabra
   
   def load_ipython_extension(ipython):
       ipython.register_magics(Abracadabra)

.. sourcecode:: bash

    $ cat example_magic/abracadabra.py
    from IPython.core.magic import (Magics, magics_class, line_magic, cell_magic)

    @magics_class
    class Abracadabra(Magics):

        @line_magic
        def abra(self, line):
            return line

        @cell_magic
        def cadabra(self, line, cell):
            return line, cell

=====================================
Introduction to IPython configuration
=====================================

.. _setting_config:

Setting configurable options
============================

Many of IPython's classes have configurable attributes (see
:doc:`options/index` for the list). These can be
configured in several ways.

Python configuration files
--------------------------

To create the blank configuration files, run::

    ipython profile create [profilename]

If you leave out the profile name, the files will be created for the
``default`` profile (see :ref:`profiles`). These will typically be located in
:file:`~/.ipython/profile_default/`, and will be named
:file:`ipython_config.py`, for historical reasons you may also find files
named with IPython prefix instead of Jupyter:
:file:`ipython_notebook_config.py`, etc. The settings in
:file:`ipython_config.py` apply to all IPython commands.

By default, configuration files are fully featured Python scripts that can
execute arbitrary code, the main usage is to set value on the configuration
object ``c`` which exist in your configuration file.

You can then configure class attributes like this::

    c.InteractiveShell.automagic = False

Be careful with spelling--incorrect names will simply be ignored, with
no error. 

To add to a collection which may have already been defined elsewhere or have
default values, you can use methods like those found on lists, dicts and
sets: append, extend, :meth:`~traitlets.config.LazyConfigValue.prepend` (like
extend, but at the front), add and update (which works both for dicts and
sets)::

    c.InteractiveShellApp.extensions.append('Cython')

.. versionadded:: 2.0
   list, dict and set methods for config values

Example configuration file
``````````````````````````

::

    # sample ipython_config.py

    c.TerminalIPythonApp.display_banner = True
    c.InteractiveShellApp.log_level = 20
    c.InteractiveShellApp.extensions = [
        'myextension'
    ]
    c.InteractiveShellApp.exec_lines = [
        'import numpy',
        'import scipy'
    ]
    c.InteractiveShellApp.exec_files = [
        'mycode.py',
        'fancy.ipy'
    ]
    c.InteractiveShell.colors = 'LightBG'
    c.InteractiveShell.xmode = 'Context'
    c.TerminalInteractiveShell.confirm_exit = False
    c.TerminalInteractiveShell.editor = 'nano'

    c.PrefilterManager.multi_line_specials = True

    c.AliasManager.user_aliases = [
     ('la', 'ls -al')
    ]

JSON Configuration files
------------------------

In case where executability of configuration can be problematic, or
configurations need to be modified programmatically, IPython also support a
limited set of functionalities via ``.json`` configuration files. 

You can defined most of the configuration options via a json object which
hierarchy represent the value you would normally set on the ``c`` object of
``.py`` configuration files. The following ``ipython_config.json`` file::

    {
        "InteractiveShell": {
            "colors": "LightBG",
        },
        "InteractiveShellApp": {
            "extensions": [
                "myextension"
            ]
        },
        "TerminalInteractiveShell": {
            "editor": "nano"
        }
    }

Is equivalent to the following ``ipython_config.py``::

    c.InteractiveShellApp.extensions = [
        'myextension'
    ]

    c.InteractiveShell.colors = 'LightBG'
    c.TerminalInteractiveShell.editor = 'nano'


Command line arguments
----------------------

Every configurable value can be set from the command line, using this
syntax::

    ipython --ClassName.attribute=value

Many frequently used options have short aliases and flags, such as
``--matplotlib`` (to integrate with a matplotlib GUI event loop) or
``--pdb`` (automatic post-mortem debugging of exceptions).

To see all of these abbreviated options, run::

    ipython --help
    jupyter notebook --help
    # etc.

Options specified at the command line, in either format, override
options set in a configuration file.

The config magic
----------------

You can also modify config from inside IPython, using a magic command::

    %config IPCompleter.greedy = True

At present, this only affects the current session - changes you make to
config are not saved anywhere. Also, some options are only read when
IPython starts, so they can't be changed like this.

.. _configure_start_ipython:

Running IPython from Python
----------------------------

If you are using :ref:`embedding` to start IPython from a normal 
python file, you can set configuration options the same way as in a 
config file by creating a traitlets config object and passing it to 
start_ipython like in the example below.

.. literalinclude:: ../../../examples/Embedding/start_ipython_config.py
    :language: python

.. _profiles:

Profiles
========

IPython can use multiple profiles, with separate configuration and
history. By default, if you don't specify a profile, IPython always runs
in the ``default`` profile. To use a new profile::

    ipython profile create foo   # create the profile foo
    ipython --profile=foo        # start IPython using the new profile

Profiles are typically stored in :ref:`ipythondir`, but you can also keep
a profile in the current working directory, for example to distribute it
with a project. To find a profile directory on the filesystem::

    ipython locate profile foo

.. _ipythondir:

The IPython directory
=====================

IPython stores its files---config, command history and extensions---in
the directory :file:`~/.ipython/` by default.

.. envvar:: IPYTHONDIR

   If set, this environment variable should be the path to a directory,
   which IPython will use for user data. IPython will create it if it
   does not exist.

.. option:: --ipython-dir=<path>

   This command line option can also be used to override the default
   IPython directory.

To see where IPython is looking for the IPython directory, use the command
``ipython locate``, or the Python function :func:`IPython.paths.get_ipython_dir`.


Systemwide configuration
========================

It can be useful to deploy systemwide ipython or ipykernel configuration
when managing environment for many users. At startup time IPython and
IPykernel will search for configuration file in multiple systemwide
locations, mainly:

  - ``/etc/ipython/``
  - ``/usr/local/etc/ipython/``

When the global install is a standalone python distribution it may also
search in distribution specific location, for example:

  - ``$ANACONDA_LOCATION/etc/ipython/``

In those locations, Terminal IPython will look for a file called
``ipython_config.py`` and ``ipython_config.json``, ipykernel will look for
``ipython_kernel_config.py`` and ``ipython_kernel.json``.

Configuration files are loaded in order and merged with configuration on
later location taking precedence on earlier locations (that is to say a user
can overwrite a systemwide configuration option).

You can see all locations in which IPython is looking for configuration files
by starting ipython in debug mode::

    $ ipython --debug -c 'exit()'

Identically with ipykernel though the command is currently blocking until
this process is killed with ``Ctrl-\``::
 
    $ python -m ipykernel --debug
=======================
Specific config details
=======================

.. _custom_prompts:

Custom Prompts
==============

.. versionchanged:: 5.0

From IPython 5, prompts are produced as a list of Pygments tokens, which are
tuples of (token_type, text). You can customise prompts by writing a method
which generates a list of tokens.

There are four kinds of prompt:

* The **in** prompt is shown before the first line of input
  (default like ``In [1]:``).
* The **continuation** prompt is shown before further lines of input
  (default like ``...:``).
* The **rewrite** prompt is shown to highlight how special syntax has been
  interpreted (default like ``----->``).
* The **out** prompt is shown before the result from evaluating the input
  (default like ``Out[1]:``).

Custom prompts are supplied together as a class. If you want to customise only
some of the prompts, inherit from :class:`IPython.terminal.prompts.Prompts`,
which defines the defaults. The required interface is like this:

.. class:: MyPrompts(shell)

   Prompt style definition. *shell* is a reference to the
   :class:`~.TerminalInteractiveShell` instance.

   .. method:: in_prompt_tokens(cli=None)
               continuation_prompt_tokens(self, cli=None, width=None)
               rewrite_prompt_tokens()
               out_prompt_tokens()

      Return the respective prompts as lists of ``(token_type, text)`` tuples.

      For continuation prompts, *width* is an integer representing the width of
      the prompt area in terminal columns.

      *cli*, where used, is the prompt_toolkit ``CommandLineInterface`` instance.
      This is mainly for compatibility with the API prompt_toolkit expects.

Here is an example Prompt class that will show the current working directory
in the input prompt:

.. code-block:: python

    from IPython.terminal.prompts import Prompts, Token
    import os

    class MyPrompt(Prompts):
         def in_prompt_tokens(self, cli=None):
             return [(Token, os.getcwd()),
                     (Token.Prompt, ' >>>')]

To set the new prompt, assign it to the ``prompts`` attribute of the IPython
shell:

.. code-block:: python

    In [2]: ip = get_ipython()
       ...: ip.prompts = MyPrompt(ip)

    /home/bob >>> # it works

See ``IPython/example/utils/cwd_prompt.py`` for an example of how to write an
extensions to customise prompts.

Inside IPython or in a startup script, you can use a custom prompts class
by setting ``get_ipython().prompts`` to an *instance* of the class.
In configuration, ``TerminalInteractiveShell.prompts_class`` may be set to
either the class object, or a string of its full importable name.

To include invisible terminal control sequences in a prompt, use
``Token.ZeroWidthEscape`` as the token type. Tokens with this type are ignored
when calculating the width.

Colours in the prompt are determined by the token types and the highlighting
style; see below for more details. The tokens used in the default prompts are
``Prompt``, ``PromptNum``, ``OutPrompt`` and ``OutPromptNum``.

.. _termcolour:

Terminal Colors
===============

.. versionchanged:: 5.0

There are two main configuration options controlling colours.

``InteractiveShell.colors`` sets the colour of tracebacks and object info (the
output from e.g. ``zip?``). It may also affect other things if the option below
is set to ``'legacy'``. It has four case-insensitive values:
``'nocolor', 'neutral', 'linux', 'lightbg'``. The default is *neutral*, which
should be legible on either dark or light terminal backgrounds. *linux* is
optimised for dark backgrounds and *lightbg* for light ones.

``TerminalInteractiveShell.highlighting_style`` determines prompt colours and
syntax highlighting. It takes the name (as a string) or class (as a subclass of
``pygments.style.Style``) of a Pygments style, or the special value ``'legacy'``
to pick a style in accordance with ``InteractiveShell.colors``.

You can see the Pygments styles available on your system by running::

    import pygments
    list(pygments.styles.get_all_styles())

Additionally, ``TerminalInteractiveShell.highlighting_style_overrides`` can override
specific styles in the highlighting. It should be a dictionary mapping Pygments
token types to strings defining the style. See `Pygments' documentation
<http://pygments.org/docs/styles/#creating-own-styles>`__ for the language used
to define styles.

Colors in the pager
-------------------

On some systems, the default pager has problems with ANSI colour codes.
To configure your default pager to allow these:

1. Set the environment PAGER variable to ``less``.
2. Set the environment LESS variable to ``-r`` (plus any other options
   you always want to pass to less by default). This tells less to
   properly interpret control sequences, which is how color
   information is given to your terminal.

.. _editors:

Editor configuration
====================

IPython can integrate with text editors in a number of different ways:

* Editors (such as `(X)Emacs`_, vim_ and TextMate_) can
  send code to IPython for execution.

* IPython's ``%edit`` magic command can open an editor of choice to edit
  a code block.

The %edit command (and its alias %ed) will invoke the editor set in your
environment as :envvar:`EDITOR`. If this variable is not set, it will default
to vi under Linux/Unix and to notepad under Windows. You may want to set this
variable properly and to a lightweight editor which doesn't take too long to
start (that is, something other than a new instance of Emacs). This way you
can edit multi-line code quickly and with the power of a real editor right
inside IPython.

You can also control the editor by setting :attr:`TerminalInteractiveShell.editor`
in :file:`ipython_config.py`.

Vim
---

Paul Ivanov's `vim-ipython <https://github.com/ivanov/vim-ipython>`_ provides
powerful IPython integration for vim.

.. _emacs:

(X)Emacs
--------

If you are a dedicated Emacs user, and want to use Emacs when IPython's
``%edit`` magic command is called you should set up the Emacs server so that
new requests are handled by the original process. This means that almost no
time is spent in handling the request (assuming an Emacs process is already
running). For this to work, you need to set your EDITOR environment variable
to 'emacsclient'. The code below, supplied by Francois Pinard, can then be
used in your :file:`.emacs` file to enable the server:

.. code-block:: common-lisp

    (defvar server-buffer-clients)
    (when (and (fboundp 'server-start) (string-equal (getenv "TERM") 'xterm))
      (server-start)
      (defun fp-kill-server-with-buffer-routine ()
        (and server-buffer-clients (server-done)))
      (add-hook 'kill-buffer-hook 'fp-kill-server-with-buffer-routine))

Thanks to the work of Alexander Schmolck and Prabhu Ramachandran,
currently (X)Emacs and IPython get along very well in other ways.

With (X)EMacs >= 24, You can enable IPython in python-mode with:

.. code-block:: common-lisp

    (require 'python)
    (setq python-shell-interpreter "ipython")

.. _`(X)Emacs`: http://www.gnu.org/software/emacs/
.. _TextMate: http://macromates.com/
.. _vim: http://www.vim.org/

.. _custom_keyboard_shortcuts:

Keyboard Shortcuts
==================

.. versionchanged:: 5.0

You can customise keyboard shortcuts for terminal IPython. Put code like this in
a :ref:`startup file <startup_files>`::

    from IPython import get_ipython
    from prompt_toolkit.enums import DEFAULT_BUFFER
    from prompt_toolkit.keys import Keys
    from prompt_toolkit.filters import HasFocus, HasSelection, ViInsertMode, EmacsInsertMode

    ip = get_ipython()
    insert_mode = ViInsertMode() | EmacsInsertMode()

    def insert_unexpected(event):
        buf = event.current_buffer
        buf.insert_text('The Spanish Inquisition')
    # Register the shortcut if IPython is using prompt_toolkit
    if getattr(ip, 'pt_app', None):
        registry = ip.pt_app.key_bindings
        registry.add_binding(Keys.ControlN,
                         filter=(HasFocus(DEFAULT_BUFFER)
                                 & ~HasSelection()
                                 & insert_mode))(insert_unexpected)


Here is a second example that bind the key sequence ``j``, ``k`` to switch to
VI input mode to ``Normal`` when in insert mode::

   from IPython import get_ipython
   from prompt_toolkit.enums import DEFAULT_BUFFER
   from prompt_toolkit.filters import HasFocus, ViInsertMode
   from prompt_toolkit.key_binding.vi_state import InputMode

   ip = get_ipython()

   def switch_to_navigation_mode(event):
      vi_state = event.cli.vi_state
      vi_state.input_mode = InputMode.NAVIGATION

   if getattr(ip, 'pt_app', None):
      registry = ip.pt_app.key_bindings
      registry.add_binding(u'j',u'k',
                           filter=(HasFocus(DEFAULT_BUFFER)
                                    & ViInsertMode()))(switch_to_navigation_mode)

For more information on filters and what you can do with the ``event`` object,
`see the prompt_toolkit docs
<https://python-prompt-toolkit.readthedocs.io/en/latest/pages/asking_for_input.html#adding-custom-key-bindings>`__.


Enter to execute
----------------

In the Terminal IPython shell – which by default uses the ``prompt_toolkit``
interface, the semantic meaning of pressing the :kbd:`Enter` key can be
ambiguous. In some case :kbd:`Enter` should execute code, and in others it
should add a new line. IPython uses heuristics to decide whether to execute or
insert a new line at cursor position. For example, if we detect that the current
code is not valid Python, then the user is likely editing code and the right
behavior is to likely to insert a new line. If the current code is a simple
statement like `ord('*')`, then the right behavior is likely to execute. Though
the exact desired semantics often varies from users to users.

As the exact behavior of :kbd:`Enter` is ambiguous, it has been special cased
to allow users to completely configure the behavior they like. Hence you can
have enter always execute code. If you prefer fancier behavior, you need to get
your hands dirty and read the ``prompt_toolkit`` and IPython documentation
though. See :ghpull:`10500`, set the
``c.TerminalInteractiveShell.handle_return`` option and get inspiration from the
following example that only auto-executes the input if it begins with a bang or
a modulo character (``!`` or ``%``). To use the following code, add it to your
IPython configuration::

    def custom_return(shell):

        """This function is required by the API. It takes a reference to
        the shell, which is the same thing `get_ipython()` evaluates to.
        This function must return a function that handles each keypress
        event. That function, named `handle` here, references `shell`
        by closure."""

        def handle(event):

            """This function is called each time `Enter` is pressed,
            and takes a reference to a Prompt Toolkit event object.
            If the current input starts with a bang or modulo, then
            the input is executed, otherwise a newline is entered,
            followed by any spaces needed to auto-indent."""

            # set up a few handy references to nested items...

            buffer = event.current_buffer
            document = buffer.document
            text = document.text

            if text.startswith('!') or text.startswith('%'): # execute the input...

                buffer.accept_action.validate_and_handle(event.cli, buffer)

            else: # insert a newline with auto-indentation...

                if document.line_count > 1: text = text[:document.cursor_position]
                indent = shell.check_complete(text)[1]
                buffer.insert_text('\n' + indent)
            
                # if you just wanted a plain newline without any indentation, you
                # could use `buffer.insert_text('\n')` instead of the lines above

        return handle

    c.TerminalInteractiveShell.handle_return = custom_return
=================
IPython shortcuts
=================

Available shortcuts in an IPython terminal.

.. warning::

  This list is automatically generated, and may not hold all available
  shortcuts. In particular, it may depend on the version of ``prompt_toolkit``
  installed during the generation of this page.


Single Filtered shortcuts
=========================

.. csv-table::
    :header: Shortcut,Filter,Description
    :widths: 30, 30, 100
    :delim: tab
    :file: single_filtered.csv


Multi Filtered shortcuts
========================

.. csv-table::
    :header: Shortcut,Filter,Description
    :widths: 30, 30, 100
    :delim: tab
    :file: multi_filtered.csv
===============
IPython options
===============

Any of the options listed here can be set in config files, at the
command line, or from inside IPython. See :ref:`setting_config` for
details.

.. toctree::

   terminal
   kernel
.. _extensions_autoreload:

==========
autoreload
==========

.. magic:: autoreload

.. automodule:: IPython.extensions.autoreload
.. _extensions_overview:

==================
IPython extensions
==================

A level above configuration are IPython extensions, Python modules which modify
the behaviour of the shell. They are referred to by an importable module name,
and can be placed anywhere you'd normally import from, or in
``.ipython/extensions/``.

Getting extensions
==================

A few important extensions are :ref:`bundled with IPython <bundled_extensions>`.
Others can be found on the `extensions index
<https://github.com/ipython/ipython/wiki/Extensions-Index>`_ on the wiki, and
the `Framework :: IPython tag <https://pypi.python.org/pypi?:action=browse&c=586>`_
on PyPI.

Extensions on PyPI can be installed using ``pip``, like any other Python package.

Using extensions
================

To load an extension while IPython is running, use the ``%load_ext`` magic:

.. sourcecode:: ipython

    In [1]: %load_ext myextension

To load it each time IPython starts, list it in your configuration file::

    c.InteractiveShellApp.extensions = [
        'myextension'
    ]

Writing extensions
==================

An IPython extension is an importable Python module that has a couple of special
functions to load and unload it. Here is a template::

    # myextension.py

    def load_ipython_extension(ipython):
        # The `ipython` argument is the currently active `InteractiveShell`
        # instance, which can be used in any way. This allows you to register
        # new magics or aliases, for example.

    def unload_ipython_extension(ipython):
        # If you want your extension to be unloadable, put that logic here.

This :func:`load_ipython_extension` function is called after your extension is
imported, and the currently active :class:`~IPython.core.interactiveshell.InteractiveShell`
instance is passed as the only argument. You can do anything you want with
IPython at that point.

:func:`load_ipython_extension` will not be called again if the user use
`%load_extension`.  The user have to explicitly ask the extension to be
reloaded (with `%reload_extension`). In case where the use ask the extension to
be reloaded, , the extension will be unloaded (with
`unload_ipython_extension`), and loaded again. 

Useful :class:`InteractiveShell` methods include :meth:`~IPython.core.interactiveshell.InteractiveShell.register_magic_function`, 
:meth:`~IPython.core.interactiveshell.InteractiveShell.push` (to add variables to the user namespace) and 
:meth:`~IPython.core.interactiveshell.InteractiveShell.drop_by_id` (to remove variables on unloading).

.. seealso::

   :ref:`defining_magics`

You can put your extension modules anywhere you want, as long as they can be
imported by Python's standard import mechanism. However, to make it easy to
write extensions, you can also put your extensions in :file:`extensions/`
within the :ref:`IPython directory <ipythondir>`. This directory is
added to :data:`sys.path` automatically.

When your extension is ready for general use, please add it to the `extensions
index <https://github.com/ipython/ipython/wiki/Extensions-Index>`_. We also
encourage you to upload it to PyPI and use the ``Framework :: IPython``
classifier, so that users can install it with standard packaging tools.

.. _bundled_extensions:

Extensions bundled with IPython
===============================

.. toctree::
   :maxdepth: 1

   autoreload
   storemagic

* ``octavemagic`` used to be bundled, but is now part of `oct2py <https://blink1073.github.io/oct2py/>`_.
  Use ``%load_ext oct2py.ipython`` to load it.
* ``rmagic`` is now part of `rpy2 <http://rpy.sourceforge.net/>`_. Use
  ``%load_ext rpy2.ipython`` to load it, and see :mod:`rpy2.ipython.rmagic` for
  details of how to use it.
* ``cythonmagic`` used to be bundled, but is now part of `cython <https://github.com/cython/cython/>`_
  Use ``%load_ext Cython`` to load it.
* ``sympyprinting`` used to be a bundled extension, but you should now use
  :func:`sympy.init_printing` instead.
.. _extensions_storemagic:

==========
storemagic
==========

.. automodule:: IPython.extensions.storemagic

.. automethod:: StoreMagics.store
=====================
 Development version
=====================

This document describes in-flight development work.

.. warning::

    Please do not edit this file by hand (doing so will likely cause merge
    conflicts for other Pull Requests). Instead, create a new file in the
    `docs/source/whatsnew/pr` folder


Released .... ...., 2019


Need to be updated:

.. toctree::
   :maxdepth: 2
   :glob:

   pr/*




.. DO NOT EDIT THIS LINE BEFORE RELEASE. FEATURE INSERTION POINT.

Backwards incompatible changes
------------------------------

.. DO NOT EDIT THIS LINE BEFORE RELEASE. INCOMPAT INSERTION POINT.
============
 6.x Series
============

.. _whatsnew650:

IPython 6.5.0
=============

Miscellaneous bug fixes and compatibility with Python 3.7.

* Autocompletion fix for modules with out ``__init__.py`` :ghpull:`11227`
* update the ``%pastebin`` magic to use ``dpaste.com`` instead og GitHub Gist
  which now requires authentication :ghpull:`11182`
* Fix crash with multiprocessing :ghpull:`11185`

.. _whatsnew640:

IPython 6.4.0
=============

Everything new in :ref:`IPython 5.7 <whatsnew570>`

* Fix display object not emitting metadata :ghpull:`11106`
* Comments failing Jedi test :ghpull:`11110` 


.. _whatsnew631:

IPython 6.3.1
=============

This is a bugfix release to switch the default completions back to IPython's
own completion machinery. We discovered some problems with the completions
from Jedi, including completing column names on pandas data frames.

You can switch the completions source with the config option
:configtrait:`Completer.use_jedi`.

.. _whatsnew630:

IPython 6.3
===========

IPython 6.3 contains all the bug fixes and features in
:ref:`IPython 5.6 <whatsnew560>`. In addition:

* A new display class :class:`IPython.display.Code` can be used to display
  syntax highlighted code in a notebook (:ghpull:`10978`).
* The :cellmagic:`html` magic now takes a ``--isolated`` option to put the
  content in an iframe (:ghpull:`10962`).
* The code to find completions using the Jedi library has had various
  adjustments. This is still a work in progress, but we hope this version has
  fewer annoyances (:ghpull:`10956`, :ghpull:`10969`, :ghpull:`10999`,
  :ghpull:`11035`, :ghpull:`11063`, :ghpull:`11065`).
* The *post* event callbacks are now always called, even when the execution failed
  (for example because of a ``SyntaxError``).
* The execution info and result objects are now made available in the
  corresponding *pre* or *post* ``*_run_cell`` :doc:`event callbacks </config/callbacks>`
  in a backward compatible manner (:ghissue:`10774` and :ghpull:`10795`).
* Performance with very long code cells (hundreds of lines) is greatly improved
  (:ghpull:`10898`). Further improvements are planned for IPython 7.

You can see all `pull requests for the 6.3 milestone
<https://github.com/ipython/ipython/pulls?utf8=%E2%9C%93&q=is%3Apr+milestone%3A6.3+is%3Aclosed>`__.

.. _whatsnew620:

IPython 6.2
===========

IPython 6.2 contains all the bugs fixes and features :ref:`available in IPython 5.5 <whatsnew550>`,
like built in progress bar support, and system-wide configuration

The following features are specific to IPython 6.2:

Function signature in completions
---------------------------------

Terminal IPython will now show the signature of the function while completing.
Only the currently highlighted function will show its signature on the line
below the completer by default. This functionality is recent, so it might be
limited; we welcome bug reports and requests for enhancements. :ghpull:`10507`

Assignments return values
-------------------------

IPython can now trigger the display hook on the last assignment of cells.
Up until 6.2 the following code wouldn't show the value of the assigned
variable::

    In[1]: xyz = "something"
    # nothing shown

You would have to actually make it the last statement::

    In [2]: xyz = "something else"
    ...   : xyz
    Out[2]: "something else"

With the option ``InteractiveShell.ast_node_interactivity='last_expr_or_assign'``
you can now do::

    In [2]: xyz = "something else"
    Out[2]: "something else"

This option can be toggled at runtime with the ``%config`` magic, and will
trigger on assignment ``a = 1``, augmented assignment ``+=``, ``-=``, ``|=`` ...
as well as type annotated assignments: ``a:int = 2``.

See :ghpull:`10598`

Recursive Call of ipdb
----------------------

Advanced users of the debugger can now correctly recursively enter ipdb. This is
thanks to ``@segevfiner`` on :ghpull:`10721`.

.. _whatsnew610:

IPython 6.1
===========

- Quotes in a filename are always escaped during tab-completion on non-Windows.
  :ghpull:`10069`

- Variables now shadow magics in autocompletion. See :ghissue:`4877` and :ghpull:`10542`.

- Added the ability to add parameters to alias_magic. For example::

    In [2]: %alias_magic hist history --params "-l 2" --line
    Created `%hist` as an alias for `%history -l 2`.

    In [3]: hist
    %alias_magic hist history --params "-l 30" --line
    %alias_magic hist history --params "-l 2" --line

  Previously it was only possible to have an alias attached to a single function,
  and you would have to pass in the given parameters every time::

    In [4]: %alias_magic hist history --line
    Created `%hist` as an alias for `%history`.

    In [5]: hist -l 2
    hist
    %alias_magic hist history --line

- To suppress log state messages, you can now either use ``%logstart -q``, pass
  ``--LoggingMagics.quiet=True`` on the command line, or set
  ``c.LoggingMagics.quiet=True`` in your configuration file.

- An additional flag ``--TerminalInteractiveShell.term_title_format`` is
  introduced to allow the user to control the format of the terminal title.  It
  is specified as a python format string, and currently the only variable it
  will format is ``{cwd}``.

- ``??``/``%pinfo2`` will now show object docstrings if the source can't be retrieved. :ghpull:`10532`
- ``IPython.display`` has gained a ``%markdown`` cell magic. :ghpull:`10563`
- ``%config`` options can now be tab completed. :ghpull:`10555`
- ``%config`` with no arguments are now unique and sorted. :ghpull:`10548`
- Completion on keyword arguments does not duplicate ``=`` sign if already present. :ghpull:`10547`
- ``%run -m <module>`` now ``<module>`` passes extra arguments to ``<module>``. :ghpull:`10546`
- completer now understand "snake case auto complete": if ``foo_bar_kittens`` is
  a valid completion, I can type ``f_b<tab>`` will complete to it. :ghpull:`10537`
- tracebacks are better standardized and will compress `/path/to/home` to `~`. :ghpull:`10515`

The following changes were also added to IPython 5.4, see :ref:`what's new in IPython 5.4 <whatsnew540>`
for more detail description:

- ``TerminalInteractiveShell`` is configurable and can be configured to
  (re)-use the readline interface.

- objects can now define a ``_repr_mimebundle_``

- Execution heuristics improve for single line statements
- ``display()`` can now return a display id to update display areas.


.. _whatsnew600:

IPython 6.0
===========

Released April 19th, 2017

IPython 6 features a major improvement in the completion machinery which is now
capable of completing non-executed code. It is also the first version of IPython
to stop compatibility with Python 2, which is still supported on the bugfix only
5.x branch. Read below for a non-exhaustive list of new features.

Make sure you have pip > 9.0 before upgrading.
You should be able to update by using:

.. code::

    pip install ipython --upgrade


.. note::

    If your pip version is greater than or equal to pip 9.0.1 you will automatically get
    the most recent version of IPython compatible with your system: on Python 2 you 
    will get the latest IPython 5.x bugfix, while in Python 3
    you will get the latest 6.x stable version.

New completion API and Interface
--------------------------------

The completer Completion API has seen an overhaul, and the new completer has
plenty of improvements both from the end users of terminal IPython and for
consumers of the API.

This new API is capable of pulling completions from :any:`jedi`, thus allowing
type inference on non-executed code. If :any:`jedi` is installed, completions like
the following are now possible without code evaluation:

    >>> data = ['Number of users', 123_456]
    ... data[0].<tab>

That is to say, IPython is now capable of inferring that `data[0]` is a string,
and will suggest completions like `.capitalize`. The completion power of IPython
will increase with new Jedi releases, and a number of bug-fixes and more completions
are already available on the development version of :any:`jedi` if you are curious.

With the help of prompt toolkit, types of completions can be shown in the
completer interface:

.. image:: ../_images/jedi_type_inference_60.png
    :alt: Jedi showing ability to do type inference
    :align: center
    :width: 400px
    :target: ../_images/jedi_type_inference_60.png

The appearance of the completer is controlled by the
``c.TerminalInteractiveShell.display_completions`` option that will show the
type differently depending on the value among ``'column'``, ``'multicolumn'``
and ``'readlinelike'``

The use of Jedi also fulfills a number of requests and fixes a number of bugs
like case-insensitive completion and completion after division operator: See
:ghpull:`10182`.

Extra patches and updates will be needed to the :mod:`ipykernel` package for
this feature to be available to other clients like Jupyter Notebook, Lab,
Nteract, Hydrogen...

The use of Jedi should be barely noticeable on recent machines, but 
can be slower on older ones.  To tweak the performance, the amount
of time given to Jedi to compute type inference can be adjusted with
``c.IPCompleter.jedi_compute_type_timeout``. The objects whose type were not
inferred will be shown as ``<unknown>``. Jedi can also be completely deactivated
by using the ``c.Completer.use_jedi=False`` option.


The old ``Completer.complete()`` API is waiting deprecation and should be
replaced replaced by ``Completer.completions()`` in the near future. Feedback on
the current state of the API and suggestions are welcome.

Python 3 only codebase
----------------------

One of the large challenges in IPython 6.0 has been the adoption of a pure
Python 3 codebase, which has led to upstream patches in pip,
pypi and warehouse to make sure Python 2 systems still upgrade to the latest
compatible Python version.

We remind our Python 2 users that IPython 5 is still compatible with Python 2.7,
still maintained and will get regular releases. Using pip 9+, upgrading IPython will
automatically upgrade to the latest version compatible with your system.

.. warning::

  If you are on a system using an older version of pip on Python 2, pip may
  still install IPython 6.0 on your system, and IPython will refuse to start.
  You can fix this by upgrading pip, and reinstalling ipython, or forcing pip to
  install an earlier version: ``pip install 'ipython<6'``

The ability to use only Python 3 on the code base of IPython brings a number
of advantages. Most of the newly written code make use of `optional function type
annotation <https://www.python.org/dev/peps/pep-0484/>`_ leading to clearer code
and better documentation.

The total size of the repository has also decreased by about 1500 lines (for the 
first time excluding the big split for 4.0). The decrease is potentially
a bit more for the sour as some documents like this one are append only and
are about 300 lines long.

The removal of the Python2/Python3 shim layer has made the code quite a lot clearer and
more idiomatic in a number of locations, and much friendlier to work with and
understand. We hope to further embrace Python 3 capabilities in the next release
cycle and introduce more of the Python 3 only idioms (yield from, kwarg only,
general unpacking) in the IPython code base, and see if we can take advantage
of these to improve user experience with better error messages and
hints.


Configurable TerminalInteractiveShell, readline interface
---------------------------------------------------------

IPython gained a new ``c.TerminalIPythonApp.interactive_shell_class`` option
that allows customizing the class used to start the terminal frontend. This
should allow a user to use custom interfaces, like reviving the former readline
interface which is now a separate package not actively maintained by the core
team. See the project to bring back the readline interface: `rlipython
<https://github.com/ipython/rlipython>`_.

This change will be backported to the IPython 5.x series.

Misc improvements
-----------------


- The :cellmagic:`capture` magic can now capture the result of a cell (from
  an expression on the last line), as well as printed and displayed output.
  :ghpull:`9851`.

- Pressing Ctrl-Z in the terminal debugger now suspends IPython, as it already
  does in the main terminal prompt.

- Autoreload can now reload ``Enum``. See :ghissue:`10232` and :ghpull:`10316`

- IPython.display has gained a :any:`GeoJSON <IPython.display.GeoJSON>` object.
  :ghpull:`10288` and :ghpull:`10253`

Functions Deprecated in 6.x Development cycle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Loading extensions from ``ipython_extension_dir`` prints a warning that this
  location is pending deprecation. This should only affect users still having
  extensions installed with ``%install_ext`` which has been deprecated since
  IPython 4.0, and removed in 5.0. Extensions still present in
  ``ipython_extension_dir`` may shadow more recently installed versions using
  pip. It is thus recommended to clean ``ipython_extension_dir`` of any
  extension now available as a package.


- ``IPython.utils.warn`` was deprecated in IPython 4.0, and has now been removed.
  instead of ``IPython.utils.warn`` inbuilt :any:`warnings` module is used.


- The function `IPython.core.oinspect.py:call_tip` is unused, was marked as
  deprecated (raising a `DeprecationWarning`) and marked for later removal.
  :ghpull:`10104`

Backward incompatible changes
------------------------------

Functions Removed in 6.x Development cycle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following functions have been removed in the
development cycle marked for Milestone 6.0.

- ``IPython/utils/process.py`` - ``is_cmd_found``
- ``IPython/utils/process.py`` - ``pycmd2argv``

- The `--deep-reload` flag and the corresponding options to inject `dreload` or
  `reload` into the interactive namespace have been removed. You have to
  explicitly import `reload` from `IPython.lib.deepreload` to use it.

- The :magic:`profile` used to print the current IPython profile, and which
  was deprecated in IPython 2.0 does now raise a `DeprecationWarning` error when
  used. It is often confused with the :magic:`prun` and the deprecation removal
  should free up the ``profile`` name in future versions.
.. _issues_list_100:

Issues closed in the 1.0 development cycle
==========================================


Issues closed in 1.2
--------------------

GitHub stats for 2013/09/09 - 2014/02/21

These lists are automatically generated, and may be incomplete or contain duplicates.

The following 13 authors contributed 84 commits.

* Benjamin Ragan-Kelley
* Daryl Herzmann
* Doug Blank
* Fernando Perez
* James Porter
* Juergen Hasch
* Julian Taylor
* Kyle Kelley
* Lawrence Fu
* Matthias Bussonnier
* Paul Ivanov
* Pascal Schetelat
* Puneeth Chaganti
* Takeshi Kanmae
* Thomas Kluyver

We closed a total of 55 issues, 38 pull requests and 17 regular issues;
this is the full list (generated with the script :file:`tools/github_stats.py`):

Pull Requests (38):

1.2.1:

* :ghpull:`4372`: Don't assume that SyntaxTB is always called with a SyntaxError
* :ghpull:`5166`: remove mktemp usage
* :ghpull:`5163`: Simplify implementation of TemporaryWorkingDirectory.
* :ghpull:`5105`: add index to format to support py2.6

1.2.0:

* :ghpull:`4972`: Work around problem in doctest discovery in Python 3.4 with PyQt
* :ghpull:`4934`: `ipython profile create` respects `--ipython-dir`
* :ghpull:`4845`: Add Origin Checking.
* :ghpull:`4928`: use importlib.machinery when available
* :ghpull:`4849`: Various unicode fixes (mostly on Windows)
* :ghpull:`4880`: set profile name from profile_dir
* :ghpull:`4908`: detect builtin docstrings in oinspect
* :ghpull:`4909`: sort dictionary keys before comparison, ordering is not guaranteed
* :ghpull:`4903`: use https for all embeds
* :ghpull:`4868`: Static path fixes
* :ghpull:`4820`: fix regex for cleaning old logs with ipcluster
* :ghpull:`4840`: Error in Session.send_raw()
* :ghpull:`4762`: whitelist alphanumeric characters for cookie_name
* :ghpull:`4748`: fix race condition in profiledir creation.
* :ghpull:`4720`: never use ssh multiplexer in tunnels
* :ghpull:`4738`: don't inject help into user_ns
* :ghpull:`4722`: allow purging local results as long as they are not outstanding
* :ghpull:`4668`: Make non-ASCII docstring unicode
* :ghpull:`4639`: Minor import fix to get qtconsole with --pylab=qt working
* :ghpull:`4453`: Play nice with App Nap
* :ghpull:`4609`: Fix bytes regex for Python 3.
* :ghpull:`4488`: fix typo in message spec doc
* :ghpull:`4346`: getpass() on Windows & Python 2 needs bytes prompt
* :ghpull:`4230`: Switch correctly to the user's default matplotlib backend after inline.
* :ghpull:`4214`: engine ID metadata should be unicode, not bytes
* :ghpull:`4232`: no highlight if no language specified
* :ghpull:`4218`: Fix display of SyntaxError when .py file is modified
* :ghpull:`4217`: avoid importing numpy at the module level
* :ghpull:`4213`: fixed dead link in examples/notebooks readme to Part 3
* :ghpull:`4183`: ESC should be handled by CM if tooltip is not on
* :ghpull:`4193`: Update for #3549: Append Firefox overflow-x fix
* :ghpull:`4205`: use TextIOWrapper when communicating with pandoc subprocess
* :ghpull:`4204`: remove some extraneous print statements from IPython.parallel
* :ghpull:`4201`: HeadingCells cannot be split or merged

1.2.1:

* :ghissue:`5101`: IPython 1.2.0: notebook fail with "500 Internal Server Error"

1.2.0:

* :ghissue:`4892`: IPython.qt test failure with python3.4
* :ghissue:`4810`: ipcluster bug in clean_logs flag
* :ghissue:`4765`: missing build script for highlight.js
* :ghissue:`4761`: ipv6 address triggers cookie exception
* :ghissue:`4721`: purge_results with jobid crashing - looking for insight
* :ghissue:`4602`: "ipcluster stop" fails after "ipcluster start --daemonize" using python3.3
* :ghissue:`3386`: Magic %paste not working in Python 3.3.2. TypeError: Type str doesn't support the buffer API
* :ghissue:`4485`: Incorrect info in "Messaging in IPython" documentation. 
* :ghissue:`4351`: /parallel/apps/launcher.py error
* :ghissue:`4334`: NotebookApp.webapp_settings static_url_prefix causes crash
* :ghissue:`4039`: Celltoolbar example issue
* :ghissue:`4256`: IPython no longer handles unicode file names 
* :ghissue:`4122`: Nbconvert [windows]: Inconsistent line endings in markdown cells exported to latex 
* :ghissue:`3819`: nbconvert add extra blank line to code block on Windows.
* :ghissue:`4203`: remove spurious print statement from parallel annoted functions
* :ghissue:`4200`: Notebook: merging a heading cell and markdown cell cannot be undone


Issues closed in 1.1
--------------------

GitHub stats for 2013/08/08 - 2013/09/09 (since 1.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

The following 25 authors contributed 337 commits.

* Benjamin Ragan-Kelley
* Bing Xia
* Bradley M. Froehle
* Brian E. Granger
* Damián Avila
* dhirschfeld
* Dražen Lučanin
* gmbecker
* Jake Vanderplas
* Jason Grout
* Jonathan Frederic
* Kevin Burke
* Kyle Kelley
* Matt Henderson
* Matthew Brett
* Matthias Bussonnier
* Pankaj Pandey
* Paul Ivanov
* rossant
* Samuel Ainsworth
* Stephan Rave
* stonebig
* Thomas Kluyver
* Yaroslav Halchenko
* Zachary Sailer


We closed a total of 76 issues, 58 pull requests and 18 regular issues;
this is the full list (generated with the script :file:`tools/github_stats.py`):

Pull Requests (58):

* :ghpull:`4188`: Allow user_ns trait to be None
* :ghpull:`4189`: always fire LOCAL_IPS.extend(PUBLIC_IPS)
* :ghpull:`4174`: various issues in markdown and rst templates
* :ghpull:`4178`: add missing data_javascript
* :ghpull:`4181`: nbconvert: Fix, sphinx template not removing new lines from headers
* :ghpull:`4043`: don't 'restore_bytes' in from_JSON
* :ghpull:`4163`: Fix for incorrect default encoding on Windows.
* :ghpull:`4136`: catch javascript errors in any output
* :ghpull:`4171`: add nbconvert config file when creating profiles
* :ghpull:`4125`: Basic exercise of `ipython [subcommand] -h` and help-all
* :ghpull:`4085`: nbconvert: Fix sphinx preprocessor date format string for Windows
* :ghpull:`4159`: don't split `.cell` and `div.cell` CSS
* :ghpull:`4158`: generate choices for `--gui` configurable from real mapping
* :ghpull:`4065`: do not include specific css in embedable one
* :ghpull:`4092`: nbconvert: Fix for unicode html headers, Windows + Python 2.x
* :ghpull:`4074`: close Client sockets if connection fails
* :ghpull:`4064`: Store default codemirror mode in only 1 place
* :ghpull:`4104`: Add way to install MathJax to a particular profile
* :ghpull:`4144`: help_end transformer shouldn't pick up ? in multiline string
* :ghpull:`4143`: update example custom.js
* :ghpull:`4142`: DOC: unwrap openssl line in public_server doc
* :ghpull:`4141`: add files with a separate `add` call in backport_pr
* :ghpull:`4137`: Restore autorestore option for storemagic
* :ghpull:`4098`: pass profile-dir instead of profile name to Kernel
* :ghpull:`4120`: support `input` in Python 2 kernels
* :ghpull:`4088`: nbconvert: Fix coalescestreams line with incorrect nesting causing strange behavior
* :ghpull:`4060`: only strip continuation prompts if regular prompts seen first
* :ghpull:`4132`: Fixed name error bug in function safe_unicode in module py3compat.
* :ghpull:`4121`: move test_kernel from IPython.zmq to IPython.kernel
* :ghpull:`4118`: ZMQ heartbeat channel: catch EINTR exceptions and continue.
* :ghpull:`4054`: use unicode for HTML export
* :ghpull:`4106`: fix a couple of default block values
* :ghpull:`4115`: Update docs on declaring a magic function
* :ghpull:`4101`: restore accidentally removed EngineError
* :ghpull:`4096`: minor docs changes
* :ghpull:`4056`: respect `pylab_import_all` when `--pylab` specified at the command-line
* :ghpull:`4091`: Make Qt console banner configurable
* :ghpull:`4086`: fix missing errno import
* :ghpull:`4030`: exclude `.git` in MANIFEST.in
* :ghpull:`4047`: Use istype() when checking if canned object is a dict
* :ghpull:`4031`: don't close_fds on Windows
* :ghpull:`4029`: bson.Binary moved
* :ghpull:`4035`: Fixed custom jinja2 templates being ignored when setting template_path
* :ghpull:`4026`: small doc fix in nbconvert
* :ghpull:`4016`: Fix IPython.start_* functions
* :ghpull:`4021`: Fix parallel.client.View map() on numpy arrays
* :ghpull:`4022`: DOC: fix links to matplotlib, notebook docs
* :ghpull:`4018`: Fix warning when running IPython.kernel tests
* :ghpull:`4019`: Test skipping without unicode paths
* :ghpull:`4008`: Transform code before %prun/%%prun runs
* :ghpull:`4014`: Fix typo in ipapp
* :ghpull:`3987`: get files list in backport_pr
* :ghpull:`3974`: nbconvert: Fix app tests on Window7 w/ Python 3.3
* :ghpull:`3978`: fix `--existing` with non-localhost IP
* :ghpull:`3939`: minor checkpoint cleanup
* :ghpull:`3981`: BF: fix nbconvert rst input prompt spacing
* :ghpull:`3960`: Don't make sphinx a dependency for importing nbconvert
* :ghpull:`3973`: logging.Formatter is not new-style in 2.6

Issues (18):

* :ghissue:`4024`: nbconvert markdown issues
* :ghissue:`4095`: Catch js error in append html in stream/pyerr
* :ghissue:`4156`: Specifying --gui=tk at the command line
* :ghissue:`3818`: nbconvert can't handle Heading with Chinese characters on Japanese Windows OS.
* :ghissue:`4134`: multi-line parser fails on ''' in comment, qtconsole and notebook.
* :ghissue:`3998`: sample custom.js needs to be updated
* :ghissue:`4078`: StoreMagic.autorestore not working in 1.0.0
* :ghissue:`3990`: Buitlin `input` doesn't work over zmq
* :ghissue:`4015`: nbconvert fails to convert all the content of a notebook
* :ghissue:`4059`: Issues with Ellipsis literal in Python 3
* :ghissue:`4103`: Wrong default argument of DirectView.clear
* :ghissue:`4100`: parallel.client.client references undefined error.EngineError
* :ghissue:`4005`: IPython.start_kernel doesn't work.
* :ghissue:`4020`: IPython parallel map fails on numpy arrays
* :ghissue:`3945`: nbconvert: commandline tests fail Win7x64 Py3.3
* :ghissue:`3977`: unable to complete remote connections for two-process 
* :ghissue:`3980`: nbconvert rst output lacks needed blank lines
* :ghissue:`3968`: TypeError: super() argument 1 must be type, not classobj (Python 2.6.6)

Issues closed in 1.0
--------------------

GitHub stats for 2012/06/30 - 2013/08/08 (since 0.13)

These lists are automatically generated, and may be incomplete or contain duplicates.

The following 155 authors contributed 4258 commits.

* Aaron Meurer
* Adam Davis
* Ahmet Bakan
* Alberto Valverde
* Allen Riddell
* Anders Hovmöller
* Andrea Bedini
* Andrew Spiers
* Andrew Vandever
* Anthony Scopatz
* Anton Akhmerov
* Anton I. Sipos
* Antony Lee
* Aron Ahmadia
* Benedikt Sauer
* Benjamin Jones
* Benjamin Ragan-Kelley
* Benjie Chen
* Boris de Laage
* Brad Reisfeld
* Bradley M. Froehle
* Brian E. Granger
* Cameron Bates
* Cavendish McKay
* chapmanb
* Chris Beaumont
* Chris Laumann
* Christoph Gohlke
* codebraker
* codespaced
* Corran Webster
* DamianHeard
* Damián Avila
* Dan Kilman
* Dan McDougall
* Danny Staple
* David Hirschfeld
* David P. Sanders
* David Warde-Farley
* David Wolever
* David Wyde
* debjan
* Diane Trout
* dkua
* Dominik Dabrowski
* Donald Curtis
* Dražen Lučanin
* drevicko
* Eric O. LEBIGOT
* Erik M. Bray
* Erik Tollerud
* Eugene Van den Bulke
* Evan Patterson
* Fernando Perez
* Francesco Montesano
* Frank Murphy
* Greg Caporaso
* Guy Haskin Fernald
* guziy
* Hans Meine
* Harry Moreno
* henryiii
* Ivan Djokic
* Jack Feser
* Jake Vanderplas
* jakobgager
* James Booth
* Jan Schulz
* Jason Grout
* Jeff Knisley
* Jens Hedegaard Nielsen
* jeremiahbuddha
* Jerry Fowler
* Jessica B. Hamrick
* Jez Ng
* John Zwinck
* Jonathan Frederic
* Jonathan Taylor
* Joon Ro
* Joseph Lansdowne
* Juergen Hasch
* Julian Taylor
* Jussi Sainio
* Jörgen Stenarson
* kevin
* klonuo
* Konrad Hinsen
* Kyle Kelley
* Lars Solberg
* Lessandro Mariano
* Mark Sienkiewicz at STScI
* Martijn Vermaat
* Martin Spacek
* Matthias Bussonnier
* Maxim Grechkin
* Maximilian Albert
* MercuryRising
* Michael Droettboom
* Michael Shuffett
* Michał Górny
* Mikhail Korobov
* mr.Shu
* Nathan Goldbaum
* ocefpaf
* Ohad Ravid
* Olivier Grisel
* Olivier Verdier
* Owen Healy
* Pankaj Pandey
* Paul Ivanov
* Pawel Jasinski
* Pietro Berkes
* Piti Ongmongkolkul
* Puneeth Chaganti
* Rich Wareham
* Richard Everson
* Rick Lupton
* Rob Young
* Robert Kern
* Robert Marchman
* Robert McGibbon
* Rui Pereira
* Rustam Safin
* Ryan May
* s8weber
* Samuel Ainsworth
* Sean Vig
* Siyu Zhang
* Skylar Saveland
* slojo404
* smithj1
* Stefan Karpinski
* Stefan van der Walt
* Steven Silvester
* Takafumi Arakaki
* Takeshi Kanmae
* tcmulcahy
* teegaar
* Thomas Kluyver
* Thomas Robitaille
* Thomas Spura
* Thomas Weißschuh
* Timothy O'Donnell
* Tom Dimiduk
* ugurthemaster
* urielshaolin
* v923z
* Valentin Haenel
* Victor Zverovich
* W. Trevor King
* y-p
* Yoav Ram
* Zbigniew Jędrzejewski-Szmek
* Zoltán Vörös


We closed a total of 1484 issues, 793 pull requests and 691 regular issues;
this is the full list (generated with the script 
:file:`tools/github_stats.py`):

Pull Requests (793):

* :ghpull:`3958`: doc update
* :ghpull:`3965`: Fix ansi color code for background yellow
* :ghpull:`3964`: Fix casing of message.
* :ghpull:`3942`: Pass on install docs
* :ghpull:`3962`: exclude IPython.lib.kernel in iptest
* :ghpull:`3961`: Longpath test fix
* :ghpull:`3905`: Remove references to 0.11 and 0.12 from config/overview.rst
* :ghpull:`3951`: nbconvert: fixed latex characters not escaped properly in nbconvert
* :ghpull:`3949`: log fatal error when PDF conversion fails
* :ghpull:`3947`: nbconvert: Make writer & post-processor aliases case insensitive.
* :ghpull:`3938`: Recompile css.
* :ghpull:`3948`: sphinx and PDF tweaks
* :ghpull:`3943`: nbconvert: Serve post-processor Windows fix
* :ghpull:`3934`: nbconvert: fix logic of verbose flag in PDF post processor
* :ghpull:`3929`: swallow enter event in rename dialog
* :ghpull:`3924`: nbconvert: Backport fixes
* :ghpull:`3925`: Replace --pylab flag with --matplotlib in usage
* :ghpull:`3910`: Added explicit error message for missing configuration arguments.
* :ghpull:`3913`: grffile to support spaces in notebook names
* :ghpull:`3918`: added check_for_tornado, closes #3916
* :ghpull:`3917`: change docs/examples refs to be just examples
* :ghpull:`3908`: what's new tweaks
* :ghpull:`3896`: two column quickhelp dialog, closes #3895
* :ghpull:`3911`: explicitly load python mode before IPython mode
* :ghpull:`3901`: don't force . relative path, fix #3897
* :ghpull:`3891`: fix #3889
* :ghpull:`3892`: Fix documentation of Kernel.stop_channels
* :ghpull:`3888`: posixify paths for Windows latex
* :ghpull:`3882`: quick fix for #3881
* :ghpull:`3877`: don't use `shell=True` in PDF export
* :ghpull:`3878`: minor template loading cleanup
* :ghpull:`3855`: nbconvert: Filter tests
* :ghpull:`3879`: finish 3870
* :ghpull:`3870`: Fix for converting notebooks that contain unicode characters.
* :ghpull:`3876`: Update parallel_winhpc.rst
* :ghpull:`3872`: removing vim-ipython, since it has it's own repo
* :ghpull:`3871`: updating docs
* :ghpull:`3873`: remove old examples
* :ghpull:`3868`: update CodeMirror component to 3.15
* :ghpull:`3865`: Escape filename for pdflatex in nbconvert
* :ghpull:`3861`: remove old external.js
* :ghpull:`3864`: add keyboard shortcut to docs
* :ghpull:`3834`: This PR fixes a few issues with nbconvert tests
* :ghpull:`3840`: prevent profile_dir from being undefined
* :ghpull:`3859`: Add "An Afternoon Hack" to docs
* :ghpull:`3854`: Catch errors filling readline history on startup
* :ghpull:`3857`: Delete extra auto
* :ghpull:`3845`: nbconvert: Serve from original build directory
* :ghpull:`3846`: Add basic logging to nbconvert
* :ghpull:`3850`: add missing store_history key to Notebook execute_requests
* :ghpull:`3844`: update payload source
* :ghpull:`3830`: mention metadata / display_data similarity in pyout spec
* :ghpull:`3848`: fix incorrect `empty-docstring`
* :ghpull:`3836`: Parse markdown correctly when mathjax is disabled
* :ghpull:`3849`: skip a failing test on windows
* :ghpull:`3828`: signature_scheme lives in Session
* :ghpull:`3831`: update nbconvert doc with new CLI
* :ghpull:`3822`: add output flag to nbconvert
* :ghpull:`3780`: Added serving the output directory if html-based format are selected.
* :ghpull:`3764`: Cleanup nbconvert templates
* :ghpull:`3829`: remove now-duplicate 'this is dev' note
* :ghpull:`3814`: add `ConsoleWidget.execute_on_complete_input` flag
* :ghpull:`3826`: try rtfd
* :ghpull:`3821`: add sphinx prolog
* :ghpull:`3817`: relax timeouts in terminal console and tests
* :ghpull:`3825`: fix more tests that fail when pandoc is missing
* :ghpull:`3824`: don't set target on internal markdown links
* :ghpull:`3816`: s/pylab/matplotlib in docs
* :ghpull:`3812`: Describe differences between start_ipython and embed
* :ghpull:`3805`: Print View has been removed
* :ghpull:`3820`: Make it clear that 1.0 is not released yet
* :ghpull:`3784`: nbconvert: Export flavors & PDF writer (ipy dev meeting)
* :ghpull:`3800`: semantic-versionify version number for non-releases
* :ghpull:`3802`: Documentation .txt to .rst
* :ghpull:`3765`: cleanup terminal console iopub handling
* :ghpull:`3720`: Fix for #3719
* :ghpull:`3787`: re-raise KeyboardInterrupt in raw_input
* :ghpull:`3770`: Organizing reveal's templates.
* :ghpull:`3751`: Use link(2) when possible in nbconvert
* :ghpull:`3792`: skip tests that require pandoc
* :ghpull:`3782`: add Importing Notebooks example
* :ghpull:`3752`: nbconvert: Add cwd to sys.path
* :ghpull:`3789`: fix raw_input in qtconsole
* :ghpull:`3756`: document the wire protocol
* :ghpull:`3749`: convert IPython syntax to Python syntax in nbconvert python template
* :ghpull:`3793`: Closes #3788
* :ghpull:`3794`: Change logo link to ipython.org
* :ghpull:`3746`: Raise a named exception when pandoc is missing
* :ghpull:`3781`: comply with the message spec in the notebook
* :ghpull:`3779`: remove bad `if logged_in` preventing new-notebook without login
* :ghpull:`3743`: remove notebook read-only view
* :ghpull:`3732`: add delay to autosave in beforeunload
* :ghpull:`3761`: Added rm_math_space to markdown cells in the basichtml.tpl to be rendered ok by mathjax after the nbconvertion.
* :ghpull:`3758`: nbconvert: Filter names cleanup
* :ghpull:`3769`: Add configurability to  tabcompletion timeout
* :ghpull:`3771`: Update px pylab test to match new output of pylab
* :ghpull:`3741`: better message when notebook format is not supported
* :ghpull:`3753`: document Ctrl-C not working in ipython kernel
* :ghpull:`3766`: handle empty metadata in pyout messages more gracefully.
* :ghpull:`3736`: my attempt to fix #3735
* :ghpull:`3759`: nbconvert: Provide a more useful error for invalid use case.
* :ghpull:`3760`: nbconvert: Allow notebook filenames without their extensions
* :ghpull:`3750`: nbconvert: Add cwd to default templates search path.
* :ghpull:`3748`: Update nbconvert docs
* :ghpull:`3734`: Nbconvert: Export extracted files into `nbname_files` subdirectory
* :ghpull:`3733`: Nicer message when pandoc is missing, closes #3730
* :ghpull:`3722`: fix two failing test in IPython.lib
* :ghpull:`3704`: Start what's new for 1.0
* :ghpull:`3705`: Complete rewrite of IPython Notebook documentation: docs/source/interactive/htmlnotebook.txt
* :ghpull:`3709`: Docs cleanup
* :ghpull:`3716`: raw_input fixes for kernel restarts
* :ghpull:`3683`: use `%matplotlib` in example notebooks
* :ghpull:`3686`: remove quarantine
* :ghpull:`3699`: svg2pdf unicode fix
* :ghpull:`3695`: fix SVG2PDF
* :ghpull:`3685`: fix Pager.detach
* :ghpull:`3675`: document new dependencies
* :ghpull:`3690`: Fixing some css minors in full_html and reveal.
* :ghpull:`3671`: nbconvert tests
* :ghpull:`3692`: Fix rename notebook - show error with invalid name
* :ghpull:`3409`: Prevent qtconsole frontend freeze on lots of output.
* :ghpull:`3660`: refocus active cell on dialog close
* :ghpull:`3598`: Statelessify mathjaxutils
* :ghpull:`3673`: enable comment/uncomment selection
* :ghpull:`3677`: remove special-case in get_home_dir for frozen dists
* :ghpull:`3674`: add CONTRIBUTING.md
* :ghpull:`3670`: use Popen command list for ipexec
* :ghpull:`3568`: pylab import adjustments
* :ghpull:`3559`: add create.Cell and delete.Cell js events
* :ghpull:`3606`: push cell magic to the head of the transformer line
* :ghpull:`3607`: NbConvert: Writers, No YAML, and stuff...
* :ghpull:`3665`: Pywin32 skips
* :ghpull:`3669`: set default client_class for QtKernelManager
* :ghpull:`3662`: add strip_encoding_cookie transformer
* :ghpull:`3641`: increase patience for slow kernel startup in tests
* :ghpull:`3651`: remove a bunch of unused `default_config_file` assignments
* :ghpull:`3630`: CSS adjustments
* :ghpull:`3645`: Don't require HistoryManager to have a shell
* :ghpull:`3643`: don't assume tested ipython is on the PATH
* :ghpull:`3654`: fix single-result AsyncResults
* :ghpull:`3601`: Markdown in heading cells (take 2)
* :ghpull:`3652`: Remove old `docs/examples`
* :ghpull:`3621`: catch any exception appending output
* :ghpull:`3585`: don't blacklist builtin names
* :ghpull:`3647`: Fix `frontend` deprecation warnings in several examples
* :ghpull:`3649`: fix AsyncResult.get_dict for single result
* :ghpull:`3648`: Fix store magic test 
* :ghpull:`3650`: Fix, config_file_name was ignored
* :ghpull:`3640`: Gcf.get_active() can return None
* :ghpull:`3571`: Added shorcuts to split cell, merge cell above and merge cell below.
* :ghpull:`3635`: Added missing slash to print-pdf call.
* :ghpull:`3487`: Drop patch for compatibility with pyreadline 1.5
* :ghpull:`3338`: Allow filename with extension in find_cmd in Windows.
* :ghpull:`3628`: Fix test for Python 3 on Windows.
* :ghpull:`3642`: Fix typo in docs
* :ghpull:`3627`: use DEFAULT_STATIC_FILES_PATH in a test instead of package dir
* :ghpull:`3624`: fix some unicode in zmqhandlers
* :ghpull:`3460`: Set calling program to UNKNOWN, when argv not in sys
* :ghpull:`3632`: Set calling program to UNKNOWN, when argv not in sys (take #2)
* :ghpull:`3629`: Use new entry point for python -m IPython
* :ghpull:`3626`: passing cell to showInPager, closes #3625
* :ghpull:`3618`: expand terminal color support
* :ghpull:`3623`: raise UsageError for unsupported GUI backends
* :ghpull:`3071`: Add magic function %drun to run code in debugger
* :ghpull:`3608`: a nicer error message when using %pylab magic
* :ghpull:`3592`: add extra_config_file
* :ghpull:`3612`: updated .mailmap
* :ghpull:`3616`: Add examples for interactive use of MPI.
* :ghpull:`3615`: fix regular expression for ANSI escapes
* :ghpull:`3586`: Corrected a typo in the format string for strftime the sphinx.py transformer of nbconvert
* :ghpull:`3611`: check for markdown no longer needed, closes #3610
* :ghpull:`3555`: Simplify caching of modules with %run
* :ghpull:`3583`: notebook small things
* :ghpull:`3594`: Fix duplicate completion in notebook
* :ghpull:`3600`: parallel: Improved logging for errors during BatchSystemLauncher.stop
* :ghpull:`3595`: Revert "allow markdown in heading cells"
* :ghpull:`3538`: add IPython.start_ipython
* :ghpull:`3562`: Allow custom nbconvert template loaders
* :ghpull:`3582`: pandoc adjustments
* :ghpull:`3560`: Remove max_msg_size
* :ghpull:`3591`: Refer to Setuptools instead of Distribute
* :ghpull:`3590`: IPython.sphinxext needs an __init__.py
* :ghpull:`3581`: Added the possibility to read a custom.css file for tweaking the final html in full_html and reveal templates.
* :ghpull:`3576`: Added support for markdown in heading cells when they are nbconverted.
* :ghpull:`3575`: tweak `run -d` message to 'continue execution'
* :ghpull:`3569`: add PYTHONSTARTUP to startup files
* :ghpull:`3567`: Trigger a single event on js app initialized
* :ghpull:`3565`: style.min.css should always exist...
* :ghpull:`3531`: allow markdown in heading cells
* :ghpull:`3577`: Simplify codemirror ipython-mode
* :ghpull:`3495`: Simplified regexp, and suggestions for clearer regexps.
* :ghpull:`3578`: Use adjustbox to specify figure size in nbconvert -> latex
* :ghpull:`3572`: Skip import irunner test on Windows.
* :ghpull:`3574`: correct static path for CM modes autoload
* :ghpull:`3558`: Add IPython.sphinxext
* :ghpull:`3561`: mention double-control-C to stop notebook server
* :ghpull:`3566`: fix event names
* :ghpull:`3564`: Remove trivial nbconvert example
* :ghpull:`3540`: allow cython cache dir to be deleted
* :ghpull:`3527`: cleanup stale, unused exceptions in parallel.error
* :ghpull:`3529`: ensure raw_input returns str in zmq shell
* :ghpull:`3541`: respect image size metadata in qtconsole
* :ghpull:`3550`: Fixing issue preventing the correct read of images by full_html and reveal exporters.
* :ghpull:`3557`: open markdown links in new tabs
* :ghpull:`3556`: remove mention of nonexistent `_margv` in macro
* :ghpull:`3552`: set overflow-x: hidden on Firefox only
* :ghpull:`3554`: Fix missing import os in latex exporter.
* :ghpull:`3546`: Don't hardcode **latex** posix paths in nbconvert
* :ghpull:`3551`: fix path prefix in nbconvert
* :ghpull:`3533`: Use a CDN to get reveal.js library.
* :ghpull:`3498`: When a notebook is written to file, name the metadata name u''.
* :ghpull:`3548`: Change to standard save icon in Notebook toolbar
* :ghpull:`3539`: Don't hardcode posix paths in nbconvert
* :ghpull:`3508`: notebook supports raw_input and %debug now
* :ghpull:`3526`: ensure 'default' is first in cluster profile list
* :ghpull:`3525`: basic timezone info
* :ghpull:`3532`: include nbconvert templates in installation
* :ghpull:`3515`: update CodeMirror component to 3.14
* :ghpull:`3513`: add 'No Checkpoints' to Revert menu
* :ghpull:`3536`: format positions are required in Python 2.6.x
* :ghpull:`3521`: Nbconvert fix, silent fail if template doesn't exist
* :ghpull:`3530`: update %store magic docstring
* :ghpull:`3528`: fix local mathjax with custom base_project_url
* :ghpull:`3518`: Clear up unused imports
* :ghpull:`3506`: %store -r restores saved aliases and directory history, as well as variables
* :ghpull:`3516`: make css highlight style configurable
* :ghpull:`3523`: Exclude frontend shim from docs build
* :ghpull:`3514`: use bootstrap `disabled` instead of `ui-state-disabled`
* :ghpull:`3520`: Added relative import of RevealExporter to __init__.py inside exporters module
* :ghpull:`3507`: fix HTML capitalization in nbconvert exporter classes
* :ghpull:`3512`: fix nbconvert filter validation
* :ghpull:`3511`: Get Tracer working after ipapi.get replaced with get_ipython
* :ghpull:`3510`: use `window.onbeforeunload=` for nav-away warning
* :ghpull:`3504`: don't use parent=self in handlers
* :ghpull:`3500`: Merge nbconvert into IPython
* :ghpull:`3478`: restore "unsaved changes" warning on unload
* :ghpull:`3493`: add a dialog when the kernel is auto-restarted
* :ghpull:`3488`: Add test suite for autoreload extension
* :ghpull:`3484`: Catch some pathological cases inside oinspect
* :ghpull:`3481`: Display R errors without Python traceback
* :ghpull:`3468`: fix `%magic` output
* :ghpull:`3430`: add parent to Configurable
* :ghpull:`3491`: Remove unexpected keyword parameter to remove_kernel
* :ghpull:`3485`: SymPy has changed its recommended way to initialize printing
* :ghpull:`3486`: Add test for non-ascii characters in docstrings
* :ghpull:`3483`: Inputtransformer: Allow classic prompts without space
* :ghpull:`3482`: Use an absolute path to iptest, because the tests are not always run from $IPYTHONDIR.
* :ghpull:`3381`: enable 2x (retina) display
* :ghpull:`3450`: Flatten IPython.frontend
* :ghpull:`3477`: pass config to subapps
* :ghpull:`3466`: Kernel fails to start when username has non-ascii characters
* :ghpull:`3465`: Add HTCondor bindings to IPython.parallel
* :ghpull:`3463`: fix typo, closes #3462
* :ghpull:`3456`: Notice for users who disable javascript
* :ghpull:`3453`: fix cell execution in firefox, closes #3447
* :ghpull:`3393`: [WIP] bootstrapify
* :ghpull:`3440`: Fix installing mathjax from downloaded file via command line
* :ghpull:`3431`: Provide means for starting the Qt console maximized and with the menu bar hidden
* :ghpull:`3425`: base IPClusterApp inherits from BaseIPythonApp
* :ghpull:`3433`: Update IPython\external\path\__init__.py
* :ghpull:`3298`: Some fixes in IPython Sphinx directive
* :ghpull:`3428`: process escapes in mathjax
* :ghpull:`3420`: thansk -> thanks
* :ghpull:`3416`: Fix doc: "principle" not "principal"
* :ghpull:`3413`: more unique filename for test
* :ghpull:`3364`: Inject requirejs in notebook and start using it.
* :ghpull:`3390`: Fix %paste with blank lines
* :ghpull:`3403`: fix creating config objects from dicts
* :ghpull:`3401`: rollback #3358
* :ghpull:`3373`: make cookie_secret configurable
* :ghpull:`3307`: switch default ws_url logic to js side
* :ghpull:`3392`: Restore anchor link on h2-h6
* :ghpull:`3369`: Use different threshold for (auto)scroll in output
* :ghpull:`3370`: normalize unicode notebook filenames
* :ghpull:`3372`: base default cookie name on request host+port
* :ghpull:`3378`: disable CodeMirror drag/drop on Safari
* :ghpull:`3358`: workaround spurious CodeMirror scrollbars
* :ghpull:`3371`: make setting the notebook dirty flag an event
* :ghpull:`3366`: remove long-dead zmq frontend.py and completer.py
* :ghpull:`3382`: cull Session digest history
* :ghpull:`3330`: Fix get_ipython_dir when $HOME is /
* :ghpull:`3319`: IPEP 13: user-expressions and user-variables
* :ghpull:`3384`: comments in tools/gitwash_dumper.py changed (''' to """)
* :ghpull:`3387`: Make submodule checks work under Python 3.
* :ghpull:`3357`: move anchor-link off of heading text
* :ghpull:`3351`: start basic tests of ipcluster Launchers
* :ghpull:`3377`: allow class.__module__ to be None
* :ghpull:`3340`: skip submodule check in package managers
* :ghpull:`3328`: decode subprocess output in launchers
* :ghpull:`3368`: Reenable bracket matching
* :ghpull:`3356`: Mpr fixes
* :ghpull:`3336`: Use new input transformation API in %time magic
* :ghpull:`3325`: Organize the JS and less files by component.
* :ghpull:`3342`: fix test_find_cmd_python
* :ghpull:`3354`: catch socket.error in utils.localinterfaces
* :ghpull:`3341`: fix default cluster count
* :ghpull:`3286`: don't use `get_ipython` from builtins in library code
* :ghpull:`3333`: notebookapp: add missing whitespace to warnings
* :ghpull:`3323`: Strip prompts even if the prompt isn't present on the first line.
* :ghpull:`3321`: Reorganize the python/server side of the notebook
* :ghpull:`3320`: define `__file__` in config files
* :ghpull:`3317`: rename `%%file` to `%%writefile`
* :ghpull:`3304`: set unlimited HWM for all relay devices
* :ghpull:`3315`: Update Sympy_printing extension load
* :ghpull:`3310`: further clarify Image docstring
* :ghpull:`3285`: load extensions in builtin trap
* :ghpull:`3308`: Speed up AsyncResult._wait_for_outputs(0)
* :ghpull:`3294`: fix callbacks as optional in js kernel.execute
* :ghpull:`3276`: Fix: "python ABS/PATH/TO/ipython.py" fails
* :ghpull:`3301`: allow python3 tests without python installed
* :ghpull:`3282`: allow view.map to work with a few more things
* :ghpull:`3284`: remove `ipython.py` entry point
* :ghpull:`3281`: fix ignored IOPub messages with no parent
* :ghpull:`3275`: improve submodule messages / git hooks
* :ghpull:`3239`: Allow "x" icon and esc key to close pager in notebook
* :ghpull:`3290`: Improved heartbeat controller to engine monitoring for long running tasks
* :ghpull:`3142`: Better error message when CWD doesn't exist on startup
* :ghpull:`3066`: Add support for relative import to %run -m (fixes #2727)
* :ghpull:`3269`: protect highlight.js against unknown languages
* :ghpull:`3267`: add missing return
* :ghpull:`3101`: use marked / highlight.js instead of pagedown and prettify
* :ghpull:`3264`: use https url for submodule
* :ghpull:`3263`: fix set_last_checkpoint when no checkpoint
* :ghpull:`3258`: Fix submodule location in setup.py
* :ghpull:`3254`: fix a few URLs from previous PR
* :ghpull:`3240`: remove js components from the repo
* :ghpull:`3158`: IPEP 15: autosave the notebook
* :ghpull:`3252`: move images out of _static folder into _images
* :ghpull:`3251`: Fix for cell magics in Qt console
* :ghpull:`3250`: Added a simple __html__() method to the HTML class
* :ghpull:`3249`: remove copy of sphinx inheritance_diagram.py
* :ghpull:`3235`: Remove the unused print notebook view
* :ghpull:`3238`: Improve the design of the tab completion UI
* :ghpull:`3242`: Make changes of Application.log_format effective
* :ghpull:`3219`: Workaround so only one CTRL-C is required for a new prompt in --gui=qt
* :ghpull:`3190`: allow formatters to specify metadata
* :ghpull:`3231`: improve discovery of public IPs
* :ghpull:`3233`: check prefixes for swallowing kernel args
* :ghpull:`3234`: Removing old autogrow JS code.
* :ghpull:`3232`: Update to CodeMirror 3 and start to ship our components
* :ghpull:`3229`: The HTML output type accidentally got removed from the OutputArea.
* :ghpull:`3228`: Typo in IPython.Parallel documentation
* :ghpull:`3226`: Text in rename dialog was way too big - making it <p>.
* :ghpull:`3225`: Removing old restuctured text handler and web service.
* :ghpull:`3222`: make BlockingKernelClient the default Client
* :ghpull:`3223`: add missing mathjax_url to new settings dict
* :ghpull:`3089`: add stdin to the notebook
* :ghpull:`3221`: Remove references to HTMLCell (dead code)
* :ghpull:`3205`: add ignored ``*args`` to HasTraits constructor
* :ghpull:`3088`: cleanup IPython handler settings
* :ghpull:`3201`: use much faster regexp for ansi coloring
* :ghpull:`3220`: avoid race condition in profile creation
* :ghpull:`3011`: IPEP 12: add KernelClient
* :ghpull:`3217`: informative error when trying to load directories
* :ghpull:`3174`: Simple class
* :ghpull:`2979`: CM configurable Take 2
* :ghpull:`3215`: Updates storemagic extension to allow for specifying variable name to load
* :ghpull:`3181`: backport If-Modified-Since fix from tornado
* :ghpull:`3200`: IFrame (VimeoVideo, ScribdDocument, ...) 
* :ghpull:`3186`: Fix small inconsistency in nbconvert: etype -> ename
* :ghpull:`3212`: Fix issue #2563, "core.profiledir.check_startup_dir() doesn't work inside py2exe'd installation"
* :ghpull:`3211`: Fix inheritance_diagram Sphinx extension for Sphinx 1.2
* :ghpull:`3208`: Update link to extensions index
* :ghpull:`3203`: Separate InputSplitter for transforming whole cells
* :ghpull:`3189`: Improve completer
* :ghpull:`3194`: finish up PR #3116
* :ghpull:`3188`: Add new keycodes
* :ghpull:`2695`: Key the root modules cache by sys.path entries.
* :ghpull:`3182`: clarify %%file docstring
* :ghpull:`3163`: BUG: Fix the set and frozenset pretty printer to handle the empty case correctly
* :ghpull:`3180`: better UsageError for cell magic with no body
* :ghpull:`3184`: Cython cache
* :ghpull:`3175`: Added missing s
* :ghpull:`3173`: Little bits of documentation cleanup
* :ghpull:`2635`: Improve Windows start menu shortcuts (#2)
* :ghpull:`3172`: Add missing import in IPython parallel magics example
* :ghpull:`3170`: default application logger shouldn't propagate
* :ghpull:`3159`: Autocompletion for zsh
* :ghpull:`3105`: move DEFAULT_STATIC_FILES_PATH to IPython.html
* :ghpull:`3144`: minor bower tweaks
* :ghpull:`3141`: Default color output for ls on OSX
* :ghpull:`3137`: fix dot syntax error in inheritance diagram
* :ghpull:`3072`: raise UnsupportedOperation on iostream.fileno()
* :ghpull:`3147`: Notebook support for a reverse proxy which handles SSL
* :ghpull:`3152`: make qtconsole size at startup configurable
* :ghpull:`3162`: adding stream kwarg to current.new_output
* :ghpull:`2981`: IPEP 10: kernel side filtering of display formats
* :ghpull:`3058`: add redirect handler for notebooks by name
* :ghpull:`3041`: support non-modules in @require
* :ghpull:`2447`: Stateful line transformers
* :ghpull:`3108`: fix some O(N) and O(N^2) operations in parallel.map
* :ghpull:`2791`: forward stdout from forked processes
* :ghpull:`3157`: use Python 3-style for pretty-printed sets
* :ghpull:`3148`: closes #3045, #3123 for tornado < version 3.0
* :ghpull:`3143`: minor heading-link tweaks
* :ghpull:`3136`: Strip useless ANSI escape codes in notebook
* :ghpull:`3126`: Prevent errors when pressing arrow keys in an empty notebook
* :ghpull:`3135`: quick dev installation instructions
* :ghpull:`2889`: Push pandas dataframes to R magic
* :ghpull:`3068`: Don't monkeypatch doctest during IPython startup.
* :ghpull:`3133`: fix argparse version check
* :ghpull:`3102`: set `spellcheck=false` in CodeCell inputarea
* :ghpull:`3064`: add anchors to heading cells
* :ghpull:`3097`: PyQt 4.10: use self._document = self.document()
* :ghpull:`3117`: propagate automagic change to shell
* :ghpull:`3118`: don't give up on weird os names
* :ghpull:`3115`: Fix example
* :ghpull:`2640`: fix quarantine/ipy_editors.py
* :ghpull:`3070`: Add info make target that was missing in old Sphinx
* :ghpull:`3082`: A few small patches to image handling
* :ghpull:`3078`: fix regular expression for detecting links in stdout
* :ghpull:`3054`: restore default behavior for automatic cluster size
* :ghpull:`3073`: fix ipython usage text
* :ghpull:`3083`: fix DisplayMagics.html docstring
* :ghpull:`3080`: noted sub_channel being renamed to iopub_channel
* :ghpull:`3079`: actually use IPKernelApp.kernel_class
* :ghpull:`3076`: Improve notebook.js documentation
* :ghpull:`3063`: add missing `%%html` magic
* :ghpull:`3075`: check for SIGUSR1 before using it, closes #3074
* :ghpull:`3051`: add width:100% to vbox for webkit / FF consistency
* :ghpull:`2999`: increase registration timeout
* :ghpull:`2997`: fix DictDB default size limit
* :ghpull:`3033`: on resume, print server info again
* :ghpull:`3062`: test double pyximport
* :ghpull:`3046`: cast kernel cwd to bytes on Python 2 on Windows
* :ghpull:`3038`: remove xml from notebook magic docstrings
* :ghpull:`3032`: fix time format to international time format
* :ghpull:`3022`: Fix test for Windows
* :ghpull:`3024`: changed instances of 'outout' to 'output' in alt texts
* :ghpull:`3013`: py3 workaround for reload in cythonmagic
* :ghpull:`2961`: time magic: shorten unnecessary output on windows
* :ghpull:`2987`: fix local files examples in markdown
* :ghpull:`2998`: fix css in .output_area pre
* :ghpull:`3003`: add $include /etc/inputrc to suggested ~/.inputrc
* :ghpull:`2957`: Refactor qt import logic. Fixes #2955
* :ghpull:`2994`: expanduser on %%file targets
* :ghpull:`2983`: fix run-all (that-> this)
* :ghpull:`2964`: fix count when testing composite error output
* :ghpull:`2967`: shows entire session history when only startsess is given
* :ghpull:`2942`: Move CM IPython theme out of codemirror folder
* :ghpull:`2929`: Cleanup cell insertion
* :ghpull:`2933`: Minordocupdate
* :ghpull:`2968`: fix notebook deletion.
* :ghpull:`2966`: Added assert msg to extract_hist_ranges()
* :ghpull:`2959`: Add command to trim the history database.
* :ghpull:`2681`: Don't enable pylab mode, when matplotlib is not importable
* :ghpull:`2901`: Fix inputhook_wx on osx
* :ghpull:`2871`: truncate potentially long CompositeErrors
* :ghpull:`2951`: use istype on lists/tuples
* :ghpull:`2946`: fix qtconsole history logic for end-of-line
* :ghpull:`2954`: fix logic for append_javascript
* :ghpull:`2941`: fix baseUrl
* :ghpull:`2903`: Specify toggle value on cell line number
* :ghpull:`2911`: display order in output area configurable
* :ghpull:`2897`: Don't rely on BaseProjectUrl data in body tag
* :ghpull:`2894`: Cm configurable
* :ghpull:`2927`: next release will be 1.0
* :ghpull:`2932`: Simplify using notebook static files from external code
* :ghpull:`2915`: added small config section to notebook docs page
* :ghpull:`2924`: safe_run_module: Silence SystemExit codes 0 and None.
* :ghpull:`2906`: Unpatch/Monkey patch CM
* :ghpull:`2921`: add menu item for undo delete cell
* :ghpull:`2917`: Don't add logging handler if one already exists.
* :ghpull:`2910`: Respect DB_IP and DB_PORT in mongodb tests
* :ghpull:`2926`: Don't die if stderr/stdout do not support set_parent() #2925
* :ghpull:`2885`: get monospace pager back
* :ghpull:`2876`: fix celltoolbar layout on FF
* :ghpull:`2904`: Skip remaining IPC test on Windows
* :ghpull:`2908`: fix last remaining KernelApp reference
* :ghpull:`2905`: fix a few remaining KernelApp/IPKernelApp changes
* :ghpull:`2900`: Don't assume test case for %time will finish in 0 time
* :ghpull:`2893`: exclude fabfile from tests
* :ghpull:`2884`: Correct import for kernelmanager on Windows
* :ghpull:`2882`: Utils cleanup
* :ghpull:`2883`: Don't call ast.fix_missing_locations unless the AST could have been modified
* :ghpull:`2855`: time(it) magic: Implement minutes/hour formatting and "%%time" cell magic
* :ghpull:`2874`: Empty cell warnings
* :ghpull:`2819`: tweak history prefix search (up/^p) in qtconsole
* :ghpull:`2868`: Import performance
* :ghpull:`2877`: minor css fixes
* :ghpull:`2880`: update examples docs with kernel move
* :ghpull:`2878`: Pass host environment on to kernel
* :ghpull:`2599`: func_kw_complete for builtin and cython with embededsignature=True using docstring
* :ghpull:`2792`: Add key "unique" to history_request protocol
* :ghpull:`2872`: fix payload keys
* :ghpull:`2869`: Fixing styling of toolbar selects on FF.
* :ghpull:`2708`: Less css
* :ghpull:`2854`: Move kernel code into IPython.kernel
* :ghpull:`2864`: Fix %run -t -N<N> TypeError
* :ghpull:`2852`: future pyzmq compatibility
* :ghpull:`2863`: whatsnew/version0.9.txt: Fix '~./ipython' -> '~/.ipython' typo
* :ghpull:`2861`: add missing KernelManager to ConsoleApp class list
* :ghpull:`2850`: Consolidate host IP detection in utils.localinterfaces
* :ghpull:`2859`: Correct docstring of ipython.py
* :ghpull:`2831`: avoid string version comparisons in external.qt
* :ghpull:`2844`: this should address the failure in #2732
* :ghpull:`2849`: utils/data: Use list comprehension for uniq_stable()
* :ghpull:`2839`: add jinja to install docs / setup.py
* :ghpull:`2841`: Miscellaneous docs fixes
* :ghpull:`2811`: Still more KernelManager cleanup
* :ghpull:`2820`: add '=' to greedy completer delims
* :ghpull:`2818`: log user tracebacks in the kernel (INFO-level)
* :ghpull:`2828`: Clean up notebook Javascript
* :ghpull:`2829`: avoid comparison error in dictdb hub history
* :ghpull:`2830`: BUG: Opening parenthesis after non-callable raises ValueError
* :ghpull:`2718`: try to fallback to pysqlite2.dbapi2 as sqlite3 in core.history
* :ghpull:`2816`: in %edit, don't save "last_call" unless last call succeeded
* :ghpull:`2817`: change ol format order
* :ghpull:`2537`: Organize example notebooks
* :ghpull:`2815`: update release/authors
* :ghpull:`2808`: improve patience for slow Hub in client tests
* :ghpull:`2812`: remove nonfunctional `-la` short arg in cython magic
* :ghpull:`2810`: remove dead utils.upgradedir
* :ghpull:`1671`: __future__ environments
* :ghpull:`2804`: skip ipc tests on Windows
* :ghpull:`2789`: Fixing styling issues with CellToolbar.
* :ghpull:`2805`: fix KeyError creating ZMQStreams in notebook
* :ghpull:`2775`: General cleanup of kernel manager code.
* :ghpull:`2340`: Initial Code to reduce parallel.Client caching
* :ghpull:`2799`: Exit code
* :ghpull:`2800`: use `type(obj) is cls` as switch when canning
* :ghpull:`2801`: Fix a breakpoint bug
* :ghpull:`2795`: Remove outdated code from extensions.autoreload
* :ghpull:`2796`: P3K: fix cookie parsing under Python 3.x (+ duplicate import is removed)
* :ghpull:`2724`: In-process kernel support (take 3)
* :ghpull:`2687`: [WIP] Metaui slideshow
* :ghpull:`2788`: Chrome frame awareness
* :ghpull:`2649`: Add version_request/reply messaging protocol
* :ghpull:`2753`: add `%%px --local` for local execution
* :ghpull:`2783`: Prefilter shouldn't touch execution_count
* :ghpull:`2333`: UI For Metadata
* :ghpull:`2396`: create a ipynbv3 json schema and a validator
* :ghpull:`2757`: check for complete pyside presence before trying to import
* :ghpull:`2782`: Allow the %run magic with '-b' to specify a file.
* :ghpull:`2778`: P3K: fix DeprecationWarning under Python 3.x 
* :ghpull:`2776`: remove non-functional View.kill method
* :ghpull:`2755`: can interactively defined classes
* :ghpull:`2774`: Removing unused code in the notebook MappingKernelManager.
* :ghpull:`2773`: Fixed minor typo causing AttributeError to be thrown.
* :ghpull:`2609`: Add 'unique' option to history_request messaging protocol
* :ghpull:`2769`: Allow shutdown when no engines are registered
* :ghpull:`2766`: Define __file__ when we %edit a real file.
* :ghpull:`2476`: allow %edit <variable> to work when interactively defined
* :ghpull:`2763`: Reset readline delimiters after loading rmagic.
* :ghpull:`2460`: Better handling of `__file__` when running scripts.
* :ghpull:`2617`: Fix for `units` argument. Adds a `res` argument.
* :ghpull:`2738`: Unicode content crashes the pager (console)
* :ghpull:`2749`: Tell Travis CI to test on Python 3.3 as well
* :ghpull:`2744`: Don't show 'try %paste' message while using magics
* :ghpull:`2728`: shift tab for tooltip
* :ghpull:`2741`: Add note to `%cython` Black-Scholes example warning of missing erf.
* :ghpull:`2743`: BUG: Octavemagic inline plots not working on Windows: Fixed
* :ghpull:`2740`: Following #2737 this error is now a name error
* :ghpull:`2737`: Rmagic: error message when moving an non-existant variable from python to R
* :ghpull:`2723`: diverse fixes for project url
* :ghpull:`2731`: %Rpush: Look for variables in the local scope first.
* :ghpull:`2544`: Infinite loop when multiple debuggers have been attached.
* :ghpull:`2726`: Add qthelp docs creation
* :ghpull:`2730`: added blockquote CSS
* :ghpull:`2729`: Fix Read the doc build, Again
* :ghpull:`2446`: [alternate 2267] Offline mathjax
* :ghpull:`2716`: remove unexisting headings level
* :ghpull:`2717`: One liner to fix debugger printing stack traces when lines of context are larger than source.
* :ghpull:`2713`: Doc bugfix: user_ns is not an attribute of Magic objects.
* :ghpull:`2690`: Fix 'import '... completion for py3 & egg files.
* :ghpull:`2691`: Document OpenMP in %%cython magic
* :ghpull:`2699`: fix jinja2 rendering for password protected notebooks
* :ghpull:`2700`: Skip notebook testing if jinja2 is not available.
* :ghpull:`2692`: Add %%cython magics to generated documentation.
* :ghpull:`2685`: Fix pretty print of types when `__module__` is not available.
* :ghpull:`2686`: Fix tox.ini
* :ghpull:`2604`: Backslashes are misinterpreted as escape-sequences by the R-interpreter.
* :ghpull:`2689`: fix error in doc (arg->kwarg) and pep-8
* :ghpull:`2683`: for downloads, replaced window.open with window.location.assign
* :ghpull:`2659`: small bugs in js are fixed
* :ghpull:`2363`: Refactor notebook templates to use Jinja2
* :ghpull:`2662`: qtconsole: wrap argument list in tooltip to match width of text body
* :ghpull:`2328`: addition of classes to generate a link or list of links from files local to the IPython HTML notebook
* :ghpull:`2668`: pylab_not_importable: Catch all exceptions, not just RuntimeErrors.
* :ghpull:`2663`: Fix issue #2660: parsing of help and version arguments
* :ghpull:`2656`: Fix irunner tests when $PYTHONSTARTUP is set
* :ghpull:`2312`: Add bracket matching to code cells in notebook
* :ghpull:`2571`: Start to document Javascript
* :ghpull:`2641`: undefinied that -> this
* :ghpull:`2638`: Fix %paste in Python 3 on Mac
* :ghpull:`2301`: Ast transfomers
* :ghpull:`2616`: Revamp API docs
* :ghpull:`2572`: Make 'Paste Above' the default paste behavior.
* :ghpull:`2574`: Fix #2244
* :ghpull:`2582`: Fix displaying history when output cache is disabled.
* :ghpull:`2591`: Fix for Issue #2584 
* :ghpull:`2526`: Don't kill paramiko tunnels when receiving ^C
* :ghpull:`2559`: Add psource, pfile, pinfo2 commands to ipdb.
* :ghpull:`2546`: use 4 Pythons to build 4 Windows installers
* :ghpull:`2561`: Fix display of plain text containing multiple carriage returns before line feed
* :ghpull:`2549`: Add a simple 'undo' for cell deletion.
* :ghpull:`2525`: Add event to kernel execution/shell reply.
* :ghpull:`2554`: Avoid stopping in ipdb until we reach the main script.
* :ghpull:`2404`: Option to limit search result in history magic command
* :ghpull:`2294`: inputhook_qt4: Use QEventLoop instead of starting up the QCoreApplication
* :ghpull:`2233`: Refactored Drag and Drop Support in Qt Console
* :ghpull:`1747`: switch between hsplit and vsplit paging (request for feedback)
* :ghpull:`2530`: Adding time offsets to the video
* :ghpull:`2542`: Allow starting IPython as `python -m IPython`.
* :ghpull:`2534`: Do not unescape backslashes in Windows (shellglob)
* :ghpull:`2517`: Improved MathJax, bug fixes
* :ghpull:`2511`: trigger default remote_profile_dir when profile_dir is set
* :ghpull:`2491`: color is supported in ironpython
* :ghpull:`2462`: Track which extensions are loaded
* :ghpull:`2464`: Locate URLs in text output and convert them to hyperlinks.
* :ghpull:`2490`: add ZMQInteractiveShell to IPEngineApp class list
* :ghpull:`2498`: Don't catch tab press when something selected
* :ghpull:`2527`: Run All Above and Run All Below
* :ghpull:`2513`: add GitHub uploads to release script
* :ghpull:`2529`: Windows aware tests for shellglob
* :ghpull:`2478`: Fix doctest_run_option_parser for Windows
* :ghpull:`2519`: clear In[ ] prompt numbers again
* :ghpull:`2467`: Clickable links
* :ghpull:`2500`: Add `encoding` attribute to `OutStream` class.
* :ghpull:`2349`: ENH: added StackExchange-style MathJax filtering
* :ghpull:`2503`: Fix traceback handling of SyntaxErrors without line numbers.
* :ghpull:`2492`: add missing 'qtconsole' extras_require
* :ghpull:`2480`: Add deprecation warnings for sympyprinting
* :ghpull:`2334`: Make the ipengine monitor the ipcontroller heartbeat and die if the ipcontroller goes down
* :ghpull:`2479`: use new _winapi instead of removed _subprocess
* :ghpull:`2474`: fix bootstrap name conflicts
* :ghpull:`2469`: Treat __init__.pyc same as __init__.py in module_list
* :ghpull:`2165`: Add -g option to %run to glob expand arguments
* :ghpull:`2468`: Tell git to ignore __pycache__ directories.
* :ghpull:`2421`: Some notebook tweaks.
* :ghpull:`2291`: Remove old plugin system
* :ghpull:`2127`: Ability to build toolbar in JS 
* :ghpull:`2445`: changes for ironpython
* :ghpull:`2420`: Pass ipython_dir to __init__() method of TerminalInteractiveShell's superclass.
* :ghpull:`2432`: Revert #1831, the `__file__` injection in safe_execfile / safe_execfile_ipy.
* :ghpull:`2216`: Autochange highlight with cell magics
* :ghpull:`1946`: Add image message handler in ZMQTerminalInteractiveShell
* :ghpull:`2424`: skip find_cmd when setting up script magics
* :ghpull:`2389`: Catch sqlite DatabaseErrors in more places when reading the history database
* :ghpull:`2395`: Don't catch ImportError when trying to unpack module functions
* :ghpull:`1868`: enable IPC transport for kernels
* :ghpull:`2437`: don't let log cleanup prevent engine start
* :ghpull:`2441`: `sys.maxsize` is the maximum length of a container.
* :ghpull:`2442`: allow iptest to be interrupted
* :ghpull:`2240`: fix message built for engine dying during task
* :ghpull:`2369`: Block until kernel termination after sending a kill signal
* :ghpull:`2439`: Py3k: Octal (0777 -> 0o777)
* :ghpull:`2326`: Detachable pager in notebook.
* :ghpull:`2377`: Fix installation of man pages in Python 3
* :ghpull:`2407`: add IPython version to message headers
* :ghpull:`2408`: Fix Issue #2366
* :ghpull:`2405`: clarify TaskScheduler.hwm doc
* :ghpull:`2399`: IndentationError display
* :ghpull:`2400`: Add scroll_to_cell(cell_number) to the notebook
* :ghpull:`2401`: unmock read-the-docs modules
* :ghpull:`2311`: always perform requested trait assignments
* :ghpull:`2393`: New option `n` to limit history search hits
* :ghpull:`2386`: Adapt inline backend to changes in matplotlib
* :ghpull:`2392`: Remove suspicious double quote
* :ghpull:`2387`: Added -L library search path to cythonmagic cell magic
* :ghpull:`2370`: qtconsole: Create a prompt newline by inserting a new block (w/o formatting)
* :ghpull:`1715`: Fix for #1688, traceback-unicode issue
* :ghpull:`2378`: use Singleton.instance() for embed() instead of manual global
* :ghpull:`2373`: fix missing imports in core.interactiveshell
* :ghpull:`2368`: remove notification widget leftover
* :ghpull:`2327`: Parallel: Support get/set of nested objects in view (e.g. dv['a.b'])
* :ghpull:`2362`: Clean up ProgressBar class in example notebook
* :ghpull:`2346`: Extra xterm identification in set_term_title
* :ghpull:`2352`: Notebook: Store the username in a cookie whose name is unique.
* :ghpull:`2358`: add backport_pr to tools
* :ghpull:`2365`: fix names of notebooks for download/save
* :ghpull:`2364`: make clients use 'location' properly (fixes #2361)
* :ghpull:`2354`: Refactor notebook templates to use Jinja2
* :ghpull:`2339`: add bash completion example
* :ghpull:`2345`: Remove references to 'version' no longer in argparse. Github issue #2343.
* :ghpull:`2347`: adjust division error message checking to account for Python 3
* :ghpull:`2305`: RemoteError._render_traceback_ calls self.render_traceback
* :ghpull:`2338`: Normalize line endings for ipexec_validate, fix for #2315.
* :ghpull:`2192`: Introduce Notification Area
* :ghpull:`2329`: Better error messages for common magic commands.
* :ghpull:`2337`: ENH: added StackExchange-style MathJax filtering
* :ghpull:`2331`: update css for qtconsole in doc
* :ghpull:`2317`: adding cluster_id to parallel.Client.__init__
* :ghpull:`2130`: Add -l option to %R magic to allow passing in of local namespace
* :ghpull:`2196`: Fix for bad command line argument to latex
* :ghpull:`2300`: bug fix: was crashing when sqlite3 is not installed
* :ghpull:`2184`: Expose store_history to execute_request messages.
* :ghpull:`2308`: Add welcome_message option to enable_pylab
* :ghpull:`2302`: Fix variable expansion on 'self'
* :ghpull:`2299`: Remove code from prefilter that duplicates functionality in inputsplitter
* :ghpull:`2295`: allow pip install from github repository directly
* :ghpull:`2280`: fix SSH passwordless check for OpenSSH
* :ghpull:`2290`: nbmanager
* :ghpull:`2288`: s/assertEquals/assertEqual (again)
* :ghpull:`2287`: Removed outdated dev docs.
* :ghpull:`2218`: Use redirect for new notebooks
* :ghpull:`2277`: nb: up/down arrow keys move to begin/end of line at top/bottom of cell
* :ghpull:`2045`: Refactoring notebook managers and adding Azure backed storage.
* :ghpull:`2271`: use display instead of send_figure in inline backend hooks
* :ghpull:`2278`: allow disabling SQLite history
* :ghpull:`2225`: Add "--annotate" option to `%%cython` magic.
* :ghpull:`2246`: serialize individual args/kwargs rather than the containers
* :ghpull:`2274`: CLN: Use name to id mapping of notebooks instead of searching.
* :ghpull:`2270`: SSHLauncher tweaks
* :ghpull:`2269`: add missing location when disambiguating controller IP
* :ghpull:`2263`: Allow docs to build on https://readthedocs.io/
* :ghpull:`2256`: Adding data publication example notebook.
* :ghpull:`2255`: better flush iopub with AsyncResults
* :ghpull:`2261`: Fix: longest_substr([]) -> ''
* :ghpull:`2260`: fix mpr again
* :ghpull:`2242`: Document globbing in `%history -g <pattern>`.
* :ghpull:`2250`: fix html in notebook example
* :ghpull:`2245`: Fix regression in embed() from pull-request #2096.
* :ghpull:`2248`: track sha of master in test_pr messages
* :ghpull:`2238`: Fast tests
* :ghpull:`2211`: add data publication message
* :ghpull:`2236`: minor test_pr tweaks
* :ghpull:`2231`: Improve Image format validation and add html width,height
* :ghpull:`2232`: Reapply monkeypatch to inspect.findsource()
* :ghpull:`2235`: remove spurious print statement from setupbase.py
* :ghpull:`2222`: adjust how canning deals with import strings
* :ghpull:`2224`: fix css typo
* :ghpull:`2223`: Custom tracebacks
* :ghpull:`2214`: use KernelApp.exec_lines/files in IPEngineApp
* :ghpull:`2199`: Wrap JS published by %%javascript in try/catch
* :ghpull:`2212`: catch errors in markdown javascript
* :ghpull:`2190`: Update code mirror 2.22 to 2.32
* :ghpull:`2200`: documentation build broken in bb429da5b
* :ghpull:`2194`: clean nan/inf in json_clean
* :ghpull:`2198`: fix mpr for earlier git version
* :ghpull:`2175`: add FileFindHandler for Notebook static files
* :ghpull:`1990`: can func_defaults
* :ghpull:`2069`: start improving serialization in parallel code
* :ghpull:`2202`: Create a unique & temporary IPYTHONDIR for each testing group.
* :ghpull:`2204`: Work around lack of os.kill in win32.
* :ghpull:`2148`: win32 iptest: Use subprocess.Popen() instead of os.system().
* :ghpull:`2179`: Pylab switch
* :ghpull:`2124`: Add an API for registering magic aliases.
* :ghpull:`2169`: ipdb: pdef, pdoc, pinfo magics all broken
* :ghpull:`2174`: Ensure consistent indentation in `%magic`.
* :ghpull:`1930`: add size-limiting to the DictDB backend
* :ghpull:`2189`: Fix IPython.lib.latextools for Python 3
* :ghpull:`2186`: removed references to h5py dependence in octave magic documentation
* :ghpull:`2183`: Include the kernel object in the event object passed to kernel events
* :ghpull:`2185`: added test for %store, fixed storemagic
* :ghpull:`2138`: Use breqn.sty in dvipng backend if possible
* :ghpull:`2182`: handle undefined param in notebooklist
* :ghpull:`1831`: fix #1814 set __file__ when running .ipy files
* :ghpull:`2051`: Add a metadata attribute to messages
* :ghpull:`1471`: simplify IPython.parallel connections and enable Controller Resume
* :ghpull:`2181`: add %%javascript, %%svg, and %%latex display magics
* :ghpull:`2116`: different images in 00_notebook-tour
* :ghpull:`2092`: %prun: Restore `stats.stream` after running `print_stream`.
* :ghpull:`2159`: show message on notebook list if server is unreachable
* :ghpull:`2176`: fix git mpr
* :ghpull:`2152`: [qtconsole] Namespace not empty at startup
* :ghpull:`2177`: remove numpy install from travis/tox scripts
* :ghpull:`2090`: New keybinding for code cell execution + cell insertion
* :ghpull:`2160`: Updating the parallel options pricing example
* :ghpull:`2168`: expand line in cell magics
* :ghpull:`2170`: Fix tab completion with IPython.embed_kernel().
* :ghpull:`2096`: embed(): Default to the future compiler flags of the calling frame.
* :ghpull:`2163`: fix 'remote_profie_dir' typo in SSH launchers
* :ghpull:`2158`: [2to3 compat ] Tuple params in func defs
* :ghpull:`2089`: Fix unittest DeprecationWarnings
* :ghpull:`2142`: Refactor test_pr.py
* :ghpull:`2140`: 2to3: Apply `has_key` fixer.
* :ghpull:`2131`: Add option append (-a) to %save
* :ghpull:`2117`: use explicit url in notebook example
* :ghpull:`2133`: Tell git that ``*.py`` files contain Python code, for use in word-diffs.
* :ghpull:`2134`: Apply 2to3 `next` fix.
* :ghpull:`2126`: ipcluster broken with any batch launcher (PBS/LSF/SGE)
* :ghpull:`2104`: Windows make file for Sphinx documentation
* :ghpull:`2074`: Make BG color of inline plot configurable
* :ghpull:`2123`: BUG: Look up the `_repr_pretty_` method on the class within the MRO rath...
* :ghpull:`2100`: [in progress] python 2 and 3 compatibility without 2to3, second try
* :ghpull:`2128`: open notebook copy in different tabs
* :ghpull:`2073`: allows password and prefix for notebook
* :ghpull:`1993`: Print View
* :ghpull:`2086`: re-aliad %ed to %edit in qtconsole
* :ghpull:`2110`: Fixes and improvements to the input splitter
* :ghpull:`2101`: fix completer deletting newline
* :ghpull:`2102`: Fix logging on interactive shell.
* :ghpull:`2088`: Fix (some) Python 3.2 ResourceWarnings
* :ghpull:`2064`: conform to pep 3110
* :ghpull:`2076`: Skip notebook 'static' dir in test suite.
* :ghpull:`2063`: Remove umlauts so py3 installations on LANG=C systems succeed.
* :ghpull:`2068`: record sysinfo in sdist
* :ghpull:`2067`: update tools/release_windows.py
* :ghpull:`2065`: Fix parentheses typo
* :ghpull:`2062`: Remove duplicates and auto-generated files from repo.
* :ghpull:`2061`: use explicit tuple in exception
* :ghpull:`2060`: change minus to \- or \(hy in manpages

Issues (691):

* :ghissue:`3940`: Install process documentation overhaul 
* :ghissue:`3946`: The PDF option for `--post` should work with lowercase 
* :ghissue:`3957`: Notebook help page broken in Firefox
* :ghissue:`3894`: nbconvert test failure
* :ghissue:`3887`: 1.0.0a1 shows blank screen in both firefox and chrome (windows 7)
* :ghissue:`3703`: `nbconvert`: Output options -- names and documentation
* :ghissue:`3931`: Tab completion not working during debugging in the notebook
* :ghissue:`3936`: Ipcluster plugin is not working with Ipython 1.0dev
* :ghissue:`3941`: IPython Notebook kernel crash on Win7x64
* :ghissue:`3926`: Ending Notebook renaming dialog with return creates new-line
* :ghissue:`3932`: Incorrect empty docstring
* :ghissue:`3928`: Passing variables to script from the workspace
* :ghissue:`3774`: Notebooks with spaces in their names breaks nbconvert latex graphics
* :ghissue:`3916`: tornado needs its own check
* :ghissue:`3915`: Link to Parallel examples "found on GitHub" broken in docs
* :ghissue:`3895`: Keyboard shortcuts box in notebook doesn't fit the screen
* :ghissue:`3912`: IPython.utils fails automated test for RC1 1.0.0
* :ghissue:`3636`: Code cell missing highlight on load
* :ghissue:`3897`: under Windows, "ipython3 nbconvert "C:/blabla/first_try.ipynb" --to latex --post PDF" POST processing action fails because of a bad parameter
* :ghissue:`3900`: python3 install syntax errors (OS X 10.8.4)
* :ghissue:`3899`: nbconvert to latex fails on notebooks with spaces in file name
* :ghissue:`3881`: Temporary Working Directory Test Fails
* :ghissue:`2750`: A way to freeze code cells in the notebook
* :ghissue:`3893`: Resize Local Image Files in Notebook doesn't work
* :ghissue:`3823`: nbconvert on windows: tex and paths
* :ghissue:`3885`: under Windows, "ipython3 nbconvert "C:/blabla/first_try.ipynb" --to latex" write "\" instead of "/" to reference file path in the .tex file
* :ghissue:`3889`: test_qt fails due to assertion error 'qt4' != 'qt'
* :ghissue:`3890`: double post, disregard this issue
* :ghissue:`3689`: nbconvert, remaining tests
* :ghissue:`3874`: Up/Down keys don't work to "Search previous command history" (besides Ctrl-p/Ctrl-n)
* :ghissue:`3853`: CodeMirror locks up in the notebook
* :ghissue:`3862`: can only connect to an ipcluster started with v1.0.0-dev (master branch) using an older ipython (v0.13.2), but cannot connect using ipython (v1.0.0-dev)
* :ghissue:`3869`: custom css not working. 
* :ghissue:`2960`: Keyboard shortcuts
* :ghissue:`3795`: ipcontroller process goes to 100% CPU, ignores connection requests
* :ghissue:`3553`: Ipython and pylab crashes in windows and canopy
* :ghissue:`3837`: Cannot set custom mathjax url, crash notebook server.
* :ghissue:`3808`: "Naming" releases ?
* :ghissue:`2431`: TypeError: must be string without null bytes, not str
* :ghissue:`3856`: `?` at end of comment causes line to execute
* :ghissue:`3731`: nbconvert: add logging for the different steps of nbconvert
* :ghissue:`3835`: Markdown cells do not render correctly when mathjax is disabled
* :ghissue:`3843`: nbconvert to rst: leftover "In[ ]"
* :ghissue:`3799`: nbconvert: Ability to specify name of output file
* :ghissue:`3726`: Document when IPython.start_ipython() should be used versus IPython.embed()
* :ghissue:`3778`: Add no more readonly view in what's new
* :ghissue:`3754`: No Print View in Notebook in 1.0dev
* :ghissue:`3798`: IPython 0.12.1 Crashes on autocompleting sqlalchemy.func.row_number properties
* :ghissue:`3811`: Opening notebook directly from the command line with multi-directory support installed
* :ghissue:`3775`: Annoying behavior when clicking on cell after execution (Ctrl+Enter)
* :ghissue:`3809`: Possible to add some bpython features?
* :ghissue:`3810`: Printing the contents of an image file messes up shell text
* :ghissue:`3702`: `nbconvert`: Default help message should be that of --help
* :ghissue:`3735`: Nbconvert 1.0.0a1 does not take into account the pdf extensions in graphs
* :ghissue:`3719`: Bad strftime format, for windows, in nbconvert exporter 
* :ghissue:`3786`: Zmq errors appearing with `Ctrl-C` in console/qtconsole
* :ghissue:`3019`: disappearing scrollbar on tooltip in Chrome 24 on Ubuntu 12.04
* :ghissue:`3785`: ipdb completely broken in Qt console
* :ghissue:`3796`: Document the meaning of milestone/issues-tags for users.
* :ghissue:`3788`: Do not auto show tooltip if docstring empty.
* :ghissue:`1366`: [Web page] No link to front page from documentation
* :ghissue:`3739`: nbconvert (to slideshow) misses some of the math in markdown cells
* :ghissue:`3768`: increase and make timeout configurable in console completion.
* :ghissue:`3724`: ipcluster only running on one cpu
* :ghissue:`1592`: better message for unsupported nbformat
* :ghissue:`2049`: Can not stop "ipython kernel" on windows
* :ghissue:`3757`: Need direct entry point to given notebook 
* :ghissue:`3745`: ImportError: cannot import name check_linecache_ipython
* :ghissue:`3701`: `nbconvert`: Final output file should be in same directory as input file
* :ghissue:`3738`: history -o works but history with -n produces identical results
* :ghissue:`3740`: error when attempting to run 'make' in docs directory
* :ghissue:`3737`: ipython nbconvert crashes with ValueError: Invalid format string.
* :ghissue:`3730`: nbconvert: unhelpful error when pandoc isn't installed
* :ghissue:`3718`: markdown cell cursor misaligned in notebook
* :ghissue:`3710`: multiple input fields for %debug in the notebook after resetting the kernel
* :ghissue:`3713`: PyCharm has problems with IPython working inside PyPy created by virtualenv
* :ghissue:`3712`: Code completion: Complete on dictionary keys
* :ghissue:`3680`: --pylab and --matplotlib flag
* :ghissue:`3698`: nbconvert: Unicode error with minus sign
* :ghissue:`3693`: nbconvert does not process SVGs into PDFs
* :ghissue:`3688`: nbconvert, figures not extracting with Python 3.x
* :ghissue:`3542`: note new dependencies in docs / setup.py
* :ghissue:`2556`: [pagedown] do not target_blank anchor link
* :ghissue:`3684`: bad message when %pylab fails due import *other* than matplotlib
* :ghissue:`3682`: ipython notebook pylab inline  import_all=False 
* :ghissue:`3596`: MathjaxUtils race condition?
* :ghissue:`1540`: Comment/uncomment selection in notebook
* :ghissue:`2702`: frozen setup: permission denied for default ipython_dir
* :ghissue:`3672`: allow_none on Number-like traits.
* :ghissue:`2411`: add CONTRIBUTING.md
* :ghissue:`481`: IPython terminal issue with Qt4Agg on XP SP3
* :ghissue:`2664`: How to preserve user variables from import clashing?
* :ghissue:`3436`: enable_pylab(import_all=False) still imports np
* :ghissue:`2630`: lib.pylabtools.figsize : NameError when using Qt4Agg backend and %pylab magic. 
* :ghissue:`3154`: Notebook: no event triggered when a Cell is created
* :ghissue:`3579`: Nbconvert: SVG are not transformed to PDF anymore
* :ghissue:`3604`: MathJax rendering problem in `%%latex` cell
* :ghissue:`3668`: AttributeError: 'BlockingKernelClient' object has no attribute 'started_channels'
* :ghissue:`3245`: SyntaxError: encoding declaration in Unicode string
* :ghissue:`3639`: %pylab inline in IPYTHON notebook throws "RuntimeError: Cannot activate multiple GUI eventloops"
* :ghissue:`3663`: frontend deprecation warnings
* :ghissue:`3661`: run -m not behaving like python -m 
* :ghissue:`3597`: re-do PR #3531 - allow markdown in Header cell
* :ghissue:`3053`: Markdown in header cells is not rendered
* :ghissue:`3655`: IPython finding its way into pasted strings. 
* :ghissue:`3620`: uncaught errors in HTML output
* :ghissue:`3646`: get_dict() error
* :ghissue:`3004`: `%load_ext rmagic` fails when legacy ipy_user_conf.py is installed (in ipython 0.13.1 / OSX 10.8)
* :ghissue:`3638`: setp() issue in ipython notebook with figure references
* :ghissue:`3634`: nbconvert reveal to pdf conversion ignores styling, prints only a single page.
* :ghissue:`1307`: Remove pyreadline workarounds, we now require pyreadline >= 1.7.1
* :ghissue:`3316`: find_cmd test failure on Windows
* :ghissue:`3494`: input() in notebook doesn't work in Python 3
* :ghissue:`3427`: Deprecate `$` as mathjax delimiter
* :ghissue:`3625`: Pager does not open from button
* :ghissue:`3149`: Miscellaneous small nbconvert feedback
* :ghissue:`3617`: 256 color escapes support
* :ghissue:`3609`: %pylab inline blows up for single process ipython
* :ghissue:`2934`: Publish the Interactive MPI Demo Notebook
* :ghissue:`3614`: ansi escapes broken in master (ls --color)
* :ghissue:`3610`: If you don't have markdown, python setup.py install says no pygments
* :ghissue:`3547`: %run modules clobber each other
* :ghissue:`3602`: import_item fails when one tries to use DottedObjectName instead of a string
* :ghissue:`3563`: Duplicate tab completions in the notebook
* :ghissue:`3599`: Problems trying to run IPython on python3 without installing...
* :ghissue:`2937`: too long completion in notebook
* :ghissue:`3479`: Write empty name for the notebooks
* :ghissue:`3505`: nbconvert: Failure in specifying user filter
* :ghissue:`1537`: think a bit about namespaces
* :ghissue:`3124`: Long multiline strings in Notebook
* :ghissue:`3464`: run -d message unclear
* :ghissue:`2706`: IPython 0.13.1 ignoring $PYTHONSTARTUP
* :ghissue:`3587`: LaTeX escaping bug in nbconvert when exporting to HTML
* :ghissue:`3213`: Long running notebook died with a coredump
* :ghissue:`3580`: Running ipython with pypy on windows
* :ghissue:`3573`: custom.js not working
* :ghissue:`3544`: IPython.lib test failure on Windows
* :ghissue:`3352`: Install Sphinx extensions
* :ghissue:`2971`: [notebook]user needs to press ctrl-c twice to stop notebook server should be put into terminal window
* :ghissue:`2413`: ipython3 qtconsole fails to install: ipython 0.13 has no such extra feature 'qtconsole' 
* :ghissue:`2618`: documentation is incorrect for install process
* :ghissue:`2595`: mac 10.8 qtconsole export history
* :ghissue:`2586`: cannot store aliases
* :ghissue:`2714`: ipython qtconsole print unittest messages in console instead his own window. 
* :ghissue:`2669`: cython magic failing to work with openmp.
* :ghissue:`3256`: Vagrant pandas instance of IPython Notebook does not respect additional plotting arguments
* :ghissue:`3010`: cython magic fail if cache dir is deleted while in session
* :ghissue:`2044`: prune unused names from parallel.error
* :ghissue:`1145`: Online help utility broken in QtConsole
* :ghissue:`3439`: Markdown links no longer open in new window (with change from pagedown to marked)
* :ghissue:`3476`:  _margv  for macros seems to be missing
* :ghissue:`3499`: Add reveal.js library (version 2.4.0) inside IPython
* :ghissue:`2771`: Wiki Migration to GitHub
* :ghissue:`2887`: ipcontroller purging some engines during connect
* :ghissue:`626`: Enable Resuming Controller
* :ghissue:`2824`: Kernel restarting after message "Kernel XXXX failed to respond to heartbeat"
* :ghissue:`2823`: %%cython magic gives ImportError: dlopen(long_file_name.so, 2): image not found
* :ghissue:`2891`: In IPython for Python 3, system site-packages comes before user site-packages
* :ghissue:`2928`: Add magic "watch" function (example)
* :ghissue:`2931`: Problem rendering pandas dataframe in  Firefox for Windows
* :ghissue:`2939`: [notebook] Figure legend not shown in inline backend if ouside the box of the axes
* :ghissue:`2972`: [notebook] in Markdown mode, press Enter key at the end of <some http link>, the next line is indented unexpectly
* :ghissue:`3069`: Instructions for installing IPython notebook on Windows
* :ghissue:`3444`: Encoding problem: cannot use if user's name is not ascii?
* :ghissue:`3335`: Reenable bracket matching
* :ghissue:`3386`: Magic %paste not working in Python 3.3.2. TypeError: Type str doesn't support the buffer API
* :ghissue:`3543`: Exception shutting down kernel from notebook dashboard (0.13.1)
* :ghissue:`3549`: Codecell size changes with selection
* :ghissue:`3445`: Adding newlines in %%latex cell
* :ghissue:`3237`: [notebook] Can't close a notebook without errors
* :ghissue:`2916`: colon invokes auto(un)indent in markdown cells
* :ghissue:`2167`: Indent and dedent in htmlnotebook
* :ghissue:`3545`: Notebook save button icon not clear
* :ghissue:`3534`: nbconvert incompatible with Windows?
* :ghissue:`3489`: Update example notebook that raw_input is allowed
* :ghissue:`3396`: Notebook checkpoint time is displayed an hour out
* :ghissue:`3261`: Empty revert to checkpoint menu if no checkpoint...
* :ghissue:`2984`: "print" magic does not work in Python 3
* :ghissue:`3524`: Issues with pyzmq and ipython on EPD update
* :ghissue:`2434`: %store magic not auto-restoring
* :ghissue:`2720`: base_url and static path
* :ghissue:`2234`: Update various low resolution graphics for retina displays
* :ghissue:`2842`: Remember passwords for pw-protected notebooks
* :ghissue:`3244`: qtconsole: ValueError('close_fds is not supported on Windows platforms if you redirect stdin/stdout/stderr',)
* :ghissue:`2215`: AsyncResult.wait(0) can hang waiting for the client to get results?
* :ghissue:`2268`: provide mean to retrieve static data path
* :ghissue:`1905`: Expose UI for worksheets within each notebook
* :ghissue:`2380`: Qt inputhook prevents modal dialog boxes from displaying
* :ghissue:`3185`: prettify on double //
* :ghissue:`2821`: Test failure: IPython.parallel.tests.test_client.test_resubmit_header
* :ghissue:`2475`: [Notebook] Line is deindented when typing eg a colon in markdown mode
* :ghissue:`2470`: Do not destroy valid notebooks
* :ghissue:`860`: Allow the standalone export of a notebook to HTML
* :ghissue:`2652`: notebook with qt backend crashes at save image location popup
* :ghissue:`1587`: Improve kernel restarting in the notebook
* :ghissue:`2710`: Saving a plot in Mac OS X backend crashes IPython
* :ghissue:`2596`: notebook "Last saved:" is misleading on file opening.
* :ghissue:`2671`: TypeError :NoneType when executed "ipython qtconsole" in windows console
* :ghissue:`2703`: Notebook scrolling breaks after pager is shown
* :ghissue:`2803`: KernelManager and KernelClient should be two separate objects
* :ghissue:`2693`: TerminalIPythonApp configuration fails without ipython_config.py
* :ghissue:`2531`: IPython 0.13.1 python 2 32-bit installer includes 64-bit ipython*.exe launchers in the scripts folder
* :ghissue:`2520`: Control-C kills port forwarding
* :ghissue:`2279`: Setting `__file__` to None breaks Mayavi import
* :ghissue:`2161`: When logged into notebook, long titles are incorrectly positioned
* :ghissue:`1292`: Notebook, Print view should not be editable...
* :ghissue:`1731`: test parallel launchers
* :ghissue:`3227`: Improve documentation of ipcontroller and possible BUG
* :ghissue:`2896`: IPController very unstable
* :ghissue:`3517`: documentation build broken in head
* :ghissue:`3522`: UnicodeDecodeError: 'ascii' codec can't decode byte on Pycharm on Windows
* :ghissue:`3448`: Please include MathJax fonts with IPython Notebook
* :ghissue:`3519`: IPython Parallel map mysteriously turns pandas Series into numpy ndarray
* :ghissue:`3345`: IPython embedded shells ask if I want to exit, but I set confirm_exit = False
* :ghissue:`3509`: IPython won't close without asking "Are you sure?" in Firefox 
* :ghissue:`3471`: Notebook jinja2/markupsafe dependencies in manual
* :ghissue:`3502`: Notebook broken in master
* :ghissue:`3302`: autoreload does not work in ipython 0.13.x, python 3.3
* :ghissue:`3475`: no warning when leaving/closing notebook on master without saved changes
* :ghissue:`3490`: No obvious feedback when kernel crashes
* :ghissue:`1912`: Move all autoreload tests to their own group
* :ghissue:`2577`: sh.py and ipython for python 3.3
* :ghissue:`3467`: %magic doesn't work
* :ghissue:`3501`: Editing markdown cells that wrap has off-by-one errors in cursor positioning
* :ghissue:`3492`: IPython for Python3
* :ghissue:`3474`: unexpected keyword argument to remove_kernel
* :ghissue:`2283`: TypeError when using '?' after a string in a %logstart session
* :ghissue:`2787`: rmagic and pandas DataFrame
* :ghissue:`2605`: Ellipsis literal triggers AttributeError
* :ghissue:`1179`: Test unicode source in pinfo
* :ghissue:`2055`: drop Python 3.1 support
* :ghissue:`2293`: IPEP 2: Input transformations
* :ghissue:`2790`: %paste and %cpaste not removing "..." lines
* :ghissue:`3480`: Testing fails because iptest.py cannot be found
* :ghissue:`2580`: will not run within PIL build directory
* :ghissue:`2797`: RMagic, Dataframe Conversion Problem 
* :ghissue:`2838`: Empty lines disappear from triple-quoted literals.
* :ghissue:`3050`: Broken link on IPython.core.display page
* :ghissue:`3473`: Config not passed down to subcommands
* :ghissue:`3462`: Setting log_format in config file results in error (and no format changes)
* :ghissue:`3311`: Notebook (occasionally) not working on windows (Sophos AV)
* :ghissue:`3461`: Cursor positioning off by a character in auto-wrapped lines
* :ghissue:`3454`:  _repr_html_ error
* :ghissue:`3457`: Space in long Paragraph Markdown cell with Chinese or Japanese
* :ghissue:`3447`: Run Cell Does not Work
* :ghissue:`1373`: Last lines in long cells are hidden
* :ghissue:`1504`: Revisit serialization in IPython.parallel
* :ghissue:`1459`: Can't connect to 2 HTTPS notebook servers on the same host
* :ghissue:`678`: Input prompt stripping broken with multiline data structures
* :ghissue:`3001`: IPython.notebook.dirty flag is not set when a cell has unsaved changes
* :ghissue:`3077`: Multiprocessing semantics in parallel.view.map
* :ghissue:`3056`: links across notebooks
* :ghissue:`3120`: Tornado 3.0
* :ghissue:`3156`: update pretty to use Python 3 style for sets
* :ghissue:`3197`: Can't escape multiple dollar signs in a markdown cell
* :ghissue:`3309`: `Image()` signature/doc improvements
* :ghissue:`3415`: Bug in IPython/external/path/__init__.py 
* :ghissue:`3446`: Feature suggestion: Download matplotlib figure to client browser
* :ghissue:`3295`: autoexported notebooks: only export explicitly marked cells
* :ghissue:`3442`: Notebook: Summary table extracted from markdown headers
* :ghissue:`3438`: Zooming notebook in chrome is broken in master 
* :ghissue:`1378`: Implement autosave in notebook
* :ghissue:`3437`: Highlighting matching parentheses
* :ghissue:`3435`: module search segfault
* :ghissue:`3424`: ipcluster --version
* :ghissue:`3434`: 0.13.2 Ipython/genutils.py doesn't exist
* :ghissue:`3426`: Feature request: Save by cell and not by line #: IPython %save magic
* :ghissue:`3412`: Non Responsive Kernel: Running a Django development server from an IPython Notebook
* :ghissue:`3408`: Save cell toolbar and slide type metadata in notebooks
* :ghissue:`3246`: %paste regression with blank lines
* :ghissue:`3404`: Weird error with $variable and grep in command line magic (!command)
* :ghissue:`3405`: Key auto-completion in dictionaries?
* :ghissue:`3259`: Codemirror linenumber css broken
* :ghissue:`3397`: Vertical text misalignment in Markdown cells
* :ghissue:`3391`: Revert #3358 once fix integrated into CM
* :ghissue:`3360`: Error 500 while saving IPython notebook
* :ghissue:`3375`: Frequent Safari/Webkit crashes
* :ghissue:`3365`: zmq frontend
* :ghissue:`2654`: User_expression issues
* :ghissue:`3389`: Store history as plain text
* :ghissue:`3388`: Ipython parallel: open TCP connection created for each result returned from engine
* :ghissue:`3385`: setup.py failure on Python 3
* :ghissue:`3376`: Setting `__module__` to None breaks pretty printing
* :ghissue:`3374`: ipython qtconsole does not display the prompt on OSX
* :ghissue:`3380`: simple call to kernel
* :ghissue:`3379`: TaskRecord key 'started' not set
* :ghissue:`3241`: notebook connection time out
* :ghissue:`3334`: magic interpreter interprets non magic commands?
* :ghissue:`3326`: python3.3: Type error when launching SGE cluster in IPython notebook
* :ghissue:`3349`: pip3 doesn't run 2to3?
* :ghissue:`3347`: Longlist support in ipdb
* :ghissue:`3343`: Make pip install / easy_install faster
* :ghissue:`3337`: git submodules broke nightly PPA builds
* :ghissue:`3206`: Copy/Paste Regression in QtConsole
* :ghissue:`3329`: Buggy linewrap in Mac OSX Terminal (Mountain Lion)
* :ghissue:`3327`: Qt version check broken
* :ghissue:`3303`: parallel tasks never finish under heavy load
* :ghissue:`1381`: '\\' for equation continuations require an extra '\' in markdown cells
* :ghissue:`3314`: Error launching IPython
* :ghissue:`3306`: Test failure when running on a Vagrant VM
* :ghissue:`3280`: IPython.utils.process.getoutput returns stderr
* :ghissue:`3299`: variables named _ or __ exhibit incorrect behavior
* :ghissue:`3196`: add an "x" or similar to htmlnotebook pager
* :ghissue:`3293`: Several 404 errors for js files Firefox
* :ghissue:`3292`: syntax highlighting in chrome on OSX 10.8.3
* :ghissue:`3288`: Latest dev version hangs on page load
* :ghissue:`3283`: ipython dev retains directory information after directory change
* :ghissue:`3279`: custom.css is not overridden in the dev IPython (1.0)
* :ghissue:`2727`: %run -m doesn't support relative imports
* :ghissue:`3268`: GFM triple backquote and unknown language
* :ghissue:`3273`: Suppressing all plot related outputs
* :ghissue:`3272`: Backspace while completing load previous page
* :ghissue:`3260`: Js error in savewidget
* :ghissue:`3247`: scrollbar in notebook when not needed?
* :ghissue:`3243`: notebook: option to view json source from browser
* :ghissue:`3265`: 404 errors when running IPython 1.0dev 
* :ghissue:`3257`: setup.py not finding submodules
* :ghissue:`3253`: Incorrect Qt and PySide version comparison
* :ghissue:`3248`: Cell magics broken in Qt console
* :ghissue:`3012`: Problems with the less based style.min.css
* :ghissue:`2390`: Image width/height don't work in embedded images
* :ghissue:`3236`: cannot set TerminalIPythonApp.log_format
* :ghissue:`3214`: notebook kernel dies if started with invalid parameter
* :ghissue:`2980`: Remove HTMLCell ?
* :ghissue:`3128`: qtconsole hangs on importing pylab (using X forwarding)
* :ghissue:`3198`: Hitting recursive depth causing all notebook pages to hang
* :ghissue:`3218`: race conditions in profile directory creation
* :ghissue:`3177`: OverflowError execption in handlers.py
* :ghissue:`2563`: core.profiledir.check_startup_dir() doesn't work inside py2exe'd installation
* :ghissue:`3207`: [Feature] folders for ipython notebook dashboard
* :ghissue:`3178`: cell magics do not work with empty lines after #2447
* :ghissue:`3204`: Default plot() colors unsuitable for red-green colorblind users
* :ghissue:`1789`: ``:\n/*foo`` turns into ``:\n*(foo)`` in triple-quoted strings.
* :ghissue:`3202`: File cell magic fails with blank lines
* :ghissue:`3199`: %%cython -a stopped working?
* :ghissue:`2688`: obsolete imports in import autocompletion
* :ghissue:`3192`: Python2, Unhandled exception, __builtin__.True = False
* :ghissue:`3179`: script magic error message loop
* :ghissue:`3009`: use XDG_CACHE_HOME for cython objects
* :ghissue:`3059`: Bugs in 00_notebook_tour example.
* :ghissue:`3104`: Integrate a javascript file manager into the notebook front end
* :ghissue:`3176`: Particular equation not rendering  (notebook)
* :ghissue:`1133`: [notebook] readonly and upload files/UI
* :ghissue:`2975`: [notebook] python file and cell toolbar
* :ghissue:`3017`: SciPy.weave broken in IPython notebook/ qtconsole 
* :ghissue:`3161`: paste macro not reading spaces correctly
* :ghissue:`2835`: %paste not working on WinXpSP3/ipython-0.13.1.py2-win32-PROPER.exe/python27
* :ghissue:`2628`: Make transformers work for lines following decorators
* :ghissue:`2612`: Multiline String containing ":\n?foo\n" confuses interpreter to replace ?foo with get_ipython().magic(u'pinfo foo')
* :ghissue:`2539`: Request: Enable cell magics inside of .ipy scripts
* :ghissue:`2507`: Multiline string does not work (includes `...`) with doctest type input in IPython notebook
* :ghissue:`2164`: Request: Line breaks in line magic command
* :ghissue:`3106`: poor parallel performance with many jobs
* :ghissue:`2438`: print inside multiprocessing crashes Ipython kernel
* :ghissue:`3155`: Bad md5 hash for package 0.13.2
* :ghissue:`3045`: [Notebook] Ipython Kernel does not start if disconnected from internet(/network?)
* :ghissue:`3146`: Using celery in python 3.3
* :ghissue:`3145`: The notebook viewer is down
* :ghissue:`2385`: grep --color not working well with notebook
* :ghissue:`3131`: Quickly install from source in a clean virtualenv?
* :ghissue:`3139`: Rolling log for ipython
* :ghissue:`3127`: notebook with pylab=inline appears to call figure.draw twice
* :ghissue:`3129`: Walking up and down the call stack
* :ghissue:`3123`: Notebook crashed if unplugged ethernet cable
* :ghissue:`3121`: NB should use normalize.css? was #3049
* :ghissue:`3087`: Disable spellchecking in notebook
* :ghissue:`3084`: ipython pyqt 4.10 incompatibilty, QTextBlockUserData
* :ghissue:`3113`: Fails to install under Jython 2.7 beta
* :ghissue:`3110`: Render of h4 headers is not correct in notebook (error in renderedhtml.css)
* :ghissue:`3109`: BUG: read_csv: dtype={'id' : np.str}: Datatype not understood
* :ghissue:`3107`: Autocompletion of object attributes in arrays
* :ghissue:`3103`: Reset locale setting in qtconsole
* :ghissue:`3090`: python3.3 Entry Point not found
* :ghissue:`3081`: UnicodeDecodeError when using Image(data="some.jpeg")
* :ghissue:`2834`: url regexp only finds one link
* :ghissue:`3091`: qtconsole breaks doctest.testmod() in Python 3.3
* :ghissue:`3074`: SIGUSR1 not available on Windows
* :ghissue:`2996`: registration::purging stalled registration high occurrence in small clusters 
* :ghissue:`3065`: diff-ability of notebooks
* :ghissue:`3067`: Crash with pygit2
* :ghissue:`3061`: Bug handling Ellipsis
* :ghissue:`3049`: NB css inconsistent behavior between ff and webkit
* :ghissue:`3039`: unicode errors when opening a new notebook
* :ghissue:`3048`: Installning ipython qtConsole should be easyer att Windows
* :ghissue:`3042`: Profile creation fails on 0.13.2 branch
* :ghissue:`3035`: docstring typo/inconsistency: mention of an xml notebook format?
* :ghissue:`3031`: HDF5 library segfault (possibly due to mismatching headers?)
* :ghissue:`2991`: In notebook importing sympy closes ipython kernel
* :ghissue:`3027`: f.__globals__ causes an error in Python 3.3
* :ghissue:`3020`: Failing test test_interactiveshell.TestAstTransform on Windows
* :ghissue:`3023`: alt text for "click to expand output" has typo in alt text
* :ghissue:`2963`: %history to print all input history of a  previous session when line range is omitted
* :ghissue:`3018`: IPython installed within virtualenv. WARNING "Please install IPython inside the virtualtenv"
* :ghissue:`2484`: Completion in Emacs *Python* buffer causes prompt to be increased.
* :ghissue:`3014`: Ctrl-C finishes notebook immediately
* :ghissue:`3007`: cython_pyximport reload broken in python3
* :ghissue:`2955`: Incompatible Qt imports when running inprocess_qtconsole
* :ghissue:`3006`: [IPython 0.13.1] The check of PyQt version is wrong
* :ghissue:`3005`: Renaming a notebook to an existing notebook name overwrites the other file
* :ghissue:`2940`: Abort trap in IPython Notebook after installing matplotlib
* :ghissue:`3000`: issue #3000
* :ghissue:`2995`: ipython_directive.py fails on multiline when prompt number < 100
* :ghissue:`2993`: File magic (%%file) does not work with paths beginning with tilde (e.g., ~/anaconda/stuff.txt)
* :ghissue:`2992`: Cell-based input for console and qt frontends?
* :ghissue:`2425`: Liaise with Spyder devs to integrate newer IPython
* :ghissue:`2986`: requesting help in a loop can damage a notebook
* :ghissue:`2978`: v1.0-dev build errors on Arch with Python 3.
* :ghissue:`2557`: [refactor] Insert_cell_at_index()
* :ghissue:`2969`: ipython command does not work in terminal
* :ghissue:`2762`: OSX wxPython (osx_cocoa, 64bit) command "%gui wx" blocks the interpreter
* :ghissue:`2956`: Silent importing of submodules differs from standard Python3.2 interpreter's behavior
* :ghissue:`2943`: Up arrow key history search gets stuck in QTConsole
* :ghissue:`2953`: using 'nonlocal' declaration in global scope causes ipython3 crash
* :ghissue:`2952`: qtconsole ignores exec_lines
* :ghissue:`2949`: ipython crashes due to atexit()
* :ghissue:`2947`: From rmagic to  an R console
* :ghissue:`2938`: docstring pane not showing in notebook
* :ghissue:`2936`: Tornado assumes invalid signature for parse_qs on Python 3.1
* :ghissue:`2935`: unable to find python after easy_install / pip install
* :ghissue:`2920`: Add undo-cell deletion menu
* :ghissue:`2914`: BUG:saving a modified .py file after loading a module kills the kernel
* :ghissue:`2925`: BUG: kernel dies if user sets sys.stderr or sys.stdout to a file object
* :ghissue:`2909`: LaTeX sometimes fails to render in markdown cells with some curly bracket + underscore combinations
* :ghissue:`2898`: Skip ipc tests on Windows
* :ghissue:`2902`: ActiveState attempt to build ipython 0.12.1 for python 3.2.2 for Mac OS failed
* :ghissue:`2899`: Test failure in IPython.core.tests.test_magic.test_time
* :ghissue:`2890`: Test failure when fabric not installed
* :ghissue:`2892`: IPython tab completion bug for paths
* :ghissue:`1340`: Allow input cells to be collapsed
* :ghissue:`2881`: ? command in notebook does not show help in Safari
* :ghissue:`2751`: %%timeit should use minutes to format running time in long running cells
* :ghissue:`2879`: When importing a module with a wrong name, ipython crashes
* :ghissue:`2862`: %%timeit should warn of empty contents
* :ghissue:`2485`: History navigation breaks in qtconsole
* :ghissue:`2785`: gevent input hook
* :ghissue:`2843`: Sliently running code in clipboard (with paste, cpaste and variants)
* :ghissue:`2784`: %run -t -N<N> error
* :ghissue:`2732`: Test failure with FileLinks class on Windows
* :ghissue:`2860`: ipython help notebook -> KeyError: 'KernelManager'
* :ghissue:`2858`: Where is the installed `ipython` script?
* :ghissue:`2856`: Edit code entered from ipython in external editor
* :ghissue:`2722`: IPC transport option not taking effect ?
* :ghissue:`2473`: Better error messages in ipengine/ipcontroller
* :ghissue:`2836`: Cannot send builtin module definitions to IP engines
* :ghissue:`2833`: Any reason not to use super() ? 
* :ghissue:`2781`: Cannot interrupt infinite loops in the notebook
* :ghissue:`2150`: clippath_demo.py in matplotlib example does not work with inline backend
* :ghissue:`2634`: Numbered list in notebook markdown cell renders with Roman numerals instead of numbers
* :ghissue:`2230`: IPython crashing during startup with "AttributeError: 'NoneType' object has no attribute 'rstrip'"
* :ghissue:`2483`: nbviewer bug? with multi-file gists
* :ghissue:`2466`: mistyping `ed -p` breaks `ed -p`
* :ghissue:`2477`: Glob expansion tests fail on Windows
* :ghissue:`2622`: doc issue: notebooks that ship with Ipython .13 are written for python 2.x
* :ghissue:`2626`: Add "Cell -> Run All Keep Going" for notebooks
* :ghissue:`1223`: Show last modification date of each notebook
* :ghissue:`2621`: user request: put link to example notebooks in Dashboard
* :ghissue:`2564`: grid blanks plots in ipython pylab inline mode (interactive)
* :ghissue:`2532`: Django shell (IPython) gives NameError on dict comprehensions
* :ghissue:`2188`: ipython crashes on ctrl-c
* :ghissue:`2391`: Request: nbformat API to load/save without changing version
* :ghissue:`2355`: Restart kernel message even though kernel is perfectly alive
* :ghissue:`2306`: Garbled input text after reverse search on Mac OS X
* :ghissue:`2297`: ipdb with separate kernel/client pushing stdout to kernel process only
* :ghissue:`2180`: Have [kernel busy] overridden only by [kernel idle]
* :ghissue:`1188`: Pylab with OSX backend keyboard focus issue and hang
* :ghissue:`2107`: test_octavemagic.py[everything] fails 
* :ghissue:`1212`: Better understand/document browser compatibility
* :ghissue:`1585`: Refactor notebook templates to use Jinja2 and make each page a separate directory
* :ghissue:`1443`: xticks scaling factor partially obscured with qtconsole and inline plotting 
* :ghissue:`1209`: can't make %result work as in doc.
* :ghissue:`1200`: IPython 0.12 Windows install fails on Vista
* :ghissue:`1127`: Interactive test scripts for Qt/nb issues
* :ghissue:`959`: Matplotlib figures hide
* :ghissue:`2071`: win32 installer issue on Windows XP
* :ghissue:`2610`: ZMQInteractiveShell.colors being ignored
* :ghissue:`2505`: Markdown Cell incorrectly highlighting after "<"
* :ghissue:`165`: Installer fails to create Start Menu entries on Windows
* :ghissue:`2356`: failing traceback in terminal ipython for first exception
* :ghissue:`2145`: Have dashboad show when server disconect
* :ghissue:`2098`: Do not crash on kernel shutdow if json file is missing
* :ghissue:`2813`: Offline MathJax is broken on 0.14dev
* :ghissue:`2807`: Test failure: IPython.parallel.tests.test_client.TestClient.test_purge_everything
* :ghissue:`2486`: Readline's history search in ipython console does not clear properly after cancellation with Ctrl+C
* :ghissue:`2709`: Cython -la doesn't work
* :ghissue:`2767`: What is IPython.utils.upgradedir ?
* :ghissue:`2210`: Placing matplotlib legend outside axis bounds causes inline display to clip it
* :ghissue:`2553`: IPython Notebooks not robust against client failures
* :ghissue:`2536`: ImageDraw in Ipython notebook not drawing lines
* :ghissue:`2264`: Feature request: Versioning messaging protocol
* :ghissue:`2589`: Creation of ~300+ MPI-spawned engines causes instability in ipcluster
* :ghissue:`2672`: notebook: inline option without pylab
* :ghissue:`2673`: Indefinite Articles & Traitlets
* :ghissue:`2705`: Notebook crashes Safari with select and drag
* :ghissue:`2721`: dreload kills ipython when it hits zmq
* :ghissue:`2806`: ipython.parallel doesn't discover globals under Python 3.3
* :ghissue:`2794`: _exit_code behaves differently in terminal vs ZMQ frontends
* :ghissue:`2793`: IPython.parallel issue with pushing pandas TimeSeries
* :ghissue:`1085`: In process kernel for Qt frontend
* :ghissue:`2760`: IndexError: list index out of range with Python 3.2
* :ghissue:`2780`: Save and load notebooks from github
* :ghissue:`2772`: AttributeError: 'Client' object has no attribute 'kill'
* :ghissue:`2754`: Fail to send class definitions from interactive session to engines namespaces
* :ghissue:`2764`: TypeError while using 'cd'
* :ghissue:`2765`: name '__file__' is not defined
* :ghissue:`2540`: Wrap tooltip if line exceeds threshold?
* :ghissue:`2394`: Startup error on ipython qtconsole (version 0.13 and 0.14-dev
* :ghissue:`2440`: IPEP 4: Python 3 Compatibility
* :ghissue:`1814`: __file__ is not defined when file end with .ipy
* :ghissue:`2759`: R magic extension interferes with tab completion
* :ghissue:`2615`: Small change needed to rmagic extension.
* :ghissue:`2748`: collapse parts of a html notebook
* :ghissue:`1661`: %paste still bugs about IndentationError and says to use %paste
* :ghissue:`2742`: Octavemagic fails to deliver inline images in IPython (on Windows)
* :ghissue:`2739`: wiki.ipython.org contaminated with prescription drug spam
* :ghissue:`2588`: Link error while executing code from cython example notebook
* :ghissue:`2550`: Rpush magic doesn't find local variables and doesn't support comma separated lists of variables
* :ghissue:`2675`: Markdown/html blockquote need css.
* :ghissue:`2419`: TerminalInteractiveShell.__init__() ignores value of ipython_dir argument
* :ghissue:`1523`: Better LaTeX printing in the qtconsole with the sympy profile
* :ghissue:`2719`: ipython fails with `pkg_resources.DistributionNotFound: ipython==0.13`
* :ghissue:`2715`: url crashes nbviewer.ipython.org
* :ghissue:`2555`: "import" module completion on MacOSX
* :ghissue:`2707`: Problem installing the new version of IPython in Windows
* :ghissue:`2696`: SymPy magic bug in IPython Notebook
* :ghissue:`2684`: pretty print broken for types created with PyType_FromSpec
* :ghissue:`2533`: rmagic breaks on Windows
* :ghissue:`2661`: Qtconsole tooltip is too wide when the function has many arguments
* :ghissue:`2679`: ipython3 qtconsole via Homebrew on Mac OS X 10.8 - pyqt/pyside import error
* :ghissue:`2646`: pylab_not_importable
* :ghissue:`2587`: cython magic pops 2 CLI windows upon execution on Windows
* :ghissue:`2660`: Certain arguments (-h, --help, --version) never passed to scripts run with ipython
* :ghissue:`2665`: Missing docs for rmagic and some other extensions
* :ghissue:`2611`: Travis wants to drop 3.1 support
* :ghissue:`2658`: Incorrect parsing of raw multiline strings
* :ghissue:`2655`: Test fails if `from __future__ import print_function` in .pythonrc.py
* :ghissue:`2651`: nonlocal with no existing variable produces too many errors
* :ghissue:`2645`: python3 is a pain (minor unicode bug)
* :ghissue:`2637`: %paste in Python 3 on Mac doesn't work
* :ghissue:`2624`: Error on launching IPython on Win 7 and Python 2.7.3
* :ghissue:`2608`: disk IO activity on cursor press
* :ghissue:`1275`: Markdown parses LaTeX math symbols as its formatting syntax in notebook
* :ghissue:`2613`: display(Math(...)) doesn't render \tau correctly
* :ghissue:`925`: Tab-completion in Qt console needn't use pager
* :ghissue:`2607`: %load_ext sympy.interactive.ipythonprinting  dammaging output
* :ghissue:`2593`: Toolbar button to open qtconsole from notebook
* :ghissue:`2602`: IPython html documentation for downloading
* :ghissue:`2598`: ipython notebook --pylab=inline replaces built-in any()
* :ghissue:`2244`: small issue: wrong printout
* :ghissue:`2590`: add easier way to execute scripts in the current directory
* :ghissue:`2581`: %hist does not work when InteractiveShell.cache_size = 0
* :ghissue:`2584`: No file COPYING
* :ghissue:`2578`: AttributeError: 'module' object has no attribute 'TestCase'
* :ghissue:`2576`: One of my notebooks won't load any more -- is there a maximum notebook size?
* :ghissue:`2560`: Notebook output is invisible when printing strings with \r\r\n line endings
* :ghissue:`2566`: if pyside partially present ipython qtconsole fails to load even if pyqt4 present
* :ghissue:`1308`: ipython qtconsole  --ssh=server --existing ... hangs
* :ghissue:`1679`: List command doesn't work in ipdb debugger the first time
* :ghissue:`2545`: pypi win32 installer creates 64bit executibles
* :ghissue:`2080`: Event loop issues with IPython 0.12 and PyQt4 (``QDialog.exec_`` and more)
* :ghissue:`2541`: Allow `python -m IPython`
* :ghissue:`2508`: subplots_adjust() does not work correctly in ipython notebook
* :ghissue:`2289`: Incorrect mathjax rendering of certain arrays of equations
* :ghissue:`2487`: Selecting and indenting
* :ghissue:`2521`: more fine-grained 'run' controls, such as 'run from here' and 'run until here'
* :ghissue:`2535`: Funny bounding box when plot with text
* :ghissue:`2523`: History not working
* :ghissue:`2514`: Issue with zooming in qtconsole
* :ghissue:`2220`: No sys.stdout.encoding in kernel based IPython
* :ghissue:`2512`: ERROR: Internal Python error in the inspect module.
* :ghissue:`2496`: Function passwd does not work in QtConsole
* :ghissue:`1453`: make engines reconnect/die when controller was restarted
* :ghissue:`2481`: ipython notebook -- clicking in a code cell's output moves the screen to the top of the code cell
* :ghissue:`2488`: Undesired plot outputs in Notebook inline mode
* :ghissue:`2482`: ipython notebook -- download may not get the latest notebook
* :ghissue:`2471`: _subprocess module removed in Python 3.3
* :ghissue:`2374`: Issues with man pages
* :ghissue:`2316`: parallel.Client.__init__ should take cluster_id kwarg
* :ghissue:`2457`: Can a R library wrapper be created with Rmagic?
* :ghissue:`1575`: Fallback frontend for console when connecting pylab=inlnie -enabled kernel?
* :ghissue:`2097`: Do not crash if history db is corrupted
* :ghissue:`2435`: ipengines fail if clean_logs enabled
* :ghissue:`2429`: Using warnings.warn() results in TypeError
* :ghissue:`2422`: Multiprocessing in ipython notebook kernel crash
* :ghissue:`2426`: ipython crashes with the following message. I do not what went wrong. Can you help me identify the problem?
* :ghissue:`2423`: Docs typo?
* :ghissue:`2257`: pip install -e fails
* :ghissue:`2418`: rmagic can't run R's read.csv on data files with NA data
* :ghissue:`2417`: HTML notebook: Backspace sometimes deletes multiple characters
* :ghissue:`2275`: notebook: "Down_Arrow" on last line of cell should move to end of line
* :ghissue:`2414`: 0.13.1 does not work with current EPD 7.3-2
* :ghissue:`2409`: there is a redundant None
* :ghissue:`2410`: Use /usr/bin/python3 instead of /usr/bin/python
* :ghissue:`2366`: Notebook Dashboard --notebook-dir and fullpath
* :ghissue:`2406`: Inability to get docstring in debugger
* :ghissue:`2398`: Show line number for IndentationErrors
* :ghissue:`2314`: HTML lists seem to interfere with the QtConsole display
* :ghissue:`1688`: unicode exception when using %run with failing script
* :ghissue:`1884`: IPython.embed changes color on error
* :ghissue:`2381`: %time doesn't work for multiline statements
* :ghissue:`1435`: Add size keywords in Image class
* :ghissue:`2372`: interactiveshell.py misses urllib and io_open imports
* :ghissue:`2371`: IPython not working
* :ghissue:`2367`: Tab expansion moves to next cell in notebook
* :ghissue:`2359`: nbviever alters the order of print and display() output
* :ghissue:`2227`: print name for IPython Notebooks has become uninformative
* :ghissue:`2361`: client doesn't use connection file's 'location' in disambiguating 'interface'
* :ghissue:`2357`: failing traceback in terminal ipython for first exception
* :ghissue:`2343`: Installing in a python 3.3b2 or python 3.3rc1 virtual environment.
* :ghissue:`2315`: Failure in test: "Test we're not loading modules on startup that we shouldn't." 
* :ghissue:`2351`: Multiple Notebook Apps: cookies not port specific, clash with each other
* :ghissue:`2350`: running unittest from qtconsole prints output to terminal
* :ghissue:`2303`:  remote tracebacks broken since 952d0d6 (PR #2223)
* :ghissue:`2330`: qtconsole does not highlight tab-completion suggestion with custom stylesheet
* :ghissue:`2325`: Parsing Tex formula fails in Notebook
* :ghissue:`2324`: Parsing Tex formula fails
* :ghissue:`1474`: Add argument to `run -n` for custom namespace
* :ghissue:`2318`: C-m n/p don't work in Markdown cells in the notebook
* :ghissue:`2309`: time.time() in ipython notebook producing impossible results
* :ghissue:`2307`: schedule tasks on newly arrived engines
* :ghissue:`2313`: Allow Notebook HTML/JS to send messages to Python code
* :ghissue:`2304`: ipengine throws KeyError: url
* :ghissue:`1878`: shell access using ! will not fill class or function scope vars
* :ghissue:`2253`: %paste does not retrieve clipboard contents under screen/tmux on OS X
* :ghissue:`1510`: Add-on (or Monkey-patch) infrastructure for HTML notebook
* :ghissue:`2273`: triple quote and %s at beginning of line with %paste
* :ghissue:`2243`: Regression in .embed()
* :ghissue:`2266`: SSH passwordless check with OpenSSH checks for the wrong thing
* :ghissue:`2217`: Change NewNotebook handler to use 30x redirect
* :ghissue:`2276`: config option for disabling history store
* :ghissue:`2239`: can't use parallel.Reference in view.map
* :ghissue:`2272`: Sympy piecewise messed up rendering
* :ghissue:`2252`: %paste throws an exception with empty clipboard
* :ghissue:`2259`: git-mpr is currently broken
* :ghissue:`2247`: Variable expansion in shell commands should work in substrings
* :ghissue:`2026`: Run 'fast' tests only
* :ghissue:`2241`: read a list of notebooks on server and bring into browser only notebook
* :ghissue:`2237`: please put python and text editor in the web only ipython
* :ghissue:`2053`: Improvements to the IPython.display.Image object
* :ghissue:`1456`: ERROR: Internal Python error in the inspect module.
* :ghissue:`2221`: Avoid importing from IPython.parallel in core
* :ghissue:`2213`: Can't trigger startup code in Engines
* :ghissue:`1464`: Strange behavior for backspace with lines ending with more than 4 spaces in notebook 
* :ghissue:`2187`: NaN in object_info_reply JSON causes parse error
* :ghissue:`214`: system command requiring administrative privileges  
* :ghissue:`2195`: Unknown option `no-edit` in git-mpr
* :ghissue:`2201`: Add documentation build to tools/test_pr.py
* :ghissue:`2205`: Command-line option for default Notebook output collapsing behavior
* :ghissue:`1927`: toggle between inline and floating figures
* :ghissue:`2171`: Can't start StarCluster after upgrading to IPython 0.13
* :ghissue:`2173`: oct2py v >= 0.3.1 doesn't need h5py anymore
* :ghissue:`2099`: storemagic needs to use self.shell
* :ghissue:`2166`: DirectView map_sync() with Lambdas Using Generators
* :ghissue:`2091`: Unable to use print_stats after %prun -r in notebook
* :ghissue:`2132`: Add fail-over for pastebin
* :ghissue:`2156`: Make it possible to install ipython without nasty gui dependencies
* :ghissue:`2154`: Scrolled long output should be off in print view by default
* :ghissue:`2162`: Tab completion does not work with IPython.embed_kernel()
* :ghissue:`2157`: IPython 0.13 / github-master cannot create logfile from scratch
* :ghissue:`2151`: missing newline when a magic is called from the qtconsole menu
* :ghissue:`2139`: 00_notebook_tour Image example broken on master
* :ghissue:`2143`: Add a %%cython_annotate magic
* :ghissue:`2135`: Running IPython from terminal
* :ghissue:`2093`: Makefile for building Sphinx documentation on Windows 
* :ghissue:`2122`: Bug in pretty printing
* :ghissue:`2120`: Notebook "Make a Copy..." keeps opening duplicates in the same tab
* :ghissue:`1997`: password cannot be used with url prefix
* :ghissue:`2129`: help/doc displayed multiple times if requested in loop
* :ghissue:`2121`: ipdb does not support input history in qtconsole
* :ghissue:`2114`: %logstart doesn't log
* :ghissue:`2085`: %ed magic fails in qtconsole
* :ghissue:`2119`: IPython fails to run on MacOS Lion 
* :ghissue:`2052`: %pylab inline magic does not work on windows
* :ghissue:`2111`: Ipython won't start on W7
* :ghissue:`2112`: Strange internal traceback
* :ghissue:`2108`: Backslash (\) at the end of the line behavior different from default Python
* :ghissue:`1425`: Ampersands can't be typed sometimes in notebook cells
* :ghissue:`1513`: Add expand/collapse support for long output elements like stdout and tracebacks
* :ghissue:`2087`: error when starting ipython
* :ghissue:`2103`: Ability to run notebook file from commandline
* :ghissue:`2082`: Qt Console output spacing
* :ghissue:`2083`: Test failures with Python 3.2 and PYTHONWARNINGS="d"
* :ghissue:`2094`: about inline
* :ghissue:`2077`: Starting IPython3 on the terminal
* :ghissue:`1760`: easy_install ipython fails on py3.2-win32
* :ghissue:`2075`: Local Mathjax install causes iptest3 error under python3
* :ghissue:`2057`: setup fails for python3 with LANG=C
* :ghissue:`2070`: shebang on Windows
* :ghissue:`2054`: sys_info missing git hash in sdists
* :ghissue:`2059`: duplicate and modified files in documentation
* :ghissue:`2056`: except-shadows-builtin osm.py:687
* :ghissue:`2058`: hyphen-used-as-minus-sign in manpages
.. _issues_list_012:

Issues closed in the 0.12 development cycle
===========================================

Issues closed in 0.12.1
-----------------------

GitHub stats for bugfix release 0.12.1 (12/28/2011-04/16/2012), backporting
pull requests from 0.13.

We closed a total of 71 issues: 44 pull requests and 27 issues; this is the
full list (generated with the script `tools/github_stats.py`).

This list is automatically generated, and may be incomplete:

Pull Requests (44):

* :ghpull:`1175`: core.completer: Clean up excessive and unused code.
* :ghpull:`1187`: misc notebook: connection file cleanup, first heartbeat, startup flush
* :ghpull:`1190`: Fix link to Chris Fonnesbeck blog post about 0.11 highlights.
* :ghpull:`1196`: docs: looks like a file path might have been accidentally pasted in the middle of a word
* :ghpull:`1206`: don't preserve fixConsole output in json
* :ghpull:`1207`: fix loadpy duplicating newlines
* :ghpull:`1213`: BUG: Minor typo in history_console_widget.py
* :ghpull:`1218`: Added -q option to %prun for suppression of the output, along with editing the dochelp string.
* :ghpull:`1222`: allow Reference as callable in map/apply
* :ghpull:`1229`: Fix display of SyntaxError in Python 3
* :ghpull:`1246`: Skip tests that require X, when importing pylab results in RuntimeError.
* :ghpull:`1253`: set auto_create flag for notebook apps
* :ghpull:`1257`: use self.kernel_manager_class in qtconsoleapp
* :ghpull:`1262`: Heartbeat no longer shares the app's Context
* :ghpull:`1283`: HeartMonitor.period should be an Integer
* :ghpull:`1284`: a fix for GH 1269
* :ghpull:`1289`: Make autoreload extension work on Python 3.
* :ghpull:`1306`: Fix %prun input parsing for escaped characters (closes #1302)
* :ghpull:`1312`: minor heartbeat tweaks
* :ghpull:`1318`: make Ctrl-D in qtconsole act same as in terminal (ready to merge)
* :ghpull:`1341`: Don't attempt to tokenize binary files for tracebacks
* :ghpull:`1353`: Save notebook as script using unicode file handle.
* :ghpull:`1363`: Fix some minor color/style config issues in the qtconsole
* :ghpull:`1364`: avoid jsonlib returning Decimal
* :ghpull:`1369`: load header with engine id when engine dies in TaskScheduler
* :ghpull:`1370`: allow draft76 websockets (Safari)
* :ghpull:`1374`: remove calls to meaningless ZMQStream.on_err
* :ghpull:`1377`: Saving non-ascii history
* :ghpull:`1396`: Fix for %tb magic.
* :ghpull:`1402`: fix symlinked /home issue for FreeBSD
* :ghpull:`1413`: get_home_dir expands symlinks, adjust test accordingly
* :ghpull:`1414`: ignore errors in shell.var_expand
* :ghpull:`1430`: Fix for tornado check for tornado < 1.1.0
* :ghpull:`1445`: Don't build sphinx docs for sdists
* :ghpull:`1463`: Fix completion when importing modules in the cwd.
* :ghpull:`1477`: fix dangling `buffer` in IPython.parallel.util
* :ghpull:`1495`: BUG: Fix pretty-printing for overzealous objects
* :ghpull:`1496`: BUG: LBYL when clearing the output history on shutdown.
* :ghpull:`1514`: DOC: Fix references to IPython.lib.pretty instead of the old location
* :ghpull:`1517`: Fix indentation bug in IPython/lib/pretty.py
* :ghpull:`1538`: store git commit hash in utils._sysinfo instead of hidden data file
* :ghpull:`1599`: Fix for %run -d in Python 3
* :ghpull:`1602`: Fix %env for Python 3
* :ghpull:`1607`: cleanup sqlitedb temporary db file after tests

Issues (27):

* :ghissue:`676`: IPython.embed() from ipython crashes twice on exit
* :ghissue:`846`: Autoreload extension doesn't work with Python 3.2
* :ghissue:`1187`: misc notebook: connection file cleanup, first heartbeat, startup flush
* :ghissue:`1191`: profile/startup files not executed with "notebook"
* :ghissue:`1197`: Interactive shell trying to: from ... import history
* :ghissue:`1198`: Kernel Has Died error in Notebook
* :ghissue:`1201`: %env magic fails with Python 3.2
* :ghissue:`1204`: double newline from %loadpy in python notebook (at least on mac)
* :ghissue:`1208`: should dv.sync_import print failed imports ?
* :ghissue:`1225`: SyntaxError display broken in Python 3
* :ghissue:`1232`: Dead kernel loop
* :ghissue:`1241`: When our debugger class is used standalone `_oh` key errors are thrown
* :ghissue:`1254`: typo in notebooklist.js breaks links
* :ghissue:`1260`: heartbeat failure on long gil-holding operation
* :ghissue:`1268`: notebook %reset magic fails with StdinNotImplementedError
* :ghissue:`1269`: Another strange input handling error
* :ghissue:`1281`: in Hub: registration_timeout must be an integer, but heartmonitor.period is CFloat
* :ghissue:`1302`: Input parsing with %prun clobbers escapes
* :ghissue:`1304`: controller/server load can disrupt heartbeat
* :ghissue:`1317`: Very slow traceback construction from Cython extension
* :ghissue:`1345`: notebook can't save unicode as script
* :ghissue:`1375`: %history -g -f file encoding issue
* :ghissue:`1401`: numpy arrays cannot be used with View.apply() in Python 3
* :ghissue:`1408`: test_get_home_dir_3 failed on Mac OS X
* :ghissue:`1412`: Input parsing issue with %prun
* :ghissue:`1421`: ipython32 %run -d breaks with NameError name 'execfile' is not defined
* :ghissue:`1484`: unhide .git_commit_info.ini


Issues closed in 0.12
---------------------

In this cycle, from August 1 to December 28 2011, we closed a total of 515
issues, 257 pull requests and 258 regular issues; this is the full list
(generated with the script `tools/github_stats.py`).

Pull requests (257):

* `1174 <https://github.com/ipython/ipython/issues/1174>`_: Remove %install_default_config and %install_profiles
* `1178 <https://github.com/ipython/ipython/issues/1178>`_: Correct string type casting in pinfo.
* `1096 <https://github.com/ipython/ipython/issues/1096>`_: Show class init and call tooltips in notebook
* `1176 <https://github.com/ipython/ipython/issues/1176>`_: Modifications to profile list
* `1173 <https://github.com/ipython/ipython/issues/1173>`_: don't load gui/pylab in console frontend
* `1168 <https://github.com/ipython/ipython/issues/1168>`_: Add --script flag as shorthand for notebook save_script option.
* `1165 <https://github.com/ipython/ipython/issues/1165>`_: encode image_tag as utf8 in [x]html export
* `1161 <https://github.com/ipython/ipython/issues/1161>`_: Allow %loadpy to load remote URLs that don't end in .py
* `1158 <https://github.com/ipython/ipython/issues/1158>`_: Add coding header when notebook exported to .py file.
* `1160 <https://github.com/ipython/ipython/issues/1160>`_: don't ignore ctrl-C during `%gui qt`
* `1159 <https://github.com/ipython/ipython/issues/1159>`_: Add encoding header to Python files downloaded from notebooks.
* `1155 <https://github.com/ipython/ipython/issues/1155>`_: minor post-execute fixes (#1154)
* `1153 <https://github.com/ipython/ipython/issues/1153>`_: Pager tearing bug
* `1152 <https://github.com/ipython/ipython/issues/1152>`_: Add support for displaying maptlotlib axes directly.
* `1079 <https://github.com/ipython/ipython/issues/1079>`_: Login/out button cleanups
* `1151 <https://github.com/ipython/ipython/issues/1151>`_: allow access to user_ns in prompt_manager
* `1120 <https://github.com/ipython/ipython/issues/1120>`_: updated vim-ipython (pending)
* `1150 <https://github.com/ipython/ipython/issues/1150>`_: BUG: Scrolling pager in vsplit on Mac OSX tears.
* `1149 <https://github.com/ipython/ipython/issues/1149>`_: #1148 (win32 arg_split)
* `1147 <https://github.com/ipython/ipython/issues/1147>`_: Put qtconsole forground when launching
* `1146 <https://github.com/ipython/ipython/issues/1146>`_: allow saving notebook.py next to notebook.ipynb
* `1128 <https://github.com/ipython/ipython/issues/1128>`_: fix pylab StartMenu item
* `1140 <https://github.com/ipython/ipython/issues/1140>`_: Namespaces for embedding
* `1132 <https://github.com/ipython/ipython/issues/1132>`_: [notebook] read-only: disable name field
* `1125 <https://github.com/ipython/ipython/issues/1125>`_: notebook : update logo
* `1135 <https://github.com/ipython/ipython/issues/1135>`_: allow customized template and static file paths for the notebook web app
* `1122 <https://github.com/ipython/ipython/issues/1122>`_: BUG: Issue #755 qt IPythonWidget.execute_file fails if filename contains...
* `1137 <https://github.com/ipython/ipython/issues/1137>`_: rename MPIExecLaunchers to MPILaunchers
* `1130 <https://github.com/ipython/ipython/issues/1130>`_: optionally ignore  shlex's ValueError in arg_split
* `1116 <https://github.com/ipython/ipython/issues/1116>`_: Shlex unicode
* `1073 <https://github.com/ipython/ipython/issues/1073>`_: Storemagic plugin
* `1143 <https://github.com/ipython/ipython/issues/1143>`_: Add post_install script to create start menu entries in Python 3
* `1138 <https://github.com/ipython/ipython/issues/1138>`_: Fix tests to work when ~/.config/ipython contains a symlink.
* `1121 <https://github.com/ipython/ipython/issues/1121>`_: Don't transform function calls on IPyAutocall objects
* `1118 <https://github.com/ipython/ipython/issues/1118>`_: protect CRLF from carriage-return action
* `1105 <https://github.com/ipython/ipython/issues/1105>`_: Fix for prompts containing newlines.
* `1126 <https://github.com/ipython/ipython/issues/1126>`_: Totally remove pager when read only (notebook)
* `1091 <https://github.com/ipython/ipython/issues/1091>`_: qtconsole : allow copy with shortcut in pager
* `1114 <https://github.com/ipython/ipython/issues/1114>`_: fix magics history in two-process ipython console
* `1113 <https://github.com/ipython/ipython/issues/1113>`_: Fixing #1112 removing failing asserts for test_carriage_return and test_beep
* `1089 <https://github.com/ipython/ipython/issues/1089>`_: Support carriage return ('\r') and beep ('\b') characters in the qtconsole
* `1108 <https://github.com/ipython/ipython/issues/1108>`_: Completer usability 2 (rebased of  pr #1082)
* `864 <https://github.com/ipython/ipython/issues/864>`_: Two-process terminal frontend (ipython core branch)
* `1082 <https://github.com/ipython/ipython/issues/1082>`_: usability and cross browser compat for completer
* `1053 <https://github.com/ipython/ipython/issues/1053>`_: minor improvements to text placement in qtconsole
* `1106 <https://github.com/ipython/ipython/issues/1106>`_: Fix display of errors in compiled code on Python 3
* `1077 <https://github.com/ipython/ipython/issues/1077>`_: allow the notebook to run without MathJax
* `1072 <https://github.com/ipython/ipython/issues/1072>`_: If object has a getdoc() method, override its normal docstring.
* `1059 <https://github.com/ipython/ipython/issues/1059>`_: Switch to simple `__IPYTHON__` global
* `1070 <https://github.com/ipython/ipython/issues/1070>`_: Execution count after SyntaxError
* `1098 <https://github.com/ipython/ipython/issues/1098>`_: notebook: config section UI
* `1101 <https://github.com/ipython/ipython/issues/1101>`_: workaround spawnb missing from pexpect.__all__
* `1097 <https://github.com/ipython/ipython/issues/1097>`_: typo, should fix #1095
* `1099 <https://github.com/ipython/ipython/issues/1099>`_: qtconsole export xhtml/utf8
* `1083 <https://github.com/ipython/ipython/issues/1083>`_: Prompts
* `1081 <https://github.com/ipython/ipython/issues/1081>`_: Fix wildcard search for updated namespaces
* `1084 <https://github.com/ipython/ipython/issues/1084>`_: write busy in notebook window title...
* `1078 <https://github.com/ipython/ipython/issues/1078>`_: PromptManager fixes
* `1064 <https://github.com/ipython/ipython/issues/1064>`_: Win32 shlex
* `1069 <https://github.com/ipython/ipython/issues/1069>`_: As you type completer, fix on Firefox
* `1039 <https://github.com/ipython/ipython/issues/1039>`_: Base of an as you type completer.
* `1065 <https://github.com/ipython/ipython/issues/1065>`_: Qtconsole fix racecondition
* `507 <https://github.com/ipython/ipython/issues/507>`_: Prompt manager
* `1056 <https://github.com/ipython/ipython/issues/1056>`_: Warning in code. qtconsole ssh -X
* `1036 <https://github.com/ipython/ipython/issues/1036>`_: Clean up javascript based on js2-mode feedback.
* `1052 <https://github.com/ipython/ipython/issues/1052>`_: Pylab fix
* `648 <https://github.com/ipython/ipython/issues/648>`_: Usermod
* `969 <https://github.com/ipython/ipython/issues/969>`_: Pexpect-u
* `1007 <https://github.com/ipython/ipython/issues/1007>`_: Fix paste/cpaste bug and refactor/cleanup that code a lot.
* `506 <https://github.com/ipython/ipython/issues/506>`_: make ENTER on a previous input field replace current input buffer
* `1040 <https://github.com/ipython/ipython/issues/1040>`_: json/jsonapi cleanup
* `1042 <https://github.com/ipython/ipython/issues/1042>`_: fix firefox (windows) break line on empty prompt number
* `1015 <https://github.com/ipython/ipython/issues/1015>`_: emacs freezes when tab is hit in ipython with latest python-mode
* `1023 <https://github.com/ipython/ipython/issues/1023>`_: flush stdout/stderr at the end of kernel init
* `956 <https://github.com/ipython/ipython/issues/956>`_: Generate "All magics..." menu live
* `1038 <https://github.com/ipython/ipython/issues/1038>`_: Notebook: don't change cell when selecting code using shift+up/down.
* `987 <https://github.com/ipython/ipython/issues/987>`_: Add Tooltip to notebook.
* `1028 <https://github.com/ipython/ipython/issues/1028>`_: Cleaner minimum version comparison 
* `998 <https://github.com/ipython/ipython/issues/998>`_: defer to stdlib for path.get_home_dir()
* `1033 <https://github.com/ipython/ipython/issues/1033>`_: update copyright to 2011/20xx-2011
* `1032 <https://github.com/ipython/ipython/issues/1032>`_: Intercept <esc> avoid closing websocket on Firefox
* `1030 <https://github.com/ipython/ipython/issues/1030>`_: use pyzmq tools where appropriate
* `1029 <https://github.com/ipython/ipython/issues/1029>`_: Restore pspersistence, including %store magic, as an extension.
* `1025 <https://github.com/ipython/ipython/issues/1025>`_: Dollar escape
* `999 <https://github.com/ipython/ipython/issues/999>`_: Fix issue #880 - more useful message to user when %paste fails
* `938 <https://github.com/ipython/ipython/issues/938>`_: changes to get ipython.el to work with the latest python-mode.el
* `1012 <https://github.com/ipython/ipython/issues/1012>`_: Add logout button.
* `1020 <https://github.com/ipython/ipython/issues/1020>`_: Dollar formatter for ! shell calls
* `1019 <https://github.com/ipython/ipython/issues/1019>`_: Use repr() to make quoted strings
* `1008 <https://github.com/ipython/ipython/issues/1008>`_: don't use crash_handler by default
* `1003 <https://github.com/ipython/ipython/issues/1003>`_: Drop consecutive duplicates when refilling readline history
* `997 <https://github.com/ipython/ipython/issues/997>`_: don't unregister interrupted post-exec functions
* `996 <https://github.com/ipython/ipython/issues/996>`_: add Integer traitlet
* `1016 <https://github.com/ipython/ipython/issues/1016>`_: Fix password hashing for Python 3
* `1014 <https://github.com/ipython/ipython/issues/1014>`_: escape minus signs in manpages
* `1013 <https://github.com/ipython/ipython/issues/1013>`_: [NumPyExampleDocstring] link was pointing to raw file
* `1011 <https://github.com/ipython/ipython/issues/1011>`_: Add hashed password support.
* `1005 <https://github.com/ipython/ipython/issues/1005>`_: Quick fix for os.system requiring str parameter
* `994 <https://github.com/ipython/ipython/issues/994>`_: Allow latex formulas in HTML output
* `955 <https://github.com/ipython/ipython/issues/955>`_: Websocket Adjustments
* `979 <https://github.com/ipython/ipython/issues/979>`_: use system_raw in terminal, even on Windows
* `989 <https://github.com/ipython/ipython/issues/989>`_: fix arguments for commands in _process_posix
* `991 <https://github.com/ipython/ipython/issues/991>`_: Show traceback, continuing to start kernel if pylab init fails
* `981 <https://github.com/ipython/ipython/issues/981>`_: Split likely multiline text when writing JSON notebooks
* `957 <https://github.com/ipython/ipython/issues/957>`_: allow change of png DPI in inline backend
* `968 <https://github.com/ipython/ipython/issues/968>`_: add wantDirectory to ipdoctest, so that directories will be checked for e
* `984 <https://github.com/ipython/ipython/issues/984>`_: Do not expose variables defined at startup to %who etc.
* `985 <https://github.com/ipython/ipython/issues/985>`_: Fixes for parallel code on Python 3
* `963 <https://github.com/ipython/ipython/issues/963>`_: disable calltips in PySide < 1.0.7 to prevent segfault
* `976 <https://github.com/ipython/ipython/issues/976>`_: Getting started on what's new
* `929 <https://github.com/ipython/ipython/issues/929>`_: Multiline history
* `964 <https://github.com/ipython/ipython/issues/964>`_: Default profile
* `961 <https://github.com/ipython/ipython/issues/961>`_: Disable the pager for the test suite
* `953 <https://github.com/ipython/ipython/issues/953>`_: Physics extension
* `950 <https://github.com/ipython/ipython/issues/950>`_: Add directory for startup files
* `940 <https://github.com/ipython/ipython/issues/940>`_: allow setting HistoryManager.hist_file with config
* `948 <https://github.com/ipython/ipython/issues/948>`_: Monkeypatch Tornado 2.1.1 so it works with Google Chrome 16.
* `916 <https://github.com/ipython/ipython/issues/916>`_: Run p ( https://github.com/ipython/ipython/pull/901 )
* `923 <https://github.com/ipython/ipython/issues/923>`_: %config magic
* `920 <https://github.com/ipython/ipython/issues/920>`_: unordered iteration of AsyncMapResults (+ a couple fixes)
* `941 <https://github.com/ipython/ipython/issues/941>`_: Follow-up to 387dcd6a, `_rl.__doc__` is `None` with pyreadline
* `931 <https://github.com/ipython/ipython/issues/931>`_: read-only notebook mode
* `921 <https://github.com/ipython/ipython/issues/921>`_: Show invalid config message on TraitErrors during init
* `815 <https://github.com/ipython/ipython/issues/815>`_: Fix #481 using custom qt4 input hook
* `936 <https://github.com/ipython/ipython/issues/936>`_: Start webbrowser in a thread.  Prevents lockup with Chrome.
* `937 <https://github.com/ipython/ipython/issues/937>`_: add dirty trick for readline import on OSX
* `913 <https://github.com/ipython/ipython/issues/913>`_: Py3 tests2
* `933 <https://github.com/ipython/ipython/issues/933>`_: Cancel in qt console closeevent should trigger event.ignore()
* `930 <https://github.com/ipython/ipython/issues/930>`_: read-only notebook mode
* `910 <https://github.com/ipython/ipython/issues/910>`_: Make import checks more explicit in %whos
* `926 <https://github.com/ipython/ipython/issues/926>`_: reincarnate -V cmdline option
* `928 <https://github.com/ipython/ipython/issues/928>`_: BUG: Set context for font size change shortcuts in ConsoleWidget
* `901 <https://github.com/ipython/ipython/issues/901>`_:   - There is a bug when running the profiler in the magic command (prun) with python3
* `912 <https://github.com/ipython/ipython/issues/912>`_: Add magic for cls on windows. Fix for #181.
* `905 <https://github.com/ipython/ipython/issues/905>`_: enable %gui/%pylab magics in the Kernel
* `909 <https://github.com/ipython/ipython/issues/909>`_: Allow IPython to run without sqlite3
* `887 <https://github.com/ipython/ipython/issues/887>`_: Qtconsole menu
* `895 <https://github.com/ipython/ipython/issues/895>`_: notebook download implies save
* `896 <https://github.com/ipython/ipython/issues/896>`_: Execfile
* `899 <https://github.com/ipython/ipython/issues/899>`_: Brian's Notebook work
* `892 <https://github.com/ipython/ipython/issues/892>`_: don't close figures every cycle with inline matplotlib backend
* `893 <https://github.com/ipython/ipython/issues/893>`_: Adding clear_output to kernel and HTML notebook
* `789 <https://github.com/ipython/ipython/issues/789>`_: Adding clear_output to kernel and HTML notebook.
* `898 <https://github.com/ipython/ipython/issues/898>`_: Don't pass unicode sys.argv with %run or `ipython script.py`
* `897 <https://github.com/ipython/ipython/issues/897>`_: Add tooltips to the notebook via 'title' attr.
* `877 <https://github.com/ipython/ipython/issues/877>`_: partial fix for issue #678
* `838 <https://github.com/ipython/ipython/issues/838>`_: reenable multiline history for terminals
* `872 <https://github.com/ipython/ipython/issues/872>`_: The constructor of Client() checks for AssertionError in validate_url to open a file instead of connection to a URL if it fails.
* `884 <https://github.com/ipython/ipython/issues/884>`_: Notebook usability fixes
* `883 <https://github.com/ipython/ipython/issues/883>`_: User notification if notebook saving fails
* `889 <https://github.com/ipython/ipython/issues/889>`_: Add drop_by_id method to shell, to remove variables added by extensions.
* `891 <https://github.com/ipython/ipython/issues/891>`_: Ability to open the notebook in a browser when it starts
* `813 <https://github.com/ipython/ipython/issues/813>`_: Create menu bar for qtconsole
* `876 <https://github.com/ipython/ipython/issues/876>`_: protect IPython from bad custom exception handlers
* `856 <https://github.com/ipython/ipython/issues/856>`_: Backgroundjobs
* `868 <https://github.com/ipython/ipython/issues/868>`_: Warn user if MathJax can't be fetched from notebook closes #744
* `878 <https://github.com/ipython/ipython/issues/878>`_: store_history=False default for run_cell
* `824 <https://github.com/ipython/ipython/issues/824>`_: History access
* `850 <https://github.com/ipython/ipython/issues/850>`_: Update codemirror to 2.15 and make the code internally more version-agnostic
* `861 <https://github.com/ipython/ipython/issues/861>`_: Fix for issue #56
* `819 <https://github.com/ipython/ipython/issues/819>`_: Adding -m option to %run, similar to -m for python interpreter.
* `855 <https://github.com/ipython/ipython/issues/855>`_: promote aliases and flags, to ensure they have priority over config files
* `862 <https://github.com/ipython/ipython/issues/862>`_: BUG: Completion widget position and pager focus.
* `847 <https://github.com/ipython/ipython/issues/847>`_: Allow connection to kernels by files
* `708 <https://github.com/ipython/ipython/issues/708>`_: Two-process terminal frontend
* `857 <https://github.com/ipython/ipython/issues/857>`_: make sdist flags work again (e.g. --manifest-only)
* `835 <https://github.com/ipython/ipython/issues/835>`_: Add Tab key to list of keys that scroll down the paging widget.
* `859 <https://github.com/ipython/ipython/issues/859>`_: Fix for issue #800
* `848 <https://github.com/ipython/ipython/issues/848>`_: Python3 setup.py install failiure
* `845 <https://github.com/ipython/ipython/issues/845>`_: Tests on Python 3
* `802 <https://github.com/ipython/ipython/issues/802>`_: DOC: extensions: add documentation for the bundled extensions
* `830 <https://github.com/ipython/ipython/issues/830>`_: contiguous stdout/stderr in notebook
* `761 <https://github.com/ipython/ipython/issues/761>`_: Windows: test runner fails if repo path (e.g. home dir) contains spaces
* `801 <https://github.com/ipython/ipython/issues/801>`_: Py3 notebook
* `809 <https://github.com/ipython/ipython/issues/809>`_: use CFRunLoop directly in `ipython kernel --pylab osx`
* `841 <https://github.com/ipython/ipython/issues/841>`_: updated old scipy.org links, other minor doc fixes
* `837 <https://github.com/ipython/ipython/issues/837>`_: remove all trailling spaces
* `834 <https://github.com/ipython/ipython/issues/834>`_: Issue https://github.com/ipython/ipython/issues/832 resolution
* `746 <https://github.com/ipython/ipython/issues/746>`_: ENH: extensions: port autoreload to current API
* `828 <https://github.com/ipython/ipython/issues/828>`_: fixed permissions (sub-modules should not be executable) + added shebang  for run_ipy_in_profiler.py
* `798 <https://github.com/ipython/ipython/issues/798>`_: pexpect & Python 3
* `804 <https://github.com/ipython/ipython/issues/804>`_: Magic 'range' crash if greater than len(input_hist)
* `821 <https://github.com/ipython/ipython/issues/821>`_: update tornado dependency to 2.1
* `807 <https://github.com/ipython/ipython/issues/807>`_: Facilitate ssh tunnel sharing by announcing ports
* `795 <https://github.com/ipython/ipython/issues/795>`_: Add cluster-id for multiple cluster instances per profile
* `742 <https://github.com/ipython/ipython/issues/742>`_: Glut
* `668 <https://github.com/ipython/ipython/issues/668>`_: Greedy completer
* `776 <https://github.com/ipython/ipython/issues/776>`_: Reworking qtconsole shortcut, add fullscreen
* `790 <https://github.com/ipython/ipython/issues/790>`_: TST: add future unicode_literals test (#786)
* `775 <https://github.com/ipython/ipython/issues/775>`_: redirect_in/redirect_out should be constrained to windows only
* `793 <https://github.com/ipython/ipython/issues/793>`_: Don't use readline in the ZMQShell
* `743 <https://github.com/ipython/ipython/issues/743>`_: Pyglet
* `774 <https://github.com/ipython/ipython/issues/774>`_: basic/initial .mailmap for nice shortlog summaries
* `770 <https://github.com/ipython/ipython/issues/770>`_: #769 (reopened)
* `784 <https://github.com/ipython/ipython/issues/784>`_: Parse user code to AST using compiler flags.
* `783 <https://github.com/ipython/ipython/issues/783>`_: always use StringIO, never cStringIO
* `782 <https://github.com/ipython/ipython/issues/782>`_: flush stdout/stderr on displayhook call
* `622 <https://github.com/ipython/ipython/issues/622>`_: Make pylab import all configurable 
* `745 <https://github.com/ipython/ipython/issues/745>`_: Don't assume history requests succeed in qtconsole
* `725 <https://github.com/ipython/ipython/issues/725>`_: don't assume cursor.selectedText() is a string
* `778 <https://github.com/ipython/ipython/issues/778>`_: don't override execfile on Python 2
* `663 <https://github.com/ipython/ipython/issues/663>`_: Python 3 compatilibility work
* `762 <https://github.com/ipython/ipython/issues/762>`_: qtconsole ipython widget's execute_file fails if filename contains spaces or quotes
* `763 <https://github.com/ipython/ipython/issues/763>`_: Set context for shortcuts in ConsoleWidget
* `722 <https://github.com/ipython/ipython/issues/722>`_: PyPy compatibility
* `757 <https://github.com/ipython/ipython/issues/757>`_: ipython.el is broken in 0.11
* `764 <https://github.com/ipython/ipython/issues/764>`_: fix "--colors=<color>" option in py-python-command-args.
* `758 <https://github.com/ipython/ipython/issues/758>`_: use ROUTER/DEALER socket names instead of XREP/XREQ
* `736 <https://github.com/ipython/ipython/issues/736>`_: enh: added authentication ability for webapp
* `748 <https://github.com/ipython/ipython/issues/748>`_: Check for tornado before running frontend.html tests.
* `754 <https://github.com/ipython/ipython/issues/754>`_: restore msg_id/msg_type aliases in top level of msg dict
* `769 <https://github.com/ipython/ipython/issues/769>`_: Don't treat bytes objects as json-safe
* `753 <https://github.com/ipython/ipython/issues/753>`_: DOC: msg['msg_type'] removed
* `766 <https://github.com/ipython/ipython/issues/766>`_: fix "--colors=<color>" option in py-python-command-args.
* `765 <https://github.com/ipython/ipython/issues/765>`_: fix "--colors=<color>" option in py-python-command-args.
* `741 <https://github.com/ipython/ipython/issues/741>`_: Run PyOs_InputHook in pager to keep plot windows interactive.
* `664 <https://github.com/ipython/ipython/issues/664>`_: Remove ipythonrc references from documentation
* `750 <https://github.com/ipython/ipython/issues/750>`_: Tiny doc fixes
* `433 <https://github.com/ipython/ipython/issues/433>`_: ZMQ terminal frontend
* `734 <https://github.com/ipython/ipython/issues/734>`_: Allow %magic argument filenames with spaces to be specified with quotes under win32
* `731 <https://github.com/ipython/ipython/issues/731>`_: respect encoding of display data from urls
* `730 <https://github.com/ipython/ipython/issues/730>`_: doc improvements for running notebook via secure protocol
* `729 <https://github.com/ipython/ipython/issues/729>`_: use null char to start markdown cell placeholder
* `727 <https://github.com/ipython/ipython/issues/727>`_: Minor fixes to the htmlnotebook
* `726 <https://github.com/ipython/ipython/issues/726>`_: use bundled argparse if system argparse is < 1.1
* `705 <https://github.com/ipython/ipython/issues/705>`_: Htmlnotebook
* `723 <https://github.com/ipython/ipython/issues/723>`_: Add 'import time' to IPython/parallel/apps/launcher.py as time.sleep is called without time being imported
* `714 <https://github.com/ipython/ipython/issues/714>`_: Install mathjax for offline use
* `718 <https://github.com/ipython/ipython/issues/718>`_: Underline keyboard shortcut characters on appropriate buttons
* `717 <https://github.com/ipython/ipython/issues/717>`_: Add source highlighting to markdown snippets
* `716 <https://github.com/ipython/ipython/issues/716>`_: update EvalFormatter to allow arbitrary expressions
* `712 <https://github.com/ipython/ipython/issues/712>`_: Reset execution counter after cache is cleared
* `713 <https://github.com/ipython/ipython/issues/713>`_: Align colons in html notebook help dialog
* `709 <https://github.com/ipython/ipython/issues/709>`_: Allow usage of '.' in notebook names
* `706 <https://github.com/ipython/ipython/issues/706>`_: Implement static publishing of HTML notebook
* `674 <https://github.com/ipython/ipython/issues/674>`_: use argparse to parse aliases & flags
* `679 <https://github.com/ipython/ipython/issues/679>`_: HistoryManager.get_session_info()
* `696 <https://github.com/ipython/ipython/issues/696>`_: Fix columnize bug, where tab completion with very long filenames would crash Qt console
* `686 <https://github.com/ipython/ipython/issues/686>`_: add ssh tunnel support to qtconsole
* `685 <https://github.com/ipython/ipython/issues/685>`_: Add SSH tunneling to engines
* `384 <https://github.com/ipython/ipython/issues/384>`_: Allow pickling objects defined interactively.
* `647 <https://github.com/ipython/ipython/issues/647>`_: My fix rpmlint
* `587 <https://github.com/ipython/ipython/issues/587>`_: don't special case for py3k+numpy
* `703 <https://github.com/ipython/ipython/issues/703>`_: make config-loading debug messages more explicit
* `699 <https://github.com/ipython/ipython/issues/699>`_: make calltips configurable in qtconsole
* `666 <https://github.com/ipython/ipython/issues/666>`_: parallel tests & extra readline escapes
* `683 <https://github.com/ipython/ipython/issues/683>`_: BF - allow nose with-doctest setting in environment
* `689 <https://github.com/ipython/ipython/issues/689>`_: Protect ipkernel from bad messages
* `702 <https://github.com/ipython/ipython/issues/702>`_: Prevent ipython.py launcher from being imported.
* `701 <https://github.com/ipython/ipython/issues/701>`_: Prevent ipython.py from being imported by accident
* `670 <https://github.com/ipython/ipython/issues/670>`_: check for writable dirs, not just existence, in utils.path
* `579 <https://github.com/ipython/ipython/issues/579>`_: Sessionwork
* `687 <https://github.com/ipython/ipython/issues/687>`_: add `ipython kernel` for starting just a kernel
* `627 <https://github.com/ipython/ipython/issues/627>`_: Qt Console history search
* `646 <https://github.com/ipython/ipython/issues/646>`_: Generate package list automatically in find_packages
* `660 <https://github.com/ipython/ipython/issues/660>`_: i658
* `659 <https://github.com/ipython/ipython/issues/659>`_: don't crash on bad config files

Regular issues (258):

* `1177 <https://github.com/ipython/ipython/issues/1177>`_: UnicodeDecodeError in py3compat from "xlrd??"
* `1094 <https://github.com/ipython/ipython/issues/1094>`_: Tooltip doesn't show constructor docstrings
* `1170 <https://github.com/ipython/ipython/issues/1170>`_: double pylab greeting with c.InteractiveShellApp.pylab = "tk" in zmqconsole
* `1166 <https://github.com/ipython/ipython/issues/1166>`_: E-mail cpaste broken
* `1164 <https://github.com/ipython/ipython/issues/1164>`_: IPython qtconsole (0.12) can't export to html with external png
* `1103 <https://github.com/ipython/ipython/issues/1103>`_: %loadpy should cut out encoding declaration
* `1156 <https://github.com/ipython/ipython/issues/1156>`_: Notebooks downloaded as Python files require a header stating the encoding
* `1157 <https://github.com/ipython/ipython/issues/1157>`_: Ctrl-C not working when GUI/pylab integration is active
* `1154 <https://github.com/ipython/ipython/issues/1154>`_: We should be less aggressive in de-registering post-execution functions
* `1134 <https://github.com/ipython/ipython/issues/1134>`_: "select-all, kill" leaves qtconsole in unusable state
* `1148 <https://github.com/ipython/ipython/issues/1148>`_: A lot of testerrors
* `803 <https://github.com/ipython/ipython/issues/803>`_: Make doctests work with Python 3
* `1119 <https://github.com/ipython/ipython/issues/1119>`_: Start menu shortcuts not created in Python 3
* `1136 <https://github.com/ipython/ipython/issues/1136>`_: The embedding machinery ignores user_ns
* `607 <https://github.com/ipython/ipython/issues/607>`_: Use the new IPython logo/font in the notebook header
* `755 <https://github.com/ipython/ipython/issues/755>`_: qtconsole ipython widget's execute_file fails if filename contains spaces or quotes
* `1115 <https://github.com/ipython/ipython/issues/1115>`_: shlex_split should return unicode
* `1109 <https://github.com/ipython/ipython/issues/1109>`_: timeit with string ending in space gives "ValueError: No closing quotation"
* `1142 <https://github.com/ipython/ipython/issues/1142>`_: Install problems
* `700 <https://github.com/ipython/ipython/issues/700>`_: Some SVG images render incorrectly in htmlnotebook
* `1117 <https://github.com/ipython/ipython/issues/1117>`_: quit() doesn't work in terminal
* `1111 <https://github.com/ipython/ipython/issues/1111>`_: ls broken after merge of #1089
* `1104 <https://github.com/ipython/ipython/issues/1104>`_: Prompt spacing weird
* `1124 <https://github.com/ipython/ipython/issues/1124>`_: Seg Fault 11 when calling PySide using "run" command
* `1088 <https://github.com/ipython/ipython/issues/1088>`_: QtConsole : can't copy from pager
* `568 <https://github.com/ipython/ipython/issues/568>`_: Test error and failure in IPython.core on windows
* `1112 <https://github.com/ipython/ipython/issues/1112>`_: testfailure in IPython.frontend on windows
* `1102 <https://github.com/ipython/ipython/issues/1102>`_: magic in IPythonDemo fails when not located at top of demo file
* `629 <https://github.com/ipython/ipython/issues/629>`_: \r and \b in qtconsole don't behave as expected
* `1080 <https://github.com/ipython/ipython/issues/1080>`_: Notebook: tab completion should close on "("
* `973 <https://github.com/ipython/ipython/issues/973>`_: Qt Console close dialog and on-top Qt Console
* `1087 <https://github.com/ipython/ipython/issues/1087>`_: QtConsole xhtml/Svg export broken ?
* `1067 <https://github.com/ipython/ipython/issues/1067>`_: Parallel test suite hangs on Python 3
* `1018 <https://github.com/ipython/ipython/issues/1018>`_: Local mathjax breaks install
* `993 <https://github.com/ipython/ipython/issues/993>`_: `raw_input` redirection to foreign kernels is extremely brittle
* `1100 <https://github.com/ipython/ipython/issues/1100>`_: ipython3 traceback unicode issue from extensions
* `1071 <https://github.com/ipython/ipython/issues/1071>`_: Large html-notebooks hang on load on a slow machine
* `89 <https://github.com/ipython/ipython/issues/89>`_: %pdoc np.ma.compress shows docstring twice
* `22 <https://github.com/ipython/ipython/issues/22>`_: Include improvements from anythingipython.el
* `633 <https://github.com/ipython/ipython/issues/633>`_: Execution count & SyntaxError
* `1095 <https://github.com/ipython/ipython/issues/1095>`_: Uncaught TypeError: Object has no method 'remove_and_cancell_tooltip'
* `1075 <https://github.com/ipython/ipython/issues/1075>`_: We're ignoring prompt customizations
* `1086 <https://github.com/ipython/ipython/issues/1086>`_: Can't open qtconsole from outside source tree
* `1076 <https://github.com/ipython/ipython/issues/1076>`_: namespace changes broke `foo.*bar*?` syntax
* `1074 <https://github.com/ipython/ipython/issues/1074>`_: pprinting old-style class objects fails (TypeError: 'tuple' object is not callable)
* `1063 <https://github.com/ipython/ipython/issues/1063>`_: IPython.utils test error due to missing unicodedata module
* `592 <https://github.com/ipython/ipython/issues/592>`_: Bug in argument parsing for %run
* `378 <https://github.com/ipython/ipython/issues/378>`_: Windows path escape issues
* `1068 <https://github.com/ipython/ipython/issues/1068>`_: Notebook tab completion broken in Firefox
* `75 <https://github.com/ipython/ipython/issues/75>`_: No tab completion after "/
* `103 <https://github.com/ipython/ipython/issues/103>`_: customizable cpaste
* `324 <https://github.com/ipython/ipython/issues/324>`_: Remove code in IPython.testing that is not being used
* `131 <https://github.com/ipython/ipython/issues/131>`_: Global variables not seen by cprofile.run()
* `851 <https://github.com/ipython/ipython/issues/851>`_: IPython shell swallows exceptions in certain circumstances
* `882 <https://github.com/ipython/ipython/issues/882>`_: ipython freezes at start if IPYTHONDIR is on an NFS mount
* `1057 <https://github.com/ipython/ipython/issues/1057>`_: Blocker: Qt console broken after "all magics" menu became dynamic
* `1027 <https://github.com/ipython/ipython/issues/1027>`_: ipython does not like white space at end of file
* `1058 <https://github.com/ipython/ipython/issues/1058>`_: New bug: Notebook asks for confirmation to leave even saved pages.
* `1061 <https://github.com/ipython/ipython/issues/1061>`_: rep (magic recall) under pypy
* `1047 <https://github.com/ipython/ipython/issues/1047>`_: Document the notebook format
* `102 <https://github.com/ipython/ipython/issues/102>`_: Properties accessed twice for classes defined interactively
* `16 <https://github.com/ipython/ipython/issues/16>`_: %store raises exception when storing compiled regex
* `67 <https://github.com/ipython/ipython/issues/67>`_: tab expansion should only take one directory level at the time
* `62 <https://github.com/ipython/ipython/issues/62>`_: Global variables undefined in interactive use of embedded ipython shell
* `57 <https://github.com/ipython/ipython/issues/57>`_: debugging with ipython does not work well outside ipython
* `38 <https://github.com/ipython/ipython/issues/38>`_: Line entry edge case error
* `980 <https://github.com/ipython/ipython/issues/980>`_: Update parallel docs for new parallel architecture
* `1017 <https://github.com/ipython/ipython/issues/1017>`_: Add small example about ipcluster/ssh startup
* `1041 <https://github.com/ipython/ipython/issues/1041>`_: Proxy Issues
* `967 <https://github.com/ipython/ipython/issues/967>`_: KernelManagers don't use zmq eventloop properly
* `1055 <https://github.com/ipython/ipython/issues/1055>`_: "All Magics" display on Ubuntu 
* `1054 <https://github.com/ipython/ipython/issues/1054>`_: ipython explodes on syntax error
* `1051 <https://github.com/ipython/ipython/issues/1051>`_: ipython3 set_next_input() failure
* `693 <https://github.com/ipython/ipython/issues/693>`_: "run -i" no longer works after %reset in terminal
* `29 <https://github.com/ipython/ipython/issues/29>`_: cPickle works in standard interpreter, but not in IPython
* `1050 <https://github.com/ipython/ipython/issues/1050>`_: ipython3 broken by commit 8bb887c8c2c447bf7
* `1048 <https://github.com/ipython/ipython/issues/1048>`_: Update docs on notebook password
* `1046 <https://github.com/ipython/ipython/issues/1046>`_: Searies of questions/issues?
* `1045 <https://github.com/ipython/ipython/issues/1045>`_: crash when exiting - previously launched embedded sub-shell
* `1043 <https://github.com/ipython/ipython/issues/1043>`_: pylab doesn't work in qtconsole
* `1044 <https://github.com/ipython/ipython/issues/1044>`_: run -p doesn't work in python 3
* `1010 <https://github.com/ipython/ipython/issues/1010>`_: emacs freezes when ipython-complete is called
* `82 <https://github.com/ipython/ipython/issues/82>`_: Update devel docs with discussion about good changelogs
* `116 <https://github.com/ipython/ipython/issues/116>`_: Update release management scipts and release.revision for git
* `1022 <https://github.com/ipython/ipython/issues/1022>`_: Pylab banner shows up with first cell to execute
* `787 <https://github.com/ipython/ipython/issues/787>`_: Keyboard selection of multiple lines in the notebook behaves inconsistently
* `1037 <https://github.com/ipython/ipython/issues/1037>`_: notepad + jsonlib: TypeError: Only whitespace may be used for indentation.
* `970 <https://github.com/ipython/ipython/issues/970>`_: Default home not writable, %HOME% does not help (windows)
* `747 <https://github.com/ipython/ipython/issues/747>`_: HOMESHARE not a good choice for "writable homedir" on Windows
* `810 <https://github.com/ipython/ipython/issues/810>`_: cleanup utils.path.get_home_dir
* `2 <https://github.com/ipython/ipython/issues/2>`_: Fix the copyright statement in source code files to be accurate
* `1031 <https://github.com/ipython/ipython/issues/1031>`_: <esc> on Firefox crash websocket
* `684 <https://github.com/ipython/ipython/issues/684>`_: %Store eliminated in configuration and magic commands in 0.11
* `1026 <https://github.com/ipython/ipython/issues/1026>`_: BUG: wrong default parameter in ask_yes_no
* `880 <https://github.com/ipython/ipython/issues/880>`_: Better error message if %paste fails
* `1024 <https://github.com/ipython/ipython/issues/1024>`_: autopx magic broken 
* `822 <https://github.com/ipython/ipython/issues/822>`_: Unicode bug in Itpl when expanding shell variables in syscalls with !
* `1009 <https://github.com/ipython/ipython/issues/1009>`_: Windows: regression in cd magic handling of paths
* `833 <https://github.com/ipython/ipython/issues/833>`_: Crash python with matplotlib and unequal length arrays
* `695 <https://github.com/ipython/ipython/issues/695>`_: Crash handler initialization is too aggressive
* `1000 <https://github.com/ipython/ipython/issues/1000>`_: Remove duplicates when refilling readline history
* `992 <https://github.com/ipython/ipython/issues/992>`_: Interrupting certain matplotlib operations leaves the inline backend 'wedged'
* `942 <https://github.com/ipython/ipython/issues/942>`_: number traits should cast if value doesn't change
* `1006 <https://github.com/ipython/ipython/issues/1006>`_: ls crashes when run on a UNC path or with non-ascii args
* `944 <https://github.com/ipython/ipython/issues/944>`_: Decide the default image format for inline figures: SVG or PNG?
* `842 <https://github.com/ipython/ipython/issues/842>`_: Python 3 on Windows (pyreadline) - expected an object with the buffer interface
* `1002 <https://github.com/ipython/ipython/issues/1002>`_: ImportError due to incorrect version checking
* `1001 <https://github.com/ipython/ipython/issues/1001>`_: Ipython "source" command?
* `954 <https://github.com/ipython/ipython/issues/954>`_: IPython embed doesn't respect namespaces
* `681 <https://github.com/ipython/ipython/issues/681>`_: pdb freezes inside qtconsole
* `698 <https://github.com/ipython/ipython/issues/698>`_: crash report "TypeError: can only concatenate list (not "unicode") to list"
* `978 <https://github.com/ipython/ipython/issues/978>`_: ipython 0.11 buffers external command output till the cmd is done
* `952 <https://github.com/ipython/ipython/issues/952>`_: Need user-facing warning in the browser if websocket connection fails
* `988 <https://github.com/ipython/ipython/issues/988>`_: Error using idlsave
* `990 <https://github.com/ipython/ipython/issues/990>`_: ipython notebook - kernel dies if matplotlib is not installed
* `752 <https://github.com/ipython/ipython/issues/752>`_: Matplotlib figures showed only once in notebook
* `54 <https://github.com/ipython/ipython/issues/54>`_: Exception hook should be optional for embedding IPython in GUIs
* `918 <https://github.com/ipython/ipython/issues/918>`_: IPython.frontend tests fail without tornado
* `986 <https://github.com/ipython/ipython/issues/986>`_: Views created with c.direct_view() fail
* `697 <https://github.com/ipython/ipython/issues/697>`_: Filter out from %who names loaded at initialization time
* `932 <https://github.com/ipython/ipython/issues/932>`_: IPython 0.11 quickref card has superfluous "%recall and"
* `982 <https://github.com/ipython/ipython/issues/982>`_: png files with executable permissions
* `914 <https://github.com/ipython/ipython/issues/914>`_: Simpler system for running code after InteractiveShell is initialised
* `911 <https://github.com/ipython/ipython/issues/911>`_: ipython crashes on startup if readline is missing
* `971 <https://github.com/ipython/ipython/issues/971>`_: bookmarks created in 0.11 are corrupt in 0.12
* `974 <https://github.com/ipython/ipython/issues/974>`_: object feature tab-completion crash
* `939 <https://github.com/ipython/ipython/issues/939>`_: ZMQShell always uses default profile
* `946 <https://github.com/ipython/ipython/issues/946>`_: Multi-tab Close action should offer option to leave all kernels alone
* `949 <https://github.com/ipython/ipython/issues/949>`_: Test suite must not require any manual interaction
* `643 <https://github.com/ipython/ipython/issues/643>`_: enable gui eventloop integration in ipkernel
* `965 <https://github.com/ipython/ipython/issues/965>`_: ipython is crashed without launch.(python3.2)
* `958 <https://github.com/ipython/ipython/issues/958>`_: Can't use os X clipboard on with qtconsole
* `962 <https://github.com/ipython/ipython/issues/962>`_: Don't require tornado in the tests
* `960 <https://github.com/ipython/ipython/issues/960>`_: crash on syntax error on Windows XP
* `934 <https://github.com/ipython/ipython/issues/934>`_: The latest ipython branch doesn't work in Chrome
* `870 <https://github.com/ipython/ipython/issues/870>`_: zmq version detection
* `943 <https://github.com/ipython/ipython/issues/943>`_: HISTIGNORE for IPython
* `947 <https://github.com/ipython/ipython/issues/947>`_: qtconsole segfaults at startup
* `903 <https://github.com/ipython/ipython/issues/903>`_: Expose a magic to control config of the inline pylab backend
* `908 <https://github.com/ipython/ipython/issues/908>`_: bad user config shouldn't crash IPython
* `935 <https://github.com/ipython/ipython/issues/935>`_: Typing `break` causes IPython to crash.
* `869 <https://github.com/ipython/ipython/issues/869>`_: Tab completion of `~/` shows no output post 0.10.x
* `904 <https://github.com/ipython/ipython/issues/904>`_: whos under pypy1.6
* `773 <https://github.com/ipython/ipython/issues/773>`_: check_security_dir() and check_pid_dir() fail on network filesystem
* `915 <https://github.com/ipython/ipython/issues/915>`_: OS X Lion Terminal.app line wrap problem
* `886 <https://github.com/ipython/ipython/issues/886>`_: Notebook kernel crash when specifying --notebook-dir on commandline
* `636 <https://github.com/ipython/ipython/issues/636>`_: debugger.py: pydb broken
* `808 <https://github.com/ipython/ipython/issues/808>`_: Ctrl+C during %reset confirm message crash Qtconsole
* `927 <https://github.com/ipython/ipython/issues/927>`_: Using return outside a function crashes ipython
* `919 <https://github.com/ipython/ipython/issues/919>`_: Pop-up segfault when moving cursor out of qtconsole window
* `181 <https://github.com/ipython/ipython/issues/181>`_: cls command does not work on windows
* `917 <https://github.com/ipython/ipython/issues/917>`_: documentation typos
* `818 <https://github.com/ipython/ipython/issues/818>`_: %run does not work with non-ascii characeters in path
* `907 <https://github.com/ipython/ipython/issues/907>`_: Errors in custom completer functions can crash IPython
* `867 <https://github.com/ipython/ipython/issues/867>`_: doc: notebook password authentication howto
* `211 <https://github.com/ipython/ipython/issues/211>`_: paste command not working
* `900 <https://github.com/ipython/ipython/issues/900>`_: Tab key should insert 4 spaces in qt console
* `513 <https://github.com/ipython/ipython/issues/513>`_: [Qt console] cannot insert new lines into console functions using tab
* `906 <https://github.com/ipython/ipython/issues/906>`_: qtconsoleapp 'parse_command_line' doen't like --existing anymore
* `638 <https://github.com/ipython/ipython/issues/638>`_: Qt console --pylab=inline and getfigs(), etc.
* `710 <https://github.com/ipython/ipython/issues/710>`_: unwanted unicode passed to args
* `436 <https://github.com/ipython/ipython/issues/436>`_: Users should see tooltips for all buttons in the notebook UI
* `207 <https://github.com/ipython/ipython/issues/207>`_: ipython crashes if atexit handler raises exception
* `692 <https://github.com/ipython/ipython/issues/692>`_: use of Tracer() when debugging works but gives error messages
* `690 <https://github.com/ipython/ipython/issues/690>`_: debugger does not print error message by default in 0.11
* `571 <https://github.com/ipython/ipython/issues/571>`_: history of multiline entries
* `749 <https://github.com/ipython/ipython/issues/749>`_: IPython.parallel test failure under Windows 7 and XP
* `890 <https://github.com/ipython/ipython/issues/890>`_: ipclusterapp.py - helep
* `885 <https://github.com/ipython/ipython/issues/885>`_: `ws-hostname` alias not recognized by notebook
* `881 <https://github.com/ipython/ipython/issues/881>`_: Missing manual.pdf?
* `744 <https://github.com/ipython/ipython/issues/744>`_: cannot create notebook in offline mode if mathjax not installed
* `865 <https://github.com/ipython/ipython/issues/865>`_: Make tracebacks from %paste show the code
* `535 <https://github.com/ipython/ipython/issues/535>`_: exception unicode handling in %run is faulty in qtconsole
* `817 <https://github.com/ipython/ipython/issues/817>`_: IPython crashed
* `799 <https://github.com/ipython/ipython/issues/799>`_: %edit magic not working on windows xp in qtconsole
* `732 <https://github.com/ipython/ipython/issues/732>`_: QTConsole wrongly promotes the index of the input line on which user presses Enter
* `662 <https://github.com/ipython/ipython/issues/662>`_: ipython test failures on Mac OS X Lion
* `650 <https://github.com/ipython/ipython/issues/650>`_: Handle bad config files better
* `829 <https://github.com/ipython/ipython/issues/829>`_: We should not insert new lines after all print statements in the notebook
* `874 <https://github.com/ipython/ipython/issues/874>`_: ipython-qtconsole: pyzmq Version Comparison
* `640 <https://github.com/ipython/ipython/issues/640>`_: matplotlib macosx windows don't respond in qtconsole
* `624 <https://github.com/ipython/ipython/issues/624>`_: ipython intermittently segfaults when figure is closed (Mac OS X)
* `871 <https://github.com/ipython/ipython/issues/871>`_: Notebook crashes if a profile is used
* `56 <https://github.com/ipython/ipython/issues/56>`_: Have %cpaste accept also Ctrl-D as a termination marker
* `849 <https://github.com/ipython/ipython/issues/849>`_: Command line options to not override profile options
* `806 <https://github.com/ipython/ipython/issues/806>`_: Provide single-port connection to kernels
* `691 <https://github.com/ipython/ipython/issues/691>`_: [wishlist] Automatically find existing kernel
* `688 <https://github.com/ipython/ipython/issues/688>`_: local security vulnerability: all ports visible to any local user.
* `866 <https://github.com/ipython/ipython/issues/866>`_: DistributionNotFound on running ipython 0.11 on Windows XP x86
* `673 <https://github.com/ipython/ipython/issues/673>`_: raw_input appears to be round-robin for qtconsole
* `863 <https://github.com/ipython/ipython/issues/863>`_: Graceful degradation when home directory not writable
* `800 <https://github.com/ipython/ipython/issues/800>`_: Timing scripts with run -t -N <N> fails on report output
* `858 <https://github.com/ipython/ipython/issues/858>`_: Typing 'continue' makes ipython0.11 crash
* `840 <https://github.com/ipython/ipython/issues/840>`_: all processes run on one CPU core
* `843 <https://github.com/ipython/ipython/issues/843>`_: "import braces" crashes ipython
* `836 <https://github.com/ipython/ipython/issues/836>`_: Strange Output after IPython Install
* `839 <https://github.com/ipython/ipython/issues/839>`_: Qtconsole segfaults when mouse exits window with active tooltip
* `827 <https://github.com/ipython/ipython/issues/827>`_: Add support for checking several limits before running task on engine
* `826 <https://github.com/ipython/ipython/issues/826>`_: Add support for creation of parallel task when no engine is running
* `832 <https://github.com/ipython/ipython/issues/832>`_: Improve error message for %logstop
* `831 <https://github.com/ipython/ipython/issues/831>`_: %logstart in read-only directory forbid any further command
* `814 <https://github.com/ipython/ipython/issues/814>`_: ipython does not start -- DistributionNotFound
* `794 <https://github.com/ipython/ipython/issues/794>`_: Allow >1 controller per profile
* `820 <https://github.com/ipython/ipython/issues/820>`_: Tab Completion feature
* `812 <https://github.com/ipython/ipython/issues/812>`_: Qt console crashes on Ubuntu 11.10
* `816 <https://github.com/ipython/ipython/issues/816>`_: Import error using Python 2.7 and dateutil2.0 No module named _thread
* `756 <https://github.com/ipython/ipython/issues/756>`_: qtconsole Windows fails to print error message for '%run nonexistent_file'
* `651 <https://github.com/ipython/ipython/issues/651>`_: Completion doesn't work on element of a list
* `617 <https://github.com/ipython/ipython/issues/617>`_: [qtconsole] %hist doesn't show anything in qtconsole
* `786 <https://github.com/ipython/ipython/issues/786>`_: from __future__ import unicode_literals does not work
* `779 <https://github.com/ipython/ipython/issues/779>`_: Using irunner from virtual evn uses systemwide ipython
* `768 <https://github.com/ipython/ipython/issues/768>`_: codepage handling of output from scripts and shellcommands are not handled properly by qtconsole
* `785 <https://github.com/ipython/ipython/issues/785>`_: Don't strip leading whitespace in repr() in notebook
* `737 <https://github.com/ipython/ipython/issues/737>`_: in pickleshare.py line52 should be "if not os.path.isdir(self.root):"?
* `738 <https://github.com/ipython/ipython/issues/738>`_: in ipthon_win_post_install.py line 38
* `777 <https://github.com/ipython/ipython/issues/777>`_: print(…, sep=…) raises SyntaxError
* `728 <https://github.com/ipython/ipython/issues/728>`_: ipcontroller crash with MPI
* `780 <https://github.com/ipython/ipython/issues/780>`_: qtconsole Out value prints before the print statements that precede it
* `632 <https://github.com/ipython/ipython/issues/632>`_: IPython Crash Report (0.10.2)
* `253 <https://github.com/ipython/ipython/issues/253>`_: Unable to install ipython on windows
* `80 <https://github.com/ipython/ipython/issues/80>`_: Split IPClusterApp into multiple Application subclasses for each subcommand
* `34 <https://github.com/ipython/ipython/issues/34>`_: non-blocking pendingResult partial results
* `739 <https://github.com/ipython/ipython/issues/739>`_: Tests fail if tornado not installed
* `719 <https://github.com/ipython/ipython/issues/719>`_: Better support Pypy
* `667 <https://github.com/ipython/ipython/issues/667>`_: qtconsole problem with default pylab profile
* `661 <https://github.com/ipython/ipython/issues/661>`_: ipythonrc referenced in magic command in 0.11
* `665 <https://github.com/ipython/ipython/issues/665>`_: Source introspection with ?? is broken
* `724 <https://github.com/ipython/ipython/issues/724>`_: crash - ipython qtconsole, %quickref
* `655 <https://github.com/ipython/ipython/issues/655>`_: ipython with qtconsole crashes
* `593 <https://github.com/ipython/ipython/issues/593>`_: HTML Notebook Prompt can be deleted . . .
* `563 <https://github.com/ipython/ipython/issues/563>`_: use argparse instead of kvloader for flags&aliases
* `751 <https://github.com/ipython/ipython/issues/751>`_: Tornado version greater than 2.0 needed for firefox 6
* `720 <https://github.com/ipython/ipython/issues/720>`_: Crash report when importing easter egg
* `740 <https://github.com/ipython/ipython/issues/740>`_: Ctrl-Enter clears line in notebook
* `772 <https://github.com/ipython/ipython/issues/772>`_: ipengine fails on Windows with "XXX lineno: 355, opcode: 0"
* `771 <https://github.com/ipython/ipython/issues/771>`_: Add python 3 tag to setup.py
* `767 <https://github.com/ipython/ipython/issues/767>`_: non-ascii in __doc__ string crashes qtconsole kernel when showing tooltip
* `733 <https://github.com/ipython/ipython/issues/733>`_: In Windows, %run fails to strip quotes from filename
* `721 <https://github.com/ipython/ipython/issues/721>`_: no completion in emacs by ipython(ipython.el)
* `669 <https://github.com/ipython/ipython/issues/669>`_: Do not accept an ipython_dir that's not writeable
* `711 <https://github.com/ipython/ipython/issues/711>`_: segfault on mac os x
* `500 <https://github.com/ipython/ipython/issues/500>`_: "RuntimeError: Cannot change input buffer during execution" in console_widget.py
* `707 <https://github.com/ipython/ipython/issues/707>`_: Copy and paste keyboard shortcuts do not work in Qt Console on OS X
* `478 <https://github.com/ipython/ipython/issues/478>`_: PyZMQ's use of memoryviews breaks reconstruction of numpy arrays
* `694 <https://github.com/ipython/ipython/issues/694>`_: Turning off callout tips in qtconsole
* `704 <https://github.com/ipython/ipython/issues/704>`_: return kills IPython
* `442 <https://github.com/ipython/ipython/issues/442>`_: Users should have intelligent autoindenting in the notebook
* `615 <https://github.com/ipython/ipython/issues/615>`_: Wireframe and implement a project dashboard page
* `614 <https://github.com/ipython/ipython/issues/614>`_: Wireframe and implement a notebook dashboard page
* `606 <https://github.com/ipython/ipython/issues/606>`_: Users should be able to use the notebook to import/export a notebook to .py or .rst
* `604 <https://github.com/ipython/ipython/issues/604>`_: A user should be able to leave a kernel running in the notebook and reconnect
* `298 <https://github.com/ipython/ipython/issues/298>`_: Users should be able to save a notebook and then later reload it
* `649 <https://github.com/ipython/ipython/issues/649>`_: ipython qtconsole (v0.11): setting "c.IPythonWidget.in_prompt = '>>> ' crashes
* `672 <https://github.com/ipython/ipython/issues/672>`_: What happened to Exit?
* `658 <https://github.com/ipython/ipython/issues/658>`_: Put the InteractiveShellApp section first in the auto-generated config files
* `656 <https://github.com/ipython/ipython/issues/656>`_: [suggestion] dependency checking for pyqt for  Windows installer
* `654 <https://github.com/ipython/ipython/issues/654>`_: broken documentation link on download page
* `653 <https://github.com/ipython/ipython/issues/653>`_: Test failures in IPython.parallel
=============
 0.10 series
=============

Release 0.10.2
==============

IPython 0.10.2 was released April 9, 2011.  This is a minor bugfix release that
preserves backward compatibility.  At this point, all IPython development
resources are focused on the 0.11 series that includes a complete architectural
restructuring of the project as well as many new capabilities, so this is
likely to be the last release of the 0.10.x series.  We have tried to fix all
major bugs in this series so that it remains a viable platform for those not
ready yet to transition to the 0.11 and newer codebase (since that will require
some porting effort, as a number of APIs have changed).

Thus, we are not opening a 0.10.3 active development branch yet, but if the
user community requires new patches and is willing to maintain/release such a
branch, we'll be happy to host it on the IPython github repositories.

Highlights of this release:

- The main one is the closing of github ticket #185, a major regression we had
  in 0.10.1 where pylab mode with GTK (or gthread) was not working correctly,
  hence plots were blocking with GTK.  Since this is the default matplotlib
  backend on Unix systems, this was a major annoyance for many users.  Many
  thanks to Paul Ivanov for helping resolve this issue.
  
- Fix IOError bug on Windows when used with -gthread.
- Work robustly if $HOME is missing from environment.
- Better POSIX support in ssh scripts (remove bash-specific idioms).
- Improved support for non-ascii characters in log files.
- Work correctly in environments where GTK can be imported but not started
  (such as a linux text console without X11).
  
For this release we merged 24 commits, contributed by the following people
(please let us know if we omitted your name and we'll gladly fix this in the
notes for the future):

* Fernando Perez
* MinRK
* Paul Ivanov
* Pieter Cristiaan de Groot
* TvrtkoM

Release 0.10.1
==============

IPython 0.10.1 was released October 11, 2010, over a year after version 0.10.
This is mostly a bugfix release, since after version 0.10 was released, the
development team's energy has been focused on the 0.11 series.  We have
nonetheless tried to backport what fixes we could into 0.10.1, as it remains
the stable series that many users have in production systems they rely on.

Since the 0.11 series changes many APIs in backwards-incompatible ways, we are
willing to continue maintaining the 0.10.x series.  We don't really have time
to actively write new code for 0.10.x, but we are happy to accept patches and
pull requests on the IPython `github site`_.  If sufficient contributions are
made that improve 0.10.1, we will roll them into future releases.  For this
purpose, we will have a branch called 0.10.2 on github, on which you can base
your contributions.

.. _github site: http://github.com/ipython

For this release, we applied approximately 60 commits totaling a diff of over
7000 lines::

    (0.10.1)amirbar[dist]> git diff --oneline rel-0.10.. | wc -l
    7296

Highlights of this release:

- The only significant new feature is that IPython's parallel computing
  machinery now supports natively the Sun Grid Engine and LSF schedulers.  This
  work was a joint contribution from Justin Riley, Satra Ghosh and Matthieu
  Brucher, who put a lot of work into it.  We also improved traceback handling
  in remote tasks, as well as providing better control for remote task IDs.

- New IPython Sphinx directive contributed by John Hunter.  You can use this
  directive to mark blocks in reStructuredText documents as containing IPython
  syntax (including figures) and the will be executed during the build:

  .. sourcecode:: ipython

      In [2]: plt.figure()  # ensure a fresh figure

      @savefig psimple.png width=4in
      In [3]: plt.plot([1,2,3])
      Out[3]: [<matplotlib.lines.Line2D object at 0x9b74d8c>]

- Various fixes to the standalone ipython-wx application.

- We now ship internally the excellent argparse library, graciously licensed
  under BSD terms by Steven Bethard.  Now (2010) that argparse has become part
  of Python 2.7 this will be less of an issue, but Steven's relicensing allowed
  us to start updating IPython to using argparse well before Python 2.7.  Many
  thanks!

- Robustness improvements so that IPython doesn't crash if the readline library
  is absent (though obviously a lot of functionality that requires readline
  will not be available).

- Improvements to tab completion in Emacs with Python 2.6.

- Logging now supports timestamps (see ``%logstart?`` for full details).

- A long-standing and quite annoying bug where parentheses would be added to
  ``print`` statements, under Python 2.5 and 2.6, was finally fixed.

- Improved handling of libreadline on Apple OSX.

- Fix ``reload`` method of IPython demos, which was broken.

- Fixes for the ipipe/ibrowse system on OSX.

- Fixes for Zope profile.

- Fix %timeit reporting when the time is longer than 1000s.

- Avoid lockups with ? or ?? in SunOS, due to a bug in termios.

- The usual assortment of miscellaneous bug fixes and small improvements.

The following people contributed to this release (please let us know if we
omitted your name and we'll gladly fix this in the notes for the future):

* Beni Cherniavsky
* Boyd Waters.
* David Warde-Farley
* Fernando Perez
* Gökhan Sever
* John Hunter
* Justin Riley
* Kiorky
* Laurent Dufrechou
* Mark E. Smith
* Matthieu Brucher
* Satrajit Ghosh
* Sebastian Busch
* Václav Šmilauer

Release 0.10
============

This release brings months of slow but steady development, and will be the last
before a major restructuring and cleanup of IPython's internals that is already
under way.  For this reason, we hope that 0.10 will be a stable and robust
release so that while users adapt to some of the API changes that will come
with the refactoring that will become IPython 0.11, they can safely use 0.10 in
all existing projects with minimal changes (if any).

IPython 0.10 is now a medium-sized project, with roughly (as reported by David
Wheeler's :command:`sloccount` utility) 40750 lines of Python code, and a diff
between 0.9.1 and this release that contains almost 28000 lines of code and
documentation.  Our documentation, in PDF format, is a 495-page long PDF
document (also available in HTML format, both generated from the same sources).

Many users and developers contributed code, features, bug reports and ideas to
this release.  Please do not hesitate in contacting us if we've failed to
acknowledge your contribution here.  In particular, for this release we have
contribution from the following people, a mix of new and regular names (in
alphabetical order by first name):

* Alexander Clausen: fix #341726.
* Brian Granger: lots of work everywhere (features, bug fixes, etc).
* Daniel Ashbrook: bug report on MemoryError during compilation, now fixed.
* Darren Dale: improvements to documentation build system, feedback, design
  ideas.
* Fernando Perez: various places.
* Gaël Varoquaux: core code, ipythonx GUI, design discussions, etc. Lots...
* John Hunter: suggestions, bug fixes, feedback.
* Jorgen Stenarson: work on many fronts, tests, fixes, win32 support, etc.
* Laurent Dufréchou: many improvements to ipython-wx standalone app.
* Lukasz Pankowski: prefilter, `%edit`, demo improvements.
* Matt Foster: TextMate support in `%edit`.
* Nathaniel Smith: fix #237073.
* Pauli Virtanen: fixes and improvements to extensions, documentation.
* Prabhu Ramachandran: improvements to `%timeit`.
* Robert Kern: several extensions.
* Sameer D'Costa: help on critical bug #269966.
* Stephan Peijnik: feedback on Debian compliance and many man pages.
* Steven Bethard: we are now shipping his :mod:`argparse` module.
* Tom Fetherston: many improvements to :mod:`IPython.demo` module.
* Ville Vainio: lots of work everywhere (features, bug fixes, etc).
* Vishal Vasta: ssh support in ipcluster.
* Walter Doerwald: work on the :mod:`IPython.ipipe` system.

Below we give an overview of new features, bug fixes and backwards-incompatible
changes.  For a detailed account of every change made, feel free to view the
project log with :command:`bzr log`.

New features
------------

* New `%paste` magic automatically extracts current contents of clipboard and
  pastes it directly, while correctly handling code that is indented or
  prepended with `>>>` or `...` python prompt markers.  A very useful new
  feature contributed by Robert Kern.

* IPython 'demos', created with the :mod:`IPython.demo` module, can now be
  created from files on disk or strings in memory.  Other fixes and
  improvements to the demo system, by Tom Fetherston.

* Added :func:`find_cmd` function to :mod:`IPython.platutils` module, to find
  commands in a cross-platform manner.

* Many improvements and fixes to Gaël Varoquaux's :command:`ipythonx`, a
  WX-based lightweight IPython instance that can be easily embedded in other WX
  applications.  These improvements have made it possible to now have an
  embedded IPython in Mayavi and other tools.

* :class:`MultiengineClient` objects now have a :meth:`benchmark` method.

* The manual now includes a full set of auto-generated API documents from the
  code sources, using Sphinx and some of our own support code.  We are now
  using the `Numpy Documentation Standard`_  for all docstrings, and we have
  tried to update as many existing ones as possible to this format.

* The new :mod:`IPython.Extensions.ipy_pretty` extension by Robert Kern
  provides configurable pretty-printing.

* Many improvements to the :command:`ipython-wx` standalone WX-based IPython
  application by Laurent Dufréchou.  It can optionally run in a thread, and
  this can be toggled at runtime (allowing the loading of Matplotlib in a
  running session without ill effects).

* IPython includes a copy of Steven Bethard's argparse_ in the
  :mod:`IPython.external` package, so we can use it internally and it is also
  available to any IPython user.  By installing it in this manner, we ensure
  zero conflicts with any system-wide installation you may already have while
  minimizing external dependencies for new users.  In IPython 0.10, We ship
  argparse version 1.0.

* An improved and much more robust test suite, that runs groups of tests in
  separate subprocesses using either Nose or Twisted's :command:`trial` runner
  to ensure proper management of Twisted-using code.  The test suite degrades
  gracefully if optional dependencies are not available, so that the
  :command:`iptest` command can be run with only Nose installed and nothing
  else.  We also have more and cleaner test decorators to better select tests
  depending on runtime conditions, do setup/teardown, etc.

* The new ipcluster now has a fully working ssh mode that should work on
  Linux, Unix and OS X.  Thanks to Vishal Vatsa for implementing this!

* The wonderful TextMate editor can now be used with %edit on OS X.  Thanks
  to Matt Foster for this patch.

* The documentation regarding parallel uses of IPython, including MPI and PBS,
  has been significantly updated and improved.

* The developer guidelines in the documentation have been updated to explain
  our workflow using :command:`bzr` and Launchpad.
  
* Fully refactored :command:`ipcluster` command line program for starting
  IPython clusters.  This new version is a complete rewrite and 1) is fully
  cross platform (we now use Twisted's process management), 2) has much
  improved performance, 3) uses subcommands for different types of clusters, 4)
  uses argparse for parsing command line options, 5) has better support for
  starting clusters using :command:`mpirun`, 6) has experimental support for
  starting engines using PBS.  It can also reuse FURL files, by appropriately
  passing options to its subcommands.  However, this new version of ipcluster
  should be considered a technology preview.  We plan on changing the API in
  significant ways before it is final.

* Full description of the security model added to the docs.

* cd completer: show bookmarks if no other completions are available.

* sh profile: easy way to give 'title' to prompt: assign to variable
  '_prompt_title'. It looks like this::
      
        [~]|1> _prompt_title = 'sudo!'
        sudo![~]|2>

* %edit: If you do '%edit pasted_block', pasted_block variable gets updated
  with new data (so repeated editing makes sense)

.. _Numpy Documentation Standard: https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt#docstring-standard

.. _argparse: http://code.google.com/p/argparse/

Bug fixes
---------

* Fix #368719, removed top-level debian/ directory to make the job of Debian
  packagers easier.
  
* Fix #291143 by including man pages contributed by Stephan Peijnik from the
  Debian project.

* Fix #358202, effectively a race condition, by properly synchronizing file
  creation at cluster startup time.

* `%timeit` now handles correctly functions that take a long time to execute
  even the first time, by not repeating them.

* Fix #239054, releasing of references after exiting.

* Fix #341726, thanks to Alexander Clausen.

* Fix #269966.  This long-standing and very difficult bug (which is actually a
  problem in Python itself) meant long-running sessions would inevitably grow
  in memory size, often with catastrophic consequences if users had large
  objects in their scripts.  Now, using `%run` repeatedly should not cause any
  memory leaks.  Special thanks to John Hunter and Sameer D'Costa for their
  help with this bug.

* Fix #295371, bug in `%history`.

* Improved support for py2exe.

* Fix #270856: IPython hangs with PyGTK

* Fix #270998: A magic with no docstring breaks the '%magic magic'

* fix #271684: -c startup commands screw up raw vs. native history

* Numerous bugs on Windows with the new ipcluster have been fixed.

* The ipengine and ipcontroller scripts now handle missing furl files
  more gracefully by giving better error messages.

* %rehashx: Aliases no longer contain dots. python3.0 binary
  will create alias python30. Fixes:
  #259716 "commands with dots in them don't work"

* %cpaste: %cpaste -r repeats the last pasted block.
  The block is assigned to pasted_block even if code
  raises exception.

* Bug #274067 'The code in get_home_dir is broken for py2exe' was
  fixed.

* Many other small bug fixes not listed here by number (see the bzr log for
  more info).
  
Backwards incompatible changes
------------------------------

* `ipykit` and related files were unmaintained and have been removed.

* The :func:`IPython.genutils.doctest_reload` does not actually call
  `reload(doctest)` anymore, as this was causing many problems with the test
  suite.  It still resets `doctest.master` to None.

* While we have not deliberately broken Python 2.4 compatibility, only minor
  testing was done with Python 2.4, while 2.5 and 2.6 were fully tested.  But
  if you encounter problems with 2.4, please do report them as bugs.

* The :command:`ipcluster` now requires a mode argument; for example to start a
  cluster on the local machine with 4 engines, you must now type::

    $ ipcluster local -n 4

* The controller now has a ``-r`` flag that needs to be used if you want to
  reuse existing furl files.  Otherwise they are deleted (the default).

* Remove ipy_leo.py. You can use :command:`easy_install ipython-extension` to
  get it.  (done to decouple it from ipython release cycle)

Issues closed in the 7.x development cycle
==========================================

Stats are not collected after version 7.17, all contribution will show up as part of the 8.0 release.

Issues closed in 7.17
---------------------

GitHub stats for 2020/06/26 - 2020/07/31 (tag: 7.16.1)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 2 issues and merged 19 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.17>`__

The following 3 authors contributed 31 commits.

* Maor Kleinberger
* Matthias Bussonnier
* Quentin Peter



Issues closed in 7.16
---------------------

GitHub stats for 2020/05/29 - 2020/06/26 (tag: 7.15.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 0 issues and merged 18 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.16>`__

The following 7 authors contributed 22 commits.

* Benjamin Ragan-Kelley
* dalthviz
* Frank Tobia
* Matthias Bussonnier
* palewire
* Paul McCarthy
* Talley Lambert


Issues closed in 7.15
---------------------

GitHub stats for 2020/05/01 - 2020/05/29 (tag: 7.14.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 1 issues and merged 29 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.15>`__

The following 6 authors contributed 31 commits.

* Blake Griffin
* Inception95
* Marcio Mazza
* Matthias Bussonnier
* Talley Lambert
* Thomas

Issues closed in 7.14
---------------------

GitHub stats for 2020/02/29 - 2020/05/01 (tag: 7.13.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 0 issues and merged 30 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.14>`__

The following 10 authors contributed 47 commits.

* Eric Wieser
* foobarbyte
* Ian Castleden
* Itamar Turner-Trauring
* Lumir Balhar
* Markus Wageringel
* Matthias Bussonnier
* Matthieu Ancellin
* Quentin Peter
* Theo Ouzhinski

Issues closed in 7.13
---------------------


GitHub stats for 2020/02/01 - 2020/02/28 (tag: 7.12.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 1 issues and merged 24 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.13>`__

The following 12 authors contributed 108 commits.

* Alex Hall
* Augusto
* Coon, Ethan T
* Daniel Hahler
* Inception95
* Itamar Turner-Trauring
* Jonas Haag
* Jonathan Slenders
* linar-jether
* Matthias Bussonnier
* Nathan Goldbaum
* Terry Davis

Issues closed in 7.12
---------------------

GitHub stats for 2020/01/01 - 2020/01/31 (tag: 7.11.1)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 2 issues and merged 14 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.12>`__

The following 11 authors contributed 48 commits.

* Augusto
* Eric Wieser
* Jeff Potter
* Mark E. Haase
* Matthias Bussonnier
* ossdev07
* ras44
* takuya fujiwara
* Terry Davis
* Thomas A Caswell
* yangyang

Issues closed in 7.11
---------------------

GitHub stats for 2019/12/01 - 2019/12/27 (tag: 7.10.1)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 4 issues and merged 36 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.11>`__

The following 16 authors contributed 114 commits.

* Augusto
* Benjamin Ragan-Kelley
* Chemss Eddine Ben Hassine
* Danny Hermes
* Dominik Miedziński
* Jonathan Feinberg
* Jonathan Slenders
* Joseph Kahn
* kousik
* Kousik Mitra
* Marc Hernandez Cabot
* Matthias Bussonnier
* Naveen Honest Raj K
* Pratyay Pandey
* Quentin Peter
* takuya fujiwara


Issues closed in 7.10.2
-----------------------


GitHub stats for 2019/12/01 - 2019/12/14 (tag: 7.10.1)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 3 issues and merged 10 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.10.2>`__

The following 3 authors contributed 11 commits.

* Jonathan Slenders
* Joseph Kahn
* Matthias Bussonnier

Issues closed in 7.10.1
-----------------------

GitHub stats for 2019/11/27 - 2019/12/01 (tag: 7.10.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 5 issues and merged 7 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.10.1>`__

The following 2 authors contributed 14 commits.

* Jonathan Slenders
* Matthias Bussonnier

Issues closed in 7.10
---------------------

GitHub stats for 2019/10/25 - 2019/11/27 (tag: 7.9.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 4 issues and merged 22 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.10>`__

The following 15 authors contributed 101 commits.

* anatoly techtonik
* Ben Lewis
* Benjamin Ragan-Kelley
* Gerrit Buss
* grey275
* Gökcen Eraslan
* Jonathan Slenders
* Joris Van den Bossche
* kousik
* Matthias Bussonnier
* Nicholas Bollweg
* Paul McCarthy
* Srinivas Reddy Thatiparthy
* Timo Kaufmann
* Tony Fast

Issues closed in 7.9
--------------------

GitHub stats for 2019/08/30 - 2019/10/25 (tag: 7.8.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 1 issues and merged 9 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.9>`__

The following 8 authors contributed 20 commits.

* Benjamin Ragan-Kelley
* Hugo
* Matthias Bussonnier
* mfh92
* Mohammad Hossein Sekhavat
* Niclas
* Vidar Tonaas Fauske
* Георгий Фролов

Issues closed in 7.8
--------------------

GitHub stats for 2019/07/26 - 2019/08/30 (tag: 7.7.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 1 issues and merged 4 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.8>`__

The following 5 authors contributed 27 commits.

* Dan Allan
* Matthias Bussonnier
* Min ho Kim
* Oscar Gustafsson
* Terry Davis

Issues closed in 7.7
--------------------

GitHub stats for 2019/07/03 - 2019/07/26 (tag: 7.6.1)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 5 issues and merged 9 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.7>`__

The following 8 authors contributed 26 commits.

* Brandon T. Willard
* juanis2112
* lllf
* Matthias Bussonnier
* Min ho Kim
* Oriol (Prodesk)
* Po-Chuan Hsieh
* techassetskris

Issues closed in 7.6
--------------------

GitHub stats for 2019/04/24 - 2019/06/28 (tag: 7.5.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 9 issues and merged 43 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.6>`__

The following 19 authors contributed 144 commits.

* Alok Singh
* Andreas
* Antony Lee
* Daniel Hahler
* Ed OBrien
* Kevin Sheppard
* Luciana da Costa Marques
* Maor Kleinberger
* Matthias Bussonnier
* Miro Hrončok
* Niclas
* Nikita Bezdolniy
* Oriol Abril
* Piers Titus van der Torren
* Pragnya Srinivasan
* Robin Gustafsson
* stonebig
* Thomas A Caswell
* zzzz-qq


Issues closed in 7.5
--------------------

GitHub stats for 2019/03/21 - 2019/04/24 (tag: 7.4.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 2 issues and merged 9 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.5>`__

The following 7 authors contributed 28 commits.

* Akshay Paropkari
* Benjamin Ragan-Kelley
* Ivan Tham
* Matthias Bussonnier
* Nick Tallant
* Sebastian Witowski
* stef-ubuntu


Issues closed in 7.4
--------------------

GitHub stats for 2019/02/18 - 2019/03/21 (tag: 7.3.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 9 issues and merged 20 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.3>`__

The following 23 authors contributed 69 commits.

* anatoly techtonik
* Benjamin Ragan-Kelley
* bnables
* Frédéric Chapoton
* Gabriel Potter
* Ian Bell
* Jake VanderPlas
* Jan S. (Milania1)
* Jesse Widner
* jsnydes
* Kyungdahm Yun
* Laurent Gautier
* Luciana da Costa Marques
* Matan Gover
* Matthias Bussonnier
* memeplex
* Mickaël Schoentgen
* Partha P. Mukherjee
* Philipp A
* Sanyam Agarwal
* Steve Nicholson
* Tony Fast
* Wouter Overmeire


Issues closed in 7.3
--------------------

GitHub stats for 2018/11/30 - 2019/02/18 (tag: 7.2.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 4 issues and merged 20 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.3>`__

The following 17 authors contributed 99 commits.

* anatoly techtonik
* Benjamin Ragan-Kelley
* Gabriel Potter
* Ian Bell
* Jake VanderPlas
* Jan S. (Milania1)
* Jesse Widner
* Kyungdahm Yun
* Laurent Gautier
* Matthias Bussonnier
* memeplex
* Mickaël Schoentgen
* Partha P. Mukherjee
* Philipp A
* Sanyam Agarwal
* Steve Nicholson
* Tony Fast

Issues closed in 7.2
--------------------

GitHub stats for 2018/10/28 - 2018/11/29 (tag: 7.1.1)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 2 issues and merged 18 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.2>`__

The following 16 authors contributed 95 commits.

* Antony Lee
* Benjamin Ragan-Kelley
* CarsonGSmith
* Chris Mentzel
* Christopher Brown
* Dan Allan
* Elliott Morgan Jobson
* is-this-valid
* kd2718
* Kevin Hess
* Martin Bergtholdt
* Matthias Bussonnier
* Nicholas Bollweg
* Pavel Karateev
* Philipp A
* Reuben Morais

Issues closed in 7.1
--------------------

GitHub stats for 2018/09/27 - 2018/10/27 (since tag: 7.0.1)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 31 issues and merged 54 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.1>`__

The following 33 authors contributed 254 commits.

* ammarmallik
* Audrey Dutcher
* Bart Skowron
* Benjamin Ragan-Kelley
* BinaryCrochet
* Chris Barker
* Christopher Moura
* Dedipyaman Das
* Dominic Kuang
* Elyashiv
* Emil Hessman
* felixzhuologist
* hongshaoyang
* Hugo
* kd2718
* kory donati
* Kory Donati
* koryd
* luciana
* luz.paz
* Massimo Santini
* Matthias Bussonnier
* Matthias Geier
* meeseeksdev[bot]
* Michael Penkov
* Mukesh Bhandarkar
* Nguyen Duy Hai
* Roy Wellington Ⅳ
* Sha Liu
* Shao Yang
* Shashank Kumar
* Tony Fast
* wim glenn


Issues closed in 7.0
--------------------

GitHub stats for 2018/07/29 - 2018/09/27 (since tag: 6.5.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 20 issues and merged 76 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A7.0>`__

The following 49 authors contributed 471 commits.

* alphaCTzo7G
* Alyssa Whitwell
* Anatol Ulrich
* apunisal
* Benjamin Ragan-Kelley
* Chaz Reid
* Christoph
* Dale Jung
* Dave Hirschfeld
* dhirschf
* Doug Latornell
* Fernando Perez
* Fred Mitchell
* Gabriel Potter
* gpotter2
* Grant Nestor
* hongshaoyang
* Hugo
* J Forde
* Jonathan Slenders
* Jörg Dietrich
* Kyle Kelley
* luz.paz
* M Pacer
* Matthew R. Scott
* Matthew Seal
* Matthias Bussonnier
* meeseeksdev[bot]
* Michael Käufl
* Olesya Baranova
* oscar6echo
* Paul Ganssle
* Paul Ivanov
* Peter Parente
* prasanth
* Shailyn javier Ortiz jimenez
* Sourav Singh
* Srinivas Reddy Thatiparthy
* Steven Silvester
* stonebig
* Subhendu Ranjan Mishra
* Takafumi Arakaki
* Thomas A Caswell
* Thomas Kluyver
* Todd
* Wei Yen
* Yarko Tymciurak
* Yutao Yuan
* Zi Chong Kao
============
 5.x Series
============

.. _whatsnew580:

IPython 5.8.0
=============

* Update inspecting function/methods for future-proofing. :ghpull:`11139`

.. _whatsnew570:

IPython 5.7
===========

* Fix IPython trying to import non-existing matplotlib backends :ghpull:`11087`
* fix for display hook not publishing object metadata :ghpull:`11101`

.. _whatsnew560:

IPython 5.6
===========

* In Python 3.6 and above, dictionaries preserve the order items were added to
  them. On these versions, IPython will display dictionaries in their native
  order, rather than sorting by the keys (:ghpull:`10958`).
* :class:`~.IPython.display.ProgressBar` can now be used as an iterator
  (:ghpull:`10813`).
* The shell object gains a :meth:`~.InteractiveShell.check_complete` method,
  to allow a smoother transition to new input processing machinery planned for
  IPython 7 (:ghpull:`11044`).
* IPython should start faster, as it no longer looks for all available pygments
  styles on startup (:ghpull:`10859`).

You can see all the PR marked for the `5.6. milestone <https://github.com/ipython/ipython/pulls?utf8=%E2%9C%93&q=is%3Apr+milestone%3A5.6+is%3Aclosed+NOT+%22Backport+PR%22+>`_,
and all the `backport versions <https://github.com/ipython/ipython/pulls?utf8=%E2%9C%93&q=is%3Apr%20milestone%3A5.6%20is%3Aclosed%20%22Backport%20PR%22%20>`__.

.. _whatsnew550:

IPython 5.5
===========

System Wide config
------------------

- IPython now looks for config files in ``{sys.prefix}/etc/ipython``
  for environment-specific configuration.
- Startup files can be found in ``/etc/ipython/startup`` or ``{sys.prefix}/etc/ipython/startup``
  in addition to the profile directory, for system-wide or env-specific startup files.

See :ghpull:`10644`

ProgressBar
-----------


IPython now has built-in support for progressbars::

    In[1]: from IPython.display import ProgressBar
    ...  : pb = ProgressBar(100)
    ...  : pb

    In[2]: pb.progress = 50

    # progress bar in cell 1 updates.

See :ghpull:`10755`


Misc
----

 - Fix ``IPython.core.display:Pretty._repr_pretty_`` had the wrong signature.
   (:ghpull:`10625`)
 - :magic:`timeit` now give a correct ``SyntaxError`` if naked ``return`` used.
   (:ghpull:`10637`)
 - Prepare the ``:ipython:`` directive to be compatible with Sphinx 1.7.
   (:ghpull:`10668`)
 - Make IPython work with OpenSSL in FIPS mode; change hash algorithm of input
   from md5 to sha1. (:ghpull:`10696`)
 - Clear breakpoints before running any script with debugger. (:ghpull:`10699`)
 - Document that :magic:`profile` is deprecated, not to be confused with :magic:`prun`. (:ghpull:`10707`)
 - Limit default number of returned completions to 500. (:ghpull:`10743`)

You can see all the PR marked for the `5.5. milestone <https://github.com/ipython/ipython/pulls?q=is%3Apr%20milestone%3A5.5%20is%3Aclosed%20NOT%20%22Backport%20PR%22>`_,
and all the `backport versions <https://github.com/ipython/ipython/pulls?utf8=%E2%9C%93&q=is%3Apr%20milestone%3A5.5%20is%3Aclosed%20%22Backport%20PR%22%20>`_.

IPython 5.4.1
=============
Released a few hours after 5.4, fix a crash when
``backports.shutil-get-terminal-size`` is not installed. :ghissue:`10629`

.. _whatsnew540:

IPython 5.4
===========

IPython 5.4-LTS is the first release of IPython after the release of the 6.x
series which is Python 3 only. It backports most of the new exposed API
additions made in IPython 6.0 and 6.1 and avoid having to write conditional
logics depending of the version of IPython.

Please upgrade to pip 9 or greater before upgrading IPython. 
Failing to do so on Python 2 may lead to a broken IPython install.

Configurable TerminalInteractiveShell
-------------------------------------

Backported from the 6.x branch as an exceptional new feature. See
:ghpull:`10373` and :ghissue:`10364`

IPython gained a new ``c.TerminalIPythonApp.interactive_shell_class`` option
that allow to customize the class used to start the terminal frontend. This
should allow user to use custom interfaces, like reviving the former readline
interface which is now a separate package not maintained by the core team.

Define ``_repr_mimebundle_``
----------------------------

Object can now define `_repr_mimebundle_` in place of multiple `_repr_*_`
methods and return a full mimebundle. This greatly simplify many implementation
and allow to publish custom mimetypes (like geojson, plotly, dataframes....).
See the ``Custom Display Logic`` example notebook for more informations.

Execution Heuristics
--------------------

The heuristic for execution in the command line interface is now more biased
toward executing for single statement. While in IPython 4.x and before a single
line would be executed when enter is pressed, IPython 5.x would insert a new
line. For single line statement this is not true anymore and if a single line is
valid Python, IPython will execute it regardless of the cursor position. Use
:kbd:`Ctrl-O` to insert a new line. :ghpull:`10489`


Implement Display IDs
---------------------

Implement display id and ability to update a given display. This should greatly
simplify a lot of code by removing the need for widgets and allow other frontend
to implement things like progress-bars.  See :ghpull:`10048`

Display function
----------------

The :func:`display() <IPython.display.display>` function is now available by
default in an IPython session, meaning users can call it on any object to see
their rich representation. This should allow for better interactivity both at
the REPL and in notebook environment.

Scripts and library that rely on display and may be run outside of IPython still
need to import the display function using ``from IPython.display import
display``. See :ghpull:`10596`


Miscs
-----

* ``_mp_main_`` is not reloaded which fixes issues with multiprocessing.
  :ghpull:`10523`
* Use user colorscheme in Pdb as well :ghpull:`10479`
* Faster shutdown. :ghpull:`10408` 
* Fix a crash in reverse search. :ghpull:`10371`
* added ``Completer.backslash_combining_completions`` boolean option to
  deactivate backslash-tab completion that may conflict with windows path. 

IPython 5.3
===========

Released on February 24th, 2017. Remarkable changes and fixes:

* Fix a bug in ``set_next_input`` leading to a crash of terminal IPython.
  :ghpull:`10231`, :ghissue:`10296`, :ghissue:`10229`
* Always wait for editor inputhook for terminal IPython :ghpull:`10239`,
  :ghpull:`10240`
* Disable ``_ipython_display_`` in terminal :ghpull:`10249`, :ghpull:`10274`
* Update terminal colors to be more visible by default on windows
  :ghpull:`10260`, :ghpull:`10238`, :ghissue:`10281`
* Add Ctrl-Z shortcut (suspend) in terminal debugger :ghpull:`10254`,
  :ghissue:`10273`
* Indent on new line by looking at the text before the cursor :ghpull:`10264`,
  :ghpull:`10275`, :ghissue:`9283`
* Update QtEventloop integration to fix some matplotlib integration issues
  :ghpull:`10201`, :ghpull:`10311`, :ghissue:`10201`
* Respect completions display style in terminal debugger :ghpull:`10305`,
  :ghpull:`10313`
* Add a config option ``TerminalInteractiveShell.extra_open_editor_shortcuts``
  to enable extra shortcuts to open the input in an editor. These are :kbd:`v`
  in vi mode, and :kbd:`C-X C-E` in emacs mode (:ghpull:`10330`).
  The :kbd:`F2` shortcut is always enabled.

IPython 5.2.2
=============

* Fix error when starting with ``IPCompleter.limit_to__all__`` configured.

IPython 5.2.1
=============

* Fix tab completion in the debugger. :ghpull:`10223`

IPython 5.2
===========

Released on January 29th, 2017. Remarkable changes and fixes:

* restore IPython's debugger to raise on quit. :ghpull:`10009`
* The configuration value ``c.TerminalInteractiveShell.highlighting_style`` can
  now directly take a class argument for custom color style. :ghpull:`9848`
* Correctly handle matplotlib figures dpi :ghpull:`9868`
* Deprecate ``-e`` flag for the ``%notebook`` magic that had no effects.
  :ghpull:`9872`
* You can now press F2 while typing at a terminal prompt to edit the contents
  in your favourite terminal editor. Set the :envvar:`EDITOR` environment
  variable to pick which editor is used. :ghpull:`9929`
* sdists will now only be ``.tar.gz`` as per upstream PyPI requirements.
  :ghpull:`9925`
* :any:`IPython.core.debugger` have gained a ``set_trace()`` method for
  convenience. :ghpull:`9947`
* The 'smart command mode' added to the debugger in 5.0 was removed, as more
  people preferred the previous behaviour. Therefore, debugger commands such as
  ``c`` will act as debugger commands even when ``c`` is defined as a variable.
  :ghpull:`10050`
* Fixes OS X event loop issues at startup, :ghpull:`10150`
* Deprecate the ``%autoindent`` magic. :ghpull:`10176`
* Emit a :any:`DeprecationWarning` when setting the deprecated
  ``limit_to_all`` option of the completer. :ghpull:`10198`
* The :cellmagic:`capture` magic can now capture the result of a cell (from an
  expression on the last line), as well as printed and displayed output.
  :ghpull:`9851`.


Changes of behavior to :any:`InteractiveShellEmbed`.

:any:`InteractiveShellEmbed` interactive behavior have changed a bit in between
5.1 and 5.2. By default ``%kill_embedded`` magic will prevent further invocation
of the current ``call location`` instead of preventing further invocation of
the current instance creation location. For most use case this will not change
much for you, though previous behavior was confusing and less consistent with
previous IPython versions.

You can now deactivate instances by using ``%kill_embedded --instance`` flag,
(or ``-i`` in short). The ``%kill_embedded`` magic also gained a
``--yes``/``-y`` option which skip confirmation step, and  ``-x``/``--exit``
which also exit the current embedded call without asking for confirmation.

See :ghpull:`10207`.



IPython 5.1
===========

* Broken ``%timeit`` on Python2 due to the use of ``__qualname__``. :ghpull:`9804`
* Restore ``%gui qt`` to create and return a ``QApplication`` if necessary. :ghpull:`9789`
* Don't set terminal title by default. :ghpull:`9801`
* Preserve indentation when inserting newlines with ``Ctrl-O``. :ghpull:`9770`
* Restore completion in debugger. :ghpull:`9785`
* Deprecate ``IPython.core.debugger.Tracer()`` in favor of simpler, newer, APIs. :ghpull:`9731`
* Restore ``NoOpContext`` context manager removed by mistake, and add `DeprecationWarning`. :ghpull:`9765`
* Add option allowing ``Prompt_toolkit`` to use 24bits colors. :ghpull:`9736`
* Fix for closing interactive matplotlib windows on OS X. :ghpull:`9854`
* An embedded interactive shell instance can be used more than once. :ghpull:`9843`
* More robust check for whether IPython is in a terminal. :ghpull:`9833`
* Better pretty-printing of dicts on PyPy. :ghpull:`9827`
* Some coloured output now looks better on dark background command prompts in Windows.
  :ghpull:`9838`
* Improved tab completion of paths on Windows . :ghpull:`9826`
* Fix tkinter event loop integration on Python 2 with ``future`` installed. :ghpull:`9824`
* Restore ``Ctrl-\`` as a shortcut to quit IPython.
* Make ``get_ipython()`` accessible when modules are imported by startup files. :ghpull:`9818`
* Add support for running directories containing a ``__main__.py`` file with the
  ``ipython`` command. :ghpull:`9813`


True Color feature
------------------

``prompt_toolkit`` uses pygments styles for syntax highlighting. By default, the
colors specified in the style are approximated using a standard 256-color
palette. ``prompt_toolkit`` also supports 24bit, a.k.a. "true", a.k.a. 16-million
color escape sequences which enable compatible terminals to display the exact
colors specified instead of an approximation. This true_color option exposes
that capability in prompt_toolkit to the IPython shell.

Here is a good source for the current state of true color support in various
terminal emulators and software projects: https://gist.github.com/XVilka/8346728



IPython 5.0
===========

Released July 7, 2016

New terminal interface
----------------------

IPython 5 features a major upgrade to the terminal interface, bringing live
syntax highlighting as you type, proper multiline editing and multiline paste,
and tab completions that don't clutter up your history.

.. image:: ../_images/ptshell_features.png
    :alt: New terminal interface features
    :align: center
    :target: ../_images/ptshell_features.png

These features are provided by the Python library `prompt_toolkit
<https://python-prompt-toolkit.readthedocs.io/en/stable/>`__, which replaces
``readline`` throughout our terminal interface.

Relying on this pure-Python, cross platform module also makes it simpler to
install IPython. We have removed dependencies on ``pyreadline`` for Windows and
``gnureadline`` for Mac.

Backwards incompatible changes
------------------------------

- The ``%install_ext`` magic function, deprecated since 4.0, has now been deleted.
  You can distribute and install extensions as packages on PyPI.
- Callbacks registered while an event is being handled will now only be called
  for subsequent events; previously they could be called for the current event.
  Similarly, callbacks removed while handling an event *will* always get that
  event. See :ghissue:`9447` and :ghpull:`9453`.
- Integration with pydb has been removed since pydb development has been stopped
  since 2012, and pydb is not installable from PyPI.
- The ``autoedit_syntax`` option has apparently been broken for many years.
  It has been removed.

New terminal interface
~~~~~~~~~~~~~~~~~~~~~~

The overhaul of the terminal interface will probably cause a range of minor
issues for existing users.
This is inevitable for such a significant change, and we've done our best to
minimise these issues.
Some changes that we're aware of, with suggestions on how to handle them:

IPython no longer uses readline configuration (``~/.inputrc``). We hope that
the functionality you want (e.g. vi input mode) will be available by configuring
IPython directly (see :doc:`/config/options/terminal`).
If something's missing, please file an issue.

The ``PromptManager`` class has been removed, and the prompt machinery simplified.
See :ref:`custom_prompts` to customise prompts with the new machinery.

:mod:`IPython.core.debugger` now provides a plainer interface.
:mod:`IPython.terminal.debugger` contains the terminal debugger using
prompt_toolkit.

There are new options to configure the colours used in syntax highlighting.
We have tried to integrate them with our classic  ``--colors`` option and
``%colors`` magic, but there's a mismatch in possibilities, so some configurations
may produce unexpected results. See :ref:`termcolour` for more information.

The new interface is not compatible with Emacs 'inferior-shell' feature. To
continue using this, add the ``--simple-prompt`` flag to the command Emacs
runs. This flag disables most IPython features, relying on Emacs to provide
things like tab completion.

Provisional Changes
-------------------

Provisional changes are experimental functionality that may, or may not, make
it into a future version of IPython, and which API may change without warnings.
Activating these features and using these API are at your own risk, and may have
security implication for your system, especially if used with the Jupyter notebook,

When running via the Jupyter notebook interfaces, or other compatible client,
you can enable rich documentation experimental functionality:

When the ``docrepr`` package is installed setting the boolean flag
``InteractiveShell.sphinxify_docstring`` to ``True``, will process the various
object through sphinx before displaying them (see the ``docrepr`` package
documentation for more information.

You need to also enable the IPython pager display rich HTML representation
using the ``InteractiveShell.enable_html_pager`` boolean configuration option.
As usual you can set these configuration options globally in your configuration
files, alternatively you can turn them on dynamically using the following
snippet:

.. code-block:: python

    ip = get_ipython()
    ip.sphinxify_docstring = True
    ip.enable_html_pager = True


You can test the effect of various combinations of the above configuration in
the Jupyter notebook, with things example like :

.. code-block:: ipython

    import numpy as np
    np.histogram?


This is part of an effort to make Documentation in Python richer and provide in
the long term if possible dynamic examples that can contain math, images,
widgets... As stated above this is nightly experimental feature with a lot of
(fun) problem to solve. We would be happy to get your feedback and expertise on
it.



Deprecated Features
-------------------

Some deprecated features are listed in this section. Don't forget to enable
``DeprecationWarning`` as an error if you are using IPython in a Continuous
Integration setup or in your testing in general:

.. code-block:: python

    import warnings
    warnings.filterwarnings('error', '.*', DeprecationWarning, module='yourmodule.*')


- ``hooks.fix_error_editor`` seems unused and is pending deprecation.
- `IPython/core/excolors.py:ExceptionColors` is  deprecated.
- `IPython.core.InteractiveShell:write()` is deprecated; use `sys.stdout` instead.
- `IPython.core.InteractiveShell:write_err()` is deprecated; use `sys.stderr` instead.
- The `formatter` keyword argument to `Inspector.info` in `IPython.core.oinspec` has no effect.
- The `global_ns` keyword argument of IPython Embed was deprecated, and has no effect. Use `module` keyword argument instead.


Known Issues:
-------------

- ``<Esc>`` Key does not dismiss the completer and does not clear the current
  buffer. This is an on purpose modification due to current technical
  limitation. Cf :ghpull:`9572`. Escape the control character which is used
  for other shortcut, and there is no practical way to distinguish. Use Ctr-G
  or Ctrl-C as an alternative.

- Cannot use ``Shift-Enter`` and ``Ctrl-Enter`` to submit code in terminal. cf
  :ghissue:`9587` and :ghissue:`9401`. In terminal there is no practical way to
  distinguish these key sequences from a normal new line return.

- ``PageUp`` and ``pageDown`` do not move through completion menu.

- Color styles might not adapt to terminal emulator themes. This will need new
  version of Pygments to be released, and can be mitigated with custom themes.
========================================
0.8 series
========================================

Release 0.8.4
=============

This was a quick release to fix an unfortunate bug that slipped into the 0.8.3
release.  The ``--twisted`` option was disabled, as it turned out to be broken
across several platforms.


Release 0.8.3
=============
  
* pydb is now disabled by default (due to %run -d problems). You can enable
  it by passing -pydb command line argument to IPython. Note that setting
  it in config file won't work.

  
Release 0.8.2
=============

* %pushd/%popd behave differently; now "pushd /foo" pushes CURRENT directory 
  and jumps to /foo. The current behaviour is closer to the documented 
  behaviour, and should not trip anyone.

  
Older releases
==============

Changes in earlier releases of IPython are described in the older file
``ChangeLog``.  Please refer to this document for details.

=============
 0.12 Series
=============

Release 0.12.1
==============

IPython 0.12.1 is a bugfix release of 0.12, pulling only bugfixes and minor
cleanup from 0.13, timed for the Ubuntu 12.04 LTS release.

See the :ref:`list of fixed issues <issues_list_012>` for specific backported issues.


Release 0.12
============

IPython 0.12 contains several major new features, as well as a large amount of
bug and regression fixes.  The 0.11 release brought with it a lot of new
functionality and major refactorings of the codebase; by and large this has
proven to be a success as the number of contributions to the project has
increased dramatically, proving that the code is now much more approachable.
But in the refactoring inevitably some bugs were introduced, and we have also
squashed many of those as well as recovered some functionality that had been
temporarily disabled due to the API changes.

The following major new features appear in this version.


An interactive browser-based Notebook with rich media support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A powerful new interface puts IPython in your browser. You can start it with
the command ``ipython notebook``:

.. figure:: ../_images/notebook_specgram.png
    :width: 400px
    :alt: The IPython notebook with embedded text, code, math and figures.
    :align: center
    :target: ../_images/notebook_specgram.png

    The new IPython notebook showing text, mathematical expressions in LaTeX,
    code, results and embedded figures created with Matplotlib.

This new interface maintains all the features of IPython you are used to, as it
is a new client that communicates with the same IPython kernels used by the
terminal and Qt console.  But the web notebook provides for a different
workflow where you can integrate, along with code execution, also text,
mathematical expressions, graphics, video, and virtually any content that a
modern browser is capable of displaying.

You can save your work sessions as documents that retain all these elements and
which can be version controlled, emailed to colleagues or saved as HTML or PDF
files for printing or publishing statically on the web.  The internal storage
format is a JSON file that can be easily manipulated for manual exporting to
other formats.

This Notebook is a major milestone for IPython, as for years we have tried to
build this kind of system.  We were inspired originally by the excellent
implementation in Mathematica, we made a number of attempts using older
technologies in earlier Summer of Code projects in 2005 (both students and
Robert Kern developed early prototypes), and in recent years we have seen the
excellent implementation offered by the `Sage <http://sagemath.org>` system.
But we continued to work on something that would be consistent with the rest of
IPython's design, and it is clear now that the effort was worth it: based on
the ZeroMQ communications architecture introduced in version 0.11, the notebook
can now retain 100% of the features of the real IPython.  But it can also
provide the rich media support and high quality Javascript libraries that were
not available in browsers even one or two years ago (such as high-quality
mathematical rendering or built-in video).

The notebook has too many useful and important features to describe in these
release notes; our documentation now contains a directory called
``examples/notebooks`` with several notebooks that illustrate various aspects
of the system.  You should start by reading those named
``00_notebook_tour.ipynb`` and ``01_notebook_introduction.ipynb`` first, and
then can proceed to read the others in any order you want.

To start the notebook server, go to a directory containing the notebooks you
want to open (or where you want to create new ones) and type::

  ipython notebook

You can see all the relevant options with::

  ipython notebook --help
  ipython notebook --help-all  # even more

and just like the Qt console, you can start the notebook server with pylab
support by using::
  
  ipython notebook --pylab

for floating matplotlib windows or::
  
  ipython notebook --pylab inline

for plotting support with automatically inlined figures.  Note that it is now
possible also to activate pylab support at runtime via ``%pylab``, so you do
not need to make this decision when starting the server.


.. _two_process_console:

Two-process terminal console
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Based on the same architecture as the notebook and the Qt console, we also have
now a terminal-based console that can connect to an external IPython kernel
(the same kernels used by the Qt console or the notebook, in fact).  While this
client behaves almost identically to the usual IPython terminal application,
this capability can be very useful to attach an interactive console to an
existing kernel that was started externally.  It lets you use the interactive
``%debug`` facilities in a notebook, for example (the web browser can't
interact directly with the debugger) or debug a third-party code where you may
have embedded an IPython kernel.

This is also something that we have wanted for a long time, and which is a
culmination (as a team effort) of the work started last year during the 2010
Google Summer of Code project.
  
Tabbed QtConsole
~~~~~~~~~~~~~~~~

The QtConsole now supports starting multiple kernels in tabs, and has a
menubar, so it looks and behaves more like a real application.  Keyboard
enthusiasts can disable the menubar with ctrl-shift-M (:ghpull:`887`).

.. figure:: ../_images/qtconsole_tabbed.png
    :width: 400px
    :alt: Tabbed IPython Qt console with embedded plots and menus.
    :align: center
    :target: ../_images/qtconsole_tabbed.png

    The improved Qt console for IPython, now with tabs to control multiple
    kernels and full menu support.


Full Python 3 compatibility
~~~~~~~~~~~~~~~~~~~~~~~~~~~

IPython can now be installed from a single codebase on Python 2 and
Python 3. The installation process for Python 3 automatically runs 2to3. The
same 'default' profile is now used for Python 2 and 3 (the previous version had
a separate 'python3' profile).

Standalone Kernel
~~~~~~~~~~~~~~~~~

The ``ipython kernel`` subcommand has been added, to allow starting a
standalone kernel, that can be used with various frontends.  You can then later
connect a Qt console or a terminal console to this kernel by typing e.g.::

  ipython qtconsole --existing

if it's the only one running, or by passing explicitly the connection
parameters (printed by the kernel at startup).


PyPy support
~~~~~~~~~~~~

The terminal interface to IPython now runs under `PyPy <http://pypy.org/>`_.
We will continue to monitor PyPy's progress, and hopefully before long at least
we'll be able to also run the notebook.  The Qt console may take longer, as Qt
is a very complex set of bindings to a huge C++ library, and that is currently
the area where PyPy still lags most behind.  But for everyday interactive use
at the terminal, with this release and PyPy 1.7, things seem to work quite well
from our admittedly limited testing.

  
Other important new features
----------------------------

* **SSH Tunnels**: In 0.11, the :mod:`IPython.parallel` Client could tunnel its
  connections to the Controller via ssh. Now, the QtConsole supports ssh tunneling,
  as do parallel engines.

* **relaxed command-line parsing**: 0.11 was released with overly-strict
  command-line parsing, preventing the ability to specify arguments with spaces,
  e.g. ``ipython --pylab qt`` or ``ipython -c "print 'hi'"``. This has
  been fixed, by using argparse. The new parsing is a strict superset of 0.11, so
  any commands in 0.11 should still work in 0.12.

* **HistoryAccessor**: The :class:`~IPython.core.history.HistoryManager` class
  for interacting with your IPython SQLite history database has been split,
  adding a parent :class:`~IPython.core.history.HistoryAccessor` class, so that
  users can write code to access and search their IPython history without being
  in an IPython session (:ghpull:`824`).

* **kernel %gui and %pylab**: The ``%gui`` and ``%pylab`` magics have been
  restored to the IPython kernel (e.g. in the qtconsole or notebook). This
  allows activation of pylab-mode, or eventloop integration after starting the
  kernel, which was unavailable in 0.11.  Unlike in the terminal, this can be
  set only once, and cannot be changed.

* **%config**: A new ``%config`` magic has been added, giving easy access to the
  IPython configuration system at runtime (:ghpull:`923`).

* **Multiline History**: Multiline readline history has been restored to the
  Terminal frontend by default (:ghpull:`838`).

* **%store**: The ``%store`` magic from earlier versions has been updated and
  re-enabled (:ref:`extensions_storemagic`; :ghpull:`1029`). To autorestore
  stored variables on startup, specify ``c.StoreMagic.autorestore = True`` in
  :file:`ipython_config.py`.


Major Bugs fixed
----------------

In this cycle, we have :ref:`closed over 500 issues <issues_list_012>`, but a
few major ones merit special mention:

* Simple configuration errors should no longer crash IPython. In 0.11, errors
  in config files, as well as invalid trait values, could crash IPython. Now,
  such errors are reported, and help is displayed.

* Certain SyntaxErrors no longer crash IPython (e.g. just typing keywords, such
  as ``return``, ``break``, etc.). See :ghissue:`704`.

* IPython path utils, such as :func:`~IPython.utils.path.get_ipython_dir` now
  check for write permissions, so IPython should function on systems where the
  default path resolution might point to a read-only location, such as
  ``HOMESHARE`` on Windows (:ghissue:`669`).

* :func:`raw_input` now works in the kernel when multiple frontends are in
  use. The request will be sent to the frontend that made the request, and an
  exception is raised if that frontend does not support stdin requests
  (e.g. the notebook) (:ghissue:`673`).

* :mod:`zmq` version detection no longer uses simple lexicographical comparison
  to check minimum version, which prevents 0.11 from working with pyzmq-2.1.10
  (:ghpull:`758`).

* A bug in PySide < 1.0.7 caused crashes on OSX when tooltips were shown
  (:ghissue:`711`). these tooltips are now disabled on old PySide
  (:ghpull:`963`).

* IPython no longer crashes when started on recent versions of Python 3 in
  Windows (:ghissue:`737`).

* Instances of classes defined interactively can now be pickled (:ghissue:`29`;
  :ghpull:`648`). Note that pickling saves a reference to the class definition,
  so unpickling the instances will only work where the class has been defined.


Backwards incompatible changes
------------------------------

* IPython connection information is no longer specified via ip/port directly,
  rather via json connection files.  These files are stored in the security
  directory, and enable us to turn on HMAC message authentication by default,
  significantly improving the security of kernels.  Various utility functions
  have been added to :mod:`IPython.lib.kernel`, for easier connecting to existing
  kernels.

* :class:`~IPython.zmq.kernelmanager.KernelManager` now has one ip, and several
  port traits, rather than several ip/port pair ``_addr`` traits. This better
  matches the rest of the code, where the ip cannot not be set separately for
  each channel.

* Custom prompts are now configured using a new class,
  :class:`~IPython.core.prompts.PromptManager`, which has traits for
  :attr:`in_template`, :attr:`in2_template` (the ``...:`` continuation prompt),
  :attr:`out_template` and :attr:`rewrite_template`. This uses Python's string
  formatting system, so you can use ``{time}`` and ``{cwd}``, although we have
  preserved the abbreviations from previous versions, e.g. ``\#`` (prompt number)
  and ``\w`` (working directory). For the list of available fields, refer to the
  source of :file:`IPython/core/prompts.py`.

* The class inheritance of the Launchers in
  :mod:`IPython.parallel.apps.launcher` used by ipcluster has changed, so that
  trait names are more consistent across batch systems. This may require a few
  renames in your config files, if you customized the command-line args for
  launching controllers and engines. The configurable names have also been
  changed to be clearer that they point to class names, and can now be
  specified by name only, rather than requiring the full import path of each
  class, e.g.::

    IPClusterEngines.engine_launcher = 'IPython.parallel.apps.launcher.MPIExecEngineSetLauncher'
    IPClusterStart.controller_launcher = 'IPython.parallel.apps.launcher.SSHControllerLauncher'

  would now be specified as::

    IPClusterEngines.engine_launcher_class = 'MPI'
    IPClusterStart.controller_launcher_class = 'SSH'

  The full path will still work, and is necessary for using custom launchers
  not in IPython's launcher module.
  
  Further, MPIExec launcher names are now prefixed with just MPI, to better match
  other batch launchers, and be generally more intuitive.  The MPIExec names are
  deprecated, but continue to work.

* For embedding a shell, note that the parameters ``user_global_ns`` and
  ``global_ns`` have been deprecated in favour of ``user_module`` and
  ``module`` respsectively.  The new parameters expect a module-like object,
  rather than a namespace dict.  The old parameters remain for backwards
  compatibility, although ``user_global_ns`` is now ignored. The ``user_ns``
  parameter works the same way as before, and calling
  :func:`~IPython.frontend.terminal.embed.embed` with no arguments still works
  as before.


Development summary and credits
-------------------------------

The previous version (IPython 0.11) was released on July 31 2011, so this
release cycle was roughly 4 1/2 months long, we closed a total of 515 issues,
257 pull requests and 258 regular issues (a :ref:`detailed list
<issues_list_012>` is available).

Many users and developers contributed code, features, bug reports and ideas to
this release.  Please do not hesitate in contacting us if we've failed to
acknowledge your contribution here.  In particular, for this release we have
had commits from the following 45 contributors, a mix of new and regular names
(in alphabetical order by first name):

* Alcides <alcides-at-do-not-span-me.com>
* Ben Edwards <bedwards-at-cs.unm.edu>
* Benjamin Ragan-Kelley <benjaminrk-at-gmail.com>
* Benjamin Thyreau <benjamin.thyreau-at-gmail.com>
* Bernardo B. Marques <bernardo.fire-at-gmail.com>
* Bernard Paulus <bprecyclebin-at-gmail.com>
* Bradley M. Froehle <brad.froehle-at-gmail.com>
* Brian E. Granger <ellisonbg-at-gmail.com>
* Christian Boos <cboos-at-bct-technology.com>
* Daniel Velkov <danielv-at-mylife.com>
* Erik Tollerud <erik.tollerud-at-gmail.com>
* Evan Patterson <epatters-at-enthought.com>
* Felix Werner <Felix.Werner-at-kit.edu>
* Fernando Perez <Fernando.Perez-at-berkeley.edu>
* Gabriel <g2p.code-at-gmail.com>
* Grahame Bowland <grahame-at-angrygoats.net>
* Hannes Schulz <schulz-at-ais.uni-bonn.de>
* Jens Hedegaard Nielsen <jenshnielsen-at-gmail.com>
* Jonathan March <jmarch-at-enthought.com>
* Jörgen Stenarson <jorgen.stenarson-at-bostream.nu>
* Julian Taylor <jtaylor.debian-at-googlemail.com>
* Kefu Chai <tchaikov-at-gmail.com>
* macgyver <neil.rabinowitz-at-merton.ox.ac.uk>
* Matt Cottingham <matt.cottingham-at-gmail.com>
* Matthew Brett <matthew.brett-at-gmail.com>
* Matthias BUSSONNIER <bussonniermatthias-at-gmail.com>
* Michael Droettboom <mdboom-at-gmail.com>
* Nicolas Rougier <Nicolas.Rougier-at-inria.fr>
* Olivier Verdier <olivier.verdier-at-gmail.com>
* Omar Andres Zapata Mesa <andresete.chaos-at-gmail.com>
* Pablo Winant <pablo.winant-at-gmail.com>
* Paul Ivanov <pivanov314-at-gmail.com>
* Pauli Virtanen <pav-at-iki.fi>
* Pete Aykroyd <aykroyd-at-gmail.com>
* Prabhu Ramachandran <prabhu-at-enthought.com>
* Puneeth Chaganti <punchagan-at-gmail.com>
* Robert Kern <robert.kern-at-gmail.com>
* Satrajit Ghosh <satra-at-mit.edu>
* Stefan van der Walt <stefan-at-sun.ac.za>
* Szabolcs Horvát <szhorvat-at-gmail.com>
* Thomas Kluyver <takowl-at-gmail.com>
* Thomas Spura <thomas.spura-at-gmail.com>
* Timo Paulssen <timonator-at-perpetuum-immobile.de>
* Valentin Haenel <valentin.haenel-at-gmx.de>
* Yaroslav Halchenko <debian-at-onerussian.com>
   
.. note::

    This list was generated with the output of
    ``git log rel-0.11..HEAD --format='* %aN <%aE>' | sed 's/@/\-at\-/' | sed 's/<>//' | sort -u``
    after some cleanup.  If you should be on this list, please add yourself.
Issues closed in the 6.x development cycle
==========================================

Issues closed in 6.3
--------------------


GitHub stats for 2017/09/15 - 2018/04/02 (tag: 6.2.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 10 issues and merged 50 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A6.3>`__

The following 35 authors contributed 253 commits.

* Anatoly Techtonik
* Antony Lee
* Benjamin Ragan-Kelley
* Corey McCandless
* Craig Citro
* Cristian Ciupitu
* David Cottrell
* David Straub
* Doug Latornell
* Fabio Niephaus
* Gergely Nagy
* Henry Fredrick Schreiner
* Hugo
* Ismael Venegas Castelló
* Ivan Gonzalez
* J Forde
* Jeremy Sikes
* Joris Van den Bossche
* Lesley Cordero
* luzpaz
* madhu94
* Matthew R. Scott
* Matthias Bussonnier
* Matthias Geier
* Olesya Baranova
* Peter Williams
* Rastislav Barlik
* Roshan Rao
* rs2
* Samuel Lelièvre
* Shailyn javier Ortiz jimenez
* Sjoerd de Vries
* Teddy Rendahl
* Thomas A Caswell
* Thomas Kluyver

Issues closed in 6.2
--------------------

GitHub stats for 2017/05/31 - 2017/09/15 (tag: 6.1.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 3 issues and merged 37 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A6.2+>`__

The following 32 authors contributed 196 commits.

* adityausathe
* Antony Lee
* Benjamin Ragan-Kelley
* Carl Smith
* Eren Halici
* Erich Spaker
* Grant Nestor
* Jean Cruypenynck
* Jeroen Demeyer
* jfbu
* jlstevens
* jus1tin
* Kyle Kelley
* M Pacer
* Marc Richter
* Marius van Niekerk
* Matthias Bussonnier
* mpacer
* Mradul Dubey
* ormung
* pepie34
* Ritesh Kadmawala
* ryan thielke
* Segev Finer
* Srinath
* Srinivas Reddy Thatiparthy
* Steven Maude
* Sudarshan Raghunathan
* Sudarshan Rangarajan
* Thomas A Caswell
* Thomas Ballinger
* Thomas Kluyver


Issues closed in 6.1
--------------------

GitHub stats for 2017/04/19 - 2017/05/30 (tag: 6.0.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 10 issues and merged 43 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A6.1+>`__

The following 26 authors contributed 116 commits.

* Alex Alekseyev
* Benjamin Ragan-Kelley
* Brian E. Granger
* Christopher C. Aycock
* Dave Willmer
* David Bradway
* ICanWaitAndFishAllDay
* Ignat Shining
* Jarrod Janssen
* Joshua Storck
* Luke Pfister
* Matthias Bussonnier
* Matti Remes
* meeseeksdev[bot]
* memeplex
* Ming Zhang
* Nick Weseman
* Paul Ivanov
* Piotr Zielinski
* ryan thielke
* sagnak
* Sang Min Park
* Srinivas Reddy Thatiparthy
* Steve Bartz
* Thomas Kluyver
* Tory Haavik


Issues closed in 6.0
--------------------

GitHub stats for 2017/04/10 - 2017/04/19 (milestone: 6.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 49 issues and merged 145 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A6.0+>`__

The following 34 authors contributed 176 commits.

* Adam Eury
* anantkaushik89
* Antonino Ingargiola
* Benjamin Ragan-Kelley
* Carol Willing
* Chilaka Ramakrishna
* chillaranand
* Denis S. Tereshchenko
* Diego Garcia
* fatData
* Fermi paradox
* fuho
* Grant Nestor
* Ian Rose
* Jeroen Demeyer
* kaushikanant
* Keshav Ramaswamy
* Matteo
* Matthias Bussonnier
* mbyt
* Michael Käufl
* michaelpacer
* Moez Bouhlel
* Pablo Galindo
* Paul Ivanov
* Piotr Przetacznik
* Rounak Banik
* sachet-mittal
* Srinivas Reddy Thatiparthy
* Tamir Bahar
* Thomas Hisch
* Thomas Kluyver
* Utkarsh Upadhyay
* Yuri Numerov
.. _issues_list_3:

Issues closed in the 3.x development cycle
==========================================


Issues closed in 3.2.1
----------------------

GitHub stats for 2015/06/22 - 2015/07/12 (since 3.2)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 1 issue and merged 3 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/milestones/3.2.1>`__

The following 5 authors contributed 9 commits.

* Benjamin Ragan-Kelley
* Matthias Bussonnier
* Nitin Dahyabhai
* Sebastiaan Mathot
* Thomas Kluyver


Issues closed in 3.2
--------------------

GitHub stats for 2015/04/03 - 2015/06/21 (since 3.1)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 7 issues and merged 30 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/milestones/3.2>`__

The following 15 authors contributed 74 commits.

* Benjamin Ragan-Kelley
* Brian Gough
* Damián Avila
* Ian Barfield
* Jason Grout
* Jeff Hussmann
* Jessica B. Hamrick
* Kyle Kelley
* Matthias Bussonnier
* Nicholas Bollweg
* Randy Lai
* Scott Sanderson
* Sylvain Corlay
* Thomas A Caswell
* Thomas Kluyver


Issues closed in 3.1
--------------------

GitHub stats for 2015/02/27 - 2015/04/03 (since 3.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 46 issues and merged 133 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/milestones/3.1>`__.

The following 33 authors contributed 344 commits:

* Abe Guerra
* Adal Chiriliuc
* Benjamin Ragan-Kelley
* Brian Drawert
* Fernando Perez
* Gareth Elston
* Gert-Ludwig Ingold
* Giuseppe Venturini
* Jakob Gager
* Jan Schulz
* Jason Grout
* Jessica B. Hamrick
* Jonathan Frederic
* Justin Tyberg
* Lorena Pantano
* mashenjun
* Mathieu
* Matthias Bussonnier
* Morten Enemark Lund
* Naveen Nathan
* Nicholas Bollweg
* onesandzeroes
* Patrick Snape
* Peter Parente
* RickWinter
* Robert Smith
* Ryan Nelson
* Scott Sanderson
* Sylvain Corlay
* Thomas Kluyver
* tmtabor
* Wieland Hoffmann
* Yuval Langer


Issues closed in 3.0
--------------------

GitHub stats for 2014/04/02 - 2015/02/13 (since 2.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 469 issues and merged 925 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/milestones/3.0>`__.

The following 155 authors contributed 5975 commits.

* A.J. Holyoake
* abalkin
* Adam Hodgen
* Adrian Price-Whelan
* Amin Bandali
* Andreas Amann
* Andrew Dawes
* Andrew Jesaitis
* Andrew Payne
* AnneTheAgile
* Aron Ahmadia
* Ben Duffield
* Benjamin ABEL
* Benjamin Ragan-Kelley
* Benjamin Schultz
* Björn Grüning
* Björn Linse
* Blake Griffith
* Boris Egorov
* Brian E. Granger
* bsvh
* Carlos Cordoba
* Cedric GESTES
* cel
* chebee7i
* Christoph Gohlke
* CJ Carey
* Cyrille Rossant
* Dale Jung
* Damián Avila
* Damon Allen
* Daniel B. Vasquez
* Daniel Rocco
* Daniel Wehner
* Dav Clark
* David Hirschfeld
* David Neto
* dexterdev
* Dimitry Kloper
* dongweiming
* Doug Blank
* drevicko
* Dustin Rodriguez
* Eric Firing
* Eric Galloway
* Erik M. Bray
* Erik Tollerud
* Ezequiel (Zac) Panepucci
* Fernando Perez
* foogunlana
* Francisco de la Peña
* George Titsworth
* Gordon Ball
* gporras
* Grzegorz Rożniecki
* Helen ST
* immerrr
* Ingolf Becker
* Jakob Gager
* James Goppert
* James Porter
* Jan Schulz
* Jason Goad
* Jason Gors
* Jason Grout
* Jason Newton
* jdavidheiser
* Jean-Christophe Jaskula
* Jeff Hemmelgarn
* Jeffrey Bush
* Jeroen Demeyer
* Jessica B. Hamrick
* Jessica Frazelle
* jhemmelg
* Jim Garrison
* Joel Nothman
* Johannes Feist
* John Stowers
* John Zwinck
* jonasc
* Jonathan Frederic
* Juergen Hasch
* Julia Evans
* Justyna Ilczuk
* Jörg Dietrich
* K.-Michael Aye
* Kalibri
* Kester Tong
* Kyle Kelley
* Kyle Rawlins
* Lev Abalkin
* Manuel Riel
* Martin Bergtholdt
* Martin Spacek
* Mateusz Paprocki
* Mathieu
* Matthias Bussonnier
* Maximilian Albert
* mbyt
* MechCoder
* Mohan Raj Rajamanickam
* mvr
* Narahari
* Nathan Goldbaum
* Nathan Heijermans
* Nathaniel J. Smith
* ncornette
* Nicholas Bollweg
* Nick White
* Nikolay Koldunov
* Nile Geisinger
* Olga Botvinnik
* Osada Paranaliyanage
* Pankaj Pandey
* Pascal Bugnion
* patricktokeeffe
* Paul Ivanov
* Peter Odding
* Peter Parente
* Peter Würtz
* Phil Elson
* Phillip Nordwall
* Pierre Gerold
* Pierre Haessig
* Raffaele De Feo
* Ramiro Gómez
* Reggie Pierce
* Remi Rampin
* Renaud Richardet
* Richard Everson
* Scott Sanderson
* Silvia Vinyes
* Simon Guillot
* Spencer Nelson
* Stefan Zimmermann
* Steve Chan
* Steven Anton
* Steven Silvester
* sunny
* Susan Tan
* Sylvain Corlay
* Tarun Gaba
* Thomas Ballinger
* Thomas Kluyver
* Thomas Robitaille
* Thomas Spura
* Tobias Oberstein
* Torsten Bittner
* unknown
* v923z
* vaibhavsagar
* W. Trevor King
* weichm
* Xiuming Chen
* Yaroslav Halchenko
* zah
.. _issues_list_5:

Issues closed in the 5.x development cycle
==========================================

Issues closed in 5.6
--------------------

GitHub stats for 2017/09/15 - 2018/04/02 (tag: 5.5.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 2 issues and merged 28 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A5.6>`__

The following 10 authors contributed 47 commits.

* Benjamin Ragan-Kelley
* Henry Fredrick Schreiner
* Joris Van den Bossche
* Matthias Bussonnier
* Mradul Dubey
* Roshan Rao
* Samuel Lelièvre
* Teddy Rendahl
* Thomas A Caswell
* Thomas Kluyver

Issues closed in 5.4
--------------------

GitHub stats for 2017/02/24 - 2017/05/30 (tag: 5.3.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 8 issues and merged 43 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A5.4+>`__

The following 11 authors contributed 64 commits.

* Benjamin Ragan-Kelley
* Carol Willing
* Kyle Kelley
* Leo Singer
* Luke Pfister
* Lumir Balhar
* Matthias Bussonnier
* meeseeksdev[bot]
* memeplex
* Thomas Kluyver
* Ximin Luo

Issues closed in 5.3
--------------------

GitHub stats for 2017/02/24 - 2017/05/30 (tag: 5.3.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 6 issues and merged 28 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A5.3+>`__

The following 11 authors contributed 53 commits.

* Benjamin Ragan-Kelley
* Carol Willing
* Justin Jent
* Kyle Kelley
* Lumir Balhar
* Matthias Bussonnier
* meeseeksdev[bot]
* Segev Finer
* Steven Maude
* Thomas A Caswell
* Thomas Kluyver


Issues closed in 5.2
--------------------

GitHub stats for 2016/08/13 - 2017/01/29 (tag: 5.1.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 30 issues and merged 74 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A5.2+>`__

The following 40 authors contributed 434 commits.

* Adam Eury
* anantkaushik89
* anatoly techtonik
* Benjamin Ragan-Kelley
* Bibo Hao
* Carl Smith
* Carol Willing
* Chilaka Ramakrishna
* Christopher Welborn
* Denis S. Tereshchenko
* Diego Garcia
* fatData
* Fermi paradox
* Fernando Perez
* fuho
* Hassan Kibirige
* Jamshed Vesuna
* Jens Hedegaard Nielsen
* Jeroen Demeyer
* kaushikanant
* Kenneth Hoste
* Keshav Ramaswamy
* Kyle Kelley
* Matteo
* Matthias Bussonnier
* mbyt
* memeplex
* Moez Bouhlel
* Pablo Galindo
* Paul Ivanov
* pietvo
* Piotr Przetacznik
* Rounak Banik
* sachet-mittal
* Srinivas Reddy Thatiparthy
* Tamir Bahar
* Thomas A Caswell
* Thomas Kluyver
* tillahoffmann
* Yuri Numerov


Issues closed in 5.1
--------------------

GitHub stats for 2016/07/08 - 2016/08/13 (tag: 5.0.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 33 issues and merged 43 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A5.1+>`__

The following 17 authors contributed 129 commits.

* Antony Lee
* Benjamin Ragan-Kelley
* Carol Willing
* Danilo J. S. Bellini
* 小明 (`dongweiming <https://github.com/dongweiming>`__)
* Fernando Perez
* Gavin Cooper
* Gil Forsyth
* Jacob Niehus
* Julian Kuhlmann
* Matthias Bussonnier
* Michael Pacer
* Nik Nyby
* Pavol Juhas
* Luke Deen Taylor
* Thomas Kluyver
* Tamir Bahar


Issues closed in 5.0
--------------------

GitHub stats for 2016/07/05 - 2016/07/07 (tag: 5.0.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 95 issues and merged 191 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A5.0+>`__

The following 27 authors contributed 229 commits.

* Adam Greenhall
* Adrian
* Antony Lee
* Benjamin Ragan-Kelley
* Carlos Cordoba
* Carol Willing
* Chris
* Craig Citro
* Dmitry Zotikov
* Fernando Perez
* Gil Forsyth
* Jason Grout
* Jonathan Frederic
* Jonathan Slenders
* Justin Zymbaluk
* Kelly Liu
* klonuo
* Matthias Bussonnier
* nvdv
* Pavol Juhas
* Pierre Gerold
* sukisuki
* Sylvain Corlay
* Thomas A Caswell
* Thomas Kluyver
* Trevor Bekolay
* Yuri Numerov
========================================
0.9 series
========================================

Release 0.9.1
=============

This release was quickly made to restore compatibility with Python 2.4, which
version 0.9 accidentally broke.  No new features were introduced, other than
some additional testing support for internal use.


Release 0.9
===========

New features
------------

* All furl files and security certificates are now put in a read-only
  directory named ~/.ipython/security.

* A single function :func:`get_ipython_dir`, in :mod:`IPython.genutils` that
  determines the user's IPython directory in a robust manner.

* Laurent's WX application has been given a top-level script called
  ipython-wx, and it has received numerous fixes. We expect this code to be
  architecturally better integrated with Gael's WX 'ipython widget' over the
  next few releases.

* The Editor synchronization work by Vivian De Smedt has been merged in.  This
  code adds a number of new editor hooks to synchronize with editors under
  Windows.

* A new, still experimental but highly functional, WX shell by Gael Varoquaux.
  This work was sponsored by Enthought, and while it's still very new, it is
  based on a more cleanly organized architecture of the various IPython
  components. We will continue to develop this over the next few releases as a
  model for GUI components that use IPython.

* Another GUI frontend, Cocoa based (Cocoa is the OSX native GUI framework),
  authored by Barry Wark.  Currently the WX and the Cocoa ones have slightly
  different internal organizations, but the whole team is working on finding
  what the right abstraction points are for a unified codebase.

* As part of the frontend work, Barry Wark also implemented an experimental
  event notification system that various ipython components can use.  In the
  next release the implications and use patterns of this system regarding the
  various GUI options will be worked out.

* IPython finally has a full test system, that can test docstrings with
  IPython-specific functionality.  There are still a few pieces missing for it
  to be widely accessible to all users (so they can run the test suite at any
  time and report problems), but it now works for the developers.  We are
  working hard on continuing to improve it, as this was probably IPython's
  major Achilles heel (the lack of proper test coverage made it effectively
  impossible to do large-scale refactoring).  The full test suite can now
  be run using the :command:`iptest` command line program.

* The notion of a task has been completely reworked.  An `ITask` interface has
  been created.  This interface defines the methods that tasks need to
  implement.  These methods are now responsible for things like submitting
  tasks and processing results.  There are two basic task types:
  :class:`IPython.kernel.task.StringTask` (this is the old `Task` object, but
  renamed) and the new :class:`IPython.kernel.task.MapTask`, which is based on
  a function.

* A new interface, :class:`IPython.kernel.mapper.IMapper` has been defined to
  standardize the idea of a `map` method.  This interface has a single `map`
  method that has the same syntax as the built-in `map`.  We have also defined
  a `mapper` factory interface that creates objects that implement
  :class:`IPython.kernel.mapper.IMapper` for different controllers.  Both the
  multiengine and task controller now have mapping capabilities.

* The parallel function capabilities have been reworks.  The major changes are
  that i) there is now an `@parallel` magic that creates parallel functions,
  ii) the syntax for multiple variable follows that of `map`, iii) both the
  multiengine and task controller now have a parallel function implementation.

* All of the parallel computing capabilities from `ipython1-dev` have been
  merged into IPython proper.  This resulted in the following new subpackages:
  :mod:`IPython.kernel`, :mod:`IPython.kernel.core`, :mod:`traitlets.config`,
  :mod:`IPython.tools` and :mod:`IPython.testing`.

* As part of merging in the `ipython1-dev` stuff, the `setup.py` script and
  friends have been completely refactored.  Now we are checking for
  dependencies using the approach that matplotlib uses.

* The documentation has been completely reorganized to accept the 
  documentation from `ipython1-dev`.

* We have switched to using Foolscap for all of our network protocols in
  :mod:`IPython.kernel`.  This gives us secure connections that are both
  encrypted and authenticated.

* We have a brand new `COPYING.txt` files that describes the IPython license
  and copyright. The biggest change is that we are putting "The IPython
  Development Team" as the copyright holder. We give more details about
  exactly what this means in this file. All developer should read this and use
  the new banner in all IPython source code files.

* sh profile: ./foo runs foo as system command, no need to do !./foo anymore

* String lists now support ``sort(field, nums = True)`` method (to easily sort
  system command output). Try it with ``a = !ls -l ; a.sort(1, nums=1)``.

* '%cpaste foo' now assigns the pasted block as string list, instead of string

* The ipcluster script now run by default with no security.  This is done
  because the main usage of the script is for starting things on localhost.
  Eventually when ipcluster is able to start things on other hosts, we will put
  security back.

* 'cd --foo' searches directory history for string foo, and jumps to that dir.
  Last part of dir name is checked first. If no matches for that are found,
  look at the whole path.

  
Bug fixes
---------

* The Windows installer has been fixed.  Now all IPython scripts have ``.bat``
  versions created.  Also, the Start Menu shortcuts have been updated.

* The colors escapes in the multiengine client are now turned off on win32 as
  they don't print correctly.

* The :mod:`IPython.kernel.scripts.ipengine` script was exec'ing
  mpi_import_statement incorrectly, which was leading the engine to crash when
  mpi was enabled.

* A few subpackages had missing ``__init__.py`` files.

* The documentation is only created if Sphinx is found.  Previously, the
  ``setup.py`` script would fail if it was missing.

* Greedy ``cd`` completion has been disabled again (it was enabled in 0.8.4) as
  it caused problems on certain platforms.
          

Backwards incompatible changes
------------------------------

* The ``clusterfile`` options of the :command:`ipcluster` command has been
  removed as it was not working and it will be replaced soon by something much
  more robust.

* The :mod:`IPython.kernel` configuration now properly find the user's
  IPython directory.

* In ipapi, the :func:`make_user_ns` function has been replaced with
  :func:`make_user_namespaces`, to support dict subclasses in namespace
  creation.

* :class:`IPython.kernel.client.Task` has been renamed
  :class:`IPython.kernel.client.StringTask` to make way for new task types.

* The keyword argument `style` has been renamed `dist` in `scatter`, `gather`
  and `map`.

* Renamed the values that the rename `dist` keyword argument can have from
  `'basic'` to `'b'`.

* IPython has a larger set of dependencies if you want all of its capabilities.
  See the `setup.py` script for details.

* The constructors for :class:`IPython.kernel.client.MultiEngineClient` and 
  :class:`IPython.kernel.client.TaskClient` no longer take the (ip,port) tuple.
  Instead they take the filename of a file that contains the FURL for that
  client.  If the FURL file is in your IPYTHONDIR, it will be found automatically
  and the constructor can be left empty.

* The asynchronous clients in :mod:`IPython.kernel.asyncclient` are now created 
  using the factory functions :func:`get_multiengine_client` and 
  :func:`get_task_client`.  These return a `Deferred` to the actual client.

* The command line options to `ipcontroller` and `ipengine` have changed to
  reflect the new Foolscap network protocol and the FURL files.  Please see the
  help for these scripts for details.

* The configuration files for the kernel have changed because of the Foolscap
  stuff.  If you were using custom config files before, you should delete them
  and regenerate new ones.

Changes merged in from IPython1
-------------------------------

New features
............

* Much improved ``setup.py`` and ``setupegg.py`` scripts.  Because Twisted and
  zope.interface are now easy installable, we can declare them as dependencies
  in our setupegg.py script.

* IPython is now compatible with Twisted 2.5.0 and 8.x.

* Added a new example of how to use :mod:`ipython1.kernel.asynclient`.

* Initial draft of a process daemon in :mod:`ipython1.daemon`.  This has not
  been merged into IPython and is still in `ipython1-dev`.
  
* The ``TaskController`` now has methods for getting the queue status.

* The ``TaskResult`` objects not have information about how long the task
  took to run.
  
* We are attaching additional attributes to exceptions ``(_ipython_*)`` that 
  we use to carry additional info around.
  
* New top-level module :mod:`asyncclient` that has asynchronous versions (that
  return deferreds) of the client classes.  This is designed to users who want
  to run their own Twisted reactor.
  
* All the clients in :mod:`client` are now based on Twisted.  This is done by 
  running the Twisted reactor in a separate thread and using the
  :func:`blockingCallFromThread` function that is in recent versions of Twisted.

* Functions can now be pushed/pulled to/from engines using
  :meth:`MultiEngineClient.push_function` and
  :meth:`MultiEngineClient.pull_function`.

* Gather/scatter are now implemented in the client to reduce the work load
  of the controller and improve performance.

* Complete rewrite of the IPython docuementation.  All of the documentation
  from the IPython website has been moved into docs/source as restructured
  text documents.  PDF and HTML documentation are being generated using 
  Sphinx.

* New developer oriented documentation: development guidelines and roadmap. 

* Traditional ``ChangeLog`` has been changed to a more useful ``changes.txt``
  file that is organized by release and is meant to provide something more
  relevant for users.

Bug fixes
.........

* Created a proper ``MANIFEST.in`` file to create source distributions.

* Fixed a bug in the ``MultiEngine`` interface.  Previously, multi-engine 
  actions were being collected with a :class:`DeferredList` with 
  ``fireononeerrback=1``.  This meant that methods were returning 
  before all engines had given their results.  This was causing extremely odd 
  bugs in certain cases. To fix this problem, we have 1) set 
  ``fireononeerrback=0`` to make sure all results (or exceptions) are in 
  before returning and 2) introduced a :exc:`CompositeError` exception 
  that wraps all of the engine exceptions.  This is a huge change as it means 
  that users will have to catch :exc:`CompositeError` rather than the actual
  exception.

Backwards incompatible changes
..............................

* All names have been renamed to conform to the lowercase_with_underscore
  convention.  This will require users to change references to all names like
  ``queueStatus`` to ``queue_status``.

* Previously, methods like :meth:`MultiEngineClient.push` and
  :meth:`MultiEngineClient.push` used ``*args`` and ``**kwargs``.  This was
  becoming a problem as we weren't able to introduce new keyword arguments into
  the API.  Now these methods simple take a dict or sequence.  This has also
  allowed us to get rid of the ``*All`` methods like :meth:`pushAll` and
  :meth:`pullAll`.  These things are now handled with the ``targets`` keyword
  argument that defaults to ``'all'``.

* The :attr:`MultiEngineClient.magicTargets` has been renamed to
  :attr:`MultiEngineClient.targets`. 

* All methods in the MultiEngine interface now accept the optional keyword
  argument ``block``.

* Renamed :class:`RemoteController` to :class:`MultiEngineClient` and 
  :class:`TaskController` to :class:`TaskClient`.

* Renamed the top-level module from :mod:`api` to :mod:`client`.

* Most methods in the multiengine interface now raise a :exc:`CompositeError`
  exception that wraps the user's exceptions, rather than just raising the raw
  user's exception.

* Changed the ``setupNS`` and ``resultNames`` in the ``Task`` class to ``push`` 
  and ``pull``.

============
 4.x Series
============

IPython 4.2
===========

IPython 4.2 (April, 2016) includes various bugfixes and improvements over 4.1.

- Fix ``ipython -i`` on errors, which was broken in 4.1.
- The delay meant to highlight deprecated commands that have moved to jupyter has been removed.
- Improve compatibility with future versions of traitlets and matplotlib.
- Use stdlib :func:`python:shutil.get_terminal_size` to measure terminal width when displaying tracebacks
  (provided by ``backports.shutil_get_terminal_size`` on Python 2).

You can see the rest `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A4.2>`__.


IPython 4.1
===========

IPython 4.1.2 (March, 2016) fixes installation issues with some versions of setuptools.

Released February, 2016. IPython 4.1 contains mostly bug fixes,
though there are a few improvements.


- IPython debugger (IPdb) now supports the number of context lines for the
  ``where`` (and ``w``) commands. The `context` keyword is also available in
  various APIs. See PR :ghpull:`9097`
- YouTube video will now show thumbnail when exported to a media that do not
  support video. (:ghpull:`9086`)
- Add warning when running `ipython <subcommand>` when subcommand is
  deprecated. `jupyter` should now be used.
- Code in `%pinfo` (also known as `??`) are now highlighter (:ghpull:`8947`)
- `%aimport` now support module completion. (:ghpull:`8884`)
- `ipdb` output is now colored ! (:ghpull:`8842`)
- Add ability to transpose columns for completion: (:ghpull:`8748`)

Many many docs improvements and bug fixes, you can see the
`list of changes <https://github.com/ipython/ipython/compare/4.0.0...4.1.0>`_

IPython 4.0
===========

Released August, 2015

IPython 4.0 is the first major release after the Big Split.
IPython no longer contains the notebook, qtconsole, etc. which have moved to
`jupyter <https://jupyter.readthedocs.io>`_.
IPython subprojects, such as `IPython.parallel <https://ipyparallel.readthedocs.io>`_ and `widgets <https://ipywidgets.readthedocs.io>`_ have moved to their own repos as well.

The following subpackages are deprecated:

- IPython.kernel (now jupyter_client and ipykernel)
- IPython.consoleapp (now jupyter_client.consoleapp)
- IPython.nbformat (now nbformat)
- IPython.nbconvert (now nbconvert)
- IPython.html (now notebook)
- IPython.parallel (now ipyparallel)
- IPython.utils.traitlets (now traitlets)
- IPython.config (now traitlets.config)
- IPython.qt (now qtconsole)
- IPython.terminal.console (now jupyter_console)

and a few other utilities.

Shims for the deprecated subpackages have been added,
so existing code should continue to work with a warning about the new home.

There are few changes to the code beyond the reorganization and some bugfixes.

IPython highlights:

- Public APIs for discovering IPython paths is moved from :mod:`IPython.utils.path` to :mod:`IPython.paths`.
  The old function locations continue to work with deprecation warnings.
- Code raising ``DeprecationWarning``
  entered by the user in an interactive session will now display the warning by
  default. See :ghpull:`8480` an :ghissue:`8478`.
- The `--deep-reload` flag and the corresponding options to inject `dreload` or
  `reload` into the interactive namespace have been deprecated, and will be
  removed in future versions. You should now explicitly import `reload` from
  `IPython.lib.deepreload` to use it.

============
 1.0 Series
============

Release 1.0.0: An Afternoon Hack
================================


IPython 1.0 requires Python ≥ 2.6.5 or ≥ 3.2.1.
It does not support Python 3.0, 3.1, or 2.5.

This is a big release.  The principal milestone is the addition of :mod:`IPython.nbconvert`,
but there has been a great deal of work improving all parts of IPython as well.

The previous version (0.13) was released on June 30, 2012,
and in this development cycle we had:

- ~12 months of work.
- ~700 pull requests merged.
- ~600 issues closed (non-pull requests).
- contributions from ~150 authors.
- ~4000 commits.

The amount of work included in this release is so large that we can only cover
here the main highlights; please see our :ref:`detailed release statistics
<issues_list_100>` for links to every issue and pull request closed on GitHub
as well as a full list of individual contributors.
It includes

Reorganization
--------------

There have been two major reorganizations in IPython 1.0:

- Added :mod:`IPython.kernel` for all kernel-related code.
  This means that :mod:`IPython.zmq` has been removed,
  and much of it is now in :mod:`IPython.kernel.zmq`,
  some of it being in the top-level :mod:`IPython.kernel`.
- We have removed the `frontend` subpackage,
  as it caused unnecessary depth.  So what was :mod:`IPython.frontend.qt`
  is now :mod:`IPython.qt`, and so on.  The one difference is that
  the notebook has been further flattened, so that
  :mod:`IPython.frontend.html.notebook` is now just `IPython.html`.
  There is a shim module, so :mod:`IPython.frontend` is still
  importable in 1.0, but there will be a warning.
- The IPython sphinx directives are now installed in :mod:`IPython.sphinx`,
  so they can be imported by other projects.


Public APIs
-----------

For the first time since 0.10 (sorry, everyone),
there is an official public API for starting IPython:

.. sourcecode:: python

    from IPython import start_ipython
    start_ipython()

This is what packages should use that start their own IPython session,
but don't actually want embedded IPython (most cases).
:func:`IPython.embed()` is used for embedding IPython into the calling namespace,
similar to calling :func:`Pdb.set_trace`, whereas :func:`start_ipython`
will start a plain IPython session, loading config and startup files as normal.

We also have added:

.. sourcecode:: python

    from IPython import get_ipython


Which is a *library* function for getting the current IPython instance,
and will return ``None`` if no IPython instance is running.
This is the official way to check whether your code is called from inside an IPython session.
If you want to check for IPython without unnecessarily importing IPython,
use this function:

.. sourcecode:: python

    def get_ipython():
        """return IPython instance if there is one, None otherwise"""
        import sys
        if "IPython" in sys.modules:
            import IPython
            return IPython.get_ipython()

Core
----

- The input transformation framework has been reworked. This fixes some corner
  cases, and adds more flexibility for projects which use IPython, like SymPy &
  SAGE. For more details, see :doc:`/config/inputtransforms`.
- Exception types can now be displayed with a custom traceback, by defining a
  ``_render_traceback_()`` method which returns a list of strings, each
  containing one line of the traceback.
- A new command, ``ipython history trim`` can be used to delete everything but
  the last 1000 entries in the history database.
- ``__file__`` is defined in both config files at load time,
  and ``.ipy`` files executed with ``%run``.
- ``%logstart`` and ``%logappend`` are no longer broken.
- Add glob expansion for ``%run``, e.g. ``%run -g script.py *.txt``.
- Expand variables (``$foo``) in Cell Magic argument line.
- By default, :command:`iptest` will exclude various slow tests.
  All tests can be run with :command:`iptest --all`.
- SQLite history can be disabled in the various cases that it does not behave well.
- ``%edit`` works on interactively defined variables.
- editor hooks have been restored from quarantine, enabling TextMate as editor,
  etc.
- The env variable PYTHONSTARTUP is respected by IPython.
- The ``%matplotlib`` magic was added, which is like the old ``%pylab`` magic,
  but it does not import anything to the interactive namespace.
  It is recommended that users switch to ``%matplotlib`` and explicit imports.
- The ``--matplotlib`` command line flag was also added. It invokes the new
  ``%matplotlib`` magic and can be used in the same way as the old ``--pylab``
  flag. You can either use it by itself as a flag (``--matplotlib``), or you
  can also pass a backend explicitly (``--matplotlib qt`` or
  ``--matplotlib=wx``, etc).


Backwards incompatible changes
******************************

- Calling :meth:`InteractiveShell.prefilter` will no longer perform static
  transformations - the processing of escaped commands such as ``%magic`` and
  ``!system``, and stripping input prompts from code blocks. This functionality
  was duplicated in :mod:`IPython.core.inputsplitter`, and the latter version
  was already what IPython relied on. A new API to transform input will be ready
  before release.
- Functions from :mod:`IPython.lib.inputhook` to control integration with GUI
  event loops are no longer exposed in the top level of :mod:`IPython.lib`.
  Code calling these should make sure to import them from
  :mod:`IPython.lib.inputhook`.
- For all kernel managers, the ``sub_channel`` attribute has been renamed to
  ``iopub_channel``.
- Users on Python versions before 2.6.6, 2.7.1 or 3.2 will now need to call
  :func:`IPython.utils.doctestreload.doctest_reload` to make doctests run 
  correctly inside IPython. Python releases since those versions are unaffected.
  For details, see :ghpull:`3068` and `Python issue 8048 <http://bugs.python.org/issue8048>`_.
- The ``InteractiveShell.cache_main_mod()`` method has been removed, and
  :meth:`~IPython.core.interactiveshell.InteractiveShell.new_main_mod` has a
  different signature, expecting a filename where earlier versions expected
  a namespace. See :ghpull:`3555` for details.
- The short-lived plugin system has been removed. Extensions are the way to go.


.. _nbconvert1:

NbConvert
---------

The major milestone for IPython 1.0 is the addition of :mod:`IPython.nbconvert` - tools for converting
IPython notebooks to various other formats.

.. warning::

    nbconvert is α-level preview code in 1.0

To use nbconvert to convert various file formats::

    ipython nbconvert --to html *.ipynb

See ``ipython nbconvert --help`` for more information.
nbconvert depends on `pandoc`_ for many of the translations to and from various formats.

.. _pandoc: http://johnmacfarlane.net/pandoc/

Notebook
--------

Major changes to the IPython Notebook in 1.0:

- The notebook is now autosaved, by default at an interval of two minutes.
  When you press 'save' or Ctrl-S, a *checkpoint* is made, in a hidden folder.
  This checkpoint can be restored, so that the autosave model is strictly safer
  than traditional save. If you change nothing about your save habits,
  you will always have a checkpoint that you have written,
  and an autosaved file that is kept up to date.
- The notebook supports :func:`raw_input` / :func:`input`, and thus also ``%debug``,
  and many other Python calls that expect user input.
- You can load custom javascript and CSS in the notebook by editing the files
  :file:`$(ipython locate profile)/static/custom/custom.{js,css}`.
- Add ``%%html``, ``%%svg``, ``%%javascript``, and ``%%latex`` cell magics
  for writing raw output in notebook cells.
- add a redirect handler and anchors on heading cells, so you can link
  across notebooks, directly to heading cells in other notebooks.
- Images support width and height metadata,
  and thereby 2x scaling (retina support).
- ``_repr_foo_`` methods can return a tuple of (data, metadata),
  where metadata is a dict containing metadata about the displayed object.
  This is used to set size, etc. for retina graphics. To enable retina matplotlib figures,
  simply set ``InlineBackend.figure_format = 'retina'`` for 2x PNG figures,
  in your :ref:`IPython config file <config_overview>` or via the ``%config`` magic.
- Add display.FileLink and FileLinks for quickly displaying HTML links to local files.
- Cells have metadata, which can be edited via cell toolbars.
  This metadata can be used by external code (e.g. reveal.js or exporters),
  when examining the notebook.
- Fix an issue parsing LaTeX in markdown cells, which required users to type ``\\\``,
  instead of ``\\``.
- Notebook templates are rendered with Jinja instead of Tornado.
- ``%%file`` has been renamed ``%%writefile`` (``%%file`` is deprecated).
- ANSI (and VT100) color parsing has been improved in both performance and
  supported values.
- The static files path can be found as ``IPython.html.DEFAULT_STATIC_FILES_PATH``,
  which may be changed by package managers.
- IPython's CSS is installed in :file:`static/css/style.min.css`
  (all style, including bootstrap), and :file:`static/css/ipython.min.css`,
  which only has IPython's own CSS. The latter file should be useful for embedding
  IPython notebooks in other pages, blogs, etc.
- The Print View has been removed. Users are encouraged to test :ref:`ipython
  nbconvert <nbconvert1>` to generate a static view.

Javascript Components
*********************

The javascript components used in the notebook have been updated significantly.

- updates to jQuery (2.0) and jQueryUI (1.10)
- Update CodeMirror to 3.14
- Twitter Bootstrap (2.3) for layout
- Font-Awesome (3.1) for icons
- highlight.js (7.3) for syntax highlighting
- marked (0.2.8) for markdown rendering
- require.js (2.1) for loading javascript

Some relevant changes that are results of this:

- Markdown cells now support GitHub-flavored Markdown (GFM),
  which includes `````python`` code blocks and tables.
- Notebook UI behaves better on more screen sizes.
- Various code cell input issues have been fixed.


Kernel
------

The kernel code has been substantially reorganized.

New features in the kernel:

- Kernels support ZeroMQ IPC transport, not just TCP
- The message protocol has added a top-level metadata field,
  used for information about messages.
- Add a `data_pub` message that functions much like `display_pub`,
  but publishes raw (usually pickled) data, rather than representations.
- Ensure that ``sys.stdout.encoding`` is defined in Kernels.
- Stdout from forked subprocesses should be forwarded to frontends (instead of crashing).

IPEP 13
*******

The KernelManager has been split into a :class:`~.KernelManager` and a :class:`~.KernelClient`.
The Manager owns a kernel and starts / signals / restarts it. There is always zero or one
KernelManager per Kernel.  Clients communicate with Kernels via zmq channels,
and there can be zero-to-many Clients connected to a Kernel at any given time.

The KernelManager now automatically restarts the kernel when it dies,
rather than requiring user input at the notebook or QtConsole UI
(which may or may not exist at restart time).

In-process kernels
******************

The Python-language frontends, particularly the Qt console, may now communicate
with in-process kernels, in addition to the traditional out-of-process
kernels. An in-process kernel permits direct access to the kernel namespace,
which is necessary in some applications. It should be understood, however, that
the in-process kernel is not robust to bad user input and will block the main
(GUI) thread while executing. Developers must decide on a case-by-case basis
whether this tradeoff is appropriate for their application.



Parallel
--------

IPython.parallel has had some refactoring as well.  
There are many improvements and fixes, but these are the major changes:

- Connections have been simplified. All ports and the serialization in use
  are written to the connection file, rather than the initial two-stage system.
- Serialization has been rewritten, fixing many bugs and dramatically improving
  performance serializing large containers.
- Load-balancing scheduler performance with large numbers of tasks has been dramatically improved.
- There should be fewer (hopefully zero) false-positives for engine failures.
- Increased compatibility with various use cases that produced serialization / argument errors
  with map, etc.
- The controller can attempt to resume operation if it has crashed,
  by passing ``ipcontroller --restore``.
- Engines can monitor the Hub heartbeat, and shutdown if the Hub disappears for too long.
- add HTCondor support in launchers


QtConsole
---------

Various fixes, including improved performance with lots of text output,
and better drag and drop support.
The initial window size of the qtconsole is now configurable via ``IPythonWidget.width``
and ``IPythonWidget.height``.

.. _issues_list_013:

Issues closed in the 0.13 development cycle
===========================================

Issues closed in 0.13
---------------------

GitHub stats since IPython 0.12 (2011/12/19 - 2012/06/30)

These lists are automatically generated, and may be incomplete or contain
duplicates.

The following 62 authors contributed 1760 commits.

* Aaron Culich
* Aaron Meurer
* Alex Kramer
* Andrew Giessel
* Andrew Straw
* André Matos
* Aron Ahmadia
* Ben Edwards
* Benjamin Ragan-Kelley
* Bradley M. Froehle
* Brandon Parsons
* Brian E. Granger
* Carlos Cordoba
* David Hirschfeld
* David Zderic
* Ernie French
* Fernando Perez
* Ian Murray
* Jason Grout
* Jens H Nielsen
* Jez Ng
* Jonathan March
* Jonathan Taylor
* Julian Taylor
* Jörgen Stenarson
* Kent Inverarity
* Marc Abramowitz
* Mark Wiebe
* Matthew Brett
* Matthias BUSSONNIER
* Michael Droettboom
* Mike Hansen
* Nathan Rice
* Pankaj Pandey
* Paul
* Paul Ivanov
* Piotr Zolnierczuk
* Piti Ongmongkolkul
* Puneeth Chaganti
* Robert Kern
* Ross Jones
* Roy Hyunjin Han
* Scott Tsai
* Skipper Seabold
* Stefan van der Walt
* Steven Johnson
* Takafumi Arakaki
* Ted Wright
* Thomas Hisch
* Thomas Kluyver
* Thomas Spura
* Thomi Richards
* Tim Couper
* Timo Paulssen
* Toby Gilham
* Tony S Yu
* W. Trevor King
* Walter Doerwald
* anatoly techtonik
* fawce
* mcelrath
* wilsaj


We closed a total of 1115 issues, 373 pull requests and 742 regular issues;
this is the full list (generated with the script 
:file:`tools/github_stats.py`):

Pull Requests (373):

* :ghpull:`1943`: add screenshot and link into releasenotes
* :ghpull:`1954`: update some example notebooks
* :ghpull:`2048`: move _encode_binary to jsonutil.encode_images
* :ghpull:`2050`: only add quotes around xunit-file on Windows
* :ghpull:`2047`: disable auto-scroll on mozilla
* :ghpull:`2015`: Fixes for %paste with special transformations
* :ghpull:`2046`: Iptest unicode
* :ghpull:`1939`: Namespaces
* :ghpull:`2042`: increase auto-scroll threshold to 100 lines
* :ghpull:`2043`: move RemoteError import to top-level
* :ghpull:`2036`: %alias_magic
* :ghpull:`1968`: Proposal of icons for .ipynb files
* :ghpull:`2037`: remove `ipython-qtconsole` gui-script
* :ghpull:`2038`: add extra clear warning to shell doc
* :ghpull:`2029`: Ship unminified js
* :ghpull:`2007`: Add custom_control and custom_page_control variables to override the Qt widgets used by qtconsole
* :ghpull:`2034`: fix&test push/pull recarrays
* :ghpull:`2028`: Reduce unhelpful information shown by pinfo
* :ghpull:`2030`: check wxPython version in inputhook
* :ghpull:`2024`: Make interactive_usage a bit more rst friendly
* :ghpull:`2031`: disable ^C^C confirmation on Windows
* :ghpull:`2027`: match stdin encoding in frontend readline test
* :ghpull:`2025`: Fix parallel test on WinXP - wait for resource cleanup.
* :ghpull:`2016`: BUG: test runner fails in Windows if filenames contain spaces.
* :ghpull:`2020`: Fix home path expansion test in Windows.
* :ghpull:`2021`: Fix Windows pathname issue in 'odd encoding' test.
* :ghpull:`2022`: don't check writability in test for get_home_dir when HOME is undefined
* :ghpull:`1996`: frontend test tweaks
* :ghpull:`2014`: relax profile regex in notebook
* :ghpull:`2012`: Mono cursor offset
* :ghpull:`2004`: Clarify generic message spec vs. Python message API in docs
* :ghpull:`2010`: notebook: Print a warning (but do not abort) if no webbrowser can be found.
* :ghpull:`2002`: Refactor %magic into a lsmagic_docs API function.
* :ghpull:`1999`: `%magic` help: display line and cell magics in alphabetical order.
* :ghpull:`1981`: Clean BG processes created by %%script on kernel exit
* :ghpull:`1994`: Fix RST misformatting.
* :ghpull:`1951`: minor notebook startup/notebook-dir adjustments
* :ghpull:`1974`: Allow path completion on notebook.
* :ghpull:`1964`: allow multiple instances of a Magic
* :ghpull:`1991`: fix _ofind attr in %page
* :ghpull:`1988`: check for active frontend in update_restart_checkbox
* :ghpull:`1979`: Add support for tox (https://tox.readthedocs.io/) and Travis CI (http://travis-ci.org/)
* :ghpull:`1970`: dblclick to restore size of images
* :ghpull:`1978`: Notebook names truncating at the first period
* :ghpull:`1825`: second attempt at scrolled long output
* :ghpull:`1934`: Cell/Worksheet metadata
* :ghpull:`1746`: Confirm restart (configuration option, and checkbox UI)
* :ghpull:`1944`: [qtconsole] take %,%% prefix into account for completion
* :ghpull:`1973`: fix another FreeBSD $HOME symlink issue
* :ghpull:`1967`: Fix psums example description in docs
* :ghpull:`1965`: fix for #1678, undo no longer clears cells
* :ghpull:`1952`: avoid duplicate "Websockets closed" dialog on ws close
* :ghpull:`1962`: Support unicode prompts
* :ghpull:`1955`: update to latest version of vim-ipython
* :ghpull:`1945`: Add --proc option to %%script
* :ghpull:`1956`: move import RemoteError after get_exc_info
* :ghpull:`1950`: Fix for copy action (Ctrl+C) when there is no pager defined in qtconsole
* :ghpull:`1948`: Fix help string for InteractiveShell.ast_node_interactivity
* :ghpull:`1942`: swallow stderr of which in utils.process.find_cmd
* :ghpull:`1940`: fix completer css on some Chrome versions
* :ghpull:`1938`: remove remaining references to deprecated XREP/XREQ names
* :ghpull:`1925`: Fix styling of superscripts and subscripts. Closes #1924.
* :ghpull:`1936`: increase duration of save messages
* :ghpull:`1937`: add %save -f
* :ghpull:`1935`: add version checking to pyreadline import test
* :ghpull:`1849`: Octave magics
* :ghpull:`1759`: github, merge PR(s) just by number(s) 
* :ghpull:`1931`: Win py3fixes
* :ghpull:`1933`: oinspect.find_file: Additional safety if file cannot be found.
* :ghpull:`1932`: Fix adding functions to CommandChainDispatcher with equal priority on Py 3
* :ghpull:`1928`: Select NoDB by default
* :ghpull:`1923`: Add IPython syntax support to the %timeit magic, in line and cell mode
* :ghpull:`1926`: Make completer recognize escaped quotes in strings.
* :ghpull:`1893`: Update Parallel Magics and Exception Display
* :ghpull:`1921`: magic_arguments: dedent but otherwise preserve indentation.
* :ghpull:`1919`: Use oinspect in CodeMagics._find_edit_target
* :ghpull:`1918`: don't warn in iptest if deathrow/quarantine are missing
* :ghpull:`1917`: Fix for %pdef on Python 3
* :ghpull:`1913`: Fix for #1428
* :ghpull:`1911`: temporarily skip autoreload tests
* :ghpull:`1909`: Fix for #1908, use os.path.normcase for safe filename comparisons
* :ghpull:`1907`: py3compat fixes for %%script and tests
* :ghpull:`1906`: ofind finds non-unique cell magics
* :ghpull:`1845`: Fixes to inspection machinery for magics
* :ghpull:`1902`: Workaround fix for gh-1632; minimal revert of gh-1424
* :ghpull:`1900`: Cython libs
* :ghpull:`1899`: add ScriptMagics to class list for generated config
* :ghpull:`1898`: minimize manpages
* :ghpull:`1897`: use glob for bad exclusion warning
* :ghpull:`1855`: %%script and %%file magics
* :ghpull:`1870`: add %%capture for capturing stdout/err
* :ghpull:`1861`: Use dvipng to format sympy.Matrix
* :ghpull:`1867`: Fix 1px margin bouncing of selected menu item.
* :ghpull:`1889`: Reconnect when the websocket connection closes unexpectedly
* :ghpull:`1886`: Fix a bug in renaming notebook
* :ghpull:`1895`: Fix error in test suite with ip.system()
* :ghpull:`1762`: add `locate` entry points
* :ghpull:`1883`: Fix vertical offset due to bold/italics, and bad browser fonts.
* :ghpull:`1875`: re-write columnize, with intermediate step.
* :ghpull:`1851`: new completer for qtconsole.
* :ghpull:`1892`: Remove suspicious quotes in interactiveshell.py
* :ghpull:`1864`: Rmagic exceptions
* :ghpull:`1829`: [notebook] don't care about leading prct in completion
* :ghpull:`1832`: Make svg, jpeg and png images resizable in notebook.
* :ghpull:`1674`: HTML Notebook carriage-return handling, take 2
* :ghpull:`1882`: Remove importlib dependency which not available in Python 2.6.
* :ghpull:`1879`: Correct stack depth for variable expansion in !system commands
* :ghpull:`1841`: [notebook] deduplicate completion results
* :ghpull:`1850`: Remove args/kwargs handling in TryNext, fix %paste error messages.
* :ghpull:`1663`: Keep line-endings in ipynb
* :ghpull:`1815`: Make : invalid in filenames in the Notebook JS code.
* :ghpull:`1819`: doc: cleanup the parallel psums example a little
* :ghpull:`1839`: External cleanup
* :ghpull:`1782`: fix Magic menu in qtconsole, split in groups
* :ghpull:`1862`: Minor bind_kernel improvements
* :ghpull:`1857`: Prevent jumping of window to input when output is clicked.
* :ghpull:`1856`: Fix 1px jumping of cells and menus in Notebook.
* :ghpull:`1852`: fix chained resubmissions
* :ghpull:`1780`: Rmagic extension
* :ghpull:`1847`: add InlineBackend to ConsoleApp class list
* :ghpull:`1836`: preserve header for resubmitted tasks
* :ghpull:`1828`: change default extension to .ipy for %save -r
* :ghpull:`1800`: Reintroduce recall
* :ghpull:`1830`: lsmagic lists magics in alphabetical order
* :ghpull:`1773`: Update SymPy profile: SymPy's latex() can now print set and frozenset
* :ghpull:`1761`: Edited documentation to use IPYTHONDIR in place of ~/.ipython
* :ghpull:`1822`: aesthetics pass on AsyncResult.display_outputs
* :ghpull:`1821`: ENTER submits the rename notebook dialog.
* :ghpull:`1820`: NotebookApp: Make the number of ports to retry user configurable.
* :ghpull:`1816`: Always use filename as the notebook name.
* :ghpull:`1813`: Add assert_in method to nose for Python 2.6
* :ghpull:`1711`: New Tooltip, New Completer and JS Refactor
* :ghpull:`1798`: a few simple fixes for docs/parallel
* :ghpull:`1812`: Ensure AsyncResult.display_outputs doesn't display empty streams
* :ghpull:`1811`: warn on nonexistent exclusions in iptest
* :ghpull:`1810`: fix for #1809, failing tests in IPython.zmq
* :ghpull:`1808`: Reposition alternate upload for firefox [need cross browser/OS/language test]
* :ghpull:`1742`: Check for custom_exceptions only once
* :ghpull:`1807`: add missing cython exclusion in iptest
* :ghpull:`1805`: Fixed a vcvarsall.bat error on win32/Py2.7 when trying to compile with m...
* :ghpull:`1739`: Dashboard improvement (necessary merge of #1658 and #1676 + fix #1492)
* :ghpull:`1770`: Cython related magic functions
* :ghpull:`1707`: Accept --gui=<...> switch in IPython qtconsole.
* :ghpull:`1797`: Fix comment which breaks Emacs syntax highlighting.
* :ghpull:`1795`: fix %gui magic
* :ghpull:`1793`: Raise repr limit for strings to 80 characters (from 30).
* :ghpull:`1794`: don't use XDG path on OS X
* :ghpull:`1792`: Unicode-aware logger
* :ghpull:`1791`: update zmqshell magics
* :ghpull:`1787`: DOC: Remove regression from qt-console docs.
* :ghpull:`1758`: test_pr, fallback on http if git protocol fail, and SSL errors...
* :ghpull:`1748`: Fix some tests for Python 3.3
* :ghpull:`1755`: test for pygments before running qt tests
* :ghpull:`1771`: Make default value of interactivity passed to run_ast_nodes configurable
* :ghpull:`1784`: restore loadpy to load
* :ghpull:`1768`: Update parallel magics
* :ghpull:`1779`: Tidy up error raising in magic decorators.
* :ghpull:`1769`: Allow cell mode timeit without setup code.
* :ghpull:`1716`: Fix for fake filenames in verbose traceback
* :ghpull:`1763`: [qtconsole] fix append_plain_html -> append_html
* :ghpull:`1732`: Refactoring of the magics system and implementation of cell magics
* :ghpull:`1630`: Merge divergent Kernel implementations
* :ghpull:`1705`: [notebook] Make pager resizable, and remember size...
* :ghpull:`1606`: Share code for %pycat and %loadpy, make %pycat aware of URLs
* :ghpull:`1757`: Open IPython notebook hyperlinks in a new window using target=_blank
* :ghpull:`1754`: Fix typo enconters->encounters
* :ghpull:`1753`: Clear window title when kernel is restarted
* :ghpull:`1449`: Fix for bug #735 : Images missing from XML/SVG export
* :ghpull:`1743`: Tooltip completer js refactor
* :ghpull:`1681`: add qt config option to clear_on_kernel_restart
* :ghpull:`1733`: Tooltip completer js refactor
* :ghpull:`1727`: terminate kernel after embed_kernel tests
* :ghpull:`1737`: add HistoryManager to ipapp class list
* :ghpull:`1686`: ENH: Open a notebook from the command line
* :ghpull:`1709`: fixes #1708, failing test in arg_split on windows
* :ghpull:`1718`: Use CRegExp trait for regular expressions.
* :ghpull:`1729`: Catch failure in repr() for %whos
* :ghpull:`1726`: use eval for command-line args instead of exec
* :ghpull:`1724`: fix scatter/gather with targets='all'
* :ghpull:`1725`: add --no-ff to git pull in test_pr
* :ghpull:`1721`: Tooltip completer js refactor
* :ghpull:`1657`: Add `wait` optional argument to `hooks.editor`
* :ghpull:`1717`: Define generic sys.ps{1,2,3}, for use by scripts.
* :ghpull:`1691`: Finish PR #1446
* :ghpull:`1710`: update MathJax CDN url for https
* :ghpull:`1713`: Make autocall regexp's configurable.
* :ghpull:`1703`: Allow TryNext to have an error message without it affecting the command chain
* :ghpull:`1714`: minor adjustments to test_pr
* :ghpull:`1704`: ensure all needed qt parts can be imported before settling for one
* :ghpull:`1706`: Mark test_push_numpy_nocopy as a known failure for Python 3
* :ghpull:`1698`: fix tooltip on token with number
* :ghpull:`1245`: pythonw py3k fixes for issue #1226
* :ghpull:`1685`: Add script to test pull request
* :ghpull:`1693`: deprecate IPYTHON_DIR in favor of IPYTHONDIR
* :ghpull:`1695`: Avoid deprecated warnings from ipython-qtconsole.desktop.
* :ghpull:`1694`: Add quote to notebook to allow it to load
* :ghpull:`1689`: Fix sys.path missing '' as first entry in `ipython kernel`.
* :ghpull:`1687`: import Binary from bson instead of pymongo
* :ghpull:`1616`: Make IPython.core.display.Image less notebook-centric
* :ghpull:`1684`: CLN: Remove redundant function definition.
* :ghpull:`1670`: Point %pastebin to gist
* :ghpull:`1669`: handle pyout messages in test_message_spec
* :ghpull:`1295`: add binary-tree engine interconnect example
* :ghpull:`1642`: Cherry-picked commits from 0.12.1 release
* :ghpull:`1659`: Handle carriage return characters ("\r") in HTML notebook output.
* :ghpull:`1656`: ensure kernels are cleaned up in embed_kernel tests
* :ghpull:`1664`: InteractiveShell.run_code: Update docstring.
* :ghpull:`1662`: Delay flushing softspace until after cell finishes
* :ghpull:`1643`: handle jpg/jpeg in the qtconsole
* :ghpull:`1652`: add patch_pyzmq() for backporting a few changes from newer pyzmq
* :ghpull:`1650`: DOC: moving files with SSH launchers
* :ghpull:`1357`: add IPython.embed_kernel() 
* :ghpull:`1640`: Finish up embed_kernel
* :ghpull:`1651`: Remove bundled Itpl module
* :ghpull:`1634`: incremental improvements to SSH launchers
* :ghpull:`1649`: move examples/test_embed into examples/tests/embed
* :ghpull:`1633`: Fix installing extension from local file on Windows
* :ghpull:`1645`: Exclude UserDict when deep reloading NumPy.
* :ghpull:`1637`: Removed a ':' which shouldn't have been there
* :ghpull:`1631`: TST: QApplication doesn't quit early enough with PySide.
* :ghpull:`1629`: evaluate a few dangling validate_message generators
* :ghpull:`1621`: clear In[] prompt numbers on "Clear All Output"
* :ghpull:`1627`: Test the Message Spec
* :ghpull:`1624`: Fixes for byte-compilation on Python 3
* :ghpull:`1615`: Add show() method to figure objects.
* :ghpull:`1625`: Fix deepreload on Python 3
* :ghpull:`1620`: pyin message now have execution_count
* :ghpull:`1457`: Update deepreload to use a rewritten knee.py. Fixes dreload(numpy).
* :ghpull:`1613`: allow map / parallel function for single-engine views
* :ghpull:`1609`: exit notebook cleanly on SIGINT, SIGTERM
* :ghpull:`1607`: cleanup sqlitedb temporary db file after tests
* :ghpull:`1608`: don't rely on timedelta.total_seconds in AsyncResult
* :ghpull:`1599`: Fix for %run -d on Python 3
* :ghpull:`1602`: Fix %env magic on Python 3.
* :ghpull:`1603`: Remove python3 profile
* :ghpull:`1604`: Exclude IPython.quarantine from installation
* :ghpull:`1600`: Specify encoding for io.open in notebook_reformat tests
* :ghpull:`1605`: Small fixes for Animation and Progress notebook
* :ghpull:`1529`: __all__ feature, improvement to dir2, and tests for both
* :ghpull:`1548`: add sugar methods/properties to AsyncResult
* :ghpull:`1535`: Fix pretty printing dispatch
* :ghpull:`1399`: Use LaTeX to print various built-in types with the SymPy printing extension
* :ghpull:`1597`: re-enter kernel.eventloop after catching SIGINT
* :ghpull:`1490`: rename plaintext cell -> raw cell
* :ghpull:`1480`: Fix %notebook magic, etc. nbformat unicode tests and fixes
* :ghpull:`1588`: Gtk3 integration with ipython works.
* :ghpull:`1595`: Examples syntax (avoid errors installing on Python 3)
* :ghpull:`1526`: Find encoding for Python files
* :ghpull:`1594`: Fix writing git commit ID to a file on build with Python 3
* :ghpull:`1556`: shallow-copy DictDB query results
* :ghpull:`1502`: small changes in response to pyflakes pass
* :ghpull:`1445`: Don't build sphinx docs for sdists
* :ghpull:`1538`: store git commit hash in utils._sysinfo instead of hidden data file
* :ghpull:`1546`: attempt to suppress exceptions in channel threads at shutdown
* :ghpull:`1559`: update tools/github_stats.py to use GitHub API v3
* :ghpull:`1563`: clear_output improvements
* :ghpull:`1560`: remove obsolete discussion of Twisted/trial from testing docs
* :ghpull:`1569`: BUG: qtconsole -- non-standard handling of \a and \b. [Fixes #1561]
* :ghpull:`1573`: BUG: Ctrl+C crashes wx pylab kernel in qtconsole.
* :ghpull:`1568`: fix PR #1567
* :ghpull:`1567`: Fix: openssh_tunnel did not parse port in `server`
* :ghpull:`1565`: fix AsyncResult.abort
* :ghpull:`1552`: use os.getcwdu in NotebookManager
* :ghpull:`1541`: display_pub flushes stdout/err
* :ghpull:`1544`: make MultiKernelManager.kernel_manager_class configurable
* :ghpull:`1517`: Fix indentation bug in IPython/lib/pretty.py
* :ghpull:`1519`: BUG: Include the name of the exception type in its pretty format.
* :ghpull:`1489`: Fix zero-copy push
* :ghpull:`1477`: fix dangling `buffer` in IPython.parallel.util
* :ghpull:`1514`: DOC: Fix references to IPython.lib.pretty instead of the old location
* :ghpull:`1481`: BUG: Improve placement of CallTipWidget
* :ghpull:`1496`: BUG: LBYL when clearing the output history on shutdown.
* :ghpull:`1508`: fix sorting profiles in clustermanager
* :ghpull:`1495`: BUG: Fix pretty-printing for overzealous objects
* :ghpull:`1472`: more general fix for #662
* :ghpull:`1483`: updated magic_history docstring
* :ghpull:`1383`: First version of cluster web service.
* :ghpull:`1398`: fix %tb after SyntaxError
* :ghpull:`1440`: Fix for failing testsuite when using --with-xml-coverage on windows.
* :ghpull:`1419`: Add %install_ext magic function.
* :ghpull:`1424`: Win32 shell interactivity
* :ghpull:`1468`: Simplify structure of a Job in the TaskScheduler
* :ghpull:`1447`: 1107 - Tab autocompletion can suggest invalid syntax
* :ghpull:`1469`: Fix typo in comment (insert space)
* :ghpull:`1463`: Fix completion when importing modules in the cwd.
* :ghpull:`1466`: Fix for issue #1437, unfriendly windows qtconsole error handling
* :ghpull:`1432`: Fix ipython directive
* :ghpull:`1465`: allow `ipython help subcommand` syntax
* :ghpull:`1416`: Conditional import of ctypes in inputhook
* :ghpull:`1462`: expedite parallel tests
* :ghpull:`1410`: Add javascript library and css stylesheet loading to JS class.
* :ghpull:`1448`: Fix for #875 Never build unicode error messages
* :ghpull:`1458`: use eval to uncan References
* :ghpull:`1450`: load mathjax from CDN via https
* :ghpull:`1451`: include heading level in JSON
* :ghpull:`1444`: Fix pyhton -> python typos
* :ghpull:`1414`: ignore errors in shell.var_expand
* :ghpull:`1430`: Fix for tornado check for tornado < 1.1.0
* :ghpull:`1413`: get_home_dir expands symlinks, adjust test accordingly
* :ghpull:`1385`: updated and prettified magic doc strings
* :ghpull:`1406`: Browser selection
* :ghpull:`1377`: Saving non-ascii history
* :ghpull:`1402`: fix symlinked /home issue for FreeBSD
* :ghpull:`1405`: Only monkeypatch xunit when the tests are run using it.
* :ghpull:`1395`: Xunit & KnownFailure
* :ghpull:`1396`: Fix for %tb magic.
* :ghpull:`1386`: Jsd3
* :ghpull:`1388`: Add simple support for running inside a virtualenv
* :ghpull:`1391`: Improve Hub/Scheduler when no engines are registered
* :ghpull:`1369`: load header with engine id when engine dies in TaskScheduler
* :ghpull:`1353`: Save notebook as script using unicode file handle.
* :ghpull:`1352`: Add '-m mod : run library module as a script' option.
* :ghpull:`1363`: Fix some minor color/style config issues in the qtconsole
* :ghpull:`1371`: Adds a quiet keyword to sync_imports
* :ghpull:`1387`: Fixing Cell menu to update cell type select box.
* :ghpull:`1296`: Wx gui example: fixes the broken example for `%gui wx`.
* :ghpull:`1372`: ipcontroller cleans up connection files unless reuse=True
* :ghpull:`1374`: remove calls to meaningless ZMQStream.on_err
* :ghpull:`1370`: allow draft76 websockets (Safari)
* :ghpull:`1368`: Ensure handler patterns are str, not unicode
* :ghpull:`1361`: Notebook bug fix branch
* :ghpull:`1364`: avoid jsonlib returning Decimal
* :ghpull:`1362`: Don't log complete contents of history replies, even in debug
* :ghpull:`1347`: fix weird magic completion in notebook
* :ghpull:`1346`: fixups for alternate URL prefix stuff
* :ghpull:`1336`: crack at making notebook.html use the layout.html template
* :ghpull:`1331`: RST and heading cells
* :ghpull:`1247`: fixes a bug causing extra newlines after comments.
* :ghpull:`1332`: notebook - allow prefixes in URL path.
* :ghpull:`1341`: Don't attempt to tokenize binary files for tracebacks
* :ghpull:`1334`: added key handler for control-s to notebook, seems to work pretty well
* :ghpull:`1338`: Fix see also in docstrings so API docs build
* :ghpull:`1335`: Notebook toolbar UI
* :ghpull:`1299`: made notebook.html extend layout.html
* :ghpull:`1318`: make Ctrl-D in qtconsole act same as in terminal (ready to merge)
* :ghpull:`1328`: Coverage
* :ghpull:`1206`: don't preserve fixConsole output in json
* :ghpull:`1330`: Add linewrapping to text cells (new feature in CodeMirror).
* :ghpull:`1309`: Inoculate clearcmd extension into %reset functionality
* :ghpull:`1327`: Updatecm2
* :ghpull:`1326`: Removing Ace edit capability.
* :ghpull:`1325`: forgotten selected_cell -> get_selected_cell
* :ghpull:`1316`: Pass subprocess test runners a suitable location for xunit output
* :ghpull:`1303`: Updatecm
* :ghpull:`1312`: minor heartbeat tweaks
* :ghpull:`1306`: Fix %prun input parsing for escaped characters (closes #1302)
* :ghpull:`1301`: New "Fix for issue #1202" based on current master.
* :ghpull:`1289`: Make autoreload extension work on Python 3.
* :ghpull:`1288`: Don't ask for confirmation when stdin isn't available
* :ghpull:`1294`: TaskScheduler.hwm default to 1 instead of 0
* :ghpull:`1283`: HeartMonitor.period should be an Integer
* :ghpull:`1264`: Aceify
* :ghpull:`1284`: a fix for GH 1269
* :ghpull:`1213`: BUG: Minor typo in history_console_widget.py
* :ghpull:`1267`: add NoDB for non-recording Hub
* :ghpull:`1222`: allow Reference as callable in map/apply
* :ghpull:`1257`: use self.kernel_manager_class in qtconsoleapp
* :ghpull:`1253`: set auto_create flag for notebook apps
* :ghpull:`1262`: Heartbeat no longer shares the app's Context
* :ghpull:`1229`: Fix display of SyntaxError in Python 3
* :ghpull:`1256`: Dewijmoize
* :ghpull:`1246`: Skip tests that require X, when importing pylab results in RuntimeError.
* :ghpull:`1211`: serve local files in notebook-dir
* :ghpull:`1224`: edit text cells on double-click instead of single-click
* :ghpull:`1187`: misc notebook: connection file cleanup, first heartbeat, startup flush
* :ghpull:`1207`: fix loadpy duplicating newlines
* :ghpull:`1129`: Unified setup.py
* :ghpull:`1199`: Reduce IPython.external.*
* :ghpull:`1218`: Added -q option to %prun for suppression of the output, along with editing the dochelp string.
* :ghpull:`1217`: Added -q option to %prun for suppression of the output, along with editing the dochelp string
* :ghpull:`1175`: core.completer: Clean up excessive and unused code.
* :ghpull:`1196`: docs: looks like a file path might have been accidentally pasted in the middle of a word
* :ghpull:`1190`: Fix link to Chris Fonnesbeck blog post about 0.11 highlights.

Issues (742):

* :ghissue:`1943`: add screenshot and link into releasenotes
* :ghissue:`1570`: [notebook] remove 'left panel' references from example.
* :ghissue:`1954`: update some example notebooks
* :ghissue:`2048`: move _encode_binary to jsonutil.encode_images
* :ghissue:`2050`: only add quotes around xunit-file on Windows
* :ghissue:`2047`: disable auto-scroll on mozilla
* :ghissue:`1258`: Magic %paste error
* :ghissue:`2015`: Fixes for %paste with special transformations
* :ghissue:`760`: Windows: test runner fails if repo path contains spaces
* :ghissue:`2046`: Iptest unicode
* :ghissue:`1939`: Namespaces
* :ghissue:`2042`: increase auto-scroll threshold to 100 lines
* :ghissue:`2043`: move RemoteError import to top-level
* :ghissue:`641`: In %magic help, remove duplicate aliases
* :ghissue:`2036`: %alias_magic
* :ghissue:`1968`: Proposal of icons for .ipynb files
* :ghissue:`825`: keyboardinterrupt crashes gtk gui when gtk.set_interactive is not available
* :ghissue:`1971`: Remove duplicate magics docs
* :ghissue:`2040`: Namespaces for cleaner public APIs
* :ghissue:`2039`: ipython parallel import exception
* :ghissue:`2035`: Getdefaultencoding test error with sympy 0.7.1_git
* :ghissue:`2037`: remove `ipython-qtconsole` gui-script
* :ghissue:`1516`: ipython-qtconsole script isn't installed for Python 2.x
* :ghissue:`1297`: "ipython -p sh" is in documentation but doesn't work
* :ghissue:`2038`: add extra clear warning to shell doc
* :ghissue:`1265`: please ship unminified js and css sources
* :ghissue:`2029`: Ship unminified js
* :ghissue:`1920`: Provide an easy way to override the Qt widget used by qtconsole
* :ghissue:`2007`: Add custom_control and custom_page_control variables to override the Qt widgets used by qtconsole
* :ghissue:`2009`: In %magic help, remove duplicate aliases
* :ghissue:`2033`: ipython parallel pushing and pulling recarrays
* :ghissue:`2034`: fix&test push/pull recarrays
* :ghissue:`2028`: Reduce unhelpful information shown by pinfo
* :ghissue:`1992`: Tab completion fails with many spaces in filename 
* :ghissue:`1885`: handle too old wx
* :ghissue:`2030`: check wxPython version in inputhook
* :ghissue:`2024`: Make interactive_usage a bit more rst friendly
* :ghissue:`2031`: disable ^C^C confirmation on Windows
* :ghissue:`2023`: Unicode test failure on OS X
* :ghissue:`2027`: match stdin encoding in frontend readline test
* :ghissue:`1901`: Windows: parallel test fails assert, leaves 14 python processes alive
* :ghissue:`2025`: Fix parallel test on WinXP - wait for resource cleanup.
* :ghissue:`1986`: Line magic function `%R` not found. (Rmagic)
* :ghissue:`1712`: test failure in ubuntu package daily build
* :ghissue:`1183`: 0.12 testsuite failures
* :ghissue:`2016`: BUG: test runner fails in Windows if filenames contain spaces.
* :ghissue:`1806`: Alternate upload methods in firefox
* :ghissue:`2019`: Windows: home directory expansion test fails
* :ghissue:`2020`: Fix home path expansion test in Windows.
* :ghissue:`2017`: Windows core test error - filename quoting
* :ghissue:`2021`: Fix Windows pathname issue in 'odd encoding' test.
* :ghissue:`1998`: call to nt.assert_true(path._writable_dir(home)) returns false in test_path.py
* :ghissue:`2022`: don't check writability in test for get_home_dir when HOME is undefined
* :ghissue:`1589`: Test failures and docs don't build on Mac OS X Lion
* :ghissue:`1996`: frontend test tweaks
* :ghissue:`2011`: Notebook server can't start cluster with hyphen-containing profile name
* :ghissue:`2014`: relax profile regex in notebook
* :ghissue:`2013`: brew install pyqt
* :ghissue:`2005`: Strange output artifacts in footer of notebook
* :ghissue:`2012`: Mono cursor offset
* :ghissue:`2004`: Clarify generic message spec vs. Python message API in docs
* :ghissue:`2006`: Don't crash when starting notebook server if runnable browser not found
* :ghissue:`2010`: notebook: Print a warning (but do not abort) if no webbrowser can be found.
* :ghissue:`2008`: pip install virtualenv
* :ghissue:`2003`: Wrong case of rmagic in docs
* :ghissue:`2002`: Refactor %magic into a lsmagic_docs API function.
* :ghissue:`2000`: kernel.js consistency with generic IPython message format.
* :ghissue:`1999`: `%magic` help: display line and cell magics in alphabetical order.
* :ghissue:`1635`: test_prun_quotes fails on Windows
* :ghissue:`1984`: Cannot restart Notebook when using `%%script --bg`
* :ghissue:`1981`: Clean BG processes created by %%script on kernel exit
* :ghissue:`1994`: Fix RST misformatting.
* :ghissue:`1949`: Introduce Notebook Magics
* :ghissue:`1985`: Kernels should start in notebook dir when manually specified
* :ghissue:`1980`: Notebook should check that --notebook-dir exists
* :ghissue:`1951`: minor notebook startup/notebook-dir adjustments
* :ghissue:`1969`: tab completion in notebook for paths not triggered
* :ghissue:`1974`: Allow path completion on notebook.
* :ghissue:`1964`: allow multiple instances of a Magic
* :ghissue:`1960`: %page not working
* :ghissue:`1991`: fix _ofind attr in %page
* :ghissue:`1982`: Shutdown qtconsole problem?
* :ghissue:`1988`: check for active frontend in update_restart_checkbox
* :ghissue:`1979`: Add support for tox (https://tox.readthedocs.io/) and Travis CI (http://travis-ci.org/)
* :ghissue:`1989`: Parallel: output of %px and %px${suffix} is inconsistent
* :ghissue:`1966`: ValueError: packer could not serialize a simple message
* :ghissue:`1987`: Notebook: MathJax offline install not recognized
* :ghissue:`1970`: dblclick to restore size of images
* :ghissue:`1983`: Notebook does not save heading level
* :ghissue:`1978`: Notebook names truncating at the first period
* :ghissue:`1553`: Limited size of output cells and provide scroll bars for such output cells
* :ghissue:`1825`: second attempt at scrolled long output
* :ghissue:`1915`: add cell-level metadata
* :ghissue:`1934`: Cell/Worksheet metadata
* :ghissue:`1746`: Confirm restart (configuration option, and checkbox UI)
* :ghissue:`1790`: Commenting function.
* :ghissue:`1767`: Tab completion problems with cell magics
* :ghissue:`1944`: [qtconsole] take %,%% prefix into account for completion
* :ghissue:`1973`: fix another FreeBSD $HOME symlink issue
* :ghissue:`1972`: Fix completion of '%tim' in the Qt console
* :ghissue:`1887`: Make it easy to resize jpeg/png images back to original size.
* :ghissue:`1967`: Fix psums example description in docs
* :ghissue:`1678`: ctrl-z clears cell output in notebook when pressed enough times
* :ghissue:`1965`: fix for #1678, undo no longer clears cells
* :ghissue:`1952`: avoid duplicate "Websockets closed" dialog on ws close
* :ghissue:`1961`: UnicodeDecodeError on directory with unicode chars in prompt
* :ghissue:`1963`: styling prompt, {color.Normal} excepts
* :ghissue:`1962`: Support unicode prompts
* :ghissue:`1959`: %page not working on qtconsole for Windows XP 32-bit
* :ghissue:`1955`: update to latest version of vim-ipython
* :ghissue:`1945`: Add --proc option to %%script
* :ghissue:`1957`: fix indentation in kernel.js
* :ghissue:`1956`: move import RemoteError after get_exc_info
* :ghissue:`1950`: Fix for copy action (Ctrl+C) when there is no pager defined in qtconsole
* :ghissue:`1948`: Fix help string for InteractiveShell.ast_node_interactivity
* :ghissue:`1941`: script magics cause terminal spam
* :ghissue:`1942`: swallow stderr of which in utils.process.find_cmd
* :ghissue:`1833`: completer draws slightly too small on Chrome
* :ghissue:`1940`: fix completer css on some Chrome versions
* :ghissue:`1938`: remove remaining references to deprecated XREP/XREQ names
* :ghissue:`1924`: HTML superscripts not shown raised in the notebook
* :ghissue:`1925`: Fix styling of superscripts and subscripts. Closes #1924.
* :ghissue:`1461`: User notification if notebook saving fails
* :ghissue:`1936`: increase duration of save messages
* :ghissue:`1542`: %save magic fails in clients without stdin if file already exists
* :ghissue:`1937`: add %save -f
* :ghissue:`1572`: pyreadline version dependency not correctly checked
* :ghissue:`1935`: add version checking to pyreadline import test
* :ghissue:`1849`: Octave magics
* :ghissue:`1759`: github, merge PR(s) just by number(s) 
* :ghissue:`1931`: Win py3fixes
* :ghissue:`1646`: Meaning of restart parameter in client.shutdown() unclear
* :ghissue:`1933`: oinspect.find_file: Additional safety if file cannot be found.
* :ghissue:`1916`: %paste doesn't work on py3
* :ghissue:`1932`: Fix adding functions to CommandChainDispatcher with equal priority on Py 3
* :ghissue:`1928`: Select NoDB by default
* :ghissue:`1923`: Add IPython syntax support to the %timeit magic, in line and cell mode
* :ghissue:`1926`: Make completer recognize escaped quotes in strings.
* :ghissue:`1929`: Ipython-qtconsole (0.12.1) hangs with Python 2.7.3, Windows 7 64 bit
* :ghissue:`1409`: [qtconsole] forward delete bring completion into current line
* :ghissue:`1922`: py3k compatibility for setupegg.py
* :ghissue:`1598`: document that sync_imports() can't handle "import foo as bar"
* :ghissue:`1893`: Update Parallel Magics and Exception Display
* :ghissue:`1890`: Docstrings for magics that use @magic_arguments are rendered wrong
* :ghissue:`1921`: magic_arguments: dedent but otherwise preserve indentation.
* :ghissue:`1919`: Use oinspect in CodeMagics._find_edit_target
* :ghissue:`1918`: don't warn in iptest if deathrow/quarantine are missing
* :ghissue:`1914`: %pdef failing on python3
* :ghissue:`1917`: Fix for %pdef on Python 3
* :ghissue:`1428`: Failing test that prun does not clobber string escapes
* :ghissue:`1913`: Fix for #1428
* :ghissue:`1911`: temporarily skip autoreload tests
* :ghissue:`1549`: autoreload extension crashes ipython
* :ghissue:`1908`: find_file errors on windows
* :ghissue:`1909`: Fix for #1908, use os.path.normcase for safe filename comparisons
* :ghissue:`1907`: py3compat fixes for %%script and tests
* :ghissue:`1904`: %%px? doesn't work, shows info for %px, general cell magic problem
* :ghissue:`1906`: ofind finds non-unique cell magics
* :ghissue:`1894`: Win64 binary install fails
* :ghissue:`1799`: Source file not found for magics
* :ghissue:`1845`: Fixes to inspection machinery for magics
* :ghissue:`1774`: Some magics seems broken
* :ghissue:`1586`: Clean up tight coupling between Notebook, CodeCell and Kernel Javascript objects
* :ghissue:`1632`: Win32 shell interactivity apparently broke qtconsole "cd" magic
* :ghissue:`1902`: Workaround fix for gh-1632; minimal revert of gh-1424
* :ghissue:`1900`: Cython libs
* :ghissue:`1503`: Cursor is offset in notebook in Chrome 17 on Linux
* :ghissue:`1426`: Qt console doesn't handle the `--gui` flag correctly.
* :ghissue:`1180`: Can't start IPython kernel in Spyder
* :ghissue:`581`: test IPython.zmq
* :ghissue:`1593`: Name embedded in notebook overrides filename
* :ghissue:`1899`: add ScriptMagics to class list for generated config
* :ghissue:`1618`: generate or minimize manpages
* :ghissue:`1898`: minimize manpages
* :ghissue:`1896`: Windows: apparently spurious warning 'Excluding nonexistent file' ... test_exampleip
* :ghissue:`1897`: use glob for bad exclusion warning
* :ghissue:`1215`: updated %quickref to show short-hand for %sc and %sx
* :ghissue:`1855`: %%script and %%file magics
* :ghissue:`1863`: Ability to silence a cell in the notebook
* :ghissue:`1870`: add %%capture for capturing stdout/err
* :ghissue:`1861`: Use dvipng to format sympy.Matrix
* :ghissue:`1867`: Fix 1px margin bouncing of selected menu item.
* :ghissue:`1889`: Reconnect when the websocket connection closes unexpectedly
* :ghissue:`1577`: If a notebook loses its network connection WebSockets won't reconnect
* :ghissue:`1886`: Fix a bug in renaming notebook
* :ghissue:`1895`: Fix error in test suite with ip.system()
* :ghissue:`1762`: add `locate` entry points
* :ghissue:`1883`: Fix vertical offset due to bold/italics, and bad browser fonts.
* :ghissue:`1875`: re-write columnize, with intermediate step.
* :ghissue:`1860`: IPython.utils.columnize sometime wrong...
* :ghissue:`1851`: new completer for qtconsole.
* :ghissue:`1892`: Remove suspicious quotes in interactiveshell.py
* :ghissue:`1854`: Class `%hierarchy` and graphiz `%%dot` magics
* :ghissue:`1827`: Sending tracebacks over ZMQ should protect against unicode failure
* :ghissue:`1864`: Rmagic exceptions
* :ghissue:`1829`: [notebook] don't care about leading prct in completion
* :ghissue:`1832`: Make svg, jpeg and png images resizable in notebook.
* :ghissue:`1674`: HTML Notebook carriage-return handling, take 2
* :ghissue:`1874`: cython_magic uses importlib, which doesn't ship with py2.6
* :ghissue:`1882`: Remove importlib dependency which not available in Python 2.6.
* :ghissue:`1878`: shell access using ! will not fill class or function scope vars
* :ghissue:`1879`: Correct stack depth for variable expansion in !system commands
* :ghissue:`1840`: New JS completer should merge completions before display
* :ghissue:`1841`: [notebook] deduplicate completion results
* :ghissue:`1736`: no good error message on missing tkinter and %paste
* :ghissue:`1741`: Display message from TryNext error in magic_paste
* :ghissue:`1850`: Remove args/kwargs handling in TryNext, fix %paste error messages.
* :ghissue:`1663`: Keep line-endings in ipynb
* :ghissue:`1872`: Matplotlib window freezes using intreractive plot in qtconsole
* :ghissue:`1869`: Improve CodeMagics._find_edit_target
* :ghissue:`1781`: Colons in notebook name causes notebook deletion without warning
* :ghissue:`1815`: Make : invalid in filenames in the Notebook JS code.
* :ghissue:`1819`: doc: cleanup the parallel psums example a little
* :ghissue:`1838`: externals cleanup
* :ghissue:`1839`: External cleanup
* :ghissue:`1782`: fix Magic menu in qtconsole, split in groups
* :ghissue:`1862`: Minor bind_kernel improvements
* :ghissue:`1859`: kernmagic during console startup
* :ghissue:`1857`: Prevent jumping of window to input when output is clicked.
* :ghissue:`1856`: Fix 1px jumping of cells and menus in Notebook.
* :ghissue:`1848`: task fails with "AssertionError: not enough buffers!" after second resubmit
* :ghissue:`1852`: fix chained resubmissions
* :ghissue:`1780`: Rmagic extension
* :ghissue:`1853`: Fix jumpy notebook behavior
* :ghissue:`1842`: task with UnmetDependency error still owned by engine
* :ghissue:`1847`: add InlineBackend to ConsoleApp class list
* :ghissue:`1846`: Exceptions within multiprocessing crash Ipython notebook kernel
* :ghissue:`1843`: Notebook does not exist and permalinks
* :ghissue:`1837`: edit magic broken in head
* :ghissue:`1834`: resubmitted tasks doesn't have same session name
* :ghissue:`1836`: preserve header for resubmitted tasks
* :ghissue:`1776`: fix magic menu in qtconsole
* :ghissue:`1828`: change default extension to .ipy for %save -r
* :ghissue:`1800`: Reintroduce recall
* :ghissue:`1671`: __future__ environments
* :ghissue:`1830`: lsmagic lists magics in alphabetical order
* :ghissue:`1835`: Use Python import in ipython profile config
* :ghissue:`1773`: Update SymPy profile: SymPy's latex() can now print set and frozenset
* :ghissue:`1761`: Edited documentation to use IPYTHONDIR in place of ~/.ipython
* :ghissue:`1772`: notebook autocomplete fail when typing number
* :ghissue:`1822`: aesthetics pass on AsyncResult.display_outputs
* :ghissue:`1460`: Redirect http to https for notebook
* :ghissue:`1287`: Refactor the notebook tab completion/tooltip
* :ghissue:`1596`: In rename dialog, <return> should submit
* :ghissue:`1821`: ENTER submits the rename notebook dialog.
* :ghissue:`1750`: Let the user disable random port selection
* :ghissue:`1820`: NotebookApp: Make the number of ports to retry user configurable.
* :ghissue:`1816`: Always use filename as the notebook name.
* :ghissue:`1775`: assert_in not present on Python 2.6
* :ghissue:`1813`: Add assert_in method to nose for Python 2.6
* :ghissue:`1498`: Add tooltip keyboard shortcuts
* :ghissue:`1711`: New Tooltip, New Completer and JS Refactor
* :ghissue:`1798`: a few simple fixes for docs/parallel
* :ghissue:`1818`: possible bug with latex / markdown
* :ghissue:`1647`: Aborted parallel tasks can't be resubmitted
* :ghissue:`1817`: Change behavior of ipython notebook --port=...
* :ghissue:`1738`: IPython.embed_kernel issues
* :ghissue:`1610`: Basic bold and italic in HTML output cells
* :ghissue:`1576`: Start and stop kernels from the notebook dashboard
* :ghissue:`1515`: impossible to shutdown notebook kernels
* :ghissue:`1812`: Ensure AsyncResult.display_outputs doesn't display empty streams
* :ghissue:`1811`: warn on nonexistent exclusions in iptest
* :ghissue:`1809`: test suite error in IPython.zmq on windows
* :ghissue:`1810`: fix for #1809, failing tests in IPython.zmq
* :ghissue:`1808`: Reposition alternate upload for firefox [need cross browser/OS/language test]
* :ghissue:`1742`: Check for custom_exceptions only once
* :ghissue:`1802`: cythonmagic tests should be skipped if Cython not available
* :ghissue:`1062`: warning message in IPython.extensions test
* :ghissue:`1807`: add missing cython exclusion in iptest
* :ghissue:`1805`: Fixed a vcvarsall.bat error on win32/Py2.7 when trying to compile with m...
* :ghissue:`1803`: MPI parallel %px bug 
* :ghissue:`1804`: Fixed a vcvarsall.bat error on win32/Py2.7 when trying to compile with mingw.
* :ghissue:`1492`: Drag target very small if IPython Dashboard has no notebooks
* :ghissue:`1562`: Offer a method other than drag-n-drop to upload notebooks
* :ghissue:`1739`: Dashboard improvement (necessary merge of #1658 and #1676 + fix #1492)
* :ghissue:`1770`: Cython related magic functions
* :ghissue:`1532`: qtconsole does not accept --gui switch
* :ghissue:`1707`: Accept --gui=<...> switch in IPython qtconsole.
* :ghissue:`1797`: Fix comment which breaks Emacs syntax highlighting.
* :ghissue:`1796`: %gui magic broken
* :ghissue:`1795`: fix %gui magic
* :ghissue:`1788`: extreme truncating of return values
* :ghissue:`1793`: Raise repr limit for strings to 80 characters (from 30).
* :ghissue:`1794`: don't use XDG path on OS X
* :ghissue:`1777`: ipython crash on wrong encoding
* :ghissue:`1792`: Unicode-aware logger
* :ghissue:`1791`: update zmqshell magics
* :ghissue:`1787`: DOC: Remove regression from qt-console docs.
* :ghissue:`1785`: IPython.utils.tests.test_process.SubProcessTestCase
* :ghissue:`1758`: test_pr, fallback on http if git protocol fail, and SSL errors...
* :ghissue:`1786`: Make notebook save failures more salient
* :ghissue:`1748`: Fix some tests for Python 3.3
* :ghissue:`1755`: test for pygments before running qt tests
* :ghissue:`1771`: Make default value of interactivity passed to run_ast_nodes configurable
* :ghissue:`1783`: part of PR #1606 (loadpy -> load) erased by magic refactoring.
* :ghissue:`1784`: restore loadpy to load
* :ghissue:`1768`: Update parallel magics
* :ghissue:`1778`: string exception in IPython/core/magic.py:232
* :ghissue:`1779`: Tidy up error raising in magic decorators.
* :ghissue:`1769`: Allow cell mode timeit without setup code.
* :ghissue:`1716`: Fix for fake filenames in verbose traceback
* :ghissue:`1763`: [qtconsole] fix append_plain_html -> append_html
* :ghissue:`1766`: Test failure in IPython.parallel
* :ghissue:`1611`: IPEP1: Cell magics and general cleanup of the Magic system
* :ghissue:`1732`: Refactoring of the magics system and implementation of cell magics
* :ghissue:`1765`: test_pr should clearn PYTHONPATH for the subprocesses
* :ghissue:`1630`: Merge divergent Kernel implementations
* :ghissue:`1705`: [notebook] Make pager resizable, and remember size...
* :ghissue:`1606`: Share code for %pycat and %loadpy, make %pycat aware of URLs
* :ghissue:`1720`: Adding interactive inline plotting to notebooks with flot
* :ghissue:`1701`: [notebook] Open HTML links in a new window by default
* :ghissue:`1757`: Open IPython notebook hyperlinks in a new window using target=_blank
* :ghissue:`1735`: Open IPython notebook hyperlinks in a new window using target=_blank
* :ghissue:`1754`: Fix typo enconters->encounters
* :ghissue:`1753`: Clear window title when kernel is restarted
* :ghissue:`735`: Images missing from XML/SVG export (for me)
* :ghissue:`1449`: Fix for bug #735 : Images missing from XML/SVG export
* :ghissue:`1752`: Reconnect Websocket when it closes unexpectedly
* :ghissue:`1751`: Reconnect Websocket when it closes unexpectedly
* :ghissue:`1749`: Load MathJax.js using HTTPS when IPython notebook server is HTTPS
* :ghissue:`1743`: Tooltip completer js refactor
* :ghissue:`1700`: A module for sending custom user messages from the kernel.
* :ghissue:`1745`: htmlnotebook: Cursor is off
* :ghissue:`1728`: ipython crash with matplotlib during picking
* :ghissue:`1681`: add qt config option to clear_on_kernel_restart
* :ghissue:`1733`: Tooltip completer js refactor
* :ghissue:`1676`: Kernel status/shutdown from dashboard
* :ghissue:`1658`: Alternate notebook upload methods
* :ghissue:`1727`: terminate kernel after embed_kernel tests
* :ghissue:`1737`: add HistoryManager to ipapp class list
* :ghissue:`945`: Open a notebook from the command line
* :ghissue:`1686`: ENH: Open a notebook from the command line
* :ghissue:`1709`: fixes #1708, failing test in arg_split on windows
* :ghissue:`1718`: Use CRegExp trait for regular expressions.
* :ghissue:`1729`: Catch failure in repr() for %whos
* :ghissue:`1726`: use eval for command-line args instead of exec
* :ghissue:`1723`: scatter/gather fail with targets='all'
* :ghissue:`1724`: fix scatter/gather with targets='all'
* :ghissue:`1725`: add --no-ff to git pull in test_pr
* :ghissue:`1722`: unicode exception when evaluating expression with non-ascii characters
* :ghissue:`1721`: Tooltip completer js refactor
* :ghissue:`1657`: Add `wait` optional argument to `hooks.editor`
* :ghissue:`123`: Define sys.ps{1,2}
* :ghissue:`1717`: Define generic sys.ps{1,2,3}, for use by scripts.
* :ghissue:`1442`: cache-size issue in qtconsole
* :ghissue:`1691`: Finish PR #1446
* :ghissue:`1446`: Fixing Issue #1442
* :ghissue:`1710`: update MathJax CDN url for https
* :ghissue:`81`: Autocall fails if first function argument begins with "-" or "+
* :ghissue:`1713`: Make autocall regexp's configurable.
* :ghissue:`211`: paste command not working
* :ghissue:`1703`: Allow TryNext to have an error message without it affecting the command chain
* :ghissue:`1714`: minor adjustments to test_pr
* :ghissue:`1509`: New tooltip for notebook
* :ghissue:`1697`: Major refactoring of the Notebook, Kernel and CodeCell JavaScript.
* :ghissue:`788`: Progress indicator in the notebook (and perhaps the Qt console)
* :ghissue:`1034`: Single process Qt console
* :ghissue:`1557`: magic function conflict while using --pylab
* :ghissue:`1476`: Pylab figure objects not properly updating
* :ghissue:`1704`: ensure all needed qt parts can be imported before settling for one
* :ghissue:`1708`: test failure in arg_split on windows
* :ghissue:`1706`: Mark test_push_numpy_nocopy as a known failure for Python 3
* :ghissue:`1696`: notebook tooltip fail on function with number
* :ghissue:`1698`: fix tooltip on token with number
* :ghissue:`1226`: Windows GUI only (pythonw) bug for IPython on Python 3.x
* :ghissue:`1245`: pythonw py3k fixes for issue #1226
* :ghissue:`1417`: Notebook Completer Class
* :ghissue:`1690`: [Bogus] Deliberately make a test fail
* :ghissue:`1685`: Add script to test pull request
* :ghissue:`1167`: Settle on a choice for $IPYTHONDIR
* :ghissue:`1693`: deprecate IPYTHON_DIR in favor of IPYTHONDIR
* :ghissue:`1672`: ipython-qtconsole.desktop is using a deprecated format
* :ghissue:`1695`: Avoid deprecated warnings from ipython-qtconsole.desktop.
* :ghissue:`1694`: Add quote to notebook to allow it to load
* :ghissue:`1240`: sys.path missing `''` as first entry when kernel launched without interface
* :ghissue:`1689`: Fix sys.path missing '' as first entry in `ipython kernel`.
* :ghissue:`1683`: Parallel controller failing with Pymongo 2.2
* :ghissue:`1687`: import Binary from bson instead of pymongo
* :ghissue:`1614`: Display Image in Qtconsole
* :ghissue:`1616`: Make IPython.core.display.Image less notebook-centric
* :ghissue:`1684`: CLN: Remove redundant function definition.
* :ghissue:`1655`: Add %open magic command to open editor in non-blocking manner
* :ghissue:`1677`: middle-click paste broken in notebook
* :ghissue:`1670`: Point %pastebin to gist
* :ghissue:`1667`: Test failure in test_message_spec
* :ghissue:`1668`: Test failure in IPython.zmq.tests.test_message_spec.test_complete "'pyout' != 'status'"
* :ghissue:`1669`: handle pyout messages in test_message_spec
* :ghissue:`1295`: add binary-tree engine interconnect example
* :ghissue:`1642`: Cherry-picked commits from 0.12.1 release
* :ghissue:`1659`: Handle carriage return characters ("\r") in HTML notebook output.
* :ghissue:`1313`: Figure out MathJax 2 support
* :ghissue:`1653`: Test failure in IPython.zmq
* :ghissue:`1656`: ensure kernels are cleaned up in embed_kernel tests
* :ghissue:`1666`: pip install ipython==dev installs version 0.8 from an old svn repo
* :ghissue:`1664`: InteractiveShell.run_code: Update docstring.
* :ghissue:`1512`: `print stuff,` should avoid newline
* :ghissue:`1662`: Delay flushing softspace until after cell finishes
* :ghissue:`1643`: handle jpg/jpeg in the qtconsole
* :ghissue:`966`: dreload fails on Windows XP with IPython 0.11 "Unexpected Error"
* :ghissue:`1500`: dreload doesn't seem to exclude numpy
* :ghissue:`1520`: kernel crash when showing tooltip (?)
* :ghissue:`1652`: add patch_pyzmq() for backporting a few changes from newer pyzmq
* :ghissue:`1650`: DOC: moving files with SSH launchers
* :ghissue:`1357`: add IPython.embed_kernel() 
* :ghissue:`1640`: Finish up embed_kernel
* :ghissue:`1651`: Remove bundled Itpl module
* :ghissue:`1634`: incremental improvements to SSH launchers
* :ghissue:`1649`: move examples/test_embed into examples/tests/embed
* :ghissue:`1171`: Recognise virtualenvs
* :ghissue:`1479`: test_extension failing in Windows
* :ghissue:`1633`: Fix installing extension from local file on Windows
* :ghissue:`1644`: Update copyright date to 2012
* :ghissue:`1636`: Test_deepreload breaks pylab irunner tests
* :ghissue:`1645`: Exclude UserDict when deep reloading NumPy.
* :ghissue:`1454`: make it possible to start engine in 'disabled' mode and 'enable' later
* :ghissue:`1641`: Escape code for the current time in PromptManager
* :ghissue:`1638`: ipython console clobbers custom sys.path
* :ghissue:`1637`: Removed a ':' which shouldn't have been there
* :ghissue:`1536`: ipython 0.12 embed shell won't run startup scripts
* :ghissue:`1628`: error: QApplication already exists in TestKillRing
* :ghissue:`1631`: TST: QApplication doesn't quit early enough with PySide.
* :ghissue:`1629`: evaluate a few dangling validate_message generators
* :ghissue:`1621`: clear In[] prompt numbers on "Clear All Output"
* :ghissue:`1627`: Test the Message Spec
* :ghissue:`1470`: SyntaxError on setup.py install with Python 3
* :ghissue:`1624`: Fixes for byte-compilation on Python 3
* :ghissue:`1612`: pylab=inline fig.show() non-existent in notebook
* :ghissue:`1615`: Add show() method to figure objects.
* :ghissue:`1622`: deepreload fails on Python 3
* :ghissue:`1625`: Fix deepreload on Python 3
* :ghissue:`1626`: Failure in new `dreload` tests under Python 3.2
* :ghissue:`1623`: IPython / matplotlib Memory error with imshow
* :ghissue:`1619`: pyin messages should have execution_count
* :ghissue:`1620`: pyin message now have execution_count
* :ghissue:`32`: dreload produces spurious traceback when numpy is involved
* :ghissue:`1457`: Update deepreload to use a rewritten knee.py. Fixes dreload(numpy).
* :ghissue:`1613`: allow map / parallel function for single-engine views
* :ghissue:`1609`: exit notebook cleanly on SIGINT, SIGTERM
* :ghissue:`1531`: Function keyword completion fails if cursor is in the middle of the complete parentheses
* :ghissue:`1607`: cleanup sqlitedb temporary db file after tests
* :ghissue:`1608`: don't rely on timedelta.total_seconds in AsyncResult
* :ghissue:`1421`: ipython32 %run -d breaks with NameError name 'execfile' is not defined
* :ghissue:`1599`: Fix for %run -d on Python 3
* :ghissue:`1201`: %env magic fails with Python 3.2
* :ghissue:`1602`: Fix %env magic on Python 3.
* :ghissue:`1603`: Remove python3 profile
* :ghissue:`1604`: Exclude IPython.quarantine from installation
* :ghissue:`1601`: Security file is not removed after shutdown by Ctrl+C or kill -INT
* :ghissue:`1600`: Specify encoding for io.open in notebook_reformat tests
* :ghissue:`1605`: Small fixes for Animation and Progress notebook
* :ghissue:`1452`: Bug fix for approval
* :ghissue:`13`: Improve robustness and debuggability of test suite
* :ghissue:`70`: IPython should prioritize __all__ during tab completion
* :ghissue:`1529`: __all__ feature, improvement to dir2, and tests for both
* :ghissue:`1475`: Custom namespace for %run
* :ghissue:`1564`: calling .abort on AsyncMapResult  results in traceback
* :ghissue:`1548`: add sugar methods/properties to AsyncResult
* :ghissue:`1535`: Fix pretty printing dispatch
* :ghissue:`1522`: Discussion: some potential Qt console refactoring
* :ghissue:`1399`: Use LaTeX to print various built-in types with the SymPy printing extension
* :ghissue:`1597`: re-enter kernel.eventloop after catching SIGINT
* :ghissue:`1490`: rename plaintext cell -> raw cell
* :ghissue:`1487`: %notebook fails in qtconsole
* :ghissue:`1545`: trailing newline not preserved in splitline ipynb
* :ghissue:`1480`: Fix %notebook magic, etc. nbformat unicode tests and fixes
* :ghissue:`1588`: Gtk3 integration with ipython works.
* :ghissue:`1595`: Examples syntax (avoid errors installing on Python 3)
* :ghissue:`1526`: Find encoding for Python files
* :ghissue:`1594`: Fix writing git commit ID to a file on build with Python 3
* :ghissue:`1556`: shallow-copy DictDB query results
* :ghissue:`1499`: various pyflakes issues
* :ghissue:`1502`: small changes in response to pyflakes pass
* :ghissue:`1445`: Don't build sphinx docs for sdists
* :ghissue:`1484`: unhide .git_commit_info.ini
* :ghissue:`1538`: store git commit hash in utils._sysinfo instead of hidden data file
* :ghissue:`1546`: attempt to suppress exceptions in channel threads at shutdown
* :ghissue:`1524`: unhide git_commit_info.ini
* :ghissue:`1559`: update tools/github_stats.py to use GitHub API v3
* :ghissue:`1563`: clear_output improvements
* :ghissue:`1558`: Ipython testing documentation still mentions twisted and trial
* :ghissue:`1560`: remove obsolete discussion of Twisted/trial from testing docs
* :ghissue:`1561`: Qtconsole - nonstandard \a and \b
* :ghissue:`1569`: BUG: qtconsole -- non-standard handling of \a and \b. [Fixes #1561]
* :ghissue:`1574`: BUG: Ctrl+C crashes wx pylab kernel in qtconsole
* :ghissue:`1573`: BUG: Ctrl+C crashes wx pylab kernel in qtconsole.
* :ghissue:`1590`: 'IPython3 qtconsole' doesn't work in Windows 7
* :ghissue:`602`: User test the html notebook
* :ghissue:`613`: Implement Namespace panel section
* :ghissue:`879`: How to handle Javascript output in the notebook
* :ghissue:`1255`: figure.show() raises an error with the inline backend
* :ghissue:`1467`: Document or bundle a git-integrated facility for stripping VCS-unfriendly binary data
* :ghissue:`1237`: Kernel status and logout button overlap
* :ghissue:`1319`: Running a cell with ctrl+Enter selects text in cell
* :ghissue:`1571`: module member autocomplete should respect __all__
* :ghissue:`1566`: ipython3 doesn't run in Win7 with Python 3.2 
* :ghissue:`1568`: fix PR #1567
* :ghissue:`1567`: Fix: openssh_tunnel did not parse port in `server`
* :ghissue:`1565`: fix AsyncResult.abort
* :ghissue:`1550`: Crash when starting notebook in a non-ascii path
* :ghissue:`1552`: use os.getcwdu in NotebookManager
* :ghissue:`1554`: wrong behavior of the all function on iterators
* :ghissue:`1541`: display_pub flushes stdout/err
* :ghissue:`1539`: Asynchrous issue when using clear_display and print x,y,z
* :ghissue:`1544`: make MultiKernelManager.kernel_manager_class configurable
* :ghissue:`1494`: Untrusted Secure Websocket broken on latest chrome dev
* :ghissue:`1521`: only install ipython-qtconsole gui script on Windows
* :ghissue:`1528`: Tab completion optionally respects __all__ (+ dir2() cleanup)
* :ghissue:`1527`: Making a progress bar work in IPython Notebook
* :ghissue:`1497`: __all__ functionality added to dir2(obj)
* :ghissue:`1518`: Pretty printing exceptions is broken
* :ghissue:`811`: Fixes for ipython unhandeled OSError exception on failure of os.getcwdu()
* :ghissue:`1517`: Fix indentation bug in IPython/lib/pretty.py
* :ghissue:`1519`: BUG: Include the name of the exception type in its pretty format.
* :ghissue:`1525`: A hack for auto-complete numpy recarray
* :ghissue:`1489`: Fix zero-copy push
* :ghissue:`1401`: numpy arrays cannot be used with View.apply() in Python 3
* :ghissue:`1477`: fix dangling `buffer` in IPython.parallel.util
* :ghissue:`1514`: DOC: Fix references to IPython.lib.pretty instead of the old location
* :ghissue:`1511`: Version comparison error ( '2.1.11' < '2.1.4' ==> True)
* :ghissue:`1506`: "Fixing" the Notebook scroll to help in visually comparing outputs
* :ghissue:`1481`: BUG: Improve placement of CallTipWidget
* :ghissue:`1241`: When our debugger class is used standalone `_oh` key errors are thrown
* :ghissue:`676`: IPython.embed() from ipython crashes twice on exit
* :ghissue:`1496`: BUG: LBYL when clearing the output history on shutdown.
* :ghissue:`1507`: python3 notebook: TypeError: unorderable types
* :ghissue:`1508`: fix sorting profiles in clustermanager
* :ghissue:`1495`: BUG: Fix pretty-printing for overzealous objects
* :ghissue:`1505`: SQLite objects created in a thread can only be used in that same thread
* :ghissue:`1482`: %history documentation out of date?
* :ghissue:`1501`: dreload doesn't seem to exclude numpy
* :ghissue:`1472`: more general fix for #662
* :ghissue:`1486`: save state of qtconsole
* :ghissue:`1485`: add history search to qtconsole
* :ghissue:`1483`: updated magic_history docstring
* :ghissue:`1383`: First version of cluster web service.
* :ghissue:`482`: test_run.test_tclass fails on Windows
* :ghissue:`1398`: fix %tb after SyntaxError
* :ghissue:`1478`: key function or lambda in sorted function doesn't find global variables
* :ghissue:`1415`: handle exit/quit/exit()/quit() variants in zmqconsole
* :ghissue:`1440`: Fix for failing testsuite when using --with-xml-coverage on windows.
* :ghissue:`1419`: Add %install_ext magic function.
* :ghissue:`1424`: Win32 shell interactivity
* :ghissue:`1434`: Controller should schedule tasks of multiple clients at the same time
* :ghissue:`1268`: notebook %reset magic fails with StdinNotImplementedError
* :ghissue:`1438`: from cherrypy import expose fails when running script form parent directory
* :ghissue:`1468`: Simplify structure of a Job in the TaskScheduler
* :ghissue:`875`: never build unicode error messages
* :ghissue:`1107`: Tab autocompletion can suggest invalid syntax
* :ghissue:`1447`: 1107 - Tab autocompletion can suggest invalid syntax
* :ghissue:`1469`: Fix typo in comment (insert space)
* :ghissue:`1463`: Fix completion when importing modules in the cwd.
* :ghissue:`1437`: unfriendly error handling with pythonw and ipython-qtconsole
* :ghissue:`1466`: Fix for issue #1437, unfriendly windows qtconsole error handling
* :ghissue:`1432`: Fix ipython directive
* :ghissue:`1465`: allow `ipython help subcommand` syntax
* :ghissue:`1394`: Wishlist: Remove hard dependency on ctypes
* :ghissue:`1416`: Conditional import of ctypes in inputhook
* :ghissue:`1462`: expedite parallel tests
* :ghissue:`1418`: Strict mode in javascript
* :ghissue:`1410`: Add javascript library and css stylesheet loading to JS class.
* :ghissue:`1427`: #922 again
* :ghissue:`1448`: Fix for #875 Never build unicode error messages
* :ghissue:`1458`: use eval to uncan References
* :ghissue:`1455`: Python3 install fails
* :ghissue:`1450`: load mathjax from CDN via https
* :ghissue:`1182`: Qtconsole, multiwindow
* :ghissue:`1439`: Notebook not storing heading celltype information
* :ghissue:`1451`: include heading level in JSON
* :ghissue:`1444`: Fix pyhton -> python typos
* :ghissue:`1412`: Input parsing issue with %prun
* :ghissue:`1414`: ignore errors in shell.var_expand
* :ghissue:`1441`: (1) Enable IPython.notebook.kernel.execute to publish display_* even it is not called with a code cell and (2) remove empty html element when execute "display_*"
* :ghissue:`1431`: Beginner Error: ipython qtconsole
* :ghissue:`1436`: "ipython-qtconsole --gui qt" hangs on 64-bit win7
* :ghissue:`1433`: websocket connection fails on Chrome
* :ghissue:`1430`: Fix for tornado check for tornado < 1.1.0
* :ghissue:`1408`: test_get_home_dir_3 failed on Mac OS X
* :ghissue:`1413`: get_home_dir expands symlinks, adjust test accordingly
* :ghissue:`1420`: fixes #922
* :ghissue:`823`: KnownFailure tests appearing as errors
* :ghissue:`1385`: updated and prettified magic doc strings
* :ghissue:`1406`: Browser selection
* :ghissue:`1411`: ipcluster starts 8 engines "successfully" but Client only finds two
* :ghissue:`1375`: %history -g -f file encoding issue
* :ghissue:`1377`: Saving non-ascii history
* :ghissue:`797`: Source introspection needs to be smarter in python 3.2
* :ghissue:`846`: Autoreload extension doesn't work with Python 3.2
* :ghissue:`1360`: IPython notebook not starting on winXP
* :ghissue:`1407`: Qtconsole segfaults on OSX when displaying some pop-up function tooltips
* :ghissue:`1402`: fix symlinked /home issue for FreeBSD
* :ghissue:`1403`: pyreadline cyclic dependency with pdb++/pdbpp module
* :ghissue:`1405`: Only monkeypatch xunit when the tests are run using it.
* :ghissue:`1404`: Feature Request: List/Dictionary tab completion
* :ghissue:`1395`: Xunit & KnownFailure
* :ghissue:`1396`: Fix for %tb magic.
* :ghissue:`1397`: Stay or leave message not working, Safari session lost.
* :ghissue:`1389`: pylab=inline inoperant through ssh tunnelling?
* :ghissue:`1386`: Jsd3
* :ghissue:`1388`: Add simple support for running inside a virtualenv
* :ghissue:`826`: Add support for creation of parallel task when no engine is running
* :ghissue:`1391`: Improve Hub/Scheduler when no engines are registered
* :ghissue:`1369`: load header with engine id when engine dies in TaskScheduler
* :ghissue:`1345`: notebook can't save unicode as script
* :ghissue:`1353`: Save notebook as script using unicode file handle.
* :ghissue:`1352`: Add '-m mod : run library module as a script' option.
* :ghissue:`1363`: Fix some minor color/style config issues in the qtconsole
* :ghissue:`1371`: Adds a quiet keyword to sync_imports
* :ghissue:`1390`: Blank screen for notebooks on Safari
* :ghissue:`1387`: Fixing Cell menu to update cell type select box.
* :ghissue:`645`: Standalone WX GUI support is broken
* :ghissue:`1296`: Wx gui example: fixes the broken example for `%gui wx`.
* :ghissue:`1254`: typo in notebooklist.js breaks links
* :ghissue:`781`: Users should be able to clone a notebook
* :ghissue:`1372`: ipcontroller cleans up connection files unless reuse=True
* :ghissue:`1374`: remove calls to meaningless ZMQStream.on_err
* :ghissue:`1382`: Update RO for Notebook
* :ghissue:`1370`: allow draft76 websockets (Safari)
* :ghissue:`1368`: Ensure handler patterns are str, not unicode
* :ghissue:`1379`: Sage link on website homepage broken
* :ghissue:`1376`: FWIW does not work with Chrome 16.0.912.77 Ubuntu 10.10
* :ghissue:`1358`: Cannot install ipython on Windows 7 64-bit
* :ghissue:`1367`: Ctrl - m  t does not toggle output in chrome
* :ghissue:`1359`: [sympyprinting] MathJax can't render \root{m}{n}
* :ghissue:`1337`: Tab in the notebook after `(` should not indent, only give a tooltip
* :ghissue:`1339`: Notebook printing broken
* :ghissue:`1344`: Ctrl + M + L does not toggle line numbering in htmlnotebook
* :ghissue:`1348`: Ctrl + M + M does not switch to markdown cell
* :ghissue:`1361`: Notebook bug fix branch
* :ghissue:`1364`: avoid jsonlib returning Decimal
* :ghissue:`1362`: Don't log complete contents of history replies, even in debug
* :ghissue:`888`: ReST support in notebooks
* :ghissue:`1205`: notebook stores HTML escaped text in the file
* :ghissue:`1351`: add IPython.embed_kernel() 
* :ghissue:`1243`: magic commands without % are not completed properly in htmlnotebook
* :ghissue:`1347`: fix weird magic completion in notebook
* :ghissue:`1355`: notebook.html extends layout.html now
* :ghissue:`1354`: min and max in the notebook
* :ghissue:`1346`: fixups for alternate URL prefix stuff
* :ghissue:`1336`: crack at making notebook.html use the layout.html template
* :ghissue:`1331`: RST and heading cells
* :ghissue:`1350`: Add '-m mod : run library module as a script' option
* :ghissue:`1247`: fixes a bug causing extra newlines after comments.
* :ghissue:`1329`: add base_url to notebook configuration options
* :ghissue:`1332`: notebook - allow prefixes in URL path.
* :ghissue:`1317`: Very slow traceback construction from Cython extension
* :ghissue:`1341`: Don't attempt to tokenize binary files for tracebacks
* :ghissue:`1300`: Cell Input collapse
* :ghissue:`1334`: added key handler for control-s to notebook, seems to work pretty well
* :ghissue:`1338`: Fix see also in docstrings so API docs build
* :ghissue:`1335`: Notebook toolbar UI
* :ghissue:`1299`: made notebook.html extend layout.html
* :ghissue:`1318`: make Ctrl-D in qtconsole act same as in terminal (ready to merge)
* :ghissue:`873`: ReST support in notebook frontend
* :ghissue:`1139`: Notebook webkit notification
* :ghissue:`1314`: Insertcell
* :ghissue:`1328`: Coverage
* :ghissue:`1206`: don't preserve fixConsole output in json
* :ghissue:`1330`: Add linewrapping to text cells (new feature in CodeMirror).
* :ghissue:`1309`: Inoculate clearcmd extension into %reset functionality
* :ghissue:`1327`: Updatecm2
* :ghissue:`1326`: Removing Ace edit capability.
* :ghissue:`1325`: forgotten selected_cell -> get_selected_cell
* :ghissue:`1316`: Pass subprocess test runners a suitable location for xunit output
* :ghissue:`1315`: Collect results from subprocess runners and spit out Xunit XML output.
* :ghissue:`1233`: Update CodeMirror to the latest version
* :ghissue:`1234`: Refactor how the notebook focuses cells
* :ghissue:`1235`: After upgrading CodeMirror check the status of some bugs
* :ghissue:`1236`: Review how select is called when notebook cells are inserted
* :ghissue:`1303`: Updatecm
* :ghissue:`1311`: Fixing CM related indentation problems.
* :ghissue:`1304`: controller/server load can disrupt heartbeat
* :ghissue:`1312`: minor heartbeat tweaks
* :ghissue:`1302`: Input parsing with %prun clobbers escapes
* :ghissue:`1306`: Fix %prun input parsing for escaped characters (closes #1302)
* :ghissue:`1251`: IPython-0.12 can't import map module on Python 3.1
* :ghissue:`1202`: Pyreadline install exclusion for 64 bit windows no longer required,  version dependency not correctly specified.
* :ghissue:`1301`: New "Fix for issue #1202" based on current master.
* :ghissue:`1242`: changed key map name to match changes to python mode
* :ghissue:`1203`: Fix for issue #1202
* :ghissue:`1289`: Make autoreload extension work on Python 3.
* :ghissue:`1263`: Different 'C-x' for shortcut, 'C-m c' not toCodeCell anymore
* :ghissue:`1259`: Replace "from (.|..) import" with absolute imports.
* :ghissue:`1278`: took a crack at making notebook.html extend layout.html
* :ghissue:`1210`: Add 'quiet' option to suppress screen output during %prun calls, edited dochelp
* :ghissue:`1288`: Don't ask for confirmation when stdin isn't available
* :ghissue:`1290`: Cell-level cut & paste overwrites multiple cells
* :ghissue:`1291`: Minor, but important fixes to cut/copy/paste.
* :ghissue:`1293`: TaskScheduler.hwm default value
* :ghissue:`1294`: TaskScheduler.hwm default to 1 instead of 0
* :ghissue:`1281`: in Hub: registration_timeout must be an integer, but heartmonitor.period is CFloat
* :ghissue:`1283`: HeartMonitor.period should be an Integer
* :ghissue:`1162`: Allow merge/split adjacent cells in notebook
* :ghissue:`1264`: Aceify
* :ghissue:`1261`: Mergesplit
* :ghissue:`1269`: Another strange input handling error
* :ghissue:`1284`: a fix for GH 1269
* :ghissue:`1232`: Dead kernel loop
* :ghissue:`1279`: ImportError: cannot import name S1 (from logging)
* :ghissue:`1276`: notebook menu item to send a KeyboardInterrupt to the kernel
* :ghissue:`1213`: BUG: Minor typo in history_console_widget.py
* :ghissue:`1248`: IPython notebook doesn't work with lastest version of tornado
* :ghissue:`1267`: add NoDB for non-recording Hub
* :ghissue:`1222`: allow Reference as callable in map/apply
* :ghissue:`1257`: use self.kernel_manager_class in qtconsoleapp
* :ghissue:`1220`: Open a new notebook while connecting to an existing kernel (opened by qtconsole or terminal or standalone)
* :ghissue:`1253`: set auto_create flag for notebook apps
* :ghissue:`1260`: heartbeat failure on long gil-holding operation
* :ghissue:`1262`: Heartbeat no longer shares the app's Context
* :ghissue:`1225`: SyntaxError display broken in Python 3
* :ghissue:`1229`: Fix display of SyntaxError in Python 3
* :ghissue:`1256`: Dewijmoize
* :ghissue:`1246`: Skip tests that require X, when importing pylab results in RuntimeError.
* :ghissue:`1250`: Wijmoize
* :ghissue:`1244`: can not imput chinese word "造" , exit right now
* :ghissue:`1194`: Adding Opera 11 as a compatible browser for ipython notebook
* :ghissue:`1198`: Kernel Has Died error in Notebook
* :ghissue:`1211`: serve local files in notebook-dir
* :ghissue:`1224`: edit text cells on double-click instead of single-click
* :ghissue:`1187`: misc notebook: connection file cleanup, first heartbeat, startup flush
* :ghissue:`1207`: fix loadpy duplicating newlines
* :ghissue:`1060`: Always save the .py file to disk next to the .ipynb
* :ghissue:`1066`: execute cell in place should preserve the current insertion-point in the notebook
* :ghissue:`1141`: "In" numbers are not invalidated when restarting kernel
* :ghissue:`1231`: pip on OSX tries to install files in /System directory.
* :ghissue:`1129`: Unified setup.py
* :ghissue:`1199`: Reduce IPython.external.*
* :ghissue:`1219`: Make all the static files path absolute.
* :ghissue:`1218`: Added -q option to %prun for suppression of the output, along with editing the dochelp string.
* :ghissue:`1217`: Added -q option to %prun for suppression of the output, along with editing the dochelp string
* :ghissue:`1216`: Pdb tab completion does not work in QtConsole
* :ghissue:`1197`: Interactive shell trying to: from ... import history
* :ghissue:`1175`: core.completer: Clean up excessive and unused code.
* :ghissue:`1208`: should dv.sync_import print failed imports ?
* :ghissue:`1186`: payloadpage.py not used by qtconsole
* :ghissue:`1204`: double newline from %loadpy in python notebook (at least on mac)
* :ghissue:`1192`: Invalid JSON data
* :ghissue:`1196`: docs: looks like a file path might have been accidentally pasted in the middle of a word
* :ghissue:`1189`: Right justify of 'in' prompt in variable prompt size configurations
* :ghissue:`1185`: ipython console not work proper with stdout...
* :ghissue:`1191`: profile/startup files not executed with "notebook"
* :ghissue:`1190`: Fix link to Chris Fonnesbeck blog post about 0.11 highlights.
* :ghissue:`1174`: Remove %install_default_config and %install_profiles
=============
 0.11 Series
=============

Release 0.11
============

IPython 0.11 is a *major* overhaul of IPython, two years in the making.  Most
of the code base has been rewritten or at least reorganized, breaking backward
compatibility with several APIs in previous versions.  It is the first major
release in two years, and probably the most significant change to IPython since
its inception.  We plan to have a relatively quick succession of releases, as
people discover new bugs and regressions.  Once we iron out any significant
bugs in this process and settle down the new APIs, this series will become
IPython 1.0.  We encourage feedback now on the core APIs, which we hope to
maintain stable during the 1.0 series.

Since the internal APIs have changed so much, projects using IPython as a
library (as opposed to end-users of the application) are the most likely to
encounter regressions or changes that break their existing use patterns.  We
will make every effort to provide updated versions of the APIs to facilitate
the transition, and we encourage you to contact us on the `development mailing
list`__ with questions and feedback.

.. __: http://mail.scipy.org/mailman/listinfo/ipython-dev

Chris Fonnesbeck recently wrote an `excellent post`__ that highlights some of
our major new features, with examples and screenshots.  We encourage you to
read it as it provides an illustrated, high-level overview complementing the
detailed feature breakdown in this document.

.. __: http://stronginference.com/post/innovations-in-ipython

A quick summary of the major changes (see below for details):

* **Standalone Qt console**: a new rich console has been added to IPython,
  started with `ipython qtconsole`.  In this application we have tried to
  retain the feel of a terminal for fast and efficient workflows, while adding
  many features that a line-oriented terminal simply can not support, such as
  inline figures, full multiline editing with syntax highlighting, graphical
  tooltips for function calls and much more.  This development was sponsored by
  `Enthought Inc.`__. See :ref:`below <qtconsole_011>` for details.

.. __: http://enthought.com

* **High-level parallel computing with ZeroMQ**. Using the same architecture
  that our Qt console is based on, we have completely rewritten our high-level
  parallel computing machinery that in prior versions used the Twisted
  networking framework.  While this change will require users to update their
  codes, the improvements in performance, memory control and internal
  consistency across our codebase convinced us it was a price worth paying.  We
  have tried to explain how to best proceed with this update, and will be happy
  to answer questions that may arise.  A full tutorial describing these
  features `was presented at SciPy'11`__, more details :ref:`below
  <parallel_011>`.

.. __: http://minrk.github.com/scipy-tutorial-2011

* **New model for GUI/plotting support in the terminal**.  Now instead of the
  various `-Xthread` flags we had before, GUI support is provided without the
  use of any threads, by directly integrating GUI event loops with Python's
  `PyOS_InputHook` API.  A new command-line flag `--gui` controls GUI support,
  and it can also be enabled after IPython startup via the new `%gui` magic.
  This requires some changes if you want to execute GUI-using scripts inside
  IPython, see :ref:`the GUI support section <gui_support>` for more details.

* **A two-process architecture.** The Qt console is the first use of a new
  model that splits IPython between a kernel process where code is executed and
  a client that handles user interaction.  We plan on also providing terminal
  and web-browser based clients using this infrastructure in future releases.
  This model allows multiple clients to interact with an IPython process
  through a :ref:`well-documented messaging protocol <messaging>` using the
  ZeroMQ networking library.

* **Refactoring.** the entire codebase has been refactored, in order to make it
  more modular and easier to contribute to.  IPython has traditionally been a
  hard project to participate because the old codebase was very monolithic.  We
  hope this (ongoing) restructuring will make it easier for new developers to
  join us.

* **Vim integration**. Vim can be configured to seamlessly control an IPython
  kernel, see the files in :file:`docs/examples/vim` for the full details.
  This work was done by Paul Ivanov, who prepared a nice `video
  demonstration`__ of the features it provides.

.. __: http://pirsquared.org/blog/2011/07/28/vim-ipython/
  
* **Integration into Microsoft Visual Studio**. Thanks to the work of the
  Microsoft `Python Tools for Visual Studio`__ team, this version of IPython
  has been integrated into Microsoft Visual Studio's Python tools open source
  plug-in.  `Details below`_

.. __: http://pytools.codeplex.com
.. _details below: ms_visual_studio_011_

* **Improved unicode support**. We closed many bugs related to unicode input.

* **Python 3**. IPython now runs on Python 3.x. See :ref:`python3_011` for
  details.

* **New profile model**. Profiles are now directories that contain all relevant
  information for that session, and thus better isolate IPython use-cases.

* **SQLite storage for history**. All history is now stored in a SQLite
  database, providing support for multiple simultaneous sessions that won't
  clobber each other as well as the ability to perform queries on all stored
  data.

* **New configuration system**. All parts of IPython are now configured via a
  mechanism inspired by the Enthought Traits library.  Any configurable element
  can have its attributes set either via files that now use real Python syntax
  or from the command-line.

* **Pasting of code with prompts**. IPython now intelligently strips out input
  prompts , be they plain Python ones (``>>>`` and ``...``) or IPython ones
  (``In [N]:`` and ``...:``).  More details :ref:`here <pasting_with_prompts>`.
  

Authors and support
-------------------

Over 60 separate authors have contributed to this release, see :ref:`below
<credits_011>` for a full list.  In particular, we want to highlight the
extremely active participation of two new core team members: Evan Patterson
implemented the Qt console, and Thomas Kluyver started with our Python 3 port
and by now has made major contributions to just about every area of IPython.

We are also grateful for the support we have received during this development
cycle from several institutions:

- `Enthought Inc`__ funded the development of our new Qt console, an effort that
  required developing major pieces of underlying infrastructure, which now
  power not only the Qt console but also our new parallel machinery.  We'd like
  to thank Eric Jones and Travis Oliphant for their support, as well as Ilan
  Schnell for his tireless work integrating and testing IPython in the
  `Enthought Python Distribution`_.

.. __: http://enthought.com
.. _Enthought Python Distribution: http://www.enthought.com/products/epd.php

- Nipy/NIH: funding via the `NiPy project`__ (NIH grant 5R01MH081909-02) helped
  us jumpstart the development of this series by restructuring the entire
  codebase two years ago in a way that would make modular development and
  testing more approachable.  Without this initial groundwork, all the new
  features we have added would have been impossible to develop.

.. __: http://nipy.org

- Sage/NSF: funding via the grant `Sage: Unifying Mathematical Software for
  Scientists, Engineers, and Mathematicians`__ (NSF grant DMS-1015114)
  supported a meeting in spring 2011 of several of the core IPython developers
  where major progress was made integrating the last key pieces leading to this
  release.

.. __: http://modular.math.washington.edu/grants/compmath09

- Microsoft's team working on `Python Tools for Visual Studio`__ developed the 
  integratron of IPython into the Python plugin for Visual Studio 2010.

.. __: http://pytools.codeplex.com

- Google Summer of Code: in 2010, we had two students developing prototypes of
  the new machinery that is now maturing in this release: `Omar Zapata`_ and
  `Gerardo Gutiérrez`_.

.. _Omar Zapata: http://ipythonzmq.blogspot.com/2010/08/ipython-zmq-status.html
.. _Gerardo Gutiérrez: http://ipythonqt.blogspot.com/2010/04/ipython-qt-interface-gsoc-2010-proposal.html>


Development summary: moving to Git and Github
---------------------------------------------

In April 2010, after `one breakage too many with bzr`__, we decided to move our
entire development process to Git and Github.com.  This has proven to be one of
the best decisions in the project's history, as the combination of git and
github have made us far, far more productive than we could be with our previous
tools.  We first converted our bzr repo to a git one without losing history,
and a few weeks later ported all open Launchpad bugs to github issues with
their comments mostly intact (modulo some formatting changes).  This ensured a
smooth transition where no development history or submitted bugs were lost.
Feel free to use our little Launchpad to Github issues `porting script`_ if you
need to make a similar transition.

.. __: http://mail.scipy.org/pipermail/ipython-dev/2010-April/005944.html
.. _porting script: https://gist.github.com/835577

These simple statistics show how much work has been done on the new release, by
comparing the current code to the last point it had in common with the 0.10
series.  A huge diff and ~2200 commits make up this cycle::

    git diff $(git merge-base 0.10.2 HEAD)  | wc -l
    288019

    git log $(git merge-base 0.10.2 HEAD)..HEAD --oneline | wc -l
    2200

Since our move to github, 511 issues were closed, 226 of which were pull
requests and 285 regular issues (:ref:`a full list with links
<issues_list_011>` is available for those interested in the details).  Github's
pull requests are a fantastic mechanism for reviewing code and building a
shared ownership of the project, and we are making enthusiastic use of it.

.. Note::

   This undercounts the number of issues closed in this development cycle,
   since we only moved to github for issue tracking in May 2010, but we have no
   way of collecting statistics on the number of issues closed in the old
   Launchpad bug tracker prior to that.

   
.. _qtconsole_011:

Qt Console
----------

IPython now ships with a Qt application that feels very much like a terminal,
but is in fact a rich GUI that runs an IPython client but supports inline
figures, saving sessions to PDF and HTML, multiline editing with syntax
highlighting, graphical calltips and much more:

.. figure:: ../_images/qtconsole.png
    :width: 400px
    :alt: IPython Qt console with embedded plots
    :align: center
    :target: ../_images/qtconsole.png

    The Qt console for IPython, using inline matplotlib plots.

We hope that many projects will embed this widget, which we've kept
deliberately very lightweight, into their own environments.  In the future we
may also offer a slightly more featureful application (with menus and other GUI
elements), but we remain committed to always shipping this easy to embed
widget.

See the `Jupyter Qt Console site <https://jupyter.org/qtconsole>`_ for a detailed
description of the console's features and use.


.. _parallel_011:

High-level parallel computing with ZeroMQ
-----------------------------------------

We have completely rewritten the Twisted-based code for high-level parallel
computing to work atop our new ZeroMQ architecture.  While we realize this will
break compatibility for a number of users, we hope to make the transition as
easy as possible with our docs, and we are convinced the change is worth it.
ZeroMQ provides us with much tighter control over memory, higher performance,
and its communications are impervious to the Python Global Interpreter Lock
because they take place in a system-level C++ thread.  The impact of the GIL in
our previous code was something we could simply not work around, given that
Twisted is itself a Python library.  So while Twisted is a very capable
framework, we think ZeroMQ fits our needs much better and we hope you will find
the change to be a significant improvement in the long run.

Our manual contains a full description of how to use IPython for parallel
computing, and the `tutorial`__ presented by Min
Ragan-Kelley at the SciPy 2011 conference provides a hands-on complement to the
reference docs.

.. __: http://minrk.github.com/scipy-tutorial-2011


Refactoring
-----------

As of this release, a significant portion of IPython has been refactored.  This
refactoring is founded on a number of new abstractions.  The main new classes
that implement these abstractions are:

* :class:`traitlets.HasTraits`.
* :class:`traitlets.config.configurable.Configurable`.
* :class:`traitlets.config.application.Application`.
* :class:`traitlets.config.loader.ConfigLoader`.
* :class:`traitlets.config.loader.Config`

We are still in the process of writing developer focused documentation about
these classes, but for now our :ref:`configuration documentation
<config_overview>` contains a high level overview of the concepts that these
classes express.

The biggest user-visible change is likely the move to using the config system
to determine the command-line arguments for IPython applications. The benefit
of this is that *all* configurable values in IPython are exposed on the
command-line, but the syntax for specifying values has changed. The gist is
that assigning values is pure Python assignment.  Simple flags exist for
commonly used options, these are always prefixed with '--'.

The IPython command-line help has the details of all the options (via
``ipython --help``), but a simple example should clarify things; the ``pylab``
flag can be used to start in pylab mode with the qt4 backend::

  ipython --pylab=qt

which is equivalent to using the fully qualified form::

  ipython --TerminalIPythonApp.pylab=qt

The long-form options can be listed via ``ipython --help-all``.


ZeroMQ architecture
-------------------

There is a new GUI framework for IPython, based on a client-server model in
which multiple clients can communicate with one IPython kernel, using the
ZeroMQ messaging framework. There is already a Qt console client, which can
be started by calling ``ipython qtconsole``. The protocol is :ref:`documented
<messaging>`.

The parallel computing framework has also been rewritten using ZMQ. The
protocol is described :ref:`here <ipyparallel:/reference/messages.md>`, and the code is in the
new :mod:`IPython.parallel` module.

.. _python3_011:

Python 3 support
----------------

A Python 3 version of IPython has been prepared. For the time being, this is
maintained separately and updated from the main codebase. Its code can be found
`here <https://github.com/ipython/ipython-py3k>`_. The parallel computing
components are not perfect on Python3, but most functionality appears to be
working.  As this work is evolving quickly, the best place to find updated
information about it is our `Python 3 wiki page`__.

.. __: http://wiki.ipython.org/index.php?title=Python_3


Unicode
-------

Entering non-ascii characters in unicode literals (``u"€ø"``) now works
properly on all platforms. However, entering these in byte/string literals
(``"€ø"``) will not work as expected on Windows (or any platform where the
terminal encoding is not UTF-8, as it typically is for Linux & Mac OS X). You
can use escape sequences (``"\xe9\x82"``) to get bytes above 128, or use
unicode literals and encode them. This is a limitation of Python 2 which we
cannot easily work around.

.. _ms_visual_studio_011:

Integration with Microsoft Visual Studio
----------------------------------------

IPython can be used as the interactive shell in the `Python plugin for
Microsoft Visual Studio`__, as seen here:

.. figure:: ../_images/ms_visual_studio.png
    :width: 500px
    :alt: IPython console embedded in Microsoft Visual Studio.
    :align: center
    :target: ../_images/ms_visual_studio.png

    IPython console embedded in Microsoft Visual Studio.

The Microsoft team developing this currently has a release candidate out using
IPython 0.11. We will continue to collaborate with them to ensure that as they
approach their final release date, the integration with IPython remains smooth.
We'd like to thank Dino Viehland and Shahrokh Mortazavi for the work they have
done towards this feature, as well as Wenming Ye for his support of our WinHPC
capabilities.

.. __: http://pytools.codeplex.com


Additional new features
-----------------------

* Added ``Bytes`` traitlet, removing ``Str``.  All 'string' traitlets should
  either be ``Unicode`` if a real string, or ``Bytes`` if a C-string. This
  removes ambiguity and helps the Python 3 transition.

* New magic ``%loadpy`` loads a python file from disk or web URL into
  the current input buffer.

* New magic ``%pastebin`` for sharing code via the 'Lodge it' pastebin.

* New magic ``%precision`` for controlling float and numpy pretty printing.

* IPython applications initiate logging, so any object can gain access to
  a the logger of the currently running Application with:

.. sourcecode:: python

    from traitlets.config.application import Application
    logger = Application.instance().log

* You can now get help on an object halfway through typing a command. For
  instance, typing ``a = zip?`` shows the details of :func:`zip`. It also
  leaves the command at the next prompt so you can carry on with it.

* The input history is now written to an SQLite database. The API for
  retrieving items from the history has also been redesigned.

* The :mod:`IPython.extensions.pretty` extension has been moved out of
  quarantine and fully updated to the new extension API.

* New magics for loading/unloading/reloading extensions have been added:
  ``%load_ext``, ``%unload_ext`` and ``%reload_ext``.

* The configuration system and configuration files are brand new. See the
  configuration system :ref:`documentation <config_index>` for more details.

* The :class:`~IPython.core.interactiveshell.InteractiveShell` class is now a
  :class:`~traitlets.config.configurable.Configurable` subclass and has traitlets
  that determine the defaults and runtime environment. The ``__init__`` method
  has also been refactored so this class can be instantiated and run without
  the old :mod:`ipmaker` module.

* The methods of :class:`~IPython.core.interactiveshell.InteractiveShell` have
  been organized into sections to make it easier to turn more sections
  of functionality into components.

* The embedded shell has been refactored into a truly standalone subclass of
  :class:`InteractiveShell` called :class:`InteractiveShellEmbed`.  All
  embedding logic has been taken out of the base class and put into the 
  embedded subclass.

* Added methods of :class:`~IPython.core.interactiveshell.InteractiveShell` to
  help it cleanup after itself. The :meth:`cleanup` method controls this. We
  couldn't do this in :meth:`__del__` because we have cycles in our object
  graph that prevent it from being called.

* Created a new module :mod:`IPython.utils.importstring` for resolving
  strings like ``foo.bar.Bar`` to the actual class.

* Completely refactored the :mod:`IPython.core.prefilter` module into
  :class:`~traitlets.config.configurable.Configurable` subclasses. Added a new
  layer into the prefilter system, called "transformations" that all new
  prefilter logic should use (rather than the older "checker/handler"
  approach).

* Aliases are now components (:mod:`IPython.core.alias`).

* New top level :func:`~IPython.frontend.terminal.embed.embed` function that can
  be called to embed IPython at any place in user's code. On the first call it
  will create an :class:`~IPython.frontend.terminal.embed.InteractiveShellEmbed`
  instance and call it. In later calls, it just calls the previously created
  :class:`~IPython.frontend.terminal.embed.InteractiveShellEmbed`.

* Created a configuration system (:mod:`traitlets.config.configurable`) that is
  based on :mod:`traitlets`. Configurables are arranged into a
  runtime containment tree (not inheritance) that i) automatically propagates
  configuration information and ii) allows singletons to discover each other in
  a loosely coupled manner. In the future all parts of IPython will be
  subclasses of :class:`~traitlets.config.configurable.Configurable`. All IPython
  developers should become familiar with the config system.

* Created a new :class:`~traitlets.config.loader.Config` for holding
  configuration information. This is a dict like class with a few extras: i)
  it supports attribute style access, ii) it has a merge function that merges
  two :class:`~traitlets.config.loader.Config` instances recursively and iii) it
  will automatically create sub-:class:`~traitlets.config.loader.Config`
  instances for attributes that start with an uppercase character.

* Created new configuration loaders in :mod:`traitlets.config.loader`. These
  loaders provide a unified loading interface for all configuration
  information including command line arguments and configuration files. We
  have two default implementations based on :mod:`argparse` and plain python
  files.  These are used to implement the new configuration system.

* Created a top-level :class:`Application` class in
  :mod:`IPython.core.application` that is designed to encapsulate the starting
  of any basic Python program. An application loads and merges all the
  configuration objects, constructs the main application, configures and
  initiates logging, and creates and configures any :class:`Configurable`
  instances and then starts the application running. An extended
  :class:`BaseIPythonApplication` class adds logic for handling the
  IPython directory as well as profiles, and all IPython entry points
  extend it.

* The :class:`Type` and :class:`Instance` traitlets now handle classes given
  as strings, like ``foo.bar.Bar``. This is needed for forward declarations.
  But, this was implemented in a careful way so that string to class
  resolution is done at a single point, when the parent
  :class:`~traitlets.HasTraitlets` is instantiated.

* :mod:`IPython.utils.ipstruct` has been refactored to be a subclass of 
  dict.  It also now has full docstrings and doctests.

* Created a Traits like implementation in :mod:`traitlets`.  This
  is a pure Python, lightweight version of a library that is similar to
  Enthought's Traits project, but has no dependencies on Enthought's code.  We
  are using this for validation, defaults and notification in our new component
  system.  Although it is not 100% API compatible with Enthought's Traits, we
  plan on moving in this direction so that eventually our implementation could
  be replaced by a (yet to exist) pure Python version of Enthought Traits.

* Added a new module :mod:`IPython.lib.inputhook` to manage the integration
  with GUI event loops using `PyOS_InputHook`.  See the docstrings in this
  module or the main IPython docs for details.

* For users, GUI event loop integration is now handled through the new
  :command:`%gui` magic command.  Type ``%gui?`` at an IPython prompt for
  documentation.

* For developers :mod:`IPython.lib.inputhook` provides a simple interface
  for managing the event loops in their interactive GUI applications.
  Examples can be found in our :file:`examples/lib` directory.

Backwards incompatible changes
------------------------------

* The Twisted-based :mod:`IPython.kernel` has been removed, and completely
  rewritten as :mod:`IPython.parallel`, using ZeroMQ.

* Profiles are now directories. Instead of a profile being a single config file,
  profiles are now self-contained directories. By default, profiles get their
  own IPython history, log files, and everything. To create a new profile, do
  ``ipython profile create <name>``.

* All IPython applications have been rewritten to use
  :class:`~traitlets.config.loader.KeyValueConfigLoader`. This means that
  command-line options have changed. Now, all configurable values are accessible
  from the command-line with the same syntax as in a configuration file.

* The command line options ``-wthread``, ``-qthread`` and
  ``-gthread`` have been removed. Use ``--gui=wx``, ``--gui=qt``, ``--gui=gtk``
  instead.

* The extension loading functions have been renamed to
  :func:`load_ipython_extension` and :func:`unload_ipython_extension`.

* :class:`~IPython.core.interactiveshell.InteractiveShell` no longer takes an
  ``embedded`` argument. Instead just use the
  :class:`~IPython.core.interactiveshell.InteractiveShellEmbed` class.

* ``__IPYTHON__`` is no longer injected into ``__builtin__``.

* :meth:`Struct.__init__` no longer takes `None` as its first argument.  It
  must be a :class:`dict` or :class:`Struct`.

* :meth:`~IPython.core.interactiveshell.InteractiveShell.ipmagic` has been
  renamed :meth:`~IPython.core.interactiveshell.InteractiveShell.magic.`

* The functions :func:`ipmagic` and :func:`ipalias` have been removed from
  :mod:`__builtins__`.

* The references to the global
  :class:`~IPython.core.interactivehell.InteractiveShell` instance (``_ip``, and
  ``__IP``) have been removed from the user's namespace. They are replaced by a
  new function called :func:`get_ipython` that returns the current
  :class:`~IPython.core.interactiveshell.InteractiveShell` instance. This
  function is injected into the user's namespace and is now the main way of
  accessing the running IPython.

* Old style configuration files :file:`ipythonrc` and :file:`ipy_user_conf.py`
  are no longer supported. Users should migrate there configuration files to
  the new format described :doc:`here </config/intro>` and
  :ref:`here <config_overview>`.

* The old IPython extension API that relied on :func:`ipapi` has been
  completely removed. The new extension API is described :ref:`here
  <extensions_overview>`.

* Support for ``qt3`` has been dropped.  Users who need this should use
  previous versions of IPython.

* Removed :mod:`shellglobals` as it was obsolete.

* Removed all the threaded shells in :mod:`IPython.core.shell`. These are no
  longer needed because of the new capabilities in
  :mod:`IPython.lib.inputhook`.

* New top-level sub-packages have been created: :mod:`IPython.core`, 
  :mod:`IPython.lib`, :mod:`IPython.utils`, :mod:`IPython.deathrow`,
  :mod:`IPython.quarantine`.  All existing top-level modules have been
  moved to appropriate sub-packages.  All internal import statements
  have been updated and tests have been added.  The build system (setup.py
  and friends) have been updated. See :doc:`/api/index` for details of these
  new sub-packages.

* :mod:`IPython.ipapi` has been moved to :mod:`IPython.core.ipapi`.
  :mod:`IPython.Shell` and :mod:`IPython.iplib` have been split and removed as
  part of the refactor.

* :mod:`Extensions` has been moved to :mod:`extensions` and all existing
  extensions have been moved to either :mod:`IPython.quarantine` or
  :mod:`IPython.deathrow`. :mod:`IPython.quarantine` contains modules that we
  plan on keeping but that need to be updated. :mod:`IPython.deathrow` contains
  modules that are either dead or that should be maintained as third party
  libraries.

* Previous IPython GUIs in :mod:`IPython.frontend` and :mod:`IPython.gui` are
  likely broken, and have been removed to :mod:`IPython.deathrow` because of the
  refactoring in the core. With proper updates, these should still work.


Known Regressions
-----------------

We do our best to improve IPython, but there are some known regressions in 0.11
relative to 0.10.2.  First of all, there are features that have yet to be
ported to the new APIs, and in order to ensure that all of the installed code
runs for our users, we have moved them to two separate directories in the
source distribution, `quarantine` and `deathrow`.  Finally, we have some other
miscellaneous regressions that we hope to fix as soon as possible.  We now
describe all of these in more detail.

Quarantine
~~~~~~~~~~

These are tools and extensions that we consider relatively easy to update to
the new classes and APIs, but that we simply haven't had time for.  Any user
who is interested in one of these is encouraged to help us by porting it and
submitting a pull request on our `development site`_.

.. _development site: http://github.com/ipython/ipython

Currently, the quarantine directory contains::

    clearcmd.py            ipy_fsops.py            ipy_signals.py
    envpersist.py          ipy_gnuglobal.py        ipy_synchronize_with.py
    ext_rescapture.py      ipy_greedycompleter.py  ipy_system_conf.py
    InterpreterExec.py     ipy_jot.py              ipy_which.py
    ipy_app_completers.py  ipy_lookfor.py          ipy_winpdb.py
    ipy_autoreload.py      ipy_profile_doctest.py  ipy_workdir.py
    ipy_completers.py      ipy_pydb.py             jobctrl.py
    ipy_editors.py         ipy_rehashdir.py        ledit.py
    ipy_exportdb.py        ipy_render.py           pspersistence.py
    ipy_extutil.py         ipy_server.py           win32clip.py

Deathrow
~~~~~~~~

These packages may be harder to update or make most sense as third-party
libraries.  Some of them are completely obsolete and have been already replaced
by better functionality (we simply haven't had the time to carefully weed them
out so they are kept here for now). Others simply require fixes to code that
the current core team may not be familiar with.  If a tool you were used to is
included here, we encourage you to contact the dev list and we can discuss
whether it makes sense to keep it in IPython (if it can be maintained).

Currently, the deathrow directory contains::

    astyle.py              ipy_defaults.py          ipy_vimserver.py
    dtutils.py             ipy_kitcfg.py            numeric_formats.py
    Gnuplot2.py            ipy_legacy.py            numutils.py
    GnuplotInteractive.py  ipy_p4.py                outputtrap.py
    GnuplotRuntime.py      ipy_profile_none.py      PhysicalQInput.py
    ibrowse.py             ipy_profile_numpy.py     PhysicalQInteractive.py
    igrid.py               ipy_profile_scipy.py     quitter.py*
    ipipe.py               ipy_profile_sh.py        scitedirector.py
    iplib.py               ipy_profile_zope.py      Shell.py
    ipy_constants.py       ipy_traits_completer.py  twshell.py


Other regressions
~~~~~~~~~~~~~~~~~

* The machinery that adds functionality to the 'sh' profile for using IPython
  as your system shell has not been updated to use the new APIs.  As a result,
  only the aesthetic (prompt) changes are still implemented. We intend to fix
  this by 0.12.  Tracked as issue 547_.

.. _547: https://github.com/ipython/ipython/issues/547

* The installation of scripts on Windows was broken without setuptools, so we
  now depend on setuptools on Windows.  We hope to fix setuptools-less
  installation, and then remove the setuptools dependency.  Issue 539_.

.. _539: https://github.com/ipython/ipython/issues/539

* The directory history `_dh` is not saved between sessions.  Issue 634_.

.. _634: https://github.com/ipython/ipython/issues/634


Removed Features
----------------

As part of the updating of IPython, we have removed a few features for the
purposes of cleaning up the codebase and interfaces.  These removals are
permanent, but for any item listed below, equivalent functionality is
available.

* The magics Exit and Quit have been dropped as ways to exit IPython. Instead,
  the lowercase forms of both work either as a bare name (``exit``) or a
  function call (``exit()``).  You can assign these to other names using
  exec_lines in the config file.


.. _credits_011:

Credits
-------

Many users and developers contributed code, features, bug reports and ideas to
this release.  Please do not hesitate in contacting us if we've failed to
acknowledge your contribution here.  In particular, for this release we have
contribution from the following people, a mix of new and regular names (in
alphabetical order by first name):

* Aenugu Sai Kiran Reddy <saikrn08-at-gmail.com>
* andy wilson <wilson.andrew.j+github-at-gmail.com>
* Antonio Cuni <antocuni>
* Barry Wark <barrywark-at-gmail.com>
* Beetoju Anuradha <anu.beethoju-at-gmail.com>
* Benjamin Ragan-Kelley <minrk-at-Mercury.local>
* Brad Reisfeld
* Brian E. Granger <ellisonbg-at-gmail.com>
* Christoph Gohlke <cgohlke-at-uci.edu>
* Cody Precord
* dan.milstein
* Darren Dale <dsdale24-at-gmail.com>
* Dav Clark <davclark-at-berkeley.edu>
* David Warde-Farley <wardefar-at-iro.umontreal.ca>
* epatters <ejpatters-at-gmail.com>
* epatters <epatters-at-caltech.edu>
* epatters <epatters-at-enthought.com>
* Eric Firing <efiring-at-hawaii.edu>
* Erik Tollerud <erik.tollerud-at-gmail.com>
* Evan Patterson <epatters-at-enthought.com>
* Fernando Perez <Fernando.Perez-at-berkeley.edu>
* Gael Varoquaux <gael.varoquaux-at-normalesup.org>
* Gerardo <muzgash-at-Muzpelheim>
* Jason Grout <jason.grout-at-drake.edu>
* John Hunter <jdh2358-at-gmail.com>
* Jens Hedegaard Nielsen <jenshnielsen-at-gmail.com>
* Johann Cohen-Tanugi <johann.cohentanugi-at-gmail.com>
* Jörgen Stenarson <jorgen.stenarson-at-bostream.nu>
* Justin Riley <justin.t.riley-at-gmail.com>
* Kiorky
* Laurent Dufrechou <laurent.dufrechou-at-gmail.com>
* Luis Pedro Coelho <lpc-at-cmu.edu>
* Mani chandra <mchandra-at-iitk.ac.in>
* Mark E. Smith
* Mark Voorhies <mark.voorhies-at-ucsf.edu>
* Martin Spacek <git-at-mspacek.mm.st>
* Michael Droettboom <mdroe-at-stsci.edu>
* MinRK <benjaminrk-at-gmail.com>
* muzuiget <muzuiget-at-gmail.com>
* Nick Tarleton <nick-at-quixey.com>
* Nicolas Rougier <Nicolas.rougier-at-inria.fr>
* Omar Andres Zapata Mesa <andresete.chaos-at-gmail.com>
* Paul Ivanov <pivanov314-at-gmail.com>
* Pauli Virtanen <pauli.virtanen-at-iki.fi>
* Prabhu Ramachandran
* Ramana <sramana9-at-gmail.com>
* Robert Kern <robert.kern-at-gmail.com>
* Sathesh Chandra <satheshchandra88-at-gmail.com>
* Satrajit Ghosh <satra-at-mit.edu>
* Sebastian Busch
* Skipper Seabold <jsseabold-at-gmail.com>
* Stefan van der Walt <bzr-at-mentat.za.net>
* Stephan Peijnik <debian-at-sp.or.at>
* Steven Bethard
* Thomas Kluyver <takowl-at-gmail.com>
* Thomas Spura <tomspur-at-fedoraproject.org>
* Tom Fetherston <tfetherston-at-aol.com>
* Tom MacWright
* tzanko
* vankayala sowjanya <hai.sowjanya-at-gmail.com>
* Vivian De Smedt <vds2212-at-VIVIAN>
* Ville M. Vainio <vivainio-at-gmail.com>
* Vishal Vatsa <vishal.vatsa-at-gmail.com>
* Vishnu S G <sgvishnu777-at-gmail.com>
* Walter Doerwald <walter-at-livinglogic.de>

.. note::

    This list was generated with the output of
    ``git log dev-0.11 HEAD --format='* %aN <%aE>' | sed 's/@/\-at\-/' | sed 's/<>//' | sort -u``
    after some cleanup.  If you should be on this list, please add yourself.
============
 8.x Series
============


IPython 8.0.1 (CVE-2022-21699)
------------------------------

IPython 8.0.1, 7.31.1 and 5.11 are security releases that change some default
values in order to prevent potential Execution with Unnecessary Privileges.

Almost all version of IPython looks for configuration and profiles in current
working directory. Since IPython was developed before pip and environments
existed it was used a convenient way to load code/packages in a project
dependant way.

In 2022, it is not necessary anymore, and can lead to confusing behavior where
for example cloning a repository and starting IPython or loading a notebook from
any Jupyter-Compatible interface that has ipython set as a kernel can lead to
code execution.


I did not find any standard way for packaged to advertise CVEs they fix, I'm
thus trying to add a ``__patched_cves__`` attribute to the IPython module that
list the CVEs that should have been fixed. This attribute is informational only
as if a executable has a flaw, this value can always be changed by an attacker.

.. code::

    In [1]: import IPython

    In [2]: IPython.__patched_cves__
    Out[2]: {'CVE-2022-21699'}

    In [3]: 'CVE-2022-21699' in IPython.__patched_cves__
    Out[3]: True

Thus starting with this version:

 - The current working directory is not searched anymore for profiles or
   configurations files.
 - Added a ``__patched_cves__`` attribute (set of strings) to IPython module that contain
   the list of fixed CVE. This is informational only.

Further details can be read on the `GitHub Advisory <https://github.com/ipython/ipython/security/advisories/GHSA-pq7m-3gw7-gq5x>`__



IPython 8.0
-----------

IPython 8.0 is bringing a large number of new features and improvements to both the
user of the terminal and of the kernel via Jupyter. The removal of compatibility
with older version of Python is also the opportunity to do a couple of
performance improvements in particular with respect to startup time.
The 8.x branch started diverging from its predecessor around IPython 7.12
(January 2020).

This release contains 250+ pull requests, in addition to many of the features
and backports that have made it to the 7.x branch. Please see the 
`8.0 milestone <https://github.com/ipython/ipython/milestone/73?closed=1>`__ for the full list of pull requests.

Please feel free to send pull requests to updates those notes after release, 
I have likely forgotten a few things reviewing 250+ PRs.

Dependencies changes/downstream packaging
-----------------------------------------

Most of our building steps have been changed to be (mostly) declarative
and follow PEP 517. We are trying to completely remove ``setup.py`` (:ghpull:`13238`) and are
looking for help to do so.

 - minimum supported ``traitlets`` version is now 5+
 - we now require ``stack_data``
 - minimal Python is now 3.8
 - ``nose`` is not a testing requirement anymore
 - ``pytest`` replaces nose.
 - ``iptest``/``iptest3`` cli entrypoints do not exists anymore. 
 - minimum officially support ``numpy`` version has been bumped, but this should
   not have much effect on packaging.


Deprecation and removal
-----------------------

We removed almost all features, arguments, functions, and modules that were
marked as deprecated between IPython 1.0 and 5.0. As a reminder, 5.0 was released
in 2016, and 1.0 in 2013. Last release of the 5 branch was 5.10.0, in May 2020.
The few remaining deprecated features we left have better deprecation warnings
or have been turned into explicit errors for better error messages.

I will use this occasion to add the following requests to anyone emitting a
deprecation warning:

 - Please add at least ``stacklevel=2`` so that the warning is emitted into the
   caller context, and not the callee one.
 - Please add **since which version** something is deprecated.

As a side note, it is much easier to conditionally compare version
numbers rather than using ``try/except`` when functionality changes with a version. 

I won't list all the removed features here, but modules like ``IPython.kernel``,
which was just a shim module around ``ipykernel`` for the past 8 years, have been
removed, and so many other similar things that pre-date the name **Jupyter**
itself.

We no longer need to add ``IPython.extensions`` to the PYTHONPATH because that is being
handled by ``load_extension``.

We are also removing ``Cythonmagic``, ``sympyprinting`` and ``rmagic`` as they are now in
other packages and no longer need to be inside IPython.


Documentation
-------------

The majority of our docstrings have now been reformatted and automatically fixed by
the experimental `Vélin <https://pypi.org/project/velin/>`_ project to conform
to numpydoc.

Type annotations
----------------

While IPython itself is highly dynamic and can't be completely typed, many of
the functions now have type annotations, and part of the codebase is now checked
by mypy.


Featured changes
----------------

Here is a features list of changes in IPython 8.0. This is of course non-exhaustive. 
Please note as well that many features have been added in the 7.x branch as well
(and hence why you want to read the 7.x what's new notes), in particular
features contributed by QuantStack (with respect to debugger protocol and Xeus
Python), as well as many debugger features that I was pleased to implement as
part of my work at QuanSight and sponsored by DE Shaw.

Traceback improvements
~~~~~~~~~~~~~~~~~~~~~~

Previously, error tracebacks for errors happening in code cells were showing a
hash, the one used for compiling the Python AST::

    In [1]: def foo():
    ...:     return 3 / 0
    ...:

    In [2]: foo()
    ---------------------------------------------------------------------------
    ZeroDivisionError                         Traceback (most recent call last)
    <ipython-input-2-c19b6d9633cf> in <module>
    ----> 1 foo()

    <ipython-input-1-1595a74c32d5> in foo()
        1 def foo():
    ----> 2     return 3 / 0
        3

    ZeroDivisionError: division by zero

The error traceback is now correctly formatted, showing the cell number in which the error happened::

    In [1]: def foo():
    ...:     return 3 / 0
    ...:

    Input In [2]: foo()
    ---------------------------------------------------------------------------
    ZeroDivisionError                         Traceback (most recent call last)
    input In [2], in <module>
    ----> 1 foo()

    Input In [1], in foo()
        1 def foo():
    ----> 2     return 3 / 0

    ZeroDivisionError: division by zero

The ``stack_data`` package has been integrated, which provides smarter information in the traceback; 
in particular it will highlight the AST node where an error occurs which can help to quickly narrow down errors.

For example in the following snippet::

    def foo(i):
        x = [[[0]]]
        return x[0][i][0]


    def bar():
        return foo(0) + foo(
            1
        ) + foo(2)


calling ``bar()`` would raise an ``IndexError`` on the return line of ``foo``,
and IPython 8.0 is capable of telling you where the index error occurs::


    IndexError
    Input In [2], in <module>
    ----> 1 bar()
            ^^^^^

    Input In [1], in bar()
          6 def bar():
    ----> 7     return foo(0) + foo(
                                ^^^^
          8         1
             ^^^^^^^^
          9     ) + foo(2)
             ^^^^

    Input In [1], in foo(i)
          1 def foo(i):
          2     x = [[[0]]]
    ----> 3     return x[0][i][0]
                       ^^^^^^^

The corresponding locations marked here with ``^`` will show up highlighted in 
the terminal and notebooks.

Finally, a colon ``::`` and line number is appended after a filename in
traceback::


    ZeroDivisionError               Traceback (most recent call last)
    File ~/error.py:4, in <module>
          1 def f():
          2     1/0
    ----> 4 f()

    File ~/error.py:2, in f()
          1 def f():
    ----> 2     1/0

Many terminals and editors have integrations enabling you to directly jump to the
relevant file/line when this syntax is used, so this small addition may have a high
impact on productivity.


Autosuggestions
~~~~~~~~~~~~~~~

Autosuggestion is a very useful feature available in `fish <https://fishshell.com/>`__, `zsh <https://en.wikipedia.org/wiki/Z_shell>`__, and `prompt-toolkit <https://python-prompt-toolkit.readthedocs.io/en/master/pages/asking_for_input.html#auto-suggestion>`__.

`Ptpython <https://github.com/prompt-toolkit/ptpython#ptpython>`__ allows users to enable this feature in
`ptpython/config.py <https://github.com/prompt-toolkit/ptpython/blob/master/examples/ptpython_config/config.py#L90>`__.

This feature allows users to accept autosuggestions with ctrl e, ctrl f,
or right arrow as described below.

1. Start ipython

.. image:: ../_images/8.0/auto_suggest_1_prompt_no_text.png

2. Run ``print("hello")``

.. image:: ../_images/8.0/auto_suggest_2_print_hello_suggest.png

3. start typing ``print`` again to see the autosuggestion

.. image:: ../_images/8.0/auto_suggest_3_print_hello_suggest.png

4. Press ``ctrl-f``, or ``ctrl-e``, or ``right-arrow`` to accept the suggestion

.. image:: ../_images/8.0/auto_suggest_4_print_hello.png

You can also complete word by word:

1. Run ``def say_hello(): print("hello")``

.. image:: ../_images/8.0/auto_suggest_second_prompt.png

2. Start typing  the first letter if ``def`` to see the autosuggestion

.. image:: ../_images/8.0/auto_suggest_d_phantom.png

3. Press ``alt-f`` (or ``escape`` followed by ``f``), to accept the first word of the suggestion

.. image:: ../_images/8.0/auto_suggest_def_phantom.png

Importantly, this feature does not interfere with tab completion:

1. After running ``def say_hello(): print("hello")``, press d

.. image:: ../_images/8.0/auto_suggest_d_phantom.png

2. Press Tab to start tab completion

.. image:: ../_images/8.0/auto_suggest_d_completions.png

3A. Press Tab again to select the first option

.. image:: ../_images/8.0/auto_suggest_def_completions.png

3B. Press ``alt f`` (``escape``, ``f``) to accept to accept the first word of the suggestion

.. image:: ../_images/8.0/auto_suggest_def_phantom.png

3C. Press ``ctrl-f`` or ``ctrl-e`` to accept the entire suggestion

.. image:: ../_images/8.0/auto_suggest_match_parens.png


Currently, autosuggestions are only shown in the emacs or vi insert editing modes:

- The ctrl e, ctrl f, and alt f shortcuts work by default in emacs mode.
- To use these shortcuts in vi insert mode, you will have to create `custom keybindings in your config.py <https://github.com/mskar/setup/commit/2892fcee46f9f80ef7788f0749edc99daccc52f4/>`__.


Show pinfo information in ipdb using "?" and "??"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In IPDB, it is now possible to show the information about an object using "?"
and "??", in much the same way that it can be done when using the IPython prompt::

    ipdb> partial?
    Init signature: partial(self, /, *args, **kwargs)
    Docstring:
    partial(func, *args, **keywords) - new function with partial application
    of the given arguments and keywords.
    File:           ~/.pyenv/versions/3.8.6/lib/python3.8/functools.py
    Type:           type
    Subclasses:

Previously, ``pinfo`` or ``pinfo2`` command had to be used for this purpose.


Autoreload 3 feature
~~~~~~~~~~~~~~~~~~~~

Example: When an IPython session is run with the 'autoreload' extension loaded,
you will now have the option '3' to select, which means the following:

    1. replicate all functionality from option 2
    2. autoload all new funcs/classes/enums/globals from the module when they are added
    3. autoload all newly imported funcs/classes/enums/globals from external modules

Try ``%autoreload 3`` in an IPython session after running ``%load_ext autoreload``.

For more information please see the following unit test : ``extensions/tests/test_autoreload.py:test_autoload_newly_added_objects``

Auto formatting with black in the CLI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If ``black`` is installed in the same environment as IPython, terminal IPython
will now *by default*  reformat the code in the CLI when possible. You can
disable this with ``--TerminalInteractiveShell.autoformatter=None``.

This feature was present in 7.x, but disabled by default.


History Range Glob feature
~~~~~~~~~~~~~~~~~~~~~~~~~~

Previously, when using ``%history``, users could specify either
a range of sessions and lines, for example:

.. code-block:: python

   ~8/1-~6/5   # see history from the first line of 8 sessions ago,
               # to the fifth line of 6 sessions ago.``

Or users could specify a glob pattern:

.. code-block:: python

   -g <pattern>  # glob ALL history for the specified pattern.

However users could *not* specify both.

If a user *did* specify both a range and a glob pattern,
then the glob pattern would be used (globbing *all* history) *and the range would be ignored*.

With this enhancement, if a user specifies both a range and a glob pattern, then the glob pattern will be applied to the specified range of history.

Don't start a multi-line cell with sunken parenthesis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From now on, IPython will not ask for the next line of input when given a single
line with more closing than opening brackets. For example, this means that if
you (mis)type ``]]`` instead of ``[]``, a ``SyntaxError`` will show up, instead of
the ``...:`` prompt continuation.

IPython shell for ipdb interact
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ipdb ``interact`` starts an IPython shell instead of Python's built-in ``code.interact()``.

Automatic Vi prompt stripping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When pasting code into IPython, it will strip the leading prompt characters if
there are any. For example, you can paste the following code into the console -
it will still work, even though each line is prefixed with prompts (`In`,
`Out`)::

    In [1]: 2 * 2 == 4
    Out[1]: True

    In [2]: print("This still works as pasted")


Previously, this was not the case for the Vi-mode prompts::

    In [1]: [ins] In [13]: 2 * 2 == 4
       ...: Out[13]: True
       ...:
      File "<ipython-input-1-727bb88eaf33>", line 1
        [ins] In [13]: 2 * 2 == 4
              ^
    SyntaxError: invalid syntax

This is now fixed, and Vi prompt prefixes - ``[ins]`` and ``[nav]`` -  are
skipped just as the normal ``In`` would be.

IPython shell can be started in the Vi mode using ``ipython --TerminalInteractiveShell.editing_mode=vi``, 
You should be able to change mode dynamically with ``%config TerminalInteractiveShell.editing_mode='vi'``

Empty History Ranges
~~~~~~~~~~~~~~~~~~~~

A number of magics that take history ranges can now be used with an empty
range. These magics are:

 * ``%save``
 * ``%load``
 * ``%pastebin``
 * ``%pycat``

Using them this way will make them take the history of the current session up
to the point of the magic call (such that the magic itself will not be
included).

Therefore it is now possible to save the whole history to a file using
``%save <filename>``, load and edit it using ``%load`` (makes for a nice usage
when followed with :kbd:`F2`), send it to `dpaste.org <http://dpast.org>`_ using
``%pastebin``, or view the whole thing syntax-highlighted with a single
``%pycat``.


Windows timing implementation: Switch to process_time
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Timing on Windows, for example with ``%%time``,  was changed from being based on ``time.perf_counter`` 
(which counted time even when the process was sleeping) to being based on ``time.process_time`` instead 
(which only counts CPU time). This brings it closer to the behavior on Linux. See :ghpull:`12984`.

Miscellaneous
~~~~~~~~~~~~~
 - Non-text formatters are not disabled in the terminal, which should simplify
   writing extensions displaying images or other mimetypes in supporting terminals.
   :ghpull:`12315`
 - It is now possible to automatically insert matching brackets in Terminal IPython using the
   ``TerminalInteractiveShell.auto_match=True`` option. :ghpull:`12586`
 - We are thinking of deprecating the current ``%%javascript`` magic in favor of a better replacement. See :ghpull:`13376`.
 - ``~`` is now expanded when part of a path in most magics :ghpull:`13385`
 - ``%/%%timeit`` magic now adds a comma every thousands to make reading a long number easier :ghpull:`13379`
 - ``"info"`` messages can now be customised to hide some fields :ghpull:`13343`
 - ``collections.UserList`` now pretty-prints :ghpull:`13320`
 - The debugger now has a persistent history, which should make it less
   annoying to retype commands :ghpull:`13246`
 - ``!pip`` ``!conda`` ``!cd`` or ``!ls`` are likely doing the wrong thing. We
   now warn users if they use one of those commands. :ghpull:`12954`
 - Make ``%precision`` work for ``numpy.float64`` type :ghpull:`12902`

Re-added support for XDG config directories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

XDG support through the years comes and goes. There is a tension between having
an identical location for configuration in all platforms versus having simple instructions. 
After initial failures a couple of years ago, IPython was modified to automatically migrate XDG
config files back into ``~/.ipython``. That migration code has now been removed.
IPython now checks the XDG locations, so if you _manually_ move your config
files to your preferred location, IPython will not move them back.


Preparing for Python 3.10
-------------------------

To prepare for Python 3.10, we have started working on removing reliance and
any dependency that is not compatible with Python 3.10. This includes migrating our
test suite to pytest and starting to remove nose. This also means that the
``iptest`` command is now gone and all testing is via pytest.

This was in large part thanks to the NumFOCUS Small Developer grant, which enabled us to
allocate \$4000 to hire `Nikita Kniazev (@Kojoley) <https://github.com/Kojoley>`_,
who did a fantastic job at updating our code base, migrating to pytest, pushing
our coverage, and fixing a large number of bugs. I highly recommend contacting
them if you need help with C++ and Python projects.

You can find all relevant issues and PRs with the SDG 2021 tag `<https://github.com/ipython/ipython/issues?q=label%3A%22Numfocus+SDG+2021%22+>`__

Removing support for older Python versions
------------------------------------------


We are removing support for Python up through 3.7, allowing internal code to use the more
efficient ``pathlib`` and to make better use of type annotations. 

.. image:: ../_images/8.0/pathlib_pathlib_everywhere.jpg
   :alt: "Meme image of Toy Story with Woody and Buzz, with the text 'pathlib, pathlib everywhere'"


We had about 34 PRs only to update some logic to update some functions from managing strings to
using Pathlib.

The completer has also seen significant updates and now makes use of newer Jedi APIs,
offering faster and more reliable tab completion.

Misc Statistics
---------------

Here are some numbers::

    7.x: 296 files, 12561 blank lines, 20282 comments, 35142 line of code.
    8.0: 252 files, 12053 blank lines, 19232 comments, 34505 line of code.

    $ git diff --stat 7.x...master  | tail -1
    340 files changed, 13399 insertions(+), 12421 deletions(-)

We have commits from 162 authors, who contributed 1916 commits in 23 month, excluding merges (to not bias toward
maintainers pushing buttons).::

   $ git shortlog  -s --no-merges  7.x...master | sort -nr
   535	Matthias Bussonnier
    86	Nikita Kniazev
    69	Blazej Michalik
    49	Samuel Gaist
    27	Itamar Turner-Trauring
    18	Spas Kalaydzhisyki
    17	Thomas Kluyver
    17	Quentin Peter
    17	James Morris
    17	Artur Svistunov
    15	Bart Skowron
    14	Alex Hall
    13	rushabh-v
    13	Terry Davis
    13	Benjamin Ragan-Kelley
     8	martinRenou
     8	farisachugthai
     7	dswij
     7	Gal B
     7	Corentin Cadiou
     6	yuji96
     6	Martin Skarzynski
     6	Justin Palmer
     6	Daniel Goldfarb
     6	Ben Greiner
     5	Sammy Al Hashemi
     5	Paul Ivanov
     5	Inception95
     5	Eyenpi
     5	Douglas Blank
     5	Coco Mishra
     5	Bibo Hao
     5	André A. Gomes
     5	Ahmed Fasih
     4	takuya fujiwara
     4	palewire
     4	Thomas A Caswell
     4	Talley Lambert
     4	Scott Sanderson
     4	Ram Rachum
     4	Nick Muoh
     4	Nathan Goldbaum
     4	Mithil Poojary
     4	Michael T
     4	Jakub Klus
     4	Ian Castleden
     4	Eli Rykoff
     4	Ashwin Vishnu
     3	谭九鼎
     3	sleeping
     3	Sylvain Corlay
     3	Peter Corke
     3	Paul Bissex
     3	Matthew Feickert
     3	Fernando Perez
     3	Eric Wieser
     3	Daniel Mietchen
     3	Aditya Sathe
     3	007vedant
     2	rchiodo
     2	nicolaslazo
     2	luttik
     2	gorogoroumaru
     2	foobarbyte
     2	bar-hen
     2	Theo Ouzhinski
     2	Strawkage
     2	Samreen Zarroug
     2	Pete Blois
     2	Meysam Azad
     2	Matthieu Ancellin
     2	Mark Schmitz
     2	Maor Kleinberger
     2	MRCWirtz
     2	Lumir Balhar
     2	Julien Rabinow
     2	Juan Luis Cano Rodríguez
     2	Joyce Er
     2	Jakub
     2	Faris A Chugthai
     2	Ethan Madden
     2	Dimitri Papadopoulos
     2	Diego Fernandez
     2	Daniel Shimon
     2	Coco Bennett
     2	Carlos Cordoba
     2	Boyuan Liu
     2	BaoGiang HoangVu
     2	Augusto
     2	Arthur Svistunov
     2	Arthur Moreira
     2	Ali Nabipour
     2	Adam Hackbarth
     1	richard
     1	linar-jether
     1	lbennett
     1	juacrumar
     1	gpotter2
     1	digitalvirtuoso
     1	dalthviz
     1	Yonatan Goldschmidt
     1	Tomasz Kłoczko
     1	Tobias Bengfort
     1	Timur Kushukov
     1	Thomas
     1	Snir Broshi
     1	Shao Yang Hong
     1	Sanjana-03
     1	Romulo Filho
     1	Rodolfo Carvalho
     1	Richard Shadrach
     1	Reilly Tucker Siemens
     1	Rakessh Roshan
     1	Piers Titus van der Torren
     1	PhanatosZou
     1	Pavel Safronov
     1	Paulo S. Costa
     1	Paul McCarthy
     1	NotWearingPants
     1	Naelson Douglas
     1	Michael Tiemann
     1	Matt Wozniski
     1	Markus Wageringel
     1	Marcus Wirtz
     1	Marcio Mazza
     1	Lumír 'Frenzy' Balhar
     1	Lightyagami1
     1	Leon Anavi
     1	LeafyLi
     1	L0uisJ0shua
     1	Kyle Cutler
     1	Krzysztof Cybulski
     1	Kevin Kirsche
     1	KIU Shueng Chuan
     1	Jonathan Slenders
     1	Jay Qi
     1	Jake VanderPlas
     1	Iwan Briquemont
     1	Hussaina Begum Nandyala
     1	Gordon Ball
     1	Gabriel Simonetto
     1	Frank Tobia
     1	Erik
     1	Elliott Sales de Andrade
     1	Daniel Hahler
     1	Dan Green-Leipciger
     1	Dan Green
     1	Damian Yurzola
     1	Coon, Ethan T
     1	Carol Willing
     1	Brian Lee
     1	Brendan Gerrity
     1	Blake Griffin
     1	Bastian Ebeling
     1	Bartosz Telenczuk
     1	Ankitsingh6299
     1	Andrew Port
     1	Andrew J. Hesford
     1	Albert Zhang
     1	Adam Johnson

This does not, of course, represent non-code contributions, for which we are also grateful.


API Changes using Frappuccino
-----------------------------

This is an experimental exhaustive API difference using `Frappuccino <https://pypi.org/project/frappuccino/>`_


The following items are new in IPython 8.0 ::

    + IPython.core.async_helpers.get_asyncio_loop()
    + IPython.core.completer.Dict
    + IPython.core.completer.Pattern
    + IPython.core.completer.Sequence
    + IPython.core.completer.__skip_doctest__
    + IPython.core.debugger.Pdb.precmd(self, line)
    + IPython.core.debugger.__skip_doctest__
    + IPython.core.display.__getattr__(name)
    + IPython.core.display.warn
    + IPython.core.display_functions
    + IPython.core.display_functions.DisplayHandle
    + IPython.core.display_functions.DisplayHandle.display(self, obj, **kwargs)
    + IPython.core.display_functions.DisplayHandle.update(self, obj, **kwargs)
    + IPython.core.display_functions.__all__
    + IPython.core.display_functions.__builtins__
    + IPython.core.display_functions.__cached__
    + IPython.core.display_functions.__doc__
    + IPython.core.display_functions.__file__
    + IPython.core.display_functions.__loader__
    + IPython.core.display_functions.__name__
    + IPython.core.display_functions.__package__
    + IPython.core.display_functions.__spec__
    + IPython.core.display_functions.b2a_hex
    + IPython.core.display_functions.clear_output(wait=False)
    + IPython.core.display_functions.display(*objs, include='None', exclude='None', metadata='None', transient='None', display_id='None', raw=False, clear=False, **kwargs)
    + IPython.core.display_functions.publish_display_data(data, metadata='None', source='<deprecated>', *, transient='None', **kwargs)
    + IPython.core.display_functions.update_display(obj, *, display_id, **kwargs)
    + IPython.core.extensions.BUILTINS_EXTS
    + IPython.core.inputtransformer2.has_sunken_brackets(tokens)
    + IPython.core.interactiveshell.Callable
    + IPython.core.interactiveshell.__annotations__
    + IPython.core.ultratb.List
    + IPython.core.ultratb.Tuple
    + IPython.lib.pretty.CallExpression
    + IPython.lib.pretty.CallExpression.factory(name)
    + IPython.lib.pretty.RawStringLiteral
    + IPython.lib.pretty.RawText
    + IPython.terminal.debugger.TerminalPdb.do_interact(self, arg)
    + IPython.terminal.embed.Set

The following items have been removed (or moved to superclass)::

    - IPython.core.application.BaseIPythonApplication.initialize_subcommand
    - IPython.core.completer.Sentinel
    - IPython.core.completer.skip_doctest
    - IPython.core.debugger.Tracer
    - IPython.core.display.DisplayHandle
    - IPython.core.display.DisplayHandle.display
    - IPython.core.display.DisplayHandle.update
    - IPython.core.display.b2a_hex
    - IPython.core.display.clear_output
    - IPython.core.display.display
    - IPython.core.display.publish_display_data
    - IPython.core.display.update_display
    - IPython.core.excolors.Deprec
    - IPython.core.excolors.ExceptionColors
    - IPython.core.history.warn
    - IPython.core.hooks.late_startup_hook
    - IPython.core.hooks.pre_run_code_hook
    - IPython.core.hooks.shutdown_hook
    - IPython.core.interactiveshell.InteractiveShell.init_deprecation_warnings
    - IPython.core.interactiveshell.InteractiveShell.init_readline
    - IPython.core.interactiveshell.InteractiveShell.write
    - IPython.core.interactiveshell.InteractiveShell.write_err
    - IPython.core.interactiveshell.get_default_colors
    - IPython.core.interactiveshell.removed_co_newlocals
    - IPython.core.magics.execution.ExecutionMagics.profile_missing_notice
    - IPython.core.magics.script.PIPE
    - IPython.core.prefilter.PrefilterManager.init_transformers
    - IPython.core.release.classifiers
    - IPython.core.release.description
    - IPython.core.release.keywords
    - IPython.core.release.long_description
    - IPython.core.release.name
    - IPython.core.release.platforms
    - IPython.core.release.url
    - IPython.core.ultratb.VerboseTB.format_records
    - IPython.core.ultratb.find_recursion
    - IPython.core.ultratb.findsource
    - IPython.core.ultratb.fix_frame_records_filenames
    - IPython.core.ultratb.inspect_error
    - IPython.core.ultratb.is_recursion_error
    - IPython.core.ultratb.with_patch_inspect
    - IPython.external.__all__
    - IPython.external.__builtins__
    - IPython.external.__cached__
    - IPython.external.__doc__
    - IPython.external.__file__
    - IPython.external.__loader__
    - IPython.external.__name__
    - IPython.external.__package__
    - IPython.external.__path__
    - IPython.external.__spec__
    - IPython.kernel.KernelConnectionInfo
    - IPython.kernel.__builtins__
    - IPython.kernel.__cached__
    - IPython.kernel.__warningregistry__
    - IPython.kernel.pkg
    - IPython.kernel.protocol_version
    - IPython.kernel.protocol_version_info
    - IPython.kernel.src
    - IPython.kernel.version_info
    - IPython.kernel.warn
    - IPython.lib.backgroundjobs
    - IPython.lib.backgroundjobs.BackgroundJobBase
    - IPython.lib.backgroundjobs.BackgroundJobBase.run
    - IPython.lib.backgroundjobs.BackgroundJobBase.traceback
    - IPython.lib.backgroundjobs.BackgroundJobExpr
    - IPython.lib.backgroundjobs.BackgroundJobExpr.call
    - IPython.lib.backgroundjobs.BackgroundJobFunc
    - IPython.lib.backgroundjobs.BackgroundJobFunc.call
    - IPython.lib.backgroundjobs.BackgroundJobManager
    - IPython.lib.backgroundjobs.BackgroundJobManager.flush
    - IPython.lib.backgroundjobs.BackgroundJobManager.new
    - IPython.lib.backgroundjobs.BackgroundJobManager.remove
    - IPython.lib.backgroundjobs.BackgroundJobManager.result
    - IPython.lib.backgroundjobs.BackgroundJobManager.status
    - IPython.lib.backgroundjobs.BackgroundJobManager.traceback
    - IPython.lib.backgroundjobs.__builtins__
    - IPython.lib.backgroundjobs.__cached__
    - IPython.lib.backgroundjobs.__doc__
    - IPython.lib.backgroundjobs.__file__
    - IPython.lib.backgroundjobs.__loader__
    - IPython.lib.backgroundjobs.__name__
    - IPython.lib.backgroundjobs.__package__
    - IPython.lib.backgroundjobs.__spec__
    - IPython.lib.kernel.__builtins__
    - IPython.lib.kernel.__cached__
    - IPython.lib.kernel.__doc__
    - IPython.lib.kernel.__file__
    - IPython.lib.kernel.__loader__
    - IPython.lib.kernel.__name__
    - IPython.lib.kernel.__package__
    - IPython.lib.kernel.__spec__
    - IPython.lib.kernel.__warningregistry__
    - IPython.paths.fs_encoding
    - IPython.terminal.debugger.DEFAULT_BUFFER
    - IPython.terminal.debugger.cursor_in_leading_ws
    - IPython.terminal.debugger.emacs_insert_mode
    - IPython.terminal.debugger.has_selection
    - IPython.terminal.debugger.vi_insert_mode
    - IPython.terminal.interactiveshell.DISPLAY_BANNER_DEPRECATED
    - IPython.terminal.ipapp.TerminalIPythonApp.parse_command_line
    - IPython.testing.test
    - IPython.utils.contexts.NoOpContext
    - IPython.utils.io.IOStream
    - IPython.utils.io.IOStream.close
    - IPython.utils.io.IOStream.write
    - IPython.utils.io.IOStream.writelines
    - IPython.utils.io.__warningregistry__
    - IPython.utils.io.atomic_writing
    - IPython.utils.io.stderr
    - IPython.utils.io.stdin
    - IPython.utils.io.stdout
    - IPython.utils.io.unicode_std_stream
    - IPython.utils.path.get_ipython_cache_dir
    - IPython.utils.path.get_ipython_dir
    - IPython.utils.path.get_ipython_module_path
    - IPython.utils.path.get_ipython_package_dir
    - IPython.utils.path.locate_profile
    - IPython.utils.path.unquote_filename
    - IPython.utils.py3compat.PY2
    - IPython.utils.py3compat.PY3
    - IPython.utils.py3compat.buffer_to_bytes
    - IPython.utils.py3compat.builtin_mod_name
    - IPython.utils.py3compat.cast_bytes
    - IPython.utils.py3compat.getcwd
    - IPython.utils.py3compat.isidentifier
    - IPython.utils.py3compat.u_format

The following signatures differ between 7.x and 8.0::

    - IPython.core.completer.IPCompleter.unicode_name_matches(self, text)
    + IPython.core.completer.IPCompleter.unicode_name_matches(text)

    - IPython.core.completer.match_dict_keys(keys, prefix, delims)
    + IPython.core.completer.match_dict_keys(keys, prefix, delims, extra_prefix='None')

    - IPython.core.interactiveshell.InteractiveShell.object_inspect_mime(self, oname, detail_level=0)
    + IPython.core.interactiveshell.InteractiveShell.object_inspect_mime(self, oname, detail_level=0, omit_sections='()')

    - IPython.core.interactiveshell.InteractiveShell.set_hook(self, name, hook, priority=50, str_key='None', re_key='None', _warn_deprecated=True)
    + IPython.core.interactiveshell.InteractiveShell.set_hook(self, name, hook, priority=50, str_key='None', re_key='None')

    - IPython.core.oinspect.Inspector.info(self, obj, oname='', formatter='None', info='None', detail_level=0)
    + IPython.core.oinspect.Inspector.info(self, obj, oname='', info='None', detail_level=0)

    - IPython.core.oinspect.Inspector.pinfo(self, obj, oname='', formatter='None', info='None', detail_level=0, enable_html_pager=True)
    + IPython.core.oinspect.Inspector.pinfo(self, obj, oname='', formatter='None', info='None', detail_level=0, enable_html_pager=True, omit_sections='()')

    - IPython.core.profiledir.ProfileDir.copy_config_file(self, config_file, path='None', overwrite=False)
    + IPython.core.profiledir.ProfileDir.copy_config_file(self, config_file, path, overwrite=False)

    - IPython.core.ultratb.VerboseTB.format_record(self, frame, file, lnum, func, lines, index)
    + IPython.core.ultratb.VerboseTB.format_record(self, frame_info)

    - IPython.terminal.embed.InteractiveShellEmbed.mainloop(self, local_ns='None', module='None', stack_depth=0, display_banner='None', global_ns='None', compile_flags='None')
    + IPython.terminal.embed.InteractiveShellEmbed.mainloop(self, local_ns='None', module='None', stack_depth=0, compile_flags='None')

    - IPython.terminal.embed.embed(**kwargs)
    + IPython.terminal.embed.embed(*, header='', compile_flags='None', **kwargs)

    - IPython.terminal.interactiveshell.TerminalInteractiveShell.interact(self, display_banner='<object object at 0xffffff>')
    + IPython.terminal.interactiveshell.TerminalInteractiveShell.interact(self)

    - IPython.terminal.interactiveshell.TerminalInteractiveShell.mainloop(self, display_banner='<object object at 0xffffff>')
    + IPython.terminal.interactiveshell.TerminalInteractiveShell.mainloop(self)

    - IPython.utils.path.get_py_filename(name, force_win32='None')
    + IPython.utils.path.get_py_filename(name)

The following are new attributes (that might be inherited)::

    + IPython.core.completer.IPCompleter.unicode_names
    + IPython.core.debugger.InterruptiblePdb.precmd
    + IPython.core.debugger.Pdb.precmd
    + IPython.core.ultratb.AutoFormattedTB.has_colors
    + IPython.core.ultratb.ColorTB.has_colors
    + IPython.core.ultratb.FormattedTB.has_colors
    + IPython.core.ultratb.ListTB.has_colors
    + IPython.core.ultratb.SyntaxTB.has_colors
    + IPython.core.ultratb.TBTools.has_colors
    + IPython.core.ultratb.VerboseTB.has_colors
    + IPython.terminal.debugger.TerminalPdb.do_interact
    + IPython.terminal.debugger.TerminalPdb.precmd

The following attribute/methods have been removed::

    - IPython.core.application.BaseIPythonApplication.deprecated_subcommands
    - IPython.core.ultratb.AutoFormattedTB.format_records
    - IPython.core.ultratb.ColorTB.format_records
    - IPython.core.ultratb.FormattedTB.format_records
    - IPython.terminal.embed.InteractiveShellEmbed.init_deprecation_warnings
    - IPython.terminal.embed.InteractiveShellEmbed.init_readline
    - IPython.terminal.embed.InteractiveShellEmbed.write
    - IPython.terminal.embed.InteractiveShellEmbed.write_err
    - IPython.terminal.interactiveshell.TerminalInteractiveShell.init_deprecation_warnings
    - IPython.terminal.interactiveshell.TerminalInteractiveShell.init_readline
    - IPython.terminal.interactiveshell.TerminalInteractiveShell.write
    - IPython.terminal.interactiveshell.TerminalInteractiveShell.write_err
    - IPython.terminal.ipapp.LocateIPythonApp.deprecated_subcommands
    - IPython.terminal.ipapp.LocateIPythonApp.initialize_subcommand
    - IPython.terminal.ipapp.TerminalIPythonApp.deprecated_subcommands
    - IPython.terminal.ipapp.TerminalIPythonApp.initialize_subcommand
.. _issues_list_011:

Issues closed in the 0.11 development cycle
===========================================

In this cycle, we closed a total of 511 issues, 226 pull requests and 285
regular issues; this is the full list (generated with the script
`tools/github_stats.py`). We should note that a few of these were made on the
0.10.x series, but we have no automatic way of filtering the issues by branch,
so this reflects all of our development over the last two years, including work
already released in 0.10.2:

Pull requests (226):

* `620 <https://github.com/ipython/ipython/issues/620>`_: Release notes and updates to GUI support docs for 0.11
* `642 <https://github.com/ipython/ipython/issues/642>`_: fix typo in docs/examples/vim/README.rst
* `631 <https://github.com/ipython/ipython/issues/631>`_: two-way vim-ipython integration
* `637 <https://github.com/ipython/ipython/issues/637>`_: print is a function, this allows to properly exit ipython
* `635 <https://github.com/ipython/ipython/issues/635>`_: support html representations in the notebook frontend
* `639 <https://github.com/ipython/ipython/issues/639>`_: Updating the credits file
* `628 <https://github.com/ipython/ipython/issues/628>`_: import pexpect from IPython.external in irunner
* `596 <https://github.com/ipython/ipython/issues/596>`_: Irunner
* `598 <https://github.com/ipython/ipython/issues/598>`_: Fix templates for CrashHandler
* `590 <https://github.com/ipython/ipython/issues/590>`_: Desktop
* `600 <https://github.com/ipython/ipython/issues/600>`_: Fix bug with non-ascii reprs inside pretty-printed lists.
* `618 <https://github.com/ipython/ipython/issues/618>`_: I617
* `599 <https://github.com/ipython/ipython/issues/599>`_: Gui Qt example and docs
* `619 <https://github.com/ipython/ipython/issues/619>`_: manpage update
* `582 <https://github.com/ipython/ipython/issues/582>`_: Updating sympy profile to match the exec_lines of isympy.
* `578 <https://github.com/ipython/ipython/issues/578>`_: Check to see if correct source for decorated functions can be displayed
* `589 <https://github.com/ipython/ipython/issues/589>`_: issue 588
* `591 <https://github.com/ipython/ipython/issues/591>`_: simulate shell expansion on %run arguments, at least tilde expansion
* `576 <https://github.com/ipython/ipython/issues/576>`_: Show message about %paste magic on an IndentationError
* `574 <https://github.com/ipython/ipython/issues/574>`_: Getcwdu
* `565 <https://github.com/ipython/ipython/issues/565>`_: don't move old config files, keep nagging the user
* `575 <https://github.com/ipython/ipython/issues/575>`_: Added more docstrings to IPython.zmq.session.
* `567 <https://github.com/ipython/ipython/issues/567>`_: fix trailing whitespace from resetting indentation
* `564 <https://github.com/ipython/ipython/issues/564>`_: Command line args in docs
* `560 <https://github.com/ipython/ipython/issues/560>`_: reorder qt support in kernel
* `561 <https://github.com/ipython/ipython/issues/561>`_: command-line suggestions
* `556 <https://github.com/ipython/ipython/issues/556>`_: qt_for_kernel: use matplotlib rcParams to decide between PyQt4 and PySide
* `557 <https://github.com/ipython/ipython/issues/557>`_: Update usage.py to newapp
* `555 <https://github.com/ipython/ipython/issues/555>`_: Rm default old config
* `552 <https://github.com/ipython/ipython/issues/552>`_: update parallel code for py3k
* `504 <https://github.com/ipython/ipython/issues/504>`_: Updating string formatting
* `551 <https://github.com/ipython/ipython/issues/551>`_: Make pylab import all configurable
* `496 <https://github.com/ipython/ipython/issues/496>`_: Qt editing keybindings
* `550 <https://github.com/ipython/ipython/issues/550>`_: Support v2 PyQt4 APIs and PySide in kernel's GUI support
* `546 <https://github.com/ipython/ipython/issues/546>`_: doc update
* `548 <https://github.com/ipython/ipython/issues/548>`_: Fix sympy profile to work with sympy 0.7.
* `542 <https://github.com/ipython/ipython/issues/542>`_: issue 440
* `533 <https://github.com/ipython/ipython/issues/533>`_: Remove unused configobj and validate libraries from externals.
* `538 <https://github.com/ipython/ipython/issues/538>`_: fix various tests on Windows
* `540 <https://github.com/ipython/ipython/issues/540>`_: support `-pylab` flag with deprecation warning
* `537 <https://github.com/ipython/ipython/issues/537>`_: Docs update
* `536 <https://github.com/ipython/ipython/issues/536>`_: `setup.py install` depends on setuptools on Windows
* `480 <https://github.com/ipython/ipython/issues/480>`_: Get help mid-command
* `462 <https://github.com/ipython/ipython/issues/462>`_: Str and Bytes traitlets
* `534 <https://github.com/ipython/ipython/issues/534>`_: Handle unicode properly in IPython.zmq.iostream
* `527 <https://github.com/ipython/ipython/issues/527>`_: ZMQ displayhook
* `526 <https://github.com/ipython/ipython/issues/526>`_: Handle asynchronous output in Qt console
* `528 <https://github.com/ipython/ipython/issues/528>`_: Do not import deprecated functions from external decorators library.
* `454 <https://github.com/ipython/ipython/issues/454>`_: New BaseIPythonApplication
* `532 <https://github.com/ipython/ipython/issues/532>`_: Zmq unicode
* `531 <https://github.com/ipython/ipython/issues/531>`_: Fix Parallel test
* `525 <https://github.com/ipython/ipython/issues/525>`_: fallback on lsof if otool not found in libedit detection
* `517 <https://github.com/ipython/ipython/issues/517>`_: Merge IPython.parallel.streamsession into IPython.zmq.session
* `521 <https://github.com/ipython/ipython/issues/521>`_: use dict.get(key) instead of dict[key] for safety from KeyErrors
* `492 <https://github.com/ipython/ipython/issues/492>`_: add QtConsoleApp using newapplication
* `485 <https://github.com/ipython/ipython/issues/485>`_: terminal IPython with newapp
* `486 <https://github.com/ipython/ipython/issues/486>`_: Use newapp in parallel code
* `511 <https://github.com/ipython/ipython/issues/511>`_: Add a new line before displaying multiline strings in the Qt console.
* `509 <https://github.com/ipython/ipython/issues/509>`_: i508
* `501 <https://github.com/ipython/ipython/issues/501>`_: ignore EINTR in channel loops
* `495 <https://github.com/ipython/ipython/issues/495>`_: Better selection of Qt bindings when QT_API is not specified
* `498 <https://github.com/ipython/ipython/issues/498>`_: Check for .pyd as extension for binary files.
* `494 <https://github.com/ipython/ipython/issues/494>`_: QtConsole zoom adjustments
* `490 <https://github.com/ipython/ipython/issues/490>`_: fix UnicodeEncodeError writing SVG string to .svg file, fixes #489
* `491 <https://github.com/ipython/ipython/issues/491>`_: add QtConsoleApp using newapplication
* `479 <https://github.com/ipython/ipython/issues/479>`_: embed() doesn't load default config
* `483 <https://github.com/ipython/ipython/issues/483>`_: Links launchpad -> github
* `419 <https://github.com/ipython/ipython/issues/419>`_: %xdel magic
* `477 <https://github.com/ipython/ipython/issues/477>`_: Add \n to lines in the log
* `459 <https://github.com/ipython/ipython/issues/459>`_: use os.system for shell.system in Terminal frontend
* `475 <https://github.com/ipython/ipython/issues/475>`_: i473
* `471 <https://github.com/ipython/ipython/issues/471>`_: Add test decorator onlyif_unicode_paths.
* `474 <https://github.com/ipython/ipython/issues/474>`_: Fix support for raw GTK and WX matplotlib backends.
* `472 <https://github.com/ipython/ipython/issues/472>`_: Kernel event loop is robust against random SIGINT.
* `460 <https://github.com/ipython/ipython/issues/460>`_: Share code for magic_edit
* `469 <https://github.com/ipython/ipython/issues/469>`_: Add exit code when running all tests with iptest.
* `464 <https://github.com/ipython/ipython/issues/464>`_: Add home directory expansion to IPYTHON_DIR environment variables.
* `455 <https://github.com/ipython/ipython/issues/455>`_: Bugfix with logger
* `448 <https://github.com/ipython/ipython/issues/448>`_: Separate out skip_doctest decorator
* `453 <https://github.com/ipython/ipython/issues/453>`_: Draft of new main BaseIPythonApplication.
* `452 <https://github.com/ipython/ipython/issues/452>`_: Use list/tuple/dict/set subclass's overridden __repr__ instead of the pretty
* `398 <https://github.com/ipython/ipython/issues/398>`_: allow toggle of svg/png inline figure format
* `381 <https://github.com/ipython/ipython/issues/381>`_: Support inline PNGs of matplotlib plots
* `413 <https://github.com/ipython/ipython/issues/413>`_: Retries and Resubmit (#411 and #412)
* `370 <https://github.com/ipython/ipython/issues/370>`_: Fixes to the display system
* `449 <https://github.com/ipython/ipython/issues/449>`_: Fix issue 447 - inspecting old-style classes.
* `423 <https://github.com/ipython/ipython/issues/423>`_: Allow type checking on elements of List,Tuple,Set traits
* `400 <https://github.com/ipython/ipython/issues/400>`_: Config5
* `421 <https://github.com/ipython/ipython/issues/421>`_: Generalise mechanism to put text at the next prompt in the Qt console.
* `443 <https://github.com/ipython/ipython/issues/443>`_: pinfo code duplication
* `429 <https://github.com/ipython/ipython/issues/429>`_: add check_pid, and handle stale PID info in ipcluster.
* `431 <https://github.com/ipython/ipython/issues/431>`_: Fix error message in test_irunner
* `427 <https://github.com/ipython/ipython/issues/427>`_: handle different SyntaxError messages in test_irunner
* `424 <https://github.com/ipython/ipython/issues/424>`_: Irunner test failure
* `430 <https://github.com/ipython/ipython/issues/430>`_: Small parallel doc typo
* `422 <https://github.com/ipython/ipython/issues/422>`_: Make ipython-qtconsole a GUI script
* `420 <https://github.com/ipython/ipython/issues/420>`_: Permit kernel std* to be redirected
* `408 <https://github.com/ipython/ipython/issues/408>`_: History request
* `388 <https://github.com/ipython/ipython/issues/388>`_: Add Emacs-style kill ring to Qt console
* `414 <https://github.com/ipython/ipython/issues/414>`_: Warn on old config files
* `415 <https://github.com/ipython/ipython/issues/415>`_: Prevent prefilter from crashing IPython
* `418 <https://github.com/ipython/ipython/issues/418>`_: Minor configuration doc fixes
* `407 <https://github.com/ipython/ipython/issues/407>`_: Update What's new documentation
* `410 <https://github.com/ipython/ipython/issues/410>`_: Install notebook frontend
* `406 <https://github.com/ipython/ipython/issues/406>`_: install IPython.zmq.gui
* `393 <https://github.com/ipython/ipython/issues/393>`_: ipdir unicode
* `397 <https://github.com/ipython/ipython/issues/397>`_: utils.io.Term.cin/out/err -> utils.io.stdin/out/err
* `389 <https://github.com/ipython/ipython/issues/389>`_: DB fixes and Scheduler HWM
* `374 <https://github.com/ipython/ipython/issues/374>`_: Various Windows-related fixes to IPython.parallel
* `362 <https://github.com/ipython/ipython/issues/362>`_: fallback on defaultencoding if filesystemencoding is None
* `382 <https://github.com/ipython/ipython/issues/382>`_: Shell's reset method clears namespace from last %run command.
* `385 <https://github.com/ipython/ipython/issues/385>`_: Update iptest exclusions (fix #375)
* `383 <https://github.com/ipython/ipython/issues/383>`_: Catch errors in querying readline which occur with pyreadline.
* `373 <https://github.com/ipython/ipython/issues/373>`_: Remove runlines etc.
* `364 <https://github.com/ipython/ipython/issues/364>`_: Single output
* `372 <https://github.com/ipython/ipython/issues/372>`_: Multiline input push
* `363 <https://github.com/ipython/ipython/issues/363>`_: Issue 125
* `361 <https://github.com/ipython/ipython/issues/361>`_: don't rely on setuptools for readline dependency check
* `349 <https://github.com/ipython/ipython/issues/349>`_: Fix %autopx magic
* `355 <https://github.com/ipython/ipython/issues/355>`_: History save thread
* `356 <https://github.com/ipython/ipython/issues/356>`_: Usability improvements to history in Qt console
* `357 <https://github.com/ipython/ipython/issues/357>`_: Exit autocall
* `353 <https://github.com/ipython/ipython/issues/353>`_: Rewrite quit()/exit()/Quit()/Exit() calls as magic
* `354 <https://github.com/ipython/ipython/issues/354>`_: Cell tweaks
* `345 <https://github.com/ipython/ipython/issues/345>`_: Attempt to address (partly) issue ipython/#342 by rewriting quit(), exit(), etc.
* `352 <https://github.com/ipython/ipython/issues/352>`_: #342: Try to recover as intelligently as possible if user calls magic().
* `346 <https://github.com/ipython/ipython/issues/346>`_: Dedent prefix bugfix + tests: #142
* `348 <https://github.com/ipython/ipython/issues/348>`_: %reset doesn't reset prompt number.
* `347 <https://github.com/ipython/ipython/issues/347>`_: Make ip.reset() work the same in interactive or non-interactive code.
* `343 <https://github.com/ipython/ipython/issues/343>`_: make readline a dependency on OSX
* `344 <https://github.com/ipython/ipython/issues/344>`_: restore auto debug behavior
* `339 <https://github.com/ipython/ipython/issues/339>`_: fix for issue 337: incorrect/phantom tooltips for magics
* `254 <https://github.com/ipython/ipython/issues/254>`_: newparallel branch (add zmq.parallel submodule)
* `334 <https://github.com/ipython/ipython/issues/334>`_: Hard reset
* `316 <https://github.com/ipython/ipython/issues/316>`_: Unicode win process
* `332 <https://github.com/ipython/ipython/issues/332>`_: AST splitter
* `325 <https://github.com/ipython/ipython/issues/325>`_: Removetwisted
* `330 <https://github.com/ipython/ipython/issues/330>`_: Magic pastebin
* `309 <https://github.com/ipython/ipython/issues/309>`_: Bug tests for GH Issues 238, 284, 306, 307. Skip module machinery if not installed. Known failures reported as 'K'
* `331 <https://github.com/ipython/ipython/issues/331>`_: Tweak config loader for PyPy compatibility.
* `319 <https://github.com/ipython/ipython/issues/319>`_: Rewrite code to restore readline history after an action
* `329 <https://github.com/ipython/ipython/issues/329>`_: Do not store file contents in history when running a .ipy file.
* `179 <https://github.com/ipython/ipython/issues/179>`_: Html notebook
* `323 <https://github.com/ipython/ipython/issues/323>`_: Add missing external.pexpect to packages
* `295 <https://github.com/ipython/ipython/issues/295>`_: Magic local scope
* `315 <https://github.com/ipython/ipython/issues/315>`_: Unicode magic args
* `310 <https://github.com/ipython/ipython/issues/310>`_: allow Unicode Command-Line options
* `313 <https://github.com/ipython/ipython/issues/313>`_: Readline shortcuts
* `311 <https://github.com/ipython/ipython/issues/311>`_: Qtconsole exit
* `312 <https://github.com/ipython/ipython/issues/312>`_: History memory
* `294 <https://github.com/ipython/ipython/issues/294>`_: Issue 290
* `292 <https://github.com/ipython/ipython/issues/292>`_: Issue 31
* `252 <https://github.com/ipython/ipython/issues/252>`_: Unicode issues
* `235 <https://github.com/ipython/ipython/issues/235>`_: Fix history magic command's bugs wrt to full history and add -O option to display full history
* `236 <https://github.com/ipython/ipython/issues/236>`_: History minus p flag
* `261 <https://github.com/ipython/ipython/issues/261>`_: Adapt magic commands to new history system.
* `282 <https://github.com/ipython/ipython/issues/282>`_: SQLite history
* `191 <https://github.com/ipython/ipython/issues/191>`_: Unbundle external libraries
* `199 <https://github.com/ipython/ipython/issues/199>`_: Magic arguments
* `204 <https://github.com/ipython/ipython/issues/204>`_: Emacs completion bugfix
* `293 <https://github.com/ipython/ipython/issues/293>`_: Issue 133
* `249 <https://github.com/ipython/ipython/issues/249>`_: Writing unicode characters to a log file. (IPython 0.10.2.git)
* `283 <https://github.com/ipython/ipython/issues/283>`_: Support for 256-color escape sequences in Qt console
* `281 <https://github.com/ipython/ipython/issues/281>`_: Refactored and improved Qt console's HTML export facility
* `237 <https://github.com/ipython/ipython/issues/237>`_: Fix185 (take two)
* `251 <https://github.com/ipython/ipython/issues/251>`_: Issue 129
* `278 <https://github.com/ipython/ipython/issues/278>`_: add basic XDG_CONFIG_HOME support
* `275 <https://github.com/ipython/ipython/issues/275>`_: inline pylab cuts off labels on log plots
* `280 <https://github.com/ipython/ipython/issues/280>`_: Add %precision magic
* `259 <https://github.com/ipython/ipython/issues/259>`_: Pyside support
* `193 <https://github.com/ipython/ipython/issues/193>`_: Make ipython cProfile-able
* `272 <https://github.com/ipython/ipython/issues/272>`_: Magic examples
* `219 <https://github.com/ipython/ipython/issues/219>`_: Doc magic pycat
* `221 <https://github.com/ipython/ipython/issues/221>`_: Doc magic alias
* `230 <https://github.com/ipython/ipython/issues/230>`_: Doc magic edit
* `224 <https://github.com/ipython/ipython/issues/224>`_: Doc magic cpaste
* `229 <https://github.com/ipython/ipython/issues/229>`_: Doc magic pdef
* `273 <https://github.com/ipython/ipython/issues/273>`_: Docs build
* `228 <https://github.com/ipython/ipython/issues/228>`_: Doc magic who
* `233 <https://github.com/ipython/ipython/issues/233>`_: Doc magic cd
* `226 <https://github.com/ipython/ipython/issues/226>`_: Doc magic pwd
* `218 <https://github.com/ipython/ipython/issues/218>`_: Doc magic history
* `231 <https://github.com/ipython/ipython/issues/231>`_: Doc magic reset
* `225 <https://github.com/ipython/ipython/issues/225>`_: Doc magic save
* `222 <https://github.com/ipython/ipython/issues/222>`_: Doc magic timeit
* `223 <https://github.com/ipython/ipython/issues/223>`_: Doc magic colors
* `203 <https://github.com/ipython/ipython/issues/203>`_: Small typos in zmq/blockingkernelmanager.py
* `227 <https://github.com/ipython/ipython/issues/227>`_: Doc magic logon
* `232 <https://github.com/ipython/ipython/issues/232>`_: Doc magic profile
* `264 <https://github.com/ipython/ipython/issues/264>`_: Kernel logging
* `220 <https://github.com/ipython/ipython/issues/220>`_: Doc magic edit
* `268 <https://github.com/ipython/ipython/issues/268>`_: PyZMQ >= 2.0.10
* `267 <https://github.com/ipython/ipython/issues/267>`_: GitHub Pages (again)
* `266 <https://github.com/ipython/ipython/issues/266>`_: OSX-specific fixes to the Qt console
* `255 <https://github.com/ipython/ipython/issues/255>`_: Gitwash typo
* `265 <https://github.com/ipython/ipython/issues/265>`_: Fix string input2
* `260 <https://github.com/ipython/ipython/issues/260>`_: Kernel crash with empty history
* `243 <https://github.com/ipython/ipython/issues/243>`_: New display system
* `242 <https://github.com/ipython/ipython/issues/242>`_: Fix terminal exit
* `250 <https://github.com/ipython/ipython/issues/250>`_: always use Session.send
* `239 <https://github.com/ipython/ipython/issues/239>`_: Makefile command & script for GitHub Pages
* `244 <https://github.com/ipython/ipython/issues/244>`_: My exit
* `234 <https://github.com/ipython/ipython/issues/234>`_: Timed history save
* `217 <https://github.com/ipython/ipython/issues/217>`_: Doc magic lsmagic
* `215 <https://github.com/ipython/ipython/issues/215>`_: History fix
* `195 <https://github.com/ipython/ipython/issues/195>`_: Formatters
* `192 <https://github.com/ipython/ipython/issues/192>`_: Ready colorize bug
* `198 <https://github.com/ipython/ipython/issues/198>`_: Windows workdir
* `174 <https://github.com/ipython/ipython/issues/174>`_: Whitespace cleanup
* `188 <https://github.com/ipython/ipython/issues/188>`_: Version info: update our version management system to use git.
* `158 <https://github.com/ipython/ipython/issues/158>`_: Ready for merge
* `187 <https://github.com/ipython/ipython/issues/187>`_: Resolved Print shortcut collision with ctrl-P emacs binding
* `183 <https://github.com/ipython/ipython/issues/183>`_: cleanup of exit/quit commands for qt console
* `184 <https://github.com/ipython/ipython/issues/184>`_: Logo added to sphinx docs
* `180 <https://github.com/ipython/ipython/issues/180>`_: Cleanup old code
* `171 <https://github.com/ipython/ipython/issues/171>`_: Expose Pygments styles as options
* `170 <https://github.com/ipython/ipython/issues/170>`_: HTML Fixes
* `172 <https://github.com/ipython/ipython/issues/172>`_: Fix del method exit test
* `164 <https://github.com/ipython/ipython/issues/164>`_: Qt frontend shutdown behavior fixes and enhancements
* `167 <https://github.com/ipython/ipython/issues/167>`_: Added HTML export
* `163 <https://github.com/ipython/ipython/issues/163>`_: Execution refactor
* `159 <https://github.com/ipython/ipython/issues/159>`_: Ipy3 preparation
* `155 <https://github.com/ipython/ipython/issues/155>`_: Ready startup fix
* `152 <https://github.com/ipython/ipython/issues/152>`_: 0.10.1 sge
* `151 <https://github.com/ipython/ipython/issues/151>`_: mk_object_info -> object_info
* `149 <https://github.com/ipython/ipython/issues/149>`_: Simple bug-fix

Regular issues (285):

* `630 <https://github.com/ipython/ipython/issues/630>`_: new.py in pwd prevents ipython from starting
* `623 <https://github.com/ipython/ipython/issues/623>`_: Execute DirectView commands while running LoadBalancedView tasks
* `437 <https://github.com/ipython/ipython/issues/437>`_: Users should have autocompletion in the notebook
* `583 <https://github.com/ipython/ipython/issues/583>`_: update manpages
* `594 <https://github.com/ipython/ipython/issues/594>`_: irunner command line options defer to file extensions
* `603 <https://github.com/ipython/ipython/issues/603>`_: Users should see colored text in tracebacks and the pager
* `597 <https://github.com/ipython/ipython/issues/597>`_: UnicodeDecodeError: 'ascii' codec can't decode byte 0xc2
* `608 <https://github.com/ipython/ipython/issues/608>`_: Organize and layout buttons in the notebook panel sections
* `609 <https://github.com/ipython/ipython/issues/609>`_: Implement controls in the Kernel panel section
* `611 <https://github.com/ipython/ipython/issues/611>`_: Add kernel status widget back to notebook
* `610 <https://github.com/ipython/ipython/issues/610>`_: Implement controls in the Cell section panel
* `612 <https://github.com/ipython/ipython/issues/612>`_: Implement Help panel section
* `621 <https://github.com/ipython/ipython/issues/621>`_: [qtconsole] on windows xp, cannot  PageUp more than once
* `616 <https://github.com/ipython/ipython/issues/616>`_: Store exit status of last command
* `605 <https://github.com/ipython/ipython/issues/605>`_: Users should be able to open different notebooks in the cwd
* `302 <https://github.com/ipython/ipython/issues/302>`_: Users should see a consistent behavior in the Out prompt in the html  notebook
* `435 <https://github.com/ipython/ipython/issues/435>`_: Notebook should not import anything by default
* `595 <https://github.com/ipython/ipython/issues/595>`_: qtconsole command issue
* `588 <https://github.com/ipython/ipython/issues/588>`_: ipython-qtconsole uses 100% CPU
* `586 <https://github.com/ipython/ipython/issues/586>`_: ? + plot() Command B0rks QTConsole Strangely
* `585 <https://github.com/ipython/ipython/issues/585>`_: %pdoc throws Errors for classes without __init__ or docstring
* `584 <https://github.com/ipython/ipython/issues/584>`_:  %pdoc throws TypeError
* `580 <https://github.com/ipython/ipython/issues/580>`_: Client instantiation AssertionError
* `569 <https://github.com/ipython/ipython/issues/569>`_: UnicodeDecodeError during startup
* `572 <https://github.com/ipython/ipython/issues/572>`_: Indented command hits error
* `573 <https://github.com/ipython/ipython/issues/573>`_: -wthread breaks indented top-level statements
* `570 <https://github.com/ipython/ipython/issues/570>`_: "--pylab inline" vs. "--pylab=inline"
* `566 <https://github.com/ipython/ipython/issues/566>`_: Can't use exec_file in config file
* `562 <https://github.com/ipython/ipython/issues/562>`_: update docs to reflect '--args=values'
* `558 <https://github.com/ipython/ipython/issues/558>`_: triple quote and %s at beginning of line
* `554 <https://github.com/ipython/ipython/issues/554>`_: Update 0.11 docs to explain Qt console and how to do a clean install
* `553 <https://github.com/ipython/ipython/issues/553>`_: embed() fails if config files not installed
* `8 <https://github.com/ipython/ipython/issues/8>`_: Ensure %gui qt works with new Mayavi and pylab
* `269 <https://github.com/ipython/ipython/issues/269>`_: Provide compatibility api for IPython.Shell().start().mainloop()
* `66 <https://github.com/ipython/ipython/issues/66>`_: Update the main What's New document to reflect work on 0.11
* `549 <https://github.com/ipython/ipython/issues/549>`_: Don't check for 'linux2' value in sys.platform
* `505 <https://github.com/ipython/ipython/issues/505>`_: Qt windows created within imported functions won't show()
* `545 <https://github.com/ipython/ipython/issues/545>`_: qtconsole ignores exec_lines
* `371 <https://github.com/ipython/ipython/issues/371>`_: segfault in qtconsole when kernel quits
* `377 <https://github.com/ipython/ipython/issues/377>`_: Failure: error (nothing to repeat)
* `544 <https://github.com/ipython/ipython/issues/544>`_: Ipython qtconsole pylab config issue.
* `543 <https://github.com/ipython/ipython/issues/543>`_: RuntimeError in completer 
* `440 <https://github.com/ipython/ipython/issues/440>`_: %run filename autocompletion "The kernel heartbeat has been inactive ... " error
* `541 <https://github.com/ipython/ipython/issues/541>`_: log_level is broken in the  ipython Application
* `369 <https://github.com/ipython/ipython/issues/369>`_: windows source install doesn't create scripts correctly
* `351 <https://github.com/ipython/ipython/issues/351>`_: Make sure that the Windows installer handles the top-level IPython scripts.
* `512 <https://github.com/ipython/ipython/issues/512>`_: Two displayhooks in zmq
* `340 <https://github.com/ipython/ipython/issues/340>`_: Make sure that the Windows HPC scheduler support is working for 0.11
* `98 <https://github.com/ipython/ipython/issues/98>`_: Should be able to get help on an object mid-command
* `529 <https://github.com/ipython/ipython/issues/529>`_: unicode problem in qtconsole for windows
* `476 <https://github.com/ipython/ipython/issues/476>`_: Separate input area in Qt Console
* `175 <https://github.com/ipython/ipython/issues/175>`_: Qt console needs configuration support
* `156 <https://github.com/ipython/ipython/issues/156>`_: Key history lost when debugging program crash
* `470 <https://github.com/ipython/ipython/issues/470>`_: decorator: uses deprecated features
* `30 <https://github.com/ipython/ipython/issues/30>`_: readline in OS X does not have correct key bindings
* `503 <https://github.com/ipython/ipython/issues/503>`_: merge IPython.parallel.streamsession and IPython.zmq.session
* `456 <https://github.com/ipython/ipython/issues/456>`_: pathname in document punctuated by dots not slashes
* `451 <https://github.com/ipython/ipython/issues/451>`_: Allow switching the default image format for inline mpl backend
* `79 <https://github.com/ipython/ipython/issues/79>`_: Implement more robust handling of config stages in Application
* `522 <https://github.com/ipython/ipython/issues/522>`_: Encoding problems
* `524 <https://github.com/ipython/ipython/issues/524>`_: otool should not be unconditionally called on osx
* `523 <https://github.com/ipython/ipython/issues/523>`_: Get profile and config file inheritance working
* `519 <https://github.com/ipython/ipython/issues/519>`_: qtconsole --pure: "TypeError: string indices must be integers, not str"
* `516 <https://github.com/ipython/ipython/issues/516>`_: qtconsole --pure: "KeyError: 'ismagic'"
* `520 <https://github.com/ipython/ipython/issues/520>`_: qtconsole --pure: "TypeError: string indices must be integers, not str"
* `450 <https://github.com/ipython/ipython/issues/450>`_: resubmitted tasks sometimes stuck as pending
* `518 <https://github.com/ipython/ipython/issues/518>`_: JSON serialization problems with ObjectId type (MongoDB)
* `178 <https://github.com/ipython/ipython/issues/178>`_: Channels should be named for their function, not their socket type
* `515 <https://github.com/ipython/ipython/issues/515>`_: [ipcluster] termination on os x
* `510 <https://github.com/ipython/ipython/issues/510>`_: qtconsole: indentation problem printing numpy arrays
* `508 <https://github.com/ipython/ipython/issues/508>`_: "AssertionError: Missing message part." in ipython-qtconsole --pure
* `499 <https://github.com/ipython/ipython/issues/499>`_: "ZMQError: Interrupted system call" when saving inline figure
* `426 <https://github.com/ipython/ipython/issues/426>`_: %edit magic fails in qtconsole
* `497 <https://github.com/ipython/ipython/issues/497>`_: Don't show info from .pyd files
* `493 <https://github.com/ipython/ipython/issues/493>`_: QFont::setPointSize: Point size <= 0 (0), must be greater than 0
* `489 <https://github.com/ipython/ipython/issues/489>`_: UnicodeEncodeError in qt.svg.save_svg
* `458 <https://github.com/ipython/ipython/issues/458>`_: embed() doesn't load default config
* `488 <https://github.com/ipython/ipython/issues/488>`_: Using IPython with RubyPython leads to problems with IPython.parallel.client.client.Client.__init()
* `401 <https://github.com/ipython/ipython/issues/401>`_: Race condition when running lbview.apply() fast multiple times in loop
* `168 <https://github.com/ipython/ipython/issues/168>`_: Scrub Launchpad links from code, docs
* `141 <https://github.com/ipython/ipython/issues/141>`_: garbage collection problem (revisited)
* `59 <https://github.com/ipython/ipython/issues/59>`_: test_magic.test_obj_del fails on win32
* `457 <https://github.com/ipython/ipython/issues/457>`_: Backgrounded Tasks not Allowed?  (but easy to slip by . . .)
* `297 <https://github.com/ipython/ipython/issues/297>`_: Shouldn't use pexpect for subprocesses in in-process terminal frontend
* `110 <https://github.com/ipython/ipython/issues/110>`_: magic to return exit status
* `473 <https://github.com/ipython/ipython/issues/473>`_: OSX readline detection fails in the debugger
* `466 <https://github.com/ipython/ipython/issues/466>`_: tests fail without unicode filename support
* `468 <https://github.com/ipython/ipython/issues/468>`_: iptest script has 0 exit code even when tests fail
* `465 <https://github.com/ipython/ipython/issues/465>`_: client.db_query() behaves different with SQLite and MongoDB
* `467 <https://github.com/ipython/ipython/issues/467>`_: magic_install_default_config test fails when there is no .ipython directory
* `463 <https://github.com/ipython/ipython/issues/463>`_: IPYTHON_DIR (and IPYTHONDIR) don't expand tilde to '~' directory
* `446 <https://github.com/ipython/ipython/issues/446>`_: Test machinery is imported at normal runtime
* `438 <https://github.com/ipython/ipython/issues/438>`_: Users should be able to use Up/Down for cell navigation
* `439 <https://github.com/ipython/ipython/issues/439>`_: Users should be able to copy notebook input and output
* `291 <https://github.com/ipython/ipython/issues/291>`_: Rename special display methods and put them lower in priority than display functions
* `447 <https://github.com/ipython/ipython/issues/447>`_: Instantiating classes without __init__ function causes kernel to crash
* `444 <https://github.com/ipython/ipython/issues/444>`_: Ctrl + t in WxIPython Causes Unexpected Behavior
* `445 <https://github.com/ipython/ipython/issues/445>`_: qt and console Based Startup Errors
* `428 <https://github.com/ipython/ipython/issues/428>`_: ipcluster doesn't handle stale pid info well
* `434 <https://github.com/ipython/ipython/issues/434>`_: 10.0.2 seg fault with rpy2
* `441 <https://github.com/ipython/ipython/issues/441>`_: Allow running a block of code in a file
* `432 <https://github.com/ipython/ipython/issues/432>`_: Silent request fails
* `409 <https://github.com/ipython/ipython/issues/409>`_: Test failure in IPython.lib
* `402 <https://github.com/ipython/ipython/issues/402>`_: History section of messaging spec is incorrect
* `88 <https://github.com/ipython/ipython/issues/88>`_: Error when inputting UTF8 CJK characters
* `366 <https://github.com/ipython/ipython/issues/366>`_: Ctrl-K should kill line and store it, so that Ctrl-y can yank it back
* `425 <https://github.com/ipython/ipython/issues/425>`_: typo in %gui magic help
* `304 <https://github.com/ipython/ipython/issues/304>`_: Persistent warnings if old configuration files exist
* `216 <https://github.com/ipython/ipython/issues/216>`_: crash of ipython when alias is used with %s and echo
* `412 <https://github.com/ipython/ipython/issues/412>`_: add support to automatic retry of tasks
* `411 <https://github.com/ipython/ipython/issues/411>`_: add support to continue tasks
* `417 <https://github.com/ipython/ipython/issues/417>`_: IPython should display things unsorted if it can't sort them
* `416 <https://github.com/ipython/ipython/issues/416>`_: wrong encode when printing unicode string
* `376 <https://github.com/ipython/ipython/issues/376>`_: Failing InputsplitterTest
* `405 <https://github.com/ipython/ipython/issues/405>`_: TraitError in traitlets.py(332) on any input
* `392 <https://github.com/ipython/ipython/issues/392>`_: UnicodeEncodeError on start
* `137 <https://github.com/ipython/ipython/issues/137>`_: sys.getfilesystemencoding return value not checked
* `300 <https://github.com/ipython/ipython/issues/300>`_: Users should be able to manage kernels and kernel sessions from the notebook UI
* `301 <https://github.com/ipython/ipython/issues/301>`_: Users should have access to working Kernel, Tabs, Edit, Help menus in the notebook
* `396 <https://github.com/ipython/ipython/issues/396>`_: cursor move triggers a lot of IO access
* `379 <https://github.com/ipython/ipython/issues/379>`_: Minor doc nit: --paging argument
* `399 <https://github.com/ipython/ipython/issues/399>`_: Add task queue limit in engine when load-balancing
* `78 <https://github.com/ipython/ipython/issues/78>`_: StringTask won't take unicode code strings
* `391 <https://github.com/ipython/ipython/issues/391>`_: MongoDB.add_record() does not work in 0.11dev
* `365 <https://github.com/ipython/ipython/issues/365>`_: newparallel on Windows
* `386 <https://github.com/ipython/ipython/issues/386>`_: FAIL: test that pushed functions have access to globals
* `387 <https://github.com/ipython/ipython/issues/387>`_: Interactively defined functions can't access user namespace
* `118 <https://github.com/ipython/ipython/issues/118>`_: Snow Leopard ipy_vimserver POLL error
* `394 <https://github.com/ipython/ipython/issues/394>`_: System escape interpreted in multi-line string
* `26 <https://github.com/ipython/ipython/issues/26>`_: find_job_cmd is too hasty to fail on Windows
* `368 <https://github.com/ipython/ipython/issues/368>`_: Installation instructions in dev docs are completely wrong
* `380 <https://github.com/ipython/ipython/issues/380>`_: qtconsole pager RST - HTML not happening consistently
* `367 <https://github.com/ipython/ipython/issues/367>`_: Qt console doesn't support ibus input method
* `375 <https://github.com/ipython/ipython/issues/375>`_: Missing libraries cause ImportError in tests
* `71 <https://github.com/ipython/ipython/issues/71>`_: temp file errors in iptest IPython.core
* `350 <https://github.com/ipython/ipython/issues/350>`_: Decide how to handle displayhook being triggered multiple times
* `360 <https://github.com/ipython/ipython/issues/360>`_: Remove `runlines` method
* `125 <https://github.com/ipython/ipython/issues/125>`_: Exec lines in config should not contribute to line numbering or history
* `20 <https://github.com/ipython/ipython/issues/20>`_: Robust readline support on OS X's builtin Python
* `147 <https://github.com/ipython/ipython/issues/147>`_: On Windows, %page is being too restrictive to split line by \r\n only
* `326 <https://github.com/ipython/ipython/issues/326>`_: Update docs and examples for parallel stuff to reflect movement away from Twisted
* `341 <https://github.com/ipython/ipython/issues/341>`_: FIx Parallel Magics for newparallel
* `338 <https://github.com/ipython/ipython/issues/338>`_: Usability improvements to Qt console
* `142 <https://github.com/ipython/ipython/issues/142>`_: unexpected auto-indenting when variables names that start with 'pass' 
* `296 <https://github.com/ipython/ipython/issues/296>`_: Automatic PDB via %pdb doesn't work
* `337 <https://github.com/ipython/ipython/issues/337>`_: exit( and quit( in Qt console produces phantom signature/docstring popup, even though quit() or exit() raises NameError
* `318 <https://github.com/ipython/ipython/issues/318>`_: %debug broken in master: invokes missing save_history() method
* `307 <https://github.com/ipython/ipython/issues/307>`_: lines ending with semicolon should not go to cache
* `104 <https://github.com/ipython/ipython/issues/104>`_: have ipengine run start-up scripts before registering with the controller
* `33 <https://github.com/ipython/ipython/issues/33>`_: The skip_doctest decorator is failing to work on Shell.MatplotlibShellBase.magic_run
* `336 <https://github.com/ipython/ipython/issues/336>`_: Missing figure development/figs/iopubfade.png for docs
* `49 <https://github.com/ipython/ipython/issues/49>`_: %clear should also delete _NN references and Out[NN] ones
* `335 <https://github.com/ipython/ipython/issues/335>`_: using setuptools installs every script twice
* `306 <https://github.com/ipython/ipython/issues/306>`_: multiline strings at end of input cause noop
* `327 <https://github.com/ipython/ipython/issues/327>`_: PyPy compatibility
* `328 <https://github.com/ipython/ipython/issues/328>`_: %run script.ipy raises "ERROR! Session/line number was not unique in database."
* `7 <https://github.com/ipython/ipython/issues/7>`_: Update the changes doc to reflect the kernel config work
* `303 <https://github.com/ipython/ipython/issues/303>`_: Users should be able to scroll a notebook w/o moving the menu/buttons
* `322 <https://github.com/ipython/ipython/issues/322>`_: Embedding an interactive IPython shell 
* `321 <https://github.com/ipython/ipython/issues/321>`_: %debug broken in master
* `287 <https://github.com/ipython/ipython/issues/287>`_: Crash when using %macros in sqlite-history branch
* `55 <https://github.com/ipython/ipython/issues/55>`_: Can't edit files whose names begin with numbers
* `284 <https://github.com/ipython/ipython/issues/284>`_: In variable no longer works in 0.11
* `92 <https://github.com/ipython/ipython/issues/92>`_: Using multiprocessing module crashes parallel IPython
* `262 <https://github.com/ipython/ipython/issues/262>`_: Fail to recover history after force-kill.
* `320 <https://github.com/ipython/ipython/issues/320>`_: Tab completing re.search objects crashes IPython
* `317 <https://github.com/ipython/ipython/issues/317>`_: IPython.kernel: parallel map issues
* `197 <https://github.com/ipython/ipython/issues/197>`_: ipython-qtconsole unicode problem in magic ls
* `305 <https://github.com/ipython/ipython/issues/305>`_: more readline shortcuts in qtconsole
* `314 <https://github.com/ipython/ipython/issues/314>`_: Multi-line, multi-block cells can't be executed.
* `308 <https://github.com/ipython/ipython/issues/308>`_: Test suite should set sqlite history to work in :memory:
* `202 <https://github.com/ipython/ipython/issues/202>`_: Matplotlib native 'MacOSX' backend broken in '-pylab' mode
* `196 <https://github.com/ipython/ipython/issues/196>`_: IPython can't deal with unicode file name.
* `25 <https://github.com/ipython/ipython/issues/25>`_: unicode bug - encoding input
* `290 <https://github.com/ipython/ipython/issues/290>`_: try/except/else clauses can't be typed, code input stops too early.
* `43 <https://github.com/ipython/ipython/issues/43>`_: Implement SSH support in ipcluster
* `6 <https://github.com/ipython/ipython/issues/6>`_: Update the Sphinx docs for the new ipcluster
* `9 <https://github.com/ipython/ipython/issues/9>`_: Getting "DeadReferenceError: Calling Stale Broker" after ipcontroller restart
* `132 <https://github.com/ipython/ipython/issues/132>`_: Ipython prevent south from working
* `27 <https://github.com/ipython/ipython/issues/27>`_: generics.complete_object broken
* `60 <https://github.com/ipython/ipython/issues/60>`_: Improve absolute import management for iptest.py
* `31 <https://github.com/ipython/ipython/issues/31>`_: Issues in magic_whos code
* `52 <https://github.com/ipython/ipython/issues/52>`_: Document testing process better
* `44 <https://github.com/ipython/ipython/issues/44>`_: Merge history from multiple sessions
* `182 <https://github.com/ipython/ipython/issues/182>`_: ipython q4thread in version 10.1 not starting properly
* `143 <https://github.com/ipython/ipython/issues/143>`_: Ipython.gui.wx.ipython_view.IPShellWidget: ignores user*_ns arguments
* `127 <https://github.com/ipython/ipython/issues/127>`_: %edit does not work on filenames consisted of pure numbers
* `126 <https://github.com/ipython/ipython/issues/126>`_: Can't transfer command line argument to script
* `28 <https://github.com/ipython/ipython/issues/28>`_: Offer finer control for initialization of input streams
* `58 <https://github.com/ipython/ipython/issues/58>`_: ipython change char '0xe9' to 4 spaces
* `68 <https://github.com/ipython/ipython/issues/68>`_: Problems with Control-C stopping ipcluster on Windows/Python2.6
* `24 <https://github.com/ipython/ipython/issues/24>`_: ipcluster does not start all the engines
* `240 <https://github.com/ipython/ipython/issues/240>`_: Incorrect method displayed in %psource
* `120 <https://github.com/ipython/ipython/issues/120>`_: inspect.getsource fails for functions defined on command line
* `212 <https://github.com/ipython/ipython/issues/212>`_: IPython ignores exceptions in the first evaulation of class attrs
* `108 <https://github.com/ipython/ipython/issues/108>`_: ipython disables python logger
* `100 <https://github.com/ipython/ipython/issues/100>`_: Overzealous introspection
* `18 <https://github.com/ipython/ipython/issues/18>`_: %cpaste freeze sync frontend
* `200 <https://github.com/ipython/ipython/issues/200>`_: Unicode error when starting ipython in a folder with non-ascii path
* `130 <https://github.com/ipython/ipython/issues/130>`_: Deadlock when importing a module that creates an IPython client
* `134 <https://github.com/ipython/ipython/issues/134>`_: multline block scrolling
* `46 <https://github.com/ipython/ipython/issues/46>`_: Input to %timeit is not preparsed
* `285 <https://github.com/ipython/ipython/issues/285>`_: ipcluster local -n 4 fails
* `205 <https://github.com/ipython/ipython/issues/205>`_: In the Qt console, Tab should insert 4 spaces when not completing
* `145 <https://github.com/ipython/ipython/issues/145>`_: Bug on MSW systems: idle can not be set as default IPython editor. Fix Suggested.
* `77 <https://github.com/ipython/ipython/issues/77>`_: ipython oops in cygwin
* `121 <https://github.com/ipython/ipython/issues/121>`_: If plot windows are closed via window controls, no more plotting is possible.
* `111 <https://github.com/ipython/ipython/issues/111>`_: Iterator version of TaskClient.map() that returns results as they become available
* `109 <https://github.com/ipython/ipython/issues/109>`_: WinHPCLauncher is a hard dependency that causes errors in the test suite
* `86 <https://github.com/ipython/ipython/issues/86>`_: Make IPython work with multiprocessing
* `15 <https://github.com/ipython/ipython/issues/15>`_: Implement SGE support in ipcluster
* `3 <https://github.com/ipython/ipython/issues/3>`_: Implement PBS support in ipcluster
* `53 <https://github.com/ipython/ipython/issues/53>`_: Internal Python error in the inspect module
* `74 <https://github.com/ipython/ipython/issues/74>`_: Manager() [from multiprocessing module] hangs ipythonx but not ipython
* `51 <https://github.com/ipython/ipython/issues/51>`_: Out not working with ipythonx
* `201 <https://github.com/ipython/ipython/issues/201>`_: use session.send throughout zmq code
* `115 <https://github.com/ipython/ipython/issues/115>`_: multiline specials not defined in 0.11 branch
* `93 <https://github.com/ipython/ipython/issues/93>`_: when looping, cursor appears at leftmost point in newline
* `133 <https://github.com/ipython/ipython/issues/133>`_: whitespace after Source introspection
* `50 <https://github.com/ipython/ipython/issues/50>`_: Ctrl-C with -gthread on Windows, causes uncaught IOError
* `65 <https://github.com/ipython/ipython/issues/65>`_: Do not use .message attributes in exceptions, deprecated in 2.6
* `76 <https://github.com/ipython/ipython/issues/76>`_: syntax error when raise is inside except process
* `107 <https://github.com/ipython/ipython/issues/107>`_: bdist_rpm causes traceback looking for a non-existant file
* `113 <https://github.com/ipython/ipython/issues/113>`_: initial magic ? (question mark) fails before wildcard
* `128 <https://github.com/ipython/ipython/issues/128>`_: Pdb instance has no attribute 'curframe'
* `139 <https://github.com/ipython/ipython/issues/139>`_: running with -pylab pollutes namespace
* `140 <https://github.com/ipython/ipython/issues/140>`_: malloc error during tab completion of numpy array member functions starting with 'c'
* `153 <https://github.com/ipython/ipython/issues/153>`_: ipy_vimserver traceback on Windows
* `154 <https://github.com/ipython/ipython/issues/154>`_: using ipython in Slicer3 show how os.environ['HOME'] is not defined
* `185 <https://github.com/ipython/ipython/issues/185>`_: show() blocks in pylab mode with ipython 0.10.1 
* `189 <https://github.com/ipython/ipython/issues/189>`_: Crash on tab completion
* `274 <https://github.com/ipython/ipython/issues/274>`_: bashism in sshx.sh
* `276 <https://github.com/ipython/ipython/issues/276>`_: Calling `sip.setapi` does not work if app has already imported from PyQt4
* `277 <https://github.com/ipython/ipython/issues/277>`_: matplotlib.image imgshow from 10.1 segfault
* `288 <https://github.com/ipython/ipython/issues/288>`_: Incorrect docstring in zmq/kernelmanager.py
* `286 <https://github.com/ipython/ipython/issues/286>`_: Fix IPython.Shell compatibility layer
* `99 <https://github.com/ipython/ipython/issues/99>`_: blank lines in history
* `129 <https://github.com/ipython/ipython/issues/129>`_: psearch: TypeError: expected string or buffer
* `190 <https://github.com/ipython/ipython/issues/190>`_: Add option to format float point output
* `246 <https://github.com/ipython/ipython/issues/246>`_: Application not conforms XDG Base Directory Specification
* `48 <https://github.com/ipython/ipython/issues/48>`_: IPython should follow the XDG Base Directory spec for configuration
* `176 <https://github.com/ipython/ipython/issues/176>`_: Make client-side history persistence readline-independent
* `279 <https://github.com/ipython/ipython/issues/279>`_: Backtraces when using ipdb do not respect -colour LightBG setting
* `119 <https://github.com/ipython/ipython/issues/119>`_: Broken type filter in magic_who_ls
* `271 <https://github.com/ipython/ipython/issues/271>`_: Intermittent problem with print output in Qt console.
* `270 <https://github.com/ipython/ipython/issues/270>`_: Small typo in IPython developer’s guide
* `166 <https://github.com/ipython/ipython/issues/166>`_: Add keyboard accelerators to Qt close dialog
* `173 <https://github.com/ipython/ipython/issues/173>`_: asymmetrical ctrl-A/ctrl-E behavior in multiline
* `45 <https://github.com/ipython/ipython/issues/45>`_: Autosave history for robustness
* `162 <https://github.com/ipython/ipython/issues/162>`_: make command history persist in ipythonqt
* `161 <https://github.com/ipython/ipython/issues/161>`_: make ipythonqt exit without dialog when exit() is called
* `263 <https://github.com/ipython/ipython/issues/263>`_: [ipython + numpy] Some test errors 
* `256 <https://github.com/ipython/ipython/issues/256>`_: reset docstring ipython 0.10 
* `258 <https://github.com/ipython/ipython/issues/258>`_: allow caching to avoid matplotlib object references
* `248 <https://github.com/ipython/ipython/issues/248>`_: Can't open and read files after upgrade from 0.10 to 0.10.0
* `247 <https://github.com/ipython/ipython/issues/247>`_: ipython + Stackless
* `245 <https://github.com/ipython/ipython/issues/245>`_: Magic save and macro missing newlines, line ranges don't match prompt numbers.
* `241 <https://github.com/ipython/ipython/issues/241>`_: "exit" hangs on terminal version of IPython
* `213 <https://github.com/ipython/ipython/issues/213>`_: ipython -pylab no longer plots interactively on 0.10.1
* `4 <https://github.com/ipython/ipython/issues/4>`_: wx frontend don't display well commands output
* `5 <https://github.com/ipython/ipython/issues/5>`_: ls command not supported in ipythonx wx frontend
* `1 <https://github.com/ipython/ipython/issues/1>`_: Document winhpcjob.py and launcher.py
* `83 <https://github.com/ipython/ipython/issues/83>`_: Usage of testing.util.DeferredTestCase should be replace with twisted.trial.unittest.TestCase
* `117 <https://github.com/ipython/ipython/issues/117>`_: Redesign how Component instances are tracked and queried
* `47 <https://github.com/ipython/ipython/issues/47>`_: IPython.kernel.client cannot be imported inside an engine
* `105 <https://github.com/ipython/ipython/issues/105>`_: Refactor the task dependencies system
* `210 <https://github.com/ipython/ipython/issues/210>`_: 0.10.1 doc mistake - New IPython Sphinx directive error
* `209 <https://github.com/ipython/ipython/issues/209>`_: can't activate IPython parallel magics
* `206 <https://github.com/ipython/ipython/issues/206>`_: Buggy linewrap in Mac OSX Terminal
* `194 <https://github.com/ipython/ipython/issues/194>`_: !sudo <command> displays password in plain text
* `186 <https://github.com/ipython/ipython/issues/186>`_: %edit issue under OS X 10.5 - IPython 0.10.1
* `11 <https://github.com/ipython/ipython/issues/11>`_: Create a daily build PPA for ipython
* `144 <https://github.com/ipython/ipython/issues/144>`_: logo missing from sphinx docs
* `181 <https://github.com/ipython/ipython/issues/181>`_: cls command does not work on windows
* `169 <https://github.com/ipython/ipython/issues/169>`_: Kernel can only be bound to localhost
* `36 <https://github.com/ipython/ipython/issues/36>`_: tab completion does not escape ()
* `177 <https://github.com/ipython/ipython/issues/177>`_: Report tracebacks of interactively entered input
* `148 <https://github.com/ipython/ipython/issues/148>`_: dictionary having multiple keys having frozenset fails to print on IPython
* `160 <https://github.com/ipython/ipython/issues/160>`_: magic_gui throws TypeError when gui magic is used
* `150 <https://github.com/ipython/ipython/issues/150>`_: History entries ending with parentheses corrupt command line on OS X 10.6.4
* `146 <https://github.com/ipython/ipython/issues/146>`_: -ipythondir - using an alternative .ipython dir for rc type stuff
* `114 <https://github.com/ipython/ipython/issues/114>`_: Interactive strings get mangled with "_ip.magic"
* `135 <https://github.com/ipython/ipython/issues/135>`_: crash on  invalid print
* `69 <https://github.com/ipython/ipython/issues/69>`_: Usage of "mycluster" profile in docs and examples
* `37 <https://github.com/ipython/ipython/issues/37>`_: Fix colors in output of ResultList on Windows
.. _issues_list_200:

Issues closed in the 2.x development cycle
==========================================

Issues closed in 2.4.1
----------------------

GitHub stats for 2014/11/01 - 2015/01/30

.. note::

    IPython 2.4.0 was released without a few of the backports listed below.
    2.4.1 has the correct patches intended for 2.4.0.

These lists are automatically generated, and may be incomplete or contain duplicates.

The following 7 authors contributed 35 commits.

* Benjamin Ragan-Kelley
* Carlos Cordoba
* Damon Allen
* Jessica B. Hamrick
* Mateusz Paprocki
* Peter Würtz
* Thomas Kluyver

We closed 10 issues and merged 6 pull requests;
this is the full list (generated with the script 
:file:`tools/github_stats.py`):

Pull Requests (10):

* :ghpull:`7106`: Changed the display order of rich output in the live notebook.
* :ghpull:`6878`: Update pygments monkeypatch for compatibility with Pygments 2.0
* :ghpull:`6778`: backport nbformat v4 to 2.x
* :ghpull:`6761`: object_info_reply field is oname, not name
* :ghpull:`6653`: Fix IPython.utils.ansispan() to ignore stray [0m
* :ghpull:`6706`: Correctly display prompt numbers that are ``None``
* :ghpull:`6634`: don't use contains in SelectWidget item_query
* :ghpull:`6593`: note how to start the qtconsole
* :ghpull:`6281`: more minor fixes to release scripts
* :ghpull:`5458`: Add support for PyQt5.

Issues (6):

* :ghissue:`7272`: qtconsole problems with pygments
* :ghissue:`7049`: Cause TypeError: 'NoneType' object is not callable in qtconsole 
* :ghissue:`6877`: Qt console doesn't work with pygments 2.0rc1
* :ghissue:`6689`: Problem with string containing two or more question marks
* :ghissue:`6702`: Cell numbering after ``ClearOutput`` preprocessor
* :ghissue:`6633`: selectwidget doesn't display 1 as a selection choice when passed in as a member of values list


Issues closed in 2.3.1
----------------------

Just one bugfix: fixed bad CRCRLF line-endings in notebooks on Windows

Pull Requests (1):

* :ghpull:`6911`: don't use text mode in mkstemp

Issues (1):

* :ghissue:`6599`: Notebook.ipynb CR+LF turned into CR+CR+LF


Issues closed in 2.3.0
----------------------

GitHub stats for 2014/08/06 - 2014/10/01

These lists are automatically generated, and may be incomplete or contain duplicates.

The following 6 authors contributed 31 commits.

* Benjamin Ragan-Kelley
* David Hirschfeld
* Eric Firing
* Jessica B. Hamrick
* Matthias Bussonnier
* Thomas Kluyver

We closed 16 issues and merged 9 pull requests;
this is the full list (generated with the script 
:file:`tools/github_stats.py`):

Pull Requests (16):

* :ghpull:`6587`: support ``%matplotlib qt5`` and ``%matplotlib nbagg``
* :ghpull:`6583`: Windows symlink test fixes
* :ghpull:`6585`: fixes :ghissue:`6473`
* :ghpull:`6581`: Properly mock winreg functions for test
* :ghpull:`6556`: Use some more informative asserts in inprocess kernel tests
* :ghpull:`6514`: Fix for copying metadata flags
* :ghpull:`6453`: Copy file metadata in atomic save
* :ghpull:`6480`: only compare host:port in Websocket.check_origin
* :ghpull:`6483`: Trim anchor link in heading cells, fixes :ghissue:`6324`
* :ghpull:`6410`: Fix relative import in appnope
* :ghpull:`6395`: update mathjax CDN url in nbconvert template
* :ghpull:`6269`: Implement atomic save
* :ghpull:`6374`: Rename ``abort_queues`` --> ``_abort_queues``
* :ghpull:`6321`: Use appnope in qt and wx gui support from the terminal; closes :ghissue:`6189`
* :ghpull:`6318`: use write_error instead of get_error_html
* :ghpull:`6303`: Fix error message when failing to load a notebook

Issues (9):

* :ghissue:`6057`: ``%matplotlib`` + qt5
* :ghissue:`6518`: Test failure in atomic save on Windows
* :ghissue:`6473`: Switching between "Raw Cell Format" and "Edit Metadata" does not work
* :ghissue:`6405`: Creating a notebook should respect directory permissions; saving should respect prior permissions
* :ghissue:`6324`: Anchors in Heading don't work.
* :ghissue:`6409`: No module named '_dummy'
* :ghissue:`6392`: Mathjax library link broken
* :ghissue:`6329`: IPython Notebook Server URL now requires "tree" at the end of the URL? (version 2.2)
* :ghissue:`6189`: ipython console freezes for increasing no of seconds in %pylab mode

Issues closed in 2.2.0
----------------------

GitHub stats for 2014/05/21 - 2014/08/06 (tag: rel-2.1.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

The following 13 authors contributed 36 commits.

* Adam Hodgen
* Benjamin Ragan-Kelley
* Björn Grüning
* Dara Adib
* Eric Galloway
* Jonathan Frederic
* Kyle Kelley
* Matthias Bussonnier
* Paul Ivanov
* Shayne Hodge
* Steven Anton
* Thomas Kluyver
* Zahari

We closed 23 issues and merged 11 pull requests;
this is the full list (generated with the script 
:file:`tools/github_stats.py`):

Pull Requests (23):

* :ghpull:`6279`: minor updates to release scripts
* :ghpull:`6273`: Upgrade default mathjax version.
* :ghpull:`6249`: always use HTTPS getting mathjax from CDN
* :ghpull:`6114`: update hmac signature comparison
* :ghpull:`6195`: Close handle on new temporary files before returning filename
* :ghpull:`6143`: pin tornado to < 4 on travis js tests
* :ghpull:`6134`: remove rackcdn https workaround for mathjax cdn
* :ghpull:`6120`: Only allow iframe embedding on same origin.
* :ghpull:`6117`: Remove / from route of TreeRedirectHandler.
* :ghpull:`6105`: only set allow_origin_pat if defined
* :ghpull:`6102`: Add newline if missing to end of script magic cell
* :ghpull:`6077`: allow unicode keys in dicts in json_clean
* :ghpull:`6061`: make CORS configurable
* :ghpull:`6081`: don’t modify dict keys while iterating through them
* :ghpull:`5803`: unify visual line handling
* :ghpull:`6005`: Changed right arrow key movement function to mirror left arrow key
* :ghpull:`6029`: add pickleutil.PICKLE_PROTOCOL
* :ghpull:`6003`: Set kernel_id before checking websocket
* :ghpull:`5994`: Fix ssh tunnel for Python3
* :ghpull:`5973`: Do not create checkpoint_dir relative to current dir
* :ghpull:`5933`: fix qt_loader import hook signature
* :ghpull:`5944`: Markdown rendering bug fix.
* :ghpull:`5917`: use shutil.move instead of os.rename

Issues (11):

* :ghissue:`6246`: Include MathJax by default or access the CDN over a secure connection
* :ghissue:`5525`: Websocket origin check fails when used with Apache WS proxy
* :ghissue:`5901`: 2 test failures in Python 3.4 in parallel group
* :ghissue:`5926`: QT console: text selection cannot be made from left to right with keyboard
* :ghissue:`5998`: use_dill does not work in Python 3.4
* :ghissue:`5964`: Traceback on Qt console exit
* :ghissue:`5787`: Error in Notebook-Generated latex (nbconvert)
* :ghissue:`5950`: qtconsole truncates help
* :ghissue:`5943`: 2.x: notebook fails to load when using HTML comments
* :ghissue:`5932`: Qt ImportDenier Does Not Adhere to PEP302
* :ghissue:`5898`: OSError when moving configuration file

Issues closed in 2.1.0
----------------------

GitHub stats for 2014/04/02 - 2014/05/21 (since 2.0.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

The following 35 authors contributed 145 commits.

* Adrian Price-Whelan
* Aron Ahmadia
* Benjamin Ragan-Kelley
* Benjamin Schultz
* Björn Linse
* Blake Griffith
* chebee7i
* Damián Avila
* Dav Clark
* dexterdev
* Erik Tollerud
* Grzegorz Rożniecki
* Jakob Gager
* jdavidheiser
* Jessica B. Hamrick
* Jim Garrison
* Jonathan Frederic
* Matthias Bussonnier
* Maximilian Albert
* Mohan Raj Rajamanickam
* ncornette
* Nikolay Koldunov
* Nile Geisinger
* Pankaj Pandey
* Paul Ivanov
* Pierre Haessig
* Raffaele De Feo
* Renaud Richardet
* Spencer Nelson
* Steve Chan
* sunny
* Susan Tan
* Thomas Kluyver
* Yaroslav Halchenko
* zah

We closed a total of 129 issues, 92 pull requests and 37 regular issues;
this is the full list (generated with the script 
:file:`tools/github_stats.py --milestone 2.1`):

Pull Requests (92):

* :ghpull:`5871`: specify encoding in msgpack.unpackb
* :ghpull:`5869`: Catch more errors from clipboard access on Windows
* :ghpull:`5866`: Make test robust against differences in line endings
* :ghpull:`5605`: Two cell toolbar fixes.
* :ghpull:`5843`: remove Firefox-specific CSS workaround
* :ghpull:`5845`: Pass Windows interrupt event to kernels as an environment variable
* :ghpull:`5835`: fix typo in v2 convert
* :ghpull:`5841`: Fix writing history with output to a file in Python 2
* :ghpull:`5842`: fix typo in nbconvert help
* :ghpull:`5846`: Fix typos in Cython example
* :ghpull:`5839`: Close graphics dev in finally clause
* :ghpull:`5837`: pass on install docs
* :ghpull:`5832`: Fixed example to work with python3
* :ghpull:`5826`: allow notebook tour instantiation to fail
* :ghpull:`5560`: Minor expansion of Cython example
* :ghpull:`5818`: interpret any exception in getcallargs as not callable
* :ghpull:`5816`: Add output to IPython directive when in verbatim mode.
* :ghpull:`5822`: Don't overwrite widget description in interact
* :ghpull:`5782`: Silence exception thrown by completer when dir() does not return a list
* :ghpull:`5807`: Drop log level to info for Qt console shutdown
* :ghpull:`5814`: Remove -i options from mv, rm and cp aliases
* :ghpull:`5812`: Fix application name when printing subcommand help.
* :ghpull:`5804`: remove an inappropriate ``!``
* :ghpull:`5805`: fix engine startup files
* :ghpull:`5806`: Don't auto-move .config/ipython if symbolic link
* :ghpull:`5716`: Add booktabs package to latex base.tplx
* :ghpull:`5669`: allows threadsafe sys.stdout.flush from background threads
* :ghpull:`5668`: allow async output on the most recent request
* :ghpull:`5768`: fix cursor keys in long lines wrapped in markdown
* :ghpull:`5788`: run cells with ``silent=True`` in ``%run nb.ipynb``
* :ghpull:`5715`: log all failed ajax API requests
* :ghpull:`5769`: Don't urlescape the text that goes into a title tag
* :ghpull:`5762`: Fix check for pickling closures
* :ghpull:`5766`: View.map with empty sequence should return empty list
* :ghpull:`5758`: Applied bug fix: using fc and ec did not properly set the figure canvas ...
* :ghpull:`5754`: Format command name into subcommand_description at run time, not import
* :ghpull:`5744`: Describe using PyPI/pip to distribute & install extensions
* :ghpull:`5712`: monkeypatch inspect.findsource only when we use it
* :ghpull:`5708`: create checkpoints dir in notebook subdirectories
* :ghpull:`5714`: log error message when API requests fail
* :ghpull:`5732`: Quick typo fix in nbformat/convert.py
* :ghpull:`5713`: Fix a NameError in IPython.parallel
* :ghpull:`5704`: Update nbconvertapp.py
* :ghpull:`5534`: cleanup some ``pre`` css inheritance
* :ghpull:`5699`: don't use common names in require decorators
* :ghpull:`5692`: Update notebook.rst fixing broken reference to notebook examples readme
* :ghpull:`5693`: Update parallel_intro.rst to fix a broken link to examples
* :ghpull:`5486`: disambiguate to location when no IPs can be determined
* :ghpull:`5574`: Remove the outdated keyboard shortcuts from notebook docs
* :ghpull:`5568`: Use ``__qualname__`` in pretty reprs for Python 3
* :ghpull:`5678`: Fix copy & paste error in docstring of ImageWidget class
* :ghpull:`5677`: Fix %bookmark -l for Python 3
* :ghpull:`5670`: nbconvert: Fix CWD imports
* :ghpull:`5647`: Mention git hooks in install documentation
* :ghpull:`5671`: Fix blank slides issue in Reveal slideshow pdf export
* :ghpull:`5657`: use 'localhost' as default for the notebook server
* :ghpull:`5584`: more semantic icons
* :ghpull:`5594`: update components with marked-0.3.2
* :ghpull:`5500`: check for Python 3.2
* :ghpull:`5582`: reset readline after running PYTHONSTARTUP
* :ghpull:`5630`: Fixed Issue :ghissue:`4012` Added Help menubar link to Github markdown doc
* :ghpull:`5613`: Fixing bug :ghissue:`5607`
* :ghpull:`5633`: Provide more help if lessc is not found.
* :ghpull:`5620`: fixed a typo in IPython.core.formatters
* :ghpull:`5619`: Fix typo in storemagic module docstring
* :ghpull:`5592`: add missing ``browser`` to notebook_aliases list
* :ghpull:`5506`: Fix ipconfig regex pattern
* :ghpull:`5581`: Fix rmagic for cells ending in comment.
* :ghpull:`5576`: only process cr if it's found
* :ghpull:`5478`: Add git-hooks install script. Update README.md
* :ghpull:`5546`: do not shutdown notebook if 'n' is part of answer
* :ghpull:`5527`: Don't remove upload items from nav tree unless explicitly requested.
* :ghpull:`5501`: remove inappropriate wheel tag override
* :ghpull:`5548`: FileNotebookManager: Use shutil.move() instead of os.rename()
* :ghpull:`5524`: never use ``for (var i in array)``
* :ghpull:`5459`: Fix interact animation page jump FF
* :ghpull:`5559`: Minor typo fix in "Cython Magics.ipynb"
* :ghpull:`5507`: Fix typo in interactive widgets examples index notebook
* :ghpull:`5554`: Make HasTraits pickleable
* :ghpull:`5535`: fix n^2 performance issue in coalesce_streams preprocessor
* :ghpull:`5522`: fix iteration over Client
* :ghpull:`5488`: Added missing require and jquery from cdn.
* :ghpull:`5516`: ENH: list generated config files in generated, and rm them upon clean
* :ghpull:`5493`: made a minor fix to one of the widget examples
* :ghpull:`5512`: Update tooltips to refer to shift-tab
* :ghpull:`5505`: Make backport_pr work on Python 3
* :ghpull:`5503`: check explicitly for 'dev' before adding the note to docs
* :ghpull:`5498`: use milestones to indicate backport
* :ghpull:`5492`: Polish whatsnew docs
* :ghpull:`5495`: Fix various broken things in docs
* :ghpull:`5496`: Exclude whatsnew/pr directory from docs builds
* :ghpull:`5489`: Fix required Python versions

Issues (37):

* :ghissue:`5364`: Horizontal scrollbar hides cell's last line on Firefox
* :ghissue:`5192`: horisontal scrollbar overlaps output or touches next cell
* :ghissue:`5840`: Third-party Windows kernels don't get interrupt signal
* :ghissue:`2412`: print history to file using qtconsole and notebook
* :ghissue:`5703`: Notebook doesn't render with "ask me every time" cookie setting in Firefox
* :ghissue:`5817`: calling mock object in IPython 2.0.0 under Python 3.4.0 raises AttributeError
* :ghissue:`5499`: Error running widgets nbconvert example
* :ghissue:`5654`: Broken links from ipython documentation
* :ghissue:`5019`: print in QT event callback doesn't show up in ipython notebook.
* :ghissue:`5800`: Only last In prompt number set ?
* :ghissue:`5801`: startup_command specified in ipengine_config.py is not executed
* :ghissue:`5690`: ipython 2.0.0 and pandoc 1.12.2.1 problem
* :ghissue:`5408`: Add checking/flushing of background output from kernel in mainloop
* :ghissue:`5407`: clearing message handlers on status=idle loses async output
* :ghissue:`5467`: Incorrect behavior of up/down keyboard arrows in code cells on wrapped lines
* :ghissue:`3085`: nicer notebook error message when lacking permissions
* :ghissue:`5765`: map_sync over empty list raises IndexError
* :ghissue:`5553`: Notebook matplotlib inline backend: can't set figure facecolor
* :ghissue:`5710`: inspect.findsource monkeypatch raises wrong exception for C extensions
* :ghissue:`5706`: Multi-Directory notebooks overwrite each other's checkpoints
* :ghissue:`5698`: can't require a function named ``f``
* :ghissue:`5569`: Keyboard shortcuts in documentation are out of date
* :ghissue:`5566`: Function name printing should use ``__qualname__`` instead of ``__name__`` (Python 3)
* :ghissue:`5676`: "bookmark -l" not working in ipython 2.0
* :ghissue:`5555`: Differentiate more clearly between Notebooks and Folders in new UI
* :ghissue:`5590`: Marked double escape 
* :ghissue:`5514`: import tab-complete fail with ipython 2.0 shell
* :ghissue:`4012`: Notebook: link to markdown formatting reference
* :ghissue:`5611`: Typo in 'storemagic' documentation
* :ghissue:`5589`: Kernel start fails when using --browser argument
* :ghissue:`5491`: Bug in Windows ipconfig ip address regular expression  
* :ghissue:`5579`: rmagic extension throws 'Error while parsing the string.' when last line is comment
* :ghissue:`5518`: Ipython2 will not open ipynb in example directory
* :ghissue:`5561`: New widget documentation has missing notebook link
* :ghissue:`5128`: Page jumping when output from widget interaction replaced
* :ghissue:`5519`: IPython.parallel.Client behavior as iterator
* :ghissue:`5510`: Tab-completion for function argument list


Issues closed in 2.0.0
----------------------


GitHub stats for 2013/08/09 - 2014/04/01 (since 1.0.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

The following 94 authors contributed 3949 commits.

* Aaron Meurer
* Abhinav Upadhyay
* Adam Riggall
* Alex Rudy
* Andrew Mark
* Angus Griffith
* Antony Lee
* Aron Ahmadia
* Arun Persaud
* Benjamin Ragan-Kelley
* Bing Xia
* Blake Griffith
* Bouke van der Bijl
* Bradley M. Froehle
* Brian E. Granger
* Carlos Cordoba
* chapmanb
* chebee7i
* Christoph Gohlke
* Christophe Pradal
* Cyrille Rossant
* Damián Avila
* Daniel B. Vasquez
* Dav Clark
* David Hirschfeld
* David P. Sanders
* David Wyde
* David Österberg
* Doug Blank
* Dražen Lučanin
* epifanio
* Fernando Perez
* Gabriel Becker
* Geert Barentsen
* Hans Meine
* Ingolf Becker
* Jake Vanderplas
* Jakob Gager
* James Porter
* Jason Grout
* Jeffrey Tratner
* Jonah Graham
* Jonathan Frederic
* Joris Van den Bossche
* Juergen Hasch
* Julian Taylor
* Katie Silverio
* Kevin Burke
* Kieran O'Mahony
* Konrad Hinsen
* Kyle Kelley
* Lawrence Fu
* Marc Molla
* Martín Gaitán
* Matt Henderson
* Matthew Brett
* Matthias Bussonnier
* Michael Droettboom
* Mike McKerns
* Nathan Goldbaum
* Pablo de Oliveira
* Pankaj Pandey
* Pascal Schetelat
* Paul Ivanov
* Paul Moore
* Pere Vilas
* Peter Davis
* Philippe Mallet-Ladeira
* Preston Holmes
* Puneeth Chaganti
* Richard Everson
* Roberto Bonvallet
* Samuel Ainsworth
* Sean Vig
* Shashi Gowda
* Skipper Seabold
* Stephan Rave
* Steve Fox
* Steven Silvester
* stonebig
* Susan Tan
* Sylvain Corlay
* Takeshi Kanmae
* Ted Drain
* Thomas A Caswell
* Thomas Kluyver
* Théophile Studer
* Volker Braun
* Wieland Hoffmann
* Yaroslav Halchenko
* Yoval P.
* Yung Siang Liau
* Zachary Sailer
* zah


We closed a total of 1121 issues, 687 pull requests and 434 regular issues;
this is the full list (generated with the script 
:file:`tools/github_stats.py`):

Pull Requests (687):

* :ghpull:`5487`: remove weird unicode space in the new copyright header
* :ghpull:`5476`: For 2.0: Fix links in Notebook Help Menu
* :ghpull:`5337`: Examples reorganization
* :ghpull:`5436`: CodeMirror shortcuts in QuickHelp
* :ghpull:`5444`: Fix numeric verification for Int and Float text widgets.
* :ghpull:`5449`: Stretch keyboard shortcut dialog
* :ghpull:`5473`: Minor corrections of git-hooks setup instructions
* :ghpull:`5471`: Add coding magic comment to nbconvert Python template
* :ghpull:`5452`: print_figure returns unicode for svg
* :ghpull:`5450`: proposal: remove codename
* :ghpull:`5462`: DOC : fixed minor error in using topological sort
* :ghpull:`5463`: make spin_thread tests more forgiving of slow VMs
* :ghpull:`5464`: Fix starting notebook server with file/directory at command line.
* :ghpull:`5453`: remove gitwash
* :ghpull:`5454`: Improve history API docs
* :ghpull:`5431`: update github_stats and gh_api for 2.0
* :ghpull:`5290`: Add dual mode JS tests
* :ghpull:`5451`: check that a handler is actually registered in ShortcutManager.handles
* :ghpull:`5447`: Add %%python2 cell magic
* :ghpull:`5439`: Point to the stable SymPy docs, not the dev docs
* :ghpull:`5437`: Install jquery-ui images
* :ghpull:`5434`: fix check for empty cells in rst template
* :ghpull:`5432`: update links in notebook help menu
* :ghpull:`5435`: Update whatsnew (notebook tour)
* :ghpull:`5433`: Document extraction of octave and R magics
* :ghpull:`5428`: Update COPYING.txt
* :ghpull:`5426`: Separate get_session_info between HistoryAccessor and HistoryManager
* :ghpull:`5419`: move prompts from margin to main column on small screens
* :ghpull:`5430`: Make sure `element` is correct in the context of displayed JS
* :ghpull:`5396`: prevent saving of partially loaded notebooks
* :ghpull:`5429`: Fix tooltip pager feature
* :ghpull:`5330`: Updates to shell reference doc
* :ghpull:`5404`: Fix broken accordion widget
* :ghpull:`5339`: Don't use fork to start the notebook in js tests
* :ghpull:`5320`: Fix for Tooltip & completer click focus bug.
* :ghpull:`5421`: Move configuration of Python test controllers into setup()
* :ghpull:`5418`: fix typo in ssh launcher send_file
* :ghpull:`5403`: remove alt-- shortcut
* :ghpull:`5389`: better log message in deprecated files/ redirect
* :ghpull:`5333`: Fix filenbmanager.list_dirs fails for Windows user profile directory
* :ghpull:`5390`: finish PR #5333
* :ghpull:`5326`: Some gardening on iptest result reporting
* :ghpull:`5375`: remove unnecessary onload hack from mathjax macro
* :ghpull:`5368`: Flexbox classes specificity fixes
* :ghpull:`5331`: fix raw_input CSS
* :ghpull:`5395`: urlencode images for rst files
* :ghpull:`5049`: update quickhelp on adding and removing shortcuts
* :ghpull:`5391`: Fix Gecko (Netscape) keyboard handling
* :ghpull:`5387`: Respect '\r' characters in nbconvert.
* :ghpull:`5399`: Revert PR #5388
* :ghpull:`5388`: Suppress output even when a comment follows ;. Fixes #4525.
* :ghpull:`5394`: nbconvert doc update
* :ghpull:`5359`: do not install less sources
* :ghpull:`5346`: give hint on where to find custom.js
* :ghpull:`5357`: catch exception in copystat
* :ghpull:`5380`: Remove DefineShortVerb... line from latex base template
* :ghpull:`5376`: elide long containers in pretty
* :ghpull:`5310`: remove raw cell placeholder on focus, closes #5238
* :ghpull:`5332`: semantic names for indicator icons
* :ghpull:`5386`: Fix import of socketserver on Python 3
* :ghpull:`5360`: remove some redundant font-family: monospace
* :ghpull:`5379`: don't instantiate Application just for default logger
* :ghpull:`5372`: Don't autoclose strings
* :ghpull:`5296`: unify keyboard shortcut and codemirror interaction
* :ghpull:`5349`: Make Hub.registration_timeout configurable
* :ghpull:`5340`: install bootstrap-tour css
* :ghpull:`5335`: Update docstring for deepreload module
* :ghpull:`5321`: Improve assignment regex to match more tuple unpacking syntax
* :ghpull:`5325`: add NotebookNotary to NotebookApp's class list
* :ghpull:`5313`: avoid loading preprocessors twice
* :ghpull:`5308`: fix HTML capitalization in Highlight2HTML
* :ghpull:`5295`: OutputArea.append_type functions are not prototype methods
* :ghpull:`5318`: Fix local import of select_figure_formats
* :ghpull:`5300`: Fix NameError: name '_rl' is not defined
* :ghpull:`5292`: focus next cell on shift+enter
* :ghpull:`5291`: debug occasional error in test_queue_status
* :ghpull:`5289`: Finishing up #5274 (widget paths fixes)
* :ghpull:`5232`: Make nbconvert html full output like notebook's html.
* :ghpull:`5288`: Correct initial state of kernel status indicator
* :ghpull:`5253`: display any output from this session in terminal console
* :ghpull:`4802`: Tour of the notebook UI (was UI elements inline with highlighting)
* :ghpull:`5285`: Update signature presentation in pinfo classes
* :ghpull:`5268`: Refactoring Notebook.command_mode
* :ghpull:`5226`: Don't run PYTHONSTARTUP file if a file or code is passed
* :ghpull:`5283`: Remove Widget.closed attribute
* :ghpull:`5279`: nbconvert: Make sure node is atleast version 0.9.12
* :ghpull:`5281`: fix a typo introduced by a rebased PR
* :ghpull:`5280`: append Firefox overflow-x fix
* :ghpull:`5277`: check that PIL can save JPEG to BytesIO
* :ghpull:`5044`: Store timestamps for modules to autoreload
* :ghpull:`5278`: Update whatsnew doc from pr files
* :ghpull:`5276`: Fix kernel restart in case connection file is deleted.
* :ghpull:`5272`: allow highlighting language to be set from notebook metadata
* :ghpull:`5158`: log refusal to serve hidden directories
* :ghpull:`5188`: New events system
* :ghpull:`5265`: Missing class def for TimeoutError
* :ghpull:`5267`: normalize unicode in notebook API tests
* :ghpull:`5076`: Refactor keyboard handling
* :ghpull:`5241`: Add some tests for utils
* :ghpull:`5261`: Don't allow edit mode up arrow to continue past index == 0
* :ghpull:`5223`: use on-load event to trigger resizable images
* :ghpull:`5252`: make one strptime call at import of jsonutil
* :ghpull:`5153`: Dashboard sorting
* :ghpull:`5169`: Allow custom header
* :ghpull:`5242`: clear _reply_content cache before using it
* :ghpull:`5194`: require latex titles to be ascii
* :ghpull:`5244`: try to avoid EADDRINUSE errors on travis
* :ghpull:`5245`: support extracted output in HTML template
* :ghpull:`5209`: make input_area css generic to cells
* :ghpull:`5246`: less %pylab, more cowbell!
* :ghpull:`4895`: Improvements to %run completions
* :ghpull:`5243`: Add Javscript to base display priority list.
* :ghpull:`5175`: Audit .html() calls take #2
* :ghpull:`5146`: Dual mode bug fixes.
* :ghpull:`5207`: Children fire event
* :ghpull:`5215`: Dashboard "Running" Tab
* :ghpull:`5240`: Remove unused IPython.nbconvert.utils.console module
* :ghpull:`5239`: Fix exclusion of tests directories from coverage reports
* :ghpull:`5203`: capture some logging/warning output in some tests
* :ghpull:`5216`: fixup positional arg handling in notebook app
* :ghpull:`5229`: get _ipython_display_ method safely
* :ghpull:`5234`: DOC : modified docs is HasTraits.traits and HasTraits.class_traits
* :ghpull:`5221`: Change widget children List to Tuple.
* :ghpull:`5231`: don't forget base_url when updating address bar in rename
* :ghpull:`5173`: Moved widget files into static/widgets/*
* :ghpull:`5222`: Unset PYTHONWARNINGS envvar before running subprocess tests.
* :ghpull:`5172`: Prevent page breaks when printing notebooks via print-view.
* :ghpull:`4985`: Add automatic Closebrackets function to Codemirror.
* :ghpull:`5220`: Make traitlets notify check more robust against classes redefining equality and bool
* :ghpull:`5197`: If there is an error comparing traitlet values when setting a trait, default to go ahead and notify of the new value.
* :ghpull:`5210`: fix pyreadline import in rlineimpl
* :ghpull:`5212`: Wrap nbconvert Markdown/Heading cells in live divs
* :ghpull:`5200`: Allow to pass option to jinja env
* :ghpull:`5202`: handle nodejs executable on debian
* :ghpull:`5112`: band-aid for completion
* :ghpull:`5187`: handle missing output metadata in nbconvert
* :ghpull:`5181`: use gnureadline on OS X
* :ghpull:`5136`: set default value from signature defaults in interact
* :ghpull:`5132`: remove application/pdf->pdf transform in javascript
* :ghpull:`5116`: reorganize who knows what about paths
* :ghpull:`5165`: Don't introspect __call__ for simple callables
* :ghpull:`5170`: Added msg_throttle sync=True widget traitlet
* :ghpull:`5191`: Translate markdown link to rst
* :ghpull:`5037`: FF Fix: alignment and scale of text widget
* :ghpull:`5179`: remove websocket url
* :ghpull:`5110`: add InlineBackend.print_figure_kwargs
* :ghpull:`5147`: Some template URL changes
* :ghpull:`5100`: remove base_kernel_url
* :ghpull:`5163`: Simplify implementation of TemporaryWorkingDirectory.
* :ghpull:`5166`: remove mktemp usage
* :ghpull:`5133`: don't use combine option on ucs package
* :ghpull:`5089`: Remove legacy azure nbmanager
* :ghpull:`5159`: remove append_json reference
* :ghpull:`5095`: handle image size metadata in nbconvert html
* :ghpull:`5156`: fix IPython typo, closes #5155
* :ghpull:`5150`: fix a link that was broken
* :ghpull:`5114`: use non-breaking space for button with no description
* :ghpull:`4778`: add APIs for installing notebook extensions
* :ghpull:`5125`: Fix the display of functions with keyword-only arguments on Python 3.
* :ghpull:`5097`: minor notebook logging changes
* :ghpull:`5047`: only validate package_data when it might be used
* :ghpull:`5121`: fix remove event in KeyboardManager.register_events
* :ghpull:`5119`: Removed 'list' view from Variable Inspector example
* :ghpull:`4925`: Notebook manager api fixes
* :ghpull:`4996`: require print_method to be a bound method
* :ghpull:`5108`: require specifying the version for gh-pages
* :ghpull:`5111`: Minor typo in docstring of IPython.parallel DirectView
* :ghpull:`5098`: mostly debugging changes for IPython.parallel
* :ghpull:`5087`: trust cells with no output
* :ghpull:`5059`: Fix incorrect `Patch` logic in widget code
* :ghpull:`5075`: More flexible box model fixes
* :ghpull:`5091`: Provide logging messages in ipcluster log when engine or controllers fail to start
* :ghpull:`5090`: Print a warning when iptest is run from the IPython source directory
* :ghpull:`5077`: flush replies when entering an eventloop
* :ghpull:`5055`: Minimal changes to import IPython from IronPython
* :ghpull:`5078`: Updating JS tests README.md
* :ghpull:`5083`: don't create js test directories unless they are being used
* :ghpull:`5062`: adjust some events in nb_roundtrip
* :ghpull:`5043`: various unicode / url fixes
* :ghpull:`5066`: remove (almost) all mentions of pylab from our examples
* :ghpull:`4977`: ensure scp destination directories exist (with mkdir -p)
* :ghpull:`5053`: Move&rename JS tests
* :ghpull:`5067`: show traceback in widget handlers
* :ghpull:`4920`: Adding PDFFormatter and kernel side handling of PDF display data
* :ghpull:`5048`: Add edit/command mode indicator
* :ghpull:`5061`: make execute button in menu bar match shift-enter
* :ghpull:`5052`: Add q to toggle the pager.
* :ghpull:`5070`: fix flex: auto
* :ghpull:`5065`: Add example of using annotations in interact
* :ghpull:`5063`: another pass on Interact example notebooks
* :ghpull:`5051`: FF Fix: code cell missing hscroll (2)
* :ghpull:`4960`: Interact/Interactive for widget
* :ghpull:`5045`: Clear timeout in multi-press keyboard shortcuts.
* :ghpull:`5060`: Change 'bind' to 'link'
* :ghpull:`5039`: Expose kernel_info method on inprocess kernel client
* :ghpull:`5058`: Fix iopubwatcher.py example script.
* :ghpull:`5035`: FF Fix: code cell missing hscroll
* :ghpull:`5040`: Polishing some docs
* :ghpull:`5001`: Add directory navigation to dashboard
* :ghpull:`5042`: Remove duplicated Channel ABC classes.
* :ghpull:`5036`: FF Fix: ext link icon same line as link text in help menu
* :ghpull:`4975`: setup.py changes for 2.0
* :ghpull:`4774`: emit event on appended element on dom
* :ghpull:`5023`: Widgets- add ability to pack and unpack arrays on JS side.
* :ghpull:`5003`: Fix pretty reprs of super() objects
* :ghpull:`4974`: make paste focus the pasted cell
* :ghpull:`5012`: Make `SelectionWidget.values` a dict
* :ghpull:`5018`: Prevent 'iptest IPython' from trying to run.
* :ghpull:`5025`: citation2latex filter (using HTMLParser)
* :ghpull:`5027`: pin lessc to 1.4
* :ghpull:`4952`: Widget test inconsistencies
* :ghpull:`5014`: Fix command mode & popup view bug
* :ghpull:`4842`: more subtle kernel indicator
* :ghpull:`5017`: Add notebook examples link to help menu.
* :ghpull:`5015`: don't write cell.trusted to disk
* :ghpull:`5007`: Update whatsnew doc from PR files
* :ghpull:`5010`: Fixes for widget alignment in FF
* :ghpull:`4901`: Add a convenience class to sync traitlet attributes
* :ghpull:`5008`: updated explanation of 'pyin' messages
* :ghpull:`5004`: Fix widget vslider spacing
* :ghpull:`4933`: Small Widget inconsistency fixes
* :ghpull:`4979`: add versioning notes to small message spec changes
* :ghpull:`4893`: add font-awesome 3.2.1
* :ghpull:`4982`: Live readout for slider widgets
* :ghpull:`4813`: make help menu a template
* :ghpull:`4939`: Embed qtconsole docs (continued)
* :ghpull:`4964`: remove shift-= merge keyboard shortcut
* :ghpull:`4504`: Allow input transformers to raise SyntaxError
* :ghpull:`4929`: Fixing various modal/focus related bugs
* :ghpull:`4971`: Fixing issues with js tests
* :ghpull:`4972`: Work around problem in doctest discovery in Python 3.4 with PyQt
* :ghpull:`4937`: pickle arrays with dtype=object
* :ghpull:`4934`: `ipython profile create` respects `--ipython-dir`
* :ghpull:`4954`: generate unicode filename
* :ghpull:`4845`: Add Origin Checking.
* :ghpull:`4916`: Fine tuning the behavior of the modal UI
* :ghpull:`4966`: Ignore sys.argv for NotebookNotary in tests
* :ghpull:`4967`: Fix typo in warning about web socket being closed
* :ghpull:`4965`: Remove mention of iplogger from setup.py
* :ghpull:`4962`: Fixed typos in quick-help text
* :ghpull:`4953`: add utils.wait_for_idle in js tests
* :ghpull:`4870`: ipython_directive, report except/warn in block and add :okexcept: :okwarning: options to suppress
* :ghpull:`4662`: Menu cleanup
* :ghpull:`4824`: sign notebooks
* :ghpull:`4943`: Docs shotgun 4
* :ghpull:`4848`: avoid import of nearby temporary with %edit
* :ghpull:`4950`: Two fixes for file upload related bugs
* :ghpull:`4927`: there shouldn't be a 'files/' prefix in FileLink[s]
* :ghpull:`4928`: use importlib.machinery when available
* :ghpull:`4949`: Remove the docscrape modules, which are part of numpydoc
* :ghpull:`4849`: Various unicode fixes (mostly on Windows)
* :ghpull:`4932`: always point py3compat.input to builtin_mod.input
* :ghpull:`4807`: Correct handling of ansi colour codes when nbconverting to latex
* :ghpull:`4922`: Python nbconvert output shouldn't have output
* :ghpull:`4912`: Skip some Windows io failures
* :ghpull:`4919`: flush output before showing tracebacks
* :ghpull:`4915`: ZMQCompleter inherits from IPCompleter
* :ghpull:`4890`: better cleanup channel FDs
* :ghpull:`4880`: set profile name from profile_dir
* :ghpull:`4853`: fix setting image height/width from metadata
* :ghpull:`4786`: Reduce spacing of heading cells
* :ghpull:`4680`: Minimal pandoc version warning
* :ghpull:`4908`: detect builtin docstrings in oinspect
* :ghpull:`4911`: Don't use `python -m package` on Windows Python 2
* :ghpull:`4909`: sort dictionary keys before comparison, ordering is not guaranteed
* :ghpull:`4374`: IPEP 23: Backbone.js Widgets
* :ghpull:`4903`: use https for all embeds
* :ghpull:`4894`: Shortcut changes
* :ghpull:`4897`: More detailed documentation about kernel_cmd
* :ghpull:`4891`: Squash a few Sphinx warnings from nbconvert.utils.lexers docstrings
* :ghpull:`4679`: JPG compression for inline pylab
* :ghpull:`4708`: Fix indent and center
* :ghpull:`4789`: fix IPython.embed
* :ghpull:`4655`: prefer marked to pandoc for markdown2html
* :ghpull:`4876`: don't show tooltip if object is not found
* :ghpull:`4873`: use 'combine' option to ucs package
* :ghpull:`4732`: Accents in notebook names and in command-line (nbconvert)
* :ghpull:`4867`: Update URL for Lawrence Hall of Science webcam image
* :ghpull:`4868`: Static path fixes
* :ghpull:`4858`: fix tb_offset when running a file
* :ghpull:`4826`: some $.html( -> $.text(
* :ghpull:`4847`: add js kernel_info request
* :ghpull:`4832`: allow NotImplementedError in formatters
* :ghpull:`4803`: BUG: fix cython magic support in ipython_directive
* :ghpull:`4865`: `build` listed twice in .gitignore. Removing one.
* :ghpull:`4851`: fix tooltip token regex for single-character names
* :ghpull:`4846`: Remove some leftover traces of irunner
* :ghpull:`4820`: fix regex for cleaning old logs with ipcluster
* :ghpull:`4844`: adjustments to notebook app logging
* :ghpull:`4840`: Error in Session.send_raw()
* :ghpull:`4819`: update CodeMirror to 3.21
* :ghpull:`4823`: Minor fixes for typos/inconsistencies in parallel docs
* :ghpull:`4811`: document code mirror tab and shift-tab
* :ghpull:`4795`: merge reveal templates
* :ghpull:`4796`: update components
* :ghpull:`4806`: Correct order of packages for unicode in nbconvert to LaTeX
* :ghpull:`4800`: Qt frontend: Handle 'aborted' prompt replies.
* :ghpull:`4794`: Compatibility fix for Python3 (Issue #4783 )
* :ghpull:`4799`: minor js test fix
* :ghpull:`4788`: warn when notebook is started in pylab mode
* :ghpull:`4772`: Notebook server info files
* :ghpull:`4797`: be conservative about kernel_info implementation
* :ghpull:`4787`: non-python kernels run python code with qtconsole
* :ghpull:`4565`: various display type validations
* :ghpull:`4703`: Math macro in jinja templates.
* :ghpull:`4781`: Fix "Source" text for the "Other Syntax" section of the "Typesetting Math" notebook
* :ghpull:`4776`: Manually document py3compat module.
* :ghpull:`4533`: propagate display metadata to all mimetypes
* :ghpull:`4785`: Replacing a for-in loop by an index loop on an array
* :ghpull:`4780`: Updating CSS for UI example.
* :ghpull:`3605`: Modal UI
* :ghpull:`4758`: Python 3.4 fixes
* :ghpull:`4735`: add some HTML error pages
* :ghpull:`4775`: Update whatsnew doc from PR files
* :ghpull:`4760`: Make examples and docs more Python 3 aware
* :ghpull:`4773`: Don't wait forever for notebook server to launch/die for tests
* :ghpull:`4768`: Qt console: Fix _prompt_pos accounting on timer flush output.
* :ghpull:`4727`: Remove Nbconvert template loading magic
* :ghpull:`4763`: Set numpydoc options to produce fewer Sphinx warnings.
* :ghpull:`4770`: always define aliases, even if empty
* :ghpull:`4766`: add `python -m` entry points for everything
* :ghpull:`4767`: remove manpages for irunner, iplogger
* :ghpull:`4751`: Added --post-serve explanation into the nbconvert docs.
* :ghpull:`4762`: whitelist alphanumeric characters for cookie_name
* :ghpull:`4625`: Deprecate %profile magic
* :ghpull:`4745`: warn on failed formatter calls
* :ghpull:`4746`: remove redundant cls alias on Windows
* :ghpull:`4749`: Fix bug in determination of public ips.
* :ghpull:`4715`: restore use of tornado static_url in templates
* :ghpull:`4748`: fix race condition in profiledir creation.
* :ghpull:`4720`: never use ssh multiplexer in tunnels
* :ghpull:`4658`: Bug fix for #4643: Regex object needs to be reset between calls in toolt...
* :ghpull:`4561`: Add Formatter.pop(type)
* :ghpull:`4712`: Docs shotgun 3
* :ghpull:`4713`: Fix saving kernel history in Python 2
* :ghpull:`4744`: don't use lazily-evaluated rc.ids in wait_for_idle
* :ghpull:`4740`: %env can't set variables
* :ghpull:`4737`: check every link when detecting virutalenv
* :ghpull:`4738`: don't inject help into user_ns
* :ghpull:`4739`: skip html nbconvert tests when their dependencies are missing
* :ghpull:`4730`: Fix stripping continuation prompts when copying from Qt console
* :ghpull:`4725`: Doc fixes
* :ghpull:`4656`: Nbconvert HTTP service
* :ghpull:`4710`: make @interactive decorator friendlier with dill
* :ghpull:`4722`: allow purging local results as long as they are not outstanding
* :ghpull:`4549`: Updated IPython console lexers.
* :ghpull:`4570`: Update IPython directive
* :ghpull:`4719`: Fix comment typo in prefilter.py
* :ghpull:`4575`: make sure to encode URL components for API requests
* :ghpull:`4718`: Fixed typo in displaypub
* :ghpull:`4716`: Remove input_prefilter hook
* :ghpull:`4691`: survive failure to bind to localhost in zmq.iostream
* :ghpull:`4696`: don't do anything if add_anchor fails
* :ghpull:`4711`: some typos in the docs
* :ghpull:`4700`: use if main block in entry points
* :ghpull:`4692`: setup.py symlink improvements
* :ghpull:`4265`: JSON configuration file
* :ghpull:`4505`: Nbconvert latex markdown images2
* :ghpull:`4608`: transparent background match ... all colors
* :ghpull:`4678`: allow ipython console to handle text/plain display
* :ghpull:`4706`: remove irunner, iplogger
* :ghpull:`4701`: Delete an old dictionary available for selecting the aligment of text.
* :ghpull:`4702`: Making reveal font-size a relative unit.
* :ghpull:`4649`: added a quiet option to %cpaste to suppress output
* :ghpull:`4690`: Option to spew subprocess streams during tests
* :ghpull:`4688`: Fixed various typos in docstrings.
* :ghpull:`4645`: CasperJs utility functions.
* :ghpull:`4670`: Stop bundling the numpydoc Sphinx extension
* :ghpull:`4675`: common IPython prefix for ModIndex
* :ghpull:`4672`: Remove unused 'attic' module
* :ghpull:`4671`: Fix docstrings in utils.text
* :ghpull:`4669`: add missing help strings to HistoryManager configurables
* :ghpull:`4668`: Make non-ASCII docstring unicode
* :ghpull:`4650`: added a note about sharing of nbconvert tempates
* :ghpull:`4646`: Fixing various output related things:
* :ghpull:`4665`: check for libedit in readline on OS X
* :ghpull:`4606`: Make running PYTHONSTARTUP optional
* :ghpull:`4654`: Fixing left padding of text cells to match that of code cells.
* :ghpull:`4306`: add raw_mimetype metadata to raw cells
* :ghpull:`4576`: Tighten up the vertical spacing on cells and make the padding of cells more consistent
* :ghpull:`4353`: Don't reset the readline completer after each prompt
* :ghpull:`4567`: Adding prompt area to non-CodeCells to indent content.
* :ghpull:`4446`: Use SVG plots in OctaveMagic by default due to lack of Ghostscript on Windows Octave
* :ghpull:`4613`: remove configurable.created
* :ghpull:`4631`: Use argument lists for command help tests
* :ghpull:`4633`: Modifies test_get_long_path_name_winr32() to allow for long path names in temp dir
* :ghpull:`4642`: Allow docs to build without PyQt installed.
* :ghpull:`4641`: Don't check for wx in the test suite.
* :ghpull:`4622`: make QtConsole Lexer configurable
* :ghpull:`4594`: Fixed #2923 Move Save Away from Cut in toolbar
* :ghpull:`4593`: don't interfere with set_next_input contents in qtconsole
* :ghpull:`4640`: Support matplotlib's Gtk3 backend in --pylab mode
* :ghpull:`4639`: Minor import fix to get qtconsole with --pylab=qt working
* :ghpull:`4637`: Fixed typo in links.txt.
* :ghpull:`4634`: Fix nbrun in notebooks with non-code cells.
* :ghpull:`4632`: Restore the ability to run tests from a function.
* :ghpull:`4624`: Fix crash when $EDITOR is non-ASCII
* :ghpull:`4453`: Play nice with App Nap
* :ghpull:`4541`: relax ipconfig matching on Windows
* :ghpull:`4552`: add pickleutil.use_dill
* :ghpull:`4590`: Font awesome for IPython slides
* :ghpull:`4589`: Inherit the width of pre code inside the input code cells.
* :ghpull:`4588`: Update reveal.js CDN to 2.5.0.
* :ghpull:`4569`: store cell toolbar preset in notebook metadata
* :ghpull:`4609`: Fix bytes regex for Python 3.
* :ghpull:`4581`: Writing unicode to stdout
* :ghpull:`4591`: Documenting codemirror shorcuts.
* :ghpull:`4607`: Tutorial doc should link to user config intro
* :ghpull:`4601`: test that rename fails with 409 if it would clobber
* :ghpull:`4599`: re-cast int/float subclasses to int/float in json_clean
* :ghpull:`4542`: new `ipython history clear` subcommand
* :ghpull:`4568`: don't use lazily-evaluated rc.ids in wait_for_idle
* :ghpull:`4572`: DOC: %profile docstring should reference %prun
* :ghpull:`4571`: no longer need 3 suffix on travis, tox
* :ghpull:`4566`: Fixing cell_type in CodeCell constructor.
* :ghpull:`4563`: Specify encoding for reading notebook file.
* :ghpull:`4452`: support notebooks in %run
* :ghpull:`4546`: fix warning condition on notebook startup
* :ghpull:`4540`: Apidocs3
* :ghpull:`4553`: Fix Python 3 handling of urllib
* :ghpull:`4543`: make hiding of initial namespace optional
* :ghpull:`4517`: send shutdown_request on exit of `ipython console`
* :ghpull:`4528`: improvements to bash completion
* :ghpull:`4532`: Hide dynamically defined metaclass base from Sphinx.
* :ghpull:`4515`: Spring Cleaning, and  Load speedup
* :ghpull:`4529`: note routing identities needed for input requests
* :ghpull:`4514`: allow restart in `%run -d`
* :ghpull:`4527`: add redirect for 1.0-style 'files/' prefix links
* :ghpull:`4526`: Allow unicode arguments to passwd_check on Python 2
* :ghpull:`4403`: Global highlight language selection.
* :ghpull:`4250`: outputarea.js: Wrap inline SVGs inside an iframe
* :ghpull:`4521`: Read wav files in binary mode
* :ghpull:`4444`: Css cleaning
* :ghpull:`4523`: Use username and password for MongoDB on ShiningPanda
* :ghpull:`4510`: Update whatsnew from PR files
* :ghpull:`4441`: add `setup.py jsversion`
* :ghpull:`4518`: Fix for race condition in url file decoding.
* :ghpull:`4497`: don't automatically unpack datetime objects in the message spec
* :ghpull:`4506`: wait for empty queues as well as load-balanced tasks
* :ghpull:`4492`: Configuration docs refresh
* :ghpull:`4508`: Fix some uses of map() in Qt console completion code.
* :ghpull:`4498`: Daemon StreamCapturer
* :ghpull:`4499`: Skip clipboard test on unix systems if headless.
* :ghpull:`4460`: Better clipboard handling, esp. with pywin32
* :ghpull:`4496`: Pass nbformat object to write call to save .py script
* :ghpull:`4466`: various pandoc latex fixes
* :ghpull:`4473`: Setup for Python 2/3
* :ghpull:`4459`: protect against broken repr in lib.pretty
* :ghpull:`4457`: Use ~/.ipython as default config directory
* :ghpull:`4489`: check realpath of env in init_virtualenv
* :ghpull:`4490`: fix possible race condition in test_await_data
* :ghpull:`4476`: Fix: Remove space added by display(JavaScript) on page reload
* :ghpull:`4398`: [Notebook] Deactivate tooltip on tab by default.
* :ghpull:`4480`: Docs shotgun 2
* :ghpull:`4488`: fix typo in message spec doc
* :ghpull:`4479`: yet another JS race condition fix
* :ghpull:`4477`: Allow incremental builds of the html_noapi docs target
* :ghpull:`4470`: Various Config object cleanups
* :ghpull:`4410`: make close-and-halt work on new tabs in Chrome
* :ghpull:`4469`: Python 3 & getcwdu
* :ghpull:`4451`: fix: allow JS test to run after shutdown test
* :ghpull:`4456`: Simplify StreamCapturer for subprocess testing
* :ghpull:`4464`: Correct description for Bytes traitlet type
* :ghpull:`4465`: Clean up MANIFEST.in
* :ghpull:`4461`: Correct TypeError message in svg2pdf
* :ghpull:`4458`: use signalstatus if exit status is undefined
* :ghpull:`4438`: Single codebase Python 3 support (again)
* :ghpull:`4198`: Version conversion, support for X to Y even if Y < X (nbformat)
* :ghpull:`4415`: More tooltips in the Notebook menu
* :ghpull:`4450`: remove monkey patch for older versions of tornado
* :ghpull:`4423`: Fix progress bar and scrolling bug.
* :ghpull:`4435`: raise 404 on not found static file
* :ghpull:`4442`: fix and add shim for change introduce by #4195
* :ghpull:`4436`: allow `require("nbextensions/extname")` to load from IPYTHONDIR/nbextensions
* :ghpull:`4437`: don't compute etags in static file handlers
* :ghpull:`4427`: notebooks should always have one checkpoint
* :ghpull:`4425`: fix js pythonisme
* :ghpull:`4195`: IPEP 21:  widget messages
* :ghpull:`4434`: Fix broken link for Dive Into Python.
* :ghpull:`4428`: bump minimum tornado version to 3.1.0
* :ghpull:`4302`: Add an Audio display class
* :ghpull:`4285`: Notebook javascript test suite using CasperJS
* :ghpull:`4420`: Allow checking for backports via milestone
* :ghpull:`4426`: set kernel cwd to notebook's directory
* :ghpull:`4389`: By default, Magics inherit from Configurable
* :ghpull:`4393`: Capture output from subprocs during test, and display on failure
* :ghpull:`4419`: define InlineBackend configurable in its own file
* :ghpull:`4303`: Multidirectory support for the Notebook
* :ghpull:`4371`: Restored ipython profile locate dir and fixed typo. (Fixes #3708).
* :ghpull:`4414`: Specify unicode type properly in rmagic
* :ghpull:`4413`: don't instantiate IPython shell as class attr
* :ghpull:`4400`: Remove 5s wait on inactivity on GUI inputhook loops
* :ghpull:`4412`: Fix traitlet _notify_trait by-ref issue
* :ghpull:`4378`: split adds new cell above, rather than below
* :ghpull:`4405`: Bring display of builtin types and functions in line with Py 2
* :ghpull:`4367`: clean up of documentation files
* :ghpull:`4401`: Provide a name of the HistorySavingThread
* :ghpull:`4384`: fix menubar height measurement
* :ghpull:`4377`: fix tooltip cancel
* :ghpull:`4293`: Factorise code in tooltip for julia monkeypatching
* :ghpull:`4292`: improve js-completer logic.
* :ghpull:`4363`: set_next_input: keep only last input when repeatedly called in a single cell
* :ghpull:`4382`: Use safe_hasattr in dir2
* :ghpull:`4379`: fix (CTRL-M -) shortcut for splitting cell in FF
* :ghpull:`4380`: Test and fixes for localinterfaces
* :ghpull:`4372`: Don't assume that SyntaxTB is always called with a SyntaxError
* :ghpull:`4342`: Return value directly from the try block and avoid a variable
* :ghpull:`4154`: Center LaTeX and figures in markdown
* :ghpull:`4311`: %load -s to load specific functions or classes
* :ghpull:`4350`: WinHPC launcher fixes
* :ghpull:`4345`: Make irunner compatible with upcoming pexpect 3.0 interface
* :ghpull:`4276`: Support container methods in config
* :ghpull:`4359`: test_pylabtools also needs to modify matplotlib.rcParamsOrig
* :ghpull:`4355`: remove hardcoded box-orient
* :ghpull:`4333`: Add Edit Notebook Metadata to Edit menu
* :ghpull:`4349`: Script to update What's New file
* :ghpull:`4348`: Call PDF viewer after latex compiling (nbconvert)
* :ghpull:`4346`: getpass() on Windows & Python 2 needs bytes prompt
* :ghpull:`4304`: use netifaces for faster IPython.utils.localinterfaces
* :ghpull:`4305`: Add even more ways to populate localinterfaces
* :ghpull:`4313`: remove strip_math_space
* :ghpull:`4325`: Some changes to improve readability.
* :ghpull:`4281`: Adjust tab completion widget if too close to bottom of page.
* :ghpull:`4347`: Remove pycolor script
* :ghpull:`4322`: Scroll to the top after change of slides in the IPython slides
* :ghpull:`4289`: Fix scrolling output (not working post clear_output changes)
* :ghpull:`4343`: Make parameters for kernel start method more general
* :ghpull:`4237`: Keywords should shadow magic functions
* :ghpull:`4338`: adjust default value of level in sync_imports
* :ghpull:`4328`: Remove unused loop variable.
* :ghpull:`4340`: fix mathjax download url to new GitHub format
* :ghpull:`4336`: use simple replacement rather than string formatting in format_kernel_cmd
* :ghpull:`4264`: catch unicode error listing profiles
* :ghpull:`4314`: catch EACCES when binding notebook app
* :ghpull:`4324`: Remove commented addthis toolbar
* :ghpull:`4327`: Use the with statement to open a file.
* :ghpull:`4318`: fix initial sys.path
* :ghpull:`4315`: Explicitly state what version of Pandoc is supported in docs/install
* :ghpull:`4316`: underscore missing on notebook_p4
* :ghpull:`4295`: Implement boundary option for load magic (#1093) 
* :ghpull:`4300`: traits defauts are strings not object
* :ghpull:`4297`: Remove an unreachable return statement.
* :ghpull:`4260`: Use subprocess for system_raw
* :ghpull:`4277`: add nbextensions
* :ghpull:`4294`: don't require tornado 3 in `--post serve`
* :ghpull:`4270`: adjust Scheduler timeout logic
* :ghpull:`4278`: add `-a` to easy_install command in libedit warning
* :ghpull:`4282`: Enable automatic line breaks in MathJax.
* :ghpull:`4279`: Fixing line-height of list items in tree view.
* :ghpull:`4253`: fixes #4039.
* :ghpull:`4131`: Add module's name argument in %%cython magic
* :ghpull:`4269`: Add mathletters option and longtable package to latex_base.tplx
* :ghpull:`4230`: Switch correctly to the user's default matplotlib backend after inline.
* :ghpull:`4271`: Hopefully fix ordering of output on ShiningPanda
* :ghpull:`4239`: more informative error message for bad serialization
* :ghpull:`4263`: Fix excludes for IPython.testing
* :ghpull:`4112`: nbconvert: Latex template refactor
* :ghpull:`4261`: Fixing a formatting error in the custom display example notebook.
* :ghpull:`4259`: Fix Windows test exclusions
* :ghpull:`4229`: Clear_output: Animation & widget related changes.
* :ghpull:`4151`: Refactor alias machinery
* :ghpull:`4153`: make timeit return an object that contains values
* :ghpull:`4258`: to-backport label is now 1.2
* :ghpull:`4242`: Allow passing extra arguments to iptest through for nose
* :ghpull:`4257`: fix unicode argv parsing
* :ghpull:`4166`: avoid executing code in utils.localinterfaces at import time
* :ghpull:`4214`: engine ID metadata should be unicode, not bytes
* :ghpull:`4232`: no highlight if no language specified
* :ghpull:`4218`: Fix display of SyntaxError when .py file is modified
* :ghpull:`4207`: add `setup.py css` command
* :ghpull:`4224`: clear previous callbacks on execute
* :ghpull:`4180`: Iptest refactoring
* :ghpull:`4105`: JS output area misaligned
* :ghpull:`4220`: Various improvements to docs formatting
* :ghpull:`4187`: Select adequate highlighter for cell magic languages
* :ghpull:`4228`: update -dev docs to reflect latest stable version
* :ghpull:`4219`: Drop bundled argparse
* :ghpull:`3851`: Adds an explicit newline for pretty-printing.
* :ghpull:`3622`: Drop fakemodule
* :ghpull:`4080`: change default behavior of database task storage
* :ghpull:`4197`: enable cython highlight in notebook
* :ghpull:`4225`: Updated docstring for core.display.Image
* :ghpull:`4175`: nbconvert: Jinjaless exporter base
* :ghpull:`4208`: Added a lightweight "htmlcore" Makefile entry
* :ghpull:`4209`: Magic doc fixes
* :ghpull:`4217`: avoid importing numpy at the module level
* :ghpull:`4213`: fixed dead link in examples/notebooks readme to Part 3
* :ghpull:`4183`: ESC should be handled by CM if tooltip is not on
* :ghpull:`4193`: Update for #3549: Append Firefox overflow-x fix
* :ghpull:`4205`: use TextIOWrapper when communicating with pandoc subprocess
* :ghpull:`4204`: remove some extraneous print statements from IPython.parallel
* :ghpull:`4201`: HeadingCells cannot be split or merged
* :ghpull:`4048`: finish up speaker-notes PR
* :ghpull:`4079`: trigger `Kernel.status_started` after websockets open
* :ghpull:`4186`: moved DummyMod to proper namespace to enable dill pickling
* :ghpull:`4190`: update version-check message in setup.py and IPython.__init__
* :ghpull:`4188`: Allow user_ns trait to be None
* :ghpull:`4189`: always fire LOCAL_IPS.extend(PUBLIC_IPS)
* :ghpull:`4174`: various issues in markdown and rst templates
* :ghpull:`4178`: add missing data_javascript
* :ghpull:`4168`: Py3 failing tests
* :ghpull:`4181`: nbconvert: Fix, sphinx template not removing new lines from headers
* :ghpull:`4043`: don't 'restore_bytes' in from_JSON
* :ghpull:`4149`: reuse more kernels in kernel tests
* :ghpull:`4163`: Fix for incorrect default encoding on Windows.
* :ghpull:`4136`: catch javascript errors in any output
* :ghpull:`4171`: add nbconvert config file when creating profiles
* :ghpull:`4172`: add ability to check what PRs should be backported in backport_pr
* :ghpull:`4167`: --fast flag for test suite!
* :ghpull:`4125`: Basic exercise of `ipython [subcommand] -h` and help-all
* :ghpull:`4085`: nbconvert: Fix sphinx preprocessor date format string for Windows
* :ghpull:`4159`: don't split `.cell` and `div.cell` CSS
* :ghpull:`4165`: Remove use of parametric tests
* :ghpull:`4158`: generate choices for `--gui` configurable from real mapping
* :ghpull:`4083`: Implement a better check for hidden values for %who etc.
* :ghpull:`4147`: Reference notebook examples, fixes #4146.
* :ghpull:`4065`: do not include specific css in embedable one
* :ghpull:`4092`: nbconvert: Fix for unicode html headers, Windows + Python 2.x
* :ghpull:`4074`: close Client sockets if connection fails
* :ghpull:`4064`: Store default codemirror mode in only 1 place
* :ghpull:`4104`: Add way to install MathJax to a particular profile
* :ghpull:`4161`: Select name when renaming a notebook
* :ghpull:`4160`: Add quotes around ".[notebook]" in readme
* :ghpull:`4144`: help_end transformer shouldn't pick up ? in multiline string
* :ghpull:`4090`: Add LaTeX citation handling to nbconvert
* :ghpull:`4143`: update example custom.js
* :ghpull:`4142`: DOC: unwrap openssl line in public_server doc
* :ghpull:`4126`: update tox.ini
* :ghpull:`4141`: add files with a separate `add` call in backport_pr
* :ghpull:`4137`: Restore autorestore option for storemagic
* :ghpull:`4098`: pass profile-dir instead of profile name to Kernel
* :ghpull:`4120`: support `input` in Python 2 kernels
* :ghpull:`4088`: nbconvert: Fix coalescestreams line with incorrect nesting causing strange behavior
* :ghpull:`4060`: only strip continuation prompts if regular prompts seen first
* :ghpull:`4132`: Fixed name error bug in function safe_unicode in module py3compat.
* :ghpull:`4121`: move test_kernel from IPython.zmq to IPython.kernel
* :ghpull:`4118`: ZMQ heartbeat channel: catch EINTR exceptions and continue.
* :ghpull:`4070`: New changes should go into pr/ folder
* :ghpull:`4054`: use unicode for HTML export
* :ghpull:`4106`: fix a couple of default block values
* :ghpull:`4107`: update parallel magic tests with capture_output API
* :ghpull:`4102`: Fix clashes between debugger tests and coverage.py
* :ghpull:`4115`: Update docs on declaring a magic function
* :ghpull:`4101`: restore accidentally removed EngineError
* :ghpull:`4096`: minor docs changes
* :ghpull:`4094`: Update target branch before backporting PR
* :ghpull:`4069`: Drop monkeypatch for pre-1.0 nose
* :ghpull:`4056`: respect `pylab_import_all` when `--pylab` specified at the command-line
* :ghpull:`4091`: Make Qt console banner configurable
* :ghpull:`4086`: fix missing errno import
* :ghpull:`4084`: Use msvcrt.getwch() for Windows pager.
* :ghpull:`4073`: rename ``post_processors`` submodule to ``postprocessors``
* :ghpull:`4075`: Update supported Python versions in tools/test_pr
* :ghpull:`4068`: minor bug fix, define 'cell' in dialog.js.
* :ghpull:`4044`: rename call methods to transform and postprocess
* :ghpull:`3744`: capture rich output as well as stdout/err in capture_output
* :ghpull:`3969`: "use strict" in most (if not all) our javascript
* :ghpull:`4030`: exclude `.git` in MANIFEST.in
* :ghpull:`4047`: Use istype() when checking if canned object is a dict
* :ghpull:`4031`: don't close_fds on Windows
* :ghpull:`4029`: bson.Binary moved
* :ghpull:`3883`: skip test on unix when x11 not available
* :ghpull:`3863`: Added working speaker notes for slides.
* :ghpull:`4035`: Fixed custom jinja2 templates being ignored when setting template_path
* :ghpull:`4002`: Drop Python 2.6 and 3.2
* :ghpull:`4026`: small doc fix in nbconvert
* :ghpull:`4016`: Fix IPython.start_* functions
* :ghpull:`4021`: Fix parallel.client.View map() on numpy arrays
* :ghpull:`4022`: DOC: fix links to matplotlib, notebook docs
* :ghpull:`4018`: Fix warning when running IPython.kernel tests
* :ghpull:`4017`: Add REPL-like printing of final/return value to %%R cell magic
* :ghpull:`4019`: Test skipping without unicode paths
* :ghpull:`4008`: Transform code before %prun/%%prun runs
* :ghpull:`4014`: Fix typo in ipapp
* :ghpull:`3997`: DOC: typos + rewording in examples/notebooks/Cell Magics.ipynb
* :ghpull:`3914`: nbconvert: Transformer tests
* :ghpull:`3987`: get files list in backport_pr
* :ghpull:`3923`: nbconvert: Writer tests
* :ghpull:`3974`: nbconvert: Fix app tests on Window7 w/ Python 3.3
* :ghpull:`3937`: make tab visible in codemirror and light red background
* :ghpull:`3933`: nbconvert: Post-processor tests
* :ghpull:`3978`: fix `--existing` with non-localhost IP
* :ghpull:`3939`: minor checkpoint cleanup
* :ghpull:`3955`: complete on % for magic in notebook
* :ghpull:`3981`: BF: fix nbconert rst input prompt spacing
* :ghpull:`3960`: Don't make sphinx a dependency for importing nbconvert
* :ghpull:`3973`: logging.Formatter is not new-style in 2.6

Issues (434):

* :ghissue:`5476`: For 2.0: Fix links in Notebook Help Menu
* :ghissue:`5337`: Examples reorganization
* :ghissue:`5436`: CodeMirror shortcuts in QuickHelp
* :ghissue:`5444`: Fix numeric verification for Int and Float text widgets.
* :ghissue:`5443`: Int and Float Widgets don't allow negative signs
* :ghissue:`5449`: Stretch keyboard shortcut dialog
* :ghissue:`5471`: Add coding magic comment to nbconvert Python template
* :ghissue:`5470`: UTF-8 Issue When Converting Notebook to a Script.
* :ghissue:`5369`: FormatterWarning for SVG matplotlib output in notebook
* :ghissue:`5460`: Can't start the notebook server specifying a notebook
* :ghissue:`2918`: CodeMirror related issues.
* :ghissue:`5431`: update github_stats and gh_api for 2.0
* :ghissue:`4887`: Add tests for modal UI
* :ghissue:`5290`: Add dual mode JS tests
* :ghissue:`5448`: Cmd+/ shortcut doesn't work in IPython master
* :ghissue:`5447`: Add %%python2 cell magic
* :ghissue:`5442`: Make a "python2" alias or rename the "python"cell magic.
* :ghissue:`2495`: non-ascii characters in the path
* :ghissue:`4554`: dictDB: Exception due to str to datetime comparission
* :ghissue:`5006`: Comm code is not run in the same context as notebook code
* :ghissue:`5118`: Weird interact behavior
* :ghissue:`5401`: Empty code cells in nbconvert rst output cause problems
* :ghissue:`5434`: fix check for empty cells in rst template
* :ghissue:`4944`: Trouble finding ipynb path in Windows 8
* :ghissue:`4605`: Change the url of Editor Shorcuts in the notebook menu.
* :ghissue:`5425`: Update COPYING.txt
* :ghissue:`5348`: BUG: HistoryAccessor.get_session_info(0) - exception
* :ghissue:`5293`: Javascript("element.append()") looks broken.
* :ghissue:`5363`: Disable saving if notebook has stopped loading
* :ghissue:`5189`: Tooltip pager mode is broken
* :ghissue:`5330`: Updates to shell reference doc
* :ghissue:`5397`: Accordion widget broken
* :ghissue:`5106`: Flexbox CSS specificity bugs
* :ghissue:`5297`: tooltip triggers focus bug
* :ghissue:`5417`: scp checking for existence of directories: directory names are incorrect
* :ghissue:`5302`: Parallel engine registration fails for slow engines
* :ghissue:`5334`: notebook's split-cell shortcut dangerous / incompatible with Neo layout (for instance)
* :ghissue:`5324`: Style of `raw_input` UI is off in notebook
* :ghissue:`5350`: Converting notebooks with spaces in their names to RST gives broken images
* :ghissue:`5049`: update quickhelp on adding and removing shortcuts
* :ghissue:`4941`: Eliminating display of intermediate stages in progress bars
* :ghissue:`5345`: nbconvert to markdown does not use backticks
* :ghissue:`5357`: catch exception in copystat
* :ghissue:`5351`: Notebook saving fails on smb share
* :ghissue:`4946`: TeX produced cannot be converted to PDF
* :ghissue:`5347`: pretty print list too slow
* :ghissue:`5238`: Raw cell placeholder is not removed when you edit the cell
* :ghissue:`5382`: Qtconsole doesn't run in Python 3
* :ghissue:`5378`: Unexpected and new conflict between PyFileConfigLoader and IPythonQtConsoleApp
* :ghissue:`4945`: Heading/cells positioning problem and cell output wrapping
* :ghissue:`5084`: Consistent approach for HTML/JS output on nbviewer
* :ghissue:`4902`: print preview does not work, custom.css not found
* :ghissue:`5336`: TypeError in bootstrap-tour.min.js
* :ghissue:`5303`: Changed Hub.registration_timeout to be a config input.
* :ghissue:`995`: Paste-able mode in terminal
* :ghissue:`5305`: Tuple unpacking for shell escape
* :ghissue:`5232`: Make nbconvert html full output like notebook's html.
* :ghissue:`5224`: Audit nbconvert HTML output
* :ghissue:`5253`: display any output from this session in terminal console
* :ghissue:`5251`: ipython console ignoring some stream messages?
* :ghissue:`4802`: Tour of the notebook UI (was UI elements inline with highlighting)
* :ghissue:`5103`: Moving Constructor definition to the top like a Function definition
* :ghissue:`5264`: Test failures on master with Anaconda
* :ghissue:`4833`: Serve /usr/share/javascript at /_sysassets/javascript/ in notebook
* :ghissue:`5071`: Prevent %pylab from clobbering interactive
* :ghissue:`5282`: Exception in widget __del__ methods in Python 3.4.
* :ghissue:`5280`: append Firefox overflow-x fix
* :ghissue:`5120`: append Firefox overflow-x fix, again
* :ghissue:`4127`: autoreload shouldn't rely on .pyc modification times
* :ghissue:`5272`: allow highlighting language to be set from notebook metadata
* :ghissue:`5050`: Notebook cells truncated with Firefox
* :ghissue:`4839`: Error in Session.send_raw()
* :ghissue:`5188`: New events system
* :ghissue:`5076`: Refactor keyboard handling
* :ghissue:`4886`: Refactor and consolidate different keyboard logic in JavaScript code
* :ghissue:`5002`: the green cell border moving forever in Chrome, when there are many code cells.
* :ghissue:`5259`: Codemirror still active in command mode
* :ghissue:`5219`: Output images appear as small thumbnails (Notebook)
* :ghissue:`4829`: Not able to connect qtconsole in Windows 8
* :ghissue:`5152`: Hide __pycache__ in dashboard directory list
* :ghissue:`5151`: Case-insesitive sort for dashboard list
* :ghissue:`4603`: Warn when overwriting a notebook with upload
* :ghissue:`4895`: Improvements to %run completions
* :ghissue:`3459`: Filename completion when run script with %run
* :ghissue:`5225`: Add JavaScript to nbconvert HTML display priority
* :ghissue:`5034`: Audit the places where we call `.html(something)`
* :ghissue:`5094`: Dancing cells in notebook
* :ghissue:`4999`: Notebook focus effects
* :ghissue:`5149`: Clicking on a TextBoxWidget in FF completely breaks dual mode.
* :ghissue:`5207`: Children fire event
* :ghissue:`5227`: display_method of objects with custom __getattr__
* :ghissue:`5236`: Cursor keys do not work to leave Markdown cell while it's being edited
* :ghissue:`5205`: Use CTuple traitlet for Widget children
* :ghissue:`5230`: notebook rename does not respect url prefix
* :ghissue:`5218`: Test failures with Python 3 and enabled warnings
* :ghissue:`5115`: Page Breaks for Print Preview Broken by display: flex - Simple CSS Fix
* :ghissue:`5024`: Make nbconvert HTML output smart about page breaking
* :ghissue:`4985`: Add automatic Closebrackets function to Codemirror.
* :ghissue:`5184`: print '\xa' crashes the interactive shell
* :ghissue:`5214`: Downloading notebook as Python (.py) fails
* :ghissue:`5211`: AttributeError: 'module' object has no attribute '_outputfile'
* :ghissue:`5206`: [CSS?] Inconsistencies in nbconvert divs and IPython Notebook divs?
* :ghissue:`5201`: node != nodejs within Debian packages
* :ghissue:`5112`: band-aid for completion
* :ghissue:`4860`: Completer As-You-Type Broken
* :ghissue:`5116`: reorganize who knows what about paths
* :ghissue:`4973`: Adding security.js with 1st attempt at is_safe
* :ghissue:`5164`: test_oinspect.test_calltip_builtin failure with python3.4
* :ghissue:`5127`: Widgets: skip intermediate callbacks during throttling
* :ghissue:`5013`: Widget alignment differs between FF and Chrome
* :ghissue:`5141`: tornado error static file
* :ghissue:`5160`: TemporaryWorkingDirectory incompatible with python3.4
* :ghissue:`5140`: WIP: %kernels magic
* :ghissue:`4987`: Widget lifecycle problems
* :ghissue:`5129`: UCS package break latex export on non-ascii 
* :ghissue:`4986`: Cell horizontal scrollbar is missing in FF but not in Chrome
* :ghissue:`4685`: nbconvert ignores image size metadata
* :ghissue:`5155`: Notebook logout button does not work (source typo)
* :ghissue:`2678`: Ctrl-m keyboard shortcut clash on Chrome OS
* :ghissue:`5113`: ButtonWidget without caption wrong height.
* :ghissue:`4778`: add APIs for installing notebook extensions
* :ghissue:`5046`: python setup.py failed vs git submodule update worked
* :ghissue:`4925`: Notebook manager api fixes
* :ghissue:`5073`: Cannot align widgets horizontally in the notebook
* :ghissue:`4996`: require print_method to be a bound method
* :ghissue:`4990`: _repr_html_ exception reporting corner case when using type(foo)
* :ghissue:`5099`: Notebook: Changing base_project_url results in failed WebSockets call
* :ghissue:`5096`: Client.map is not fault tolerant
* :ghissue:`4997`: Inconsistent %matplotlib qt behavior
* :ghissue:`5041`: Remove more .html(...) calls.
* :ghissue:`5078`: Updating JS tests README.md
* :ghissue:`4977`: ensure scp destination directories exist (with mkdir -p)
* :ghissue:`3411`: ipython parallel: scp failure.
* :ghissue:`5064`: Errors during interact display at the terminal, not anywhere in the notebook
* :ghissue:`4921`: Add PDF formatter and handling
* :ghissue:`4920`: Adding PDFFormatter and kernel side handling of PDF display data
* :ghissue:`5048`: Add edit/command mode indicator
* :ghissue:`4889`: Add UI element for indicating command/edit modes
* :ghissue:`5052`: Add q to toggle the pager.
* :ghissue:`5000`: Closing pager with keyboard in modal UI
* :ghissue:`5069`: Box model changes broke the Keyboard Shortcuts help modal
* :ghissue:`4960`: Interact/Interactive for widget
* :ghissue:`4883`: Implement interact/interactive for widgets
* :ghissue:`5038`: Fix multiple press keyboard events
* :ghissue:`5054`: UnicodeDecodeError: 'ascii' codec can't decode byte 0xc6 in position 1: ordinal not in range(128)
* :ghissue:`5031`: Bug during integration of IPython console in Qt application
* :ghissue:`5057`: iopubwatcher.py example is broken.
* :ghissue:`4747`: Add event for output_area adding an output
* :ghissue:`5001`: Add directory navigation to dashboard
* :ghissue:`5016`: Help menu external-link icons break layout in FF
* :ghissue:`4885`: Modal UI behavior changes
* :ghissue:`5009`: notebook signatures don't work
* :ghissue:`4975`: setup.py changes for 2.0
* :ghissue:`4774`: emit event on appended element on dom
* :ghissue:`5020`: Python Lists translated to javascript objects in widgets
* :ghissue:`5003`: Fix pretty reprs of super() objects
* :ghissue:`5012`: Make `SelectionWidget.values` a dict
* :ghissue:`4961`: Bug when constructing a selection widget with both values and labels
* :ghissue:`4283`: A `<` in a markdown cell strips cell content when converting to latex
* :ghissue:`4006`: iptest IPython broken
* :ghissue:`4251`: & escaped to &amp; in tex ?
* :ghissue:`5027`: pin lessc to 1.4
* :ghissue:`4323`: Take 2: citation2latex filter (using HTMLParser)
* :ghissue:`4196`: Printing notebook from browser gives 1-page truncated output
* :ghissue:`4842`: more subtle kernel indicator
* :ghissue:`4057`: No path to notebook examples from Help menu
* :ghissue:`5015`: don't write cell.trusted to disk
* :ghissue:`4617`: Changed url link in Help dropdown menu.
* :ghissue:`4976`: Container widget layout broken on Firefox
* :ghissue:`4981`: Vertical slider layout broken
* :ghissue:`4793`: Message spec changes related to `clear_output`
* :ghissue:`4982`: Live readout for slider widgets
* :ghissue:`4813`: make help menu a template
* :ghissue:`4989`: Filename tab completion completely broken
* :ghissue:`1380`: Tab should insert 4 spaces in # comment lines
* :ghissue:`2888`: spaces vs tabs
* :ghissue:`1193`: Allow resizing figures in notebook
* :ghissue:`4504`: Allow input transformers to raise SyntaxError
* :ghissue:`4697`: Problems with height after toggling header and toolbar...
* :ghissue:`4951`: TextWidget to code cell command mode bug.
* :ghissue:`4809`: Arbitrary scrolling (jumping) in clicks in modal UI for notebook
* :ghissue:`4971`: Fixing issues with js tests
* :ghissue:`4972`: Work around problem in doctest discovery in Python 3.4 with PyQt
* :ghissue:`4892`: IPython.qt test failure with python3.4
* :ghissue:`4863`: BUG: cannot create an OBJECT array from memory buffer
* :ghissue:`4704`: Subcommand `profile` ignores --ipython-dir 
* :ghissue:`4845`: Add Origin Checking.
* :ghissue:`4870`: ipython_directive, report except/warn in block and add :okexcept: :okwarning: options to suppress
* :ghissue:`4956`: Shift-Enter does not move to next cell
* :ghissue:`4662`: Menu cleanup
* :ghissue:`4824`: sign notebooks
* :ghissue:`4848`: avoid import of nearby temporary with %edit
* :ghissue:`4731`: %edit files mistakenly import modules in /tmp
* :ghissue:`4950`: Two fixes for file upload related bugs
* :ghissue:`4871`: Notebook upload fails after Delete
* :ghissue:`4825`: File Upload URL set incorrectly
* :ghissue:`3867`: display.FileLinks should work in the exported html verion of a notebook
* :ghissue:`4948`: reveal: ipython css overrides reveal themes
* :ghissue:`4947`: reveal: slides that are too big?
* :ghissue:`4051`: Test failures with Python 3 and enabled warnings
* :ghissue:`3633`: outstanding issues over in ipython/nbconvert repo
* :ghissue:`4087`: Sympy printing in the example notebook
* :ghissue:`4627`: Document various QtConsole embedding approaches.
* :ghissue:`4849`: Various unicode fixes (mostly on Windows)
* :ghissue:`3653`: autocompletion in "from package import <tab>"
* :ghissue:`4583`: overwrite? prompt gets EOFError in 2 process
* :ghissue:`4807`: Correct handling of ansi colour codes when nbconverting to latex
* :ghissue:`4611`: Document how to compile .less files in dev docs.
* :ghissue:`4618`: "Editor Shortcuts" link is broken in help menu dropdown notebook
* :ghissue:`4522`: DeprecationWarning: the sets module is deprecated
* :ghissue:`4368`: No symlink from ipython to ipython3 when inside a python3 virtualenv
* :ghissue:`4234`: Math without $$ doesn't show up when converted to slides
* :ghissue:`4194`: config.TerminalIPythonApp.nosep does not work
* :ghissue:`1491`: prefilter not called for multi-line notebook cells
* :ghissue:`4001`: Windows IPython executable /scripts/ipython not working
* :ghissue:`3959`: think more carefully about text wrapping in nbconvert
* :ghissue:`4907`: Test for traceback depth fails on Windows
* :ghissue:`4906`: Test for IPython.embed() fails on Windows
* :ghissue:`4912`: Skip some Windows io failures
* :ghissue:`3700`: stdout/stderr should be flushed printing exception output... 
* :ghissue:`1181`: greedy completer bug in terminal console
* :ghissue:`2032`: check for a few places we should be using DEFAULT_ENCODING
* :ghissue:`4882`: Too many files open when starting and stopping kernel repeatedly
* :ghissue:`4880`: set profile name from profile_dir
* :ghissue:`4238`: parallel.Client() not using profile that notebook was run with?
* :ghissue:`4853`: fix setting image height/width from metadata
* :ghissue:`4786`: Reduce spacing of heading cells
* :ghissue:`4680`: Minimal pandoc version warning
* :ghissue:`3707`: nbconvert: Remove IPython magic commands from --format="python" output
* :ghissue:`4130`: PDF figures as links from png or svg figures
* :ghissue:`3919`: Allow --profile to be passed a dir.
* :ghissue:`2136`: Handle hard newlines in pretty printer
* :ghissue:`4790`: Notebook modal UI: "merge cell below" key binding, `shift+=`, does not work with some keyboard layouts
* :ghissue:`4884`: Keyboard shortcut changes
* :ghissue:`1184`: slow handling of keyboard input
* :ghissue:`4913`: Mathjax, Markdown, tex, env* and italic
* :ghissue:`3972`: nbconvert: Template output testing
* :ghissue:`4903`: use https for all embeds
* :ghissue:`4874`: --debug does not work if you set .kernel_cmd
* :ghissue:`4679`: JPG compression for inline pylab
* :ghissue:`4708`: Fix indent and center
* :ghissue:`4789`: fix IPython.embed
* :ghissue:`4759`: Application._load_config_files log parameter default fails
* :ghissue:`3153`: docs / file menu: explain how to exit the notebook
* :ghissue:`4791`: Did updates to ipython_directive bork support for cython magic snippets?
* :ghissue:`4385`: "Part 4 - Markdown Cells.ipynb" nbviewer example seems not well referenced in current online documentation page https://ipython.org/ipython-doc/stable/interactive/notebook.htm
* :ghissue:`4655`: prefer marked to pandoc for markdown2html
* :ghissue:`3441`: Fix focus related problems in the notebook
* :ghissue:`3402`: Feature Request: Save As (latex, html,..etc) as a menu option in Notebook rather than explicit need to invoke nbconvert
* :ghissue:`3224`: Revisit layout of notebook area
* :ghissue:`2746`: rerunning a cell with long output (exception) scrolls to much (html notebook)
* :ghissue:`2667`: can't save opened notebook if accidentally delete the notebook in the dashboard
* :ghissue:`3026`: Reporting errors from _repr_<type>_ methods
* :ghissue:`1844`: Notebook does not exist and permalinks
* :ghissue:`2450`: [closed PR] Prevent jumping of window to input when output is clicked.
* :ghissue:`3166`: IPEP 16: Notebook multi directory dashboard and URL mapping
* :ghissue:`3691`: Slight misalignment of Notebook menu bar with focus box
* :ghissue:`4875`: Empty tooltip with `object_found = false` still being shown
* :ghissue:`4432`: The SSL cert for the MathJax CDN is invalid and URL is not protocol agnostic
* :ghissue:`2633`: Help text should leave current cell active
* :ghissue:`3976`: DOC: Pandas link on the notebook help menu?
* :ghissue:`4082`: /new handler redirect cached by browser
* :ghissue:`4298`: Slow ipython --pylab and ipython notebook startup
* :ghissue:`4545`: %store magic not working
* :ghissue:`4610`: toolbar UI enhancements
* :ghissue:`4782`: New modal UI
* :ghissue:`4732`: Accents in notebook names and in command-line (nbconvert)
* :ghissue:`4752`: link broken in docs/examples
* :ghissue:`4835`: running ipython on python files adds an extra traceback frame
* :ghissue:`4792`: repr_html exception warning on qtconsole with pandas  #4745 
* :ghissue:`4834`: function tooltip issues
* :ghissue:`4808`: Docstrings in Notebook not displayed properly and introspection
* :ghissue:`4846`: Remove some leftover traces of irunner
* :ghissue:`4810`: ipcluster bug in clean_logs flag
* :ghissue:`4812`: update CodeMirror for the notebook
* :ghissue:`671`: add migration guide for old IPython config
* :ghissue:`4783`: ipython 2dev  under windows / (win)python 3.3 experiment
* :ghissue:`4772`: Notebook server info files
* :ghissue:`4765`: missing build script for highlight.js
* :ghissue:`4787`: non-python kernels run python code with qtconsole
* :ghissue:`4703`: Math macro in jinja templates.
* :ghissue:`4595`: ipython notebook XSS vulnerable
* :ghissue:`4776`: Manually document py3compat module.
* :ghissue:`4686`: For-in loop on an array in cell.js
* :ghissue:`3605`: Modal UI
* :ghissue:`4769`: Ipython 2.0 will not startup on py27 on windows
* :ghissue:`4482`: reveal.js converter not including CDN by default?
* :ghissue:`4761`: ipv6 address triggers cookie exception
* :ghissue:`4580`: rename or remove %profile magic
* :ghissue:`4643`: Docstring does not open properly
* :ghissue:`4714`: Static URLs are not auto-versioned
* :ghissue:`2573`: document code mirror keyboard shortcuts
* :ghissue:`4717`: hang in parallel.Client when using SSHAgent
* :ghissue:`4544`: Clarify the requirement for pyreadline on Windows
* :ghissue:`3451`: revisit REST /new handler to avoid systematic crawling.
* :ghissue:`2922`: File => Save as '.py' saves magic as code 
* :ghissue:`4728`: Copy/Paste stripping broken in version > 0.13.x in QTConsole
* :ghissue:`4539`: Nbconvert: Latex to PDF conversion fails on notebooks with accented letters
* :ghissue:`4721`: purge_results with jobid crashing - looking for insight
* :ghissue:`4620`: Notebook with ? in title defies autosave, renaming and deletion.
* :ghissue:`4574`: Hash character in notebook name breaks a lot of things
* :ghissue:`4709`: input_prefilter hook not called
* :ghissue:`1680`: qtconsole should support --no-banner and custom banner
* :ghissue:`4689`: IOStream IP address configurable
* :ghissue:`4698`: Missing "if __name__ == '__main__':" check in /usr/bin/ipython
* :ghissue:`4191`: NBConvert: markdown inline and locally referenced files have incorrect file location for latex 
* :ghissue:`2865`: %%!? does not display the shell execute docstring
* :ghissue:`1551`: Notebook should be saved before printing
* :ghissue:`4612`: remove `Configurable.created` ?
* :ghissue:`4629`: Lots of tests fail due to space in sys.executable
* :ghissue:`4644`: Fixed URLs for notebooks
* :ghissue:`4621`: IPython 1.1.0 Qtconsole syntax highlighting highlights python 2 only built-ins when using python 3
* :ghissue:`2923`: Move Delete Button Away from Save Button in the HTML notebook toolbar
* :ghissue:`4615`: UnicodeDecodeError
* :ghissue:`4431`: ipython slow in os x mavericks?
* :ghissue:`4538`: DOC: document how to change ipcontroller-engine.json in case controller was started with --ip="*"
* :ghissue:`4551`: Serialize methods and closures
* :ghissue:`4081`: [Nbconvert][reveal] link to font awesome ?
* :ghissue:`4602`: "ipcluster stop" fails after "ipcluster start --daemonize" using python3.3
* :ghissue:`4578`: NBconvert fails with unicode errors when `--stdout` and file redirection is specified and HTML entities are present
* :ghissue:`4600`: Renaming new notebook to an exist name silently deletes the old one
* :ghissue:`4598`: Qtconsole docstring pop-up fails on method containing defaulted enum argument
* :ghissue:`951`: Remove Tornado monkeypatch
* :ghissue:`4564`: Notebook save failure
* :ghissue:`4562`: nbconvert: Default encoding problem on OS X
* :ghissue:`1675`: add file_to_run=file.ipynb capability to the notebook
* :ghissue:`4516`: `ipython console` doesn't send a `shutdown_request`
* :ghissue:`3043`: can't restart pdb session in ipython
* :ghissue:`4524`: Fix bug with non ascii passwords in notebook login
* :ghissue:`1866`: problems rendering an SVG?
* :ghissue:`4520`: unicode error when trying Audio('data/Bach Cello Suite #3.wav') 
* :ghissue:`4493`: Qtconsole cannot print an ISO8601 date at nanosecond precision
* :ghissue:`4502`: intermittent parallel test failure test_purge_everything 
* :ghissue:`4495`: firefox 25.0: notebooks report "Notebook save failed", .py script save fails, but .ipynb save succeeds
* :ghissue:`4245`: nbconvert latex: code highlighting causes error
* :ghissue:`4486`: Test for whether inside virtualenv does not work if directory is symlinked
* :ghissue:`4485`: Incorrect info in "Messaging in IPython" documentation. 
* :ghissue:`4447`: Ipcontroller broken in current HEAD on windows
* :ghissue:`4241`: Audio display object
* :ghissue:`4463`: Error on empty c.Session.key
* :ghissue:`4454`: UnicodeDecodeError when starting Ipython notebook on a directory containing a file with a non-ascii character
* :ghissue:`3801`: Autocompletion: Fix issue #3723 -- ordering of completions for magic commands and variables with same name
* :ghissue:`3723`: Code completion: 'matplotlib' and '%matplotlib'
* :ghissue:`4396`: Always checkpoint al least once ?
* :ghissue:`2524`: [Notebook] Clear kernel queue
* :ghissue:`2292`: Client side tests for the notebook
* :ghissue:`4424`: Dealing with images in multidirectory environment
* :ghissue:`4388`: Make writing configurable magics easier
* :ghissue:`852`: Notebook should be saved before downloading
* :ghissue:`3708`: ipython profile locate should also work
* :ghissue:`1349`: `?` may generate hundreds of cell 
* :ghissue:`4381`: Using hasattr for trait_names instead of just looking for it directly/using __dir__?
* :ghissue:`4361`: Crash Ultratraceback/ session history
* :ghissue:`3044`: IPython notebook autocomplete for filename string converts multiple spaces to a single space
* :ghissue:`3346`: Up arrow history search shows duplicates in Qtconsole
* :ghissue:`3496`: Fix import errors when running tests from the source directory
* :ghissue:`4114`: If default profile doesn't exist, can't install mathjax to any location
* :ghissue:`4335`: TestPylabSwitch.test_qt fails
* :ghissue:`4291`: serve like option for nbconvert --to latex
* :ghissue:`1824`: Exception before prompting for password during ssh connection
* :ghissue:`4309`: Error in nbconvert - closing </code> tag is not inserted in HTML under some circumstances
* :ghissue:`4351`: /parallel/apps/launcher.py error
* :ghissue:`3603`: Upcoming issues with nbconvert
* :ghissue:`4296`: sync_imports() fails in python 3.3
* :ghissue:`4339`: local mathjax install doesn't work
* :ghissue:`4334`: NotebookApp.webapp_settings static_url_prefix causes crash
* :ghissue:`4308`: Error when use "ipython notebook" in win7 64 with python2.7.3 64.
* :ghissue:`4317`: Relative imports broken in the notebook (Windows)
* :ghissue:`3658`: Saving Notebook clears "Kernel Busy" status from the page and titlebar
* :ghissue:`4312`: Link broken on ipython-doc stable
* :ghissue:`1093`: Add boundary options to %load
* :ghissue:`3619`: Multi-dir webservice design
* :ghissue:`4299`: Nbconvert, default_preprocessors to list of dotted name not list of obj
* :ghissue:`3210`: IPython.parallel tests seem to hang on ShiningPanda
* :ghissue:`4280`: MathJax Automatic Line Breaking
* :ghissue:`4039`: Celltoolbar example issue
* :ghissue:`4247`: nbconvert --to latex: error when converting greek letter
* :ghissue:`4273`: %%capture not capturing rich objects like plots (IPython 1.1.0)
* :ghissue:`3866`: Vertical offsets in LaTeX output for nbconvert
* :ghissue:`3631`: xkcd mode for the IPython notebook
* :ghissue:`4243`: Test exclusions not working on Windows
* :ghissue:`4256`: IPython no longer handles unicode file names 
* :ghissue:`3656`: Audio displayobject
* :ghissue:`4223`: Double output on Ctrl-enter-enter
* :ghissue:`4184`: nbconvert: use r pygmentize backend when highlighting "%%R" cells 
* :ghissue:`3851`: Adds an explicit newline for pretty-printing.
* :ghissue:`3622`: Drop fakemodule
* :ghissue:`4122`: Nbconvert [windows]: Inconsistent line endings in markdown cells exported to latex 
* :ghissue:`3819`: nbconvert add extra blank line to code block on Windows.
* :ghissue:`4203`: remove spurious print statement from parallel annoted functions
* :ghissue:`4200`: Notebook: merging a heading cell and markdown cell cannot be undone
* :ghissue:`3747`: ipynb -> ipynb transformer
* :ghissue:`4024`: nbconvert markdown issues
* :ghissue:`3903`: on Windows, 'ipython3 nbconvert "C:/blabla/first_try.ipynb" --to slides' gives an unexpected result, and '--post serve' fails
* :ghissue:`4095`: Catch js error in append html in stream/pyerr
* :ghissue:`1880`: Add parallelism to test_pr 
* :ghissue:`4085`: nbconvert: Fix sphinx preprocessor date format string for Windows
* :ghissue:`4156`: Specifying --gui=tk at the command line
* :ghissue:`4146`: Having to prepend 'files/' to markdown image paths is confusing 
* :ghissue:`3818`: nbconvert can't handle Heading with Chinese characters on Japanese Windows OS.
* :ghissue:`4134`: multi-line parser fails on ''' in comment, qtconsole and notebook.
* :ghissue:`3998`: sample custom.js needs to be updated
* :ghissue:`4078`: StoreMagic.autorestore not working in 1.0.0
* :ghissue:`3990`: Buitlin `input` doesn't work over zmq
* :ghissue:`4015`: nbconvert fails to convert all the content of a notebook
* :ghissue:`4059`: Issues with Ellipsis literal in Python 3
* :ghissue:`2310`: "ZMQError: Interrupted system call" from RichIPythonWidget
* :ghissue:`3807`: qtconsole ipython 0.13.2 - html/xhtml export fails
* :ghissue:`4103`: Wrong default argument of DirectView.clear
* :ghissue:`4100`: parallel.client.client references undefined error.EngineError
* :ghissue:`484`: Drop nosepatch
* :ghissue:`3350`: Added longlist support in ipdb.
* :ghissue:`1591`: Keying 'q' doesn't quit the interactive help in Wins7
* :ghissue:`40`: The tests in test_process fail under Windows
* :ghissue:`3744`: capture rich output as well as stdout/err in capture_output
* :ghissue:`3742`: %%capture to grab rich display outputs
* :ghissue:`3863`: Added working speaker notes for slides.
* :ghissue:`4013`: Iptest fails in dual python installation
* :ghissue:`4005`: IPython.start_kernel doesn't work.
* :ghissue:`4020`: IPython parallel map fails on numpy arrays
* :ghissue:`3914`: nbconvert: Transformer tests
* :ghissue:`3923`: nbconvert: Writer tests
* :ghissue:`3945`: nbconvert: commandline tests fail Win7x64 Py3.3
* :ghissue:`3937`: make tab visible in codemirror and light red background
* :ghissue:`3935`: No feedback for mixed tabs and spaces
* :ghissue:`3933`: nbconvert: Post-processor tests
* :ghissue:`3977`: unable to complete remote connections for two-process 
* :ghissue:`3939`: minor checkpoint cleanup
* :ghissue:`3955`: complete on % for magic in notebook
* :ghissue:`3954`: all magics should be listed when completing on %
* :ghissue:`3980`: nbconvert rst output lacks needed blank lines
* :ghissue:`3968`: TypeError: super() argument 1 must be type, not classobj (Python 2.6.6)
* :ghissue:`3880`: nbconvert: R&D remaining tests
* :ghissue:`2440`: IPEP 4: Python 3 Compatibility
Issues closed in the 8.x development cycle
==========================================

GitHub stats for 2022/01/05 - 2022/01/12 (tag: 8.0.0rc1)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 26 issues and merged 307 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A8.0>`__

The following 99 authors contributed 372 commits.

* 007vedant
* Adam Hackbarth
* Aditya Sathe
* Ahmed Fasih
* Albert Zhang
* Alex Hall
* Andrew Port
* Ankitsingh6299
* Arthur Moreira
* Ashwin Vishnu
* Augusto
* BaoGiang HoangVu
* bar-hen
* Bart Skowron
* Bartosz Telenczuk
* Bastian Ebeling
* Benjamin Ragan-Kelley
* Blazej Michalik
* blois
* Boyuan Liu
* Brendan Gerrity
* Carol Willing
* Coco Bennett
* Coco Mishra
* Corentin Cadiou
* Daniel Goldfarb
* Daniel Mietchen
* Daniel Shimon
* digitalvirtuoso
* Dimitri Papadopoulos
* dswij
* Eric Wieser
* Erik
* Ethan Madden
* Faris A Chugthai
* farisachugthai
* Gal B
* gorogoroumaru
* Hussaina Begum Nandyala
* Inception95
* Iwan Briquemont
* Jake VanderPlas
* Jakub Klus
* James Morris
* Jay Qi
* Jeroen Bédorf
* Joyce Er
* juacrumar
* Juan Luis Cano Rodríguez
* Julien Rabinow
* Justin Palmer
* Krzysztof Cybulski
* L0uisJ0shua
* lbennett
* LeafyLi
* Lightyagami1
* Lumir Balhar
* Mark Schmitz
* Martin Skarzynski
* martinRenou
* Matt Wozniski
* Matthias Bussonnier
* Meysam Azad
* Michael T
* Michael Tiemann
* Naelson Douglas
* Nathan Goldbaum
* Nick Muoh
* nicolaslazo
* Nikita Kniazev
* NotWearingPants
* Paul Ivanov
* Paulo S. Costa
* Pete Blois
* Peter Corke
* PhanatosZou
* Piers Titus van der Torren
* Rakessh Roshan
* Ram Rachum
* rchiodo
* Reilly Tucker Siemens
* Romulo Filho
* rushabh-v
* Sammy Al Hashemi
* Samreen Zarroug
* Samuel Gaist
* Sanjana-03
* Scott Sanderson
* skalaydzhiyski
* sleeping
* Snir Broshi
* Spas Kalaydzhisyki
* Sylvain Corlay
* Terry Davis
* Timur Kushukov
* Tobias Bengfort
* Tomasz Kłoczko
* Yonatan Goldschmidt
* 谭九鼎
============
 7.x Series
============


.. _version 7.31:

IPython 7.31
============

IPython 7.31 brings a couple of backports and fixes from the 8.0 branches,
it is likely one of the last releases of the 7.x series, as 8.0 will probably be released
between this release and what would have been 7.32.

Please test 8.0 beta/rc releases in addition to this release.

This Releases:
 - Backport some fixes for Python 3.10 (:ghpull:`13412`)
 - use full-alpha transparency on dvipng rendered LaTeX (:ghpull:`13372`)

Many thanks to all the contributors to this release. You can find all individual
contributions to this milestone `on github
<https://github.com/ipython/ipython/milestone/95>`__.

Thanks as well to the `D. E. Shaw group <https://deshaw.com/>`__ for sponsoring
work on IPython and related libraries.


.. _version 7.30:

IPython 7.30
============

IPython 7.30 fixes a couple of bugs introduce in previous releases (in
particular with respect to path handling), and introduce a few features and
improvements:

Notably we will highlight :ghpull:`13267` "Document that ``%run`` can execute
notebooks and ipy scripts.", which is the first commit of Fernando Pérez since
mid 2016 (IPython 5.1). If you are new to IPython, Fernando created IPython in
2001. The other most recent contribution of Fernando to IPython itself was
May 2018, by reviewing and merging PRs. I want to note that Fernando is still
active but mostly as a mentor and leader of the whole Jupyter organisation, but
we're still happy to see him contribute code !

:ghpull:`13290` "Use sphinxify (if available) in object_inspect_mime path"
should allow richer Repr of docstrings when using jupyterlab inspector.

:ghpull:`13311` make the debugger use ``ThreadPoolExecutor`` for debugger cmdloop.
This should fix some issues/infinite loop, but let us know if you come across
any regressions. In particular this fixes issues with `kmaork/madbg <https://github.com/kmaork/madbg>`_,
a remote debugger for IPython.

Note that this is likely the ante-penultimate release of IPython 7.x as a stable
branch, as I hope to release IPython 8.0 as well as IPython 7.31 next
month/early 2022.

IPython 8.0 will drop support for Python 3.7, removed nose as a dependency, and
7.x will only get critical bug fixes with 8.x becoming the new stable. This will
not be possible without `NumFOCUS Small Development Grants
<https://numfocus.org/programs/small-development-grants>`_ Which allowed us to
hire `Nikita Kniazev <https://github.com/Kojoley>`_ who provide Python and C++
help and contracting work.


Many thanks to all the contributors to this release. You can find all individual
contributions to this milestone `on github
<https://github.com/ipython/ipython/milestone/94?closed=1>`__.

Thanks as well to the `D. E. Shaw group <https://deshaw.com/>`__ for sponsoring
work on IPython and related libraries.


.. _version 7.29:

IPython 7.29
============


IPython 7.29 brings a couple of new functionalities to IPython and a number of bugfixes.
It is one of the largest recent release, relatively speaking, with close to 15 Pull Requests.


 - fix an issue where base64 was returned instead of bytes when showing figures :ghpull:`13162`
 - fix compatibility with PyQt6, PySide 6 :ghpull:`13172`. This may be of
   interest if you are running on Apple Silicon as only qt6.2+ is natively
   compatible.
 - fix matplotlib qtagg eventloop :ghpull:`13179`
 - Multiple docs fixes, typos, ... etc.
 - Debugger will now exit by default on SigInt :ghpull:`13218`, this will be
   useful in notebook/lab if you forgot to exit the debugger. "Interrupt Kernel"
   will now exist the debugger.

It give Pdb the ability to skip code in decorators. If functions contain a
special value names ``__debuggerskip__ = True|False``, the function will not be
stepped into, and Pdb will step into lower frames only if the value is set to
``False``. The exact behavior is still likely to have corner cases and will be
refined in subsequent releases. Feedback welcome. See the debugger module
documentation for more info. Thanks to the `D. E. Shaw
group <https://deshaw.com/>`__ for funding this feature.

The main branch of IPython is receiving a number of changes as we received a
`NumFOCUS SDG <https://numfocus.org/programs/small-development-grants>`__
($4800), to help us finish replacing ``nose`` by ``pytest``, and make IPython
future proof with an 8.0 release.


Many thanks to all the contributors to this release. You can find all individual
contributions to this milestone `on github
<https://github.com/ipython/ipython/milestone/93>`__.

Thanks as well to the `D. E. Shaw group <https://deshaw.com/>`__ for sponsoring
work on IPython and related libraries.


.. _version 7.28:

IPython 7.28
============


IPython 7.28 is again a minor release that mostly bring bugfixes, and couple of
improvement. Many thanks to MrMino, who again did all the work this month, and
made a number of documentation improvements.

Here is a non-exhaustive list of changes,

Fixes:

 - async with doesn't allow newlines :ghpull:`13090`
 - Dynamically changing to vi mode via %config magic) :ghpull:`13091`

Virtualenv handling fixes:

 - init_virtualenv now uses Pathlib :ghpull:`12548`
 - Fix Improper path comparison of virtualenv directories :ghpull:`13140`
 - Fix virtual environment user warning for lower case pathes :ghpull:`13094`
 - Adapt to all sorts of drive names for cygwin :ghpull:`13153`

New Features:

 - enable autoplay in embed YouTube player :ghpull:`13133`

 Documentation:

 - Fix formatting for the core.interactiveshell documentation :ghpull:`13118`
 - Fix broken ipyparallel's refs :ghpull:`13138`
 - Improve formatting of %time documentation :ghpull:`13125`
 - Reword the YouTubeVideo autoplay WN :ghpull:`13147`


Highlighted features
--------------------


``YouTubeVideo`` autoplay and the ability to add extra attributes to ``IFrame``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can add any extra attributes to the ``<iframe>`` tag using the new
``extras`` argument in the ``IFrame`` class. For example::

    In [1]: from IPython.display import IFrame

    In [2]: IFrame(src="src", width=300, height=300, extras=['loading="eager"'])

The above cells will result in the following HTML code being displayed in a
notebook::

    <iframe
        width="300"
        height="300"
        src="src"
        frameborder="0"
        allowfullscreen
        loading="eager"
    ></iframe>

Related to the above, the ``YouTubeVideo`` class now takes an
``allow_autoplay`` flag, which sets up the iframe of the embedded YouTube video
such that it allows autoplay.

.. note::
    Whether this works depends on the autoplay policy of the browser rendering
    the HTML allowing it. It also could get blocked by some browser extensions.

Try it out!
::

    In [1]: from IPython.display import YouTubeVideo

    In [2]: YouTubeVideo("dQw4w9WgXcQ", allow_autoplay=True)



Thanks
------

Many thanks to all the contributors to this release. You can find all individual
contributions to this milestone `on github
<https://github.com/ipython/ipython/milestone/92>`__.

Thanks as well to the `D. E. Shaw group <https://deshaw.com/>`__ for sponsoring
work on IPython and related libraries.


.. _version 7.27:

IPython 7.27
============

IPython 7.27 is a minor release that fixes a couple of issues and compatibility.

- Add support for GTK4 :ghpull:`131011`
- Add support for Qt6 :ghpull:`13085`
- Fix an issue with pip magic on windows :ghpull:`13093`

Thanks
------

Many thanks to all the contributors to this release. You can find all individual
contributions to this milestone `on github
<https://github.com/ipython/ipython/milestone/91>`__.

Thanks as well to the `D. E. Shaw group <https://deshaw.com/>`__ for sponsoring
work on IPython and related libraries.

.. _version 7.26:

IPython 7.26
============

IPython 7.26 is a minor release that fixes a couple of issues, updates in API
and Copyright/Licenses issues around various part of the codebase.

We'll highlight `this issue <https://github.com/ipython/ipython/issues/13039>`
pointing out we were including and refereeing to code from Stack Overflow which
was CC-BY-SA, hence incompatible with the BSD license of IPython. This lead us
to a rewriting of the corresponding logic which in our case was done in a more
efficient way (in our case we were searching string prefixes instead of full
strings).

You will notice also a number of documentation improvements and cleanup.

Of particular interest are the following Pull-requests:


 - The IPython directive now uses Sphinx logging for warnings. :ghpull:`13030`.
 - Add expiry days option to pastebin magic and change http protocol to https.
   :ghpull:`13056`
 - Make Ipython.utils.timing work with jupyterlite :ghpull:`13050`.

Pastebin magic expiry days option
---------------------------------

The Pastebin magic now has ``-e`` option to determine 
the number of days for paste expiration. For example
the paste that created with ``%pastebin -e 20 1`` magic will
be available for next 20 days.





Thanks
------

Many thanks to all the contributors to this release and in particular MrMino who
is doing most of the work those days. You can find all individual contributions
to this milestone `on github <https://github.com/ipython/ipython/milestone/90>`__.

Thanks as well to the `D. E. Shaw group <https://deshaw.com/>`__ for sponsoring
work on IPython and related libraries.


.. _version 7.25:

IPython 7.25
============

IPython 7.25 is a minor release that contains a single bugfix, which is highly
recommended for all users of ipdb, ipython debugger %debug magic and similar.

Issuing commands like ``where`` from within the debugger would reset the
local variables changes made by the user. It is interesting to look at the root
cause of the issue as accessing an attribute (``frame.f_locals``) would trigger
this side effects.

Thanks in particular to the patience from the reporters at D.E. Shaw for their
initial bug report that was due to a similar coding oversight in an extension,
and who took time to debug and narrow down the problem.

Thanks
------

Many thanks to all the contributors to this release you can find all individual
contributions to this milestone `on github <https://github.com/ipython/ipython/milestone/89>`__.

Thanks as well to the `D. E. Shaw group <https://deshaw.com/>`__ for sponsoring
work on IPython and related libraries.


.. _version 7.24:

IPython 7.24
============

Third release of IPython for 2021, mostly containing bug fixes. A couple of not
typical updates:

Misc
----


 - Fix an issue where ``%recall`` would both succeeded and print an error message
   it failed. :ghpull:`12952`
 - Drop support for NumPy 1.16 – practically has no effect beyond indicating in
   package metadata that we do not support it. :ghpull:`12937`

Debugger improvements
---------------------

The debugger (and ``%debug`` magic) have been improved and can skip or hide frames
originating from files that are not writable to the user, as these are less
likely to be the source of errors, or be part of system files this can be a useful
addition when debugging long errors.

In addition to the global ``skip_hidden True|False`` command, the debugger has
gained finer grained control of predicates as to whether to a frame should be
considered hidden. So far 3 predicates are available :

  - ``tbhide``: frames containing the local variable ``__tracebackhide__`` set to
    True.
  - ``readonly``: frames originating from readonly files, set to False.
  - ``ipython_internal``: frames that are likely to be from IPython internal
    code, set to True.

You can toggle individual predicates during a session with

.. code-block::

   ipdb> skip_predicates readonly True

Read-only files will now be considered hidden frames.


You can call ``skip_predicates`` without arguments to see the states of current
predicates:

.. code-block::

    ipdb> skip_predicates
    current predicates:
        tbhide : True
        readonly : False
        ipython_internal : True

If all predicates are set to ``False``,  ``skip_hidden`` will practically have
no effect. We attempt to warn you when all predicates are False.

Note that the ``readonly`` predicate may increase disk access as we check for
file access permission for all frames on many command invocation, but is usually
cached by operating systems. Let us know if you encounter any issues.

As the IPython debugger does not use the traitlets infrastructure for
configuration, by editing your ``.pdbrc`` files and appending commands you would
like to be executed just before entering the interactive prompt. For example:


.. code::

    # file : ~/.pdbrc
    skip_predicates readonly True
    skip_predicates tbhide False

Will hide read only frames by default and show frames marked with
``__tracebackhide__``.




Thanks
------

Many thanks to all the contributors to this release you can find all individual
contributions to this milestone `on github <https://github.com/ipython/ipython/milestone/87>`__.

Thanks as well to the `D. E. Shaw group <https://deshaw.com/>`__ for sponsoring
work on IPython and related libraries, in particular above mentioned
improvements to the debugger.




.. _version 7.23:

IPython 7.23 and 7.23.1
=======================


Third release of IPython for 2021, mostly containing bug fixes. A couple of not
typical updates:

 - We moved to GitHub actions away from Travis-CI, the transition may not be
   100% complete (not testing on nightly anymore), but as we ran out of
   Travis-Ci hours on the IPython organisation that was a necessary step.
   :ghpull:`12900`.

 - We have a new dependency: ``matplotlib-inline``, which try to extract
   matplotlib inline backend specific behavior. It is available on PyPI and
   conda-forge thus should not be a problem to upgrade to this version. If you
   are a package maintainer that might be an extra dependency to package first.
   :ghpull:`12817` (IPython 7.23.1 fix a typo that made this change fail)

In the addition/new feature category, ``display()`` now have a ``clear=True``
option to clear the display if any further outputs arrives, allowing users to
avoid having to use ``clear_output()`` directly. :ghpull:`12823`.

In bug fixes category, this release fix an issue when printing tracebacks
containing Unicode characters :ghpull:`12758`.

In code cleanup category :ghpull:`12932` remove usage of some deprecated
functionality for compatibility with Python 3.10.



Thanks
------

Many thanks to all the contributors to this release you can find all individual
contributions to this milestone `on github <https://github.com/ipython/ipython/milestone/86>`__.
In particular MrMino for responding to almost all new issues, and triaging many
of the old ones, as well as takluyver, minrk, willingc for reacting quikly when
we ran out of CI Hours.

Thanks as well to organisations, QuantStack (martinRenou and SylvainCorlay) for
extracting matplotlib inline backend into its own package, and the `D. E. Shaw group
<https://deshaw.com/>`__ for sponsoring work on IPython and related libraries.


.. _version 7.22:

IPython 7.22
============

Second release of IPython for 2021, mostly containing bug fixes. Here is a quick
rundown of the few changes.

- Fix some ``sys.excepthook`` shenanigan when embedding with qt, recommended if
  you – for example – use `napari <https://napari.org>`__. :ghpull:`12842`.
- Fix bug when using the new ipdb ``%context`` magic :ghpull:`12844`
- Couples of deprecation cleanup :ghpull:`12868`
- Update for new dpast.com api if you use the ``%pastbin`` magic. :ghpull:`12712`
- Remove support for numpy before 1.16. :ghpull:`12836`


Thanks
------

We have a new team member that you should see more often on the IPython
repository, Błażej Michalik (@MrMino) have been doing regular contributions to
IPython, and spent time replying to many issues and guiding new users to the
codebase; they now have triage permissions to the IPython repository and we'll
work toward giving them more permission in the future.

Many thanks to all the contributors to this release you can find all individual
contributions to this milestone `on github <https://github.com/ipython/ipython/milestone/84>`__.

Thanks as well to organisations, QuantStack for working on debugger
compatibility for Xeus_python, and the `D. E. Shaw group
<https://deshaw.com/>`__ for sponsoring work on IPython and related libraries.

.. _version 721:

IPython 7.21
============

IPython 7.21 is the first release we have back on schedule of one release every
month; it contains a number of minor fixes and improvements, notably, the new
context command for ipdb


New "context" command in ipdb
-----------------------------

It is now possible to change the number of lines shown in the backtrace
information in ipdb using "context" command. :ghpull:`12826`

(thanks @MrMino, there are other improvement from them on master).

Other notable changes in IPython 7.21
-------------------------------------

- Fix some issues on new osx-arm64 :ghpull:`12804`, :ghpull:`12807`. 
- Compatibility with Xeus-Python for debugger protocol, :ghpull:`12809`
- Misc docs fixes for compatibility and uniformity with Numpydoc.
  :ghpull:`12824`


Thanks
------

Many thanks to all the contributors to this release you can find all individual
contribution to this milestone `on github <https://github.com/ipython/ipython/milestone/83>`__.


.. _version 720:

IPython 7.20
============

IPython 7.20 is the accumulation of 3 month of work on IPython, spacing between
IPython release have been increased from the usual once a month for various
reason.

   - Mainly as I'm too busy and the effectively sole maintainer, and
   - Second because not much changes happened before mid December.

The main driver for this release was the new version of Jedi 0.18 breaking API;
which was taken care of in the master branch early in 2020 but not in 7.x as I
though that by now 8.0 would be out.

The inclusion of a resolver in pip did not help and actually made things worse.
If usually I would have simply pinned Jedi to ``<0.18``; this is not a solution
anymore as now pip is free to install Jedi 0.18, and downgrade IPython.

I'll do my best to keep the regular release, but as the 8.0-dev branch and 7.x
are starting to diverge this is becoming difficult in particular with my limited
time, so if you have any cycles to spare I'll appreciate your help to respond to
issues and pushing 8.0 forward.

Here are thus some of the changes for IPython 7.20.

  - Support for PyQt5 >= 5.11 :ghpull:`12715`
  - ``%reset`` remove imports more agressively :ghpull:`12718`
  - fix the ``%conda`` magic :ghpull:`12739`
  - compatibility with Jedi 0.18, and bump minimum Jedi version. :ghpull:`12793`


.. _version 719:

IPython 7.19
============

IPython 7.19 accumulative two month of works, bug fixes and improvements, there
was exceptionally no release last month.

  - Fix to restore the ability to specify more than one extension using command
    line flags when using traitlets 5.0 :ghpull:`12543`
  - Docs docs formatting that make the install commands work on zsh
    :ghpull:`12587`
  - Always display the last frame in tracebacks even if hidden with
    ``__tracebackhide__`` :ghpull:`12601`
  - Avoid an issue where a callback can be registered multiple times.
    :ghpull:`12625`
  - Avoid an issue in debugger mode where frames changes could be lost.
    :ghpull:`12627`

  - Never hide the frames that invoke a debugger, even if marked as hidden by
    ``__tracebackhide__`` :ghpull:`12631`
  - Fix calling the debugger in a recursive manner :ghpull:`12659`


A number of code changes have landed on master and we are getting close to
enough new features and codebase improvement that a 8.0 start to make sens.
For downstream packages, please start working on migrating downstream testing
away from iptest and using pytest, as nose will not work on Python 3.10 and we
will likely start removing it as a dependency for testing.

.. _version 718:

IPython 7.18
============

IPython 7.18 is a minor release that mostly contains bugfixes.

 - ``CRLF`` is now handled by magics my default; solving some issues due to copy
   pasting on windows. :ghpull:`12475`

 - Requiring pexpect ``>=4.3`` as we are Python 3.7+ only and earlier version of
   pexpect will be incompatible. :ghpull:`12510`

 - Minimum jedi version is now 0.16. :ghpull:`12488`



.. _version 717:

IPython 7.17
============

IPython 7.17 brings a couple of new improvements to API and a couple of user
facing changes to make the terminal experience more user friendly.

:ghpull:`12407` introduces the ability to pass extra argument to the IPython
debugger class; this is to help a new project from ``kmaork``
(https://github.com/kmaork/madbg) to feature a fully remote debugger.

:ghpull:`12410` finally remove support for 3.6, while the codebase is still
technically compatible; IPython will not install on Python 3.6.

lots of work on the debugger and hidden frames from ``@impact27`` in
:ghpull:`12437`, :ghpull:`12445`, :ghpull:`12460` and in particular
:ghpull:`12453` which make the debug magic more robust at handling spaces.

Biggest API addition is code transformation which is done before code execution;
IPython allows a number of hooks to catch non-valid Python syntax (magic, prompt
stripping...etc). Transformers are usually called many time; typically:

 - When trying to figure out whether the code is complete and valid (should we
   insert a new line or execute ?)
 - During actual code execution pass before giving the code to Python's
   ``exec``.

This lead to issues when transformer might have had side effects; or do external
queries. Starting with IPython 7.17 you can expect your transformer to be called
less time.

Input transformers are now called only once in the execution path of
`InteractiveShell`, allowing to register transformer that potentially have side
effects (note that this is not recommended). Internal methods `should_run_async`, and
`run_cell_async` now take a recommended optional `transformed_cell`, and
`preprocessing_exc_tuple` parameters that will become mandatory at some point in
the future; that is to say cells need to be explicitly transformed to be valid
Python syntax ahead of trying to run them. :ghpull:`12440`;

``input_transformers`` can now also have an attribute ``has_side_effects`` set
to `True`, when this attribute is present; this  will prevent the transformers
from being ran when IPython is trying to guess whether the user input is
complete. Note that this may means you will need to explicitly execute in some
case where your transformations are now not ran; but will not affect users with
no custom extensions.


API Changes
-----------

Change of API and exposed objects automatically detected using `frappuccino
<https://pypi.org/project/frappuccino/>`_


 The following items are new since 7.16.0::

     + IPython.core.interactiveshell.InteractiveShell.get_local_scope(self, stack_depth)

 The following signatures differ since 7.16.0::

     - IPython.core.interactiveshell.InteractiveShell.run_cell_async(self, raw_cell, store_history=False, silent=False, shell_futures=True)
     + IPython.core.interactiveshell.InteractiveShell.run_cell_async(self, raw_cell, store_history=False, silent=False, shell_futures=True, *, transformed_cell=None, preprocessing_exc_tuple=None)

     - IPython.core.interactiveshell.InteractiveShell.should_run_async(self, raw_cell)
     + IPython.core.interactiveshell.InteractiveShell.should_run_async(self, raw_cell, *, transformed_cell=None, preprocessing_exc_tuple=None)

     - IPython.terminal.debugger.TerminalPdb.pt_init(self)
     + IPython.terminal.debugger.TerminalPdb.pt_init(self, pt_session_options=None)

This method was added::

     + IPython.core.interactiveshell.InteractiveShell.get_local_scope

Which is now also present on subclasses::

     + IPython.terminal.embed.InteractiveShellEmbed.get_local_scope
     + IPython.terminal.interactiveshell.TerminalInteractiveShell.get_local_scope


.. _version 716:

IPython 7.16.1, 7.16.2
======================

IPython 7.16.1 was release immediately after 7.16.0 to fix a conda packaging issue.
The source is identical to 7.16.0 but the file permissions in the tar are different.

IPython 7.16.2 pins jedi dependency to "<=0.17.2" which should prevent some
issues for users still on python 3.6. This may not be sufficient as pip may
still allow to downgrade IPython.

Compatibility with Jedi > 0.17.2 was not added as this would have meant bumping
the minimal version to >0.16.

IPython 7.16
============


The default traceback mode will now skip frames that are marked with
``__tracebackhide__ = True`` and show how many traceback frames have been
skipped. This can be toggled by using :magic:`xmode` with the ``--show`` or
``--hide`` attribute. It will have no effect on non verbose traceback modes.

The ipython debugger also now understands ``__tracebackhide__`` as well and will
skip hidden frames when displaying. Movement up and down the stack will skip the
hidden frames and will show how many frames were hidden. Internal IPython frames
are also now hidden by default. The behavior can be changed with the
``skip_hidden`` while in the debugger, command and accepts "yes", "no", "true"
and "false" case insensitive parameters.


Misc Noticeable changes:
------------------------

- Exceptions are now (re)raised when running notebooks via the :magic:`%run`, helping to catch issues in workflows and
  pipelines. :ghpull:`12301`
- Fix inputhook for qt 5.15.0 :ghpull:`12355`
- Fix wx inputhook :ghpull:`12375`
- Add handling for malformed pathext env var (Windows) :ghpull:`12367`
- use $SHELL in system_piped :ghpull:`12360` for uniform behavior with
  ipykernel.

Reproducible Build
------------------

IPython 7.15 reproducible build did not work, so we try again this month
:ghpull:`12358`.


API Changes
-----------

Change of API and exposed objects automatically detected using `frappuccino
<https://pypi.org/project/frappuccino/>`_ (still in beta):


The following items are new and mostly related to understanding ``__tracebackhide__``::

    + IPython.core.debugger.Pdb.do_down(self, arg)
    + IPython.core.debugger.Pdb.do_skip_hidden(self, arg)
    + IPython.core.debugger.Pdb.do_up(self, arg)
    + IPython.core.debugger.Pdb.hidden_frames(self, stack)
    + IPython.core.debugger.Pdb.stop_here(self, frame)


The following items have been removed::

    - IPython.core.debugger.Pdb.new_do_down
    - IPython.core.debugger.Pdb.new_do_up

Those were implementation details.


.. _version 715:

IPython 7.15
============

IPython 7.15 brings a number of bug fixes and user facing improvements.

Misc Noticeable changes:
------------------------

 - Long completion name have better elision in terminal :ghpull:`12284`
 - I've started to test on Python 3.9 :ghpull:`12307` and fix some errors.
 - Hi DPI scaling of figures when using qt eventloop :ghpull:`12314`
 - Document the ability to have systemwide configuration for IPython.
   :ghpull:`12328`
 - Fix issues with input autoformatting :ghpull:`12336`
 - ``IPython.core.debugger.Pdb`` is now interruptible (:ghpull:`12168`, in 7.14
   but forgotten in release notes)
 - Video HTML attributes (:ghpull:`12212`, in 7.14 but forgotten in release
   notes)

Reproducible Build
------------------

Starting with IPython 7.15, I am attempting to provide reproducible builds,
that is to say you should be able from the source tree to generate an sdist
and wheel that are identical byte for byte with the publish version on PyPI.

I've only tested on a couple of machines so far and the process is relatively
straightforward, so this mean that IPython not only have a deterministic build
process, but also I have either removed, or put under control all effects of
the build environments on the final artifact.  I encourage you to attempt the
build process on your machine as documented in :ref:`core_developer_guide`
and let me know if you do not obtain an identical artifact.

While reproducible builds is critical to check that the supply chain of (open
source) software has not been compromised, it can also help to speedup many
of the build processes in large environment (conda, apt...) by allowing
better caching of intermediate build steps.

Learn more on `<https://reproducible-builds.org/>`_. `Reflections on trusting
trust <https://dl.acm.org/doi/10.1145/358198.358210>`_ is also one of the
cornerstone and recommended reads on this subject.

.. note::

   The build commit from which the sdist is generated is also `signed
   <https://en.wikipedia.org/wiki/Digital_signature>`_, so you should be able to
   check it has not been compromised, and the git repository is a `merkle-tree
   <https://en.wikipedia.org/wiki/Merkle_tree>`_, you can check the consistency
   with `git-fsck <https://git-scm.com/docs/git-fsck>`_ which you likely `want
   to enable by default
   <https://gist.github.com/mbbx6spp/14b86437e794bffb4120>`_.

NEP29: Last version to support Python 3.6
-----------------------------------------

IPython 7.15 will be the Last IPython version to officially support Python
3.6, as stated by `NumPy Enhancement Proposal 29
<https://numpy.org/neps/nep-0029-deprecation_policy.html>`_. Starting with
next minor version of IPython I may stop testing on Python 3.6 and may stop
publishing release artifacts that install on Python 3.6

Highlighted features
--------------------

Highlighted features are not new, but seem to not be widely known, this
section will help you discover in more narrative form what you can do with
IPython.

Increase Tab Completion Menu Height
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In terminal IPython it is possible to increase the hight of the tab-completion
menu. To do so set the value of
:configtrait:`TerminalInteractiveShell.space_for_menu`, this will reserve more
space at the bottom of the screen for various kind of menus in IPython including
tab completion and searching in history. 

Autoformat Code in the terminal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have a preferred code formatter, you can configure IPython to
reformat your code. Set the value of
:configtrait:`TerminalInteractiveShell.autoformatter` to for example ``'black'``
and IPython will auto format your code when possible.


.. _version 714:

IPython 7.14
============

IPython  7.14 is a minor release that fix a couple of bugs and prepare
compatibility with new or future versions of some libraries. 

Important changes:
------------------

 - Fix compatibility with Sphinx 3+ :ghpull:`12235`
 - Remove deprecated matplotlib parameter usage, compatibility with matplotlib
   3.3+ :`122250`

Misc Changes
------------

 - set ``.py`` extension when editing current buffer in vi/emacs. :ghpull:`12167`
 - support for unicode identifiers in ``?``/``??`` :ghpull:`12208`
 - add extra options to the ``Video`` Rich objects :ghpull:`12212`
 - add pretty-printing to ``SimpleNamespace`` :ghpull:`12230`

IPython.core.debugger.Pdb is now interruptible
----------------------------------------------

A ``KeyboardInterrupt`` will now interrupt IPython's extended debugger, in order to make Jupyter able to interrupt it. (:ghpull:`12168`)

Video HTML attributes
---------------------

Add an option to `IPython.display.Video` to change the attributes of the HTML display of the video (:ghpull:`12212`)


Pending deprecated imports
--------------------------

Many object present in ``IPython.core.display`` are there for internal use only,
and should  already been imported from ``IPython.display`` by users and external
libraries. Trying to import those from ``IPython.core.display`` is still possible
but will trigger a
deprecation warning in later versions of IPython and will become errors in the
future. 

This will simplify compatibility with other Python kernels (like Xeus-Python),
and simplify code base. 




.. _version 713:

IPython 7.13
============

IPython 7.13 is the final release of the 7.x branch since master is diverging
toward an 8.0. Exiting new features have already been merged in 8.0 and will
not be available on the 7.x branch. All the changes below have been backported
from the master branch.


 - Fix inability to run PDB when inside an event loop :ghpull:`12141`
 - Fix ability to interrupt some processes on windows :ghpull:`12137`
 - Fix debugger shortcuts :ghpull:`12132`
 - improve tab completion when inside a string by removing irrelevant elements :ghpull:`12128`
 - Fix display of filename tab completion when the path is long :ghpull:`12122`
 - Many removal of Python 2 specific code path :ghpull:`12110`
 - displaying wav files do not require NumPy anymore, and is 5x to 30x faster :ghpull:`12113`

See the list of all closed issues and pull request on `github
<https://github.com/ipython/ipython/pulls?q=is%3Aclosed+milestone%3A7.13>`_.

.. _version 712:

IPython 7.12
============

IPython 7.12 is a minor update that mostly brings code cleanup, removal of
longtime deprecated function and a couple update to documentation cleanup as well.

Notable changes are the following:

 - Exit non-zero when ipython is given a file path to run that doesn't exist :ghpull:`12074`
 - Test PR on ARM64 with Travis-CI :ghpull:`12073`
 - Update CI to work with latest Pytest :ghpull:`12086`
 - Add infrastructure to run ipykernel eventloop via trio :ghpull:`12097`
 - Support git blame ignore revs :ghpull:`12091`
 - Start multi-line ``__repr__`` s on their own line :ghpull:`12099`

.. _version 7111:

IPython 7.11.1
==============

A couple of deprecated functions (no-op) have been reintroduces in py3compat as
Cython was still relying on them, and will be removed in a couple of versions.

.. _version 711:

IPython 7.11
============

IPython 7.11 received a couple of compatibility fixes and code cleanup.

A number of function in the ``py3compat`` have been removed; a number of types
in the IPython code base are now non-ambiguous and now always ``unicode``
instead of ``Union[Unicode,bytes]``; many of the relevant code path have thus
been simplified/cleaned and types annotation added.

IPython support several verbosity level from exceptions. ``xmode plain`` now
support chained exceptions. :ghpull:`11999`

We are starting to remove ``shell=True`` in some usages of subprocess. While not directly
a security issue (as IPython is made to run arbitrary code anyway) it is not good
practice and we'd like to show the example. :ghissue:`12023`. This discussion
was started by ``@mschwager`` thanks to a new auditing tool they are working on
with duo-labs (`dlint <https://github.com/duo-labs/dlint>`_).

Work around some bugs in Python 3.9 tokenizer :ghpull:`12057`

IPython will now print its version after a crash. :ghpull:`11986`

This is likely the last release from the 7.x series that will see new feature.
The master branch will soon accept large code changes and thrilling new
features; the 7.x branch will only start to accept critical bug fixes, and
update dependencies.

.. _version 7102:

IPython 7.10.2
==============

IPython 7.10.2 fix a couple of extra incompatibility between IPython, ipdb,
asyncio and Prompt Toolkit 3.

.. _version 7101:

IPython 7.10.1
==============

IPython 7.10.1 fix a couple of incompatibilities with Prompt toolkit 3 (please
update Prompt toolkit to 3.0.2 at least), and fixes some interaction with
headless IPython.

.. _version 7100:

IPython 7.10.0
==============

IPython 7.10 is the first double digit minor release in the  last decade, and
first since the release of IPython 1.0, previous double digit minor release was
in August 2009.

We've been trying to give you regular release on the last Friday of every month
for a guaranty of rapid access to bug fixes and new features.

Unlike the previous first few releases that have seen only a couple of code
changes, 7.10 bring a number of changes, new features and bugfixes.

Stop Support for Python 3.5 – Adopt NEP 29
------------------------------------------

IPython has decided to follow the informational `NEP 29
<https://numpy.org/neps/nep-0029-deprecation_policy.html>`_ which layout a clear
policy as to which version of (C)Python and NumPy are supported.

We thus dropped support for Python 3.5, and cleaned up a number of code path
that were Python-version dependant. If you are on 3.5 or earlier pip should
automatically give you the latest compatible version of IPython so you do not
need to pin to a given version.

Support for Prompt Toolkit 3.0
------------------------------

Prompt Toolkit 3.0 was release a week before IPython 7.10 and introduces a few
breaking changes. We believe IPython 7.10 should be compatible with both Prompt
Toolkit 2.x and 3.x, though it has not been extensively tested with 3.x so
please report any issues.


Prompt Rendering Performance improvements
-----------------------------------------

Pull Request :ghpull:`11933` introduced an optimisation in the prompt rendering
logic that should decrease the resource usage of IPython when using the
_default_ configuration but could potentially introduce a regression of
functionalities if you are using a custom prompt.

We know assume if you haven't changed the default keybindings that the prompt
**will not change** during the duration of your input – which is for example
not true when using vi insert mode that switches between `[ins]` and `[nor]`
for the current mode.

If you are experiencing any issue let us know.

Code autoformatting
-------------------

The IPython terminal can now auto format your code just before entering a new
line or executing a command. To do so use the
``--TerminalInteractiveShell.autoformatter`` option and set it to ``'black'``;
if black is installed IPython will use black to format your code when possible.

IPython cannot always properly format your code; in particular it will
auto formatting with *black* will only work if:

   - Your code does not contains magics or special python syntax.

   - There is no code after your cursor.

The Black API is also still in motion; so this may not work with all versions of
black.

It should be possible to register custom formatter, though the API is till in
flux.

Arbitrary Mimetypes Handing in Terminal (Aka inline images in terminal)
-----------------------------------------------------------------------

When using IPython terminal it is now possible to register function to handle
arbitrary mimetypes. While rendering non-text based representation was possible in
many jupyter frontend; it was not possible in terminal IPython, as usually
terminal are limited to displaying text. As many terminal these days provide
escape sequences to display non-text; bringing this loved feature to IPython CLI
made a lot of sens. This functionality will not only allow inline images; but
allow opening of external program; for example ``mplayer`` to "display" sound
files.

So far only the hooks necessary for this are in place, but no default mime
renderers added; so inline images will only be available via extensions. We will
progressively enable these features by default in the next few releases, and
contribution is welcomed.

We welcome any feedback on the API. See :ref:`shell_mimerenderer` for more
informations.

This is originally based on work form in :ghpull:`10610` from @stephanh42
started over two years ago, and still a lot need to be done.

MISC
----

 - Completions can define their own ordering :ghpull:`11855`
 - Enable Plotting in the same cell than the one that import matplotlib
   :ghpull:`11916`
 - Allow to store and restore multiple variables at once :ghpull:`11930`

You can see `all pull-requests <https://github.com/ipython/ipython/pulls?q=is%3Apr+milestone%3A7.10+is%3Aclosed>`_ for this release.

API Changes
-----------

Change of API and exposed objects automatically detected using `frappuccino <https://pypi.org/project/frappuccino/>`_ (still in beta):

The following items are new in IPython 7.10::

    + IPython.terminal.shortcuts.reformat_text_before_cursor(buffer, document, shell)
    + IPython.terminal.interactiveshell.PTK3
    + IPython.terminal.interactiveshell.black_reformat_handler(text_before_cursor)
    + IPython.terminal.prompts.RichPromptDisplayHook.write_format_data(self, format_dict, md_dict='None')

The following items have been removed in 7.10::

    - IPython.lib.pretty.DICT_IS_ORDERED

The following signatures differ between versions::

    - IPython.extensions.storemagic.restore_aliases(ip)
    + IPython.extensions.storemagic.restore_aliases(ip, alias='None')

Special Thanks
--------------

 - @stephanh42 who started the work on inline images in terminal 2 years ago
 - @augustogoulart who spent a lot of time triaging issues and responding to
   users.
 - @con-f-use who is my (@Carreau) first sponsor on GitHub, as a reminder if you
   like IPython, Jupyter and many other library of the SciPy stack you can
   donate to numfocus.org non profit

.. _version 790:

IPython 7.9.0
=============

IPython 7.9 is a small release with a couple of improvement and bug fixes.

 - Xterm terminal title should be restored on exit :ghpull:`11910`
 - special variables ``_``,``__``, ``___`` are not set anymore when cache size
   is 0 or less.  :ghpull:`11877`
 - Autoreload should have regained some speed by using a new heuristic logic to
   find all objects needing reload. This should avoid large objects traversal
   like pandas dataframes. :ghpull:`11876`
 - Get ready for Python 4. :ghpull:`11874`
 - `%env` Magic now has heuristic to hide potentially sensitive values :ghpull:`11896`

This is a small release despite a number of Pull Request Pending that need to
be reviewed/worked on. Many of the core developers have been busy outside of
IPython/Jupyter and we thanks all contributor for their patience; we'll work on
these as soon as we have time.


.. _version780:

IPython 7.8.0
=============

IPython 7.8.0 contain a few bugfix and 2 new APIs:

 - Enable changing the font color for LaTeX rendering :ghpull:`11840`
 - and Re-Expose some PDB API (see below)

Expose Pdb API
--------------

Expose the built-in ``pdb.Pdb`` API. ``Pdb`` constructor arguments are generically
exposed, regardless of python version.
Newly exposed arguments:

- ``skip`` - Python 3.1+
- ``nosiginnt`` - Python 3.2+
- ``readrc`` - Python 3.6+

Try it out::

    from IPython.terminal.debugger import TerminalPdb
    pdb = TerminalPdb(skip=["skipthismodule"])


See :ghpull:`11840`

.. _version770:

IPython 7.7.0
=============

IPython 7.7.0 contain multiple bug fixes and documentation updates; Here are a
few of the outstanding issue fixed:

   - Fix a bug introduced in 7.6 where the ``%matplotlib`` magic would fail on
     previously acceptable arguments :ghpull:`11814`.
   - Fix the manage location on freebsd :ghpull:`11808`.
   - Fix error message about aliases after ``%reset`` call in ipykernel
     :ghpull:`11806`
   - Fix Duplication completions in emacs :ghpull:`11803`

We are planning to adopt `NEP29 <https://github.com/numpy/numpy/pull/14086>`_
(still currently in draft) which may make this minor version of IPython the
last one to support Python 3.5 and will make the code base more aggressive
toward removing compatibility with older versions of Python.

GitHub now support to give only "Triage" permissions to users; if you'd like to
help close stale issues and labels issues please reach to us with your GitHub
Username and we'll add you to the triage team. It is a great way to start
contributing and a path toward getting commit rights.

.. _version761:

IPython 7.6.1
=============

IPython 7.6.1 contain a critical bugfix in the ``%timeit`` magic, which would
crash on some inputs as a side effect of :ghpull:`11716`. See :ghpull:`11812`


.. _whatsnew760:

IPython 7.6.0
=============

IPython 7.6.0 contains a couple of bug fixes and number of small features
additions as well as some compatibility with the current development version of
Python 3.8.

   - Add a ``-l`` option to :magic:`psearch` to list the available search
     types. :ghpull:`11672`
   - Support ``PathLike`` for ``DisplayObject`` and ``Image``. :ghpull:`11764`
   - Configurability of timeout in the test suite for slow platforms.
     :ghpull:`11756`
   - Accept any casing for matplotlib backend. :ghpull:`121748`
   - Properly skip test that requires numpy to be installed :ghpull:`11723`
   - More support for Python 3.8 and positional only arguments (pep570)
     :ghpull:`11720`
   - Unicode names for the completion are loaded lazily on first use which
     should decrease startup time. :ghpull:`11693`
   - Autoreload now update the types of reloaded objects; this for example allow
     pickling of reloaded objects. :ghpull:`11644`
   - Fix a bug where ``%%time`` magic would suppress cell output. :ghpull:`11716`


Prepare migration to pytest (instead of nose) for testing
---------------------------------------------------------

Most of the work between 7.5 and 7.6 was to prepare the migration from our
testing framework to pytest. Most of the test suite should now work by simply
issuing ``pytest`` from the root of the repository.

The migration to pytest is just at its beginning. Many of our test still rely
on IPython-specific plugins for nose using pytest (doctest using IPython syntax
is one example of this where test appear as "passing", while no code has been
ran). Many test also need to be updated like ``yield-test`` to be properly
parametrized tests.

Migration to pytest allowed me to discover a number of issues in our test
suite; which was hiding a number of subtle issues – or not actually running
some of the tests in our test suite – I have thus corrected many of those; like
improperly closed resources; or used of deprecated features. I also made use of
the ``pytest --durations=...`` to find some of our slowest test and speed them
up (our test suite can now be up to 10% faster). Pytest as also a variety of
plugins and flags which will make the code quality of IPython and the testing
experience better.

Misc
----

We skipped the release of 7.6 at the end of May, but will attempt to get back
on schedule. We are starting to think about making introducing backward
incompatible change and start the 8.0 series.

Special Thanks to Gabriel (@gpotter2 on GitHub), who among other took care many
of the remaining task for 7.4 and 7.5, like updating the website.

.. _whatsnew750:

IPython 7.5.0
=============

IPython 7.5.0 consist mostly of bug-fixes, and documentation updates, with one
minor new feature. The `Audio` display element can now be assigned an element
id when displayed in browser. See :ghpull:`11670`

The major outstanding bug fix correct a change of behavior that was introduce
in 7.4.0 where some cell magics would not be able to access or modify global
scope when using the ``@needs_local_scope`` decorator. This was typically
encountered with the ``%%time`` and ``%%timeit`` magics. See :ghissue:`11659`
and :ghpull:`11698`.

.. _whatsnew740:

IPython 7.4.0
=============

Unicode name completions
------------------------

Previously, we provided completion for a unicode name with its relative symbol.
With this, now IPython provides complete suggestions to unicode name symbols.

As on the PR, if user types ``\LAT<tab>``, IPython provides a list of
possible completions. In this case, it would be something like::

   'LATIN CAPITAL LETTER A',
   'LATIN CAPITAL LETTER B',
   'LATIN CAPITAL LETTER C',
   'LATIN CAPITAL LETTER D',
   ....

This help to type unicode character that do not have short latex aliases, and
have long unicode names. for example ``Ͱ``, ``\GREEK CAPITAL LETTER HETA``.

This feature was contributed by Luciana Marques :ghpull:`11583`.

Make audio normalization optional
---------------------------------

Added 'normalize' argument to `IPython.display.Audio`. This argument applies
when audio data is given as an array of samples. The default of `normalize=True`
preserves prior behavior of normalizing the audio to the maximum possible range.
Setting to `False` disables normalization.


Miscellaneous
-------------

 - Fix improper acceptation of ``return`` outside of functions. :ghpull:`11641`.
 - Fixed PyQt 5.11 backwards incompatibility causing sip import failure.
   :ghpull:`11613`.
 - Fix Bug where ``type?`` would crash IPython. :ghpull:`1608`.
 - Allow to apply ``@needs_local_scope`` to cell magics for convenience.
   :ghpull:`11542`.

.. _whatsnew730:

IPython 7.3.0
=============

.. _whatsnew720:

IPython 7.3.0 bring several bug fixes and small improvements that you will
described bellow. 

The biggest change to this release is the implementation of the ``%conda`` and
``%pip`` magics, that will attempt to install packages in the **current
environment**. You may still need to restart your interpreter or kernel for the
change to be taken into account, but it should simplify installation of packages
into remote environment. Installing using pip/conda from the command line is
still the prefer method.

The ``%pip`` magic was already present, but was only printing a warning; now it
will actually forward commands to pip. 

Misc bug fixes and improvements:

 - Compatibility with Python 3.8.
 - Do not expand shell variable in execution magics, and added the
   ``no_var_expand`` decorator for magic requiring a similar functionality
   :ghpull:`11516`
 - Add ``%pip`` and ``%conda`` magic :ghpull:`11524`
 - Re-initialize posix aliases after a ``%reset`` :ghpull:`11528`
 - Allow the IPython command line to run ``*.ipynb`` files :ghpull:`11529`

IPython 7.2.0
=============

IPython 7.2.0 brings minor bugfixes, improvements, and new configuration options:

 - Fix a bug preventing PySide2 GUI integration from working :ghpull:`11464`
 - Run CI on Mac OS ! :ghpull:`11471`
 - Fix IPython "Demo" mode. :ghpull:`11498`
 - Fix ``%run`` magic  with path in name :ghpull:`11499`
 - Fix: add CWD to sys.path *after* stdlib :ghpull:`11502`
 - Better rendering of signatures, especially long ones. :ghpull:`11505`
 - Re-enable jedi by default if it's installed :ghpull:`11506`
 - Add New ``minimal`` exception reporting mode (useful for educational purpose). See :ghpull:`11509`


Added ability to show subclasses when using pinfo and other utilities
---------------------------------------------------------------------

When using ``?``/``??`` on a class, IPython will now list the first 10 subclasses.

Special Thanks to Chris Mentzel of the Moore Foundation for this feature. Chris
is one of the people who played a critical role in IPython/Jupyter getting
funding.

We are grateful for all the help Chris has given us over the years,
and we're now proud to have code contributed by Chris in IPython.

OSMagics.cd_force_quiet configuration option
--------------------------------------------

You can set this option to force the %cd magic to behave as if ``-q`` was passed:
::

    In [1]: cd /
    /

    In [2]: %config OSMagics.cd_force_quiet = True

    In [3]: cd /tmp

    In [4]:

See :ghpull:`11491`

In vi editing mode, whether the prompt includes the current vi mode can now be configured
-----------------------------------------------------------------------------------------

Set the ``TerminalInteractiveShell.prompt_includes_vi_mode`` to a boolean value
(default: True) to control this feature. See :ghpull:`11492`

.. _whatsnew710:

IPython 7.1.0
=============

IPython 7.1.0 is the first minor release after 7.0.0 and mostly brings fixes to
new features, internal refactoring, and fixes for regressions that happened during the 6.x->7.x
transition. It also brings **Compatibility with Python 3.7.1**, as we're
unwillingly relying on a bug in CPython.

New Core Dev:

 - We welcome Jonathan Slenders to the commiters. Jonathan has done a fantastic
   work on prompt_toolkit, and we'd like to recognise his impact by giving him
   commit rights. :ghissue:`11397`

Notable Changes

 - Major update of "latex to unicode" tab completion map (see below)

Notable New Features:

 - Restore functionality and documentation of the **sphinx directive**, which
   is now stricter (fail on error by daefault), has new configuration options,
   has a brand new documentation page :ref:`ipython_directive` (which needs
   some cleanup). It is also now *tested* so we hope to have less regressions.
   :ghpull:`11402`

 - ``IPython.display.Video`` now supports ``width`` and ``height`` arguments,
   allowing a custom width and height to be set instead of using the video's
   width and height. :ghpull:`11353`

 - Warn when using ``HTML('<iframe>')`` instead of ``IFrame`` :ghpull:`11350`

 - Allow Dynamic switching of editing mode between vi/emacs and show
   normal/input mode in prompt when using vi. :ghpull:`11390`. Use ``%config
   TerminalInteractiveShell.editing_mode = 'vi'`` or ``%config
   TerminalInteractiveShell.editing_mode = 'emacs'`` to dynamically switch
   between modes.


Notable Fixes:

 - Fix entering of **multi-line blocks in terminal** IPython, and various
   crashes in the new input transformation machinery :ghpull:`11354`,
   :ghpull:`11356`, :ghpull:`11358`. These also fix a **Compatibility bug
   with Python 3.7.1**.

 - Fix moving through generator stack in ipdb :ghpull:`11266`

 - %Magic command arguments now support quoting. :ghpull:`11330`

 - Re-add ``rprint`` and ``rprinte`` aliases. :ghpull:`11331`

 - Remove implicit dependency on ``ipython_genutils`` :ghpull:`11317`

 - Make ``nonlocal`` raise ``SyntaxError`` instead of silently failing in async
   mode. :ghpull:`11382`

 - Fix mishandling of magics and ``= !`` assignment just after a dedent in
   nested code blocks :ghpull:`11418`

 - Fix instructions for custom shortcuts :ghpull:`11426`


Notable Internals improvements:

 - Use of ``os.scandir`` (Python 3 only) to speed up some file system operations.
   :ghpull:`11365`

 - use ``perf_counter`` instead of ``clock`` for more precise
   timing results with ``%time`` :ghpull:`11376`

Many thanks to all the contributors and in particular to ``bartskowron`` and
``tonyfast`` who handled some pretty complicated bugs in the input machinery. We
had a number of first time contributors and maybe hacktoberfest participants that
made significant contributions and helped us free some time to focus on more
complicated bugs.

You
can see all the closed issues and Merged PR, new features and fixes `here
<https://github.com/ipython/ipython/issues?utf8=%E2%9C%93&q=+is%3Aclosed+milestone%3A7.1+>`_.

Unicode Completion update
-------------------------

In IPython 7.1 the Unicode completion map has been updated and synchronized with
the Julia language.

Added and removed character characters:

 ``\jmath`` (``ȷ``), ``\\underleftrightarrow`` (U+034D, combining) have been
 added, while ``\\textasciicaron`` have been removed

Some sequences have seen their prefix removed:

 - 6 characters ``\text...<tab>`` should now be inputed with ``\...<tab>`` directly,
 - 45 characters ``\Elz...<tab>`` should now be inputed with ``\...<tab>`` directly,
 - 65 characters ``\B...<tab>`` should now be inputed with ``\...<tab>`` directly,
 - 450 characters ``\m...<tab>`` should now be inputed with ``\...<tab>`` directly,

Some sequences have seen their prefix shortened:

 - 5 characters ``\mitBbb...<tab>`` should now be inputed with ``\bbi...<tab>`` directly,
 - 52 characters ``\mit...<tab>`` should now be inputed with ``\i...<tab>`` directly,
 - 216 characters ``\mbfit...<tab>`` should now be inputed with ``\bi...<tab>`` directly,
 - 222 characters ``\mbf...<tab>`` should now be inputed with ``\b...<tab>`` directly,

A couple of characters had their sequence simplified:

 - ``ð``, type ``\dh<tab>``, instead of ``\eth<tab>``
 - ``ħ``, type ``\hbar<tab>``, instead of ``\Elzxh<tab>``
 - ``ɸ``, type ``\ltphi<tab>``, instead of ``\textphi<tab>``
 - ``ϴ``, type ``\varTheta<tab>``, instead of ``\textTheta<tab>``
 - ``ℇ``, type ``\eulermascheroni<tab>``, instead of ``\Eulerconst<tab>``
 - ``ℎ``, type ``\planck<tab>``, instead of ``\Planckconst<tab>``

 - U+0336 (COMBINING LONG STROKE OVERLAY), type ``\strike<tab>``, instead of ``\Elzbar<tab>``.

A couple of sequences have been updated:

 - ``\varepsilon`` now gives ``ɛ`` (GREEK SMALL LETTER EPSILON) instead of ``ε`` (GREEK LUNATE EPSILON SYMBOL),
 - ``\underbar`` now gives U+0331 (COMBINING MACRON BELOW) instead of U+0332 (COMBINING LOW LINE).


.. _whatsnew700:

IPython 7.0.0
=============

Released Thursday September 27th, 2018

IPython 7 includes major feature improvements.
This is also the second major version of IPython to support only
Python 3 – starting at Python 3.4. Python 2 is still community-supported
on the bugfix only 5.x branch, but we remind you that Python 2 "end of life"
is on Jan 1st 2020.

We were able to backport bug fixes to the 5.x branch thanks to our backport bot which
backported more than `70 Pull-Requests
<https://github.com/ipython/ipython/pulls?page=3&q=is%3Apr+sort%3Aupdated-desc+author%3Aapp%2Fmeeseeksdev++5.x&utf8=%E2%9C%93>`_, but there are still many PRs that required manual work. This is an area of the project where you can easily contribute by looking for `PRs that still need manual backport <https://github.com/ipython/ipython/issues?q=label%3A%22Still+Needs+Manual+Backport%22+is%3Aclosed+sort%3Aupdated-desc>`_

The IPython 6.x branch will likely not see any further release unless critical
bugs are found.

Make sure you have pip > 9.0 before upgrading. You should be able to update by running:

.. code::

    pip install ipython --upgrade

.. only:: ipydev

  If you are trying to install or update an ``alpha``, ``beta``, or ``rc``
  version, use pip ``--pre`` flag.

  .. code::

      pip install ipython --upgrade --pre


Or, if you have conda installed: 

.. code::
   
   conda install ipython



Prompt Toolkit 2.0
------------------

IPython 7.0+ now uses ``prompt_toolkit 2.0``. If you still need to use an earlier
``prompt_toolkit`` version, you may need to pin IPython to ``<7.0``.

Autowait: Asynchronous REPL
---------------------------

Staring with IPython 7.0 on Python 3.6+, IPython can automatically ``await``
top level code. You should not need to access an event loop or runner
yourself. To learn more, read the :ref:`autoawait` section of our docs, see
:ghpull:`11265`, or try the following code::

    Python 3.6.0
    Type 'copyright', 'credits' or 'license' for more information
    IPython 7.0.0 -- An enhanced Interactive Python. Type '?' for help.

    In [1]: import aiohttp
       ...: result = aiohttp.get('https://api.github.com')

    In [2]: response = await result
    <pause for a few 100s ms>

    In [3]: await response.json()
    Out[3]:
    {'authorizations_url': 'https://api.github.com/authorizations',
     'code_search_url': 'https://api.github.com/search/code?q={query}{&page,per_page,sort,order}',
    ...
    }

.. note::

   Async integration is experimental code, behavior may change or be removed
   between Python and IPython versions without warnings.

Integration is by default with `asyncio`, but other libraries can be configured --
like ``curio`` or ``trio`` -- to improve concurrency in the REPL::

    In [1]: %autoawait trio

    In [2]: import trio

    In [3]: async def child(i):
       ...:     print("   child %s goes to sleep"%i)
       ...:     await trio.sleep(2)
       ...:     print("   child %s wakes up"%i)

    In [4]: print('parent start')
       ...: async with trio.open_nursery() as n:
       ...:     for i in range(3):
       ...:         n.spawn(child, i)
       ...: print('parent end')
    parent start
       child 2 goes to sleep
       child 0 goes to sleep
       child 1 goes to sleep
       <about 2 seconds pause>
       child 2 wakes up
       child 1 wakes up
       child 0 wakes up
    parent end

See :ref:`autoawait` for more information.


Asynchronous code in a Notebook interface or any other frontend using the
Jupyter Protocol will require further updates to the IPykernel package.

Non-Asynchronous code
~~~~~~~~~~~~~~~~~~~~~

As the internal API of IPython is now asynchronous, IPython needs to run under
an event loop. In order to allow many workflows, (like using the :magic:`%run`
magic, or copy-pasting code that explicitly starts/stop event loop), when
top-level code is detected as not being asynchronous, IPython code is advanced
via a pseudo-synchronous runner, and may not advance pending tasks.

Change to Nested Embed
~~~~~~~~~~~~~~~~~~~~~~

The introduction of the ability to run async code had some effect on the
``IPython.embed()`` API. By default, embed will not allow you to run asynchronous
code unless an event loop is specified.

Effects on Magics
~~~~~~~~~~~~~~~~~

Some magics will not work with async until they're updated.
Contributions welcome.

Expected Future changes
~~~~~~~~~~~~~~~~~~~~~~~

We expect more internal but public IPython functions to become ``async``, and
will likely end up having a persistent event loop while IPython is running.

Thanks
~~~~~~

This release took more than a year in the making.
The code was rebased a number of
times; leading to commit authorship that may have been lost in the final
Pull-Request. Huge thanks to many people for contribution, discussion, code,
documentation, use-cases: dalejung, danielballan, ellisonbg, fperez, gnestor,
minrk, njsmith, pganssle, tacaswell, takluyver , vidartf ... And many others.


Autoreload Improvement
----------------------

The magic :magic:`%autoreload 2 <autoreload>` now captures new methods added to
classes. Earlier, only methods existing as of the initial import were being
tracked and updated.  

This new feature helps dual environment development - Jupyter+IDE - where the
code gradually moves from notebook cells to package files as it gets
structured.

**Example**: An instance of the class ``MyClass`` will be able to access the
method ``cube()`` after it is uncommented and the file ``file1.py`` is saved on
disk.


.. code::

   # notebook

   from mymodule import MyClass
   first = MyClass(5)

.. code::

   # mymodule/file1.py

   class MyClass:

       def __init__(self, a=10):
           self.a = a

       def square(self):
           print('compute square')
           return self.a*self.a

       # def cube(self):
       #     print('compute cube')
       #     return self.a*self.a*self.a




Misc
----

The autoindent feature that was deprecated in 5.x was re-enabled and
un-deprecated in :ghpull:`11257`

Make :magic:`%run -n -i ... <run>` work correctly. Earlier, if :magic:`%run` was
passed both arguments, ``-n`` would be silently ignored. See :ghpull:`10308`


The :cellmagic:`%%script` (as well as :cellmagic:`%%bash`,
:cellmagic:`%%ruby`... ) cell magics now raise by default if the return code of
the given code is non-zero (thus halting execution of further cells in a
notebook). The behavior can be disable by passing the ``--no-raise-error`` flag.


Deprecations
------------

A couple of unused functions and methods have been deprecated and will be removed
in future versions:

  - ``IPython.utils.io.raw_print_err``
  - ``IPython.utils.io.raw_print``

  
Backwards incompatible changes
------------------------------

* The API for transforming input before it is parsed as Python code has been
  completely redesigned: any custom input transformations will need to be
  rewritten. See :doc:`/config/inputtransforms` for details of the new API.
============
 2.x Series
============

Release 2.4
===========

January, 2014

.. note::

    Some of the patches marked for 2.4 were left out of 2.4.0.
    Please use 2.4.1.

- backport read support for nbformat v4 from IPython 3
- support for PyQt5 in the kernel (not QtConsole)
- support for Pygments 2.0

For more information on what fixes have been backported to 2.4,
see our :ref:`detailed release info <issues_list_200>`.


Release 2.3.1
=============

November, 2014

- Fix CRCRLF line-ending bug in notebooks on Windows

For more information on what fixes have been backported to 2.3.1,
see our :ref:`detailed release info <issues_list_200>`.

Release 2.3.0
=============

October, 2014

- improve qt5 support
- prevent notebook data loss with atomic writes

For more information on what fixes have been backported to 2.3,
see our :ref:`detailed release info <issues_list_200>`.

Release 2.2.0
=============

August, 2014

- Add CORS configuration

For more information on what fixes have been backported to 2.2,
see our :ref:`detailed release info <issues_list_200>`.

Release 2.1.0
=============

May, 2014

IPython 2.1 is the first bugfix release for 2.0.
For more information on what fixes have been backported to 2.1,
see our :ref:`detailed release info
<issues_list_200>`.


Release 2.0.0
=============

April, 2014

IPython 2.0 requires Python ≥ 2.7.2 or ≥ 3.3.0.
It does not support Python 3.0, 3.1, 3.2, 2.5, or 2.6.

The principal milestones of 2.0 are:

- interactive widgets for the notebook
- directory navigation in the notebook dashboard
- persistent URLs for notebooks
- a new modal user interface in the notebook
- a security model for notebooks

Contribution summary since IPython 1.0 in August, 2013:

- ~8 months of work
- ~650 pull requests merged
- ~400 issues closed (non-pull requests)
- contributions from ~100 authors
- ~4000 commits

The amount of work included in this release is so large that we can only cover
here the main highlights; please see our :ref:`detailed release statistics
<issues_list_200>` for links to every issue and pull request closed on GitHub
as well as a full list of individual contributors.

New stuff in the IPython notebook
---------------------------------

Directory navigation
********************

.. image:: /_images/2.0/treeview.png
    :width: 392px
    :alt: Directory navigation
    :align: center

The IPython notebook dashboard allows navigation into subdirectories.
URLs are persistent based on the notebook's path and name,
so no more random UUID URLs.

Serving local files no longer needs the ``files/`` prefix.
Relative links across notebooks and other files should work just as if notebooks were regular HTML files.

Interactive widgets
*******************

.. image:: /_images/2.0/widgets.png
    :width: 392px
    :alt: Interactive widgets
    :align: center

IPython 2.0 adds :mod:`IPython.html.widgets`, for manipulating
Python objects in the kernel with GUI controls in the notebook.
IPython comes with a few built-in widgets for simple data types,
and an API designed for developers to build more complex widgets.
See the `widget docs`_ for more information.

.. _widget docs: http://nbviewer.ipython.org/github/ipython/ipython/blob/2.x/examples/Interactive%20Widgets/Index.ipynb


Modal user interface
********************

The notebook has added separate Edit and Command modes,
allowing easier keyboard commands and making keyboard shortcut customization possible.
See the new `User Interface notebook`_ for more information.

.. _User Interface Notebook: http://nbviewer.ipython.org/github/ipython/ipython/blob/2.x/examples/Notebook/User%20Interface.ipynb


You can familiarize yourself with the updated notebook user interface, including an
explanation of Edit and Command modes, by going through the short guided tour
which can be started from the Help menu.

.. image:: /_images/2.0/user-interface.png
    :width: 392px
    :alt: Interface tour
    :align: center


Security
********

2.0 introduces a security model for notebooks,
to prevent untrusted code from executing on users' behalf when notebooks open.
A quick summary of the model:

- Trust is determined by signing notebooks.
- Untrusted HTML output is sanitized.
- Untrusted Javascript is never executed.
- HTML and Javascript in Markdown are never trusted.

Dashboard "Running" tab
***********************

.. image:: /_images/2.0/running-crop.png
    :width: 392px
    :alt: Running tab
    :align: center

The dashboard now has a "Running" tab which shows all of the running notebooks.

Single codebase Python 3 support
--------------------------------

IPython previously supported Python 3 by running 2to3 during setup. We
have now switched to a single codebase which runs natively on Python 2.7
and 3.3.

For notes on how to maintain this, see :doc:`/development/pycompat`.

Selecting matplotlib figure formats
-----------------------------------

Deprecate single-format ``InlineBackend.figure_format``
configurable in favor of ``InlineBackend.figure_formats``,
which is a set, supporting multiple simultaneous figure formats (e.g. png, pdf).

This is available at runtime with the new API function :func:`IPython.display.set_matplotlib_formats`.

clear_output changes
--------------------

* There is no longer a 500ms delay when calling ``clear_output``.
* The ability to clear stderr and stdout individually was removed.
* A new ``wait`` flag that prevents ``clear_output`` from being executed until new
  output is available.  This eliminates animation flickering by allowing the
  user to double buffer the output.
* The output div height is remembered when the ``wait=True`` flag is used.

Extending configurable containers
---------------------------------

Some configurable traits are containers (list, dict, set)
Config objects now support calling ``extend``, ``update``, ``insert``, etc.
on traits in config files, which will ultimately result in calling
those methods on the original object.

The effect being that you can now add to containers without having to copy/paste
the initial value::

    c = get_config()
    c.InlineBackend.rc.update({ 'figure.figsize' : (6, 4) })

Changes to hidden namespace on startup
--------------------------------------

Previously, all names declared in code run at startup
(startup files, ``ipython -i script.py``, etc.)
were added to the hidden namespace, which hides the names from tools like ``%whos``.
There are two changes to this behavior:

1. Scripts run on the command-line ``ipython -i script.py``now behave the same as if they were
   passed to ``%run``, so their variables are never hidden.
2. A boolean config flag ``InteractiveShellApp.hide_initial_ns`` has been added to optionally
   disable the hidden behavior altogether. The default behavior is unchanged.

Using dill to expand serialization support
------------------------------------------

The new function :func:`~IPython.utils.pickleutil.use_dill` allows
dill to extend serialization support in :mod:`IPython.parallel` (closures, etc.).
A :meth:`DirectView.use_dill` convenience method was also added, to enable dill
locally and on all engines with one call.

New IPython console lexer
-------------------------

The IPython console lexer has been rewritten and now supports tracebacks
and customized input/output prompts. See the :ref:`new lexer docs <console_lexer>`
for details.

DisplayFormatter changes
------------------------

There was no official way to query or remove callbacks in the Formatter API.
To remedy this, the following methods are added to :class:`BaseFormatter`:

- ``lookup(instance)`` - return appropriate callback or a given object
- ``lookup_by_type(type_or_str)`` - return appropriate callback for a given type or ``'mod.name'`` type string
- ``pop(type_or_str)`` - remove a type (by type or string).
  Pass a second argument to avoid KeyError (like dict).

All of the above methods raise a KeyError if no match is found.

And the following methods are changed:

- ``for_type(type_or_str)`` - behaves the same as before, only adding support for ``'mod.name'``
  type strings in addition to plain types. This removes the need for ``for_type_by_name()``,
  but it remains for backward compatibility.

Formatters can now raise NotImplementedError in addition to returning None
to indicate that they cannot format a given object.

Exceptions and Warnings
***********************

Exceptions are no longer silenced when formatters fail.
Instead, these are turned into a :class:`~IPython.core.formatters.FormatterWarning`.
A FormatterWarning will also be issued if a formatter returns data of an invalid type
(e.g. an integer for 'image/png').


Other changes
-------------

* `%%capture` cell magic now captures the rich display output, not just
  stdout/stderr

* In notebook, Showing tooltip on tab has been disables to avoid conflict with
  completion, Shift-Tab could still be used to invoke tooltip when inside
  function signature and/or on selection.

* ``object_info_request`` has been replaced by ``object_info`` for consistency in the javascript API.
  ``object_info`` is a simpler interface to register callback that is incompatible with ``object_info_request``.

* Previous versions of IPython on Linux would use the XDG config directory,
  creating :file:`~/.config/ipython` by default. We have decided to go
  back to :file:`~/.ipython` for consistency among systems. IPython will
  issue a warning if it finds the XDG location, and will move it to the new
  location if there isn't already a directory there.

* Equations, images and tables are now centered in Markdown cells.
* Multiline equations are now centered in output areas; single line equations
  remain left justified.

* IPython config objects can be loaded from and serialized to JSON.
  JSON config file have the same base name as their ``.py`` counterpart,
  and will be loaded with higher priority if found.

* bash completion updated with support for all ipython subcommands and flags, including nbconvert

* ``ipython history trim``: added ``--keep=<N>`` as an alias for the more verbose
  ``--HistoryTrim.keep=<N>``
* New ``ipython history clear`` subcommand, which is the same as the newly supported
  ``ipython history trim --keep=0``

* You can now run notebooks in an interactive session via ``%run notebook.ipynb``.

* Print preview is back in the notebook menus, along with options to
  download the open notebook in various formats. This is powered by
  nbconvert.

* :exc:`~IPython.nbconvert.utils.pandoc.PandocMissing` exceptions will be
  raised if Pandoc is unavailable, and warnings will be printed if the version
  found is too old. The recommended Pandoc version for use with nbconvert is
  1.12.1.

* The InlineBackend.figure_format now supports JPEG output if PIL/Pillow is available.

* Input transformers (see :doc:`/config/inputtransforms`) may now raise
  :exc:`SyntaxError` if they determine that input is invalid. The input
  transformation machinery in IPython will handle displaying the exception to
  the user and resetting state.

* Calling ``container.show()`` on javascript display is deprecated and will
  trigger errors on future IPython notebook versions. ``container`` now show
  itself as soon as non-empty

* Added ``InlineBackend.print_figure_kwargs`` to allow passing keyword arguments
  to matplotlib's ``Canvas.print_figure``. This can be used to change the value of
  ``bbox_inches``, which is 'tight' by default, or set the quality of JPEG figures.

* A new callback system has been introduced. For details, see :doc:`/config/callbacks`.

* jQuery and require.js are loaded from CDNs in the default HTML template,
  so javascript is available in static HTML export (e.g. nbviewer).

Backwards incompatible changes
------------------------------

* Python 2.6 and 3.2 are no longer supported: the minimum required
  Python versions are now 2.7 and 3.3.
* The Transformer classes have been renamed to Preprocessor in nbconvert and
  their ``call`` methods have been renamed to ``preprocess``.
* The ``call`` methods of nbconvert post-processors have been renamed to
  ``postprocess``.

* The module ``IPython.core.fakemodule`` has been removed.

* The alias system has been reimplemented to use magic functions. There should be little
  visible difference while automagics are enabled, as they are by default, but parts of the
  :class:`~IPython.core.alias.AliasManager` API have been removed.

* We fixed an issue with switching between matplotlib inline and GUI backends,
  but the fix requires matplotlib 1.1 or newer.  So from now on, we consider
  matplotlib 1.1 to be the minimally supported version for IPython. Older
  versions for the most part will work, but we make no guarantees about it.

* The :command:`pycolor` command has been removed. We recommend the much more capable
  :command:`pygmentize` command from the `Pygments <http://pygments.org/>`_ project.
  If you need to keep the exact output of :command:`pycolor`, you can still use
  ``python -m IPython.utils.PyColorize foo.py``.

* :mod:`IPython.lib.irunner` and its command-line entry point have been removed.
  It had fallen out of use long ago.

* The ``input_prefilter`` hook has been removed, as it was never
  actually used by the code. The input transformer system offers much
  more powerful APIs to work with input code. See
  :doc:`/config/inputtransforms` for details.

* :class:`IPython.core.inputsplitter.IPythonInputSplitter` no longer has a method
  ``source_raw_reset()``, but gains :meth:`~IPython.core.inputsplitter.IPythonInputSplitter.raw_reset`
  instead. Use of ``source_raw_reset`` can be replaced with::

      raw = isp.source_raw
      transformed = isp.source_reset()

* The Azure notebook manager was removed as it was no longer compatible with the notebook storage scheme.

* Simplifying configurable URLs

  - base_project_url is renamed to base_url (base_project_url is kept as a deprecated alias, for now)
  - base_kernel_url configurable is removed (use base_url)
  - websocket_url configurable is removed (use base_url)
.. Developers should add in this file, during each release cycle, information
.. about important changes they've made, in a summary format that's meant for
.. end users.  For each release we normally have three sections: features,  bug
.. fixes and api breakage.
.. Please remember to credit the authors of the contributions by name,
.. especially when they are new users or developers who do not regularly
.. participate  in IPython's development.

.. _whatsnew_index:

=====================
What's new in IPython
=====================

..
    this will appear in the docs if we are not releasing a version (ie if
    `_version_extra` in release.py is an empty string)

.. only:: ipydev

   Development version in-progress features:
   
   .. toctree::

      development


This section documents the changes that have been made in various versions of
IPython. Users should consult these pages to learn about new features, bug
fixes and backwards incompatibilities. Developers should summarize the
development work they do here in a user friendly format.

.. toctree::
   :maxdepth: 1

   version8
   github-stats-8
   version7
   github-stats-7
   version6
   github-stats-6
   version5
   github-stats-5
   version4
   github-stats-4
   version3
   github-stats-3
   version3_widget_migration
   version2.0
   github-stats-2.0
   version1.0
   github-stats-1.0
   version0.13
   github-stats-0.13
   version0.12
   github-stats-0.12
   version0.11
   github-stats-0.11
   version0.10
   version0.9
   version0.8

..
   this makes a hidden toctree that keeps sphinx from complaining about
   documents included nowhere when building docs for stable
   We place it at the end as it will still be reachable via prev/next links.
   
.. only:: ipystable

   .. toctree::
      :hidden:

      development
Migrating Widgets to IPython 3
==============================

Upgrading Notebooks
-------------------

1. The first thing you'll notice when upgrading an IPython 2.0 widget
   notebook to IPython 3.0 is the "Notebook converted" dialog. Click
   "ok".
2. All of the widgets distributed with IPython have been renamed. The
   "Widget" suffix was removed from the end of the class name. i.e.
   ``ButtonWidget`` is now ``Button``.
3. ``ContainerWidget`` was renamed to ``Box``.
4. ``PopupWidget`` was removed from IPython, because bootstrapjs was 
   problematic (creates global variables, etc.). If you use the
   ``PopupWidget``, try using a ``Box`` widget instead. If your notebook
   can't live without the popup functionality, subclass the ``Box``
   widget (both in Python and JS) and use JQuery UI's ``draggable()``
   and ``resizable()`` methods to mimic the behavior.
5. ``add_class`` and ``remove_class`` were removed. More often than not
   a new attribute exists on the widget that allows you to achieve the
   same explicitly. i.e. the ``Button`` widget now has a
   ``button_style`` attribute which you can set to 'primary', 'success',
   'info', 'warning', 'danger', or '' instead of using ``add_class`` to
   add the bootstrap class. ``VBox`` and ``HBox`` classes (flexible
   ``Box`` subclasses) were added that allow you to avoid using
   ``add_class`` and ``remove_class`` to make flexible box model
   layouts. As a last resort, if you can't find a built in attribute for
   the class you want to use, a new ``_dom_classes`` list trait was
   added, which combines ``add_class`` and ``remove_class`` into one
   stateful list.
6. ``set_css`` and ``get_css`` were removed in favor of explicit style
   attributes - ``visible``, ``width``, ``height``, ``padding``,
   ``margin``, ``color``, ``background_color``, ``border_color``,
   ``border_width``, ``border_radius``, ``border_style``,
   ``font_style``, ``font_weight``, ``font_size``, and ``font_family``
   are a few. If you can't find a trait to see the css attribute you
   need, you can, in order of preference, (A) subclass to create your
   own custom widget, (B) use CSS and the ``_dom_classes`` trait to set
   ``_dom_classes``, or (C) use the ``_css`` dictionary to set CSS
   styling like ``set_css`` and ``get_css``.
7. For selection widgets, such as ``Dropdown``, the ``values`` argument
   has been renamed to ``options``.

Upgrading Custom Widgets
------------------------

Javascript
~~~~~~~~~~

1. If you are distributing your widget and decide to use the deferred
   loading technique (preferred), you can remove all references to the
   WidgetManager and the register model/view calls (see the Python
   section below for more information).
2. In 2.0 require.js was used incorrectly, that has been fixed and now
   loading works more like Python's import. Requiring
   ``widgets/js/widget`` doesn't import the ``WidgetManager`` class,
   instead it imports a dictionary that exposes the classes within that
   module:

   .. code:: javascript

       {
       'WidgetModel': WidgetModel,
       'WidgetView': WidgetView,
       'DOMWidgetView': DOMWidgetView,
       'ViewList': ViewList,
       }

   If you decide to continue to use the widget registry (by registering
   your widgets with the manager), you can import a dictionary with a
   handle to the WidgetManager class by requiring
   ``widgets/js/manager``. Doing so will import:

   .. code:: javascript

       {'WidgetManager': WidgetManager}

3. Don't rely on the ``IPython`` namespace for anything. To inherit from
   the DOMWidgetView, WidgetView, or WidgetModel, require
   ``widgets/js/widget`` as ``widget``. If you were inheriting from
   DOMWidgetView, and the code looked like this:

   .. code:: javascript

       IPython.DOMWidgetView.extend({...})

   It would become this:

   .. code:: javascript

       widget.DOMWidgetView.extend({...})

4. Custom models are encouraged. When possible, it's recommended to move
   your code into a custom model, so actions are performed 1 time,
   instead of N times where N is the number of displayed views.

Python
~~~~~~

Generally, custom widget Python code can remain unchanged. If you
distribute your custom widget, you may be using ``display`` and
``Javascript`` to publish the widget's Javascript to the front-end. That
is no longer the recommended way of distributing widget Javascript.
Instead have the user install the Javascript to his/her nbextension
directory or their profile's static directory. Then use the new
``_view_module`` and ``_model_module`` traitlets in combination with
``_view_name`` and ``_model_name`` to instruct require.js on how to load
the widget's Javascript. The Javascript is then loaded when the widget
is used for the first time.

Details
-------

Asynchronous
~~~~~~~~~~~~

In the IPython 2.x series the only way to register custom widget views
and models was to use the registry in the widget manager. Unfortunately,
using this method made distributing and running custom widgets difficult. The widget
maintainer had to either use the rich display framework to push the
widget's Javascript to the notebook or instruct the users to install the
Javascript by hand in a custom profile. With the first method, the
maintainer would have to be careful about when the Javascript was pushed
to the front-end. If the Javascript was pushed on Python widget
``import``, the widgets wouldn't work after page refresh. This is
because refreshing the page does not restart the kernel, and the Python
``import`` statement only runs once in a given kernel instance (unless
you reload the Python modules, which isn't straight forward). This meant
the maintainer would have to have a separate ``push_js()`` method that
the user would have to call after importing the widget's Python code.

Our solution was to add support for loading widget views and models
using require.js paths. Thus the comm and widget frameworks now support
lazy loading. To do so, everything had to be converted to asynchronous
code. HTML5 promises are used to accomplish that
(`#6818 <https://github.com/ipython/ipython/pull/6818>`__,
`#6914 <https://github.com/ipython/ipython/pull/6914>`__).

Symmetry
~~~~~~~~

In IPython 3.0, widgets can be instantiated from the front-end
(`#6664 <https://github.com/ipython/ipython/pull/6664>`__). On top of
this, a widget persistence API was added
(`#7163 <https://github.com/ipython/ipython/pull/7163>`__,
`#7227 <https://github.com/ipython/ipython/pull/7227>`__). With the
widget persistence API, you can persist your widget instances using
Javascript. This makes it easy to persist your widgets to your notebook
document (with a small amount of custom JS). By default, the widgets are
persisted to your web browsers local storage which makes them reappear
when your refresh the page.

Smaller Changes
~~~~~~~~~~~~~~~

-  Latex math is supported in widget ``description``\ s
   (`#5937 <https://github.com/ipython/ipython/pull/5937>`__).
-  Widgets can be display more than once within a single container
   widget (`#5963 <https://github.com/ipython/ipython/pull/5963>`__,
   `#6990 <https://github.com/ipython/ipython/pull/6990>`__).
-  ``FloatRangeSlider`` and ``IntRangeSlider`` were added
   (`#6050 <https://github.com/ipython/ipython/pull/6050>`__).
-  "Widget" was removed from the ends of all of the widget class names
   (`#6125 <https://github.com/ipython/ipython/pull/6125>`__).
-  ``ContainerWidget`` was renamed to ``Box``
   (`#6125 <https://github.com/ipython/ipython/pull/6125>`__).
-  ``HBox`` and ``VBox`` widgets were added
   (`#6125 <https://github.com/ipython/ipython/pull/6125>`__).
-  ``add\_class`` and ``remove\_class`` were removed in favor of a
   ``_dom_classes`` list
   (`#6235 <https://github.com/ipython/ipython/pull/6235>`__).
-  ``get\_css`` and ``set\_css`` were removed in favor of explicit
   traits for widget styling
   (`#6235 <https://github.com/ipython/ipython/pull/6235>`__).
-  ``jslink`` and ``jsdlink`` were added
   (`#6454 <https://github.com/ipython/ipython/pull/6454>`__,
   `#7468 <https://github.com/ipython/ipython/pull/7468>`__).
-  An ``Output`` widget was added, which allows you to ``print`` and
   ``display`` within widgets
   (`#6670 <https://github.com/ipython/ipython/pull/6670>`__).
-  ``PopupWidget`` was removed
   (`#7341 <https://github.com/ipython/ipython/pull/7341>`__).
-  A visual cue was added for widgets with 'dead' comms
   (`#7227 <https://github.com/ipython/ipython/pull/7227>`__).
-  A ``SelectMultiple`` widget was added (a ``Select`` widget that
   allows multiple things to be selected at once)
   (`#6890 <https://github.com/ipython/ipython/pull/6890>`__).
-  A class was added to help manage children views
   (`#6990 <https://github.com/ipython/ipython/pull/6990>`__).
-  A warning was added that shows on widget import because it's expected
   that the API will change again by IPython 4.0. This warning can be
   suppressed (`#7107 <https://github.com/ipython/ipython/pull/7107>`__,
   `#7200 <https://github.com/ipython/ipython/pull/7200>`__,
   `#7201 <https://github.com/ipython/ipython/pull/7201>`__,
   `#7204 <https://github.com/ipython/ipython/pull/7204>`__).

Comm and Widget PR Index
------------------------

Here is a chronological list of PRs affecting the widget and comm frameworks for IPython 3.0. Note that later PRs may revert changes
made in earlier PRs:

- Add placeholder attribute to text widgets
  `#5652 <https://github.com/ipython/ipython/pull/5652>`__
- Add latex support in widget labels,
  `#5937 <https://github.com/ipython/ipython/pull/5937>`__
- Allow widgets to display more than once within container widgets.
  `#5963 <https://github.com/ipython/ipython/pull/5963>`__
- use require.js,
  `#5980 <https://github.com/ipython/ipython/pull/5980>`__
- Range widgets
  `#6050 <https://github.com/ipython/ipython/pull/6050>`__
- Interact on\_demand option
  `#6051 <https://github.com/ipython/ipython/pull/6051>`__
- Allow text input on slider widgets
  `#6106 <https://github.com/ipython/ipython/pull/6106>`__
- support binary buffers in comm messages
  `#6110 <https://github.com/ipython/ipython/pull/6110>`__
- Embrace the flexible box model in the widgets
  `#6125 <https://github.com/ipython/ipython/pull/6125>`__
- Widget trait serialization
  `#6128 <https://github.com/ipython/ipython/pull/6128>`__
- Make Container widgets take children as the first positional
  argument `#6153 <https://github.com/ipython/ipython/pull/6153>`__
- once-displayed
  `#6168 <https://github.com/ipython/ipython/pull/6168>`__
- Validate slider value, when limits change
  `#6171 <https://github.com/ipython/ipython/pull/6171>`__
- Unregistering comms in Comm Manager
  `#6216 <https://github.com/ipython/ipython/pull/6216>`__
- Add EventfulList and EventfulDict trait types.
  `#6228 <https://github.com/ipython/ipython/pull/6228>`__
- Remove add/remove\_class and set/get\_css.
  `#6235 <https://github.com/ipython/ipython/pull/6235>`__
- avoid unregistering widget model twice
  `#6250 <https://github.com/ipython/ipython/pull/6250>`__
- Widget property lock should compare json states, not python states
  `#6332 <https://github.com/ipython/ipython/pull/6332>`__
- Strip the IPY\_MODEL\_ prefix from widget IDs before referencing
  them. `#6377 <https://github.com/ipython/ipython/pull/6377>`__
- "event" is not defined error in Firefox
  `#6437 <https://github.com/ipython/ipython/pull/6437>`__
- Javascript link
  `#6454 <https://github.com/ipython/ipython/pull/6454>`__
- Bulk update of widget attributes
  `#6463 <https://github.com/ipython/ipython/pull/6463>`__
- Creating a widget registry on the Python side.
  `#6493 <https://github.com/ipython/ipython/pull/6493>`__
- Allow widget views to be loaded from require modules
  `#6494 <https://github.com/ipython/ipython/pull/6494>`__
- Fix Issue #6530
  `#6532 <https://github.com/ipython/ipython/pull/6532>`__
- Make comm manager (mostly) independent of InteractiveShell
  `#6540 <https://github.com/ipython/ipython/pull/6540>`__
- Add semantic classes to top-level containers for single widgets
  `#6609 <https://github.com/ipython/ipython/pull/6609>`__
- Selection Widgets: forcing 'value' to be in 'values'
  `#6617 <https://github.com/ipython/ipython/pull/6617>`__
- Allow widgets to be constructed from Javascript
  `#6664 <https://github.com/ipython/ipython/pull/6664>`__
- Output widget
  `#6670 <https://github.com/ipython/ipython/pull/6670>`__
- Minor change in widgets.less to fix alignment issue
  `#6681 <https://github.com/ipython/ipython/pull/6681>`__
- Make Selection widgets respect values order.
  `#6747 <https://github.com/ipython/ipython/pull/6747>`__
- Widget persistence API
  `#6789 <https://github.com/ipython/ipython/pull/6789>`__
- Add promises to the widget framework.
  `#6818 <https://github.com/ipython/ipython/pull/6818>`__
- SelectMultiple widget
  `#6890 <https://github.com/ipython/ipython/pull/6890>`__
- Tooltip on toggle button
  `#6923 <https://github.com/ipython/ipython/pull/6923>`__
- Allow empty text box \*while typing\* for numeric widgets
  `#6943 <https://github.com/ipython/ipython/pull/6943>`__
- Ignore failure of widget MathJax typesetting
  `#6948 <https://github.com/ipython/ipython/pull/6948>`__
- Refactor the do\_diff and manual child view lists into a separate
  ViewList object
  `#6990 <https://github.com/ipython/ipython/pull/6990>`__
- Add warning to widget namespace import.
  `#7107 <https://github.com/ipython/ipython/pull/7107>`__
- lazy load widgets
  `#7120 <https://github.com/ipython/ipython/pull/7120>`__
- Fix padding of widgets.
  `#7139 <https://github.com/ipython/ipython/pull/7139>`__
- Persist widgets across page refresh
  `#7163 <https://github.com/ipython/ipython/pull/7163>`__
- Make the widget experimental error a real python warning
  `#7200 <https://github.com/ipython/ipython/pull/7200>`__
- Make the widget error message shorter and more understandable.
  `#7201 <https://github.com/ipython/ipython/pull/7201>`__
- Make the widget warning brief and easy to filter
  `#7204 <https://github.com/ipython/ipython/pull/7204>`__
- Add visual cue for widgets with dead comms
  `#7227 <https://github.com/ipython/ipython/pull/7227>`__
- Widget values as positional arguments
  `#7260 <https://github.com/ipython/ipython/pull/7260>`__
- Remove the popup widget
  `#7341 <https://github.com/ipython/ipython/pull/7341>`__
- document and validate link, dlink
  `#7468 <https://github.com/ipython/ipython/pull/7468>`__
- Document interact 5637
  `#7525 <https://github.com/ipython/ipython/pull/7525>`__
- Update some broken examples of using widgets
  `#7547 <https://github.com/ipython/ipython/pull/7547>`__
- Use Output widget with Interact
  `#7554 <https://github.com/ipython/ipython/pull/7554>`__
- don't send empty execute\_result messages
  `#7560 <https://github.com/ipython/ipython/pull/7560>`__
- Validation on the python side
  `#7602 <https://github.com/ipython/ipython/pull/7602>`__
- only show prompt overlay if there's a prompt
  `#7661 <https://github.com/ipython/ipython/pull/7661>`__
- Allow predictate to be used for comparison in selection widgets
  `#7674 <https://github.com/ipython/ipython/pull/7674>`__
- Fix widget view persistence.
  `#7680 <https://github.com/ipython/ipython/pull/7680>`__
- Revert "Use Output widget with Interact"
  `#7703 <https://github.com/ipython/ipython/pull/7703>`__
.. _issues_list_4:

Issues closed in the 4.x development cycle
==========================================


Issues closed in 4.2
--------------------

GitHub stats for 2015/02/02 - 2016/04/20 (since 4.1)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 10 issues and merged 22 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A4.2+>`__

The following 10 authors contributed 27 commits.

* Benjamin Ragan-Kelley
* Carlos Cordoba
* Gökhan Karabulut
* Jonas Rauber
* Matthias Bussonnier
* Paul Ivanov
* Sebastian Bank
* Thomas A Caswell
* Thomas Kluyver
* Vincent Woo


Issues closed in 4.1
--------------------

GitHub stats for 2015/08/12 - 2016/02/02 (since 4.0.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 60 issues and merged 148 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/issues?q=milestone%3A4.1+>`__

The following 52 authors contributed 468 commits.

* Aaron Meurer
* Alexandre Avanian
* Anthony Sottile
* Antony Lee
* Arthur Loder
* Ben Kasel
* Ben Rousch
* Benjamin Ragan-Kelley
* bollwyvl
* Carol Willing
* Christopher Roach
* Douglas La Rocca
* Fairly
* Fernando Perez
* Frank Sachsenheim
* Guillaume DOUMENC
* Gábor Luk
* Hoyt Koepke
* Ivan Timokhin
* Jacob Niehus
* JamshedVesuna
* Jan Schulz
* Jan-Philip Gehrcke
* jc
* Jessica B. Hamrick
* jferrara
* John Bohannon
* John Kirkham
* Jonathan Frederic
* Kyle Kelley
* Lev Givon
* Lilian Besson
* lingxz
* Matthias Bussonnier
* memeplex
* Michael Droettboom
* naught101
* Peter Waller
* Pierre Gerold
* Rémy Léone
* Scott Sanderson
* Shanzhuo Zhang
* Sylvain Corlay
* Tayfun Sen
* Thomas A Caswell
* Thomas Ballinger
* Thomas Kluyver
* Vincent Legoll
* Wouter Bolsterlee
* xconverge
* Yuri Numerov
* Zachary Pincus


Issues closed in 4.0
--------------------


GitHub stats for 2015/02/27 - 2015/08/11 (since 3.0)

These lists are automatically generated, and may be incomplete or contain duplicates.

We closed 35 issues and merged 125 pull requests.
The full list can be seen `on GitHub <https://github.com/ipython/ipython/milestone/21>`__

The following 69 authors contributed 1186 commits.

* Abe Guerra
* Adal Chiriliuc
* Alexander Belopolsky
* Andrew Murray
* Antonio Russo
* Benjamin Ragan-Kelley
* Björn Linse
* Brian Drawert
* chebee7i
* Daniel Rocco
* Donny Winston
* Drekin
* Erik Hvatum
* Fernando Perez
* Francisco de la Peña
* Frazer McLean
* Gareth Elston
* Gert-Ludwig Ingold
* Giuseppe Venturini
* Ian Barfield
* Ivan Pozdeev
* Jakob Gager
* Jan Schulz
* Jason Grout
* Jeff Hussmann
* Jessica B. Hamrick
* Joe Borg
* Joel Nothman
* Johan Forsberg
* Jonathan Frederic
* Justin Tyberg
* Koen van Besien
* Kyle Kelley
* Lorena Pantano
* Lucretiel
* Marin Gilles
* mashenjun
* Mathieu
* Matthias Bussonnier
* Merlijn van Deen
* Mikhail Korobov
* Naveen Nathan
* Nicholas Bollweg
* nottaanibot
* Omer Katz
* onesandzeroes
* Patrick Snape
* patter001
* Peter Parente
* Pietro Battiston
* RickWinter
* Robert Smith
* Ryan Nelson
* Scott Sanderson
* Sebastiaan Mathot
* Sylvain Corlay
* thethomask
* Thomas A Caswell
* Thomas Adriaan Hellinger
* Thomas Kluyver
* Tianhui Michael Li
* tmtabor
* unknown
* Victor Ramirez
* Volker Braun
* Wieland Hoffmann
* Yuval Langer
* Zoltán Vörös
* Élie Michel
=============
 0.13 Series
=============

Release 0.13
============

IPython 0.13 contains several major new features, as well as a large amount of
bug and regression fixes.  The previous version (0.12) was released on December
19 2011, and in this development cycle we had:

- ~6 months of work.
- 373 pull requests merged.
- 742 issues closed (non-pull requests).
- contributions from 62 authors.
- 1760 commits.
- a diff of 114226 lines.

The amount of work included in this release is so large, that we can only cover
here the main highlights; please see our :ref:`detailed release statistics
<issues_list_013>` for links to every issue and pull request closed on GitHub
as well as a full list of individual contributors.


Major Notebook improvements: new user interface and more
--------------------------------------------------------

The IPython Notebook, which has proven since its release to be wildly popular,
has seen a massive amount of work in this release cycle, leading to a
significantly improved user experience as well as many new features.

The first user-visible change is a reorganization of the user interface; the
left panel has been removed and was replaced by a real menu system and a
toolbar with icons.  Both the toolbar and the header above the menu can be
collapsed to leave an unobstructed working area:

.. image:: ../_images/ipy_013_notebook_spectrogram.png
    :width: 460px
    :alt: New user interface for Notebook
    :align: center
    :target: ../_images/ipy_013_notebook_spectrogram.png

The notebook handles very long outputs much better than before (this was a
serious usability issue when running processes that generated massive amounts
of output).  Now, in the presence of outputs longer than ~100 lines, the
notebook will automatically collapse to a scrollable area and the entire left
part of this area controls the display: one click in this area will expand the
output region completely, and a double-click will hide it completely.  This
figure shows both the scrolled and hidden modes:

.. image:: ../_images/ipy_013_notebook_long_out.png
    :width: 460px
    :alt: Scrolling and hiding of long output in the notebook.
    :align: center
    :target: ../_images/ipy_013_notebook_long_out.png

.. note::

   The auto-folding of long outputs is disabled in Firefox due to bugs in its
   scrolling behavior.  See :ghpull:`2047` for details.

Uploading notebooks to the dashboard is now easier: in addition to drag and
drop (which can be finicky sometimes), you can now click on the upload text and
use a regular file dialog box to select notebooks to upload. Furthermore, the
notebook dashboard now auto-refreshes its contents and offers buttons to shut
down any running kernels (:ghpull:`1739`):

.. image:: ../_images/ipy_013_dashboard.png
    :width: 460px
    :alt: Improved dashboard
    :align: center
    :target: ../_images/ipy_013_dashboard.png


Cluster management
~~~~~~~~~~~~~~~~~~

The notebook dashboard can now also start and stop clusters, thanks to a new
tab in the dashboard user interface:

.. image:: ../_images/ipy_013_dashboard_cluster.png
    :width: 460px
    :alt: Cluster management from the notebook dashboard
    :align: center
    :target: ../_images/ipy_013_dashboard_cluster.png

This interface allows, for each profile you have configured, to start and stop
a cluster (and optionally override the default number of engines corresponding
to that configuration).  While this hides all error reporting, once you have a
configuration that you know works smoothly, it is a very convenient interface
for controlling your parallel resources.


New notebook format
~~~~~~~~~~~~~~~~~~~

The notebooks saved now use version 3 of our format, which supports heading
levels as well as the concept of 'raw' text cells that are not rendered as
Markdown.  These will be useful with converters_ we are developing, to pass raw
markup (say LaTeX).  That conversion code is still under heavy development and
not quite ready for prime time, but we welcome help on this front so that we
can merge it for full production use as soon as possible.

.. _converters: https://github.com/ipython/nbconvert

.. note::

   v3 notebooks can *not* be read by older versions of IPython, but we provide
   a `simple script`_ that you can use in case you need to export a v3
   notebook to share with a v2 user.

.. _simple script: https://gist.github.com/1935808


JavaScript refactoring
~~~~~~~~~~~~~~~~~~~~~~
  
All the client-side JavaScript has been decoupled to ease reuse of parts of the
machinery without having to build a full-blown notebook. This will make it much
easier to communicate with an IPython kernel from existing web pages and to
integrate single cells into other sites, without loading the full notebook
document-like UI. :ghpull:`1711`.
    
This refactoring also enables the possibility of writing dynamic javascript
widgets that are returned from Python code and that present an interactive view
to the user, with callbacks in Javascript executing calls to the Kernel.  This
will enable many interactive elements to be added by users in notebooks.

An example of this capability has been provided as a proof of concept in
:file:`examples/widgets` that lets you directly communicate with one or more
parallel engines, acting as a mini-console for parallel debugging and
introspection.

    
Improved tooltips
~~~~~~~~~~~~~~~~~

The object tooltips have gained some new functionality. By pressing tab several
times, you can expand them to see more of a docstring, keep them visible as you
fill in a function's parameters, or transfer the information to the pager at the
bottom of the screen. For the details, look at the example notebook
:file:`01_notebook_introduction.ipynb`.

.. figure:: ../_images/ipy_013_notebook_tooltip.png
    :width: 460px
    :alt: Improved tooltips in the notebook.
    :align: center
    :target: ../_images/ipy_013_notebook_tooltip.png

    The new notebook tooltips.

Other improvements to the Notebook
----------------------------------

These are some other notable small improvements to the notebook, in addition to
many bug fixes and minor changes to add polish and robustness throughout:

* The notebook pager (the area at the bottom) is now Resizable by dragging its
  divider handle, a feature that had been requested many times by just about
  anyone who had used the notebook system.  :ghpull:`1705`.

* It is now possible to open notebooks directly from the command line; for
  example: ``ipython notebook path/`` will automatically set ``path/`` as the
  notebook directory, and ``ipython notebook path/foo.ipynb`` will further
  start with the ``foo.ipynb`` notebook opened.  :ghpull:`1686`.
  
* If a notebook directory is specified with ``--notebook-dir`` (or with the
  corresponding configuration flag ``NotebookManager.notebook_dir``), all
  kernels start in this directory.

* Fix codemirror clearing of cells with ``Ctrl-Z``; :ghpull:`1965`.
  
* Text (markdown) cells now line wrap correctly in the notebook, making them
  much easier to edit :ghpull:`1330`.
  
* PNG and JPEG figures returned from plots can be interactively resized in the
  notebook, by dragging them from their lower left corner. :ghpull:`1832`.

* Clear ``In []`` prompt numbers on "Clear All Output".  For more
  version-control-friendly ``.ipynb`` files, we now strip all prompt numbers
  when doing a "Clear all output".  This reduces the amount of noise in
  commit-to-commit diffs that would otherwise show the (highly variable) prompt
  number changes. :ghpull:`1621`.

* The notebook server now requires *two* consecutive ``Ctrl-C`` within 5
  seconds (or an interactive confirmation) to terminate operation.  This makes
  it less likely that you will accidentally kill a long-running server by
  typing ``Ctrl-C`` in the wrong terminal.  :ghpull:`1609`.

* Using ``Ctrl-S`` (or ``Cmd-S`` on a Mac) actually saves the notebook rather
  than providing the fairly useless browser html save dialog.  :ghpull:`1334`.
  
* Allow accessing local files from the notebook (in urls), by serving any local
  file as the url ``files/<relativepath>``.  This makes it possible to, for
  example, embed local images in a notebook.  :ghpull:`1211`.

      
Cell magics
-----------

We have completely refactored the magic system, finally moving the magic
objects to standalone, independent objects instead of being the mixin class
we'd had since the beginning of IPython (:ghpull:`1732`).  Now, a separate base
class is provided in :class:`IPython.core.magic.Magics` that users can subclass
to create their own magics.  Decorators are also provided to create magics from
simple functions without the need for object orientation.  Please see the
:ref:`magic` docs for further details.

All builtin magics now exist in a few subclasses that group together related
functionality, and the new :mod:`IPython.core.magics` package has been created
to organize this into smaller files.
    
This cleanup was the last major piece of deep refactoring needed from the
original 2001 codebase.
    
We have also introduced a new type of magic function, prefixed with `%%`
instead of `%`, which operates at the whole-cell level.  A cell magic receives
two arguments: the line it is called on (like a line magic) and the body of the
cell below it.
    
Cell magics are most natural in the notebook, but they also work in the
terminal and qt console, with the usual approach of using a blank line to
signal cell termination.
    
For example, to time the execution of several statements::

    %%timeit x = 0   # setup
    for i in range(100000):
        x += i**2

This is particularly useful to integrate code in another language, and cell
magics already exist for shell scripts, Cython, R and Octave. Using ``%%script
/usr/bin/foo``, you can run a cell in any interpreter that accepts code via
stdin.

Another handy cell magic makes it easy to write short text files: ``%%file
~/save/to/here.txt``.

The following cell magics are now included by default; all those that use
special interpreters (Perl, Ruby, bash, etc.) assume you have the requisite
interpreter installed:

* ``%%!``: run cell body with the underlying OS shell; this is similar to
  prefixing every line in the cell with ``!``.
  
* ``%%bash``: run cell body under bash.
  
* ``%%capture``: capture the output of the code in the cell (and stderr as
  well).  Useful to run codes that produce too much output that you don't even
  want scrolled.
  
* ``%%file``: save cell body as a file.
  
* ``%%perl``: run cell body using Perl.
  
* ``%%prun``: run cell body with profiler (cell extension of ``%prun``).
  
* ``%%python3``: run cell body using Python 3.
  
* ``%%ruby``: run cell body using Ruby.
  
* ``%%script``: run cell body with the script specified in the first line.
  
* ``%%sh``: run cell body using sh.
  
* ``%%sx``: run cell with system shell and capture process output (cell
  extension of ``%sx``).
  
* ``%%system``: run cell with system shell (``%%!`` is an alias to this).
  
* ``%%timeit``: time the execution of the cell (extension of ``%timeit``).

This is what some of the script-related magics look like in action:

.. image:: ../_images/ipy_013_notebook_script_cells.png
    :width: 460px
    :alt: Cluster management from the notebook dashboard
    :align: center
    :target: ../_images/ipy_013_notebook_script_cells.png
  
In addition, we have also a number of :ref:`extensions <extensions_overview>`
that provide specialized magics.  These typically require additional software
to run and must be manually loaded via ``%load_ext <extension name>``, but are
extremely useful.  The following extensions are provided:

**Cython magics** (extension ``cythonmagic``)
    This extension provides magics to automatically build and compile Python
    extension modules using the Cython_ language. You must install Cython
    separately, as well as a C compiler, for this to work.  The examples
    directory in the source distribution ships with a full notebook
    demonstrating these capabilities:

.. image:: ../_images/ipy_013_notebook_cythonmagic.png
    :width: 460px
    :alt: Cython magic
    :align: center
    :target: ../_images/ipy_013_notebook_cythonmagic.png

.. _cython: http://cython.org

**Octave magics** (extension ``octavemagic``)
    This extension provides several magics that support calling code written in
    the Octave_ language for numerical computing.  You can execute single-lines
    or whole blocks of Octave code, capture both output and figures inline
    (just like matplotlib plots), and have variables automatically converted
    between the two languages.  To use this extension, you must have Octave
    installed as well as the oct2py_ package.  The examples
    directory in the source distribution ships with a full notebook
    demonstrating these capabilities:

.. image:: ../_images/ipy_013_notebook_octavemagic.png
    :width: 460px
    :alt: Octave magic
    :align: center
    :target: ../_images/ipy_013_notebook_octavemagic.png

.. _octave: http://www.gnu.org/software/octave
.. _oct2py: http://pypi.python.org/pypi/oct2py

**R magics** (extension ``rmagic``)
    This extension provides several magics that support calling code written in
    the R_ language for statistical data analysis.  You can execute
    single-lines or whole blocks of R code, capture both output and figures
    inline (just like matplotlib plots), and have variables automatically
    converted between the two languages.  To use this extension, you must have
    R installed as well as the rpy2_ package that bridges Python and R.  The
    examples directory in the source distribution ships with a full notebook
    demonstrating these capabilities:

.. image:: ../_images/ipy_013_notebook_rmagic.png
    :width: 460px
    :alt: R magic
    :align: center
    :target: ../_images/ipy_013_notebook_rmagic.png

.. _R: http://www.r-project.org
.. _rpy2: http://rpy.sourceforge.net/rpy2.html


Tab completer improvements
--------------------------

Useful tab-completion based on live inspection of objects is one of the most
popular features of IPython. To make this process even more user-friendly, the
completers of both the Qt console and the Notebook have been reworked.

The Qt console comes with a new ncurses-like tab completer, activated by
default, which lets you cycle through the available completions by pressing tab,
or select a completion with the arrow keys (:ghpull:`1851`).

.. figure:: ../_images/ipy_013_qtconsole_completer.png
    :width: 460px
    :alt: ncurses-like completer, with highlighted selection.
    :align: center
    :target: ../_images/ipy_013_qtconsole_completer.png

    The new improved Qt console's ncurses-like completer allows to easily
    navigate thought long list of completions.

In the notebook, completions are now sourced both from object introspection and
analysis of surrounding code, so limited completions can be offered for
variables defined in the current cell, or while the kernel is busy 
(:ghpull:`1711`).


We have implemented a new configurable flag to control tab completion on
modules that provide the ``__all__`` attribute::

  IPCompleter.limit_to__all__= Boolean

This instructs the completer to honor ``__all__`` for the completion.
Specifically, when completing on ``object.<tab>``, if True: only those names
in ``obj.__all__`` will be included.  When False [default]: the ``__all__``
attribute is ignored. :ghpull:`1529`.


Improvements to the Qt console
------------------------------

The Qt console continues to receive improvements and refinements, despite the
fact that it is by now a fairly mature and robust component.  Lots of small
polish has gone into it, here are a few highlights:

* A number of changes were made to the underlying code for easier integration
  into other projects such as Spyder_ (:ghpull:`2007`, :ghpull:`2024`).

* Improved menus with a new Magic menu that is organized by magic groups (this
  was made possible by the reorganization of the magic system
  internals). :ghpull:`1782`.

* Allow for restarting kernels without clearing the qtconsole, while leaving a
  visible indication that the kernel has restarted. :ghpull:`1681`.
  
* Allow the native display of jpeg images in the qtconsole. :ghpull:`1643`.

.. _spyder: https://code.google.com/p/spyderlib


  
Parallel
--------

The parallel tools have been improved and fine-tuned on multiple fronts.  Now,
the creation of an :class:`IPython.parallel.Client` object automatically
activates a line and cell magic function ``px`` that sends its code to all the
engines. Further magics can be easily created with the :meth:`.Client.activate`
method, to conveniently execute code on any subset of engines. :ghpull:`1893`.

The ``%%px`` cell magic can also be given an optional targets argument, as well
as a ``--out`` argument for storing its output.

A new magic has also been added, ``%pxconfig``, that lets you configure various
defaults of the parallel magics.  As usual, type  ``%pxconfig?`` for details.

The exception reporting in parallel contexts has been improved to be easier to
read.  Now, IPython directly reports the remote exceptions without showing any
of the internal execution parts:

.. image::  ../_images/ipy_013_par_tb.png
    :width: 460px
    :alt: Improved parallel exceptions.
    :align: center
    :target: ../_images/ipy_013_par_tb.png

The parallel tools now default to using ``NoDB`` as the storage backend for
intermediate results.  This means that the default usage case will have a
significantly reduced memory footprint, though certain advanced features are
not available with this backend.

The parallel magics now display all output, so you can do parallel plotting or
other actions with complex display.  The ``px`` magic has now both line and cell
modes, and in cell mode finer control has been added about how to collate
output from multiple engines. :ghpull:`1768`.

There have also been incremental improvements to the SSH launchers:
    
* add to_send/fetch steps for moving connection files around.
  
* add SSHProxyEngineSetLauncher, for invoking to `ipcluster engines` on a
  remote host. This can be used to start a set of engines via PBS/SGE/MPI
  *remotely*.
    
This makes the SSHLauncher usable on machines without shared filesystems.

A number of 'sugar' methods/properties were added to AsyncResult that are
quite useful (:ghpull:`1548`) for everday work:
    
    * ``ar.wall_time`` = received - submitted
    * ``ar.serial_time`` = sum of serial computation time
    * ``ar.elapsed`` = time since submission (wall_time if done)
    * ``ar.progress`` = (int) number of sub-tasks that have completed
    * ``len(ar)`` = # of tasks
    * ``ar.wait_interactive()``: prints progress
    
Added :meth:`.Client.spin_thread` / :meth:`~.Client.stop_spin_thread` for
running spin in a background thread, to keep zmq queue clear.  This can be used
to ensure that timing information is as accurate as possible (at the cost of
having a background thread active).

Set TaskScheduler.hwm default to 1 instead of 0.  1 has more
predictable/intuitive behavior, if often slower, and thus a more logical
default.  Users whose workloads require maximum throughput and are largely
homogeneous in time per task can make the optimization themselves, but now the
behavior will be less surprising to new users. :ghpull:`1294`.


Kernel/Engine unification
-------------------------

This is mostly work 'under the hood', but it is actually a *major* achievement
for the project that has deep implications in the long term: at last, we have
unified the main object that executes as the user's interactive shell (which we
refer to as the *IPython kernel*) with the objects that run in all the worker
nodes of the parallel computing facilities (the *IPython engines*).  Ever since
the first implementation of IPython's parallel code back in 2006, we had wanted
to have these two roles be played by the same machinery, but a number of
technical reasons had prevented that from being true.

In this release we have now merged them, and this has a number of important
consequences:

* It is now possible to connect any of our clients (qtconsole or terminal
  console) to any individual parallel engine, with the *exact* behavior of
  working at a 'regular' IPython console/qtconsole.  This makes debugging,
  plotting, etc. in parallel scenarios vastly easier.

* Parallel engines can always execute arbitrary 'IPython code', that is, code
  that has magics, shell extensions, etc.  In combination with the ``%%px``
  magics, it is thus extremely natural for example to send to all engines a
  block of Cython or R code to be executed via the new Cython and R magics. For
  example, this snippet would send the R block to all active engines in a
  cluster::

    %%px
    %%R
    ... R code goes here
  
* It is possible to embed not only an interactive shell with the
  :func:`IPython.embed` call as always, but now you can also embed a *kernel*
  with :func:`IPython.embed_kernel()`.  Embedding an IPython kernel in an
  application is useful when you want to use :func:`IPython.embed` but don't
  have a terminal attached on stdin and stdout.

* The new :func:`IPython.parallel.bind_kernel` allows you to promote Engines to
  listening Kernels, and connect QtConsoles to an Engine and debug it
  directly.

In addition, having a single core object through our entire architecture also
makes the project conceptually cleaner, easier to maintain and more robust.
This took a lot of work to get in place, but we are thrilled to have this major
piece of architecture finally where we'd always wanted it to be.
  

Official Public API
-------------------

We have begun organizing our API for easier public use, with an eye towards an
official IPython 1.0 release which will firmly maintain this API compatible for
its entire lifecycle.  There is now an :mod:`IPython.display` module that
aggregates all display routines, and the :mod:`traitlets.config` namespace has
all public configuration tools.  We will continue improving our public API
layout so that users only need to import names one level deeper than the main
``IPython`` package to access all public namespaces.


IPython notebook file icons
---------------------------

The directory ``docs/resources`` in the source distribution contains SVG and
PNG versions of our file icons, as well as an ``Info.plist.example`` file with
instructions to install them on Mac OSX.  This is a first draft of our icons,
and we encourage contributions from users with graphic talent to improve them
in the future.

	  
New top-level `locate` command
------------------------------

Add `locate` entry points; these would be useful for quickly locating IPython
directories and profiles from other (non-Python) applications. :ghpull:`1762`.
    
Examples::
    
    $> ipython locate
    /Users/me/.ipython
  
    $> ipython locate profile foo
    /Users/me/.ipython/profile_foo
  
    $> ipython locate profile
    /Users/me/.ipython/profile_default
  
    $> ipython locate profile dne
    [ProfileLocate] Profile u'dne' not found.

	
Other new features and improvements
-----------------------------------

* **%install_ext**: A new magic function to install an IPython extension from
  a URL. E.g. ``%install_ext
  https://bitbucket.org/birkenfeld/ipython-physics/raw/default/physics.py``.

* The ``%loadpy`` magic is no longer restricted to Python files, and has been
  renamed ``%load``. The old name remains as an alias.

* New command line arguments will help external programs find IPython folders:
  ``ipython locate`` finds the user's IPython directory, and ``ipython locate
  profile foo`` finds the folder for the 'foo' profile (if it exists).

* The :envvar:`IPYTHON_DIR` environment variable, introduced in the Great
  Reorganization of 0.11 and existing only in versions 0.11-0.13, has been
  deprecated. As described in :ghpull:`1167`, the complexity and confusion of
  migrating to this variable is not worth the aesthetic improvement. Please use
  the historical :envvar:`IPYTHONDIR` environment variable instead.

* The default value of *interactivity* passed from
  :meth:`~IPython.core.interactiveshell.InteractiveShell.run_cell` to
  :meth:`~IPython.core.interactiveshell.InteractiveShell.run_ast_nodes`
  is now configurable.

* New ``%alias_magic`` function to conveniently create aliases of existing
  magics, if you prefer to have shorter names for personal use.

* We ship unminified versions of the JavaScript libraries we use, to better
  comply with Debian's packaging policies.

* Simplify the information presented by ``obj?/obj??`` to eliminate a few
  redundant fields when possible.  :ghpull:`2038`.

* Improved continuous integration for IPython.  We now have automated test runs
  on `Shining Panda <https://jenkins.shiningpanda.com/ipython>`_ and `Travis-CI
  <http://travis-ci.org/#!/ipython/ipython>`_, as well as `Tox support
  <http://tox.testrun.org>`_.

* The `vim-ipython`_ functionality (externally developed) has been updated to
  the latest version.

.. _vim-ipython: https://github.com/ivanov/vim-ipython

* The ``%save`` magic now has a ``-f`` flag to force overwriting, which makes
  it much more usable in the notebook where it is not possible to reply to
  interactive questions from the kernel. :ghpull:`1937`.

* Use dvipng to format sympy.Matrix, enabling display of matrices in the Qt
  console with the sympy printing extension. :ghpull:`1861`.

* Our messaging protocol now has a reasonable test suite, helping ensure that
  we don't accidentally deviate from the spec and possibly break third-party
  applications that may have been using it.  We encourage users to contribute
  more stringent tests to this part of the test suite.  :ghpull:`1627`.

* Use LaTeX to display, on output, various built-in types with the SymPy
  printing extension. :ghpull:`1399`.

* Add Gtk3 event loop integration and example. :ghpull:`1588`.

* ``clear_output`` improvements, which allow things like progress bars and other
  simple animations to work well in the notebook (:ghpull:`1563`):
    
    * `clear_output()` clears the line, even in terminal IPython, the QtConsole
      and plain Python as well, by printing `\r` to streams.
    
    * `clear_output()` avoids the flicker in the notebook by adding a delay,
      and firing immediately upon the next actual display message.
    
    * `display_javascript` hides its `output_area` element, so using display to
      run a bunch of javascript doesn't result in ever-growing vertical space.

* Add simple support for running inside a virtualenv.  While this doesn't
  supplant proper installation (as users should do), it helps ad-hoc calling of
  IPython from inside a virtualenv. :ghpull:`1388`.

  
Major Bugs fixed
----------------

In this cycle, we have :ref:`closed over 740 issues <issues_list_013>`, but a
few major ones merit special mention:

* The ``%pastebin`` magic has been updated to point to gist.github.com, since
  unfortunately http://paste.pocoo.org has closed down. We also added a -d flag
  for the user to provide a gist description string. :ghpull:`1670`.

* Fix ``%paste`` that would reject certain valid inputs. :ghpull:`1258`.

* Fix sending and receiving of Numpy structured arrays (those with composite
  dtypes, often used as recarrays). :ghpull:`2034`.

* Reconnect when the websocket connection closes unexpectedly. :ghpull:`1577`.

* Fix truncated representation of objects in the debugger by showing at least
  80 characters' worth of information.  :ghpull:`1793`.

* Fix logger to be Unicode-aware: logging could crash ipython if there was
  unicode in the input. :ghpull:`1792`.

* Fix images missing from XML/SVG export in the Qt console. :ghpull:`1449`.

* Fix deepreload on Python 3. :ghpull:`1625`, as well as having a much cleaner
  and more robust implementation of deepreload in general. :ghpull:`1457`.


Backwards incompatible changes
------------------------------

* The exception :exc:`IPython.core.error.TryNext` previously accepted
  arguments and keyword arguments to be passed to the next implementation
  of the hook. This feature was removed as it made error message propagation
  difficult and violated the principle of loose coupling.
============
 3.x Series
============

IPython 3.2.3
=============

Fixes compatibility with Python 3.4.4.

IPython 3.2.2
=============

Address vulnerabilities when files have maliciously crafted filenames (CVE-2015-6938),
or vulnerability when opening text files with malicious binary content (CVE pending).

Users are **strongly** encouraged to upgrade immediately.
There are also a few small unicode and nbconvert-related fixes.


IPython 3.2.1
=============

IPython 3.2.1 is a small bugfix release, primarily for cross-site security fixes in the notebook.
Users are **strongly** encouraged to upgrade immediately.
There are also a few small unicode and nbconvert-related fixes.

See :ref:`issues_list_3` for details.


IPython 3.2
===========

IPython 3.2 contains important security fixes. Users are **strongly** encouraged to upgrade immediately.

Highlights:

- Address cross-site scripting vulnerabilities CVE-2015-4706, CVE-2015-4707
- A security improvement that set the secure attribute to login cookie to prevent them to be sent over http
- Revert the face color of matplotlib axes in the inline backend to not be transparent.
- Enable mathjax safe mode by default
- Fix XSS vulnerability in JSON error messages
- Various widget-related fixes

See :ref:`issues_list_3` for details.


IPython 3.1
===========

Released April 3, 2015

The first 3.x bugfix release, with 33 contributors and 344 commits.
This primarily includes bugfixes to notebook layout and focus problems.


Highlights:

- Various focus jumping and scrolling fixes in the notebook.
- Various message ordering and widget fixes in the notebook.
- Images in markdown and output are confined to the notebook width.
  An `.unconfined` CSS class is added to disable this behavior per-image.
  The resize handle on output images is removed.
- Improved ordering of tooltip content for Python functions, putting the signature at the top.
- Fix UnicodeErrors when displaying some objects with unicode reprs on Python 2.
- Set the kernel's working directory to the notebook directory when running ``nbconvert --execute``,
  so that behavior matches the live notebook.
- Allow setting custom SSL options for the tornado server with ``NotebookApp.ssl_options``,
  and protect against POODLE with default settings by disabling SSLv3.
- Fix memory leak in the IPython.parallel Controller on Python 3.


See :ref:`issues_list_3` for details.


Release 3.0
===========

Released February 27, 2015

This is a really big release. Over 150 contributors, and almost 6000 commits in a bit under a year.
Support for languages other than Python is greatly improved,
notebook UI has been significantly redesigned,
and a lot of improvement has happened in the experimental interactive widgets.
The message protocol and document format have both been updated,
while maintaining better compatibility with previous versions than prior updates.
The notebook webapp now enables editing of any text file, and even
a web-based terminal (on Unix platforms).

3.x will be the last monolithic release of IPython,
as the next release cycle will see the growing project split into its Python-specific and language-agnostic components.
Language-agnostic projects (notebook, qtconsole, etc.) will move under the umbrella of the new Project Jupyter name,
while Python-specific projects (interactive Python shell, Python kernel, IPython.parallel)
will remain under IPython, and be split into a few smaller packages.
To reflect this, IPython is in a bit of a transition state.
The logo on the notebook is now the Jupyter logo.
When installing kernels system-wide, they go in a `jupyter` directory.
We are going to do our best to ease this transition for users and developers.

Big changes are ahead.


Using different kernels
-----------------------

.. image:: ../_images/kernel_selector_screenshot.png
   :alt: Screenshot of 'new' dropdown showing different kernels
   :align: center

You can now choose a kernel for a notebook within the user interface, rather
than starting up a separate notebook server for each kernel you want to use. The
syntax highlighting adapts to match the language you're working in.

Information about the kernel is stored in the notebook file, so when you open a
notebook, it will automatically start the correct kernel.

It is also easier to use the Qt console and the terminal console with other
kernels, using the --kernel flag::

    ipython qtconsole --kernel bash
    ipython console --kernel bash

    # To list available kernels
    ipython kernelspec list

Kernel authors should see :ref:`kernelspecs` for how to register their kernels
with IPython so that these mechanisms work.

Typing unicode identifiers
--------------------------

.. image:: /_images/unicode_completion.png

Complex expressions can be much cleaner when written with a wider choice of
characters. Python 3 allows unicode identifiers, and IPython 3 makes it easier
to type those, using a feature from Julia. Type a backslash followed by a LaTeX
style short name, such as ``\alpha``. Press tab, and it will turn into α.

Widget migration guide
----------------------
The widget framework has a lot of backwards incompatible changes.
For information about migrating widget notebooks and custom widgets to 3.0 refer
to the :doc:`widget migration guide<version3_widget_migration>`.

Other new features
------------------

* :class:`~.TextWidget` and :class:`~.TextareaWidget` objects now include a
  ``placeholder`` attribute, for displaying placeholder text before the
  user has typed anything.

* The :magic:`load` magic can now find the source for objects in the user namespace.
  To enable searching the namespace, use the ``-n`` option.

  .. sourcecode:: ipython

      In [1]: %load -n my_module.some_function

* :class:`~.DirectView` objects have a new :meth:`~.DirectView.use_cloudpickle`
  method, which works like ``view.use_dill()``, but causes the ``cloudpickle``
  module from PiCloud's `cloud`__ library to be used rather than dill or the
  builtin pickle module.

  __ https://pypi.python.org/pypi/cloud

* Added a .ipynb exporter to nbconvert.  It can be used by passing `--to notebook`
  as a commandline argument to nbconvert.

* New nbconvert preprocessor called :class:`~.ClearOutputPreprocessor`. This
  clears the output from IPython notebooks.

* New preprocessor for nbconvert that executes all the code cells in a notebook.
  To run a notebook and save its output in a new notebook::

      ipython nbconvert InputNotebook --ExecutePreprocessor.enabled=True --to notebook --output Executed

* Consecutive stream (stdout/stderr) output is merged into a single output
  in the notebook document.
  Previously, all output messages were preserved as separate output fields in the JSON.
  Now, the same merge is applied to the stored output as the displayed output,
  improving document load time for notebooks with many small outputs.

* ``NotebookApp.webapp_settings`` is deprecated and replaced with
  the more informatively named ``NotebookApp.tornado_settings``.

* Using :magic:`timeit` prints warnings if there is at least a 4x difference in timings
  between the slowest and fastest runs, since this might meant that the multiple
  runs are not independent of one another.

* It's now possible to provide mechanisms to integrate IPython with other event
  loops, in addition to the ones we already support. This lets you run GUI code
  in IPython with an interactive prompt, and to embed the IPython
  kernel in GUI applications. See :doc:`/config/eventloops` for details. As part
  of this, the direct ``enable_*`` and ``disable_*`` functions for various GUIs
  in :mod:`IPython.lib.inputhook` have been deprecated in favour of
  :meth:`~.InputHookManager.enable_gui` and :meth:`~.InputHookManager.disable_gui`.

* A ``ScrollManager`` was added to the notebook.  The ``ScrollManager`` controls how the notebook document is scrolled using keyboard.  Users can inherit from the ``ScrollManager`` or ``TargetScrollManager`` to customize how their notebook scrolls.  The default ``ScrollManager`` is the ``SlideScrollManager``, which tries to scroll to the nearest slide or sub-slide cell.

* The function :func:`~IPython.html.widgets.interaction.interact_manual` has been
  added which behaves similarly to :func:`~IPython.html.widgets.interaction.interact`,
  but adds a button to explicitly run the interacted-with function, rather than
  doing it automatically for every change of the parameter widgets. This should
  be useful for long-running functions.

* The ``%cython`` magic is now part of the Cython module. Use `%load_ext Cython` with a version of Cython >= 0.21 to have access to the magic now.

* The Notebook application now offers integrated terminals on Unix platforms,
  intended for when it is used on a remote server. To enable these, install
  the ``terminado`` Python package.

* The Notebook application can now edit any plain text files, via a full-page CodeMirror instance.

* Setting the default highlighting language for nbconvert with the config option
  ``NbConvertBase.default_language`` is deprecated. Nbconvert now respects
  metadata stored in the :ref:`kernel spec <kernelspecs>`.

* IPython can now be configured systemwide, with files in :file:`/etc/ipython`
  or :file:`/usr/local/etc/ipython` on Unix systems,
  or :file:`{%PROGRAMDATA%}\\ipython` on Windows.

* Added support for configurable user-supplied `Jinja
  <http://jinja.pocoo.org/>`_ HTML templates for the notebook.  Paths to
  directories containing template files can be specified via
  ``NotebookApp.extra_template_paths``.  User-supplied template directories
  searched first by the notebook, making it possible to replace existing
  templates with your own files.

  For example, to replace the notebook's built-in ``error.html`` with your own,
  create a directory like ``/home/my_templates`` and put your override template
  at ``/home/my_templates/error.html``.  To start the notebook with your custom
  error page enabled, you would run::

      ipython notebook '--extra_template_paths=["/home/my_templates/"]'

  It's also possible to override a template while also `inheriting
  <http://jinja.pocoo.org/docs/dev/templates/#template-inheritance>`_ from that
  template, by prepending ``templates/`` to the ``{% extends %}`` target of
  your child template.  This is useful when you only want to override a
  specific block of a template.  For example, to add additional CSS to the
  built-in ``error.html``, you might create an override that looks like::

    {% extends "templates/error.html" %}

    {% block stylesheet %}
    {{super()}}
    <style type="text/css">
      /* My Awesome CSS */
    </style>
    {% endblock %}

* Added a widget persistence API.  This allows you to persist your notebooks interactive widgets.
  Two levels of control are provided:
  1. Higher level- ``WidgetManager.set_state_callbacks`` allows you to register callbacks for loading and saving widget state.  The callbacks you register are automatically called when necessary.
  2. Lower level- the ``WidgetManager`` Javascript class now has ``get_state`` and ``set_state`` methods that allow you to get and set the state of the widget runtime.

  Example code for persisting your widget state to session data::

    %%javascript
    require(['widgets/js/manager'], function(manager) {
        manager.WidgetManager.set_state_callbacks(function() { // Load
            return JSON.parse(sessionStorage.widgets_state || '{}');
        }, function(state) { // Save
            sessionStorage.widgets_state = JSON.stringify(state);
        });
    });

* Enhanced support for :magic:`env` magic.  As before, :magic:`env` with no
  arguments displays all environment variables and values.  Additionally,
  :magic:`env` can be used to get or set individual environment variables. To
  display an individual value, use the `%env var` syntax. To set a value, use
  `env var val` or `env var=val`. Python value expansion using `$` works as usual.


Backwards incompatible changes
------------------------------

* The :ref:`message protocol <messaging>` has been updated from version 4 to version 5.
  Adapters are included, so IPython frontends can still talk to kernels that
  implement protocol version 4.

* The notebook format has been updated from version 3 to version 4.
  Read-only support for v4 notebooks has been backported to IPython 2.4.
  Notable changes:
  
  * heading cells are removed in favor or markdown headings
  * notebook outputs and output messages are more consistent with each other
  * use :func:`IPython.nbformat.read` and :func:`~IPython.nbformat.write`
    to read and write notebook files
    instead of the deprecated :mod:`IPython.nbformat.current` APIs.

  You can downgrade a notebook to v3 via ``nbconvert``::
  
      ipython nbconvert --to notebook --nbformat 3 <notebook>
  
  which will create :file:`notebook.v3.ipynb`, a copy of the notebook in v3 format.

* :func:`IPython.core.oinspect.getsource` call specification has changed:

  * `oname` keyword argument has been added for property source formatting
  * `is_binary` keyword argument has been dropped, passing ``True`` had
    previously short-circuited the function to return ``None`` unconditionally

* Removed the octavemagic extension: it is now available as ``oct2py.ipython``.

* Creating PDFs with LaTeX no longer uses a post processor.
  Use `nbconvert --to pdf` instead of `nbconvert --to latex --post pdf`.

* Used https://github.com/jdfreder/bootstrap2to3 to migrate the Notebook to Bootstrap 3.

  Additional changes:

  - Set `.tab-content .row` `0px;` left and right margin (bootstrap default is `-15px;`)
  - Removed `height: @btn_mini_height;` from `.list_header>div, .list_item>div` in `tree.less`
  - Set `#header` div `margin-bottom: 0px;`
  - Set `#menus` to `float: left;`
  - Set `#maintoolbar .navbar-text` to `float: none;`
  - Added no-padding convenience class.
  - Set border of #maintoolbar to 0px

* Accessing the `container` DOM object when displaying javascript has been
  deprecated in IPython 2.0 in favor of accessing `element`. Starting with
  IPython 3.0 trying to access `container` will raise an error in browser
  javascript console.

* ``IPython.utils.py3compat.open`` was removed: :func:`io.open` provides all
  the same functionality.

* The NotebookManager and ``/api/notebooks`` service has been replaced by
  a more generic ContentsManager and ``/api/contents`` service,
  which supports all kinds of files.
* The Dashboard now lists all files, not just notebooks and directories.
* The ``--script`` hook for saving notebooks to Python scripts is removed,
  use :samp:`ipython nbconvert --to python {notebook}` instead.

* The ``rmagic`` extension is deprecated, as it is now part of rpy2. See
  :mod:`rpy2.ipython.rmagic`.

* :meth:`~.KernelManager.start_kernel` and :meth:`~.KernelManager.format_kernel_cmd`
  no longer accept a ``executable`` parameter. Use the kernelspec machinery instead.

* The widget classes have been renamed from `*Widget` to `*`.  The old names are
  still functional, but are deprecated.  i.e. `IntSliderWidget` has been renamed
  to `IntSlider`.
* The ContainerWidget was renamed to Box and no longer defaults as a flexible
  box in the web browser.  A new FlexBox widget was added, which allows you to
  use the flexible box model.

* The notebook now uses a single websocket at `/kernels/<kernel-id>/channels` instead of separate
  `/kernels/<kernel-id>/{shell|iopub|stdin}` channels. Messages on each channel are identified by a
  `channel` key in the message dict, for both send and recv.


Content Security Policy
```````````````````````

The Content Security Policy is a web standard for adding a layer of security to
detect and mitigate certain classes of attacks, including Cross Site Scripting
(XSS) and data injection attacks. This was introduced into the notebook to
ensure that the IPython Notebook and its APIs (by default) can only be embedded
in an iframe on the same origin.

Override ``headers['Content-Security-Policy']`` within your notebook
configuration to extend for alternate domains and security settings.::

    c.NotebookApp.tornado_settings = {
        'headers': {
            'Content-Security-Policy': "frame-ancestors 'self'"
        }
    }

Example policies::

    Content-Security-Policy: default-src 'self' https://*.jupyter.org

Matches embeddings on any subdomain of jupyter.org, so long as they are served
over SSL.

There is a `report-uri <https://developer.mozilla.org/en-US/docs/Web/HTTP/Headers/Content-Security-Policy/report-uri>`_ endpoint available for logging CSP violations, located at
``/api/security/csp-report``. To use it, set ``report-uri`` as part of the CSP::

    c.NotebookApp.tornado_settings = {
        'headers': {
            'Content-Security-Policy': "frame-ancestors 'self'; report-uri /api/security/csp-report"
        }
    }

It simply provides the CSP report as a warning in IPython's logs. The default
CSP sets this report-uri relative to the ``base_url`` (not shown above).

For a more thorough and accurate guide on Content Security Policies, check out
`MDN's Using Content Security Policy <https://developer.mozilla.org/en-US/docs/Web/Security/CSP/Using_Content_Security_Policy>`_ for more examples.
Incompatible change switch to perl
----------------------------------

Document which filename start with ``incompat-`` will be gathers in their own
incompatibility section.

Starting with IPython 42, only perl code execution is allowed. See :ghpull:`42`
Antigravity feature
===================

Example new antigravity feature. Try ``import antigravity`` in a Python 3
console.
.. _execution_semantics:

Execution of cells in the IPython kernel
========================================

When IPython kernel receives `execute_request <https://jupyter-client.readthedocs.io/en/latest/messaging.html#execute>`_
with user code, it processes the message in the following phases:

1. Fire the ``pre_execute`` event.
2. Fire the ``pre_run_cell`` event unless silent is ``True``.
3. Execute ``run_cell`` method to preprocess ``code``, compile and run it, see below for details.
4. If execution succeeds, expressions in ``user_expressions`` are computed.
   This ensures that any error in the expressions don't affect the main code execution.
5. Fire the ``post_execute`` event.
6. Fire the ``post_run_cell`` event unless silent is ``True``.

.. seealso::

    :doc:`/config/callbacks`


Running user ``code``
=====================

First, the ``code`` cell is transformed to expand ``%magic`` and ``!system``
commands by ``IPython.core.inputtransformer2``. Then expanded cell is compiled
using standard Python :func:`compile` function and executed.

Python :func:`compile` function provides ``mode`` argument to select one
of three ways of compiling code:

*single*
  Valid for a single interactive statement (though the source can contain
  multiple lines, such as a for loop).  When compiled in this mode, the
  generated bytecode contains special instructions that trigger the calling of
  :func:`sys.displayhook` for any expression in the block that returns a value.
  This means that a single statement can actually produce multiple calls to
  :func:`sys.displayhook`, if for example it contains a loop where each
  iteration computes an unassigned expression would generate 10 calls::

      for i in range(10):
          i**2

*exec*
  An arbitrary amount of source code, this is how modules are compiled.
  :func:`sys.displayhook` is *never* implicitly called.

*eval*
  A single expression that returns a value.  :func:`sys.displayhook` is *never*
  implicitly called.


The ``code`` field is split into individual blocks each of which is valid for
execution in 'single' mode, and then:

- If there is only a single block: it is executed in 'single' mode.

- If there is more than one block:

  * if the last block is a single line long, run all but the last in 'exec' mode
    and the very last one in 'single' mode.  This makes it easy to type simple
    expressions at the end to see computed values.

  * if the last block is no more than two lines long, run all but the last in
    'exec' mode and the very last one in 'single' mode.  This makes it easy to
    type simple expressions at the end to see computed values.  - otherwise
    (last one is also multiline), run all in 'exec' mode

  * otherwise (last block is also multiline), run all in 'exec' mode as a single
    unit.


.. _config_overview:

============================================
Overview of the IPython configuration system
============================================

This section describes the IPython configuration system. This is based on
:mod:`traitlets.config`; see that documentation for more information
about the overall architecture.

Configuration file location
===========================

So where should you put your configuration files? IPython uses "profiles" for
configuration, and by default, all profiles will be stored in the so called
"IPython directory". The location of this directory is determined by the
following algorithm:

* If the ``ipython-dir`` command line flag is given, its value is used.

* If not, the value returned by :func:`IPython.paths.get_ipython_dir`
  is used. This function will first look at the :envvar:`IPYTHONDIR`
  environment variable and then default to :file:`~/.ipython`.
  Historical support for the :envvar:`IPYTHON_DIR` environment variable will
  be removed in a future release.

For most users, the configuration directory will be :file:`~/.ipython`.

Previous versions of IPython on Linux would use the XDG config directory,
creating :file:`~/.config/ipython` by default. We have decided to go
back to :file:`~/.ipython` for consistency among systems. IPython will
issue a warning if it finds the XDG location, and will move it to the new
location if there isn't already a directory there.

Once the location of the IPython directory has been determined, you need to know
which profile you are using. For users with a single configuration, this will
simply be 'default', and will be located in
:file:`<IPYTHONDIR>/profile_default`.

The next thing you need to know is what to call your configuration file. The
basic idea is that each application has its own default configuration filename.
The default named used by the :command:`ipython` command line program is
:file:`ipython_config.py`, and *all* IPython applications will use this file.
The IPython kernel will load its own config file *after*
:file:`ipython_config.py`. To load a particular configuration file instead of
the default, the name can be overridden by the ``config_file`` command line
flag.

To generate the default configuration files, do::

    $ ipython profile create

and you will have a default :file:`ipython_config.py` in your IPython directory
under :file:`profile_default`.

.. note::

    IPython configuration options are case sensitive, and IPython cannot
    catch misnamed keys or invalid values.
    
    By default IPython will also ignore any invalid configuration files. 

.. versionadded:: 5.0

    IPython can be configured to abort in case of invalid configuration file.
    To do so set the environment variable ``IPYTHON_SUPPRESS_CONFIG_ERRORS`` to
    `'1'` or `'true'`


Locating these files
--------------------

From the command-line, you can quickly locate the IPYTHONDIR or a specific
profile with:

.. sourcecode:: bash

    $ ipython locate
    /home/you/.ipython
    
    $ ipython locate profile foo
    /home/you/.ipython/profile_foo

These map to the utility functions: :func:`IPython.paths.get_ipython_dir`
and :func:`IPython.paths.locate_profile` respectively.


.. _profiles_dev:

Profiles
========

A profile is a directory containing configuration and runtime files, such as
logs, connection info for the parallel apps, and your IPython command history.

The idea is that users often want to maintain a set of configuration files for
different purposes: one for doing numerical computing with NumPy and SciPy and
another for doing symbolic computing with SymPy. Profiles make it easy to keep a
separate configuration files, logs, and histories for each of these purposes.

Let's start by showing how a profile is used:

.. code-block:: bash

    $ ipython --profile=sympy

This tells the :command:`ipython` command line program to get its configuration
from the "sympy" profile. The file names for various profiles do not change. The
only difference is that profiles are named in a special way. In the case above,
the "sympy" profile means looking for :file:`ipython_config.py` in :file:`<IPYTHONDIR>/profile_sympy`.

The general pattern is this: simply create a new profile with:

.. code-block:: bash

    $ ipython profile create <name>

which adds a directory called ``profile_<name>`` to your IPython directory. Then
you can load this profile by adding ``--profile=<name>`` to your command line
options. Profiles are supported by all IPython applications.

IPython extends the config loader for Python files so that you can inherit
config from another profile. To do this, use a line like this in your Python
config file:

.. sourcecode:: python

    load_subconfig('ipython_config.py', profile='default')
:orphan:

Writing code for Python 2 and 3
===============================

.. module:: IPython.utils.py3compat
   :synopsis: Python 2 & 3 compatibility helpers


IPython 6 requires Python 3, so our compatibility module
``IPython.utils.py3compat`` is deprecated and will be removed in a future
version. In most cases, we recommend you use the `six module
<https://six.readthedocs.io/>`__ to support compatible code. This is widely
used by other projects, so it is familiar to many developers and thoroughly
battle-tested.

Our ``py3compat`` module provided some more specific unicode conversions than
those offered by ``six``. If you want to use these, copy them into your own code
from IPython 5.x. Do not rely on importing them from IPython, as the module may
be removed in the future.

.. seealso::

   `Porting Python 2 code to Python 3 <https://docs.python.org/3/howto/pyporting.html>`_
     Official information in the Python docs.

   `Python-Modernize <https://python-modernize.readthedocs.io/en/latest/>`_
     A tool which helps make code compatible with Python 3.

   `Python-Future <https://python-future.org/>`_
     Another compatibility tool, which focuses on writing code for Python 3 and
     making it work on Python 2.
:orphan:

==============================================
Connection Diagrams of The IPython ZMQ Cluster
==============================================

IPython parallel has moved to ipyparallel -
see :ref:`ipyparallel:/reference/connections.md` for the documentation.
Making simple Python wrapper kernels
====================================

.. versionadded:: 3.0

You can now re-use the kernel machinery in IPython to easily make new kernels.
This is useful for languages that have Python bindings, such as `Octave
<http://www.gnu.org/software/octave/>`_ (via
`Oct2Py <http://blink1073.github.io/oct2py/>`_), or languages
where the REPL can be controlled in a tty using `pexpect <https://pexpect.readthedocs.io/en/latest/>`_,
such as bash.

.. seealso::

   `bash_kernel <https://github.com/takluyver/bash_kernel>`_
     A simple kernel for bash, written using this machinery

Required steps
--------------

Subclass :class:`ipykernel.kernelbase.Kernel`, and implement the
following methods and attributes:

.. class:: MyKernel

   .. attribute:: implementation
                  implementation_version
                  language
                  language_version
                  banner
    
     Information for :ref:`msging_kernel_info` replies. 'Implementation' refers
     to the kernel (e.g. IPython), and 'language' refers to the language it
     interprets (e.g. Python). The 'banner' is displayed to the user in console
     UIs before the first prompt. All of these values are strings.

   .. attribute:: language_info

     Language information for :ref:`msging_kernel_info` replies, in a dictionary.
     This should contain the key ``mimetype`` with the mimetype of code in the
     target language (e.g. ``'text/x-python'``), and ``file_extension`` (e.g.
     ``'py'``).
     It may also contain keys ``codemirror_mode`` and ``pygments_lexer`` if they
     need to differ from :attr:`language`.

     Other keys may be added to this later.

   .. method:: do_execute(code, silent, store_history=True, user_expressions=None, allow_stdin=False)
   
     Execute user code.
     
     :param str code: The code to be executed.
     :param bool silent: Whether to display output.
     :param bool store_history: Whether to record this code in history and
         increase the execution count. If silent is True, this is implicitly
         False.
     :param dict user_expressions: Mapping of names to expressions to evaluate
         after the code has run. You can ignore this if you need to.
     :param bool allow_stdin: Whether the frontend can provide input on request
         (e.g. for Python's :func:`raw_input`).
     
     Your method should return a dict containing the fields described in
     :ref:`execution_results`. To display output, it can send messages
     using :meth:`~ipykernel.kernelbase.Kernel.send_response`.
     See :doc:`messaging` for details of the different message types.

To launch your kernel, add this at the end of your module::

    if __name__ == '__main__':
        from ipykernel.kernelapp import IPKernelApp
        IPKernelApp.launch_instance(kernel_class=MyKernel)

Example
-------

``echokernel.py`` will simply echo any input it's given to stdout::

    from ipykernel.kernelbase import Kernel

    class EchoKernel(Kernel):
        implementation = 'Echo'
        implementation_version = '1.0'
        language = 'no-op'
        language_version = '0.1'
        language_info = {'mimetype': 'text/plain'}
        banner = "Echo kernel - as useful as a parrot"

        def do_execute(self, code, silent, store_history=True, user_expressions=None,
                       allow_stdin=False):
            if not silent:
                stream_content = {'name': 'stdout', 'text': code}
                self.send_response(self.iopub_socket, 'stream', stream_content)

            return {'status': 'ok',
                    # The base class increments the execution count
                    'execution_count': self.execution_count,
                    'payload': [],
                    'user_expressions': {},
                   }

    if __name__ == '__main__':
        from ipykernel.kernelapp import IPKernelApp
        IPKernelApp.launch_instance(kernel_class=EchoKernel)

Here's the Kernel spec ``kernel.json`` file for this::

    {"argv":["python","-m","echokernel", "-f", "{connection_file}"],
     "display_name":"Echo"
    }


Optional steps
--------------

You can override a number of other methods to improve the functionality of your
kernel. All of these methods should return a dictionary as described in the
relevant section of the :doc:`messaging spec <messaging>`.

.. class:: MyBetterKernel

   .. method:: do_complete(code, cusor_pos)

     Code completion
     
     :param str code: The code already present
     :param int cursor_pos: The position in the code where completion is requested
     
     .. seealso::
     
        :ref:`msging_completion` messages

   .. method:: do_inspect(code, cusor_pos, detail_level=0)

     Object introspection
     
     :param str code: The code
     :param int cursor_pos: The position in the code where introspection is requested
     :param int detail_level: 0 or 1 for more or less detail. In IPython, 1 gets
         the source code.
     
     .. seealso::
     
        :ref:`msging_inspection` messages

   .. method:: do_history(hist_access_type, output, raw, session=None, start=None, stop=None, n=None, pattern=None, unique=False)

     History access. Only the relevant parameters for the type of history
     request concerned will be passed, so your method definition must have defaults
     for all the arguments shown with defaults here.

     .. seealso::
     
        :ref:`msging_history` messages

   .. method:: do_is_complete(code)
   
     Is code entered in a console-like interface complete and ready to execute,
     or should a continuation prompt be shown?
     
     :param str code: The code entered so far - possibly multiple lines
     
     .. seealso::
     
        :ref:`msging_is_complete` messages

   .. method:: do_shutdown(restart)

     Shutdown the kernel. You only need to handle your own clean up - the kernel
     machinery will take care of cleaning up its own things before stopping.
     
     :param bool restart: Whether the kernel will be started again afterwards
     
     .. seealso::
     
        :ref:`msging_shutdown` messages
:orphan:

==========================
Making kernels for IPython
==========================

Kernels are now part of Jupyter - see
:ref:`jupyterclient:kernels` for the documentation.
How IPython works
=================

Terminal IPython
----------------

When you type ``ipython``, you get the original IPython interface, running in
the terminal. It does something like this::

    while True:
        code = input(">>> ")
        exec(code)

Of course, it's much more complex, because it has to deal with multi-line
code, tab completion using :mod:`readline`, magic commands, and so on. But the
model is like that: prompt the user for some code, and when they've entered it,
exec it in the same process. This model is often called a REPL, or
Read-Eval-Print-Loop.

The IPython Kernel
------------------

All the other interfaces—the Notebook, the Qt console, ``ipython console`` in
the terminal, and third party interfaces—use the IPython Kernel. This is a
separate process which is responsible for running user code, and things like
computing possible completions. Frontends communicate with it using JSON
messages sent over `ZeroMQ <http://zeromq.org/>`_ sockets; the protocol they use is described in
:ref:`jupyterclient:messaging`.

The core execution machinery for the kernel is shared with terminal IPython:

.. image:: figs/ipy_kernel_and_terminal.png

A kernel process can be connected to more than one frontend simultaneously. In
this case, the different frontends will have access to the same variables.

.. TODO: Diagram illustrating this?

This design was intended to allow easy development of different frontends based
on the same kernel, but it also made it possible to support new languages in the
same frontends, by developing kernels in those languages, and we are refining
IPython to make that more practical.

Today, there are two ways to develop a kernel for another language. Wrapper
kernels reuse the communications machinery from IPython, and implement only the
core execution part. Native kernels implement execution and communications in
the target language:

.. image:: figs/other_kernels.png

Wrapper kernels are easier to write quickly for languages that have good Python
wrappers, like `octave_kernel <https://pypi.python.org/pypi/octave_kernel>`_, or
languages where it's impractical to implement the communications machinery, like
`bash_kernel <https://pypi.python.org/pypi/bash_kernel>`_. Native kernels are
likely to be better maintained by the community using them, like
`IJulia <https://github.com/JuliaLang/IJulia.jl>`_ or `IHaskell <https://github.com/gibiansky/IHaskell>`_.

.. seealso::

   :ref:`jupyterclient:kernels`
   
   :doc:`wrapperkernels`

:orphan:

================================
Messaging for Parallel Computing
================================

IPython parallel has moved to ipyparallel -
see :ref:`ipyparallel:/reference/messages.md` for the documentation.
=========================
IPython GUI Support Notes
=========================

IPython allows GUI event loops to be run in an interactive IPython session.
This is done using Python's PyOS_InputHook hook which Python calls
when the :func:`raw_input` function is called and waiting for user input.
IPython has versions of this hook for wx, pyqt4 and pygtk.

When a GUI program is used interactively within IPython, the event loop of
the GUI should *not* be started. This is because, the PyOS_Inputhook itself
is responsible for iterating the GUI event loop.

IPython has facilities for installing the needed input hook for each GUI 
toolkit and for creating the needed main GUI application object. Usually,
these main application objects should be created only once and for some
GUI toolkits, special options have to be passed to the application object
to enable it to function properly in IPython.

We need to answer the following questions:

* Who is responsible for creating the main GUI application object, IPython
  or third parties (matplotlib, enthought.traits, etc.)?

* What is the proper way for third party code to detect if a GUI application
  object has already been created?  If one has been created, how should
  the existing instance be retrieved?

* In a GUI application object has been created, how should third party code
  detect if the GUI event loop is running. It is not sufficient to call the
  relevant function methods in the GUI toolkits (like ``IsMainLoopRunning``)
  because those don't know if the GUI event loop is running through the
  input hook.

* We might need a way for third party code to determine if it is running
  in IPython or not.  Currently, the only way of running GUI code in IPython
  is by using the input hook, but eventually, GUI based versions of IPython
  will allow the GUI event loop in the more traditional manner. We will need
  a way for third party code to distinguish between these two cases.

Here is some sample code I have been using to debug this issue::

    from matplotlib import pyplot as plt

    from enthought.traits import api as traits

    class Foo(traits.HasTraits):
        a = traits.Float()

    f = Foo()
    f.configure_traits()

    plt.plot(range(10))


.. _developer_guide:

=====================================================
Developer's guide for third party tools and libraries
=====================================================

.. important::

    This guide contains information for developers of third party tools and
    libraries that use IPython. Alternatively, documentation for core
    **IPython** development can be found in the :doc:`../coredev/index`.

.. toctree::
   :maxdepth: 1

   how_ipython_works
   wrapperkernels
   execution
   lexer
   config
   inputhook_app
:orphan:

Messaging in IPython
====================

The message specification is now part of Jupyter - see
:ref:`jupyterclient:messaging` for the documentation.
.. _console_lexer:

New IPython Console Lexer
-------------------------

.. versionadded:: 2.0.0

The IPython console lexer has been rewritten and now supports tracebacks
and customized input/output prompts. An entire suite of lexers is now
available at :mod:`IPython.lib.lexers`. These include:

IPythonLexer & IPython3Lexer
  Lexers for pure IPython (python + magic/shell commands)

IPythonPartialTracebackLexer & IPythonTracebackLexer
  Supports 2.x and 3.x via the keyword `python3`. The partial traceback
  lexer reads everything but the Python code appearing in a traceback.
  The full lexer combines the partial lexer with an IPython lexer.

IPythonConsoleLexer
  A lexer for IPython console sessions, with support for tracebacks.
  Supports 2.x and 3.x via the keyword `python3`.

IPyLexer
  A friendly lexer which examines the first line of text and from it,
  decides whether to use an IPython lexer or an IPython console lexer.
  Supports 2.x and 3.x via the keyword `python3`.

Previously, the :class:`IPythonConsoleLexer` class was available at
:mod:`IPython.sphinxext.ipython_console_hightlight`.  It was inserted
into Pygments' list of available lexers under the name `ipython`.  It should
be mentioned that this name is inaccurate, since an IPython console session
is not the same as IPython code (which itself is a superset of the Python
language).

Now, the Sphinx extension inserts two console lexers into Pygments' list of
available lexers. Both are IPyLexer instances under the names: `ipython` and
`ipython3`. Although the names can be confusing (as mentioned above), their
continued use is, in part, to maintain backwards compatibility and to
aid typical usage. If a project needs to make Pygments aware of more than just
the IPyLexer class, then one should not make the IPyLexer class available under
the name `ipython` and use `ipy` or some other non-conflicting value.

Code blocks such as:

.. code-block:: rst

    .. code-block:: ipython

        In [1]: 2**2
        Out[1]: 4

will continue to work as before, but now, they will also properly highlight
tracebacks.  For pure IPython code, the same lexer will also work:

.. code-block:: rst

    .. code-block:: ipython

        x = ''.join(map(str, range(10)))
        !echo $x

Since the first line of the block did not begin with a standard IPython console
prompt, the entire block is assumed to consist of IPython code instead.
.. _history:

=======
History
=======

Origins
=======

IPython was starting in 2001 by Fernando Perez while he was a graduate student
at the University of Colorado, Boulder. IPython as we know it today grew out
of the following three projects:

* ipython by Fernando Pérez. Fernando began using Python and ipython began as
  an outgrowth of his desire for things like Mathematica-style prompts, access
  to previous output (again like Mathematica's % syntax) and a flexible
  configuration system (something better than :envvar:`PYTHONSTARTUP`).
* IPP by Janko Hauser. Very well organized, great usability. Had
  an old help system. IPP was used as the "container" code into
  which Fernando added the functionality from ipython and LazyPython.
* LazyPython by Nathan Gray. Simple but very powerful. The quick
  syntax (auto parens, auto quotes) and verbose/colored tracebacks
  were all taken from here.

Here is how Fernando describes the early history of IPython:

    When I found out about IPP and LazyPython I tried to join all three
    into a unified system. I thought this could provide a very nice
    working environment, both for regular programming and scientific
    computing: shell-like features, IDL/Matlab numerics, Mathematica-type
    prompt history and great object introspection and help facilities. I
    think it worked reasonably well, though it was a lot more work than I
    had initially planned.

.. _about_index:

=============
About IPython
=============

.. toctree::
   :maxdepth: 1

   history
   license_and_copyright

.. _tips:

=====================
IPython Tips & Tricks
=====================

The `IPython cookbook
<https://github.com/ipython/ipython/wiki?path=Cookbook>`_ details more things
you can do with IPython.

.. This is not in the current version:


Embed IPython in your programs
------------------------------

A few lines of code are enough to load a complete IPython inside your own
programs, giving you the ability to work with your data interactively after
automatic processing has been completed. See :ref:`the embedding section <embedding>`.

Run doctests
------------

Run your doctests from within IPython for development and debugging. The
special ``%doctest_mode`` command toggles a mode where the prompt, output and
exceptions display matches as closely as possible that of the default Python
interpreter. In addition, this mode allows you to directly paste in code that
contains leading '>>>' prompts, even if they have extra leading whitespace
(as is common in doctest files). This combined with the ``%hist -t`` call to
see your translated history allows for an easy doctest workflow, where you
can go from doctest to interactive execution to pasting into valid Python code
as needed.

Use IPython to present interactive demos
----------------------------------------

Use the :class:`IPython.lib.demo.Demo` class to load any Python script as an interactive
demo. With a minimal amount of simple markup, you can control the execution of
the script, stopping as needed. See :ref:`here <interactive_demos>` for more.

Suppress output
---------------

Put a ';' at the end of a line to suppress the printing of output. This is
useful when doing calculations which generate long output you are not
interested in seeing. It also keeps the object out of the output cache, so if
you're working with large temporary objects, they'll be released from memory sooner.

Lightweight 'version control'
-----------------------------

When you call ``%edit`` with no arguments, IPython opens an empty editor
with a temporary file, and it returns the contents of your editing
session as a string variable. Thanks to IPython's output caching
mechanism, this is automatically stored::

    In [1]: %edit

    IPython will make a temporary file named: /tmp/ipython_edit_yR-HCN.py

    Editing... done. Executing edited code...

    hello - this is a temporary file

    Out[1]: "print('hello - this is a temporary file')\n"

Now, if you call ``%edit -p``, IPython tries to open an editor with the
same data as the last time you used %edit. So if you haven't used %edit
in the meantime, this same contents will reopen; however, it will be
done in a new file. This means that if you make changes and you later
want to find an old version, you can always retrieve it by using its
output number, via '%edit _NN', where NN is the number of the output
prompt.

Continuing with the example above, this should illustrate this idea::

    In [2]: edit -p

    IPython will make a temporary file named: /tmp/ipython_edit_nA09Qk.py

    Editing... done. Executing edited code...

    hello - now I made some changes

    Out[2]: "print('hello - now I made some changes')\n"

    In [3]: edit _1

    IPython will make a temporary file named: /tmp/ipython_edit_gy6-zD.py

    Editing... done. Executing edited code...

    hello - this is a temporary file

    IPython version control at work :)

    Out[3]: "print('hello - this is a temporary file')\nprint('IPython version control at work :)')\n"


This section was written after a contribution by Alexander Belchenko on
the IPython user list.

=================
IPython reference
=================

.. _command_line_options:

Command-line usage
==================

You start IPython with the command::

    $ ipython [options] files

If invoked with no options, it executes the file and exits, passing the
remaining arguments to the script, just as if you had specified the same
command with python. You may need to specify `--` before args to be passed
to the script, to prevent IPython from attempting to parse them.
If you add the ``-i`` flag, it drops you into the interpreter while still
acknowledging any options you may have set in your ``ipython_config.py``. This
behavior is different from standard Python, which when called as python ``-i``
will only execute one file and ignore your configuration setup.

Please note that some of the configuration options are not available at the
command line, simply because they are not practical here. Look into your
configuration files for details on those. There are separate configuration files
for each profile, and the files look like :file:`ipython_config.py` or
:file:`ipython_config_{frontendname}.py`.  Profile directories look like
:file:`profile_{profilename}` and are typically installed in the
:envvar:`IPYTHONDIR` directory, which defaults to :file:`$HOME/.ipython`. For
Windows users, :envvar:`HOME` resolves to :file:`C:\\Users\\{YourUserName}` in
most instances.

Command-line Options
--------------------

To see the options IPython accepts, use ``ipython --help`` (and you probably
should run the output through a pager such as ``ipython --help | less`` for
more convenient reading).  This shows all the options that have a single-word
alias to control them, but IPython lets you configure all of its objects from
the command-line by passing the full class name and a corresponding value; type
``ipython --help-all`` to see this full list.  For example::

    $ ipython --help-all
    <...snip...>
    --matplotlib=<CaselessStrEnum> (InteractiveShellApp.matplotlib)
        Default: None
        Choices: ['auto', 'gtk', 'gtk3', 'gtk4', 'inline', 'nbagg', 'notebook', 'osx', 'qt', 'qt4', 'qt5', 'tk', 'wx']
        Configure matplotlib for interactive use with the default matplotlib
        backend.
    <...snip...>


Indicate that the following::

   $ ipython --matplotlib qt


is equivalent to::

   $ ipython --InteractiveShellApp.matplotlib='qt'

Note that in the second form, you *must* use the equal sign, as the expression
is evaluated as an actual Python assignment.  While in the above example the
short form is more convenient, only the most common options have a short form,
while any configurable variable in IPython can be set at the command-line by
using the long form.  This long form is the same syntax used in the
configuration files, if you want to set these options permanently.


Interactive use
===============

IPython is meant to work as a drop-in replacement for the standard interactive
interpreter. As such, any code which is valid python should execute normally
under IPython (cases where this is not true should be reported as bugs). It
does, however, offer many features which are not available at a standard python
prompt. What follows is a list of these.


Caution for Windows users
-------------------------

Windows, unfortunately, uses the ``\`` character as a path separator. This is a
terrible choice, because ``\`` also represents the escape character in most
modern programming languages, including Python. For this reason, using '/'
character is recommended if you have problems with ``\``.  However, in Windows
commands '/' flags options, so you can not use it for the root directory. This
means that paths beginning at the root must be typed in a contrived manner
like: ``%copy \opt/foo/bar.txt \tmp``

.. _magic:

Magic command system
--------------------

IPython will treat any line whose first character is a % as a special
call to a 'magic' function. These allow you to control the behavior of
IPython itself, plus a lot of system-type features. They are all
prefixed with a % character, but parameters are given without
parentheses or quotes.

Lines that begin with ``%%`` signal a *cell magic*: they take as arguments not
only the rest of the current line, but all lines below them as well, in the
current execution block.  Cell magics can in fact make arbitrary modifications
to the input they receive, which need not even be valid Python code at all.
They receive the whole block as a single string.

As a line magic example, the :magic:`cd` magic works just like the OS command of
the same name::

      In [8]: %cd
      /home/fperez

The following uses the builtin :magic:`timeit` in cell mode::

  In [10]: %%timeit x = range(10000)
      ...: min(x)
      ...: max(x)
      ...:
  1000 loops, best of 3: 438 us per loop

In this case, ``x = range(10000)`` is called as the line argument, and the
block with ``min(x)`` and ``max(x)`` is called as the cell body.  The
:magic:`timeit` magic receives both.

If you have 'automagic' enabled (as it is by default), you don't need to type in
the single ``%`` explicitly for line magics; IPython will scan its internal
list of magic functions and call one if it exists. With automagic on you can
then just type ``cd mydir`` to go to directory 'mydir'::

      In [9]: cd mydir
      /home/fperez/mydir

Cell magics *always* require an explicit ``%%`` prefix, automagic
calling only works for line magics.

The automagic system has the lowest possible precedence in name searches, so
you can freely use variables with the same names as magic commands. If a magic
command is 'shadowed' by a variable, you will need the explicit ``%`` prefix to
use it:

.. sourcecode:: ipython

    In [1]: cd ipython     # %cd is called by automagic
    /home/fperez/ipython

    In [2]: cd=1 	   # now cd is just a variable

    In [3]: cd .. 	   # and doesn't work as a function anymore
    File "<ipython-input-3-9fedb3aff56c>", line 1
      cd ..
          ^
    SyntaxError: invalid syntax


    In [4]: %cd .. 	   # but %cd always works
    /home/fperez

    In [5]: del cd     # if you remove the cd variable, automagic works again

    In [6]: cd ipython

    /home/fperez/ipython

Line magics, if they return a value, can be assigned to a variable using the
syntax ``l = %sx ls`` (which in this particular case returns the result of `ls`
as a python list). See :ref:`below <manual_capture>` for more information.

Type ``%magic`` for more information, including a list of all available magic
functions at any time and their docstrings. You can also type
``%magic_function_name?`` (see :ref:`below <dynamic_object_info>` for
information on the '?' system) to get information about any particular magic
function you are interested in.

The API documentation for the :mod:`IPython.core.magic` module contains the full
docstrings of all currently available magic commands.

.. seealso::

   :doc:`magics`
     A list of the line and cell magics available in IPython by default

   :ref:`defining_magics`
     How to define and register additional magic functions


Access to the standard Python help
----------------------------------

Simply type ``help()`` to access Python's standard help system. You can
also type ``help(object)`` for information about a given object, or
``help('keyword')`` for information on a keyword. You may need to configure your
PYTHONDOCS environment variable for this feature to work correctly.

.. _dynamic_object_info:

Dynamic object information
--------------------------

Typing ``?word`` or ``word?`` prints detailed information about an object. If
certain strings in the object are too long (e.g. function signatures) they get
snipped in the center for brevity. This system gives access variable types and
values, docstrings, function prototypes and other useful information.

If the information will not fit in the terminal, it is displayed in a pager
(``less`` if available, otherwise a basic internal pager).

Typing ``??word`` or ``word??`` gives access to the full information, including
the source code where possible. Long strings are not snipped.

The following magic functions are particularly useful for gathering
information about your working environment:

    * :magic:`pdoc` **<object>**: Print (or run through a pager if too long) the
      docstring for an object. If the given object is a class, it will
      print both the class and the constructor docstrings.
    * :magic:`pdef` **<object>**: Print the call signature for any callable
      object. If the object is a class, print the constructor information.
    * :magic:`psource` **<object>**: Print (or run through a pager if too long)
      the source code for an object.
    * :magic:`pfile` **<object>**: Show the entire source file where an object was
      defined via a pager, opening it at the line where the object
      definition begins.
    * :magic:`who`/:magic:`whos`: These functions give information about identifiers
      you have defined interactively (not things you loaded or defined
      in your configuration files). %who just prints a list of
      identifiers and %whos prints a table with some basic details about
      each identifier.

The dynamic object information functions (?/??, ``%pdoc``,
``%pfile``, ``%pdef``, ``%psource``) work on object attributes, as well as
directly on variables. For example, after doing ``import os``, you can use
``os.path.abspath??``.


Command line completion
+++++++++++++++++++++++

At any time, hitting TAB will complete any available python commands or
variable names, and show you a list of the possible completions if
there's no unambiguous one. It will also complete filenames in the
current directory if no python names match what you've typed so far.


Search command history
++++++++++++++++++++++

IPython provides two ways for searching through previous input and thus
reduce the need for repetitive typing:

   1. Start typing, and then use the up and down arrow keys (or :kbd:`Ctrl-p`
      and :kbd:`Ctrl-n`) to search through only the history items that match
      what you've typed so far.
   2. Hit :kbd:`Ctrl-r`: to open a search prompt. Begin typing and the system
      searches your history for lines that contain what you've typed so
      far, completing as much as it can.

IPython will save your input history when it leaves and reload it next
time you restart it. By default, the history file is named
:file:`.ipython/profile_{name}/history.sqlite`.

Autoindent
++++++++++

Starting with 5.0, IPython uses `prompt_toolkit` in place of ``readline``,
it thus can recognize lines ending in ':' and indent the next line,
while also un-indenting automatically after 'raise' or 'return',
and support real multi-line editing as well as syntactic coloration
during edition.

This feature does not use the ``readline`` library anymore, so it will
not honor your :file:`~/.inputrc` configuration (or whatever
file your :envvar:`INPUTRC` environment variable points to).

In particular if you want to change the input mode to ``vi``, you will need to
set the ``TerminalInteractiveShell.editing_mode`` configuration  option of IPython.

Session logging and restoring
-----------------------------

You can log all input from a session either by starting IPython with the
command line switch ``--logfile=foo.py`` (see :ref:`here <command_line_options>`)
or by activating the logging at any moment with the magic function :magic:`logstart`.

Log files can later be reloaded by running them as scripts and IPython
will attempt to 'replay' the log by executing all the lines in it, thus
restoring the state of a previous session. This feature is not quite
perfect, but can still be useful in many cases.

The log files can also be used as a way to have a permanent record of
any code you wrote while experimenting. Log files are regular text files
which you can later open in your favorite text editor to extract code or
to 'clean them up' before using them to replay a session.

The :magic:`logstart` function for activating logging in mid-session is used as
follows::

    %logstart [log_name [log_mode]]

If no name is given, it defaults to a file named 'ipython_log.py' in your
current working directory, in 'rotate' mode (see below).

'%logstart name' saves to file 'name' in 'backup' mode. It saves your
history up to that point and then continues logging.

%logstart takes a second optional parameter: logging mode. This can be
one of (note that the modes are given unquoted):

    * [over:] overwrite existing log_name.
    * [backup:] rename (if exists) to log_name~ and start log_name.
    * [append:] well, that says it.
    * [rotate:] create rotating logs log_name.1~, log_name.2~, etc.

Adding the '-o' flag to '%logstart' magic (as in '%logstart -o [log_name [log_mode]]')
will also include output from iPython in the log file.

The :magic:`logoff` and :magic:`logon` functions allow you to temporarily stop and
resume logging to a file which had previously been started with
%logstart. They will fail (with an explanation) if you try to use them
before logging has been started.

.. _system_shell_access:

System shell access
-------------------

Any input line beginning with a ``!`` character is passed verbatim (minus
the ``!``, of course) to the underlying operating system. For example,
typing ``!ls`` will run 'ls' in the current directory.

.. _manual_capture:

Manual capture of command output and magic output
-------------------------------------------------

You can assign the result of a system command to a Python variable with the
syntax ``myfiles = !ls``. Similarly, the result of a magic (as long as it returns
a value) can be assigned to a variable.  For example, the syntax ``myfiles = %sx ls``
is equivalent to the above system command example (the :magic:`sx` magic runs a shell command
and captures the output).  Each of these gets machine
readable output from stdout (e.g. without colours), and splits on newlines. To
explicitly get this sort of output without assigning to a variable, use two
exclamation marks (``!!ls``) or the :magic:`sx` magic command without an assignment.
(However, ``!!`` commands cannot be assigned to a variable.)

The captured list in this example has some convenience features. ``myfiles.n`` or ``myfiles.s``
returns a string delimited by newlines or spaces, respectively. ``myfiles.p``
produces `path objects <http://pypi.python.org/pypi/path.py>`_ from the list items.
See :ref:`string_lists` for details.

IPython also allows you to expand the value of python variables when
making system calls. Wrap variables or expressions in {braces}::

    In [1]: pyvar = 'Hello world'
    In [2]: !echo "A python variable: {pyvar}"
    A python variable: Hello world
    In [3]: import math
    In [4]: x = 8
    In [5]: !echo {math.factorial(x)}
    40320

For simple cases, you can alternatively prepend $ to a variable name::

    In [6]: !echo $sys.argv
    [/home/fperez/usr/bin/ipython]
    In [7]: !echo "A system variable: $$HOME"  # Use $$ for literal $
    A system variable: /home/fperez

Note that `$$` is used to represent a literal `$`.

System command aliases
----------------------

The :magic:`alias` magic function allows you to define magic functions which are in fact
system shell commands. These aliases can have parameters.

``%alias alias_name cmd`` defines 'alias_name' as an alias for 'cmd'

Then, typing ``alias_name params`` will execute the system command 'cmd
params' (from your underlying operating system).

You can also define aliases with parameters using ``%s`` specifiers (one per
parameter). The following example defines the parts function as an
alias to the command ``echo first %s second %s`` where each ``%s`` will be
replaced by a positional parameter to the call to %parts::

    In [1]: %alias parts echo first %s second %s
    In [2]: parts A B
    first A second B
    In [3]: parts A
    ERROR: Alias <parts> requires 2 arguments, 1 given.

If called with no parameters, :magic:`alias` prints the table of currently
defined aliases.

The :magic:`rehashx` magic allows you to load your entire $PATH as
ipython aliases. See its docstring for further details.


.. _dreload:

Recursive reload
----------------

The :mod:`IPython.lib.deepreload` module allows you to recursively reload a
module: changes made to any of its dependencies will be reloaded without
having to exit. To start using it, do::

    from IPython.lib.deepreload import reload as dreload


Verbose and colored exception traceback printouts
-------------------------------------------------

IPython provides the option to see very detailed exception tracebacks,
which can be especially useful when debugging large programs. You can
run any Python file with the %run function to benefit from these
detailed tracebacks. Furthermore, both normal and verbose tracebacks can
be colored (if your terminal supports it) which makes them much easier
to parse visually.

See the magic :magic:`xmode` and :magic:`colors` functions for details.

These features are basically a terminal version of Ka-Ping Yee's cgitb
module, now part of the standard Python library.


.. _input_caching:

Input caching system
--------------------

IPython offers numbered prompts (In/Out) with input and output caching
(also referred to as 'input history'). All input is saved and can be
retrieved as variables (besides the usual arrow key recall), in
addition to the :magic:`rep` magic command that brings a history entry
up for editing on the next command line.

The following variables always exist:

* ``_i``, ``_ii``, ``_iii``: store previous, next previous and next-next
  previous inputs.

* ``In``, ``_ih`` : a list of all inputs; ``_ih[n]`` is the input from line
  ``n``. If you overwrite In with a variable of your own, you can remake the
  assignment to the internal list with a simple ``In=_ih``.

Additionally, global variables named ``_i<n>`` are dynamically created (``<n>``
being the prompt counter), so ``_i<n> == _ih[<n>] == In[<n>]``.

For example, what you typed at prompt 14 is available as ``_i14``, ``_ih[14]``
and ``In[14]``.

This allows you to easily cut and paste multi line interactive prompts
by printing them out: they print like a clean string, without prompt
characters. You can also manipulate them like regular variables (they
are strings), modify or exec them.

You can also re-execute multiple lines of input easily by using the magic
:magic:`rerun` or :magic:`macro` functions. The macro system also allows you to
re-execute previous lines which include magic function calls (which require
special processing). Type %macro? for more details on the macro system.

A history function :magic:`history` allows you to see any part of your input
history by printing a range of the _i variables.

You can also search ('grep') through your history by typing
``%hist -g somestring``. This is handy for searching for URLs, IP addresses,
etc. You can bring history entries listed by '%hist -g' up for editing
with the %recall command, or run them immediately with :magic:`rerun`.

.. _output_caching:

Output caching system
---------------------

For output that is returned from actions, a system similar to the input
cache exists but using _ instead of _i. Only actions that produce a
result (NOT assignments, for example) are cached. If you are familiar
with Mathematica, IPython's _ variables behave exactly like
Mathematica's % variables.

The following variables always exist:

    * [_] (a single underscore): stores previous output, like Python's
      default interpreter.
    * [__] (two underscores): next previous.
    * [___] (three underscores): next-next previous.

Additionally, global variables named _<n> are dynamically created (<n>
being the prompt counter), such that the result of output <n> is always
available as _<n> (don't use the angle brackets, just the number, e.g.
``_21``).

These variables are also stored in a global dictionary (not a
list, since it only has entries for lines which returned a result)
available under the names _oh and Out (similar to _ih and In). So the
output from line 12 can be obtained as ``_12``, ``Out[12]`` or ``_oh[12]``. If you
accidentally overwrite the Out variable you can recover it by typing
``Out=_oh`` at the prompt.

This system obviously can potentially put heavy memory demands on your
system, since it prevents Python's garbage collector from removing any
previously computed results. You can control how many results are kept
in memory with the configuration option ``InteractiveShell.cache_size``.
If you set it to 0, output caching is disabled. You can also use the :magic:`reset`
and :magic:`xdel` magics to clear large items from memory.

Directory history
-----------------

Your history of visited directories is kept in the global list _dh, and
the magic :magic:`cd` command can be used to go to any entry in that list. The
:magic:`dhist` command allows you to view this history. Do ``cd -<TAB>`` to
conveniently view the directory history.


Automatic parentheses and quotes
--------------------------------

These features were adapted from Nathan Gray's LazyPython. They are
meant to allow less typing for common situations.

Callable objects (i.e. functions, methods, etc) can be invoked like this
(notice the commas between the arguments)::

    In [1]: callable_ob arg1, arg2, arg3
    ------> callable_ob(arg1, arg2, arg3)

.. note::
   This feature is disabled by default. To enable it, use the ``%autocall``
   magic command. The commands below with special prefixes will always work,
   however.

You can force automatic parentheses by using '/' as the first character
of a line. For example::

    In [2]: /globals # becomes 'globals()'

Note that the '/' MUST be the first character on the line! This won't work::

    In [3]: print /globals # syntax error

In most cases the automatic algorithm should work, so you should rarely
need to explicitly invoke /. One notable exception is if you are trying
to call a function with a list of tuples as arguments (the parenthesis
will confuse IPython)::

    In [4]: zip (1,2,3),(4,5,6) # won't work

but this will work::

    In [5]: /zip (1,2,3),(4,5,6)
    ------> zip ((1,2,3),(4,5,6))
    Out[5]: [(1, 4), (2, 5), (3, 6)]

IPython tells you that it has altered your command line by displaying
the new command line preceded by ``--->``.

You can force automatic quoting of a function's arguments by using ``,``
or ``;`` as the first character of a line. For example::

    In [1]: ,my_function /home/me  # becomes my_function("/home/me")

If you use ';' the whole argument is quoted as a single string, while ',' splits
on whitespace::

    In [2]: ,my_function a b c    # becomes my_function("a","b","c")

    In [3]: ;my_function a b c    # becomes my_function("a b c")

Note that the ',' or ';' MUST be the first character on the line! This
won't work::

    In [4]: x = ,my_function /home/me # syntax error

IPython as your default Python environment
==========================================

Python honors the environment variable :envvar:`PYTHONSTARTUP` and will
execute at startup the file referenced by this variable. If you put the
following code at the end of that file, then IPython will be your working
environment anytime you start Python::

    import os, IPython
    os.environ['PYTHONSTARTUP'] = ''  # Prevent running this again
    IPython.start_ipython()
    raise SystemExit

The ``raise SystemExit`` is needed to exit Python when
it finishes, otherwise you'll be back at the normal Python ``>>>``
prompt.

This is probably useful to developers who manage multiple Python
versions and don't want to have correspondingly multiple IPython
versions. Note that in this mode, there is no way to pass IPython any
command-line options, as those are trapped first by Python itself.

.. _Embedding:

Embedding IPython
=================

You can start a regular IPython session with

.. sourcecode:: python

    import IPython
    IPython.start_ipython(argv=[])

at any point in your program.  This will load IPython configuration,
startup files, and everything, just as if it were a normal IPython session.
For information on setting configuration options when running IPython from
python, see :ref:`configure_start_ipython`.

It is also possible to embed an IPython shell in a namespace in your Python
code. This allows you to evaluate dynamically the state of your code, operate
with your variables, analyze them, etc. For example, if you run the following
code snippet::

  import IPython

  a = 42
  IPython.embed()

and within the IPython shell, you reassign `a` to `23` to do further testing of 
some sort, you can then exit::

  >>> IPython.embed()
  Python 3.6.2 (default, Jul 17 2017, 16:44:45) 
  Type 'copyright', 'credits' or 'license' for more information
  IPython 6.2.0.dev -- An enhanced Interactive Python. Type '?' for help.

  In [1]: a = 23

  In [2]: exit()

Once you exit and print `a`, the value 23 will be shown::


  In: print(a)
  23

It's important to note that the code run in the embedded IPython shell will 
*not* change the state of your code and variables, **unless** the shell is 
contained within the global namespace. In the above example, `a` is changed 
because this is true.

To further exemplify this, consider the following example::

  import IPython
  def do():
      a = 42
      print(a)
      IPython.embed()
      print(a)

Now if call the function and complete the state changes as we did above, the
value `42` will be printed. Again, this is because it's not in the global
namespace:: 

  do()

Running a file with the above code can lead to the following session::

  >>> do()
  42
  Python 3.6.2 (default, Jul 17 2017, 16:44:45) 
  Type 'copyright', 'credits' or 'license' for more information
  IPython 6.2.0.dev -- An enhanced Interactive Python. Type '?' for help.

  In [1]: a = 23

  In [2]: exit()
  42

.. note::

  At present, embedding IPython cannot be done from inside IPython.
  Run the code samples below outside IPython.

This feature allows you to easily have a fully functional python
environment for doing object introspection anywhere in your code with a
simple function call. In some cases a simple print statement is enough,
but if you need to do more detailed analysis of a code fragment this
feature can be very valuable.

It can also be useful in scientific computing situations where it is
common to need to do some automatic, computationally intensive part and
then stop to look at data, plots, etc.
Opening an IPython instance will give you full access to your data and
functions, and you can resume program execution once you are done with
the interactive part (perhaps to stop again later, as many times as
needed).

The following code snippet is the bare minimum you need to include in
your Python programs for this to work (detailed examples follow later)::

    from IPython import embed

    embed() # this call anywhere in your program will start IPython

You can also embed an IPython *kernel*, for use with qtconsole, etc. via
``IPython.embed_kernel()``. This should work the same way, but you can
connect an external frontend (``ipython qtconsole`` or ``ipython console``),
rather than interacting with it in the terminal.

You can run embedded instances even in code which is itself being run at
the IPython interactive prompt with '%run <filename>'. Since it's easy
to get lost as to where you are (in your top-level IPython or in your
embedded one), it's a good idea in such cases to set the in/out prompts
to something different for the embedded instances. The code examples
below illustrate this.

You can also have multiple IPython instances in your program and open
them separately, for example with different options for data
presentation. If you close and open the same instance multiple times,
its prompt counters simply continue from each execution to the next.

Please look at the docstrings in the :mod:`~IPython.frontend.terminal.embed`
module for more details on the use of this system.

The following sample file illustrating how to use the embedding
functionality is provided in the examples directory as embed_class_long.py.
It should be fairly self-explanatory:

.. literalinclude:: ../../../examples/Embedding/embed_class_long.py
    :language: python

Once you understand how the system functions, you can use the following
code fragments in your programs which are ready for cut and paste:

.. literalinclude:: ../../../examples/Embedding/embed_class_short.py
    :language: python

Using the Python debugger (pdb)
===============================

Running entire programs via pdb
-------------------------------

pdb, the Python debugger, is a powerful interactive debugger which
allows you to step through code, set breakpoints, watch variables,
etc.  IPython makes it very easy to start any script under the control
of pdb, regardless of whether you have wrapped it into a 'main()'
function or not. For this, simply type ``%run -d myscript`` at an
IPython prompt. See the :magic:`run` command's documentation for more details, including
how to control where pdb will stop execution first.

For more information on the use of the pdb debugger, see :ref:`debugger-commands`
in the Python documentation.

IPython extends the debugger with a few useful additions, like coloring of
tracebacks. The debugger will adopt the color scheme selected for IPython.

The ``where`` command has also been extended to take as argument the number of
context line to show. This allows to a many line of context on shallow stack trace:

.. code::

    In [5]: def foo(x):
    ...:     1
    ...:     2
    ...:     3
    ...:     return 1/x+foo(x-1)
    ...:     5
    ...:     6
    ...:     7
    ...:

    In[6]: foo(1)
    # ...
    ipdb> where 8
    <ipython-input-6-9e45007b2b59>(1)<module>
    ----> 1 foo(1)

    <ipython-input-5-7baadc3d1465>(5)foo()
        1 def foo(x):
        2     1
        3     2
        4     3
    ----> 5     return 1/x+foo(x-1)
        6     5
        7     6
        8     7

    > <ipython-input-5-7baadc3d1465>(5)foo()
        1 def foo(x):
        2     1
        3     2
        4     3
    ----> 5     return 1/x+foo(x-1)
        6     5
        7     6
        8     7


And less context on shallower Stack Trace:

.. code::

    ipdb> where 1
    <ipython-input-13-afa180a57233>(1)<module>
    ----> 1 foo(7)

    <ipython-input-5-7baadc3d1465>(5)foo()
    ----> 5     return 1/x+foo(x-1)

    <ipython-input-5-7baadc3d1465>(5)foo()
    ----> 5     return 1/x+foo(x-1)

    <ipython-input-5-7baadc3d1465>(5)foo()
    ----> 5     return 1/x+foo(x-1)

    <ipython-input-5-7baadc3d1465>(5)foo()
    ----> 5     return 1/x+foo(x-1)


Post-mortem debugging
---------------------

Going into a debugger when an exception occurs can be
extremely useful in order to find the origin of subtle bugs, because pdb
opens up at the point in your code which triggered the exception, and
while your program is at this point 'dead', all the data is still
available and you can walk up and down the stack frame and understand
the origin of the problem.

You can use the :magic:`debug` magic after an exception has occurred to start
post-mortem debugging. IPython can also call debugger every time your code
triggers an uncaught exception. This feature can be toggled with the :magic:`pdb` magic
command, or you can start IPython with the ``--pdb`` option.

For a post-mortem debugger in your programs outside IPython,
put the following lines toward the top of your 'main' routine::

    import sys
    from IPython.core import ultratb
    sys.excepthook = ultratb.FormattedTB(mode='Verbose',
    color_scheme='Linux', call_pdb=1)

The mode keyword can be either 'Verbose' or 'Plain', giving either very
detailed or normal tracebacks respectively. The color_scheme keyword can
be one of 'NoColor', 'Linux' (default) or 'LightBG'. These are the same
options which can be set in IPython with ``--colors`` and ``--xmode``.

This will give any of your programs detailed, colored tracebacks with
automatic invocation of pdb.

.. _pasting_with_prompts:

Pasting of code starting with Python or IPython prompts
=======================================================

IPython is smart enough to filter out input prompts, be they plain Python ones
(``>>>`` and ``...``) or IPython ones (``In [N]:`` and ``...:``).  You can
therefore copy and paste from existing interactive sessions without worry.

The following is a 'screenshot' of how things work, copying an example from the
standard Python tutorial::

    In [1]: >>> # Fibonacci series:

    In [2]: ... # the sum of two elements defines the next

    In [3]: ... a, b = 0, 1

    In [4]: >>> while b < 10:
       ...:     ...     print(b)
       ...:     ...     a, b = b, a+b
       ...:
    1
    1
    2
    3
    5
    8

And pasting from IPython sessions works equally well::

    In [1]: In [5]: def f(x):
       ...:        ...:     "A simple function"
       ...:        ...:     return x**2
       ...:    ...:

    In [2]: f(3)
    Out[2]: 9

.. _gui_support:

GUI event loop support
======================

IPython has excellent support for working interactively with Graphical User
Interface (GUI) toolkits, such as wxPython, PyQt4/PySide, PyGTK and Tk. This is
implemented by running the toolkit's event loop while IPython is waiting for
input.

For users, enabling GUI event loop integration is simple.  You simple use the
:magic:`gui` magic as follows::

    %gui [GUINAME]

With no arguments, ``%gui`` removes all GUI support.  Valid ``GUINAME``
arguments include ``wx``, ``qt``, ``qt5``, ``gtk``, ``gtk3`` ``gtk4``, and
``tk``.

Thus, to use wxPython interactively and create a running :class:`wx.App`
object, do::

    %gui wx

You can also start IPython with an event loop set up using the `--gui`
flag::

    $ ipython --gui=qt

For information on IPython's matplotlib_ integration (and the ``matplotlib``
mode) see :ref:`this section <matplotlib_support>`.

For developers that want to integrate additional event loops with IPython, see
:doc:`/config/eventloops`.

When running inside IPython with an integrated event loop, a GUI application
should *not* start its own event loop. This means that applications that are
meant to be used both
in IPython and as standalone apps need to have special code to detects how the
application is being run. We highly recommend using IPython's support for this.
Since the details vary slightly between toolkits, we point you to the various
examples in our source directory :file:`examples/IPython Kernel/gui/` that
demonstrate these capabilities.

PyQt and PySide
---------------

.. attempt at explanation of the complete mess that is Qt support

When you use ``--gui=qt`` or ``--matplotlib=qt``, IPython can work with either
PyQt4 or PySide.  There are three options for configuration here, because
PyQt4 has two APIs for QString and QVariant: v1, which is the default on
Python 2, and the more natural v2, which is the only API supported by PySide.
v2 is also the default for PyQt4 on Python 3.  IPython's code for the QtConsole
uses v2, but you can still use any interface in your code, since the
Qt frontend is in a different process.

The default will be to import PyQt4 without configuration of the APIs, thus
matching what most applications would expect. It will fall back to PySide if
PyQt4 is unavailable.

If specified, IPython will respect the environment variable ``QT_API`` used
by ETS.  ETS 4.0 also works with both PyQt4 and PySide, but it requires
PyQt4 to use its v2 API.  So if ``QT_API=pyside`` PySide will be used,
and if ``QT_API=pyqt`` then PyQt4 will be used *with the v2 API* for
QString and QVariant, so ETS codes like MayaVi will also work with IPython.

If you launch IPython in matplotlib mode with ``ipython --matplotlib=qt``,
then IPython will ask matplotlib which Qt library to use (only if QT_API is
*not set*), via the 'backend.qt4' rcParam.  If matplotlib is version 1.0.1 or
older, then IPython will always use PyQt4 without setting the v2 APIs, since
neither v2 PyQt nor PySide work.

.. warning::

    Note that this means for ETS 4 to work with PyQt4, ``QT_API`` *must* be set
    to work with IPython's qt integration, because otherwise PyQt4 will be
    loaded in an incompatible mode.

    It also means that you must *not* have ``QT_API`` set if you want to
    use ``--gui=qt`` with code that requires PyQt4 API v1.


.. _matplotlib_support:

Plotting with matplotlib
========================

matplotlib_ provides high quality 2D and 3D plotting for Python. matplotlib_
can produce plots on screen using a variety of GUI toolkits, including Tk,
PyGTK, PyQt4 and wxPython. It also provides a number of commands useful for
scientific computing, all with a syntax compatible with that of the popular
Matlab program.

To start IPython with matplotlib support, use the ``--matplotlib`` switch. If
IPython is already running, you can run the :magic:`matplotlib` magic.  If no
arguments are given, IPython will automatically detect your choice of
matplotlib backend.  You can also request a specific backend with
``%matplotlib backend``, where ``backend`` must be one of: 'tk', 'qt', 'wx',
'gtk', 'osx'.  In the web notebook and Qt console, 'inline' is also a valid
backend value, which produces static figures inlined inside the application
window instead of matplotlib's interactive figures that live in separate
windows.

.. _interactive_demos:

Interactive demos with IPython
==============================

IPython ships with a basic system for running scripts interactively in
sections, useful when presenting code to audiences. A few tags embedded
in comments (so that the script remains valid Python code) divide a file
into separate blocks, and the demo can be run one block at a time, with
IPython printing (with syntax highlighting) the block before executing
it, and returning to the interactive prompt after each block. The
interactive namespace is updated after each block is run with the
contents of the demo's namespace.

This allows you to show a piece of code, run it and then execute
interactively commands based on the variables just created. Once you
want to continue, you simply execute the next block of the demo. The
following listing shows the markup necessary for dividing a script into
sections for execution as a demo:

.. literalinclude:: ../../../examples/IPython Kernel/example-demo.py
    :language: python

In order to run a file as a demo, you must first make a Demo object out
of it. If the file is named myscript.py, the following code will make a
demo::

    from IPython.lib.demo import Demo

    mydemo = Demo('myscript.py')

This creates the mydemo object, whose blocks you run one at a time by
simply calling the object with no arguments. Then call it to run each step
of the demo::

    mydemo()

Demo objects can be
restarted, you can move forward or back skipping blocks, re-execute the
last block, etc. See the :mod:`IPython.lib.demo` module and the
:class:`~IPython.lib.demo.Demo` class for details.

Limitations: These demos are limited to
fairly simple uses. In particular, you cannot break up sections within
indented code (loops, if statements, function definitions, etc.)
Supporting something like this would basically require tracking the
internal execution state of the Python interpreter, so only top-level
divisions are allowed. If you want to be able to open an IPython
instance at an arbitrary point in a program, you can use IPython's
:ref:`embedding facilities <Embedding>`.

.. include:: ../links.txt
=================
Python vs IPython
=================

This document is meant to highlight the main differences between the Python
language and what are the specific constructs you can do only in IPython.

Unless expressed otherwise all of the constructs you will see here will raise a
``SyntaxError`` if run in a pure Python shell, or if executing in a Python
script.

Each of these features is described more in detail in the further parts of the documentation.


Quick overview:
===============


All the following constructs are valid IPython syntax:

.. code-block:: ipython

    In [1]: ?

.. code-block:: ipython

    In [1]: ?object


.. code-block:: ipython

    In [1]: object?

.. code-block:: ipython

    In [1]: *pattern*?

.. code-block:: ipython

    In [1]: %shell like --syntax

.. code-block:: ipython

    In [1]: !ls

.. code-block:: ipython

    In [1]: my_files = !ls ~/
    In [1]: for i, file in enumerate(my_files):
       ...:     raw = !echo $file
       ...:     !echo {file[0].upper()} $raw


.. code-block:: ipython

    In [1]: %%perl magic --function
       ...: @months = ("July", "August", "September");
       ...: print $months[0];


Each of these constructs is compiled by IPython into valid python code and will
do most of the time what you expect it will do. Let's see each of these examples
in more detail.


Accessing help
==============

As IPython is mostly an interactive shell, the question mark is a simple
shortcut to get help. A question mark alone will bring up the IPython help:

.. code-block:: ipython

    In [1]: ?

    IPython -- An enhanced Interactive Python
    =========================================

    IPython offers a combination of convenient shell features, special commands
    and a history mechanism for both input (command history) and output (results
    caching, similar to Mathematica). It is intended to be a fully compatible
    replacement for the standard Python interpreter, while offering vastly
    improved functionality and flexibility.

    At your system command line, type 'ipython -h' to see the command line
    options available. This document only describes interactive features.

    MAIN FEATURES
    -------------
    ...

A single question mark before or after an object available in the current
namespace will show help relative to this object:

.. code-block:: ipython

    In [6]: object?
    Docstring: The most base type
    Type:      type


A double question mark will try to pull out more information about the object,
and if possible display the python source code of this object.

.. code-block:: ipython

    In[1]: import collections
    In[2]: collections.Counter??

    Init signature: collections.Counter(*args, **kwds)
    Source:
    class Counter(dict):
        '''Dict subclass for counting hashable items.  Sometimes called a bag
        or multiset.  Elements are stored as dictionary keys and their counts
        are stored as dictionary values.

        >>> c = Counter('abcdeabcdabcaba')  # count elements from a string

        >>> c.most_common(3)                # three most common elements
        [('a', 5), ('b', 4), ('c', 3)]
        >>> sorted(c)                       # list all unique elements
        ['a', 'b', 'c', 'd', 'e']
        >>> ''.join(sorted(c.elements()))   # list elements with repetitions
        'aaaaabbbbcccdde'
        ...



If you are looking for an object, the use of wildcards ``*`` in conjunction
with a question mark will allow you to search the current namespace for objects with
matching names:

.. code-block:: ipython

    In [24]: *int*?
    FloatingPointError
    int
    print


Shell Assignment
================


When doing interactive computing it is a common need to access the underlying shell.
This is doable through the use of the exclamation mark ``!`` (or bang).

This allows to execute simple commands when present in beginning of the line:

.. code-block:: ipython

    In[1]: !pwd
    /User/home/

Change directory:

.. code-block:: ipython

    In[1]: !cd /var/etc

Or edit file:

.. code-block:: ipython

    In[1]: !mvim myfile.txt


The line after the bang can call any program installed in the underlying
shell, and support variable expansion in the form of ``$variable`` or ``{variable}``.
The later form of expansion supports arbitrary python expressions:

.. code-block:: ipython

    In[1]: file = 'myfile.txt'

    In[2]: !mv $file {file.upper()}


The bang (``!``) can also be present on the right hand side of an assignment, just
after the equal sign, or separated from it by a white space. In this case the
standard output of the command after the bang will be split out into lines
in a list-like object and assigned to the left hand side.

This allows you, for example, to put the list of files of the current working directory in a variable:

.. code-block:: ipython

    In[1]: my_files = !ls


You can combine the different possibilities in for loops, conditions, functions...:

.. code-block:: ipython

    my_files = !ls ~/
    for i, file in enumerate(my_files):
        raw = !echo $backup $file
        !cp $file {file.split('.')[0] + '.bak'}


Magics
------

Magic functions (magics) are often present in the form of shell-like syntax, but they are
python functions under the hood. The syntax and assignment possibilities are
similar to the one with the bang (``!``) syntax, but with more flexibility and
power. Magic functions start with a percent sign (``%``) or double percent signs (``%%``).

A magic call with a single percent sign will act only on one line:

.. code-block:: ipython

    In[1]: %xmode
    Exception reporting mode: Verbose

Magics support assignment:

.. code-block:: ipython

    In [1]: results = %timeit -r1 -n1 -o list(range(1000))
    1 loops, best of 1: 21.1 µs per loop

    In [2]: results
    Out[2]: <TimeitResult : 1 loops, best of 1: 21.1 µs per loop>

Magics with double percent signs (``%%``) can spread over multiple lines, but they do not support assignments:

.. code-block:: ipython

    In[1]: %%bash
    ...  : echo "My shell is:" $SHELL
    ...  : echo "My disk usage is:"
    ...  : df -h
    My shell is: /usr/local/bin/bash
    My disk usage is:
    Filesystem      Size   Used  Avail Capacity  iused   ifree %iused  Mounted on
    /dev/disk1     233Gi  216Gi   16Gi    94% 56788108 4190706   93%   /
    devfs          190Ki  190Ki    0Bi   100%      656       0  100%   /dev
    map -hosts       0Bi    0Bi    0Bi   100%        0       0  100%   /net
    map auto_home    0Bi    0Bi    0Bi   100%        0       0  100%   /hom
.. _autoawait:

Asynchronous in REPL: Autoawait
===============================

.. note::

   This feature is experimental and behavior can change between python and
   IPython version without prior deprecation.

Starting with IPython 7.0, and when using Python 3.6 and above, IPython offer the
ability to run asynchronous code from the REPL. Constructs which are
:exc:`SyntaxError` s in the Python REPL can be used seamlessly in IPython.

The examples given here are for terminal IPython, running async code in a
notebook interface or any other frontend using the Jupyter protocol needs
IPykernel version 5.0 or above. The details of how async code runs in IPykernel
will differ between IPython, IPykernel and their versions.

When a supported library is used, IPython will automatically allow Futures and
Coroutines in the REPL to be ``await`` ed. This will happen if an :ref:`await
<await>` (or any other async constructs like async-with, async-for) is used at
top level scope, or if any structure valid only in `async def
<https://docs.python.org/3/reference/compound_stmts.html#async-def>`_ function
context are present. For example, the following being a syntax error in the
Python REPL::

    Python 3.6.0 
    [GCC 4.2.1]
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import aiohttp
    >>> session = aiohttp.ClientSession()
    >>> result = session.get('https://api.github.com')
    >>> response = await result
      File "<stdin>", line 1
        response = await result
                              ^
    SyntaxError: invalid syntax

Should behave as expected in the IPython REPL::

    Python 3.6.0
    Type 'copyright', 'credits' or 'license' for more information
    IPython 7.0.0 -- An enhanced Interactive Python. Type '?' for help.

    In [1]: import aiohttp
       ...: session = aiohttp.ClientSession()
       ...: result = session.get('https://api.github.com')

    In [2]: response = await result
    <pause for a few 100s ms>

    In [3]: await response.json()
    Out[3]:
    {'authorizations_url': 'https://api.github.com/authorizations',
     'code_search_url': 'https://api.github.com/search/code?q={query}...',
    ...
    }


You can use the ``c.InteractiveShell.autoawait`` configuration option and set it
to :any:`False` to deactivate automatic wrapping of asynchronous code. You can
also use the :magic:`%autoawait` magic to toggle the behavior at runtime::

    In [1]: %autoawait False

    In [2]: %autoawait
    IPython autoawait is `Off`, and set to use `asyncio`



By default IPython will assume integration with Python's provided
:mod:`asyncio`, but integration with other libraries is provided. In particular
we provide experimental integration with the ``curio`` and ``trio`` library.

You can switch the current integration by using the
``c.InteractiveShell.loop_runner`` option or the ``autoawait <name
integration>`` magic.

For example::

    In [1]: %autoawait trio

    In [2]: import trio

    In [3]: async def child(i):
       ...:     print("   child %s goes to sleep"%i)
       ...:     await trio.sleep(2)
       ...:     print("   child %s wakes up"%i)

    In [4]: print('parent start')
       ...: async with trio.open_nursery() as n:
       ...:     for i in range(5):
       ...:         n.spawn(child, i)
       ...: print('parent end')
    parent start
       child 2 goes to sleep
       child 0 goes to sleep
       child 3 goes to sleep
       child 1 goes to sleep
       child 4 goes to sleep
       <about 2 seconds pause>
       child 2 wakes up
       child 1 wakes up
       child 0 wakes up
       child 3 wakes up
       child 4 wakes up
    parent end


In the above example, ``async with`` at top level scope is a syntax error in
Python.

Using this mode can have unexpected consequences if used in interaction with
other features of IPython and various registered extensions. In particular if
you are a direct or indirect user of the AST transformers, these may not apply
to your code.

When using command line IPython, the default loop (or runner) does not process
in the background, so top level asynchronous code must finish for the REPL to
allow you to enter more code. As with usual Python semantics, the awaitables are
started only when awaited for the first time. That is to say, in first example,
no network request is done between ``In[1]`` and ``In[2]``.


Effects on IPython.embed()
--------------------------

IPython core being asynchronous, the use of ``IPython.embed()`` will now require
a loop to run. By default IPython will use a fake coroutine runner which should
allow ``IPython.embed()`` to be nested. Though this will prevent usage of the
:magic:`%autoawait` feature when using IPython embed. 

You can set a coroutine runner explicitly for ``embed()`` if you want to run
asynchronous code, though the exact behavior is undefined.

Effects on Magics
-----------------

A couple of magics (``%%timeit``, ``%timeit``, ``%%time``, ``%%prun``) have not
yet been updated to work with asynchronous code and will raise syntax errors
when trying to use top-level ``await``. We welcome any contribution to help fix
those, and extra cases we haven't caught yet. We hope for better support in Core
Python for top-level Async code.

Internals
---------

As running asynchronous code is not supported in interactive REPL (as of Python
3.7) we have to rely to a number of complex workarounds and heuristics to allow
this to happen. It is interesting to understand how this works in order to
comprehend potential bugs, or provide a custom runner.

Among the many approaches that are at our disposition, we find only one that
suited out need. Under the hood we use the code object from a async-def function
and run it in global namespace after modifying it to not create a new
``locals()`` scope::

    async def inner_async():
        locals().update(**global_namespace)
        #
        # here is user code
        #
        return last_user_statement
    codeobj = modify(inner_async.__code__)
    coroutine = eval(codeobj, user_ns)
    display(loop_runner(coroutine))



The first thing you'll notice is that unlike classical ``exec``, there is only
one namespace. Second, user code runs in a function scope, and not a module
scope.

On top of the above there are significant modification to the AST of
``function``, and ``loop_runner`` can be arbitrary complex. So there is a
significant overhead to this kind of code.

By default the generated coroutine function will be consumed by Asyncio's
``loop_runner = asyncio.get_evenloop().run_until_complete()`` method if
``async`` mode is deemed necessary, otherwise the coroutine will just be
exhausted in a simple runner. It is possible, though, to change the default
runner.

A loop runner is a *synchronous*  function responsible from running a coroutine
object.

The runner is responsible for ensuring that ``coroutine`` runs to completion,
and it should return the result of executing the coroutine. Let's write a
runner for ``trio`` that print a message when used as an exercise, ``trio`` is
special as it usually prefers to run a function object and make a coroutine by
itself, we can get around this limitation by wrapping it in an async-def without
parameters and passing this value to ``trio``::


    In [1]: import trio
       ...: from types import CoroutineType
       ...:
       ...: def trio_runner(coro:CoroutineType):
       ...:     print('running asynchronous code')
       ...:     async def corowrap(coro):
       ...:         return await coro
       ...:     return trio.run(corowrap, coro)

We can set it up by passing it to ``%autoawait``::

    In [2]: %autoawait trio_runner

    In [3]: async def async_hello(name):
       ...:     await trio.sleep(1)
       ...:     print(f'Hello {name} world !')
       ...:     await trio.sleep(1)

    In [4]: await async_hello('async')
    running asynchronous code
    Hello async world !


Asynchronous programming in python (and in particular in the REPL) is still a
relatively young subject. We expect some code to not behave as you expect, so
feel free to contribute improvements to this codebase and give us feedback.

We invite you to thoroughly test this feature and report any unexpected behavior
as well as propose any improvement.

Using Autoawait in a notebook (IPykernel)
-----------------------------------------

Update ipykernel to version 5.0 or greater::

   pip install ipykernel ipython --upgrade
   # or
   conda install ipykernel ipython --upgrade

This should automatically enable :magic:`autoawait` integration. Unlike
terminal IPython, all code runs on ``asyncio`` eventloop, so creating a loop by
hand will not work, including with magics like :magic:`%run` or other
frameworks that create the eventloop themselves. In cases like these you can
try to use projects like `nest_asyncio
<https://github.com/erdewit/nest_asyncio>`_ and follow `this discussion
<https://github.com/jupyter/notebook/issues/3397#issuecomment-419386811>`_

Difference between terminal IPython and IPykernel
-------------------------------------------------

The exact asynchronous code running behavior varies between Terminal IPython and
IPykernel. The root cause of this behavior is due to IPykernel having a
*persistent* `asyncio` loop running, while Terminal IPython starts and stops a
loop for each code block. This can lead to surprising behavior in some cases if
you are used to manipulating asyncio loop yourself, see for example
:ghissue:`11303` for a longer discussion but here are some of the astonishing
cases.

This behavior is an implementation detail, and should not be relied upon. It can
change without warnings in future versions of IPython.

In terminal IPython a loop is started for each code blocks only if there is top
level async code::

   $ ipython
   In [1]: import asyncio
      ...: asyncio.get_event_loop()
   Out[1]: <_UnixSelectorEventLoop running=False closed=False debug=False>

   In [2]:

   In [2]: import asyncio
      ...: await asyncio.sleep(0)
      ...: asyncio.get_event_loop()
   Out[2]: <_UnixSelectorEventLoop running=True closed=False debug=False>

See that ``running`` is ``True`` only in the case were we ``await sleep()``

In a Notebook, with ipykernel the asyncio eventloop is always running::

   $ jupyter notebook
   In [1]: import asyncio
      ...: loop1 = asyncio.get_event_loop()
      ...: loop1
   Out[1]: <_UnixSelectorEventLoop running=True closed=False debug=False>

   In [2]: loop2 = asyncio.get_event_loop()
      ...: loop2
   Out[2]: <_UnixSelectorEventLoop running=True closed=False debug=False>

   In [3]: loop1 is loop2
   Out[3]: True

In Terminal IPython background tasks are only processed while the foreground
task is running, if and only if the foreground task is async::

   $ ipython
   In [1]: import asyncio
      ...:
      ...: async def repeat(msg, n):
      ...:     for i in range(n):
      ...:         print(f"{msg} {i}")
      ...:         await asyncio.sleep(1)
      ...:     return f"{msg} done"
      ...:
      ...: asyncio.ensure_future(repeat("background", 10))
   Out[1]: <Task pending coro=<repeat() running at <ipython-input-1-02d0ef250fe7>:3>>

   In [2]: await asyncio.sleep(3)
   background 0
   background 1
   background 2
   background 3

   In [3]: import time
   ...: time.sleep(5)

   In [4]: await asyncio.sleep(3)
   background 4
   background 5
   background 6g

In a Notebook, QtConsole, or any other frontend using IPykernel, background
tasks should behave as expected.
=======================
Built-in magic commands
=======================

.. note::

    To Jupyter users: Magics are specific to and provided by the IPython kernel.
    Whether Magics are available on a kernel is a decision that is made by
    the kernel developer on a per-kernel basis. To work properly, Magics must
    use a syntax element which is not valid in the underlying language. For
    example, the IPython kernel uses the `%` syntax element for Magics as `%`
    is not a valid unary operator in Python. However, `%` might have meaning in
    other languages.

Here is the help auto-generated from the docstrings of all the available Magics
functions that IPython ships with. 

You can create and register your own Magics with IPython. You can find many user
defined Magics on `PyPI <https://pypi.io>`_. Feel free to publish your own and
use the ``Framework :: IPython`` trove classifier. 


.. include:: magics-generated.txt
.. _tutorial:

======================
Introducing IPython
======================

You don't need to know anything beyond Python to start using IPython – just type
commands as you would at the standard Python prompt. But IPython can do much
more than the standard prompt. Some key features are described here. For more
information, check the :ref:`tips page <tips>`, or look at examples in the
`IPython cookbook <https://github.com/ipython/ipython/wiki/Cookbook%3A-Index>`_.

If you haven't done that yet see :ref:`how to install ipython <install>`.

If you've never used Python before, you might want to look at `the official
tutorial <http://docs.python.org/tutorial/>`_ or an alternative, `Dive into
Python <https://www.diveinto.org/python3/table-of-contents.html>`_.

Start IPython by issuing the ``ipython`` command from your shell, you should be
greeted by the following::

    Python 3.6.0
    Type 'copyright', 'credits' or 'license' for more information
    IPython 6.0.0.dev -- An enhanced Interactive Python. Type '?' for help.

    In [1]:


Unlike the Python REPL, you will see that the input prompt is ``In [N]:``
instead of ``>>>``. The number ``N`` in the prompt will be used later in this
tutorial but should usually not impact the computation.

You should be able to type single line expressions and press enter to evaluate
them. If an expression is incomplete, IPython will automatically detect this and
add a new line when you press :kbd:`Enter` instead of executing right away.

Feel free to explore multi-line text input. Unlike many other REPLs, with
IPython you can use the up and down arrow keys when editing multi-line
code blocks.

Here is an example of a longer interaction with the IPython REPL,
which we often refer to as an IPython *session* ::

    In [1]: print('Hello IPython')
    Hello IPython

    In [2]: 21 * 2
    Out[2]: 42

    In [3]: def say_hello(name):
       ...:     print('Hello {name}'.format(name=name))
       ...:

We won't get into details right now, but you may notice a few differences to
the standard Python REPL. First, your code should be syntax-highlighted as you
type. Second, you will see that some results will have an ``Out[N]:`` prompt,
while some other do not. We'll come to this later.

Depending on the exact command you are typing you might realize that sometimes
:kbd:`Enter` will add a new line, and sometimes it will execute the current
statement. IPython tries to guess what you are doing, so most of the time you
should not have to care. Though if by any chance IPython does not the right
thing you can force execution of the current code block by pressing in sequence
:kbd:`Esc` and :kbd:`Enter`. You can also force the insertion of a new line at
the position of the cursor by using :kbd:`Ctrl-o`.

The four most helpful commands
==============================

The four most helpful commands, as well as their brief description, is shown
to you in a banner, every time you start IPython:

==========    =========================================================
command       description
==========    =========================================================
?             Introduction and overview of IPython's features.
%quickref     Quick reference.
help          Python's own help system.
object?       Details about 'object', use 'object??' for extra details.
==========    =========================================================

Tab completion
==============

Tab completion, especially for attributes, is a convenient way to explore the
structure of any object you're dealing with. Simply type ``object_name.<TAB>``
to view the object's attributes. Besides Python objects and keywords, tab
completion also works on file and directory names.

Starting with IPython 6.0, if ``jedi`` is installed, IPython will try to pull
completions from Jedi as well. This allows to not only inspect currently
existing objects, but also to infer completion statically without executing
code. There is nothing particular needed to get this to work, simply use tab
completion on more complex expressions like the following::

    >>> data = ['Number of users', 123456]
    ... data[0].<tab>

IPython and Jedi will be able to infer that ``data[0]`` is actually a string
and should show relevant completions like ``upper()``, ``lower()`` and other
string methods. You can use the :kbd:`Tab` key to cycle through completions,
and while a completion is highlighted, its type will be shown as well.
When the type of the completion is a function, the completer will also show the
signature of the function when highlighted.

Exploring your objects
======================

Typing ``object_name?`` will print all sorts of details about any object,
including docstrings, function definition lines (for call arguments) and
constructor details for classes. To get specific information on an object, you
can use the magic commands ``%pdoc``, ``%pdef``, ``%psource`` and ``%pfile``

.. _magics_explained:

Magic functions
===============

IPython has a set of predefined 'magic functions' that you can call with a
command line style syntax.  There are two kinds of magics, line-oriented and
cell-oriented.  **Line magics** are prefixed with the ``%`` character and work
much like OS command-line calls: they get as an argument the rest of the line,
where arguments are passed without parentheses or quotes. **Lines magics** can
return results and can be used in the right hand side of an assignment.  **Cell
magics** are prefixed with a double ``%%``, and they are functions that get as
an argument not only the rest of the line, but also the lines below it in a
separate argument.

Magics are useful as convenient functions where Python syntax is not the most
natural one, or when one want to embed invalid python syntax in their work flow. 

The following examples show how to call the built-in :magic:`timeit` magic, both
in line and cell mode::

      In [1]: %timeit range(1000)
      100000 loops, best of 3: 7.76 us per loop

      In [2]: %%timeit x = range(10000)
      ...: max(x)
      ...: 
      1000 loops, best of 3: 223 us per loop

The built-in magics include:

- Functions that work with code: :magic:`run`, :magic:`edit`, :magic:`save`,
  :magic:`macro`, :magic:`recall`, etc.

- Functions which affect the shell: :magic:`colors`, :magic:`xmode`,
  :magic:`automagic`, etc.

- Other functions such as :magic:`reset`, :magic:`timeit`,
  :cellmagic:`writefile`, :magic:`load`, or :magic:`paste`.

You can always call magics using the ``%`` prefix, and if you're calling a line
magic on a line by itself, as long as the identifier is not defined in your
namespace, you can omit even that::

    run thescript.py

You can toggle this behavior by running the :magic:`automagic` magic.  Cell
magics must always have the ``%%`` prefix.

A more detailed explanation of the magic system can be obtained by calling
``%magic``, and for more details on any magic function, call ``%somemagic?`` to
read its docstring. To see all the available magic functions, call
``%lsmagic``.

.. seealso::

    The :ref:`magic` section of the documentation goes more in depth into how
    the magics works and how to define your own, and :doc:`magics` for a list of
    built-in magics.

    `Cell magics`_ example notebook

Running and Editing
-------------------

The :magic:`run` magic command allows you to run any python script and load all
of its data directly into the interactive namespace. Since the file is re-read
from disk each time, changes you make to it are reflected immediately (unlike
imported modules, which have to be specifically reloaded). IPython also includes
:ref:`dreload <dreload>`, a recursive reload function.

``%run`` has special flags for timing the execution of your scripts (-t), or
for running them under the control of either Python's pdb debugger (-d) or
profiler (-p).

The :magic:`edit` command gives a reasonable approximation of multi-line editing,
by invoking your favorite editor on the spot. IPython will execute the
code you type in there as if it were typed interactively. Note that for
:magic:`edit` to work, the call to startup your editor has to be a blocking
call. In a GUI environment, your editor likely will have such an option.

Debugging
---------

After an exception occurs, you can call :magic:`debug` to jump into the Python
debugger (pdb) and examine the problem. Alternatively, if you call :magic:`pdb`,
IPython will automatically start the debugger on any uncaught exception. You can
print variables, see code, execute statements and even walk up and down the call
stack to track down the true source of the problem. This can be an efficient way
to develop and debug code, in many cases eliminating the need for print
statements or external debugging tools.

You can also step through a program from the beginning by calling
``%run -d theprogram.py``.

History
=======

IPython stores both the commands you enter, and the results it produces. You
can easily go through previous commands with the up- and down-arrow keys, or
access your history in more sophisticated ways.

Input and output history are kept in variables called ``In`` and ``Out``, keyed
by the prompt numbers, e.g. ``In[4]``. The last three objects in output history
are also kept in variables named ``_``, ``__`` and ``___``.

You can use the ``%history`` magic function to examine past input and output.
Input history from previous sessions is saved in a database, and IPython can be
configured to save output history.

Several other magic functions can use your input history, including ``%edit``, 
``%rerun``, ``%recall``, ``%macro``, ``%save`` and ``%pastebin``. You can use a
standard format to refer to lines::

    %pastebin 3 18-20 ~1/1-5
    
This will take line 3 and lines 18 to 20 from the current session, and lines
1-5 from the previous session.

System shell commands
=====================

To run any command at the system shell, simply prefix it with ``!``, e.g.::

    !ping www.bbc.co.uk
    
You can capture the output into a Python list, e.g.: ``files = !ls``. To pass
the values of Python variables or expressions to system commands, prefix them
with $: ``!grep -rF $pattern ipython/*`` or wrap in `{braces}`. See :ref:`our
shell section <system_shell_access>` for more details.

Define your own system aliases
------------------------------

It's convenient to have aliases to the system commands you use most often. This
allows you to work seamlessly from inside IPython with the same commands you are
used to in your system shell. IPython comes with some pre-defined aliases and a
complete system for changing directories, both via a stack (see :magic:`pushd`,
:magic:`popd` and :magic:`dhist`) and via direct :magic:`cd`. The latter keeps a
history of visited directories and allows you to go to any previously visited
one.


Configuration
=============

Much of IPython can be tweaked through :doc:`configuration </config/intro>`.
To get started, use the command ``ipython profile create`` to produce the
default config files. These will be placed in
:file:`~/.ipython/profile_default`, and contain comments explaining
what the various options do.

Profiles allow you to use IPython for different tasks, keeping separate config
files and history for each one. More details in :ref:`the profiles section
<profiles>`.

.. _startup_files:

Startup Files
-------------

If you want some code to be run at the beginning of every IPython session, the
easiest way is to add Python (.py) or IPython (.ipy) scripts to your
:file:`profile_default/startup/` directory. Files here will be executed as soon
as the IPython shell is constructed, before any other code or scripts you have
specified. The files will be run in order of their names, so you can control the
ordering with prefixes, like ``10-myimports.py``.

.. include:: ../links.txt
========
Tutorial
========

This section of IPython documentation will walk you through most of the IPython
functionality. You do not need to have any deep knowledge of Python to read this
tutorial, though some sections might make slightly more sense if you have already
done some work in the classic Python REPL.

.. note::
    
    Some part of this documentation are more than a decade old so might be out
    of date, we welcome any report of inaccuracy, and Pull Requests that make
    that up to date.

.. toctree::
   :maxdepth: 2
   :hidden:

   tutorial
   plotting
   reference
   shell
   autoawait
   tips
   python-ipython-diff
   magics

.. seealso::

    `A Qt Console for Jupyter <https://jupyter.org/qtconsole/>`__
    `The Jupyter Notebook <https://jupyter-notebook.readthedocs.io/en/latest/>`__
.. _plotting:

Rich Outputs
------------

One of the main feature of IPython when used as a kernel is its ability to
show rich output. This means that object that can be representing as image,
sounds, animation, (etc...) can be shown this way if the frontend support it.

In order for this to be possible, you need to use the ``display()`` function,
that should be available by default on IPython 5.4+ and 6.1+, or that you can
import with ``from IPython.display import display``. Then use ``display(<your
object>)`` instead of ``print()``, and if possible your object will be displayed
with a richer representation. In the terminal of course, there won't be much
difference as object are most of the time represented by text, but in notebook
and similar interface you will get richer outputs.


Plotting
--------

.. note::

    Starting with IPython 5.0 and matplotlib 2.0 you can avoid the use of
    IPython's specific magic and use
    ``matplotlib.pyplot.ion()``/``matplotlib.pyplot.ioff()`` which have the
    advantages of working outside of IPython as well.


One major feature of the IPython kernel is the ability to display plots that 
are the output of running code cells. The IPython kernel is designed to work 
seamlessly with the matplotlib_ plotting library to provide this functionality.

To set this up, before any plotting or import of matplotlib is performed you
must execute the ``%matplotlib``  :ref:`magic command <magics_explained>`. This
performs the necessary behind-the-scenes setup for IPython to work correctly
hand in hand with ``matplotlib``; it does *not*, however, actually execute any
Python ``import`` commands, that is, no names are added to the namespace.

If the ``%matplotlib`` magic is called without an argument, the
output of a plotting command is displayed using the default ``matplotlib``
backend in a separate window. Alternatively, the backend can be explicitly
requested using, for example::

  %matplotlib gtk

A particularly interesting backend, provided by IPython, is the ``inline``
backend.  This is available only for the Jupyter Notebook and the
Jupyter QtConsole.  It can be invoked as follows::

  %matplotlib inline

With this backend, the output of plotting commands is displayed *inline* within
frontends like the Jupyter notebook, directly below the code cell that produced
it. The resulting plots will then also be stored in the notebook document.

.. seealso::

    `Plotting with Matplotlib`_  example notebook


The matplotlib_ library also ships with ``%matplotlib notebook`` command that
allows interactive figures if your environment allows it.

See the matplotlib_ documentation for more information. 

.. include:: ../links.txt
.. _ipython_as_shell:

.. note::

   This page has been kept for historical reason. You most likely want to use
   `Xonsh <https://xon.sh/>`__ instead of this.


=========================
IPython as a system shell
=========================



Overview
========

It is possible to adapt IPython for system shell usage. In the past, IPython
shipped a special 'sh' profile for this purpose, but it had been quarantined
since 0.11 release, and in 1.0 it was removed altogether. Nevertheless, much
of this section relies on machinery which does not require a custom profile.

You can set up your own 'sh' :ref:`profile <Profiles>` to be different from
the default profile such that:

 * Prompt shows the current directory (see `Prompt customization`_)
 * Make system commands directly available (in alias table) by running the
   ``%rehashx`` magic. If you install new programs along your PATH, you might
   want to run ``%rehashx`` to update the alias table
 * turn ``%autocall`` to full mode


Environment variables
=====================

Rather than manipulating os.environ directly, you may like to use the magic
`%env` command.  With no arguments, this displays all environment variables
and values.  To get the value of a specific variable, use `%env var`.  To set
the value of a specific variable, use `%env foo bar`, `%env foo=bar`.  By
default values are considered to be strings so quoting them is unnecessary.
However, Python variables are expanded as usual in the magic command, so
`%env foo=$bar` means "set the environment variable foo to the value of the
Python variable `bar`".

Aliases
=======

Once you run ``%rehashx``, all of your $PATH has been loaded as IPython aliases,
so you should be able to type any normal system command and have it executed.
See ``%alias?``  and ``%unalias?`` for details on the alias facilities. See also
``%rehashx?`` for details on the mechanism used to load $PATH.

.. warning::

   See info at the top of the page. You most likely want to use
   `Xonsh <https://xon.sh/>`__ instead of this.

Directory management
====================

Since each command passed by IPython to the underlying system is executed
in a subshell which exits immediately, you can NOT use !cd to navigate
the filesystem.

IPython provides its own builtin ``%cd`` magic command to move in the
filesystem (the % is not required with automagic on). It also maintains
a list of visited directories (use ``%dhist`` to see it) and allows direct
switching to any of them. Type ``cd?`` for more details.

``%pushd``, ``%popd`` and ``%dirs`` are provided for directory stack handling.


Prompt customization
====================

See :ref:`custom_prompts`.


.. _string_lists:

String lists
============

String lists (IPython.utils.text.SList) are handy way to process output
from system commands. They are produced by ``var = !cmd`` syntax.

First, we acquire the output of 'ls -l'::

    [Q:doc/examples]|2> lines = !ls -l
     ==
    ['total 23',
     '-rw-rw-rw- 1 ville None 1163 Sep 30  2006 example-demo.py',
     '-rw-rw-rw- 1 ville None 1927 Sep 30  2006 example-embed-short.py',
     '-rwxrwxrwx 1 ville None 4606 Sep  1 17:15 example-embed.py',
     '-rwxrwxrwx 1 ville None 1017 Sep 30  2006 example-gnuplot.py',
     '-rwxrwxrwx 1 ville None  339 Jun 11 18:01 extension.py',
     '-rwxrwxrwx 1 ville None  113 Dec 20  2006 seteditor.py',
     '-rwxrwxrwx 1 ville None  245 Dec 12  2006 seteditor.pyc']

Now, let's take a look at the contents of 'lines' (the first number is
the list element number)::

    [Q:doc/examples]|3> lines
                    <3> SList (.p, .n, .l, .s, .grep(), .fields() available). Value:

    0: total 23
    1: -rw-rw-rw- 1 ville None 1163 Sep 30  2006 example-demo.py
    2: -rw-rw-rw- 1 ville None 1927 Sep 30  2006 example-embed-short.py
    3: -rwxrwxrwx 1 ville None 4606 Sep  1 17:15 example-embed.py
    4: -rwxrwxrwx 1 ville None 1017 Sep 30  2006 example-gnuplot.py
    5: -rwxrwxrwx 1 ville None  339 Jun 11 18:01 extension.py
    6: -rwxrwxrwx 1 ville None  113 Dec 20  2006 seteditor.py
    7: -rwxrwxrwx 1 ville None  245 Dec 12  2006 seteditor.pyc

Now, let's filter out the 'embed' lines::

    [Q:doc/examples]|4> l2 = lines.grep('embed',prune=1)
    [Q:doc/examples]|5> l2
                    <5> SList (.p, .n, .l, .s, .grep(), .fields() available). Value:

    0: total 23
    1: -rw-rw-rw- 1 ville None 1163 Sep 30  2006 example-demo.py
    2: -rwxrwxrwx 1 ville None 1017 Sep 30  2006 example-gnuplot.py
    3: -rwxrwxrwx 1 ville None  339 Jun 11 18:01 extension.py
    4: -rwxrwxrwx 1 ville None  113 Dec 20  2006 seteditor.py
    5: -rwxrwxrwx 1 ville None  245 Dec 12  2006 seteditor.pyc

Now, we want strings having just file names and permissions::

    [Q:doc/examples]|6> l2.fields(8,0)
                    <6> SList (.p, .n, .l, .s, .grep(), .fields() available). Value:

    0: total
    1: example-demo.py -rw-rw-rw-
    2: example-gnuplot.py -rwxrwxrwx
    3: extension.py -rwxrwxrwx
    4: seteditor.py -rwxrwxrwx
    5: seteditor.pyc -rwxrwxrwx

Note how the line with 'total' does not raise IndexError.

If you want to split these (yielding lists), call fields() without
arguments::

    [Q:doc/examples]|7> _.fields()
                    <7>
    [['total'],
     ['example-demo.py', '-rw-rw-rw-'],
     ['example-gnuplot.py', '-rwxrwxrwx'],
     ['extension.py', '-rwxrwxrwx'],
     ['seteditor.py', '-rwxrwxrwx'],
     ['seteditor.pyc', '-rwxrwxrwx']]

If you want to pass these separated with spaces to a command (typical
for lists if files), use the .s property::


    [Q:doc/examples]|13> files = l2.fields(8).s
    [Q:doc/examples]|14> files
                    <14> 'example-demo.py example-gnuplot.py extension.py seteditor.py seteditor.pyc'
    [Q:doc/examples]|15> ls $files
    example-demo.py  example-gnuplot.py  extension.py  seteditor.py  seteditor.pyc

SLists are inherited from normal Python lists, so every list method is
available::

    [Q:doc/examples]|21> lines.append('hey')


Real world example: remove all files outside version control
------------------------------------------------------------

First, capture output of "hg status"::

    [Q:/ipython]|28> out = !hg status
     ==
    ['M IPython\\extensions\\ipy_kitcfg.py',
     'M IPython\\extensions\\ipy_rehashdir.py',
    ...
     '? build\\lib\\IPython\\Debugger.py',
     '? build\\lib\\IPython\\extensions\\InterpreterExec.py',
     '? build\\lib\\IPython\\extensions\\InterpreterPasteInput.py',
    ...

(lines starting with ? are not under version control).

::

    [Q:/ipython]|35> junk = out.grep(r'^\?').fields(1)
    [Q:/ipython]|36> junk
                <36> SList (.p, .n, .l, .s, .grep(), .fields() availab
    ...
    10: build\bdist.win32\winexe\temp\_ctypes.py
    11: build\bdist.win32\winexe\temp\_hashlib.py
    12: build\bdist.win32\winexe\temp\_socket.py

Now we can just remove these files by doing 'rm $junk.s'.

The .n, .s, .p properties
-------------------------

Properties of `SList <https://ipython.readthedocs.io/en/stable/api/generated/IPython.utils.text.html?highlight=SList#IPython.utils.text.SList>`_ wrapper
provide a convenient ways to use contained text in different formats:

* ``.n`` returns (original) string with lines separated by a newline
* ``.s`` returns string with lines separated by single space (for
  convenient passing to system commands)
* ``.p`` returns list of "path" objects from detected file names

.. error::

   You went too far scroll back up. You most likely want to use
   `Xonsh <https://xon.sh/>`__ instead of this.

.. _core_developer_guide:

=================================
Guide for IPython core Developers
=================================

This guide documents the development of IPython itself.  Alternatively,
developers of third party tools and libraries that use IPython should see the
:doc:`../development/index`.


For instructions on how to make a developer install see :ref:`devinstall`.

Backporting Pull requests
=========================

All pull requests should usually be made against ``master``, if a Pull Request
need to be backported to an earlier release; then it should be tagged with the
correct ``milestone``.

If you tag a pull request with a milestone **before** merging the pull request,
and the base ref is ``master``, then our backport bot should automatically create
a corresponding pull-request that backport on the correct branch.

If you have write access to the IPython repository you can also just mention the
**backport bot** to do the work for you. The bot is evolving so instructions may
be different. At the time of this writing you can use::

    @meeseeksdev[bot] backport [to] <branchname>

The bot will attempt to backport the current pull-request and issue a PR if
possible. 

.. note::

    The ``@`` and ``[bot]`` when mentioning the bot should be optional and can
    be omitted.

If the pull request cannot be automatically backported, the bot should tell you
so on the PR and apply a "Need manual backport" tag to the origin PR.

.. _release_process:

IPython release process
=======================

This document contains the process that is used to create an IPython release.

Conveniently, the ``release`` script in the ``tools`` directory of the ``IPython``
repository automates most of the release process. This document serves as a
handy reminder and checklist for the release manager.

During the release process, you might need the extra following dependencies:

 - ``keyring`` to access your GitHub authentication tokens
 - ``graphviz`` to generate some graphs in the documentation
 - ``ghpro`` to generate the stats

Make sure you have all the required dependencies to run the tests as well.

You can try to ``source tools/release_helper.sh`` when releasing via bash, it 
should guide you through most of the process.


1. Set Environment variables
----------------------------

Set environment variables to document previous release tag, current
release milestone, current release version, and git tag.

These variables may be used later to copy/paste as answers to the script
questions instead of typing the appropriate command when the time comes. These
variables are not used by the scripts directly; therefore, there is no need to
``export`` them. The format for bash is as follows, but note that these values
are just an example valid only for the 5.0 release; you'll need to update them
for the release you are actually making::

    PREV_RELEASE=4.2.1
    MILESTONE=5.0
    VERSION=5.0.0
    BRANCH=master

For `reproducibility of builds <https://reproducible-builds.org/specs/source-date-epoch/>`_,
we recommend setting ``SOURCE_DATE_EPOCH`` prior to running the build; record the used value
of ``SOURCE_DATE_EPOCH`` as it may not be available from build artifact. You
should be able to use ``date +%s`` to get a formatted timestamp::

    SOURCE_DATE_EPOCH=$(date +%s)


2. Create GitHub stats and finish release note
----------------------------------------------

.. note::

    This step is optional if making a Beta or RC release.

.. note::

    Before generating the GitHub stats, verify that all closed issues and pull
    requests have `appropriate milestones
    <https://github.com/ipython/ipython/wiki/Dev:-GitHub-workflow#milestones>`_.
    `This search
    <https://github.com/ipython/ipython/issues?q=is%3Aclosed+no%3Amilestone+is%3Aissue>`_
    should return no results before creating the GitHub stats.

If a major release:

    - merge any pull request notes into what's new::

          python tools/update_whatsnew.py

    - update ``docs/source/whatsnew/development.rst``, to ensure it covers
      the major release features

    - move the contents of ``development.rst`` to ``versionX.rst`` where ``X`` is
      the numerical release version

    - generate summary of GitHub contributions, which can be done with::

          python tools/github_stats.py --milestone $MILESTONE > stats.rst

      which may need some manual cleanup of ``stats.rst``. Add the cleaned
      ``stats.rst`` results to ``docs/source/whatsnew/github-stats-X.rst``
      where ``X`` is the numerical release version (don't forget to add it to
      the git repository as well). If creating a major release, make a new
      ``github-stats-X.rst`` file; if creating a minor release, the content
      from ``stats.rst`` may simply be added to the top of an existing
      ``github-stats-X.rst`` file.

    - Edit ``docs/source/whatsnew/index.rst`` to list the new ``github-stats-X``
      file you just created.

    - You do not need to temporarily remove the first entry called
      ``development``, nor re-add it after the release, it will automatically be
      hidden when releasing a stable version of IPython (if ``_version_extra``
      in ``release.py`` is an empty string.

      Make sure that the stats file has a header or it won't be rendered in
      the final documentation.

To find duplicates and update `.mailmap`, use::

    git log --format="%aN <%aE>" $PREV_RELEASE... | sort -u -f

If a minor release you might need to do some of the above points manually, and
forward port the changes.

3. Make sure the repository is clean
------------------------------------

of any file that could be problematic.
   Remove all non-tracked files with:

   .. code::

       git clean -xfdi

   This will ask for confirmation before removing all untracked files. Make
   sure the ``dist/`` folder is clean to avoid any stale builds from
   previous build attempts.


4. Update the release version number
------------------------------------

Edit ``IPython/core/release.py`` to have the current version.

in particular, update version number and ``_version_extra`` content in
``IPython/core/release.py``.

Step 5 will validate your changes automatically, but you might still want to
make sure the version number matches pep440.

In particular, ``rc`` and ``beta`` are not separated by ``.`` or the ``sdist``
and ``bdist`` will appear as different releases. For example, a valid version
number for a release candidate (rc) release is: ``1.3rc1``. Notice that there
is no separator between the '3' and the 'r'. Check the environment variable
``$VERSION`` as well.

You will likely just have to modify/comment/uncomment one of the lines setting
``_version_extra``


5. Run the `tools/build_release` script
---------------------------------------

Running ``tools/build_release`` does all the file checking and building that
the real release script will do. This makes test installations, checks that
the build procedure runs OK, and tests other steps in the release process.

The ``build_release`` script will in particular verify that the version number
match PEP 440, in order to avoid surprise at the time of build upload.

We encourage creating a test build of the docs as well. 

6. Create and push the new tag
------------------------------

Commit the changes to release.py::

    git commit -am "release $VERSION" -S
    git push origin $BRANCH

(omit the ``-S`` if you are no signing the package)

Create and push the tag::

    git tag -am "release $VERSION" "$VERSION" -s
    git push origin $VERSION

(omit the ``-s`` if you are no signing the package)

Update release.py back to ``x.y-dev`` or ``x.y-maint`` commit and push::

    git commit -am "back to development" -S
    git push origin $BRANCH

(omit the ``-S`` if you are no signing the package)

Now checkout the tag we just made::

    git checkout $VERSION

7. Run the release script
-------------------------

Run the ``release`` script, this step requires having a current wheel, Python
>=3.4 and Python 2.7.::

    ./tools/release

This makes the tarballs and wheels, and puts them under the ``dist/``
folder. Be sure to test the ``wheels``  and the ``sdist`` locally before
uploading them to PyPI. We do not use an universal wheel as each wheel
installs an ``ipython2`` or ``ipython3`` script, depending on the version of
Python it is built for. Using an universal wheel would prevent this.

Check the shasum of files with::

    shasum -a 256 dist/*

and takes notes of them you might need them to update the conda-forge recipes.
Rerun the command and check the hash have not changed::

    ./tools/release
    shasum -a 256 dist/*

Use the following to actually upload the result of the build::

    ./tools/release upload

It should posts them to ``archive.ipython.org`` and to PyPI.

PyPI/Warehouse will automatically hide previous releases. If you are uploading
a non-stable version, make sure to log-in to PyPI and un-hide previous version.


8. Draft a short release announcement
-------------------------------------

The announcement should include:

- release highlights
- a link to the html version of the *What's new* section of the documentation
- a link to upgrade or installation tips (if necessary)

Post the announcement to the mailing list and or blog, and link from Twitter.

.. note::

    If you are doing a RC or Beta, you can likely skip the next steps.

9. Update milestones on GitHub
-------------------------------

These steps will bring milestones up to date:

- close the just released milestone
- open a new milestone for the next release (x, y+1), if the milestone doesn't
  exist already

10. Update the IPython website
------------------------------

The IPython website should document the new release:

- add release announcement (news, announcements)
- update current version and download links
- update links on the documentation page (especially if a major release)

11. Update readthedocs
----------------------

Make sure to update readthedocs and set the latest tag as stable, as well as
checking that previous release is still building under its own tag.

12. Update the Conda-Forge feedstock
------------------------------------

Follow the instructions on `the repository <https://github.com/conda-forge/ipython-feedstock>`_

13. Celebrate!
--------------

Celebrate the release and please thank the contributors for their work. Great
job!



Old Documentation
=================

Out of date documentation is still available and have been kept for archival purposes.

.. note::

  Developers documentation used to be on the IPython wiki, but are now out of
  date. The wiki is though still available for historical reasons: `Old IPython
  GitHub Wiki.  <https://github.com/ipython/ipython/wiki/Dev:-Index>`_
.. _kernel_install:

Installing the IPython kernel
=============================

.. seealso::

   :ref:`Installing Jupyter <jupyter:install>`
     The IPython kernel is the Python execution backend for Jupyter.

The Jupyter Notebook and other frontends automatically ensure that the IPython kernel is available.
However, if you want to use a kernel with a different version of Python, or in a virtualenv or conda environment,
you'll need to install that manually.

Kernels for Python 2 and 3
--------------------------

If you're running Jupyter on Python 3, you can set up a Python 2 kernel after
checking your version of pip is greater than 9.0::

    python2 -m pip --version

Then install with ::

    python2 -m pip install ipykernel
    python2 -m ipykernel install --user

Or using conda, create a Python 2 environment::

    conda create -n ipykernel_py2 python=2 ipykernel
    source activate ipykernel_py2    # On Windows, remove the word 'source'
    python -m ipykernel install --user

.. note::

    IPython 6.0 stopped support for Python 2, so
    installing IPython on Python 2 will give you an older version (5.x series).

If you're running Jupyter on Python 2 and want to set up a Python 3 kernel,
follow the same steps, replacing ``2`` with ``3``.

The last command installs a :ref:`kernel spec <jupyterclient:kernelspecs>` file
for the current python installation. Kernel spec files are JSON files, which
can be viewed and changed with a normal text editor.

.. _multiple_kernel_install:

Kernels for different environments
----------------------------------

If you want to have multiple IPython kernels for different virtualenvs or conda
environments, you will need to specify unique names for the kernelspecs.

Make sure you have ipykernel installed in your environment. If you are using
``pip`` to install ``ipykernel`` in a conda env, make sure ``pip`` is
installed:

.. sourcecode:: bash

    source activate myenv
    conda install pip
    conda install ipykernel # or pip install ipykernel

For example, using conda environments, install a ``Python (myenv)`` Kernel in a first
environment:

.. sourcecode:: bash

    source activate myenv
    python -m ipykernel install --user --name myenv --display-name "Python (myenv)"

And in a second environment, after making sure ipykernel is installed in it:

.. sourcecode:: bash

    source activate other-env
    python -m ipykernel install --user --name other-env --display-name "Python (other-env)"

The ``--name`` value is used by Jupyter internally. These commands will overwrite
any existing kernel with the same name. ``--display-name`` is what you see in
the notebook menus.

Using virtualenv or conda envs, you can make your IPython kernel in one env available to Jupyter in a different env. To do so, run ipykernel install from the kernel's env, with --prefix pointing to the Jupyter env:

.. sourcecode:: bash

    /path/to/kernel/env/bin/python -m ipykernel install --prefix=/path/to/jupyter/env --name 'python-my-env'

Note that this command will create a new configuration for the kernel in one of the preferred location (see ``jupyter --paths`` command for more details):

* system-wide (e.g. /usr/local/share),
* in Jupyter's env (sys.prefix/share),
* per-user (~/.local/share or ~/Library/share)

If you want to edit the kernelspec before installing it, you can do so in two steps.
First, ask IPython to write its spec to a temporary location:

.. sourcecode:: bash

    ipython kernel install --prefix /tmp

edit the files in /tmp/share/jupyter/kernels/python3 to your liking, then when you are ready, tell Jupyter to install it (this will copy the files into a place Jupyter will look):

.. sourcecode:: bash

    jupyter kernelspec install /tmp/share/jupyter/kernels/python3
.. _install_index:

============
Installation
============

.. toctree::
   :maxdepth: 3
   :hidden:


   install
   kernel_install



This sections will guide you through :ref:`installing IPython itself <install>`, and
installing :ref:`kernels for Jupyter <kernel_install>` if you wish to work with
multiple version of Python, or multiple environments.


Quick install reminder
~~~~~~~~~~~~~~~~~~~~~~

Here is a quick reminder of the commands needed for installation if you are
already familiar with IPython and are just searching to refresh your memory:

Install IPython:

.. code-block:: bash

    $ pip install ipython


Install and register an IPython kernel with Jupyter:


.. code-block:: bash

    $ python -m pip install ipykernel

    $ python -m ipykernel install [--user] [--name <machine-readable-name>] [--display-name <"User Friendly Name">]

for more help see 

.. code-block:: bash
    
    $ python -m ipykernel install  --help
    


.. seealso::

   `Installing Jupyter <https://jupyter.readthedocs.io/en/latest/install.html>`__
     The Notebook, nbconvert, and many other former pieces of IPython are now
     part of Project Jupyter.


.. _install:

Installing IPython
==================


IPython 6 requires Python ≥ 3.3. IPython 5.x can be installed on Python 2.


Quick Install 
-------------

With ``pip`` already installed :

.. code-block:: bash

    $ pip install ipython

This installs IPython as well as its dependencies.

If you want to use IPython with notebooks or the Qt console, you should also
install Jupyter ``pip install jupyter``.



Overview
--------

This document describes in detail the steps required to install IPython. For a
few quick ways to get started with package managers or full Python
distributions, see `the install page <https://ipython.org/install.html>`_ of the
IPython website.

Please let us know if you have problems installing IPython or any of its
dependencies.

IPython and most dependencies should be installed via :command:`pip`.
In many scenarios, this is the simplest method of installing Python packages.
More information about :mod:`pip` can be found on
`its PyPI page <https://pip.pypa.io>`__.


More general information about installing Python packages can be found in
`Python's documentation <http://docs.python.org>`_.

.. _dependencies:

Dependencies
~~~~~~~~~~~~

IPython relies on a number of other Python packages. Installing using a package
manager like pip or conda will ensure the necessary packages are installed.
Manual installation without dependencies is possible, but not recommended.
The dependencies can be viewed with package manager commands,
such as :command:`pip show ipython` or :command:`conda info ipython`.


Installing IPython itself
~~~~~~~~~~~~~~~~~~~~~~~~~

IPython requires several dependencies to work correctly, it is not recommended
to install IPython and all its dependencies manually as this can be quite long
and troublesome. You should use the python package manager ``pip``.

Installation using pip
~~~~~~~~~~~~~~~~~~~~~~

Make sure you have the latest version of :mod:`pip` (the Python package
manager) installed. If you do not, head to `Pip documentation
<https://pip.pypa.io/en/stable/installing/>`_ and install :mod:`pip` first.

The quickest way to get up and running with IPython is to install it with pip:

.. code-block:: bash

    $ pip install ipython

That's it.


Installation from source
~~~~~~~~~~~~~~~~~~~~~~~~

To install IPython from source,
grab the latest stable tarball of IPython `from PyPI
<https://pypi.python.org/pypi/ipython>`__.  Then do the following:

.. code-block:: bash

    tar -xzf ipython-5.1.0.tar.gz
    cd ipython-5.1.0
    # The [test] extra ensures test dependencies are installed too:
    pip install '.[test]'

Do not invoke ``setup.py`` directly as this can have undesirable consequences
for further upgrades. We do not recommend using ``easy_install`` either.

If you are installing to a location (like ``/usr/local``) that requires higher
permissions, you may need to run the last command with :command:`sudo`. You can
also install in user specific location by using the ``--user`` flag in
conjunction with pip.

To run IPython's test suite, use the :command:`pytest` command:

.. code-block:: bash

    $ pytest

.. _devinstall:

Installing the development version
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is also possible to install the development version of IPython from our
`Git <http://git-scm.com/>`_ source code repository.  To do this you will
need to have Git installed on your system.  


Then do:

.. code-block:: bash

    $ git clone https://github.com/ipython/ipython.git
    $ cd ipython
    $ pip install -e '.[test]'

The :command:`pip install -e .` command allows users and developers to follow
the development branch as it changes by creating links in the right places and
installing the command line scripts to the appropriate locations.

Then, if you want to update your IPython at any time, do:

.. code-block:: bash

    $ git pull

If the dependencies or entrypoints have changed, you may have to run

.. code-block:: bash

    $ pip install -e .

again, but this is infrequent.
.. _api-index:

###################
  The IPython API
###################

.. only:: html

   :Release: |version|
   :Date: |today|

.. include:: generated/gen.txt
