# Astropy Project Governance

Please visit our website to learn more about the [Astropy Team](http://www.astropy.org/team.html).
All Astropy community members are expected to abide by the [Astropy Project Code of Conduct](http://www.astropy.org/code_of_conduct.html).
Contributing to Astropy
=======================

Reporting Issues
----------------

When opening an issue to report a problem, please try to provide a minimal code
example that reproduces the issue along with details of the operating
system and the Python, NumPy, and `astropy` versions you are using.

Contributing Code and Documentation
-----------------------------------

So you are interested in contributing to the Astropy Project?  Excellent!
We love contributions! Astropy is open source, built on open source, and
we'd love to have you hang out in our community.

Anti Imposter Syndrome Reassurance
----------------------------------

We want your help. No, really.

There may be a little voice inside your head that is telling you that you're not
ready to be an open source contributor; that your skills aren't nearly good
enough to contribute. What could you possibly offer a project like this one?

We assure you - the little voice in your head is wrong. If you can write code or
documentation, you can contribute code to open source.
Contributing to open source projects is a fantastic way to advance one's coding
and open source workflow skills. Writing perfect code isn't the measure of a good
developer (that would disqualify all of us!); it's trying to create something,
making mistakes, and learning from those mistakes. That's how we all improve,
and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can
help out by writing documentation, tests, or even giving feedback about the
project (and yes - that includes giving feedback about the contribution
process). Some of these contributions may be the most valuable to the project as
a whole, because you're coming to the project with fresh eyes, so you can see
the errors and assumptions that seasoned contributors have glossed over.

Note: This text was originally written by
[Adrienne Lowe](https://github.com/adriennefriend) for a
[PyCon talk](https://www.youtube.com/watch?v=6Uj746j9Heo), and was adapted by
Astropy based on its use in the README file for the
[MetPy project](https://github.com/Unidata/MetPy).

How to Contribute, Best Practices
---------------------------------

Most contributions to Astropy are done via [pull
requests](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests)
from GitHub users' forks of the [astropy
repository](https://github.com/astropy/astropy). If you are new to this
style of development, you will want to read over our [development
workflow](https://docs.astropy.org/en/latest/development/workflow/development_workflow.html).

You may also/instead be interested in contributing to an
[astropy affiliated package](https://www.astropy.org/affiliated/).
Affiliated packages are astronomy-related software packages that are not a part
of the `astropy` core package, but build on it for more specialized applications
and follow the Astropy guidelines for reuse, interoperability, and interfacing.
Each affiliated package has its own developers/maintainers and its own specific
guidelines for contributions, so be sure to read their docs.

Once you open a pull request (which should be opened against the ``main``
branch, not against any of the other branches), please make sure to
include the following:

- **Code**: the code you are adding, which should follow
  our [coding guidelines](https://docs.astropy.org/en/latest/development/codeguide.html) as much as possible.

- **Tests**: these are usually tests to ensure code that previously
  failed now works (regression tests), or tests that cover as much as possible
  of the new functionality to make sure it does not break in the future and
  also returns consistent results on all platforms (since we run these tests on
  many platforms/configurations). For more information about how to write
  tests, see our [testing guidelines](https://docs.astropy.org/en/latest/development/testguide.html).

- **Documentation**: if you are adding new functionality, be sure to include a
  description in the main documentation (in ``docs/``). Again, we have some
  detailed [documentation guidelines](https://docs.astropy.org/en/latest/development/docguide.html) to help you out.

- **Performance improvements**: if you are making changes that impact `astropy`
  performance, consider adding a performance benchmark in the
  [astropy-benchmarks](https://github.com/astropy/astropy-benchmarks)
  repository. You can find out more about how to do this
  [in the README for that repository](https://github.com/astropy/astropy-benchmarks#contributing-benchmarks).

- **Changelog entry**: whether you are fixing a bug or adding new
  functionality, you should add a changelog fragment in the ``docs/changes/``
  directory. See ``docs/changes/README.rst`` for some guidance on the creation
  of this file.

  If you are opening a pull request you may not know
  the PR number yet, but you can add it once the pull request is open. If you
  are not sure where to put the changelog entry, wait until a maintainer
  has reviewed your PR and assigned it to a milestone.

  You do not need to include a changelog entry for fixes to bugs introduced in
  the developer version and therefore are not present in the stable releases. In
  general you do not need to include a changelog entry for minor documentation
  or test updates. Only user-visible changes (new features/API changes, fixed
  issues) need to be mentioned. If in doubt, ask the core maintainer reviewing
  your changes.

Checklist for Contributed Code
------------------------------

Before being merged, a pull request for a new feature will be reviewed to see if
it meets the following requirements. If you are unsure about how to meet all of these
requirements, please submit the PR and ask for help and/or guidance. An Astropy
maintainer will collaborate with you to make sure that the pull request meets the
requirements for inclusion in the package:

**Scientific Quality** (when applicable)
  * Is the submission relevant to astronomy?
  * Are references included to the origin source for the algorithm?
  * Does the code perform as expected?
  * Has the code been tested against previously existing implementations?

**Code Quality**
  * Are the [coding guidelines](https://docs.astropy.org/en/latest/development/codeguide.html) followed?
  * Is the code compatible with Python >=3.8?
  * Are there dependencies other than the `astropy` core, the Python Standard
    Library, and NumPy 1.18.0 or later?
    * Is the package importable even if the C-extensions are not built?
    * Are additional dependencies handled appropriately?
    * Do functions that require additional dependencies raise an `ImportError`
      if they are not present?

**Testing**
  * Are the [testing guidelines](https://docs.astropy.org/en/latest/development/testguide.html) followed?
  * Are the inputs to the functions sufficiently tested?
  * Are there tests for any exceptions raised?
  * Are there tests for the expected performance?
  * Are the sources for the tests documented?
  * Have tests that require an [optional dependency](https://docs.astropy.org/en/latest/development/testguide.html#tests-requiring-optional-dependencies)
    been marked as such?
  * Does ``tox -e test`` run without failures?

**Documentation**
  * Are the [documentation guidelines](https://docs.astropy.org/en/latest/development/docguide.html) followed?
  * Is there a docstring in [numpydoc format](https://numpydoc.readthedocs.io/en/latest/format.html) in the function describing:
    * What the code does?
    * The format of the inputs of the function?
    * The format of the outputs of the function?
    * References to the original algorithms?
    * Any exceptions which are raised?
    * An example of running the code?
  * Is there any information needed to be added to the docs to describe the
    function?
  * Does the documentation build without errors or warnings?

**License**
  * Is the `astropy` license included at the top of the file?
  * Are there any conflicts with this code and existing codes?

**Astropy requirements**
  * Do all the GitHub Actions and CircleCI tests pass? If not, are they allowed to fail?
  * If applicable, has an entry been added into the changelog?
  * Can you check out the pull request and repeat the examples and tests?

Other Tips
----------

- Behind the scenes, we conduct a number of tests or checks with new pull requests.
  This is a technique that is called continuous integration, and we use GitHub Actions
  and CircleCI. To prevent the automated tests from running, you can add ``[ci skip]``
  to your commit message. This is useful if your PR is a work in progress (WIP) and
  you are not yet ready for the tests to run. For example:

      $ git commit -m "WIP widget [ci skip]"

  - If you already made the commit without including this string, you can edit
    your existing commit message by running:

        $ git commit --amend

- If your commit makes substantial changes to the documentation but none of
  those changes include code snippets, then you can use ``[ci skip]``,
  which will skip all CI except RTD, where the documentation is built.

- When contributing trivial documentation fixes (i.e., fixes to typos, spelling,
  grammar) that don't contain any special markup and are not associated with
  code changes, please include the string ``[ci skip]`` in your commit
  message.

      $ git commit -m "Fixed typo [ci skip]"

- ``[ci skip]`` and ``[skip ci]`` are the same and can be used interchangeably.
<!-- This comments are hidden when you submit the pull request,
so you do not need to remove them! -->

<!-- Please be sure to check out our contributing guidelines,
https://github.com/astropy/astropy/blob/main/CONTRIBUTING.md .
Please be sure to check out our code of conduct,
https://github.com/astropy/astropy/blob/main/CODE_OF_CONDUCT.md . -->

<!-- If you are new or need to be re-acquainted with Astropy
contributing workflow, please see
http://docs.astropy.org/en/latest/development/workflow/development_workflow.html .
There is even a practical example at
https://docs.astropy.org/en/latest/development/workflow/git_edit_workflow_examples.html#astropy-fix-example . -->

<!-- Astropy coding style guidelines can be found here:
https://docs.astropy.org/en/latest/development/codeguide.html#coding-style-conventions
Our testing infrastructure enforces to follow a subset of the PEP8 to be
followed. You can check locally whether your changes have followed these by
running the following command:

tox -e codestyle

-->

<!-- Please just have a quick search on GitHub to see if a similar
pull request has already been posted.
We have old closed pull requests that might provide useful code or ideas
that directly tie in with your pull request. -->

<!-- We have several automatic features that run when a pull request is open.
They can appear daunting but do not worry because maintainers will help
you navigate them, if necessary. -->

### Description
<!-- Provide a general description of what your pull request does.
Complete the following sentence and add relevant details as you see fit. -->

<!-- In addition please ensure that the pull request title is descriptive
and allows maintainers to infer the applicable subpackage(s). -->

<!-- READ THIS FOR MANUAL BACKPORT FROM A MAINTAINER:
Apply "skip-basebranch-check" label **before** you open the PR! -->

This pull request is to address ...

<!-- If the pull request closes any open issues you can add this.
If you replace <Issue Number> with a number, GitHub will automatically link it.
If this pull request is unrelated to any issues, please remove
the following line. -->

Fixes #<Issue Number>

### Checklist for package maintainer(s)
<!-- This section is to be filled by package maintainer(s) who will
review this pull request. -->

This checklist is meant to remind the package maintainer(s) who will review this pull request of some common things to look for. This list is not exhaustive.

- [ ] Do the proposed changes actually accomplish desired goals?
- [ ] Do the proposed changes follow the [Astropy coding guidelines](https://docs.astropy.org/en/latest/development/codeguide.html)?
- [ ] Are tests added/updated as required? If so, do they follow the [Astropy testing guidelines](https://docs.astropy.org/en/latest/development/testguide.html)?
- [ ] Are docs added/updated as required? If so, do they follow the [Astropy documentation guidelines](https://docs.astropy.org/en/latest/development/docguide.html#astropy-documentation-rules-and-guidelines)?
- [ ] Is rebase and/or squash necessary? If so, please provide the author with appropriate instructions. Also see ["When to rebase and squash commits"](https://docs.astropy.org/en/latest/development/when_to_rebase.html).
- [ ] Did the CI pass? If no, are the failures related? If you need to run daily and weekly cron jobs as part of the PR, please apply the `Extra CI` label.
- [ ] Is a change log needed? If yes, did the change log check pass? If no, add the `no-changelog-entry-needed` label. If this is a manual backport, use the `skip-changelog-checks` label unless special changelog handling is necessary.
- [ ] Is a milestone set? Milestone must be set but `astropy-bot` check might be missing; do not let the green checkmark fool you.
- [ ] At the time of adding the milestone, if the milestone set requires a backport to release branch(es), apply the appropriate `backport-X.Y.x` label(s) *before* merge.
---
name: Bug report
about: Create a report describing unexpected or incorrect behavior in astropy.
labels: Bug
---

<!-- This comments are hidden when you submit the issue,
so you do not need to remove them! -->

<!-- Please be sure to check out our contributing guidelines,
https://github.com/astropy/astropy/blob/main/CONTRIBUTING.md .
Please be sure to check out our code of conduct,
https://github.com/astropy/astropy/blob/main/CODE_OF_CONDUCT.md . -->

<!-- Please have a search on our GitHub repository to see if a similar
issue has already been posted.
If a similar issue is closed, have a quick look to see if you are satisfied
by the resolution.
If not please go ahead and open an issue! -->

<!-- Please check that the development version still produces the same bug.
You can install development version with
pip install git+https://github.com/astropy/astropy
command. -->

### Description
<!-- Provide a general description of the bug. -->

### Expected behavior
<!-- What did you expect to happen. -->

### Actual behavior
<!-- What actually happened. -->
<!-- Was the output confusing or poorly described? -->

### Steps to Reproduce
<!-- Ideally a code example could be provided so we can run it ourselves. -->
<!-- If you are pasting code, use triple backticks (```) around
your code snippet. -->
<!-- If necessary, sanitize your screen output to be pasted so you do not
reveal secrets like tokens and passwords. -->

1. [First Step]
2. [Second Step]
3. [and so on...]

```python
# Put your Python code snippet here.
```

### System Details
<!-- Even if you do not think this is necessary, it is useful information for the maintainers.
Please run the following snippet and paste the output below:
import platform; print(platform.platform())
import sys; print("Python", sys.version)
import numpy; print("Numpy", numpy.__version__)
import erfa; print("pyerfa", erfa.__version__)
import astropy; print("astropy", astropy.__version__)
import scipy; print("Scipy", scipy.__version__)
import matplotlib; print("Matplotlib", matplotlib.__version__)
-->
---
name: Feature request
about: Suggest an idea to improve astropy
labels: "Feature Request"
---

<!-- This comments are hidden when you submit the issue,
so you do not need to remove them! -->

<!-- Please be sure to check out our contributing guidelines,
https://github.com/astropy/astropy/blob/main/CONTRIBUTING.md .
Please be sure to check out our code of conduct,
https://github.com/astropy/astropy/blob/main/CODE_OF_CONDUCT.md . -->

<!-- Please have a search on our GitHub repository to see if a similar
issue has already been posted.
If a similar issue is closed, have a quick look to see if you are satisfied
by the resolution.
If not please go ahead and open an issue! -->

### Description
<!-- Provide a general description of the feature you would like. -->
<!-- If you want to, you can suggest a draft design or API. -->
<!-- This way we have a deeper discussion on the feature. -->


### Additional context
<!-- Add any other context or screenshots about the feature request here. -->
<!-- This part is optional. -->
[![Travis CI Build Status](https://travis-ci.org/libexpat/libexpat.svg?branch=master)](https://travis-ci.org/libexpat/libexpat)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/libexpat/libexpat?svg=true)](https://ci.appveyor.com/project/libexpat/libexpat)
[![Packaging status](https://repology.org/badge/tiny-repos/expat.svg)](https://repology.org/metapackage/expat/versions)


# Expat, Release 2.2.9

This is Expat, a C library for parsing XML, started by
[James Clark](https://en.wikipedia.org/wiki/James_Clark_(programmer)) in 1997.
Expat is a stream-oriented XML parser.  This means that you register
handlers with the parser before starting the parse.  These handlers
are called when the parser discovers the associated structures in the
document being parsed.  A start tag is an example of the kind of
structures for which you may register handlers.

Expat supports the following compilers:
- GNU GCC >=4.5
- LLVM Clang >=3.5
- Microsoft Visual Studio >=8.0/2005

Windows users should use the
[`expat_win32` package](https://sourceforge.io/projects/expat/files/expat_win32/),
which includes both precompiled libraries and executables, and source code for
developers.

Expat is [free software](https://www.gnu.org/philosophy/free-sw.en.html).
You may copy, distribute, and modify it under the terms of the License
contained in the file
[`COPYING`](https://github.com/libexpat/libexpat/blob/master/expat/COPYING)
distributed with this package.
This license is the same as the MIT/X Consortium license.

If you are building Expat from a check-out from the
[Git repository](https://github.com/libexpat/libexpat/),
you need to run a script that generates the configure script using the
GNU autoconf and libtool tools.  To do this, you need to have
autoconf 2.58 or newer. Run the script like this:

```console
./buildconf.sh
```

Once this has been done, follow the same instructions as for building
from a source distribution.

To build Expat from a source distribution, you first run the
configuration shell script in the top level distribution directory:

```console
./configure
```

There are many options which you may provide to configure (which you
can discover by running configure with the `--help` option).  But the
one of most interest is the one that sets the installation directory.
By default, the configure script will set things up to install
libexpat into `/usr/local/lib`, `expat.h` into `/usr/local/include`, and
`xmlwf` into `/usr/local/bin`.  If, for example, you'd prefer to install
into `/home/me/mystuff/lib`, `/home/me/mystuff/include`, and
`/home/me/mystuff/bin`, you can tell `configure` about that with:

```console
./configure --prefix=/home/me/mystuff
```

Another interesting option is to enable 64-bit integer support for
line and column numbers and the over-all byte index:

```console
./configure CPPFLAGS=-DXML_LARGE_SIZE
```

However, such a modification would be a breaking change to the ABI
and is therefore not recommended for general use &mdash; e.g. as part of
a Linux distribution &mdash; but rather for builds with special requirements.

After running the configure script, the `make` command will build
things and `make install` will install things into their proper
location.  Have a look at the `Makefile` to learn about additional
`make` options.  Note that you need to have write permission into
the directories into which things will be installed.

If you are interested in building Expat to provide document
information in UTF-16 encoding rather than the default UTF-8, follow
these instructions (after having run `make distclean`).
Please note that we configure with `--without-xmlwf` as xmlwf does not
support this mode of compilation (yet):

1. Mass-patch `Makefile.am` files to use `libexpatw.la` for a library name:
   <br/>
   `find -name Makefile.am -exec sed
       -e 's,libexpat\.la,libexpatw.la,'
       -e 's,libexpat_la,libexpatw_la,'
       -i {} +`

1. Run `automake` to re-write `Makefile.in` files:<br/>
   `automake`

1. For UTF-16 output as unsigned short (and version/error strings as char),
   run:<br/>
   `./configure CPPFLAGS=-DXML_UNICODE --without-xmlwf`<br/>
   For UTF-16 output as `wchar_t` (incl. version/error strings), run:<br/>
   `./configure CFLAGS="-g -O2 -fshort-wchar" CPPFLAGS=-DXML_UNICODE_WCHAR_T
       --without-xmlwf`
   <br/>Note: The latter requires libc compiled with `-fshort-wchar`, as well.

1. Run `make` (which excludes xmlwf).

1. Run `make install` (again, excludes xmlwf).

Using `DESTDIR` is supported.  It works as follows:

```console
make install DESTDIR=/path/to/image
```

overrides the in-makefile set `DESTDIR`, because variable-setting priority is

1. commandline
1. in-makefile
1. environment

Note: This only applies to the Expat library itself, building UTF-16 versions
of xmlwf and the tests is currently not supported.

When using Expat with a project using autoconf for configuration, you
can use the probing macro in `conftools/expat.m4` to determine how to
include Expat.  See the comments at the top of that file for more
information.

A reference manual is available in the file `doc/reference.html` in this
distribution.


The CMake build system is still *experimental* and will replace the primary
build system based on GNU Autotools at some point when it is ready.
For an idea of the available (non-advanced) options for building with CMake:

```console
# rm -f CMakeCache.txt ; cmake -D_EXPAT_HELP=ON -LH . | grep -B1 ':.*=' | sed 's,^--$,,'
// Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel ...
CMAKE_BUILD_TYPE:STRING=

// Install path prefix, prepended onto install directories.
CMAKE_INSTALL_PREFIX:PATH=/usr/local

// Path to a program.
DOCBOOK_TO_MAN:FILEPATH=/usr/bin/docbook2x-man

// build man page for xmlwf
EXPAT_BUILD_DOCS:BOOL=ON

// build the examples for expat library
EXPAT_BUILD_EXAMPLES:BOOL=ON

// build fuzzers for the expat library
EXPAT_BUILD_FUZZERS:BOOL=OFF

// build the tests for expat library
EXPAT_BUILD_TESTS:BOOL=ON

// build the xmlwf tool for expat library
EXPAT_BUILD_TOOLS:BOOL=ON

// Character type to use (char|ushort|wchar_t) [default=char]
EXPAT_CHAR_TYPE:STRING=char

// install expat files in cmake install target
EXPAT_ENABLE_INSTALL:BOOL=ON

// Use /MT flag (static CRT) when compiling in MSVC
EXPAT_MSVC_STATIC_CRT:BOOL=OFF

// build a shared expat library
EXPAT_SHARED_LIBS:BOOL=ON

// Treat all compiler warnings as errors
EXPAT_WARNINGS_AS_ERRORS:BOOL=OFF

// Make use of getrandom function (ON|OFF|AUTO) [default=AUTO]
EXPAT_WITH_GETRANDOM:STRING=AUTO

// utilize libbsd (for arc4random_buf)
EXPAT_WITH_LIBBSD:BOOL=OFF

// Make use of syscall SYS_getrandom (ON|OFF|AUTO) [default=AUTO]
EXPAT_WITH_SYS_GETRANDOM:STRING=AUTO
```
