# <a href="https://numpy.org/"><img alt="NumPy" src="/branding/logo/primary/numpylogo.svg" height="60"></a>

<!--[![Azure Pipelines](https://dev.azure.com/numpy/numpy/_apis/build/status/numpy.numpy?branchName=main)](-->
<!--https://dev.azure.com/numpy/numpy/_build/latest?definitionId=1?branchName=main)-->
<!--[![Actions build_test](https://github.com/numpy/numpy/actions/workflows/build_test.yml/badge.svg)](-->
<!--https://github.com/numpy/numpy/actions/workflows/build_test.yml)-->
<!--[![TravisCI](https://app.travis-ci.com/numpy/numpy.svg?branch=main)](-->
<!--https://app.travis-ci.com/numpy/numpy)-->
<!--[![CircleCI](https://img.shields.io/circleci/project/github/numpy/numpy/main.svg?label=CircleCI)](-->
<!--https://circleci.com/gh/numpy/numpy)-->
<!--[![Codecov](https://codecov.io/gh/numpy/numpy/branch/main/graph/badge.svg)](-->
<!--https://codecov.io/gh/numpy/numpy)-->

[![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](
https://numfocus.org)
[![PyPI Downloads](https://img.shields.io/pypi/dm/numpy.svg?label=PyPI%20downloads)](
https://pypi.org/project/numpy/)
[![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/numpy.svg?label=Conda%20downloads)](
https://anaconda.org/conda-forge/numpy)
[![Stack Overflow](https://img.shields.io/badge/stackoverflow-Ask%20questions-blue.svg)](
https://stackoverflow.com/questions/tagged/numpy)
[![Nature Paper](https://img.shields.io/badge/DOI-10.1038%2Fs41592--019--0686--2-blue)](
https://doi.org/10.1038/s41586-020-2649-2)

NumPy is the fundamental package for scientific computing with Python.

- **Website:** https://www.numpy.org
- **Documentation:** https://numpy.org/doc
- **Mailing list:** https://mail.python.org/mailman/listinfo/numpy-discussion
- **Source code:** https://github.com/numpy/numpy
- **Contributing:** https://www.numpy.org/devdocs/dev/index.html
- **Bug reports:** https://github.com/numpy/numpy/issues
- **Report a security vulnerability:** https://tidelift.com/docs/security

It provides:

- a powerful N-dimensional array object
- sophisticated (broadcasting) functions
- tools for integrating C/C++ and Fortran code
- useful linear algebra, Fourier transform, and random number capabilities

Testing:

NumPy requires `pytest` and `hypothesis`.  Tests can then be run after installation with:

    python -c 'import numpy; numpy.test()'

Code of Conduct
----------------------

NumPy is a community-driven open source project developed by a diverse group of
[contributors](https://numpy.org/teams/). The NumPy leadership has made a strong
commitment to creating an open, inclusive, and positive community. Please read the
[NumPy Code of Conduct](https://numpy.org/code-of-conduct/) for guidance on how to interact
with others in a way that makes our community thrive.

Call for Contributions
----------------------

The NumPy project welcomes your expertise and enthusiasm!

Small improvements or fixes are always appreciated; issues labeled as ["good
first issue"](https://github.com/numpy/numpy/labels/good%20first%20issue)
may be a good starting point. If you are considering larger contributions
to the source code, please contact us through the [mailing
list](https://mail.python.org/mailman/listinfo/numpy-discussion) first.

Writing code isn’t the only way to contribute to NumPy. You can also:
- review pull requests
- help us stay on top of new and old issues
- develop tutorials, presentations, and other educational materials
- maintain and improve [our website](https://github.com/numpy/numpy.org)
- develop graphic design for our brand assets and promotional materials
- translate website content
- help with outreach and onboard new contributors
- write grant proposals and help with other fundraising efforts

For more information about the ways you can contribute to NumPy, visit [our website](https://numpy.org/contribute/). 
If you’re unsure where to start or how your skills fit in, reach out! You can
ask on the mailing list or here, on GitHub, by opening a new issue or leaving a
comment on a relevant issue that is already open.

Our preferred channels of communication are all public, but if you’d like to
speak to us in private first, contact our community coordinators at
numpy-team@googlegroups.com or on Slack (write numpy-team@googlegroups.com for
an invitation).

We also have a biweekly community call, details of which are announced on the
mailing list. You are very welcome to join.

If you are new to contributing to open source, [this
guide](https://opensource.guide/how-to-contribute/) helps explain why, what,
and how to successfully get involved.
PocketFFT
---------

This is a heavily modified implementation of FFTPack [1,2], with the following
advantages:

- strictly C99 compliant
- more accurate twiddle factor computation
- very fast plan generation
- worst case complexity for transform sizes with large prime factors is
  `N*log(N)`, because Bluestein's algorithm [3] is used for these cases.


Some code details
-----------------

Twiddle factor computation:

- making use of symmetries to reduce number of sin/cos evaluations
- all angles are reduced to the range `[0; pi/4]` for higher accuracy
- an adapted implementation of `sincospi()` is used, which actually computes
  `sin(x)` and `(cos(x)-1)`.
- if `n` sin/cos pairs are required, the adjusted `sincospi()` is only called
  `2*sqrt(n)` times; the remaining values are obtained by evaluating the
  angle addition theorems in a numerically accurate way.

Parallel invocation:

- Plans only contain read-only data; all temporary arrays are allocated and
  deallocated during an individual FFT execution. This means that a single plan
  can be used in several threads at the same time.

Efficient codelets are available for the factors:

- 2, 3, 4, 5, 7, 11 for complex-valued FFTs
- 2, 3, 4, 5 for real-valued FFTs

Larger prime factors are handled by somewhat less efficient, generic routines.

For lengths with very large prime factors, Bluestein's algorithm is used, and
instead of an FFT of length `n`, a convolution of length `n2 >= 2*n-1`
is performed, where `n2` is chosen to be highly composite.


[1] Swarztrauber, P. 1982, Vectorizing the Fast Fourier Transforms
    (New York: Academic Press), 51
[2] https://www.netlib.org/fftpack/
[3] https://en.wikipedia.org/wiki/Chirp_Z-transform
NumPy has a Code of Conduct, please see: https://numpy.org/code-of-conduct
# Contributing to numpy

## Reporting issues

When reporting issues please include as much detail as possible about your
operating system, numpy version and python version. Whenever possible, please
also include a brief, self-contained code example that demonstrates the problem.

If you are reporting a segfault please include a GDB traceback, which you can
generate by following
[these instructions.](https://github.com/numpy/numpy/blob/main/doc/source/dev/development_environment.rst#debugging)

## Contributing code

Thanks for your interest in contributing code to numpy!

+ If this is your first time contributing to a project on GitHub, please read
through our
[guide to contributing to numpy](https://numpy.org/devdocs/dev/index.html)
+ If you have contributed to other projects on GitHub you can go straight to our
[development workflow](https://numpy.org/devdocs/dev/development_workflow.html)

Either way, please be sure to follow our
[convention for commit messages](https://numpy.org/devdocs/dev/development_workflow.html#writing-the-commit-message).

If you are writing new C code, please follow the style described in
``doc/C_STYLE_GUIDE``.

Suggested ways to work on your development version (compile and run
the tests without interfering with system packages) are described in
``doc/source/dev/development_environment.rst``.

### A note on feature enhancements/API changes

If you are interested in adding a new feature to NumPy, consider
submitting your feature proposal to the [mailing list][mail], 
which is the preferred forum for discussing new features and
API changes.

[mail]: https://mail.python.org/mailman/listinfo/numpy-discussion
<!--         ----------------------------------------------------------------
                MAKE SURE YOUR PR GETS THE ATTENTION IT DESERVES!
                ----------------------------------------------------------------

*  FORMAT IT RIGHT:
      http://www.numpy.org/devdocs/dev/development_workflow.html#writing-the-commit-message

*  IF IT'S A NEW FEATURE OR API CHANGE, TEST THE WATERS:
      http://www.numpy.org/devdocs/dev/development_workflow.html#get-the-mailing-list-s-opinion

*  HIT ALL THE GUIDELINES:
      https://numpy.org/devdocs/dev/index.html#guidelines

*  WHAT TO DO IF WE HAVEN'T GOTTEN BACK TO YOU:
      http://www.numpy.org/devdocs/dev/development_workflow.html#getting-your-pr-reviewed
-->
# NumPy Logo Guidelines
These guidelines are meant to help keep the NumPy logo consistent and recognizable across all its uses. They also provide a common language for referring to the logos and their components.

The primary logo is the horizontal option (logomark and text next to each other) and the secondary logo is the stacked version (logomark over text). I’ve also provided the logomark on its own (meaning it doesn’t have text). When in doubt, it’s preferable to use primary or secondary options over the logomark alone.

## Color
The full color options are a combo of two shades of blue, rgb(77, 171, 207) and rgb(77, 119, 207), while light options are rgb(255, 255, 255) and dark options are rgb(1, 50, 67).

Whenever possible, use the full color logos. One color logos (light or dark) are to be used when full color will not have enough contrast, usually when logos must be on colored backgrounds.

## Minimum Size
Please do not make the primary logo smaller than 50px wide, secondary logo smaller than 35px wide, or logomark smaller than 20px wide.

## Logo Integrity
A few other notes to keep in mind when using the logo:
- Make sure to scale the logo proportionally.
- Maintain a good amount of space around the logo. Don’t let it overlap with text, images, or other elements.
- Do not try and recreate or modify the logo. For example, do not use the logomark and then try to write NumPy in another font.
Note that since Python 3.6 the builtin tracemalloc module can be used to
track allocations inside numpy.
Numpy places its CPU memory allocations into the `np.lib.tracemalloc_domain`
domain.
See https://docs.python.org/3/library/tracemalloc.html.

The tool that used to be here has been deprecated.
