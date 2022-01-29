
pypocketfft version
-------------------

SciPy currently vendors [pypocketfft] at:

    commit bf2c431c21213b7c5e23c2f542009b0bd3ec1445
    Merge: 8234f5c 9c252b2
    Author: Martin Reinecke <martin@mpa-garching.mpg.de>
    Date:   Tue Oct 20 15:16:14 2020 +0200

        Merge branch 'good_size_keywords' into 'master'

        Handle keyword args in good_size

        See merge request mtr/pypocketfft!41

pypocketfft: https://gitlab.mpcdf.mpg.de/mtr/pypocketfft
pypocketfft
===========

This package provides Fast Fourier, trigonometric, and Hartley transforms with a
simple Python interface.

The central algorithms are derived from Paul Swarztrauber's FFTPACK code
(http://www.netlib.org/fftpack).

Features
--------
- supports fully complex and half-complex (i.e., complex-to-real and
  real-to-complex) FFTs, discrete sine/cosine transforms, and Hartley transforms
- achieves very high accuracy for all transforms
- supports multidimensional arrays and selection of the axes to be transformed
- supports single, double, and long double precision
- makes use of CPU vector instructions when performing 2-D and higher-dimensional
  transforms
- supports prime-length transforms without degrading to O(N**2) performance
- has optional OpenMP support for multidimensional transforms
# AMOS

A Portable Package for Bessel Functions of a Complex Argument
and Nonnegative Order

This algorithm is a package of subroutines for computing Bessel
functions and Airy functions.  The routines are updated
versions of those routines found in TOMS algorithm 644.

## Disclaimer

```
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                 ISSUED BY SANDIA LABORATORIES,
*                   A PRIME CONTRACTOR TO THE
*               UNITED STATES DEPARTMENT OF ENERGY
* * * * * * * * * * * * * *  NOTICE   * * * * * * * * * * * * * * *
* THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
* UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
* UNITED STATES DEPARTMENT OF ENERGY, NOR ANY OF THEIR
* EMPLOYEES, NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR
* EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
* LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS
* OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
* DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
* PRIVATELY OWNED RIGHTS.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* THIS CODE HAS BEEN APPROVED FOR UNLIMITED RELEASE.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
```
This contains a few hacks to re-use BOOST C++ library
[tests data](https://github.com/boostorg/math/tree/develop/test) for
validating scipy.special functions.

convert.py
----------

This script parses the BOOST data and writes the data into a CVS text file.
For each data file, a subdirectory is created, as one boost data file (.ipp)
can have several test sets.

From the root of this repo
```
# Clone boostorg repo
git clone --depth=1 https://github.com/boostorg/math.git boostmath

# Remove existing data
rm -rf scipy/special/tests/data/boost/*

# Run the coverter script (potentially also update exclude regexes)
python scipy/special/utils/convert.py

# Verify all the new files are used in test_data.py
git diff --stat HEAD | grep "Bin 0"
git diff HEAD -- scipy/special/tests/test_data.py
```

It may be desirable to remove whitespace only changes (see
[#12357](https://github.com/scipy/scipy/pull/12357)) for instructions on that.
SciPy Code of Conduct
======

You can read our Code of Conduct by following [this link](../doc/source/dev/conduct/code_of_conduct.rst). 


Alternatively, you can find it under `scipy/doc/source/dev/conduct/code_of_conduct.rst`. 

<!-- 
Thanks for contributing a pull request! Please ensure that
your PR satisfies the checklist before submitting:
http://scipy.github.io/devdocs/dev/contributor/development_workflow.html#checklist-before-submitting-a-pr

Also, please name and describe your PR as you would write a
commit message:
http://scipy.github.io/devdocs/dev/contributor/development_workflow.html#writing-the-commit-message

Note that we are a team of volunteers; we appreciate your
patience during the review process.

Again, thanks for contributing!
-->

#### Reference issue
<!--Example: Closes gh-WXYZ.-->

#### What does this implement/fix?
<!--Please explain your changes.-->

#### Additional information
<!--Any additional information you think is important.-->
# SciPy Documentation

## How to build the docs

To build the html docs for local development, SciPy itself needs to be built so your
environment needs to be set up for that.  For details on that, see the
[Contributor Guide](http://scipy.github.io/devdocs/dev/contributor/contributor_toc.html#development-environment)).

Also ensure to initialize and update submodules (this pulls in the SciPy Sphinx
theme and `numpydoc`):
```
git submodule update --init
```

Now to build both SciPy itself and the docs, use:
```
python3 runtests.py --doc html
```

Alternatively, if you prefer to build SciPy and the docs separately rather
than use `runtests.py`:
```
python setup.py develop  # in the root of the repo
cd doc && make html-scipyorg
```

In case the SciPy version found by the above command is different from that of the
latest commit in the repo, you will see a message like:
```
installed scipy 5fd20ec1aa != current repo git version '35fd20ec1a'
```

This indicates that you're likely picking up the wrong SciPy install, check
with `python -c "import scipy; print(scipy.__file__)"`.

If the build is successful, you can open it in your browser with `make show`
(which will open `build/html-scipyorg/index.html`).


## Building pdf docs

To build the pdf docs, which requires a LaTeX install and can be more fiddly
to get to work, replace the doc build commands in the section above with:
```
python3 runtests.py --doc latex
```
or:
```
make latex
```

That will use Sphinx to generate the LaTeX sources. To then produce a pdf,
navigate to `doc/build/latex/` and run:
```
make all-pdf
```

That will produce a file `scipy-ref.pdf` in `build/latex/`.


## Building documentation for a release

For building all the documentation artifacts for a release, run:
```
make dist
```

This will build SciPy in-place (to ensure the version is correct), build html
and pdf docs as well as create a zip archive of the html docs that can easily
be redistributed.


## Layout of the docs in this repository

- `source` is where most of the content lives.
  - `dev` contains the contributor and developer guides as well as the governance
    docs and the code of conduct.
  - `tutorial` contains all tutorial content.
- `release` contains the release notes. Note that those normally should not be
  updated as part of a PR; we keep releases notes for the upcoming releases
  on the wiki of the main SciPy repo.