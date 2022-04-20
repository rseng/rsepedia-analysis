# Release 0.9.8
- Fixed an issue with unused variables appearing in common terms in jacobians.
- Fixes a bug printing large integers, in some environments: Large numbers such as 8.034023767017109e+27 whch were actually ints would be printed as an the int 8034023767017108950029959168 and then the C++ compiler would complain as that's more than the maximum int value. This has been fixed by only printing ints as ints if they are `MIN_INT < number < MAX_INT` else we print them as a float (and we get the sceintific notation).
- Fixed tests to pass with sympy 1.10 and required latest cellmlmanip, which also workes with sympy1.10. Updated sympy requirement to be >=1.9, < 1.11

# Release 0.9.5
- Corrected a type in the generated output for `--rush-larsen-c`
- Updated the help text to no longer imply that CellML files can be converted from remote URIs (a local file is required)
- Fixed required sympy version to be < 1.9 since the jacobean generation has in sympy1.10 changed and makes chaste_codegen tests fail
 
# Release 0.9.4
- This version drops python 3.5 support. The reason for this is python 3.5 is end of life and the chase project will soon be dropping support.
- Corrected Rush-Larsen output for `--rush-larsen-c`

# Release 0.9.3
- Performance upgrade for `--rush-larsen` using caching on linearity checking.

# Release 0.9.2
- Corrected a typo in command line argument `--skip-ingularity-fixes` renaming it to `--skip-singularity-fixes`

# Release 0.9.1
- Added RushLarsen translators for allowing output as c code (`--rush-larsen-c`) and labview (`--rush-larsen-labview`), These are for people wanting to generate generic c or labview code when using chaste_codegen as a stand-alone utility and are not used by chaste

# Release 0.9.0
- Updated to the latest version of the Web Lab Ontology.
- Better error messages when working with transitively finding terms, and encountering subjects that don't point to an existing node in the model.

# Release 0.8.0
- Added `--use-model-factory` flag to allow code to be added for models to self-register to the model factory of the ApPredict peoject.
- Renamed backward euler cells to be in line with others using `BackwardEulerOpt` / `BackwardEuler` instead of `BackwardEuler` and `BackwardEulerNoLot`.

# Release 0.7.0
- The singularities fixing code from release 0.6.3 has moved to the latest cellmlmanip release (0.3.0). You may see small differences in generated code, due to singularities now being fixed prior to any unit conversion taking place. These should not cause any differences in results when run with chaste.
- Improved capacitance unit checking.

# Release 0.6.3
- Fixed issue with backward euler opt models where unused state variables appeard in ComputeResidual
- Fixed memory leek on generated models within lookup table interpolation on nan/inf values in singularities.
- Removed chaste warnings from backward euler template for cases where inf/nan would end up in the table. Intead an error is thrown when multiple (more than 2) such warnings happen.

- Added an automatic fix for removable singularities in GHK equations (which can be switched off with --skip-ingularity-fixes).
  The process looks for equations of any of the following forms, where U is a function of V:
  - `U / (exp(U) - 1.0)`
  - `U / (1.0 - exp(U))`
  - `(exp(U) - 1.0) / U`
  - `(1.0 - exp(U)) / U`  
  It replaces these with a piecewise 1e-7 either side of U==0 drawing a stright line in the region.
  For example `(V + 5)/(exp(V + 5) - 1)` becomes `((fabs(-V - 5.0000000000000000) < fabs(-4.9999999000000000 / 2 - -5.0000001000000000 / 2)) ? ((-5.0000001000000000 + 5.0) / (-1.0 + exp(-5.0000001000000000 + 5.0)) + (--5.0000001000000000 + V) * ((-4.9999999000000000 + 5.0) / (-1.0 + exp(-4.9999999000000000 + 5.0)) - (-5.0000001000000000 + 5.0) / (-1.0 + exp(-5.0000001000000000 + 5.0))) / (--5.0000001000000000 - 4.9999999000000000)) : ((5.0 + V) / (-1.0 + exp(5.0 + V))))`

  See for more details appendix B in: Johnstone, R. H. (2018). Uncertainty characterisation in action potential modelling for cardiac drug safety [PhD thesis]. University of Oxford. https://ora.ox.ac.uk/objects/uuid:0a28829c-828d-4641-bfb0-11193ef47195
- For lookup tables prevented expressions of the form 1 / A from being added, instead adding A. 1 / A type expressions were causing issues, when A is 0. While many cases have a piecewise to prevent hissing this case, lookup table interpolation might cause issues.
- Updated the way BackwardEuler models are calculated, to allow the jacobian to be taken into consideration for lookup tables.
- Fixed a bug with BackwardEuler models where jacobian common term equations (e.g. var_x0) ended up in lookup tables.

# Release 0.5.4
- Fixed sympy deprecation warning when using sympy 1.7 and bumped cellmlmanip recuirement up to ensure sympy 1.7 compatibility
- Improved support for secondary trigonometric functions such as sec and acoth.
- When used with Cellmlmanip version 0.2.2+ an improved printing of devisions is used. For example `1 / (1/cos(x))` now gets rendered as `cos(x)` whereas previously it would be `1 / 1 / cos(x)` giving the incorrect result. An side effcet of the change is that powers of formulas get extra brackets e.g. `pow((1 / x), 2)`.
- chaste_codegen uses placeholder functions for some common maths functions like exp, in order to delay evaluation till the point where the code is written. For programmers using chaste_codegen as a library, there now is a function called `subs_math_func_placeholders` which can be used on any sympy expression to substitute these placeholders.
  E.g.
  ```
  >> expr
  2.0 * exp_(V)
  >> subs_math_func_placeholders(expr)
  2.0 * exp(V)```

# Release 0.5.3
- Added an additional error messages if cellml files can't be loaded and a warning if a lookup table is specified for a tag not present in the model.

# Release 0.5.2
- Error messages have been improved, especially for errors caused by invalid or missing metadata.

# Release 0.5.1
- Corrected usage mentioned in the readme

# Release 0.5.0
- Now implements lookup tables.
- Multiplication equations x * y * y now have 1.0 terms removed in a way that works more generically.

# Release 0.4.1
- Now outputs sqrt(x) instead of pow(x, 0.5).

# Release 0.4.0
- This release explicitly adds versions for dependencies, rather than leaving it up to the cellmlmanip and Jinja2 packages. Versions are semi-strict allowing for minor updates (which should not break compatibility).
- This release contains a bug fix with regards to the stimulus sign, which was calculated incorrectly when the stimulus equation had been changed due to unit conversions.
- The model type displayed in the header of generated files has been fixed and now displays correctly (instead of always displaying 'normal')

# Release 0.3.0
- This release includes the required ontology ttl files in the release itself.
- Removed test data from release package to save space.

# Release 0.2.0
- The command line interface now allows generating multiple model types in one go.
- The command line interface now has a --show-outputs option
- Some command line argument names have changed to more closely match what pycml used to use, command line arguments are not backwards compatible with release 0.1.0
- An issue with generating modifiers that are also parameters has been fixed, by generating the modifier in situ if it doesn't have a defining equation.
- An issue with GetIIonic where the sign of the equations was incorrect has been fixed by analysing the equation before any unit conversion has taken place.

# Release 0.1.0
Initial release of chaste_codegen
## Description
<!--- Describe your changes in detail -->

## Motivation and Context
<!--- Why is this change required? What problem does it solve? -->
<!--- If it fixes an open issue, please link to the issue here, using the 'fixes #<issue>' syntax. -->

## Types of changes
<!--- What types of changes does your code introduce? Put an `x` in all the boxes that apply: -->
- [ ] Bug fix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to change)

## Documentation checklist
<!--- Go over all the following points, and put an `x` in all the boxes that apply. -->
- [ ] I have updated all documentation in the code where necessary.
- [ ] I have checked spelling in all (new) comments and documentation.
- [ ] I have added a note to RELEASE.md if relevant (new feature, breaking change, or notable bug fix).
- [ ] I have updated version & citation.txt & citation.cff version.

## Testing
- [X] Testing is done automatically and codecov shows test coverage
- [ ] This cannot be tested automatically <!-- describe how it has been tested)-->

# Contributing to Chaste codegen

## Installation

Users install `chaste_codegen` using pip.
For instructions on installing as a user, see [README.md](README.md)

Developers should:

1. Clone the repository
2. Create a virtual environment, using e.g. virtualenv or conda. Make sure to use Python3 (e.g. `$ virtualenv venv -p python3`). (If you are on windows you might need to install virtual env first with `pip install virtualenv`. Make sure your python3 installation is in your path.)
3. Activate the environment (e.g. `$ source venv/bin/activate`). (On Windows, virtualenv creates a batch file to activate the virtualenv: `\path\to\env\Scripts\activate`)
4. Install the developer requirements into the virtual environment: `pip install -r dev-requirements/dev.txt`
5. Run the tests: `$ python -m pytest`.

### Requirements

There are two lists of requirements.

1. User requirements, specified in `setup.py`.
2. Developer requirements, specified in the `dev-requirements` directory.

User requirements specify minimum versions of each dependency, and can be used in an existing Python ecosystem.
Developer requirements specify "pinned" versions of each dependency, and should be installed inside a virtual environment.
Having pinned versions ensures different developers get consistent results.
Continuous integration testing happens with user requirements.

Using a virtualenv for development also helps you notice when requirements are missing from the developer requirements.
Similarly, on CI testing anything not listed in the user requirements will cause the tests to fail.

**User requirements** are listed in `setup.py`.
They are divided into base dependencies (listed in `install_requires`) and `test` dependencies (listed in `extras_require[test]`).
Users install these requirements automatically when they `pip install chaste_codegen`.

**Developer requirements** are listed in `base.in`, `test.in`, and `dev.in` (where `test` requires `base`, while `dev` requires `base` and `test`).
These `.in` files are compiled into `.txt` files for `pip` using [pip-tools](https://pypi.org/project/pip-tools/).
To compile them, use the Makefile in the `dev-requirements` folder, by simply typing `make`.

To install the developer requirements into your virtualenv, make sure you've activated the virtual environment, and then:

```sh
$ pip install -r dev-requirements/dev.txt
```


## Coding style guidelines

We follow the [PEP8 recommendations](https://www.python.org/dev/peps/pep-0008/) for coding style, and use [flake8](http://flake8.pycqa.org/en/latest/) to check our PEP8 adherence. To run, type

```sh
$ flake8
```

### Python version

Python 3.5+


## Testing

We're using [pytest](https://docs.pytest.org/en/latest/). To run, type

```sh
$ python -m pytest
```


## Documentation

Every method and every class should have a [docstring](https://www.python.org/dev/peps/pep-0257/) that describes in plain terms what it does, and what the expected input and output is.

These docstrings can be fairly simple, but can also make use of [reStructuredText](http://docutils.sourceforge.net/docs/user/rst/quickref.html), a markup language designed specifically for writing [technical documentation](https://en.wikipedia.org/wiki/ReStructuredText). For example, you can link to other classes and methods by writing ```:class:`myokit.Model` ``` and  ```:meth:`run()` ```.


## Infrastructure & configuration files

### Visual Studio Code Development Environment
Visual Studio Code (not to be confused with Visual Studio) is a cross-platform free lightweight IDE from Microsoft.
To use it for ``chaste_codegen`` development you will need to install the [Python extension](https://marketplace.visualstudio.com/items?itemName=ms-python.python).
To be able to run and debug tests you will also need to install the [pytest extension](https://code.visualstudio.com/docs/python/testing).
In your settings you should also make it use the correct vritual environment. This can be done via the UI but it is a bit hit and miss. 
However there is a settings file in the `.vscode` folder in the ``chaste_codegen`` folder called `settings.json`. Note that `.vscode` is included in our `.gitignore` so any changes will not be tracked or committed.
Below is an example `settings.json`; replace `<python_path>` with the folder that contains `python.exe` for your virtual environment (in Windows it is in `\Scripts` and on linux in `/bin `).
`<project_path>` is just the root folder of the ``chaste_codegen`` code.
```
{
    "python.pythonPath" : "<python_path>",
    "python.venvPath": "${workspaceFolder}",
    "python.venvFolders": [
        "<project_path>"
    ],

    "python.testing.unittestEnabled": false,
    "python.testing.nosetestsEnabled": false,
    "python.testing.pytestEnabled": true,
    "python.testing.pytestArgs": ["--rootdir=.", "--verbose"]
}
```

### Git setup
The tests contain a large amount of reference files. When reference files are updated it's a common practice to regenerate them all (after they have been tested with chaste). Often only a few will have changes.
In order to hide reference files for which only the timestamps have changed, please set up your git environment as follows.
```
git config --global filter.strip_gen_time.clean "sed 's;^//! on .*;//! on (date omitted as unimportant);'"
git config --global filter.strip_gen_time.smudge cat
```
See https://git-scm.com/book/en/v2/Customizing-Git-Git-Attributes#_keyword_expansion for full details.
If you are on windows you may need to install sed for windows and add it to your path: http://gnuwin32.sourceforge.net/packages/sed.htm
![workflow](https://github.com/ModellingWebLab/chaste-codegen/actions/workflows/pytest.yml/badge.svg) [![Documentation Status](https://readthedocs.org/projects/chaste-codegen/badge/?version=latest)](https://chaste-codegen.readthedocs.io/en/latest/?badge=latest) [![codecov](https://codecov.io/gh/ModellingWebLab/chaste-codegen/branch/master/graph/badge.svg)](https://codecov.io/gh/ModellingWebLab/chaste-codegen)

# Code generation for cardiac Chaste

The `chaste_codegen` module takes [CellML](https://www.cellml.org/) models as input, via [cellmlmanip](https://github.com/ModellingWebLab/cellmlmanip) to read and manipulate them, then uses templating to generate C++ code.

The [jinja2](http://jinja.pocoo.org/) templating engine is used.

## Installing 
We recommend installing chaste_codegen in a virtual environment (using virtualenv or conda)

Users install `chaste_codegen` using pip.

```
pip install chaste_codegen
```

## Using `chaste_codegen`
After installation, chaste_codegen can be called using the `chaste_codegen` command:
```
usage: chaste_codegen [-h] [--version] [--normal] [--cvode]
                      [--cvode-data-clamp] [--backward-euler] [--rush-larsen]
                      [--grl1] [--grl2] [--rush-larsen-labview]
                      [--rush-larsen-c] [-j] [-o OUTFILE]
                      [--output-dir OUTPUT_DIR] [--show-outputs] [-c CLS_NAME]
                      [-q] [--skip-ingularity-fixes] [-y] [--opt] [-m]
                      [--lookup-table <metadata tag> min max step]
                      [--use-model-factory]
                      cellml_file
chaste_codegen: error: the following arguments are required: cellml_file

```

For more information about the available options call
`chaste_codegen -h` or see the [CodeGenerationFromCellML guide](https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides/CodeGenerationFromCellML) 


## Release notes
For release notes see [RELEASE.md](./RELEASE.md)


## Documentation
API documentation explaining how to use cellmlmanip can be found on [readthedocs](https://chaste-codegen.readthedocs.io/en/latest/)


## Contributing
For guidelines on contributing to `chaste_codegen`, please see [CONTRIBUTING.md](CONTRIBUTING.md).
.. Root of all chaste_codegen docs

.. _GitHub: https://github.com/ModellingWebLab/chaste-codegen


Welcome to the chaste_codegen documentation
===========================================

chaste_codegen is hosted on GitHub_, where you can find the code and installation instructions.

This page provides the *API*, or *developer documentation* for ``chaste_codegen``.

API documentation
=================
.. automodapi:: chaste_codegen
   :no-inheritance-diagram:
 