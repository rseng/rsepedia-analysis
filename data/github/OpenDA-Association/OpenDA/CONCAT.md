
# Contributing to OpenDA

This document briefly describes the policy for contributing changes to OpenDA.
First we would like to emphasize that anyone can create a local copy and make your own
local modifications of OpenDA using the normal commands of github and git. Only when you want to
share your code with the OpenDA community this document applies.

## Share your code with the OpenDA community

For OpenDA we distinguish two methods for contributing:
1. Create a fork on github and then a pull request. When you contribute for the first time, then this is probably what you want.
2. Ask a core-developer of OpenDA to create a branch where you can push you contributions. If you do not know who the core-developers of OpenDA are then you can contact one at info@openda.org

## Repository structure

As a general rule all code changes are first committed to a branch. These branches are usually personal, so you can work here in your without interference of other developments. There are a few special branches, with the name master or starting with release. These names are reserved and should not be used for development.

## Step by step

### Fork scenario
1. Create a fork on github
2. Clone your newly created fork `git clone https://github.com/<your account>/OpenDA.git`
3. Create a branch with a name that describes your intended developments. `git checkout -b my_branch`
4. Write and test your code in your branch.
5. Commit to your local copy and push to your fork of OpenDA
   * `git add bla`
   * `git commit -m "my description"`
   * `git push origin my_branch`
6. Once finished or more regularly: merge changes in OpenDA into your work
  * Update your fork of OpenDA on github
  * `git fetch origin master` to update your local copy (from your fork)
  * `git merge master`
7. Push changes to the server again `git push origin my_branch` (to your fork)
8. Create a pull request on github

### Directly on main repository scenario
1. Clone your OpenDA repository `git clone https://github.com/OpenDA-Association/OpenDA.git`
3. Create a branch with a name that describes your intended developments.  `git checkout -b my_branch`
4. Write and test your code in your branch.
5. Commit to local copy and push to your fork of OpenDA
  * `git add bla`
  * `git commit -m "my description"`
  * `git push origin my_branch`
6. Once finished or more regularly: merge changes in OpenDA into your work
  * `git fetch origin master` to update your local copy
  * `git merge master`
7. Push changes to OpenDA on github again `git pus origin my_branch`
8. Create a pull request on github

## Questions

Any questions? You can post them on the [Forum](https://sourceforge.net/p/openda/discussion/?source=navbar) or mail us at info@openda.org
[![Java CI](https://github.com/OpenDA-Association/OpenDA/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/OpenDA-Association/OpenDA/actions/workflows/ci.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/465833e082b54b279105a280b36c75b8)](https://www.codacy.com/gh/OpenDA-Association/OpenDA/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=OpenDA-Association/OpenDA&amp;utm_campaign=Badge_Grade)

# OpenDA

OpenDA is an open interface standard for (and free implementation of) a set of tools to quickly implement data-assimilation and calibration for arbitrary numerical models. OpenDA wants to stimulate the use of data-assimilation and calibration by lowering the implementation costs and enhancing the exchange of software among researchers and end-users.
A model that conforms to the OpenDA standard can use all the tools that are available in OpenDA. This allows experimentation with data-assimilation/calibration methods without the need for extensive programming. Reversely, developers of data-assimilation/calibration software that make their implementations compatible with the OpenDA interface will make their new methods usable for all OpenDA users (either for free or on a commercial basis).
OpenDA has been designed for high performance. Hence, even large-scale models can use it. Also, OpenDA allows users to optimize the interaction between their model and the data-assimilation/calibration methods. Hence, data-assimilation with OpenDA can be as efficient as with custom-made implementations of data-assimilation methods.
OpenDA is an Open Source project. Contributions are welcome from anyone wishing to participate in the further development of the OpenDA toolset.

## Features of OpenDA

Data-assimilation methods

- Ensemble KF (EnKF)
- Ensemble SquareRoot KF (EnSR)
- Steady State KF
- Particle Filter
- 3DVar
- DudEnKF (still under research)
- DudEnSR (still under research)

Parameter estimation (calibration) methods:

- Dud
- Sparse Dud
- Simplex
- Powell
- Gridded full search
- Shuffled Comples Evolution (SCE)
- Generalized Likelihood Uncertainty Estimation (GLUE)
- (L)BFGS
- Conjugate Gradient: Fleetjer-Reeves, Polak-Ribiere, Steepest Descent
- Uncertainty Analaysis methods
- GLUE
- DELSA

Language interfaces

- C/C++
- Java
- Fortran77/90

These files are part of the OpenDA software. For more information see our website at
http://www.openda.org


This module is essentially empty. It contains no useful code or data, but serves as
a template for creating new modules. One can create a new module by:

1. make a copy of the template 

	> svn cp module_template my_module

2. change module name in several files:

    - `module.properties`
	- `java/unit_test_info.txt`
	- `unit_test_info.txt`

3. create your own useful content

Other material found in this template

There is an ant build file `build.xml` to build java code. It copies the resulting
jar-file to the directory bin. 

You can also add external java resources in `java/resources/` and external native libraries in `bin_external/`.
	
This template contains a small java class and unit-test. It is small but fully functional.

    java/src/org/openda/NothingUseful.java
    java/test/org/openda/NothingUsefulTest.java

There is also an external test case in the directory tests

    tests/run_test.sh
D-Flow Flexible Mesh is not part of OpenDA and should be installed separately on your computer. 

All tests start the D-Flow Flexible Mesh executable with a shell or bat script located in the directory *./test_dir/stochModel/bin*.

There are 2 scripts in this directory: 
- start_dimr.bat. Use this script if you have installed Delft3D FM Suite where DIMR (Delft Integrated Model Runner)is used to start the dflowfm.dll
- start_dimr.sh. Use this script if you have a Linux installation.
 
If D-Flow Flexible Mesh cannot be started by OpenDA, you probably need to check and modify the installation path or other details in one of these scripts. 

The configuration of the examples in this directory was updated and tested on Windows10 and CentOS Linux version 7:
- OpenDA 3.0.0
- Delft3D Flexible Mesh Suite 2020.02 HM
- dimr 2.18.05.73189

For other distributions or platforms the configuration might not be correct.

	Status: 
	Tests that are enabled in run_all_tests.bat run correctly. 
	Tests that are commented out do not run correctly. WORK IN PROGRESS.


## build executable

Dependencies

```bat
    pip install pyinstaller
```

First run the default option to create a one folder build
```bat
    pyinstaller src/reactive_pollution_model.py
```

To create a single executable use
```bat
    pyinstaller --onefile src/reactive_pollution_model.py
```

The resulting folder/executable are located in the `dist` folder.
Note that when running the option with `--onefile` after running the one folder build you need to clean the `dist` folder manually. 