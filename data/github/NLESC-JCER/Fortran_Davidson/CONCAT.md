# Unreleased

# Version 0.2.0

## Changed

* Improved variables names.
* Removed unnecessary ritz_vectors calculations.
* Optimized the DPR for the free version [#34](https://github.com/NLESC-JCER/Fortran_Davidson/issues/34).
* Removed recomputation of the residues [#34](https://github.com/NLESC-JCER/Fortran_Davidson/issues/34).

# Version 0.1.0

## Changed

* Update the whole projection matrix after adding some correction vector. Replace the block update schema.

## Fixed

* Fixed [several bugs](https://github.com/NLESC-JCER/Fortran_Davidson/issues/29) in the matrix free implementation.

# Version 0.0.5

## New

 * Select the initial orthonormal basis set based on the lowest diagonal elements of the matrix

# Version 0.0.4

## changed
 *  For the Matrix-free implementation, the function that computes the target matrix on the fly,
 was replaced by another function representing the action of a matrix-operator over a vector

# Version 0.0.3

## Changed
 * split the `dense` and `matrix` free into modules
 * Moved the `lapack` calls to their own module
 * Moved the array utilities to their own module
 
## Deleted
 * Removed the `DPR` vs `GJD` benchmark

# Version 0.0.2

## Added
 * Algorithm for the [Generalized eigenvalue problem](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix#Generalized_eigenvalue_problem) using the Davidson method.
 * Matrix free implementation of the Davidson method using the **DPR** correction.
 * Tests for the lapack wrappers.
 
## Changed
 * Append the orthonormalized correction vector to the previous subspace by computing only the
 matrix elements involved with the new correction vector.
 
 * Used as maximum dimension of the subspace 10 times the number of requested eigenvalues.

# Version 0.0.1

## Dependencies
 
 * Fortran compiler >= 6.0 (supporting `submodules` and `move_alloc`)
 * [lapack](http://www.netlib.org/lapack/) or [MKL](https://software.intel.com/en-us/mkl).
 * [Cmake](https://cmake.org/) >= 3.0
[![Build Status](https://travis-ci.org/NLESC-JCER/Fortran_Davidson.svg?branch=master)](https://travis-ci.org/NLESC-JCER/Fortran_Davidson) [![Pipeline](https://github.com/NLESC-JCER/Fortran_Davidson/workflows/FortranCI/badge.svg)](https://github.com/NLESC-JCER/Fortran_Davidson/actions) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3637964.svg)](https://doi.org/10.5281/zenodo.3637964)

Davidson Eigensolver
===================
This package contains a Modern Fortran implementation of the *Davidson diagonalization algorithms*.
Different schemas are available to compute the correction.

Available correction methods are:
 * **DPR**: Diagonal-Preconditioned-Residue
 * **GJD**: Generalized Jacobi Davidson


### Note:
The Davidson method is suitable for **diagonal-dominant symmetric matrices**, that are quite common
in certain scientific problems like [electronic structure](https://en.wikipedia.org/wiki/Electronic_structure). The Davidson method could be not practical
for other kinds of symmetric matrices.

Usage
-----
The following program call the `eigensolver` *subroutine* from the `davidson` module and computes
the lowest 3 eigenvalues and corresponding eigenvectors, using the *GJD* method with a tolerance
of `1e-8` and `100` maximum iteration.
```fortran
program main
  use numeric_kinds, only: dp
  use davidson, only: generalized_eigensolver, generate_diagonal_dominant
 
  implicit none

  integer, parameter :: dim = 50
  integer, parameter :: lowest = 3
  real(dp), dimension(dim, dim) :: matrix, second_matrix
  real(dp), dimension(lowest) :: eigenvalues
  real(dp), dimension(dim, lowest) :: eigenvectors
  real(dp) :: tolerance
  integer:: max_dim_subspace, max_iterations, lowest

  matrix = generate_diagonal_dominant(dim, 1d-4)
  second_matrix = generate_diagonal_dominant(dim, 1d-4, 1.0_dp)
  max_iterations = 1000
  max_dim_subspace = 20
  tolerance = 1d-8
  call generalized_eigensolver(matrix, eigenvalues, eigenvectors, lowest, "GJD", max_iterations, &
       tolerance, final_iterations, max_dim_subspace, second_matrix)
  print *, eigenvalues
  print *, eigenvectors

end program main
```
The helper  `generate_diagonal_dominant` function generates a diagonal dominant
matrix with entries to the diagonal close to row number `(i=1, number_of_rows)`
and random number of the order `1e-4` on the off-diagonal entries.

**Variables**:
 * `matrix` (*in*) matrix to diagonalize
 * `eigenvalues` (*out*) resulting eigenvalues
 * `eigenvectors` (*out*) resulting eigenvectors
 * `lowest`(*in*) number of eigenvalues to compute
 * `method`(*in*) Either "DPR" or "GJD"
 * `max_iterations`(*in*) maximum number of iterations
 * `tolerance`(*in*) Numerical tolerance for convergence
 * `final_iterations`(*output*) returns the number of iterations that were needed to converge
 * `max_dim_subspace`(*in*, *optional*) Dimension of the subspace of search
 * `second_matrix`(*in*, optional) Optional matrix to compute the generalized eigenvalue problem
 
### References:
 * [Davidson diagonalization method and its applications to electronic structure calculations](https://www.semanticscholar.org/paper/DAVIDSON-DIAGONALIZATION-METHOD-AND-ITS-APPLICATION-Liao/5811eaf768d1a006f505dfe24f329874a679ba59)
 * [Numerical Methods for Large Eigenvalue Problem](https://doi.org/10.1137/1.9781611970739)

Installation and Testing
------------------------

To compile execute:
```
cmake -H. -Bbuild && cmake --build build
```

To use another compiler (e.g. ifort):
```
cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=ifort && cmake --build build
```

To run the test:
```
cmake -H. -Bbuild -DENABLE_TEST=ON && cmake --build build
cd build && ctest -V
```

To Debug compile as:
```
cmake -H. -Bbuild  -DCMAKE_BUILD_TYPE=Debug && cmake --build build
```

Dependencies
------------
This packages assumes that you have installed the following packages:
 * A Fortran compiler >=  version 6.0 
 * [CMake](https://cmake.org/)
 * [Lapack](http://www.netlib.org/lapack/)
	
Optionally, If an [MKL](https://software.intel.com/en-us/mkl) library is available the package will try to find it.
 Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance, race,
religion, or sexual identity and orientation.

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
reported by contacting the project team at f.zapata@esciencecenter.nl. All
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
# Contributing guidelines

We welcome any kind of contribution to our software, from simple comment or question to a full fledged [pull request](https://help.github.com/articles/about-pull-requests/). Please read and follow our [Code of Conduct](CODE_OF_CONDUCT.md).

A contribution can be one of the following cases:

1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

## You have a question

1. use the search functionality [here](https://github.com/NLESC-JCER/Fortran_Davidson/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue;
1. apply the "Question" label; apply other labels when relevant.

## You think you may have found a bug

1. use the search functionality [here]() to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the [SHA hashcode](https://help.github.com/articles/autolinked-references-and-urls/#commit-shas) of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
1. apply relevant labels to the newly created issue.

## You want to make some kind of change to the code base

1. (**important**) announce your plan to the rest of the community _before you start working_. This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
1. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/) and [here](https://help.github.com/articles/syncing-a-fork/));
1. make sure the existing tests still work by running ``python setup.py test``;
1. add your own tests (if necessary);
1. update or expand the documentation;
1. [push](http://rogerdudler.github.io/git-guide/) your feature branch to (your fork of) the Fortran_Davidson repository on GitHub;
1. create the pull request, e.g. following the instructions [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
