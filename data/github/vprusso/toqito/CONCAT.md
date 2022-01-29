# Changelog

## 0.0.5

- Fix: Bug in `swap.py`. Added test to cover bug found.

- Fix: Bug in `ppt_distinguishability.py` that prevented checking for higher
  dimensions. Added test to cover bug found.

- Fix: Bug in `ppt_distinguishability.py` that prevented checking 2-dimensional
  cases.

- Fix: Bug in `state_distinguishability.py`. Value returned was being multiplied
  by an unnecessary factor. Added other tests that would have caught this.

- Fix: Bug in `partial_trace.py`. Added test to cover bug found.

- Fix: Bug in `partial_transpose.py` for non-square matrices. Adding test cases to
  cover this previously failing case.
 
- Fix: Bug in `purity.py`. Squaring matrix was using `**2` instead of 
  `np.linalg.matrix_power` 

- Fix: Bug in `random_state_vector.py`. Added tests to cover.

- Fix: Bug in `schmidt_rank.py`. Added tests to cover.

- Feature: Added in `symmetric_extension_hierarchy.py`. The function is a
  hierarchy of semidefinite programs that eventually converge to the separable
  value for distinguishing an ensemble of quantum states. The first level is
  just the standard PPT distinguishability probability. Computing higher levels
  of this hierarchy can become intractable quite quickly.

- Feature: Added in ability to perform both minimum-error and unambiguous state
  discrimination in `state_distinguishability.py`. Adding additional tests to
  `test_state_distinguishability.py` to cover this extra feature.

- Feature: Added `majorizes.py`; a function to determine whether a vector or matrix
  majorizes another. This is used as a criterion to determine whether a quantum
  state can be converted to another via LOCC. Added `test_majorizes.py` for unit
  testing.

- Feature: Added `perfect_matchings.py`; a function that calculates all the 
  perfect matchings of a given list of objects. Added 
  `test_perfect_matchings.py` for unit testing.

- Feature: Added `breuer.py` under `states/`; a specific bound-entangled state of 
  interest. Added in `test_breuer.py` for unit testing.

- Feature: Added `brauer.py` under `states/`. Added in `test_brauer.py` for unit 
  testing.

- Feature: Added `entanglement_of_formation.py` under `state_props/`. Added in 
  `test_entanglement_of_formation.py` for unit testing.

- Feature: Added `l1_norm_coherence.py` under `state_props/`. Added 
  `test_l1_norm_coherence.py` for unit testing.

- Enhancement: Adding further tests for `symmetric_projection.py`.

- Enhancement: Adding further tests for `ppt_distinguishability.py`.

- Enhancement: Adding further tests for `reduction.py`

- Enhancement: Consolidating `conclusive_state_exclusion.py` and 
  `unambiguous_state_exclusion.py` to just `state_exclusion.py`. Adding in
  optional parameter to specify method of exclusion. Consolidated tests in
  `test_state_exclusion.py`.

## 0.0.6

 - Feature: Added `in_separable_ball.py` under `state_props/`. Knowing whether a
 density matrix (quantum state) is contained in the separable ball centered at
 the maximally-mixed state is useful for separability testing. Added 
 `test_in_separable_ball.py` for unit testing.

- Feature: Added in ability to perform both minimum-error and unambiguous state
  discrimination in `ppt_distinguishability.py`. Adding additional tests to
  `ppt_state_distinguishability.py` to cover this extra feature.

- Feature: Added in ability to compute both primal and dual optimization 
  problems in `ppt_distinguishability.py`. This gives the user the capability
  to obtain the measurement operators and also to use the computationally
  efficiency of the dual problem to make quicker numeric calculations.
 
- Feature: Added `unvec.py` under `matrix_ops/`. This feature is the inverse of
  pre-existing "vec" operation found in `vec.py`. That is, it allows one to take
  a vector and transform it to a (possibly square) matrix.

- Feature: Added `choi_to_kraus.py` under `channel_ops/`. This feature allows
  the user to convert a Choi matrix to a list of Kraus operators. This is the 
  inverse of the existing `kraus_to_choi.py` functionality that `toqito` 
  currently offers.

- Feature: Add `is_mutually_orthogonal.py` under `state_props/`. This feature
  allows the user to determine whether a given set of states (vectors) are
  mutually orthogonal with each other.

- Feature: Add `has_symmetric_extension.py` under `state_props/`. This feature 
  is very useful for determining whether a given state is entangled by checking
  whether there exists a symmetric extension for that state.

- Feature: Add `cvx_kron.py` under `helper/`. This feature allows one to compute
  the Kronecker product between two entities, where either two may be of type
  `np.ndarray` or a `cvxpy` expression.
 
- Feature: Add `is_identity.py` to `state_props/`. This feature allows one to check
  if a given matrix is equal to the identity matrix. 
 
- Enhancement: Previously ignore tests for `channel_props/` are now being run as
  part of the testing suite.

- Enhancement: More robust error checking and adding tolerance arguments for
  various matrix properties found in `matrix_props/`.

- Doc-Fix: Trailing back-tick quotes removed in `li_norm_coherence.py`.
 
- Fix: Bug in `symmetric_extension_hierarchy.py` for `level=1` case. Added test
  to cover this bug.

- Fix: Bug in `states/tile.py` produced 3-dimensional vectors when they should
  have in fact been 9-dimensional vectors.

- Fix: Warning in testing suite for `test_symmetric_projection.py` for 
  previously using the deprecated `numpy.matrix` class has since been fixed and
  the warning hence resolved.

- Fix: The `states/isotropic.py` file no longer requires the use of the 
  deprecated `numpy.matrix` class.

## 0.0.7

- Feature: Added `is_unital.py` under `channel_props/`. This feature allows
  the user to determine whether a given channel (specified by either its Choi
  matrix or by its Kraus operators) is unital.
 
- Feature: Added `is_trace_preserving.py` under `channel_props/`. This feature
  allows the user to determine whether a given channel (specified by either its
  Choi matrix or by its Kraus operators) is trace-preserving.

- Feature: Added `log_negativity.py` under `state_props`. This feature
  allows the user to calculate the log negativty of a quantum state.

- Feature: Added `channel_fidelity.py` under `channel_metrics/`. This 
  feature allows the user to calculate the channel fidelity between the Choi 
  representations of two quantum channels.

- Feature: Added `is_idempotent.py` under `matrix_props/`. This feature
  allows the user to determine whether a given matrix is idempotent.

- Enhancement: Adding `rtol` and `atol` tolerance parameters for 
  `is_herm_preserving.py`.

- Enhancement: Improving speed of calculating the classical value of nonlocal
  game. This enhancement is taken from QETLAB which was inspired by pre-print
  arxiv:2005.13418.

- Enhancement: Parallel repetitions for the classical value of a nonlocal game
  is now supported.

- Fix: The `partial_channel.py` function has been enhanced to deal with 
  completely positive maps specified by Kraus operators as input.

- Fix: The GHZ state now supports either dimension or parameter `1` instead 
  of it previously being `2`

## 1.0.0

- Fix: Various documentation fixes and updates.

- Fix: Index error for unambiguous quantum state distinguishability.
 
## 1.0.1

- Feature: Added `bures_distance.py` under `state_metrics/`. This feature
  allows the user to determine whether a given matrix is idempotent.

- Feature: Added `singlet.py` under `states/`. This feature allows one to yield
  a singlet state of dimension `dim`.

- Feature: Added `is_quantum_channel` under `channel_props`. This feature
  allows one to check whether a given input provided as either a list of Kraus
  operators, or, a Choi matrix constitutes a valid quantum channel.

- Fix: Permute systems had a bug where if the `inv_perm` option in 
  `permute_systems.py` was selected, the standard permutation was calculated 
  (not the inverse permutation). Further unit tests are included to catch
  similar failures.

- Fix: The `partial_transpose.py` function did not accurately calculate the 
  partial transpose on matrices of certain dimension. The fix for 
  `permute_systems.py` fixes the issue with `partial_transpose.py`. Further unit
  tests are included to catch similar failures.

- Fix: The `partial_trace.py` function was not accurately calculating the 
  partial trace when the argument was specified as a list of dimensions for 
  certain cases. This has been fixed and further test cases have been included 
  to prevent this from occurring.

- Fix: The `swap.py` function was not accurately swapping on all sub-systems.
  Further unit tests are included to catch similar failures.

- Fix: The `hadamard.py` function was not yielding Hadamard matrices of proper
  size and value. Fixed and added tests to cover this case.

- Fix: The `schmidt_decomposition.py` function was taking an incorrect argument 
  into the SVD function. Fixed and added further tests cases to cover.

- Fix: The `is_product_vector.py` was making use of the 
  `schmidt_decomposition.py` function incorrectly. Fixed and added further 
  test cases.

- Enhancement: Adding ability for `schmidt_rank.py` function to process not 
  just vectors, but also matrices. Adding in unit tests to cover this case.

- Enhancement: Adding the ability for `is_product.py` function to process not
  just vectors, but also matrices. Adding in unit tests to cover this case.

- Enhancement: Simplified code for `nonlocal_game.py` and 
  `extended_nonlocal_game.py`

- Enhancement: Some general documentation clean-up. 

## 1.0.2

- Feature: Added `dual_channel` under `channel_ops`. This feature
  allows one to calculate the dual channel of a given input provided as either a 
  list of Kraus operators, or, a Choi matrix that represents a valid quantum 
  channel.
  
- Fix: Thanks to Jake Xuereb for a fix in the `concurrence` fucntion that was 
  incorrectly not taking the transpose.# ![logo](./docs/figures/logo.svg "logo")

(Theory of Quantum Information Toolkit)

The `toqito` package is an open source Python library for studying various
objects in quantum information, namely, states, channels, and measurements.

Specifically, `toqito` focuses on providing numerical tools to study problems
pertaining to entanglement theory, nonlocal games, matrix analysis, and other
aspects of quantum information that are often associated with computer science.

`toqito` aims to fill the needs of quantum information researchers who want
numerical and computational tools for manipulating quantum states,
measurements, and channels. It can also be used as a tool to enhance the
experience of students and instructors in classes pertaining to quantum
information.


[![build status](http://img.shields.io/travis/vprusso/toqito.svg?style=plastic)](https://travis-ci.org/vprusso/toqito)
[![doc status](https://readthedocs.org/projects/toqito/badge/?version=latest&style=plastic)](https://toqito.readthedocs.io/en/latest/)
[![codecov](https://codecov.io/gh/vprusso/toqito/branch/master/graph/badge.svg?style=plastic)](https://codecov.io/gh/vprusso/toqito)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4743211.svg)](https://doi.org/10.5281/zenodo.4743211)
[![Unitary Fund](https://img.shields.io/badge/Supported%20By-UNITARY%20FUND-brightgreen.svg?style=plastic)](http://unitary.fund)

## Installing

The preferred way to install the `toqito` package is via `pip`:

```
pip install toqito
```

Alternatively, to install, you may also run the following command from the
top-level package directory.

```
python setup.py install
```

## Using

Full documentation along with specific examples and tutorials are provided
here: [https://toqito.readthedocs.io/](https://toqito.readthedocs.io/).

More information can also be found on the following
[toqito homepage](https://vprusso.github.io/toqito/).

## Testing

The `pytest` module is used for testing. To run the suite of tests for `toqito`,
run the following command in the root directory of this project.

```
pytest --cov-report term-missing --cov=toqito tests/
```

## Citing

You can cite `toqito` using the following DOI:
10.5281/zenodo.4743211


If you are using the `toqito` software package in research work, please include
an explicit mention of `toqito` in your publication. Something along the lines
of:

```
To solve problem "X" we used `toqito`; a package for studying certain
aspects of quantum information.
```

A BibTeX entry that you can use to cite `toqito` is provided here:

```bib
@misc{toqito,
   author       = {Vincent Russo},
   title        = {toqito: A {P}ython toolkit for quantum information, version 1.0.0},
   howpublished = {\url{https://github.com/vprusso/toqito}},
   month        = May,
   year         = 2021,
   doi          = {10.5281/zenodo.4743211}
 }
```

## Contributing

All contributions, bug reports, bug fixes, documentation improvements,
enhancements, and ideas are welcome.

A detailed overview on how to contribute can be found in the
[contributing guide](https://github.com/vprusso/toqito/blob/master/.github/CONTRIBUTING.md).

## License

[MIT License](http://opensource.org/licenses/mit-license.php>)
---
title: 'toqito -- Theory of quantum information toolkit: A Python package for studying quantum information'
tags:
  - Python
  - quantum information
  - quantum computing
  - entanglement
  - nonlocal games
  - matrix analysis
authors:
  - name: Vincent Russo
    affiliation: "1, 2"
affiliations:
  - name: ISARA Corporation
    index: 1    
  - name: Modellicity Inc.
    index: 2
date: 21 Mar 2021
bibliography: paper.bib

---

# Summary

`toqito` is an open source library for studying various objects in quantum
information, namely: states, channels, and measurements. `toqito` focuses on
providing numerical tools to study problems pertaining to: entanglement theory,
nonlocal games, matrix analysis, and other aspects of quantum information that
are often associated with computer science. 

# Statement of Need

While there are many outstanding feature-rich Python packages to study quantum
information, they are often focused on physics applications 
[@johansson2013qutip; @killoran2019strawberry; @steiger2018projectq].
Other excellent software offerings that are closer in scope to `toqito`, such
as `QETLAB` [@johnston2016qetlab], are written in non-opensource languages and
therefore require the users to have access to costly licenses.

`toqito` possesses functions for fundamental operations including the partial
trace, partial transpose, and others. `toqito` also makes use of the
`cvxpy` [@diamond2016cvxpy] convex optimization module to solve various
semidefinite programs that pertain to problems arising in the study of nonlocal
games, state discrimination, and other problems in quantum information.
`toqito` provides the ability to either directly calculate or estimate the
classical and quantum values of nonlocal games. `toqito` also provides numerous
functions for performing operations on and for determining properties of
quantum states, channels, and measurements. `toqito` provides utilities for
exploring measures of entanglement and properties of entangled states. Support
for generating random quantum states and measurement operators is also
provided. 

The `toqito` module is supported for Python 3.7 and makes use of many of the more
modern features of the language including f-strings, type hinting, and others.
`toqito` is available on GitHub (https://github.com/vprusso/toqito) and can be
installed via pip (https://pypi.org/project/toqito/). Further information of
features and uses can be found on the documentation page
(https://toqito.readthedocs.io/en/latest/).

# Acknowledgements
This research is supported by the Unitary Fund [@zeng2019unitary].

# References

## Expected Behavior


## Actual Behavior


## Steps to Reproduce the Problem

1.
1.
1.

## Specifications

- Version:
- Platform:
- Subsystem:# How to contribute

We welcome contributions from external contributors, and this document
describes how to merge code changes into this `toqito`. 

## Getting Started

*    Make sure you have a [GitHub account](https://github.com/signup/free).
*    [Fork](https://help.github.com/articles/fork-a-repo/) this repository on GitHub.
*    On your local machine,
     [clone](https://help.github.com/articles/cloning-a-repository/) your fork of
     the repository.

## Making Changes

*    Add some really awesome code to your local fork.  It's usually a 
     [good idea](http://blog.jasonmeridth.com/posts/do-not-issue-pull-requests-from-your-master-branch/)
     to make changes on a 
     [branch](https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/)
     with the branch name relating to the feature you are going to add.
*    When you are ready for others to examine and comment on your new feature,
     navigate to your fork of `toqito` on GitHub and open a 
     [pull request](https://help.github.com/articles/using-pull-requests/) (PR). Note that
     after you launch a PR from one of your fork's branches, all
     subsequent commits to that branch will be added to the open pull request
     automatically.  Each commit added to the PR will be validated for
     mergability, compilation and test suite compliance; the results of these tests
     will be visible on the PR page.
*    If you're providing a new feature, you must add test cases and documentation.
*    When the code is ready to go, make sure you run the test suite using pytest.
*    When you're ready to be considered for merging, check the "Ready to go"
     box on the PR page to let the `toqito` devs know that the changes are complete.
     The code will not be merged until this box is checked, the continuous
     integration returns check marks,
     and the primary developer approves the reviews.

# Additional Resources

*    [General GitHub documentation](https://help.github.com/)
*    [PR best practices](http://codeinthehole.com/writing/pull-requests-and-other-good-practices-for-teams-using-github/)
*    [A guide to contributing to software packages](http://www.contribution-guide.org)
*    [Thinkful PR example](http://www.thinkful.com/learn/github-pull-request-tutorial/#Time-to-Submit-Your-First-PR)

## Code of Conduct

### Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of 
experience, nationality, personal appearance, race, religion, or sexual identity
and orientation.

### Our Standards

Examples of behavior that contributes to creating a positive environment
include:

*    Using welcoming and inclusive language
*    Being respectful of differing viewpoints and experiences
*    Gracefully accepting constructive criticism
*    Focusing on what is best for the community
*    Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

*    The use of sexualized language or imagery and unwelcome sexual attention or 
     advances
*    Trolling, insulting/derogatory comments, and personal or political attacks
*    Public or private harassment
*    Publishing others' private information, such as a physical or electronic
     address, without explicit permission
*    Other conduct which could reasonably be considered inappropriate in a
     professional setting

### Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

### Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

### Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at [vincentrusso1@gmail.com]. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

### Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
## Description
Provide a brief description of the PR's purpose here.

## Todos
Notable points that this PR has either accomplished or will accomplish.
  -  [ ] TODO 1

## Questions
-  [ ] Question1

## Status
-  [ ] Ready to go