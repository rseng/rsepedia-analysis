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
-  [ ] Ready to goQuantum state distinguishability
=================================

In this tutorial we are going to cover the problem of *quantum state
distinguishability* (sometimes analogously referred to as quantum state
discrimination). We are going to briefly describe the problem setting and then
describe how one may use :code:`toqito` to calculate the optimal probability
with which this problem can be solved when given access to certain
measurements.

Further information beyond the scope of this tutorial can be found in the text
[tWatrousQI]_ as well as the course [tSikoraSDP]_.

The state distinguishability problem
-------------------------------------

The quantum state distinguishability problem is phrased as follows.

1. Alice possesses an ensemble of :math:`n` quantum states:

    .. math::
        \begin{equation}
            \eta = \left( (p_0, \rho_0), \ldots, (p_n, \rho_n)  \right),
        \end{equation}

where :math:`p_i` is the probability with which state :math:`\rho_i` is
selected from the ensemble. Alice picks :math:`\rho_i` with probability
:math:`p_i` from her ensemble and sends :math:`\rho_i` to Bob.

2. Bob receives :math:`\rho_i`. Both Alice and Bob are aware of how the
   ensemble is defined but he does *not* know what index :math:`i`
   corresponding to the state :math:`\rho_i` he receives from Alice is.

3. Bob wants to guess which of the states from the ensemble he was given. In
   order to do so, he may measure :math:`\rho_i` to guess the index :math:`i`
   for which the state in the ensemble corresponds.

This setting is depicted in the following figure.

.. figure:: figures/quantum_state_distinguish.svg
   :alt: quantum state distinguishability
   :align: center

   The quantum state distinguishability setting.

Depending on the sets of measurements that Alice and Bob are allowed to use,
the optimal probability of distinguishing a given set of states is characterized
by the following image.

.. figure:: figures/measurement_inclusions.svg
   :width: 200
   :alt: measurement hierarchy
   :align: center

   Measurement hierarchy.

That is,the probability that Alice and Bob are able to distinguish using PPT
measurements is a natural upper bound on the optimal probability of
distinguishing via separable measurements.

In general:

* LOCC: These are difficult objects to handle mathematically; difficult to
  design protocols for and difficult to provide bounds on their power.

* Separable: Separable measurements have a nicer structure than LOCC.
  Unfortunately, optimizing over separable measurements in NP-hard.

* PPT: PPT measurements offer a nice structure and there exists efficient
  techniques that allow one to optimize over the set of PPT measurements via
  semidefinite programming.

Optimal probability of distinguishing a quantum state
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The optimal probability with which Bob can distinguish the state he is given
may be obtained by solving the following semidefinite program (SDP).

.. math::
    \begin{align*}
        \text{maximize:} \quad & \sum_{i=0}^n p_i \langle M_i,
        \rho_i \rangle \\
        \text{subject to:} \quad & \sum_{i=0}^n M_i = \mathbb{I}_{\mathcal{X}},\\
                                 & M_i \in \text{Pos}(\mathcal{X})
    \end{align*}

This optimization problem is solved in :code:`toqito` to obtain the optimal
probability with which Bob can distinguish state :math:`\rho_i`.

To illustrate how we can phrase and solve this problem in :code:`toqito`,
consider the following example. Assume Alice has an ensemble of quantum states

.. math::
    \eta = \{ (1/2, \rho_0), (1/2, \rho_1) \}

such that 

.. math::
    \rho_0 = | 0 \rangle \langle 0 | = \begin{pmatrix}
                1 & 0 \\
                0 & 0
             \end{pmatrix} \quad \text{and} \quad
    \rho_1 = | 1 \rangle \langle 1 | = \begin{pmatrix}
                0 & 0 \\
                0 & 1
             \end{pmatrix}


These states are completely orthogonal to each other, and it is known that Bob
can optimally distinguish the state he is given perfectly, i.e. with probability
:math:`1`.

Using :code:`toqito`, we can calculate this probability directly as follows:

.. code-block:: python

    >>> from toqito.states import basis
    >>> from toqito.state_opt import state_distinguishability
    >>> 
    >>> # Define the standard basis |0> and |1>
    >>> e_0, e_1 = basis(2, 0), basis(2, 1)
    >>>
    >>> # Define the corresponding density matrices of |0> and |1> 
    >>> # given as |0><0| and |1><1|, respectively.
    >>> e_00 = e_0 * e_0.conj().T
    >>> e_11 = e_1 * e_1.conj().T
    >>>
    >>> # Define a list of states and a corresponding list of 
    >>> # probabilities with which those states are selected.
    >>> states = [e_00, e_11] 
    >>> probs = [1/2, 1/2]
    >>>
    >>> # Calculate the probability with which Bob can 
    >>> # distinguish the state he is provided.
    >>> state_distinguishability(states, probs)
    1.000

Specifying similar state distinguishability problems can be done so using this
general pattern.

Optimal probability of distinguishing a state via PPT measurements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We may consider the quantum state distinguishability scenario under somewhat
different and more limited set of circumstances. Specifically, we may want to
ask the same question but restrict to enforcing that in order to determine the
state that Bob is provided, he is limited to using a certain class of
measurement. There are a wider class of measurements with respect to the ones
we considered in the previous example referred to as PPT (positive partial
transpose).

The problem of state distinguishability with respect to PPT measurements can
also be framed as an SDP and was initially presented in this manner in
[tCosentino13]_

.. math::

    \begin{equation}
        \begin{aligned}
            \text{minimize:} \quad & \frac{1}{k} \text{Tr}(Y) \\
            \text{subject to:} \quad & Y \geq \text{T}_{\mathcal{A}}
                                      (\rho_j), \quad j = 1, \ldots, k, \\
                                     & Y \in \text{Herm}(\mathcal{A} \otimes
                                      \mathcal{B}).
        \end{aligned}
    \end{equation}

Using :code:`toqito`, we can determine the optimal probability for Bob to
distinguish a given state from an ensemble if he is only given access to PPT
measurements.

Consider the following Bell states

.. math::
    \begin{equation}
        \begin{aligned}
            | \psi_0 \rangle = \frac{|00\rangle + |11\rangle}{\sqrt{2}}, &\quad
            | \psi_1 \rangle = \frac{|01\rangle + |10\rangle}{\sqrt{2}}, \\
            | \psi_2 \rangle = \frac{|01\rangle - |10\rangle}{\sqrt{2}}, &\quad
            | \psi_3 \rangle = \frac{|00\rangle - |11\rangle}{\sqrt{2}}.
        \end{aligned}
    \end{equation}

It was shown in [tCosentino13]_ and later extended in [tCR13]_ that for the following set of states

.. math::
    \begin{equation}
        \begin{aligned}
            \rho_1^{(2)} &= |\psi_0 \rangle | \psi_0 \rangle \langle \psi_0 | \langle \psi_0 |, \\
            \rho_2^{(2)} &= |\psi_1 \rangle | \psi_3 \rangle \langle \psi_1 | \langle \psi_3 |, \\
            \rho_3^{(2)} &= |\psi_2 \rangle | \psi_3 \rangle \langle \psi_2 | \langle \psi_3 |, \\
            \rho_4^{(2)} &= |\psi_3 \rangle | \psi_3 \rangle \langle \psi_3 | \langle \psi_3 |, \\
        \end{aligned}
    \end{equation}

that the optimal probability of distinguishing via a PPT measurement should yield
:math:`7/8 \approx 0.875`.

This ensemble of states and some of its properties with respect to
distinguishability were initially considered in [tYDY12]_. In :code:`toqito`,
we can calculate the probability with which Bob can distinguish these states
via PPT measurements in the following manner.

.. code-block:: python

    >>> import numpy as np
    >>> from toqito.states import bell
    >>> from toqito.state_opt import ppt_distinguishability
    >>> # Bell vectors:
    >>> psi_0 = bell(0)
    >>> psi_1 = bell(2)
    >>> psi_2 = bell(3)
    >>> psi_3 = bell(1)
    >>>
    >>> # YDY vectors from [tYDY12]_:
    >>> x_1 = np.kron(psi_0, psi_0)
    >>> x_2 = np.kron(psi_1, psi_3)
    >>> x_3 = np.kron(psi_2, psi_3)
    >>> x_4 = np.kron(psi_3, psi_3)
    >>>
    >>> # YDY density matrices:
    >>> rho_1 = x_1 * x_1.conj().T
    >>> rho_2 = x_2 * x_2.conj().T
    >>> rho_3 = x_3 * x_3.conj().T
    >>> rho_4 = x_4 * x_4.conj().T
    >>>
    >>> states = [rho_1, rho_2, rho_3, rho_4]
    >>> probs = [1 / 4, 1 / 4, 1 / 4, 1 / 4]
    >>> ppt_distinguishability(states, probs)
    0.875

Probability of distinguishing a state via separable measurements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As previously mentioned, optimizing over the set of separable measurements is
NP-hard. However, there does exist a hierarchy of semidefinite programs which
eventually does converge to the separable value. This hierarchy is based off
the notion of symmetric extensions. More information about this hierarchy of
SDPs can be found here [tNav08]_.

References
------------------------------
.. [tWatrousQI] Watrous, John
    "The theory of quantum information"
    Section: "A semidefinite program for optimal measurements"
    Cambridge University Press, 2018

.. [tNav08] Navascués, Miguel.
    "Pure state estimation and the characterization of entanglement."
    Physical review letters 100.7 (2008): 070503.
    https://arxiv.org/abs/0707.4398

.. [tSikoraSDP] Sikora, Jamie
    "Semidefinite programming in quantum theory (lecture series)"
    Lecture 2: Semidefinite programs for nice problems and popular functions
    Perimeter Institute for Theoretical Physics, 2019

.. [tCosentino13] Cosentino, Alessandro,
    "Positive-partial-transpose-indistinguishable states via semidefinite programming",
    Physical Review A 87.1 (2013): 012321.
    https://arxiv.org/abs/1205.1031

.. [tCR13] Cosentino, Alessandro and Russo, Vincent
    "Small sets of locally indistinguishable orthogonal maximally entangled states",
    Quantum Information & Computation, Volume 14, 
    https://arxiv.org/abs/1307.3232

.. [tYDY12] Yu, Nengkun, Runyao Duan, and Mingsheng Ying.
    "Four locally indistinguishable ququad-ququad orthogonal
    maximally entangled states."
    Physical review letters 109.2 (2012): 020506.
    https://arxiv.org/abs/1107.3224
Introductory Tutorial
======================

This tutorial will illustrate the basics of how to use :code:`toqito`. This
will cover how to instantiate and use the fundamental objects that
:code:`toqito` provides; namely quantum states, channels, and measurements.

This is a user guide for :code:`toqito` and is not meant to serve as an
introduction to quantum information. For introductory material on quantum
information, please consult "Quantum Information and Quantum Computation" by
Nielsen and Chuang or the freely available lecture notes `"Introduction to
Quantum Computing"
<https://cs.uwaterloo.ca/~watrous/LectureNotes/CPSC519.Winter2006/all.pdf)>`_
by John Watrous.

More advanced tutorials can be found on the `tutorials page
<https://toqito.readthedocs.io/en/latest/tutorials.html>`_.

This tutorial assumes you have :code:`toqito` installed on your machine. If you
do not, please consult the `installation instructions
<https://toqito.readthedocs.io/en/latest/install.html>`_.

States
------

A *quantum state* is a density operator

.. math::
    \rho \in \text{D}(\mathcal{X})

where :math:`\mathcal{X}` is a complex Euclidean space and where
:math:`\text{D}(\cdot)` represents the set of density matrices, that is, the
set of matrices that are positive semidefinite with trace equal to :math:`1`.

Quantum States
^^^^^^^^^^^^^^

A complete overview of the scope of quantum states can be found
`here <https://toqito.readthedocs.io/en/latest/states.html#quantum-states>`_

The standard basis bra vectors given as :math:`|0\rangle` and :math:`|1\rangle` where

.. math::
    | 0 \rangle = [1, 0]^{\text{T}} \quad \text{and} \quad | 1 \rangle = [0, 1]^{\text{T}}

can be defined in :code:`toqito` as such

.. code-block:: python

    >>> from toqito.states import basis
    >>> # |0>
    >>> basis(2, 0)
    [[1]
    [0]]
    >>> # |1>
    >>> basis(2, 1)
    [[0]
    [1]]

One may define one of the four Bell states written as

.. math::
    u_0 = \frac{1}{\sqrt{2}} \left(| 00 \rangle + | 11 \rangle \right)

using :code:`toqito` as

.. code-block:: python

    >>> import numpy as np
    >>> e_0, e_1 = basis(2, 0), basis(2, 1)
    >>> u_0 = 1/np.sqrt(2) * (np.kron(e_0, e_0) + np.kron(e_1, e_1))
    [[0.70710678],
     [0.        ],
     [0.        ],
     [0.70710678]]

The corresponding density operator of :math:`u_0` can be obtained from

.. math::
    \rho_0 = u_0 u_0^* = \frac{1}{2} 
    \begin{pmatrix} 
        1 & 0 & 0 & 1 \\
        0 & 0 & 0 & 0 \\
        0 & 0 & 0 & 0 \\
        1 & 0 & 0 & 1
    \end{pmatrix} \in \text{D}(\mathcal{X}).

In :code:`toqito`, that can be obtained as 

.. code-block:: python
    
    >>> rho_0 = u_0 * u_0.conj().T
     [0.  0.  0.  0. ]
     [0.  0.  0.  0. ]
     [0.5 0.  0.  0.5]]

Alternatively, we may leverage the :code:`bell` function in :code:`toqito` to
generate all four Bell states defined as

.. math::
    \begin{equation}
        \begin{aligned}
            u_0 = \frac{1}{\sqrt{2}} \left(| 00 \rangle + | 11 \rangle \right), &\quad 
            u_1 = \frac{1}{\sqrt{2}} \left(| 00 \rangle - | 11 \rangle \right), \\
            u_2 = \frac{1}{\sqrt{2}} \left(| 01 \rangle + | 10 \rangle \right), &\quad
            u_3 = \frac{1}{\sqrt{2}} \left(| 01 \rangle - | 10 \rangle \right),
        \end{aligned}
    \end{equation}

in a more concise manner as 

.. code-block:: python

    >>> from toqito.states import bell
    >>> import numpy as np
    >>> bell(0)
    [[0.70710678],
     [0.        ],
     [0.        ],
     [0.70710678]]

The Bell states constitute one such well-known class of quantum states. There
are many other classes of states that are widely used in the field of quantum
information. For instance, the GHZ state 

.. math::
    | GHZ \rangle = \frac{1}{\sqrt{2}} \left( | 000 \rangle + | 111 \rangle \right)

is a well-known 3-qubit quantum state. We can invoke this using :code:`toqito` as

.. code-block:: python

    >>> from toqito.states import ghz
    >>> ghz(2, 3).toarray()
    [[0.70710678],
     [0.        ],
     [0.        ],
     [0.        ],
     [0.        ],
     [0.        ],
     [0.        ],
     [0.70710678]]

While the 3-qubit form of the GHZ state is arguably the most notable, it is
possible to define a generalized GHZ state

.. math::
    | GHZ_n \rangle = \frac{1}{\sqrt{n}} \left( | 0 \rangle^{\otimes n} + | 1
    \rangle^{\otimes n} \right).

This generalized state may be obtained in :code:`toqito` as well. For instance,
here is the GHZ state :math:`\mathbb{C}^{4^{\otimes 7}}` as 

.. math::
    \frac{1}{\sqrt{30}} \left(| 0000000 \rangle + 2| 1111111 \rangle + 3|
    2222222 \rangle + 4| 3333333\rangle \right).

Using :code:`toqito`, we can see this generates the appropriate generalized GHZ
state.

.. code-block:: python

    >>> from toqito.states import ghz
    >>> ghz(4, 7, np.array([1, 2, 3, 4]) / np.sqrt(30)).toarray()
    [[0.18257419],
     [0.        ],
     [0.        ],
     ...,
     [0.        ],
     [0.        ],
     [0.73029674]])

Properties of Quantum States
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given a quantum state, it is often useful to be able to determine certain
*properties* of the state.

For instance, we can check if a quantum state is pure, that is, if the density
matrix that describes the state has rank 1.

Any one of the Bell states serve as an example of a pure state

.. code-block:: python

    >>> from toqito.states import bell
    >>> from toqito.state_props import is_pure
    >>> rho = bell(0) * bell(0).conj().T
    >>> is_pure(rho)
    True

Another property that is useful is whether a given state is PPT (positive
partial transpose), that is, whether the state remains positive after taking
the partial transpose of the state.

For quantum states consisting of shared systems of either dimension :math:`2
\otimes 2` or :math:`2 \otimes 3`, the notion of whether a state is PPT serves
as a method to determine whether a given quantum state is entangled or
separable.

As an example, any one of the Bell states constitute a canonical maximally
entangled state over :math:`2 \otimes 2` and therefore should not satisfy the
PPT criterion.

.. code-block:: python

    >>> from toqito.states import bell
    >>> from toqito.state_props import is_ppt
    >>> rho = bell(2) * bell(2).conj().T
    >>> is_ppt(rho)
    False

As we can see, the PPT criterion is :code:`False` for an entangled state in
:math:`2 \otimes 2`.

Further properties that one can check via :code:`toqito` may be found `on this page
<https://toqito.readthedocs.io/en/latest/states.html#properties-of-quantum-states>`_.


Operations on Quantum States
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(Coming soon).

Distance Metrics for Quantum States
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given two quantum states, it is often useful to have some way in which to
quantify how similar or different one state is from another.

One well known metric is the *fidelity* function defined for two quantum
states. For two states :math:`\rho` and :math:`\sigma`, one defines the
fidelity between :math:`\rho` and :math:`\sigma` as 

.. math::
    || \sqrt{\rho} \sqrt{\sigma} ||_1,

where :math:`|| \cdot ||_1` denotes the trace norm. 

The fidelity fucntion yields a value between :math:`0` and :math:`1`, with
:math:`0` representing the scenarion where :math:`\rho` and :math:`\sigma` are
as different as can be and where a value of :math:`1` indicates a scenario
where :math:`\rho` and :math:`\sigma` are identical.

Let us consider an example in :code:`toqito` where we wish to calculate the
fidelity function between quantum states that happen to be identitcal.

.. code-block:: python

    >>> from toqito.states import bell
    >>> from toqito.state_metrics import fidelity
    >>>
    >>> # Define two identical density operators.
    >>> rho = bell(0)*bell(0).conj().T
    >>> sigma = bell(0)*bell(0).conj().T
    >>> 
    >>> # Calculate the fidelity between `rho` and `sigma`
    >>> fidelity(rho, sigma)
    1

There are a number of other metrics one can compute on two density matrices
including the trace norm, trace distance. These and others are also available
in :code:`toqito`. For a full list of distance metrics one can compute on
quatum states, consult the docs.

Channels
--------

A *quantum channel* can be defined as a completely positive and trace
preserving linear map.

More formally, let :math:`\mathcal{X}` and :math:`\mathcal{Y}` represent
complex Euclidean spaces and let :math:`\text{L}(\cdot)` represent the set of
linear operators. Then a quantum channel, :math:`\Phi` is defined as

.. math::
    \Phi: \text{L}(\mathcal{X}) \rightarrow \text{L}(\mathcal{Y})

such that :math:`\Phi` is completely positive and trace preserving.

Quantum Channels
^^^^^^^^^^^^^^^^
(Coming soon).

Properties of Quantum Channels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
(Coming soon).

Operations on Quantum Channels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
(Coming soon).

Measurements
------------

A *measurement* can be defined as a function

.. math::
    \mu: \Sigma \rightarrow \text{Pos}(\mathcal{X})

satisfying

.. math::
    \sum_{a \in \Sigma} \mu(a) = \mathbb{I}_{\mathcal{X}}

where :math:`\Sigma` represents a set of measurement outcomes and where
:math:`\mu(a)` represents the measurement operator associated with outcome
:math:`a \in \Sigma`.

States
=====================

A *quantum state* is a density operator

.. math::
    \rho \in \text{D}(\mathcal{X})

where :math:`\mathcal{X}` is a complex Euclidean space and where
:math:`\text{D}(\cdot)` represents the set of density matrices.

Distance Metrics for Quantum States
-----------------------------------

.. toctree::

.. autosummary::
   :toctree: _autosummary

    toqito.state_metrics.fidelity
    toqito.state_metrics.helstrom_holevo
    toqito.state_metrics.hilbert_schmidt
    toqito.state_metrics.sub_fidelity
    toqito.state_metrics.trace_distance
    toqito.state_metrics.trace_norm
    toqito.state_metrics.bures_distance
    toqito.state_metrics.matsumoto_fidelity

Optimizations over Quantum States
-----------------------------------------------

.. toctree::

.. autosummary::
   :toctree: _autosummary

    toqito.state_opt.optimal_clone
    toqito.state_opt.ppt_distinguishability
    toqito.state_opt.state_distinguishability
    toqito.state_opt.state_exclusion
    toqito.state_opt.state_helper
    toqito.state_opt.symmetric_extension_hierarchy

Operations on Quantum States
----------------------------

.. toctree::

.. autosummary::
   :toctree: _autosummary

    toqito.state_ops.pure_to_mixed
    toqito.state_ops.schmidt_decomposition

Properties of Quantum States
----------------------------

.. toctree::

.. autosummary::
   :toctree: _autosummary

    toqito.state_props.concurrence
    toqito.state_props.entanglement_of_formation
    toqito.state_props.has_symmetric_extension
    toqito.state_props.in_separable_ball
    toqito.state_props.is_ensemble
    toqito.state_props.is_mixed
    toqito.state_props.is_mutually_orthogonal
    toqito.state_props.is_mutually_unbiased_basis
    toqito.state_props.is_ppt
    toqito.state_props.is_product
    toqito.state_props.is_pure
    toqito.state_props.l1_norm_coherence
    toqito.state_props.log_negativity
    toqito.state_props.negativity
    toqito.state_props.purity
    toqito.state_props.schmidt_rank
    toqito.state_props.von_neumann_entropy

Quantum States
--------------

.. toctree::

.. autosummary::
   :toctree: _autosummary

    toqito.states.basis
    toqito.states.bell
    toqito.states.brauer
    toqito.states.breuer
    toqito.states.chessboard
    toqito.states.domino
    toqito.states.gen_bell
    toqito.states.ghz
    toqito.states.gisin
    toqito.states.horodecki
    toqito.states.isotropic
    toqito.states.max_entangled
    toqito.states.max_mixed
    toqito.states.tile
    toqito.states.w_state
    toqito.states.werner
    toqito.states.singlet
Nonlocal games
================================================================

In this tutorial, we are going to cover the notion of a *nonlocal game*; a
mathematical framework that abstractly models a physical system. The simplest
instance of a nonlocal game involves two players, Alice and Bob, who are not
allowed to communicate with each other once the game has started and who play
cooperatively against an adversary referred to as the referee.

A primary challenge that arises when studying these games is to determine the
maximum probability with which Alice and Bob are able to achieve a winning
outcome. 

This probability is highly dependent on the type of *strategy* that Alice and
Bob use in the game. A *classical strategy* is one in which Alice and Bob have
access to classical resources. The best that Alice and Bob can do using a
classical strategy is known as the *classical value* of the game. Similarly, a
*quantum strategy* is one in which Alice and Bob have access to quantum
resources. The best that Alice and Bob can do using a quantum strategy is known
as the *quantum value* of the game.

Calculating the classical value of a game is NP-hard as we need to perform a
brute-force check to see which strategy yields the classical value of the game. 

Using :code:`toqito`, we will be able to specify a nonlocal game and be able to
directly calculate the classical value and also place lower bounds on the
quantum value.

Further information beyond the scope of this tutorial on nonlocal games can be
found in [tCHTW04]_. Further information on the lower bound technique can be found in
[tLD07]_.

Two-player nonlocal games
--------------------------

A *two-player nonlocal game* consists of players that we give the names *Alice*
and *Bob*:

.. figure:: figures/alice_and_bob.svg
   :alt: nonlocal game
   :align: center

   The players: Alice and Bob.

Alice and Bob are in separate locations and cannot communicate once the game
begins. Prior to the game however, Alice and Bob are free to communicate with
each other. In addition to the players, there is also another party in this
game that is referred to as the *referee*.

.. figure:: figures/referee.svg
   :alt: nonlocal game
   :align: center

   The referee.

Alice and Bob want to play in a cooperative fashion against the referee.

Now that we have set the stage with respect to the actors and actresses we will
encounter in this game, let us see how the game is actually played.

.. figure:: figures/nonlocal_game.svg
   :alt: nonlocal game
   :align: center

   A two-player nonlocal game.

A nonlocal game unfolds in the following manner.

1. The referee randomly generates questions denoted as :math:`x` and :math:`y`.
   The referee sends the question :math:`x` to Alice and the question :math:`y`
   to Bob. The referee also keeps a copy of :math:`x` and :math:`y` for
   reference.

2. Alice and Bob each receive their respective questions. They are then each
   expected to respond to their questions with answers that we denote as
   :math:`a` and :math:`b`. Alice sends :math:`a` to the referee, and Bob sends
   :math:`b`.

3. When the referee receives :math:`a` and :math:`b` from Alice and Bob, the
   referee evaluates a particular function that is predicated on the questions
   :math:`x` and :math:`y` as well as the answers :math:`a` and :math:`b`. The
   outcome of this function is either :math:`0` or :math:`1`, where an outcome
   of :math:`0` indicates a loss for Alice and Bob and an outcome of :math:`1`
   indicates a win for Alice and Bob.

Alice and Bob's goal in the above game is to get the function in Step-3 to
output a :math:`1`, or equivalently, to indicate a winning outcome. This type
of game is referred to as a *nonlocal game*.

Classical and Quantum Strategies
--------------------------------

Now that we have the framework for a nonlocal game, we can consider the
player's *strategy*; how the players play the game given access to certain
resources. There are a number of strategies that the players can use, but for
simplicity, we will restrict our attention to two types of strategies.

1. *Classical strategies*: The players answer the questions in a deterministic
   manner.

2. *Quantum strategies*: The players make use of quantum resources in the form
   of a shared quantum state and respective sets of measurements.

Classical strategies
^^^^^^^^^^^^^^^^^^^^

A *classical strategy* for a nonlocal game is one where the players
deterministically produce an output for every possible combination of inputs
they may receive in the game. The *classical value* of a nonlocal game is the
maximum probability achieved by the players over all classical strategies. For
a nonlocal game, :math:`G`, we use :math:`\omega(G)` to represent the classical
value of :math:`G`.

One question you may have is whether a classical strategy can be improved by
introducing randomness. If the players randomly select their answers, is it
possible for them to do potentially better than if they had just played
deterministically? As it happens, probabilistic classical strategies cannot
perform any better than deterministic classical strategies.

There is therefore no loss in generality in restricting our analysis of
classical strategies to deterministic ones and it is assumed that when we use
the term classical strategy that we implicitly mean a classical strategy that
is played deterministically.

Quantum strategies
^^^^^^^^^^^^^^^^^^

A *quantum strategy* for a nonlocal game is one where the players prepare a
quantum state prior to the start of the game along with respective sets of
measurements that they apply to their respective portions of the shared state
during the game based on the questions they receive to generate their answers.
The *quantum value* of a nonlocal game is the maximum probability achieved by
the players over all quantum strategies. For a nonlocal game, :math:`G`, we use
:math:`\omega^*(G)` to represent the quantum value of :math:`G`.

.. figure:: figures/nonlocal_game_quantum_strategy.svg
   :alt: nonlocal game quantum strategy
   :align: center

   A two-player nonlocal game invoking a quantum strategy.

Let us describe the high-level steps for how Alice and Bob play using a quantum
strategy.

1. Alice and Bob prepare a state :math:`\sigma \in \text{D}(\mathcal{U} \otimes
   \mathcal{V})` prior to the start of the game. We use :math:`\textsf{U}`` and
   :math:`\textsf{V}` to denote the respective registers of spaces :math:`\textsf{U}`
   and :math:`\textsf{V}`.

2. The referee sends question :math:`x \in \Sigma_A` to Alice and :math:`y \in
   \Sigma_B` to Bob. 

3. Alice and Bob perform a *measurement* on their system. The outcome of this
   measurement yields their answers :math:`a \in \Gamma_A` and :math:`b \in
   \Gamma_B`. Specifically, Alice and Bob have collections of measurements

.. math::
    \begin{equation}
        \begin{aligned}
            \{ A_a^x : a \in \Gamma_{\text{A}} \} \subset \text{Pos}(\mathcal{U})
            \quad \text{and} \quad 
            \{ B_b^y : b \in \Gamma_{\text{B}} \} \subset \text{Pos}(\mathcal{V}),
        \end{aligned}
    \end{equation}

such that the measurements satisfy

.. math::
    \begin{equation}
        \begin{aligned}
            \sum_{a \in \Gamma_A} A_a^x = \mathbb{I}_{\mathcal{U}}
            \quad \text{and} \quad 
            \sum_{b \in \Gamma_B} B_b^y = \mathbb{I}_{\mathcal{V}}
        \end{aligned}
    \end{equation}

4. The referee determines whether Alice and Bob win or lose, based on the
   questions :math:`x` and :math:`y` as well as the answers :math:`a` and
   :math:`b`. 

For certain games, the probability that the players obtain a winning outcome is
higher if they use a quantum strategy as opposed to a classical one. This
striking separation is one primary motivation to study nonlocal games, as it
provides examples of tasks that benefit from the manipulation of quantum
information. 

Calculating the classical value
-------------------------------
(Coming soon)

Calculating the quantum value
------------------------------

The ability to calculate the quantum value for an arbitrary nonlocal game is a
highly non-trivial task. Indeed, the quantum value is only known in special
cases for certain nonlocal games.

For an arbitrary nonlocal game, there exist approaches that place upper and
lower bounds on the quantum value. The lower bound approach is calculated using
the technique of semidefinite programming [tLD07]_. While this method is efficient
to carry out, it does not guarantee convergence to the quantum value (although
in certain cases, it is attained).

The primary idea of this approach is to note that fixing the measurements on one
system yields the optimal measurements of the other system via an SDP. The
algorithm proceeds in an iterative manner between two SDPs. In the first SDP, we
assume that Bob's measurements are fixed, and Alice's measurements are to be
optimized over. In the second SDP, we take Alice's optimized measurements from
the first SDP and now optimize over Bob's measurements. This method is repeated
until the quantum value reaches a desired numerical precision.

For completeness, the first SDP where we fix Bob's measurements and optimize
over Alice's measurements is given as SDP-1.

.. math::

    \begin{equation}
        \begin{aligned}
            \textbf{SDP-1:} \quad & \\
            \text{maximize:} \quad & \sum_{(x,y \in \Sigma)} \pi(x,y)
                                     \sum_{(a,b) \in \Gamma}
                                     V(a,b|x,y)
                                     \langle B_b^y, A_a^x \rangle \\
            \text{subject to:} \quad & \sum_{a \in \Gamma_{\mathsf{A}}} =
                                       \tau, \qquad \qquad
                                       \forall x \in \Sigma_{\mathsf{A}}, \\
                               \quad & A_a^x \in \text{Pos}(\mathcal{A}),
                                       \qquad
                                       \forall x \in \Sigma_{\mathsf{A}}, \
                                       \forall a \in \Gamma_{\mathsf{A}}, \\
                                       & \tau \in \text{D}(\mathcal{A}).
        \end{aligned}
    \end{equation}

Similarly, the second SDP where we fix Alice's measurements and optimize over
Bob's measurements is given as SDP-2.

.. math::

    \begin{equation}
        \begin{aligned}
            \textbf{SDP-2:} \quad & \\
            \text{maximize:} \quad & \sum_{(x,y \in \Sigma)} \pi(x,y)
                                     \sum_{(a,b) \in \Gamma} V(a,b|x,y)
                                     \langle B_b^y, A_a^x \rangle \\
            \text{subject to:} \quad & \sum_{b \in \Gamma_{\mathsf{B}}} =
                                       \mathbb{I}, \qquad \qquad
                                       \forall y \in \Sigma_{\mathsf{B}}, \\
                               \quad & B_b^y \in \text{Pos}(\mathcal{B}),
                               \qquad \forall y \in \Sigma_{\mathsf{B}}, \
                               \forall b \in \Gamma_{\mathsf{B}}.
        \end{aligned}
    \end{equation}


Lower bounding the quantum value in `toqito`
---------------------------------------------

The :code:`toqito` software implements both of these optimization problems using
the :code:`cvxpy` library. We see-saw between the two SDPs until the value we
obtain reaches a specific precision threshold.

As we are not guaranteed to obtain the true quantum value of a given nonlocal
game as this approach can get stuck in a local minimum, the :code:`toqito`
function allows the user to specify an :code:`iters` argument that runs the
see-saw approach a number of times and then returns the highest of the values
obtained.

Example: Lower bounding the quantum value of the CHSH game
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let us consider calculating the lower bound on the quantum value of the CHSH
game.

.. note::
    As the CHSH game is a subtype of nonlocal game referred to as an XOR game,
    we do not necessarily need to resort to this lower bound technique as there
    exists a specific SDP formulation that one can use to directly compute the
    quantum value of an XOR game. More information on how one defines the CHSH
    game as well as this method to directly calculate the quantum value of an
    XOR game is provided in `"Calculating the Quantum and Classical Value of a Two-Player XOR Game" <https://toqito.readthedocs.io/en/latest/tutorials.xor_quantum_value.html>`_

We will use the CHSH game here as an illustrative example as we already know
what the optimal quantum value should be.

The first step is to use :code:`numpy` to encode a matrix that encapsulates the
probabilities with which the questions are asked to Alice and Bob. As defined in
the CHSH game, each of the four pairs
:math:`\{(0, 0), (0, 1), (1, 0), (1, 1)\}` are all equally likely. We encode
this in the matrix as follows.

.. code-block:: python

    >>> # Creating the probability matrix.
    >>> import numpy as np
    >>> prob_mat = np.array([[1 / 4, 1 / 4], [1 / 4, 1 / 4]])

Next, we want to loop through all possible combinations of question and answer
pairs and populate the :math:(a, b, x, y)^{th}` entry of that matrix with a
:math:`1` in the event that the winning condition is satisfied. Otherwise, if
the winning condition is not satisfied for that particular choice of
:math:`a, b, x,` and :math:`y`, we place a :math:`0` at that position.

The following code performs this operation and places the appropriate entries
in this matrix into the :code:`pred_mat` variable.

.. code-block:: python

    >>> # Creating the predicate matrix.
    >>> import numpy as np
    >>> num_alice_inputs, num_alice_outputs = 2, 2
    >>> num_bob_inputs, num_bob_outputs = 2, 2
    >>>
    >>> pred_mat = np.zeros(
    >>>     (num_alice_outputs, num_bob_outputs, num_alice_inputs, num_bob_inputs)
    >>> )
    >>>
    >>> for a_alice in range(num_alice_outputs):
    >>>     for b_bob in range(num_bob_outputs):
    >>>         for x_alice in range(num_alice_inputs):
    >>>             for y_bob in range(num_bob_inputs):
    >>>                 if a_alice ^ b_bob == x_alice * y_bob:
    >>>                     pred_mat[a_alice, b_bob, x_alice, y_bob] = 1
    >>> print(pred_mat)
    [[[[1. 1.]
       [1. 0.]]

      [[0. 0.]
       [0. 1.]]]


     [[[0. 0.]
       [0. 1.]]

      [[1. 1.]
       [1. 0.]]]]

Now that we have both :code:`prob_mat` and :code:`pred_mat` defined, we can
use :code:`toqito` to determine the lower bound on the quantum value.

.. code-block:: python

    >>> from toqito.nonlocal_games.nonlocal_game import NonlocalGame
    >>> chsh = NonlocalGame(prob_mat, pred_mat)
    >>> chsh.quantum_value_lower_bound()
    0.8535539268303678

In this case, we can see that the quantum value of the CHSH game is in fact
attained as :math:`\cos^2(\pi/8) \approx 0.85355`.

The FFL game
-------------

The *FFL (Fortnow, Feige, Lovasz) game* is a nonlocal game specified as
follows.

.. math::
    \begin{equation}
        \begin{aligned} 
            &\pi(0, 0) = \frac{1}{3}, \quad 
             \pi(0, 1) = \frac{1}{3}, \quad 
             \pi(1, 0) = \frac{1}{3}, \quad
             \pi(1, 1) = 0, \\ 
            &(x,y) \in \Sigma_A \times \Sigma_B, \qquad \text{and} \qquad (a, b) \in \Gamma_A \times \Gamma_B,
        \end{aligned}
    \end{equation}

where

    .. math::
        \begin{equation}
            \Sigma_A = \{0, 1\}, \quad \Sigma_B = \{0, 1\}, \quad \Gamma_A =
            \{0,1\}, \quad \text{and} \quad \Gamma_B = \{0, 1\}.
        \end{equation}

Alice and Bob win the FFL game if and only if the following equation is
satisfied

    .. math::
        \begin{equation}
        a \lor x = b \lor y.
        \end{equation}

It is well-known that both the classical and quantum value of this nonlocal
game is :math:`2/3` [tCHTW04]_. We can verify this fact using :code:`toqito`.
The following example encodes the FFL game. We then calculate the classical
value and calculate lower bounds on the quantum value of the FFL game.

.. code-block:: python

    >>> import numpy as np
    >>> from toqito.nonlocal_games.nonlocal_game import NonlocalGame
    >>>
    >>> # Specify the number of inputs, and number of outputs.
    >>> num_alice_in, num_alice_out = 2, 2
    >>> num_bob_in, num_bob_out = 2, 2
    >>> 
    >>> # Define the probability matrix of the FFL game.
    >>> prob_mat = np.array([[1/3, 1/3], [1/3, 0]])
    >>>
    >>>
    >>> # Define the predicate matrix of the FFL game.
    >>> pred_mat = np.zeros((num_alice_out, num_bob_out, num_alice_in, num_bob_in))
    >>> for a_alice in range(num_alice_out):
    >>>     for b_bob in range(num_bob_out):
    >>>         for x_alice in range(num_alice_in):
    >>>             for y_bob in range(num_bob_in):
    >>>                 if (a_alice or x_alice) != (b_bob or y_bob):
    >>>                     pred_mat[a_alice, b_bob, x_alice, y_bob] = 1
    >>> # Define the FFL game object.
    >>> ffl = NonlocalGame(prob_mat, pred_mat)
    >>> ffl.classical_value()
    0.6666666666666666
    >>> ffl.quantum_value_lower_bound()
    0.6666857549041076

In this case, we obtained the correct quantum value of :math:`2/3`, however,
the lower bound technique is not guaranteed to converge to the true quantum
value in general.

Parallel repetitions of nonlocal games
--------------------------------------
(Coming soon).

References
------------------------------

.. [tCHTW04] Cleve, Richard, Hoyer, Peter, Toner, Benjamin, and Watrous, John
    "Consequences and limits of nonlocal strategies"
    Computational Complexity 2004. Proceedings. 19th IEEE Annual Conference.
    https://arxiv.org/abs/quant-ph/0404076

.. [tLD07] Liang, Yeong-Cherng, and Andrew C. Doherty.
    "Bounds on quantum correlations in Bell-inequality experiments."
    Physical Review A 75.4 (2007): 042103.
    https://arxiv.org/abs/quant-ph/0608128
Nonlocal Games
=====================

A number of nonlocal game-related functions.

XOR Games
---------------------------------------------

.. automodule:: toqito.nonlocal_games.xor_game
   :members:
   :undoc-members:
   :show-inheritance:

Nonlocal Games
---------------

.. automodule:: toqito.nonlocal_games.nonlocal_game
   :members:
   :undoc-members:
   :show-inheritance:

Extended Nonlocal Games
------------------------

.. automodule:: toqito.nonlocal_games.extended_nonlocal_game
   :members:
   :undoc-members:
   :show-inheritance:

Quantum Hedging Protocols
-------------------------

.. automodule:: toqito.nonlocal_games.quantum_hedging
   :members:
   :undoc-members:
   :show-inheritance:
Superdense Coding
==================

Coming soon.Random Objects
==============

A number of functions for generating random objects.

Random
-------

.. toctree::

.. autosummary::
   :toctree: _autosummary

    toqito.random.random_density_matrix
    toqito.random.random_ginibre
    toqito.random.random_povm
    toqito.random.random_state_vector
    toqito.random.random_unitary
Tutorials
=========

Tutorials for :code:`toqito`.

Nonlocal Games
--------------------------

.. toctree::
   :maxdepth: 5

   tutorials.nonlocal_games

XOR Games
--------------------------

.. toctree::
   :maxdepth: 5

   tutorials.xor_quantum_value

Extended Nonlocal Games
--------------------------

.. toctree::
   :maxdepth: 5

   tutorials.extended_nonlocal_games

Quantum State Distinguishability
--------------------------------

.. toctree::
   :maxdepth: 5

   tutorials.state_distinguishability

Quantum State Exclusion
-----------------------

.. toctree::
   :maxdepth: 5

   tutorials.state_exclusion

Quantum Communication
-----------------------------------------------------

.. toctree::
   :maxdepth: 5

   tutorials.superdense_coding
Extended nonlocal games
==========================

In this tutorial, we will define the concept of an *extended nonlocal game*.
Extended nonlocal games are a more general abstraction of nonlocal games
wherein the referee, who previously only provided questions and answers to the
players, now share a state with the players and is able to perform a
measurement on that shared state. 

Every extended nonlocal game has a *value* associated to it. Analogously to
nonlocal games, this value is a quantity that dictates how well the players can
perform a task in the extended nonlocal game model when given access to certain
resources. We will be using :code:`toqito` to calculate these quantities.

We will also look at existing results in the literature on these values and be
able to replicate them using :code:`toqito`. Much of the written content in
this tutorial will be directly taken from [tRusso17]_.

Extended nonlocal games have a natural physical interpretation in the setting
of tripartite steering [tCSAN15]_ and in device-independent quantum scenarios [tTFKW13]_. For
more information on extended nonlocal games, please refer to [tJMRW16]_ and
[tRusso17]_.

The extended nonlocal game model
--------------------------------

An *extended nonlocal game* is similar to a nonlocal game in the sense that it
is a cooperative game played between two players Alice and Bob against a
referee. The game begins much like a nonlocal game, with the referee selecting
and sending a pair of questions :math:`(x,y)` according to a fixed probability
distribution. Once Alice and Bob receive :math:`x` and :math:`y`, they respond
with respective answers :math:`a` and :math:`b`. Unlike a nonlocal game, the
outcome of an extended nonlocal game is determined by measurements performed by
the referee on its share of the state initially provided to it by Alice and
Bob. 

.. figure:: figures/extended_nonlocal_game.svg
   :alt: extened nonlocal game
   :align: center

   An extended nonlocal game.

Specifically, Alice and Bob's winning probability is determined by
collections of measurements, :math:`V(a,b|x,y) \in \text{Pos}(\mathcal{R})`,
where :math:`\mathcal{R} = \mathbb{C}^m` is a complex Euclidean space with
:math:`m` denoting the dimension of the referee's quantum system--so if Alice
and Bob's response :math:`(a,b)` to the question pair :math:`(x,y)` leaves the
referee's system in the quantum state

.. math::
    \sigma_{a,b}^{x,y} \in \text{D}(\mathcal{R}),

then their winning and losing probabilities are given by

.. math::
    \left\langle V(a,b|x,y), \sigma_{a,b}^{x,y} \right\rangle 
    \quad \text{and} \quad 
    \left\langle \mathbb{I} - V(a,b|x,y), \sigma_{a,b}^{x,y} \right\rangle.

Strategies for extended nonlocal games
---------------------------------------

An extended nonlocal game :math:`G` is defined by a pair :math:`(\pi, V)`,
where :math:`\pi` is a probability distribution of the form

.. math::
    \pi : \Sigma_A \times \Sigma_B \rightarrow [0, 1]

on the Cartesian product of two alphabets :math:`\Sigma_A` and
:math:`\Sigma_B`, and :math:`V` is a function of the form

.. math::
    V : \Gamma_A \times \Gamma_B \times \Sigma_A \times \Sigma_B \rightarrow \text{Pos}(\mathcal{R})

for :math:`\Sigma_A` and :math:`\Sigma_B` as above, :math:`\Gamma_A` and
:math:`\Gamma_B` being alphabets, and :math:`\mathcal{R}` refers to the
referee's space. Just as in the case for nonlocal games, we shall use the
convention that

.. math::
    \Sigma = \Sigma_A \times \Sigma_B \quad \text{and} \quad \Gamma = \Gamma_A \times \Gamma_B

to denote the respective sets of questions asked to Alice and Bob and the sets
of answers sent from Alice and Bob to the referee.

When analyzing a strategy for Alice and Bob, it may be convenient to define a
function

.. math::
    K : \Gamma_A \times \Gamma_B \times \Sigma_A \times \Sigma_B \rightarrow \text{Pos}(\mathcal{R}).

We can represent Alice and Bob's winning probability for an extended nonlocal
game as

.. math::
    \sum_{(x,y) \in \Sigma} \pi(x,y) \sum_{(a,b) \in \Gamma} \left\langle V(a,b|x,y), K(a,b|x,y) \right\rangle.

Standard quantum strategies for extended nonlocal games
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A *standard quantum strategy* for an extended nonlocal game consists of
finite-dimensional complex Euclidean spaces :math:`\mathcal{U}` for Alice and
:math:`\mathcal{V}` for Bob, a quantum state :math:`\sigma \in
\text{D}(\mathcal{U} \otimes \mathcal{R} \otimes \mathcal{V})`, and two
collections of measurements

.. math::
    \{ A_a^x : a \in \Gamma_A \} \subset \text{Pos}(\mathcal{U})
    \quad \text{and} \quad
    \{ B_b^y : b \in \Gamma_B \} \subset \text{Pos}(\mathcal{V}),

for each :math:`x \in \Sigma_A` and :math:`y \in \Sigma_B` respectively. As
usual, the measurement operators satisfy the constraint that

.. math::
    \sum_{a \in \Gamma_A} A_a^x = \mathbb{I}_{\mathcal{U}} 
    \quad \text{and} \quad
    \sum_{b \in \Gamma_B} B_b^y = \mathbb{I}_{\mathcal{V}},

for each :math:`x \in \Sigma_A` and :math:`y \in \Sigma_B`.

When the game is played, Alice and Bob present the referee with a quantum
system so that the three parties share the state :math:`\sigma \in
\text{D}(\mathcal{U} \otimes \mathcal{R} \otimes \mathcal{V})`. The referee
selects questions :math:`(x,y) \in \Sigma` according to the distribution
:math:`\pi` that is known to all participants in the game.

The referee then sends :math:`x` to Alice and :math:`y` to Bob. At this point,
Alice and Bob make measurements on their respective portions of the state
:math:`\sigma` using their measurement operators to yield an outcome to send
back to the referee. Specifically, Alice measures her portion of the state
:math:`\sigma` with respect to her set of measurement operators :math:`\{A_a^x
: a \in \Gamma_A\}`, and sends the result :math:`a \in \Gamma_A` of this
measurement to the referee. Likewise, Bob measures his portion of the state
:math:`\sigma` with respect to his measurement operators 
:math:`\{B_b^y : b \in \Gamma_B\}` to yield the outcome :math:`b \in \Gamma_B`,
that is then sent back to the referee.

At the end of the protocol, the referee measures its quantum system with
respect to the measurement :math:`\{V(a,b|x,y), \mathbb{I}-V(a,b|x,y)\}`.

The winning probability for such a strategy in this game :math:`G = (\pi,V)` is
given by

.. math::
    \sum_{(x,y) \in \Sigma} \pi(x,y) \sum_{(a,b) \in \Gamma}
    \left \langle A_a^x \otimes V(a,b|x,y) \otimes B_b^y,
    \sigma
    \right \rangle.

For a given extended nonlocal game :math:`G = (\pi,V)`, we write
:math:`\omega^*(G)` to denote the *standard quantum value* of :math:`G`, which
is the supremum value of Alice and Bob's winning probability over all standard
quantum strategies for :math:`G`.

Unentangled strategies for extended nonlocal games
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An *unentangled strategy* for an extended nonlocal game is simply a standard
quantum strategy for which the state :math:`\sigma \in \text{D}(\mathcal{U}
\otimes \mathcal{R} \otimes \mathcal{V})` initially prepared by Alice and Bob
is fully separable.

Any unentangled strategy is equivalent to a strategy where Alice and Bob store
only classical information after the referee's quantum system has been provided
to it.

For a given extended nonlocal game :math:`G = (\pi, V)` we write
:math:`\omega(G)` to denote the *unentangled value* of :math:`G`, which is the
supremum value for Alice and Bob's winning probability in :math:`G` over all
unentangled strategies. The unentangled value of any extended nonlocal game,
:math:`G`, may be written as

.. math::
    \omega(G) = \max_{f, g}
    \lVert
    \sum_{(x,y) \in \Sigma} \pi(x,y)
    V(f(x), g(y)|x, y)
    \rVert

where the maximum is over all functions :math:`f : \Sigma_A \rightarrow
\Gamma_A` and :math:`g : \Sigma_B \rightarrow \Gamma_B`.

Non-signaling strategies for extended nonlocal games
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A *non-signaling strategy* for an extended nonlocal game consists of a function

.. math::
    K : \Gamma_A \times \Gamma_B \times \Sigma_A \times \Sigma_B \rightarrow \text{Pos}(\mathcal{R})

such that

.. math::
    \sum_{a \in \Gamma_A} K(a,b|x,y) = \rho_b^y \quad \text{and} \quad \sum_{b \in \Gamma_B} K(a,b|x,y) = \sigma_a^x,

for all :math:`x \in \Sigma_A` and :math:`y \in \Sigma_B` where
:math:`\{\rho_b^y : y \in \Sigma_B, b \in \Gamma_B\}` and :math:`\{\sigma_a^x:
x \in \Sigma_A, a \in \Gamma_A\}` are collections of operators satisfying

.. math::
    \sum_{a \in \Gamma_A} \sigma_a^x = \tau = \sum_{b \in \Gamma_B} \rho_b^y,

for every choice of :math:`x \in \Sigma_A` and :math:`y \in \Sigma_B` and where
:math:`\tau \in \text{D}(\mathcal{R})` is a density operator.

For any extended nonlocal game, :math:`G = (\pi, V)`, the winning probability
for a non-signaling strategy is given by

.. math::
    \sum_{(x,y) \in \Sigma} \pi(x,y) \sum_{(a,b) \in \Gamma} \left\langle V(a,b|x,y) K(a,b|x,y) \right\rangle.

We denote the *non-signaling value* of :math:`G` as :math:`\omega_{ns}(G)`
which is the supremum value of the winning probability of :math:`G` taken over
all non-signaling strategies for Alice and Bob.

Relationships between different strategies and values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For an extended nonlocal game, :math:`G`, the values have the following relationship:


.. note::
    .. math::
        0 \leq \omega(G) \leq \omega^*(G) \leq \omega_{ns}(G) \leq 1.

Example: The BB84 extended nonlocal game
-----------------------------------------

The *BB84 extended nonlocal game* is defined as follows. Let :math:`\Sigma_A =
\Sigma_B = \Gamma_A = \Gamma_B = \{0,1\}`, define

.. math::
    \begin{equation}
        \begin{aligned}
            V(0,0|0,0) = \begin{pmatrix}
                            1 & 0 \\
                            0 & 0
                         \end{pmatrix}, &\quad
            V(1,1|0,0) = \begin{pmatrix}
                            0 & 0 \\
                            0 & 1
                         \end{pmatrix}, \\
            V(0,0|1,1) = \frac{1}{2}\begin{pmatrix}
                            1 & 1 \\
                            1 & 1
                         \end{pmatrix}, &\quad
            V(1,1|1,1) = \frac{1}{2}\begin{pmatrix}
                            1 & -1 \\
                            -1 & 1
                         \end{pmatrix},
        \end{aligned}
    \end{equation}

define 

.. math::
    V(a,b|x,y) = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}

for all :math:`a \not= b` or :math:`x \not= y`, define :math:`\pi(0,0) =
\pi(1,1) = 1/2`, and define :math:`\pi(x,y) = 0` if :math:`x \not=y`.

We can encode the BB84 game, :math:`G_{BB84} = (\pi, V)`, in :code:`numpy`
arrays where :code:`prob_mat` corresponds to the probability distribution
:math:`\pi` and where :code:`pred_mat` corresponds to the operator :math:`V`. 

.. code-block:: python
    
    >>> """Define the BB84 extended nonlocal game."""
    >>> import numpy as np
    >>> from toqito.states import basis
    >>>
    >>> # The basis: {|0>, |1>}:
    >>> e_0, e_1 = basis(2, 0), basis(2, 1)
    >>>
    >>> # The basis: {|+>, |->}:
    >>> e_p = (e_0 + e_1) / np.sqrt(2)
    >>> e_m = (e_0 - e_1) / np.sqrt(2)
    >>>
    >>> # The dimension of referee's measurement operators:
    >>> dim = 2
    >>> # The number of outputs for Alice and Bob:
    >>> a_out, b_out = 2, 2
    >>> # The number of inputs for Alice and Bob:
    >>> a_in, b_in = 2, 2
    >>> 
    >>> # Define the predicate matrix V(a,b|x,y) \in Pos(R)
    >>> bb84_pred_mat = np.zeros([dim, dim, a_out, b_out, a_in, b_in])
    >>>
    >>> # V(0,0|0,0) = |0><0|
    >>> bb84_pred_mat[:, :, 0, 0, 0, 0] = e_0 * e_0.conj().T
    >>> # V(1,1|0,0) = |1><1|
    >>> bb84_pred_mat[:, :, 1, 1, 0, 0] = e_1 * e_1.conj().T
    >>> # V(0,0|1,1) = |+><+|
    >>> bb84_pred_mat[:, :, 0, 0, 1, 1] = e_p * e_p.conj().T
    >>> # V(1,1|1,1) = |-><-|
    >>> bb84_pred_mat[:, :, 1, 1, 1, 1] = e_m * e_m.conj().T
    >>>
    >>> # The probability matrix encode \pi(0,0) = \pi(1,1) = 1/2
    >>> bb84_prob_mat = 1/2*np.identity(2)

The unentangled value of the BB84 extended nonlocal game
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It was shown in [tTFKW13]_ and [tJMRW16]_ that

.. math::
    \omega(G_{BB84}) = \cos^2(\pi/8).

This can be verified in :code:`toqito` as follows.

.. code-block:: python

    >>> """Calculate the unentangled value of the BB84 extended nonlocal game."""
    >>> from toqito.nonlocal_games.extended_nonlocal_game import ExtendedNonlocalGame
    >>> 
    >>> # Define an ExtendedNonlocalGame object based on the BB84 game.
    >>> bb84 = ExtendedNonlocalGame(bb84_prob_mat, bb84_pred_mat)
    >>> 
    >>> # The unentangled value is cos(pi/8)**2 \approx 0.85356
    >>> bb84.unentangled_value() 
    0.8535533905544173

The BB84 game also exhibits strong parallel repetition. We can specify how many
parallel repetitions for :code:`toqito` to run. The example below provides an
example of two parallel repetitions for the BB84 game.

.. code-block:: python

    >>> """The unentangled value of BB84 under parallel repetition."""
    >>> from toqito.nonlocal_games.extended_nonlocal_game import ExtendedNonlocalGame
    >>> 
    >>> # Define the bb84 game for two parallel repetitions.
    >>> bb84_2_reps = ExtendedNonlocalGame(bb84_prob_mat, bb84_pred_mat, 2)
    >>> 
    >>> # The unentangled value for two parallel repetitions is cos(pi/8)**4 \approx 0.72855
    >>> bb84_2_reps.unentangled_value() 
    0.7285533940730632

It was shown in [tJMRW16]_ that the BB84 game possesses the property of strong
parallel repetition. That is,

.. math::
    \omega(G_{BB84}^r) = \omega(G_{BB84})^r

for any integer :math:`r`. 

The standard quantum value of the BB84 extended nonlocal game
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can calculate lower bounds on the standard quantum value of the BB84 game
using :code:`toqito` as well.

.. code-block:: python

    >>> """Calculate lower bounds on the standard quantum value of the BB84 extended nonlocal game."""
    >>> from toqito.nonlocal_games.extended_nonlocal_game import ExtendedNonlocalGame
    >>> 
    >>> # Define an ExtendedNonlocalGame object based on the BB84 game.
    >>> bb84_lb = ExtendedNonlocalGame(bb84_prob_mat, bb84_pred_mat)
    >>> 
    >>> # The standard quantum value is cos(pi/8)**2 \approx 0.85356
    >>> bb84_lb.quantum_value_lower_bound()
    0.8535533236834885

From [tJMRW16]_, it is known that :math:`\omega(G_{BB84}) =
\omega^*(G_{BB84})`, however, if we did not know this beforehand, we could
attempt to calculate upper bounds on the standard quantum value. 

There are a few methods to do this, but one easy way is to simply calculate the
non-signaling value of the game as this provides a natural upper bound on the
standard quantum value. Typically, this bound is not tight and usually not all
that useful in providing tight upper bounds on the standard quantum value,
however, in this case, it will prove to be useful.

The non-signaling value of the BB84 extended nonlocal game
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using :code:`toqito`, we can see that :math:`\omega_{ns}(G) = \cos^2(\pi/8)`.

.. code-block:: python

    >>> """Calculate the non-signaling value of the BB84 extended nonlocal game."""
    >>> from toqito.nonlocal_games.extended_nonlocal_game import ExtendedNonlocalGame
    >>> 
    >>> # Define an ExtendedNonlocalGame object based on the BB84 game.
    >>> bb84 = ExtendedNonlocalGame(bb84_prob_mat, bb84_pred_mat)
    >>> 
    >>> # The non-signaling value is cos(pi/8)**2 \approx 0.85356
    >>> bb84.nonsignaling_value() 
    0.853486975032519

So we have the relationship that

.. math::
    \omega(G_{BB84}) = \omega^*(G_{BB84}) = \omega_{ns}(G_{BB84}) = \cos^2(\pi/8).

It turns out that strong parallel repetition does *not* hold in the
non-signaling scenario for the BB84 game. This was shown in [tRusso17]_, and we
can observe this by the following snippet.

.. code-block:: python

    >>> """The non-signaling value of BB84 under parallel repetition."""
    >>> from toqito.nonlocal_games.extended_nonlocal_game import ExtendedNonlocalGame
    >>> 
    >>> # Define the bb84 game for two parallel repetitions.
    >>> bb84_2_reps = ExtendedNonlocalGame(bb84_prob_mat, bb84_pred_mat, 2)
    >>> 
    >>> # The non-signaling value for two parallel repetitions is cos(pi/8)**4 \approx 0.73825
    >>> bb84_2_reps.nonsignaling_value() 
    0.7382545498689419

Note that :math:`0.73825 \geq \cos(\pi/8)^4 \approx 0.72855` and therefore we
have that

.. math::
    \omega_{ns}(G^r_{BB84}) \not= \omega_{ns}(G_{BB84})^r

for :math:`r = 2`.

Example: The CHSH extended nonlocal game
-----------------------------------------

Let us now define another extended nonlocal game, :math:`G_{CHSH}`.

Let :math:`\Sigma_A = \Sigma_B = \Gamma_A = \Gamma_B = \{0,1\}`, define a
collection of measurements :math:`\{V(a,b|x,y) : a \in \Gamma_A, b \in
\Gamma_B, x \in \Sigma_A, y \in \Sigma_B\} \subset \text{Pos}(\mathcal{R})`
such that

.. math::
    \begin{equation}
        \begin{aligned}
            V(0,0|0,0) = V(0,0|0,1) = V(0,0|1,0) = \begin{pmatrix}
                                                    1 & 0 \\
                                                    0 & 0
                                                   \end{pmatrix}, \\
            V(1,1|0,0) = V(1,1|0,1) = V(1,1|1,0) = \begin{pmatrix}
                                                    0 & 0 \\
                                                    0 & 1
                                                   \end{pmatrix}, \\
            V(0,1|1,1) = \frac{1}{2}\begin{pmatrix}
                                        1 & 1 \\
                                        1 & 1
                                    \end{pmatrix}, \\
            V(1,0|1,1) = \frac{1}{2} \begin{pmatrix}
                                        1 & -1 \\
                                        -1 & 1
                                     \end{pmatrix},
        \end{aligned}
    \end{equation}

define 

.. math::
    V(a,b|x,y) = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}

for all :math:`a \oplus b \not= x \land y`, and define :math:`\pi(0,0) =
\pi(0,1) = \pi(1,0) = \pi(1,1) = 1/4`.

In the event that :math:`a \oplus b \not= x \land y`, the referee's measurement
corresponds to the zero matrix. If instead it happens that :math:`a \oplus b =
x \land y`, the referee then proceeds to measure with respect to one of the
measurement operators. This winning condition is reminiscent of the standard
CHSH nonlocal game.

We can encode :math:`G_{CHSH}` in a similar way using :code:`numpy` arrays as
we did for :math:`G_{BB84}`.

.. code-block:: python

    >>> """Define the CHSH extended nonlocal game."""
    >>> import numpy as np
    >>>
    >>> # The dimension of referee's measurement operators:
    >>> dim = 2
    >>> # The number of outputs for Alice and Bob:
    >>> a_out, b_out = 2, 2
    >>> # The number of inputs for Alice and Bob:
    >>> a_in, b_in = 2, 2
    >>> 
    >>> # Define the predicate matrix V(a,b|x,y) \in Pos(R)
    >>> chsh_pred_mat = np.zeros([dim, dim, a_out, b_out, a_in, b_in])
    >>>
    >>> # V(0,0|0,0) = V(0,0|0,1) = V(0,0|1,0).
    >>> chsh_pred_mat[:, :, 0, 0, 0, 0] = np.array([[1, 0], [0, 0]])
    >>> chsh_pred_mat[:, :, 0, 0, 0, 1] = np.array([[1, 0], [0, 0]])
    >>> chsh_pred_mat[:, :, 0, 0, 1, 0] = np.array([[1, 0], [0, 0]])
    >>>
    >>> # V(1,1|0,0) = V(1,1|0,1) = V(1,1|1,0).
    >>> chsh_pred_mat[:, :, 1, 1, 0, 0] = np.array([[0, 0], [0, 1]])
    >>> chsh_pred_mat[:, :, 1, 1, 0, 1] = np.array([[0, 0], [0, 1]])
    >>> chsh_pred_mat[:, :, 1, 1, 1, 0] = np.array([[0, 0], [0, 1]])
    >>>
    >>> # V(0,1|1,1)
    >>> chsh_pred_mat[:, :, 0, 1, 1, 1] = 1/2 * np.array([[1, 1], [1, 1]])
    >>>
    >>> # V(1,0|1,1)
    >>> chsh_pred_mat[:, :, 1, 0, 1, 1] = 1/2 * np.array([[1, -1], [-1, 1]])
    >>>
    >>> # The probability matrix encode \pi(0,0) = \pi(0,1) = \pi(1,0) = \pi(1,1) = 1/4.
    >>> chsh_prob_mat = np.array([[1/4, 1/4], [1/4, 1/4]])


Example: The unentangled value of the CHSH extended nonlocal game
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Similar to what we did for the BB84 extended nonlocal game, we can also compute
the unentangled value of :math:`G_{CHSH}`.

.. code-block:: python

    >>> """Calculate the unentangled value of the CHSH extended nonlocal game."""
    >>> from toqito.nonlocal_games.extended_nonlocal_game import ExtendedNonlocalGame
    >>> 
    >>> # Define an ExtendedNonlocalGame object based on the CHSH game.
    >>> chsh = ExtendedNonlocalGame(chsh_prob_mat, chsh_pred_mat)
    >>> 
    >>> # The unentangled value is 3/4 = 0.75
    >>> chsh.unentangled_value() 
    0.7499999999992315

We can also run multiple repetitions of :math:`G_{CHSH}`.

.. code-block:: python

    >>> """The unentangled value of CHSH under parallel repetition."""
    >>> from toqito.nonlocal_games.extended_nonlocal_game import ExtendedNonlocalGame
    >>> 
    >>> # Define the CHSH game for two parallel repetitions.
    >>> chsh_2_reps = ExtendedNonlocalGame(chsh_prob_mat, chsh_pred_mat, 2)
    >>> 
    >>> # The unentangled value for two parallel repetitions is (3/4)**2 \approx 0.5625
    >>> chsh_2_reps.unentangled_value() 
    0.5625000000002018

Note that strong parallel repetition holds as

.. math::
    \omega(G_{CHSH})^2 = \omega(G_{CHSH}^2) = \left(\frac{3}{4}\right)^2.

Example: The non-signaling value of the CHSH extended nonlocal game
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To obtain an upper bound for :math:`G_{CHSH}`, we can calculate the
non-signaling value.

.. code-block:: python

    >>> """Calculate the non-signaling value of the CHSH extended nonlocal game."""
    >>> from toqito.nonlocal_games.extended_nonlocal_game import ExtendedNonlocalGame
    >>> 
    >>> # Define an ExtendedNonlocalGame object based on the CHSH game.
    >>> chsh = ExtendedNonlocalGame(chsh_prob_mat, chsh_pred_mat)
    >>> 
    >>> # The non-signaling value is 3/4 = 0.75
    >>> chsh.nonsignaling_value() 
    0.7500002249607216

As we know that :math:`\omega(G_{CHSH}) = \omega_{ns}(G_{CHSH}) = 3/4` and that

.. math::
    \omega(G) \leq \omega^*(G) \leq \omega_{ns}(G)

for any extended nonlocal game, :math:`G`, we may also conclude that
:math:`\omega^*(G) = 3/4`.

Note the SCS convex optimization solver will generate a large number of warnings of the form

```
WARN: A->p (column pointers) not strictly increasing
```

This is a known issue, and while it does not appear to impact the correctness
of the results, it is an outstanding issue for the :code:`toqito` project.

Example: An extended nonlocal game with quantum advantage
----------------------------------------------------------

So far, we have only seen examples of extended nonlocal games where the
standard quantum and unentangled values are equal. Here we'll see an example of
an extended nonlocal game where the standard quantum value is *strictly higher*
than the unentangled value.


Example: A monogamy-of-entanglement game with mutually unbiased bases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let :math:`\zeta = \exp(\frac{2 \pi i}{3})` and consider the following four
mutually unbiased bases:

.. math::
    \begin{equation}\label{eq:MUB43}
    \begin{aligned}
      \mathcal{B}_0 &= \left\{ e_0,\: e_1,\: e_2 \right\}, \\
      \mathcal{B}_1 &= \left\{ \frac{e_0 + e_1 + e_2}{\sqrt{3}},\:
      \frac{e_0 + \zeta^2 e_1 + \zeta e_2}{\sqrt{3}},\:
      \frac{e_0 + \zeta e_1 + \zeta^2 e_2}{\sqrt{3}} \right\}, \\
      \mathcal{B}_2 &= \left\{ \frac{e_0 + e_1 + \zeta e_2}{\sqrt{3}},\:
      \frac{e_0 + \zeta^2 e_1 + \zeta^2 e_2}{\sqrt{3}},\:
      \frac{e_0 + \zeta e_1 + e_2}{\sqrt{3}} \right\}, \\
      \mathcal{B}_3 &= \left\{ \frac{e_0 + e_1 + \zeta^2 e_2}{\sqrt{3}},\:
      \frac{e_0 + \zeta^2 e_1 + e_2}{\sqrt{3}},\:
      \frac{e_0 + \zeta e_1 + \zeta e_2}{\sqrt{3}} \right\}.
    \end{aligned}
    \end{equation} 

Define an extended nonlocal game :math:`G_{MUB} = (\pi,R)` so that

.. math::

 		\pi(0) = \pi(1) = \pi(2) = \pi(3) = \frac{1}{4}

and :math:`R` is such that

.. math::
 		{ R(0|x), R(1|x), R(2|x) }

represents a measurement with respect to the basis :math:`\mathcal{B}_x`, for
each :math:`x \in \{0,1,2,3\}`.

Taking the description of :math:`G_{MUB}`, we can encode this as follows.

.. code-block:: python

    >>> """Define the monogamy-of-entanglement game defined by MUBs."""
    >>>  prob_mat = 1 / 4 * np.identity(4)
    >>>
    >>>  dim = 3
    >>>  e_0, e_1, e_2 = basis(dim, 0), basis(dim, 1), basis(dim, 2)
    >>>
    >>>  eta = np.exp((2 * np.pi * 1j) / dim)
    >>>  mub_0 = [e_0, e_1, e_2]
    >>>  mub_1 = [
    >>>      (e_0 + e_1 + e_2) / np.sqrt(3),
    >>>      (e_0 + eta ** 2 * e_1 + eta * e_2) / np.sqrt(3),
    >>>      (e_0 + eta * e_1 + eta ** 2 * e_2) / np.sqrt(3),
    >>>  ]
    >>>  mub_2 = [
    >>>      (e_0 + e_1 + eta * e_2) / np.sqrt(3),
    >>>      (e_0 + eta ** 2 * e_1 + eta ** 2 * e_2) / np.sqrt(3),
    >>>      (e_0 + eta * e_1 + e_2) / np.sqrt(3),
    >>>  ]
    >>>  mub_3 = [
    >>>      (e_0 + e_1 + eta ** 2 * e_2) / np.sqrt(3),
    >>>      (e_0 + eta ** 2 * e_1 + e_2) / np.sqrt(3),
    >>>      (e_0 + eta * e_1 + eta * e_2) / np.sqrt(3),
    >>>  ]
    >>>
    >>>  # List of measurements defined from mutually unbiased basis.
    >>>  mubs = [mub_0, mub_1, mub_2, mub_3]
    >>> 
    >>>  num_in = 4
    >>>  num_out = 3
    >>>  pred_mat = np.zeros([dim, dim, num_out, num_out, num_in, num_in], dtype=complex)
    >>>
    >>>  pred_mat[:, :, 0, 0, 0, 0] = mubs[0][0] * mubs[0][0].conj().T
    >>>  pred_mat[:, :, 1, 1, 0, 0] = mubs[0][1] * mubs[0][1].conj().T
    >>>  pred_mat[:, :, 2, 2, 0, 0] = mubs[0][2] * mubs[0][2].conj().T
    >>>
    >>>  pred_mat[:, :, 0, 0, 1, 1] = mubs[1][0] * mubs[1][0].conj().T
    >>>  pred_mat[:, :, 1, 1, 1, 1] = mubs[1][1] * mubs[1][1].conj().T
    >>>  pred_mat[:, :, 2, 2, 1, 1] = mubs[1][2] * mubs[1][2].conj().T
    >>>
    >>>  pred_mat[:, :, 0, 0, 2, 2] = mubs[2][0] * mubs[2][0].conj().T
    >>>  pred_mat[:, :, 1, 1, 2, 2] = mubs[2][1] * mubs[2][1].conj().T
    >>>  pred_mat[:, :, 2, 2, 2, 2] = mubs[2][2] * mubs[2][2].conj().T
    >>>
    >>>  pred_mat[:, :, 0, 0, 3, 3] = mubs[3][0] * mubs[3][0].conj().T
    >>>  pred_mat[:, :, 1, 1, 3, 3] = mubs[3][1] * mubs[3][1].conj().T
    >>>  pred_mat[:, :, 2, 2, 3, 3] = mubs[3][2] * mubs[3][2].conj().T

Now that we have encoded :math:`G_{MUB}`, we can calculate the unentangled value.

.. code-block:: python

    >>> g_mub = ExtendedNonlocalGame(prob_mat, pred_mat)
    >>> unent_val = g_mub.unentangled_value()
    >>> unent_val
    0.6545084973280103

That is, we have that 

.. math::

    \omega(G_{MUB}) = \frac{3 + \sqrt{5}}{8} \approx 0.65409.

However, if we attempt to run a lower bound on the standard quantum value, we
obtain.

.. code-block:: python

    >>> q_val = g_mub.quantum_value()
    >>> q_val
    0.660931321341278

Note that as we are calculating a lower bound, it is possible that a value this
high will not be obtained, or in other words, the algorithm can get stuck in a
local maximum that prevents it from finding the global maximum.

It is uncertain what the optimal standard quantum strategy is for this game,
but the value of such a strategy is bounded as follows

.. math::

    2/3 \geq \omega^*(G) \geq 0.6609.

For further information on the :math:`G_{MUB}` game, consult [tRusso17]_.

References
------------------------------

.. [tJMRW16] Johnston, Nathaniel, Mittal, Rajat, Russo, Vincent, Watrous, John
    "Extended non-local games and monogamy-of-entanglement games"
    Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences 472.2189 (2016),
    https://arxiv.org/abs/1510.02083

.. [tCSAN15] Cavalcanti, Daniel, Skrzypczyk, Paul, Aguilar, Gregory, Nery, Ranieri
    "Detection of entanglement in asymmetric quantum networks and multipartite quantum steering"
    Nature Communications, 6(7941), 2015
    https://arxiv.org/abs/1412.7730

.. [tTFKW13] Tomamichel, Marco, Fehr, Serge, Kaniewski, Jkedrzej, and Wehner, Stephanie.
    "A Monogamy-of-Entanglement Game With Applications to Device-Independent Quantum Cryptography"
    New Journal of Physics 15.10 (2013): 103002,
    https://arxiv.org/abs/1210.4359

.. [tRusso17] Russo, Vincent
    "Extended nonlocal games"
    https://arxiv.org/abs/1704.07375

Quantum state exclusion
=======================

In this tutorial, we are going to cover the problem of *quantum state
exclusion*. We are going to briefly describe the problem setting and then
describe how one may use :code:`toqito` to calculate the optimal probability
with which this problem can be solved for a number of different scenarios.

Quantum state exclusion is very closely related to the problem of quantum state
distinguishability. It may be useful to consult the following tutorial that
covers quantum state distinguishability:

* `Quantum State Distinguishability <https://toqito.readthedocs.io/en/latest/tutorials.state_distinguishability.html>`_

Further information beyond the scope of this tutorial can be found in the text
[tPBR12]_ as well as the course [tBJOP14]_.


The state exclusion problem
---------------------------

The quantum state exclusion problem is phrased as follows.

1. Alice possesses an ensemble of :math:`n` quantum states:

    .. math::
        \begin{equation}
            \eta = \left( (p_0, \rho_0), \ldots, (p_n, \rho_n)  \right),
        \end{equation}

where :math:`p_i` is the probability with which state :math:`\rho_i` is
selected from the ensemble. Alice picks :math:`\rho_i` with probability
:math:`p_i` from her ensemble and sends :math:`\rho_i` to Bob.

2. Bob receives :math:`\rho_i`. Both Alice and Bob are aware of how the
   ensemble is defined but he does *not* know what index :math:`i`
   corresponding to the state :math:`\rho_i` he receives from Alice is.

3. Bob wants to guess which of the states from the ensemble he was *not* given.
   In order to do so, he may measure :math:`\rho_i` to guess the index
   :math:`i` for which the state in the ensemble corresponds.

This setting is depicted in the following figure.

.. figure:: figures/quantum_state_distinguish.svg
   :alt: quantum state exclusion
   :align: center

   The quantum state exclusion setting.

.. note::
    The primary difference between the quantum state distinguishability
    scenario and the quantum state exclusion scenario is that in the former,
    Bob want to guess which state he was given, and in the latter, Bob wants to
    guess which state he was *not* given.

Optimal probability of conclusively excluding a quantum state
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



Optimal probability of unambiguously excluding a quantum state
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

References
------------------------------

.. [tPBR12] Pusey, Matthew, Barret, Jonathan, and Rudolph, Terry
    "On the reality of the quantum state"
    Nature Physics 8.6 (2012): 475-478.
    arXiv:1111.3328

.. [tBJOP14] Bandyopadhyay, Somshubhro, Jain, Rahul, Oppenheim, Jonathan, Perry, Christopher
    "Conclusive exclusion of quantum states"
    Physical Review A 89.2 (2014): 022336.
    arXiv:1306.4683

Linear Algebra
==============

A number of linear algebra-related functions.

Matrices
-----------

.. toctree::

.. autosummary::
   :toctree: _autosummary

    toqito.matrices.clock
    toqito.matrices.cnot
    toqito.matrices.fourier
    toqito.matrices.gell_mann
    toqito.matrices.gen_gell_mann
    toqito.matrices.gen_pauli
    toqito.matrices.hadamard
    toqito.matrices.iden
    toqito.matrices.pauli
    toqito.matrices.shift

Operations on Matrices and Vectors
----------------------------------

.. toctree::

.. autosummary::
   :toctree: _autosummary

    toqito.matrix_ops.tensor
    toqito.matrix_ops.unvec
    toqito.matrix_ops.vec

Properties of Matrices and Vectors
----------------------------------

.. toctree::

.. autosummary::
   :toctree: _autosummary

    toqito.matrix_props.is_commuting
    toqito.matrix_props.is_density
    toqito.matrix_props.is_diagonal
    toqito.matrix_props.is_hermitian
    toqito.matrix_props.is_identity
    toqito.matrix_props.is_idempotent
    toqito.matrix_props.is_normal
    toqito.matrix_props.is_permutation
    toqito.matrix_props.is_positive_definite
    toqito.matrix_props.is_projection
    toqito.matrix_props.is_positive_semidefinite
    toqito.matrix_props.is_square
    toqito.matrix_props.is_symmetric
    toqito.matrix_props.is_unitary
    toqito.matrix_props.majorizes
Getting started
===============

Installing
^^^^^^^^^^

1. Ensure you have Python 3.7 or greater installed on your machine.

2. Consider using a `virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtualenv/>`_.


3. The preferred way to install the :code:`toqito` package is via :code:`pip`.

.. code-block:: bash

    pip install toqito

Alternatively, to install, you may also run the following command from the
top-level package directory.

.. code-block:: bash

    python setup.py install

Installing BLAS/LAPACK
^^^^^^^^^^^^^^^^^^^^^^

The :code:`toqito` module makes heavy use of the :code:`cvxpy` module for solving various
convex optimization problems that naturally arise for certain problems in
quantum information. The installation instructions for :code:`cvxpy` may be found on
the project's `installation page <https://www.cvxpy.org/install/index.html>`_.

As a dependency for many of the solvers, you will need to ensure you have the
BLAS and LAPACK mathematical libraries installed on your machine. If you have
:code:`numpy` working on your machine, it is likely that you already have these
libraries on your machine. The :code:`cvxpy` module provides many different solvers
to select from for solving SDPs. We tend to use the
`SCS <https://github.com/cvxgrp/scs>`_ solver. Ensure that you have the :code:`scs`
Python module installed and built for your machine.

Testing
^^^^^^^

The :code:`pytest` module is used for testing. In order to run and :code:`pytest`, you will need to ensure it is
installed on your machine. Consult the `pytest <https://docs.pytest.org/en/latest/>`_ website for more information. To
run the suite of tests for :code:`toqito`, run the following command in the root directory of this project:

.. code-block:: bash

    pytest --cov-report term-missing --cov=toqito tests/

Contributing
^^^^^^^^^^^^

All contributions, bug reports, bug fixes, documentation improvements,
enhancements, and ideas are welcome.

A detailed overview on how to contribute can be found in the
`contributing guide <https://github.com/vprusso/toqito/blob/master/.github/CONTRIBUTING.md>`_.

Reporting Issues
^^^^^^^^^^^^^^^^

Please report any issues you encounter via the
`issue template <https://github.com/vprusso/toqito/blob/master/.github/ISSUE_TEMPLATE.md>`_.

Citing
^^^^^^

You can cite :code:`toqito` using the following DOI: `10.5281/zenodo.4743211 <https://zenodo.org/record/4743211>`_.

If you are using the :code:`toqito` software package in research work, please
include an explicit mention of :code:`toqito` in your publication. Something
along the lines of:

    To solve problem "X" we used `toqito`; a package for studying certain
    aspects of quantum information.

A BibTeX entry that you can use to cite :code:`toqito` is provided here:

.. code-block:: bash

    @misc{toqito,
       author       = {Vincent Russo},
       title        = {toqito: A {P}ython toolkit for quantum information, version 1.0.0},
       howpublished = {\url{https://github.com/vprusso/toqito}},
       month        = Mar,
       year         = 2021,
       doi          = {10.5281/zenodo.4743211}
     }
Measurements
=====================

A *measurement* can be defined as a function

.. math::
    \mu : \Sigma \rightarrow \text{Pos}(\mathcal{X})

satisfying

.. math::
    \sum_{a \in \Sigma} \mu(a) = \mathbb{I}_{\mathcal{X}}

where :math:`\Sigma` represents a set of measurement outcomes and
where :math:`\mu(a)` represents the measurement operator associated
with outcome :math:`a \in \Sigma`.

Operations on Measurements
--------------------------

.. toctree::

.. autosummary::
   :toctree: _autosummary

    toqito.measurement_ops.measure


Properties of Measurements
--------------------------

.. toctree::

.. autosummary::
   :toctree: _autosummary

    toqito.measurement_props.is_povm
.. figure:: figures/logo.svg
   :alt: toqito logo
   :align: center

toqito
======

The :code:`toqito` package is an open source library for studying various
objects in quantum information, namely, states, channels, and measurements.
toqito focuses on providing numerical tools to study problems pertaining to
entanglement theory, nonlocal games, and other aspects of quantum information
that are often associated with computer science.

.. image:: http://img.shields.io/travis/vprusso/toqito.svg?style=plastic
   :alt: Build Status
   :target: https://travis-ci.org/vprusso/toqito
.. image:: https://readthedocs.org/projects/toqito/badge/?version=latest&style=plastic
   :alt: Documentation
   :target: https://toqito.readthedocs.io/en/latest/
.. image:: https://codecov.io/gh/vprusso/toqito/branch/master/graph/badge.svg?style=plastic
   :alt: Codecov
   :target: https://codecov.io/gh/vprusso/toqito
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4743211.svg?style=plastic
   :alt: DOI
   :target: https://doi.org/10.5281/zenodo.4743211
.. image:: https://img.shields.io/badge/Supported%20By-UNITARY%20FUND-brightgreen.svg?style=plastic
   :alt: Unitary Fund
   :target: http://unitary.fund


User Documentation
------------------

.. toctree::
    :maxdepth: 2

    getting_started
    intro_tutorial

.. toctree::
    :maxdepth: 2
    :caption: Fundamental Objects

    states
    channels
    measurements

.. toctree::
    :maxdepth: 2
    :caption: Mathematics

    matrices
    perms

.. toctree::
    :maxdepth: 2
    :caption: Random Objects

    random

.. toctree::
    :maxdepth: 2
    :caption: Nonlocal Games

    nonlocal_games


.. toctree::
    :maxdepth: 2
    :caption: Tutorials

    tutorials




Permutations
============

A number of permutations-related functions.

Permutations
-----------------------------------

.. toctree::

.. autosummary::
   :toctree: _autosummary

    toqito.perms.antisymmetric_projection
    toqito.perms.perfect_matchings
    toqito.perms.perm_sign
    toqito.perms.permutation_operator
    toqito.perms.permute_systems
    toqito.perms.swap
    toqito.perms.swap_operator
    toqito.perms.symmetric_projection
    toqito.perms.unique_perms
Channels
=========

A *quantum channel* can be defined as a completely positive and trace preserving
linear map.

More formally, let :math:`\mathcal{X}` and :math:`\mathcal{Y}` represent complex
Euclidean spaces and let :math:`\text{L}(\cdot)` represent the set of linear
operators. Then a quantum channel, :math:`\Phi` is defined as

.. math::
    \Phi : \text{L}(\mathcal{X}) \rightarrow \text{L}(\mathcal{Y})

such that :math:`\Phi` is completely positive and trace preserving.

Distance Metrics for Quantum Channels
-------------------------------------

.. toctree::

.. autosummary::
   :toctree: _autosummary

    toqito.channel_metrics.channel_fidelity

Quantum Channels
----------------

.. toctree::

.. autosummary::
   :toctree: _autosummary

    toqito.channels.choi
    toqito.channels.dephasing
    toqito.channels.depolarizing
    toqito.channels.partial_trace
    toqito.channels.partial_transpose
    toqito.channels.realignment
    toqito.channels.reduction

Operations on Quantum Channels
------------------------------

.. toctree::

.. autosummary::
   :toctree: _autosummary

    toqito.channel_ops.apply_channel
    toqito.channel_ops.choi_to_kraus
    toqito.channel_ops.kraus_to_choi
    toqito.channel_ops.partial_channel
    toqito.channel_ops.dual_channel

Properties of Quantum Channels
------------------------------

.. toctree::

.. autosummary::
   :toctree: _autosummary

    toqito.channel_props.is_completely_positive
    toqito.channel_props.is_herm_preserving
    toqito.channel_props.is_positive
    toqito.channel_props.is_trace_preserving
    toqito.channel_props.is_unital
    toqito.channel_props.choi_rank
    toqito.channel_props.is_quantum_channel
    toqito.channel_props.is_unitary
Calculating the quantum and classical value of a two-player XOR game
=====================================================================

In this tutorial, we will cover the concept of an *XOR game*. We will also
showcase how the :code:`toqito` software package can be used to calculate the
classical and quantum value of a given XOR game.

For readers who are already familiar with XOR games and who simply want to see
how to use :code:`toqito` to study these objects, they are welcome to consult the
documentation page, and more specifically the function `xor\_game\_value
<https://toqito.readthedocs.io/en/latest/nonlocal_games.xor_games.html>`_.

Further information beyond the scope of this tutorial on the notion of XOR
games along with the method of computing their quantum value may be found in
[tCSUU08]_.

Two-player XOR games
--------------------

A two-player XOR game is a nonlocal game in which the winning condition is
predicated on an XOR function. For more information on the more general class
of nonlocal games along with how one defines classical and quantum strategies
for these games, please refer to the example:

* `Lower Bounds on the Quantum Value of a Two-Player Nonlocal Game
  <https://toqito.readthedocs.io/en/latest/tutorials.nonlocal_quantum_lower_bound.html>`_

.. note::
    It is *not* known how to directly compute the quantum value of an arbitrary
    nonlocal game. For the subset of XOR games, it turns out that it is
    possible to directly calculate the quantum value by solving a semidefinite
    program. The :code:`toqito` package obtains the quantum value of an XOR game
    in this manner.

The rest of this tutorial is concerned with analyzing specific XOR games.

The CHSH game
-------------

The *CHSH game* is a two-player XOR game with the following probability
distribution and question and answer sets.

.. math::
    \begin{equation}
        \begin{aligned} \pi(x,y) = \frac{1}{4}, \qquad (x,y) \in
                        \Sigma_A \times
            \Sigma_B, \qquad \text{and} \qquad (a, b) \in \Gamma_A \times
            \Gamma_B,
        \end{aligned}
    \end{equation}

where

.. math::
    \begin{equation}
        \Sigma_A = \{0, 1\}, \quad \Sigma_B = \{0, 1\}, \quad \Gamma_A =
        \{0,1\}, \quad \text{and} \quad \Gamma_B = \{0, 1\}.
    \end{equation}

Alice and Bob win the CHSH game if and only if the following equation is
satisfied

.. math::
    \begin{equation}
        a \oplus b = x y.
    \end{equation}

Recall that :math:`\oplus` refers to the XOR operation. 

For each question scenario, the following table provides what the winning
condition must be equal to for each question tuple to induce a winning outcome.

.. table::
    :align: center

    +-------------+-------------+----------------------+
    | :math:`x`   | :math:`y`   |  :math:`a \oplus b`  |
    +=============+=============+======================+
    | :math:`0`   | :math:`0`   | :math:`0`            |
    +-------------+-------------+----------------------+
    | :math:`0`   | :math:`1`   | :math:`0`            |
    +-------------+-------------+----------------------+
    | :math:`1`   | :math:`0`   | :math:`0`            |
    +-------------+-------------+----------------------+
    | :math:`1`   | :math:`1`   | :math:`1`            |
    +-------------+-------------+----------------------+

In order to specify an XOR game in :code:`toqito`, we will define two matrices:

    * :code:`prob_mat`: A matrix whose :math:`(x, y)^{th}` entry corresponds to
      the probability that Alice receives question :math:`x` and Bob receives
      question :math:`y`.

    * :code:`pred_mat`: A matrix whose :math:`(x, y)^{th}` entry corresponds to
      the winning choice of :math:`a` and :math:`b` when Alice receives
      :math:`x` and Bob receives :math:`y` from the referee.

For the CHSH game, the `prob_mat` and `pred_mat` variables are defined as follows.

.. code-block:: python

    >>> import numpy as np
    >>> prob_mat = np.array([[1/4, 1/4],
    >>>                      [1/4, 1/4]])
    >>> pred_mat = np.array([[0, 0],
    >>>                      [0, 1]])

That is, the :code:`prob_mat` matrix encapsulates that each question pair
:math:`\{(0,0), (0, 1), (1, 0), (1, 1)\}` is equally likely. 

The :code:`pred_mat` matrix indicates what the winning outcome of Alice and Bob
should be. For instance, :code:`pred_mat[0][0] = 0` describes the scenario where
Alice and Bob both receive :math:`0` as input. As we want to satisfy the
winning condition :math:`x \land y = a \oplus b`, we must have that :math:`a
\oplus b = 0` to satisfy the case when both :math:`x` and :math:`y` are equal
to zero. A similar logic can be followed to populate the remaining entries of
the :code:`pred_mat` variable.

We will use both of the :code:`prob_mat` and :code:`pred_mat` variables in the
coming subsections to make use of the :code:`toqito` package to compute both the
classical and quantum value of the CHSH game.

A classical strategy for the CHSH game
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can begin by asking; is it possible for Alice and Bob to win for every
single question pair they receive with certainty? If Alice and Bob use a
classical strategy, the answer to this question is "no". To see why, consider
the following equations:

.. math::
    \begin{equation}
        \begin{aligned}
            a_0 \oplus b_0 = 0, \quad a_0 \oplus b_1 = 0, \\
            a_1 \oplus b_0 = 0, \quad a_1 \oplus b_1 = 1.
        \end{aligned}
    \end{equation}

In the above equation, :math:`a_x` is Alice's answer in the event that she
receives question :math:`x` from the referee for :math:`x \in \Sigma_A`.
Similarly, :math:`b_y` is Bob's answer when Bob receives question :math:`y`
from the referee for :math:`y \in \Sigma_B`. These equations express the
winning conditions that Alice and Bob must satisfy in order to perfectly win
the CHSH game. That is, if it's possible to satisfy all of these equations
simultaneously, it's not possible for them to lose. 

One could perform a brute-force check to see that there is no possible way for
Alice and Bob to simultaneously satisfy all four equations. The best they can
do is satisfy three out of the four equations 

.. math::
    \begin{equation}
        \begin{aligned}
            a_0 \oplus b_0 = 0, \quad a_0 \oplus b_1 = 0, \\
            a_1 \oplus b_0 = 0.
        \end{aligned}
    \end{equation}

They can achieve this if they either have answers :math:`a_0 = b_0 = a_1 = b_1
= 0` or :math:`a_0 = b_0 = a_1 = b_1 = 1`.

Since it is not possible to satisfy all four equations, but it is possible to
satisfy three out of the four equations, the classical value of the CHSH game
is :math:`3/4`, or stated in an equivalent way

.. math::
    \begin{equation}
        \omega(G_{CHSH}) = 3/4 = 0.75.
    \end{equation}

We can verify this by making use of :code:`toqito` to compute the classical
value of the CHSH game.


.. code-block:: python

    >>> from toqito.nonlocal_games.xor_game import XORGame
    >>> chsh = XORGame(prob_mat, pred_mat)
    >>> chsh.classical_value()
    0.75

A quantum strategy for the CHSH game
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

What is very intriguing about the CHSH game is that it is an example of a
nonlocal game where the players can do *strictly better* if they make use of a
quantum strategy instead of a classical one. The quantum strategy that allows
the players to do strictly better is composed of the following shared state and
sets of measurements.

* State: The players prepare and share the state: 

    .. math::
        \begin{equation}
            | \psi \rangle = \frac{1}{\sqrt{2}}
            \left(| 00 \rangle + | 11 \rangle \right).
        \end{equation}

* Measurements: The players measure with respect to the following basis
    
    .. math::
        \begin{equation}
            \begin{aligned}
                | \phi_0 \rangle &= \cos(\theta)|0 \rangle + \sin(\theta)|1 \rangle, \\
                | \phi_1 \rangle &= -\sin(\theta)|0 \rangle + \cos(\theta)|1 \rangle,
            \end{aligned}
        \end{equation}

such that

* If :math:`x = 0` Alice sets :math:`\theta = 0`.
  Otherwise, if :math:`x = 1`, Alice sets :math:`\theta = \pi/4`.

* If :math:`y = 0` Bob sets :math:`\theta = \pi/8`.
  Otherwise, if :math:`y = 1`, Bob sets :math:`\theta = -\pi/8`.

We can now analyze how well this particular quantum strategy performs by
analyzing what occurs in each of the four possible scenarios. For brevity, we
will just analyze the first case, but analyzing the remaining cases follows a
similar analysis.

* Case: :math:`x = 0, y = 0`: 

In this case, Alice and Bob win if :math:`a = b = 0` or if :math:`a = b = 1`.
Alice receives question :math:`x` and selects her measurements constructed from
the basis as specified in the strategy.

.. math::
    \begin{equation}
        A_0^0 = | \phi_0 \rangle \langle \phi_0 |
        \quad \text{and} \quad
        A_1^0 = | \phi_1 \rangle \langle \phi_1 |
    \end{equation}

where 

.. math::
    \begin{equation}
        \begin{aligned}
            | \phi_0 \rangle &= \cos(0)| 0 \rangle + \sin(0)| 1 \rangle, \\
            | \phi_1 \rangle &= -\sin(0)| 0 \rangle + \cos(0)| 1 \rangle.
        \end{aligned}
    \end{equation}

In a similar way, since Bob receives question :math:`y = 0`, he selects his
measurements from the basis

.. math::
    \begin{equation}
        \begin{aligned}
            | \phi_0 \rangle &= \cos(\pi/8)| 0 \rangle + \sin(\pi/8)| 1 \rangle, \\
            | \phi_1 \rangle &= -\sin(\pi/8)| 0 \rangle + \cos(\pi/8)| 1 \rangle.
        \end{aligned}
    \end{equation}

where the measurement operators themselves are defined as

.. math::
    \begin{equation}
        B_0^0 = | \phi_0 \rangle
        \quad \text{and} \quad
        B_1^0 = | \phi_1 \rangle \langle \phi_1 |
    \end{equation}.

Using these measurements, we can calculate the probability that Alice and Bob
win on the inputs :math:`x = 0` and :math:`y = 0` as

.. math::
    \begin{equation}
        p(a, b|0, 0) = \langle \psi | A_0^0 \otimes B_0^0 | \psi \rangle + 
                       \langle \psi | A_1^0 \otimes B_1^0 | \psi \rangle.
    \end{equation}

Calculating the above equation and normalizing by a factor of :math:`1/4`, we
obtain the value of :math:`\cos^2(\pi/8)`. Calculating the remaining three
cases of :math:`(x = 0, y = 1), (x = 1, y = 0)`, and :math:`(x = 1, y = 1)`
follow a similar analysis.

We can see that using this quantum strategy the players win the CHSH game with
a probability of :math:`\cos^2(\pi/8) \approx 0.85355`, which is quite a bit
better than the best classical strategy yielding a probability of :math:`3/4`
to win. As it turns out, the winning probability :math:`\cos^2(\pi/8)` using a
quantum strategy is optimal, which we can represent as
:math:`\omega^*(G_{CHSH}) = \cos^2(\pi/8)`.

We can calculate the quantum value of the CHSH game using :code:`toqito` as
follows:

.. code-block:: python

    >>> chsh.quantum_value()
    0.8535533885683664

For reference, the complete code to calculate both the classical and quantum
values of the CHSH game is provided below.

.. code-block:: python

    >>> import numpy as np
    >>> from toqito.nonlocal_games.xor_game import XORGame
    >>> prob_mat = np.array([[1/4, 1/4],
    >>>                      [1/4, 1/4]])
    >>> pred_mat = np.array([[0, 0],
    >>>                      [0, 1]])
    >>> chsh = XORGame(prob_mat, pred_mat)
    >>> chsh.classical_value()
    0.75
    >>> chsh.quantum_value()
    0.8535533885683664

The odd cycle game
------------------

The *odd cycle game* is another two-player XOR game with the following question and answer sets

.. math::
    \begin{equation}
        \begin{aligned} 
            \Sigma_{A} = \Sigma_B = \mathbb{Z}_n \qquad \text{and} \qquad \Gamma_A = \Gamma_B = \{0, 1\},
        \end{aligned}
    \end{equation}

where :math:`\pi` is the uniform probability distribution over the question set.

As an example, we can specify the odd cycle game for :math:`n=5` and calculate
the classical and quantum values of this game.

.. code-block:: python

    >>> import numpy as np
    >>> from toqito.nonlocal_games.xor_game import XORGame
    >>>
    >>> # Define the probability matrix.
    >>> prob_mat = np.array([
    >>>    [0.1, 0.1, 0, 0, 0],
    >>>    [0, 0.1, 0.1, 0, 0],
    >>>    [0, 0, 0.1, 0.1, 0],
    >>>    [0, 0, 0, 0.1, 0.1],
    >>>    [0.1, 0, 0, 0, 0.1]])
    >>>
    >>> # Define the predicate matrix.
    >>> pred_mat = np.array([
    >>>    [0, 1, 0, 0, 0],
    >>>    [0, 0, 1, 0, 0],
    >>>    [0, 0, 0, 1, 0],
    >>>    [0, 0, 0, 0, 1],
    >>>    [1, 0, 0, 0, 0]])
    >>>
    >>> # Compute the classical and quantum values.
    >>> odd_cycle = XORGame(prob_mat, pred_mat)
    >>> odd_cycle.classical_value()
    0.9
    >>> odd_cycle.quantum_value()
    0.9755282544736033

Note that the odd cycle game is another example of an XOR game where the
players are able to win with a strictly higher probability if they adopt a
quantum strategy. For a general XOR game, Alice and Bob may perform equally
well whether they adopt either a quantum or classical strategy. It holds that
the quantum value for any XOR game is a natural upper bound on the classical
value. That is, for an XOR game, :math:`G`, it holds that

.. math::
    \omega(G) \leq \omega^*(G),

for every XOR game :math:`G`.
    

References
------------------------------

.. [tCSUU08] Cleve, Richard, Slofstra, William, Unger, Falk, and Upadhyay, Sarvagya
    "Perfect parallel repetition theorem for quantum XOR proof systems"
    Computational Complexity 17.2 (2008): 282-299.
    https://arxiv.org/abs/quant-ph/0608146

﻿toqito.channels.partial\_transpose
==================================

.. currentmodule:: toqito.channels

.. autofunction:: partial_transpose﻿toqito.perms.permutation\_operator
==================================

.. currentmodule:: toqito.perms

.. autofunction:: permutation_operator﻿toqito.states.basis
===================

.. currentmodule:: toqito.states

.. autofunction:: basis﻿toqito.state\_metrics.sub\_fidelity
===================================

.. currentmodule:: toqito.state_metrics

.. autofunction:: sub_fidelity﻿toqito.state\_props.concurrence
===============================

.. currentmodule:: toqito.state_props

.. autofunction:: concurrence﻿toqito.matrix\_props.is\_identity
=================================

.. currentmodule:: toqito.matrix_props

.. autofunction:: is_identity﻿toqito.matrix\_props.is\_square
===============================

.. currentmodule:: toqito.matrix_props

.. autofunction:: is_square﻿toqito.states.w\_state
======================

.. currentmodule:: toqito.states

.. autofunction:: w_state﻿toqito.state\_opt.state\_helper
===============================

.. automodule:: toqito.state_opt.state_helper

   
   
   

   
   
   

   
   
   

   
   
   



﻿toqito.measurement\_props.is\_povm
==================================

.. currentmodule:: toqito.measurement_props

.. autofunction:: is_povm﻿toqito.state\_opt.ppt\_distinguishability
=========================================

.. currentmodule:: toqito.state_opt

.. autofunction:: ppt_distinguishability﻿toqito.perms.antisymmetric\_projection
======================================

.. currentmodule:: toqito.perms

.. autofunction:: antisymmetric_projection﻿toqito.states.singlet
=====================

.. currentmodule:: toqito.states

.. autofunction:: singlet﻿toqito.matrices.pauli
=====================

.. currentmodule:: toqito.matrices

.. autofunction:: pauli﻿toqito.channel\_ops.kraus\_to\_choi
===================================

.. currentmodule:: toqito.channel_ops

.. autofunction:: kraus_to_choi﻿toqito.state\_props.log\_negativity
===================================

.. currentmodule:: toqito.state_props

.. autofunction:: log_negativity﻿toqito.measurement\_ops.measure
===============================

.. currentmodule:: toqito.measurement_ops

.. autofunction:: measure﻿toqito.state\_opt.optimal\_clone
================================

.. currentmodule:: toqito.state_opt

.. autofunction:: optimal_clone﻿toqito.random.random\_density\_matrix
=====================================

.. currentmodule:: toqito.random

.. autofunction:: random_density_matrix﻿toqito.channels.reduction
=========================

.. currentmodule:: toqito.channels

.. autofunction:: reduction﻿toqito.state\_opt.state\_distinguishability
===========================================

.. currentmodule:: toqito.state_opt

.. autofunction:: state_distinguishability﻿toqito.matrix\_props.is\_density
================================

.. currentmodule:: toqito.matrix_props

.. autofunction:: is_density﻿toqito.matrix\_props.majorizes
==============================

.. currentmodule:: toqito.matrix_props

.. autofunction:: majorizes﻿toqito.matrix\_ops.unvec
========================

.. currentmodule:: toqito.matrix_ops

.. autofunction:: unvec﻿toqito.random.random\_povm
==========================

.. currentmodule:: toqito.random

.. autofunction:: random_povm﻿toqito.state\_props.entanglement\_of\_formation
===============================================

.. currentmodule:: toqito.state_props

.. autofunction:: entanglement_of_formation﻿toqito.states.brauer
====================

.. currentmodule:: toqito.states

.. autofunction:: brauer﻿toqito.channel\_props.choi\_rank
================================

.. currentmodule:: toqito.channel_props

.. autofunction:: choi_rank﻿toqito.channel\_props.is\_trace\_preserving
===========================================

.. currentmodule:: toqito.channel_props

.. autofunction:: is_trace_preserving﻿toqito.matrix\_props.is\_commuting
==================================

.. currentmodule:: toqito.matrix_props

.. autofunction:: is_commuting﻿toqito.random.random\_ginibre
=============================

.. currentmodule:: toqito.random

.. autofunction:: random_ginibre﻿toqito.matrix\_props.is\_diagonal
=================================

.. currentmodule:: toqito.matrix_props

.. autofunction:: is_diagonal﻿toqito.perms.perfect\_matchings
===============================

.. currentmodule:: toqito.perms

.. autofunction:: perfect_matchings﻿toqito.states.gisin
===================

.. currentmodule:: toqito.states

.. autofunction:: gisin﻿toqito.channels.choi
====================

.. currentmodule:: toqito.channels

.. autofunction:: choi﻿toqito.state\_props.is\_ppt
===========================

.. currentmodule:: toqito.state_props

.. autofunction:: is_ppt﻿toqito.matrix\_props.is\_projection
===================================

.. currentmodule:: toqito.matrix_props

.. autofunction:: is_projection﻿toqito.matrices.shift
=====================

.. currentmodule:: toqito.matrices

.. autofunction:: shift﻿toqito.states.gen\_bell
=======================

.. currentmodule:: toqito.states

.. autofunction:: gen_bell﻿toqito.channels.realignment
===========================

.. currentmodule:: toqito.channels

.. autofunction:: realignment﻿toqito.state\_props.is\_ensemble
================================

.. currentmodule:: toqito.state_props

.. autofunction:: is_ensemble﻿toqito.states.bell
==================

.. currentmodule:: toqito.states

.. autofunction:: bell﻿toqito.perms.perm\_sign
=======================

.. currentmodule:: toqito.perms

.. autofunction:: perm_sign﻿toqito.channel\_ops.choi\_to\_kraus
===================================

.. currentmodule:: toqito.channel_ops

.. autofunction:: choi_to_kraus﻿toqito.matrix\_props.is\_symmetric
==================================

.. currentmodule:: toqito.matrix_props

.. autofunction:: is_symmetric﻿toqito.channel\_props.is\_completely\_positive
==============================================

.. currentmodule:: toqito.channel_props

.. autofunction:: is_completely_positive﻿toqito.matrix\_props.is\_positive\_definite
===========================================

.. currentmodule:: toqito.matrix_props

.. autofunction:: is_positive_definite﻿toqito.states.max\_mixed
========================

.. currentmodule:: toqito.states

.. autofunction:: max_mixed﻿toqito.state\_props.is\_mixed
=============================

.. currentmodule:: toqito.state_props

.. autofunction:: is_mixed﻿toqito.state\_props.negativity
==============================

.. currentmodule:: toqito.state_props

.. autofunction:: negativity﻿toqito.state\_metrics.hilbert\_schmidt
======================================

.. currentmodule:: toqito.state_metrics

.. autofunction:: hilbert_schmidt﻿toqito.states.ghz
=================

.. currentmodule:: toqito.states

.. autofunction:: ghz﻿toqito.states.domino
====================

.. currentmodule:: toqito.states

.. autofunction:: domino﻿toqito.states.isotropic
=======================

.. currentmodule:: toqito.states

.. autofunction:: isotropic﻿toqito.matrices.iden
====================

.. currentmodule:: toqito.matrices

.. autofunction:: iden﻿toqito.matrix\_props.is\_unitary
================================

.. currentmodule:: toqito.matrix_props

.. autofunction:: is_unitary﻿toqito.perms.swap\_operator
===========================

.. currentmodule:: toqito.perms

.. autofunction:: swap_operator﻿toqito.matrix\_ops.tensor
=========================

.. currentmodule:: toqito.matrix_ops

.. autofunction:: tensor﻿toqito.matrices.cnot
====================

.. currentmodule:: toqito.matrices

.. autofunction:: cnot﻿toqito.states.werner
====================

.. currentmodule:: toqito.states

.. autofunction:: werner﻿toqito.state\_props.l1\_norm\_coherence
=======================================

.. currentmodule:: toqito.state_props

.. autofunction:: l1_norm_coherence﻿toqito.states.breuer
====================

.. currentmodule:: toqito.states

.. autofunction:: breuer﻿toqito.matrix\_props.is\_idempotent
===================================

.. currentmodule:: toqito.matrix_props

.. autofunction:: is_idempotent﻿toqito.channel\_props.is\_positive
==================================

.. currentmodule:: toqito.channel_props

.. autofunction:: is_positive﻿toqito.state\_metrics.trace\_distance
=====================================

.. currentmodule:: toqito.state_metrics

.. autofunction:: trace_distance﻿toqito.state\_metrics.matsumoto\_fidelity
=========================================

.. currentmodule:: toqito.state_metrics

.. autofunction:: matsumoto_fidelity﻿toqito.matrix\_ops.vec
======================

.. currentmodule:: toqito.matrix_ops

.. autofunction:: vec﻿toqito.state\_ops.schmidt\_decomposition
========================================

.. currentmodule:: toqito.state_ops

.. autofunction:: schmidt_decomposition﻿toqito.channel\_props.is\_unitary
=================================

.. currentmodule:: toqito.channel_props

.. autofunction:: is_unitary﻿toqito.channels.dephasing
=========================

.. currentmodule:: toqito.channels

.. autofunction:: dephasing﻿toqito.state\_props.is\_mutually\_unbiased\_basis
=================================================

.. currentmodule:: toqito.state_props

.. autofunction:: is_mutually_unbiased_basis﻿toqito.channel\_ops.apply\_channel
==================================

.. currentmodule:: toqito.channel_ops

.. autofunction:: apply_channel﻿toqito.channel\_props.is\_unital
================================

.. currentmodule:: toqito.channel_props

.. autofunction:: is_unital﻿toqito.matrices.gen\_pauli
==========================

.. currentmodule:: toqito.matrices

.. autofunction:: gen_pauli﻿toqito.matrices.gen\_gell\_mann
===============================

.. currentmodule:: toqito.matrices

.. autofunction:: gen_gell_mann﻿toqito.matrices.gell\_mann
==========================

.. currentmodule:: toqito.matrices

.. autofunction:: gell_mann﻿toqito.state\_props.is\_pure
============================

.. currentmodule:: toqito.state_props

.. autofunction:: is_pure﻿toqito.state\_opt.state\_exclusion
==================================

.. currentmodule:: toqito.state_opt

.. autofunction:: state_exclusion﻿toqito.perms.swap
=================

.. currentmodule:: toqito.perms

.. autofunction:: swap﻿toqito.states.chessboard
========================

.. currentmodule:: toqito.states

.. autofunction:: chessboard﻿toqito.channel\_ops.partial\_channel
====================================

.. currentmodule:: toqito.channel_ops

.. autofunction:: partial_channel﻿toqito.perms.unique\_perms
==========================

.. currentmodule:: toqito.perms

.. autofunction:: unique_perms﻿toqito.state\_metrics.fidelity
==============================

.. currentmodule:: toqito.state_metrics

.. autofunction:: fidelity﻿toqito.state\_ops.pure\_to\_mixed
=================================

.. currentmodule:: toqito.state_ops

.. autofunction:: pure_to_mixed﻿toqito.random.random\_state\_vector
===================================

.. currentmodule:: toqito.random

.. autofunction:: random_state_vector﻿toqito.state\_props.in\_separable\_ball
=======================================

.. currentmodule:: toqito.state_props

.. autofunction:: in_separable_ball﻿toqito.matrix\_props.is\_normal
===============================

.. currentmodule:: toqito.matrix_props

.. autofunction:: is_normal﻿toqito.state\_props.has\_symmetric\_extension
=============================================

.. currentmodule:: toqito.state_props

.. autofunction:: has_symmetric_extension﻿toqito.channels.partial\_trace
==============================

.. currentmodule:: toqito.channels

.. autofunction:: partial_trace﻿toqito.states.tile
==================

.. currentmodule:: toqito.states

.. autofunction:: tile﻿toqito.channel\_ops.dual\_channel
=================================

.. currentmodule:: toqito.channel_ops

.. autofunction:: dual_channel﻿toqito.state\_props.purity
==========================

.. currentmodule:: toqito.state_props

.. autofunction:: purity﻿toqito.perms.permute\_systems
=============================

.. currentmodule:: toqito.perms

.. autofunction:: permute_systems﻿toqito.channel\_props.is\_quantum\_channel
==========================================

.. currentmodule:: toqito.channel_props

.. autofunction:: is_quantum_channel﻿toqito.channels.depolarizing
============================

.. currentmodule:: toqito.channels

.. autofunction:: depolarizing﻿toqito.state\_metrics.helstrom\_holevo
======================================

.. currentmodule:: toqito.state_metrics

.. autofunction:: helstrom_holevo﻿toqito.state\_props.is\_product
===============================

.. currentmodule:: toqito.state_props

.. autofunction:: is_product﻿toqito.state\_metrics.bures\_distance
=====================================

.. currentmodule:: toqito.state_metrics

.. autofunction:: bures_distance﻿toqito.state\_opt.symmetric\_extension\_hierarchy
=================================================

.. currentmodule:: toqito.state_opt

.. autofunction:: symmetric_extension_hierarchy﻿toqito.perms.symmetric\_projection
==================================

.. currentmodule:: toqito.perms

.. autofunction:: symmetric_projection﻿toqito.states.horodecki
=======================

.. currentmodule:: toqito.states

.. autofunction:: horodecki﻿toqito.random.random\_unitary
=============================

.. currentmodule:: toqito.random

.. autofunction:: random_unitary﻿toqito.state\_props.is\_mutually\_orthogonal
============================================

.. currentmodule:: toqito.state_props

.. autofunction:: is_mutually_orthogonal﻿toqito.states.max\_entangled
============================

.. currentmodule:: toqito.states

.. autofunction:: max_entangled﻿toqito.state\_props.schmidt\_rank
=================================

.. currentmodule:: toqito.state_props

.. autofunction:: schmidt_rank﻿toqito.channel\_props.is\_herm\_preserving
==========================================

.. currentmodule:: toqito.channel_props

.. autofunction:: is_herm_preserving﻿toqito.state\_props.von\_neumann\_entropy
=========================================

.. currentmodule:: toqito.state_props

.. autofunction:: von_neumann_entropy﻿toqito.channel\_metrics.channel\_fidelity
=========================================

.. currentmodule:: toqito.channel_metrics

.. autofunction:: channel_fidelity﻿toqito.matrices.fourier
=======================

.. currentmodule:: toqito.matrices

.. autofunction:: fourier﻿toqito.matrix\_props.is\_hermitian
==================================

.. currentmodule:: toqito.matrix_props

.. autofunction:: is_hermitian﻿toqito.matrix\_props.is\_permutation
====================================

.. currentmodule:: toqito.matrix_props

.. autofunction:: is_permutation﻿toqito.state\_metrics.trace\_norm
=================================

.. currentmodule:: toqito.state_metrics

.. autofunction:: trace_norm﻿toqito.matrix\_props.is\_positive\_semidefinite
===============================================

.. currentmodule:: toqito.matrix_props

.. autofunction:: is_positive_semidefinite﻿toqito.matrices.clock
=====================

.. currentmodule:: toqito.matrices

.. autofunction:: clock﻿toqito.matrices.hadamard
========================

.. currentmodule:: toqito.matrices

.. autofunction:: hadamard