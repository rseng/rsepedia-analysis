<a href="http://mc-stan.org">
<img src="https://raw.githubusercontent.com/stan-dev/logos/master/logo.png" width=200 alt="Stan Logo"/>
</a>

<b>Stan</b> is a C++ package providing

* full Bayesian inference using the No-U-Turn sampler (NUTS), a variant of Hamiltonian Monte Carlo (HMC),
* approximate Bayesian inference using automatic differentiation variational inference (ADVI), and
* penalized maximum likelihood estimation (MLE) using L-BFGS optimization.

It is built on top of the [Stan Math library](https://github.com/stan-dev/math), which provides

* a full first- and higher-order automatic differentiation library based on C++ template overloads, and
* a supporting fully-templated matrix, linear algebra, and probability special function library.

There are interfaces available in R, Python, MATLAB, Julia, Stata, Mathematica, and for the command line.

[![DOI](https://zenodo.org/badge/19868/stan-dev/stan.svg)](https://zenodo.org/badge/latestdoi/19868/stan-dev/stan)

Home Page
---------
Stan's home page, with links to everything you'll need to use Stan is:

[http://mc-stan.org/](http://mc-stan.org/)

Interfaces
----------
There are separate repositories in the stan-dev GitHub organization for the interfaces, higher-level libraries and lower-level libraries.  

Source Repository
-----------------
Stan's source-code repository is hosted here on GitHub.

Licensing
---------
The Stan math library, core Stan code, and CmdStan are licensed under new BSD. RStan and PyStan are licensed under GPLv3, with other interfaces having other open-source licenses.

Note that the Stan math library depends on the Intel TBB library which is licensed under the Apache 2.0 license. This dependency implies an additional restriction as compared to the new BSD lincense alone. The Apache 2.0 license is incompatible with GPL-2 licensed code if distributed as a unitary binary. You may refer to the Apache 2.0 evaluation page on the [Stan Math wiki](https://github.com/stan-dev/math/wiki/Apache-2.0-License-Evaluation).
#### Summary:
Please provide a short couple sentence summary.


#### Description:
Describe the issue as clearly as possible.


#### Reproducible Steps:
Please report steps to reproduce the issue. If it's not possible to reproduce, please include a description of how you discovered the issue.

If you have a reproducible example, please include it.


#### Current Output:
The current output. Knowing what is the current behavior is useful.


#### Expected Output:
Describe what you expect the output to be. Knowing the correct behavior is also very useful.


#### Additional Information:
Provide any additional information here.

#### Current Version:
v2.28.2
## Contributing Code

Start here: [Introduction to Stan for New Developers](https://github.com/stan-dev/stan/wiki/Introduction-to-Stan-for-New-Developers)

The instructions for contributors including the Git process, the testing process, and the Stan coding standards, see the top-level Wiki:

* [Stan Developer Process](https://github.com/stan-dev/stan/wiki)

## Stan Manual Issues

If you would like to make a suggestion for fixing a typo or unclear part of the manual, please make a comment on the currently open next manual issue (it will be numbered by current release followed by `++`) and we will fix it in the next release after the current release:

* [Next Manual Issue](https://github.com/stan-dev/stan/issues?utf8=✓&q=is%3Aopen%20label%3ABug%20%22next%20manual%22%20label%3ADocumentation)

This is easier for us than having to deal with a typo pull request, though we will gladly take those, too.

## Submitting an Issue

If you want to report a bug or request a feature, the right place to do this is on each individual project's issue tracker.  Please first search to see if the issue has already been reported, and if you need to report a new issue, click on the `New` button in the upper right corner.

If you suspect a mathematical function is returning the wrong results or want to request a new function be added to Stan:

* [Stan Math Library](https://github.com/stan-dev/math/issues)

If you suspect a problem with the Stan language, the fitting algorithms, or with I/O, or if you want to see a new feature in all of the interfaces:

* [Stan Language, Algorithms, and Services](https://github.com/stan-dev/stan/issues)

If you suspect a problem or would like to see a new feature in one of the interfaces, choose the interface:

* [CmdStan Interface](https://github.com/stan-dev/cmdstan/issues)
* [RStan](https://github.com/stan-dev/rstan/issues)
* [PyStan](https://github.com/stan-dev/pystan/issues)
* [Stan.jl](https://github.com/goedman/Stan.jl/issues)
* [StataStan](https://github.com/stan-dev/statastan/issues)
* [MathematicaStan](https://github.com/stan-dev/MathematicaStan/issues)

If you suspect a problem or would like to request a feature in one of our higher-level packages:

* [RStanArm](https://github.com/stan-dev/rstanarm/issues)
* [ShinyStan](https://github.com/stan-dev/shinystan/issues)
* [BayesPlot](https://github.com/stan-dev/bayesplot/issues)
* [Loo](https://github.com/stan-dev/loo/issues)

If you have a problem with Stan's emacs mode or would like to request a feature:

* [Stan Emacs Mode](https://github.com/stan-dev/stan-mode/issues)

If you suspect a problem or would like to see something added to the web site (http://mc-stan.org)[http://mc-stan.org]:

* [Stan Web Pages](https://github.com/stan-dev/stan-dev.github.io/issues)

#### Submission Checklist

- [ ] Run unit tests: `./runTests.py src/test/unit`
- [ ] Run cpplint: `make cpplint`
- [ ] Declare copyright holder and open-source license: see below

#### Summary

#### Intended Effect

#### How to Verify

#### Side Effects

#### Documentation

#### Copyright and Licensing

Please list the copyright holder for the work you are submitting (this will be you or your assignee, such as a university or company):



By submitting this pull request, the copyright holder is agreeing to license the submitted work under the following licenses:
- Code: BSD 3-clause (https://opensource.org/licenses/BSD-3-Clause)
- Documentation: CC-BY 4.0 (https://creativecommons.org/licenses/by/4.0/)
