![JuMP logo](https://jump.dev/JuMP.jl/dev/assets/logo-with-text-background.svg "JuMP logo")
---

[![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](https://numfocus.org)

JuMP is a domain-specific modeling language for [mathematical optimization](https://en.wikipedia.org/wiki/Mathematical_optimization)
embedded in [Julia](https://julialang.org/). You can find out more about us by
visiting [jump.dev](https://jump.dev).


**Latest Release**: [![version](https://juliahub.com/docs/JuMP/DmXqY/0.22.2/version.svg)](https://juliahub.com/ui/Packages/JuMP/DmXqY/0.22.2) (`release-0.22` branch):
  * Installation via the Julia package manager:
    * `import Pkg; Pkg.add("JuMP")`
  * Get help:
    * Read the [Documentation](https://jump.dev/JuMP.jl/stable/)
    * Ask a question on the [Community forum]
  * Testing status:
    * Github Actions: [![Build Status](https://github.com/jump-dev/JuMP.jl/workflows/CI/badge.svg?branch=release-0.22)](https://github.com/jump-dev/JuMP.jl/actions?query=workflow%3ACI)
  * [![deps](https://juliahub.com/docs/JuMP/deps.svg)](https://juliahub.com/ui/Packages/JuMP/DmXqY?t=2)

**Development version** (`master` branch):
  * Installation via the Julia package manager:
    * `import Pkg; Pkg.add(Pkg.PackageSpec(name="JuMP", rev="master"))`
  * Get help:
    * Read the [Documentation](https://jump.dev/JuMP.jl/dev/)
    * Join the [Developer chatroom](https://gitter.im/JuliaOpt/JuMP-dev)
    * Read the [NEWS](https://github.com/jump-dev/JuMP.jl/tree/master/NEWS.md)
  * Testing status:
    * Github Actions: [![Build Status](https://github.com/jump-dev/JuMP.jl/workflows/CI/badge.svg?branch=master)](https://github.com/jump-dev/JuMP.jl/actions?query=workflow%3ACI)
    * Test coverage: [![codecov](https://codecov.io/gh/jump-dev/JuMP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/JuMP.jl)

## Need help?

Use the [Community forum] to search for answers to previously asked questions,
or ask a new question.

The post [Please read: make it easier to help you](https://discourse.julialang.org/t/please-read-make-it-easier-to-help-you/14757),
describes the best practices for asking a question.

[Community forum]: https://discourse.julialang.org/c/domain/opt

## Bug reports

Please report any issues via the Github [issue tracker]. All types of issues
are welcome and encouraged; this includes bug reports, documentation typos,
feature requests, etc.

[issue tracker]: https://github.com/jump-dev/JuMP.jl/issues

## Citing JuMP

If you find JuMP useful in your work, we kindly request that you cite the
following paper ([pdf](https://mlubin.github.io/pdf/jump-sirev.pdf)):

```bibtex
@article{DunningHuchetteLubin2017,
    author = {Iain Dunning and Joey Huchette and Miles Lubin},
    title = {JuMP: A Modeling Language for Mathematical Optimization},
    journal = {SIAM Review},
    volume = {59},
    number = {2},
    pages = {295-320},
    year = {2017},
    doi = {10.1137/15M1020575},
}
```

For an earlier work where we presented a prototype implementation of JuMP, see
[here](https://dx.doi.org/10.1287/ijoc.2014.0623):

```bibtex
@article{LubinDunningIJOC,
    author = {Miles Lubin and Iain Dunning},
    title = {Computing in Operations Research Using Julia},
    journal = {INFORMS Journal on Computing},
    volume = {27},
    number = {2},
    pages = {238-248},
    year = {2015},
    doi = {10.1287/ijoc.2014.0623},
}
```

A preprint of this paper is [freely available](https://arxiv.org/abs/1312.1431).

---

![NumFOCUS logo](https://jump.dev/JuMP.jl/dev/assets/numfocus-logo.png)

JuMP is a Sponsored Project of NumFOCUS, a 501(c)(3) nonprofit charity in the
United States. NumFOCUS provides JuMP with fiscal, legal, and administrative
support to help ensure the health and sustainability of the project. Visit
[numfocus.org](https://numfocus.org) for more information.

You can support JuMP by [donating](https://numfocus.salsalabs.org/donate-to-jump/index.html).

Donations to JuMP are managed by NumFOCUS. For donors in the United States,
your gift is tax-deductible to the extent provided by law. As with any donation,
you should consult with your tax adviser about your particular tax situation.

JuMP's largest expense is the annual JuMP-dev workshop. Donations will help us
provide travel support for JuMP-dev attendees and take advantage of other
opportunities that arise to support JuMP development.
# JuMP release notes

Release notes for JuMP have been moved to the documentation.

See https://jump.dev/JuMP.jl/latest/release_notes/
# JuMP community standards

Welcome! JuMP is part of the Julia community, and we uphold its [community
standards](https://julialang.org/community/standards/). We repeat the key points
below with some minor specializations for JuMP.

The JuMP community is committed to maintaining a welcoming, civil, and
constructive environment. We expect the following standards to be observed and
upheld by all participants in any community forum (Discourse, GitHub, gitter,
etc.).

## Be respectful and inclusive.

Please do not use overtly sexual language or imagery, and do not attack anyone
based on any aspect of personal identity, including gender, sexuality, religion,
ethnicity, race, age or ability. Keep in mind that what you write in public
forums is read by many people who don't know you personally, so please refrain
from making prejudiced or sexual jokes and comments – even ones that you might
consider acceptable in private. Ask yourself if a comment or statement might
make someone feel unwelcomed or like an outsider.

## Give credit.

All participants in the JuMP community are expected to respect copyright laws
and ethical attribution standards. This applies to both code and written
materials, such as documentation or blog posts. Materials that violate the law,
are plagiaristic, or ethically dubious in some way will be removed from
officially-maintained lists of resources.

## Get involved.

The JuMP community is built on a foundation of reciprocity and collaboration. Be
aware that most community members contribute on a voluntary basis, so ideas and
bug reports are ok, but demands are not. Along these lines, it's important to
keep in mind that, unless stated otherwise, JuMP's interfaces to commercial
solvers are maintained by volunteers who find them useful. Pull requests are
always welcomed. Please check out our [contributing guide](https://github.com/jump-dev/JuMP.jl/blob/master/CONTRIBUTING.md).
You can also reach out in the [developer chatroom](https://gitter.im/JuliaOpt/JuMP-dev)
with questions about how to get started.

## Any concerns?

If you have a conflict or concern that requires resolution, please contact the
[Julia Community Stewards](https://julialang.org/community/stewards/).
# How to contribute to JuMP

The contribution guide has been moved to the documentation.

See https://jump.dev/JuMP.jl/latest/developers/contributing
---
name: Bug report
about: Help us track down bugs in JuMP

---

Welcome to JuMP!

Please read the following before posting a new bug report:

- If you have a question or are unsure if the behavior you're experiencing is a bug, please search or post to our Discourse site: https://discourse.julialang.org/c/domain/opt. Questions posted to Discourse have broader visibility and are likely to be answered more quickly than issues filed here.

- If you're experiencing a bug that is solver-specific (e.g., it only happens when you use Gurobi), the issue is best raised at the relevant solver-wrapper repository (e.g., Gurobi.jl). If you are unsure, you should raise the problem on Discourse first.

- If you are reasonably confident your issue is a bug in JuMP, this is the right place. Be sure to include as much relevant information as possible, including a minimal reproducible example. See https://help.github.com/articles/basic-writing-and-formatting-syntax/ for background on how to format text and code on GitHub issues.

Thanks for contributing to JuMP!
---
name: Feature request
about: Suggest an idea for JuMP

---

Welcome to JuMP! You can use this Github issue to suggest a new feature for JuMP.

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.

Thanks for contributing to JuMP!
# JuMP Documentation README

**The documentation currently requires Julia 1.6 to build.**

JuMP's documentation is written with [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl).

## Initial setup

To build the documentation, you need to do a series of initialization steps.
However, you only need to do this once!

First, you will need a local copy of JuMP. If you don't have one already, run:
```
$ julia -e 'import Pkg; Pkg.develop("JuMP")'
```

This will create a copy of JuMP at `~/.julia/dev/JuMP`. (On Windows, this will
be located at `C:\\Users\\<your_user_name>\\.julia\\dev\\JuMP`.) Open a terminal,
and `cd` to that directory:
```
$ cd ~/.julia/dev/JuMP
```

The next step is to setup the `docs` environment.
```
$ julia --project=docs -e 'import Pkg; Pkg.instantiate(); Pkg.develop(Pkg.PackageSpec(path="."))'
```

Now you're ready to build the documentation.

## Building the docs

Build the docs as follows:
```
$ cd ~/.julia/dev/JuMP
$ julia --project=docs docs/make.jl
```

The compiled documents can be viewed at `~/.julia/dev/JuMP/docs/build/index.html`.
```@raw html
<img class="display-light-only" src="assets/logo-with-text.svg" alt="JuMP logo"/>
<img class="display-dark-only" src="assets/logo-dark-with-text.svg" alt="JuMP logo"/>
```

```@meta
# These comments do not display in the HTML output.
# See https://github.com/JuliaDocs/Documenter.jl/issues/674.
```

# Introduction

Welcome to the documentation for JuMP!

## What is JuMP?

[JuMP](https://github.com/jump-dev/JuMP.jl) is a domain-specific modeling
language for [mathematical optimization](https://en.wikipedia.org/wiki/Mathematical_optimization)
embedded in [Julia](https://julialang.org/). It currently supports a number of
open-source and commercial solvers for a variety of problem classes, including
linear, mixed-integer, second-order conic, semidefinite, and nonlinear
programming.

!!! tip
    If you aren't sure if you should use JuMP, read [Should I use JuMP?](@ref).

## Resources for getting started

There are few ways to get started with JuMP:

* Read the [Installation Guide](@ref).
* Read the introductory tutorials [Getting started with Julia](@ref) and
  [Getting started with JuMP](@ref).
* Browse some of our modeling tutorials, including classics such as
  [The diet problem](@ref), or the [Maximum likelihood estimation](@ref) problem
  using nonlinear programming.

!!! tip
    Need help? Join the [community forum](https://discourse.julialang.org/c/domain/opt/13)
    to search for answers to commonly asked questions.

    Before asking a question, make sure to read the post [make it easier to help you](https://discourse.julialang.org/t/psa-make-it-easier-to-help-you/14757),
    which contains a number of tips on how to ask a good question.

## How the documentation is structured

Having a high-level overview of how this documentation is structured will help
you know where to look for certain things.

* **Tutorials** contain worked examples of solving problems with JuMP. Start
  here if you are new to JuMP, or you have a particular problem class you want
  to model.

* The **Manual** contains short code-snippets that explain how to achieve
  specific tasks in JuMP. Look here if you want to know how to achieve a
  particular task, such as how to [Delete a variable](@ref delete_a_variable) or
  how to [Modify an objective coefficient](@ref).

* The **API Reference** contains a complete list of the functions you can use in
  JuMP. Look here if you want to know how to use a particular function.

* The **Background information** section contains background reading material to
  provide context to JuMP. Look here if you want an understanding of what JuMP
  is and why we created it, rather than how to use it.

* The **Developer docs** section contains information for people contributing to
  JuMP development or writing JuMP extensions. Don't worry about this section if
  you are using JuMP to formulate and solve problems as a user.

* The **MathOptInterface** section is a self-contained copy of the documentation
  for MathOptInterface. Look here for functions and constants beginning with
  `MOI.`, as well as for general information on how MathOptInterface works.

## Citing JuMP

If you find JuMP useful in your work, we kindly request that you cite the
following paper ([pdf](https://mlubin.github.io/pdf/jump-sirev.pdf)):

``` sourceCode
@article{DunningHuchetteLubin2017,
author = {Iain Dunning and Joey Huchette and Miles Lubin},
title = {JuMP: A Modeling Language for Mathematical Optimization},
journal = {SIAM Review},
volume = {59},
number = {2},
pages = {295-320},
year = {2017},
doi = {10.1137/15M1020575},
}
```

For an earlier work where we presented a prototype implementation of JuMP, see
[here](https://dx.doi.org/10.1287/ijoc.2014.0623):

``` sourceCode
@article{LubinDunningIJOC,
author = {Miles Lubin and Iain Dunning},
title = {Computing in Operations Research Using Julia},
journal = {INFORMS Journal on Computing},
volume = {27},
number = {2},
pages = {238-248},
year = {2015},
doi = {10.1287/ijoc.2014.0623},
}
```

A preprint of this paper is [freely available](https://arxiv.org/abs/1312.1431).

## NumFOCUS

![NumFOCUS logo](assets/numfocus-logo.png)

JuMP is a Sponsored Project of NumFOCUS, a 501(c)(3) nonprofit charity in the
United States. NumFOCUS provides JuMP with fiscal, legal, and administrative
support to help ensure the health and sustainability of the project. Visit
[numfocus.org](https://numfocus.org) for more information.

You can support JuMP by [donating](https://numfocus.salsalabs.org/donate-to-jump/index.html).

Donations to JuMP are managed by NumFOCUS. For donors in the United States,
your gift is tax-deductible to the extent provided by law. As with any donation,
you should consult with your tax adviser about your particular tax situation.

JuMP's largest expense is the annual JuMP-dev workshop. Donations will help us
provide travel support for JuMP-dev attendees and take advantage of other
opportunities that arise to support JuMP development.
# Installation Guide

This guide explains how to install Julia and  JuMP. If you have installation
troubles, read the [Common installation issues](@ref) section below.

## Install Julia

JuMP is a package for [Julia](https://julialang.org). To use JuMP, first
[download and install](https://julialang.org/downloads/) Julia.

!!! tip
    If you  are new  to Julia, read our [Getting started with Julia](@ref)
    tutorial.

### Which version should I pick?

You can install the "Current stable release" or the "Long-term support (LTS)
release".

 * The "Current stable release" is the latest release of Julia. It has access to
   newer features, and is likely faster.
 * The "Long-term support release" is an older version of Julia that has
   continued to receive bug and security fixes. However, it may not have the
   latest features or performance improvements.

For most users, you should install the "Current stable release", and whenever
Julia releases a new version of the current stable release, you should update
your version of Julia. Note that any code you write on one version of the
current stable release will continue to work on all subsequent releases.

For users in restricted software environments (e.g., your enterprise IT controls
what software you can install), you may be better off installing the long-term
support release because you will not have to update Julia as frequently.

## Install JuMP

From Julia, JuMP is installed using the built-in package manager:
```julia
import Pkg
Pkg.add("JuMP")
```

!!! tip
    We recommend you create a Pkg _environment_ for each project you use JuMP
    for, instead of adding lots of packages to the global environment. The
    [Pkg manager documentation](https://julialang.github.io/Pkg.jl/v1/environments/)
    has more information on this topic.

When we release a new version of JuMP, you can update with:
```julia
import Pkg
Pkg.update("JuMP")
```

## Install a solver

JuMP depends on solvers to solve optimization problems. Therefore, you will need
to install one before you can solve problems with JuMP.

Install a solver using the Julia package manager, replacing `"Clp"` by the
Julia package name as appropriate.
```julia
import Pkg
Pkg.add("Clp")
```

Once installed, you can use Clp as a solver with JuMP as follows, using
[`set_optimizer_attributes`](@ref) to set solver-specific options:
```julia
using JuMP
using Clp
model = Model(Clp.Optimizer)
set_optimizer_attributes(model, "LogLevel" => 1, "PrimalTolerance" => 1e-7)
```

!!! note
    Most packages follow the `ModuleName.Optimizer` naming convention, but
    exceptions may exist. See the README of the Julia package's GitHub
    repository for more details on how to use a particular solver, including any
    solver-specific options.

## Supported solvers

Most solvers are not written in Julia, and some require commercial licenses to
use, so installation is often more complex.
  * If a solver has `Manual` in the `Installation` column, the solver requires a
    manual installation step, such as downloading and installing a binary, or
    obtaining a commercial license. Consult the README of the relevant Julia
    package for more information.
  * If the solver has `Manualᴹ` in the `Installation` column, the solver
    requires an installation of [MATLAB](https://www.mathworks.com/products/matlab.html).
  * If the `Installation` column is missing an entry, installing the Julia
    package will download and install any relevant solver binaries
    automatically, and you shouldn't need to do anything other than `Pkg.add`.

Solvers with a missing entry in the `Julia Package` column are written in Julia.
The link in the `Solver` column is the corresponding Julia package.

| Solver                                                                         | Julia Package                                                                    | Installation | License | Supports             |
| ------------------------------------------------------------------------------ | -------------------------------------------------------------------------------- | ------------ | ------- | ---------------------|
| [Alpine.jl](https://github.com/lanl-ansi/Alpine.jl)                            |                                                                                  |        | Triad NS | (MI)NLP                   |
| [Artelys Knitro](https://www.artelys.com/knitro)                               | [KNITRO.jl](https://github.com/jump-dev/KNITRO.jl)                               | Manual | Comm.    | (MI)LP, (MI)SOCP, (MI)NLP |
| [BARON](http://minlp.com/baron)                                                | [BARON.jl](https://github.com/joehuchette/BARON.jl)                              | Manual | Comm.    | (MI)NLP                   |
| [Bonmin](http://github.com/coin-or/Bonmin)                                     | [AmplNLWriter.jl](https://github.com/jump-dev/AmplNLWriter.jl)                   |        | EPL      | (MI)NLP                   |
| [Cbc](https://github.com/coin-or/Cbc)                                          | [Cbc.jl](https://github.com/jump-dev/Cbc.jl)                                     |        | EPL      | (MI)LP                    |
| [CDCS](https://github.com/oxfordcontrol/CDCS)                                  | [CDCS.jl](https://github.com/oxfordcontrol/CDCS.jl)                              | Manualᴹ | GPL     | LP, SOCP, SDP             |
| [CDD](https://github.com/cddlib/cddlib)                                        | [CDDLib.jl](https://github.com/JuliaPolyhedra/CDDLib.jl)                         |        | GPL      | LP                        |
| [Clp](https://github.com/coin-or/Clp)                                          | [Clp.jl](https://github.com/jump-dev/Clp.jl)                                     |        | EPL      | LP                        |
| [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl)                          |                                                                                  |        | Apache   | LP, QP, SOCP, SDP         |
| [Couenne](http://github.com/coin-or/Couenne)                                   | [AmplNLWriter.jl](https://github.com/jump-dev/AmplNLWriter.jl)                   |        | EPL      | (MI)NLP                   |
| [CPLEX](https://www.ibm.com/analytics/cplex-optimizer/)                        | [CPLEX.jl](https://github.com/jump-dev/CPLEX.jl)                                 | Manual | Comm.    | (MI)LP, (MI)SOCP          |
| [CSDP](https://github.com/coin-or/Csdp)                                        | [CSDP.jl](https://github.com/jump-dev/CSDP.jl)                                   |        | EPL      | LP, SDP                   |
| [EAGO.jl](https://github.com/psorlab/EAGO.jl)                                  |                                                                                  |        | MIT | NLP                    |
| [ECOS](https://github.com/ifa-ethz/ecos)                                       | [ECOS.jl](https://github.com/jump-dev/ECOS.jl)                                   |        | GPL      | LP, SOCP                  |
| [FICO Xpress](https://www.fico.com/en/products/fico-xpress-optimization-suite) | [Xpress.jl](https://github.com/jump-dev/Xpress.jl)                               | Manual | Comm.    | (MI)LP, (MI)SOCP          |
| [GLPK](http://www.gnu.org/software/glpk/)                                      | [GLPK.jl](https://github.com/jump-dev/GLPK.jl)                                   |        | GPL      | (MI)LP                    |
| [Gurobi](https://gurobi.com)                                                   | [Gurobi.jl](https://github.com/jump-dev/Gurobi.jl)                               | Manual | Comm.    | (MI)LP, (MI)SOCP          |
| [HiGHS](https://github.com/ERGO-Code/HiGHS)                                    | [HiGHS.jl](https://github.com/jump-dev/HiGHS.jl)                                 |        |MIT       | LP                        |
| [Hypatia.jl](https://github.com/chriscoey/Hypatia.jl)                          |                                                                                  |        | MIT      | LP, SOCP, SDP             |
| [Ipopt](https://github.com/coin-or/Ipopt)                                      | [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl)                                 |        | EPL      | LP, QP, NLP               |
| [Juniper.jl](https://github.com/lanl-ansi/Juniper.jl)                          |                                                                                  |        | MIT      | (MI)SOCP, (MI)NLP         |
| [MadNLP.jl](https://github.com/sshin23/MadNLP.jl)                              |                                                                                  |        | MIT      | LP, QP, NLP               |
| [MOSEK](https://www.mosek.com/)                                                | [MosekTools.jl](https://github.com/jump-dev/MosekTools.jl)                       | Manual | Comm.    | (MI)LP, (MI)SOCP, SDP     |
| [NLopt](https://github.com/stevengj/nlopt)                                     | [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl)                                 |        | GPL      | LP, QP, NLP               |
| [OSQP](https://osqp.org/)                                                      | [OSQP.jl](https://github.com/oxfordcontrol/OSQP.jl)                              |        | Apache   | LP, QP                    |
| [PATH](http://pages.cs.wisc.edu/~ferris/path.html)                             | [PATHSolver.jl](https://github.com/chkwon/PATHSolver.jl)                 |        | MIT      | MCP                       |
| [Pavito.jl](https://github.com/jump-dev/Pavito.jl)                             |                                                                                  |        | MPL-2    | (MI)NLP                   |
| [ProxSDP.jl](https://github.com/mariohsouto/ProxSDP.jl)                        |                                                                                  |        | MIT      | LP, SOCP, SDP             |
| [SCIP](https://scipopt.org/)                                                   | [SCIP.jl](https://github.com/scipopt/SCIP.jl)                            |        | ZIB      | (MI)LP, (MI)NLP           |
| [SCS](https://github.com/cvxgrp/scs)                                           | [SCS.jl](https://github.com/jump-dev/SCS.jl)                                     |        | MIT      | LP, SOCP, SDP             |
| [SDPA](http://sdpa.sourceforge.net/)                                           | [SDPA.jl](https://github.com/jump-dev/SDPA.jl), [SDPAFamily.jl](https://github.com/ericphanson/SDPAFamily.jl) |  | GPL | LP, SDP |
| [SDPNAL](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/)               | [SDPNAL.jl](https://github.com/jump-dev/SDPNAL.jl)                               | Manualᴹ | CC BY-SA | LP, SDP                  |
| [SDPT3](https://blog.nus.edu.sg/mattohkc/softwares/sdpt3/)                     | [SDPT3.jl](https://github.com/jump-dev/SDPT3.jl)                                 | Manualᴹ | GPL      | LP, SOCP, SDP            |
| [SeDuMi](http://sedumi.ie.lehigh.edu/)                                         | [SeDuMi.jl](https://github.com/jump-dev/SeDuMi.jl)                               | Manualᴹ | GPL      | LP, SOCP, SDP            |
| [Tulip.jl](https://github.com/ds4dm/Tulip.jl)                                  |                                                                                  |        | MPL-2     | LP                       |

Where:
- LP = Linear programming
- QP = Quadratic programming
- SOCP = Second-order conic programming (including problems with convex
  quadratic constraints or objective)
- MCP = Mixed-complementarity programming
- NLP = Nonlinear programming
- SDP = Semidefinite programming
- (MI)XXX = Mixed-integer equivalent of problem type `XXX`

!!! note
    Developed a solver or solver wrapper? This table is open for new
    contributions! Start by making a pull request to edit the [installation.md](https://github.com/jump-dev/JuMP.jl/blob/master/docs/src/installation.md)
    file.

!!! note
    Developing a solver or solver wrapper? See [Models](@ref jump_models) and the
    [MathOptInterface docs](https://jump.dev/MathOptInterface.jl/stable/) for
    more details on how JuMP interacts with solvers. Please get in touch via the
    [Developer Chatroom](https://jump.dev/pages/governance/#developer-chatroom)
    with any questions about connecting new solvers with JuMP.

## AMPL-based solvers

Use [AmplNLWriter](https://github.com/jump-dev/AmplNLWriter.jl) to access
solvers that support the [nl format](https://en.wikipedia.org/wiki/Nl_(format)).

Some solvers, such as [Bonmin](https://github.com/coin-or/Bonmin) and
[Couenne](https://github.com/coin-or/Couenne) can be installed via the Julia
package manager. Others need to be manually installed.

Consult the AMPL documentation for a [complete list of supported solvers](https://ampl.com/products/solvers/all-solvers-for-ampl/).

## GAMS-based solvers

Use [GAMS.jl](https://github.com/GAMS-dev/gams.jl) to access solvers available
through [GAMS](https://www.gams.com). Such solvers include:
[AlphaECP](https://www.gams.com/latest/docs/S_ALPHAECP.html),
[Antigone](https://www.gams.com/latest/docs/S_ANTIGONE.html),
[BARON](https://www.gams.com/latest/docs/S_BARON.html),
[CONOPT](https://www.gams.com/latest/docs/S_CONOPT.html),
[Couenne](https://www.gams.com/latest/docs/S_COUENNE.html),
[LocalSolver](https://www.gams.com/latest/docs/S_LOCALSOLVER.html),
[PATHNLP](https://www.gams.com/latest/docs/S_PATHNLP.html),
[SHOT](https://www.gams.com/latest/docs/S_SHOT.html),
[SNOPT](https://www.gams.com/latest/docs/S_SNOPT.html),
[SoPlex](https://www.gams.com/latest/docs/S_SOPLEX.html).
See a complete list [here](https://www.gams.com/latest/docs/S_MAIN.html).

!!! note
    [GAMS.jl](https://github.com/GAMS-dev/gams.jl) requires an installation of
    the commercial software [GAMS](https://www.gams.com) for which a
    [free community license](https://www.gams.com/latest/docs/UG_License.html#GAMS_Community_Licenses)
    exists.

## NEOS-based solvers

Use [NEOSServer.jl](https://github.com/odow/NEOSServer.jl) to access solvers
available through the [NEOS Server](https://neos-server.org).

## Previously supported solvers

The following solvers were compatible with JuMP up to release 0.18 but are
not yet compatible with the latest version because they do not implement the
new MathOptInterface API:

- [Pajarito](https://github.com/JuliaOpt/Pajarito.jl)

Please join the [Developer Chatroom](https://jump.dev/pages/governance/#developer-chatroom)
if you have interest in reviving a previously supported solver.

## Common installation issues

!!! tip
    When in doubt, run `import Pkg; Pkg.update()` to see if updating your
    packages fixes the issue. Remember you will need to exit Julia and start a
    new session for the changes to take effect.


### Check the version of your packages

Each package is versioned with a [three-part number](https://semver.org) of the
form `vX.Y.Z`. You can check which versions you have installed with
`import Pkg; Pkg.status()`.

This should almost always be the most-recent release. You can check the releases
of a package by going to the relevant GitHub page, and navigating to the
"releases" page. For example, the list of JuMP releases is available at:
[https://github.com/jump-dev/JuMP.jl/releases](https://github.com/jump-dev/JuMP.jl/releases).

If you post on the [community forum](https://discourse.julialang.org/c/domain/opt/13),
please include the output of `Pkg.status()`!

### Unsatisfiable requirements detected

Did you get an error like `Unsatisfiable requirements detected for package JuMP`?
The Pkg documentation has a [section on how to understand and manage these conflicts](https://julialang.github.io/Pkg.jl/v1/managing-packages/#conflicts).

### Installing new packages can make JuMP downgrade to an earlier version

Another common complaint is that after adding a new package, code that
previously worked no longer works.

This usually happens because the new package is not compatible with the latest
version of JuMP. Therefore, the package manager rolls-back JuMP to an earlier
version! Here's an example.

First, we add JuMP:
```julia
(jump_example) pkg> add JuMP
  Resolving package versions...
Updating `~/jump_example/Project.toml`
  [4076af6c] + JuMP v0.21.5
Updating `~/jump_example/Manifest.toml`
  ... lines omitted ...
```
The `+ JuMP v0.21.5` line indicates that JuMP has been added at version
`0.21.5`. However, watch what happens when we add [JuMPeR](https://github.com/iainnz/JuMPeR.jl):
```julia
(jump_example) pkg> add JuMPeR
  Resolving package versions...
Updating `~/jump_example/Project.toml`
  [4076af6c] ↓ JuMP v0.21.5 ⇒ v0.18.6
  [707a9f91] + JuMPeR v0.6.0
Updating `~/jump_example/Manifest.toml`
  ... lines omitted ...
```
JuMPeR gets added at version `0.6.0` (`+ JuMPeR v0.6.0`), but JuMP gets
downgraded from `0.21.5` to `0.18.6` (`↓ JuMP v0.21.5 ⇒ v0.18.6`)! The reason
for this is that JuMPeR doesn't support a version of JuMP newer than `0.18.6`.

!!! tip
    Pay careful attention to the output of the package manager when adding new
    packages, especially when you see a package being downgraded!
# Should I use JuMP?

JuMP is an [algebraic modeling language](@ref algebraic-modeling-language) for
mathematical optimization written in the [Julia language](https://julialang.org).

This page explains when you should consider using JuMP, and importantly, when
you should _not_ use JuMP.

## When should I use JuMP?

You should use JuMP if you have a constrained optimization problem for which you
can formulate:
 * a set of decision variables
 * a scalar objective function
 * a set of constraints.

Key reasons to use JuMP include:

 - User friendliness
   - JuMP has syntax that mimics natural mathematical expressions. (See the
     section on [algebraic modeling languages](@ref algebraic-modeling-language).)
 - Speed
   - Benchmarking has shown that JuMP can create problems at similar speeds to
     special-purpose modeling languages such as [AMPL](https://ampl.com/).
   - JuMP communicates with most solvers in memory, avoiding the need to write
     intermediary files.
 - Solver independence
   - JuMP uses a generic solver-independent interface provided by the
     [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl)
     package, making it easy to change between a number of open-source and
     commercial optimization software packages ("solvers"). The
     [Supported solvers](@ref) section contains a table of the currently
     supported solvers.
 - Access to advanced algorithmic techniques
   - JuMP supports efficient _in-memory_ re-solves of linear programs, which
     previously required using solver-specific or low-level C++ libraries.
   - JuMP provides access to solver-independent and solver-dependent
     [Callbacks](@ref callbacks_manual).
 - Ease of embedding
   - JuMP itself is written purely in Julia. Solvers are the only binary
     dependencies.
   - Automated install of many solver dependencies.
     - JuMP provides automatic installation of many open-source solvers. This is
       different to modeling languages in Python which require you to download
       and install a solver yourself.
   - Being embedded in a general-purpose programming language makes it easy to
     solve optimization problems as part of a larger workflow (e.g., inside a
     simulation, behind a web server, or as a subproblem in a decomposition
     algorithm).
     - As a trade-off, JuMP's syntax is constrained by the syntax available in
       Julia.
   - JuMP is [MPL](https://www.mozilla.org/MPL/2.0/) licensed, meaning that it
     can be embedded in commercial software that complies with the terms of the
     license.

## When should I not use JuMP?

JuMP supports a broad range of optimization classes. However, there are still
some that it doesn't support, or that are better supported by other software
packages.

### You want to optimize a complicated Julia function

Packages in Julia compose well. It's common for people to pick two unrelated
packages and use them in conjunction to create novel behavior. JuMP isn't one of
those packages.

If you want to optimize a ordinary differential equation from
[DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)
or tune a neural network from [Flux.jl](https://github.com/FluxML/Flux.jl),
consider using other packages such as:
 * [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl)
 * [GalacticOptim.jl](https://github.com/SciML/GalacticOptim.jl)
 * [Nonconvex.jl](https://github.com/JuliaNonconvex/Nonconvex.jl)

### Black-box, derivative free, or unconstrained optimization

JuMP does support nonlinear programs with constraints and objectives containing
user-defined functions. However, the functions must be automatically
differentiable, or need to provide explicit derivatives. (See
[User-defined Functions](@ref) for more information.)

If your function is a black-box that is non-differentiable (e.g., the output of
a simulation written in C++), JuMP is not the right tool for the job. This also
applies if you want to use a derivative free method.

Even if your problem is differentiable, if it is unconstrained there is limited
benefit (and downsides in the form of more overhead) to using JuMP over tools
which are only concerned with function minimization.

Alternatives to consider are:
 * [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl)
 * [GalacticOptim.jl](https://github.com/SciML/GalacticOptim.jl)
 * [NLopt.jl](https://github.com/JuliaOpt/NLopt.jl)

### Multiobjective programs

If your problem has more than one objective, JuMP is not the right tool for the
job. However, [we're working on fixing this!](https://github.com/jump-dev/JuMP.jl/issues/2099).

Alternatives to consider are:
 * [vOptGeneric.jl](https://github.com/vOptSolver/vOptGeneric.jl)

### Disciplined convex programming

JuMP does not support [disciplined convex programming (DCP)](https://dcp.stanford.edu).

Alternatives to consider are:
 * [Convex.jl](https://github.com/jump-dev/Convex.jl)

!!! note
    `Convex.jl` is also built on MathOptInterface, and shares the same set of
    underlying solvers. However, you input problems differently, and Convex.jl
    checks that the problem is DCP.

### Stochastic programming

JuMP requires deterministic input data.

If you have stochastic input data, consider using a JuMP extension such as:
 * [InfiniteOpt.jl](https://github.com/pulsipher/InfiniteOpt.jl)
 * [StochasticPrograms.jl](https://github.com/martinbiel/StochasticPrograms.jl)
 * [SDDP.jl](https://github.com/odow/SDDP.jl)
# Release notes

## Version 0.22.2 (January 10, 2021)

- New features:
  - The function `all_nl_constraints` now returns all nonlinear constraints
    in a model
  - `start_value` and `set_start_value` can now be used to get and set the
    primal start for constraint references
  - Plural macros now return a tuple containing the elements that were defined
    instead of `nothing`
  - Anonymous variables are now printed as `_[i]` where `i` is the index of the
    variable instead of `noname`. Calling `name(x)` still returns `""` so this
    is non-breaking.
- Bug fixes:
  - Fixed handling of `min` and `max` in nonlinear expressions
  - CartesianIndex is no longer allowed as a key for DenseAxisArrays.
- Documentation, maintenance:
  - Improved the performance of GenericAffExpr
  - Added a tutorial on the Travelling Salesperson Problem
  - Added a tutorial on querying the Hessian of a nonlinear program
  - Added documentation on using custom solver binaries.

## Version 0.22.1 (November 29, 2021)

- New features:
 * Export `OptimizationSense` enum, with instances: `MIN_SENSE`, `MAX_SENSE`,
   and `FEASIBILITY_SENSE`
 * Add `Base.isempty(::Model)` to match `Base.empty(::Model)`
- Bug fixes:
 * Fix bug in container with tuples as indices
 * Fix bug in `set_time_limit_sec`
- Documentation, maintenance
 * Add tutorial "Design patterns for larger models"
 * Remove release notes section from PDF
 * General edits of the documentation and error messages

## Version 0.22.0 (November 10, 2021)

**JuMP v0.22 is a breaking release**

### Breaking changes

JuMP 0.22 contains a number of breaking changes. However, these should be
invisible for the majority of users. You will mostly encounter these breaking
changes if you: wrote a JuMP extension, accessed `backend(model)`, or called
 `@SDconstraint`.

The breaking changes are as follows:

 * MathOptInterface has been updated to v0.10.4. For users who have interacted
   with the MOI backend, this contains a large number of breaking changes. Read
   the [MathOptInterface release notes](https://jump.dev/MathOptInterface.jl/v0.10/release_notes/#v0.10.0-(September-6,-2021))
   for more details.
 * The `bridge_constraints` keyword argument to `Model` and `set_optimizer` has
   been renamed `add_bridges` to reflect that more thing were bridged than just
   constraints.
 * The `backend(model)` field now contains a concrete instance of a
   `MOI.Utilities.CachingOptimizer` instead of one with an abstractly typed
   optimizer field. In most cases, this will lead to improved performance.
   However, calling `set_optimizer` after `backend` invalidates the old
   backend. For example:
   ```julia
   model = Model()
   b = backend(model)
   set_optimizer(model, GLPK.Optimizer)
   @variable(model, x)
   # b is not updated with `x`! Get a new b by calling `backend` again.
   new_b = backend(model)
   ```
 * All usages of `@SDconstraint` are deprecated. The new syntax is
   `@constraint(model, X >= Y, PSDCone())`.
 * Creating a `DenseAxisArray` with a `Number` as an axis will now display a
   warning. This catches a common error in which users write
   `@variable(model, x[length(S)])` instead of
   `@variable(model, x[1:length(S)])`.
 * The `caching_mode` argument to `Model`, e.g.,
   `Model(caching_mode = MOIU.MANUAL)` mode has been removed. For more control
   over the optimizer, use `direct_model` instead.
 * The previously deprecated `lp_objective_perturbation_range` and
   `lp_rhs_perturbation_range` functions have been removed. Use
   `lp_sensitivity_report` instead.
 * The `.m` fields of `NonlinearExpression` and `NonlinearParameter` have been
   renamed to `.model`.
 * Infinite variable bounds are now ignored. Thus, `@variable(model, x <= Inf)`
   will show `has_upper_bound(x) == false`. Previously, these bounds were passed
   through to the solvers which caused numerical issues for solvers expecting
   finite bounds.
 * The `variable_type` and `constraint_type` functions were removed. This should
   only affect users who previously wrote JuMP extensions. The functions can be
   deleted without consequence.
 * The internal functions `moi_mode`, `moi_bridge_constraints`,
   `moi_add_constraint`, and `moi_add_to_function_constant` are no longer
   exported.
 * The un-used method `Containers.generate_container` has been deleted.
 * The `Containers` API has been refactored, and `_build_ref_sets` is now
   public as `Containers.build_ref_sets`.
 * The `parse_constraint_` methods for extending `@constraint` at parse time
   have been refactored in a breaking way. Consult the Extensions documentation
   for more details and examples.

### New features

 * The `TerminationStatusCode` and `ResultStatusCode` enums are now exported
   by JuMP. Prefer `termination_status(model) == OPTIMAL` instead of
   `== MOI.OPTIMAL`, although the `MOI.` prefix way still works.
 * Copy a `x::DenseAxisArray` to an `Array` by calling `Array(x)`.
 * `NonlinearExpression` is now a subtype of `AbstractJuMPScalar`
 * Constraints such as `@constraint(model, x + 1 in MOI.Integer())` are now
   supported.
 * `primal_feasibility_report` now accepts a function as the first argument.
 * Scalar variables `@variable(model, x[1:2] in MOI.Integer())` creates two
   variables, both of which are constrained to be in the set `MOI.Integer`.
 * Conic constraints can now be specified as inequalities under a different
   partial ordering. So `@constraint(model, x - y in MOI.Nonnegatives())` can
   now be written as `@constraint(model, x >= y, MOI.Nonnegatives())`.
 * Names are now set for vectorized constraints.

### Documentation, maintenance and performance

 * The documentation is now available as a PDF.
 * The documentation now includes a full copy of the MathOptInterface
   documentation to make it easy to link concepts between the docs. (The
   MathOptInterface documentation has also been significantly improved.)
 * The documentation contains a large number of improvements and clarifications
   on a range of topics. Thanks to @sshin23, @DilumAluthge, and @jlwether.
 * The documentation is now built with Julia 1.6 instead of 1.0.
 * Various error messages have been improved to be more readable.
 * Fixed a performance issue when `show` was called on a `SparseAxisArray` with
   a large number of elements.
 * Fixed a bug displaying barrier and simplex iterations in `solution_summary`.
 * Fixed a bug by implementing `hash` for `DenseAxisArray` and
   `SparseAxisArray`.
 * Names are now only set if the solver supports them. Previously, this
   prevented solvers such as Ipopt from being used with `direct_model`.
 * `MutableArithmetics.Zero` is converted into a `0.0` before being returned to
   the user. Previously, some calls to `@expression` would return the
   undocumented `MutableArithmetics.Zero()` object. One example is summing over
   an empty set `@expression(model, sum(x[i] for i in 1:0))`. You will now get
   `0.0` instead.
 * `AffExpr` and `QuadExpr` can now be used with `== 0` instead of `iszero`.
   This fixes a number of issues relating to Julia standard libraries such as
   `LinearAlgebra` and `SparseArrays`.
 * Fixed a bug when registering a user-defined function with splatting.

## Version 0.21.10 (September 4, 2021)

For a detailed list of the closed issues and pull requests from this release,
see the [tag notes](https://github.com/jump-dev/JuMP.jl/releases/tag/v0.21.10).
A summary of changes are as follows:

- New features:
  * Add `add_NL_expression`
  * `add_NL_xxx` functions now support `AffExpr` and `QuadExpr` as terms
- Documentation, maintenance and performance:
  * Fix bug in `solution_summary`
  * Fix bug in `relax_integrality`
  * Improve error message in `lp_sensitivity_report`

## Version 0.21.9 (August 1, 2021)

For a detailed list of the closed issues and pull requests from this release,
see the [tag notes](https://github.com/jump-dev/JuMP.jl/releases/tag/v0.21.9).
A summary of changes are as follows:

- New features:
  * Containers now support arbitrary container types by passing the type to the
    `container` keyword and overloading `Containers.container`.
  * `is_valid` now supports nonlinear constraints
  * Added `unsafe_backend` for querying the inner-most optimizer of a JuMP
    model.
  * Nonlinear parameters now support the plural `@NLparameters` macro.
  * Containers (e.g., `DenseAxisArray`) can now be used in vector-valued
    constraints.
- Documentation, maintenance and performance:
  * Various improvements to the documentation.

## Version 0.21.8 (May 8, 2021)

For a detailed list of the closed issues and pull requests from this release,
see the [tag notes](https://github.com/jump-dev/JuMP.jl/releases/tag/v0.21.8).
A summary of changes are as follows:

- New features:
  * The `@constraint` macro is now extendable in the same way as `@variable`.
  * `AffExpr` and `QuadExpr` can now be used in nonlinear macros.
- Bug fixes:
  * Fixed a bug in `lp_sensitivity_report`.
  * Fixed an inference issue when creating empty `SparseAxisArray`s.

## Version 0.21.7 (April 12, 2021)

For a detailed list of the closed issues and pull requests from this release,
see the [tag notes](https://github.com/jump-dev/JuMP.jl/releases/tag/v0.21.7).
A summary of changes are as follows:

- New features:
  * Added `primal_feasibility_report`, which can be used to check whether a
    primal point satisfies primal feasibility.
  * Added `coefficient`, which returns the coefficient associated with a
    variable in affine and quadratic expressions.
  * Added `copy_conflict`, which returns the IIS of an infeasible model.
  * Added `solution_summary`, which returns (and prints) a struct containing a
    summary of the solution.
  * Allow `AbstractVector` in vector constraints instead of just `Vector`.
  * Added `latex_formulation(model)` which returns an object representing the
    latex formulation of a model. Use `print(latex_formulation(model))` to print
    the formulation as a string.
  * User-defined functions in nonlinear expressions are now automatically
    registered to aid quick model prototyping. However, a warning is printed to
    encourage the manual registration.
  * DenseAxisArray's now support broadcasting over multiple arrays.
  * Container indices can now be iterators of `Base.SizeUnknown`.
- Bug fixes:
  * Fixed bug in `rad2deg` and `deg2rad` in nonlinear expressions.
  * Fixed a MethodError bug in `Containers` when forcing container type.
  * Allow partial slicing of a DenseAxisArray, resolving an issue from 2014!
  * Fixed a bug printing variable names in IJulia.
  * Ending an IJulia cell with `model` now prints a summary of the model (like
    in the REPL) not the latex formulation. Use `print(model)` to print the latex
    formulation.
  * Fixed a bug when copying models containing nested arrays.
- Documentation, performance improvements, and general maintenance:
  * Tutorials are now part of the documentation, and more refactoring has taken
    place.
  * Added JuliaFormatter added as a code formatter.
  * Added some precompilation statements to reduce initial latency.
  * Various improvements to error messages to make them more helpful.
  * Improved performance of `value(::NonlinearExpression)`.
  * Improved performance of `fix(::VariableRef)`.

## Version 0.21.6 (January 29, 2021)

For a detailed list of the closed issues and pull requests from this release,
see the [tag notes](https://github.com/jump-dev/JuMP.jl/releases/tag/v0.21.6).
A summary of changes are as follows:

- New features:
  * Added support for skew symmetric variables via
    `@variable(model, X[1:2, 1:2] in SkewSymmetricMatrixSpace())`.
  * `lp_sensitivity_report` has been added which significantly improves the
    performance of querying the sensitivity summary of an LP.
    `lp_objective_perturbation_range` and `lp_rhs_perturbation_range` are
    deprecated.
  * Dual warm-starts are now supported with `set_dual_start_value` and
    `dual_start_value`.
  * `∈` (`\in<tab>`) can now be used in macros instead of `=` or `in`.
  * Use `haskey(model::Model, key::Symbol)` to check if a name `key` is
    registered in a model.
  * Added `unregister(model::Model, key::Symbol)` to unregister a name `key`
    from `model`.
  * Added `callback_node_status` for use in callbacks.
  * Added `print_bridge_graph` to visualize the bridging graph generated by
    MathOptInterface.
  * Improved error message for containers with duplicate indices.
- Bug fixes:
  * Various fixes to pass tests on Julia 1.6.
  * Fixed a bug in the printing of nonlinear expressions in IJulia.
  * Fixed a bug when nonlinear expressions are passed to user-defined functions.
  * Some internal functions that were previously exported are now no longer
    exported.
  * Fixed a bug when relaxing a fixed binary variable.
  * Fixed a `StackOverflowError` that occurred when `SparseAxisArray`s had a
    large number of elements.
  * Removed an unnecessary type assertion in `list_of_constraint_types`.
  * Fixed a bug when copying models with registered expressions.
- Documentation and general maintenance:
  * The documentation has been significantly overhauled. It now has distinct
    sections for the manual, API reference, and examples. The existing examples
    in `/examples` have now been moved to `/docs/src/examples` and rewritten
    using `Literate.jl`, and they are now included in the documentation.
  * JuliaFormatter has been applied to most of the codebase. This will continue
    to roll out over time, as we fix upstream issues in the formatter, and will
    eventually become compulsory.
  * The root cause of a large number of method invalidations has been resolved.
  * We switched continuous integration from Travis and Appveyor to GitHub
    Actions.

## Version 0.21.5 (September 18, 2020)

For a detailed list of the closed issues and pull requests from this release,
see the [tag notes](https://github.com/jump-dev/JuMP.jl/releases/tag/v0.21.5).
A summary of changes are as follows:

- Fix deprecation warnings
- Throw `DimensionMismatch` for incompatibly sized functions and sets
- Unify treatment of `keys(x)` on JuMP containers

## Version 0.21.4 (September 14, 2020)

For a detailed list of the closed issues and pull requests from this release,
see the [tag notes](https://github.com/jump-dev/JuMP.jl/releases/tag/v0.21.4).
A summary of changes are as follows:

- New features:
  * Add debug info when adding unsupported constraints
  * Add `relax_integrality` for solving continuous relaxation
  * Allow querying constraint conflicts
- Bug fixes:
  * Dispatch on `Real` for `MOI.submit`
  * Implement `copy` for `CustomSet` in tests
  * Don't export private macros
  * Fix invalid assertion in nonlinear
- Error if constraint has `NaN` right-hand side
- Improve speed of tests
  * Lots of work modularizing files in `/test`
- Improve line numbers in macro error messages
- Print nonlinear subexpressions
- Various documentation updates
- Dependency updates:
  * Datastructures 0.18
  * MathOptFormat v0.5
  * Prep for MathOptInterface 0.9.15

## Version 0.21.3 (June 18, 2020)

- Added Special Order Sets (SOS1 and SOS2) to JuMP with default weights to ease
  the creation of such constraints (#2212).
- Added functions `simplex_iterations`, `barrier_iterations` and `node_count`
  (#2201).
- Added function `reduced_cost` (#2205).
- Implemented `callback_value` for affine and quadratic expressions (#2231).
- Support `MutableArithmetics.Zero` in objective and constraints (#2219).
- Documentation improvements:
  * Mention tutorials in the docs (#2223).
  * Update COIN-OR links (#2242).
  * Explicit link to the documentation of `MOI.FileFormats` (#2253).
  * Typo fixes (#2261).
- Containers improvements:
  * Fix `Base.map` for `DenseAxisArray` (#2235).
  * Throw `BoundsError` if number of indices is incorrect for `DenseAxisArray`
    and `SparseAxisArray` (#2240).
- Extensibility improvements:
  * Implement a `set_objective` method fallback that redirects to
    `set_objective_sense` and `set_objective_function` (#2247).
  * Add `parse_constraint` method with arbitrary number of arguments (#2051).
  * Add `parse_constraint_expr` and `parse_constraint_head` (#2228).

## Version 0.21.2 (April 2, 2020)

- Added `relative_gap()` to access `MOI.RelativeGap()` attribute (#2199).
- Documentation fixes:
  * Added link to source for docstrings in the documentation (#2207).
  * Added docstring for `@variables` macro (#2216).
  * Typo fixes (#2177, #2184, #2182).
- Implementation of methods for Base functions:
  * Implemented `Base.empty!` for `JuMP.Model` (#2198).
  * Implemented `Base.conj` for JuMP scalar types (#2209).
- Bug fixes:
  * Fixed sum of expression with scalar product in macro (#2178).
  * Fixed writing of nonlinear models to MathOptFormat (#2181).
  * Fixed construction of empty SparseAxisArray (#2179).
  * Fixed constraint with zero function (#2188).

## Version 0.21.1 (Feb 18, 2020)

- Improved the clarity of the `with_optimizer` deprecation warning.

## Version 0.21 (Feb 16, 2020)

Breaking changes:

- Deprecated `with_optimizer` (#2090, #2084, #2141). You can replace
  `with_optimizer` by either nothing, `optimizer_with_attributes` or a closure:
  * replace `with_optimizer(Ipopt.Optimizer)` by `Ipopt.Optimizer`.
  * replace `with_optimizer(Ipopt.Optimizer, max_cpu_time=60.0)`
    by `optimizer_with_attributes(Ipopt.Optimizer, "max_cpu_time" => 60.0)`.
  * replace `with_optimizer(Gurobi.Optimizer, env)` by `() -> Gurobi.Optimizer(env)`.
  * replace `with_optimizer(Gurobi.Optimizer, env, Presolve=0)`
    by `optimizer_with_attributes(() -> Gurobi.Optimizer(env), "Presolve" => 0)`.

  alternatively to `optimizer_with_attributes`, you can also set the attributes
  separately with `set_optimizer_attribute`.
- Renamed `set_parameter` and `set_parameters` to `set_optimizer_attribute` and
  `set_optimizer_attributes` (#2150).
- Broadcast should now be explicit inside macros. `@SDconstraint(model, x >= 1)`
  and `@constraint(model, x + 1 in SecondOrderCone())` now throw an error
  instead of broadcasting `1` along the dimension of `x` (#2107).
- `@SDconstraint(model, x >= 0)` is now equivalent to `@constraint(model, x in PSDCone())`
  instead of `@constraint(model, (x .- 0) in PSDCone())` (#2107).
- The macros now create the containers with `map` instead of `for` loops,
  as a consequence, containers created by `@expression` can now have any element
  type and containers of constraint references now have concrete element types
  when possible. This fixes a long-standing issue where `@expression` could
  only be used to generate a collection of linear expressions. Now it works for
  quadratic expressions as well (#2070).
- Calling `deepcopy(::AbstractModel)` now throws an error.
- The constraint name is now printed in the model string (#2108).

New features:

- Added support for solver-independent and solver-specific callbacks (#2101).
- Added `write_to_file` and `read_from_file`, supported formats are CBF, LP,
  MathOptFormat, MPS and SDPA (#2114).
- Added support for complementarity constraints (#2132).
- Added support for indicator constraints (#2092).
- Added support for querying multiple solutions with the `result` keyword (#2100).
- Added support for constraining variables on creation (#2128).
- Added method `delete` that deletes a vector of variables at once if it is
  supported by the underlying solver (#2135).
- The arithmetic between JuMP expression has be refactored into the
  MutableArithmetics package (#2107).
- Improved error on complex values in NLP (#1978).
- Added an example of column generation (#2010).

Bug fixes:

- Incorrect coefficients generated when using Symmetric variables (#2102)

## Version 0.20.1 (Oct 18, 2019)

- Add sections on `@variables` and `@constraints` in the documentation (#2062).
- Fixed product of sparse matrices for Julia v1.3 (#2063).
- Added `set_objective_coefficient` to modify the coefficient of a linear term
  of the objective function (#2008).
- Added `set_time_limit_sec`, `unset_time_limit_sec` and `time_limit_sec` to set
  and query the time limit for the solver in seconds (#2053).

## Version 0.20.0 (Aug 24, 2019)

- Documentation updates.
- Numerous bug fixes.
- Better error messages (#1977, #1978, #1997, #2017).
- Performance improvements (#1947, #2032).
- Added LP sensitivity summary functions `lp_objective_perturbation_range`
  and `lp_rhs_perturbation_range` (#1917).
- Added functions `dual_objective_value`, `raw_status` and `set_parameter`.
- Added function `set_objective_coefficient` to modify the coefficient of
  a linear term of the objective (#2008).
- Added functions `set_normalized_rhs`, `normalized_rhs`, and
  `add_to_function_constant` to modify and get the constant part
  of a constraint (#1935, #1960).
- Added functions `set_normalized_coefficient` and `normalized_coefficient`
  to modify and get the coefficient of a linear term of a constraint
  (#1935, #1960).
- Numerous other improvements in MOI 0.9, see the `NEWS.md` file of MOI for more
  details.

## Version 0.19.2 (June 8, 2019)

- Fix a bug in derivatives that could arise in models with nested nonlinear
  subexpressions.

## Version 0.19.1 (May 12, 2019)

- Usability and performance improvements.
- Bug fixes.

## Version 0.19.0 (February 15, 2019)

**JuMP 0.19 contains significant breaking changes.**

Breaking changes:

- JuMP's abstraction layer for communicating with solvers changed from
  [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl) (MPB) to
  [MathOptInterface](https://github.com/JuliaOpt/MathOptInterface.jl)
  (MOI). MOI addresses many longstanding design issues. (See @mlubin's
  [slides](https://www.juliaopt.org/meetings/bordeaux2018/lubin.pdf) from
  JuMP-dev 2018.) JuMP 0.19 is compatible only with solvers that have been
  updated for MOI. See the
  [installation guide](https://www.juliaopt.org/JuMP.jl/dev/installation/)
  for a list of solvers that have and have not yet been updated.

- Most solvers have been renamed to `PackageName.Optimizer`. For example,
  `GurobiSolver()` is now `Gurobi.Optimizer`.

- Solvers are no longer added to a model via `Model(solver = XXX(kwargs...))`.
  Instead use `Model(with_optimizer(XXX, kwargs...))`. For example, `Model(with_optimizer(Gurobi.Optimizer, OutputFlag=0))`.

- JuMP containers (e.g., the objects returned by `@variable`) have been
  redesigned. `Containers.SparseAxisArray` replaces `JuMPDict`, `JuMPArray` was
  rewritten (inspired by `AxisArrays`) and renamed `Containers.DenseAxisArray`,
  and you can now request a container type with the `container=` keyword to the
  macros. See the corresponding
  [documentation](https://www.juliaopt.org/JuMP.jl/dev/variables/#Variable-containers-1)
  for more details.

- The statuses returned by solvers have changed. See the possible status
  values
  [here](https://www.juliaopt.org/MathOptInterface.jl/stable/apireference.html#Termination-Status-1).
  The MOI statuses are much richer than the MPB statuses and can be used to
  distinguish between previously indistinguishable cases (e.g. did the solver
  have a feasible solution when it stopped because of the time limit?).

- Starting values are separate from result values. Use `value` to query
  the value of a variable in a solution. Use `start_value` and `set_start_value`
  to get and set an initial starting point provided to the solver. The solutions
  from previous solves are no longer automatically set as the starting points
  for the next solve.

- The data structures for affine and quadratic expressions `AffExpr` and
  `QuadExpr` have changed. Internally, terms are stored in dictionaries instead
  of lists. Duplicate coefficients can no longer exist. Accessors and iteration
  methods have changed.

- `JuMPNLPEvaluator` no longer includes the linear and quadratic parts of the
  model in the evaluation calls. These are now handled separately to allow NLP
  solvers that support various types of constraints.

- JuMP solver-independent callbacks have been replaced by solver-specific
  callbacks. See your favorite solver for more details. (See the note below: No
  solver-specific callbacks are implemented yet.)

- The `norm()` syntax is no longer recognized inside macros. Use the
  `SecondOrderCone()` set instead.

- JuMP no longer performs automatic transformation between special quadratic
  forms and second-order cone constraints. Support for these
  constraint classes depends on the solver.

- The symbols `:Min` and `:Max` are no longer used as optimization senses.
  Instead, JuMP uses the `OptimizationSense` enum from MathOptInterface.
  `@objective(model, Max, ...)`, `@objective(model, Min, ...)`,
  `@NLobjective(model, Max, ...)`, and `@objective(model, Min, ...)` remain
  valid, but `@objective(m, :Max, ...)` is no longer accepted.

- The sign conventions for duals has changed in some cases for consistency with
  conic duality (see the
  [documentation](https://www.juliaopt.org/MathOptInterface.jl/v0.6.2/apimanual.html#Duals-1)).
  The `shadow_price` helper method returns duals with signs that match
  conventional LP interpretations of dual values as sensitivities of the
  objective value to relaxations of constraints.

- `@constraintref` is no longer defined. Instead, create the appropriate
  container to hold constraint references manually. For example,
  ```julia
  constraints = Dict() # Optionally, specify types for improved performance.
  for i in 1:N
    constraints[i] = @constraint(model, ...)
  end
  ```

- The `lowerbound`, `upperbound`, and `basename` keyword arguments to the `@variable`
  macro have been renamed to `lower_bound`, `upper_bound`, and `base_name`,
  for consistency with JuMP's new
  [style recommendations](https://www.juliaopt.org/JuMP.jl/dev/style/).

- We rely on broadcasting syntax to apply accessors to collections of
  variables, e.g., `value.(x)` instead of `getvalue(x)` for collections. (Use
  `value(x)` when `x` is a scalar object.)

New features:

- Splatting (like `f(x...)`) is recognized in restricted settings in nonlinear
  expressions.

- Support for deleting constraints and variables.

- The documentation has been completely rewritten using docstrings and
  Documenter.

- Support for modeling mixed conic and quadratic models (e.g., conic models
  with quadratic objectives and bi-linear matrix inequalities).

- Significantly improved support for modeling new types of constraints and for
  extending JuMP's macros.

- Support for providing dual warm starts.

- Improved support for accessing solver-specific attributes (e.g., the
  irreducible inconsistent subsystem).

- Explicit control of whether symmetry-enforcing constraints are added to PSD
  constraints.

- Support for modeling exponential cones.

- Significant improvements in internal code quality and testing.

- Style and naming guidelines.

- Direct mode and manual mode provide explicit control over when copies of a
  model are stored or regenerated. See the corresponding
  [documentation](https://www.juliaopt.org/JuMP.jl/dev/solvers/).

There are known regressions from JuMP 0.18 that will be addressed in a future
release (0.19.x or later):

- Performance regressions in model generation
  ([issue](https://github.com/JuliaOpt/JuMP.jl/issues/1403)). Please file an
  issue anyway if you notice a significant performance regression. We have
  plans to address a number of performance issues, but we might not be aware of
  all of them.

- Fast incremental NLP solves are not yet reimplemented
  ([issue](https://github.com/JuliaOpt/JuMP.jl/issues/1185)).

- We do not yet have an implementation of solver-specific callbacks.

- The column generation syntax in `@variable` has been removed (i.e., the
  `objective`, `coefficients`, and `inconstraints` keyword arguments). Support
  for column generation will be re-introduced in a future release.

- The ability to solve the continuous relaxation (i.e. via
  `solve(model; relaxation = true)`) is not yet reimplemented ([issue](https://github.com/JuliaOpt/JuMP.jl/issues/1611)).

## Version 0.18.5 (December 1, 2018)

   * Support views in some derivative evaluation functions.
   * Improved compatibility with PackageCompiler.

## Version 0.18.4 (October 8, 2018)

   * Fix a bug in model printing on Julia 0.7 and 1.0.

## Version 0.18.3 (October 1, 2018)

   * Add support for Julia v1.0 (Thanks @ExpandingMan)
   * Fix matrix expressions with quadratic functions (#1508)

## Version 0.18.2 (June 10, 2018)

   * Fix a bug in second-order derivatives when expressions are present (#1319)
   * Fix a bug in `@constraintref` (#1330)

## Version 0.18.1 (April 9, 2018)

   * Fix for nested tuple destructuring (#1193)
   * Preserve internal model when relaxation=true (#1209)
   * Minor bug fixes and updates for example

## Version 0.18.0 (July 27, 2017)

   * Drop support for Julia 0.5.
   * Update for ForwardDiff 0.5.
   * Minor bug fixes.

## Version 0.17.1 (June 9, 2017)

   * Use of `constructconstraint!` in `@SDconstraint`.
   * Minor bug fixes.

## Version 0.17.0 (May 27, 2017)

   * **Breaking change**: Mixing quadratic and conic constraints is no longer supported.
   * **Breaking change**: The ``getvariable`` and ``getconstraint`` functions are replaced by indexing on the corresponding symbol. For instance, to access the variable with name ``x``, one should now write ``m[:x]`` instead of ``getvariable(m, :x)``. As a consequence, creating a variable and constraint with the same name now triggers a warning, and accessing one of them afterwards throws an error. This change is breaking only in the latter case.
   * Addition of the ``getobjectivebound`` function that mirrors the functionality of the MathProgBase ``getobjbound`` function except that it takes into account transformations performed by JuMP.
   * Minor bug fixes.

The following changes are primarily of interest to developers of JuMP extensions:

   * The new syntax ``@constraint(model, expr in Cone)`` creates the constraint ensuring that ``expr`` is inside ``Cone``. The ``Cone`` argument is passed to ``constructconstraint!`` which enables the call to the dispatched to an extension.
   * The ``@variable`` macro now calls ``constructvariable!`` instead of directly calling the ``Variable`` constructor. Extra arguments and keyword arguments passed to ``@variable`` are passed to ``constructvariable!`` which enables the call to be dispatched to an extension.
   * Refactor the internal function ``conicdata`` (used build the MathProgBase conic model) into smaller subfunctions to make these parts reusable by extensions.

## Version 0.16.2 (March 28, 2017)

   * Minor bug fixes and printing tweaks
   * Address deprecation warnings for Julia 0.6

## Version 0.16.1 (March 7, 2017)

   * Better support for ``AbstractArray`` in JuMP (Thanks @tkoolen)
   * Minor bug fixes

## Version 0.16.0 (February 23, 2017)

   * **Breaking change**: JuMP no longer has a mechanism for selecting solvers by default (the previous mechanism was flawed and incompatible with Julia 0.6). Not specifying a solver before calling ``solve()`` will result in an error.
   * **Breaking change**: User-defined functions are no longer global. The first argument to ``JuMP.register`` is now a JuMP ``Model`` object within whose scope the function will be registered. Calling ``JuMP.register`` without a ``Model`` now produces an error.
   * **Breaking change**: Use the new ``JuMP.fix`` method to fix a variable to a value or to update the value to which a variable is fixed. Calling ``setvalue`` on a fixed variable now results in an error in order to avoid silent behavior changes. (Thanks @joaquimg)
   * Nonlinear expressions now print out similarly to linear/quadratic expressions (useful for debugging!)
   * New ``category`` keyword to ``@variable``. Used for specifying categories of anonymous variables.
   * Compatibility with Julia 0.6-dev.
   * Minor fixes and improvements (Thanks @cossio, @ccoffrin, @blegat)

## Version 0.15.1 (January 31, 2017)

  * Bugfix for ``@LinearConstraints`` and friends

## Version 0.15.0 (December 22, 2016)

  * Julia 0.5.0 is the minimum required version for this release.
  * Document support for BARON solver
  * Enable info callbacks in more states than before, e.g. for recording solutions.
    New ``when`` argument to ``addinfocallback`` ([#814](https://github.com/JuliaOpt/JuMP.jl/pull/814), thanks @yeesian)
  * Improved support for anonymous variables. This includes new warnings for potentially confusing use of the traditional non-anonymous syntax:
    * When multiple variables in a model are given the same name
    * When non-symbols are used as names, e.g., ``@variable(m, x[1][1:N])``
  * Improvements in iterating over JuMP containers ([#836](https://github.com/JuliaOpt/JuMP.jl/pull/836), thanks @IssamT)
  * Support for writing variable names in .lp file output (Thanks @leethargo)
  * Support for querying duals to SDP problems (Thanks @blegat)
  * The comprehension syntax with curly braces ``sum{}``, ``prod{}``, and ``norm2{}`` has been deprecated
    in favor of Julia's native comprehension syntax ``sum()``, ``prod()`` and ``norm()`` as previously announced.
    (For early adopters of the new syntax, ``norm2()`` was renamed to ``norm()`` without deprecation.)
  * Unit tests rewritten to use Base.Test instead of FactCheck
  * Improved support for operations with matrices of JuMP types (Thanks @ExpandingMan)
  * The syntax to halt a solver from inside a callback has changed from ``throw(CallbackAbort())`` to ``return JuMP.StopTheSolver``
  * Minor bug fixes

## Version 0.14.2 (December 12, 2016)

  * Allow singeton anonymous variables (includes bugfix)

## Version 0.14.1 (September 12, 2016)

  * More consistent handling of states in informational callbacks,
    includes a new ``when`` parameter to ``addinfocallback`` for
    specifying in which state an informational callback should be called.

## Version 0.14.0 (August 7, 2016)

  * Compatibility with Julia 0.5 and ForwardDiff 0.2
  * Support for "anonymous" variables, constraints, expressions, and parameters, e.g.,
    ``x = @variable(m, [1:N])`` instead of ``@variable(m, x[1:N])``
  * Support for retrieving constraints from a model by name via ``getconstraint``
  * ``@NLconstraint`` now returns constraint references (as expected).
  * Support for vectorized expressions within lazy constraints
  * On Julia 0.5, parse new comprehension syntax ``sum(x[i] for i in 1:N if isodd(i))``
    instead of ``sum{ x[i], i in 1:N; isodd(i) }``. The old syntax with curly
    braces will be deprecated in JuMP 0.15.
  * Now possible to provide nonlinear expressions as "raw" Julia ``Expr`` objects
    instead of using JuMP's nonlinear macros. This input format is useful for
    programmatically generated expressions.
  * ``s/Mathematical Programming/Mathematical Optimization/``
  * Support for local cuts (Thanks to @madanim, Mehdi Madani)
  * Document Xpress interface developed by @joaquimg, Joaquim Dias Garcia
  * Minor bug and deprecation fixes (Thanks @odow, @jrevels)


## Version 0.13.2 (May 16, 2016)

  * Compatibility update for MathProgBase

## Version 0.13.1 (May 3, 2016)

  * Fix broken deprecation for ``registerNLfunction``.

## Version 0.13.0 (April 29, 2016)

  * Most exported methods and macros have been renamed to avoid camelCase. See the list of changes [here](https://github.com/JuliaOpt/JuMP.jl/blob/e53d0db67cde2a4b80d0c1281f4b49eb0128a1f5/src/deprecated.jl#L30). There is a 1-1 mapping from the old names to the new, and it is safe to simply replace the names to update existing models.
  * Specify variable lower/upper bounds in ``@variable`` using the ``lowerbound`` and ``upperbound`` keyword arguments.
  * Change name printed for variable using the ``basename`` keyword argument to ``@variable``.
  * New ``@variables`` macro allows multiline declaration of groups of variables.
  * A number of solver methods previously available only through MathProgBase are now exposed directly in JuMP. The fix was [recorded](https://youtu.be/qF1lZPJ3a5A) live!
  * Compatibility fixes with Julia 0.5.
  * The "end" indexing syntax is no longer supported within JuMPArrays which do not use 1-based indexing until upstream issues are resolved, see [here](https://github.com/JuliaOpt/JuMP.jl/issues/730).

## Version 0.12.2 (March 9, 2016)

  * Small fixes for nonlinear optimization

## Version 0.12.1 (March 1, 2016)

  * Fix a regression in slicing for JuMPArrays (when not using 1-based indexing)

## Version 0.12.0 (February 27, 2016)

  * The automatic differentiation functionality has been completely rewritten with a number of user-facing changes:
      - ``@defExpr`` and ``@defNLExpr`` now take the model as the first argument. The previous one-argument version of ``@defExpr`` is deprecated; all expressions should be named. E.g., replace ``@defExpr(2x+y)`` with ``@defExpr(jump_model, my_expr, 2x+y)``.
      - JuMP no longer uses Julia's variable binding rules for efficiently re-solving a sequence of nonlinear models. Instead, we have introduced nonlinear parameters. This is a breaking change, so we have added a warning message when we detect models that may depend on the old behavior.
      - Support for user-defined functions integrated within nonlinear JuMP expressions.
  * Replaced iteration over ``AffExpr`` with ``Number``-like scalar iteration; previous iteration behavior is now available via ``linearterms(::AffExpr)``.
  * Stopping the solver via ``throw(CallbackAbort())`` from a callback no longer triggers an exception. Instead, ``solve()`` returns ``UserLimit`` status.
  * ``getDual()`` now works for conic problems (Thanks @emreyamangil.)

## Version 0.11.3 (February 4, 2016)

  * Bug-fix for problems with quadratic objectives and semidefinite constraints

## Version 0.11.2 (January 14, 2016)

  * Compatibility update for Mosek

## Version 0.11.1 (December 1, 2015)

  * Remove usage of `@compat` in tests.
  * Fix updating quadratic objectives for nonlinear models.

## Version 0.11.0 (November 30, 2015)

  * Julia 0.4.0 is the minimum required version for this release.
  * Fix for scoping semantics of index variables in sum{}. Index variables no longer leak into the surrounding scope.
  * Addition of the ``solve(m::Model, relaxation=true)`` keyword argument to solve the standard continuous realaxation of model ``m``
  * The ``getConstraintBounds()`` method allows access to the lower and upper bounds of all constraints in a (nonlinear) model.
  * Update for breaking changes in MathProgBase

## Version 0.10.3 (November 20, 2015)

  * Fix a rare error when parsing quadratic expressions
  * Fix ``Variable()`` constructor with default arguments
  * Detect unrecognized keywords in ``solve()``

## Version 0.10.2 (September 28, 2015)

  * Fix for deprecation warnings

## Version 0.10.1 (September 3, 2015)

  * Fixes for ambiguity warnings.
  * Fix for breaking change in precompilation syntax in Julia 0.4-pre

## Version 0.10.0 (August 31, 2015)

  * Support (on Julia 0.4 and later) for conditions in indexing ``@defVar`` and ``@addConstraint`` constructs, e.g. ``@defVar(m, x[i=1:5,j=1:5; i+j >= 3])``
  * Support for vectorized operations on Variables and expressions. See the documentation for details.
  * New ``getVar()`` method to access variables in a model by name
  * Support for semidefinite programming.
  * Dual solutions are now available for general nonlinear problems. You may call ``getDual`` on a reference object for a nonlinear constraint, and ``getDual`` on a variable object for Lagrange multipliers from active bounds.
  * Introduce warnings for two common performance traps: too many calls to ``getValue()`` on a collection of variables and use of the ``+`` operator in a loop to sum expressions.
  * Second-order cone constraints can be written directly with the ``norm()`` and ``norm2{}`` syntax.
  * Implement MathProgBase interface for querying Hessian-vector products.
  * Iteration over ``JuMPContainer``s is deprecated; instead, use the ``keys`` and ``values`` functions, and ``zip(keys(d),values(d))`` for the old behavior.
  * ``@defVar`` returns ``Array{Variable,N}`` when each of ``N`` index sets are of the form ``1:nᵢ``.
  * Module precompilation: on Julia 0.4 and later, ``using JuMP`` is now much faster.

## Version 0.9.3 (August 11, 2015)

  * Fixes for FactCheck testing on julia v0.4.

## Version 0.9.2 (June 27, 2015)

  * Fix bug in @addConstraints.

## Version 0.9.1 (April 25, 2015)

  * Fix for Julia 0.4-dev.
  * Small infrastructure improvements for extensions.

## Version 0.9.0 (April 18, 2015)

  * Comparison operators for constructing constraints (e.g. ``2x >= 1``) have been deprecated. Instead, construct the constraints explicitly in
    the ``@addConstraint`` macro to add them to the model, or in the ``@LinearConstraint`` macro to create a stand-alone linear constraint instance.
  * ``getValue()`` method implemented to compute the value of a nonlinear subexpression
  * JuMP is now released under the Mozilla Public License version 2.0 (was previously LGPL). MPL is a copyleft license which is less restrictive than LGPL, especially for embedding JuMP within other applications.
  * A number of performance improvements in ReverseDiffSparse for computing derivatives.
  * ``MathProgBase.getsolvetime(m)`` now returns the solution time reported by the solver, if available. (Thanks @odow, Oscar Dowson)
  * Formatting fix for LP format output. (Thanks @sbebo, Leonardo Taccari).

## Version 0.8.0 (February 17, 2015)

  * Nonlinear subexpressions now supported with the ``@defNLExpr`` macro.
  * SCS supported for solving second-order conic problems.
  * ``setXXXCallback`` family deprecated in favor of ``addXXXCallback``.
  * Multiple callbacks of the same type can be registered.
  * Added support for informational callbacks via ``addInfoCallback``.
  * A ``CallbackAbort`` exception can be thrown from callback to safely exit optimization.

## Version 0.7.4 (February 4, 2015)

  * Reduced costs and linear constraint duals are now accessible when quadratic constraints are present.
  * Two-sided nonlinear constraints are supported.
  * Methods for accessing the number of variables and constraints in a model are renamed.
  * New default procedure for setting initial values in nonlinear optimization: project zero onto the variable bounds.
  * Small bug fixes.


## Version 0.7.3 (January 14, 2015)

  * Fix a method ambiguity conflict with Compose.jl (cosmetic fix)

## Version 0.7.2 (January 9, 2015)

  * Fix a bug in ``sum(::JuMPDict)``
  * Added the ``setCategory`` function to change a variables category (e.g. continuous or binary)
  after construction, and ``getCategory`` to retrieve the variable category.

## Version 0.7.1 (January 2, 2015)

  * Fix a bug in parsing linear expressions in macros. Affects only Julia 0.4 and later.

## Version 0.7.0 (December 29, 2014)

### Linear/quadratic/conic programming

  * **Breaking change**: The syntax for column-wise model generation has been changed to use keyword arguments in ``@defVar``.
  * On Julia 0.4 and later, variables and coefficients may be multiplied in any order within macros. That is, variable*coefficient is now valid syntax.
  * ECOS supported for solving second-order conic problems.

### Nonlinear programming

  * Support for skipping model generation when solving a sequence of nonlinear models with changing data.
  * Fix a memory leak when solving a sequence of nonlinear models.
  * The ``@addNLConstraint`` macro now supports the three-argument version to define sets of nonlinear constraints.
  * KNITRO supported as a nonlinear solver.
  * Speed improvements for model generation.
  * The ``@addNLConstraints`` macro supports adding multiple (groups of) constraints at once. Syntax is similar to ``@addConstraints``.
  * Discrete variables allowed in nonlinear problems for solvers which support them (currently only KNITRO).

### General

  * Starting values for variables may now be specified with ``@defVar(m, x, start=value)``.
  * The ``setSolver`` function allows users to change the solver subsequent to model creation.
  * Support for "fixed" variables via the ``@defVar(m, x == 1)`` syntax.
  * Unit tests rewritten to use FactCheck.jl, improved testing across solvers.

## Version 0.6.3 (October 19, 2014)

  * Fix a bug in multiplying two AffExpr objects.


## Version 0.6.2 (October 11, 2014)

  * Further improvements and bug fixes for printing.
  * Fixed a bug in ``@defExpr``.
  * Support for accessing expression graphs through the MathProgBase NLP interface.

## Version 0.6.1 (September 19, 2014)

  * Improvements and bug fixes for printing.

## Version 0.6.0 (September 9, 2014)

  * Julia 0.3.0 is the minimum required version for this release.
  * ``buildInternalModel(m::Model)`` added to build solver-level model in memory without optimizing.
  * Deprecate ``load_model_only`` keyword argument to ``solve``.
  * Add groups of constraints with ``@addConstraints`` macro.
  * Unicode operators now supported, including ``∑`` for ``sum``, ``∏`` for ``prod``, and ``≤``/``≥``
  * Quadratic constraints supported in ``@addConstraint`` macro.
  * Quadratic objectives supported in ``@setObjective`` macro.
  * MathProgBase solver-independent interface replaces Ipopt-specific interface for nonlinear problems
    - **Breaking change**: ``IpoptOptions`` no longer supported to specify solver options, use ``m = Model(solver=IpoptSolver(options...))`` instead.
  * New solver interfaces: ECOS, NLopt, and nonlinear support for MOSEK
  * New option to control whether the lazy constraint callback is executed at each node in the B&B tree or just when feasible solutions are found
  * Add support for semicontinuous and semi-integer variables for those solvers that support them.
  * Add support for index dependencies (e.g. triangular indexing) in ``@defVar``, ``@addConstraint``, and ``@defExpr`` (e.g. ``@defVar(m, x[i=1:10,j=i:10])``).
    - This required some changes to the internal structure of JuMP containers, which may break code that explicitly stored ``JuMPDict`` objects.

## Version 0.5.8 (September 24, 2014)

  * Fix a bug with specifying solvers (affects Julia 0.2 only)

## Version 0.5.7 (September 5, 2014)

  * Fix a bug in printing models

## Version 0.5.6 (September 2, 2014)
  * Add support for semicontinuous and semi-integer variables for those solvers that support them.
    - **Breaking change**: Syntax for ``Variable()`` constructor has changed (use of this interface remains discouraged)
  * Update for breaking changes in MathProgBase

## Version 0.5.5 (July 6, 2014)

  * Fix bug with problem modification: adding variables that did not appear in existing constraints or objective.

## Version 0.5.4 (June 19, 2014)

  * Update for breaking change in MathProgBase which reduces loading times for ``using JuMP``
  * Fix error when MIPs not solved to optimality


## Version 0.5.3 (May 21, 2014)

  * Update for breaking change in ReverseDiffSparse

## Version 0.5.2 (May 9, 2014)

  * Fix compatibility with Julia 0.3 prerelease


## Version 0.5.1 (May 5, 2014)

  * Fix a bug in coefficient handling inside lazy constraints and user cuts

## Version 0.5.0 (May 2, 2014)

  * Support for nonlinear optimization with exact, sparse second-order derivatives automatically computed. Ipopt is currently the only solver supported.
  * ``getValue`` for ``AffExpr`` and ``QuadExpr``
  * **Breaking change**: ``getSolverModel`` replaced by ``getInternalModel``, which returns the internal MathProgBase-level model
  * Groups of constraints can be specified with ``@addConstraint`` (see documentation for details). This is not a breaking change.
  * ``dot(::JuMPDict{Variable},::JuMPDict{Variable})`` now returns the corresponding quadratic expression.

## Version 0.4.1 (March 24, 2014)

  * Fix bug where change in objective sense was ignored when re-solving a model.
  * Fix issue with handling zero coefficients in AffExpr.

## Version 0.4.0 (March 10, 2014)

  * Support for SOS1 and SOS2 constraints.
  * Solver-independent callback for user heuristics.
  * ``dot`` and ``sum`` implemented for ``JuMPDict`` objects. Now you can say ``@addConstraint(m, dot(a,x) <= b)``.
  * Developers: support for extensions to JuMP. See definition of Model in ``src/JuMP.jl`` for more details.
  * Option to construct the low-level model before optimizing.

## Version 0.3.2 (February 17, 2014)

 * Improved model printing
   - Preliminary support for IJulia output

## Version 0.3.1 (January 30, 2014)

 * Documentation updates
   - Support for MOSEK
   - CPLEXLink renamed to CPLEX

## Version 0.3.0 (January 21, 2014)

 * Unbounded/infeasibility rays: getValue() will return the corresponding
   components of an unbounded ray when a model is unbounded, if supported
   by the selected solver. getDual() will return an infeasibility ray (Farkas proof)
   if a model is infeasible and the selected solver supports this feature.
 * Solver-independent callbacks for user generated cuts.
 * Use new interface for solver-independent QCQP.
 * ``setlazycallback`` renamed to ``setLazyCallback`` for consistency.

## Version 0.2.0 (December 15, 2013)

  * **Breaking change**: Objective sense is specified in setObjective
    instead of in the Model constructor.
  * **Breaking change**: ``lpsolver`` and ``mipsolver`` merged into
    single ``solver`` option.
  * Problem modification with efficient LP restarts and MIP warm-starts.
  * Relatedly, column-wise modeling now supported.
  * Solver-independent callbacks supported. Currently we support only
    a "lazy constraint" callback, which works with Gurobi, CPLEX, and GLPK.
    More callbacks coming soon.

## Version 0.1.2 (November 16, 2013)

  * Bug fixes for printing, improved error messages.
  * Allow ``AffExpr`` to be used in macros; e.g.,
    ``ex = y + z; @addConstraint(m, x + 2*ex <= 3)``

## Version 0.1.1 (October 23, 2013)

  * Update for solver specification API changes in MathProgBase.

## Version 0.1.0 (October 3, 2013)

  * Initial public release.
# How to use a custom binary

Many solvers are not written in Julia, but instead in languages like C or C++.
JuMP interacts with these solvers through binary dependencies.

For many open-source solvers, we automatically install the appropriate binary
when you run `Pkg.add("Solver")`. For example, `Pkg.add("ECOS")` will also
install the ECOS binary.

This page explains how this installation works, and how you can use a custom
binary.

!!! compat
    These instructions require Julia 1.6 or later.

## Background

Each solver that JuMP supports is structured as a Julia package. For example,
the interface for the [ECOS](https://github.com/embotech/ecos) solver is
provided by the [ECOS.jl](https://github.com/jump-dev/ECOS.jl) package.

!!! tip
    This page uses the example of ECOS.jl because it is simple to compile. Other
    solvers follow similar conventions. For example, the interface to the Clp
    solver is provided by Clp.jl.

The ECOS.jl package provides an interface between the C API of ECOS and
MathOptInterface. However, it does not handle the installation of the solver
binary; that is the job for a JLL package.

A JLL is a Julia package that wraps a pre-compiled binary.
Binaries are built using [Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil)
(for example, [ECOS](https://github.com/JuliaPackaging/Yggdrasil/blob/master/E/ECOS/build_tarballs.jl))
and hosted in the [JuliaBinaryWrappers](https://github.com/JuliaBinaryWrappers)
GitHub repository (for example, [ECOS_jll.jl](https://github.com/JuliaBinaryWrappers/ECOS_jll.jl)).

JLL packages contain little code. Their only job is to `dlopen` a dynamic
library, along with any dependencies.

JLL packages manage their binary dependencies using [Julia's artifact system](https://pkgdocs.julialang.org/v1/artifacts/).
Each JLL package has an `Artifacts.toml` file which describes where to find each
binary artifact for each different platform that it might be installed on. Here
is the [Artifacts.toml file for ECOS_jll.jl](https://github.com/JuliaBinaryWrappers/ECOS_jll.jl/blob/main/Artifacts.toml).

The binaries installed by the JLL package should be sufficient for most users.
In rare cases, however, you may require a custom binary. The two main reasons to
use a custom binary are:

 * You want a binary with custom compilation settings (for example, debugging)
 * You want a binary with a set of dependencies that are not available on
   Yggdrasil (for example, a commerial solver like Gurobi or CPLEX).

The following sections explain how to replace the binaries provided by a JLL
package with the custom ones you have compiled. As a reminder, we use ECOS as an
example for simplicity, but the steps are the same for other solvers.

## [Explore the JLL you want to override](@id jll_structure)

The first step is to explore the structure and filenames of the JLL package we
want to override.

Find the location of the files using `.artifact_dir`:
```julia
julia> using ECOS_jll

julia> ECOS_jll.artifact_dir
"/Users/oscar/.julia/artifacts/2addb75332eff5a1657b46bb6bf30d2410bc7ecf"
```

!!! tip
    This path may be different on other machines.

Here is what it contains:
```julia
julia> readdir(ECOS_jll.artifact_dir)
4-element Vector{String}:
 "include"
 "lib"
 "logs"
 "share"

julia> readdir(joinpath(ECOS_jll.artifact_dir, "lib"))
1-element Vector{String}:
 "libecos.dylib"
```

Other solvers may have a `bin` directory containing executables. To use a custom
binary of ECOS, we need to replace `/lib/libecos.dylib` with our custom binary.

## [Compile a custom binary](@id compile_ecos)

The next step is to compile a custom binary. Because ECOS is written in C with
no dependencies, this is easy to do if you have a C compiler:
```julia
oscar@Oscars-MBP jll_example % git clone https://github.com/embotech/ecos.git
[... lines omitted ...]
oscar@Oscars-MBP jll_example % cd ecos
oscar@Oscars-MBP ecos % make shared
[... many lines omitted...]
oscar@Oscars-MBP ecos % mkdir lib
oscar@Oscars-MBP ecos % cp libecos.dylib lib
```

!!! warning
    Compiling custom solver binaries is an advanced operation. Due to the
    complexities of compiling various solvers, the JuMP community is unable to
    help you diagnose and fix compilation issues.

After this compilation step, we now have a folder `/tmp/jll_example/ecos`
that contains `lib` and `include` directories with the same files as `ECOS_jll`:
```julia
julia> readdir(joinpath("ecos", "lib"))
1-element Vector{String}:
 "libecos.dylib"
```

## Overriding a single library

To override the `libecos` library, we need to know what `ECOS_jll` calls it. (In
most cases, it will also be `libecos`, but not always.)

There are two ways you can check.

 1. Check the bottom of the JLL's GitHub README. For example,
    [ECOS_jll](https://github.com/JuliaBinaryWrappers/ECOS_jll.jl#products) has
    a single `LibraryProduct` called `libecos`.
 2. Type `ECOS_jll.` and the press the `[TAB]` key twice to auto-complete
    available options:
    ```julia
    julia> ECOS_jll.
    LIBPATH           PATH_list          best_wrapper       get_libecos_path   libecos_handle
    LIBPATH_list      __init__           dev_jll            is_available       libecos_path
    PATH              artifact_dir       find_artifact_dir  libecos
    ```
    Here you can see there is `libecos`, and more usefully for us,
    `libecos_path`.

Once you know the name of the variable to override (the one that ends in
`_path`), use  [Preferences.jl](https://github.com/JuliaPackaging/Preferences.jl)
to specify a new path:
```julia
using Preferences
set_preferences!(
    "LocalPreferences.toml",
    "ECOS_jll",
    "libecos_path" => "/tmp/jll_example/ecos/lib/libecos"
)
```

This will create a file in your current directory called `LocalPreferences.toml`
with the contents:
```julia
[ECOS_jll]
libecos_path = "/tmp/jll_example/ecos/lib/libecos"
```

Now if you restart Julia, you will see:
```julia
julia> using ECOS_jll

julia> ECOS_jll.libecos
"/tmp/jll_example/ecos/lib/libecos"
```

To go back to using the default library, just delete the `LocalPreferences.toml`
file.

## Overriding an entire artifact

Sometimes a solver may provide a number of libraries and executables, and
specifying the path for each of the becomes tedious. In this case, we can use
Julia's `Override.toml` to replace an entire artifact.

Overriding an entire artifact requires you to replicate the structure and
contents of the JLL package that we [explored above](@ref jll_structure).

In most cases you need only reproduce the `include`, `lib`, and `bin`
directories (if they exist). You can safely ignore any `logs` or `share`
directories. Take careful note of what files each directory contains and what
they are called.

For our ECOS example, we already reproduced the structure when we
[compiled ECOS](@ref compile_ecos).

So, now we need to tell Julia to use our custom installation instead of
the default. We can do this by making an override file at
`~/.julia/artifacts/Overrides.toml`.

`Overrides.toml` has the following content:
```julia
# Override for ECOS_jll
2addb75332eff5a1657b46bb6bf30d2410bc7ecf = "/tmp/jll_example/ecos"
```
where `2addb75332eff5a1657b46bb6bf30d2410bc7ecf` is the folder from the original
`ECOS_jll.artifact_dir` and `"/tmp/jll_example/ecos"` is the location of our new
installation. Replace these as appropriate for your system.

If you restart Julia after creating the override file, you will see:
```julia
julia> using ECOS_jll

julia> ECOS_jll.artifact_dir
"/tmp/jll_example/ecos"
```
Now when we use ECOS it will use our custom binary.

## Using Cbc with a custom binary

As a second example, we demonstrate how to use
[Cbc.jl](https://github.com/jump-dev/Cbc.jl) with a custom binary.

### Explore the JLL you want to override

First, let's check where `Cbc_jll` is installed:
```julia
julia> using Cbc_jll

julia> Cbc_jll.artifact_dir
"/Users/oscar/.julia/artifacts/e481bc81db5e229ba1f52b2b4bd57484204b1b06"

julia> readdir(Cbc_jll.artifact_dir)
5-element Vector{String}:
 "bin"
 "include"
 "lib"
 "logs"
 "share"

julia> readdir(joinpath(Cbc_jll.artifact_dir, "bin"))
1-element Vector{String}:
 "cbc"

julia> readdir(joinpath(Cbc_jll.artifact_dir, "lib"))
10-element Vector{String}:
 "libCbc.3.10.5.dylib"
 "libCbc.3.dylib"
 "libCbc.dylib"
 "libCbcSolver.3.10.5.dylib"
 "libCbcSolver.3.dylib"
 "libCbcSolver.dylib"
 "libOsiCbc.3.10.5.dylib"
 "libOsiCbc.3.dylib"
 "libOsiCbc.dylib"
 "pkgconfig"
```

### Compile a custom binary

Next, we need to compile Cbc. Cbc can be difficult to compile (it has a lot of
dependencies), but for macOS users there is a homebrew recipe:
```
(base) oscar@Oscars-MBP jll_example % brew install cbc
[ ... lines omitted ... ]
(base) oscar@Oscars-MBP jll_example % brew list cbc
/usr/local/Cellar/cbc/2.10.5/bin/cbc
/usr/local/Cellar/cbc/2.10.5/include/cbc/ (76 files)
/usr/local/Cellar/cbc/2.10.5/lib/libCbc.3.10.5.dylib
/usr/local/Cellar/cbc/2.10.5/lib/libCbcSolver.3.10.5.dylib
/usr/local/Cellar/cbc/2.10.5/lib/libOsiCbc.3.10.5.dylib
/usr/local/Cellar/cbc/2.10.5/lib/pkgconfig/ (2 files)
/usr/local/Cellar/cbc/2.10.5/lib/ (6 other files)
/usr/local/Cellar/cbc/2.10.5/share/cbc/ (59 files)
/usr/local/Cellar/cbc/2.10.5/share/coin/ (4 files)
```

### Override single libraries

To use `Preferences.jl` to override specific libraries we first check the names
of each library in `Cbc_jll`:
```julia
julia> Cbc_jll.
LIBPATH               cbc                    get_libcbcsolver_path  libOsiCbc_path
LIBPATH_list          cbc_path               is_available           libcbcsolver
PATH                  dev_jll                libCbc                 libcbcsolver_handle
PATH_list             find_artifact_dir      libCbc_handle          libcbcsolver_path
__init__              get_cbc_path           libCbc_path
artifact_dir          get_libCbc_path        libOsiCbc
best_wrapper          get_libOsiCbc_path     libOsiCbc_handle
```

Then we add the following to `LocalPreferences.toml`:
```julia
[Cbc_jll]
cbc_path = "/usr/local/Cellar/cbc/2.10.5/bin/cbc"
libCbc_path = "/usr/local/Cellar/cbc/2.10.5/lib/libCbc.3.10.5"
libOsiCbc_path = "/usr/local/Cellar/cbc/2.10.5/lib/libOsiCbc.3.10.5"
libcbcsolver_path = "/usr/local/Cellar/cbc/2.10.5/lib/libCbcSolver.3.10.5"
```

!!! info
    Note that capitalization matters, so `libcbcsolver_path` corresponds to
    `libCbcSolver.3.10.5`.

### Override entire artifact

To use the homebrew install as our custom binary we add the following to
`~/.julia/artifacts/Overrides.toml`:
```julia
# Override for Cbc_jll
e481bc81db5e229ba1f52b2b4bd57484204b1b06 = "/usr/local/Cellar/cbc/2.10.5"
```
```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# [Extensions](@id extensions_manual)

JuMP provides a variety of ways to extend the basic modeling functionality.

!!! tip
    This documentation in this section is still a work-in-progress. The best
    place to look for ideas and help when writing a new JuMP extension are
    existing JuMP extensions. Examples include:
     * [BilevelJuMP.jl](https://github.com/joaquimg/BilevelJuMP.jl)
     * [Coluna.jl](https://github.com/atoptima/Coluna.jl)
     * [InfiniteOpt.jl](https://github.com/pulsipher/InfiniteOpt.jl)
     * [Plasmo.jl](https://github.com/zavalab/Plasmo.jl)
     * [PolyJuMP.jl](https://github.com/jump-dev/PolyJuMP.jl)
     * [SDDP.jl](https://github.com/odow/SDDP.jl)
     * [StochasticPrograms.jl](https://github.com/martinbiel/StochasticPrograms.jl)
     * [SumOfSquares.jl](https://github.com/jump-dev/SumOfSquares.jl)
     * [vOptGeneric.jl](https://github.com/vOptSolver/vOptGeneric.jl)

## Define a new set

To define a new set for JuMP, subtype `MOI.AbstractScalarSet` or
`MOI.AbstractVectorSet` and implement `Base.copy` for the set. That's it!

```jldoctest define_new_set
struct _NewVectorSet <: MOI.AbstractVectorSet
    dimension::Int
end
Base.copy(x::_NewVectorSet) = x

model = Model()
@variable(model, x[1:2])
@constraint(model, x in _NewVectorSet(2))

# output

[x[1], x[2]] ∈ _NewVectorSet(2)
```

However, for vector-sets, this requires the user to specify the dimension
argument to their set, even though we could infer it from the length of `x`!

You can make a more user-friendly set by subtyping [`AbstractVectorSet`](@ref)
and implementing [`moi_set`](@ref).

```jldoctest define_new_set
struct NewVectorSet <: JuMP.AbstractVectorSet end
JuMP.moi_set(::NewVectorSet, dim::Int) = _NewVectorSet(dim)

model = Model()
@variable(model, x[1:2])
@constraint(model, x in NewVectorSet())

# output

[x[1], x[2]] ∈ _NewVectorSet(2)
```

## Extend [`@variable`](@ref)

Just as `Bin` and `Int` create binary and integer variables, you can extend
the [`@variable`](@ref) macro to create new types of variables. Here is an
explanation by example, where we create a `AddTwice` type, that creates a tuple
of two JuMP variables instead of a single variable.

First, create a new struct. This can be anything. Our struct holds a
[`VariableInfo`](@ref) object that stores bound information, and whether the
variable is binary or integer.
```jldoctest new_variable
julia> struct AddTwice
           info::JuMP.VariableInfo
       end
```

Second, implement [`build_variable`](@ref), which takes `::Type{AddTwice}` as
an argument, and returns an instance of `AddTwice`. Note that you can also
receive keyword arguments.
```jldoctest new_variable
julia> function JuMP.build_variable(
           _err::Function,
           info::JuMP.VariableInfo,
           ::Type{AddTwice};
           kwargs...
       )
           println("Can also use $kwargs here.")
           return AddTwice(info)
       end
```

Third, implement [`add_variable`](@ref), which takes the instance of `AddTwice`
from the previous step, and returns something. Typically, you will want to call
[`add_variable`](@ref) here. For example, our `AddTwice` call is going to add
two JuMP variables.
```jldoctest new_variable
julia> function JuMP.add_variable(
           model::JuMP.Model,
           duplicate::AddTwice,
           name::String,
       )
           a = JuMP.add_variable(
               model,
               JuMP.ScalarVariable(duplicate.info),
               name * "_a",
            )
           b = JuMP.add_variable(
               model,
               JuMP.ScalarVariable(duplicate.info),
               name * "_b",
            )
           return (a, b)
       end
```

Now `AddTwice` can be passed to [`@variable`](@ref) similar to `Bin` or `Int`.
However, now it adds two variables instead of one!
```jldoctest new_variable
julia> model = Model();

julia> @variable(model, x[i=1:2], AddTwice, kw=i)
Can also use Base.Iterators.Pairs(:kw => 1) here.
Can also use Base.Iterators.Pairs(:kw => 2) here.
2-element Vector{Tuple{VariableRef, VariableRef}}:
 (x[1]_a, x[1]_b)
 (x[2]_a, x[2]_b)

julia> num_variables(model)
4

julia> first(x[1])
x[1]_a

julia> last(x[2])
x[2]_b
```

## Extend [`@constraint`](@ref)

The [`@constraint`](@ref) macro has three steps that can be intercepted and
extended: parse time, build time, and add time.

### Parse

To extend the [`@constraint`](@ref) macro at parse time, implement one of the
following methods:

 * [`parse_constraint_head`](@ref)
 * [`parse_constraint_call`](@ref)

!!! warning
    Extending the constraint macro at parse time is an advanced operation and
    has the potential to interfere with existing JuMP syntax. Please discuss
    with the [developer chatroom](https://gitter.im/JuliaOpt/jump-dev) before
    publishing any code that implements these methods.

[`parse_constraint_head`](@ref) should be implemented to intercept an expression
based on the `.head` field of `Base.Expr`. For example:
```jldoctest
julia> using JuMP

julia> const MutableArithmetics = JuMP._MA;

julia> model = Model(); @variable(model, x);

julia> function JuMP.parse_constraint_head(
           _error::Function,
           ::Val{:(:=)},
           lhs,
           rhs,
       )
           println("Rewriting := as ==")
           new_lhs, parse_code = MutableArithmetics.rewrite(lhs)
           build_code = :(
               build_constraint($(_error), $(new_lhs), MOI.EqualTo($(rhs)))
           )
           return false, parse_code, build_code
       end

julia> @constraint(model, x + x := 1.0)
Rewriting := as ==
2 x = 1.0
```

[`parse_constraint_call`](@ref) should be implemented to intercept an expression
of the form `Expr(:call, op, args...)`. For example:
```jldoctest
julia> using JuMP

julia> const MutableArithmetics = JuMP._MA;

julia> model = Model(); @variable(model, x);

julia> function JuMP.parse_constraint_call(
           _error::Function,
           is_vectorized::Bool,
           ::Val{:my_equal_to},
           lhs,
           rhs,
       )
           println("Rewriting my_equal_to to ==")
           new_lhs, parse_code = MutableArithmetics.rewrite(lhs)
           build_code = if is_vectorized
               :(build_constraint($(_error), $(new_lhs), MOI.EqualTo($(rhs)))
           )
           else
               :(build_constraint.($(_error), $(new_lhs), MOI.EqualTo($(rhs))))
           end
           return parse_code, build_code
       end

julia> @constraint(model, my_equal_to(x + x, 1.0))
Rewriting my_equal_to to ==
2 x = 1.0
```

!!! tip
    When parsing a constraint you can recurse into sub-constraint (e.g., the
    `{expr}` in `z => {x <= 1}`) by calling [`parse_constraint`](@ref).

### Build

To extend the [`@constraint`](@ref) macro at build time, implement a new
[`build_constraint`](@ref) method.

This may mean implementing a method for a specific function or set created at
parse time, or it may mean implementing a method which handles additional
positional arguments.

[`build_constraint`](@ref) must return an [`AbstractConstraint`](@ref), which
can either be an [`AbstractConstraint`](@ref) already supported by JuMP, e.g., `ScalarConstraint` or `VectorConstraint`, or a custom
[`AbstractConstraint`](@ref) with a corresponding [`add_constraint`](@ref)
method (see [Add](@ref extension_add_constraint)).

!!! tip
    The easiest way to extend [`@constraint`](@ref) is via an additional
    positional argument to [`build_constraint`](@ref).

Here is an example of adding extra arguments to [`build_constraint`](@ref):
```jldoctest
julia> model = Model(); @variable(model, x);

julia> struct MyConstrType end

julia> function JuMP.build_constraint(
            _error::Function,
            f::JuMP.GenericAffExpr,
            set::MOI.EqualTo,
            extra::Type{MyConstrType};
            d = 0,
       )
            new_set = MOI.LessThan(set.value + d)
            return JuMP.build_constraint(_error, f, new_set)
       end

julia> @constraint(model, my_con, x == 0, MyConstrType, d = 2)
my_con : x ≤ 2.0
```

!!! note
    Only a single positional argument can be given to a particular constraint.
    Extensions that seek to pass multiple arguments (e.g., `Foo` and `Bar`)
    should combine them into one argument type (e.g., `FooBar`).

### [Add](@id extension_add_constraint)

[`build_constraint`](@ref) returns an [`AbstractConstraint`](@ref) object. To
extend [`@constraint`](@ref) at add time, define a subtype of
[`AbstractConstraint`](@ref), implement [`build_constraint`](@ref) to return an
instance of the new type, and then implement [`add_constraint`](@ref).

Here is an example:
```jldoctest
julia> model = Model(); @variable(model, x);

julia> struct MyTag
           name::String
       end

julia> struct MyConstraint{S} <: AbstractConstraint
           name::String
           f::AffExpr
           s::S
       end

julia> function JuMP.build_constraint(
            _error::Function,
            f::AffExpr,
            set::MOI.AbstractScalarSet,
            extra::MyTag,
       )
            return MyConstraint(extra.name, f, set)
       end

julia> function JuMP.add_constraint(
            model::Model,
            con::MyConstraint,
            name::String,
       )
            return add_constraint(
                model,
                ScalarConstraint(con.f, con.s),
                "$(con.name)[$(name)]",
            )
       end

julia> @constraint(model, my_con, 2x <= 1, MyTag("my_prefix"))
my_prefix[my_con] : 2 x - 1 ≤ 0.0
```

## The extension dictionary

Every JuMP model has a field `.ext::Dict{Symbol,Any}` that can be used by
extensions. This is useful if your extensions to [`@variable`](@ref) and
[`@constraint`](@ref) need to store information between calls.

The most common way to initialize a model with information in the `.ext`
dictionary is to provide a new constructor:
```jldoctest
julia> function MyModel()
           model = Model()
           model.ext[:MyModel] = 1
           return model
       end
MyModel (generic function with 1 method)

julia> model = MyModel()
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> model.ext
Dict{Symbol, Any} with 1 entry:
  :MyModel => 1
```

If you define extension data, implement [`copy_extension_data`](@ref)
to support [`copy_model`](@ref).

## Defining new JuMP models

If extending individual calls to [`@variable`](@ref) and [`@constraint`](@ref)
is not sufficient, it is possible to implement a new model via a subtype of
[`AbstractModel`](@ref). You can also define new [`AbstractVariableRef`](@ref)s
to create different types of JuMP variables.

!!! warning
    Extending JuMP in this manner is an advanced operation. We strongly
    encourage you to consider how you can use the methods mentioned in the
    previous sections to achieve your aims instead of defining new model and
    variable types. Consult the [developer chatroom](https://gitter.im/JuliaOpt/jump-dev)
    _before_ starting work on this.

If you define new types, you will need to implement a considerable number of
methods, and doing so will require a detailed understanding of the JuMP
internals. Therefore, the list of methods to implement is currently
undocumented.

The easiest way to extend JuMP by defining a new model type is to follow an
existing example. A simple example to follow is the [JuMPExtension module](https://github.com/jump-dev/JuMP.jl/blob/master/test/JuMPExtension.jl)
in the JuMP test suite. The best example of an external JuMP extension that
implements an [`AbstractModel`](@ref) is [InfiniteOpt.jl](https://github.com/pulsipher/InfiniteOpt.jl).

# How to contribute to JuMP

Welcome! This document explains some ways you can contribute to JuMP.

## Code of Conduct

This project and everyone participating in it is governed by the
[JuMP Code of Conduct](https://github.com/jump-dev/JuMP.jl/blob/master/CODE_OF_CONDUCT.md).
By participating, you are expected to uphold this code.

## Join the community forum

First up, join the [community forum](https://discourse.julialang.org/c/domain/opt).

The forum is a good place to ask questions about how to use JuMP. You can also
use the forum to discuss possible feature requests and bugs before raising a
GitHub issue (more on this below).

Aside from asking questions, the easiest way you can contribute to JuMP is to
help answer questions on the forum!

## Join the developer chatroom

If you're interested in contributing code to JuMP, the next place to join is the
[developer chatroom](https://gitter.im/JuliaOpt/JuMP-dev). Let us know what you
have in mind, and we can point you in the right direction.

## Improve the documentation

Chances are, if you asked (or answered) a question on the community forum, then
it is a sign that the [documentation](https://jump.dev/JuMP.jl/dev/) could be
improved. Moreover, since it is your question, you are probably the best-placed
person to improve it!

The docs are written in Markdown and are built using
[Documenter.jl](https://github.com/JuliaDocs/Documenter.jl).
You can find the source of all the docs
[here](https://github.com/jump-dev/JuMP.jl/tree/master/docs).

If your change is small (like fixing typos, or one or two sentence corrections),
the easiest way to do this is via GitHub's online editor. (GitHub has
[help](https://help.github.com/articles/editing-files-in-another-user-s-repository/)
on how to do this.)

If your change is larger, or touches multiple files, you will need to make the
change locally and then use Git to submit a pull request. (See
[Contribute code to JuMP](@ref) below for more on this.)

!!! tip
    If you need any help, come join the
    [developer chatroom](https://gitter.im/JuliaOpt/JuMP-dev) and we will walk
    you through the process.

## File a bug report

Another way to contribute to JuMP is to file
[bug reports](https://github.com/jump-dev/JuMP.jl/issues/new?template=bug_report.md).

Make sure you read the info in the box where you write the body of the issue
before posting. You can also find a copy of that info
[here](https://github.com/jump-dev/JuMP.jl/blob/master/.github/ISSUE_TEMPLATE/bug_report.md).

!!! tip
    If you're unsure whether you have a real bug, post on the
    [community forum](https://discourse.julialang.org/c/domain/opt)
    first. Someone will either help you fix the problem, or let you know the
    most appropriate place to open a bug report.

## Contribute code to JuMP

Finally, you can also contribute code to JuMP!

!!! warning
    If you do not have experience with Git, GitHub, and Julia development, the
    first steps can be a little daunting. However, there are lots of tutorials
    available online, including these for:
     * [GitHub](https://guides.github.com/activities/hello-world/)
     * [Git and GitHub](https://try.github.io/)
     * [Git](https://git-scm.com/book/en/v2)
     * [Julia package development](https://docs.julialang.org/en/v1/stdlib/Pkg/#Developing-packages-1)
    If you need any help, come join the [developer chatroom](https://gitter.im/JuliaOpt/JuMP-dev)
    and we will walk you through the process.

Once you are familiar with Git and GitHub, the workflow for contributing code to
JuMP is similar to the following:

**Step 1: decide what to work on**

The first step is to find an [open issue](https://github.com/jump-dev/JuMP.jl/issues)
(or open a new one) for the problem you want to solve. Then, _before_ spending
too much time on it, discuss what you are planning to do in the issue to see if
other contributors are fine with your proposed changes. Getting feedback early can
improve code quality, and avoid time spent writing code that does not get merged into
JuMP.

!!! tip
    At this point, remember to be patient and polite; you may get a _lot_ of
    comments on your issue! However, do not be afraid! Comments mean that people are
    willing to help you improve the code that you are contributing to JuMP.

**Step 2: fork JuMP**

Go to [https://github.com/jump-dev/JuMP.jl](https://github.com/jump-dev/JuMP.jl)
and click the "Fork" button in the top-right corner. This will create a copy of
JuMP under your GitHub account.

**Step 3: install JuMP locally**

Open Julia and run:
```julia
] dev JuMP
```
This will download the JuMP Git repository to `~/.julia/dev/JuMP`. If you're on
Windows, this will be `C:\\Users\\<my_name>\\.julia\\dev\\JuMP`.

!!! warning
    `] command` means "first type `]` to enter the Julia pkg mode, then type the
    rest. Don't copy-paste the code directly.

**Step 4: checkout a new branch**

!!! note
    In the following, replace any instance of `GITHUB_ACCOUNT` with your GitHub
    user name.

The next step is to checkout a development branch. In a terminal (or command
prompt on Windows), run:
```
$ cd ~/.julia/dev/JuMP

$ git remote add GITHUB_ACCOUNT https://github.com/GITHUB_ACCOUNT/JuMP.jl.git

$ git checkout master

$ git pull

$ git checkout -b my_new_branch
```

!!! tip
    Lines starting with `$` mean "run these in a terminal (command prompt on
    Windows)."

**Step 5: make changes**

Now make any changes to the source code inside the `~/.julia/dev/JuMP`
directory.

Make sure you:
 * Follow the [Style guide](@ref) and run [JuliaFormatter](@ref)
 * Add tests and documentation for any changes or new features

!!! tip
    When you change the source code, you'll need to restart Julia for the
    changes to take effect. This is a pain, so install
    [Revise.jl](https://github.com/timholy/Revise.jl).

**Step 6a: test your code changes**

To test that your changes work, run the JuMP test-suite by opening Julia and
running:
```julia
cd("~/.julia/dev/JuMP")
] activate .
] test
```

!!! warning
    Running the tests might take a long time (~10--15 minutes).

!!! tip
    If you're using Revise.jl, you can also run the tests by calling `include`:
    ```julia
    include("test/runtests.jl")
    ```
    This can be faster if you want to re-run the tests multiple times.

**Step 6b: test your documentation changes**

Open Julia, then run:
```julia
cd("~/.julia/dev/JuMP/docs")
] activate .
include("src/make.jl")
```

!!! warning
    Building the documentation might take a long time (~10 minutes).

!!! tip
    If there's a problem with the tests that you don't know how to fix, don't
    worry. Continue to step 5, and one of the JuMP contributors will comment
    on your pull request telling you how to fix things.

**Step 7: make a pull request**

Once you've made changes, you're ready to push the changes to GitHub. Run:
```
$ cd ~/.julia/dev/JuMP

$ git add .

$ git commit -m "A descriptive message of the changes"

$ git push -u GITHUB_ACCOUNT my_new_branch
```

Then go to [https://github.com/jump-dev/JuMP.jl](https://github.com/jump-dev/JuMP.jl)
and follow the instructions that pop up to open a pull request.

**Step 8: respond to comments**

At this point, remember to be patient and polite; you may get a _lot_ of
comments on your pull request! However, do not be afraid! A lot of comments
means that people are willing to help you improve the code that you are
contributing to JuMP.

To respond to the comments, go back to step 5, make any changes, test the
changes in step 6, and then make a new commit in step 7. Your PR will
automatically update.

**Step 9: cleaning up**

Once the PR is merged, clean-up your Git repository ready for the
next contribution!
```
$ cd ~/.julia/dev/JuMP

$ git checkout master

$ git pull
```

!!! note
    If you have suggestions to improve this guide, please make a pull request!
    It's particularly helpful if you do this after your first pull request
    because you'll know all the parts that could be explained better.

Thanks for contributing to JuMP!
# Style guide and design principles

## Style guide

This section describes the coding style rules that apply to JuMP code and that
we recommend for JuMP models and surrounding Julia code. The motivations for
a style guide include:

- conveying best practices for writing readable and maintainable code
- reducing the amount of time spent on
  [bike-shedding](https://en.wikipedia.org/wiki/Law_of_triviality) by
  establishing basic naming and formatting conventions
- lowering the barrier for new contributors by codifying the existing practices
  (e.g., you can be more confident your code will pass review if you follow the style guide)

In some cases, the JuMP style guide diverges from the
[Julia style guide](https://docs.julialang.org/en/v1.0.0/manual/style-guide/).
All such cases will be explicitly noted and justified.

The JuMP style guide adopts many recommendations from the
[Google style guides](https://github.com/google/styleguide).

!!! info
    The style guide is always a work in progress, and not all JuMP code
    follows the rules. When modifying JuMP, please fix the style violations
    of the surrounding code (i.e., leave the code tidier than when you
    started). If large changes are needed, consider separating them into
    another PR.

### JuliaFormatter

JuMP uses [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) as
an autoformatting tool.

We use the options contained in [`.JuliaFormatter.toml`](https://github.com/jump-dev/JuMP.jl/blob/master/.JuliaFormatter.toml).

To format code, `cd` to the JuMP directory, then run:
```julia
] add JuliaFormatter@0.13.2
using JuliaFormatter
format("src")
format("test")
```

!!! info
    A continuous integration check verifies that all PRs made to JuMP have
    passed the formatter.

The following sections outline extra style guide points that are not fixed
automatically by JuliaFormatter.

### Whitespace

For conciseness, never use more than one blank line within a function, and never
begin a function with a blank line.

Bad:
```julia
function foo(x)
    y = 2 * x


    return y
end

function foo(x)

    y = 2 * x
    return y
end
```

### Juxtaposed multiplication

Only use juxtaposed multiplication when the right-hand side is a symbol.

Good:
```julia
2x  # Acceptable if there are space constraints.
2 * x  # This is preferred if space is not an issue.
2 * (x + 1)
```

Bad:
```julia
2(x + 1)
```

### Empty vectors

For a type `T`, `T[]` and `Vector{T}()` are equivalent ways to create an
empty vector with element type `T`. Prefer `T[]` because it is more concise.

### Comments

For non-native speakers and for general clarity, comments in code must be proper
English sentences with appropriate punctuation.

Good:
```julia
# This is a comment demonstrating a good comment.
```

Bad:
```julia
# a bad comment
```

### JuMP macro syntax

For consistency, always use parentheses.

Good:
```julia
@variable(model, x >= 0)
```

Bad:
```julia
@variable model x >= 0
```

For consistency, always use `constant * variable` as opposed to
`variable * constant`. This makes it easier to read models in
ambiguous cases like `a * x`.

Good:
```julia
a = 4
@constraint(model, 3 * x <= 1)
@constraint(model, a * x <= 1)
```

Bad:
```julia
a = 4
@constraint(model, x * 3 <= 1)
@constraint(model, x * a <= 1)
```

In order to reduce boilerplate code, prefer the plural form of macros over lots
of repeated calls to singular forms.

Good:
```julia
@variables(model, begin
    x >= 0
    y >= 1
    z <= 2
end)
```

Bad:
```julia
@variable(model, x >= 0)
@variable(model, y >= 1)
@variable(model, z <= 2)
```

An exception is made for calls with many keyword arguments, since these need to
be enclosed in parentheses in order to parse properly.

Acceptable:
```julia
@variable(model, x >= 0, start = 0.0, base_name = "my_x")
@variable(model, y >= 1, start = 2.0)
@variable(model, z <= 2, start = -1.0)
```

Also acceptable:
```julia
@variables(model, begin
    x >= 0, (start = 0.0, base_name = "my_x")
    y >= 1, (start = 2.0)
    z <= 2, (start = -1.0)
end)
```

While we always use `in` for `for`-loops, it is acceptable to use `=` in the
container declarations of JuMP macros.

Okay:
```julia
@variable(model, x[i=1:3])
```
Also okay:
```julia
@variable(model, x[i in 1:3])
```

### Naming

```julia
module SomeModule end
function some_function end
const SOME_CONSTANT = ...
struct SomeStruct
  some_field::SomeType
end
@enum SomeEnum ENUM_VALUE_A ENUM_VALUE_B
some_local_variable = ...
some_file.jl # Except for ModuleName.jl.
```

### Exported and non-exported names

Begin private module level functions and constants with an underscore. All other
objects in the scope of a module should be exported. (See JuMP.jl for an example
of how to do this.)

Names beginning with an underscore should only be used for distinguishing
between exported (public) and non-exported (private) objects. Therefore, never
begin the name of a local variable with an underscore.

```julia
module MyModule

export public_function, PUBLIC_CONSTANT

function _private_function()
    local_variable = 1
    return
end

function public_function end

const _PRIVATE_CONSTANT = 3.14159
const PUBLIC_CONSTANT = 1.41421

end
```

### Use of underscores within names

The Julia style guide recommends avoiding underscores "when readable", for
example, `haskey`, `isequal`, `remotecall`, and `remotecall_fetch`. This
convention creates the potential for unnecessary bikeshedding and also forces
the user to recall the presence/absence of an underscore, e.g., "was that
argument named `basename` or `base_name`?". For consistency, *always use
underscores* in variable names and function names to separate words.

### Use of `!`

Julia has a convention of appending `!` to a function name if the function
modifies its arguments. We recommend to:

- Omit `!` when the name itself makes it clear that modification is taking
  place, e.g., `add_constraint` and `set_name`. We depart from the Julia style
  guide because `!` does not provide a reader with any additional information
  in this case, and adherence to this convention is not uniform even in base
  Julia itself (consider `Base.println` and `Base.finalize`).
- Use `!` in all other cases. In particular it can be used to distinguish
  between modifying and non-modifying variants of the same function like `scale`
  and `scale!`.

Note that `!` is *not* a self-documenting feature because it is still
ambiguous which arguments are modified when multiple arguments are present.
Be sure to document which arguments are modified in the method's docstring.

See also the Julia style guide recommendations for
[ordering of function arguments](https://docs.julialang.org/en/v1/manual/style-guide/#Write-functions-with-argument-ordering-similar-to-Julia-Base).

### Abbreviations

Abbreviate names to make the code more readable, not to save typing.
Don't arbitrarily delete letters from a word to abbreviate it (e.g., `indx`).
Use abbreviations consistently within a body of code (e.g., do not mix
`con` and `constr`, `idx` and `indx`).

Common abbreviations:

- `num` for `number`
- `con` for `constraint`

### No one-letter variable names

Where possible, avoid one-letter variable names.

Use `model = Model()` instead of `m = Model()`

Exceptions are made for indices in loops.

### User-facing `MethodError`

Specifying argument types for methods is mostly optional in Julia, which means
that it's possible to find out that you are working with unexpected types deep in
the call chain. Avoid this situation or handle it with a helpful error message.
*A user should see a `MethodError` only for methods that they called directly.*

Bad:
```julia
_internal_function(x::Integer) = x + 1
# The user sees a MethodError for _internal_function when calling
# public_function("a string"). This is not very helpful.
public_function(x) = _internal_function(x)
```

Good:
```julia
_internal_function(x::Integer) = x + 1
# The user sees a MethodError for public_function when calling
# public_function("a string"). This is easy to understand.
public_function(x::Integer) = _internal_function(x)
```

If it is hard to provide an error message at the top of the call chain,
then the following pattern is also ok:
```julia
_internal_function(x::Integer) = x + 1
function _internal_function(x)
    error(
        "Internal error. This probably means that you called " *
        "public_function() with the wrong type.",
    )
end
public_function(x) = _internal_function(x)
```

### `@enum` vs. `Symbol`

The `@enum` macro lets you define types with a finite number of values that
are explicitly enumerated (like `enum` in C/C++). `Symbol`s are lightweight
strings that are used to represent identifiers in Julia (for example, `:x`).

`@enum` provides type safety and can have docstrings attached to explain the
possible values. Use `@enum`s when applicable, e.g., for reporting statuses.
Use strings to provide long-form additional information like error messages.

Use of `Symbol` should typically be reserved for identifiers, e.g., for lookup
in the JuMP model (`model[:my_variable]`).

### `using` vs. `import`

`using ModuleName` brings all symbols exported by the module `ModuleName`
into scope, while `import ModuleName` brings only the module itself into scope.
(See the Julia
[manual](https://docs.julialang.org/en/v1/manual/modules/#modules-1)) for
examples and more details.

For the same reason that `from <module> import *` is not recommended in python
([PEP 8](https://www.python.org/dev/peps/pep-0008/#imports)), avoid
`using ModuleName` except in throw-away scripts or at the REPL. The `using`
statement makes it harder to track where symbols come from and exposes the code
to ambiguities when two modules export the same symbol.

Prefer `using ModuleName: x, p` to `import ModuleName.x, ModuleName.p` and
`import MyModule: x, p` because the `import` versions allow method extension
without qualifying with the module name.

Similarly, `using ModuleName: ModuleName` is an acceptable substitute for
`import ModuleName`, because it does not bring all symbols exported by
`ModuleName` into scope. However, we prefer `import ModuleName` for consistency.

## Documentation

This section describes the writing style that should be used when writing
documentation for JuMP (and supporting packages).

We can recommend the documentation style guides by [Divio](https://www.divio.com/blog/documentation/),
[Google](https://developers.google.com/style/), and [Write the Docs](https://www.writethedocs.org/guide/)
as general reading for those writing documentation. This guide delegates a
thorough handling of the topic to those guides and instead elaborates on the
points more specific to Julia and documentation that use [Documenter](https://github.com/JuliaDocs/Documenter.jl).

 - Be concise
 - Use lists instead of long sentences
 - Use numbered lists when describing a sequence, e.g., (1) do X, (2) then Y
 - Use bullet points when the items are not ordered
 - Example code should be covered by doctests
 - When a word is a Julia symbol and not an English word, enclose it with
   backticks. In addition, if it has a docstring in this doc add a link using
   `@ref`. If it is a plural, add the "s" after the closing backtick. For example,
   ```
   [`VariableRef`](@ref)s
   ```
 - Use [`@meta`](https://juliadocs.github.io/Documenter.jl/v0.21/man/syntax/#@meta-block-1)
   blocks for TODOs and other comments that shouldn't be visible to readers.
   For example,
   ````@markdown
   ```@meta
   # TODO: Mention also X, Y, and Z.
   ```
   ````
### Docstrings

- Every exported object needs a docstring
- All examples in docstrings should be [`jldoctests`](https://juliadocs.github.io/Documenter.jl/stable/man/doctests/)
- Always use complete English sentences with proper punctuation
- Do not terminate lists with punctuation (e.g., as in this doc)

Here is an example:
````julia
"""
    signature(args; kwargs...)

Short sentence describing the function.

Optional: add a slightly longer paragraph describing the function.

## Notes

 - List any notes that the user should be aware of

## Examples

```jldoctest
julia> 1 + 1
2
```
"""
````

## Testing

Use a module to encapsulate tests, and structure all tests as functions. This
avoids leaking local variables between tests.

Here is a basic skeleton:
```julia
module TestPkg

using Test

_helper_function() = 2

function test_addition()
    @test 1 + 1 == _helper_function()
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

end # TestPkg

TestPkg.runtests()
```

Break the tests into multiple files, with one module per file, so that subsets
of the codebase can be tested by calling `include` with the relevant file.

## Design principles

TODO: How to structure and test large JuMP models, libraries that use JuMP.

For how to write a solver, see MOI.
# Development roadmap

This page is not JuMP documentation *per se* but are notes for the JuMP
community. The JuMP developers have compiled this roadmap document to
share their plans and goals. Contributions to roadmap issues are especially
invited.

## JuMP 1.0

JuMP 1.0 will be ready to release roughly when all of these tasks are completed.
Some but not all of these tasks are summarized in the
[JuMP 1.0 milestone](https://github.com/jump-dev/JuMP.jl/milestone/12).

- Create a website for JuMP (**Done**: [jump.dev](https://jump.dev))
- Deprecate the JuliaOpt organization and move repositories to the
  [JuMP-dev](https://github.com/JuMP-dev) organization (**Done**)
- Address major regressions from JuMP 0.18
  - Performance ([#1403](https://github.com/jump-dev/JuMP.jl/issues/1403),
                 [#1654](https://github.com/jump-dev/JuMP.jl/issues/1654),
                 [#1607](https://github.com/jump-dev/JuMP.jl/issues/1607))
  - Callbacks (**Done**: see `examples/callbacks.jl`)
  - Column generation syntax (**Done**: see `examples/cutting_stock_column_generation.jl`)
  - Support for second-order cones in Gurobi, CPLEX, and Xpress (**Done**)
- Fix issues that we promised MOI would fix
  - Checking feasibility of solutions (**Done**: [#2466](https://github.com/jump-dev/JuMP.jl/pull/2466))
  - Accessing IIS (**Done**: see [Conflicts](@ref))
  - Accessing multiple results from solvers (**Done**: [Gurobi#392](https://github.com/jump-dev/Gurobi.jl/pull/392))
  - Dual warm-starts (**Done**: [#2214](https://github.com/jump-dev/JuMP.jl/pull/2214))
- Address "easy" usability issues
  - Line numbers in error messages (**Done**: [#2276](https://github.com/jump-dev/JuMP.jl/pull/2276))
  - LP sensitivity summary (**Done**: see [Sensitivity analysis for LP](@ref))
  - Inferred element types for collections in macros (**Done**: [#2070](https://github.com/jump-dev/JuMP.jl/pull/2070))
  - Expose solver-independent options from JuMP (**Done**: see [`set_silent`](@ref) etc.)
- Improve the documentation ([#1062](https://github.com/jump-dev/JuMP.jl/issues/1062))
  - Separate how-to, concept explanation, and technical reference following the
    [Divio recommendations](https://www.divio.com/blog/documentation/) (**Done**)
  - Fully integrate [JuMPTutorials](https://github.com/jump-dev/JuMPTutorials.jl)
    with JuMP's documentation (**Done**)
- Developer experience
  - Get JuMP's unit tests running faster. See [#1745](https://github.com/jump-dev/JuMP.jl/pull/1745). (**Done**)
- All solvers should complete the transition to MOI (**Done**)
- Provide packages for installing Bonmin and Couenne (**Done**)
- [MathOptFormat](https://github.com/odow/MathOptFormat.jl) 1.0 (**Done**)

## MOI 1.0

```@meta
# TODO: List MOI 1.0 items here.
```

## Beyond JuMP 1.0

```@meta
# TODO: Copy over list of items not tied to JuMP 1.0. These should have more
# elaborate explanations so that potential contributors know what we mean,
# i.e., a few sentences each or a link to a document/issue.
```
```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP, GLPK
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# [Algebraic modeling languages](@id algebraic-modeling-language)

JuMP is an algebraic modeling language for mathematical optimization written in
the [Julia language](https://julialang.org). In this page, we explain what an
algebraic modeling language actually is.

## What is an algebraic modeling language?

If you have taken a class in mixed-integer linear programming, you will have
seen a formulation like:
```math
\begin{aligned}
\min \;     & c^\top x \\
\text{s.t.} & A x = b  \\
            & x \ge 0  \\
            & x_i \in \mathbb{Z}, \quad \forall i \in \mathcal{I}
\end{aligned}
```
where `c`, `A`, and `b` are appropriately sized vectors and matrices of data,
and $\mathcal{I}$ denotes the set of variables that are integer.

Solvers expect problems in a _standard form_ like this because it limits the
types of constraints that they need to consider. This makes writing a solver
much easier.

!!! info "What is a solver?"
    A solver is a software package that computes solutions to one or more
    classes of problems.

    For example, [GLPK](https://www.gnu.org/software/glpk/) is a solver for
    linear programming (LP) and mixed integer programming (MIP) problems. It
    incorporates algorithms such as the simplex method and the interior-point
    method.

    JuMP currently supports a number of open-source and commercial solvers,
    which can be viewed in the [Supported-solvers](@ref) table.


Despite the textbook view of a linear program, you probably formulated problems
algebraically like so:
```math
\begin{aligned}
\max \;     & \sum\limits_{i = 1}^n c_i x_i                   \\
\text{s.t.} & \sum\limits_{i = 1}^n w_i x_i \le b             \\
            & x_i \ge 0 \quad \forall i = 1,\ldots,n          \\
            & x_i \in \mathbb{Z} \quad \forall i = 1,\ldots,n.
\end{aligned}
```
!!! info
    Do you recognize this formulation? It's the knapsack problem.

Users prefer to write problems in _algebraic form_ because it is more
convenient. For example, we used $\le b$, even though the standard form
only supported constraints of the form $Ax = b$.

We could convert our knapsack problem into the standard form by adding a new
slack variable $x_0$:
```math
\begin{aligned}
\max \;     & \sum\limits_{i = 1}^n c_i x_i            \\
\text{s.t.} & x_0 + \sum\limits_{i = 1}^n w_i x_i = b  \\
            & x_i \ge 0 \quad \forall i = 0,\ldots,n   \\
            & x_i \in \mathbb{Z} \quad \forall i = 1,\ldots,n.
\end{aligned}
```
However, as models get more complicated, this manual conversion becomes more and
more error-prone.

An algebraic modeling language is a tool that simplifies the translation between
the algebraic form of the modeler, and the standard form of the solver.

Each algebraic modeling language has two main parts:

 1. A domain specific language for the user to write down problems in algebraic
    form.
 2. A converter from the algebraic form into a standard form supported by the
    solver (and back again).

Part 2 is less trivial than it might seem, because each solver has a unique
application programming interface (API) and data structure for representing
optimization models and obtaining results.

JuMP uses the
[MathOptInterface.jl](https://github.com/jump-dev/MathOptInterface.jl)
package to abstract these differences between solvers.

### What is MathOptInterface?

MathOptInterface (MOI) is an abstraction layer designed to provide an
interface to mathematical optimization solvers so that users do not need to
understand multiple solver-specific APIs. MOI can be used directly, or through
a higher-level modeling interface like JuMP.

There are three main parts to MathOptInterface:

 1. A solver-independent API that abstracts concepts such as adding and deleting
    variables and constraints, setting and getting parameters, and querying
    results. For more information on the MathOptInterface API, read the
    [documentation](@ref moi_documentation).

 2. An automatic rewriting system based on equivalent formulations of a
    constraint. For more information on this rewriting system, read the
    [LazyBridgeOptimizer](@ref) section of the manual, and our
    [paper on arXiv](https://arxiv.org/abs/2002.03447).

 3. Utilities for managing how and when models are copied to solvers. For more
    information on this, read the [CachingOptimizer](@ref) section of the
    manual.

## From user to solver

This section provides a brief summary of the steps that happen in order to
translate the model that the user writes into a model that the solver
understands.

### Step I: writing in algebraic form

JuMP provides the first part of an algebraic modeling language using the
[`@variable`](@ref), [`@objective`](@ref), and [`@constraint`](@ref) macros.

For example, here's how we write the knapsack problem in JuMP:
```jldoctest
julia> using JuMP, GLPK

julia> function algebraic_knapsack(c, w, b)
           n = length(c)
           model = Model(GLPK.Optimizer)
           @variable(model, x[1:n] >= 0, Int)
           @objective(model, Max, sum(c[i] * x[i] for i = 1:n))
           @constraint(model, sum(w[i] * x[i] for i = 1:n) <= b)
           optimize!(model)
           return value.(x)
       end
algebraic_knapsack (generic function with 1 method)

julia> algebraic_knapsack([1, 2], [0.5, 0.5], 1.25)
2-element Vector{Float64}:
 0.0
 2.0
```
This formulation is compact, and it closely matches the algebraic formulation of
the model we wrote out above.

### Step II: algebraic to functional

For the next step, JuMP's macros re-write the variables and constraints into a
functional form. Here's what the JuMP code looks like after this step:
```jldoctest
julia> using JuMP, GLPK

julia> function nonalgebraic_knapsack(c, w, b)
           n = length(c)
           model = Model(GLPK.Optimizer)
           x = [VariableRef(model) for i = 1:n]
           for i = 1:n
               set_lower_bound(x[i], 0)
               set_integer(x[i])
               set_name(x[i], "x[$i]")
           end
           obj = AffExpr(0.0)
           for i = 1:n
               add_to_expression!(obj, c[i], x[i])
           end
           set_objective(model, MAX_SENSE, obj)
           lhs = AffExpr(0.0)
           for i = 1:n
               add_to_expression!(lhs, w[i], x[i])
           end
           con = build_constraint(error, lhs, MOI.LessThan(b))
           add_constraint(model, con)
           optimize!(model)
           return value.(x)
       end
nonalgebraic_knapsack (generic function with 1 method)

julia> nonalgebraic_knapsack([1, 2], [0.5, 0.5], 1.25)
2-element Vector{Float64}:
 0.0
 2.0
```

Hopefully you agree that the macro version is much easier to read!

### Part III: JuMP to MathOptInterface

In the third step, JuMP converts the functional form of the problem, i.e.,
`nonalgebraic_knapsack`, into the MathOptInterface API:
```jldoctest
julia> using MathOptInterface, GLPK

julia> const MOI = MathOptInterface;

julia> function mathoptinterface_knapsack(optimizer, c, w, b)
           n = length(c)
           model = MOI.instantiate(optimizer)
           x = MOI.add_variables(model, n)
           for i in 1:n
               MOI.add_constraint(model, x[i], MOI.GreaterThan(0.0))
               MOI.add_constraint(model, x[i], MOI.Integer())
               MOI.set(model, MOI.VariableName(), x[i], "x[$i]")
           end
           MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
           obj = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(c, x), 0.0)
           MOI.set(model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
           MOI.add_constraint(
               model,
               MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(w, x), 0.0),
               MOI.LessThan(b),
           )
           MOI.optimize!(model)
           return MOI.get.(model, MOI.VariablePrimal(), x)
       end
mathoptinterface_knapsack (generic function with 1 method)

julia> mathoptinterface_knapsack(GLPK.Optimizer, [1.0, 2.0], [0.5, 0.5], 1.25)
2-element Vector{Float64}:
 0.0
 2.0
```
The code is becoming more verbose and looking less like the mathematical
formulation that we started with.

### Step IV: MathOptInterface to GLPK

As a final step, the [GLPK.jl](https://github.com/jump-dev/GLPK.jl) package
converts the MathOptInterface form, i.e., `mathoptinterface_knapsack`, into a
GLPK-specific API:
```jldoctest
julia> using GLPK

julia> function glpk_knapsack(c, w, b)
           n = length(c)
           model = glp_create_prob()
           glp_add_cols(model, n)
           for i in 1:n
               glp_set_col_bnds(model, i, GLP_LO, 0.0, GLP_DBL_MAX)
               glp_set_col_kind(model, i, GLP_IV)
               glp_set_col_name(model, i, "x[$i]")
           end
           glp_set_obj_dir(model, GLP_MAX)
           for i in 1:n
               glp_set_obj_coef(model, i, c[i])
           end
           glp_set_obj_coef(model, 0, 0.0)
           glp_add_rows(model, 1)
           glp_set_mat_row(
               model,
               1,
               length(w),
               Ref(Cint.(1:n), 0),
               Ref(w, 0),
           )
           glp_set_row_bnds(model, 1, GLP_UP, -GLP_DBL_MAX, b)
           simplex_options = glp_smcp()
           glp_init_smcp(simplex_options)
           simplex_options.msg_lev = GLP_MSG_ERR
           glp_simplex(model, simplex_options)
           options = glp_iocp()
           glp_init_iocp(options)
           options.msg_lev = GLP_MSG_ERR
           glp_intopt(model, options)
           x = glp_mip_col_val.(model, 1:n)
           glp_delete_prob(model)
           return x
       end
glpk_knapsack (generic function with 1 method)

julia> glpk_knapsack([1.0, 2.0], [0.5, 0.5], 1.25)
2-element Vector{Float64}:
 0.0
 2.0
```
We've now gone from a algebraic model that looked identical to the mathematical
model we started with, to a verbose function that uses GLPK-specific
functionality.

The difference between `algebraic_knapsack` and `glpk_knapsack` highlights the
benefit that algebraic modeling languages provide to users. Moreover, if we used
a different solver, the solver-specific function would be entirely different. A
key benefit of an algebraic modeling language is that you can change the solver
without needing to rewrite the model.
```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# [Solver-independent Callbacks](@id callbacks_manual)

Many mixed-integer (linear, conic, and nonlinear) programming solvers offer
the ability to modify the solve process. Examples include changing branching
decisions in branch-and-bound, adding custom cutting planes, providing custom
heuristics to find feasible solutions, or implementing on-demand separators to
add new constraints only when they are violated by the current solution (also
known as lazy constraints).

While historically this functionality has been limited to solver-specific
interfaces, JuMP provides solver-independent support for three types of
callbacks:

 1. lazy constraints
 2. user-cuts
 3. heuristic solutions

## Available solvers

Solver-independent callback support is limited to a few solvers. This includes
[CPLEX](https://github.com/jump-dev/CPLEX.jl),
[GLPK](https://github.com/jump-dev/GLPK.jl),
[Gurobi](https://github.com/jump-dev/Gurobi.jl), and
[Xpress](https://github.com/jump-dev/Xpress.jl).

!!! warning
    While JuMP provides a solver-independent way of accessing callbacks, you
    should not assume that you will see identical behavior when running the same
    code on different solvers. For example, some solvers may ignore user-cuts
    for various reasons, while other solvers may add every user-cut. Read the
    underlying solver's callback documentation to understand details specific to
    each solver.

!!! tip
    This page discusses solver-_independent_ callbacks. However, each solver
    listed above also provides a solver-_dependent_ callback to provide access
    to the full range of solver-specific features. Consult the solver's README
    for an example of how to use the solver-dependent callback. This will
    require you to understand the C interface of the solver.

## Things you can and cannot do during solver-independent callbacks

There is a limited range of things you can do during a callback. Only use the
functions and macros explicitly stated in this page of the documentation, or in
the [Callbacks tutorial](@ref callbacks_tutorial).

Using any other part of the JuMP API (for example, adding a constraint with [`@constraint`](@ref)
or modifying a variable bound with [`set_lower_bound`](@ref)) is undefined
behavior, and your solver may throw an error, return an incorrect solution, or
result in a segfault that aborts Julia.

In each of the three solver-independent callbacks, there are two things you may
query:
 - [`callback_node_status`](@ref) returns an [`MOI.CallbackNodeStatusCode`](@ref)
   enum indicating if the current primal solution is integer feasible.
 - [`callback_value`](@ref) returns the current primal solution of a variable.

If you need to query any other information, use a solver-dependent callback
instead. Each solver supporting a solver-dependent callback has information on
how to use it in the README of their GitHub repository.

If you want to modify the problem in a callback, you _must_ use a lazy
constraint.

!!! warning
    You can only set each callback once. Calling `set` twice will over-write
    the earlier callback. In addition, if you use a solver-independent
    callback, you cannot set a solver-dependent callback.

## Lazy constraints

Lazy constraints are useful when the full set of constraints is too large to
explicitly include in the initial formulation. When a MIP solver reaches a new
solution, for example with a heuristic or by solving a problem at a node in
the branch-and-bound tree, it will give the user the chance to provide
constraints that would make the current solution infeasible. For some more
information about lazy constraints, see this [blog post by Paul Rubin](https://orinanobworld.blogspot.com/2012/08/user-cuts-versus-lazy-constraints.html).

A lazy constraint callback can be set using the following syntax:

```julia
model = Model(GLPK.Optimizer)
@variable(model, x <= 10, Int)
@objective(model, Max, x)
function my_callback_function(cb_data)
    status = callback_node_status(cb_data, model)
    if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        # `callback_value(cb_data, x)` is not integer (to some tolerance).
        # If, for example, your lazy constraint generator requires an
        # integer-feasible primal solution, you can add a `return` here.
        return
    elseif status == MOI.CALLBACK_NODE_STATUS_INTEGER
        # `callback_value(cb_data, x)` is integer (to some tolerance).
    else
        @assert status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
        # `callback_value(cb_data, x)` might be fractional or integer.
    end
    x_val = callback_value(cb_data, x)
    if x_val > 2 + 1e-6
        con = @build_constraint(x <= 2)
        MOI.submit(model, MOI.LazyConstraint(cb_data), con)
    end
end
MOI.set(model, MOI.LazyConstraintCallback(), my_callback_function)
```

!!! info
    The lazy constraint callback _may_ be called at fractional or integer
    nodes in the branch-and-bound tree. There is no guarantee that the
    callback is called at _every_ primal solution.

!!! warning
    Only add a lazy constraint if your primal solution violates the constraint.
    Adding the lazy constraint irrespective of feasibility may result in the
    solver returning an incorrect solution, or lead to many constraints being
    added, slowing down the solution process.
    ```julia
    model = Model(GLPK.Optimizer)
    @variable(model, x <= 10, Int)
    @objective(model, Max, x)
    function bad_callback_function(cb_data)
        # Don't do this!
        con = @build_constraint(x <= 2)
        MOI.submit(model, MOI.LazyConstraint(cb_data), con)
    end
    function good_callback_function(cb_data)
        if callback_value(x) > 2
            con = @build_constraint(x <= 2)
            MOI.submit(model, MOI.LazyConstraint(cb_data), con)
        end
    end
    MOI.set(model, MOI.LazyConstraintCallback(), good_callback_function)
    ```

!!! warning
    During the solve, a solver may visit a point that was cut off by a previous
    lazy constraint, for example, because the earlier lazy constraint was removed
    during presolve. However, the solver will not stop until it reaches a
    solution that satisfies all added lazy constraints.

## User cuts

User cuts, or simply cuts, provide a way for the user to tighten the LP
relaxation using problem-specific knowledge that the solver cannot or is
unable to infer from the model. Just like with lazy constraints, when a MIP
solver reaches a new node in the branch-and-bound tree, it will give the user
the chance to provide cuts to make the current relaxed (fractional) solution
infeasible in the hopes of obtaining an integer solution. For more details
about the difference between user cuts and lazy constraints see the
aforementioned [blog post](https://orinanobworld.blogspot.com/2012/08/user-cuts-versus-lazy-constraints.html).

A user-cut callback can be set using the following syntax:

```julia
model = Model(GLPK.Optimizer)
@variable(model, x <= 10.5, Int)
@objective(model, Max, x)
function my_callback_function(cb_data)
    x_val = callback_value(cb_data, x)
    con = @build_constraint(x <= floor(x_val))
    MOI.submit(model, MOI.UserCut(cb_data), con)
end
MOI.set(model, MOI.UserCutCallback(), my_callback_function)
```

!!! warning
    User cuts must not change the set of integer feasible solutions.
    Equivalently, user cuts can only remove fractional solutions. If you add a
    cut that removes an integer solution (even one that is not optimal), the
    solver may return an incorrect solution.

!!! info
    The user-cut callback _may_ be called at fractional nodes in the
    branch-and-bound tree. There is no guarantee that the callback is called
    at _every_ fractional primal solution.

## Heuristic solutions

Integer programming solvers frequently include heuristics that run at the
nodes of the branch-and-bound tree. They aim to find integer solutions quicker
than plain branch-and-bound would to tighten the bound, allowing us to fathom
nodes quicker and to tighten the integrality gap.

Some heuristics take integer solutions and explore their "local neighborhood"
(for example, flipping binary variables, fix some variables and solve a smaller MILP)
and others take fractional solutions and attempt to round them in an
intelligent way.

You may want to add a heuristic of your own if you have some special insight
into the problem structure that the solver is not aware of, for example, you can
consistently take fractional solutions and intelligently guess integer
solutions from them.

A heuristic solution callback can be set using the following syntax:

```julia
model = Model(GLPK.Optimizer)
@variable(model, x <= 10.5, Int)
@objective(model, Max, x)
function my_callback_function(cb_data)
    x_val = callback_value(cb_data, x)
    status = MOI.submit(
        model, MOI.HeuristicSolution(cb_data), [x], [floor(Int, x_val)]
    )
    println("I submitted a heuristic solution, and the status was: ", status)
end
MOI.set(model, MOI.HeuristicCallback(), my_callback_function)
```

The third argument to `submit` is a vector of JuMP variables, and the
fourth argument is a vector of values corresponding to each variable.

`MOI.submit` returns an enum that depends on whether the solver accepted the
solution. The possible return codes are:

 - `MOI.HEURISTIC_SOLUTION_ACCEPTED`
 - `MOI.HEURISTIC_SOLUTION_REJECTED`
 - `MOI.HEURISTIC_SOLUTION_UNKNOWN`

!!! warning
    Some solvers may accept partial solutions. Others require a feasible integer
    solution for every variable. If in doubt, provide a complete solution.

!!! info
    The heuristic solution callback _may_ be called at fractional nodes in the
    branch-and-bound tree. There is no guarantee that the callback is called
    at _every_ fractional primal solution.
```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP, GLPK, SCS
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# [Models](@id jump_models)

JuMP models are the fundamental building block that we use to construct
optimization problems. They hold things like the variables and constraints, as
well as which solver to use and even solution information.

!!! info
    JuMP uses "optimizer" as a synonym for "solver." Our convention is to use
    "solver" to refer to the underlying software, and use "optimizer" to refer
    to the Julia object that wraps the solver. For example, `GLPK` is a solver,
    and `GLPK.Optimizer` is an optimizer.

!!! tip
    See [Supported solvers](@ref) for a list of available solvers.

## Create a model

Create a model by passing an optimizer to [`Model`](@ref):
```jldoctest
julia> model = Model(GLPK.Optimizer)
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: GLPK
```

If you don't know which optimizer you will be using at creation time, create a
model without an optimizer, and then call [`set_optimizer`](@ref) at any time
prior to [`optimize!`](@ref):
```jldoctest
julia> model = Model()
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> set_optimizer(model, GLPK.Optimizer)
```

!!! tip
    Don't know what the fields `Model mode` and `CachingOptimizer state` mean?
    Read the [Backends](@ref) section.

### What is the difference?

For most models, there is no difference between passing the optimizer to
[`Model`](@ref), and calling [`set_optimizer`](@ref).

However, if an optimizer does not support a constraint in the model, the timing
of when an error will be thrown can differ:

 * If you pass an optimizer, an error will be thrown when you try to add the
   constraint.
 * If you call [`set_optimizer`](@ref), an error will be thrown when you try to
   solve the model via [`optimize!`](@ref).

Therefore, most users should pass an optimizer to [`Model`](@ref) because it
provides the earliest warning that your solver is not suitable for the model you
are trying to build. However, if you are modifying a problem by adding and
deleting different constraint types, you may need to use
[`set_optimizer`](@ref). See [Switching optimizer for the relaxed problem](@ref)
for an example of when this is useful.

### Reducing time-to-first-solve latency

By default, JuMP uses [bridges](@ref LazyBridgeOptimizer) to reformulate the
model you are building into an equivalent model supported by the solver.

However, if your model is already supported by the solver, bridges add latency
(read [The "time-to-first-solve" issue](@ref)). This is particularly noticeable
for small models.

To reduce the "time-to-first-solve", try passing `add_bridges = false`.
```jldoctest
julia> model = Model(GLPK.Optimizer; add_bridges = false);
```
or
```jldoctest
julia> model = Model();

julia> set_optimizer(model, GLPK.Optimizer; add_bridges = false)
```

However, be wary! If your model and solver combination needs bridges, an error
will be thrown:
```jldoctest
julia> model = Model(SCS.Optimizer; add_bridges = false);


julia> @variable(model, x)
x

julia> @constraint(model, 2x <= 1)
ERROR: Constraints of type MathOptInterface.ScalarAffineFunction{Float64}-in-MathOptInterface.LessThan{Float64} are not supported by the solver.

If you expected the solver to support your problem, you may have an error in your formulation. Otherwise, consider using a different solver.

The list of available solvers, along with the problem types they support, is available at https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers.
[...]
```

### Solvers which expect environments

Some solvers accept (or require) positional arguments such as a license
environment or a path to a binary executable. For these solvers, you can pass
a function to [`Model`](@ref) which takes zero arguments and returns an instance
of the optimizer.

A common use-case for this is passing an environment to Gurobi:
```julia
julia> grb_env = Gurobi.Env();

julia> model = Model(() -> Gurobi.Optimizer(grb_env))
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: Gurobi
```

## [Solver options](@id solver_options)

JuMP uses "attribute" as a synonym for "option." Use
[`optimizer_with_attributes`](@ref) to create an optimizer with some attributes
initialized:
```jldoctest
julia> model = Model(optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => 0))
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: GLPK
```

Alternatively, use [`set_optimizer_attribute`](@ref) to set an attribute after
the model has been created:
```jldoctest
julia> model = Model(GLPK.Optimizer);

julia> set_optimizer_attribute(model, "msg_lev", 0)

julia> get_optimizer_attribute(model, "msg_lev")
0
```

## Print the model

By default, `show(model)` will print a summary of the problem:
```jldoctest model_print
julia> model = Model(); @variable(model, x >= 0); @objective(model, Max, x);

julia> model
A JuMP Model
Maximization problem with:
Variable: 1
Objective function type: VariableRef
`VariableRef`-in-`MathOptInterface.GreaterThan{Float64}`: 1 constraint
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
Names registered in the model: x
```

Use `print` to print the formulation of the model (in IJulia, this will render
as LaTeX.
```jldoctest model_print
julia> print(model)
Max x
Subject to
 x ≥ 0.0
```

!!! warning
    This format is specific to JuMP and may change in any future release. It is
    not intended to be an instance format. To write the model to a file, use
    [`write_to_file`](@ref) instead.


Use [`latex_formulation`](@ref) to display the model in LaTeX form.
```jldoctest model_print
julia> latex_formulation(model)
$$ \begin{aligned}
\max\quad & x\\
\text{Subject to} \quad & x \geq 0.0\\
\end{aligned} $$
```

In IJulia (and Documenter), ending a cell in with [`latex_formulation`](@ref)
will render the model in LaTeX!
```@example
using JuMP                # hide
model = Model()           # hide
@variable(model, x >= 0)  # hide
@objective(model, Max, x) # hide
latex_formulation(model)
```

## Turn off output

Use [`set_silent`](@ref) and [`unset_silent`](@ref) to disable or enable
printing output from the solver.
```jldoctest
julia> model = Model(GLPK.Optimizer);

julia> set_silent(model)

julia> unset_silent(model)
```

!!! tip
    Most solvers will also have a [solver-specific option](@ref solver_options)
    to provide finer-grained control over the output. Consult their README's for
    details.

## Set a time limit

Use [`set_time_limit_sec`](@ref), [`unset_time_limit_sec`](@ref), and
[`time_limit_sec`](@ref) to manage time limits.
```jldoctest
julia> model = Model(GLPK.Optimizer);

julia> set_time_limit_sec(model, 60.0)


julia> time_limit_sec(model)
60.0

julia> unset_time_limit_sec(model)

julia> time_limit_sec(model)
2.147483647e6
```

!!! info
    Some solvers do not support time limits. In these cases, an error will be
    thrown.

## Write a model to file

JuMP can write models to a variety of file-formats using [`write_to_file`](@ref)
and [`Base.write`](@ref).

For most common file formats, the file type will be detected from the extension.
For example, here is how to write an MPS file:
```jldoctest file_formats; setup=:(model = Model())
julia> write_to_file(model, "model.mps")
```

To write to a specific `io::IO`, use [`Base.write`](@ref). Specify the file type
by passing a [`MOI.FileFormats.FileFormat`](@ref) enum.
```jldoctest file_formats; setup=:(model = Model(); io = IOBuffer())
julia> write(io, model; format = MOI.FileFormats.FORMAT_MPS)
```

## Read a model from file

JuMP models can be created from file formats using [`read_from_file`](@ref) and
[`Base.read`](@ref).

```jldoctest file_formats
julia> model = read_from_file("model.mps")
A JuMP Model
Minimization problem with:
Variables: 0
Objective function type: AffExpr
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> seekstart(io);

julia> model2 = read(io, Model; format = MOI.FileFormats.FORMAT_MPS)
A JuMP Model
Minimization problem with:
Variables: 0
Objective function type: AffExpr
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
```

!!! note
    Because file formats do not serialize the containers of JuMP variables and
    constraints, the names in the model will _not_ be registered. Therefore, you
    cannot access named variables and constraints via `model[:x]`. Instead, use
    [`variable_by_name`](@ref) or [`constraint_by_name`](@ref) to access
    specific variables or constraints.

## Relax integrality

Use [`relax_integrality`](@ref) to remove any integrality constraints from the
model, such as integer and binary restrictions on variables.
[`relax_integrality`](@ref) returns a function that can be later called with
zero arguments to re-add the removed constraints:
```jldoctest
julia> model = Model();

julia> @variable(model, x, Int)
x

julia> num_constraints(model, VariableRef, MOI.Integer)
1

julia> undo = relax_integrality(model);

julia> num_constraints(model, VariableRef, MOI.Integer)
0

julia> undo()

julia> num_constraints(model, VariableRef, MOI.Integer)
1
```

### Switching optimizer for the relaxed problem

A common reason for relaxing integrality is to compute dual variables of the
relaxed problem. However, some mixed-integer linear solvers (for example, Cbc) do not
return dual solutions, even if the problem does not have integrality
restrictions.

Therefore, after [`relax_integrality`](@ref) you should call
[`set_optimizer`](@ref) with a solver that does support dual solutions, such as
Clp.

For example, instead of:
```julia
using JuMP, Cbc
model = Model(Cbc.Optimizer)
@variable(model, x, Int)
undo = relax_integrality(model)
optimize!(model)
reduced_cost(x)  # Errors
```
do:
```julia
using JuMP, Cbc, Clp
model = Model(Cbc.Optimizer)
@variable(model, x, Int)
undo = relax_integrality(model)
set_optimizer(model, Clp.Optimizer)
optimize!(model)
reduced_cost(x)  # Works
```

## Backends

!!! info
    This section discusses advanced features of JuMP. For new users, you may
    want to skip this section. You don't need to know how JuMP manages problems
    behind the scenes to create and solve JuMP models.

A JuMP [`Model`](@ref) is a thin layer around a *backend* of type
[`MOI.ModelLike`](@ref) that stores the optimization problem and acts as the
optimization solver.

However, if you construct a model like `Model(GLPK.Optimizer)`, the backend is
not a `GLPK.Optimizer`, but a more complicated object.

From JuMP, the MOI backend can be accessed using the [`backend`](@ref) function.
Let's see what the [`backend`](@ref) of a JuMP [`Model`](@ref) is:
```jldoctest models_backends
julia> model = Model(GLPK.Optimizer);

julia> b = backend(model)
MOIU.CachingOptimizer{MOIB.LazyBridgeOptimizer{GLPK.Optimizer}, MOIU.UniversalFallback{MOIU.Model{Float64}}}
in state EMPTY_OPTIMIZER
in mode AUTOMATIC
with model cache MOIU.UniversalFallback{MOIU.Model{Float64}}
  fallback for MOIU.Model{Float64}
with optimizer MOIB.LazyBridgeOptimizer{GLPK.Optimizer}
  with 0 variable bridges
  with 0 constraint bridges
  with 0 objective bridges
  with inner model A GLPK model
```

Uh oh! Even though we passed a `GLPK.Optimizer`, the backend is a much more
complicated object.

### CachingOptimizer

A `MOIU.CachingOptimizer` is a layer that abstracts the difference between
solvers that support incremental modification (for example, they support adding
variables one-by-one), and solvers that require the entire problem in a single
API call (for example, they only accept the `A`, `b` and `c` matrices of a linear
program).

It has two parts:

 1. A cache, where the model can be built and modified incrementally
    ```jldoctest models_backends
    julia> b.model_cache
    MOIU.UniversalFallback{MOIU.Model{Float64}}
    fallback for MOIU.Model{Float64}
    ```
 2. An optimizer, which is used to solve the problem
    ```jldoctest models_backends
    julia> b.optimizer
    MOIB.LazyBridgeOptimizer{GLPK.Optimizer}
    with 0 variable bridges
    with 0 constraint bridges
    with 0 objective bridges
    with inner model A GLPK model
    ```

!!! info
    The [LazyBridgeOptimizer](@ref) section explains what a
    `LazyBridgeOptimizer` is.

The `CachingOptimizer` has logic to decide when to copy the problem from the
cache to the optimizer, and when it can efficiently update the optimizer
in-place.

A `CachingOptimizer` may be in one of three possible states:

* `NO_OPTIMIZER`: The CachingOptimizer does not have any optimizer.
* `EMPTY_OPTIMIZER`: The CachingOptimizer has an empty optimizer, and it is not
  synchronized with the cached model.
* `ATTACHED_OPTIMIZER`: The CachingOptimizer has an optimizer, and it is
  synchronized with the cached model.

A `CachingOptimizer` has two modes of operation:

* `AUTOMATIC`: The `CachingOptimizer` changes its state when necessary. For
  example, [`optimize!`](@ref) will automatically call `attach_optimizer` (an
  optimizer must have been previously set). Attempting to add a constraint or
  perform a modification not supported by the optimizer results in a drop to
  `EMPTY_OPTIMIZER` mode.
* `MANUAL`: The user must change the state of the `CachingOptimizer` using
  [`MOIU.reset_optimizer(::JuMP.Model)`](@ref),
  [`MOIU.drop_optimizer(::JuMP.Model)`](@ref), and
  [`MOIU.attach_optimizer(::JuMP.Model)`](@ref). Attempting to perform
  an operation in the incorrect state results in an error.

By default [`Model`](@ref) will create a `CachingOptimizer` in `AUTOMATIC` mode.

### LazyBridgeOptimizer

The second layer that JuMP applies automatically is a `LazyBridgeOptimizer`. A
`LazyBridgeOptimizer` is an MOI layer that attempts to transform the problem
from the formulation provided by the user into an equivalent problem supported
by the solver. This may involve adding new variables and constraints to the
optimizer. The transformations are selected from a set of known recipes called
_bridges_.

A common example of a bridge is one that splits an interval constraint like
`@constraint(model, 1 <= x + y <= 2)` into two constraints,
`@constraint(model, x + y >= 1)` and `@constraint(model, x + y <= 2)`.

Use the `add_bridges = false` keyword to remove the bridging layer:
```jldoctest
julia> model = Model(GLPK.Optimizer; add_bridges = false)
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: GLPK

julia> backend(model)
MOIU.CachingOptimizer{GLPK.Optimizer, MOIU.UniversalFallback{MOIU.Model{Float64}}}
in state EMPTY_OPTIMIZER
in mode AUTOMATIC
with model cache MOIU.UniversalFallback{MOIU.Model{Float64}}
  fallback for MOIU.Model{Float64}
with optimizer A GLPK model
```

### Unsafe backend

In some advanced use-cases, it is necessary to work with the inner optimization
model directly. To access this model, use [`unsafe_backend`](@ref):
```jldoctest models_backends
julia> backend(model)
MOIU.CachingOptimizer{MOIB.LazyBridgeOptimizer{GLPK.Optimizer}, MOIU.UniversalFallback{MOIU.Model{Float64}}}
in state EMPTY_OPTIMIZER
in mode AUTOMATIC
with model cache MOIU.UniversalFallback{MOIU.Model{Float64}}
  fallback for MOIU.Model{Float64}
with optimizer MOIB.LazyBridgeOptimizer{GLPK.Optimizer}
  with 0 variable bridges
  with 0 constraint bridges
  with 0 objective bridges
  with inner model A GLPK model

julia> unsafe_backend(model)
A GLPK model
```

!!! warning
    [`backend`](@ref) and [`unsafe_backend`](@ref) are advanced routines. Read
    their docstrings to understand the caveats of their usage, and only call
    them if you wish to access low-level solver-specific functions.

## Direct mode

Using a `CachingOptimizer` results in an additional copy of the model being
stored by JuMP in the `.model_cache` field. To avoid this overhead, create a
JuMP model using [`direct_model`](@ref):
```jldoctest direct_mode
julia> model = direct_model(GLPK.Optimizer())
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: DIRECT
Solver name: GLPK
```

!!! warning
    Solvers that do not support incremental modification do not support
    `direct_model`. An error will be thrown, telling you to use a
    `CachingOptimizer` instead.

The benefit of using [`direct_model`](@ref) is that there are no extra layers
(for example, `Cachingoptimizer` or `LazyBridgeOptimizer`) between `model` and the
provided optimizer:
```jldoctest direct_mode
julia> backend(model)
A GLPK model
```

A downside of direct mode is that there is no bridging layer. Therefore, only
constraints which are natively supported by the solver are supported. For
example, `GLPK.jl` does not implement constraints of the form `l <= a' x <= u`.
```julia direct_mode
julia> @variable(model, x[1:2]);

julia> @constraint(model, 1 <= x[1] + x[2] <= 2)
ERROR: Constraints of type MathOptInterface.ScalarAffineFunction{Float64}-in-MathOptInterface.Interval{Float64} are not supported by the solver.

If you expected the solver to support your problem, you may have an error in your formulation. Otherwise, consider using a different solver.

The list of available solvers, along with the problem types they support, is available at https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers.
[...]
```
```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP, GLPK
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# [Solutions](@id jump_solutions)

This section of the manual describes how to access a solved solution to a
problem. It uses the following model as an example:
```jldoctest solutions
model = Model(GLPK.Optimizer)
@variable(model, x >= 0)
@variable(model, y[[:a, :b]] <= 1)
@objective(model, Max, -12x - 20y[:a])
@expression(model, my_expr, 6x + 8y[:a])
@constraint(model, my_expr >= 100)
@constraint(model, c1, 7x + 12y[:a] >= 120)
optimize!(model)
print(model)

# output

Max -12 x - 20 y[a]
Subject to
 6 x + 8 y[a] ≥ 100.0
 c1 : 7 x + 12 y[a] ≥ 120.0
 x ≥ 0.0
 y[a] ≤ 1.0
 y[b] ≤ 1.0
```

## Solutions summary

[`solution_summary`](@ref) can be used for checking the summary of the optimization solutions.

```jldoctest solutions; filter=r"[0-9]+.[0-9]+"
julia> solution_summary(model)
* Solver : GLPK

* Status
  Termination status : OPTIMAL
  Primal status      : FEASIBLE_POINT
  Dual status        : FEASIBLE_POINT
  Message from the solver:
  "Solution is optimal"

* Candidate solution
  Objective value      : -205.14285714285714
  Objective bound      : Inf
  Dual objective value : -205.1428571428571

* Work counters
  Solve time (sec)   : 0.00008

julia> solution_summary(model, verbose=true)
* Solver : GLPK

* Status
  Termination status : OPTIMAL
  Primal status      : FEASIBLE_POINT
  Dual status        : FEASIBLE_POINT
  Result count       : 1
  Has duals          : true
  Message from the solver:
  "Solution is optimal"

* Candidate solution
  Objective value      : -205.14285714285714
  Objective bound      : Inf
  Dual objective value : -205.1428571428571
  Primal solution :
    x : 15.428571428571429
    y[a] : 1.0
    y[b] : 1.0
  Dual solution :
    c1 : 1.7142857142857142

* Work counters
  Solve time (sec)   : 0.00008
```

## Why did the solver stop?

Use[`termination_status`](@ref) to understand why the solver stopped.

```jldoctest solutions
julia> termination_status(model)
OPTIMAL::TerminationStatusCode = 1
```

The [`MOI.TerminationStatusCode`](@ref) enum describes the full list of statuses
that could be returned.

Common return values include `OPTIMAL`, `LOCALLY_SOLVED`, `INFEASIBLE`,
`DUAL_INFEASIBLE`, and `TIME_LIMIT`.

!!! info
    A return status of `OPTIMAL` means the solver found (and proved) a
    globally optimal solution. A return status of `LOCALLY_SOLVED` means the
    solver found a locally optimal solution (which may also be globally
    optimal, but it could not prove so).

!!! warning
    A return status of `DUAL_INFEASIBLE` does not guarantee that the primal
    is unbounded. When the dual is infeasible, the primal is unbounded if
    there exists a feasible primal solution.

Use [`raw_status`](@ref) to get a solver-specific string explaining why the
optimization stopped:
```jldoctest solutions
julia> raw_status(model)
"Solution is optimal"
```

## Primal solutions

### Primal solution status

Use [`primal_status`](@ref) to return an [`MOI.ResultStatusCode`](@ref) enum
describing the status of the primal solution.
```jldoctest solutions
julia> primal_status(model)
FEASIBLE_POINT::ResultStatusCode = 1
```
Other common returns are `NO_SOLUTION`, and `INFEASIBILITY_CERTIFICATE`.
The first means that the solver doesn't have a solution to return, and the
second means that the primal solution is a certificate of dual infeasibility (a
primal unbounded ray).

You can also use [`has_values`](@ref), which returns `true` if there is a
solution that can be queried, and `false` otherwise.
```jldoctest solutions
julia> has_values(model)
true
```

### Objective values

The objective value of a solved problem can be obtained via
[`objective_value`](@ref). The best known bound on the optimal objective
value can be obtained via [`objective_bound`](@ref). If the solver supports it,
the value of the dual objective can be obtained via
[`dual_objective_value`](@ref).

```jldoctest solutions
julia> objective_value(model)
-205.14285714285714

julia> objective_bound(model)  # GLPK only implements objective bound for MIPs
Inf

julia> dual_objective_value(model)
-205.1428571428571
```

### Primal solution values

If the solver has a primal solution to return, use [`value`](@ref) to access it:
```julia solutions
julia> value(x)
15.428571428571429
```

Broadcast [`value`](@ref) over containers:
```julia solutions
julia> value.(y)
1-dimensional DenseAxisArray{Float64,1,...} with index sets:
    Dimension 1, Symbol[:a, :b]
And data, a 2-element Array{Float64,1}:
 1.0
 1.0
```

[`value`](@ref) also works on expressions:
```jldoctest solutions
julia> value(my_expr)
100.57142857142857
```
and constraints:
```jldoctest solutions
julia> value(c1)
120.0
```
!!! info
    Calling [`value`](@ref) on a constraint returns the constraint function
    evaluated at the solution.

## Dual solutions

### Dual solution status

Use [`dual_status`](@ref) to return an [`MOI.ResultStatusCode`](@ref) enum
describing the status of the dual solution.
```jldoctest solutions
julia> dual_status(model)
FEASIBLE_POINT::ResultStatusCode = 1
```
Other common returns are `NO_SOLUTION`, and `INFEASIBILITY_CERTIFICATE`.
The first means that the solver doesn't have a solution to return, and the
second means that the dual solution is a certificate of primal infeasibility (a
dual unbounded ray).

You can also use [`has_duals`](@ref), which returns `true` if there is a
solution that can be queried, and `false` otherwise.
```jldoctest solutions
julia> has_duals(model)
true
```

### Dual solution values

If the solver has a dual solution to return, use [`dual`](@ref) to access it:
```julia solutions
julia> dual(c1)
1.7142857142857142
```

Query the duals of variable bounds using [`LowerBoundRef`](@ref),
[`UpperBoundRef`](@ref), and [`FixRef`](@ref):
```julia solutions
julia> dual(LowerBoundRef(x))
0.0

julia> dual.(UpperBoundRef.(y))
1-dimensional DenseAxisArray{Float64,1,...} with index sets:
    Dimension 1, Symbol[:a, :b]
And data, a 2-element Array{Float64,1}:
 -0.5714285714285694
  0.0
```

!!! warning
    JuMP's definition of duality is independent of the objective sense. That is,
    the sign of feasible duals associated with a constraint depends on the
    direction of the constraint and not whether the problem is maximization or
    minimization. **This is a different convention from linear programming
    duality in some common textbooks.** If you have a linear program, and you
    want the textbook definition, you probably want to use [`shadow_price`](@ref)
    and [`reduced_cost`](@ref) instead.

```julia solutions
julia> shadow_price(c1)
1.7142857142857142

julia> reduced_cost(x)
0.0

julia> reduced_cost.(y)
1-dimensional DenseAxisArray{Float64,1,...} with index sets:
    Dimension 1, Symbol[:a, :b]
And data, a 2-element Array{Float64,1}:
  0.5714285714285694
 -0.0
```

## Recommended workflow

The recommended workflow for solving a model and querying the solution is
something like the following:
```jldoctest solutions
if termination_status(model) == OPTIMAL
    println("Solution is optimal")
elseif termination_status(model) == TIME_LIMIT && has_values(model)
    println("Solution is suboptimal due to a time limit, but a primal solution is available")
else
    error("The model was not solved correctly.")
end
println("  objective value = ", objective_value(model))
if primal_status(model) == FEASIBLE_POINT
    println("  primal solution: x = ", value(x))
end
if dual_status(model) == FEASIBLE_POINT
    println("  dual solution: c1 = ", dual(c1))
end

# output

Solution is optimal
  objective value = -205.14285714285714
  primal solution: x = 15.428571428571429
  dual solution: c1 = 1.7142857142857142
```

## OptimizeNotCalled errors

Modifying a model after calling [`optimize!`](@ref) will reset the model into
the `MOI.OPTIMIZE_NOT_CALLED` state. If you attempt to query solution
information, an `OptimizeNotCalled` error will be thrown.

If you are iteratively querying solution information and modifying a model,
query all the results first, then modify the problem.

For example, instead of:
```julia
model = Model(GLPK.Optimizer)
@variable(model, x >= 0)
optimize!(model)
set_lower_bound(x, 1)  # This will modify the model
x_val = value(x)       # This will fail because the model has been modified
set_start_value(x, x_val)
```
do
```julia
model = Model(GLPK.Optimizer)
@variable(model, x >= 0)
optimize!(model)
x_val = value(x)
set_lower_bound(x, 1)
set_start_value(x, x_val)
```

```@meta
# TODO: How to accurately measure the solve time.
```

## Accessing attributes

[MathOptInterface](@ref moi_documentation) defines many model attributes that
can be queried. Some attributes can be directly accessed by getter functions.
These include:
- [`solve_time`](@ref)
- [`relative_gap`](@ref)
- [`simplex_iterations`](@ref)
- [`barrier_iterations`](@ref)
- [`node_count`](@ref)

## Sensitivity analysis for LP

Given an LP problem and an optimal solution corresponding to a basis, we can
question how much an objective coefficient or standard form right-hand side
coefficient (c.f., [`normalized_rhs`](@ref)) can change without violating primal
or dual feasibility of the basic solution.

Note that not all solvers compute the basis, and for sensitivity analysis, the
solver interface must implement `MOI.ConstraintBasisStatus`.

To give a simple example, we could analyze the sensitivity of the optimal
solution to the following (non-degenerate) LP problem:

```jldoctest solutions_sensitivity
model = Model(GLPK.Optimizer)
@variable(model, x[1:2])
set_lower_bound(x[2], -0.5)
set_upper_bound(x[2], 0.5)
@constraint(model, c1, x[1] + x[2] <= 1)
@constraint(model, c2, x[1] - x[2] <= 1)
@objective(model, Max, x[1])
print(model)

# output

Max x[1]
Subject to
 c1 : x[1] + x[2] ≤ 1.0
 c2 : x[1] - x[2] ≤ 1.0
 x[2] ≥ -0.5
 x[2] ≤ 0.5
```

To analyze the sensitivity of the problem we could check the allowed
perturbation ranges of, for example, the cost coefficients and the right-hand side
coefficient of the constraint `c1` as follows:

```jldoctest solutions_sensitivity
julia> optimize!(model)

julia> value.(x)
2-element Vector{Float64}:
 1.0
 0.0

julia> report = lp_sensitivity_report(model);

julia> x1_lo, x1_hi = report[x[1]]
(-1.0, Inf)

julia> println("The objective coefficient of x[1] could decrease by $(x1_lo) or increase by $(x1_hi).")
The objective coefficient of x[1] could decrease by -1.0 or increase by Inf.

julia> x2_lo, x2_hi = report[x[2]]
(-1.0, 1.0)

julia> println("The objective coefficient of x[2] could decrease by $(x2_lo) or increase by $(x2_hi).")
The objective coefficient of x[2] could decrease by -1.0 or increase by 1.0.

julia> c_lo, c_hi = report[c1]
(-1.0, 1.0)

julia> println("The RHS of c1 could decrease by $(c_lo) or increase by $(c_hi).")
The RHS of c1 could decrease by -1.0 or increase by 1.0.
```

The range associated with a variable is the range of the allowed perturbation of
the corresponding objective coefficient. Note that the current primal solution
remains optimal within this range; however the corresponding dual solution might
change since a cost coefficient is perturbed. Similarly, the range associated
with a constraint is the range of the allowed perturbation of the corresponding
right-hand side coefficient. In this range the current dual solution remains
optimal, but the optimal primal solution might change.

If the problem is degenerate, there are multiple optimal bases and hence these
ranges might not be as intuitive and seem too narrow, for example, a larger cost
coefficient perturbation might not invalidate the optimality of the current
primal solution. Moreover, if a problem is degenerate, due to finite precision,
it can happen that, for example, a perturbation seems to invalidate a basis even though
it doesn't (again providing too narrow ranges). To prevent this, increase the
`atol` keyword argument to [`lp_sensitivity_report`](@ref). Note that this might
make the ranges too wide for numerically challenging instances. Thus, do not
blindly trust these ranges, especially not for highly degenerate or numerically
unstable instances.

## Conflicts

When the model you input is infeasible, some solvers can help you find the
cause of this infeasibility by offering a conflict, that is, a subset of the
constraints that create this infeasibility. Depending on the solver,
this can also be called an IIS (irreducible inconsistent subsystem).

The function [`compute_conflict!`](@ref) is used to trigger the computation of
a conflict. Once this process is finished, the attribute
[`MOI.ConflictStatus`](@ref) returns a [`MOI.ConflictStatusCode`](@ref).

If there is a conflict, you can query from each constraint whether it
participates in the conflict or not using the attribute
[`MOI.ConstraintConflictStatus`](@ref), which returns a
[`MOI.ConflictParticipationStatusCode`](@ref).

To create a new model containing only the constraints that participate in the
conflict, use [`copy_conflict`](@ref). It may be helpful to write this model
to a file for easier debugging using [`write_to_file`](@ref).

For instance, this is how you can use this functionality:

```julia
using JuMP
model = Model() # You must use a solver that supports conflict refining/IIS
# computation, like CPLEX or Gurobi
# for example, using Gurobi; model = Model(Gurobi.Optimizer)
@variable(model, x >= 0)
@constraint(model, c1, x >= 2)
@constraint(model, c2, x <= 1)
optimize!(model)

# termination_status(model) will likely be INFEASIBLE,
# depending on the solver

compute_conflict!(model)
if MOI.get(model, MOI.ConflictStatus()) != MOI.CONFLICT_FOUND
    error("No conflict could be found for an infeasible model.")
end

# Both constraints participate in the conflict.
MOI.get(model, MOI.ConstraintConflictStatus(), c1)
MOI.get(model, MOI.ConstraintConflictStatus(), c2)

# Get a copy of the model with only the constraints in the conflict.
new_model, reference_map = copy_conflict(model)
```

Conflicting constraints can be collected in a list and printed
as follows:

```julia
conflict_constraint_list = ConstraintRef[]
for (F, S) in list_of_constraint_types(model)
    for con in all_constraints(model, F, S)
        if MOI.get(model, MOI.ConstraintConflictStatus(), con) == MOI.IN_CONFLICT
            push!(conflict_constraint_list, con)
            println(con)
        end
    end
end
```

## Multiple solutions

Some solvers support returning multiple solutions. You can check how many
solutions are available to query using [`result_count`](@ref).

Functions for querying the solutions, for example, [`primal_status`](@ref) and
[`value`](@ref), all take an additional keyword argument `result` which can be
used to specify which result to return.

!!! warning
    Even if [`termination_status`](@ref) is `OPTIMAL`, some of the returned
    solutions may be suboptimal! However, if the solver found at least one
    optimal solution, then `result = 1` will always return an optimal solution.
    Use [`objective_value`](@ref) to assess the quality of the remaining
    solutions.

```julia
using JuMP
model = Model()
@variable(model, x[1:10] >= 0)
# ... other constraints ...
optimize!(model)

if termination_status(model) != OPTIMAL
    error("The model was not solved correctly.")
end

an_optimal_solution = value.(x; result = 1)
optimal_objective = objective_value(model; result = 1)
for i in 2:result_count(model)
    @assert has_values(model; result = i)
    println("Solution $(i) = ", value.(x; result = i))
    obj = objective_value(model; result = i)
    println("Objective $(i) = ", obj)
    if isapprox(obj, optimal_objective; atol = 1e-8)
        print("Solution $(i) is also optimal!")
    end
end
```

## Checking feasibility of solutions

To check the feasibility of a primal solution, use
[`primal_feasibility_report`](@ref), which takes a `model`, a dictionary mapping
each variable to a primal solution value (defaults to the last solved solution),
and a tolerance `atol` (defaults to `0.0`).

The function returns a dictionary which maps the infeasible constraint
references to the distance between the primal value of the constraint and the
nearest point in the corresponding set. A point is classed as infeasible if the
distance is greater than the supplied tolerance `atol`.

```@meta
# Add a filter here because the output of the dictionary is not ordered, and
# changes in printing order will cause the doctest to fail.
```
```jldoctest feasibility; filter=[r"x.+?\=\> 0.1", r"c1.+? \=\> 0.01"]
julia> model = Model(GLPK.Optimizer);

julia> @variable(model, x >= 1, Int);

julia> @variable(model, y);

julia> @constraint(model, c1, x + y <= 1.95);

julia> point = Dict(x => 1.9, y => 0.06);

julia> primal_feasibility_report(model, point)
Dict{Any, Float64} with 2 entries:
  x integer         => 0.1
  c1 : x + y ≤ 1.95 => 0.01

julia> primal_feasibility_report(model, point; atol = 0.02)
Dict{Any, Float64} with 1 entry:
  x integer => 0.1
```

If the point is feasible, an empty dictionary is returned:
```jldoctest feasibility
julia> primal_feasibility_report(model, Dict(x => 1.0, y => 0.0))
Dict{Any, Float64}()
```

To use the primal solution from a solve, omit the `point` argument:
```jldoctest feasibility
julia> optimize!(model)

julia> primal_feasibility_report(model)
Dict{Any, Float64}()
```

Pass `skip_mising = true` to skip constraints which contain variables that are
not in `point`:
```jldoctest feasibility
julia> primal_feasibility_report(model, Dict(x => 2.1); skip_missing = true)
Dict{Any, Float64} with 1 entry:
  x integer => 0.1
```

You can also use the functional form, where the first argument is a function
that maps variables to their primal values:
```jldoctest feasibility
julia> optimize!(model)

julia> primal_feasibility_report(v -> value(v), model)
Dict{Any, Float64}()
```
```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
```

# [Variables](@id jump_variables)

The term *variable* in mathematical optimization has many meanings. For example,
*optimization* variables (also called decision variables) are the unknowns ``x``
that we are solving for in the problem:
```math
\begin{align}
    & \min_{x \in \mathbb{R}^n} & f_0(x)
    \\
    & \;\;\text{s.t.} & f_i(x) & \in \mathcal{S}_i & i = 1 \ldots m
\end{align}
```

To complicate things, Julia uses *variable* to mean a binding between a name and
a value. For example, in the statement:
```jldoctest
julia> x = 1
1
```
`x` is a variable that stores the value `1`.

JuMP uses *variable* in a third way, to mean an instance of the
[`VariableRef`](@ref) struct. JuMP variables are the link between Julia and the
optimization variables inside a JuMP model.

This page explains how to create and manage JuMP variables in a variety of
contexts.

## Create a variable

Create variables using the [`@variable`](@ref) macro. When creating a variable,
you can also specify variable bounds:
```jldoctest variables_2
model = Model()
@variable(model, x_free)
@variable(model, x_lower >= 0)
@variable(model, x_upper <= 1)
@variable(model, 2 <= x_interval <= 3)
@variable(model, x_fixed == 4)
print(model)

# output

Feasibility
Subject to
 x_fixed = 4.0
 x_lower ≥ 0.0
 x_interval ≥ 2.0
 x_upper ≤ 1.0
 x_interval ≤ 3.0
```

!!! warning
    When creating a variable with a single lower- or upper-bound, and the
    value of the bound is not a numeric literal (for example, `1` or `1.0`), the
    name of the variable _must_ appear on the left-hand side. Putting the name
    on the right-hand side is an error. For example, to create a variable `x`:
    ```julia
    a = 1
    @variable(model, x >= 1)      # ✓ Okay
    @variable(model, 1.0 <= x)    # ✓ Okay
    @variable(model, x >= a)      # ✓ Okay
    @variable(model, a <= x)      # × Not okay
    @variable(model, x >= 1 / 2)  # ✓ Okay
    @variable(model, 1 / 2 <= x)  # × Not okay
    ```

### Containers of variables

The [`@variable`](@ref) macro also supports creating collections of JuMP
variables. We'll cover some brief syntax here; read the [Variable containers](@ref)
section for more details.

You can create arrays of JuMP variables:
```jldoctest; setup=:(model = Model())
julia> @variable(model, x[1:2, 1:2])
2×2 Matrix{VariableRef}:
 x[1,1]  x[1,2]
 x[2,1]  x[2,2]

julia> x[1, 2]
x[1,2]
```

Index sets can be named, and bounds can depend on those names:
```jldoctest; setup=:(model = Model())
julia> @variable(model, sqrt(i) <= x[i = 1:3] <= i^2)
3-element Vector{VariableRef}:
 x[1]
 x[2]
 x[3]

julia> x[2]
x[2]
```

Sets can be any Julia type that supports iteration:
```jldoctest; setup=:(model = Model())
julia> @variable(model, x[i = 2:3, j = 1:2:3, ["red", "blue"]] >= 0)
3-dimensional DenseAxisArray{VariableRef,3,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, 1:2:3
    Dimension 3, ["red", "blue"]
And data, a 2×2×2 Array{VariableRef, 3}:
[:, :, "red"] =
 x[2,1,red]  x[2,3,red]
 x[3,1,red]  x[3,3,red]

[:, :, "blue"] =
 x[2,1,blue]  x[2,3,blue]
 x[3,1,blue]  x[3,3,blue]

julia> x[2, 1, "red"]
x[2,1,red]
```

Sets can depend upon previous indices:
```jldoctest; setup=:(model = Model())
julia> @variable(model, u[i = 1:2, j = i:3])
JuMP.Containers.SparseAxisArray{VariableRef, 2, Tuple{Int64, Int64}} with 5 entries:
  [1, 1]  =  u[1,1]
  [1, 2]  =  u[1,2]
  [1, 3]  =  u[1,3]
  [2, 2]  =  u[2,2]
  [2, 3]  =  u[2,3]
```
and we can filter elements in the sets using the `;` syntax:
```jldoctest; setup=:(model = Model())
julia> @variable(model, v[i = 1:9; mod(i, 3) == 0])
JuMP.Containers.SparseAxisArray{VariableRef, 1, Tuple{Int64}} with 3 entries:
  [3]  =  v[3]
  [6]  =  v[6]
  [9]  =  v[9]
```

## Registered variables

When you create variables, JuMP registers them inside the model using their
corresponding symbol. Get a registered name using `model[:key]`:
```jldoctest
julia> model = Model()
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> @variable(model, x)
x

julia> model
A JuMP Model
Feasibility problem with:
Variable: 1
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
Names registered in the model: x

julia> model[:x] === x
true
```

Registered names are most useful when you start to write larger models and
want to break up the model construction into functions:
```jldoctest
julia> function set_objective(model::Model)
           @objective(model, Min, 2 * model[:my_x] + 1)
           return
       end
set_objective (generic function with 1 method)

julia> model = Model();

julia> @variable(model, my_x);

julia> set_objective(model)

julia> print(model)
Min 2 my_x + 1
Subject to
```

## [Anonymous variables](@id anonymous_variables)

To reduce the likelihood of accidental bugs, and because JuMP registers
variables inside a model, creating two variables with the same name is an error:
```julia
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, x)
ERROR: An object of name x is already attached to this model. If this
    is intended, consider using the anonymous construction syntax, e.g.,
    `x = @variable(model, [1:N], ...)` where the name of the object does
    not appear inside the macro.

    Alternatively, use `unregister(model, :x)` to first unregister
    the existing name from the model. Note that this will not delete the
    object; it will just remove the reference at `model[:x]`.
[...]
```

A common reason for encountering this error is adding variables in a loop.

As a work-around, JuMP provides *anonymous* variables. Create a scalar valued
anonymous variable by omitting the name argument:
```jldoctest; setup=:(model=Model())
julia> x = @variable(model)
_[1]
```
Anonymous variables get printed as an underscore followed by a unique index of
the variable.

!!! warning
    The index of the variable may not correspond to the column of the variable
    in the solver!

Create a container of anonymous JuMP variables by dropping the name in front of
the `[`:
```jldoctest; setup=:(model=Model())
julia> y = @variable(model, [1:2])
2-element Vector{VariableRef}:
 _[1]
 _[2]
```

The `<=` and `>=` short-hand cannot be used to set bounds on scalar-valued
anonymous JuMP variables. Instead, use the `lower_bound` and `upper_bound`
keywords:
```jldoctest; setup=:(model=Model())
julia> x_lower = @variable(model, lower_bound = 1.0)
_[1]

julia> x_upper = @variable(model, upper_bound = 2.0)
_[2]

julia> x_interval = @variable(model, lower_bound = 3.0, upper_bound = 4.0)
_[3]
```

## Variable names

In addition to the symbol that variables are registered with, JuMP variables
have a `String` name that is used for printing and writing to file formats.

Get and set the name of a variable using [`name`](@ref) and [`set_name`](@ref):
```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> name(x)
"x"

julia> set_name(x, "my_x_name")

julia> x
my_x_name
```

Override the default choice of name using the `base_name` keyword:
```jldoctest variable_name_vector
julia> model = Model();

julia> @variable(model, x[i=1:2], base_name = "my_var")
2-element Vector{VariableRef}:
 my_var[1]
 my_var[2]
```

Note that names apply to each element of the container, not to the container of
variables:
```jldoctest variable_name_vector
julia> name(x[1])
"my_var[1]"

julia> set_name(x[1], "my_x")

julia> x
2-element Vector{VariableRef}:
 my_x
 my_var[2]
```

### Retrieve a variable by name

Retrieve a variable from a model using [`variable_by_name`](@ref):
```jldoctest variable_name_vector
julia> variable_by_name(model, "my_x")
my_x
```
If the name is not present, `nothing` will be returned:
```jldoctest variable_name_vector
julia> variable_by_name(model, "bad_name")
```

You can only look up individual variables using [`variable_by_name`](@ref).
Something like this will not work:
```jldoctest
julia> model = Model();

julia> @variable(model, [i = 1:2], base_name = "my_var")
2-element Vector{VariableRef}:
 my_var[1]
 my_var[2]

julia> variable_by_name(model, "my_var")
```

To look up a collection of variables, do not use [`variable_by_name`](@ref).
Instead, register them using the `model[:key] = value` syntax:
```jldoctest
julia> model = Model();

julia> model[:x] = @variable(model, [i = 1:2], base_name = "my_var")
2-element Vector{VariableRef}:
 my_var[1]
 my_var[2]

julia> model[:x]
2-element Vector{VariableRef}:
 my_var[1]
 my_var[2]
```

## String names, symbolic names, and bindings

It's common for new users to experience confusion relating to JuMP variables.
Part of the problem is the overloaded use of "variable" in mathematical
optimization, along with the difference between the name that a variable is
registered under and the `String` name used for printing.

Here's a summary of the differences:

 * JuMP variables are created using [`@variable`](@ref).
 * JuMP variables can be named or anonymous.
 * Named JuMP variables have the form `@variable(model, x)`. For named
   variables:
   * The `String` name of the variable is set to `"x"`.
   * A Julia variable `x` is created that binds `x` to  the JuMP variable.
   * The name `:x` is registered as a key in the model with the value `x`.
 * Anonymous JuMP variables have the form `x = @variable(model)`. For anonymous
   variables:
   * The `String` name of the variable is set to `""`. When printed, this is
     replaced with `"_[i]"` where `i` is the index of the variable.
   * You control the name of the Julia variable used as the binding.
   * No name is registered as a key in the model.
 * The `base_name` keyword can override the `String` name of the variable.
 * You can manually register names in the model via `model[:key] = value`

Here's an example that should make things clearer:
```jldoctest
julia> model = Model();

julia> x_binding = @variable(model, base_name = "x")
x

julia> model
A JuMP Model
Feasibility problem with:
Variable: 1
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> x
ERROR: UndefVarError: x not defined

julia> x_binding
x

julia> name(x_binding)
"x"

julia> model[:x_register] = x_binding
x

julia> model
A JuMP Model
Feasibility problem with:
Variable: 1
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
Names registered in the model: x_register

julia> model[:x_register]
x

julia> model[:x_register] === x_binding
true

julia> x
ERROR: UndefVarError: x not defined
```

## Create, delete, and modify variable bounds

Query whether a variable has a bound using [`has_lower_bound`](@ref),
[`has_upper_bound`](@ref), and [`is_fixed`](@ref):
```jldoctest variables_2
julia> has_lower_bound(x_free)
false

julia> has_upper_bound(x_upper)
true

julia> is_fixed(x_fixed)
true
```

If a variable has a particular bound, query the value of it using
[`lower_bound`](@ref), [`upper_bound`](@ref), and [`fix_value`](@ref):
```jldoctest variables_2
julia> lower_bound(x_interval)
2.0

julia> upper_bound(x_interval)
3.0

julia> fix_value(x_fixed)
4.0
```

Querying the value of a bound that does not exist will result in an error.

Delete variable bounds using [`delete_lower_bound`](@ref),
[`delete_upper_bound`](@ref), and [`unfix`](@ref):
```jldoctest variables_2
julia> delete_lower_bound(x_lower)

julia> has_lower_bound(x_lower)
false

julia> delete_upper_bound(x_upper)

julia> has_upper_bound(x_upper)
false

julia> unfix(x_fixed)

julia> is_fixed(x_fixed)
false
```

Set or update variable bounds using [`set_lower_bound`](@ref),
[`set_upper_bound`](@ref), and [`fix`](@ref):
```jldoctest variables_2
julia> set_lower_bound(x_lower, 1.1)

julia> set_upper_bound(x_upper, 2.1)

julia> fix(x_fixed, 4.1)
```

Fixing a variable with existing bounds will throw an error. To delete the bounds
prior to fixing, use `fix(variable, value; force = true)`.

```jldoctest
julia> model = Model();

julia> @variable(model, x >= 1)
x

julia> fix(x, 2)
ERROR: Unable to fix x to 2 because it has existing variable bounds. Consider calling `JuMP.fix(variable, value; force=true)` which will delete existing bounds before fixing the variable.

julia> fix(x, 2; force = true)


julia> fix_value(x)
2.0
```

!!! tip
    Use [`fix`](@ref) instead of `@constraint(model, x == 2)`. The former
    modifies variable bounds, while the latter adds a new linear constraint to
    the problem.

## Binary variables

Binary variables are constrained to the set ``x \in \{0, 1\}``.

Create a binary variable by passing `Bin` as an optional positional argument:
```jldoctest variables_binary; setup=:(model=Model())
julia> @variable(model, x, Bin)
x
```

Check if a variable is binary using [`is_binary`](@ref):
```jldoctest variables_binary
julia> is_binary(x)
true
```
Delete a binary constraint using [`unset_binary`](@ref):
```jldoctest variables_binary
julia> unset_binary(x)

julia> is_binary(x)
false
```

Binary variables can also be created by setting the `binary` keyword to `true`:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x, binary=true)
x
```
or by using [`set_binary`](@ref):
```jldoctest; setup=:(model=Model())
julia> @variable(model, x)
x

julia> set_binary(x)
```

## Integer variables

Integer variables are constrained to the set ``x \in \mathbb{Z}``.

Create an integer variable by passing `Int` as an optional positional argument:
```jldoctest variables_integer; setup=:(model=Model())
julia> @variable(model, x, Int)
x
```

Check if a variable is integer using [`is_integer`](@ref):
```jldoctest variables_integer
julia> is_integer(x)
true
```
Delete an integer constraint using [`unset_integer`](@ref).
```jldoctest variables_integer
julia> unset_integer(x)

julia> is_integer(x)
false
```

Integer variables can also be created by setting the `integer` keyword to
`true`:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x, integer=true)
x
```
or by using [`set_integer`](@ref):
```jldoctest; setup=:(model=Model())
julia> @variable(model, x)
x

julia> set_integer(x)
```

!!! tip
    The [`relax_integrality`](@ref) function relaxes all integrality constraints
    in the model, returning a function that can be called to undo the operation
    later on.

## Start values

There are two ways to provide a primal starting solution (also called MIP-start
or a warmstart) for each variable:

 - using the `start` keyword in the [`@variable`](@ref) macro
 - using [`set_start_value`](@ref)

The starting value of a variable can be queried using [`start_value`](@ref). If
no start value has been set, [`start_value`](@ref) will return `nothing`.

```jldoctest variables_start; setup=:(model=Model())
julia> @variable(model, x)
x

julia> start_value(x)

julia> @variable(model, y, start = 1)
y

julia> start_value(y)
1.0

julia> set_start_value(y, 2)

julia> start_value(y)
2.0
```

!!! warning
    Some solvers do not support start values. If a solver does not support start
    values, an `MathOptInterface.UnsupportedAttribute{MathOptInterface.VariablePrimalStart}`
    error will be thrown.

!!! tip
    To set the optimal solution from a previous solve as a new starting value,
    use [`all_variables`](@ref) to get a vector of all the variables in the
    model, then run:
    ```julia
    x = all_variables(model)
    x_solution = value.(x)
    set_start_value.(x, x_solution)
    ```

## [Delete a variable](@id delete_a_variable)

Use [`delete`](@ref) to delete a variable from a model. Use [`is_valid`](@ref)
to check if a variable belongs to a model and has not been deleted.

```jldoctest variables_delete; setup=:(model=Model())
julia> @variable(model, x)
x

julia> is_valid(model, x)
true

julia> delete(model, x)

julia> is_valid(model, x)
false
```

Deleting a variable does not unregister the corresponding name from the model.
Therefore, creating a new variable of the same name will throw an error:
```jldoctest variables_delete
julia> @variable(model, x)
ERROR: An object of name x is already attached to this model. If this
    is intended, consider using the anonymous construction syntax, e.g.,
    `x = @variable(model, [1:N], ...)` where the name of the object does
    not appear inside the macro.

    Alternatively, use `unregister(model, :x)` to first unregister
    the existing name from the model. Note that this will not delete the
    object; it will just remove the reference at `model[:x]`.
[...]
```

After calling [`delete`](@ref), call [`unregister`](@ref) to remove the symbolic
reference:
```jldoctest variables_delete
julia> unregister(model, :x)

julia> @variable(model, x)
x
```

!!! info
    [`delete`](@ref) does not automatically [`unregister`](@ref) because we do
    not distinguish between names that are automatically registered by JuMP
    macros and names that are manually registered by the user by setting values
    in [`object_dictionary`](@ref). In addition, deleting a variable and then
    adding a new variable of the same name is an easy way to introduce bugs into
    your code.

## Variable containers

JuMP provides a mechanism for creating collections of variables in three types
of data structures, which we refer to as *containers*.

The three types are `Array`s, `DenseAxisArray`s, and `SparseAxisArray`s. We
explain each of these in the following.

!!! tip
    You can read more about containers in the [Containers](@ref) section.

### Arrays

We have already seen the creation of an array of JuMP variables with the
`x[1:2]` syntax. This can be extended to create multi-dimensional
arrays of JuMP variables. For example:
```jldoctest variables_arrays; setup=:(model=Model())
julia> @variable(model, x[1:2, 1:2])
2×2 Matrix{VariableRef}:
 x[1,1]  x[1,2]
 x[2,1]  x[2,2]
```

Arrays of JuMP variables can be indexed and sliced as follows:
```jldoctest variables_arrays
julia> x[1, 2]
x[1,2]

julia> x[2, :]
2-element Vector{VariableRef}:
 x[2,1]
 x[2,2]
```

Variable bounds can depend upon the indices:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[i=1:2, j=1:2] >= 2i + j)
2×2 Matrix{VariableRef}:
 x[1,1]  x[1,2]
 x[2,1]  x[2,2]

julia> lower_bound.(x)
2×2 Matrix{Float64}:
 3.0  4.0
 5.0  6.0
```

JuMP will form an `Array` of JuMP variables when it can determine at compile
time that the indices are one-based integer ranges. Therefore `x[1:b]` will
create an `Array` of JuMP variables, but `x[a:b]` will not. If JuMP cannot
determine that the indices are one-based integer ranges (for example, in the case of
`x[a:b]`), JuMP will create a `DenseAxisArray` instead.

### [DenseAxisArrays](@id variable_jump_arrays)

We often want to create arrays where the indices are not one-based integer
ranges. For example, we may want to create a variable indexed by the name of a
product or a location. The syntax is the same as that above, except with an
arbitrary vector as an index as opposed to a one-based range. The biggest
difference is that instead of returning an `Array` of JuMP variables, JuMP will
return a `DenseAxisArray`. For example:
```jldoctest variables_jump_arrays; setup=:(model=Model())
julia> @variable(model, x[1:2, [:A,:B]])
2-dimensional DenseAxisArray{VariableRef,2,...} with index sets:
    Dimension 1, Base.OneTo(2)
    Dimension 2, [:A, :B]
And data, a 2×2 Matrix{VariableRef}:
 x[1,A]  x[1,B]
 x[2,A]  x[2,B]
```

DenseAxisArrays can be indexed and sliced as follows:
```jldoctest variables_jump_arrays
julia> x[1, :A]
x[1,A]

julia> x[2, :]
1-dimensional DenseAxisArray{VariableRef,1,...} with index sets:
    Dimension 1, [:A, :B]
And data, a 2-element Vector{VariableRef}:
 x[2,A]
 x[2,B]
```

Bounds can depend upon indices:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[i=2:3, j=1:2:3] >= 0.5i + j)
2-dimensional DenseAxisArray{VariableRef,2,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, 1:2:3
And data, a 2×2 Matrix{VariableRef}:
 x[2,1]  x[2,3]
 x[3,1]  x[3,3]

julia> lower_bound.(x)
2-dimensional DenseAxisArray{Float64,2,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, 1:2:3
And data, a 2×2 Matrix{Float64}:
 2.0  4.0
 2.5  4.5
```

### [SparseAxisArrays](@id variable_sparseaxisarrays)

The third container type that JuMP natively supports is `SparseAxisArray`.
These arrays are created when the indices do not form a rectangular set.
For example, this applies when indices have a dependence upon previous
indices (called *triangular indexing*). JuMP supports this as follows:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[i=1:2, j=i:2])
JuMP.Containers.SparseAxisArray{VariableRef, 2, Tuple{Int64, Int64}} with 3 entries:
  [1, 1]  =  x[1,1]
  [1, 2]  =  x[1,2]
  [2, 2]  =  x[2,2]
```

We can also conditionally create variables via a JuMP-specific syntax. This
syntax appends a comparison check that depends upon the named indices and is
separated from the indices by a semi-colon (`;`). For example:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[i=1:4; mod(i, 2)==0])
JuMP.Containers.SparseAxisArray{VariableRef, 1, Tuple{Int64}} with 2 entries:
  [2]  =  x[2]
  [4]  =  x[4]
```

#### Performance considerations

When using the semi-colon as a filter, JuMP iterates over *all* indices and
evaluates the conditional for each combination. If there are many index
dimensions and a large amount of sparsity, this can be inefficient.

For example:
```jldoctest; setup=:(model=Model()), filter=r"[0-9\.]+ seconds.+"
julia> N = 10
10

julia> S = [(1, 1, 1), (N, N, N)]
2-element Vector{Tuple{Int64, Int64, Int64}}:
 (1, 1, 1)
 (10, 10, 10)

julia> @time @variable(model, x1[i=1:N, j=1:N, k=1:N; (i, j, k) in S])
  0.203861 seconds (392.22 k allocations: 23.977 MiB, 99.10% compilation time)
JuMP.Containers.SparseAxisArray{VariableRef, 3, Tuple{Int64, Int64, Int64}} with 2 entries:
  [1, 1, 1   ]  =  x1[1,1,1]
  [10, 10, 10]  =  x1[10,10,10]

julia> @time @variable(model, x2[S])
  0.045407 seconds (65.24 k allocations: 3.771 MiB, 99.15% compilation time)
1-dimensional DenseAxisArray{VariableRef,1,...} with index sets:
    Dimension 1, [(1, 1, 1), (10, 10, 10)]
And data, a 2-element Vector{VariableRef}:
 x2[(1, 1, 1)]
 x2[(10, 10, 10)]
```

The first option is slower because it is equivalent to:
```julia
x1 = Dict()
for i in 1:N
    for j in 1:N
        for k in 1:N
            if (i, j, k) in S
                x1[i, j, k] = @variable(model)
            end
        end
    end
end
```
If performance is a concern, explicitly construct the set of indices instead of
using the filtering syntax.

### [Forcing the container type](@id variable_forcing)

When creating a container of JuMP variables, JuMP will attempt to choose the
tightest container type that can store the JuMP variables. Thus, it will prefer
to create an Array before a DenseAxisArray and a DenseAxisArray before a
SparseAxisArray. However, because this happens at compile time, JuMP does not
always make the best choice. To illustrate this, consider the following example:
```jldoctest variable_force_container; setup=:(model=Model())
julia> A = 1:2
1:2

julia> @variable(model, x[A])
1-dimensional DenseAxisArray{VariableRef,1,...} with index sets:
    Dimension 1, 1:2
And data, a 2-element Vector{VariableRef}:
 x[1]
 x[2]
```
Since the value (and type) of `A` is unknown at parsing time, JuMP is unable to
infer that `A` is a one-based integer range. Therefore, JuMP creates a
`DenseAxisArray`, even though it could store these two variables in a standard
one-dimensional `Array`.

We can share our knowledge that it is possible to store these JuMP variables as
an array by setting the `container` keyword:
```jldoctest variable_force_container
julia> @variable(model, y[A], container=Array)
2-element Vector{VariableRef}:
 y[1]
 y[2]
```
JuMP now creates a vector of JuMP variables instead of a DenseAxisArray.
Choosing an invalid container type will throw an error.

### User-defined containers

In addition to the built-in container types, you can create your own collections
of JuMP variables.

!!! tip
    This is a point that users often overlook: you are not restricted to the
    built-in container types in JuMP.

For example, the following code creates a dictionary with symmetric matrices as
the values:
```jldoctest; setup=:(model=Model())
julia> variables = Dict{Symbol,Array{VariableRef,2}}(
           key => @variable(model, [1:2, 1:2], Symmetric, base_name = "$(key)")
           for key in [:A, :B]
       )
Dict{Symbol, Matrix{VariableRef}} with 2 entries:
  :A => [A[1,1] A[1,2]; A[1,2] A[2,2]]
  :B => [B[1,1] B[1,2]; B[1,2] B[2,2]]
```

Another common scenario is a request to add variables to existing containers,
for example:
```julia
using JuMP
model = Model()
@variable(model, x[1:2] >= 0)
# Later I want to add
@variable(model, x[3:4] >= 0)
```
This is not possible with the built-in JuMP container types. However, you can
use regular Julia types instead:
```jldoctest
model = Model()
x = model[:x] = @variable(model, [1:2], lower_bound = 0, base_name = "x")
append!(x, @variable(model, [1:2], lower_bound = 0, base_name = "y"))
model[:x]

# output

4-element Vector{VariableRef}:
 x[1]
 x[2]
 y[1]
 y[2]
```

## Semidefinite variables

A square symmetric matrix ``X`` is positive semidefinite if all eigenvalues are
nonnegative.

Declare a matrix of JuMP variables to be positive semidefinite by passing `PSD`
as an optional positional argument:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[1:2, 1:2], PSD)
2×2 LinearAlgebra.Symmetric{VariableRef, Matrix{VariableRef}}:
 x[1,1]  x[1,2]
 x[1,2]  x[2,2]
```

!!! note
    `x` must be a square 2-dimensional `Array` of JuMP variables; it cannot be a
    DenseAxisArray or a SparseAxisArray.

## Symmetric variables

Declare a square matrix of JuMP variables to be symmetric by passing
`Symmetric`  as an optional positional argument:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[1:2, 1:2], Symmetric)
2×2 LinearAlgebra.Symmetric{VariableRef, Matrix{VariableRef}}:
 x[1,1]  x[1,2]
 x[1,2]  x[2,2]
```

## [The `@variables` macro](@id variables)

If you have many [`@variable`](@ref) calls, JuMP provides the macro
[`@variables`](@ref) that can improve readability:
```jldoctest; setup=:(model=Model())
julia> @variables(model, begin
           x
           y[i=1:2] >= i, (start = i, base_name = "Y_$i")
           z, Bin
       end)
(x, VariableRef[Y_1[1], Y_2[2]], z)

julia> print(model)
Feasibility
Subject to
 Y_1[1] ≥ 1.0
 Y_2[2] ≥ 2.0
 z binary
```
The [`@variables`](@ref) macro returns a tuple of the variables that were
defined.

!!! note
    Keyword arguments must be contained within parentheses.

## [Variables constrained on creation](@id jump_variables_on_creation)

All uses of the [`@variable`](@ref) macro documented so far translate into
separate calls for variable creation and the adding of any bound or integrality
constraints.

For example, `@variable(model, x >= 0, Int)`, is equivalent to:
```julia
@variable(model, x)
set_lower_bound(x, 0.0)
set_integer(x)
```
Importantly, the bound and integrality constraints are added _after_ the
variable has been created.

However, some solvers require a set specifying the variable domain to be given
when the variable is first created. We say that these variables are
_constrained on creation_.

Use `in` within [`@variable`](@ref) to access the special syntax for
constraining variables on creation.

For example, the following creates a vector of variables that belong to the
[`SecondOrderCone`](@ref):
```jldoctest constrained_variables; setup=:(model=Model())
julia> @variable(model, y[1:3] in SecondOrderCone())
3-element Vector{VariableRef}:
 y[1]
 y[2]
 y[3]
```

For contrast, the standard syntax is as follows:
```jldoctest constrained_variables
julia> @variable(model, x[1:3])
3-element Vector{VariableRef}:
 x[1]
 x[2]
 x[3]

julia> @constraint(model, x in SecondOrderCone())
[x[1], x[2], x[3]] ∈ MathOptInterface.SecondOrderCone(3)
```

An alternate syntax to `x in Set` is to use the `set` keyword of
[`@variable`](@ref). This is most useful when creating anonymous variables:
```julia
x = @variable(model, [1:3], set = SecondOrderCone())
```

!!! note
    You cannot delete the constraint associated with a variable constrained on
    creation.

### Example: positive semidefinite variables

An alternative to the syntax in [Semidefinite variables](@ref), declare a matrix
of JuMP variables to be positive semidefinite using [`PSDCone`](@ref):
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[1:2, 1:2] in PSDCone())
2×2 LinearAlgebra.Symmetric{VariableRef, Matrix{VariableRef}}:
 x[1,1]  x[1,2]
 x[1,2]  x[2,2]
```

### Example: symmetric variables

As an alternative to the syntax in [Symmetric variables](@ref), declare a matrix
of JuMP variables to be symmetric using [`SymMatrixSpace`](@ref):
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[1:2, 1:2] in SymMatrixSpace())
2×2 LinearAlgebra.Symmetric{VariableRef, Matrix{VariableRef}}:
 x[1,1]  x[1,2]
 x[1,2]  x[2,2]
```

### Example: skew-symmetric variables

Declare a matrix of JuMP variables to be skew-symmetric using
[`SkewSymmetricMatrixSpace`](@ref):
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[1:2, 1:2] in SkewSymmetricMatrixSpace())
2×2 Matrix{AffExpr}:
 0        x[1,2]
 -x[1,2]  0
```
Note that even though `x` is 2 by 2, only one decision variable is added to
`model`; the remaining elements in `x` are linear transformations of the single
variable.

### Why use variables constrained on creation?

For most users, it does not matter if you use the constrained on creation
syntax. Therefore, use whatever syntax you find most convenient.

However, if you use [`direct_model`](@ref), you may be forced to use the
constrained on creation syntax.

The technical difference between variables constrained on creation and the
standard JuMP syntax is that variables constrained on creation calls
`MOI.add_constrained_variables`, while the standard JuMP syntax calls
`MOI.add_variables` and then `MOI.add_constraint`.

Consult the implementation of solver package you are using to see if your solver
requires `MOI.add_constrained_variables`.
```@meta
DocTestSetup = quote
    using JuMP
end
```

# Containers

JuMP provides specialized containers similar to [`AxisArrays`](https://github.com/JuliaArrays/AxisArrays.jl)
that enable multi-dimensional arrays with non-integer indices.

These containers are created automatically by JuMP's macros. Each macro has the
same basic syntax:
```julia
@macroname(model, name[key1=index1, index2; optional_condition], other stuff)
```

The containers are generated by the
`name[key1=index1, index2; optional_condition]` syntax. Everything else is
specific to the particular macro.

Containers can be named, for example, `name[key=index]`, or unnamed, for example,
`[key=index]`. We call unnamed containers _anonymous_.

We call the bits inside the square brackets and before the `;` the _index sets_.
The index sets can be named, for example, `[i = 1:4]`, or they can be unnamed, for example,
`[1:4]`.

We call the bit inside the square brackets and after the `;` the _condition_.
Conditions are optional.

In addition to the standard JuMP macros like [`@variable`](@ref) and
[`@constraint`](@ref), which construct containers of variables and constraints
respectively, you can use [`Containers.@container`](@ref) to construct
containers with arbitrary elements.

We will use this macro to explain the three types of containers that are
natively supported by JuMP: `Array`,
[`Containers.DenseAxisArray`](@ref), and [`Containers.SparseAxisArray`](@ref).

## Array

An `Array` is created when the index sets are rectangular and the index sets are
of the form `1:n`.
```jldoctest containers_array
julia> Containers.@container(x[i = 1:2, j = 1:3], (i, j))
2×3 Matrix{Tuple{Int64, Int64}}:
 (1, 1)  (1, 2)  (1, 3)
 (2, 1)  (2, 2)  (2, 3)
```
The result is a normal Julia `Array`, so you can do all the usual things.

### Slicing

Arrays can be sliced
```jldoctest containers_array
julia> x[:, 1]
2-element Vector{Tuple{Int64, Int64}}:
 (1, 1)
 (2, 1)

julia> x[2, :]
3-element Vector{Tuple{Int64, Int64}}:
 (2, 1)
 (2, 2)
 (2, 3)
```

### Looping

Use `eachindex` to loop over the elements:
```jldoctest containers_array
julia> for key in eachindex(x)
           println(x[key])
       end
(1, 1)
(2, 1)
(1, 2)
(2, 2)
(1, 3)
(2, 3)
```

### Get the index sets

Use `axes` to obtain the index sets:
```jldoctest containers_array
julia> axes(x)
(Base.OneTo(2), Base.OneTo(3))
```

### Broadcasting

Broadcasting over an Array returns an Array
```jldoctest containers_array
julia> swap(x::Tuple) = (last(x), first(x))
swap (generic function with 1 method)

julia> swap.(x)
2×3 Matrix{Tuple{Int64, Int64}}:
 (1, 1)  (2, 1)  (3, 1)
 (1, 2)  (2, 2)  (3, 2)
```

## DenseAxisArray

A [`Containers.DenseAxisArray`](@ref) is created when the index sets are
rectangular, but not of the form `1:n`. The index sets can be of any type.

```jldoctest containers_dense
julia> x = Containers.@container([i = 1:2, j = [:A, :B]], (i, j))
2-dimensional DenseAxisArray{Tuple{Int64, Symbol},2,...} with index sets:
    Dimension 1, Base.OneTo(2)
    Dimension 2, [:A, :B]
And data, a 2×2 Matrix{Tuple{Int64, Symbol}}:
 (1, :A)  (1, :B)
 (2, :A)  (2, :B)
```

### Slicing

DenseAxisArrays can be sliced
```jldoctest containers_dense
julia> x[:, :A]
1-dimensional DenseAxisArray{Tuple{Int64, Symbol},1,...} with index sets:
    Dimension 1, Base.OneTo(2)
And data, a 2-element Vector{Tuple{Int64, Symbol}}:
 (1, :A)
 (2, :A)

julia> x[1, :]
1-dimensional DenseAxisArray{Tuple{Int64, Symbol},1,...} with index sets:
    Dimension 1, [:A, :B]
And data, a 2-element Vector{Tuple{Int64, Symbol}}:
 (1, :A)
 (1, :B)
```

### Looping

Use `eachindex` to loop over the elements:
```jldoctest containers_dense
julia> for key in eachindex(x)
           println(x[key])
       end
(1, :A)
(2, :A)
(1, :B)
(2, :B)
```

### Get the index sets

Use `axes` to obtain the index sets:
```jldoctest containers_dense
julia> axes(x)
(Base.OneTo(2), [:A, :B])
```

### Broadcasting

Broadcasting over a DenseAxisArray returns a DenseAxisArray
```jldoctest containers_dense
julia> swap(x::Tuple) = (last(x), first(x))
swap (generic function with 1 method)

julia> swap.(x)
2-dimensional DenseAxisArray{Tuple{Symbol, Int64},2,...} with index sets:
    Dimension 1, Base.OneTo(2)
    Dimension 2, [:A, :B]
And data, a 2×2 Matrix{Tuple{Symbol, Int64}}:
 (:A, 1)  (:B, 1)
 (:A, 2)  (:B, 2)
```

### Access internal data

Use `Array(x)` to copy the internal data array into a new `Array`:
```jldoctest containers_dense
julia> Array(x)
2×2 Matrix{Tuple{Int64, Symbol}}:
 (1, :A)  (1, :B)
 (2, :A)  (2, :B)
```

To access the internal data without a copy, use `x.data`.
```jldoctest containers_dense
julia> x.data
2×2 Matrix{Tuple{Int64, Symbol}}:
 (1, :A)  (1, :B)
 (2, :A)  (2, :B)
```

## SparseAxisArray

A [`Containers.SparseAxisArray`](@ref) is created when the index sets are
non-rectangular. This occurs in two circumstances:

An index depends on a prior index:
```jldoctest containers_sparse
julia> Containers.@container([i = 1:2, j = i:2], (i, j))
JuMP.Containers.SparseAxisArray{Tuple{Int64, Int64}, 2, Tuple{Int64, Int64}} with 3 entries:
  [1, 1]  =  (1, 1)
  [1, 2]  =  (1, 2)
  [2, 2]  =  (2, 2)
```

The `[indices; condition]` syntax is used:
```jldoctest containers_sparse
julia> x = Containers.@container([i = 1:3, j = [:A, :B]; i > 1 && j == :B], (i, j))
JuMP.Containers.SparseAxisArray{Tuple{Int64, Symbol}, 2, Tuple{Int64, Symbol}} with 2 entries:
  [2, B]  =  (2, :B)
  [3, B]  =  (3, :B)
```
Here we have the index sets `i = 1:3, j = [:A, :B]`, followed by `;`, and then a
condition, which evaluates to `true` or `false`: `i > 1 && j == :B`.

### Slicing

```@meta
# TODO: This is included so we know to update the documentation when this is fixed.
```

Slicing is not supported.
```jldoctest containers_sparse
julia> x[:, :B]
ERROR: ArgumentError: Indexing with `:` is not supported by Containers.SparseAxisArray
[...]
```

### Looping

Use `eachindex` to loop over the elements:
```jldoctest containers_sparse
julia> for key in eachindex(x)
           println(x[key])
       end
(2, :B)
(3, :B)
```

### Broadcasting

Broadcasting over a SparseAxisArray returns a SparseAxisArray

```jldoctest containers_sparse
julia> swap(x::Tuple) = (last(x), first(x))
swap (generic function with 1 method)

julia> swap.(x)
JuMP.Containers.SparseAxisArray{Tuple{Symbol, Int64}, 2, Tuple{Int64, Symbol}} with 2 entries:
  [2, B]  =  (:B, 2)
  [3, B]  =  (:B, 3)
```

## Forcing the container type

Pass `container = T` to use `T` as the container. For example:
```jldoctest; filter=r"\([1-2], [1-2]\) \=\> [2-4]"
julia> Containers.@container([i = 1:2, j = 1:2], i + j, container = Array)
2×2 Matrix{Int64}:
 2  3
 3  4

julia> Containers.@container([i = 1:2, j = 1:2], i + j, container = Dict)
Dict{Tuple{Int64, Int64}, Int64} with 4 entries:
  (1, 2) => 3
  (1, 1) => 2
  (2, 2) => 4
  (2, 1) => 3
```
You can also pass `DenseAxisArray` or `SparseAxisArray`.

## How different container types are chosen

If the compiler can prove _at compile time_ that the index sets are rectangular,
and indexed by a compact set of integers that start at `1`,
[`Containers.@container`](@ref) will return an array. This is the case if your
index sets are visible to the macro as `1:n`:
```jldoctest
julia> Containers.@container([i=1:3, j=1:5], i + j)
3×5 Matrix{Int64}:
 2  3  4  5  6
 3  4  5  6  7
 4  5  6  7  8
```
or an instance of `Base.OneTo`:
```jldoctest; setup=:(model = Model())
julia> set = Base.OneTo(3)
Base.OneTo(3)

julia> Containers.@container([i=set, j=1:5], i + j)
3×5 Matrix{Int64}:
 2  3  4  5  6
 3  4  5  6  7
 4  5  6  7  8
```

If the compiler can prove that the index set is rectangular, but not necessarily
of the form `1:n` at compile time, then a [`Containers.DenseAxisArray`](@ref)
will be constructed instead:
```jldoctest
julia> set = 1:3
1:3

julia> Containers.@container([i=set, j=1:5], i + j)
2-dimensional DenseAxisArray{Int64,2,...} with index sets:
    Dimension 1, 1:3
    Dimension 2, Base.OneTo(5)
And data, a 3×5 Matrix{Int64}:
 2  3  4  5  6
 3  4  5  6  7
 4  5  6  7  8
```

!!! info
    What happened here? Although we know that `set` contains `1:3`, at compile
    time the `typeof(set)` is a `UnitRange{Int}`. Therefore, Julia can't prove
    that the range starts at `1` (it only finds this out at runtime), and it
    defaults to a  `DenseAxisArray`. The case where we explicitly wrote
    `i = 1:3` worked because the macro can "see" the `1` at compile time.

However, if you know that the indices do form an `Array`, you can force the
container type with `container = Array`:
```jldoctest
julia> set = 1:3
1:3

julia> Containers.@container([i=set, j=1:5], i + j, container = Array)
3×5 Matrix{Int64}:
 2  3  4  5  6
 3  4  5  6  7
 4  5  6  7  8
```

Here's another example with something similar:
```jldoctest
julia> a = 1
1

julia> Containers.@container([i=a:3, j=1:5], i + j)
2-dimensional DenseAxisArray{Int64,2,...} with index sets:
    Dimension 1, 1:3
    Dimension 2, Base.OneTo(5)
And data, a 3×5 Matrix{Int64}:
 2  3  4  5  6
 3  4  5  6  7
 4  5  6  7  8

julia> Containers.@container([i=1:a, j=1:5], i + j)
1×5 Matrix{Int64}:
 2  3  4  5  6
```

Finally, if the compiler cannot prove that the index set is rectangular, a
[`Containers.SparseAxisArray`](@ref) will be created.

This occurs when some indices depend on a previous one:
```jldoctest
julia> Containers.@container([i=1:3, j=1:i], i + j)
JuMP.Containers.SparseAxisArray{Int64, 2, Tuple{Int64, Int64}} with 6 entries:
  [1, 1]  =  2
  [2, 1]  =  3
  [2, 2]  =  4
  [3, 1]  =  4
  [3, 2]  =  5
  [3, 3]  =  6
```
or if there is a condition on the index sets:
```jldoctest
julia> Containers.@container([i = 1:5; isodd(i)], i^2)
JuMP.Containers.SparseAxisArray{Int64, 1, Tuple{Int64}} with 3 entries:
  [1]  =  1
  [3]  =  9
  [5]  =  25
```

The condition can depend on multiple indices, the only requirement is that it is
an expression that returns `true` or `false`:
```jldoctest
julia> condition(i, j) = isodd(i) && iseven(j)
condition (generic function with 1 method)

julia> Containers.@container([i = 1:2, j = 1:4; condition(i, j)], i + j)
JuMP.Containers.SparseAxisArray{Int64, 2, Tuple{Int64, Int64}} with 2 entries:
  [1, 2]  =  3
  [1, 4]  =  5
```
```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# Nonlinear Modeling

JuMP has support for general smooth nonlinear (convex and nonconvex)
optimization problems. JuMP is able to provide exact, sparse second-order
derivatives to solvers. This information can improve solver accuracy and
performance.

There are three main changes to solve nonlinear programs in JuMP.
 * Use [`@NLobjective`](@ref) instead of [`@objective`](@ref)
 * Use [`@NLconstraint`](@ref) instead of [`@constraint`](@ref)
 * Use [`@NLexpression`](@ref) instead of [`@expression`](@ref)

!!! info
    There are some restrictions on what syntax you can use in the `@NLxxx`
    macros. Make sure to read the [Syntax notes](@ref).

## Set a nonlinear objective

Use [`@NLobjective`](@ref) to set a nonlinear objective.

```jldoctest; setup=:(model = Model(); @variable(model, x[1:2]))
julia> @NLobjective(model, Min, exp(x[1]) - sqrt(x[2]))
```
## Add a nonlinear constraint

Use [`@NLconstraint`](@ref) to add a nonlinear constraint.

```jldoctest; setup=:(model = Model(); @variable(model, x[1:2]))
julia> @NLconstraint(model, exp(x[1]) <= 1)
exp(x[1]) - 1.0 ≤ 0

julia> @NLconstraint(model, [i = 1:2], x[i]^i >= i)
2-element Vector{NonlinearConstraintRef{ScalarShape}}:
 x[1] ^ 1.0 - 1.0 ≥ 0
 x[2] ^ 2.0 - 2.0 ≥ 0

julia> @NLconstraint(model, con[i = 1:2], prod(x[j] for j = 1:i) == i)
2-element Vector{NonlinearConstraintRef{ScalarShape}}:
 (*)(x[1]) - 1.0 = 0
 x[1] * x[2] - 2.0 = 0
```

!!! info
    You can only create nonlinear constraints with `<=`, `>=`, and `==`.
    More general `Nonlinear`-in-`Set` constraints are not supported.

## Create a nonlinear expression

Use [`@NLexpression`](@ref) to create nonlinear expression objects. The syntax
is identical to [`@expression`](@ref), except that the expression can contain
nonlinear terms.

```jldoctest nl_expression; setup=:(model = Model(); @variable(model, x[1:2]))
julia> expr = @NLexpression(model, exp(x[1]) + sqrt(x[2]))
subexpression[1]: exp(x[1]) + sqrt(x[2])

julia> my_anon_expr = @NLexpression(model, [i = 1:2], sin(x[i]))
2-element Vector{NonlinearExpression}:
 subexpression[2]: sin(x[1])
 subexpression[3]: sin(x[2])

julia> @NLexpression(model, my_expr[i = 1:2], sin(x[i]))
2-element Vector{NonlinearExpression}:
 subexpression[4]: sin(x[1])
 subexpression[5]: sin(x[2])
```

Nonlinear expression can be used in [`@NLobjective`](@ref), [`@NLconstraint`](@ref),
and even nested in other [`@NLexpression`](@ref)s.

```jldoctest nl_expression
julia> @NLobjective(model, Min, expr^2 + 1)

julia> @NLconstraint(model, [i = 1:2], my_expr[i] <= i)
2-element Vector{NonlinearConstraintRef{ScalarShape}}:
 subexpression[4] - 1.0 ≤ 0
 subexpression[5] - 2.0 ≤ 0

julia> @NLexpression(model, nested[i = 1:2], sin(my_expr[i]))
2-element Vector{NonlinearExpression}:
 subexpression[6]: sin(subexpression[4])
 subexpression[7]: sin(subexpression[5])
```

## Create a nonlinear parameter

For nonlinear models only, JuMP offers a syntax for explicit "parameter" objects,
which are constants in the model that can be efficiently updated between solves.

Nonlinear parameters are declared by using the [`@NLparameter`](@ref) macro
and may be indexed by arbitrary sets analogously to JuMP variables and
expressions.

The initial value of the parameter must be provided on the right-hand side of
the `==` sign.

```jldoctest nonlinear_parameters; setup=:(model = Model(); @variable(model, x))
julia> @NLparameter(model, p[i = 1:2] == i)
2-element Vector{NonlinearParameter}:
 parameter[1] == 1.0
 parameter[2] == 2.0
```

Create anonymous parameters using the `value` keyword:
```jldoctest nonlinear_parameters
julia> anon_parameter = @NLparameter(model, value = 1)
parameter[3] == 1.0
```

!!! info
    A parameter is not an optimization variable. It must be fixed to a value with
    `==`. If you want a parameter that is `<=` or `>=`, create a variable instead
    using [`@variable`](@ref).

Use [`value`](@ref) and [`set_value`](@ref) to query or update the value of a
parameter.

```jldoctest nonlinear_parameters
julia> value.(p)
2-element Vector{Float64}:
 1.0
 2.0

julia> set_value(p[2], 3.0)
3.0

julia> value.(p)
2-element Vector{Float64}:
 1.0
 3.0
```

Nonlinear parameters can be used *within nonlinear macros* only:

```jldoctest nonlinear_parameters
julia> @objective(model, Max, p[1] * x)
ERROR: MethodError: no method matching *(::NonlinearParameter, ::VariableRef)
[...]

julia> @NLobjective(model, Max, p[1] * x)

julia> @expression(model, my_expr, p[1] * x^2)
ERROR: MethodError: no method matching *(::NonlinearParameter, ::QuadExpr)
Closest candidates are:
[...]

julia> @NLexpression(model, my_nl_expr, p[1] * x^2)
subexpression[1]: parameter[1] * x ^ 2.0
```

### When to use a parameter

Nonlinear parameters are useful when solving nonlinear models in a sequence:

```@example
using JuMP, Ipopt
model = Model(Ipopt.Optimizer)
set_silent(model)
@variable(model, z)
@NLparameter(model, x == 1.0)
@NLobjective(model, Min, (z - x)^2)
optimize!(model)
@show value(z) # Equals 1.0.

# Now, update the value of x to solve a different problem.
set_value(x, 5.0)
optimize!(model)
@show value(z) # Equals 5.0
nothing #hide
```

!!! info
    Using nonlinear parameters can be faster than creating a new model from
    scratch with updated data because JuMP is able to avoid repeating a number
    of steps in processing the model before handing it off to the solver.

## Syntax notes

The syntax accepted in nonlinear macros is more restricted than the syntax
for linear and quadratic macros. We note some important points below.

### No operator overloading

There is no operator overloading provided to build up nonlinear expressions.
For example, if `x` is a JuMP variable, the code `3x` will return an
`AffExpr` object that can be used inside of future expressions and linear
constraints. However, the code `sin(x)` is an error. All nonlinear
expressions must be inside of macros.

```jldoctest; setup=:(model = Model(); @variable(model, x))
julia> expr = sin(x) + 1
ERROR: sin is not defined for type AbstractVariableRef. Are you trying to build a nonlinear problem? Make sure you use @NLconstraint/@NLobjective.
[...]

julia> expr = @NLexpression(model, sin(x) + 1)
subexpression[1]: sin(x) + 1.0
```

### Scalar operations only

Except for the splatting syntax discussed below, all expressions
must be simple scalar operations. You cannot use `dot`, matrix-vector products,
vector slices, etc.
```jldoctest nlp_scalar_only; setup=:(model = Model(); @variable(model, x[1:2]); @variable(model, y); c = [1, 2])
julia> @NLobjective(model, Min, c' * x + 3y)
ERROR: Unexpected array [1 2] in nonlinear expression. Nonlinear expressions may contain only scalar expressions.
[...]
```

Translate vector operations into explicit `sum()` operations:
```jldoctest nlp_scalar_only
julia> @NLobjective(model, Min, sum(c[i] * x[i] for i = 1:2) + 3y)
```

Or use an [`@expression`](@ref):
```jldoctest nlp_scalar_only
julia> @expression(model, expr, c' * x)
x[1] + 2 x[2]

julia> @NLobjective(model, Min, expr + 3y)

```

### Splatting

The [splatting operator](https://docs.julialang.org/en/v1/manual/faq/#...-splits-one-argument-into-many-different-arguments-in-function-calls-1)
  `...` is recognized in a very restricted setting for expanding function
  arguments. The expression splatted can be *only* a symbol. More complex
  expressions are not recognized.

```jldoctest; filter=r"≤|<="
julia> model = Model();

julia> @variable(model, x[1:3]);

julia> @NLconstraint(model, *(x...) <= 1.0)
x[1] * x[2] * x[3] - 1.0 ≤ 0

julia> @NLconstraint(model, *((x / 2)...) <= 0.0)
ERROR: LoadError: Unexpected expression in (*)(x / 2...). JuMP supports splatting only symbols. For example, x... is ok, but (x + 1)..., [x; y]... and g(f(y)...) are not.
```

## User-defined Functions

JuMP's library of recognized univariate functions is derived from the
[Calculus.jl](https://github.com/johnmyleswhite/Calculus.jl) package.

In addition to this list of functions, it is possible to register custom
*user-defined* nonlinear functions.

!!! tip
    User-defined functions can be used anywhere in [`@NLobjective`](@ref),
    [`@NLconstraint`](@ref), and [`@NLexpression`](@ref).

!!! tip
    JuMP will attempt to automatically register functions it detects in your
    nonlinear expressions, which usually means manually registering
    a function is not needed. Two exceptions are if you want to provide custom
    derivatives, or if the function is not available in the scope of the
    nonlinear expression.

!!! warning
    User-defined functions must return a scalar output. For a work-around, see
    [User-defined functions with vector outputs](@ref).

### Automatic differentiation

JuMP does not support black-box optimization, so all user-defined functions must
provide derivatives in some form.

Fortunately, JuMP supports **automatic differentiation of user-defined functions**,
a feature to our knowledge not available in any comparable modeling systems.

!!! info
    Automatic differentiation is *not* finite differencing. JuMP's automatically
    computed derivatives are not subject to approximation error.

JuMP uses [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) to
perform automatic differentiation; see the ForwardDiff.jl
[documentation](https://www.juliadiff.org/ForwardDiff.jl/v0.10.2/user/limitations.html)
for a description of how to write a function suitable for automatic
differentiation.

!!! warning
    Get an error like `No method matching Float64(::ForwardDiff.Dual)`? Read
    this section, and see the guidelines at [ForwardDiff.jl](https://www.juliadiff.org/ForwardDiff.jl/release-0.10/user/limitations.html).

The most common error is that your user-defined function is not generic with
respect to the number type, that is, don't assume that the input to the function
is `Float64`.
```julia
f(x::Float64) = 2 * x  # This will not work.
f(x::Real)    = 2 * x  # This is good.
f(x)          = 2 * x  # This is also good.
```

Another reason you may encounter this error is if you create arrays inside
your function which are `Float64`.
```julia
function bad_f(x...)
    y = zeros(length(x))  # This constructs an array of `Float64`!
    for i = 1:length(x)
        y[i] = x[i]^i
    end
    return sum(y)
end

function good_f(x::T...) where {T<:Real}
    y = zeros(T, length(x))  # Construct an array of type `T` instead!
    for i = 1:length(x)
        y[i] = x[i]^i
    end
    return sum(y)
end
```

### Register a function

To register a user-defined function with derivatives computed by
automatic differentiation, use the [`register`](@ref) method as in the following
example:

```@example
using JuMP #hide
square(x) = x^2
f(x, y) = (x - 1)^2 + (y - 2)^2

model = Model()

register(model, :square, 1, square; autodiff = true)
register(model, :my_f, 2, f; autodiff = true)

@variable(model, x[1:2] >= 0.5)
@NLobjective(model, Min, my_f(x[1], square(x[2])))
```

The above code creates a JuMP model with the objective function
`(x[1] - 1)^2 + (x[2]^2 - 2)^2`. The arguments to [`register`](@ref) are:
 1. The model for which the functions are registered.
 2. A Julia symbol object which serves as the name of the user-defined function
    in JuMP expressions.
 3. The number of input arguments that the function takes.
 4. The Julia method which computes the function
 5. A flag to instruct JuMP to compute exact gradients automatically.

!!! tip
    The symbol `:my_f` doesn't have to match the name of the function `f`.
    However, it's more readable if it does. Make sure you use `my_f`
    and not `f` in the macros.

!!! warning
    If you use multi-variate user-defined functions, JuMP will disable
    second-derivative information. This can lead to significant slow-downs in
    some cases. Only use a user-defined function if you cannot write out the
    expression algebraically in the macro.

!!! warning
    User-defined functions cannot be re-registered and will not update if you
    modify the underlying Julia function. If you want to change a user-defined
    function between solves, rebuild the model or use a different name. To use
    a different name programmatically, see [Raw expression input](@ref).

### Register a function and gradient

Forward-mode automatic differentiation as implemented by ForwardDiff.jl has a
computational cost that scales linearly with the number of input dimensions. As
such, it is not the most efficient way to compute gradients of user-defined
functions if the number of input arguments is large. In this case, users may
want to provide their own routines for evaluating gradients.

#### Univariate functions

For univariate functions, the gradient function `∇f` returns a number that
represents the first-order derivative:
```@example
using JuMP #hide
f(x) = x^2
∇f(x) = 2x
model = Model()
register(model, :my_square, 1, f, ∇f; autodiff = true)
@variable(model, x >= 0)
@NLobjective(model, Min, my_square(x))
```
If `autodiff = true`, JuMP will use automatic differentiation to compute the
hessian.

#### Multivariate functions

For multivariate functions, the gradient function `∇f` must take a gradient
vector as the first argument that is filled in-place:
```@example
using JuMP #hide
f(x, y) = (x - 1)^2 + (y - 2)^2
function ∇f(g::AbstractVector{T}, x::T, y::T) where {T}
    g[1] = 2 * (x - 1)
    g[2] = 2 * (y - 2)
    return
end

model = Model()
register(model, :my_square, 2, f, ∇f)
@variable(model, x[1:2] >= 0)
@NLobjective(model, Min, my_square(x[1], x[2]))
```

!!! warning
    Make sure the first argument to `∇f` supports an `AbstractVector`, and do
    not assume the input is `Float64`.

### Register a function, gradient, and hessian

!!! warning
    The ability to explicitly register a hessian is only available for
    univariate functions.

Instead of automatically differentiating the hessian, you can instead pass a
function which returns a number representing the second-order derivative.

```@example
using JuMP #hide
f(x) = x^2
∇f(x) = 2x
∇²f(x) = 2
model = Model()
register(model, :my_square, 1, f, ∇f, ∇²f)
@variable(model, x >= 0)
@NLobjective(model, Min, my_square(x))
```

### User-defined functions with vector inputs

User-defined functions which take vectors as input arguments (for example,
`f(x::Vector)`) are *not* supported. Instead, use Julia's splatting syntax to
create a function with scalar arguments. For example, instead of
```julia
f(x::Vector) = sum(x[i]^i for i in 1:length(x))
```
define:
```julia
f(x...) = sum(x[i]^i for i in 1:length(x))
```

This function `f` can be used in a JuMP model as follows:
```@example
using JuMP #hide
model = Model()
@variable(model, x[1:5] >= 0)
f(x...) = sum(x[i]^i for i in 1:length(x))
register(model, :f, 5, f; autodiff = true)
@NLobjective(model, Min, f(x...))
```

!!! tip
    Make sure to read the syntax restrictions of [Splatting](@ref).

## Factors affecting solution time

The execution time when solving a nonlinear programming problem can be divided
into two parts, the time spent in the optimization algorithm (the solver) and
the time spent evaluating the nonlinear functions and corresponding derivatives.
Ipopt explicitly displays these two timings in its output, for example:

```
Total CPU secs in IPOPT (w/o function evaluations)   =      7.412
Total CPU secs in NLP function evaluations           =      2.083
```

For Ipopt in particular, one can improve the performance by installing advanced
sparse linear algebra packages, see [Installation Guide](@ref). For other
solvers, see their respective documentation for performance tips.

The function evaluation time, on the other hand, is the responsibility of the
modeling language. JuMP computes derivatives by using reverse-mode automatic
differentiation with graph coloring methods for exploiting sparsity of the
Hessian matrix [^1]. As a conservative bound, JuMP's performance here currently
may be expected to be within a factor of 5 of AMPL's.

## Querying derivatives from a JuMP model

For some advanced use cases, one may want to directly query the derivatives of a
JuMP model instead of handing the problem off to a solver.
Internally, JuMP implements the [`MOI.AbstractNLPEvaluator`](@ref) interface. To
obtain an NLP evaluator object from a JuMP model, use [`NLPEvaluator`](@ref).
[`index`](@ref) returns the [`MOI.VariableIndex`](@ref) corresponding to a JuMP
variable. `MOI.VariableIndex` itself is a type-safe wrapper for `Int64` (stored
in the `.value` field.)

For example:

```jldoctest derivatives
raw_index(v::MOI.VariableIndex) = v.value
model = Model()
@variable(model, x)
@variable(model, y)
@NLobjective(model, Min, sin(x) + sin(y))
values = zeros(2)
x_index = raw_index(JuMP.index(x))
y_index = raw_index(JuMP.index(y))
values[x_index] = 2.0
values[y_index] = 3.0
d = NLPEvaluator(model)
MOI.initialize(d, [:Grad])
MOI.eval_objective(d, values) # == sin(2.0) + sin(3.0)

# output
1.0504174348855488
```

```jldoctest derivatives
∇f = zeros(2)
MOI.eval_objective_gradient(d, ∇f, values)
(∇f[x_index], ∇f[y_index]) # == (cos(2.0), cos(3.0))

# output
(-0.4161468365471424, -0.9899924966004454)
```

Only nonlinear constraints (those added with [`@NLconstraint`](@ref)), and
nonlinear objectives (added with [`@NLobjective`](@ref)) exist in the scope of
the [`NLPEvaluator`](@ref).

The [`NLPEvaluator`](@ref) *does not evaluate derivatives of linear or quadratic
constraints or objectives*.

The [`index`](@ref) method applied to a nonlinear constraint reference object
returns its index as a [`NonlinearConstraintIndex`](@ref). The `.value` field of
[`NonlinearConstraintIndex`](@ref) stores the raw integer index. For example:

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @NLconstraint(model, cons1, sin(x) <= 1);

julia> @NLconstraint(model, cons2, x + 5 == 10);

julia> typeof(cons1)
NonlinearConstraintRef{ScalarShape} (alias for ConstraintRef{Model, NonlinearConstraintIndex, ScalarShape})

julia> index(cons1)
NonlinearConstraintIndex(1)

julia> index(cons2)
NonlinearConstraintIndex(2)
```

```@meta
# TODO: Provide a link for how to access the linear and quadratic parts of the
# model.
```

Note that for one-sided nonlinear constraints, JuMP subtracts any values on the
right-hand side when computing expressions. In other words, one-sided nonlinear
constraints are always transformed to have a right-hand side of zero.

This method of querying derivatives directly from a JuMP model is convenient for
interacting with the model in a structured way, for example, for accessing derivatives
of specific variables. For example, in statistical maximum likelihood estimation
problems, one is often interested in the Hessian matrix at the optimal solution,
which can be queried using the [`NLPEvaluator`](@ref).

## Raw expression input

!!! warning
    This section requires advanced knowledge of Julia's `Expr`. You should read
    the [Expressions and evaluation](https://docs.julialang.org/en/v1/manual/metaprogramming/#Expressions-and-evaluation)
    section of the Julia documentation first.

In addition to the [`@NLexpression`](@ref), [`@NLobjective`](@ref) and
[`@NLconstraint`](@ref) macros, it is also possible to provide Julia `Expr`
objects directly by using [`add_NL_expression`](@ref),
[`set_NL_objective`](@ref) and [`add_NL_constraint`](@ref).

This input form may be useful if the expressions are generated programmatically.

### Add a nonlinear expression

Use [`add_NL_expression`](@ref) to add a nonlinear expression to the model.

```jldoctest; setup=:(using JuMP; model = Model())
julia> @variable(model, x)
x

julia> expr = :($(x) + sin($(x)^2))
:(x + sin(x ^ 2))

julia> expr_ref = add_NL_expression(model, expr)
subexpression[1]: x + sin(x ^ 2.0)
```
This is equivalent to
```jldoctest; setup=:(using JuMP; model = Model(); @variable(model, x))
julia> expr_ref = @NLexpression(model, x + sin(x^2))
subexpression[1]: x + sin(x ^ 2.0)
```

!!! note
    You must interpolate the variables directly into the expression `expr`.

### Set the objective function

Use [`set_NL_objective`](@ref) to set a nonlinear objective.

```jldoctest; setup=:(using JuMP; model = Model(); @variable(model, x))
julia> expr = :($(x) + $(x)^2)
:(x + x ^ 2)

julia> set_NL_objective(model, MIN_SENSE, expr)
```
This is equivalent to
```jldoctest; setup=:(using JuMP; model = Model(); @variable(model, x))
julia> @NLobjective(model, Min, x + x^2)
```

!!! note
    You must use `MIN_SENSE` or `MAX_SENSE` instead of `Min` and `Max`.

### Add a constraint

Use [`add_NL_constraint`](@ref) to add a nonlinear constraint.

```jldoctest; setup=:(using JuMP; model = Model(); @variable(model, x))
julia> expr = :($(x) + $(x)^2)
:(x + x ^ 2)

julia> add_NL_constraint(model, :($(expr) <= 1))
(x + x ^ 2.0) - 1.0 ≤ 0
```

This is equivalent to
```jldoctest; setup=:(using JuMP; model = Model(); @variable(model, x))
julia> @NLconstraint(model, Min, x + x^2 <= 1)
(x + x ^ 2.0) - 1.0 ≤ 0
```

### More complicated examples

Raw expression input is most useful when the expressions are generated
programmatically, often in conjunction with user-defined functions.

As an example, we construct a model with the nonlinear constraints `f(x) <= 1`,
where `f(x) = x^2` and `f(x) = sin(x)^2`:
```jldoctest
julia> function main(functions::Vector{Function})
           model = Model()
           @variable(model, x)
           for (i, f) in enumerate(functions)
               f_sym = Symbol("f_$(i)")
               register(model, f_sym, 1, f; autodiff = true)
               add_NL_constraint(model, :($(f_sym)($(x)) <= 1))
           end
           print(model)
           return
       end
main (generic function with 1 method)

julia> main([x -> x^2, x -> sin(x)^2])
Feasibility
Subject to
 f_1(x) - 1.0 ≤ 0
 f_2(x) - 1.0 ≤ 0
```

As another example, we construct a model with the constraint
`x^2 + sin(x)^2 <= 1`:
```jldoctest
julia> function main(functions::Vector{Function})
           model = Model()
           @variable(model, x)
           expr = Expr(:call, :+)
           for (i, f) in enumerate(functions)
               f_sym = Symbol("f_$(i)")
               register(model, f_sym, 1, f; autodiff = true)
               push!(expr.args, :($(f_sym)($(x))))
           end
           add_NL_constraint(model, :($(expr) <= 1))
           print(model)
           return
       end
main (generic function with 1 method)

julia> main([x -> x^2, x -> sin(x)^2])
Feasibility
Subject to
 (f_1(x) + f_2(x)) - 1.0 ≤ 0
```

[^1]: Dunning, Huchette, and Lubin, "JuMP: A Modeling Language for Mathematical Optimization", SIAM Review, [PDF](https://mlubin.github.io/pdf/jump-sirev.pdf).
```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# Objectives

This page describes macros and functions related to linear and quadratic
objective functions only, unless otherwise indicated. For nonlinear objective
functions, see [Nonlinear Modeling](@ref).

## Set a linear objective

Use the [`@objective`](@ref) macro to set a linear objective function.

Use `Min` to create a minimization objective:
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Min, 2x + 1)
2 x + 1
```

Use `Max` to create a maximization objective:
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Max, 2x + 1)
2 x + 1
```

## Set a quadratic objective

Use the [`@objective`](@ref) macro to set a quadratic objective function.

Use `^2` to have a variable squared:
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Min, x^2 + 2x + 1)
x² + 2 x + 1
```

You can also have bilinear terms between variables:
```jldoctest; setup = :(model=Model())
julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> @objective(model, Max, x * y + x + y)
x*y + x + y
```

## Query the objective function

Use [`objective_function`](@ref) to return the current objective function.
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Min, 2x + 1)
2 x + 1

julia> objective_function(model)
2 x + 1
```

## Evaluate the objective function at a point

Use [`value`](@ref) to evaluate an objective function at a point specifying values for variables.

```jldoctest; setup = :(model=Model())
julia> @variable(model, x[1:2]);

julia> @objective(model, Min, 2x[1]^2 + x[1] + 0.5*x[2])
2 x[1]² + x[1] + 0.5 x[2]

julia> f = objective_function(model)
2 x[1]² + x[1] + 0.5 x[2]

julia> point = Dict(x[1] => 2.0, x[2] => 1.0);

julia> value(z -> point[z], f)
10.5
```

## Query the objective sense

Use [`objective_sense`](@ref) to return the current objective sense.
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Min, 2x + 1)
2 x + 1

julia> objective_sense(model)
MIN_SENSE::OptimizationSense = 0
```

## Modify an objective

To modify an objective, call [`@objective`](@ref) with the new objective
function.
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Min, 2x)
2 x

julia> @objective(model, Max, -2x)
-2 x
```

Alternatively, use [`set_objective_function`](@ref).

```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Min, 2x)
2 x

julia> new_objective = @expression(model, -2 * x)
-2 x

julia> set_objective_function(model, new_objective)
```

## Modify an objective coefficient

Use [`set_objective_coefficient`](@ref) to modify an objective coefficient.
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Min, 2x)
2 x

julia> set_objective_coefficient(model, x, 3)

julia> objective_function(model)
3 x
```

!!! info
    There is no way to modify the coefficient of a quadratic term. Set a new
    objective instead.

## Modify the objective sense

Use [`set_objective_sense`](@ref) to modify the objective sense.
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Min, 2x)
2 x

julia> objective_sense(model)
MIN_SENSE::OptimizationSense = 0

julia> set_objective_sense(model, MAX_SENSE);

julia> objective_sense(model)
MAX_SENSE::OptimizationSense = 1
```

Alternatively, call [`@objective`](@ref) and pass the existing objective
function.
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @objective(model, Min, 2x)
2 x

julia> @objective(model, Max, objective_function(model))
2 x
```
```@meta
DocTestSetup = quote
    using JuMP
end
```

# Expressions

JuMP has three types of expressions: affine, quadratic, and nonlinear. These
expressions can be inserted into constraints or into the objective. This is
particularly useful if an expression is used in multiple places in the model.

## Affine expressions

There are four ways of constructing an affine expression in JuMP: with the
[`@expression`](@ref) macro, with operator overloading, with the `AffExpr`
constructor, and with [`add_to_expression!`](@ref).

### Macros

The recommended way to create an affine expression is via the
[`@expression`](@ref) macro.

```jldoctest affine_macro
model = Model()
@variable(model, x)
@variable(model, y)
ex = @expression(model, 2x + y - 1)

# output

2 x + y - 1
```

This expression can be used in the objective or added to a constraint. For
example:
```jldoctest affine_macro
@objective(model, Min, 2 * ex - 1)
objective_function(model)

# output

4 x + 2 y - 3
```

Just like variables and constraints, named expressions can also be created. For
example
```jldoctest
model = Model()
@variable(model, x[i = 1:3])
@expression(model, expr[i = 1:3], i * sum(x[j] for j in i:3))
expr

# output

3-element Vector{AffExpr}:
 x[1] + x[2] + x[3]
 2 x[2] + 2 x[3]
 3 x[3]
```

!!! tip
    You can read more about containers in the [Containers](@ref) section.

### Operator overloading

Expressions can also be created without macros. However, note that in some
cases, this can be much slower that constructing an expression using macros.

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
ex = 2x + y - 1

# output

2 x + y - 1
```

### Constructors

A third way to create an affine expression is by the `AffExpr` constructor. The
first argument is the constant term, and the remaining arguments are
variable-coefficient pairs.

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
ex = AffExpr(-1.0, x => 2.0, y => 1.0)

# output

2 x + y - 1
```

### `add_to_expression!`

The fourth way to create an affine expression is by using
[`add_to_expression!`](@ref). Compared to the operator overloading method,
this approach is faster because it avoids constructing temporary objects.
The [`@expression`](@ref) macro uses [`add_to_expression!`](@ref)
behind-the-scenes.
```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
ex = AffExpr(-1.0)
add_to_expression!(ex, 2.0, x)
add_to_expression!(ex, 1.0, y)

# output

2 x + y - 1
```

!!! warning
    Read the section [Initializing arrays](@ref) for some cases to be careful
    about when using [`add_to_expression!`](@ref).

### Removing zero terms

Use [`drop_zeros!`](@ref) to remove terms from an affine expression with a `0`
coefficient.

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @expression(model, ex, x + 1 - x)
0 x + 1

julia> drop_zeros!(ex)

julia> ex
1
```

### Coefficients

Use [`coefficient`](@ref) to return the coefficient associated with a variable
in an affine expression.

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> @expression(model, ex, 2x + 1)
2 x + 1

julia> coefficient(ex, x)
2.0

julia> coefficient(ex, y)
0.0
```

## Quadratic expressions

Like affine expressions, there are four ways of constructing a quadratic
expression in JuMP: macros, operator overloading, constructors, and
[`add_to_expression!`](@ref).

### Macros

The [`@expression`](@ref) macro can be used to create quadratic expressions by
including quadratic terms.

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
ex = @expression(model, x^2 + 2 * x * y + y^2 + x + y - 1)

# output

x² + 2 y*x + y² + x + y - 1
```

### Operator overloading

Operator overloading can also be used to create quadratic expressions. The same
performance warning (discussed in the affine expression section) applies.

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
ex = x^2 + 2 * x * y + y^2 + x + y - 1

# output

x² + 2 x*y + y² + x + y - 1
```

### Constructors

Quadratic expressions can also be created using the `QuadExpr` constructor. The
first argument is an affine expression, and the remaining arguments are pairs,
where the first term is a `JuMP.UnorderedPair` and the second term is the
coefficient.

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
aff_expr = AffExpr(-1.0, x => 1.0, y => 1.0)
quad_expr = QuadExpr(aff_expr, UnorderedPair(x, x) => 1.0,
                     UnorderedPair(x, y) => 2.0, UnorderedPair(y, y) => 1.0)

# output

x² + 2 x*y + y² + x + y - 1
```

### `add_to_expression!`

Finally, [`add_to_expression!`](@ref) can also be used to add quadratic terms.

```jldoctest
model = Model()
@variable(model, x)
@variable(model, y)
ex = QuadExpr(x + y - 1.0)
add_to_expression!(ex, 1.0, x, x)
add_to_expression!(ex, 2.0, x, y)
add_to_expression!(ex, 1.0, y, y)

# output

x² + 2 x*y + y² + x + y - 1
```

!!! warning
    Read the section [Initializing arrays](@ref) for some cases to be careful
    about when using [`add_to_expression!`](@ref).

### Removing zero terms

Use [`drop_zeros!`](@ref) to remove terms from a quadratic expression with a `0`
coefficient.

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @expression(model, ex, x^2 + x + 1 - x^2)
0 x² + x + 1

julia> drop_zeros!(ex)

julia> ex
x + 1
```

### Coefficients

Use [`coefficient`](@ref) to return the coefficient associated with a pair of variables
in a quadratic expression.

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> @expression(model, ex, 2*x*y + 3*x)
2 x*y + 3 x

julia> coefficient(ex, x, y)
2.0

julia> coefficient(ex, x, x)
0.0

julia> coefficient(ex, y, x)
2.0

julia> coefficient(ex, x)
3.0
```

## Nonlinear expressions

Nonlinear expressions can be constructed only using the [`@NLexpression`](@ref)
macro and can be used only in [`@NLobjective`](@ref), [`@NLconstraint`](@ref),
and other [`@NLexpression`](@ref)s. Moreover, quadratic and affine expressions
cannot be used in the nonlinear macros. For more details, see the [Nonlinear
Modeling](@ref) section.

## Initializing arrays

JuMP implements `zero(AffExpr)` and `one(AffExpr)` to support various functions
in `LinearAlgebra` (for example, accessing the off-diagonal of a `Diagonal`
matrix).
```jldoctest
julia> zero(AffExpr)
0

julia> one(AffExpr)
1
```

However, this can result in a subtle bug if you call
[`add_to_expression!`](@ref) or the [MutableArithmetics API](https://github.com/jump-dev/MutableArithmetics.jl)
on an element created by `zeros` or `ones`:
```jldoctest
julia> x = zeros(AffExpr, 2)
2-element Vector{AffExpr}:
 0
 0

julia> add_to_expression!(x[1], 1.1)
1.1

julia> x
2-element Vector{AffExpr}:
 1.1
 1.1
```

Notice how we modified `x[1]`, but we also changed `x[2]`!

This happened because `zeros(AffExpr, 2)` calls `zero(AffExpr)` once to obtain a
zero element, and then creates an appropriately sized array filled with the same
element.

This also happens with broadcasting calls containing a conversion of `0` or `1`:
```jldoctest
julia> x = Vector{AffExpr}(undef, 2)
2-element Vector{AffExpr}:
 #undef
 #undef

julia> x .= 0
2-element Vector{AffExpr}:
 0
 0

julia> add_to_expression!(x[1], 1.1)
1.1

julia> x
2-element Vector{AffExpr}:
 1.1
 1.1
```

The recommended way to create an array of empty expressions is as follows:
```jldoctest
julia> x = Vector{AffExpr}(undef, 2)
2-element Vector{AffExpr}:
 #undef
 #undef

julia> for i in eachindex(x)
           x[i] = AffExpr(0.0)
       end

julia> add_to_expression!(x[1], 1.1)
1.1

julia> x
2-element Vector{AffExpr}:
 1.1
 0
```

Alternatively, use non-mutating operation to avoid updating `x[1]` in-place:
```jldoctest
julia> x = zeros(AffExpr, 2)
2-element Vector{AffExpr}:
 0
 0

julia> x[1] += 1.1
1.1

julia> x
2-element Vector{AffExpr}:
 1.1
 0
```
Note that for large expressions this will be slower due to the allocation of
additional temporary objects.
```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
    import GLPK
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# [Constraints](@id jump_constraints)

JuMP is based on the [MathOptInterface (MOI) API](@ref moi_documentation).
Because of this, JuMP uses the following standard form to represent problems:
```math
\begin{align}
    & \min_{x \in \mathbb{R}^n} & f_0(x)
    \\
    & \;\;\text{s.t.} & f_i(x) & \in \mathcal{S}_i & i = 1 \ldots m
\end{align}
```
Each constraint, ``f_i(x) \in \mathcal{S}_i``, is composed of a function and a
set. For example, instead of calling ``a^\top x \le b`` a
*less-than-or-equal-to* constraint, we say that it is a
*scalar-affine-in-less-than* constraint, where the function ``a^\top x`` belongs
to the *less-than* set ``(-\infty, b]``. We use the shorthand
*function-in-set* to refer to constraints composed of different types of
functions and sets.

This page explains how to write various types of constraints in JuMP. For
nonlinear constraints, see [Nonlinear Modeling](@ref) instead.

## Add a constraint

Add a constraint to a JuMP model using the [`@constraint`](@ref) macro. The
syntax to use depends on the type of constraint you wish to add.

### Add a linear constraint

Create linear constraints using the [`@constraint`](@ref) macro:
```jldoctest
model = Model()
@variable(model, x[1:3])
@constraint(model, c1, sum(x) <= 1)
@constraint(model, c2, x[1] + 2 * x[3] >= 2)
@constraint(model, c3, sum(i * x[i] for i in 1:3) == 3)
@constraint(model, c4, 4 <= 2 * x[2] <= 5)
print(model)

# output

Feasibility
Subject to
 c3 : x[1] + 2 x[2] + 3 x[3] = 3.0
 c2 : x[1] + 2 x[3] ≥ 2.0
 c1 : x[1] + x[2] + x[3] ≤ 1.0
 c4 : 2 x[2] ∈ [4.0, 5.0]
```

### Normalization

JuMP normalizes constraints by moving all of the terms containing variables to
the left-hand side and all of the constant terms to the right-hand side. Thus,
we get:
```jldoctest; setup=:(model=Model(); @variable(model, x))
julia> @constraint(model, c, 2x + 1 <= 4x + 4)
c : -2 x <= 3.0
```

### [Add a quadratic constraint](@id quad_constraints)

In addition to affine functions, JuMP also supports constraints with quadratic
terms. For example:
```jldoctest con_quadratic; setup=:(model=Model())
julia> @variable(model, x[i=1:2])
2-element Vector{VariableRef}:
 x[1]
 x[2]

julia> @variable(model, t >= 0)
t

julia> @constraint(model, my_q, x[1]^2 + x[2]^2 <= t^2)
my_q : x[1]² + x[2]² - t² <= 0.0
```

!!! tip
    Because solvers can take advantage of the knowledge that a constraint is
    quadratic, prefer adding quadratic constraints using [`@constraint`](@ref),
    rather than [`@NLconstraint`](@ref).

### Vectorized constraints

You can also add constraints to JuMP using vectorized linear algebra. For
example:
```jldoctest con_vector; setup=:(model = Model())
julia> @variable(model, x[i=1:2])
2-element Vector{VariableRef}:
 x[1]
 x[2]

julia> A = [1 2; 3 4]
2×2 Matrix{Int64}:
 1  2
 3  4

julia> b = [5, 6]
2-element Vector{Int64}:
 5
 6

julia> @constraint(model, con, A * x .== b)
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.EqualTo{Float64}}, ScalarShape}}:
 con : x[1] + 2 x[2] = 5.0
 con : 3 x[1] + 4 x[2] = 6.0
```

!!! note
    Make sure to use [Julia's dot syntax](https://docs.julialang.org/en/v1/manual/functions/index.html#man-vectorized-1)
    in front of the comparison operators (for example, `.==`, `.>=`, and `.<=`). If you
    use a comparison without the dot, an error will be thrown.

### Containers of constraints

The [`@constraint`](@ref) macro supports creating collections of constraints.
We'll cover some brief syntax here; read the [Constraint containers](@ref)
section for more details:

Create arrays of constraints:
```jldoctest; setup=:(model=Model(); @variable(model, x[1:3]))
julia> @constraint(model, c[i=1:3], x[i] <= i^2)
3-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 c[1] : x[1] ≤ 1.0
 c[2] : x[2] ≤ 4.0
 c[3] : x[3] ≤ 9.0

julia> c[2]
c[2] : x[2] ≤ 4.0
```

Sets can be any Julia type that supports iteration:
```jldoctest; setup=:(model=Model(); @variable(model, x[1:3]))
julia> @constraint(model, c[i=2:3, ["red", "blue"]], x[i] <= i^2)
2-dimensional DenseAxisArray{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape},2,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, ["red", "blue"]
And data, a 2×2 Matrix{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 c[2,red] : x[2] ≤ 4.0  c[2,blue] : x[2] ≤ 4.0
 c[3,red] : x[3] ≤ 9.0  c[3,blue] : x[3] ≤ 9.0

julia> c[2, "red"]
c[2,red] : x[2] ≤ 4.0
```

Sets can depend upon previous indices:
```jldoctest; setup=:(model=Model(); @variable(model, x[1:3]))
julia> @constraint(model, c[i=1:3, j=i:3], x[i] <= j)
JuMP.Containers.SparseAxisArray{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}, 2, Tuple{Int64, Int64}} with 6 entries:
  [1, 1]  =  c[1,1] : x[1] ≤ 1.0
  [1, 2]  =  c[1,2] : x[1] ≤ 2.0
  [1, 3]  =  c[1,3] : x[1] ≤ 3.0
  [2, 2]  =  c[2,2] : x[2] ≤ 2.0
  [2, 3]  =  c[2,3] : x[2] ≤ 3.0
  [3, 3]  =  c[3,3] : x[3] ≤ 3.0
```
and you can filter elements in the sets using the `;` syntax:
```jldoctest; setup=:(model=Model(); @variable(model, x[1:9]))
julia> @constraint(model, c[i=1:9; mod(i, 3) == 0], x[i] <= i)
JuMP.Containers.SparseAxisArray{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}, 1, Tuple{Int64}} with 3 entries:
  [3]  =  c[3] : x[3] ≤ 3.0
  [6]  =  c[6] : x[6] ≤ 6.0
  [9]  =  c[9] : x[9] ≤ 9.0
```

## Registered constraints

When you create constraints, JuMP registers them inside the model using their
corresponding symbol. Get a registered name using `model[:key]`:
```jldoctest
julia> model = Model()
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> @variable(model, x)
x

julia> @constraint(model, my_c, 2x <= 1)
my_c : 2 x ≤ 1.0

julia> model
A JuMP Model
Feasibility problem with:
Variable: 1
`AffExpr`-in-`MathOptInterface.LessThan{Float64}`: 1 constraint
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
Names registered in the model: my_c, x

julia> model[:my_c] === my_c
true
```

## Anonymous constraints

To reduce the likelihood of accidental bugs, and because JuMP registers
constraints inside a model, creating two constraints with the same name is an
error:
```julia
julia> model = Model();

julia> @variable(model, x)
x

julia> @constraint(model, c, 2x <= 1)
c : 2 x <= 1.0

julia> @constraint(model, c, 2x <= 1)
ERROR: An object of name c is already attached to this model. If this
    is intended, consider using the anonymous construction syntax, e.g.,
    `x = @variable(model, [1:N], ...)` where the name of the object does
    not appear inside the macro.

    Alternatively, use `unregister(model, :c)` to first unregister
    the existing name from the model. Note that this will not delete the
    object; it will just remove the reference at `model[:c]`.
[...]
```

A common reason for encountering this error is adding constraints in a loop.

As a work-around, JuMP provides *anonymous* constraints. Create an anonymous
constraint by omitting the name argument:
```jldoctest; setup=:(model=Model(); @variable(model, x))
julia> c = @constraint(model, 2x <= 1)
2 x <= 1.0
```

Create a container of anonymous constraints by dropping the name in front of
the `[`:
```jldoctest; setup=:(model=Model(); @variable(model, x[1:3]))
julia> c = @constraint(model, [i = 1:3], x[i] <= i)
3-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 x[1] ≤ 1.0
 x[2] ≤ 2.0
 x[3] ≤ 3.0
```

## Constraint names

In addition to the symbol that constraints are registered with, constraints have
a `String` name that is used for printing and writing to file formats.

Get and set the name of a constraint using [`name(::JuMP.ConstraintRef)`](@ref)
and [`set_name(::JuMP.ConstraintRef, ::String)`](@ref):
```jldoctest
julia> model = Model(); @variable(model, x);

julia> @constraint(model, con, x <= 1)
con : x <= 1.0

julia> name(con)
"con"

julia> set_name(con, "my_con_name")

julia> con
my_con_name : x <= 1.0
```

Override the default choice of name using the `base_name` keyword:
```jldoctest constraint_name_vector
julia> model = Model(); @variable(model, x);

julia> con = @constraint(model, [i=1:2], x <= i, base_name = "my_con")
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 my_con[1] : x ≤ 1.0
 my_con[2] : x ≤ 2.0
```

Note that names apply to each element of the container, not to the container of
constraints:
```jldoctest constraint_name_vector
julia> name(con[1])
"my_con[1]"

julia> set_name(con[1], "c")

julia> con
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 c : x ≤ 1.0
 my_con[2] : x ≤ 2.0
```

### Retrieve a constraint by name

Retrieve a constraint from a model using [`constraint_by_name`](@ref):
```jldoctest constraint_name_vector
julia> constraint_by_name(model, "c")
c : x ≤ 1.0
```

If the name is not present, `nothing` will be returned:
```jldoctest constraint_name_vector
julia> constraint_by_name(model, "bad_name")
```

You can only look up individual constraints using [`constraint_by_name`](@ref).
Something like this will not work:
```jldoctest
julia> model = Model(); @variable(model, x);

julia> con = @constraint(model, [i=1:2], x <= i, base_name = "my_con")
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 my_con[1] : x ≤ 1.0
 my_con[2] : x ≤ 2.0

julia> constraint_by_name(model, "my_con")
```

To look up a collection of constraints, do not use [`constraint_by_name`](@ref).
Instead, register them using the `model[:key] = value` syntax:
```jldoctest
julia> model = Model(); @variable(model, x);

julia> model[:con] = @constraint(model, [i=1:2], x <= i, base_name = "my_con")
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 my_con[1] : x ≤ 1.0
 my_con[2] : x ≤ 2.0

julia> model[:con]
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 my_con[1] : x ≤ 1.0
 my_con[2] : x ≤ 2.0
```

## String names, symbolic names, and bindings

It's common for new users to experience confusion relating to constraints.
Part of the problem is the difference between the name that a constraint is
registered under and the `String` name used for printing.

Here's a summary of the differences:

 * Constraints are created using [`@constraint`](@ref).
 * Constraints can be named or anonymous.
 * Named constraints have the form `@constraint(model, c, expr)`. For named
   constraints:
   * The `String` name of the constraint is set to `"c"`.
   * A Julia variable `c` is created that binds `c` to  the JuMP constraint.
   * The name `:c` is registered as a key in the model with the value `c`.
 * Anonymous constraints have the form `c = @constraint(model, expr)`. For
   anonymous constraints:
   * The `String` name of the constraint is set to `""`.
   * You control the name of the Julia variable used as the binding.
   * No name is registered as a key in the model.
 * The `base_name` keyword can override the `String` name of the constraint.
 * You can manually register names in the model via `model[:key] = value`.

Here's an example of the differences:
```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> c_binding = @constraint(model, 2x <= 1, base_name = "c")
c : 2 x <= 1.0

julia> model
A JuMP Model
Feasibility problem with:
Variable: 1
`AffExpr`-in-`MathOptInterface.LessThan{Float64}`: 1 constraint
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
Names registered in the model: x

julia> c
ERROR: UndefVarError: c not defined

julia> c_binding
c : 2 x <= 1.0

julia> name(c_binding)
"c"

julia> model[:c_register] = c_binding
c : 2 x <= 1.0

julia> model
A JuMP Model
Feasibility problem with:
Variable: 1
`AffExpr`-in-`MathOptInterface.LessThan{Float64}`: 1 constraint
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
Names registered in the model: c_register, x

julia> model[:c_register]
c : 2 x <= 1.0

julia> model[:c_register] === c_binding
true

julia> c
ERROR: UndefVarError: c not defined
```

## The `@constraints` macro

If you have many [`@constraint`](@ref) calls, use the [`@constraints`](@ref)
macro to improve readability:
```jldoctest; setup=:(model=Model(); @variable(model, x))
julia> @constraints(model, begin
           2x <= 1
           c, x >= -1
       end)
(2 x ≤ 1.0, c : x ≥ -1.0)

julia> print(model)
Feasibility
Subject to
 c : x ≥ -1.0
 2 x ≤ 1.0
```
The [`@constraints`](@ref) macro returns a tuple of the constraints that were
defined.

## [Duality](@id constraint_duality)

JuMP adopts the notion of [conic duality from MathOptInterface](@ref Duality).
For linear programs, a feasible dual on a `>=` constraint is nonnegative and a
feasible dual on a `<=` constraint is nonpositive. If the constraint is an
equality constraint, it depends on which direction is binding.

!!! warning
    JuMP's definition of duality is independent of the objective sense. That is,
    the sign of feasible duals associated with a constraint depends on the
    direction of the constraint and not whether the problem is maximization or
    minimization. **This is a different convention from linear programming
    duality in some common textbooks.** If you have a linear program, and you
    want the textbook definition, you probably want to use [`shadow_price`](@ref)
    and [`reduced_cost`](@ref) instead.

The dual value associated with a constraint in the most recent solution can be
accessed using the [`dual`](@ref) function. Use [`has_duals`](@ref) to check if
the model has a dual solution available to query. For example:
```jldoctest con_duality
julia> model = Model(GLPK.Optimizer);

julia> @variable(model, x)
x

julia> @constraint(model, con, x <= 1)
con : x <= 1.0

julia> @objective(model, Min, -2x)
-2 x

julia> has_duals(model)
false

julia> optimize!(model)

julia> has_duals(model)
true

julia> dual(con)
-2.0

julia> @objective(model, Max, 2x)
2 x

julia> optimize!(model)

julia> dual(con)
-2.0
```

To help users who may be less familiar with conic duality, JuMP provides
[`shadow_price`](@ref), which returns a value that can be interpreted as the
improvement in the objective in response to an infinitesimal relaxation (on the
scale of one unit) in the right-hand side of the constraint.
[`shadow_price`](@ref) can be used only on linear constraints with a `<=`, `>=`,
or `==` comparison operator.

In the example above, `dual(con)` returned `-2.0` regardless of the optimization
sense. However, in the second case when the optimization sense is `Max`,
[`shadow_price`](@ref) returns:
```jldoctest con_duality
julia> shadow_price(con)
2.0
```

### Duals of variable bounds

To query the dual variables associated with a variable bound, first obtain a
constraint reference using one of [`UpperBoundRef`](@ref),
[`LowerBoundRef`](@ref), or [`FixRef`](@ref), and then call [`dual`](@ref) on
the returned constraint reference. The [`reduced_cost`](@ref) function may
simplify this process as it returns the shadow price of an active bound of
a variable (or zero, if no active bound exists).
```jldoctest
julia> model = Model(GLPK.Optimizer);

julia> @variable(model, x <= 1)
x

julia> @objective(model, Min, -2x)
-2 x

julia> optimize!(model)

julia> dual(UpperBoundRef(x))
-2.0

julia> reduced_cost(x)
-2.0
```

## Modify a constant term

This section explains how to modify the constant term in a constraint. There are
multiple ways to achieve this goal; we explain three options.

### Option 1: change the right-hand side

Use [`set_normalized_rhs`](@ref) to modify the right-hand side (constant)
term of a linear or quadratic  constraint. Use [`normalized_rhs`](@ref) to query
the right-hand side term.
```jldoctest; setup = :(model = Model(); @variable(model, x))
julia> @constraint(model, con, 2x <= 1)
con : 2 x <= 1.0

julia> set_normalized_rhs(con, 3)

julia> con
con : 2 x <= 3.0

julia> normalized_rhs(con)
3.0
```

!!! warning
    [`set_normalized_rhs`](@ref) sets the right-hand side term of the
    normalized constraint. See [Normalization](@ref) for more details.

### Option 2: use fixed variables

If constraints are complicated, for example, they are composed of a number of
components, each of which has a constant term, then it may be difficult to
calculate what the right-hand side term is in the standard form.

For this situation, JuMP includes the ability to *fix* variables to a
value using the [`fix`](@ref) function. Fixing a variable sets its lower
and upper bound to the same value. Thus, changes in a constant term can be
simulated by adding a new variable and fixing it to different values. Here is
an example:
```jldoctest; setup = :(model = Model(); @variable(model, x))
julia> @variable(model, const_term)
const_term

julia> @constraint(model, con, 2x <= const_term + 1)
con : 2 x - const_term <= 1.0

julia> fix(const_term, 1.0)
```
The constraint `con` is now equivalent to `2x <= 2`.

!!! warning
    Fixed variables are not replaced with constants when communicating the
    problem to a solver. Therefore, even though `const_term` is fixed, it is
    still a decision variable, and so `const_term * x` is bilinear.

### Option 3: modify the function's constant term

The third option is to use [`add_to_function_constant`](@ref). The constant
given is added to the function of a `func`-in-`set` constraint. In the following
example, adding `2` to the function has the effect of removing `2` to the
right-hand side:
```jldoctest con_add; setup = :(model = Model(); @variable(model, x))
julia> @constraint(model, con, 2x <= 1)
con : 2 x <= 1.0

julia> add_to_function_constant(con, 2)

julia> con
con : 2 x <= -1.0

julia> normalized_rhs(con)
-1.0
```

In the case of interval constraints, the constant is removed from each bound:
```jldoctest con_add_interval; setup = :(model = Model(); @variable(model, x))
julia> @constraint(model, con, 0 <= 2x + 1 <= 2)
con : 2 x ∈ [-1.0, 1.0]

julia> add_to_function_constant(con, 3)

julia> con
con : 2 x ∈ [-4.0, -2.0]
```

## Modify a variable coefficient

To modify the coefficients for a linear term (modifying the coefficient of a
quadratic term is not supported) in a constraint, use
[`set_normalized_coefficient`](@ref). To query the current coefficient, use
[`normalized_coefficient`](@ref).
```jldoctest; setup = :(model = Model(); @variable(model, x[1:2]))
julia> @constraint(model, con, 2x[1] + x[2] <= 1)
con : 2 x[1] + x[2] ≤ 1.0

julia> set_normalized_coefficient(con, x[2], 0)

julia> con
con : 2 x[1] ≤ 1.0

julia> normalized_coefficient(con, x[2])
0.0
```

!!! warning
    [`set_normalized_coefficient`](@ref) sets the coefficient of the normalized
    constraint. See [Normalization](@ref) for more details.

## Delete a constraint

Use [`delete`](@ref) to delete a constraint from a model. Use [`is_valid`](@ref)
to check if a constraint belongs to a model and has not been deleted.

```jldoctest constraints_delete; setup = :(model=Model(); @variable(model, x))
julia> @constraint(model, con, 2x <= 1)
con : 2 x <= 1.0

julia> is_valid(model, con)
true

julia> delete(model, con)

julia> is_valid(model, con)
false
```

Deleting a constraint does not unregister the symbolic reference from the model.
Therefore, creating a new constraint of the same name will throw an error:
```jldoctest constraints_delete
julia> @constraint(model, con, 2x <= 1)
ERROR: An object of name con is already attached to this model. If this
    is intended, consider using the anonymous construction syntax, e.g.,
    `x = @variable(model, [1:N], ...)` where the name of the object does
    not appear inside the macro.

    Alternatively, use `unregister(model, :con)` to first unregister
    the existing name from the model. Note that this will not delete the
    object; it will just remove the reference at `model[:con]`.
[...]
```

After calling [`delete`](@ref), call [`unregister`](@ref) to remove the symbolic
reference:
```jldoctest constraints_delete
julia> unregister(model, :con)

julia> @constraint(model, con, 2x <= 1)
con : 2 x <= 1.0
```

!!! info
    [`delete`](@ref) does not automatically [`unregister`](@ref) because we do
    not distinguish between names that are automatically registered by JuMP
    macros, and names that are manually registered by the user by setting values
    in [`object_dictionary`](@ref). In addition, deleting a constraint and then
    adding a new constraint of the same name is an easy way to introduce bugs
    into your code.

## Start values

Provide a starting value (also called warmstart) for a constraint's dual using
[`set_dual_start_value`](@ref).

The start value of a constraint's dual can be queried using
[`dual_start_value`](@ref). If no start value has been set,
[`dual_start_value`](@ref) will return `nothing`.

```jldoctest constraint_dual_start; setup=:(model=Model())
julia> @variable(model, x)
x

julia> @constraint(model, con, x >= 10)
con : x ≥ 10.0

julia> dual_start_value(con)

julia> set_dual_start_value(con, 2)

julia> dual_start_value(con)
2.0
```

Vector-valued constraints require a vector warmstart:
```jldoctest constraint_dual_start_vector; setup=:(model=Model())
julia> @variable(model, x[1:3])
3-element Vector{VariableRef}:
 x[1]
 x[2]
 x[3]

julia> @constraint(model, con, x in SecondOrderCone())
con : [x[1], x[2], x[3]] in MathOptInterface.SecondOrderCone(3)

julia> dual_start_value(con)

julia> set_dual_start_value(con, [1.0, 2.0, 3.0])

julia> dual_start_value(con)
3-element Vector{Float64}:
 1.0
 2.0
 3.0
```

To take the dual solution from the last solve and use it as the starting point
for a new solve, use:
```julia
for (F, S) in list_of_constraint_types(model)
    for con in all_constraints(model, F, S)
        set_dual_start_value(con, dual(con))
    end
end
```

!!! note
    Some constraints might not have well defined duals, hence you might need to
    filter `(F, S)` pairs.

## Constraint containers

Like [Variable containers](@ref), JuMP provides a mechanism for building groups
of constraints compactly. References to these groups of constraints are returned
in *containers*. Three types of constraint containers are supported: `Array`s,
`DenseAxisArray`s, and `SparseAxisArray`s. We explain each of these in the
following.

!!! tip
    You can read more about containers in the [Containers](@ref) section.

### [Arrays](@id constraint_arrays)

One way of adding a group of constraints compactly is the following:
```jldoctest constraint_arrays; setup=:(model=Model(); @variable(model, x))
julia> @constraint(model, con[i = 1:3], i * x <= i + 1)
3-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 con[1] : x ≤ 2.0
 con[2] : 2 x ≤ 3.0
 con[3] : 3 x ≤ 4.0
```
JuMP returns references to the three constraints in an `Array` that is bound to
the Julia variable `con`. This array can be accessed and sliced as you would
with any Julia array:
```jldoctest constraint_arrays
julia> con[1]
con[1] : x <= 2.0

julia> con[2:3]
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 con[2] : 2 x ≤ 3.0
 con[3] : 3 x ≤ 4.0
```

Anonymous containers can also be constructed by dropping the name (for example, `con`)
before the square brackets:
```jldoctest constraint_arrays
julia> con = @constraint(model, [i = 1:2], i * x <= i + 1)
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 x ≤ 2.0
 2 x ≤ 3.0
```

Just like [`@variable`](@ref), JuMP will form an `Array` of constraints when it
can determine at parse time that the indices are one-based integer ranges.
Therefore `con[1:b]` will create an `Array`, but `con[a:b]` will not. A special
case is `con[Base.OneTo(n)]` which will produce an `Array`. If JuMP cannot
determine that the indices are one-based integer ranges (for example, in the case of
`con[a:b]`), JuMP will create a `DenseAxisArray` instead.

### DenseAxisArrays

The syntax for constructing a [`DenseAxisArray`](@ref Containers.DenseAxisArray)
of constraints is very similar to the
[syntax for constructing](@ref variable_jump_arrays) a `DenseAxisArray` of
variables.

```jldoctest constraint_jumparrays; setup=:(model=Model(); @variable(model, x))
julia> @constraint(model, con[i = 1:2, j = 2:3], i * x <= j + 1)
2-dimensional DenseAxisArray{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape},2,...} with index sets:
    Dimension 1, Base.OneTo(2)
    Dimension 2, 2:3
And data, a 2×2 Matrix{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 con[1,2] : x ≤ 3.0    con[1,3] : x ≤ 4.0
 con[2,2] : 2 x ≤ 3.0  con[2,3] : 2 x ≤ 4.0
```

### SparseAxisArrays

The syntax for constructing a
[`SparseAxisArray`](@ref Containers.SparseAxisArray) of constraints is very
similar to the [syntax for constructing](@ref variable_sparseaxisarrays) a
`SparseAxisArray` of variables.

```jldoctest constraint_jumparrays; setup=:(model=Model(); @variable(model, x))
julia> @constraint(model, con[i = 1:2, j = 1:2; i != j], i * x <= j + 1)
JuMP.Containers.SparseAxisArray{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}, 2, Tuple{Int64, Int64}} with 2 entries:
  [1, 2]  =  con[1,2] : x ≤ 3.0
  [2, 1]  =  con[2,1] : 2 x ≤ 2.0
```

!!! warning
    If you have many index dimensions and a large amount of sparsity, read
    [Performance considerations](@ref).

### Forcing the container type

When creating a container of constraints, JuMP will attempt to choose the
tightest container type that can store the constraints. However, because this
happens at parse time, it does not always make the best choice. Just like in
[`@variable`](@ref), you can force the type of container using the `container`
keyword. For syntax and the reason behind this, take a look at the
[variable docs](@ref variable_forcing).

## Accessing constraints from a model

Query the types of function-in-set constraints in a model using
[`list_of_constraint_types`](@ref):
```jldoctest con_access
julia> model = Model();

julia> @variable(model, x[i=1:2] >= i, Int);

julia> @constraint(model, x[1] + x[2] <= 1);

julia> list_of_constraint_types(model)
3-element Vector{Tuple{Type, Type}}:
 (AffExpr, MathOptInterface.LessThan{Float64})
 (VariableRef, MathOptInterface.GreaterThan{Float64})
 (VariableRef, MathOptInterface.Integer)
```

For a given combination of function and set type, use
[`num_constraints`](@ref) to access the number of constraints and
[`all_constraints`](@ref) to access a list of their references:
```jldoctest con_access
julia> num_constraints(model, VariableRef, MOI.Integer)
2

julia> cons = all_constraints(model, VariableRef, MOI.Integer)
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.VariableIndex, MathOptInterface.Integer}, ScalarShape}}:
 x[1] integer
 x[2] integer
```

Use [`constraint_object`](@ref) to get an instance of an
[`AbstractConstraint`](@ref) object that stores the constraint data:
```jldoctest con_access
julia> con = constraint_object(cons[1])
ScalarConstraint{VariableRef, MathOptInterface.Integer}(x[1], MathOptInterface.Integer())

julia> con.func
x[1]

julia> con.set
MathOptInterface.Integer()
```

## MathOptInterface constraints

Because JuMP is based on MathOptInterface, you can add any constraints supported
by MathOptInterface using the function-in-set syntax. For a list of supported
functions and sets, read [Standard form problem](@ref).

!!! note
    We use `MOI` as an alias for the `MathOptInterface` module. This alias is
    defined by `using JuMP`. You may also define it in your code as follows:
    ```julia
    import MathOptInterface
    const MOI = MathOptInterface
    ```

For example, the following two constraints are equivalent:
```jldoctest moi; setup=:(model=Model(); @variable(model, x[1:3]))
julia> @constraint(model, 2 * x[1] <= 1)
2 x[1] ≤ 1.0

julia> @constraint(model, 2 * x[1] in MOI.LessThan(1.0))
2 x[1] ≤ 1.0
```

You can also use any set defined by MathOptInterface:
```jldoctest moi
julia> @constraint(model, x - [1; 2; 3] in MOI.Nonnegatives(3))
[x[1] - 1, x[2] - 2, x[3] - 3] ∈ MathOptInterface.Nonnegatives(3)

julia> @constraint(model, x in MOI.ExponentialCone())
[x[1], x[2], x[3]] ∈ MathOptInterface.ExponentialCone()
```

!!! info
    Similar to how JuMP defines the `<=` and `>=` syntax as a convenience way to
    specify [`MOI.LessThan`](@ref) and [`MOI.GreaterThan`](@ref) constraints,
    the remaining sections in this page describe functions and syntax that have
    been added for the convenience of common modeling situations.

## Set inequality syntax

For modeling convenience, the syntax `@constraint(model, x >= y, Set())` is
short-hand for `@constraint(model, x - y in Set())`. Therefore, the following
calls are equivalent:
```jldoctest set_inequality
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> y = [0.5, 0.75];

julia> @constraint(model, x >= y, MOI.Nonnegatives(2))
[x[1] - 0.5, x[2] - 0.75] ∈ MathOptInterface.Nonnegatives(2)

julia> @constraint(model, y <= x, MOI.Nonnegatives(2))
[x[1] - 0.5, x[2] - 0.75] ∈ MathOptInterface.Nonnegatives(2)

julia> @constraint(model, x - y in MOI.Nonnegatives(2))
[x[1] - 0.5, x[2] - 0.75] ∈ MathOptInterface.Nonnegatives(2)
```

Non-zero constants are not supported in this syntax:
```jldoctest set_inequality
julia> @constraint(model, x >= 1, MOI.Nonnegatives(2))
ERROR: Operation `sub_mul` between `Vector{VariableRef}` and `Int64` is not allowed. You should use broadcast.
Stacktrace:
[...]
```
Use instead:
```jldoctest set_inequality
julia> @constraint(model, x .- 1 >= 0, MOI.Nonnegatives(2))
[x[1] - 1, x[2] - 1] ∈ MathOptInterface.Nonnegatives(2)
```

## Second-order cone constraints

A [`SecondOrderCone`](@ref) constrains the variables `t` and `x` to the set:
```math
||x||_2 \le t,
```
and ``t \ge 0``. It can be added as follows:
```jldoctest
julia> model = Model();

julia> @variable(model, t)
t

julia> @variable(model, x[1:2])
2-element Vector{VariableRef}:
 x[1]
 x[2]

julia> @constraint(model, [t; x] in SecondOrderCone())
[t, x[1], x[2]] ∈ MathOptInterface.SecondOrderCone(3)
```

## Rotated second-order cone constraints

A [`RotatedSecondOrderCone`](@ref) constrains the variables `t`, `u`, and `x`
to the set:
```math
||x||_2^2 \le 2 t \cdot u
```
and ``t, u \ge 0``. It can be added as follows:
```jldoctest
julia> model = Model();

julia> @variable(model, t)
t

julia> @variable(model, u)
u

julia> @variable(model, x[1:2])
2-element Vector{VariableRef}:
 x[1]
 x[2]

julia> @constraint(model, [t; u; x] in RotatedSecondOrderCone())
[t, u, x[1], x[2]] ∈ MathOptInterface.RotatedSecondOrderCone(4)
```

## Semi-integer and semi-continuous variables

Semi-continuous variables are constrained to the set
``x \in \{0\} \cup [l, u]``.

Create a semi-continuous variable using the `MOI.Semicontinuous` set:
```jldoctest; setup = :(model = Model(); @variable(model, x))
julia> @constraint(model, x in MOI.Semicontinuous(1.5, 3.5))
x in MathOptInterface.Semicontinuous{Float64}(1.5, 3.5)
```

Semi-integer variables  are constrained to the set
``x \\in \{0\} \cup \{l, l+1, \dots, u\}``.

Create a semi-integer variable using the `MOI.Semiinteger` set:
```jldoctest; setup = :(model = Model(); @variable(model, x))
julia> @constraint(model, x in MOI.Semiinteger(1.0, 3.0))
x in MathOptInterface.Semiinteger{Float64}(1.0, 3.0)
```

## Special Ordered Sets of Type 1

In a Special Ordered Set of Type 1 (often denoted SOS-I or SOS1), at most one
element can take a non-zero value.

Construct SOS-I constraints using the [`SOS1`](@ref) set:
```jldoctest con_sos; setup=:(model = Model())
julia> @variable(model, x[1:3])
3-element Vector{VariableRef}:
 x[1]
 x[2]
 x[3]

julia> @constraint(model, x in SOS1())
[x[1], x[2], x[3]] in MathOptInterface.SOS1{Float64}([1.0, 2.0, 3.0])
```

Although not required for feasibility, solvers can benefit from an ordering of
the variables (for example, the variables represent different factories to build, at
most one factory can be built, and the factories can be ordered according to
cost). To induce an ordering, a vector of weights can be provided, and the
variables are ordered according to their corresponding weight.

For example, in the constraint:
```jldoctest con_sos
julia> @constraint(model, x in SOS1([3.1, 1.2, 2.3]))
[x[1], x[2], x[3]] in MathOptInterface.SOS1{Float64}([3.1, 1.2, 2.3])
```
the variables `x` have precedence `x[2]`, `x[3]`, `x[1]`.

## Special Ordered Sets of Type 2

In a Special Ordered Set of Type 2 (SOS-II), at most two elements can be
non-zero, and if there are two non-zeros, they must be consecutive according to
the ordering induced by a weight vector.

Construct SOS-II constraints using the [`SOS2`](@ref) set:
```jldoctest con_sos
julia> @constraint(model, x in SOS2([3.0, 1.0, 2.0]))
[x[1], x[2], x[3]] in MathOptInterface.SOS2{Float64}([3.0, 1.0, 2.0])
```
The possible non-zero pairs are (`x[1]`, `x[3]`) and (`x[2]`, `x[3]`):

If the weight vector is omitted, JuMP induces an ordering from `1:length(x)`:
```jldoctest con_sos
julia> @constraint(model, x in SOS2())
[x[1], x[2], x[3]] in MathOptInterface.SOS2{Float64}([1.0, 2.0, 3.0])
```

## Indicator constraints

Indicator constraints consist of a binary variable and a linear constraint. The
constraint holds when the binary variable takes the value `1`. The constraint
may or may not hold when the binary variable takes the value `0`.

To enforce the constraint `x + y <= 1` when the binary variable `a` is `1`, use:
```jldoctest indicator; setup=:(model = Model())
julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> @variable(model, a, Bin)
a

julia> @constraint(model, a => {x + y <= 1})
a => {x + y ≤ 1.0}
```

If the constraint must hold when `a` is zero, add `!` or `¬` before the binary
variable;
```jldoctest indicator
julia> @constraint(model, !a => {x + y <= 1})
!a => {x + y ≤ 1.0}
```

## Semidefinite constraints

To constrain a matrix to be positive semidefinite (PSD), use [`PSDCone`](@ref):
```jldoctest con_psd; setup=:(model = Model())
julia> @variable(model, X[1:2, 1:2])
2×2 Matrix{VariableRef}:
 X[1,1]  X[1,2]
 X[2,1]  X[2,2]

julia> @constraint(model, X >= 0, PSDCone())
[X[1,1]  X[1,2];
 X[2,1]  X[2,2]] ∈ PSDCone()
```

!!! tip
    Where possible, prefer constructing a matrix of
    [Semidefinite variables](@ref) using the [`@variable`](@ref) macro, rather
    than adding a constraint like `@constraint(model, X >= 0, PSDCone())`. In
    some solvers, adding the constraint via [`@constraint`](@ref) is less
    efficient, and can result in additional intermediate variables and
    constraints being added to the model.

The inequality `X >= Y` between two square matrices `X` and `Y` is understood as
constraining `X - Y` to be positive semidefinite.
```jldoctest con_psd
julia> Y = [1 2; 2 1]
2×2 Matrix{Int64}:
 1  2
 2  1

julia> @constraint(model, X >= Y, PSDCone())
[X[1,1] - 1  X[1,2] - 2;
 X[2,1] - 2  X[2,2] - 1] ∈ PSDCone()
```

### Symmetry

Solvers supporting PSD constraints usually expect to be given a matrix that
is *symbolically* symmetric, that is, for which the expression in corresponding
off-diagonal entries are the same. In our example, the expressions of entries
`(1, 2)` and `(2, 1)` are respectively `X[1,2] - 2` and `X[2,1] - 2` which are
different.

To bridge the gap between the constraint modeled and what the solver
expects, solvers may add an equality constraint `X[1,2] - 2 == X[2,1] - 2` to
force symmetry. Use `LinearAlgebra.Symmetric` to explicitly tell the solver that
the matrix is symmetric:
```jldoctest con_psd
julia> import LinearAlgebra

julia> Z = [X[1, 1] X[1, 2]; X[1, 2] X[2, 2]]
2×2 Matrix{VariableRef}:
 X[1,1]  X[1,2]
 X[1,2]  X[2,2]

julia> @constraint(model, LinearAlgebra.Symmetric(Z) >= 0, PSDCone())
[X[1,1]  X[1,2];
 X[1,2]  X[2,2]] ∈ PSDCone()
```

Note that the lower triangular entries are ignored even if they are
different so use it with caution:
```jldoctest con_psd
julia> @constraint(model, LinearAlgebra.Symmetric(X) >= 0, PSDCone())
[X[1,1]  X[1,2];
 X[1,2]  X[2,2]] ∈ PSDCone()
```
(Note the `(2, 1)` element of the constraint is `X[1,2]`, not `X[2,1]`.)

## Complementarity constraints

A mixed complementarity constraint `F(x) ⟂ x` consists of finding `x` in the
interval `[lb, ub]`, such that the following holds:

- `F(x) == 0` if `lb < x < ub`
- `F(x) >= 0` if `lb == x`
- `F(x) <= 0` if `x == ub`

JuMP supports mixed complementarity constraints via `complements(F(x), x)` or
`F(x) ⟂ x` in the [`@constraint`](@ref) macro. The interval set `[lb, ub]` is
obtained from the variable bounds on `x`.

For example, to define the problem `2x - 1 ⟂ x` with `x ∈ [0, ∞)`, do:
```jldoctest complementarity; setup=:(model=Model())
julia> @variable(model, x >= 0)
x

julia> @constraint(model, 2x - 1 ⟂ x)
[2 x - 1, x] ∈ MathOptInterface.Complements(2)
```
This problem has a unique solution at `x = 0.5`.

The perp operator `⟂` can be entered in most editors (and the Julia REPL) by
typing `\perp<tab>`.

An alternative approach that does not require the `⟂` symbol uses the
`complements` function as follows:
```jldoctest complementarity
julia> @constraint(model, complements(2x - 1, x))
[2 x - 1, x] ∈ MathOptInterface.Complements(2)
```

In both cases, the mapping `F(x)` is supplied as the first argument, and the
matching variable `x` is supplied as the second.

Vector-valued complementarity constraints are also supported:
```jldoctest complementarity
julia> @variable(model, -2 <= y[1:2] <= 2)
2-element Vector{VariableRef}:
 y[1]
 y[2]

julia> M = [1 2; 3 4]
2×2 Matrix{Int64}:
 1  2
 3  4

julia> q = [5, 6]
2-element Vector{Int64}:
 5
 6

julia> @constraint(model, M * y + q ⟂ y)
[y[1] + 2 y[2] + 5, 3 y[1] + 4 y[2] + 6, y[1], y[2]] ∈ MathOptInterface.Complements(4)
```
# [Callbacks](@id CallbackAPI)

More information can be found in the [Callbacks](@ref callbacks_manual) section
of the manual.

## Macros

```@docs
@build_constraint
```

## Callback variable primal

```@docs
callback_value
```

## Callback node status

```@docs
callback_node_status
```
# [Extensions](@id ExtensionAPI)

More information can be found in the [Extensions](@ref extensions_manual)
section of the manual.

## Define a new set

```@docs
AbstractVectorSet
```

## Extend `@variable`

```@docs
ScalarVariable
VariableInfo
add_variable
build_variable
```

## Extend `@constraint`

```@docs
build_constraint
add_constraint
AbstractShape
shape
reshape_vector
reshape_set
dual_shape
ScalarShape
VectorShape
SquareMatrixShape
SymmetricMatrixShape
operator_to_set
parse_constraint
parse_constraint_head
parse_constraint_call
```
# [Objectives](@id ObjectiveAPI)

More information can be found in the [Objectives](@ref) section of the manual.

## Objective functions

```@docs
@objective
objective_function
set_objective_function
set_objective_coefficient
set_objective
objective_function_type
objective_function_string
show_objective_function_summary
```

## Objective sense

```@docs
objective_sense
set_objective_sense
```
# [Models](@id ModelAPI)

More information can be found in the [Models](@ref jump_models) section of
the manual.

## Constructors

```@docs
Model
direct_model
```

## Enums

```@docs
ModelMode
AUTOMATIC
MANUAL
DIRECT
```

## Basic functions

```@docs
backend
unsafe_backend
name(::AbstractModel)
solver_name
Base.empty!(::Model)
Base.isempty(::Model)
mode
object_dictionary
unregister
latex_formulation
```

## Working with attributes

```@docs
set_optimizer
optimizer_with_attributes
get_optimizer_attribute
set_optimizer_attribute
set_optimizer_attributes
set_silent
unset_silent
set_time_limit_sec
unset_time_limit_sec
time_limit_sec
```

## Copying

```@docs
ReferenceMap
copy_model
copy_extension_data
Base.copy(::AbstractModel)
```
## I/O

```@docs
write_to_file
Base.write(::IO, ::Model; ::MOI.FileFormats.FileFormat)
read_from_file
Base.read(::IO, ::Type{Model}; ::MOI.FileFormats.FileFormat)
```

## Caching Optimizer

```@docs
MOIU.reset_optimizer(::JuMP.Model)
MOIU.drop_optimizer(::JuMP.Model)
MOIU.attach_optimizer(::JuMP.Model)
```

## Bridge tools

```@docs
bridge_constraints
print_bridge_graph
```

## Extension tools

```@docs
AbstractModel
operator_warn
error_if_direct_mode
```
# [Solutions](@id SolutionAPI)

More information can be found in the [Solutions](@ref jump_solutions) section of
the manual.

## Basic utilities

```@docs
JuMP.optimize!
NoOptimizer
OptimizeNotCalled
solution_summary
```

## Termination status

```@docs
termination_status
raw_status
result_count
```

## Primal solutions

```@docs
primal_status
has_values
value
```

## Dual solutions

```@docs
dual_status
has_duals
dual
shadow_price
reduced_cost
```

## Basic attributes

```@docs
objective_value
objective_bound
dual_objective_value
solve_time
relative_gap
simplex_iterations
barrier_iterations
node_count
```

## [Conflicts](@id ref_conflicts)

```@docs
compute_conflict!
copy_conflict
```

## Sensitivity

```@docs
lp_sensitivity_report
SensitivityReport
```

## Feasibility

```@docs
primal_feasibility_report
```
# [Variables](@id VariableAPI)

More information can be found in the [Variables](@ref jump_variables) section of
the manual.

## Macros

```@docs
@variable
@variables
```

## Basic utilities

```@docs
VariableRef
num_variables
all_variables
owner_model
index(::VariableRef)
optimizer_index(::VariableRef)
check_belongs_to_model
VariableNotOwned
VariableConstrainedOnCreation
VariablesConstrainedOnCreation
```

## Names

```@docs
name(::JuMP.VariableRef)
set_name(::JuMP.VariableRef, ::String)
variable_by_name
```

## Start values

```@docs
set_start_value
start_value
```

## Lower bounds

```@docs
has_lower_bound
lower_bound
set_lower_bound
delete_lower_bound
LowerBoundRef
```

## Upper bounds

```@docs
has_upper_bound
upper_bound
set_upper_bound
delete_upper_bound
UpperBoundRef
```

## Fixed bounds

```@docs
is_fixed
fix_value
fix
unfix
FixRef
```

## Integer variables

```@docs
is_integer
set_integer
unset_integer
IntegerRef
```

## Binary variables

```@docs
is_binary
set_binary
unset_binary
BinaryRef
```

## Integrality utilities

```@docs
relax_integrality
```

## Extensions

```@docs
AbstractVariable
AbstractVariableRef
parse_one_operator_variable
```
# [Containers](@id ContainersAPI)

More information can be found in the [Containers](@ref) section of the manual.

```@docs
Containers
Containers.DenseAxisArray
Containers.SparseAxisArray
Containers.container
Containers.default_container
Containers.@container
Containers.VectorizedProductIterator
Containers.vectorized_product
Containers.NestedIterator
Containers.nested
```

For advanced users, the following functions are provided to aid the writing of
macros that use the container functionality.
```@docs
Containers.build_ref_sets
Containers.container_code
```
# [Nonlinear Modeling](@id NonlinearAPI)

More information can be found in the [Nonlinear Modeling](@ref) section of the
manual.

## [Constraints](@id ref_nl_constraints)

```@docs
@NLconstraint
@NLconstraints
NonlinearConstraintIndex
num_nl_constraints
add_NL_constraint
all_nl_constraints
```

## [Expressions](@id ref_nl_expressions)

```@docs
@NLexpression
@NLexpressions
NonlinearExpression
add_NL_expression
```

## [Objectives](@id ref_nl_objectives)

```@docs
@NLobjective
set_NL_objective
```

## [Parameters](@id ref_nl_parameters)

```@docs
@NLparameter
@NLparameters
NonlinearParameter
value(::JuMP.NonlinearParameter)
set_value(::JuMP.NonlinearParameter, ::Number)
```

## User-defined functions

```@docs
register
```

## Derivatives

```@docs
NLPEvaluator
```
# [Expressions](@id ExpressionAPI)

More information can be found in the [Expressions](@ref) section of the manual.


## Macros

```@docs
@expression
@expressions
```

## Affine expressions

```@docs
GenericAffExpr
AffExpr
linear_terms
```

## Quadratic expressions


```@docs
GenericQuadExpr
QuadExpr
UnorderedPair
quad_terms
```

## Utilities and modifications

```@docs
constant
coefficient
isequal_canonical
add_to_expression!
drop_zeros!
map_coefficients
map_coefficients_inplace!
```

## JuMP-to-MOI converters

```@docs
variable_ref_type
jump_function
jump_function_type
moi_function
moi_function_type
```
# [Constraints](@id ConstraintAPI)

More information can be found in the [Constraints](@ref jump_constraints)
section of the manual.

## Macros

```@docs
@constraint
@constraints
ConstraintRef
AbstractConstraint
ScalarConstraint
VectorConstraint
```

## Names

```@docs
name(::ConstraintRef{Model,<:JuMP._MOICON})
set_name(::ConstraintRef{Model,<:JuMP._MOICON}, ::String)
constraint_by_name
```

## Modification

```@docs
normalized_coefficient
set_normalized_coefficient

normalized_rhs
set_normalized_rhs

add_to_function_constant
```

## Deletion

```@docs
JuMP.delete
is_valid
ConstraintNotOwned
```

## Query constraints

```@docs
list_of_constraint_types
all_constraints
num_constraints
index(::ConstraintRef)
optimizer_index(::ConstraintRef{Model})
constraint_object
```

## Start values

```@docs
set_dual_start_value
dual_start_value
```

## Special sets

```@docs
SecondOrderCone
RotatedSecondOrderCone
PSDCone
SOS1
SOS2
SkewSymmetricMatrixSpace
SkewSymmetricMatrixShape
SymMatrixSpace
moi_set
```

## Printing

```@docs
function_string
constraints_string
in_set_string
show_constraints_summary
```
# Introduction

[Linear programs (LPs)](https://en.wikipedia.org/wiki/Linear_programming) are a
fundamental class of optimization problems of the form:
```math
\begin{align}
    \min_{x \in \mathbb{R}^n} & \sum\limits_{i=1}^n c_i x_i \\
    \;\;\text{s.t.} & l_j \le \sum\limits_{i=1}^n a_{ij} x_i \le u_j & j = 1 \ldots m \\
    & l_i \le x_i \le u_i & i = 1 \ldots n.
\end{align}
```
The most important thing to note is that all terms are of the form
`coefficient * variable`, and that there are no nonlinear terms or
multiplications between variables.

Mixed-integer linear programs (MILPs) are extensions of linear programs in which
some (or all) of the decision variables take discrete values.

## How to choose a solver

Almost all solvers support linear programs; look for "LP" in the list of
[Supported solvers](@ref). However, fewer solvers support mixed-integer linear
programs. Solvers supporting discrete variables start with "(MI)" in the list of
[Supported solvers](@ref).

## How these tutorials are structured

Having a high-level overview of how this part of the documentation is structured
will help you know where to look for certain things.

 * The following tutorials are worked examples that present a problem in words,
   then formulate it in mathematics, and then solve it in JuMP. This usually
   involves some sort of visualization of the solution. Start here if you are
   new to JuMP.
   * [The diet problem](@ref)
   * [The cannery problem](@ref)
   * [The facility location problem](@ref)
   * [Financial modeling problems](@ref)
   * [Network flow problems](@ref)
   * [N-Queens](@ref)
   * [Sudoku](@ref)
 * The [Tips and tricks](@ref linear_tips_and_tricks) tutorial contains a number
   of helpful reformulations and tricks you can use when modeling linear
   programs. Look here if you are stuck trying to formulate a problem as a
   linear program.
 * The [Callbacks](@ref callbacks_tutorial) tutorial explains how to write a
   variety of solver-independent callbacks. Look here if you want to write a
   callback.
 * The remaining tutorials are less verbose and styled in the form of short code
   examples. These tutorials have less explanation, but may contain useful
   code snippets, particularly if they are similar to a problem you are trying
   to solve.
# Introduction

[Nonlinear programs (NLPs)](https://en.wikipedia.org/wiki/Nonlinear_programming)
are a class of optimization problems in which some of the constraints or the
objective function are nonlinear:
```math
\begin{align}
    \min_{x \in \mathbb{R}^n} & f_0(x) \\
    \;\;\text{s.t.} & l_j \le f_j(x) \le u_j & j = 1 \ldots m \\
    & l_i \le x_i \le u_i & i = 1 \ldots n.
\end{align}
```

Mixed-integer nonlinear linear programs (MINLPs) are extensions of nonlinear
programs in which some (or all) of the decision variables take discrete values.

## How to choose a solver

JuMP supports a range of nonlinear solvers; look for "NLP" in the list
of [Supported solvers](@ref). However, very few solvers support mixed-integer
nonlinear linear programs. Solvers supporting discrete variables start with
"(MI)" in the list of [Supported solvers](@ref).

If the only nonlinearities in your model are quadratic terms (that is,
multiplication between two decision variables), you can also use second-order
cone solvers, which are indicated by "SOCP." In most cases, these solvers are
restricted to convex quadratic problems and will error if you pass a nonconvex
quadratic function; however, Gurobi has the ability to solve nonconvex quadratic
terms.

## How these tutorials are structured

Having a high-level overview of how this part of the documentation is structured
will help you know where to look for certain things.

 * The following tutorials are worked examples that present a problem in words,
   then formulate it in mathematics, and then solve it in JuMP. This usually
   involves some sort of visualization of the solution. Start here if you are
   new to JuMP.
   * [Rocket Control](@ref)
   * [Optimal control for a Space Shuttle reentry trajectory](@ref)
   * [Quadratic portfolio optimization](@ref)
 * The [Tips and tricks](@ref nonlinear_tips_and_tricks) tutorial contains a
   number of helpful reformulations and tricks you can use when modeling
   nonlinear programs. Look here if you are stuck trying to formulate a problem
   as a nonlinear program.
 * The [Computing Hessians](@ref) is an advanced tutorial which explains how to
   compute the Hessian of the Lagrangian of a nonlinear program. This is useful
   only in particular cases.
 * The remaining tutorials are less verbose and styled in the form of short code
   examples. These tutorials have less explanation, but may contain useful
   code snippets, particularly if they are similar to a problem you are trying
   to solve.
# Introduction

[Conic programs](https://en.wikipedia.org/wiki/Conic_optimization) are a class
of convex nonlinear optimization problems which use cones to represent the
nonlinearities. They have the form:
```math
\begin{align}
    & \min_{x \in \mathbb{R}^n} & f_0(x) \\
    & \;\;\text{s.t.} & f_j(x) \in \mathcal{S}_j & \;\; j = 1 \ldots m
\end{align}
```

Mixed-integer conic programs (MICPs) are extensions of conic programs in which
some (or all) of the decision variables take discrete values.

## How to choose a solver

JuMP supports a range of conic solvers, although support differs on what types
of cones each solver supports. In the list of [Supported solvers](@ref), "SOCP"
denotes solvers supporting second-order cones and "SDP" denotes solvers
supporting semidefinite cones. In addition, solvers such as SCS and Mosek have
support for the exponential cone. Moreover, due to the bridging system in
MathOptInterface, many of these solvers support a much wider range of exotic
cones than they natively support. Solvers supporting discrete variables start
with "(MI)" in the list of [Supported solvers](@ref).

## How these tutorials are structured

Having a high-level overview of how this part of the documentation is structured
will help you know where to look for certain things.

 * The following tutorials are worked examples that present a problem in words,
   then formulate it in mathematics, and then solve it in JuMP. This usually
   involves some sort of visualization of the solution. Start here if you are
   new to JuMP.
   * [Experiment design](@ref)
   * [Logistic regression](@ref)
 * The [Tips and tricks](@ref conic_tips_and_tricks) tutorial contains a
   number of helpful reformulations and tricks you can use when modeling
   conic programs. Look here if you are stuck trying to formulate a problem
   as a conic program.
 * The remaining tutorials are less verbose and styled in the form of short code
   examples. These tutorials have less explanation, but may contain useful
   code snippets, particularly if they are similar to a problem you are trying
   to solve.
# Introduction

The purpose of these "Getting started" tutorials is to teach new users the
basics of Julia and JuMP.

## How these tutorials are structured

Having a high-level overview of how this part of the documentation is structured
will help you know where to look for certain things.

 * The "Getting started with ..." tutorials are basic introductions to different
   aspects of JuMP and Julia. If you are new to JuMP and Julia, start by reading
   them in the following order:
   * [Getting started with Julia](@ref)
   * [Getting started with JuMP](@ref)
   * [Getting started with sets and indexing](@ref)
   * [Getting started with data and plotting](@ref)
 * Julia has a reputation for being "fast." Unfortunately, it is also easy to
   write _slow_ Julia code. [Performance tips](@ref) contains a number of
   important tips on how to improve the performance of models you write in JuMP.
 * [Design patterns for larger models](@ref) is a more advanced tutorial
   that is aimed at users writing large JuMP models. It's in the "Getting
   started" section to give you an early preview of how JuMP makes it easy to
   structure larger models. If you are new to JuMP you may want to skip
   or briefly skim this tutorial, and come back to it once you have written a
   few JuMP models.
