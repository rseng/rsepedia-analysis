## Future Release (v1.4.1) Date XX.XX.XX

## Release (v1.4.0) Date 22.02.12

* PR #91 migrated to CMake (Issues #85, #64, #26)
* PR #95 added proper coverage behavior to CMake (Issues #94, #93, #92)
* PR #99 remove .gcda and .gcno files on `make scrub` (Issues #84)
* PR #98 Include Geant4 comparison in readthedocs and example-usecases/ (Issues #86, #75)

## Release (v1.3.3) Date 22.02.11

* Hotfix 1.3.3 to add Geant4 comparison plots to master so docs render correctly before next release. 

## Release (v1.3.2) Date 22.01.08 

* PR #77 extended testing to more input files (Issue #35)
* PR #78 updated documentation to systematicall include units of all parameters (Issue #74)
* PR #79 updated documentation to clarify a requirement when installing ROOT dependency on WSL (Issue #76)
* PR #80 improved the testing on the included example notebook (Issue #65)

## Release (v1.3.1) Date 21.12.30 

* PR #62 added fail-on-warning flag for readthedocs (Issue #41)
* PR #68 used breathe to bridge doxygen documentation to readthedocs (Issue #47)
* PR #69 added code coverage and more testing (Issue #49, #44)

## Release (v1.3.0) Date 21.12.27

* PR #54 Move MersenneTwister to built-in C++ function. (Issue #27)
	* The software now requires C++ 11 or newer.
* PR #55 fixed sphinx documentation warnings. (Issue #46)
* PR #56 Travis now checks if files are empty. (Issue #50)
* PR #57 adds testing for c++1--c++17 under ROOT 6.24.08
* PR #58 Improvements to documentation:
	* CONTRIBUTING is no longer cluttered with contents of other documents such as templates and the code of conduct.
	* Syntax for verbosity flags was clarified. (See comments in issue #48)
	* An error in the lifetime units in the documentation was fixed. (All code was already consistent in this regard. Issue #53)
	* CONTRIBUTING now outlines a workflow for changes to be made. (Issue #34)
	* Added arXiv badge. (See issue #37)
	* Include python packages necessary for building docs in the provided yml file.
	* Ignore doc build output and remove build output from repository.
	* Other formatting improvements.
* PR #61 fixes regex to extend to OSX (issue #29)

## Release (v1.2.4) Date 21.12.15

* PR #51 Hotfix to remove the bug introduced by verbosity (solves issue #12)

## Release (v1.2.3) Date 21.11.19

* Hotfix to add submitted tag for JOSS

## Release (v1.2.2) Date 21.11.14

* Hotfix to fix various typos

## Release (v1.2.1) Date 21.11.14

* Hotfix to update release notes. 

## Release (v1.2.0) Date 21.11.13

* PR #38 adds documentation on readthedocs (solves issue #23)
* PR #40 fixes up example (solves issue #7)


## Release (v1.1.3) Date 21.11.06

* PR #16 to lower default verbosity of realizeCascades (solves issue #12; bug introduced see v1.2.4)
* PR #20 adds more clarity to README documentation at top level (solves issue #8)

## Release (v1.1.1) Date 21.10.23

* PR #34 Hotfix to fix Zenodo DOI

## Release (v1.1.0) Date: 21.10.23

* PR #18 and #25 implementing JOSS paper updates (solves issue #11)
* PR #21 correcting the ROOT citation (solves issue #4)
* PR #19 updating the documentation on the mechanism for contributing in CONTRIBUTING.md (solves issue #9)
* PR #28 fix the c++17 bug for the Mersenne Twister in clang on OSX 10.14.6 (issue #27) --> had to re-open this issue PR #28 did NOT solve it
* PR #30 Got a Zenodo DOI (issue #24) 
* PR #17 Got merged accidentally with PR #30 and it partially addressed issue #10 
* PR #31 Got merged (closes out issue #10)
## Template
The following template is not required, but if you do not use it, please be sure to include all answers to all of the questions in some other way.

**Does your pull request resolve or partially resolve an issue?** 
Yes / No.

**If Yes, which issue?** 

**Does your pull request implement code improvements?**
Yes / No.

**Does your pull request implement any breaking changes?**
Yes / No.

**If breaking changes are implemented, please describe:**

**Testing:**  
This pull request:
[ ] Alters the existing CI in some way.
[ ] Adds a new step to the CI.
[ ] Does not introduce any features that the CI could not already test.
[ ] Is not accompanied by necessary CI changes due to some limitation described below. (Please also describe how new features can be manually tested.)

See the `README.md` file in the `nrCascadeSim/tests` directory for more instructions on how to make tests. 
[![Build Status](https://app.travis-ci.com/villano-lab/nrCascadeSim.svg?branch=master)](https://app.travis-ci.com/villano-lab/nrCascadeSim)
[![Documentation Status](https://readthedocs.org/projects/nrcascadesim/badge/?version=latest)](https://nrcascadesim.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![arXiv: 2104.02742](https://img.shields.io/badge/arXiv-2104.02742-orange.svg?style=flat)](https://arxiv.org/abs/2104.02742)
[![codecov](https://codecov.io/gh/villano-lab/nrCascadeSim/branch/master/graph/badge.svg?token=Q6XPU6LPPL)](https://codecov.io/gh/villano-lab/nrCascadeSim)

# nrCascadeSim - a tool for generating nuclear recoil spectra resulting from neutron capture

<!---[![status](https://joss.theoj.org/papers/fd8076268036956d3bf08193c4fc2db9/status.svg)](https://joss.theoj.org/papers/fd8076268036956d3bf08193c4fc2db9)-->
[![status](https://joss.theoj.org/papers/d69ced49c5c17fdbf637e0747d815deb/status.svg)](https://joss.theoj.org/papers/d69ced49c5c17fdbf637e0747d815deb)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5579857.svg)](https://doi.org/10.5281/zenodo.5579857) <br/>

A C/C++ set of executables and library for computing the nuclear recoil spectra left behind by capturing neutrons when all the capture gammas escape the medium. 

<img src="paper/SorVsLin_fig.png" width="500">

You can find more detailed documentation of *nrCascadeSim* [here](https://nrcascadesim.readthedocs.io/en/latest/index.html). The documentation contains e.g. a [guide to get started](https://nrcascadesim.readthedocs.io/en/latest/01_Getting_Started.html).

## CITATION

If you decide to use this code, or if you want to add a reference to it, please cite the latest archived version,

> Villano, A.N., Harris, K., Brown, S. , 2021, nrCascadeSim - A tool for generating nuclear recoil spectra resulting from neutron capture [Code, v1.4.0] [[DOI:10.5281/zenodo.5579857]](https://zenodo.org/record/5579857).

```
@software{nrcascadesim,
  author = {Villano, A.N. and Harris, K. and Brown S.},
  title = {{nrCascadeSim - A tool for generating nuclear recoil spectra resulting from neutron capture [Code, v1.4.0]}},
  year         = {2022},
  publisher    = {Zenodo},
  version      = {v1.4.0},
  doi          = {DOI:10.5281/zenodo.5579857},
  url          = {https://doi.org/10.5281/zenodo.5579857},
  howpublished={The code can be found under \url{https://github.com/villano-lab}.}
}
```

## VERSION HISTORY

- 12.02.2022: Release of [version 1.4.0](https://github.com/villano-lab/nrCascadeSim/releases/tag/v1.4.0)
- 11.02.2022: Release of [version 1.3.3](https://github.com/villano-lab/nrCascadeSim/releases/tag/v1.3.3)
- 08.01.2022: Release of [version 1.3.2](https://github.com/villano-lab/nrCascadeSim/releases/tag/v1.3.2)
- 30.12.2021: Release of [version 1.3.1](https://github.com/villano-lab/nrCascadeSim/releases/tag/v1.3.1)
- 27.12.2021: Release of [version 1.3.0](https://github.com/villano-lab/nrCascadeSim/releases/tag/v1.3.0)
- 15.12.2021: Release of [version 1.2.4](https://github.com/villano-lab/nrCascadeSim/releases/tag/v1.2.4)
- 19.11.2021: Release of [version 1.2.3](https://github.com/villano-lab/nrCascadeSim/releases/tag/v1.2.3)
- 13.11.2021: Release of [version 1.2.0](https://github.com/villano-lab/nrCascadeSim/releases/tag/v1.2.0)
- 06.11.2021: Release of [version 1.1.3](https://github.com/villano-lab/nrCascadeSim/releases/tag/v1.1.3)

## AUTHORS & CONTACT

The authors of *nrCascadeSim* are [A.N. Villano](https://github.com/villaa), [K. Harris](https://github.com/gerudo7), and S. Brown.

For questions, support, bug reports, or other suggestions, please open an [issue](https://github.com/villano-lab/nrCascadeSim/issues).


## LICENSE

This project is licensed under the MIT License - see the LICENSE file.

<!---## Ionization

Ionization assumes the [Lindhard](https://www.osti.gov/biblio/4701226) model:

  Y = k*g(&epsilon;)/(1+kg(&epsilon;))  
  g(&epsilon;) = a\*&epsilon;<sup>&gamma;</sup> + b\*&epsilon;<sup>&omega;</sup> + &epsilon;  
  &epsilon;(E<sub>r</sub>) = 11.5E<sub>r</sub>\[keV\]Z<sup>-7/3</sup>

Using the accepted value for Silicon (*k = 0.143*) or Germanium (*k = 0.159*), whichever is
appropriate; *a = 3*; *b = 0.7*; *&gamma; = 0.15*; and *&omega; = 0.6*.

*Last updated 08 Jan, 2022, v1.2.3*
#Release Checklist

Please be sure to do the following when making a release.

1. Update release page of readthedocs documentation `docs/source/08_Release_History.rst`.
2. Update citation page of readthedocs `docs/source/07_Citations.rst`. (in 3 places--find em' all!)
3. Update the front-facing README in the Version history and citations sections. (in 4 places including version history -- find em' all!)
5. Update release notes `RELEASENOTES.md`.
6. Update `docs/source/index.rst` if badges changed.
7. Restore fail-on-warning on .readthedocs.yaml if it was turned off.
8. Make a release on github.
9. Be sure the `stable` build of readthedocs points to the new release.
10. Be sure to create a version on readthedocs of the new release. 
11. If the version is a patch, deactivate the docs version for the previous patch of the same minor version. (Only one docs version for each minor version should be active at a time.)
12. Be sure codecov website is switched to default to master branch.
13. Update Zenodo.
The following code of conduct must be followed in all places associated with this repository itself. This includes issues, pull requests, and all files contained within the repository. The code of conduct applies to all contributors and community members.

In addition to following GitHub's terms of service and any laws applicable to you, please:

* Do your best to properly attribute any copied, vendored, or otherwise "borrowed" code. In particular, make sure to provide LICENSE files when necessary.
* Be courteous to other contributors/community members.
* Do not make major changes to LICENSE or CONTRIBUTING.md
    * Clarifications and spelling/grammar fixes are fine; modifying the rules is not.
* Avoid vulgar or offensive language.
* Stay on-topic; all discussion, including issues and pull requests, must relate to this repository in some way, and should not focus on something other than the code and its modification.
* Do not add or include malicious code/malware of any kind, including but not limited to ransomware, adware, bloatware, and spyware.

Failure to comply with this code of conduct could result in having your access to the community restricted, such as having offending posts removed or being barred from further submissions.

You can report code of conduct violations to a maintainer (@gerudo7 or @villaa) or the @villano-lab group via direct message or email. If you see a violation of GitHub's ToS or of a local or federal law, please report it to the appropriate authorities first. Thanks!
# Workflow

Below is an outline of our general workflow. It is primarily intended for our software team, 
but outside contributors may wish to draw from it as well.

1. Using gitflow locally, start a branch and make any desired changes.
2. Update the RELEASENOTES.md document and push the branch to GitHub.
3. Create a PR that has a base branch of develop (for features and bugfixes) or master (for releases and hotfixes).
4. If the checks pass and we're happy with the release, merge the PR. **Do not delete the branch on GitHub yet!**
5. On the local copy, run `git flow finish` for the branch you started. (At this point, it is safe to delete the branch on GitHub.)
    * If you merged back to develop, this is the last step. If you merged back to master:
6. Push the local copy of develop.
7. Push tags.
8. Create a release on GitHub from the version tag you just added.
8. Download a tar.gz and update the [Zenodo entry](https://zenodo.org/record/5579858#.YXQsrHlMHOQ) with the new version.

# Standards

All features of `nrCascadeSim` should be tested with Travis-CI wherever possible.  This is done
through the top-level `.travis.yml` file.  You can find Travis-CI's docs
[here](https://docs.travis-ci.com/). Instructions for adding tests can be found in `nrCascadeSim/tests/README.md`.

Addition of features should be accompanied by some form of documentation outlining the use of the
new features. Similarly, if a feature is altered, documentation should be adjusted to match.

Variable and function names should relate to what the variable or function is doing.  
**Unacceptable:** `func` - Name does not describe what the function does at all.  
**Acceptable:** `calculate` - The name indicates a category (that the function calculates something), but it is still vague.  
**Best:** `integrate` - The name clearly indicates what the function does (take an integral).

*(Spot something that's not up to our standards? Submit an issue or pull request!)*

## C/C++ Standards

Keep code neatly organized:
* Leave space before and after multi-line function definitions
* Indent code within braces, one indent per nested brace
* Avoid long lines of code where it's possible to break the code up into multiple lines instead.

Example:

```C++
int n = 0;

int multiply(float a, float b){
    return a*b;
};

void main(int n) {
    if (1 == 1){
        while (n < 3){
            n += 1;
            multiply(1,1);
        }
    }
};
```

## Documentation Standards

All functions accessible to the user should have at least one example provided, possibly more if
the usage is complex or varies significantly.

All variables and options available to the user should be clearly defined.

Documentation files should, at the end of the file, note the date corresponding to the last time
they were updated as well as the relevant version number.


# Pull Requests

Pull Requests (PRs) are created in order to submit to the owner(s) of the repository some code for
consideration. Typically that code will address some issue or improve the code in some way, we
should be clear about how we expect PRs to improve the code in our contributing documentation.
When creating the pull request you have to supply a comparison branch.  When submitting PRs,
please default to submitting to the `develop` branch and never submit directly to the `master`
branch.  This allows us to correctly manage versioning.  If you submit a pull request directly to
the `master` branch, you will be asked to change the target to `develop` (or another applicable
branch).

In order for us to be consistent with GitFlow we should adhere to the following:

**PRs submitted by the development team:** we should locally create a feature branch using the
gitflow protocol on our local machine, and then push that branch. That branch should then be
selected as the "comparison" branch for the PR. Further for the merger to be compatible with
gitflow we should define the base branch as "develop." So the steps are:

1. create a feature branch locally make some changes and push it to the remote (GitHub)
2. open a pull request with base branch `develop` and comparison branch the feature you just created `feature/XXX`
3. When you're done committing (and pushing!) to the feature branch push the button on GitHub to merge the PR back--it will merge it to develop 
4. Delete the branch on github and your local machine and add notes to the upcoming release
5. the feature will be released when the code team does the next release. 

**PR requests submitted from outside our development team:** are very similar to those from the
development team, but the team won't have access to or control over the feature branch created. It
would be created by a fork of the repository. So it looks like this:

1. create a fork of the repository with a branch dedicated to the issue (could be the local `master` we can't enforce any naming conventions there). 
2. open a PR with base branch `develop` and the comparison branch the branch on the fork you just created. 
3. When you're done committing alert the development team in the PR by using the @villaa or other tags. 
4. This will be merged back by the development team if the criteria for code improvement are met. 
5. the feature will be released when the code team does the next release. 

All PRs will be automatically by Travis-CI.  Please note whether you updated the CI or
whether no change was needed.  If for some reason a new, untested feature is implemented, but you
are unable to implement the necessary CI, explain why and how it can be manually tested.

## Release Documentation

When a PR is accepted it will be staged for a release. We will make note of all the currently
staged changes in the RELEASENOTES.md file. It is helpful, but not necessary to put a short
description under the `Next Release` section briefly describing the work done in a PR.  
A template is provided in `pull_request_template.md`.

**Other information:**
Anything else you want to say.

# GitHub Issues

Issues fall into three categories:
* Bug report
* Feature request
* Documentation issue

When submitting issues, please be specific. We can only resolve a bug if we know what about the
program isn't working, implement a feature if we know what aspect is being improved, and clarify
documentation if we know what part is unclear.

Below are outlines for determining what your issue qualifies as. When
submitting an issue, please specify which of these three categories you think it belongs in. We
understand that the three categories can overlap, so don't worry too much if you aren't sure if
the category you chose is appropriate. When creating an issue, you will be given the option to 
select a template; these are just to
help people know what to write, and their use is not strictly required (although it may help us
address the issue faster).

## Bug report

When submitting a bug report, please make sure to include any information necessary for us to
reproduce the bug. If we can't reproduce it, it will be much harder to diagnose and solve the
issue.

An issue is a bug report if:
* The code does not run or only partially runs.
* The code does not build.
* A command associated with the code fails to run despite matching the documentation.
* The code takes an inordinately long amount of time to run, taking hardware into account.
* The code gives no output either to a binary file, a log file, or the terminal, when it should be giving some kind of output there.
* A command that should give consistent output gives different output each time.
* The result of a command directly contradicts what the documentation says should occur.

An issue is not a bug report if:
* The code does not interface with an environment that the documentation does not specify it will interface with. (Feature request)
* The code is missing the ability to do something you think it should be able to do, but the documentation does not specify it is capable of. (Feature request)
* The documentation is unclear but the code does not give results that directly contradict it. (Documentation issue)

## Feature request

An issue is a feature request if:
* You are requesting for the code to interface in a new way that it currently does not, such as a new command or argument.
* You are proposing a particular way to increase the speed of the code.
* You are pointing out where the code could be more user-friendly.
* You are otherwise requesting for the code to do something it is not yet written to do.

An issue is not a feature request if:
* It does not affect the code, only the documentation. (Documentation issue)
* It is to fix unexpected behavior. (Bug report)
* You are providing the feature you are requesting. (Pull request)


## Documentation

An issue is a documentation issue if:
* A command mentioned in the documentation cannot be found.
* You would like an example and there is no similar example in the documentation.
* There is a part of the documentation you are asking for us to clarify.
* There are spelling and grammar errors, missing images, or broken links which you do not know how best to fix with a pull request.

An issue is not a documentation issue if:
* You provide fixed wording, spelling, and/or grammar for all issues you point out. (Pull request)
* The code attempts to run but fails. (Bug report)
* You are looking for a way to do something that you do not know exists and is not mentioned in the documentation. (Feature request)

*Last updated 20 December, 2021, v1.2.4*
# nrCascadeSim - a tool for generating nuclear recoil spectra resulting from neutron capture

This directory contains the object files for the `nrCascadeSim` library. 

## Separate Files

* `libcap.so` -- shared object for the whole `nrCascadeSim` libraries.
* `cascadeProd.o` -- cascade simulating functions
* `edepmath.o` -- mathematical functions
* `isotope_info.o` -- information like neutron separation energies for specific isotopes
* `lindhard.o` -- functions to implement ionization yield
* `rootUtil.o` --  `ROOT` file writing
* `weisskopf.o` -- Weisskopf lifetime estimate functions
The format of the files in this director are:

Each line represents an individual cascade realization for a capture. Each number in the line is defined as follows:

(total fraction): the absolute fraction that this cascade represents in this capture sample

(symbol for post-capture nucleus): like 71Ge represents the final nucleus after the capture

(A value for post-capture nucleus): This is redundant for now

(vector of the energy levels): these are the energy levels AFTER the neutron separation energy.

(vector of level lifetimes): these are the lifetimes in attoseconds of each intermediate level

Here is an example:


0.000671	29Si	29	[4840.4	0.0]	[3500	100000000000000]

Here, this cascade is meant to be simulated with an overall fraction of 0.07%. The capture happens
on 28-Si and proceeds with two emitted gammas. First a gamma is emitted by the difference
between the neutron separation energy and the 4840.4 keV. The next gamma is emitted between the
4840.4 level and the ground state. 

## A Note on Estimation in Levelfiles

If the lifetime of the state was not available, then the Weisskopf estimate for an electric dipole was used.
This generally gives shorter lifetimes than found experimentally.
The dipole is chosen because it is the lowest allowed multipolarity for the transitions calculated
because this leads to underestimation instead of overestimation of the half-lives.
Calculating a lifetime that is too large would eliminate aspects of the spectrum if it prevents decay-in-flight,
while calculating one that is too small will increase the probability of decay-in-flight,
changing the probabilities of some energy readings occurring,
but will not affect what energies can be read in the first place.

[Weisskopf, V. F. Radiative Transition Probabilities in Nuclei. September 1951. Physical Review. Vol83 I5 P1073--1073.](https://link.aps.org/doi/10.1103/PhysRev.83.1073)This directory contains files made for the Geant4 implementation. 
No guarantee is made for the accuracy of these files.---
title: '`nrCascadeSim` - A simulation tool for nuclear recoil cascades resulting from neutron capture'
bibliography: references.bib
tags:
  - C++
  - Simulation
  - Nuclear Physics
authors:
  - name: A.N. Villano
    affiliation: 1
    orcid: 0000-0002-3893-7259
  - name: Kitty Harris
    affiliation: 1
    orcid: 0000-0001-5406-8367
  - name: Staci Brown
    affiliation: 2
affiliations:
 - name: Department of Physics, University of Colorado Denver, Denver CO 80217, USA
   index: 1
 - name: Department of Applied Mathematics & Statistics, University of New Mexico, Albuquerque NM 87131, USA
   index: 2
date: 16 October 2021
#nocite: '@*'
---

# Summary

Neutron capture-induced nuclear recoils have emerged as an important tool for detector
calibrations in direct dark matter detection and coherent elastic neutrino-nucleus scattering
(CE$\mathrm{\nu}$NS).

`nrCascadeSim` is a C++ command-line tool for generating simulation data for energy deposits
resulting from neutron capture on pure materials. Presently, capture events within silicon,
germanium, neon, and argon are supported. While the software was developed for solid state
detector calibration, it can be used for any application which requires simulated neutron
capture-induced nuclear recoil data.

A "cascade" occurs when a neutron becomes part of a nucleus.  The neutron can be captured to one
of many discrete energy levels, or states; if the energy level is nonzero (not the ground state),
then the state will eventually change so that it is zero.  This can happen either all at once or
in multiple steps &mdash; that is, the captured neutron may go from its state to the ground state,
or it may go to another state with lower energy that is not the ground state (provided that one
exists).  The cascade refers to the particular "path" of energy levels that a captured neutron
takes to get to the ground state from the neutron separation energy. Currently, the code assumes
that the neutrons that enter the nuclear system have zero kinetic energy; this is a good
approximation for thermal neutrons because 0.0254\ eV (the average kinetic energy of a thermal
neutron) is small compared to most nuclear recoil energy scales, making it negligible.

`nrCascadeSim` models many of these cascades at once and saves the energies along with other
useful data to a single file. The output file is a `ROOT` file [@ROOT]. 



# Models Used

When modeling deposits from neutron capture events, we want to look at the recoil of the atom as a
result of these cascades.  To determine how much energy is deposited, we must track how much the
atom slows down between steps of the cascade as well as how each nuclear state change affects the
atom's kinetic energy.  `nrCascadeSim` assumes a constant deceleration that results from the atom
colliding with other nearby electrons and nuclei. This means that it must simulate, along with the
steps of the cascade, the time between each state &mdash; to calculate how much the atom slows
down. And it must also simulate the angle between the atom's momentum before a decay and the
momentum boost (gamma ray) resulting from the decay &mdash; to calculate the resulting momenta.
The time between steps is simulated as an exponential random variable based on the state's
half-life, and the angle is simulated as having an isotropic distribution.  Cascade selection is
weighted by isotope abundance [@abundances;@nudat2] and cross-section as well as the probability
of the energy level.  In existing levelfiles, energy levels are derived from [@Ge] for germanium
and from [@Si] for silicon.

The above process models the recoil energies, and the output gives both the total recoil energy
for a cascade as well as the energy per step.  For some applications, this may be the desired
output, or the user may already have a particular process they will use for converting this energy
to what they wish to measure.  However, we also include, for convenience, the ionization yield and
ionization energy of these recoils. Ionization yield is a fraction that, when multiplied by the
energy, gives the ionization energy, and ionization energy is the amount of energy that would be
read out if an otherwise equivalent electron recoil were to occur. This calculation is useful
because many solid-state detectors read out the ionization energy for nuclear recoils. This
ionization yield assumes the Lindhard model [@lindhard].

Figure \ref{LindvSor_fig} compares the normalized frequencies of ionization energies from the
Lindhard [@lindhard] model with the Sorensen [@sorensen] yield model, which is applied after the
simulation using Python, and applies detector resolution models to both. This figure demonstrates
one example of user-applied analysis utilizing the energy deposits at each step instead of the
ionization energy.

![An overlaid histogram showing an example use case in which points are generated and then multiple yield models and resolutions are applied.  The "Small Res (1/5)" histograms have Gaussians with 1/5 of the width of their counterparts. \label{LindvSor_fig}](SorVsLin_fig.pdf)

# Statement of Need

The goal of this software is to simplify the computation of the nuclear recoil spectrum following
neutron capture for a variety of applications.  These include nuclear recoil calibrations for dark
matter direct detection and coherent neutrino detection (CE$\mathrm{\nu}$NS). In these cases as the
particle detection has become more sensitive (detectors having a lower energy threshold) it is now
possible to use the capture-induced nuclear recoil events for detector calibrations. Additionally, 
thermalized neutrons will provide large backgrounds that have heretofore not been modeled. The key
roadblock to studying these scenarios is the complexity of calculating the nuclear recoil
spectrum. 

`nrCascadeSim` addresses this need by allowing users to generate nuclear recoil simulations that
reflect a variety of single-element detector setups. The energy levels that the recoiling nuclei
may pass between and their respective lifetimes are customizable, and multiple isotopes of the
same element can be present within the same simulation. Pre-defined energy level files exist for
silicon and germanium, which take into account the natural abundance data of each isotope in
[@abundances] and [@nudat2].  Output values include energy deposits at each step along each
individual cascade, total kinetic energy deposits, and ionization energy deposits. 


# State of the Field

While there are tools, such as the open-source GEANT4 [@Geant4] framework, that allow users to simulate
neutron capture, existing tools are not built specifically for neutron capture-based nuclear
recoils as `nrCascadeSim` is and therefore use some underlying assumptions that `nrCascadeSim`
does not. The main approximation often used in GEANT4 that we avoid in `nrCascadeSim` is that all
recoils decay directly to the ground state. While this works for some applications, it is
necessary to be more precise when an accurate spectrum of neutron capture-based recoils is needed
for analyses such as calibration or background subtraction. Figure \ref{G4comp} shows a
comparison for the energy deposits produced by Geant4 for natural silicon compared with those
produced by `nrCascadeSim`. The figure does not include any instrumentation resolution and shows a
highly prominent peak around 1.25 keV recoil energy (coming from capture on $^{29}$Si directly to
the ground state) whereas the `nrCascadeSim` shows another direct-to-ground contribution (from
capture on $^{28}$Si) at around 1.0 keV recoil energy and generally far more "spread out" recoils
coming from two- or more step cascades. 

![An overlaid histogram showing how the Geant4 `v10.7.3` energy deposits compare with those from
`nrCascadeSim` for natural silicon. \label{G4comp}](Silicon_comparison.pdf)

Recently, the power of the neutron capture-induced events has been acknowledged in the
CE$\mathrm{\nu}$NS field [@crab]. That initial study, however, used the FIFRELIN code
[@PhysRevC.82.054616], which was originally developed for modeling fission fragments and has been
updated to use statistical models of gamma emission for the purpose of modeling fission-fragment
deexcitation [@FIFRELIN].  `nrCascadeSim` takes the complementary approach of beginning with small
to medium-sized nuclei and modeling the cascades in more exact detail.  The goal is for the code
to be extended to heavier nuclei but still using this detailed approach.     

# Acknowledgments

This material is based upon work supported by the U.S. Department of Energy, Office of Science,
Office of High Energy Physics (HEP) under Award Number DE-SC0021364.

# References
# Tests

This directory houses the testing suite for the code. Each script in this directory is designed to
be run and return 0 on success. When a new test is added it should be added to the Travis-CI
execution list by modifying the file `.travis.yml` in the top-level of this repository.

## `realizeCascades`

The main executable for this code is `realizeCascades` and is simply tested by providing identical
input and comparing the `md5` sum of the `ROOT` output file against the expected value. This
happens inside the script `test_realizeCascades` in this directory.    

# Creating Tests

There are two types of tests that can be created: tests based off of a shell script, and tests
composed of Python scripts. The Python tests use the environment in the top directory
`nrCascadeSim_env.yml`. 

## Script Tests

Place script tests in the `nrCascadeSim/tests` (this) folder. Insert the execution of
the script in the `script:` section of the file `nrCascadeSim/.travis.yml`
after the command `cd tests`.

## Python Tests

Python tests are allowed in the `nrCascadeSim/example-usecase` directory. In
that directory modify the script `test_local.py` to include testing of all of
the Python functions/Jupyter notebooks that are necessary to test the
functionality.  

## C++ Tests

C++ tests are supplied in the `bin` subdirectory 
to test the coverage offered by examples provided on readthedocs.
They can be built from the top-level directory using the command `make tests`;
they are not normal usecases, so they are not made by default.
They should be runnable from this (the `tests`) directory 
(you will need to supply the path to the binaries).
# Example Usecase

The primary purpose of this directory is to house components of `Yields_and_Resolutions.ipynb`, 
which serves as an example of how one might use the outputs of `nrCascadeSim`. 
In this example, the yield of the simulated data is taken, 
a mathematical representation of the detector's resolution is applied,
and the resulting spectra are plotted.

## Instructions

*These instructions are provided in the notebook as well if you would prefer to see them in context.*

The notebook calls a pre-generated file `data/file.root`, which will re-generate the provided image exactly if left unchanged.
If you would like to see the full functionality of nrCascadeSim, however, you can replace this file or generate a new one and alter the path in the library.
You can replace the file by running the script `generate.sh` or by the command below from the repository's top level:
    `bin/realizeCascades -n 10000 -o test-example/data/file.root levelfiles/Si28_ngam_all_cascades.txt`

## Files

* `Yields_and_Resolutions.ipynb` - the notebook serving as the example described above.
* `Z_G4Comparison.ipynb` - a notebook comparing data output with nrCascadeSim to data output with Geant4 in the context of applying yield models.
* `environment.yml` - a conda environment yml file including everything needed to run this example.
* `generate.sh` - a short shell script for running a new simulation using nrCascadeSim. Note that this will overwrite the provided file.
* `rpune.py` - a python script for taking in Geant4 data and removing information that is not necessary for Z_G4Comparison.ipynb (see notebook).
* `requirements.txt` - a pip requirements text file including all pip packages needed to run this example, as an alternative to the yml file.
* `SorVsLin.svg` - a static image of one particular run of the simulation and analysis.
* `standard.mplstyle` - formatting parameters read by the plotting library.
* `test_local.py` - script for CI to ensure that the notebook is running without errors.
* `data` directory - contains some precompiled data for the notebook for those who do not wish to generate their own.
* `plots` directory - contains plots generated by notebooks.
* `python` directory - contains libraries used for the data processing.This directory contains scripts that you may find useful.

## Scripts
`scale.sh` - Multiply [something] by a factor.  

## Example Commands
`combinecommand.txt` - Combine two levelfiles in an organized way without duplicates. (I think - Though I only see one file being called here?)  
(Alt description) - After copying the contents of one levelfile into another, use this command to sort entries and remove duplicates.  
`cascade_finder_commands.txt` - Collection of commands for... something.---
name: Feature request
about: Suggest an idea for improving or expanding this project.
title: ''
labels: enhancement
assignees: villaa

---

I am submitting a **feature request**.

**The feature I am requesting is for:**
make / realizeCascades / CI / specific file / Something new / etc.

**I am requesting:**
A completely new feature / An improvement on an existing feature / etc.

**Ideas for implementation:**
(Optional)

**Other Information:**
Anything else you want to say.
---
name: Documentation Issue
about: Point out issues with documentation.
title: ''
labels: documentation
assignees: gerudo7

---

I am submitting a **documentation issue**.

**The file(s) in question is/are:**
README.md / CONTRIBUTING.md / LICENSE / etc.

**The problem is in the following category/categories:**
Clarity / Examples / Broken links and images / Typos, spelling, and grammar / Undocumented Information / Out-of-date / Other

**Description of the problem:**
Describe whatever is wrong with the documentation or could otherwise be improved.
---
name: Bug report
about: Create a report to help us find and fix bugs.
title: ''
labels: bug
assignees: gerudo7

---

I am submitting a **bug report**.

**This bug occurs in:**
make / realizeCascades / CI / specific file / etc.

**Expected behavior:**
____ should ____.

**Current behavior:**
____ instead does ____.

**Steps to reproduce:**
1. Do thing
2. Do thing
3. Result

**Other Information:**
Anything else you want to say.

**Relevant Output:**
Provide a log file, text from terminal, "No output", etc.
====================
6. Geant4 Comparison
====================

------
Geant4
------

`Geant4  <https://geant4.web.cern.ch/>`_ is a widely popular open-source toolkit for particle physics. 
Unlike nrCascadeSim, Geant4 is a large toolkit with a broad range of applications 
and less focus on any one particular application. 

^^^^^^^^^^
Comparison
^^^^^^^^^^

*Note: You can find the full code for this comparison in:*
`example-usecase/Z_G4Comparison.ipynb <https://github.com/villano-lab/nrCascadeSim/blob/master/example-usecase/Z_G4Comparison.ipynb>`_.

Most notably, Geant4 seems to be much more likely to generate single-step events, 
in which all of the energy from the event is deposited at once, than nrCascadeSim.
Under-generation of events with multiple energy deposits is problematic for modeling detectors with non-linear yield models -
that is, for modeling detectors where the energy read out is not a linear function fo the energy deposited,
making it necessary to take the yielded energy from each individual deposit separately.
Where nrCascadeSim generated single-step events under 10% of the time for silicon, germanium, neon, and argon,
Geant4 generated single-step events over 90% of the time for these elements.

Below are plots comparing the output of nrCascadeSim to the output of Geant4. 
Each plot has data generated by both programs for a single element.
The energies plotted are the total energy deposits for a single event - 
if an event had multiple steps, then the energies for each step were summed before adding them to the plot.
The histograms are weighted to show the portion of events falling in each energy bin.

.. image:: https://raw.githubusercontent.com/villano-lab/nrCascadeSim/master/example-usecase/plots/Silicon_comparison.png
   :width: 500

.. image:: https://raw.githubusercontent.com/villano-lab/nrCascadeSim/master/example-usecase/plots/Germanium_comparison.png
   :width: 500

.. image:: https://raw.githubusercontent.com/villano-lab/nrCascadeSim/master/example-usecase/plots/Neon_comparison.png
   :width: 500

.. image:: https://raw.githubusercontent.com/villano-lab/nrCascadeSim/master/example-usecase/plots/Argon_comparison.png
   :width: 500
====================================
4. Examples: Putting it all together
====================================

These examples are in the `tests/bin` directory and can be compiled using `make tests` from teh top level.

---------------------------
Reading a Level Input File
---------------------------

This example simply reads a level input file called `inputfile.txt` of the correct format. It then
prints out the contents of that file in a format showing the properties of each cascade to
standard output. 


.. code-block:: C 

  //library commands
  #include "cascadeProd.h"
  #include "lindhard.h"
  #include "weisskopf.h"
  #include "isotope_info.h"
  
  //ROOT stuff
  #include "rootUtil.h"

  //get the file name
  string filename="inputfile.txt";

  //read the contents of a file into a cli object
  int numc;
  bool success=false;
  cli *cascadeFile = readCascadeDistributionFile(numc,filenames[i],success);

  //print the info that was read in
  for(int i=0;i<numc;i++){
      cout << "Cascade ID: " << i+1 << "/" << numc << endl;
      cout << "Fraction of this cascade: " << cascadeFile[i].frac << endl;
      cout << "Neutron separation: " << cascadeFile[i].Sn << endl;
      cout << "Mass number: " << cascadeFile[i].A << endl;
      cout << "Number of steps: " << cascadeFile[i].n << endl;
      cout << endl;
      cout << "Energy Levels (keV)\t|\ttau (fs)" << endl;
      cout << "------------------------------------------------" << endl;
      cout << setfill('0') << setw(5) << setprecision(5);
      cout << "      *****       " << "\t \t" << " ***** " << endl;
      for(int j=0;j<cascadeFile[i].n;j++){
        cout << "      "<< cascadeFile[i].Elev[j] << "       " << "\t \t" << " "<< cascadeFile[i].taus[j] << " " << endl;
   }

   return 0;
          

----------------------------
Printing Isotope Information
----------------------------

This example prints all of the information inside of the isotope database to standard output.

.. code-block:: C 

  //library commands
  #include "cascadeProd.h"
  #include "lindhard.h"
  #include "weisskopf.h"
  #include "isotope_info.h"

  //ROOT stuff
  #include "rootUtil.h"

  //print out all isotope information in database
  int main(){

      listStuff();

      return 0;

  }

----------------------------------
Fetching an Ionization Yield Model
----------------------------------

This example computes the expected ionization based on the standard Lindhard model
[Lindhard1963]_. The computation is done for a slowing germanium ion between 100 eV and 50 eV. The
result is printed to standard output.

.. code-block:: C 

  //library commands
  #include "cascadeProd.h"
  #include "lindhard.h"
  #include "weisskopf.h"
  #include "isotope_info.h"
  #include <iostream>
  //ROOT stuff
  #include "rootUtil.h"

  //seed an MT random number with 1
  std::mt19937 *mtrand = new std::mt19937(1);

  //get the ionization deposited between 100 and 50 eV
  // double *ionization;
  double E0 = 100; //eV
  double E1 = 50; //eV
  double *ionization = geIonizationInRange_k(E0,E1,0.159,mtrand); //k-value for Germanium (accepted)

  //print the ionization
  int main(){
      std::cout << "Ionization Energy: " << ionization[0] << " eV " << endl;
      std::cout << "Ionization Pairs: " << ionization[1] << " eV " << endl;
  return 0;
  }

--------------------------------------------------
Generating a Single Cascade Realization and Saving
--------------------------------------------------

This example reads in a cascade input file called `inputfile.txt`, realizes approximately 10,000
total cascades in the distribution governed by the input file, and then ports the outputs to a ROOT file
named `output.root`. Each cascade that is realized is printed to standard output showing how many
events were realized for each particular type of cascade. 

.. code-block:: C 

  //library commands
  #include "cascadeProd.h"
  #include "lindhard.h"
  #include "weisskopf.h"
  #include "isotope_info.h"
  #include <iostream>
  //ROOT stuff
  #include "rootUtil.h"

  //get the file name
  string filenames="inputfile.txt";

  //read the contents of a file into a cli object
  int numc;
  bool success=false;
  cli *cascadeFile = readCascadeDistributionFile(numc,filenames,success);

  //get a root file and make
  TFile *f = TFile::Open("output.root","recreate");
  TTree *t = new TTree("cascade","cascade");

  //random number
  std::mt19937 *mtrand = new std::mt19937(1);

  //calculate the cascades
    int main(){
    int num = 10000;
    for(int k=0;k<numc;k++){
      int nrealize = num*cascadeFile[k].frac;
      std::cout << "Realizing " << nrealize << " events of cascade ID " << cascadeFile[k].cid << endl;
      cri *cascade_data;
      cascade_data = Cascade(nrealize,cascadeFile[k].cid,cascadeFile[k].Sn,cascadeFile[k].n,cascadeFile[k].Elev,cascadeFile[k].taus,cascadeFile[k].A,mtrand);
      std::cout << "Cascade realization " << k << " success: " << addToNRTTree(t,nrealize,cascade_data,cascadeFile[k]) << endl;

      freecriarray(nrealize,cascade_data);
    }
    freecliarray(numc,cascadeFile);

    //write the ROOT file
    t->Write("",TObject::kOverwrite);
    f->Close();

    return 0;
  }
    
========================================
2. Executables of *nrCascadeSim*
========================================

Using the library functions defined in *nrCascadeSim* through the shared-object library
`libncap.so` users can define their own executables for a task or use the pre-defined executables.


------------------------------------------------
Using *realizeCascades* command-line executable 
------------------------------------------------

To get a quick list of the expected inputs and flags one can always run `realizeCascades --help`.
The result will be:

.. code-block:: bash 

  Usage:  (null) options [ inputfile(s) ]
    -d, --seed          <integer>      Seed for random numbers
    -h, --help                         Print usage 
    -n, --numgen        <number>       Number of traces to generate
    -o, --outfile       <filename>     Name the output file
    -s, --silent                       Silent, no standard out
    -v, --verbose       <level>        Print verbose messages at level <level>.
                                       Currently must use `--verbose=<level>` or `-v<level>` - no spaces.
    -V, --version                      Print version and exit
    -l, --log           <filename>     Log additional output to the specified file.
                                       If this option is not used, no logging will occur.


The `realizeCascades` command will run the simulation a specified number of times for a given
input file.  

Note that `realizeCascades` must be run from the `nrCascadeSim/bin` directory unless the user has
added it to the path like by doing `sudo make install`.  Also note that `ROOT` must be present in
the current environment for this command to work.

^^^^^^^^^
Arguments
^^^^^^^^^

""""""""""""""""""
Required Arguments
""""""""""""""""""

All three of these arguments are required:
* ``-n, --numgen`` specifies the total number of cascade events to be simulated. (example: `-n 100000` to simulate one hundred thousand events.)
* ``-o, --outfile`` specifies the location of the output file. (example: `-o ~/output.root` to output to a file `output.root` in the home directory.)
* The main argument (no prefix) specifies the input file. (example: `levelfiles/Si28_ngam_all_cascades_rfmt_sorted.txt` to call a levelfile with all cascades for 28Si available.)

This makes the full example:

.. code-block:: bash 

  realizeCascades -n 100000 -o ~/output.root levelfiles/Si28_ngam_all_cascades_rfmt_sorted.txt


to simulate 100000 events for 28Si and output them to a file in the home directory.

""""""""""""""""""
Optional Arguments
""""""""""""""""""

* ``-h, --help`` display the help menu. This overrides other options. Help will be displayed and program will exit. 
* ``-s, --silent`` silent. Nothing will be printed to stdout.
* ``-v, --verbose`` verbosity level. Default to 1 where only the random seed is printed. Max level is currently 2 where a lot of level/simulation information is printed. Currently must use ``--verbose=<level>`` or ``-v<level>`` - no spaces.
* ``-V, --version`` version. Print the version tag of the code and exit.  
* ``-l, --log`` log file. Specify a file to print the output to.  


^^^^^^^^^^^^^^^^^^
Reproducible Files
^^^^^^^^^^^^^^^^^^

The ``-d, --seed`` flag for the seed will result in files with consistent data. 

However, additional binary data may result in checksums being different despite the data being the
same.  If you want a reproducible file that can be compared to another by an md5 checksum, append
to the output file's name:  `?reproducible=fixedname`  (you will either need to put the filename
in quotes or escape the `?` character).  Example: `-o "output.root?reproducible=fixedname"`  This
surpresses various forms of metadata that result in changes to the binary even for the same data
(see `ROOT` page on reproducible-files_ ).

.. _reproducible-files: https://root.cern.ch/doc/master/classTFile.html#ad0377adf2f3d88da1a1f77256a140d60 

Unfortunately, the binary still seems to be influenced by the environment it is generated in,
so at present there is no md5sum to compare to that will work across all devices.

-------------------------------------
*realizeCascades* cascade input file 
-------------------------------------

^^^^^^^^^^^^^^^^^^^^^^^^
Levelfile (Input) Format
^^^^^^^^^^^^^^^^^^^^^^^^

The levelfile is a singular text file read by the program using regular expressions.  Each row in
a levelfile corresponds to one possible cascade, which should include a relative weight for the
probability of the cascade's occurrence.  While it is helpful to create columns that are easy for
the user to read, columns can be delineated by any number of spaces.

The general format of one row of an input file is:

.. code-block:: bash
   
   weight isotope A [..,E2,E1,E0] [..,tau2,tau1,inf]

Each portion of this row is described in the table below.

.. list-table:: Input Row 
   :widths: 25 25 50
   :header-rows: 1

   * - Name
     - Format
     - Description
   * - `weight`
     - `#.##` or `#e+/-##`
     - The probability of this cascade occuring, normalized to unity with all other cascades. This variable includes a weight for the isotope's abundance compared to other isotopes listed within the same levelfile. (If only one isotope is present within the levelfile, the abundance weight is not needed.) Weights can be given in decimal form or scientific notation (e.g. 0.000671 or 6.71e-04). 
   * - `isotope`
     - `##Xx`
     - The isotope of the nucleus *after* capture. (For example, if 28Si is present, it will become 29Si, so 29Si should be listed.) This should be formatted as two numbers, one capital letter, and one lower-case letter (e.g. 29Si, 74Ge).
   * - `A`
     - `##`
     - The number of particles in the nucleus after capture. This should match the first two digits of `isotope`. For example, if `isotope` is 72Ge, `A` should be 72.
   * - `energies` 
     - `[... E2 E1 0]`
     - An ordered list of the energy levels traversed (keV), including the ground state (0 keV), separated by spaces. These should be in the decreasing order, the order in which the nucleus will go through the states. Do not include the separation energy to account for the initial unbound state before capture; this is already assumed.
   * - `lifetimes`
     - `[... tau2 tau1 inf]`
     - An ordered list of the lifetimes of the energy levels traversed (fs), separated by spaces. It must be the same length as the list of energies, and the lifetimes should be in the same order as the energies. The last entry is `100000000000000.0` (1e+14 fs, or 100 ms), which is effectively infinite on the timescale of the simulation, to indicate that the state is stable at the ground state.

Note: for the lifetimes one can also use strings representing multipolarity (like `w(M1)`,
`w(E1)`, etc. for magnetic or electric dipole transitions respectively) to instruct the program to
use the corresponding Weisskopf estimate [Weisskopf1951]_. These estimates are not very accurate
and are known to be systematically low. 

^^^^^^^^^^
On Weights
^^^^^^^^^^

The sum of the probabilities must be less than or equal to one in order for the simulation to 
work properly. If the sum is less than one, the simulation may skip generating some points in 
the output &mdash; for example, when requesting 100 entries, if the total probability is 0.95, 
one would expect 95 entries on average &mdash; but the input cascades will still be at the 
correct proportions with respect to one another. If the sum is greater than one, the simulation 
may not reach certain cascades at all--for instance, if a file has 12 cascades, and the 
probabilities of the first 10 add up to 1, then the last two will never be generated.

"""""""""""""""""""""""""""""""""""
An example for calculating weights:
"""""""""""""""""""""""""""""""""""

A silicon detector has three isotopes, which become 29Si, 30Si, and 31Si after capture.  The
abundances within the detector are 60%, 30%, and 10%, respectively.  Each has three possible
cascades we want to model, which we list below in our (incomplete) draft of the levelfile:

.. code-block:: bash

   weight? 29Si 29 [0]         [100000000000000.0]
   weight? 29Si 29 [5000 0]    [0.84   100000000000000.0]
   weight? 29Si 29 [3000 0]    [0.5    100000000000000.0]
   weight? 30Si 30 [0]         [100000000000000.0]
   weight? 30Si 30 [4000 0]    [1      100000000000000.0]
   weight? 30Si 30 [2000 0]    [0.15   100000000000000.0]
   weight? 31Si 31 [0]         [100000000000000.0]
   weight? 31Si 31 [4999 0]    [0.15   100000000000000.0]
   weight? 31Si 31 [540  0]    [.954   100000000000000.0]


Let's say the probabilities of the cascade occurring **within the respective isotopes** are as below:


+-----------+-------------+-------+---------------+---------------+
| **29Si**  | Cascade     | `[0]` | `[5000    0]` | `[3000    0]` | 
+-----------+-------------+-------+---------------+---------------+
|           | Probability | 0.35  |    0.5        |    0.15       |
+-----------+-------------+-------+---------------+---------------+

+-----------+-------------+-------+---------------+---------------+
| **30Si**  | Cascade     | `[0]` | `[4000    0]` | `[2000    0]` | 
+-----------+-------------+-------+---------------+---------------+
|           | Probability | 0.8   |    0.1        |    0.1        |
+-----------+-------------+-------+---------------+---------------+

+-----------+-------------+-------+---------------+---------------+
| **31Si**  | Cascade     | `[0]` | `[4999    0]` | `[540    0]`  | 
+-----------+-------------+-------+---------------+---------------+
|           | Probability | 0.2   |    0.3        |    0.5        |
+-----------+-------------+-------+---------------+---------------+


Then the relative probabilities **within the simulation** are:

+-----------+-------------+-------+---------------+---------------+
| **29Si**  | Cascade     | `[0]` | `[5000    0]` | `[3000    0]` | 
+-----------+-------------+-------+---------------+---------------+
|           | Probability | 0.21  |    0.3        |    0.09       |
+-----------+-------------+-------+---------------+---------------+

+-----------+-------------+-------+---------------+---------------+
| **30Si**  | Cascade     | `[0]` | `[4000    0]` | `[2000    0]` | 
+-----------+-------------+-------+---------------+---------------+
|           | Probability | 0.24  |    0.03       |    0.03       |
+-----------+-------------+-------+---------------+---------------+

+-----------+-------------+-------+---------------+---------------+
| **31Si**  | Cascade     | `[0]` | `[4999    0]` | `[540    0]`  | 
+-----------+-------------+-------+---------------+---------------+
|           | Probability | 0.0   |    0.03       |    0.05       |
+-----------+-------------+-------+---------------+---------------+


Making our completed levelfile:

.. code-block:: bash

   0.21    29Si 29 [0]         [100000000000000.0]
   0.30    29Si 29 [5000 0]    [0.84   100000000000000.0]
   0.09    29Si 29 [3000 0]    [0.5    100000000000000.0]
   0.24    30Si 30 [0]         [100000000000000.0]
   0.03    30Si 30 [4000 0]    [1      100000000000000.0]
   0.03    30Si 30 [2000 0]    [0.15   100000000000000.0]
   0.02    31Si 31 [0]         [100000000000000.0]
   0.03    31Si 31 [4999 0]    [0.15   100000000000000.0]
   0.05    31Si 31 [540  0]    [.954   100000000000000.0]

^^^^^^^^^^^^^^^^^^^^^^^^^
On Energies and Lifetimes
^^^^^^^^^^^^^^^^^^^^^^^^^

In the following levelfile row, the nth lifetime entry corresponds to the nth energy level entry.

.. code-block:: bash

   0.30    29Si 29 [5000 4000 3000 2000 1000 0]    [0.84 0.95 1.35 0.03 0.11 100000000000000.0]

Therefore, the program reads this as:

+---------------+-----------+-----------+-----------+-----------+-----------+
| Energy level: | 5000 keV  | 4000 keV  | 3000 keV  | 2000 keV  | 1000 keV  |
+---------------+-----------+-----------+-----------+-----------+-----------+
| **Lifetime:** | 0.84 fs   | 0.95 fs   | 1.35 fs   | 0.03 fs   | 0.11 fs   |
+---------------+-----------+-----------+-----------+-----------+-----------+


-------------------------------------
*realizeCascades* cascade output file 
-------------------------------------

Note: ROOT_ is needed to open these files.

.. _ROOT: https://root.cern/install/

A file that contains the separate NR deposits, along with their Ionization deposits.  In addition
all of the exiting gamma energies and times should be listed.

The output files are `*.root` files and therefore cannot be read as text.
Instead, they need to be imported to a program to be read out.
One straightforward way of reading these files is with python and the
`uproot <https://pypi.org/project/uproot/>`_ package.

The `*.root` files store information in a tree-like structure. The top-most key in the output
files will be `cascade` (there are no other top-level keys). Beneath this, the following keys
exist:  

.. list-table:: Output Structure 
   :widths: 25 25 25 50
   :header-rows: 1
   
   * -  `Name`  
     -  `Shape`       
     -  **Units** 
     -  Description 
   * -  `n`  
     -  `1D Array`    
     -   N/A       
     -  Array denoting the number of energy levels in a given cascade. This includes intermediate levels and the ground state.
   * -  `cid`    
     -  `1D Array`    
     -   N/A       
     -   Array of cascade IDs. The cascade ID is the number of the row in the levelfile which contains the cascade used. These count starting from zero.
   * -  `Elev` 
     -  `Jagged Array` 
     -   **keV**   
     - Array of energy level inputs. Each entry is an array of size `n`.
   * - `taus`
     - `Jagged Array`
     -  **femto-sec (fs)**  
     -  Array of lifetime inputs. Each entry is an array of size `n`.
   * -  `delE`   
     -  `Jagged Array`
     -   **eV**    
     -   Array of energy deposits between energy levels. Each entry is an array of size `n - 1`. It contains the individual energy deposits, not the total energy deposit. If using a custom nonlinear ionization model, these are the best to operate on.
   * -  `I` 
     - `Jagged Array`
     -  None   
     -  Array containing the ionization calculations for each energy deposit. Each entry is an array of size `n - 1`. This ionization is given in terms of a number of charges.
   * -  `Ei` 
     -  `Jagged Array`
     -   **eV** 
     -  Array of calculated ionization energy per step. These energies are conversions of `delE` to ionization energies. Each entry is an array of size `n - 1` containing the individual ionization energies. The Lindhard model is used here.
   * -  `time`  
     -  `Jagged Array`
     -   **fs**  
     -   Array of the time spent at each energy level. Each entry is an array of size `n` containing individual times.
   * -   `Eg`  
     -   `Jagged Array`
     -    **MeV**   
     -    Array of gamma energies. Each entry is an array of gamma energies, corresponding to an energy deposit.

The ordering of values in the arrays are consistent; that is, the nth entry of `n` corresponds to
the nth entry of `cid`, the nth entry of `Elev`, and so on.  The length of each main array should
be equal to the number of simulations; that is, if running 10000 events, `n` and `cid` will have
lengths of 10000 and the jagged arrays will have first dimensions of length 10000.

.. image:: https://raw.githubusercontent.com/villano-lab/nrCascadeSim/master/output_structure.svg 
   :width: 750 
.. The three most important abstract base classes of *obscura* are

.. #. ``DM_Particle``
.. #. ``DM_Distribution``
.. #. ``DM_Detector``

.. We will discuss the interface each of these classes provide in more detail.
.. But first we take a look at the detection targets in direct DM search experiments, namely nuclei, bound electrons in atoms, and bound electrons in crystals.
=====================
Citing *nrCascadeSim*
=====================

-----------
How to cite
-----------

If you decide to use this code, or if you want to add a reference to it, please cite the latest archived version,

    Villano, A.N., Harris, K., Brown, S., 2021, nrCascadeSim - A tool for generating nuclear recoil spectra resulting from neutron capture [Code, v1.4.0] [DOI:10.5281/zenodo.5579857]

.. raw:: html

	<details>
	<summary><a>Bibtex entry.</a></summary>
 
.. code-block::

    @software{nrcascadesim,
    author = {Villano, A.N. and Harris, K. and Brown S.},
    title = {{nrCascadeSim - A tool for generating nuclear recoil spectra resulting from neutron capture [Code, v1.4.0]}},
    year         = {2022},
    publisher    = {Zenodo},
    version      = {v1.4.0},
    doi          = {DOI:10.5281/zenodo.5579857},
    url          = {https://doi.org/10.5281/zenodo.5579857},
    howpublished={The code can be found under \url{https://github.com/villano-lab}.}
    }

.. raw:: html

	</details>


---------------------------------------------------
Research and research software using *nrCascadeSim*
---------------------------------------------------

The library *nrCascadeSim* has been applied to obtain the scientific results of the following papers

#. **First observation of isolated nuclear recoils following neutron capture**
  
  A.N. Villano, M. Fritts, N. Mast, S. Brown, P. Cushman, K. Harris, V. Mandic 

  .. image:: https://img.shields.io/badge/arXiv-2110.02751-B31B1B.svg
      :target: https://arxiv.org/abs/2110.02751
      :alt: [arXiv:2110.02751]


.. Here is a list of research software using *nrCascadeSim*:

.. #. Emken, T., 2021, `Dark Matter Simulation Code for Underground Scatterings - Sun Edition (DaMaSCUS-SUN) <https://github.com/temken/DaMaSCUS-SUN>`_ Astrophysics Source Code Library, record `[ascl:2102.018] <https://ascl.net/2102.018>`_, `[DOI:10.5281/zenodo.4559874] <https://zenodo.org/record/4559874>`_

.. .. image:: https://github.com/temken/obscura/actions/workflows/main.yml/badge.svg?branch=master
..   :target: https://github.com/temken/obscura/actions/workflows/main.yml
..   :alt: Build Status
.. .. image:: https://codecov.io/gh/temken/obscura/branch/master/graph/badge.svg?token=1Pe1QMcngr
..   :target: https://codecov.io/gh/temken/obscura
..   :alt: Code Coverage 
.. image:: https://app.travis-ci.com/villano-lab/nrCascadeSim.svg?branch=master 
   :target: https://app.travis-ci.com/villano-lab/nrCascadeSim
   :alt: Build Status 
.. image:: https://readthedocs.org/projects/nrcascadesim/badge/?version=latest
   :target: https://nrcascadesim.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License: MIT
.. image:: https://img.shields.io/badge/arXiv-2104.02742-orange.svg?style=flat
   :target: https://arxiv.org/abs/2104.02742
   :alt: arXiv: 2104.02742
.. image:: https://codecov.io/gh/villano-lab/nrCascadeSim/branch/master/graph/badge.svg?token=Q6XPU6LPPL
   :target: https://codecov.io/gh/villano-lab/nrCascadeSim
    

============================================================================================
*nrCascadeSim* - a tool for generating nuclear recoil spectra resulting from neutron capture
============================================================================================

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5579857.svg
   :target: https://doi.org/10.5281/zenodo.5579857
   :alt: DOI
.. image:: https://joss.theoj.org/papers/d69ced49c5c17fdbf637e0747d815deb/status.svg
   :target: https://joss.theoj.org/papers/d69ced49c5c17fdbf637e0747d815deb
   :alt: JOSS paper

The purpose of this code is to simulate energy deposits due to cascading of energy levels
following neutron capture.  This code was written for use in nuclear recoil calibration for dark
matter detectors, but may be useful in other particle physics applications as well, including
coherent elastic neutrino nucleus scattering (CEνNS).  Currently, we use a constant
acceleration model for the atom slowing down and calculation of the ionization energy.  We also
use the Lindhard model for calculating the ionization, but the output is complete enough to allow
the user to choose their ionization yield model after simulation.  The code currently supports
Neon, Argon, Silicon, and Germanium cascades slowing down in a lattice of like material.

The documentation does not contain a review of the physics implemented in the library.

If you want to contribute to `nrCascadeSim`, please check out the `contribution guidelines
<https://github.com/villano-lab/nrCascadeSim/blob/master/CONTRIBUTING.md>`_.

.. image:: https://raw.githubusercontent.com/villano-lab/nrCascadeSim/master/paper/SorVsLin_fig.png 
   :width: 500

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   01_Getting_Started
   02_Executables
   03_Doxy_Libraries
   04_Examples
   05_Documentation
   06_Geant4_Comparison
   07_Citations
   08_Release_History
   09_License
   10_Contact
   References

===============
Release history
===============

* 12.02.2022: Release of `v1.4.0 <https://github.com/villano-lab/nrCascadeSim/releases/tag/v1.4.0>`_
* 11.02.2022: Release of `v1.3.3 <https://github.com/villano-lab/nrCascadeSim/releases/tag/v1.3.3>`_
* 08.01.2022: Release of `v1.3.2 <https://github.com/villano-lab/nrCascadeSim/releases/tag/v1.3.2>`_
* 30.12.2021: Release of `v1.3.1 <https://github.com/villano-lab/nrCascadeSim/releases/tag/v1.3.1>`_
* 27.12.2021: Release of `v1.3.0 <https://github.com/villano-lab/nrCascadeSim/releases/tag/v1.3.0>`_
* 15.12.2021: Release of `v1.2.4 <https://github.com/villano-lab/nrCascadeSim/releases/tag/v1.2.4>`_
* 19.11.2021: Release of `v1.2.3 <https://github.com/villano-lab/nrCascadeSim/releases/tag/v1.2.3>`_
* 13.11.2021: Release of `v1.2.0 <https://github.com/villano-lab/nrCascadeSim/releases/tag/v1.2.0>`_
* 06.11.2021: Release of `v1.1.3 <https://github.com/villano-lab/nrCascadeSim/releases/tag/v1.1.3>`_
==================
1. Getting started
==================

------------
Installation
------------

*nrCascadeSim* is designed to run in a Unix-based system and is tested via Travis-CI_ using the
Xenial_ (Ubuntu 16.04), Bionic_ (Ubuntu 18.04), and Focal_ (Ubuntu 20.04) distributions, all on
x86_64 architecture. It is also tested on Mac OSX Mojave (10.14), Catalina (10.15), and Big Sur
(11.5) via Travis-CI_. 

Because of a combination of the operating systems tested and the versions of `ROOT` (see below)
the Travis-CI_ tests three C++ compiler standards: `c++11` (Xenial_, Mojave, Catalina, Big Sur);
`c++14` (Bionic_); and `c++17` (Focal_). 

.. _Travis-CI: https://app.travis-ci.com/github/villano-lab/nrCascadeSim
.. _Xenial: https://docs.travis-ci.com/user/reference/xenial/ 
.. _Bionic: https://docs.travis-ci.com/user/reference/bionic/ 
.. _Focal:  https://docs.travis-ci.com/user/reference/focal/

^^^^^^^^^^^^
Dependencies
^^^^^^^^^^^^

""""""""""""""""""""""""""""""""""""
1. `ROOT <https://root.cern/>`_
""""""""""""""""""""""""""""""""""""

Travis-CI_ tests two versions of `ROOT`: 6.24.06 and 6.20.00. It is assumed (though not explicitly
tested) that the versions in between those two versions will also work. 

To install `ROOT` please follow the instructions on the `CERN  <https://root.cern/install/>`_
website.

It is intended that *nrCascadeSim* is compatible with all versions; known compatibility with
version 6. 

On Linux machines, you can also install `ROOT` via a `pre-packaged binary
<https://root.cern/install/#download-a-pre-compiled-binary-distribution>`_ run::

	wget https://root.cern/download/root_v6.24.02.Linux-ubuntu20-x86_64-gcc9.3.tar.gz
	tar -xzvf root_v6.24.02.Linux-ubuntu20-x86_64-gcc9.3.tar.gz
	source root/bin/thisroot.sh # also available: thisroot.{csh,fish,bat}

**Note:** You will need to find the appropriate ROOT link for your Linux distribution. 
The one provided above is used for Focal. Bionic uses ``root_v6.24.06.Linux-ubuntu18-x86_64-gcc7.5.tar.gz``
and Xenial uses ``root_v6.24.06.Linux-ubuntu20-x86_64-gcc9.3.tar.gz``.

**If you are using WSL,** you will need to `install the conda package <https://root.cern/install#conda>`_ (recommended) 
or `build ROOT from source <https://root.cern/install#build-from-source>`_. 
Using a pre-compiled binary or installing from most package managers will not work for WSL users.

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
2. `gcc <https://gcc.gnu.org/>`_
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. To install *gcc* on a Mac, we can use `homebrew <https://brew.sh/>`_ ::

..	brew install libconfig

You almost certainly have a version of `gcc` already and *nrCascadeSim* is compatible with version
4.4.7 or newer.

On Linux machines, you can build `gcc` via the `apt-get` manager::

	sudo apt-get install gcc

If, for some reason, you need to use a C++ standard older than 11, please use nrCascadeSim v1.2.2 or older.


^^^^^^^^^^^^^^^^
Download & Build
^^^^^^^^^^^^^^^^

The `nrCascadeSim` source code can be downloaded by cloning this `git repository
<https://github.com/villano-lab/nrCascadeSim>`_: ::

   git clone https://github.com/villano-lab/nrCascadeSim.git
   cd nrCascadeSim 

The code is compiled and the executable and library is built by `CMake <https://cmake.org/>`_. To build, choose a build directory `/path/to/build/` and do.::

        mkdir -p /path/to/build
        cd /path/to/build
	cmake -DCMAKE_INSTALL_PREFIX=/path/to/install /path/to/source	
	make

If everything worked well, the executable and library file are created in the build directory as::

	realizeCascades
	regexPlayground
	lib/libncap.so	

And various `.o` files. If you would like to build the testing suite, use::

        make tests

To clean the installation use::

	make clean

To install to the directory `/path/to/install` use::

        make install

From the build directory. Then, if you want to uninstall use::

        make uninstall
       



---------------------------------------
Using *nrCascadeSim* command-line tools
---------------------------------------

*nrCascadeSim* contains pre-built executables built on the library functions that can be used as
command-line tools. See Section 2 for list of the executables and inputs and outputs. 

.. .. warning::

..	The import of these parameters via libconfig is very case-sensitive. A float parameter has to be set to e.g. *1.0*, and **not** just *1*.

..  .. raw:: html

..	<details>
.. 	<summary><a>The full configuration file</a></summary>
 
 
.. .. raw:: html

..	</details>

----------------------------------
Using *nrCascadeSim* as a library
----------------------------------

If we want to use *nrCascadeSim* functions in an external code, we can do so and import it as a library.
We recommend to do this inside your Make build, where the *nrCascadeSim* libraries can be linked
with the `-fPIC` and `-lncap` flags.

Use the following headers when including the library functions:


.. code-block:: c++

  //library commands
  #include "cascadeProd.h"
  #include "lindhard.h"
  #include "weisskopf.h"
  #include "isotope_info.h"
  #include "MersenneTwister.h"
  
  //ROOT stuff
  #include "rootUtil.h"


Using the `ROOT` utilities  will also require having `ROOT` installed and adding `root-config
--cflags --glibs` to the compilation line in your Make file. 


.. As an instructional example `this repository <https://github.com/temken/template_cpp_cmake_obscura>`_ contains a C++ project template built with CMake that imports and uses the *obscura* library.
========================================
3. The libraries of *nrCascadeSim*
========================================

The functionality of the executables of *nrCascadeSim* such as `realizeCascades` is constructed
from a set of library functions whose code is stored in the `libncap.so` file and whose prototypes
are stored in several header files. Below we list the different header files with the internal
functions, data structures, and their uses. 

---------------
`cascadeProd.h`
---------------

This header has the prototypes for the functions that generate the statistical realizations of
each cascade by Monte-Carlo simulation. The structures and functions related to our data
structures are prototyped below:

.. doxygenstruct:: cli
   :members:

.. doxygenstruct:: cri
   :members:

There are also some utility functions that are used for reading the cascade input files, and they
are prototyped below. 

.. doxygenfunction:: readCascadeDistributionFile 
   :project: nrCascadeSim

.. doxygenfunction:: interpretDbl 
   :project: nrCascadeSim

.. doxygenfunction:: interpretSn 
   :project: nrCascadeSim

.. doxygenfunction:: interpretWeisskopf 
   :project: nrCascadeSim

.. doxygenfunction:: interpretElevVector 
   :project: nrCascadeSim

.. doxygenfunction:: interpretTauVector 
   :project: nrCascadeSim

.. doxygenfunction:: vsplit 
   :project: nrCascadeSim

The functions below provide the functionality to calculate various details of the
atom/ion trajectories for the supported elements: germanium, silicon, argon, neon. At this time
there are separate functions for each of the supported elements; this is meant to be unified in
the future in order to support a wider range of elements. For now we always use
constant-acceleration `S2` stopping. `S2` refers to the parameter from the Lindhard paper [Lindhard1963]_. 

.. doxygenfunction:: Cascade 
   :project: nrCascadeSim

.. doxygenfunction:: geCascade 
   :project: nrCascadeSim

.. doxygenfunction:: geDecay 
   :project: nrCascadeSim

.. doxygenfunction:: geStop 
   :project: nrCascadeSim

.. doxygenfunction:: rgeS2 
   :project: nrCascadeSim

.. doxygenfunction:: vgeS2 
   :project: nrCascadeSim

.. doxygenfunction:: vgeS2func 
   :project: nrCascadeSim

.. doxygenfunction:: siCascade 
   :project: nrCascadeSim

.. doxygenfunction:: siDecay 
   :project: nrCascadeSim

.. doxygenfunction:: siStop 
   :project: nrCascadeSim

.. doxygenfunction:: rsiS2 
   :project: nrCascadeSim

.. doxygenfunction:: vsiS2 
   :project: nrCascadeSim

.. doxygenfunction:: vsiS2func 
   :project: nrCascadeSim

.. doxygenfunction:: arCascade 
   :project: nrCascadeSim

.. doxygenfunction:: arDecay 
   :project: nrCascadeSim

.. doxygenfunction:: arStop 
   :project: nrCascadeSim

.. doxygenfunction:: rarS2 
   :project: nrCascadeSim

.. doxygenfunction:: varS2 
   :project: nrCascadeSim

.. doxygenfunction:: varS2func 
   :project: nrCascadeSim

.. doxygenfunction:: neCascade 
   :project: nrCascadeSim

.. doxygenfunction:: neDecay 
   :project: nrCascadeSim

.. doxygenfunction:: neStop 
   :project: nrCascadeSim

.. doxygenfunction:: rneS2 
   :project: nrCascadeSim

.. doxygenfunction:: vneS2 
   :project: nrCascadeSim

.. doxygenfunction:: vneS2func 
   :project: nrCascadeSim

---------------
`edepmath.h`
---------------

In this header is contained prototypes for supporting mathematical functions. Most of the content
are functions to assist with drawing numbers from specific probability
distributions.  


.. doxygenfile:: edepmath.h
   :project: nrCascadeSim

---------------
`lindhard.h`
---------------

In this header is contained prototypes for functions to furnish simple representations of the
Lindhard ionization model [Lindhard1963]_. They generally help return the ionization yield fraction given at
a particular starting energy (in **eV**). There are also specified functions to return the
ionization for an atom slowing down from one starting energy to another (as would happen in one
step of the cascade). Again, as in `cascadeProd.h` there are separate functions for each isotope
currently and this is intended to be unified in the future.   


.. doxygenfile:: lindhard.h
   :project: nrCascadeSim

---------------
`weisskopf.h`
---------------

In this header is contained a prototype for obtaining the Weisskopf decay-time estimate
[Weisskopf1951]_ for a gamma decay of a certain energy (in **MeV**) and certain multipolarity
(like `M1`, `E1`, etc.).


.. doxygenfile:: weisskopf.h
   :project: nrCascadeSim

----------------
`isotope_info.h`
----------------

In this header is contained prototypes for getting various isotope information. In the future this
should be replaced with a more robust API to a database to get all of this information. For now,
the information needed is hard-coded in the library. 

.. doxygenfile:: isotope_info.h
   :project: nrCascadeSim


-------------------
`rootUtil.h`
-------------------

In this header is contained prototypes for interfacing with the `ROOT` [ROOT1997]_ system. This is
only for the  writing of the output file.  

.. doxygenfile:: rootUtil.h
   :project: nrCascadeSim
==================
References
==================

.. .. [ref] author, *title*, `journal <>`_, `[arXiv:xxxx] <https://arxiv.org/abs/xxxx>`_
.. https://journals.aps.org/pr/abstract/10.1103/PhysRev.83.1073 Weisskopf estimates
.. .. [Weisskopf1951] V.F. Weisskopf, *Radiative transition probabilities in nuclei*, `Phys.Rev. 83 (1951) 1073 <https://journals.aps.org/pr/abstract/10.1103/PhysRev.83.1073>`_, `[arXiv:0203002] <https://arxiv.org/abs/0203002>`_.
.. [Lindhard1963] J. Lindhard, V. Nielsen, M. Scharff, and P.V. Thomsen, *Integral equations governing radiation effects (notes on atomic collisions III)*, `Kgl. Danske Videnskab., Selskab. Mat. Fys. Medd. 10 (1963) <https://www.osti.gov/biblio/4701226-integral-equations-governing-radiation-effects-notes-atomic-collisions-iii>`_
.. [ROOT1997] Rene Brun and Fons Rademakers , *ROOT  - An object oriented data analysis framework*, `Nucl. Inst. & Meth. in Phys. Res. A 389 (1997); Release v6.22/00 <https://doi.org/10.5281/zenodo.3895852>`_
.. [Weisskopf1951] V.F. Weisskopf, *Radiative transition probabilities in nuclei*, `Phys.Rev. 83 (1951) 1073 <https://journals.aps.org/pr/abstract/10.1103/PhysRev.83.1073>`_
=================
Contact & Support
=================

The authors of *nrCascadeSim* are `A.N. Villano <https://github.com/villaa>`_, `K. Harris
<https://github.com/gerudo7>`_, and S. Brown.

For questions, support, bug reports, or other suggestions, please open an issue
on `github <https://github.com/villano-lab/nrCascadeSim/issues>`_.
=================================
5. Documentation and Sphinx 
=================================

This documentation is hosted via `readthedocs <https://nrcascadesim.readthedocs.io/en/latest/>`_
but can also be built using `sphinx
<https://www.sphinx-doc.org/en/master/tutorial/getting-started.html>`_ in a stand-alone way. 

To build the documentation do the following starting at the root directory of the `repository
<https://github.com/villano-lab/nrCascadeSim>`_. 

.. code-block:: bash

   pip install sphinx
   sphinx-build --version
   sphinx-build -b html docs/source/ docs/build/html

This will build html documentation in `docs/build/html/index.html` which you can view with your
browser. 

After making the documentation there will be a `Makefile` inside the `docs/` directory that can be
used for further building. You can use the following commands, then, to re-build the `html`
documentation, or build the `epub` documentation, or `latex` and `pdf` documentation.

.. code-block:: bash

   cd docs/
   make html
   make epub
   make latexpdf

.. **protoSENSEI@surface**
.. ^^^^^^^^^^^^^^^^^^^^^^^
.. 
.. * **SENSEI: First Direct-Detection Constraints on sub-GeV Dark Matter from a Surface Run**
..   
..   SENSEI Collaboration (Michael Crisler et al.)
.. 
..   .. image:: https://img.shields.io/badge/Phys.Rev.Lett.-121(2018)no.6-255773.svg
..       :target: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.061803
..       :alt: Phys.Rev.Lett. 121 (2018) no.6, 061803
..   .. image:: https://img.shields.io/badge/arXiv-1804.00088-B31B1B.svg
..       :target: https://arxiv.org/abs/1804.00088
..       :alt: [arXiv:1804.00088]


