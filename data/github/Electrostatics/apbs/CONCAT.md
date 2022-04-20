# Contributing to APBS

We are excited to have you contribute to APBS!
If you have general questions about the code, please see the [documentation](https://apbs.readthedocs.io/en/latest/) and/or start a thread in the [discussion forum](https://github.com/Electrostatics/apbs/discussions). 
Once you have familiarity with the code, you're ready to start contributing:

## Find or create a problem

Find a problem/feature that needs to be resolved in the code [issues](https://github.com/Electrostatics/apbs/issues).
If the problem you're trying to solve doesn't exist, create a Github issue to resolve some or all of the problem or to add features.
For example, you could
  * Fix a bug
  * Add a new capability with 3 parts - create issue for each part it possible
  * In general, make problem/feature small enough to finish it within a week or two

## Create a branch

Create a git branch using the convention of `github_user`/`issue_#` where `github_user` is your user name and `issue_#` is the issue number from step 2 above.

## Create and pass tests

Create a test that replicates the problem/feature and fails and show how your fix results in a working test.
Work on the code until your test, as well as all previous tests, pass.

## Do not submit messy code

Run your code through formatting (e.g., [psf/black](https://github.com/psf/black)) and linting (e.g., [pylint](https://www.pylint.org/) and [flake8](https://flake8.pycqa.org/en/latest/)) tools appropriate to the language for your code.

## Submit a pull request

1. Commit your changes to Git and push your branch to the Github repo.
2. Create a [pull request](https://github.com/Electrostatics/apbs/compare?expand=1) and add reviewers (e.g., `intendo` or `sobolevnrm`) to the request.
[Reference the issue you are trying to fix](https://docs.github.com/en/github/managing-your-work-on-github/linking-a-pull-request-to-an-issue) in your pull request.
3. If the pull request passes [Github Actions](https://github.com/features/actions) and peer review, then the branch will be merged with `main`, your branch will be deleted, and the development team will be very grateful for your contribution!

Thank you for considering to contribute to our code!
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

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
reported by contacting the project team at nathanandrewbaker@gmail.com. All
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

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
APBS - Adaptive Poisson-Boltzmann Solver
========================================

[![Documentation Status](https://readthedocs.org/projects/apbs/badge/?version=latest)](https://apbs.readthedocs.io/en/latest/?badge=latest)
[![Github Action Build Status](https://github.com/Electrostatics/apbs/workflows/Build/badge.svg)](https://github.com/Electrostatics/apbs/actions)

This repository contains the APBS software.

For more information about APBS, please see

* Home page:  http://www.poissonboltzmann.org/
* Documentation: http://apbs.readthedocs.io
APBS 3.2 CHANGELOG
==================

These are notes for APBS version 3.2
------------------------------------

* Binary releases may be found on [GitHub](https://github.com/Electrostatics/apbs/releases) and on [SourceForge](http://sourceforge.net/projects/apbs/files/apbs).

### New Features

* Poisson-Boltzmann Analytical Method (PYGBE, see [ADD_AUTHORS](ADD_LINK)) integrated with APBS. PYGBE can be installed from the "slic" branch of the [PYGE Repository](https://github.com/matiasmartineza/pygbe.git).

* Poisson-Boltzmann Analytical Method (PBAM, see [Lotan & Head-Gordon](http://pubs.acs.org/doi/full/10.1021/ct050263p)) and Semi-Analytical Method (PBSAM, see [Yap & Head-Gordon](http://pubs.acs.org/doi/abs/10.1021/ct100145f)) integrated with APBS. PBSAM is currently only available in the Linux and OS X distributions.
    - Examples are located with the APBS examples in the pbam/ and pbsam/ directories.
    - More information and documentation may be found in the [PBAM](http://www.poissonboltzmann.org/external_contributions/extern-pbam/) and [PBSAM](http://www.poissonboltzmann.org/external_contributions/extern-pbsam/) sections of the APBS website.
* Tree-Code Accelerated Boundary Integral Poisson-Boltzmann Method (TABI-PB) integrated with APBS.(See [Geng & Krasny](http://www.sciencedirect.com/science/article/pii/S0021999113002404))
    - Examples are located with the APBS examples in the bem/, bem-pKa/, and bem-binding-energies/ folders
    - Included NanoShaper alternative to MSMS.
    - More information and documentation may be found in the [Contributions](http://www.poissonboltzmann.org/external_contributions/extern-tabi/) section of the APBS website
* Added binary DX format support to the appropriate APBS tools.
* Test suite amended and expanded.
* Removed hard-coded limitation to number of grid points used to determine surface accessibility.
* Moved [PDB2PQR](https://github.com/Electrostatics/pdb2pqr) components to it's own repository which can be installed separately.

### Known Bugs / Limitations

* FETK not building in windows due to C standard restrictions in the Microsoft compiler implementation.

### Minor Updates

* PB(S)AM now requires the key work 'pos' for the term argument.
* PB(S)AM 'surf' keyword has been replaced with the 'usemesh' keyword.
* PB(S)AM 'salt' keyword has been replaced with the 'ion' keyword.
* PB(S)AM dynamics parameters are no longer accepted in the ELEC section.
* PB(S)AM now has only one type of ELEC method: pb(s)am_auto.
* PB(S)AM 'gridpts' keyword has been replaced with 'dime' keyword.
* PB(S)AM 'dx' and '3dmap' keywords are deprecated to use the 'write' one instead.
* BEM mesh keyword now requires method names instead of just integer values.
* GEOFLOW ELEC type has been change from 'geoflow-auto' to 'geoflow'.
* Fixed miscellaneous Windows build issues.
* Update the build configurations for the Pythons libraries.
* Removed unused files that no longer worked (e.g. Makefile in example directories)
* Updated Sphinx configuration to build API documentation

### Notes

* The following are included in APBS as Git submodules:
- Geometric Flow ([link](https://github.com/Electrostatics/geoflow_c/tree/e8ce510a670e0b7f3501e72be6141fc20328f947))
- FETk ([link](https://github.com/Electrostatics/FETK/tree/0c6fdeabe8929acea7481cb1480b5706b343b7e0))
- PBAM/PBSAM ([link](https://github.com/davas301/pb_solvers/tree/4805cbec02b30e9bae927f03ac2fecd3217c4dad))
- TABI-PB ([link](https://github.com/lwwilson1/TABIPB/tree/941eff91acd4153a06764e34d29b633c6e3b980f))
APBS Release Procedure
-----------------------
 1. Change Version Number
	 - [ ] Edit [VERSION]([https://github.com/Electrostatics/apbs/blob/main/apbs/VERSION)
		Increment the value after the comment block, which is in the form
	     M_m_u
	     Where:
		 - M is the Major version - increment if there are breaking changes or dropping support for previous features
		 - m is the Minor version - increment for new features added
		 - u is the Micro version - increment for adding small changes like new tests or fixing small bugs

 2. Update the Releases document
	 - [ ] Edit [docs/releases.rst]([https://github.com/Electrostatics/apbs/blob/main/apbs/docs/releases.rst)
	   - Document major/minor changes for this release
   
 3. Update License info
	   - [ ] Update license dates and information in source files
	   - In apbs/src edit all .c source files and all .h header files, update dates
  
 4. Create a Pull Request (PR)
     - [ ] Create a new [Pull Request](https://github.com/Electrostatics/apbs/pulls)
		 - Base branch should be `release`
		 - Source branch should be `main`
		 - Briefly describe the changes included

 5. Check tests
     - Go to the [Actions](https://github.com/Electrostatics/apbs/actions) tab in GitHub
     - Tests are performed for three target platforms:
       - Ubuntu
	   - MacOSX
	   - Windows
     - [ ] Ensure that the builds and associated tests were successful
	 - [ ] Ensure that the use tests were successful
	 - [ ] Ensure that the build artifacts were uploaded to the action

 6. Merge the PR
	 - [ ] Ensure that the [Release](https://github.com/Electrostatics/apbs/releases) is correctly created
	 - [ ] Ensure that the builds and associated tests were successful
	 - [ ] Ensure that the use tests were successful
	 - [ ] Ensure that the build artifacts were uploaded to the Release

 7. Update http://www.poissonboltzmann.org/apbs/release-history with new release information.

 8. Pat yourself on the back for a job well done!
# APBS Tools and Utilities

Please see [tools and utilities](https://apbs.readthedocs.io/en/latest/using/tools.html)
External Packages
===================

If you want to update the submodule to the 
latest commit available, that isn't done with

`git submodule update`

To update the commit, go to the submodue directory:
`externals/submodule/` and do the following:

1. git checkout master
2. git pull

This will update it within the submodule. Then you can

`cd ../../`

And add changes to the APBS repo.
APBS validation and test cases
==============================

Please see [tests cases](https://apbs.readthedocs.io/en/latest/using/tests.html)# APBS version validation and test cases

## APBS examples and test cases

This directory serves as the root directory for the APBS test suite.
In each directory you will find example input files to use with APBS and a README file displaying the results for different versions of APBS.

Executing <code>make test</code> in each directory will run the examples for that directory and log the results to <code>TESTRESULTS.log</code>.
Executing <code>make test</code> from the root examples directory will run all the tests listed below.
Tests will either pass, pass with rounding error (within 10<sup>-9</sup> of the expected result), or fail outright.

| Example | README file | Source | Description | 
| ---- | ---- | ---- | ---- | 
| Actin dimer (actin-dimer) | [actin-dimer/README.md](actin-dimer/README.md) | Dave Sept | Calculate binding energies for actin dimers. This is an example of a large biomolecule binding energy calculation that often requires parallel focusing. | 
| Alkane nonpolar solvation energies (alkanes) | [alkanes/README.md](alkanes/README.md) | Nathan Baker, Jason Wagoner | Calculate nonpolar solvation energies for various alkanes.  Taken from Wagoner JA, Baker NA. Assessing implicit models for nonpolar mean solvation forces: the importance of dispersion and volume terms. [Proc Natl Acad Sci USA, 103, 8331-8336, 2006.](http://dx.doi.org/10.1073/pnas.0600118103) |
| Born ion (born) | [born/README.md](born/README.md) | Nathan Baker | Calculate solvation energies for ions of various sizes and compare to the analytical results. |
| FKBP (FKBP) | [FKBP/README.md](FKBP/README.md) | Jung-Hsin Lin | Binding of various (small) ligands to FKBP.  Analogous to HCA binding case (except it works). |
| HCA ligand binding (hca-bind) | [hca-bind/README.md](hca-bind/README.md) | UHBD | Calculate the binding of a small molecule (acetazolamide) to a medium-sized protein (human carbonic anhydrase). |
| Acetic acid ionization (ionize) | [ionize/README.md](ionize/README.md) | UHBD | Calculate electrostatic contributions to the ionization energy of acetic acid. | 
| Ion-ion PMF (ion-pmf) | [ion-pmf/README.md](ion-pmf/README.md) | Nathan Baker | Calculate solvation energies and solvation force components for ion pairs. |
| Ion-protein interaction energies (ion-protein) | [ion-protein/README.md](ion-protein/README.md) | Dave Sept | Calculate polar energy of placing an ion near a macromolecule. |
| PKA-balanol binding (pka-lig) | [pka-lig/README.md](pka-lig/README.md) | Chung Wong | Calculate binding energies of a ligand to protein kinase A. |
| PKA-balanol binding (pka-lig)/UHDB | [pka-lig/UHDB/readme.md](pka-lig/UHDB/readme.md)| Chun Wong | Shows the calculations done using molecular surface based and van der Waals dielectric definitions. |
| Coulomb's law (point-pmf) | [point-pmf/README.md](point-pmf/README.md) | Nathan Baker | See how well we do reproducing Coulomb's law. |
| Methanol solvation (solv) | [solv/README.md](solv/README.md) | UHBD | Calculate the solvation energies of methanol and methoxide. | 
| Protein-RNA interactions (protein-rna) | [protein-rna/README.md](protein-rna/README.md) | David Draper | Calculate the salt dependence of protein interactions with box B RNA hairpin. |
| Geometric flow solvation model | [geoflow/README.md](geoflow/README.md) | Elizabeth Jurrus | Calculate the dielectric interface profile across the solute-solvent boundary in a thermodynamically sef-consistent fashion. |
| Ion Binding to DNA Duplexes Using SMPBE| [smpbe/readme.md](smpbe/reamdme.md) | | Calculate PBE taking into account the finite ion size. |
| Focusing Membrane Boundary Condition| [membrane/readme.md](membrane/readme.md) | | Solve the PBE with a single atom using focusing membrane boundary conditions. |
| NMR Structure of the RNA binding Domain | [bem/readme.md](bem/readme.md) | | Calculate the solvation complex using the boundary element method as is implemented in APBS. |
| Born Ion | [opal/README.md](opal/README.md) | Nathan Baker | The Born ion is a canonical electrostatic's test case for which there is an analytical solution. This example examines the solvation free energy. |
| Miscellaneous | [misc/README.md](misc/README.md) | | A collection of pqr files of molecules that have interesting potentials. |README for Born APBS examples
=============================

This is the canonical electrostatics test case: Born ion. A non-polarizable ion with a single embedded point charge; has an analytical solution for the potential. We examine the solvation free energy as a function of ionic radius.

Please see apbs.in for details on the particular solvation energy calculations. Analytical results are given in pmf.dat.

This example was contributed by Nathan Baker.

Input File|Description|APBS Version|Results (kJ/mol)|Analytical (kJ/mol)
---|---|---|---|---
[apbs-mol-auto.in](apbs-mol-auto.in)|Sequential, 3 A sphere, 3-level focusing to 0.188 A, srfm mol|**1.5**|**-229.7740**|-230.62
|||1.4.2|-229.774
|||1.4.1|-229.7736|
|||1.4|-229.7736<sup>[3](#3)</sup>
|||1.3|-229.7735
|||1.2.1|-229.7735
|||1.2|-229.7735
|||1.1.0|-229.7735
|||1.0.0|-229.7735
|||0.5.1|-229.7735
|||0.5.0|-229.7735
|||0.4.0|-229.7735<sup>[1](#1)</sup>
|||0.3.2|-229.7248
|||0.3.1|-229.7248
|||0.3.0|-229.7248
|||0.2.6|-229.7248
|||0.2.5|-229.7248
|||0.2.4|-227.1859
|||0.2.3|-227.1589
|||0.2.2|-227.186
|||0.2.1|-227.19
|||0.2.0|-227.19
|||0.1.8|-227.19
[apbs-smol-auto.in](apbs-smol-auto.in)|Sequential, 3 A sphere, 3-level focusing to 0.188 A, srfm smol|**1.5**|**-229.012**|-230.62
|||1.4.2|-229.0124
|||1.4.1|-229.0124
|||1.4|-229.0124
|||1.3|-229.0124
|||1.2.1|-229.0124
|||1.2|-229.0124<sup>[2](#2)</sup>
|||1.1.0|-229.0123
|||1.0.0|-229.0123
|||0.5.1|-229.0123
|||0.5.0|-229.0123
|||0.4.0|-229.0123
[apbs-mol-parallel.in](apbs-mol-parallel.in)|Parallel with 4 processors, 3 A sphere, focusing to 0.103 A, srfm mol|**1.5**|**-230.492**|-230.62
|||1.4.2|-230.492
|||1.4.1|-230.4918<sup>[4](#4)</sup>|
|||1.4|-230.4919<sup>[3](#3)</sup>
|||1.3|-230.4918
|||1.2.1|-230.4918
|||1.2|-230.4918<sup>[2](#2)</sup>
|||1.1.0|-230.4916
|||1.0.0|-230.4916
|||0.5.1|-230.4916
|||0.5.0|-230.4916
|||0.4.0|-230.4916
|||0.2.1|-230.77
[apbs-smol-parallel.in](apbs-smol-parallel.in)|Parallel with 4 processors, 3 A sphere, focusing to 0.103 A, srfm smol|**1.5**|**-229.387**|-230.62
|||1.4.2|-229.387
|||1.4.1|-229.3871
|||1.4|-229.3871
|||1.3|-229.3871
|||1.2.1|-229.3871
|||1.2|-229.3871<sup>[2](#2)</sup>
|||1.1.0|-229.3872
|||1.0.0|-229.3872
|||0.5.1|-229.3872
|||0.5.0|-229.3872
|||0.4.0|-229.3872<sup>[1](#1)</sup>
|||0.3.2|-226.3529
|||0.3.1|-226.3529
|||0.3.0|-229.5849
|||0.2.6|-229.5849
|||0.2.5|-229.5849
|||0.2.4|-226.2276
|||0.2.3|-226.2276
|||0.2.2|-226.2276
|||0.2.0|-226.228
|||0.1.8|-226.23
[apbs-smol-parallel.in](apbs-mol-fem.in)|Finite Element Method, 3 A sphere, 3-level focusing to 0.188 A, srfm mol|**1.4.1-binary**|**-231.9550**|-230.62
[apbs-smol-parallel.in](apbs-smol-fem.in)|Finite Element Method, 3 A sphere, 3-level focusing to 0.188 A, srfm smol|**1.4.1-binary**|**-230.9760**|-230.62

<a name=1></a><sup>1</sup> The discrepancy in values between versions 0.4.0 and 0.3.2 is most likely due to three factors:

-   A bug fix in Vacc_molAcc which removed spurious regions of high internal dielectric values
-   A switch in the algorithm used to compute the dielectric smoothing for srfm smol.
-   The addition of the Vacc sphere density (sdens keyword) as a variable and a change in the default sdens value from 3.0 to 10.0

<a name=2></a><sup>2</sup> APBS 1.2 has switched the multigrid smoothing algorithm from standard Gauss-Seidel to Gauss-Seidel red/black in order to facilitate parallelization. This switch has caused small differences in individual calculation energies which, when combined to the final answer, create larger errors (up to 0.04%). These errors can be reduced by resetting the APBS error tolerance to 1e-9 or smaller values. For a more detailed explanation, please see the APBS FAQ, [here](http://www.poissonboltzmann.org/docs/apbs-faq/#sources error calculation).

<a name=3></a><sup>3</sup> The discrepancy in values between versions 1.3 and 1.4 is most likely due to the following factor(s):

-   Translation of contrib/pmgZ library from FORTRAN to C
-   Differences in numerical implementations between FORTRAN and C compilers result in small round-off discrepencies
-   Small margins due to these round-off discrepencies acumulate in the computations

<a name=4></a><sup>4</sup> The discrepancy in the result between versions 1.4 and 1.4.1-binary is most likely due to a reporting error.

Please see the ChangeLog or the [APBS website](http://www.poissonboltzmann.org/) for more information.


README for TABI-PB (Boundary Element Method) Binding Energy Examples
====================================================================

The example input files included in this folder uses a boundary element approach called
TABI-PB to solve the PBE. BEMs have the characteristic that only the boundary of the
domain has to be discretized. This is particularly useful for problems in which the data
of interest is at the boundary of the solution.

This directory contains three example .in files:
        1. 1d30.in
        2. 1d30_monomer1.in
        3. 1d30_monomer2.in

These files provide an example to demonstrate the calculation of binding energy on 1d30.
More details are available on the APBS website contributions section.

Input File| APBS Version| Result (kJ/mol)| Expected (kJ/mol)
---|---|---|---
[1d30.in](1d30.in)| **3.0**| **-22113.098**| **-22113.098**
[1d30_monomer1.in](1d30_monomer1.in)| **3.0**| **-26225.275**| **-26225.275**
[1d30_monomer2.in](1d30_monomer2.in)| **3.0**| **-779.948**| **-779.948**
The example input files in this directory calculate nonpolar solvation energies for alkanes based on the protocol described in Wagoner JA, Baker NA. Assessing implicit models for nonpolar mean solvation forces: the importance of dispersion and volume terms. Proc Natl Acad Sci USA, 103, 8331-8336, 2006. [(http://dx.doi.org/10.1073/pnas.0600118103)](http://dx.doi.org/10.1073/pnas.0600118103).

This example was contributed by Nathan Baker and Jason Wagoner.


APBS Version|Alkane|SASA (Å<sup>2</sup>)|SASA energy (kJ/mol)|SAV (Å<sup>3</sup>)|SAV energy (kJ/mol)|WCA energy (kJ/mol)|Total nonpolar solvation energy (kJ/mol)
---|---|---|---|---|---|---|---
**3.0**|2-methylbutane|214.202|1.82072|253.665|60.7274|-48.1507|1.439739455792E+01
||butane|193.855|1.64777|217.863|52.1564|-41.7207|1.208346456826E+01
||cyclohexane|221.799|1.88529|267.435|64.0239|-52.3691|1.354016672221E+01
||cyclopentane|193.638|1.64593|217.998|52.1887|-44.471|9.363673200142E+00
||ethane|139.427|1.18513|140.346|33.5988|-25.3612|9.422717598546E+00
||hexane|250.291|2.12748|298.053|71.3539|-57.0807|1.640068943201E+01
||isobutane|192.744|1.63832|218.943|52.415|-40.8218|1.323144287435E+01
||methane|105.42|0.896066|95.985|22.9788|-15.9805|7.894367190329E+00
||neopentane|210.755|1.79141|251.127|60.1198|-47.4149|1.449633815052E+01
||pentane|222.524|1.89145|258.93|61.9878|-49.4003|1.447900211546E+01
||propane|170.391|1.44832|183.573|43.9474|-33.4721|1.192358496286E+01


APBS Version|Alkane|SASA (Å<sup>2</sup>)|SASA energy (kJ/mol)|SAV (Å<sup>3</sup>)|SAV energy (kJ/mol)|WCA energy (kJ/mol)|Total nonpolar solvation energy (kJ/mol)
---|---|---|---|---|---|---|---
**1.5**|2-methylbutane|214.202|1.82072|253.665|60.7274|-48.1507|1.439739455792E+01
||butane|193.855|1.64777|217.863|52.1564|-41.7207|1.208346456826E+01
||cyclohexane|221.799|1.88529|267.435|64.0239|-52.3691|1.354016672221E+01
||cyclopentane|193.638|1.64593|217.998|52.1887|-44.471|9.363673200142E+00
||ethane|139.427|1.18513|140.346|33.5988|-25.3612|9.422717598546E+00
||hexane|250.291|2.12748|298.053|71.3539|-57.0807|1.640068943201E+01
||isobutane|192.744|1.63832|218.943|52.415|-40.8218|1.323144287435E+01
||methane|105.42|0.896066|95.985|22.9788|-15.9805|7.894367190329E+00
||neopentane|210.755|1.79141|251.127|60.1198|-47.4149|1.449633815052E+01
||pentane|222.524|1.89145|258.93|61.9878|-49.4003|1.447900211546E+01
||propane|170.391|1.44832|183.573|43.9474|-33.4721|1.192358496286E+01


APBS Version|Alkane|SASA (Å<sup>2</sup>)|SASA energy (kJ/mol)|SAV (Å<sup>3</sup>)|SAV energy (kJ/mol)|WCA energy (kJ/mol)|Total nonpolar solvation energy (kJ/mol)
---|---|---|---|---|---|---|---
**1.4.2**|2-methylbutane|214.202|1.82072|253.665|60.7274|-48.1507|1.439740000000e+01
||butane|193.855|1.64777|217.863|52.1564|-41.7207|1.208350000000e+01
||cyclohexane|221.799|1.88529|267.435|64.0239|-52.3691|1.354020000000e+01
||cyclopentane|193.638|1.64593|217.998|52.1887|-44.471|9.363670000000e+00
||ethane|139.427|1.18513|140.346|33.5988|-25.3612|9.422720000000e+00
||hexane|250.291|2.12748|298.053|71.3539|-57.0807|1.640070000000e+01
||isobutane|192.744|1.63832|218.943|52.415|-40.8218|1.323140000000e+01
||methane|105.42|0.896066|95.985|22.9788|-15.9805|7.894370000000e+00
||neopentane|210.755|1.79141|251.127|60.1198|-47.4149|1.449630000000e+01
||pentane|222.524|1.89145|258.93|61.9878|-49.4003|1.447900000000e+01
||propane|170.391|1.44832|183.573|43.9474|-33.4721|1.192360000000e+01


APBS Version|Alkane|SASA (Å<sup>2</sup>)|SASA energy (kJ/mol)|SAV (Å<sup>3</sup>)|SAV energy (kJ/mol)|WCA energy (kJ/mol)|Total nonpolar solvation energy (kJ/mol)
---|---|---|---|---|---|---|---
**1.4.1-binary**|2-methylbutane|214.202|1.82072|253.665|60.7274|-48.1507|1.439739455792E+01
||butane|193.855|1.64777|217.863|52.1564|-41.7207|1.208346456826E+01
||cyclohexane|221.799|1.88529|267.435|64.0239|-52.3691|1.354016672221E+01
||cyclopentane|193.638|1.64593|217.998|52.1887|-44.471|9.363673200142E+00
||ethane|139.427|1.18513|140.346|33.5988|-25.3612|9.422717598546E+00
||hexane|250.291|2.12748|298.053|71.3539|-57.0807|1.640068943201E+01
||isobutane|192.744|1.63832|218.943|52.415|-40.8218|1.323144287435E+01
||methane|105.42|0.896066|95.985|22.9788|-15.9805|7.894367190329E+00
||neopentane|210.755|1.79141|251.127|60.1198|-47.4149|1.449633815052E+01
||pentane|222.524|1.89145|258.93|61.9878|-49.4003|1.447900211546E+01
||propane|170.391|1.44832|183.573|43.9474|-33.4721|1.192358496286E+01


APBS Version|Alkane|SASA (Å<sup>2</sup>)|SASA energy (kJ/mol)|SAV (Å<sup>3</sup>)|SAV energy (kJ/mol)|WCA energy (kJ/mol)|Total nonpolar solvation energy (kJ/mol)
---|---|---|---|---|---|---|---
**1.4**|2-methylbutane|214.202|1.82072|253.665|60.7274|-48.1507|1.439739455792E+01
||butane|193.855|1.64777|217.863|52.1564|-41.7207|1.208346456826E+01
||cyclohexane|221.799|1.88529|267.435|64.0239|-52.3691|1.354016672221E+01
||cyclopentane|193.638|1.64593|217.998|52.1887|-44.471|9.363673200142E+00
||ethane|139.427|1.18513|140.346|33.5988|-25.3612|9.422717598546E+00
||hexane|250.291|2.12748|298.053|71.3539|-57.0807|1.640068943201E+01
||isobutane|192.744|1.63832|218.943|52.415|-40.8218|1.323144287435E+01
||methane|105.42|0.896066|95.985|22.9788|-15.9805|7.894367190329E+00
||neopentane|210.755|1.79141|251.127|60.1198|-47.4149|1.449633815052E+01
||pentane|222.524|1.89145|258.93|61.9878|-49.4003|1.447900211546E+01
||propane|170.391|1.44832|183.573|43.9474|-33.4721|1.192358496286E+01

APBS Version|Alkane|SASA (Å<sup>2</sup>)|SASA energy (kJ/mol)|SAV (Å<sup>3</sup>)|SAV energy (kJ/mol)|WCA energy (kJ/mol)|Total nonpolar solvation energy (kJ/mol)
---|---|---|---|---|---|---|---
**1.3**|2-methylbutane|214.202|1.82072|253.906|60.7852|-48.3035|1.439739455792E+01
||butane|193.855|1.64777|218.119|52.2176|-41.6443|1.208346456826E+01
||cyclohexane|221.799|1.88529|267.165|63.9593|-52.4787|1.354016672221E+01
||cyclopentane|193.638|1.64593|218.412|52.2877|-44.3607|9.363673200142E+00
||ethane|139.427|1.18513|140.692|33.6818|-25.4382|9.422717598546E+00
||hexane|250.291|2.12748|297.703|71.2701|-57.1544|1.640068943201E+01
||isobutane|192.744|1.63832|219.153|52.4652|-40.8617|1.323144287435E+01
||methane|105.42|0.896066|95.4242|22.8446|-15.9414|7.894367190329E+00
||neopentane|210.755|1.79141|251.314|60.1647|-47.4807|1.449633815052E+01
||pentane|222.524|259.149|1.89145|62.0403|-49.456|1.447900211546E+01
||propane|170.391|182.703|1.44832|43.739|-33.4629|1.192358496286E+01

APBS Version|Alkane|SASA (Å<sup>2</sup>)|SASA energy (kJ/mol)|SAV (Å<sup>3</sup>)|SAV energy (kJ/mol)|WCA energy (kJ/mol)|Total nonpolar solvation energy (kJ/mol)
---|---|---|---|---|---|---|---
1.2.1|2-methylbutane|214.202|1.82072|253.906|60.7852|-48.3035|1.439739455792E+01
||butane|193.855|1.64777|218.119|52.2176|-41.6443|1.208346456826E+01
||cyclohexane|221.799|1.88529|267.165|63.9593|-52.4787|1.354016672221E+01
||cyclopentane|193.638|1.64593|218.412|52.2877|-44.3607|9.363673200142E+00
||ethane|139.427|1.18513|140.692|33.6818|-25.4382|9.422717598546E+00
||hexane|250.291|2.12748|297.703|71.2701|-57.1544|1.640068943201E+01
||isobutane|192.744|1.63832|219.153|52.4652|-40.8617|1.323144287435E+01
||methane|105.42|0.896066|95.4242|22.8446|-15.9414|7.894367190329E+00
||neopentane|210.755|1.79141|251.314|60.1647|-47.4807|1.449633815052E+01
||pentane|222.524|259.149|1.89145|62.0403|-49.456|1.447900211546E+01
||propane|170.391|182.703|1.44832|43.739|-33.4629|1.192358496286E+01

APBS Version|Alkane|SASA (Å<sup>2</sup>)|SASA energy (kJ/mol)|SAV (Å<sup>3</sup>)|SAV energy (kJ/mol)|WCA energy (kJ/mol)|Total nonpolar solvation energy (kJ/mol)
---|---|---|---|---|---|---|---
1.2|2-methylbutane|214.202|1.82072|253.906|60.7852|-48.3035|1.439739455792E+01
||butane|193.855|1.64777|218.119|52.2176|-41.6443|1.208346456826E+01
||cyclohexane|221.799|1.88529|267.165|63.9593|-52.4787|1.354016672221E+01
||cyclopentane|193.638|1.64593|218.412|52.2877|-44.3607|9.363673200142E+00
||ethane|139.427|1.18513|140.692|33.6818|-25.4382|9.422717598546E+00
||hexane|250.291|2.12748|297.703|71.2701|-57.1544|1.640068943201E+01
||isobutane|192.744|1.63832|219.153|52.4652|-40.8617|1.323144287435E+01
||methane|105.42|0.896066|95.4242|22.8446|-15.9414|7.894367190329E+00
||neopentane|210.755|1.79141|251.314|60.1647|-47.4807|1.449633815052E+01
||pentane|222.524|259.149|1.89145|62.0403|-49.456|1.447900211546E+01
||propane|170.391|182.703|1.44832|43.739|-33.4629|1.192358496286E+01

APBS Version|Alkane|SASA (Å<sup>2</sup>)|SASA energy (kJ/mol)|SAV (Å<sup>3</sup>)|SAV energy (kJ/mol)|WCA energy (kJ/mol)|Total nonpolar solvation energy (kJ/mol)
---|---|---|---|---|---|---|---
1.1.0|2-methylbutane|214.202|1.82072|253.906|60.7852|-48.3035|1.439739455792E+01
||butane|193.855|1.64777|218.119|52.2176|-41.6443|1.208346456826E+01
||cyclohexane|221.799|1.88529|267.165|63.9593|-52.4787|1.354016672221E+01
||cyclopentane|193.638|1.64593|218.412|52.2877|-44.3607|9.363673200142E+00
||ethane|139.427|1.18513|140.692|33.6818|-25.4382|9.422717598546E+00
||hexane|250.291|2.12748|297.703|71.2701|-57.1544|1.640068943201E+01
||isobutane|192.744|1.63832|219.153|52.4652|-40.8617|1.323144287435E+01
||methane|105.42|0.896066|95.4242|22.8446|-15.9414|7.894367190329E+00
||neopentane|210.755|1.79141|251.314|60.1647|-47.4807|1.449633815052E+01
||pentane|222.524|259.149|1.89145|62.0403|-49.456|1.447900211546E+01
||propane|170.391|182.703|1.44832|43.739|-33.4629|1.192358496286E+01

APBS Version|Alkane|SASA (Å<sup>2</sup>)|SASA energy (kJ/mol)|SAV (Å<sup>3</sup>)|SAV energy (kJ/mol)|WCA energy (kJ/mol)|Total nonpolar solvation energy (kJ/mol)
---|---|---|---|---|---|---|---
1.0.0|2-methylbutane|214.202|1.82072|253.906|60.7852|-48.3035|1.439739455792E+01
||butane|193.855|1.64777|218.119|52.2176|-41.6443|1.208346456826E+01
||cyclohexane|221.799|1.88529|267.165|63.9593|-52.4787|1.354016672221E+01
||cyclopentane|193.638|1.64593|218.412|52.2877|-44.3607|9.363673200142E+00
||ethane|139.427|1.18513|140.692|33.6818|-25.4382|9.422717598546E+00
||hexane|250.291|2.12748|297.703|71.2701|-57.1544|1.640068943201E+01
||isobutane|192.744|1.63832|219.153|52.4652|-40.8617|1.323144287435E+01
||methane|105.42|0.896066|95.4242|22.8446|-15.9414|7.894367190329E+00
||neopentane|210.755|1.79141|251.314|60.1647|-47.4807|1.449633815052E+01
||pentane|222.524|259.149|1.89145|62.0403|-49.456|1.447900211546E+01
||propane|170.391|182.703|1.44832|43.739|-33.4629|1.192358496286E+01

APBS Version|Alkane|SASA (Å<sup>2</sup>)|SASA energy (kJ/mol)|SAV (Å<sup>3</sup>)|SAV energy (kJ/mol)|WCA energy (kJ/mol)|Total nonpolar solvation energy (kJ/mol)
---|---|---|---|---|---|---|---
0.5.1<sup>[1](#1)</sup>|2-methylbutane|214.202|1.82072|253.906|60.7852|-48.3035|1.439739455792E+01
||butane|193.855|1.64777|218.119|52.2176|-41.6443|1.208346456826E+01
||cyclohexane|221.799|1.88529|267.165|63.9593|-52.4787|1.354016672221E+01
||cyclopentane|193.638|1.64593|218.412|52.2877|-44.3607|9.363673200142E+00
||ethane|139.427|1.18513|140.692|33.6818|-25.4382|9.422717598546E+00
||hexane|250.291|2.12748|297.703|71.2701|-57.1544|1.640068943201E+01
||isobutane|192.744|1.63832|219.153|52.4652|-40.8617|1.323144287435E+01
||methane|105.42|0.896066|95.4242|22.8446|-15.9414|7.894367190329E+00
||neopentane|210.755|1.79141|251.314|60.1647|-47.4807|1.449633815052E+01
||pentane|222.524|259.149|1.89145|62.0403|-49.456|1.447900211546E+01
||propane|170.391|182.703|1.44832|43.739|-33.4629|1.192358496286E+01

APBS Version|Alkane|SASA (Å<sup>2</sup>)|SASA energy (kJ/mol)|SAV (Å<sup>3</sup>)|SAV energy (kJ/mol)|WCA energy (kJ/mol)|Total nonpolar solvation energy (kJ/mol)
---|---|---|---|---|---|---|---
0.5.0|2-methylbutane|214.202|1.82072|253.906|60.7852|-48.3035|1.430239579640E+01
||butane|193.855|1.64777|218.119|52.2176|-41.6443|1.222110127537E+01
||cyclohexane|221.799|1.88529|267.165|63.9593|-52.4787|1.336586748209E+01
||cyclopentane|193.638|1.64593|218.412|52.2877|-44.3607|9.572911268235E+00
||ethane|139.427|1.18513|140.692|33.6818|-25.4382|9.428701470984E+00
||hexane|250.291|2.12748|297.703|71.2701|-57.1544|1.624316652259E+01
||isobutane|192.744|1.63832|219.153|52.4652|-40.8617|1.324178842307E+01
||methane|105.42|0.896066|95.4242|22.8446|-15.9414|7.799212389992E+00
||neopentane|210.755|1.79141|251.314|60.1647|-47.4807|1.447540747648E+01
||pentane|222.524|259.149|1.89145|62.0403|-49.456|1.447574303821E+01
||propane|170.391|182.703|1.44832|43.739|-33.4629|1.172438897305E+01

<a name=1></a><sup>1</sup> The discrepancy in values between versions 0.5.1 and 0.5.0 is most likely due to the following factor(s):

-   Removal of the dime keyword and implementation of variable grid lengths

Please see the ChangeLog or the [APBS website](http://www.poissonboltzmann.org/) for more information.
README for pka-lig APBS examples
================================

The example input files in this directory calculate the binding energes of a ligand to protein kinase A.

This example was contributed by Chung Wong.

Input File|Description|APBS Version|Results (kJ/mol)|UHBD (kJ/mol)
---|---|---|---|---
[apbs-mol-vdw.in](apbs-mol-vdw.in)|2-level focusing to 0.250 A spacing, VdW surface, srfm mol|**1.5**|**8.08352**|8.876
|||1.4.2|8.08352
|||1.4.1|8.0835
|||1.4|8.0835
|||1.3|8.0835
|||1.2.1|8.0835
|||1.2|8.0835<sup>[4](#4)</sup>
|||1.1.0|8.0858
|||1.0.0|8.0858
|||0.5.1|8.0858<sup>[3](#3)</sup>
|||0.5.0|8.0640
||||0.4.0|8.0640
[apbs-smol-vdw.in](apbs-smol-vdw.in)|2-level focusing to 0.250 A spacing, VdW surface, srfm smol|**1.5**|**20.9630**|8.876
|||1.4.2|20.9630
|||1.4.1|20.9630
|||1.4|20.9630
|||1.3|20.9630
|||1.2.1|20.9630
|||1.2|20.9630<sup>[4](#4)</sup>
|||1.1.0|20.9628
|||1.0.0|20.9628
|||0.5.1|20.9628<sup>[3](#2)</sup>
|||0.5.0|20.9542
|||0.4.0|20.9542<sup>[2](#2)</sup>
|||0.3.2|8.0640<sup>[1](#1)</sup>
|||0.3.1|6.6465
|||0.3.0|6.6465
|||0.2.6|6.6465
|||0.2.5|6.6465
|||0.2.4|6.6465
|||0.2.3|6.6465
|||0.2.2|6.6465
|||0.2.1|6.647
|||0.2.0|6.647
|||0.1.8|6.65
[apbs-mol-surf.in](apbs-mol-surf.in)|2-level focusing to 0.250 A spacing, molecular surface, srfm mol|**1.5**|**119.2610**|86.50
|||1.4.2|119.2610
|||1.4.1|119.2608
|||1.4|119.2608
|||1.3|119.2608
|||1.2.1|119.2608
|||1.2|119.2608<sup>[4](#4)</sup>
|||1.1.0|119.2607
|||1.0.0|119.2607
|||0.5.1|119.2607<sup>[3](#3)</sup>
|||0.5.0|119.2347
|||0.4.0|119.2347
[apbs-smol-surf.in](apbs-smol-surf.in)|2-level focusing to 0.250 A spacing, molecular surface, srfm smol|**1.5**|**108.8770**|86.50
|||1.4.2|108.8770
|||1.4.1|108.8773
|||1.4|108.8773<sup>[5](#5)</sup>
|||1.3|108.8748
|||1.2.1|108.8748
|||1.2|108.8748<sup>[4](#4)</sup>
|||1.1.0|108.8773
|||1.0.0|108.8773
|||0.5.1|108.8773<sup>[3](#3)</sup>
|||0.5.0|108.8540
|||0.4.0|108.8540<sup>[2](#2)</sup>
|||0.3.2|94.8705<sup>[1](#1)</sup>
|||0.3.1|97.0147
|||0.3.0|97.0147
|||0.2.6|97.0147
|||0.2.5|97.0147
|||0.2.4|97.0147
|||0.2.3|97.0147
|||0.2.2|97.0147
|||0.2.1|97.015
|||0.2.0|97.015
|||0.1.8|97.01

<a name=1></a><sup>1</sup> The grid dimensions (dime) changed from 65\^3 to 97\^3 in the 0.3.2 release to give a finer mesh.

<a name=2></a><sup>2</sup> The discrepancy in values between versions 0.4.0 and 0.3.2 is most likely due to three factors:

-   A bug fix in Vacc\_molAcc which removed spurious regions of high internal dielectric values
-   A switch in the algorithm used to compute the dielectric smoothing for srfm smol.
-   The addition of the Vacc sphere density (sdens keyword) as a variable and a change in the default sdens value from 3.0 to 10.0

<a name=3></a><sup>3</sup> The discrepancy in values between versions 0.5.1 and 0.5.0 is most likely due to the following factor(s):

-   Bug fix regarding multipole behavior for neutral proteins

<a name=4></a><sup>4</sup> APBS 1.2 has switched the multigrid smoothing algorithm from standard Gauss-Seidel to Gauss-Seidel red/black in order to facilitate parallelization. This switch has caused small differences in individual calculation energies which, when combined to the final answer, create larger errors (up to 0.04%). These errors can be reduced by resetting the APBS error tolerance to 1e-9 or smaller values. For a more detailed explanation, please see the APBS FAQ, [here](http://www.poissonboltzmann.org/docs/apbs-faq/#sources error calculation).

<a name=5></a><sup>5</sup> The discrepancy in values between versions 1.3 and 1.4 is most likely due to the following factor(s):

-   Translation of contrib/pmgZ library from FORTRAN to C
-   Differences in numerical implementations between FORTRAN and C compilers result in small round-off discrepencies
-   Small margins due to these round-off discrepencies acumulate in the computations

Please see the ChangeLog or the [APBS website](http://www.poissonboltzmann.org/) for more information.


README for protein-RNA APBS example
===================================

This is example is taken directly from

> García-García C, Draper DE. Electrostatic interactions in a peptide-RNA complex. Journal of Molecular Biology. 331 (1), 75-88, 2003. <http://dx.doi.org/10.1016/S0022-2836(03)00615-6>

with the minimized PDB files kindly provided by David Draper. It uses the change in binding free energy with ionic strength to estimate the change in number of ions "bound" (diffusively) to RNA upon protein binding.

A more comprehensive walkthrough of this example is available [online](http://www.poissonboltzmann.org/examples/Protein-Rna_Tutorial/).

The APBS results agree reasonably well with Draper's calculations; there are some differences in parameterization, surface definitions, etc. that lead to small differents in the results.

To run this example, make sure that the apbs and dxmath binaries are in your path; e.g.,

    export PATH=/path/to/apbs/executable:$PATH
    export PATH=/path/to/tools/dxmath/executable:$PATH

<center>
|APBS Version|Input file|APBS Result|
|---|---|---|
|3.0   |apbs-0.025.in|8.674116429351E+01
|      |apbs-0.050.in| 9.606836713867E+01
|      |apbs-0.075.in| 1.011537214883E+02
|      |apbs-0.100.in| 1.046142116108E+02
|      |apbs-0.125.in| 1.072226817611E+02
|      |apbs-0.150.in| 1.093084123761E+02
|      |apbs-0.175.in| 1.110412443877E+02
|      |apbs-0.200.in| 1.125199716537E+02
|      |apbs-0.225.in| 1.138070465620E+02
|      |apbs-0.250.in| 1.149444369078E+02
|      |apbs-0.275.in| 1.159616972338E+02
|      |apbs-0.300.in| 1.168804254687E+02
|      |apbs-0.325.in| 1.177168854907E+02
|      |apbs-0.400.in| 1.198456038802E+02
|      |apbs-0.500.in| 1.220607673699E+02
|      |apbs-0.600.in| 1.238080564885E+02
|      |apbs-0.700.in| 1.252364090877E+02
|      |apbs-0.800.in| 1.264340604647E+02
|1.5   |apbs-0.025.in|8.674120000000e+01
|      |apbs-0.050.in|9.606840000000e+01
|      |apbs-0.075.in|1.011540000000e+02
|      |apbs-0.100.in|1.046140000000e+02
|      |apbs-0.125.in|1.072230000000e+02
|      |apbs-0.150.in|1.093080000000e+02
|      |apbs-0.175.in|1.110410000000e+02
|      |apbs-0.200.in|1.125200000000e+02
|      |apbs-0.225.in|1.138070000000e+02
|      |apbs-0.250.in|1.149440000000e+02
|      |apbs-0.275.in|1.159620000000e+02
|      |apbs-0.300.in|1.168800000000e+02
|      |apbs-0.325.in|1.177170000000e+02
|      |apbs-0.400.in|1.198460000000e+02
|      |apbs-0.500.in|1.220610000000e+02
|      |apbs-0.600.in|1.238080000000e+02
|      |apbs-0.700.in|1.252360000000e+02
|      |apbs-0.800.in|1.264340000000e+02
|1.4.2 |apbs-0.025.in|8.674120000000e+01
|      |apbs-0.050.in|9.606840000000e+01
|      |apbs-0.075.in|1.011540000000e+02
|      |apbs-0.100.in|1.046140000000e+02
|      |apbs-0.125.in|1.072230000000e+02
|      |apbs-0.150.in|1.093080000000e+02
|      |apbs-0.175.in|1.110410000000e+02
|      |apbs-0.200.in|1.125200000000e+02
|      |apbs-0.225.in|1.138070000000e+02
|      |apbs-0.250.in|1.149440000000e+02
|      |apbs-0.275.in|1.159620000000e+02
|      |apbs-0.300.in|1.168800000000e+02
|      |apbs-0.325.in|1.177170000000e+02
|      |apbs-0.400.in|1.198460000000e+02
|      |apbs-0.500.in|1.220610000000e+02
|      |apbs-0.600.in|1.238080000000e+02
|      |apbs-0.700.in|1.252360000000e+02
|      |apbs-0.800.in|1.264340000000e+02
</center>
<!---
%%%%%%
Commented this out since is the result of running all the example. Reather, I am putting the results of the test in here. JB.
%%%%%%
Input file|Description|APBS version|APBS results||Draper PB results||Draper experimental results||
---|---|---|---|---|---|---|---|---
||||n = -(d Δ G)/(d ln [KCl])/RT|-d( Δ G)/(d log<sub>10</sub> [KCl]) (kcal)|n = -(d Δ G)/(d ln [KCl])/RT|-d( Δ G)/(d log<sub>10</sub> [KCl]) (kcal)|n = -(d Δ G)/(d ln [KCl])/RT|-d( Δ G)/(d log<sub>10</sub> [KCl]) (kcal)
'make all'|Run a series of binding energy calculations at different ionic strengths|**1.4.2**|**-(4.52831 ± 0.0758878)**|**6.1561 ± 0.109612**|-(4.3 ± 0.2)|5.9 ± 0.2|-(4.4 ± 0.2)|6.0 ± 0.2
          |                                                                        |1.0.0           |-(4.52 ± 0.08)|6.2 ± 0.1
--->
Please see the ChangeLog or the [APBS website](http://www.poissonboltzmann.org/) for more information.
APBS ion pmf test case
======================

This is the two-ion test case for force evaluation from the Im et al. *Comput. Phys. Commun.* **111**, 59--75 (1998) paper. The output of these run(s) are solvation energies for ion pairs (convert to total energy by adding Coulombic forces) and solvation force components. Molecule 0 is fixed (at x=-3) and molecule 1 was moved.

To run this test case, simply type `./runme.sh`, making sure that the apbs binary is in your path.

Here's the output you can look for (all forces in kJ/mol/Å):

sasa - Solvent-accessible surface area force (scaled by the surface tension,    gamma). **NOTE THAT THIS FORCE HAS CHANGED!** We're now using SASA forces instead of the spline-surface square gradient forces used in previous examples (and by Im et al).

qf - Reaction field force

db - Dielectric boundary force

ib - Ionic boundary force (osmotic pressure)

Here are the force results visually estimated from the Im et al paper (see above) with forces in kJ/mol/Å with one ion fixed at x = -3.0 Å.

Mol 1 location|QF||DB||IB||
---|---|---|---|---|---|---
||Mol 0|Mol 1|Mol 0|Mol 1|Mol 0|Mol 1
-3.00|0.0E+00|0.0E+00|0.0E+00|0.0E+00|0|0
-2.50|8.4E+00|-1.3E+01|1.1E+02|-1.1E+02|0|0
-2.00|3.8E+01|-3.4E+01|8.2E+01|-8.2E+01|0|0
-1.50|6.3E+01|-5.5E+01|4.2E+01|-4.0E+01|0|0
-1.00|8.8E+01|8.8E+01|4.2E+00|-4.2E+01|0|0
-0.50|1.0E+02|-1.0E+02|-1.2E+01|1.6E+01|0|0
0.00|1.1E+02|0|-2.1E+01|2.1E+01|0|0
0.50|9.6E+01|-9.8E+01|-1.1E+01|4.2E+00|0|0
1.00|8.4E+01|-8.8E+01|4.2E+00|0.0E+00|0|0
1.50|6.7E+01|-6.5E+01|6.3E+00|4.2E+00|0|0
2.00|5.3E+01|-5.3E+01|6.3E+01|6.3E+01|0|0
2.50|4.6E+01|-4.8E+01|6.3E+01|-2.1E+01|0|0
3.00|3.8E+01|-3.8E+01|6.3E+00|-6.3E+00|0|0

Here are the results from APBS 3.0

    polarforces
     qf  0  8.398642197666E+01  -1.324564548552E-05  -1.613435632529E-05
     ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
     db  0  6.148357059184E+00  -2.667517425897E-05  -5.378919678211E-05
     tot 0  9.013477903584E+01  -3.992081974449E-05  -6.992355310740E-05
     qf  1  -8.466423642736E+01  -1.836748045969E-05  -3.062224428458E-05
     ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
     db  1  2.882739230548E+00  -3.149946357588E-05  -6.291495506459E-05
     tot 1  -8.178149719681E+01  -4.986694403557E-05  -9.353719934917E-05
     tot all 8.353281839029E+00  -8.978776378007E-05  -1.634607524566E-04

    apolarforces
     sasa  0  -1.099776974333E+01  0.000000000000E+00  0.000000000000E+00
     sav   0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
     wca   0  -5.826577086437E-01  -2.766670515801E-05  -2.766670515838E-05
     tot   0  -1.158042745197E+01  -2.766670515801E-05  -2.766670515838E-05
     sasa  1  1.111862435589E+01  0.000000000000E+00  0.000000000000E+00
     sav   1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
     wca   1  5.826560767576E-01  -2.767485007141E-05  -2.767485007183E-05
     tot   1  1.170128043265E+01  -2.767485007141E-05  -2.767485007183E-05
     tot all  1.208529806779E-01  -5.534155522943E-05  -5.534155523021E-05

Here are the results from APBS 1.5

    x = 1.000

    polarforces
     qf  0  8.398642197666E+01  -1.324564548552E-05  -1.613435632529E-05
     ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
     db  0  6.148357059184E+00  -2.667517425897E-05  -5.378919678211E-05
     tot 0  9.013477903584E+01  -3.992081974449E-05  -6.992355310740E-05
     qf  1  -8.466423642736E+01  -1.836748045969E-05  -3.062224428458E-05
     ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
     db  1  2.882739230548E+00  -3.149946357588E-05  -6.291495506459E-05
     tot 1  -8.178149719681E+01  -4.986694403557E-05  -9.353719934917E-05
     tot all 8.353281839029E+00  -8.978776378007E-05  -1.634607524566E-04

    apolarforces
     sasa  0  -1.099776974333E+01  0.000000000000E+00  0.000000000000E+00
     sav   0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
     wca   0  -5.826577086437E-01  -2.766670515801E-05  -2.766670515838E-05
     tot   0  -1.158042745197E+01  -2.766670515801E-05  -2.766670515838E-05
     sasa  1  1.111862435589E+01  0.000000000000E+00  0.000000000000E+00
     sav   1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
     wca   1  5.826560767576E-01  -2.767485007141E-05  -2.767485007183E-05
     tot   1  1.170128043265E+01  -2.767485007141E-05  -2.767485007183E-05
     tot all  1.208529806779E-01  -5.534155522943E-05  -5.534155523021E-05

Here are the results from APBS 1.4.2:

    x = 1.000

    polarforces
      qf  0  8.398642197666E+01  -1.324564548552E-05  -1.613435632529E-05
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      db  0  6.148357059184E+00  -2.667517425897E-05  -5.378919678211E-05
      tot 0  9.013477903584E+01  -3.992081974449E-05  -6.992355310740E-05
      qf  1  -8.466423642736E+01  -1.836748045969E-05  -3.062224428458E-05
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      db  1  2.882739230548E+00  -3.149946357588E-05  -6.291495506459E-05
      tot 1  -8.178149719681E+01  -4.986694403557E-05  -9.353719934917E-05
      tot all 8.353281839029E+00  -8.978776378007E-05  -1.634607524566E-04

    apolarforces
      sasa  0  -1.099776974333E+01  0.000000000000E+00  0.000000000000E+00
      sav   0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      wca   0  -5.826577086437E-01  -2.766670515801E-05  -2.766670515838E-05
      tot   0  -1.158042745197E+01  -2.766670515801E-05  -2.766670515838E-05
      sasa  1  1.111862435589E+01  0.000000000000E+00  0.000000000000E+00
      sav   1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      wca   1  5.826560767576E-01  -2.767485007141E-05  -2.767485007183E-05
      tot   1  1.170128043265E+01  -2.767485007141E-05  -2.767485007183E-05
      tot all  1.208529806779E-01  -5.534155522943E-05  -5.534155523021E-05

Here are the results from APBS 1.4.1-binary:

    x = 1.000

    polarforces
      qf  0  8.398642197665E+01  -1.324563987321E-05  -1.613435751204E-05
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      db  0  6.148357059184E+00  -2.667517427189E-05  -5.378919673813E-05
      tot 0  9.013477903583E+01  -3.992081414510E-05  -6.992355425017E-05
      qf  1  -8.466423642736E+01  -1.836748529478E-05  -3.062224653038E-05
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      db  1  2.882739230548E+00  -3.149946357422E-05  -6.291495504390E-05
      tot 1  -8.178149719681E+01  -4.986694886900E-05  -9.353720157428E-05
      tot all 8.353281839026E+00  -8.978776301409E-05  -1.634607558245E-04

    apolarforces
      sasa  0  -1.099776974333E+01  0.000000000000E+00  0.000000000000E+00
      sav   0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      wca   0  -5.826577086437E-01  -2.766670515801E-05  -2.766670515838E-05
      tot   0  -1.158042745197E+01  -2.766670515801E-05  -2.766670515838E-05
      sasa  1  1.111862435589E+01  0.000000000000E+00  0.000000000000E+00
      sav   1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      wca   1  5.826560767576E-01  -2.767485007143E-05  -2.767485007183E-05
      tot   1  1.170128043265E+01  -2.767485007143E-05  -2.767485007183E-05
      tot all  1.208529806779E-01  -5.534155522944E-05  -5.534155523021E-05

      Global net ELEC energy = -1.125191605614E+03 kJ/mol
      Global net APOL energy = -1.365308087297E+01 kJ/mol

Here are the results from APBS 1.4:

    x = 1.000

    polarforces
      qf  0   8.398642197666E+01  -1.324564266726E-05  -1.613435892778E-05
      ib  0   0.000000000000E+00   0.000000000000E+00   0.000000000000E+00
      db  0   6.148357059184E+00  -2.667517424343E-05  -5.378919673334E-05
      tot 0   9.013477903584E+01  -3.992081691068E-05  -6.992355566112E-05
      qf  1  -8.466423642735E+01  -1.836748431279E-05  -3.062224860445E-05
      ib  1   0.000000000000E+00   0.000000000000E+00   0.000000000000E+00
      db  1   2.882739230548E+00  -3.149946358193E-05  -6.291495505310E-05
      tot 1  -8.178149719681E+01  -4.986694789471E-05  -9.353720365755E-05
      tot all 8.353281839034E+00  -8.978776480540E-05  -1.634607593187E-04

    apolarforces
      sasa  0  -1.099776974333E+01   0.000000000000E+00   0.000000000000E+00
      sav   0   0.000000000000E+00   0.000000000000E+00   0.000000000000E+00
      wca   0  -5.826577086437E-01  -2.766670515801E-05  -2.766670515838E-05
      tot   0  -1.158042745197E+01  -2.766670515801E-05  -2.766670515838E-05
      sasa  1   1.111862435589E+01   0.000000000000E+00   0.000000000000E+00
      sav   1   0.000000000000E+00   0.000000000000E+00   0.000000000000E+00
      wca   1   5.826560767576E-01  -2.767485007143E-05  -2.767485007183E-05
      tot   1   1.170128043265E+01  -2.767485007143E-05  -2.767485007183E-05
      tot all   1.208529806779E-01  -5.534155522944E-05  -5.534155523021E-05

      Global net ELEC energy = -1.125191605614E+03 kJ/mol

Here are the results from APBS 1.3:

    x = 1.000

    polarforces
      qf  0  8.398640174378E+01  9.818388346045E-06  9.795992548983E-06
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      db  0  6.148361411050E+00  2.004819537885E-05  3.920937677211E-05
      tot 0  9.013476315483E+01  2.986658372490E-05  4.900536932109E-05
      qf  1  -8.466424236101E+01 -9.207441185028E-06 -2.534007494961E-05
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      db  1  2.882747163943E+00  1.731862478386E-05  3.267273837455E-05
      tot 1  -8.178149519706E+01 8.111183598834E-06  7.332663424942E-06
      tot all 8.353267957769E+00 3.797776732373E-05  5.633803274603E-05

    apolarforces
      sasa  0  -1.099776974333E+01  0.000000000000E+00  0.000000000000E+00
      sav   0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      wca   0  -5.826577086437E-01  -2.766670515801E-05  -2.766670515838E-05
      tot   0  -1.158042745197E+01  -2.766670515801E-05  -2.766670515838E-05
      sasa  1  1.111862435589E+01  0.000000000000E+00  0.000000000000E+00
      sav   1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      wca   1  5.826560767576E-01  -2.767485007141E-05  -2.767485007183E-05
      tot   1  1.170128043265E+01  -2.767485007141E-05  -2.767485007183E-05
      tot all  1.208529806779E-01  -5.534155522943E-05  -5.534155523021E-05

      Global net ELEC energy = -1.125191531404E+03 kJ/mol

Here are the results from APBS 1.2.1:

    x = 1.000

    polarforces
      qf  0  8.398640174377E+01  9.818386932509E-06  9.795992245139E-06
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      db  0  6.148361411050E+00  2.004819530432E-05  3.920937680245E-05
      tot 0  9.013476315482E+01  2.986658223683E-05  4.900536904759E-05
      qf  1  -8.466424236099E+01 -9.207442514897E-06 -2.534008066980E-05
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      db  1  2.882747163943E+00  1.731862480409E-05  3.267273838847E-05
      tot 1  -8.178149519705E+01 8.111182289197E-06  7.332657718672E-06
      tot all 8.353267957774E+00 3.797776452602E-05  5.633802676626E-05

    apolarforces
      sasa  0  -1.099776974333E+01  0.000000000000E+00  0.000000000000E+00
      sav   0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      wca   0  -5.826577086437E-01  -2.766670515801E-05  -2.766670515838E-05
      tot   0  -1.158042745197E+01  -2.766670515801E-05  -2.766670515838E-05
      sasa  1  1.111862435589E+01  0.000000000000E+00  0.000000000000E+00
      sav   1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      wca   1  5.826560767576E-01  -2.767485007143E-05  -2.767485007183E-05
      tot   1  1.170128043265E+01  -2.767485007143E-05  -2.767485007183E-05
      tot all  1.208529806779E-01  -5.534155522944E-05  -5.534155523021E-05

      Global net ELEC energy = -1.125191531404E+03 kJ/mol

Here are the results from APBS 1.2:

    x = 1.000

    polarforces
      qf  0  8.398640174377E+01  9.818386932509E-06  9.795992245139E-06
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      db  0  6.148361411050E+00  2.004819530432E-05  3.920937680245E-05
      tot 0  9.013476315482E+01  2.986658223683E-05  4.900536904759E-05
      qf  1  -8.466424236099E+01 -9.207442514897E-06 -2.534008066980E-05
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      db  1  2.882747163943E+00  1.731862480409E-05  3.267273838847E-05
      tot 1  -8.178149519705E+01 8.111182289197E-06  7.332657718672E-06
      tot all 8.353267957774E+00 3.797776452602E-05  5.633802676626E-05

    apolarforces
      sasa  0  -1.099776974333E+01  0.000000000000E+00  0.000000000000E+00
      sav   0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      wca   0  -5.826577086437E-01  -2.766670515801E-05  -2.766670515838E-05
      tot   0  -1.158042745197E+01  -2.766670515801E-05  -2.766670515838E-05
      sasa  1  1.111862435589E+01  0.000000000000E+00  0.000000000000E+00
      sav   1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      wca   1  5.826560767576E-01  -2.767485007143E-05  -2.767485007183E-05
      tot   1  1.170128043265E+01  -2.767485007143E-05  -2.767485007183E-05
      tot all  1.208529806779E-01  -5.534155522944E-05  -5.534155523021E-05

      Global net ELEC energy = -1.125191531404E+03 kJ/mol

Here are the results from APBS 1.1.0:

    x = 1.000

    polarforces
      qf  0  8.398636786524E+01  3.277928951109E-05  2.248755327954E-05
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      db  0  6.148315264050E+00  6.212090029364E-05  9.765261473306E-05
      tot 0  9.013468312929E+01  9.490018980473E-05  1.201401680126E-04
      qf  1  -8.466453877446E+01  1.260694683844E-05  -2.688056989348E-05
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      db  1  2.882714236549E+00  9.382342500375E-05  1.225207937768E-04
      tot 1  -8.178182453791E+01  1.064303718422E-04  9.564022388330E-05
      tot all 8.352858591379E+00  2.013305616469E-04  2.157803918959E-04

    apolarforces
      sasa  0  -1.099776974333E+01  0.000000000000E+00  0.000000000000E+00
      sav   0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      wca   0  -5.826577086437E-01  -2.766670515801E-05  -2.766670515838E-05
      tot   0  -1.158042745197E+01  -2.766670515801E-05  -2.766670515838E-05
      sasa  1  1.111862435589E+01  0.000000000000E+00  0.000000000000E+00
      sav   1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      wca   1  5.826560767576E-01  -2.767485007143E-05  -2.767485007183E-05
      tot   1  1.170128043265E+01  -2.767485007143E-05  -2.767485007183E-05
      tot all  1.208529806779E-01  -5.534155522944E-05  -5.534155523021E-05

      Global net ELEC energy = -1.125192402906E+03 kJ/mol

Here are the results from APBS 0.5.0:

    x = -3.000
      sasa 0 2.014e-01 0.000e+00 0.000e+00
      sasa 1 2.014e-01 0.000e+00 0.000e+00
      sasa  0  2.014243542734E-01  0.000000000000E+00  0.000000000000E+00
      qf  0  2.402694548782E-01  1.193684463475E-04  9.707135192876E-05
      qf  1  2.402694548782E-01  1.193684463475E-04  9.707135192876E-05
      db  0  -1.040061396733E+00  5.059871689052E-05  1.004179762749E-04
      db  1  -1.040061396733E+00  5.059871689052E-05  1.004179762749E-04
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      Global net ELEC energy = -1.494109084160E+03 kJ/mol
      Global net APOL energy = 1.073537020853E-01 kJ/mol
    x = -2.500
      sasa 0 -1.047e+01 0.000e+00 0.000e+00
      sasa 1 1.088e+01 0.000e+00 0.000e+00
      sasa  0  -1.047406642222E+01  0.000000000000E+00  0.000000000000E+00
      qf  0  9.670776686682E+00  3.067299706303E-05  8.525779585775E-06
      qf  1  -9.776393142038E+00  1.160701475502E-04  9.737743502613E-05
      db  0  1.128766544950E+02  4.216378707828E-05  1.156092840465E-04
      db  1  -1.018562119621E+02  1.332582044569E-04  2.007546725043E-04
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      Global net ELEC energy = -1.453172513724E+03 kJ/mol
      Global net APOL energy = 1.127544530116E-01 kJ/mol
    x = -2.000
      sasa 0 -1.088e+01 0.000e+00 0.000e+00
      sasa 1 1.092e+01 0.000e+00 0.000e+00
      sasa  0  -1.087691513076E+01  0.000000000000E+00  0.000000000000E+00
      qf  0  3.600129707630E+01  1.673257382577E-05  -3.172493445854E-06
      qf  1  -3.698431978276E+01  1.362593318370E-04  1.348136933229E-04
      db  0  7.641808295612E+01  8.426812747154E-05  1.501550525155E-04
      db  1  -7.630995607470E+01  1.930747292853E-04  2.577078768374E-04
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      Global net ELEC energy = -1.395864905308E+03 kJ/mol
      Global net APOL energy = 1.202783777713E-01 kJ/mol
    x = -1.50
      sasa 0 -1.055e+01 0.000e+00 0.000e+00
      sasa 1 1.064e+01 0.000e+00 0.000e+00
      sasa  0  -1.055463616393E+01  0.000000000000E+00  0.000000000000E+00
      qf  0  6.263390243427E+01  -1.039272251661E-05  -3.653632695350E-05
      qf  1  -6.274351718559E+01  6.171199346571E-05  3.641783990404E-05
      db  0  4.141531409141E+01  8.939020846417E-05  1.464873997274E-04
      db  1  -3.980350645690E+01  1.618242786613E-04  2.096874862673E-04
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      Global net ELEC energy = -1.343413508063E+03 kJ/mol
      Global net APOL energy = 1.278197055953E-01 kJ/mol
    x = -1.00
      sasa 0 -1.104e+01 0.000e+00 0.000e+00
      sasa 1 1.080e+01 0.000e+00 0.000e+00
      sasa  0  -1.103805461418E+01  0.000000000000E+00  0.000000000000E+00
      qf  0  8.602515079485E+01  1.679958402518E-05  -1.286070388760E-06
      qf  1  -8.589474499290E+01  4.904393981795E-05  2.948346215562E-05
      db  0  5.146620466700E+00  1.263326043428E-04  1.972642536915E-04
      db  1  -2.462461594184E+00  1.419028113702E-04  2.025478366045E-04
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      Global net ELEC energy = -1.297718358292E+03 kJ/mol
      Global net APOL energy = 1.354538497618E-01 kJ/mol
    x = -0.500
      sasa 0 -1.088e+01 0.000e+00 0.000e+00
      sasa 1 1.080e+01 0.000e+00 0.000e+00
      sasa  0  -1.087691513076E+01  0.000000000000E+00  0.000000000000E+00
      qf  0  1.011157426881E+02  3.495539842345E-05  2.650020071179E-05
      qf  1  -1.018411965567E+02  2.379958570171E-05  7.021000073448E-06
      db  0  -1.147809320720E+01  1.908780322530E-04  2.832722451268E-04
      db  1  2.221256509495E+01  1.621554833141E-04  2.462408961005E-04
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      Global net ELEC energy = -1.253921668205E+03 kJ/mol
      Global net APOL energy = 1.431402031209E-01 kJ/mol
    x = 0.000
      sasa 0 -1.047e+01 0.000e+00 0.000e+00
      sasa 1 1.055e+01 0.000e+00 0.000e+00
      sasa  0  -1.047406642222E+01  0.000000000000E+00  0.000000000000E+00
      qf  0  1.076485044405E+02  4.150254188932E-05  2.584327997210E-05
      qf  1  -1.076452387179E+02  1.116736665075E-05  -1.918731685710E-05
      db  0  -2.164756414609E+01  1.200915665417E-04  1.684440777983E-04
      db  1  2.928788039760E+01  1.264642983339E-04  1.736212475138E-04
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      Global net ELEC energy = -1.211966318316E+03 kJ/mol
      Global net APOL energy = 1.508265564799E-01 kJ/mol
    x = 0.500
      sasa 0 -1.072e+01 0.000e+00 0.000e+00
      sasa 1 1.059e+01 0.000e+00 0.000e+00
      sasa  0  -1.071577564735E+01  0.000000000000E+00  0.000000000000E+00
      qf  0  1.027354991386E+02  2.555484206199E-05  1.063932308940E-05
      qf  1  -1.027864740990E+02  4.081203889223E-05  2.562624208050E-05
      db  0  -1.024263914547E+01  6.237801032056E-05  9.677805192977E-05
      db  1  4.694564975253E+00  7.269351372102E-05  1.062464666717E-04
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      Global net ELEC energy = -1.169410509835E+03 kJ/mol
      Global net APOL energy = 1.585419149460E-01 kJ/mol
    x = 1.000
      sasa 0 -1.100e+01 0.000e+00 0.000e+00
      sasa 1 1.112e+01 0.000e+00 0.000e+00
      sasa  0  -1.099776974333E+01  0.000000000000E+00  0.000000000000E+00
      qf  0  8.398636786523E+01  3.277928825388E-05  2.248755650954E-05
      qf  1  -8.466453877447E+01  1.260694492290E-05  -2.688056972615E-05
      db  0  6.148315264050E+00  6.212090033654E-05  9.765261477356E-05
      db  1  2.882714236549E+00  9.382342501820E-05  1.225207938255E-04
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      Global net ELEC energy = -1.125192402906E+03 kJ/mol
      Global net APOL energy = 1.661992631981E-01 kJ/mol
    x = 1.500
      sasa 0 -1.068e+01 0.000e+00 0.000e+00
      sasa 1 1.072e+01 0.000e+00 0.000e+00
      sasa  0  -1.067549077649E+01  0.000000000000E+00  0.000000000000E+00
      qf  0  6.661368503641E+01  3.283120005903E-05  2.336264942404E-05
      qf  1  -6.687857770544E+01  3.794378201242E-07  -4.537653753875E-05
      db  0  7.100467586189E+00  6.377840082367E-05  1.004771191938E-04
      db  1  3.082008432733E+00  1.019851332113E-04  1.318438583362E-04
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      Global net ELEC energy = -1.087409226682E+03 kJ/mol
      Global net APOL energy = 1.739378257498E-01 kJ/mol
    x = 2.000
      sasa 0 -1.092e+01 0.000e+00 0.000e+00
      sasa 1 1.096e+01 0.000e+00 0.000e+00
      sasa  0  -1.091720000162E+01  0.000000000000E+00  0.000000000000E+00
      qf  0  5.395122212118E+01  3.402365136823E-05  2.463252404579E-05
      qf  1  -5.435238595073E+01  2.345850272885E-05  5.498720638168E-06
      db  0  7.099496363901E+00  6.338076881104E-05  9.960971260094E-05
      db  1  1.822551374452E+00  7.422787937764E-05  1.097869885768E-04
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      Global net ELEC energy = -1.056724540553E+03 kJ/mol
      Global net APOL energy = 1.814733525524E-01 kJ/mol
    x = 2.500
      sasa 0 -1.072e+01 0.000e+00 0.000e+00
      sasa 1 1.068e+01 0.000e+00 0.000e+00
      sasa  0  -1.071577564735E+01  0.000000000000E+00  0.000000000000E+00
      qf  0  4.454222472665E+01  3.558462567387E-05  2.660491501258E-05
      qf  1  -4.504040092314E+01  -1.044092908553E-05  -3.622045970037E-05
      db  0  7.053724444849E+00  6.155505900578E-05  9.694085460358E-05
      db  1  -2.309382970894E+00  6.271891825021E-05  9.531078126617E-05
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      Global net ELEC energy = -1.032072744366E+03 kJ/mol
      Global net APOL energy = 1.893047314465E-01 kJ/mol
    x = 3.000
      sasa 0 -1.092e+01 0.000e+00 0.000e+00
      sasa 1 1.051e+01 0.000e+00 0.000e+00
      sasa  0  -1.091720000162E+01  0.000000000000E+00  0.000000000000E+00
      qf  0  3.737527031761E+01  3.343902445084E-05  2.403559379423E-05
      qf  1  -3.737539015853E+01  4.414820479892E-05  3.441749305654E-05
      db  0  7.008341960925E+00  6.385969683875E-05  1.003016820956E-04
      db  1  -7.008312180867E+00  7.315905619890E-05  1.101572582696E-04
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      Global net ELEC energy = -1.011312542835E+03 kJ/mol
      Global net APOL energy = 1.969446766344E-01 kJ/mol
    x = 3.500
      sasa 0 -1.064e+01 0.000e+00 0.000e+00
      sasa 1 1.068e+01 0.000e+00 0.000e+00
      sasa  0  -1.063520590564E+01  0.000000000000E+00  0.000000000000E+00
      qf  0  3.179181222228E+01  3.472564738899E-05  2.566914250463E-05
      qf  1  -3.259379322155E+01  1.216664204331E-05  -2.497181119826E-05
      db  0  6.969782295514E+00  6.257908926389E-05  9.840967875397E-05
      db  1  7.904107869827E+00  8.728726024045E-05  1.158752820808E-04
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      Global net ELEC energy = -9.936580311110E+02 kJ/mol
      Global net APOL energy = 2.046252289721E-01 kJ/mol
    x = 4.000
      sasa 0 0.000e+00 0.000e+00 0.000e+00
      sasa 1 0.000e+00 0.000e+00 0.000e+00
      sasa  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      qf  0  2.736006693935E+01  3.271235623929E-05  2.341325695619E-05
      qf  1  -2.776033811952E+01  8.396915386866E-06  -3.668786316582E-05
      db  0  6.937279812724E+00  6.412225118915E-05  1.008451446546E-04
      db  1  1.489454481196E+00  1.060308624524E-04  1.379944549983E-04
      ib  0  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      ib  1  0.000000000000E+00  0.000000000000E+00  0.000000000000E+00
      Global net ELEC energy = -9.786249660218E+02 kJ/mol
      Global net APOL energy = 2.091848317949E-01 kJ/mol
README for ion-protein APBS examples
====================================

The example input files in this directory calculate the energy of placing an ion near a macromolecule.

This example was contributed by Dave Sept.

Input File|Description|APBS Version|Results (kJ/mol)|UHBD (kJ/mol)
---|---|---|---|---|
[apbs-mol-pdiel2.in](apbs-mol-pdiel2.in)|0.53 A resolution, pdie 2, srfm mol|**3.0**|**15.5916**|**23.58**
|||1.5|15.5916
|||1.4.2|15.5916
|||1.4.1|15.5916
|||1.4|15.5916
|||1.3|15.5916<sup>[3](#3)</sup>
|||1.2.1|15.5819
|||1.2|15.5819
|||1.1.0|15.5819<sup>[2](#2)</sup>
|||1.0.0|15.5916
|||0.5.1|15.5916
|||0.5.0|15.5916
|||0.4.0|15.5916
[apbs-smol-pdiel2.in](apbs-smol-pdiel2.in)|0.53 A resolution, pdie 2, srfm smol|**3.0**|**23.5554**|**23.58**
|||1.5|23.5554
|||1.4.2|23.5554
|||1.4.1|23.5554
|||1.4|23.5554
|||1.3|23.5554<sup>[3](#3)</sup>
|||1.2.1|23.5458
|||1.2|23.5458
|||1.1.0|23.5458<sup>[2](#2)</sup>
|||1.0.0|23.5554
|||0.5.1|23.5554
|||0.5.0|23.5554
|||0.4.0|23.5554<sup>[1](#1)</sup>
|||0.3.2|21.4763
|||0.3.1|19.8794
|||0.3.0|19.8794
|||0.2.6|19.8794
|||0.2.5|19.8794
|||0.2.4|19.8794
|||0.2.3|19.8652
|||0.2.2|21.4530
|||0.2.1|21.453
|||0.2.0|21.453
|||0.1.8|21.45
[apbs-mol-pdiel12.in](apbs-mol-pdiel12.in)|0.53 A resolution, pdie 12, srfm mol|**3.0**|**18.0272**|**23.58**
|||1.5|18.0272
|||1.4.2|18.0272
|||1.4.1|18.0272
|||1.4|18.0272
|||1.3|18.0272<sup>[3](#3)</sup>
|||1.2.1|18.0176
|||1.2|18.0176
|||1.1.0|18.0176<sup>[2](#2)</sup>
|||1.0.0|18.0272
|||0.5.1|18.0272
|||0.5.0|18.0272
|||0.4.0|18.0272
[apbs-smol-pdiel12.in](apbs-smol-pdiel12.in)|0.53 A resolution, pdie 12, srfm smol|**3.0**|**19.2825**|**23.58**
|||1.5|19.2825
|||1.4.2|19.2825
|||1.4.1|19.2825
|||1.4|19.2825
|||1.3|19.2825<sup>[3](#3)</sup>
|||1.2.1|19.2728
|||1.2|19.2728
|||1.1.0|19.2728<sup>[2](#2)</sup>
|||1.0.0|19.2825
|||0.5.1|19.2825
|||0.5.0|19.2825
|||0.4.0|19.2825<sup>[1](#1)</sup>
|||0.3.2|21.4763
|||0.3.1|18.9205
|||0.3.0|17.4207
|||0.2.6|17.4207
|||0.2.5|17.4207
|||0.2.4|17.4207
|||0.2.3|17.4049
|||0.2.2|18.8953
|||0.2.1|18.895
|||0.2.0|18.895
|||0.1.8|18.90

<a name=1></a><sup>1</sup> The discrepancy in values between versions 0.4.0 and 0.3.2 is most likely due to three factors:

-   A bug fix in Vacc\_molAcc which removed spurious regions of high internal dielectric values
-   A switch in the algorithm used to compute the dielectric smoothing for srfm smol.
-   The addition of the Vacc sphere density (sdens keyword) as a variable and a change in the default sdens value from 3.0 to 10.0

<a name=2></a><sup>2</sup> The discrepancy in values between versions 1.0.0 and 1.1.0 is due to a bugfix in the implementation of the boundary conditions. This bug introduces a very small error (generally less than 1%) the calculated results. This error is most prominent when the molecule substantially overlaps the boundary (e.g., in the current example) and is often symptomatic of insufficiently-large problem domains.

<a name=3></a><sup>3</sup> The discrepancy in values between versions 1.2.1 and 1.3 is most likely due to the following factor(s):

-   Fixed a bug in Vpmg.c which causes zero potential values on boundaries in non-focusing calculations.

Please see the ChangeLog or the [APBS website](http://www.poissonboltzmann.org/) for more information.
README for geometric flow APBS examples
=============================

The example input files in this directory use a differential
geometry-based geometric flow solvation model to describes a smooth
dielectric interface profile across the solvent-solute boundary in a
thermodynamically self-consistent fashion.  The main parameters of the
model are the solute/solvent dielectric coefficients, solvent pressure on
the solute, microscopic surface tension, solvent density, and molecular
force-field parameters.

More information can be found here in [Chen et al.](http://www.ncbi.nlm.nih.gov/pubmed/20938489) and [Thomas et
al](http://www.ncbi.nlm.nih.gov/pubmed/23212974).

The parameters used in the input files (imidazol and glycerol) came from Thomas et al.

| Input File  | Description                 | APBS Version     | Global Net ELEC Energy | Global Net APOL Energy |
|-------------|-----------------------------|------------------|------------------------|------------------------|
| imidazol.in | Serial, mdh boundary condts.| **3.0** | **-10.3022** | **0.541742** |
|||1.5|-10.3022|0.541742
|||1.4.2|-10.3022|0.541742
README for TABI-PB (Boundary Element Method) pKa Examples
=========================================================

The example input files included in this folder uses a boundary element approach called
TABI-PB to solve the PBE. BEMs have the characteristic that only the boundary of the
domain has to be discretized. This is particularly useful for problems in which the data
of interest is at the boundary of the solution.

This directory contains five example .in files:
        1. 2LZT-ASH66.in
        2. 2LZT-ASP66.in
        3. 2LZT-noASH66.in
        4. 2LZT-noASP66.in
        5. ASH66.in
        6. ASP66.in

These files can be used to demonstrate TABI-PB pKa calculations, documented at:
http://www.poissonboltzmann.org/examples/Lysozyme_pKa_example/

File Input| APBS Version| Result (kJ/mol)| Expected (kJ/mol)
---|---|---|---
[ASH66.in](ASH66.in)| **3.0**| **-4.026**| **-4.026**
[2LTZ-ASH66.in](2LTZ-ASH66.in)| **3.0**| **-369.020**| **-369.020**
[2LTZ-noASH66.in](2LTZ-noASH66.in)| **3.0**| **-367.795**| **-367.795**
README for solv APBS examples
=============================

The example input files in this directory calculate the solvation energies of methanol and methoxide.

The source for this example is UHBD.


Input File                          | Description | APBS Version | Methanol Results (kJ/mol) | Methoxide Results (kJ/mol) |Difference (kJ/mol)| Methanol UHBD (kJ/mol)|Methoxide UHBD (kJ/mol)|Difference UHBD (kJ/mol)
---------------------------|---------|-------------|----------|----|-----------|------|----------|--------------
[apbs-mol.in](apbs-mol.in)| Focusing to 0.25 A, srfm mol| **3.0** | **-36.2486** | **-390.4122**|**-354.1635**|**-35.595**|**-390.023**|**-354.424**
| | |1.5   |-36.2486|-390.4120|-354.1640
| | |1.4.2 |-36.2486|-390.4120|-354.1640
| | |1.4.1 |-36.2486|-390.4122|-354.1635
| | |1.4   |-36.2486|-390.4122<sup>[3](#3)</sup>|-354.1635
| | |1.3 | -36.2486| -390.4121| -354.1635
| | |1.2.1 | -36.2486| -390.4121| -354.1635
| | |1.2<sup>[2](#2)</sup> |-36.2486| -390.4121| -354.1635
| | |1.1.0 |-36.2486| -390.4119	| -354.1632
| | |1.0.0 |-36.2486| -390.4119	| -354.1632
| | |0.5.1 |-36.2486| -390.4119	| -354.1632
| | |0.5.0 |-36.2486| -390.4119	| -354.1632
| | |0.4.0 |-36.2486| -390.4119	| -354.1632
[apbs-smol.in](apbs-smol.in) | Focusing to 0.25 A, srfm smol | **3.0** | **-37.5759** | **-391.2388**| **-353.6629**| **-35.595**| **-390.023**| **-354.424**
| | |1.5|-37.5759|-391.2390|-353.6630
| | |1.4.2|-37.5759|-391.2390|-353.6630
| | |1.4.1|-37.5759|-391.2389|-353.6629
| | |1.4 |-37.5759<sup>[3](#3)</sup>|-391.2389<sup>[3](#3)</sup>| -353.6629
| | |1.3 |-37.5760| -391.2388| -353.6629
| | |1.2.1 |-37.5760| -391.2388| -353.6629
| | |1.2<sup>[2](#2)</sup> |-37.5760| -391.2388| -353.6629
| | |1.1.0 |-37.5760| -391.2388| -353.6627
| | |1.0.0 |-37.5760| -391.2388| -353.6627
| | |0.5.1 |-37.5760| -391.2388| -353.6627
| | |0.5.0 |-37.5760| -391.2388| -353.6627
| | |0.4.0<sup>[1](#1)</sup> |-37.5760| -391.2388| -353.6627
| | |0.3.2 |-36.2486| -390.4119| -354.1632
| | |0.3.1 |-36.2486| -390.4119| -354.1632
| | |0.3.0 |-36.2486| -390.4119| -354.1632
| | |0.2.6 |-36.2486| -390.4119| -354.1632
| | |0.2.5 |-36.2486| -390.4119| -354.1632
| | |0.2.4 |-36.2486| -390.4119| -354.1632
| | |0.2.3 |-36.2223| -391.7995| -355.5771
| | |0.2.2 |-36.2223| -391.7995| -355.5771
| | |0.2.1 |-36.222| -391.800| -355.577
| | |0.2.0 |-36.222| -391.800| -355.577
| | |0.1.8 |-36.222| -391.800| -355.577




<a name=1></a><sup>1</sup> The discrepancy in values between versions 0.4.0 and 0.3.2 is most likely due to three factors:
A bug fix in Vacc_molAcc which removed spurious regions of high internal dielectric values
A switch in the algorithm used to compute the dielectric smoothing for srfm smol
The addition of the Vacc sphere density (sdens keyword) as a variable and a change in the default sdens value from 3.0 to 10.0

<a name=2></a><sup>2</sup> APBS 1.2 has switched the multigrid smoothing algorithm from standard Gauss-Seidel to Gauss-Seidel red/black in order to facilitate parallelization. This switch has caused small differences in individual calculation energies which, when combined to the final answer, create larger errors (up to 0.04%). These errors can be reduced by resetting the APBS error tolerance to 1e-9 or smaller values. For a more detailed explanation, please see here.

<a name=3></a><sup>3</sup> The discrepancy in values between versions 1.3 and 1.4 is most likely due to the following factor(s):

Translation of contrib/pmgZ library from FORTRAN to C
Differences in numerical implementations between FORTRAN and C compilers result in small round-off discrepencies
Small margins due to these round-off discrepencies acumulate in the computations

Please see the ChangeLog or the [APBS website](http://www.poissonboltzmann.org/) for more information.
README for PBAM APBS examples
=============================

The example input files in this directory uses an analytical 
solution to the linear Poisson-Boltzmann equation.
The main parameters of the model are the solute/solvent dielectric 
coefficients, system temperature and salt concentration.

More information can be found here in [Cite!]() and 
[Lotan and Head-Gordon](http://pubs.acs.org/doi/full/10.1021/ct050263p).

**toy_energyforce.in**|**APBS Version**|**Result Energy Mol 1**|**Expected Energy Mol 1**|**Result Force Mol 1**|**Expected Force Mol 1**
---|---|---|---|---|---
---|1.5|-4.965018552290E+01|-4.965019E+01|(-7.018534E-16, -1.074613E+01, -5.196457E-16)|(-7.018534E-16, -1.074613E+01, -5.196457E-16)

**toy_energyforce.in**|**APBS Version**|**Result Energy Mol 2**|**Expected Energy Mol 2**|**Result Force Mol 2**|**Expected Force Mol 2**
---|---|---|---|---|---
---|1.5|-4.965018552290E+01|-4.965019E+01|(7.018535E-16, 1.074613E+01, 5.137566E-16)|(7.018535E-16, 1.074613E+01, 5.137566E-16)

README for FKBP APBS examples
=============================

The example input files in this directory simulate the binding of various (small) ligands to FKBP. Analogous to HCA binding case (except it works).

In order to calculate solvation energy upon binding you will need to take the results from these input files and subtract from them the results obtained from the `coulomb` utility found at `apbs/tools/manip/coulomb`. The values returned from this utility are:

-   1d7h-dmso: -15.0930 kJ/mol (analytical value -15.103 kJ/mol)
-   1d7i-dss: -11.9670 kJ/mol (analytical value -11.975 kJ/mol)

This example was contributed by Jung-Hsin Lin.

Input File|Description|APBS Version|Results (kJ/mol)|UHBD (kJ/mol)
---|---|---|---|---
[1d7h-dmso/apbs-mol.in](1d7h-dmso/apbs-mol.in)|1d7h-dmso, 2-level focusing to 0.225 A, VdW surface, srfm mol|**3.0**|**15.0081**|**19.097**
|||1.5|15.0081
|||1.4.2|15.0081
|||1.4.1|15.0081
|||1.4|15.0081<sup>[4](#4)</sup>
|||1.3|15.0077
|||1.2.1|15.0077<sup>[3](#3)</sup>
|||1.2|15.0087<sup>[2](#2)</sup>
|||1.1.0|15.0089
|||1.0.0|15.0089
|||0.5.1|15.0089
|||0.5.0|15.0089
|||0.4.0|15.0089
[1d7h-dmso/apbs-smol.in](1d7h-dmso/apbs-smol.in)|1d7h-dmso, 2-level focusing to 0.225 A, VdW surface, srfm smol|**3.0**|**16.2445**|**19.097**
|||1.5|16.2445
|||1.4.2|16.2445
|||1.4.1|16.2445
|||1.4|16.2445<sup>[4](#4)</sup>
|||1.3|16.2446
|||1.2.1|16.2446<sup>[3](#3)</sup>
|||1.2|16.2456<sup>[2](#2)</sup>
|||1.1.0|16.2458
|||1.0.0|16.2458
|||0.5.1|16.2458
|||0.5.0|16.2458
|||0.4.0|16.2458<sup>[1](#1)</sup>
|||0.3.2|15.0089
|||0.3.1|15.0089
|||0.3.0|15.0089
|||0.2.6|15.0089
|||0.2.5|15.0089
|||0.2.4|15.0089
|||0.2.3|15.0097
|||0.2.2|14.5886
|||0.2.1|14.589
|||0.2.0|14.589
|||0.1.8|14.591
[1d7i-dss/apbs-mol.in](1d7i-dss/apbs-mol.in)|1d7i-dss, 2-level focusing to 0.225 A, VdW surface, srfm mol|**3.0**|**14.4250**|**16.231**
|||1.5|14.4250
|||1.4.2|14.4250
|||1.4.1|14.4250
|||1.4|14.4250
|||1.3|14.4250
|||1.2.1|14.4250<sup>[3](#3)</sup>
|||1.2|14.4253<sup>[2](#2)</sup>
|||1.1.0|14.4254
|||1.0.0|14.4254
|||0.5.1|14.4254
|||0.5.0|14.4254
|||0.4.0|14.4254
[1d7i-dss/apbs-smol.in](1d7i-dss/apbs-smol.in)|1d7i-dss, 2-level focusing to 0.225 A, VdW surface, srfm smol|**3.0**|**15.4515**|**16.231**
|||1.5|15.4515
|||1.4.2|15.4515
|||1.4.1|15.4515
|||1.4|15.4515
|||1.3|15.4515
|||1.2.1|15.4515<sup>[3](#3)</sup>
|||1.2|15.4517
|||1.1.0|15.4517
|||1.0.0|15.4517
|||0.5.1|15.4517
|||0.5.0|15.4517
|||0.4.0|15.4517<sup>[1](#1)</sup>
|||0.3.2|14.4254
|||0.3.1|14.4254
|||0.3.0|14.4254
|||0.2.6|14.4254
|||0.2.5|14.4254
|||0.2.4|14.4254
|||0.2.3|14.4254
|||0.2.2|14.3865
|||0.2.1|14.387
|||0.2.0|14.387
|||0.1.8|15.210

<a name=1></a><sup>1</sup> The discrepancy in values between versions 0.4.0 and 0.3.2 is most likely due to three factors:

-   A bug fix in Vacc\_molAcc which removed spurious regions of high internal dielectric values
-   A switch in the algorithm used to compute the dielectric smoothing for srfm smol
-   The addition of the Vacc sphere density (sdens keyword) as a variable and a change in the default sdens value from 3.0 to 10.0

<a name=2></a><sup>2</sup> APBS 1.2 has switched the multigrid smoothing algorithm from standard Gauss-Seidel to Gauss-Seidel red/black in order to facilitate parallelization. This switch has caused small differences in individual calculation energies which, when combined to the final answer, create larger errors (up to 0.04%). These errors can be reduced by resetting the APBS error tolerance to 1e-9 or smaller values. For a more detailed explanation, please see the APBS FAQ, [here](http://www.poissonboltzmann.org/docs/apbs-faq/#sources error calculation).

<a name=3></a><sup>3</sup> The discrepancy in values between versions 1.2 and 1.2.1 is most likely due to the following factor(s):

-   Fixed a bug in Vpmg\_fillcoCoefMolIon which causes npbe based calculations to return very large energies

<a name=4></a><sup>4</sup> The discrepancy in values between versions 1.3 and 1.4 is most likely due to the following factor(s):

-   Translation of contrib/pmgZ library from FORTRAN to C
-   Differences in numerical implementations between FORTRAN and C compilers result in small round-off discrepencies
-   Small margins due to these round-off discrepencies acumulate in the computations

Please see the ChangeLog or the [APBS website](http://www.poissonboltzmann.org/) for more information.
README for TABI-PB (Boundary Element Method) Examples
=====================================================

The example input files included in this folder uses a boundary element approach called
TABI-PB to solve the PBE. BEMs have the characteristic that only the boundary of the
domain has to be discretized. This is particularly useful for problems in which the data
of interest is at the boundary of the solution.

This directory contains four example .in files.

Two examples calculate surface potentials for 1a63 in a 0.15 M salt solution:
        1. 1a63_NanoShaper_SES.in uses NanoShaper to create an SES triangulation.
        2. 1a63_NanoShaper_Skin.in uses NanoShaper to create a Skin surface triangulation.

Two examples calculate surface potentials for 451c in a 0.15 M salt solution:
        1. 451c_order1.in uses a 1st order Taylor series expansion for the treecode.
        2. 451c_order5.in uses a 5th order Taylor series expansion for the treecode.

Additionally, more pqr files are available in the test_proteins directory.
More details on TABI-PB are available on the APBS website contributions section.

Input File| APBS Version| Result (kJ/mol) | Expected (kJ/mol)
---|---|---|---
[451c_order1.in](451c_order1.in)| **3.0**| **-4907.441**| -4907.443
| | 1.5| -4907.455|
[451c_order5.in](451c_order5.in)| **3.0**| **-4920.112**| -4920.133
| | 1.5| -4920.091|
README for PBSAM-gly Examples
=============================

This method calculates the energies and forces for inter-molecule interactions. For more information please visit the [PB-SAM](http://www.poissonboltzmann.org/external_contributions/extern-pbsam) website.

**gly_energyforce.in**|**APBS Version**|**Result Mol 1**|**Expected Mol 1**|**Result Force Mol 1**|**Expected Force Mol 1**
---|---|---|---|---|---
---|1.5|6.176608555569E-05|6.176609E-05|(-5.376354E-04, -4.238467E-04, -2.339671E-04)|(-5.376354E-04, -4.238467E-04, -2.339671E-04)

**gly_energyforce.in**|**APBS Version**|**Result Mol 2**|**Expected Mol 2**|**Result Force Mol 2**|**Expected Force Mol 2**
---|---|---|---|---|---
---|1.5|6.210593742215E-05|6.210594E-05|(5.351512E-04, 4.459663E-04, 2.411459E-04)|(5.351512E-04, 4.459663E-04, 2.411459E-04)README for point-pmf APBS examples
==================================

The example input files in this directory see how well we do reproducing Coulomb's law.

Be sure to edit the runme.sh script to point to the APBS binary.

This example was contributed by Nathan Baker.

Input File|Description|APBS Version|Results (kJ/mol)||||Analytical (kJ/mol)||||
---|---|---|---|---|---|---|---|---|---|---
||||1 A Dist|2 A Dist|3 A Dist|4 A Dist|1 A Dist|2 A Dist|3 A Dist|4 A Dist
[apbs.in](apbs.in)|Focusing to 0.21 A, srfm spl2|**1.5**|**18.3082**|**8.90669**|**5.9096**|**4.43014**|17.686|8.843|5.89533|4.4215
|||1.4.2|18.3082|8.90669|5.9096|4.43014
|||1.4.1|18.3082|8.9067|5.9096|4.4301
|||1.4|18.3082|8.9067|5.9096|4.4301
|||1.3|18.3082|8.9067|5.9096|4.4301
|||1.2.1|18.3082|8.9067|5.9096|4.4301
|||1.2|18.3082|8.9067|5.9096|4.4301
|||1.1.0|18.3082|8.9067|5.9096|4.4301
|||1.0.0|18.3082|8.9067|5.9096|4.4301
|||0.5.1|18.3082|8.9067|5.9096|4.4301
|||0.5.0|18.3082|8.9067|5.9096|4.4301
|||0.4.0|18.3082|8.9067|5.9096|4.4301
|||0.3.2|18.3082|8.9067|5.9096|4.4301
|||0.3.1|18.3082|8.9067|5.9096|4.4301
|||0.3.0|18.3082|8.9067|5.9096|4.4301
|||0.2.6|18.3082|8.9067|5.9096|4.4301
|||0.2.5|18.3082|8.9067|5.9096|4.4301
|||0.2.4|18.3082|8.9067|5.9096|4.4301
|||0.2.3|18.3082|8.9067|5.9096|4.4301
|||0.2.2|18.3082|8.9067|5.9096|4.4301
|||0.2.1|18.308|8.907|5.910|4.430
|||0.2.0|18.31|8.91|5.91|4.43
|||0.1.8|18.308|8.907|5.910|4.430

README for ionize APBS examples
===============================

The example input files in this directory calculate electrostatic contributions to the ionization energy of acetic acid.

Input File|Description|APBS Version|Results (kJ/mol)||||UHBD (kJ/mol)||||
---|---|---|---|---|---|---|---|---|---|---
||||Acetic Acid|Acetate|Proton|Ionization Energy|Acetic Acid|Acetate|Proton|Ionization Energy
[apbs-mol.in](apbs-mol.in)|3-level focusing to 0.094 A, srfm mol|**3.0**|**-22.6788**|**-199.746**|**-297.46**|**-474.527**|**-22.22**|**-198.04**|**-295.79**|**-471.61**
|||1.5|-22.6788|-199.746|-297.46|-474.527|
|||1.4.2|-22.6788|-199.746|-297.46|-474.527
|||1.4.1|-22.6788|-199.7463|-297.4598|-474.5273
|||1.4|-22.6788|-199.7463|-297.4598|-474.5273
|||1.3|-22.6788|-199.7463|-297.4598|-474.5273
|||1.2.1|-22.6788|-199.7463|-297.4598|-474.5273
|||1.2<sup>[1](#1)</sup>|-22.6788|-199.7463|-297.4598|-474.5273
|||1.1.0|-22.6787|-199.7462|-297.4599|-474.5273
|||1.0.0|-22.6787|-199.7462|-297.4599|-474.5273
|||0.5.1|-22.6787|-199.7462|-297.4599|-474.5273
|||0.5.0|-22.6787|-199.7462|-297.4599|-474.5273
|||0.4.0|-22.6787|-199.7462|-297.4599|-474.5273
|||0.3.2|-22.6787|-199.7462|-297.4599|-474.5273
|||0.3.1|-22.6787|-199.7462|-297.4599|-474.5273
|||0.3.0|-22.6787|-199.7462|-297.4599|-474.5273
|||0.2.6|-22.6787|-199.7462|-297.4599|-474.5273
|||0.2.5|-22.6787|-199.7462|-297.4599|-474.5273
|||0.2.4|-22.6787|-199.7462|-297.4599|-474.5273
|||0.2.3|-22.629|-195.3135|-290.8751|-463.5617
|||0.2.2|-22.628|-195.3051|-290.8591|-463.5373
|||0.2.1|-22.627|-195.305|-290.859|-463.537
|||0.2.0|-22.627|-195.305|-290.859|-463.537
|||0.1.8|-22.63|-195.31|-290.86|-463.54
[apbs-smol.in](apbs-smol.in)|3-level focusing to 0.094 A, srfm smol|**3.0**|**-22.3305**|**-198.488**|**-295.967**|**-472.125**|**-22.22**|**-198.04**|**-295.79**|**-471.61**
|||1.5|-22.3305|-198.488|-295.967|-472.125
|||1.4.2|-22.3305|-198.488|-295.967|-472.125
|||1.4.1|-22.3305|-198.4883|-295.9669|-472.1247
|||1.4|-22.3305|-198.4883|-295.9669|-472.1247
|||1.3|-22.3305|-198.4883|-295.9669|-472.1247
|||1.2.1|-22.3305|-198.4883|-295.9669|-472.1247
|||1.2<sup>[1](#1)</sup>|-22.3305|-198.4883|-295.9669|-472.1247
|||1.1.0|-22.3304|-198.4881|-295.9670|-472.1247
|||1.0.0|-22.3304|-198.4881|-295.9670|-472.1247
|||0.5.1|-22.3304|-198.4881|-295.9670|-472.1247
|||0.5.0|-22.3304|-198.4881|-295.9670|-472.1247
|||0.4.0|-22.3304|-198.4881|-295.9670|-472.1247

<a name=1></a><sup>1</sup> APBS 1.2 has switched the multigrid smoothing algorithm from standard Gauss-Seidel to Gauss-Seidel red/black in order to facilitate parallelization. This switch has caused small differences in individual calculation energies which, when combined to the final answer, create larger errors (up to 0.04%). These errors can be reduced by resetting the APBS error tolerance to 1e-9 or smaller values. For a more detailed explanation, please see the APBS FAQ, [here](http://www.poissonboltzmann.org/docs/apbs-faq/#sources error calculation).

Please see the ChangeLog or the [APBS website](http://www.poissonboltzmann.org/) for more information.
README for Actin-Dimer APBS examples
====================================

The example input files in this directory calculate binding energies for actin dimers. This is an example of a large biomolecule binding energy calculation that often requires parallel focusing.

This example was contributed by Dave Sept.

Input File                          | Description | APBS Version | Results (kJ/mol) | UHBD (kJ/mol)
------------------------------------|-------------|--------------|------------------|--------------
[apbs-mol-auto.in](apbs-mol-auto.in)| Sequential, 2-level focusing to ≤ 0.725 A, NPBE, srfm mol| **3.0** | **104.868** | 106.7 (1.00 A res., NPBE)
| | |1.5 | 104.868 |
| | |1.4.2 |104.868
| | |1.4.1 |104.8683
| | |1.4   |104.8683
| | |1.3 | 104.8683<sup>[8](#8)</sup>
| | |1.2.1 | 104.867
| | |1.2 |104.867
| | |1.1.0 |104.867<sup>[5](#5)</sup>
| | |1.0.0 |104.868
| | |0.5.1 |104.868<sup>[3](#3)</sup>
| | |0.5.0 | 105.0338<sup>[2](#2)</sup>
| | |0.4.0 |104.8895
[apbs-smol-auto.in](apbs-smol-auto.in) | Sequential, 2-level focusing to ≤ 0.725 A, NPBE, srfm smol | **3.0** | **109.5841** | 106.7 (1.00 A res., NPBE)
| | |1.5   |109.5841
| | |1.4.2 |109.5841
| | |1.4.1 |109.5841
| | |1.4   |109.5841
| | |1.3   |109.5841<sup>[8](#8)</sup>
| | |1.2.1 |109.5829
| | |1.2   |109.5829
| | |1.1.0 |109.5829<sup>[5](#5)</sup>
| | |1.0.0 |109.5841
| | |0.5.1 |109.5841<sup>[3](#3)</sup>
| | |0.5.0 |109.7518<sup>[2](#2)</sup>
| | |0.4.0 |109.6043<sup>[1](#1)</sup>
| | |0.3.2 |90.8704
| | |0.3.1 |88.6101
| | |0.3.0 |88.6101
| | |0.2.6 |88.6101
| | |0.2.5 |88.6101
| | |0.2.4 |88.6101
| | |0.2.3 |88.6064
| | |0.2.2 |90.829
| | |0.2.1 |90.829
| | |0.2.0 |90.829
| | |0.1.8 |90.84
[apbs-mol-parallel.in](apbs-mol-parallel.in) |Parallel with 8 processors, focusing to \~0.9 A, LPBE, srfm mol |**1.5**|**98.1746**|106.7 (1.00 A res., NPBE)
| | |1.4.2 |98.1746
| | |1.4.1 |98.1746
| | |1.4   |98.1746
| | |1.3   |98.1746<sup>[8](#8)</sup>
| | |1.2.1 |98.1733<sup>[7](#7)</sup>
| | |1.2   |98.1635<sup>[6](#6)</sup>
| | |1.1.0 |98.1630<sup>[5](#5)</sup>
| | |1.0.0 |98.1643<sup>[4](#4)</sup>
| | |0.5.1 |98.1654<sup>[3](#3)</sup>
| | |0.5.0 |98.3530<sup>[2](#2)</sup>
| | |0.4.0 |98.1834
[apbs-smol-parallel.in](apbs-smol-parallel.in)|Parallel with 8 processors, focusing to \~0.9 A, LPBE, srfm smol|**1.5**|**115.542**|106.7 (1.00 A res., NPBE)
| | |1.4.2 |115.542
| | |1.4.1 |115.5421
| | |1.4   |115.5421<sup>[9](#9)</sup>
| | |1.3   |115.5422<sup>[8](#8)</sup>
| | |1.2.1 |115.5409<sup>[7](#7)</sup>
| | |1.2   |115.5563<sup>[6](#6)</sup>
| | |1.1.0 |115.5560<sup>[5](#5)</sup>
| | |1.0.0 |115.5573<sup>[4](#4)</sup>
| | |0.5.1 |115.5584<sup>[3](#3)</sup>
| | |0.5.0 |115.7492<sup>[2](#2)</sup>
| | |0.4.0 |115.5751<sup>[1](#1)</sup>
| | |0.3.2 |87.1121
| | |0.3.1 |87.1121
| | |0.3.0 |90.2573
| | |0.2.6 |90.2573
| | |0.2.5 |90.2573
| | |0.2.4 |90.2573
| | |0.2.3 |90.2543
| | |0.2.2 |91.9450
| | |0.2.1 |91.945
| | |0.2.0 |91.939
| | |0.1.8 |91.67

<a name=1></a><sup>1</sup> The discrepancy in values between versions 0.4.0 and 0.3.2 is most likely due to three factors:

-   A bug fix in Vacc\_molAcc which removed spurious regions of high internal dielectric values
-   A switch in the algorithm used to compute the dielectric smoothing for srfm smol.
-   The addition of the Vacc sphere density (sdens keyword) as a variable and a change in the default sdens value from 3.0 to 10.0

<a name=2></a><sup>2</sup> The discrepancy in values between versions 0.5.0 and 0.4.0 is most likely due to the following factor(s):

-   A change in the autofocusing routine for APBS

<a name=3></a><sup>3</sup> The discrepancy in values between versions 0.5.1 and 0.5.0 is most likely due to the following factor(s):

-   Bug fix regarding multipole behavior for neutral proteins

<a name=4></a><sup>4</sup> The discrepancy in values between versions 0.5.1 and 1.0.0 was due to the execution of the previous APBS tests on a PowerPC platform with the XLC/XLF compilers. Running with binaries compiled with gcc/gfortran or the Intel compilers gives identical results between versions 0.5.1 and 1.0.0.

<a name=5></a><sup>5</sup> The discrepancy in values between versions 1.0.0 and 1.1.0 is due to a bugfix in the implementation of the boundary conditions. This bug introduces a very small error (generally less than 1%) the calculated results.

<a name=6></a><sup>6</sup> APBS 1.2 has switched the multigrid smoothing algorithm from standard Gauss-Seidel to Gauss-Seidel red/black in order to facilitate parallelization. This switch has caused small differences in individual calculation energies which, when combined to the final answer, create larger errors (up to 0.04%). These errors can be reduced by resetting the APBS error tolerance to 1e-9 or smaller values. For a more detailed explanation, please see the FAQ, [here](http://www.poissonboltzmann.org/docs/apbs-faq/#sources error calculation).

<a name=7></a><sup>7</sup> The discrepancy in values between versions 1.2 and 1.2.1 is most likely due to the following factor(s):

-   Fixed a bug in Vpmg\_fillcoCoefMolIon which causes npbe based calculations to return very large energies

<a name=8></a><sup>8</sup> The discrepancy in values between versions 1.2.1 and 1.3 is most likely due to the following factor(s):

-   Fixed a bug in Vpmg.c which causes zero potential values on boundaries in non-focusing calculations.

<a name=9></a><sup>9</sup> The discrepancy in values between versions 1.3 and 1.4 is most likely due to the following factor(s):

-   Translation of contrib/pmgZ library from FORTRAN to C
-   Differences in numerical implementations between FORTRAN and C compilers result in small round-off discrepencies
-   Small margins due to these round-off discrepencies acumulate in the computations

Please see the ChangeLog or the [APBS website](http://www.poissonboltzmann.org/) for more information.
README file for the Misc folder.
==================================================================

A collection of PQR files of molecules with interesting potentials

|                        |                                      |
|------------------------|--------------------------------------|
| [fas2.pqr](fas2.pqr)   | Fasciculin-2 (1FSC)                  |
| [mache.pqr](mache.pqr) | Mouse acetylcholinesterase (1MAH)    |
| [achbp.pqr](achbp.pqr) | Acetylcholine binding protein (1I9B) |

README for hca-bind APBS examples
=================================

The example input files in this directory calculate the binding of a small molecule (acetazolamide) to a medium-sized protein (human carbonic anhydrase).

The UHBD calculations where performed using a van der Waals surface definition for the dielectric coefficient. This is simulated in the APBS input files by setting srad to 0.0.

Input File|Description|APBS Version|Results (kJ/mol)|UHBD (kJ/mol)
---|---|---|---|---
[apbs-mol.in](apbs-mol.in)|2-level focusing to 0.225 A, VdW surface, srfm mol|**3.0**|**-52.4648**|**-70.00**
|||1.5|-52.4648
|||1.4.2|-52.4648
|||1.4.1|-52.4648<sup>[6](#6)</sup>
|||1.4|-51.4648<sup>[5](#5)</sup>
|||1.3|-52.4647
|||1.2.1|-52.4647
|||1.2|-52.4647<sup>[4](#4)</sup>
|||1.1.0|-52.4669
|||1.0.0|-52.4669
|||0.5.1|-52.4669<sup>[3](#3)</sup>
|||0.5.0|-52.1062<sup>[2](#2)</sup>
|||0.4.0|-52.4414
[apbs-smol.in](apbs-smol.in)|2-level focusing to 0.225 A, VdW surface, srfm smol|**3.0**|**-54.0598**|**-70.00**
|||1.5|-54.0598
|||1.4.2|-54.0598
|||1.4.1|-54.0598
|||1.3|-54.0598
|||1.2.1|-54.0598
|||1.2|-54.0598<sup>[4](#4)</sup>
|||1.1.0|-54.0587
|||1.0.0|-54.0587
|||0.5.1|-54.0587<sup>[3](#3)</sup>
|||0.5.0|-54.7039<sup>[2](#2)</sup>
|||0.4.0|-54.0393<sup>[1](#1)</sup>
|||0.3.2|-57.1192
|||0.3.1|-57.1192
|||0.3.0|-57.1192
|||0.2.6|-57.1192
|||0.2.5|-57.1192
|||0.2.4|-57.1192
|||0.2.3|-57.1123
|||0.2.2|-57.1123
|||0.2.1|-57.112
|||0.2.0|-57.711
|||0.1.8|-58.51

<a name=1></a><sup>1</sup> The discrepancy in values between versions 0.4.0 and 0.3.2 is most likely due to three factors:

-   A bug fix in Vacc\_molAcc which removed spurious regions of high internal dielectric values
-   A switch in the algorithm used to compute the dielectric smoothing for srfm smol
-   The addition of the Vacc sphere density (sdens keyword) as a variable and a change in the default sdens value from 3.0 to 10.0

<a name=2></a><sup>2</sup> The discrepancy in values between versions 0.5.0 and 0.4.0 is most likely due to the following factor(s):

-   A change in the autofocusing routine for APBS

<a name=3></a><sup>3</sup> The discrepancy in values between versions 0.5.1 and 0.5.0 is most likely due to the following factor(s):

-   Bug fix regarding multipole behavior for neutral proteins

<a name=4></a><sup>4</sup> APBS 1.2 has switched the multigrid smoothing algorithm from standard Gauss-Seidel to Gauss-Seidel red/black in order to facilitate parallelization. This switch has caused small differences in individual calculation energies which, when combined to the final answer, create larger errors (up to 0.04%). These errors can be reduced by resetting the APBS error tolerance to 1e-9 or smaller values. For a more detailed explanation, please see the APBS FAQ, [here](http://www.poissonboltzmann.org/docs/apbs-faq/#sources error calculation).

<a name=5></a><sup>5</sup> The discrepancy in values between versions 1.3 and 1.4 is most likely due to the following factor(s):

-   Translation of contrib/pmgZ library from FORTRAN to C
-   Differences in numerical implementations between FORTRAN and C compilers result in small round-off discrepencies
-   Small margins due to these round-off discrepencies acumulate in the computations

<a name=6></a><sup>6</sup> The discrepancy in the result between versions 1.4 and 1.4.1-binary is most likely due to a reporting error.

Please see the ChangeLog or the [APBS website](http://www.poissonboltzmann.org/) for more information.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Smartphone (please complete the following information):**
 - Device: [e.g. iPhone6]
 - OS: [e.g. iOS8.1]
 - Browser [e.g. stock browser, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
.. _GitHub issues: https://github.com/Electrostatics/apbs/issues

============
Getting help
============

-------------
GitHub issues
-------------

Our preferred mechanism for user questions and feedback is via `GitHub issues`_.
We monitor these issues daily and usually respond within a few days.

-------------
Announcements
-------------

Announcements about updates to the APBS-PDB2PQR software and related news are available through our `mailing list <http://us11.campaign-archive1.com/home/?u=a5808042b2b3ea90ee3603cd8&id=28701e36f0>`_; please `register for updates <http://eepurl.com/by4eQr>`_.

----------------------
Contacting the authors
----------------------

If all else fails, feel free to contact nathanandrewbaker@gmail.com.
===============
Further reading
===============

-------------------------
General solvation reviews
-------------------------

* Baker NA. Poisson-Boltzmann methods for biomolecular electrostatics. Methods in Enzymology, 383, 94-118, 2004. http://www.sciencedirect.com/science/article/pii/S0076687904830052
* Baker NA, McCammon JA. Electrostatic interactions. Structural Bioinformatics. Weissig H, Bourne PE, eds., 2005. http://dx.doi.org/10.1002/0471721204.ch21
* Baker NA. Biomolecular applications of Poisson-Boltzmann methods. Reviews in Computational Chemistry. Lipkowitz KB, Larter R, Cundari TR., 21, 2005. http://dx.doi.org/10.1002/0471720895.ch5
* Baker NA. Improving implicit solvent simulations: a Poisson-centric view. Curr Opin Struct Biol, 15, 137-43, 2005. http://dx.doi.org/10.1016/j.sbi.2005.02.001
* Baker NA, Bashford D, Case DA. Implicit solvent electrostatics in biomolecular simulation. New Algorithms for Macromolecular Simulation. Leimkuhler B, Chipot C, Elber R, Laaksonen A, Mark A, Schlick T, Schutte C, Skeel R, eds., 2006. http://dx.doi.org/10.1007/3-540-31618-3_15
* Dong F, Olsen B, Baker NA. Computational Methods for Biomolecular Electrostatics. Methods in Cell Biology: Biophysical Tools for Biologists, 84, 843-870, 2008. http://www.sciencedirect.com/science/article/pii/S0091679X0784026X
* Ren P, Chun J, Thomas DG, Schnieders MJ, Marucho M, Zhang J, Baker NA. Biomolecular electrostatics and solvation: a computational perspective. Quart Rev Biophys, 45 (4), 427-491, 2012. http://dx.doi.org/10.1017/S003358351200011X

---------------------
APBS parallel solvers
---------------------

* Baker NA, Sept D, Joseph S, Holst MJ, McCammon JA. Electrostatics of nanosystems: application to microtubules and the ribosome. Proc Natl Acad Sci USA, 98, 10037-41, 2001. http://dx.doi.org/10.1073/pnas.181342398
* Baker NA, Sept D, Holst MJ, McCammon JA. The adaptive multilevel finite element solution of the Poisson-Boltzmann equation on massively parallel computers. IBM J Res Devel, 45, 427-38, 2001. http://dx.doi.org/10.1147/rd.453.0427

---------------------
APBS multigrid solver
---------------------

* M. Holst, Adaptive numerical treatment of elliptic systems on manifolds. Advances in Computational Mathematics 15, 139-191, 2001 http://dx.doi.org/10.1023/A:1014246117321
* M. Holst and F. Saied, Numerical solution of the nonlinear Poisson-Boltzmann equation: Developing more robust and efficient methods. J. Comput. Chem. 16, 337-364, 1995.
* M. Holst and F. Saied, Multigrid solution of the Poisson-Boltzmann equation. J. Comput. Chem. 14, 105-113, 1993.

--------------------------
APBS finite element solver
--------------------------

* Holst M, Baker NA, Wang F. Adaptive multilevel finite element solution of the Poisson-Boltzmann equation I: algorithms and examples. J Comput Chem, 21, 1319-42, 2000. http://bit.ly/1goFAFE
* Baker N, Holst M, Wang F. Adaptive multilevel finite element solution of the Poisson-Boltzmann equation II: refinement schemes based on solvent accessible surfaces. J Comput Chem, 21, 1343-52, 2000. http://bit.ly/1dNSP8l

--------------------------
APBS geometric flow solver
--------------------------

* Chen Z, Baker NA, Wei GW. Differential geometry based solvation model I: Eulerian formulation, J Comput Phys, 229, 8231-58, 2010. http://dx.doi.org/10.1016/j.jcp.2010.06.036
* Chen Z, Baker NA, Wei GW. Differential geometry based solvation model II: Lagrangian formulation. J Math Biol, 63, 1139-1200, 2011. http://dx.doi.org/10.1007/s00285-011-0402-z
* Chen Z, Zhao S, Chun J, Thomas DG, Baker NA, Wei GW. Variational approach for nonpolar solvation analysis. Journal of Chemical Physics, 137, 084101, 2012. http://dx.doi.org/10.1063/1.4745084
* Thomas DG, Chun J, Chen Z, Wei G, Baker NA. Parameterization of a Geometric flow implicit solvation model. Journal of Computational Chemistry, 34, 687-95, 2013. http://dx.doi.org/10.1002/jcc.23181
* Daily M, Chun J, Heredia-Langner A, Baker NA. Origin of parameter degeneracy and molecular shape relationships in geometric-flow calculations of solvation free energies. J Chem Phys, 139, 204108, 2013. http://dx.doi.org/10.1063/1.4832900

-------------------------------
TABI-PB boundary element solver
-------------------------------

* Geng W, Krasny R. A treecode-accelerated boundary integral Poisson–Boltzmann solver for electrostatics of solvated biomolecules, J Comput Phys, 247, 62-78, 2013. https://doi.org/10.1016/j.jcp.2013.03.056

-----------------------------------------------------------
Structural bioinformatics based on electrostatic properties
-----------------------------------------------------------

* Zhang X, Bajaj CL, Kwon B, Dolinsky TJ, Nielsen JE, Baker NA. Application of new multi-resolution methods for the comparison of biomolecular electrostatic properties in the absence of global structural similarity. Multiscale Model Simul, 5, 1196-213, 2006. http://dx.doi.org/10.1137/050647670
* Chakraborty S, Rao BJ, Baker N, Ásgeirsson B. Structural phylogeny by profile extraction and multiple superimposition using electrostatic congruence as a discriminator. Intrinsically Disordered Proteins, 1 (1), e25463, 2013. https://www.landesbioscience.com/journals/idp/article/25463/

-------------------
Other fun with APBS
-------------------

* Wagoner JA, Baker NA. Assessing implicit models for nonpolar mean solvation forces: the importance of dispersion and volume terms. Proc Natl Acad Sci USA, 103, 8331-6, 2006. http://dx.doi.org/10.1073/pnas.0600118103
* Swanson JMJ, Wagoner JA, Baker NA, McCammon JA. Optimizing the Poisson dielectric boundary with explicit solvent forces and energies: lessons learned with atom-centered dielectric functions. J Chem Theory Comput, 3, 170-83, 2007. http://dx.doi.org/10.1021/ct600216k
* Schnieders MJ, Baker NA, Ren P, Ponder JW. Polarizable Atomic Multipole Solutes in a Poisson-Boltzmann Continuum. J Chem Phys, 126, 124114, 2007. http://dx.doi.org/10.1063/1.2714528
* Callenberg KM, Choudhary OP, de Forest GL, Gohara DW, Baker NA, Grabe M. APBSmem: A graphical interface for electrostatic calculations at the membrane. PLoS ONE, 5, e12722, 2010. http://dx.doi.org/10.1371/journal.pone.0012722
* Unni S, Huang Y, Hanson RM, Tobias M, Krishnan S, Li WW, Nielsen JE, Baker NA. Web servers and services for electrostatics calculations with APBS and PDB2PQR. J Comput Chem, 32 (7), 1488-1491, 2011. http://dx.doi.org/10.1002/jcc.21720
* Konecny R, Baker NA, McCammon JA. iAPBS: a programming interface to the adaptive Poisson–Boltzmann solver. Computational Science and Discovery, 5, 015005, 2012. http://dx.doi.org/10.1088/1749-4699/5/1/015005
* Jurrus E, Engel D, Star K, Monson K, Brandi J, Felberg LE, Brookes DH, Wilson L, Chen J, Liles K, Chun M, Li P, Gohara DW, Dolinsky T, Konecny R, Koes DR, Nielsen JE, Head-Gordon T, Geng W, Krasny R, Wei G-W, Holst MJ, McCammon JA, Baker NA. Improvements to the APBS biomolecular solvation software suite. Protein Sci, 27 (1), 112-128, 2018. https://doi.org/10.1002/pro.3280
* Laureanti J, Brandi J, Offor E, Engel D, Rallo R, Ginovska B, Martinez X, Baaden M, Baker NA. Visualizing biomolecular electrostatics in virtual reality with UnityMol‐APBS. Protein Sci, 29 (1), 237-246, 2020. https://doi.org/10.1002/pro.3773.. APBS documentation master file, created by
   sphinx-quickstart on Fri Jul 17 09:09:16 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

APBS - Adaptive Poisson-Boltzmann Solver
========================================

:Release: |release|
:Date: |today|

========
Overview
========

An understanding of electrostatic interactions is essential for the study of biomolecular processes.
The structures of proteins and other biopolymers are being determined at an increasing rate through structural genomics and other efforts while specific linkages of these biopolymers in cellular pathways or supramolecular assemblages are being detected by genetic and proteomic studies.
To integrate this information in physical models for drug discovery or other applications requires the ability to evaluate the energetic interactions within and between biopolymers.
Among the various components of molecular energetics, solvation properties and electrostatic interactions are of special importance due to the long range of these interactions and the substantial charges of typical biopolymer components.

APBS solves the equations of continuum electrostatics for large biomolecular assemblages.
This software was designed “from the ground up” using modern design principles to ensure its ability to interface with other computational packages and evolve as methods and applications change over time.
The APBS code is accompanied by extensive documentation for both users and programmers and is supported by a variety of utilities for preparing calculations and analyzing results.
Finally, the free, open-source APBS license ensures its accessibility to the entire biomedical community.

========
Contents
========

.. toctree::
   :maxdepth: 2

   getting/index
   using/index
   background
   supporting
   help
   reading
   formats/index
   api/index
   releases
   todo

------------------
Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
===============
Supporting APBS
===============

--------------------------
Please register as a user!
--------------------------

Please help ensure continued support for APBS-PDB2PQR by `registering your use of our software <http://eepurl.com/by4eQr>`_.

-------------------
Citing our software
-------------------

If you use APBS in your research, please cite one or more of the following papers:

* Jurrus E, Engel D, Star K, Monson K, Brandi J, Felberg LE, Brookes DH, Wilson L, Chen J, Liles K, Chun M, Li P, Gohara DW, Dolinsky T, Konecny R, Koes DR, Nielsen JE, Head-Gordon T, Geng W, Krasny R, Wei G-W, Holst MJ, McCammon JA, Baker NA. Improvements to the APBS biomolecular solvation software suite. Protein Sci, 27 (1), 112-128, 2018. https://doi.org/10.1002/pro.3280
* R. Bank and M. Holst, A New Paradigm for Parallel Adaptive Meshing Algorithms. SIAM Review 45, 291-323, 2003. http://epubs.siam.org/doi/abs/10.1137/S003614450342061
* Baker NA, Sept D, Joseph S, Holst MJ, McCammon JA. Electrostatics of nanosystems: application to microtubules and the ribosome. Proc. Natl. Acad. Sci. USA 98, 10037-10041 2001. http://www.pnas.org/content/98/18/10037
* M. Holst, Adaptive numerical treatment of elliptic systems on manifolds. Advances in Computational Mathematics 15, 139-191, 2001 http://dx.doi.org/10.1023/A:1014246117321
* M. Holst and F. Saied, Numerical solution of the nonlinear Poisson-Boltzmann equation: Developing more robust and efficient methods. J. Comput. Chem. 16, 337-364, 1995.
* M. Holst and F. Saied, Multigrid solution of the Poisson-Boltzmann equation. J. Comput. Chem. 14, 105-113, 1993.

------------------------
Supporting organizations
------------------------

The PDB2PQR authors would like to give special thanks to the supporting organizations behind the APBS and PDB2PQR software:

`National Institutes of Health <http://www.nih.gov>`_
  Primary source of funding for APBS via grant GM069702
`National Biomedical Computation Resource <http://nbcr.ucsd.edu/>`_
  Deployment and computational resources support from 2002 to 2020
`National Partnership for Advanced Computational Infrastructure <http://www.hipersoft.rice.edu/npaci/>`_
  Funding and computational resources
`Washington University in St. Louis <http://biochem.wustl.edu>`_
  Start-up fundingSolvation model background
==========================

----------------
Solvation models
----------------

Electrostatic and solvation models can be roughly divided into two classes ([Warshel2006]_, [Roux1999]_, [Ren2012]_) explicit solvent models that treat the solvent in atomic detail and implicit solvent models that generally replace the explicit solvent with a dielectric continuum.
Each method has its strengths and weaknesses.
While explicit solvent models offer some of the highest levels of detail, they generally require extensive sampling to converge properties of interest.
On the other hand, implicit solvent models trade detail and some accuracy for the “pre-equilibration” of solvent degrees of freedom and elimination of sampling for these degrees of freedom. Implicit solvent methods are popular for a variety of biomedical research problems.

The polar solvation energy is generally associated with a difference in charging free energies in vacuum and solvent.
A variety of implicit solvent models are available to biomedical researchers to describe polar solvation; however, the most widely-used methods are currently the Generalized Born (GB) and Poisson-Boltzmann (PB) models.
GB and related methods are very fast heuristic models for estimating the polar solvation energies of biomolecular structures and therefore are often used in high-throughput applications such as molecular dynamics simulations.
PB methods can be formally derived from more detailed theories and offer a somewhat slower, but often more accurate, method for evaluating polar solvation properties and often serve as the basis for parameterization and testing of GB methods.
Finally, unlike most GB methods, PB models provide a global solution for the electrostatic potential and field within and around a biomolecule, therefore making them uniquely suited to visualization and other structural analyses, diffusion simulations, and a number of other methods which require global electrostatic properties.

The PB equation ([Fogolari2002]_, [Lamm2003]_, [Grochowski2007]_, [Baker2005]_) is a nonlinear elliptic partial differential equation of the form shown in the figure below which is solved for the electrostatic potential.
The coefficients of this equation are directly related to the molecular structure of the system under consideration.
PB theory is approximate and, as a result, has several well-known limitations which can affect its accuracy ([Grochowski2007]_, [Netz2000]_), particularly for strongly charged systems or high salt concentrations.
However, despite these limitations, PB methods are still very important for biomolecular structural analysis, modeling, and simulation.
Furthermore, these limitations are currently being addressed through new implicit solvent models and hybrid treatments which extend the applicability of PB theory while preserving some of its computational efficiency.
There are currently examples of both types of treatments which leverage APBS ([Azuara2006]_, [Chu2007]_, [Vitalis2004]_).

.. image:: /media/pb-schematic.png

PB methods provide polar solvation energies and therefore must be complemented by non-polar solvation models to provide a complete view of biomolecular solvent-solute interactions. non-polar solvation is generally associated with the insertion of the uncharged solute into solvent. There are many non-polar solvation models available; however, work by Levy et al. [Levy2003]_ as well as our own research [Wagoner2006]_ has demonstrated the importance of non-polar implicit solvent models which include treatment of attractive solute-solvent dispersion terms.
This model has been implemented in APBS and can also be easily transformed into simpler popular non-polar models (e.g., solvent-accessible surface area).
While this model can be used separately from PB to analyze non-polar contributions to solvation energy, its preferred use is coupled to the PB equation through a geometric flow model [Chen2010]_ which treats polar and non-polar interactions in the same framework and reduces the number of user-specified empirical parameters.

.. _errors:

----------------------------
Caveats and sources of error
----------------------------

^^^^^^^^^^^
Model error
^^^^^^^^^^^

When performing solvation calculations using APBS, it is important to keep in mind that you are using an approximate model for solvation.
Therefore, your answers may contain errors related to approximations in the model.
Many review articles have covered the nature of these approximations, we will stress the highlights below.

""""""""""""""""""""""""""
Linear dielectric response
""""""""""""""""""""""""""

The Poisson-Boltzmann equation models the solvent as a dielectric continuum that responds linearly to all applied fields.
In particular, under this model, very strong fields can induce unrealistically strong polarization in the dielectric media representing the solvent and/or the solute interior.
However, molecular solvents or solutes cannot support an infinite amount of polarization: they are limited by their density, their finite dipole moments, and their finite degree of electronic polarizability.
Therefore, the continuum model assumption of linear dielectric response can break down in situations with strong electric fields; e.g., around nucleic acids or very highly-charged proteins.

"""""""""""""""""""""""""
Local dielectric response
"""""""""""""""""""""""""

The Poisson-Boltzmann equation models the solvent as a dielectric continuum that also responds locally to all applied fields. 
n other words, under this model, the local polarization at a point x is only dependent on the field at point x.
However, molecular solvents and solutes clearly don't obey this assumption: the variety of covalent, steric, and other non-bonded intra- and inter-molecular interactions ensures that the polarization at point x is dependent on solute-field interactions in a non-vanishing neighborhood around x.
One way to limit the impact of this flawed assumption, is to model solute response as "explicitly" as possible in your continuum electrostatics problems.
In other words, rather than relying upon the continuum model to reproduce conformational relaxation or response in your solute, model such response in detail through molecular simulations or other conformational sampling.

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Ambiguity of dielectric interfaces and coefficient values
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Violation of the assumptions of linear and local dielectric response in real molecular systems leads to serious ambiguity in the definition of the dielectric coefficient in the Poisson-Boltzmann equation.
In particular, while the values for bulk solvent (i.e., far away from the solute) response are well-defined, all other values of the dielectric coefficient are ambiguous.
In general, continuum models assume a constant low-dielectric value inside the solute and the bulk solvent value outside the solute.
This assumption creates tremendous sensitivity of calculation results on the placement of the dielectric interface (usually determined by solute atomic radii) and the specific value of the internal solute dielectric.
In general, errors arising from this assumption can be minimized by using internal dielectric values that are consistent with the solute atomic radii parameterization.

""""""""""""""""""""""""""""""""""""""""""""""""""
No specific ion-solvent or ion-solute interactions
""""""""""""""""""""""""""""""""""""""""""""""""""

Most Poisson-Boltzmann models assume that ions do not interact directly with the solvent: they are charges embedded in the same dielectric material as the bulk solvent.
This assumption implies that ions experience no "desolvation" penalty as they interact with the solute surface.
Additionally, most Poisson-Boltzmann models assume that ions interaction with the solute only through electrostatic and hard-sphere steric potentials.
However, this assumption neglects some of the subtlety of ion-protein interactions; in particular, dispersive interactions that can possibly lead to some degree of ion specificity.

"""""""""""""""""""""""
Mean field ion behavior
"""""""""""""""""""""""

Finally, the Poisson-Boltzmann model is a "mean field" description of ionic solutions.
This means that ions only experience the average influence of other ions in the system; the model neglects fluctuations in the ionic atmosphere and correlations between the ions in solution.
Such correlations and fluctuations can be very important at high ionic charge densities; e.g., for multivalent ions, high ion concentrations, or the high-density ionic regions near highly-charged biomolecules.

^^^^^^^^^^^^^^^^^^^^
Parameter set errors
^^^^^^^^^^^^^^^^^^^^

.. todo::

   Under construction; please see https://arxiv.org/abs/1705.10035 for an initial discussion.
   Saved as issue https://github.com/Electrostatics/apbs/issues/481 

^^^^^^^^^^^^^^^^^^^^^^
Structure-based errors
^^^^^^^^^^^^^^^^^^^^^^

Electrostatics calculations can be very sensitive to errors in the structure, including:

* Misplaced atoms or sidechains

* Missing regions of biomolecular structure

* Incorrect titration state assignments

Of these errors, incorrect titration states are the most common and, often, the most problematic.
The software package PDB2PQR was created to minimize all of the above problems and we recommend its use to "pre-process" structures before electrostatics calculations.

^^^^^^^^^^^^^^^^^^^^
Discretization error
^^^^^^^^^^^^^^^^^^^^

The Poisson-Boltzmann partial differential equation must be discretized in order to be solved on a computer.
APBS discretizes the equation in spacing by evaluating the problem coefficients and solving for the electrostatic potential on a set of grid (finite difference) or mesh (finite element) points.
However, this discretization is an approximation to the actual, continuously-specified problem coefficients.
Coarser discretization of coefficients and the solution reduce the overall accuracy and introduce errors into the final potential and calculated energies.

It is very important to evaluate the sensitivity of your calculated energies to the grid spacings and lengths.
In general, it is a good idea to scan a range of grid spacings and lengths before starting a problem and choose the largest problem domain with the smallest grid spacing that gives consistent results (e.g., results that don't change as you further reduce the grid spacing).

^^^^^^^^^^^^^^^^^^^^^^^^^^
Solver and round-off error
^^^^^^^^^^^^^^^^^^^^^^^^^^

APBS uses iterative solvers to solve the nonlinear algebraic equations resulting from the discretized Poisson-Boltzmann equation.
Iterative solvers obtain solutions to algebraic equations which are accurate within a specified error tolerance.
Current versions of APBS use a fixed error tolerance of 10\ :sup:`-6` which implies approximately 1 part per million root-mean-squared error in calculated potentials.
Such error tolerances have been empirically observed to give good accuracy in the calculated energies obtained with APBS. 

However, it is important to note that the error in potential does not necessarily directly relate to the error in the energies calculated by APBS.
In particular, most meaningful energies are calculated as differences between energies from several calculations.
While the accuracy of each separate energy can be related to the solver error tolerance, the energy difference can only be loosely bounded by the error tolerance.

This issue is illustrated in the protein kinase ligand binding example provided with APBS as ``pka-lig`` and analyzed below.
This example demonstrates that, while the errors for each calculation remain small, the overall error in the computed energy can be very large; particularly when two different methods are compared.

.. list-table:: Sensitivity of PB energies to iterative solver error tolerance (APBS 1.2)
   :header-rows: 1

   * - Error tolerance
     - Protein energy
     - Protein energy relative error (with respect to 10\ :sup:`-12` tolerance)
     - Ligand energy
     - Ligand energy relative error (with respect to 10\ :sup:`-12` tolerance)
     - Complex energy
     - Complex energy relative error (with respect to 10\ :sup:`-12` tolerance)
     - Binding energy
     - Binding energy relative error (with respect to 10\ :sup:`-12` tolerance)
   * - 1.00E-06
     - 3.01E+05
     - 2.47E-08
     - 1.05E+04
     - 1.42E-08
     - 3.11E+05
     - 2.45E-08
     - 8.08E+00
     - 7.75E-06
   * - 1.00E-09
     - 3.01E+05
     - 3.19E-11
     - 1.05E+04
     - 1.71E-11
     - 3.11E+05
     - 2.45E-08
     - 8.08E+00
     - 2.48E-09
   * - 1.00E-12
     - 3.01E+05
     - 0.00E+00
     - 1.05E+04
     - 0.00E+00
     - 3.11E+05
     - 0.00E+00
     - 8.08E+00
     - 0.00E+00

---------------
Further reading
---------------

.. [Azuara2006] Azuara C, Lindahl E, Koehl P, Orland H, and Delarue M, PDB_Hydro: incorporating dipolar solvents with variable density in the Poisson-Boltzmann treatment of macromolecule electrostatics. Nucleic Acids Research, 2006. 34: p. W38-W42.

.. [Baker2005] Baker NA, Biomolecular Applications of Poisson-Boltzmann Methods, in Reviews in Computational Chemistry, Lipkowitz KB, Larter R, and Cundari TR, Editors. 2005, John Wiley and Sons.

.. [Chen2010] Chen Z, Baker NA, Wei GW. Differential geometry based solvation model I: Eulerian formulation, J Comput Phys, 229, 8231-58, 2010.

.. [Chu2007] Chu VB, Bai Y, Lipfert J, Herschlag D, and Doniach S, Evaluation of Ion Binding to DNA Duplexes Using a Size-Modified Poisson-Boltzmann Theory. Biophysical Journal, 2007. 93(9): p. 3202-9.

.. [Fogolari2002] Fogolari F, Brigo A, and Molinari H, The Poisson-Boltzmann equation for biomolecular electrostatics: a tool for structural biology. Journal of Molecular Recognition, 2002. 15(6): p. 377-92.

.. [Grochowski2007] Grochowski P, lstrok A, and Trylska J, Continuum molecular electrostatics, salt effects and counterion binding. A review of the Poisson-Boltzmann theory and its modifications. Biopolymers, 2007. 89(2): p. 93-113.

.. [Lamm2003] Lamm G, The Poisson-Boltzmann Equation, in Reviews in Computational Chemistry, Lipkowitz KB, Larter R, and Cundari TR, Editors. 2003, John Wiley and Sons, Inc. p. 147-366.

.. [Levy2003] Levy RM, Zhang LY, Gallicchio E, and Felts AK, On the nonpolar hydration free energy of proteins: surface area and continuum solvent models for the solute-solvent interaction energy. Journal of the American Chemical Society, 2003. 125(31): p. 9523-30.

.. [Netz2000] Netz RR and Orland H, Beyond Poisson-Boltzmann: Fluctuation effects and correlation functions. European Physical Journal E, 2000. 1(2-3): p. 203-14.

.. [Ren2012] Ren P, Chun J, Thomas DG, Schnieders M, Marucho M, Zhang J, Baker NA, Biomolecular electrostatics and solvation: a computational perspective. Quarterly Reviews of Biophysics, 2012. 45(4): p. 427-491.

.. [Roux1999] Roux B and Simonson T, Implicit solvent models. Biophysical Chemistry, 1999. 78(1-2): p. 1-20.

.. [Vitalis2004] Vitalis A, Baker NA, McCammon JA, ISIM: A program for grand canonical Monte Carlo simulations of the ionic environment of biomolecules, Molecular Simulation, 2004, 30(1), 45-61.

.. [Wagoner2006] Wagoner JA and Baker NA, Assessing implicit models for nonpolar mean solvation forces: the importance of dispersion and volume terms. Proceedings of the National Academy of Sciences of the United States of America, 2006. 103(22): p. 8331-6.

.. [Warshel2006] Warshel A, Sharma PK, Kato M, and Parson WW, Modeling electrostatic effects in proteins. Biochimica et Biophysica Acta (BBA) - Proteins & Proteomics, 2006. 1764(11): p. 1647-76.

===============
Release history
===============


---------------------
APBS X.X.X ()
---------------------

* Binary releases may be found in `GitHub releases <https://github.com/Electrostatics/apbs/releases>`_.

^^^^^^^^^^^^
New Features
^^^^^^^^^^^^

* 

^^^^^^^^^^^^^^^^^^^^^^^^^^
Known Bugs and Limitations
^^^^^^^^^^^^^^^^^^^^^^^^^^

(see v3.4.0)

* 

^^^^^^^^^^^^^
Minor Updates
^^^^^^^^^^^^^

* Updated Windows build to Windows 2022; Visual Studio updated to v17, 2022
* Specified Ubuntu build (20.04) and Mac build (11)
* Fixed CMake syntax error around Python detection
* Statically link Python in the Mac and Ubuntu builds

^^^^^
Notes
^^^^^

* The following are included in APBS:

  * `Geometric Flow <https://github.com/Electrostatics/geoflow_c/tree/39d53269c084f1dc1caa71de95dca77f19da739e>`_
  * `FETk <https://github.com/Electrostatics/FETK/tree/8c2b67fe587336ba73f77573f13e31ecb1a5a7f9>`_
  * `PBAM/PBSAM <https://github.com/Electrostatics/pb_solvers/tree/d3ba994d7ec2b2cad5b3e843784c7cb9f41ace37>`_
  * `TABI-PB <https://github.com/Treecodes/TABI-PB/tree/fe1c237b057418fed48535db125394607040d9de>`_


---------------------
APBS 3.4.0 (Jan 2022)
---------------------

* Binary releases may be found in `GitHub releases <https://github.com/Electrostatics/apbs/releases>`_.

^^^^^^^^^^^^
New Features
^^^^^^^^^^^^

* Revamped build system
* Most submodule switched to using CMake's FetchContent
* FETK is now required; currently using v1.9.2
* Automatic release processes implemented
* Cross-platform builds performed on GitHub Actions
* Pre-compiled binaries posted to each Release
* Binaries are currently single-threaded (no OpenMP)

^^^^^^^^^^^^^^^^^^^^^^^^^^
Known Bugs and Limitations
^^^^^^^^^^^^^^^^^^^^^^^^^^

* Automated build is only single threaded
* pb_solvers has been disabled; requires further development for full integration
* Support for the PYGBE solver is minimal so far; user needs to install PYGBE (e.g. via `pip install pygbe`)
* Docker build was failing during tests and needs to be fixed

^^^^^
Notes
^^^^^

* The following are included in APBS:

  * `Geometric Flow <https://github.com/Electrostatics/geoflow_c/tree/39d53269c084f1dc1caa71de95dca77f19da739e>`_
  * `FETk <https://github.com/Electrostatics/FETK/tree/8c2b67fe587336ba73f77573f13e31ecb1a5a7f9>`_
  * `PBAM/PBSAM <https://github.com/Electrostatics/pb_solvers/tree/d3ba994d7ec2b2cad5b3e843784c7cb9f41ace37>`_
  * `TABI-PB <https://github.com/Treecodes/TABI-PB/tree/fe1c237b057418fed48535db125394607040d9de>`_


-------------------
APBS 3.0 (May 2020)
-------------------

* Binary releases may be found in `GitHub releases <https://github.com/Electrostatics/apbs/releases>`_ and on `SourceForge <http://sourceforge.net/projects/apbs/files/apbs>`_.

^^^^^^^^^^^^
New Features
^^^^^^^^^^^^

* Poisson-Boltzmann Analytical Method (PBAM, see `Lotan & Head-Gordon <http://pubs.acs.org/doi/full/10.1021/ct050263p>`_) and Semi-Analytical Method (PBSAM, see `Yap & Head-Gordon <http://pubs.acs.org/doi/abs/10.1021/ct100145f>`_) integrated with APBS. PBSAM is currently only available in the Linux and OS X distributions.

  * Examples are located with the APBS examples in the pbam/ and pbsam/ directories.

* Tree-Code Accelerated Boundary Integral Poisson-Boltzmann Method (TABI-PB) integrated with APBS (See `Geng & Krasny <http://www.sciencedirect.com/science/article/pii/S0021999113002404>`_).

  * Examples are located with the APBS examples in the bem/, bem-pKa/, and bem-binding-energies/ folders
  * Included NanoShaper alternative to MSMS.
  * More information and documentation may be found in the `Contributions <http://www.poissonboltzmann.org/external_contributions/extern-tabi/>`_ section of the APBS website

* Added binary DX format support to the appropriate APBS tools.
* Test suite amended and expanded.
* Removed hard-coded limitation to number of grid points used to determine surface accessibility.

^^^^^^^^^^^^^^^^^^^^^^^^^^
Known Bugs and Limitations
^^^^^^^^^^^^^^^^^^^^^^^^^^

* PBSAM not building in windows due to C standard restrictions in the Microsoft compiler implementation.

^^^^^^^^^^^^^
Minor Updates
^^^^^^^^^^^^^

* PB(S)AM now requires the key work 'pos' for the term argument.
* PB(S)AM 'surf' keyword has been replaced with the 'usemesh' keyword.
* PB(S)AM 'salt' keyword has been replaced with the 'ion' keyword.
* PB(S)AM dynamics parameters are no longer accepted in the ELEC section.
* PB(S)AM now has only one type of ELEC method: pb(s)am_auto.
* PB(S)AM 'gridpts' keyword has been replaced with 'dime' keyword.
* PB(S)AM 'dx' and '3dmap' keywords are deprecated to use the 'write' one instead.
* BEM mesh keyword now requires method names instead of just integer values.
* GEOFLOW ELEC type has been change from 'geoflow-auto' to 'geoflow'.
* Fixed miscellaneous Windows build issues.
* Update the build configurations for the Pythons libraries.

^^^^^
Notes
^^^^^

* The following are included in APBS as Git submodules:

  * `Geometric Flow <https://github.com/Electrostatics/geoflow_c/tree/e8ce510a670e0b7f3501e72be6141fc20328f947>`_
  * `FETk <https://github.com/Electrostatics/FETK/tree/0c6fdeabe8929acea7481cb1480b5706b343b7e0>`_
  * `PBAM/PBSAM <https://github.com/davas301/pb_solvers/tree/4805cbec02b30e9bae927f03ac2fecd3217c4dad>`_
  * `TABI-PB <https://github.com/lwwilson1/TABIPB/tree/941eff91acd4153a06764e34d29b633c6e3b980f>`_


-------------------
APBS 1.5 (Oct 2016)
-------------------

Dearest APBS users,

We are pleased to announce the latest release of APBS 1.5. The latest version of APBS includes several notable features and bug fixes. This release includes the addition of Poisson-Boltzmann Analytical-Method (PB-AM), Poisson-Boltzmann Semi-Analytical Method (PB-SAM) and the Treecode-Accelerated Boundary Integral Poisson-Boltzmann method (TABI). Additionally, we have made improvements to the build system and the system tests, as well as miscellaneous bug fixes. A full change log may be found `here <https://github.com/Electrostatics/apbs/blob/master/apbs/doc/ChangeLog.md>`_.

For help with installation, building, or running APBS, please visit https://gitter.im/Electrostatics/help.

We thank you for your continued support of APBS.

Sincerely,

The APBS Development Team

-----------------------
APBS 1.4.2.1 (Jan 2016)
-----------------------

^^^^^^^^^^^^
New features
^^^^^^^^^^^^

* Poisson-Boltzmann Semi-Analytical Method (PB-SAM) packaged and built with APBS.
* New Geometric flow API and improvements in speed.
* Support for BinaryDX file format.
* SOR solver added for mg-auto input file option.
* DXMath improvements.
* Test suit improvements:

  * APBS build in Travis-CI
  * Geometric Flow tests added.
  * Protein RNA tests enabled.
  * Intermediate results testing.

* Example READMEs onverted to markdown and updated with latest results. 

^^^^^^^^^
Bug fixes
^^^^^^^^^

* OpenMPI (mg-para) functionality restored.
* Fixed parsing PQR files that contained records other than ATOM and HETATM.
* Geometric Flow boundary indexing bug fixed.
* Build fixes:

  * Out of source CMake build are again working.
  * Python library may be built.
  * CentOS 5 binary builds for glibc compatibility.
  * Pull requests merged.

* Removed irrelevant warning messages.

^^^^^
Notes
^^^^^

The following packages are treated as submodules in APBS:

* Geometric Flow has been moved to its own `repository <https://github.com/Electrostatics/geoflow_c>`_.
* FETk has been `cloned <https://github.com/Electrostatics/FETK>`_ so that we could effect updates.
* PB-SAM lives here:  https://github.com/Electrostatics/PB-SAM

Added a `chat feature <https://gitter.im/Electrostatics/help>`_ for users.

^^^^^^^^^^
Known bugs
^^^^^^^^^^

* Travis CI Linux builds are breaking because Geometric Flow relies on C++11 and Travis boxen have an old GCC that does not support C++11. This also and issue for CentOS 5.
* BEM is temporarily disabled due to build issues.
* Geometric Flow build is currently broken on Windows using Visual Studio.

-----------------------
APBS 1.4.2.0 (Jan 2016)
-----------------------

^^^^^^^^^^^^^
Binary builds
^^^^^^^^^^^^^

Binary releases are available.

^^^^^^^^^^^^
New features
^^^^^^^^^^^^

* Poisson-Boltzmann Semi-Analytical Method (PB-SAM) packaged and build with APBS.
* New Geometric flow API and improvements: https://github.com/Electrostatics/apbs/issues/235
* Support for BinaryDX file format: https://github.com/Electrostatics/apbs/issues/216
* SOR solver added for mg-auto input file option.
* DXMath improvements https://github.com/Electrostatics/apbs/issues/168 and https://github.com/Electrostatics/apbs/issues/216
* Test suite improvements:

  * APBS build in Travis-CI
  * Geometric Flow test added.
  * Protein RNA test enabled https://github.com/Electrostatics/apbs/issues/149
  * Intermediate result testing https://github.com/Electrostatics/apbs/issues/64

* Example READMEs converted to markdown and updated with latest results.

^^^^^^^^^
Bug fixes
^^^^^^^^^

* OpenMPI (mg-para) functionality restored: https://github.com/Electrostatics/apbs/issues/190
* Fized parsing PQR files that contained records other than ATOM and HETATM: https://github.com/Electrostatics/apbs/issues/77 and https://github.com/Electrostatics/apbs/issues/214
* Geometrix Flow boundary indexing bug fixed.
* Build fixes:

  * Out of source CMake build are again working.
  * Python library may be built:  https://github.com/Electrostatics/apbs/issues/372
  * CentOS 5 binary builds for glibc compability.
  * Pull requests merged.

*  Removed irrelevant warning messages: https://github.com/Electrostatics/apbs/issues/378

^^^^^
Notes
^^^^^

* The following packages are treated as submodules in APBS:

  * Geometric Flow has been moved to its own repository:  https://github.com/Electrostatics/geoflow_c/
  * FETk has been cloned: https://github.com/Electrostatics/FETK/
  * PB-SAM lives here:  https://github.com/Electrostatics/PB-SAM/

* Added chat feature at https://gitter.im/Electrostatics/help/ for users. 

^^^^^^^^^^
Known bugs
^^^^^^^^^^

* Travis-CI Linux builds are breaking because Geometric Flow relies on C++11 and Travis boxen have an old GCC that does not support C++11. This is also an issue for CentOS 5.
* BEM is temporarily disabled due to build issues.
* Geometric Flow build is currently broken on Windows using Visual Studio.

---------------------
APBS 1.4.1 (Aug 2014)
---------------------

^^^^^^^
Summary
^^^^^^^

We are pleased to announced the release of APBS 1.4.1. This was primarily a bug fix release; however, we have added a few features we'd like to hightlight below.
We would like to also highlight our new website, still located at: http://www.poissonboltzmann.org. This site is also hosted at GitHub and we hope that the new organization will make it easier for people to find the content they need. While we are still in the process of migrating some remaining content, we have added links to the previous page when needed.
Thank you for your continuing support of APBS. As always, please use our mailing list to send up questions or comments about our software.

^^^^^^^^^^^^^^^^
Detailed changes
^^^^^^^^^^^^^^^^

* Multigrid bug fix for volumes with large problem domain.
* We have added a preliminary implementation of geometric flow.
* Finite element method support has been re-enabled.
* Migration of the APBS source tree to `GitHub <http://github.com/Electrostatics/apbs>`_ for better collaboration, issue tracking, and source code management.
* Improved test suite.

---------------------
APBS 1.4.0 (Jul 2012)
---------------------

^^^^^^^
Summary
^^^^^^^

We are pleased to announce the release of APBS 1.4.0. This version of APBS includes a massive rewrite to eliminate FORTRAN from the software code base to improve portability and facilitate planned optimization and parallelization activities. A more detailed list of changes is provided below.
Starting with this release, we have created separate installation packages for the APBS binaries, examples, and programming documentation. This change is in response to user requests and recognition of the large size of the examples and documentation directories.

^^^^^^^^^^^^^^^^
Detailed changes
^^^^^^^^^^^^^^^^


* Removed FORTRAN dependency from APBS
* Direct line by line translation of all source from contrib/pmgZ
* Functions replaced and tested incrementally to ensure code congruence
* Created new subfolder src/pmgC for translated pmg library
* Created new macros for 2d, 3d matrix access
* In src/generic/apbs/vmatrix.h
* Simulate native FORTRAN 2 and 3 dimensional arrays
* Use 1-indexed, column-major ordering
* Allowed direct 1-1 translation from FORTRAN to ensurre code congruence
* Added additional debugging and output macros to src/generic/apbs/vhal.h
* Added message, error message, assertion, warning, and abort macros
* Macro behavior modified by the --enable-debug flag for configure
* Non-error messages directed to stderr in debug, io.mc otherwise
* All error messages are directed to stdout
* In debug mode, verbose location information is provided
* Added additional flags to configure
* --with-fetk replaces FETK_INCLUDE, FETK_LIBRARY environment flags
* --with-efence enables compiling with electric fence library
* --enable-debug eliminates compiling optimization and includes line no info
* ---enable-profiling adds profiling information and sets --enable-debug
* --enable-verbose-debug prints lots of function specific information

-------------------
APBS 1.3 (Oct 2010)
-------------------

^^^^^^^^^^^^
New features
^^^^^^^^^^^^

* Added in new read and write binary (gz) commands. Can read gzipped DX files directly.
* Added new write format to output the atomic potentials to a flat file (see atompot)
* Added new functionality for using a previously calculated potential map for a new calculation.
* Added a new program for converting delphi potential maps to OpenDX format. tools/mesh/del2dx
* Updated Doxygen manual with call/caller graphs.  Replaced HTML with PDF.
* Added tools/matlab/solver with simple Matlab LPBE solver for prototyping, teaching, etc.
* Deprecated APBS XML output format.
* Deprecated nlev keyword.
* Added etol keyword, which allows user-defined error tolerance in LPBE and NPBE calculations (default errtol value is 1.0e-6).
* Added more explanatory error messages for the case in which parm keyword is missing from APBS input file for apolar calculations.
* Added a polar and apolor forces calculation example to examples/born/ .
* Added warning messages for users who try to compile APBS with --enable-tinker flag and run APBS stand-alone execution.
* Switched default Opal service urls from sccne.wustl.edu to NBCR.
* Added a sanity check in routines.c: 'bcfl map' in the input file requires 'usemap pot' statement in the input file as well.
* Introduced Vpmgp_size() routine to replace F77MGSZ call in vpmg.c
* Updated test results for APBS-1.3 release.
    
^^^^^^^^^
Bug fixes
^^^^^^^^^

* Modified Vpmg_dbForce with some grid checking code provided by Matteo Rotter.
* Fixed a bug in psize.py per Michael Lerner's suggestion. The old version of psize.py gives wrong cglen and fglen results in special cases (e.g., all y coordinates are negative values).
* Fixed a bug in examples/scripts/checkforces.sh: the condition for "Passed with rounding error" is abs(difference) < errortol, not the other way around.
* Fixed the help string in ApbsClient.py .
* Fixed a bug in Vacc_atomdSASA(): the atom SASA needs to be reset to zero displacement after finite melement methods.
* Fixed a bug in Vpmg_dbForce(): the initialization of rtot should appear before it is used.
* Fixed a bug in initAPOL(): center should be initialized before used.
* Fixed a bug in routines.c: eliminated spurious "Invalid data type for writing!" and "Invalid format for writing!" from outputs with "write atompot" statement in the input file.
* Fixed a bug in vpmg.c: fixed zero potential value problem on eges and corners in non-focusing calculations.

---------------------
APBS 1.2.1 (Dec 2009)
---------------------

^^^^^^^^^
Bug fixes
^^^^^^^^^

* Added in warning into focusFillBound if there is a large value detected in setting the boundary conditions during a focusing calculation
* Added in a check and abort in Vpmg_qmEnergy if chopped values are detected. This occurs under certain conditions for NPBE calculations where focusing cuts into a low-dielectric regions.
* Fixed a bug in Vpmg_MolIon that causes npbe based calculations to return very large energies.  This occurs under certain conditions for NPBE calculations where focusing cuts into a low-dielectric regions.

---------------------
APBS 1.2.0 (Oct 2009)
---------------------

^^^^^^^^^^^^
New features
^^^^^^^^^^^^

* Updated NBCR opal service urls from http://ws.nbcr.net/opal/... to http://ws.nbcr.net/opal2/... 
* Increased the number of allowed write statements from 10 to 20
* Updated inputgen.py with --potdx and --istrng options added, original modification code provided by Miguel Ortiz-Lombardía
* Added more information on PQR file parsing failures
* Added in support for OpenMP calculations for multiprocessor machines.
* Changed default Opal service from http://ws.nbcr.net/opal2/services/APBS_1.1.0 to http://sccne.wustl.edu:8082/opal2/services/apbs-1.2

^^^^^^^^^^^^^
Modifications
^^^^^^^^^^^^^

* Applied Robert Konecny's patch to bin/routines.h (no need to include apbscfg.h in routines.h)

^^^^^^^^^
Bug fixes
^^^^^^^^^

* Added a remove_Valist function in Python wrapper files, to fix a memory leak problem in pdb2pka
* Fixed a bug in smooth.c: bandwidth iband, jband and kband (in grid units) should be positive integers
* Fixed a bug in psize.py: for a pqr file with no ATOM entries but only HETATM entries in it, inputgen.py should still create an APBS input file with reasonable grid lengths
* Fixed a bug in Vgrid_integrate: weight w should return to 1.0 after every i, j or k loop is finished
* Fixed a bug in routines.c, now runGB.py and main.py in tools/python/ should be working again instead of producing segfault
* Fixed a few bugs in ApbsClient.py.in related to custom-defined APBS Opal service urls, now it should be OK to use custom-defined APBS Opal service urls for PDB2PQR web server installations

---------------------
APBS 1.1.0 (Mar 2009)
---------------------

^^^^^^^^^^^^
New features
^^^^^^^^^^^^

* Moved APBS user guide and tutorial to MediaWiki
* Added in support for OpenMPI for parallel calculations
* Added in command line support for Opal job submissions (Code by Samir Unni)
* Allowed pathname containing spaces in input file, as long as the whole pathname is in quotes ("")
* Documented 'make test' and related features

^^^^^^^^^^^^^
Modifications
^^^^^^^^^^^^^

* Modified the function bcCalc to march through the data array linearly when setting boundary conditions. This removes duplication of grid points on the edge of the array and corners.
* Clarified documentation on the IDs assigned to input maps, PQRs, parameter files, etc.
* pdated tutorial to warn against spaces in APBS working directory path in VMD; updated user guide to warn against spaces in APBS installation path on Windows
* 'make test' has been reconfigured to run before issuing make install (can be run from top directory)
* Removed tools/visualization/vmd from tools directory in lieu of built-in support in VMD
* Path lengths can now be larger than 80 characters
* Expanded authorship list
* Added in 'make test-opal' as a post install test (run from the examples install directory)
* Added additional concentrations to protein-rna test case to better encompass experimental conditions used by Garcia-Garcia and Draper; this improves agreement with the published data

^^^^^^^^^
Bug fixes
^^^^^^^^^

* Fixed typos in User Guide (ion keyword) and clarified SMPBE keyword usage
* Fixed typo in User Guide (writemat: poission -> poisson)
* Updated psize.py with Robert's patch to fix inconsistent assignment of fine grid numbers in some (very) rare cases
* Fixed bug with boundary condition assignment.  This could potentially affect all calculations; however, probably has limited impact:  many test cases gave identical results after the bug fix; the largest change in value was < 0.07%.

---------------------
APBS 1.0.0 (Apr 2008)
---------------------

^^^^^^^^^^^^
New features
^^^^^^^^^^^^


* Changed license to New BSD style open source license (see http://www.opensource.org/licenses/bsd-license.php) for more information
* Added in a feature limited version of PMG (Aqua) that reduces the memory footprint of an APBS run by 2-fold
* Modified several routines to boost the speed of APBS calculations by approximately 10% when combined with the low memory version of APBS
* Simplified parameter input for ION and SMPBE keywords (key-value pairs) 
* Examples and documentation for size-modified PBE code (Vincent Chu et al)
* Added in "fast" compile-time option that uses optimized parameters for multigrid calculations
* mg-dummy calculations can be run with any number (n>3) of grid points
* Updated PMG license to LGPL
* Added per-atom SASA information output from APOLAR calculations
* Added per-atom verbosity to APOLAR calculation outputs
* Ability to read-in MCSF-format finite element meshes (e.g., as produced by Holst group GAMER software)
* Updated installation instructions in user guide
* Updated inputgen.py to write out the electrostatic potential for APBS input file.

^^^^^^^^^
Bug fixes
^^^^^^^^^

* Updated tools/python/apbslib* for new NOsh functionality
* Clarified ELEC/DIME and ELEC/PDIME documentation
* Added more transparent warnings/error messages about path lengths which exceed the 80-character limit
* Fixed small typo in user guide in installation instructions
* Fixed memory leaks throughout the APBS code
* Fixed NOsh_parseREAD errors for input files with \r line feeds.
* Fixed a variable setting error in the test examples
* Fixed a bug where memory usage is reported incorrectly for large allocations on 64-bit systems
* Added DTRSV to APBS-supplied BLAS to satisfy FEtk SuperLU dependency
* Fixed a small bug in routines.c to print out uncharged molecule id
* Limited calculation of forces when surface maps are read in 

---------------------
APBS 0.5.1 (Jul 2007)
---------------------

^^^^^^^^^^^^
New features
^^^^^^^^^^^^

* Replaced APOLAR->glen and APOLAR->dime keywords with APOLAR->grid
* Deprecated mergedx. Added mergedx2
    
    * mergedx2 takes the bounding box that a user wishes to calculate a map for, as well as a resolution of the output map. An output map meeting those specifications is calculated and store.
    
* Added pKa tutorial
* Added warning about strange grid settings (MGparm)
* Fixed a bug in vpmg.c that occured when a user supplied a dielectric map with the ionic strength set to zero, causing the map to not be used.
* Removed deprecated (as of 0.5.0) tools/manip/acc in lieu of new APOLAR APBS features
* Added enumerations for return codes, new PBE solver (SMPBE) and linear/ nonlinear types
* Added in code for Size-Modified PBE (SMPBE)


^^^^^^^^^^^^^^^^^^^^^^^^^
Bug fixes and API changes
^^^^^^^^^^^^^^^^^^^^^^^^^

* Fixed buffer over-run problem
* Fixed case inconsistency with inputgen.py and psize.py scripts which caused problems with some versions of Python
* Fixed bug wherein 'bcfl sdh' behaved essentially like 'bcfl zero'.  Now we have the correct behavior:  'bcfl sdh' behaves like 'bcfl mdh' thanks to the multipole code added by Mike Schnieders.  Interestingly, this bug didn't have a major on the large-molecule test cases/examples provided by APBS but did affect the small molecule solvation energies.  Thanks to Bradley Scott Perrin for reporting this bug.
* Added support for chain IDs in noinput.py
* Fixed bug in noinput.py where REMARK lines would cause the script to fail.

---------------------
APBS 0.5.0 (Jan 2007)
---------------------

^^^^^^^^^^^^
New features
^^^^^^^^^^^^

* Significantly streamlined the configure/build/install procedure:
    
    * Most common compiler/library options now detected by default
    * MALOC is now included as a "plugin" to simplify installation and compatibility issue
    
* Added new APOLAR section to input file and updated documentation -- this function renders tools/manip/acc obsolete.
* Added support for standard one-character chain IDs in PQR files. 
* Added a new "spl4" charge method (chgm) option to support a quintic B-spline discretization (thanks to Michael Schnieders).
* Updated psize.py
* Added a new "spl4" ion-accessibility coefficient model (srfm) option that uses a 7th order polynomial. This option provides the higher order continuity necessary for stable force calculations with atomic multipole force fields (thanks to Michael Schnieders).
* Modified the "sdh" boundary condition (bcfl) option to include dipoles and quadrupoles.  Well-converged APBS calculations won't change with the dipole and quadrupole molecular moments included in the boundary potential estimate, but calculations run with the boundary close to the solute should give an improved result (thanks to Michael Schnieders). 
* Updated documentation to reflect new iAPBS features (NAMD support)
* Added Gemstone example to the tutorial
* New example demonstrating salt dependence of protein-RNA interactions.
* Added code to allow for an interface with TINKER (thanks to Michael Schnieders).
* The Python wrappers are now disabled by default.  They can be compiled by passing the --enable-python flag to the configure script.  This will allow users to attempt to compile the wrappers on various systems as desired.
* Added XML support for reading in parameter files when using PDB files as input.  New XML files can be found in tools/conversion/param/vparam.
* Added XML support for reading "PQR" files in XML format.
* Added support for command line --version and --help flags. 
* Added support for XML output options via the --output-file and  --output-format flags.
* Updated runme script in ion-pmf example to use environmental variable for APBS path
* Modified the license to allow exceptions for packaging APBS binaries with several visualization programs.  PMG license modifed as well.
* Added a DONEUMANN macro to vfetk.c to test FEM problems with all Neumann boundary conditions (e.g., membranes).
* Added Vpmg_splineSelect to select the correct Normalization method with either cubic or quintic (7th order polynomial) spline methods.
* Modified the selection criteria in Vpmg_qfForce, Vpmg_ibForce and Vpmg_dbnpForce for use with the new spline based (spl4) method. 
* Added ion-pmf to the make test suite.
* Updated splash screen to include new PMG acknowledgment
* Added runGB.py and readGB.py to the Python utilities, which calculate solvation energy based on APBS parameterized Generalized Born model.
* Updated authorship and tool documentation
* Deprecated ELEC->gamma keyword in lieu of APOLAR->gamma

^^^^^^^^^^^^^^^^^^^^^^^^^
Bug fixes and API changes
^^^^^^^^^^^^^^^^^^^^^^^^^

* Cleanup of documentation, new Gemstone example
* Clarified usage of dime in mg-para ELEC statements
* Massive cleanup of NOsh, standardizing molecule and calculation IDs and making the serial focusing procedure more robust
* Removed MGparm partOlap* data members; the parallel focusing centering is now done entirely within NOsh
* Updated the user manual and tutorial
* Updated psize.py to use centers and radii when calculating grid sizes (thanks to John Mongan)
* Fixed problems with FEM-based NPBE, LPBE, and LRPBE calculations
* Fixed a minor bug in the configure script that prevented MPI libraries from being found when using certain compilers.
* Updated acinclude.m4, aclocal.m4, config/* for new version (1.9) of automake and compatibility with new MALOC
* Fixed a bug where reading in a file in PDB format had atom IDs starting  at 1 rather than 0, leading to a segmentation fault.
* Fixed a bug in mypde.f where double precision values were initialized with single precision number (causing multiplication errors).
* Fixed a bug in the FEM code. Now calls the npbe solver works properly with FEtk 1.40
* Modified the FEMParm struct to contain a new variable pkey, which is  required for selecting the correct path in AM_Refine

---------------------
APBS 0.4.0 (Dec 2005)
---------------------

^^^^^^^^^^^^
New features
^^^^^^^^^^^^

* New version of the 'acc' program available.
* Added additional verbosity to APBS output.
* Added tools/python/vgrid to the autoconf script. The directory compiles with the rest of the Python utilities and is used for manipulating dx files.
* Modified the tools/python/noinput.py example to show the ability to get and print energy and force vectors directly into Python arrays.
* Added dx2uhbd tool to tools/mesh for converting from dx format to UHBD format (Thanks to Robert Konecny)
* Added ability of tools/manip/inputgen.py to split a single mg-para APBS input file into multiple asynchronous input files.
* Modified inputgen.py to be more flexible for developers wishing to directly interface with APBS.
* Added Vclist cell list class to replace internal hash table in Vacc
* Modified Vacc class to use Vclist, including changes to the Vacc interface (and required changes throughout the code)
* Consolidated Vpmg_ctor and Vpmg_ctorFocus into Vpmg_ctor
* Consolidated vpmg.c, vpmg-force.c, vpmg-energy.c, vpmg-setup.c
* Added autoconf support for compilation on the MinGW32 Windows Environment
* Added autoconf support (with Python) for Mac OS 10.4 (Tiger)
* Added the function Vpmg_solveLaplace to solve homogeneous versions of Poisson's equation using Laplacian eigenfunctions.
* Modified the dielectric smoothing algorithm (srfm smol) to a 9 point method based on Bruccoleri, et al.  J Comput Chem 18 268-276 (1997).  NOTE:  This is a faster and more flexible smoothing method.  However, when combined with the the molecular surface bugfixes listed below, this change has the potential to make the srfm smol method give very different results from what was calculated in APBS 0.3.2.  Users who need backwards compatibility are encouraged to use the spline based smoothing method (srfm spl2) or the molecular surface without smoothing (srfm mol).
* Added new 'sdens' input keyword to allow user to control the sphere density used in Vacc.  This became necessary due to the Vacc_molAcc bug fix listed below.  Only applies to srfm mol and srfm smol.
* Made the examples directory documentation much more streamlined.
* Added tests for examples directory.  Users can now issue a "make test" in the desired directory to compare local results with expected results. Also includes timing results for tests for comparison between installations.

^^^^^^^^^
Bug fixes
^^^^^^^^^

* Fixed a bug in Vpmg_qmEnergy to remove a spurious coefficient of z_i^2 from the energy calculation.  This generated incorrect results for multivalent ions (but then again, the validity of the NPBE is questionable for multivalents...)  (Big thanks to Vincent Chu)
* Fixed a bug in vacc.c where atoms with radii less than 1A were not considered instead of atoms with no radii.
* Fixed error in tools/mesh/dx2mol.c (Thanks to Fred Damberger)
* Fixed floating point error which resulted in improper grid spacings for some cases.
* Fixed a bug in Vacc_molAcc which generates spurious regions of high internal dielectric for molecular surface-based dielectric definitions.  These regions were very small and apparently affected energies by 1-2% (when used with the 'srfm mol'; the 'srfm smol' can potentially give larger deviations).  The new version of the molecular surface is actually faster (requires 50-70% of the time for most cases) but we should all be using the spline surface anyway -- right? (Thanks to John Mongan and Jessica Swanson for finding this bug).
* Fixed a bug in vpmg.c that caused an assertion error when writing out laplacian maps (Thanks to Vincent Chu).
* Ensured Vpmg::ccf was always re-initialized (in the case where the Vpmg object is being re-used).
* Removed a spurious error estimation in finite element calculations.
* Clarified the role of ccf and other variables in mypde.f and vpmg.c by expanding/revising the inline comments.

---------------------
APBS 0.3.2 (Nov 2004)
---------------------

^^^^^^^^^^^^
New features
^^^^^^^^^^^^

* Updated tutorial with more mg-auto examples
* Updated apbs.spec file for generating RPMs on more platforms.
* Added new Python wrapper to tools/python directory showing how to run APBS without PQR and .in inputs.
* Python wrappers are now configured to compile on more architectures/ from more compilers.
* Updated tools/conversion/pdb2pqr to a new version (0.1.0) of PDB2PQR, which now can handle nucleic acids, rebuild missing heavy atoms, add hydrogens, and perform some optimization.

^^^^^^^^^
Bug fixes
^^^^^^^^^

* The dimensions of the fine grids in the pka-lig example calculations were increased to give more reliable results (albeit ones which don't agree with the reported UHBD values as well).
* hz in mgparse.c causes name clash with AIX environmental variable; fixed.
* Fixed documentation to state that using a kappa map does not ignore ELEC ION statements.
* Added a stability fix for printing charge densities for LPBE-type calculations.
* Fixed a bug in NPBE calculations which led to incorrect charge densities and slightly modified total energies.
* Modified the origin when creating UHBD grids to match standard UHBD format.
* Fixed VASSERT error caused by rounding error when reading in dx grid files.

---------------------
APBS 0.3.1 (Apr 2004)
---------------------

^^^^^^^^^^^^
New features
^^^^^^^^^^^^

* New APBS tutorial
* New :file:`tools/python/vgrid/mergedx.py` script to merge dx files generated from parallel APBS runs back into a single dx file.

^^^^^^^^^
Bug fixes
^^^^^^^^^

* Fixed bug in parallel calculations where atoms or points on a border between two processors were not included.  Modified setup algorithm for parallel calculations to allow partitions in order to obtain grid points and spacing from the global grid information.
* Modified extEnergy function to work with parallel calculations, giving better accuracy.

---------------------
APBS 0.3.0 (Feb 2004)
---------------------

^^^^
News
^^^^

APBS is now supported by the NIH via NIGMS grant GM69702-01.

^^^^^^^^^^^^^^^^^^^^^^^^^
Changes that affect users
^^^^^^^^^^^^^^^^^^^^^^^^^

* New version of the documentation
* New directory structure in tools/
* Finished fe-manual mode for ELEC calculations -- this is the adaptive finite element solver
* Added documetnation for fe-manual
* New apbs/tools/manip/inputgen.py script to automatically generate input APBS files from PQR data
* Added new asynchronous mode in mg-para parallel calculations to enable running on-demand and/or limited resources
* Added new script (tools/manip/async.sh) to convert mg-para calculations in mg-async calculations
* Added following aliases for some of the more obscure parameters in the input files:

  * chgm 0 ==> chgm spl0
  * chgm 1 ==> chgm spl2
  * srfm 0 ==> srfm mol
  * srfm 1 ==> srfm smol
  * srfm 2 ==> srfm spl2
  * bcfl 0 ==> bcfl zero
  * bcfl 1 ==> bcfl sdh
  * bcfl 2 ==> bcfl mdh
  * bcfl 4 ==> bcfl focus
  * calcenergy 0 ==> calcenergy no
  * calcenergy 1 ==> calcenergy total
  * calcenergy 2 ==> calcenergy comps
  * calcforce 0 ==> calcforce no
  * calcforce 1 ==> calcforce total
  * calcforce 2 ==> calcforce comps

* Example input files have been updated to reflect this change. NOTE: the code is backward-compliant; i.e., old input files WILL still work.
* Added new READ options "PARM" and "MOL PDB", see documentation for more information. These options allow users to use unparameterized PDB files together with a parameter database.
* Updated the documentation
* Now include support for chain IDs and other optional fields in PQR/PDB files
* Added support for parsing PDB files
* Renamed:

* amber2charmm -> amber2charmm.sh
* pdb2pqr -> pdb2pqr.awk
* qcd2pqr -> qcd2pqr.awk

* Added a new Python-based pdb2pqr (tools/conversion/pdb2pqr) script that allows users to choose parameters from different forcefields.
* Updated Python wrappers (tools/python) and added the python directory to autoconf/automake.
* Reformatted examples/README.html for readability.

^^^^^^^^^
Bug fixes
^^^^^^^^^

* Fixed bug in PQR parsing that can cause PDB/PQR files to be mis-read when they contain residues with numbers in their names (Thanks to Robert Konecny and Joanna Trylska)
* Fixed bug when writing out number/charge density: unrealistic densities reported inside iVdW surface.
* Fixed bug in VMD read_dx utility
* Invalid map IDs now result in an error message instead of a core dump (thanks to Marco Berrera)
* Modified mechanism for cubic-spline output, fixing a bug associated with zero-radius atoms
* Fixed omission of srfm in sections of documentation (thanks to Sameer Varma)
* Made autoconf/automake configure setup more robust on Solaris 8 platforms (thanks to Ben Carrington)
   
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Changes that affect developers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* New docuemtnation setup
* New tools/ directory structure
* Changed Vgreen interface and improved efficiency
* Changed Vopot interface to support multiple grids
* Added several norm and seminorm functions to Vgrid class
* Altered --with-blas syntax in configure scripts and removed --disable-blas
* Documented high-level frontend routines
* Cool new class and header-file dependency graphs courtesy of Doxygen and Graphviz
* Added substantial mypde.c-based functionality to Vfetk
* Moved chgm from PBEparm to MGparm
* Minor changes to Vfetk: removed genIcos and added genCube
* FEM solution of RPBE working again (see test/reg-fem) and is probably more up-to-date than test/fem
* Updated API documentation
* Changed many NOsh, FEMparm, MGparm variables to enums
* Changes to Valist and Vatom classes
* Fixed minor bugs in documentation formatting
* Made Vopot more robust
* Created Vparam class for parameter file parsing
* Added vparam* parameter database flat files to tools/conversion/param

---------------------
APBS 0.2.6 (Jan 2003)
---------------------

* Changed license to GPL
* Made a few routines compliant with ANSI X3.159-1989 by removing snprintf (compliant with ISO/IEC 9899:1999).  This is basically for the sake of OSF support.

---------------------
APBS 0.2.5 (Nov 2002)
---------------------

* Improved consistency between energies evaluated with "chgm 0" and "chgm 1"
* Made charge-field energy evaluation consistent for user-supplied charge maps
* Added new psize.py script courtesy of Todd Dolinsky.
* Updated list of APBS-related tools in User Guide.
* Added RPM capabilities courtesy of Steve Bond.
* Removed annoying excess verbosity from Vgrid.
* Updated Blue Horizon compilation instructions (thanks to Robert Konecny and Giri Chukkapalli)
* Updated autoconf/automake/libtool setup and added --disable-tools option

---------------------
APBS 0.2.4 (Oct 2002)
---------------------

* Fixed bug which set one of the  z-boundaries to zero for "bcfl 1".  This can perturb results for systems where the grid boundaries are particularly close to the biomolecule.  While this is an embarassing bug, most systems using settings suggested by the psize script appear largely unaffected (see examples/README.html).  Thanks to Michael Grabe for finding this bug (Michael, you can stop finding bugs now...)
* Updated VMD scripts to agree with the current OpenDX output format
* A COMMENT:  As far as I can tell, the current version of OpenDX-formatted output (same as version 0.2.3) is fully compliant with the OpenDX standards (see, for example,  http://opendx.npaci.edu/docs/html/pages/usrgu065.htm#HDRXMPLES).   However, I realize this is different than the format for previous versions and would encourage all users to update their APBS-based applications to accomodate these changes.  The best solution would be for all downstream applications to use the APBS Vgrid class (see http://agave.wustl.edu/apbs/doc/api/html/group__Vgrid.html) to manipulate the data since this class should remain internally consistent between releases.  Finally, I would love to have some OpenDX guru who uses APBS to contact me so I can solidfy the data ouput format of APBS.  I'm about ready to permanently switch to another format if I can't reach a consensus with OpenDX...

---------------------
APBS 0.2.3 (Oct 2002)
---------------------

* Fixed bugs in salt-dependent Helmholtz/nonlinear term of PBE affecting both LPBE and NPBE calculations.  While this bug fix only changes most energies by < 2 kJ/mol, it is recommended that all users upgrade.  Many thanks to Michael Grabe for finding and carefully documenting this bug!
* A parameter (chgm) has been added which controls the charge discretization method used.  Therefore, this version contains substantial changes in both the API and input file syntax.  Specifically:
    
    * PBEparm has two new members (chgm, setchgm)
    * Vpmg_fillco requires another argument
    * Vpmg\_*Force functions require additional arguments
    * Input files must now contain the keyword "chgm #" where # is an integer
    * Please see the documentation for more information.
    
* Fixed problems with "slicing" off chunks of the mesh during I/O of focused calculations
* Updated authors list
* New CHARMM parameters -- Robert Konecny
* Created enumerations for common surface and charge discretization methods
* Added Vmgrid class to support easy manipulation of nested grid data
* Added more verbosity to error with NPBE forces
* Added working Python wrappers -- Todd Dolinksy
* Modified VMD scripts read_dx and loadstuff.vmd

---------------------
APBS 0.2.2 (Aug 2002)
---------------------

* There were several other changes along the way... I lost track.
* Changed coordinate indexing in some energy calculations
* Updated documentation to reflect recent changes on Blue Horizon
* Improved speed of problem setup BUT NOW RESTRICT use of input coefficient maps (see documentation)
* Updated documentation, placing particular emphasis on use of Intel compilers and vendor BLAS on Intel Linux systems
* Fixed bug for nonpolar force evaluation in Vpmg_dbnpForce
* Removed MG test scripts; use :file:`bin/*.c` for templates/testing
* Made main driver code completely memory-leak free (i.e., if you wanted to wrap it and call it repeatedly -- Thanks to Robert Konecny)
* Fixed main driver code for compatibility with SGI compilers (Thanks to Fabrice Leclerc)
* Made focused evaluation of forces more sensible.
* Added 'print force' statement
* Fixed bug in OpenDX input/output (OpenDX documentation is lousy!)

---------------------
APBS 0.2.1 (Apr 2002)
---------------------

This version requires the latest version of MALOC to work properly!

* Syntax changes
    
    * The writepot and writeacc keywords have been generalized and new I/O features have been added.  The syntax is now:
        
        * write pot dx potential
        * write smol dx surface
        * etc.  Please see the User's Manual for more information
        
    * The read keywords has been generalized and new I/O features have been added which support the use of pre-calculated coefficient grids, etc.  The correct syntax for reading in a molecule is now "read mol pqr mol.pqr end"; please see the User's Manual for more information.
    * The "mg" keyword is no longer supported; all input files should use "mg-manual" or one of the other alternatives instead.
    
* A change in the behavior of the "calcenergy" keyword; passing an argument of 2 to this keyword now prints out per-atom energies in addition to the energy component information
* A new option has been added to tools/manip/acc to give per-atom solvent-accessible surface area contributions
* A new option has been added to tools/manip/coulomb to give per-atom electrostatic energies
* Added tools/mesh/dxmath for performing arithmetic on grid-based data (i.e., adding potential results from two calculations, etc.)
* Added tools/mesh/uhbd_asc2bin for converting UHBD-format grid files from ASCII to binary (contributed by Dave Sept)
* Improvement of VMD visualization scripts (contributed by Dave Sept)
* The API has changed significantly; please see the Programmer's Manual.
* Working (but still experimental) Python wrappers for major APBS functions.
* More flexible installation capabilities (pointed out by Steve Bond)
* Added ability to use vendor-supplied BLAS
* Brought up-to-date with new MALOC

---------------------
APBS 0.2.0 (Mar 2002)
---------------------

This version is a public (beta) release candidate and includes:

* Slight modification of the user and programmer's guides
* Scripts for visualization of potential results in VMD (Contributed by Dave Sept)
* Corrections to some of the example input files
* A few additional API features

This release requires a new version of MALOC. 

---------------------
APBS 0.1.8 (Jan 2002)
---------------------

This version is a public (beta) release candidate and includes the following bug-fixes:

* Added warning to parallel focusing 
* Added several test cases and validated the current version of the code for all but one (see examples/README.html)
* Fixed atom partitioning bug and external energy evaluation during focusing
* Added new program for converting OpenDX format files to MOLMOL (by Jung-Hsin Lin)

You should definitely upgrade, the previous versions may produce unreliable results.

---------------------
APBS 0.1.7 (Dec 2001)
---------------------

This version is a public (beta) release candidate and includes the following bug-fixes:

* Fixed I/O for potential in UHBD format (thanks, Richard!)
* Re-arranged garbage collection routines in driver code
* Improved FORTRAN/C interfaces
* Re-configured autoconf/libtool setup for more accurate library version number reporting

---------------------
APBS 0.1.6 (Nov 2001)
---------------------

This version is a public (beta) release candidate and includes the following bug-fixes and features:

* Fixed printf formatting in UHBD potential output
* Added input file support for parallel focusing
* Fixed small bug in parsing writeacc syntax (thanks, Dave)
* Added output file support for parallel focusing
* Changed some documentation

You need to download a new version of MALOC for this release.   

---------------------
APBS 0.1.5 (Oct 2001)
---------------------

This version features minor bug fixes and several new features:

* Fixed shift in center of geometry for OpenDX I/O
* Made energy evaluation more robust when using NPBE
* Rearrangments of files and modified compilation behavior
* Input file support for ion species of varying valency and concentration
* Input file support incorrect nlev/dime combinations; APBS now finds acceptable settings near to the user's requested values
* "Automatic focusing".  Users now simply specify the physical parameters (temperature, dielectric, etc.), the coarse and fine grid lengths and centers, and APBS calculates the rest

---------------------
APBS 0.1.4 (Sep 2001)
---------------------

This version features major bug fixes introduced in the 0.1.3 release:

* Chain ID support has been **removed** from the PDB/PQR parser (if anyone has a nice, flexible PDB parser they'd like to contribute to the code, I'd appreciate it)
* Configure script has been made compatible with OSF
* Bug fix in disabling FEtk-specific header files

---------------------
APBS 0.1.3 (Sep 2001)
---------------------

This version features a few improvements in scripts, PDB parsing flexibility, and portability, including:

* Dave Sept upgraded the psize and shift scripts to allow more flexibility in PDB formats.
* Chain ID support has been added to the PDB/PQR parser
* Removed -g from compiler flags during linking of C and FORTAN under OSF (thanks to Dagmar Floeck and Julie Mitchell for help debugging this problem)

---------------------
APBS 0.1.2 (Sep 2001)
---------------------

This version is mainly designed to increase portability by switching to libtool for library creation and linking.
Of course, it also contains a few bug fixes.
Highlights include:

* Changes to the User Manual
* Addition of a Programmer's Manual
* Various FEtk-related things (no particular impact to the user)
* Improvements to the test systems
* Change in the format for printing energies
* Change in directory structure
* Fixed centering bug in main driver (only impacted I/o)
* Fixed error message bug in VPMG class
* Fixed grid length bug (popped up during sanity checks) in VPMG class
* Switched to libtool for linking
* Note that Compaq Tru64 Alpha users may still experience problems while compiling due to some strangess with linking C and FORTRAN objects.

---------------------
APBS 0.1.1 (Aug 2001)
---------------------

I am slightly less pleased to announce the first bug-fix for APBS, version 0.1.1.
This fixes compilation problems that popped up for several folks,
including:

* Syntax errors with non-GNU compilers
* Errors in the installation instructions
* Installation of binary in machine-specific directory

---------------------
APBS 0.1.0 (Aug 2001)
---------------------

I am pleased to announce the availability of a pre-beta version of the Adaptive Poisson-Boltzmann Solver (APBS) code to selected research groups.
APBS is new software designed to solve the Poisson-Boltzmann equation for very large biomolecular systems.
For more information, please visit the APBS web site at http://mccammon.ucsd.edu/apbs.

This release is designed to allow interested users to get familiar with the code. 
It is not currently fully functional; it only provides for the sequential multigrid (Cartesian mesh) solution of the linearized and nonlinear Poisson-Boltzmann equation.
User-friendly parallel support will be incorporated into the next release.
Other limitations that may impact its immediate usefulness are:

* No finite element support.  This is awaiting the public release of the Holst group's FEtk library.
* Somewhat inefficient coefficient evaluation (i.e., problem setup).  This should be fixed in the next release or two.

Rather than serving as a production code, this release represents a request for help in breaking the software and finding its deficiencies
before a public beta.

If you are interested in testing this early release, please go to http://wasabi.ucsd.edu/~nbaker/apbs/download/.
Since this is not a public release of APBS, you will need to enter the user-name "apbs-beta" and the password "q94p$fa!" for access to this site.
Once there, please follow the instructions to download and install APBS.

If you are not interested in trying out this early release, but would like to stay informed about subsequent versions of APBS, please consider subscribing to the APBS announcements mailing list by sending the message "subscribe apbs-announce" to majordomo@mccammon.ucsd.edu.

Thank you for your time and interest in the APBS software.

""""""""""""""""""""""""""
Documentation "to-do" list
""""""""""""""""""""""""""
.. todolist::
.. _apbsxmlparm:

APBS XML parameter format
=========================

This parameter file format has the following form:

.. code-block:: xml

   <ffname>
      <residue>
          <name>resname</name>
          <atom>
              <name>atomname</name>
              <charge>atomcharge</charge>
              <radius>atomradius</radius>
              <epsilon>atomepsilon</epsilon>
          </atom>
          ...
      </residue>
      ...
   </ffname>

The variables in this example are:

``ffname``
  The name of the forcefield. This is the root element of the XML file.

``resname``
  A string giving the residue name, as provided in the PDB file to be parameterized.

``atomname``
  A string giving the atom name, as provided in the PDB file to be parameterized.

``atomcharge``
  A float giving the atomic charge (in electrons).

``atomradius``
  A float giving the atomic Radius (in Å).

``atomepsilon``
  A float giving the Lennard-Jones well depth :math:`\epsilon` (in kJ/mol).
  This is used for the calculation of WCA energies in apolar solvation energies and forces.
  We assume that the Lennard-Jones potential is defined in the "AMBER style".. _pqr:

PQR molecular structure format
==============================

This format is a modification of the PDB format which allows users to add charge and radius parameters to existing PDB data while keeping it in a format amenable to visualization with standard molecular graphics programs.
The origins of the PQR format are somewhat uncertain, but has been used by several computational biology software programs, including MEAD and AutoDock.
UHBD uses a very similar format called QCD.

APBS reads very loosely-formatted PQR files: all fields are whitespace-delimited rather than the strict column formatting mandated by the PDB format.
This more liberal formatting allows coordinates which are larger/smaller than ± 999 Å.
APBS reads data on a per-line basis from PQR files using the following format:::

  Field_name Atom_number Atom_name Residue_name Chain_ID Residue_number X Y Z Charge Radius

where the whitespace is the most important feature of this format.
The fields are:

``Field_name``
  A string which specifies the type of PQR entry and should either be ATOM or HETATM in order to be parsed by APBS.

``Atom_number``
  An integer which provides the atom index.

``Atom_name``
  A string which provides the atom name.

``Residue_name``
  A string which provides the residue name.

``Chain_ID``
  An optional string which provides the chain ID of the atom.
  Note that chain ID support is a new feature of APBS 0.5.0 and later versions.

``Residue_number``
  An integer which provides the residue index.

``X Y Z``
  3 floats which provide the atomic coordinates (in Å)

``Charge``
  A float which provides the atomic charge (in electrons).

``Radius``
  A float which provides the atomic radius (in Å).

Clearly, this format can deviate wildly from PDB due to the use of whitespaces rather than specific column widths and alignments.
This deviation can be particularly significant when large coordinate values are used.
However, in order to maintain compatibility with most molecular graphics programs, the PDB2PQR program and the utilities provided with APBS attempt to preserve the PDB format as much as possible.
.. _pdb:

PDB molecular structure format
==============================

The PDB file format is described in detail in the `Protein Data Bank documentation <http://www.rcsb.org/pdb/static.do?p=file_formats/pdb/index.html>`_.
Harwell-Boeing matrix format
============================

This is the sparse matrix output format used by APBS for analyses of the matrix operators which are constructed during PB solution.
This format was implemented so matrix operators could by decomposed with SuperLU and ARPACK but this also serves as a useful general mechanism for sparse matrix input and output.
.. _mcsf:

MCSF mesh format
================

APBS reads and writes meshes in the `FEtk MCSF format <http://fetk.org/codes/mc/>`.
.. _opendx:

OpenDX scalar data format
=========================

We output most discretized scalar data (e.g., potential, accessibility, etc.) from APBS in the data format used by the OpenDX software package.
The OpenDX data format is very flexible; the following sections describe the application of this format for APBS multigrid and finite element datasets.

The multigrid data format has the following form:

.. code-block:: bash

   object 1 class gridpositions counts nx ny nz
   origin xmin ymin zmin
   delta hx 0.0 0.0
   delta 0.0 hy 0.0 
   delta 0.0 0.0 hz
   object 2 class gridconnections counts nx ny nz
   object 3 class array type double rank 0 items n data follows
   u(0,0,0) u(0,0,1) u(0,0,2)
   ...
   u(0,0,nz-3) u(0,0,nz-2) u(0,0,nz-1)
   u(0,1,0) u(0,1,1) u(0,1,2)
   ...
   u(0,1,nz-3) u(0,1,nz-2) u(0,1,nz-1)
   ...
   u(0,ny-1,nz-3) u(0,ny-1,nz-2) u(0,ny-1,nz-1)
   u(1,0,0) u(1,0,1) u(1,0,2)
   ...
   attribute "dep" string "positions"
   object "regular positions regular connections" class field
   component "positions" value 1
   component "connections" value 2
   component "data" value 3`

The variables in this format include:

``nx ny nz``
  The number of grid points in the x-, y-, and z-directions

``xmin ymin zmin``
  The coordinates of the grid lower corner

``hx hy hz``
  The grid spacings in the x-, y-, and z-directions.

``n``
  The total number of grid points; :math:`n = nx * ny * nz`

``u(*,*,*)``
  The data values, ordered with the z-index increasing most quickly, followed by the y-index, and then the x-index.

For finite element solutions, the OpenDX format takes the following form:

.. code-block:: bash

   object 1 class array type float rank 1 shape 3 items N
   v1x v1y v1z
   v2x v2y v2z
   ...
   vNx vNy vNz
   object 2 class array type int rank 1 shape 4 items M
   s1a s1b s1c s1d
   s2a s2b s2c s2d
   ...
   sMa sMb sMc sMd
   attribute "element type" string "tetrahedra"
   object 3 class array type float rank 0 items N
   u1
   u2
   ...
   uN
   attribute "dep" string "positions"
   object "irregular positions irregular connections" class field
   component "positions" value 1
   component "connections" value 2
   component "data" value 3
   end

where the variables in this format are:

``N``
  Number of vertices

``vix viy viz``
  Coordinates of vertex i

``M``
  Number of simplices

``sia sib sic sid``
  IDs of vertices in simplex i

``ui``
  Data value associated with vertex i
File formats
============

----------------------------
Mesh and scalar data formats
----------------------------

.. toctree:: 
   :maxdepth: 1

   opendx
   mcsf
   uhbd

---------------------------
Molecular structure formats
---------------------------

.. toctree::
   :maxdepth: 1

   pqr
   pdb
   xml-struct

--------------
Matrix formats
--------------

.. toctree:: 
   :maxdepth: 1

   harwell

=================
Parameter formats
=================

.. toctree::
   :maxdepth: 1

   apbs-xml-parm
   flat-parm

.. _xmlstruct:

XML molecular structure format
==============================

The XML structure format was designed as a light-weight alternative to remediate some of the shortcomings of the flat-file format.
By use of XML, issues related to extra fields in the file or columns merging together can easily be remedied.
Additionally, APBS will only parse the necessary information from the XML file and will ignore all other information, so users wishing to store extra data related to a residue or atom can do so inline without affecting APBS.

This data format has the following form:

.. code-block:: xml

   <roottag>
      <residue>
          <atom>
              <x>x</x>
              <y>y</y>
              <z>z</z>
              <charge>charge</charge>
              <radius>radius</radius>
          </atom>
          ...
      </residue>
      ...
   </roottag>

The variables in this example are:

``roottag``
  This is the root element of the XML file. The value is not important to APBS - APBS simply checks that it is closed at the end of the file.

``x y z``
  A float giving the {x, y, z}-coordinate of the atom in Å.

``charge``
  A float giving the atomic charge (in electrons).

``atomradius``
  A float giving the atomic Radius (in Å).

.. note::

   Yes, we probably should have used `PDBML <http://pdbml.pdb.org/>`_ instead.
.. _apbsflatparm:

APBS flat-file parameter format
===============================

This parameter file format is a series of lines of the form:

.. code-block:: bash

   Residue_name Atom_name Charge Radius Epsilon

where the whitespaces are important and denote separation between the fields.
The fields here are:

``Residue_name``
  A string giving the residue name, as provided in the PDB file to be parametrized.

``Atom_name``
  A string giving the atom name, as provided in the PDB file to be parametrized.

``Charge``
  A float giving the atomic charge (in electrons).

``Radius``
  A float giving the atomic radius (in Å).

``Epsilon``
  A float giving the Lennard-Jones well depth (epsilon, in kJ/mol).
  This is used for the calculation of WCA energies in apolar solvation energies and forces.
  We assume that the Lennard-Jones potential is defined in the "AMBER style"
UHBD scalar data format
=======================

We also support scalar data output in the legacy "UHBD format" for use with programs such as `UHBD <http://browndye.ucsd.edu/>`_ and `SDA <https://mcm.h-its.org/sda/>`_.

UHBD data is written in the format:

.. code-block:: c

    /* Write out the header */
    Vio_printf(sock, "%72s\n", title);
    Vio_printf(sock, "%12.5e%12.5e%7d%7d%7d%7d%7d\n", 1.0, 0.0, -1, 0,
      nz, 1, nz);
    Vio_printf(sock, "%7d%7d%7d%12.5e%12.5e%12.5e%12.5e\n", nx, ny, nz,
      hx, (xmin-hx), (ymin-hx), (zmin-hx));
    Vio_printf(sock, "%12.5e%12.5e%12.5e%12.5e\n", 0.0, 0.0, 0.0, 0.0);
    Vio_printf(sock, "%12.5e%12.5e%7d%7d", 0.0, 0.0, 0, 0);

    /* Write out the entries */
    icol = 0;
    for (k=0; k<nz; k++) {
        Vio_printf(sock, "\n%7d%7d%7d\n", k+1, thee->nx, thee->ny);
        icol = 0;
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {
                u = k*(nx)*(ny)+j*(nx)+i;
                icol++;
                Vio_printf(sock, " %12.5e", thee->data[u]);
                if (icol == 6) {
                    icol = 0;
                    Vio_printf(sock, "\n");
                }
            }
        }
    }
.. _GitHub repository: https://github.com/Electrostatics/apbs
.. _GitHub releases: https://github.com/Electrostatics/apbs/releases

=============================
How to build APBS from source
=============================

These instructions assume that you have downloaded the source code from `GitHub releases`_.

.. caution:: We do not recommend cloning directly from the head of the master branch because it is typically under development and could be unstable. Unless you really know what you are doing, we advise you to skip the next step.

-------------------------------
Get source directly from Github
-------------------------------

Here are the commands to get the source directly from our `GitHub repository`_, 

.. code:: bash

   git clone https://github.com/Electrostatics/apbs
   cd apbs

-----------------
Shortcut to build
-----------------

There is a script that is used to build APBS in the Github Actions. You may want to use the file, :file:`.build.sh`, as a template for building APBS.

-----------------
Import submodules
-----------------

*As of v3.4.0:* Submodules are only used for Pybind11, so this step is only required if building the Python interface.

We are using Git submodules to manage various pieces of code.  To build the master branch, after cloning it, you will need to do the following from within the top of the source directory:

.. code:: bash

   git submodule init
   git submodule update

------------
Set up CMake
------------

From the top of the source directory, the basic commands for configuring the APBS build for CMake are

.. code:: bash

   mkdir build
   cd build
   # NOTE: This will be you $APBS_BUILD_DIR
   export APBS_BUILD_DIR=`echo $(PWD)`
   cmake ..

To see all the options you can run:

.. code:: bash

   cd $APBS_BUILD_DIR
   ccmake ..

Additional features can be built using the flags described below.

^^^^^^^^^^^^^^
Geometric flow
^^^^^^^^^^^^^^

If you want to use the geometric flow implementation, when invoking CMake, set :makevar:`ENABLE_GEOFLOW` to ``ON``; e.g.,

.. code:: bash

   cd $APBS_BUILD_DIR
   cmake -DENABLE_GEOFLOW=ON ..

^^^^^^^^^^^
Using PB-AM
^^^^^^^^^^^

If you want to use the Poisson-Boltzmann Analytical Method developed by the Teresa Head-Gordon lab, set the CMake variable :makevar:`ENABLE_PBAM` to ``ON``.

.. warning::

   PB-AM currently runs on OS X or Linux only.

.. code:: bash

   cd $APBS_BUILD_DIR
   cmake -DENABLE_PBAM=ON ..

^^^^^^^^^^^^^
Using TABI-PB
^^^^^^^^^^^^^

If you want to use the Treecode-Accelerated Boundary Integral method (TABI-PB) developed by Robert Krasny and Weihua Geng, set the CMake variable :makevar:`ENABLE_BEM` to ``ON``.

TABI-PB requires the use of a molecular surface mesh generation software to create a surface representation of the molecule.
By default, TABI-PB uses NanoShaper to generate an SES or Skin surface.
See `TABI-PB documentation <https://github.com/Treecodes/TABI-PB>`_ for details on choosing NanoShaper.
When TABI-PB runs, it will attempt to generate a surface mesh by looking in your path for the mesh generation executable.
A user can obtain the appropriate executable using the steps described below. The user then must place these executables in their path.

.. code:: bash

   cd $APBS_BUILD_DIR
   cmake -DENABLE_BEM=ON ..

"""""""""""""""""""""""""""""
Getting NanoShaper executable
"""""""""""""""""""""""""""""

Surface meshing software executables are currently pre-built for OS X, Linux, and Windows and can be installed via CMake.
The executables will be placed in the :file:`bin` directory of your build.

NanoShaper is a molecular surface mesh generation software package developed by W. Rocchia and S. Decherchi.

.. code:: bash

   cd $APBS_BUILD_DIR
   cmake -DGET_NanoShaper=ON ..

^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Using finite element support
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*As of v3.4.0:* The Finite Element Toolkit, FETK, is required for building APBS.  
You can set the version of FETK used with the CMake variable :makevar:`FETK_VERSION`.
That variable will be set to a working default if not manually set.

.. code:: bash

   cd $APBS_BUILD_DIR
   cmake -DENABLE_FETK_=ON -DFETK_VERSION=v1.9.2

For advanced users, you can use a version of FETK other than a released version by setting ``FETK_VERSION``
to the desired git commit hash instead of a version number:

.. code:: bash

   cd $APBS_BUILD_DIR
   cmake -DENABLE_FETK=ON -DFETK_VERSION=[git hash]


^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Enabling APBS Python support
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

APBS Python support requires a local installation of `SWIG <http://www.swig.org/>`_.

Assuming SWIG is installed, APBS Python support can be enabled by setting the CMake variable :makevar:`ENABLE_PYTHON` to ``ON``.
If you are on Linux you will also need to set the CMake variable :makevar:`BUILD_SHARED_LIBS` to ``OFF``.

.. code:: bash

   cd $APBS_BUILD_DIR
   cmake -DENABLE_PYTHON=ON ..

---------------------------
Building the code - minimal
---------------------------

Assuming the Cmake command completed successfully, APBS can be built with

.. code:: bash

   cd $APBS_BUILD_DIR
   # Run cmake with the options you prefer:
   make -j

----------------------------
Building the code - advanced
----------------------------

.. code:: bash

   export INSTALL_DIR=$SOME_DIR/apbs
   export PATH=$INSTALL_DIR/bin:$PATH
   # NOTE: In case you need to debug the source code:
   # export RELEASE_TYPE=Debug
   export RELEASE_TYPE=Release
   # NOTE: If cmake or make fail, save yourself and make sure your remove
   #       everything including the build directory. This code base uses
   #       many older autoconf based projects that do not know how to save
   #       state or recover from partial builds. If cmake or make fail, you
   #       should figure out how to fix it and then remove everything and
   #       try again.
   rmdir $APBS_BUILD_DIR
   mkdir -p $APBS_BUILD_DIR
   cd $APBS_BUILD_DIR
   # NOTE: In case you need to debug cmake, use verbose debug/trace mode:
   # cmake -S .. -B $BUILD_DIR --trace-source=../CMakeLists.txt --trace-expand \
   cmake                                        \
      -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR       \
      -DCMAKE_BUILD_TYPE=$RELEASE_TYPE          \
      -DENABLE_GEOFLOW=ON                       \
      -DENABLE_BEM=ON                           \
      -DENABLE_FETK=ON                          \
      -DFETK_VERSION=[version]                  \
      -DENABLE_OPENMP=ON                        \
      -DENABLE_PBAM=ON                          \
      -DENABLE_PBSAM=ON                         \
      -DENABLE_PYTHON=ON                        \
      -DENABLE_TESTS=ON                         \
      -DBUILD_SHARED_LIBS=ON                    \
      ..
   make -j

------------
Testing APBS
------------

.. code:: bash

   cd $APBS_BUILD_DIR
   # NOTE: Assuming you have already built APBS
   # NOTE: So that the apbs and optional NanoShaper binaries are in the path:
   export PATH="$APBS_BUILD_DIR/bin:$PATH"
   ctest -C Release --output-on-failure

---------------
Installing APBS
---------------

.. code:: bash

   export INSTALL_DIR="Some directory - default is /usr/local"
   cd $APBS_BUILD_DIR
   cmake                                  \
      -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
      # NOTE: Add cmake options that you used during the Build APBS section
   ..
   make -j install
.. _registering: http://eepurl.com/by4eQr
.. _GitHub releases: https://github.com/Electrostatics/apbs/releases
.. _Visual C++ Redistributable Package:  https://aka.ms/vs/17/release/vc_redist.x86.exe

============
Getting APBS
============

.. note::

   *Before you begin!* APBS funding is dependent on your help for continued development and support. Please `register <http://eepurl.com/by4eQr>`_ before using the software so we can accurately report the number of users to our funding agencies.

-----------
Web servers
-----------

Most functionality is available through our online web servers.

The web server offers a simple way to use both APBS and PDB2PQR without the need to download and install additional programs.

After `registering`_, please visit http://server.poissonboltzmann.org/ to access the web server.

-------------------------------------
Installing from pre-compiled binaries
-------------------------------------

The best way to install APBS is via download of a pre-compiled binary from `GitHub releases`_ (after `registering`_, of course).

^^^^^^^^^^^^
Requirements
^^^^^^^^^^^^

The pre-compiled binaries include nearly all dependencies, so the only requirement is Python:

* All platforms
  * Python 3 (tested with 3.9)

.. caution:: 

  On Linux you may need to set your LD_LIBRARY_PATH and PATH environment variables:
  For example, in bash with APBS installed in $HOME/apbs, you would set the environment variables in your .bashrc like:
  .. code-block:: bash

    export LD_LIBRARY_PATH=$HOME/apbs/lib:${LD_LIBRARY_PATH}
    export PATH=$HOME/apbs/bin:${PATH}

^^^^^^^^^^^^^^^^^^
What's in the box?
^^^^^^^^^^^^^^^^^^

The binary distributions typically provide the following contents:

bin
  contains the main APBS executable
include
  header files for building software that calls APBS
lib
  libraries for building software that calls APBS
share/apbs/docs
  the APBS documentation
share/apbs/examples
  APBS examples
share/apbs/tests
  the APBS test suite
share/apbs/tools
  useful programs to help process APBS input and output


---------------------------
Installing from source code
---------------------------

Those who enjoy an adventure can download the source code from `GitHub releases`_ and install from source code following the directions at the link below:

.. toctree::
   :maxdepth: 2

   source

------------------------
Current platform support
------------------------

+------------+---------+---------+------------+------+-------+------+---------+-------------+
| OS         | PYTHON  | GEOFLOW | BEM,       | FETK | PBSAM | PBAM | PYTHON  | SHARED_LIBS |
|            | VERSION |         | NanoShaper |      |       |      | SUPPORT |             |
+============+=========+=========+============+======+=======+======+=========+=============+
| Ubuntu     | 3.9     | Yes     | Yes        | Yes  | Yes   | Yes  | Yes     | No          |
+------------+---------+---------+------------+------+-------+------+---------+-------------+
| MacOSX     | 3.9     | Yes     | Yes        | Yes  | Yes   | Yes  | Yes     | No          |
+------------+---------+---------+------------+------+-------+------+---------+-------------+
| Windows 10 | 3.9     | Yes     | Yes        | Yes  | Yes   | Yes  | No      | No          |
+------------+---------+---------+------------+------+-------+------+---------+-------------+
.. _api-label:

=============
API reference
=============

.. currentmodule:: apbs
.. module:: apbs

The :program:`apbs` command provides a command-line interface to APBS's functionality.
It is built on classes and functions in the :mod:`apbs` module.
The API (application programming interface) of :mod:`apbs` is documented here for developers who might want to directly use the APBS code.

.. Note::

   The API is still changing and there is currently no guarantee that
   it will remain stable between minor releases.

.. todo:: Add API documentation for all new Python routines
.. _units:

The following are preferred APBS units.

=======
Numbers
=======

Many quantities are scaled by Avogadro's number, the number of atoms in a mole (mol), :math:`N_A = 6.02214076 \times 10^{23} \, \text{mol}^{-1}`.

======
Length
======

The preferred unit of length is the ångström (Å), equal to 10\ :sup:`-10` meters.

======
Volume
======

The preferred unit of volume is Å\ :sup:`3`, equal to 10\ :sup:`-27` liters.

=========================
Density and concentration
=========================

The preferred unit of density is number per Å\ :sup:`3`, corresponding to a
concentration of approximately 1660.5391 mol L\ :sup:`-1` or molar (M).

===========
Temperature
===========

The preferred unit of temperature is Kelvin (K).

======
Charge
======

The preferred unit of charge is :math:`e_c = 1.602176634 \times 10^{−19} \, \text{C}`.
The following number is often useful:  :math:`N_A e_c = 9.64853321233 \times 10^4 \, \text{C} \, \text{mol}^{-1}`.

======
Energy
======

The preferred unit of energy is :math:`k_B \, T` or :math:`R \,T` where

* Boltzmann's constant: :math:`k_B = 1.38064852 \times 10^{-23} \, \text{J} \, \text{K}^{-1}`
* Gas constant: :math:`R = N_A k_B = 8.31446261815324 \text{J} \, \text{K}^{-1} \, \text{mol}^{-1}`
* Temperature: :math:`T`

If :math:`T \approx 298 \, \text{K}`, then :math:`R\, T \approx 2.49 \, \text{kJ}`.

===============
Surface tension
===============

The preferred unit for surface tension is kJ mol\ :sup:`-1` Å\ :sup:`-2`.
Values for the surface tension of water in these models often range from 0.105 to 0.301 kJ mol\ :sup:`-1` Å\ :sup:`-2`. [#Sharp]_
However, these values can vary significantly depending on the model used. [#Thomas]_ [#Wagoner]_

========
Pressure
========

The preferred unit for pressure is kJ mol\ :sup:`-1` Å\ :sup:`-3`.
Values for the surface tension of water in these models vary significantly depending on the model used (e.g., between 0.0004 and 0.146 kJ mol\ :sup:`-1` Å\ :sup:`-3`). [#Thomas]_ [#Wagoner]_

=======================
Electrostatic potential
=======================

The preferred unit of electrostatic potential is :math:`k_B \, T \, e_c^{-1}` or :math:`R \, T \, e_c^{-1}`.
If :math:`T \approx 298 \, \text{K}`, then :math:`k_B \, T \, e_c^{-1} = R \, T \, e_c^{-1} \approx 0.0256 \, \text{J} \, \text{C}^{-1} = 25.6 \, \text{mV}`.

.. [#Sharp] Sharp KA, Nicholls A, Fine RF, Honig B. Reconciling the magnitude of the microscopic and macroscopic hydrophobic effects. Science, 252, 106-109, 1991. DOI:`10.1126/science.2011744 <http://doi.org/10.1126/science.2011744>`_.

.. [#Thomas] Thomas DG, Chun J, Zhen C, Wei GW, Baker NA. Parameterization of a geometric flow implicit solvation model. J Comput Chem, 34, 687-695, 2013. DOI:`10.1002/jcc.23181 <https://doi.org/10.1002/jcc.23181>`_.

.. [#Wagoner] Wagoner JA and Baker NA.  Assessing implicit models for nonpolar mean solvation forces: The importance of dispersion and volume terms. Proc Natl Acad Sci USA, 103, 8331-9336, 2006. DOI:`10.1073/pnas.0600118103 <https://doi.org/10.1073/pnas.0600118103>`_.==============================
APBS validation and test cases
==============================

This directory serves as the root directory for the APBS test suite.  This
directory contains python source files used for testing an input file
containing the input files used by the apbs executable and the expected results for each test case.

The default input file is called `test_cases.cfg`, and the main testing program is called
`apbs_tester.py`.

-----
Usage
-----

.. note::

  It is important that you run the tests from the command line and that you run them from within the `INSTALL_DIR/share/apbs/tests` directory.


A usage description for `apbs_tester.py` can be obtained by running:

.. code-block:: console

   python3 apbs_tester.py -h

-------------
Test Sections
-------------

The sections of the test file, :file:`test_cases.cfg`, follow the following format::

  [Some-Target_Test]
  input_dir     : ../path/to/some-example
  some-forces   : forces
  some-input    : * * 1.0E+01 2.0E+02

where:

* The first element in brackets :mailheader:`[Some-Target_Test]` describes the name of the *target_test* section.

* After the first element, the remaining elements are *property*/*value* pairs

* The first property is the :makevar:`input_dir`.
  This is the location of all input files referenced in other properties.

* A property has a *name* that is also the basename of the input file concatenated with the file extension, :file:`{input-file}.in`.

* The property *name* will also be used for the output from :file:`apbs some-input.in` to create the output file, :file:`some-input.out`

* If the *value* of the property is :makevar:`forces` the test will calcalate forces.

* If the *value* of the property is a list of floating point numbers, the values are expected outputs.

* If a :makevar:`*` is used in place of a floating point number, the output will be ignored. Some test cases have multiple outputs.
  The test function parses each of these, but if a :makevar:`*` is used, the output will be ignored in testing.
  Most often, the first outputs are intermediate values followed by a final output, and the test case is only concerned with the final output.

--------
Examples
--------

The following will run the `apbs_tester.py` specifying the path to the apbs executable and using the geoflow section of the test_cases.cfg file:

.. code-block:: console

   python3 apbs_tester.py -e ${INSTALL_DIR}/bin/apbs -c test_cases.cfg -t geoflow

where :makevar:`INSTALL_DIR` is the path to your installation directory.
==========
Using APBS
==========

.. _PDB ID: https://www.rcsb.org/pages/help/advancedsearch/pdbIDs

.. note::

   *Before you begin!* PDB2PQR funding is dependent on your help for continued development and support. Please `register <http://eepurl.com/by4eQr>`_ before using the software so we can accurately report the number of users to our funding agencies.

APBS is often used together with the `PDB2PQR software <https://github.com/Electrostatics/pdb2pqr>`_; e.g., ,in the following type of workflow

#. Start with a `PDB ID`_ or locally generated PDB file (see :doc:`/formats/pdb`).
#. Assign titration states and parameters with :program:`pdb2pqr` to convert the protein and ligands to PQR format (see :doc:`/formats/pqr`).
#. Perform electrostatics calculations with :program:`apbs` (can be done from within the `PDB2PQR web server <web-server>`_).
#. Visualize results from within PDB2PQR web server or with :ref:`other-software`.

--------------
Web server use
--------------

Most users will use PDB2PQR through `the web server <http://server.poissonboltzmann.org/>`_ (after `registering <http://eepurl.com/by4eQr>`_, of course).
However, it is also possible to install local versions of PDB2PQR and run these through the command line.

----------------
Command line use
----------------

.. code-block:: bash

   apbs [options] input-file

where the list of ``[options]`` can be obtained by running APBS with the ``--help`` option.
The input file format is described below.

-----------------
Input file syntax
-----------------

APBS has a :ref:`new input file format <new_input_format>` that accepts `YAML- <http://yaml.org>`_ or `JSON- <http://json.org>`_ format input.
The :ref:`old APBS input format <old_input_format>` has been deprecated but will continue to be supported for the next few releases.

.. toctree::
   :maxdepth: 2
   :caption: New and old input file formats

   input/new/index
   input/old/index

.. _examples:

--------
Examples
--------

.. _PDB2PQR:  https://github.com/Electrostatics/pdb2pqr
.. _APBS-PDB2PQR web interface:  http://server.poissonboltzmann.org/

APBS examples start with a PQR file (e.g., generated by PDB2PQR_).
Some of these examples can be performed through the `APBS-PDB2PQR web interface`_ but most require a command-line APBS program.

.. toctree::
   :maxdepth: 2

   examples/solvation-energies
   examples/binding-energies
   examples/salt-linkage
   examples/visualization-pymol
   examples/visualization-unitymol

--------------------
Tests and validation
--------------------

APBS is distributed with testing tools and validation examples.

.. toctree::
   :maxdepth: 1

   tests

-------------------
Tools and utilities
-------------------

APBS is distributed with utilities designed to simplify typical tasks associated with electrostatics calculations.

.. toctree::
   :maxdepth: 1

   tools

.. _other-software:

--------------
Other software
--------------

A variety of other software can be used to visualize and process the results of PDB2PQR and APBS calculations.

^^^^^^^^^^^^^^^^^^^^^^
Visualization software
^^^^^^^^^^^^^^^^^^^^^^

Examples of visualization software that work with output from PDB2PQR and APBS:

* `PyMOL <https://pymol.org/>`_
* `VMD <https://www.ks.uiuc.edu/Research/vmd/>`_
* `Chimera <https://www.cgl.ucsf.edu/chimera/>`_
* `PMV <http://mgltools.scripps.edu/packages/pmv>`_

^^^^^^^^^^^^^^^^^^^^
Dynamics simulations
^^^^^^^^^^^^^^^^^^^^

As an example of PDB2PQR and APBS integration with molecular mechanics sofware, the `iAPBS <https://mccammon.ucsd.edu/iapbs/>`_ library was developed to facilitate the integration of APBS with other molecular simulation packages.
This library has enabled the integration of APBS with several molecular dynamics packages, including `NAMD <http://www.ks.uiuc.edu/Research/namd/>`_, `AMBER <http://ambermd.org/>`_, and `CHARMM <https://www.charmm.org/charmm/>`_.

APBS is also used directly by Brownian dynamics software such as `SDA <https://mcm.h-its.org/sda/>`_ and `BrownDye <http://browndye.ucsd.edu/>`_.

--------------
Notes on units
--------------

APBS and PDB2PQR use a few different sets of units, explained in the following sections:

.. toctree::
   :maxdepth: 2

   units--------------------
Conversion utilities
--------------------

.. todo::  Update this documentation with the *new APBS syntax* (see :ref:`new_input_format`).

^^^^^^^^^^^^^^^
amber2charmm.sh
^^^^^^^^^^^^^^^

A bash script for converting AMBER atom names to CHARMM names.
Found in :file:`tools/conversion`

^^^^^^
del2dx
^^^^^^

Converts DelPhi-format map files (electrostatic potential, etc.) to APBS OpenDX format.
Found in :file:`tools/mesh`

^^^^^^
dx2mol
^^^^^^

Converts OpenDX format map files to MolMol format.
Found in :file:`tools/mesh`

^^^^^^^
dx2uhbd
^^^^^^^

Converts OpenDX format map files to UHBD format.
Found in :file:`tools/mesh`

^^^^^^^^^^^
qcd2pqr.awk
^^^^^^^^^^^

An awk script for converting from UHBD QCD format to PQR format.

----------------------
Benchmarking utilities
----------------------

^^^^^^^^^
benchmark
^^^^^^^^^

Benchmark file I/O for reading/writing scalar data.
Found in :file:`tools/mesh`

^^^^^^^^^^^^
uhbd_asc2bin
^^^^^^^^^^^^

Converts UHBD ASCII-format files to binary format.
Found in :file:`tools/mesh`

^^^^^^^^^^^^^^^^
WHATIF2AMBER.sed
^^^^^^^^^^^^^^^^

A sed script for converting WHATIF atoms names to the AMBER naming scheme.
Found in :file:`tools/conversion`

----------------------------
Setup and analysis utilities
----------------------------

^^^^^^^^
analysis
^^^^^^^^

Calculates various metrics from input scalar data.
Found in :file:`tools/mesh`

^^^^
born
^^^^

Calculate generalized Born forces and energies.
Found in :file:`tools/manip`

^^^^^^^
coulomb
^^^^^^^

Calculate Coulomb forces and energies.
Found in :file:`tools/manip`

.. _dxmath:

^^^^^^
dxmath
^^^^^^

Performs simple arithmetic operations with Cartesian grid data.  
This program takes as input a file with operations specified in a stack-based (RPN) manner.
For example, a command file which adds grid1 and grid2, multiplies the result by 5.3, adds grid4, subtracts 99.3 from the whole thing, and writes the result on grid5 would have the form:

.. code-block:: mathematica
   
   grid1
   grid2 +
   5.3 *
   grid4 +
   99.3 -
   grid5 =

The file names, scalar values, and operations must be separated by tabs, line breaks, or white space.
Comments can be included between the character # and a new line (in the usual shell script fashion).
Found in :file:`tools/mesh`

^^^^^^^^^^^
inputgen.py
^^^^^^^^^^^

Create an APBS input file using :ref:`psize` data.
Found in :file:`tools/manip`

^^^^^^^^^^^^^^^^^^^^
mergedx and mergedx2
^^^^^^^^^^^^^^^^^^^^

Combine multiple OpenDX files into a single resampled file.
:program:`mergedx2` can perform a number of grid manipulation operations, including:

* Combining multiple OpenDX map files
* Resampling of one or more OpenDX map files (for example to alter the grid spacing of separate OpenDX files for further manipulation)
* Extracting a subregion of an existing OpenDX map file.

Found in :file:`tools/mesh`

^^^^^^
mgmesh
^^^^^^

Prints out acceptable combinations of :doc:`input/old/elec/nlev` and :doc:`input/old/elec/dime` for multigrd calculations.
Found in :file:`tools/mesh`

^^^^^^^^^^
multivalue
^^^^^^^^^^

This program evaluates OpenDX scalar data at a series of user-specified points and returns the value of the data at each point.
Found in :file:`tools/mesh`

.. _psize:

^^^^^^^^
psize.py
^^^^^^^^

Suggest grid sizes and spacings for APBS given an input molecule.
Found in :file:`tools/manip`

^^^^^^^^^^
similarity
^^^^^^^^^^

Computes similarity between two scalar grid datasets.
Found in :file:`tools/mesh`

^^^^^^
smooth
^^^^^^

Convolve grid data with various filters.
Found in :file:`tools/mesh`

.. _print:

PRINT input file section
========================

.. note::  

   This section has been ported to the *new APBS syntax* (see :ref:`new_input_format`).
   See :ref:`process_new_input` for more information.

This is a very simple section that allows linear combinations of calculated properties to be written to standard output.
The syntax of this section is:

.. code-block:: bash

    PRINT {what} [id op id op...] END

The first mandatory argument is ``what``, the quantity to manipulate or print.
This variable is a string that can assume the following values:

``elecEnergy``
  Print electrostatic energies as calculated with an earlier :ref:`elec` :ref:`calcenergy` command.
``elecForce``
  Print electrostatic forces as calculated with an earlier :ref:`elec` :ref:`calcforce` command.
``apolEnergy``
  Print apolar energies as calculated with an earlier :ref:`apolar` :ref:`calcenergy` command.
``apolForce``
  Print electrostatic forces as calculated with an earlier :ref:`apolar` :ref:`calcforce` command.

The next arguments are a series of ``id op id op id op ... id`` commands where every ``id`` is immediately followed by an ``op`` and another ``id``.

``id``
  This is a variable string or integer denoting the ID of a particular :ref:`elec` or :ref:`apolar` calculations.
  String values of ``id`` correspond to the optional "names" that can be assigned to :ref:`elec` or :ref:`apolar` calculations.
  Integer values of id are assumed to corresponding to the sequentially-assigned integer IDs for :ref:`elec` or :ref:`apolar` calculations.
  These IDs start at 1 and are incremented (independently) for each new :ref:`elec` or :ref:`apolar` calculation.
``op``
  Specify the arithmetic operation (``+`` for addition and ``-`` for subtraction) to be performed on the calculated quantities

For example:

.. code-block:: python

   # Energy change due to binding
   print energy complex - ligand - protein end
   # Energy change due to solvation
   print energy solvated - reference end
   # Solvation energy change due to binding
   print energy complex_solv - complex_ref - ligand_solv + ligand_ref - protein_solv + protein_ref end

.. _old_input_format:

=====================
Old APBS input format
=====================

In the old format, APBS input files are loosely-formatted files which contain information about the input, parameters, and output for each calculation.

These files are whitespace- or linefeed-delimited.
Comments can be added to the input files via the ``#`` character; all text between the ``#`` and the end of the line is not parsed by APBS.
If pathnames used in the input file contain spaces, then the entire pathname must be enclosed in quotes.
For example, if you wanted to refer to the file :file:`foo` which resides in a directory with spaces in its name, then you should refer to :file:`foo` as :file:`"/path with spaces/foo"`.
Specific examples of APBS input are provided in :ref:`examples`.

APBS input files contain three basic sections which can be repeated any number of times:

* :ref:`read_old_input`:
  Section for specifying data-reading input.
  For the *new APBS syntax*, see :ref:`read_new_input`.
* :ref:`elec`:
  Section for specifying polar solvation (electrostatics) calculation parameters.
  For the *new APBS syntax*, see :ref:`calculate_new_input`.
* :ref:`apolar`:
  Section for specifying apolar solvation calculation parameters.
  For the *new APBS syntax*, see :ref:`nonpolar_new_input`.
* :ref:`print`:
  Section for specifying summary output.
  For the *new APBS syntax*, see :ref:`process_new_input`.

The APBS input file is constructed from these sections in the following format:

.. code-block:: bash
   
   READ
   ...
   END
   
   ELEC
   ...
   END
   
   APOLAR
   ...
   END
   
   PRINT
   ...
   END
   
   QUIT

These sections can occur in any order and can be repeated any number of times.
However, the sections are inter-dependent.
For example, PRINT requires ELEC and/or APOLAR while ELEC requires one or more READ sections.
Sections can also be repeated; several READ statements may be used to load molecules and multiple ELEC or APOLAR sections would specify various electrostatics calculations on one or more molecules.

Each section has the following syntax:

.. code-block:: bash
   
   SECTION [name <id>]

where the optional ``name`` argument allows the user to include a string to identify the section.
In the absence of this argument, sections are assigned numerical IDs.

.. toctree::
   :maxdepth: 1
   :caption: Input file sections

   read
   elec/index
   apolar/index
   print
.. _read_old_input:

READ input file section
=======================

.. note::  

   This section has been ported to the *new APBS syntax* (see :ref:`new_input_format`).
   See :ref:`read_new_input` for more information.

The READ block of an APBS input file has the following general format:

.. code-block:: bash
   
   READ
       [ keywords... ]
   END

where ``keywords`` is or more of the keywords described below (the line breaks and indentation are for clarity; only whitespace is necessary).

.. note::
   One of these sections must be present for every molecule involved in the APBS calculation.
   Molecule and "map" IDs are assigned implicitly assigned for each molecule/map read, based on order and starting at 1 and incremented independently for each input type.
   In other words, each input PQR file is assigned an ID 1, 2, 3, ...; each input dielectric map is assigned an independent ID 1, 2, 3, ...; etc.

------
charge
------

This command allows APBS to read the fixed (molecular) charge density function mapped to a mesh.
The inputs are maps of charge densities; these values have units of e\ :sub:`c` Å\ :sup:`-3`, where e\ :sub:`c` is the electron charge.
In general, this command will read charge-maps written by :ref:`elec` :ref:`old_write` commands.
The syntax of this command is:

.. code-block:: bash

   READ charge {format} {path} END 

``format``
  Specify the format of the charge map.
  Acceptable values include:

  ``dx``
    :ref:`opendx`

  ``gz``
    gzipped (zlib) compressed :ref:`opendx`.
    Files can be read directly in compressed form.

``path``
  The location of the charge map file.


----
diel
----

This command allows APBS to read the dielectric function mapped to 3 meshes shifted by one-half grid spacing in the x, y, and z directions.
The inputs are maps of dielectric variables between the solvent and biomolecular dielectric constants; these values are unitless.
In general, this command will read dielectric maps written by by :ref:`elec` :ref:`old_write` commands.
The syntax of this command is:

.. code-block:: bash

   READ diel {format} {path-x} {path-y} {path-z} END

``format``
  The format of the dielectric map.

  ``dx``
    :ref:`opendx`
  
  ``gz``
    gzipped (zlib) compressed :ref:`opendx`.
    Files can be read directly in compressed form.

``path-x``
  The location of the x-shifted dielectric map file.

``path-y``
  The location of the y-shifted dielectric map file.

``path-z`` The location of the z-shifted dielectric map file.

.. note::

   If you choose this option and have a non-zero ionic strength, you must also include a READ kappa_ statement.

-----
kappa
-----

This command allows APBS to read the ion-accessibility function mapped to a mesh.
The inputs are maps of ion accessibility values which range between 0 and the build Debye-Hückel screening parameter; these values have units of Å\ :sup:`-2`.
In general, this command will read kappa-maps written by by :ref:`elec` :ref:`old_write` commands.
The syntax of this command is:

.. code-block:: bash

   READ kappa {format} {path} END

``format``
  Specify the format of the charge map.
  Acceptable values include:

  ``dx``
    :ref:`opendx`
  
  ``gz``
    gzipped (zlib) compressed :ref:`opendx`.
    Files can be read directly in compressed form.

``path``
  The location of the map file.


.. note::

   If you choose this option, you must also include a read diel statement.

---
mol
---

This command specifies the molecular data to be read into APBS.
The syntax is

.. code-block:: bash

   READ mol {format} {path} END

``format``
  The format of the input data.

  ``pqr``
    Specify that molecular data is in :ref:`PQR format <pqr>`.
  
  ``pdb``
    Specify that molecular data is in pseudo-PDB format.
    If this type of structure file is used, then a parameter file must also be specified with a READ parm_ statement to provide charge and radius parameters for the biomolecule's atoms.

``path``
  The location of the molecular data file.

----
parm
----

This command specifies the charge and radius data to be used with pseudo-PDB-format molecule files.
The syntax is:

.. code-block:: bash

   READ parm {format} {path} END

``format``
  The format of the parameter file.

  ``flat``
    Specify that the parameter file is in :ref:`APBS flat-file parameter format <apbsflatparm>`.

  ``xml``
    Specify that the parameter file is in :ref:`APBS XML parameter format <apbsxmlparm>`

``path``
  The location of the parameter data file.

.. note::
   
   APBS provides a few example files as part of the source code distribution.
   Currently, example files only contain the polar parameters that can also be assigned more easily through the PDB2PQR software.

---
pot
---

This command allows APBS to read the electrostatic potential mapped to a mesh.
The inputs are maps of the electrostatic potential from a previous calculation.
In general, this command will read potential-maps written by by :ref:`elec` :ref:`old_write` commands.
The syntax of this command is:

.. code-block:: bash

   READ pot {format} {path} END

``format``
  Specify the format of the charge map.
  Acceptable values include:

  ``dx``
    :ref:`opendx`

  ``gz``
    gzipped (zlib) compressed :ref:`opendx`.
    Files can be read directly in compressed form.

``path``
  The location of the map file.

.. note::
   
   To use this functionality you must set the :ref:`bcfl` keyword to ``map``.
   See also: :ref:`usemap`.
.. _swin:

swin
====

.. currentmodule::  apbs.input_file.calculate

.. note::


   Some instances of this keyword have been moved to the *new APBS syntax* (see :ref:`new_input_format`):

   * For finite difference calculations, see :func:`finite_difference.FiniteDifference.surface_spline_window`
   * For finite element calculations, see :func:`finite_element.FiniteElement.surface_spline_window`

   .. todo::  move other instances of this keyword to the new syntax

Specify the size of the support (i.e., the rate of change) for spline-based surface definitions (see :ref:`elecsrfm`).
The syntax is:

.. code-block:: bash
   
   swin {win}

where ``win`` is a floating point number for the spline window width (in Å).
Usually 0.3 Å.

Note that, per the analysis of Nina, Im, and Roux (`article <http://dx.doi.org/10.1016/S0301-4622(98)00236-1>`_)</a>, the force field parameters (radii) generally need to be adjusted if the spline window is changed.
.. _calcenergy:

calcenergy
==========

.. currentmodule:: apbs.input_file.calculate

.. note::  

   Some versions of this command have been ported to the *new APBS syntax* (see :ref:`new_input_format`):


   * Nonpolar calculations:  see :func:`nonpolar.Nonpolar.calculate_energy` for more information.
   * Polar calculations:

     * Boundary elements:  see :func:`boundary_element.BoundaryElement.calculate_energy`
     * Finite difference:  see :func:`finite_difference.FiniteDifference.calculate_energy`
     * Finite element:  see :func:`finite_element.FiniteElement.calculate_energy`

This optional keyword controls energy output from an apolar solvation calculation.
The syntax is:

.. code-block:: bash

   calcenergy <flag>

where ``flag`` is a string denoting what type of energy to calculate:

``no``
  (Deprecated) Don't calculate any energies.
``total``
  Calculate and return total apolar energy for the entire molecule.
``comps``
  Calculate and return total apolar energy for the entire molecule as well as the energy components for each atom.

.. note::
   This option must be used consistently (with the same ``flag`` value) for all calculations that will appear in subsequent :ref:`print` statements.
.. _temp:

temp
====

.. note::  

   Some versions of this command have been ported to the *new APBS syntax* (see :ref:`new_input_format`):

   .. currentmodule:: apbs.input_file.calculate

   * Boundary-element Poisson-Boltzmann polar calculations: See :func:`boundary_element.BoundaryElement.temperature`.
   * Nonpolar calculations: See :func:`nonpolar.Nonpolar.temperature`.
   * Finite-difference Poisson-Boltzmann polar calculations: See :func:`finite_difference.FiniteDifference.temperature`.
   * Finite-element Poisson-Boltzmann polar calculations: See :func:`finite_element.FiniteElement.temperature`.

   .. todo:: add other uses to new syntax

This keyword specifies the temperature for the calculation.
The syntax is:

.. code-block:: bash

   temp {T}

where ``T`` is the floating point value of the temperature for calculation (in K).
.. _grid:

grid
====

.. note::  

   Some versions of this command have been ported to the *new APBS syntax* (see :ref:`new_input_format`):


   * Nonpolar calculations:
      .. currentmodule:: apbs.input_file.calculate.nonpolar

      See :func:`Nonpolar.grid_spacings` for more information.

Specify the grid spacings for multigrid and volume integral calculations.
This value may be different in each direction.
The syntax is:

.. code-block:: bash

   grid {hx hy hz}

where ``hx hy hz`` are the (floating point) grid spacings in the x-, y-, and z-directions (respectively) in Å.
.. _srad:

srad
====

.. note::  

   .. currentmodule:: apbs.input_file.calculate

   Some versions of this command have been ported to the *new APBS syntax* (see :ref:`new_input_format`):

   * Boundary-element Poisson-Boltzmann polar calculations:  see :func:`boundary_element.Mesh.solvent_radius` for more information.
   * Finite-difference Poisson-Boltzmann polar calculations:  see :func:`finite_difference.FiniteDifference.solvent_radius` for more information.
   * Finite-element Poisson-Boltzmann polar calculations:  see :func:`finite_element.FiniteElement.solvent_radius` for more information.
   * Nonpolar calculations:  see :func:`nonpolar.Nonpolar.solvent_radius` for more information.

This keyword specifies the radius of the solvent molecules; this parameter is used to define various solvent-related surfaces and volumes (see :ref:`elecsrfm`).
This value is usually set to 1.4 Å for a water-like molecular surface and set to 0 Å for a van der Waals surface.
The syntax is:

.. code-block:: bash

   srad {radius}

where ``radius`` is the floating point value of the solvent radius (in Å).
This keyword is ignored for ``srfm spl2`` (see :ref:`elecsrfm`).

.. _mol:

mol
===

.. note::  

   Some versions of this command have been ported to the *new APBS syntax* (see :ref:`new_input_format`):

   .. currentmodule::  apbs.input_file.calculate

   * Nonpolar calculations: see :func:`nonpolar.Nonpolar.molecule` for more information.
   * Boundary-element Poisson-Boltzmann calculations: see :func:`boundary_element.BoundaryElement.molecule` for more information.
   * Finite-difference Poisson-Boltzmann calculations: see :func:`finite_difference.FiniteDifference.molecule` for more information.
   * Finite-element Poisson-Boltzmann calculations: see :func:`finite_element.FiniteElement.molecule` for more information.

.. todo:: port for other calculation types

This term specifies the molecule for which the calculation is to be performed.
The syntax is:

.. code-block:: bash
   
   mol {id}
   

where ``id`` is the integer ID of the molecule for which the apolar calculation is to be performed.
The molecule IDs are based on the order in which molecules are read by ``READ mol`` statements (see :ref:`read_old_input`), starting from 1.
.. _sdens:

sdens
=====

.. currentmodule:: apbs.input_file.calculate

.. note::  

   Some versions of this command have been ported to the *new APBS syntax* (see :ref:`new_input_format`):

   * Boundary-element Poisson-Boltzmann polar calculations:  see :func:`boundary_element.Mesh.surface_density` for more information.
   * Nonpolar calculations:  see :func:`nonpolar.Nonpolar.surface_density` for more information.

This keyword specifies the number of quadrature points per Å\ :superscript:`2` to use in calculation surface terms (e.g., molecular surface, solvent accessible surface).
This keyword is ignored when :ref:`srad` is 0.0 (e.g., for van der Waals surfaces) or when :ref:`elecsrfm` is ``spl2`` (e.g., for spline surfaces).
The syntax is:

.. code-block:: bash

   sdens {density}

where ``density`` is a floating point number indicating the number of grid points per Å\ :superscript:`-2`.
A typical value is 10.0.

.. note::
   There is a strong correlation between the value used for the sphere density, the accuracy of the results, and the APBS calculation time.
.. _bconc:

bconc
=====

.. note::  

   Some versions of this command have been ported to the *new APBS syntax* (see :ref:`new_input_format`):


   * Nonpolar calculations:
      .. currentmodule:: apbs.input_file.calculate.nonpolar

      See :func:`Nonpolar.solvent_density` for more information.

This keyword specifies the bulk solvent density.
This coefficient multiplies the integral term of the apolar model discussed above and can be set to zero to eliminate integral contributions to the apolar solvation calculation.
The syntax is:

.. code-block:: bash

   bconc <density>

where ``density`` is a floating point number giving the bulk solvent density in Å\ :superscript:`-3`.
.. _calcforce:

calcforce
=========

.. currentmodule:: apbs.input_file.calculate

.. note::  

   Some versions of this command have been ported to the *new APBS syntax* (see :ref:`new_input_format`):


   * Nonpolar calculations:  see :func:`nonpolar.Nonpolar.calculate_forces` for more information.
   * Polar calculations:

     * Boundary element:  see :func:`boundary_element.BoundaryElement.calculate_forces`
     * Finite difference:  see :func:`finite_difference.FiniteDifference.calculate_forces`
     * Finite element:  see :func:`finite_element.FiniteElement.calculate_forces`

This optional keyword controls energy output from an apolar solvation calculation.
The syntax is:

.. code-block:: bash

   calcforce {flag}

where ``flag`` is a text string that specifies the types of force values to be returned:

``no``
  (Deprecated) don't calculate any forces.
``total``
  Calculate and return total electrostatic and apolar forces for the entire molecule.
``comps``
  Calculate and return total electrostatic and apolar forces for the entire molecule as well as force components for each atom.

The possible outputs from calcforce are:

``tot {n}``
  total force for atom *n*
``qf {n}``
  fixed charge force for atom *n*
``db {n}``
  dielectric boundary force for atom *n*
``ib {n}``
  ionic boundary force for atom *n*

The values will be printed in three columns which correspond to the x, y, and z components of the force vector.

.. note::
   This option must be used consistently (with the same ``flag`` value) for all calculations that will appear in subsequent :ref:`print` statements.
.. _dpos:

dpos
====

.. currentmodule:: apbs.input_file.calculate.nonpolar

.. note::  

   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`):  see :func:`Nonpolar.displacement` for more information.

This is the displacement used for finite-difference-based calculations of surface area derivatives.
I know, this is a terrible way to calculate surface area derivatives -- we're working on replacing it with an analytic version.
In the meantime, please use this parameter with caution.
If anyone has code for a better method, please share!

The syntax is:

.. code-block:: bash

   dpos {displacement}

where ``displacement`` is a floating point number indicating the finite difference displacement for force (surface area derivative) calculations in units of Å.

.. warning::
   This parameter is very dependent on ``sdens`` (see :doc:`../generic/sdens`); e.g., smaller values of ``dpos`` require larger values of ``sdens``.
.. _apolar:

APOLAR input file section
=========================

.. note::  

   This section has been ported to the *new APBS syntax* (see :ref:`new_input_format`).
   See :ref:`nonpolar_new_input` for more information.

This section is the main component for apolar solvation calculations in APBS runs.
There may be several APOLAR sections, operating on different molecules or using different parameters for multiple runs on the same molecule.
The syntax of this section is:

.. code-block:: bash

   APOLAR [name id]
     <keywords...>
   END

The first (optional) argument is:

.. code-block:: bash

   name <id>

where ``id`` is a unique string which can be assigned to the calculation to facilitate later operations (particularly in the :doc:`../print` statements).
The ``keywords...`` describing the parameters of the apolar calculation are discussed in more detail below:

.. toctree::
   :maxdepth: 2
   :caption: APOLAR keywords:

   ../generic/bconc
   ../generic/calcenergy
   ../generic/calcforce
   dpos
   gamma
   ../generic/grid
   ../generic/mol
   press
   ../generic/sdens
   ../generic/srad
   srfm
   ../generic/swin
   ../generic/temp

APBS apolar calculations follow the very generic framework described in  Wagoner JA, Baker NA. Assessing implicit models for nonpolar mean solvation forces: the importance of dispersion and volume terms. Proc Natl Acad Sci USA, 103, 8331-8336, 2006. doi:`10.1073/pnas.0600118103 <http://dx.doi.org/10.1073/pnas.0600118103>`_.

Nonpolar solvation potentials of mean force (energies) are calculated according to:

.. math::

   {W}^{(\mathrm{np})}(x) = \gamma A(x) + pV(x) + \bar \rho \sum^N_{i=1} \int _{\Omega} u_i^{(\mathrm{att})} (x_i, y) \theta (x,y) \, \mathrm{d}y 

and mean nonpolar solvation forces are calculated according to:

.. math::

   \mathbf{F}_i^{(\mathrm{np})}(x) = -\gamma \frac{\partial A (x)}{\partial x_i} - p \int _{\Gamma _i (x)} \frac{y-x_i}{\lVert y - x_i \rVert} \, \mathrm{d}y - \bar \rho \sum _{i=1}^N \int _{\Omega} \frac{\partial u_i^{(\mathrm{att})}(x_i,y)}{\partial x_i} \theta (x,y) \, \mathrm{d}y 

In these equations, :math:`\gamma` is the repulsive (hard sphere) solvent surface tension (see :ref:`gamma`), *A* is the conformation-dependent solute surface area (see :ref:`srad` and :ref:`apolarsrfm` keywords), *p* is the repulsive (hard sphere) solvent pressure (see :ref:`press` keyword), *V* is the conformation-dependent solute volume (see :ref:`srad` and :ref:`apolarsrfm` keywords), :math:`\rho` (see :ref:`bconc` keywords) is the bulk solvent density, and the integral involves the attractive portion (defined in a Weeks-Chandler-Andersen sense) of the Lennard-Jones interactions between the solute and the solvent integrated over the region of the problem domain outside the solute volume *V*.
Lennard-Jones parameters are taken from APBS parameter files as read in through an APBS input file READ statement (see :ref:`read_old_input`).

.. note::

   The above expressions can easily be reduced to simpler apolar solvation formalisms by setting one or more of the coefficients to zero through the keywords.

.. warning::

   All APOLAR calculations require a parameter file which contains Lennard-Jones radius and well-depth parameters for all the atoms in the solute PDB.
   This parameter file must also contain radius and well-depth parameters for water (specifically: residue "WAT" and atom "OW").
   Complete parameter files for protein and nucleic acid parameters are not currently available; we prefer geometric flow calculations (coupled polar and apolar components) rather than this model.
   .. _gamma:

gamma
=====

.. currentmodule:: apbs.input_file.calculate.nonpolar

.. note::  

   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`):  see :func:`Nonpolar.surface_tension` for more information.

This keyword specifies the surface tension coefficient for apolar solvation models.

.. code-block:: bash

   gamma { value }

where ``value`` is a floating point number designating the surface tension in units of kJ mol\ :superscript:`-1` Å\ :superscript:`-2`.
This term can be set to zero to eliminate the :abbr:`SASA (solvent-accessible surface area)` contributions to the apolar solvation calculations.

.. todo::

   Resolve unit confusion with geometric flow :ref:`gamma` keyword.
   https://github.com/Electrostatics/apbs/issues/490.. _press:

press
=====

.. currentmodule:: apbs.input_file.calculate.nonpolar

.. note::  

   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`):  see :func:`Nonpolar.pressue` for more information.

This term specifies the solvent pressure in kJ mol\ :superscript:`-1` Å\ :superscript:`-3`.
This coefficient multiplies the volume term of the apolar model and can be set to zero to eliminate volume contributions to the apolar solvation calculation.
The syntax is:

.. code-block:: bash

   press {value}

where ``value`` is the floating point value of the pressure coefficient in kJ mol\ :superscript:`-1` Å\ :superscript:`-3`.

.. todo::

   Resolve unit confusion with geometric flow ``press`` keyword and the apolar :ref:`press` keyword.
   Documented in https://github.com/Electrostatics/apbs/issues/499
.. _apolarsrfm:

srfm (apolar)
=============

.. currentmodule:: apbs.input_file.calculate.nonpolar

.. note::  

   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`):  see :func:`Nonpolar.surface_method` for more information.

This keyword specifies the model used to construct the solvent-related surface and volume.
The syntax is:

.. code-block:: bash

   srfm {flag}

where ``flag`` is a string that specifies the model used for surface and volume.
Acceptable values of flag include:

``sacc``
  Solvent-accessible (also called "probe-inflated") surface and volume.

  
.. _usemap:

usemap
======

.. note::

   .. currentmodule::  apbs.input_file.calculate.generic

   Some instances of this keyword have been ported to the *new APBS syntax* (see :ref:`new_input_format`); see :class:`UseMap`

  .. todo::  Port other uses to new syntax.

Specify pre-calculated coefficient maps to be used in the Poisson-Boltzmann calculation.
These must have been input via an earlier READ statement (see :ref:`read_old_input`).

The syntax for this command is:

.. code-block:: bash
   
   usemap {type} {id}

where the mandatory keywords are:

``type``
  A string that specifies the type of pre-calculated map to be read in:

  ``diel``
    Dielectric function map (as read by :ref:`read_old_input` ``diel``); this causes the :ref:`pdie`, :ref:`sdie`, :ref:`srad`, :ref:`swin`, and :ref:`elecsrfm` parameters and the radii of the biomolecular atoms to be ignored when computing dielectric maps for the Poisson-Boltzmann equation.
    Note that the :ref:`pdie` and :ref:`sdie` values are still used for some boundary condition calculations as specified by :ref:`bcfl`.
  ``kappa``
    Mobile ion-accessibility function map (as read by :ref:`read_old_input` ``kappa``); this causes the :ref:`swin` and :ref:`elecsrfm` parameters and the radii of the biomolecular atoms to be ignored when computing mobile ion values for the Poisson-Boltzmann equation.
    The :ref:`ion` parameter is not ignored and will still be used.
  ``charge``
    Charge distribution map (as read by :ref:`read_old_input` ``charge``); this causes the :ref:`chgm` parameter and the charges of the biomolecular atoms to be ignored when assembling the fixed charge distribution for the Poisson-Boltzmann equation.
  ``pot``
    Potential map (as read by :ref:`read_old_input` ``pot``); this option requires setting :ref:`bcfl` to ``map``.

``id``
  As described in the READ command documentation (see :ref:`read_old_input`), this integer ID specifies the particular map read in with READ.
  These IDs are assigned sequentially, starting from 1, and incremented independently for each map type read by APBS.
  In other words, a calculation that uses two PQR files, one parameter file, three charge maps, and four dielectric maps would have PQR files with IDs 1-2, a parameter file with ID 1, charge maps with IDs 1-3, and dielectric maps with IDs 1-4.

.. _ntraj:

ntraj
=====

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

Specify the number of Brownian Dynamic trajectories desired for the PB-(S)AM run.

.. code-block:: bash
   
   ntraj {traj}

where ``traj`` is an integer of the number of desired trajectories.
.. _tree_order:

tree_order
==========

.. currentmodule::  apbs.input_file.calculate.boundary_element

.. note::

   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`); see :func:`TABIParameters.tree_order`.

TABI-PB parameter that specifies the order of the treecode multipole expansion.
The syntax is:

.. code-block:: bash

   tree_order {order}

where ``order`` is an integer that indicates the Taylor expansion order.
Users can adjust the order for different accuracy. 
A typical choice for this parameter is 3.
.. _grid2d:

grid2d
======

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

Specify the filename and location of a 2D cross sectional potential to be written to.

.. code-block:: bash

   grid2d {filename} {axis} {axis_value}

``filename``
  String for the name of the 2D grid to be printed out

``axis``
  String of either x, y, or z, for which cartesian axis the grid will be computed along

``axis_value``
  A floating point number of the position along ``axis`` that will be used.

.. note::

   Multiple 2D files can be printed out with 1 PB-AM run. Just specify them with more grid2d flags.

.. todo::
   
   The PB-(S)AM ``grid2d`` keyword should not exist; please replace it ASAP with the :ref:`old_write` command.
   Documented in https://github.com/Electrostatics/apbs/issues/493

.. _3dmap:

3dmap
=====

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

Specify the name of the file into which the potential surface on the coarse-grain molecule surface will be printed.

.. code-block:: bash
   
   3dmap {filename}

where ``filename`` is a string for the name of the file where a 3D grid will be printed out.

.. todo::
   
   The PB-(S)AM ``3dmap`` keyword should not exist; please replace it ASAP with the :ref:`old_write` command.
   Documented this todo as https://github.com/Electrostatics/apbs/issues/482
.. _targetNum:

targetNum
=========

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

Specify the target number of vertices in the initial finite element mesh for :ref:`femanual` calculations.
The syntax is:

.. code-block:: bash

   targetNum { num }

where ``num`` is an integer denoting the target number of vertices in initial mesh.
Initial refinement will continue until this number is reached or the the longest edge of every simplex is below :ref:`targetRes`.
.. _old_write:

write
=====

.. currentmodule::  apbs.input_file.calculate

.. note::
   
   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`); see :class:`generic.WriteMap` and specific functions for different calculations:

   * Finite differences; see :func:`finite_difference.write_atom_potentials`.
   * Finite elements; see :func:`finite_elements.write_atom_potentials`.

This controls the output of scalar data calculated during the Poisson-Boltzmann run.
This keyword can be repeated several times to provide various types of data output from APBS.
The syntax is:

.. code-block:: bash

   write {type} {format} {stem}

``type``
  A string indicating what type of data to output:

  ``charge``
    Write out the biomolecular charge distribution in units of e\ :sub:`c` (electron charge) per Å\ :sup:`3` (multigrid only).
  ``pot``
    Write out the electrostatic potential over the entire problem domain in units of k\ :sub:`b` T e\ :sub:`c`\ :sup:`-1` (multigrid and finite element), where

    k\ :sub:`b`
      Boltzmann's constant:  1.3806504 × 10\ :sup:`−23` J K\ :sup:`-1`

    T
      The temperature of your calculation in K

    e\ :sub:`c`
      is the charge of an electron:  1.60217646 × 10\ :sup:`-19` C

    As an example, if you ran your calculation at 300 K, then the potential would be written out as multiples of
    k\ :sub:`b` T e\ :sub:`c`\ :sup:`-1` = (1.3806504 × 10\ :sup:`−23` J K\ :sup:`-1`) × (300 K) × (1.60217646 × 10\ :sup:`-19` C)\ :sup:`-1` = (4.1419512 × 10\ :sup:`-21` J) × (6.241509752 × 10\ :sup:`18` C\ :sup:`-1`) = 25.85202 mV

  ``atompot``
    Write out the electrostatic potential at each atom location in units of k\ :sub:`b` T e\ :sub:`c`\ :sup:`-1` (multigrid and finite element).
  ``smol``
    Write out the solvent accessibility defined by the molecular surface definition (see :ref:`elecsrfm` ``smol``).
    Values are unitless and range from 0 (inaccessible) to 1 (accessible). (multigrid and finite element).
  ``sspl``
    Write out the spline-based solvent accessibility (see :ref:`elecsrfm` ``spl2``).
    Values are unitless and range from 0 (inaccessible) to 1 (accessible) (multigrid and finite element)
  ``vdw``
    Write out the van der Waals-based solvent accessibility (see :ref:`elecsrfm` ``smol`` with :ref:`srad` 0.0).
    Values are unitless and range from 0 (inaccessible) to 1 (accessible). (multigrid and finite element)
  ``ivdw``
    Write out the inflated van der Waals-based ion accessibility (see :ref:`elecsrfm` ``smol``).
    Values are unitless and range from 0 (inaccessible) to 1 (accessible). (multigrid and finite element)
  ``lap``
    Write out the Laplacian of the potential :math:`\nabla^2 \phi` in units of k\ :sub:`B` T e\ :sub:`c`\ :sup:`-1` Å\ :sup:`-2`  (multigrid only).
  ``edens``
    Write out the "energy density" :math:`-\nabla \cdot \epsilon \nabla \phi` in units of k\ :sub:`B` T e\ :sub:`c`\ :sup:`-1` Å\ :sup:`-2`  (multigrid only).
  ``ndens``
    Write out the total mobile ion number density for all ion species in units of M (multigrid only).
    The output is calculated according to the formula (for nonlinear PB calculations):  :math:`\rho(x) = \sum_i^N {\bar{\rho}_i e^{-q_i\phi(x) - V_i (x)}}`, where *N* is the number of ion species, :math:`\bar{\rho}_i` is the bulk density of ion species *i*, :math:`q_i` is the charge of ion species *i*, :math:`\phi(x)` is the electrostatic potential, and :math:`V_i` is the solute-ion interaction potential for species *i*.
  ``qdens``
    Write out the total mobile ion charge density for all ion species in units of e\ :sub:`c` M (multigrid only).
    The output is calculated according to the formula (for nonlinear PB calculations):  :math:`\rho(x) = \sum_i^N {\bar{\rho}_i q_i e^{-q_i\phi(x) - V_i (x)}}`, where *N* is the number of ion species, :math:`\bar{\rho}_i` is the bulk density of ion species *i*, :math:`q_i` is the charge of ion species *i*, :math:`\phi(x)` is the electrostatic potential, and :math:`V_i` is the solute-ion interaction potential for species *i*.
  ``dielx`` or ``diely`` or ``dielz``
    Write out the dielectric map shifted by 1/2 grid spacing in the {x, y, z}-direction (see :ref:`read_old_input` ``diel``).
    The values are unitless (multigrid only).

``format``
  A string that specifies the format for writing out the data:

  ``dx``
    Write out data in :doc:`/formats/opendx`.
    This is the preferred format for APBS I/O. (multigrid and finite element).

  ``avs``
    Write out data in AVS UCD format. (finite element only).

  ``uhbd``
    Write out data in :doc:`/formats/uhbd`. (multigrid only).

  ``gz``
    Write out :doc:`/formats/opendx` in gzipped (zlib) compatible format.
    Appends .dx.gz to the filename.

  ``flat``
    Write out data as a plain text file. (multigrid and finite element).

``stem``
  A string that specifies the path for the output; files are written to :file:`stem.{XYZ}`, where ``XYZ`` is determined by the file format (and processor rank for parallel calculations).
  If the pathname contains spaces, then it must be surrounded by double quotes.
.. _akeyPRE:

akeyPRE
=======

.. currentmodule:: apbs.input_file.calculate.finite_element

.. note::  This command has been ported to :func:`FiniteElement.a_priori_refinement`.

Specifies how the initial finite element mesh should be constructed (from refinement of a very coarse 8-tetrahedron mesh prior to the solve-estimate-refine iteration in :ref:`femanual` finite element calculations.
The syntax is:

.. code-block:: bash

   akeyPRE {key}

where ``key`` is a text string that specifies the method used to guide initial refinement and takes one of the values:

``unif``
  Uniform refinement
``geom``
  Geometry-based refinement at molecular surfaces and charges

.. _randorient:

randorient
==========

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

Flag to indicate that the molecules should have a random orientation in subsequent PB-(S)AM calculations.
.. _writemat:

writemat
========

.. note::

   This command is deprecated and will not ported to the *new APBS syntax* (see :ref:`new_input_format`).

This controls the output of the mathematical operators in the Poisson-Boltzmann equation as matrices in Harwell-Boeing matrix format (multigrid only).
The syntax is:

.. code-block:: bash
   
   writemat {type} {stem}

where

``type``
  A string that indicates what type of operator to output.

  ``poisson``
    Write out the Poisson operator :math:`-\nabla \cdot \epsilon \nabla`.

``stem``
  A string that specifies the path for the output.
.. _glen:

glen
====

.. currentmodule:: apbs.input_file.calculate.finite_difference

.. note::  

   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`); see :class:`GridDimensions` for more information.

   
Specify the mesh domain lengths for multigrid :ref:`mgmanual` calculations.
These lengths may be different in each direction.
The syntax is:

.. code-block:: bash
   
   glen {xlen ylen zlen}

where ```xlen ylen zlen`` are the (floating point) grid lengths in the x-, y-, and z-directions (respectively) in Å.
.. _gridpts:

gridpts
=======

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

Specify the number of gridpoints in each cartesian dimension.

.. code-block:: bash
   
   gridpts {pts}

where ``pts`` is a integer number indicating the number of grid points.

.. todo::
   
   The PB-(S)AM ``gridpts`` keyword should not exist; it's duplicative of the existing :ref:`dime` keyword!
   Documented in https://github.com/Electrostatics/apbs/issues/494

.. _async:

async
=====

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

An optional keyword to perform an asynchronous parallel focusing Poisson-Boltzmann equation.
The syntax is

.. code-block:: bash

   async {rank}

where ``rank`` is the integer ID of the particular processor to masquerade as.
Processor IDs range from *0* to *N-1*, where *N* is the total number of processors in the run (see :ref:`pdime`).
Processor IDs are related to their position in the overall grid by :math:`p = nx ny k + nx j + i`  where :math:`nx` is the number of processors in the x-direction, :math:`ny` is the number of processors in the y-direction, :math:`nz` is the number of processors in the z-direction, :math:`i` is the index of the processor in the x-direction, :math:`j` is the index of the processor in the y-direction, :math:`k` is the index of the processor in the z-direction, and :math:`p` is the overall rank of the processor.

.. _npbe:

npbe
====

.. currentmodule:: apbs.input_file.calculate

.. note::  

   Some aspects of this command have been moved to the *new APBS syntax* (see :ref:`new_input_format`): 

   * Finite difference:  see :func:`finite_difference.FiniteDifference.equation` for more information.

   * Finite element:  see :func:`finite_element.FiniteElement.equation` for more information.

.. todo::  port for other types of calculations.

Specifies that the nonlinear (full) Poisson-Boltzmann equation should be solved.

.. note::

   The options :ref:`lpbe`, :ref:`npbe`, :ref:`lrpbe`, :ref:`nrpbe` are mutually exclusive.
.. _mgmanual:

mg-manual
=========

.. currentmodule:: apbs.input_file.calculate.finite_difference

.. note::

   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`); see :class:`Manual` and :func:`FiniteDifference.boundary_condition`.

Manually-configured finite differnece multigrid Poisson-Boltzmann calculations.

This is a standard single-point multigrid PBE calculation without focusing or additional refinement.
The ``mg-manual`` calculation offers the most control of parameters to the user.
Several of these calculations can be strung together to perform focusing calculations by judicious choice of the :ref:`bcfl` flag; however, the setup of the focusing is not automated as it is in :ref:`mgauto` and :ref:`mgpara` calculations and therefore this command should primarily be used by more experienced users.

.. toctree::
   :maxdepth: 2
   :caption: ELEC mg-manual keywords:

   bcfl
   ../generic/calcenergy
   ../generic/calcforce
   chgm
   dime
   etol
   gcent
   glen
   ../generic/grid
   ion
   lpbe
   lrpbe
   ../generic/mol
   nlev
   npbe
   pdie
   ../generic/sdens
   sdie
   ../generic/srad
   srfm
   ../generic/swin
   ../generic/temp
   usemap
   write
   writemat
.. _maxsolve:

maxsolve
========

.. currentmodule::  apbs.input_file.calculate.finite_element

.. note::  This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`); see :func:`FiniteElement.maximum_refinement_iterations`. 

Specify the number of times to perform the solve-estimate-refine iteration of the finite element solver (:ref:`femanual`).
The syntax is:

.. code-block:: bash
   
   maxsolve { num }

where `num` is an integer indicating the desired maximum number of iterations.
.. _runname:

runname
=======

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

Specify the output name for the PB-(S)AM calculation.

..code-block:: bash

   runname {name}

where ``name`` is a string indicating the prefix for all PB-(S)AM output files.
.. _gcent:

gcent
======

.. currentmodule::  apbs.input_file.calculate.finite_difference

.. note::  This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`); see :class:`GridCenter` for more information.

Specify the center of the grid based on a molecule's center or absolute coordinates :ref:`mgmanual` multigrid calculations.
The syntax is:

.. code-block:: bash
   
   fgcent { mol id | xcent ycent zcent }

where a user can specify **either**

``mol {id}``
  Center the grid on molecule with integer ID id; as assigned in the READ section (see :ref:`read_old_input`) of the input file.
  Molecule IDs are assigned in the order they are read, starting at 1.

**or** the user can specify

``xcent ycent zcent``
  Center the grids on the coordinates (floating point numbers in Å) at which the grid is centered.
  Based on the input molecule PDB coordinate frame.

.. _nlev:

nlev
====

.. note::

   ..currentmodule:: apbs.input_file.calculate.finite_difference

   This command has been eliminated in the *new APBS syntax* (see :ref:`new_input_format`); see :class:`GridDimensions` for more information.


Specify the depth of the multilevel hierarchy used in the :ref:`mgmanual` multigrid solver.
See :ref:`dime` for a discussion of how nlev relates to grid dimensions.
The syntax is:

.. code-block:: bash
   
   nlev {lev}

where ``lev`` is an integer indicating the desired depth of the multigrid hierarchy.

.. _geoflowauto:

geoflow-auto
============

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

To increase the accuracy of our implicit solvent modeling, we have implemented a differential geometry based geometric flow
solvation model `(Thomas, 2013) <https://www.ncbi.nlm.nih.gov/pubmed/23212974>`_.
In this model, polar and nonpolar solvation free energies are coupled and the solvent-solute boundary is determined in a self-consistent manner.
Relevant references are provided in :doc:`Recommended reading </reading>`.
This section provides a brief overview of the method.

The solutions for the electrostatic potential :math:`\phi` and the characteristic function :math:`S` (related to the solvent density) are obtained by minimizing a free energy functional that includes both polar and nonpolar solvation energy terms.
Minimization of the functional with respect to :math:`\phi` gives the Poisson-Boltzmann equation with a dielectric coefficient :math:`\epsilon` has the solute value :math:`\epsilon_m` where :math:`S = 1` and the solvent value :math:`\epsilon_s` where :math:`S = 0`.
Minimization of the free energy functional with respect to :math:`S` gives

.. math::

   -\nabla\cdot\left(\gamma\frac{\nabla S}{\parallel\nabla S\parallel}\right)+p-\rho_0U^{att}+\rho_m\phi - \frac{1}{2}\epsilon_m\mid\nabla\phi\mid^2+\frac{1}{2}\epsilon_s\mid\nabla\phi\mid^2=0 

where :math:`\gamma` is the microscopic surface tension, :math:`p` is the hydrostatic pressure, and :math:`U^{att}` is the attractive portion of the van der Waals dispersion interaction between the solute and the solvent.

Keywords for this calculation type include:

.. toctree::
   :maxdepth: 2
   :caption: ELEC geoflow-auto keywords:

   bcfl
   ../generic/bconc
   etol
   gamma-geoflow
   lpbe
   ../generic/mol
   pdie
   press-geoflow
   sdie
   vdwdisp

.. warning::

   Although the ``ion`` and ``lpbe`` keywords will be accepted in the geoflow-auto calculation, the treatment of salt is not currently implemented in APBS geometric flow.

.. todo::
   
   Add LPBE/NPBE support to geometric flow or remove the ``ion`` and ``lpbe`` keywords.
   Documented in https://github.com/Electrostatics/apbs/issues/491

.. todo::
   
   If there's only one mode, then we can change the keyword from ``geoflow-auto`` to just ``geoflow``.
   Documented in https://github.com/Electrostatics/apbs/issues/492

.. _salt:

salt
====

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

Specify the monovalent salt concentration of the system, in molar. This is usually a value between 0.00 to 0.15.

.. code-block:: bash
   
   salt {saltConc}

where ``saltConc`` is the floating point value of the monovalent salt concentration in molar.

.. todo::

   The PB-(S)AM ``salt`` keyword should be eradicated and replaced with the :ref:`ion` keyword.
   Documented in https://github.com/Electrostatics/apbs/issues/501
   
.. _xyz:

xyz
===

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

For each molecule in the system and for each trajectory, specify a xyz file for the starting position of that molecule.
The syntax is:

.. code-block:: bash
   
   xyz {molecule_id} {filename}

``molecule_id``
  An integer (starting at 1) of the molecule index from the READ  section

``filename``
  The name of the file for the xyz coordinates of the molecule center for a given trajectory.
  The trajectories for a given molecule should be ordered sequentially in the ELEC section.

.. todo::
   
   It would be nice to incorporate the ``xyz`` functionality into the :ref:`read_old_input` block.
   Documented in https://github.com/Electrostatics/apbs/issues/505
.. _vdwdisp:

vdwdisp
=======

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

Specify whether the attractive van der Waals contribution to the geometric flow potential is on or off.

.. code-block:: bash

    vdwdisp { flag }

where ``flag`` is 0 (vdw off) or 1 (vdw on).
.. _imat:

imat
====

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

This keyword can be used to load in the surface integral matrices previously generated by PB-SAM named as :file:`mol{m}sph{s}.bin` for molecule ID *s* and coarse-grained sphere *s* (see :ref:`pbamauto` for more information).
The syntax is:

.. code-block:: bash
   
   imat {prefix}

where ``prefix`` is the filename prefix :file:`mol{m}sph`.
The :file:`{s}.bin` will be appended during the program run.

.. todo::

   It would be better to generalize the :ref:`read_old_input` section of the input file rather than use the ``imat`` command.
   This command also needs to be cleaned up -- it's too fragile.
   Documented in https://github.com/Electrostatics/apbs/issues/495.. _pdie:

pdie
====

.. currentmodule::  apbs.input_file.calculate

.. note:: 

   Some versions of this command have been ported to the *new APBS syntax* (see :ref:`new_input_format`):

   * Boundary element Poisson-Boltzmann calculations: see :func:`boundary_element.BoundaryElement.solute_dielectric` for more information.
   * Finite difference Poisson-Boltzmann calculations: see :func:`finite_difference.FiniteDifference.solute_dielectric` for more information.
   * Finite element Poisson-Boltzmann calculations: see :func:`finite_element.FiniteElement.solute_dielectric` for more information.

.. todo:: port for other calculation types

Specify the dielectric constant of the solute molecule.
The syntax is:

.. code-block:: bash

   pdie {diel}

where ``diel`` is the floating point value of the unitless biomolecular dielectric constant.
This is usually a value between 2 to 20, where lower values consider only electronic polarization and higher values consider additional polarization due to intramolecular motion.
The dielectric value must be :math:`\ge 1`.
.. _pbamauto:

pbam-auto
=========

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

PB-AM is an analytical solution to the linearized Poisson-Boltzmann equation for multiple spherical objects of arbitrary charge distribution in an ionic solution.
More details on the method are available in `Lotan, Head-Gordon (2006) <http://pubs.acs.org/doi/full/10.1021/ct050263p>`_.
The physical calculations are uses to perform various actions on a system of molecules such as calculation of energies, forces, torques, electrostatic potentials, and Brownian dynamics schemes.
This fast method coarse-grains all molecules of the system into single spheres large enough to contain all molecule atoms.

.. todo::

   If there's only one mode to PBAM, let's call it ``pbam`` instead of ``pbam-auto``.
   Documented in https://github.com/Electrostatics/apbs/issues/498

The current implementation of PB-AM in APBS includes:

* Calculation of energies, forces and torques
* Calculation of electrostatic potentials
* Brownian dynamics simulations

Keywords for this calculation type include:

.. toctree::
   :maxdepth: 2
   :caption: ELEC pbam-auto keywords:

   3dmap
   diff
   dx
   grid2d
   gridpts
   ../generic/mol
   ntraj
   pbc
   pdie
   randorient
   runname
   runtype
   salt
   sdie
   ../generic/temp
   term
   termcombine
   units
   xyz

======================
Background information
======================

PB-AM is an analytical solution to the linearized Poisson-Boltzmann equation for multiple spherical objects of arbitrary charge distribution in an ionic solution.
The solution can be reduced to a simple system of equations as follows:

.. math::

   A = \Gamma \cdot (\Delta \cdot T \cdot A + E) 

Where :math:`A^{(i)}` represents the effective multipole expansion of the charge distributions of molecule :math:`i`.
:math:`E^{(i)}` is the free charge distribution of molecule :math:`i`.
:math:`\Gamma` is a dielectric boundary-crossing operator, :math:`\Delta` is a cavity polarization operator, :math:`T` an operator that transforms the multipole expansion to a local coordinate frame.
:math:`A^{(i)}` is solved for through an iterative SCF method.


From the above formulation, computation of the interaction energy :math:`\Omega^{(i)}` for molecule :math:`i`, is given as follows:

.. math::

   \Omega^{(i)}=\frac{1}{\epsilon_s} \left \langle \sum_{j \ne i}^N  T \cdot A^{(j)} ,  A^{(i)} \right \rangle

where :math:`\langle M, N \rangle` denotes the inner product.
Forces can be obtained from

.. math::

   \textbf{F}^{(i)} = \nabla_i \Omega^{(i)}=\frac{1}{\epsilon_s} \left[ \langle \nabla_i \,T \cdot A^{(i)} ,  A^{(i)} \rangle +  \langle T \cdot A^{(i)} ,   \nabla_i \, A^{(i)} \rangle \right]

.. _chgm:

chgm
====

.. currentmodule:: apbs.input_file.calculate

.. note::  

   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`):
   
   * For finite difference, see :func:`finite_difference.FiniteDifference.charge_discretization` for more information.
   * For finite element, see :func:`finite_element.FiniteElement.charge_discretization` for more information.

Specify the method by which the biomolecular point charges (i.e., Dirac delta functions) by which charges are mapped to the grid for a multigrid (:ref:`mgmanual`, :ref:`mgauto`, :ref:`mgpara`) Poisson-Boltzmann calculation.
As we are attempting to model delta functions, the support (domain) of these discretized charge distributions is always strongly dependent on the grid spacing.
The syntax is:

.. code-block:: bash

   chgm {flag}

``flag`` is a text string that specifies the type of discretization:

``spl0``
  Traditional trilinear interpolation (linear splines).
  The charge is mapped onto the nearest-neighbor grid points.
  Resulting potentials are very sensitive to grid spacing, length, and position.
``spl2``
  Cubic B-spline discretization.
  The charge is mapped onto the nearest- and next-nearest-neighbor grid points.
  Resulting potentials are somewhat less sensitive (than ``spl0``) to grid spacing, length, and position.
``spl4``
  Quintic B-spline discretization.
  Similar to ``spl2``, except the charge/multipole is additionally mapped to include next-next-nearest neighbors (125 grid points receive charge density).

.. _pbsamauto:

pbsam-auto
==========

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

PB-SAM is a semi-analytical solution to the linearized Poisson-Boltzmann equation for multiple molecules of arbitrary charge distribution in an ionic solution.
The solution is an extension of the :ref:`analytical method <pbamauto>`, leveraging fast-multipole methods as well as boundary elements.
Each molecule is coarse-grained as a system of overlapping spheres, whose surface charges are represented by multipole expansions.
For details on the method, please see `Yap, Head-Gordon (2010) <http://pubs.acs.org/doi/abs/10.1021/ct100145f>`_ and `Yap, Head-Gordon (2013) <http://pubs.acs.org/doi/abs/10.1021/ct400048q>`_.

.. todo::

   If there's only one mode to PBAM, let's call it ``pbsam`` instead of ``pbsam-auto``.

The current implementation of PB-SAM in APBS includes:

* Calculation of energies, forces and torques
* Calculation of electrostatic potentials
* Brownian dynamics simulations

Keywords for this calculation type include:

.. toctree::
   :maxdepth: 2
   :caption: ELEC pbsam-auto keywords:

   3dmap
   diff
   dx
   exp
   grid2d
   imat
   ntraj
   pbc
   pdie
   randorient
   runname
   runtype
   salt
   sdie
   surf
   ../generic/temp
   term
   termcombine
   tolsp
   units
   xyz

======================
Background information
======================

PB-SAM is a semi-analytical solution to the linearized Poisson-Boltzmann equation for multiple molecules of arbitrary charge distribution in an ionic solution.
The solution is an extension of the analytical method, leveraging Fast-Multipole methods as well as boundary elements.
Each molecule is coarse-grained as a system of overlapping spheres, whose surface charges are represented by the multipole expansions :math:`H^{(i)}` and :math:`F^{(i)}`.
To solve for the potential, the following interactions are considered:

* Intra-molecular interactions between overlapping spheres are treated numerically
* Intra-molecular interactions between non-overlapping spheres are treated analytically
* Inter-molecular interactions between spheres on different molecules

With these interactions, the multipole expansions are solved with an iterative SCF method, briefly given as

.. math::

   H^{(i,k)} &= I_{E}^{(i,k)} \cdot \left ( H^{(i,k)} + F^{(i,k)} + T \cdot H^{(j,l)} \right ) \\
   F^{(i,k)} &= I_{E}^{(i,k)} \cdot \left ( H^{(i,k)} + F^{(i,k)} + T \cdot F^{(j,l)} \right )

Where :math:`H^{(i)}` and :math`F^{(i)}` are multipole expansions, :math:`I_{E}^{(i,k)}` is the exposed surface integral matrix for sphere :math:`k` of molecule :math:`i`, and :math:`T` is an operator that transforms the multipole expansion to a local coordinate frame.

From the above formulation, computation of the interaction energy :math:`\Omega^{(i)}` for molecule :math:`i`, is given as a sum of all the interactions of spheres :math:`k` within it with all external spheres (in a simplified form) as follows:

.. math::

   \Omega^{(i)} = \frac{1}{\epsilon_s} \left \langle \sum_{k \, in\, i} \sum_{j \ne i}^N \sum_{l\, in \, j}  T \cdot H^{(j,l)} ,  H^{(i,k)} \right \rangle

where :math:`\langle  M, N \rangle` denotes the inner product.

When energy is computed, forces follow as:

.. math::

   \textbf{F}^{(i)} = \nabla_i \Omega^{(i)}=\frac{1}{\epsilon_s} [ \langle \nabla_i \,T \cdot H^{(j,l)} ,  H^{(i,k)} \rangle +  \langle T \cdot H^{(j,l)},   \nabla_i \, H^{(i,k)} \rangle

The method to calculate the torque is discussed in `Yap, Head-Gordon (2013) <http://pubs.acs.org/doi/abs/10.1021/ct400048q>`_.

============
PB-SAM files
============

-------------------
Vertex/surface file
-------------------

As part of the coarse-graining process a definition of the molecular surface is necessary.

-----------------------
Coarse-grained PQR file
-----------------------

The coarse-graining process will produce a new PQR file :file:`mol{#}_cg.pqr` that contains the original PQR concatenated with coarse-graining spherical centers.
The number `#` refers to the order the file was read during the :doc:`../read` statements.

---------------------------
IMAT: surface integral file
---------------------------

The surface integrals are computed for the boundary element part of PB-SAM.
Their calculation can be quite time-consuming, so the first time they are computed for a system, they are saved to the working directory with the name :file:`mol{m}sph{s}.bin``.
The *m* in :file:`mol{m}sph{s}.bin`` is the ordered ID of the molecule from the PQR section.
The *s* in :file:`mol{m}sph{s}.bin`` refers to coarse-grained sphere *s* of the molecule.

-------------------------
Multipole expansion files
-------------------------

Much like the IMAT files, the expansion files are files generated from self-polarization that are useful and time-saving methods for running a system of full-mutual polarziation on many molecules.
If no expansion path is provided, the program will calculate self-polarization for each type of molecule in the system and save files of the form :file:`mol{m}{H,F}.{s}.exp`, where *m* is the molecule ID, *H* and *F* refer to the respective expansion (see above), and *s* is the coarse-grained sphere number.
.. _ion:

ion
===

.. note::

   .. currentmodule:: apbs.input_file.calculate

   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`); see:
   
   * General information about this object:  :class:`generic.MobileIons` for more information.
   * Boundary element implementation:  :func:`boundary_element.BoundaryElement.ions`.
   * Finite difference implementation:  :func:`finite_difference.FiniteDifference.ions`.
   * Finite element implementation:  :func:`finite_element.FiniteElement.ions`.

Specify the bulk concentrations of mobile ion species present in the system.
This command can be repeated as necessary to specify multiple types of ions; however, only the largest ionic radius is used to determine the ion-accessibility function.
The total bulk system of ions must be electroneutral which means the charge densities/concentrations of positive and negative ions must be equal.
The syntax is:

.. code-block:: bash

   ion charge {charge} conc {conc} radius {radius}

where

``charge``
  Mobile ion species charge (floating point number in ec)
``conc``
  Mobile ion species concentration (floating point number in M)
``radius``
  Mobile ion species radius (floating point number in Å)

.. _domainLength:

domainLength
============

.. currentmodule:: apbs.input_file.calculate.finite_element

.. note::  This command has been ported to the *new APBS syntax* (see :func:`FiniteElement.domain_length`).

Specify the rectangular finite element mesh domain lengths for :ref:`femanual` finite element calculations.
This length may be different in each direction.
If the :ref:`usemesh` keyword is included, then this command is ignored.
The syntax is:

.. code-block:: bash

   domainLength {xlen ylen zlen}

where the parameters ``xlen ylen zlen`` are floating point numbers that specify the mesh lengths in the x-, y-, and z-directions (respectively) in units of Å.

.. _pbc:

pbc
===

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

This keyword is used to indicate if 3D periodic boundary conditions (PBCs) will be used in a PB-(S)AM calculation.
If used, a box length must also be specified, in Ångstroms.

.. code-block:: bash
   
   pbc {boxlength}

where ``boxlength`` is the floating point value of the box length in Ångstroms.

.. note::

   The box is centered at the origin (0, 0, 0).
   The code assumes a minimum image convention, so it only includes the closest image of the neighboring molecules.
   For this convention to always be preserved, the periodic box is assumed to be large enough such that the electrostatic forces are sufficiently attenuated beyond one boxlength.
   Generally, the program assumes a mutual polarization cutoff of 100 Å for the mutual polarization, so if the boxlength is shorter, the cutoff will be reduced to boxlength/2.
.. _lrpbe:

lrpbe
=====

.. currentmodule:: apbs.input_file.calculate

.. note::  

   Some aspects of this command have been moved to the *new APBS syntax* (see :ref:`new_input_format`): 

   * Finite difference:  see :func:`finite_difference.FiniteDifference.equation` for more information.

   * Finite element:  see :func:`finite_element.FiniteElement.equation` for more information.

.. todo::  port for other types of calculations.


Specifies that the linear form of the regularized Poisson-Boltzmann equation (RPBE) should be solved.
The regularized PBE equation replaces the point charge distribution with the corresponding Green's function.
As a result of this replacement, the solution corresponds to the reaction field instead of the total potential; the total potential can be recovered by adding the appropriate Coulombic terms to the solution.
Likewise, this equation immediately yields the solvation energy without the need for reference calculations.

.. note::

   The options :ref:`lpbe`, :ref:`npbe`, :ref:`lrpbe`, :ref:`nrpbe` are mutually exclusive... _mgpara:

mg-para
=======

.. currentmodule:: apbs.input_file.calculate.finite_difference

.. note::

   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`); see :func:`Focus.parallel` and :func:`FiniteDifference.calculation_type`.


Automatically-configured parallel focusing multigrid Poisson-Boltzmann calculations.

This calculation closely resembles :ref:`mgauto` in syntax.
However, it is designed to perform electrostatics calculations on systems in a parallel focusing fashion.

.. toctree::
   :maxdepth: 2
   :caption: ELEC mg-para keywords:

   async
   bcfl
   ../generic/calcenergy
   ../generic/calcforce
   cgcent
   cglen
   chgm
   dime
   etol
   fgcent
   fglen
   ion
   lpbe
   lrpbe
   ../generic/mol
   npbe
   ofrac
   pdie
   pdime
   ../generic/sdens
   sdie
   ../generic/srad
   srfm
   ../generic/swin
   ../generic/temp
   usemap
   write
   writemat
.. _targetRes: 

targetRes
=========

.. currentmodule::  apbs.input_file.calculate.finite_element

.. note::
   
   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`); see :func:`FiniteElement.initial_mesh_resolution`.

Specify the target resolution of the simplices in a finite element mesh (:ref:`femanual`).
The syntax is:

.. code-block:: bash
   
   targetRes { res }

where ``res`` is a floating point number denoting the target resolution for longest edges of simplices in mesh (in Å).
Refinement will continue until the longest edge of every simplex is below this value or the number of vertices reaches :ref:`targetNum`.
.. _ofrac:

ofrac
=====

.. currentmodule::  apbs.input_file.calculate.finite_difference

.. note::  This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`); see :func:`ParallelFocus.overlap_fraction`.

Specify the amount of overlap to include between the individual processors meshes in a parallel focusing calculation (:ref:`mgpara`).
The syntax is:

.. code-block:: bash

   ofrac {frac}

where ``frac`` is a floating point value between 0.0 and 1.0 denoting the amount of overlap between processors.
Empirical evidence suggests that an value of 0.1 is sufficient to generate stable energies.
However, this value may not be sufficient to generate stable forces and/or good quality isocontours.
For example, the following table illustrates the change in energies and visual artifacts in isocontours as a function of ``ofrac`` values for a small peptide (2PHK:B).

.. list-table:: Sensitivity of 2PHK:B solvation energy calculations to ofrac values.
   :widths: auto
   :header-rows: 1

   * - ``ofrac`` value
     - Energy (kJ/mol)
     - Visual artifact in isocontour?
   * - 0.05
     - 342.79
     - No
   * - 0.06
     - 342.00
     - No
   * - 0.07
     - 341.12
     - Yes
   * - 0.08
     - 341.14
     - Yes
   * - 0.09
     - 342.02
     - Yes
   * - 0.10
     - 340.84
     - Yes
   * - 0.11
     - 339.67
     - No
   * - 0.12
     - 341.10
     - No
   * - 0.13
     - 341.10
     - No
   * - 0.14
     - 341.32
     - No
   * - 0.15
     - 341.54
     - No

In general, larger ``ofrac`` values will reduce the parallel efficiency but will improve the accuracy.

For broad spatial support of the splines, every charge included in partition needs to be at least 1 grid space (:ref:`chgm` ``spl0``), 2 grid spaces (:ref:`chgm` ``spl2``), or 3 grid spaces (:ref:`chgm` ``spl4``) away from the partition boundary.

.. _termcombine:

termcombine
===========

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

Combine multiple PB-(S)AM Brownian dynamics trajectory termination conditions (see :doc:`term`) via logic operators

.. code-block:: bash
   
   termcombine {op}

where ``op`` is either the string ``or`` or ``and``.
If  ``and`` is selected, all listed termination conditions must be fulfilled before the simulation ends.
If ``or`` is selected, only one needs to be satisfied to complete the simulation.
.. _pygbe:

pygbe
============

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

.. todo::  @jbardhan should add references to papers and any other details about PYGBE

Keywords for this calculation type include:

.. toctree::
   :maxdepth: 2
   :caption: ELEC pygbe keywords:
.. _cglen:

cglen
=====

.. currentmodule:: apbs.input_file.calculate.finite_difference

.. note::  

   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`):  see :func:`Focus.coarse_grid_dimensions` for more information.

Specify the length of the coarse grid (in a focusing calculation) for an automatic multigrid (:ref:`mgauto`, :ref:`mgpara`) Poisson-Boltzmann calculation.
This may be different in each direction.

.. code-block:: bash

   cglen {xlen ylen zlen}

This is the starting mesh, so it should be large enough to complete enclose the biomolecule and ensure that the chosen boundary condition (see :ref:`bcfl`) is appropriate.

``xlen ylen zlen``
  Grid lengths (floating point numbers) in the x-, y-, and z-directions in Å.
.. _femanual:

fe-manual
=========

.. currentmodule:: apbs.input_file.calculate.finite_element

.. note::  This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`); see :class:`FiniteElement` for more information.

Manually-configured adaptive finite element Poisson-Boltzmann calculations.

This is a single-point PBE calculation performed by our adaptive finite element PBE solver.
It requires that APBS be linked to the Michael Holst group `FEtk finite element library <http://www.fetk.org>`_ during compilation.
The finite element solver uses a "solve-estimate-refine" cycle.
Specifically, starting from an initial mesh, it performs the following iteration:

#. solve the problem with the current mesh
#. estimate the error in the solution
#. adaptively refine the mesh to reduce the error

This iteration is repeated until a global error tolerance is reached.

Keywords for this calculation type include:

.. toctree::
   :maxdepth: 2
   :caption: ELEC fe-manual keywords:

   akeyPRE
   akeySOLVE
   bcfl
   ../generic/calcenergy
   ../generic/calcforce
   chgm
   domainLength
   ekey
   etol
   ion
   lpbe
   lrpbe
   maxsolve
   maxvert
   ../generic/mol
   npbe
   nrpbe
   pdie
   ../generic/sdens
   sdie
   ../generic/srad
   srfm
   ../generic/swin
   targetNum
   targetRes
   ../generic/temp
   usemesh
   write


.. note::

   The finite element methods are currently most useful for a select set of problems which can benefit from adaptive refinement of the solution.
   Furthermore, this implementation is experimental.
   In general, the sequential and parallel focusing multigrid methods offer the most efficient solution of the PBE for most systems.
.. _units_pbsam:

units
=====

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

Specify the units for energy/force/potential output in PB-(S)AM calculations:

.. code-block:: bash
   
   units {flag}

where ``flag`` specifies the unit system:

``kcalmol``
  kcal/mol

``jmol``
  J/mol

``kT``
  kT

Force units will be energy units/Angstrom and potential units will be energy units/electron.

.. todo::

   It would be great to use the same units everywhere in APBS.
   Documented in https://github.com/Electrostatics/apbs/issues/485
.. _maxvert:

maxvert
=======

.. currentmodule::  apbs.input_file.calculate.finite_element

.. note::  This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`); see :func:`FiniteElement.maximum_vertices`. 

Specify the maximum number of vertices to allow during solve-estimate-refine cycle of finite element solver (:ref:`femanual`).
This places a limit on the memory that can be used by the solver.
The syntax is:

.. code-block:: bash
   
   maxvert { num }

where ``num`` is an integer indicating the maximum number of vertices.
.. _mesh:

mesh
====

.. currentmodule::  apbs.input_file.calculate.boundary_element

.. note::

   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`); see :func:`Mesh.surface_method` and :func:`Mesh.software`.


TABI-PB parameter that specifies the meshing software used to generate surface mesh.
The syntax is:

.. code-block:: bash

   mesh {flag}

where ``flag`` is an integer indicating the meshing software to be used:

.. _NanoShaper: https://www.electrostaticszone.eu/downloads

0
  Formerly used for msms, no longer supported.
1
  SES implementation in NanoShaper_
2
  Skin surface implementation in NanoShaper_

Note that the executable NanoShaper_ must be included in your path to use them.
.. _surf:

surf
====

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

This keyword can be used to load in the MSMS vertex file for coarse-graining (see :ref:`pbsamauto`)
The syntax is:

.. code-block:: bash

   surf {prefix}

where ``prefix`` refers to the filename :file:{prefix}.vert`.

.. todo::
   
   The PB-SAM ``surf`` command is redundant with and should be replaced by the existing :ref:`usemesh` command.
   Documented in https://github.com/Electrostatics/apbs/issues/502
.. _cgcent:

cgcent
======

.. currentmodule:: apbs.input_file.calculate.finite_difference

.. note::  

   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`):  see :func:`Focus.coarse_grid_center` for more information.

This keyword controls electrostatic energy output from a Poisson-Boltzmann calculation
The syntax is:

.. code-block:: bash

   cgcent { mol id | xcent ycent zcent }

The arguments for this keyword are **either**

``mol id``
  Center the grid on molecule with integer ID ``id``; as assigned in the ``READ`` section with a ``READ mol`` command (see :ref:`read_old_input`)

**or**

``xcent ycent zcent``
  Center the grid on the (floating point) coordinates (in Å) at which the grid is centered.
  Based on the PDB coordinate frame.

.. _elec:

ELEC input file section
=======================

.. note::

   Some versions of this command have been ported to the *new APBS syntax* (see :ref:`new_input_format`):

   * :ref:`boundary_element`
   * :ref:`finite_difference`
   * :ref:`finite_element`

.. todo::  port other versions of command

The ELEC block of an APBS input file is used for polar solvation (electrostatics) calculations and has the following syntax:

.. code-block:: bash
   
   ELEC [ name {id} ]
        {type}
        {keywords...}
   END

The optional ``id`` variable is a simple string that allows ELEC statements to be named.
Since numerous ELEC blocks may appear in an APBS input file, it can be difficult to keep track of them all.
It is possible to assign an optional name (string) to each ELEC block to simplify the organizational process.

The ``type`` command defines the types of ELEC calculation to be performed and includes:

* Finite difference multigrid calculations with `PMG <http://www.fetk.org>`_.

  * :ref:`mgauto`
  * :ref:`mgpara`
  * :ref:`mgmanual`

* `Geometric flow solvation <https://www.ncbi.nlm.nih.gov/pubmed/23212974>`_ finite difference calculations

  * :ref:`geoflowauto`

* Boundary element method calculations with `TABI-PB <https://doi.org/10.1016/j.jcp.2013.03.056>`_.

  * :ref:`tabi`

* Analytic and semi-analytic Poisson-Boltzmann approximations

  * :ref:`pbamauto`
  * :ref:`pbsamauto`

* Finite element calculations with `FEtk <http://www.fetk.org>`_.

  * :ref:`femanual`

* No-op modes for generating coefficient maps

  * :ref:`mgdummy`

Finally, the ``keywords`` are calculation-specific commands that customize the particular type of calculation.
This section is the main component for polar solvation calculations in APBS runs.
There may be several ELEC sections, operating on different molecules or using different parameters for multiple runs on the same molecule.
The order of the ELEC statement can matter since certain types of boundary conditions (:ref:`bcfl`) can require information about previous calculations.


.. toctree::
   :maxdepth: 1
   :caption: Calculation type keywords

   tabi
   fe-manual
   geoflow-auto
   mg-auto
   mg-manual
   mg-para
   mg-dummy
   pbam-auto
   pbsam-auto
.. _tabi:

tabi
==========

.. currentmodule::  apbs.input_file.calculate.boundary_element

.. note::
   
   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`); see :class:`BoundaryElement` and :class:`TABIParameters`.

This mode uses the TABI-PB integral equation software from Geng and Krasny to solve the linearized Poisson-Boltmzann equation. 
Boundary element methods offer the ability to focus numerical effort on a much smaller region of the problem domain:  the interface between the molecule and the solvent.
In this method, two coupled integral equations defined on the solute-solvent boundary define a mathematical relationship between the electrostatic surface potential and its normal derivative with a set of integral kernels consisting of Coulomb and screened Coulomb potentials with their normal derivatives.
The boundary element method requires a surface triangulation, generated by a program such as `NanoShaper <https://www.electrostaticszone.eu/downloads>`_, on which to discretize the integral equations.

For more information, see the Geng & Krasny `2013 J Comput Phys paper <https://doi.org/10.1016/j.jcp.2013.03.056>`_.

.. toctree::
   :maxdepth: 1
   :caption: ELEC tabi keywords

   ion
   mac
   mesh
   ../generic/mol
   outdata
   pdie
   ../generic/sdens
   sdie
   ../generic/srad
   ../generic/temp
   tree_n0
   tree_order


======================
Background information
======================

The Treecode-Accelerated Boundary Integral Poisson-Boltzmann solver (TABI-PB; `Geng, 2013 <http://www.sciencedirect.com/science/article/pii/S0021999113002404>`_ calculates electrostatics of solvated biomolecules governed by the linearized Poisson-Boltzmann equation.
It uses a well-posed boundary integral Poisson-Boltzmann formulation to ensure rapid convergence.
In addition, a fast treecode algorithm for the screened Coulomb potential `(Li, 2009) <http://www.sciencedirect.com/science/article/pii/S0021999109000916>`_ is applied to speed up the matrix-vector products in each GMRES iteration.
The molecular surfaces, which divide the entire domain into solute region and solvent region, are generated by `NanoShaper <https://www.electrostaticszone.eu/downloads>`_.

-----------------
TABI-PB algorithm
-----------------

The coupled integral equations derived from the linearized Poisson-Boltzmann equation are

.. math::
   
   \frac{1}{2}(1+\epsilon)\phi({x})-\int_\Gamma(K_1({x},{y})\frac{\partial{\phi({y})}}{\partial{v}}+K_2({x},{y})\phi({y}))dS_{y} &=S_1({x}) \\
   \frac{1}{2}(1+\frac{1}{\epsilon})\frac{\partial{\phi({x})}}{\partial{v}}-\int_\Gamma(K_3({x},{y})\frac{\partial{\phi({y})}}{\partial{v}}+K_4({x},{y})\phi({y}))dS_{y} &=S_2({x}), {x}\in\Gamma

for the surface potential :math:`\phi`, and its normal derivative :math:`\frac{\partial\phi}{\partial v}` on the surface :math`\Gamma`.
The kernels :math:`K_{1,2,3,4}` are linear combinations of the Coulomb and screened Coulomb potentials:

.. math::
   
   G_0({x},{y}) &= \frac{1}{4\pi |{x}-{y}|} \\
   G_{\kappa}({x},{y}) &= \frac{e^{-\kappa|{x}-{y}|}}{4\pi |{x}-{y}|}

and their first and second derivatives.

The sums in the discretized form of the integral equations above have the form of *N*-body interactions,

.. math::
   
   V_i = \sum_{j=1,j \neq i}^{N} q_j G({x}_i,{x}_j), i=1,\ldots,N 

where :math:`G` is the screened Coulomb potential kernel, :math:`{x}_i, {x}_j` are the centroids of the triangles, and :math:`q_j` is the charge at :math:`{x}_j`.
The particles (centroids of the triangles) are divided into a hierarchy of clusters having a tree structure.
The treecode replaces the :math:`\mathcal{O}(N^2)` particle-particle interactions by :math:`\mathcal{O}(N \log N)` particle-cluster interactions and TABI-PB utilizes this feature efficiently.

======
Output
======

The TABI-PB code produces an output file called :file:`surface_potential.dat` containing:

* number of nodes, number of triangles
* node index, vertices, normal vector, surface potential (kJ mol\ :sup:`-1` e\ :sub:`c`\ :sup:`-1`), surface potential normal derivatives (kJ mol\ :sup:`-1` e\ :sub:`c`\ :sup:`-1` A\ :sup:`-1`)
* connectivity data for surface triangulation

The format is given below:

.. code-block:: bash
   
   num_node num_triangle
   node_index x y z norm_x norm_y norm_z phi norm_phi
   (et cetera)
   node_index1 node_index2 node_index3

The TABI-PB code prints the free energy of solvation and Coulombic free energy in kJ/mol, along with some other information such as CPU time and the GMRES residuals at each step.

Additionally, TABI-PB can optionally output a VTK polygonal data file containing color mappable potentials and normal derivatives of potentials on the faces and vertices of the mesh.
The VTK file can be visualized using `ParaView <https://www.paraview.org/>`_.

.. _dime:

dime
====

.. currentmodule:: apbs.input_file.calculate.finite_difference

.. note::  

   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`); see :class:`GridDimensions` for more information.


Specifies the number of grid points per processor for grid-based discretization.
The syntax is:

.. code-block:: bash
   
   dime {nx ny nz}

For :ref:`mgmanual` calculations, the arguments are dependent on the choice of :ref:`nlev` by the formula: :math:`n = c 2^{l + 1} + 1` where *n* is the dime argument, *c* is a non-zero integer, *l* is the :ref:`nlev` value.
The most common values for grid dimensions are 65, 97, 129, and 161 (they can be different in each direction); these are all compatible with a :ref:`nlev` value of 4.
If you happen to pick a "bad" value for the dimensions (i.e., mismatch with :ref:`nlev`), the APBS code will adjust the specified :ref:`dime` downwards to more appropriate values.
This means that "bad" values will typically result in lower resolution/accuracy calculations!
The arguments for this keyword are:

``nx ny nz``
  The (integer) number of grid points in the x-, y-, and z-directions, respectively.

.. note::
   dime should be interpreted as the number of grid points per processor for all calculations, including :ref:`mgpara`.
   This interpretation helps manage the amount of memory per-processor - generally the limiting resource for most calculations.

.. _sdie:

sdie
====

.. currentmodule::  apbs.input_file.calculate

.. note:: 

   Some versions of this command have been ported to the *new APBS syntax* (see :ref:`new_input_format`):

   * Boundary element Poisson-Boltzmann calculations: see :func:`boundary_element.BoundaryElement.solvent_dielectric` for more information.
   * Finite difference Poisson-Boltzmann calculations: see :func:`finite_difference.FiniteDifference.solvent_dielectric` for more information.
   * Finite element Poisson-Boltzmann calculations: see :func:`finite_element.FiniteElement.solvent_dielectric` for more information.

.. todo:: port for other calculation types

Specify the dielectric constant of the solvent.
The syntax is:

.. code-block:: bash
   
   sdie {diel}

where ``diel`` is a floating point number representing the solvent dielectric constant (unitless).
This number must be :math:`\ge 1`.
Bulk water at biologically-relevant temperatures is usually modeled with a dielectric constant of 78-80.
.. _term:

term
====

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

Specify a termination condition for a PB-(S)AM Brownian dynamics trajectory.
The syntax is:

.. code-block:: bash
   
   term {type} {options}

where the ``options`` are determined by the ``type`` as follows:

``contact {file}``
  Termination based on molecular contact conditions.
  ``file`` is a string for the contact file filename.
  The contact file has a list formatted as follows:  ``moltype1 at1 moltype2 at2 dist`` where ``moltype1``  and ``moltype2``  are indices of the molecular types, ``at1`` is the index of an atom from the first molecular type, ``at2`` is the index of an atom from the second molecular type, and ``dist`` is the maximum distance between the two atoms that defines the contact.
  ``pad`` is distance criterion that will be checked in the case that the true atom contact distance may not be fulfilled.

  .. note::

     Sometimes these distances cannot be reached due to the assumption in this model that the molecule is spherical.
     If this is the case, the atom positions are transformed to the molecule surface and surface points are compared to the pad distance.

``{pos} {val} {molecule}``
  Specify a position termination condition for a given molecule.
  where ``pos`` is one of the following options: ``x<=, x>=, y<=, y>=, z<=, z>=, r<=, r>=``.
  ``val`` is the value along the given axis to check against.
  ``molecule`` is the molecule index (1 based) according to the order of molecules listed in the ``READ`` section that this condition applies to.
  This command can be understood as:  "Terminate the simulation when molecule ``molecule`` fulfills the condition ``pos`` ``val``".

  .. todo::

     Add a constant keyword (e.g., like ``position``) before the ``{pos}`` argument of ``term``.
     Documented in https://github.com/Electrostatics/apbs/issues/503

``time {val}``
  Specify a time termination condition where ``val`` is a floating point number for the trajectory time limit (in picoseconds).
.. _mac:

mac
===

.. currentmodule:: apbs.input_file.calculate.boundary_element

.. note::

   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`); see :func:`TABIParameters.multipole_acceptance_criterion`.

TABI-PB parameter, multipole acceptance criterion (MAC), that controls distance ratio at which the method uses direct summation or Taylor approximation (a particle-cluster interaction) to calculate the integral kernels.
The syntax is:

.. code-block:: bash

   mac {theta}

where ``theta`` is a floating-point number from 0 to 1 controlling the distance ratio.
This multipole acceptance criterion (MAC) is :math:`\frac{r_c}{R}\leqslant \theta`, where :math:`r_c` is the cluster radius, and :math:`R` is the distance of the particle to the cluster center.
If the above relationship is satisfied, the Taylor approximation will be used instead of direct summation.
A typical value for this parameter is 0.8.
.. _usemesh:

usemesh
=======

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

Specify the external finite element mesh to be used in the finite element Poisson-Boltzmann calculation (:ref:`femanual`).
These must have been input via an earlier READ mesh statement (see :ref:`read_old_input`).
The syntax is:

.. code-block:: bash

   usemesh {id}

where ``id`` is an integer ID specifying the particular map read in with :ref:`read_old_input`.
These IDs are assigned sequentially, starting from 1, and incremented independently for each mesh read by APBS.
.. _tolsp:

tolsp
=====

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

This is an undocumented parameter from the PB-SAM code that does something with the coarseness of the molecule representation.
The PB-SAM authors recommend a value of 2.5.

.. todo::
   
   We need documentation for ``tolsp``.
   Documented in https://github.com/Electrostatics/apbs/issues/504
.. _ekey:

ekey
====

.. currentmodule::  apbs.input_file.calculate.finite_element

.. note::  This command has been ported to the *new APBS syntax* (see :func:`FiniteElement.error_based_refinement`).

Specify the method used to determine the error tolerance in the solve-estimate-refine iterations of the finite element solver (:ref:`femanual`).
The syntax is:

.. code-block:: bash

   ekey { flag }

where ``flag`` is a text string that determines the method for error calculation.

``simp``
  Per-simplex error limit
``global``
  Global (whole domain) error limit
``frac``
  Fraction of simplices you'd like to see refined at each iteration
.. _fgcent:

fgcent
======

.. currentmodule:: apbs.input_file.calculate.finite_difference

.. note::  

   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`):  see :func:`Focus.fine_grid_center` for more information.


Specify the center of the fine grid (in a focusing calculation) based on a molecule's center or absolute coordinates for :ref:`mgpara` and :ref:`mgauto` multigrid calculations.
The syntax is:

.. code-block
   
   fgcent { mol id | xcent ycent zcent }

where a user can specify **either**

``mol {id}``
  Center the grid on molecule with integer ID id; as assigned in the READ section (see :ref:`read_old_input`) of the input file.
  Molecule IDs are assigned in the order they are read, starting at 1.

**or** the user can specify

``xcent ycent zcent``
  Center the grids on the coordinates (floating point numbers in Å) at which the grid is centered.
  Based on the input molecule PDB coordinate frame.

.. _pdime:

pdime
=====

.. currentmodule::  apbs.input_file.calculate.finite_difference

.. note::  This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`); see :func:`ParallelFocus.processor_array`.


Specify the processor array to be used in a parallel focusing (:ref:`mgpara`) calculation.
The syntax is:

.. code-block:: bash
   
   pdime {npx npy npz}

where ``npx npy npz`` are the integer number of processors to be used in the x-, y- and z-directions of the system.
The product ``npx × npy × npz`` should be less than or equal to the total number of processors with which APBS was invoked (usually via mpirun).
If more processors are provided at invocation than actually used during the run, the extra processors are not used in the calculation.
The processors are tiled across the domain in a Cartesian fashion with a specified amount of overlap (see :ref:`ofrac`) between each processor to ensure continuity of the solution.
Each processor's subdomain will contain the number of grid points specified by the dime keyword.
For broad spatial support of the splines, every charge included in partition needs to be at least 1 grid space (:ref:`chgm` ``spl0``), 2 grid spaces (:ref:`chgm` ``spl2``), or 3 grid spaces (:ref:`chgm` ``spl4``) away from the partition boundary.
.. _nrpbe:

nrpbe
=====

.. currentmodule:: apbs.input_file.calculate

.. note::  

   Some aspects of this command have been moved to the *new APBS syntax* (see :ref:`new_input_format`): 

   * Finite difference:  see :func:`finite_difference.FiniteDifference.equation` for more information.

   * Finite element:  see :func:`finite_element.FiniteElement.equation` for more information.

.. todo::  port for other types of calculations.

Specifies that the nonlinear form of the regularized Poisson-Boltzmann equation (RPBE) should be solved.
The regularized PBE equation replaces the point charge distribution with the corresponding Green's function.
As a result of this replacement, the solution corresponds to the reaction field instead of the total potential; the total potential can be recovered by adding the appropriate Coulombic terms to the solution.
Likewise, this equation immediately yields the solvation energy without the need for reference calculations.

.. note::

   The options :ref:`lpbe`, :ref:`npbe`, :ref:`lrpbe`, :ref:`nrpbe` are mutually exclusive.
   
.. note::

   This functionality is only available with FEM-based solvers.
.. _tree_n0:

tree_n0
=======

.. currentmodule::  apbs.input_file.calculate.boundary_element

.. note::  This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`); see :func:`TABIParameters.maximum_particles`.

TABI-PB parameter that specifies the maximum number of particles in a treecode leaf.
This controls leaf size in the process of building the tree structure.
The syntax is:

.. code-block:: bash

   tree_n0 {max_number}

where ``max_number`` is an integer.
A typical value for this parameter is 500.
.. _dx:

dx
==

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).


Specify the name of the file into which the potential will be printed.

.. code-block:: bash
   
   dx {filename}

where ``filename`` is a string for the name of the file where an OpenDX file will be printed out.

.. todo::
   
   The PB-(S)AM ``dx`` keyword should not exist; please replace it ASAP with the :ref:`old_write` command.
   Documented in https://github.com/Electrostatics/apbs/issues/488
.. _runtype:

runtype
=======

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

Indicate what type of calculation you would like to run with the PB-(S)AM model.

.. code-block:: bash
   
   runtype {type}

where ``type`` is the type of calculation to be perfomed:

``energyforce``
  Compute and print out the interaction energies, forces and torques on each molecule. 

``electrostatics``
  Print the electrostatic potential of points in the system.

``dynamics``
  Perform a Brownian Dynamics simulation, using forces and torques generated from the PB-(S)AM model.
  The calculation of force and torque has been integrated into a Brownian dynamics scheme that is detailed in `Yap EH, Head-Gordon TL (2013) <http://pubs.acs.org/doi/abs/10.1021/ct400048q>`_
  This option will generate a series of files of the form

  :file:`dyn_{toy}.pqr`
    The starting configuration of the system for the first trajectory

  :file:`dyn_{toy}.stat`
    A file that prints how each trajectory was terminated and the time that this occurred at.

  :file:`dyn_{toy}_traj.xyz`
    A VMD-readable xyz file for the trajectory of ``traj``.

  :file:`dyn_toy_traj.dat`
    A file with positions, forces and torques for the system.

  .. todo::

     The dynamics part of the PB-(S)AM code should be moved out of the ``ELEC`` section.
     Documented in https://github.com/Electrostatics/apbs/issues/500
gamma
=====

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

This keyword specifies the surface tension coefficient for apolar solvation models.

.. code-block:: bash

   gamma { value }

where ``value`` is a floating point number designating the surface tension in units of kcal mol\ :superscript:`-1` Å\ :superscript:`-2`.
This term can be set to zero to eliminate the :abbr:`SASA (solvent-accessible surface area)` contributions to the apolar solvation calculations.

.. warning::

   *Either* this documentation is incorrect *or* the implementation needs to be changed to use kJ mol\ :superscript:`-1` Å\ :superscript:`-2` instead of kcal.

.. todo::

   Resolve unit confusion with geometric flow :ref:`gamma` keyword.
   https://github.com/Electrostatics/apbs/issues/490.. _mgdummy:

mg-dummy
========

.. currentmodule::  apbs.input_file.calculate.finite_difference

.. note::
   
   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`); see :func:`FiniteDifference.noop` for more information.


Not a Poisson-Boltzmann calculation.
Many calculations of surface and charge distribution properties which do not require solution of the PBE.

This type of calculation allows users to write out dielectric, ion-accessibility, and charge distribution, and other types of maps that depend solely on biomolecular geometry.
Since these maps depend only on geometry, they can be written out without actually solving the PB equation. 

.. toctree::
   :maxdepth: 2
   :caption: ELEC mg-auto keywords:

   bcfl
   chgm
   dime
   gcent
   glen
   ../generic/grid
   ion
   lpbe
   lrpbe
   ../generic/mol
   npbe
   pdie
   ../generic/sdens
   sdie
   ../generic/srad
   srfm
   ../generic/swin
   ../generic/temp
   write
.. _exp:

exp
===

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

This keyword can be used to load in the expansion matrices from files.
They will have been previously generated, and will be named :file:`mol{m}.{H, F}.[s].exp` (see :ref:`pbamauto` for more information).
The syntax is:

.. code-block:: bash
   
   exp {prefix}

where ``prefix`` is the filename prefix :file:`mol{m}sph`.
The *H* or *F* and :file:`{s}.bin` will be appended during the program run.

.. todo::

   It would be better to generalize the :ref:`read_old_input` section of the input file rather than use the ``exp`` command.
   This command also needs to be cleaned up -- it's too fragile.
   Documented at https://github.com/Electrostatics/apbs/issues/489.. _elecsrfm:

srfm (elec)
===========

.. currentmodule::  apbs.input_file.calculate

.. note::

   Some uses of this command have been ported to the *new APBS syntax* (see :ref:`new_input_format`):
   
   * Finite differences; see :func:`finite_difference.FiniteDifference.surface_method` for more information
   * Finite elements; see :func:`finite_element.FiniteElement.surface_method` for more information

Specify the model used to construct the dielectric and ion-accessibility coefficients.
The syntax for this command is:

.. code-block:: bash

   srfm {flag}

where ``flag`` is a string describing the coefficient model:

``mol``
  The dielectric coefficient is defined based on a molecular surface definition.
  The problem domain is divided into two spaces.
  The "free volume" space is defined by the union of solvent-sized spheres (see :ref:`srad`) which do not overlap with biomolecular atoms.
  This free volume is assigned bulk solvent dielectric values.
  The complement of this space is assigned biomolecular dielectric values.
  With a non-zero solvent radius (srad), this choice of coefficient corresponds to the traditional definition used for PB calculations.
  When the solvent radius is set to zero, this corresponds to a van der Waals surface definition.
  The ion-accessibility coefficient is defined by an "inflated" van der Waals model.
  Specifically, the radius of each biomolecular atom is increased by the radius of the ion species (as specified with the :ref:`ion` keyword).
  The problem domain is then divided into two spaces.
  The space inside the union of these inflated atomic spheres is assigned an ion-accessibility value of 0; the complement space is assigned bulk ion accessibility values.

``smol``
  The dielectric and ion-accessibility coefficients are defined as for mol (see above).
  However, they are then "smoothed" by a 9-point harmonic averaging to somewhat reduce sensitivity to the grid setup as described by Bruccoleri et al. J Comput Chem 18 268-276, 1997 (`10.1007/s00214-007-0397-0 <http://dx.doi.org/10.1007/s00214-007-0397-0>`_).

``spl2``
  The dielectric and ion-accessibility coefficients are defined by a cubic-spline surface as described by Im et al, Comp Phys Commun 111 (1-3) 59-75, 1998 (`10.1016/S0010-4655(98)00016-2 <https://doi.org/10.1016/S0010-4655(98)00016-2>`_).
  The width of the dielectric interface is controlled by the :ref:`swin` parameter.
  These spline-based surface definitions are very stable with respect to grid parameters and therefore ideal for calculating forces.
  However, they require substantial reparameterization of the force field; interested users should consult Nina et al, Biophys Chem 78 (1-2) 89-96, 1999 (`10.1016/S0301-4622(98)00236-1 <http://dx.doi.org/10.1016/S0301-4622(98)00236-1>`_).
  Additionally, these surfaces can generate unphysical results with non-zero ionic strengths; this is an on-going area of development.

``spl4``
  The dielectric and ion-accessibility coefficients are defined by a 7th order polynomial.
  This surface definition has characteristics similar to spl2, but provides higher order continuity necessary for stable force calculations with atomic multipole force fields (up to quadrupole).
press
=====

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

This term specifies the solvent pressure in kJ mol\ :superscript:`-1` Å\ :superscript:`-3`.
This coefficient multiplies the volume term of the apolar model and can be set to zero to eliminate volume contributions to the apolar solvation calculation.
The syntax is:

.. code-block:: bash

   press {value}

where ``value`` is the floating point value of the pressure coefficient in kJ mol\ :superscript:`-1` Å\ :superscript:`-3`.

.. warning::

   *Either* this documentation is incorrect *or* the implementation needs to be changed to use kJ mol\ :superscript:`-1` Å\ :superscript:`-3` instead of kcal.

.. todo::

   Resolve unit confusion with geometric flow ``press`` keyword and the apolar :ref:`press` keyword.
   Documented in https://github.com/Electrostatics/apbs/issues/499
.. _etol:

etol
====

.. currentmodule:: apbs.input_file.calculate

.. note::  

   Some versions of this command have been ported to the *new APBS syntax* (see :ref:`new_input_format`):

   * For a boundary element calculation, see :func:`boundary_element.BoundaryElement.error_tolerance`.
   * For a finite difference calculation, see :func:`finite_difference.FiniteDifference.error_tolerance`.
   * For a finite element calculation, see :func:`finite_element.FiniteElement.error_tolerance`.

.. todo::  add documentation links for other instances.

Specifies the tolerance for iterations of the partial differntial equation solvers:
The syntax is:

.. code-block:: bash
   
   etol { tol }

where ``tol`` is the (floating point) numerical value for the error tolerance.

For finite difference solvers, this keyword is optional and is intended for :ref:`mgmanual`, :ref:`mgauto`, and :ref:`mgpara` calculation types.

For finite element solvers, this keyword specifies the tolerance for error-based adaptive refinement during the solve-estimate-refine iterations of the finite element solver (:ref:`femanual`), where ``tol`` is the (floating point) numerical value for the error tolerance.
.. _outdata:

outdata
=======

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

TABI-PB parameter that specifies the file type for printing the output data.
The syntax is:

.. code-block:: bash
   
   outdata {flag}

where ``flag`` is an integer indicating the output file types:

0
  .dat format
1
  Both the .dat format and a VTK polygonal data file that can be visualized in the ParaView software.
  The VTK file contains color mappable potentials and normal derivatives of potentials on the faces and vertices of the mesh.


.. todo::

   The integer flag values for ``mesh`` should really be replaced by human-readable strings.
.. _fglen:

fglen
=====


.. currentmodule:: apbs.input_file.calculate.finite_difference

.. note::  

   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`):  see :func:`Focus.fine_grid_dimensions` for more information.


Specifies the fine mesh domain lengths in a multigrid focusing calculation (:ref:`mgpara` or :ref:`mgauto`); this may be different in each direction.
The syntax is:

.. code-block:: bash

   fglen {xlen ylen zlen}

This should enclose the region of interest in the molecule.
The arguments to this command are:

``xlen ylen zlen``
  Grid lengths (floating point numbers) in the x-, y-, and z-directions in Å.

.. _mgauto:

mg-auto
=======

.. currentmodule:: apbs.input_file.calculate.finite_difference

.. note::

   This command has been ported to the *new APBS syntax* (see :ref:`new_input_format`).
   See :class:`FiniteDifference` and :class:`Focus` for more information.

Automatically configured finite difference Poisson-Boltzmann calculations.

This multigrid calculation automatically sets up and performs a string of single-point PBE calculations to "focus" on a region of interest (binding site, etc.) in a system.
It is basically an automated version of :ref:`mgmanual` designed for easier use.
Most users should use this version of ELEC.

Focusing is a method for solving the Poisson-Boltzmann equation in a finite difference setting.
Some of the earliest references to this method are from Gilson and Honig [#Gilson]_.
The method starts by solving the equation on a coarse grid (i.e., few grid points) with large dimensions (i.e., grid lengths).
The solution on this coarse grid is then used to set the Dirichlet boundary condition values for a smaller problem domain -- and therefore a finer grid -- surrounding the region of interest.
The finer grid spacing in the smaller problem domain often provides greater accuracy in the solution.

The following keywords are present in mg-auto ELEC blocks; all keywords are required unless otherwise noted.

.. note::

   During focusing calculations, you may encounter the message "WARNING! Unusually large potential values detected on the focusing boundary!" for some highly charged systems based on location of the focusing boundary.
   First, you should determine if you received any other warning or error messages as part of this calculation, particularly those referring to exceeded number of iterations or error tolerance (:ref:`etol`). 
   Next, you should check if the calculation converged to a reasonable answer.
   In particular, you should check sensitivity to the grid spacing by making small changes to the grid lengths (via the :ref:`fglen` parameter) and see if the changes in energies are correspondingly small.
   If so, then this warning can be safely ignored.


.. toctree::
   :maxdepth: 2
   :caption: ELEC mg-auto keywords:

   bcfl
   ../generic/calcenergy
   ../generic/calcforce
   cgcent
   cglen
   chgm
   dime
   etol
   fgcent
   fglen
   ion
   lpbe
   lrpbe
   ../generic/mol
   npbe
   pdie
   ../generic/sdens
   sdie
   ../generic/srad
   srfm
   ../generic/swin
   ../generic/temp
   usemap
   write
   writemat

.. [#Gilson] Gilson MK and Honig BH, Calculation of electrostatic potentials in an enzyme active site. Nature, 1987. 330(6143): p. 84-6. DOI:`10.1038/330084a0 <http://dx.doi.org/10.1038/330084a0>`_.. _lpbe:

lpbe
====

.. currentmodule:: apbs.input_file.calculate

.. note::  

   Some aspects of this command have been moved to the *new APBS syntax* (see :ref:`new_input_format`): 

   * Finite difference:  see :func:`finite_difference.FiniteDifference.equation` for more information.

   * Finite element:  see :func:`finite_element.FiniteElement.equation` for more information.

.. todo::  port for other types of calculations.


Specifies that the linearized Poisson-Boltzmann equation should be solved.

.. note::

   The options :ref:`lpbe`, :ref:`npbe`, :ref:`lrpbe`, :ref:`nrpbe` are mutually exclusive... _bcfl:

bcfl
====

.. currentmodule::  apbs.input_file.calculate

.. note::  

   Some versions of this command have been ported to the *new APBS syntax* (see :ref:`new_input_format`):
   
   * Finite difference boundary conditions; see :func:`finite_difference.FiniteDifference.boundary_condition`.
   * Finite element boundary conditions; see :func:`finite_element.FiniteElement.boundary_condition`.

Specifies the type of boundary conditions used to solve the Poisson-Boltzmann equation.
The syntax is:

.. code-block:: bash
   
   bcfl {flag}

where ``flag`` is a text string that identifies the type of conditions to be used.

``zero``
  "Zero" boundary condition. Dirichlet conditions where the potential at the boundary is set to zero.
  This condition is not commonly used and can result in large errors if used inappropriately.
``sdh``
  "Single Debye-Hückel" boundary condition.
  Dirichlet condition where the potential at the boundary is set to the values prescribed by a Debye-Hückel model for a single sphere with a point charge, dipole, and quadrupole.
  The sphere radius in this model is set to the radius of the biomolecule and the sphere charge, dipole, and quadrupole are set to the total moments of the protein.
  This condition works best when the boundary is sufficiently far from the biomolecule.
``mdh``
  "Multiple Debye-Hückel" boundary condition.
  Dirichlet condition where the potential at the boundary is set to the values prescribed by a Debye-Hückel model for a multiple, non-interacting spheres with a point charges.
  The radii of the non-interacting spheres are set to the atomic radii of and the sphere charges are set to the atomic charges.
  This condition works better than sdh for closer boundaries but can be very slow for large biomolecules.
``focus``
  "Focusing" boundary condition.
  Dirichlet condition where the potential at the boundary is set to the values computed by the previous (usually lower-resolution) PB calculation.
  This is **only** used in sequential focusing performed manually in :ref:`mgmanual` calculations.
  All of the boundary points should lie within the domain of the previous calculation for best accuracy; if any boundary points lie outside, their values are computed using single Debye-Hückel boundary conditions (see above).
``map``
  Specifying map allows a previously calculated potential map to be used in a new focusing calculation.
  A typical scenario is using the same coarse grid for multiple focusing calculations.
  A potential map can be written once from a coarse grid calculation, then used in subsequent runs to bypass the need to recalculate the coarse grid.
  See the READ keyword pot (see :ref:`read_old_input`) and the attached example files for its use.


.. _diff:

diff
====

.. todo::  This command has not yet been ported to the *new APBS syntax* (see :ref:`new_input_format`).

Specify the diffusion coefficients for each molecule in the system for a PB-(S)AM Brownian dynamics calculation.

.. code-block:: bash

   diff {type} {dTrans} {dRot}

``type``
  a string indicating the molecule dynamics type
  
  ``stat``
    Stationary.

  ``rot``
    Object is fixed but rotates

  ``move``
    Object moves and rotates.

``dTrans``
  Translational diffusion coefficient in units of Å\ :sup:`2` ps\ :sup:`-1`.
  Used only with the ``move`` keyword.

``dRot``
  Rotational diffusion coefficient.
  Used with the ``move`` and ``rot`` keywords.

.. todo::
   
   What are the units for ``dRot``?
   Documented as https://github.com/Electrostatics/apbs/issues/486

.. note::

   The order of these keywords is expected to be identical to the order of the molecules in the READ section.

.. todo::
   
   Add a ``mol id`` flag rather than have an implicit ordering of the ``diff`` keywords.
   Documented in https://github.com/Electrostatics/apbs/issues/487
.. _akeySOLVE:

akeySOLVE
=========

.. note::  This command has been eliminated in the *new APBS syntax* (see :ref:`new_input_format`) and will not be ported.

Specifies how the the finite element mesh should be adaptively subdivided during the solve-estimate-refine iterations of a :ref:`femanual` finite element calculation.
The syntax is:

.. code-block:: bash

   akeySOLVE {key}

where ``key`` is a text string that specifies the method used to guide adaptive refinement:

``resi``
  Residual-based a *posteriori* refinement.

==============================
Data processing input file API
==============================

.. automodule::  apbs.input_file.process
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:.. _new_input_format:

=================================
YAML- and JSON-format input files
=================================

In its new input format, APBS accepts either `JSON <json.org>`_- or `YAML <yaml.org>`_-format input files.

These input files consist of the following keywords and objects:

.. toctree::
   :maxdepth: 3

   read
   calculate/index
   process

.. todo:: finish other required and optional sections

^^^^^^^^^^^^^^^^^^^^^^^^^^
Input file class structure
^^^^^^^^^^^^^^^^^^^^^^^^^^

The input file object parsing and validation follows the basic pattern implemented in :class:`apbs.input_file.InputFile` (see below).
This class should serve as a template for adding new input file sections.

.. automodule::  apbs.input_file
   :members:
   :undoc-members:
   :show-inheritance:.. _read_new_input:

==========================================
Data loading input file section (required)
==========================================

.. currentmodule:: apbs.input_file.read

This required section is denoted by the keyword ``read`` and is described in :class:`Read`.
The section includes objects indexed by the following keywords:

* ``molecules``:  a list of molecule input objects of class :class:`Molecule`
* ``potential maps``:  a list of electrostatic potential map input objects of class :class:`Map`
* ``charge density maps``:  a list of charge density map input objects of class :class:`Map`
* ``ion accessibility maps``:  a list of ion accessibility map input objects of class :class:`Map`
* ``dielectric maps``:  a list of dielectric map input objects of class :class:`DielectricMapGroup`
* ``parameters``:  a list of parameter files of class :class:`Parameter`

--------------------------------------------------------------------------------

---------------------------
Data loading input file API
---------------------------

.. automodule::  apbs.input_file.read
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:
.. _process_new_input:

=====================================
Results processing section (optional)
=====================================


This optional section is denoted by the keyword ``process`` and includes lists of arithmetic operation objects (see :class:`apbs.input_file.process.Operation`) indexed by the following keywords:

* ``sums``:  a list weighted sum operations
* ``products``:  a list of weighted product operations
* ``exps``:  a list of element-wise exponentiation operations

The syntax for these objects is described in :class:`apbs.input_file.process.Process`.

See also:

.. toctree::
   :maxdepth: 2

   process_api
.. _finite_difference:

================================================
Finite-difference Poisson-Boltzmann calculations
================================================

.. todo::

   Make this more user-friendly by:

   * Adding introductory text with contents
   * Only showing inherited members in API documentation
   * Only showing undocumented members in API documentation

.. automodule::  apbs.input_file.calculate.finite_difference
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:
.. _finite_element:

================================================
Finite-element Poisson-Boltzmann calculations
================================================

.. todo::  FINISH

.. todo::

   Make this more user-friendly by:

   * Adding introductory text with contents
   * Only showing inherited members in API documentation
   * Only showing undocumented members in API documentation

.. automodule::  apbs.input_file.calculate.finite_element
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:
.. _boundary_element:

================================================
Boundary-element Poisson-Boltzmann calculations
================================================

This input file section documents boundary element method calculations with `TABI-PB <https://doi.org/10.1016/j.jcp.2013.03.056>`_.

.. todo::  Add PyGBE

.. todo::

   Make this more user-friendly by:

   * Adding introductory text with contents
   * Only showing inherited members in API documentation
   * Only showing undocumented members in API documentation

.. automodule::  apbs.input_file.calculate.boundary_element
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:
.. _calculate_new_input:

=========================================
Calculation input file section (required)
=========================================

This required section is denoted by the keyword ``calculate`` and includes a list of objects of the type :class:`apbs.input_file.calculate.Calculate`.

.. todo:: improve documentation with outline and/or example.

-----------------
Calculation types
-----------------

.. toctree::
   :maxdepth: 2

   boundary_element
   finite_difference
   finite_element
   nonpolar

-------------
Calculate API
-------------

.. automodule::  apbs.input_file.calculate
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:

-------------------------------------------------------

.. automodule::  apbs.input_file.calculate.generic
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:.. _nonpolar_new_input:

================================
Nonpolar grid-based calculations
================================

.. todo::  Provide overview and/or separate API documentation from main docs.

.. automodule::  apbs.input_file.calculate.nonpolar.Nonpolar
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:
==================
Solvation energies
==================

.. todo::  Update this documentation with the *new APBS syntax* (see :ref:`new_input_format`).

Solvation energies are usually decomposed into a free energy cycle as shown in the free energy cycle below.
Note that such solvation energies often performed on fixed conformations; as such, they are more correctly called "potentials of mean force".
More details on using APBS for the polar and nonpolar portions of such a cycle are given in the following sections.

.. figure:: /media/apbs_sol_eng.png

Our model solvation free energy cycle illustrating several steps:

1. The solvation energy to be calculated.
2. Charging of the solute in solution (e.g., inhomogeneous dielectric, ions present).
3. Introduction of attractive solute-solvent dispersive interaction interactions (e.g., an integral of Weeks-Chandler-Andersen interactions over the solvent-accessible volume).
4. Introduction of repulsive solute-solvent interaction (e.g., cavity formation).
5. Basically a null step although it could be used to offset unwanted energies added in Steps 3 and 4 above.
6. Charging of the solute in a vacuum or homogeneous dielectric environment in the absence of mobile ions.

---------------
Polar solvation
---------------

The full free energy cycle is usually decomposed into polar and nonpolar parts.
The polar portion is usually represented by the charging energies in Steps 2 and 6:

.. math::

   \Delta_p G = \Delta_2 G - \Delta_6 G 

Energies returned from APBS electrostatics calculations are charging free energies.
Therefore, to calculate the polar contribution to the solvation free energy, we simply need to setup two calculations corresponding to Steps 2 and 6 in the free energy cycle.
Note that the electrostatic charging free energies returned by APBS include self-interaction terms.
These are the energies of a charge distribution interacting with itself.
Such self-interaction energies are typically very large and extremely sensitive to the problem discretization (grid spacing, location, etc.).
Therefore, it is very important that the two calculations in Steps 2 and 6 are performed with identical grid spacings, lengths, and centers, in order to ensure appropriate matching (or "cancellation") of self-energy terms.

^^^^^^^^
Born ion
^^^^^^^^

One of the canonical examples for polar solvation is the Born ion: a nonpolarizable sphere with a single charge at its center surrounded by an aqueous medium.
Consider the transfer of a non-polarizable ion between two dielectrics.
In the initial state, the dielectric coefficient inside and outside the ion is :math:`\epsilon\_{\mathrm {in}}`, and in the final state, the dielectric coefficient inside the ion is :math:`\epsilon\_{\mathrm {in}}` and the dielectric coefficient outside the ion is :math:`\epsilon\_{\mathrm {in}}`.
In the absence of external ions, the polar solvation energy of this transfer for this system is given by:

.. math::
   
   \Delta\_p G\_{\mathrm{Born}}= \frac{q^2}{8\pi\epsilon\_0 a}\left (\frac{1}{\epsilon\_{\mathrm {out}}}-\frac{1}{\epsilon\_{\mathrm {in}}}\right)

where q is the ion charge, a is the ion radius, and the two ε variables denote the two dielectric coefficients.
This model assumes zero ionic strength.

Note that, in the case of transferring an ion from vacuum, where :math:`\epsilon\_{\mathrm {in}} = 1`, the expression becomes

.. math::

   \Delta\_p G\_{\mathrm{Born}}= \frac{q^2}{8\pi\epsilon\_0 a}\left (\frac{1}{\epsilon\_{\mathrm {out}}}-1\right)

We can setup a PQR file for the Born ion for use with APBS with the contents:

.. code-block:: bash

   REMARK  This is an ion with a 3 A radius and a +1 e charge
   ATOM      1   I  ION     1 0.000   0.000   0.000  1.00 3.00

We're interested in performing two APBS calculations for the charging free energies in homogeneous and heterogeneous dielectric coefficients.
We'll assume the internal dielectric coefficient is 1 (e.g., a vacuum) and the external dielectric coefficient is 78.54 (e.g., water).
For these settings, the polar Born ion solvation energy expression has the form

.. math::

   \Delta_p G_{\mathrm{Born}} = -691.85 \biggl( \frac{z^2}{R} \biggr) \mathrm {kJ \, A/mol}

where z is the ion charge in electrons and R is the ion size in Å.

This solvation energy calculation can be setup in APBS with the following input file:

.. code-block:: bash

   # READ IN MOLECULES
   read
     mol pqr born.pqr
   end
   elec name solv # Electrostatics calculation on the solvated state
     mg-manual # Specify the mode for APBS to run
     dime 97 97 97 # The grid dimensions
     nlev 4 # Multigrid level parameter
     grid 0.33 0.33 0.33 # Grid spacing
     gcent mol 1 # Center the grid on molecule 1
     mol 1 # Perform the calculation on molecule 1
     lpbe # Solve the linearized Poisson-Boltzmann equation
     bcfl mdh # Use all multipole moments when calculating the potential
     pdie 1.0 # Solute dielectric
     sdie 78.54 # Solvent dielectric
     chgm spl2 # Spline-based discretization of the delta functions
     srfm mol # Molecular surface definition
     srad 1.4 # Solvent probe radius (for molecular surface)
     swin 0.3 # Solvent surface spline window (not used here)
     sdens 10.0 # Sphere density for accessibility object
     temp 298.15 # Temperature
     calcenergy total # Calculate energies
     calcforce no # Do not calculate forces
   end
   elec name ref # Calculate potential for reference (vacuum) state
     mg-manual
     dime 97 97 97
     nlev 4
     grid 0.33 0.33 0.33
     gcent mol 1
     mol 1
     lpbe
     bcfl mdh
     pdie 1.0
     sdie 1.0
     chgm spl2
     srfm mol
     srad 1.4
     swin 0.3
     sdens 10.0
     temp 298.15
     calcenergy total
     calcforce no
   end
   # Calculate solvation energy
   print energy solv - ref end
   quit

Running this example with a recent version of APBS should give an answer of -229.59 kJ/mol which is in good agreement with the -230.62 kJ/mol predicted by the analytic formula above.

.. note::

   The Born example above can be easily generalized to other polar solvation energy calculations.
   For example, ions could be added to the solv ELEC, dielectric constants could be modified, surface definitions could be changed (in both ELEC sections!), or more complicated molecules could be examined.
   Many of the examples included with APBS also demonstrate solvation energy calculations.

.. note::

   As molecules get larger, it is important to examine the sensitivity of the calculated polar solvation energies with respect to grid spacings and dimensions.

----------------
Apolar solvation
----------------

Referring back to the solvation free energy cycle, the nonpolar solvation free energy is usually represented by the energy changes in Steps 3 through 5:

.. math::

   \Delta_n G = (\Delta_3 G - \Delta_5 G) + \Delta_4 G


where Step 4 represents the energy of creating a cavity in solution and Steps 3-5 is the energy associated with dispersive interactions between the solute and solvent.
There are many possible choices for modeling this nonpolar solvation process.
APBS implements a relatively general model described by `Wagoner and Baker (2006) <http://www.pnas.org/content/103/22/8331>`_ and references therein.
The implementation and invocation of this model is described in more in the :ref:`apolar` documentation.
Our basic model for the cavity creation term (Step 4) is motivated by scaled particle theory and has the form

.. math::

   \Delta_4 G = pV + \gamma A

where :math:`p` is the solvent pressure (:ref:`press` keyword), :math:`V` is the solute volume, :math:`\gamma` is the solvent surface tension (:ref:`gamma` keyword), and :math:`A` is the solute surface area.

Our basic model for the dispersion terms (Steps 3 and 5) follow a Weeks-Chandler-Anderson framework as proposed by `Levy et al (2002) <http://onlinelibrary.wiley.com/doi/10.1002/jcc.10045/abstract>`_:

.. math::

   \Delta_3 G - \Delta_5 G = \overset{-} \rho \int_\omega u^{(att)}(y)\theta(y)dy

where :math:`\overline{\rho}` is the bulk solvent density (:ref:`bconc` keyword), :math:`\Omega` is the problem domain, :math:`u^{\mathrm{(att)}}(y)` is the attractive dispersion interaction between the solute and the solvent at point y with dispersive Lennard-Jones parameters specified in APBS parameter files, and :math:`\theta(y)` describes the solvent accessibility of point y.

The ability to independently adjust :ref:`press`, :ref:`gamma`, and :ref:`bconc` means that the general nonpolar solvation model presented above can be easily adapted to other popular nonpolar solvation models.
For example, setting :ref:`press` and :ref:`bconc` to zero yields a typical solvent-accessible surface area model.

=====================================
Protein-RNA binding linked equilibria
=====================================

.. todo::  Update this documentation with the *new APBS syntax* (see :ref:`new_input_format`).

Before reading this example, please review :ref:`errors` for relevant caveats.

------------
Introduction
------------

This example is taken from `a paper by García-García and Draper <http://dx.doi.org/10.1016/S0022-2836\(03\)00615-6>`_.
Special thanks to `David Draper <http://pmcb.jhu.edu/inactive%20pages/draper-profile.html>`_ who provided the PDB files.
This example explores the electrostatic contributions to the binding interaction between a 22-residue α-helical peptide of protein λ with the "box B" RNA hairpin structure.
In particular, this example uses nonlinear Poisson-Boltzmann equation calculations to look at the non-specific screening effects of monovalent salt on the peptide-RNA complex.
García-García and Draper isolated the contribution of KCl concentration to the binding of the folded peptide with the folded RNA hairpin and determined a fairly linear relationship between the binding free energy :math:`\Delta_{b} G` and the logarithm of the KCl concentration which yields 

.. math::

   \frac{\partial\Delta_{b}G}{\partial\log_{10}[{\rm KCl}]} = {6.0 \pm 0.2 ~ } {\rm kcal/mol}

This slope can be used to determine the number  of KCl ions linked to the binding equilibrium through the expression

.. math::

   n = -\frac{\partial \Delta_b G}{{RT} \partial \log_{10}[{\rm KCl}]} = {-4.52 \pm 0.08~ } {\rm kcal/mol}

where :math:`RT` is the thermal energy, to determine :math:`n = -4.4 \pm 0.2` for the RNA-peptide binding equilibrium.
:math:`RT` is equal to :math:`kT * N_a` where :math:`kT` is the product of the Boltzmann constant :math:`k` (equal to the gas constant :math:`R/N_a`), and the temperature :math:`T` (at STP it is 298.15 K) and :math:`N_a` is Avogadro's constant.
Thus :math:`RT` is equal to

.. math::
   
   R ~ ({\mathrm{Joules}}/{\mathrm{Kelvin}}) * T~({\mathrm {Kelvin}}) * N_a~({\mathrm {mols}}) * {1~\mathrm{kJ}}/{1000~\mathrm J}

which roughly equals

.. math::

   (1.38 \times 10^{-23}) \times (6.022 \times 10^{23}) \times (298.15)/(1000)

which is approximately 2.479 kJ/mol or 0.593 kcal/mol.

García-García and Draper used nonlinear Poisson-Boltzmann equation calculations to estimate the electrostatic contributions to the binding free energy as a function of the monovalent salt concentration.
As :ref:`discussed elsewhere <errors>`, the Poisson-Boltzmann equation is only able to describe non-specific interactions of ions with solutes, including the effects of ion size and charge but otherwise ignoring the important differences between ionic species.
Interestingly (and perhaps surprisingly), they find excellent agreement between the experimental binding energy dependence on KCl and their Poisson-Boltzmann calculations with equivalent concentrations of monovalent ions.
This agreement strongly suggests that the binding of RNA and the peptide is primarily determined by electrostatic interactions.
It also suggests that the primary interaction of the KCl with this system is through non-specific screening interactions.
The García-García and Draper nonlinear Poisson-Boltzmann equation calculations gave:

.. math::

   \frac{\partial\Delta_{b}G}{\partial\log_{10}[{\rm KCl}]} = {5.9 \pm 0.2 ~ } {\rm kcal/mol}
 
and :math:`n = -4.3 \pm 0.2` for KCl linkage to the RNA-peptide binding equilibrium.

-------------------
APBS implementation
-------------------

This example follows the calculations from their paper.

The PQR files are included in the :file:`examples/protein-rna/` directory of the apbs repository.
This directory also includes a :file:`template.txt` file that serves as a template for the APBS input files with ``IONSTR`` as a placeholder for the ionic strength.
This file is also shown here:

.. code-block:: bash

   read  
     mol pqr model_outNB.pqr
     mol pqr model_outNpep.pqr
     mol pqr model_outBoxB19.pqr
   end
   elec name complex
     mg-auto
     dime 65 97 129
     cglen 45.3322 54.9498 82.2633
     fglen 45.3322 52.3234 68.3902
     cgcent mol 1
     fgcent mol 1
     mol 1
     npbe
     bcfl sdh
     pdie 4.0
     ion charge 1 conc IONSTR radius 2.0
     ion charge -1 conc IONSTR radius 2.0
     sdie 80.0
     srfm mol
     chgm spl2
     sdens 10.00
     srad 1.40
     swin 0.30
     temp 298.15
     calcenergy total
     calcforce no
     write qdens dx qdens-complex-IONSTR
     write ndens dx ndens-complex-IONSTR
   end
   elec name peptide
     mg-auto
     dime 65 97 129
     cglen 45.3322 54.9498 82.2633
     fglen 45.3322 52.3234 68.3902
     cgcent mol 1
     fgcent mol 1
     mol 2
     npbe
     bcfl sdh
     pdie 4.0
     sdie 80.0 
     ion charge 1 conc IONSTR radius 2.0 
     ion charge -1 conc IONSTR radius 2.0 
     srfm mol 
     chgm spl2 
     sdens 10.00 
     srad 1.40 
     swin 0.30 
     temp 298.15 
     calcenergy total 
     calcforce no 
     write qdens dx qdens-peptide-IONSTR 
     write ndens dx ndens-peptide-IONSTR 
   end 
   elec name rna 
     mg-auto 
     dime 65 97 129 
     cglen 45.3322 54.9498 82.2633 
     fglen 45.3322 52.3234 68.3902 
     cgcent mol 1 
     fgcent mol 1 
     mol 3 
     npbe 
     bcfl sdh 
     pdie 4.0 
     sdie 80.0 
     ion charge 1 conc IONSTR radius 2.0 
     ion charge -1 conc IONSTR radius 2.0 
     srfm mol 
     chgm spl2 
     sdens 10.00 
     srad 1.40 
     swin 0.30 
     temp 298.15 
     calcenergy total 
     calcforce no 
     write qdens dx qdens-rna-IONSTR 
     write ndens dx ndens-rna-IONSTR 
   end
   print elecEnergy complex - peptide - rna end 
   quit

As used in the template file, the READ command, our calculation will have three parts:  

* Calculation of the total electrostatic energy (including self-interaction energies) of the peptide-RNA complex. This calculation is named complex in the input file.  
* Calculation of the total electrostatic energy (including self-interaction energies) of the peptide. This calculation is named peptide in the input file.  
* Calculation of the total electrostatic energy (including self-interaction energies) of the RNA. This calculation is named rna in the input file.  

The calculations themselves will not be overly demanding, since we will use relatively coarse grids.
This grid coarseness has a significant impact on the absolute electrostatic binding energy we obtain from this particular calculation: the calculated energy isn't converged with respect to grid spacing.
However, the overall slope of binding energy with respect to monovalent ion concentration is rather insensitive with respect to the grid spacing, allowing us to save computational time and effort during the calculations.
The calculation will conclude with a :doc:`/using/input/old/print` command which will combine the total energies from the three parts to obtain our approximate absolute electrostatic binding energy for the complex at 0.225 M monovalent salt concentration.
It is very important to note that this absolute energy no meaning in isolation for several reasons:  

* It is not converged with respect to grid spacing  
* It does not contain other very important non-electrostatic aspects of the binding energy which are important for the measured affinity  

``IONSTR`` is a placeholder that represents the ion concentration for the APBS calculation.

You will also have to create a :file:`dxmath.txt` file which contains the following.

.. code-block:: bash

   qdens-complex-IONSTR.dx
   qdens-pep-IONSTR.dx -
   qdens-rna-IONSTR.dx -
   qdens-diff-IONSTR.dx = 

:ref:`dxmath` will subtract the dx maps of the individual peptide and RNA from the overall structure (and prints to the :file:`qdens-diff-IONSTR.dx` file.

----------------------
Automation with Python
----------------------

We have provided Python scripts :file:`apbs_{win, unix}_dx.py` that run the necessary APBS calculations and analyze the results.
When you run these programs, you need to be in the same directory as ``template.txt`` and ``dxmath.txt``.
This script will create all the input files for the tests as well as run apbs and dxmath on your :file:`template.txt` and :file:`dxmath.txt` files.
Most of the syntax fills in the ion concentrations in the template file, and the call commands actually run the calculations on each input.

-------------
Visualization
-------------

The :file:`qdens-diff-0.225.dx` file produced by the script can be viewed in PyMOL or another visualization program to give something similar to the following imaged which show the difference in charge density before and after binding.

.. image:: /media/rna-qdens-pymol.jpg

.. image:: /media/rna-qdens-vmd.jpg

========================
VIsualization with PyMOL
========================

The `PyMOL <http://www.pymol.org/>`_ molecular graphics software package can both run APBS and visualize resulting electrostatic potentials. 
Below are instructions for performing a basic demonstration of how to go from a PDB entry to a plot of structure and potential in PyMOL using APBS.

------------------------
Run the APBS calculation
------------------------

* Load your PQR file you created into PyMOL (:guilabel:`File → Open...`) and choose your favorite graphical representation of the molecular structure.

* Go to :guilabel:`Plugin → APBS Tools...` to open the APBS calculation plugin.

* Under the :guilabel:`Main` tab of the PyMOL APBS Tools window, select :guilabel:`Use another PQR` and either browse to (via the :guilabel:`Choose Externally Generated PQR` button) or input the path to your PQR file. This step is necessary to ensure you use the radii and charges assigned by PDB2PQR.

* Under the :guilabel:`APBS Location` tab of the PyMOL APBS Tools window, either browse to (via the APBS binary location: button) or input the path to your local APBS binary. It is not necessary to provide a path to the APBS :file:`psize.py` binary for most biomolecules.

* Under the :guilabel:`Temporary File Locations` tab of the PyMOL APBS Tools window, customize the locations of the various temporary files created during the run. This can be useful if you want to save the generated files for later use.

* Under the :guilabel:`Configuration` tab of the PyMOL APBS Tools window, press :guilabel:`Set grid` to set the grid spacings. The default values are usually sufficient for all but the most highly charged biomolecules.

* Under the :guilabel:`Configuration` tab of the PyMOL APBS Tools window, customize the remaining parameters; the defaults are usually OK.

  .. note::

     0.150 M concentrations for the +1 and −1 ion species are often useful to ensure that electrostatic properties are not overly exaggerated.

* Under the :guilabel:`Configuration` tab of the PyMOL APBS Tools window, press the Run :guilabel:`APBS button` to start the APBS calculation. Depending on the speed of your computer, this could take a few minutes. The :guilabel:`Run APBS` button will become unselected when the calculation is finished.

---------------------
Visualize the results
---------------------

Before proceeding, you must load the electrostatic potential data into PyMOL. Under the :guilabel:`Visualization` tab of the PyMOL APBS Tools window, press the :guilabel:`Update` button.

^^^^^^^^^^^^^^^^^^^^^^^^^
Electrostatic isocontours
^^^^^^^^^^^^^^^^^^^^^^^^^

PyMOL makes this step very easy: adjust the positive and negative "Contour" fields to the desired values (usually ±1, ±5, or ±10 kT/e)
and press the :guilabel:`Positive Isosurface`, :guilabel:`Negative Isosurface`, and :guilabel:`Show buttons`.

At this point, you probably have a figure that looks something like the image below.

.. figure:: /media/fas2-iso-pymol.png

   ±1 kT/e electrostatic potential isocontours of FAS2 in PyMOL

If the colors are not as you expect, you can change the colors of the objects iso_neg and iso_pos in the main menu.
By convention (for electrostatics in chemistry), red is negative (think oxygen atoms in carboxyl groups) and blue positive (think nitrogen atoms in amines).

^^^^^^^^^^^^^^^^^^
Surface potentials
^^^^^^^^^^^^^^^^^^

If you haven't already, hide the isocontours by pressing the :guilabel:`Positive Isosurface`, :guilabel:`Negative Isosurface`, and :guilabel:`Hide` buttons.
The surface potential is also straightforward to visualize.
Set the "Low" and "High"values to the desired values (usually ±1, ±5, or ±10 kT/e) at which the surface colors are clamped at red (-) or blue (+).
Check the "Solvent accessible surface" and "Color by potential on sol. acc. surf." buttons to plot the potential on the solvent-accessible (probe-inflated or Lee-Richards) surface.
Press the :guilabel:`Molecular Surface` :guilabel:`Show` button to load the surface potential.

.. figure:: /media/fas2-surf-pymol.png

   ±5 kT/e electrostatic potential of FAS2 in PyMOL plotted on the solvent-accessible surface.

The solvent-accessible surface tends to reveal more global features of the surface potential.
Tighter surfaces (e.g., van der Waals and molecular or Connolly surfaces) provide more information about the shape of the biomolecule but otherwise tend to simply map atomic surface charges onto the biomolecular surface.
PyMOL can simultaneously provide geometric information (from the molecular surface) and useful electrostatic potential information (from the solvent-accessible surface).
To visualize the molecule in this way, simply uncheck the "Solvent accessible surface"box and check the "Color by potential on sol. acc. surf." box on the :guilabel:`Visualization` tab.
=============================
Virtual reality with UnityMol
=============================

Molecular visualization software packages provide the ability for users to explore the 3D representations molecular structures and properties.
Typical user interaction is limited to panning, zooming, and rotating the molecule using a mouse and keyboard while viewing on a standard computing monitor.
These techniques support a pseudo 3-dimensional view of a molecule to understand its structure but lack the true depth perception people are used to with stereoscopic vision in the real world.

New advancements in virtual reality (VR) technologies has resulted in lower costs and systems that are easier to use to many consumers.
Compared to past VR hardware, these new systems have several key advancements including lower latency, higher frame rates, and improved resolution.
Additionally, these systems are equipped with better optics and motion tracking and a more robust software ecosystem.

We are extending the visualization capabilities for APBS through the incorporation of a VR device with molecular rendering software.
We are currently experimenting with the HTC Vive, which allows a person to walk around a 15' by 15' physical space while wearing a head mounted display.
Precise head movements are matched in virtual reality with no noticeable latency.
Additionally, the HTC Vive controllers are motion tracked with millimeter precision and provide a valuable method for interacting with virtual objects.
We have enabled VR using the HTC Vive in the `UnityMol molecular visualization software <http://www.baaden.ibpc.fr/umol/>`_ (created by Baaden, et al.) and incorporated electrostatic surface data (see figure below and a `YouTube video <https://www.youtube.com/watch?v=Xxb3W8jnnp8&t=21s>`_).
New viewing capabilities now include walking around, grabbing (using the motion controllers), and scaling (gestures) of molecules.
We are actively working with Dr. Baaden and his group to determine the best use of interaction techniques for users to interact with molecular models through his software.

.. figure:: /media/1fas_VR.png

   View of UnityMol form the monitor as it is being used in VR with controllers.

For future work, we would like to further extend UnityMol in the HTC Vive to include natural user interactions for viewing multiple molecules, vary the electrostatic results from APBS, and change molecular attributes.
We envision this tool will also enable virtual collaboration for participant in different locations.
Each participant will be able to view, gesture and interact with the same data in the same VR space.
Finally, we would like to explore the use of VR for research related to docking of different molecules.

--------------------
Getting the software
--------------------

#. Download :file:`UnityMol-APBS-PS.zip` from `SourceForge <https://sourceforge.net/projects/unitymol-apbs/>`_.
#. Unzip :file:`UnityMol-APBS-PS.zip`; the resulting folder contains :file:`UnityMol-APBS.zip` and :file:`APBS-PDB2PQR.zip`.
#. :file:`Unzip UnityMol-APBS.zip`; the resulting folder contains :program:`UnityMol.exe`.
#. Optionally unzip :file:`APBS-PDB2PQR.zip` into :file:`C:\` to generate three directories :file:`apbs` (containing :program:`apbs` executable), :file:`pdb2pqr` (containing :program:`pdb2pqr` executable), and :file:`OutputFiles`. Alternatively, these executables can be downloaded and installed separately.

------------------
Using the software
------------------

Launch UnityMol.exe :program:`UnityMol.exe` to start the VR visualization.
The user interface is illustrated below.

.. figure:: /media/UI.png

   UnityMol-APBS user interface for PDB2PQR and APBS.
   (A) The main UnityMolAPBS menu; orange box highlights the two buttons used to open the APBS and PDB2PQR tools.
   (B) The main menu for interactions with APBS and PDB2PQR.
   Blue boxes show the buttons to launch PDB2PQR and APBS executables, green boxes show the location of the options used for producing the image in below, and the purple boxes highlight the two input fields required to use custom force fields and custom residue names.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Acetylcholinesterase example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The example illustrates the VR vizualization of the electrostatic surface potentials and electrostatic field lines of *Torpedo californica* acetylcholinesterase (AChE).

#. Download :file:`5ei5.pdb` from https://www.rcsb.org/structure/5EI5
#. Open UnityMol-APBS (VR or desktop)
#. Load :file:`5ei5.pdb` file
#. Open the :guilabel:`PDB2PQR panel`
#. Choose :guilabel:`options` (examples below) or run the default (default force field is AMBER)

   * :guilabel:`apbs-input` generates input file necessary for APBS
   * :guilabel:`drop-water` removes explicit water molecules from structure
   * :guilabel:`summary` writes atom names and sequence to a new file
   * :guilabel:`salt` writes salt bridge donor and acceptor atoms to a new file
   * :guilabel:`hbond` writes hydrogen bonding donors and acceptors to a new file. The resulting :file:`.hbond` and :file:`.salt` files can be loaded as a new selection in UnityMol-APBS

#. Select :guilabel:`all(5EI5)` and run PDB2PQR
#. :file:`5ei5X.pqr` is written to a file and is immediately loaded for the user.
#. Select :file:`all(5EI5)` and run APBS
#. :file:`5ei5X.dx` is created and loaded into the selection :guilabel:`all(5EI5X)` automatically
#. Select the :guilabel:`+` button on the :guilabel:`all(5EI5X)` selection tab, then select :guilabel:`surface`
#. Select :guilabel:`color by charge`
#. Select the :guilabel:`+` button on the :guilabel:`all(5EI5X)` selection tab, then select :guilabel:`field lines`

As a result of these steps, you should see a figure similar to the following.

.. figure:: /media/AChE.png

   Electrostatic surface potential and field lines of Torpedo californica AChE (PDB ID 5EI5) with bound alkylene-linked bis-tacrine.
   (A) Electrostatic field lines and protein secondary structure shown with alpha helices (yellow), beta sheets (blue), and random coils (white).
   Residues Tyr70, Trp84, Trp279, and Phe330 are shown interacting with alkylene-linked bis-tacrine via hydrogen bonding and π-π stacking interactions.
   The red oval highlights the potential gradient.
   (B) AChE surface model with field lines and mapped electrostatic surface potentials shown with neutral, negative, and positive charges in white, red, and blue, respectively.
   Field lines are calculated from a gradient (value of 0.2) and depicted with the starting points in red and the ending points in blue.
   The orientation is the same in Figures A and B, where the alkylene-linked bis-tacrine can be seen occupying the catalytic gorge.
   The white circle highlights the potential gradient exiting the catalytic gorge.
================
Binding energies
================

In general, implicit solvent models are used to calculation the contribution of solvation to binding free energies.
Additional binding free energy contributions (molecular mechanics energies, entropic changes, etc.) must be calculated separately and are not discussed in this tutorial.

-----------------
Free energy cycle
-----------------

Our framework for calculating solvation contributions to binding free energies is shown in the figure below:

.. image:: /media/apbs_bind_eng.png

This binding free energy cycle illustrates binding in terms of transfer free energies from a homogeneous dielectric environment (where interactions are described by Coulomb's law) to an inhomogeneous dielectric environment with differing internal (green) and external (cyan) dielectric constants.
The binding (dissociation) free energy is depicted in Step 3.
The binding free energy is given by

.. math::

   \Delta_b G = -\Delta_3 G =\Delta_4 G-\Delta_1 G-\Delta_2 G.

The following sections provide more detail on calculating individual terms of this equation.

---------------------------
Binding energy calculations
---------------------------

The most general method for calculating binding free energies divides the binding process up into solvation :math:`\Delta\Delta_s G` and Coulombic :math:`\Delta\Delta_c G` components:

.. math::

   \Delta\Delta_b G = \Delta\Delta_s G + \Delta\Delta_c G.

As mentioned above, this framework neglects the numerous other mechanical and entropic components actually involved in the binding process.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Solvation contribution to binding
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If we're just interested in calculating the solvation contributions to binding (steps 4 and 2 in the binding free energy cycle), then we simply need to follow the instructions from the :doc:`solvation-energies` section for the complex and isolated components.
The solvation energy contribution to the binding is then given by

.. math::

   \Delta\Delta_s G = \Delta_4 G - \Delta_2 G = \Delta_s G_{cmpx} - \Delta_s G_{mol1} - \Delta_s G_{mol2}

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Coulombic contribution to binding
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To complete the binding free energy cycle, we need to add intermolecular Coulombic contributions to the solvation energy change upon binding to get the total electrostatic/solvent contribution to the binding free energy.
In particular, we're interested in the change in Coulombic electrostatic energy upon binding, as given by

.. math::

   \Delta\Delta_c G = -\Delta_1 G =  \Delta_c G_{cmpx} - \Delta_c G_{mol1} - \Delta_c G_{mol2}

Each of the  quantities in this equation is the sum of pairwise Coulombic interactions between all atoms in the molecule (or complex) for a particular uniform dielectric.
In order to combine these Coulombic binding energies with the solvation energies described above, we need to make sure consistent dielectric constants are used.
In particular, Coulombic interactions should be calculated using the same uniform dielectric constant as the reference state of the solvation energy above.
For example, if solvation energies are calculated for transferring a protein from a homogeneous medium with uniform dielectric of  to an inhomogeneous medium with internal dielectric :math:`\epsilon_u` and external dielectric :math:`\epsilon_v`, then Coulombic energies should be calculated using a dielectric of :math:`\epsilon_u`.
The APBS accessory program :file:`tools/manip/coulomb` was created to help with the calculation of these analytic individual per-molecule Coulombic energies.
Given a PQR file as input, the :file:`tools/manip/coulomb` program calculates Coulombic energies for a vacuum dielectric (e.g., a uniform dielectric of 1).
If the reference dielectric is :math:`\epsilon_u`, then all energies returned by :file:`tools/manip/coulomb` need to be divided by :math:`\epsilon_u`.


^^^^^^^^^^^^^^
Other examples
^^^^^^^^^^^^^^

Several binding energy examples are distributed in the :file:`examples` directory with APBS.
                                                                                
  855
  26.7783817  34.2836311  35.8978118  27.3693644  34.3718208  35.0835203
  27.2500667  34.7693554  36.6473199  25.8929529  34.7253322  35.6951552
  26.6142677  32.8526458  36.2218570  26.9920881  32.6775049  37.2292078
  25.1527424  32.3891512  36.1524727  24.7082096  32.7360603  35.2196642
  25.1273377  31.2994731  36.1655793  24.3049146  32.8821039  37.3246880
  24.7634809  32.5370443  38.2513778  24.2917562  33.9720034  37.3321755
  22.5957266  32.2768716  37.2620219  22.1174790  32.5393270  38.9903550
  22.2734255  33.5822510  39.2663413  21.0657470  32.2850244  39.1221404
  22.7220461  31.9010765  39.6347610  27.4439657  32.0086624  35.2704124
  27.5767313  32.3497001  34.0964173  27.9571120  30.8827680  35.7619373
  27.8051102  30.6567495  36.7345223  28.6839266  29.8844547  34.9691881
  29.0578363  30.3619189  34.0635358  29.8994105  29.3512604  35.7301216
  29.5740909  28.8453618  36.6391591  30.7648198  28.4008434  34.9005667
  31.1260229  28.9082641  34.0060090  31.6170725  28.0709099  35.4947263
  30.1905548  27.5214691  34.6088628  30.7043910  30.4532250  36.0682029
  31.4355731  30.1164310  36.5912347  27.7498859  28.7377794  34.5993209
  27.1886688  28.0700506  35.4689926  27.5785547  28.5206774  33.3022718
  28.1483523  29.0754754  32.6796305  26.7763967  27.4582392  32.7046408
  26.1410471  26.9971418  33.4608498  25.8915034  28.0456919  31.5962632
  26.5296664  28.4561116  30.8137143  25.3067584  27.2410244  31.1504989
  24.9486447  29.1295637  32.0722153  25.4218136  30.4464892  32.2253573
  26.4436585  30.6890719  31.9734906  24.5791852  31.4459362  32.7372461
  24.9549762  32.4496740  32.8699215  23.2520761  31.1243481  33.0837521
  22.4451419  32.0546990  33.6451210  22.9129987  32.8874300  33.7417952
  22.7586866  29.8226497  32.8818972  21.7213371  29.6539273  33.1304436
  23.6129545  28.8167399  32.3874202  23.2682523  27.8021493  32.2528735
  27.6911236  26.3898119  32.1168863  28.8417785  26.6694184  31.7793292
  27.1535326  25.1903420  31.9121434  26.1952892  25.0716376  32.2084634
  27.8170870  24.0624248  31.2563012  28.8148670  24.3611533  30.9349421
  27.9758199  22.9389919  32.2921372  28.4881634  23.3629805  33.1557593
  26.9953008  22.5899979  32.6160613  28.8042951  21.7515886  31.7837119
  28.2333186  21.1916165  31.0430641  29.7134368  22.1227026  31.3106439
  29.2215063  20.8249301  32.9307475  30.1993662  20.4061486  32.6932447
  29.3408543  21.3885077  33.8561262  28.1969588  19.7062725  33.1344315
  27.2575070  20.1383101  33.4793692  28.0174021  19.2120764  32.1795820
  28.6848319  18.7076225  34.1037622  28.9182831  19.1344975  34.9888579
  27.9929875  17.9908405  34.2702658  29.5085221  18.2506535  33.7392034
  27.0438278  23.6415144  30.0092541  25.8225721  23.7714851  29.9554462
  27.7638201  23.1137727  29.0313784  28.7642451  23.0935584  29.1686470
  27.2704467  22.4755328  27.8237725  26.1824811  22.4248084  27.8663673
  27.6700300  23.3148588  26.6008180  27.2058510  24.2979661  26.6791340
  28.7520661  23.4461789  26.6089360  27.2602544  22.6673942  25.2659117
  27.6332051  21.6442155  25.2207666  25.7398141  22.6406576  25.1153293
  25.3469489  23.6445752  25.2760337  25.4683363  22.2958715  24.1175794
  25.3011825  21.9574839  25.8425625  27.8714595  23.4430001  24.1000711
  28.9578218  23.3940636  24.1751010  27.5624267  22.9945739  23.1559073
  27.5620972  24.4876504  24.1316152  27.8490726  21.0603937  27.7495624
  29.0667338  20.8813022  27.7268846  26.9684755  20.0664724  27.6921302
  25.9954081  20.3240297  27.7753424  27.2875488  18.6643002  27.4204995
  28.3578216  18.5099327  27.5579326  26.5561924  17.7230037  28.4032866
  25.4810903  17.8347480  28.2628510  26.9392315  16.2622244  28.1001374
  28.0066159  16.1146556  28.2647508  26.3785275  15.5814508  28.7406217
  26.6978401  16.0050245  27.0687692  26.9138722  18.0969299  29.8632860
  27.9938423  18.0323266  29.9962786  26.6165044  19.1283324  30.0524024
  26.2334344  17.2341981  30.9327726  26.6120257  16.2126232  30.8992501
  26.4487545  17.6455598  31.9189037  25.1544735  17.2340583  30.7776738
  26.9515814  18.3755811  25.9546024  25.7839270  18.3564389  25.5604791
  27.9818222  18.1965683  25.1316858  28.9019401  18.1862241  25.5481267
  27.8618442  17.8227596  23.7259626  26.9361644  18.2361255  23.3254241
  29.0302195  18.4025399  22.9036897  29.9504654  17.8850099  23.1746970
  28.8360958  18.1677707  21.8570904  29.2299717  19.9283646  23.0529108
  28.3826393  20.3635808  23.5826559  30.5109809  20.2205236  23.8416748
  31.3781987  19.8451003  23.2984439  30.6191043  21.2921254  24.0089698
  30.4643517  19.7350717  24.8164769  29.3125857  20.6151017  21.6871869
  28.3622450  20.4971256  21.1666563  29.5038367  21.6807265  21.8135071
  30.1045919  20.1697929  21.0851215  27.8107603  16.2932459  23.5986593
  28.7449998  15.6000959  24.0060211  26.7370058  15.7754164  23.0060920
  25.9792324  16.4079619  22.7919860  26.5981573  14.4086385  22.5106084
  27.5506493  13.8939105  22.6365331  25.5521132  13.6449574  23.3415088
  25.8271809  13.6783483  24.3956881  24.5772946  14.1164457  23.2168004
  25.4781141  12.1858053  22.9158545  26.4856113  11.5248220  22.7053064
  24.2972587  11.6348537  22.7630588  24.2845960  10.6887671  22.4098327
  23.4554023  12.1711918  22.9168205  26.3077790  14.4214015  20.9971890
  25.2614612  13.9703101  20.5231955  27.2505678  14.9576959  20.2237093
  28.1029317  15.2912275  20.6507619  27.2621845  14.8192424  18.7745442
  26.2507812  14.9286581  18.3832791  27.8861576  15.6018830  18.3430273
  27.8207553  13.4608655  18.3470875  28.4529737  12.7330358  19.1186967
  27.6235284  13.1284148  17.0760510  27.1266932  13.8145216  16.5259763
  28.1347666  11.9251925  16.4081836  27.8247129  11.0387237  16.9615781
  27.5222157  11.8690321  14.9998863  27.6648719  12.8354748  14.5163749
  28.0377696  11.1152089  14.4047792  26.0201465  11.5214585  15.0465290
  25.9078955  10.4373273  15.0324937  25.5772630  11.8868679  15.9730708
  25.2218247  12.1336235  13.8898168  24.1778663  11.8381534  13.9948727
  25.2796720  13.2201232  13.9557194  25.7408470  11.6756739  12.5268020
  26.7919058  11.9503199  12.4373079  25.6678718  10.5898645  12.4647156
  24.9779853  12.2923429  11.4243109  25.0379155  13.2998810  11.4618173
  25.3366572  11.9976803  10.5272925  24.0019806  12.0359825  11.4670255
  29.6650608  11.9146805  16.3766307  30.2783536  10.8924312  16.6765658
  30.2641788  13.0694956  16.0927276  29.6452134  13.8257516  15.8376539
  31.7161545  13.2867739  15.9707338  32.2203413  12.3389571  15.7824756
  31.9847284  14.2157372  14.7808998  31.4973288  15.1744924  14.9581529
  33.4621913  14.4642951  14.4833580  33.9667307  13.5192343  14.2823787
  33.5533680  15.1159781  13.6143471  33.9404291  14.9594480  15.3285074
  31.4208465  13.6202278  13.6411104  31.4809275  14.2665381  12.9337862
  32.3188698  13.9098940  17.2249468  33.3774288  13.4804613  17.6830716
  31.6569387  14.9336418  17.7706665  30.7628340  15.1855713  17.3741263
  32.1675654  15.7892353  18.8364549  33.2148496  15.5359931  19.0014132
  32.0789503  17.2618043  18.3891901  32.1695021  17.3030130  17.3037033
  31.0746638  17.6185124  18.6178753  33.1269047  18.2549047  18.9482320
  33.8371852  18.4738466  18.1509169  32.4251487  19.5543603  19.3339759
  31.8305161  19.3759559  20.2298402  33.1620118  20.3314818  19.5372420
  31.7780533  19.8830863  18.5208244  33.9319601  17.8195938  20.1776709
  34.5426960  16.9469669  19.9463477  34.6106515  18.6246567  20.4594130
  33.2716671  17.6116638  21.0196157  31.3827701  15.5272174  20.1250387
  30.1729304  15.7616851  20.1989707  32.0830373  15.0234601  21.1369837
  33.0748009  14.8954181  20.9949775  31.5619984  14.7435848  22.4722784
  30.5676260  15.1777692  22.5759704  31.4840820  13.2331415  22.7093613
  32.4700323  12.7972440  22.5478695  31.2002783  13.0552583  23.7466764
  30.4654235  12.5448934  21.7983565  29.4697317  12.9347174  22.0098555
  30.7156997  12.7394251  20.7554900  30.4954418  11.0356123  22.0566010
  31.4905717  10.6464863  21.8409601  30.2774006  10.8509652  23.1084896
  29.4629511  10.2909817  21.2121492  29.4294516   9.2529878  21.5433049
  28.4830635  10.7350225  21.3876178  29.7838121  10.3404761  19.7714578
  30.7171179  10.0027400  19.5841718  29.1294584   9.7873300  19.2366658
  29.7066587  11.2899194  19.4357225  32.4729240  15.3727935  23.5068725
  33.6964646  15.3518829  23.3576812  31.8793423  15.9099864  24.5593286
  30.8703131  15.9028461  24.6031147  32.6262338  16.5670659  25.6143291
  33.0375029  15.8110409  26.2832687  33.4496168  17.1344962  25.1805463
  31.7831313  17.5162352  26.4391642  30.5504636  17.4917071  26.4238971
  32.4980874  18.3248812  27.2030927  33.5036759  18.2951612  27.1135620
  31.9664733  19.1726510  28.2491310  30.9275898  19.4120689  28.0222667
  32.0157474  18.3763488  29.5615379  31.6182580  17.3770027  29.3840857
  33.0487239  18.2828488  29.8967681  31.1858289  19.0386844  30.6529369
  31.6891788  19.9508646  30.9734171  30.2171835  19.3200584  30.2398058
  30.9606673  18.1018300  31.8465786  30.0451323  17.2535931  31.7925125
  31.5942708  18.2829286  32.9098385  32.7799871  20.4694145  28.3069328
  33.9972138  20.4495630  28.5071743  32.1083439  21.5960702  28.1102427
  31.1161394  21.5366048  27.9308250  32.6431904  22.9526545  28.2565141
  33.6001718  22.9189948  28.7771242  32.8457511  23.6558722  26.9046236
  33.0398118  24.7137168  27.0819341  34.0340999  23.1042464  26.1273605
  33.8598641  22.0618861  25.8606599  34.1826522  23.6962393  25.2243154
  34.9301391  23.1786380  26.7435787  31.6912402  23.5451861  26.1076933
  31.8283047  24.1041733  25.3393503  31.7017773  23.7917037  29.1017831
  30.5663553  23.4083673  29.3847411  32.1698341  24.9645510  29.5028827
  33.1116798  25.2354624  29.2585456  31.4110796  25.9318236  30.2864859
  30.3654013  25.6258740  30.3187260  31.9543818  25.9989221  31.7240411
  31.3358065  26.6705549  32.3193061  32.0195324  24.6402111  32.4291275
  32.7353969  23.9867892  31.9304505  32.3401437  24.7765689  33.4619688
  31.0371395  24.1680748  32.4206851  33.2660496  26.5029300  31.6788060
  33.5631335  26.5869869  32.5878447  31.4729288  27.3071914  29.6240538
  32.2734071  27.5520812  28.7166100  30.6270866  28.2239402  30.0782979
  29.9955774  27.9473841  30.8164505  30.6918035  29.6466411  29.7474638
  31.7280011  29.9012295  29.5247099  29.8503056  29.9611255  28.5028443
  30.0352680  29.1955390  27.7492815  28.3417066  30.0582127  28.7512367
  28.1095067  30.9179118  29.3798423  27.8213988  30.1700112  27.7999434
  27.9904050  29.1508212  29.2424386  30.2850940  31.1955243  27.9929954
  29.7818491  31.3641721  27.1930015  30.2523922  30.4984303  30.9310851
  29.5540692  30.0182194  31.8225759  30.6097767  31.7763956  30.9060896
  31.0915264  32.0783077  30.0712758  30.1508765  32.8163309  31.8238775
  29.6533846  32.3506473  32.6746638  31.3642398  33.6150135  32.3362682
  32.2641000  33.0097481  32.2265584  31.5086516  34.5105853  31.7318660
  31.2595011  33.9997911  33.8176770  31.1893011  33.0833752  34.4036516
  32.1868762  34.4967527  34.1026204  30.0754700  34.9201042  34.1487614
  29.8164129  35.8817001  33.3910633  29.4008949  34.6970654  35.1821592
  29.1531247  33.7184371  31.0825746  29.4765000  34.2666848  30.0245491
  27.9321020  33.8489210  31.5984199  27.7300734  33.3796466  32.4696903
  26.9041010  34.7275324  31.0497433  27.3902433  35.5354475  30.5028154
  26.0470645  33.9243988  30.0635143  25.4874737  33.1619306  30.6053089
  25.3472227  34.5853697  29.5522089  26.6844214  33.4499553  29.3173278
  26.0351212  35.3704870  32.1433009  25.8974012  34.8493636  33.2541926
  25.3983169  36.4892487  31.7902239  25.6008504  36.8449196  30.8668592
  24.4027327  37.2055014  32.6075635  24.8592288  37.4561857  33.5650967
  23.9291587  38.5147507  31.9352265  22.9385038  38.7593510  32.3186121
  24.8527487  39.6716918  32.3188224  25.8702610  39.4783818  31.9791758
  24.4872533  40.5961130  31.8715780  24.8495360  39.7918856  33.4021857
  23.8278426  38.4437149  30.4012519  23.1832955  37.6177115  30.1006765
  23.3955755  39.3711735  30.0256106  24.8143720  38.3245811  29.9532194
  23.1813244  36.3481944  32.9307591  22.7335912  36.3188021  34.0758969
  22.6457547  35.6461497  31.9337825  23.0867543  35.6653796  31.0253902
  21.4829881  34.7721319  32.0518268  21.4623504  34.3683113  33.0640240
  20.1848454  35.5859583  31.8759389  19.3388465  34.9513764  32.1399523
  20.1978458  36.4093938  32.5900781  19.9326951  36.1578961  30.4769016
  20.4822783  35.6486768  29.4782607  19.1309595  37.1124957  30.3648392
  21.5718997  33.5708866  31.0931131  22.4860426  33.4620453  30.2656220
  20.6157519  32.6510364  31.2301121  19.9263624  32.7887054  31.9552986
  20.5357768  31.4449137  30.4203929  21.4878054  30.9185247  30.4894225
  19.4392560  30.5465485  31.0025623  18.4703034  31.0405832  30.9306450
  19.4028084  29.6075827  30.4501542  19.6515642  30.3316654  32.0498913
  20.2884786  31.7597635  28.9367394  20.8479426  31.0889348  28.0721868
  19.4972972  32.7916684  28.6227013  19.1344878  33.3603673  29.3743428
  19.2351248  33.1989221  27.2455448  18.7725501  32.3544877  26.7345244
  18.2311863  34.3591884  27.2341717  18.6450019  35.2325837  27.7380923
  17.9912816  34.6274385  26.2052432  17.3130859  34.0575583  27.7383944
  20.5363114  33.5267037  26.4890508  20.7433030  33.0340181  25.3825237
  21.4353357  34.3130130  27.0711869  21.1822367  34.7290427  27.9560852
  22.7326510  34.6629395  26.4827056  22.5594011  34.9985362  25.4602729
  23.3563532  35.8407700  27.2438721  23.5489161  35.5580123  28.2788093
  24.6356525  36.3676655  26.6004317  24.4219134  36.7024004  25.5853487
  25.0153460  37.2104058  27.1781721  25.3974014  35.5886174  26.5718629
  22.4486275  36.9189733  27.2100885  22.7539844  37.5855131  27.8298922
  23.6692867  33.4503807  26.4090201  24.2634042  33.2172410  25.3566633
  23.7589231  32.6219386  27.4590207  23.2361037  32.8442965  28.2941044
  24.5337575  31.3734899  27.4317771  25.5759055  31.6275793  27.2380259
  24.4622387  30.6893471  28.8029771  23.4338846  30.4136036  29.0364646
  25.0755911  29.7884290  28.7856472  24.8421421  31.3529153  29.5797616
  24.0865014  30.4172583  26.3057953  24.9242850  29.8444618  25.6104735
  22.7760666  30.2926783  26.0743679  22.1516928  30.7466962  26.7256505
  22.1794054  29.5284168  24.9765112  22.4447822  28.4790987  25.1057979
  20.6433515  29.6566264  25.0097798  20.2797578  29.5177348  26.0278845
  20.3734113  30.6635818  24.6915729  19.9555075  28.6150816  24.1062260
  20.6489807  28.2909250  23.3302273  19.7174788  27.7307116  24.6972699
  18.6810641  29.1049108  23.4081206  18.0845912  30.1332249  23.8109368
  18.2936658  28.4472048  22.4144393  22.7148865  29.9904427  23.6138361
  23.1295879  29.1615263  22.8081062  22.7334405  31.3074937  23.3579861
  22.4533975  31.9033746  24.1239323  23.2170930  31.9203864  22.1077517
  22.6673951  31.4996546  21.2657605  23.0199378  33.4546264  22.1269984
  23.0766534  33.8407597  23.1447717  23.8522052  33.9003379  21.5821273
  21.7451016  33.9799041  21.4410707  21.9537868  35.0002344  21.1192319
  21.5292363  33.3946135  20.5472189  20.4988755  34.0560990  22.3267498
  20.7595746  34.5251820  23.2755405  19.7630838  34.6911532  21.8333715
  19.8632759  32.6899439  22.5635407  19.5479955  32.2593612  21.6130677
  20.5980532  32.0243678  23.0164437  18.7141613  32.8107797  23.4746733
  17.9853101  33.3986823  23.0961079  18.3280319  31.8945315  23.6523408
  19.0334456  33.1789764  24.3593664  24.6880380  31.5693432  21.8663803
  25.0486360  31.0968966  20.7883252  25.5202019  31.7719388  22.8877164
  25.0966799  32.1601385  23.7184452  26.9553156  31.4536764  22.9171383
  27.4554084  32.0386326  22.1451652  27.5481395  31.8634995  24.2852052
  26.9608173  31.4108027  25.0840838  29.0053773  31.4231423  24.4609851
  29.6122873  31.8086372  23.6416893  29.3955602  31.8092883  25.4026070
  29.0717595  30.3353321  24.4814491  27.5097697  33.3898880  24.4641655
  26.4917826  33.7690155  24.3745698  27.8814547  33.6536552  25.4543144
  28.1294524  33.8736949  23.7091346  27.2008179  29.9734821  22.5890711
  27.9733690  29.6545163  21.6859140  26.5250364  29.0568492  23.2831520
  25.9142773  29.3850636  24.0175634  26.6467985  27.6105946  23.0921379
  27.7072994  27.3600890  23.0664656  25.9884440  26.8899207  24.2793155
  24.9963251  27.3148402  24.4318531  25.8365513  25.8446184  24.0103901
  26.7460393  26.9079020  25.6027382  28.1576502  26.8799063  25.6539097
  28.7624782  26.9009941  24.7594154  28.8229567  26.7819187  26.8875676
  29.9022661  26.7483812  26.9094291  28.0886216  26.6977784  28.0812684
  28.6007702  26.5989359  29.0269172  26.6854462  26.7413551  28.0418604
  26.1181842  26.6709592  28.9581818  26.0185886  26.8480455  26.8076580
  24.9389395  26.8625982  26.7827399  26.0698047  27.1039758  21.7594318
  26.6616857  26.2175799  21.1474825  24.9590031  27.6581502  21.2635766
  24.4834945  28.3313090  21.8474376  24.4177567  27.4055811  19.9178873
  24.1952910  26.3442003  19.8081553  23.1302704  28.2180693  19.7062457
  23.2643796  29.2146154  20.1268679  22.9583151  28.3521103  18.6382554
  21.8827024  27.5837077  20.3411078  22.1523438  27.0327642  21.2421558
  21.2144364  28.3944116  20.6314877  21.1315753  26.6507867  19.3774480
  21.2002108  27.0391234  18.3612565  21.5795130  25.6571540  19.3827011
  19.6409954  26.5816631  19.7208009  19.2268981  27.5863558  19.6354334
  19.1372711  25.9601632  18.9804494  19.3742901  26.0651925  21.0790448
  19.7586300  26.6803150  21.7819307  18.3791472  26.0677621  21.2517530
  19.7195047  25.1237565  21.2003157  25.4521179  27.7459631  18.8504628
  25.7605691  26.8964825  18.0183512  26.0253887  28.9492861  18.9065304
  25.6945394  29.5799436  19.6227267  27.1186907  29.3892383  18.0382037
  26.7605397  29.3598501  17.0091312  27.4728317  30.8501520  18.3653101
  26.5851728  31.4626177  18.2068042  27.7658557  30.9397621  19.4113357
  28.6051037  31.3833300  17.4730037  29.5357723  30.8752007  17.7254862
  28.3621521  31.1655703  16.4329699  28.8291368  32.8867128  17.6232957
  28.7126746  33.4659910  18.6947572  29.1736163  33.5783832  16.5590065
  29.3248190  34.5686341  16.6878706  29.2820997  33.1285057  15.6613188
  28.3300968  28.4447489  18.1246988  28.8295049  28.0042622  17.0921612
  28.7592307  28.0790102  19.3365623  28.3273629  28.5314428  20.1296012
  29.8190358  27.1001666  19.6022289  30.7532825  27.4936937  19.2015847
  29.9850052  26.9328438  21.1214725  30.2388192  27.8975338  21.5609815
  29.0413439  26.6107050  21.5616938  31.0503326  25.9276430  21.5031240
  32.4092752  26.2792827  21.4163383  32.6874725  27.2755446  21.1056504
  33.4053353  25.3335100  21.7224465  34.4478795  25.6059487  21.6490731
  33.0359256  24.0284952  22.1138940  33.9771329  23.0967959  22.4030765
  34.8596666  23.4638260  22.3130440  31.6761340  23.6783450  22.2059718
  31.4271195  22.6694773  22.5001719  30.6833707  24.6277498  21.9008867
  29.6399212  24.3548929  21.9581430  29.5622740  25.7452190  18.9182317
  30.4157177  25.2403564  18.1894369  28.3805132  25.1624742  19.1265163
  27.7473335  25.6261210  19.7623302  27.9124799  23.9367874  18.4865979
  28.5728397  23.1186673  18.7742319  26.5030913  23.6387149  19.0241062
  25.8278871  24.4619546  18.7909346  26.1019624  22.7292594  18.5769029
  26.5401977  23.5102254  20.1059073  27.9467485  24.0524662  16.9497438
  28.5462507  23.2172606  16.2709617  27.3493852  25.1207732  16.4134589
  26.9473262  25.7746506  17.0698999  27.2673958  25.4425333  14.9882631
  26.7360088  24.6339443  14.4862790  26.4519073  26.7439423  14.8059396
  26.6567998  27.4472238  15.6130348  26.7613270  27.2241773  13.8776576
  24.9488688  26.5073904  14.7047200  24.3282299  26.7691335  13.6850010
  24.3073002  25.9957502  15.7286879  23.3125514  25.8507789  15.6311735
  24.7924074  25.8700603  16.6056157  28.6619596  25.5124931  14.3316034
  28.8994930  24.8645725  13.3119936  29.6013378  26.2494110  14.9283979
  29.3254898  26.7581112  15.7562255  30.9835533  26.4150526  14.4554209
  30.9600029  26.7381002  13.4146076  31.6930583  27.5078050  15.2810207
  31.5537499  27.3009706  16.3421668  32.7623726  27.4456487  15.0790040
  31.2557446  28.9502593  14.9677768  30.4733954  29.1812479  14.0180015
  31.7538275  29.8897851  15.6323367  31.8048417  25.1099533  14.4931357
  32.7571227  24.9604055  13.7267229  31.4235495  24.1434043  15.3344577
  30.6378577  24.3541942  15.9331361  32.0357530  22.8152226  15.4472859
  33.0134144  22.8379099  14.9658656  32.2800430  22.5018723  16.9371498
  31.3568498  22.6168648  17.5051673  32.6065667  21.4656898  17.0248437
  33.3914128  23.3534583  17.5387607  34.5656194  23.1479558  17.2723298
  33.0842909  24.3191353  18.3702195  33.8291704  24.8787543  18.7600975
  32.1119620  24.5484642  18.5188646  31.2488596  21.7093339  14.6989267
  31.5588159  20.5261149  14.8457612  30.2519016  22.0663052  13.8782939
  30.0345846  23.0508488  13.8185021  29.5164101  21.1467913  12.9962097
  29.0637117  21.7239059  12.1898887  30.2163224  20.4437005  12.5445640
  28.4002689  20.3365874  13.6700173  27.8974027  19.3806942  13.0781368
  28.0227878  20.6898783  14.9000278  28.4535742  21.5188625  15.2839093
  26.9617058  20.0602306  15.6979922  26.9242914  19.0005109  15.4457035
  27.2747231  20.1764851  17.2052364  27.2601779  21.2266990  17.4966252
  26.2441686  19.4147865  18.0492691  26.2276456  18.3645202  17.7580238
  26.5047930  19.4865483  19.1051755  25.2510174  19.8427594  17.9130843
  28.6589194  19.6061106  17.5533662  29.4411295  20.1848510  17.0621171
  28.8246337  19.6621464  18.6292412  28.7268926  18.5663052  17.2335206
  25.5987881  20.6917800  15.3766625  25.4478535  21.9139560  15.3668161
  24.5954876  19.8531828  15.1259913  24.8424469  18.8740406  15.1047748
  23.2128944  20.1981080  14.7728453  22.9724944  21.1793342  15.1821929
  23.0989006  20.2734122  13.2382748  23.8949912  20.9148966  12.8602377
  23.2649197  19.2783662  12.8254381  21.7502707  20.8045816  12.7262094
  20.8757232  21.2242250  13.5154411  21.5072357  20.7910103  11.4986710
  22.2329728  19.1750021  15.3831058  21.8222221  18.2141005  14.7264303
  21.8950094  19.3540080  16.6639765  22.2372996  20.1838440  17.1268880
  21.1139648  18.4047106  17.4658468  20.7326672  17.6003848  16.8366598
  21.7743053  17.9588581  18.2096483  19.9227285  19.0243981  18.2016179
  19.6543970  20.2235887  18.1146596  19.2035478  18.1858993  18.9419595
  19.5511221  17.2395042  19.0022648  18.1675098  18.5773145  19.8988791
  17.5518217  19.3673762  19.4689546  17.2730260  17.3760722  20.2370197
  17.8950579  16.5205298  20.5002951  16.6505404  17.6269924  21.0959358
  16.3576810  17.0185347  19.0628948  15.7511372  17.8921550  18.8239862
  16.9621588  16.7780346  18.1882841  15.4471167  15.8282370  19.3833422
  15.9640535  14.7354292  19.7168296  14.2067076  15.9699221  19.2768095
  18.7993009  19.1087056  21.1847277  19.7553355  18.5235311  21.6998899
  18.2388571  20.1953935  21.7155747  17.4230587  20.5839755  21.2642946
  18.7024382  20.8451567  22.9386380  19.6798782  20.4435793  23.2057596
  18.8472314  22.3569655  22.7133457  17.9141510  22.7499429  22.3094630
  19.0085904  22.8401723  23.6769779  19.9712571  22.7613143  21.8190797
  19.9810732  22.6270680  20.4774611  19.1758450  22.2127611  19.8889295
  21.1972242  23.0380659  19.9846783  21.4466650  22.9713298  19.0082044
  22.0359046  23.4647306  20.9852539  23.3493531  23.9454061  20.9819058
  23.8950363  24.0089794  20.0521358  23.9326198  24.3291089  22.1994525
  24.9558925  24.6738907  22.2203288  23.1854297  24.2531708  23.3887839
  23.6203085  24.5765945  24.3228978  21.8763010  23.7381685  23.3797122
  21.3260218  23.6355257  24.3033278  21.2684244  23.3255132  22.1771490
  17.7494591  20.6107113  24.1064673  16.5320149  20.5056847  23.9366262
  18.3106608  20.6550727  25.3079926  19.3184346  20.7060646  25.3524001
  17.6010980  20.8536078  26.5728853  16.6008110  21.2367515  26.3708429
  17.4592701  19.5429785  27.3676710  16.7508709  19.7050524  28.1801236
  16.9677341  18.3555676  26.5362334  17.6828794  18.1193310  25.7482792
  16.8493170  17.4796042  27.1739606  16.0056752  18.5993747  26.0855812
  18.6948937  19.1864911  27.9311603  18.6559227  18.2495608  28.1367199
  18.3495129  21.8934848  27.4081013  19.5372973  22.1535277  27.1844903
  17.6696298  22.4703988  28.3960118  16.6962947  22.2297321  28.5179209
  18.2790013  23.2731477  29.4500419  19.3432071  23.0432020  29.5010076
  18.1342471  24.7699679  29.1475714  18.6687565  25.0002686  28.2259107
  17.0815895  25.0062087  28.9919500  18.6901697  25.6374059  30.2598335
  20.0665015  25.5924367  30.5467548  20.7225723  24.9880824  29.9379483
  20.5953345  26.3267399  31.6204315  21.6544195  26.2789272  31.8263570
  19.7362839  27.1147628  32.4103345  20.2499026  27.8035632  33.4560542
  21.2073257  27.7373487  33.4782914  18.3574937  27.1776568  32.1196540
  17.7122520  27.7946572  32.7274685  17.8338983  26.4282586  31.0496232
  16.7739791  26.4522765  30.8435751  17.6560236  22.9131158  30.8002902
  16.4323923  22.8464505  30.9339546  18.5058406  22.6838920  31.7956624
  19.4910090  22.7314138  31.5781346  18.1537635  22.3665844  33.1740469
  17.0679760  22.3711604  33.2702149  18.6270012  20.9456396  33.5127447
  18.2443640  20.2602331  32.7564466  19.7158161  20.9039077  33.4817995
  18.1353604  20.4641476  34.8792083  17.2908622  21.1449914  35.5051496
  18.5456549  19.3600961  35.3047843  18.7077735  23.4371725  34.1275169
  19.8701782  23.4108182  34.5369076  17.8676111  24.4152168  34.4696793
  16.9430332  24.3851166  34.0641684  18.2033792  25.5309205  35.3643892
  19.0355543  26.0847321  34.9298929  16.9937953  26.4692868  35.4526508
  16.6724368  26.7529310  34.4504774  16.1698852  25.9266513  35.9161749
  17.2893806  27.7334114  36.2672847  18.2006512  28.5109705  35.9004300
  16.5964894  27.9565757  37.2846131  18.6337014  25.0760963  36.7719531
  19.4311146  25.7528024  37.4226708  18.1640176  23.9072912  37.2262053
  17.5700457  23.3751713  36.6063466  18.5181834  23.3076475  38.5139725
  18.2848560  24.0260960  39.2997850  17.6464592  22.0647343  38.7307555
  17.8919464  21.3023242  37.9914462  17.8283384  21.6574503  39.7253374
  16.5918019  22.3254752  38.6420678  20.0164142  22.9610608  38.6358668
  20.5445255  22.8862465  39.7460756  20.7097783  22.7884329  37.5068687
  20.1782567  22.7809174  36.6480710  22.1445818  22.4681500  37.4013280
  22.6019034  22.5273089  38.3889591  22.3376237  21.0351625  36.8880207
  23.3964109  20.7802859  36.9345304  21.5515365  19.9838705  37.6772974
  20.4808978  20.1042838  37.5116825  21.8419246  18.9849134  37.3519292
  21.7622009  20.0849649  38.7419784  21.9057023  20.9740552  35.5518835
  21.9554506  20.0526574  35.2869979  22.9191768  23.4367364  36.4922794
  24.1441921  23.3465950  36.4006502  22.2222123  24.3858621  35.8456905
  21.2221192  24.3182180  35.9695759  22.7092980  25.3660345  34.8542084
  21.8307445  25.8866056  34.4732106  23.6279648  26.4053937  35.5327066
  24.5319501  25.9030737  35.8771256  23.9328132  27.1527522  34.8001355
  22.9979549  27.1334104  36.7309644  22.6658884  26.4113839  37.4770311
  23.7660586  27.7571823  37.1881939  21.8236707  28.0297484  36.3272281
  22.1725468  28.7757750  35.6131476  21.0307584  27.4284942  35.8823987
  21.2613238  28.7241405  37.5637411  20.7640600  27.9719916  38.1762827
  22.0791250  29.1553097  38.1411448  20.2971998  29.7764103  37.1947568
  19.5497451  29.3708681  36.6498240  19.8679187  30.1884892  38.0108657
  20.7251759  30.5026632  36.6384067  23.3425747  24.7231452  33.6098679
  24.2209946  25.3122313  32.9751217  22.8678651  23.5401390  33.2304345
  22.1085344  23.1564303  33.7748447  23.4204054  22.7401881  32.1313567
  24.3798558  23.1636453  31.8344101  23.7057208  21.3080911  32.5872897
  22.7903788  20.8542358  32.9672344  24.3014389  20.4139862  31.5037426
  25.2057056  20.8733712  31.1044824  24.5448665  19.4407755  31.9299922
  23.5806207  20.2641999  30.6999642  24.6811154  21.3599424  33.5960815
  24.5589567  20.5709873  34.1292542  22.4998094  22.7390090  30.9198515
  21.3210682  22.3915901  30.9974508  23.0741357  23.1021796  29.7817496
  24.0550775  23.3368910  29.8348809  22.5913832  22.8019406  28.4463600
  21.5019069  22.7711068  28.4369071  23.0896600  23.8902317  27.4857449
  24.1769572  23.9461632  27.5393012  22.8307729  23.5998791  26.4675237
  22.5173508  25.2642618  27.7609499  23.0891370  26.0893931  28.7492444
  23.9460687  25.7453063  29.3093728  22.5298467  27.3484587  29.0281281
  22.9593276  27.9804259  29.7913785  21.3808493  27.7703879  28.3390454
  20.9214165  28.7179015  28.5788108  20.8060379  26.9460999  27.3578723
  19.9055955  27.2613127  26.8516742  21.3810921  25.7007213  27.0578135
  20.9301100  25.0689747  26.3068142  23.1493540  21.4498213  28.0001055
  24.2866437  21.0896949  28.3148213  22.3946270  20.7325880  27.1830286
  21.4606050  21.0504472  26.9668787  22.9075499  19.6044421  26.4055074
  23.9969367  19.6290718  26.4307233  22.4665165  18.2415806  26.9817689
  23.1509067  17.4751740  26.6178998  22.4480974  18.1772522  28.5154409
  21.6936004  18.8557276  28.9136815  22.2099311  17.1657986  28.8445405
  23.4217949  18.4600024  28.9153873  21.1742360  17.9058495  26.5475087
  20.8531191  17.2125711  27.1287759  22.4746669  19.7865246  24.9522828
  21.4328599  20.3889977  24.6834820  23.2719518  19.2770603  24.0152717
  24.1275661  18.8425213  24.3303094  22.9386304  19.1443439  22.5905998
  21.8664595  19.3059156  22.4792533  23.6399974  20.2081217  21.7189391
  23.3171285  21.1879626  22.0707082  25.1734610  20.1696820  21.7907142
  25.5503458  19.2163829  21.4202247  25.5854271  20.9723596  21.1791015
  25.5012026  20.3109246  22.8205968  23.2166745  20.0763896  20.2477450
  22.1304535  20.1277610  20.1725473  23.6446912  20.8919190  19.6649028
  23.5615661  19.1304729  19.8301304  23.2553747  17.7200336  22.1389656
  24.3241283  17.1896226  22.4502643  22.3254917  17.1013593  21.4149647
  21.4628676  17.6047904  21.2645414  22.3595884  15.6816069  21.0442942
  23.3783989  15.3088538  21.1495950  21.4731595  14.8682543  21.9896199
  20.4503684  15.2432523  21.9515779  21.4870927  13.3780201  21.6633635
  22.5192046  13.0380237  21.5787706  20.9663839  12.8241189  22.4445548
  20.9811926  13.1967441  20.7150477  21.9711120  14.9899714  23.2983789
  21.2499425  14.7597657  23.8887563  21.9010108  15.4667311  19.6064694
  20.7588635  15.7711215  19.2561441  22.7903876  14.9310123  18.7729397
  23.6751363  14.6601329  19.1779650  22.5200181  14.5031920  17.3881706
  21.7582275  15.1501256  16.9532524  23.7939469  14.6313662  16.5411761
  24.5580685  13.9645913  16.9407798  23.5619524  14.3185900  15.5231128
  24.3429807  16.0636439  16.4980876  23.5413670  16.7325480  16.1848754
  24.6613347  16.3700898  17.4945540  25.5295043  16.1914246  15.5354721
  26.4943499  15.3935083  15.6374988  25.5021590  17.1081692  14.6812464
  21.9692006  13.0748325  17.3021193  22.4166723  12.1959000  18.0719668
  21.0838157  12.8419272  16.4522577
                                                                                
  855
  26.7783817  34.2836311  35.8978118  27.3693644  34.3718208  35.0835203
  27.2500667  34.7693554  36.6473199  25.8929529  34.7253322  35.6951552
  26.6142677  32.8526458  36.2218570  26.9920881  32.6775049  37.2292078
  25.1527424  32.3891512  36.1524727  24.7082096  32.7360603  35.2196642
  25.1273377  31.2994731  36.1655793  24.3049146  32.8821039  37.3246880
  24.7634809  32.5370443  38.2513778  24.2917562  33.9720034  37.3321755
  22.5957266  32.2768716  37.2620219  22.1174790  32.5393270  38.9903550
  22.2734255  33.5822510  39.2663413  21.0657470  32.2850244  39.1221404
  22.7220461  31.9010765  39.6347610  27.4439657  32.0086624  35.2704124
  27.5767313  32.3497001  34.0964173  27.9571120  30.8827680  35.7619373
  27.8051102  30.6567495  36.7345223  28.6839266  29.8844547  34.9691881
  29.0578363  30.3619189  34.0635358  29.8994105  29.3512604  35.7301216
  29.5740909  28.8453618  36.6391591  30.7648198  28.4008434  34.9005667
  31.1260229  28.9082641  34.0060090  31.6170725  28.0709099  35.4947263
  30.1905548  27.5214691  34.6088628  30.7043910  30.4532250  36.0682029
  31.4355731  30.1164310  36.5912347  27.7498859  28.7377794  34.5993209
  27.1886688  28.0700506  35.4689926  27.5785547  28.5206774  33.3022718
  28.1483523  29.0754754  32.6796305  26.7763967  27.4582392  32.7046408
  26.1410471  26.9971418  33.4608498  25.8915034  28.0456919  31.5962632
  26.5296664  28.4561116  30.8137143  25.3067584  27.2410244  31.1504989
  24.9486447  29.1295637  32.0722153  25.4218136  30.4464892  32.2253573
  26.4436585  30.6890719  31.9734906  24.5791852  31.4459362  32.7372461
  24.9549762  32.4496740  32.8699215  23.2520761  31.1243481  33.0837521
  22.4451419  32.0546990  33.6451210  22.9129987  32.8874300  33.7417952
  22.7586866  29.8226497  32.8818972  21.7213371  29.6539273  33.1304436
  23.6129545  28.8167399  32.3874202  23.2682523  27.8021493  32.2528735
  27.6911236  26.3898119  32.1168863  28.8417785  26.6694184  31.7793292
  27.1535326  25.1903420  31.9121434  26.1952892  25.0716376  32.2084634
  27.8170870  24.0624248  31.2563012  28.8148670  24.3611533  30.9349421
  27.9758199  22.9389919  32.2921372  28.4881634  23.3629805  33.1557593
  26.9953008  22.5899979  32.6160613  28.8042951  21.7515886  31.7837119
  28.2333186  21.1916165  31.0430641  29.7134368  22.1227026  31.3106439
  29.2215063  20.8249301  32.9307475  30.1993662  20.4061486  32.6932447
  29.3408543  21.3885077  33.8561262  28.1969588  19.7062725  33.1344315
  27.2575070  20.1383101  33.4793692  28.0174021  19.2120764  32.1795820
  28.6848319  18.7076225  34.1037622  28.9182831  19.1344975  34.9888579
  27.9929875  17.9908405  34.2702658  29.5085221  18.2506535  33.7392034
  27.0438278  23.6415144  30.0092541  25.8225721  23.7714851  29.9554462
  27.7638201  23.1137727  29.0313784  28.7642451  23.0935584  29.1686470
  27.2704467  22.4755328  27.8237725  26.1824811  22.4248084  27.8663673
  27.6700300  23.3148588  26.6008180  27.2058510  24.2979661  26.6791340
  28.7520661  23.4461789  26.6089360  27.2602544  22.6673942  25.2659117
  27.6332051  21.6442155  25.2207666  25.7398141  22.6406576  25.1153293
  25.3469489  23.6445752  25.2760337  25.4683363  22.2958715  24.1175794
  25.3011825  21.9574839  25.8425625  27.8714595  23.4430001  24.1000711
  28.9578218  23.3940636  24.1751010  27.5624267  22.9945739  23.1559073
  27.5620972  24.4876504  24.1316152  27.8490726  21.0603937  27.7495624
  29.0667338  20.8813022  27.7268846  26.9684755  20.0664724  27.6921302
  25.9954081  20.3240297  27.7753424  27.2875488  18.6643002  27.4204995
  28.3578216  18.5099327  27.5579326  26.5561924  17.7230037  28.4032866
  25.4810903  17.8347480  28.2628510  26.9392315  16.2622244  28.1001374
  28.0066159  16.1146556  28.2647508  26.3785275  15.5814508  28.7406217
  26.6978401  16.0050245  27.0687692  26.9138722  18.0969299  29.8632860
  27.9938423  18.0323266  29.9962786  26.6165044  19.1283324  30.0524024
  26.2334344  17.2341981  30.9327726  26.6120257  16.2126232  30.8992501
  26.4487545  17.6455598  31.9189037  25.1544735  17.2340583  30.7776738
  26.9515814  18.3755811  25.9546024  25.7839270  18.3564389  25.5604791
  27.9818222  18.1965683  25.1316858  28.9019401  18.1862241  25.5481267
  27.8618442  17.8227596  23.7259626  26.9361644  18.2361255  23.3254241
  29.0302195  18.4025399  22.9036897  29.9504654  17.8850099  23.1746970
  28.8360958  18.1677707  21.8570904  29.2299717  19.9283646  23.0529108
  28.3826393  20.3635808  23.5826559  30.5109809  20.2205236  23.8416748
  31.3781987  19.8451003  23.2984439  30.6191043  21.2921254  24.0089698
  30.4643517  19.7350717  24.8164769  29.3125857  20.6151017  21.6871869
  28.3622450  20.4971256  21.1666563  29.5038367  21.6807265  21.8135071
  30.1045919  20.1697929  21.0851215  27.8107603  16.2932459  23.5986593
  28.7449998  15.6000959  24.0060211  26.7370058  15.7754164  23.0060920
  25.9792324  16.4079619  22.7919860  26.5981573  14.4086385  22.5106084
  27.5506493  13.8939105  22.6365331  25.5521132  13.6449574  23.3415088
  25.8271809  13.6783483  24.3956881  24.5772946  14.1164457  23.2168004
  25.4781141  12.1858053  22.9158545  26.4856113  11.5248220  22.7053064
  24.2972587  11.6348537  22.7630588  24.2845960  10.6887671  22.4098327
  23.4554023  12.1711918  22.9168205  26.3077790  14.4214015  20.9971890
  25.2614612  13.9703101  20.5231955  27.2505678  14.9576959  20.2237093
  28.1029317  15.2912275  20.6507619  27.2621845  14.8192424  18.7745442
  26.2507812  14.9286581  18.3832791  27.8861576  15.6018830  18.3430273
  27.8207553  13.4608655  18.3470875  28.4529737  12.7330358  19.1186967
  27.6235284  13.1284148  17.0760510  27.1266932  13.8145216  16.5259763
  28.1347666  11.9251925  16.4081836  27.8247129  11.0387237  16.9615781
  27.5222157  11.8690321  14.9998863  27.6648719  12.8354748  14.5163749
  28.0377696  11.1152089  14.4047792  26.0201465  11.5214585  15.0465290
  25.9078955  10.4373273  15.0324937  25.5772630  11.8868679  15.9730708
  25.2218247  12.1336235  13.8898168  24.1778663  11.8381534  13.9948727
  25.2796720  13.2201232  13.9557194  25.7408470  11.6756739  12.5268020
  26.7919058  11.9503199  12.4373079  25.6678718  10.5898645  12.4647156
  24.9779853  12.2923429  11.4243109  25.0379155  13.2998810  11.4618173
  25.3366572  11.9976803  10.5272925  24.0019806  12.0359825  11.4670255
  29.6650608  11.9146805  16.3766307  30.2783536  10.8924312  16.6765658
  30.2641788  13.0694956  16.0927276  29.6452134  13.8257516  15.8376539
  31.7161545  13.2867739  15.9707338  32.2203413  12.3389571  15.7824756
  31.9847284  14.2157372  14.7808998  31.4973288  15.1744924  14.9581529
  33.4621913  14.4642951  14.4833580  33.9667307  13.5192343  14.2823787
  33.5533680  15.1159781  13.6143471  33.9404291  14.9594480  15.3285074
  31.4208465  13.6202278  13.6411104  31.4809275  14.2665381  12.9337862
  32.3188698  13.9098940  17.2249468  33.3774288  13.4804613  17.6830716
  31.6569387  14.9336418  17.7706665  30.7628340  15.1855713  17.3741263
  32.1675654  15.7892353  18.8364549  33.2148496  15.5359931  19.0014132
  32.0789503  17.2618043  18.3891901  32.1695021  17.3030130  17.3037033
  31.0746638  17.6185124  18.6178753  33.1269047  18.2549047  18.9482320
  33.8371852  18.4738466  18.1509169  32.4251487  19.5543603  19.3339759
  31.8305161  19.3759559  20.2298402  33.1620118  20.3314818  19.5372420
  31.7780533  19.8830863  18.5208244  33.9319601  17.8195938  20.1776709
  34.5426960  16.9469669  19.9463477  34.6106515  18.6246567  20.4594130
  33.2716671  17.6116638  21.0196157  31.3827701  15.5272174  20.1250387
  30.1729304  15.7616851  20.1989707  32.0830373  15.0234601  21.1369837
  33.0748009  14.8954181  20.9949775  31.5619984  14.7435848  22.4722784
  30.5676260  15.1777692  22.5759704  31.4840820  13.2331415  22.7093613
  32.4700323  12.7972440  22.5478695  31.2002783  13.0552583  23.7466764
  30.4654235  12.5448934  21.7983565  29.4697317  12.9347174  22.0098555
  30.7156997  12.7394251  20.7554900  30.4954418  11.0356123  22.0566010
  31.4905717  10.6464863  21.8409601  30.2774006  10.8509652  23.1084896
  29.4629511  10.2909817  21.2121492  29.4294516   9.2529878  21.5433049
  28.4830635  10.7350225  21.3876178  29.7838121  10.3404761  19.7714578
  30.7171179  10.0027400  19.5841718  29.1294584   9.7873300  19.2366658
  29.7066587  11.2899194  19.4357225  32.4729240  15.3727935  23.5068725
  33.6964646  15.3518829  23.3576812  31.8793423  15.9099864  24.5593286
  30.8703131  15.9028461  24.6031147  32.6262338  16.5670659  25.6143291
  33.0375029  15.8110409  26.2832687  33.4496168  17.1344962  25.1805463
  31.7831313  17.5162352  26.4391642  30.5504636  17.4917071  26.4238971
  32.4980874  18.3248812  27.2030927  33.5036759  18.2951612  27.1135620
  31.9664733  19.1726510  28.2491310  30.9275898  19.4120689  28.0222667
  32.0157474  18.3763488  29.5615379  31.6182580  17.3770027  29.3840857
  33.0487239  18.2828488  29.8967681  31.1858289  19.0386844  30.6529369
  31.6891788  19.9508646  30.9734171  30.2171835  19.3200584  30.2398058
  30.9606673  18.1018300  31.8465786  30.0451323  17.2535931  31.7925125
  31.5942708  18.2829286  32.9098385  32.7799871  20.4694145  28.3069328
  33.9972138  20.4495630  28.5071743  32.1083439  21.5960702  28.1102427
  31.1161394  21.5366048  27.9308250  32.6431904  22.9526545  28.2565141
  33.6001718  22.9189948  28.7771242  32.8457511  23.6558722  26.9046236
  33.0398118  24.7137168  27.0819341  34.0340999  23.1042464  26.1273605
  33.8598641  22.0618861  25.8606599  34.1826522  23.6962393  25.2243154
  34.9301391  23.1786380  26.7435787  31.6912402  23.5451861  26.1076933
  31.8283047  24.1041733  25.3393503  31.7017773  23.7917037  29.1017831
  30.5663553  23.4083673  29.3847411  32.1698341  24.9645510  29.5028827
  33.1116798  25.2354624  29.2585456  31.4110796  25.9318236  30.2864859
  30.3654013  25.6258740  30.3187260  31.9543818  25.9989221  31.7240411
  31.3358065  26.6705549  32.3193061  32.0195324  24.6402111  32.4291275
  32.7353969  23.9867892  31.9304505  32.3401437  24.7765689  33.4619688
  31.0371395  24.1680748  32.4206851  33.2660496  26.5029300  31.6788060
  33.5631335  26.5869869  32.5878447  31.4729288  27.3071914  29.6240538
  32.2734071  27.5520812  28.7166100  30.6270866  28.2239402  30.0782979
  29.9955774  27.9473841  30.8164505  30.6918035  29.6466411  29.7474638
  31.7280011  29.9012295  29.5247099  29.8503056  29.9611255  28.5028443
  30.0352680  29.1955390  27.7492815  28.3417066  30.0582127  28.7512367
  28.1095067  30.9179118  29.3798423  27.8213988  30.1700112  27.7999434
  27.9904050  29.1508212  29.2424386  30.2850940  31.1955243  27.9929954
  29.7818491  31.3641721  27.1930015  30.2523922  30.4984303  30.9310851
  29.5540692  30.0182194  31.8225759  30.6097767  31.7763956  30.9060896
  31.0915264  32.0783077  30.0712758  30.1508765  32.8163309  31.8238775
  29.6533846  32.3506473  32.6746638  31.3642398  33.6150135  32.3362682
  32.2641000  33.0097481  32.2265584  31.5086516  34.5105853  31.7318660
  31.2595011  33.9997911  33.8176770  31.1893011  33.0833752  34.4036516
  32.1868762  34.4967527  34.1026204  30.0754700  34.9201042  34.1487614
  29.8164129  35.8817001  33.3910633  29.4008949  34.6970654  35.1821592
  29.1531247  33.7184371  31.0825746  29.4765000  34.2666848  30.0245491
  27.9321020  33.8489210  31.5984199  27.7300734  33.3796466  32.4696903
  26.9041010  34.7275324  31.0497433  27.3902433  35.5354475  30.5028154
  26.0470645  33.9243988  30.0635143  25.4874737  33.1619306  30.6053089
  25.3472227  34.5853697  29.5522089  26.6844214  33.4499553  29.3173278
  26.0351212  35.3704870  32.1433009  25.8974012  34.8493636  33.2541926
  25.3983169  36.4892487  31.7902239  25.6008504  36.8449196  30.8668592
  24.4027327  37.2055014  32.6075635  24.8592288  37.4561857  33.5650967
  23.9291587  38.5147507  31.9352265  22.9385038  38.7593510  32.3186121
  24.8527487  39.6716918  32.3188224  25.8702610  39.4783818  31.9791758
  24.4872533  40.5961130  31.8715780  24.8495360  39.7918856  33.4021857
  23.8278426  38.4437149  30.4012519  23.1832955  37.6177115  30.1006765
  23.3955755  39.3711735  30.0256106  24.8143720  38.3245811  29.9532194
  23.1813244  36.3481944  32.9307591  22.7335912  36.3188021  34.0758969
  22.6457547  35.6461497  31.9337825  23.0867543  35.6653796  31.0253902
  21.4829881  34.7721319  32.0518268  21.4623504  34.3683113  33.0640240
  20.1848454  35.5859583  31.8759389  19.3388465  34.9513764  32.1399523
  20.1978458  36.4093938  32.5900781  19.9326951  36.1578961  30.4769016
  20.4822783  35.6486768  29.4782607  19.1309595  37.1124957  30.3648392
  21.5718997  33.5708866  31.0931131  22.4860426  33.4620453  30.2656220
  20.6157519  32.6510364  31.2301121  19.9263624  32.7887054  31.9552986
  20.5357768  31.4449137  30.4203929  21.4878054  30.9185247  30.4894225
  19.4392560  30.5465485  31.0025623  18.4703034  31.0405832  30.9306450
  19.4028084  29.6075827  30.4501542  19.6515642  30.3316654  32.0498913
  20.2884786  31.7597635  28.9367394  20.8479426  31.0889348  28.0721868
  19.4972972  32.7916684  28.6227013  19.1344878  33.3603673  29.3743428
  19.2351248  33.1989221  27.2455448  18.7725501  32.3544877  26.7345244
  18.2311863  34.3591884  27.2341717  18.6450019  35.2325837  27.7380923
  17.9912816  34.6274385  26.2052432  17.3130859  34.0575583  27.7383944
  20.5363114  33.5267037  26.4890508  20.7433030  33.0340181  25.3825237
  21.4353357  34.3130130  27.0711869  21.1822367  34.7290427  27.9560852
  22.7326510  34.6629395  26.4827056  22.5594011  34.9985362  25.4602729
  23.3563532  35.8407700  27.2438721  23.5489161  35.5580123  28.2788093
  24.6356525  36.3676655  26.6004317  24.4219134  36.7024004  25.5853487
  25.0153460  37.2104058  27.1781721  25.3974014  35.5886174  26.5718629
  22.4486275  36.9189733  27.2100885  22.7539844  37.5855131  27.8298922
  23.6692867  33.4503807  26.4090201  24.2634042  33.2172410  25.3566633
  23.7589231  32.6219386  27.4590207  23.2361037  32.8442965  28.2941044
  24.5337575  31.3734899  27.4317771  25.5759055  31.6275793  27.2380259
  24.4622387  30.6893471  28.8029771  23.4338846  30.4136036  29.0364646
  25.0755911  29.7884290  28.7856472  24.8421421  31.3529153  29.5797616
  24.0865014  30.4172583  26.3057953  24.9242850  29.8444618  25.6104735
  22.7760666  30.2926783  26.0743679  22.1516928  30.7466962  26.7256505
  22.1794054  29.5284168  24.9765112  22.4447822  28.4790987  25.1057979
  20.6433515  29.6566264  25.0097798  20.2797578  29.5177348  26.0278845
  20.3734113  30.6635818  24.6915729  19.9555075  28.6150816  24.1062260
  20.6489807  28.2909250  23.3302273  19.7174788  27.7307116  24.6972699
  18.6810641  29.1049108  23.4081206  18.0845912  30.1332249  23.8109368
  18.2936658  28.4472048  22.4144393  22.7148865  29.9904427  23.6138361
  23.1295879  29.1615263  22.8081062  22.7334405  31.3074937  23.3579861
  22.4533975  31.9033746  24.1239323  23.2170930  31.9203864  22.1077517
  22.6673951  31.4996546  21.2657605  23.0199378  33.4546264  22.1269984
  23.0766534  33.8407597  23.1447717  23.8522052  33.9003379  21.5821273
  21.7451016  33.9799041  21.4410707  21.9537868  35.0002344  21.1192319
  21.5292363  33.3946135  20.5472189  20.4988755  34.0560990  22.3267498
  20.7595746  34.5251820  23.2755405  19.7630838  34.6911532  21.8333715
  19.8632759  32.6899439  22.5635407  19.5479955  32.2593612  21.6130677
  20.5980532  32.0243678  23.0164437  18.7141613  32.8107797  23.4746733
  17.9853101  33.3986823  23.0961079  18.3280319  31.8945315  23.6523408
  19.0334456  33.1789764  24.3593664  24.6880380  31.5693432  21.8663803
  25.0486360  31.0968966  20.7883252  25.5202019  31.7719388  22.8877164
  25.0966799  32.1601385  23.7184452  26.9553156  31.4536764  22.9171383
  27.4554084  32.0386326  22.1451652  27.5481395  31.8634995  24.2852052
  26.9608173  31.4108027  25.0840838  29.0053773  31.4231423  24.4609851
  29.6122873  31.8086372  23.6416893  29.3955602  31.8092883  25.4026070
  29.0717595  30.3353321  24.4814491  27.5097697  33.3898880  24.4641655
  26.4917826  33.7690155  24.3745698  27.8814547  33.6536552  25.4543144
  28.1294524  33.8736949  23.7091346  27.2008179  29.9734821  22.5890711
  27.9733690  29.6545163  21.6859140  26.5250364  29.0568492  23.2831520
  25.9142773  29.3850636  24.0175634  26.6467985  27.6105946  23.0921379
  27.7072994  27.3600890  23.0664656  25.9884440  26.8899207  24.2793155
  24.9963251  27.3148402  24.4318531  25.8365513  25.8446184  24.0103901
  26.7460393  26.9079020  25.6027382  28.1576502  26.8799063  25.6539097
  28.7624782  26.9009941  24.7594154  28.8229567  26.7819187  26.8875676
  29.9022661  26.7483812  26.9094291  28.0886216  26.6977784  28.0812684
  28.6007702  26.5989359  29.0269172  26.6854462  26.7413551  28.0418604
  26.1181842  26.6709592  28.9581818  26.0185886  26.8480455  26.8076580
  24.9389395  26.8625982  26.7827399  26.0698047  27.1039758  21.7594318
  26.6616857  26.2175799  21.1474825  24.9590031  27.6581502  21.2635766
  24.4834945  28.3313090  21.8474376  24.4177567  27.4055811  19.9178873
  24.1952910  26.3442003  19.8081553  23.1302704  28.2180693  19.7062457
  23.2643796  29.2146154  20.1268679  22.9583151  28.3521103  18.6382554
  21.8827024  27.5837077  20.3411078  22.1523438  27.0327642  21.2421558
  21.2144364  28.3944116  20.6314877  21.1315753  26.6507867  19.3774480
  21.2002108  27.0391234  18.3612565  21.5795130  25.6571540  19.3827011
  19.6409954  26.5816631  19.7208009  19.2268981  27.5863558  19.6354334
  19.1372711  25.9601632  18.9804494  19.3742901  26.0651925  21.0790448
  19.7586300  26.6803150  21.7819307  18.3791472  26.0677621  21.2517530
  19.7195047  25.1237565  21.2003157  25.4521179  27.7459631  18.8504628
  25.7605691  26.8964825  18.0183512  26.0253887  28.9492861  18.9065304
  25.6945394  29.5799436  19.6227267  27.1186907  29.3892383  18.0382037
  26.7605397  29.3598501  17.0091312  27.4728317  30.8501520  18.3653101
  26.5851728  31.4626177  18.2068042  27.7658557  30.9397621  19.4113357
  28.6051037  31.3833300  17.4730037  29.5357723  30.8752007  17.7254862
  28.3621521  31.1655703  16.4329699  28.8291368  32.8867128  17.6232957
  28.7126746  33.4659910  18.6947572  29.1736163  33.5783832  16.5590065
  29.3248190  34.5686341  16.6878706  29.2820997  33.1285057  15.6613188
  28.3300968  28.4447489  18.1246988  28.8295049  28.0042622  17.0921612
  28.7592307  28.0790102  19.3365623  28.3273629  28.5314428  20.1296012
  29.8190358  27.1001666  19.6022289  30.7532825  27.4936937  19.2015847
  29.9850052  26.9328438  21.1214725  30.2388192  27.8975338  21.5609815
  29.0413439  26.6107050  21.5616938  31.0503326  25.9276430  21.5031240
  32.4092752  26.2792827  21.4163383  32.6874725  27.2755446  21.1056504
  33.4053353  25.3335100  21.7224465  34.4478795  25.6059487  21.6490731
  33.0359256  24.0284952  22.1138940  33.9771329  23.0967959  22.4030765
  34.8596666  23.4638260  22.3130440  31.6761340  23.6783450  22.2059718
  31.4271195  22.6694773  22.5001719  30.6833707  24.6277498  21.9008867
  29.6399212  24.3548929  21.9581430  29.5622740  25.7452190  18.9182317
  30.4157177  25.2403564  18.1894369  28.3805132  25.1624742  19.1265163
  27.7473335  25.6261210  19.7623302  27.9124799  23.9367874  18.4865979
  28.5728397  23.1186673  18.7742319  26.5030913  23.6387149  19.0241062
  25.8278871  24.4619546  18.7909346  26.1019624  22.7292594  18.5769029
  26.5401977  23.5102254  20.1059073  27.9467485  24.0524662  16.9497438
  28.5462507  23.2172606  16.2709617  27.3493852  25.1207732  16.4134589
  26.9473262  25.7746506  17.0698999  27.2673958  25.4425333  14.9882631
  26.7360088  24.6339443  14.4862790  26.4519073  26.7439423  14.8059396
  26.6567998  27.4472238  15.6130348  26.7613270  27.2241773  13.8776576
  24.9488688  26.5073904  14.7047200  24.3282299  26.7691335  13.6850010
  24.3073002  25.9957502  15.7286879  23.3125514  25.8507789  15.6311735
  24.7924074  25.8700603  16.6056157  28.6619596  25.5124931  14.3316034
  28.8994930  24.8645725  13.3119936  29.6013378  26.2494110  14.9283979
  29.3254898  26.7581112  15.7562255  30.9835533  26.4150526  14.4554209
  30.9600029  26.7381002  13.4146076  31.6930583  27.5078050  15.2810207
  31.5537499  27.3009706  16.3421668  32.7623726  27.4456487  15.0790040
  31.2557446  28.9502593  14.9677768  30.4733954  29.1812479  14.0180015
  31.7538275  29.8897851  15.6323367  31.8048417  25.1099533  14.4931357
  32.7571227  24.9604055  13.7267229  31.4235495  24.1434043  15.3344577
  30.6378577  24.3541942  15.9331361  32.0357530  22.8152226  15.4472859
  33.0134144  22.8379099  14.9658656  32.2800430  22.5018723  16.9371498
  31.3568498  22.6168648  17.5051673  32.6065667  21.4656898  17.0248437
  33.3914128  23.3534583  17.5387607  34.5656194  23.1479558  17.2723298
  33.0842909  24.3191353  18.3702195  33.8291704  24.8787543  18.7600975
  32.1119620  24.5484642  18.5188646  31.2488596  21.7093339  14.6989267
  31.5588159  20.5261149  14.8457612  30.2519016  22.0663052  13.8782939
  30.0345846  23.0508488  13.8185021  29.5164101  21.1467913  12.9962097
  29.0637117  21.7239059  12.1898887  30.2163224  20.4437005  12.5445640
  28.4002689  20.3365874  13.6700173  27.8974027  19.3806942  13.0781368
  28.0227878  20.6898783  14.9000278  28.4535742  21.5188625  15.2839093
  26.9617058  20.0602306  15.6979922  26.9242914  19.0005109  15.4457035
  27.2747231  20.1764851  17.2052364  27.2601779  21.2266990  17.4966252
  26.2441686  19.4147865  18.0492691  26.2276456  18.3645202  17.7580238
  26.5047930  19.4865483  19.1051755  25.2510174  19.8427594  17.9130843
  28.6589194  19.6061106  17.5533662  29.4411295  20.1848510  17.0621171
  28.8246337  19.6621464  18.6292412  28.7268926  18.5663052  17.2335206
  25.5987881  20.6917800  15.3766625  25.4478535  21.9139560  15.3668161
  24.5954876  19.8531828  15.1259913  24.8424469  18.8740406  15.1047748
  23.2128944  20.1981080  14.7728453  22.9724944  21.1793342  15.1821929
  23.0989006  20.2734122  13.2382748  23.8949912  20.9148966  12.8602377
  23.2649197  19.2783662  12.8254381  21.7502707  20.8045816  12.7262094
  20.8757232  21.2242250  13.5154411  21.5072357  20.7910103  11.4986710
  22.2329728  19.1750021  15.3831058  21.8222221  18.2141005  14.7264303
  21.8950094  19.3540080  16.6639765  22.2372996  20.1838440  17.1268880
  21.1139648  18.4047106  17.4658468  20.7326672  17.6003848  16.8366598
  21.7743053  17.9588581  18.2096483  19.9227285  19.0243981  18.2016179
  19.6543970  20.2235887  18.1146596  19.2035478  18.1858993  18.9419595
  19.5511221  17.2395042  19.0022648  18.1675098  18.5773145  19.8988791
  17.5518217  19.3673762  19.4689546  17.2730260  17.3760722  20.2370197
  17.8950579  16.5205298  20.5002951  16.6505404  17.6269924  21.0959358
  16.3576810  17.0185347  19.0628948  15.7511372  17.8921550  18.8239862
  16.9621588  16.7780346  18.1882841  15.4471167  15.8282370  19.3833422
  15.9640535  14.7354292  19.7168296  14.2067076  15.9699221  19.2768095
  18.7993009  19.1087056  21.1847277  19.7553355  18.5235311  21.6998899
  18.2388571  20.1953935  21.7155747  17.4230587  20.5839755  21.2642946
  18.7024382  20.8451567  22.9386380  19.6798782  20.4435793  23.2057596
  18.8472314  22.3569655  22.7133457  17.9141510  22.7499429  22.3094630
  19.0085904  22.8401723  23.6769779  19.9712571  22.7613143  21.8190797
  19.9810732  22.6270680  20.4774611  19.1758450  22.2127611  19.8889295
  21.1972242  23.0380659  19.9846783  21.4466650  22.9713298  19.0082044
  22.0359046  23.4647306  20.9852539  23.3493531  23.9454061  20.9819058
  23.8950363  24.0089794  20.0521358  23.9326198  24.3291089  22.1994525
  24.9558925  24.6738907  22.2203288  23.1854297  24.2531708  23.3887839
  23.6203085  24.5765945  24.3228978  21.8763010  23.7381685  23.3797122
  21.3260218  23.6355257  24.3033278  21.2684244  23.3255132  22.1771490
  17.7494591  20.6107113  24.1064673  16.5320149  20.5056847  23.9366262
  18.3106608  20.6550727  25.3079926  19.3184346  20.7060646  25.3524001
  17.6010980  20.8536078  26.5728853  16.6008110  21.2367515  26.3708429
  17.4592701  19.5429785  27.3676710  16.7508709  19.7050524  28.1801236
  16.9677341  18.3555676  26.5362334  17.6828794  18.1193310  25.7482792
  16.8493170  17.4796042  27.1739606  16.0056752  18.5993747  26.0855812
  18.6948937  19.1864911  27.9311603  18.6559227  18.2495608  28.1367199
  18.3495129  21.8934848  27.4081013  19.5372973  22.1535277  27.1844903
  17.6696298  22.4703988  28.3960118  16.6962947  22.2297321  28.5179209
  18.2790013  23.2731477  29.4500419  19.3432071  23.0432020  29.5010076
  18.1342471  24.7699679  29.1475714  18.6687565  25.0002686  28.2259107
  17.0815895  25.0062087  28.9919500  18.6901697  25.6374059  30.2598335
  20.0665015  25.5924367  30.5467548  20.7225723  24.9880824  29.9379483
  20.5953345  26.3267399  31.6204315  21.6544195  26.2789272  31.8263570
  19.7362839  27.1147628  32.4103345  20.2499026  27.8035632  33.4560542
  21.2073257  27.7373487  33.4782914  18.3574937  27.1776568  32.1196540
  17.7122520  27.7946572  32.7274685  17.8338983  26.4282586  31.0496232
  16.7739791  26.4522765  30.8435751  17.6560236  22.9131158  30.8002902
  16.4323923  22.8464505  30.9339546  18.5058406  22.6838920  31.7956624
  19.4910090  22.7314138  31.5781346  18.1537635  22.3665844  33.1740469
  17.0679760  22.3711604  33.2702149  18.6270012  20.9456396  33.5127447
  18.2443640  20.2602331  32.7564466  19.7158161  20.9039077  33.4817995
  18.1353604  20.4641476  34.8792083  17.2908622  21.1449914  35.5051496
  18.5456549  19.3600961  35.3047843  18.7077735  23.4371725  34.1275169
  19.8701782  23.4108182  34.5369076  17.8676111  24.4152168  34.4696793
  16.9430332  24.3851166  34.0641684  18.2033792  25.5309205  35.3643892
  19.0355543  26.0847321  34.9298929  16.9937953  26.4692868  35.4526508
  16.6724368  26.7529310  34.4504774  16.1698852  25.9266513  35.9161749
  17.2893806  27.7334114  36.2672847  18.2006512  28.5109705  35.9004300
  16.5964894  27.9565757  37.2846131  18.6337014  25.0760963  36.7719531
  19.4311146  25.7528024  37.4226708  18.1640176  23.9072912  37.2262053
  17.5700457  23.3751713  36.6063466  18.5181834  23.3076475  38.5139725
  18.2848560  24.0260960  39.2997850  17.6464592  22.0647343  38.7307555
  17.8919464  21.3023242  37.9914462  17.8283384  21.6574503  39.7253374
  16.5918019  22.3254752  38.6420678  20.0164142  22.9610608  38.6358668
  20.5445255  22.8862465  39.7460756  20.7097783  22.7884329  37.5068687
  20.1782567  22.7809174  36.6480710  22.1445818  22.4681500  37.4013280
  22.6019034  22.5273089  38.3889591  22.3376237  21.0351625  36.8880207
  23.3964109  20.7802859  36.9345304  21.5515365  19.9838705  37.6772974
  20.4808978  20.1042838  37.5116825  21.8419246  18.9849134  37.3519292
  21.7622009  20.0849649  38.7419784  21.9057023  20.9740552  35.5518835
  21.9554506  20.0526574  35.2869979  22.9191768  23.4367364  36.4922794
  24.1441921  23.3465950  36.4006502  22.2222123  24.3858621  35.8456905
  21.2221192  24.3182180  35.9695759  22.7092980  25.3660345  34.8542084
  21.8307445  25.8866056  34.4732106  23.6279648  26.4053937  35.5327066
  24.5319501  25.9030737  35.8771256  23.9328132  27.1527522  34.8001355
  22.9979549  27.1334104  36.7309644  22.6658884  26.4113839  37.4770311
  23.7660586  27.7571823  37.1881939  21.8236707  28.0297484  36.3272281
  22.1725468  28.7757750  35.6131476  21.0307584  27.4284942  35.8823987
  21.2613238  28.7241405  37.5637411  20.7640600  27.9719916  38.1762827
  22.0791250  29.1553097  38.1411448  20.2971998  29.7764103  37.1947568
  19.5497451  29.3708681  36.6498240  19.8679187  30.1884892  38.0108657
  20.7251759  30.5026632  36.6384067  23.3425747  24.7231452  33.6098679
  24.2209946  25.3122313  32.9751217  22.8678651  23.5401390  33.2304345
  22.1085344  23.1564303  33.7748447  23.4204054  22.7401881  32.1313567
  24.3798558  23.1636453  31.8344101  23.7057208  21.3080911  32.5872897
  22.7903788  20.8542358  32.9672344  24.3014389  20.4139862  31.5037426
  25.2057056  20.8733712  31.1044824  24.5448665  19.4407755  31.9299922
  23.5806207  20.2641999  30.6999642  24.6811154  21.3599424  33.5960815
  24.5589567  20.5709873  34.1292542  22.4998094  22.7390090  30.9198515
  21.3210682  22.3915901  30.9974508  23.0741357  23.1021796  29.7817496
  24.0550775  23.3368910  29.8348809  22.5913832  22.8019406  28.4463600
  21.5019069  22.7711068  28.4369071  23.0896600  23.8902317  27.4857449
  24.1769572  23.9461632  27.5393012  22.8307729  23.5998791  26.4675237
  22.5173508  25.2642618  27.7609499  23.0891370  26.0893931  28.7492444
  23.9460687  25.7453063  29.3093728  22.5298467  27.3484587  29.0281281
  22.9593276  27.9804259  29.7913785  21.3808493  27.7703879  28.3390454
  20.9214165  28.7179015  28.5788108  20.8060379  26.9460999  27.3578723
  19.9055955  27.2613127  26.8516742  21.3810921  25.7007213  27.0578135
  20.9301100  25.0689747  26.3068142  23.1493540  21.4498213  28.0001055
  24.2866437  21.0896949  28.3148213  22.3946270  20.7325880  27.1830286
  21.4606050  21.0504472  26.9668787  22.9075499  19.6044421  26.4055074
  23.9969367  19.6290718  26.4307233  22.4665165  18.2415806  26.9817689
  23.1509067  17.4751740  26.6178998  22.4480974  18.1772522  28.5154409
  21.6936004  18.8557276  28.9136815  22.2099311  17.1657986  28.8445405
  23.4217949  18.4600024  28.9153873  21.1742360  17.9058495  26.5475087
  20.8531191  17.2125711  27.1287759  22.4746669  19.7865246  24.9522828
  21.4328599  20.3889977  24.6834820  23.2719518  19.2770603  24.0152717
  24.1275661  18.8425213  24.3303094  22.9386304  19.1443439  22.5905998
  21.8664595  19.3059156  22.4792533  23.6399974  20.2081217  21.7189391
  23.3171285  21.1879626  22.0707082  25.1734610  20.1696820  21.7907142
  25.5503458  19.2163829  21.4202247  25.5854271  20.9723596  21.1791015
  25.5012026  20.3109246  22.8205968  23.2166745  20.0763896  20.2477450
  22.1304535  20.1277610  20.1725473  23.6446912  20.8919190  19.6649028
  23.5615661  19.1304729  19.8301304  23.2553747  17.7200336  22.1389656
  24.3241283  17.1896226  22.4502643  22.3254917  17.1013593  21.4149647
  21.4628676  17.6047904  21.2645414  22.3595884  15.6816069  21.0442942
  23.3783989  15.3088538  21.1495950  21.4731595  14.8682543  21.9896199
  20.4503684  15.2432523  21.9515779  21.4870927  13.3780201  21.6633635
  22.5192046  13.0380237  21.5787706  20.9663839  12.8241189  22.4445548
  20.9811926  13.1967441  20.7150477  21.9711120  14.9899714  23.2983789
  21.2499425  14.7597657  23.8887563  21.9010108  15.4667311  19.6064694
  20.7588635  15.7711215  19.2561441  22.7903876  14.9310123  18.7729397
  23.6751363  14.6601329  19.1779650  22.5200181  14.5031920  17.3881706
  21.7582275  15.1501256  16.9532524  23.7939469  14.6313662  16.5411761
  24.5580685  13.9645913  16.9407798  23.5619524  14.3185900  15.5231128
  24.3429807  16.0636439  16.4980876  23.5413670  16.7325480  16.1848754
  24.6613347  16.3700898  17.4945540  25.5295043  16.1914246  15.5354721
  26.4943499  15.3935083  15.6374988  25.5021590  17.1081692  14.6812464
  21.9692006  13.0748325  17.3021193  22.4166723  12.1959000  18.0719668
  21.0838157  12.8419272  16.4522577
