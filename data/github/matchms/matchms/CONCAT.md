# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Updated and extended plotting functionality, now located in `matchms.plotting`.
Contains three plot types: `plot_spectrum()` or `spectrum.plot()`, `plot_spectra_mirror()` or `spectrum.plot_against()` and `plot_spectra_array()` [#303](https://github.com/matchms/matchms/pull/303)

### Changed

- `Spectrum` objects got an update of the basic spectrum plots `spectrum.plot()` [#303](https://github.com/matchms/matchms/pull/303)
- `require_precursor_mz()` filter will now also discard nonsensical m/z values < 10.0 (value can be adapted by user) [#309](https://github.com/matchms/matchms/pull/309)

### Fixed

- Updated to new url for `load_from_usi` function (old link was broken) [#310](https://github.com/matchms/matchms/pull/310)

## [0.12.0] - 2022-01-18

### Added

- peak comments (as an `mz: comment` dictionary) are now part of metadata and can be addressed via a `Spectrum()` object `peak_comments` property [#284](https://github.com/matchms/matchms/pull/284)
- peak comments are dynamically updated whenever the respective peaks are changed [#277](https://github.com/matchms/matchms/pull/277)

### Changed

- Major refactoring of unit test layout now using a spectrum builder pattern [#261](https://github.com/matchms/matchms/pull/261)
- Spikes object now has different getitem method that allows to extract specific peaks as mz/intensity pair (or array) [#291](https://github.com/matchms/matchms/pull/291)
- `add_parent_mass()` filter now better handles existing entries (including fields "parent_mass", "exact_mass" and "parentmass") [#292](https://github.com/matchms/matchms/pull/292)
- minor improvement of compound name cleaning in `derive_adduct_from_name()` filter [#280](https://github.com/matchms/matchms/pull/280)
- `save_as_msp()` now writes peak comments (if present) to the output file [#277](https://github.com/matchms/matchms/pull/277)
- `load_from_msp()` now also reads peak comments [#277](https://github.com/matchms/matchms/pull/277)

### Fixed

- able to handle spectra containg empty/zero intensities [#289](https://github.com/matchms/matchms/pull/289)

## [0.11.0] - 2021-12-16

## Added

- better, more flexible string handling of `ModifiedCosine` [#275](https://github.com/matchms/matchms/pull/275)
- matchms logger, replacing all former `print` statments to better control logging output [#271](https://github.com/matchms/matchms/pull/271)
- `add_logging_to_file()`, `set_matchms_logger_level()`, `reset_matchms_logger()` functions to adapt logging output to user needs [#271](https://github.com/matchms/matchms/pull/271)

## Changed

- `save_as_msp()` can now also write to files with other than ".msp" extensions such as ".dat" [#276](https://github.com/matchms/matchms/pull/276)
- refactored `add_precursor_mz`, including better logging [#275](https://github.com/matchms/matchms/pull/275)

## [0.10.0] - 2021-11-21

### Added

- `Spectrum()` objects now also allows generating hashes, e.g. `hash(spectrum)` [#259](https://github.com/matchms/matchms/pull/259)
- `Spectrum()` objects can generate `.spectrum_hash()` and `.metadata_hash()` to track changes to peaks or metadata [#259](https://github.com/matchms/matchms/pull/259)
- `load_from_mgf()` now accepts both a path to a mgf file or a file-like object from a preloaded MGF file [#258](https://github.com/matchms/matchms/pull/258)
- `add_retention` filters with function `add_retention_time()` and `add_retention_index()` [#265](https://github.com/matchms/matchms/pull/265)

### Changed

- Code linting triggered by pylint update [#257](https://github.com/matchms/matchms/pull/257)
- Refactored `add_parent_mass()` filter can now also handle missing charge entries (if ionmode is known) [#252](https://github.com/matchms/matchms/pull/252)

## [0.9.2] - 2021-07-20

### Added

- Support for Python 3.9 [#240](https://github.com/matchms/matchms/issues/240)

### Changed

- Use `bool` instead of `numpy.bool` [#245](https://github.com/matchms/matchms/pull/245)

## [0.9.1] - 2021-06-16

### Fixed

- Correctly handle charge=0 entries in `add_parent_mass` filter [#236](https://github.com/matchms/matchms/pull/236)
- Reordered written metadata in MSP export for compatability with MS-FINDER & MS-DIAL [#230](https://github.com/matchms/matchms/pull/230)
- Update README.rst to fix fstring-quote python example [#226](https://github.com/matchms/matchms/pull/226)

## [0.9.0] - 2021-05-06

### Added

- new `matchms.networking` module which allows to build and export graphs from `scores` objects [#198](https://github.com/matchms/matchms/pull/198)
- Expand list of known negative ionmode adducts and conversion rules [#213](https://github.com/matchms/matchms/pull/213)
- `.to_numpy` method for Spikes class which allows to run `spectrum.peaks.to_numpy` [#214](https://github.com/matchms/matchms/issues/214)
- `save_as_msp()` function to export spectrums to .msp file [#215](https://github.com/matchms/matchms/pull/215)

### Changed

- `add_precursor_mz()` filter now also checks for metadata in keys `precursormz` and `precursor_mass` [#223](https://github.com/matchms/matchms/pull/223)
- `load_from_msp()` now handles .msp files containing multiple peaks per line separated by `;` [#221](https://github.com/matchms/matchms/pull/221)
- `add_parent_mass()` now includes `overwrite_existing_entry` option (default is False) [#225](https://github.com/matchms/matchms/pull/225)

### Fixed

- `add_parent_mass()` filter now makes consistent use of cleaned adducts [#225](https://github.com/matchms/matchms/pull/225)

## [0.8.2] - 2021-03-08

### Added

- Added filter function 'require_precursor_mz' and added 1 assert function in 'ModifiedCosine' [#191](https://github.com/matchms/matchms/pull/191)

- `make_charge_int()` to convert charge field to integer [#184](https://github.com/matchms/matchms/issues/184)

### Changed

- now deprecated: `make_charge_scalar()`, use `make_charge_int()` instead [#183](https://github.com/matchms/matchms/pull/183)

### Fixed

- Make `load_from_msp` work with different whitespaces [#192](https://github.com/matchms/matchms/issues/192)
- Very minor bugs in `add_parent_mass` [#188](https://github.com/matchms/matchms/pull/188)

## [0.8.1] - 2021-02-19

### Fixed

- Add package data to pypi tar.gz file (to fix Bioconda package) [#179](https://github.com/matchms/matchms/pull/179)

## [0.8.0] - 2021-02-16

### Added

- helper functions to clean adduct strings, `clean_adduct()` [#170](https://github.com/matchms/matchms/pull/170)

### Changed

- more thorough adduct cleaning effecting `derive_adduct_from_name()` and `derive_ionmode()` [#171](https://github.com/matchms/matchms/issues/171)
- significant expansion of `add_parent_mass()` filter to take known adduct properties into account [#170](https://github.com/matchms/matchms/pull/170)

## Fixed

- too unspecific formula detection (and removal) from given compound names in `derive_formula_from_name` [#172](https://github.com/matchms/matchms/issues/172)
- no longer ignore n_max setting in `reduce_to_number_of_peaks` filter [#177](https://github.com/matchms/matchms/issues/177)

## [0.7.0] - 2021-01-04

### Added

- `scores_by_query` and `scores_by reference` now accept sort=True to return sorted scores [#153](https://github.com/matchms/matchms/pull/153)

### Changed

- `Scores.scores` is now returning a structured array [#153](https://github.com/matchms/matchms/pull/153)

### Fixed

- Minor bug in `add_precursor_mz` [#161](https://github.com/matchms/matchms/pull/161)
- Minor bug in `Spectrum` class (missing metadata deepcopy) [#153](https://github.com/matchms/matchms/pull/153)
- Minor bug in `Spectrum` class (__eq__ method was not working with numpy arrays in metadata) [#153](https://github.com/matchms/matchms/pull/153)

## [0.6.2] - 2020-12-03

### Changed

- Considerable performance improvement for CosineGreedy and CosineHungarian [#159](https://github.com/matchms/matchms/pull/159)

## [0.6.1] - 2020-11-26

### Added

- PrecursorMzMatch for deriving precursor m/z matches within a given tolerance [#156](https://github.com/matchms/matchms/pull/156)

### Changed

- Raise error for improper use of reduce_to_number_of_peaks filter [#151](https://github.com/matchms/matchms/pull/151)
- Renamed ParentmassMatch to ParentMassMatch [#156](https://github.com/matchms/matchms/pull/156)

### Fixed

- Fix minor issue with msp importer to avoid failing with unknown characters [#151](https://github.com/matchms/matchms/pull/151)

## [0.6.0] - 2020-09-14

### Added

- Four new peak filtering functions [#119](https://github.com/matchms/matchms/pull/119)
- score_by_reference and score_by_query methods to Scores [#142](https://github.com/matchms/matchms/pull/142)
- is_symmetric option to speed up all-vs-all type score calculation [#59](https://github.com/matchms/matchms/issues/59)
- Support for Python 3.8 [#145](https://github.com/matchms/matchms/pull/145)

### Changed

- Refactor similarity scores to be instances of BaseSimilarity class [#135](https://github.com/matchms/matchms/issues/135)
- Marked Scores.calculate() method as deprecated [#135](https://github.com/matchms/matchms/issues/135)

### Removed

- calculate_parallel function [#135](https://github.com/matchms/matchms/issues/135)
- Scores.calculate_parallel method [#135](https://github.com/matchms/matchms/issues/135)
- similarity.FingerprintSimilarityParallel class (now part of similarity.FingerprintSimilarity) [#135](https://github.com/matchms/matchms/issues/135)
- similarity.ParentmassMatchParallel class (now part of similarity.ParentmassMatch) [#135](https://github.com/matchms/matchms/issues/135)

## [0.5.2] - 2020-08-26

### Changed

- Revision of JOSS manuscript [#137](https://github.com/matchms/matchms/pull/137)

## [0.5.1] - 2020-08-19

### Added

- Basic submodule documentation and more code examples [#128](https://github.com/matchms/matchms/pull/128)

### Changed

- Extended, updated, and corrected documentation for filter functions [#118](https://github.com/matchms/matchms/pull/118)

## [0.5.0] - 2020-08-05

### Added

- Read mzML and mzXML files to create Spectrum objects from it [#110](https://github.com/matchms/matchms/pull/110)
- Read msp files to create Spectrum objects from it [#102](https://github.com/matchms/matchms/pull/102)
- Peak weighting option for CosineGreedy and ModifiedCosine score [#96](https://github.com/matchms/matchms/issues/96)
- Peak weighting option for CosineHungarian score [#112](https://github.com/matchms/matchms/pull/112)
- Similarity score based on comparing parent masses [#79](https://github.com/matchms/matchms/pull/79)
- Method for instantiating a spectrum from the metabolomics USI [#93](https://github.com/matchms/matchms/pull/93)

### Changed

- CosineGreedy function is now numba based [#86](https://github.com/matchms/matchms/pull/86)
- Extended readthedocs documentation [#82](https://github.com/matchms/matchms/issues/82)

### Fixed

- Incorrect denominator for cosine score normalization [#98](https://github.com/matchms/matchms/pull/98)

## [0.4.0] - 2020-06-11

### Added

- Filter add_fingerprint to derive molecular fingerprints [#42](https://github.com/matchms/matchms/issues/42)
- Similarity scores based on molecular fingerprints [#42](https://github.com/matchms/matchms/issues/42)
- Add extensive compound name cleaning and harmonization [#23](https://github.com/matchms/matchms/issues/23)
- Faster cosine score implementation using numba [#29](https://github.com/matchms/matchms/issues/29)
- Cosine score based on Hungarian algorithm [#40](https://github.com/matchms/matchms/pull/40)
- Modified cosine score [#26](https://github.com/matchms/matchms/issues/26)
- Import and export of spectrums from json files [#15](https://github.com/matchms/matchms/issues/15)
- Doc strings for many methods [#49](https://github.com/matchms/matchms/issues/49)
- Examples in doc strings which are tested on CI [#49](https://github.com/matchms/matchms/issues/49)

### Changed

- normalize_intensities filter now also normalizes losses [#69](https://github.com/matchms/matchms/issues/69)

### Removed

## [0.3.4] - 2020-05-29

### Changed

- Fix verify step in conda publish workflow
- Fixed mixed up loss intensity order. [#20](https://github.com/matchms/matchms/issues/20)

## [0.3.3] - 2020-05-27

### Added

- Build workflow runs the tests after installing the package [#47](https://github.com/matchms/matchms/pull/47)

### Changed

- tests were removed from the package (see setup.py) [#47](https://github.com/matchms/matchms/pull/47)

## [0.3.2] - 2020-05-26

### Added

- Workflow improvements
  - Use artifacts in build workflow
  - List artifact folder in build workflow

### Changed

- Workflow improvements [#244](https://github.com/matchms/matchms-backup/pull/244)
  - merge anaconda and python build workflows
  - fix conda package install command in build workflow
  - publish only on ubuntu machine
  - update workflow names
  - test conda packages on windows and unix separately
  - install conda package generated by the workflow
  - split workflows into multiple parts
  - use default settings for conda action
- data folder is handled by setup.py but not meta.yml

### Removed

- remove python build badge [#244](https://github.com/matchms/matchms-backup/pull/244)
- Moved ``spec2vec`` similarity related functionality from ``matchms`` to [iomega/spec2vec](https://github.com/iomega/spec2vec)
- removed build step in build workflow
- removed conda build scripts: conda/build.sh and conda/bld.bat
- removed conda/condarc.yml
- removed conda_build_config.yaml
- removed testing from publish workflow

## [0.3.1] - 2020-05-19

### Added

- improve conda package [#225](https://github.com/matchms/matchms/pull/225)
  - Build scripts for Windows and Unix(MacOS and Linux) systems
  - verify conda package after uploading to anaconda repository by installing it
  - conda package also includes `matchms/data` folder

### Changed

- conda package fixes [#223](https://github.com/matchms/matchms/pull/223)
  - move conda receipe to conda folder
  - fix conda package installation issue
  - add extra import tests for conda package
  - add instructions to build conda package locally
  - automatically find matchms package in setup.py
  - update developer instructions
  - increase verbosity while packaging
  - skip builds for Python 2.X
  - more flexible package versions
  - add deployment requirements to meta.yml
- verify conda package [#225](https://github.com/matchms/matchms/pull/225)
  - use conda/environment.yml when building the package
- split anaconda workflow [#225](https://github.com/matchms/matchms/pull/225)
  - conda build: tests conda packages on every push and pull request
  - conda publish: publish and test conda package on release
  - update the developer instructions
  - move conda receipe to conda folder

## [0.3.0] - 2020-05-13

### Added

- Spectrum, Scores class, save_to_mgf, load_from_mgf, normalize_intensities, calculate_scores [#66](https://github.com/matchms/matchms/pull/66) [#67](https://github.com/matchms/matchms/pull/67) [#103](https://github.com/matchms/matchms/pull/103) [#108](https://github.com/matchms/matchms/pull/108) [#113](https://github.com/matchms/matchms/pull/113) [#115](https://github.com/matchms/matchms/pull/115) [#151](https://github.com/matchms/matchms/pull/151) [#152](https://github.com/matchms/matchms/pull/152) [#121](https://github.com/matchms/matchms/pull/121) [#154](https://github.com/matchms/matchms/pull/154) [#134](https://github.com/matchms/matchms/pull/134) [#159](https://github.com/matchms/matchms/pull/159) [#161](https://github.com/matchms/matchms/pull/161) [#198](https://github.com/matchms/matchms/pull/198)
- Spikes class [#150](https://github.com/matchms/matchms/pull/150) [#167](https://github.com/matchms/matchms/pull/167)
- Anaconda package [#70](https://github.com/matchms/matchms/pull/70) [#68](https://github.com/matchms/matchms/pull/68) [#181](https://github.com/matchms/matchms/pull/181)
- Sonarcloud [#80](https://github.com/matchms/matchms/pull/80) [#79](https://github.com/matchms/matchms/pull/79) [#149](https://github.com/matchms/matchms/pull/149) [#169](https://github.com/matchms/matchms/pull/169)
- Normalization filter [#83](https://github.com/matchms/matchms/pull/83)
- SpeciesString filter [#181](https://github.com/matchms/matchms/pull/181)
- Select by relative intensity filter [#98](https://github.com/matchms/matchms/pull/98)
- Select-by capability based on mz and intensity [#87](https://github.com/matchms/matchms/pull/87)
- Default filters [#97](https://github.com/matchms/matchms/pull/97)
- integration test [#89](https://github.com/matchms/matchms/pull/89) [#147](https://github.com/matchms/matchms/pull/147) [#156](https://github.com/matchms/matchms/pull/156) [#194](https://github.com/matchms/matchms/pull/194)
- cosine greedy similarity function [#112](https://github.com/matchms/matchms/pull/112)
- parent mass filter [#116](https://github.com/matchms/matchms/pull/116) [#122](https://github.com/matchms/matchms/pull/122) [#158](https://github.com/matchms/matchms/pull/158)
- require_minimum_number_of_peaks filter [#131](https://github.com/matchms/matchms/pull/131) [#155](https://github.com/matchms/matchms/pull/155)
- reduce_to_number_of_peaks filter [#209](https://github.com/matchms/matchms/pull/209)
- inchi filters [#145](https://github.com/matchms/matchms/pull/145) [#127](https://github.com/matchms/matchms/pull/127) [#181](https://github.com/matchms/matchms/pull/181)
- losses [#160](https://github.com/matchms/matchms/pull/160)
- vesion string checks [#185](https://github.com/matchms/matchms/pull/185)
- Spec2Vec [#183](https://github.com/matchms/matchms/pull/183) [#165](https://github.com/matchms/matchms/pull/165)
- functions to verify inchies [#181](https://github.com/matchms/matchms/pull/181) [#180](https://github.com/matchms/matchms/pull/180)
- documentation using radthedocs [#196](https://github.com/matchms/matchms/pull/196) [#197](https://github.com/matchms/matchms/pull/197)
- build status badges [#174](https://github.com/matchms/matchms/pull/174)
- vectorize spec2vec [#206](https://github.com/matchms/matchms/pull/206)

### Changed

- Seperate filters [#97](https://github.com/matchms/matchms/pull/97)
- Translate filter steps to new structure (interpret charge and ionmode) [#73](https://github.com/matchms/matchms/pull/73)
- filters returning a new spectrum [#100](https://github.com/matchms/matchms/pull/100)
- Flowchart diagram [#135](https://github.com/matchms/matchms/pull/135)
- numpy usage [#191](https://github.com/matchms/matchms/pull/191)
- consistency of the import statements [#189](https://github.com/matchms/matchms/pull/189)

## [0.2.0] - 2020-04-03

### Added

- Anaconda actions

## [0.1.0] - 2020-03-19

### Added

- This is the initial version of Spec2Vec from https://github.com/iomega/Spec2Vec
- (later splitted into matchms + spec2vec)

[Unreleased]: https://github.com/matchms/matchms/compare/0.12.0...HEAD
[0.12.0]: https://github.com/matchms/matchms/compare/0.11.0...0.12.0
[0.11.0]: https://github.com/matchms/matchms/compare/0.10.0...0.11.0
[0.10.0]: https://github.com/matchms/matchms/compare/0.9.2...0.10.0
[0.9.2]: https://github.com/matchms/matchms/compare/0.9.0...0.9.2
[0.9.1]: https://github.com/matchms/matchms/compare/0.9.0...0.9.1
[0.9.0]: https://github.com/matchms/matchms/compare/0.8.2...0.9.0
[0.8.2]: https://github.com/matchms/matchms/compare/0.8.1...0.8.2
[0.8.1]: https://github.com/matchms/matchms/compare/0.8.0...0.8.1
[0.8.0]: https://github.com/matchms/matchms/compare/0.7.0...0.8.0
[0.7.0]: https://github.com/matchms/matchms/compare/0.6.2...0.7.0
[0.6.2]: https://github.com/matchms/matchms/compare/0.6.1...0.6.2
[0.6.1]: https://github.com/matchms/matchms/compare/0.6.0...0.6.1
[0.6.0]: https://github.com/matchms/matchms/compare/0.5.2...0.6.0
[0.5.2]: https://github.com/matchms/matchms/compare/0.5.1...0.5.2
[0.5.1]: https://github.com/matchms/matchms/compare/0.5.0...0.5.1
[0.5.0]: https://github.com/matchms/matchms/compare/0.4.0...0.5.0
[0.4.0]: https://github.com/matchms/matchms/compare/0.3.4...0.4.0
[0.3.4]: https://github.com/matchms/matchms/compare/0.3.3...0.3.4
[0.3.3]: https://github.com/matchms/matchms/compare/0.3.2...0.3.3
[0.3.2]: https://github.com/matchms/matchms/compare/0.3.1...0.3.2
[0.3.1]: https://github.com/matchms/matchms/compare/0.3.0...0.3.1
[0.3.0]: https://github.com/matchms/matchms/compare/0.2.0...0.3.0
[0.2.0]: https://github.com/matchms/matchms/compare/0.1.0...0.2.0
[0.1.0]: https://github.com/matchms/matchms/releases/tag/0.1.0
# Contributing guidelines

We welcome any kind of contribution to our software, from simple comment or question to a full fledged [pull request](https://help.github.com/articles/about-pull-requests/). Please read and follow our [Code of Conduct](CODE_OF_CONDUCT.rst).

A contribution can be one of the following cases:

1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation);
1. you want to make a new release of the code base.

The sections below outline the steps in each case.

## You have a question

1. use the search functionality [here](https://github.com/matchms/matchms/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue;
1. apply the "Question" label; apply other labels when relevant.

## You think you may have found a bug

1. use the search functionality [here](https://github.com/matchms/matchms/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the [SHA hashcode](https://help.github.com/articles/autolinked-references-and-urls/#commit-shas) of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
1. apply relevant labels to the newly created issue.

## You want to make some kind of change to the code base

1. (**important**) announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
1. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/) and [here](https://help.github.com/articles/syncing-a-fork/));
1. make sure the existing tests still work by running ``pytest``;
1. add your own tests (if necessary);
1. update or expand the documentation;
1. update the `CHANGELOG.md` file with change;
1. [push](http://rogerdudler.github.io/git-guide/>) your feature branch to (your fork of) the matchms repository on GitHub;
1. create the pull request, e.g. following the instructions [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.

## You want to make a new release of the code base

To create release you need write permission on the repository.

1. Check author list in `citation.cff` and `.zenodo.json` files
1. Bump the version using `bump2version <major|minor|patch>`. For example, `bump2version major` will increase major version numbers everywhere its needed (code, meta, etc.) in the repo.
1. Update the `CHANGELOG.md` to include changes made
1. Goto [GitHub release page](https://github.com/matchms/matchms/releases)
1. Press draft a new release button
1. Fill version, title and description field
1. Press the Publish Release button
1. Wait until [PyPi publish workflow](https://github.com/matchms/matchms/actions/workflows/CI_publish_pypi.yml) has completed
1. Verify new release is on [PyPi](https://pypi.org/project/matchms/#history)
1. Wait until new release is also on Bioconda (https://anaconda.org/bioconda/matchms) via a automaticly created PR on [bioconda recipes repo](https://github.com/bioconda/bioconda-recipes/pulls?q=is%3Apr+is%3Aopen+matchms)
1. Test matchms from bioconda by manually running [Conda verify](https://github.com/matchms/matchms/actions/workflows/conda_verify.yml) workflow

Also a Zenodo entry will be made for the release with its own DOI.
---
title: matchms - processing and similarity evaluation of mass spectrometry data.
tags:
  - Python
  - mass spectrometry
  - metadata cleaning
  - data processing
  - similarity measures
  - metabolomics

authors:
  - name: Florian Huber
    orcid: 0000-0002-3535-9406
    affiliation: 1
  - name: Stefan Verhoeven
    orcid: 0000-0002-5821-2060
    affiliation: 1
  - name: Christiaan Meijer
    orcid: 0000-0002-5529-5761
    affiliation: 1
  - name: Hanno Spreeuw
    orcid: 0000-0002-5057-0322
    affiliation: 1
  - name: Efraín Manuel Villanueva Castilla
    orcid: 0000-0001-7665-3575
    affiliation: 2
  - name: Cunliang Geng
    orcid: 0000-0002-1409-8358
    affiliation: 1
  - name: Justin J. J. van der Hooft
    orcid: 0000-0002-9340-5511
    affiliation: 3
  - name: Simon Rogers
    orcid: 0000-0003-3578-4477
    affiliation: 2
  - name: Adam Belloum
    orcid: 0000-0001-6306-6937
    affiliation: 1
  - name: Faruk Diblen
    orcid: 0000-0002-0989-929X
    affiliation: 1
  - name: Jurriaan H. Spaaks
    orcid: 0000-0002-7064-4069
    affiliation: 1

affiliations:
 - name: Netherlands eScience Center, Science Park 140, 1098XG Amsterdam, The Netherlands
   index: 1
 - name: School of Computing Science, University of Glasgow, Glasgow, United Kingdom
   index: 2
 - name: Bioinformatics Group, Plant Sciences Group, University of Wageningen, Wageningen, the Netherlands
   index: 3
date: 16 June 2020
bibliography: paper.bib

---

# Summary

Mass spectrometry data is at the heart of numerous applications in the biomedical and life sciences.
With growing use of high-throughput techniques, researchers need to analyze larger and more complex datasets. In particular through joint effort in the research community, fragmentation mass spectrometry datasets are growing in size and number.
Platforms such as MassBank [@horai_massbank_2010], GNPS [@Wang2016] or MetaboLights [@haug_metabolights_2020] serve as an open-access hub for sharing of raw, processed, or annotated fragmentation mass spectrometry data.
Without suitable tools, however, exploitation of such datasets remains overly challenging. 
In particular, large collected datasets contain data acquired using different instruments and measurement conditions, and can further contain a significant fraction of inconsistent, wrongly labeled, or incorrect metadata (annotations).

``matchms`` is an open-source Python package to import, process, clean, and compare mass spectrometry data (MS/MS) (see \autoref{fig:flowchart}).
It allows to implement and run an easy-to-follow, easy-to-reproduce workflow from raw mass spectra to pre- and post-processed spectral data. 
Raw data can be imported from the commonly used formats msp, mzML [@martens_mzmlcommunity_2011], mzXML, MGF (mzML, mzXML, MGF file importers are built on top of pyteomics [@levitsky_pyteomics_2019;@goloborodko_pyteomicspython_2013], as well as from JSON files (as provided by GNPS), but also via Universal Spectrum Identifiers (USI) [@wang_interactive_2020]. Further data formats or more extensive options regarding metadata parsing can best be handled by using pyteomics [@levitsky_pyteomics_2019] or pymzml [@kosters_pymzml_2018].
``matchms`` contains numerous metadata cleaning and harmonizing filter functions that can easily be stacked to construct a desired pipeline (\autoref{fig:filtering}), which can also easily be extended by custom functions wherever needed. Available filters include extensive cleaning, correcting, checking of key metadata fields such as compound name, structure annotations (InChI, SMILES, InChIKey), ionmode, adduct, or charge.
Many of the provided metadata cleaning filters were designed for handling and improving GNPS-style MGF or JSON datasets. For future versions, however, we aim to further extend this to other commonly used public databases.

![Flowchart of ``matchms`` workflow. Reference and query spectrums are filtered using the same set of set filters (here: filter A and filter B). Once filtered, every reference spectrum is compared to every query spectrum using the ``matchms.Scores`` object. \label{fig:flowchart}](flowchart_matchms.png)

Current Python tools for working with MS/MS data include pyOpenMS [@rost_pyopenms_2014], a wrapper for OpenMS [@rost_openms_2016] with a strong focus on processing and filtering of raw mass spectral data. 
pyOpenMS has a wide range of peak processing functions which can be used to further complement a ``matchms`` filtering pipeline.
Another, more lightweight and native Python package with a focus on spectra visualization is ``spectrum_utils`` [@bittremieux_spectrum_utils_2020].
``matchms`` focuses on comparing and linking large number of mass spectra. Many of its built-in filters are aimed at handling large mass spectra datasets from common public data libraries such as GNPS.

``matchms`` provides functions to derive different similarity scores between spectra. Those include the established spectra-based measures of the cosine score or modified cosine score [@watrous_mass_2012].
The package also offers fast implementations of common similarity measures (Dice, Jaccard, Cosine) that can be used to compute similarity scores between molecular fingerprints (rdkit, morgan1, morgan2, morgan3, all implemented using rdkit [@rdkit]).
``matchms`` facilitates easily deriving similarity measures between large number of spectra at comparably fast speed due to score implementations based on NumPy [@van_der_walt_numpy_2011], SciPy [@2020SciPy-NMeth], and Numba [@lam_numba_2015]. Additional similarity measures can easily be added using the ``matchms`` API. 
The provided API also allows to quickly compare, sort, and inspect query versus reference spectra using either the included similarity scores or added custom measures.
The API was designed to be easily extensible so that users can add their own filters for spectra processing, or their own similarity functions for spectral comparisons.
The present set of filters and similarity functions was mostly geared towards smaller molecules and natural compounds, but it could easily be extended by functions specific to larger peptides or proteins.

``matchms`` is freely accessible either as conda package (https://anaconda.org/nlesc/matchms), or in form of source-code on GitHub (https://github.com/matchms/matchms). For further code examples and documentation see https://matchms.readthedocs.io/en/latest/.
All main functions are covered by tests and continuous integration to offer reliable functionality.
We explicitly value future contributions from a mass spectrometry interested community and hope that ``matchms`` can serve as a reliable and accessible entry point for handling complex mass spectrometry datasets using Python. 


# Example workflow
A typical workflow with ``matchms`` looks as indicated in \autoref{fig:flowchart}, or as described in the following code example.
```python
from matchms.importing import load_from_mgf
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms import calculate_scores
from matchms.similarity import CosineGreedy

# Read spectrums from a MGF formatted file
file = load_from_mgf("all_your_spectrums.mgf")

# Apply filters to clean and enhance each spectrum
spectrums = []
for spectrum in file:
    spectrum = default_filters(spectrum)
    spectrum = normalize_intensities(spectrum)
    spectrums.append(spectrum)

# Calculate Cosine similarity scores between all spectrums
scores = calculate_scores(references=spectrums,
                          queries=spectrums,
                          similarity_function=CosineGreedy())

# Print the calculated scores for each spectrum pair
for score in scores:
    (reference, query, score, n_matching) = score
    # Ignore scores between same spectrum and
    # pairs which have less than 20 peaks in common
    if reference is not query and n_matching >= 20:
        print(f"Reference scan id: {reference.metadata['scans']}")
        print(f"Query scan id: {query.metadata['scans']}")
        print(f"Score: {score:.4f}")
        print(f"Number of matching peaks: {n_matching}")
        print("----------------------------")
```

![``matchms`` provides a range of filter functions to process spectrum peaks and metadata. Filters can easily be stacked and combined to build a desired pipeline. The API also makes it easy to extend customer pipelines by adding own filter functions. \label{fig:filtering}](filtering_sketch.png)

# Processing spectrum peaks and plotting
``matchms`` provides numerous filters to process mass spectra peaks. Below a simple example to remove low intensity peaks from a spectrum (\autoref{fig:peak_filtering}).
```python
from matchms.filtering import require_minimum_number_of_peaks
from matchms.filtering import select_by_mz
from matchms.filtering import select_by_relative_intensity

def process_peaks(s):
    s = select_by_mz(s, mz_from=0, mz_to=1000)
    s = select_by_relative_intensity(s, intensity_from=0.001)
    s = require_minimum_number_of_peaks(s, n_required=10)
    return s

# Apply processing steps to spectra (here to a single "spectrum_raw")
spectrum_processed = process_peaks(spectrum_raw)

# Plot raw spectrum (all and zoomed in)
spectrum_raw.plot()
spectrum_raw.plot(intensity_to=0.02)

# Plot processed spectrum (all and zoomed in)
spectrum_processed.plot()
spectrum_processed.plot(intensity_to=0.02)
```

![Example of ``matchms`` peak filtering applied to an actual spectrum using ``select_by_relative_intensity`` to remove peaks of low relative intensity. Spectra are plotted using the provided ``spectrum.plot()`` function. \label{fig:peak_filtering}](peak_filtering.png)


# References
---
name: Custom issue template
about: Describe this issue template's purpose here.
title: ''
labels: ''
assignees: ''

---


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
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. 

**Describe the solution you'd like**
A clear and concise description of what would be needed, what it has to do.

**Any good starting points?**
Do you know of any alternative solutions or do you have some helping code/links at hand?

**Additional context**
Add any other context or screenshots about the feature request here.

**Scientific reference**
If this request is about adding a method to matchms (new similarity measure, new filter etc.), please point us to a scientific reference (if possible).
###############################################################################
Contributor Covenant Code of Conduct
###############################################################################

Our Pledge
**********

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance, race,
religion, or sexual identity and orientation.

Our Standards
*************

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

Our Responsibilities
********************

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

Scope
*****

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

Enforcement
***********

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at generalization@esciencecenter.nl. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

Attribution
***********

This Code of Conduct is adapted from the `Contributor Covenant <https://www.contributor-covenant.org>`_, version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html
.. image:: readthedocs/_static/matchms_header.png
   :target: readthedocs/_static/matchms.png
   :align: left
   :alt: matchms

Matchms is an open-source Python package to import, process, clean, and compare mass spectrometry data (MS/MS). It allows to implement and run an easy-to-follow, easy-to-reproduce workflow from raw mass spectra to pre- and post-processed spectral data. Spectral data can be imported from common formats such mzML, mzXML, msp, metabolomics-USI, MGF, or json (e.g. GNPS-syle json files). Matchms then provides filters for metadata cleaning and checking, as well as for basic peak filtering. Finally, matchms was build to import and apply different similarity measures to compare large amounts of spectra. This includes common Cosine scores, but can also easily be extended by custom measures. Example for spectrum similarity measures that were designed to work in matchms are `Spec2Vec <https://github.com/iomega/spec2vec>`_ and `MS2DeepScore <https://github.com/matchms/ms2deepscore>`_.

If you use matchms in your research, please cite the following software paper:  

F Huber, S. Verhoeven, C. Meijer, H. Spreeuw, E. M. Villanueva Castilla, C. Geng, J.J.J. van der Hooft, S. Rogers, A. Belloum, F. Diblen, J.H. Spaaks, (2020). matchms - processing and similarity evaluation of mass spectrometry data. Journal of Open Source Software, 5(52), 2411, https://doi.org/10.21105/joss.02411

.. list-table::
   :widths: 25 25
   :header-rows: 1

   * - 
     - Badges
   * - **fair-software.nl recommendations**
     - 
   * - \1. Code repository
     - |GitHub Badge|
   * - \2. License
     - |License Badge|
   * - \3. Community Registry
     - |Conda Badge| |Pypi Badge| |Research Software Directory Badge|
   * - \4. Enable Citation
     - |JOSS Badge| |Zenodo Badge|
   * - \5. Checklists
     - |CII Best Practices Badge| |Howfairis Badge|
   * - **Code quality checks**
     -
   * - Continuous integration
     - |CI Build|
   * - Documentation
     - |ReadTheDocs Badge|
   * - Code Quality
     - |Sonarcloud Quality Gate Badge| |Sonarcloud Coverage Badge|


.. |GitHub Badge| image:: https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue
   :target: https://github.com/matchms/matchms
   :alt: GitHub Badge

.. |License Badge| image:: https://img.shields.io/github/license/matchms/matchms
   :target: https://github.com/matchms/matchms
   :alt: License Badge

.. |Conda Badge| image:: https://anaconda.org/bioconda/matchms/badges/version.svg
   :target: https://anaconda.org/bioconda/matchms
   :alt: Conda Badge

.. |Pypi Badge| image:: https://img.shields.io/pypi/v/matchms?color=blue
   :target: https://pypi.org/project/matchms/
   :alt: Pypi Badge

.. |Research Software Directory Badge| image:: https://img.shields.io/badge/rsd-matchms-00a3e3.svg
   :target: https://www.research-software.nl/software/matchms
   :alt: Research Software Directory Badge

.. |Zenodo Badge| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3859772.svg
   :target: https://doi.org/10.5281/zenodo.3859772
   :alt: Zenodo Badge

.. |JOSS Badge| image:: https://joss.theoj.org/papers/10.21105/joss.02411/status.svg
   :target: https://doi.org/10.21105/joss.02411
   :alt: JOSS Badge

.. |CII Best Practices Badge| image:: https://bestpractices.coreinfrastructure.org/projects/3792/badge
   :target: https://bestpractices.coreinfrastructure.org/projects/3792
   :alt: CII Best Practices Badge

.. |Howfairis Badge| image:: https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green
   :target: https://fair-software.eu
   :alt: Howfairis badge

.. |CI Build| image:: https://github.com/matchms/matchms/actions/workflows/CI_build.yml/badge.svg
    :alt: Continuous integration workflow
    :target: https://github.com/matchms/matchms/actions/workflows/CI_build.yml

.. |ReadTheDocs Badge| image:: https://readthedocs.org/projects/matchms/badge/?version=latest
    :alt: Documentation Status
    :scale: 100%
    :target: https://matchms.readthedocs.io/en/latest/?badge=latest

.. |Sonarcloud Quality Gate Badge| image:: https://sonarcloud.io/api/project_badges/measure?project=matchms_matchms&metric=alert_status
   :target: https://sonarcloud.io/dashboard?id=matchms_matchms
   :alt: Sonarcloud Quality Gate

.. |Sonarcloud Coverage Badge| image:: https://sonarcloud.io/api/project_badges/measure?project=matchms_matchms&metric=coverage
   :target: https://sonarcloud.io/component_measures?id=matchms_matchms&metric=Coverage&view=list
   :alt: Sonarcloud Coverage


**********************************
Latest changes (matchms >= 0.11.0)
**********************************

Matchms now allows proper logging. Matchms functions and method report unexpected or 
undesired behavior as logging WARNING, and additional information as INFO.
The default logging level is set to WARNING. If you want to output additional
logging messages, you can lower the logging level to INFO using set_matchms_logger_level:

.. code-block:: python

    from matchms import set_matchms_logger_level

    set_matchms_logger_level("INFO")

If you want to suppress logging warnings, you can also raise the logging level
to ERROR by:

.. code-block:: python

    set_matchms_logger_level("ERROR")

To write logging entries to a local file, you can do the following:

.. code-block:: python

    from matchms.logging_functions import add_logging_to_file

    add_logging_to_file("sample.log", loglevel="INFO")

If you want to write the logging messages to a local file while silencing the
stream of such messages, you can do the following:

.. code-block:: python

    from matchms.logging_functions import add_logging_to_file

    add_logging_to_file("sample.log", loglevel="INFO",
                        remove_stream_handlers=True)


**********************************
Latest changes (matchms >= 0.10.0)
**********************************

- 2 new filters in ``matchms.filtering``: ``add_retention_time()`` and ``add_retention_index()``, to consistently add retention time/index to the spectrum metadata
- Hashes! ``Spectrum``-objects now allow to compute different hashes:

.. code-block:: python

    from matchms.importing import load_from_mgf

    # Read spectrums from a MGF formatted file, for other formats see https://matchms.readthedocs.io/en/latest/api/matchms.importing.html 
    spectrums = list(load_from_mgf("tests/pesticides.mgf"))
    
    # Spectrum hashes are generated based on MS/MS peak m/z and intensities
    # Those will change if any processing step affects the peaks.
    spectrum_hashes = [s.spectrum_hash() for s in spectrums]
    # Metadata hashes are generated based on the spectrum metadata
    # Those will change if any processing step affects the metadata.
    metadata_hashes = [s.metadata_hash() for s in spectrums]
    # `hash(spectrum)` will return a hash that is a combination of spectrum and metadata hash
    # Those will hence change if any processing step affects the peaks and/or the metadata.
    hashes = [hash(s) for s in spectrums]

***********************
Documentation for users
***********************
For more extensive documentation `see our readthedocs <https://matchms.readthedocs.io/en/latest/>`_ and our `matchms introduction tutorial <https://blog.esciencecenter.nl/build-your-own-mass-spectrometry-analysis-pipeline-in-python-using-matchms-part-i-d96c718c68ee>`_.

Installation
============

Prerequisites:  

- Python 3.7, 3.8 or 3.9
- Anaconda (recommended)

We recommend installing matchms from Anaconda Cloud with

.. code-block:: console

  # install matchms in a new virtual environment to avoid dependency clashes
  conda create --name matchms python=3.8
  conda activate matchms
  conda install --channel bioconda --channel conda-forge matchms

Alternatively, matchms can also be installed using ``pip`` but users will then either have to install ``rdkit`` on their own or won't be able to use the entire functionality. Without ``rdkit`` installed several filter functions related to processing and cleaning chemical metadata will not run.
To install matchms with ``pip`` simply run

.. code-block:: console

  pip install matchms

matchms ecosystem -> additional functionalities
===============================================

Matchms functionalities can be complemented by additional packages.  
To date we are aware of:

+ `Spec2Vec <https://github.com/iomega/spec2vec>`_ an alternative machine-learning spectral similarity score that can simply be installed by `pip install spec2vec` and be imported as `from spec2vec import Spec2Vec` following the same API as the scores in `matchms.similarity`.

+ `MS2DeepScore <https://github.com/matchms/ms2deepscore>`_ a supervised, deep-learning based spectral similarity score that can simply be installed by `pip install ms2deepscore` and be imported as `from ms2deepscore import MS2DeepScore` following the same API as the scores in `matchms.similarity`.

+ `matchmsextras <https://github.com/matchms/matchmsextras>`_ which contains additional functions to create networks based on spectral similarities, to run spectrum searchers against `PubChem`, or additional plotting methods.

+ `memo <https://github.com/mandelbrot-project/memo>`_ a method allowing a Retention Time (RT) agnostic alignment of metabolomics samples using the fragmentation spectra (MS2) of their consituents.

*(if you know of any other packages that are fully compatible with matchms, let us know!)*

Introduction
============

To get started with matchms, we recommend following our `matchms introduction tutorial <https://blog.esciencecenter.nl/build-your-own-mass-spectrometry-analysis-pipeline-in-python-using-matchms-part-i-d96c718c68ee>`_.

Alternatively, here below is a small example of using matchms to calculate the Cosine score between mass Spectrums in the `tests/pesticides.mgf <https://github.com/matchms/matchms/blob/master/tests/pesticides.mgf>`_ file.

.. code-block:: python

    from matchms.importing import load_from_mgf
    from matchms.filtering import default_filters
    from matchms.filtering import normalize_intensities
    from matchms import calculate_scores
    from matchms.similarity import CosineGreedy

    # Read spectrums from a MGF formatted file, for other formats see https://matchms.readthedocs.io/en/latest/api/matchms.importing.html 
    file = load_from_mgf("tests/pesticides.mgf")

    # Apply filters to clean and enhance each spectrum
    spectrums = []
    for spectrum in file:
        # Apply default filter to standardize ion mode, correct charge and more.
        # Default filter is fully explained at https://matchms.readthedocs.io/en/latest/api/matchms.filtering.html .
        spectrum = default_filters(spectrum)
        # Scale peak intensities to maximum of 1
        spectrum = normalize_intensities(spectrum)
        spectrums.append(spectrum)

    # Calculate Cosine similarity scores between all spectrums
    # For other similarity score methods see https://matchms.readthedocs.io/en/latest/api/matchms.similarity.html .
    scores = calculate_scores(references=spectrums,
                              queries=spectrums,
                              similarity_function=CosineGreedy())

    # Print the calculated scores for each spectrum pair
    for score in scores:
        (reference, query, score) = score
        # Ignore scores between same spectrum and
        # pairs which have less than 20 peaks in common
        if reference is not query and score["matches"] >= 20:
            print(f"Reference scan id: {reference.metadata['scans']}")
            print(f"Query scan id: {query.metadata['scans']}")
            print(f"Score: {score['score']:.4f}")
            print(f"Number of matching peaks: {score['matches']}")
            print("----------------------------")

Glossary of terms
=================

.. list-table::
   :header-rows: 1

   * - Term
     - Description
   * - adduct / addition product
     - During ionization in a mass spectrometer, the molecules of the injected compound break apart
       into fragments. When fragments combine into a new compound, this is known as an addition
       product, or adduct.  `Wikipedia <https://en.wikipedia.org/wiki/Adduct>`__
   * - GNPS
     - Knowledge base for sharing of mass spectrometry data (`link <https://gnps.ucsd.edu/ProteoSAFe/static/gnps-splash.jsp>`__).
   * - InChI / :code:`INCHI`
     - InChI is short for International Chemical Identifier. InChIs are useful
       in retrieving information associated with a certain molecule from a
       database.
   * - InChIKey / InChI key / :code:`INCHIKEY`
     - An identifier for molecules. For example, the InChI key for carbon
       dioxide is :code:`InChIKey=CURLTUGMZLYLDI-UHFFFAOYSA-N` (yes, it
       includes the substring :code:`InChIKey=`).
   * - MGF File / Mascot Generic Format
     - A plan ASCII file format to store peak list data from a mass spectrometry experiment. Links: `matrixscience.com <http://www.matrixscience.com/help/data_file_help.html#GEN>`__,
       `fiehnlab.ucdavis.edu <https://fiehnlab.ucdavis.edu/projects/lipidblast/mgf-files>`__.
   * - parent mass / :code:`parent_mass`
     - Actual mass (in Dalton) of the original compound prior to fragmentation.
       It can be recalculated from the precursor m/z by taking
       into account the charge state and proton/electron masses.
   * - precursor m/z / :code:`precursor_mz`
     - Mass-to-charge ratio of the compound targeted for fragmentation.
   * - SMILES
     - A line notation for describing the structure of chemical species using
       short ASCII strings. For example, water is encoded as :code:`O[H]O`,
       carbon dioxide is encoded as :code:`O=C=O`, etc. SMILES-encoded species may be converted to InChIKey `using a resolver like this one <https://cactus.nci.nih.gov/chemical/structure>`__. The Wikipedia entry for SMILES is `here <https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system>`__.


****************************
Documentation for developers
****************************

Installation
============

To install matchms, do:

.. code-block:: console

  git clone https://github.com/matchms/matchms.git
  cd matchms
  conda create --name matchms-dev python=3.8
  conda activate matchms-dev
  # Install rdkit using conda, rest of dependencies can be installed with pip
  conda install -c conda-forge rdkit
  python -m pip install --upgrade pip
  pip install --editable .[dev]

Run the linter with:

.. code-block:: console

  prospector

Automatically fix incorrectly sorted imports:

.. code-block:: console

  isort --recursive .

Files will be changed in place and need to be committed manually.

Run tests (including coverage) with:

.. code-block:: console

  pytest


Conda package
=============

The conda packaging is handled by a `recipe at Bioconda <https://github.com/bioconda/bioconda-recipes/blob/master/recipes/matchms/meta.yaml>`_.

Publishing to PyPI will trigger the creation of a `pull request on the bioconda recipes repository <https://github.com/bioconda/bioconda-recipes/pulls?q=is%3Apr+is%3Aopen+matchms>`_
Once the PR is merged the new version of matchms will appear on `https://anaconda.org/bioconda/matchms <https://anaconda.org/bioconda/matchms>`_

Flowchart
=========

.. figure:: paper/flowchart_matchms.png
  :width: 400
  :alt: Flowchart
  
  Flowchart of matchms workflow. Reference and query spectrums are filtered using the same
  set of set filters (here: filter A and filter B). Once filtered, every reference spectrum is compared to
  every query spectrum using the matchms.Scores object.

Contributing
============

If you want to contribute to the development of matchms,
have a look at the `contribution guidelines <CONTRIBUTING.md>`_.

*******
License
*******

Copyright (c) 2021, Netherlands eScience Center

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*******
Credits
*******

This package was created with `Cookiecutter
<https://github.com/audreyr/cookiecutter>`_ and the `NLeSC/python-template
<https://github.com/NLeSC/python-template>`_.
.. matchms documentation master file, created by
   sphinx-quickstart on Tue Apr  7 09:16:44 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to matchms's documentation!
===================================

Matchms is an open-access Python package to import, process, clean, and compare mass spectrometry data (MS/MS). It allows to implement and run an easy-to-follow, easy-to-reproduce workflow from raw mass spectra to pre- and post-processed spectral data.

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   API <api/matchms.rst>

Introduction
============

Matchms was designed to easily build custom spectra processing pipelines and to compute spectra similarities (see flowchart). Spectral data can be imported from common formats such mzML, mzXML, msp, metabolomics-USI, MGF, or json (e.g. GNPS-syle json files). Matchms then provides filters for metadata cleaning and checking, as well as for basic peak filtering. Finally, matchms was build to import and apply different similarity measures to compare large amounts of spectra. This includes common Cosine scores, but can also easily be extended by custom measures.

.. image:: _static/flowchart_matchms.png
  :width: 800
  :alt: matchms workflow illustration

Installation
============

Prerequisites:

- Python 3.7, 3.8, or 3.9
- Anaconda

Install matchms from Anaconda Cloud with

.. code-block:: console

  # install matchms in a new virtual environment to avoid dependency clashes
  conda create --name matchms python=3.8
  conda activate matchms
  conda install --channel nlesc --channel bioconda --channel conda-forge matchms

Example
=======

Below is a small example of using matchms to calculate the Cosine score between mass Spectrums in the `tests/pesticides.mgf <https://github.com/matchms/matchms/blob/master/tests/pesticides.mgf>`_ file.

.. testcode::

    from matchms.importing import load_from_mgf
    from matchms.filtering import default_filters
    from matchms.filtering import normalize_intensities
    from matchms import calculate_scores
    from matchms.similarity import CosineGreedy

    # Read spectrums from a MGF formatted file, for other formats see https://matchms.readthedocs.io/en/latest/api/matchms.importing.html
    file = load_from_mgf("../tests/pesticides.mgf")

    # Apply filters to clean and enhance each spectrum
    spectrums = []
    for spectrum in file:
        # Apply default filter to standardize ion mode, correct charge and more.
        # Default filter is fully explained at https://matchms.readthedocs.io/en/latest/api/matchms.filtering.html .
        spectrum = default_filters(spectrum)
        # Scale peak intensities to maximum of 1
        spectrum = normalize_intensities(spectrum)
        spectrums.append(spectrum)

    # Calculate Cosine similarity scores between all spectrums
    # For other similarity score methods see https://matchms.readthedocs.io/en/latest/api/matchms.similarity.html .
    # Because references and queries are here the same spectra, we can set is_symmetric=True
    scores = calculate_scores(references=spectrums,
                              queries=spectrums,
                              similarity_function=CosineGreedy(),
                              is_symmetric=True)

    # This computed all-vs-all similarity scores, the array of which can be accessed as scores.scores
    print(f"Size of matrix of computed similarities: {scores.scores.shape}")

    # Matchms allows to get the best matches for any query using scores_by_query
    query = spectrums[15]  # just an example
    best_matches = scores.scores_by_query(query, sort=True)

    # Print the calculated scores for each spectrum pair
    for (reference, score) in best_matches[:10]:
        # Ignore scores between same spectrum
        if reference is not query:
            print(f"Reference scan id: {reference.metadata['scans']}")
            print(f"Query scan id: {query.metadata['scans']}")
            print(f"Score: {score['score']:.4f}")
            print(f"Number of matching peaks: {score['matches']}")
            print("----------------------------")

Should output

.. testoutput::

    Size of matrix of computed similarities: (76, 76)
    Reference scan id: 613
    Query scan id: 2161
    Score: 0.8646
    Number of matching peaks: 14
    ----------------------------
    Reference scan id: 603
    Query scan id: 2161
    Score: 0.8237
    Number of matching peaks: 14
    ----------------------------
    ...

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
