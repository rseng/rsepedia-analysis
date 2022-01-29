# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## Added

- Now supports Python 3.9 (including CI test runs) [#40](https://github.com/iomega/spec2vec/issues/40)

## [0.5.0] - 2021-06-18

## Changed

- Spec2Vec is now using gensim >= 4.0.0 [#62](https://github.com/iomega/spec2vec/pull/62)

## [0.4.0] - 2021-02-10

## Changed

- refactored `Spec2Vec` to now accept `Spectrum` or `SpectrumDocument` as input [#51](https://github.com/iomega/spec2vec/issues/51)

## Fixed

- updated and fixed code examples  [#51](https://github.com/iomega/spec2vec/issues/51)
- updated and fixed attribute typing [#51](https://github.com/iomega/spec2vec/issues/51)

## [0.3.4] - 2021-02-10

### Changed

- update required numba version to >=0.51 to avoid issues between numba and numpy [#55](https://github.com/iomega/spec2vec/pull/55)

## [0.3.3] - 2021-02-09

### Added

- Metadata getter method for `SpectrumDocument` [#50](https://github.com/iomega/spec2vec/pull/50)
- Implement `is_symmetric=True` option for `Spec2Vec.matrix` method [#53](https://github.com/iomega/spec2vec/pull/53)

### Changed

- Change default for `n_decimals` parameter from 1 to 2 [#50](https://github.com/iomega/spec2vec/pull/50)

## [0.3.2] - 2020-12-03

### Changed

- Add optional progress bar for spec2vec.matrix() calculations (default is False) [#43](https://github.com/iomega/spec2vec/pull/43)

## [0.3.1] - 2020-09-23

### Changed

- Implement faster, numba-based cosine similarity function [#29](https://github.com/iomega/spec2vec/pull/29)

## [0.3.0] - 2020-09-16

### Added

- Support for Python 3.8 [#35](https://github.com/iomega/spec2vec/pull/35)

### Changed

- Refactored Spec2Vec class to provide .pair() and .matrix() methods [#35](https://github.com/iomega/spec2vec/pull/35)

### Removed

- Spec2VecParallel (is now included as Spec2Vec.matrix()) [#35](https://github.com/iomega/spec2vec/pull/35)

## [0.2.0] - 2020-06-18

### Added

- Wrapper for training a gensim word2vec model [#13](https://github.com/iomega/spec2vec/tree/13-gensim-wrapper)
- Basic logger for word2vec model training [#11](https://github.com/iomega/spec2vec/issues/11)

### Changed

- Extend spec2vec similarity calculation to handle missing words [#9](https://github.com/iomega/spec2vec/issues/9)
- Extend documentation and given code examples [#15](https://github.com/iomega/spec2vec/issues/15)
- Updated the integration test to work with matchms 0.4.0 [#7](https://github.com/iomega/spec2vec/issues/7)

## [0.1.0] - 2020-06-02

### Added

- Matchms as dependency [#4](https://github.com/iomega/spec2vec/pull/4)
- Bump2version config

### Changed

- Splitted spec2vec from [matchms]. See (https://github.com/matchms/matchms) [#1](https://github.com/iomega/spec2vec/pull/1) [#4](https://github.com/iomega/spec2vec/pull/4)
  - Updated packaging related configuration
  - Update the GH Actions workflows
  - Updated the documentation
  - Updated the badges
  - Updated the integration and unit tests
  - Zenodo metadata
  
### Fixed

### Removed

- Fossa configuration
- Flowchart

[Unreleased]: https://github.com/iomega/spec2vec/compare/0.5.0...HEAD
[0.5.0]: https://github.com/iomega/spec2vec/compare/0.4.0...0.5.0
[0.4.0]: https://github.com/iomega/spec2vec/compare/0.3.4...0.4.0
[0.3.4]: https://github.com/iomega/spec2vec/compare/0.3.3...0.3.4
[0.3.3]: https://github.com/iomega/spec2vec/compare/0.3.2...0.3.3
[0.3.2]: https://github.com/iomega/spec2vec/compare/0.3.1...0.3.2
[0.3.1]: https://github.com/iomega/spec2vec/compare/0.3.0...0.3.1
[0.3.0]: https://github.com/iomega/spec2vec/compare/0.2.0...0.3.0
[0.2.0]: https://github.com/iomega/spec2vec/compare/0.1.0...0.2.0
[0.1.0]: https://github.com/iomega/spec2vec/releases/tag/0.1.0
# Contributing guidelines

We welcome any kind of contribution to our software, from simple comment or question to a full fledged [pull request](https://help.github.com/articles/about-pull-requests/). Please read and follow our [Code of Conduct](CODE_OF_CONDUCT.rst).

A contribution can be one of the following cases:

1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation);
1. you want to make a new release of the code base.

The sections below outline the steps in each case.

## You have a question

1. use the search functionality [here](https://github.com/iomega/spec2vec/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue;
1. apply the "Question" label; apply other labels when relevant.

## You think you may have found a bug

1. use the search functionality [here](https://github.com/iomega/spec2vec/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the [SHA hashcode](https://help.github.com/articles/autolinked-references-and-urls/#commit-shas) of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
1. apply relevant labels to the newly created issue.

## You want to make some kind of change to the code base

1. (**important**) announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
1. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/) and [here](https://help.github.com/articles/syncing-a-fork/));
1. make sure the existing tests still work by running ``python setup.py test``;
1. add your own tests (if necessary);
1. update or expand the documentation;
1. update the `CHANGELOG.md` file with change;
1. [push](http://rogerdudler.github.io/git-guide/>) your feature branch to (your fork of) the spec2vec repository on GitHub;
1. create the pull request, e.g. following the instructions [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.

## You want to make a new release of the code base

To create release you need write permission on the repository.

1. Check author list in `citation.cff` and `.zenodo.json` files
1. Bump the version using `bump2version <major|minor|patch>`. For example, `bump2version major` will increase major version numbers everywhere its needed (code, meta, etc.) in the repo.
1. Update the `CHANGELOG.md` to include changes made
1. Goto [GitHub release page](https://github.com/iomega/spec2vec/releases)
1. Press draft a new release button
1. Fill version, title and description field
1. Press the Publish Release button

A GitHub action will run which will publish the new version to [anaconda](https://anaconda.org/nlesc/spec2vec).
Also a Zenodo entry will be made for the release with its own DOI.
