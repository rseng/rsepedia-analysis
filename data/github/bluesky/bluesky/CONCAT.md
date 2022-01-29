# Bluesky — An Experiment Specification & Orchestration Engine

[![Build Status](https://img.shields.io/github/workflow/status/bluesky/bluesky/Unit%20Tests)](https://github.com/bluesky/bluesky/actions?query=workflow%3A%22Unit+Tests%22+branch%3Amaster)
[![PyPI](https://img.shields.io/pypi/v/bluesky)](https://pypi.org/project/bluesky/)
[![conda-forge](https://img.shields.io/conda/vn/conda-forge/bluesky)](https://anaconda.org/conda-forge/bluesky)

The Bluesky Python Package is an experiment specification and orchestration engine. 
- Specify the logic of an experiment in a high-level, hardware-abstracted way.
- First-class support for adaptive feedback between analysis and acquisition.
- Data is emitted in a streaming fashion in standard Python data structures.
- Pause/resume, robust error handling, and rich metadata capture are built in.

[**Bluesky Documentation**](http://blueskyproject.io/bluesky).

The Bluesky Project enables experimental science at the lab-bench or facility scale. It is a collection of Python libraries that are co-developed but independently useful and may be adopted *a la carte*.

[**Bluesky Project Documentation**](http://blueskyproject.io).

<!--- Provide a general summary of the issue in the Title above -->

## Expected Behavior
<!--- If you're describing a bug, tell us what should happen -->
<!--- If you're suggesting a change/improvement, tell us how it should work -->

## Current Behavior
<!--- If describing a bug, tell us what happens instead of the expected behavior -->
<!--- If suggesting a change/improvement, explain the difference from current behavior -->

## Possible Solution
<!--- Not obligatory, but suggest a fix/reason for the bug, -->
<!--- or ideas how to implement the addition or change -->

## Steps to Reproduce (for bugs)
<!--- Provide a link to a live example, or an unambiguous set of steps to -->
<!--- reproduce this bug. Include code to reproduce, if relevant -->
1.
2.
3.

## Context
<!--- How has this issue affected you? What are you trying to accomplish? -->
<!--- Providing context helps us come up with a solution that is most useful in the real world -->

## Your Environment
<!--- Include as many relevant details about the environment you experienced the bug in -->
# Contributing

## Getting Started

* Make sure you have a [GitHub account](https://github.com/signup/free)
* Submit a ticket for your issue, assuming one does not already exist.
  * Clearly describe the issue including steps to reproduce when it is a bug.
  * Make sure you fill in the earliest version that you know has the issue.
* Fork the repository on GitHub


## Making Changes

* Create a topic branch from where you want to base your work.
  * This is usually the master branch.
  * Only target release branches if you are certain your fix must be on that
    branch.
  * To quickly create a topic branch based on master; `git checkout -b
    fix/master/my_contribution master`. Please avoid working directly on the
    `master` branch.
* Make commits of logical units.
* Check for unnecessary whitespace with `git diff --check` before committing.
* Make sure your commit messages are in the proper format (see below)
* Make sure you have added the necessary tests for your changes.
* Run _all_ the tests to assure nothing else was accidentally broken.

### Writing the commit message

Commit messages should be clear and follow a few basic rules. Example:

```
ENH: add functionality X to bluesky.<submodule>.

The first line of the commit message starts with a capitalized acronym
(options listed below) indicating what type of commit this is.  Then a blank
line, then more text if needed.  Lines shouldn't be longer than 72
characters.  If the commit is related to a ticket, indicate that with
"See #3456", "See ticket 3456", "Closes #3456" or similar.
```

Describing the motivation for a change, the nature of a bug for bug fixes 
or some details on what an enhancement does are also good to include in a 
commit message. Messages should be understandable without looking at the code 
changes. 

Standard acronyms to start the commit message with are:
```
API: an (incompatible) API change
BLD: change related to building numpy
BUG: bug fix
CI : continuous integration
DEP: deprecate something, or remove a deprecated object
DEV: development tool or utility
DOC: documentation
ENH: enhancement
MNT: maintenance commit (refactoring, typos, etc.)
REV: revert an earlier commit
STY: style fix (whitespace, PEP8)
TST: addition or modification of tests
REL: related to releases
```
## The Pull Request

* Now push to your fork
* Submit a [pull request](https://help.github.com/articles/using-pull-requests) to this branch. This is a start to the conversation.

At this point you're waiting on us. We like to at least comment on pull requests within three business days 
(and, typically, one business day). We may suggest some changes or improvements or alternatives.

Hints to make the integration of your changes easy (and happen faster):
- Keep your pull requests small
- Don't forget your unit tests
- All algorithms need documentation, don't forget the .rst file
- Don't take changes requests to change your code personally
<!--- Provide a general summary of your changes in the Title above -->

## Description
<!--- Describe your changes in detail -->

## Motivation and Context
<!--- Why is this change required? What problem does it solve? -->
<!--- If it fixes an open issue, please link to the issue here. -->

## How Has This Been Tested?
<!--- Please describe in detail how you tested your changes. -->
<!--- Include details of your testing environment, and the tests you ran to -->
<!--- see how your change affects other areas of the code, etc. -->

<!--
## Screenshots (if appropriate):
-->
