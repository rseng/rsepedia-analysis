# Cantera Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, caste, color, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
conduct@cantera.org.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
[https://www.contributor-covenant.org/version/2/0/code_of_conduct.html][v2.0].

Community Impact Guidelines were inspired by
[Mozilla's code of conduct enforcement ladder][Mozilla CoC].

For answers to common questions about this code of conduct, see the FAQ at
[https://www.contributor-covenant.org/faq][FAQ]. Translations are available
at [https://www.contributor-covenant.org/translations][translations].

[homepage]: https://www.contributor-covenant.org
[v2.0]: https://www.contributor-covenant.org/version/2/0/code_of_conduct.html
[Mozilla CoC]: https://github.com/mozilla/diversity
[FAQ]: https://www.contributor-covenant.org/faq
[translations]: https://www.contributor-covenant.org/translations
# Contributing to Cantera

* For significant changes, please start a discussion on the Cantera
  Users' Group or create an issue on the [Cantera/enhancements](https://github.com/Cantera/enhancements/issues/new/choose) repository
  on GitHub to plan your modifications so that they can be implemented
  efficiently and in a way that doesn't conflict with any other planned
  future development
* Fork the `Cantera/cantera` repository on Github
* Clone your new repository or add it as a remote to an existing repository
* Check out the existing `main` branch, then start a new feature branch for
  your work
* When making changes, write code that is consistent with the surrounding code
  (see the [style guidelines](#style-guidelines) below)
* Add tests for any new features that you are implementing to either the
  GoogleTest-based test suite or the Python test suite.
* Add examples that highlight new capabilities, or update existing
  examples to make use of new features.
* As you make changes, commit them to your feature branch
  * Configure Git with your name and e-mail address before making any commits
  * Use descriptive commit messages (summary line of no more than 72 characters,
    followed by a blank line and a more detailed summary, if any)
  * Make related changes in a single commit, and unrelated changes in separate
    commits
  * Make sure that your commits do not include any undesired files, e.g., files
    produced as part of the build process or other temporary files.
  * Use Git's history-rewriting features (i.e., `git rebase -i`; see
    https://help.github.com/articles/about-git-rebase/) to organize your commits
    and squash "fixup" commits and reversions.
  * Do not merge your branch with `main`. If needed, you should rebase your branch
    onto the most recent `HEAD` commit of `main`.
  * Periodically run the test suite (`scons test`) to make sure that your
    changes are not causing any test failures.
* Push the changes on your new feature branch to your forked copy of the
  `Cantera/cantera` repository on GitHub.

* Submit a Pull Request on Github, from your forked copy. Check the results
  of the continuous-integration tests run using GitHub Actions and resolve
  any issues that arise.
* Additional discussion of good Git & Github workflow is provided at
  http://matplotlib.org/devel/gitwash/development_workflow.html and
  https://docs.scipy.org/doc/numpy-1.15.0/dev/gitwash/development_workflow.html
* Cantera is licensed under a [BSD
  license](https://github.com/Cantera/cantera/blob/main/License.txt) which
  allows others to freely modify the code, and if your Pull Request is accepted,
  then that code will be release under this license as well. The copyright for
  Cantera is held collectively by the contributors. If you have made a
  significant contribution, please add your name to the `AUTHORS` file.

# Style Guidelines

* Try to follow the style of surrounding code, and use variable names that
  follow existing patterns. Pay attention to indentation and spacing.
* Configure your editor to use 4 spaces per indentation level, and **never to
  use tabs**.
* Avoid introducing trailing whitespace
* Limit line lengths to 88 characters when possible
* Write comments to explain non-obvious operations
* Use whitespaces to improve code readability (examples: after commas; before and
  after mathematical operators (`+`/`-`/`*`/`/` except `^`), binary operators
  (`&&`/`||`/...), and comparisons (`<`/`>`/`==`/...); before and after equality
  signs `=` unless used for the assignment of a default parameter)
* Do not go out of your way to change formatting in otherwise unmodified code

## C++

* All classes, member variables, and methods should have Doxygen-style comments
  (e.g., comment lines starting with `//!` or comment blocks starting with `/*!`)
* Avoid defining non-trivial functions in header files
* Header files should include an 'include guard'
* Protected and private member variable names are generally prefixed with
  `m_`. For most classes, member variables should not be public.
* Class names use `InitialCapsNames`
* Methods use `camelCaseNames`
* Do not indent the contents of namespaces
* Code should follow the C++11 standard, with minimum required compiler versions
  GCC 4.8, Clang 3.4, MSVC 14.0 (2015) and Intel 15.0.
* Avoid manual memory management (i.e. `new` and `delete`), preferring to use
  standard library containers, as well as `std::unique_ptr` and
  `std::shared_ptr` when dynamic allocation is required.
* Portions of Boost which are "header only" may be used. If possible, include
  Boost header files only within .cpp files rather than other header files to
  avoid unnecessary increases in compilation time. Boost should not be added
  to the public interface unless its existence and use is optional. This keeps
  the number of dependencies low for users of Cantera. In these cases,
  `CANTERA_API_NO_BOOST` should be used to conditionally remove Boost dependencies.
* While Cantera does not specifically follow these rules, the following style
  guides are useful references for possible style choices and the rationales behind them.
  * The Google C++ Style Guide: https://google.github.io/styleguide/cppguide.html
  * http://geosoft.no/development/cppstyle.html
* For any new code, do *not* use the `doublereal` and `integer` typedefs for the
  basic types `double` and `int`, but also do not go out of your way to change
  uses of these in otherwise unmodified code.

## Python

* Style generally follows PEP8 (https://www.python.org/dev/peps/pep-0008/)
* Code in `.py` and `.pyx` files needs to be written to work with Python 3
* The minimum Python version that Cantera supports is Python 3.6, so code should only use features added in Python 3.6 or earlier
* Please use double quotes in all new Python code
# How to get support

> This project has a [Code of Conduct](https://github.com/Cantera/cantera/blob/main/CODE_OF_CONDUCT.md).
> By interacting with this repository, organisation, or community you agree to
> abide by its terms.

For **help**, **support** and **questions** please create a post on the
**[Cantera Users' Group](https://groups.google.com/group/cantera-users)**.
Any discussion of Cantera functionality such as how to use certain function
calls, syntax problems, input files, etc. should be directed to the Users' Group.

Further, the **[Cantera Gitter Chat](https://gitter.im/Cantera/Lobby)** is an
infrequently monitored chat room that can be used to discuss tangentially-related
topics such as how to model the underlying physics of a problem, share cool
applications that you have developed, etc.

Please **_do not_** raise an issue on GitHub unless it is a bug report or a
feature request. Issues that do not fall into these categories will be closed.
If you're not sure, please make a post on the
[Users' Group](https://groups.google.com/group/cantera-users) and someone will
be able to help you out.

## Documentation

The [documentation](https://cantera.org/documentation)
offers a number of starting points:

- [Python tutorial](https://cantera.org/tutorials/python-tutorial.html)
- [Application Examples in Python (Jupyter)](https://github.com/Cantera/cantera-jupyter#cantera-jupyter)
- [A guide to Cantera's input file format](https://cantera.org/tutorials/input-files.html)
- [Information about the Cantera community](https://cantera.org/community.html)

Documentation for the [development version of
Cantera](https://cantera.org/documentation/dev-docs.html) is also available.

## Contributions

See [`CONTRIBUTING.md`](https://github.com/Cantera/cantera/blob/main/CONTRIBUTING.md) on how to contribute.
<!-- Thanks for contributing code! Please include a description of your change and check your pull request against the list below. For further questions, refer to the contributing guide (https://github.com/Cantera/cantera/blob/main/CONTRIBUTING.md). -->

**Changes proposed in this pull request**

<!-- Provide a clear and concise description of changes and/or features introduced in this pull request. -->

-
-
-

**If applicable, fill in the issue number this pull request is fixing**

<!-- Issues with issue number '<issue>' are referenced as #<issue>. To link to an issue in the enhancements repository, use Cantera/enhancements#<issue>. -->

Closes #

**If applicable, provide an example illustrating new features this pull request is introducing**

<!-- A minimal, complete, and reproducible example demonstrating features introduced by this pull request. See https://stackoverflow.com/help/minimal-reproducible-example for additional suggestions on how to create such an example. -->

**Checklist**

- [ ] The pull request includes a clear description of this code change
- [ ] Commit messages have short titles and reference relevant issues
- [ ] Build passes (`scons build` & `scons test`) and unit tests address code coverage
- [ ] Style & formatting of contributed code follows [contributing guidelines](https://github.com/Cantera/cantera/blob/main/CONTRIBUTING.md)
- [ ] The pull request is ready for review
---
name: Bug report
about: Report reproducible software issues so we can improve
title: ''
labels: ''
assignees: ''
---

<!-- Please fill in the following information to report a problem with Cantera. If you have a question about using Cantera, please post it on our Google Users' Group (https://groups.google.com/forum/#!forum/cantera-users). Feature enhancements should be discussed in the dedicated Cantera enhancements repository (https://github.com/Cantera/enhancements/new/choose) -->

**Problem description**

<!-- A clear and concise description of what the bug is. -->

**Steps to reproduce**

<!-- A minimal, complete, and reproducible example demonstrating the problem. See https://stackoverflow.com/help/minimal-reproducible-example for additional suggestions on how to create such an example. -->

1. Open '...'
2. Run '....'
3. See error '....'

**Behavior**

<!-- Describe the result of executing the above steps, and how this differs from what you expect to happen. -->

**System information**

- Cantera version: [for example, 2.5.0 or the git commit hash]
- OS: [for example, Windows 10]
- Python/MATLAB/other software versions:

**Attachments**

<!-- If applicable, attach scripts and/or input files to help explain your problem. Please do *not* attach screenshots of code or terminal output. -->

**Additional context**

<!-- Add any other context about the problem here. -->
---
name: Feature request
about: Suggest a new feature to enhance Cantera's capabilities
title: ''
labels: ''
assignees: ''
---

Feature requests have been moved to
[Cantera/enhancements](https://github.com/Cantera/enhancements/issues/new/choose) and should be
opened there. Thank you for your suggestions!
