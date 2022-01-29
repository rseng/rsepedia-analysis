# Security Policy

## Supported Versions

The following table lists versions and whether they are supported. Security
vulnerability reports will be accepted and acted upon for all supported
versions.

| Version | Supported          |
| ------- | ------------------ |
| 3.5.x   | :white_check_mark: |
| 3.4.x   | :white_check_mark: |
| 3.3.x   | :x:                |
| < 3.3   | :x:                |


## Reporting a Vulnerability

If you have found a security vulnerability, in order to keep it confidential,
please do not report an issue on GitHub.

Please email us details of the vulnerability at matplotlib@numfocus.org;
include a description and proof-of-concept that is [short and
self-contained](http://www.sscce.org/).

You should expect a response within a week of your email. Depending on the
severity of the issue, this may require some time to draft an immediate bugfix
release. Less severe issues may be held until the next release.

We do not award bounties for security vulnerabilities.

# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
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
matplotlib@numfocus.org (monitored by Thomas Caswell (\@tacaswell) and Hannah
Aizenman (\@story645)) or a report can be made using the [NumFOCUS Code of Conduct report form][numfocus form]. 
If community leaders cannot come to a resolution about enforcement, reports will be escalated to the NumFocus Code of Conduct committee (conduct@numfocus.org).
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

[numfocus form]: https://numfocus.typeform.com/to/ynjGdT

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
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
# Test project for matplotlib sphinx extensions

A tiny sphinx project from ``sphinx-quickstart`` with all default answers.
Please refer to the [developers guide](https://matplotlib.org/devel/index.html).
## PR Summary

## PR Checklist

<!-- Please mark any checkboxes that do not apply to this PR as [N/A]. -->
**Tests and Styling**
- [ ] Has pytest style unit tests (and `pytest` passes).
- [ ] Is [Flake 8](https://flake8.pycqa.org/en/latest/) compliant (install `flake8-docstrings` and run `flake8 --docstring-convention=all`).

**Documentation**
- [ ] New features are documented, with examples if plot related.
- [ ] New features have an entry in `doc/users/next_whats_new/` (follow instructions in README.rst there).
- [ ] API changes documented in `doc/api/next_api_changes/` (follow instructions in README.rst there).
- [ ] Documentation is sphinx and numpydoc compliant (the docs should [build](https://matplotlib.org/devel/documenting_mpl.html#building-the-docs) without error).

<!--
Thank you so much for your PR!  To help us review your contribution, please
consider the following points:

- A development guide is available at https://matplotlib.org/devdocs/devel/index.html.

- Help with git and github is available at
  https://matplotlib.org/devel/gitwash/development_workflow.html.

- Do not create the PR out of main, but out of a separate branch.

- The PR title should summarize the changes, for example "Raise ValueError on
  non-numeric input to set_xlim".  Avoid non-descriptive titles such as
  "Addresses issue #8576".

- The summary should provide at least 1-2 sentences describing the pull request
  in detail (Why is this change required?  What problem does it solve?) and
  link to any relevant issues.

- If you are contributing fixes to docstrings, please pay attention to
  http://matplotlib.org/devel/documenting_mpl.html#formatting.  In particular,
  note the difference between using single backquotes, double backquotes, and
  asterisks in the markup.

We understand that PRs can sometimes be overwhelming, especially as the
reviews start coming in.  Please let us know if the reviews are unclear or
the recommended next step seems overly demanding, if you would like help in
addressing a reviewer's comments, or if you have been waiting too long to hear
back on your PR.
-->
