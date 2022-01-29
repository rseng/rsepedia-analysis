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
reported by contacting the project team at j.m.hoch@uu.nl. All
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
# How to contribute

This python-package is a first outcome of an interdisciplinary project aimed at understanding the complex interplay between conflict and climate and environment.
As such, the presented code and functionalities can only be seen as a first step towards a fully-fledged model.
We therefore strongly encourage other users to contribute to this project!

## General notes

When contributing to this repository, please first discuss the change you wish to make via issue, email, or any other method with the owners of this repository before making a change.

Please note we have a code of conduct, please follow it in all your interactions with the project, the project owners, and users of the project.

## Getting Started

* Make sure you have a GitHub account.
* Fork the repository on GitHub.

## Making Changes

* Create a topic branch from where you want to base your work.
  * This is usually the dev branch.
  * Only target release branches if you are certain your fix must be on that
    branch.
  * To quickly create a topic branch based on master, run `git checkout -b
    fix/dev/my_contribution master`. Please avoid working directly on the
    `dev` (or `master`) branch.
* Make commits of logical and atomic units. Write a [good commit message][commit]!
* Make sure you have added the necessary tests for your changes.

## Submitting Changes

* Push your changes to a topic branch in your fork of the repository.
* Submit a pull request to the repository.
* The core team looks at pull requests as soon as possible, but no maximum waiting time can be given here.
* After feedback has been given we expect responses within two weeks. After two
  weeks we may close the pull request if it isn't showing any activity.

[commit]: http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html---
title: 'CoPro: a data-driven modelling framework for conflict risk projections'
tags:
  - Python
  - climate change
  - projections
  - conflict
  - climate security
  - water
  - risk
authors:
  - name: Jannis M. Hoch^[corresponding author]
    orcid: 0000-0003-3570-6436
    affiliation: 1
  - name: Sophie de Bruin
    orcid: 0000-0003-3429-349X
    affiliation: "1, 2"
  - name: Niko Wanders
    orcid: 0000-0002-7102-5454
    affiliation: 1
affiliations:
 - name: Department of Physical Geography, Utrecht University, Utrecht, the Netherlands
   index: 1
 - name: PBL Netherlands Environmental Assessment Agency, the Hague, the Netherlands
   index: 2
date: 18 February 2021
bibliography: bibliography.bib
---

# Summary

Climate change and environmental degradation are increasingly recognized as factors that can contribute to conflict risk under specific conditions.
In light of predicted shifts in climate patterns and the potentially resulting battle for increasingly scarce resources, it is widely acknowledged that there is an actual risk of increased armed conflict. To efficiently plan and implement adaptation and mitigation measures, it is key to first obtain an understanding of conflict drivers and spatial conflict risk distribution. And second, conflict risk needs to be projected to a given point in the future to be able to prepare accordingly. With CoPro, building and running models investigating the interplay between conflict and climate is made easier. By means of a clear workflow, maps of conflict risk for today as well as the future can be produced. Despite the structured workflow, CoPro caters for a variety of settings and input data, thereby capturing the multitude of facets of the climate-environment-conflict nexus.

# Statement of need 

There is increasing consensus that climate change can exacerbate the risk of (armed) conflict [@koubi2019climate; @mach2019climate]. Nevertheless, making (operational) projections of conflict risk is still challenging due to several reasons [@cederman2017predicting]. Building upon recent, similar approaches to use data-driven models [@colaresi2017robot] and statistical approaches [@witmer2017subnational; @hegre2016forecasting], CoPro is a novel, fully open, and extensible Python-model facilitating the set-up, execution, and evaluation of machine-learning models predicting conflict risk. CoPro provides a structured workflow including pre- and post-processing tools, making it accessible to all levels of experience. Such a user-friendly tool is needed not only to integrate the different disciplines, but also to extend the modeling approach with new insights and data - after all, the established links between climate and societal factors with conflict are still weak [@koubi2019climate; @mach2019climate]. In addition to scholarly explorations of the inter-dependencies and the importance of various conflict drivers, model output such as maps of spatially-disaggregated projected conflict risk can be an invaluable input to inform the decision-making process in affected regions.

Since conflicts are of all times and not limited to specific regions or countries, CoPro is designed with user-flexibility in mind. Therefore, the number and variables provided to the model is not specified, allowing for bespoke model designs. Depending on the modeling exercise and data used, several machine-learning models and pre-processing algorithms are available in CoPro. In its current form, the supervised learning techniques support vector classifier, k-neighbors classifier, and random-forest classifier are implemented. Catering for different model designs is of added value because of the non-linear and sometimes irrational - 'law-breaking' [@cederman2017predicting] - nature of conflicts. On top of that, the analyses can be run at any spatial scale, allowing for better identification of sub-national drivers of conflict risk. After all, conflict onset and conflicts are often limited to specific areas where driving factors coincide. 

Since the replicability of scientific results is important when developing forecast and projection models [@hegre2017introduction], CoPro produces reproducible output using transparent models. Hence, by making model code openly available and by including dedicated features in the model, we hope to advance the existing body of tools developed to project conflict risk.

These functionalities altogether make CoPro suited for both 'quick-and-dirty' and in-depth analyses of the relative importances of climate, environmental, and societal drivers as well as for assessments how conflict risk can change both in time and space.

# Acknowledgements
This research was supported by a Pathways to Sustainability Acceleration Grant from the Utrecht University.
We kindly acknowledge the valuable contributions from all partners at PBL, PRIO (Peace Research Institute Oslo), Uppsala University, and Utrecht University.

# References
