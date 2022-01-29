
We welcome people who want to make contributions to Numba, big or small!
Even simple documentation improvements are encouraged.

# Asking questions

Numba has a [discourse forum](https://numba.discourse.group/) for longer/more
involved questions and an IRC channel on
[gitter.im](https://gitter.im/numba/numba) for quick questions and interactive
help.

# Ways to help:

There's lots of ways to help improve Numba, some of these require creating code
changes, see **contributing patches** below.

## Quick things:

* Answer a question asked on [discourse](https://numba.discourse.group/) or
  [gitter.im](https://gitter.im/numba/numba).
* Review a page of documentation, check it makes sense, that it's clear and
  still relevant, that the examples are present, good and working. Fix anything
  that needs updating in a pull request.
* Make a file that is not `flake8` compliant meet the standard, a list of all
  failing files is in the `exclude` section of the [`.flake8` config](https://github.com/numba/numba/blob/master/.flake8),
  then create a pull request with the change.

## More involved things:

* Review a pull request, you don't need to be a compiler engineer to do an
  initial review of a pull request. It's incredibly helpful to have pull
  requests go through a review to just make sure the code change is well formed,
  documented, efficient and clear. Further, if the code is fixing a bug, making
  sure that tests are present demonstrating it is fixed! Look out for PRs with
  the [`needs initial review`](https://github.com/numba/numba/labels/needs%20initial%20review)
  label.
* Work on fixing or implementing something in the code base, there are a lot of
  [`good first issue's`](https://github.com/numba/numba/labels/good%20first%20issue)
  and [`good second issue's`](https://github.com/numba/numba/labels/good%20first%20issue).
  For implementing new features/functionality, the extension API is the best
  thing to use and a guide to using `@overload` in particular is
  [here](https://numba.pydata.org/numba-doc/dev/extending/overloading-guide.html)
  and the API documentation is [here](https://numba.pydata.org/numba-doc/latest/extending/high-level.html#implementing-functions).

## Contributing patches

Please fork the Numba repository on Github, and create a new branch
containing your work.  When you are done, open a pull request.

# Further reading

Please read the [contributing guide](
https://numba.pydata.org/numba-doc/dev/developer/contributing.html).
<!--

Thanks for wanting to contribute to Numba :)

First, if you need some help or want to chat to the core developers, please
visit https://gitter.im/numba/numba for real time chat or post to the Numba
forum https://numba.discourse.group/.

Here's some guidelines to help the review process go smoothly.

0. Please write a description in this text box of the changes that are being
   made.

1. Please ensure that you have written units tests for the changes made/features
   added.

2. If you are closing an issue please use one of the automatic closing words as
   noted here: https://help.github.com/articles/closing-issues-using-keywords/

3. If your pull request is not ready for review but you want to make use of the
   continuous integration testing facilities here, please click the arrow besides
   "Create Pull Request" and choose "Create Draft Pull Request".
   When it's ready for review, you can click the button "ready to review" near
   the end of the pull request
   (besides "This pull request is still a work in progress".)
   The maintainers will then be automatically notified to review it.

4. Once review has taken place please do not add features or make changes out of
   the scope of those requested by the reviewer (doing this just add delays as
   already reviewed code ends up having to be re-reviewed/it is hard to tell
   what is new etc!). Further, please do not rebase your branch on master/force
   push/rewrite history, doing any of these causes the context of any comments
   made by reviewers to be lost. If conflicts occur against master they should
   be resolved by merging master into the branch used for making the pull
   request.

Many thanks in advance for your cooperation!

-->
---
name: Feature Request
about: Tell us about something in the Python language/NumPy you'd like Numba to support. Not for asking general questions - see below.

---

---

<!--

Thanks for opening an issue! To help the Numba team handle your information
efficiently, please first ensure that there is no other issue present that
already describes the issue you have
(search at https://github.com/numba/numba/issues?&q=is%3Aissue).

-->

## Feature request

<!--

Please include details of the feature you would like to see, why you would
like to see it/the use case.

-->
---
name: First Release Candidate Checklist (maintainer only)
about: Checklist template for the first release of every series
title: Numba X.Y.Zrc1 Checklist (FIXME)
labels: task

---


## Numba X.Y.Z

* [ ] Merge to master.
    - [ ] "remaining Pull-Requests from milestone".
* [ ] Review deprecation schedule and notices. Make PRs if need be.
* [ ] Merge change log changes.
    - [ ] "PR with changelog entries".
* [ ] Create X.Y release branch.
* [ ] Pin llvmlite to `>=0.A.0rc1,<0.A+1.0`.
* [ ] Pin NumPy if needed
* [ ] Pin tbb if needed
* [ ] Annotated tag X.Y.Zrc1 on release branch.
* [ ] Build and upload conda packages on buildfarm (check "upload").
* [ ] Build wheels (`$PYTHON_VERSIONS`) on the buildfarm.
* [ ] Verify packages uploaded to Anaconda Cloud and move to `numba/label/main`.
* [ ] Upload wheels and sdist to PyPI (upload from `ci_artifacts`).
* [ ] Verify wheels for all platforms arrived on PyPi.
* [ ] Initialize and verify ReadTheDocs build.
* [ ] Clean up `ci_artifacts`.
* [ ] Send RC announcement email / post announcement to discourse group.
* [ ] Post link to Twitter.

### Post Release:

* [ ] Tag X.Y+1.0dev0 to start new development cycle on `master`.
* [ ] Update llvmlite dependency spec to match next version via PR to `master`.
* [ ] Update release checklist template with any additional bullet points that
      may have arisen during the release.
* [ ] Close milestone (and then close this release issue).
---
name: Bug Report
about: Report a bug. Not for asking general questions - see below.

---

<!--

Thanks for opening an issue! To help the Numba team handle your information
efficiently, please first ensure that there is no other issue present that
already describes the issue you have
(search at https://github.com/numba/numba/issues?&q=is%3Aissue).

-->

## Reporting a bug

<!--

Before submitting a bug report please ensure that you can check off these boxes:

-->

- [ ] I have tried using the latest released version of Numba (most recent is
 visible in the change log (https://github.com/numba/numba/blob/master/CHANGE_LOG).
- [ ] I have included a self contained code sample to reproduce the problem.
  i.e. it's possible to run as 'python bug.py'.

<!--

Please include details of the bug here, including, if applicable, what you
expected to happen!

-->
---
name: Subsequent Release Candidate Checklist (maintainer only)
about: Checklist template for all subsequent releases (RC 2-N, FINAL and PATCH) of every series
title: Numba X.Y.Zrc1 Checklist (FIXME)
labels: task

---


## numba X.Y.Z

* [ ] Cherry-pick items from the X.Y.Z milestone into a PR.
* [ ] Approve change log modifications and cherry-pick.
* [ ] Merge change log modifications and cherry-picks to X.Y release branch.
  * [ ] https://github.com/numba/numba/pull/XXXX
* [ ] Review, merge and check execution of release notebook. (FINAL ONLY)
* [ ] Annotated tag X.Y.Z on release branch (no `v` prefix).
* [ ] Build and upload conda packages on buildfarm (check `upload`).
* [ ] Verify packages uploaded to Anaconda Cloud and move to
  `numba/label/main`.
* [ ] Build wheels (`$PYTHON_VERSIONS`) on the buildfarm.
* [ ] Upload wheels and sdist to PyPI (upload from `ci_artifacts`).
* [ ] Verify wheels for all platforms arrived on PyPi.
* [ ] Verify ReadTheDocs build.
* [ ] Send RC/FINAL announcement email / post announcement to discourse group.
* [ ] Post link to Twitter.
* [ ] Post link to python-announce-list@python.org.

### Post release

* [ ] Clean up `ci_artifacts` by moving files to subdirectories
* [ ] Update release checklist template.
* [ ] Ping Anaconda Distro team to trigger a build for `defaults` (FINAL ONLY).
* [ ] Create a release on Github at https://github.com/numba/numba/releases (FINAL ONLY).
* [ ] Close milestone (and then close this release issue).
# DAG Roadmap

This directory includes a representation of the Numba roadmap in the form of a
DAG.  We have done this to enable a highly granular display of enhancements to
Numba that also shows the relationships between these tasks. Many tasks have
prerequisites, and we've found that issue trackers, Kanban boards, and
time-bucketed roadmap documentation all fail to represent this information in
different ways.

## Requirements

```
conda install jinja2 python-graphviz pyyaml
```

## Usage

```
./render.py -o dagmap.html dagmap.yaml
```

The generated HTML file will look for `jquery.graphviz.svg.js` in the same
directory.

## Updating the DAG

Copy one of the existing tasks and edit:
  * `label`: text appears on the node.  Embed `\n` for line breaks.
  * `id`: Referenced to indicate a dependency
  * `description`: Shown in the tooltip.  Automatically word-wrapped.
  * `depends_on`: Optional list of task IDs which this task depends on.

The `style` section of the file is not used yet.

## Notes

The HTML rendering of the graph is based on a slightly modified version of
(jquery.graphviz.svg)[https://github.com/mountainstorm/jquery.graphviz.svg/].
Its license is:
```
Copyright (c) 2015 Mountainstorm
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```