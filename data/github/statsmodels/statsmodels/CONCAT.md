Release Notes
=============

The list of changes for each statsmodels release can be found [here](https://www.statsmodels.org/devel/release/index.html). Full details are available in the [commit logs](https://github.com/statsmodels/statsmodels).
- [ ] closes #xxxx
- [ ] tests added / passed. 
- [ ] code/documentation is well formatted.  
- [ ] properly formatted commit message. See 
      [NumPy's guide](https://docs.scipy.org/doc/numpy-1.15.1/dev/gitwash/development_workflow.html#writing-the-commit-message). 

<details>


**Notes**:

* It is essential that you add a test when making code changes. Tests are not 
  needed for doc changes.
* When adding a new function, test values should usually be verified in another package (e.g., R/SAS/Stata).
* When fixing a bug, you must add a test that would produce the bug in main and
  then show that it is fixed with the new code.
* New code additions must be well formatted. Changes should pass flake8. If on Linux or OSX, you can
  verify you changes are well formatted by running 
  ```
  git diff upstream/main -u -- "*.py" | flake8 --diff --isolated
  ```
  assuming `flake8` is installed. This command is also available on Windows 
  using the Windows System for Linux once `flake8` is installed in the 
  local Linux environment. While passing this test is not required, it is good practice and it help 
  improve code quality in `statsmodels`.
* Docstring additions must render correctly, including escapes and LaTeX.

</details>
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

#### Describe the bug

[A clear and concise description of what the bug is. This should explain **why** the current behaviour is a problem and why the expected output is a better solution.]

#### Code Sample, a copy-pastable example if possible


```python
# Your code here that produces the bug
# This example should be self-contained, and so not rely on external data.
# It should run in a fresh ipython session, and so include all relevant imports.
```
<details>

**Note**: As you can see, there are many issues on our GitHub tracker, so it is very possible that your issue has been posted before. Please check first before submitting so that we do not have to handle and close duplicates.

**Note**: Please be sure you are using the latest released version of `statsmodels`, or a recent build of `main`. If your problem has been fixed in an unreleased version, you might be able to use `main` until a new release occurs. 

**Note**: If you are using a released version, have you verified that the bug exists in the main branch of this repository? It helps the limited resources if we know problems exist in the current main branch so that they do not need to check whether the code sample produces a bug in the next release.

</details>


If the issue has not been resolved, please file it in the issue tracker.

#### Expected Output

A clear and concise description of what you expected to happen.

#### Output of ``import statsmodels.api as sm; sm.show_versions()``

<details>

[paste the output of ``import statsmodels.api as sm; sm.show_versions()`` here below this line]

</details>
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

#### Is your feature request related to a problem? Please describe
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

#### Describe the solution you'd like
A clear and concise description of what you want to happen.

#### Describe alternatives you have considered
A clear and concise description of any alternative solutions or features you have considered.

#### Additional context
Add any other context about the feature request here.This directory holds files that were once part of statsmodels but
are no longer maintained.  They are retained here in order to have their
git histories readily available, but should *not* be considered usable.
# Continuous Integration Tools

These scripts are used to implement Continuous Integration on Travis, AppVeyor
and Azure.  They should not be removed without careful consideration of the 
consequences.
# Documentation Documentation

We use a combination of sphinx and Jupyter notebooks for the documentation.
Jupyter notebooks should be used for longer, self-contained examples demonstrating
a topic.
Sphinx is nice because we get the tables of contents and API documentation.

## Build Process

Building the docs requires a few additional dependencies. You can get most
of these with

```bash

   pip install -e .[docs]

```

From the root of the project.
Some of the examples rely on `rpy2` to execute R code from the notebooks.
It's not included in the setup requires since it's known to be difficult to
install.

To generate the HTML docs, run ``make html`` from the ``docs`` directory.
This executes a few distinct builds

1. datasets
2. notebooks
3. sphinx

# Notebook Builds

We're using `nbconvert` to execute the notebooks, and then convert them
to HTML. The conversion is handled by `statsmodels/tools/nbgenerate.py`.
The default python kernel (embedded in the notebook) is `python3`.
You need at least `nbconvert==4.2.0` to specify a non-default kernel,
which can be passed in the Makefile.
The code in this folder is based on the Federal Reserve Bank of New York code
found at https://github.com/FRBNY-TimeSeriesAnalysis/Nowcasting, which was
downloaded as of commit 19f365cab8269e3aac3faa11ad091d6e913c5c43. Only the
files from that repository which were required for generating the test results
are included here.

In additionm the following files from the original package have been modified
(use git diff against the above repository to see the changes)

- functions/dfm.m
- functions/update_nowcast.m

The following files are not a part of the original package:

- test_DFM_blocks.m
- test_DFM.m
- test_news_blocks.m
- test_news.m
- test_spec_blocks.xls
- test_spec.xls
