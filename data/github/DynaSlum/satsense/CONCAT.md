# Changelog

## version 0.9

- Added a large amount of documentation
  - Available at https://satsense.readthedocs.io/en/latest/
  - Includes:
    - Installation instructions
    - Example notebook for feature extraction
    - API Documentation and docstrings

- Bug fixes:
  - Histogram of Gradients
    - fixed the absolute sine difference calculation
  - Fixed the padding around the data when splitting the generator.
  - Fixed the generation of windows

- Development:
  - Added automated versioning
  - Increased code maintainability

## version 0.8
- Initial release
- Features included:
  - Histogram of Gradients
  - Pantex
  - NDVI
    - also available:
    - RgNDVI (Red-green based)
    - RbNDVI (Red-blue based)
    - NDSI (Snow Cover Index)
    - NDWI (Water Cover Index)
    - WVSI (Soil Cover Index)
  - Lacunarity
  - SIFT
  - TextonContributions are very welcome. Please make sure there is a github issue
associated with with every pull request. Creating an issue is also a good
way to propose/discuss new features or get help with using satsense.

# Installation for development

Please follow the installation instructions on
[readthedocs](https://satsense.readthedocs.io/en/latest/installation.html)
to get started.

# Testing

Please add unit tests for the code you are writing (e.g. when fixing a bug, implement
a test that demonstrates the bug is fixed). You can run the unit tests locally
with the command

```python
python setup.py test
```

# Coding style

Please make sure your code is formatted according to
[PEP8](https://www.python.org/dev/peps/pep-0008/) and docstrings are written
according to [PEP257](https://www.python.org/dev/peps/pep-0257/). Publicly visible
functions should have
[numpy style docstrings](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html).

Please autoformat your code with the following commands before making a pull request:

```bash
isort satsense/my_file.py
yapf -i satsense/my_file.py
```

Please use prospector to check that your code meets our standards:

```bash
prospector satsense/my_file.py
```

# Pull requests

Please create a pull request early, to keep other developers informed of what you're doing.
Limit the amount of work in a pull request to fixing a single bug or adding a single new feature.
Make sure the unit tests on Travis pass and review the comments by Codacy (click the Travis/Codacy
buttons below your pull request). Note that Codacy occasionally reports false positives, ask if in
doubt.

# Documentation

All public functions should have
[numpy style docstrings](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html).
You can build the documentation locally by running

```bash
python setup.py build_sphinx
```

Use

```bash
python setup.py build_sphinx -Ea
```

to build everying from scratch. Please check that there are no warnings.

## Converting Notebooks for documentation

If you update the notebooks please update their counterparts in the doc folder by using `jupyter nbconvert`

From the root of the project:
```bash
jupyter nbconvert --to rst notebooks/**/*.ipynb --output-dir=doc/notebooks/
```

# Creating a release

Make sure to update the version number and release date in CITATION.cff.
---
name: User Story
about: Create a scrum user story
title: ''
labels: ''
assignees: ''

---

**As a** <!-- user -->
**I want to** <!-- what -->
**because** <!-- why -->
# Notebooks accompanying the satsense Python package

## FeatureExtraction
Notebooks demonstrating how to extract features from satellite images using
satsense.

## Performance
Notebooks demonstrating how to use the performance metrics with satsense and
the utility functions to convert ground truth data to and from mask files and
shapefiles.

## Classification

Notebooks demonstrating classification using features calculated with Satsense.
