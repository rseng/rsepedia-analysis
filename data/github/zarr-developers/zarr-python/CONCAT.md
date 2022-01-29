<div align="center">
  <img src="https://raw.githubusercontent.com/zarr-developers/community/master/logos/logo2.png"><br>
</div>

# Zarr

<table>
<tr>
  <td>Latest Release</td>
  <td>
    <a href="https://pypi.org/project/zarr/">
    <img src="https://badge.fury.io/py/zarr.svg" alt="latest release" />
    </a>
  </td>
</tr>
  <td></td>
  <td>
    <a href="https://anaconda.org/anaconda/zarr/">
    <img src="https://anaconda.org/conda-forge/zarr/badges/version.svg" alt="latest release" />
    </a>
</td>
</tr>
<tr>
  <td>Package Status</td>
  <td>
		<a href="https://pypi.org/project/zarr/">
		<img src="https://img.shields.io/pypi/status/zarr.svg" alt="status" />
		</a>
  </td>
</tr>
<tr>
  <td>License</td>
  <td>
    <a href="https://github.com/zarr-developers/zarr-python/blob/master/LICENSE">
    <img src="https://img.shields.io/pypi/l/zarr.svg" alt="license" />
    </a>
</td>
</tr>
<tr>
  <td>Build Status</td>
  <td>
    <a href="https://travis-ci.org/zarr-developers/zarr-python">
    <img src="https://travis-ci.org/zarr-developers/zarr-python.svg?branch=master" alt="travis build status" />
    </a>
  </td>
</tr>
<tr>
  <td>Coverage</td>
  <td>
    <a href="https://codecov.io/gh/zarr-developers/zarr-python">
    <img src="https://codecov.io/gh/zarr-developers/zarr-python/branch/master/graph/badge.svg"/ alt="coverage">
    </a>
  </td>
</tr>
<tr>
  <td>Downloads</td>
  <td>
    <a href="https://zarr.readthedocs.io">
    <img src="https://pepy.tech/badge/zarr" alt="pypi downloads" />
    </a>
  </td>
</tr>
<tr>
	<td>Gitter</td>
	<td>
		<a href="https://gitter.im/zarr-developers/community">
		<img src="https://badges.gitter.im/zarr-developers/community.svg" />
		</a>
	</td>
</tr>
<tr>
	<td>Citation</td>
	<td>
		<a href="https://doi.org/10.5281/zenodo.3773450">
			<img src="https://zenodo.org/badge/DOI/10.5281/zenodo.3773450.svg" alt="DOI">
		</a>
	</td>
</tr>

</table>

## What is it?

Zarr is a Python package providing an implementation of compressed, chunked, N-dimensional arrays, designed for use in parallel computing. See the [documentation](https://zarr.readthedocs.io) for more information.

## Main Features

- [**Create**](https://zarr.readthedocs.io/en/stable/tutorial.html#creating-an-array) N-dimensional arrays with any NumPy `dtype`.
- [**Chunk arrays**](https://zarr.readthedocs.io/en/stable/tutorial.html#chunk-optimizations) along any dimension.
- [**Compress**](https://zarr.readthedocs.io/en/stable/tutorial.html#compressors) and/or filter chunks using any NumCodecs codec.
- [**Store arrays**](https://zarr.readthedocs.io/en/stable/tutorial.html#tutorial-storage) in memory, on disk, inside a zip file, on S3, etc...
- [**Read**](https://zarr.readthedocs.io/en/stable/tutorial.html#reading-and-writing-data) an array [**concurrently**](https://zarr.readthedocs.io/en/stable/tutorial.html#parallel-computing-and-synchronization) from multiple threads or processes.
- Write to an array concurrently from multiple threads or processes.
- Organize arrays into hierarchies via [**groups**](https://zarr.readthedocs.io/en/stable/tutorial.html#groups).

## Where to get it

Zarr can be installed from PyPI using `pip`:

```bash
pip install zarr
```

or via `conda`:

```bash
conda install -c conda-forge zarr
```

For more details, including how to install from source, see the [installation documentation](https://zarr.readthedocs.io/en/stable/#installation).
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, gender identity and expression, level of experience, nationality, personal appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable behavior and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviors that they deem inappropriate, threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. Examples of representing a project or community include using an official project e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event. Representation of a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at zarr.conduct@gmail.com. The project team will review and investigate all complaints, and will respond in a way that it deems appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
For bug reports, please follow the template below. For enhancement proposals, feel free
to use whatever template makes sense (major new features should be discussed in the
Zarr specifications repository https://github.com/zarr-developers/zarr-specs).

#### Minimal, reproducible code sample, a copy-pastable example if possible

```python
# Your code here

```

#### Problem description

Explain why the current behavior is a problem, what the expected output/behaviour 
is, and why the expected output/behaviour is a better solution.

#### Version and installation information

Please provide the following:

* Value of ``zarr.__version__``
* Value of ``numcodecs.__version__``
* Version of Python interpreter
* Operating system (Linux/Windows/Mac)
* How Zarr was installed (e.g., "using pip into virtual environment", or "using conda")

Also, if you think it might be relevant, please provide the output from ``pip freeze`` or
``conda env export`` depending on which was used to install Zarr.
Contributing
============

Please see the [project documentation](http://zarr.readthedocs.io/en/stable/contributing.html) for information about contributing to Zarr.

[Description of PR]

TODO:
* [ ] Add unit tests and/or doctests in docstrings
* [ ] Add docstrings and API docs for any new/modified user-facing classes and functions
* [ ] New/modified features documented in docs/tutorial.rst
* [ ] Changes documented in docs/release.rst
* [ ] GitHub Actions have all passed
* [ ] Test coverage is 100% (Codecov passes)
