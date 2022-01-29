SWIFTsimIO
==========

![Build Status](https://github.com/swiftsim/swiftsimio/actions/workflows/pytest.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/swiftsimio/badge/?version=latest)](https://swiftsimio.readthedocs.io/en/latest/?badge=latest)
[![JOSS Status](https://joss.theoj.org/papers/e85c85f49b99389d98f9b6d81f090331/status.svg)](https://joss.theoj.org/papers/e85c85f49b99389d98f9b6d81f090331)


The SWIFT astrophysical simulation code (http://swift.dur.ac.uk) is used
widely. There exists many ways of reading the data from SWIFT, which outputs
HDF5 files. These range from reading directly using `h5py` to using a complex
system such as `yt`; however these either are unsatisfactory (e.g. a lack of
unit information in reading HDF5), or too complex for most use-cases. 
`swiftsimio` provides an object-oriented API to read (dynamically) data
from SWIFT.

Full documentation is available at [ReadTheDocs](http://swiftsimio.readthedocs.org).

Getting set up with `swiftsimio` is easy; it (by design) has very few
requirements. There are a number of optional packages that you can install
to make the experience better and these are recommended.


Requirements
------------

This requires `python` `v3.6.0` or higher. Unfortunately it is not
possible to support `swiftsimio` on versions of python lower than this.
It is important that you upgrade if you are still a `python2` user.

### Python packages


+ `numpy`, required for the core numerical routines.
+ `h5py`, required to read data from the SWIFT HDF5 output files.
+ `unyt`, required for symbolic unit calculations (depends on sympy`).

### Optional packages


+ `numba`, highly recommended should you wish to use the in-built visualisation
  tools.
+ `scipy`, required if you wish to generate smoothing lengths for particle types
  that do not store this variable in the snapshots (e.g. dark matter)
+ `tqdm`, required for progress bars for some long-running tasks. If not installed
  no progress bar will be shown.
+ `py-sphviewer`, if you wish to use our integration with this visualisation
  code.


Installing
----------

`swiftsimio` can be installed using the python packaging manager, `pip`,
or any other packaging manager that you wish to use:

`pip install swiftsimio`


Citing
------

Please cite `swiftsimio` using the JOSS [paper](https://joss.theoj.org/papers/10.21105/joss.02430):

```bibtex
@article{Borrow2020,
  doi = {10.21105/joss.02430},
  url = {https://doi.org/10.21105/joss.02430},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {52},
  pages = {2430},
  author = {Josh Borrow and Alexei Borrisov},
  title = {swiftsimio: A Python library for reading SWIFT data},
  journal = {Journal of Open Source Software}
}
```

If you use any of the subsampled projection backends, we ask that you cite our
relevant SPHERIC [paper](https://arxiv.org/abs/2106.05281). Note that citing
the arXiv version here is recommended as the ADS cannot track conference
proceedings well.

```bibtex
@article{Borrow2021
  title={Projecting SPH Particles in Adaptive Environments}, 
  author={Josh Borrow and Ashley J. Kelly},
  year={2021},
  eprint={2106.05281},
  archivePrefix={arXiv},
  primaryClass={astro-ph.GA}
}
```
Contributing to SWIFTsimIO
==========================

Contributions for SWIFTsimIO should come through our GitHub repository,
available at https://github.com/swiftsim/swiftsimio.

Contributions are always welcome, but you should make sure of the following:

+ Your contributions pass all unit tests (you can check this with `pytest`)
+ Your contributions add unit tests for new functionality
+ Your contributions are formatted with the `black` formatter (see `format.sh`)
+ Your contributions are documented fully under `/docs`.

You should also abide by the following code of conduct:

### Code of Conduct

The community of participants in open source Astronomy projects is made up of
members from around the globe with a diverse set of skills, personalities,
and experiences. It is through these differences that our community
experiences success and continued growth. We expect everyone in our community
to follow these guidelines when interacting with others both inside and
outside of our community. Our goal is to keep ours a positive, inclusive,
successful, and growing community.

As members of the community,

+ We pledge to treat all people with respect and provide a harassment- and
  bullying-free environment, regardless of sex, sexual orientation and/or
  gender identity, disability, physical appearance, body size, race,
  nationality, ethnicity, and religion. In particular, sexual language and
  imagery, sexist, racist, or otherwise exclusionary jokes are not appropriate.
+ We pledge to respect the work of others by recognizing
  acknowledgement/citation requests of original authors. As authors, we pledge
  to be explicit about how we want our own work to be cited or acknowledged.
+ We pledge to welcome those interested in joining the community, and realize
  that including people with a variety of opinions and backgrounds will only
  serve to enrich our community. In particular, discussions relating to
  pros/cons of various technologies, programming languages, and so on are
  welcome, but these should be done with respect, taking proactive measure to
  ensure that all participants are heard and feel confident that they can
  freely express their opinions.
+ We pledge to welcome questions and answer them respectfully, paying
  particular attention to those new to the community. We pledge to provide
  respectful criticisms and feedback in forums, especially in discussion
  threads resulting from code contributions.
+ We pledge to be conscientious of the perceptions of the wider community and
  to respond to criticism respectfully. We will strive to model behaviours that
  encourage productive debate and disagreement, both within our community and
  where we are criticized. We will treat those outside our community with the
  same respect as people within our community.
+ We pledge to help the entire community follow the code of conduct, and to
  not remain silent when we see violations of the code of conduct. We will
  take action when members of our community violate this code such as
  contacting joshua.borrow@durham.ac.uk with the subject line SWIFTsimIO Code
  of Conduct (all emails sent in this fashion will be treated with the
  strictest confidence) or talking privately with the person.
+ This code of conduct applies to all community situations online and
  offline, including mailing lists, forums, social media, conferences,
  meetings, associated social events, and one-to-one interactions.

Any related activity or project organized by members of the SWIFTsimIO
community, including affiliated packages, are welcome to have their own codes
of conduct, but agree to also abide by the present code of conduct.

Parts of this code of conduct have been adapted from the PSF code of conduct and
the Astropy code of conduct: https://www.astropy.org/code_of_conduct.html.