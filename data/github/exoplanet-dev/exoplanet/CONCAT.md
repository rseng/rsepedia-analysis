<p align="center">
  <img width="240" src="https://raw.githubusercontent.com/exoplanet-dev/exoplanet/main/docs/_static/logo.png">
  <br><br>
  <a href="https://github.com/exoplanet-dev/exoplanet/actions/workflows/tests.yml">
    <img src="https://github.com/exoplanet-dev/exoplanet/actions/workflows/tests.yml/badge.svg" alt="Tests">
  </a>
  <a href="https://docs.exoplanet.codes">
    <img src="https://readthedocs.org/projects/exoplanet/badge/?version=latest" alt="Docs">
  </a>
  <a href="https://coveralls.io/github/exoplanet-dev/exoplanet?branch=main">
    <img src="https://coveralls.io/repos/github/exoplanet-dev/exoplanet/badge.svg?branch=main" alt="Coverage">
  </a>
</p>

# exoplanet

Fast & scalable MCMC for all your exoplanet needs! _exoplanet_ is a toolkit for
probabilistic modeling of time series data in astronomy with a focus on
observations of [exoplanets](https://en.wikipedia.org/wiki/Exoplanet), using
[PyMC3](https://docs.pymc.io). _PyMC3_ is a flexible and high-performance model
building language and inference engine that scales well to problems with a large
number of parameters. _exoplanet_ extends _PyMC3_'s language to support many of
the custom functions and distributions required when fitting exoplanet datasets.

Read the full documentation at [docs.exoplanet.codes](https://docs.exoplanet.codes).

## Installation

The quickest way to get started is to use [pip](https://pip.pypa.io):

```bash
python -m pip install exoplanet
```

Note that you will need Python (>=3.6) installed for this to work, but then this
will install all the required dependencies.

Check out the [main installation documentation
page](https://docs.exoplanet.codes/en/latest/user/install/) for more options.

## Usage

Check out the tutorials and API docs on [the docs
page](https://docs.exoplanet.codes) for example usage and much more info. You
can also find more in-depth examples on the [exoplanet case studies
page](https://gallery.exoplanet.codes).

## Contributing

_exoplanet_ is an open source project and we would love it if you wanted to
contribute. Check out [the developer
documentation](https://docs.exoplanet.codes/en/latest/user/dev/) for more info
about getting started.
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
reported by contacting the project team at foreman.mackey@gmail.com. All
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
## Contributing to exoplanet

### Reporting issues

If you find a bug or other unexpected behavior while using `exoplanet`, open an issue
on the [GitHub repository](https://github.com/exoplanet-dev/exoplanet/issues) and we
will try to respond and (hopefully) solve the problem in a timely manner. Similarly,
if you have a feature request or question about the library, the best place to post
those is currently on GitHub as an issue, but that is likely change if the user
community grows. If you report an issue, please give the details needed to reproduce
the problem (version of exoplanet, its dependencies, and your platform) and a small
standalone piece of code that demonstrates the problem clearly.


### Contributing code

We welcome contributions to the codebase of all scales from typo fixes to new features,
but if you would like to add a substantial feature, it would be a good idea to first
open an issue that describes your plan so that we can discuss in advance. More details
about the process for contributing can be found in the developer docs at:
https://docs.exoplanet.codes/en/stable/user/dev/.
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
This section should include a simple, standalone code snippet that demonstrates the issue.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Your setup (please complete the following information):**
 - Version of exoplanet:
 - Operating system:
 - Python version & installation method (pip, conda, etc.):

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: enhancement
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
---
title: "`exoplanet`: Gradient-based probabilistic inference for exoplanet data & other astronomical time series"
tags:
  - Python
  - astronomy
authors:
  - name: Daniel Foreman-Mackey
    orcid: 0000-0002-9328-5652
    affiliation: 1
  - name: Rodrigo Luger
    orcid: 0000-0002-0296-3826
    affiliation: "1,2"
  - name: Eric Agol
    orcid: 0000-0002-0802-9145
    affiliation: "3,2"
  - name: Thomas Barclay
    orcid: 0000-0001-7139-2724
    affiliation: 4
  - name: Luke G. Bouma
    orcid: 0000-0002-0514-5538
    affiliation: 5
  - name: Timothy D. Brandt
    orcid: 0000-0003-2630-8073
    affiliation: 6
  - name: Ian Czekala
    orcid: 0000-0002-1483-8811
    affiliation: "7,8,9,10"
  - name: Trevor J. David
    orcid: 0000-0001-6534-6246
    affiliation: "1,11"
  - name: Jiayin Dong
    orcid: 0000-0002-3610-6953
    affiliation: "7,8"
  - name: Emily A. Gilbert
    orcid: 0000-0002-0388-8004
    affiliation: 12
  - name: Tyler A. Gordon
    orcid: 0000-0001-5253-1987
    affiliation: 3
  - name: Christina Hedges
    orcid: 0000-0002-3385-8391
    affiliation: "13,14"
  - name: Daniel R. Hey
    orcid: 0000-0003-3244-5357
    affiliation: "15,16"
  - name: Brett M. Morris
    orcid: 0000-0003-2528-3409
    affiliation: 17
  - name: Adrian M. Price-Whelan
    orcid: 0000-0003-0872-7098
    affiliation: 1
  - name: Arjun B. Savel
    orcid: 0000-0002-2454-768X
    affiliation: 18
affiliations:
  - name: Center for Computational Astrophysics, Flatiron Institute, New York, NY, USA
    index: 1
  - name: Virtual Planetary Laboratory, University of Washington, Seattle, WA, USA
    index: 2
  - name: Department of Astronomy, University of Washington, University of Washington, Seattle, WA, USA
    index: 3
  - name: Center for Space Sciences and Technology, University of Maryland, Baltimore County, Baltimore, MD, USA
    index: 4
  - name: Department of Astrophysical Sciences, Princeton University, Princeton, NJ, USA
    index: 5
  - name: Department of Physics, University of California, Santa Barbara, Santa Barbara, CA, USA
    index: 6
  - name: Department of Astronomy and Astrophysics, The Pennsylvania State University, University Park, PA, USA
    index: 7
  - name: Center for Exoplanets and Habitable Worlds, The Pennsylvania State University, University Park, PA, USA
    index: 8
  - name: Center for Astrostatistics, The Pennsylvania State University, University Park, PA, USA
    index: 9
  - name: Institute for Computational and Data Sciences, The Pennsylvania State University, University Park, PA, USA
    index: 10
  - name: Department of Astrophysics, American Museum of Natural History, New York, NY, USA
    index: 11
  - name: Department of Astronomy and Astrophysics, University of Chicago, Chicago, IL, USA
    index: 12
  - name: NASA Ames Research Center, Moffett Field, CA, USA
    index: 13
  - name: Bay Area Environmental Research Institute, Moffett Field, CA, USA
    index: 14
  - name: Sydney Institute for Astronomy, School of Physics, University of Sydney, Camperdown, New South Wales, Australia
    index: 15
  - name: Stellar Astrophysics Centre, Department of Physics and Astronomy, Aarhus University, Aarhus, Denmark
    index: 16
  - name: Center for Space and Habitability, University of Bern, Bern, Switzerland
    index: 17
  - name: Department of Astronomy, University of Maryland, College Park, MD, USA
    index: 18
date: 23 April 2021
bibliography: paper.bib
---

# Summary

`exoplanet` is a toolkit for probabilistic modeling of astronomical time series
data, with a focus on observations of exoplanets, using `PyMC3` [@pymc3].
`PyMC3` is a flexible and high-performance model-building language and inference
engine that scales well to problems with a large number of parameters.
`exoplanet` extends `PyMC3`’s modeling language to support many of the custom
functions and probability distributions required when fitting exoplanet datasets
or other astronomical time series.

While it has been used for other applications, such as the study of stellar
variability [e.g., @gillen20; @medina20], the primary purpose of `exoplanet` is
the characterization of exoplanets [e.g., @gilbert20; @plavchan20] or multiple-star
systems [e.g., @czekala21] using time-series photometry, astrometry, and/or
radial velocity. In particular, the typical use case would be to use one or more
of these datasets to place constraints on the physical and orbital parameters of
the system, such as planet mass or orbital period, while simultaneously taking
into account the effects of stellar variability.

# Statement of need

Time-domain astronomy is a priority of the observational astronomical community,
with huge survey datasets currently available and more forthcoming. Within this
research domain, there is significant investment into the discovery and
characterization of exoplanets, planets orbiting stars other than our Sun. These
datasets are large (on the scale of hundreds of thousands of observations per
star from space-based observatories such as _Kepler_ and _TESS_), and the
research questions are becoming more ambitious (in terms of both the
computational cost of the physical models and the flexibility of these models).
The packages in the _exoplanet_ ecosystem are designed to enable rigorous
probabilistic inference with these large datasets and high-dimensional models by
providing a high-performance and well-tested infrastructure for integrating
these models with modern modeling frameworks such as `PyMC3`. Since its initial
release at the end of 2018, `exoplanet` has been widely used, with 64
citations of the Zenodo record [@zenodo] so far.

# The _exoplanet_ software ecosystem

Besides the primary `exoplanet` package, the _exoplanet_ ecosystem of projects
includes several other libraries. This paper describes, and is the primary
reference for, this full suite of packages. The following provides a short
description of each library within this ecosystem and discusses how they are
related.

- `exoplanet`[^exoplanet] is the primary library, and it includes implementations
  of many special functions required for exoplanet data analysis. These include
  the spherical geometry for computing orbits, some exoplanet-specific
  distributions for eccentricity [@kipping13b; @vaneylen19] and limb darkening
  [@kipping13], and exposure-time integrated limb-darkened transit light curves.
- `exoplanet-core`[^exoplanet-core] provides efficient, well-tested, and
  differentiable implementations of all of the exoplanet-specific operations
  that must be compiled for performance. These include an efficient solver for
  Kepler's equation [based on the algorithm proposed by @raposo17] and limb
  darkened transit light curves [@agol20]. Besides the implementation for
  `PyMC3`, `exoplanet-core` includes implementations in `numpy` [@numpy] and
  `jax` [@jax].
- `celerite2`[^celerite2], is an updated implementation of the _celerite_
  algorithm[^celerite] [@foremanmackey17; @foremanmackey18] for scalable
  Gaussian Process regression for time series data. Like `exoplanet-core`,
  `celerite2` includes support for `numpy`, `jax`, and `PyMC3`, as well as some
  recent generalizations of the _celerite_ algorithm [@gordon20].
- `pymc3-ext`[^pymc3-ext], includes a set of helper functions to make `PyMC3`
  more amenable to the typical astronomical data analysis workflow. For example,
  it provides a tuning schedule for `PyMC3`'s sampler [based on the method used
  by the `Stan` project and described by @carpenter17] that provides better
  performance on models with correlated parameters.
- `rebound-pymc3`[^rebound-pymc3] provides an interface between _REBOUND_
  [@rein12], _REBOUNDx_ [@tamayo20], and `PyMC3` to enable inference with full
  N-body orbit integration.

# Documentation & case studies

The main documentation page for the _exoplanet_ libraries lives at
[docs.exoplanet.codes](https://docs.exoplanet.codes) where it is hosted on
[ReadTheDocs](https://readthedocs.org). The tutorials included with the
documentation are automatically executed on every push or pull request to the
GitHub repository, with the goal of ensuring that the tutorials are always
compatible with the current version of the code. The `celerite2` project has its
own documentation page at
[celerite2.readthedocs.io](https://celerite2.readthedocs.io), with tutorials
that are similarly automatically executed.

Alongside these documentation pages, there is a parallel "Case Studies" website
at [gallery.exoplanet.codes](https://gallery.exoplanet.codes) that includes more
detailed example use cases for `exoplanet` and the other libraries described
here. Like the tutorials on the documentation page, these case studies are
automatically executed using GitHub Actions, but at lower cadence (once a week
and when a new release of the `exoplanet` library is made) since the runtime is
much longer. \autoref{fig:figure} shows the results of two example case studies
demonstrating some of the potential use cases of the `exoplanet` software
ecosystem.

![Some examples of datasets fit using `exoplanet`. The full analyses behind
these examples are available on the "Case Studies" page as Jupyter notebooks.
(left) A fit to the light curves of a transiting exoplanet observed by two
different space-based photometric surveys: Kepler and TESS. (right) The phase-folded
radial velocity time series for an exoplanet observed from different
observatories with different instruments, fit simultaneously using `exoplanet`.
\label{fig:figure}](figures/figure.png)

# Similar tools

There is a rich ecosystem of tooling available for inference with models such as
the ones supported by `exoplanet`. Each of these tools has its own set of
strengths and limitations and we will not make a detailed comparison here, but
it is worth listing some of these tools and situating `exoplanet` in this
context.

Some of the most popular tools in this space include (and note that this is far
from a comprehensive list!) `EXOFAST` [@eastman13; @eastman19], `radvel`
[@fulton18], `juliet` [@espinoza19], `exostriker` [@trifonov19], `PYANETI`
[@barragan19], `allesfitter` [@guenther20], and `orbitize` [@blunt20]. Similar
tools also exist for modeling observations of eclipsing binary systems,
including `JKTEBOP` [@southworth04], `eb` [@irwin11], and `PHOEBE` [@conroy20].
These packages all focus on providing a high-level interface for designing
models and then executing a fit. In contrast, `exoplanet` is designed to be lower
level and more conceptually similar to tools like `batman` [@kreidberg15],
`PyTransit` [@parviainen15], `ldtk` [@parviainen15b], `ellc` [@maxted16],
`starry` [@luger19], or `Limbdark.jl` [@agol20], which provide the building
blocks for evaluating the models required for inference with exoplanet datasets.
In fact, several of the higher-level packages listed above include these
lower-level libraries as dependencies, and our hope is that `exoplanet` could
provide the backend for future high-level libraries.

As emphasized in the title of this paper, the main selling point of `exoplanet`
when compared to other tools in this space is that it supports differentiation
of all components of the model and is designed to integrate seamlessly with the
`aesara` [@aesara] automatic differentiation framework used by `PyMC3`. It is
worth noting that `aesara` was previously known as `Theano` [@theano], so these
names are sometimes used interchangeably in the `PyMC3` or `exoplanet`
documentation[^theano-aesara]. This allows the use of modern inference
algorithms such as No U-Turn Sampling [@hoffman14] or Automatic Differentiation
Variational Inference [@kucukelbir17]. These algorithms can have some
computational and conceptual advantages over inference methods that do not use
gradients, especially for high-dimensional models. The computation of gradients
is also useful for model optimization; this is necessary when, say, searching
for new exoplanets, mapping out degeneracies or multiple modes of a posterior,
or estimating uncertainties from a Hessian. Care has been taken to provide
gradients which are numerically stable, and more accurate and faster to evaluate
than finite-difference gradients.

# Acknowledgements

We would like to thank the Astronomical Data Group at Flatiron for listening to
every iteration of this project and for providing great feedback every step of
the way.

This research was partially conducted during the _Exostar19_ program at the
_Kavli Institute for Theoretical Physics_ at UC Santa Barbara, which was
supported in part by the National Science Foundation under Grant No. NSF
PHY-1748958.

Besides the software cited above, `exoplanet` is also built on top of `ArviZ`
[@arviz] and `AstroPy` [@astropy13; @astropy18].

# References

[^exoplanet]: [https://github.com/exoplanet-dev/exoplanet](https://github.com/exoplanet-dev/exoplanet)
[^exoplanet-core]: [https://github.com/exoplanet-dev/exoplanet-core](https://github.com/exoplanet-dev/exoplanet-core)
[^celerite2]: [https://celerite2.readthedocs.io](https://celerite2.readthedocs.io)
[^celerite]: [https://celerite.readthedocs.io](https://celerite.readthedocs.io)
[^pymc3-ext]: [https://github.com/exoplanet-dev/pymc3-ext](https://github.com/exoplanet-dev/pymc3-ext)
[^rebound-pymc3]: [https://github.com/exoplanet-dev/rebound-pymc3](https://github.com/exoplanet-dev/rebound-pymc3)
[^theano-aesara]: More information about this distinction is available at [https://docs.exoplanet.codes/en/stable/user/theano/](https://docs.exoplanet.codes/en/stable/user/theano/)
---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.6
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Light travel time delay

```{code-cell}
:nbsphinx: hidden

import exoplanet

exoplanet.utils.docs_setup()
print(f"exoplanet.__version__ = '{exoplanet.__version__}'")
```

A simple example showing the effects of light travel time delay on an edge-on planet in an orbit similar to that of Earth.

```{code-cell}
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
import exoplanet as xo
```

Let's instantiate an orbit corresponding to a massless Earth orbiting the Sun in a circular orbit, viewed edge-on:

```{code-cell}
orbit = xo.orbits.KeplerianOrbit(
    period=365.25,
    t0=0.0,
    incl=0.5 * np.pi,
    ecc=0.0,
    omega=0.0,
    Omega=0.0,
    m_planet=0.0,
    m_star=1.0,
    r_star=1.0,
)
```

In the observer's coordinate system, the Earth will be seen to transit the Sun when its `x` position is zero. We specified the reference transit time to be `t0 = 0`, so let's check that this is indeed the case. Below we compute the time closest to the point where `x = 0`:

```{code-cell}
t = np.linspace(-0.01, 0.01, 9999)
x, y, z = orbit.get_planet_position(t)
t0_obs1 = t[np.argmin(np.abs(x.eval()))]
print("{:.4f}".format(t0_obs1))
```

No surprises here! But what happens if we add light travel time delay into the mix?

```{code-cell}
x, y, z = orbit.get_planet_position(t, light_delay=True)
t0_obs2 = t[np.argmin(np.abs(x.eval()))]
print("{:.4f}".format(t0_obs2))
```

Recall that the time is in days by default. Let's see what it is in minutes:

```{code-cell}
print("{:.2f}".format((t0_obs2 * u.day).to(u.minute)))
```

So when we account for light travel time delay, the Earth appears to transit the Sun *early* by about 8.3 minutes. In `exoplanet`, the reference time `t0` is specified relative to the origin of the coordinate system (in this case, the position of the Sun). You can check that this is in fact the [light travel time from the Sun to the Earth](https://image.gsfc.nasa.gov/poetry/venus/q89.html).

+++

Since transits happen early, you might guess that the opposite is true for secondary eclipses: they appear to occur late. Let's check this. Absent any time delay, the center of the eclipse should again occur when `x = 0`, half a period after the transit:

```{code-cell}
t = np.linspace(
    0.5 * orbit.period.eval() - 0.01, 0.5 * orbit.period.eval() + 0.01, 9999
)
x, y, z = orbit.get_planet_position(t)
t0_obs1 = t[np.argmin(np.abs(x.eval()))]
print("{:.4f}".format(t0_obs1))
```

```{code-cell}
x, y, z = orbit.get_planet_position(t, light_delay=True)
t0_obs2 = t[np.argmin(np.abs(x.eval()))]
print("{:.4f}".format(t0_obs2))
```

The difference between those two is

```{code-cell}
print("{:.2f}".format(((t0_obs2 - t0_obs1) * u.day).to(u.minute)))
```

as expected.

+++

As an example, here's a comparison between a transit light curve that includes light travel delay and one that doesnt't:

```{code-cell}
# Region around transit
t = np.linspace(-0.5, 0.5, 20000)

# No light delay
light_curve1 = (
    xo.LimbDarkLightCurve(0.5, 0.25)
    .get_light_curve(orbit=orbit, r=0.1, t=t)
    .eval()
)

# With light delay
light_curve2 = (
    xo.LimbDarkLightCurve(0.5, 0.25)
    .get_light_curve(orbit=orbit, r=0.1, t=t, light_delay=True)
    .eval()
)

# Plot!
fig = plt.figure(figsize=(8, 3))
plt.plot(t, light_curve1, label="no delay")
plt.plot(t, light_curve2, label="with delay")
plt.xlabel("time [days]")
plt.ylabel("relative flux")
plt.legend(fontsize=10, loc="lower right")
_ = plt.title(
    "Light delay causes transits to occur 8.32 minutes early", fontsize=14
)
```
---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.6
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# About these tutorials

+++

This and the following tutorials are automatically executed with every change of the code to make sure that they are always up to date with the code.
As a result, they are designed to require only a relatively small amount of computation time; when using `exoplanet` for research you will probably find that your runtimes are longer.
For more in-depth tutorials with real-world applications and real data, check out the [Case Studies page](https://gallery.exoplanet.codes).

At the top of each tutorial, you'll find a cell like the following that indicates the version of `exoplanet` that was used to generate the tutorial:

```{code-cell}
import exoplanet

exoplanet.utils.docs_setup()
print(f"exoplanet.__version__ = '{exoplanet.__version__}'")
```

That cell also includes a call to the `exoplanet.utils.docs_setup` function that will squash some warnings (these are generally caused by Theano/Aesara; see {ref}`theano` for more info) and set up our `matplotlib` style.

To exectute a tutorial on your own, you can click on the buttons at the top right or this page to launch the notebook using [Binder](https://mybinder.org) or download the `.ipynb` file directly.

```{code-cell}

```
---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.6
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

(autodiff)=

# Automatic differentation & gradient-based inference

```{code-cell}
import exoplanet

exoplanet.utils.docs_setup()
print(f"exoplanet.__version__ = '{exoplanet.__version__}'")
```

The major selling point of `exoplanet` compared to other similar libraries is that it integrates with the `PyMC3` probabilistic modeling framework.
`PyMC3` offers gradient-based inference algorithms (more on these below) that can be more computationally efficient (per effective sample) than other tools commonly used for probabilistic inference in astrophysics.
When I am writing this, gradient-based inference methodology is not widely used in astro because it can be difficult to compute the relevant gradients (derivatives of your log probability function with respect to the parameters).
`PyMC3` (and many other frameworks like it) handle this issue using [automatic differentiation (AD)](https://en.wikipedia.org/wiki/Automatic_differentiation), a method (or collection of methods) for automatically propagating derivatives through your code.

It's beyond the scope of this tutorial to go into too many details about AD and most users of `exoplanet` shouldn't need to interact with this too much, but this should at least give you a little taste of the kinds of things AD can do for you and demonstrate how this translates into efficient inference with probabilistic models.
The main thing that I want to emphasize here is that AD is not the same as *symbolic differentiation* (it's not going to provide you with a mathematical expression for your gradients), but it's also not the same as numerical methods like [finite difference](https://en.wikipedia.org/wiki/Finite_difference).
Using AD to evaluate the gradients of your model will generally be faster, more efficient, and more numerically stable than alternatives, but there are always exceptions to any rule.
There are times when providing your AD framework with a custom implementation and/or differentation rule for a particular function is beneficial in terms of cost and stability.
`exoplanet` is designed to provide these custom implementations only where it is useful (e.g. solving Kepler's equation or evaluating limb-darkened light curves) and then rely on the existing AD toolkit elsewhere.

## Automatic differentation in Theano/Aesara

Let's start by giving a quick taste for how AD works in the context of `PyMC3` and `exoplanet`.
`PyMC3` is built on top of a library called `Aesara` (formerly `Theano`, see {ref}`theano`) that was designed for AD and fast linear algebra in the context of neural networks.
The learning curve can be steep because your code ends up looking a little different from the usual way you would write Python code.
(If you're looking for something more "Pythonic", [JAX](https://github.com/google/jax) might be a good alternative, but for our purposes we'll stick with `Aesara`.)
In `Aesara`, you start be defining the *relationships* between variables and then ask to evaluate a variable using this "graph".
This might seem a little counterintuitive and you hopefully won't have to worry about this too much, but it can be a handy workflow because `Aesara` will compile and optimize your model in the background (generally making it faster).
This structure allows us to evaluate the derivative of some variable with respect to another.

Here's a simple example where we define a relationship between two variables `x` and `y`:

$$
y = \exp\left[\sin\left(\frac{2\,\pi\,x}{3}\right)\right]
$$

and then calculate the derivative using AD.
For comparison, the symbolic derivative is:

$$
\frac{\mathrm{d}y}{\mathrm{d}x} = \frac{2\,\pi}{3}\,\exp\left[\sin\left(\frac{2\,\pi\,x}{3}\right)\right]\,\cos\left(\frac{2\,\pi\,x}{3}\right)
$$

```{code-cell}
import numpy as np
import matplotlib.pyplot as plt
from aesara_theano_fallback import aesara
import aesara_theano_fallback.tensor as at

# Define the relationship between two variables x and y=f(x)
x_ = at.vector("x")
y_ = at.exp(at.sin(2 * np.pi * x_ / 3.0))

# Here aesara will compile a function to evaluate y
func = aesara.function([x_], y_)

# Request the gradient of y with respect to x
# Note the call to `sum`. This is a bit of a cheat since by
# default Aesara only computes gradients of a *scalar* function
grad = aesara.function([x_], aesara.grad(at.sum(y_), x_))

# Plot the function and its derivative
x = np.linspace(-3, 3, 500)
plt.plot(x, func(x), "k", lw=1, label="f(x)")
plt.plot(x, grad(x), "C0", label="df(x)/dx; AD", lw=3)

# Overplot the symbolic derivative for comparison
plt.plot(
    x,
    2 * np.pi * func(x) * np.cos(2 * np.pi * x / 3.0) / 3.0,
    "C1",
    lw=1,
    label="df(x)/dx; symbolic",
)

plt.legend()
plt.title("a silly function and its derivative")
plt.ylabel("f(x); df(x)/dx")
_ = plt.xlabel("x")
```

This example is obviously pretty artificial, but I think that you can imagine how something like this would start to come in handy when your models get more complicated.
In particular, I think that you'll probably pretty regularly find yourself experimenting with different choices of parameters and it would be a real pain to be required to re-write all your derivative code for every new choice of parameterization.

## Gradient-based inference

Now that we have a way of easily evaluating the derivatives of scalar functions, we have a wide array of new inference algorithms at our disposal.
For example, if we use Aesara to evaluate our log probability function, we could pass the value and gradient to [`scipy.optimize.minimize`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html) in order to find the maximum likelihood or maximum a posteriori parameters.
But, the point of `exoplanet` is that it integrates with `PyMC3` which provides implementations of various gradient-based Markov chain Monte Carlo (MCMC) methods like [Hamiltonian Monte Carlo](https://en.wikipedia.org/wiki/Hamiltonian_Monte_Carlo) and [No U-Turn Sampling](https://arxiv.org/abs/1111.4246).
`PyMC3` also supports gradient-based methods for [variational inference](https://en.wikipedia.org/wiki/Variational_Bayesian_methods), but (so far) `exoplanet` had mostly been used for MCMC.

Take a look at {ref}`intro-to-pymc3` for some specific introductory examples of using `PyMC3` for MCMC, but the general motivation for using `PyMC3` (instead of [emcee](https://emcee.readthedocs.io) or [dynesty](https://dynesty.readthedocs.io), for example) is that it can be far more computationally efficient (per effective sample), especially when your model has more than a couple of parameters.
We're not going to do a systematic study of the relative performance of different methods because that is fraught and can be problem specific, but let's look at another very artificial example where we compare the performance of sampling a 50-dimensional Gaussian using `PyMC3` and [`emcee`](https://emcee.readthedocs.io).

```{code-cell}
import time
import emcee
import logging
import pymc3 as pm
import arviz as az
import pymc3_ext as pmx

# Make PyMC3 less noisy just for this example!
logger = logging.getLogger("pymc3")
logger.setLevel(logging.ERROR)

# First define and sample the model using PyMC3
with pm.Model() as model:
    pm.Normal("x", shape=50)

    # Run a fake chain just to get everything to compile
    pm.sample(
        tune=1, draws=1, chains=1, return_inferencedata=True, progressbar=False
    )

    # Run and time the sampling
    start = time.time()
    pymc_trace = pm.sample(
        tune=1000,
        draws=1000,
        chains=2,
        cores=1,
        return_inferencedata=True,
        progressbar=False,
    )
    pymc_time = time.time() - start

    # Compute the cost per sample
    pymc_per_eff = pymc_time / np.mean(az.summary(pymc_trace)["ess_bulk"])

# Then sample the same model using emcee
func = pmx.get_theano_function_for_var(model.logpt, model=model)


def logp(theta):
    return func(theta)


init = np.random.randn(124, model.ndim)
nwalkers, ndim = init.shape
sampler = emcee.EnsembleSampler(nwalkers, ndim, logp)

start = time.time()
init = sampler.run_mcmc(init, 1000)
sampler.reset()
sampler.run_mcmc(init, 3000)
emcee_time = time.time() - start

emcee_trace = az.from_emcee(sampler=sampler)
emcee_per_eff = emcee_time / np.mean(az.summary(emcee_trace)["ess_bulk"])
```

After running these chains, we can compare the relative performance (remember that this is a very artificial example so caveats, small print, yadda, etc.).
To compare the performance of different MCMC samplers, the best metric is the computational cost per effective sample (for the parameter/function that you most care about).
In this case, we'll look at the total runtime (including burn-in) divided by an estaimte of the the mean effective sample size across all 50 dimensions (lower is better):

```{code-cell}
print(
    f"For PyMC3, the runtime per effective sample is: {pymc_per_eff * 1e3:.2f} ms"
)
print(
    f"For emcee, the runtime per effective sample is: {emcee_per_eff * 1e3:.2f} ms"
)
```

So, in this case, PyMC3 is a few order of magnitude more efficient.
As your models get more complex (in terms of computational cost, number of parameters, and the geometry of your parameter space), your mileage might vary, but my experience is that using `PyMC3` is almost always a good idea if you can!

:::{tip}
If you find that your model is not sampling well, this is almost always caused by **parameterization issues** and you'll probably get drastically better performance with some reparameterization.
See {ref}`reparameterization` for some tips and tricks.
:::

```{code-cell}

```
---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.6
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

(reparameterization)=

# Reparameterization

+++

One thing that you'll often find when using `exoplanet` (and `PyMC3` or other libraries for gradient-based inference) is that performance of your inference can be quite sensitive to the parameterization of your problem.
For example, you might see warnings from `PyMC3` telling you to consider reparameterizing your model (because of divergences) or (worse!) you might find that sampling is unreasonably slow or the number of effective samples is very small.
These issues can almost always be solved by investigating any sources of non-linear degeneracies between parameters in your model or heavy tails in the posterior density.

[The Stan User's Guide](https://mc-stan.org/docs/2_27/stan-users-guide/reparameterization-section.html) includes some general advice about reparameterization that can be a useful place to start, but I've found that the best choices tend to be pretty problem specific.
Throughout these tutorials and the [Case Studies](https://gallery.exoplanet.codes) we have tried to select sensible parameterizations (sometimes after many experiments) and comment wherever we've chosen something non-standard.

Some of the trickiest parameters are things like orbital eccentricity, inclination, and other angles.
In the case of eccentricity, you can sometimes reparameterize in terms of an observable (such as transit duration, see the "Quick fits for TESS light curves" case study).
For angles it can sometimes be better to parameterize in terms of sums or differences of pairs of angles (see the example at {ref}`data-and-models/astrometry` for a demo).

As I learn more general advice for reparameterization of `exoplanet` models, I'll try to keep this page updated, but in the meantime, feel free to start [a "discussion" on the GitHub repository](https://github.com/exoplanet-dev/exoplanet/discussions) if you discover a particularly good or tricky parameterization problem.

```{code-cell}

```
---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.6
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

(data-and-models)=

# Data & models

```{code-cell}
import exoplanet

exoplanet.utils.docs_setup()
print(f"exoplanet.__version__ = '{exoplanet.__version__}'")
```

This tutorial provides an overview of the kinds of datasets that `exoplanet` has been designed to model.
It's worth noting that this is certainly not a comprehensive list and `exoplanet` is designed to be *low-level* library that can be built on top of to support new use cases.

:::{tip}
For real world use cases and examples, check out the [Case Studies](https://gallery.exoplanet.codes) page.
:::

+++

## Overview

One primary interface to `exoplanet` is via the {class}`exoplanet.orbits.KeplerianOrbit` object (and its generalizations like {class}`exoplanet.orbits.TTVOrbit`).
The `KeplerianOrbit` object defines a set of independent (but possibly massive) bodies orbiting a central body on Keplerian trajectories.
These orbits are parameterized by orbital elements, and one of the primary uses of `exoplanet` is the measurement of these orbital elements given observations.

Given a `KeplerianOrbit`, users of `exoplanet` can evaluate things like the positions and velocities of all the bodies as a function of time.
For example, here's how you could define an orbit with a single body and plot the orbits:

```{code-cell}
import numpy as np
import matplotlib.pyplot as plt
import exoplanet as xo

# Define the orbit
orbit = xo.orbits.KeplerianOrbit(
    period=10.0,  # All times are in days
    t0=0.5,
    incl=0.5 * np.pi,  # All angles are in radians
    ecc=0.3,
    omega=-2.5,
    Omega=1.2,
    m_planet=0.05,  # All masses and distances are in Solar units
)

# Get the position and velocity as a function of time
t = np.linspace(0, 20, 5000)
x, y, z = orbit.get_relative_position(t)
vx, vy, vz = orbit.get_star_velocity(t)

# Plot the coordinates
# Note the use of `.eval()` throughout since the coordinates are all
# Aesara/Theano objects
fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
ax = axes[0]
ax.plot(t, x.eval(), label="x")
ax.plot(t, y.eval(), label="y")
ax.plot(t, z.eval(), label="z")
ax.set_ylabel("position of orbiting body [$R_*$]")
ax.legend(fontsize=10, loc=1)

ax = axes[1]
ax.plot(t, vx.eval(), label="$v_x$")
ax.plot(t, vy.eval(), label="$v_y$")
ax.plot(t, vz.eval(), label="$v_z$")
ax.set_xlim(t.min(), t.max())
ax.set_xlabel("time [days]")
ax.set_ylabel("velocity of central [$R_*$/day]")
_ = ax.legend(fontsize=10, loc=1)
```

The key feture of `exoplanet` is that all of the parameters to a `KeplerianOrbit` can be `PyMC3` variables.
This means that these elements are now something that you can *infer*.
For example, if we want to fit for the orbital period, we can define a `PyMC3` model like the following:

```{code-cell}
import pymc3 as pm

with pm.Model():
    log_period = pm.Normal("log_period", mu=np.log(10), sigma=2.0)
    orbit = xo.orbits.KeplerianOrbit(
        period=pm.math.exp(log_period),  # ...
    )

    # Define the rest of you model using `orbit`...
```

In the following sections, we will go through some of the specific ways that you might use this `orbit` to define a model, but first some orbital conventions!

+++

## Orbital conventions

These orbits are specified with respect to a single central body (generally the most massive body) and then a system of non-interacting bodies orbiting the central.
This is a good parameterization for exoplanetary systems and binary stars, but it is sometimes not sufficient for systems with multiple massive bodies where the interactions will be important to the dynamics.

We follow a set of internally consistent orbital conventions that also attempts to respect as many of the established conventions of the stellar and exoplanetary fields as possible.
The location of a body on a Keplerian orbit with respect to the center of mass is given in the perifocal plane by

$$
\boldsymbol{r} =
  \left [
    \begin{array}{c}
      x\\
      y\\
      z\\
    \end{array}
  \right ].
$$

By construction, the orbit of the body in the perifocal frame is constrained entirely to the $x -y$ plane as shown by this figure:

![Orbital geometry](../_static/orbit3D.png)

The range of orbital convention choices mainly stems from how the location of the body in the perifocal frame is converted to the observer frame,

$$
\boldsymbol{R} =
  \left [
    \begin{array}{c}
      X\\
      Y\\
      Z\\
    \end{array}
  \right ].
$$

We choose to orient the $\hat{\boldsymbol{X}}$ unit vector in the north direction, $\hat{\boldsymbol{Y}}$ in the east direction, and $\hat{\boldsymbol{Z}}$ towards the observer. Under this convention, the inclination of an orbit is defined by the dot-product of the orbital angular momentum vector with $\hat{\boldsymbol{Z}}$, conforming to the common astrometric convention that a zero inclination ($i = 0$) orbit is face-on, with the test body orbiting in a counter-clockwise manner around the primary body.

In the stellar and exoplanetary fields, there is less agreement on which  segment of the line of nodes is defined as the *ascending* node.
We choose to define the ascending node as the point where the orbiting body crosses the plane of the sky moving *away* from the observer (i.e., crossing from $Z > 0$ to $Z< 0$).
This convention has historically been the choice of the visual binary field; the opposite convention occasionally appears in exoplanet and planetary studies.

To implement the transformation from perifocal frame to observer frame, we consider the rotation matrices

$$
  \boldsymbol{P}_x(\phi) = \left [
  \begin{array}{ccc}
    1 & 0 & 0 \\
    0 & \cos \phi & - \sin \phi \\
    0 & \sin \phi & \cos \phi \\
    \end{array}\right]
$$

$$
  \boldsymbol{P}_z (\phi) = \left [
  \begin{array}{ccc}
    \cos \phi & - \sin \phi & 0\\
    \sin \phi & \cos \phi & 0 \\
    0 & 0 & 1 \\
    \end{array}\right].
$$

which result in a *clockwise* rotation of the axes, as defined using the right hand rule.

This means when we look down the $z$-axis, for a positive angle $\phi$, it would be as if the $x$ and $y$ axes rotated clockwise.
In order to find out what defines counter-clockwise when considering the other rotations, we look to the right hand rule and cross products of the axes unit vectors.
Since $\hat{\boldsymbol{x}} \times \hat{\boldsymbol{y}} = \hat{\boldsymbol{z}}$, when looking down the $z$ axis the direction of the $x$-axis towards the $y$-axis defines counter clockwise.
Similarly, we have $\hat{\boldsymbol{y}} \times \hat{\boldsymbol{z}} = \hat{\boldsymbol{x}}$, and $\hat{\boldsymbol{z}} \times \hat{\boldsymbol{x}} = \hat{\boldsymbol{y}}$.

To convert a position $\boldsymbol{r}$ in the perifocal frame to the observer frame $\boldsymbol{R}$ (referencing the figure above), three rotations are required:
(i) a rotation about the $z$-axis through an angle $\omega$ so that the $x$-axis coincides with the line of nodes at the ascending node,
(ii) a rotation about the $x$-axis through an angle ($-i$) so that the two planes are coincident, and finally
(iii) a rotation about the $z$-axis through an angle $\Omega$. Applying these rotation matrices yields

$$
  \boldsymbol{R} =
  \boldsymbol{P}_z (\Omega) \boldsymbol{P}_x(-i) \boldsymbol{P}_z(\omega)
  \boldsymbol{r}.
$$

As a function of true anomaly $f$, the position in the observer frame is given by

$$
  \begin{array}{lc}
    X =& r [ \cos \Omega (\cos \omega \cos f - \sin \omega \sin f)  - \sin \Omega  \cos i (\sin \omega \cos f + \cos \omega \sin f) ] \\
    Y =& r [ \sin \Omega (\cos \omega \cos f - \sin \omega \sin f) + \cos \Omega \cos i(\sin \omega \cos f + \cos \omega \sin f) ] \\
    Z =& - r \sin i (\sin \omega \cos f + \cos \omega \sin f).\\
\end{array}
$$

Simplifying the equations using the sum angle identities, we find

$$
  \begin{array}{lc}
    X =& r (\cos \Omega \cos(\omega + f) - \sin(\Omega) \sin(\omega + f) \cos(i)) \\
    Y =& r (\sin \Omega \cos(\omega + f) + \cos(\Omega) \sin(\omega + f) \cos(i)) \\
    Z =& - r \sin(\omega + f) \sin(i).\\
\end{array}
$$

The central body is defined by any two of radius $R_\star$, mass $M_\star$, or density $\rho_\star$.
The third physical parameter can always be computed using the other two so `exoplanet` will throw an error if three are given (even if they are numerically self-consistent).

The orbits are then typically defined by (check out the docstring {class}`exoplanet.orbits.KeplerianOrbit` for all the options):

* the period $P$ or semi-major axis $a$ of the orbit,
* the time of conjunction $t_0$ or time of periastron passage $t_\mathrm{p}$,
* the planet's radius $R_\mathrm{pl}$ and mass $M_\mathrm{pl}$,
* the eccentricity $e$, argument of periastron $\omega$ and the longitude of ascending node $\Omega$, and
* the inclination of the orbital plane $i$ or the impact parameter $b$.

A `KeplerianOrbit` object in `exoplanet` will always be fully specified.
For example, if the orbit is defined using the period $P$, the semi-major axis $a$ will be computed using Kepler's third law.

+++

(data-and-models/rvs)=

## Radial velocities

A {class}`exoplanet.orbits.KeplerianOrbit` can be used to compute the expected radial velocity time series for a given set of parameters.
One typical parameterization for a radial velocity fit would look something like this:

```{code-cell}
import arviz as az
import pymc3_ext as pmx
import aesara_theano_fallback.tensor as tt

# Create a dummy dataset
random = np.random.default_rng(1234)
t_plot = np.linspace(0, 50, 500)
t = np.sort(random.uniform(0, 50, 20))
rv_err = 0.5 + np.zeros_like(t)
rv_obs = 5.0 * np.sin(2 * np.pi * t / 10.0) + np.sqrt(
    0.05 ** 2 + rv_err ** 2
) * random.normal(size=len(t))

with pm.Model():

    # Period, semi-amplitude, and eccentricity
    log_period = pm.Normal("log_period", mu=np.log(10.0), sigma=1.0)
    period = pm.Deterministic("period", tt.exp(log_period))
    log_semiamp = pm.Normal("log_semiamp", mu=np.log(5.0), sd=2.0)
    semiamp = pm.Deterministic("semiamp", tt.exp(log_semiamp))
    ecc = pm.Uniform("ecc", lower=0, upper=1)

    # At low eccentricity, omega and the phase of periastron (phi) are
    # correlated so it can be best to fit in (omega ± phi) / 2
    plus = pmx.Angle("plus")
    minus = pmx.Angle("minus")
    phi = pm.Deterministic("phi", plus + minus)
    omega = pm.Deterministic("omega", plus - minus)

    # For non-zero eccentricity, it can sometimes be better to use
    # sqrt(e)*sin(omega) and sqrt(e)*cos(omega) as your parameters:
    #
    #     ecs = pmx.UnitDisk("ecs", testval=0.01 * np.ones(2))
    #     ecc = pm.Deterministic("ecc", tt.sum(ecs ** 2, axis=0))
    #     omega = pm.Deterministic("omega", tt.arctan2(ecs[1], ecs[0]))
    #     phi = pmx.Angle("phi")

    # Jitter & the system mean velocity offset
    log_jitter = pm.Normal("log_jitter", mu=np.log(0.05), sd=5.0)
    zero_point = pm.Normal("zero_point", mu=0, sd=10.0)

    # Then we define the orbit
    tperi = pm.Deterministic("tperi", period * phi / (2 * np.pi))
    orbit = xo.orbits.KeplerianOrbit(
        period=period, t_periastron=tperi, ecc=ecc, omega=omega
    )

    # And define the RV model
    rv_model = zero_point + orbit.get_radial_velocity(t, K=semiamp)

    # Finally add in the observation model
    err = tt.sqrt(rv_err ** 2 + tt.exp(2 * log_jitter))
    pm.Normal("obs", mu=rv_model, sigma=rv_err, observed=rv_obs)

    # We'll also track the model just for plotting purposes
    pm.Deterministic(
        "rv_plot", zero_point + orbit.get_radial_velocity(t_plot, K=semiamp)
    )

    soln = pmx.optimize(vars=[plus, minus, ecc])
    soln = pmx.optimize(soln)
    trace = pmx.sample(
        tune=1000,
        draws=1000,
        cores=2,
        chains=2,
        start=soln,
        return_inferencedata=True,
    )

# Plot the results
rv_plot = trace.posterior["rv_plot"].values
q16, q50, q84 = np.percentile(rv_plot, [16, 50, 84], axis=(0, 1))
plt.errorbar(t, rv_obs, yerr=rv_err, fmt="+k", label="data")
plt.plot(t_plot, q50)
plt.fill_between(t_plot, q16, q84, alpha=0.3, label="posterior")
plt.xlim(0, 50)
plt.legend(fontsize=12, loc=2)
plt.xlabel("time [days]")
plt.ylabel("radial velocity [m/s]")

az.summary(trace, var_names=["^(?!rv_plot).*"], filter_vars="regex")
```

Instead of using `semiamp` as a parameter, we could have passed `m_planet` and `incl` as parameters to `KeplerianOrbit` and fit the physical orbit.
Similarly, if we we doing a joint fit with a transit light curve (see {ref}`data-and-models/combining` below) we could have parameterized in terms of `t0` (the time of a reference transit) instead of the phase of periastron.

As mentioned in the comments above and discussed in {ref}`reparameterization`, the performance of these fits can be problem specific and depend on the choice of parameters.
Therefore it is often worthwhile experimenting with different parameterizations, but this can be a good place to start for radial velocity-only fits.

+++

(data-and-models/astrometry)=

## Astrometry

Astrometric observations usually consist of measurements of the separation and position angle of the secondary star (or directly imaged exoplanet), relative to the primary star as a function of time, but `exoplanet` could also be used to model the motion of the center of light for and unresolved orbit.
The typical {class}`exoplanet.orbits.KeplerianOrbit` definition for and astrometric dataset will be similar to a radial velocity fit:

```{code-cell}
random = np.random.default_rng(5678)
t_plot = np.linspace(0, 22 * 365.25, 500)
t = np.sort(random.uniform(0.0, 22 * 365.25, 45))
rho_err = random.uniform(0.05, 0.1, len(t))
theta_err = random.uniform(0.05, 0.1, len(t))

with pm.Model():

    # Period, semi-major axis, eccentricity, and t0
    log_period = pm.Normal("log_period", mu=np.log(25.0 * 365.25), sigma=1.0)
    period = pm.Deterministic("period", tt.exp(log_period))
    log_a = pm.Normal("log_a", mu=np.log(0.3), sd=2.0)
    a = pm.Deterministic("a", tt.exp(log_a))
    ecc = pm.Uniform("ecc", lower=0, upper=1)
    tperi = pm.Normal("tperi", mu=3500.0, sigma=1000.0)

    # For astrometric orbits, a good choice of parameterization can be
    # (Omega ± omega) / 2
    plus = pmx.Angle("plus")
    minus = pmx.Angle("minus")
    Omega = pm.Deterministic("Omega", plus + minus)
    omega = pm.Deterministic("omega", plus - minus)

    # We'll use a uniform prior on cos(incl)
    cos_incl = pm.Uniform("cos_incl", lower=-1.0, upper=1.0, testval=0.3)
    incl = pm.Deterministic("incl", tt.arccos(cos_incl))

    # Then we define the orbit
    orbit = xo.orbits.KeplerianOrbit(
        a=a,
        t_periastron=tperi,
        period=period,
        incl=incl,
        ecc=ecc,
        omega=omega,
        Omega=Omega,
    )

    # And define the RV model
    rho_model, theta_model = orbit.get_relative_angles(t)

    # ================================================== #
    # Simulate data from the model for testing           #
    # You should remove the following lines in your code #
    # ================================================== #
    rho_obs, theta_obs = pmx.eval_in_model([rho_model, theta_model])
    rho_obs += rho_err * random.normal(size=len(t))
    theta_obs += theta_err * random.normal(size=len(t))
    # =============== end simulated data =============== #

    # Define the observation model this is simple for rho:
    pm.Normal("rho_obs", mu=rho_model, sd=rho_err, observed=rho_obs)

    # But we want to be cognizant of the fact that theta wraps so the following
    # is equivalent to
    #
    #   pm.Normal("obs_theta", mu=theta_model, observed=theta_obs, sd=theta_err)
    #
    # but takes into account the wrapping
    theta_diff = tt.arctan2(
        tt.sin(theta_model - theta_obs), tt.cos(theta_model - theta_obs)
    )
    pm.Normal("theta_obs", mu=theta_diff, sd=theta_err, observed=0.0)

    # We'll also track the model just for plotting purposes
    rho_plot, theta_plot = orbit.get_relative_angles(t_plot)
    pm.Deterministic("rho_plot", rho_plot)
    pm.Deterministic("theta_plot", theta_plot)

    trace = pmx.sample(
        tune=1000,
        draws=1000,
        cores=2,
        chains=2,
        target_accept=0.95,
        return_inferencedata=True,
    )

# Plot the results
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

q16, q50, q84 = np.percentile(
    trace.posterior["rho_plot"].values, [16, 50, 84], axis=(0, 1)
)
ax1.errorbar(t, rho_obs, yerr=rho_err, fmt="+k", label="data")
ax1.plot(t_plot, q50)
ax1.fill_between(t_plot, q16, q84, alpha=0.3, label="posterior")
ax1.set_ylabel(r"$\rho$ [arcsec]")

q16, q50, q84 = np.percentile(
    trace.posterior["theta_plot"].values, [16, 50, 84], axis=(0, 1)
)
ax2.errorbar(t, theta_obs, yerr=rho_err, fmt="+k", label="data")
ax2.plot(t_plot, q50)
ax2.fill_between(t_plot, q16, q84, alpha=0.3, label="posterior")
ax2.set_xlim(t_plot.min(), t_plot.max())
ax2.legend(fontsize=12, loc=2)
ax2.set_ylabel(r"$\theta$ [radians]")
ax2.set_xlabel("time [days]")

# Compute the convergence stats
az.summary(trace, var_names=["^(?!.*_plot).*"], filter_vars="regex")
```

In this example, the units of `a` are arcseconds, but the {func}`exoplanet.orbits.KeplerianOrbit.get_relative_angles` function accepts a parallax argument if you have a constraint on the parallax of the system.
If you provide parallax as an argument, `a` should be provided in the usual units of Solar radii.

+++

(data-and-models/transits)=

## Transits, occultations, and eclipses

`exoplanet` has built in support for evaluating quadratically limb darkened light curve models using the algorithm from [Agol et al. (2020)](https://arxiv.org/abs/1908.03222).
If you need flexible surface models or higher order limb darkening, check out the [`starry` package](https://starry.readthedocs.io) which also integrates with `PyMC3`.

Transit and occultation modeling is one of the primary applications of `exoplanet` so there are quite a few options (including transit timing variations, detached eclipsing binary modeling, and much more) that are highlighted on the [Case Studies](https://gallery.exoplanet.codes) page.
But, a bread-and-butter transit model implemented in `exoplanet` might look something like the following:

```{code-cell}
random = np.random.default_rng(123)
num_transits = 4
t = np.arange(0, 35, 0.02)
yerr = 5e-4

with pm.Model():

    # The baseline flux
    mean = pm.Normal("mean", mu=0.0, sigma=1.0)

    # Often the best parameterization is actually in terms of the
    # times of two reference transits rather than t0 and period
    t0 = pm.Normal("t0", mu=4.35, sigma=1.0)
    t1 = pm.Normal("t1", mu=33.2, sigma=1.0)
    period = pm.Deterministic("period", (t1 - t0) / num_transits)

    # The Kipping (2013) parameterization for quadratic limb darkening
    # paramters
    u = xo.distributions.QuadLimbDark("u", testval=np.array([0.3, 0.2]))

    # The radius ratio and impact parameter; these parameters can
    # introduce pretty serious covariances and are ripe for
    # reparameterization
    log_r = pm.Normal("log_r", mu=np.log(0.04), sigma=2.0)
    r = pm.Deterministic("r", tt.exp(log_r))
    b = xo.distributions.ImpactParameter("b", ror=r, testval=0.35)

    # Set up a Keplerian orbit for the planets
    orbit = xo.orbits.KeplerianOrbit(period=period, t0=t0, b=b)

    # Compute the model light curve; note that we index using `[:, 0]`
    # since `get_light_curve` returns an object with the shape
    # `(n_times, n_planets)`
    light_curve = (
        xo.LimbDarkLightCurve(u[0], u[1]).get_light_curve(
            orbit=orbit, r=r, t=t
        )[:, 0]
        + mean
    )

    # Here we track the value of the model light curve for plotting
    # purposes
    pm.Deterministic("light_curve", light_curve)

    # ================================================== #
    # Simulate data from the model for testing           #
    # You should remove the following lines in your code #
    # ================================================== #
    y = pmx.eval_in_model(light_curve)
    y += yerr * random.normal(size=len(y))
    # =============== end simulated data =============== #

    # The likelihood function assuming known Gaussian uncertainty
    pm.Normal("obs", mu=light_curve, sd=yerr, observed=y)

    trace = pmx.sample(
        tune=1000,
        draws=1000,
        cores=2,
        chains=2,
        return_inferencedata=True,
    )

# Plot the results
q16, q50, q84 = np.percentile(
    trace.posterior["light_curve"].values, [16, 50, 84], axis=(0, 1)
)
plt.plot(t, y, ".k", ms=2, label="data")
plt.plot(t, q50)
plt.fill_between(t, q16, q84, alpha=0.3, label="posterior")
plt.xlim(0.0, 35)
plt.legend(fontsize=12, loc=3)
plt.xlabel("time [days]")
plt.ylabel("relative flux")

# Compute the convergence stats
az.summary(trace, var_names=["^(?!light_curve).*"], filter_vars="regex")
```

(data-and-models/combining)=

## Combining datasets

Since `exoplanet` is built on top of `PyMC3`, it has the capacity to support essentially arbitrariliy complicated models.
This means that you can share parameters or fit multiple datasets however you want.
We won't go into too many details about this here, but you can see some examples on the [Case Studies](https://gallery.exoplanet.codes) of joint transit/radial velocity fits, or inferences based on datasets from multiple instruments.

The basic idea behind how to simultaneously fit multiple datasets is that you'll need to define the relationship between the datasets (perhaps share parameters or a shared `KeplerianOrbit`) and then you'll include the likelihood for each dataset in your model.
In all of the above examples, we had a line like

```python
pm.Normal("obs", mu=rv_model, sigma=rv_err, observed=rv_obs)
```

in all of our models.
This defines a Gaussian likelihood conditioned on the `observed` data.
To combine datasets, you can simply add multiple lines like this (one for each dataset) and, behind the scenes, `PyMC3` will multiply these likelihoods (or actually add their logarithms) as it should.

For more concrete examples, check out the [Case Studies](https://gallery.exoplanet.codes) and (if that's not sufficient) feel free to start [a "discussion" on the GitHub repository](https://github.com/exoplanet-dev/exoplanet/discussions) asking for help.

```{code-cell}

```
---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.6
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Citing exoplanet & its dependencies

+++

The _exoplanet_ package is mostly just glue that connects many other ideas and software.
In a situation like this, it can be easy to forget about the important infrastructure upon which our science is built.
In order to make sure that you can easily give credit where credit is due, we have tried to make it as painless as possible to work out which citations are expected for a model fit using _exoplanet_ by including a {func}`exoplanet.citations.get_citations_for_model` function that introspects the current PyMC3 model and constructs a list of citations for the functions used in that model.

For example, you might compute a quadratically limb darkened light curve using `starry` (via the {class}`exoplanet.LimbDarkLightCurve` class):

```{code-cell}
import pymc3 as pm
import exoplanet as xo

with pm.Model() as model:
    u = xo.distributions.QuadLimbDark("u")
    orbit = xo.orbits.KeplerianOrbit(period=10.0)
    light_curve = xo.LimbDarkLightCurve(u[0], u[1])
    transit = light_curve.get_light_curve(r=0.1, orbit=orbit, t=[0.0, 0.1])

    txt, bib = xo.citations.get_citations_for_model()
```

The {func}`exoplanet.citations.get_citations_for_model` function would generate an acknowledgement that cites:

- [PyMC3](https://docs.pymc.io/#citing-pymc3): for the inference engine and modeling framework,
- [Theano/Aesara](https://aesara.readthedocs.io/en/latest/citation.html): for the numerical infrastructure,
- [AstroPy](http://www.astropy.org/acknowledging.html): for units and constants,
- [Kipping (2013)](https://arxiv.org/abs/1308.0009): for the reparameterization of the limb darkening parameters for a quadratic law, and
- [Luger, et al. (2018)](https://arxiv.org/abs/1810.06559): for the light curve calculation.

The first output from {func}`exoplanet.citations.get_citations_for_model` gives the acknowledgement text:

```{code-cell}
print(txt)
```

And the second output is a string with BibTeX entries for each of the citations in the acknowledgement text:

```{code-cell}
print(bib.split("\n\n")[0] + "\n\n...")
```

```{code-cell}

```
---
jupytext:
  encoding: '# -*- coding: utf-8 -*-'
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.6
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

(intro-to-pymc3)=

# A quick intro to PyMC3

```{code-cell}
import exoplanet

exoplanet.utils.docs_setup()
print(f"exoplanet.__version__ = '{exoplanet.__version__}'")
```

Gradient-based inference methods (like [Hamiltonian Monte Carlo (HMC)](https://en.wikipedia.org/wiki/Hamiltonian_Monte_Carlo)) haven't been widely used in astrophysics, but they are the standard methods for probabilistic inference using Markov chain Monte Carlo (MCMC) in many other fields.
*exoplanet* is designed to provide the building blocks for fitting many exoplanet datasets using this technology, and this tutorial presents some of the basic features of the [PyMC3](https://docs.pymc.io/) modeling language and inference engine.
The [documentation for PyMC3](https://docs.pymc.io/) includes many other tutorials that you should check out to get more familiar with the features that are available.

In this tutorial, we will go through two simple examples of fitting some data using PyMC3.
The first is the classic fitting a line to data with unknown error bars, and the second is a more relevant example where we fit a radial velocity model to the public radial velocity observations of [51 Peg](https://en.wikipedia.org/wiki/51_Pegasi).
You can read more about fitting lines to data [in the bible of line fitting](https://arxiv.org/abs/1008.4686) and you can see another example of fitting the 51 Peg data using HMC (this time using [Stan](http://mc-stan.org)) [here](https://dfm.io/posts/stan-c++/).

## Hello world (AKA fitting a line to data)

My standard intro to a new modeling language or inference framework is to fit a line to data.
So. Let's do that with PyMC3.

To start, we'll generate some fake data using a linear model.
Feel free to change the random number seed to try out a different dataset.

```{code-cell}
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(42)

true_m = 0.5
true_b = -1.3
true_logs = np.log(0.3)

x = np.sort(np.random.uniform(0, 5, 50))
y = true_b + true_m * x + np.exp(true_logs) * np.random.randn(len(x))

plt.plot(x, y, ".k")
plt.ylim(-2, 2)
plt.xlabel("x")
_ = plt.ylabel("y")
```

To fit a model to these data, our model will have 3 parameters: the slope $m$, the intercept $b$, and the log of the uncertainty $\log(\sigma)$.
To start, let's choose broad uniform priors on these parameters:

$$
\begin{eqnarray}
p(m) &=& \left\{\begin{array}{ll}
1/10 & \mathrm{if}\,-5 < m < 5 \\
0 & \mathrm{otherwise} \\
\end{array}\right. \\
p(b) &=& \left\{\begin{array}{ll}
1/10 & \mathrm{if}\,-5 < b < 5 \\
0 & \mathrm{otherwise} \\
\end{array}\right. \\
p(\log(\sigma)) &=& \left\{\begin{array}{ll}
1/10 & \mathrm{if}\,-5 < \log(\sigma) < 5 \\
0 & \mathrm{otherwise} \\
\end{array}\right.
\end{eqnarray}
$$

Then, the log-likelihood function will be

$$
\log p(\{y_n\}\,|\,m,\,b,\,\log(\sigma)) = -\frac{1}{2}\sum_{n=1}^N \left[\frac{(y_n - m\,x_n - b)^2}{\sigma^2} + \log(2\,\pi\,\sigma^2)\right]
$$

[**Note:** the second normalization term is needed in this model because we are fitting for $\sigma$ and the second term is *not* a constant.]

Another way of writing this model that might not be familiar is the following:

$$
\begin{eqnarray}
m &\sim& \mathrm{Uniform}(-5,\,5) \\
b &\sim& \mathrm{Uniform}(-5,\,5) \\
\log(\sigma) &\sim& \mathrm{Uniform}(-5,\,5) \\
y_n &\sim& \mathrm{Normal}(m\,x_n+b,\,\sigma)
\end{eqnarray}
$$

This is the way that a model like this is often defined in statistics and it will be useful when we implement out model in PyMC3 so take a moment to make sure that you understand the notation.

Now, let's implement this model in PyMC3.
The documentation for the distributions available in PyMC3's modeling language can be [found here](https://docs.pymc.io/api/distributions/continuous.html) and these will come in handy as you go on to write your own models.

```{code-cell}
import pymc3 as pm

with pm.Model() as model:

    # Define the priors on each parameter:
    m = pm.Uniform("m", lower=-5, upper=5)
    b = pm.Uniform("b", lower=-5, upper=5)
    logs = pm.Uniform("logs", lower=-5, upper=5)

    # Define the likelihood. A few comments:
    #  1. For mathematical operations like "exp", you can't use
    #     numpy. Instead, use the mathematical operations defined
    #     in "pm.math".
    #  2. To condition on data, you use the "observed" keyword
    #     argument to any distribution. In this case, we want to
    #     use the "Normal" distribution (look up the docs for
    #     this).
    pm.Normal("obs", mu=m * x + b, sd=pm.math.exp(logs), observed=y)

    # This is how you will sample the model. Take a look at the
    # docs to see that other parameters that are available.
    trace = pm.sample(
        draws=1000, tune=1000, chains=2, cores=2, return_inferencedata=True
    )
```

Now since we now have samples, let's make some diagnostic plots.
The first plot to look at is the "traceplot" implemented in PyMC3.
In this plot, you'll see the marginalized distribution for each parameter on the left and the trace plot (parameter value as a function of step number) on the right.
In each panel, you should see two lines with different colors.
These are the results of different independent chains and if the results are substantially different in the different chains then there is probably something going wrong.

```{code-cell}
import arviz as az

_ = az.plot_trace(trace, var_names=["m", "b", "logs"])
```

It's also good to quantify that "looking substantially different" argument.
This is implemented in PyMC3 as the "summary" function.
In this table, some of the key columns to look at are `n_eff` and `Rhat`.

1. `n_eff` shows an estimate of the number of effective (or independent) samples for that parameter. In this case, `n_eff` should probably be around 500 per chain (there should have been 2 chains run).
2. `Rhat` shows the [Gelman–Rubin statistic](https://docs.pymc.io/api/diagnostics.html#pymc3.diagnostics.gelman_rubin) and it should be close to 1.

```{code-cell}
az.summary(trace, var_names=["m", "b", "logs"])
```

The last diagnostic plot that we'll make here is the [corner plot made using corner.py](https://corner.readthedocs.io).

```{code-cell}
import corner

_ = corner.corner(
    trace,
    truths=dict(m=true_m, b=true_b, logs=true_logs),
)
```

**Extra credit:** Here are a few suggestions for things to try out while getting more familiar with PyMC3:

1. Try initializing the parameters using the `testval` argument to the distributions. Does this improve performance in this case? It will substantially improve performance in more complicated examples.
2. Try changing the priors on the parameters. For example, try the "uninformative" prior [recommended by Jake VanderPlas on his blog](http://jakevdp.github.io/blog/2014/06/14/frequentism-and-bayesianism-4-bayesian-in-python/#Prior-on-Slope-and-Intercept).
3. What happens as you substantially increase or decrease the simulated noise? Does the performance change significantly? Why?

## A more realistic example: radial velocity exoplanets

While the above example was cute, it doesn't really fully exploit the power of PyMC3 and it doesn't really show some of the real issues that you will face when you use PyMC3 as an astronomer.
To get a better sense of how you might use PyMC3 in Real Life™, let's take a look at a more realistic example: fitting a Keplerian orbit to radial velocity observations.

One of the key aspects of this problem that I want to highlight is the fact that PyMC3 (and the underlying model building framework [Theano (recently renamed to Aesara)](https://aesara.readthedocs.io/)) don't have out-of-the-box support for the root-finding that is required to solve Kepler's equation.
As part of the process of computing a Keplerian RV model, we must solve the equation:

$$
M = E - e\,\sin E
$$

for the eccentric anomaly $E$ given some mean anomaly $M$ and eccentricity $e$.
There are commonly accepted methods of solving this equation using [Newton's method](https://en.wikipedia.org/wiki/Newton%27s_method), but if we want to expose that to PyMC3, we have to define a [custom Theano operation](https://aesara.readthedocs.io/en/latest/extending/index.html) with a custom gradient.
I won't go into the details of the math (because [I blogged about it](https://dfm.io/posts/stan-c++/)) and I won't go into the details of the implementation.
So, for this tutorial, we'll use the custom Kepler solver that is implemented as part of *exoplanet* and fit the publicly available radial velocity observations of the famous exoplanetary system 51 Peg using PyMC3.

First, we need to download the data from the exoplanet archive:

```{code-cell}
import requests
import pandas as pd
import matplotlib.pyplot as plt

# Download the dataset from the Exoplanet Archive:
url = "https://exoplanetarchive.ipac.caltech.edu/data/ExoData/0113/0113357/data/UID_0113357_RVC_001.tbl"
r = requests.get(url)
if r.status_code != requests.codes.ok:
    r.raise_for_status()
data = np.array(
    [
        l.split()
        for l in r.text.splitlines()
        if not l.startswith("\\") and not l.startswith("|")
    ],
    dtype=float,
)
t, rv, rv_err = data.T
t -= np.mean(t)

# Plot the observations "folded" on the published period:
# Butler et al. (2006) https://arxiv.org/abs/astro-ph/0607493
lit_period = 4.230785
plt.errorbar(
    (t % lit_period) / lit_period, rv, yerr=rv_err, fmt=".k", capsize=0
)
plt.xlim(0, 1)
plt.ylim(-110, 110)
plt.annotate(
    "period = {0:.6f} days".format(lit_period),
    xy=(1, 0),
    xycoords="axes fraction",
    xytext=(-5, 5),
    textcoords="offset points",
    ha="right",
    va="bottom",
    fontsize=12,
)
plt.ylabel("radial velocity [m/s]")
_ = plt.xlabel("phase")
```

Now, here's the implementation of a radial velocity model in PyMC3.
Some of this will look familiar after the Hello World example, but things are a bit more complicated now.
Take a minute to take a look through this and see if you can follow it.
There's a lot going on, so I want to point out a few things to pay attention to:

1. All of the mathematical operations (for example `exp` and `sqrt`) are being performed using Theano/Aesara (see {ref}`theano`) instead of NumPy.
2. All of the parameters have initial guesses provided. This is an example where this makes a big difference because some of the parameters (like period) are very tightly constrained.
3. Some of the lines are wrapped in `Deterministic` distributions. This can be useful because it allows us to track values as the chain progresses even if they're not parameters. For example, after sampling, we will have a sample for `bkg` (the background RV trend) for each step in the chain. This can be especially useful for making plots of the results.
4. Similarly, at the end of the model definition, we compute the RV curve for a single orbit on a fine grid. This can be very useful for diagnosing fits gone wrong.
5. For parameters that specify angles (like $\omega$, called `w` in the model below), it can be inefficient to sample in the angle directly because of the fact that the value wraps around at $2\pi$. Instead, it can be better to sample the unit vector specified by the angle or as a parameter in a unit disk, when combined with eccentricity. In practice, this can be achieved by sampling a 2-vector from an isotropic Gaussian and normalizing the components by the norm. These are implemented in the [pymc3-ext](https://github.com/exoplanet-dev/pymc3-ext) package.

```{code-cell}
import pymc3_ext as pmx
import aesara_theano_fallback.tensor as tt

import exoplanet as xo

with pm.Model() as model:

    # Parameters
    logK = pm.Uniform(
        "logK",
        lower=0,
        upper=np.log(200),
        testval=np.log(0.5 * (np.max(rv) - np.min(rv))),
    )
    logP = pm.Uniform(
        "logP", lower=0, upper=np.log(10), testval=np.log(lit_period)
    )
    phi = pm.Uniform("phi", lower=0, upper=2 * np.pi, testval=0.1)

    # Parameterize the eccentricity using:
    #  h = sqrt(e) * sin(w)
    #  k = sqrt(e) * cos(w)
    hk = pmx.UnitDisk("hk", testval=np.array([0.01, 0.01]))
    e = pm.Deterministic("e", hk[0] ** 2 + hk[1] ** 2)
    w = pm.Deterministic("w", tt.arctan2(hk[1], hk[0]))

    rv0 = pm.Normal("rv0", mu=0.0, sd=10.0, testval=0.0)
    rvtrend = pm.Normal("rvtrend", mu=0.0, sd=10.0, testval=0.0)

    # Deterministic transformations
    n = 2 * np.pi * tt.exp(-logP)
    P = pm.Deterministic("P", tt.exp(logP))
    K = pm.Deterministic("K", tt.exp(logK))
    cosw = tt.cos(w)
    sinw = tt.sin(w)
    t0 = (phi + w) / n

    # The RV model
    bkg = pm.Deterministic("bkg", rv0 + rvtrend * t / 365.25)
    M = n * t - (phi + w)

    # This is the line that uses the custom Kepler solver
    f = xo.orbits.get_true_anomaly(M, e + tt.zeros_like(M))
    rvmodel = pm.Deterministic(
        "rvmodel", bkg + K * (cosw * (tt.cos(f) + e) - sinw * tt.sin(f))
    )

    # Condition on the observations
    pm.Normal("obs", mu=rvmodel, sd=rv_err, observed=rv)

    # Compute the phased RV signal
    phase = np.linspace(0, 1, 500)
    M_pred = 2 * np.pi * phase - (phi + w)
    f_pred = xo.orbits.get_true_anomaly(M_pred, e + tt.zeros_like(M_pred))
    rvphase = pm.Deterministic(
        "rvphase", K * (cosw * (tt.cos(f_pred) + e) - sinw * tt.sin(f_pred))
    )
```

In this case, I've found that it is useful to first optimize the parameters to find the "maximum a posteriori" (MAP) parameters and then start the sampler from there.
This is useful here because MCMC is not designed to *find* the maximum of the posterior; it's just meant to sample the shape of the posterior.
The performance of all MCMC methods can be really bad when the initialization isn't good (especially when some parameters are very well constrained).
To find the maximum a posteriori parameters using PyMC3, you can use the `optimize` function from [pymc3-ext](https://github.com/exoplanet-dev/pymc3-ext):

```{code-cell}
with model:
    map_params = pmx.optimize()
```

Let's make a plot to check that this initialization looks reasonable.
In the top plot, we're looking at the RV observations as a function of time with the initial guess for the long-term trend overplotted in blue.
In the lower panel, we plot the "folded" curve where we have wrapped the observations onto the best-fit period and the prediction for a single overplotted in orange.
If this doesn't look good, try adjusting the initial guesses for the parameters and see if you can get a better fit.

**Exercise:** Try changing the initial guesses for the parameters (as specified by the `testval` argument) and see how sensitive the results are to these values. Are there some parameters that are less important? Why is this?

```{code-cell}
fig, axes = plt.subplots(2, 1, figsize=(8, 8))

period = map_params["P"]

ax = axes[0]
ax.errorbar(t, rv, yerr=rv_err, fmt=".k")
ax.plot(t, map_params["bkg"], color="C0", lw=1)
ax.set_ylim(-110, 110)
ax.set_ylabel("radial velocity [m/s]")
ax.set_xlabel("time [days]")

ax = axes[1]
ax.errorbar(t % period, rv - map_params["bkg"], yerr=rv_err, fmt=".k")
ax.plot(phase * period, map_params["rvphase"], color="C1", lw=1)
ax.set_ylim(-110, 110)
ax.set_ylabel("radial velocity [m/s]")
ax.set_xlabel("phase [days]")

plt.tight_layout()
```

Now let's sample the posterior starting from our MAP estimate (here we're using the sampling routine from [pymc3-ext](https://github.com/exoplanet-dev/pymc3-ext) which wraps the PyMC3 function with some better defaults).

```{code-cell}
with model:
    trace = pmx.sample(
        draws=1000,
        tune=1000,
        start=map_params,
        chains=2,
        cores=2,
        target_accept=0.95,
        return_inferencedata=True,
    )
```

As above, it's always a good idea to take a look at the summary statistics for the chain.
If everything went as planned, there should be more than 1000 effective samples per chain and the Rhat values should be close to 1.
(Not too bad for about 30 seconds of run time!)

```{code-cell}
az.summary(
    trace,
    var_names=["logK", "logP", "phi", "e", "w", "rv0", "rvtrend"],
)
```

Similarly, we can make the corner plot again for this model.

```{code-cell}
_ = corner.corner(trace, var_names=["K", "P", "e", "w"])
```

Finally, the last plot that we'll make here is of the posterior predictive density.
In this case, this means that we want to look at the distribution of predicted models that are consistent with the data.
As above, the top plot shows the raw observations as black error bars and the RV trend model is overplotted in blue.
But, this time, the blue line is actually composed of 25 lines that are samples from the posterior over trends that are consistent with the data.
In the bottom panel, the orange lines indicate the same 25 posterior samples for the RV curve of one orbit.

```{code-cell}
fig, axes = plt.subplots(2, 1, figsize=(8, 8))

period = map_params["P"]

ax = axes[0]
ax.errorbar(t, rv, yerr=rv_err, fmt=".k")
ax.set_ylabel("radial velocity [m/s]")
ax.set_xlabel("time [days]")

ax = axes[1]
ax.errorbar(t % period, rv - map_params["bkg"], yerr=rv_err, fmt=".k")
ax.set_ylabel("radial velocity [m/s]")
ax.set_xlabel("phase [days]")

bkg = trace.posterior["bkg"].values
rvphase = trace.posterior["rvphase"].values

for ind in np.random.randint(np.prod(bkg.shape[:2]), size=25):
    i = np.unravel_index(ind, bkg.shape[:2])
    axes[0].plot(t, bkg[i], color="C0", lw=1, alpha=0.3)
    axes[1].plot(phase * period, rvphase[i], color="C1", lw=1, alpha=0.3)

axes[0].set_ylim(-110, 110)
axes[1].set_ylim(-110, 110)

plt.tight_layout()
```
