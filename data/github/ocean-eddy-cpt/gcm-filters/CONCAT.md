The contributor guidelines for GCM-Filters [can be found in our online documentation](https://gcm-filters.readthedocs.io/en/latest/how_to_contribute.html)
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

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at nora.loose@gmail.com. The project team will review and investigate all complaints, and will respond in a way that it deems appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
## GCM Filters

[![Tests](https://github.com/ocean-eddy-cpt/gcm-filters/workflows/Tests/badge.svg)](https://github.com/ocean-eddy-cpt/gcm-filters/actions?query=workflow%3ATests)
[![codecov](https://codecov.io/gh/ocean-eddy-cpt/gcm-filters/branch/master/graph/badge.svg?token=ZKRiulYe68)](https://codecov.io/gh/ocean-eddy-cpt/gcm-filters)
[![pre-commit](https://github.com/ocean-eddy-cpt/gcm-filters/workflows/pre-commit/badge.svg)](https://github.com/ocean-eddy-cpt/gcm-filters/actions?query=workflow%3Apre-commit)
[![Documentation Status](https://readthedocs.org/projects/gcm-filters/badge/?version=latest)](https://gcm-filters.readthedocs.io/en/latest/?badge=latest)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/gcm_filters.svg)](https://anaconda.org/conda-forge/gcm_filters)
[![PyPI version](https://badge.fury.io/py/gcm-filters.svg)](https://badge.fury.io/py/gcm-filters)
[![Downloads](https://pepy.tech/badge/gcm-filters)](https://pepy.tech/project/gcm-filters)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03947/status.svg)](https://doi.org/10.21105/joss.03947)

GCM-Filters: Diffusion-based Spatial Filtering of Gridded Data

### Description

**GCM-Filters** is a python package that performs spatial filtering analysis in a flexible and efficient way.
The GCM-Filters algorithm applies a discrete Laplacian to smooth a field through an iterative process that resembles diffusion ([Grooms et al., 2021](https://doi.org/10.1029/2021MS002552)).
The package can be used for either gridded observational data or gridded data that is produced by General Circulation Models (GCMs) of ocean, weather, and climate.
Such GCM data come on complex curvilinear grids, whose geometry is respected by the GCM-Filters Laplacians.
Through integration with [dask](https://dask.org/), GCM-Filters enables parallel, out-of-core filter analysis on both CPUs and GPUs.

### Installation

GCM-Filters can be installed using conda:
```shell
conda install -c conda-forge gcm_filters
```

GCM-Filters can also be installed with pip:
```shell
pip install gcm_filters
```

### Getting Started

To learn how to use GCM-Filters for your data, visit the [GCM-Filters documentation](https://gcm-filters.readthedocs.io/).


### Binder Demo

Click the button below to run an interactive demo of GCM-Filters in Binder:

[![badge](https://img.shields.io/badge/launch-binder-579aca.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAFkAAABZCAMAAABi1XidAAAB8lBMVEX///9XmsrmZYH1olJXmsr1olJXmsrmZYH1olJXmsr1olJXmsrmZYH1olL1olJXmsr1olJXmsrmZYH1olL1olJXmsrmZYH1olJXmsr1olL1olJXmsrmZYH1olL1olJXmsrmZYH1olL1olL0nFf1olJXmsrmZYH1olJXmsq8dZb1olJXmsrmZYH1olJXmspXmspXmsr1olL1olJXmsrmZYH1olJXmsr1olL1olJXmsrmZYH1olL1olLeaIVXmsrmZYH1olL1olL1olJXmsrmZYH1olLna31Xmsr1olJXmsr1olJXmsrmZYH1olLqoVr1olJXmsr1olJXmsrmZYH1olL1olKkfaPobXvviGabgadXmsqThKuofKHmZ4Dobnr1olJXmsr1olJXmspXmsr1olJXmsrfZ4TuhWn1olL1olJXmsqBi7X1olJXmspZmslbmMhbmsdemsVfl8ZgmsNim8Jpk8F0m7R4m7F5nLB6jbh7jbiDirOEibOGnKaMhq+PnaCVg6qWg6qegKaff6WhnpKofKGtnomxeZy3noG6dZi+n3vCcpPDcpPGn3bLb4/Mb47UbIrVa4rYoGjdaIbeaIXhoWHmZYHobXvpcHjqdHXreHLroVrsfG/uhGnuh2bwj2Hxk17yl1vzmljzm1j0nlX1olL3AJXWAAAAbXRSTlMAEBAQHx8gICAuLjAwMDw9PUBAQEpQUFBXV1hgYGBkcHBwcXl8gICAgoiIkJCQlJicnJ2goKCmqK+wsLC4usDAwMjP0NDQ1NbW3Nzg4ODi5+3v8PDw8/T09PX29vb39/f5+fr7+/z8/Pz9/v7+zczCxgAABC5JREFUeAHN1ul3k0UUBvCb1CTVpmpaitAGSLSpSuKCLWpbTKNJFGlcSMAFF63iUmRccNG6gLbuxkXU66JAUef/9LSpmXnyLr3T5AO/rzl5zj137p136BISy44fKJXuGN/d19PUfYeO67Znqtf2KH33Id1psXoFdW30sPZ1sMvs2D060AHqws4FHeJojLZqnw53cmfvg+XR8mC0OEjuxrXEkX5ydeVJLVIlV0e10PXk5k7dYeHu7Cj1j+49uKg7uLU61tGLw1lq27ugQYlclHC4bgv7VQ+TAyj5Zc/UjsPvs1sd5cWryWObtvWT2EPa4rtnWW3JkpjggEpbOsPr7F7EyNewtpBIslA7p43HCsnwooXTEc3UmPmCNn5lrqTJxy6nRmcavGZVt/3Da2pD5NHvsOHJCrdc1G2r3DITpU7yic7w/7Rxnjc0kt5GC4djiv2Sz3Fb2iEZg41/ddsFDoyuYrIkmFehz0HR2thPgQqMyQYb2OtB0WxsZ3BeG3+wpRb1vzl2UYBog8FfGhttFKjtAclnZYrRo9ryG9uG/FZQU4AEg8ZE9LjGMzTmqKXPLnlWVnIlQQTvxJf8ip7VgjZjyVPrjw1te5otM7RmP7xm+sK2Gv9I8Gi++BRbEkR9EBw8zRUcKxwp73xkaLiqQb+kGduJTNHG72zcW9LoJgqQxpP3/Tj//c3yB0tqzaml05/+orHLksVO+95kX7/7qgJvnjlrfr2Ggsyx0eoy9uPzN5SPd86aXggOsEKW2Prz7du3VID3/tzs/sSRs2w7ovVHKtjrX2pd7ZMlTxAYfBAL9jiDwfLkq55Tm7ifhMlTGPyCAs7RFRhn47JnlcB9RM5T97ASuZXIcVNuUDIndpDbdsfrqsOppeXl5Y+XVKdjFCTh+zGaVuj0d9zy05PPK3QzBamxdwtTCrzyg/2Rvf2EstUjordGwa/kx9mSJLr8mLLtCW8HHGJc2R5hS219IiF6PnTusOqcMl57gm0Z8kanKMAQg0qSyuZfn7zItsbGyO9QlnxY0eCuD1XL2ys/MsrQhltE7Ug0uFOzufJFE2PxBo/YAx8XPPdDwWN0MrDRYIZF0mSMKCNHgaIVFoBbNoLJ7tEQDKxGF0kcLQimojCZopv0OkNOyWCCg9XMVAi7ARJzQdM2QUh0gmBozjc3Skg6dSBRqDGYSUOu66Zg+I2fNZs/M3/f/Grl/XnyF1Gw3VKCez0PN5IUfFLqvgUN4C0qNqYs5YhPL+aVZYDE4IpUk57oSFnJm4FyCqqOE0jhY2SMyLFoo56zyo6becOS5UVDdj7Vih0zp+tcMhwRpBeLyqtIjlJKAIZSbI8SGSF3k0pA3mR5tHuwPFoa7N7reoq2bqCsAk1HqCu5uvI1n6JuRXI+S1Mco54YmYTwcn6Aeic+kssXi8XpXC4V3t7/ADuTNKaQJdScAAAAAElFTkSuQmCC)](https://mybinder.org/v2/gh/pangeo-data/pangeo-docker-images/812971e?urlpath=git-pull%3Frepo%3Dhttps%253A%252F%252Fgithub.com%252Focean-eddy-cpt%252Fgcm-filters%26urlpath%3Dlab%252Ftree%252Fgcm-filters%252Fdocs%252Fexamples%26branch%3Dmaster)


## Get in touch

Report bugs, suggest features or view the source code on [GitHub](https://github.com/ocean-eddy-cpt/gcm-filters).


## License and copyright

GCM-Filters is licensed under version 3 of the Gnu Lesser General Public License.

Development occurs on GitHub at <https://github.com/ocean-eddy-cpt/gcm-filters>.


## How to cite GCM-Filters

If you are using GCM-Filters and would like to cite it in academic publications, we
would certainly appreciate it. We recommend two citations.

- Loose et al., (2022). GCM-Filters: A Python Package for Diffusion-based Spatial Filtering of Gridded Data.
  Journal of Open Source Software, 7(70), 3947, https://doi.org/10.21105/joss.03947

  Here’s an example of a BibTeX entry:

  ```shell
  @article{Loose2022,
     author = {Nora Loose and Ryan Abernathey and Ian Grooms and Julius Busecke and Arthur Guillaumin and Elizabeth Yankovsky and Gustavo Marques and Jacob Steinberg and Andrew Slavin Ross and Hemant Khatri and Scott Bachman and Laure Zanna and Paige Martin},
     title = {GCM-Filters: A Python Package for Diffusion-based Spatial Filtering of Gridded Data},
     journal = {Journal of Open Source Software},
     volume = {7},
     number = {70},
     pages = {3947},
     doi = {10.21105/joss.03947},
     url = {https://doi.org/10.21105/joss.03947},
     year = {2022},
     publisher = {The Open Journal},
  }
  ```


- Grooms et al., (2021). Diffusion-Based Smoothers for Spatial Filtering of Gridded Geophysical Data.
  Journal of Advances in Modeling Earth Systems, 13, e2021MS002552, https://doi.org/10.1029/2021MS002552

  Here’s an example of a BibTeX entry:

  ```shell
  @article{Grooms2021,
     author = {Grooms, I. and Loose, N. and Abernathey, R. and Steinberg, J. M. and Bachman, S. D. and Marques, G. and Guillaumin, A. P. and Yankovsky, E.},
     title = {Diffusion-Based Smoothers for Spatial Filtering of Gridded Geophysical Data},
     journal = {Journal of Advances in Modeling Earth Systems},
     volume = {13},
     number = {9},
     pages = {e2021MS002552},
     keywords = {spatial filtering, coarse graining, data analysis},
     doi = {https://doi.org/10.1029/2021MS002552},
     url = {https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2021MS002552},
     year = {2021}
  }
  ```
---
title: 'GCM-Filters: A Python Package for Diffusion-based Spatial Filtering of Gridded Data'
tags:
  - Python
  - ocean modeling
  - climate modeling
  - fluid dynamics
authors:
  - name: Nora Loose
    orcid: 0000-0002-3684-9634
    affiliation: 1
  - name: Ryan Abernathey
    orcid: 0000-0001-5999-4917
    affiliation: 2
  - name: Ian Grooms
    orcid: 0000-0002-4678-7203
    affiliation: 1
  - name: Julius Busecke
    orcid: 0000-0001-8571-865X
    affiliation: 2
  - name: Arthur Guillaumin
    orcid: 0000-0003-1571-4228
    affiliation: 3
  - name: Elizabeth Yankovsky
    orcid: 0000-0003-3612-549X
    affiliation: 3
  - name: Gustavo Marques
    orcid: 0000-0001-7238-0290
    affiliation: 4
  - name: Jacob Steinberg
    orcid: 0000-0002-2609-6405
    affiliation: 5
  - name: Andrew Slavin Ross
    orcid: 0000-0002-2368-6979
    affiliation: 3
  - name: Hemant Khatri
    orcid: 0000-0001-6559-9059
    affiliation: 6
  - name: Scott Bachman
    orcid: 0000-0002-6479-4300
    affiliation: 4
  - name: Laure Zanna
    orcid: 0000-0002-8472-4828
    affiliation: 3
  - name: Paige Martin
    orcid: 0000-0003-3538-633X
    affiliation: "2, 7"
affiliations:
  - name: Department of Applied Mathematics, University of Colorado Boulder, Boulder, CO, USA
    index: 1
  - name: Lamont-Doherty Earth Observatory, Columbia University, New York, NY, USA
    index: 2
  - name: Courant Institute of Mathematical Sciences, New York University, New York, NY, USA
    index: 3
  - name: Climate and Global Dynamics Division, National Center for Atmospheric Research, Boulder, CO, USA
    index: 4
  - name: Woods Hole Oceanographic Institution, Woods Hole, MA, USA
    index: 5
  - name: Earth, Ocean and Ecological Sciences, University of Liverpool, UK
    index: 6
  - name: Research School of Earth Sciences, Australian National University, Canberra, Australia
    index: 7

date: 1 November 2021
bibliography: paper.bib

---

# Summary

`GCM-Filters` is a python package that allows scientists to perform spatial filtering analysis in an easy, flexible and efficient way. The package implements the filtering method based on the discrete Laplacian operator that was introduced by @grooms2021diffusion. The filtering algorithm is analogous to smoothing via diffusion; hence the name *diffusion-based filters*. `GCM-Filters` can be used with either gridded observational data or gridded data that is produced by General Circulation Models (GCMs) of ocean, weather, and climate. Spatial filtering of observational or GCM data is a common analysis method in the Earth Sciences, for example to study oceanic and atmospheric motions at different spatial scales or to develop subgrid-scale parameterizations for ocean models.

`GCM-Filters` provides filters that are highly configurable, with the goal to be useful for a wide range of scientific applications. The user has different options for selecting the filter scale and filter shape.
The filter scale can be defined in several ways: a fixed length scale (e.g., 100 km), a scale tied to a model grid scale (e.g., 1$^\circ$), or a scale tied to a varying dynamical scale (e.g., the Rossby radius of deformation). As an example, \autoref{fig1} shows unfiltered and filtered relative vorticity, where the filter scale is set to a model grid scale of 4$^\circ$. `GCM-Filters` also allows for anisotropic, i.e., direction-dependent, filtering.
Finally, the filter shape -- currently: either Gaussian or Taper -- determines how sharply the filter separates scales above and below the target filter scale.

![(Left) Snapshot of unfiltered surface relative vorticity  $\zeta = \partial_x v - \partial_y u$ from a global 0.1$^\circ$ simulation with MOM6 [@adcroft2019MOM6]. (Right) Relative vorticity filtered to 4$^\circ$, obtained by applying `GCM-Filters` to the field $\zeta$ on the left. The plots are made with `matplotlib` [@Hunter2007] and `cartopy` [@Cartopy].\label{fig1}](filtered_vorticity.png){ width=100% }

# Statement of Need

Spatial filtering is commonly used as a scientific tool for analyzing gridded data. An example of an existing spatial filtering tool in python is the `ndimage.gaussian_filter` function in `SciPy` [@2020SciPy-NMeth], implemented as a sequence of convolution filters. While being a valuable tool for image processing (or blurring), `SciPy`'s Gaussian filter is of limited use for GCM data; it assumes a regular and rectangular Cartesian grid, employs a simple boundary condition, and the definitions of filter scale and shape have little or no flexibility. The python package `GCM-Filters` is specificially designed to filter GCM data, and seeks to solve a number of challenges for the user:

1. GCM data comes on irregular curvilinear grids with spatially varying grid-cell geometry.
2. Continental boundaries require careful and special treatment when filtering ocean GCM output.
3. Earth Science applications benefit from configurable filters, where the definition of filter scale and shape is flexible.
4. GCM output is often too large to process in memory, requiring distributed and / or delayed execution.

The `GCM-Filters` algorithm [@grooms2021diffusion] applies a discrete Laplacian to smooth a field through an iterative process that resembles diffusion. The discrete Laplacian takes into account the varying grid-cell geometry and uses a no-flux boundary condition, mimicking how diffusion is internally implemented in GCMs. The no-flux boundary conditions ensures that the filter preserves the integral: $\int_{\Omega} \bar{f}(x,y) \,dA = \int_{\Omega} f (x,y)\, dA$, where $f$ is the original field, $\bar{f}$ the filtered field, and $\Omega$ the ocean domain. Conservation of the integral is a desirable filter property for many physical quantities, for example energy or ocean salinity. More details on the filter properties can be found in @grooms2021diffusion.

An important goal of `GCM-Filters` is to enable computationally efficient filtering. The user can employ `GCM-Filters` on either CPUs or GPUs, with `NumPy` [@harris2020array] or `CuPy` [@cupy2017learningsys] input data. `GCM-Filters` leverages `Dask` [@dask] and `Xarray` [@hoyer2017xarray] to support filtering of larger-than-memory datasets and computational flexibility.

# Usage

The main `GCM-Filters` class that the user will interface with is the `gcm_filters.Filter` object. When creating a filter object, the user specifies how they want to smooth their data, including the desired filter shape and filter scale. At this stage, the user also picks the grid type that matches their GCM data, given a predefined list of grid types. Each grid type has an associated discrete Laplacian, and requires different *grid variables* that the user must provide (the latter are usually available to the user as part of the GCM output). Currently, `GCM-Filters` provides a number of different grid types and associated discrete Laplacians:

* Grid types with **scalar Laplacians** that can be used for filtering scalar fields, for example temperature or vorticity (see \autoref{fig1}). The currently implemented grid types are compatible with different ocean GCM grids including MOM5 [@mom5], MOM6 [@adcroft2019MOM6] and the POP2 [@pop2] tripole grid.
* Grid types with **vector Laplacians** that can be used for filtering vector fields, for example horizontal velocity $(u,v)$. The currently implemented grid type is compatible with ocean GCM grids that use an Arakawa C-grid convention; examples include MOM6 [@adcroft2019MOM6] and the MITgcm [@mitgcm].

Atmospheric model grids are not yet supported, but could be implemented in `GCM-Filters`. Users are encouraged to contribute more grid types and Laplacians via pull requests.
While we are excited to share `GCM-Filters` at version `0.2.3`, we plan to continue improving and maintaining the package for the long run and welcome new contributors from the broader community.

# Acknowledgements

This work was supported by the National Science Foundation grants OCE 1912302, OCE 1912325, OCE 1912332, OCE 1912420, GEO 1912357, and the NOAA grant CVP NA19OAR4310364.
Busecke received support from the Gordon and Betty Moore Foundation.
This research is supported in part by the generosity of Eric and Wendy Schmidt by recommendation of Schmidt Futures, as part of its Virtual Earth System Research Institute (VESRI).

# References
---
name: Feature request
about: Suggest an idea for this project
title: "[FEATURE]"
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
name: Documentation request
about: Help us make the documentation better
title: "[DOC]"
labels: documentation
assignees: ''

---

**What were you looking for/trying to do?**:

**What was missing from the docs or not clear enough?**:

**How could the docs be improved? Additional examples, different wording, etc.**
---
name: Usage/General question
about: A general question about how to use GCM-Filters
title: ''
labels: question
assignees: ''

---
---
name: Bug report
about: Create a report to help us improve
title: "[BUG]"
labels: bug
assignees: ''

---

<!-- Please include a self-contained copy-pastable example that generates the issue if possible.

Please be concise with code posted. See guidelines below on how to provide a good bug report:

- Craft Minimal Bug Reports: http://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports
- Minimal Complete Verifiable Examples: https://stackoverflow.com/help/mcve

Bug reports that follow these guidelines are easier to diagnose, and so are often handled much more quickly.
-->

**What happened**:

**What you expected to happen**:

**Minimal Complete Verifiable Example (MCVE)**:

```python
# Put your MCVE code here
```

**Anything else we need to know?**:

**GCM-Filters version and environment**:

<details>
<!-- To get the version number of GCM-Filters do `import gcm_filters; print(gcm_filters.__version__)`-->

<!-- Paste the output of `conda list` from your shell here. -->

</details>
How to cite GCM-Filters
=======================

If you are using GCM-Filters and would like to cite it in academic publications, we
would certainly appreciate it. We recommend two citations.

     - Loose et al., (2022). GCM-Filters: A Python Package for Diffusion-based Spatial Filtering of Gridded Data.
       Journal of Open Source Software, 7(70), 3947, https://doi.org/10.21105/joss.03947

       Here’s an example of a BibTeX entry::

           @article{Loose2022,
              author = {Nora Loose and Ryan Abernathey and Ian Grooms and Julius Busecke and Arthur Guillaumin and Elizabeth Yankovsky and Gustavo Marques and Jacob Steinberg and Andrew Slavin Ross and Hemant Khatri and Scott Bachman and Laure Zanna and Paige Martin},
              title = {GCM-Filters: A Python Package for Diffusion-based Spatial Filtering of Gridded Data},
              journal = {Journal of Open Source Software},
              volume = {7},
              number = {70},
              pages = {3947},
              doi = {10.21105/joss.03947},
              url = {https://doi.org/10.21105/joss.03947},
              year = {2022},
              publisher = {The Open Journal},
           }

     - Grooms et al., (2021). Diffusion-Based Smoothers for Spatial Filtering of Gridded Geophysical Data.
       Journal of Advances in Modeling Earth Systems, 13, e2021MS002552, https://doi.org/10.1029/2021MS002552

       Here’s an example of a BibTeX entry::

           @article{Grooms2021,
              author = {Grooms, I. and Loose, N. and Abernathey, R. and Steinberg, J. M. and Bachman, S. D. and Marques, G. and Guillaumin, A. P. and Yankovsky, E.},
              title = {Diffusion-Based Smoothers for Spatial Filtering of Gridded Geophysical Data},
              journal = {Journal of Advances in Modeling Earth Systems},
              volume = {13},
              number = {9},
              pages = {e2021MS002552},
              keywords = {spatial filtering, coarse graining, data analysis},
              doi = {https://doi.org/10.1029/2021MS002552},
              url = {https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2021MS002552},
              year = {2021}
           }

Contributor Guide
=================

Step-by-step
------------
This guide will take you through the necessary steps in order to contribute code to the repository.


1. Setup
^^^^^^^^
If you are not super familiar with the terminology of `forking`, `pull request` etc, here is a `git tutorial <https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork>`_ to get you started.

Fork the main repository, and clone your fork onto the machine of your choice.

Make sure to install and activate the environment with::

   conda create -n gcm-filters-env -c conda-forge --file requirements.txt --file requirements-dev.txt

If you find that `conda` takes a long time, try to install `mamba` with ``conda install mamba`` and then
and do::

   mamba create -n gcm-filters-env -c conda-forge --file requirements.txt --file requirements-dev.txt


and finally activate the environment::

   conda activate gcm-filters-env

Before you can efficiently test your code you need to install this package in the now activated environment::

   pip install -e . --no-deps

For the linting also install `pre-commit <https://pre-commit.com>`_ now::

   pre-commit install

2. Write your code
^^^^^^^^^^^^^^^^^^
This one is pretty obvious 😁

3. Write and run tests
^^^^^^^^^^^^^^^^^^^^^^

Now figure out a good way to test the functionality of your code and run the test suite with::

   pytest

You will likely have to iterate here several times until all tests pass.

4. Linting and Docstrings
^^^^^^^^^^^^^^^^^^^^^^^^^
Once your tests pass, you want to make sure that all the code is formatted properly and has docstrings.

Run all the linters with::

   pre-commit run --all-files

Some things will automatically be reformatted, others need manual fixes. Follow the instructions in the terminal
until all checks pass.


Once you got everything to pass, you can stage and commit your changes and push them to the remote github repository.

How to change the documentation
-------------------------------

In order to build the documentation locally you should build and activate the docs environment::

   mamba env create -f docs/environment.yml

   conda activate gcm-filters-docs

Then navigate to the docs folder and build the docs locally with::

   cd docs

   make html

Once that is done you can open the created html files in `docs/_build/index.html` with your webbrowser::

   open _build/index.html

You can then edit and rebuild the docs until you are satisfied and submit a PR as usual.
Factoring the Gaussian Filter
=============================

This section provides background on applying a Gaussian filter with a large filter scale by repeatedly applying Gaussian filters with a smaller filter scale.
This was not discussed in `Grooms et al. (2021) <https://doi.org/10.1029/2021MS002552>`_, so this section provides extra detail.

The :py:class:`gcm_filters.Filter` class has an argument ``n_iterations`` that will automatically split a Gaussian filter into a sequence of Gaussian filters with smaller scales, each of which is less sensitive to roundoff errors. If a user is encountering instability with the standard Gaussian filter, the user can try setting ``n_iterations`` to an integer greater than 1. In general ``n_iterations`` should be set to the smallest value that avoids numerical instability. The user should input their desired filter scale, and then the code will automatically select a smaller filter scale for each of the constituent filters in such a way that the final result achieves the filter scale input by the user.

.. ipython:: python

    factored_gaussian_filter_x2 = gcm_filters.Filter(
        filter_scale=40,
        dx_min=1,
        n_iterations=2,  # number of constituent filters
        filter_shape=gcm_filters.FilterShape.GAUSSIAN,
        grid_type=gcm_filters.GridType.REGULAR,
    )
    factored_gaussian_filter_x2

.. ipython:: python

    factored_gaussian_filter_x3 = gcm_filters.Filter(
        filter_scale=40,
        dx_min=1,
        n_iterations=3,  # number of constituent filters
        filter_shape=gcm_filters.FilterShape.GAUSSIAN,
        grid_type=gcm_filters.GridType.REGULAR,
    )
    factored_gaussian_filter_x3


.. note:: When ``n_iterations`` is greater than 1, ``n_steps`` is the number of steps in a single small-scale Gaussian filter. The total number of steps taken by the algorithm is the product of ``n_iterations`` and ``n_steps``. If ``n_steps`` is not set by the user, the code automatically changes ``n_steps`` to a default value that gives an accurate approximation to the small-scale Gaussian filters that comprise the factored filter. In the example above, the code determined a smaller ``n_steps`` for the second filter compared to the first filter because the filter scale for each of the constituent filters is smaller.

Applying the Taper filter repeatedly is not equivalent to applying it once with a different scale, so this method only works for the Gaussian filter.

Mathematical Background
-----------------------

Per equation (13) of `Grooms et al. (2021) <https://doi.org/10.1029/2021MS002552>`_, the data to be filtered can be thought of as a vector :math:`\vec{f}` that can be expanded in an orthonormal basis of eigenfunctions of the discrete Laplacian

.. math:: \vec{f} = \sum_i \hat{f}_i\vec{q}_i

where :math:`\vec{q}_i` are the orthonormal basis vectors.
Each basis vector is associated with a wavenumber :math:`k_i` and a length scale :math:`2\pi/k_i`.
The filtered data can also be expanded in the same basis

.. math:: \sum_i \hat{f}_ip(k_i^2)\vec{q}_i

The polynomial :math:`p` approximates the target filter.
If you filter the data :math:`N` times using the same filter, then the result has the following expansion

.. math:: \sum_i \hat{f}_i(p(k_i^2))^N\vec{q}_i

For a Gaussian filter, the target filter is

.. math:: p(k^2) \approx g_L(k) = \text{exp}\left\{-\frac{L^2}{24}k^2\right\}

where :math:`L` is the ``filter_scale``.
The Gaussian target has the nice property that :math:`g_L^N = g_{\sqrt{N}L}`, i.e. if you apply the Gaussian filter :math:`N` times with scale :math:`L`, the result is the same as if you applied the Gaussian filter once with scale :math:`\sqrt{N}L`.
The way the code uses this is that instead of filtering once with scale :math:`L`, it has the option to filter :math:`N` times with scale :math:`L/\sqrt{N}`.

Inexact Equivalence
-------------------

Unfortunately, because we're using a polynomial approximation rather than an exact Gaussian filter, filtering once with scale :math:`L` is not exactly the same as filtering :math:`N` times with scale :math:`L/\sqrt{N}`.
The difference between the filtered field obtained using a single Gaussian versus :math:`N` factored Gaussians is exactly

.. math:: \sum_i \hat{f}_i(p_1(k_i^2) - (p_N(k^2))^N)\vec{q}_i

where :math:`p_1(k^2)` is the polynomial that approximates :math:`g_{L}(k)` and :math:`p_N(k^2)` is the polynomial that approximates :math:`g_{L/\sqrt{N}}(k)`.
We can bound the 2-norm of this error as follows.
For :math:`p_1` we have the error

.. math:: g_L(k) = p_1(k^2) + e_1(k)
    :label: unfactored-error

and for :math:`p_N` we have the error

.. math:: g_{L/\sqrt{N}}(k) = p_N(k^2) + e_N(k)

The aforementioned property of the Gaussian implies that

.. math:: g_{L}(k) = (p_N(k^2) + e_N(k))^N = (p_N(k^2))^N + N (p_N(k^2))^{N-1} e_N(k) + \ldots + (e_N(k))^N
    :label: factored-error

(using the binomial expansion.)
Subtracting :eq:`factored-error` from :eq:`unfactored-error` gives us an expression for the difference between the polynomial approximation with scale :math:`L` and the factored approximation using :math:`N` filters each with scale :math:`L/\sqrt{N}`:

.. math:: p_1(k^2) - (p_N(k^2))^N = - e_1(k) + N p_N(k^2)^{N-1} e_N(k) + \ldots + e_N(k)^N \sim - e_1(k) + N p_N(k^2)^{N-1} e_N(k)

where the last expression is in the limit of small errors :math:`|e_1(k)|` and :math:`|e_N(k)|` with :math:`N` fixed.
The difference in the two filtered fields is thus

.. math:: \sum_i \hat{f}_i(p_1(k_i^2) - (p_N(k_i^2))^N)\vec{q}_i\sim\sum_i \hat{f}_i(- e_1(k_i) + N p_N(k_i^2)^{N-1} e_N(k_i))\vec{q}_i

and the squared norm of this asymptotic approximation is exactly

.. math :: \sum_i \hat{f}_i^2(- e_1(k) + N p_N(k^2)^{N-1} e_N(k))^2.

The default choice of ``n_steps`` implies that :math:`|e_1(k)|` and :math:`|e_N(k)|` are both less than about 0.01, and the approximating polynomial is approximately bounded between 0 and 1.
Together these imply that

.. math :: (- e_1(k) + N p_N(k^2)^{N-1} e_N(k))^2 < 0.0001 (1+N)^2

The squared norm of the difference in the filtered fields is thus approximately bounded by

.. math :: 0.0001 (1+N)^2 \sum_i \hat{f}_i^2 = 0.0001(1+N)^2\|\vec{f}\|^2

The norm of the difference divided by the norm of the unfiltered field is thus approximately bounded by :math:`0.01(1+N)`.
This is why :math:`N` should be chosen as small as possible while avoiding numerical instability: as :math:`N` increases the difference between applying the filter once vs :math:`N` times increases.

Closing Comments
----------------

Note that the same ideas can be used to bound the norm of the difference between the filtered field that would be obtained using the exact filter :math:`g`, and the filtered field obtained using the polynomial approximation with :math:`N=1`.
In this case the analysis is simpler and the result is that the norm of the difference divided by the norm of the unfiltered field is bounded by 0.01.
Since this doesn't rely on factoring the filter, this bound is true for both the Gaussian and Taper filters.
Basic Filtering
================

The Filter Object
------------------


The core object in GCM-Filters is the :py:class:`gcm_filters.Filter` object. Its full documentation below enumerates all possible options.

.. autofunction:: gcm_filters.Filter


Details related to ``filter_scale``, ``filter_shape``, ``transition_width``, and ``n_steps`` can be found in the :doc:`theory`.
The following sections explain the options for ``grid_type`` and ``grid_vars`` in more detail.

Grid types
----------

To define a filter, we need to pick a grid and associated Laplacian that matches our data.
The currently implemented grid types are:

.. ipython:: python

    import gcm_filters
    list(gcm_filters.GridType)

This list will grow as we implement more Laplacians.

The following table provides an overview of these different grid type options: what grid they are suitable for, whether they handle land (i.e., continental boundaries), what boundary condition the Laplacian operators use, and whether they come with a scalar or vector Laplacian. You can also find links to example usages.

+--------------------------------+-----------------------------------------------------------------+--------------+--------------------+------------------+------------------------------------------+
| ``GridType``                   | Grid                                                            | Handles land | Boundary condition | Laplacian type   | Example                                  |
+================================+=================================================================+==============+====================+==================+==========================================+
| ``REGULAR``                    | Cartesian grid                                                  | no           | periodic           | Scalar Laplacian |                                          |
+--------------------------------+-----------------------------------------------------------------+--------------+--------------------+------------------+------------------------------------------+
| ``REGULAR_WITH_LAND``          | Cartesian grid                                                  | yes          | periodic           | Scalar Laplacian | see below                                |
+--------------------------------+-----------------------------------------------------------------+--------------+--------------------+------------------+------------------------------------------+
| ``IRREGULAR_WITH_LAND``        | locally orthogonal grid                                         | yes          | periodic           | Scalar Laplacian | :doc:`examples/example_filter_types`;    |
|                                |                                                                 |              |                    |                  | :doc:`examples/example_tripole_grid`     |
+--------------------------------+-----------------------------------------------------------------+--------------+--------------------+------------------+------------------------------------------+
| ``MOM5U``                      | Velocity-point on Arakawa B-Grid                                | yes          | periodic           | Scalar Laplacian |                                          |
+--------------------------------+-----------------------------------------------------------------+--------------+--------------------+------------------+------------------------------------------+
| ``MOM5T``                      | Tracer-point on Arakawa B-Grid                                  | yes          | periodic           | Scalar Laplacian |                                          |
+--------------------------------+-----------------------------------------------------------------+--------------+--------------------+------------------+------------------------------------------+
| ``TRIPOLAR_POP_WITH_LAND``     | locally orthogonal grid                                         | yes          | tripole            | Scalar Laplacian | :doc:`examples/example_tripole_grid`     |
+--------------------------------+-----------------------------------------------------------------+--------------+--------------------+------------------+------------------------------------------+
| ``VECTOR_C_GRID``              | `Arakawa C-Grid <https://en.wikipedia.org/wiki/Arakawa_grids>`_ | yes          | periodic           | Vector Laplacian | :doc:`examples/example_vector_laplacian` |
+--------------------------------+-----------------------------------------------------------------+--------------+--------------------+------------------+------------------------------------------+

Grid types with scalar Laplacians can be used for filtering scalar fields (such as temperature), and grid types with vector Laplacians can be used for filtering vector fields (such as velocity).

Grid types for simple fixed factor filtering
++++++++++++++++++++++++++++++++++++++++++++

The remaining grid types are for a special type of filtering: **simple fixed factor filtering** to achieve a fixed *coarsening* factor (see also the :doc:`theory`). If you specify one of the following grid types for your data, ``gcm_filters`` will internally transform your original (locally orthogonal) grid to a uniform Cartesian grid with `dx = dy = 1`, and perform fixed factor filtering on the uniform grid. After this is done, ``gcm_filters`` transforms the filtered field back to your original grid.
In practice, this coordinate transformation is achieved by area weighting and deweighting (see :doc:`theory`). This is why the following grid types have the suffix ``AREA_WEIGHTED``.

+-----------------------------------------------+-------------------------+--------------+--------------------+------------------+--------------------------------------+
| ``GridType``                                  | Grid                    | Handles land | Boundary condition | Laplacian type   | Example                              |
+===============================================+=========================+==============+====================+==================+======================================+
| ``REGULAR_AREA_WEIGHTED``                     | locally orthogonal grid | no           | periodic           | Scalar Laplacian |                                      |
+-----------------------------------------------+-------------------------+--------------+--------------------+------------------+--------------------------------------+
| ``REGULAR_WITH_LAND_AREA_WEIGHTED``           | locally orthogonal grid | yes          | periodic           | Scalar Laplacian | :doc:`examples/example_filter_types`;|
|                                               |                         |              |                    |                  | :doc:`examples/example_tripole_grid` |
+-----------------------------------------------+-------------------------+--------------+--------------------+------------------+--------------------------------------+
| ``TRIPOLAR_REGULAR_WITH_LAND_AREA_WEIGHTED``  | locally orthogonal grid | yes          | tripole            | Scalar Laplacian | :doc:`examples/example_tripole_grid` |
+-----------------------------------------------+-------------------------+--------------+--------------------+------------------+--------------------------------------+


Grid variables
--------------

Each grid type from the above two tables has different *grid variables* that must be provided as Xarray `DataArrays <http://xarray.pydata.org/en/stable/user-guide/data-structures.html>`_. For example, let's assume we are on a Cartesian grid (with uniform grid spacing equal to 1), and we want to use the grid type ``REGULAR_WITH_LAND``. To find out what the required grid variables for this grid type are, we can use this utility function.

.. ipython:: python

    gcm_filters.required_grid_vars(gcm_filters.GridType.REGULAR_WITH_LAND)

``wet_mask`` is a binary array representing the topography on our grid. Here the convention is that the array is 1 in the ocean (“wet points”) and 0 on land (“dry points”).

.. ipython:: python

    import numpy as np
    import xarray as xr

    ny, nx = (128, 256)

    mask_data = np.ones((ny, nx))
    mask_data[(ny // 4):(3 * ny // 4), (nx // 4):(3 * nx // 4)] = 0
    wet_mask = xr.DataArray(mask_data, dims=['y', 'x'])

.. ipython:: python
    :okwarning:

    @savefig wet_mask.png
    wet_mask.plot()

We have made a big island.


.. note:: Some more complicated grid types require more grid variables.
    The units for these variables should be *consistent*, but no specific system of units is required.
    For example, if grid cell edge lengths are defined using kilometers, then the filter scale and ``dx_min`` should also be defined using kilometers, and the grid cell areas should be defined in square kilometers.

Creating the Filter Object
--------------------------

We create a filter object as follows.

.. ipython:: python

    filter = gcm_filters.Filter(
        filter_scale=4,
        dx_min=1,
        filter_shape=gcm_filters.FilterShape.TAPER,
        grid_type=gcm_filters.GridType.REGULAR_WITH_LAND,
        grid_vars={'wet_mask': wet_mask}
    )
    filter

The string representation for the filter object in the last line includes some of the parameters it was initiliazed with, to help us keep track of what we are doing.
We have created a Taper filter that will filter our data by a fixed factor of 4.

Applying the Filter
-------------------

We can now apply the filter object that we created above to some data. Let's create a random 3D cube of data that matches our grid.

.. ipython:: python

    nt = 10
    data = np.random.rand(nt, ny, nx)
    da = xr.DataArray(data, dims=['time', 'y', 'x'])
    da

We now mask our data with the ``wet_mask``.

.. ipython:: python

   da_masked = da.where(wet_mask)

.. ipython:: python
    :okwarning:

    @savefig data.png
    da_masked.isel(time=0).plot()

Now that we have some data, we can apply our filter. We need to specify which dimension names to apply the filter over. In this case, it is ``y``, ``x``.

.. warning:: The dimension order matters! Since some filters deal
    with anisotropic grids, the latitude / y dimension must appear first
    in order to obtain the correct result. That is not an issue for this simple
    (isotropic) toy example but needs to be kept in mind for applications on
    real GCM grids.

.. ipython:: python

    %time da_filtered = filter.apply(da_masked, dims=['y', 'x'])

.. ipython:: python

    da_filtered

Let's visualize what the filter did.

.. ipython:: python
    :okwarning:

    @savefig data_filtered.png
    da_filtered.isel(time=0).plot()


Using Dask
-----------

Up to now, we have filtered *eagerly*; when we called ``.apply``, the results were computed immediately and stored in memory.
``GCM-Filters`` is also designed to work seamlessly with Dask array inputs. With `dask <https://dask.org/>`_, we can filter *lazily*, deferring the filter computations and possibly executing them in parallel.
We can do this with our synthetic data by converting them to dask.

.. ipython:: python
    :okwarning:

    da_dask = da_masked.chunk({'time': 2})
    da_dask

We now filter our data lazily.

.. ipython:: python
    :okwarning:

    da_filtered_lazy = filter.apply(da_dask, dims=['y', 'x'])
    da_filtered_lazy

Nothing has actually been computed yet.
We can trigger computation as follows:

.. ipython:: python

    %time da_filtered_computed = da_filtered_lazy.compute()

Here we got only a very modest speedup because our example data are too small. For bigger data, the performance benefit will be more evident.
.. gcm-filters documentation master file, created by
   sphinx-quickstart on Tue Jan 12 09:24:23 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

GCM-Filters: Diffusion-based Spatial Filtering of Gridded Data
===============================================================

.. image:: https://github.com/ocean-eddy-cpt/gcm-filters/workflows/Tests/badge.svg
   :target: https://github.com/ocean-eddy-cpt/gcm-filters/actions?query=workflow%3ATests

.. image:: https://codecov.io/gh/ocean-eddy-cpt/gcm-filters/branch/master/graph/badge.svg?token=ZKRiulYe68
   :target: https://codecov.io/gh/ocean-eddy-cpt/gcm-filters

.. image:: https://img.shields.io/conda/vn/conda-forge/gcm_filters.svg
   :target: https://anaconda.org/conda-forge/gcm_filters

.. image:: https://badge.fury.io/py/gcm-filters.svg
   :target: https://badge.fury.io/py/gcm-filters

.. image:: https://pepy.tech/badge/gcm-filters
   :target: https://pepy.tech/project/gcm-filters

.. image:: https://joss.theoj.org/papers/10.21105/joss.03947/status.svg
   :target: https://doi.org/10.21105/joss.03947

|
**GCM-Filters** is a python package that performs spatial filtering analysis in a flexible and efficient way.
The GCM-Filters algorithm applies a discrete Laplacian to smooth a field through an iterative process that resembles diffusion (see :doc:`theory` or `Grooms et al., 2021 <https://doi.org/10.1029/2021MS002552>`_).
The package can be used for either gridded observational data or gridded data that is produced by General Circulation Models (GCMs) of ocean, weather, and climate.
Such GCM data come on complex curvilinear grids, whose geometry is respected by the GCM-Filters Laplacians.
Through integration with `dask <https://dask.org/>`_, GCM-Filters enables parallel, out-of-core filter analysis on both CPUs and GPUs.

Getting Started
----------------

.. toctree::
   :maxdepth: 1

   installation
   theory
   basic_filtering
   gpu
   examples/example_filter_types
   examples/example_tripole_grid
   examples/example_vector_laplacian
   factored_gaussian
   examples/example_numerical_instability
   examples/example_satellite_observations

References
----------

.. toctree::
   :maxdepth: 1

   AMS 2022 Talk <https://noraloose.github.io/ams2022-talk/>
   api
   how_to_contribute
   how_to_cite


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _reStructuredText: https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html
Installation
------------

Requirements
^^^^^^^^^^^^

GCM-Filters is compatible with python 3. It requires xarray_, numpy_, and dask_.

Installation from conda forge
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

GCM-Filters can be installed via conda forge::

    conda install -c conda-forge gcm_filters


Installation from pip
^^^^^^^^^^^^^^^^^^^^^

GCM-Filters can also be installed with pip::

    pip install gcm_filters

This will install the latest release from
`pypi <https://pypi.python.org/pypi>`_.


Installation from GitHub
^^^^^^^^^^^^^^^^^^^^^^^^

GCM-Filters is under active development. To obtain the latest development version,
you may clone the `source repository <https://github.com/ocean-eddy-cpt/gcm-filters>`_
and install it::

    git clone https://github.com/ocean-eddy-cpt/gcm-filters.git
    cd gcm-filters
    python setup.py install

or simply::

    pip install git+https://github.com/ocean-eddy-cpt/gcm-filters.git

Users are encouraged to `fork <https://help.github.com/articles/fork-a-repo/>`_
GCM-Filters and submit issues_ and `pull requests`_.

.. _dask: http://dask.pydata.org
.. _numpy: https://numpy.org
.. _xarray: http://xarray.pydata.org
.. _issues: https://github.com/ocean-eddy-cpt/gcm-filters/issues
.. _`pull requests`: https://github.com/ocean-eddy-cpt/gcm-filters/pulls


How to run the example notebooks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to run the example notebooks in this documentation, you will need a few extra dependencies that you can install via::

   conda env create -f docs/environment.yml

   conda activate gcm-filters-docs
API Reference
=============

Filter Class
------------

.. autoclass:: gcm_filters.Filter
   :members:


Kernels
-------

.. automodule:: gcm_filters.kernels
   :members:
   :show-inheritance:
Filter Theory
=============

.. ipython:: python
    :suppress:

    import gcm_filters
    import numpy as np

The theory behind ``gcm-filters`` is described here at a high level.
For a more detailed treatment, see `Grooms et al. (2021) <https://doi.org/10.1029/2021MS002552>`_.

Filter Scale and Shape
----------------------

Any low-pass spatial filter should have a target length scale such that the filtered field keeps the part of the signal with length scales larger than the target length scale, and smoothes out smaller scales. In the context of this package the target length scale is called ``filter_scale``.

A spatial filter can also have a *shape* that determines how sharply it separates scales above and below the target length scale.
The filter shape can be thought of in terms of the kernel of a convolution filter

.. math:: \bar{f} = \int G(x - x')f(x') dx'

where :math:`f` is the function being filtered, :math:`G` is the filter kernel, and :math:`x'` is a dummy integration variable. (We note, however, that our filter is not *exactly* the same as a convolution filter. So our filter with a Gaussian target does not exactly produce the same results as a convolution against a Gaussian kernel on the sphere.)

This package currently has two filter shapes: ``GAUSSIAN`` and ``TAPER``.

.. ipython:: python

    list(gcm_filters.FilterShape)

For the ``GAUSSIAN`` filter the ``filter_scale`` equals :math:`\sqrt{12}\times` the standard deviation of the Gaussian.
\I.e. if you want to use a Gaussian filter with standard deviation L, then you should set ``filter_scale`` equal to L :math:`\times\sqrt{12}`.
This strange-seeming choice makes the Gaussian kernel have the same effect on large scales as a boxcar filter of width ``filter_scale``.
Thus ``filter_scale`` can be thought of as the "coarse grid scale" of the filtered field.

We can create a Gaussian filter as follows.

.. ipython:: python

    gaussian_filter = gcm_filters.Filter(
        filter_scale=4,
        dx_min=1,
        filter_shape=gcm_filters.FilterShape.GAUSSIAN,
        grid_type=gcm_filters.GridType.REGULAR,
    )
    gaussian_filter

Once the filter has been constructed, the method ``plot_shape`` can be used to plot the shape of the target filter and the approximate filter.

.. ipython:: python
    :okwarning:

    @savefig gaussian_shape.png
    gaussian_filter.plot_shape()

The distinction between the target filter and the approximate filter will be discussed below.

.. note:: ``plot_shape`` does not plot the shape of the filter *kernel*. Instead, it plots the frequency response of the filter for each wavenumber :math:`k`.
    In other words, the plot shows how the filter attenuates different scales in the data.
    Length scales are related to wavenumbers by :math:`\ell = 2\pi/k`.
    The filter leaves large scales unchanged, so the plot shows values close to 1 for small :math:`k`.
    The filter damps out small scales, so the plots shows values close to 0 for large :math:`k`.

The definition of the ``TAPER`` filter is more complex, but the ``filter_scale`` has the same meaning: it corresponds to the width of a qualitatively-similar boxcar filter.

.. ipython:: python

    taper_filter = gcm_filters.Filter(
        filter_scale=4,
        dx_min=1,
        filter_shape=gcm_filters.FilterShape.TAPER,
        grid_type=gcm_filters.GridType.REGULAR,
    )
    taper_filter

.. ipython:: python
    :okwarning:

    @savefig taper_shape.png
    taper_filter.plot_shape()


The plot above shows that the ``TAPER`` filter is more scale-selective than the Gaussian filter; it does a better job of leaving scales larger than the filter scale unchanged, and removing scales smaller than the filter scale.
The drawbacks of the ``TAPER`` filter are that it requires higher computational cost for the same filter scale (due to a higher number of necessary filter steps, see below), and it can produce negative values for the filtered field even when the unfiltered field is positive.

The Taper filter has a tunable parameter ``transition_width`` that controls how sharply the filter separates scales above and below the filter scale.
``transition_width`` = 1 would be the same as a complete *projection* onto the large scales, leaving the small scales completely zeroed out.
This would require a very high computational cost, and is not at all recommended!
The default is ``transition_width`` = :math:`\pi`.
Larger values for ``transition_width`` reduce the cost and the likelihood of producing negative values from positive data, but make the filter less scale-selective. In the example below, we choose ``transition_width`` = :math:`2\pi`.

.. ipython:: python

    wider_taper_filter = gcm_filters.Filter(
        filter_scale=4,
        dx_min=1,
        filter_shape=gcm_filters.FilterShape.TAPER,
        transition_width=2*np.pi,
        grid_type=gcm_filters.GridType.REGULAR,
    )
    wider_taper_filter

.. ipython:: python
    :okwarning:

    @savefig wider_taper_shape.png
    wider_taper_filter.plot_shape()

.. note:: The Taper filter is similar to the `Lanczos filter <https://journals.ametsoc.org/view/journals/apme/18/8/1520-0450_1979_018_1016_lfioat_2_0_co_2.xml>`_.
    Both are 1 for a range of large scales and 0 for a range of small scales, with a transition in between.
    The difference is in the transition region: in the transition region the Lanczos filter is straight line connecting 1 and 0, while the Taper filter is a smoother cubic.
    The Lanczos filter is typically described in terms of its "half-power cutoff wavelength"; the Taper filter can be similarly described.
    The half-power cutoff wavelength for the Taper filter with a ``filter_scale`` of :math:`L` and a ``transition_width`` of :math:`X` is :math:`2LX/(X+1)`.


Filter Steps
------------

The filter goes through several steps to produce the final filtered field.
There are two different kinds of steps: *Laplacian* and *Biharmonic* steps.
At each Laplacian step, the filtered field is updated using the following formula

.. math:: \bar{f} \leftarrow \bar{f} + \frac{1}{s_{j}}\Delta \bar{f}

The filtered field is initialized to :math:`\bar{f}=f` and :math:`\Delta` denotes a discrete Laplacian.
At each Biharmonic step, the filtered field is updated using

.. math:: \bar{f}\leftarrow \bar{f}+\frac{2R\{s_j\}}{|s_j|^2} \Delta\bar{f} + \frac{1}{|s_j|^2}\Delta^2\bar{f}

where :math:`R\{\cdot\}` denotes the real part of a complex number.

The total number of steps, ``n_steps``, and the values of :math:`s_j` are automatically selected by the code to produce the desired filter scale and shape.
If the filter scale is much larger than the grid scale, many steps are required.
Also, the Taper filter requires more steps than the Gaussian filter for the same ``filter_scale``; in the above examples the Taper filters required ``n_steps`` = 16, but the Gaussian filter only ``n_steps`` = 5.

The code allows users to set their own ``n_steps``.
Biharmonic steps are counted as 2 steps because their cost is approximately twice as much as a Laplacian step.
So with ``n_steps`` = 3 you might get one Laplacian plus one biharmonic step, or three Laplacian steps.
(The user cannot choose how ``n_steps`` is split between Laplacian and Biharmonic steps; that split is set internally in the code.)

For example, the user might want to use a smaller number of steps to reduce the cost. The caveat is that the accuracy will be reduced, so the filter might not act as expected: it may not have the right shape or the right length scale. To illustrate this, we create a new filter with a smaller number of steps than the default ``n_steps`` = 16, and plot the result.

.. ipython:: python
    :okwarning:

    taper_filter_8steps = gcm_filters.Filter(
        filter_scale=4,
        dx_min=1,
        filter_shape=gcm_filters.FilterShape.TAPER,
        n_steps=8,
        grid_type=gcm_filters.GridType.REGULAR,
    )
    taper_filter_8steps

.. ipython:: python
    :okwarning:

    @savefig taper_8steps_shape.png
    taper_filter_8steps.plot_shape()


The example above shows that using ``n_steps`` = 8 still yields a very accurate approximation of the target filter, at half the cost of the default. The main drawback in this example is that the filter slightly *amplifies* large scales, which also implies that it will not conserve variance.

The example below shows what happens with ``n_steps`` = 4.

.. ipython:: python
    :okwarning:

    taper_filter_4steps = gcm_filters.Filter(
        filter_scale=4,
        dx_min=1,
        filter_shape=gcm_filters.FilterShape.TAPER,
        n_steps=4,
        grid_type=gcm_filters.GridType.REGULAR,
    )
    taper_filter_4steps

.. ipython:: python
    :okwarning:

    @savefig taper_4steps_shape.png
    taper_filter_4steps.plot_shape()


.. warning::

    For this example of a Taper filter with a filter factor of 4, ``n_steps = 4`` is simply not enough to get a good approximation of the target filter. The ``taper_filter_4steps`` object created here will still "work" but it will not behave as expected; specifically, it will smooth more than expected - it will act like a filter with a larger filter scale.

The minimum number of steps is 3; if ``n_steps`` is not set by the user, or if it is set to a value less than 3, the code automatically changes ``n_steps`` to the default value.


Numerical Stability
-------------------

When the filter scale is much larger than the grid scale the filter can become unstable to roundoff errors.
The usual manifestation of these roundoff errors is high-amplitude small-scale noise in the filtered field.
(This problem is worse for the Taper filter than the Gaussian filter.)

.. tip::
    In such cases, the user has a few options to try to regain stability.

    1. If the data being filtered is single-precision, it might help to promote it to double precision (or higher) before filtering.
    2. If a user is encountering instability with the standard Gaussian filter, the user can try setting ``n_iterations`` to an integer greater than 1. Read more about this in :doc:`factored_gaussian`.
    3. The user can also try reducing ``n_steps``, but must not reduce it too much or the resulting filter will not behave as expected.
    4. Users might elect to *coarsen* their data before filtering, i.e. to reduce the resolution of the input data before applying the filter. This has the effect of increasing the grid size, and thus decreasing the gap between the filter scale and the grid scale.
    5. The final option is simply to use a different approach to filtering, not based on ``gcm-filters``.

:doc:`examples/example_numerical_instability` has an example of numerical instability, as well as examples of avoiding the instability by increasing the precision and coarsening the data.

Spatially-Varying Filter Scale
------------------------------

In the foregoing discussion the filter scale is fixed over the physical domain.
It is possible to vary the filter scale over the domain by introducing a *diffusivity* :math:`\kappa`.
(This diffusivity is nondimensional.)
The Laplacian steps are altered to

.. math:: \bar{f} \leftarrow \bar{f} + \frac{1}{s_{j}}\nabla\cdot(\kappa\nabla \bar{f})

and the Biharmonic steps are similarly altered by replacing :math:`\Delta` with :math:`\nabla\cdot(\kappa\nabla)`.
With :math:`\kappa` the *local* filter scale is :math:`\sqrt{\kappa}\times` ``filter_scale``.
For reasons given in `Grooms et al. (2021) <https://doi.org/10.1029/2021MS002552>`_, we require :math:`\kappa\le 1`, and at least one place in the domain where :math:`\kappa = 1`.
Thus, when using variable :math:`\kappa`, ``filter_scale`` sets the *largest* filter scale in the domain and the local filter scale can be reduced by making :math:`\kappa<1`.

Suppose, for example, that you want the local filter scale to be :math:`L(x,y)`.
You can achieve this in ``gcm-filters`` as follows.

1. Set ``filter_scale`` equal to the maximum of :math:`L(x,y)` over the domain. (Call this value :math:`L_{max}`).
2. Set :math:`\kappa` equal to :math:`L(x,y)^2/L_{max}^2`.

:doc:`examples/example_filter_types` has examples of filtering with spatially-varying filter scale.

Anisotropic Filtering
---------------------

It is possible to have different filter scales in different directions, and to have both the scales and directions vary over the domain.
This is achieved by replacing :math:`\kappa` in the previous section with a :math:`2\times2` symmetric and positive definite matrix (for a 2D domain), i.e. replacing :math:`\Delta` with :math:`\nabla\cdot(\mathbf{K}\nabla)`.
``gcm-filters`` currently only supports diagonal :math:`\mathbf{K}`, i.e. the principal axes of the anisotropic filter are aligned with the grid, so that the user only inputs one :math:`\kappa` for each grid direction, rather than a full :math:`2\times2` matrix.
Just like in the previous section, we require that each of these two :math:`\kappa` be less than or equal to 1, and the interpretation is also the same: the local filter scale in a particular direction is :math:`\sqrt{\kappa}\times` ``filter_scale``.

Suppose, for example, that you want to filter with a scale of 60 in the grid-x direction and a scale of 30 in the grid-y direction.
Then you would set ``filter_scale`` =  60, with :math:`\kappa_x = 1` to get a filter scale of 60 in the grid-x direction.
Next, to get a filter scale of 30 in the grid-y direction you would set :math:`\kappa_y=1/4`.

The :doc:`examples/example_filter_types` has examples of anisotropic filtering.

.. _Fixed factor filtering:

Fixed Factor Filtering
----------------------

:doc:`examples/example_filter_types` also shows methods designed specifically for the case where the user wants to set the local filter scale equal to a multiple :math:`m` of the local grid scale to achieve a fixed *coarsening* factor.
This can be achieved using the anisotropic diffusion described in the previous section.

An alternative way to achieve filtering with fixed coarsening factor :math:`m` is what we refer to as **simple fixed factor filtering**. This method is somewhat ad hoc, and *not* equivalent to fixed factor filtering via anisotropic diffusion. On the upside, simple fixed factor filtering is often significantly faster and yields very similar results in practice, as seen in :doc:`examples/example_filter_types`. The code handles simple fixed factor filtering as follows:

1. It multiplies the unfiltered data by the local grid cell area.
2. It applies a filter with ``filter_scale`` = :math:`m` *as if* the grid scale were uniform.
3. It divides the resulting field by the local grid cell area.

The first step is essentially a coordinate transformation where your original (locally orthogonal) grid is transformed to a uniform Cartesian grid with :math:`dx = dy = 1`. The third step is the reverse coordinate transformation.

.. note:: The three steps above are handled internally by ``gcm-filters`` if the user chooses one of the following grid types:

   * ``REGULAR_AREA_WEIGHTED``
   * ``REGULAR_WITH_LAND_AREA_WEIGHTED``
   * ``TRIPOLAR_REGULAR_WITH_LAND_AREA_WEIGHTED``

   together with ``filter_scale`` = :math:`m` and ``dx_min`` = 1. (For simple fixed factor filtering, only ``dx_min`` on the transformed uniform grid matters; and here we have ``dx_min`` = 1). Read more about the different grid types in :doc:`basic_filtering`.


Filtering Vectors
-----------------

In Cartesian geometry the Laplacian of a vector field can be obtained by taking the Laplacian of each component of the vector field, so vector fields can be filtered as described in the foregoing sections.
On smooth manifolds, the Laplacian of a vector field is not the same as the Laplacian of each component of the vector field.
Users may wish to use a **vector Laplacian** to filter vector fields.
The filter is constructed in exactly the same way; the only difference is in how the Laplacian is defined.
Rather than taking a scalar field and returning a scalar field, the vector Laplacian takes a vector field as input and returns a vector field.
To distinguish this from the scalar Laplacian, we refer to the filter based on a scalar Laplacian as a *diffusion-based* filter and the filter based on a vector Laplacian as a *viscosity-based* filter.
:doc:`examples/example_vector_laplacian` has examples of viscosity-based filtering.
