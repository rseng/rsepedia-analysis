About the Copyright Holders
===========================

*   Copyright (c) 2008-2011 AQR Capital Management, LLC

    AQR Capital Management began pandas development in 2008. Development was
    led by Wes McKinney. AQR released the source under this license in 2009.
*   Copyright (c) 2011-2012, Lambda Foundry, Inc.

    Wes is now an employee of Lambda Foundry, and remains the pandas project
    lead.
*   Copyright (c) 2011-2012, PyData Development Team

    The PyData Development Team is the collection of developers of the PyData
    project. This includes all of the PyData sub-projects, including pandas. The
    core team that coordinates development on GitHub can be found here:
    https://github.com/pydata.

Full credits for pandas contributors can be found in the documentation.

Our Copyright Policy
====================

PyData uses a shared copyright model. Each contributor maintains copyright
over their contributions to PyData. However, it is important to note that
these contributions are typically only changes to the repositories. Thus,
the PyData source code, in its entirety, is not the copyright of any single
person or institution. Instead, it is the collective copyright of the
entire PyData Development Team. If individual contributors want to maintain
a record of what changes/contributions they have specific copyright on,
they should indicate their copyright in the commit message of the change
when they commit the change to one of the PyData repositories.

With this in mind, the following banner should be used in any source code
file to indicate the copyright and license terms:

```
#-----------------------------------------------------------------------------
# Copyright (c) 2012, PyData Development Team
# All rights reserved.
#
# Distributed under the terms of the BSD Simplified License.
#
# The full license is in the LICENSE file, distributed with this software.
#-----------------------------------------------------------------------------
```

Other licenses can be found in the LICENSES directory.

License
=======

pandas is distributed under a 3-clause ("Simplified" or "New") BSD
license. Parts of NumPy, SciPy, numpydoc, bottleneck, which all have
BSD-compatible licenses, are included. Their licenses follow the pandas
license.
<div align="center">
  <img src="https://pandas.pydata.org/static/img/pandas.svg"><br>
</div>

-----------------

# pandas: powerful Python data analysis toolkit
[![PyPI Latest Release](https://img.shields.io/pypi/v/pandas.svg)](https://pypi.org/project/pandas/)
[![Conda Latest Release](https://anaconda.org/conda-forge/pandas/badges/version.svg)](https://anaconda.org/anaconda/pandas/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3509134.svg)](https://doi.org/10.5281/zenodo.3509134)
[![Package Status](https://img.shields.io/pypi/status/pandas.svg)](https://pypi.org/project/pandas/)
[![License](https://img.shields.io/pypi/l/pandas.svg)](https://github.com/pandas-dev/pandas/blob/main/LICENSE)
[![Azure Build Status](https://dev.azure.com/pandas-dev/pandas/_apis/build/status/pandas-dev.pandas?branch=main)](https://dev.azure.com/pandas-dev/pandas/_build/latest?definitionId=1&branch=main)
[![Coverage](https://codecov.io/github/pandas-dev/pandas/coverage.svg?branch=main)](https://codecov.io/gh/pandas-dev/pandas)
[![Downloads](https://static.pepy.tech/personalized-badge/pandas?period=month&units=international_system&left_color=black&right_color=orange&left_text=PyPI%20downloads%20per%20month)](https://pepy.tech/project/pandas)
[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/pydata/pandas)
[![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](https://numfocus.org)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)

## What is it?

**pandas** is a Python package that provides fast, flexible, and expressive data
structures designed to make working with "relational" or "labeled" data both
easy and intuitive. It aims to be the fundamental high-level building block for
doing practical, **real world** data analysis in Python. Additionally, it has
the broader goal of becoming **the most powerful and flexible open source data
analysis / manipulation tool available in any language**. It is already well on
its way towards this goal.

## Main Features
Here are just a few of the things that pandas does well:

  - Easy handling of [**missing data**][missing-data] (represented as
    `NaN`, `NA`, or `NaT`) in floating point as well as non-floating point data
  - Size mutability: columns can be [**inserted and
    deleted**][insertion-deletion] from DataFrame and higher dimensional
    objects
  - Automatic and explicit [**data alignment**][alignment]: objects can
    be explicitly aligned to a set of labels, or the user can simply
    ignore the labels and let `Series`, `DataFrame`, etc. automatically
    align the data for you in computations
  - Powerful, flexible [**group by**][groupby] functionality to perform
    split-apply-combine operations on data sets, for both aggregating
    and transforming data
  - Make it [**easy to convert**][conversion] ragged,
    differently-indexed data in other Python and NumPy data structures
    into DataFrame objects
  - Intelligent label-based [**slicing**][slicing], [**fancy
    indexing**][fancy-indexing], and [**subsetting**][subsetting] of
    large data sets
  - Intuitive [**merging**][merging] and [**joining**][joining] data
    sets
  - Flexible [**reshaping**][reshape] and [**pivoting**][pivot-table] of
    data sets
  - [**Hierarchical**][mi] labeling of axes (possible to have multiple
    labels per tick)
  - Robust IO tools for loading data from [**flat files**][flat-files]
    (CSV and delimited), [**Excel files**][excel], [**databases**][db],
    and saving/loading data from the ultrafast [**HDF5 format**][hdfstore]
  - [**Time series**][timeseries]-specific functionality: date range
    generation and frequency conversion, moving window statistics,
    date shifting and lagging


   [missing-data]: https://pandas.pydata.org/pandas-docs/stable/user_guide/missing_data.html
   [insertion-deletion]: https://pandas.pydata.org/pandas-docs/stable/user_guide/dsintro.html#column-selection-addition-deletion
   [alignment]: https://pandas.pydata.org/pandas-docs/stable/user_guide/dsintro.html?highlight=alignment#intro-to-data-structures
   [groupby]: https://pandas.pydata.org/pandas-docs/stable/user_guide/groupby.html#group-by-split-apply-combine
   [conversion]: https://pandas.pydata.org/pandas-docs/stable/user_guide/dsintro.html#dataframe
   [slicing]: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#slicing-ranges
   [fancy-indexing]: https://pandas.pydata.org/pandas-docs/stable/user_guide/advanced.html#advanced
   [subsetting]: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#boolean-indexing
   [merging]: https://pandas.pydata.org/pandas-docs/stable/user_guide/merging.html#database-style-dataframe-or-named-series-joining-merging
   [joining]: https://pandas.pydata.org/pandas-docs/stable/user_guide/merging.html#joining-on-index
   [reshape]: https://pandas.pydata.org/pandas-docs/stable/user_guide/reshaping.html
   [pivot-table]: https://pandas.pydata.org/pandas-docs/stable/user_guide/reshaping.html
   [mi]: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#hierarchical-indexing-multiindex
   [flat-files]: https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#csv-text-files
   [excel]: https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#excel-files
   [db]: https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#sql-queries
   [hdfstore]: https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#hdf5-pytables
   [timeseries]: https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#time-series-date-functionality

## Where to get it
The source code is currently hosted on GitHub at:
https://github.com/pandas-dev/pandas

Binary installers for the latest released version are available at the [Python
Package Index (PyPI)](https://pypi.org/project/pandas) and on [Conda](https://docs.conda.io/en/latest/).

```sh
# conda
conda install pandas
```

```sh
# or PyPI
pip install pandas
```

## Dependencies
- [NumPy - Adds support for large, multi-dimensional arrays, matrices and high-level mathematical functions to operate on these arrays](https://www.numpy.org)
- [python-dateutil - Provides powerful extensions to the standard datetime module](https://dateutil.readthedocs.io/en/stable/index.html)
- [pytz - Brings the Olson tz database into Python which allows accurate and cross platform timezone calculations](https://github.com/stub42/pytz)

See the [full installation instructions](https://pandas.pydata.org/pandas-docs/stable/install.html#dependencies) for minimum supported versions of required, recommended and optional dependencies.

## Installation from sources
To install pandas from source you need [Cython](https://cython.org/) in addition to the normal
dependencies above. Cython can be installed from PyPI:

```sh
pip install cython
```

In the `pandas` directory (same one where you found this file after
cloning the git repo), execute:

```sh
python setup.py install
```

or for installing in [development mode](https://pip.pypa.io/en/latest/cli/pip_install/#install-editable):


```sh
python -m pip install -e . --no-build-isolation --no-use-pep517
```

If you have `make`, you can also use `make develop` to run the same command.

or alternatively

```sh
python setup.py develop
```

See the full instructions for [installing from source](https://pandas.pydata.org/pandas-docs/stable/install.html#installing-from-source).

## License
[BSD 3](LICENSE)

## Documentation
The official documentation is hosted on PyData.org: https://pandas.pydata.org/pandas-docs/stable

## Background
Work on ``pandas`` started at [AQR](https://www.aqr.com/) (a quantitative hedge fund) in 2008 and
has been under active development since then.

## Getting Help

For usage questions, the best place to go to is [StackOverflow](https://stackoverflow.com/questions/tagged/pandas).
Further, general questions and discussions can also take place on the [pydata mailing list](https://groups.google.com/forum/?fromgroups#!forum/pydata).

## Discussion and Development
Most development discussions take place on GitHub in this repo. Further, the [pandas-dev mailing list](https://mail.python.org/mailman/listinfo/pandas-dev) can also be used for specialized discussions or design issues, and a [Gitter channel](https://gitter.im/pydata/pandas) is available for quick development related questions.

## Contributing to pandas [![Open Source Helpers](https://www.codetriage.com/pandas-dev/pandas/badges/users.svg)](https://www.codetriage.com/pandas-dev/pandas)

All contributions, bug reports, bug fixes, documentation improvements, enhancements, and ideas are welcome.

A detailed overview on how to contribute can be found in the **[contributing guide](https://pandas.pydata.org/docs/dev/development/contributing.html)**.

If you are simply looking to start working with the pandas codebase, navigate to the [GitHub "issues" tab](https://github.com/pandas-dev/pandas/issues) and start looking through interesting issues. There are a number of issues listed under [Docs](https://github.com/pandas-dev/pandas/issues?labels=Docs&sort=updated&state=open) and [good first issue](https://github.com/pandas-dev/pandas/issues?labels=good+first+issue&sort=updated&state=open) where you could start out.

You can also triage issues which may include reproducing bug reports, or asking for vital information such as version numbers or reproduction instructions. If you would like to start triaging issues, one easy way to get started is to [subscribe to pandas on CodeTriage](https://www.codetriage.com/pandas-dev/pandas).

Or maybe through using pandas you have an idea of your own or are looking for something in the documentation and thinking ‘this can be improved’...you can do something about it!

Feel free to ask questions on the [mailing list](https://groups.google.com/forum/?fromgroups#!forum/pydata) or on [Gitter](https://gitter.im/pydata/pandas).

As contributors and maintainers to this project, you are expected to abide by pandas' code of conduct. More information can be found at: [Contributor Code of Conduct](https://github.com/pandas-dev/pandas/blob/main/.github/CODE_OF_CONDUCT.md)
Release Notes
=============

The list of changes to Pandas between each release can be found
[here](https://pandas.pydata.org/pandas-docs/stable/whatsnew/index.html). For full
details, see the commit logs at https://github.com/pandas-dev/pandas.
Directory containing the pandas website (hosted at https://pandas.pydata.org).

The website sources are in `web/pandas/`, which also include a `config.yml` file
containing the settings to build the website. The website is generated with the
command `./pandas_web.py pandas`. See `./pandas_web.py --help` and the header of
the script for more information and options.

After building the website, to navigate it, it is needed to access the web using
an http server (a not open the local files with the browser, since the links and
the image sources are absolute to where they are served from). The easiest way
to run an http server locally is to run `python -m http.server` from the
`web/build/` directory.
# Donate to pandas

<div id="salsalabs-donate-container">
</div>
<script type="text/javascript"
        src="https://default.salsalabs.org/api/widget/template/4ba4e328-1855-47c8-9a89-63e4757d2151/?tId=salsalabs-donate-container">
</script>

_pandas_ is a Sponsored Project of [NumFOCUS](https://numfocus.org/), a 501(c)(3) nonprofit charity in the United States.
NumFOCUS provides _pandas_ with fiscal, legal, and administrative support to help ensure the
health and sustainability of the project. Visit numfocus.org for more information.

Donations to _pandas_ are managed by NumFOCUS. For donors in the United States, your gift is tax-deductible
to the extent provided by law. As with any donation, you should consult with your tax adviser about your particular tax situation.
# Try pandas online

<section>
    <pre data-executable>
import pandas
fibonacci = pandas.Series([1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144])
fibonacci.sum()
    </pre>
    <script src="https://combinatronics.com/ines/juniper/v0.1.0/dist/juniper.min.js"></script>
    <script>new Juniper({ repo: 'datapythonista/pandas-web' })</script>
</section>

## Interactive tutorials

You can also try _pandas_ on [Binder](https://mybinder.org/) for one of the next topics:

- Exploratory analysis of US presidents
- Preprocessing the Titanic dataset to train a machine learning model
- Forecasting the stock market

_(links will be added soon)_
# Contribute to pandas

_pandas_ is and will always be **free**. To make the development sustainable, we need _pandas_ users, corporate
and individual, to support the development by providing their time and money.

You can find more information about current developers in the [team page](about/team.html),
and about current sponsors in the [sponsors page](about/sponsors.html).

<section>
    <div class="container mt-5">
      <div class="row text-center">
        <div class="col-md-4">
          <span class="fa-stack fa-4x">
            <i class="fas fa-circle fa-stack-2x pink"></i>
            <i class="fas fa-building fa-stack-1x fa-inverse"></i>
          </span>
          <h4 class="service-heading mt-3 fw-bold blue">Corporate support</h4>
          <p class="text-muted">
            pandas depends on companies and institutions using the software to support its development. Hiring
            people to work on pandas, or letting existing employees to contribute to the
            software. Or sponsoring pandas with funds, so the project can hire people to
            progress on the <a href="about/roadmap.html">pandas roadmap</a>.
          </p>
          <p>More information in the <a href="about/sponsors.html">sponsors page</a></p>
        </div>
        <div class="col-md-4">
          <span class="fa-stack fa-4x">
            <i class="fas fa-circle fa-stack-2x pink"></i>
            <i class="fas fa-users fa-stack-1x fa-inverse"></i>
          </span>
          <h4 class="service-heading mt-3 fw-bold blue">Individual contributors</h4>
          <p class="text-muted">
            pandas is mostly developed by volunteers. All kind of contributions are welcome,
            such as contributions to the code, to the website (including graphical designers),
            to the documentation (including translators) and others. There are tasks for all
            levels, including beginners.
          </p>
          <p>More information in the <a href="{{ base_url }}/docs/development/index.html">contributing page</a></p>
        </div>
        <div class="col-md-4">
          <span class="fa-stack fa-4x">
            <i class="fas fa-circle fa-stack-2x pink"></i>
            <i class="fas fa-dollar-sign fa-stack-1x fa-inverse"></i>
          </span>
          <h4 class="service-heading mt-3 fw-bold blue">Donations</h4>
          <p class="text-muted">
            Individual donations are appreciated, and are used for things like the project
            infrastructure, travel expenses for our volunteer contributors to attend
            the in-person sprints, or to give small grants to develop features.
          </p>
          <p>Make your donation in the <a href="donate.html">donate page</a></p>
        </div>
      </div>
    </div>
</section>
# Getting started

## Installation instructions

The next steps provides the easiest and recommended way to set up your
environment to use pandas. Other installation options can be found in
the [advanced installation page]({{ base_url}}/docs/getting_started/install.html).

1. Download [Anaconda](https://www.anaconda.com/distribution/) for your operating system and
   the latest Python version, run the installer, and follow the steps. Please note:

    - It is not needed (and discouraged) to install Anaconda as root or administrator.
    - When asked if you wish to initialize Anaconda3, answer yes.
    - Restart the terminal after completing the installation.

    Detailed instructions on how to install Anaconda can be found in the
    [Anaconda documentation](https://docs.anaconda.com/anaconda/install/).

2. In the Anaconda prompt (or terminal in Linux or MacOS), start JupyterLab:

    <img class="img-fluid" alt="" src="{{ base_url }}/static/img/install/anaconda_prompt.png"/>

3. In JupyterLab, create a new (Python 3) notebook:

    <img class="img-fluid" alt="" src="{{ base_url }}/static/img/install/jupyterlab_home.png"/>

4. In the first cell of the notebook, you can import pandas and check the version with:

    <img class="img-fluid" alt="" src="{{ base_url }}/static/img/install/pandas_import_and_version.png"/>

5. Now you are ready to use pandas, and you can write your code in the next cells.

## Tutorials

You can learn more about pandas in the [tutorials]({{ base_url }}/docs/getting_started/intro_tutorials/),
and more about JupyterLab in the
[JupyterLab documentation](https://jupyterlab.readthedocs.io/en/stable/user/interface.html).

## Books

The book we recommend to learn pandas is [Python for Data Analysis](https://amzn.to/2KI5JJw),
by [Wes McKinney](https://wesmckinney.com/), creator of pandas.

<a href="https://amzn.to/2KI5JJw">
    <img alt="Python for Data Analysis" src="{{ base_url }}/static/img/pydata_book.gif"/>
</a>

## Videos

<iframe width="560" height="315" frameborder="0"
src="https://www.youtube.com/embed/_T8LGqJtuGc"
allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
allowfullscreen></iframe>

## Cheat sheet

[pandas cheat sheet](https://pandas.pydata.org/Pandas_Cheat_Sheet.pdf)
# About pandas

## History of development

In 2008, _pandas_ development began at [AQR Capital Management](https://www.aqr.com).
By the end of 2009 it had been [open sourced](https://en.wikipedia.org/wiki/Open_source),
and is actively supported today by a community of like-minded individuals around the world who
contribute their valuable time and energy to help make open source _pandas_
possible. Thank you to [all of our contributors](team.html).

Since 2015, _pandas_ is a [NumFOCUS sponsored project](https://numfocus.org/sponsored-projects).
This will help ensure the success of development of _pandas_ as a world-class open-source project.

### Timeline

- **2008**:  Development of _pandas_ started
- **2009**: _pandas_ becomes open source
- **2012**: First edition of _Python for Data Analysis_ is published
- **2015**: _pandas_ becomes a [NumFOCUS sponsored project](https://numfocus.org/sponsored-projects)
- **2018**: First in-person core developer sprint

## Library Highlights

- A fast and efficient **DataFrame** object for data manipulation with
  integrated indexing;

- Tools for **reading and writing data** between in-memory data structures and
  different formats: CSV and text files, Microsoft Excel, SQL databases, and
  the fast HDF5 format;

- Intelligent **data alignment** and integrated handling of **missing data**:
  gain automatic label-based alignment in computations and easily manipulate
  messy data into an orderly form;

- Flexible **reshaping** and pivoting of data sets;

- Intelligent label-based **slicing**, **fancy indexing**, and **subsetting**
  of large data sets;

- Columns can be inserted and deleted from data structures for **size
  mutability**;

- Aggregating or transforming data with a powerful **group by** engine
  allowing split-apply-combine operations on data sets;

- High performance **merging and joining** of data sets;

- **Hierarchical axis indexing** provides an intuitive way of working with
  high-dimensional data in a lower-dimensional data structure;

- **Time series**-functionality: date range generation and frequency
  conversion, moving window statistics, date shifting and lagging.
  Even create domain-specific time offsets and join time
  series without losing data;

- Highly **optimized for performance**, with critical code paths written in
  [Cython](http://www.cython.org/) or C.

- Python with *pandas* is in use in a wide variety of **academic and
  commercial** domains, including Finance, Neuroscience, Economics,
  Statistics, Advertising, Web Analytics, and more.

## Mission

_pandas_ aims to be the fundamental high-level building block for doing practical,
real world data analysis in Python.
Additionally, it has the broader goal of becoming the most powerful and flexible
open source data analysis / manipulation tool available in any language.

## Vision

A world where data analytics and manipulation software is:

- Accessible to everyone
- Free for users to use and modify
- Flexible
- Powerful
- Easy to use
- Fast

## Values

Is in the core of _pandas_ to be respectful and welcoming with everybody,
users, contributors and the broader community. Regardless of level of experience,
gender, gender identity and expression, sexual orientation, disability,
personal appearance, body size, race, ethnicity, age, religion, or nationality.
# Citing and logo

## Citing pandas

If you use _pandas_ for a scientific publication, we would appreciate citations to the published software and the
following paper:

- [pandas on Zenodo](https://zenodo.org/record/3715232#.XoqFyC2ZOL8),
   Please find us on Zenodo and replace with the citation for the version you are using. You can replace the full author
   list from there with "The pandas development team" like in the example below.

        @software{reback2020pandas,
            author       = {The pandas development team},
            title        = {pandas-dev/pandas: Pandas},
            month        = feb,
            year         = 2020,
            publisher    = {Zenodo},
            version      = {latest},
            doi          = {10.5281/zenodo.3509134},
            url          = {https://doi.org/10.5281/zenodo.3509134}
        }

- [Data structures for statistical computing in python](https://conference.scipy.org/proceedings/scipy2010/pdfs/mckinney.pdf),
   McKinney, Proceedings of the 9th Python in Science Conference, Volume 445, 2010.

        @InProceedings{ mckinney-proc-scipy-2010,
          author    = { {W}es {M}c{K}inney },
          title     = { {D}ata {S}tructures for {S}tatistical {C}omputing in {P}ython },
          booktitle = { {P}roceedings of the 9th {P}ython in {S}cience {C}onference },
          pages     = { 56 - 61 },
          year      = { 2010 },
          editor    = { {S}t\'efan van der {W}alt and {J}arrod {M}illman },
          doi       = { 10.25080/Majora-92bf1922-00a }
        }

## Brand and logo

When using the project name _pandas_, please use it in lower case, even at the beginning of a sentence.

The official logos of _pandas_ are:

### Primary logo

<table class="table logo">
    <tr>
        <td>
            <img alt="" src="{{ base_url }}/static/img/pandas.svg"/>
        </td>
        <td style="background-color: #150458">
            <img alt="" src="{{ base_url }}/static/img/pandas_white.svg"/>
        </td>
    </tr>
</table>

### Secondary logo

<table class="table logo">
    <tr>
        <td>
            <img alt="" src="{{ base_url }}/static/img/pandas_secondary.svg"/>
        </td>
        <td style="background-color: #150458">
            <img alt="" src="{{ base_url }}/static/img/pandas_secondary_white.svg"/>
        </td>
    </tr>
</table>

### Logo mark

<table class="table logo">
    <tr>
        <td>
            <img alt="" src="{{ base_url }}/static/img/pandas_mark.svg"/>
        </td>
        <td style="background-color: #150458">
            <img alt="" src="{{ base_url }}/static/img/pandas_mark_white.svg"/>
        </td>
    </tr>
</table>

### Logo usage

The pandas logo is available in full color and white accent.
The full color logo should only appear against white backgrounds.
The white accent logo should go against contrasting color background.

When using the logo, please follow the next directives:

- Primary logo should never be seen under 1 inch in size for printing and 72px for web
- The secondary logo should never be seen under 0.75 inch in size for printing and 55px for web
- Leave enough margin around the logo (leave the height of the logo in the top, bottom and both sides)
- Do not distort the logo by changing its proportions
- Do not place text or other elements on top of the logo

### Colors

<table class="table">
    <tr>
        <td style="text-align: center;">
            <svg xmlns="http://www.w3.org/2000/svg" width="100" height="100">
                <circle cx="50" cy="50" r="50" fill="#150458"/>
            </svg>
            <br/>
            <b style="color: #150458;">Blue</b><br/>
            RGB: R21 G4 B88<br/>
            HEX: #150458
        </td>
        <td style="text-align: center;">
            <svg xmlns="http://www.w3.org/2000/svg" width="100" height="100">
                <circle cx="50" cy="50" r="50" fill="#ffca00"/>
            </svg>
            <br/>
            <b style="color: #150458;">Yellow</b><br/>
            RGB: R255 G202 B0<br/>
            HEX: #FFCA00
        </td>
        <td style="text-align: center;">
            <svg xmlns="http://www.w3.org/2000/svg" width="100" height="100">
                <circle cx="50" cy="50" r="50" fill="#e70488"/>
            </svg>
            <br/>
            <b style="color: #150458;">Pink</b><br/>
            RGB: R231 G4 B136<br/>
            HEX: #E70488
        </td>
    </tr>
</table>
# Sponsors

## NumFOCUS

![](https://numfocus.org/wp-content/uploads/2018/01/optNumFocus_LRG.png)

_pandas_ is a Sponsored Project of [NumFOCUS](https://numfocus.org/), a 501(c)(3) nonprofit charity in the United States.
NumFOCUS provides _pandas_ with fiscal, legal, and administrative support to help ensure the
health and sustainability of the project. Visit numfocus.org for more information.

Donations to _pandas_ are managed by NumFOCUS. For donors in the United States, your gift is tax-deductible
to the extent provided by law. As with any donation, you should consult with your tax adviser about your particular tax situation.

## Become a sponsor

As a free and open source project, _pandas_ relies on the support of the community of users for its development.
If you work for an organization that uses and benefits from _pandas_, please consider supporting pandas. There
are different ways, such as employing people to work on pandas, funding the project, or becoming a
[NumFOCUS sponsor](https://numfocus.org/sponsors) to support the broader ecosystem. Please contact us at
[admin@numfocus.org](mailto:admin@numfocus.org) to discuss.

## Institutional partners

Institutional partners are companies and universities that support the project by employing contributors.
Current institutional partners include:

<ul>
    {% for company in sponsors.active if company.kind == "partner" %}
        <li><a href="{{ company.url }}">{{ company.name }}</a>: {{ company.description }}</li>
    {% endfor %}
</ul>

## Sponsors

Sponsors are organizations that provide funding for pandas. Current sponsors include:

<ul>
    {% for company in sponsors.active if company.kind == "regular" %}
        <li><a href="{{ company.url }}">{{ company.name }}</a>: {{ company.description }}</li>
    {% endfor %}
</ul>

## In-kind sponsors

In-kind sponsors are organizations that support pandas development with goods or services.
Current in-kind sponsors include:

<ul>
    {% for company in sponsors.inkind %}
        <li><a href="{{ company.url }}">{{ company.name }}</a>: {{ company.description }}</li>
    {% endfor %}
</ul>

## Past institutional partners

<ul>
    {% for company in sponsors.past if company.kind == "partner" %}
        <li><a href="{{ company.url }}">{{ company.name }}</a></li>
    {% endfor %}
</ul>
# Team

## Contributors

_pandas_ is made with love by more than [2,000 volunteer contributors](https://github.com/pandas-dev/pandas/graphs/contributors).

If you want to support pandas development, you can find information in the [donations page](../donate.html).

## Maintainers

<div class="card-group maintainers">
    {% for person in maintainers.people %}
        <div class="card">
            <img class="card-img-top" alt="" src="{{ person.avatar_url }}"/>
            <div class="card-body">
                <h6 class="card-title">
                    {% if person.blog %}
                        <a href="{{ person.blog }}">
                            {{ person.name or person.login }}
                        </a>
                    {% else %}
                        {{ person.name or person.login }}
                    {% endif %}
                </h6>
                <p class="card-text small"><a href="{{ person.html_url }}">{{ person.login }}</a></p>
            </div>
        </div>
    {% endfor %}
</div>

## Diversity and Inclusion

> _pandas_ expressly welcomes and encourages contributions from anyone who faces under-representation, discrimination in the technology industry
> or anyone willing to increase the diversity of our team.
> We have identified visible gaps and obstacles in sustaining diversity and inclusion in the open-source communities and we are proactive in increasing
> the diversity of our team.
> We have a [code of conduct](../community/coc.html) to ensure a friendly and welcoming environment.
> Please send an email to [pandas-code-of-conduct-committee](mailto:pandas-coc@googlegroups.com), if you think we can do a
> better job at achieving this goal.

## Governance

Wes McKinney is the Benevolent Dictator for Life (BDFL).

The project governance is available in the [project governance documents](https://github.com/pandas-dev/pandas-governance).

## Code of conduct committee

<ul>
    {% for person in maintainers.coc %}
        <li>{{ person }}</li>
    {% endfor %}
</ul>

## NumFOCUS committee

<ul>
    {% for person in maintainers.numfocus %}
        <li>{{ person }}</li>
    {% endfor %}
</ul>

## Emeritus maintainers

<ul>
    {% for person in maintainers.emeritus %}
        <li>{{ person }}</li>
    {% endfor %}
</ul>
# Roadmap

This page provides an overview of the major themes in pandas'
development. Each of these items requires a relatively large amount of
effort to implement. These may be achieved more quickly with dedicated
funding or interest from contributors.

An item being on the roadmap does not mean that it will *necessarily*
happen, even with unlimited funding. During the implementation period we
may discover issues preventing the adoption of the feature.

Additionally, an item *not* being on the roadmap does not exclude it
from inclusion in pandas. The roadmap is intended for larger,
fundamental changes to the project that are likely to take months or
years of developer time. Smaller-scoped items will continue to be
tracked on our [issue tracker](https://github.com/pandas-dev/pandas/issues).

See [Roadmap evolution](#roadmap-evolution) for proposing
changes to this document.

## Extensibility

Pandas `extending.extension-types` allow
for extending NumPy types with custom data types and array storage.
Pandas uses extension types internally, and provides an interface for
3rd-party libraries to define their own custom data types.

Many parts of pandas still unintentionally convert data to a NumPy
array. These problems are especially pronounced for nested data.

We'd like to improve the handling of extension arrays throughout the
library, making their behavior more consistent with the handling of
NumPy arrays. We'll do this by cleaning up pandas' internals and
adding new methods to the extension array interface.

## String data type

Currently, pandas stores text data in an `object` -dtype NumPy array.
The current implementation has two primary drawbacks: First, `object`
-dtype is not specific to strings: any Python object can be stored in an
`object` -dtype array, not just strings. Second: this is not efficient.
The NumPy memory model isn't especially well-suited to variable width
text data.

To solve the first issue, we propose a new extension type for string
data. This will initially be opt-in, with users explicitly requesting
`dtype="string"`. The array backing this string dtype may initially be
the current implementation: an `object` -dtype NumPy array of Python
strings.

To solve the second issue (performance), we'll explore alternative
in-memory array libraries (for example, Apache Arrow). As part of the
work, we may need to implement certain operations expected by pandas
users (for example the algorithm used in, `Series.str.upper`). That work
may be done outside of pandas.

## Apache Arrow interoperability

[Apache Arrow](https://arrow.apache.org) is a cross-language development
platform for in-memory data. The Arrow logical types are closely aligned
with typical pandas use cases.

We'd like to provide better-integrated support for Arrow memory and
data types within pandas. This will let us take advantage of its I/O
capabilities and provide for better interoperability with other
languages and libraries using Arrow.

## Block manager rewrite

We'd like to replace pandas current internal data structures (a
collection of 1 or 2-D arrays) with a simpler collection of 1-D arrays.

Pandas internal data model is quite complex. A DataFrame is made up of
one or more 2-dimensional "blocks", with one or more blocks per dtype.
This collection of 2-D arrays is managed by the BlockManager.

The primary benefit of the BlockManager is improved performance on
certain operations (construction from a 2D array, binary operations,
reductions across the columns), especially for wide DataFrames. However,
the BlockManager substantially increases the complexity and maintenance
burden of pandas.

By replacing the BlockManager we hope to achieve

-   Substantially simpler code
-   Easier extensibility with new logical types
-   Better user control over memory use and layout
-   Improved micro-performance
-   Option to provide a C / Cython API to pandas' internals

See [these design
documents](https://dev.pandas.io/pandas2/internal-architecture.html#removal-of-blockmanager-new-dataframe-internals)
for more.

## Decoupling of indexing and internals

The code for getting and setting values in pandas' data structures
needs refactoring. In particular, we must clearly separate code that
converts keys (e.g., the argument to `DataFrame.loc`) to positions from
code that uses these positions to get or set values. This is related to
the proposed BlockManager rewrite. Currently, the BlockManager sometimes
uses label-based, rather than position-based, indexing. We propose that
it should only work with positional indexing, and the translation of
keys to positions should be entirely done at a higher level.

Indexing is a complicated API with many subtleties. This refactor will
require care and attention. More details are discussed at
<https://github.com/pandas-dev/pandas/wiki/(Tentative)-rules-for-restructuring-indexing-code>

## Numba-accelerated operations

[Numba](https://numba.pydata.org) is a JIT compiler for Python code.
We'd like to provide ways for users to apply their own Numba-jitted
functions where pandas accepts user-defined functions (for example,
`Series.apply`,
`DataFrame.apply`,
`DataFrame.applymap`, and in groupby and
window contexts). This will improve the performance of
user-defined-functions in these operations by staying within compiled
code.

## Documentation improvements

We'd like to improve the content, structure, and presentation of the
pandas documentation. Some specific goals include

-   Overhaul the HTML theme with a modern, responsive design
    (`15556`)
-   Improve the "Getting Started" documentation, designing and writing
    learning paths for users different backgrounds (e.g. brand new to
    programming, familiar with other languages like R, already familiar
    with Python).
-   Improve the overall organization of the documentation and specific
    subsections of the documentation to make navigation and finding
    content easier.

## Performance monitoring

Pandas uses [airspeed velocity](https://asv.readthedocs.io/en/stable/)
to monitor for performance regressions. ASV itself is a fabulous tool,
but requires some additional work to be integrated into an open source
project's workflow.

The [asv-runner](https://github.com/asv-runner) organization, currently
made up of pandas maintainers, provides tools built on top of ASV. We
have a physical machine for running a number of project's benchmarks,
and tools managing the benchmark runs and reporting on results.

We'd like to fund improvements and maintenance of these tools to

-   Be more stable. Currently, they're maintained on the nights and
    weekends when a maintainer has free time.
-   Tune the system for benchmarks to improve stability, following
    <https://pyperf.readthedocs.io/en/latest/system.html>
-   Build a GitHub bot to request ASV runs *before* a PR is merged.
    Currently, the benchmarks are only run nightly.

## Roadmap Evolution

Pandas continues to evolve. The direction is primarily determined by
community interest. Everyone is welcome to review existing items on the
roadmap and to propose a new item.

Each item on the roadmap should be a short summary of a larger design
proposal. The proposal should include

1.  Short summary of the changes, which would be appropriate for
    inclusion in the roadmap if accepted.
2.  Motivation for the changes.
3.  An explanation of why the change is in scope for pandas.
4.  Detailed design: Preferably with example-usage (even if not
    implemented yet) and API documentation
5.  API Change: Any API changes that may result from the proposal.

That proposal may then be submitted as a GitHub issue, where the pandas
maintainers can review and comment on the design. The [pandas mailing
list](https://mail.python.org/mailman/listinfo/pandas-dev) should be
notified of the proposal.

When there's agreement that an implementation would be welcome, the
roadmap should be updated to include the summary and a link to the
discussion issue.
# Code of conduct

As contributors and maintainers of this project, and in the interest of
fostering an open and welcoming community, we pledge to respect all people who
contribute through reporting issues, posting feature requests, updating
documentation, submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free
experience for everyone, regardless of level of experience, gender, gender
identity and expression, sexual orientation, disability, personal appearance,
body size, race, ethnicity, age, religion, or nationality.

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery
* Personal attacks
* Trolling or insulting/derogatory comments
* Public or private harassment
* Publishing other's private information, such as physical or electronic
  addresses, without explicit permission
* Other unethical or unprofessional conduct

Furthermore, we encourage inclusive behavior - for example,
please don’t say “hey guys!” but “hey everyone!”.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

By adopting this Code of Conduct, project maintainers commit themselves to
fairly and consistently applying these principles to every aspect of managing
this project. Project maintainers who do not follow or enforce the Code of
Conduct may be permanently removed from the project team.

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community.

A working group of community members is committed to promptly addressing any
reported issues. The working group is made up of pandas contributors and users.
Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the working group by e-mail (pandas-coc@googlegroups.com).
Messages sent to this e-mail address will not be publicly visible but only to
the working group members. The working group currently includes

<ul>
    {% for person in maintainers.coc %}
    <li>{{ person }}</li>
    {% endfor %}
</ul>

All complaints will be reviewed and investigated and will result in a response
that is deemed necessary and appropriate to the circumstances. Maintainers are
obligated to maintain confidentiality with regard to the reporter of an
incident.

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 1.3.0, available at
[https://www.contributor-covenant.org/version/1/3/0/][version],
and the [Swift Code of Conduct][swift].

[homepage]: https://www.contributor-covenant.org
[version]: https://www.contributor-covenant.org/version/1/3/0/
[swift]: https://swift.org/community/#code-of-conduct
# Ecosystem

Increasingly, packages are being built on top of pandas to address
specific needs in data preparation, analysis and visualization. This is
encouraging because it means pandas is not only helping users to handle
their data tasks but also that it provides a better starting point for
developers to build powerful and more focused data tools. The creation
of libraries that complement pandas' functionality also allows pandas
development to remain focused around its original requirements.

This is an inexhaustive list of projects that build on pandas in order
to provide tools in the PyData space. For a list of projects that depend
on pandas, see the [libraries.io usage page for
pandas](https://libraries.io/pypi/pandas/usage) or [search pypi for
pandas](https://pypi.org/search/?q=pandas).

We'd like to make it easier for users to find these projects, if you
know of other substantial projects that you feel should be on this list,
please let us know.

## Statistics and machine learning

### [Statsmodels](https://www.statsmodels.org/)

Statsmodels is the prominent Python "statistics and econometrics
library" and it has a long-standing special relationship with pandas.
Statsmodels provides powerful statistics, econometrics, analysis and
modeling functionality that is out of pandas' scope. Statsmodels
leverages pandas objects as the underlying data container for
computation.

### [sklearn-pandas](https://github.com/paulgb/sklearn-pandas)

Use pandas DataFrames in your [scikit-learn](https://scikit-learn.org/)
ML pipeline.

### [Featuretools](https://github.com/alteryx/featuretools/)

Featuretools is a Python library for automated feature engineering built
on top of pandas. It excels at transforming temporal and relational
datasets into feature matrices for machine learning using reusable
feature engineering "primitives". Users can contribute their own
primitives in Python and share them with the rest of the community.

### [Compose](https://github.com/alteryx/compose)

Compose is a machine learning tool for labeling data and prediction engineering.
It allows you to structure the labeling process by parameterizing
prediction problems and transforming time-driven relational data into
target values with cutoff times that can be used for supervised learning.

## Visualization

### [Altair](https://altair-viz.github.io/)

Altair is a declarative statistical visualization library for Python.
With Altair, you can spend more time understanding your data and its
meaning. Altair's API is simple, friendly and consistent and built on
top of the powerful Vega-Lite JSON specification. This elegant
simplicity produces beautiful and effective visualizations with a
minimal amount of code. Altair works with Pandas DataFrames.

### [Bokeh](https://bokeh.pydata.org)

Bokeh is a Python interactive visualization library for large datasets
that natively uses the latest web technologies. Its goal is to provide
elegant, concise construction of novel graphics in the style of
Protovis/D3, while delivering high-performance interactivity over large
data to thin clients.

[Pandas-Bokeh](https://github.com/PatrikHlobil/Pandas-Bokeh) provides a
high level API for Bokeh that can be loaded as a native Pandas plotting
backend via

```
pd.set_option("plotting.backend", "pandas_bokeh")
```

It is very similar to the matplotlib plotting backend, but provides
interactive web-based charts and maps.

### [seaborn](https://seaborn.pydata.org)

Seaborn is a Python visualization library based on
[matplotlib](https://matplotlib.org). It provides a high-level,
dataset-oriented interface for creating attractive statistical graphics.
The plotting functions in seaborn understand pandas objects and leverage
pandas grouping operations internally to support concise specification
of complex visualizations. Seaborn also goes beyond matplotlib and
pandas with the option to perform statistical estimation while plotting,
aggregating across observations and visualizing the fit of statistical
models to emphasize patterns in a dataset.

### [plotnine](https://github.com/has2k1/plotnine/)

Hadley Wickham's [ggplot2](https://ggplot2.tidyverse.org/) is a
foundational exploratory visualization package for the R language. Based
on ["The Grammar of
Graphics"](https://www.cs.uic.edu/~wilkinson/TheGrammarOfGraphics/GOG.html)
it provides a powerful, declarative and extremely general way to
generate bespoke plots of any kind of data.
Various implementations to other languages are available.
A good implementation for Python users is [has2k1/plotnine](https://github.com/has2k1/plotnine/).

### [IPython Vega](https://github.com/vega/ipyvega)

[IPython Vega](https://github.com/vega/ipyvega) leverages [Vega](https://github.com/vega/vega) to create plots within Jupyter Notebook.

### [Plotly](https://plot.ly/python)

[Plotly's](https://plot.ly/) [Python API](https://plot.ly/python/)
enables interactive figures and web shareability. Maps, 2D, 3D, and
live-streaming graphs are rendered with WebGL and
[D3.js](https://d3js.org/). The library supports plotting directly from
a pandas DataFrame and cloud-based collaboration. Users of [matplotlib,
ggplot for Python, and
Seaborn](https://plot.ly/python/matplotlib-to-plotly-tutorial/) can
convert figures into interactive web-based plots. Plots can be drawn in
[IPython Notebooks](https://plot.ly/ipython-notebooks/) , edited with R
or MATLAB, modified in a GUI, or embedded in apps and dashboards. Plotly
is free for unlimited sharing, and has
[cloud](https://plot.ly/product/plans/),
[offline](https://plot.ly/python/offline/), or
[on-premise](https://plot.ly/product/enterprise/) accounts for private
use.

### [QtPandas](https://github.com/draperjames/qtpandas)

Spun off from the main pandas library, the
[qtpandas](https://github.com/draperjames/qtpandas) library enables
DataFrame visualization and manipulation in PyQt4 and PySide
applications.

## IDE

### [IPython](https://ipython.org/documentation.html)

IPython is an interactive command shell and distributed computing
environment. IPython tab completion works with Pandas methods and also
attributes like DataFrame columns.

### [Jupyter Notebook / Jupyter Lab](https://jupyter.org)

Jupyter Notebook is a web application for creating Jupyter notebooks. A
Jupyter notebook is a JSON document containing an ordered list of
input/output cells which can contain code, text, mathematics, plots and
rich media. Jupyter notebooks can be converted to a number of open
standard output formats (HTML, HTML presentation slides, LaTeX, PDF,
ReStructuredText, Markdown, Python) through 'Download As' in the web
interface and `jupyter convert` in a shell.

Pandas DataFrames implement `_repr_html_`and `_repr_latex` methods which
are utilized by Jupyter Notebook for displaying (abbreviated) HTML or
LaTeX tables. LaTeX output is properly escaped. (Note: HTML tables may
or may not be compatible with non-HTML Jupyter output formats.)

See `Options and Settings <options>` and
`Available Options <options.available>`
for pandas `display.` settings.

### [quantopian/qgrid](https://github.com/quantopian/qgrid)

qgrid is "an interactive grid for sorting and filtering DataFrames in
IPython Notebook" built with SlickGrid.

### [Spyder](https://www.spyder-ide.org/)

Spyder is a cross-platform PyQt-based IDE combining the editing,
analysis, debugging and profiling functionality of a software
development tool with the data exploration, interactive execution, deep
inspection and rich visualization capabilities of a scientific
environment like MATLAB or Rstudio.

Its [Variable
Explorer](https://docs.spyder-ide.org/variableexplorer.html) allows
users to view, manipulate and edit pandas `Index`, `Series`, and
`DataFrame` objects like a "spreadsheet", including copying and
modifying values, sorting, displaying a "heatmap", converting data
types and more. Pandas objects can also be renamed, duplicated, new
columns added, copyed/pasted to/from the clipboard (as TSV), and
saved/loaded to/from a file. Spyder can also import data from a variety
of plain text and binary files or the clipboard into a new pandas
DataFrame via a sophisticated import wizard.

Most pandas classes, methods and data attributes can be autocompleted in
Spyder's [Editor](https://docs.spyder-ide.org/editor.html) and [IPython
Console](https://docs.spyder-ide.org/ipythonconsole.html), and Spyder's
[Help pane](https://docs.spyder-ide.org/help.html) can retrieve and
render Numpydoc documentation on pandas objects in rich text with Sphinx
both automatically and on-demand.

## API

### [pandas-datareader](https://github.com/pydata/pandas-datareader)

`pandas-datareader` is a remote data access library for pandas
(PyPI:`pandas-datareader`). It is based on functionality that was
located in `pandas.io.data` and `pandas.io.wb` but was split off in
v0.19. See more in the [pandas-datareader
docs](https://pandas-datareader.readthedocs.io/en/latest/):

The following data feeds are available:

- Google Finance
- Tiingo
- Morningstar
- IEX
- Robinhood
- Enigma
- Quandl
- FRED
- Fama/French
- World Bank
- OECD
- Eurostat
- TSP Fund Data
- Nasdaq Trader Symbol Definitions
- Stooq Index Data
- MOEX Data

### [quandl/Python](https://github.com/quandl/Python)

Quandl API for Python wraps the Quandl REST API to return Pandas
DataFrames with timeseries indexes.

### [pydatastream](https://github.com/vfilimonov/pydatastream)

PyDatastream is a Python interface to the [Thomson Dataworks Enterprise
(DWE/Datastream)](http://dataworks.thomson.com/Dataworks/Enterprise/1.0/)
SOAP API to return indexed Pandas DataFrames with financial data. This
package requires valid credentials for this API (non free).

### [pandaSDMX](https://pandasdmx.readthedocs.io)

pandaSDMX is a library to retrieve and acquire statistical data and
metadata disseminated in [SDMX](https://www.sdmx.org) 2.1, an
ISO-standard widely used by institutions such as statistics offices,
central banks, and international organisations. pandaSDMX can expose
datasets and related structural metadata including data flows,
code-lists, and data structure definitions as pandas Series or
MultiIndexed DataFrames.

### [fredapi](https://github.com/mortada/fredapi)

fredapi is a Python interface to the [Federal Reserve Economic Data
(FRED)](https://fred.stlouisfed.org/) provided by the Federal Reserve
Bank of St. Louis. It works with both the FRED database and ALFRED
database that contains point-in-time data (i.e. historic data
revisions). fredapi provides a wrapper in Python to the FRED HTTP API,
and also provides several convenient methods for parsing and analyzing
point-in-time data from ALFRED. fredapi makes use of pandas and returns
data in a Series or DataFrame. This module requires a FRED API key that
you can obtain for free on the FRED website.

## Domain specific

### [Geopandas](https://github.com/kjordahl/geopandas)

Geopandas extends pandas data objects to include geographic information
which support geometric operations. If your work entails maps and
geographical coordinates, and you love pandas, you should take a close
look at Geopandas.

### [xarray](https://github.com/pydata/xarray)

xarray brings the labeled data power of pandas to the physical sciences
by providing N-dimensional variants of the core pandas data structures.
It aims to provide a pandas-like and pandas-compatible toolkit for
analytics on multi-dimensional arrays, rather than the tabular data for
which pandas excels.

## Out-of-core

### [Blaze](https://blaze.pydata.org/)

Blaze provides a standard API for doing computations with various
in-memory and on-disk backends: NumPy, Pandas, SQLAlchemy, MongoDB,
PyTables, PySpark.

### [Dask](https://dask.readthedocs.io/en/latest/)

Dask is a flexible parallel computing library for analytics. Dask
provides a familiar `DataFrame` interface for out-of-core, parallel and
distributed computing.

### [Dask-ML](https://dask-ml.readthedocs.io/en/latest/)

Dask-ML enables parallel and distributed machine learning using Dask
alongside existing machine learning libraries like Scikit-Learn,
XGBoost, and TensorFlow.

### [Koalas](https://koalas.readthedocs.io/en/latest/)

Koalas provides a familiar pandas DataFrame interface on top of Apache
Spark. It enables users to leverage multi-cores on one machine or a
cluster of machines to speed up or scale their DataFrame code.

### [Odo](http://odo.pydata.org)

Odo provides a uniform API for moving data between different formats. It
uses pandas own `read_csv` for CSV IO and leverages many existing
packages such as PyTables, h5py, and pymongo to move data between non
pandas formats. Its graph based approach is also extensible by end users
for custom formats that may be too specific for the core of odo.

### [Ray](https://ray.readthedocs.io/en/latest/pandas_on_ray.html)

Pandas on Ray is an early stage DataFrame library that wraps Pandas and
transparently distributes the data and computation. The user does not
need to know how many cores their system has, nor do they need to
specify how to distribute the data. In fact, users can continue using
their previous Pandas notebooks while experiencing a considerable
speedup from Pandas on Ray, even on a single machine. Only a
modification of the import statement is needed, as we demonstrate below.
Once you've changed your import statement, you're ready to use Pandas on
Ray just like you would Pandas.

```
# import pandas as pd
import ray.dataframe as pd
```

### [Vaex](https://docs.vaex.io/)

Increasingly, packages are being built on top of pandas to address
specific needs in data preparation, analysis and visualization. Vaex is
a python library for Out-of-Core DataFrames (similar to Pandas), to
visualize and explore big tabular datasets. It can calculate statistics
such as mean, sum, count, standard deviation etc, on an N-dimensional
grid up to a billion (10^9^) objects/rows per second. Visualization is
done using histograms, density plots and 3d volume rendering, allowing
interactive exploration of big data. Vaex uses memory mapping, zero
memory copy policy and lazy computations for best performance (no memory
wasted).

- ``vaex.from_pandas``
- ``vaex.to_pandas_df``

## Data cleaning and validation

### [pyjanitor](https://github.com/ericmjl/pyjanitor/)

Pyjanitor provides a clean API for cleaning data, using method chaining.

### [Engarde](https://engarde.readthedocs.io/en/latest/)

Engarde is a lightweight library used to explicitly state your
assumptions about your datasets and check that they're *actually* true.

## Extension data types

Pandas provides an interface for defining
`extension types <extending.extension-types>` to extend NumPy's type system. The following libraries
implement that interface to provide types not found in NumPy or pandas,
which work well with pandas' data containers.

### [cyberpandas](https://cyberpandas.readthedocs.io/en/latest)

Cyberpandas provides an extension type for storing arrays of IP
Addresses. These arrays can be stored inside pandas' Series and
DataFrame.

### [Pandas-Genomics](https://pandas-genomics.readthedocs.io/en/latest/)

Pandas-Genomics provides an extension type and extension array for working
 with genomics data.  It also includes `genomics` accessors for many useful properties
 and methods related to QC and analysis of genomics data.

### [Pint-Pandas](https://github.com/hgrecco/pint-pandas)

Pint-Pandas provides an extension type for storing numeric arrays with units.
These arrays can be stored inside pandas' Series and DataFrame. Operations
between Series and DataFrame columns which use pint's extension array are then
units aware.

## Accessors

A directory of projects providing
`extension accessors <extending.register-accessors>`. This is for users to discover new accessors and for library
authors to coordinate on the namespace.

  | Library                                                              | Accessor   |  Classes              |
  | ---------------------------------------------------------------------|------------|-----------------------|
  | [cyberpandas](https://cyberpandas.readthedocs.io/en/latest)          | `ip`       | `Series`              |
  | [pdvega](https://altair-viz.github.io/pdvega/)                       | `vgplot`   | `Series`, `DataFrame` |
  | [pandas-genomics](https://pandas-genomics.readthedocs.io/en/latest/) | `genomics` | `Series`, `DataFrame` |
  | [pandas_path](https://github.com/drivendataorg/pandas-path/)         | `path`     | `Index`, `Series`     |
  | [pint-pandas](https://github.com/hgrecco/pint-pandas)                | `pint`     | `Series`, `DataFrame` |
  | [composeml](https://github.com/alteryx/compose)                      | `slice`    | `DataFrame`           |
  | [woodwork](https://github.com/alteryx/woodwork)                      | `slice`    | `Series`, `DataFrame` |
## Development tools

### [pandas-stubs](https://github.com/VirtusLab/pandas-stubs)

While pandas repository is partially typed, the package itself doesn't expose this information for external use.
Install pandas-stubs to enable basic type coverage of pandas API.

Learn more by reading through these issues [14468](https://github.com/pandas-dev/pandas/issues/14468),
[26766](https://github.com/pandas-dev/pandas/issues/26766), [28142](https://github.com/pandas-dev/pandas/issues/28142).

See installation and usage instructions on the [github page](https://github.com/VirtusLab/pandas-stubs).
Title: pandas extension arrays
Date: 2019-01-04

# pandas extension arrays

Extensibility was a major theme in pandas development over the last couple of
releases. This post introduces the pandas extension array interface: the
motivation behind it and how it might affect you as a pandas user. Finally, we
look at how extension arrays may shape the future of pandas.

Extension Arrays are just one of the changes in pandas 0.24.0. See the
[whatsnew][whatsnew] for a full changelog.

## The Motivation

Pandas is built on top of NumPy. You could roughly define a Series as a wrapper
around a NumPy array, and a DataFrame as a collection of Series with a shared
index. That's not entirely correct for several reasons, but I want to focus on
the "wrapper around a NumPy array" part. It'd be more correct to say "wrapper
around an array-like object".

Pandas mostly uses NumPy's builtin data representation; we've restricted it in
places and extended it in others. For example, pandas' early users cared greatly
about timezone-aware datetimes, which NumPy doesn't support. So pandas
internally defined a `DatetimeTZ` dtype (which mimics a NumPy dtype), and
allowed you to use that dtype in `Index`, `Series`, and as a column in a
`DataFrame`. That dtype carried around the tzinfo, but wasn't itself a valid
NumPy dtype.

As another example, consider `Categorical`. This actually composes *two* arrays:
one for the `categories` and one for the `codes`. But it can be stored in a
`DataFrame` like any other column.

Each of these extension types pandas added is useful on its own, but carries a
high maintenance cost. Large sections of the codebase need to be aware of how to
handle a NumPy array or one of these other kinds of special arrays. This made
adding new extension types to pandas very difficult.

Anaconda, Inc. had a client who regularly dealt with datasets with IP addresses.
They wondered if it made sense to add an [IPArray][IPArray] to pandas. In the
end, we didn't think it passed the cost-benefit test for inclusion in pandas
*itself*, but we were interested in defining an interface for third-party
extensions to pandas. Any object implementing this interface would be allowed in
pandas. I was able to write [cyberpandas][cyberpandas] outside of pandas, but it
feels like using any other dtype built into pandas.

## The Current State

As of pandas 0.24.0, all of pandas' internal extension arrays (Categorical,
Datetime with Timezone, Period, Interval, and Sparse) are now built on top of
the ExtensionArray interface. Users shouldn't notice many changes. The main
thing you'll notice is that things are cast to `object` dtype in fewer places,
meaning your code will run faster and your types will be more stable. This
includes storing `Period` and `Interval` data in `Series` (which were previously
cast to object dtype).

Additionally, we'll be able to add *new* extension arrays with relative ease.
For example, 0.24.0 (optionally) solved one of pandas longest-standing pain
points: missing values casting integer-dtype values to float.


```python
>>> int_ser = pd.Series([1, 2], index=[0, 2])
>>> int_ser
0    1
2    2
dtype: int64

>>> int_ser.reindex([0, 1, 2])
0    1.0
1    NaN
2    2.0
dtype: float64
```

With the new [IntegerArray][IntegerArray] and nullable integer dtypes, we can
natively represent integer data with missing values.

```python
>>> int_ser = pd.Series([1, 2], index=[0, 2], dtype=pd.Int64Dtype())
>>> int_ser
0    1
2    2
dtype: Int64

>>> int_ser.reindex([0, 1, 2])
0      1
1    NaN
2      2
dtype: Int64
```

One thing it does slightly change how you should access the raw (unlabeled)
arrays stored inside a Series or Index, which is occasionally useful. Perhaps
the method you're calling only works with NumPy arrays, or perhaps you want to
disable automatic alignment.

In the past, you'd hear things like "Use `.values` to extract the NumPy array
from a Series or DataFrame." If it were a good resource, they'd tell you that's
not *entirely* true, since there are some exceptions. I'd like to delve into
those exceptions.

The fundamental problem with `.values` is that it serves two purposes:

1. Extracting the array backing a Series, Index, or DataFrame
2. Converting the Series, Index, or DataFrame to a NumPy array

As we saw above, the "array" backing a Series or Index might not be a NumPy
array, it may instead be an extension array (from pandas or a third-party
library). For example, consider `Categorical`,

```python
>>> cat = pd.Categorical(['a', 'b', 'a'], categories=['a', 'b', 'c'])
>>> ser = pd.Series(cat)
>>> ser
0    a
1    b
2    a
dtype: category
Categories (3, object): ['a', 'b', 'c']

>>> ser.values
[a, b, a]
Categories (3, object): ['a', 'b', 'c']
```

In this case `.values` is a Categorical, not a NumPy array. For period-dtype
data, `.values` returns a NumPy array of `Period` objects, which is expensive to
create. For timezone-aware data, `.values` converts to UTC and *drops* the
timezone info. These kind of surprises (different types, or expensive or lossy
conversions) stem from trying to shoehorn these extension arrays into a NumPy
array. But the entire point of an extension array is for representing data NumPy
*can't* natively represent.

To solve the `.values` problem, we've split its roles into two dedicated methods:

1. Use `.array` to get a zero-copy reference to the underlying data
2. Use `.to_numpy()` to get a (potentially expensive, lossy) NumPy array of the
   data.

So with our Categorical example,

```python
>>> ser.array
[a, b, a]
Categories (3, object): ['a', 'b', 'c']

>>> ser.to_numpy()
array(['a', 'b', 'a'], dtype=object)
```

To summarize:

- `.array` will *always* be a an ExtensionArray, and is always a zero-copy
   reference back to the data.
- `.to_numpy()` is *always* a NumPy array, so you can reliably call
   ndarray-specific methods on it.

You shouldn't ever need `.values` anymore.

## Possible Future Paths

Extension Arrays open up quite a few exciting opportunities. Currently, pandas
represents string data using Python objects in a NumPy array, which is slow.
Libraries like [Apache Arrow][arrow] provide native support for variable-length
strings, and the [Fletcher][fletcher] library provides pandas extension arrays
for Arrow arrays. It will allow [GeoPandas][geopandas] to store geometry data
more efficiently. Pandas (or third-party libraries) will be able to support
nested data, data with units, geo data, GPU arrays. Keep an eye on the
[pandas ecosystem][eco] page, which will keep track of third-party extension
arrays. It's an exciting time for pandas development.

## Other Thoughts

I'd like to emphasize that this is an *interface*, and not a concrete array
implementation. We are *not* reimplementing NumPy here in pandas. Rather, this
is a way to take any array-like data structure (one or more NumPy arrays, an
Apache Arrow array, a CuPy array) and place it inside a DataFrame. I think
getting pandas out of the array business, and instead thinking about
higher-level tabular data things, is a healthy development for the project.

This works perfectly with NumPy's [`__array_ufunc__`][ufunc] protocol and
[NEP-18][nep18]. You'll be able to use the familiar NumPy API on objects that
aren't backed by NumPy memory.

## Upgrade

These new goodies are all available in the recently released pandas 0.24.

conda:

    conda install -c conda-forge pandas

pip:

    pip install --upgrade pandas

As always, we're happy to hear feedback on the [mailing list][ml],
[@pandas-dev][twitter], or [issue tracker][tracker].

Thanks to the many contributors, maintainers, and [institutional
partners][partners] involved in the pandas community.


[IPArray]: https://github.com/pandas-dev/pandas/issues/18767
[cyberpandas]: https://cyberpandas.readthedocs.io
[IntegerArray]: http://pandas.pydata.org/pandas-docs/version/0.24/reference/api/pandas.arrays.IntegerArray.html
[fletcher]: https://github.com/xhochy/fletcher
[arrow]: https://arrow.apache.org
[ufunc]: https://numpy.org/neps/nep-0013-ufunc-overrides.html
[nep18]: https://www.numpy.org/neps/nep-0018-array-function-protocol.html
[ml]: https://mail.python.org/mailman/listinfo/pandas-dev
[twitter]: https://twitter.com/pandas_dev
[tracker]: https://github.com/pandas-dev/pandas/issues
[partners]: https://github.com/pandas-dev/pandas-governance/blob/master/people.md
[eco]: http://pandas.pydata.org/pandas-docs/stable/ecosystem.html#extension-data-types
[whatsnew]: http://pandas.pydata.org/pandas-docs/version/0.24/whatsnew/v0.24.0.html
[geopandas]: https://github.com/geopandas/geopandas
Title: 2019 pandas user survey
Date: 2019-08-22

<style type="text/css">
table td {
    background: none;
}

table tr.even td {
    background: none;
}

table {
	text-shadow: none;
}

</style>

# 2019 pandas user survey

Pandas recently conducted a user survey to help guide future development.
Thanks to everyone who participated! This post presents the high-level results.

This analysis and the raw data can be found [on GitHub](https://github.com/pandas-dev/pandas-user-surveys) and run on Binder

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pandas-dev/pandas-user-surveys/master?filepath=2019.ipynb)


We had about 1250 repsonses over the 15 days we ran the survey in the summer of 2019.

## About the Respondents

There was a fair amount of representation across pandas experience and frequeny of use, though the majority of respondents are on the more experienced side.



![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_4_0.png)




![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_5_0.png)


We included a few questions that were also asked in the [Python Developers Survey](https://www.jetbrains.com/research/python-developers-survey-2018/) so we could compare Pandas' population to Python's.

90% of our respondents use Python as a primary language (compared with 84% from the PSF survey).





    Yes    90.67%
    No      9.33%
    Name: Is Python your main language?, dtype: object



Windows users are well represented (see [Steve Dower's talk](https://www.youtube.com/watch?v=uoI57uMdDD4) on this topic).





    Linux      61.57%
    Windows    60.21%
    MacOS      42.75%
    Name: What Operating Systems do you use?, dtype: object



For environment isolation, [conda](https://conda.io/en/latest/) was the most popular.




![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_13_0.png)


Most repondents are Python 3 only.





    3        92.39%
    2 & 3     6.80%
    2         0.81%
    Name: Python 2 or 3?, dtype: object



## Pandas APIs

It can be hard for open source projects to know what features are actually used. We asked a few questions to get an idea.

CSV and Excel are (for better or worse) the most popular formats.



![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_18_0.png)


In preperation for a possible refactor of pandas internals, we wanted to get a sense for
how common wide (100s of columns or more) DataFrames are.



![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_20_0.png)


Pandas is slowly growing new exentension types. Categoricals are the most popular,
and the nullable integer type is already almost as popular as datetime with timezone.



![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_22_0.png)


More and better examples seem to be a high-priority development item.
Pandas recently received a NumFOCUS grant to improve our documentation,
which we're using to write tutorial-style documentation, which should help
meet this need.



![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_24_0.png)


We also asked about specific, commonly-requested features.



![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_26_0.png)


Of these, the clear standout is "scaling" to large datasets. A couple observations:

1. Perhaps pandas' documentation should do a better job of promoting libraries that provide scalable dataframes (like [Dask](https://dask.org), [vaex](https://dask.org), and [modin](https://modin.readthedocs.io/en/latest/))
2. Memory efficiency (perhaps from a native string data type, fewer internal copies, etc.) is a valuable goal.

After that, the next-most critical improvement is integer missing values. Those were actually added in [Pandas 0.24](https://pandas.pydata.org/pandas-docs/stable/whatsnew/v0.24.0.html#optional-integer-na-support), but they're not the default, and there's still some incompatibilites with the rest of pandas API.

Pandas is a less conservative library than, say, NumPy. We're approaching 1.0, but on the way we've made many deprecations and some outright API breaking changes. Fortunately, most people are OK with the tradeoff.





    Yes    94.89%
    No      5.11%
    Name: Is Pandas stable enough for you?, dtype: object



There's a perception (which is shared by many of the pandas maintainers) that the pandas API is too large. To measure that, we asked whether users thought that pandas' API was too large, too small, or just right.



![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_31_0.png)


Finally, we asked for an overall satisfaction with the library, from 1 (not very unsatisfied) to 5 (very satisfied).



![png]({{ base_url }}/static/img/blog/2019-user-survey/2019_33_0.png)


Most people are very satisfied. The average response is 4.39. I look forward to tracking this number over time.

If you're analyzing the raw data, be sure to share the results with us [@pandas_dev](https://twitter.com/pandas_dev).
Title: pandas 1.0
Date: 2020-01-29

# pandas 1.0

Today pandas celebrates its 1.0.0 release. In many ways this is just a normal release with a host of new features, performance improvements, and bug fixes, which are documented in our [release notes](https://pandas.pydata.org/pandas-docs/version/1.0.0/whatsnew/v1.0.0.html). But it’s also something a bit more — a milestone for the project beyond just the commits. We wanted to take some time to reflect on where we've been and where we're going.

## Reflections

The world of scientific Python has changed a lot since pandas was started.  In 2011, [the ecosystem was fragmented](https://wesmckinney.com/blog/a-roadmap-for-rich-scientific-data-structures-in-python/): a standard *rich* data structure for statistics and data science had yet to emerge. This echos a similar story for NumPy, which consolidated array efforts that were [previously fragmented](https://numpy.org/old_array_packages.html).

Over the subsequent years, pandas emerged as a *de facto* standard. It’s used by data scientists and analysts and as a data structure for other libraries to build on top of. StackOverflow [cited pandas](https://stackoverflow.blog/2017/09/14/python-growing-quickly/) as one of the reasons for Python being the fastest growing major programming language.

![Growth of pandas](https://149351115.v2.pressablecdn.com/wp-content/uploads/2017/09/related_tags_over_time-1-1000x1000.png)

Today, the ecosystem is in another phase of exploration.
Several new DataFrame implementations are cropping up to fill needs not met by pandas.
We're [working with those projects](https://datapythonista.me/blog/dataframe-summit-at-euroscipy.html) to establish shared standards and semantics for rich data structures.

## Community and Project Health

This release cycle is the first to involve any kind of grant funding for pandas. [Pandas received funding](https://chanzuckerberg.com/eoss/proposals/) as part of the CZI’s [*Essential Open Source Software for Science*](https://medium.com/@cziscience/the-invisible-foundations-of-biomedicine-4ab7f8d4f5dd) [program](https://medium.com/@cziscience/the-invisible-foundations-of-biomedicine-4ab7f8d4f5dd). The pandas project relies overwhelmingly on volunteer contributors. These volunteer contributions are shepherded and augmented by some maintainers who are given time from their employers — our [institutional partners](https://github.com/pandas-dev/pandas-governance/blob/master/people.md#institutional-partners). The largest work item in our grant award was library maintenance, which specifically includes working with community members to address our large backlog of open issues and pull requests.

While a “1.0.0” version might seem arbitrary or anti-climactic (given that pandas as a codebase is nearly 12 years old), we see it as a symbolic milestone celebrating the growth of our core developer team and depth of our contributor base.  Few open source projects are ever truly “done” and pandas is no different. We recognize the essential role that pandas now occupies, and we intend to continue to evolve the project and adapt to the needs of the world’s data wranglers.

## Going Forward

Our [roadmap](https://pandas.pydata.org/pandas-docs/version/1.0.0/development/roadmap.html) contains an up-to-date listing of where we see the project heading over the next few years.
Needless to say, there's still plenty to do.

Check out the [release notes](https://pandas.pydata.org/pandas-docs/version/1.0.0/whatsnew/v1.0.0.html) and visit the [installation page](https://pandas.pydata.org/pandas-docs/version/1.0.0/getting_started/install.html) for instructions on updating to pandas 1.0.
To report a security vulnerability to pandas, please go to https://tidelift.com/security and see the instructions there.
# Contributor Code of Conduct

As contributors and maintainers of this project, and in the interest of
fostering an open and welcoming community, we pledge to respect all people who
contribute through reporting issues, posting feature requests, updating
documentation, submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free
experience for everyone, regardless of level of experience, gender, gender
identity and expression, sexual orientation, disability, personal appearance,
body size, race, ethnicity, age, religion, or nationality.

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery
* Personal attacks
* Trolling or insulting/derogatory comments
* Public or private harassment
* Publishing other's private information, such as physical or electronic
  addresses, without explicit permission
* Other unethical or unprofessional conduct

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

By adopting this Code of Conduct, project maintainers commit themselves to
fairly and consistently applying these principles to every aspect of managing
this project. Project maintainers who do not follow or enforce the Code of
Conduct may be permanently removed from the project team.

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community.

A working group of community members is committed to promptly addressing any
reported issues. The working group is made up of pandas contributors and users.
Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the working group by e-mail (pandas-coc@googlegroups.com).
Messages sent to this e-mail address will not be publicly visible but only to
the working group members. The working group currently includes

- Safia Abdalla
- Tom Augspurger
- Joris Van den Bossche
- Camille Scott
- Nathaniel Smith

All complaints will be reviewed and investigated and will result in a response
that is deemed necessary and appropriate to the circumstances. Maintainers are
obligated to maintain confidentiality with regard to the reporter of an
incident.

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 1.3.0, available at
[https://www.contributor-covenant.org/version/1/3/0/][version],
and the [Swift Code of Conduct][swift].

[homepage]: https://www.contributor-covenant.org
[version]: https://www.contributor-covenant.org/version/1/3/0/
[swift]: https://swift.org/community/#code-of-conduct
# Contributing to pandas

A detailed overview on how to contribute can be found in the **[contributing guide](https://pandas.pydata.org/docs/dev/development/contributing.html)**.
- [ ] closes #xxxx
- [ ] tests added / passed
- [ ] Ensure all linting tests pass, see [here](https://pandas.pydata.org/pandas-docs/dev/development/contributing_codebase.html#pre-commit) for how to run them
- [ ] whatsnew entry
---

name: Feature Request
about: Suggest an idea for pandas
title: "ENH:"
labels: "Enhancement, Needs Triage"

---

#### Is your feature request related to a problem?

[this should provide a description of what the problem is, e.g. "I wish I could use pandas to do [...]"]

#### Describe the solution you'd like

[this should provide a description of the feature request, e.g. "`DataFrame.foo` should get a new parameter `bar` that [...]", try to write a docstring for the desired feature]

#### API breaking implications

[this should provide a description of how this feature will affect the API]

#### Describe alternatives you've considered

[this should provide a description of any alternative solutions or features you've considered]

#### Additional context

[add any other context, code examples, or references to existing implementations about the feature request here]

```python
# Your code here, if applicable

```
