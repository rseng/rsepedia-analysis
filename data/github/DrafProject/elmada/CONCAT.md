<img src="https://github.com/DrafProject/elmada/raw/main/doc/images/elmada_logo.svg" width="450" alt="elmada logo">

---

# elmada: Dynamic electricity carbon emission factors and prices for Europe

**Status:**
[![PyPI](https://img.shields.io/pypi/v/elmada?color=success&label=pypi%20package)](https://pypi.python.org/pypi/elmada)
[![CI](https://github.com/DrafProject/elmada/actions/workflows/CI.yml/badge.svg)](https://github.com/DrafProject/elmada/actions/workflows/CI.yml)
[![CI with conda](https://github.com/DrafProject/elmada/actions/workflows/CI_conda.yml/badge.svg)](https://github.com/DrafProject/elmada/actions/workflows/CI_conda.yml)
[![codecov](https://codecov.io/gh/DrafProject/elmada/branch/main/graph/badge.svg?token=EOKKJG48A9)](https://codecov.io/gh/DrafProject/elmada)

**Usage:**
[![python](https://img.shields.io/badge/python-3.7_|_3.8_|_3.9-blue?logo=python&logoColor=white)](https://github.com/DrafProject/elmada)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![status](https://joss.theoj.org/papers/10.21105/joss.03625/status.svg)][JOSS paper]

**Contribution:**
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1)](https://pycqa.github.io/isort/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Gitter](https://badges.gitter.im/DrafProject/elmada.svg)](https://gitter.im/DrafProject/elmada)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.0-4baaaa.svg)](CODE_OF_CONDUCT.md)

The open-source Python package **elmada** provides electricity carbon emission factors and wholesale prices for European countries.
The target group includes modelers of distributed energy hubs who need **el**ectricity **ma**rket **da**ta (short: **elmada**), e.g., to evaluate the environmental effect of demand response.
**elmada** is part of the [Draf Project] but can be used as a standalone package.

<img src="https://github.com/DrafProject/elmada/raw/main/doc/images/elmada_scheme_scribble.svg" width="650" alt="Elmada scheme scribble">

## Features

* __Dynamic electricity Carbon Emission Factors (CEFs)__ are calculated depending on country and year in up to quarter-hourly resolution.
There are two types of CEFs: __Grid Mix Emission Factors (XEFs)__ and __Marginal Emission Factors (MEFs)__.
While XEFs reflect the carbon footprint of an electricity use (attributional approach), MEFs estimate the carbon impact (consequential approach) of a change in electricity demand (Learn more in the [white paper][CEFWhitepaper] from Tomorrow and WattTime).
Choose between
  * __XEFs__ from fuel type-specific [ENTSO-E] electricity generation data only for Germany (`XEF_EP`),
  * and __XEFs__ & __MEFs__ from merit order based simulations for [30 European Countries][Europe30] (`XEF_PP`, `XEF_PWL`, `MEF_PP`, `MEF_PWL`).
  The according Power Plant method (`PP`) and Piecewise Linear method (`PWL`) are described in the open-access [Applied Energy paper].
  The data used depend on the method chosen, see [scheme below](#cef-scheme).

* __Wholesale electricity prices__ are provided for European countries. You can choose between the real historical [ENTSO-E] data (`hist_EP`) or the simulation results of the `PP` / `PWL` method.

* Other interesting market data such as merit order lists & plots, fuel-specific generation data, or power plant lists are provided as a by-product of the CEF calculations.

## Methodology

With the `XEF_EP` method, XEFs are calculated by multiplying the share matrix *S* (fuel type specific share of electricity generation per time step from [ENTSO-E]) with the intensity vector *ε*  (fuel type specific life cycle carbon emission intensities from [Tranberg.2019]):

<img src="https://render.githubusercontent.com/render/math?math=\mathrm{XEF}^\mathrm{EP}_{t} = S_{t,f}\cdot\varepsilon_f">

The methods `PP`, `PWL`, and `PWLv` are explained in the [Applied Energy paper]. Here is an overview:
 <!-- Converted from pptx via https://convertio.co/ -->
 <img src="https://github.com/DrafProject/elmada/raw/main/doc/images/scheme_CEF_calculation.svg" id='cef-scheme' width="900" alt="scheme_CEF_calculation">

# Data

## Geographic scope

In `elmada`, two-letter country codes ([ISO 3166-1 alpha-2](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2)) are used.

The countries supported by `elmada` can be seen in the map below which is the output of `elmada.plots.cef_country_map(year=2020, method="XEF_EP")`.

<img src="doc/images/cef_country_map.svg" width="600" alt="cef_country_map">

In the [Usage section](#usage) they are referred to as Europe30.
They include:

* 20 countries analyzed in the [Applied Energy paper]: AT, BE, CZ, DE, DK, ES, FI, FR, GB, GR, HU, IE, IT, LT, NL, PL, PT, RO, RS, SI
* 8 countries with only [one reported fossil fuel type][APENsupplPage8]: BA, CH, EE, LV, ME, MK, NO, SE
* 2 countries where installed generation capacity data for 2019 were only available after the publication of the [Applied Energy paper]: BG, SK

## Data modes

You can use **elmada** in two data modes which can be set with `elmada.set_mode(mode=<MODE>)`:

* `mode="safe"` (default):
  * Pre-cached data for 4 years and 20 countries are used. The data are described in the [Applied Energy paper].
  * The years are 2017 to 2020 and the countries AT, BE, CZ, DE, DK, ES, FI, FR, GB, GR, HU, IE, IT, LT, NL, PL, PT, RO, RS, SI.
  * The data is available in the space-saving and quick-to-read [Parquet format] under [.../safe_cache].
* `mode="live"`:
  * Up-to-date data are retrieved on demand and are cached to an OS-specific directory, see `elmada.paths.CACHE_DIR`. A symbolic link to it can be conveniently created by executing `elmada.make_symlink_to_cache()`.
  * Available years are 2017 until the present.
  * Slow due to API requests.
  * Requires valid API keys of ENTSO-E, Morph, Quandl, see [table below](#data-sources).

## Data sources

| Description | Local data location | Source | Channel | Involved in |
|-|-|-|-|-|
| Generation time series & installed generation capacities | [.../safe_cache] or `CACHE_DIR` | [ENTSO-E] | 🔌 on-demand-retrieval via [EntsoePandasClient] (requires valid [ENTSO-E API key]) | CEFs via `EP`, `PP`, `PWL`, `PWLv` |
| Carbon prices (EUA)| [.../safe_cache] or `CACHE_DIR` | [Sandbag] & [ICE] | 🔌 on-demand-retrieval via [Quandl] (requires valid [Quandl API key]) | CEFs via `PP`, `PWL`, `PWLv` |
| Share of CCGT among gas power plants | [.../safe_cache] or `CACHE_DIR` | [GEO] | 🔌 on-demand-download via [Morph] (requires valid [Morph API key])| CEFs via `PWL`, `PWLv` |
| (Average) fossil power plants sizes | [.../safe_cache] or `CACHE_DIR` | [GEO] | 🔌 on-demand-scraping via [BeautifulSoup4] | CEFs via `PWL`, `PWLv` |
| German fossil power plant list with efficiencies | [.../safe_cache] or `CACHE_DIR` | [OPSD] | 🔌 on-demand-download from [here][opsd_download] | CEFs via `PP`, `PWL`, `PWLv` |
| Transmission & distribution losses | [.../worldbank] | [Worldbank] | 💾 manual download from [here][wb] | CEFs via `PP`, `PWL`, `PWLv` |
| Fuel prices for 2015 (+ trends) | [.../from_other.py] (+ [.../destatis]) | [Konstantin.2017] (+ [DESTATIS]) | 🔢 hard-coded values (+ 💾 manual download from [here][destatis_download]) | CEFs via `PP`, `PWL`, `PWLv` |
| Fuel type-specific carbon emission intensities | [.../from_other.py] & [.../tranberg] | [Quaschning] & [Tranberg.2019] | 🔢 hard-coded values | CEFs via `EP`, `PP`, `PWL`, `PWLv` |

## Time zones

The data is in local time since the [Draf Project] focuses on the modeling of individual local energy hubs.
Standard time is used i.e. daylight saving time is ignored.
Also see [this table](https://github.com/DrafProject/marginal-emission-factors/blob/main/README.md#time-zones) of the time zones used.

# Installation

## Using `pip`

```sh
python -m pip install elmada
```

NOTE: Read [here](https://snarky.ca/why-you-should-use-python-m-pip/) why you should use `python -m pip` instead of `pip`.

## From source using conda

For a conda environment including a full editable **elmada** version do the following steps.

Clone the source repository:

```sh
git clone https://github.com/DrafProject/elmada.git
cd elmada
```

Create an conda environment based on `environment.yml` and install an editable local **elmada** version:

```sh
conda env create
```

Activate the environment:

```sh
conda activate elmada
```

## From source without using conda

### For Unix

```sh
git clone https://github.com/DrafProject/elmada.git
cd elmada
python3 -m venv env
source env/bin/activate
python -m pip install -e .[dev]
```

### For Windows

```sh
git clone https://github.com/DrafProject/elmada.git
cd elmada
py -m venv env
.\env\Scripts\activate
py -m pip install -e .[dev]
```

# Tests

This should always work:

```sh
pytest -m="not apikey"
```

This works only if API keys are set as described [below](#optional-set-your-api-keys-and-go-live-mode):

```sh
pytest
```

# Usage

```py
import elmada
```

## OPTIONAL: Set your API keys and go live mode

```py
elmada.set_api_keys(entsoe="YOUR_ENTSOE_KEY", morph="YOUR_MORPH_KEY", quandl="YOUR_QUANDL_KEY")
# NOTE: API keys are stored in an OS-dependent config directory for later use.

elmada.set_mode("live")
```

## Carbon Emission factors

```py
elmada.get_emissions(year=2019, country="DE", method="MEF_PWL", freq="60min", use_datetime=True)
```

... returns marginal emission factors calculated by the `PWL` method with hourly datetime index:

```sh
2019-01-01 00:00:00     990.103492
2019-01-01 01:00:00     959.758367
                          ...
2019-12-31 22:00:00    1064.122146
2019-12-31 23:00:00    1049.852079
Freq: 60T, Name: MEFs, Length: 8760, dtype: float64
```

The `method` argument of `get_emissions()` takes strings that consists of two parts seperated by an underscore.
The first part is the type of emission factor: grid mix emission factors (`XEF`) or marginal emission factors (`MEF`).
The second part determines the calculation method: power plant method (`PP`), piecewise linear method (`PWL`),  or piecewise linear method in validation mode (`PWLv`).

The first part can be omitted (`_PP`, `_PWL`, `_PWLv`) to return a DataFrame that includes additional information.

```py
elmada.get_emissions(year=2019, country="DE", method="_PWL")
```

... returns all output from the PWL method:

```sh
      residual_load  total_load marginal_fuel  efficiency  marginal_cost         MEFs        XEFs
0          21115.00    51609.75       lignite    0.378432      40.889230   990.103492  204.730151
1          18919.50    51154.50       lignite    0.390397      39.636039   959.758367  164.716687
...             ...         ...           ...         ...            ...          ...         ...
8758       27116.00    41652.00       lignite    0.352109      43.946047  1064.122146  388.542911
8759       25437.75    39262.75       lignite    0.356895      43.356723  1049.852079  376.009477
[8760 rows x 7 columns]
```

Additionally, XEFs can be calculated from historic fuel type-specific generation data (`XEF_EP`).

Here is an overview of valid `method` argument values:

| `method` | Return type | Return values | Restriction |
| --: | -- | -- | -- |
| `XEF_PP` | Series | XEFs using PP method | DE |
| `XEF_PWL` | Series | XEFs using PWL method | [Europe30] |
| `XEF_PWLv` | Series | XEFs using PWLv method | DE |
| `MEF_PP` | Series | MEFs from PP method | DE |
| `MEF_PWL` | Series | MEFs using PWL method | [Europe30] |
| `MEF_PWLv` | Series | MEFs using PWLv method | DE |
| `_PP` | Dataframe | extended data for PP method | DE |
| `_PWL` | Dataframe | extended data for PWL method | [Europe30] |
| `_PWLv` | Dataframe | extended data for PWLv method | DE |
| `XEF_EP` | Series | XEFs using fuel type-specific generation data from [ENTSO-E] | [Europe30] |

You can plot the carbon emission factors with

```py
elmada.plots.cefs_scatter(year=2019, country="DE", method="MEF_PP")
```

<img src="https://github.com/DrafProject/elmada/raw/main/doc/images/cefs_scatter.png" width="600" alt="CEFs">

## Wholesale prices

```py
elmada.get_prices(year=2019, country="DE", method="hist_EP")
```

```sh
0       28.32
1       10.07
        ...  
8758    38.88
8759    37.39
Length: 8760, dtype: float64
```

Possible values for the `method` argument of `get_prices()` are:

| `method` | Description | Restriction |
| --: | -- | -- |
| `PP` | Using the power plant method | DE |
| `PWL` | Using piecewise linear method | [Europe30] |
| `PWLv` | Using piecewise linear method in validation mode | DE |
| `hist_EP` | Using historic [ENTSO-E] data | [Europe30] without BA, ME, MK|
| `hist_SM` | Using historic [Smard] data | used only as backup for DE, 2015 and 2018 |

## Merit order

```py
elmada.plots.merit_order(year=2019, country="DE", method="PP")
```

... plots the merit order:

<img src="https://github.com/DrafProject/elmada/raw/main/doc/images/merit_order.svg" width="600" alt="merit_order">

```py
elmada.get_merit_order(year=2019, country="DE", method="PP")
```

... returns the merit order as DataFrame with detailed information on individual power plant blocks.

## Pre-processed data

The following table describes additional `elmada` functions that provide pre-processed data.
Keyword arguments are for example `kw = dict(year=2019, freq="60min", country="DE")`.

| `elmada.` function call | Return type (Dimensions) | Return value | Usage in `elmada` | Used within |
| -- | -- | -- | -- | -- |
| `get_el_national_generation(**kw)` | DataFrame (time, fuel type) | National electricity generation | Share matrix *S* | `XEF_EP` method |
| `get_el_national_generation(**kw).sum(axis=1)` | Series (time) | Total national electricity generation | Proxy for the total load | XEFs calculations |
| `get_residual_load(**kw)` | Series (time) | Conventional national generation | Proxy for the residual load (see [scheme above](#methodology)) | `PP`, `PWL` and `PWLv`|

# Contributing

Contributions in any form are welcome! To contribute changes, please have a look at our [contributing guidelines](CONTRIBUTING.md).

In short:

1. Fork the project and create a feature branch to work on in your fork (`git checkout -b new-feature`).
1. Commit your changes to the feature branch and push the branch to GitHub (`git push origin my-new-feature`).
1. On GitHub, create a new pull request from the feature branch.

# Citing elmada

If you use **elmada** for academic work please cite this paper published in the Journal for Open Source Software:

[![status](https://joss.theoj.org/papers/10.21105/joss.03625/status.svg)][JOSS paper]

```bibtex
@article{Fleschutz2021,
  title = {elmada: Dynamic electricity carbon emission factors and prices for Europe},
  author = {Markus Fleschutz and Michael D. Murphy},
  journal = {Journal of Open Source Software},
  publisher = {The Open Journal},
  year = {2021},
  volume = {6},
  number = {66},
  pages = {3625},
  doi = {10.21105/joss.03625},
}
```

If you use the PP or PWL method, please also cite the open-access [Applied Energy paper]:

[![APEN](https://img.shields.io/badge/AppliedEnergy-10.1016/j.apenergy.2021.117040-brightgreen)][Applied Energy paper]

```bibtex
@article{Fleschutz2021b,
  title = {The effect of price-based demand response on carbon emissions in European electricity markets: The importance of adequate carbon prices},
  author = {Markus Fleschutz and Markus Bohlayer and Marco Braun and Gregor Henze and Michael D. Murphy},
  journal = {Applied Energy},
  year = {2021},
  volume = {295},
  issn = {0306-2619},
  pages = {117040},
  doi = {10.1016/j.apenergy.2021.117040},
}
```

# License

Copyright (c) 2021 Markus Fleschutz

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

<!-- SOURCES -->
[.../destatis]: elmada/data/raw/destatis
[.../from_other.py]: elmada/from_other.py
[.../safe_cache]: elmada/data/safe_cache
[.../tranberg]: elmada/data/raw/tranberg
[.../worldbank]: elmada/data/raw/worldbank
[APENsupplPage8]: https://ars.els-cdn.com/content/image/1-s2.0-S0306261921004992-mmc1.pdf#page=8
[Applied Energy paper]: https://doi.org/10.1016/j.apenergy.2021.117040
[BeautifulSoup4]: https://pypi.org/project/beautifulsoup4
[destatis_download]: https://www.destatis.de/DE/Themen/Wirtschaft/Preise/Publikationen/Energiepreise/energiepreisentwicklung-xlsx-5619001.xlsx?__blob=publicationFile
[DESTATIS]: https://www.destatis.de
[Draf Project]: https://github.com/DrafProject
[ENTSO-E API key]: https://transparency.entsoe.eu/content/static_content/Static%20content/web%20api/Guide.html
[ENTSO-E]: https://transparency.entsoe.eu/
[EntsoePandasClient]: https://github.com/EnergieID/entsoe-py#EntsoePandasClient
[Europe30]: #geographic-scope
[GEO]: http://globalenergyobservatory.org
[ICE]: https://www.theice.com
[JOSS paper]: https://doi.org/10.21105/joss.03625
[Konstantin.2017]: https://doi.org/10.1007/978-3-662-49823-1
[Morph API key]: https://morph.io/documentation/api
[Morph]: https://morph.io
[opsd_download]: https://data.open-power-system-data.org/conventional_power_plants/latest
[OPSD]: https://open-power-system-data.org
[Parquet format]: https://parquet.apache.org
[Quandl API key]: https://docs.quandl.com/docs#section-authentication
[Quandl]: https://www.quandl.com
[Quaschning]: https://www.volker-quaschning.de/datserv/CO2-spez/index_e.ph
[Sandbag]: https://sandbag.org.uk/carbon-price-viewer
[Smard]: https://www.smard.de/en
[Tranberg.2019]: https://doi.org/10.1016/j.esr.2019.100367
[wb]: https://databank.worldbank.org/reports.aspx?source=2&series=EG.ELC.LOSS.ZS
[CEFWhitepaper]: https://docplayer.net/217796110-A-vision-for-how-ambitious-organizations-can-accurately-measure-electricity-emissions-to-take-genuine-action.html
[Worldbank]: https://databank.worldbank.org/reports.aspx?source=2&series=EG.ELC.LOSS.ZS
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, caste, color, religion, or sexual identity
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
[INSERT CONTACT METHOD].
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

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
[https://www.contributor-covenant.org/version/2/0/code_of_conduct.html][v2.0].

Community Impact Guidelines were inspired by
[Mozilla's code of conduct enforcement ladder][Mozilla CoC].

For answers to common questions about this code of conduct, see the FAQ at
[https://www.contributor-covenant.org/faq][FAQ]. Translations are available
at [https://www.contributor-covenant.org/translations][translations].

[homepage]: https://www.contributor-covenant.org
[v2.0]: https://www.contributor-covenant.org/version/2/0/code_of_conduct.html
[Mozilla CoC]: https://github.com/mozilla/diversity
[FAQ]: https://www.contributor-covenant.org/faq
[translations]: https://www.contributor-covenant.org/translations
# How to contribute

To guarantee the further evolution of elmada in the long term, we depend on the support of volunteer developers.

Some of the resources to look at if you're interested in contributing:
* [Join us on Gitter to chat!](https://gitter.im/DrafProject/elmada)

## Licensing

By contributing to elmada, e.g. through opening a pull request or submitting a patch, you represent that your contributions are your own original work and that you have the right to license them, and you agree that your contributions are licensed under the LGPL 3 license.

## Submitting bug reports

[Open an issue on GitHub](https://github.com/DrafProject/elmada/issues/new) to report bugs or other problems.

## Submitting changes

To contribute changes:

1. Fork the project on GitHub
1. Create a feature branch to work on in your fork (``git checkout -b new-fix-or-feature``)
1. Commit your changes to the feature branch after running black to format your code
1. Push the branch to GitHub (``git push origin new-fix-or-feature``)
1. On GitHub, create a new [pull request](https://github.com/DrafProject/elmada/pull/new/master) from the feature branch

### Pull requests

Before submitting a pull request, check whether you have:

* Added or updated documentation for your changes
* Added tests if you implemented new functionality

When opening a pull request, please provide a clear summary of your changes!

### Commit messages

Please try to write clear commit messages. One-line messages are fine for small changes, but bigger changes should look like this:

    A brief summary of the commit

    A paragraph or bullet-point list describing what changed and its impact,
    covering as many lines as needed.

## Testing

We have existing test coverage for the key functionality of elmada.

All tests are in the ``elmada/tests`` directory and use [pytest](https://docs.pytest.org/en/latest/).

Our test coverage is not perfect. An easy way to contribute code is to work on better tests.

## Coding conventions

Start reading our code and you'll get the hang of it.

We mostly follow the official [Style Guide for Python Code (PEP8)](https://www.python.org/dev/peps/pep-0008/).

We have chosen to use the uncompromising code formatter, [`black`](https://github.com/psf/black/).
If run from the root directory of this repo, `pyproject.toml` should ensure the line lengths are restricted to 100.
The philosophy behind using black is to have uniform style throughout the project dictated by code.
Since `black` is designed to minimise diffs, and make patches more human readable, this also makes code reviews more efficient.

## Attribution

The layout and content of this document is based on the contribution guidelines of the projects [OpenGovernment](https://github.com/opengovernment/opengovernment/blob/master/CONTRIBUTING.md) and [Calliope](https://github.com/calliope-project/calliope/blob/master/CONTRIBUTING.md).
---
title: 'elmada: Dynamic electricity carbon emission factors and prices for Europe'
tags:
  - energy
  - electricity market
  - carbon emissions
  - marginal emissions
  - python
authors:
 - name: Markus Fleschutz
   orcid: 0000-0002-8516-9635
   affiliation: 1, 2
 - name: Michael D. Murphy
   orcid: 0000-0002-4269-2581
   affiliation: 1
affiliations:
 - name: Department of Process, Energy and Transport Engineering, Munster Technological University
   index: 1
 - name: Institute of Refrigeration, Air-Conditioning, and Environmental Engineering, Karlsruhe University of Applied Sciences
   index: 2
date: 23 July 2021
bibliography: paper.bib
---

# Summary

The expansion of intermittent renewable energy sources such as solar and wind requires increased operational flexibility in electricity systems.
Energy system models at the scale of individual decentral energy hubs can help decision-makers of energy hubs such as city quarters or industrial sites evaluate the cost and carbon emission saving potentials of their flexibility.
For national scale models, the carbon emissions of the electricity supply system are endogenously determined.
However, low-level models (at the scale of decentral energy hubs) need this information as input.
Since specific carbon emissions of national electricity supply systems fluctuate hourly, the usage of dynamic (i.e. at least hourly resolved) carbon emission factors (CEFs) is essential [@Prina2020].

`elmada` is an easy-to-use open-source Python package designed to provide dynamic electricity CEFs and prices for European countries.
The target group includes modelers of distributed energy hubs who need electricity market data.
This is where the name **elmada** comes from: **el**ectricity **ma**rket **da**ta.
`elmada` is developed in the open on GitHub [@ElmadaGitHub].
Each release is archived on Zenodo [@ElmadaZenodo].

# Statement of Need

Dynamic CEFs are important for the environmental assessment of electricity supply in not fully decarbonized energy systems.
To the best of the authors' knowledge, `elmada` is the first free and open-source Python interface for dynamic CEFs in Europe.
This makes `elmada` an important complement to existing commercial services.

At the moment, there are two main commercial services that provide an Application Programming Interface (API) for historical dynamic CEFs: the `electricityMap` API [@ElecMapApi] and the Automated Emissions Reduction from WattTime [@WattTime].
The `electricityMap` is maintained by Tomorrow, a startup based in Denmark, and WattTime is a nonprofit organization in the USA.
However, both focus on real-time CEFs as incentive signals for demand response answering the question "How clean is my electricity right now?".
We elaborate more on `electricityMap` here, as they originate in Europe, which is also the focus of `elmada`.
The services of WattTime are broadly similar.

`electricityMap` is a software project that visualizes the carbon emission intensity linked to the generation and consumption of electricity on a global choropleth map.
Additionally, the `electricityMap` API provides historical, real-time (current hour), forecast, and since recently also marginal data. The calculation methods consider international energy exchanges and the fact that the list of data sources is curated by Tomorrow (the company behind it) makes it save-to-use as a live incentive signal e.g. for carbon-based demand response applications.
However, the use of `electricityMap` API requires a data-dependent payment even for the historic data, so it is not free of charge.

There are two types of dynamic CEFs:

* grid-mix emission factors (XEFs), which represent the emission intensity based on the current generation mix of the electricity system,
* and marginal emission factors (MEFs), which quantify the emission intensity of the generators likely to react to a marginal system change.

Currently, there is no multi-national solution for modelers of decentral energy hubs searching for free historical hourly CEFs (in particular MEFs).
This gap often leads to the usage of yearly average CEFs, which are potentially misleading [@Hawkes2010].
We close this gap by providing the conveniently installable Python package `elmada` that calculates XEFs and MEFs in hourly (or higher) resolution for 30 European countries for free.
`elmada` provides modelers a no-regret (free) entry point to European dynamic CEFs leaving it open to the modeler to later switch to a paid service that generate more accurate CEFs, e.g. through advanced methods such as the consideration of cross-border energy flows through flow tracing [@Tranberg2019].

# Functionality

`elmada` calculates both types of dynamic CEFs: XEFs and MEFs.
MEFs are more challenging to approximate than XEFs since MEFs require the identification of the marginal power plants per time step.
In `elmada`, this is done through a merit order simulation within the power plant (PP) and piecewise linear (PWL) method described in [@Fleschutz2021].

Also, historical and simulated day-ahead electricity market prices are provided.
They can be used either for the economic evaluation of electricity demands or to model the incentive signal of price-based demand response.

Currently, `elmada` provides data for 30 European countries and for each year since 2017.
`elmada` works mainly with data from the ENTSO-E Transparency Platform [@ENTSOE].

# Current and Future Usage

So far, `elmada` has been used in a study where MEFs, XEFs and the results of load shift simulations based on them are compared across 20 European countries [@Fleschutz2021].
In ongoing research, `elmada` is used to quantify the costs and emission-saving potentials that arise from the exploitation of existing and future flexibility in decentral energy hubs.

We hope that `elmada` reduces the difficulty associated with the use of dynamic CEFs and prices in the modeling of decentral energy systems.

# Acknowledgements

The author acknowledges funding by the MTU Risam scholarship scheme and the German Federal Ministry of the Environment, Nature Conservation and Nuclear Safety (BMU) via the project WIN4Climate (No. 03KF0094 A) as part of the National Climate Initiative.

# References
