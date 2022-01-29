<!--
SPDX-FileCopyrightText: 2020 Wolfgang Traylor <wolfgang.traylor@senckenberg.de>

SPDX-License-Identifier: CC-BY-4.0
-->

# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog][] by Olivier Lacan, and this project adheres to [Semantic Versioning][].

[Keep a Changelog]: <https://keepachangelog.com/en/1.0.0/>
[Semantic Versioning]: <https://semver.org/spec/v2.0.0.html>

## [1.1.3] - 2021-12-14

### Added
- Option to reset date in `Fauna::World` - e.g. to start simulation from the beginning in the next location/gridcell

### Changed
- Move quickstart guide from Doxygen docs to `README.md` [#52]
- Update PlantUML from 1.2019.09 to 1.2021.15

### Fixed
- Restore library interface for backward compatibility: `Fauna::World::World()` and `Fauna::World::is_activated()`
- Exception if started without HFTs
- Now HFT table files won’t be created if there are no HFTs defined.
- Too low digestibility values in output, due to wrong weights in aggregation.
- CI doxygen generation failing with error "sh: 1: /usr/sbin/dot: not found". `dot` is now disabled.

## [1.1.2] - 2021-07-27

### Added
- `demo_results.Rmd` now only tries to plot those files that are available [#39]

### Fixed
- Allow `simulation.one_hft_per_habitat` only for `herbivore_type = "Cohort"`
- Issue that only one HFT got created and outputted if `simulation.one_hft_per_habitat = true`

## [1.1.1] - 2021-07-26

### Fixed
- Implement `simulation.one_hft_per_habitat` option

## [1.1.0] - 2021-06-15

### Added
- Code coverage badge with [codecov.io](https://codecov.io/gh/wtraylor/modular_megafauna_model/)
- Software metadata in [codemeta format](https://codemeta.github.io/)
- Add parameter sanity checks: mortality must not exceed reproduction; expenditure must not exceed intake [#45]
- New text table output: `individual_density` and `body_fat` [#43]
- The parameter `hft.breeding_season.start` can now take a month name (alternative to the Julian day).

### Changed
- Rename output file `mass_density_per_hft` to `mass_density` [#44]

### Fixed
- When configuring CMake with `BUILD_DEMO_SIMULATOR=ON` but `BUILD_TESTING=OFF`, the `megafauna.toml` instruction file was not copied even though `run_demo` needs it.
- Valgrind memory check in CI didn’t fail
- Doxygen issues with ReadTheDocs [#34] [#41] [#42]

## [1.0.3] - 2021-05-10

### Fixed
- Anabolism efficiency was forced to be lower than catabolism efficiency, but for no apparent reason.

## [1.0.2] - 2021-05-05

### Changed
- Increase demo simulation years so we can see some population ups and downs.

### Fixed
- The unit tests were broken after the last release.
- cpptoml wouldn’t compile with latest `g++` (11.1.0).

## [1.0.1] - 2021-04-29

### Fixed
- `DigestiveLimit::FixedFraction` scaled like `DigestiveLimit::Allometric`. Fixed fraction is now fixed and not scaling. As a result, the newborns of the example herbivore (`example.toml`) now die if their intake is set to a fixed fraction. Therefore the example uses the allometric limit.
- The exponent of `Hft::digestion_allometric` was wrong. It was set to the exponent of `Hft::expenditure_basal_rate`, but it must be 1 minus that because the digestion-limited intake is *relative* to body mass.
- If `DigestiveLimit::Allometric` was used, the calculated intake rates were far too low because the fraction was not multiplied with body mass.

### Added
- Check code format in continuous integration.

## [1.0.0] - 2021-04-26

### Added
- R scripts and LibreOffice document to reproduce figures in `docs/images/`. [#31]
- Minor additions to docs.
- Integration for [readthedocs.org](https://readthedocs.org)
- CI check for Git tag matching version in CMakeLists.txt

### Changed
- Use sward density for `Fauna::HalfMaxIntake` (instead of whole-habitat grass density). This requires the vegetation model to provide correct “FPC” values (fractional area covered by grass in the habitat). [#12]
- Use modern `rmarkdown` package to render demo results.

### Fixed
- Demo simulator has now monthly ambient air temperature.
- Replace copyrighted figure `docs/images/thermoregulation.png` with my own. [#32]

### Removed
- Snow depth. It is not used anywhere. [#14]

## [0.6.0] - 2021-04-08

### Added
- Invitation to Matrix room in README.

### Fixed
- Add licensing information. [#22]

### Removed
- Nitrogen retention by herbivores. The function `Habitat::add_excreted_nitrogen()` has been removed. The vegetation model must now handle recycling of nitrogen itself. This simplifies the code and reduces the number of parameters. Compare [#15]

## [0.5.5] - 2020-08-04

### Fixed
- Floating point imprecision. [#29]
- Too much forage got eaten in one day in the demo simulator. [#29]

## [0.5.4] - 2020-08-03

### Fixed
- Misleading error message if TOML parameter has wrong type. [#25]
- Group parameter `expenditure.components` cannot be re-used. [#28]
- No error if an HFT group is defined twice. [#27]
- No instruction file check that catabolism must be smaller than anabolism coefficient. [#26]

## [0.5.3] - 2020-06-05

### Fixed
- Crash when reading instruction file without parameter `hft.digestion.i_g_1992_ijk`.

### Changed
- Single values in TOML file can be parsed like a single-element array.

## [0.5.2] - 2020-06-02

### Fixed
- Crash if not all HFTs were established immediately at the start of the start of the simulation.

## [0.5.1] - 2020-05-29

### Added
- Executable `megafauna_insfile_linter` for checking if a TOML instruction file is valid.

### Fixed
- Herbivores are now only established if the `do_herbivores==true` flag is passed into MMM.
- Parameter `hft.digestion.i_g_1992_ijk` is not mandatory anymore.
- Fix compiling errors under GCC 10.1.0 about not finding `std` exception.

## [0.5.0] - 2020-04-10

### Added
- Revised and expanded Doxygen page on model description/discussion.
- Parameter `hft.expenditure.fmr_multiplier`. [#9]

### Changed
- Output `eaten_nitrogen_per_ind` is now in milligram, not kilogram, per day and individual.
- Remove `Fauna::Parameters` from the library’s include interface.
- The Doxygen documentation of the library’s include interface can be parsed without errors.
- Reproductive success at parturition is now based on body condition at day of conception. [#10]
- Scale maximum daily intake (DMI) allometrically from fraction of adult mass. [#8]
    - Parameter `hft.digestion.allometric.coefficient` changed to `.value_male_adult`.
- Expenditure component `Allometric` is now `BasalMetabolicRate` and `FieldMetabolicRate`. [#9]
    - Parameter `hft.expenditure.allometric` changed to `.basal_rate`.

### Removed
- Individual mode. [#7]

## [0.4.0] - 2020-03-27

### Added
- New output tables: `available_forage`, `eaten_forage_per_ind`, `eaten_nitrogen_per_ind`.

## [0.3.1] - 2020-03-26

### Fixed
- `Fauna::World` can be constructed without an instruction file. [#4]
- Correctly plot demo results with `demo_results.Rmd` if output is daily or decadal.
- Fix “First day before last” error in decadal output. [#3]

## [0.3.0] - 2020-03-17
### Added
- Net energy content model: `NetEnergyModel::GrossEnergyFraction`
    - Parameter `forage.gross_energy`
    - Parameter `hft.digestion.digestibility_multiplier`
    - Parameter `hft.digestion.k_fat`
    - Parameter `hft.digestion.k_maintenance`
    - Parameter `hft.digestion.me_coefficient`
- Model description of forage energy and digestion.
- Parameter `hft.body_mass.empty`.
- Parameter `hft.body_fat.catabolism_efficiency`.

### Changed
- Fractional body fat now refers to the empty body (i.e. without ingesta, blood, etc.)  and not the live body mass.
- HFT is made optional. New implementations of `HerbivoreInterface` don’t require an HFT. [#5]

### Fixed
- Order of HFT names in text table output is now guaranteed.
- Compiler flags specific to GCC are removed. [#1]
- Unknown TOML parameters/keys now issue an error. [#2]
- Parse TOML parameter `hft.body_fat.gross_energy`.
- Run directory is now created in script `tools/run_valgrind_memcheck`.

### Removed
- The old “default” net energy model `get_net_energy_content_default()`.
- Parameter `hft.digestion.efficiency` (less efficient digestion of hindgut
  fermenters is now in `hft.digestion.digestibility_multiplier`).

## [0.2.0] - 2019-12-06
### Added
- New instruction file parameters, which were constants before:
    - `hft.digestion.anabolism_coefficient` (formerly `Fauna::FatMassEnergyBudget::FACTOR_ANABOLISM`)
    - `hft.digestion.catabolism_coefficient` (formerly `Fauna::FatMassEnergyBudget::FACTOR_CATABOLISM`)
    - `hft.digestion.efficiency` (formerly `Fauna::DIGESTION_EFFICIENCY_HINDGUTS`)
    - `hft.digestion.i_g_1992_ijk` (formerly constants in the function object `Fauna::GetDigestiveLimitIlliusGordon1992`)
    - `simulation.metabolizable_energy` (formerly `Fauna::ME_COEFFICIENT_GRASS`)
    - `hft.reproduction.logistic.growth_rate` and `.midpoint` (formerly constant parameters in `Fauna::ReprIlliusOConnor2000`).
- Default HFT groups “ruminants” and “hindguts” (replacing parameter `hft.digestion.type`)
- Simple bash script `run_demo` produced in the build folder to execute demo simulator with results in one command.

### Changed
- Replaced `std::map` with `std::array` in `Fauna::ForageValues` to improve speed.
- Various little performance improvements.
- Turned `Fauna::GetNetEnergyContent` interface (strategy design pattern) to a function.
- Turned `Fauna::GetDigestiveLimitIlliusGordon1992` into the function `Fauna::get_digestive_limit_illius_gordon_1992()`.
- Renamed TOML parameter `hft.foraging.net_energy_model` to `hft.digestion.net_energy_model`
- Renamed `Fauna::ReproductionModel::IlliusOConnor2000` to `Logistic`.

### Fixed
- The parameter `hft.mortality.factors` was not parsed.

### Removed
- Instruction file parameter `hft.digestion.type` and the corresponding `Fauna::DigestionType` and `Fauna::Hft::digestion_type`.

## [0.1.0] - 2019-11-07
### Added
- Herbivores in cohort and individual mode.
    - Energy expenditure components:
        - Allometric
        - Based on [Taylor et al. (1981)][]
        - Thermoregulation
    - Diet composer: only grass
    - Reproduction models:
        - Constant annual reproduction rate
        - Based on [Illius & O’Connor (2000)][]
        - Linear relationship with body condition
    - Mortality factors:
        - Constant annual background mortality
        - Death at end of lifespan
        - Starvation mortality based on [Illius & O’Connor (2000)][]
        - Starvation at a threshold value of body condition
    - Foraging limits:
        - Functional response based on [Illius & O’Connor (2000)][]
        - General Holling Type II functional response
    - Net energy in forage: formula used by [Illius & O’Connor (2000)][]
- Continuous integration (CI) for GitLab.
- Output in tab-separated text tables.
- TOML instruction file reader.
- Demo simulator with simple logistic grass growth.

[Illius & O’Connor (2000)]: <https://doi.org/10.2307/3800911>
[Taylor et al. (1981)]: <https://doi.org/10.1017/S0003356100040617>

[Unreleased]: https://github.com/wtraylor/modular_megafauna_model/compare/1.1.3...develop
[1.1.3]: https://github.com/wtraylor/modular_megafauna_model/compare/1.1.2...1.1.3
[1.1.2]: https://github.com/wtraylor/modular_megafauna_model/compare/1.1.1...1.1.2
[1.1.1]: https://github.com/wtraylor/modular_megafauna_model/compare/1.1.0...1.1.1
[1.1.0]: https://github.com/wtraylor/modular_megafauna_model/compare/1.0.3...1.1.0
[1.0.3]: https://github.com/wtraylor/modular_megafauna_model/compare/1.0.2...1.0.3
[1.0.2]: https://github.com/wtraylor/modular_megafauna_model/compare/1.0.1...1.0.2
[1.0.1]: https://github.com/wtraylor/modular_megafauna_model/compare/1.0.0...1.0.1
[1.0.0]: https://github.com/wtraylor/modular_megafauna_model/compare/0.6.0...1.0.0
[0.6.0]: https://github.com/wtraylor/modular_megafauna_model/compare/0.5.5...0.6.0
[0.5.5]: https://github.com/wtraylor/modular_megafauna_model/compare/0.5.4...0.5.5
[0.5.4]: https://github.com/wtraylor/modular_megafauna_model/compare/0.5.3...0.5.4
[0.5.3]: https://github.com/wtraylor/modular_megafauna_model/compare/0.5.2...0.5.3
[0.5.2]: https://github.com/wtraylor/modular_megafauna_model/compare/0.5.1...0.5.2
[0.5.1]: https://github.com/wtraylor/modular_megafauna_model/compare/0.5.0...0.5.1
[0.5.0]: https://github.com/wtraylor/modular_megafauna_model/compare/0.4.0...0.5.0
[0.4.0]: https://github.com/wtraylor/modular_megafauna_model/compare/0.3.1...0.4.0
[0.3.1]: https://github.com/wtraylor/modular_megafauna_model/compare/0.3.0...0.3.1
[0.3.0]: https://github.com/wtraylor/modular_megafauna_model/compare/0.2.0...0.3.0
[0.2.0]: https://github.com/wtraylor/modular_megafauna_model/compare/0.1.0...0.2.0
[0.1.0]: https://github.com/wtraylor/modular_megafauna_model/releases/tag/0.1.0

[#1]: https://github.com/wtraylor/modular_megafauna_model/issues/1
[#2]: https://github.com/wtraylor/modular_megafauna_model/issues/2
[#3]: https://github.com/wtraylor/modular_megafauna_model/issues/3
[#4]: https://github.com/wtraylor/modular_megafauna_model/issues/4
[#5]: https://github.com/wtraylor/modular_megafauna_model/issues/5
[#6]: https://github.com/wtraylor/modular_megafauna_model/issues/6
[#7]: https://github.com/wtraylor/modular_megafauna_model/issues/7
[#8]: https://github.com/wtraylor/modular_megafauna_model/issues/8
[#9]: https://github.com/wtraylor/modular_megafauna_model/issues/9
[#10]: https://github.com/wtraylor/modular_megafauna_model/issues/10
[#12]: https://github.com/wtraylor/modular_megafauna_model/issues/12
[#14]: https://github.com/wtraylor/modular_megafauna_model/issues/14
[#15]: https://github.com/wtraylor/modular_megafauna_model/issues/15
[#22]: https://github.com/wtraylor/modular_megafauna_model/issues/22
[#25]: https://github.com/wtraylor/modular_megafauna_model/issues/25
[#26]: https://github.com/wtraylor/modular_megafauna_model/issues/26
[#28]: https://github.com/wtraylor/modular_megafauna_model/issues/28
[#29]: https://github.com/wtraylor/modular_megafauna_model/issues/29
[#31]: https://github.com/wtraylor/modular_megafauna_model/issues/31
[#32]: https://github.com/wtraylor/modular_megafauna_model/issues/32
[#34]: https://github.com/wtraylor/modular_megafauna_model/issues/34
[#39]: https://github.com/wtraylor/modular_megafauna_model/issues/39
[#41]: https://github.com/wtraylor/modular_megafauna_model/issues/41
[#43]: https://github.com/wtraylor/modular_megafauna_model/issues/43
[#44]: https://github.com/wtraylor/modular_megafauna_model/issues/44
[#45]: https://github.com/wtraylor/modular_megafauna_model/issues/45
[#52]: https://github.com/wtraylor/modular_megafauna_model/issues/52
<!--
SPDX-FileCopyrightText: 2020 Wolfgang Traylor <wolfgang.traylor@senckenberg.de>

SPDX-License-Identifier: CC-BY-4.0
-->

Modular Megafauna Model
=======================

[![JOSS](https://joss.theoj.org/papers/10.21105/joss.03631/status.svg)](https://doi.org/10.21105/joss.03631)
[![DOI](docs/images/zenodo_doi.svg)](https://zenodo.org/badge/latestdoi/228426088)
[![REUSE-compliant](docs/images/reuse-compliant.svg)][REUSE]
[![Documentation Status](https://readthedocs.org/projects/modular-megafauna-model/badge/?version=latest)](https://modular-megafauna-model.readthedocs.io/en/latest/?badge=latest)
[![Pipeline Status](https://gitlab.com/wtraylor/modular_megafauna_model/badges/master/pipeline.svg)](https://gitlab.com/wtraylor/modular_megafauna_model/-/commits/master)
[![codecov](https://codecov.io/gh/wtraylor/modular_megafauna_model/branch/master/graph/badge.svg)](https://codecov.io/gh/wtraylor/modular_megafauna_model)

[![LGPL logo](docs/images/lgpl.svg)](https://choosealicense.com/licenses/lgpl-3.0/)

[REUSE]: https://reuse.software



Overview
--------

The Modular Megafauna Model (MMM) is a dynamic, mechanistic, and process-based abstraction of large herbivore populations through space and time.
Herbivores feed, grow, reproduce, and die on a daily simulation cycle.
The amount of forage as well as relevant environment variables need to be supplied by the calling program, for instance by a dynamic global vegetation model (DGVM).

The megafauna model itself is spatially inexplicit, but operates in different singular spatial units, which can be mapped to spatially meaningful cells by the calling framework.
Likewise, the megafauna model is ignorant of calendar dates, but is only aware of the Julian day of the year.

As the name suggests, modularity is a primary design goal.
Processes and model components can be switched on and off, replaced, and expanded.
Since this C++ library is Free Software, the scientific community is encouraged to use, study, change, and redistribute the source code.

Come join the MMM user and developer room on Matrix: <https://matrix.to/#/!rnevkLtJTORmvyzFHD:matrix.org?via=matrix.org>

Project documentation: <https://modular-megafauna-model.readthedocs.io/>

Repository Structure
--------------------

This project follows the [Pitchfork Layout](https://github.com/vector-of-bool/pitchfork).

- `docs/`: Doxygen documentation.
- `examples/`: Exemplary instruction files.
- `external/`: Embedded external projects, which remain unmodified.
- `include/`: Public API header files.
- `LICENSES/`: Folder with licenses, compliant with [REUSE][].
- `src/`: Source and (private) header files of the project. Subdirectories correspond to C++ namespaces.
- `tests/`: Unit tests and test scripts.
- `tools/`: Different helper tools for the developer.

Quickstart
----------

### Compile and Run Unit Tests

To check if the modular megafauna model works by itself correctly, you should compile it and run the unit tests.
You will need [Cmake](https://cmake.org) (3.10 or later) and a C++ compiler supporting the C++11 standard.
It should work well with GCC, the [GNU C++ compiler](https://gcc.gnu.org) (4.8 or later).

Open a Unix shell in the root of this repository and run:

```sh
mkdir -p build
cd build
cmake -DBUILD_TESTING=ON ..
make megafauna_unit_tests
./megafauna_unit_tests
```

If the unit tests all pass successfully, you are good to go!

### Run the Demo Simulator

With the CMake option `BUILD_DEMO_SIMULATOR=ON` you can compile the demo simulator.
This independent program is a very simple grass simulator that hosts the megafauna library.
It should work out of the box with the instructions files in the `examples/` folder.

Again, open a terminal in the repository root folder and run:

```sh
mkdir -p build
cd build
cmake -DBUILD_DEMO_SIMULATOR=ON ..
make megafauna_demo_simulator
./megafauna_demo_simulator "../examples/megafauna.toml" "../examples/demo_simulation.toml"
```

Congratulations, you have run your first simulation!

The output files are tab-separated text (`.tsv`) files.
The simulator does not want to have you lose your results and will refuse to overwrite existing output files.
So if you want to run it again, you will first need to remove the previously created output tables: `rm *.tsv`

A very simple [RMarkdown](https://rmarkdown.rstudio.com/) file to visualize the demo output is provided in the `build` folder.
It is called `demo_results.Rmd`.
You can open and render it in [RStudio](https://www.rstudio.com/).

Alternatively, you can compile the RMarkdown file into an HTML file without RStudio.
You will need to have [R](https://www.r-project.org/) (version 3, 4, or later) and [Pandoc](https://pandoc.org) installed.

Open a terminal in the `build` directory and execute `R` to get into an interactive R session.
In the R console, you execute the following commands:

```r
install.packages("rmarkdown")
library(rmarkdown)
render("demo_results.Rmd")
quit(save = "no")
```

This should have produced the HTML file `demo_results.html`, which you can open in your web browser.
The results should show a rising population curve, which means that the example herbivore has been able to survive and reproduce.
You can compare it to the expected output (“artifacts”) produced in the last Continuous Integration (CI) run on [GitLab][gitlab_ci_result] (download “demo_simulation:archive” through the three-dots menu).

[gitlab_ci_result]: <https://gitlab.com/wtraylor/modular_megafauna_model/-/pipelines?page=1&scope=all&ref=master&status=success>

Feel free to play with the MMM parameters in `megafauna.toml` and the Demo Simulator parameters in `demo_simulation.toml`.

### Customize the Instruction File

The Modular Megafauna Model library requires a single instruction file to set general simulation parameters and herbivore (HFT = Herbivore Functional Type) parameters.
The instruction file is in [TOML v0.5][] format; see there for a description of the syntax.
String parameters are generally case insensitive.

[TOML v0.5]: https://github.com/toml-lang/toml/tree/v0.5.0

All possible options are listed in the example file under [examples/megafauna.toml](examples/megafauna.toml).
You will find all parameters explained in detail in the Doxygen documentation of the two classes [Fauna::Parameters][] and [Fauna::Hft][].
The parameter names and cascaded categories are kept as consistent as possible between the C++ classes and the TOML file.

[Fauna::Parameters]: <https://modular-megafauna-model.readthedocs.io/en/latest/structFauna_1_1Parameters.html>
[Fauna::Hft]: <https://modular-megafauna-model.readthedocs.io/en/latest/classFauna_1_1Hft.html>

Note that the example HFT in [examples/megafauna.toml](examples/megafauna.toml) is intentionally fictional.
The parameter values do fall in a range realistic for an ungulate grazer, but for your study you will want to parametrize your organism of interest based on physiological data from the literature.

Both HFTs and HFT groups are represented as arrays of tables in the TOML syntax.
You can define any number of HFTs, but they need to have unique names.
Any HFT can be assigned to any number of groups to inherit parameters from that group.
However, you cannot cascade groups, and an HFT cannot inherit from another HFT.

In general, the user is forced to specify all parameters.
This way one instruction file is always complete and self-explanatory.
Default values in the C++ code might change between model versions.
However, the model should still yield the same results with the same instruction file.

Compiling the Documentation
---------------------------

You can find the automatically compiled documentation of the latest release here:
<https://modular-megafauna-model.readthedocs.io/>

As a bare minimum, you will need to have [CMake](https://cmake.org) (version 3.10 or higher) and [Doxygen](https://www.doxygen.nl) (version 1.8.16 or higher) installed.

Open a Unix shell (terminal) in the root directory of the megafauna library, and execute the following lines.
On Windows, you can use the [Windows Subsystem for Linux]() <!--TODO-->(Windows 10 or higher) or try to compile it with the CMake GUI and/or an IDE.

```bash
mkdir build
cd build
cmake -DBUILD_DOC=ON ..
make megafauna_docs
```

Don’t worry if warning messages appear. Usually, most of the documentation
will be fine.
Now open the created file `docs/index.html` in a web browser.

### Optional Build Requirements for the Documentation
- Java Runtime Environment (JRE) and [Graphviz](www.graphviz.org) to compile [PlantUML](http://plantuml.com) diagrams. See [here](http://plantuml.com/graphviz-dot) for details.
- [LaTeX](www.latex-project.org) to render mathematical formulas offline, and [BibTeX](www.bibtex.org) for the bibliography.

Existing Integrations
---------------------

Originally this megafauna model was developed for the dynamic global vegetation model [LPJ-GUESS](http://iis4.nateko.lu.se/lpj-guess/).
On the Lund subversion server there exists a branch `megafauna` that integrates this library into LPJ-GUESS.
The educational version of LPJ-GUESS does *not* integrate MMM.
LPJ-GUESS is proprietary and closed-source.
Please contact the maintainers of LPJ-GUESS to kindly ask for access.

Other dynamic vegetation models can include the megafauna model as an external library, too.
Learn more in the Doxygen documentation.

Changing the Codebase
---------------------

Flexibility and extensibility were high design goals for developing the megafauna model.
Hopefully you will find it possible to implement the necessary code changes/extensions for your particular research questions.
You will need basic skills with Git and C++ (C++11 standard) in order to contribute.

On the index/main page of the Doxygen documentation you will be directed to the resources you need to contribute.
Please also read through the file [CONTRIBUTING.md](CONTRIBUTING.md).

Continuous Integration (CI) runs through gitlab.com in this mirror repository: <https://gitlab.com/wtraylor/modular_megafauna_model>

Note that for running the model, you don’t need to change the source code.
Most parameters can be set in the instruction file.

Bugs and Issues
---------------

Known bugs and improvements are collected in the [GitHub issue tracker](https://github.com/wtraylor/modular_megafauna_model/issues/).
If you discover a new bug, please use the issue tracker to [report](https://github.com/wtraylor/modular_megafauna_model/issues/new) it.
If you can fix it yourself, fork this repository, fix the bug, and request a pull into the main repository.
Compare [Understanding the GitHub flow](https://guides.github.com/introduction/flow/).

Authors
-------

- Wolfgang Traylor (wolfgang.traylor@senckenberg.de) ![ORCID](docs/images/orcid.png) <https://orcid.org/0000-0002-4813-1072>, Senckenberg Biodiversity and Climate Research Centre ([SBiK-F][])

[SBiK-F]: <https://www.senckenberg.de/en/institutes/sbik-f/>

Similar Projects
----------------

- Dangal, Shree R. S., Hanqin Tian, Chaoqun Lu, Wei Ren, Shufen Pan, Jia Yang, Nicola Di Cosmo, and Amy Hessl. 2017. “Integrating Herbivore Population Dynamics into a Global Land Biosphere Model: Plugging Animals into the Earth System.” Journal of Advances in Modeling Earth Systems 9 (8): 2920–45. <https://doi.org/10.1002/2016MS000904>.
- Illius, A. W., and T. G. O’Connor. 2000. “Resource Heterogeneity and Ungulate Population Dynamics.” Oikos 89 (2): 283–94. <https://doi.org/10.1034/j.1600-0706.2000.890209.x>.
- Pachzelt, Adrian, Anja Rammig, Steven Higgins, and Thomas Hickler. 2013. “Coupling a Physiological Grazer Population Model with a Generalized Model for Vegetation Dynamics.” Ecological Modelling 263: 92–102. <https://doi.org/http://dx.doi.org/10.1016/j.ecolmodel.2013.04.025>.
- Zhu, Dan, Philippe Ciais, Jinfeng Chang, Gerhard Krinner, Shushi Peng, Nicolas Viovy, Josep Peñuelas, and Sergey Zimov. 2018. “The Large Mean Body Size of Mammalian Herbivores Explains the Productivity Paradox During the Last Glacial Maximum.” Nature Ecology & Evolution. <https://doi.org/10.1038/s41559-018-0481-y>.

License
-------

- This project follows the [REUSE][] standard:
    - Every file has its copyright/license information either in a comment at the top or in a separate text file with the extension `.license`.
    - All license texts can be found in the directory `LICENSES/`.
    - Project information and licenses for Git submodules can be found in the text file `.reuse/dep5`.
- The Modular Megafauna Model is Free Software under the [GNU Lesser General Public License v3.0 or later][lgpl].
- The software documentation, the accompanying images, and configuration files are licensed under the [Creative Commons Attribution 4.0 International][cc-by-4.0].

[cc-by-4.0]: https://creativecommons.org/licenses/by/4.0
[lgpl]: https://www.gnu.org/licenses/lgpl-3.0-standalone.html

### Third-Party Work

- The [Catch2](https://github.com/catchorg/Catch2) test framework (`tests/catch.hpp`) is licensed under the [Boost Software License](http://www.boost.org/LICENSE_1_0.txt).

- The [CMake-codecov](https://github.com/RWTH-HPC/CMake-codecov) scripts in `external/CMake-codecov/` are licensed under the BSD-3-Clause License.

- The version of the [PlantUML](http://plantuml.com) file (`docs/plantuml.jar`), which is used to render UML diagrams in the Doxygen documentation, is under the [MIT license](http://opensource.org/licenses/MIT).

- The [cpptoml](https://github.com/skystrife/cpptoml) library (`external/cpptoml/`) is licensed under the [MIT license](http://opensource.org/licenses/MIT).
<!--
SPDX-FileCopyrightText: 2020 Wolfgang Traylor <wolfgang.traylor@senckenberg.de>

SPDX-License-Identifier: CC-BY-4.0
-->

# Contributing to the Codebase

This document is a guide for developers who want to contribute to the megafauna library.
It explains the organization of the software project and sets guidelines for code style, structure, and format.
This document is only about *syntax* only.

After contributing something, don’t forget to add your name to:

- the file header in a new line starting with `SPDX-FileCopyrightText: ...` (following the [REUSE][] standard),
- the “Authors” section in the `README.md`, and
- the list of authors in the citation files [CITATION.cff](CITATION.cff) and [codemeta.json](codemeta.json).

Note that the authors list in the [Zenodo archive][] is automatically derived from the contributors in the Git history.

[REUSE]: https://reuse.software
[Zenodo archive]: <https://zenodo.org/badge/latestdoi/228426088>

## Table of Contents

<!-- vim-markdown-toc GFM -->

* [Version Control](#version-control)
    * [Branches](#branches)
    * [Release Versioning](#release-versioning)
        * [Checklist for merging into `master`](#checklist-for-merging-into-master)
    * [Commit Messages](#commit-messages)
    * [Continuous Integration](#continuous-integration)
* [Licensing](#licensing)
* [Coding Guidelines](#coding-guidelines)
    * [Repository Structure](#repository-structure)
    * [Code Format](#code-format)
        * [Naming Code Elements](#naming-code-elements)
        * [Ordering](#ordering)
        * [File Header](#file-header)
    * [Unit Tests](#unit-tests)
    * [Code Checkers](#code-checkers)
    * [Doxygen Documentation](#doxygen-documentation)
        * [Markdown](#markdown)
        * [BibTeX Bibliography](#bibtex-bibliography)

<!-- vim-markdown-toc -->

## Version Control

### Branches
- The `master` branch is reserved for stable releases, tagged with the version numbers.
- This repository follows Vincent Driessen’s [Successful Git Branching Model](https://nvie.com/posts/a-successful-git-branching-model/).
- If you are new to Git branching, check out this tutorial: [Learn Git Branching](https://learngitbranching.js.org/)

### Release Versioning
This project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html):

- The pattern of realease tas is `MAJOR.MINOR.PATCH` (e.g., `1.9.1` or `1.12.0`).
- If the new version cannot read old instruction files anymore or breaks the library interface, increment `MAJOR`.
- If the new version introduces a new feature, but still interoperates like the old version, increment `MINOR`.
- If the new version only fixes bugs, extends or amends the documentation or refactors code, increment `PATCH`.

#### Prereleases
Prerelease versions stay on the `develop` branch.
You can create tags, but don’t merge into the `master` branch.

The tag name should have some meaningful and numbered appendix to the release version it is moving towards.
For example, if the prerelease prepares version `1.3.0`, the tag may be `1.3.0-alpha.1` in the early development phase.
For final testing you may use `1.3.0-beta.1`, `1.3.0-beta.2`, and so on.

Don’t write the prerelease version into `CMakeLists.txt`, `codemeta.json`, and `CITATION.cff`.
That would be difficult to maintain because if you don’t revert the version number in these files to `0.0.0` (which indicates development) immediately and just continue to commit, you end up with different commits containing the same version number in the metadata.
The version numbers in `CMakeLists.txt`, `codemeta.json`, and `CITATION.cff` are really only for actual releases, which are citeable and documented on Read the Docs.

#### Checklist for merging into `master`

Each merge into the `master` branch is a release and should have a Git tag.

1. List your changes in `CHANGELOG.md`, following the formatting guidelines there.
    - Rename the “Unreleased” section to the to-be-released version in `CHANGELOG.md`.
    - At the bottom of the file, add the URL for the release, comparing it to the previous version. Orient yourself by the existing link URLs.
2. Set the new version in `CMakeLists.txt` under `VERSION`.
3. Update metadata files:
    - Set the version and the `date-released:` field in `CITATION.cff`. The date format is `YYYY-MM-DD`.
    - Set the `"version":` and `"dateModified":` fields in `codemeta.json`.
4. Now do the merge: `git switch master && git merge --no-ff develop`
5. Create a new release on GitHub, which will trigger [Zenodo](https://zenodo.org) to archive the code and mint a DOI.
    - The release and the tag description should summarize the changes (which you can copy-paste from `CHANGELOG.md`.
    - The name of the tag and the release is just the exact version, e.g. `0.1.2`.
6. Fast-forward the `develop` branch: `git switch develop && git merge --ff master`
7. Your first commit in `develop` resets everything so that it cannot be confused with a released version:
    - Set `VERSION 0.0.0` in `CMakeLists.txt`.
    - Set `version: 0.0.0` in `CITATION.cff`, and empty the `date-released:` field.
    - In `codemeta.json` set `"version": "0.0.0"` and `"dateModified": ""`.
    - Prepare the `[Unreleased]` section in `CHANGELOG.md`. Update the URL for `[Unreleased]` in the link list of the bottom; it should compare `develop` with the latest release.
8. Check that Zenodo and Read the Docs have received the latest version:
    - Zenodo: <https://doi.org/10.5281/zenodo.4710254>
    - Read the Docs: <https://modular-megafauna-model.readthedocs.io/en/latest/>
9. If applicable: Close the [Milestone][] for this release on GitHub.
10. Announce the release in the [Matrix channel][].

[Matrix channel]: <https://matrix.to/#/!rnevkLtJTORmvyzFHD:matrix.org?via=matrix.org>
[Milestone]: <https://github.com/wtraylor/modular_megafauna_model/milestones>

### Commit Messages
Follow Chris Beams’ guide for crafting your Git commit messages: [How to Write a Git Commit Message](https://chris.beams.io/posts/git-commit/)

> 1. Separate subject from body with a blank line
> 1. Limit the subject line to 50 characters
> 1. Capitalize the subject line
> 1. Do not end the subject line with a period
> 1. Use the imperative mood in the subject line
> 1. Wrap the body at 72 characters
> 1. Use the body to explain what and why vs. how

### Continuous Integration
The repository should always be in a valid state.
A number of tests are defined in the file `.gitlab-ci.yml`.
This file works with GitLab Continuous Integration (CI).

The CI script also runs [Valgrind](https://valgrind.org) memory check.
With the bash script `tools/run_valgrind_memcheck` you can execute a memory check manually on your local machine.
**Always make sure contributions to the codebase don’t have memory leaks.**

## Licensing
- Familiarize yourself with the REUSE standard in this tutorial: <https://reuse.software/tutorial/>
- When you create a new file, add a REUSE license header with the same license as similar files in the project.
- When you contribute to a file, add yourself as a copyright holder to the REUSE license header.
- When you create a commit with Git, use the `-s/--signoff` flag in order to sign the [Developer Certificate of Origin][DCO]. This way you certify that you wrote or otherwise have the right to submit the code you’re contributing to the project.
    - Just come into the habit of writing `git commit -s`.

[DCO]: https://developercertificate.org/

## Coding Guidelines

### Repository Structure
This project follows the [Pitchfork Layout](https://github.com/vector-of-bool/pitchfork) for C++ projects.
Here is a summary of the relevant parts:

- **Namespace Folders:** The `src/` directory has subfolders reflecting the namespaces of the contained components. To minimize the danger of name collision, the header **include guards** contain the namespace hierarchy also, e.g. `FAUNA_OUTPUT_HABITAT_DATA_H`.
- **Separate Header Placement:** Header (`*.h`) and source (`*.cpp`) files are kept together in `src/` if they are *private* (not part of the library interface). *Public* headers are placed in `include/` while their corresponding source files remain in `src/`.
- **Coherence:** A header and corresponding source file (= *physical component*)contain code for *logical component.* If in doubt, rather air on the side of granularity and create several individual components.
- **Merged Test Placement:** Any logical/physical unit has its unit test in a file in the same folder with the same file name, but with the suffix `.test.cpp`.

### Code Format
The C++ code layout follows the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).
No worries, you don’t need to read everything.
Just see that your IDE or text editor auto-formats your code using [Clang](http://clang.llvm.org/).
The default Google layout is exported in the [.clang-format](.clang-format) file, ready to be read by `clang`.

This repository has a [.editorconfig](.editorconfig) file to define indentation rule for different file types.
Please install the plugin for your text editor if available: [editorconfig.org/](https://editorconfig.org/)

#### Naming Code Elements

- **Files** are always lower-case with underscores to separate words.
    - Header files end with `.h`, source files with `.cpp`, and the corresponding unit test files with `.test.cpp`.
    - If a file only contains one class, name the file like the class.
    - If a file contains several classes, use a plural like `net_energy_models.h`.
    - If a file contains a collection of functionality, use an abstract grouping noun, e.g. `stochasticity.h` or `nitrogen.h`.
- **Classes** are named in CamelCase with upper-case first letter, e.g. `MyExampleClass`. Don’t repeat the namespace in the class name (avoid something like `Output::OutputDataClass`).
    - **Enum** types are like classes.
- **Functions** are imperative verbs with underscores, e.g. `create_new_herbivores()`.
- Global **constants** as well as static const member and function variables are all-uppercase with underscores, e.g. `MY_GLOBAL_CONSTANT`.
    - C++11-style **enum class elements** don’t have global scope and thus don’t require a prefix. Since the shouting tone of all-uppercase names is distracting, just use CamelCase for the enum members, e.g. `OutputInterval::Annual`.
- **Namespaces** are short and lower-case with first letter capitalized, e.g. `Fauna`.

#### Ordering
An example class definition in a header file:
```cpp
class MyExampleClass{
 public: // -> public members first
  /// Constructor
  MyExampleClass(); // -> constructors are always first

  /// Destructor
  ~MyExampleClass(); // -> destructor comes next

  // -> The following member functions in alphabetical order:
  /// Create new herbivore instance.
  Herbivore* create_herbivore();

  /// Get the type of herbivore.
  HerbivoreType get_herbivore_type();

  /// A public global constant in this class.
  static const int MY_CONSTANT = 10;

 private:
  /// Helper function to perform calculations.
  void my_private_function();

  // -> Finally the private member variables NOT in alphabetical order, but in
  // the order of desired initialization.
  HerbivoreType herbi_type;
  int var = 10;
};
```

In the corresponding source file, *all* function definitions (both private & public) are in alphabetical order, except for the constructors and destructor, which come first.

If there is more than one class in the header file, separate their function definitions in blocks with big comment captions, for example like this:
```cpp
//============================================================
// HftPopulationsMap
//============================================================
```

If there are functions local to this file, put them in an anonymous namespace before any other definitions:
```cpp
namespace {
  int my_local_function(){
  // ...
  };
}
```

#### File Header
Begin each `.h` or `.cpp` file with a doxygen header containing a brief description.
The description will appear in the file list of the generated doxygen documentation.
Ususally the brief description will be the same for a `.h` and its `.cpp` file.
If it is only one class in the header file, you can also copy the `\brief` description from that class.
Here is an example:

<!--TODO: REUSE header-->

```cpp
/**
 * \file
 * \brief Management classes of herbivore populations.
 * \copyright LGPL-3.0-or-later
 * \date <current year>
 */
```

We omit the `\author` field because it might be difficult to keep track of all authors who have contributed.
(That’s what version control is for.)
Instead, all contributors of the project shall be listed collectively in the “Authors” section of the `README.md`.

The `\date` field is only relevant for the copyright. (Use Git to see when the file has been changed.)

### Unit Tests
Make sure to write a unit test for every logical component.
If you create a `.cpp` file, there should most likely also be a corresponding `.test.cpp` file that checks the public functions of the class or classes.

Unit tests use the [Catch2](https://github.com/catchorg/Catch2) framework in the [single header](https://raw.githubusercontent.com/catchorg/Catch2/master/single_include/catch2/catch.hpp) distribution.
(The `tests/catch.hpp` file can be updated from time to time.)
To run the unit tests after building the megafauna library, run `./megafauna_unit_tests` in the build directory.
To disable building the unit tests, you can call `cmake -DBUILD_TESTING=OFF /path/to/repo`.

### Code Checkers

Run [cppclean](https://github.com/myint/cppclean) on the code to find unnecessary `#include`s, unnecessary functions, and a lot more.
Execute the helper script `./tools/cppclean.sh` in the Bash.

### Doxygen Documentation

- The documentation is completely in English, preferably with American spelling.
- Images are placed in `docs/images/`. If the figure was plotted with a script, save the script in the same folder and with the same file name as the image.

#### Markdown
- Overview pages are written in [Markdown](http://www.doxygen.nl/manual/markdown.html) in `*.md`-files in the `docs/` folder.
- Make a new line after each sentence (and perhaps after logical sentence structures): [Inner Sentence Rule](https://cirosantilli.com/markdown-style-guide#option-wrap-inner-sentence).

#### BibTeX Bibliography
[BibTeX](www.bibtex.org) is used for bibliographic references: [docs/bibliography.bib](docs/bibliography.bib).
The Doxygen command `\cite` is used for that.
This makes browsing the Doxygen documentation easier.

In general you should not need to put any references to scientific publications in *comments* in the source code.
Better you explain everything in a narrative form in the Doxygen *documentation* and use the `\cite` command for that.
If you do cite in source code comments, make sure that the reference is uniquely identifiable in `bibliography.bib`.

Write the in-text citation in the doxygen documentation (either in a C++ file or in a Markdown document) in the [APA format](https://www.mendeley.com/guides/apa-citation-guide#2_In_Text) followed by the `\cite` reference:
```plain
Illius & O’Connor (2000) \cite illius2000resource states that ...
Blaxter (1989, p. 123) \cite blaxter1989energy states that ...
(McDonald et al. 2010, p. 123 \cite mcdonald2010animal)
```

Use [bibsort](http://ftp.math.utah.edu/pub/bibsort/) to sort the
bibliography entries by label.
Or you can do it manually, too.
A sorted bibliography is easier to read and better to maintain with
version control software.

- Be parsimonious! Some fields are not handled by Doxygen.
- Include the abstract whenever possible.
- Don’t use `journaltitle` and `date`. They are not recognized by Doxygen. Use `year`, `month`, etc. instead of `date`.
- Use `{...}` brackets instead of `""`.
- Have equal signs (`=`) line up vertically (for prettiness).
- BibTeX identifiers (inspired by the BibTeX export of [Google Scholar](https://scholar.google.com)):
`<author><year><firstword>` (all lowercase and without delimiter)
    - author: Family name of first author as it would be cited (including van/von/…)
    - year: Publication year.
    - firstword: First word of the title excluding ‘the’, ‘a’, ‘an’, ‘of’, ‘is’, ‘are’, ‘were’, and the like. Hyphens, dashes, apostrophes, and slashes within the first (compound) word are simply omitted.
    - If the above produces non-unique IDs, use the second word, or (if even that fails) the third.
# cpptoml
A header-only library for parsing [TOML][toml] configuration files.

Targets: [TOML v0.5.0][currver] as of August 2018.

This includes support for the new DateTime format, inline tables,
multi-line basic and raw strings, digit separators, hexadecimal integers,
octal integers, binary integers, and float special values.

Alternatives:
- [Boost.toml][boost.toml] is a C++ implementation of a TOML parser using
  the Boost library. As of writing, it supports v0.5.0 as well.
- [ctoml][ctoml] is a C++11 implementation of a TOML parser, but only
  supports v0.2.0.
- [libtoml][libtoml] is a C implementation of a TOML parser, which can be
  linked to from your C++ programs easily. As of April 2016, it supports
  v0.4.0.
- [tinytoml][tinytoml] is a C++11 implementation of a TOML parser, which
  also supports v0.4.0 as of November 2015.

## Build Status
[![Build Status](https://travis-ci.org/skystrife/cpptoml.svg?branch=master)](https://travis-ci.org/skystrife/cpptoml)

## Test Results

From [the toml-test suite][toml-test]:

```
126 passed, 0 failed
```

We also currently maintain (but hopefully not indefinitely!) a [fork of the
toml-test suite][toml-test-fork] that adds tests for features and
clarifications that have been added to the TOML spec more recently than
toml-test has been updated. We pass every test there.

```
148 passed, 0 failed
```

# Compilation
Requires a well conforming C++11 compiler. On OSX this means clang++ with
libc++ and libc++abi (the default clang installed with XCode's command line
tools is sufficient).

On Linux, you should be able to use g++ >= 4.8.x, or clang++ with libc++
and libc++abi (if your package manager supplies this; most don't).

Compiling the examples can be done with cmake:

```
mkdir build
cd build
cmake ../
make
```

# Example Usage
To parse a configuration file from a file, you can do the following:

```cpp
auto config = cpptoml::parse_file("config.toml");
```

`parse_file()` returns a (shared pointer to a) `cpptoml::table`, which you
can then query. It will throw an instance of `cpptoml::parse_exception` in
the event that the file failed to parse, and the exception message should
contain the line number the error occurred as well as a description of the
error.

## Obtaining Basic Values
You can find basic values like so:

```cpp
auto val = config->get_as<int64_t>("my-int");
// val is a cpptoml::option<int64_t>

if (val)
{
    // *val is the integer value for the key "my-int"
}
else
{
    // "my-int" either did not exist or was not an integer
}
```

To simplify things, you can specify default a default value using the
`value_or` function on the `option`:

```cpp
auto baz = config->get_as<double>("baz").value_or(0.5);
// baz is now the double value for key "baz", if it exists, or 0.5 otherwise
```

cpptoml has extended support for dates and times beyond the TOML v0.4.0
spec. Specifically, it supports

- Local Date (`local_date`), which simply represents a date and lacks any time
  information, e.g. `1980-08-02`;
- Local Time (`local_time`), which simply represents a time and lacks any
  date or zone information, e.g. `12:10:03.001`;
- Local Date-time (`local_datetime`), which represents a date and a time,
  but lacks zone information, e.g. `1980-08-02T12:10:03.001`;
- and Offset Date-time (`offset_datetime`), which represents a date, a
  time, and timezone information, e.g. `1980-08-02T12:10:03.001-07:00`

Here are the fields of the date/time objects in cpptoml:

- year (`local_date`, `local_datetime`, `offset_datetime`)
- month (`local_date`, `local_datetime`, `offset_datetime`)
- day (`local_date`, `local_datetime`, `offset_datetime`)
- hour (`local_time`, `local_datetime`, `offset_datetime`)
- minute (`local_time`, `local_datetime`, `offset_datetime`)
- second (`local_time`, `local_datetime`, `offset_datetime`)
- microsecond (`local_time`, `local_datetime`, `offset_datetime`)
- hour\_offset (`offset_datetime`)
- minute\_offset (`offset_datetime`)

There are convenience functions `cpptoml::offset_datetime::from_zoned()` and
`cpptoml::offset_datetime::from_utc()` to convert `struct tm`s to
`cpptoml::offset_datetime`s.

## Nested Tables
If you want to look up things in nested tables, there are two ways of doing
this. Suppose you have the following structure:

```toml
[first-table]
key1 = 0.1
key2 = 1284

[first-table.inner]
key3 = "hello world"
```

Here's an idiomatic way of obtaining all three keys' values:

```cpp
auto config = cpptoml::parse_file("config.toml");
auto key1 = config->get_qualified_as<double>("first-table.key1");
auto key2 = config->get_qualified_as<int>("first-table.key2");
auto key3 = config->get_qualified_as<std::string>("first-table.inner.key3");
```

(Note that, because the TOML spec allows for "." to occur in a table name,
you won't *always* be able to do this for any nested key, but in practice
you should be fine.)

A slightly more verbose way of getting them would be to first obtain the
individual tables, and then query those individual tables for their keys
like so:

```cpp
auto config = cpptoml::parse_file("config.toml");

auto first = config->get_table("first-table");
auto key1 = first->get_as<double>("key1");
auto key2 = first->get_as<int>("key2");

auto inner = first->get_table("inner");
auto key3 = inner->get_as<std::string>("key3");
```

The function `get_table_qualified` also exists, so obtaining the inner
table could be written as

```cpp
auto inner2 = config->get_table_qualified("first-table.inner");
```

## Arrays of Values
Suppose you had a configuration file like the following:

```toml
arr = [1, 2, 3, 4, 5]
mixed-arr = [[1, 2, 3, 4, 5], ["hello", "world"], [0.1, 1.1, 2.1]]
```

To obtain an array of values, you can do the following:

```cpp
auto config = cpptoml::parse_file("config.toml");

auto vals = config->get_array_of<int64_t>("arr");
// vals is a cpptoml::option<std::vector<int64_t>>

for (const auto& val : *vals)
{
    // val is an int64_t
}
```

`get_array_of` will return an `option<vector<T>>`, which will be empty if
the key does not exist, is not of the array type, or contains values that
are not of type `T`.

For nested arrays, it looks like the following:

```cpp
auto nested = config->get_array_of<cpptoml::array>("mixed-arr");

auto ints = (*nested)[0]->get_array_of<int64_t>();
// ints is a cpptoml::option<std::vector<int64_t>>

auto strings = (*nested)[1]->get_array_of<std::string>();
auto doubles = (*nested)[2]->get_array_of<double>();
```

There is also a `get_qualified_array_of` for simplifying arrays located
inside nested tables.

## Arrays of Tables
Suppose you had a configuration file like the following:

```toml
[[table-array]]
key1 = "hello"

[[table-array]]
key1 = "can you hear me"
```

Arrays of tables are represented as a separate type in `cpptoml`. They can
be obtained like so:

```cpp
auto config = cpptoml::parse_file("config.toml");

auto tarr = config->get_table_array("table-array");

for (const auto& table : *tarr)
{
    // *table is a cpptoml::table
    auto key1 = table->get_as<std::string>("key1");
}
```

## More Examples
You can look at the files files `parse.cpp`, `parse_stdin.cpp`, and
`build_toml.cpp` in the root directory for some more examples.

`parse_stdin.cpp` shows how to use the visitor pattern to traverse an
entire `cpptoml::table` for serialization.

`build_toml.cpp` shows how to construct a TOML representation in-memory and
then serialize it to a stream.

[currver]: https://github.com/toml-lang/toml/blob/master/versions/en/toml-v0.5.0.md
[toml]: https://github.com/toml-lang/toml
[toml-test]: https://github.com/BurntSushi/toml-test
[toml-test-fork]: https://github.com/skystrife/toml-test
[ctoml]: https://github.com/evilncrazy/ctoml
[libtoml]: https://github.com/ajwans/libtoml
[tinytoml]: https://github.com/mayah/tinytoml
[boost.toml]: https://github.com/ToruNiina/Boost.toml
# Tutor {#page_tutor}
<!--
SPDX-FileCopyrightText: 2020 Wolfgang Traylor <wolfgang.traylor@senckenberg.de>

SPDX-License-Identifier: CC-BY-4.0
-->

\brief Instructions how expand the code base for your own needs.

\tableofcontents

## Herbivores Tutorials {#sec_tutor_herbivores}

### How to add a new herbivore class {#sec_new_herbivore_class}

The model design allows a complete substitution of the herbivore class.
If you want to implement a completely new model behaviour, you can derive your new class from \ref Fauna::HerbivoreInterface and write it from scratch.
If you want to build upon the base functionality, derive it from \ref Fauna::HerbivoreBase.

Then, derive a new class from \ref Fauna::PopulationInterface to manage and construct your object instances.
In \ref Fauna::WorldConstructor::create_populations(), create all instances of that population class for one habitat.

@startuml "Relationships for a new herbivore type."
	!include diagrams.iuml!new_herbivore_type
@enduml

### How to add a new energy expenditure component {#sec_new_expenditure_component}
- Add a new enum entry in \ref Fauna::ExpenditureComponent.
- TOML instruction file: Add new possible string value for the HFT parameter `expenditure.components` in \ref Fauna::InsfileReader::read_hft(); include it in the error message.
- Implement your algorithm as a free function or a class. See \ref expenditure_components.h for examples.
- Call your model in \ref Fauna::HerbivoreBase::get_todays_expenditure().
- Update the UML diagram in Section \ref sec_herbivorebase.

### How to add a new foraging limit {#sec_new_foraging_limit}
A foraging limit constrains the daily uptake of forage mass by a herbivore individual.
Foraging limits are implemented as functors (without using the [strategy design pattern](\ref sec_strategy), though).
Which ones are activated is defined by `foraging.limits` in \ref Fauna::Hft.
They are called in \ref Fauna::GetForageDemands::get_max_foraging().

- Add a new enum entry in \ref Fauna::ForagingLimit.
- TOML instruction file: Add a new possible string value for the HFT parameter `foraging.limits` in \ref Fauna::InsfileReader::read_hft()
- Implement your foraging limit (preferably as a function object in the file \ref foraging_limits.h, but you can do as you wish).
Make sure that an exception is thrown if it is called with an unknown forage type.
- Call your implementation in \ref Fauna::GetForageDemands::get_max_foraging().
- Update the UML diagram in \ref sec_herbivorebase.

### How to add a new digestive limit {#sec_new_digestive_limit}
On top of foraging limitations, the daily forage uptake can be constrained by maximum digestive throughput.
The implementation is almost parallel to a [foraging limit](\ref sec_new_foraging_limit).

- Add a new enum entry in \ref Fauna::DigestiveLimit.
- Read the new value for `digestion.limit` from the TOML instruction file in \ref Fauna::InsfileReader::read_hft().
- Implement your digestive limit algorithm as a free function or an object. If it is not much code, put it in \ref foraging_limits.h, otherwise create a new file for it.
Make sure that an exception is thrown if it is called with an unknown forage type.
- Call your implementation in \ref Fauna::GetForageDemands::get_max_digestion().
- Update the UML diagram in \ref sec_herbivorebase.

### How to add a new reproduction model {#sec_new_reproduction_model}
A reproduction model defines the offspring per female individual for each simulation day.

- Create a new enum entry in \ref Fauna::ReproductionModel.
- Read the new value for `reproduction.model` from the instruction file in \ref Fauna::InsfileReader::read_hft().
- Create your class or function in \ref reproduction_models.h or in a separate file.
- Call your model in \ref Fauna::HerbivoreBase::get_todays_offspring_proportion().
- Update the UML diagram in \ref sec_herbivorebase.

### How to add a new diet composer {#sec_new_diet_composer}
In a scenario with multiple forage types, the herbivore decides what to include in its diet.
This decision is modelled by an implementation of a so called “diet composer model”: \ref Fauna::DietComposer.
You can implement your own model as a new class or a simple function; just call it in \ref Fauna::GetForageDemands::get_diet_composition().

- Create a new enum entry in \ref Fauna::DietComposer.
- Read the new value for `foraging.diet_composer` in \ref Fauna::InsfileReader::read_hft().
- Call your model in \ref Fauna::GetForageDemands::get_diet_composition().
- Update the UML diagram in \ref sec_herbivorebase.

### How to add a new mortality factor {#sec_new_mortality_factor}
Any death event of an herbivore is modelled by a mortality factor.
Whether you want to have herbivores die by for instance disease, drought, or predators, you should create a new mortality factor.

- Create a new enum entry in \ref Fauna::MortalityFactor.
- Parse the new possible value for the set `mortality.factors` in \ref Fauna::InsfileReader::read_hft().
- Implement your mortality model as a function or class in \ref mortality_factors.h or in a separate file (if it’s more complex).
- Call the mortality factor in \ref Fauna::HerbivoreBase::apply_mortality_factors_today().
- Update the UML diagram in \ref sec_herbivorebase.

## Forage Tutorials {#sec_tutor_forage}

### How to add a new forage type {#sec_new_forage_type}

- Create new enum entry in \ref Fauna::ForageType and add it to \ref Fauna::FORAGE_TYPES by expanding the initializing function `get_all_forage_types()`, which is declared in local namespace in \ref forage_types.cpp.

- Add a short name without blanks in \ref Fauna::get_forage_type_name().

- Derive new class from \ref Fauna::ForageBase.
	+ Implement a `merge()` method, like \ref Fauna::GrassForage::merge().

- Add a new member variable in \ref Fauna::HabitatForage of that class.
	+ Add it in \ref Fauna::HabitatForage::operator[]().
	+ Call your `merge()` function in \ref Fauna::HabitatForage::merge().

- Adjust the implementation of \ref Fauna::Habitat::get_available_forage() and \ref Fauna::Habitat::remove_eaten_forage() in your vegetation model.

- For \ref Fauna::NetEnergyModel::GrossEnergyFraction you will need to define a value for the parameter `forage.gross_energy.new_forage_type` in the TOML file.

- Herbivores
	+ Check all foraging and digestion limits (\ref foraging_limits.h) whether they need to be expanded.
	+ Probably you will need to implement [a new diet composer](\ref sec_new_diet_composer) or adjust existing ones.

- Test Simulations
	+ If you want to use your forage type in the demo simulations, expand \ref Fauna::Demo::SimpleHabitat by a new growth model (analoguous to \ref Fauna::Demo::LogisticGrass).
	Also update the UML diagram in the class documentation of `SimpleHabitat`.

@startuml "Relationships for a new forage type."
	!include diagrams.iuml!new_forage_type
@enduml

### How to change forage net energy content {#sec_change_netenergy}

- Implement your model in a function or class in the file \ref net_energy_models.h.
- Add a new enum item in \ref Fauna::NetEnergyModel.
- Parse that value for the HFT parameter `digestion.net_energy_model` in \ref Fauna::InsfileReader::read_hft().
- Execute your function or class in \ref Fauna::HerbivoreBase::get_net_energy_content().
- Update the UML diagram in \ref sec_herbivorebase.

### How to add a new forage distribution algorithm {#sec_new_forage_distribution}

- Derive a new class from \ref Fauna::DistributeForage and implement your algorithm. If the algorithm is not too big, put it in \ref forage_distribution_algorithms.h, otherwise create a new set of files.
- Return a reference to an object of your class in \ref Fauna::WorldConstructor::create_distribute_forage() if it is selected in the parameters.
- Add a new enum entry in \ref Fauna::ForageDistributionAlgorithm.
- Parse your new value from the parameter `simulation.forage_distribution` in \ref Fauna::InsfileReader::read_hft().

## Parameters Tutorials {#sec_tutor_parameters}

### How to add a new global parameter {#sec_new_global_parameter}

- Declare your parameter as a member variable in \ref Fauna::Parameters and initialize it with a valid default value. The `enum class` type for Enum parameters should be declared in \ref parameters.h.
- Write a validity check in \ref Fauna::Parameters::is_valid().
- Parse the parameter from the TOML instruction file in \ref Fauna::InsfileReader. General parameters are parsed in Fauna::InsfileReader::read_table_simulation().
    + If you are creating a whole new set of parameters, it might make sense to group them in a TOML table. Existing tables are `output` and `simulation`. You should then write a new private member function similar to \ref Fauna::InsfileReader::read_table_simulation().

### How to add a new HFT parameter {#sec_new_hft_parameter}

- Declare a new member variable in \ref Fauna::Hft and initialize it with a valid default value. The `enum class` type for Enum parameters should be declared in \ref hft.h.
- Write a validity check in \ref Fauna::Hft::is_valid().

\see \ref sec_design_parameters

## Output Tutorials {#sec_output_tutor}

### How to add a new output variable {#sec_new_output_variable}
After separating the megafauna model from LPJ-GUESS into its own library only a minimal set of variables are made available for output.
If your variable of interest is already present in \ref Fauna::Output::HabitatData or \ref Fauna::Output::HerbivoreData, then you can skip the first step in the following tutorial.

- Extend the appropriate container: \ref Fauna::Output::HabitatData or \ref Fauna::Output::HerbivoreData by a new member variable and initialize it with zero.
    + Add it in the `reset()` function.
	+ Implement average building in the appropriate `merge()` function.
	+ For herbivore data, you need to add it to \ref Fauna::Output::HerbivoreData::create_datapoint(). If your value is *per individual*, you will need to weight the value by individual density; if it is *per area* or *per habitat*, you can calculate the sum.
	+ Assign a value to the variable somewhere in daily simulation.

- Write the variable in \ref Fauna::Output::TextTableWriter.
    + Add a selector for your new output file as a boolean member variable in \ref Fauna::Output::TextTableWriterOptions. Pay attention to place it in the right Doxygen group.
    + Parse the name of the selector in \ref Fauna::InsfileReader::read_table_output_text_tables() and add it there to the valid options in the error message.
    + Optional: Consider adding it in the example output in the file `examples/megafauna.toml` and plotting it in `tools/demo_simulator/demo_results.Rmd`. In that case the output file needs to be listed as an artifact in `.gitlab-ci.yml` under the job “demo_simulation”.
    + Add an output file stream (`std::ofstream`) for your variable as a private member variable in \ref Fauna::Output::TextTableWriter.
    + In the constructor \ref Fauna::Output::TextTableWriter::TextTableWriter(), add your new output file to the list of file streams if it is selected in the options.
    + Initialize the column captions of your new file in \ref Fauna::Output::TextTableWriter::write_captions().
    + Write the data to the file in \ref Fauna::Output::TextTableWriter::write_datapoint().

\note If you want to add a variable that is not *per herbivore mass*, you would have to use mass density as weight.

### How to write output to another format {#sec_new_output_writer}
The default output format are very simple tab-separated plaintext tables.
The class \ref Fauna::Output::TextTableWriter is responsible for this format.
However, you can also replace that output format with another one, for instance writing to a NetCDF file or forwarding it to another program or library.

- Derive a new class from \ref Fauna::Output::WriterInterface.
- Add a new enum entry to \ref Fauna::OutputFormat.
- Parse the new option in \ref Fauna::InsfileReader::read_table_output().
- Create a new instance of your writer class in the constructor \ref Fauna::World::World() if your enum entry is selected.

\see \ref sec_design_output design

-------------------------------------------------

\copyright <a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/80x15.png" /></a> This software documentation is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
\author Wolfgang Traylor, Senckenberg BiK-F
\date 2019
# Object-oriented Design Concepts {#page_object_orientation}
<!--
SPDX-FileCopyrightText: 2020 Wolfgang Traylor <wolfgang.traylor@senckenberg.de>

SPDX-License-Identifier: CC-BY-4.0
-->

\brief Some programming principles and paradigms in object-oriented software engineering.

\tableofcontents

In the megafauna model a couple of object-oriented design patterns were employed that are explained here along with general concepts of object-oriented programming.

\todo Add explanation of the PIMPL idiom (compilation firewall). Compare [How to implement the pimpl idiom by using unique_ptr - Fluent C++](https://www.fluentcpp.com/2017/09/22/make-pimpl-using-unique_ptr/).

## Good Programming Practice {#sec_good_practice}

### Information Hiding {#sec_information_hiding}
Information or data hiding means that any parts that may be subject to later changes in the software design must not be accessible from other modules or from clients.
Only the very necessary access may be granted in well-defined, minimal interfaces.

Always assume the worst: If access is given to any other part of the software, that part may change it in unpredictable ways!

Declaring class members `private` (**encapsulation**) is one way to data hiding.

### Rule of Three {#sec_rule_of_three}
If a class explicitely defines at least one of the following methods, it should most likely also define the other ones:

- Destructor
- Copy Constructor
- Copy Assignment Operator

\todo In C++11 the **Rule of Five** becomes relevant!

## S-O-L-I-D Design principles ## {#sec_design_solid}

### Single Responsibility Principle {#sec_single_responsibility}
A class should have only a single responsibility:
A class should have only one reason to change.

### Open/Closed Principle {#sec_open_closed}
A class/module/function should be open for extension, but closed for modification.

### Liskov’s Substitution Principle {#sec_liskov_substitution}
Objects in a program should be replaceable with instances of their subtypes without altering the correctness of that program.

### Interface Segregation Principle {#sec_interface_segregation}
Many client-specific interfaces are better than one general-purpose interface.

### Dependency Inversion Principle {#sec_dependency_inversion}
1. High-level modules should not depend on low-level modules. Both should depend on abstractions.
2. Abstractions should not depend on details. Details should depend on abstractions.

## Inversion of Control {#sec_inversion_of_control}
The design principle of inversion of control is also called *Hollywood Principle:* “Don’t call us, we’ll call you.”
An architecture following that principle is built around a generic framework which directs the flow of control by delegating tasks to various, interchangeable submodules.
This approach makes the system more modular and extensible.

@startuml
	hide members
	hide methods
	Framework ..> Client1 : <<create & call>>
	Framework ..> Client2 : <<create & call>>
@enduml

Inversion of control is related to the [Dependency Inversion Principle](\ref sec_dependency_inversion), which differs in that it is concerned about the relationship between high-level and low-level modules rather than the framework design.

### Dependency Injection
Dependency Injection is a major technique to implement inversion of control:
One object (the framework) supplies the dependencies for another object (the client), which then become a part of the client’s state.
This is a direct alternative to global objects/variables.

Two kinds of dependency injection are used
1. **Setter Injection:** A client receives its dependency *after* being constructed, via a setter method. This is dangerous, because until initialization through the setter method, the client might be in an invalid state.
2. **Constructor Injection:** Any object of the client class receives its dependency in the constructor.

For example: The \ref Fauna::HftList object is not a global variable, but is instead being passed down from \ref Fauna::World to other classes that need it.

In the case of \ref Fauna::Habitat, the \ref Fauna::World is oblivious to the concrete realization of the interface, which makes it possible to substitute parts of the model without changing the framework.

\attention Never use global variables. They make unit tests impossible and violate the principles of modularity and compartmentalization.

## Design Patterns

### Singleton {#sec_singleton}
A class is called *singleton* if it permits only one global instantiation in the program.
This object-oriented approach has advantages over global variables because it is more generally flexible and the time of instantiation is flexible.
However, it is only justified for a class that is at the very top of the management hierarchy in the software architecture.

The basic implementation is as follows:

    class MySingleton{
    public:
     static MySingleton& get_instance(){
 	     static MySingleton global_instance; // creates the instance on first call
 	     return global_instance;
     }
    private:
     MySingleton(){}                     // Constructor hidden from the outside
     MySingleton(MySingleton const&);    // deleted copy constructor
     void operator=(MySingleton const&); // deleted assignment constructor
    }

To get access to the instance or trigger the initial instantiation, use `MySingleton::get_instance();`.

Wrapping global variables and free functions into a singleton class is good, but it is better to *avoid singletons all together* and instead follow the principle of [Inversion of Control](\ref sec_inversion_of_control).

\note The only instance where the Singleton pattern is used is in the demo simulator: \ref Fauna::Demo::Framework.

### Strategy {#sec_strategy}
The strategy design pattern defines a set of interchangable algorithms, each of which are encapsulated in a class that is derived from one abstract interface parent class.
Thanks to C++ polymorphism, the used algorithm can change during runtime.
Here is a basic implementation:

    struct Strategy {
    	virtual operator()(const int param) const = 0;
    };

    struct StrategyOne: public Strategy {
    	virtual operator()(const int param) const{ /*...*/ };
    };

    struct StrategyTwo: public Strategy {
    	virtual operator()(const int param) const{ /*...*/ };
    };

@startuml
	hide members
	hide methods
	interface Strategy
	Strategy <|-- StrategyOne
	Strategy <|-- StrategyTwo
	Client --> Strategy : <<use>>
@enduml

\anchor sec_functors
In this example, the classes are implemented as *function objects* (**functors**) because they overload the `operator()`.
In contrast to simple functions, they could also hold state variables and make use of class inheritance.
(They should still be light-weight, though.)
Their implemented algorithm can be called by using the object like a function:

    Strategy* do_algorithm = new StrategyOne; // just substitute the instantiated class here

    int i = 1;
    do_algorithm(i); // this line does not need to change when the strategy is replaced.

**Naming Conventions:**
Obviously, the class name (and the names of their instances and pointers) should be verbs.
The class name should of course be capitalized.

Strictly speaking, the strategy pattern aims to offer the possibility to substitute an algorithm during the lifetime of a client object.
In the megafauna model that is usually not the case, but rather it is used as a means of [dependency injection](sec_dependency_inversion).

### Facade {#sec_facade}
A facade class presents a simple interface to interact with another module or subsystem of the software.
The complexity of the subsystem is hidden behind the public functions of the facade class.
The subsystem doesn’t know about the facade.

@startuml "Structure of the facade design pattern."
	hide members
	hide methods
	client ..> facade : <<use>>
	facade ..> subsystem.class1 : <<call>>
	facade ..> subsystem.class2 : <<call>>
	facade ..> subsystem.class3 : <<call>>
@enduml

-------------------------------------------------

\copyright <a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/80x15.png" /></a> This software documentation is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
\author Wolfgang Traylor, Senckenberg BiK-F
\date 2019
# Recipes {#page_recipes}
<!--
SPDX-FileCopyrightText: 2020 Wolfgang Traylor <wolfgang.traylor@senckenberg.de>

SPDX-License-Identifier: CC-BY-4.0
-->

\brief Instruction file snippets and templates.

\tableofcontents

## Constant Population {#sec_recipe_constant_population}

Prescribing a constant herbivore density to the model may be used to simulate the effects of different stocking regimes.
Reproduction should then be disabled.
The herbivores may be immortal so that the stocking rate remains constant even if the animals are starving.

This is a minimal example herbivore instruction file for cattle:

```toml
[[hft]]
name = "Cattle"
[hft.body_fat]
deviation          = 0.125 # body condition
maximum            = 0.3   # kg/kg
maximum_daily_gain = 0.0   # kg/kg/day
[hft.body_mass]
female = 400 # kg
male   = 400 # kg
[hft.digestion]
type = "Ruminant"
limit = "IlliusGordon1992"
[hft.establishment]
age_range = { first = 1, last = 1 } # years
density   = 1.0 # ind/km²
[hft.expenditure]
components = [ "Taylor1981" ]
[hft.foraging]
diet_composer           = "PureGrazer"
limits                  = [ "IlliusOConnor2000" ]
half_max_intake_density = 40 # gDM/m²
[hft.reproduction]
model = "None"
```

The parameter values are mostly taken from Illius & O’Connor (2000)\cite illius2000resource.

With this instruction file, one immortal cohort of adult body mass will be established.
The mass density may fluctuate because of fat mass, but individual density will remain constant.

The `establishment.density` is probably the most crucial parameter in this setup.
Currently, the density can only prescribed globally and is constant.

-------------------------------------------------

\copyright <a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/80x15.png" /></a> This software documentation is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
\author Wolfgang Traylor, Senckenberg BiK-F
\date 2019
Welcome Page
============
<!--
SPDX-FileCopyrightText: 2020 Wolfgang Traylor <wolfgang.traylor@senckenberg.de>

SPDX-License-Identifier: CC-BY-4.0
-->

This source code documentation is a guide for users and developers of the megafauna library.
Before you continue, make sure you have read the `README.md` in the root of this repository in order to gain a general overview.

This documentation contains auto-generated source code documentation and in addition a few pages that generally describe the model and the source code architecture.
You can find all of these manually written pages under “Related Pages” in the upper navigation bar.

A good first read for users and developers is the \ref page_model.

Throughout this documentation you will find diagrams in UML (Unified Modeling Language), mostly class diagrams.
If you are not familiar with UML at all, consider looking at one of the many introductions you can find online.

## Contributing

If you are considering to contribute to the source code, make yourself familiar with the \ref page_design as well as the \ref page_object_orientation.
Please observe the [Open-closed Principle](\ref sec_open_closed) wherever possible: Only change a class *if necessary.* Try to extend it instead or—even better—write a new class with your functionality.
In \ref page_tutor you might find already a step-by-step instruction for your use case.
In `CONTRIBUTING.md` you will find all the details about what your code should look like and how you can commit your changes to the repository.

Please keep the documentation pages \ref page_design and \ref page_model updated.
The UML diagrams are rendered by PlantUML and collected all in one file: `docs/diagrams.iuml`.
You will find references for the PlantUML syntax on <http://www.plantuml.com>.

Come into the practice of **test-driven development**.
Write a test for each one of your new classes!
Check all of its public methods and see what example test scenarios it has to fulfill.
Also check if exceptions are thrown correctly.

A scientifically important paradigm of this library is that **no parameters are hard-coded.**
Any biologically relevant parameter value is subject to uncertainty.
In order to be included in a sensitivity/uncertainty analysis, the parameter must be in the instruction file.

## Integration into a Vegetation Model

If you want to introduce the megafauna library to a vegetation model, you should study the code of the demo simulator (`tools/demo_simulator`) and its documentation on the page \ref page_demo_simulator.

Since the MMM is licensed under the [LGPLv3][LGPL] (or later), you can use it in differently licensed (even proprietary) vegetation models, as long as you …

1. give prominent notice about the megafauna library and its license,
2. attach the LGPL license, and
3. convey the library source code.

But see the [LGPL][] § 4 for all the legal details on “Combined Work.”

[LGPL]: https://www.gnu.org/licenses/lgpl-3.0.html

-------------------------------------------------

\copyright <a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/80x15.png" /></a> This software documentation is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
\author Wolfgang Traylor, Senckenberg BiK-F
\date 2019
# Ecological Model Discussion {#page_model}
<!--
SPDX-FileCopyrightText: 2020 Wolfgang Traylor <wolfgang.traylor@senckenberg.de>

SPDX-License-Identifier: CC-BY-4.0
-->

\brief Scientific background of the Modular Megafauna Model.

\tableofcontents

## Introduction

This document explains the design choices for the modules and concepts in the megafauna model from a scientific rather than a programmatical angle.
It discusses the different submodels in the framework: what their assumptions are, how to use them, and how to combine them.

Because of the modular nature of the software, you can “generate” a range of models with different combinations and configurations of submodels.
Therefore this documentation page must not be seen as a “model description”, but rather as a loose discussion of available model components.
It is a *living* document that ought to be expanded when new features are introduced to the library.

The first section, \ref sec_model_goals, shall help you clarify the direction of your modeling project.
The section \ref sec_basic_model_concepts introduces the simulation framework of the Modular Megafauna Model: where it is flexible and where it is constrained.
The following sections describe conceptual elements of the \ref sec_dynamic_populations, which are currently represented by the herbivore type “cohort,” and conclude with some lessons learned from emerging model behavior.
At the end of this document, you will find a list of \ref sec_symbols_and_abbreviations and a remark on the choice of \ref sec_units_of_measurement in the model.

Some aspects of the model can only be evaluated in the context of the connected vegetation model.
For [LPJ-GUESS](http://iis4.nateko.lu.se/lpj-guess/) you will find those aspects in the megafauna doxygen page of the LPJ-GUESS repository.

## Model Goals {#sec_model_goals}
Before you start your modeling project, you should have your model **goal** defined.
The formulation of the model goal paints in broad strokes a picture of the direction you want to take.
In the next step, the definition of a model **purpose** will help you convert the goal statement into an **objective** statement.
The objective is more concrete than the goal and can be further refined into **model specifications,** which serve as reference points for evaluating the output of your model candidates and whether you have reached your goal.
For a more more detailed discussion of these terms see Overton (1990)\cite overton1990strategy.

The Modular Megafauna Model, by its flexible nature, shall serve as broad of a range of model goals as possible.
Its modularity is supposed to help in the iterative modeling process.
You can combine different submodels to create a number of different intermediate developmental models.
They can then be evaluated against your model specifications.
When composing your model, be wary of the **“complexity paradox”:** The more complex (i.e. “realistic”) your model is, the more uncertain it is and the more difficult it is to know if your tests of the model are meaningful (Oreskes 2003\cite oreskes2003role).

The only goal given by the software architecture of the Modular Megafauna Model can be stated thus:
*to simulate herbivore–vegetation dynamics over time.*
By selecting model components, the goal becomes more specific, and so the herbivore type “cohort” (see Section \ref sec_herbivore_cohorts) has the goal
*to dynamically simulate herbivore population densities as they emerge from basic mechanistic processes.*
Herbivore densities are not prescribed, but emerge **bottom-up.**
Based on that, project-specific model objective statements can be formulated, for instance, “to simulate bison populations in the North American prairie for the *purpose* of estimating potential pre-European carrying capacity.”
The model components must then be selected or newly implemented to meet that objective.
In North America, for example snow cover may play a role and should be included in the model.
If the objective were to simulate wildebeest in the Serengeti, snow wouldn’t play a role and should be excluded from the model.

Therefore, before you start working with the Modular Megafauna Model, gain as much clarity as possible about your model goals, purposes, objectives, and specifications.
You might even want to consider a **preregistration,** see for example Nosek et al. (2018)\cite nosek2018preregistration, Lee et al. (2019)\cite lee2019robust, and Dirnagl (2020)\cite dirnagl2020preregistration.

## Basic Model Concepts {#sec_basic_model_concepts}

The world of the megafauna model is comprised of **simulation units.**
Each such unit consists of a **habitat** and the herbivore **populations** inhabiting it.
The habitat must be implemented by the outside vegetation model.
Output can be spatially aggregated over several habitats by assigning multiple habitats to the same **aggregation unit.**
Again, this is done by the habitat implementation of the vegetation model.
What kind of **herbivores** populate the world can be defined by the user.
**Cohorts** are the first herbivore type implemented.

For the megafauna model, all habitats are of equal size and spatially inexplicit.
Spatial meaning must be given by the outside vegetation model; for example LPJ-GUESS maps each grid cell in the longitude/latitude raster to one aggregation unit.
Currently, there is no interaction between habitats, that is herbivores cannot move from one to the other.
Each habitat can be thought of as a little homogeneous capsulated world of forage and herbivores.

\image html images/model_entities.svg "Basic model entities in the Modular Megafauna Model. Aggregation units (disjunct sets of simulation units) are only needed to aggregate output spatially. Cohorts are only one possible implementation of herbivores."

The simulations run with **daily time steps.** The predecessor model by Adrian Pachzelt \cite pachzelt2013coupling operated on a monthly schedule and was thus much faster.
However, the Modular Megafauna Model should be applicable on different spatial and temporal scales.
LPJ-GUESS simulates vegetation in a daily schedule, and so naturally the attached herbivores should be treated the
same.
Moreover, there has been no formal analysis how much a coarser temporal resolution affects the model outcome.
So it seemed better to air on the side of a finer resolution.

The vegetation model makes plant biomass available to herbivores as dry matter forage mass per area.
The megafauna model works with a hard-coded set of **forage types,** for example “grass.”
It is up to the vegetation model to map its own representation of edible biomass to forage types, for instance by aggregating all graminoids, herbs, and forbs to one “grass” forage mass.
The advantage of this approach is that the herbivore model is highly decoupled from the vegetation model.
Placing forage types as hard-coded entities at the core of the herbivore model makes it easy to design herbivore diet preferences, foraging, and digestion around them.
New forage types can be implemented as necessary.

## Dynamic Populations {#sec_dynamic_populations}

### Herbivore Cohorts {#sec_herbivore_cohorts}

A cohort represents all herbivores born in one year.
The state variables of a cohort represent an average over all individual animals within the cohort.
The main reason to not simulate individuals is to save computational resources.

These are the state variables for each herbivore object:

- Age
- Sex
- Current energy need
- Fat mass

Each herbivore cohort is assigned a **herbivore functional type (HFT).**
The HFT can be interpreted as representing a species, a guild, or a trophic level.
An HFT is simply a user-defined, constant set of parameters defining physiology, life history, and everything else.
Each **cohort population** contains all cohorts of one HFT in a particular habitat.
Therefore, the maximum number of cohorts within one population is given by the HFT life span in years times two, for the two sexes.

Offspring of large herbivores usually shows an even sex ratio.
Most model processes don’t differentiate between males and females, only body size and age of maturity have sex-specific parameters.
During the model design, it seemed advisable to at least set the basis for gender differentiation because some large herbivores do show pronounced sexual dimorphism not only in size (e.g. bison or proboscideans) but also in diet and behavior (e.g. elephants, Shannon et al. (2013)\cite shannon2013diet, and steppe bison, Guthrie (1990)\cite guthrie1990frozen).

Each year, the newborn animals of all reproductive cohorts of one HFT in one habitat are combined into one new cohort.
This effectively eliminates the connection between parents and offspring.
Therefore it is not possible to implement an exchange of energy through lactation between parents and young directly.

### Body Mass and Composition {#sec_body_mass_and_composition}

The user-defined **live body mass** of simulated herbivores is the sum of blood, gut contents (ingesta), structural (fat-free) mass and deposited body fat.
It is very important to realize that the **body fat** that the model works with is the fraction of fat in the empty body.
The **empty body mass** is the live body mass minus blood, ingesta, hair, and antlers/horns.
The body fat is total lipid content, which is also known as ether extract, free lipid content, or crude fat (Hyvönen, 1996 \cite hyvonen1996approach).
This is different from the mass of suet and organ fat because the fat tissue also contains water.
The term **lean body mass** is live body mass minus all fat mass (compare quote from Blaxter, 1989 \cite blaxter1989energy in Section \ref sec_fat_as_energy_storage).
Note that lean body mass includes digesta, hair, antlers, and blood.

\image html images/body_composition.svg "Body composition in the Modular Megafauna Model."

Calder (1996, p. 14) \cite calder1996function discusses the question of variable ingesta load:

> Should the total body mass used in allometry include gut contents, a major
> source of variability (but representing mass that the animal must be designed
> to support and carry), or should gut contents be subtracted from live mass?
> Disallowing gut contents is not practical in studies wherein the animals are
> not, or should not be sacrificed.

In this line of argument, the gut contents are always included the body mass values given in the megafauna model.
It is designed for large herbivores, in particular extinct ones, and their body mass is most commonly given as a total live weight.
Live body mass is easy to measure.
Most allometric regressions are based on it.

Technically, the empty body mass is different from the “ingesta-free mass” in the literature because ingesta-free mass usually includes hair.
For less furry animals the two can be considered approximately equal, though.

### Fat as Energy Storage {#sec_fat_as_energy_storage}

The variable amount of body fat, which serves as energy reserves, is a critical component of the herbivore simulations.
As Blaxter (1989, p. 51) \cite blaxter1989energy explains, the ingesta-free animal body can be viewed as composed of fat and fat-free mass:

> Schematically the body can be regarded as consisting of two components – fat
> and non-fat. […] The non-fat material consists of water, the minerals of bone
> and soft tissue, carbohydrate, nitrogen-containing compounds, and, in the
> living animal, the contents of the digestive tract. The non-fat component of
> the body is usually referred to as the *lean body mass* or *fat-free mass.*
> Many studies, commencing with those of Murray (1922) and embracing a wide
> range of adult species, have shown that the chemical composition of the
> fat-free body is approximately constant. The wide range of composition of
> animals is largely, but not entirely, due to variation in the proportion of
> fat.

When defining the fractional body fat parameter for an herbivore, you should not rely on measurements of weight loss of starving or fattening animals.
In such data it is difficult to disentangle the contributions of changing fat mass, gut contents, water content, and fat-free mass (e.g. Reimers et al., 1982 \cite reimers1982body).

The model ignores the contribution of catabolizing protein altogether.
This is a simplification as several studies have shown that the contribution of mobilized protein can play a considerable role in meeting energy requirements.
Reimers et al. (1982, pp. 1813, 1819) \cite reimers1982body observe that reindeer lost about 31% of their body protein during winter.
Torbit et al. (1988) \cite torbit1988calibration calculate and discuss the energetic value of catabolizing protein.
Parker et al. (1993) \cite parker1993seasonal note that Sitka black-tailed deer lost 10–15% of their protein reserves during winter.
Nontheless fat is unquestionably by far the most important energy reserve.
Fluctuating body fat is convenient to model as an energy pool that can be filled to a maximum and emptied to zero.
The interactions for protein synthesis and depletion are far more complex and less studied.

Any forage energy that is ingested beyond maintenance needs is converted to body fat.
Section \ref sec_energy_content details the efficiency of fat anabolism.
When forage intake is not enough to meet energy needs, fat reserves are catabolized, and the energy is directly available as net energy to balance any energy “debts” on a daily basis.
Illius & O’Connor (2000) \cite illius2000resource assumed an efficiency of 100% to mobilize fat reserves in cattle when they directly converted the combustion (gross) energy of fat tissue to net energy.
Armstrong & Robertson (2000) \cite armstrong2000energetics use a factor of 80% in their sheep model, citing a 1990 publication by the Australian Standing Committee on Agriculture (SCA) \cite corbett1990feeding.
The megafauna model allows the catabolism efficiency to be specified by the user.
This efficiency factor is multiplied with the fat gross energy to derive the net energy gain (MJ) from burning one kg of body fat.

### Ontogenetic Growth

The growth curve is currently linear.
The body mass of a cohort that hasn’t reached physical maturity yet is calculated as a linear interpolation between the user-defined neonate body mass and adult body mass.
However, the growth curve should be sigmoid, compare Price (1985, pp. 187–190)\cite price1985growth and this quote by Blaxter (1989)\cite blaxter1989energy, p. 242f:

> Growth in weight is characteristically sigmoid; it accelerates during a short
> initial period and then declines until, as maturity approaches, it approaches
> zero. A large number of different functions have been used to describe this
> relationship between weight and time. Virtually all of them state that dW/dt,
> the rate of change in weight with time (or *rate of growth*) is a function of
> weight at the time that dW/dt is measured. Such functions include the
> logistic equation, the Gompertz function and the Bertalanffy function. These
> are derived algebraically in most textbooks of biomathematics (see Causto
> 1977). They were generalised by F.J. Richards (1959) and have been reviewed
> and critically analysed by Parts (1982).

### Energy Budget

\image html images/energy_budget.svg "Model of energy budget for a ruminant or hindgut fermenter. Modified after Minson (1990), Fig. 5.1."

#### Energy Expenditure

The allometric scaling of metabolic rate in the megafauna model does not differentiate between inter- and intraspecific scaling.
This model assumption might need to be re-evaluated in the future.
As Makarieva et al. \cite makarieva2009comment point out:

> However, it has repeatedly been observed that young animals have elevated
> metabolic rates compared with what is predicted for their body mass from
> interspecific scaling (5–8).

Compare also Glazier (2005) \cite glazier2005beyond for a discussion on intraspecific metabolic scaling.

#### Energy Content of Forage {#sec_energy_content}

The model for energy content in herbivore forage presented here is based on the partitioning of metabolizable energy.
A historical overview of the model framework is given by Ferrell & Oltjen (2008) \cite ferrell2008asas.
Its conceptual shortcomings and difficulties in practical methodology are summarized by Birkett & de Lange (2001) \cite birkett2001limitations.

The diagram on the energy budget shows how energy from the forage is used by an herbivore: **Gross energy** (\f$GE\f$) is the heat that could be produced from complete combustion of the feedstuff.
From that, the part which is not excreted in feces is the **digestible energy** (\f$DE\f$).
Some proportion of it is then lost to urine and gas production, but the rest is **metabolizable energy** (\f$ME\f$).
After deducing now the losses due to heat increment, the remaining **net energy** (\f$NE\f$) is effectively utilizable for all physiological processes.

Gross energy depends only on the physical properties of the forage and measured in a combustion chamber.
It is therefore independent of the animal.
McDonald et al. (2010, p. 259)\cite mcdonald2010animal provide an overview of gross energy in different feedstuffs for livestock: It typically ranges between 18 to 20 MJ/kgDM.
The measurements by Golley (1961) \cite golley1961energy suggest that there is some seasonal variation in gross energy of leaves.
However, the model assumes a constant value.

The proportional **dry-matter digestibility** (\f$DMD\f$) of the forage is a central variable in the model.
It measures the fraction of the gross energy that is usable by the animal.
The rest gets excreted in the feces because it is undigestible fiber: protected cellulose and hemicellulose, silica, and cutin.
Agricultural research has shown that the digestibility is closely correlated with metabolizable energy and net energy (Minson, 1990, p. 7 \cite minson1990forage).
Therefore digestibility is modeled as the *one* indicator for forage quality, i.e. forage energy density, and must be given by the vegetation model.
This neglects any other effects on the digestibility, like interactions of different forages or effects of the individual animal on the digestibility.
Digestibility is best measured *in vivo* in the rumen of a living ruminant, but there exist various indirect methods with reliable conversions.
For an overview see Minson (1990) \cite minson1990forage and McDonald (2010) \cite mcdonald2010animal.
In order to precisely define the model interface: Formulas in the megafauna model assume *in vivo* digestibility of livestock ruminants.

Research has focused mainly on the digestion of ruminant livestock, and so the megafauna model works primarily with the well-established formulas for ruminants.
Digestibility is also defined for ruminants.
To account for the less efficient digestion of hindgut fermenters, the user can define a **digestibility multiplier** to convert ruminant digestibility to hindgut digestibility.
This approach is taken from by and Pachzelt et al. (2015) \cite pachzelt2015potential, who cite Illius & Gordon (1992) \cite illius1992modelling

Ruminants typically lose a relatively constant fraction of about 19% of digestible energy in urine and methane (López et al. 2000 \cite lopez2000prediction, McDonald et al. 2010 \cite mcdonald2010animal, p. 258).
The difference between cattle and sheep is very small here (McDonald et al. 2010, p. 260).
McDonald et al. (2010, p. 258) specify that 11–13 percent of digestible energy is lost as methane.
The 19% loss to urine and gases is often expressed as the ratio of metabolizable energy to digestible energy, \f$ME/DE=0.81\f$.
This ratio is also known as the **metabolizable energy coefficient** (e.g. in Robbins, 1983\cite robbins1983wildlife).With a gross energy of about 19 MJ/kg, metabolizable energy in the digestible fraction of the forage is then about 15–16 MJ/kg.
Various herbivore models work with these numbers, for instance: Givens et al. (1989)\cite givens1989digestibility, Illius and Gorden (1991)\cite illius1991prediction, Parker et al. (1991) \cite parker1996foraging, Illius and Gordon (1999)\cite illius1999scaling, Smallegange and Brinsting (2002)\cite smallegange2002food.

\warning
Some publication, like Minson (1990)\cite minson1990forage, use the term “metabolizability of energy” or “metabolizable energy coefficient” to refer to the \f$ME/GE\f$ ratio: the metabolizable fraction of the *gross* energy.
This includes fecal losses and has the dry-matter digestibility already calculated in.
However, the modular megafauna model works with explicit digestibility values and the ME/DE ratio.
You could divide the \f$ME/GE\f$ ratio by the fractional digestibility to get \f$ME/DE\f$.

A unitless **net energy coefficient** (\f$k\f$) defines the efficiency of using the metabolizable energy for meeting maintenance energy needs, i.e. for converting metabolizable energy content to **net energy** content (\f$NE\f$) of the forage.
(In Robbins (1983)\cite robbins1983wildlife it is called *NEC*.)
Many livestock models differentiate between different k values to reflect different conversion efficiencies: for meeting maintenance needs (\f$k_m\f$), for growth and fattening (\f$k_f\f$), and for lactation (\f$k_l\f$) (Blaxter 1989, p. 254ff\cite blaxter1989energy; Minson 1990\cite minson1990forage, p. 151).

In the Modular Megafauna Model, the energy budget calculates with the “currency” net energy, which broadly represents the available oxidizable metabolic fuels—glucose and fatty acid.
Basal and field metabolic rate and other energy expenditures are directly “paid” with net energy.
Therefore the efficiency factor \f$k_m\f$ is used to convert from metabolizable energy to net energy.
Body fat is anabolized from metabolizable forage energy with the efficiency factor \f$k_f\f$.

\image html images/retention_over_intake.svg "Model of energy retention in an herbivore. km and kf denote the slope of the line, i.e. the efficiency of utilizing metabolizable energy. When fed maintenance requirements, the animal will neither gain nor lose weight. Below that point it will starve (i.e. catabolize reserves) and above it will build reserves (i.e. anabolize fat). After McDonald et al. (2010), Fig. 11.5."

Feeding trials have shown that the net energy coefficient can linearly depend on the metabolizable energy content of the forage (Robbins, 1983, p. 296f; Minson, 1990, p. 93, 155).
However, this effect seems to be mostly related to very high levels of feeding and by pelleting the feed.
In this model, the energy coefficients \f$k_m\f$ and \f$k_f\f$ are assumed to be constant.

\remark
Internally the model converts first from metabolizable energy to net energy to pay energy expenditures.
If there is excess net energy, this gets converted *afterwards* to body fat.
The amount of net energy required to build up one kilogram of body fat is given by the product of fat gross energy content, \f$k_m\f$ and \f$k_f^{-1}\f$.
Note that Illius & O’Connor (2000) \cite illius2000resource probably took the same approach when they specify an anabolism coefficient of 54.6 MJ/kg, citing Blaxter (1989) \cite blaxter1989energy.
Namely, 54.6 MJ/kg is the product of 39 MJ/kg, \f$k_m=0.70\f$, and the inverse of \f$k_f=0.50\f$; probably Illius & O’Connor (2000) took the latter two figures from Table 12.1 on page 259 in Blaxter (1989) for oxen on an “average diet.”

In summary:
Net energy content, \f$NE\f$ in MJ/kgDM, depends on variable dry-matter digestibility, \f$DMD\f$, as the key variable.
Gross energy content, \f$GE\f$, is user-specified for each forage type.
Only the digestible fraction of the gross energy in dry matter is counted as digestible energy, \f$DE\f$.
How much metabolizable energy can be extracted from the digested part of the forage is species-specific and defined by the user as the metabolizable energy coefficient or \f$ME/DE\f$ ratio.
A user-defined factor, \f$k_m\f$, defines how efficient the metabolizable energy is used to meet net energy needs for maintenance and other activities.
The factor \f$k_f\f$ denotes the efficiency for converting from \f$ME\f$ to body fat (anabolism).
The net energy content is given by:

\f[
NE = ME * k_m = DE * \frac{ME}{DE} * k_m = GE * DMD * \frac{ME}{DE} * k_m
\f]

#### Thermoregulation by Conductance {#sec_thermoregulation}

This model of thermoregulation is often called the **Scholander-Irving model** and was published in two seminal papers in 1950: Scholander et al. (1950a) \cite scholander1950adaptation and Scholander et al. (1950b) \cite scholander1950heat.
The more detailed implementation is taken from Peters (1983) \cite peters1983ecological.

Homeothermic animals have extra energy costs to maintain their body core temperature.
Through basal metabolism and other ways of energy burning, heat is already passively created.
Thermoregulatory costs arise when the ambient temperature drops below the *lower critical temperature*: the passive heat from thermoneutral metabolism is not counterbalance heat loss to the environment.
The rate of heat loss depends on the *thermal conductance* of the whole animal (energy flow per temperature difference), which in turn depends on the *thermal conductivity* (energy flow per temperature difference and per thickness) of fur and skin and the body surface.
Conductance is the inverse of resistance or insulation, and conductivity is the inverse of resistivity.

- \f$T_{crit}\f$: Lower critical temperature [°C].
- \f$T_{core}\f$: Body core temperature [°C].
- \f$T_{air}\f$: Ambient air temperature [°C].
- \f$E_{neu}\f$: Thermoneutral metabolic rate [MJ/ind/day]
- \f$C\f$: Whole-body thermal conductance [W/ind].
- \f$\Phi\f$: Heat loss [MJ/ind/day]

\f[
T_{crit} = T_{core} - \frac{E_{neu}}{C}
\f]

\f[
\Phi = C * max(T_{crit} - T_{air}, 0)
\f]

\image html thermoregulation.svg "Schematic description of the effects of external temperature on the metabolic rate in homeotherms. After Peters 1983, Fig. 5.6"

\note
In its current form, the model only considers costs when temperatures are too low.
Overheating effects are not implemented since the model was developed with the focus on Arctic megafauna.

The critical parameter for thermoregulatory expenditure is the **whole-body conductance:** the rate of heat flow per difference between core and air temperature (W/°C).
The conductance can be approximated from the average conductivity and the body surface.
Conductivity is the inverse of insulation: it is the heat flow per temperature difference per area.

Body surface in m² scales roughly as \f$0.09*BM^{0.66}\f$ (Hudson & White 1985\cite hudson1985bioenergetics).

#### Foraging {#sec_foraging}

An herbivore’s daily dry matter intake can be limited by any of a number of factors, as illustrated by the following figure:

@startuml "Levels of herbivore intake constraints. What and how much of the available forage an herbivore ingests is limited by a cascade of internal and external factors."
	!include diagrams.iuml!intake_limit_levels
@enduml

The **functional response** is the intake rate as a function of available forage (Holling 1959a\cite holling1959components) and corresponds to the “foraging limit” in the megafauna model.
Grazers are generally thought to have a Type II functional response (Owen-Smith 2002\cite owensmith2002metaphysiological): their intake rate quickly increases with increasing food abundance towards an asymptotic maximum.
This maximum is the rate at which an herbivore could *theoretically* ingest forage if there were no constraints of forage abundance, digestive capacity, or metabolic requirements.

The Type II functional response is commonly expressed as a hyperbolically saturating function (Holling 1959b\cite holling1959some), which is also known as [“Michaelis–Menten” function](https://en.wikipedia.org/wiki/Michaelis%E2%80%93Menten_kinetics) in other contexts.
The critical parameter of this function is the forage density at which the intake rate reaches half of its maximum.
Owen-Smith (2002)\cite owensmith2002metaphysiological calls it \f$v_{1/2}\f$; Illius and O’Connor (2000) call it \f$\beta\f$.

\note
Illius and O’Connor (2000)\cite illius2000resource and subsequent similar models (e.g. Pachzelt et al. 2013\cite pachzelt2013coupling) use an empirically derived half-maximum intake density, \f$\beta\f$, to constrain intake.
However, they set the *daily digestive limit* as the asymptotic maximum.
This can successfully create a density dependence effect: when herbivore densities rise and forage becomes scarce, the intake rate *gradually* decreases.
However, it is not congruent with the original empirical measurement of \f$\beta\f$, which assumes a *short-term* intake rate irrespective of digestive capacity.

\image html images/functional_response_types.svg "The different types of functional responses. This image is in the Public Domain."

Following a seminal publication by Spalinger and Hobbs (1992)\cite spalinger1992mechanisms, a lot of work has been done to model the functional response of grazers mechanistically (e.g.
Illius and Fitzgibbon 1994\cite illius1994costs;
Bradbury et al. 1996\cite bradbury1996relationship;
Fortin et al. 2002\cite fortin2002temporal;
Hobbs et al. 2003\cite hobbs2003challenges;
Fortin et al. 2004\cite fortin2004multitasking;
Robinson and Merrill 2012\cite robinson2012influence).
In the short term, grazers may be limited by encounter rate (moving between forage patches) or handling rate (cropping and chewing), which means they might not be able to ingest the amount of forage they desire even though there is enough forage available.
However, in the long term (i.e. over days and months), the intake of grazers is most likely to be digestion-limited.
That is because they can compensate encounter and handling limitation by increasing their daily foraging time.
In other words, the functional response curve increases so sharply with even low forage density that it does not play a major role in the long term.

### Mortality {#sec_mortality}

Cohorts are deleted when they have reached their **life span.**
The population numbers are gradually decreased with a user-specified **annual background mortality** and by **starvation.**

#### Death of Starvation

In the process of starvation, different fat depots are mobilized in a typical sequence: rump fat, subcutaneous fat, visceral fat, and, finally, marrow fat (Hanks, 2004/1981 \cite hanks2004characterization).
When the energy reserves of an animal are exhausted, it will die of starvation.
Body fat, i.e. lipid in the ingesta-free body, is then zero.
Different studies have found that the carcasses of large herbivores that have starved to death contain virtually no body fat anymore (Reimers et al., 1982 \cite reimers1982body; Depperschmidt et al., 1987 \cite depperschmidt1987body), but chemical analysis of fat content in carcass samples can be imprecise (Depperschmidt et al., 1987).

### Population Dynamics {#sec_population_dynamics}

#### Population Stability
Several qualities of the model result in a propensity towards extreme instability of the simulated herbivore populations, with exponential irruptions followed by sudden crashes (sometimes to extinction).
Here are some of the reasons for these “boom–bust cycles”:

- If there is **no regular seasonal die-off,** there is effectively **no early density dependence** effect. Herbivores just reproduce exponentially until the population crashes completely. In winter or the dry season, the fat storage should drop so far that a fraction of the population dies because then density dependence effects occur at those times when there is not enough forage in the vegetation period to completely fill all animals’ fat storage. Another controlling factor can be **variable climate:** long or harsh winters/dry seasons can also regulate the population and prevent uninterrupted exponential growth. See for example Stewart et al. (2005)\cite stewart2005densitydependent for a discussion of the roles of summer versus winter in regulating populations of large mammals.
- There is **no resource heterogeneity,** only grass. If there was more heterogeneity with low-quality and high-quality forage, populations might be more stable (Owen-Smith 2004\cite owensmith2004functional).
- There is **no movement** between habitats. That is, a crashing population has no way to escape extinction by moving somewhere else. The situation in each habitat is comparable to an island population (e.g. Klein 1968\cite klein1968reindeer).
- **Annual grass allocation** in the case of MMM being coupled with LPJ-GUESS: Since LPJ-GUESS allocates one year’s NPP at the end of the year (Dec. 31st), herbivores will have to starve until the end of the year when all forage is eaten.

For further discussion and brainstorming see also these excerpts from Wolfgang Traylor’s lab notebook on [Open Science Framework](https://osf.io/):

- [2017-09 Why do Populations Crash?](https://osf.io/df69m/)
- [2017-10 Christmas Present](https://osf.io/fj74e/)
- [2017-12 Coexistence and Stability](https://osf.io/bajku/)
- [2018-02 Closeup on a Crash](https://osf.io/5pnyg/)
- [2018-06 Fat, Stability, and Movement](https://osf.io/8swth/)

#### Minimum Density Threshold {#sec_minimum_density_threshold}
The parameter \ref Fauna::Hft::mortality_minimum_density_threshold defines at which point a dwindling population (sum of all cohorts) may be considered dead.
It is an arbitrary and artificial, but critical “tuning parameter” for realistic model performance.
Possible re-establishment only happens if all cohorts are dead within one habitat.

It is important to keep this parameter low enough for slowly breeding and long-lived animals because otherwise they may die out after establishment:
After establishment, the background mortality continually diminishes the adult cohorts, and after some years the total population (all cohorts together) my drop below the `minimum_density_threshold` before reproduction could compensate.

On the other hand, the `minimum_density_threshold` should not be set *too* low as this would result in extremely thin “ghost” populations that are effectively preventing re-establishment.

#### Species Coexistence {#sec_coexistence}

The classical competitive exclusion principle predicts that no two species can coexist in the long term if they each solely depend on one shared resource (Hardin 1960\cite hardin1960competitive).
One species will inevitably outcompete the other one.
Though there are indeed ecological mechanisms that can facilitate coexistence with a shared resource (Chesson 2000\cite chesson2000mechanisms), the parameter space for this to happen in a model is usually very narrow (e.g. van Langevelde et al. 2008\cite vanlangevelde2008intantaneous).

In order to simply avoid competition among different HFTs, the option \ref Fauna::Parameters::one_hft_per_habitat can be enabled: Each HFT exists on its own, without any interaction with other species.
With that option enabled, all HFTs should each be assigned to the same number of habitats.
It is the responsibility of the host application (the vegetation model) to ensure that the number of habitats is an integer multiple of the HFT count.

## Symbols and Abbreviations {#sec_symbols_and_abbreviations}

- \f$bf\f$  = Current body fat as fraction of lipids per empty body [frac.]
- \f$BM\f$  = Live body mass [kg/ind]
- \f$DM\f$  = Dry matter
- \f$DMD\f$ = Dry-matter digestibility [frac.]
- \f$eb\f$  = Empty body fraction [frac.]
- \f$FM\f$  = Fat mass of an individual [kg/ind]
- \f$GE\f$  = Gross energy, also known as heat of combustion or calorific value [MJ/kgDM]
- \f$k\f$   = Net energy coefficient, efficiency of converting ME to usable energy
    - Subscript \f$_f\f$ = conversion to fat gross energy
    - Subscript \f$_m\f$ = conversion to NE, to meet needs of maintenance and field metabolic rate
    - Subscript \f$_p\f$ = conversion to protein gross energy
- \f$ME\f$  = Metabolizable energy content in forage [MJ/kgDM]
- \f$MRT\f$ = Mean retention time [hours]
- \f$NE\f$  = Net energy content in forage, usable for meeting energy requirements [MJ/kgDM]
- \f$SM\f$  = Structural mass [kg/ind]
- General subscripts:
    - \f$_{ad}\f$    = Adult
    - \f$_{birth}\f$ = At birth/for neonates
    - \f$_{max}\f$   = Maximum (body fat, fat mass, reproduction rate, …)

Note that “mass” and “weight” are used interchangeably.

## Units of Measurement {#sec_units_of_measurement}
- All forage values (e.g. available grass biomass, consumed forage) are *dry matter mass* in kilograms (`DMkg`).
- Any forage per area (e.g. forage in a habitat) is `kgDM/km²`.
- Herbivore-related mass values (e.g. body mass, fat mass) are also `kg`, but live mass (see \ref sec_body_mass_and_composition).
- Population densities of herbivores are either in `kg/km²` or `ind/km²` (“ind” = “individuals”).
- Digestibility values are interpreted as in-vitro digestibility (see \ref sec_energy_content).

\remark
The units of measurement were primarily chosen in a way to yield numbers broadly around zero for calculation.
Floating point operations are most precise then, and the values can be printed in a text output table.
When post-processing the output, you can convert to your units of choice, e.g. `ind/ha` instead of `ind/km²`.

-------------------------------------------------

\copyright <a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/80x15.png" /></a> This software documentation is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
\author Wolfgang Traylor, Senckenberg BiK-F
\date 2019
# Software Design {#page_design}
<!--
SPDX-FileCopyrightText: 2020 Wolfgang Traylor <wolfgang.traylor@senckenberg.de>

SPDX-License-Identifier: CC-BY-4.0
-->

\brief Notes on the software design of the herbivore model from a programmer’s perspective.

\tableofcontents

\todo Encourage to do refactoring of the architecture if it supports modularity and flexibility.

## Overview {#sec_design_overview}

The megafauna model aims to apply principles of object-oriented programming as much as possible (see the page [Object-oriented Programming](\ref page_object_orientation)).
Its architecture is modular and extensible.
Each part can be tested in unit tests.

The following UML diagram shows through which interfaces the megafauna model interacts with other components:

@startuml "Component diagram of the basic interactions of the megafauna model."
	!include diagrams.iuml!basic_components
@enduml

The Modular Megafauna Model was originally within the code base of the LPJ-GUESS ecosystem model.
In 2019 it was separated into its own library for the following reasons:

- Model is **reusable** with other vegetation models.
- Model can be better **tested** in isolation and with continuous integration (CI).
- More **freedom** to structure the repository and its documentation.
- The model can be **licensed** and **distributed** independently from LPJ-GUESS.
- **Memory leaks** can be found more easily since they are not mingled with LPJ-GUESS.

### Simulation Design

The basic simulation design is simple:

- **Herbivores** (\ref Fauna::HerbivoreInterface) are independent entities that interact with their environment.
- Each herbivore lives in a **habitat** (\ref Fauna::Habitat), grouped in **populations** (\ref Fauna::PopulationInterface).
  In LPJ-GUESS this would correspond to a plant `Individual` growing in a `Patch`.
- **Fauna::World** is the framework class running the simulation.

@startuml "Most important classes in the megafauna model."
	!include diagrams.iuml!important_classes
@enduml

All interactions between herbivores and their environment happen through \ref Fauna::Habitat.
The herbivores don’t feed themselves and don’t have any direct connection to the habitat.
Instead, they are assigned their forage.
With each simulated day (\ref Fauna::HerbivoreInterface::simulate_day()), it is calculated, how much forage they would like to consume (\ref Fauna::HerbivoreInterface::get_forage_demands()).
Based on all forage demands, the forage is distributed with a user-specified forage distribution algorithm (\ref Fauna::DistributeForage).
Then they can eat their portion (\ref Fauna::HerbivoreInterface::eat()).
This approach follows the [Inversion of Control Principle](\ref sec_inversion_of_control).

Similarly, the Habitat does not interact with the herbivores either.
It does not even *know* about the herbivore populations, as it is capsuled in \ref Fauna::SimulationUnit.

## Forage Classes {#sec_design_forage_classes}

The model is designed to make implementation of multiple types of forage (like grass, browse, moss, etc.) easy.
Each forage type is listed in \ref Fauna::ForageType.
The global constant \ref Fauna::FORAGE_TYPES holds all entries of this enum so that it’s easy to iterate over them.

The template class \ref Fauna::ForageValues serves as a multi-purpose container for any forage-specific values.
Many arithmetic operators are defined to perform calculations over all forage types at once.
For any specific use of the class, a semantic `typedef` is defined, e.g. \ref Fauna::ForageMass or \ref Fauna::Digestibility.
This helps to directly see in the code what a variable contains.

Forage types need to have specific model properties.
Grass, for instance, has the property *sward density,* which would not make sense for leaves of trees.
Therefore, a second set of forage classes is defined with one class for each forage type.
All these classes inherit from \ref Fauna::ForageBase.

Any forage properties are defined by the habitat implementation in \ref Fauna::Habitat::get_available_forage().
They can be used for example in algorithms of:

- forage distribution (\ref Fauna::DistributeForage),
- diet composition (\ref Fauna::GetForageDemands::get_diet_composition),
- digestion limits (\ref Fauna::GetForageDemands::get_max_digestion), or
- foraging limits (\ref Fauna::GetForageDemands::get_max_foraging).

@startuml "Forage classes in the megafauna model."
	!include diagrams.iuml!forage_classes
@enduml

\see \ref sec_new_forage_type

## The Herbivore {#sec_design_the_herbivore}

The simulation framework can operate with any class that implements \ref Fauna::HerbivoreInterface (compare \ref sec_liskov_substitution).
Which class to choose is defined by the instruction file parameter \ref Fauna::Parameters::herbivore_type.

Currently, only one herbivore class is implemented: \ref Fauna::HerbivoreCohort.
The herbivore model performs calculations generally *per area* and not per individual.
The area size of a habitat is undefined.

@startuml "Class diagram of the default herbivore class: Fauna::HerbivoreCohort."
	!include diagrams.iuml!herbivore_classes
@enduml

\see \ref sec_new_herbivore_class

### HerbivoreBase {#sec_herbivorebase}
The herbivore class itself can be seen as a mere framework (compare \ref sec_inversion_of_control) that integrates various components:

- Each herbivore is of one **Herbivore Functional Type (HFT):** \ref Fauna::Hft
- The herbivore’s own **energy budget**: \ref Fauna::FatmassEnergyBudget.
- Its **energy needs**, defined by \ref Fauna::Hft::expenditure_components.
The herbivore object is self-responsible to call the implementation of the given expenditure models.
(A strategy pattern would not work here as different expenditure models need to know different variables.)
- How much the herbivore **is able to digest** is limited by a single algorithm defined in \ref Fauna::Hft::digestion_limit.
- How much the herbivore **is able to forage** can be constrained by various factors which are defined as a set of \ref Fauna::Hft::foraging_limits.
- The **diet composition** (i.e. feeding preferences in a scenario with multiple forage types) is controlled by a the model selected in \ref Fauna::Hft::foraging_diet_composer, whose implementation should be called in \ref Fauna::GetForageDemands::get_diet_composition().
- How much **net energy** the herbivore is able to gain from feeding on forage is calculated by a selected net energy model: \ref Fauna::NetEnergyModel.
(given by [constructor injection](\ref sec_inversion_of_control)).
- **Death** of herbivores is controlled by a set of \ref Fauna::Hft::mortality_factors.
For a cohort that means that the density is proportionally reduced.
The corresponding population objects will release dead herbivore objects automatically.

@startuml "Model components around Fauna::HerbivoreBase. Each component is selected by an HFT enum parameter by the user through the instruction file. The herbivore class then creates/calls the appropriate classes and functions."
	!include diagrams.iuml!herbivorebase_compartments
@enduml

### Populations {#sec_populations}
Each herbivore class needs a specific population class, implementing \ref Fauna::PopulationInterface, which manages a list of class instances of the same HFT.
Each habitat (\ref Fauna::Habitat) is populated by herbivores.
The class \ref Fauna::SimulationUnit contains a habitat and the herbivore populations (\ref Fauna::PopulationList).

A herbivore population instantiates new herbivore objects in the function \ref Fauna::PopulationInterface::establish().
For cohort herbivores, there is a simple helper class to construct new objects: \ref Fauna::CreateHerbivoreCohort.
The `establish()` function is called by the simulation framework (\ref Fauna::World).
In this design, the framework is only responsible for triggering the spawning of herbivores.
How the reproduce and die is managed by the herbivore class itself, and the corresponding population and creator class.

@startuml "Herbivore population classes."
	!include diagrams.iuml!population_classes
@enduml

## Error Handling {#sec_design_error_handling}

### Exceptions {#sec_design_exceptions}
The library uses the C++ standard library exceptions defined in `<stdexcept>`.
All exceptions are derived from `std::exception`:

@startuml "Standard library exceptions used in the megafauna library."
	!include diagrams.iuml!exception_classes
@enduml

Any function that potentially *creates* an exception declares that in its doxygen description.

\warning Beware that any function—unless documented otherwise—will not catch exceptions from calls to other functions.
Therefore, even if a function does not announce a potential exception throw in its documentation, exceptions created in other functions can arise.

Exceptions are used…:
- …to check if the TOML instruction file is correct.
- …to check if parameters in public methods are valid.
- …to check the validity of variables coming from outside of the herbivory module where there are no contracts defined and ensured.

You throw an exception (in this case class `std::invalid_argument`) like this:

    ```cpp
    if (/*error occurs/*)
        throw std::invalid_argument("My error message");
    ```

Each class makes no assumptions about the simulation framework (e.g. that parameters have been checked), but solely relies on the class contracts in the code documentation.

The megafauna library will never exit the program or print messages to STDERR or STDOUT.
No exceptions arising in the library will be handled.
The host program is responsible to handle exceptions from the library and to inform the user and halt the program.

Exception messages in the megafauna library should start with the fully qualified function name (including namespace) of the function that is creating the exception.

\remark
If you debug with [`gdb`](https://www.gnu.org/software/gdb) and want to backtrace an exception, use the command `catch throw`.
That forces gdb to stop at an exception, and then you can use the command `backtrace` to see the function stack.

### Assertions {#sec_design_assertions}
At appropriate places, `assert()` is called (defined in the standard library header `<cassert>`/`assert.h`).
`assert()` calls are only expanded by the compiler if compilation happens for DEBUG mode; in RELEASE, they are completely ignored.

Assertions are used…:
- …within non-public methods to check within-class functionality.
- …to verify the result of an algorithm within a function.
- …in code regions that might be expanded later: An assert call serves as a reminder for the developer to implement all necessary dependencies.

## Parameters {#sec_design_parameters}

All user-defined simulation parameters are contained in the two classes \ref Fauna::Hft and \ref Fauna::Parameters.
All parameters must be constant within one simulation run.
Since some classes work with pointers to the classes \ref Fauna::Hft and \ref Fauna::Parameters, all objects of these classes must not be moved in memory.
For that reason, `std::shared_ptr` is used in \ref Fauna::HftList.

The host program only passes the path to the [TOML](https://github.com/toml-lang/toml) instruction file to the class \ref Fauna::World.
The parameters are parsed by the megafauna library independently, using [cpptoml](https://github.com/skystrife/cpptoml).
This is done by the class \ref Fauna::InsfileReader.

The TOML format is chosen because it is easy to read for humans and easy to parse with a free library.
Many people might be already familiar with similar syntax from DOS ini files.
The table arrays are particularly useful for defining an HFT set.
Alternative formats don’t have all these advantages.
The YAML standard is a bit too complex for our purpose, but would have been another good candidate.
JSON and XML are not so easy for humans to read and edit.

Following the [Inversion of Control](\ref sec_inversion_of_control) principle, as few classes as possible have direct access to the classes that hold the parameters.
These classes play the role of the “framework”: They call any client classes _only_ with the very necessary parameters instead of the complete \ref Fauna::Hft or \ref Fauna::Parameters objects.
The following diagram gives an overview:

@startuml "Classes that have direct access to the parameter-holding classes Fauna::Hft and Fauna::Parameters."
	!include diagrams.iuml!parameters_access
@enduml

## Output {#sec_design_output}

### Output Classes {#sec_design_output_classes}

Output classes within the herbivory module are collected in the namespace \ref Fauna::Output.
- The two structs \ref Fauna::Output::HabitatData and \ref Fauna::Output::HerbivoreData are simple data containers.
- The struct \ref Fauna::Output::CombinedData represents one datapoint (‘tupel’/‘observation’) of all output variables in space and time.

@startuml "Output classes of the herbivory module."
	!include diagrams.iuml!output_classes
@enduml

There are three levels of data aggregation:

1) Each day in \ref Fauna::World::simulate_day(), a new set of output data (\ref Fauna::Output::CombinedData) is created for each simulation unit.
For this, the habitat data is taken as is, but the herbivore data is aggregated per HFT (see \ref Fauna::Output::HerbivoreData::create_datapoint()).
This level of aggregation is **spatial within one habitat**.
Sums and averages are calculated.
Any variables *per habitat* or *per area* are summed, for instance herbivore densities.
Variables *per individual* are averaged, using individual density as weight.

2) The second level of aggregation happens also in \ref Fauna::World::simulate_day().
The datapoint for that day is added to the temporal average in the \ref Fauna::SimulationUnit object.
This level of aggregation is therefore **temporal across days**.

3) The third level of aggregation takes place in \ref Fauna::Output::Aggregator.
Here, the accumulated temporal averages from the simulation units are combined in spatial aggregation units (\ref Fauna::Output::Datapoint::aggregation_unit).
This level of aggregation is therefore **spatial across habitats**.

For the latter two aggregation levels the function \ref Fauna::Output::CombinedData::merge() is used.

\note
All time-dependent variables are always **per day.**
For example, there is no such thing like *forage eaten in one year.*
This way, all variables can be aggregated using the same algorithm, whether they are time-independent (like *individual density*) or represent a time-dependent rate (like *mortality* or *eaten forage*).


#### Pros and Cons of the Output Design {#sec_output_prosandcons}

The pros of this design:

- Simplicity: There are only few classes.
- Separation of concerns: Each class (herbivores and habitats) is self-responsible for managing its own output, and the output data containers are self-responsible for aggregating their data.
- Diversity of data structures: There is no restriction in regards to data type for new member variables in the output containers (as long as they can be merged).

The cons of this design:

- Strong coupling: The output module is highly dependent on the data structure of the output containers.
- Rigidity of data containers: Ideally, the containers should be oblivious to the details of the data they hold.
- Lack of modularity: A submodule of, e.g. \ref Fauna::HerbivoreBase cannot easily deliver its own output variable.
- Cumbersome extensibility: New output variables need to be introduced in various places (see \ref sec_new_output_variable).
That is a violation of the [Open/Closed Principle](\ref sec_open_closed).
- Any variable that is specific to a submodule or interface implementation (e.g. `bodyfat` is specific to \ref Fauna::HerbivoreBase) will produce undefined values if that submodule is not active.
The user is then responsible to interpret them as invalid or disable their output.
So far, there is no check of congruency between [parameters](\ref Fauna::Parameters)/[HFT settings](\ref Fauna::Hft) and the selection of output variables in the output module.

-------------------------------------------------

\copyright <a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/80x15.png" /></a> This software documentation is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
\author Wolfgang Traylor, Senckenberg BiK-F
\date 2019
# Demo Simulator {#page_demo_simulator}
<!--
SPDX-FileCopyrightText: 2020 Wolfgang Traylor <wolfgang.traylor@senckenberg.de>

SPDX-License-Identifier: CC-BY-4.0
-->

\brief Introduction to the demo vegetation model that uses the megafauna library.

\tableofcontents

## Overview

The demo simulator shall demonstrate how to integrate the megafauna library.
Moreover, it serves as a testing framework to run the megafauna model with as little overhead as possible and in a controlled environment.

The demo simulator classes are all in the namespace \ref Fauna::Demo.
The central class running the program is \ref Fauna::Demo::Framework.
It employs \ref Fauna::World to execute the megafauna simulation.

The class \ref Fauna::Demo::SimpleHabitat implements a very basic vegetation model that can be parametrized with custom parameters in the instruction file.
Only this one kind of vegetation model is implemented.
Grass growth with a logistic growth function:

\image html images/logistic_growth.svg "Logistic growth of demo grass model."

The `SimpleHabitat` class corresponds to the LPJ-GUESS `Patch`.

Each “habitat group” can be considered a list of \ref Fauna::Demo::SimpleHabitat objects.
The “habitat group” corresponds conceptually to the LPJ-GUESS `Gridcell`.

Each aggregation unit (habitat group) comprises several habitats.
Since all habitats have the same properties and there is no stochasticity, the output from all aggregation units will look the same.

@startuml "Class diagram of the megafauna demo simulator."
	!include diagrams.iuml!demo_simulator_classes
@enduml

### Parameters

An example instruction file is provided in `examples/demo_simulation.toml`.
It is completely separate from the megafauna library instruction file,
It emulates the scenario of the metaphysiological model by Norman Owen-Smith \cite owensmith2002metaphysiological during growing season.

-------------------------------------------------

\copyright <a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/80x15.png" /></a> This software documentation is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
\author Wolfgang Traylor, Senckenberg BiK-F
\date 2019
