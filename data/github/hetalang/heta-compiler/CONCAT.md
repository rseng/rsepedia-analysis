# Change Log

## 0.6.10

- add `format: Dot` support
- do not support units syntax like this `{units: mL2}`
- fix bug: Simbio event without period
- Simbio storage name: nameless => nameless_ns
- use NaNMath in Julia export

## 0.6.9 - multispace

- remove support of `DSwitcher`, `CSwitcher` in `{format: SLV}`
- support of `CSwitcher`, `DSwitcher` in `{format: DBSolve}`
- support active/inactive events for `format: Matlab` 
- multispace export for `format: Matlab`
- multispace export for `format: Mrgsolve`, `format: DBSolve`, `format: SLV`
- multispace export for `format: SBML`, `format: Simbio`
- `spaceFilter` prop for `format: HetaCode`
- remove unnecessary rules from `format: Julia`

## 0.6.8

- `TimeSwitcher`, `DSwitcher`, `CSwitcher` support in mrgsolve
- support `#setScenario` and `Scenario` (no export)
- sbml export: proper sequence of listOf, remove garbage
- remove `SimSolver` export

## 0.6.7 - ready for JOSS

- error messages for unsupported switchers in SLV, DBSolve, mrgsolve
- update reserved words list with "default" (to support mrgsolve)
- extend api documentation
- export only concrete namespaces in Simbio

## 0.6.6

- remove `@SimpleTask` Class
- export to SBML L3 + timeUnits from `@TimeScale`
- replace "markdown" package by "markdown-it"
- support HetaCode export
- fix error with SLV events export
- fix log file creation
- fix error message in case of circular assignments
- support `piecewise` function in SBML, Matlab, Simbio, Julia
- support `piecewise` function in SBML module
- remove unsupported period in Simbio message

## 0.6.5

- Added `JuliaExport` format
- Add `@StopSwitcher` class as experimental
- Remove mathjs-translate deps
- Fix bug in comments inside dictionary
- ban `<-` syntax in Process expression
- remove npm dependance
- prettify internal errors
- skip units check if errors on previous stages

## 0.6.4

- `@TimeScale` component and time terms checking
- support periodic `@TimeEvent` in SBML and Simbio
- bug fix: temporally remove support of `{active: false}` events in Matlab
- bug fix: `powTransform` in SLV
- `atStart` and` {active: true}` in SimSolver
- fix renaming of function names

## 0.6.3

- check expressions for (=) sign
- bug fix: JS error in case of circular refs

## 0.6.2

- update to support SimSolver v0.3

## 0.6.1

- fix bug with units checking and dimensionless

## 0.6.0 - supports Heta standard of v0.4

- add #defineUnit instead of @UnitDef class
- checking terms for @Compartment, @Species, @Reaction
- checking units consistency for all @Record's assignments
- options.unitsCheck in declaration file
- supporting of dimensionless units in format [], 1, (1e-3)
- excel sheets numeration from 0 not from 1
- optional "id" property for #export
- advanced units checking for #export {format: Simbio, ...}
- support @_Switcher {active: false}
- fix support of logical operators: and, or, not, xor
- add `reversible` property (default `true`) to `@Process`

## 0.5.18

- Message for user to install the latest version
- Stop export if compilation errors exist
- support of ternary operator in Simbio

## 0.5.17

- replace SimSolverPlatform by Platform function in SimSolver export
- output property in @Record
- support `output` prop for SimSolver
- support `output` prop for DBSolve, SLV
- support `output` prop for Mrgsolve
- support `output` prop for Matlab
- bug fix in XLSX module

## 0.5.16

- fix bugs in sbml-module: unary minus in <plus>

## 0.5.15

- add scope prop in julia's events
- use > instead of >= in SBML and Simbio events to support run at start
- add `atStart` prop to `_Switcher`
- extend julia format
- rename #export Julia to SimSolver
- add properties to @Const: scale, lower, upper
- sbml-module: use dimensionless for simplified units
- sbml-module: support base units
- minor bug fixes

## 0.5.14

- update matlab run
- ternary operators support in Matlab
- remove ifg0 support

## 0.5.13 - SimSolver support

- add default tasks to Julia export
- notify in limitation in SBML module: StoichiometryMath, EventWithDelay, FastReaction
- fix Julia export for static and dynamic when switches compartment
- add --distDir, --metaDir, --juliaOnly to CLI options
- add "options.juliaOnly" to declaration
- week unit support in Matlab

## 0.5.12

- fix transformation in Julia: nthRoot, exponentiale, zero dynamic, ode_ priority, factorial
- add bind() method to `DSwitcher`
- export `DSwitcher` to formats: SLV, DBSolve, Matlab, Simbiology, SBML
- multispace Julia export
- use `trigger` instead of `condition` in `CSwitcher`
- export `CSwitcher` to formats: SLV, DBSolve, Matlab, Simbiology, SBML

## 0.5.11

- fix bug when `trigger` is empty
- add messages about unsupported SBML features
- add `logFormat` builder option for saving logs in JSON
- fix slv/dbsolve multiple `stop`

## 0.5.10 - skipped

## 0.5.9

- remove `repeatCount` prop from `@TimeSwitcher`
- support of multiple `@TimeSwitcher` in Matlab
- prettify code in SLV/DBSolve with `groupConstBy` prop
- fix error in DBSolve event target

## 0.5.8

- add spaceFilter for Mrgsolve, Julia
- version selection in #export {format: SBML}: support for L2V3, L2V4, L2V5
- add support @TimeSwitcher as event in SimBiology (instead of doses)
- fix bug with "not a deep clone in Expression"
- fix bug with empty period in Matlab
- include correct description of TimeSwitcher in Julia
- update structure of Julia format

## 0.5.7

- add spaceFilter for SLV, DBSolve, Simbio, Matlab, SBML
- fix UnitsExpr: empty units prop, dot in string
- output all dynamic in SLV, DBSolve

## 0.5.6 - multispace export

- default export of `units` as UnitExpr in `@UnitDef`
- multispace for JSON, YAML, XLSX: `spaceFilter` property
- catch error when XLSX file is busy

## 0.5.5 - sbml export

- fixes in unit conversion for SBML export
- support time symbol in SBML export
- remove name from SBML assignments

## 0.5.4

- add AnotherXLSX export
- fix rounding problems in units with prefixes
- fix name of container in Simbio export

## 0.5.3

- Add Export of format DBSolve
- Fix errors in unit rebase
- Add DSwitcher class
- Rename CondSwitcher to CSwitcher
- draft DSwitcher in Julia export

## 0.5.2

- rewrites Heta logs: all logs are stored in container
- fix Matlab export without event
- SBML module: support speciesType
- pretty Unit.toHTML
- faster clone() method for components
- support multiplier in Unit's hash
- minute as qsp unit

## 0.5.1 - IRT ready

- support {free: true} in #export DBSolve
- update _Component.clone() method to exclude cloning of namespace
- fixes to sbmlParse() function for clear math expressions
- "src/browser.js" is an entry point for browser apps and webpack
- method references() to _Component
- add "omit" property to #export formats: JSON, YAML, XLSX
- exclude "fs" and "path" libs from core code for browser support
- updates to nunjucks templates for easy usage in browser apps
- remove specific JS errors from Container and components
- multiple dependencies updates

## 0.5.0 - first public

- corresponds to Heta standard v0.2.4, see <https://hetalang.github.io/>
# Export formats

Following [Heta specifications](specifications/) exporting to different formats can be done by `#export` action. The following formats are implemented in Heta compiler.

- [JSON](#json)
- [YAML](#yaml)
- [SLV](#slv)
- [DBSolve](#dbsolve)
- [SBML](#sbml)
- [Simbio](#simbio)
- [Mrgsolve](#mrgsolve)
- [XLSX](#xlsx)
- [Julia](#julia)
- [Matlab](#matlab)

See also [Features support table](#features-support)

The general format for all export actions is the following:
```heta
#export {
    format: JSON, // or other supported formats, required
    filepath: path/to/output, // Relative or absolute path to generated directory or file
    spaceFilter: [nameless, another] // [] if set, namespaces out of the list will be skipped.
    ... // other options
};
```

## JSON

Export to [JSON structure](https://www.json.org/) (array) storing the content of whole platform or selected namespaces (see spaceFilter option).

### Properties

| property | type | required | default | ref | description | 
| ---------|------|----------|---------|-----|-------------|
| omit | string[] | | | | Array of properties paths to exclude from output. |
| noUnitsExpr | boolean | | false | | If `false` or not set all units will be written in format of UnitsExpr. If `true` all unit will be written in Unit array format. |

### Output files

**[filepath].json** : all content created for selected space.

**Example**

```heta
#export {
    format: JSON,
    filepath: output, // save result in file "dist/output.json"
    omit: [aux.wiki], // omit aux.wiki properties from components
    noUnitsExpr: false, // save units in format UnitsExpr
    spaceFilter: [ nameless, another ]
};
```

## YAML

Export to [YAML structure](https://yaml.org/) (array) representing the content of namespace.

### Properties

All options is the same as for [JSON format](#json).

### Output files

**[filepath].yml** : all content created for selected space.

**Example**

```heta
#export {
    format: YAML,
    filepath: output, // save result in file "dist/output.json"
    omit: [aux.wiki], // omit aux.wiki properties from components
    noUnitsExpr: false, // save units in format UnitsExpr
    spaceFilter: [ nameless, another ]
};
```

## SLV

Export to SLV format which is the model format for [DBSolveOptimum](http://insysbio.com/en/software/db-solve-optimum).

### Properties

| property | type | required | default | ref | description | 
| ---------|------|----------|---------|-----|-------------|
| eventsOff | boolean | | | | if `eventsOff = true` the switchers will not be exported to DBSolve events. |
| powTransform | "keep" / "operator" / "function" | | "keep" | | This is option describing if the transformation of x^y and pow(x, y) is required. |
| groupConstBy | string/path | | `tags[0]` | | How to group const in Initial Values of DBSolve file. Should be written in format of JSON path |

### Output files

**[filepath]/[namespace].slv** : model created based on namespace which can be opened by DBSolveOptimum.

### Known restrictions

- `Compartment` which changes in time may result in wrong ODE.
- `CSwitcher` and `DSwitcher` are not supported.
- Initialization of `Record` by expression does not work: `x1 .= k1 * A` (not supported).
- `Infinity`, `-Infinity`, `NaN` values is not supported
- boolean operators like `and`, `or`, etc. are not supported

**Example**

```heta
#export {
    format: SLV,
    filepath: model, // save results in file "dist/model.slv"
    spaceFilter: nameless, // namespace used for model generation
    eventsOff: false, // all switchers will be transformed to DBSolve events
    powTransform: keep, // use x^y and pow(x, y) without changes
    groupConstBy: "tags[1]" // use the second tag
};
```

## DBSolve

Export to DBSolve format which is the model format for [DBSolveOptimum](http://insysbio.com/en/software/db-solve-optimum).

This is the updated version of SLV export format which supports compartment volumes changed in time and initializing records by arbitrary expressions.

### Properties

| property | type | required | default | ref | description | 
| ---------|------|----------|---------|-----|-------------|
| powTransform | "keep" / "operator" / "function" | | "keep" | | This is option describing if the transformation of x^y and pow(x, y) is required. |
| groupConstBy | string/path | | `tags[0]` | | How to group const in Initial Values of DBSolve file. Should be written in format of JSON path |

### Output files

**[filepath]/[namespace].slv** : model created based on namespace which can be opened by DBSolveOptimum.

### Known restrictions

- `Infinity`, `-Infinity`, `NaN` values is not supported
- boolean operators like `and`, `or`, etc. are not supported

**Example**

```heta
#export {
    format: DBSolve,
    filepath: model, // save results in file "dist/model.slv"
    spaceFilter: nameless, // namespace used for model generation
    powTransform: keep // use x^y and pow(x, y) without changes
};
```

## SBML

Export to [SBML format](http://sbml.org/Main_Page).

### Properties

| property | type | required | default | ref | description | 
| ---------|------|----------|---------|-----|-------------|
| version | string | | L2V4 | | SBML version in format: `L2V4`. Possible values are `L2V3`, `L2V4`, `L2V5`, `L3V1`, `L3V2` |

### Output files

**[filepath]/[namespace].xml** : SBML formatted model

**Example:**

```heta
#export {
    format: SBML,
    filepath: model, // save results in file "dist/model.xml"
    spaceFilter: nameless, // namespace used for model generation
    version: L2V4 // Level 2 Version 4
};
```

## Simbio

Export to [Simbiology](https://www.mathworks.com/products/simbiology.html)/Matlab code (m files). The code can be run to create simbiology project.

### Properties

-

### Output files

**[filepath]/[namespace].m** : Code which can be run in Matlab environment to generate Simbio model.
**[filepath]/fun.m** : Auxilary mathematical functions to support Simbio code. This code should be placed in the same directory as simbio project.

**Example:**
```heta
#export {
    format: Simbio,
    filepath: model, // save results in directory "dist/model"
    spaceFilter: nameless // namespace used for model generation
};
```

## Mrgsolve

Export to [mrgsolve](http://mrgsolve.github.io/) model format (cpp file).

### Properties

- 

### Output files

**[filepath]/[namespace].cpp** : Code which can be run in mrgsolve environment.
**[filepath]/run.cpp** : Code for fast run.

### Known restrictions

- `CSwitcher` is not supported.
- `DSwitcher` is not supported.
- Initialization by MathExpr is not supported. Do not use `S1 .= x * y`.

**Example:**

```heta
#export {
    format: Mrgsolve,
    filepath: model, // save results in file "dist/model.cpp"
    spaceFilter: nameless // namespace used for model generation
};
```

## XLSX

Creation of Excel file (.xlsx) which contains components of namespace.

### Properties

| property | type | required | default | ref | description | 
| ---------|------|----------|---------|-----|-------------|
| omitRows | number | | | | If set this creates empty rows in output sheets. |
| omit | string[] | | | | Array of properties paths to exclude from output. |
| splitByClass | boolean | | | | If `true` the components will be splitted by class and saved as several sheets: one sheet per a class. |

### Output files

**[filepath].xlsx** : File which can be opened in Excel.

**Example:**

```heta
#export {
    format: XLSX,
    filepath: output, // save result in file "dist/output.xlsx"
    spaceFilter: nameless, // output all from nameless namespace
    omitRows: 5, // include 5 empty rows between header and the first line
    omit: [aux.wiki], // omit aux.wiki properties from components
    splitByClass: true // split classed to different sheets
};
```

## Julia

Creation of Julia files (.jl).

### Properties

- 

### Output files

**[filepath]/model.jl** : File storing model code.
**[filepath]/run.jl** : Code to run model.

**Example:**

```heta
#export {
    format: Julia,
    filepath: julia_code, // save result in directory "dist/julia_code"
    spaceFilter: nameless // create model based on nameless namespace
};
```

## Matlab

Creation of Matlab files (.m) which represent ODE and code to run ODE.

### Properties

- 

### Output files

**[filepath]/model.m** : File storing model code.
**[filepath]/param.m** : storing constants initialization
**[filepath]/run.m** : Code to run model.

### Known restrictions

- `CSwitcher` is not supported.

**Example:**

```heta
#export {
    format: Matlab,
    filepath: matlab_code, // save result in directory "dist/matlab_code"
    spaceFilter: nameless // create model based on nameless namespace
};
```

## Features support

*na* means "not applicable"

| | SLV | DBSolve | Julia | Mrgsolve/R | Matlab | Simbio/Matlab | SBML | JSON, YAML | XLSX |
|--|--|--|--|--|--|--|--|--|--|
|`@UnitDef` class                      |na |na |na |na |na |+ |+ |+ |+ 
|`@TimeSwitcher` class                 |+  |+  |+  |+  |+  |+ |+ |+ |+
|`@TimeSwitcher {start: 6}`                              |+ |+ |+ |+ |+ |+ |+ |+ |+
|`@TimeSwitcher {start: 0}`                              |+ |+ |+ |+ |+ |- |+ |+ |+
|`@TimeSwitcher {start: time_start}` with ref to `@Const`|+ |+ |+ |+ |+ |+ |+ |+ |+
|`@TimeSwitcher {period: 12}` infinite repeat            |+ |+ |+ |+ |+ |+ |+ |+ |+
|`@TimeSwitcher {stop: 120}` stop time for repeat        |+ |+ |+ |+ |+ |+ |+ |+ |+
|`@CSwitcher` class                                      |- |+ |+ |+ |+ |+ |+ |+ |+
|`@CSwitcher` with interpolation                         |- |- |+ |- |+ |+ |na|na|na
|`@DSwitcher` class                                      |- |+ |+ |+ |+ |+ |+ |+ |+
|`@DSwitcher` with interpolation                         |- |- |+ |- |+ |+ |na|na|na
|MathExpr: arithmetic functions                          |+ |+ |+ |+ |+ |+ |+ |+ |+
|MathExpr: boolean operators                             |- |- |+ |+ |+ |+ |+ |+ |+
|MathExpr: ternary operator                              |+ |+ |+ |- |+ |+ |+ |+ |+
|MathExpr: `piecewise` function                          |- |- |+ |- |+ |+ |+ |+ |+
|MathExpr: `e`, `pi`                                     |+ |+ |+ |+ |+ |+ |+ |+ |+
|MathExpr: `Infinity`, `NaN`                             |- |- |+ |+ |+ |+ |+ |+ |+
|Const: `Infinity`, `NaN`                                |- |- |+ |+ |+ |+ |+ |+ |+
|`@Scenario` support                                     |- |- |- |- |- |- |- |+ |+
# API references

This is a documentation for the JavaScript interface of __Heta compiler__.

For users guidance and CLI references see the main [documentation site](https://hetalang.github.io/#/heta-compiler/).

*Under developments. If you have questions contact the developers directly.*

## Main classes

- [Builder]{@link Builder}
- [Container]{@link Container}
- [Namespace]{@link Namespace}

## Modules

- [ModuleSystem]{@link ModuleSystem}
- [_Module]{@link _Module}

## Elements

- [Top]{@link Top}
- [Component]{@link Component}
- [Const]{@link Const}
- [Unit]{@link Unit}
- [Expression]{@link Expression}

## Auxiliary

- [Logger]{@link Logger}
- [Transport]{@link Transport}
- [JSONTransport]{@link JSONTransport}
- [StdoutTransport]{@link StdoutTransport}
- [StringTransport]{@link StringTransport}
[![Heta project](https://img.shields.io/badge/%CD%B1-Heta_project-blue)](https://hetalang.github.io/)
[![GitHub issues](https://img.shields.io/github/issues/hetalang/heta-compiler.svg)](https://GitHub.com/hetalang/heta-compiler/issues/)
[![Autotests](https://github.com/hetalang/heta-compiler/workflows/Autotests/badge.svg)](https://github.com/hetalang/heta-compiler/actions)
[![Coverage Status](https://coveralls.io/repos/github/hetalang/heta-compiler/badge.svg?branch=master)](https://coveralls.io/github/hetalang/heta-compiler?branch=master)
[![GitHub npm](https://img.shields.io/npm/v/heta-compiler/latest.svg)](https://www.npmjs.com/package/heta-compiler)
[![Documentation](https://img.shields.io/website?down_color=yellow&label=Documentation&up_color=green&url=https%3A%2F%2Fhetalang.github.io%2F#%2Fheta-compiler%2F)](https://hetalang.github.io/#/heta-compiler/)
[![status](https://joss.theoj.org/papers/ebff76c368d3adb720afe414ef6b29fb/status.svg)](https://joss.theoj.org/papers/ebff76c368d3adb720afe414ef6b29fb)
[![GitHub license](https://img.shields.io/github/license/hetalang/heta-compiler.svg)](https://github.com/hetalang/heta-compiler/blob/master/LICENSE)

# Heta compiler

**Heta compiler** is a software tool for the compilation of Heta-based QSP modeling platforms. Heta compiler can also be used as a JavaScript/Node package to develop modeling tools.

To read the full documentation, visit the Heta project homepage: <https://hetalang.github.io/#/heta-compiler/>.

## Table of contents

- [Introduction](#introduction)
- [How to cite](#how-to-cite)
- [Installation](#installation)
- [Supported tools](#supported-tools)
- [Usage of command line interface](#usage-of-command-line-interface)
- [Usage in NodeJS packages](#usage-in-nodejs-packages)
- [Known issues and limitations](#known-issues-and-limitations)
- [Getting help](#getting-help)
- [Contribute](#contribute)
- [License](#license)
- [Authors and history](#authors-and-history)

## Introduction

**Heta compiler** is a tool for the development of Quantitative Systems Pharmacology and Systems Biology platforms. It allows combining modules written in different formats like: [Heta language code](https://hetalang.github.io/#/specifications/), Excel sheets, [JSON](https://en.wikipedia.org/wiki/JSON)/[YAML](https://en.wikipedia.org/wiki/YAML) formatted structures, [SBML](http://sbml.org/) and transforming them into the dynamical model/models of different formats.

Quantitative Systems Pharmacology (QSP) is a discipline that uses mathematical computer models to characterize biological systems, disease processes and drug pharmacology. QSP typically deals with mechanism-based dynamical models described by ODE systems. Sometimes the modeling systems includes hundred or thousand of components and developed by a research group involving people with different expertise.

Heta compiler can be used as the framework for a QSP modeling project of any size and complexity. It can be easily integrated with existed infrastructure, workflows or used as a part of the CI/CD strategy. The pre-formulated requirements of Heta compiler are:

- storing the QSP models and data in integrated infrastructure;
- support iterative modeling platform updates (continuous development approach);
- support of models written in human-readable text and table formats;
- export models and data to different popular formats on the fly.

## How to cite

Metelkin, E., (2021). Heta compiler: a software tool for the development of large-scale QSP models and compilation into simulation formats. __Journal of Open Source Software, 6(67), 3708__, [DOI: 10.21105/joss.03708](https://doi.org/10.21105/joss.03708)

## Installation

[NodeJS](https://nodejs.org/en/) must be installed prior to Heta compiler installation. Currently **NodeJS v8/v10** are recommended.

The next steps should be taken using console (shell): **cmd**, **PowerShell**, **sh**, **bash** depending on your operating system.

1. Check Node version. It should be >= 8.0.0.
    ```bash
    node -v
    # v8.0.0 or newer
    ```

2. The latest stable version of Heta compiler can be installed from npm
    ```bash
    npm i -g heta-compiler
    ```
    **OR** The development version can be installed directly from GitHub
    ```bash
    npm i -g git+https://github.com/hetalang/heta-compiler.git
    ```

## Supported tools

>for more information see [export formats](export-formats)

Heta compiler was created to support exporting to different popular modeling formats.
One of the main development effort is to extend a list of supporting formats and allow people to have the same results working in different tools.
The current version supports the following formats:

- DBSolveOptimum .SLV files [link](http://insysbio.com/en/software/db-solve-optimum)
- SBML L2V4 [link](http://sbml.org/)
- mrgsolve .CPP files [link](https://mrgsolve.github.io/user_guide/)
- Simbiology/Matlab .M files [link](https://www.mathworks.com/products/simbiology.html)
- Matlab describing ODEs file [link](https://www.mathworks.com/help/matlab/ordinary-differential-equations.html)
- Julia format
- JSON formatted file
- YAML formatted file
- Excel sheets

## Usage of command line interface

Heta compiler comes with a built-in CLI which can be used to compile files from the command line.

>To learn more about options, see [CLI references](./cli-references)

The following is the example where we create a Heta module and compile it into SBML format. For example you want to create platform in directory "/path/to/my-platform" (target directory)

1. Create Heta file: *index.heta* in the target directory with the content:
    ```heta
    comp1 @Compartment;
    s1 @Species { compartment: comp1 };
    r1 @Reaction { actors: s1 => };

    comp1 .= 1;
    s1 .= 10;
    r1 := k1*s1*comp1;
    k1 @Const = 1e-2;

    #export {
        format: SBML,
        filepath: model
    };
    ```

2. Be sure you are in the target directory, use command `cd /path/to/my-platform` or similar if not. Compile the platform:
    ```bash
    heta build
    ```
    Heta builder takes "index.heta" file (module) as default, reads it and transforms to SBML file as declared in *index.heta*.

3. See results of compilation in directory /path/to/my-platform/**dist**.

>If you would like to load the platform form several files using `include` statement inside "index.heta", see [specifications](https://hetalang.github.io/#/specifications/include).

## Creating a Heta platform template

Platform can be structured using a prepared template of folders and pre-constructed embedded files.
Heta compiler provides the `heta init` tool for creating a such a modeling platform template.

The tool creates a draft platform including supplementary files and directories including the `platform.js` file, and files for __git__ repository.

>For more information see the [CLI references](./cli-references?id=quotheta-initquot-command) documentation.

```
heta init
$ heta init
Creating a template platform in directory: "Y:\draft"...
? Platform id (template) draft
? Platform id draft
? Platform notes (platform notes)
? Platform notes platform notes
? Platform version (v0.1.0)
? Platform version v0.1.0
? Platform license (UNLICENSED)
? Platform license UNLICENSED
? Set options (y/N)
? Set options No
? Select file types (Use arrow keys)
> heta
  heta+xlsx
  heta+xlsx extended
  xlsx
  json
  yaml
```

## Usage in NodeJS packages

Heta compiler has been written in NodeJS environment and can be used as a package for browser or server-side tools and applications.

> To learn more more, see [API docs](https://hetalang.github.io/heta-compiler/dev)
(under development).

```javascript
const { Container } = require('heta-compiler');

// platform code in Q-array format
let qArr = [
    { class: 'Compartment', id: 'comp1', assignments: {start_: '1'} },
    { class: 'Species', id: 's1', compartment: 'comp1', assignments: {start_: '10'} },
    { class: 'Reaction', id: 'r1', actors: 's1 =>', assignments: {ode_: 'k1*s1*comp1'} },
    { class: 'Const', id: 'r1', actors: 's1 =>', num: 1e-2 },
    { action: 'export', format: 'SBML', filepath: 'model' }
];

// compilation
let c = (new Container)
    .loadMany(qArr)
    .knitMany();
// get export element
let output = c.exportStorage[0]
    .make();

// print sbml code to console
console.log(output[0].content);

// check errors
console.log(c.hetaErrors());
```

## Known issues and limitations

To see a list of the supported format features, go to [features support](export-formats#features-support) table.

The tool is under active development so there are a lot of features to implement. To help us prioritize them write an [issue](https://github.com/hetalang/heta-compiler/issues).

## Getting help

 - Read Heta documentation on <https://hetalang.github.io/>
 - Use [Gitter Chatroom](https://gitter.im/hetalang/community?utm_source=readme).
 - Use [Issue Tracker](https://github.com/hetalang/heta-compiler/issues)

## Contribute

- [Source Code](https://github.com/hetalang/heta-compiler)
- [Issue Tracker](https://github.com/hetalang/heta-compiler/issues)
- See also contributing in [Heta project](https://hetalang.github.io/#/CONTRIBUTING)

## License

Licensed under the Apache License, Version 2.0. See the [LICENSE](./LICENSE) text.

## Authors and history

The original author of the project is [Evgeny Metelkin](https://github.com/metelkin). The tool was inspired by the idea that large scale dynamical systems used in QSP and SB require the specific tool which allows writing model code in unified formats and transforming them depending on one's needs: to database-like format or ODEs. Working with large models should be as easy as with the small ones.

- The initial prototype 0.1.x was developed in 2017 and named as **qs3p** (quantitative systems pharmacology programming platform). It was used in several [InSysBio LLC](http://insysbio.com) projects including [IRT](https://irt.insysbio.com/) and **Alzheimer disease consortium**.

- The next versions of **qs3p-js** used the updated format of platform components and a new approach for storing them. A set of new exporting formats was supported. The current version supports Heta code including actions, modules, namespaces. It was used as the main infrastructure for the development of the large- and middle-scale QSP platforms developed in the framework of InSysBio services.

- In 2020 the tool was renamed to **Heta compiler** and published as a Free Open Source project on [GitHub](https://GitHub.com/hetalang/heta-compiler) under Apache 2.0 license. Since then Heta compiler has been developed in the framework of [Heta project](https://hetalang.github.io/).

Copyright 2019-2021, InSysBio LLC# Migrate from v0.5 to v0.6

Heta compiler of **version 0.6** follows the [Heta standard](/specifications/) of **version 0.4**.

The new Heta specifications has several incompatibilities. To use the newer Heta compiler you should make some updates to thee platform files. To see the other updates see [change log](./CHANGELOG).

You don't need it for newly created platform.

## Update platform

1. **Check** that your platform can be build without errors in the current builder v0.5.x.

    *If you use Git you should commit the latest changes before updating formats.*

1. Install the latest version of **Heta compiler**.

    ```bash
    npm install -g heta-compiler
    ```

1. If you use declaration file **platform.json** update the property `builderVersion: ^0.6.0`.

    ```json
    {
        "builderVersion": "^0.6.0",
        "id": "my-platform",
        "notes": "platform notes",
        "version": "v0.1.0",
        ...
    }
    ```

1. If you use **qsp-units.heta** substitute it by the downloaded file [qsp-units.heta](https://raw.githubusercontent.com/hetalang/heta-compiler/master/bin/init/qsp-units.heta ':target=_blank :download')

1. Update the `@UnitDef` instances by `#defineUnit` action.

    ```heta
    // version 0.5
    Da @UnitDef { units: g/mole };  
    ```

    ```heta
    // version 0.6
    Da #defineUnit { units: g/mole };
    ```

1. Now the numeration of sheets in module of XLSX type starts from zero. Update the code if you use them.

    **include statement for xlsx type**

    ```heta
    // version 0.5
    include table.xlsx type xlsx with { sheet: 1, omitRows: 3 }
    include table.xlsx type xlsx with { sheet: 2, omitRows: 3 }
    ```

    ```heta
    // version 0.6
    include table.xlsx type xlsx with { sheet: 0, omitRows: 3 }
    include table.xlsx type xlsx with { sheet: 1, omitRows: 3 }
    ```

1. Build and check errors.
# Command Line Interface

*This page describe how to work with Heta compiler from the console (shell).*

If `heta` command is not available check your Heta compiler [installation](./README) and content of system paths (`PATH` variable in Windows).

## Table of contents

- ["heta" command](#quothetaquot-command)
- ["heta build" command](#quotheta-buildquot-command)
- ["heta init" command](#quotheta-initquot-command)
- ["heta help" command](#quotheta-helpquot-command)
- [Declaration file format](#declaration-file-format)

## "heta" command

`heta` is the prefix command for working with the tools. Writing the command alone prints the information about the tool and available options. `heta help` does the same.

```
$ heta
Usage: heta [options] [command]

Command line utilities for working with Heta compiler
  version: 0.4.31
  license: Apache-2.0

Options:
  -v, --version  output the version number
  -h, --help     output usage information

Commands:
  build [dir]    Compile Heta-based platform and create set of export files.
  init [dir]     Create template platform files in directory
  help [cmd]     display help for [cmd]
```

## "heta build" command

`heta build` runs the compilation of the platform.
It uses the main source file (index) as an initial point to compile the platform.

The default run of `heta build` (no options set, no configuration file) will do the following:

1. Looking for **index.heta** in parent working directory of shell.
2. Running parsing of index file as module of type "heta" and all files (modules) mentioned by `include` statement inside `index.heta`.
4. Creation of export files declared by `#export` actions to **dist/** directory.
5. If there are compiling errors the file **build.log** will be created in working directory.

### Running build with CLI options

CLI options allow setting specific options for build command. Use `heta build -h` to see the list of options.

At the end of the command line you can set **[dir]** path which will be used as a working directory (WD) of Heta compiler run instead of shell working directory. Absolute and relative path are possible here. If the path not set the shell WD will be used as WD of Heta.

List of `heta build` options:

| option | type | default | description |
|--|--|--|--|
| --source | \<string\> | index.heta | Path to main heta module. This allows using another name and path of index Heta module. |
| --type | \<string\> | heta | Type of source file. This option allows to select type of module which will be applied for parsing. Available values: heta/xlsx/json/yaml/sbml. |
| --debug | | | Working in debugging mode. All parsed files will be saved in JSON files in **meta** directory. |
| --units-check | | | If set all records will be checked for units consistency. |
| --skip-export | | | If set no export files will be created. |
| --julia-only | | | Run in Julia supporting mode: skip declared exports, add default export to Julia. |
| --dist-dir | \<string\> | |  Set export directory path, where to store exported files. |
| --meta-dir | \<string\> | |  Set meta directory path. |
| --log-mode | string | error | The rule in which case the log file should be created. Possible values are: never/error/always |
| -d, --declaration | string | platform | The filepath to declaration file (see below) without extension. The command will search the declaration file based on option trying a set of extensions: .json/.json5/.yml. |

#### Example 1

Let's our shell working directory is **/path/to/my-platform/**. Our Heta module is located in **src/table.xlsx** subdirectory and has a type xlsx (Excel sheet). To run compilation and save export files you should use the command.
```
heta build --source src/table.xlsx --type xlsx
```

#### Example 2

Run compilation without exporting files using "index.heta" as entry point.
```
heta build --skip-export
```

#### Example 3

Declaring working directory (WD) inside command line. 
```
heta build y:/my-platform
```

### Running build with declaration file

Declaration is a file in specific format which is located in working directory of modeling platform. As default it has a name **platform.json** and it is JSON formatted. It has two purposes:
- It annotates the developed modeling platform by some specific properties like id, notes, constibutors, repository, etc.
- To customize compiler's behavior: files location, outputs, etc. It can be used instead of CLI options.

>To use the arbitrary name of declaration file use `--declaration [filename]` option.

The declaration file can be in one of the following formats:  [JSON](https://en.wikipedia.org/wiki/JSON), [JSON5](https://json5.org/), [YAML](https://en.wikipedia.org/wiki/YAML) with the same schema.

>It is a good idea to start the model development from creation of **platform.json** file. If declaration file is set you have not to use additional options in `heta build`.

To create a draft declaration file use ["heta init" command](#"heta-init"-command). To learn all properties of declaration file see [declaration file format](#declaration-file-format).

#### Example

The following declaration file changes the default **dist** directory and displays only Heta error in console.

```json
{
  "id": "test",
  "options": {
    "distDir": "output",
    "logLevel": "error"
  }
}
```

## "heta init" command

`heta init` creates template files for QSP platform. Running the command without options will create template in current working directory. You can set another directory for creating template using optional path **[dir]** at the end of command line.

After running the command a developer should answer to series of questions about the initialized platform in prompt mode. To create the default platform use `--silent` option.

| option | type | description |
|--|--|--|
| -f, --force || This option allows rewriting the existed files and directories. |
| -s, --silent || Run initialization in silent mode with default options. |

#### Example

The current working directory is "Y:\my-folder". We are going to create a new heta project creating subdirectory "test".

```
$ heta init -f test
Creating a template platform in directory: "Y:\my-folder\test"...
? Platform id test-platform       
? Platform notes Temporal platform to test initialization
? Platform version v0.1.0
? Platform license UNLICENSED
? Set options No
? Select file types xlsx
Platform template is created.
{
  ...
}
DONE.
```

If we check the content of the created "test" directory we will see the following:

**src/index.heta** : the default heta file containing platform template. [index.heta](https://raw.githubusercontent.com/hetalang/heta-compiler/master/bin/init/index.heta ':target=_blank :download')

**src/qsp-units.heta** : pre-defined list of units. You didn't have to use it. [qsp-units.heta](https://raw.githubusercontent.com/hetalang/heta-compiler/master/bin/init/qsp-units.heta ':target=_blank :download')

**platform.json** : default declaration file.

**.gitattributes, .gitignore** : templates to help working with [Git](https://git-scm.com/)


## "heta help" command

To obtain a list of options and command description one can use command help.
It can be done by two ways:
```
heta help <command>
```
or
```
heta <command> -h
```

## Declaration file format

There are properties in declaration file which do not change compilation process. They can be used for annotation of a developed QSP platform for summarizing annotation and auxilary information.

| option | type | CLI option | default value | description |
|--|--|--|--|--|
| id | string ||| This is an unique identifier of modeling platform. Do not use spaces. *Annotation element*  |
| notes | string ||| Put a description in it. This helps people discover your package. *Annotation element.*|
| version | string ||| Version of the platform. Substantial changes to the platform should come along with changes to the version. It is recommended to follow [semver](https://semver.org/) rules. |
| keywords | string[] ||| Array of keywords for possible indexing platform. This helps people discover your package. *Annotation element.* |
| homepage | string ||| The URL to the page supporting the developed platform. *Annotation element.* |
| repository | object || {} | Specify the place where your code lives. This is helpful for people who want to contribute. |
| repository.type | string ||| Type of repository, for example: "git", "svn". *Annotation element.* |
| repository.url | string ||| The URL of source repository, for example: "https://github.com/insysbio/heta-case-mini.git". *Annotation element.* |
| license | string ||| Short udentifier under which license the platform is distributed. It is important especially for Open source platfroms. If you’re using a common license such as BSD-2-Clause or MIT, add a current [SPDX license identifier](https://spdx.org/licenses/). *Annotation element.* |
| private | boolean ||| Set true if a platform must not be shared in public repositories. *Annotation element.* |
| contributors | string[] ||| Array of authors and contributors. Please follow the format "Albert Einstein <albert.einstein@gmail.com> (https://einstein.org/cv)" *Annotation element.* |
| builderVersion | string ||| The required version of Heta compiler which should compile the code. This prevents running the old compiler for the updated QSP platforms. The string must follow the rules of [semantic versioning](https://docs.npmjs.com/about-semantic-versioning). See also [semantic versioning calculator](https://semver.npmjs.com/). |
| importModule | object | | {} | Container for the description of index module. |
| importModule.source | string | --source | index.heta | Path to index heta module. Absolute and relative filepaths are applicable. Example: "src/table.xlsx" |
| importModule.type | string | --type | heta | Type of source file. This option set type of module which will be applied for parsing. Available values: heta/xlsx/json/yaml/sbml. |
| options | object | | {} | A set of compiler options. |
| options.logMode | string | --log-mode | error | The rule in which case the log file should be created. Possible values are: never/error/always. |
| options.logPath | string | | build.log | Filepath where the log file should be created. |
| options.logLevel | string | | info | When parsing the compiler prints the messages to the shell. Here you can set a level of printing messages. Possible values: "info", "warn", "error". For example if you set "warn", only warnings and errors will be printed. |
| options.logFormat | string | | `string` | The format of saving logs to file. The default value is `string` which corresponds the format similar to console. Full list of options is : `string`, `json`.|
| options.unitsCheck | boolean | --units-check | false | If `true` all Record will be checked for units consistency. |
| options.skipExport | boolean | --skip-export | false | If `true` no export files will be created. |
| options.juliaOnly | boolean | --julia-only | false | If `true` the compilation will run Julia supporting mode. |
| options.distDir | string | --dist-dir | dist | At default all export files are created inside **dist** directory. The option can set the another target for storing outputs. |
| options.debug | boolean | --debug | false | Working in debugging mode. All parsed modules will be saved in JSON files in meta directory. |
| options.metaDir | string | --meta-dir | meta | If `options.debug` is set as `true` this option changes the target directory for meta files. |
| options.exitWithoutError | boolean | | false | If there are errors in compilation the `heta build` command return status 1 to console (which means error). If you set true this will return 0. This can be helpful for using autotesting and CI/CD automatization. |

Using neither declaration file nor CLI options is equivalent to the following declaration:
```json
{
    "builderVersion": "*",
    "options": {
        "logMode": "error",
        "logPath": "build.log",
        "logLevel": "info",
        "logFormat": "string",
        "distDir": "dist",
        "metaDir": "meta",
        "debug": false,
        "unitsCheck": false,
        "skipExport": false,
        "juliaOnly": false,
        "exitWithoutError": false
    },
    "importModule": {
        "source": "index.heta",
        "type": "heta"
    }
}
```
# TODO

## modules:

[x] xlsx
[x] json
[x] yaml
[x] sbml
[ ] markdown -> Page
[ ] csv

## exports

[x] DBSolve & SLV (DBSolve)
[x] JSON + YAML
[x] SBML L2
[x] SBML L3
[x] mrgsolve (R)
[x] simbio (Matlab)
[x] xlsx (Heta)
[x] another xlsx
[x] matlab
[x] Julia
[x] Heta-code (Heta)
[ ] csv
[ ] rxode (R)
[ ] dat (DBSolve)
[ ] ModelingToolkit (Julia)
[ ] ODEs in markdown/latex/ascii
[x] DOT language / Graphviz

## bugs

- calculate units for pow function

## features

- remove unnecessary rules
- checking units for diff eq
- check unit consistency for Species: amount/area if compartment is area
- AnyUnit for zero numbers
- checking legal functions inside Expressions and its arguments
- highlight multiline comments in Heta dictionary and array (with/without comma)
- `#move`, `#moveNS`
- atStart to exports: Matlab, DBSolve

## ideas

- check file format for modules
- syntax highlight in web
- add "ignoreCompartment" property in Species
- do not translate base units in SBML export like second => _second
- automatic creation of modifiers in SBML
- avoid insert for existed elements: get warning or #forceInsert 
- `@Dose` class to use with simbiology/mrgsolve/nonmem doses
- `heta update` => `npm i heta-compiler`
- support null for properties: highlight, parse, heta standard
- stoichiometry as `@Const` and `@Record`
- #defineFunction + function checking
- updating properties with `one::s1.assignments.start_ 5.5;`
- remove `isAmount`, `compartment` properties from `@Reaction`

## remove lodash

https://medium.com/swlh/are-we-ready-to-replace-lodash-60cd651f6c58
https://github.com/you-dont-need/You-Dont-Need-Lodash-Underscore

- merge
- defaults
- sortBy
- has
- get
- set
- flatten => flat
- omit
- pick
- times
- trim => trim
- cloneDeep !
- unique
- mapValues
- groupBy
- fromPairs
- drop
- cloneDeepWith
- intersection

### Dose class

```heta
dose1 @Dose {
  target: A,
  amount: 100,
  start: 0,
  period: 12,
  repeatCount: 4,
  rate: 0.1, // for injection
  duration: 1
};
dose2 @Dose {
  target: A,
  amount: dose_amount,
  start: start1,
  period: period1,
  repeatCount: 4,
  rate: rate1 // for injection
};
```
---
title: 'Heta compiler: a software tool for the development of large-scale QSP models and compilation into simulation formats'
tags:
  - systems biology
  - quantitative systems pharmacology
  - mathematical modeling
  - pharmacometrics
  - format conversion
  - sbml
  - heta
authors:
  - name: Evgeny Metelkin
    orcid: 0000-0001-6612-5373
    affiliation: 1
affiliations:
 - name: InSysBio LLC
   index: 1
date: 20 August 2021
bibliography: paper.bib

---

# Summary

Today mathematical modeling is becoming more and more popular in biomedicine and drug development. __Quantitative systems pharmacology__ (QSP), a relatively new research discipline, is devoted to complex models describing organisms, diseases, and drug dynamics. Designing these models presents a set of challenging methodological problems like managing a huge amount of data, dealing with large-scale models, time-consuming calculations, etc. __Heta compiler__ is a small and fast software tool written in JavaScript which manages infrastructure for QSP modeling projects. The purpose of the tool is to build and integrate QSP platform modules, to check their completeness and consistency, and then to compile everything into runnable code that will be executed in simulation software. A user can apply a command-line interface to run the model building process. Alternatively, Heta compiler can be used as a package for developing web-based applications or be integrated with simulation software.

# Statement of need

The large and still growing Systems Biology (SB) and Systems Pharmacology modeling communities utilize a variety of software tools for simulation and data analysis [@Stephanou2018; @Mentre2020; @Knight-Schrijver2016]. Usually, the modelers solve the algebraic-differential equations or perform parameters identification or sensitivity analysis. While being useful for tackling specific problems, each software tool often has no user-friendly way for routine operations like step-by-step model creation and maintenance. Furthermore, different tools have their own internal model format which cannot be reused.

This paper presents Heta compiler which provides a convenient and flexible way for the development of dynamic large-scale models based on the __Heta language__ code. The compiler translates the source modeling code into a variety of formats to be run in simulation software tools. Heta compiler also provides information on errors in a model which can be used to debug.

This tool is an effort to resolve the typical problems in a QSP project by creating a controllable working environment.
The pre-formulated requirements are:  

-	store QSP models and data in integrated infrastructure, 
-	support iterative platform updates, 
-	support of models written in human-readable formats as well as in tables, 
-	help for model code reuse and sharing,
- provide interface for storing several models in a single platform,
-	export models and data to different popular formats so it can be used out-of-the-box.

# Heta formats

`Heta compiler` has been evolving alongside the Heta language [@metelkin2019] specification. Heta is a series of human-readable and writable formats for QSP and Systems Biology projects: Heta code, table representation, JSON, and YAML notation. Heta describes dynamic models in the process-description format i.e., as interacting components that describe volumes, concentrations, amounts, rates. On the other side, it was designed to be easily transformed into ODEs or other formats.

The standardization of process-description modeling notation was also pursued in formats like SBML, CellML, Antimony. However the Heta standard can be distinguished by the specific features:  

-	Human-readable/writable code that can be used for model development or modification.
-	Easy code parsing.
-	Modularity: QSP/SB platform can be subdivided into several files and spaces for better project management.
- Multiple interchangeable representation: human-readable code, tables, JSON, YAML.
-	Reusability: modeling platforms should be easily extended for other projects.
-	Reach annotation capabilities for better model code revision.
-	Simple transformation to popular modeling formats or general-purpose ODEs.
- Support of translation from/to SBML [@Hucka2003].

## Example

Code in \autoref{fig:model-code} is an example of the Heta code describing a simple one-compartment model. The metabolic scheme of the model can be found in \autoref{fig:model-scheme}.

![Model code in Heta format: `index.heta` file. \label{fig:model-code}](model-code.png){ width=60% }

![One compartment model with two metabolites and one reaction.\label{fig:model-scheme}](model-scheme.png){ width=60% }

When the size of a model code is large it is recommended to subdivide it into modules but this model is small and can be placed into a single file, e.g `index.heta`. To build the platform with `Heta compiler` one can run the compilation with the following command in a command terminal: `heta build`.

# Features overview

`Heta compiler` includes the parser of the Heta formats and supports all features of the [Heta specifications](https://hetalang.github.io/#/specifications/) of version 0.4.1. 
It was designed to support exporting to different popular modeling formats. The current version supports the following formats: 

-	DBSolveOptimum
-	SBML of levels 2 and 3
-	mrgsolve
-	Simbiology
-	Matlab describing ODEs file
-	Julia language code
-	JSON/YAML
-	Excel sheets

`Heta compiler` can work in two modes: as a command-line tool for model development or as a library to be incorporated into third-party applications. 
The source code is written in pure JavaScript and can be run in the Node environment or in a browser.
It can be used for both: server-side and front-end applications.

To use `Heta compiler` in a modeling project a user should follow the specific formats and agreements.
Project files i.e. model code, datasets, figures, etc. should be stored in the same directory.
This directory typically includes an optional `platform.json` declaration file that stores the supplementary information of a platform as well as specific parameters for platform compilation.
The alternative way to set the options of the compiler is to use command-line arguments.
The list of them can be shown with `heta build -h` command in a shell.

![A typical workflow of `heta compiler` in a modeling project.\label{fig:workflow}](workflow.png)

# Results and discussion

`Heta compiler` can be used as the framework for a QSP modeling project of any size and complexity.  Currently, it is applied for the development and maintenance of a variety of commercial and open-source modeling projects [@faah; @covid].
`Heta compiler` has also been used for the development of web applications like the [Immune Response Template](https://irt.insysbio.com/) navigator and "PK/RO simulator" R-Shiny application [@mAb-app].

The Heta-based formats are friendly for version control systems like Git and SVN because of the modular structure and the text-based representation.
Heta compiler can easily be integrated with existing modeling infrastructure, workflows or used as a part of a CI/CD workflow.

`Heta compiler` is a part of the Heta project which is an initiative for the development of full-cycle infrastructure for modeling in pharmacology and biology: <https://hetalang.github.io>.

# References
# JOSS paper

This folder stores files for the publication in JOSS.
# title
ldalklas
kalkdlsa

ksaldklsa
l;adsl;lds
