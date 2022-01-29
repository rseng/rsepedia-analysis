# Change Log for Dash
All notable changes to `dash` will be documented in this file.
This project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- [#1876](https://github.com/plotly/dash/pull/1876) Delays finalizing `Dash.config` attributes not used in the constructor until `init_app()`.
- [#1869](https://github.com/plotly/dash/pull/1869), [#1873](https://github.com/plotly/dash/pull/1873) Upgrade Plotly.js to v2.8.3. This includes:
  - [Feature release 2.5.0](https://github.com/plotly/plotly.js/releases/tag/v2.5.0):
    - 3D traces are now compatible with `no-unsafe-eval` CSP rules.
  - [Feature release 2.6.0](https://github.com/plotly/plotly.js/releases/tag/v2.5.0):
    - Add `smith` subplots and `scattersmith` traces, for drawing Smith charts.
  - [Feature release 2.7.0](https://github.com/plotly/plotly.js/releases/tag/v2.5.0):
    - Add text data for `histogram` traces.
    - Fix an interaction between `uirevision` and `autorange` that pops up in some cases of mixed clientside / serverside figure generation.
  - [Feature release 2.8.0](https://github.com/plotly/plotly.js/releases/tag/v2.5.0):
    - Add horizontal colorbars.
    - Add text data on `heatmap` and related trace types.
    - Control legend group title fonts.
  - Patch releases [2.5.1](https://github.com/plotly/plotly.js/releases/tag/v2.5.1), [2.6.1](https://github.com/plotly/plotly.js/releases/tag/v2.6.1), [2.6.2](https://github.com/plotly/plotly.js/releases/tag/v2.6.2), [2.6.3](https://github.com/plotly/plotly.js/releases/tag/v2.6.3), [2.6.4](https://github.com/plotly/plotly.js/releases/tag/v2.6.4), [2.8.1](https://github.com/plotly/plotly.js/releases/tag/v2.8.1), [2.8.2](https://github.com/plotly/plotly.js/releases/tag/v2.8.2), and [2.8.3](https://github.com/plotly/plotly.js/releases/tag/v2.8.3) containing bugfixes.
  - This PR also upgrades various other dependencies of dash renderer and component suites.

- [#1745](https://github.com/plotly/dash/pull/1745):
    Improve our `extras_require`: there are now five options here, each with a well-defined role:
    - `dash[dev]`: for developing and building dash components.
    - `dash[testing]`: for using the `pytest` plugins in the `dash.testing` module
    - `dash[diskcache]`: required if you use `DiskcacheLongCallbackManager`
    - `dash[celery]`: required if you use `CeleryLongCallbackManager`
    - `dash[ci]`: mainly for internal use, these are additional requirements for the Dash CI tests, exposed for other component libraries to use a matching configuration.

### Added
- [#1883](https://github.com/plotly/dash/pull/1883) in DataTable added `page_current` to `persisted_props` as requested in [#1860](https://github.com/plotly/dash/issues/1860)



- [#1763](https://github.com/plotly/dash/pull/1763):
    ## Dash and Dash Renderer

    - `Input`, `State`, and `Output` now accept components instead of ID strings and Dash `callback` will auto-generate the component's ID under-the-hood if not supplied. This allows usage like:

    ```python
    my_input = dcc.Input()
    my_output = html.Div()
    app.layout = html.Div([my_input, my_output])

    @dash.callback(Output(my_output, 'children'), Input(my_input, 'value'))
    def update(value):
        return f'You have entered {value}'
    ```

    Or, if using Python >=3.8 you can use the `:=` walrus operator:
    ```python
    app.layout = html.Div([
        my_input := dcc.Input(),
        my_output := html.Div()
    ])

    @dash.callback(Output(my_output, 'children'), Input(my_input, 'value'))
    def update(value):
        return f'You have entered {value}'
    ```

  [#1894](https://github.com/plotly/dash/pull/1894) restricted this feature so auto-generated IDs are not allowed if the app uses `dash_snapshots` (a Dash Enterprise package) or if the component uses `persistence`, as this can create confusing errors. Callback definitions can still reference components in these cases, but those components must have explicit IDs.

    ## Dash Core Components

    ### Rearranged Keyword Arguments & Flexible Types
    **`Dropdown`, `RadioItem`, and `Checklist`**
    - Rearranged Keyword Arguments - `options` & `value` are now the first two keywords which means they can be supplied as positional arguments without the keyword. Supplying the keywords (`options=` and `value=`) is still supported.
    - Flexible Types - `options` can be supplied in two new forms:
      1. An array of `string|number|bool` where `label` and `value` are equal to the items in the list.
      2. A dictionary where the keys and values set as `value` and `label` respectively.

    Before:

    ```python
    dcc.Dropdown(
        options=[
            {'label': 'New York', 'value': 'New York'},
            {'label': 'Montreal', 'value': 'Montreal'},
        ],
        value='New York'
    )
    ```

    or

    ```python
    dcc.Dropdown(
        options=[
            {'label': 'New York', 'value': 'NYC'},
            {'label': 'Montreal', 'value': 'MTL'},
        ],
        value='New York'
    )
    ```

    After:

    ```python
    dcc.Dropdown(['New York', 'Montreal'], 'New York')
    ```

    Or

    ```python
    dcc.Dropdown({'NYC': 'New York', 'MTL': 'Montreal'}, 'New York')
    ```

    **`RangeSlider` & `Slider`**
    - Rearranged Keyword Arugments - `min`, `max`, and `step` are now the first three keyword arguments which means they can be supplied as positional arguments without the keyword.
    - Flexible Types
      - `step` will be calculated implicitly if not given.
      - `marks` will be auto generated if not given. It will use `min` and `max` and will respect `step` if supplied. Auto generated marks labels are SI unit formatted. Around 5 human-readable marks will be created.
      - To remove the Slider's marks, set `marks=None`.

    Before:

    ```python
    dcc.Slider(marks={1: 2, 2: 2, 3: 3})
    ```

    After:

    ```python
    dcc.Slider(min=1, max=3, step=1)
    ```

    Or equivalently:

    ```python
    dcc.Slider(1, 3, 1)
    ```

    Step can also be omitted and the `Slider` will attempt to create a nice, human readable  step with SI units and around 5 marks:

    ```python
    dcc.Slider(0, 100)
    ```

    The SI units and ranges supported in `marks` are:
    * `µ` - micro, 10⁻⁶
    * `m` - milli, 10⁻³
    * `​` (none) - 10⁰
    * `k` - kilo, 10³
    * `M` - mega, 10⁶
    * `G` - giga, 10⁹
    * `T` - tera, 10¹²
    * `P` - peta, 10¹⁵
    * `E` - exa, 10¹⁸

    _Ranges below 10µ are not supported by the Slider. This is a bug: https://github.com/plotly/dash/issues/1766_

    **`DataTable`**

    - Rearranged Keyword Arguments - `data` and `columns` the first twokeyword arguments which means they can be supplied as positional arguments without the keyword.
    - Inferred Properties - If `columns` isn't supplied then it is extracted from the the first row in `data`

    Before:

    ```python
    dash_table.DataTable(data=df.to_dict('records'), columns=[{'name': i, 'id': i} for i in df.columns])
    ```

    After:

    ```python
    dash_table.DataTable(data=df.to_dict('records'))
    ```

    ### New Component Properties

    **`Checklist` & `RadioItems`**

    - A new property `inline` appends `display: inline-block` to `labelStyle`.

    ```python
    dcc.Checklist(inline=True)
    ```

### Fixed
- [#1879](https://github.com/plotly/dash/pull/1879) Delete redundancy in pattern-matching callback implementation, specifically when `ALL` and `MATCH` wildcards are used together. This patch was submitted by an anonymous Dash Enterprise customer. Many thanks!

- [#1858](https://github.com/plotly/dash/pull/1858) Support `mini-css-extract-plugin` Webpack plugin with `@plotly/webpack-dash-dynamic-import` node package - used by components to support dash async chunks. Updated dependencies of other `@plotly` node packages.

- [#1836](https://github.com/plotly/dash/pull/1836) Fix `__all__` in dcc and table for extras: dcc download helpers and table format helpers. This also restores this functionality to the obsolete top-level packages `dash_core_components` and `dash_table`.

- [#1822](https://github.com/plotly/dash/pull/1822) Remove Radium from renderer dependencies, as part of investigating React 17 support.

- [#1779](https://github.com/plotly/dash/pull/1779):
    - Clean up our handling of serialization problems, including fixing `orjson` for Python 3.6
    - Added the ability for `dash.testing` `percy_snapshot` methods to choose widths to generate.

- [#1778](https://github.com/plotly/dash/pull/1778) DataTable: Fix React warnings stating
  that each child in a list should have a unique "key" prop

- [#1895](https://github.com/plotly/dash/pull/1895) Support debug=True if native namespace-packages are present

## [2.0.0] - 2021-08-03

## Dash and Dash Renderer

### Added
- [#1702](https://github.com/plotly/dash/pull/1702) Added a new `@app.long_callback` decorator to support callback functions that take a long time to run. See the PR and documentation for more information.
- [#1514](https://github.com/plotly/dash/pull/1514) Perform json encoding using the active plotly JSON engine.  This will default to the faster orjson encoder if the `orjson` package is installed.
- [#1736](https://github.com/plotly/dash/pull/1736) Add support for `request_refresh_jwt` hook and retry requests that used expired JWT tokens.

### Changed
- [#1679](https://github.com/plotly/dash/pull/1679) Restructure `dash`, `dash-core-components`, `dash-html-components`, and `dash-table` into a singular monorepo and move component packages into `dash`. This change makes the component modules available for import within the `dash` namespace, and simplifies the import pattern for a Dash app. From a development standpoint, all future changes to component modules will be made within the `components` directory, and relevant packages updated with the `dash-update-components` CLI command.
- [#1707](https://github.com/plotly/dash/pull/1707) Change the default value of the `compress` argument to the `dash.Dash` constructor to `False`. This change reduces CPU usage, and was made in recognition of the fact that many deployment platforms (e.g. Dash Enterprise) already apply their own compression. If deploying to an environment that does not already provide compression, the Dash 1 behavior may be restored by adding `compress=True` to the `dash.Dash` constructor.
- [#1734](https://github.com/plotly/dash/pull/1734) Added `npm run build` script to simplify build process involving `dash-renderer` and subcomponent libraries within `dash`.

### Fixed
- [#1857](https://github.com/plotly/dash/pull/1857) Fixed a regression with `dcc.Slider` and `dcc.RangeSlider` where steps were not being set to marks if None was passed as the prop argument.  Added a check to set the min and max based on the range of marks if they are not explicitly defined (for more info, see [#1843](https://github.com/plotly/dash/issues/1843) and [#1851](https://github.com/plotly/dash/issues/1843)).


## Dash Core Components
### Added

- [#1729](https://github.com/plotly/dash/pull/1729) Include F#, C#, and MATLAB in markdown code highlighting, for the upcoming .NET and MATLAB flavors of dash.

- [#1735](https://github.com/plotly/dash/pull/1735) Upgrade Plotly.js to v2.4.2. This includes:
  - [Feature release 2.3.0](https://github.com/plotly/plotly.js/releases/tag/v2.3.0):
    - More number formatting options due to `d3-format` upgrade.
    - Many new `geo` projections.
    - Improved rendering and performance of `scattergl`, `splom` and `parcoords` traces.
  - [Feature release 2.4.0](https://github.com/plotly/plotly.js/releases/tag/v2.4.0):
    - `legend.groupclick`
    - `bbox` of hover items in event data, to support custom dash-driven hover effects
  - Patch releases [2.3.1](https://github.com/plotly/plotly.js/releases/tag/v2.3.1), [2.4.1](https://github.com/plotly/plotly.js/releases/tag/v2.4.1), and [2.4.2](https://github.com/plotly/plotly.js/releases/tag/v2.4.2) containing various bug fixes.

- [#1735](https://github.com/plotly/dash/pull/1735) New `dcc.Tooltip` component. This is particularly useful for rich hover information on `dcc.Graph` charts, using the `bbox` information included in the event data in plotly.js v2.4.0

## Dash Table
### Added

- [#1729](https://github.com/plotly/dash/pull/1729) Include F#, C#, and MATLAB in markdown code highlighting, for the upcoming .NET and MATLAB flavors of dash.

## Dash HTML Components
### Removed

- [#1734](https://github.com/plotly/dash/pull/1734) Removed the following obsolete `html` elements - `<command>`, `<element>`, `<isindex>`, `<listing>`, `<multicol>`, `<nextid>`. These are obsolete and had been previously removed from the reference table.

## [1.21.0] - 2021-07-09

## Dash and Dash Renderer
### Added
- [#1675](https://github.com/plotly/dash/pull/1675) Add new `Dash` constructor argument `extra_hot_reload_paths`. This allows you to re-initialize the Python code of the app when non-Python files change, if you know that these files impact the app.

### Changed
- [#1675](https://github.com/plotly/dash/pull/1675) Remove the constraint that `requests_pathname_prefix` ends with `routes_pathname_prefix`. When you are serving your app behind a reverse proxy that rewrites URLs that constraint needs to be violated.
- [#1611](https://github.com/plotly/dash/pull/1611) and [#1685](https://github.com/plotly/dash/pull/1685) Package dash-renderer artifacts and dependencies with Dash, and source renderer resources from within Dash.
- [#1567](https://github.com/plotly/dash/pull/1567) Julia component generator puts components into `src/jl` - fixes an issue on case-insensitive filesystems when the component name and module name match (modulo case) and no prefix is used. Also reduces JS/Julia clutter in the overloaded `src` directory.

### Fixed
- [#1664](https://github.com/plotly/dash/pull/1664) Fix [#1649](https://github.com/plotly/dash/issues/1649), makes the devtools readable with a dark theme.
- [#1640](https://github.com/plotly/dash/pull/1640) Fix [#1475](https://github.com/plotly/dash/issues/1475), missing `timing_information` after certain modifications to Flask behavior

## Dash Core Components
### Fixed

- [#963](https://github.com/plotly/dash-core-components/pull/963) Fixes [#885](https://github.com/plotly/dash-core-components/issues/885)

  This applies the fix from [#878](https://github.com/plotly/dash-core-components/pull/878) to the RangeSlider.
  It not only fixes the bug where the tooltips were visible when slider was not, but it also reduces the lag in the
  tooltip when the slider handles are moved.

### Updated
- [#939](https://github.com/plotly/dash-core-components/pull/939) Upgrade Plotly.js to v2.2.1. Note that this is a major version upgrade to Plotly.js, however we are not treating this as a breaking change for DCC as the majority of breaking changes in Plotly.js do not affect the Dash API. The one exception is that several trace types that have long been deprecated are removed entirely.
  - [Major release 2.0.0](https://github.com/plotly/plotly.js/releases/tag/v2.0.0):
    - Stop exporting d3 as `Plotly.d3`, and remove many other deep pieces of the public API. This does not affect the `dcc.Graph` component, but if you make use of `Plotly` from the global scope in some other way you may be affected.
    - Drop the deprecated trace types `contourgl` and `area`, as well as legacy pre-`scatterpolar` polar attributes `bar.r`, `bar.t`, `scatter.r`, `scatter.t`, `layout.radialaxis`, `layout.angularaxis`. Use `scatterpolar`, `barpolar`, and `polar` subplots instead.
    - `heatmapgl` and `pointcloud` trace types, and the `transform` attribute are deprecated, and will be removed in a future release.
    - Increase CSP safety by removing function constructors. 3D plots still use function constructors, but if you place one of the non-3D bundles (including the new `strict` bundle) in your `assets` folder you will have no function constructors.
    - Remove "Aa" text in legends.
    - Default `hovermode` to "closest".
    - Default `textposition` to "auto" in `bar` traces. If you previously used the `bar.text` attribute for hover only, you will need to explicitly set `textposition="none"`.
    - Add `bar.marker.pattern`, `image.zsmooth`, and various other features and bugfixes.
  - [Feature release 2.1.0](https://github.com/plotly/plotly.js/releases/tag/v2.1.0):
    - New `icicle` trace type.
    - New `legendrank` trace attribute.
    - Several other additions and bug fixes.
  - [Feature release 2.2.0](https://github.com/plotly/plotly.js/releases/tag/v2.2.0):
    - Legend group titles
    - Half-year directive (`%h`) for date formatting
    - Several other bug fixes and performance improvements
  - [Patch release 2.2.1](https://github.com/plotly/plotly.js/releases/tag/v2.2.1) containing a security fix.

### Added
- [#932](https://github.com/plotly/dash-core-components/pull/932) Adds a new copy to clipboard component.
- [#948](https://github.com/plotly/dash-core-components/pull/948)] Adds `disabled_days` prop to `DatePickerRange` and `DatePickerSingle` components. With this prop you can specify days that should be made unselectable in the date picker, in addition to those that fall outside of the range specified by `min_date_allowed` and `max_date_allowed`.

### Changed
- [#972](https://github.com/plotly/dash-core-components/pull/972) Updated R package vignettes and `dash-info.yaml` to regenerate examples without attaching now-deprecated core component packages (`dashHtmlComponents`, `dashCoreComponents`, or `dashTable`).

## Dash HTML Components
### Changed
- [#194](https://github.com/plotly/dash-html-components/pull/194) Updated dependencies and build process
- [#190](https://github.com/plotly/dash-core-components/pull/190) Updated R package vignettes and `dash-info.yaml` to regenerate examples without attaching now-deprecated core component packages (`dashHtmlComponents`, `dashCoreComponents`, or `dashTable`).

## Dash Table
### Fixed
- [#907](https://github.com/plotly/dash-table/pull/907)
  - Fix a bug where pagination did not work or was not visible. [#834](https://github.com/plotly/dash-table/issues/834)
  - Fix a bug where if you are on a page that no longer exists after the data is updated, no data is displayed. [#892](https://github.com/plotly/dash-table/issues/892)


### Added
- [#916](https://github.com/plotly/dash-table/pull/916)
  - Added `html` option to `markdown_options` prop. This enables the use of html tags in markdown text.

- [#545](https://github.com/plotly/dash-table/issues/545)
    - Case insensitive filtering
    - New props: `filter_options` - to control case of all filters, `columns.filter_options` - to control filter case for each column
    - New operators: `i=`, `ieq`, `i>=`, `ige`, `i>`, `igt`, `i<=`, `ile`, `i<`, `ilt`, `i!=`, `ine`, `icontains` - for case-insensitive filtering, `s=`, `seq`, `s>=`, `sge`, `s>`, `sgt`, `s<=`, `sle`, `s<`, `slt`, `s!=`, `sne`, `scontains` - to force case-sensitive filtering on case-insensitive columns

### Changed
- [#918](https://github.com/plotly/dash-core-components/pull/918) Updated all dependencies. In particular the `highlight.js` upgrade changes code highlighting in markdown: we have long used their "github" style, this has been updated to more closely match current github styles.
- [#901](https://github.com/plotly/dash-core-components/pull/901) Updated R package `dash-info.yaml` to regenerate example without attaching now-deprecated core component packages (`dashHtmlComponents`, `dashCoreComponents`, or `dashTable`).


## [1.20.0] - 2021-04-08

## Dash and Dash Renderer
### Changed
- [#1531](https://github.com/plotly/dash/pull/1531) Update the format of the docstrings to make them easier to read in the reference pages of Dash Docs and in the console. This also addresses [#1205](https://github.com/plotly/dash/issues/1205)
- [#1553](https://github.com/plotly/dash/pull/1553) Increase the z-index of the Dash error menu from 1001 to 1100 in order to make sure it appears above Bootstrap components.

### Fixed
- [#1546](https://github.com/plotly/dash/pull/1546) Validate callback request `outputs` vs `output` to avoid a perceived security issue.

## Dash Core Components
### Added
- [#863](https://github.com/plotly/dash-core-components/pull/863) Adds a new `Download` component. Along with this several utility functions are added to help construct the appropriate data format:
  - `dcc.send_file` - send a file from disk
  - `dcc.send_data_frame` - send a `DataFrame`, using one of its writer methods
  - `dcc.send_bytes` - send a bytestring or the result of a bytestring writer
  - `dcc.send_string` - send a string or the result of a string writer

### Changed
- [#923](https://github.com/plotly/dash-core-components/pull/923)
  Set `autoComplete` to off in `dcc.Dropdown`. This fixes [#808](https://github.com/plotly/dash-core-components/issues/808)

### Fixed
- [#930](https://github.com/plotly/dash-core-components/pull/930) Fixed a bug [#867](https://github.com/plotly/dash-core-components/issues/867) with `DatePickerRange` that would sometimes shift the allowed dates by one day.
- [#934](https://github.com/plotly/dash-core-components/pull/934) Fixed a bug in `EnhancedTab` component that ignored `disabled_className` property

## Dash HTML Components
### Fixed
- [#179](https://github.com/plotly/dash-html-components/pull/179) - Fixes [#77](https://github.com/plotly/dash-html-components/issues/77) Added `allow` and `referrerPolicy` properties to `html.Iframe`

- [#178](https://github.com/plotly/dash-html-components/pull/178) - Fix [#161](https://github.com/plotly/dash-html-components/issues/161) <object> `data` property, and fix [#129](https://github.com/plotly/dash-html-components/issues/129) obsolete, deprecated, and discouraged elements. No elements were removed, but comments were added to the documentation about these elements detailing their limitations.

## Dash Table
### Changed
- [#862](https://github.com/plotly/dash-table/pull/862) - update docstrings per https://github.com/plotly/dash/issues/1205
- [#878](https://github.com/plotly/dash-table/pull/878) - update build process to use Webpack 5 and other latest dependencies

## [1.19.0] - 2021-01-19

## Dash and Dash Renderer
### Added
- [#1508](https://github.com/plotly/dash/pull/1508) Fix [#1403](https://github.com/plotly/dash/issues/1403): Adds an x button
to close the error messages box.
- [#1525](https://github.com/plotly/dash/pull/1525) Adds support for callbacks which have overlapping inputs and outputs. Combined with `dash.callback_context` this addresses many use cases which require circular callbacks.

### Changed
- [#1503](https://github.com/plotly/dash/pull/1506) Fix [#1466](https://github.com/plotly/dash/issues/1466): loosen `dash[testing]` requirements for easier integration in external projects. This PR also bumps many `dash[dev]` requirements.

### Fixed
- [#1530](https://github.com/plotly/dash/pull/1530) Dedent error messages more carefully.
- [#1527](https://github.com/plotly/dash/issues/1527) 🐛 `get_asset_url` now pulls from an external source if `assets_external_path` is set.
  - updated `_add_assets_resource` to build asset urls the same way as `get_asset_url`.
  - updated doc string for `assets_external_path` Dash argument to be more clear that it will always be joined with the `assets_url_path` argument when determining the url to an external asset.
- [#1493](https://github.com/plotly/dash/pull/1493) Fix [#1143](https://github.com/plotly/dash/issues/1143), a bug where having a file with one of several common names (test.py, code.py, org.py, etc) that imports a dash component package would make `import dash` fail with a cryptic error message asking whether you have a file named "dash.py"

## Dash Core Components
### Fixed
- [#905](https://github.com/plotly/dash-core-components/pull/905) Make sure the `figure` prop of `dcc.Graph` receives updates from user interactions in the graph, by using the same `layout` object as provided in the prop rather than cloning it. Fixes [#879](https://github.com/plotly/dash-core-components/issues/879).
- [#903](https://github.com/plotly/dash-core-components/pull/903) Part of fixing dash import bug https://github.com/plotly/dash/issues/1143

### Updated
- [#911](https://github.com/plotly/dash-core-components/pull/911), [#906](https://github.com/plotly/dash-core-components/pull/906)
  - Upgraded Plotly.js to [1.58.4](https://github.com/plotly/plotly.js/releases/tag/v1.58.4)
    - Patch Release [1.58.4](https://github.com/plotly/plotly.js/releases/tag/v1.58.4)
    - Patch Release [1.58.3](https://github.com/plotly/plotly.js/releases/tag/v1.58.3)

### Added
- [#888](https://github.com/plotly/dash-core-components/pull/888) Adds a `drag_value` prop to `dcc.Slider`to be able to fire callbacks from dragging and releasing the slider.

## Dash HTML Components
### Fixed
- [#169](https://github.com/plotly/dash-html-components/pull/169) - part of fixing dash import bug https://github.com/plotly/dash/issues/1143

## Dash Table
### Fixed
- [#854](https://github.com/plotly/dash-table/pull/854) - part of fixing dash import bug https://github.com/plotly/dash/issues/1143

## [1.18.1] - 2020-12-09

## [1.18.0] - 2020-12-07

## [1.17.0] - 2020-10-29
### Changed
- [#1442](https://github.com/plotly/dash/pull/1442) Update from React 16.13.0 to 16.14.0
### Fixed
- [#1434](https://github.com/plotly/dash/pull/1434) Fix [#1432](https://github.com/plotly/dash/issues/1432) for Julia to import non-core component packages without possible errors.

### Changed
- [#1448](https://github.com/plotly/dash/pull/1448) Provide a hint in the callback error when the user forgot to make `app.callback(...)` a decorator.

## [1.16.3] - 2020-10-07
### Fixed
- [#1426](https://github.com/plotly/dash/pull/1426) Fix a regression caused by `flask-compress==1.6.0` causing performance degradation on server requests

## [1.16.2] - 2020-09-25
### Fixed
- [#1415](https://github.com/plotly/dash/pull/1415) Fix a regression with some layouts callbacks involving dcc.Tabs, not yet loaded dash_table.DataTable and dcc.Graph to not be called
- [#1416](https://github.com/plotly/dash/pull/1416) Make callback graph more robust for complex apps and some specific props (`width` in particular) that previously caused errors.

## [1.16.1] - 2020-09-16
### Changed
- [#1376](https://github.com/plotly/dash/pull/1376) Extends the `getTransform` logic in the renderer to handle `persistenceTransforms` for both nested and non-nested persisted props. This was used to to fix [dcc#700](https://github.com/plotly/dash-core-components/issues/700) in conjunction with [dcc#854](https://github.com/plotly/dash-core-components/pull/854) by using persistenceTransforms to strip the time part of the datetime so that datepickers can persist when defined in callbacks.

### Fixed
- [#1408](https://github.com/plotly/dash/pull/1408) Fixes a bug where the callback graph layout would reset whenever a callback fired, losing user-initiated layout changes ([#1402](https://github.com/plotly/dash/issues/1402)) or creating a new force layout ([#1401](https://github.com/plotly/dash/issues/1401))

## [1.16.0] - 2020-09-03
### Added
- [#1371](https://github.com/plotly/dash/pull/1371) You can now get [CSP `script-src` hashes](https://developer.mozilla.org/en-US/docs/Web/HTTP/Headers/Content-Security-Policy/script-src) of all added inline scripts by calling `app.csp_hashes()` (both Dash internal inline scripts, and those added with `app.clientside_callback`) .

### Changed
- [#1385](https://github.com/plotly/dash/pull/1385) Closes [#1350](https://github.com/plotly/dash/issues/1350) and fixes a previously undefined callback behavior when multiple elements are stacked on top of one another and their `n_clicks` props are used as inputs of the same callback. The callback will now trigger once with all the triggered `n_clicks` props changes.
- [#1179](https://github.com/plotly/dash/pull/1179) New and improved callback graph in the debug menu. Now based on Cytoscape for much more interactivity, plus callback profiling including number of calls, fine-grained time information, bytes sent and received, and more. You can even add custom timing information on the server with `callback_context.record_timing(name, seconds)`

### Fixed
- [#1384](https://github.com/plotly/dash/pull/1384) Fixed a bug introduced by [#1180](https://github.com/plotly/dash/pull/1180) breaking use of `prevent_initial_call` as a positional arg in callback definitions

## [1.15.0] - 2020-08-25
### Added
- [#1355](https://github.com/plotly/dash/pull/1355) Removed redundant log message and consolidated logger initialization. You can now control the log level - for example suppress informational messages from Dash with `app.logger.setLevel(logging.WARNING)`.
- [#1253](https://github.com/plotly/dash/pull/1253), [#1377](https://github.com/plotly/dash/pull/1377) Added experimental `--jl-prefix` option to `dash-generate-components`, optionally generates Julia version of components and corresponding Julia package

### Changed
- [#1180](https://github.com/plotly/dash/pull/1180) and [#1375](https://github.com/plotly/dash/pull/1375) `Input`, `Output`, and `State` in callback definitions don't need to be in lists. You still need to provide `Output` items first, then `Input` items, then `State`, and the list form is still supported. In particular, if you want to return a single output item wrapped in a length-1 list, you should still wrap the `Output` in a list. This can be useful for procedurally-generated callbacks.
- [#1368](https://github.com/plotly/dash/pull/1368) Updated pytest to v6.0.1. To avoid deprecation warnings, this also updated pytest-sugar to 0.9.4 and pytest-mock to 3.2.0. The pytest-mock update only effects python >= 3.0. Pytest-mock remains pinned at 2.0.0 for python == 2.7.

## [1.14.0] - 2020-07-27
### Added
- [#1343](https://github.com/plotly/dash/pull/1343) Add `title` parameter to set the
document title. This is the recommended alternative to setting app.title or overriding
the index HTML.
- [#1315](https://github.com/plotly/dash/pull/1315) Add `update_title` parameter to set or disable the "Updating...." document title during updates. Closes [#856](https://github.com/plotly/dash/issues/856) and [#732](https://github.com/plotly/dash/issues/732)

## [1.13.4] - 2020-06-25
### Fixed
- [#1310](https://github.com/plotly/dash/pull/1310) Fix a regression since 1.13.0 preventing more than one loading state from being shown at a time.

## [1.13.3] - 2020-06-19

## [1.13.2] - 2020-06-18
### Fixed
- [#1305](https://github.com/plotly/dash/issues/1305)
    - Fix regression that causes crash when `FLASK_ENV` is modified during app execution
    - Fix regression that caused tests using `_wait_for_callbacks` to fail

## [1.13.1] - 2020-06-17

## [1.13.0] - 2020-06-17
### Added
- [#1289](https://github.com/plotly/dash/pull/1289) Supports `DASH_PROXY` env var to tell `app.run_server` to report the correct URL to view your app, when it's being proxied. Throws an error if the proxy is incompatible with the host and port you've given the server.
- [#1240](https://github.com/plotly/dash/pull/1240) Adds `callback_context` to clientside callbacks (e.g. `dash_clientside.callback_context.triggered`). Supports `triggered`, `inputs`, `inputs_list`, `states`, and `states_list`, all of which closely resemble their serverside cousins.

### Changed
- [#1237](https://github.com/plotly/dash/pull/1237) Closes [#920](https://github.com/plotly/dash/issues/920): Converts hot reload fetch failures into a server status indicator showing whether the latest fetch succeeded or failed. Callback fetch failures still appear as errors but have a clearer message.
- [#1254](https://github.com/plotly/dash/pull/1254) Modifies the callback chain implementation and improves performance for apps with a lot of components

### Fixed
- [#1255](https://github.com/plotly/dash/pull/1255) Hard hot reload targets only the current window, not the top - so if your app is in an iframe you will only reload the app
- [#1249](https://github.com/plotly/dash/pull/1249) Fixes [#919](https://github.com/plotly/dash/issues/919) so `dash.testing` is compatible with more `pytest` plugins, particularly `pytest-flake8` and `pytest-black`.
- [#1248](https://github.com/plotly/dash/pull/1248) Fixes [#1245](https://github.com/plotly/dash/issues/1245), so you can use prop persistence with components that have dict IDs, ie for pattern-matching callbacks.
- [#1185](https://github.com/plotly/dash/pull/1185) Sort asset directories, same as we sort files inside those directories. This way if you need your assets loaded in a certain order, you can add prefixes to subdirectory names and enforce that order.
- [#1288](https://github.com/plotly/dash/pull/1288) Closes [#1285](https://github.com/plotly/dash/issues/1285): Debug=True should work in the __main__ module.

## [1.12.0] - 2020-05-05
### Added
- [#1228](https://github.com/plotly/dash/pull/1228) Adds control over firing callbacks on page (or layout chunk) load. Individual callbacks can have their initial calls disabled in their definition `@app.callback(..., prevent_initial_call=True)` and similar for `app.clientside_callback`. The app-wide default can also be changed with `app=Dash(prevent_initial_callbacks=True)`, then individual callbacks may disable this behavior.
- [#1201](https://github.com/plotly/dash/pull/1201) New attribute `app.validation_layout` allows you to create a multi-page app without `suppress_callback_exceptions=True` or layout function tricks. Set this to a component layout containing the superset of all IDs on all pages in your app.
- [#1078](https://github.com/plotly/dash/pull/1078) Permit usage of arbitrary file extensions for assets within component libraries

### Fixed
- [#1224](https://github.com/plotly/dash/pull/1224) Fixes [#1223](https://github.com/plotly/dash/issues/1223), a very specific situation in which initial callbacks will not fire.
- [#1220](https://github.com/plotly/dash/pull/1220) Fixes [#1216](https://github.com/plotly/dash/issues/1216), a set of related issues about pattern-matching callbacks with `ALL` wildcards in their `Output` which would fail if no components matched the pattern.
- [#1212](https://github.com/plotly/dash/pull/1212) Fixes [#1200](https://github.com/plotly/dash/issues/1200) - prior to Dash 1.11, if none of the inputs to a callback were on the page, it was not an error. This was, and is now again, treated as though the callback raised PreventUpdate. The one exception to this is with pattern-matching callbacks, when every Input uses a multi-value wildcard (ALL or ALLSMALLER), and every Output is on the page. In that case the callback fires as usual.
- [#1201](https://github.com/plotly/dash/pull/1201) Fixes [#1193](https://github.com/plotly/dash/issues/1193) - prior to Dash 1.11, you could use `flask.has_request_context() == False` inside an `app.layout` function to provide a special layout containing all IDs for validation purposes in a multi-page app. Dash 1.11 broke this when we moved most of this validation into the renderer. This change makes it work again.

## [1.11.0] - 2020-04-10
### Added
- [#1103](https://github.com/plotly/dash/pull/1103) Pattern-matching IDs and callbacks. Component IDs can be dictionaries, and callbacks can reference patterns of components, using three different wildcards: `ALL`, `MATCH`, and `ALLSMALLER`, available from `dash.dependencies`. This lets you create components on demand, and have callbacks respond to any and all of them. To help with this, `dash.callback_context` gets three new entries: `outputs_list`, `inputs_list`, and `states_list`, which contain all the ids, properties, and except for the outputs, the property values from all matched components.
- [#1103](https://github.com/plotly/dash/pull/1103) `dash.testing` option `--pause`: after opening the dash app in a test, will invoke `pdb` for live debugging of both Javascript and Python. Use with a single test case like `pytest -k cbwc001 --pause`.

### Changed
- [#1103](https://github.com/plotly/dash/pull/1103) Multiple changes to the callback pipeline:
  - `dash.callback_context.triggered` now does NOT reflect any initial values, and DOES reflect EVERY value which has been changed either by activity in the app or as a result of a previous callback. That means that the initial call of a callback with no prerequisite callbacks will list nothing as triggering. For backward compatibility, we continue to provide a length-1 list for `triggered`, but its `id` and `property` are blank strings, and `bool(triggered)` is `False`.
  - A user interaction which returns the same property value as was previously present will not trigger the component to re-render, nor trigger callbacks using that property as an input.
  - Callback validation is now mostly done in the browser, rather than in Python. A few things - mostly type validation, like ensuring IDs are strings or dicts and properties are strings - are still done in Python, but most others, like ensuring outputs are unique, inputs and outputs don't overlap, and (if desired) that IDs are present in the layout, are done in the browser. This means you can define callbacks BEFORE the layout and still validate IDs to the layout; and while developing an app, most errors in callback definitions will not halt the app.

### Fixed
- [#1103](https://github.com/plotly/dash/pull/1103) Fixed multiple bugs with chained callbacks either not triggering, inconsistently triggering, or triggering multiple times. This includes: [#635](https://github.com/plotly/dash/issues/635), [#832](https://github.com/plotly/dash/issues/832), [#1053](https://github.com/plotly/dash/issues/1053), [#1071](https://github.com/plotly/dash/issues/1071), and [#1084](https://github.com/plotly/dash/issues/1084). Also fixed [#1105](https://github.com/plotly/dash/issues/1105): async components that aren't rendered by the page (for example in a background Tab) would block the app from executing callbacks.

## [1.10.0] - 2020-04-01
### Added
- [#1134](https://github.com/plotly/dash/pull/1134) Allow `dash.run_server()` host and port parameters to be set with environment variables HOST & PORT, respectively

### Changed
- [#1145](https://github.com/plotly/dash/pull/1145) Update from React 16.8.6 to 16.13.0

### Fixed
- [#1142](https://github.com/plotly/dash/pull/1142) [Persistence](https://dash.plot.ly/persistence): Also persist 0, empty string etc

## [1.9.1] - 2020-02-27
### Added
- [#1133](https://github.com/plotly/dash/pull/1133) Allow the `compress` config variable to be set with an environment variable with DASH_COMPRESS=FALSE

## [1.9.0] - 2020-02-04
### Fixed
- [#1080](https://github.com/plotly/dash/pull/1080) Handle case where dash fails to load when used inside an iframe with a sandbox attribute that only has allow-scripts

## [1.8.0] - 2020-01-14
### Added
- [#1073](https://github.com/plotly/dash/pull/1073) Two new functions to simplify usage handling URLs and pathnames: `app.get_relative_path` & `app.trim_relative_path`.
These functions are particularly useful for apps deployed on Dash Enterprise where the apps served under a URL prefix (the app name) which is unlike apps served on localhost:8050.
    - `app.get_relative_path` returns a path with the config setting `requests_pathname_prefix` prefixed. Use `app.get_relative_path` anywhere you would provide a relative pathname, like `dcc.Link(href=app.relative_path('/page-2'))` or even as an alternative to `app.get_asset_url` with e.g. `html.Img(src=app.get_relative_path('/assets/logo.png'))`.
    - `app.trim_relative_path` a path with `requests_pathname_prefix` and leading & trailing
    slashes stripped from it. Use this function in callbacks that deal with `dcc.Location` `pathname`
    routing.
    Example usage:
    ```python
    app.layout = html.Div([
        dcc.Location(id='url'),
        html.Div(id='content')
    ])
    @app.callback(Output('content', 'children'), [Input('url', 'pathname')])
    def display_content(path):
        page_name = app.strip_relative_path(path)
        if not page_name:  # None or ''
            return html.Div([
                html.Img(src=app.get_relative_path('/assets/logo.png')),
                dcc.Link(href=app.get_relative_path('/page-1')),
                dcc.Link(href=app.get_relative_path('/page-2')),
            ])
        elif page_name == 'page-1':
            return chapters.page_1
        if page_name == "page-2":
            return chapters.page_2
    ```

### Changed
- [#1035](https://github.com/plotly/dash/pull/1035) Simplify our build process.
- [#1074](https://github.com/plotly/dash/pull/1074) Error messages when providing an incorrect property to a component have been improved: they now specify the component type, library, version, and ID (if available).

### Fixed
- [#1037](https://github.com/plotly/dash/pull/1037) Fix no_update test to allow copies, such as those stored and retrieved from a cache.

## [1.7.0] - 2019-11-27
### Added
- [#967](https://github.com/plotly/dash/pull/967) Add support for defining
clientside JavaScript callbacks via inline strings.
- [#1020](https://github.com/plotly/dash/pull/1020) Allow `visit_and_snapshot` API in `dash.testing.browser` to stay on the page so you can run other checks.

### Changed
- [#1026](https://github.com/plotly/dash/pull/1026) Better error message when you forget to wrap multiple `children` in an array, and they get passed to other props.

### Fixed
- [#1018](https://github.com/plotly/dash/pull/1006) Fix the `dash.testing` **stop** API with process application runner in Python2. Use `kill()` instead of `communicate()` to avoid hanging.
- [#1027](https://github.com/plotly/dash/pull/1027) Fix bug with renderer callback lock never resolving with non-rendered async component using the asyncDecorator

## [1.6.1] - 2019-11-14
### Fixed
- [#1006](https://github.com/plotly/dash/pull/1006) Fix IE11 / ES5 compatibility and validation issues
- [#1006](https://github.com/plotly/dash/pull/1006) Fix bug with renderer wrapper component TreeContainer to prevent useless re-renders
- [#1001](https://github.com/plotly/dash/pull/1001)
  - Fix and improve the `clear_input()` API in `dash.testing`, so it's more robust handling react `input`.
  - make the `percy_snapshot()` API more robust, and the timeout of `wait_for_callbacks` (if set to True) will not fail the snapshot execution, but logged as potential error.

## [1.6.0] - 2019-11-04
### Fixed
- [#999](https://github.com/plotly/dash/pull/999) Fix fingerprint for component suites with `metadata` in version.
- [#983](https://github.com/plotly/dash/pull/983) Fix the assets loading issues when dashR application runner is handling with an app defined by string chunk.

## [1.5.1] - 2019-10-29
### Fixed
- [#987](https://github.com/plotly/dash/pull/987) Fix cache string handling for component suites with nested folders in their packages.
- [#986](https://github.com/plotly/dash/pull/986) Fix a bug with evaluation of `_force_eager_loading` when application is loaded with gunicorn

## [1.5.0] - 2019-10-29
### Added
- [#964](https://github.com/plotly/dash/pull/964) Adds support for preventing updates in clientside functions.
  - Reject all updates with `throw window.dash_clientside.PreventUpdate;`
  - Reject a single output by returning `window.dash_clientside.no_update`
- [#899](https://github.com/plotly/dash/pull/899) Add support for async dependencies and components
- [#973](https://github.com/plotly/dash/pull/973) Adds support for resource caching and adds a fallback caching mechanism through etag

### Fixed
- [#974](https://github.com/plotly/dash/pull/974) Fix and improve a percy snapshot behavior issue we found in dash-docs testing. It adds a flag `wait_for_callbacks` to ensure that, in the context of a dash app testing, the percy snapshot action will happen only after all callbacks get fired.

## [1.4.1] - 2019-10-17
### Fixed
- [#969](https://github.com/plotly/dash/pull/969) Fix warnings emitted by react devtools coming from our own devtools components.

## [1.4.0] - 2019-10-08
### Added
- [#948](https://github.com/plotly/dash/pull/948) Support setting working directory for R apps run using the `dashr` fixture, primarily useful for tests with assets. `dashr.start_server` supports a `cwd` argument to set an explicit working directory, and has smarter defaults when it's omitted: if `app` is a path to an R script, uses the directory of that path; if `app` is a string, uses the directory the test file itself is in.
- [#944](https://github.com/plotly/dash/pull/944)
  - Relevant `dash.testing` methods can now be called with either an element or a CSS selector: `select_dcc_dropdown`, `multiple_click`, `clear_input`, `zoom_in_graph_by_ratio`, `click_at_coord_fractions`.
  - Three new `dash.testing` methods: `clear_local_storage`, `clear_session_storage`, and `clear_storage` (to clear both together)
- [#937](https://github.com/plotly/dash/pull/937) `dash.testing` adds two APIs `zoom_in_graph_by_ratio` and `click_at_coord_fractions` about advanced interactions using mouse `ActionChain`
- [#938](https://github.com/plotly/dash/issues/938) Add debugging traces to dash backend about serving component suites, to verify the installed packages whenever in doubt.

### Fixed
- [#944](https://github.com/plotly/dash/pull/944) Fix a bug with persistence being toggled on/off on an existing component.

## [1.3.1] - 2019-09-19
### Changed
- Bump dash-core-components version from 1.2.0 to [1.2.1](https://github.com/plotly/dash-core-components/blob/master/CHANGELOG.md#120---2019-09-19)

## [1.3.0] - 2019-09-17
### Added
- [#923](https://github.com/plotly/dash/pull/923) Add one configuration `--percy-assets` in `pytest` to specify extra application assets path if needed.

- [#918](https://github.com/plotly/dash/pull/918) Add `wait_for_element_by_id` and `visit_and_snapshot` APIs in browser, add `raw_command` option (with higher priority than the default waitress one) and optional `start_timeout` argument to handle large applications within the process runner.

- [#903](https://github.com/plotly/dash/pull/903) Persistence: enable props edited by the user to persist across recreating the component or reloading the page. Components need to define three new props: `persistence`, `persisted_props`, and `persistence_type` as described in the lead comment of `src/persistence.js`. App developers then enable this behavior by, in the simplest case, setting `persistence: true` on the component. First use case is table, see [dash-table#566](https://github.com/plotly/dash-table/pull/566)

### Changed
- Bump dash-table version from 4.2.0 to [4.3.0](https://github.com/plotly/dash-table/blob/master/CHANGELOG.md#430---2019-09-17)
- Bump dash-core-components version from 1.1.2 to [1.2.0](https://github.com/plotly/dash-core-components/blob/master/CHANGELOG.md#120---2019-09-17)
- Bump dash-renderer version from 1.0.1 to [1.1.0](https://github.com/plotly/dash/blob/master/dash-renderer/CHANGELOG.md#110---2019-09-17)

### Fixed
- [#915](https://github.com/plotly/dash/issues/915) Fix `dash-generate-components` on Windows.
- [#829](https://github.com/plotly/dash/issues/829) Fix the `--remote` pytest argument which was not effective in the code, adding a new argument `--remote-url` to support the selenium grid usage in the cloud.
- [#910](https://github.com/plotly/dash/pull/910) Reduce the dash-renderer packages size on **PyPI** about 55% by removing the source maps. To do more advanced debugging, the source maps needs to be generated from source code with `npm run build:local` and pip install in editable mode, i.e. `pip install -e .`

## [1.2.0] - 2019-08-27
### Added
- [#860](https://github.com/plotly/dash/pull/860) Add a new arg `dev_tools_prune_errors` to `app.run_server` and `app.enable_dev_tools`. Default `True`, tracebacks only include user code and below. Set it `False` for the previous behavior showing all the Dash and Flask parts of the stack.

### Changed
- Bump dash-table version from 4.1.0 to [4.2.0](https://github.com/plotly/dash-table/blob/master/CHANGELOG.md#420---2019-08-27)
- Bump dash-core-components version from 1.1.1 to [1.1.2](https://github.com/plotly/dash-core-components/blob/master/CHANGELOG.md#112---2019-08-27)
- Bump dash-html-components version from 1.0.0 to [1.0.1](https://github.com/plotly/dash-html-components/blob/master/CHANGELOG.md#101---2019-08-27)
- Bump dash-renderer version from 1.0.0 to [1.0.1](https://github.com/plotly/dash/blob/dev/dash-renderer/CHANGELOG.md#101---2019-08-27)

### Fixed
- [#874](https://github.com/plotly/dash/pull/874) Clean all the binary assets in dash-renderer, add tool to build all the required bundles from fresh source code to avoid confusion of the assets and improve the release process. Fixes [#868](https://github.com/plotly/dash/pull/868) and [#734](https://github.com/plotly/dash/pull/734)

## [1.1.1] - 2019-08-06
### Changed
- Bump dash-core-components version from 1.1.0 to [1.1.1](https://github.com/plotly/dash-core-components/blob/master/CHANGELOG.md#111---2019-08-06)

## [1.1.0] - 2019-08-05
### Added
- [#827](https://github.com/plotly/dash/pull/827) Add support for dashR testing to the `dash.testing` pytest framework.

### Changed
- Bump dash-table version from 4.0.2 to [4.1.0](https://github.com/plotly/dash-table/blob/master/CHANGELOG.md#410---2019-08-05)
- Bump dash-core-components version from 1.0.0 to [1.1.0](https://github.com/plotly/dash-core-components/blob/master/CHANGELOG.md#110---2019-08-05)

## [1.0.2] - 2019-07-15
### Changed
- Bump dash-table version from 4.0.1 to [4.0.2](https://github.com/plotly/dash-table/blob/master/CHANGELOG.md#402---2019-07-15)

### Fixed
- [#821](https://github.com/plotly/dash/pull/821) Fix a bug with callback error reporting, [#791](https://github.com/plotly/dash/issues/791).

## [1.0.1] - 2019-07-09
### Changed
- 💥 [#808](https://github.com/plotly/dash/pull/808) Remove strong `dash.testing` dependencies per community feedback. Testing users should do `pip install dash[testing]` afterwards.

- [#805](https://github.com/plotly/dash/pull/805) Add headless mode for dash.testing, add `pytest_setup_options` hook for full configuration of `WebDriver Options`.

- Bump dash-table version from 4.0.0 to [4.0.1](https://github.com/plotly/dash-table/blob/master/CHANGELOG.md#401---2019-07-09)

## [1.0.0] - 2019-06-20
### Changed
- 💥 [#761](https://github.com/plotly/dash/pull/761) Several breaking changes to the `dash.Dash` API:
  - Remove two obsolete constructor kwargs: `static_folder` and `components_cache_max_age`
  - Remove the misspelled `supress_callback_exceptions` fallback
  - Remove the unused `resources.config.infer_from_layout`
  - Revamp `app.config`: ALL constructor args are now stored in `config`, with three exceptions: `server`, `index_string`, and `plugins`. None of these are stored in any other instance attributes anymore.
  - Change `hot_reload_interval` from msec to seconds, for consistency with `hot_reload_watch_interval`
  - When called from `enable_dev_tools`, `debug=True` by default. It's still `False` by default from `run_server`.

- ✨ [#744](https://github.com/plotly/dash/pull/744) Introducing Dash Testing (`dash.testing`) - read the full tutorial at <https://dash.plotly.com/testing>.

- [#753](https://github.com/plotly/dash/pull/753) `Component` no longer inherits `MutableMapping`, so `values`, `keys`, and more are no longer methods. Fixes an issue reported in [dcc#440](https://github.com/plotly/dash-core-components/issues/440) where components with certain prop names defined but not provided would cause a failure to render. During component generation we now disallow all props with leading underscores or matching a few remaining reserved words: `UNDEFINED`, `REQUIRED`, `to_plotly_json`, `available_properties`, and `available_wildcard_properties`.

- [#739](https://github.com/plotly/dash/pull/739) Allow the Flask app to be provided to Dash after object initialization. This allows users to define Dash layouts etc when using the app factory pattern, or any other pattern that inhibits access to the app object. This broadly complies with the flask extension API, allowing Dash to be considered as a Flask extension where it needs to be.

- [#774](https://github.com/plotly/dash/pull/774) Allow the Flask app to set the Dash app name if the name is not provided by users.

- [#722](https://github.com/plotly/dash/pull/722) Assets are served locally by default. Both JS scripts and CSS files are affected. This improves robustness and flexibility in numerous situations, but in certain cases initial loading could be slowed. To restore the previous CDN serving, set `app.scripts.config.serve_locally = False` (and similarly with `app.css`, but this is generally less important).

- [#724](https://github.com/plotly/dash/pull/724), [renderer#175](https://github.com/plotly/dash-renderer/pull/175) Undo/redo toolbar is removed by default, you can enable it with `app=Dash(show_undo_redo=true)`. The CSS hack `._dash-undo-redo:{display:none;}` is no longer needed

- 💥 [#709](https://github.com/plotly/dash/pull/709) Merge the `dash-renderer` project into the main dash repo to simplify feature dev workflow. We will keep the [deprecated one](https://github.com/plotly/dash-renderer) for archive purpose.

## [0.43.0] - 2019-05-15
### Changed
- Bump dash-core-components version from 0.47.0 to [0.48.0](https://github.com/plotly/dash-core-components/blob/master/CHANGELOG.md#0480---2019-05-15)
- Bump dash-renderer version from 0.23.0 to [0.24.0](https://github.com/plotly/dash-renderer/blob/master/CHANGELOG.md#0240---2019-05-15)
- Bump dash-table version from 3.6.0 to [3.7.0](https://github.com/plotly/dash-table/blob/master/CHANGELOG.md#370---2019-05-15)

### Fixed
- [renderer#170](https://github.com/plotly/dash-renderer/pull/170) Fix regression on handling PreventUpdate (204 NO CONTENT)

## [0.42.0] - 2019-04-25
### Added
- [#687](https://github.com/plotly/dash/pull/687), [renderer#100](https://github.com/plotly/dash-renderer/pull/100) Dev Tools support. A new UI in the application that automatically display JavaScript & Python error messages, validates your component's properties, and displays a graph of your callback's dependencies. Only enabled in debug mode. Turn this on and off with two new config flags in `app.run_server`:
  - `dev_tools_props_check` - turn on/off property validation.
  - `dev_tools_ui` - turn on/off the UI.

### Fixed
- [renderer#148](https://github.com/plotly/dash-renderer/issues/148) Fix regression for `children=0` case.

## [0.41.0] - 2019-04-10
### Added
- [#672](https://github.com/plotly/dash/pull/672), [renderer#143](https://github.com/plotly/dash-renderer/pull/143) Support for "Clientside Callbacks" - an escape hatch to execute your callbacks in JavaScript instead of Python
- [#676](https://github.com/plotly/dash/pull/676) Add `dev_tools_ui` config flag in `app.run_server` (serialized in `<script id="_dash-config" type="application/json">`) to display or hide the forthcoming Dev Tools UI in Dash's front-end (dash-renderer).
- [#680](https://github.com/plotly/dash/pull/680) Partial updates: leave some multi-output updates unchanged while updating others

### Removed
- [renderer#145](https://github.com/plotly/dash-renderer/pull/145) Remove `dash_renderer._set_react_version` support for 15.4.2 and 16.2.0

### Changed
- Bump dash-core-components version from 0.45.0 to [0.46.0](https://github.com/plotly/dash-core-components/blob/master/CHANGELOG.md#0460---2019-04-10)
- [renderer#145](https://github.com/plotly/dash-renderer/pull/145) Update from React 15.4.2 to React 16.8.6

## [0.40.0] - 2019-03-25
### Changed
- Bump dash-core-components version from 0.44.0 to [0.45.0](https://github.com/plotly/dash-core-components/blob/master/CHANGELOG.md#0450---2019-03-25)
- Bump dash-html-components version from 0.14.0 to [0.15.0](https://github.com/plotly/dash-html-components/blob/master/CHANGELOG.md#0150---2019-03-25)
- [renderer#140](https://github.com/plotly/dash-renderer/pull/140), [renderer#126](https://github.com/plotly/dash-renderer/pull/126) Optimize rendering, and always assign `setProps` to components even with no callbacks to use it.

## [0.39.0] - 2019-03-04
### Added
- [#436](https://github.com/plotly/dash/pull/436) Allow multiple outputs from a single callback.
- [#367](https://github.com/plotly/dash/pull/367) Support custom JavaScript hooks to modify callback payloads and responses.
- [#623](https://github.com/plotly/dash/pull/623) Modify the flask response with custom cookies or headers, using `dash.callback_context.response`.
- [renderer#93](https://github.com/plotly/dash-renderer/pull/93) Loading states API

### Changed
- Bump dash-core-components version from 0.43.1 to [0.44.0](https://github.com/plotly/dash-core-components/blob/master/CHANGELOG.md#0440---2019-03-04)
- Bump dash-html-components version from 0.13.5 to [0.14.0](https://github.com/plotly/dash-html-components/blob/master/CHANGELOG.md#0140---2019-03-04)
- Bump dash-table version from 3.5.0 to [3.6.0](https://github.com/plotly/dash-table/blob/master/CHANGELOG.md#360---2019-03-04)

## [0.38.0] - 2019-02-25
### Added
- [#603](https://github.com/plotly/dash/pull/603) Add components libraries js/css distribution to hot reload watch.
- [#608](https://github.com/plotly/dash/pull/608), [renderer#124](https://github.com/plotly/dash-renderer/pull/124) Callback context:
  - Know which inputs caused a callback to fire: `dash.callback_context.triggered`
  - Input/State values by name `dash.callback_context.states.get('btn.n_clicks')`

### Changed
- Bump dash-table version from 3.4.0 to [3.5.0](https://github.com/plotly/dash-table/blob/master/CHANGELOG.md#350---2019-02-25)
- Bump dash-renderer version from 0.18.0 to [0.19.0](https://github.com/plotly/dash-renderer/blob/master/CHANGELOG.md#0190---2019-02-25)

### Fixed
- Fix missing indentation for generated metadata.json [#600](https://github.com/plotly/dash/issues/600)
- Fix missing component prop docstring error [#598](https://github.com/plotly/dash/issues/598)
- [#492](https://github.com/plotly/dash/pull/492) Move `__repr__` to base component instead of being generated.
- [#605](https://github.com/plotly/dash/pull/605) Raise exception when same input & output are used in a callback

## [0.37.0] - 2019-02-11
### Removed
- [renderer#118](https://github.com/plotly/dash-renderer/pull/118) Removed redux logger for the dev.

### Changed
- [#565](https://github.com/plotly/dash/pull/565) Add core libraries as version locked dependencies
- Bump dash-table version from 3.3.0 to [3.4.0](https://github.com/plotly/dash-table/blob/master/CHANGELOG.md#340---2019-02-08)
- Bump dash-renderer version from 0.17.0 to [0.18.0](https://github.com/plotly/dash-renderer/blob/master/CHANGELOG.md#0180---2019-02-11)
- Bump dash-core-components version from 0.43.0 to [0.43.1](https://github.com/plotly/dash-core-components/blob/master/CHANGELOG.md#0431---2019-02-11)

### Fixed
- [#563](https://github.com/plotly/dash/pull/563) Fix collections.abc deprecation warning for Python 3.8

## [0.36.0] - 2019-01-25
### Removed
- [#550](https://github.com/plotly/dash/pull/550), [renderer#114](https://github.com/plotly/dash-renderer/pull/114) Remove support for `Event` system. Use event properties instead, for example the `n_clicks` property instead of the `click` event, see [#531](https://github.com/plotly/dash/issues/531). `dash_renderer` MUST be upgraded to >=0.17.0 together with this, and it is recommended to update `dash_core_components` to >=0.43.0 and `dash_html_components` to >=0.14.0.

## [0.35.3] - 2019-01-23
### Changed
- [#547](https://github.com/plotly/dash/pull/547)
  - `assets_folder` argument now defaults to 'assets'
  - The assets folder is now always relative to the given root path of `name` argument, the default of `__main__` will get the `cwd`.
  - No longer coerce the name argument from the server if the server argument is provided.

### Fixed
- [#547](https://github.com/plotly/dash/pull/547)
  - Asset blueprint takes routes prefix into it's static path.
  - Asset url path no longer strip routes from requests.
- [#548](https://github.com/plotly/dash/pull/548) Remove print statement from PreventUpdate error handler.
- [#524](https://github.com/plotly/dash/pull/524) Removed ComponentRegistry dist cache.

## [0.35.2] - 2019-01-11
### Fixed
- [#522](https://github.com/plotly/dash/pull/522) Fix typo in some exception names
- [renderer#110](https://github.com/plotly/dash-renderer/pull/110)
  - Keep the config store state on soft reload.
  - AppProvider returns `Loading...` if no configs as before [renderer#108](https://github.com/plotly/dash-renderer/pull/108).

## 0.35.1 - 2018-12-27
### Fixed
- [#518](https://github.com/plotly/dash/pull/518) Always skip `dynamic` resources from index resources collection.

## 0.35.0 - 2018-12-18
### Added
- [#483](https://github.com/plotly/dash/pull/483) Experimental `--r-prefix` option to `dash-generate-components`, optionally generates R version of components and corresponding R package.

## 0.34.0 - 2018-12-17
### Added
- [#490](https://github.com/plotly/dash/pull/490) Add `--ignore` option to `dash-generate-components`, defaults to `^_`.

### Removed
- [renderer#108](https://github.com/plotly/dash-renderer/pull/108) Unused login api and Authentication component

### Fixed
- Add `key` to rendered components, fixing [renderer#379](https://github.com/plotly/dash-core-components/issues/379)

## 0.33.0 - 2018-12-10
### Added
- [#487](https://github.com/plotly/dash/pull/487) Add specific Dash exception types to replace generic exceptions (`InvalidIndexException`, `DependencyException`, `ResourceException`)

## 0.32.2 - 2018-12-09
### Fixed
- [#485](https://github.com/plotly/dash/pull/485) Fix typo in missing events/inputs error message

## 0.32.1 - 2018-12-07
### Changed
- [#484](https://github.com/plotly/dash/pull/484) Mute dash related missing props docstring from extract-meta warnings

## 0.32.0 - 2018-12-07
### Added
- [#478](https://github.com/plotly/dash/pull/478), [renderer#104](https://github.com/plotly/dash-renderer/issues/104) Support for .map file extension and dynamic (on demand) loading
- [renderer#107](https://github.com/plotly/dash-renderer/pull/107) [Redux devtools](https://github.com/zalmoxisus/redux-devtools-extension) support

## 0.31.1 - 2018-11-29
### Fixed
- [#473](https://github.com/plotly/dash/pull/473) Fix `_imports_.py` indentation generation.

## 0.31.0 - 2018-11-29
### Added
- [#451](https://github.com/plotly/dash/pull/451) Combine `extract-meta` and Python component files generation in a cli

### Fixed
- Fix a bug [renderer#66](https://github.com/plotly/dash-renderer/issues/66) in the ON_PROP_CHANGE callback where history was not correctly set when acting on more than one component. In particular, the 'undo' button should now work as expected.

## 0.30.0 - 2018-11-14
### Added
- [#362](https://github.com/plotly/dash/pull/362), [renderer#73](https://github.com/plotly/dash-renderer/pull/73) Hot reloading from the browser.
- Silence routes logging with `dev_tools_silence_routes_logging`.

## 0.29.0 - 2018-11-06
### Added
- [#444](https://github.com/plotly/dash/pull/444) Add component namespaces registry, collect the resources needed by component library when they are imported instead of crawling the layout.

## 0.28.7 - 2018-11-05
### Fixed
- [#450](https://github.com/plotly/dash/pull/450) Use the same prop name black list for component generation in all supported Python versions. Closes [#361](https://github.com/plotly/dash/issues/361).

## 0.28.6 - 2018-11-05
### Fixed
- [#443](https://github.com/plotly/dash/pull/443) `Dash.registered_paths` changed to a `collections.defaultdict(set)`, was appending the same package paths on every index.

## 0.28.5 - 2018-10-18
### Fixed
- [#431](https://github.com/plotly/dash/pull/431) Replace windows endline when generating components class docstrings.

## 0.28.4 - 2018-10-18
### Fixed
- [#430](https://github.com/plotly/dash/pull/430) Fix `Component.traverse()` and `Component.traverse_with_paths()` for components with `children` of type `tuple`, not just `list`.

## 0.28.3 - 2018-10-17
### Fixed
- [#418](https://github.com/plotly/dash/pull/418) Fix http-equiv typo
- Include missing polyfills to restore Internet Explorer support, restore whatwg-fetch [renderer#87](https://github.com/plotly/dash-renderer/issues/87)

## 0.28.2 - 2018-10-05
### Changed
- [#377](https://github.com/plotly/dash/pull/377) Move `add_url` function definition out of `Dash.__init__`

## 0.28.1 - 2018-09-26
### Fixed
- [#407](https://github.com/plotly/dash/pull/407) Missing favicon package_data from setup.py

## 0.28.0 - 2018-09-26
### Added
- [#406](https://github.com/plotly/dash/pull/406) Default favicon for dash apps.
- Bust the cache of the assets favicon.

### Fixed
- [#403](https://github.com/plotly/dash/pull/403) Remove the first and last blank lines from the HTML index string.

## 0.27.0 - 2018-09-20
### Added
- [#369](https://github.com/plotly/dash/pull/369), [renderer#77](https://github.com/plotly/dash-renderer/pull/77) Allow serving dev bundles from the components suite, enable with `app.run_server(dev_tools_serve_dev_bundles=True)`

### Fixed
- [#350](https://github.com/plotly/dash/pull/350) Use HTML5 syntax for the meta tag

## 0.26.6 - 2018-09-19
### Fixed
- [#387](https://github.com/plotly/dash/pull/387) Add `Cache-Control` headers to files served by `Dash.serve_component_suites`, and time modified query string to collected components suites resources.
- [#394](https://github.com/plotly/dash/pull/394) Add `InvalidResourceError` error and a Flask error handler so unregistered paths in `serve_component_suites` return a 404 instead of 500.

## 0.26.5 - 2018-09-10
### Fixed
- [#374](https://github.com/plotly/dash/pull/374) Fix `get_asset_url` with a different `assets_url_path`.

## 0.26.4 - 2018-08-28
### Fixed
- Set `url_base_pathname` to `None` in `Dash.__init__`. Fix [#364](https://github.com/plotly/dash/issues/364)

## 0.26.3 - 2018-08-27
### Added
- `Dash.get_asset_url` will give the prefixed url for the asset file.

### Fixed
- [#351](https://github.com/plotly/dash/pull/351) Prefix assets files with `requests_pathname_prefix`.

## 0.26.2 - 2018-08-26
### Fixed
- [#343](https://github.com/plotly/dash/pull/343) Only create the assets blueprint once for apps that provide the same flask instance to multiple dash instances.

## 0.26.1 - 2018-08-26
### Fixed
- [#336](https://github.com/plotly/dash/pull/336) Fix bug in `_validate_layout` which would not let a user set `app.layout` to be a function that returns a layout [(fixes #334)](https://github.com/plotly/dash/issues/334).

## 0.26.0 - 2018-08-20
### Added
- [#318](https://github.com/plotly/dash/pull/318) Add `assets_ignore` init keyword, regex filter for the assets files.

## 0.25.1 - 2018-08-20
### Fixed
- [#335](https://github.com/plotly/dash/pull/335) Ensure CSS/JS external resources are loaded before the assets.

## 0.25.0 - 2018-08-14
### Added
- [#322](https://github.com/plotly/dash/pull/322) Take config values from init or environ variables (Prefixed with `DASH_`).

### Fixed
- Take `requests_pathname_prefix` config when creating scripts tags.
- `requests/routes_pathname_prefix` must start and end with `/`.
- `requests_pathname_prefix` must end with `routes_pathname_prefix`. If you supplied both `requests` and `routes` pathname before this update, make sure `requests_pathname_prefix` ends with the same value as `routes_pathname_prefix`.
- `url_base_pathname` sets both `requests/routes` pathname, cannot supply it with either `requests` or `routes` pathname prefixes.

## 0.24.2 - 2018-08-13
### Fixed
- [#320](https://github.com/plotly/dash/pull/320) Disallow duplicate component ids in the initial layout.

## 0.24.1 - 2018-08-10
### Fixed
- Fix bug [#321](https://github.com/plotly/dash/issues/321) where importing Dash components with no props would result in an error.
- Fix a bug in 0.23.1 where importing components with arguments that are Python keywords could cause an error. In particular, this fixes `dash-html-components` with Python 3.7.

## 0.24.0 - 2018-08-10
### Added
- [#319](https://github.com/plotly/dash/pull/309) Add a modified time query string to assets included in the index in order to bust the cache.

## 0.23.1 - 2018-08-02
### Added
- [#316](https://github.com/plotly/dash/pull/316) Add `ie-compat` meta tag to the index by default.
- [#305](https://github.com/plotly/dash/pull/305) Add `external_script` and `external_css` keywords to dash `__init__`.
- Dash components are now generated at build-time and then imported rather than generated when a module is imported. This should reduce the time it takes to import Dash component libraries, and makes Dash compatible with IDEs.

## 0.22.1 - 2018-08-01
### Fixed
- [#273](https://github.com/plotly/dash/pull/273) Raise a more informative error if a non-JSON-serializable value is returned from a callback.

## 0.22.0 - 2018-07-25
### Added
- [#286](https://github.com/plotly/dash/pull/286) Asset files & index customization.
- [#294](https://github.com/plotly/dash/pull/294) Raise an error if there is no layout present when the server is run.
- [renderer#55](https://github.com/plotly/dash-renderer/pull/55) Add `_dash-error` class to the "Error loading layout" and "Error loading dependencies" messages.

### Fixed
- Attempting to render a `Boolean` value to the page no longer crashes the app.
- [renderer#57](https://github.com/plotly/dash-renderer/issues/57) If a callback references an `id` which does not exist in the DOM tree at the time it is executed, throw a more informative front-end exception.
- [renderer#54](https://github.com/plotly/dash-renderer/pull/54) Previously, if a component called `updateProps` with multiple properties, Dash would fire the callback multiple times (once for each property). Now the callback only fires once.

## 0.21.1 - 2018-04-10
### Added
- [#237](https://github.com/plotly/dash/pull/237) Support `aria-*` and `data-*` attributes in all dash html components. These new keywords can be added using a dictionary expansion, e.g. `html.Div(id="my-div", **{"data-toggle": "toggled", "aria-toggled": "true"})`
- [renderer#45](https://github.com/plotly/dash-renderer/pull/45) Allow user to choose between React versions '15.4.2' and '16.2.0':
```python
import dash_renderer

# Set the react version before setting up the Dash application
dash_renderer._set_react_version('16.2.0')

app = dash.Dash(...)
```

### Fixed
- [renderer#50](https://github.com/plotly/dash-renderer/pull/50) Update MANIFEST.in to include `react` and `react-dom` bundles for development mode

## 0.21.0 - 2018-02-21
### Added
- [#207](https://github.com/plotly/dash/pull/207) Support React components using [Flow](https://flow.org/en/docs/react/) types. `component_loader` now has the following behavior to create docstrings as determined in discussion in [#187](https://github.com/plotly/dash/issues/187):
  1. If a Dash component has `PropTypes`-generated typing, the docstring uses the `PropTypes`, _regardless of whether the component also has Flow types (current behavior)._
  2. Otherwise if a Dash component has Flow types but _not `PropTypes`_, the docstring now uses the objects generated by `react-docgen` from the Flow types.

### Fixed
- [renderer#42](https://github.com/plotly/dash-renderer/pull/42) Fix [renderer#41](https://github.com/plotly/dash-renderer/issues/41) and [renderer#44](https://github.com/plotly/dash-renderer/issues/44).
  - In some cases, during initialization, callbacks may fired multiple times instead of just once. This only happens in certain scenarios where outputs have overlapping inputs and those inputs are leaves (they don't have any inputs of their own).
  - If an output component is returned from a callback and its inputs were _not_ returned from the same input (i.e. they were already visible), then the callback to update the output would not fire. This has now been fixed. A common scenario where this app structure exists is within a Tabbed app, where there are global controls that update each tab's contents and the tab's callback just displays new output containers.

## 0.20.0 - 2018-01-19
### Added
- [#190](https://github.com/plotly/dash/pull/190) `exceptions.PreventUpdate` can be raised inside a callback to prevent the callback from updating the app. See <https://community.plotly.com/t/improving-handling-of-aborted-callbacks/7536/2>.

### Removed
- Removes logging from redux middleware from production build based on process.env.NODE_ENV.

### Changed
- Many pylint style fixes: [#163](https://github.com/plotly/dash/pull/163), [#164](https://github.com/plotly/dash/pull/164), [#165](https://github.com/plotly/dash/pull/165), [#166](https://github.com/plotly/dash/pull/166), [#167](https://github.com/plotly/dash/pull/167), [#168](https://github.com/plotly/dash/pull/168), [#169](https://github.com/plotly/dash/pull/169), [#172](https://github.com/plotly/dash/pull/172), [#173](https://github.com/plotly/dash/pull/173), [#181](https://github.com/plotly/dash/pull/181), [#185](https://github.com/plotly/dash/pull/185), [#186](https://github.com/plotly/dash/pull/186), [#193](https://github.com/plotly/dash/pull/193)
- [#184](https://github.com/plotly/dash/pull/184) New integration test framework.
- [#174](https://github.com/plotly/dash/pull/174) Submodules are now imported into the `dash` namespace for better IDE completion.

## 0.19.0 - 2017-10-16
### Changed
- 🔒 Remove CSRF protection measures. CSRF-style attacks are not relevant to Dash apps. Dash's API uses `POST` requests with content type `application/json` which are not susceptible to unwanted requests from 3rd party sites. See [#141](https://github.com/plotly/dash/issues/141).
- 🔒 `app.server.secret_key` is no longer required since CSRF protection was removed. Setting `app.server.secret_key` was difficult to document and a very common source of confusion, so it's great that users won't get bitten by this anymore :tada:
- 🐞 [renderer#22](https://github.com/plotly/dash-renderer/pull/22), [renderer#28](https://github.com/plotly/dash-renderer/pull/28) Previously, old requests could override new requests if their response was longer than the new one. This caused subtle bugs when apps are deployed on multiple processes or threads with component callbacks that update at varying rates like urls. Originally reported in [#133](https://github.com/plotly/dash/issues/133). This fix should also improve performance when many updates happen at once as outdated requests will get dropped instead of updating the UI. Performance issue with the first PR reported in [renderer#27](https://github.com/plotly/dash-renderer/issues/27) and fixed in the second PR.
- [renderer#21](https://github.com/plotly/dash-renderer/pull/21) Fix an issue where a callback would be fired excessively. Previously, the callback would be called as many times as it had inputs. Now, it is called less.

## 0.18.3 - 2017-09-08
### Added
- `app.config` is now a `dict` instead of a class. You can set config variables with `app.config['suppress_callback_exceptions'] = True` now. The previous class-based syntax (e.g. `app.config.suppress_callback_exceptions`) has been maintained for backwards compatibility.
- 🐌 Experimental behaviour for a customizable "loading state". When a callback is in motion, Dash now appends a `<div class="_dash-loading-callback"/>` to the DOM. Users can style this element using custom CSS to display loading screen overlays. This feature is in alpha, we may remove it at any time.

### Fixed
- Fix a bug from 0.18.2 that removed the ability for dash to serve the app on any route besides `/`.
- Fix a bug from 0.18.0 with the new config variables when used in a multi-app setting, causing config to be shared across apps. Originally reported in <https://community.plotly.com/t/flask-endpoint-error/5691/7>
- Rename config setting `supress_callback_exceptions` to `suppress_callback_exceptions`. The original spelling is kept for backward compatibility.
- 🐞 (renderer) Fix a bug where Dash would fire updates for each parent of a grandchild node that shared the same grandparent. Originally reported in <https://community.plotly.com/t/specifying-dependency-tree-traversal/5080/5>
- 🐞 (renderer) Fix a bug where the document title that displays "Updating..." wouldn't change if the callback raised an Exception. Now it will be removed on any response, even a failure.

## 0.18.2 - 2017-09-07
### Added
- [#70](https://github.com/plotly/dash/pull/70) 🔧 Add an `endpoint` to each of the URLs to allow for multiple routes.

## 0.18.1 - 2017-09-07
### Fixed
- [#128](https://github.com/plotly/dash/pull/128) 🐛 If `app.layout` is a function, then it used to be called excessively. Now it is called just once on startup and just once on page load.

## 0.18.0 - 2017-09-07
### Changed
- 🔒 Remove the `/static/` folder and endpoint that is implicitly initialized by flask. This is too implicit for my comfort level: I worry that users will not be aware that their files in their `static` folder are accessible
- ⚡️ Remove all API calls to the Plotly API (<https://api.plotly.com/>), the authentication endpoints and decorators, and the associated `filename`, `sharing` and `app_url` arguments. This was never documented or officially supported. Authentication has been moved to the [`dash-auth` package](https://github.com/plotly/dash-auth).
- [#107](https://github.com/plotly/dash/pull/107) ✏️ Sort prop names in exception messages.

### Added
- 🔧 Add two new `config` variables: `routes_pathname_prefix` and `requests_pathname_prefix` to provide more flexibility for API routing when Dash apps are run behind proxy servers. `routes_pathname_prefix` is a prefix applied to the backend routes and `requests_pathname_prefix` prefixed in requests made by Dash's front-end. `dash-renderer==0.8.0rc3` uses these endpoints.
- [#112](https://github.com/plotly/dash/pull/112) 🔧 Add `id` to `KeyError` exceptions in components.

### Fixed
- ✏️ Fix a typo in an exception.
- 🔧 Replaced all illegal characters in environment variables.

### 🔧 Maintenance
- 📝 Update README.md
- ✅ Fix CircleCI tests. Note that the [`dash-renderer`](https://github.com/plotly/dash-renderer) contains the bulk of the integration tests.
- 💄 Flake8 fixes and tests (fixes [#99](https://github.com/plotly/dash/issues/99))
- ✨ Add this CHANGELOG.md.

## 0.17.3 - 2017-06-22
✨ This is the initial open-source release of Dash.
# How to make a new Dash back end

Dash as a framework purposely generally does as much of its work as it can in the front end, both as a way to maximize performance (limit the demand on servers, limit network latency) and as a way to make it easier to create new back ends. Nevertheless, the Python back end is generally the canonical implementation, both because it’s the first and because it lives in the same repo as the front end, https://github.com/plotly/dash

In order to make a new back end then, the primary goal is to replicate the functionality of the Python version. This doesn’t mean we need to match the Python syntax; each back end should implement features in a way that’s natural to users of that language, but all else equal keeping this syntax and naming close to the Python version will help with documentation and with cross-language usage, for example to allow Python users (who at least at first will outnumber users of other languages by a large margin) to help other users on the community forum.

## Dash Components

Fundamentally, Dash components are React components. Take a look at the contents of https://github.com/plotly/dash-core-components/tree/master/dash_core_components - in each component repo we generate assets that need to be served to the browser:
- one main JavaScript bundle (`dash_core_components.js`)
- maybe some sub-bundles to load asynchronously (`async-*.js`)
- normally a sourcemap for each bundle (`*.js.map`)

We also include a couple of files that are only used to generate the components:
- `metadata.json`: a structured description of each component and its react props and prop types
- `package-info.json`: a copy of the package.json file from the repo, mainly useful for grabbing the version number

Then we also have Python files. Most of these are generated directly from `metadata.json` - for example `Checklist.py` is a class definition corresponding to the Checklist React component, and if you look inside it you’ll see it inherits from the `dash.development.base_component.Component` class, it has a docstring listing all props and their types in Python notation (instead of “array of objects”, it says “list of dicts”), and it has a constructor that ensures you can only create it with appropriate props. Then there are some Python files that are NOT generated, but copied from https://github.com/plotly/dash-core-components/tree/master/dash_core_components_base - `__init__.py` collects all the component classes and explicitly lists all the browser assets in `_js_dist` and possibly `_css_dist`.

Each back end must have a way to generate its own wrappers for these React components. The Python wrappers are generated source code files as described above, created by [`_py_components_generation.py`](https://github.com/plotly/dash/blob/dev/dash/development/_py_components_generation.py) - other languages may choose to either generate files, or create precompiled objects of some sort, but the key requirements are:
- Provide a natural way for users of this language to create components as data structures, keeping in mind that components may be nested inside the children prop of some other components, and most props are optional.
- To the extent that we can help users with built-in documentation and IDE auto-completion, we should try to do that.
- When requested by the framework, the component must serialize to JSON. This looks like:
  `{"namespace": "dash_core_components", "type": "Checklist", "props": {...}}`
- When a component package is included in the program (ie, during the equivalent of the Python statement `import dash_core_components`) the package must make itself known to the framework, so that the framework knows to include the main JavaScript bundle in the HTML for the page and knows where to find all the other JS and related files. Notice how in `__init__.py` in `_js_dist` some files are marked `"async": True|"eager"|"lazy"` or `"dynamic": True`, but the main bundle `dash_core_components.js` has neither. This full complexity is not needed in a new back end, just know you can ignore `"async": "eager"` files and treat all of the others with flags as not to include in the initial HTML but to be served later if requested.
- Some packages also have CSS files that must be loaded for its components to look and function correctly. In Python we collect these in `_css_dist` and register them with Dash during package import, the same way we do with `_js_dist`.
- Note that we CANNOT wait until a component from the package has been instantiated to alert the framework that this package is in use. It’s very common for the initial page layout to not include components from all packages. This MUST happen earlier.
- The Python artifacts are generated in the same repo as the JavaScript source files. Currently this is also the case for R and Julia, but we’re moving away from that model. New back ends should copy these JavaScript bundles to a new repo, and generate whatever other artifacts they need in this new repo.

## Dash server routes

There is a relatively small set of routes (urls/url patterns) that a Dash server must be prepared to serve. A good way to get a feel for them is to load https://dash.plotly.com/ and look at the page structure (“Elements” tab in Chrome devtools) and the network requests that are made on page load and their responses (“Network” tab). You can see all of the route patterns if you look at the main [`dash.dash` file](https://github.com/plotly/dash/blob/dev/dash/dash.py) and search for `self._add_url` - plus the one `self.server.register_blueprint` call above it. These routes are:
- `""` and anything not caught by other routes listed below (see https://dash.plotly.com/urls): the “index” route, serving the HTML page skeleton. The Python back end implements a bunch of customization options for this, but for a first version these are unnecessary. See the template given in [`_default_index`](https://github.com/plotly/dash/blob/357f22167d40ef00c92ff165aa6df23c622799f6/dash/dash.py#L58-L74) for the basic structure and pieces that need to be included, most importantly `{%config}` that includes info like whether to display in-browser devtools, and `{%scripts}` and `{%renderer}` that load the necessary `<script>` tags and initialize the dash renderer.
- `_dash-component-suites/<package_name>/<path>`: serve the JavaScript and related files given by each component package, as described above. Note that we include an extra cache-busting portion in each filename. In the Python version this is the version number and unix timestamp of the file as reported by the filesystem
- `_dash-layout"`: Gives the JSON structure of the initial Dash components to render on the page (`app.layout` in Python)
- `_dash-dependencies`: A JSON array of all callback functions defined in the app.
- `_dash-update-component`: Handles all server-side callbacks. Responds in JSON. Note that `children` outputs can contain Dash components, which must be JSON encoded the same way as in `_dash-layout`.
- `_reload-hash`: once hot reloading is implemented and enabled (normally only during development), this route tells the page when something has changed on the server and the page should refresh.
- `_favicon.ico`: the favicon, the title bar icon for the page, normally the Plotly logo.
- `assets/<path>`: static files to use on the page, potentially nested. CSS and JS files in this directory should be included in the HTML index page, and any file in this directory should be served if requested by the page. In Python this is handled by the `register_blueprint` call.

## Choosing a web server

You'll need to choose an HTTP(S) server for the programming language that you're adding to the Dash framework. For example, here are the HTTP(S) servers currently uses for the existing known distributions of Dash:

- [Dash Python](https://github.com/plotly/dash) uses [Flask](https://github.com/pallets/flask/) (though there is also a [community-maintained version that uses Django](https://github.com/GibbsConsulting/django-plotly-dash))
- [Dash Julia](https://github.com/plotly/dash.jl) uses [HTTP.jl](https://github.com/JuliaWeb/HTTP.jl)
- [Dash.NET](https://github.com/plotly/dash.net) uses [Giraffe](https://github.com/giraffe-fsharp/Giraffe)
- [Dash.R](https://github.com/plotly/dashr) uses [Fiery](https://github.com/thomasp85/fiery)

## Dash App Objects

In Python, users first create a dash object:
```py
from dash import Dash
app = Dash(...)
```
Then they set the layout to some nested Dash components:
```py
app.layout = html.Div(...)
```
Then they add callbacks:
```py
@app.callback(Output(...), Input(...), Input(...), State(...), ...)
def my_callback(input1, input2, state):
    <do stuff>
    return my_output_value

app.clientside_callback(
    “<JavaScript function as a string or reference>”,
    Output(...), Input(...), Input(...), State(...), ...
)
```
And finally they run the server
```py
app.run_server(...)
```

Any new back end needs to provide all of that functionality, one way or another.

The `Dash()` constructor has lots of options. They’re all listed and documented [here]( https://github.com/plotly/dash/blob/357f22167d40ef00c92ff165aa6df23c622799f6/dash/dash.py#L113-L253) - some are Python-specific (`name`, `server`, `plugins`), others should eventually be replicated but many can be left out of the first proof-of-concept implementation.

Similarly the `app.run_server` method has a lot of options, listed [here](https://github.com/plotly/dash/blob/357f22167d40ef00c92ff165aa6df23c622799f6/dash/dash.py#L1596-L1671) - again some are Python-specific (`flask_run_options`) and others can be added over time. Notice that many of these are simply passed on to `self.enable_dev_tools` - that’s because in Python the flask `server.run` command (called at the end of `app.run_server`) is normally only used in development, in production a more powerful server such as gunicorn is used, but a user may still want to enable devtools using a production server. You’re the expert on the new back end language, only you know if such a split makes sense. We don't want to write our own web server framework, you should be able to choose an existing one. Ideally this server should be easy to install on all major operating systems / platforms where your language is used, and for production scalability it should be able to run multiple workers or otherwise make full use of a multi-core machine. If it has the ability in development mode to create friendly error messages with stack traces to help debugging, that's a plus, but if not we can probably build that ourselves. 

## Contributing your own backends

We're very happy to work with anyone who is interested in making their own back end for Dash! Please get in touch through the [community forum](https://community.plotly.com/c/dash/) if you have started such an endeavor. In our experience, adding a new language to the Dash framework is not a weekend project, though we are constantly working to make this easier. At the moment, please anticipate a few months of work (including writing documentation). So far, mature Dash back ends exist for:

- [Python](https://github.com/plotly/dash)
- [Julia](https://github.com/plotly/dash.jl)
- [F#/.NET](https://github.com/plotly/dash.net)
- [R](https://github.com/plotly/dashr)

MATLAB, Scala, and SAS are examples of other scientific programming languages that could one day benefit from Dash's low-code, front end framework for creating AI and data science apps.

Happy coding!
# Dash

[![CircleCI](https://img.shields.io/circleci/project/github/plotly/dash/master.svg)](https://circleci.com/gh/plotly/dash)
[![GitHub](https://img.shields.io/github/license/plotly/dash.svg?color=dark-green)](https://github.com/plotly/dash/blob/master/LICENSE)
[![PyPI](https://img.shields.io/pypi/v/dash.svg?color=dark-green)](https://pypi.org/project/dash/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/dash.svg?color=dark-green)](https://pypi.org/project/dash/)
[![GitHub commit activity](https://img.shields.io/github/commit-activity/y/plotly/dash.svg?color=dark-green)](https://github.com/plotly/dash/graphs/contributors)
[![LGTM Alerts](https://img.shields.io/lgtm/alerts/g/plotly/dash.svg)](https://lgtm.com/projects/g/plotly/dash/alerts)
[![LGTM Grade](https://img.shields.io/lgtm/grade/python/g/plotly/dash.svg)](https://lgtm.com/projects/g/plotly/dash/context:python)

#### *Dash is the most downloaded, trusted Python framework for building ML & data science web apps*.

Built on top of [Plotly.js](https://github.com/plotly/plotly.js), [React](https://reactjs.org/) and [Flask](https://palletsprojects.com/p/flask/), Dash ties modern UI elements like dropdowns, sliders, and graphs directly to your analytical Python code. Read [our tutorial](https://dash.plotly.com/getting-started) (proudly crafted ❤️ with Dash itself).

- [Docs](https://dash.plotly.com/getting-started): Create your first Dash app in under 5 minutes

- [dash.gallery](https://dash.gallery): Dash app gallery with Python & R code

### Dash App Examples

| Dash App | Description |
|--- | :---: |
|![Sample Dash App](https://user-images.githubusercontent.com/1280389/30086128-9bb4a28e-9267-11e7-8fe4-bbac7d53f2b0.gif) | Here’s a simple example of a Dash App that ties a Dropdown to a Plotly Graph. As the user selects a value in the Dropdown, the application code dynamically exports data from Google Finance into a Pandas DataFrame. This app was written in just **43** lines of code ([view the source](https://gist.github.com/chriddyp/3d2454905d8f01886d651f207e2419f0)). |
|![Crossfiltering Dash App](https://user-images.githubusercontent.com/1280389/30086123-97c58bde-9267-11e7-98a0-7f626de5199a.gif)|Dash app code is declarative and reactive, which makes it easy to build complex apps that contain many interactive elements. Here’s an example with 5 inputs, 3 outputs, and cross filtering. This app was composed in just 160 lines of code, all of which were Python.|
|![Dash App with Mapbox map showing walmart store openings](https://user-images.githubusercontent.com/1280389/30086299-768509d0-9268-11e7-8e6b-626ac9ca512c.gif)| Dash uses [Plotly.js](https://github.com/plotly/plotly.js) for charting. About 50 chart types are supported, including maps. |
|![Financial report](https://github.com/plotly/dash-docs/blob/516f80c417051406210b94ea23a6d3b6cd84d146/assets/images/gallery/dash-financial-report.gif)| Dash isn't just for dashboards. You have full control over the look and feel of your applications. Here's a Dash App that's styled to look like a PDF report. |

To learn more about Dash, read the [extensive announcement letter](https://medium.com/@plotlygraphs/introducing-dash-5ecf7191b503) or [jump in with the user guide](https://plotly.com/dash).

### Dash OSS & Dash Enterprise

With Dash Open Source, Dash apps run on your local laptop or workstation, but cannot be easily accessed by others in your organization.

Scale up with Dash Enterprise when your Dash app is ready for department or company-wide consumption. Or, launch your initiative with Dash Enterprise from the start to unlock developer productivity gains and hands-on acceleration from Plotly's team.

ML Ops Features: A one-stop shop for ML Ops: Horizontally scalable hosting, deployment, and authentication for your Dash apps. No IT or DevOps required.
- [**App manager**](https://plotly.com/dash/app-manager/) Deploy & manage Dash apps without needing IT or a DevOps team. App Manager gives you point & click control over all aspects of your Dash deployments.
- [**Kubernetes scaling**](https://plotly.com/dash/kubernetes/) Ensure high availability of Dash apps and scale horizontally with Dash Enterprise’s Kubernetes architecture. No IT or Helm required.
- [**No code auth**](https://plotly.com/dash/authentication/) Control Dash app access in a few clicks. Dash Enterprise supports LDAP, AD, PKI, Okta, SAML, OpenID Connect, OAuth, SSO, and simple email authentication.
- [**Job Queue**](https://plotly.com/dash/job-queue/) The Job Queue is the key to building scalable Dash apps. Move heavy computation from synchronous Dash callbacks to the Job Queue for asynchronous background processing.

Low-Code Features: Low-code Dash app capabilities that supercharge developer productivity.
- [**Design Kit**](https://plotly.com/dash/design-kit/) Design like a pro without writing a line of CSS. Easily arrange, style, brand, and customize your Dash apps.
- [**Snapshot Engine**](https://plotly.com/dash/snapshot-engine/) Save & share Dash app views as links or PDFs. Or, run a Python job through Dash and have Snapshot Engine email a report when the job is done.
- [**Dashboard Toolkit**](https://plotly.com/dash/toolkit/) Drag & drop layouts, chart editing, and crossfilter for your Dash apps.
- [**Embedding**](https://plotly.com/dash/embedding/) Natively embed Dash apps in an existing web application or website without the use of IFrames.

Enterprise AI Features: Everything that your data science team needs to rapidly deliver AI/ML research and business initiatives.
- [**AI App Marketplace**](https://plotly.com/dash/ai-and-ml-templates/) Dash Enterprise ships with dozens of Dash app templates for business problems where AI/ML is having the greatest impact.
- [**Big Data for Pything**](https://plotly.com/dash/big-data-for-python/) Connect to Python's most popular big data back ends: Dask, Databricks, NVIDIA RAPIDS, Snowflake, Postgres, Vaex, and more.
- [**GPU & Dask Acceleration**](https://plotly.com/dash/gpu-dask-acceleration/) Dash Enterprise puts Python’s most popular HPC stack for GPU and parallel CPU computing in the hands of business users.
- [**Data Science Workspaces**](https://plotly.com/dash/workspaces/) Be productive from Day 1. Write and execute Python, R, & Julia code from Dash Enterprise's onboard code editor.

See [https://plotly.com/contact-us/](https://plotly.com/contact-us/) to get in touch.

![image](https://images.prismic.io/plotly-marketing-website/493eec39-8467-4610-b9d0-d6ad3ea61423_Dash+Open+source%2BDash+enterprise2-01.jpg?auto=compress,format)
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

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at chris@plotly.com. The project team will review and investigate all complaints, and will respond in a way that it deems appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [https://contributor-covenant.org/version/1/4][version]

[homepage]: https://contributor-covenant.org
[version]: https://contributor-covenant.org/version/1/4/
# Contributor Guide

## Getting Started

Glad that you decided to make your contribution in Dash, to set up your development environment, run the following commands:

```bash
# in your working directory
$ git clone git@github.com:plotly/dash.git
$ cd dash
$ python3 -m venv .venv/dev
# activate the virtualenv
# on windows `.venv\dev\scripts\activate`
# on some linux / mac environments, use `.` instead of `source`
$ source .venv/dev/bin/activate
# install dash and dependencies
$ pip install -e .[testing,dev]  # in some shells you need \ to escape []
$ npm install
# this script will build the dash-core-components, dash-html-components, dash-table,
# and renderer bundles; this will build all bundles from source code in their
# respective directories. The only true source of npm version is defined
# in package.json for each package.
$ npm run build  # runs `renderer build` and `npm build` in dcc, html, table
# build and install components used in tests
$ npm run setup-tests.py # or npm run setup-tests.R
# you should see dash points to a local source repo
$ pip list | grep dash
```

### Dash-Renderer Beginner Guide

`Dash Renderer` began as a separate repository. It was merged into the main `Dash` repository as part of the 1.0 release. It is the common frontend for all Dash backends (**R** and **Python**), and manages React Component layout and backend event handling.

If you want to contribute or simply dig deeper into Dash, we encourage you to play and taste it. This is the most efficient way to learn and understand everything under the hood.

For contributors with a primarily **Python** or **R** background, this section might help you understand more details about developing and debugging in JavaScript world.

As of Dash 1.2, the renderer bundle and its peer dependencies can be packed and generated from the source code. The `dash-renderer\package.json` file is the one version of the truth for dash renderer version and npm dependencies. A build tool `renderer`, which is a tiny Python script installed by Dash as a command-line tool, has a few commands which can be run from within the `dash/dash-renderer` directory:

1. `renderer clean` deletes all the previously generated assets by this same tool.
2. `renderer npm` installs all the npm modules using this `package.json` files. Note that the `package-lock.json` file is the computed reference product for the versions defined with tilde(~) or caret(^) syntax in npm.
3. `renderer bundles` parses the locked version JSON, copies all the peer dependencies into dash_renderer folder, bundles the renderer assets, and generates an `__init__.py` to map all the resources. There are also a list of helpful `scripts` property defined in `package.json` you might need to do some handy tasks like linting, syntax format with prettier, etc.
4. `renderer digest` computes the content hash of each asset in `dash_renderer` folder, prints out the result in logs, and dumps into a JSON file `digest.json`. Use this when you have a doubt about the current assets in `dash_renderer`, and compare it with previous result in one shot by this command.
5. `renderer build` runs 1, 2, 3, 4 in sequence as a complete build process from scratch.
6. `renderer build local` runs the same order as in 5 and also generates source maps for debugging purposes.

When a change in renderer code doesn't reflect in your browser as expected, this could be: confused bundle generation, caching issue in a browser, Python package not in `editable` mode, etc. The new tool reduces the risk of bundle assets by adding the digest to help compare asset changes.

### Development of `dash-core-components`, `dash-html-components`, and `dash_table`

Specific details on making changes and contributing to `dcc`, `html`, and `dash_table` can be found within their respective sub-directories in the `components` directory. Once changes have been made in the specific directories, the `dash-update-components` command line tool can be used to update the build artifacts and dependencies of the respective packages within Dash. For example, if a change has been made to `dash-core-components`, use `dash-update-components "dash-core-components"` to move the build artifacts to Dash. By default, this is set to update `all` packages.

## Python 2 And 3 Compatibility

Writing Python 2/3 compatible code might be a challenging task for contributors used to working on one particular version, especially new learners who start directly with Python 3.

We use `python-future` as our tool to mainly write Python 3 code and make it back-compatible to Python 2.7 (the only Python 2 version Dash supports). Please refer to [this list of idioms](https://python-future.org/compatible_idioms.html "https://python-future.org/compatible_idioms.html") for more details on working with `python-future`.

## Git

Use the [GitHub flow](https://guides.github.com/introduction/flow/) when proposing contributions to this repository (i.e. create a feature branch and submit a PR against the default branch).

### Organize your commits

For pull request with notable file changes or a big feature development, we highly recommend to organize the commits in a logical manner, so it

- makes a code review experience much more pleasant
- facilitates a possible cherry picking with granular commits

*an intuitive [example](https://github.com/plotly/dash-core-components/pull/548) is worth a thousand words.*

#### Git Desktop

Git command veterans might argue that a simple terminal and a cherry switch keyboard is the most elegant solution. But in general, a desktop tool makes the task easier.

1. <https://www.gitkraken.com/git-client>
2. <https://desktop.github.com/>

### Emoji

Plotlyers love to use emoji as an effective communication medium for:

**Commit Messages**

Emojis make the commit messages :cherry_blossom:. If you have no idea about what to add ? Here is a nice [cheatsheet](https://gitmoji.carloscuesta.me/) and just be creative!

**Code Review Comments**

- :dancer: `:dancer:` - used to indicate you can merge! Equivalent to GitHub's :squirrel:
- :cow2: `:cow2:` cow tip - minor coding style or code flow point
- :tiger2: `:tiger2:` testing tiger - something needs more tests, or tests need to be improved
- :snake: `:snake:` security snake - known or suspected security flaw
- :goat: `:goat:` grammar goat
- :smile_cat: `:smile_cat:` happy cat - for bits of code that you really like!
- :dolls: `:dolls:` documentation dolls
- :pill: `:pill:` performance enhancing pill
- :hocho: `:hocho:` removal of large chunks of code (obsolete stuff, or feature removals)
- :bug: `:bug:` - a bug of some kind. 8 legged or 6. Sometimes poisonous.
- :camel: :palm_tree: `:camel:` `:palm_tree:` - The Don't Repeat Yourself (DRY) camel or palm tree.
- :space_invader: `:space_invader:` - Too much space or too little.
- :spaghetti: `:spaghetti:` - copy-pasta, used to signal code that was copy-pasted without being updated

### Coding Style

We use `flake8`, `pylint`, and [`black`](https://black.readthedocs.io/en/stable/) for linting. please refer to the relevant steps in `.circleci/config.yml`.

## Tests

Our tests use Google Chrome via Selenium. You will need to install [ChromeDriver](https://chromedriver.chromium.org/getting-started) matching the version of Chrome installed on your system. Here are some helpful tips for [Mac](https://www.kenst.com/2015/03/installing-chromedriver-on-mac-osx/) and [Windows](http://jonathansoma.com/lede/foundations-2018/classes/selenium/selenium-windows-install/).

We use [pytest](https://docs.pytest.org/en/latest/) as our test automation framework, plus [jest](https://jestjs.io/) for a few renderer unit tests. You can `npm run test` to run them all, but this command simply runs `pytest` with no arguments, then `cd dash-renderer && npm run test` for the renderer unit tests.

Most of the time, however, you will want to just run a few relevant tests and let CI run the whole suite. `pytest` lets you specify a directory or file to run tests from (eg `pytest tests/unit`) or a part of the test case name using `-k` - for example `pytest -k cbcx004` will run a single test, or `pytest -k cbcx` will run that whole file. See the [testing tutorial](https://dash.plotly.com/testing) to learn about the test case ID convention we use.

### Unit Tests

For simple API changes, please add adequate unit tests under `/tests/unit`

Note: *You might find out that we have more integration tests than unit tests in Dash. This doesn't mean unit tests are not important, the [test pyramid](https://martinfowler.com/articles/practical-test-pyramid.html) is still valid. Dash project has its unique traits which needs more integration coverage than typical software project, another reason was that dash was a quick prototype crafted by chris in a lovely montreal summer.*

### Integration Tests

We introduced the `dash.testing` feature in [Dash 1.0](https://community.plotly.com/t/announcing-dash-testing/24868). It makes writing a Dash integration test much easier. Please read the [tutorial](https://dash.plotly.com/testing) and add relevant integration tests with any new features or bug fixes.

## Financial Contributions

Dash, and many of Plotly's open source products, have been funded through direct sponsorship by companies. [Get in touch] about funding feature additions, consulting, or custom app development.

[Dash Core Components]: https://dash.plotly.com/dash-core-components
[Dash HTML Components]: https://github.com/plotly/dash-html-components
[write your own components]: https://dash.plotly.com/plugins
[Dash Component Boilerplate]: https://github.com/plotly/dash-component-boilerplate
[issues]: https://github.com/plotly/dash-core-components/issues
[GitHub flow]: https://guides.github.com/introduction/flow/
[semantic versioning]: https://semver.org/
[Dash Community Forum]: https://community.plotly.com/c/dash
[Get in touch]: https://plotly.com/products/consulting-and-oem
[Documentation]: https://github.com/orgs/plotly/projects/8
[Dash Docs]: https://github.com/plotly/dash-docs
*Start with a description of this PR. Then edit the list below to the items that make sense for your PR scope, and check off the boxes as you go!*

## Contributor Checklist

- [ ] I have broken down my PR scope into the following TODO tasks
   -  [ ] task 1
   -  [ ] task 2
- [ ] I have run the tests locally and they passed. (refer to testing section in [contributing](https://github.com/plotly/dash/blob/master/CONTRIBUTING.md))
- [ ] I have added tests, or extended existing tests, to cover any new features or bugs fixed in this PR

### optionals

- [ ] I have added entry in the `CHANGELOG.md`
- [ ] If this PR needs a follow-up in **dash docs**, **community thread**, I have mentioned the relevant URLS as follows
    -  [ ] this GitHub [#PR number]() updates the dash docs
    -  [ ] here is the show and tell thread in Plotly Dash community
---
name: Bug report
about: Create a report to help us improve
title: "[BUG]"
labels: 'Type: Bug'
assignees: ''

---

Thank you so much for helping improve the quality of Dash!

We do our best to catch bugs during the release process, but we rely on your help to find the ones that slip through.


**Describe your context**
Please provide us your environment, so we can easily reproduce the issue.

-  replace the result of `pip list | grep dash` below
```
dash                 0.42.0
dash-core-components 0.47.0
dash-html-components 0.16.0
dash-renderer        0.23.0
dash-table           3.6.0
```
-  if frontend related, tell us your Browser, Version and OS

    - OS: [e.g. iOS]
    - Browser [e.g. chrome, safari]
    - Version [e.g. 22]

**Describe the bug**

A clear and concise description of what the bug is.

**Expected behavior**

A clear and concise description of what you expected to happen.

**Screenshots**

If applicable, add screenshots or screen recording to help explain your problem.
---
name: Feature request
about: Suggest an idea for this project
title: "[Feature Request]"
labels: 'Type: Enhancement'
assignees: ''

---

Thanks so much for your interest in Dash!

Before posting an issue here, please check the Dash [community forum](https://community.plotly.com/c/dash) to see if the topic has already been discussed. The community forum is also great for implementation questions. When in doubt, please feel free to just post the issue here :)

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
# Change Log for dash-table
### NOTE: as of v2.0, changes in dash-table are all being recorded in the main dash changelog.
### This file is kept only for historical purposes.
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).


## [4.12.0] - 2021-07-09
### Fixed
- [#907](https://github.com/plotly/dash-table/pull/907)
  - Fix a bug where pagination did not work or was not visible. [#834](https://github.com/plotly/dash-table/issues/834)
  - Fix a bug where if you are on a page that no longer exists after the data is updated, no data is displayed. [#892](https://github.com/plotly/dash-table/issues/892)


### Added
- [#916](https://github.com/plotly/dash-table/pull/916)
  - Added `html` option to `markdown_options` prop. This enables the use of html tags in markdown text.

- [#545](https://github.com/plotly/dash-table/issues/545)
    - Case insensitive filtering
    - New props: `filter_options` - to control case of all filters, `columns.filter_options` - to control filter case for each column
    - New operators: `i=`, `ieq`, `i>=`, `ige`, `i>`, `igt`, `i<=`, `ile`, `i<`, `ilt`, `i!=`, `ine`, `icontains` - for case-insensitive filtering, `s=`, `seq`, `s>=`, `sge`, `s>`, `sgt`, `s<=`, `sle`, `s<`, `slt`, `s!=`, `sne`, `scontains` - to force case-sensitive filtering on case-insensitive columns

### Changed
- [#918](https://github.com/plotly/dash-core-components/pull/918) Updated all dependencies. In particular the `highlight.js` upgrade changes code highlighting in markdown: we have long used their "github" style, this has been updated to more closely match current github styles.
- [#901](https://github.com/plotly/dash-core-components/pull/901) Updated R package `dash-info.yaml` to regenerate example without attaching now-deprecated core component packages (`dashHtmlComponents`, `dashCoreComponents`, or `dashTable`).

## [4.11.3] - 2021-04-08
### Changed
- [#862](https://github.com/plotly/dash-table/pull/862) - update docstrings per https://github.com/plotly/dash/issues/1205
- [#878](https://github.com/plotly/dash-table/pull/878) - update build process to use Webpack 5 and other latest dependencies

## [4.11.2] - 2021-01-19
### Fixed
- [#854](https://github.com/plotly/dash-table/pull/854) - part of fixing dash import bug https://github.com/plotly/dash/issues/1143

## [4.11.1] - 2020-12-07
### Fixed
- [#844](https://github.com/plotly/dash-table/pull/844) Fix a bug where the table is using classes that are styled by Bootstrap

## [4.11.0] - 2020-10-29
### Fixed
- [#841](https://github.com/plotly/dash-table/pull/841)
  - Fix prop-types regression causing console errors in browser devtools
  - Fix syntax highlighting regression for Markdown cells
- [#842](https://github.com/plotly/dash-table/pull/842) Fix a regression introduced with [#722](https://github.com/plotly/dash-table/pull/722) causing the tooltips to be misaligned with respect to their parent cell and incompletely addressed in [#817](https://github.com/plotly/dash-table/pull/817)

### Added
- [#841](https://github.com/plotly/dash-table/pull/841) Add Julia syntax highlighting support for Markdown cells
- [#831](https://github.com/plotly/dash-table/pull/831) Add the `tooltip_header` prop and add nested prop `use_with` (with values: `header`, `data`, `both`) to the `tooltip` prop to configure header cell tooltips

## [4.10.1] - 2020-09-03
-Dash.jl Julia component generation

## [4.10.0] - 2020-08-25
### Added
- [#820](https://github.com/plotly/dash-table/pull/820) Add support for Dash.jl Julia built components

### Fixed
- [#817](https://github.com/plotly/dash-table/pull/817) Fix a regression introduced with [#722](https://github.com/plotly/dash-table/pull/722) causing the tooltips to be misaligned with respect to their parent cell
- [#818](https://github.com/plotly/dash-table/pull/818) Fix a regression causing copy/paste not to work when selecting a range of cells with Shift + mouse click
- [#819](https://github.com/plotly/dash-table/pull/819) Fix pagination `page_current` and `page_count` fields to accommodate larger numbers

## [4.9.0] - 2020-07-27
### Added
- [#808](https://github.com/plotly/dash-table/pull/808)Fix a regression introduced with [#787](https://github.com/plotly/dash-table/pull/787) making it impossible to open markdown links in the current tab.
    - Adds a new `markdown_options` property that supports:
        - `link_target` nested prop with values `_blank`, `_parent`, `_self`, `_top` or an arbitrary string (default: `_blank`)

### Fixed
- [#806](https://github.com/plotly/dash-table/pull/806) Fix a bug where fixed rows a misaligned after navigating or editing cells [#803](https://github.com/plotly/dash-table/issues/803)
- [#809](https://github.com/plotly/dash-table/pull/809) Fix a bug where a scrollbar flickers on table render [#801](https://github.com/plotly/dash-table/issues/801)

## [4.8.1] - 2020-06-19
### Fixed
- [#798](https://github.com/plotly/dash-table/pull/798) Fix a bug where headers are not aligned with columns after an update [#797](https://github.com/plotly/dash-table/issues/797)

## [4.8.0] - 2020-06-17
### Added
- [#787](https://github.com/plotly/dash-table/pull/787) Add `cell_selectable` property to allow/disallow cell selection

### Changed
- [#787](https://github.com/plotly/dash-table/pull/787)
    - Clicking on a link in a Markdown cell now requires a single click instead of two
    - Links in Markdown cells now open a new tab (target="_blank")

### Fixed
- [#785](https://github.com/plotly/dash-table/pull/785) Fix a bug where the table does not refresh correctly if a property was previously missing
- [#793](https://github.com/plotly/dash-table/pull/793)
    - Fix a bug where headers are not aligned with columns with fixed_rows [#777](https://github.com/plotly/dash-table/issues/777)
    - Fix a regression where headers don't scroll horizontally with fixed_rows [#780](https://github.com/plotly/dash-table/issues/780)

## [4.7.0] - 2020-05-05
### Added
- [#729](https://github.com/plotly/dash-table/pull/729) Improve conditional styling
    - `style_data_conditional`: Add support for `row_index` and `column_id` array of values
    - `style_header_conditional`: Add support for `header_index` and `column_id` array of values
    - `style_filter_conditional`: Add support for `column_id` array of values
    - `style_cell_conditional`: Add support for `column_id` array of values
    - `style_data_conditional`: Add new conditions `state: 'active'|'selected'` to customize selected and active cell styles

### Fixed
- [#722](https://github.com/plotly/dash-table/pull/722) Fix a bug where row height is misaligned when using fixed_columns and/or fixed_rows
- [#728](https://github.com/plotly/dash-table/pull/728) Fix copy/paste on readonly cells
- [#724](https://github.com/plotly/dash-table/pull/724) Fix `active_cell` docstring: clarify optional nature of the `row_id` nested prop
- [#732](https://github.com/plotly/dash-table/pull/732) Fix a bug where opening a dropdown scrolled the table down its last row
- [#731](https://github.com/plotly/dash-table/pull/731) Fix a bug where `data=None` and `columns=None` caused the table to throw an error
- [#766](https://github.com/plotly/dash-table/pull/766) Sanitize table `id` for stylesheet injection (fixes usage with Pattern-Matching callbacks)

## Changed
- [#758](https://github.com/plotly/dash-table/pull/758) Improve error message for invalid filter queries

## [4.6.2] - 2020-04-01
### Changed
- [#713](https://github.com/plotly/dash-table/pull/713) Update from React 16.8.6 to 16.13.0

## [4.6.1] - 2020-02-27
### Added
- [#711](https://github.com/plotly/dash-table/pull/711) Added R examples to package help

### Changed
- [#704](https://github.com/plotly/dash-table/pull/704) Renamed async modules with hyphen `-` instead of tilde `~`

## [4.6.0] - 2020-01-14
### Added
- [#606](https://github.com/plotly/dash-table/pull/606) Add markdown support for table cells. Cells will be rendered as markdown if the column `presentation` is specified as `markdown`.
    - Add highlight.js for syntax highlighting. If `window.hljs` is specified, that will be used for highlighting instead.

### Fixed
- [#670](https://github.com/plotly/dash-table/pull/670) Fix a bug where `derived_filter_query_structure` was not getting updated properly
- [#677](https://github.com/plotly/dash-table/pull/677) Fix a bug where the table fails to load when used inside an iframe with a sandbox attribute that only has allow-scripts
- [#665](https://github.com/plotly/dash-table/pull/665) Fix a bug in Firefox where the dropdown cells height is incorrect

## [4.5.1] - 2019-11-14
### Fixed
- [#637](https://github.com/plotly/dash-table/pull/637) Fix multiple issues
  - Fix IE11 compatibility issues and add ES5 compatibility and validation
  - Fix a bug with `loading_state` being handled incorrectly, causing the table to steal focus

## [4.5.0] - 2019-10-29
### Changed
- [#554](https://github.com/plotly/dash-table/pull/554) Async loading of `xlsx` library on export

## [4.4.1] - 2019-10-17
### Fixed
- [#618](https://github.com/plotly/dash-table/issues/618) Fix a bug with keyboard navigation not working correctly in certain circumstances when the table contains `readonly` columns.
- [#206](https://github.com/plotly/dash-table/issues/206) Fix a bug with copy/paste to and from column filters not working.
- [#561](https://github.com/plotly/dash-table/issues/561) Fix an incorrect React PureComponent usage causing warnings in DevTools.
- [#611](https://github.com/plotly/dash-table/issues/611) Fix a bug with copy/paste causing hidden columns to be removed from the table

## [4.4.0] - 2019-10-08
### Added
[#546](https://github.com/plotly/dash-table/issues/546)
- New prop `export_columns` that takes values `all` or `visible` (default). This prop controls the columns used during export

[#597](https://github.com/plotly/dash-table/issues/597)
- Add `is blank` unary operator. Returns true for `undefined`, `null` and `''`.

[#299](https://github.com/plotly/dash-table/issues/299)
- New prop `page_count` that sets the maximum number of pages that are
  accessible via the pagination menu when using backend pagination.

### Changed
[#598](https://github.com/plotly/dash-table/issues/598)
- Allow values with whitespaces in column filters

[#580](https://github.com/plotly/dash-table/issues/580)
- Change pagination menu button UI to use arrow icons instead of plain
  buttons
- Move pagination menu to bottom-right of table
- Include go-to-first and go-to-last buttons
- Include current-page and total-pages display in pagination menu
- Include input box for user to navigate directly to a page

### Fixed

[#460](https://github.com/plotly/dash-table/issues/460)
- The `datestartswith` relational operator now supports number comparison
- Fixed a bug where the implicit operator for columns was `equal` instead of the expected default for the column type

[#546](https://github.com/plotly/dash-table/issues/546)
- Visible columns are used correctly for both header and data rows

[#563](https://github.com/plotly/dash-table/issues/563)
- Fixed a bug where any string beginning with a relational operator was being interpreted as that operator being applied to the rest of the string (e.g., "lens" was interpreted as "<=ns")

[#591](https://github.com/plotly/dash-table/issues/591)
- Fixed row and column selection when multiple tables are present

[#600](https://github.com/plotly/dash-table/issues/600)
- Fixed reconciliation when validation default value is `0` (number)
- Apply reconciliation value when deleting cells, if possible

## [4.3.0] - 2019-09-17
### Added
[#566](https://github.com/plotly/dash-table/pull/566)
- Support persisting user edits when the component or the page is reloaded. New props are `persistence`, `persistence_type`, and `persisted_props`. Set `persistence` to a truthy value to enable, the other two modify persistence behavior. See [plotly/dash#903](https://github.com/plotly/dash/pull/903) for more details.

[#319](https://github.com/plotly/dash-table/issues/319)
- New 'loading_state' prop that contains information about which prop, if any, is being computed.
- Table no longer allows for editing while the `data` prop is loading.

### Fixed
[#578](https://github.com/plotly/dash-table/pull/578)
- Fix [#576](https://github.com/plotly/dash-table/issues/576), editing column names or deleting columns while other columns are hidden causing the hidden columns to be lost.
- Fix an unreported bug that clicking "Cancel" at the column name edit prompt would clear the name, rather than leaving it unchanged as it should.

[#569](https://github.com/plotly/dash-table/issues/569), [#544](https://github.com/plotly/dash-table/issues/544)
- Allow empty strings in all `filter_query` (e.g filter_query: '{colA} eq ""')

[#567](https://github.com/plotly/dash-table/issues/567)
- Add support for missing `border-radius` in style_** props
- Fix table's inner vs. outer container styling

[#18](https://github.com/plotly/dash-table/issues/18)
- Fix row selection vertical and horizontal alignment

[#103](https://github.com/plotly/dash-table/issues/103)
- Simplify usage for multi-line cells and ellipsis. The cell's content now inherits the value of
`white-space`, `overflow` and `text-overflow` from its parent, making it possible to style
multi-line & ellipsis with `style_data` and other style props.

[#583](https://github.com/plotly/dash-table/issues/583)
- Fix regression when editing the content of a cell in a scrolled virtualized table

[#539](https://github.com/plotly/dash-table/issues/539)
- Fix bug where boolean values are not showing up in the table

## [4.2.0] - 2019-08-27
### Added
[#317](https://github.com/plotly/dash-table/issues/317)
- New `column.selectable` nested prop that displays a selection checkbox or radio button in the column.
- New `column_selectable` prop to choose whether columns can be selected or not, and whether a single or
    multiple selections can be in effect at the same time.
- New `selected_columns` prop that contains the list of visible and hidden columns that are currently selected
- New `derived_viewport_selected_columns` that contains the list of visible columns that are currently selected
    This prop is read-only. Use `selected_columns` in callbacks instead.

### Fixed
[#533](https://github.com/plotly/dash-table/issues/533)
- Fixed problem clearing one column shifting everything to the left and
leaving the last column blank
- Add merge_duplicate_headers prop to correct `export_format: display` behavior.
[#549](https://github.com/plotly/dash-table/issues/549)
- Fixed renaming of single-row headers in the GUI

## [4.1.0] - 2019-08-05
### Added
[#314](https://github.com/plotly/dash-table/issues/314)
- New `column.hideable` flag that displays an "eye" action icon in the column
    Accepts a boolean, array of booleans, 'last' or 'first'. Clicking on the "eye" will add the column to the `hidden_columns` prop.
    `hidden_columns` can be added back through the Columns toggle menu whether they are hideable or not.
- New accepted values for `column.clearable`, `column.deletable` and `column.renamable`
    These props now also accept 'last' and 'first'.
    - 'last' will display the action only on the last row of the headers
    - 'first' will display the action only on the first row of the headers

[#313](https://github.com/plotly/dash-table/issues/313)
- Ability to export table as csv or xlsx file.

[#497](https://github.com/plotly/dash-table/pull/497)
- New `column.clearable` flag that displays a "eraser" action in the column
    Accepts a boolean or array of booleans for multi-line headers.
    Clicking a merged column's "eraser" will clear all related columns.

    - Clearing column(s) will remove the appropriate data props from each datum
    row of `data`.
    - Additionally clearing the column will reset the filter for the affected column(s)

[#318](https://github.com/plotly/dash-table/issues/318)
- Headers are included when copying from the table to different
tabs and elsewhere. They are ignored when copying from the table onto itself and
between two tables within the same tab.

### Changed
[#497](https://github.com/plotly/dash-table/pull/497)
- Like for clearing above, deleting through the "trash" action will also
reset the filter for the affected column(s)

### Fixed
[#524](https://github.com/plotly/dash-table/issues/524)
- Fixed readonly dropdown cells content (display label, not value)

[#259](https://github.com/plotly/dash-table/issues/259)
- Fixed columns `sticky` on Safari

[#491](https://github.com/plotly/dash-table/issues/491)
- Fixed inconsistent behaviors when editing cell headers

[#521](https://github.com/plotly/dash-table/pull/521)
- Fixed white line artifacts when rendering the table with browser zoom different from 100%

## [4.0.2] - 2019-07-15
### Fixed
[#489](https://github.com/plotly/dash-table/issues/489)
- Add `fill_width` prop to replace `content_style` prop removed in the [4.0 API rework](https://github.com/plotly/dash-table/pull/446)

## [4.0.1] - 2019-07-09
### Changed
[#488](https://github.com/plotly/dash-table/pull/488)
- Update table build for use as a library and make consistent with other Dash repos

## [4.0.0] - 2019-06-20
### Changed
[#446](https://github.com/plotly/dash-table/pull/446)
- Table API rework
#### NEW
    - `column.sort_as_null`: Allows sorting behavior customization.
        Accepts an array of string, number or booleans.

#### REMOVED
    - `column.clearable`: Allows clearing the value of a dropdown cell.
        Removed in favor of `dropdown_**` `clearable` nested property.
    - `column.hidden`: Allows hiding column
        Removed. Stay tuned by following https://github.com/plotly/dash-table/issues/314.
    - `column.options`
        Removed. Redundant with `dropdown`.
    - `content_style`
        Removed. Deemed unnecessary. NOTE - This was added back in 4.0.2 under the name the "fill_width" property name.
    - `pagination_settings`
        Replaced by two props `page_current` and `page_size`.

#### RENAMED
    - `column_static_tooltip`
        Renamed to `tooltip`.
    - `column_conditional_tooltips`
        Renamed to `tooltip_conditional`.
    - `filter`
        Renamed to `filter_query`.
    - `sort_type`
        Renamed to `sort_mode`.
    - `derived_filter_structure`
        Renamed to `derived_filter_query_structure`.

#### MODIFIED
    - `column.deletable`: Allows column deletion.
        Now accepts a boolean or an array of booleans (for multi-line headers).
        For example, if there are multiple headers and you want the second header row to be deletable, this would be `[False, True]`.
    - `column.editable_name`: Allows column renaming.
        Renamed to `column.renamable`
        Now accepts a boolean or an array of booleans (for multi-line headers).
        For example, if there are multiple headers and you want the second row's header's name to be editable, this would be `[False, True]`.
    - `column.id`
        Now accepts `string` only -- `number` column ids can be casted to string.
    - `n_fixed_columns`: Will fix columns to the left.
        Renamed to `fixed_columns`
        Now accepts an object { headers: boolean, data: number } instead of a number.
        { headers: true } determines the number of columns to fix automatically. For example, if the rows are selectable or deletable, { headers: true } would fix those columns automatically. If { headers: true, data: 2 }, it would fix the first two data columns in addition to the selectable and deletable if visible.
    - `n_fixed_rows`: Will fix rows to the top.
        Renamed to `fixed_rows`
        Now accepts an object { headers: boolean, data: number } instead of a number.
        { headers: true } determines the number of rows to fix automatically (i.e. if there are multiple headers, it will fix all of them as well as the filter row).
        { headers: true, data: 2} would fix all of the header rows as well as the first 2 data rows.
    -  `pagination_mode`
        Renamed to `page_action`.
        `'fe'` is now `'native'`, `'be'` is now `'custom'`, and `false` is now '`none'`
    -  `column_static_dropdown`
        Renamed to `dropdown`.
        Now an object with each entry referring to a Column ID. Each nested prop expects.
        `clearable` and `options`.
    - `column_conditional_dropdowns`
        Renamed to `dropdown_conditional`.
        `condition` changed to the same `if` nested prop used by styles.
        `dropdown` renamed to `options`.
    - `dropdown_properties`
        Renamed to `dropdown_data`.
        Matches the `data` structure.
    - `tooltips`
        Renamed to `tooltip_data`.
        Matches the `data` structure.
    - `filtering`
        Renamed to `filter_action`.
    - `sorting`
        Renamed to `sort_action`.
    - `sorting_treat_empty_string_as_none`
        Renamed to `sort_as_null`.
        Now accepts an array of string, number or booleans that can be ignored during sort.
        Table-level prop for the `column.sort_as_null` column nested prop.
    - `style_data_conditional`
        Renamed `filter` to `filter_query`.

### Added
[#320](https://github.com/plotly/dash-table/issues/320)
- Ability to conditionally format columns if editing is disabled.

[#456](https://github.com/plotly/dash-table/issues/456)
- Support for dash-table is now available for R users of Dash.

### Fixed
[#434](https://github.com/plotly/dash-table/issues/434)
- Fix CSS borders properties overwrite style_* borders properties.

[#435](https://github.com/plotly/dash-table/issues/435)
- selected_cells background color is set through styling pipeline / derivations.

## [3.7.0] - 2019-05-15
### Added
[#397](https://github.com/plotly/dash-table/pull/397), [#410](https://github.com/plotly/dash-table/pull/410)
- Improve filtering syntax and capabilities
    - new field syntax `{myField}`
    - short form by-column filter
        - implicit column and default operator based on column type
            - Text and Any columns default to `contains`
            - Numeric columns default to `eq`
            - Date columns default to `datestartswith`
        - implicit column (e.g `ne "value"` becomes `{my-column} ne "value"`)
    - new `contains` relational operator for strings
    - new `datestartswith` relational operator for dates
    - new `eq` behavior (will attempt to convert and compare numeric values if possible)
    - new readonly `derived_filter_structure` prop exposing the query structure in a programmatically friendlier way

[#412](https://github.com/plotly/dash-table/pull/412)
- Add support for row IDs, based on the `'id'` attribute of each row of `data`
    - IDs will not be displayed unless there is a column with `id='id'`
    - `active_cell`, `start_cell`, `end_cell`, and items in `selected_cells` contain row and column IDs: All are now dicts  `{'row', 'column', 'row_id' and 'column_id'}` rather than arrays `[row, column]`.
    - Added new props mirroring all existing row indices props:
        - `selected_row_ids` mirrors `selected_rows`
        - `derived_viewport_row_ids` mirrors `derived_viewport_indices`
        - `derived_virtual_row_ids` mirrors `derived_virtual_indices`
        - `derived_viewport_selected_row_ids` mirrors `derived_viewport_selected_rows`
        - `derived_virtual_selected_row_ids` mirrors `derived_virtual_selected_rows`

[#424](https://github.com/plotly/dash-table/pull/424)
- Customizable cell borders through `style_**` props
    - cell borders now no longer use `box-shadow` and use `border` instead
    - Supports CSS shorthands:
        border, border_bottom, border_left, border_right, border_top
    - style_** props will ignore the following CSS rules:
        border_bottom_color, border_bottom_left_radius, border_bottom_right_radius, border_bottom_style, border_bottom_width, border_collapse, border_color, border_corner_shape, border_image_source, border_image_width, border_left_color, border_left_style, border_left_width, border_right_color, border_right_style, border_right_width, border_spacing, border_style, border_top_color, border_top_left_radius, border_top_right_radius, border_top_style, border_top_width, border_width
    - Styles priority:
        1. Props priority in decreasing order
            style_data_conditional
            style_data
            style_filter_conditional
            style_filter
            style_header_conditional
            style_header
            style_cell_conditional
            style_cell
        2. Within each props, higher index rules win over lower index rules
        3. Previously applied styles of equal priority win over later ones (applied top to bottom, left to right)

### Changed
[#397](https://github.com/plotly/dash-table/pull/397)
- Rename `filtering_settings` to `filter`

[#417](https://github.com/plotly/dash-table/pull/417)
- Rename `sorting_settings` to `sort_by`

[#412](https://github.com/plotly/dash-table/pull/412)
- `active_cell` and `selected_cells` items are dicts `{'row', 'column', 'row_id' and 'column_id'}` instead of arrays `[row, column]`

## [3.6.0] - 2019-03-04
### Fixed
[#189](https://github.com/plotly/dash-table/issues/189)
- Added `format` nested prop to columns
    - Applied to columns with `type=numeric` (more to come)
    - Uses [d3-format](https://github.com/d3/d3-format) under the hood
    - `format.locale` for localization configuration
    - `format.prefix` for SI prefix configuration
    - `format.specifier` for formatting configuration
    - `format.separate_4digits` to configure grouping behavior for numbers with 4 digits or less
    - Python helpers (dash_table.FormatTemplate)
- Added `locale_format` prop to table (default localization configuration, merged with column.format.locale)

[#387](https://github.com/plotly/dash-core/issues/387)
- Fix filtering conditions using floats

## [3.5.0] - 2019-02-25
### Added
[#342](https://github.com/plotly/dash-core/issues/342)
- Added `column_type` condition to style `if`; allows applying styles based on the type of the column for props
    - `style_cell_conditional`
    - `style_data_conditional`
    - `style_filter_conditional`
    - `style_header_conditional`

### Fixed
[#347](https://github.com/plotly/dash-core/issues/347)
- Fixed table behavior when `filtering_settings` is updated through a callback

[#322](https://github.com/plotly/dash-core/issues/322)
- Added LICENSE file to Python distributable

[#342](https://github.com/plotly/dash-core/issues/342)
- Added already supported `filter` nested prop / condition to `style_data_conditional` props definition

## [3.4.0] - 2019-02-08
### Added
[#364](https://github.com/plotly/dash-table/pull/364)
- Added the `datetime` data type

### Changed
[#224](https://github.com/plotly/dash-table/issues/224)
- Added support for unquoted column id with
  - letters, numbers, [-+:.]
- Added support for single and double quoted column id with arbitrary name

### Fixed
[#365](https://github.com/plotly/dash-table/issues/365)
- Incorrect tooltip behavior if cell is in a fixed row or column

- Incorrect default value for `column_static_tooltip` changed from [] to {}

## [3.3.0] - 2019-02-01
### Added
[#307](https://github.com/plotly/dash-core/issues/307)
- Added tooltip_delay and tooltip_duration props to tweak table's tooltips display behavior
- Added tooltips, column_static_tooltip, column_conditional_tooltips to define the tooltip
applicable to a certain cell in the table with nested props delay and duration to override
table's default behavior

## [3.2.0] - 2019-01-25
### Added
[#297](https://github.com/plotly/dash-core/issues/297)
- Added column.validation nested prop to tweak coercion and validation behavior
    - allow_null (boolean): [numeric, text] Allow null/undefined/NaN value
    - default (any): [numeric, text] Default value to use on validation/coercion failure
- Added on user-initiated data change processing (column.on_change.action)
    - Coerce: As Validation + attempts to convert the user-provided value into the destination type
    - None: Accept the user-provided value without verification
    - Validation: Check if the user-provided value is of the correct type/format/etc.
- Added on user-initiated data change failure processing (column.on_change.failure)
    This comes into effect after the column.on_change.action has failed
    - Accept: Accept the user-provide value
    - Default Applies the value defined in column.validation.default
    - Reject: Confirms the failure
### Changed
[#297](https://github.com/plotly/dash-core/issues/297)
- Moved column.type `dropdown` to column.presentation=`dropdown`

## [3.1.12] - 2019-01-11
### Fixed
- Regression, misaligned header [#324](https://github.com/plotly/dash-core/issues/324)
### Maintenance
- Test with head of both Dash v0.x and Dash v1.x [#20](https://github.com/plotly/dash-core/issues/20)

## [3.1.11] - 2018-12-10
### Fixed
- Selection, navigation, copy from readonly cell [#276](https://github.com/plotly/dash-table/issues/276)

## [3.1.10] - 2018-12-10
### Removed
- Deprecated nested property 'displayed_pages' from 'pagination_settings' [#275](https://github.com/plotly/dash-table/issues/275)

## [3.1.9] - 2018-12-06
### Added
- Source map [#284](https://github.com/plotly/dash-table/issues/284)
    Related Dash issue [#480](https://github.com/plotly/dash/issues/480)
### Changed
- Refactoring in preparation for data types [#280](https://github.com/plotly/dash-table/issues/280)

## [3.1.8] - 2018-12-04
### Added
- Virtualization [#234](https://github.com/plotly/dash-table/issues/234)
### Fixed
- Linting correctly applied to all code [#254](https://github.com/plotly/dash-table/issues/254)
### Changed
- Update dependencies [#278](https://github.com/plotly/dash-table/pull/278)
- Update dependencies [#274](https://github.com/plotly/dash-table/pull/274)
- Update dependencies [#251](https://github.com/plotly/dash-table/pull/251)


## [3.1.7] - 2018-11-19
### Fixed
- Visual offset with vertical scroll [#216](https://github.com/plotly/dash-table/issues/216)

## [3.1.6] - 2018-11-15
### Changed
- Generate python components classes for IDE support [#243](https://github.com/plotly/dash-table/pull/243)

## [3.1.5] - 2018-11-09
### Fixed
- Fix python package regression [#235](https://github.com/plotly/dash-table/issues/235)

## [3.1.4] - 2018-11-06
### Added
- New derived props for `selected_rows` [#147](https://github.com/plotly/dash-table/issues/147)
- Package library the UMD way [#212](https://github.com/plotly/dash-table/issues/212)

## [3.1.3] - 2018-11-05
### Fixed
- Fix load in IE 11 [#217](https://github.com/plotly/dash-table/issues/217)

## [3.1.2] - 2018-11-02
### Fixed
The version in the package didn't get updated.

## [3.1.1] - 2018-11-02
### Fixed
The remote URL path for the bundle was incorrect.

## [3.1.0] - 2018-11-02
- 3.1.0 (Alpha) Release of the Dash Table

## [3.1.0-rc21] - 2018-11-02
### Fixed
- Fire submit when pressing enter key in `IsolatedInput`. Fixes issue [#194](https://github.com/plotly/dash-table/issues/194)

## [3.1.0-rc20] - 2018-11-01
### Fixed
- Fix performance degradation on load [#208](https://github.com/plotly/dash-table/pull/208) [#200](https://github.com/plotly/dash-table/pull/200) [#198](https://github.com/plotly/dash-table/issues/198)

## [3.1.0-rc19] - 2018-11-01
### Changed
- Change default styles [#193](https://github.com/plotly/dash-table/pull/193) [#150](https://github.com/plotly/dash-table/issues/150)
    - prop `content_style` defaults to 'grow' instead of 'fit'
    - prop `style_table` width nested property defaults to '100%' if not provided
- Change cell styling and filter display [#196](https://github.com/plotly/dash-table/pull/196) [#150](https://github.com/plotly/dash-table/issues/150)
    - uneditable cells can be clicked & navigated, the mouse cursor is the default one
    - filter inputs have a placeholder that is visible on hover and focus
    - first filter input placeholder is always visible

## [3.1.0-rc18] - 2018-10-31
### Changed
- Rename table component to DataTable [#187](https://github.com/plotly/dash-table/pull/187) [#154](https://github.com/plotly/dash-table/issues/154)

## [3.1.0-rc17] - 2018-10-31
### Fixed
- Fix install on Linux [#184](https://github.com/plotly/dash-table/pull/184) [#137](https://github.com/plotly/dash-table/issues/137)

## [3.1.0-rc16] - 2018-10-30
### Fixed
- Fix copy/paste behavior when copying rows larger than data [#180](https://github.com/plotly/dash-table/pull/180) [#142](https://github.com/plotly/dash-table/issues/142)

## [3.1.0-rc15] - 2018-10-30
### Changed
- Column 'editable' prop takes precedence over table 'editable' prop[#182](https://github.com/plotly/dash-table/pull/182) [#175](https://github.com/plotly/dash-table/issues/175)

## [3.1.0-rc14] - 2018-10-30
### Changed
- Rename sorting_settings columnId -> column_id [#183](https://github.com/plotly/dash-table/pull/183) [#171](https://github.com/plotly/dash-table/issues/171)

## [3.1.0-rc13] - 2018-10-30
### Changed
- Allow keyboard navigation on focused input [#172](https://github.com/plotly/dash-table/pull/172) [#141](https://github.com/plotly/dash-table/issues/141) [#143](https://github.com/plotly/dash-table/issues/143)

## [3.1.0-rc12] - 2018-10-29
### Changed
- Rename selected_cell -> selected_cells [#181](https://github.com/plotly/dash-table/pull/181) [#177](https://github.com/plotly/dash-table/issues/177)

## [3.1.0-rc11]
### Fixed
- Fix regressions linked to the style_as_list_view feature / prop [#179](https://github.com/plotly/dash-table/pull/179)

## [3.1.0-rc10]
### Changed
- Improved props docstrings [#163](https://github.com/plotly/dash-table/issues/163)

## [3.1.0-rc9]
### Changed
- Sort ascending on first click [#164](https://github.com/plotly/dash-table/pull/164) [#118](https://github.com/plotly/dash-table/issues/118)
    - Flips icons displayed so that they are pointing up on ascending and down on descending.

## [3.1.0-rc8]
### Changed
- Improve props typing [#158](https://github.com/plotly/dash-table/pulls/158) [#143](https://github.com/plotly/dash-table/issues/143)

## [3.1.0-rc7]
### Changed
- Make table id optional (generate random one if needed) [#155](https://github.com/plotly/dash-table/pulls/155) [#139](https://github.com/plotly/dash-table/issues/139)

## [3.1.0-rc6]
### Added
- Styling API refactoring

    - Remove column width / maxWidth / minWidth
    - Rename property table_style to css

    Cell: All table cells
    Data: All data driven cells (no operations, headers, filters)
    Filter: All basic filter cells
    Header: All header cells

    Priority
    Data: style_data_conditional > style_data > style_cell_conditional > style_cell
    Filter: style_filter_conditional > style_filter > style_cell_conditional > style_cell
    Header: style_header_conditional > style_header > style_cell_conditional > style_cell

    Merge Logic
    Only properties defined at a higher priority level will override properties
    defined at lower priority levels. For example if A is applied then B, A+B will be..

    A = {
        background_color: 'floralwhite',
        color: 'red',
        font_type: 'monospace',
        width: 100
    }

    B = {
        color: 'black',
        font_size: 22
    }

    A+B = {
        background_color: 'floralwhite', // from A, not overridden
        color: 'black', // from B, A overridden
        font_size: 22, // from B
        font_type: 'monospace', // from A
        width: 100 // from A
    }

    - Add new property style_table of form
        { ...CSSProperties }
    - Add new property style_cell of form
        { ...CSSProperties }
    - Add new property style_data of form
        { ...CSSProperties }
    - Add new property style_filter of form
        { ...CSSProperties }
    - Add new property style_header of form
        { ...CSSProperties }
    - Add new property style_cell_conditional of form
        [{
            if: { column_id: string | number },
            ...CSSProperties
        }]
    - Add new property style_data_conditional of form
        [{
            if: { column_id: string | number, filter: string, row_index: number | 'odd' | 'even' },
            ...CSSProperties
        }]
    - Add new property style_filter_conditional of form
        [{
            if: { column_id: string | number },
            ...CSSProperties
        }]
    - Add new property style_header_conditional of form
        [{
            if: { column_id: string | number, header_index: number | 'odd' | 'even' },
            ...CSSProperties
        }]
    - All CSSProperties are supported in kebab-case, camelCase and snake_case

### Changed
- Renaming 'dataframe' props to 'data'
    dataframe -> data
    dataframe_previous -> data_previous
    dataframe_timestamp -> data_timestamp
    derived_virtual_dataframe -> derived_virtual_data
    derived_viewport_dataframe -> derived_viewport_data

## [3.1.0-rc5]
### Fixed
- Tests and fixes for editable/readonly [#134](https://github.com/plotly/dash-table/pull/134) [#132](https://github.com/plotly/dash-table/issues/132)

## [3.1.0-rc4]
### Added
Columns width percentage and default (fit to content) support
- Added prop content_style that takes values 'fit' or 'grow' (Default='fit')
- Added width percentage support
- Modified default column behavior from fixed width to 'fit content'
- Modified width, min-width, max-width interaction on columns

Width percentage

    Columns can now accept '%' width, minWidth, maxWidth

    For the percentages to have meaning, the dash-table must be forced to have a width and the content of the dash-table must be forced to grow to fill the available space made available by the container (by default the table is only as big as it needs to be).

    To use percentage-based column widths, add:

    - content style
        content_style='grow'

    - table style (example)
        table_style=[{ selector: '.dash-spreadsheet', rule: 'width: 100%; max-width: 100%' }]

    - column with %-based width
        columns=[{
            id: 'column',
            width: '40%'
        }]

Default column width

    Columns now default to 'fit to content' when no width is defined
    Note: If pagination is used or the dataframe modified, the column width will be re-evaluated on each modification.

Interaction between width, min-width and max-width

    Column min-width and max-width do not default to width value is not defined.

## [3.1.0-rc3]
### Fixed
- Miscellaneous fixes for pagination, virtual df and viewport df [#112](https://github.com/plotly/dash-table/pull/112)

## [3.1.0-rc2]
### Changed
- 1 new internal facing derived/controlled property:
    columns: Columns -> columns: VisibleColumns
    Gets rid of conditional processing for hidden columns in the cell and header factories as well as in navigation/selection handlers

Clean up column offsets
- A bunch of offsets were introduced to the table in the previous development cycle (2.x -> 3.0). Turns out these offsets are neither useful or necessary


Validate compatibility of filtering, sorting, pagination

External facing classes and attributes
- (ATTRIBUTE) data-dash-column=<columnId>

    .dash-cell,
    .dash-header {
        &[data-dash-column='ticker'] {
            // styling
        }
    }

- (CLASS) dash-cell
- (CLASS) dash-header

- (CLASS) dash-delete-cell
- (CLASS) dash-delete-header
- (CLASS) dash-select-cell
- (CLASS) dash-select-header

- (CLASS) dash-cell-value

- (CLASS) dash-freeze-left
- (CLASS) dash-freeze-top
- (CLASS) dash-spreadsheet
- (CLASS) dash-spreadsheet-container
- (CLASS) dash-spreadsheet-inner

## [3.1.0-rc1]
### Added
Version 3.1 of the Dash-Table builds upon the 3.0 table and solidifies the external facing API of the table
- introducing the notion of derived properties
- virtual and viewport dataframe and indices for more flexibility
- code refactoring to simplify and improve the existing implementation / prepare for the future
- documentation of the API and table features
- additional e2e, integration and unit tests for a more mature development platform

Derived Properties
- Derived properties are new to 3.1
- They are readonly properties that represent a transform from multiple 'first-class' properties of the component.
- For example, `derived_viewport_dataframe` is a readonly view based on `(dataframe, filtering params, sorting params, pagination params) --> derived_viewport_dataframe`

Derived properties allow the component to expose complex state that can be useful for a Dash Server developer but without introducing dual states, a situation where multiple properties may represent the same state within the component, making it necessary to reconcile them on each prop update.

Virtual and Viewport Dataframe
- 4 new external facing derived properties and 4 internal facing controlled properties that represent:
    1. the filtered and sorted dataframe and the indices mapping
    2. the filtered, sorted and paginated dataframe and the indices mapping

    - `derived_viewport_dataframe`
    - `derived_viewport_indices`
    - `derived_virtual_dataframe`
    - `derived_virtual_indices`

    In the event where sorting, filtering or pagination is done on the Dash server, it is possible that some or all derived dataframes will be equal to the dataframe prop.

## [3.0.0-rc22]
### Fixed
- Fix regression for user select
    Sorting arrow will no longer highlight.

## [3.0.0-rc21]
### Changed
- Improve performance when the user clicks outside of the table [#104](https://github.com/plotly/dash-table/pull/104)
    Clicking outside of the table was setting the table's `is_focused` property.
    Setting component properties in Dash can be expensive: it can cause the
    entire app to re-render.
    Now, clicking outside the table will update the component more efficiently,
    prevent excessive application re-renders.

## [3.0.0-rc20]
### Fixed
- Fix incorrect border around table cells when not filled [#102](https://github.com/plotly/dash-table/pull/102) [#101](https://github.com/plotly/dash-table/issues/101)
    Table styling has been changed for frozen rows and columns. Default styling change from:

    - frozen rows: { height: 500px } to { height: fit-content, max-height: 500px }
    - frozen columns: { width: 500px } to { width: fit-content, max-width: 500px }

## [3.0.0-rc19]
### Fixed
- Fix dropdown position & behavior on scroll [#96](https://github.com/plotly/dash-table/issues/96)
    Limitation: The dropdown in fixed columns behaves differently from the dropdown in the non-fixed portion of the table. Because of layers of overflow & positioning, the dropdown does not show outside of the table is instead part of it. Opening the dropdown in bottom rows will require scrolling vs. displaying on top of the table.

## [3.0.0-rc18]
### Added
Basic Filtering & Preparation work for advaced filtering
- Additional filtering_type prop that can take value 'basic' (or eventually 'advanced')
    This prop defines whether the user is presented with the UI to filter by column or with complex expressions
    The default value is 'basic'

    Note: The filtering row counts against n_fixed_rows

- Additional filtering_types prop that takes an array of values with valid values 'basic' (and eventually 'advanced')
    This prop defines what type of filtering are available to the user
    The default value is ['basic']

    Note: This value needs to be consistent with `filtering_type`

## [3.0.0-rc17]
### Added
- Make buttons non-selectable [#105](https://github.com/plotly/dash-table/pull/105) [#91](https://github.com/plotly/dash-table/issues/91)

## [3.0.0-rc16]
### Fixed
- Fix keyboard navigation after copy [#90](https://github.com/plotly/dash-table/pull/90) [#49](https://github.com/plotly/dash-table/issues/49)

## [3.0.0-rc15]
### Fixed
- Fix global copy paste regression [#87](https://github.com/plotly/dash-table/pull/87) [#75](https://github.com/plotly/dash-table/issues/75)
- Fix data paste [#87](https://github.com/plotly/dash-table/pull/87) [#88](https://github.com/plotly/dash-table/issues/88)

## [3.0.0-rc14]
### Fixed
- Empty dropdown setting value regression fix [#85](https://github.com/plotly/dash-table/pull/85) [#83](https://github.com/plotly/dash-table/issues/83)

## RC13 - Modify click & sequential click behavior
### Changed
- Partial implementation of new click & sequential click behavior [#79](https://github.com/plotly/dash-table/pull/79) [#77](https://github.com/plotly/dash-table/issues/77)
    - First click selects the cell's content and will cause user input to override the cell content.
    - Second click into the cell will remove the selection and position the cursor accordingly.

## [3.0.0-rc12]
### Fixed
- Fix border style [#78](https://github.com/plotly/dash-table/pull/78) [#68](https://github.com/plotly/dash-table/issues/68)
- Fix selection [#78](https://github.com/plotly/dash-table/pull/78) [#73](https://github.com/plotly/dash-table/issues/73)
- Fix row offset [#78](https://github.com/plotly/dash-table/pull/78) [#76](https://github.com/plotly/dash-table/issues/76)

## [3.0.0-rc11]
### Fixed
- Fix copy/paste regression [#70](https://github.com/plotly/dash-table/pull/70) [#64](https://github.com/plotly/dash-table/issues/64)
- Fix click/blur regression [#70](https://github.com/plotly/dash-table/pull/70) [#65](https://github.com/plotly/dash-table/issues/65)
- Fix delete regression [#70](https://github.com/plotly/dash-table/pull/70) [#67](https://github.com/plotly/dash-table/issues/67)

## [3.0.0-rc10]
### Fixed
- Fix double click regression [#63](https://github.com/plotly/dash-table/pull/63)

## [3.0.0-rc9]
### Added
Treat empty strings as none
- sorting_treat_empty_string_as_none takes value True or False
    Overrides sorting default behavior to consider empty strings as a nully value.
    Note: This is a stopgag prop, full implementation of sorting overrides will most probably deprecate it.
    Default value is False.

## [3.0.0-rc8]
### Fixed
setProps bug fix
- Fixing initialization issue where the FE wrongly believes it's working in DEV mode

## [3.0.0-rc7]
### Added
Additional `sorting_type` prop that can take value 'multi' or 'single'
This prop defines whether the user can sort based on multiple columns or can only sort by one column at a time. The default value is 'single'.

## [3.0.0-rc6]
### Added
Consolidating virtualization, sorting, filtering
- First steps to make sorting work from both FE and BE *
- and consistent with Virtualization settings *

New Props
- sorting -> ['fe', 'be', true, false] (default: false) -- replaces `sortable` prop
- sorting_settings -> array of { field, ascending } -- replaces `sort` prop
- `virtual_dataframe` (READONLY)
- `virtual_dataframe_indices` (READONLY; not officially supported yet -- IN DEVELOPMENT)

virtual_dataframe vs. dataframe
- the virtual dataframe is the content of the viewport for the user (e.g. user has a 10k rows dataframe with FE/250 lines paging, on 1st page -> the `virtual_dataframe` contains items [0,250[ of the dataframe); the dataframe still contains 10k items
- 10k rows, no paging, sorting and filtering -> the virtual dataframe contains items visible in the viewport, in the visible order; the dataframe still contains 10k items
- if the user modifies a cell, the dataframe and the virtual_dataframe are updated with the new data

### Deprecated
- sortable
- sort
- dataframe behavior on sort (see below)

## [3.0.0-rc5]
### Added
New props for Conditional Style, Conditional Dropdown, Filter
- filtering -> ['fe', 'be', true, false] (default: false)
- filtering_settings -> AST query string (default: '')
- column_conditional_dropdowns
- column_static_dropdown
- column_conditional_styles
- column_static_style
- row_conditional_styles
- row_static_style

### Deprecated
- column style
- column options
- dropdown_properties prop

## [3.0.0-rc4]
### Added
Version 3.0 of the Dash-Table expands vastly on the capability of the 2.x table and provides features:
- visually freezing rows and/or columns
- filtering in either FE or BE, basic filtering UI
- sorting in either FE or BE, basic sorting UI
- pagination in either FE or BE, basic pagination UI
- performance optimizations
- basic coverage through e2e, integration and unit tests

Virtualization
- See v_be_page_usage.py and v_fe_page_usage.py for FE and BE usage scenarios.
- virtual_dataframe and virtual_dataframe_indices are exposed and expected to be *readonly*. Setting them from the BE will have no impact on the FE display.
FE Virtualization
- BE is not expected to update the dataframe when the virtualization settings are updated.

BE Virtualization
- BE is expected to update the dataframe when the virtualization settings are updated.

Freeze Top Rows (Limitations)
- the table styling is forced to { table-layout: fixed; width: 0 !important; } to ensure the frozen section and the rest of the table stay in sync (width-wise); this means that the width of the table is only driven by the width of the columns (default width is 200px)
- can't freeze rows and columns at the same time

Freeze Left Columns (Limitations)
- performance is highly impacted if the table is in a scrollable container as the frozen columns position has to be recalculated on each scroll event; impact is minimal up to 50-100 items and makes the table difficult to use with 250-500 items
- can't freeze rows and columns at the same time
- when using merged headers, make sure that the number of fixed columns respects the merged headers, otherwise there will be some unresolved visual bugs/artifacts
- rows are assumed to all have the same height

Deletable Columns (Limitations)
- there might be unintended side-effects if used with BE virtualization (the act of deleting a column / columns modifies the dataframe)

Performance Improvements
- Table now renders and navigates faster
- Typing in cell does not modify dataframe until focus is lost / edit is confirmed ("enter" or "tab)

### Deprecated
- prop "update_on_unfocus" has been removed
# Dash Table

An interactive `DataTable` for [Dash](https://dash.plotly.com/).

:point_right: [Documentation](https://dash.plotly.com/datatable)

## Quickstart

```
pip install dash
```

```python
from dash import Dash, dash_table
import pandas as pd

df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/solar.csv')

app = Dash(__name__)

app.layout = dash_table.DataTable(
    id='table',
    columns=[{"name": i, "id": i} for i in df.columns],
    data=df.to_dict('records'),
)

if __name__ == '__main__':
    app.run_server(debug=True)
```

![Interactive Dash DataTable](https://user-images.githubusercontent.com/1280389/47935912-67187080-deb2-11e8-8936-34b0c99b518f.png)

## Background

Dash DataTable is an interactive table component designed for viewing, editing, and exploring large datasets.

DataTable is rendered with standard, semantic HTML `<table/>` markup, which makes it accessible, responsive, and easy to style.

This component was written from scratch in React.js and Typescript specifically for the Dash community. Its API was designed to be ergonomic and its behavior is completely customizable through its properties.

DataTable was designed with a featureset that allows that Dash users to create complex, spreadsheet driven applications with no compromises. We're excited to continue to work with users and companies that [invest in DataTable's future](https://plotly.com/products/consulting-and-oem/).

Note: DataTable is currently supported in Chrome, Firefox, Safari, Edge (version 15+), and Internet Explorer 11.

Share your DataTable Dash apps on the [community forum](https://community.plotly.com/t/show-and-tell-community-thread/7554)!

## Contributing

See [CONTRIBUTING.md](https://github.com/plotly/dash-table/blob/master/CONTRIBUTING.md)
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

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at chris@plotly.com. The project team will review and investigate all complaints, and will respond in a way that it deems appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
# Change Log for dash-html-components
### NOTE: as of v2.0, changes in dash-html-components are all being recorded in the main dash changelog.
### This file is kept only for historical purposes.
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [1.1.4] - 2021-07-09

### Changed
- [#194](https://github.com/plotly/dash-html-components/pull/194) Updated dependencies and build process
- [#190](https://github.com/plotly/dash-core-components/pull/190) Updated R package vignettes and `dash-info.yaml` to regenerate examples without attaching now-deprecated core component packages (`dashHtmlComponents`, `dashCoreComponents`, or `dashTable`).

## [1.1.3] - 2021-04-08
### Fixed
- [#179](https://github.com/plotly/dash-html-components/pull/179) - Fixes [#77](https://github.com/plotly/dash-html-components/issues/77) Added `allow` and `referrerPolicy` properties to `html.Iframe`

- [#178](https://github.com/plotly/dash-html-components/pull/178) - Fix [#161](https://github.com/plotly/dash-html-components/issues/161) <object> `data` property, and fix [#129](https://github.com/plotly/dash-html-components/issues/129) obsolete, deprecated, and discouraged elements. No elements were removed, but comments were added to the documentation about these elements detailing their limitations.

## [1.1.2] - 2021-01-19
### Fixed
- [#169](https://github.com/plotly/dash-html-components/pull/169) - part of fixing dash import bug https://github.com/plotly/dash/issues/1143

## [1.1.1] - 2020-09-03
- Dash.jl Julia component generation

## [1.1.0] - 2020-08-25
### Added
- [#165](https://github.com/plotly/dash-html-components/pull/165) Add support for Dash.jl Julia component generation.

## [1.0.3] - 2020-04-01
### Updated
- Update generated props
- Update generated R artifacts

## [1.0.2] - 2019-11-14
### Fixed
- [#143](https://github.com/plotly/dash-html-components/pull/143) Fix IE11 compatibility issues and ES5 compatibility and validation

## [1.0.1] - 2019-08-27
### Updated
- Generated documentation

## [1.0.0] - 2019-06-20
### Added
- [#119](https://github.com/plotly/dash-html-components/pull/119)
    - Added `formEncType`, `formMethod`, `formTarget` attributes to Button
    - Added `autoComplete` attribute to Select

### Removed
- [#113](https://github.com/plotly/dash-html-components/pull/113) Removed `version.py` - use `__version__` in the main namespace.

## [0.16.0] - 2019-04-25
### Fixed
- [#110](https://github.com/plotly/dash-html-components/pull/110), [#111](https://github.com/plotly/dash-html-components/pull/111) Improved the property definitions in advance of the Dev Tools property validation.
    In particular:
    - Boolean properties like `hidden` accept a bool or a case insensitive string with the same name (e.g. `'hidden'` or `'HIDDEN'`)
    - Numeric properties like `rows`, `max`, `min` allow a stringified number or a number

### Added
- Added `formNoValidate` & `inputMode` properties.

## [0.15.0] - 2019-03-25
### Changed
- Remove undefined `setProps` handling [#103](https://github.com/plotly/dash-html-components/pull/103)

## [0.14.0] - 2019-03-04
### Added
- Added `data-dash-is-loading` attribute to all components, that holds the new `loading_state.is_loading` prop.

## [0.13.5] - 2019-01-11
### Changed
- Added `.idea`, `tests`, `dist`, `.circleci` to npmignore.
- Added repository url and long_description to setup.py
- Merged in `dashHtmlComponents` R package and updated to 0.13.5

### Removed
- Removed click events - these have been obsolete since 0.7.0 [#89](https://github.com/plotly/dash-html-components/pull/89)

## [0.13.4] - 2018-12-17
### Fixed
- Fix build from wrong dash version.

## [0.13.3] - 2018-12-17
### Fixed
- `n_clicks`/`n_clicks_timestamp` PropType changed from invalid `integer` to `number`.
- omit `n_clicks`/`n_clicks_timestamp` from wrapped element props.

## [0.13.2] - 2018-09-21
### Fixed
- Fixes Python3.7 incompatibility with `0.13.0` and `0.13.1`.

### Changed
- Regenerated classes with Python3.7 to remove `async` keyword.

## [0.13.1] - 2018-09-20
### Fixed
- Renamed `__init__.py` external_path to dash_html_components.min.js

## [0.13.0] - 2018-09-20
### Added
- Unminified dev bundle support. [#64](https://github.com/plotly/dash-html-components/pull/64)

## Unreleased

## [0.12.0] - 2018-06-01
### Changed
- `dash_html_components/__init__.py` now imports from Python class files rather than generating classes at runtime,
adding support for IDE auto complete etc.

## [0.11.0] - 2018-06-01
### Added
- A `n_clicks_timestamp` property was added to all of the components. This property represents the date that the element was clicked on and can be used to determine _which element was clicked on_ in callbacks with multiple elements. This is considered a stop-gap solution: ultimately we'll want a solution that works for _all_ properties across all components, not just the `n_clicks` property. https://github.com/plotly/dash-html-components/pull/45

## [0.10.1] - 2018-04-29
### Added
- `aria-*` and `data-*` attributes are now supported in all dash html components [#40](https://github.com/plotly/dash-html-components/pull/40)
    These new keywords can be added using a dictionary expansion, e.g.
    ```
    html.Div(id="my-div", **{"data-toggle": "toggled", "aria-toggled": "true"})
    ```
- The `role` attribute was added to all components
- The `autoComplete` property was added to `textarea`

## [0.10.0] - 2018-04-03
### Added
- Previously, if a user named their app file `dash.py`, an unhelpful error
message would be raised. Now, `import dash_html_components` will check if
the user has a file named `dash.py` and warn the users appropriately.
https://github.com/plotly/dash-html-components/pull/39

## [0.9.0] - 2018-02-11
### Changed
- Moved `PropTypes` import from using `react` to using `prop-types` package to support using React 16+ in `dash-renderer`

### Added
- Added `Picture` and `Base` components
- Added `muted` property to `Audio` component

## [0.8.0] - 2017-09-29
### Added
- A `key` property has been added to every component. See https://reactjs.org/docs/lists-and-keys.html for more about this attribute.

## [0.7.0] - 2017-07-18
### Added
- A `n_clicks` property has been added to every component that gets incremented automatically when the element has been clicked on

## [0.2.3] - 2016-07-20
### Fixed
- `style` propType is now correctly set to object, not string

## [0.2.2] - 2016-07-17
### Fixed
- Issue with component metadata path in pypi package

## [0.2.0] - 2016-07-07
### Added
- Fix issues with attribute casing

## 0.1.0 - 2016-06-28
- Initial release

[0.2.3]: https://github.com/plotly/dash-html-components/compare/v0.2.2...v0.2.3
[0.2.2]: https://github.com/plotly/dash-html-components/compare/v0.2.0...v0.2.2
[0.2.0]: https://github.com/plotly/dash-html-components/compare/v0.1.0...v0.2.0
# dash-html-components

Vanilla HTML components for [Dash][]

### Install dependencies

1. Create a virtual env and activate.
    ```
    $ virtualenv venv
    $ venv/bin/activate
    ```
    _Note: venv\Scripts\activate for Windows_

2. Install Python packages required to build components.
    ```
    $ pip install -r dev-requirements.txt
    ```
3. Generate components and install npm packages
    ```
    $ npm install
    ```

### Generating HTML Components

The components in `src/components`, as well as the export index in
`src/index.js` are programmatically generated from element definitions in
`scripts/`. To regenerate:


```sh
$ npm run generate-components
```
The list of attributes is regenerated by scraping the
[MDN HTML attribute reference][].

_Note: This step will have already been done for you when you ran `npm install`_

### Development

#### Testing your components in Dash

1. Watch for changes

        $ npm run build:watch

2. Install module locally (after every change)

        # Generate metadata, and build the JavaScript bundle
        $ npm run install-local

        # Now you're done. For subsequent changes, if you've got `npm run build:watch`
        $ python setup.py install

3. Run the Dash layout you want to test

        # Import dash_html_components to your layout, then run it:
        $ python my_dash_layout.py

#### Installing Python package locally

Before publishing to PyPi, you can test installing the module locally:

```sh
# Install in `site-packages` on your machine
$ npm run install-local
```

#### Uninstalling Python package locally

```sh
$ npm run uninstall-local
```

## Contributing

See the [contributing guide](CONTRIBUTING.md) for guidelines on contributing to this project.


### Create a production build and publish:

1. Build your code:
    ```
    $ npm run build
    ```
2. Create a Python tarball
    ```
    $ python setup.py sdist
    ```
    This distribution tarball will get generated in the `dist/` folder

3. Test your tarball by copying it into a new environment and installing it locally:
    ```
    $ pip install dash-html-components-<new-version>.tar.gz
    ```

4. If it works, then you can publish the component to NPM and PyPI:
    1. Publish on PyPI
        ```
        $ twine upload dist/*
        ```
    2. Cleanup the dist folder (optional)
        ```
        $ rm -rf dist
        ```
    3. Publish on NPM (Optional if chosen False in `publish_on_npm`)
        ```
        $ npm publish
        ```
        _Publishing your component to NPM will make the JavaScript bundles available on the unpkg CDN. By default, Dash servers the component library's CSS and JS from the remote unpkg CDN, so if you haven't published the component package to NPM you'll need to set the `serve_locally` flags to `True` (unless you choose `False` on `publish_on_npm`). We will eventually make `serve_locally=True` the default, [follow our progress in this issue](https://github.com/plotly/dash/issues/284)._

[Dash]: https://plotly.com/dash
[MDN HTML attribute reference]: https://developer.mozilla.org/en-US/docs/Web/HTML/Attributes
[NPM package authors]: https://www.npmjs.com/package/dash-html-components/access
[PyPi]: https://pypi.python.org/pypi
# Contributing to dash-html-components

## Getting Started

Refer to the [readme](README.md) for installation and development workflow instructions.

## Contributions

[Dash HTML Components][] consist of generic HTML5 elements based on the [MDN spec][] (also see the [W3 index of elements][]). For more complex UI components, see the [Dash Core Components][]. Contributions are welcome! This repository's open [issues][] are a good place to start.

## Making a Contribution
_For larger features, your contribution will have a higher likelihood of getting merged if you create an issue to discuss the changes that you'd like to make before you create a pull request._

1. Create a pull request.
2. After a review has been done and your changes have been approved, they will be merged and included in a future release of Dash.
3. If significant enough, you have created an issue about documenting the new feature or change and you have added it to the [dash-docs](https://github.com/plotly/dash-docs) project.

## Coding Style

Please lint any additions to react components with `npm run lint`. Rules defined in [.eslintrc](.eslintrc) are inherited from [`dash-components-archetype`](https://github.com/plotly/dash-components-archetype)'s [eslintrc-react.json][]

## Pull Request Guidelines

Use the [GitHub flow][] when proposing contributions to this repository (i.e. create a feature branch and submit a PR against the master branch).

## Financial Contributions

If your company wishes to sponsor development of open source dash components, please [get in touch][].

[Dash HTML Components]: https://dash.plotly.com/dash-html-components
[MDN spec]: https://developer.mozilla.org/en-US/docs/Web/HTML/Element
[W3 index of elements]: https://dev.w3.org/html5/html-author/#index-of-elements
[Dash Core Components]: https://github.com/plotly/dash-core-components
[issues]: https://github.com/plotly/dash-html-components/issues
[GitHub flow]: https://guides.github.com/introduction/flow/
[eslintrc-react.json]: https://github.com/plotly/dash-components-archetype/blob/master/config/eslint/eslintrc-react.json
[get in touch]: https://plotly.com/products/consulting-and-oem

# Change Log for dash-core-components
### NOTE: as of v2.0, changes in dash-core-component are all being recorded in the main dash changelog.
### This file is kept only for historical purposes.
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [1.17.1] - 2021-07-12

### Fixed

- Removed unnecessary Julia files from npm package


## [1.17.0] - 2021-07-09

### Fixed

- [#963](https://github.com/plotly/dash-core-components/pull/963) Fixes [#885](https://github.com/plotly/dash-core-components/issues/885)

  This applies the fix from [#878](https://github.com/plotly/dash-core-components/pull/878) to the RangeSlider.
  It not only fixes the bug where the tooltips were visible when slider was not, but it also reduces the lag in the
  tooltip when the slider handles are moved.

### Updated
- [#939](https://github.com/plotly/dash-core-components/pull/939) Upgrade Plotly.js to v2.2.1. Note that this is a major version upgrade to Plotly.js, however we are not treating this as a breaking change for DCC as the majority of breaking changes in Plotly.js do not affect the Dash API. The one exception is that several trace types that have long been deprecated are removed entirely.
  - [Major release 2.0.0](https://github.com/plotly/plotly.js/releases/tag/v2.0.0):
    - Stop exporting d3 as `Plotly.d3`, and remove many other deep pieces of the public API. This does not affect the `dcc.Graph` component, but if you make use of `Plotly` from the global scope in some other way you may be affected.
    - Drop the deprecated trace types `contourgl` and `area`, as well as legacy pre-`scatterpolar` polar attributes `bar.r`, `bar.t`, `scatter.r`, `scatter.t`, `layout.radialaxis`, `layout.angularaxis`. Use `scatterpolar`, `barpolar`, and `polar` subplots instead.
    - `heatmapgl` and `pointcloud` trace types, and the `transform` attribute are deprecated, and will be removed in a future release.
    - Increase CSP safety by removing function constructors. 3D plots still use function constructors, but if you place one of the non-3D bundles (including the new `strict` bundle) in your `assets` folder you will have no function constructors.
    - Remove "Aa" text in legends.
    - Default `hovermode` to "closest".
    - Default `textposition` to "auto" in `bar` traces. If you previously used the `bar.text` attribute for hover only, you will need to explicitly set `textposition="none"`.
    - Add `bar.marker.pattern`, `image.zsmooth`, and various other features and bugfixes.
  - [Feature release 2.1.0](https://github.com/plotly/plotly.js/releases/tag/v2.1.0):
    - New `icicle` trace type.
    - New `legendrank` trace attribute.
    - Several other additions and bug fixes.
  - [Feature release 2.2.0](https://github.com/plotly/plotly.js/releases/tag/v2.2.0):
    - Legend group titles
    - Half-year directive (`%h`) for date formatting
    - Several other bug fixes and performance improvements
  - [Patch release 2.2.1](https://github.com/plotly/plotly.js/releases/tag/v2.2.1) containing a security fix.

### Added
- [#932](https://github.com/plotly/dash-core-components/pull/932) Adds a new copy to clipboard component.
- [#948](https://github.com/plotly/dash-core-components/pull/948)] Adds `disabled_days` prop to `DatePickerRange` and `DatePickerSingle` components. With this prop you can specify days that should be made unselectable in the date picker, in addition to those that fall outside of the range specified by `min_date_allowed` and `max_date_allowed`.

### Changed
- [#972](https://github.com/plotly/dash-core-components/pull/972) Updated R package vignettes and `dash-info.yaml` to regenerate examples without attaching now-deprecated core component packages (`dashHtmlComponents`, `dashCoreComponents`, or `dashTable`).

## [1.16.0] - 2021-04-08
### Added
- [#863](https://github.com/plotly/dash-core-components/pull/863) Adds a new `Download` component. Along with this several utility functions are added to help construct the appropriate data format:
  - `dcc.send_file` - send a file from disk
  - `dcc.send_data_frame` - send a `DataFrame`, using one of its writer methods
  - `dcc.send_bytes` - send a bytestring or the result of a bytestring writer
  - `dcc.send_string` - send a string or the result of a string writer

### Changed
- [#923](https://github.com/plotly/dash-core-components/pull/923)
  Set `autoComplete` to off in `dcc.Dropdown`. This fixes [#808](https://github.com/plotly/dash-core-components/issues/808)

### Fixed
- [#930](https://github.com/plotly/dash-core-components/pull/930) Fixed a bug [#867](https://github.com/plotly/dash-core-components/issues/867) with `DatePickerRange` that would sometimes shift the allowed dates by one day.
- [#934](https://github.com/plotly/dash-core-components/pull/934) Fixed a bug in `EnhancedTab` component that ignored `disabled_className` property

## [1.15.0] - 2021-01-19
### Fixed
- [#905](https://github.com/plotly/dash-core-components/pull/905) Make sure the `figure` prop of `dcc.Graph` receives updates from user interactions in the graph, by using the same `layout` object as provided in the prop rather than cloning it. Fixes [#879](https://github.com/plotly/dash-core-components/issues/879).
- [#903](https://github.com/plotly/dash-core-components/pull/903) Part of fixing dash import bug https://github.com/plotly/dash/issues/1143

### Updated
- [#911](https://github.com/plotly/dash-core-components/pull/911), [#906](https://github.com/plotly/dash-core-components/pull/906)
  - Upgraded Plotly.js to [1.58.4](https://github.com/plotly/plotly.js/releases/tag/v1.58.4)
    - Patch Release [1.58.4](https://github.com/plotly/plotly.js/releases/tag/v1.58.4)
    - Patch Release [1.58.3](https://github.com/plotly/plotly.js/releases/tag/v1.58.3)

### Added
- [#888](https://github.com/plotly/dash-core-components/pull/888) Adds a `drag_value` prop to `dcc.Slider`to be able to fire callbacks from dragging and releasing the slider.

## [1.14.1] - 2020-12-09
### Updated
- [#898](https://github.com/plotly/dash-core-components/pull/898)
    - Patch Release [1.58.2](https://github.com/plotly/plotly.js/releases/tag/v1.58.2)

## [1.14.0] - 2020-12-07
### Updated
- [#889](https://github.com/plotly/dash-core-components/pull/889), [#893](https://github.com/plotly/dash-core-components/pull/893)
  - Upgraded Plotly.js to [1.58.1](https://github.com/plotly/plotly.js/releases/tag/v1.58.1)
    - Patch Release [1.58.1](https://github.com/plotly/plotly.js/releases/tag/v1.58.1)
    - [Feature release of Plotly.js 1.58.0](https://github.com/plotly/plotly.js/releases/tag/v1.58.0) which:
      - Add `ticklabelposition` attribute to cartesian axes and colorbars [#5275](https://github.com/plotly/plotly.js/pull/5275)
      - Add "strict" `autotypenumbers` to axes and `layout` [#5240](https://github.com/plotly/plotly.js/pull/5240)
      - Add `itemwidth` to legends [#5212](https://github.com/plotly/plotly.js/pull/5212)
      - Add `root.color` attribute to `sunburst` and `treemap` traces [#5232](https://github.com/plotly/plotly.js/pull/5232), [#5245](https://github.com/plotly/plotly.js/pull/5245)
      - Enable fast image rendering for all linear axes [#5307](https://github.com/plotly/plotly.js/pull/5307)
      - Rework matches and scaleanchor so they work together [#5287](https://github.com/plotly/plotly.js/pull/5287)

## [1.13.0] - 2020-10-29
### Added
- [#871](https://github.com/plotly/dash-core-components/pull/871) Add Julia syntax highlighting support for dcc.Markdown

### Fixed
- [#878](https://github.com/plotly/dash-core-components/pull/878)
  - Fixed [#751](https://github.com/plotly/dash-core-components/issues/751), a bug that causes `dcc.Slider` and `dcc.RangerSlider` tooltips to be visible even if the slider component isn't visible (e.g. overflow),

### Updated
- [#875](https://github.com/plotly/dash-core-components/pull/875)
  - Upgraded Plotly.js to [1.57.1](https://github.com/plotly/plotly.js/releases/tag/v1.57.1)
    - Patch release [1.57.1](https://github.com/plotly/plotly.js/releases/tag/v1.57.1)
    - [Feature release of Plotly.js 1.57.0](https://github.com/plotly/plotly.js/releases/tag/v1.57.0) which:
      - Add "domain" axis references in layout `images`, `shapes` and `annotations` [#5014](https://github.com/plotly/plotly.js/pull/5014)
      - Add `rotation` attribute to `sunburst` traces [#5171](https://github.com/plotly/plotly.js/pull/5171), [#5201](https://github.com/plotly/plotly.js/pull/5201)
      - Add computed margins in "full-json" export [#5203](https://github.com/plotly/plotly.js/pull/5203)
    - [Feature release of Plotly.js 1.56.0](https://github.com/plotly/plotly.js/releases/tag/v1.56.0) which:
        - Introduce period positioning attributes on date axes in various cartesian traces [#5074](https://github.com/plotly/plotly.js/pull/5074), [#5175](https://github.com/plotly/plotly.js/pull/5175)
        - Add minexponent attribute to improve control over SI prefixes in axis tick labels [#5121](https://github.com/plotly/plotly.js/pull/5121),
        - Add sort attribute to sunburst and treemap traces to disable automatic sort [#5164](https://github.com/plotly/plotly.js/pull/5164)
        - Handle rgba colors in colorscale of surface traces [#5166](https://github.com/plotly/plotly.js/pull/5166)
    - Patch release [1.55.2](https://github.com/plotly/plotly.js/releases/tag/v1.55.2)

## [1.12.1] - 2020-09-16
### Fixed
- [#854](https://github.com/plotly/dash-core-components/pull/854) Used `persistenceTransforms` to strip the time part of the datetime in the persited props of DatePickerSingle (date) and DatePickerRange (end_date, start_date), fixing [dcc#700](https://github.com/plotly/dash-core-components/issues/700).

### Added
- [#850](https://github.com/plotly/dash-core-components/pull/850) Add property `prependData` to `Graph` to support `Plotly.prependTraces`
  + refactored the existing `extendTraces` API to be a single `mergeTraces` API that can handle both `prepend` as well as `extend`.

### Updated
- [#864](https://github.com/plotly/dash-core-components/pull/864) Upgraded Plotly.js to [1.55.2](https://github.com/plotly/plotly.js/releases/tag/v1.55.2)

## [1.12.0] - 2020-09-03
### Updated
- [#858](https://github.com/plotly/dash-core-components/pull/858)
  - Upgraded Plotly.js to [1.55.1](https://github.com/plotly/plotly.js/releases/tag/v1.55.1)
  - Patch release [1.55.1](https://github.com/plotly/plotly.js/releases/tag/v1.55.1)
  - [Feature release of Plotly.js 1.55.0](https://github.com/plotly/plotly.js/releases/tag/v1.55.0) which:
    - Introduce "period" `ticklabelmode` on cartesian date axes [#4993](https://github.com/plotly/plotly.js/pull/4993), [#5055](https://github.com/plotly/plotly.js/pull/5055), [#5060](https://github.com/plotly/plotly.js/pull/5060), [#5065](https://github.com/plotly/plotly.js/pull/5065), [#5088](https://github.com/plotly/plotly.js/pull/5088), [#5089](https://github.com/plotly/plotly.js/pull/5089)
    - Add new formatting options for weeks and quarters [#5026](https://github.com/plotly/plotly.js/pull/5026)
    - Add `source` attribute to `image` traces for fast rendering [#5075](https://github.com/plotly/plotly.js/pull/5075)
    - Add `zsmooth` attribute for discrete `heatmapgl` traces [#4953](https://github.com/plotly/plotly.js/pull/4953)
    - Add horizontal and vertical markers for arrow charts [#5010](https://github.com/plotly/plotly.js/pull/5010)
    - Add touch support to `rangeslider` [#5025](https://github.com/plotly/plotly.js/pull/5025)

## [1.11.0] - 2020-08-25
### Added
- [#851](https://github.com/plotly/dash-core-components/pull/851) Add support for Dash.jl Julia built components
- [#840](https://github.com/plotly/dash-core-components/pull/840) Add styling properties to `dcc.Loading` component
  + `parent_className`: Add CSS class for the outermost `dcc.Loading` parent div DOM node
  + `parent_style`: Add CSS style property for the outermost `dcc.Loading` parent div DOM node
  + provides a workaround for the previous behaviour the of `className` property, which changed in [#740](https://github.com/plotly/dash-core-components/pull/740). `parent_className` (or inline styles in `parent_style`) now allow CSS rules to be applied to the outermost `dcc.Loading` div, which is no longer covered by `className` on loading completion as of Dash Core Components `>= 1.9.1` (Dash `>= 1.11.0`).

## [1.10.2] - 2020-07-27
- [#835](https://github.com/plotly/dash-core-components/pull/835)
  - Upgraded Plotly.js to [1.54.7](https://github.com/plotly/plotly.js/releases/tag/v1.54.7)
    - Patch release [1.54.7](https://github.com/plotly/plotly.js/releases/tag/v1.54.7)
    - Patch release [1.54.6](https://github.com/plotly/plotly.js/releases/tag/v1.54.6)
    - Patch release [1.54.5](https://github.com/plotly/plotly.js/releases/tag/v1.54.5)
    - Patch release [1.54.4](https://github.com/plotly/plotly.js/releases/tag/v1.54.4)

## [1.10.1] - 2020-06-17
### Updated
- [#824](https://github.com/plotly/dash-core-components/pull/824)
  - Upgraded plotly.js to [1.54.3](https://github.com/plotly/plotly.js/releases/tag/v1.54.3)
    - Patch release [1.54.3](https://github.com/plotly/plotly.js/releases/tag/v1.54.3)
    - Patch release [1.54.2](https://github.com/plotly/plotly.js/releases/tag/v1.54.2)

## [1.10.0] - 2020-05-05
### Changed
- [#793](https://github.com/plotly/dash-core-components/pull/793) Added title key (i.e. HTML `title` attribute) to option dicts in `dcc.Dropdown` `options[]` list property.

### Fixed
- [#792](https://github.com/plotly/dash-core-components/pull/792) Improved the robustness of `dcc.Store` components, fixing [#456](https://github.com/plotly/dash-core-components/issues/456) whereby persistent storage could become corrupted, and fixing lifecycle issues that prevented adding `Store` components to the page after initial loading.
- [#790](https://github.com/plotly/dash-core-components/pull/790) Fixed bug where the dcc.Dropdown dropdown was hidden by the dash_table.DataTable fixed rows and columns.

### Updated
- [#800](https://github.com/plotly/dash-core-components/pull/800)
  - Upgraded plotly.js to [1.54.1](https://github.com/plotly/plotly.js/releases/tag/v1.54.1)
  - [Feature release of Plotly.js 1.54.0](https://github.com/plotly/plotly.js/releases/tag/v1.54.0) which:
    - Introduces new drag modes "drawline", "drawrect", "drawcircle", "drawopenpath", "drawclosedpath", adds optional modebar buttons for drawing & removing new shapes inside cartesian subplots, adds newshape and activeshape attributes to layout, and adds editable and fillrule attributes to layout.shapes[#4775](https://github.com/plotly/plotly.js/pull/4775)
  - Add angle and allowoverlap attributes to marker of scattermapbox traces[#4575](https://github.com/plotly/plotly.js/pull/4575), [#4794](https://github.com/plotly/plotly.js/pull/4794)
  - Also contains various other fixes

## [1.9.1] - 2020-04-10
### Changed
- [#740](https://github.com/plotly/dash-core-components/pull/740) Keep components that are loading in the DOM, but not visible, as opposed to removing them entirely. This will ensure that the size of the component's container does not shrink or expand when the component goes into the loading state.

### Fixed
- [#740](https://github.com/plotly/dash-core-components/pull/740) Fixed bug in which mapbox `uirevision` was not behaving when inside a `dcc.Loading` component

## [1.9.0] - 2020-04-01
### Changed
- [#766](https://github.com/plotly/dash-core-components/pull/766) Update from React 16.8.6 to 16.13.0
- [#768](https://github.com/plotly/dash-core-components/pull/768) Added title property to dcc.Link
- [#776](https://github.com/plotly/dash-core-components/pull/776) Update dcc.Link to set href as children if children not defined. Makes href a required prop as well.
- [#767](https://github.com/plotly/dash-core-components/pull/767) Updated dcc.Link to respond to click modifiers, and added a target prop.
- [#774](https://github.com/plotly/dash-core-components/pull/774) Fixed dcc.Location firing callbacks for wrong property.
- [772](https://github.com/plotly/dash-core-components/pull/772) Modified dcc.Link to work with absolute paths if refresh=True.

### Updated
- [#784](https://github.com/plotly/dash-core-components/pull/784)
  - [Feature release of Plotly.js 1.53.0](https://github.com/plotly/plotly.js/releases/tag/v1.53.0) which contains:
    - `rangebreaks` on date axes [#4614](https://github.com/plotly/plotly.js/pull/4614)
    - (x|y) unified `hovermode` [#4620](https://github.com/plotly/plotly.js/pull/4620)
    - "hovered data" mode to `spikesnap` [#4665](https://github.com/plotly/plotly.js/pull/4665)
    - "full-json" export format to `Plotly.toImage` and `Plotly.dowloadImage` [#4593](https://github.com/plotly/plotly.js/pull/4593)
    - node.customdata and link.customdata in `sankey` traces [#4621](https://github.com/plotly/plotly.js/pull/4621)
    - `opacityscale` for `surface` traces [#4480](https://github.com/plotly/plotly.js/pull/4480)

## [1.8.1] -2020-02-27
### Added
- [#760](https://github.com/plotly/dash-core-components/pull/760) Added R examples to package help

### Changed
- [#762](https://github.com/plotly/dash-core-components/pull/762) Renamed async modules with hyphen `-` instead of tilde `~`

## [1.8.0] - 2020-02-04
### Changed
- [#743](https://github.com/plotly/dash-core-components/pull/743) Location component now emits an event on URL path update from Link component
- [#739](https://github.com/plotly/dash-core-components/pull/739) Async Slider and RangeSlider
- [#729](https://github.com/plotly/dash-core-components/pull/729) Handle case where dcc fails to load when used inside an iframe with a sandbox attribute that only has allow-scripts

### Fixed
- [#730](https://github.com/plotly/dash-core-components/pull/730) Fixed bug in which input components with type `number` did not correctly update their values.
- [#731](https://github.com/plotly/dash-core-components/pull/731) Fixed bug where non-clearable dropdowns could still be cleared by typing backspace

### Updated
- [#747](https://github.com/plotly/dash-core-components/pull/747)
  - Upgrade plotly.js to [1.52.2](https://github.com/plotly/plotly.js/releases/tag/v1.52.2)

## [1.7.1] - 2020-01-15 (JS-only)
### Fixed
- [#734](https://github.com/plotly/dash-core-components/pull/734) Fix JS-facing release bug where `Plotly.js` was listed in `devDependencies` instead of `dependencies`

## [1.7.0] - 2020-01-14
### Added
- [#711](https://github.com/plotly/dash-core-components/pull/711) Added support for `dcc.Link` (dccLink) and nested `dcc.Markdown` (dccMarkdown) react components inside of `dcc.Markdown`
- [#706](https://github.com/plotly/dash-core-components/pull/706)
  - Added new `responsive` property that overrides the underlying Plotly.js graph responsiveness from Dash-land
  - Added responsiveness on graph parent element resize (previously only worked on window.resize)
  - Added new `dash-graph--pending` class to dcc.Graph, present while resizing, (re-)rendering, loading

### Changed
- [#723](https://github.com/plotly/dash-core-components/pull/723) Changed npm package content to allow source code inclusion from other projects
- [#725](https://github.com/plotly/dash-core-components/pull/725) Improve async graph performance by parallelizing resource fetching instead of fetching sequentially
- [#720](https://github.com/plotly/dash-core-components/pull/720) `highlight.js` is now bundled into the package, and no longer sets the `window.hljs` variable. Similarly to how `plotly.js` is handled, it is overridden by a user-provided version if one exists.

### Updated
- [#732](https://github.com/plotly/dash-core-components/pull/732)
  - Upgraded plotly.js to [1.52.1](https://github.com/plotly/plotly.js/releases/tag/v1.52.1)
  - [Feature release 1.52.0](https://github.com/plotly/plotly.js/releases/tag/v1.52.0) which contains:
    - Enable loading locale bundles before plotly.js bundles [#4453](https://github.com/plotly/plotly.js/pull/4453)
    - `ko` localization [#4315](https://github.com/plotly/plotly.js/pull/4315)
  - Patch release [1.52.1](https://github.com/plotly/plotly.js/releases/tag/v1.52.1) containing several bug fixes.
- [#706](https://github.com/plotly/dash-core-components/pull/706)
  - Upgraded plotly.js to [1.51.3](https://github.com/plotly/plotly.js/releases/tag/v1.51.3)

## [1.6.0] - 2019-11-27
### Updated
- Upgraded plotly.js to 1.51.2 [#708](https://github.com/plotly/dash-core-components/pull/708)
  - Patch release [1.51.2](https://github.com/plotly/plotly.js/releases/tag/v1.51.2) containing several bug fixes.

### Changed
- [#695](https://github.com/plotly/dash-core-components/pull/695) Improvements to Slider and RangeSlider
  - Marks outside of the range specified by `min` and `max` are now omitted when the slider renders.
  - Padding is now dependent on the orientation (vertical or horizontal), and whether or not tooltips are always displayed.
  - The whitespace is now preserved for `marks` labels.

### Added
- [#695](https://github.com/plotly/dash-core-components/pull/695) Added new property `verticalHeight` to Slider and RangeSlider, to allow the user to specify the height (in px) of vertical sliders. This defaults to `400`.

## [1.5.1] - 2019-11-14
### Fixed
- [#696](https://github.com/plotly/dash-core-components/pull/696) Fix IE11 compatibility issues and ES5 compatibility and validation

### Changed
- [#687](https://github.com/plotly/dash-core-components/pull/687/) Use `start_date`, `min_date_allowed`, `end_date`, or `max_date_allowed` for the initial visible month if the value of the parameter `initial_visible_month` is not supplied.

## [1.5.0] - 2019-11-04
### Added
- [#692](https://github.com/plotly/dash-core-components/pull/692) Async DatePickerSingle, DatePickerRange, Dropdown, Markdown, Upload components

### Updated
- [#693](https://github.com/plotly/dash-core-components/pull/693) Upgraded plotly.js to 1.51.1
  - [Feature release 1.51.0](https://github.com/plotly/plotly.js/releases/tag/v1.51.0) which contains:
    - A new `image` trace type to display 3- or 4-channel color images as data
    - `automargin` for `pie` charts for better readability when labeling lots of small slices
    - Toggle-type `updatemenus`
    - `zh-CN` localization
    - And various other small features and bug fixes
  - Patch release [1.51.1](https://github.com/plotly/plotly.js/releases/tag/v1.51.1) containing several bug fixes.

## [1.4.0] - 2019-10-29
### Added
- [#616](https://github.com/plotly/dash-core-components/pull/616) Async Graph and Plotly.js

## [1.3.1] - 2019-10-17
### Updated
- Upgraded plotly.js to 1.50.1 [#681](https://github.com/plotly/dash-core-components/issues/681)
  - Patch release [1.50.1](https://github.com/plotly/plotly.js/releases/tag/v1.50.1) containing several bug fixes.

### Fixed
- [#681](https://github.com/plotly/dash-core-components/issues/681) Fix a bug with the dcc.Graph component logging errors in certain circumstances when nested inside a dcc.Loading component

## [1.3.0] - 2019-10-08
### Added
- Added `search_value` prop to `Dropdown`, for server-side options loading/filtering. [#660](https://github.com/plotly/dash-core-components/pull/660)

### Updated
- Upgraded plotly.js to 1.50.0 [#675](https://github.com/plotly/dash-core-components/pull/675)
  - [Feature release 1.50.0](https://github.com/plotly/plotly.js/releases/tag/v1.50.0) which contains:
    - A new `treemap` trace type for display of hierarchical data.
    - `texttemplate` support for all traces with on-graph text, and custom date formatting for templated on-graph and hover text.
    - Transitions (animation) for `bar` charts.
    - Numerous other performance improvements, features, and bug fixes.
  - Patch release [1.49.5](https://github.com/plotly/plotly.js/releases/tag/v1.49.5) containing several bug fixes.

## [1.2.1] - 2019-09-19
### Fixed
- Fix regression in DatePickerRange, DatePickerSingle, Input
[#652](https://github.com/plotly/dash-core-components/issues/652)


## [1.2.0] - 2019-09-17
### Added
- Added support for persistence of user-edited props to value-input components: `Checklist`, `DatePickerRange`, `DatePickerSingle`, `Dropdown`, `Input`, `RadioItems`, `RangeSlider`, `Slider`, `Tabs`, and `Textarea`. New props are `persistence`, `persistence_type`, and `persisted_props`. Set `persistence` to a truthy value to enable, the other two modify persistence behavior. See [plotly/dash#903](https://github.com/plotly/dash/pull/903) for more details. [#646](https://github.com/plotly/dash-core-components/pull/646)

### Fixed
- Fixed `Slider` and `RangeSlider` components with `tooltip.always_visible` [#640](https://github.com/plotly/dash-core-components/issues/640)

- Fixed an infinite loop problem when `Graph` is wrapped by `Loading` component [#608](https://github.com/plotly/dash-core-components/issues/608)

## [1.1.2] - 2019-08-27
### Fixed
- Fixed problems with `Graph` components leaking events and being recreated multiple times if declared with no ID [#604](https://github.com/plotly/dash-core-components/pull/604)

- Fixed problem with `DatePickerRange` component about `clearable` not working [#614](https://github.com/plotly/dash-core-components/issues/614) and [#594](https://github.com/plotly/dash-core-components/issues/594)

### Updated
- Upgraded plotly.js to 1.49.4 [#612](https://github.com/plotly/dash-core-components/issues/612)
  - Patch releases [1.49.4](https://github.com/plotly/plotly.js/releases/tag/v1.49.4), [1.49.3](https://github.com/plotly/plotly.js/releases/tag/v1.49.3), [1.49.2](https://github.com/plotly/plotly.js/releases/tag/v1.49.2)


## [1.1.1] - 2019-08-06
### Updated
- Upgraded plotly.js to 1.49.1 [#595](https://github.com/plotly/dash-core-components/issues/595)
  - Patch release [1.49.1](https://github.com/plotly/plotly.js/releases/tag/v1.49.1)

## [1.1.0] - 2019-08-05
### Changed
- Fixed inconsistent behavior of `input` with `type=number` [#580](https://github.com/plotly/dash-core-components/pull/580)

### Updated
- Upgraded plotly.js to 1.49.0 [#589](https://github.com/plotly/dash-core-components/pull/589)
  - [Feature release 1.49.0](https://github.com/plotly/plotly.js/releases/tag/v1.49.0) which contains:
    - New `indicator` trace type for gauge and KPI displays.
    - Lots of tile map improvements: `choroplethmapbox` and `densitymapbox` trace types, numerous `style` options for `mapbox` subplots that do not require a Mapbox access token, and more.
    - Various bug fixes and smaller improvements.

## [1.0.0] - 2019-06-20
### Added
- `Markdown` components support code highlighting - no need to switch to `SyntaxHighlighter`, which has been removed. Use triple backticks, with the opening backticks followed by the language name or abbreviation. [#562](https://github.com/plotly/dash-core-components/pull/562) Supported languages:
    - Bash
    - CSS
    - HTTP
    - JavaScript
    - Python
    - JSON
    - Markdown
    - HTML, XML
    - R
    - Ruby
    - SQL
    - Shell Session
    - YAML
- Added a `dedent` prop to `Markdown` components, and enabled it by default - removing all matching leading whitespace from every line that has any non-whitespace content. You can disable this with `dedent=False`. [#569](https://github.com/plotly/dash-core-components/pull/569)
- Ability to add tooltips to `Slider` and `RangeSlider`, which can be visible always or on hover. Tooltips also take a position argument. [#564](https://github.com/plotly/dash-core-components/pull/564)

### Fixed
- Fixed `min_date_allowed` and `max_date_allowed` bug in `DatePickerRange` [#551](https://github.com/plotly/dash-core-components/issues/551)
- Fixed unwanted `resize()` calls on unmounted `Graph`s [#534](https://github.com/plotly/dash-core-components/issues/534)
- Fixed `tab--disabled` CSS class issue in `Tab` component with custom styling [#568](https://github.com/plotly/dash-core-components/pull/568)

### Changed
- Changed `dcc.Checklist` prop `values` to `value`, to match all the other input components [#558](https://github.com/plotly/dash-core-components/pull/558). Also improved prop types for `Dropdown` and `RadioItems` `value` props to consistently accept both strings and numbers.

### Removed
- 💥 Removed the `SyntaxHighlighter` component. This is now built into `Markdown` [#562](https://github.com/plotly/dash-core-components/pull/562).
- Removed the `containerProps` prop in `Markdown` - after the refactor of [#562](https://github.com/plotly/dash-core-components/pull/562), its function is served by the `id`, `className`, and `style` props. [#569](https://github.com/plotly/dash-core-components/pull/569)
- Removed `version.py` - use `__version__` in the main namespace instead. [#555](https://github.com/plotly/dash-core-components/pull/555)

### Updated
- Upgraded plotly.js to 1.48.3 [#571](https://github.com/plotly/dash-core-components/pull/571)
  - [Feature release 1.48.0](https://github.com/plotly/plotly.js/releases/tag/v1.48.0) which contains:
    - New `funnel` and `funnelarea` trace types
    - Shared color axes and colorbars
    - Sorting cartesian axes by the value on the opposite axis
    - Improvements to `bar` & `waterfall` text, legend clicking, histogram binning, hover text, and more
  - Patch releases [1.48.3](https://github.com/plotly/plotly.js/releases/tag/v1.48.3), [1.48.2](https://github.com/plotly/plotly.js/releases/tag/v1.48.2), [1.48.1](https://github.com/plotly/plotly.js/releases/tag/v1.48.1), [1.47.4](https://github.com/plotly/plotly.js/releases/tag/v1.47.4), [1.47.3](https://github.com/plotly/plotly.js/releases/tag/v1.47.3), [1.47.2](https://github.com/plotly/plotly.js/releases/tag/v1.47.2), [1.47.1](https://github.com/plotly/plotly.js/releases/tag/v1.47.1) containing numerous bug fixes

## [0.48.0] - 2019-05-15
### Added
- `figure` prop in `dcc.Graph` now accepts a `frames` key
- Improved the `Dropdown` options description for dash-docs [#547](https://github.com/plotly/dash-core-components/pull/547)
- Added `optionHeight` prop to `Dropdown` [#552](https://github.com/plotly/dash-core-components/pull/552)
- Merged in R `dashCoreComponents` package and updated to 0.48.0

### Removed
- Removed unused `key` prop from `dcc.ConfirmDialog`

## [0.47.0] - 2019-04-25
### Fixed
- Fixed style regression in DatePickerSingle and DatePickerRange [#518](https://github.com/plotly/dash-core-components/issues/518)
- **Breaking** - In `dcc.Input`, fixed several HTML properties that weren't properly camel cased and therefore were not actually applied in the DOM [#523](https://github.com/plotly/dash-core-components/pull/523/):
    - `autocomplete` is now `autoComplete`
    - `autofocus` is now `autoFocus`
    - `inputmode` is now `inputMode`
    - `maxlength` is now `maxLength`
    - `minlength` is now `minLength`
- Improved property definitions of various components in anticipation of upcoming component validation (https://github.com/plotly/dash-renderer/pull/100/). [#523](https://github.com/plotly/dash-core-components/pull/523/)
- **Breaking** - `n_blur_timestamp` & `n_submit_timestamp` in `Input` & `Textarea` is now a number instead of a date object/string. This matches the form of `n_clicks_timestamp` as used in `dash_html_components`. [#523](https://github.com/plotly/dash-core-components/pull/523/)
- Fixed an issue with `ConfirmDialog` that could display multiple confirm popups instead of a single popup in certain contexts. [#523](https://github.com/plotly/dash-core-components/pull/523/)
- `dcc.Markdown` `containerProps` is now applied to the component's container again. This was broken in 0.45.0's `react-markdown` upgrade.

### Changed
- `dcc.Interval` will reset its timer when re-enabled. Previously, if the `dcc.Interval` was disabled, it's "clock would keep running": when it was reenabled, it would fire at the next `interval` on the previous clock's schedule, rather than `interval` milliseconds later. For example, previously the schedule might look like this:
    ```
    0 seconds: interval started with `interval=5000`
    5 seconds: `n_intervals=1`
    7 seconds: callback sets `disabled=True`
    10 seconds: interval continues to run, but doesn't fire an update
    13 seconds: callback sets `disabled=False`
    15 seconds: interval fires an update: `n_intervals=2`
    ```

    Now, it will look like this:
    ```
    0 seconds: interval started with `interval=5000`
    5 seconds: `n_intervals=1`
    7 seconds: callback sets `disabled=True` - interval stops
    13 seconds: callback sets `disabled=False` - clock resets
    18 seconds: interval fires an update: `n_intervals=2`
    ```

## [0.46.0] - 2019-04-10
### Added
- `extendData` prop for `Graph` component. This feeds `Plotly.extendTraces` for incremental data updates. [#461](https://github.com/plotly/dash-core-components/pull/461)

### Changed
[#508](https://github.com/plotly/dash-core-components/pull/508)
- Upgrade from React 15.4.2 to 16.8.6
- Upgrade from react-date 12.3.0 to 20.1.0

### Fixed
- Fix unnecessary `loading_state` prop for `Input` component. [#498](https://github.com/plotly/dash-core-components/issues/498)
- Ensure `DatePickerSingle` callbacks fire with cleared dates. [#511](https://github.com/plotly/dash-core-components/pull/511)
- Fixes incorrect default values for `config` prop of `Graph`. [#515](https://github.com/plotly/dash-core-components/pull/515)

### Updated
- Upgraded plotly.js to 1.47.0 [#516](https://github.com/plotly/dash-core-components/pull/516)
  - [Feature release 1.47.0](https://github.com/plotly/plotly.js/releases/tag/v1.47.0) which contains:
    - New `volume` gl3d trace type
    - Interactive node grouping for Sankey diagrams, using box or lasso selection
    - Add way for Plotly.toImage and Plotly.downloadImage to export images with current graph width/height by passing width/height option as null
    - Improvements to hover labels, legends, and more
  - [Feature release 1.46.0](https://github.com/plotly/plotly.js/releases/tag/v1.46.0) which contains:
    - New `waterfall` trace type
    - New `sunburst` trace type
    - Implement connectgaps on surface traces
    - Implement hovertemplate for box and violin points
    - Display hover labels above modebar, ensuring that the hover labels are always visible within the graph div
  - Patch releases [1.46.1](https://github.com/plotly/plotly.js/releases/tag/v1.46.1), [1.45.3](https://github.com/plotly/plotly.js/releases/tag/v1.45.3), [1.45.2](https://github.com/plotly/plotly.js/releases/tag/v1.45.2), and [1.45.1](https://github.com/plotly/plotly.js/releases/tag/v1.45.1) containing numerous bug fixes

## [0.45.0] - 2019-03-25
### Added
- `restyleData` prop for `Graph` component [#483](https://github.com/plotly/dash-core-components/pull/483)

### Changed
- `dcc.Markdown` now uses GitHub-flavored markdown instead of CommonMark markdown. This was done as part of upgrading `react-markdown` third party library from 2.4.5 to [4.0.6](https://github.com/rexxars/react-markdown/blob/master/CHANGELOG.md#406---2019-01-04). Compare the differences in these online editors: [CommonMark editor](https://spec.commonmark.org/dingus/), [GitHub Markdown Editor](https://rexxars.github.io/react-markdown/). Notable changes:
    - A line break is now needed between paragraphs and lists.
    - Links are automatically rendered.
    - Many more features are supported like tables and strikethrough.

### Fixed
- Fix Vertical Slider regression [#479](https://github.com/plotly/dash/issues/479)
- Fix Slider regression [#485](https://github.com/plotly/dash/issues/485)

## [0.44.0] - 2019-03-04
### Added
- Loading component [#267](https://github.com/plotly/dash/issues/267)

### Updated
- Upgraded plotly.js to 1.45.0 [#470](https://github.com/plotly/dash-core-components/pull/470)
  - [Feature release 1.45.0](https://github.com/plotly/plotly.js/releases/tag/v1.45.0) which contains:
     - Sankey diagram improvements including circular networks, node grouping, and concentration colorscales
     - Matching cartesian axes
     - Better bar, box, and violin alignment control
     - Orthographic 3D projections
     - Hovertemplate support for more trace types, including all 3D traces
     - And many other features and bug fixes
  - Patch release [1.44.4](https://github.com/plotly/plotly.js/releases/tag/v1.44.4) containing numerous bug fixes

## [0.43.1] - 2019-02-11
### Updated
- Upgraded plotly.js to 1.44.3 [#458](https://github.com/plotly/dash-core-components/pull/458)
  - Patch releases [1.44.2](https://github.com/plotly/plotly.js/releases/tag/v1.44.2) and [1.44.3](https://github.com/plotly/plotly.js/releases/tag/v1.44.3) containing numerous bug fixes

## [0.43.0] - 2019-01-25
### Added
- Added event props `n_blur` and `n_clicks` - along with `n_blur_timestamp` and `n_clicks_timestamp` - in `Textarea` components, to maintain the functionality lost by removing the `click` and `blur` events. All other events were already covered by existing props. [#444](https://github.com/plotly/dash-core-components/pull/444)
- Merged in R `dashCoreComponents` package and updated to 0.43.0

### Fixed
- Fix dynamically disabling and enabling `Interval` components [#436](https://github.com/plotly/dash-core-components/pull/436)
- Clear date in DatePickerSingle and DatePickerRange [#434](https://github.com/plotly/dash-core-components/issues/434)

### Removed
- Removed `Event` system - see https://github.com/plotly/dash/issues/531 for details. [#444](https://github.com/plotly/dash-core-components/pull/444)

### Updated
- Upgraded plotly.js to 1.44.1 [#445](https://github.com/plotly/dash-core-components/pull/445)
  - [Feature release 1.44.0](https://github.com/plotly/plotly.js/releases/tag/v1.44.0) which contains:
    - A new `isosurface` gl3d trace type
    - Animated transitions via `Plotly.react` using `layout.transitions`
    - `hovertemplate` support in many more trace types
    - And many other features and bug fixes
  - Patch releases [1.44.1](https://github.com/plotly/plotly.js/releases/tag/v1.44.1) and [1.43.2](https://github.com/plotly/plotly.js/releases/tag/v1.43.2) containing numerous bug fixes

## [0.42.1] - 2019-01-07
### Fixed
- Fix `dcc.Store` type changes [#427](https://github.com/plotly/dash-core-components/pull/427)

## [0.42.0] - 2018-12-27
### Fixed
- Fix `dcc.Store` null values in list causing an infinite loop [#424](https://github.com/plotly/dash-core-components/pull/424)

### Updated
- Upgraded plotly.js to 1.43.1 [#423](https://github.com/plotly/dash-core-components/pull/423). This includes:
  - [Feature release 1.43.0](https://github.com/plotly/plotly.js/releases/tag/v1.43.0), which contains:
    - `multicategory` axis type for hierarchical categories
    - `uirevision` attributes to control persistence of user-driven changes to the graph
    - `hovertemplate` for more control of hover labels for certain trace types
    - Alignment options for titles and legend items
    - And many other features and bug fixes
  - Patch releases [1.43.1](https://github.com/plotly/plotly.js/releases/tag/v1.43.1), [1.42.5](https://github.com/plotly/plotly.js/releases/tag/v1.42.5), [1.42.4](https://github.com/plotly/plotly.js/releases/tag/v1.42.4), and [1.42.3](https://github.com/plotly/plotly.js/releases/tag/v1.42.3) containing numerous bug fixes

## [0.41.0] - 2018-12-11
### Added
- `dangerously_allow_html` prop for Markdown component for allowing HTML.

## [0.40.5] - 2018-12-11
### Fixed
- Fix typos in DatePickerSingle props [#361](https://github.com/plotly/dash-core-components/pull/361)

## [0.40.4] - 2018-12-10
### Fixed
- Add map files to manifest [#413](https://github.com/plotly/dash-core-components/pull/413)

## [0.40.3] - 2018-12-07
### Added
- Source map [#404](https://github.com/plotly/dash-core-components/issues/404)
    Related Dash issue [#480](https://github.com/plotly/dash/issues/480)

## [0.40.2] - 2018-12-04
### Fixed
- Put Input value set in onBlur/onSubmit under a debounce check [#384](https://github.com/plotly/dash-core-components/pull/384)

## [0.40.1] - 2018-12-04
### Fixed
- Fixed issue [#390](https://github.com/plotly/dash-core-components/issues/390) by providing better styles for vertical Tabs.

## [0.40.0] - 2018-11-28
### Added
- Add Logout button (dash-deployment-server authentication integration) [#388](https://github.com/plotly/dash-core-components/pull/388)

## [0.39.0] - 2018-11-12
### Changed
- Updated `react` and `react-dom` to version `^16.6.1`
- Updated `react-docgen` to `^2.21.0`
- Updated `react-select-fast-filter-options` to `^0.2.3`
- Updated `react-virtualized-select` to `^3.1.3`
- Upgraded `babel` and dependencies to `7.1.5`
- Upgraded `enzyme` and dependencies to `3.7.0`
- Removed `react-select` because it's unused - we're using `react-virtualized-select` instead.

## [0.38.1] - 2018-11-14
### Fixed
- The issue [#115](https://github.com/plotly/dash-core-components/issues/115)
with datepicker, which didn't appear when the date value is `None` was fixed.
By default, current month would appear in such cases for
`DatePickerRange` and `DatePickerSingle` components.
See pull request https://github.com/plotly/dash-core-components/pull/201.
- Refactored the way the Graph component would generate an unique id if none provided.
- Default CSS imported via `style-loader` is now placed at top, so that user supplied CSS can overwrite it, fixes [#380](https://github.com/plotly/dash-core-components/issues/380)

## [0.38.0] - 2018-11-07
### Fixed
- Changed the way the default CSS files for some components are loaded to being loaded with webpack instead of as dependencies.

## [0.37.2] - 2018-11-07
### Changed
- Updated `react-select` to version `2.1.0`

## [0.37.1] - 2018-11-07
### Added
- Added `clickannotation` event to `dcc.Graph`.
See https://github.com/plotly/dash-core-components/pull/182.

## [0.37.0] - 2018-11-04
### Fixed
- Some Input props weren't being picked up by React. Changed:
    - `autocomplete` to `autoComplete`
    - `autofocus` to `autoFocus`
    - `inputmode` to `inputMode`
    - `maxlength` to `maxLength`
    - `minlength` to `minLength`
### Added
- Unit tests for `Input` component.
- New `debounce` prop for `Input` component that determines if the input should update to a Dash server immediately, or on 'Enter' key. Fixes [#169](https://github.com/plotly/dash-core-components/issues/169)
### Changed
- `min` and `max` prop now won't update the input when values are lower than or greater than `min` and `max` respectively. Fixes[#173](https://github.com/plotly/dash-core-components/issues/173)
- `step` prop can now be a `number` and is therefor set correctly on the corresponding `<input/>` tag. Fixes [#292](https://github.com/plotly/dash-core-components/issues/292)

## [0.36.0] - 2018-11-01
### Fixed
- The `npm start` command now runs the Demo app again [#346](https://github.com/plotly/dash-core-components/issues/346)

## [0.36.0] - 2018-10-31
### Updated
- Upgraded plotly.js to 1.42.2 [#354](https://github.com/plotly/dash-core-components/pull/354). This includes:
  - [Feature release 1.42.0](https://github.com/plotly/plotly.js/releases/tag/v1.42.0) which contains:
    - A new trace type `parcats` (parallel categories)
    - `modebar` styling
    - `pie` trace titles
    - And many other features and bug fixes
  - Patch releases [1.42.2](https://github.com/plotly/plotly.js/releases/tag/v1.42.2) and [1.42.1](https://github.com/plotly/plotly.js/releases/tag/v1.42.1) containing numerous bug fixes

## [0.35.2] - 2018-10-30
### Fixed
- Fix Input not used in callbacks resetting the value on updates. [#350](https://github.com/plotly/dash-core-components/pull/350)

## [0.35.1] - 2018-10-29
### Fixed
- Fix Dropdown options prop docstring typo [#328](https://github.com/plotly/dash-core-components/pull/328/files)
- Fix Input readonly prop -> readOnly [#348](https://github.com/plotly/dash-core-components/pull/348)

## [0.35.0] - 2018-10-22
### Added
- n_blur/n_submit and timestamps props added to the Input component. [#326](https://github.com/plotly/dash-core-components/pull/326)

## [0.34.0] - 2018-10-17
### Added
- `npm run test-unit` will run new Jest+Enzyme unit tests
- Unit tests for Tabs component
### Fixed
- Fixed bug in Tabs component where value was resetting if using callback-less mode [#331](https://github.com/plotly/dash-core-components/issues/331)
- Fixed bug with default Tabs value not being set to children's Tab value (if it's set)
- Fixed bug where Tabs.children.props wheren't being selected properly, related to [#84](https://github.com/plotly/dash-renderer/issues/84)

## [0.33.1] -- 2018-10-17
### Fixed
- Fix Store component nested data [#333](https://github.com/plotly/dash-core-components/pull/333)

## [0.33.0] -- 2018-10-04
### Fixed
- Patched plotly.js version 1.41.3 [#320](https://github.com/plotly/dash-core-components/pull/320). This includes patch releases [1.41.3](https://github.com/plotly/plotly.js/releases/tag/v1.41.3), [1.41.2](https://github.com/plotly/plotly.js/releases/tag/v1.41.2), and [1.41.1](https://github.com/plotly/plotly.js/releases/tag/v1.41.1) containing numerous bug fixes.

## [0.32.0] - 2018-10-2
### Added
- Added Store component [#248](https://github.com/plotly/dash-core-components/pull/248)


## [0.31.0] - 2018-09-21
### Changed
- Updated NPM scripts:
  - `test` now runs Selenium integration tests
  - `format` runs Prettier formatter
  - There are new `build` scripts, most notably `build:watch` runs a watcher and rebuilds upon changes
  - There's a new `publish-all` script that publishes to NPM and PyPi
### Fixed
- The `start` script will now run the `Demo` application

## [0.30.2] - 2018-09-21
### Fixed
- Fixed regression in Graph component where it wouldn't resize correctly [#256](https://github.com/plotly/dash-core-components/issues/256)

## [0.30.1] - 2018-09-20
### Fixed
- Renamed `__init__.py` external_path to dash_core_components.min.js

## [0.30.0] - 2018-09-20
### Added
- Unminified dev bundle support. [#293](https://github.com/plotly/dash-core-components/pull/293)

## [0.29.0] -- 2018-09-13
### Updated
- Upgraded Plotly.js, the underlying library behind the `dash_core_components.Graph` component, to version 1.41.0 [#302](https://github.com/plotly/dash-core-components/pull/302). Feature release [1.41.0](https://github.com/plotly/plotly.js/releases/tag/v1.41.0) adds:
  - Stacked area charts, as part of `scatter` traces
  - A new `barpolar` trace type for polar "bar" traces. This replaces and improves the deprecated `area` trace type.
  - A `responsive` plot config option
  - And various other features and bug fixes.

  Many of these features were funded directly by companies that rely on this library. If your organization or company would like to sponsor particular features or bug fixes in these open source libraries, please reach out: http://plotly.com/products/consulting-and-oem

## [0.28.3] - 2018-09-07
### Changed
- The `Interval` component's `max_interval` prop can now be used to stop/restart the interval. Fixes [#266](https://github.com/plotly/dash-core-components/issues/266)
- The `Graph` component's `id` is now not required to be set.
### Fixed
- Fixed bug where Graph would resize randomly when rerendered, for example in a dcc.Tabs component.

## [0.28.2] - 2018-09-06
### Fixed
- Fixed bug in Tabs component where initial tab content wasn't rendering, [#282](https://github.com/plotly/dash-core-components/issues/282)
- Fixed bug in Tabs component where no default Tab is selected if Tabs.value is empty

## [0.28.1] - 2018-08-29
### Changed
- `candlestick` and `OHLC` charts are now plotted using the `Plotly.react` method instead of the `Plotly.newPlot` method.
### Fixed
- Fix bug where front-end error was thrown when setting `Graph.figure = {}` (fixes [#260]).

## [0.28.0]
### Updated
- Upgraded plotly.js to version 1.40.1 [#280](https://github.com/plotly/dash-core-components/pull/280). This includes [Feature release 1.40.0](https://github.com/plotly/plotly.js/releases/tag/v1.40.0) and patch releases [1.40.1](https://github.com/plotly/plotly.js/releases/tag/v1.40.1), [1.39.4](https://github.com/plotly/plotly.js/releases/tag/v1.39.4), [1.39.3](https://github.com/plotly/plotly.js/releases/tag/v1.39.3), and [1.39.2](https://github.com/plotly/plotly.js/releases/tag/v1.39.2), containing:
  - `contour` trace legend items
  - `piecolorway` and `extendpiecolors` attributes for more control over `pie` colors
  - And many other features and bug fixes

## [0.27.2]
### Fixed
- `Tabs.children` can now be undefined, so you can update them dynamically. [#265](https://github.com/plotly/dash-core-components/issues/265)
- `Tabs` Callback-less version no longer has the 2nd tab selected by default [#262](https://github.com/plotly/dash-core-components/issues/262)
- Fixes bug with nested `Tabs` [#273](https://github.com/plotly/dash-core-components/issues/273) and [#272](https://github.com/plotly/dash-core-components/issues/272)

## [0.27.1]
### Fixed
- `ConfirmDialogProvider` can now be used without a callback. [#241](https://github.com/plotly/dash-core-components/pull/241)
- `ConfirmDialog`, only fire `submit` when `submit` is clicked. [#242](https://github.com/plotly/dash-core-components/issues/242) fixed in [#241](https://github.com/plotly/dash-core-components/pull/241)

## [0.27.0]
### Changed
- `dash_core_components/__init__.py` now imports from python class files rather than generating classes at runtime,
adding support for IDE autocomplete etc.

## [0.26.0]
### Added
- New Tabs and Tab components! [#213](https://github.com/plotly/dash-core-components/pull/213#pullrequestreview-135893345)

## [0.25.1]
### Fixed
- `__init__` version formatting for unpkg.

## [0.25.0]
### Added
- `ConfirmDialog` and `ConfirmDialogProvider` components [#211](https://github.com/plotly/dash-core-components/pull/211)

## [0.24.1]
### Fixed
- Improved DatePickerRange, fixing issues [#209](https://github.com/plotly/dash-core-components/issues/209) and [#152](https://github.com/plotly/dash-core-components/issues/152)
- Link component now is a proper <a> tag so you can right click on it, and will scroll back to top. Fixes [#99](https://github.com/plotly/dash-core-components/issues/99), implemented in [#215](https://github.com/plotly/dash-core-components/pull/215)
- Added `max_interval` prop to `Interval` component, fixing issue [#222](https://github.com/plotly/dash-core-components/issues/222)


## [0.24.0]
### Updated
- Upgraded plotly.js, the underlying library behind the
`dash_core_components.Graph` component, to version 1.39.1 [#228](https://github.com/plotly/dash-core-components/pull/228). This includes:
  - [Feature release 1.39.0](https://github.com/plotly/plotly.js/releases/tag/v1.39.0), which contains:
    - `layout.template` capability, along with `Plotly.makeTemplate` and `Plotly.validateTemplate`.
    - A new 3D `streamtube` trace type
    - `scattergl` on-graph text
    - `polar` polygon grids
    - And many other features and bug fixes
  - Patch releases [1.39.1](https://github.com/plotly/plotly.js/releases/tag/v1.39.1), [1.38.3](https://github.com/plotly/plotly.js/releases/tag/v1.38.3), [1.38.2](https://github.com/plotly/plotly.js/releases/tag/v1.38.2), and [1.38.1](https://github.com/plotly/plotly.js/releases/tag/v1.38.1) containing numerous bug fixes.
  - Many of these features were funded directly by companies that rely on this library. If your organization or company would like to sponsor particular features or bug fixes in these open source libraries, please reach out: http://plotly.com/products/consulting-and-oem

## [0.23.0]
### Updated
- Upgraded plotly.js to version 1.38.0 [#207](https://github.com/plotly/dash-core-components/pull/207). This includes:
  - Feature releases [1.38.0](https://github.com/plotly/plotly.js/releases/tag/v1.38.0), [1.37.0](https://github.com/plotly/plotly.js/releases/tag/v1.37.0), and [1.36.0](https://github.com/plotly/plotly.js/releases/tag/v1.36.0), which contain:
    - A new 3D `cone` trace type to visualize vector fields
    - A new `splom` (aka Scatter PLOt Matrix) trace type
    - Performance improvements for `splom` and other multi-subplot graphs
    - `plotly_legendclick` and `plotly_legenddoubleclick` events
    - Multi-selection and click-to-select for `parcoords` axes
    - Fixed-size shapes
    - And many other features and bug fixes
  - Patch releases [1.37.1](https://github.com/plotly/plotly.js/releases/tag/v1.37.1) and [1.36.1](https://github.com/plotly/plotly.js/releases/tag/v1.36.1) containing numerous bug fixes

## [0.22.2] - 2018-05-22
### Fixed
- `dcc.Input` component now handles `disabled=False` property.
- Broken sourcemaps for debugging.
### Added
- Testing configuration for CHROMEPATH and SERVER_PROCESSES

## [0.22.1] - 2018-04-09
### Fixed
- Various bugs with the `ohlc` and `candlestick` chart type in the `dcc.Graph`
component were fixed. See https://github.com/plotly/dash-core-components/pull/184.

## [0.22.0] - 2018-04-03
### Added
- Previously, if a user named their app file `dash.py`, an unhelpful error
message would be raised. Now, `import dash_core_components` will check if
the user has a file named `dash.py` and warn the users appropriately.
https://github.com/plotly/dash-core-components/pull/177

## [0.21.1] - 2018-03-28
### Fixed
- In some cases, frequently multi-page apps, the `dcc.Graph` interactive properties
will stop working (`selectedData`, `hoverData`, `relayoutData`). This should be fixed now. https://github.com/plotly/dash-core-components/pull/178
- `dcc.Graph` will now resize after it it plotted for the first time. This should fix issues
where the `dcc.Graph` component was not fitting to the size of its container. https://github.com/plotly/dash-core-components/pull/178

## [0.21.0] - 2018-03-12
### Updated
- Upgraded Plotly.js, the underlying library behind the
`dash_core_components.Graph` component, to version 1.35.2 [#174](https://github.com/plotly/dash-core-components/pull/174). This includes [Feature release 1.35.0](https://github.com/plotly/plotly.js/releases/tag/v1.35.0) and patch releases [1.35.2](https://github.com/plotly/plotly.js/releases/tag/v1.35.2), [1.35.1](https://github.com/plotly/plotly.js/releases/tag/v1.35.1), which contain:
  - An `automargin` attribute for cartesian axes with long tick labels
  - Support for typed arrays as data inputs
  - layout `grid` attribute for easy subplot generation
  - And numerous other features and bug fixes

## [0.20.2] - 2018-03-05
### Fixed
- The `selectedData`, `clickData`, and `hoverData` callbacks were being attached without being
removed every time the graph was updated. They are now removed and reattached. #172

## [0.20.1] - 2018-03-01
### Fixed
- The `serve_locally` was broken - the Plotly.js bundle wasn't being served correctly.

## [0.20.0] - 2018-03-01
### Updated
- Upgraded plotly.js, the underlying library behind the
`dash_core_components.Graph` component, to version 1.34.0 [#170](https://github.com/plotly/dash-core-components/pull/170). [Feature release 1.34.0](https://github.com/plotly/plotly.js/releases/tag/v1.34.0) contains:
  - `Plotly.react`, a new do-it-all API method that creates and updates graphs using the same API signature
  - Constraint-type contours in `contour` traces
  - `notched` `box` traces
  - Localization machinery for auto-formatted date axis ticks
  - And many other features and bug fixes

## [0.19.0] - 2018-02-11
### Changed
- `PropTypes` now uses `prop-types` package instead of `React` to support move to React 16+

## [0.18.1] - 2017-01-25
### Fixed
- Patched plotly.js to version 1.33.1 [#151](https://github.com/plotly/dash-core-components/pull/151). [Patch release 1.33.1](https://github.com/plotly/plotly.js/releases/tag/v1.33.1) includes numerous bug fixes.

## [0.18.0] - 2017-01-19
### Updated
- Upgraded Plotly.js, the underlying library behind the `dash_core_components.Graph` component, to [version 1.33.0](https://github.com/plotly/plotly.js/releases/tag/v1.33.0). This was a huge release! Here are some of the new features that are available:
  - Completely rewritten `scattergl` trace type using `regl` [plotly.js#2258](https://github.com/plotly/plotly.js/pull/2258)
  - Completely rewritten polar chart renderer accompanied by new `scatterpolar` and `scatterpolargl` trace types [plotly.js#2200](https://github.com/plotly/plotly.js/pull/2200). The old polar trace types - `scatter` with `(r, t)` coordinates, bar with `(r, t)` coordinates, and `area` - are now deprecated.
  - Add the ability to draw layout images and layout shapes on subplot with `scattergl` traces [plotly.js#2258](https://github.com/plotly/plotly.js/pull/2258)
  - Add `fill` capabilities to `scattergl` traces [plotly.js#2258](https://github.com/plotly/plotly.js/pull/2258)
  - Add `spikedistance`, `hoverdistance` and `spikesnap` for more customizable spikes and hover behavior on cartesian subplots [plotly.js#2247](https://github.com/plotly/plotly.js/pull/2247)
  - Add official Spanish translation (locale `es`) [plotly.js#2249](https://github.com/plotly/plotly.js/pull/2249)
  - Add official French translation (locale `fr`) [plotly.js#2252](https://github.com/plotly/plotly.js/pull/2252)
  - And numerous other features and bug fixes
  - Many of these features were funded directly by companies that rely on this library. If your organization or company would like to sponsor particular features or bug fixes in these open source libraries, please reach out: http://plotly.com/products/consulting-and-oem

## [0.17.1] - 2017-01-18
### Fixed
- Previously, if `None` is supplied to `SyntaxHighlighter` or `Markdown`, the
component would not render and the app would break. This is problematic because
if `children` isn't supplied (as done in the case for when you are updating that
property from a callback), `None` is the default property. Fixes https://github.com/plotly/dash-core-components/issues/147. This bug was introduced in
v0.15.4.

## [0.17.0] - 2017-01-11
### Added
- The `dcc.Graph` component now includes `pointNumbers` inside `selectedData`
and `hoverData` if the chart type is a `histogram`, `histogram2d`, or `histogram2dcontour`.

## [0.16.0] - 2017-01-11
### Updated
- Upgraded Plotly.js, the underlying library behind the `dash_core_components.Graph` component, to version 1.32.0 [#143](https://github.com/plotly/dash-core-components/pull/143). [Feature release 1.32.0](https://github.com/plotly/plotly.js/releases/tag/v1.32.0) includes:
  - Localization machinery, including an official German translation (locale `de`)
  - A new `violin` trace type
  - Selection improvements: `selected` and `unselected` attribute containers to customize selection states, Support for multi-selections, persistent selections.
  - `colorway` attribute to customize the trace-to-trace color sequence
  - Add `tickformatstops` to set tick format per cartesian axis range
  - And many more features and bug fixes

## [0.15.5] - 2017-01-08
### Fixed
- The `dash_core_components.Location` and `dash_core_components.Link` properties
should now work on Internet Explorer.
Thanks to @nedned for suggesting a solution.
Fixes https://github.com/plotly/dash-core-components/pull/113

## [0.15.4] - 2017-12-21
### Changed
- The `dash_core_components.Location` component now supports `hash`,
`href`, and `search` in addition to the already supported `pathname`
(mimicking the `window.location` API). `href` can be used to handle
`pathname`, `hash`, and `search` in aggregate, or each can be manipulated
independently.
- The `children` property of `dash_core_components.Markdown` and
`dash_core_components.SyntaxHighlighter` now accepts an
array of strings (previously it *had* to be a string). Now,
if an array is provided, it is collapsed into a string with line
breaks (see #134).

## [0.15.3] - 2017-12-11
### Fixed
- Patched plotly.js to version 1.31.2 [#125](https://github.com/plotly/dash-core-components/pull/125), including patch releases [1.31.2](https://github.com/plotly/plotly.js/releases/tag/v1.31.2) and [1.31.1](https://github.com/plotly/plotly.js/releases/tag/v1.31.1) with numerous bug fixes.


## [0.15.2] - 2017-11-24
### :sweat_smile: Added
- The `Interval` component has a new property: `n_intervals`. This is an
integer that increases every time that the interval passes. This allows you
to use the `Interval` component without using the `events=[Event(...)]` pattern
inside the callback.

This is similar to the `n_clicks` property of the `dash_html_components`
components.
This was the last use case for `events=[Event(...)]` inside the
`dash_core_components` library. Ultimately, we may be able to deprecate this
pattern.

### Changed
- The `dash_core_components.Input(type='number')` component actually converts
the values to floats or integers, instead of passing the numbers back as strings.
https://github.com/plotly/dash-core-components/pull/100
Big thanks to community contributor @Madhu94!

### Fixed
- The `disable_click` property in the `dcc.Upload` component now works.
https://github.com/plotly/dash-core-components/pull/106.
Big thanks to community contributor @Akronix!
- Several properties in several components had the wrong `propTypes`.
This has been fixed, improving the documentation for the Dash python classes
(and removing warnings in JS development).
Big thanks to community contributor @Akronix!

## [0.15.1] - 2017-11-23
### Fixed
- Attempt to fix the JS builds from 0.15.0 but actually nothing changed.

## [0.15.0] - 2017-11-19
- Bad build. See 0.15.2 for the correct build

## [0.14.0] - 2017-10-17
### :sparkles: Added
- An `Upload` component! :tada: See [https://plotly.com/dash/dash-core-components/upload](https://plotly.com/dash/dash-core-components/upload) for docs.

## [0.13.0] - 2017-10-05
### Updated
- Upgraded plotly.js to version 1.31.0. This includes:
  - Two huge feature releases by the plotly.js team! :clap: - [1.31.0](https://github.com/plotly/plotly.js/releases/tag/v1.31.0) and [1.30.0](https://github.com/plotly/plotly.js/releases/tag/v1.30.0), which contain:
    - A new `table` trace type
    - `geo.center` making geo views fully reproducible
    - Extend lasso and select-box drag modes to `bar`, `histogram`, `scattergeo`, and `choropleth` trace types
    - `aggregate` transforms
    - And many other features and bug fixes
  - Patch release [1.30.1](https://github.com/plotly/plotly.js/releases/tag/v1.30.1) with bug fixes

## [0.12.7] - 2017-09-26
### :bug: Fixed
- Fixed issues related to updating the `max_date_allowed` property of `DatePickerSingle` and `DatePickerRange`  programmatically through callbacks
- Clicking on the end date in the `DatePickerRange` will now open up the calendar to the end date (https://github.com/plotly/dash-core-components/issues/80)

### Maintenance
- Cleaned up `DatePickerSingle` and `DatePickerRange`

## [0.12.6] - 2017-09-11
### :bug: Fixed
-  Non-ascii characters, like chinese characters, are now supported as
   search strings in the `dcc.Dropdown` component (https://github.com/plotly/dash-core-components/pull/75)

## [0.12.5] - 2017-09-11
### :bug: Fixed
-  The `Interval` component was constantly resetting its interval on every update. Initially reported in https://community.plotly.com/t/multiple-interval-object-in-a-single-page/5699/3
- Removed the used `label` property from the `Slider` component
- Provide a more descriptive documentation for the `marks` property of the `Slider` component

### :stars: Added
- A `disabled` property on the `Interval` component will disable the interval component from firing its updates.

## [0.12.4] - 2017-08-18
### Added
- Added `className` and `style` properties to the parent `div`s of the `Checklist`, `Dropdown`, `Graph` and `RadioItems` component. As requested in https://github.com/plotly/dash-core-components/issues/57, solved in https://github.com/plotly/dash-core-components/pull/60

## [0.12.3] - 2017-08-17
### Fixed
- Previously, the `max_date_allowed` could not be selected. This issue has been fixed, issue first reported in https://community.plotly.com/t/solved-datepicker-in-dash/4816/10

## [0.12.2] - 2017-08-10
### Fixed
- Previously, when the `options` of a `dcc.Dropdown` would change, the options would no longer be searchable. That has been fixed. Issue was originally reported in https://community.plotly.com/t/dropdown-not-searching-values-when-typing/5323/3

## [0.12.1] - 2017-08-09
### Fixed
- Disabled portal settings on `dcc.DatePickerSingle` and `dcc.DatePickerRange` when `vertical=True`. `with_portal` and `with_full_screen_portal` will only apply if `vertical=False`.

## [0.12.0] - 2017-08-09
### Added
- Added two new date picker components: `dcc.DatePickerSingle` and `dcc.DatePickerRange`

## [0.11.1] - 2017-08-07
### Fixed
- Added support for all of the valid HTML attributes of the `Input` component.
- Added support for a few more `type` values of the `Input` component. The
  full list of valid types are 'text', 'number', 'password', 'email', 'range', 'search', 'tel', 'url', 'hidden'.
  Note that type values that don't have cross-browser support are not included (such as `datetime`)

## [0.11.0] - 2017-08-04
### Added
- The `Dropdown` component renders `options` much, much faster. It can render 50,000 options (client-side) without crashing! This fixes https://github.com/plotly/dash/issues/103

## [0.10.0] - 2017-08-03
### Updated
- Upgraded plotly.js to version 1.29.3, including feature releases [1.29.0](https://github.com/plotly/plotly.js/releases/tag/v1.29.0) and [1.28.0](https://github.com/plotly/plotly.js/releases/tag/v1.28.0), and patch releases [1.29.3](https://github.com/plotly/plotly.js/releases/tag/v1.29.3), [1.29.2](https://github.com/plotly/plotly.js/releases/tag/v1.29.2), [1.29.1](https://github.com/plotly/plotly.js/releases/tag/v1.29.1), [1.28.2](https://github.com/plotly/plotly.js/releases/tag/v1.28.2), [1.28.1](https://github.com/plotly/plotly.js/releases/tag/v1.28.1), and [1.27.1](https://github.com/plotly/plotly.js/releases/tag/v1.27.1) This includes TONS of fixes and improvements, notably:
    - Add touch interactions to cartesian, gl2d and ternary subplots including for
    select and lasso drag modes
    - Add support for contour line labels in contour and contourcarpet traces
    - Add support for select and lasso drag modes on scattermapbox traces
    - Add reset view and toggle hover mode bar buttons to mapbox subplots
    - Add support for array marker.opacity settings in scattermapbox traces
    - Add namelength layout and trace attribute to control the trace name's
    visible length in hover labels
    - Add cliponaxis attribute to scatter and scatterternary traces to allow
    markers and text nodes to be displayed above their subplot's axes
    - Add axis layer attribute with 'above traces' and 'below traces' values

## [0.9.0] - 2017-07-28
### Added
- A `config` property of the `Graph` component that exposes the [plotly.js config properties](https://plotly.com/javascript/configuration-options/). Here's an example that hides 2 buttons and makes the elements in the graph "editable":
```
import dash
import dash_core_components as dcc
import dash_html_components as html

app = dash.Dash()

app.layout = html.Div([
    dcc.Graph(
        id='my-graph',
        figure={'data': [{'x': [1, 2, 3]}]},
        config={'editable': True, 'modeBarButtonsToRemove': ['pan2d', 'lasso2d']}
    )
])

if __name__ == '__main__':
    app.run_server(debug=True)
```

## [0.8.0] - 2017-07-27
### Added
- A new `Textarea` component for displaying the simple Textarea HTML element. The content of the `Textarea` is controlled through the `value` property:
```
dcc.Textarea(id='my-text-area' value='''
SELECT * FROM MY_TABLES
LIMIT 10;
''')
```

## [0.7.1] - 2017-07-24
### Fixed
- Clearing a Graph selection box sets the `selectedData` value to `None` (`null` in JavaScript). Before, it didn't change the `selectedData` property, preventing the user and the Dash developer from clearing selections. Fixes https://github.com/plotly/dash/issues/97, thanks to @pmbaumgartner for reporting.


## [0.7.0] - 2017-07-20
### Added
- The `clearable` property to the `Dropdown`, which toggles on and off the "x" on the side of the dropdown that clears the current selection.
- The `searchable` property to the `Dropdown`, which toggles on and off whether the `Dropdown` is searchable.

### Fixed
- Clicking on the little `x` on the side of the Dropdown to clear the currently selected value didn't work. Now it does. If `multi=false`, then `null` (or Python's `None`) is set. If `multi=True`, then `[]` is set.

## [0.6.0] - 2017-07-18
### Added
- The `Slider` and the `RangeSlider` component can update when the user finishes dragging the slider rather than just while they drag. The default behaviour has remained the same (updates while dragging) but you can toggle that the updates only get fired on "mouse up" by setting `updatemode` to `'mouseup'` (`'drag'` is the default).
- A `Link` and `Location` were added. `Location` represents the address bar of the web browser and `Link` provides a way to modify the address bar without refreshing the page. Combined, these two components can be used to create a "single page app" with multiple URLs. That is, apps that have multiple URLs but surfing between the different pages doesn't trigger a full page refresh like it would with traditional links.
- Previously, if callback functions weren't supplied to a component, it wouldn't update. This caused a lot of confusion: users would create a simple layout without any callbacks and then wonder why the sliders wouldn't slide or the text inputs wouldn't update. Now, all of the components manage their own state and their appearance will update regardless of whether Dash has assigned a callback to them.

## [0.5.3] - 2017-07-03
### Added
- A `range` object is now included in the `selectedData` object that specifies
  that dimensions of the selected region.
- A `lassoPoints` object is now included in the `selectedData` object that
  provides coordinates of the lassoed region.

## [0.5.2] - 2017-07-03
### Added
- A new property `clear_on_unhover` on the `Graph` component will clear the
  `hoverData` property when the user "unhovers" from a point if True. If False,
  then the `hoverData` property will be equal to the data from the last point
  that was hovered over. The default is False.
# Dash Core Components

This package provides the core React component suite for [Dash][].

[![CircleCI](https://circleci.com/gh/plotly/dash-core-components.svg?style=svg)](https://circleci.com/gh/plotly/dash-core-components)

## Development

This package is part of `dash`, and if you install `dash` in development mode with extras as below, you can develop in this portion as well.
From the root of the `dash` repo:

```bash
# It's recommended to install your python packages in a virtualenv
# As of dash 2.0, python 3 is required
$ python -m venv venv && . venv/bin/activate

# make sure dash is installed with dev and testing dependencies
$ pip install -e .[dev,testing]  # in some shells you need \ to escape []

# run the build process - this will build all of dash, including dcc
$ npm i && npm run build

# install dcc in editable mode
$ pip install -e .
```

### Code quality and tests

### To run integration tests (test_integration.py)
You can run the Selenium integration tests with the
```sh
npm test
```

### Testing your components in Dash
1. Run the build watcher by running
        $ npm run build:watch

2. Run the dash layout you want to test

        # Import dash_core_components to your layout, then run it:
        $ python my_dash_layout.py

## Dash Component Boilerplate

See the [dash-component-boilerplate](https://github.com/plotly/dash-component-boilerplate) repo for more information.

[Dash]: https://plotly.com/dash
[Dash Component Boilerplate]: (https://github.com/plotly/dash-component-boilerplate)
[NPM package authors]: https://www.npmjs.com/package/dash-core-components/access
[PyPi]: https://pypi.python.org/pypi


## Big Thanks
Cross-browser Testing Powered by [![image](https://user-images.githubusercontent.com/1394467/64290307-e4c66600-cf33-11e9-85a1-12c82230a597.png)](https://saucelabs.com)
# Contributing to dash-core-components

## Getting Started

Refer to the [readme](README.md) for installation and development instructions.

## Contributions

[Dash Core Components][] consist of pluggable components for creating interactive user interfaces. For generic HTML5 elements, see [Dash HTML Components][]. Contributions are welcome! This repository's open [issues][] are a good place to start. Another way to contribute is to [write your own components][] using for instance the [component boilerplate](https://github.com/plotly/dash-component-boilerplate).

## Coding Style

Please lint any additions to react components with `npm run lint`. Rules defined in [.eslintrc](.eslintrc) are inherited from [`dash-components-archetype`](https://github.com/plotly/dash-components-archetype)'s [eslintrc-react.json][]

## Pull Request Guidelines

Use the [GitHub flow][] when proposing contributions to this repository (i.e. create a feature branch and submit a PR against the master branch).

## Making a Contribution
_For larger features, your contribution will have a higher likelihood of getting merged if you create an issue to discuss the changes that you'd like to make before you create a pull request._

1. Create a pull request.
2. After a review has been done and your changes have been approved, they will be merged and included in a future release of Dash.
3. If significant enough, you have created an issue about documenting the new feature or change and you have added it to the [dash-docs](https://github.com/plotly/dash-docs) project.

## Running the Tests

In order to run the tests, you first need to have built the JavaScript
`dash_core_components` library. You will need to build the library again if
you've pulled from upstream otherwise you may be running with an out of date
`bundle.js`. See the instructions for building `bundle.js` in the [Testing
Locally](README.md#testing-locally) section of README.md.

## Updating Plotly.js

1. Update the version of `plotly.js` in package.json. Always use an exact version without "^" or "~"
2. Run `npm install` followed by `npm run build`, the Plotly.js artifact will be copied and bundled over as required
4. Update `CHANGELOG.md` with links to the releases and a description of the changes. The message should state (see the existing `CHANGELOG.md` for examples):
    * If you're only bumping the patch level, the heading is "Fixed" and the text starts "Patched plotly.js". Otherwise the heading is "Updated" and the text starts "Upgraded plotly.js"
    * The new plotly.js version number, and the PR in which this was done
    * All major or minor versions included, with links to their release pages and a summary of the major new features in each. If there are multiple minor/major releases included, be sure to look at all of their release notes to construct the summary. Call minor versions "feature" versions for the benefit of users not steeped in semver terminology.
    * All patch versions included, with links to their release pages and a note that these fix bugs
5. When bumping the dcc version, a plotly.js patch/minor/major constitutes a dcc patch/minor/major respectively as well.

## Financial Contributions

If your company wishes to sponsor development of open source dash components, please [get in touch][].

[Dash Core Components]: https://dash.plotly.com/dash-core-components
[Dash HTML Components]: https://github.com/plotly/dash-html-components
[write your own components]: https://dash.plotly.com/plugins
[Dash Components Archetype]: https://github.com/plotly/dash-components-archetype
[issues]: https://github.com/plotly/dash-core-components/issues
[GitHub flow]: https://guides.github.com/introduction/flow/
[eslintrc-react.json]: https://github.com/plotly/dash-components-archetype/blob/master/config/eslint/eslintrc-react.json
[contributors]: https://github.com/plotly/dash-core-components/graphs/contributors
[semantic versioning]: https://semver.org/
[Dash Community Forum]: https://community.plotly.com/c/dash
[Confirmation Modal component]: https://github.com/plotly/dash-core-components/pull/211#issue-195280462
[Confirmation Modal announcement]: https://community.plotly.com/t/announcing-dash-confirmation-modal-feedback-welcome/11627
[get in touch]: https://plotly.com/products/consulting-and-oem
Empty# @plotly/prettier-config-dash

Shareable Prettier configuration for Dash projects.# @plotly/eslint-config-dash

Shareable linting configurations for Dash projects.EmptyEmpty