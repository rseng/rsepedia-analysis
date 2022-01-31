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

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at xarray-core-team@googlegroups.com. The project team will review and investigate all complaints, and will respond in a way that it deems appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
Xarray's contributor guidelines [can be found in our online documentation](http://xarray.pydata.org/en/stable/contributing.html)
# How to issue an xarray release in 16 easy steps

Time required: about an hour.

These instructions assume that `upstream` refers to the main repository:

```sh
$ git remote -v
{...}
upstream        https://github.com/pydata/xarray (fetch)
upstream        https://github.com/pydata/xarray (push)
```

<!-- markdownlint-disable MD031 -->

 1. Ensure your main branch is synced to upstream:
     ```sh
     git switch main
     git pull upstream main
     ```
 2. Confirm there are no commits on stable that are not yet merged
    ([ref](https://github.com/pydata/xarray/pull/4440)):
     ```sh
     git merge upstream/stable
     ```
 3. Add a list of contributors with:
    ```sh
    git log "$(git tag --sort=v:refname | tail -1).." --format=%aN | sort -u | perl -pe 's/\n/$1, /'
    ```
    This will return the number of contributors:
    ```sh
    git log "$(git tag --sort=v:refname | tail -1).." --format=%aN | sort -u | wc -l
    ```
 4. Write a release summary: ~50 words describing the high level features. This
    will be used in the release emails, tweets, GitHub release notes, etc.
 5. Look over whats-new.rst and the docs. Make sure "What's New" is complete
    (check the date!) and add the release summary at the top.
    Things to watch out for:
    - Important new features should be highlighted towards the top.
    - Function/method references should include links to the API docs.
    - Sometimes notes get added in the wrong section of whats-new, typically
      due to a bad merge. Check for these before a release by using git diff,
      e.g., `git diff v{0.X.Y-1} whats-new.rst` where {0.X.Y-1} is the previous
      release.
 6. Open a PR with the release summary and whatsnew changes; in particular the
    release headline should get feedback from the team on what's important to include.
 7. After merging, again ensure your main branch is synced to upstream:
     ```sh
     git pull upstream main
     ```
 8. If you have any doubts, run the full test suite one final time!
      ```sh
      pytest
      ```
 9. Check that the ReadTheDocs build is passing.
10. Issue the release on GitHub. Click on "Draft a new release" at
    <https://github.com/pydata/xarray/releases>. Type in the version number (with a "v")
    and paste the release summary in the notes.
11. This should automatically trigger an upload of the new build to PyPI via GitHub Actions.
    Check this has run [here](https://github.com/pydata/xarray/actions/workflows/pypi-release.yaml),
    and that the version number you expect is displayed [on PyPI](https://pypi.org/project/xarray/)
12. Update the stable branch (used by ReadTheDocs) and switch back to main:
     ```sh
      git switch stable
      git rebase main
      git push --force upstream stable
      git switch main
     ```
    You may need to first fetch it with `git fetch upstream`,
    and check out a local version with `git checkout -b stable upstream/stable`.

    It's OK to force push to `stable` if necessary. (We also update the stable
    branch with `git cherry-pick` for documentation only fixes that apply the
    current released version.)
13. Add a section for the next release {0.X.Y+1} to doc/whats-new.rst:
     ```rst
     .. _whats-new.0.X.Y+1:

     v0.X.Y+1 (unreleased)
     ---------------------

     New Features
     ~~~~~~~~~~~~


     Breaking changes
     ~~~~~~~~~~~~~~~~


     Deprecations
     ~~~~~~~~~~~~


     Bug fixes
     ~~~~~~~~~


     Documentation
     ~~~~~~~~~~~~~


     Internal Changes
     ~~~~~~~~~~~~~~~~

     ```
14. Commit your changes and push to main again:
      ```sh
      git commit -am 'New whatsnew section'
      git push upstream main
      ```
    You're done pushing to main!

15. Update the docs. Login to <https://readthedocs.org/projects/xray/versions/>
    and switch your new release tag (at the bottom) from "Inactive" to "Active".
    It should now build automatically.
16. Issue the release announcement to mailing lists & Twitter. For bug fix releases, I
    usually only email xarray@googlegroups.com. For major/feature releases, I will email a broader
    list (no more than once every 3-6 months):
      - pydata@googlegroups.com
      - xarray@googlegroups.com
      - numpy-discussion@scipy.org
      - scipy-user@scipy.org
      - pyaos@lists.johnny-lin.com

    Google search will turn up examples of prior release announcements (look for
    "ANN xarray").
    Some of these groups require you to be subscribed in order to email them.

<!-- markdownlint-enable MD013 -->

## Note on version numbering

We follow a rough approximation of semantic version. Only major releases (0.X.0)
should include breaking changes. Minor releases (0.X.Y) are for bug fixes and
backwards compatible new features, but if a sufficient number of new features
have arrived we will issue a major release even if there are no compatibility
breaks.

Once the project reaches a sufficient level of maturity for a 1.0.0 release, we
intend to follow semantic versioning more strictly.
<!-- Feel free to remove check-list items aren't relevant to your change -->

- [ ] Closes #xxxx
- [ ] Tests added
- [ ] User visible changes (including notable bug fixes) are documented in `whats-new.rst`
- [ ] New functions/methods are listed in `api.rst`
# Property-based tests using Hypothesis

This directory contains property-based tests using a library
called [Hypothesis](https://github.com/HypothesisWorks/hypothesis-python).

The property tests for xarray are a work in progress - more are always welcome.
They are stored in a separate directory because they tend to run more examples
and thus take longer, and so that local development can run a test suite
without needing to `pip install hypothesis`.

## Hang on, "property-based" tests?

Instead of making assertions about operations on a particular piece of
data, you use Hypothesis to describe a *kind* of data, then make assertions
that should hold for *any* example of this kind.

For example: "given a 2d ndarray of dtype uint8 `arr`,
`xr.DataArray(arr).plot.imshow()` never raises an exception".

Hypothesis will then try many random examples, and report a minimised
failing input for each error it finds.
[See the docs for more info.](https://hypothesis.readthedocs.io/en/master/)
# Benchmark CI

<!-- Author: @jaimergp -->
<!-- Last updated: 2021.07.06 -->
<!-- Describes the work done as part of https://github.com/scikit-image/scikit-image/pull/5424 -->

## How it works

The `asv` suite can be run for any PR on GitHub Actions (check workflow `.github/workflows/benchmarks.yml`) by adding a `run-benchmark` label to said PR. This will trigger a job that will run the benchmarking suite for the current PR head (merged commit) against the PR base (usually `main`).

We use `asv continuous` to run the job, which runs a relative performance measurement. This means that there's no state to be saved and that regressions are only caught in terms of performance ratio (absolute numbers are available but they are not useful since we do not use stable hardware over time). `asv continuous` will:

* Compile `scikit-image` for _both_ commits. We use `ccache` to speed up the process, and `mamba` is used to create the build environments.
* Run the benchmark suite for both commits, _twice_  (since `processes=2` by default).
* Generate a report table with performance ratios:
    * `ratio=1.0` -> performance didn't change.
    * `ratio<1.0` -> PR made it slower.
    * `ratio>1.0` -> PR made it faster.

Due to the sensitivity of the test, we cannot guarantee that false positives are not produced. In practice, values between `(0.7, 1.5)` are to be considered part of the measurement noise. When in doubt, running the benchmark suite one more time will provide more information about the test being a false positive or not.

## Running the benchmarks on GitHub Actions

1. On a PR, add the label `run-benchmark`.
2. The CI job will be started. Checks will appear in the usual dashboard panel above the comment box.
3. If more commits are added, the label checks will be grouped with the last commit checks _before_ you added the label.
4. Alternatively, you can always go to the `Actions` tab in the repo and [filter for `workflow:Benchmark`](https://github.com/scikit-image/scikit-image/actions?query=workflow%3ABenchmark). Your username will be assigned to the `actor` field, so you can also filter the results with that if you need it.

## The artifacts

The CI job will also generate an artifact. This is the `.asv/results` directory compressed in a zip file. Its contents include:

* `fv-xxxxx-xx/`. A directory for the machine that ran the suite. It contains three files:
    * `<baseline>.json`, `<contender>.json`: the benchmark results for each commit, with stats.
    * `machine.json`: details about the hardware.
* `benchmarks.json`: metadata about the current benchmark suite.
* `benchmarks.log`: the CI logs for this run.
* This README.

## Re-running the analysis

Although the CI logs should be enough to get an idea of what happened (check the table at the end), one can use `asv` to run the analysis routines again.

1. Uncompress the artifact contents in the repo, under `.asv/results`. This is, you should see `.asv/results/benchmarks.log`, not `.asv/results/something_else/benchmarks.log`. Write down the machine directory name for later.
2. Run `asv show` to see your available results. You will see something like this:

```
$> asv show

Commits with results:

Machine    : Jaimes-MBP
Environment: conda-py3.9-cython-numpy1.20-scipy

    00875e67

Machine    : fv-az95-499
Environment: conda-py3.7-cython-numpy1.17-pooch-scipy

    8db28f02
    3a305096
```

3. We are interested in the commits for `fv-az95-499` (the CI machine for this run). We can compare them with `asv compare` and some extra options. `--sort ratio` will show largest ratios first, instead of alphabetical order. `--split` will produce three tables: improved, worsened, no changes. `--factor 1.5` tells `asv` to only complain if deviations are above a 1.5 ratio. `-m` is used to indicate the machine ID (use the one you wrote down in step 1). Finally, specify your commit hashes: baseline first, then contender!

```
$> asv compare --sort ratio --split --factor 1.5 -m fv-az95-499 8db28f02 3a305096

Benchmarks that have stayed the same:

       before           after         ratio
     [8db28f02]       [3a305096]
     <ci-benchmark-check~9^2>
              n/a              n/a      n/a  benchmark_restoration.RollingBall.time_rollingball_ndim
      1.23±0.04ms       1.37±0.1ms     1.12  benchmark_transform_warp.WarpSuite.time_to_float64(<class 'numpy.float64'>, 128, 3)
       5.07±0.1μs       5.59±0.4μs     1.10  benchmark_transform_warp.ResizeLocalMeanSuite.time_resize_local_mean(<class 'numpy.float32'>, (192, 192, 192), (192, 192, 192))
      1.23±0.02ms       1.33±0.1ms     1.08  benchmark_transform_warp.WarpSuite.time_same_type(<class 'numpy.float32'>, 128, 3)
       9.45±0.2ms       10.1±0.5ms     1.07  benchmark_rank.Rank3DSuite.time_3d_filters('majority', (32, 32, 32))
       23.0±0.9ms         24.6±1ms     1.07  benchmark_interpolation.InterpolationResize.time_resize((80, 80, 80), 0, 'symmetric', <class 'numpy.float64'>, True)
         38.7±1ms         41.1±1ms     1.06  benchmark_transform_warp.ResizeLocalMeanSuite.time_resize_local_mean(<class 'numpy.float32'>, (2048, 2048), (192, 192, 192))
       4.97±0.2μs       5.24±0.2μs     1.05  benchmark_transform_warp.ResizeLocalMeanSuite.time_resize_local_mean(<class 'numpy.float32'>, (2048, 2048), (2048, 2048))
       4.21±0.2ms       4.42±0.3ms     1.05  benchmark_rank.Rank3DSuite.time_3d_filters('gradient', (32, 32, 32))

...
```

If you want more details on a specific test, you can use `asv show`. Use `-b pattern` to filter which tests to show, and then specify a commit hash to inspect:

```
$> asv show -b time_to_float64 8db28f02

Commit: 8db28f02 <ci-benchmark-check~9^2>

benchmark_transform_warp.WarpSuite.time_to_float64 [fv-az95-499/conda-py3.7-cython-numpy1.17-pooch-scipy]
  ok
  =============== ============= ========== ============= ========== ============ ========== ============ ========== ============
  --                                                                N / order
  --------------- --------------------------------------------------------------------------------------------------------------
      dtype_in       128 / 0     128 / 1      128 / 3     1024 / 0    1024 / 1    1024 / 3    4096 / 0    4096 / 1    4096 / 3
  =============== ============= ========== ============= ========== ============ ========== ============ ========== ============
    numpy.uint8    2.56±0.09ms   523±30μs   1.28±0.05ms   130±3ms     28.7±2ms    81.9±3ms   2.42±0.01s   659±5ms    1.48±0.01s
    numpy.uint16   2.48±0.03ms   530±10μs   1.28±0.02ms   130±1ms    30.4±0.7ms   81.1±2ms    2.44±0s     653±3ms    1.47±0.02s
   numpy.float32    2.59±0.1ms   518±20μs   1.27±0.01ms   127±3ms     26.6±1ms    74.8±2ms   2.50±0.01s   546±10ms   1.33±0.02s
   numpy.float64   2.48±0.04ms   513±50μs   1.23±0.04ms   134±3ms     30.7±2ms    85.4±2ms   2.55±0.01s   632±4ms    1.45±0.01s
  =============== ============= ========== ============= ========== ============ ========== ============ ========== ============
  started: 2021-07-06 06:14:36, duration: 1.99m
```

## Other details

### Skipping slow or demanding tests

To minimize the time required to run the full suite, we trimmed the parameter matrix in some cases and, in others, directly skipped tests that ran for too long or require too much memory. Unlike `pytest`, `asv` does not have a notion of marks. However, you can `raise NotImplementedError` in the setup step to skip a test. In that vein, a new private function is defined at `benchmarks.__init__`: `_skip_slow`. This will check if the `ASV_SKIP_SLOW` environment variable has been defined. If set to `1`, it will raise `NotImplementedError` and skip the test. To implement this behavior in other tests, you can add the following attribute:

```python
from . import _skip_slow  # this function is defined in benchmarks.__init__

def time_something_slow():
    pass

time_something.setup = _skip_slow
```
# Proposal: Xarray flexible indexes refactoring

Current status: https://github.com/pydata/xarray/projects/1

## 1. Data Model

Indexes are used in Xarray to extract data from Xarray objects using coordinate labels instead of using integer array indices. Although the indexes used in an Xarray object can be accessed (or built on-the-fly) via public methods like `to_index()` or properties like `indexes`, those are mainly used internally.

The goal of this project is to make those indexes 1st-class citizens of Xarray's data model. As such, indexes should clearly be separated from Xarray coordinates with the following relationships:

- Index -> Coordinate: one-to-many
- Coordinate -> Index: one-to-zero-or-one

An index may be built from one or more coordinates. However, each coordinate must relate to one index at most. Additionally, a coordinate may not be tied to any index.

The order in which multiple coordinates relate to an index may matter. For example, Scikit-Learn's [`BallTree`](https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.BallTree.html#sklearn.neighbors.BallTree) index with the Haversine metric requires providing latitude and longitude values in that specific order. As another example, the order in which levels are defined in a `pandas.MultiIndex` may affect its lexsort depth (see [MultiIndex sorting](https://pandas.pydata.org/pandas-docs/stable/user_guide/advanced.html#sorting-a-multiindex)).

Xarray's current data model has the same index-coordinate relationships than stated above, although this assumes that multi-index "virtual" coordinates are counted as coordinates (we can consider them as such, with some constraints). More importantly, This refactoring would turn the current one-to-one relationship between a dimension and an index into a many-to-many relationship, which would overcome some current limitations.

For example, we might want to select data along a dimension which has several coordinates:

```python
>>> da
<xarray.DataArray (river_profile: 100)>
array([...])
Coordinates:
  * drainage_area  (river_profile) float64 ...
  * chi            (river_profile) float64 ...
```

In this example, `chi` is a transformation of the `drainage_area` variable that is often used in geomorphology. We'd like to select data along the river profile using either `da.sel(drainage_area=...)` or `da.sel(chi=...)` but that's not currently possible. We could rename the `river_profile` dimension to one of the coordinates, then use `sel` with that coordinate, then call `swap_dims` if we want to use `sel` with the other coordinate, but that's not ideal. We could also build a `pandas.MultiIndex` from `drainage_area` and `chi`, but that's not optimal (there's no hierarchical relationship between these two coordinates).

Let's take another example:

```python
>>> da
<xarray.DataArray (x: 200, y: 100)>
array([[...], [...]])
Coordinates:
  * lon      (x, y) float64 ...
  * lat      (x, y) float64 ...
  * x        (x) float64 ...
  * y        (y) float64 ...
```

This refactoring would allow creating a geographic index for `lat` and `lon` and two simple indexes for `x` and `y` such that we could select data with either `da.sel(lon=..., lat=...)` or  `da.sel(x=..., y=...)`.

Refactoring the dimension -> index one-to-one relationship into many-to-many would also introduce some issues that we'll need to address, e.g., ambiguous cases like `da.sel(chi=..., drainage_area=...)` where multiple indexes may potentially return inconsistent positional indexers along a dimension.

## 2. Proposed API changes

### 2.1 Index wrapper classes

Every index that is used to select data from Xarray objects should inherit from a base class, e.g., `XarrayIndex`, that provides some common API. `XarrayIndex` subclasses would generally consist of thin wrappers around existing index classes such as `pandas.Index`, `pandas.MultiIndex`, `scipy.spatial.KDTree`, etc.

There is a variety of features that an xarray index wrapper may or may not support:

- 1-dimensional vs. 2-dimensional vs. n-dimensional coordinate (e.g., `pandas.Index` only supports 1-dimensional coordinates while a geographic index could be built from n-dimensional coordinates)
- built from a single vs multiple coordinate(s) (e.g., `pandas.Index` is built from one coordinate, `pandas.MultiIndex` may be built from an arbitrary number of coordinates and a geographic index would typically require two latitude/longitude coordinates)
- in-memory vs. out-of-core (dask) index data/coordinates (vs. other array backends)
- range-based vs. point-wise selection
- exact vs. inexact lookups

Whether or not a `XarrayIndex` subclass supports each of the features listed above should be either declared explicitly via a common API or left to the implementation. An `XarrayIndex` subclass may encapsulate more than one underlying object used to perform the actual indexing. Such "meta" index would typically support a range of features among those mentioned above and would automatically select the optimal index object for a given indexing operation.

An `XarrayIndex` subclass must/should/may implement the following properties/methods:

- a `from_coords` class method that creates a new index wrapper instance from one or more Dataset/DataArray coordinates (+ some options)
- a `query` method that takes label-based indexers as argument (+ some options) and that returns the corresponding position-based indexers
- an `indexes` property to access the underlying index object(s) wrapped by the `XarrayIndex` subclass
- a `data` property to access index's data and map it to coordinate data (see [Section 4](#4-indexvariable))
- a `__getitem__()` implementation to propagate the index through DataArray/Dataset indexing operations
- `equals()`, `union()` and `intersection()` methods for data alignment (see [Section 2.6](#26-using-indexes-for-data-alignment))
- Xarray coordinate getters (see [Section 2.2.4](#224-implicit-coodinates))
- a method that may return a new index and that will be called when one of the corresponding coordinates is dropped from the Dataset/DataArray (multi-coordinate indexes)
- `encode()`/`decode()` methods that would allow storage-agnostic serialization and fast-path reconstruction of the underlying index object(s) (see [Section 2.8](#28-index-encoding))
- one or more "non-standard" methods or properties that could be leveraged in Xarray 3rd-party extensions like Dataset/DataArray accessors (see [Section 2.7](#27-using-indexes-for-other-purposes))

The `XarrayIndex` API has still to be defined in detail.

Xarray should provide a minimal set of built-in index wrappers (this could be reduced to the indexes currently supported in Xarray, i.e., `pandas.Index` and `pandas.MultiIndex`). Other index wrappers may be implemented in 3rd-party libraries (recommended). The `XarrayIndex` base class should be part of Xarray's public API.

#### 2.1.1 Index discoverability

For better discoverability of Xarray-compatible indexes, Xarray could provide some mechanism to register new index wrappers, e.g., something like [xoak's `IndexRegistry`](https://xoak.readthedocs.io/en/latest/_api_generated/xoak.IndexRegistry.html#xoak.IndexRegistry) or [numcodec's registry](https://numcodecs.readthedocs.io/en/stable/registry.html).

Additionally (or alternatively), new index wrappers may be registered via entry points as is already the case for storage backends and maybe other backends (plotting) in the future.

Registering new indexes either via a custom registry or via entry points should be optional. Xarray should also allow providing `XarrayIndex` subclasses in its API (Dataset/DataArray constructors, `set_index()`, etc.).

### 2.2 Explicit vs. implicit index creation

#### 2.2.1 Dataset/DataArray's `indexes` constructor argument

The new `indexes` argument of Dataset/DataArray constructors may be used to specify which kind of index to bind to which coordinate(s). It would consist of a mapping where, for each item, the key is one coordinate name (or a sequence of coordinate names) that must be given in `coords` and the value is the type of the index to build from this (these) coordinate(s):

```python
>>> da = xr.DataArray(
...     data=[[275.2, 273.5], [270.8, 278.6]],
...     dims=('x', 'y'),
...     coords={
...         'lat': (('x', 'y'), [[45.6, 46.5], [50.2, 51.6]]),
...         'lon': (('x', 'y'), [[5.7, 10.5], [6.2, 12.8]]),
...     },
...     indexes={('lat', 'lon'): SpatialIndex},
... )
<xarray.DataArray (x: 2, y: 2)>
array([[275.2, 273.5],
       [270.8, 278.6]])
Coordinates:
  * lat      (x, y) float64 45.6 46.5 50.2 51.6
  * lon      (x, y) float64 5.7 10.5 6.2 12.8
```

More formally, `indexes` would accept `Mapping[CoordinateNames, IndexSpec]` where:

- `CoordinateNames = Union[CoordinateName, Tuple[CoordinateName, ...]]` and `CoordinateName = Hashable`
- `IndexSpec = Union[Type[XarrayIndex], Tuple[Type[XarrayIndex], Dict[str, Any]], XarrayIndex]`, so that index instances or index classes + build options could be also passed

Currently index objects like `pandas.MultiIndex` can be passed directly to `coords`, which in this specific case results in the implicit creation of virtual coordinates. With the new `indexes` argument this behavior may become even more confusing than it currently is. For the sake of clarity, it would be appropriate to eventually drop support for this specific behavior and treat any given mapping value given in `coords` as an array that can be wrapped into an Xarray variable, i.e., in the case of a multi-index:

```python
>>> xr.DataArray([1.0, 2.0], dims='x', coords={'x': midx})
<xarray.DataArray (x: 2)>
array([1., 2.])
Coordinates:
    x        (x) object ('a', 0) ('b', 1)
```

A possible, more explicit solution to reuse a `pandas.MultiIndex` in a DataArray/Dataset with levels exposed as coordinates is proposed in [Section 2.2.4](#224-implicit-coordinates).

#### 2.2.2 Dataset/DataArray's `set_index` method

New indexes may also be built from existing sets of coordinates or variables in a Dataset/DataArray using the `.set_index()` method.

The [current signature](http://xarray.pydata.org/en/stable/generated/xarray.DataArray.set_index.html#xarray.DataArray.set_index) of `.set_index()` is tailored to `pandas.MultiIndex` and tied to the concept of a dimension-index. It is therefore hardly reusable as-is in the context of flexible indexes proposed here.

The new signature may look like one of these:

- A. `.set_index(coords: CoordinateNames, index: Union[XarrayIndex, Type[XarrayIndex]], **index_kwargs)`: one index is set at a time, index construction options may be passed as keyword arguments
- B. `.set_index(indexes: Mapping[CoordinateNames, Union[Type[XarrayIndex], Tuple[Type[XarrayIndex], Dict[str, Any]]]])`: multiple indexes may be set at a time from a mapping of coordinate or variable name(s) as keys and `XarrayIndex` subclasses (maybe with a dict of build options) as values. If variable names are given as keys of they will be promoted as coordinates

Option A looks simple and elegant but significantly departs from the current signature. Option B is more consistent with the Dataset/DataArray constructor signature proposed in the previous section and would be easier to adopt in parallel with the current signature that we could still support through some depreciation cycle.

The `append` parameter of the current `.set_index()` is specific to `pandas.MultiIndex`. With option B we could still support it, although we might want to either drop it or move it to the index construction options in the future.

#### 2.2.3 Implicit default indexes

In general explicit index creation should be preferred over implicit index creation. However, there is a majority of cases where basic `pandas.Index` objects could be built and used as indexes for 1-dimensional coordinates. For convenience, Xarray should automatically build such indexes for the coordinates where no index has been explicitly assigned in the Dataset/DataArray constructor or when indexes have been reset / dropped.

For which coordinates?

- A. only 1D coordinates with a name matching their dimension name
- B. all 1D coordinates

When to create it?

- A. each time when a new Dataset/DataArray is created
- B. only when we need it (i.e., when calling `.sel()` or `indexes`)

Options A and A are what Xarray currently does and may be the best choice considering that indexes could possibly be invalidated by coordinate mutation.

Besides `pandas.Index`, other indexes currently supported in Xarray like `CFTimeIndex` could be built depending on the coordinate data type.

#### 2.2.4 Implicit coordinates

Like for the indexes, explicit coordinate creation should be preferred over implicit coordinate creation. However, there may be some situations where we would like to keep creating coordinates implicitly for backwards compatibility.

For example, it is currently possible to pass a `pandas.MulitIndex` object as a coordinate to the Dataset/DataArray constructor:

```python
>>> midx = pd.MultiIndex.from_arrays([['a', 'b'], [0, 1]], names=['lvl1', 'lvl2'])
>>> da = xr.DataArray([1.0, 2.0], dims='x', coords={'x': midx})
>>> da
<xarray.DataArray (x: 2)>
array([1., 2.])
Coordinates:
  * x        (x) MultiIndex
  - lvl1     (x) object 'a' 'b'
  - lvl2     (x) int64 0 1
```

In that case, virtual coordinates are created for each level of the multi-index. After the index refactoring, these coordinates would become real coordinates bound to the multi-index.

In the example above a coordinate is also created for the `x` dimension:

```python
>>> da.x
<xarray.DataArray 'x' (x: 2)>
array([('a', 0), ('b', 1)], dtype=object)
Coordinates:
  * x        (x) MultiIndex
  - lvl1     (x) object 'a' 'b'
  - lvl2     (x) int64 0 1
```

With the new proposed data model, this wouldn't be a requirement anymore: there is no concept of a dimension-index. However, some users might still rely on the `x` coordinate so we could still (temporarily) support it for backwards compatibility.

Besides `pandas.MultiIndex`, there may be other situations where we would like to reuse an existing index in a new Dataset/DataArray (e.g., when the index is very expensive to build), and which might require implicit creation of one or more coordinates.

The example given here is quite confusing, though: this is not an easily predictable behavior. We could entirely avoid the implicit creation of coordinates, e.g., using a helper function that generates coordinate + index dictionaries that we could then pass directly to the DataArray/Dataset constructor:

```python
>>> coords_dict, index_dict = create_coords_from_index(midx, dims='x', include_dim_coord=True)
>>> coords_dict
{'x': <xarray.Variable (x: 2)>
 array([('a', 0), ('b', 1)], dtype=object),
 'lvl1': <xarray.Variable (x: 2)>
 array(['a', 'b'], dtype=object),
 'lvl2': <xarray.Variable (x: 2)>
 array([0, 1])}
>>> index_dict
{('lvl1', 'lvl2'): midx}
>>> xr.DataArray([1.0, 2.0], dims='x', coords=coords_dict, indexes=index_dict)
<xarray.DataArray (x: 2)>
array([1., 2.])
Coordinates:
    x        (x) object ('a', 0) ('b', 1)
  * lvl1     (x) object 'a' 'b'
  * lvl2     (x) int64 0 1
```

### 2.2.5 Immutable indexes

Some underlying indexes might be mutable (e.g., a tree-based index structure that allows dynamic addition of data points) while other indexes like `pandas.Index` aren't. To keep things simple, it is probably better to continue considering all indexes in Xarray as immutable (as well as their corresponding coordinates, see [Section 2.4.1](#241-mutable-coordinates)).

### 2.3 Index access

#### 2.3.1 Dataset/DataArray's `indexes` property

The `indexes` property would allow easy access to all the indexes used in a Dataset/DataArray. It would return a `Dict[CoordinateName, XarrayIndex]` for easy index lookup from coordinate name.

#### 2.3.2 Additional Dataset/DataArray properties or methods

In some cases the format returned by the `indexes` property would not be the best (e.g, it may return duplicate index instances as values). For convenience, we could add one more property / method to get the indexes in the desired format if needed.

### 2.4 Propagate indexes through operations

#### 2.4.1 Mutable coordinates

Dataset/DataArray coordinates may be replaced (`__setitem__`) or dropped (`__delitem__`) in-place, which may invalidate some of the indexes. A drastic though probably reasonable solution in this case would be to simply drop all indexes bound to those replaced/dropped coordinates. For the case where a 1D basic coordinate that corresponds to a dimension is added/replaced, we could automatically generate a new index (see [Section 2.2.4](#224-implicit-indexes)).

We must also ensure that coordinates having a bound index are immutable, e.g., still wrap them into `IndexVariable` objects (even though the `IndexVariable` class might change substantially after this refactoring).

#### 2.4.2 New Dataset/DataArray with updated coordinates

Xarray provides a variety of Dataset/DataArray operations affecting the coordinates and where simply dropping the index(es) is not desirable. For example:

- multi-coordinate indexes could be reduced to single coordinate indexes
  - like in `.reset_index()` or `.sel()` applied on a subset of the levels of a `pandas.MultiIndex` and that internally call `MultiIndex.droplevel` and `MultiIndex.get_loc_level`, respectively
- indexes may be indexed themselves
  - like `pandas.Index` implements `__getitem__()`
  - when indexing their corresponding coordinate(s), e.g., via `.sel()` or `.isel()`, those indexes should be indexed too
  - this might not be supported by all Xarray indexes, though
- some indexes that can't be indexed could still be automatically (re)built in the new Dataset/DataArray
  - like for example building a new `KDTree` index from the selection of a subset of an initial collection of data points
  - this is not always desirable, though, as indexes may be expensive to build
  - a more reasonable option would be to explicitly re-build the index, e.g., using `.set_index()`
- Dataset/DataArray operations involving alignment (see [Section 2.6](#26-using-indexes-for-data-alignment))

### 2.5 Using indexes for data selection

One main use of indexes is label-based data selection using the DataArray/Dataset `.sel()` method. This refactoring would introduce a number of API changes that could go through some depreciation cycles:

- the keys of the mapping given to `indexers` (or the names of `indexer_kwargs`) would not correspond to only dimension names but could be the name of any coordinate that has an index
- for a `pandas.MultiIndex`, if no dimension-coordinate is created by default (see [Section 2.2.4](#224-implicit-coordinates)), providing dict-like objects as indexers should be depreciated
- there should be the possibility to provide additional options to the indexes that support specific selection features (e.g., Scikit-learn's `BallTree`'s `dualtree` query option to boost performance).
  - the best API is not trivial here, since `.sel()` may accept indexers passed to several indexes (which should still be supported for convenience and compatibility), and indexes may have similar options with different semantics
  - we could introduce a new parameter like `index_options: Dict[XarrayIndex, Dict[str, Any]]` to pass options grouped by index
- the `method` and `tolerance` parameters are specific to `pandas.Index` and would not be supported by all indexes: probably best is to eventually pass those arguments as `index_options`
- the list valid indexer types might be extended in order to support new ways of indexing data, e.g., unordered selection of all points within a given range
  - alternatively, we could reuse existing indexer types with different semantics depending on the index, e.g., using `slice(min, max, None)` for unordered range selection

With the new data model proposed here, an ambiguous situation may occur when indexers are given for several coordinates that share the same dimension but not the same index, e.g., from the example in [Section 1](#1-data-model):

```python
da.sel(x=..., y=..., lat=..., lon=...)
```

The easiest solution for this situation would be to raise an error. Alternatively, we could introduce a new parameter to specify how to combine the resulting integer indexers (i.e., union vs intersection), although this could already be achieved by chaining `.sel()` calls or combining `.sel()` with `.merge()` (it may or may not be straightforward).

### 2.6 Using indexes for data alignment

Another main use if indexes is data alignment in various operations. Some considerations regarding alignment and flexible indexes:

- support for alignment should probably be optional for an `XarrayIndex` subclass.
  - like `pandas.Index`, the index wrapper classes that support it should implement `.equals()`, `.union()` and/or `.intersection()`
  - support might be partial if that makes sense (outer, inner, left, right, exact...).
  - index equality might involve more than just the labels: for example a spatial index might be used to check if the coordinate system (CRS) is identical for two sets of coordinates
  - some indexes might implement inexact alignment, like in [#4489](https://github.com/pydata/xarray/pull/4489) or a `KDTree` index that selects nearest-neighbors within a given tolerance
  - alignment may be "multi-dimensional", i.e., the `KDTree` example above vs. dimensions aligned independently of each other
- we need to decide what to do when one dimension has more than one index that supports alignment
  - we should probably raise unless the user explicitly specify which index to use for the alignment
- we need to decide what to do when one dimension has one or more index(es) but none support alignment
  - either we raise or we fail back (silently) to alignment based on dimension size
- for inexact alignment, the tolerance threshold might be given when building the index and/or when performing the alignment
- are there cases where we want a specific index to perform alignment and another index to perform selection?
  - it would be tricky to support that unless we allow multiple indexes per coordinate
  - alternatively, underlying indexes could be picked internally in a "meta" index for one operation or another, although the risk is to eventually have to deal with an explosion of index wrapper classes with different meta indexes for each combination that we'd like to use.

### 2.7 Using indexes for other purposes

Xarray also provides a number of Dataset/DataArray methods where indexes are used in various ways, e.g.,

- `resample` (`CFTimeIndex` and a `DatetimeIntervalIndex`)
- `DatetimeAccessor` & `TimedeltaAccessor` properties (`CFTimeIndex` and a `DatetimeIntervalIndex`)
- `interp` & `interpolate_na`,
   - with `IntervalIndex`, these become regridding operations. Should we support hooks for these operations?
- `differentiate`, `integrate`, `polyfit`
   - raise an error if not a "simple" 1D index?
- `pad`
- `coarsen` has to make choices about output index labels.
- `sortby`
- `stack`/`unstack`
- plotting
    - `plot.pcolormesh` "infers" interval breaks along axes, which are really inferred `bounds` for the appropriate indexes.
    - `plot.step` again uses `bounds`. In fact, we may even want `step` to be the default 1D plotting function if the axis has `bounds` attached.

It would be reasonable to first restrict those methods to the indexes that are currently available in Xarray, and maybe extend the `XarrayIndex` API later upon request when the opportunity arises.

Conversely, nothing should prevent implementing "non-standard" API in 3rd-party `XarrayIndex` subclasses that could be used in DataArray/Dataset extensions (accessors). For example, we might want to reuse a `KDTree` index to compute k-nearest neighbors (returning a DataArray/Dataset with a new dimension) and/or the distances to the nearest neighbors (returning a DataArray/Dataset with a new data variable).

### 2.8 Index encoding

Indexes don't need to be directly serializable since we could (re)build them from their corresponding coordinate(s). However, it would be useful if some indexes could be encoded/decoded to/from a set of arrays that would allow optimized reconstruction and/or storage, e.g.,

- `pandas.MultiIndex` -> `index.levels` and `index.codes`
- Scikit-learn's `KDTree` and `BallTree` that use an array-based representation of an immutable tree structure

## 3. Index representation in DataArray/Dataset's `repr`

Since indexes would become 1st class citizen of Xarray's data model, they deserve their own section in Dataset/DataArray `repr` that could look like:

```
<xarray.DataArray (x: 2, y: 2)>
array([[5.4, 7.8],
       [6.2, 4.7]])
Coordinates:
  * lon      (x, y) float64 10.2 15.2 12.6 17.6
  * lat      (x, y) float64 40.2 45.6 42.2 47.6
  * x        (x) float64 200.0 400.0
  * y        (y) float64 800.0 1e+03
Indexes:
  lat, lon     <SpatialIndex coords=(lat, lon) dims=(x, y)>
  x            <PandasIndexWrapper>
  y            <PandasIndexWrapper>
```

To keep the `repr` compact, we could:

- consolidate entries that map to the same index object, and have an short inline repr for `XarrayIndex` object
- collapse the index section by default in the HTML `repr`
- maybe omit all trivial indexes for 1D coordinates that match the dimension name

## 4. `IndexVariable`

`IndexVariable` is currently used to wrap a `pandas.Index` as a variable, which would not be relevant after this refactoring since it is aimed at decoupling indexes and variables.

We'll probably need to move elsewhere some of the features implemented in `IndexVariable` to:

- ensure that all coordinates with an index are immutable (see [Section 2.4.1](#241-mutable-coordinates))
  - do not set values directly, do not (re)chunk (even though it may be already chunked), do not load, do not convert to sparse/dense, etc.
- directly reuse index's data when that's possible
  - in the case of a `pandas.Index`, it makes little sense to have duplicate data (e.g., as a NumPy array) for its corresponding coordinate
- convert a variable into a `pandas.Index` using `.to_index()` (for backwards compatibility).

Other `IndexVariable` API like `level_names` and `get_level_variable()` would not useful anymore: it is specific to how we currently deal with `pandas.MultiIndex` and virtual "level" coordinates in Xarray.

## 5. Chunked coordinates and/or indexers

We could take opportunity of this refactoring to better leverage chunked coordinates (and/or chunked indexers for data selection). There's two ways to enable it:

A. support for chunked coordinates is left to the index
B. support for chunked coordinates is index agnostic and is implemented in Xarray

As an example for B, [xoak](https://github.com/ESM-VFC/xoak) supports building an index for each chunk, which is coupled with a two-step data selection process (cross-index queries + brute force "reduction" look-up). There is an example [here](https://xoak.readthedocs.io/en/latest/examples/dask_support.html). This may be tedious to generalize this to other kinds of operations, though. Xoak's Dask support is rather experimental, not super stable (it's quite hard to control index replication and data transfer between Dask workers with the default settings), and depends on whether indexes are thread-safe and/or serializable.

Option A may be more reasonable for now.

## 6. Coordinate duck arrays

Another opportunity of this refactoring is support for duck arrays as index coordinates. Decoupling coordinates and indexes would *de-facto* enable it.

However, support for duck arrays in index-based operations such as data selection or alignment would probably require some protocol extension, e.g.,

```python
class MyDuckArray:
    ...

    def _sel_(self, indexer):
        """Prepare the label-based indexer to conform to this coordinate array."""
        ...
        return new_indexer

    ...
```

For example, a `pint` array would implement `_sel_` to perform indexer unit conversion or raise, warn, or just pass the indexer through if it has no units.
xarray: N-D labeled arrays and datasets
=======================================

.. image:: https://github.com/pydata/xarray/workflows/CI/badge.svg?branch=main
   :target: https://github.com/pydata/xarray/actions?query=workflow%3ACI
.. image:: https://codecov.io/gh/pydata/xarray/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/pydata/xarray
.. image:: https://readthedocs.org/projects/xray/badge/?version=latest
   :target: https://xarray.pydata.org/
.. image:: https://img.shields.io/badge/benchmarked%20by-asv-green.svg?style=flat
  :target: https://pandas.pydata.org/speed/xarray/
.. image:: https://img.shields.io/pypi/v/xarray.svg
   :target: https://pypi.python.org/pypi/xarray/
.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/python/black
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.598201.svg
   :target: https://doi.org/10.5281/zenodo.598201


**xarray** (formerly **xray**) is an open source project and Python package
that makes working with labelled multi-dimensional arrays simple,
efficient, and fun!

Xarray introduces labels in the form of dimensions, coordinates and
attributes on top of raw NumPy_-like arrays, which allows for a more
intuitive, more concise, and less error-prone developer experience.
The package includes a large and growing library of domain-agnostic functions
for advanced analytics and visualization with these data structures.

Xarray was inspired by and borrows heavily from pandas_, the popular data
analysis package focused on labelled tabular data.
It is particularly tailored to working with netCDF_ files, which were the
source of xarray's data model, and integrates tightly with dask_ for parallel
computing.

.. _NumPy: https://www.numpy.org
.. _pandas: https://pandas.pydata.org
.. _dask: https://dask.org
.. _netCDF: https://www.unidata.ucar.edu/software/netcdf

Why xarray?
-----------

Multi-dimensional (a.k.a. N-dimensional, ND) arrays (sometimes called
"tensors") are an essential part of computational science.
They are encountered in a wide range of fields, including physics, astronomy,
geoscience, bioinformatics, engineering, finance, and deep learning.
In Python, NumPy_ provides the fundamental data structure and API for
working with raw ND arrays.
However, real-world datasets are usually more than just raw numbers;
they have labels which encode information about how the array values map
to locations in space, time, etc.

Xarray doesn't just keep track of labels on arrays -- it uses them to provide a
powerful and concise interface. For example:

-  Apply operations over dimensions by name: ``x.sum('time')``.
-  Select values by label instead of integer location:
   ``x.loc['2014-01-01']`` or ``x.sel(time='2014-01-01')``.
-  Mathematical operations (e.g., ``x - y``) vectorize across multiple
   dimensions (array broadcasting) based on dimension names, not shape.
-  Flexible split-apply-combine operations with groupby:
   ``x.groupby('time.dayofyear').mean()``.
-  Database like alignment based on coordinate labels that smoothly
   handles missing values: ``x, y = xr.align(x, y, join='outer')``.
-  Keep track of arbitrary metadata in the form of a Python dictionary:
   ``x.attrs``.

Documentation
-------------

Learn more about xarray in its official documentation at https://xarray.pydata.org/

Contributing
------------

You can find information about contributing to xarray at our `Contributing page <https://xarray.pydata.org/en/latest/contributing.html#>`_.

Get in touch
------------

- Ask usage questions ("How do I?") on `StackOverflow`_.
- Report bugs, suggest features or view the source code `on GitHub`_.
- For less well defined questions or ideas, or to announce other projects of
  interest to xarray users, use the `mailing list`_.

.. _StackOverFlow: https://stackoverflow.com/questions/tagged/python-xarray
.. _mailing list: https://groups.google.com/forum/#!forum/xarray
.. _on GitHub: https://github.com/pydata/xarray

NumFOCUS
--------

.. image:: https://numfocus.org/wp-content/uploads/2017/07/NumFocus_LRG.png
   :scale: 25 %
   :target: https://numfocus.org/

Xarray is a fiscally sponsored project of NumFOCUS_, a nonprofit dedicated
to supporting the open source scientific computing community. If you like
Xarray and want to support our mission, please consider making a donation_
to support our efforts.

.. _donation: https://numfocus.salsalabs.org/donate-to-xarray/

History
-------

Xarray is an evolution of an internal tool developed at `The Climate
Corporation`__. It was originally written by Climate Corp researchers Stephan
Hoyer, Alex Kleeman and Eugene Brevdo and was released as open source in
May 2014. The project was renamed from "xray" in January 2016. Xarray became a
fiscally sponsored project of NumFOCUS_ in August 2018.

__ http://climate.com/
.. _NumFOCUS: https://numfocus.org

License
-------

Copyright 2014-2019, xarray Developers

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Xarray bundles portions of pandas, NumPy and Seaborn, all of which are available
under a "3-clause BSD" license:
- pandas: setup.py, xarray/util/print_versions.py
- NumPy: xarray/core/npcompat.py
- Seaborn: _determine_cmap_params in xarray/core/plot/utils.py

Xarray also bundles portions of CPython, which is available under the "Python
Software Foundation License" in xarray/core/pycompat.py.

Xarray uses icons from the icomoon package (free version), which is
available under the "CC BY 4.0" license.

The full text of these licenses are included in the licenses directory.
.. _contributing:

**********************
Contributing to xarray
**********************


.. note::

  Large parts of this document came from the `Pandas Contributing
  Guide <http://pandas.pydata.org/pandas-docs/stable/contributing.html>`_.

Where to start?
===============

All contributions, bug reports, bug fixes, documentation improvements,
enhancements, and ideas are welcome.

If you are brand new to *xarray* or open-source development, we recommend going
through the `GitHub "issues" tab <https://github.com/pydata/xarray/issues>`_
to find issues that interest you. There are a number of issues listed under
`Documentation <https://github.com/pydata/xarray/issues?q=is%3Aissue+is%3Aopen+label%3Adocumentation>`_
and `good first issue
<https://github.com/pydata/xarray/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22>`_
where you could start out. Once you've found an interesting issue, you can
return here to get your development environment setup.

Feel free to ask questions on the `mailing list
<https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!forum/xarray>`_.

.. _contributing.bug_reports:

Bug reports and enhancement requests
====================================

Bug reports are an important part of making *xarray* more stable. Having a complete bug
report will allow others to reproduce the bug and provide insight into fixing. See
`this stackoverflow article <https://stackoverflow.com/help/mcve>`_ for tips on
writing a good bug report.

Trying out the bug-producing code on the *main* branch is often a worthwhile exercise
to confirm that the bug still exists. It is also worth searching existing bug reports and
pull requests to see if the issue has already been reported and/or fixed.

Bug reports must:

#. Include a short, self-contained Python snippet reproducing the problem.
   You can format the code nicely by using `GitHub Flavored Markdown
   <http://github.github.com/github-flavored-markdown/>`_::

      ```python
      import xarray as xr
      ds = xr.Dataset(...)

      ...
      ```

#. Include the full version string of *xarray* and its dependencies. You can use the
   built in function::

      ```python
      import xarray as xr
      xr.show_versions()

      ...
      ```

#. Explain why the current behavior is wrong/not desired and what you expect instead.

The issue will then show up to the *xarray* community and be open to comments/ideas
from others.

.. _contributing.github:

Working with the code
=====================

Now that you have an issue you want to fix, enhancement to add, or documentation
to improve, you need to learn how to work with GitHub and the *xarray* code base.

.. _contributing.version_control:

Version control, Git, and GitHub
--------------------------------

To the new user, working with Git is one of the more daunting aspects of contributing
to *xarray*.  It can very quickly become overwhelming, but sticking to the guidelines
below will help keep the process straightforward and mostly trouble free.  As always,
if you are having difficulties please feel free to ask for help.

The code is hosted on `GitHub <https://www.github.com/pydata/xarray>`_. To
contribute you will need to sign up for a `free GitHub account
<https://github.com/signup/free>`_. We use `Git <http://git-scm.com/>`_ for
version control to allow many people to work together on the project.

Some great resources for learning Git:

* the `GitHub help pages <http://help.github.com/>`_.
* the `NumPy's documentation <http://docs.scipy.org/doc/numpy/dev/index.html>`_.
* Matthew Brett's `Pydagogue <http://matthew-brett.github.io/pydagogue/>`_.

Getting started with Git
------------------------

`GitHub has instructions <http://help.github.com/set-up-git-redirect>`__ for installing git,
setting up your SSH key, and configuring git.  All these steps need to be completed before
you can work seamlessly between your local repository and GitHub.

.. _contributing.forking:

Forking
-------

You will need your own fork to work on the code. Go to the `xarray project
page <https://github.com/pydata/xarray>`_ and hit the ``Fork`` button. You will
want to clone your fork to your machine::

    git clone https://github.com/your-user-name/xarray.git
    cd xarray
    git remote add upstream https://github.com/pydata/xarray.git

This creates the directory `xarray` and connects your repository to
the upstream (main project) *xarray* repository.

.. _contributing.dev_env:

Creating a development environment
----------------------------------

To test out code changes, you'll need to build *xarray* from source, which
requires a Python environment. If you're making documentation changes, you can
skip to :ref:`contributing.documentation` but you won't be able to build the
documentation locally before pushing your changes.

.. _contributiong.dev_python:

Creating a Python Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before starting any development, you'll need to create an isolated xarray
development environment:

- Install either `Anaconda <https://www.anaconda.com/download/>`_ or `miniconda
  <https://conda.io/miniconda.html>`_
- Make sure your conda is up to date (``conda update conda``)
- Make sure that you have :ref:`cloned the repository <contributing.forking>`
- ``cd`` to the *xarray* source directory

We'll now kick off a two-step process:

1. Install the build dependencies
2. Build and install xarray

.. code-block:: sh

   # Create and activate the build environment
   conda create -c conda-forge -n xarray-tests python=3.8

   # This is for Linux and MacOS
   conda env update -f ci/requirements/environment.yml

   # On windows, use environment-windows.yml instead
   conda env update -f ci/requirements/environment-windows.yml

   conda activate xarray-tests

   # or with older versions of Anaconda:
   source activate xarray-tests

   # Build and install xarray
   pip install -e .

At this point you should be able to import *xarray* from your locally
built version:

.. code-block:: sh

   $ python  # start an interpreter
   >>> import xarray
   >>> xarray.__version__
   '0.10.0+dev46.g015daca'

This will create the new environment, and not touch any of your existing environments,
nor any existing Python installation.

To view your environments::

      conda info -e

To return to your root environment::

      conda deactivate

See the full conda docs `here <http://conda.pydata.org/docs>`__.

Creating a branch
-----------------

You want your ``main`` branch to reflect only production-ready code, so create a
feature branch before making your changes. For example::

    git branch shiny-new-feature
    git checkout shiny-new-feature

The above can be simplified to::

    git checkout -b shiny-new-feature

This changes your working directory to the shiny-new-feature branch.  Keep any
changes in this branch specific to one bug or feature so it is clear
what the branch brings to *xarray*. You can have many "shiny-new-features"
and switch in between them using the ``git checkout`` command.

To update this branch, you need to retrieve the changes from the ``main`` branch::

    git fetch upstream
    git merge upstream/main

This will combine your commits with the latest *xarray* git ``main``.  If this
leads to merge conflicts, you must resolve these before submitting your pull
request.  If you have uncommitted changes, you will need to ``git stash`` them
prior to updating.  This will effectively store your changes, which can be
reapplied after updating.

.. _contributing.documentation:

Contributing to the documentation
=================================

If you're not the developer type, contributing to the documentation is still of
huge value. You don't even have to be an expert on *xarray* to do so! In fact,
there are sections of the docs that are worse off after being written by
experts. If something in the docs doesn't make sense to you, updating the
relevant section after you figure it out is a great way to ensure it will help
the next person.

.. contents:: Documentation:
   :local:


About the *xarray* documentation
--------------------------------

The documentation is written in **reStructuredText**, which is almost like writing
in plain English, and built using `Sphinx <http://sphinx-doc.org/>`__. The
Sphinx Documentation has an excellent `introduction to reST
<http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`__. Review the Sphinx docs to perform more
complex changes to the documentation as well.

Some other important things to know about the docs:

- The *xarray* documentation consists of two parts: the docstrings in the code
  itself and the docs in this folder ``xarray/doc/``.

  The docstrings are meant to provide a clear explanation of the usage of the
  individual functions, while the documentation in this folder consists of
  tutorial-like overviews per topic together with some other information
  (what's new, installation, etc).

- The docstrings follow the **NumPy Docstring Standard**, which is used widely
  in the Scientific Python community. This standard specifies the format of
  the different sections of the docstring. See `this document
  <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`_
  for a detailed explanation, or look at some of the existing functions to
  extend it in a similar manner.

- The tutorials make heavy use of the `ipython directive
  <http://matplotlib.org/sampledoc/ipython_directive.html>`_ sphinx extension.
  This directive lets you put code in the documentation which will be run
  during the doc build. For example:

  .. code:: rst

      .. ipython:: python

          x = 2
          x ** 3

  will be rendered as::

      In [1]: x = 2

      In [2]: x ** 3
      Out[2]: 8

  Almost all code examples in the docs are run (and the output saved) during the
  doc build. This approach means that code examples will always be up to date,
  but it does make building the docs a bit more complex.

- Our API documentation in ``doc/api.rst`` houses the auto-generated
  documentation from the docstrings. For classes, there are a few subtleties
  around controlling which methods and attributes have pages auto-generated.

  Every method should be included in a ``toctree`` in ``api.rst``, else Sphinx
  will emit a warning.


How to build the *xarray* documentation
---------------------------------------

Requirements
~~~~~~~~~~~~
Make sure to follow the instructions on :ref:`creating a development environment above <contributing.dev_env>`, but
to build the docs you need to use the environment file ``ci/requirements/doc.yml``.

.. code-block:: sh

    # Create and activate the docs environment
    conda env create -f ci/requirements/doc.yml
    conda activate xarray-docs

    # or with older versions of Anaconda:
    source activate xarray-docs

    # Build and install xarray
    pip install -e .

Building the documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~

To build the documentation run::

    cd doc/
    make html

Then you can find the HTML output in the folder ``xarray/doc/_build/html/``.

The first time you build the docs, it will take quite a while because it has to run
all the code examples and build all the generated docstring pages. In subsequent
evocations, Sphinx will try to only build the pages that have been modified.

If you want to do a full clean build, do::

    make clean
    make html

.. _contributing.code:

Contributing to the code base
=============================

.. contents:: Code Base:
   :local:

Code standards
--------------

Writing good code is not just about what you write. It is also about *how* you
write it. During :ref:`Continuous Integration <contributing.ci>` testing, several
tools will be run to check your code for stylistic errors.
Generating any warnings will cause the test to fail.
Thus, good style is a requirement for submitting code to *xarray*.

In addition, because a lot of people use our library, it is important that we
do not make sudden changes to the code that could have the potential to break
a lot of user code as a result, that is, we need it to be as *backwards compatible*
as possible to avoid mass breakages.

Code Formatting
~~~~~~~~~~~~~~~

xarray uses several tools to ensure a consistent code format throughout the project:

- `Black <https://black.readthedocs.io/en/stable/>`_ for standardized
  code formatting
- `blackdoc <https://blackdoc.readthedocs.io/en/stable/>`_ for
  standardized code formatting in documentation
- `Flake8 <http://flake8.pycqa.org/en/latest/>`_ for general code quality
- `isort <https://github.com/timothycrosley/isort>`_ for standardized order in imports.
  See also `flake8-isort <https://github.com/gforcada/flake8-isort>`_.
- `mypy <http://mypy-lang.org/>`_ for static type checking on `type hints
  <https://docs.python.org/3/library/typing.html>`_

We highly recommend that you setup `pre-commit hooks <https://pre-commit.com/>`_
to automatically run all the above tools every time you make a git commit. This
can be done by running::

   pre-commit install

from the root of the xarray repository. You can skip the pre-commit checks
with ``git commit --no-verify``.


Backwards Compatibility
~~~~~~~~~~~~~~~~~~~~~~~

Please try to maintain backwards compatibility. *xarray* has a growing number of users with
lots of existing code, so don't break it if at all possible.  If you think breakage is
required, clearly state why as part of the pull request.

Be especially careful when changing function and method signatures, because any change
may require a deprecation warning. For example, if your pull request means that the
argument ``old_arg`` to ``func`` is no longer valid, instead of simply raising an error if
a user passes ``old_arg``, we would instead catch it:

.. code-block:: python

    def func(new_arg, old_arg=None):
        if old_arg is not None:
            from warnings import warn

            warn(
                "`old_arg` has been deprecated, and in the future will raise an error."
                "Please use `new_arg` from now on.",
                DeprecationWarning,
            )

            # Still do what the user intended here

This temporary check would then be removed in a subsequent version of xarray.
This process of first warning users before actually breaking their code is known as a
"deprecation cycle", and makes changes significantly easier to handle both for users
of xarray, and for developers of other libraries that depend on xarray.


.. _contributing.ci:

Testing With Continuous Integration
-----------------------------------

The *xarray* test suite runs automatically the
`GitHub Actions <https://docs.github.com/en/free-pro-team@latest/actions>`__,
continuous integration service, once your pull request is submitted.

A pull-request will be considered for merging when you have an all 'green' build. If any
tests are failing, then you will get a red 'X', where you can click through to see the
individual failed tests. This is an example of a green build.

.. image:: _static/ci.png

.. note::

   Each time you push to your PR branch, a new run of the tests will be
   triggered on the CI. If they haven't already finished, tests for any older
   commits on the same branch will be automatically cancelled.

.. _contributing.tdd:


Test-driven development/code writing
------------------------------------

*xarray* is serious about testing and strongly encourages contributors to embrace
`test-driven development (TDD) <http://en.wikipedia.org/wiki/Test-driven_development>`_.
This development process "relies on the repetition of a very short development cycle:
first the developer writes an (initially failing) automated test case that defines a desired
improvement or new function, then produces the minimum amount of code to pass that test."
So, before actually writing any code, you should write your tests.  Often the test can be
taken from the original GitHub issue.  However, it is always worth considering additional
use cases and writing corresponding tests.

Adding tests is one of the most common requests after code is pushed to *xarray*.  Therefore,
it is worth getting in the habit of writing tests ahead of time so that this is never an issue.

Like many packages, *xarray* uses `pytest
<http://doc.pytest.org/en/latest/>`_ and the convenient
extensions in `numpy.testing
<http://docs.scipy.org/doc/numpy/reference/routines.testing.html>`_.

Writing tests
~~~~~~~~~~~~~

All tests should go into the ``tests`` subdirectory of the specific package.
This folder contains many current examples of tests, and we suggest looking to these for
inspiration.  If your test requires working with files or
network connectivity, there is more information on the `testing page
<https://github.com/pydata/xarray/wiki/Testing>`_ of the wiki.

The ``xarray.testing`` module has many special ``assert`` functions that
make it easier to make statements about whether DataArray or Dataset objects are
equivalent. The easiest way to verify that your code is correct is to
explicitly construct the result you expect, then compare the actual result to
the expected correct result::

    def test_constructor_from_0d():
        expected = Dataset({None: ([], 0)})[None]
        actual = DataArray(0)
        assert_identical(expected, actual)

Transitioning to ``pytest``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

*xarray* existing test structure is *mostly* classed based, meaning that you will
typically find tests wrapped in a class.

.. code-block:: python

    class TestReallyCoolFeature:
        ...

Going forward, we are moving to a more *functional* style using the
`pytest <http://doc.pytest.org/en/latest/>`__ framework, which offers a richer
testing framework that will facilitate testing and developing. Thus, instead of
writing test classes, we will write test functions like this:

.. code-block:: python

    def test_really_cool_feature():
        ...

Using ``pytest``
~~~~~~~~~~~~~~~~

Here is an example of a self-contained set of tests that illustrate multiple
features that we like to use.

- functional style: tests are like ``test_*`` and *only* take arguments that are either
  fixtures or parameters
- ``pytest.mark`` can be used to set metadata on test functions, e.g. ``skip`` or ``xfail``.
- using ``parametrize``: allow testing of multiple cases
- to set a mark on a parameter, ``pytest.param(..., marks=...)`` syntax should be used
- ``fixture``, code for object construction, on a per-test basis
- using bare ``assert`` for scalars and truth-testing
- ``assert_equal`` and ``assert_identical`` from the ``xarray.testing`` module for xarray object comparisons.
- the typical pattern of constructing an ``expected`` and comparing versus the ``result``

We would name this file ``test_cool_feature.py`` and put in an appropriate place in the
``xarray/tests/`` structure.

.. TODO: confirm that this actually works

.. code-block:: python

    import pytest
    import numpy as np
    import xarray as xr
    from xarray.testing import assert_equal


    @pytest.mark.parametrize("dtype", ["int8", "int16", "int32", "int64"])
    def test_dtypes(dtype):
        assert str(np.dtype(dtype)) == dtype


    @pytest.mark.parametrize(
        "dtype",
        [
            "float32",
            pytest.param("int16", marks=pytest.mark.skip),
            pytest.param(
                "int32", marks=pytest.mark.xfail(reason="to show how it works")
            ),
        ],
    )
    def test_mark(dtype):
        assert str(np.dtype(dtype)) == "float32"


    @pytest.fixture
    def dataarray():
        return xr.DataArray([1, 2, 3])


    @pytest.fixture(params=["int8", "int16", "int32", "int64"])
    def dtype(request):
        return request.param


    def test_series(dataarray, dtype):
        result = dataarray.astype(dtype)
        assert result.dtype == dtype

        expected = xr.DataArray(np.array([1, 2, 3], dtype=dtype))
        assert_equal(result, expected)



A test run of this yields

.. code-block:: shell

   ((xarray) $ pytest test_cool_feature.py -v
    =============================== test session starts ================================
    platform darwin -- Python 3.6.4, pytest-3.2.1, py-1.4.34, pluggy-0.4.0 --
    cachedir: ../../.cache
    plugins: cov-2.5.1, hypothesis-3.23.0
    collected 11 items

    test_cool_feature.py::test_dtypes[int8] PASSED
    test_cool_feature.py::test_dtypes[int16] PASSED
    test_cool_feature.py::test_dtypes[int32] PASSED
    test_cool_feature.py::test_dtypes[int64] PASSED
    test_cool_feature.py::test_mark[float32] PASSED
    test_cool_feature.py::test_mark[int16] SKIPPED
    test_cool_feature.py::test_mark[int32] xfail
    test_cool_feature.py::test_series[int8] PASSED
    test_cool_feature.py::test_series[int16] PASSED
    test_cool_feature.py::test_series[int32] PASSED
    test_cool_feature.py::test_series[int64] PASSED

    ================== 9 passed, 1 skipped, 1 xfailed in 1.83 seconds ==================

Tests that we have ``parametrized`` are now accessible via the test name, for
example we could run these with ``-k int8`` to sub-select *only* those tests
which match ``int8``.


.. code-block:: shell

   ((xarray) bash-3.2$ pytest  test_cool_feature.py  -v -k int8
   =========================== test session starts ===========================
   platform darwin -- Python 3.6.2, pytest-3.2.1, py-1.4.31, pluggy-0.4.0
   collected 11 items

   test_cool_feature.py::test_dtypes[int8] PASSED
   test_cool_feature.py::test_series[int8] PASSED


Running the test suite
----------------------

The tests can then be run directly inside your Git clone (without having to
install *xarray*) by typing::

    pytest xarray

The tests suite is exhaustive and takes a few minutes.  Often it is
worth running only a subset of tests first around your changes before running the
entire suite.

The easiest way to do this is with::

    pytest xarray/path/to/test.py -k regex_matching_test_name

Or with one of the following constructs::

    pytest xarray/tests/[test-module].py
    pytest xarray/tests/[test-module].py::[TestClass]
    pytest xarray/tests/[test-module].py::[TestClass]::[test_method]

Using `pytest-xdist <https://pypi.python.org/pypi/pytest-xdist>`_, one can
speed up local testing on multicore machines, by running pytest with the optional -n argument::

    pytest xarray -n 4

This can significantly reduce the time it takes to locally run tests before
submitting a pull request.

For more, see the `pytest <http://doc.pytest.org/en/latest/>`_ documentation.

Running the performance test suite
----------------------------------

Performance matters and it is worth considering whether your code has introduced
performance regressions.  *xarray* is starting to write a suite of benchmarking tests
using `asv <https://github.com/spacetelescope/asv>`__
to enable easy monitoring of the performance of critical *xarray* operations.
These benchmarks are all found in the ``xarray/asv_bench`` directory.  asv
supports both python2 and python3.

To use all features of asv, you will need either ``conda`` or
``virtualenv``. For more details please check the `asv installation
webpage <https://asv.readthedocs.io/en/latest/installing.html>`_.

To install asv::

    pip install git+https://github.com/spacetelescope/asv

If you need to run a benchmark, change your directory to ``asv_bench/`` and run::

    asv continuous -f 1.1 upstream/main HEAD

You can replace ``HEAD`` with the name of the branch you are working on,
and report benchmarks that changed by more than 10%.
The command uses ``conda`` by default for creating the benchmark
environments. If you want to use virtualenv instead, write::

    asv continuous -f 1.1 -E virtualenv upstream/main HEAD

The ``-E virtualenv`` option should be added to all ``asv`` commands
that run benchmarks. The default value is defined in ``asv.conf.json``.

Running the full benchmark suite can take up to one hour and use up a few GBs of RAM.
Usually it is sufficient to paste only a subset of the results into the pull
request to show that the committed changes do not cause unexpected performance
regressions.  You can run specific benchmarks using the ``-b`` flag, which
takes a regular expression.  For example, this will only run tests from a
``xarray/asv_bench/benchmarks/groupby.py`` file::

    asv continuous -f 1.1 upstream/main HEAD -b ^groupby

If you want to only run a specific group of tests from a file, you can do it
using ``.`` as a separator. For example::

    asv continuous -f 1.1 upstream/main HEAD -b groupby.GroupByMethods

will only run the ``GroupByMethods`` benchmark defined in ``groupby.py``.

You can also run the benchmark suite using the version of *xarray*
already installed in your current Python environment. This can be
useful if you do not have ``virtualenv`` or ``conda``, or are using the
``setup.py develop`` approach discussed above; for the in-place build
you need to set ``PYTHONPATH``, e.g.
``PYTHONPATH="$PWD/.." asv [remaining arguments]``.
You can run benchmarks using an existing Python
environment by::

    asv run -e -E existing

or, to use a specific Python interpreter,::

    asv run -e -E existing:python3.6

This will display stderr from the benchmarks, and use your local
``python`` that comes from your ``$PATH``.

Information on how to write a benchmark and how to use asv can be found in the
`asv documentation <https://asv.readthedocs.io/en/latest/writing_benchmarks.html>`_.

..
   TODO: uncomment once we have a working setup
         see https://github.com/pydata/xarray/pull/5066

   The *xarray* benchmarking suite is run remotely and the results are
   available `here <http://pandas.pydata.org/speed/xarray/>`_.

Documenting your code
---------------------

Changes should be reflected in the release notes located in ``doc/whats-new.rst``.
This file contains an ongoing change log for each release.  Add an entry to this file to
document your fix, enhancement or (unavoidable) breaking change.  Make sure to include the
GitHub issue number when adding your entry (using ``:issue:`1234```, where ``1234`` is the
issue/pull request number).

If your code is an enhancement, it is most likely necessary to add usage
examples to the existing documentation.  This can be done following the section
regarding documentation :ref:`above <contributing.documentation>`.

Contributing your changes to *xarray*
=====================================

Committing your code
--------------------

Keep style fixes to a separate commit to make your pull request more readable.

Once you've made changes, you can see them by typing::

    git status

If you have created a new file, it is not being tracked by git. Add it by typing::

    git add path/to/file-to-be-added.py

Doing 'git status' again should give something like::

    # On branch shiny-new-feature
    #
    #       modified:   /relative/path/to/file-you-added.py
    #

The following defines how a commit message should be structured:

    * A subject line with `< 72` chars.
    * One blank line.
    * Optionally, a commit message body.

Please reference the relevant GitHub issues in your commit message using ``GH1234`` or
``#1234``.  Either style is fine, but the former is generally preferred.

Now you can commit your changes in your local repository::

    git commit -m

Pushing your changes
--------------------

When you want your changes to appear publicly on your GitHub page, push your
forked feature branch's commits::

    git push origin shiny-new-feature

Here ``origin`` is the default name given to your remote repository on GitHub.
You can see the remote repositories::

    git remote -v

If you added the upstream repository as described above you will see something
like::

    origin  git@github.com:yourname/xarray.git (fetch)
    origin  git@github.com:yourname/xarray.git (push)
    upstream        git://github.com/pydata/xarray.git (fetch)
    upstream        git://github.com/pydata/xarray.git (push)

Now your code is on GitHub, but it is not yet a part of the *xarray* project.  For that to
happen, a pull request needs to be submitted on GitHub.

Review your code
----------------

When you're ready to ask for a code review, file a pull request. Before you do, once
again make sure that you have followed all the guidelines outlined in this document
regarding code style, tests, performance tests, and documentation. You should also
double check your branch changes against the branch it was based on:

#. Navigate to your repository on GitHub -- https://github.com/your-user-name/xarray
#. Click on ``Branches``
#. Click on the ``Compare`` button for your feature branch
#. Select the ``base`` and ``compare`` branches, if necessary. This will be ``main`` and
   ``shiny-new-feature``, respectively.

Finally, make the pull request
------------------------------

If everything looks good, you are ready to make a pull request.  A pull request is how
code from a local repository becomes available to the GitHub community and can be looked
at and eventually merged into the ``main`` version.  This pull request and its associated
changes will eventually be committed to the ``main`` branch and available in the next
release.  To submit a pull request:

#. Navigate to your repository on GitHub
#. Click on the ``Pull Request`` button
#. You can then click on ``Commits`` and ``Files Changed`` to make sure everything looks
   okay one last time
#. Write a description of your changes in the ``Preview Discussion`` tab
#. Click ``Send Pull Request``.

This request then goes to the repository maintainers, and they will review
the code. If you need to make more changes, you can make them in
your branch, add them to a new commit, push them to GitHub, and the pull request
will automatically be updated.  Pushing them to GitHub again is done by::

    git push origin shiny-new-feature

This will automatically update your pull request with the latest code and restart the
:ref:`Continuous Integration <contributing.ci>` tests.


Delete your merged branch (optional)
------------------------------------

Once your feature branch is accepted into upstream, you'll probably want to get rid of
the branch. First, update your ``main`` branch to check that the merge was successful::

    git fetch upstream
    git checkout main
    git merge upstream/main

Then you can do::

    git branch -D shiny-new-feature

You need to use a upper-case ``-D`` because the branch was squashed into a
single commit before merging. Be careful with this because ``git`` won't warn
you if you accidentally delete an unmerged branch.

If you didn't delete your branch using GitHub's interface, then it will still exist on
GitHub. To delete it there do::

    git push origin --delete shiny-new-feature


PR checklist
------------

- **Properly comment and document your code.** See `"Documenting your code" <https://xarray.pydata.org/en/stable/contributing.html#documenting-your-code>`_.
- **Test that the documentation builds correctly** by typing ``make html`` in the ``doc`` directory. This is not strictly necessary, but this may be easier than waiting for CI to catch a mistake. See `"Contributing to the documentation" <https://xarray.pydata.org/en/stable/contributing.html#contributing-to-the-documentation>`_.
- **Test your code**.

    - Write new tests if needed. See `"Test-driven development/code writing" <https://xarray.pydata.org/en/stable/contributing.html#test-driven-development-code-writing>`_.
    - Test the code using `Pytest <http://doc.pytest.org/en/latest/>`_. Running all tests (type ``pytest`` in the root directory) takes a while, so feel free to only run the tests you think are needed based on your PR (example: ``pytest xarray/tests/test_dataarray.py``). CI will catch any failing tests.
    - By default, the upstream dev CI is disabled on pull request and push events. You can override this behavior per commit by adding a <tt>[test-upstream]</tt> tag to the first line of the commit message. For documentation-only commits, you can skip the CI per commit by adding a "[skip-ci]" tag to the first line of the commit message.

- **Properly format your code** and verify that it passes the formatting guidelines set by `Black <https://black.readthedocs.io/en/stable/>`_ and `Flake8 <http://flake8.pycqa.org/en/latest/>`_. See `"Code formatting" <https://xarray.pydata.org/en/stablcontributing.html#code-formatting>`_. You can use `pre-commit <https://pre-commit.com/>`_ to run these automatically on each commit.

    - Run ``pre-commit run --all-files`` in the root directory. This may modify some files. Confirm and commit any formatting changes.

- **Push your code and** `create a PR on GitHub <https://help.github.com/en/articles/creating-a-pull-request>`_.
- **Use a helpful title for your pull request** by summarizing the main contributions rather than using the latest commit message. If the PR addresses an `issue <https://github.com/pydata/xarray/issues>`_, please `reference it <https://help.github.com/en/articles/autolinked-references-and-urls>`_.
.. currentmodule:: xarray

What's New
==========

.. ipython:: python
    :suppress:

    import numpy as np
    import pandas as pd
    import xarray as xray
    import xarray
    import xarray as xr

    np.random.seed(123456)

.. _whats-new.0.21.1:

v0.21.1 (unreleased)
--------------------

New Features
~~~~~~~~~~~~


Breaking changes
~~~~~~~~~~~~~~~~


Deprecations
~~~~~~~~~~~~


Bug fixes
~~~~~~~~~


Documentation
~~~~~~~~~~~~~


Internal Changes
~~~~~~~~~~~~~~~~


.. _whats-new.0.21.0:

v0.21.0 (27 January 2022)
-------------------------

Many thanks to the 20 contributors to the v0.21.0 release!

Abel Aoun, Anderson Banihirwe, Ant Gib, Chris Roat, Cindy Chiao,
Deepak Cherian, Dominik Stańczak, Fabian Hofmann, Illviljan, Jody Klymak, Joseph
K Aicher, Mark Harfouche, Mathias Hauser, Matthew Roeschke, Maximilian Roos,
Michael Delgado, Pascal Bourgault, Pierre, Ray Bell, Romain Caneill, Tim Heap,
Tom Nicholas, Zeb Nicholls, joseph nowak, keewis.


New Features
~~~~~~~~~~~~
- New top-level function :py:func:`cross`. (:issue:`3279`, :pull:`5365`).
  By `Jimmy Westling <https://github.com/illviljan>`_.
- ``keep_attrs`` support for :py:func:`where` (:issue:`4141`, :issue:`4682`, :pull:`4687`).
  By `Justus Magin <https://github.com/keewis>`_.
- Enable the limit option for dask array in the following methods :py:meth:`DataArray.ffill`, :py:meth:`DataArray.bfill`, :py:meth:`Dataset.ffill` and :py:meth:`Dataset.bfill` (:issue:`6112`)
  By `Joseph Nowak <https://github.com/josephnowak>`_.

Breaking changes
~~~~~~~~~~~~~~~~
- Rely on matplotlib's default datetime converters instead of pandas' (:issue:`6102`, :pull:`6109`).
  By `Jimmy Westling <https://github.com/illviljan>`_.
- Improve repr readability when there are a large number of dimensions in datasets or dataarrays by
  wrapping the text once the maximum display width has been exceeded. (:issue:`5546`, :pull:`5662`)
  By `Jimmy Westling <https://github.com/illviljan>`_.


Deprecations
~~~~~~~~~~~~
- Removed the lock kwarg from the zarr and pydap backends, completing the deprecation cycle started in :issue:`5256`.
  By `Tom Nicholas <https://github.com/TomNicholas>`_.
- Support for ``python 3.7`` has been dropped. (:pull:`5892`)
  By `Jimmy Westling <https://github.com/illviljan>`_.


Bug fixes
~~~~~~~~~
- Preserve chunks when creating a :py:class:`DataArray` from another :py:class:`DataArray`
  (:pull:`5984`). By `Fabian Hofmann <https://github.com/FabianHofmann>`_.
- Properly support :py:meth:`DataArray.ffill`, :py:meth:`DataArray.bfill`, :py:meth:`Dataset.ffill` and :py:meth:`Dataset.bfill` along chunked dimensions (:issue:`6112`).
  By `Joseph Nowak <https://github.com/josephnowak>`_.

- Subclasses of ``byte`` and ``str`` (e.g. ``np.str_`` and ``np.bytes_``) will now serialise to disk rather than raising a ``ValueError: unsupported dtype for netCDF4 variable: object`` as they did previously (:pull:`5264`).
  By `Zeb Nicholls <https://github.com/znicholls>`_.

- Fix applying function with non-xarray arguments using :py:func:`xr.map_blocks`.
  By `Cindy Chiao <https://github.com/tcchiao>`_.

- No longer raise an error for an all-nan-but-one argument to
  :py:meth:`DataArray.interpolate_na` when using `method='nearest'` (:issue:`5994`, :pull:`6144`).
  By `Michael Delgado <https://github.com/delgadom>`_.
- `dt.season <https://xarray.pydata.org/en/stable/generated/xarray.DataArray.dt.season.html>`_  can now handle NaN and NaT.  (:pull:`5876`).
  By `Pierre Loicq <https://github.com/pierreloicq>`_.
- Determination of zarr chunks handles empty lists for encoding chunks or variable chunks that occurs in certain cirumstances (:pull:`5526`). By `Chris Roat <https://github.com/chrisroat>`_.

Internal Changes
~~~~~~~~~~~~~~~~

- Replace ``distutils.version`` with ``packaging.version``  (:issue:`6092`).
  By `Mathias Hauser <https://github.com/mathause>`_.

- Removed internal checks for ``pd.Panel`` (:issue:`6145`).
  By `Matthew Roeschke <https://github.com/mroeschke>`_.

- Add ``pyupgrade`` pre-commit hook (:pull:`6152`).
  By `Maximilian Roos <https://github.com/max-sixty>`_.

.. _whats-new.0.20.2:

v0.20.2 (9 December 2021)
-------------------------

This is a bugfix release to resolve (:issue:`3391`, :issue:`5715`). It also
includes performance improvements in unstacking to a ``sparse`` array and a
number of documentation improvements.

Many thanks to the 20 contributors:

Aaron Spring, Alexandre Poux, Deepak Cherian, Enrico Minack, Fabien Maussion,
Giacomo Caria, Gijom, Guillaume Maze, Illviljan, Joe Hamman, Joseph Hardin, Kai
Mühlbauer, Matt Henderson, Maximilian Roos, Michael Delgado, Robert Gieseke,
Sebastian Weigand and Stephan Hoyer.


Breaking changes
~~~~~~~~~~~~~~~~
- Use complex nan when interpolating complex values out of bounds by default (instead of real nan) (:pull:`6019`).
  By `Alexandre Poux <https://github.com/pums974>`_.

Performance
~~~~~~~~~~~

- Significantly faster unstacking to a ``sparse`` array. :pull:`5577`
  By `Deepak Cherian <https://github.com/dcherian>`_.

Bug fixes
~~~~~~~~~
- :py:func:`xr.map_blocks` and :py:func:`xr.corr` now work when dask is not installed (:issue:`3391`, :issue:`5715`, :pull:`5731`).
  By `Gijom <https://github.com/Gijom>`_.
- Fix plot.line crash for data of shape ``(1, N)`` in _title_for_slice on format_item (:pull:`5948`).
  By `Sebastian Weigand <https://github.com/s-weigand>`_.
- Fix a regression in the removal of duplicate backend entrypoints (:issue:`5944`, :pull:`5959`)
  By `Kai Mühlbauer <https://github.com/kmuehlbauer>`_.
- Fix an issue that datasets from being saved when time variables with units that ``cftime`` can parse but pandas can not were present (:pull:`6049`).
  By `Tim Heap <https://github.com/mx-moth>`_.

Documentation
~~~~~~~~~~~~~

- Better examples in docstrings for groupby and resampling reductions (:pull:`5871`).
  By `Deepak Cherian <https://github.com/dcherian>`_,
  `Maximilian Roos <https://github.com/max-sixty>`_,
  `Jimmy Westling <https://github.com/illviljan>`_ .
- Add list-like possibility for tolerance parameter in the reindex functions.
  By `Antoine Gibek <https://github.com/antscloud>`_,

Internal Changes
~~~~~~~~~~~~~~~~

- Use ``importlib`` to replace functionality of ``pkg_resources`` in
  backend plugins tests. (:pull:`5959`).
  By `Kai Mühlbauer <https://github.com/kmuehlbauer>`_.


.. _whats-new.0.20.1:

v0.20.1 (5 November 2021)
-------------------------

This is a bugfix release to fix :issue:`5930`.

Bug fixes
~~~~~~~~~
- Fix a regression in the detection of the backend entrypoints (:issue:`5930`, :pull:`5931`)
  By `Justus Magin <https://github.com/keewis>`_.

Documentation
~~~~~~~~~~~~~

- Significant improvements to  :ref:`api`. By `Deepak Cherian <https://github.com/dcherian>`_.

.. _whats-new.0.20.0:

v0.20.0 (1 November 2021)
-------------------------

This release brings improved support for pint arrays, methods for weighted standard deviation, variance,
and sum of squares, the option to disable the use of the bottleneck library, significantly improved performance of
unstack, as well as many bugfixes and internal changes.

Many thanks to the 40 contributors to this release!:

Aaron Spring, Akio Taniguchi, Alan D. Snow, arfy slowy, Benoit Bovy, Christian Jauvin, crusaderky, Deepak Cherian,
Giacomo Caria, Illviljan, James Bourbeau, Joe Hamman, Joseph K Aicher, Julien Herzen, Kai Mühlbauer,
keewis, lusewell, Martin K. Scherer, Mathias Hauser, Max Grover, Maxime Liquet, Maximilian Roos, Mike Taves, Nathan Lis,
pmav99, Pushkar Kopparla, Ray Bell, Rio McMahon, Scott Staniewicz, Spencer Clark, Stefan Bender, Taher Chegini,
Thomas Nicholas, Tomas Chor, Tom Augspurger, Victor Negîrneac, Zachary Blackwood, Zachary Moon, and Zeb Nicholls.

New Features
~~~~~~~~~~~~
- Add ``std``, ``var``,  ``sum_of_squares`` to :py:class:`~core.weighted.DatasetWeighted` and :py:class:`~core.weighted.DataArrayWeighted`.
  By `Christian Jauvin <https://github.com/cjauvin>`_.
- Added a :py:func:`get_options` method to xarray's root namespace (:issue:`5698`, :pull:`5716`)
  By `Pushkar Kopparla <https://github.com/pkopparla>`_.
- Xarray now does a better job rendering variable names that are long LaTeX sequences when plotting (:issue:`5681`, :pull:`5682`).
  By `Tomas Chor <https://github.com/tomchor>`_.
- Add an option (``"use_bottleneck"``) to disable the use of ``bottleneck`` using :py:func:`set_options` (:pull:`5560`)
  By `Justus Magin <https://github.com/keewis>`_.
- Added ``**kwargs`` argument to :py:meth:`open_rasterio` to access overviews (:issue:`3269`).
  By `Pushkar Kopparla <https://github.com/pkopparla>`_.
- Added ``storage_options`` argument to :py:meth:`to_zarr` (:issue:`5601`, :pull:`5615`).
  By `Ray Bell <https://github.com/raybellwaves>`_, `Zachary Blackwood <https://github.com/blackary>`_ and
  `Nathan Lis <https://github.com/wxman22>`_.
- Added calendar utilities :py:func:`DataArray.convert_calendar`, :py:func:`DataArray.interp_calendar`, :py:func:`date_range`, :py:func:`date_range_like` and :py:attr:`DataArray.dt.calendar` (:issue:`5155`, :pull:`5233`).
  By `Pascal Bourgault <https://github.com/aulemahal>`_.
- Histogram plots are set with a title displaying the scalar coords if any, similarly to the other plots (:issue:`5791`, :pull:`5792`).
  By `Maxime Liquet <https://github.com/maximlt>`_.
- Slice plots display the coords units in the same way as x/y/colorbar labels (:pull:`5847`).
  By `Victor Negîrneac <https://github.com/caenrigen>`_.
- Added a new :py:attr:`Dataset.chunksizes`, :py:attr:`DataArray.chunksizes`, and :py:attr:`Variable.chunksizes`
  property, which will always return a mapping from dimension names to chunking pattern along that dimension,
  regardless of whether the object is a Dataset, DataArray, or Variable. (:issue:`5846`, :pull:`5900`)
  By `Tom Nicholas <https://github.com/TomNicholas>`_.

Breaking changes
~~~~~~~~~~~~~~~~
- The minimum versions of some dependencies were changed:

  =============== ====== ====
  Package         Old    New
  =============== ====== ====
  cftime          1.1    1.2
  dask            2.15   2.30
  distributed     2.15   2.30
  lxml            4.5    4.6
  matplotlib-base 3.2    3.3
  numba           0.49   0.51
  numpy           1.17   1.18
  pandas          1.0    1.1
  pint            0.15   0.16
  scipy           1.4    1.5
  seaborn         0.10   0.11
  sparse          0.8    0.11
  toolz           0.10   0.11
  zarr            2.4    2.5
  =============== ====== ====

- The ``__repr__`` of a :py:class:`xarray.Dataset`'s ``coords`` and ``data_vars``
  ignore ``xarray.set_option(display_max_rows=...)`` and show the full output
  when called directly as, e.g., ``ds.data_vars`` or ``print(ds.data_vars)``
  (:issue:`5545`, :pull:`5580`).
  By `Stefan Bender <https://github.com/st-bender>`_.

Deprecations
~~~~~~~~~~~~

- Deprecate :py:func:`open_rasterio` (:issue:`4697`, :pull:`5808`).
  By `Alan Snow <https://github.com/snowman2>`_.
- Set the default argument for `roll_coords` to `False` for :py:meth:`DataArray.roll`
  and :py:meth:`Dataset.roll`. (:pull:`5653`)
  By `Tom Nicholas <https://github.com/TomNicholas>`_.
- :py:meth:`xarray.open_mfdataset` will now error instead of warn when a value for ``concat_dim`` is
  passed alongside ``combine='by_coords'``.
  By `Tom Nicholas <https://github.com/TomNicholas>`_.

Bug fixes
~~~~~~~~~

- Fix ZeroDivisionError from saving dask array with empty dimension (:issue: `5741`).
  By `Joseph K Aicher <https://github.com/jaicher>`_.
- Fixed performance bug where ``cftime`` import attempted within various core operations if ``cftime`` not
  installed (:pull:`5640`).
  By `Luke Sewell <https://github.com/lusewell>`_
- Fixed bug when combining named DataArrays using :py:func:`combine_by_coords`. (:pull:`5834`).
  By `Tom Nicholas <https://github.com/TomNicholas>`_.
- When a custom engine was used in :py:func:`~xarray.open_dataset` the engine
  wasn't initialized properly, causing missing argument errors or inconsistent
  method signatures. (:pull:`5684`)
  By `Jimmy Westling <https://github.com/illviljan>`_.
- Numbers are properly formatted in a plot's title (:issue:`5788`, :pull:`5789`).
  By `Maxime Liquet <https://github.com/maximlt>`_.
- Faceted plots will no longer raise a `pint.UnitStrippedWarning` when a `pint.Quantity` array is plotted,
  and will correctly display the units of the data in the colorbar (if there is one) (:pull:`5886`).
  By `Tom Nicholas <https://github.com/TomNicholas>`_.
- With backends, check for path-like objects rather than ``pathlib.Path``
  type, use ``os.fspath`` (:pull:`5879`).
  By `Mike Taves <https://github.com/mwtoews>`_.
- ``open_mfdataset()`` now accepts a single ``pathlib.Path`` object (:issue: `5881`).
  By `Panos Mavrogiorgos <https://github.com/pmav99>`_.
- Improved performance of :py:meth:`Dataset.unstack` (:pull:`5906`). By `Tom Augspurger <https://github.com/TomAugspurger>`_.

Documentation
~~~~~~~~~~~~~

- Users are instructed to try ``use_cftime=True`` if a ``TypeError`` occurs when combining datasets and one of the types involved is a subclass of ``cftime.datetime`` (:pull:`5776`).
  By `Zeb Nicholls <https://github.com/znicholls>`_.
- A clearer error is now raised if a user attempts to assign a Dataset to a single key of
  another Dataset. (:pull:`5839`)
  By `Tom Nicholas <https://github.com/TomNicholas>`_.

Internal Changes
~~~~~~~~~~~~~~~~

- Explicit indexes refactor: avoid ``len(index)`` in ``map_blocks`` (:pull:`5670`).
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Explicit indexes refactor: decouple ``xarray.Index``` from ``xarray.Variable`` (:pull:`5636`).
  By `Benoit Bovy <https://github.com/benbovy>`_.
- Fix ``Mapping`` argument typing to allow mypy to pass on ``str`` keys (:pull:`5690`).
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- Annotate many of our tests, and fix some of the resulting typing errors. This will
  also mean our typing annotations are tested as part of CI. (:pull:`5728`).
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- Improve the performance of reprs for large datasets or dataarrays. (:pull:`5661`)
  By `Jimmy Westling <https://github.com/illviljan>`_.
- Use isort's `float_to_top` config. (:pull:`5695`).
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- Remove use of the deprecated ``kind`` argument in
  :py:meth:`pandas.Index.get_slice_bound` inside :py:class:`xarray.CFTimeIndex`
  tests (:pull:`5723`).  By `Spencer Clark <https://github.com/spencerkclark>`_.
- Refactor `xarray.core.duck_array_ops` to no longer special-case dispatching to
  dask versions of functions when acting on dask arrays, instead relying numpy
  and dask's adherence to NEP-18 to dispatch automatically. (:pull:`5571`)
  By `Tom Nicholas <https://github.com/TomNicholas>`_.
- Add an ASV benchmark CI and improve performance of the benchmarks (:pull:`5796`)
  By `Jimmy Westling <https://github.com/illviljan>`_.
- Use ``importlib`` to replace functionality of ``pkg_resources`` such
  as version setting and loading of resources. (:pull:`5845`).
  By `Martin K. Scherer <https://github.com/marscher>`_.


.. _whats-new.0.19.0:

v0.19.0 (23 July 2021)
----------------------

This release brings improvements to plotting of categorical data, the ability to specify how attributes
are combined in xarray operations, a new high-level :py:func:`unify_chunks` function, as well as various
deprecations, bug fixes, and minor improvements.


Many thanks to the 29 contributors to this release!:

Andrew Williams, Augustus, Aureliana Barghini, Benoit Bovy, crusaderky, Deepak Cherian, ellesmith88,
Elliott Sales de Andrade, Giacomo Caria, github-actions[bot], Illviljan, Joeperdefloep, joooeey, Julia Kent,
Julius Busecke, keewis, Mathias Hauser, Matthias Göbel, Mattia Almansi, Maximilian Roos, Peter Andreas Entschev,
Ray Bell, Sander, Santiago Soler, Sebastian, Spencer Clark, Stephan Hoyer, Thomas Hirtz, Thomas Nicholas.

New Features
~~~~~~~~~~~~
- Allow passing argument ``missing_dims`` to :py:meth:`Variable.transpose` and :py:meth:`Dataset.transpose`
  (:issue:`5550`, :pull:`5586`)
  By `Giacomo Caria <https://github.com/gcaria>`_.
- Allow passing a dictionary as coords to a :py:class:`DataArray` (:issue:`5527`,
  reverts :pull:`1539`, which had deprecated this due to python's inconsistent ordering in earlier versions).
  By `Sander van Rijn <https://github.com/sjvrijn>`_.
- Added :py:meth:`Dataset.coarsen.construct`, :py:meth:`DataArray.coarsen.construct` (:issue:`5454`, :pull:`5475`).
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Xarray now uses consolidated metadata by default when writing and reading Zarr
  stores (:issue:`5251`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- New top-level function :py:func:`unify_chunks`.
  By `Mattia Almansi <https://github.com/malmans2>`_.
- Allow assigning values to a subset of a dataset using positional or label-based
  indexing (:issue:`3015`, :pull:`5362`).
  By `Matthias Göbel <https://github.com/matzegoebel>`_.
- Attempting to reduce a weighted object over missing dimensions now raises an error (:pull:`5362`).
  By `Mattia Almansi <https://github.com/malmans2>`_.
- Add ``.sum`` to :py:meth:`~xarray.DataArray.rolling_exp` and
  :py:meth:`~xarray.Dataset.rolling_exp` for exponentially weighted rolling
  sums. These require numbagg 0.2.1;
  (:pull:`5178`).
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- :py:func:`xarray.cov` and :py:func:`xarray.corr` now lazily check for missing
  values if inputs are dask arrays (:issue:`4804`, :pull:`5284`).
  By `Andrew Williams <https://github.com/AndrewWilliams3142>`_.
- Attempting to ``concat`` list of elements that are not all ``Dataset`` or all ``DataArray`` now raises an error (:issue:`5051`, :pull:`5425`).
  By `Thomas Hirtz <https://github.com/thomashirtz>`_.
- allow passing a function to ``combine_attrs`` (:pull:`4896`).
  By `Justus Magin <https://github.com/keewis>`_.
- Allow plotting categorical data (:pull:`5464`).
  By `Jimmy Westling <https://github.com/illviljan>`_.
- Allow removal of the coordinate attribute ``coordinates`` on variables by setting ``.attrs['coordinates']= None``
  (:issue:`5510`).
  By `Elle Smith <https://github.com/ellesmith88>`_.
- Added :py:meth:`DataArray.to_numpy`, :py:meth:`DataArray.as_numpy`, and :py:meth:`Dataset.as_numpy`. (:pull:`5568`).
  By `Tom Nicholas <https://github.com/TomNicholas>`_.
- Units in plot labels are now automatically inferred from wrapped :py:meth:`pint.Quantity` arrays. (:pull:`5561`).
  By `Tom Nicholas <https://github.com/TomNicholas>`_.

Breaking changes
~~~~~~~~~~~~~~~~

- The default ``mode`` for :py:meth:`Dataset.to_zarr` when ``region`` is set
  has changed to the new ``mode="r+"``, which only allows for overriding
  pre-existing array values. This is a safer default than the prior ``mode="a"``,
  and allows for higher performance writes (:pull:`5252`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- The main parameter to :py:func:`combine_by_coords` is renamed to `data_objects` instead
  of `datasets` so anyone calling this method using a named parameter will need to update
  the name accordingly (:issue:`3248`, :pull:`4696`).
  By `Augustus Ijams <https://github.com/aijams>`_.

Deprecations
~~~~~~~~~~~~

- Removed the deprecated ``dim`` kwarg to :py:func:`DataArray.integrate` (:pull:`5630`)
- Removed the deprecated ``keep_attrs`` kwarg to :py:func:`DataArray.rolling` (:pull:`5630`)
- Removed the deprecated ``keep_attrs`` kwarg to :py:func:`DataArray.coarsen` (:pull:`5630`)
- Completed deprecation of passing an ``xarray.DataArray`` to :py:func:`Variable` - will now raise a ``TypeError`` (:pull:`5630`)

Bug fixes
~~~~~~~~~
- Fix a minor incompatibility between partial datetime string indexing with a
  :py:class:`CFTimeIndex` and upcoming pandas version 1.3.0 (:issue:`5356`,
  :pull:`5359`).
  By `Spencer Clark <https://github.com/spencerkclark>`_.
- Fix 1-level multi-index incorrectly converted to single index (:issue:`5384`,
  :pull:`5385`).
  By `Benoit Bovy <https://github.com/benbovy>`_.
- Don't cast a duck array in a coordinate to :py:class:`numpy.ndarray` in
  :py:meth:`DataArray.differentiate` (:pull:`5408`)
  By `Justus Magin <https://github.com/keewis>`_.
- Fix the ``repr`` of :py:class:`Variable` objects with ``display_expand_data=True``
  (:pull:`5406`)
  By `Justus Magin <https://github.com/keewis>`_.
- Plotting a pcolormesh with ``xscale="log"`` and/or ``yscale="log"`` works as
  expected after improving the way the interval breaks are generated (:issue:`5333`).
  By `Santiago Soler <https://github.com/santisoler>`_
- :py:func:`combine_by_coords` can now handle combining a list of unnamed
  ``DataArray`` as input (:issue:`3248`, :pull:`4696`).
  By `Augustus Ijams <https://github.com/aijams>`_.


Internal Changes
~~~~~~~~~~~~~~~~
- Run CI on the first & last python versions supported only; currently 3.7 & 3.9.
  (:pull:`5433`)
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- Publish test results & timings on each PR.
  (:pull:`5537`)
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- Explicit indexes refactor: add a ``xarray.Index.query()`` method in which
  one may eventually provide a custom implementation of label-based data
  selection (not ready yet for public use). Also refactor the internal,
  pandas-specific implementation into ``PandasIndex.query()`` and
  ``PandasMultiIndex.query()`` (:pull:`5322`).
  By `Benoit Bovy <https://github.com/benbovy>`_.

.. _whats-new.0.18.2:

v0.18.2 (19 May 2021)
---------------------

This release reverts a regression in xarray's unstacking of dask-backed arrays.

.. _whats-new.0.18.1:

v0.18.1 (18 May 2021)
---------------------

This release is intended as a small patch release to be compatible with the new
2021.5.0 ``dask.distributed`` release. It also includes a new
``drop_duplicates`` method, some documentation improvements, the beginnings of
our internal Index refactoring, and some bug fixes.

Thank you to all 16 contributors!

Anderson Banihirwe, Andrew, Benoit Bovy, Brewster Malevich, Giacomo Caria,
Illviljan, James Bourbeau, Keewis, Maximilian Roos, Ravin Kumar, Stephan Hoyer,
Thomas Nicholas, Tom Nicholas, Zachary Moon.

New Features
~~~~~~~~~~~~
- Implement :py:meth:`DataArray.drop_duplicates`
  to remove duplicate dimension values (:pull:`5239`).
  By `Andrew Huang <https://github.com/ahuang11>`_.
- Allow passing ``combine_attrs`` strategy names to the ``keep_attrs`` parameter of
  :py:func:`apply_ufunc` (:pull:`5041`)
  By `Justus Magin <https://github.com/keewis>`_.
- :py:meth:`Dataset.interp` now allows interpolation with non-numerical datatypes,
  such as booleans, instead of dropping them. (:issue:`4761` :pull:`5008`).
  By `Jimmy Westling <https://github.com/illviljan>`_.
- Raise more informative error when decoding time variables with invalid reference dates.
  (:issue:`5199`, :pull:`5288`). By `Giacomo Caria <https://github.com/gcaria>`_.


Bug fixes
~~~~~~~~~
- Opening netCDF files from a path that doesn't end in ``.nc`` without supplying
  an explicit ``engine`` works again (:issue:`5295`), fixing a bug introduced in
  0.18.0.
  By `Stephan Hoyer <https://github.com/shoyer>`_

Documentation
~~~~~~~~~~~~~
- Clean up and enhance docstrings for the :py:class:`DataArray.plot` and ``Dataset.plot.*``
  families of methods (:pull:`5285`).
  By `Zach Moon <https://github.com/zmoon>`_.

- Explanation of deprecation cycles and how to implement them added to contributors
  guide. (:pull:`5289`)
  By `Tom Nicholas <https://github.com/TomNicholas>`_.


Internal Changes
~~~~~~~~~~~~~~~~

- Explicit indexes refactor: add an ``xarray.Index`` base class and
  ``Dataset.xindexes`` / ``DataArray.xindexes`` properties. Also rename
  ``PandasIndexAdapter`` to ``PandasIndex``, which now inherits from
  ``xarray.Index`` (:pull:`5102`).
  By `Benoit Bovy <https://github.com/benbovy>`_.
- Replace ``SortedKeysDict`` with python's ``dict``, given dicts are now ordered.
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- Updated the release guide for developers. Now accounts for actions that are automated via github
  actions. (:pull:`5274`).
  By `Tom Nicholas <https://github.com/TomNicholas>`_.

.. _whats-new.0.18.0:

v0.18.0 (6 May 2021)
--------------------

This release brings a few important performance improvements, a wide range of
usability upgrades, lots of bug fixes, and some new features. These include
a plugin API to add backend engines, a new theme for the documentation,
curve fitting methods, and several new plotting functions.

Many thanks to the 38 contributors to this release: Aaron Spring, Alessandro Amici,
Alex Marandon, Alistair Miles, Ana Paula Krelling, Anderson Banihirwe, Aureliana Barghini,
Baudouin Raoult, Benoit Bovy, Blair Bonnett, David Trémouilles, Deepak Cherian,
Gabriel Medeiros Abrahão, Giacomo Caria, Hauke Schulz, Illviljan, Mathias Hauser, Matthias Bussonnier,
Mattia Almansi, Maximilian Roos, Ray Bell, Richard Kleijn, Ryan Abernathey, Sam Levang, Spencer Clark,
Spencer Jones, Tammas Loughran, Tobias Kölling, Todd, Tom Nicholas, Tom White, Victor Negîrneac,
Xianxiang Li, Zeb Nicholls, crusaderky, dschwoerer, johnomotani, keewis


New Features
~~~~~~~~~~~~

- apply ``combine_attrs`` on data variables and coordinate variables when concatenating
  and merging datasets and dataarrays (:pull:`4902`).
  By `Justus Magin <https://github.com/keewis>`_.
- Add :py:meth:`Dataset.to_pandas` (:pull:`5247`)
  By `Giacomo Caria <https://github.com/gcaria>`_.
- Add :py:meth:`DataArray.plot.surface` which wraps matplotlib's `plot_surface` to make
  surface plots (:issue:`2235` :issue:`5084` :pull:`5101`).
  By `John Omotani <https://github.com/johnomotani>`_.
- Allow passing multiple arrays to :py:meth:`Dataset.__setitem__` (:pull:`5216`).
  By `Giacomo Caria <https://github.com/gcaria>`_.
- Add 'cumulative' option to :py:meth:`Dataset.integrate` and
  :py:meth:`DataArray.integrate` so that result is a cumulative integral, like
  :py:func:`scipy.integrate.cumulative_trapezoidal` (:pull:`5153`).
  By `John Omotani <https://github.com/johnomotani>`_.
- Add ``safe_chunks`` option to :py:meth:`Dataset.to_zarr` which allows overriding
  checks made to ensure Dask and Zarr chunk compatibility (:issue:`5056`).
  By `Ryan Abernathey <https://github.com/rabernat>`_
- Add :py:meth:`Dataset.query` and :py:meth:`DataArray.query` which enable indexing
  of datasets and data arrays by evaluating query expressions against the values of the
  data variables (:pull:`4984`).
  By `Alistair Miles <https://github.com/alimanfoo>`_.
- Allow passing ``combine_attrs`` to :py:meth:`Dataset.merge` (:pull:`4895`).
  By `Justus Magin <https://github.com/keewis>`_.
- Support for `dask.graph_manipulation
  <https://docs.dask.org/en/latest/graph_manipulation.html>`_ (requires dask >=2021.3)
  By `Guido Imperiale <https://github.com/crusaderky>`_
- Add :py:meth:`Dataset.plot.streamplot` for streamplot plots with :py:class:`Dataset`
  variables (:pull:`5003`).
  By `John Omotani <https://github.com/johnomotani>`_.
- Many of the arguments for the :py:attr:`DataArray.str` methods now support
  providing an array-like input. In this case, the array provided to the
  arguments is broadcast against the original array and applied elementwise.
- :py:attr:`DataArray.str` now supports ``+``, ``*``, and ``%`` operators. These
  behave the same as they do for :py:class:`str`, except that they follow
  array broadcasting rules.
- A large number of new :py:attr:`DataArray.str` methods were implemented,
  :py:meth:`DataArray.str.casefold`, :py:meth:`DataArray.str.cat`,
  :py:meth:`DataArray.str.extract`, :py:meth:`DataArray.str.extractall`,
  :py:meth:`DataArray.str.findall`, :py:meth:`DataArray.str.format`,
  :py:meth:`DataArray.str.get_dummies`, :py:meth:`DataArray.str.islower`,
  :py:meth:`DataArray.str.join`, :py:meth:`DataArray.str.normalize`,
  :py:meth:`DataArray.str.partition`, :py:meth:`DataArray.str.rpartition`,
  :py:meth:`DataArray.str.rsplit`, and  :py:meth:`DataArray.str.split`.
  A number of these methods allow for splitting or joining the strings in an
  array. (:issue:`4622`)
  By `Todd Jennings <https://github.com/toddrjen>`_
- Thanks to the new pluggable backend infrastructure external packages may now
  use the ``xarray.backends`` entry point to register additional engines to be used in
  :py:func:`open_dataset`, see the documentation in :ref:`add_a_backend`
  (:issue:`4309`, :issue:`4803`, :pull:`4989`, :pull:`4810` and many others).
  The backend refactor has been sponsored with the "Essential Open Source Software for Science"
  grant from the `Chan Zuckerberg Initiative <https://chanzuckerberg.com>`_ and
  developed by `B-Open <https://www.bopen.eu>`_.
  By `Aureliana Barghini <https://github.com/aurghs>`_ and `Alessandro Amici <https://github.com/alexamici>`_.
- :py:attr:`~core.accessor_dt.DatetimeAccessor.date` added (:issue:`4983`, :pull:`4994`).
  By `Hauke Schulz <https://github.com/observingClouds>`_.
- Implement ``__getitem__`` for both :py:class:`~core.groupby.DatasetGroupBy` and
  :py:class:`~core.groupby.DataArrayGroupBy`, inspired by pandas'
  :py:meth:`~pandas.core.groupby.GroupBy.get_group`.
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Switch the tutorial functions to use `pooch <https://github.com/fatiando/pooch>`_
  (which is now a optional dependency) and add :py:func:`tutorial.open_rasterio` as a
  way to open example rasterio files (:issue:`3986`, :pull:`4102`, :pull:`5074`).
  By `Justus Magin <https://github.com/keewis>`_.
- Add typing information to unary and binary arithmetic operators operating on
  :py:class:`Dataset`, :py:class:`DataArray`, :py:class:`Variable`,
  :py:class:`~core.groupby.DatasetGroupBy` or
  :py:class:`~core.groupby.DataArrayGroupBy` (:pull:`4904`).
  By `Richard Kleijn <https://github.com/rhkleijn>`_.
- Add a ``combine_attrs`` parameter to :py:func:`open_mfdataset` (:pull:`4971`).
  By `Justus Magin <https://github.com/keewis>`_.
- Enable passing arrays with a subset of dimensions to
  :py:meth:`DataArray.clip` & :py:meth:`Dataset.clip`; these methods now use
  :py:func:`xarray.apply_ufunc`; (:pull:`5184`).
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- Disable the `cfgrib` backend if the `eccodes` library is not installed (:pull:`5083`).
  By `Baudouin Raoult <https://github.com/b8raoult>`_.
- Added :py:meth:`DataArray.curvefit` and :py:meth:`Dataset.curvefit` for general curve fitting applications. (:issue:`4300`, :pull:`4849`)
  By `Sam Levang <https://github.com/slevang>`_.
- Add options to control expand/collapse of sections in display of Dataset and
  DataArray. The function :py:func:`set_options` now takes keyword arguments
  ``display_expand_attrs``, ``display_expand_coords``, ``display_expand_data``,
  ``display_expand_data_vars``, all of which can be one of ``True`` to always
  expand, ``False`` to always collapse, or ``default`` to expand unless over a
  pre-defined limit (:pull:`5126`).
  By `Tom White <https://github.com/tomwhite>`_.
- Significant speedups in :py:meth:`Dataset.interp` and :py:meth:`DataArray.interp`.
  (:issue:`4739`, :pull:`4740`).
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Prevent passing `concat_dim` to :py:func:`xarray.open_mfdataset` when
  `combine='by_coords'` is specified, which should never have been possible (as
  :py:func:`xarray.combine_by_coords` has no `concat_dim` argument to pass to).
  Also removes unneeded internal reordering of datasets in
  :py:func:`xarray.open_mfdataset` when `combine='by_coords'` is specified.
  Fixes (:issue:`5230`).
  By `Tom Nicholas <https://github.com/TomNicholas>`_.
- Implement ``__setitem__`` for ``xarray.core.indexing.DaskIndexingAdapter`` if
  dask version supports item assignment. (:issue:`5171`, :pull:`5174`)
  By `Tammas Loughran <https://github.com/tammasloughran>`_.

Breaking changes
~~~~~~~~~~~~~~~~
- The minimum versions of some dependencies were changed:

  ============ ====== ====
  Package      Old    New
  ============ ====== ====
  boto3        1.12   1.13
  cftime       1.0    1.1
  dask         2.11   2.15
  distributed  2.11   2.15
  matplotlib   3.1    3.2
  numba        0.48   0.49
  ============ ====== ====

- :py:func:`open_dataset` and :py:func:`open_dataarray` now accept only the first argument
  as positional, all others need to be passed are keyword arguments. This is part of the
  refactor to support external backends (:issue:`4309`, :pull:`4989`).
  By `Alessandro Amici <https://github.com/alexamici>`_.
- Functions that are identities for 0d data return the unchanged data
  if axis is empty. This ensures that Datasets where some variables do
  not have the averaged dimensions are not accidentially changed
  (:issue:`4885`, :pull:`5207`).
  By `David Schwörer <https://github.com/dschwoerer>`_.
- :py:attr:`DataArray.coarsen` and :py:attr:`Dataset.coarsen` no longer support passing ``keep_attrs``
  via its constructor. Pass ``keep_attrs`` via the applied function, i.e. use
  ``ds.coarsen(...).mean(keep_attrs=False)`` instead of ``ds.coarsen(..., keep_attrs=False).mean()``.
  Further, coarsen now keeps attributes per default (:pull:`5227`).
  By `Mathias Hauser <https://github.com/mathause>`_.
- switch the default of the :py:func:`merge` ``combine_attrs`` parameter to
  ``"override"``. This will keep the current behavior for merging the ``attrs`` of
  variables but stop dropping the ``attrs`` of the main objects (:pull:`4902`).
  By `Justus Magin <https://github.com/keewis>`_.

Deprecations
~~~~~~~~~~~~

- Warn when passing `concat_dim` to :py:func:`xarray.open_mfdataset` when
  `combine='by_coords'` is specified, which should never have been possible (as
  :py:func:`xarray.combine_by_coords` has no `concat_dim` argument to pass to).
  Also removes unneeded internal reordering of datasets in
  :py:func:`xarray.open_mfdataset` when `combine='by_coords'` is specified.
  Fixes (:issue:`5230`), via (:pull:`5231`, :pull:`5255`).
  By `Tom Nicholas <https://github.com/TomNicholas>`_.
- The `lock` keyword argument to :py:func:`open_dataset` and :py:func:`open_dataarray` is now
  a backend specific option. It will give a warning if passed to a backend that doesn't support it
  instead of being silently ignored. From the next version it will raise an error.
  This is part of the refactor to support external backends (:issue:`5073`).
  By `Tom Nicholas <https://github.com/TomNicholas>`_ and `Alessandro Amici <https://github.com/alexamici>`_.


Bug fixes
~~~~~~~~~
- Properly support :py:meth:`DataArray.ffill`, :py:meth:`DataArray.bfill`, :py:meth:`Dataset.ffill`, :py:meth:`Dataset.bfill` along chunked dimensions.
  (:issue:`2699`).
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Fix 2d plot failure for certain combinations of dimensions when `x` is 1d and `y` is
  2d (:issue:`5097`, :pull:`5099`).
  By `John Omotani <https://github.com/johnomotani>`_.
- Ensure standard calendar times encoded with large values (i.e. greater than
  approximately 292 years), can be decoded correctly without silently overflowing
  (:pull:`5050`).  This was a regression in xarray 0.17.0.
  By `Zeb Nicholls <https://github.com/znicholls>`_.
- Added support for `numpy.bool_` attributes in roundtrips using `h5netcdf` engine with `invalid_netcdf=True` [which casts `bool`s to `numpy.bool_`] (:issue:`4981`, :pull:`4986`).
  By `Victor Negîrneac <https://github.com/caenrigen>`_.
- Don't allow passing ``axis`` to :py:meth:`Dataset.reduce` methods (:issue:`3510`, :pull:`4940`).
  By `Justus Magin <https://github.com/keewis>`_.
- Decode values as signed if attribute `_Unsigned = "false"` (:issue:`4954`)
  By `Tobias Kölling <https://github.com/d70-t>`_.
- Keep coords attributes when interpolating when the indexer is not a Variable. (:issue:`4239`, :issue:`4839` :pull:`5031`)
  By `Jimmy Westling <https://github.com/illviljan>`_.
- Ensure standard calendar dates encoded with a calendar attribute with some or
  all uppercase letters can be decoded or encoded to or from
  ``np.datetime64[ns]`` dates with or without ``cftime`` installed
  (:issue:`5093`, :pull:`5180`).
  By `Spencer Clark <https://github.com/spencerkclark>`_.
- Warn on passing ``keep_attrs`` to ``resample`` and ``rolling_exp`` as they are ignored, pass ``keep_attrs``
  to the applied function instead (:pull:`5265`).
  By `Mathias Hauser <https://github.com/mathause>`_.

Documentation
~~~~~~~~~~~~~
- New section on :ref:`add_a_backend` in the "Internals" chapter aimed to backend developers
  (:issue:`4803`, :pull:`4810`).
  By `Aureliana Barghini <https://github.com/aurghs>`_.
- Add :py:meth:`Dataset.polyfit` and :py:meth:`DataArray.polyfit` under "See also" in
  the docstrings of :py:meth:`Dataset.polyfit` and :py:meth:`DataArray.polyfit`
  (:issue:`5016`, :pull:`5020`).
  By `Aaron Spring <https://github.com/aaronspring>`_.
- New sphinx theme & rearrangement of the docs (:pull:`4835`).
  By `Anderson Banihirwe <https://github.com/andersy005>`_.

Internal Changes
~~~~~~~~~~~~~~~~
- Enable displaying mypy error codes and ignore only specific error codes using
  ``# type: ignore[error-code]`` (:pull:`5096`).
  By `Mathias Hauser <https://github.com/mathause>`_.
- Replace uses of ``raises_regex`` with the more standard
  ``pytest.raises(Exception, match="foo")``;
  (:pull:`5188`), (:pull:`5191`).
  By `Maximilian Roos <https://github.com/max-sixty>`_.

.. _whats-new.0.17.0:

v0.17.0 (24 Feb 2021)
---------------------

This release brings a few important performance improvements, a wide range of
usability upgrades, lots of bug fixes, and some new features. These include
better ``cftime`` support, a new quiver plot, better ``unstack`` performance,
more efficient memory use in rolling operations, and some python packaging
improvements. We also have a few documentation improvements (and more planned!).

Many thanks to the 36 contributors to this release: Alessandro Amici, Anderson
Banihirwe, Aureliana Barghini, Ayrton Bourn, Benjamin Bean, Blair Bonnett, Chun
Ho Chow, DWesl, Daniel Mesejo-León, Deepak Cherian, Eric Keenan, Illviljan, Jens
Hedegaard Nielsen, Jody Klymak, Julien Seguinot, Julius Busecke, Kai Mühlbauer,
Leif Denby, Martin Durant, Mathias Hauser, Maximilian Roos, Michael Mann, Ray
Bell, RichardScottOZ, Spencer Clark, Tim Gates, Tom Nicholas, Yunus Sevinchan,
alexamici, aurghs, crusaderky, dcherian, ghislainp, keewis, rhkleijn

Breaking changes
~~~~~~~~~~~~~~~~
- xarray no longer supports python 3.6

  The minimum version policy was changed to also apply to projects with irregular
  releases. As a result, the minimum versions of some dependencies have changed:

  ============ ====== ====
  Package      Old    New
  ============ ====== ====
  Python       3.6    3.7
  setuptools   38.4   40.4
  numpy        1.15   1.17
  pandas       0.25   1.0
  dask         2.9    2.11
  distributed  2.9    2.11
  bottleneck   1.2    1.3
  h5netcdf     0.7    0.8
  iris         2.2    2.4
  netcdf4      1.4    1.5
  pseudonetcdf 3.0    3.1
  rasterio     1.0    1.1
  scipy        1.3    1.4
  seaborn      0.9    0.10
  zarr         2.3    2.4
  ============ ====== ====

  (:issue:`4688`, :pull:`4720`, :pull:`4907`, :pull:`4942`)
- As a result of :pull:`4684` the default units encoding for
  datetime-like values (``np.datetime64[ns]`` or ``cftime.datetime``) will now
  always be set such that ``int64`` values can be used.  In the past, no units
  finer than "seconds" were chosen, which would sometimes mean that ``float64``
  values were required, which would lead to inaccurate I/O round-trips.
- Variables referred to in attributes like ``bounds`` and ``grid_mapping``
  can be set as coordinate variables. These attributes are moved to
  :py:attr:`DataArray.encoding` from :py:attr:`DataArray.attrs`. This behaviour
  is controlled by the ``decode_coords`` kwarg to :py:func:`open_dataset` and
  :py:func:`open_mfdataset`.  The full list of decoded attributes is in
  :ref:`weather-climate` (:pull:`2844`, :issue:`3689`)
- As a result of :pull:`4911` the output from calling :py:meth:`DataArray.sum`
  or :py:meth:`DataArray.prod` on an integer array with ``skipna=True`` and a
  non-None value for ``min_count`` will now be a float array rather than an
  integer array.

Deprecations
~~~~~~~~~~~~

- ``dim`` argument to :py:meth:`DataArray.integrate` is being deprecated in
  favour of a ``coord`` argument, for consistency with :py:meth:`Dataset.integrate`.
  For now using ``dim`` issues a ``FutureWarning``. It will be removed in
  version 0.19.0 (:pull:`3993`).
  By `Tom Nicholas <https://github.com/TomNicholas>`_.
- Deprecated ``autoclose`` kwargs from :py:func:`open_dataset` are removed (:pull:`4725`).
  By `Aureliana Barghini <https://github.com/aurghs>`_.
- the return value of :py:meth:`Dataset.update` is being deprecated to make it work more
  like :py:meth:`dict.update`. It will be removed in version 0.19.0 (:pull:`4932`).
  By `Justus Magin <https://github.com/keewis>`_.

New Features
~~~~~~~~~~~~
- :py:meth:`~xarray.cftime_range` and :py:meth:`DataArray.resample` now support
  millisecond (``"L"`` or ``"ms"``) and microsecond (``"U"`` or ``"us"``) frequencies
  for ``cftime.datetime`` coordinates (:issue:`4097`, :pull:`4758`).
  By `Spencer Clark <https://github.com/spencerkclark>`_.
- Significantly higher ``unstack`` performance on numpy-backed arrays which
  contain missing values; 8x faster than previous versions in our benchmark, and
  now 2x faster than pandas (:pull:`4746`).
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- Add :py:meth:`Dataset.plot.quiver` for quiver plots with :py:class:`Dataset` variables.
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Add ``"drop_conflicts"`` to the strategies supported by the ``combine_attrs`` kwarg
  (:issue:`4749`, :pull:`4827`).
  By `Justus Magin <https://github.com/keewis>`_.
- Allow installing from git archives (:pull:`4897`).
  By `Justus Magin <https://github.com/keewis>`_.
- :py:class:`~core.rolling.DataArrayCoarsen` and :py:class:`~core.rolling.DatasetCoarsen`
  now implement a ``reduce`` method, enabling coarsening operations with custom
  reduction functions (:issue:`3741`, :pull:`4939`).
  By `Spencer Clark <https://github.com/spencerkclark>`_.
- Most rolling operations use significantly less memory. (:issue:`4325`).
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Add :py:meth:`Dataset.drop_isel` and :py:meth:`DataArray.drop_isel`
  (:issue:`4658`, :pull:`4819`).
  By `Daniel Mesejo <https://github.com/mesejo>`_.
- Xarray now leverages updates as of cftime version 1.4.1, which enable exact I/O
  roundtripping of ``cftime.datetime`` objects (:pull:`4758`).
  By `Spencer Clark <https://github.com/spencerkclark>`_.
- :py:func:`open_dataset` and :py:func:`open_mfdataset` now accept ``fsspec`` URLs
  (including globs for the latter) for ``engine="zarr"``, and so allow reading from
  many remote and other file systems (:pull:`4461`)
  By `Martin Durant <https://github.com/martindurant>`_
- :py:meth:`DataArray.swap_dims` & :py:meth:`Dataset.swap_dims` now accept dims
  in the form of kwargs as well as a dict, like most similar methods.
  By `Maximilian Roos <https://github.com/max-sixty>`_.

Bug fixes
~~~~~~~~~
- Use specific type checks in ``xarray.core.variable.as_compatible_data`` instead of
  blanket access to ``values`` attribute (:issue:`2097`)
  By `Yunus Sevinchan <https://github.com/blsqr>`_.
- :py:meth:`DataArray.resample` and :py:meth:`Dataset.resample` do not trigger
  computations anymore if :py:meth:`Dataset.weighted` or
  :py:meth:`DataArray.weighted` are applied (:issue:`4625`, :pull:`4668`). By
  `Julius Busecke <https://github.com/jbusecke>`_.
- :py:func:`merge` with ``combine_attrs='override'`` makes a copy of the attrs
  (:issue:`4627`).
- By default, when possible, xarray will now always use values of
  type ``int64`` when encoding and decoding ``numpy.datetime64[ns]`` datetimes.  This
  ensures that maximum precision and accuracy are maintained in the round-tripping
  process (:issue:`4045`, :pull:`4684`). It also enables encoding and decoding standard
  calendar dates with time units of nanoseconds (:pull:`4400`).
  By `Spencer Clark <https://github.com/spencerkclark>`_ and `Mark Harfouche
  <http://github.com/hmaarrfk>`_.
- :py:meth:`DataArray.astype`, :py:meth:`Dataset.astype` and :py:meth:`Variable.astype` support
  the ``order`` and ``subok`` parameters again. This fixes a regression introduced in version 0.16.1
  (:issue:`4644`, :pull:`4683`).
  By `Richard Kleijn <https://github.com/rhkleijn>`_ .
- Remove dictionary unpacking when using ``.loc`` to avoid collision with ``.sel`` parameters (:pull:`4695`).
  By `Anderson Banihirwe <https://github.com/andersy005>`_.
- Fix the legend created by :py:meth:`Dataset.plot.scatter` (:issue:`4641`, :pull:`4723`).
  By `Justus Magin <https://github.com/keewis>`_.
- Fix a crash in orthogonal indexing on geographic coordinates with ``engine='cfgrib'``
  (:issue:`4733` :pull:`4737`).
  By `Alessandro Amici <https://github.com/alexamici>`_.
- Coordinates with dtype ``str`` or ``bytes`` now retain their dtype on many operations,
  e.g. ``reindex``, ``align``, ``concat``, ``assign``, previously they were cast to an object dtype
  (:issue:`2658` and :issue:`4543`).
  By `Mathias Hauser <https://github.com/mathause>`_.
- Limit number of data rows when printing large datasets. (:issue:`4736`, :pull:`4750`).
  By `Jimmy Westling <https://github.com/illviljan>`_.
- Add ``missing_dims`` parameter to transpose (:issue:`4647`, :pull:`4767`).
  By `Daniel Mesejo <https://github.com/mesejo>`_.
- Resolve intervals before appending other metadata to labels when plotting (:issue:`4322`, :pull:`4794`).
  By `Justus Magin <https://github.com/keewis>`_.
- Fix regression when decoding a variable with a ``scale_factor`` and ``add_offset`` given
  as a list of length one (:issue:`4631`).
  By `Mathias Hauser <https://github.com/mathause>`_.
- Expand user directory paths (e.g. ``~/``) in :py:func:`open_mfdataset` and
  :py:meth:`Dataset.to_zarr` (:issue:`4783`, :pull:`4795`).
  By `Julien Seguinot <https://github.com/juseg>`_.
- Raise DeprecationWarning when trying to typecast a tuple containing a :py:class:`DataArray`.
  User now prompted to first call `.data` on it (:issue:`4483`).
  By `Chun Ho Chow <https://github.com/chunhochow>`_.
- Ensure that :py:meth:`Dataset.interp` raises ``ValueError`` when interpolating
  outside coordinate range and ``bounds_error=True`` (:issue:`4854`,
  :pull:`4855`).
  By `Leif Denby <https://github.com/leifdenby>`_.
- Fix time encoding bug associated with using cftime versions greater than
  1.4.0 with xarray (:issue:`4870`, :pull:`4871`).
  By `Spencer Clark <https://github.com/spencerkclark>`_.
- Stop :py:meth:`DataArray.sum` and :py:meth:`DataArray.prod` computing lazy
  arrays when called with a ``min_count`` parameter (:issue:`4898`, :pull:`4911`).
  By `Blair Bonnett <https://github.com/bcbnz>`_.
- Fix bug preventing the ``min_count`` parameter to :py:meth:`DataArray.sum` and
  :py:meth:`DataArray.prod` working correctly when calculating over all axes of
  a float64 array (:issue:`4898`, :pull:`4911`).
  By `Blair Bonnett <https://github.com/bcbnz>`_.
- Fix decoding of vlen strings using h5py versions greater than 3.0.0 with h5netcdf backend (:issue:`4570`, :pull:`4893`).
  By `Kai Mühlbauer <https://github.com/kmuehlbauer>`_.
- Allow converting :py:class:`Dataset` or :py:class:`DataArray` objects with a ``MultiIndex``
  and at least one other dimension to a ``pandas`` object (:issue:`3008`, :pull:`4442`).
  By `ghislainp <https://github.com/ghislainp>`_.

Documentation
~~~~~~~~~~~~~
- Add information about requirements for accessor classes (:issue:`2788`, :pull:`4657`).
  By `Justus Magin <https://github.com/keewis>`_.
- Start a list of external I/O integrating with ``xarray`` (:issue:`683`, :pull:`4566`).
  By `Justus Magin <https://github.com/keewis>`_.
- Add concat examples and improve combining documentation (:issue:`4620`, :pull:`4645`).
  By `Ray Bell <https://github.com/raybellwaves>`_ and
  `Justus Magin <https://github.com/keewis>`_.
- explicitly mention that :py:meth:`Dataset.update` updates inplace (:issue:`2951`, :pull:`4932`).
  By `Justus Magin <https://github.com/keewis>`_.
- Added docs on vectorized indexing (:pull:`4711`).
  By `Eric Keenan <https://github.com/EricKeenan>`_.

Internal Changes
~~~~~~~~~~~~~~~~
- Speed up of the continuous integration tests on azure.

  - Switched to mamba and use matplotlib-base for a faster installation of all dependencies (:pull:`4672`).
  - Use ``pytest.mark.skip`` instead of ``pytest.mark.xfail`` for some tests that can currently not
    succeed (:pull:`4685`).
  - Run the tests in parallel using pytest-xdist (:pull:`4694`).

  By `Justus Magin <https://github.com/keewis>`_ and `Mathias Hauser <https://github.com/mathause>`_.
- Use ``pyproject.toml`` instead of the ``setup_requires`` option for
  ``setuptools`` (:pull:`4897`).
  By `Justus Magin <https://github.com/keewis>`_.
- Replace all usages of ``assert x.identical(y)`` with ``assert_identical(x,  y)``
  for clearer error messages (:pull:`4752`).
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- Speed up attribute style access (e.g. ``ds.somevar`` instead of ``ds["somevar"]``) and
  tab completion in IPython (:issue:`4741`, :pull:`4742`).
  By `Richard Kleijn <https://github.com/rhkleijn>`_.
- Added the ``set_close`` method to ``Dataset`` and ``DataArray`` for backends
  to specify how to voluntary release all resources. (:pull:`#4809`)
  By `Alessandro Amici <https://github.com/alexamici>`_.
- Update type hints to work with numpy v1.20 (:pull:`4878`).
  By `Mathias Hauser <https://github.com/mathause>`_.
- Ensure warnings cannot be turned into exceptions in :py:func:`testing.assert_equal` and
  the other ``assert_*`` functions (:pull:`4864`).
  By `Mathias Hauser <https://github.com/mathause>`_.
- Performance improvement when constructing DataArrays. Significantly speeds up
  repr for Datasets with large number of variables.
  By `Deepak Cherian <https://github.com/dcherian>`_.

.. _whats-new.0.16.2:

v0.16.2 (30 Nov 2020)
---------------------

This release brings the ability to write to limited regions of ``zarr`` files,
open zarr files with :py:func:`open_dataset` and :py:func:`open_mfdataset`,
increased support for propagating ``attrs`` using the ``keep_attrs`` flag, as
well as numerous bugfixes and documentation improvements.

Many thanks to the 31 contributors who contributed to this release: Aaron
Spring, Akio Taniguchi, Aleksandar Jelenak, alexamici, Alexandre Poux, Anderson
Banihirwe, Andrew Pauling, Ashwin Vishnu, aurghs, Brian Ward, Caleb, crusaderky,
Dan Nowacki, darikg, David Brochart, David Huard, Deepak Cherian, Dion Häfner,
Gerardo Rivera, Gerrit Holl, Illviljan, inakleinbottle, Jacob Tomlinson, James
A. Bednar, jenssss, Joe Hamman, johnomotani, Joris Van den Bossche, Julia Kent,
Julius Busecke, Kai Mühlbauer, keewis, Keisuke Fujii, Kyle Cranmer, Luke
Volpatti, Mathias Hauser, Maximilian Roos, Michaël Defferrard, Michal
Baumgartner, Nick R. Papior, Pascal Bourgault, Peter Hausamann, PGijsbers, Ray
Bell, Romain Martinez, rpgoldman, Russell Manser, Sahid Velji, Samnan Rahee,
Sander, Spencer Clark, Stephan Hoyer, Thomas Zilio, Tobias Kölling, Tom
Augspurger, Wei Ji, Yash Saboo, Zeb Nicholls,

Deprecations
~~~~~~~~~~~~

- :py:attr:`~core.accessor_dt.DatetimeAccessor.weekofyear` and :py:attr:`~core.accessor_dt.DatetimeAccessor.week`
  have been deprecated. Use ``DataArray.dt.isocalendar().week``
  instead (:pull:`4534`). By `Mathias Hauser <https://github.com/mathause>`_.
  `Maximilian Roos <https://github.com/max-sixty>`_, and `Spencer Clark <https://github.com/spencerkclark>`_.
- :py:attr:`DataArray.rolling` and :py:attr:`Dataset.rolling` no longer support passing ``keep_attrs``
  via its constructor. Pass ``keep_attrs`` via the applied function, i.e. use
  ``ds.rolling(...).mean(keep_attrs=False)`` instead of ``ds.rolling(..., keep_attrs=False).mean()``
  Rolling operations now keep their attributes per default (:pull:`4510`).
  By `Mathias Hauser <https://github.com/mathause>`_.

New Features
~~~~~~~~~~~~

- :py:func:`open_dataset` and :py:func:`open_mfdataset`
  now works with ``engine="zarr"`` (:issue:`3668`, :pull:`4003`, :pull:`4187`).
  By `Miguel Jimenez <https://github.com/Mikejmnez>`_ and `Wei Ji Leong <https://github.com/weiji14>`_.
- Unary & binary operations follow the ``keep_attrs`` flag (:issue:`3490`, :issue:`4065`, :issue:`3433`, :issue:`3595`, :pull:`4195`).
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Added :py:meth:`~core.accessor_dt.DatetimeAccessor.isocalendar()` that returns a Dataset
  with year, week, and weekday calculated according to the ISO 8601 calendar. Requires
  pandas version 1.1.0 or greater (:pull:`4534`). By `Mathias Hauser <https://github.com/mathause>`_,
  `Maximilian Roos <https://github.com/max-sixty>`_, and `Spencer Clark <https://github.com/spencerkclark>`_.
- :py:meth:`Dataset.to_zarr` now supports a ``region`` keyword for writing to
  limited regions of existing Zarr stores (:pull:`4035`).
  See :ref:`io.zarr.appending` for full details.
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- Added typehints in :py:func:`align` to reflect that the same type received in ``objects`` arg will be returned (:pull:`4522`).
  By `Michal Baumgartner <https://github.com/m1so>`_.
- :py:meth:`Dataset.weighted` and :py:meth:`DataArray.weighted` are now executing value checks lazily if weights are provided as dask arrays (:issue:`4541`, :pull:`4559`).
  By `Julius Busecke <https://github.com/jbusecke>`_.
- Added the ``keep_attrs`` keyword to ``rolling_exp.mean()``; it now keeps attributes
  per default. By `Mathias Hauser <https://github.com/mathause>`_ (:pull:`4592`).
- Added ``freq`` as property to :py:class:`CFTimeIndex` and into the
  ``CFTimeIndex.repr``. (:issue:`2416`, :pull:`4597`)
  By `Aaron Spring <https://github.com/aaronspring>`_.

Bug fixes
~~~~~~~~~

- Fix bug where reference times without padded years (e.g. ``since 1-1-1``) would lose their units when
  being passed by ``encode_cf_datetime`` (:issue:`4422`, :pull:`4506`). Such units are ambiguous
  about which digit represents the years (is it YMD or DMY?). Now, if such formatting is encountered,
  it is assumed that the first digit is the years, they are padded appropriately (to e.g. ``since 0001-1-1``)
  and a warning that this assumption is being made is issued. Previously, without ``cftime``, such times
  would be silently parsed incorrectly (at least based on the CF conventions) e.g. "since 1-1-1" would
  be parsed (via ``pandas`` and ``dateutil``) to ``since 2001-1-1``.
  By `Zeb Nicholls <https://github.com/znicholls>`_.
- Fix :py:meth:`DataArray.plot.step`. By `Deepak Cherian <https://github.com/dcherian>`_.
- Fix bug where reading a scalar value from a NetCDF file opened with the ``h5netcdf`` backend would raise a ``ValueError`` when ``decode_cf=True`` (:issue:`4471`, :pull:`4485`).
  By `Gerrit Holl <https://github.com/gerritholl>`_.
- Fix bug where datetime64 times are silently changed to incorrect values if they are outside the valid date range for ns precision when provided in some other units (:issue:`4427`, :pull:`4454`).
  By `Andrew Pauling <https://github.com/andrewpauling>`_
- Fix silently overwriting the ``engine`` key when passing :py:func:`open_dataset` a file object
  to an incompatible netCDF (:issue:`4457`). Now incompatible combinations of files and engines raise
  an exception instead. By `Alessandro Amici <https://github.com/alexamici>`_.
- The ``min_count`` argument to :py:meth:`DataArray.sum()` and :py:meth:`DataArray.prod()`
  is now ignored when not applicable, i.e. when ``skipna=False`` or when ``skipna=None``
  and the dtype does not have a missing value (:issue:`4352`).
  By `Mathias Hauser <https://github.com/mathause>`_.
- :py:func:`combine_by_coords` now raises an informative error when passing coordinates
  with differing calendars (:issue:`4495`). By `Mathias Hauser <https://github.com/mathause>`_.
- :py:attr:`DataArray.rolling` and :py:attr:`Dataset.rolling` now also keep the attributes and names of of (wrapped)
  ``DataArray`` objects, previously only the global attributes were retained (:issue:`4497`, :pull:`4510`).
  By `Mathias Hauser <https://github.com/mathause>`_.
- Improve performance where reading small slices from huge dimensions was slower than necessary (:pull:`4560`). By `Dion Häfner <https://github.com/dionhaefner>`_.
- Fix bug where ``dask_gufunc_kwargs`` was silently changed in :py:func:`apply_ufunc` (:pull:`4576`). By `Kai Mühlbauer <https://github.com/kmuehlbauer>`_.

Documentation
~~~~~~~~~~~~~
- document the API not supported with duck arrays (:pull:`4530`).
  By `Justus Magin <https://github.com/keewis>`_.
- Mention the possibility to pass functions to :py:meth:`Dataset.where` or
  :py:meth:`DataArray.where` in the parameter documentation (:issue:`4223`, :pull:`4613`).
  By `Justus Magin <https://github.com/keewis>`_.
- Update the docstring of :py:class:`DataArray` and :py:class:`Dataset`.
  (:pull:`4532`);
  By `Jimmy Westling <https://github.com/illviljan>`_.
- Raise a more informative error when :py:meth:`DataArray.to_dataframe` is
  is called on a scalar, (:issue:`4228`);
  By `Pieter Gijsbers <https://github.com/pgijsbers>`_.
- Fix grammar and typos in the :doc:`contributing` guide (:pull:`4545`).
  By `Sahid Velji <https://github.com/sahidvelji>`_.
- Fix grammar and typos in the :doc:`user-guide/io` guide (:pull:`4553`).
  By `Sahid Velji <https://github.com/sahidvelji>`_.
- Update link to NumPy docstring standard in the :doc:`contributing` guide (:pull:`4558`).
  By `Sahid Velji <https://github.com/sahidvelji>`_.
- Add docstrings to ``isnull`` and ``notnull``, and fix the displayed signature
  (:issue:`2760`, :pull:`4618`).
  By `Justus Magin <https://github.com/keewis>`_.

Internal Changes
~~~~~~~~~~~~~~~~

- Optional dependencies can be installed along with xarray by specifying
  extras as ``pip install "xarray[extra]"`` where ``extra`` can be one of ``io``,
  ``accel``, ``parallel``, ``viz`` and ``complete``. See docs for updated
  :ref:`installation instructions <installation-instructions>`.
  (:issue:`2888`, :pull:`4480`).
  By `Ashwin Vishnu <https://github.com/ashwinvis>`_, `Justus Magin
  <https://github.com/keewis>`_ and `Mathias Hauser
  <https://github.com/mathause>`_.
- Removed stray spaces that stem from black removing new lines (:pull:`4504`).
  By `Mathias Hauser <https://github.com/mathause>`_.
- Ensure tests are not skipped in the ``py38-all-but-dask`` test environment
  (:issue:`4509`). By `Mathias Hauser <https://github.com/mathause>`_.
- Ignore select numpy warnings around missing values, where xarray handles
  the values appropriately, (:pull:`4536`);
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- Replace the internal use of ``pd.Index.__or__`` and ``pd.Index.__and__`` with ``pd.Index.union``
  and ``pd.Index.intersection`` as they will stop working as set operations in the future
  (:issue:`4565`). By `Mathias Hauser <https://github.com/mathause>`_.
- Add GitHub action for running nightly tests against upstream dependencies (:pull:`4583`).
  By `Anderson Banihirwe <https://github.com/andersy005>`_.
- Ensure all figures are closed properly in plot tests (:pull:`4600`).
  By `Yash Saboo <https://github.com/yashsaboo>`_, `Nirupam K N
  <https://github.com/Nirupamkn>`_ and `Mathias Hauser
  <https://github.com/mathause>`_.

.. _whats-new.0.16.1:

v0.16.1 (2020-09-20)
---------------------

This patch release fixes an incompatibility with a recent pandas change, which
was causing an issue indexing with a ``datetime64``. It also includes
improvements to ``rolling``, ``to_dataframe``, ``cov`` & ``corr`` methods and
bug fixes. Our documentation has a number of improvements, including fixing all
doctests and confirming their accuracy on every commit.

Many thanks to the 36 contributors who contributed to this release:

Aaron Spring, Akio Taniguchi, Aleksandar Jelenak, Alexandre Poux,
Caleb, Dan Nowacki, Deepak Cherian, Gerardo Rivera, Jacob Tomlinson, James A.
Bednar, Joe Hamman, Julia Kent, Kai Mühlbauer, Keisuke Fujii, Mathias Hauser,
Maximilian Roos, Nick R. Papior, Pascal Bourgault, Peter Hausamann, Romain
Martinez, Russell Manser, Samnan Rahee, Sander, Spencer Clark, Stephan Hoyer,
Thomas Zilio, Tobias Kölling, Tom Augspurger, alexamici, crusaderky, darikg,
inakleinbottle, jenssss, johnomotani, keewis, and rpgoldman.

Breaking changes
~~~~~~~~~~~~~~~~

- :py:meth:`DataArray.astype` and :py:meth:`Dataset.astype` now preserve attributes. Keep the
  old behavior by passing `keep_attrs=False` (:issue:`2049`, :pull:`4314`).
  By `Dan Nowacki <https://github.com/dnowacki-usgs>`_ and `Gabriel Joel Mitchell <https://github.com/gajomi>`_.

New Features
~~~~~~~~~~~~

- :py:meth:`~xarray.DataArray.rolling` and :py:meth:`~xarray.Dataset.rolling`
  now accept more than 1 dimension. (:pull:`4219`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- :py:meth:`~xarray.DataArray.to_dataframe` and :py:meth:`~xarray.Dataset.to_dataframe`
  now accept a ``dim_order`` parameter allowing to specify the resulting dataframe's
  dimensions order (:issue:`4331`, :pull:`4333`).
  By `Thomas Zilio <https://github.com/thomas-z>`_.
- Support multiple outputs in :py:func:`xarray.apply_ufunc` when using
  ``dask='parallelized'``. (:issue:`1815`, :pull:`4060`).
  By `Kai Mühlbauer <https://github.com/kmuehlbauer>`_.
- ``min_count`` can be supplied to reductions such as ``.sum`` when specifying
  multiple dimension to reduce over; (:pull:`4356`).
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- :py:func:`xarray.cov` and :py:func:`xarray.corr` now handle missing values; (:pull:`4351`).
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- Add support for parsing datetime strings formatted following the default
  string representation of cftime objects, i.e. YYYY-MM-DD hh:mm:ss, in
  partial datetime string indexing, as well as :py:meth:`~xarray.cftime_range`
  (:issue:`4337`). By `Spencer Clark <https://github.com/spencerkclark>`_.
- Build ``CFTimeIndex.__repr__`` explicitly as :py:class:`pandas.Index`. Add ``calendar`` as a new
  property for :py:class:`CFTimeIndex` and show ``calendar`` and ``length`` in
  ``CFTimeIndex.__repr__`` (:issue:`2416`, :pull:`4092`)
  By `Aaron Spring <https://github.com/aaronspring>`_.
- Use a wrapped array's ``_repr_inline_`` method to construct the collapsed ``repr``
  of :py:class:`DataArray` and :py:class:`Dataset` objects and
  document the new method in :doc:`internals/index`. (:pull:`4248`).
  By `Justus Magin <https://github.com/keewis>`_.
- Allow per-variable fill values in most functions. (:pull:`4237`).
  By `Justus Magin <https://github.com/keewis>`_.
- Expose ``use_cftime`` option in :py:func:`~xarray.open_zarr` (:issue:`2886`, :pull:`3229`)
  By `Samnan Rahee <https://github.com/Geektrovert>`_ and `Anderson Banihirwe <https://github.com/andersy005>`_.

Bug fixes
~~~~~~~~~

- Fix indexing with datetime64 scalars with pandas 1.1 (:issue:`4283`).
  By `Stephan Hoyer <https://github.com/shoyer>`_ and
  `Justus Magin <https://github.com/keewis>`_.
- Variables which are chunked using dask only along some dimensions can be chunked while storing with zarr along previously
  unchunked dimensions (:pull:`4312`) By `Tobias Kölling <https://github.com/d70-t>`_.
- Fixed a bug in backend caused by basic installation of Dask (:issue:`4164`, :pull:`4318`)
  `Sam Morley <https://github.com/inakleinbottle>`_.
- Fixed a few bugs with :py:meth:`Dataset.polyfit` when encountering deficient matrix ranks (:issue:`4190`, :pull:`4193`). By `Pascal Bourgault <https://github.com/aulemahal>`_.
- Fixed inconsistencies between docstring and functionality for :py:meth:`DataArray.str.get`
  and :py:meth:`DataArray.str.wrap` (:issue:`4334`). By `Mathias Hauser <https://github.com/mathause>`_.
- Fixed overflow issue causing incorrect results in computing means of :py:class:`cftime.datetime`
  arrays (:issue:`4341`). By `Spencer Clark <https://github.com/spencerkclark>`_.
- Fixed :py:meth:`Dataset.coarsen`, :py:meth:`DataArray.coarsen` dropping attributes on original object (:issue:`4120`, :pull:`4360`). By `Julia Kent <https://github.com/jukent>`_.
- fix the signature of the plot methods. (:pull:`4359`) By `Justus Magin <https://github.com/keewis>`_.
- Fix :py:func:`xarray.apply_ufunc` with ``vectorize=True`` and ``exclude_dims`` (:issue:`3890`).
  By `Mathias Hauser <https://github.com/mathause>`_.
- Fix `KeyError` when doing linear interpolation to an nd `DataArray`
  that contains NaNs (:pull:`4233`).
  By `Jens Svensmark <https://github.com/jenssss>`_
- Fix incorrect legend labels for :py:meth:`Dataset.plot.scatter` (:issue:`4126`).
  By `Peter Hausamann <https://github.com/phausamann>`_.
- Fix ``dask.optimize`` on ``DataArray`` producing an invalid Dask task graph (:issue:`3698`)
  By `Tom Augspurger <https://github.com/TomAugspurger>`_
- Fix ``pip install .`` when no ``.git`` directory exists; namely when the xarray source
  directory has been rsync'ed by PyCharm Professional for a remote deployment over SSH.
  By `Guido Imperiale <https://github.com/crusaderky>`_
- Preserve dimension and coordinate order during :py:func:`xarray.concat` (:issue:`2811`, :issue:`4072`, :pull:`4419`).
  By `Kai Mühlbauer <https://github.com/kmuehlbauer>`_.
- Avoid relying on :py:class:`set` objects for the ordering of the coordinates (:pull:`4409`)
  By `Justus Magin <https://github.com/keewis>`_.

Documentation
~~~~~~~~~~~~~

- Update the docstring of :py:meth:`DataArray.copy` to remove incorrect mention of 'dataset' (:issue:`3606`)
  By `Sander van Rijn <https://github.com/sjvrijn>`_.
- Removed skipna argument from :py:meth:`DataArray.count`, :py:meth:`DataArray.any`, :py:meth:`DataArray.all`. (:issue:`755`)
  By `Sander van Rijn <https://github.com/sjvrijn>`_
- Update the contributing guide to use merges instead of rebasing and state
  that we squash-merge. (:pull:`4355`). By `Justus Magin <https://github.com/keewis>`_.
- Make sure the examples from the docstrings actually work (:pull:`4408`).
  By `Justus Magin <https://github.com/keewis>`_.
- Updated Vectorized Indexing to a clearer example.
  By `Maximilian Roos <https://github.com/max-sixty>`_

Internal Changes
~~~~~~~~~~~~~~~~

- Fixed all doctests and enabled their running in CI.
  By `Justus Magin <https://github.com/keewis>`_.
- Relaxed the :ref:`mindeps_policy` to support:

  - all versions of setuptools released in the last 42 months (but no older than 38.4)
  - all versions of dask and dask.distributed released in the last 12 months (but no
    older than 2.9)
  - all versions of other packages released in the last 12 months

  All are up from 6 months (:issue:`4295`)
  `Guido Imperiale <https://github.com/crusaderky>`_.
- Use :py:func:`dask.array.apply_gufunc <dask.array.gufunc.apply_gufunc>` instead of
  :py:func:`dask.array.blockwise` in :py:func:`xarray.apply_ufunc` when using
  ``dask='parallelized'``. (:pull:`4060`, :pull:`4391`, :pull:`4392`)
  By `Kai Mühlbauer <https://github.com/kmuehlbauer>`_.
- Align ``mypy`` versions to ``0.782`` across ``requirements`` and
  ``.pre-commit-config.yml`` files. (:pull:`4390`)
  By `Maximilian Roos <https://github.com/max-sixty>`_
- Only load resource files when running inside a Jupyter Notebook
  (:issue:`4294`) By `Guido Imperiale <https://github.com/crusaderky>`_
- Silenced most ``numpy`` warnings such as ``Mean of empty slice``. (:pull:`4369`)
  By `Maximilian Roos <https://github.com/max-sixty>`_
- Enable type checking for :py:func:`concat` (:issue:`4238`)
  By `Mathias Hauser <https://github.com/mathause>`_.
- Updated plot functions for matplotlib version 3.3 and silenced warnings in the
  plot tests (:pull:`4365`). By `Mathias Hauser <https://github.com/mathause>`_.
- Versions in ``pre-commit.yaml`` are now pinned, to reduce the chances of
  conflicting versions. (:pull:`4388`)
  By `Maximilian Roos <https://github.com/max-sixty>`_



.. _whats-new.0.16.0:

v0.16.0 (2020-07-11)
---------------------

This release adds `xarray.cov` & `xarray.corr` for covariance & correlation
respectively; the `idxmax` & `idxmin` methods, the `polyfit` method &
`xarray.polyval` for fitting polynomials, as well as a number of documentation
improvements, other features, and bug fixes. Many thanks to all 44 contributors
who contributed to this release:

Akio Taniguchi, Andrew Williams, Aurélien Ponte, Benoit Bovy, Dave Cole, David
Brochart, Deepak Cherian, Elliott Sales de Andrade, Etienne Combrisson, Hossein
Madadi, Huite, Joe Hamman, Kai Mühlbauer, Keisuke Fujii, Maik Riechert, Marek
Jacob, Mathias Hauser, Matthieu Ancellin, Maximilian Roos, Noah D Brenowitz,
Oriol Abril, Pascal Bourgault, Phillip Butcher, Prajjwal Nijhara, Ray Bell, Ryan
Abernathey, Ryan May, Spencer Clark, Spencer Hill, Srijan Saurav, Stephan Hoyer,
Taher Chegini, Todd, Tom Nicholas, Yohai Bar Sinai, Yunus Sevinchan,
arabidopsis, aurghs, clausmichele, dmey, johnomotani, keewis, raphael dussin,
risebell

Breaking changes
~~~~~~~~~~~~~~~~

- Minimum supported versions for the following packages have changed: ``dask >=2.9``,
  ``distributed>=2.9``.
  By `Deepak Cherian <https://github.com/dcherian>`_
- ``groupby`` operations will restore coord dimension order. Pass ``restore_coord_dims=False``
  to revert to previous behavior.
- :meth:`DataArray.transpose` will now transpose coordinates by default.
  Pass ``transpose_coords=False`` to revert to previous behaviour.
  By `Maximilian Roos <https://github.com/max-sixty>`_
- Alternate draw styles for :py:meth:`plot.step` must be passed using the
  ``drawstyle`` (or ``ds``) keyword argument, instead of the ``linestyle`` (or
  ``ls``) keyword argument, in line with the `upstream change in Matplotlib
  <https://matplotlib.org/api/prev_api_changes/api_changes_3.1.0.html#passing-a-line2d-s-drawstyle-together-with-the-linestyle-is-deprecated>`_.
  (:pull:`3274`)
  By `Elliott Sales de Andrade <https://github.com/QuLogic>`_
- The old ``auto_combine`` function has now been removed in
  favour of the :py:func:`combine_by_coords` and
  :py:func:`combine_nested` functions. This also means that
  the default behaviour of :py:func:`open_mfdataset` has changed to use
  ``combine='by_coords'`` as the default argument value. (:issue:`2616`, :pull:`3926`)
  By `Tom Nicholas <https://github.com/TomNicholas>`_.
- The ``DataArray`` and ``Variable`` HTML reprs now expand the data section by
  default (:issue:`4176`)
  By `Stephan Hoyer <https://github.com/shoyer>`_.

New Features
~~~~~~~~~~~~
- :py:meth:`DataArray.argmin` and :py:meth:`DataArray.argmax` now support
  sequences of 'dim' arguments, and if a sequence is passed return a dict
  (which can be passed to :py:meth:`DataArray.isel` to get the value of the minimum) of
  the indices for each dimension of the minimum or maximum of a DataArray.
  (:pull:`3936`)
  By `John Omotani <https://github.com/johnomotani>`_, thanks to `Keisuke Fujii
  <https://github.com/fujiisoup>`_ for work in :pull:`1469`.
- Added :py:func:`xarray.cov` and :py:func:`xarray.corr` (:issue:`3784`, :pull:`3550`, :pull:`4089`).
  By `Andrew Williams <https://github.com/AndrewWilliams3142>`_ and `Robin Beer <https://github.com/r-beer>`_.
- Implement :py:meth:`DataArray.idxmax`, :py:meth:`DataArray.idxmin`,
  :py:meth:`Dataset.idxmax`, :py:meth:`Dataset.idxmin`.  (:issue:`60`, :pull:`3871`)
  By `Todd Jennings <https://github.com/toddrjen>`_
- Added :py:meth:`DataArray.polyfit` and :py:func:`xarray.polyval` for fitting
  polynomials. (:issue:`3349`, :pull:`3733`, :pull:`4099`)
  By `Pascal Bourgault <https://github.com/aulemahal>`_.
- Added :py:meth:`xarray.infer_freq` for extending frequency inferring to CFTime indexes and data (:pull:`4033`).
  By `Pascal Bourgault <https://github.com/aulemahal>`_.
- ``chunks='auto'`` is now supported in the ``chunks`` argument of
  :py:meth:`Dataset.chunk`. (:issue:`4055`)
  By `Andrew Williams <https://github.com/AndrewWilliams3142>`_
- Control over attributes of result in :py:func:`merge`, :py:func:`concat`,
  :py:func:`combine_by_coords` and :py:func:`combine_nested` using
  combine_attrs keyword argument. (:issue:`3865`, :pull:`3877`)
  By `John Omotani <https://github.com/johnomotani>`_
- `missing_dims` argument to :py:meth:`Dataset.isel`,
  :py:meth:`DataArray.isel` and :py:meth:`Variable.isel` to allow replacing
  the exception when a dimension passed to ``isel`` is not present with a
  warning, or just ignore the dimension. (:issue:`3866`, :pull:`3923`)
  By `John Omotani <https://github.com/johnomotani>`_
- Support dask handling for :py:meth:`DataArray.idxmax`, :py:meth:`DataArray.idxmin`,
  :py:meth:`Dataset.idxmax`, :py:meth:`Dataset.idxmin`.  (:pull:`3922`, :pull:`4135`)
  By `Kai Mühlbauer <https://github.com/kmuehlbauer>`_ and `Pascal Bourgault <https://github.com/aulemahal>`_.
- More support for unit aware arrays with pint (:pull:`3643`, :pull:`3975`, :pull:`4163`)
  By `Justus Magin <https://github.com/keewis>`_.
- Support overriding existing variables in ``to_zarr()`` with ``mode='a'`` even
  without ``append_dim``, as long as dimension sizes do not change.
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- Allow plotting of boolean arrays. (:pull:`3766`)
  By `Marek Jacob <https://github.com/MeraX>`_
- Enable using MultiIndex levels as coordinates in 1D and 2D plots (:issue:`3927`).
  By `Mathias Hauser <https://github.com/mathause>`_.
- A ``days_in_month`` accessor for :py:class:`xarray.CFTimeIndex`, analogous to
  the ``days_in_month`` accessor for a :py:class:`pandas.DatetimeIndex`, which
  returns the days in the month each datetime in the index.  Now days in month
  weights for both standard and non-standard calendars can be obtained using
  the :py:class:`~core.accessor_dt.DatetimeAccessor` (:pull:`3935`).  This
  feature requires cftime version 1.1.0 or greater.  By
  `Spencer Clark <https://github.com/spencerkclark>`_.
- For the netCDF3 backend, added dtype coercions for unsigned integer types.
  (:issue:`4014`, :pull:`4018`)
  By `Yunus Sevinchan <https://github.com/blsqr>`_
- :py:meth:`map_blocks` now accepts a ``template`` kwarg. This allows use cases
  where the result of a computation could not be inferred automatically.
  By `Deepak Cherian <https://github.com/dcherian>`_
- :py:meth:`map_blocks` can now handle dask-backed xarray objects in ``args``. (:pull:`3818`)
  By `Deepak Cherian <https://github.com/dcherian>`_
- Add keyword ``decode_timedelta`` to :py:func:`xarray.open_dataset`,
  (:py:func:`xarray.open_dataarray`, :py:func:`xarray.open_dataarray`,
  :py:func:`xarray.decode_cf`) that allows to disable/enable the decoding of timedeltas
  independently of time decoding (:issue:`1621`)
  `Aureliana Barghini <https://github.com/aurghs>`_

Enhancements
~~~~~~~~~~~~
- Performance improvement of :py:meth:`DataArray.interp` and :py:func:`Dataset.interp`
  We performs independant interpolation sequentially rather than interpolating in
  one large multidimensional space. (:issue:`2223`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- :py:meth:`DataArray.interp` now support interpolations over chunked dimensions (:pull:`4155`). By `Alexandre Poux <https://github.com/pums974>`_.
- Major performance improvement for :py:meth:`Dataset.from_dataframe` when the
  dataframe has a MultiIndex (:pull:`4184`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.
  - :py:meth:`DataArray.reset_index` and :py:meth:`Dataset.reset_index` now keep
  coordinate attributes (:pull:`4103`). By `Oriol Abril <https://github.com/OriolAbril>`_.
- Axes kwargs such as ``facecolor`` can now be passed to :py:meth:`DataArray.plot` in ``subplot_kws``.
  This works for both single axes plots and FacetGrid plots.
  By `Raphael Dussin <https://github.com/raphaeldussin>`_.
- Array items with long string reprs are now limited to a
  reasonable width (:pull:`3900`)
  By `Maximilian Roos <https://github.com/max-sixty>`_
- Large arrays whose numpy reprs would have greater than 40 lines are now
  limited to a reasonable length.
  (:pull:`3905`)
  By `Maximilian Roos <https://github.com/max-sixty>`_

Bug fixes
~~~~~~~~~
- Fix errors combining attrs in :py:func:`open_mfdataset` (:issue:`4009`, :pull:`4173`)
  By `John Omotani <https://github.com/johnomotani>`_
- If groupby receives a ``DataArray`` with name=None, assign a default name (:issue:`158`)
  By `Phil Butcher <https://github.com/pjbutcher>`_.
- Support dark mode in VS code (:issue:`4024`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- Fix bug when converting multiindexed pandas objects to sparse xarray objects. (:issue:`4019`)
  By `Deepak Cherian <https://github.com/dcherian>`_.
- ``ValueError`` is raised when ``fill_value`` is not a scalar in :py:meth:`full_like`. (:issue:`3977`)
  By `Huite Bootsma <https://github.com/huite>`_.
- Fix wrong order in converting a ``pd.Series`` with a MultiIndex to ``DataArray``.
  (:issue:`3951`, :issue:`4186`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_ and `Stephan Hoyer <https://github.com/shoyer>`_.
- Fix renaming of coords when one or more stacked coords is not in
  sorted order during stack+groupby+apply operations. (:issue:`3287`,
  :pull:`3906`) By `Spencer Hill <https://github.com/spencerahill>`_
- Fix a regression where deleting a coordinate from a copied :py:class:`DataArray`
  can affect the original :py:class:`DataArray`.  (:issue:`3899`, :pull:`3871`)
  By `Todd Jennings <https://github.com/toddrjen>`_
- Fix :py:class:`~xarray.plot.FacetGrid` plots with a single contour. (:issue:`3569`, :pull:`3915`).
  By `Deepak Cherian <https://github.com/dcherian>`_
- Use divergent colormap if ``levels`` spans 0. (:issue:`3524`)
  By `Deepak Cherian <https://github.com/dcherian>`_
- Fix :py:class:`~xarray.plot.FacetGrid` when ``vmin == vmax``. (:issue:`3734`)
  By `Deepak Cherian <https://github.com/dcherian>`_
- Fix plotting when ``levels`` is a scalar and ``norm`` is provided. (:issue:`3735`)
  By `Deepak Cherian <https://github.com/dcherian>`_
- Fix bug where plotting line plots with 2D coordinates depended on dimension
  order. (:issue:`3933`)
  By `Tom Nicholas <https://github.com/TomNicholas>`_.
- Fix ``RasterioDeprecationWarning`` when using a ``vrt`` in ``open_rasterio``. (:issue:`3964`)
  By `Taher Chegini <https://github.com/cheginit>`_.
- Fix ``AttributeError`` on displaying a :py:class:`Variable`
  in a notebook context. (:issue:`3972`, :pull:`3973`)
  By `Ian Castleden <https://github.com/arabidopsis>`_.
- Fix bug causing :py:meth:`DataArray.interpolate_na` to always drop attributes,
  and added `keep_attrs` argument. (:issue:`3968`)
  By `Tom Nicholas <https://github.com/TomNicholas>`_.
- Fix bug in time parsing failing to fall back to cftime. This was causing time
  variables with a time unit of `'msecs'` to fail to parse. (:pull:`3998`)
  By `Ryan May <https://github.com/dopplershift>`_.
- Fix weighted mean when passing boolean weights (:issue:`4074`).
  By `Mathias Hauser <https://github.com/mathause>`_.
- Fix html repr in untrusted notebooks: fallback to plain text repr. (:pull:`4053`)
  By `Benoit Bovy <https://github.com/benbovy>`_.
- Fix :py:meth:`DataArray.to_unstacked_dataset` for single-dimension variables. (:issue:`4049`)
  By `Deepak Cherian <https://github.com/dcherian>`_
- Fix :py:func:`open_rasterio` for ``WarpedVRT`` with specified ``src_crs``. (:pull:`4104`)
  By `Dave Cole <https://github.com/dtpc>`_.

Documentation
~~~~~~~~~~~~~
- update the docstring of :py:meth:`DataArray.assign_coords` : clarify how to
  add a new coordinate to an existing dimension and illustrative example
  (:issue:`3952`, :pull:`3958`) By
  `Etienne Combrisson <https://github.com/EtienneCmb>`_.
- update the docstring of :py:meth:`Dataset.diff` and
  :py:meth:`DataArray.diff` so it does document the ``dim``
  parameter as required. (:issue:`1040`, :pull:`3909`)
  By `Justus Magin <https://github.com/keewis>`_.
- Updated :doc:`Calculating Seasonal Averages from Timeseries of Monthly Means
  <examples/monthly-means>` example notebook to take advantage of the new
  ``days_in_month`` accessor for :py:class:`xarray.CFTimeIndex`
  (:pull:`3935`). By `Spencer Clark <https://github.com/spencerkclark>`_.
- Updated the list of current core developers. (:issue:`3892`)
  By `Tom Nicholas <https://github.com/TomNicholas>`_.
- Add example for multi-dimensional extrapolation and note different behavior
  of ``kwargs`` in :py:meth:`Dataset.interp` and :py:meth:`DataArray.interp`
  for 1-d and n-d interpolation (:pull:`3956`).
  By `Matthias Riße <https://github.com/risebell>`_.
- Apply ``black`` to all the code in the documentation (:pull:`4012`)
  By `Justus Magin <https://github.com/keewis>`_.
- Narrative documentation now describes :py:meth:`map_blocks`: :ref:`dask.automatic-parallelization`.
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Document ``.plot``, ``.dt``, ``.str`` accessors the way they are called. (:issue:`3625`, :pull:`3988`)
  By `Justus Magin <https://github.com/keewis>`_.
- Add documentation for the parameters and return values of :py:meth:`DataArray.sel`.
  By `Justus Magin <https://github.com/keewis>`_.

Internal Changes
~~~~~~~~~~~~~~~~
- Raise more informative error messages for chunk size conflicts when writing to zarr files.
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Run the ``isort`` pre-commit hook only on python source files
  and update the ``flake8`` version. (:issue:`3750`, :pull:`3711`)
  By `Justus Magin <https://github.com/keewis>`_.
- Add `blackdoc <https://blackdoc.readthedocs.io>`_ to the list of
  checkers for development. (:pull:`4177`)
  By `Justus Magin <https://github.com/keewis>`_.
- Add a CI job that runs the tests with every optional dependency
  except ``dask``. (:issue:`3794`, :pull:`3919`)
  By `Justus Magin <https://github.com/keewis>`_.
- Use ``async`` / ``await`` for the asynchronous distributed
  tests. (:issue:`3987`, :pull:`3989`)
  By `Justus Magin <https://github.com/keewis>`_.
- Various internal code clean-ups (:pull:`4026`,  :pull:`4038`).
  By `Prajjwal Nijhara <https://github.com/pnijhara>`_.

.. _whats-new.0.15.1:

v0.15.1 (23 Mar 2020)
---------------------

This release brings many new features such as :py:meth:`Dataset.weighted` methods for weighted array
reductions, a new jupyter repr by default, and the start of units integration with pint. There's also
the usual batch of usability improvements, documentation additions, and bug fixes.

Breaking changes
~~~~~~~~~~~~~~~~

- Raise an error when assigning to the ``.values`` or ``.data`` attribute of
  dimension coordinates i.e. ``IndexVariable`` objects. This has been broken since
  v0.12.0. Please use :py:meth:`DataArray.assign_coords` or :py:meth:`Dataset.assign_coords`
  instead. (:issue:`3470`, :pull:`3862`)
  By `Deepak Cherian <https://github.com/dcherian>`_

New Features
~~~~~~~~~~~~

- Weighted array reductions are now supported via the new :py:meth:`DataArray.weighted`
  and :py:meth:`Dataset.weighted` methods. See :ref:`comput.weighted`. (:issue:`422`, :pull:`2922`).
  By `Mathias Hauser <https://github.com/mathause>`_.
- The new jupyter notebook repr (``Dataset._repr_html_`` and
  ``DataArray._repr_html_``) (introduced in 0.14.1) is now on by default. To
  disable, use ``xarray.set_options(display_style="text")``.
  By `Julia Signell <https://github.com/jsignell>`_.
- Added support for :py:class:`pandas.DatetimeIndex`-style rounding of
  ``cftime.datetime`` objects directly via a :py:class:`CFTimeIndex` or via the
  :py:class:`~core.accessor_dt.DatetimeAccessor`.
  By `Spencer Clark <https://github.com/spencerkclark>`_
- Support new h5netcdf backend keyword `phony_dims` (available from h5netcdf
  v0.8.0 for :py:class:`~xarray.backends.H5NetCDFStore`.
  By `Kai Mühlbauer <https://github.com/kmuehlbauer>`_.
- Add partial support for unit aware arrays with pint. (:pull:`3706`, :pull:`3611`)
  By `Justus Magin <https://github.com/keewis>`_.
- :py:meth:`Dataset.groupby` and :py:meth:`DataArray.groupby` now raise a
  `TypeError` on multiple string arguments. Receiving multiple string arguments
  often means a user is attempting to pass multiple dimensions as separate
  arguments and should instead pass a single list of dimensions.
  (:pull:`3802`)
  By `Maximilian Roos <https://github.com/max-sixty>`_
- :py:func:`map_blocks` can now apply functions that add new unindexed dimensions.
  By `Deepak Cherian <https://github.com/dcherian>`_
- An ellipsis (``...``) is now supported in the ``dims`` argument of
  :py:meth:`Dataset.stack` and :py:meth:`DataArray.stack`, meaning all
  unlisted dimensions, similar to its meaning in :py:meth:`DataArray.transpose`.
  (:pull:`3826`)
  By `Maximilian Roos <https://github.com/max-sixty>`_
- :py:meth:`Dataset.where` and :py:meth:`DataArray.where` accept a lambda as a
  first argument, which is then called on the input; replicating pandas' behavior.
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- ``skipna`` is available in :py:meth:`Dataset.quantile`, :py:meth:`DataArray.quantile`,
  :py:meth:`core.groupby.DatasetGroupBy.quantile`, :py:meth:`core.groupby.DataArrayGroupBy.quantile`
  (:issue:`3843`, :pull:`3844`)
  By `Aaron Spring <https://github.com/aaronspring>`_.
- Add a diff summary for `testing.assert_allclose`. (:issue:`3617`, :pull:`3847`)
  By `Justus Magin <https://github.com/keewis>`_.

Bug fixes
~~~~~~~~~

- Fix :py:meth:`Dataset.interp` when indexing array shares coordinates with the
  indexed variable (:issue:`3252`).
  By `David Huard <https://github.com/huard>`_.
- Fix recombination of groups in :py:meth:`Dataset.groupby` and
  :py:meth:`DataArray.groupby` when performing an operation that changes the
  size of the groups along the grouped dimension. By `Eric Jansen
  <https://github.com/ej81>`_.
- Fix use of multi-index with categorical values (:issue:`3674`).
  By `Matthieu Ancellin <https://github.com/mancellin>`_.
- Fix alignment with ``join="override"`` when some dimensions are unindexed. (:issue:`3681`).
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Fix :py:meth:`Dataset.swap_dims` and :py:meth:`DataArray.swap_dims` producing
  index with name reflecting the previous dimension name instead of the new one
  (:issue:`3748`, :pull:`3752`). By `Joseph K Aicher
  <https://github.com/jaicher>`_.
- Use ``dask_array_type`` instead of ``dask_array.Array`` for type
  checking. (:issue:`3779`, :pull:`3787`)
  By `Justus Magin <https://github.com/keewis>`_.
- :py:func:`concat` can now handle coordinate variables only present in one of
  the objects to be concatenated when ``coords="different"``.
  By `Deepak Cherian <https://github.com/dcherian>`_.
- xarray now respects the over, under and bad colors if set on a provided colormap.
  (:issue:`3590`, :pull:`3601`)
  By `johnomotani <https://github.com/johnomotani>`_.
- ``coarsen`` and ``rolling`` now respect ``xr.set_options(keep_attrs=True)``
  to preserve attributes. :py:meth:`Dataset.coarsen` accepts a keyword
  argument ``keep_attrs`` to change this setting. (:issue:`3376`,
  :pull:`3801`) By `Andrew Thomas <https://github.com/amcnicho>`_.
- Delete associated indexes when deleting coordinate variables. (:issue:`3746`).
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Fix :py:meth:`Dataset.to_zarr` when using ``append_dim`` and ``group``
  simultaneously. (:issue:`3170`). By `Matthias Meyer <https://github.com/niowniow>`_.
- Fix html repr on :py:class:`Dataset` with non-string keys (:pull:`3807`).
  By `Maximilian Roos <https://github.com/max-sixty>`_.

Documentation
~~~~~~~~~~~~~

- Fix documentation of :py:class:`DataArray` removing the deprecated mention
  that when omitted, `dims` are inferred from a `coords`-dict. (:pull:`3821`)
  By `Sander van Rijn <https://github.com/sjvrijn>`_.
- Improve the :py:func:`where` docstring.
  By `Maximilian Roos <https://github.com/max-sixty>`_
- Update the installation instructions: only explicitly list recommended dependencies
  (:issue:`3756`).
  By `Mathias Hauser <https://github.com/mathause>`_.

Internal Changes
~~~~~~~~~~~~~~~~

- Remove the internal ``import_seaborn`` function which handled the deprecation of
  the ``seaborn.apionly`` entry point (:issue:`3747`).
  By `Mathias Hauser <https://github.com/mathause>`_.
- Don't test pint integration in combination with datetime objects. (:issue:`3778`, :pull:`3788`)
  By `Justus Magin <https://github.com/keewis>`_.
- Change test_open_mfdataset_list_attr to only run with dask installed
  (:issue:`3777`, :pull:`3780`).
  By `Bruno Pagani <https://github.com/ArchangeGabriel>`_.
- Preserve the ability to index with ``method="nearest"`` with a
  :py:class:`CFTimeIndex` with pandas versions greater than 1.0.1
  (:issue:`3751`). By `Spencer Clark <https://github.com/spencerkclark>`_.
- Greater flexibility and improved test coverage of subtracting various types
  of objects from a :py:class:`CFTimeIndex`. By `Spencer Clark
  <https://github.com/spencerkclark>`_.
- Update Azure CI MacOS image, given pending removal.
  By `Maximilian Roos <https://github.com/max-sixty>`_
- Remove xfails for scipy 1.0.1 for tests that append to netCDF files (:pull:`3805`).
  By `Mathias Hauser <https://github.com/mathause>`_.
- Remove conversion to ``pandas.Panel``, given its removal in pandas
  in favor of xarray's objects.
  By `Maximilian Roos <https://github.com/max-sixty>`_

.. _whats-new.0.15.0:


v0.15.0 (30 Jan 2020)
---------------------

This release brings many improvements to xarray's documentation: our examples are now binderized notebooks (`click here <https://mybinder.org/v2/gh/pydata/xarray/main?urlpath=lab/tree/doc/examples/weather-data.ipynb>`_)
and we have new example notebooks from our SciPy 2019 sprint (many thanks to our contributors!).

This release also features many API improvements such as a new
:py:class:`~core.accessor_dt.TimedeltaAccessor` and support for :py:class:`CFTimeIndex` in
:py:meth:`~DataArray.interpolate_na`); as well as many bug fixes.

Breaking changes
~~~~~~~~~~~~~~~~
- Bumped minimum tested versions for dependencies:

  - numpy 1.15
  - pandas 0.25
  - dask 2.2
  - distributed 2.2
  - scipy 1.3

- Remove ``compat`` and ``encoding`` kwargs from ``DataArray``, which
  have been deprecated since 0.12. (:pull:`3650`).
  Instead, specify the ``encoding`` kwarg when writing to disk or set
  the :py:attr:`DataArray.encoding` attribute directly.
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- :py:func:`xarray.dot`, :py:meth:`DataArray.dot`, and the ``@`` operator now
  use ``align="inner"`` (except when ``xarray.set_options(arithmetic_join="exact")``;
  :issue:`3694`) by `Mathias Hauser <https://github.com/mathause>`_.

New Features
~~~~~~~~~~~~
- Implement :py:meth:`DataArray.pad` and :py:meth:`Dataset.pad`. (:issue:`2605`, :pull:`3596`).
  By `Mark Boer <https://github.com/mark-boer>`_.
- :py:meth:`DataArray.sel` and :py:meth:`Dataset.sel` now support :py:class:`pandas.CategoricalIndex`. (:issue:`3669`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- Support using an existing, opened h5netcdf ``File`` with
  :py:class:`~xarray.backends.H5NetCDFStore`. This permits creating an
  :py:class:`~xarray.Dataset` from a h5netcdf ``File`` that has been opened
  using other means (:issue:`3618`).
  By `Kai Mühlbauer <https://github.com/kmuehlbauer>`_.
- Implement ``median`` and ``nanmedian`` for dask arrays. This works by rechunking
  to a single chunk along all reduction axes. (:issue:`2999`).
  By `Deepak Cherian <https://github.com/dcherian>`_.
- :py:func:`~xarray.concat` now preserves attributes from the first Variable.
  (:issue:`2575`, :issue:`2060`, :issue:`1614`)
  By `Deepak Cherian <https://github.com/dcherian>`_.
- :py:meth:`Dataset.quantile`, :py:meth:`DataArray.quantile` and ``GroupBy.quantile``
  now work with dask Variables.
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Added the ``count`` reduction method to both :py:class:`~core.rolling.DatasetCoarsen`
  and :py:class:`~core.rolling.DataArrayCoarsen` objects. (:pull:`3500`)
  By `Deepak Cherian <https://github.com/dcherian>`_
- Add ``meta`` kwarg to :py:func:`~xarray.apply_ufunc`;
  this is passed on to :py:func:`dask.array.blockwise`. (:pull:`3660`)
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Add ``attrs_file`` option in :py:func:`~xarray.open_mfdataset` to choose the
  source file for global attributes in a multi-file dataset (:issue:`2382`,
  :pull:`3498`). By `Julien Seguinot <https://github.com/juseg>`_.
- :py:meth:`Dataset.swap_dims` and :py:meth:`DataArray.swap_dims`
  now allow swapping to dimension names that don't exist yet. (:pull:`3636`)
  By `Justus Magin <https://github.com/keewis>`_.
- Extend :py:class:`~core.accessor_dt.DatetimeAccessor` properties
  and support ``.dt`` accessor for timedeltas
  via :py:class:`~core.accessor_dt.TimedeltaAccessor` (:pull:`3612`)
  By `Anderson Banihirwe <https://github.com/andersy005>`_.
- Improvements to interpolating along time axes (:issue:`3641`, :pull:`3631`).
  By `David Huard <https://github.com/huard>`_.

  - Support :py:class:`CFTimeIndex` in :py:meth:`DataArray.interpolate_na`
  - define 1970-01-01 as the default offset for the interpolation index for both
    :py:class:`pandas.DatetimeIndex` and :py:class:`CFTimeIndex`,
  - use microseconds in the conversion from timedelta objects to floats to avoid
    overflow errors.

Bug fixes
~~~~~~~~~
- Applying a user-defined function that adds new dimensions using :py:func:`apply_ufunc`
  and ``vectorize=True`` now works with ``dask > 2.0``. (:issue:`3574`, :pull:`3660`).
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Fix :py:meth:`~xarray.combine_by_coords` to allow for combining incomplete
  hypercubes of Datasets (:issue:`3648`).  By `Ian Bolliger
  <https://github.com/bolliger32>`_.
- Fix :py:func:`~xarray.combine_by_coords` when combining cftime coordinates
  which span long time intervals (:issue:`3535`).  By `Spencer Clark
  <https://github.com/spencerkclark>`_.
- Fix plotting with transposed 2D non-dimensional coordinates. (:issue:`3138`, :pull:`3441`)
  By `Deepak Cherian <https://github.com/dcherian>`_.
- :py:meth:`plot.FacetGrid.set_titles` can now replace existing row titles of a
  :py:class:`~xarray.plot.FacetGrid` plot. In addition :py:class:`~xarray.plot.FacetGrid` gained
  two new attributes: :py:attr:`~xarray.plot.FacetGrid.col_labels` and
  :py:attr:`~xarray.plot.FacetGrid.row_labels` contain :py:class:`matplotlib.text.Text` handles for both column and
  row labels. These can be used to manually change the labels.
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Fix issue with Dask-backed datasets raising a ``KeyError`` on some computations involving :py:func:`map_blocks` (:pull:`3598`).
  By `Tom Augspurger <https://github.com/TomAugspurger>`_.
- Ensure :py:meth:`Dataset.quantile`, :py:meth:`DataArray.quantile` issue the correct error
  when ``q`` is out of bounds (:issue:`3634`) by `Mathias Hauser <https://github.com/mathause>`_.
- Fix regression in xarray 0.14.1 that prevented encoding times with certain
  ``dtype``, ``_FillValue``, and ``missing_value`` encodings (:issue:`3624`).
  By `Spencer Clark <https://github.com/spencerkclark>`_
- Raise an error when trying to use :py:meth:`Dataset.rename_dims` to
  rename to an existing name (:issue:`3438`, :pull:`3645`)
  By `Justus Magin <https://github.com/keewis>`_.
- :py:meth:`Dataset.rename`, :py:meth:`DataArray.rename` now check for conflicts with
  MultiIndex level names.
- :py:meth:`Dataset.merge` no longer fails when passed a :py:class:`DataArray` instead of a :py:class:`Dataset`.
  By `Tom Nicholas <https://github.com/TomNicholas>`_.
- Fix a regression in :py:meth:`Dataset.drop`: allow passing any
  iterable when dropping variables (:issue:`3552`, :pull:`3693`)
  By `Justus Magin <https://github.com/keewis>`_.
- Fixed errors emitted by ``mypy --strict`` in modules that import xarray.
  (:issue:`3695`) by `Guido Imperiale <https://github.com/crusaderky>`_.
- Allow plotting of binned coordinates on the y axis in :py:meth:`plot.line`
  and :py:meth:`plot.step` plots (:issue:`3571`,
  :pull:`3685`) by `Julien Seguinot <https://github.com/juseg>`_.
- setuptools is now marked as a dependency of xarray
  (:pull:`3628`) by `Richard Höchenberger <https://github.com/hoechenberger>`_.

Documentation
~~~~~~~~~~~~~
- Switch doc examples to use `nbsphinx <https://nbsphinx.readthedocs.io>`_ and replace
  ``sphinx_gallery`` scripts with Jupyter notebooks. (:pull:`3105`, :pull:`3106`, :pull:`3121`)
  By `Ryan Abernathey <https://github.com/rabernat>`_.
- Added :doc:`example notebook <examples/ROMS_ocean_model>` demonstrating use of xarray with
  Regional Ocean Modeling System (ROMS) ocean hydrodynamic model output. (:pull:`3116`)
  By `Robert Hetland <https://github.com/hetland>`_.
- Added :doc:`example notebook <examples/ERA5-GRIB-example>` demonstrating the visualization of
  ERA5 GRIB data. (:pull:`3199`)
  By `Zach Bruick <https://github.com/zbruick>`_ and
  `Stephan Siemen <https://github.com/StephanSiemen>`_.
- Added examples for :py:meth:`DataArray.quantile`, :py:meth:`Dataset.quantile` and
  ``GroupBy.quantile``. (:pull:`3576`)
  By `Justus Magin <https://github.com/keewis>`_.
- Add new :doc:`example notebook <examples/apply_ufunc_vectorize_1d>` example notebook demonstrating
  vectorization of a 1D function using :py:func:`apply_ufunc` , dask and numba.
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Added example for :py:func:`~xarray.map_blocks`. (:pull:`3667`)
  By `Riley X. Brady <https://github.com/bradyrx>`_.

Internal Changes
~~~~~~~~~~~~~~~~
- Make sure dask names change when rechunking by different chunk sizes. Conversely, make sure they
  stay the same when rechunking by the same chunk size. (:issue:`3350`)
  By `Deepak Cherian <https://github.com/dcherian>`_.
- 2x to 5x speed boost (on small arrays) for :py:meth:`Dataset.isel`,
  :py:meth:`DataArray.isel`, and :py:meth:`DataArray.__getitem__` when indexing by int,
  slice, list of int, scalar ndarray, or 1-dimensional ndarray.
  (:pull:`3533`) by `Guido Imperiale <https://github.com/crusaderky>`_.
- Removed internal method ``Dataset._from_vars_and_coord_names``,
  which was dominated by ``Dataset._construct_direct``. (:pull:`3565`)
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- Replaced versioneer with setuptools-scm. Moved contents of setup.py to setup.cfg.
  Removed pytest-runner from setup.py, as per deprecation notice on the pytest-runner
  project. (:pull:`3714`) by `Guido Imperiale <https://github.com/crusaderky>`_.
- Use of isort is now enforced by CI.
  (:pull:`3721`) by `Guido Imperiale <https://github.com/crusaderky>`_


.. _whats-new.0.14.1:

v0.14.1 (19 Nov 2019)
---------------------

Breaking changes
~~~~~~~~~~~~~~~~

- Broken compatibility with ``cftime < 1.0.3`` . By `Deepak Cherian <https://github.com/dcherian>`_.

  .. warning::

    cftime version 1.0.4 is broken
    (`cftime/126 <https://github.com/Unidata/cftime/issues/126>`_);
    please use version 1.0.4.2 instead.

- All leftover support for dates from non-standard calendars through ``netcdftime``, the
  module included in versions of netCDF4 prior to 1.4 that eventually became the
  `cftime <https://github.com/Unidata/cftime/>`_ package, has been removed in favor of relying solely on
  the standalone ``cftime`` package (:pull:`3450`).
  By `Spencer Clark <https://github.com/spencerkclark>`_.

New Features
~~~~~~~~~~~~
- Added the ``sparse`` option to :py:meth:`~xarray.DataArray.unstack`,
  :py:meth:`~xarray.Dataset.unstack`, :py:meth:`~xarray.DataArray.reindex`,
  :py:meth:`~xarray.Dataset.reindex` (:issue:`3518`).
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- Added the ``fill_value`` option to :py:meth:`DataArray.unstack` and
  :py:meth:`Dataset.unstack` (:issue:`3518`, :pull:`3541`).
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- Added the ``max_gap`` kwarg to :py:meth:`~xarray.DataArray.interpolate_na` and
  :py:meth:`~xarray.Dataset.interpolate_na`. This controls the maximum size of the data
  gap that will be filled by interpolation. By `Deepak Cherian <https://github.com/dcherian>`_.
- Added :py:meth:`Dataset.drop_sel` & :py:meth:`DataArray.drop_sel` for dropping labels.
  :py:meth:`Dataset.drop_vars` & :py:meth:`DataArray.drop_vars` have been added for
  dropping variables (including coordinates). The existing :py:meth:`Dataset.drop` &
  :py:meth:`DataArray.drop` methods remain as a backward compatible
  option for dropping either labels or variables, but using the more specific methods is encouraged.
  (:pull:`3475`)
  By `Maximilian Roos <https://github.com/max-sixty>`_
- Added :py:meth:`Dataset.map` & ``GroupBy.map`` & ``Resample.map`` for
  mapping / applying a function over each item in the collection, reflecting the widely used
  and least surprising name for this operation.
  The existing ``apply`` methods remain for backward compatibility, though using the ``map``
  methods is encouraged.
  (:pull:`3459`)
  By `Maximilian Roos <https://github.com/max-sixty>`_
- :py:meth:`Dataset.transpose` and :py:meth:`DataArray.transpose` now support an ellipsis (``...``)
  to represent all 'other' dimensions. For example, to move one dimension to the front,
  use ``.transpose('x', ...)``. (:pull:`3421`)
  By `Maximilian Roos <https://github.com/max-sixty>`_
- Changed ``xr.ALL_DIMS`` to equal python's ``Ellipsis`` (``...``), and changed internal usages to use
  ``...`` directly. As before, you can use this to instruct a ``groupby`` operation
  to reduce over all dimensions. While we have no plans to remove ``xr.ALL_DIMS``, we suggest
  using ``...``. (:pull:`3418`)
  By `Maximilian Roos <https://github.com/max-sixty>`_
- :py:func:`xarray.dot`, and :py:meth:`DataArray.dot` now support the
  ``dims=...`` option to sum over the union of dimensions of all input arrays
  (:issue:`3423`) by `Mathias Hauser <https://github.com/mathause>`_.
- Added new ``Dataset._repr_html_`` and ``DataArray._repr_html_`` to improve
  representation of objects in Jupyter. By default this feature is turned off
  for now. Enable it with ``xarray.set_options(display_style="html")``.
  (:pull:`3425`) by `Benoit Bovy <https://github.com/benbovy>`_ and
  `Julia Signell <https://github.com/jsignell>`_.
- Implement `dask deterministic hashing
  <https://docs.dask.org/en/latest/custom-collections.html#deterministic-hashing>`_
  for xarray objects. Note that xarray objects with a dask.array backend already used
  deterministic hashing in previous releases; this change implements it when whole
  xarray objects are embedded in a dask graph, e.g. when :py:meth:`DataArray.map_blocks` is
  invoked. (:issue:`3378`, :pull:`3446`, :pull:`3515`)
  By `Deepak Cherian <https://github.com/dcherian>`_ and
  `Guido Imperiale <https://github.com/crusaderky>`_.
- Add the documented-but-missing :py:meth:`~core.groupby.DatasetGroupBy.quantile`.
- xarray now respects the ``DataArray.encoding["coordinates"]`` attribute when writing to disk.
  See :ref:`io.coordinates` for more. (:issue:`3351`, :pull:`3487`)
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Add the documented-but-missing :py:meth:`~core.groupby.DatasetGroupBy.quantile`.
  (:issue:`3525`, :pull:`3527`). By `Justus Magin <https://github.com/keewis>`_.

Bug fixes
~~~~~~~~~
- Ensure an index of type ``CFTimeIndex`` is not converted to a ``DatetimeIndex`` when
  calling :py:meth:`Dataset.rename`, :py:meth:`Dataset.rename_dims` and :py:meth:`Dataset.rename_vars`.
  By `Mathias Hauser <https://github.com/mathause>`_. (:issue:`3522`).
- Fix a bug in :py:meth:`DataArray.set_index` in case that an existing dimension becomes a level
  variable of MultiIndex. (:pull:`3520`). By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- Harmonize ``_FillValue``, ``missing_value`` during encoding and decoding steps. (:pull:`3502`)
  By `Anderson Banihirwe <https://github.com/andersy005>`_.
- Fix regression introduced in v0.14.0 that would cause a crash if dask is installed
  but cloudpickle isn't (:issue:`3401`) by `Rhys Doyle <https://github.com/rdoyle45>`_
- Fix grouping over variables with NaNs. (:issue:`2383`, :pull:`3406`).
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Make alignment and concatenation significantly more efficient by using dask names to compare dask
  objects prior to comparing values after computation. This change makes it more convenient to carry
  around large non-dimensional coordinate variables backed by dask arrays. Existing workarounds involving
  ``reset_coords(drop=True)`` should now be unnecessary in most cases.
  (:issue:`3068`, :issue:`3311`, :issue:`3454`, :pull:`3453`).
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Add support for cftime>=1.0.4. By `Anderson Banihirwe <https://github.com/andersy005>`_.
- Rolling reduction operations no longer compute dask arrays by default. (:issue:`3161`).
  In addition, the ``allow_lazy`` kwarg to ``reduce`` is deprecated.
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Fix ``GroupBy.reduce`` when reducing over multiple dimensions.
  (:issue:`3402`). By `Deepak Cherian <https://github.com/dcherian>`_
- Allow appending datetime and bool data variables to zarr stores.
  (:issue:`3480`). By `Akihiro Matsukawa <https://github.com/amatsukawa>`_.
- Add support for numpy >=1.18 (); bugfix mean() on datetime64 arrays on dask backend
  (:issue:`3409`, :pull:`3537`). By `Guido Imperiale <https://github.com/crusaderky>`_.
- Add support for pandas >=0.26 (:issue:`3440`).
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Add support for pseudonetcdf >=3.1 (:pull:`3485`).
  By `Barron Henderson <https://github.com/barronh>`_.

Documentation
~~~~~~~~~~~~~
- Fix leap year condition in `monthly means example <http://xarray.pydata.org/en/stable/examples/monthly-means.html>`_.
  By `Mickaël Lalande <https://github.com/mickaellalande>`_.
- Fix the documentation of :py:meth:`DataArray.resample` and
  :py:meth:`Dataset.resample`,  explicitly stating that a
  datetime-like dimension is required. (:pull:`3400`)
  By `Justus Magin <https://github.com/keewis>`_.
- Update the :ref:`terminology` page to address multidimensional coordinates. (:pull:`3410`)
  By `Jon Thielen <https://github.com/jthielen>`_.
- Fix the documentation of :py:meth:`Dataset.integrate` and
  :py:meth:`DataArray.integrate` and add an example to
  :py:meth:`Dataset.integrate`. (:pull:`3469`)
  By `Justus Magin <https://github.com/keewis>`_.

Internal Changes
~~~~~~~~~~~~~~~~

- Added integration tests against `pint <https://pint.readthedocs.io/>`_.
  (:pull:`3238`, :pull:`3447`, :pull:`3493`, :pull:`3508`)
  by `Justus Magin <https://github.com/keewis>`_.

  .. note::

    At the moment of writing, these tests *as well as the ability to use pint in general*
    require `a highly experimental version of pint
    <https://github.com/andrewgsavage/pint/pull/6>`_ (install with
    ``pip install git+https://github.com/andrewgsavage/pint.git@refs/pull/6/head)``.
    Even with it, interaction with non-numpy array libraries, e.g. dask or sparse, is broken.

- Use Python 3.6 idioms throughout the codebase. (:pull:`3419`)
  By `Maximilian Roos <https://github.com/max-sixty>`_

- Run basic CI tests on Python 3.8. (:pull:`3477`)
  By `Maximilian Roos <https://github.com/max-sixty>`_

- Enable type checking on default sentinel values (:pull:`3472`)
  By `Maximilian Roos <https://github.com/max-sixty>`_

- Add ``Variable._replace`` for simpler replacing of a subset of attributes (:pull:`3472`)
  By `Maximilian Roos <https://github.com/max-sixty>`_

.. _whats-new.0.14.0:

v0.14.0 (14 Oct 2019)
---------------------

Breaking changes
~~~~~~~~~~~~~~~~
- This release introduces a rolling policy for minimum dependency versions:
  :ref:`mindeps_policy`.

  Several minimum versions have been increased:

  ============ ================== ====
  Package      Old                New
  ============ ================== ====
  Python       3.5.3              3.6
  numpy        1.12               1.14
  pandas       0.19.2             0.24
  dask         0.16 (tested: 2.4) 1.2
  bottleneck   1.1 (tested: 1.2)  1.2
  matplotlib   1.5 (tested: 3.1)  3.1
  ============ ================== ====

  Obsolete patch versions (x.y.Z) are not tested anymore.
  The oldest supported versions of all optional dependencies are now covered by
  automated tests (before, only the very latest versions were tested).

  (:issue:`3222`, :issue:`3293`, :issue:`3340`, :issue:`3346`, :issue:`3358`).
  By `Guido Imperiale <https://github.com/crusaderky>`_.

- Dropped the ``drop=False`` optional parameter from :py:meth:`Variable.isel`.
  It was unused and doesn't make sense for a Variable. (:pull:`3375`).
  By `Guido Imperiale <https://github.com/crusaderky>`_.

- Remove internal usage of :py:class:`collections.OrderedDict`. After dropping support for
  Python <=3.5, most uses of ``OrderedDict`` in xarray were no longer necessary. We
  have removed the internal use of the ``OrderedDict`` in favor of Python's builtin
  ``dict`` object which is now ordered itself. This change will be most obvious when
  interacting with the ``attrs`` property on Dataset and DataArray objects.
  (:issue:`3380`, :pull:`3389`). By `Joe Hamman <https://github.com/jhamman>`_.

New functions/methods
~~~~~~~~~~~~~~~~~~~~~

- Added :py:func:`~xarray.map_blocks`, modeled after :py:func:`dask.array.map_blocks`.
  Also added :py:meth:`Dataset.unify_chunks`, :py:meth:`DataArray.unify_chunks` and
  :py:meth:`testing.assert_chunks_equal`. (:pull:`3276`).
  By `Deepak Cherian <https://github.com/dcherian>`_ and
  `Guido Imperiale <https://github.com/crusaderky>`_.

Enhancements
~~~~~~~~~~~~

- ``core.groupby.GroupBy`` enhancements. By `Deepak Cherian <https://github.com/dcherian>`_.

  - Added a repr (:pull:`3344`). Example::

      >>> da.groupby("time.season")
      DataArrayGroupBy, grouped over 'season'
      4 groups with labels 'DJF', 'JJA', 'MAM', 'SON'

  - Added a ``GroupBy.dims`` property that mirrors the dimensions
    of each group (:issue:`3344`).

- Speed up :py:meth:`Dataset.isel` up to 33% and :py:meth:`DataArray.isel` up to 25% for small
  arrays (:issue:`2799`, :pull:`3375`). By
  `Guido Imperiale <https://github.com/crusaderky>`_.

Bug fixes
~~~~~~~~~
- Reintroduce support for :mod:`weakref` (broken in v0.13.0). Support has been
  reinstated for :py:class:`~xarray.DataArray` and :py:class:`~xarray.Dataset` objects only.
  Internal xarray objects remain unaddressable by weakref in order to save memory
  (:issue:`3317`). By `Guido Imperiale <https://github.com/crusaderky>`_.
- Line plots with the ``x`` or ``y`` argument set to a 1D non-dimensional coord
  now plot the correct data for 2D DataArrays
  (:issue:`3334`). By `Tom Nicholas <https://github.com/TomNicholas>`_.
- Make :py:func:`~xarray.concat` more robust when merging variables present in some datasets but
  not others (:issue:`508`). By `Deepak Cherian <https://github.com/dcherian>`_.
- The default behaviour of reducing across all dimensions for
  :py:class:`~xarray.core.groupby.DataArrayGroupBy` objects has now been properly removed
  as was done for :py:class:`~xarray.core.groupby.DatasetGroupBy` in 0.13.0 (:issue:`3337`).
  Use ``xarray.ALL_DIMS`` if you need to replicate previous behaviour.
  Also raise nicer error message when no groups are created (:issue:`1764`).
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Fix error in concatenating unlabeled dimensions (:pull:`3362`).
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Warn if the ``dim`` kwarg is passed to rolling operations. This is redundant since a dimension is
  specified when the :py:class:`~core.rolling.DatasetRolling` or :py:class:`~core.rolling.DataArrayRolling` object is created.
  (:pull:`3362`). By `Deepak Cherian <https://github.com/dcherian>`_.

Documentation
~~~~~~~~~~~~~

- Created a glossary of important xarray terms (:issue:`2410`, :pull:`3352`).
  By `Gregory Gundersen <https://github.com/gwgundersen>`_.
- Created a "How do I..." section (:ref:`howdoi`) for solutions to common questions. (:pull:`3357`).
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Add examples for :py:meth:`Dataset.swap_dims` and :py:meth:`DataArray.swap_dims`
  (pull:`3331`, pull:`3331`). By `Justus Magin <https://github.com/keewis>`_.
- Add examples for :py:meth:`align`, :py:meth:`merge`, :py:meth:`combine_by_coords`,
  :py:meth:`full_like`, :py:meth:`zeros_like`, :py:meth:`ones_like`, :py:meth:`Dataset.pipe`,
  :py:meth:`Dataset.assign`, :py:meth:`Dataset.reindex`, :py:meth:`Dataset.fillna` (:pull:`3328`).
  By `Anderson Banihirwe <https://github.com/andersy005>`_.
- Fixed documentation to clean up an unwanted file created in ``ipython`` example
  (:pull:`3353`). By `Gregory Gundersen <https://github.com/gwgundersen>`_.

.. _whats-new.0.13.0:

v0.13.0 (17 Sep 2019)
---------------------

This release includes many exciting changes: wrapping of
`NEP18 <https://www.numpy.org/neps/nep-0018-array-function-protocol.html>`_ compliant
numpy-like arrays; new :py:meth:`~Dataset.plot.scatter` plotting method that can scatter
two ``DataArrays`` in a ``Dataset`` against each other; support for converting pandas
DataFrames to xarray objects that wrap ``pydata/sparse``; and more!

Breaking changes
~~~~~~~~~~~~~~~~

- This release increases the minimum required Python version from 3.5.0 to 3.5.3
  (:issue:`3089`). By `Guido Imperiale <https://github.com/crusaderky>`_.
- The ``isel_points`` and ``sel_points`` methods are removed, having been deprecated
  since v0.10.0. These are redundant with the ``isel`` / ``sel`` methods.
  See :ref:`vectorized_indexing` for the details
  By `Maximilian Roos <https://github.com/max-sixty>`_
- The ``inplace`` kwarg for public methods now raises an error, having been deprecated
  since v0.11.0.
  By `Maximilian Roos <https://github.com/max-sixty>`_
- :py:func:`~xarray.concat` now requires the ``dim`` argument. Its ``indexers``, ``mode``
  and ``concat_over`` kwargs have now been removed.
  By `Deepak Cherian <https://github.com/dcherian>`_
- Passing a list of colors in ``cmap`` will now raise an error, having been deprecated since
  v0.6.1.
- Most xarray objects now define ``__slots__``. This reduces overall RAM usage by ~22%
  (not counting the underlying numpy buffers); on CPython 3.7/x64, a trivial DataArray
  has gone down from 1.9kB to 1.5kB.

  Caveats:

  - Pickle streams produced by older versions of xarray can't be loaded using this
    release, and vice versa.
  - Any user code that was accessing the ``__dict__`` attribute of
    xarray objects will break. The best practice to attach custom metadata to xarray
    objects is to use the ``attrs`` dictionary.
  - Any user code that defines custom subclasses of xarray classes must now explicitly
    define ``__slots__`` itself. Subclasses that don't add any attributes must state so
    by defining ``__slots__ = ()`` right after the class header.
    Omitting ``__slots__`` will now cause a ``FutureWarning`` to be logged, and will raise an
    error in a later release.

  (:issue:`3250`) by `Guido Imperiale <https://github.com/crusaderky>`_.
- The default dimension for :py:meth:`Dataset.groupby`, :py:meth:`Dataset.resample`,
  :py:meth:`DataArray.groupby` and :py:meth:`DataArray.resample` reductions is now the
  grouping or resampling dimension.
- :py:meth:`DataArray.to_dataset` requires ``name`` to be passed as a kwarg (previously ambiguous
  positional arguments were deprecated)
- Reindexing with variables of a different dimension now raise an error (previously deprecated)
- ``xarray.broadcast_array`` is removed (previously deprecated in favor of
  :py:func:`~xarray.broadcast`)
- ``Variable.expand_dims`` is removed (previously deprecated in favor of
  :py:meth:`Variable.set_dims`)

New functions/methods
~~~~~~~~~~~~~~~~~~~~~

- xarray can now wrap around any
  `NEP18 <https://www.numpy.org/neps/nep-0018-array-function-protocol.html>`_ compliant
  numpy-like library (important: read notes about ``NUMPY_EXPERIMENTAL_ARRAY_FUNCTION`` in
  the above link). Added explicit test coverage for
  `sparse <https://github.com/pydata/sparse>`_. (:issue:`3117`, :issue:`3202`).
  This requires `sparse>=0.8.0`. By `Nezar Abdennur <https://github.com/nvictus>`_
  and `Guido Imperiale <https://github.com/crusaderky>`_.

- :py:meth:`~Dataset.from_dataframe` and :py:meth:`~DataArray.from_series` now
  support ``sparse=True`` for converting pandas objects into xarray objects
  wrapping sparse arrays. This is particularly useful with sparsely populated
  hierarchical indexes. (:issue:`3206`)
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- The xarray package is now discoverable by mypy (although typing hints coverage is not
  complete yet). mypy type checking is now enforced by CI. Libraries that depend on
  xarray and use mypy can now remove from their setup.cfg the lines::

    [mypy-xarray]
    ignore_missing_imports = True

  (:issue:`2877`, :issue:`3088`, :issue:`3090`, :issue:`3112`, :issue:`3117`,
  :issue:`3207`)
  By `Guido Imperiale <https://github.com/crusaderky>`_
  and `Maximilian Roos <https://github.com/max-sixty>`_.

- Added :py:meth:`DataArray.broadcast_like` and :py:meth:`Dataset.broadcast_like`.
  By `Deepak Cherian <https://github.com/dcherian>`_ and `David Mertz
  <https://github.com/DavidMertz>`_.

- Dataset plotting API for visualizing dependencies between two DataArrays!
  Currently only :py:meth:`Dataset.plot.scatter` is implemented.
  By `Yohai Bar Sinai <https://github.com/yohai>`_ and `Deepak Cherian <https://github.com/dcherian>`_

- Added :py:meth:`DataArray.head`, :py:meth:`DataArray.tail` and :py:meth:`DataArray.thin`;
  as well as :py:meth:`Dataset.head`, :py:meth:`Dataset.tail` and :py:meth:`Dataset.thin` methods.
  (:issue:`319`) By `Gerardo Rivera <https://github.com/dangomelon>`_.

Enhancements
~~~~~~~~~~~~

- Multiple enhancements to :py:func:`~xarray.concat` and :py:func:`~xarray.open_mfdataset`.
  By `Deepak Cherian <https://github.com/dcherian>`_

  - Added ``compat='override'``. When merging, this option picks the variable from the first dataset
    and skips all comparisons.

  - Added ``join='override'``. When aligning, this only checks that index sizes are equal among objects
    and skips checking indexes for equality.

  - :py:func:`~xarray.concat` and :py:func:`~xarray.open_mfdataset` now support the ``join`` kwarg.
    It is passed down to :py:func:`~xarray.align`.

  - :py:func:`~xarray.concat` now calls :py:func:`~xarray.merge` on variables that are not concatenated
    (i.e. variables without ``concat_dim`` when ``data_vars`` or ``coords`` are ``"minimal"``).
    :py:func:`~xarray.concat` passes its new ``compat`` kwarg down to :py:func:`~xarray.merge`.
    (:issue:`2064`)

  Users can avoid a common bottleneck when using :py:func:`~xarray.open_mfdataset` on a large number of
  files with variables that are known to be aligned and some of which need not be concatenated.
  Slow equality comparisons can now be avoided, for e.g.::

    data = xr.open_mfdataset(files, concat_dim='time', data_vars='minimal',
                             coords='minimal', compat='override', join='override')

- In :py:meth:`~xarray.Dataset.to_zarr`, passing ``mode`` is not mandatory if
  ``append_dim`` is set, as it will automatically be set to ``'a'`` internally.
  By `David Brochart <https://github.com/davidbrochart>`_.

- Added the ability to initialize an empty or full DataArray
  with a single value. (:issue:`277`)
  By `Gerardo Rivera <https://github.com/dangomelon>`_.

- :py:func:`~xarray.Dataset.to_netcdf()` now supports the ``invalid_netcdf`` kwarg when used
  with ``engine="h5netcdf"``. It is passed to ``h5netcdf.File``.
  By `Ulrich Herter <https://github.com/ulijh>`_.

- ``xarray.Dataset.drop`` now supports keyword arguments; dropping index
  labels by using both ``dim`` and ``labels`` or using a
  :py:class:`~core.coordinates.DataArrayCoordinates` object are deprecated (:issue:`2910`).
  By `Gregory Gundersen <https://github.com/gwgundersen>`_.

- Added examples of :py:meth:`Dataset.set_index` and
  :py:meth:`DataArray.set_index`, as well are more specific error messages
  when the user passes invalid arguments (:issue:`3176`).
  By `Gregory Gundersen <https://github.com/gwgundersen>`_.

- :py:meth:`Dataset.filter_by_attrs` now filters the coordinates as well as the variables.
  By `Spencer Jones <https://github.com/cspencerjones>`_.

Bug fixes
~~~~~~~~~

- Improve "missing dimensions" error message for :py:func:`~xarray.apply_ufunc`
  (:issue:`2078`).
  By `Rick Russotto <https://github.com/rdrussotto>`_.
- :py:meth:`~xarray.DataArray.assign_coords` now supports dictionary arguments
  (:issue:`3231`).
  By `Gregory Gundersen <https://github.com/gwgundersen>`_.
- Fix regression introduced in v0.12.2 where ``copy(deep=True)`` would convert
  unicode indices to dtype=object (:issue:`3094`).
  By `Guido Imperiale <https://github.com/crusaderky>`_.
- Improved error handling and documentation for `.expand_dims()`
  read-only view.
- Fix tests for big-endian systems (:issue:`3125`).
  By `Graham Inggs <https://github.com/ginggs>`_.
- XFAIL several tests which are expected to fail on ARM systems
  due to a ``datetime`` issue in NumPy (:issue:`2334`).
  By `Graham Inggs <https://github.com/ginggs>`_.
- Fix KeyError that arises when using .sel method with float values
  different from coords float type (:issue:`3137`).
  By `Hasan Ahmad <https://github.com/HasanAhmadQ7>`_.
- Fixed bug in ``combine_by_coords()`` causing a `ValueError` if the input had
  an unused dimension with coordinates which were not monotonic (:issue:`3150`).
  By `Tom Nicholas <https://github.com/TomNicholas>`_.
- Fixed crash when applying ``distributed.Client.compute()`` to a DataArray
  (:issue:`3171`). By `Guido Imperiale <https://github.com/crusaderky>`_.
- Better error message when using groupby on an empty DataArray (:issue:`3037`).
  By `Hasan Ahmad <https://github.com/HasanAhmadQ7>`_.
- Fix error that arises when using open_mfdataset on a series of netcdf files
  having differing values for a variable attribute of type list. (:issue:`3034`)
  By `Hasan Ahmad <https://github.com/HasanAhmadQ7>`_.
- Prevent :py:meth:`~xarray.DataArray.argmax` and :py:meth:`~xarray.DataArray.argmin` from calling
  dask compute (:issue:`3237`). By `Ulrich Herter <https://github.com/ulijh>`_.
- Plots in 2 dimensions (pcolormesh, contour) now allow to specify levels as numpy
  array (:issue:`3284`). By `Mathias Hauser <https://github.com/mathause>`_.
- Fixed bug in :meth:`DataArray.quantile` failing to keep attributes when
  `keep_attrs` was True (:issue:`3304`). By `David Huard <https://github.com/huard>`_.

Documentation
~~~~~~~~~~~~~

- Created a `PR checklist <https://xarray.pydata.org/en/stable/contributing.html/contributing.html#pr-checklist>`_
  as a quick reference for tasks before creating a new PR
  or pushing new commits.
  By `Gregory Gundersen <https://github.com/gwgundersen>`_.

- Fixed documentation to clean up unwanted files created in ``ipython`` examples
  (:issue:`3227`).
  By `Gregory Gundersen <https://github.com/gwgundersen>`_.

.. _whats-new.0.12.3:

v0.12.3 (10 July 2019)
----------------------

New functions/methods
~~~~~~~~~~~~~~~~~~~~~

- New methods :py:meth:`Dataset.to_stacked_array` and
  :py:meth:`DataArray.to_unstacked_dataset` for reshaping Datasets of variables
  with different dimensions
  (:issue:`1317`).
  This is useful for feeding data from xarray into machine learning models,
  as described in :ref:`reshape.stacking_different`.
  By `Noah Brenowitz <https://github.com/nbren12>`_.

Enhancements
~~~~~~~~~~~~

- Support for renaming ``Dataset`` variables and dimensions independently
  with :py:meth:`~Dataset.rename_vars` and :py:meth:`~Dataset.rename_dims`
  (:issue:`3026`).
  By `Julia Kent <https://github.com/jukent>`_.

- Add ``scales``, ``offsets``, ``units`` and ``descriptions``
  attributes to :py:class:`~xarray.DataArray` returned by
  :py:func:`~xarray.open_rasterio`. (:issue:`3013`)
  By `Erle Carrara <https://github.com/ecarrara>`_.

Bug fixes
~~~~~~~~~

- Resolved deprecation warnings from newer versions of matplotlib and dask.
- Compatibility fixes for the upcoming pandas 0.25 and NumPy 1.17 releases.
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- Fix summaries for multiindex coordinates (:issue:`3079`).
  By `Jonas Hörsch <https://github.com/coroa>`_.
- Fix HDF5 error that could arise when reading multiple groups from a file at
  once (:issue:`2954`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

.. _whats-new.0.12.2:

v0.12.2 (29 June 2019)
----------------------

New functions/methods
~~~~~~~~~~~~~~~~~~~~~

- Two new functions, :py:func:`~xarray.combine_nested` and
  :py:func:`~xarray.combine_by_coords`, allow for combining datasets along any
  number of dimensions, instead of the one-dimensional list of datasets
  supported by :py:func:`~xarray.concat`.

  The new ``combine_nested`` will accept the datasets as a nested
  list-of-lists, and combine by applying a series of concat and merge
  operations. The new ``combine_by_coords`` instead uses the dimension
  coordinates of datasets to order them.

  :py:func:`~xarray.open_mfdataset` can use either ``combine_nested`` or
  ``combine_by_coords`` to combine datasets along multiple dimensions, by
  specifying the argument ``combine='nested'`` or ``combine='by_coords'``.

  The older function ``auto_combine`` has been deprecated,
  because its functionality has been subsumed by the new functions.
  To avoid FutureWarnings switch to using ``combine_nested`` or
  ``combine_by_coords``, (or set the ``combine`` argument in
  ``open_mfdataset``). (:issue:`2159`)
  By `Tom Nicholas <https://github.com/TomNicholas>`_.

- :py:meth:`~xarray.DataArray.rolling_exp` and
  :py:meth:`~xarray.Dataset.rolling_exp` added, similar to pandas'
  ``pd.DataFrame.ewm`` method. Calling ``.mean`` on the resulting object
  will return an exponentially weighted moving average.
  By `Maximilian Roos <https://github.com/max-sixty>`_.

- New :py:func:`DataArray.str <core.accessor_str.StringAccessor>` for string
  related manipulations, based on ``pandas.Series.str``.
  By `0x0L <https://github.com/0x0L>`_.

- Added ``strftime`` method to ``.dt`` accessor, making it simpler to hand a
  datetime ``DataArray`` to other code expecting formatted dates and times.
  (:issue:`2090`). :py:meth:`~xarray.CFTimeIndex.strftime` is also now
  available on :py:class:`CFTimeIndex`.
  By `Alan Brammer <https://github.com/abrammer>`_ and
  `Ryan May <https://github.com/dopplershift>`_.

- ``GroupBy.quantile`` is now a method of ``GroupBy``
  objects  (:issue:`3018`).
  By `David Huard <https://github.com/huard>`_.

- Argument and return types are added to most methods on ``DataArray`` and
  ``Dataset``, allowing static type checking both within xarray and external
  libraries. Type checking with `mypy <http://mypy-lang.org/>`_ is enabled in
  CI (though not required yet).
  By `Guido Imperiale <https://github.com/crusaderky>`_
  and `Maximilian Roos <https://github.com/max-sixty>`_.

Enhancements to existing functionality
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Add ``keepdims`` argument for reduce operations (:issue:`2170`)
  By `Scott Wales <https://github.com/ScottWales>`_.
- Enable ``@`` operator for DataArray. This is equivalent to :py:meth:`DataArray.dot`
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- Add ``fill_value`` argument for reindex, align, and merge operations
  to enable custom fill values. (:issue:`2876`)
  By `Zach Griffith <https://github.com/zdgriffith>`_.
- :py:meth:`DataArray.transpose` now accepts a keyword argument
  ``transpose_coords`` which enables transposition of coordinates in the
  same way as :py:meth:`Dataset.transpose`. :py:meth:`DataArray.groupby`
  :py:meth:`DataArray.groupby_bins`, and :py:meth:`DataArray.resample` now
  accept a keyword argument ``restore_coord_dims`` which keeps the order
  of the dimensions of multi-dimensional coordinates intact (:issue:`1856`).
  By `Peter Hausamann <https://github.com/phausamann>`_.
- Clean up Python 2 compatibility in code (:issue:`2950`)
  By `Guido Imperiale <https://github.com/crusaderky>`_.
- Better warning message when supplying invalid objects to ``xr.merge``
  (:issue:`2948`).  By `Mathias Hauser <https://github.com/mathause>`_.
- Add ``errors`` keyword argument to ``Dataset.drop`` and :py:meth:`Dataset.drop_dims`
  that allows ignoring errors if a passed label or dimension is not in the dataset
  (:issue:`2994`).
  By `Andrew Ross <https://github.com/andrew-c-ross>`_.

IO related enhancements
~~~~~~~~~~~~~~~~~~~~~~~

- Implement :py:func:`~xarray.load_dataset` and
  :py:func:`~xarray.load_dataarray` as alternatives to
  :py:func:`~xarray.open_dataset` and :py:func:`~xarray.open_dataarray` to
  open, load into memory, and close files, returning the Dataset or DataArray.
  These functions are helpful for avoiding file-lock errors when trying to
  write to files opened using ``open_dataset()`` or ``open_dataarray()``.
  (:issue:`2887`)
  By `Dan Nowacki <https://github.com/dnowacki-usgs>`_.
- It is now possible to extend existing :ref:`io.zarr` datasets, by using
  ``mode='a'`` and the new ``append_dim`` argument in
  :py:meth:`~xarray.Dataset.to_zarr`.
  By `Jendrik Jördening <https://github.com/jendrikjoe>`_,
  `David Brochart <https://github.com/davidbrochart>`_,
  `Ryan Abernathey <https://github.com/rabernat>`_ and
  `Shikhar Goenka <https://github.com/shikharsg>`_.
- ``xr.open_zarr`` now accepts manually specified chunks with the ``chunks=``
  parameter. ``auto_chunk=True`` is equivalent to ``chunks='auto'`` for
  backwards compatibility. The ``overwrite_encoded_chunks`` parameter is
  added to remove the original zarr chunk encoding.
  By `Lily Wang <https://github.com/lilyminium>`_.
- netCDF chunksizes are now only dropped when original_shape is different,
  not when it isn't found. (:issue:`2207`)
  By `Karel van de Plassche <https://github.com/Karel-van-de-Plassche>`_.
- Character arrays' character dimension name decoding and encoding handled by
  ``var.encoding['char_dim_name']`` (:issue:`2895`)
  By `James McCreight <https://github.com/jmccreight>`_.
- open_rasterio() now supports rasterio.vrt.WarpedVRT with custom transform,
  width and height (:issue:`2864`).
  By `Julien Michel <https://github.com/jmichel-otb>`_.

Bug fixes
~~~~~~~~~

- Rolling operations on xarray objects containing dask arrays could silently
  compute the incorrect result or use large amounts of memory (:issue:`2940`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- Don't set encoding attributes on bounds variables when writing to netCDF.
  (:issue:`2921`)
  By `Deepak Cherian <https://github.com/dcherian>`_.
- NetCDF4 output: variables with unlimited dimensions must be chunked (not
  contiguous) on output. (:issue:`1849`)
  By `James McCreight <https://github.com/jmccreight>`_.
- indexing with an empty list creates an object with zero-length axis (:issue:`2882`)
  By `Mayeul d'Avezac <https://github.com/mdavezac>`_.
- Return correct count for scalar datetime64 arrays (:issue:`2770`)
  By `Dan Nowacki <https://github.com/dnowacki-usgs>`_.
- Fixed max, min exception when applied to a multiIndex (:issue:`2923`)
  By `Ian Castleden <https://github.com/arabidopsis>`_
- A deep copy deep-copies the coords (:issue:`1463`)
  By `Martin Pletcher <https://github.com/pletchm>`_.
- Increased support for `missing_value` (:issue:`2871`)
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Removed usages of `pytest.config`, which is deprecated (:issue:`2988`)
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- Fixed performance issues with cftime installed (:issue:`3000`)
  By `0x0L <https://github.com/0x0L>`_.
- Replace incorrect usages of `message` in pytest assertions
  with `match` (:issue:`3011`)
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- Add explicit pytest markers, now required by pytest
  (:issue:`3032`).
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- Test suite fixes for newer versions of pytest (:issue:`3011`, :issue:`3032`).
  By `Maximilian Roos <https://github.com/max-sixty>`_
  and `Stephan Hoyer <https://github.com/shoyer>`_.

.. _whats-new.0.12.1:

v0.12.1 (4 April 2019)
----------------------

Enhancements
~~~~~~~~~~~~

- Allow ``expand_dims`` method to support inserting/broadcasting dimensions
  with size > 1. (:issue:`2710`)
  By `Martin Pletcher <https://github.com/pletchm>`_.

Bug fixes
~~~~~~~~~

- Dataset.copy(deep=True) now creates a deep copy of the attrs (:issue:`2835`).
  By `Andras Gefferth <https://github.com/kefirbandi>`_.
- Fix incorrect ``indexes`` resulting from various ``Dataset`` operations
  (e.g., ``swap_dims``, ``isel``, ``reindex``, ``[]``) (:issue:`2842`,
  :issue:`2856`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

.. _whats-new.0.12.0:

v0.12.0 (15 March 2019)
-----------------------

Highlights include:

- Removed support for Python 2. This is the first version of xarray that is
  Python 3 only!
- New :py:meth:`~xarray.DataArray.coarsen` and
  :py:meth:`~xarray.DataArray.integrate` methods. See :ref:`comput.coarsen`
  and :ref:`compute.using_coordinates` for details.
- Many improvements to cftime support. See below for details.

Deprecations
~~~~~~~~~~~~

- The ``compat`` argument to ``Dataset`` and the ``encoding`` argument to
  ``DataArray`` are deprecated and will be removed in a future release.
  (:issue:`1188`)
  By `Maximilian Roos <https://github.com/max-sixty>`_.

cftime related enhancements
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Resampling of standard and non-standard calendars indexed by
  :py:class:`~xarray.CFTimeIndex` is now possible. (:issue:`2191`).
  By `Jwen Fai Low <https://github.com/jwenfai>`_ and
  `Spencer Clark <https://github.com/spencerkclark>`_.

- Taking the mean of arrays of :py:class:`cftime.datetime` objects, and
  by extension, use of :py:meth:`~xarray.DataArray.coarsen` with
  :py:class:`cftime.datetime` coordinates is now possible. By `Spencer Clark
  <https://github.com/spencerkclark>`_.

- Internal plotting now supports ``cftime.datetime`` objects as time series.
  (:issue:`2164`)
  By `Julius Busecke <https://github.com/jbusecke>`_ and
  `Spencer Clark <https://github.com/spencerkclark>`_.

- :py:meth:`~xarray.cftime_range` now supports QuarterBegin and QuarterEnd offsets (:issue:`2663`).
  By `Jwen Fai Low <https://github.com/jwenfai>`_

- :py:meth:`~xarray.open_dataset` now accepts a ``use_cftime`` argument, which
  can be used to require that ``cftime.datetime`` objects are always used, or
  never used when decoding dates encoded with a standard calendar.  This can be
  used to ensure consistent date types are returned when using
  :py:meth:`~xarray.open_mfdataset` (:issue:`1263`) and/or to silence
  serialization warnings raised if dates from a standard calendar are found to
  be outside the :py:class:`pandas.Timestamp`-valid range (:issue:`2754`).  By
  `Spencer Clark <https://github.com/spencerkclark>`_.

- :py:meth:`pandas.Series.dropna` is now supported for a
  :py:class:`pandas.Series` indexed by a :py:class:`~xarray.CFTimeIndex`
  (:issue:`2688`). By `Spencer Clark <https://github.com/spencerkclark>`_.

Other enhancements
~~~~~~~~~~~~~~~~~~

- Added ability to open netcdf4/hdf5 file-like objects with ``open_dataset``.
  Requires (h5netcdf>0.7 and h5py>2.9.0). (:issue:`2781`)
  By `Scott Henderson <https://github.com/scottyhq>`_
- Add ``data=False`` option to ``to_dict()`` methods. (:issue:`2656`)
  By `Ryan Abernathey <https://github.com/rabernat>`_
- :py:meth:`DataArray.coarsen` and
  :py:meth:`Dataset.coarsen` are newly added.
  See :ref:`comput.coarsen` for details.
  (:issue:`2525`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- Upsampling an array via interpolation with resample is now dask-compatible,
  as long as the array is not chunked along the resampling dimension.
  By `Spencer Clark <https://github.com/spencerkclark>`_.
- :py:func:`xarray.testing.assert_equal` and
  :py:func:`xarray.testing.assert_identical` now provide a more detailed
  report showing what exactly differs between the two objects (dimensions /
  coordinates / variables / attributes)  (:issue:`1507`).
  By `Benoit Bovy <https://github.com/benbovy>`_.
- Add ``tolerance`` option to ``resample()`` methods ``bfill``, ``pad``,
  ``nearest``. (:issue:`2695`)
  By `Hauke Schulz <https://github.com/observingClouds>`_.
- :py:meth:`DataArray.integrate` and
  :py:meth:`Dataset.integrate` are newly added.
  See :ref:`compute.using_coordinates` for the detail.
  (:issue:`1332`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- Added :py:meth:`~xarray.Dataset.drop_dims` (:issue:`1949`).
  By `Kevin Squire <https://github.com/kmsquire>`_.

Bug fixes
~~~~~~~~~

- Silenced warnings that appear when using pandas 0.24.
  By `Stephan Hoyer <https://github.com/shoyer>`_
- Interpolating via resample now internally specifies ``bounds_error=False``
  as an argument to ``scipy.interpolate.interp1d``, allowing for interpolation
  from higher frequencies to lower frequencies.  Datapoints outside the bounds
  of the original time coordinate are now filled with NaN (:issue:`2197`). By
  `Spencer Clark <https://github.com/spencerkclark>`_.
- Line plots with the ``x`` argument set to a non-dimensional coord now plot
  the correct data for 1D DataArrays.
  (:issue:`2725`). By `Tom Nicholas <https://github.com/TomNicholas>`_.
- Subtracting a scalar ``cftime.datetime`` object from a
  :py:class:`CFTimeIndex` now results in a :py:class:`pandas.TimedeltaIndex`
  instead of raising a ``TypeError`` (:issue:`2671`).  By `Spencer Clark
  <https://github.com/spencerkclark>`_.
- backend_kwargs are no longer ignored when using open_dataset with pynio engine
  (:issue:'2380')
  By `Jonathan Joyce <https://github.com/jonmjoyce>`_.
- Fix ``open_rasterio`` creating a WKT CRS instead of PROJ.4 with
  ``rasterio`` 1.0.14+ (:issue:`2715`).
  By `David Hoese <https://github.com/djhoese>`_.
- Masking data arrays with :py:meth:`xarray.DataArray.where` now returns an
  array with the name of the original masked array (:issue:`2748` and :issue:`2457`).
  By `Yohai Bar-Sinai <https://github.com/yohai>`_.
- Fixed error when trying to reduce a DataArray using a function which does not
  require an axis argument. (:issue:`2768`)
  By `Tom Nicholas <https://github.com/TomNicholas>`_.
- Concatenating a sequence of :py:class:`~xarray.DataArray` with varying names
  sets the name of the output array to ``None``, instead of the name of the
  first input array. If the names are the same it sets the name to that,
  instead to the name of the first DataArray in the list as it did before.
  (:issue:`2775`). By `Tom Nicholas <https://github.com/TomNicholas>`_.

- Per the `CF conventions section on calendars
  <http://cfconventions.org/cf-conventions/cf-conventions.html#calendar>`_,
  specifying ``'standard'`` as the calendar type in
  :py:meth:`~xarray.cftime_range` now correctly refers to the ``'gregorian'``
  calendar instead of the ``'proleptic_gregorian'`` calendar (:issue:`2761`).

.. _whats-new.0.11.3:

v0.11.3 (26 January 2019)
-------------------------

Bug fixes
~~~~~~~~~

- Saving files with times encoded with reference dates with timezones
  (e.g. '2000-01-01T00:00:00-05:00') no longer raises an error
  (:issue:`2649`).  By `Spencer Clark <https://github.com/spencerkclark>`_.
- Fixed performance regression with ``open_mfdataset`` (:issue:`2662`).
  By `Tom Nicholas <https://github.com/TomNicholas>`_.
- Fixed supplying an explicit dimension in the ``concat_dim`` argument to
  to ``open_mfdataset`` (:issue:`2647`).
  By `Ben Root <https://github.com/WeatherGod>`_.

.. _whats-new.0.11.2:

v0.11.2 (2 January 2019)
------------------------

Removes inadvertently introduced setup dependency on pytest-runner
(:issue:`2641`). Otherwise, this release is exactly equivalent to 0.11.1.

.. warning::

  This is the last xarray release that will support Python 2.7. Future releases
  will be Python 3 only, but older versions of xarray will always be available
  for Python 2.7 users. For the more details, see:

  - `Xarray Github issue discussing dropping Python 2 <https://github.com/pydata/xarray/issues/1829>`__
  - `Python 3 Statement <http://www.python3statement.org/>`__
  - `Tips on porting to Python 3 <https://docs.python.org/3/howto/pyporting.html>`__

.. _whats-new.0.11.1:

v0.11.1 (29 December 2018)
--------------------------

This minor release includes a number of enhancements and bug fixes, and two
(slightly) breaking changes.

Breaking changes
~~~~~~~~~~~~~~~~

- Minimum rasterio version increased from 0.36 to 1.0 (for ``open_rasterio``)
- Time bounds variables are now also decoded according to CF conventions
  (:issue:`2565`). The previous behavior was to decode them only if they
  had specific time attributes, now these attributes are copied
  automatically from the corresponding time coordinate. This might
  break downstream code that was relying on these variables to be
  brake downstream code that was relying on these variables to be
  not decoded.
  By `Fabien Maussion <https://github.com/fmaussion>`_.

Enhancements
~~~~~~~~~~~~

- Ability to read and write consolidated metadata in zarr stores (:issue:`2558`).
  By `Ryan Abernathey <https://github.com/rabernat>`_.
- :py:class:`CFTimeIndex` uses slicing for string indexing when possible (like
  :py:class:`pandas.DatetimeIndex`), which avoids unnecessary copies.
  By `Stephan Hoyer <https://github.com/shoyer>`_
- Enable passing ``rasterio.io.DatasetReader`` or ``rasterio.vrt.WarpedVRT`` to
  ``open_rasterio`` instead of file path string. Allows for in-memory
  reprojection, see  (:issue:`2588`).
  By `Scott Henderson <https://github.com/scottyhq>`_.
- Like :py:class:`pandas.DatetimeIndex`, :py:class:`CFTimeIndex` now supports
  "dayofyear" and "dayofweek" accessors (:issue:`2597`).  Note this requires a
  version of cftime greater than 1.0.2.  By `Spencer Clark
  <https://github.com/spencerkclark>`_.
- The option ``'warn_for_unclosed_files'`` (False by default) has been added to
  allow users to enable a warning when files opened by xarray are deallocated
  but were not explicitly closed. This is mostly useful for debugging; we
  recommend enabling it in your test suites if you use xarray for IO.
  By `Stephan Hoyer <https://github.com/shoyer>`_
- Support Dask ``HighLevelGraphs`` by `Matthew Rocklin <https://github.com/mrocklin>`_.
- :py:meth:`DataArray.resample` and :py:meth:`Dataset.resample` now supports the
  ``loffset`` kwarg just like pandas.
  By `Deepak Cherian <https://github.com/dcherian>`_
- Datasets are now guaranteed to have a ``'source'`` encoding, so the source
  file name is always stored (:issue:`2550`).
  By `Tom Nicholas <https://github.com/TomNicholas>`_.
- The ``apply`` methods for ``DatasetGroupBy``, ``DataArrayGroupBy``,
  ``DatasetResample`` and ``DataArrayResample`` now support passing positional
  arguments to the applied function as a tuple to the ``args`` argument.
  By `Matti Eskelinen <https://github.com/maaleske>`_.
- 0d slices of ndarrays are now obtained directly through indexing, rather than
  extracting and wrapping a scalar, avoiding unnecessary copying. By `Daniel
  Wennberg <https://github.com/danielwe>`_.
- Added support for ``fill_value`` with
  :py:meth:`~xarray.DataArray.shift` and :py:meth:`~xarray.Dataset.shift`
  By `Maximilian Roos <https://github.com/max-sixty>`_

Bug fixes
~~~~~~~~~

- Ensure files are automatically closed, if possible, when no longer referenced
  by a Python variable (:issue:`2560`).
  By `Stephan Hoyer <https://github.com/shoyer>`_
- Fixed possible race conditions when reading/writing to disk in parallel
  (:issue:`2595`).
  By `Stephan Hoyer <https://github.com/shoyer>`_
- Fix h5netcdf saving scalars with filters or chunks (:issue:`2563`).
  By `Martin Raspaud <https://github.com/mraspaud>`_.
- Fix parsing of ``_Unsigned`` attribute set by OPENDAP servers. (:issue:`2583`).
  By `Deepak Cherian <https://github.com/dcherian>`_
- Fix failure in time encoding when exporting to netCDF with versions of pandas
  less than 0.21.1 (:issue:`2623`).  By `Spencer Clark
  <https://github.com/spencerkclark>`_.
- Fix MultiIndex selection to update label and level (:issue:`2619`).
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

.. _whats-new.0.11.0:

v0.11.0 (7 November 2018)
-------------------------

Breaking changes
~~~~~~~~~~~~~~~~

- Finished deprecations (changed behavior with this release):

  - ``Dataset.T`` has been removed as a shortcut for :py:meth:`Dataset.transpose`.
    Call :py:meth:`Dataset.transpose` directly instead.
  - Iterating over a ``Dataset`` now includes only data variables, not coordinates.
    Similarily, calling ``len`` and ``bool`` on a ``Dataset`` now
    includes only data variables.
  - ``DataArray.__contains__`` (used by Python's ``in`` operator) now checks
    array data, not coordinates.
  - The old resample syntax from before xarray 0.10, e.g.,
    ``data.resample('1D', dim='time', how='mean')``, is no longer supported will
    raise an error in most cases. You need to use the new resample syntax
    instead, e.g., ``data.resample(time='1D').mean()`` or
    ``data.resample({'time': '1D'}).mean()``.


- New deprecations (behavior will be changed in xarray 0.12):

  - Reduction of :py:meth:`DataArray.groupby` and :py:meth:`DataArray.resample`
    without dimension argument will change in the next release.
    Now we warn a FutureWarning.
    By `Keisuke Fujii <https://github.com/fujiisoup>`_.
  - The ``inplace`` kwarg of a number of `DataArray` and `Dataset` methods is being
    deprecated and will be removed in the next release.
    By `Deepak Cherian <https://github.com/dcherian>`_.


- Refactored storage backends:

  - Xarray's storage backends now automatically open and close files when
    necessary, rather than requiring opening a file with ``autoclose=True``. A
    global least-recently-used cache is used to store open files; the default
    limit of 128 open files should suffice in most cases, but can be adjusted if
    necessary with
    ``xarray.set_options(file_cache_maxsize=...)``. The ``autoclose`` argument
    to ``open_dataset`` and related functions has been deprecated and is now a
    no-op.

    This change, along with an internal refactor of xarray's storage backends,
    should significantly improve performance when reading and writing
    netCDF files with Dask, especially when working with many files or using
    Dask Distributed. By `Stephan Hoyer <https://github.com/shoyer>`_


- Support for non-standard calendars used in climate science:

  - Xarray will now always use :py:class:`cftime.datetime` objects, rather
    than by default trying to coerce them into ``np.datetime64[ns]`` objects.
    A :py:class:`~xarray.CFTimeIndex` will be used for indexing along time
    coordinates in these cases.
  - A new method :py:meth:`~xarray.CFTimeIndex.to_datetimeindex` has been added
    to aid in converting from a  :py:class:`~xarray.CFTimeIndex` to a
    :py:class:`pandas.DatetimeIndex` for the remaining use-cases where
    using a :py:class:`~xarray.CFTimeIndex` is still a limitation (e.g. for
    resample or plotting).
  - Setting the ``enable_cftimeindex`` option is now a no-op and emits a
    ``FutureWarning``.

Enhancements
~~~~~~~~~~~~

- :py:meth:`xarray.DataArray.plot.line` can now accept multidimensional
  coordinate variables as input. `hue` must be a dimension name in this case.
  (:issue:`2407`)
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Added support for Python 3.7. (:issue:`2271`).
  By `Joe Hamman <https://github.com/jhamman>`_.
- Added support for plotting data with `pandas.Interval` coordinates, such as those
  created by :py:meth:`~xarray.DataArray.groupby_bins`
  By `Maximilian Maahn <https://github.com/maahn>`_.
- Added :py:meth:`~xarray.CFTimeIndex.shift` for shifting the values of a
  CFTimeIndex by a specified frequency. (:issue:`2244`).
  By `Spencer Clark <https://github.com/spencerkclark>`_.
- Added support for using ``cftime.datetime`` coordinates with
  :py:meth:`~xarray.DataArray.differentiate`,
  :py:meth:`~xarray.Dataset.differentiate`,
  :py:meth:`~xarray.DataArray.interp`, and
  :py:meth:`~xarray.Dataset.interp`.
  By `Spencer Clark <https://github.com/spencerkclark>`_
- There is now a global option to either always keep or always discard
  dataset and dataarray attrs upon operations. The option is set with
  ``xarray.set_options(keep_attrs=True)``, and the default is to use the old
  behaviour.
  By `Tom Nicholas <https://github.com/TomNicholas>`_.
- Added a new backend for the GRIB file format based on ECMWF *cfgrib*
  python driver and *ecCodes* C-library. (:issue:`2475`)
  By `Alessandro Amici <https://github.com/alexamici>`_,
  sponsored by `ECMWF <https://github.com/ecmwf>`_.
- Resample now supports a dictionary mapping from dimension to frequency as
  its first argument, e.g., ``data.resample({'time': '1D'}).mean()``. This is
  consistent with other xarray functions that accept either dictionaries or
  keyword arguments. By `Stephan Hoyer <https://github.com/shoyer>`_.

- The preferred way to access tutorial data is now to load it lazily with
  :py:meth:`xarray.tutorial.open_dataset`.
  :py:meth:`xarray.tutorial.load_dataset` calls `Dataset.load()` prior
  to returning (and is now deprecated). This was changed in order to facilitate
  using tutorial datasets with dask.
  By `Joe Hamman <https://github.com/jhamman>`_.
- ``DataArray`` can now use ``xr.set_option(keep_attrs=True)`` and retain attributes in binary operations,
  such as (``+, -, * ,/``). Default behaviour is unchanged (*Attributes will be dismissed*). By `Michael Blaschek <https://github.com/MBlaschek>`_

Bug fixes
~~~~~~~~~

- ``FacetGrid`` now properly uses the ``cbar_kwargs`` keyword argument.
  (:issue:`1504`, :issue:`1717`)
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Addition and subtraction operators used with a CFTimeIndex now preserve the
  index's type. (:issue:`2244`).
  By `Spencer Clark <https://github.com/spencerkclark>`_.
- We now properly handle arrays of ``datetime.datetime`` and ``datetime.timedelta``
  provided as coordinates. (:issue:`2512`)
  By `Deepak Cherian <https://github.com/dcherian>`_.
- ``xarray.DataArray.roll`` correctly handles multidimensional arrays.
  (:issue:`2445`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- ``xarray.plot()`` now properly accepts a ``norm`` argument and does not override
  the norm's ``vmin`` and ``vmax``. (:issue:`2381`)
  By `Deepak Cherian <https://github.com/dcherian>`_.
- ``xarray.DataArray.std()`` now correctly accepts ``ddof`` keyword argument.
  (:issue:`2240`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- Restore matplotlib's default of plotting dashed negative contours when
  a single color is passed to ``DataArray.contour()`` e.g. ``colors='k'``.
  By `Deepak Cherian <https://github.com/dcherian>`_.


- Fix a bug that caused some indexing operations on arrays opened with
  ``open_rasterio`` to error (:issue:`2454`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- Subtracting one CFTimeIndex from another now returns a
  ``pandas.TimedeltaIndex``, analogous to the behavior for DatetimeIndexes
  (:issue:`2484`).  By `Spencer Clark <https://github.com/spencerkclark>`_.
- Adding a TimedeltaIndex to, or subtracting a TimedeltaIndex from a
  CFTimeIndex is now allowed (:issue:`2484`).
  By `Spencer Clark <https://github.com/spencerkclark>`_.
- Avoid use of Dask's deprecated ``get=`` parameter in tests
  by `Matthew Rocklin <https://github.com/mrocklin>`_.
- An ``OverflowError`` is now accurately raised and caught during the
  encoding process if a reference date is used that is so distant that
  the dates must be encoded using cftime rather than NumPy (:issue:`2272`).
  By `Spencer Clark <https://github.com/spencerkclark>`_.

- Chunked datasets can now roundtrip to Zarr storage continually
  with `to_zarr` and ``open_zarr`` (:issue:`2300`).
  By `Lily Wang <https://github.com/lilyminium>`_.

.. _whats-new.0.10.9:

v0.10.9 (21 September 2018)
---------------------------

This minor release contains a number of backwards compatible enhancements.

Announcements of note:

- Xarray is now a NumFOCUS fiscally sponsored project! Read
  `the announcement <https://numfocus.org/blog/xarray-joins-numfocus-sponsored-projects>`_
  for more details.
- We have a new :doc:`roadmap` that outlines our future development plans.

- ``Dataset.apply`` now properly documents the way `func` is called.
  By `Matti Eskelinen <https://github.com/maaleske>`_.

Enhancements
~~~~~~~~~~~~

- :py:meth:`~xarray.DataArray.differentiate` and
  :py:meth:`~xarray.Dataset.differentiate` are newly added.
  (:issue:`1332`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- Default colormap for sequential and divergent data can now be set via
  :py:func:`~xarray.set_options()`
  (:issue:`2394`)
  By `Julius Busecke <https://github.com/jbusecke>`_.

- min_count option is newly supported in :py:meth:`~xarray.DataArray.sum`,
  :py:meth:`~xarray.DataArray.prod` and :py:meth:`~xarray.Dataset.sum`, and
  :py:meth:`~xarray.Dataset.prod`.
  (:issue:`2230`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- :py:func:`~plot.plot()` now accepts the kwargs
  ``xscale, yscale, xlim, ylim, xticks, yticks`` just like pandas. Also ``xincrease=False, yincrease=False`` now use matplotlib's axis inverting methods instead of setting limits.
  By `Deepak Cherian <https://github.com/dcherian>`_. (:issue:`2224`)

- DataArray coordinates and Dataset coordinates and data variables are
  now displayed as `a b ... y z` rather than `a b c d ...`.
  (:issue:`1186`)
  By `Seth P <https://github.com/seth-p>`_.
- A new CFTimeIndex-enabled :py:func:`cftime_range` function for use in
  generating dates from standard or non-standard calendars.  By `Spencer Clark
  <https://github.com/spencerkclark>`_.

- When interpolating over a ``datetime64`` axis, you can now provide a datetime string instead of a ``datetime64`` object. E.g. ``da.interp(time='1991-02-01')``
  (:issue:`2284`)
  By `Deepak Cherian <https://github.com/dcherian>`_.

- A clear error message is now displayed if a ``set`` or ``dict`` is passed in place of an array
  (:issue:`2331`)
  By `Maximilian Roos <https://github.com/max-sixty>`_.

- Applying ``unstack`` to a large DataArray or Dataset is now much faster if the MultiIndex has not been modified after stacking the indices.
  (:issue:`1560`)
  By `Maximilian Maahn <https://github.com/maahn>`_.

- You can now control whether or not to offset the coordinates when using
  the ``roll`` method and the current behavior, coordinates rolled by default,
  raises a deprecation warning unless explicitly setting the keyword argument.
  (:issue:`1875`)
  By `Andrew Huang <https://github.com/ahuang11>`_.

- You can now call ``unstack`` without arguments to unstack every MultiIndex in a DataArray or Dataset.
  By `Julia Signell <https://github.com/jsignell>`_.

- Added the ability to pass a data kwarg to ``copy`` to create a new object with the
  same metadata as the original object but using new values.
  By `Julia Signell <https://github.com/jsignell>`_.

Bug fixes
~~~~~~~~~

- ``xarray.plot.imshow()`` correctly uses the ``origin`` argument.
  (:issue:`2379`)
  By `Deepak Cherian <https://github.com/dcherian>`_.

- Fixed ``DataArray.to_iris()`` failure while creating ``DimCoord`` by
  falling back to creating ``AuxCoord``. Fixed dependency on ``var_name``
  attribute being set.
  (:issue:`2201`)
  By `Thomas Voigt <https://github.com/tv3141>`_.
- Fixed a bug in ``zarr`` backend which prevented use with datasets with
  invalid chunk size encoding after reading from an existing store
  (:issue:`2278`).
  By `Joe Hamman <https://github.com/jhamman>`_.

- Tests can be run in parallel with pytest-xdist
  By `Tony Tung <https://github.com/ttung>`_.

- Follow up the renamings in dask; from dask.ghost to dask.overlap
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- Now raises a ValueError when there is a conflict between dimension names and
  level names of MultiIndex. (:issue:`2299`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- Follow up the renamings in dask; from dask.ghost to dask.overlap
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- Now :py:func:`~xarray.apply_ufunc` raises a ValueError when the size of
  ``input_core_dims`` is inconsistent with the number of arguments.
  (:issue:`2341`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- Fixed ``Dataset.filter_by_attrs()`` behavior not matching ``netCDF4.Dataset.get_variables_by_attributes()``.
  When more than one ``key=value`` is passed into ``Dataset.filter_by_attrs()`` it will now return a Dataset with variables which pass
  all the filters.
  (:issue:`2315`)
  By `Andrew Barna <https://github.com/docotak>`_.

.. _whats-new.0.10.8:

v0.10.8 (18 July 2018)
----------------------

Breaking changes
~~~~~~~~~~~~~~~~

- Xarray no longer supports python 3.4. Additionally, the minimum supported
  versions of the following dependencies has been updated and/or clarified:

  - pandas: 0.18 -> 0.19
  - NumPy: 1.11 -> 1.12
  - Dask: 0.9 -> 0.16
  - Matplotlib: unspecified -> 1.5

  (:issue:`2204`). By `Joe Hamman <https://github.com/jhamman>`_.

Enhancements
~~~~~~~~~~~~

- :py:meth:`~xarray.DataArray.interp_like` and
  :py:meth:`~xarray.Dataset.interp_like` methods are newly added.
  (:issue:`2218`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- Added support for curvilinear and unstructured generic grids
  to :py:meth:`~xarray.DataArray.to_cdms2` and
  :py:meth:`~xarray.DataArray.from_cdms2` (:issue:`2262`).
  By `Stephane Raynaud <https://github.com/stefraynaud>`_.

Bug fixes
~~~~~~~~~

- Fixed a bug in ``zarr`` backend which prevented use with datasets with
  incomplete chunks in multiple dimensions (:issue:`2225`).
  By `Joe Hamman <https://github.com/jhamman>`_.

- Fixed a bug in :py:meth:`~Dataset.to_netcdf` which prevented writing
  datasets when the arrays had different chunk sizes (:issue:`2254`).
  By `Mike Neish <https://github.com/neishm>`_.

- Fixed masking during the conversion to cdms2 objects by
  :py:meth:`~xarray.DataArray.to_cdms2` (:issue:`2262`).
  By `Stephane Raynaud <https://github.com/stefraynaud>`_.

- Fixed a bug in 2D plots which incorrectly raised an error when 2D coordinates
  weren't monotonic (:issue:`2250`).
  By `Fabien Maussion <https://github.com/fmaussion>`_.

- Fixed warning raised in :py:meth:`~Dataset.to_netcdf` due to deprecation of
  `effective_get` in dask (:issue:`2238`).
  By `Joe Hamman <https://github.com/jhamman>`_.

.. _whats-new.0.10.7:

v0.10.7 (7 June 2018)
---------------------

Enhancements
~~~~~~~~~~~~

- Plot labels now make use of metadata that follow CF conventions
  (:issue:`2135`).
  By `Deepak Cherian <https://github.com/dcherian>`_ and `Ryan Abernathey <https://github.com/rabernat>`_.

- Line plots now support facetting with ``row`` and ``col`` arguments
  (:issue:`2107`).
  By `Yohai Bar Sinai <https://github.com/yohai>`_.

- :py:meth:`~xarray.DataArray.interp` and :py:meth:`~xarray.Dataset.interp`
  methods are newly added.
  See :ref:`interp` for the detail.
  (:issue:`2079`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

Bug fixes
~~~~~~~~~

- Fixed a bug in ``rasterio`` backend which prevented use with ``distributed``.
  The ``rasterio`` backend now returns pickleable objects (:issue:`2021`).
  By `Joe Hamman <https://github.com/jhamman>`_.

.. _whats-new.0.10.6:

v0.10.6 (31 May 2018)
---------------------

The minor release includes a number of bug-fixes and backwards compatible
enhancements.

Enhancements
~~~~~~~~~~~~

- New PseudoNetCDF backend for many Atmospheric data formats including
  GEOS-Chem, CAMx, NOAA arlpacked bit and many others. See
  :ref:`io.PseudoNetCDF` for more details.
  By `Barron Henderson <https://github.com/barronh>`_.

- The :py:class:`Dataset` constructor now aligns :py:class:`DataArray`
  arguments in ``data_vars`` to indexes set explicitly in ``coords``,
  where previously an error would be raised.
  (:issue:`674`)
  By `Maximilian Roos <https://github.com/max-sixty>`_.

- :py:meth:`~DataArray.sel`, :py:meth:`~DataArray.isel` & :py:meth:`~DataArray.reindex`,
  (and their :py:class:`Dataset` counterparts) now support supplying a ``dict``
  as a first argument, as an alternative to the existing approach
  of supplying `kwargs`. This allows for more robust behavior
  of dimension names which conflict with other keyword names, or are
  not strings.
  By `Maximilian Roos <https://github.com/max-sixty>`_.

- :py:meth:`~DataArray.rename` now supports supplying ``**kwargs``, as an
  alternative to the existing approach of supplying a ``dict`` as the
  first argument.
  By `Maximilian Roos <https://github.com/max-sixty>`_.

- :py:meth:`~DataArray.cumsum` and :py:meth:`~DataArray.cumprod` now support
  aggregation over multiple dimensions at the same time. This is the default
  behavior when dimensions are not specified (previously this raised an error).
  By `Stephan Hoyer <https://github.com/shoyer>`_

- :py:meth:`DataArray.dot` and :py:func:`dot` are partly supported with older
  dask<0.17.4. (related to :issue:`2203`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- Xarray now uses `Versioneer <https://github.com/warner/python-versioneer>`__
  to manage its version strings. (:issue:`1300`).
  By `Joe Hamman <https://github.com/jhamman>`_.

Bug fixes
~~~~~~~~~

- Fixed a regression in 0.10.4, where explicitly specifying ``dtype='S1'`` or
  ``dtype=str`` in ``encoding`` with ``to_netcdf()`` raised an error
  (:issue:`2149`).
  `Stephan Hoyer <https://github.com/shoyer>`_

- :py:func:`apply_ufunc` now directly validates output variables
  (:issue:`1931`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- Fixed a bug where ``to_netcdf(..., unlimited_dims='bar')`` yielded NetCDF
  files with spurious 0-length dimensions (i.e. ``b``, ``a``, and ``r``)
  (:issue:`2134`).
  By `Joe Hamman <https://github.com/jhamman>`_.

- Removed spurious warnings with ``Dataset.update(Dataset)`` (:issue:`2161`)
  and ``array.equals(array)`` when ``array`` contains ``NaT`` (:issue:`2162`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- Aggregations with :py:meth:`Dataset.reduce` (including ``mean``, ``sum``,
  etc) no longer drop unrelated coordinates (:issue:`1470`). Also fixed a
  bug where non-scalar data-variables that did not include the aggregation
  dimension were improperly skipped.
  By `Stephan Hoyer <https://github.com/shoyer>`_

- Fix :meth:`~DataArray.stack` with non-unique coordinates on pandas 0.23
  (:issue:`2160`).
  By `Stephan Hoyer <https://github.com/shoyer>`_

- Selecting data indexed by a length-1 ``CFTimeIndex`` with a slice of strings
  now behaves as it does when using a length-1 ``DatetimeIndex`` (i.e. it no
  longer falsely returns an empty array when the slice includes the value in
  the index) (:issue:`2165`).
  By `Spencer Clark <https://github.com/spencerkclark>`_.

- Fix ``DataArray.groupby().reduce()`` mutating coordinates on the input array
  when grouping over dimension coordinates with duplicated entries
  (:issue:`2153`).
  By `Stephan Hoyer <https://github.com/shoyer>`_

- Fix ``Dataset.to_netcdf()`` cannot create group with ``engine="h5netcdf"``
  (:issue:`2177`).
  By `Stephan Hoyer <https://github.com/shoyer>`_

.. _whats-new.0.10.4:

v0.10.4 (16 May 2018)
----------------------

The minor release includes a number of bug-fixes and backwards compatible
enhancements. A highlight is ``CFTimeIndex``, which offers support for
non-standard calendars used in climate modeling.

Documentation
~~~~~~~~~~~~~

- New FAQ entry, :ref:`ecosystem`.
  By `Deepak Cherian <https://github.com/dcherian>`_.
- :ref:`assigning_values` now includes examples on how to select and assign
  values to a :py:class:`~xarray.DataArray` with ``.loc``.
  By `Chiara Lepore <https://github.com/chiaral>`_.

Enhancements
~~~~~~~~~~~~

- Add an option for using a ``CFTimeIndex`` for indexing times with
  non-standard calendars and/or outside the Timestamp-valid range; this index
  enables a subset of the functionality of a standard
  ``pandas.DatetimeIndex``.
  See :ref:`CFTimeIndex` for full details.
  (:issue:`789`, :issue:`1084`, :issue:`1252`)
  By `Spencer Clark <https://github.com/spencerkclark>`_ with help from
  `Stephan Hoyer <https://github.com/shoyer>`_.
- Allow for serialization of ``cftime.datetime`` objects (:issue:`789`,
  :issue:`1084`, :issue:`2008`, :issue:`1252`) using the standalone ``cftime``
  library.
  By `Spencer Clark <https://github.com/spencerkclark>`_.
- Support writing lists of strings as netCDF attributes (:issue:`2044`).
  By `Dan Nowacki <https://github.com/dnowacki-usgs>`_.
- :py:meth:`~xarray.Dataset.to_netcdf` with ``engine='h5netcdf'`` now accepts h5py
  encoding settings ``compression`` and ``compression_opts``, along with the
  NetCDF4-Python style settings ``gzip=True`` and ``complevel``.
  This allows using any compression plugin installed in hdf5, e.g. LZF
  (:issue:`1536`). By `Guido Imperiale <https://github.com/crusaderky>`_.
- :py:meth:`~xarray.dot` on dask-backed data will now call :func:`dask.array.einsum`.
  This greatly boosts speed and allows chunking on the core dims.
  The function now requires dask >= 0.17.3 to work on dask-backed data
  (:issue:`2074`). By `Guido Imperiale <https://github.com/crusaderky>`_.
- ``plot.line()`` learned new kwargs: ``xincrease``, ``yincrease`` that change
  the direction of the respective axes.
  By `Deepak Cherian <https://github.com/dcherian>`_.

- Added the ``parallel`` option to :py:func:`open_mfdataset`. This option uses
  ``dask.delayed`` to parallelize the open and preprocessing steps within
  ``open_mfdataset``. This is expected to provide performance improvements when
  opening many files, particularly when used in conjunction with dask's
  multiprocessing or distributed schedulers (:issue:`1981`).
  By `Joe Hamman <https://github.com/jhamman>`_.

- New ``compute`` option in :py:meth:`~xarray.Dataset.to_netcdf`,
  :py:meth:`~xarray.Dataset.to_zarr`, and :py:func:`~xarray.save_mfdataset` to
  allow for the lazy computation of netCDF and zarr stores. This feature is
  currently only supported by the netCDF4 and zarr backends. (:issue:`1784`).
  By `Joe Hamman <https://github.com/jhamman>`_.


Bug fixes
~~~~~~~~~

- ``ValueError`` is raised when coordinates with the wrong size are assigned to
  a :py:class:`DataArray`. (:issue:`2112`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- Fixed a bug in :py:meth:`~xarray.DataArray.rolling` with bottleneck. Also,
  fixed a bug in rolling an integer dask array. (:issue:`2113`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- Fixed a bug where `keep_attrs=True` flag was neglected if
  :py:func:`apply_ufunc` was used with :py:class:`Variable`. (:issue:`2114`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- When assigning a :py:class:`DataArray` to :py:class:`Dataset`, any conflicted
  non-dimensional coordinates of the DataArray are now dropped.
  (:issue:`2068`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- Better error handling in ``open_mfdataset`` (:issue:`2077`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- ``plot.line()`` does not call ``autofmt_xdate()`` anymore. Instead it changes
  the rotation and horizontal alignment of labels without removing the x-axes of
  any other subplots in the figure (if any).
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Colorbar limits are now determined by excluding ±Infs too.
  By `Deepak Cherian <https://github.com/dcherian>`_.
  By `Joe Hamman <https://github.com/jhamman>`_.
- Fixed ``to_iris`` to maintain lazy dask array after conversion (:issue:`2046`).
  By `Alex Hilson <https://github.com/AlexHilson>`_ and `Stephan Hoyer <https://github.com/shoyer>`_.

.. _whats-new.0.10.3:

v0.10.3 (13 April 2018)
------------------------

The minor release includes a number of bug-fixes and backwards compatible enhancements.

Enhancements
~~~~~~~~~~~~

- :py:meth:`~xarray.DataArray.isin` and :py:meth:`~xarray.Dataset.isin` methods,
  which test each value in the array for whether it is contained in the
  supplied list, returning a bool array. See :ref:`selecting values with isin`
  for full details. Similar to the ``np.isin`` function.
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- Some speed improvement to construct :py:class:`~xarray.core.rolling.DataArrayRolling`
  object (:issue:`1993`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- Handle variables with different values for ``missing_value`` and
  ``_FillValue`` by masking values for both attributes; previously this
  resulted in a ``ValueError``. (:issue:`2016`)
  By `Ryan May <https://github.com/dopplershift>`_.

Bug fixes
~~~~~~~~~

- Fixed ``decode_cf`` function to operate lazily on dask arrays
  (:issue:`1372`). By `Ryan Abernathey <https://github.com/rabernat>`_.
- Fixed labeled indexing with slice bounds given by xarray objects with
  datetime64 or timedelta64 dtypes (:issue:`1240`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- Attempting to convert an xarray.Dataset into a numpy array now raises an
  informative error message.
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- Fixed a bug in decode_cf_datetime where ``int32`` arrays weren't parsed
  correctly (:issue:`2002`).
  By `Fabien Maussion <https://github.com/fmaussion>`_.
- When calling `xr.auto_combine()` or `xr.open_mfdataset()` with a `concat_dim`,
  the resulting dataset will have that one-element dimension (it was
  silently dropped, previously) (:issue:`1988`).
  By `Ben Root <https://github.com/WeatherGod>`_.

.. _whats-new.0.10.2:

v0.10.2 (13 March 2018)
-----------------------

The minor release includes a number of bug-fixes and enhancements, along with
one possibly **backwards incompatible change**.

Backwards incompatible changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- The addition of ``__array_ufunc__`` for xarray objects (see below) means that
  NumPy `ufunc methods`_ (e.g., ``np.add.reduce``) that previously worked on
  ``xarray.DataArray`` objects by converting them into NumPy arrays will now
  raise ``NotImplementedError`` instead. In all cases, the work-around is
  simple: convert your objects explicitly into NumPy arrays before calling the
  ufunc (e.g., with ``.values``).

.. _ufunc methods: https://docs.scipy.org/doc/numpy/reference/ufuncs.html#methods

Enhancements
~~~~~~~~~~~~

- Added :py:func:`~xarray.dot`, equivalent to :py:func:`numpy.einsum`.
  Also, :py:func:`~xarray.DataArray.dot` now supports ``dims`` option,
  which specifies the dimensions to sum over.
  (:issue:`1951`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- Support for writing xarray datasets to netCDF files (netcdf4 backend only)
  when using the `dask.distributed <https://distributed.readthedocs.io>`_
  scheduler (:issue:`1464`).
  By `Joe Hamman <https://github.com/jhamman>`_.

- Support lazy vectorized-indexing. After this change, flexible indexing such
  as orthogonal/vectorized indexing, becomes possible for all the backend
  arrays. Also, lazy ``transpose`` is now also supported. (:issue:`1897`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- Implemented NumPy's ``__array_ufunc__`` protocol for all xarray objects
  (:issue:`1617`). This enables using NumPy ufuncs directly on
  ``xarray.Dataset`` objects with recent versions of NumPy (v1.13 and newer):

  .. ipython:: python

      ds = xr.Dataset({"a": 1})
      np.sin(ds)

  This obliviates the need for the ``xarray.ufuncs`` module, which will be
  deprecated in the future when xarray drops support for older versions of
  NumPy. By `Stephan Hoyer <https://github.com/shoyer>`_.

- Improve :py:func:`~xarray.DataArray.rolling` logic.
  :py:func:`~xarray.core.rolling.DataArrayRolling` object now supports
  :py:func:`~xarray.core.rolling.DataArrayRolling.construct` method that returns a view
  of the DataArray / Dataset object with the rolling-window dimension added
  to the last axis. This enables more flexible operation, such as strided
  rolling, windowed rolling, ND-rolling, short-time FFT and convolution.
  (:issue:`1831`, :issue:`1142`, :issue:`819`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- :py:func:`~plot.line()` learned to make plots with data on x-axis if so specified. (:issue:`575`)
  By `Deepak Cherian <https://github.com/dcherian>`_.

Bug fixes
~~~~~~~~~

- Raise an informative error message when using ``apply_ufunc`` with numpy
  v1.11 (:issue:`1956`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- Fix the precision drop after indexing datetime64 arrays (:issue:`1932`).
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- Silenced irrelevant warnings issued by ``open_rasterio`` (:issue:`1964`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- Fix kwarg `colors` clashing with auto-inferred `cmap` (:issue:`1461`)
  By `Deepak Cherian <https://github.com/dcherian>`_.
- Fix :py:func:`~xarray.plot.imshow` error when passed an RGB array with
  size one in a spatial dimension.
  By `Zac Hatfield-Dodds <https://github.com/Zac-HD>`_.

.. _whats-new.0.10.1:

v0.10.1 (25 February 2018)
--------------------------

The minor release includes a number of bug-fixes and backwards compatible enhancements.

Documentation
~~~~~~~~~~~~~

- Added a new guide on :ref:`contributing` (:issue:`640`)
  By `Joe Hamman <https://github.com/jhamman>`_.
- Added apply_ufunc example to :ref:`/examples/weather-data.ipynb#Toy-weather-data` (:issue:`1844`).
  By `Liam Brannigan <https://github.com/braaannigan>`_.
- New entry `Why don’t aggregations return Python scalars?` in the
  :doc:`getting-started-guide/faq` (:issue:`1726`).
  By `0x0L <https://github.com/0x0L>`_.

Enhancements
~~~~~~~~~~~~
**New functions and methods**:

- Added :py:meth:`DataArray.to_iris` and
  :py:meth:`DataArray.from_iris` for
  converting data arrays to and from Iris_ Cubes with the same data and coordinates
  (:issue:`621` and :issue:`37`).
  By `Neil Parley <https://github.com/nparley>`_ and `Duncan Watson-Parris <https://github.com/duncanwp>`_.
- Experimental support for using `Zarr`_ as storage layer for xarray
  (:issue:`1223`).
  By `Ryan Abernathey <https://github.com/rabernat>`_ and
  `Joe Hamman <https://github.com/jhamman>`_.
- New :py:meth:`~xarray.DataArray.rank` on arrays and datasets. Requires
  bottleneck (:issue:`1731`).
  By `0x0L <https://github.com/0x0L>`_.
- ``.dt`` accessor can now ceil, floor and round timestamps to specified frequency.
  By `Deepak Cherian <https://github.com/dcherian>`_.

**Plotting enhancements**:

- :func:`xarray.plot.imshow` now handles RGB and RGBA images.
  Saturation can be adjusted with ``vmin`` and ``vmax``, or with ``robust=True``.
  By `Zac Hatfield-Dodds <https://github.com/Zac-HD>`_.
- :py:func:`~plot.contourf()` learned to contour 2D variables that have both a
  1D coordinate (e.g. time) and a 2D coordinate (e.g. depth as a function of
  time) (:issue:`1737`).
  By `Deepak Cherian <https://github.com/dcherian>`_.
- :py:func:`~plot.plot()` rotates x-axis ticks if x-axis is time.
  By `Deepak Cherian <https://github.com/dcherian>`_.
- :py:func:`~plot.line()` can draw multiple lines if provided with a
  2D variable.
  By `Deepak Cherian <https://github.com/dcherian>`_.

**Other enhancements**:

- Reduce methods such as :py:func:`DataArray.sum()` now handles object-type array.

  .. ipython:: python

      da = xr.DataArray(np.array([True, False, np.nan], dtype=object), dims="x")
      da.sum()

  (:issue:`1866`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- Reduce methods such as :py:func:`DataArray.sum()` now accepts ``dtype``
  arguments. (:issue:`1838`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- Added nodatavals attribute to DataArray when using :py:func:`~xarray.open_rasterio`. (:issue:`1736`).
  By `Alan Snow <https://github.com/snowman2>`_.
- Use ``pandas.Grouper`` class in xarray resample methods rather than the
  deprecated ``pandas.TimeGrouper`` class (:issue:`1766`).
  By `Joe Hamman <https://github.com/jhamman>`_.
- Experimental support for parsing ENVI metadata to coordinates and attributes
  in :py:func:`xarray.open_rasterio`.
  By `Matti Eskelinen <https://github.com/maaleske>`_.
- Reduce memory usage when decoding a variable with a scale_factor, by
  converting 8-bit and 16-bit integers to float32 instead of float64
  (:pull:`1840`), and keeping float16 and float32 as float32 (:issue:`1842`).
  Correspondingly, encoded variables may also be saved with a smaller dtype.
  By `Zac Hatfield-Dodds <https://github.com/Zac-HD>`_.
- Speed of reindexing/alignment with dask array is orders of magnitude faster
  when inserting missing values  (:issue:`1847`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- Fix ``axis`` keyword ignored when applying ``np.squeeze`` to ``DataArray`` (:issue:`1487`).
  By `Florian Pinault <https://github.com/floriankrb>`_.
- ``netcdf4-python`` has moved the its time handling in the ``netcdftime`` module to
  a standalone package (`netcdftime`_). As such, xarray now considers `netcdftime`_
  an optional dependency. One benefit of this change is that it allows for
  encoding/decoding of datetimes with non-standard calendars without the
  ``netcdf4-python`` dependency (:issue:`1084`).
  By `Joe Hamman <https://github.com/jhamman>`_.

.. _Zarr: http://zarr.readthedocs.io/

.. _Iris: http://scitools.org.uk/iris

.. _netcdftime: https://unidata.github.io/netcdftime

**New functions/methods**

- New :py:meth:`~xarray.DataArray.rank` on arrays and datasets. Requires
  bottleneck (:issue:`1731`).
  By `0x0L <https://github.com/0x0L>`_.

Bug fixes
~~~~~~~~~
- Rolling aggregation with ``center=True`` option now gives the same result
  with pandas including the last element (:issue:`1046`).
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- Support indexing with a 0d-np.ndarray (:issue:`1921`).
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- Added warning in api.py of a netCDF4 bug that occurs when
  the filepath has 88 characters (:issue:`1745`).
  By `Liam Brannigan <https://github.com/braaannigan>`_.
- Fixed encoding of multi-dimensional coordinates in
  :py:meth:`~Dataset.to_netcdf` (:issue:`1763`).
  By `Mike Neish <https://github.com/neishm>`_.
- Fixed chunking with non-file-based rasterio datasets (:issue:`1816`) and
  refactored rasterio test suite.
  By `Ryan Abernathey <https://github.com/rabernat>`_
- Bug fix in open_dataset(engine='pydap') (:issue:`1775`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- Bug fix in vectorized assignment  (:issue:`1743`, :issue:`1744`).
  Now item assignment to :py:meth:`~DataArray.__setitem__` checks
- Bug fix in vectorized assignment  (:issue:`1743`, :issue:`1744`).
  Now item assignment to :py:meth:`DataArray.__setitem__` checks
  coordinates of target, destination and keys. If there are any conflict among
  these coordinates, ``IndexError`` will be raised.
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.
- Properly point ``DataArray.__dask_scheduler__`` to
  ``dask.threaded.get``.  By `Matthew Rocklin <https://github.com/mrocklin>`_.
- Bug fixes in :py:meth:`DataArray.plot.imshow`: all-NaN arrays and arrays
  with size one in some dimension can now be plotted, which is good for
  exploring satellite imagery (:issue:`1780`).
  By `Zac Hatfield-Dodds <https://github.com/Zac-HD>`_.
- Fixed ``UnboundLocalError`` when opening netCDF file (:issue:`1781`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- The ``variables``, ``attrs``, and ``dimensions`` properties have been
  deprecated as part of a bug fix addressing an issue where backends were
  unintentionally loading the datastores data and attributes repeatedly during
  writes (:issue:`1798`).
  By `Joe Hamman <https://github.com/jhamman>`_.
- Compatibility fixes to plotting module for NumPy 1.14 and pandas 0.22
  (:issue:`1813`).
  By `Joe Hamman <https://github.com/jhamman>`_.
- Bug fix in encoding coordinates with ``{'_FillValue': None}`` in netCDF
  metadata (:issue:`1865`).
  By `Chris Roth <https://github.com/czr137>`_.
- Fix indexing with lists for arrays loaded from netCDF files with
  ``engine='h5netcdf`` (:issue:`1864`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- Corrected a bug with incorrect coordinates for non-georeferenced geotiff
  files (:issue:`1686`). Internally, we now use the rasterio coordinate
  transform tool instead of doing the computations ourselves. A
  ``parse_coordinates`` kwarg has beed added to :py:func:`~open_rasterio`
  (set to ``True`` per default).
  By `Fabien Maussion <https://github.com/fmaussion>`_.
- The colors of discrete colormaps are now the same regardless if `seaborn`
  is installed or not (:issue:`1896`).
  By `Fabien Maussion <https://github.com/fmaussion>`_.
- Fixed dtype promotion rules in :py:func:`where` and :py:func:`concat` to
  match pandas (:issue:`1847`). A combination of strings/numbers or
  unicode/bytes now promote to object dtype, instead of strings or unicode.
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- Fixed bug where :py:meth:`~xarray.DataArray.isnull` was loading data
  stored as dask arrays (:issue:`1937`).
  By `Joe Hamman <https://github.com/jhamman>`_.

.. _whats-new.0.10.0:

v0.10.0 (20 November 2017)
--------------------------

This is a major release that includes bug fixes, new features and a few
backwards incompatible changes. Highlights include:

- Indexing now supports broadcasting over dimensions, similar to NumPy's
  vectorized indexing (but better!).
- :py:meth:`~DataArray.resample` has a new groupby-like API like pandas.
- :py:func:`~xarray.apply_ufunc` facilitates wrapping and parallelizing
  functions written for NumPy arrays.
- Performance improvements, particularly for dask and :py:func:`open_mfdataset`.

Breaking changes
~~~~~~~~~~~~~~~~

- xarray now supports a form of vectorized indexing with broadcasting, where
  the result of indexing depends on dimensions of indexers,
  e.g., ``array.sel(x=ind)`` with ``ind.dims == ('y',)``. Alignment between
  coordinates on indexed and indexing objects is also now enforced.
  Due to these changes, existing uses of xarray objects to index other xarray
  objects will break in some cases.

  The new indexing API is much more powerful, supporting outer, diagonal and
  vectorized indexing in a single interface.
  The ``isel_points`` and ``sel_points`` methods are deprecated, since they are
  now redundant with the ``isel`` / ``sel`` methods.
  See :ref:`vectorized_indexing` for the details (:issue:`1444`,
  :issue:`1436`).
  By `Keisuke Fujii <https://github.com/fujiisoup>`_ and
  `Stephan Hoyer <https://github.com/shoyer>`_.

- A new resampling interface to match pandas' groupby-like API was added to
  :py:meth:`Dataset.resample` and :py:meth:`DataArray.resample`
  (:issue:`1272`). :ref:`Timeseries resampling <resampling>` is
  fully supported for data with arbitrary dimensions as is both downsampling
  and upsampling (including linear, quadratic, cubic, and spline interpolation).

  Old syntax:

  .. ipython::
    :verbatim:

    In [1]: ds.resample("24H", dim="time", how="max")
    Out[1]:
    <xarray.Dataset>
    [...]

  New syntax:

  .. ipython::
    :verbatim:

    In [1]: ds.resample(time="24H").max()
    Out[1]:
    <xarray.Dataset>
    [...]

  Note that both versions are currently supported, but using the old syntax will
  produce a warning encouraging users to adopt the new syntax.
  By `Daniel Rothenberg <https://github.com/darothen>`_.

- Calling ``repr()`` or printing xarray objects at the command line or in a
  Jupyter Notebook will not longer automatically compute dask variables or
  load data on arrays lazily loaded from disk (:issue:`1522`).
  By `Guido Imperiale <https://github.com/crusaderky>`_.

- Supplying ``coords`` as a dictionary to the ``DataArray`` constructor without
  also supplying an explicit ``dims`` argument is no longer supported. This
  behavior was deprecated in version 0.9 but will now raise an error
  (:issue:`727`).

- Several existing features have been deprecated and will change to new
  behavior in xarray v0.11. If you use any of them with xarray v0.10, you
  should see a ``FutureWarning`` that describes how to update your code:

  - ``Dataset.T`` has been deprecated an alias for ``Dataset.transpose()``
    (:issue:`1232`). In the next major version of xarray, it will provide short-
    cut lookup for variables or attributes with name ``'T'``.
  - ``DataArray.__contains__`` (e.g., ``key in data_array``) currently checks
    for membership in ``DataArray.coords``. In the next major version of
    xarray, it will check membership in the array data found in
    ``DataArray.values`` instead (:issue:`1267`).
  - Direct iteration over and counting a ``Dataset`` (e.g., ``[k for k in ds]``,
    ``ds.keys()``, ``ds.values()``, ``len(ds)`` and ``if ds``) currently
    includes all variables, both data and coordinates. For improved usability
    and consistency with pandas, in the next major version of xarray these will
    change to only include data variables (:issue:`884`). Use ``ds.variables``,
    ``ds.data_vars`` or ``ds.coords`` as alternatives.

- Changes to minimum versions of dependencies:

  - Old numpy < 1.11 and pandas < 0.18 are no longer supported (:issue:`1512`).
    By `Keisuke Fujii <https://github.com/fujiisoup>`_.
  - The minimum supported version bottleneck has increased to 1.1
    (:issue:`1279`).
    By `Joe Hamman <https://github.com/jhamman>`_.

Enhancements
~~~~~~~~~~~~

**New functions/methods**

- New helper function :py:func:`~xarray.apply_ufunc` for wrapping functions
  written to work on NumPy arrays to support labels on xarray objects
  (:issue:`770`). ``apply_ufunc`` also support automatic parallelization for
  many functions with dask. See :ref:`comput.wrapping-custom` and
  :ref:`dask.automatic-parallelization` for details.
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- Added new method :py:meth:`Dataset.to_dask_dataframe`, convert a dataset into
  a dask dataframe.
  This allows lazy loading of data from a dataset containing dask arrays (:issue:`1462`).
  By `James Munroe <https://github.com/jmunroe>`_.

- New function :py:func:`~xarray.where` for conditionally switching between
  values in xarray objects, like :py:func:`numpy.where`:

  .. ipython::
    :verbatim:

    In [1]: import xarray as xr

    In [2]: arr = xr.DataArray([[1, 2, 3], [4, 5, 6]], dims=("x", "y"))

    In [3]: xr.where(arr % 2, "even", "odd")
    Out[3]:
    <xarray.DataArray (x: 2, y: 3)>
    array([['even', 'odd', 'even'],
           ['odd', 'even', 'odd']],
          dtype='<U4')
    Dimensions without coordinates: x, y

  Equivalently, the :py:meth:`~xarray.Dataset.where` method also now supports
  the ``other`` argument, for filling with a value other than ``NaN``
  (:issue:`576`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- Added :py:func:`~xarray.show_versions` function to aid in debugging
  (:issue:`1485`).
  By `Joe Hamman <https://github.com/jhamman>`_.

**Performance improvements**

- :py:func:`~xarray.concat` was computing variables that aren't in memory
  (e.g. dask-based) multiple times; :py:func:`~xarray.open_mfdataset`
  was loading them multiple times from disk. Now, both functions will instead
  load them at most once and, if they do, store them in memory in the
  concatenated array/dataset (:issue:`1521`).
  By `Guido Imperiale <https://github.com/crusaderky>`_.

- Speed-up (x 100) of ``xarray.conventions.decode_cf_datetime``.
  By `Christian Chwala <https://github.com/cchwala>`_.

**IO related improvements**

- Unicode strings (``str`` on Python 3) are now round-tripped successfully even
  when written as character arrays (e.g., as netCDF3 files or when using
  ``engine='scipy'``) (:issue:`1638`). This is controlled by the ``_Encoding``
  attribute convention, which is also understood directly by the netCDF4-Python
  interface. See :ref:`io.string-encoding` for full details.
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- Support for ``data_vars`` and ``coords`` keywords from
  :py:func:`~xarray.concat` added to :py:func:`~xarray.open_mfdataset`
  (:issue:`438`). Using these keyword arguments can significantly reduce
  memory usage and increase speed.
  By `Oleksandr Huziy <https://github.com/guziy>`_.

- Support for :py:class:`pathlib.Path` objects added to
  :py:func:`~xarray.open_dataset`, :py:func:`~xarray.open_mfdataset`,
  ``xarray.to_netcdf``, and :py:func:`~xarray.save_mfdataset`
  (:issue:`799`):

  .. ipython::
    :verbatim:

    In [2]: from pathlib import Path  # In Python 2, use pathlib2!

    In [3]: data_dir = Path("data/")

    In [4]: one_file = data_dir / "dta_for_month_01.nc"

    In [5]: xr.open_dataset(one_file)
    Out[5]:
    <xarray.Dataset>
    [...]

  By `Willi Rath <https://github.com/willirath>`_.

- You can now explicitly disable any default ``_FillValue`` (``NaN`` for
  floating point values) by passing the encoding ``{'_FillValue': None}``
  (:issue:`1598`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- More attributes available in :py:attr:`~xarray.Dataset.attrs` dictionary when
  raster files are opened with :py:func:`~xarray.open_rasterio`.
  By `Greg Brener <https://github.com/gbrener>`_.

- Support for NetCDF files using an ``_Unsigned`` attribute to indicate that a
  a signed integer data type should be interpreted as unsigned bytes
  (:issue:`1444`).
  By `Eric Bruning <https://github.com/deeplycloudy>`_.

- Support using an existing, opened netCDF4 ``Dataset`` with
  :py:class:`~xarray.backends.NetCDF4DataStore`. This permits creating an
  :py:class:`~xarray.Dataset` from a netCDF4 ``Dataset`` that has been opened using
  other means (:issue:`1459`).
  By `Ryan May <https://github.com/dopplershift>`_.

- Changed :py:class:`~xarray.backends.PydapDataStore` to take a Pydap dataset.
  This permits opening Opendap datasets that require authentication, by
  instantiating a Pydap dataset with a session object. Also added
  :py:meth:`xarray.backends.PydapDataStore.open` which takes a url and session
  object (:issue:`1068`).
  By `Philip Graae <https://github.com/mrpgraae>`_.

- Support reading and writing unlimited dimensions with h5netcdf (:issue:`1636`).
  By `Joe Hamman <https://github.com/jhamman>`_.

**Other improvements**

- Added ``_ipython_key_completions_`` to xarray objects, to enable
  autocompletion for dictionary-like access in IPython, e.g.,
  ``ds['tem`` + tab -> ``ds['temperature']`` (:issue:`1628`).
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- Support passing keyword arguments to ``load``, ``compute``, and ``persist``
  methods. Any keyword arguments supplied to these methods are passed on to
  the corresponding dask function (:issue:`1523`).
  By `Joe Hamman <https://github.com/jhamman>`_.

- Encoding attributes are now preserved when xarray objects are concatenated.
  The encoding is copied from the first object  (:issue:`1297`).
  By `Joe Hamman <https://github.com/jhamman>`_ and
  `Gerrit Holl <https://github.com/gerritholl>`_.

- Support applying rolling window operations using bottleneck's moving window
  functions on data stored as dask arrays (:issue:`1279`).
  By `Joe Hamman <https://github.com/jhamman>`_.

- Experimental support for the Dask collection interface (:issue:`1674`).
  By `Matthew Rocklin <https://github.com/mrocklin>`_.

Bug fixes
~~~~~~~~~

- Suppress ``RuntimeWarning`` issued by ``numpy`` for "invalid value comparisons"
  (e.g. ``NaN``). Xarray now behaves similarly to pandas in its treatment of
  binary and unary operations on objects with NaNs (:issue:`1657`).
  By `Joe Hamman <https://github.com/jhamman>`_.

- Unsigned int support for reduce methods with ``skipna=True``
  (:issue:`1562`).
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- Fixes to ensure xarray works properly with pandas 0.21:

  - Fix :py:meth:`~xarray.DataArray.isnull` method (:issue:`1549`).
  - :py:meth:`~xarray.DataArray.to_series` and
    :py:meth:`~xarray.Dataset.to_dataframe` should not return a ``pandas.MultiIndex``
    for 1D data (:issue:`1548`).
  - Fix plotting with datetime64 axis labels (:issue:`1661`).

  By `Stephan Hoyer <https://github.com/shoyer>`_.

- :py:func:`~xarray.open_rasterio` method now shifts the rasterio
  coordinates so that they are centered in each pixel (:issue:`1468`).
  By `Greg Brener <https://github.com/gbrener>`_.

- :py:meth:`~xarray.Dataset.rename` method now doesn't throw errors
  if some ``Variable`` is renamed to the same name as another ``Variable``
  as long as that other ``Variable`` is also renamed (:issue:`1477`). This
  method now does throw when two ``Variables`` would end up with the same name
  after the rename (since one of them would get overwritten in this case).
  By `Prakhar Goel <https://github.com/newt0311>`_.

- Fix :py:func:`xarray.testing.assert_allclose` to actually use ``atol`` and
  ``rtol`` arguments when called on ``DataArray`` objects (:issue:`1488`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- xarray ``quantile`` methods now properly raise a ``TypeError`` when applied to
  objects with data stored as ``dask`` arrays (:issue:`1529`).
  By `Joe Hamman <https://github.com/jhamman>`_.

- Fix positional indexing to allow the use of unsigned integers (:issue:`1405`).
  By `Joe Hamman <https://github.com/jhamman>`_ and
  `Gerrit Holl <https://github.com/gerritholl>`_.

- Creating a :py:class:`Dataset` now raises ``MergeError`` if a coordinate
  shares a name with a dimension but is comprised of arbitrary dimensions
  (:issue:`1120`).
  By `Joe Hamman <https://github.com/jhamman>`_.

- :py:func:`~xarray.open_rasterio` method now skips rasterio's ``crs``
  attribute if its value is ``None`` (:issue:`1520`).
  By `Leevi Annala <https://github.com/leevei>`_.

- Fix :py:func:`xarray.DataArray.to_netcdf` to return bytes when no path is
  provided (:issue:`1410`).
  By `Joe Hamman <https://github.com/jhamman>`_.

- Fix :py:func:`xarray.save_mfdataset` to properly raise an informative error
  when objects other than  ``Dataset`` are provided (:issue:`1555`).
  By `Joe Hamman <https://github.com/jhamman>`_.

- :py:func:`xarray.Dataset.copy` would not preserve the encoding property
  (:issue:`1586`).
  By `Guido Imperiale <https://github.com/crusaderky>`_.

- :py:func:`xarray.concat` would eagerly load dask variables into memory if
  the first argument was a numpy variable (:issue:`1588`).
  By `Guido Imperiale <https://github.com/crusaderky>`_.

- Fix bug in :py:meth:`~xarray.Dataset.to_netcdf` when writing in append mode
  (:issue:`1215`).
  By `Joe Hamman <https://github.com/jhamman>`_.

- Fix ``netCDF4`` backend to properly roundtrip the ``shuffle`` encoding option
  (:issue:`1606`).
  By `Joe Hamman <https://github.com/jhamman>`_.

- Fix bug when using ``pytest`` class decorators to skiping certain unittests.
  The previous behavior unintentionally causing additional tests to be skipped
  (:issue:`1531`). By `Joe Hamman <https://github.com/jhamman>`_.

- Fix pynio backend for upcoming release of pynio with Python 3 support
  (:issue:`1611`). By `Ben Hillman <https://github/brhillman>`_.

- Fix ``seaborn`` import warning for Seaborn versions 0.8 and newer when the
  ``apionly`` module was deprecated.
  (:issue:`1633`). By `Joe Hamman <https://github.com/jhamman>`_.

- Fix COMPAT: MultiIndex checking is fragile
  (:issue:`1833`). By `Florian Pinault <https://github.com/floriankrb>`_.

- Fix ``rasterio`` backend for Rasterio versions 1.0alpha10 and newer.
  (:issue:`1641`). By `Chris Holden <https://github.com/ceholden>`_.

Bug fixes after rc1
~~~~~~~~~~~~~~~~~~~

- Suppress warning in IPython autocompletion, related to the deprecation
  of ``.T`` attributes (:issue:`1675`).
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- Fix a bug in lazily-indexing netCDF array. (:issue:`1688`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- (Internal bug) MemoryCachedArray now supports the orthogonal indexing.
  Also made some internal cleanups around array wrappers (:issue:`1429`).
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- (Internal bug) MemoryCachedArray now always wraps ``np.ndarray`` by
  ``NumpyIndexingAdapter``. (:issue:`1694`)
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- Fix importing xarray when running Python with ``-OO`` (:issue:`1706`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- Saving a netCDF file with a coordinates with a spaces in its names now raises
  an appropriate warning (:issue:`1689`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- Fix two bugs that were preventing dask arrays from being specified as
  coordinates in the DataArray constructor (:issue:`1684`).
  By `Joe Hamman <https://github.com/jhamman>`_.

- Fixed ``apply_ufunc`` with ``dask='parallelized'`` for scalar arguments
  (:issue:`1697`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- Fix "Chunksize cannot exceed dimension size" error when writing netCDF4 files
  loaded from disk (:issue:`1225`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- Validate the shape of coordinates with names matching dimensions in the
  DataArray constructor (:issue:`1709`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- Raise ``NotImplementedError`` when attempting to save a MultiIndex to a
  netCDF file (:issue:`1547`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- Remove netCDF dependency from rasterio backend tests.
  By `Matti Eskelinen <https://github.com/maaleske>`_

Bug fixes after rc2
~~~~~~~~~~~~~~~~~~~

- Fixed unexpected behavior in ``Dataset.set_index()`` and
  ``DataArray.set_index()`` introduced by pandas 0.21.0. Setting a new
  index with a single variable resulted in 1-level
  ``pandas.MultiIndex`` instead of a simple ``pandas.Index``
  (:issue:`1722`).  By `Benoit Bovy <https://github.com/benbovy>`_.

- Fixed unexpected memory loading of backend arrays after ``print``.
  (:issue:`1720`).  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

.. _whats-new.0.9.6:

v0.9.6 (8 June 2017)
--------------------

This release includes a number of backwards compatible enhancements and bug
fixes.

Enhancements
~~~~~~~~~~~~

- New :py:meth:`~xarray.Dataset.sortby` method to ``Dataset`` and ``DataArray``
  that enable sorting along dimensions (:issue:`967`).
  See :ref:`the docs <reshape.sort>` for examples.
  By `Chun-Wei Yuan <https://github.com/chunweiyuan>`_ and
  `Kyle Heuton <https://github.com/kheuton>`_.

- Add ``.dt`` accessor to DataArrays for computing datetime-like properties
  for the values they contain, similar to ``pandas.Series`` (:issue:`358`).
  By `Daniel Rothenberg <https://github.com/darothen>`_.

- Renamed internal dask arrays created by ``open_dataset`` to match new dask
  conventions (:issue:`1343`).
  By `Ryan Abernathey <https://github.com/rabernat>`_.

- :py:meth:`~xarray.as_variable` is now part of the public API (:issue:`1303`).
  By `Benoit Bovy <https://github.com/benbovy>`_.

- :py:func:`~xarray.align` now supports ``join='exact'``, which raises
  an error instead of aligning when indexes to be aligned are not equal.
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- New function :py:func:`~xarray.open_rasterio` for opening raster files with
  the `rasterio <https://rasterio.readthedocs.io/en/latest/>`_ library.
  See :ref:`the docs <io.rasterio>` for details.
  By `Joe Hamman <https://github.com/jhamman>`_,
  `Nic Wayand <https://github.com/NicWayand>`_ and
  `Fabien Maussion <https://github.com/fmaussion>`_

Bug fixes
~~~~~~~~~

- Fix error from repeated indexing of datasets loaded from disk (:issue:`1374`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- Fix a bug where ``.isel_points`` wrongly assigns unselected coordinate to
  ``data_vars``.
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- Tutorial datasets are now checked against a reference MD5 sum to confirm
  successful download (:issue:`1392`). By `Matthew Gidden
  <https://github.com/gidden>`_.

- ``DataArray.chunk()`` now accepts dask specific kwargs like
  ``Dataset.chunk()`` does. By `Fabien Maussion <https://github.com/fmaussion>`_.

- Support for ``engine='pydap'`` with recent releases of Pydap (3.2.2+),
  including on Python 3 (:issue:`1174`).

Documentation
~~~~~~~~~~~~~

- A new `gallery <http://xarray.pydata.org/en/latest/auto_gallery/index.html>`_
  allows to add interactive examples to the documentation.
  By `Fabien Maussion <https://github.com/fmaussion>`_.

Testing
~~~~~~~

- Fix test suite failure caused by changes to ``pandas.cut`` function
  (:issue:`1386`).
  By `Ryan Abernathey <https://github.com/rabernat>`_.

- Enhanced tests suite by use of ``@network`` decorator, which is
  controlled via ``--run-network-tests`` command line argument
  to ``py.test`` (:issue:`1393`).
  By `Matthew Gidden <https://github.com/gidden>`_.

.. _whats-new.0.9.5:

v0.9.5 (17 April, 2017)
-----------------------

Remove an inadvertently introduced print statement.

.. _whats-new.0.9.3:

v0.9.3 (16 April, 2017)
-----------------------

This minor release includes bug-fixes and backwards compatible enhancements.

Enhancements
~~~~~~~~~~~~

- New :py:meth:`~xarray.DataArray.persist` method to Datasets and DataArrays to
  enable persisting data in distributed memory when using Dask (:issue:`1344`).
  By `Matthew Rocklin <https://github.com/mrocklin>`_.

- New :py:meth:`~xarray.DataArray.expand_dims` method for ``DataArray`` and
  ``Dataset`` (:issue:`1326`).
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

Bug fixes
~~~~~~~~~

- Fix ``.where()`` with ``drop=True`` when arguments do not have indexes
  (:issue:`1350`). This bug, introduced in v0.9, resulted in xarray producing
  incorrect results in some cases.
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- Fixed writing to file-like objects with :py:meth:`~xarray.Dataset.to_netcdf`
  (:issue:`1320`).
  `Stephan Hoyer <https://github.com/shoyer>`_.

- Fixed explicitly setting ``engine='scipy'`` with ``to_netcdf`` when not
  providing a path (:issue:`1321`).
  `Stephan Hoyer <https://github.com/shoyer>`_.

- Fixed open_dataarray does not pass properly its parameters to open_dataset
  (:issue:`1359`).
  `Stephan Hoyer <https://github.com/shoyer>`_.

- Ensure test suite works when runs from an installed version of xarray
  (:issue:`1336`). Use ``@pytest.mark.slow`` instead of a custom flag to mark
  slow tests.
  By `Stephan Hoyer <https://github.com/shoyer>`_

.. _whats-new.0.9.2:

v0.9.2 (2 April 2017)
---------------------

The minor release includes bug-fixes and backwards compatible enhancements.

Enhancements
~~~~~~~~~~~~

- ``rolling`` on Dataset is now supported (:issue:`859`).

- ``.rolling()`` on Dataset is now supported (:issue:`859`).
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- When bottleneck version 1.1 or later is installed, use bottleneck for rolling
  ``var``, ``argmin``, ``argmax``, and ``rank`` computations. Also, rolling
  median now accepts a ``min_periods`` argument (:issue:`1276`).
  By `Joe Hamman <https://github.com/jhamman>`_.

- When ``.plot()`` is called on a 2D DataArray and only one dimension is
  specified with ``x=`` or ``y=``, the other dimension is now guessed
  (:issue:`1291`).
  By `Vincent Noel <https://github.com/vnoel>`_.

- Added new method :py:meth:`~Dataset.assign_attrs` to ``DataArray`` and
  ``Dataset``, a chained-method compatible implementation of the
  ``dict.update`` method on attrs (:issue:`1281`).
  By `Henry S. Harrison <https://hsharrison.github.io>`_.

- Added new ``autoclose=True`` argument to
  :py:func:`~xarray.open_mfdataset` to explicitly close opened files when not in
  use to prevent occurrence of an OS Error related to too many open files
  (:issue:`1198`).
  Note, the default is ``autoclose=False``, which is consistent with
  previous xarray behavior.
  By `Phillip J. Wolfram <https://github.com/pwolfram>`_.

- The ``repr()`` of ``Dataset`` and ``DataArray`` attributes uses a similar
  format to coordinates and variables, with vertically aligned entries
  truncated to fit on a single line (:issue:`1319`).  Hopefully this will stop
  people writing ``data.attrs = {}`` and discarding metadata in notebooks for
  the sake of cleaner output.  The full metadata is still available as
  ``data.attrs``.
  By `Zac Hatfield-Dodds <https://github.com/Zac-HD>`_.

- Enhanced tests suite by use of ``@slow`` and ``@flaky`` decorators, which are
  controlled via ``--run-flaky`` and ``--skip-slow`` command line arguments
  to ``py.test`` (:issue:`1336`).
  By `Stephan Hoyer <https://github.com/shoyer>`_ and
  `Phillip J. Wolfram <https://github.com/pwolfram>`_.

- New aggregation on rolling objects :py:meth:`~core.rolling.DataArrayRolling.count`
  which providing a rolling count of valid values (:issue:`1138`).

Bug fixes
~~~~~~~~~
- Rolling operations now keep preserve original dimension order (:issue:`1125`).
  By `Keisuke Fujii <https://github.com/fujiisoup>`_.

- Fixed ``sel`` with ``method='nearest'`` on Python 2.7 and 64-bit Windows
  (:issue:`1140`).
  `Stephan Hoyer <https://github.com/shoyer>`_.

- Fixed ``where`` with ``drop='True'`` for empty masks (:issue:`1341`).
  By `Stephan Hoyer <https://github.com/shoyer>`_ and
  `Phillip J. Wolfram <https://github.com/pwolfram>`_.

.. _whats-new.0.9.1:

v0.9.1 (30 January 2017)
------------------------

Renamed the "Unindexed dimensions" section in the ``Dataset`` and
``DataArray`` repr (added in v0.9.0) to "Dimensions without coordinates"
(:issue:`1199`).

.. _whats-new.0.9.0:

v0.9.0 (25 January 2017)
------------------------

This major release includes five months worth of enhancements and bug fixes from
24 contributors, including some significant changes that are not fully backwards
compatible. Highlights include:

- Coordinates are now *optional* in the xarray data model, even for dimensions.
- Changes to caching, lazy loading and pickling to improve xarray's experience
  for parallel computing.
- Improvements for accessing and manipulating ``pandas.MultiIndex`` levels.
- Many new methods and functions, including
  :py:meth:`~DataArray.quantile`,
  :py:meth:`~DataArray.cumsum`,
  :py:meth:`~DataArray.cumprod`
  :py:attr:`~DataArray.combine_first`
  :py:meth:`~DataArray.set_index`,
  :py:meth:`~DataArray.reset_index`,
  :py:meth:`~DataArray.reorder_levels`,
  :py:func:`~xarray.full_like`,
  :py:func:`~xarray.zeros_like`,
  :py:func:`~xarray.ones_like`
  :py:func:`~xarray.open_dataarray`,
  :py:meth:`~DataArray.compute`,
  :py:meth:`Dataset.info`,
  :py:func:`testing.assert_equal`,
  :py:func:`testing.assert_identical`, and
  :py:func:`testing.assert_allclose`.

Breaking changes
~~~~~~~~~~~~~~~~

- Index coordinates for each dimensions are now optional, and no longer created
  by default :issue:`1017`. You can identify such dimensions without coordinates
  by their appearance in list of "Dimensions without coordinates" in the
  ``Dataset`` or ``DataArray`` repr:

  .. ipython::
    :verbatim:

    In [1]: xr.Dataset({"foo": (("x", "y"), [[1, 2]])})
    Out[1]:
    <xarray.Dataset>
    Dimensions:  (x: 1, y: 2)
    Dimensions without coordinates: x, y
    Data variables:
        foo      (x, y) int64 1 2

  This has a number of implications:

  - :py:func:`~align` and :py:meth:`~Dataset.reindex` can now error, if
    dimensions labels are missing and dimensions have different sizes.
  - Because pandas does not support missing indexes, methods such as
    ``to_dataframe``/``from_dataframe`` and ``stack``/``unstack`` no longer
    roundtrip faithfully on all inputs. Use :py:meth:`~Dataset.reset_index` to
    remove undesired indexes.
  - ``Dataset.__delitem__`` and :py:meth:`~Dataset.drop` no longer delete/drop
    variables that have dimensions matching a deleted/dropped variable.
  - ``DataArray.coords.__delitem__`` is now allowed on variables matching
    dimension names.
  - ``.sel`` and ``.loc`` now handle indexing along a dimension without
    coordinate labels by doing integer based indexing. See
    :ref:`indexing.missing_coordinates` for an example.
  - :py:attr:`~Dataset.indexes` is no longer guaranteed to include all
    dimensions names as keys. The new method :py:meth:`~Dataset.get_index` has
    been added to get an index for a dimension guaranteed, falling back to
    produce a default ``RangeIndex`` if necessary.

- The default behavior of ``merge`` is now ``compat='no_conflicts'``, so some
  merges will now succeed in cases that previously raised
  ``xarray.MergeError``. Set ``compat='broadcast_equals'`` to restore the
  previous default. See :ref:`combining.no_conflicts` for more details.

- Reading :py:attr:`~DataArray.values` no longer always caches values in a NumPy
  array :issue:`1128`. Caching of ``.values`` on variables read from netCDF
  files on disk is still the default when :py:func:`open_dataset` is called with
  ``cache=True``.
  By `Guido Imperiale <https://github.com/crusaderky>`_ and
  `Stephan Hoyer <https://github.com/shoyer>`_.
- Pickling a ``Dataset`` or ``DataArray`` linked to a file on disk no longer
  caches its values into memory before pickling (:issue:`1128`). Instead, pickle
  stores file paths and restores objects by reopening file references. This
  enables preliminary, experimental use of xarray for opening files with
  `dask.distributed <https://distributed.readthedocs.io>`_.
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- Coordinates used to index a dimension are now loaded eagerly into
  :py:class:`pandas.Index` objects, instead of loading the values lazily.
  By `Guido Imperiale <https://github.com/crusaderky>`_.
- Automatic levels for 2d plots are now guaranteed to land on ``vmin`` and
  ``vmax`` when these kwargs are explicitly provided (:issue:`1191`). The
  automated level selection logic also slightly changed.
  By `Fabien Maussion <https://github.com/fmaussion>`_.

- ``DataArray.rename()`` behavior changed to strictly change the ``DataArray.name``
  if called with string argument, or strictly change coordinate names if called with
  dict-like argument.
  By `Markus Gonser <https://github.com/magonser>`_.

- By default ``to_netcdf()`` add a ``_FillValue = NaN`` attributes to float types.
  By `Frederic Laliberte <https://github.com/laliberte>`_.

- ``repr`` on ``DataArray`` objects uses an shortened display for NumPy array
  data that is less likely to overflow onto multiple pages (:issue:`1207`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- xarray no longer supports python 3.3, versions of dask prior to v0.9.0,
  or versions of bottleneck prior to v1.0.

Deprecations
~~~~~~~~~~~~

- Renamed the ``Coordinate`` class from xarray's low level API to
  :py:class:`~xarray.IndexVariable`. ``Variable.to_variable`` and
  ``Variable.to_coord`` have been renamed to
  :py:meth:`~xarray.Variable.to_base_variable` and
  :py:meth:`~xarray.Variable.to_index_variable`.
- Deprecated supplying ``coords`` as a dictionary to the ``DataArray``
  constructor without also supplying an explicit ``dims`` argument. The old
  behavior encouraged relying on the iteration order of dictionaries, which is
  a bad practice (:issue:`727`).
- Removed a number of methods deprecated since v0.7.0 or earlier:
  ``load_data``, ``vars``, ``drop_vars``, ``dump``, ``dumps`` and the
  ``variables`` keyword argument to ``Dataset``.
- Removed the dummy module that enabled ``import xray``.

Enhancements
~~~~~~~~~~~~

- Added new method :py:meth:`~DataArray.combine_first` to ``DataArray`` and
  ``Dataset``, based on the pandas method of the same name (see :ref:`combine`).
  By `Chun-Wei Yuan <https://github.com/chunweiyuan>`_.

- Added the ability to change default automatic alignment (arithmetic_join="inner")
  for binary operations via :py:func:`~xarray.set_options()`
  (see :ref:`math automatic alignment`).
  By `Chun-Wei Yuan <https://github.com/chunweiyuan>`_.

- Add checking of ``attr`` names and values when saving to netCDF, raising useful
  error messages if they are invalid. (:issue:`911`).
  By `Robin Wilson <https://github.com/robintw>`_.
- Added ability to save ``DataArray`` objects directly to netCDF files using
  :py:meth:`~xarray.DataArray.to_netcdf`, and to load directly from netCDF files
  using :py:func:`~xarray.open_dataarray` (:issue:`915`). These remove the need
  to convert a ``DataArray`` to a ``Dataset`` before saving as a netCDF file,
  and deals with names to ensure a perfect 'roundtrip' capability.
  By `Robin Wilson <https://github.com/robintw>`_.
- Multi-index levels are now accessible as "virtual" coordinate variables,
  e.g., ``ds['time']`` can pull out the ``'time'`` level of a multi-index
  (see :ref:`coordinates`). ``sel`` also accepts providing multi-index levels
  as keyword arguments, e.g., ``ds.sel(time='2000-01')``
  (see :ref:`multi-level indexing`).
  By `Benoit Bovy <https://github.com/benbovy>`_.
- Added ``set_index``, ``reset_index`` and ``reorder_levels`` methods to
  easily create and manipulate (multi-)indexes (see :ref:`reshape.set_index`).
  By `Benoit Bovy <https://github.com/benbovy>`_.
- Added the ``compat`` option ``'no_conflicts'`` to ``merge``, allowing the
  combination of xarray objects with disjoint (:issue:`742`) or
  overlapping (:issue:`835`) coordinates as long as all present data agrees.
  By `Johnnie Gray <https://github.com/jcmgray>`_. See
  :ref:`combining.no_conflicts` for more details.
- It is now possible to set ``concat_dim=None`` explicitly in
  :py:func:`~xarray.open_mfdataset` to disable inferring a dimension along
  which to concatenate.
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- Added methods :py:meth:`DataArray.compute`, :py:meth:`Dataset.compute`, and
  :py:meth:`Variable.compute` as a non-mutating alternative to
  :py:meth:`~DataArray.load`.
  By `Guido Imperiale <https://github.com/crusaderky>`_.
- Adds DataArray and Dataset methods :py:meth:`~xarray.DataArray.cumsum` and
  :py:meth:`~xarray.DataArray.cumprod`.  By `Phillip J. Wolfram
  <https://github.com/pwolfram>`_.

- New properties :py:attr:`Dataset.sizes` and :py:attr:`DataArray.sizes` for
  providing consistent access to dimension length on both ``Dataset`` and
  ``DataArray`` (:issue:`921`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- New keyword argument ``drop=True`` for :py:meth:`~DataArray.sel`,
  :py:meth:`~DataArray.isel` and :py:meth:`~DataArray.squeeze` for dropping
  scalar coordinates that arise from indexing.
  ``DataArray`` (:issue:`242`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- New top-level functions :py:func:`~xarray.full_like`,
  :py:func:`~xarray.zeros_like`, and :py:func:`~xarray.ones_like`
  By `Guido Imperiale <https://github.com/crusaderky>`_.
- Overriding a preexisting attribute with
  :py:func:`~xarray.register_dataset_accessor` or
  :py:func:`~xarray.register_dataarray_accessor` now issues a warning instead of
  raising an error (:issue:`1082`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- Options for axes sharing between subplots are exposed to
  :py:class:`~xarray.plot.FacetGrid` and :py:func:`~xarray.plot.plot`, so axes
  sharing can be disabled for polar plots.
  By `Bas Hoonhout <https://github.com/hoonhout>`_.
- New utility functions :py:func:`~xarray.testing.assert_equal`,
  :py:func:`~xarray.testing.assert_identical`, and
  :py:func:`~xarray.testing.assert_allclose` for asserting relationships
  between xarray objects, designed for use in a pytest test suite.
- ``figsize``, ``size`` and ``aspect`` plot arguments are now supported for all
  plots (:issue:`897`). See :ref:`plotting.figsize` for more details.
  By `Stephan Hoyer <https://github.com/shoyer>`_ and
  `Fabien Maussion <https://github.com/fmaussion>`_.
- New :py:meth:`~Dataset.info` method to summarize ``Dataset`` variables
  and attributes. The method prints to a buffer (e.g. ``stdout``) with output
  similar to what the command line utility ``ncdump -h`` produces (:issue:`1150`).
  By `Joe Hamman <https://github.com/jhamman>`_.
- Added the ability write unlimited netCDF dimensions with the ``scipy`` and
  ``netcdf4`` backends via the new ``xray.Dataset.encoding`` attribute
  or via the ``unlimited_dims`` argument to ``xray.Dataset.to_netcdf``.
  By `Joe Hamman <https://github.com/jhamman>`_.
- New :py:meth:`~DataArray.quantile` method to calculate quantiles from
  DataArray objects (:issue:`1187`).
  By `Joe Hamman <https://github.com/jhamman>`_.


Bug fixes
~~~~~~~~~
- ``groupby_bins`` now restores empty bins by default (:issue:`1019`).
  By `Ryan Abernathey <https://github.com/rabernat>`_.

- Fix issues for dates outside the valid range of pandas timestamps
  (:issue:`975`). By `Mathias Hauser <https://github.com/mathause>`_.

- Unstacking produced flipped array after stacking decreasing coordinate values
  (:issue:`980`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- Setting ``dtype`` via the ``encoding`` parameter of ``to_netcdf`` failed if
  the encoded dtype was the same as the dtype of the original array
  (:issue:`873`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- Fix issues with variables where both attributes ``_FillValue`` and
  ``missing_value`` are set to ``NaN`` (:issue:`997`).
  By `Marco Zühlke <https://github.com/mzuehlke>`_.

- ``.where()`` and ``.fillna()`` now preserve attributes (:issue:`1009`).
  By `Fabien Maussion <https://github.com/fmaussion>`_.

- Applying :py:func:`broadcast()` to an xarray object based on the dask backend
  won't accidentally convert the array from dask to numpy anymore (:issue:`978`).
  By `Guido Imperiale <https://github.com/crusaderky>`_.

- ``Dataset.concat()`` now preserves variables order (:issue:`1027`).
  By `Fabien Maussion <https://github.com/fmaussion>`_.

- Fixed an issue with pcolormesh (:issue:`781`). A new
  ``infer_intervals`` keyword gives control on whether the cell intervals
  should be computed or not.
  By `Fabien Maussion <https://github.com/fmaussion>`_.

- Grouping over an dimension with non-unique values with ``groupby`` gives
  correct groups.
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- Fixed accessing coordinate variables with non-string names from ``.coords``.
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- :py:meth:`~xarray.DataArray.rename` now simultaneously renames the array and
  any coordinate with the same name, when supplied via a :py:class:`dict`
  (:issue:`1116`).
  By `Yves Delley <https://github.com/burnpanck>`_.

- Fixed sub-optimal performance in certain operations with object arrays (:issue:`1121`).
  By `Yves Delley <https://github.com/burnpanck>`_.

- Fix ``.groupby(group)`` when ``group`` has datetime dtype (:issue:`1132`).
  By `Jonas Sølvsteen <https://github.com/j08lue>`_.

- Fixed a bug with facetgrid (the ``norm`` keyword was ignored, :issue:`1159`).
  By `Fabien Maussion <https://github.com/fmaussion>`_.

- Resolved a concurrency bug that could cause Python to crash when
  simultaneously reading and writing netCDF4 files with dask (:issue:`1172`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- Fix to make ``.copy()`` actually copy dask arrays, which will be relevant for
  future releases of dask in which dask arrays will be mutable (:issue:`1180`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

- Fix opening NetCDF files with multi-dimensional time variables
  (:issue:`1229`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.

Performance improvements
~~~~~~~~~~~~~~~~~~~~~~~~

- ``xarray.Dataset.isel_points`` and ``xarray.Dataset.sel_points`` now
  use vectorised indexing in numpy and dask (:issue:`1161`), which can
  result in several orders of magnitude speedup.
  By `Jonathan Chambers <https://github.com/mangecoeur>`_.

.. _whats-new.0.8.2:

v0.8.2 (18 August 2016)
-----------------------

This release includes a number of bug fixes and minor enhancements.

Breaking changes
~~~~~~~~~~~~~~~~

- :py:func:`~xarray.broadcast` and :py:func:`~xarray.concat` now auto-align
  inputs, using ``join=outer``. Previously, these functions raised
  ``ValueError`` for non-aligned inputs.
  By `Guido Imperiale <https://github.com/crusaderky>`_.

Enhancements
~~~~~~~~~~~~

- New documentation on :ref:`panel transition`. By
  `Maximilian Roos <https://github.com/max-sixty>`_.
- New ``Dataset`` and ``DataArray`` methods :py:meth:`~xarray.Dataset.to_dict`
  and :py:meth:`~xarray.Dataset.from_dict` to allow easy conversion between
  dictionaries and xarray objects (:issue:`432`). See
  :ref:`dictionary IO<dictionary io>` for more details.
  By `Julia Signell <https://github.com/jsignell>`_.
- Added ``exclude`` and ``indexes`` optional parameters to :py:func:`~xarray.align`,
  and ``exclude`` optional parameter to :py:func:`~xarray.broadcast`.
  By `Guido Imperiale <https://github.com/crusaderky>`_.
- Better error message when assigning variables without dimensions
  (:issue:`971`). By `Stephan Hoyer <https://github.com/shoyer>`_.
- Better error message when reindex/align fails due to duplicate index values
  (:issue:`956`). By `Stephan Hoyer <https://github.com/shoyer>`_.

Bug fixes
~~~~~~~~~

- Ensure xarray works with h5netcdf v0.3.0 for arrays with ``dtype=str``
  (:issue:`953`). By `Stephan Hoyer <https://github.com/shoyer>`_.
- ``Dataset.__dir__()`` (i.e. the method python calls to get autocomplete
  options) failed if one of the dataset's keys was not a string (:issue:`852`).
  By `Maximilian Roos <https://github.com/max-sixty>`_.
- ``Dataset`` constructor can now take arbitrary objects as values
  (:issue:`647`). By `Maximilian Roos <https://github.com/max-sixty>`_.
- Clarified ``copy`` argument for :py:meth:`~xarray.DataArray.reindex` and
  :py:func:`~xarray.align`, which now consistently always return new xarray
  objects (:issue:`927`).
- Fix ``open_mfdataset`` with ``engine='pynio'`` (:issue:`936`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- ``groupby_bins`` sorted bin labels as strings (:issue:`952`).
  By `Stephan Hoyer <https://github.com/shoyer>`_.
- Fix bug introduced by v0.8.0 that broke assignment to datasets when both the
  left and right side have the same non-unique index values (:issue:`956`).

.. _whats-new.0.8.1:

v0.8.1 (5 August 2016)
----------------------

Bug fixes
~~~~~~~~~

- Fix bug in v0.8.0 that broke assignment to Datasets with non-unique
  indexes (:issue:`943`). By `Stephan Hoyer <https://github.com/shoyer>`_.

.. _whats-new.0.8.0:

v0.8.0 (2 August 2016)
----------------------

This release includes four months of new features and bug fixes, including
several breaking changes.

.. _v0.8.0.breaking:

Breaking changes
~~~~~~~~~~~~~~~~

- Dropped support for Python 2.6 (:issue:`855`).
- Indexing on multi-index now drop levels, which is consistent with pandas.
  It also changes the name of the dimension / coordinate when the multi-index is
  reduced to a single index (:issue:`802`).
- Contour plots no longer add a colorbar per default (:issue:`866`). Filled
  contour plots are unchanged.
- ``DataArray.values`` and ``.data`` now always returns an NumPy array-like
  object, even for 0-dimensional arrays with object dtype (:issue:`867`).
  Previously, ``.values`` returned native Python objects in such cases. To
  convert the values of scalar arrays to Python objects, use the ``.item()``
  method.

Enhancements
~~~~~~~~~~~~

- Groupby operations now support grouping over multidimensional variables. A new
  method called :py:meth:`~xarray.Dataset.groupby_bins` has also been added to
  allow users to specify bins for grouping. The new features are described in
  :ref:`groupby.multidim` and :ref:`/examples/multidimensional-coords.ipynb`.
  By `Ryan Abernathey <https://github.com/rabernat>`_.

- DataArray and Dataset method :py:meth:`where` now supports a ``drop=True``
  option that clips coordinate elements that are fully masked.  By
  `Phillip J. Wolfram <https://github.com/pwolfram>`_.

- New top level :py:func:`merge` function allows for combining variables from
  any number of ``Dataset`` and/or ``DataArray`` variables. See :ref:`merge`
  for more details. By `Stephan Hoyer <https://github.com/shoyer>`_.

- :py:meth:`DataArray.resample` and :py:meth:`Dataset.resample` now support the
  ``keep_attrs=False`` option that determines whether variable and dataset
  attributes are retained in the resampled object. By
  `Jeremy McGibbon <https://github.com/mcgibbon>`_.

- Better multi-index support in :py:meth:`DataArray.sel`,
  :py:meth:`DataArray.loc`, :py:meth:`Dataset.sel` and
  :py:meth:`Dataset.loc`, which now behave more closely to pandas and
  which also accept dictionaries for indexing based on given level names
  and labels (see :ref:`multi-level indexing`).
  By `Benoit Bovy <https://github.com/benbovy>`_.

- New (experimental) decorators :py:func:`~xarray.register_dataset_accessor` and
  :py:func:`~xarray.register_dataarray_accessor` for registering custom xarray
  extensions without subclassing. They are described in the new documentation
  page on :ref:`internals`. By `Stephan Hoyer <https://github.com/shoyer>`_.

- Round trip boolean datatypes. Previously, writing boolean datatypes to netCDF
  formats would raise an error since netCDF does not have a `bool` datatype.
  This feature reads/writes a `dtype` attribute to boolean variables in netCDF
  files. By `Joe Hamman <https://github.com/jhamman>`_.

- 2D plotting methods now have two new keywords (`cbar_ax` and `cbar_kwargs`),
  allowing more control on the colorbar (:issue:`872`).
  By `Fabien Maussion <https://github.com/fmaussion>`_.

- New Dataset method :py:meth:`Dataset.filter_by_attrs`, akin to
  ``netCDF4.Dataset.get_variables_by_attributes``, to easily filter
  data variables using its attributes.
  `Filipe Fernandes <https://github.com/ocefpaf>`_.

Bug fixes
~~~~~~~~~

- Attributes were being retained by default for some resampling
  operations when they should not. With the ``keep_attrs=False`` option, they
  will no longer be retained by default. This may be backwards-incompatible
  with some scripts, but the attributes may be kept by adding the
  ``keep_attrs=True`` option. By
  `Jeremy McGibbon <https://github.com/mcgibbon>`_.

- Concatenating xarray objects along an axis with a MultiIndex or PeriodIndex
  preserves the nature of the index (:issue:`875`). By
  `Stephan Hoyer <https://github.com/shoyer>`_.

- Fixed bug in arithmetic operations on DataArray objects whose dimensions
  are numpy structured arrays or recarrays :issue:`861`, :issue:`837`. By
  `Maciek Swat <https://github.com/maciekswat>`_.

- ``decode_cf_timedelta`` now accepts arrays with ``ndim`` >1 (:issue:`842`).
   This fixes issue :issue:`665`.
   `Filipe Fernandes <https://github.com/ocefpaf>`_.

- Fix a bug where `xarray.ufuncs` that take two arguments would incorrectly
  use to numpy functions instead of dask.array functions (:issue:`876`). By
  `Stephan Hoyer <https://github.com/shoyer>`_.

- Support for pickling functions from  ``xarray.ufuncs`` (:issue:`901`). By
  `Stephan Hoyer <https://github.com/shoyer>`_.

- ``Variable.copy(deep=True)`` no longer converts MultiIndex into a base Index
  (:issue:`769`). By `Benoit Bovy <https://github.com/benbovy>`_.

- Fixes for groupby on dimensions with a multi-index (:issue:`867`). By
  `Stephan Hoyer <https://github.com/shoyer>`_.

- Fix printing datasets with unicode attributes on Python 2 (:issue:`892`). By
  `Stephan Hoyer <https://github.com/shoyer>`_.

- Fixed incorrect test for dask version (:issue:`891`). By
  `Stephan Hoyer <https://github.com/shoyer>`_.

- Fixed `dim` argument for `isel_points`/`sel_points` when a `pandas.Index` is
  passed. By `Stephan Hoyer <https://github.com/shoyer>`_.

- :py:func:`~xarray.plot.contour` now plots the correct number of contours
  (:issue:`866`). By `Fabien Maussion <https://github.com/fmaussion>`_.

.. _whats-new.0.7.2:

v0.7.2 (13 March 2016)
----------------------

This release includes two new, entirely backwards compatible features and
several bug fixes.

Enhancements
~~~~~~~~~~~~

- New DataArray method :py:meth:`DataArray.dot` for calculating the dot
  product of two DataArrays along shared dimensions. By
  `Dean Pospisil <https://github.com/deanpospisil>`_.

- Rolling window operations on DataArray objects are now supported via a new
  :py:meth:`DataArray.rolling` method. For example:

  .. ipython::
    :verbatim:

    In [1]: import xarray as xr
       ...: import numpy as np

    In [2]: arr = xr.DataArray(np.arange(0, 7.5, 0.5).reshape(3, 5), dims=("x", "y"))

    In [3]: arr
    Out[3]:
    <xarray.DataArray (x: 3, y: 5)>
    array([[ 0. ,  0.5,  1. ,  1.5,  2. ],
           [ 2.5,  3. ,  3.5,  4. ,  4.5],
           [ 5. ,  5.5,  6. ,  6.5,  7. ]])
    Coordinates:
      * x        (x) int64 0 1 2
      * y        (y) int64 0 1 2 3 4

    In [4]: arr.rolling(y=3, min_periods=2).mean()
    Out[4]:
    <xarray.DataArray (x: 3, y: 5)>
    array([[  nan,  0.25,  0.5 ,  1.  ,  1.5 ],
           [  nan,  2.75,  3.  ,  3.5 ,  4.  ],
           [  nan,  5.25,  5.5 ,  6.  ,  6.5 ]])
    Coordinates:
      * x        (x) int64 0 1 2
      * y        (y) int64 0 1 2 3 4

  See :ref:`comput.rolling` for more details. By
  `Joe Hamman <https://github.com/jhamman>`_.

Bug fixes
~~~~~~~~~

- Fixed an issue where plots using pcolormesh and Cartopy axes were being distorted
  by the inference of the axis interval breaks. This change chooses not to modify
  the coordinate variables when the axes have the attribute ``projection``, allowing
  Cartopy to handle the extent of pcolormesh plots (:issue:`781`). By
  `Joe Hamman <https://github.com/jhamman>`_.

- 2D plots now better handle additional coordinates which are not ``DataArray``
  dimensions (:issue:`788`). By `Fabien Maussion <https://github.com/fmaussion>`_.


.. _whats-new.0.7.1:

v0.7.1 (16 February 2016)
-------------------------

This is a bug fix release that includes two small, backwards compatible enhancements.
We recommend that all users upgrade.

Enhancements
~~~~~~~~~~~~

- Numerical operations now return empty objects on no overlapping labels rather
  than raising ``ValueError`` (:issue:`739`).
- :py:class:`~pandas.Series` is now supported as valid input to the ``Dataset``
  constructor (:issue:`740`).

Bug fixes
~~~~~~~~~

- Restore checks for shape consistency between data and coordinates in the
  DataArray constructor (:issue:`758`).
- Single dimension variables no longer transpose as part of a broader
  ``.transpose``. This behavior was causing ``pandas.PeriodIndex`` dimensions
  to lose their type (:issue:`749`)
- :py:class:`~xarray.Dataset` labels remain as their native type on ``.to_dataset``.
  Previously they were coerced to strings (:issue:`745`)
- Fixed a bug where replacing a ``DataArray`` index coordinate would improperly
  align the coordinate (:issue:`725`).
- ``DataArray.reindex_like`` now maintains the dtype of complex numbers when
  reindexing leads to NaN values (:issue:`738`).
- ``Dataset.rename`` and ``DataArray.rename`` support the old and new names
  being the same (:issue:`724`).
- Fix :py:meth:`~xarray.Dataset.from_dataframe` for DataFrames with Categorical
  column and a MultiIndex index (:issue:`737`).
- Fixes to ensure xarray works properly after the upcoming pandas v0.18 and
  NumPy v1.11 releases.

Acknowledgments
~~~~~~~~~~~~~~~

The following individuals contributed to this release:

- Edward Richards
- Maximilian Roos
- Rafael Guedes
- Spencer Hill
- Stephan Hoyer

.. _whats-new.0.7.0:

v0.7.0 (21 January 2016)
------------------------

This major release includes redesign of :py:class:`~xarray.DataArray`
internals, as well as new methods for reshaping, rolling and shifting
data. It includes preliminary support for :py:class:`pandas.MultiIndex`,
as well as a number of other features and bug fixes, several of which
offer improved compatibility with pandas.

New name
~~~~~~~~

The project formerly known as "xray" is now "xarray", pronounced "x-array"!
This avoids a namespace conflict with the entire field of x-ray science. Renaming
our project seemed like the right thing to do, especially because some
scientists who work with actual x-rays are interested in using this project in
their work. Thanks for your understanding and patience in this transition. You
can now find our documentation and code repository at new URLs:

- http://xarray.pydata.org
- http://github.com/pydata/xarray/

To ease the transition, we have simultaneously released v0.7.0 of both
``xray`` and ``xarray`` on the Python Package Index. These packages are
identical. For now, ``import xray`` still works, except it issues a
deprecation warning. This will be the last xray release. Going forward, we
recommend switching your import statements to ``import xarray as xr``.

.. _v0.7.0.breaking:

Breaking changes
~~~~~~~~~~~~~~~~

- The internal data model used by ``xray.DataArray`` has been
  rewritten to fix several outstanding issues (:issue:`367`, :issue:`634`,
  `this stackoverflow report`_). Internally, ``DataArray`` is now implemented
  in terms of ``._variable`` and ``._coords`` attributes instead of holding
  variables in a ``Dataset`` object.

  This refactor ensures that if a DataArray has the
  same name as one of its coordinates, the array and the coordinate no longer
  share the same data.

  In practice, this means that creating a DataArray with the same ``name`` as
  one of its dimensions no longer automatically uses that array to label the
  corresponding coordinate. You will now need to provide coordinate labels
  explicitly. Here's the old behavior:

  .. ipython::
    :verbatim:

    In [2]: xray.DataArray([4, 5, 6], dims="x", name="x")
    Out[2]:
    <xray.DataArray 'x' (x: 3)>
    array([4, 5, 6])
    Coordinates:
      * x        (x) int64 4 5 6

  and the new behavior (compare the values of the ``x`` coordinate):

  .. ipython::
    :verbatim:

    In [2]: xray.DataArray([4, 5, 6], dims="x", name="x")
    Out[2]:
    <xray.DataArray 'x' (x: 3)>
    array([4, 5, 6])
    Coordinates:
      * x        (x) int64 0 1 2

- It is no longer possible to convert a DataArray to a Dataset with
  ``xray.DataArray.to_dataset`` if it is unnamed. This will now
  raise ``ValueError``. If the array is unnamed, you need to supply the
  ``name`` argument.

.. _this stackoverflow report: http://stackoverflow.com/questions/33158558/python-xray-extract-first-and-last-time-value-within-each-month-of-a-timeseries

Enhancements
~~~~~~~~~~~~

- Basic support for :py:class:`~pandas.MultiIndex` coordinates on xray objects, including
  indexing, :py:meth:`~DataArray.stack` and :py:meth:`~DataArray.unstack`:

  .. ipython::
    :verbatim:

    In [7]: df = pd.DataFrame({"foo": range(3), "x": ["a", "b", "b"], "y": [0, 0, 1]})

    In [8]: s = df.set_index(["x", "y"])["foo"]

    In [12]: arr = xray.DataArray(s, dims="z")

    In [13]: arr
    Out[13]:
    <xray.DataArray 'foo' (z: 3)>
    array([0, 1, 2])
    Coordinates:
      * z        (z) object ('a', 0) ('b', 0) ('b', 1)

    In [19]: arr.indexes["z"]
    Out[19]:
    MultiIndex(levels=[[u'a', u'b'], [0, 1]],
               labels=[[0, 1, 1], [0, 0, 1]],
               names=[u'x', u'y'])

    In [14]: arr.unstack("z")
    Out[14]:
    <xray.DataArray 'foo' (x: 2, y: 2)>
    array([[  0.,  nan],
           [  1.,   2.]])
    Coordinates:
      * x        (x) object 'a' 'b'
      * y        (y) int64 0 1

    In [26]: arr.unstack("z").stack(z=("x", "y"))
    Out[26]:
    <xray.DataArray 'foo' (z: 4)>
    array([  0.,  nan,   1.,   2.])
    Coordinates:
      * z        (z) object ('a', 0) ('a', 1) ('b', 0) ('b', 1)

  See :ref:`reshape.stack` for more details.

  .. warning::

      xray's MultiIndex support is still experimental, and we have a long to-
      do list of desired additions (:issue:`719`), including better display of
      multi-index levels when printing a ``Dataset``, and support for saving
      datasets with a MultiIndex to a netCDF file. User contributions in this
      area would be greatly appreciated.

- Support for reading GRIB, HDF4 and other file formats via PyNIO_. See
  :ref:`io.pynio` for more details.
- Better error message when a variable is supplied with the same name as
  one of its dimensions.
- Plotting: more control on colormap parameters (:issue:`642`). ``vmin`` and
  ``vmax`` will not be silently ignored anymore. Setting ``center=False``
  prevents automatic selection of a divergent colormap.
- New ``xray.Dataset.shift`` and ``xray.Dataset.roll`` methods
  for shifting/rotating datasets or arrays along a dimension:

  .. ipython:: python
      :okwarning:

      array = xray.DataArray([5, 6, 7, 8], dims="x")
      array.shift(x=2)
      array.roll(x=2)

  Notice that ``shift`` moves data independently of coordinates, but ``roll``
  moves both data and coordinates.
- Assigning a ``pandas`` object directly as a ``Dataset`` variable is now permitted. Its
  index names correspond to the ``dims`` of the ``Dataset``, and its data is aligned.
- Passing a :py:class:`pandas.DataFrame` or ``pandas.Panel`` to a Dataset constructor
  is now permitted.
- New function ``xray.broadcast`` for explicitly broadcasting
  ``DataArray`` and ``Dataset`` objects against each other. For example:

  .. ipython:: python

      a = xray.DataArray([1, 2, 3], dims="x")
      b = xray.DataArray([5, 6], dims="y")
      a
      b
      a2, b2 = xray.broadcast(a, b)
      a2
      b2

.. _PyNIO: https://www.pyngl.ucar.edu/Nio.shtml

Bug fixes
~~~~~~~~~

- Fixes for several issues found on ``DataArray`` objects with the same name
  as one of their coordinates (see :ref:`v0.7.0.breaking` for more details).
- ``DataArray.to_masked_array`` always returns masked array with mask being an
  array (not a scalar value) (:issue:`684`)
- Allows for (imperfect) repr of Coords when underlying index is PeriodIndex (:issue:`645`).
- Fixes for several issues found on ``DataArray`` objects with the same name
  as one of their coordinates (see :ref:`v0.7.0.breaking` for more details).
- Attempting to assign a ``Dataset`` or ``DataArray`` variable/attribute using
  attribute-style syntax (e.g., ``ds.foo = 42``) now raises an error rather
  than silently failing (:issue:`656`, :issue:`714`).
- You can now pass pandas objects with non-numpy dtypes (e.g., ``categorical``
  or ``datetime64`` with a timezone) into xray without an error
  (:issue:`716`).

Acknowledgments
~~~~~~~~~~~~~~~

The following individuals contributed to this release:

- Antony Lee
- Fabien Maussion
- Joe Hamman
- Maximilian Roos
- Stephan Hoyer
- Takeshi Kanmae
- femtotrader

v0.6.1 (21 October 2015)
------------------------

This release contains a number of bug and compatibility fixes, as well
as enhancements to plotting, indexing and writing files to disk.

Note that the minimum required version of dask for use with xray is now
version 0.6.

API Changes
~~~~~~~~~~~

- The handling of colormaps and discrete color lists for 2D plots in
  ``xray.DataArray.plot`` was changed to provide more compatibility
  with matplotlib's ``contour`` and ``contourf`` functions (:issue:`538`).
  Now discrete lists of colors should be specified using ``colors`` keyword,
  rather than ``cmap``.

Enhancements
~~~~~~~~~~~~

- Faceted plotting through ``xray.plot.FacetGrid`` and the
  ``xray.plot.plot`` method. See :ref:`plotting.faceting` for more details
  and examples.
- ``xray.Dataset.sel`` and ``xray.Dataset.reindex`` now support
  the ``tolerance`` argument for controlling nearest-neighbor selection
  (:issue:`629`):

  .. ipython::
    :verbatim:

    In [5]: array = xray.DataArray([1, 2, 3], dims="x")

    In [6]: array.reindex(x=[0.9, 1.5], method="nearest", tolerance=0.2)
    Out[6]:
    <xray.DataArray (x: 2)>
    array([  2.,  nan])
    Coordinates:
      * x        (x) float64 0.9 1.5

  This feature requires pandas v0.17 or newer.
- New ``encoding`` argument in ``xray.Dataset.to_netcdf`` for writing
  netCDF files with compression, as described in the new documentation
  section on :ref:`io.netcdf.writing_encoded`.
- Add ``xray.Dataset.real`` and ``xray.Dataset.imag``
  attributes to Dataset and DataArray (:issue:`553`).
- More informative error message with ``xray.Dataset.from_dataframe``
  if the frame has duplicate columns.
- xray now uses deterministic names for dask arrays it creates or opens from
  disk. This allows xray users to take advantage of dask's nascent support for
  caching intermediate computation results. See :issue:`555` for an example.

Bug fixes
~~~~~~~~~

- Forwards compatibility with the latest pandas release (v0.17.0). We were
  using some internal pandas routines for datetime conversion, which
  unfortunately have now changed upstream (:issue:`569`).
- Aggregation functions now correctly skip ``NaN`` for data for ``complex128``
  dtype (:issue:`554`).
- Fixed indexing 0d arrays with unicode dtype (:issue:`568`).
- ``xray.DataArray.name`` and Dataset keys must be a string or None to
  be written to netCDF (:issue:`533`).
- ``xray.DataArray.where`` now uses dask instead of numpy if either the
  array or ``other`` is a dask array. Previously, if ``other`` was a numpy array
  the method was evaluated eagerly.
- Global attributes are now handled more consistently when loading remote
  datasets using ``engine='pydap'`` (:issue:`574`).
- It is now possible to assign to the ``.data`` attribute of DataArray objects.
- ``coordinates`` attribute is now kept in the encoding dictionary after
  decoding (:issue:`610`).
- Compatibility with numpy 1.10 (:issue:`617`).

Acknowledgments
~~~~~~~~~~~~~~~

The following individuals contributed to this release:

- Ryan Abernathey
- Pete Cable
- Clark Fitzgerald
- Joe Hamman
- Stephan Hoyer
- Scott Sinclair

v0.6.0 (21 August 2015)
-----------------------

This release includes numerous bug fixes and enhancements. Highlights
include the introduction of a plotting module and the new Dataset and DataArray
methods ``xray.Dataset.isel_points``, ``xray.Dataset.sel_points``,
``xray.Dataset.where`` and ``xray.Dataset.diff``. There are no
breaking changes from v0.5.2.

Enhancements
~~~~~~~~~~~~

- Plotting methods have been implemented on DataArray objects
  ``xray.DataArray.plot`` through integration with matplotlib
  (:issue:`185`). For an introduction, see :ref:`plotting`.
- Variables in netCDF files with multiple missing values are now decoded as NaN
  after issuing a warning if open_dataset is called with mask_and_scale=True.
- We clarified our rules for when the result from an xray operation is a copy
  vs. a view (see :ref:`copies_vs_views` for more details).
- Dataset variables are now written to netCDF files in order of appearance
  when using the netcdf4 backend (:issue:`479`).

- Added ``xray.Dataset.isel_points`` and ``xray.Dataset.sel_points``
  to support pointwise indexing of Datasets and DataArrays (:issue:`475`).

  .. ipython::
    :verbatim:

    In [1]: da = xray.DataArray(
       ...:     np.arange(56).reshape((7, 8)),
       ...:     coords={"x": list("abcdefg"), "y": 10 * np.arange(8)},
       ...:     dims=["x", "y"],
       ...: )

    In [2]: da
    Out[2]:
    <xray.DataArray (x: 7, y: 8)>
    array([[ 0,  1,  2,  3,  4,  5,  6,  7],
           [ 8,  9, 10, 11, 12, 13, 14, 15],
           [16, 17, 18, 19, 20, 21, 22, 23],
           [24, 25, 26, 27, 28, 29, 30, 31],
           [32, 33, 34, 35, 36, 37, 38, 39],
           [40, 41, 42, 43, 44, 45, 46, 47],
           [48, 49, 50, 51, 52, 53, 54, 55]])
    Coordinates:
    * y        (y) int64 0 10 20 30 40 50 60 70
    * x        (x) |S1 'a' 'b' 'c' 'd' 'e' 'f' 'g'

    # we can index by position along each dimension
    In [3]: da.isel_points(x=[0, 1, 6], y=[0, 1, 0], dim="points")
    Out[3]:
    <xray.DataArray (points: 3)>
    array([ 0,  9, 48])
    Coordinates:
        y        (points) int64 0 10 0
        x        (points) |S1 'a' 'b' 'g'
      * points   (points) int64 0 1 2

    # or equivalently by label
    In [9]: da.sel_points(x=["a", "b", "g"], y=[0, 10, 0], dim="points")
    Out[9]:
    <xray.DataArray (points: 3)>
    array([ 0,  9, 48])
    Coordinates:
        y        (points) int64 0 10 0
        x        (points) |S1 'a' 'b' 'g'
      * points   (points) int64 0 1 2

- New ``xray.Dataset.where`` method for masking xray objects according
  to some criteria. This works particularly well with multi-dimensional data:

  .. ipython:: python

      ds = xray.Dataset(coords={"x": range(100), "y": range(100)})
      ds["distance"] = np.sqrt(ds.x ** 2 + ds.y ** 2)

      @savefig where_example.png width=4in height=4in
      ds.distance.where(ds.distance < 100).plot()

- Added new methods ``xray.DataArray.diff`` and ``xray.Dataset.diff``
  for finite difference calculations along a given axis.

- New ``xray.DataArray.to_masked_array`` convenience method for
  returning a numpy.ma.MaskedArray.

  .. ipython:: python

      da = xray.DataArray(np.random.random_sample(size=(5, 4)))
      da.where(da < 0.5)
      da.where(da < 0.5).to_masked_array(copy=True)

- Added new flag "drop_variables" to ``xray.open_dataset`` for
  excluding variables from being parsed. This may be useful to drop
  variables with problems or inconsistent values.

Bug fixes
~~~~~~~~~

- Fixed aggregation functions (e.g., sum and mean) on big-endian arrays when
  bottleneck is installed (:issue:`489`).
- Dataset aggregation functions dropped variables with unsigned integer dtype
  (:issue:`505`).
- ``.any()`` and ``.all()`` were not lazy when used on xray objects containing
  dask arrays.
- Fixed an error when attempting to saving datetime64 variables to netCDF
  files when the first element is ``NaT`` (:issue:`528`).
- Fix pickle on DataArray objects (:issue:`515`).
- Fixed unnecessary coercion of float64 to float32 when using netcdf3 and
  netcdf4_classic formats (:issue:`526`).

v0.5.2 (16 July 2015)
---------------------

This release contains bug fixes, several additional options for opening and
saving netCDF files, and a backwards incompatible rewrite of the advanced
options for ``xray.concat``.

Backwards incompatible changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- The optional arguments ``concat_over`` and ``mode`` in ``xray.concat`` have
  been removed and replaced by ``data_vars`` and ``coords``. The new arguments are both
  more easily understood and more robustly implemented, and allowed us to fix a bug
  where ``concat`` accidentally loaded data into memory. If you set values for
  these optional arguments manually, you will need to update your code. The default
  behavior should be unchanged.

Enhancements
~~~~~~~~~~~~

- ``xray.open_mfdataset`` now supports a ``preprocess`` argument for
  preprocessing datasets prior to concatenaton. This is useful if datasets
  cannot be otherwise merged automatically, e.g., if the original datasets
  have conflicting index coordinates (:issue:`443`).
- ``xray.open_dataset`` and ``xray.open_mfdataset`` now use a
  global thread lock by default for reading from netCDF files with dask. This
  avoids possible segmentation faults for reading from netCDF4 files when HDF5
  is not configured properly for concurrent access (:issue:`444`).
- Added support for serializing arrays of complex numbers with `engine='h5netcdf'`.
- The new ``xray.save_mfdataset`` function allows for saving multiple
  datasets to disk simultaneously. This is useful when processing large datasets
  with dask.array. For example, to save a dataset too big to fit into memory
  to one file per year, we could write:

  .. ipython::
    :verbatim:

    In [1]: years, datasets = zip(*ds.groupby("time.year"))

    In [2]: paths = ["%s.nc" % y for y in years]

    In [3]: xray.save_mfdataset(datasets, paths)

Bug fixes
~~~~~~~~~

- Fixed ``min``, ``max``, ``argmin`` and ``argmax`` for arrays with string or
  unicode types (:issue:`453`).
- ``xray.open_dataset`` and ``xray.open_mfdataset`` support
  supplying chunks as a single integer.
- Fixed a bug in serializing scalar datetime variable to netCDF.
- Fixed a bug that could occur in serialization of 0-dimensional integer arrays.
- Fixed a bug where concatenating DataArrays was not always lazy (:issue:`464`).
- When reading datasets with h5netcdf, bytes attributes are decoded to strings.
  This allows conventions decoding to work properly on Python 3 (:issue:`451`).

v0.5.1 (15 June 2015)
---------------------

This minor release fixes a few bugs and an inconsistency with pandas. It also
adds the ``pipe`` method, copied from pandas.

Enhancements
~~~~~~~~~~~~

- Added ``xray.Dataset.pipe``, replicating the `new pandas method`_ in version
  0.16.2. See :ref:`transforming datasets` for more details.
- ``xray.Dataset.assign`` and ``xray.Dataset.assign_coords``
  now assign new variables in sorted (alphabetical) order, mirroring the
  behavior in pandas. Previously, the order was arbitrary.

.. _new pandas method: http://pandas.pydata.org/pandas-docs/version/0.16.2/whatsnew.html#pipe

Bug fixes
~~~~~~~~~

- ``xray.concat`` fails in an edge case involving identical coordinate variables (:issue:`425`)
- We now decode variables loaded from netCDF3 files with the scipy engine using native
  endianness (:issue:`416`). This resolves an issue when aggregating these arrays with
  bottleneck installed.

v0.5 (1 June 2015)
------------------

Highlights
~~~~~~~~~~

The headline feature in this release is experimental support for out-of-core
computing (data that doesn't fit into memory) with :doc:`user-guide/dask`. This includes a new
top-level function ``xray.open_mfdataset`` that makes it easy to open
a collection of netCDF (using dask) as a single ``xray.Dataset`` object. For
more on dask, read the `blog post introducing xray + dask`_ and the new
documentation section :doc:`user-guide/dask`.

.. _blog post introducing xray + dask: https://www.anaconda.com/blog/developer-blog/xray-dask-out-core-labeled-arrays-python/

Dask makes it possible to harness parallelism and manipulate gigantic datasets
with xray. It is currently an optional dependency, but it may become required
in the future.

Backwards incompatible changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- The logic used for choosing which variables are concatenated with
  ``xray.concat`` has changed. Previously, by default any variables
  which were equal across a dimension were not concatenated. This lead to some
  surprising behavior, where the behavior of groupby and concat operations
  could depend on runtime values (:issue:`268`). For example:

  .. ipython::
    :verbatim:

    In [1]: ds = xray.Dataset({"x": 0})

    In [2]: xray.concat([ds, ds], dim="y")
    Out[2]:
    <xray.Dataset>
    Dimensions:  ()
    Coordinates:
        *empty*
    Data variables:
        x        int64 0

  Now, the default always concatenates data variables:

  .. ipython:: python
      :suppress:

      ds = xray.Dataset({"x": 0})

  .. ipython:: python

      xray.concat([ds, ds], dim="y")

  To obtain the old behavior, supply the argument ``concat_over=[]``.

Enhancements
~~~~~~~~~~~~

- New ``xray.Dataset.to_array`` and enhanced
  ``xray.DataArray.to_dataset`` methods make it easy to switch back
  and forth between arrays and datasets:

  .. ipython:: python

      ds = xray.Dataset(
          {"a": 1, "b": ("x", [1, 2, 3])},
          coords={"c": 42},
          attrs={"Conventions": "None"},
      )
      ds.to_array()
      ds.to_array().to_dataset(dim="variable")

- New ``xray.Dataset.fillna`` method to fill missing values, modeled
  off the pandas method of the same name:

  .. ipython:: python

      array = xray.DataArray([np.nan, 1, np.nan, 3], dims="x")
      array.fillna(0)

  ``fillna`` works on both ``Dataset`` and ``DataArray`` objects, and uses
  index based alignment and broadcasting like standard binary operations. It
  also can be applied by group, as illustrated in
  :ref:`/examples/weather-data.ipynb#Fill-missing-values-with-climatology`.
- New ``xray.Dataset.assign`` and ``xray.Dataset.assign_coords``
  methods patterned off the new :py:meth:`DataFrame.assign <pandas.DataFrame.assign>`
  method in pandas:

  .. ipython:: python

      ds = xray.Dataset({"y": ("x", [1, 2, 3])})
      ds.assign(z=lambda ds: ds.y ** 2)
      ds.assign_coords(z=("x", ["a", "b", "c"]))

  These methods return a new Dataset (or DataArray) with updated data or
  coordinate variables.
- ``xray.Dataset.sel`` now supports the ``method`` parameter, which works
  like the parameter of the same name on ``xray.Dataset.reindex``. It
  provides a simple interface for doing nearest-neighbor interpolation:

  .. use verbatim because I can't seem to install pandas 0.16.1 on RTD :(

  .. ipython::
      :verbatim:

      In [12]: ds.sel(x=1.1, method="nearest")
      Out[12]:
      <xray.Dataset>
      Dimensions:  ()
      Coordinates:
          x        int64 1
      Data variables:
          y        int64 2

      In [13]: ds.sel(x=[1.1, 2.1], method="pad")
      Out[13]:
      <xray.Dataset>
      Dimensions:  (x: 2)
      Coordinates:
        * x        (x) int64 1 2
      Data variables:
          y        (x) int64 2 3

  See :ref:`nearest neighbor lookups` for more details.
- You can now control the underlying backend used for accessing remote
  datasets (via OPeNDAP) by specifying ``engine='netcdf4'`` or
  ``engine='pydap'``.
- xray now provides experimental support for reading and writing netCDF4 files directly
  via `h5py`_ with the `h5netcdf`_ package, avoiding the netCDF4-Python package. You
  will need to install h5netcdf and specify ``engine='h5netcdf'`` to try this
  feature.
- Accessing data from remote datasets now has retrying logic (with exponential
  backoff) that should make it robust to occasional bad responses from DAP
  servers.
- You can control the width of the Dataset repr with ``xray.set_options``.
  It can be used either as a context manager, in which case the default is restored
  outside the context:

  .. ipython:: python

      ds = xray.Dataset({"x": np.arange(1000)})
      with xray.set_options(display_width=40):
          print(ds)

  Or to set a global option:

  .. ipython::
      :verbatim:

      In [1]: xray.set_options(display_width=80)

  The default value for the ``display_width`` option is 80.

.. _h5py: http://www.h5py.org/
.. _h5netcdf: https://github.com/shoyer/h5netcdf

Deprecations
~~~~~~~~~~~~

- The method ``load_data()`` has been renamed to the more succinct
  ``xray.Dataset.load``.

v0.4.1 (18 March 2015)
----------------------

The release contains bug fixes and several new features. All changes should be
fully backwards compatible.

Enhancements
~~~~~~~~~~~~

- New documentation sections on :ref:`time-series` and
  :ref:`combining multiple files`.
- ``xray.Dataset.resample`` lets you resample a dataset or data array to
  a new temporal resolution. The syntax is the `same as pandas`_, except you
  need to supply the time dimension explicitly:

  .. ipython:: python
      :verbatim:

      time = pd.date_range("2000-01-01", freq="6H", periods=10)
      array = xray.DataArray(np.arange(10), [("time", time)])
      array.resample("1D", dim="time")

  You can specify how to do the resampling with the ``how`` argument and other
  options such as ``closed`` and ``label`` let you control labeling:

  .. ipython:: python
      :verbatim:

      array.resample("1D", dim="time", how="sum", label="right")

  If the desired temporal resolution is higher than the original data
  (upsampling), xray will insert missing values:

  .. ipython:: python
      :verbatim:

      array.resample("3H", "time")

- ``first`` and ``last`` methods on groupby objects let you take the first or
  last examples from each group along the grouped axis:

  .. ipython:: python
      :verbatim:

      array.groupby("time.day").first()

  These methods combine well with ``resample``:

  .. ipython:: python
      :verbatim:

      array.resample("1D", dim="time", how="first")


- ``xray.Dataset.swap_dims`` allows for easily swapping one dimension
  out for another:

  .. ipython:: python

      ds = xray.Dataset({"x": range(3), "y": ("x", list("abc"))})
      ds
      ds.swap_dims({"x": "y"})

  This was possible in earlier versions of xray, but required some contortions.
- ``xray.open_dataset`` and ``xray.Dataset.to_netcdf`` now
  accept an ``engine`` argument to explicitly select which underlying library
  (netcdf4 or scipy) is used for reading/writing a netCDF file.

.. _same as pandas: http://pandas.pydata.org/pandas-docs/stable/timeseries.html#up-and-downsampling

Bug fixes
~~~~~~~~~

- Fixed a bug where data netCDF variables read from disk with
  ``engine='scipy'`` could still be associated with the file on disk, even
  after closing the file (:issue:`341`). This manifested itself in warnings
  about mmapped arrays and segmentation faults (if the data was accessed).
- Silenced spurious warnings about all-NaN slices when using nan-aware
  aggregation methods (:issue:`344`).
- Dataset aggregations with ``keep_attrs=True`` now preserve attributes on
  data variables, not just the dataset itself.
- Tests for xray now pass when run on Windows (:issue:`360`).
- Fixed a regression in v0.4 where saving to netCDF could fail with the error
  ``ValueError: could not automatically determine time units``.

v0.4 (2 March, 2015)
--------------------

This is one of the biggest releases yet for xray: it includes some major
changes that may break existing code, along with the usual collection of minor
enhancements and bug fixes. On the plus side, this release includes all
hitherto planned breaking changes, so the upgrade path for xray should be
smoother going forward.

Breaking changes
~~~~~~~~~~~~~~~~

- We now automatically align index labels in arithmetic, dataset construction,
  merging and updating. This means the need for manually invoking methods like
  ``xray.align`` and ``xray.Dataset.reindex_like`` should be
  vastly reduced.

  :ref:`For arithmetic<math automatic alignment>`, we align
  based on the **intersection** of labels:

  .. ipython:: python

      lhs = xray.DataArray([1, 2, 3], [("x", [0, 1, 2])])
      rhs = xray.DataArray([2, 3, 4], [("x", [1, 2, 3])])
      lhs + rhs

  :ref:`For dataset construction and merging<merge>`, we align based on the
  **union** of labels:

  .. ipython:: python

      xray.Dataset({"foo": lhs, "bar": rhs})

  :ref:`For update and __setitem__<update>`, we align based on the **original**
  object:

  .. ipython:: python

      lhs.coords["rhs"] = rhs
      lhs

- Aggregations like ``mean`` or ``median`` now skip missing values by default:

  .. ipython:: python

      xray.DataArray([1, 2, np.nan, 3]).mean()

  You can turn this behavior off by supplying the keyword argument
  ``skipna=False``.

  These operations are lightning fast thanks to integration with bottleneck_,
  which is a new optional dependency for xray (numpy is used if bottleneck is
  not installed).
- Scalar coordinates no longer conflict with constant arrays with the same
  value (e.g., in arithmetic, merging datasets and concat), even if they have
  different shape (:issue:`243`). For example, the coordinate ``c`` here
  persists through arithmetic, even though it has different shapes on each
  DataArray:

  .. ipython:: python

      a = xray.DataArray([1, 2], coords={"c": 0}, dims="x")
      b = xray.DataArray([1, 2], coords={"c": ("x", [0, 0])}, dims="x")
      (a + b).coords

  This functionality can be controlled through the ``compat`` option, which
  has also been added to the ``xray.Dataset`` constructor.
- Datetime shortcuts such as ``'time.month'`` now return a ``DataArray`` with
  the name ``'month'``, not ``'time.month'`` (:issue:`345`). This makes it
  easier to index the resulting arrays when they are used with ``groupby``:

  .. ipython:: python

      time = xray.DataArray(
          pd.date_range("2000-01-01", periods=365), dims="time", name="time"
      )
      counts = time.groupby("time.month").count()
      counts.sel(month=2)

  Previously, you would need to use something like
  ``counts.sel(**{'time.month': 2}})``, which is much more awkward.
- The ``season`` datetime shortcut now returns an array of string labels
  such `'DJF'`:

  .. ipython:: python

      ds = xray.Dataset({"t": pd.date_range("2000-01-01", periods=12, freq="M")})
      ds["t.season"]

  Previously, it returned numbered seasons 1 through 4.
- We have updated our use of the terms of "coordinates" and "variables". What
  were known in previous versions of xray as "coordinates" and "variables" are
  now referred to throughout the documentation as "coordinate variables" and
  "data variables". This brings xray in closer alignment to `CF Conventions`_.
  The only visible change besides the documentation is that ``Dataset.vars``
  has been renamed ``Dataset.data_vars``.
- You will need to update your code if you have been ignoring deprecation
  warnings: methods and attributes that were deprecated in xray v0.3 or earlier
  (e.g., ``dimensions``, ``attributes```) have gone away.

.. _bottleneck: https://github.com/pydata/bottleneck

Enhancements
~~~~~~~~~~~~

- Support for ``xray.Dataset.reindex`` with a fill method. This
  provides a useful shortcut for upsampling:

  .. ipython:: python

      data = xray.DataArray([1, 2, 3], [("x", range(3))])
      data.reindex(x=[0.5, 1, 1.5, 2, 2.5], method="pad")

  This will be especially useful once pandas 0.16 is released, at which point
  xray will immediately support reindexing with
  `method='nearest' <https://github.com/pydata/pandas/pull/9258>`_.
- Use functions that return generic ndarrays with DataArray.groupby.apply and
  Dataset.apply (:issue:`327` and :issue:`329`). Thanks Jeff Gerard!
- Consolidated the functionality of ``dumps`` (writing a dataset to a netCDF3
  bytestring) into ``xray.Dataset.to_netcdf`` (:issue:`333`).
- ``xray.Dataset.to_netcdf`` now supports writing to groups in netCDF4
  files (:issue:`333`). It also finally has a full docstring -- you should read
  it!
- ``xray.open_dataset`` and ``xray.Dataset.to_netcdf`` now
  work on netCDF3 files when netcdf4-python is not installed as long as scipy
  is available (:issue:`333`).
- The new ``xray.Dataset.drop`` and ``xray.DataArray.drop`` methods
  makes it easy to drop explicitly listed variables or index labels:

  .. ipython:: python
      :okwarning:

      # drop variables
      ds = xray.Dataset({"x": 0, "y": 1})
      ds.drop("x")

      # drop index labels
      arr = xray.DataArray([1, 2, 3], coords=[("x", list("abc"))])
      arr.drop(["a", "c"], dim="x")

- ``xray.Dataset.broadcast_equals`` has been added to correspond to
  the new ``compat`` option.
- Long attributes are now truncated at 500 characters when printing a dataset
  (:issue:`338`). This should make things more convenient for working with
  datasets interactively.
- Added a new documentation example, :ref:`/examples/monthly-means.ipynb`. Thanks Joe
  Hamman!

Bug fixes
~~~~~~~~~

- Several bug fixes related to decoding time units from netCDF files
  (:issue:`316`, :issue:`330`). Thanks Stefan Pfenninger!
- xray no longer requires ``decode_coords=False`` when reading datasets with
  unparseable coordinate attributes (:issue:`308`).
- Fixed ``DataArray.loc`` indexing with ``...`` (:issue:`318`).
- Fixed an edge case that resulting in an error when reindexing
  multi-dimensional variables (:issue:`315`).
- Slicing with negative step sizes (:issue:`312`).
- Invalid conversion of string arrays to numeric dtype (:issue:`305`).
- Fixed``repr()`` on dataset objects with non-standard dates (:issue:`347`).

Deprecations
~~~~~~~~~~~~

- ``dump`` and ``dumps`` have been deprecated in favor of
  ``xray.Dataset.to_netcdf``.
- ``drop_vars`` has been deprecated in favor of ``xray.Dataset.drop``.

Future plans
~~~~~~~~~~~~

The biggest feature I'm excited about working toward in the immediate future
is supporting out-of-core operations in xray using Dask_, a part of the Blaze_
project. For a preview of using Dask with weather data, read
`this blog post`_ by Matthew Rocklin. See :issue:`328` for more details.

.. _Dask: http://dask.pydata.org
.. _Blaze: http://blaze.pydata.org
.. _this blog post: http://matthewrocklin.com/blog/work/2015/02/13/Towards-OOC-Slicing-and-Stacking/

v0.3.2 (23 December, 2014)
--------------------------

This release focused on bug-fixes, speedups and resolving some niggling
inconsistencies.

There are a few cases where the behavior of xray differs from the previous
version. However, I expect that in almost all cases your code will continue to
run unmodified.

.. warning::

    xray now requires pandas v0.15.0 or later. This was necessary for
    supporting TimedeltaIndex without too many painful hacks.

Backwards incompatible changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Arrays of :py:class:`datetime.datetime` objects are now automatically cast to
  ``datetime64[ns]`` arrays when stored in an xray object, using machinery
  borrowed from pandas:

  .. ipython:: python

      from datetime import datetime

      xray.Dataset({"t": [datetime(2000, 1, 1)]})

- xray now has support (including serialization to netCDF) for
  :py:class:`~pandas.TimedeltaIndex`. :py:class:`datetime.timedelta` objects
  are thus accordingly cast to ``timedelta64[ns]`` objects when appropriate.
- Masked arrays are now properly coerced to use ``NaN`` as a sentinel value
  (:issue:`259`).

Enhancements
~~~~~~~~~~~~

- Due to popular demand, we have added experimental attribute style access as
  a shortcut for dataset variables, coordinates and attributes:

  .. ipython:: python

      ds = xray.Dataset({"tmin": ([], 25, {"units": "celsius"})})
      ds.tmin.units

  Tab-completion for these variables should work in editors such as IPython.
  However, setting variables or attributes in this fashion is not yet
  supported because there are some unresolved ambiguities (:issue:`300`).
- You can now use a dictionary for indexing with labeled dimensions. This
  provides a safe way to do assignment with labeled dimensions:

  .. ipython:: python

      array = xray.DataArray(np.zeros(5), dims=["x"])
      array[dict(x=slice(3))] = 1
      array

- Non-index coordinates can now be faithfully written to and restored from
  netCDF files. This is done according to CF conventions when possible by
  using the ``coordinates`` attribute on a data variable. When not possible,
  xray defines a global ``coordinates`` attribute.
- Preliminary support for converting ``xray.DataArray`` objects to and from
  CDAT_ ``cdms2`` variables.
- We sped up any operation that involves creating a new Dataset or DataArray
  (e.g., indexing, aggregation, arithmetic) by a factor of 30 to 50%. The full
  speed up requires cyordereddict_ to be installed.

.. _CDAT: http://uvcdat.llnl.gov/
.. _cyordereddict: https://github.com/shoyer/cyordereddict

Bug fixes
~~~~~~~~~

- Fix for ``to_dataframe()`` with 0d string/object coordinates (:issue:`287`)
- Fix for ``to_netcdf`` with 0d string variable (:issue:`284`)
- Fix writing datetime64 arrays to netcdf if NaT is present (:issue:`270`)
- Fix align silently upcasts data arrays when NaNs are inserted (:issue:`264`)

Future plans
~~~~~~~~~~~~

- I am contemplating switching to the terms "coordinate variables" and "data
  variables" instead of the (currently used) "coordinates" and "variables",
  following their use in `CF Conventions`_ (:issue:`293`). This would mostly
  have implications for the documentation, but I would also change the
  ``Dataset`` attribute ``vars`` to ``data``.
- I no longer certain that automatic label alignment for arithmetic would be a
  good idea for xray -- it is a feature from pandas that I have not missed
  (:issue:`186`).
- The main API breakage that I *do* anticipate in the next release is finally
  making all aggregation operations skip missing values by default
  (:issue:`130`). I'm pretty sick of writing ``ds.reduce(np.nanmean, 'time')``.
- The next version of xray (0.4) will remove deprecated features and aliases
  whose use currently raises a warning.

If you have opinions about any of these anticipated changes, I would love to
hear them -- please add a note to any of the referenced GitHub issues.

.. _CF Conventions: http://cfconventions.org/Data/cf-conventions/cf-conventions-1.6/build/cf-conventions.html

v0.3.1 (22 October, 2014)
-------------------------

This is mostly a bug-fix release to make xray compatible with the latest
release of pandas (v0.15).

We added several features to better support working with missing values and
exporting xray objects to pandas. We also reorganized the internal API for
serializing and deserializing datasets, but this change should be almost
entirely transparent to users.

Other than breaking the experimental DataStore API, there should be no
backwards incompatible changes.

New features
~~~~~~~~~~~~

- Added ``xray.Dataset.count`` and ``xray.Dataset.dropna``
  methods, copied from pandas, for working with missing values (:issue:`247`,
  :issue:`58`).
- Added ``xray.DataArray.to_pandas`` for
  converting a data array into the pandas object with the same dimensionality
  (1D to Series, 2D to DataFrame, etc.) (:issue:`255`).
- Support for reading gzipped netCDF3 files (:issue:`239`).
- Reduced memory usage when writing netCDF files (:issue:`251`).
- 'missing_value' is now supported as an alias for the '_FillValue' attribute
  on netCDF variables (:issue:`245`).
- Trivial indexes, equivalent to ``range(n)`` where ``n`` is the length of the
  dimension, are no longer written to disk (:issue:`245`).

Bug fixes
~~~~~~~~~

- Compatibility fixes for pandas v0.15 (:issue:`262`).
- Fixes for display and indexing of ``NaT`` (not-a-time) (:issue:`238`,
  :issue:`240`)
- Fix slicing by label was an argument is a data array (:issue:`250`).
- Test data is now shipped with the source distribution (:issue:`253`).
- Ensure order does not matter when doing arithmetic with scalar data arrays
  (:issue:`254`).
- Order of dimensions preserved with ``DataArray.to_dataframe`` (:issue:`260`).

v0.3 (21 September 2014)
------------------------

New features
~~~~~~~~~~~~

- **Revamped coordinates**: "coordinates" now refer to all arrays that are not
  used to index a dimension. Coordinates are intended to allow for keeping track
  of arrays of metadata that describe the grid on which the points in "variable"
  arrays lie. They are preserved (when unambiguous) even though mathematical
  operations.
- **Dataset math** ``xray.Dataset`` objects now support all arithmetic
  operations directly. Dataset-array operations map across all dataset
  variables; dataset-dataset operations act on each pair of variables with the
  same name.
- **GroupBy math**: This provides a convenient shortcut for normalizing by the
  average value of a group.
- The dataset ``__repr__`` method has been entirely overhauled; dataset
  objects now show their values when printed.
- You can now index a dataset with a list of variables to return a new dataset:
  ``ds[['foo', 'bar']]``.

Backwards incompatible changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``Dataset.__eq__`` and ``Dataset.__ne__`` are now element-wise operations
  instead of comparing all values to obtain a single boolean. Use the method
  ``xray.Dataset.equals`` instead.

Deprecations
~~~~~~~~~~~~

- ``Dataset.noncoords`` is deprecated: use ``Dataset.vars`` instead.
- ``Dataset.select_vars`` deprecated: index a ``Dataset`` with a list of
  variable names instead.
- ``DataArray.select_vars`` and ``DataArray.drop_vars`` deprecated: use
  ``xray.DataArray.reset_coords`` instead.

v0.2 (14 August 2014)
---------------------

This is major release that includes some new features and quite a few bug
fixes. Here are the highlights:

- There is now a direct constructor for ``DataArray`` objects, which makes it
  possible to create a DataArray without using a Dataset. This is highlighted
  in the refreshed ``tutorial``.
- You can perform aggregation operations like ``mean`` directly on
  ``xray.Dataset`` objects, thanks to Joe Hamman. These aggregation
  methods also worked on grouped datasets.
- xray now works on Python 2.6, thanks to Anna Kuznetsova.
- A number of methods and attributes were given more sensible (usually shorter)
  names: ``labeled`` -> ``sel``,  ``indexed`` -> ``isel``, ``select`` ->
  ``select_vars``, ``unselect`` -> ``drop_vars``, ``dimensions`` -> ``dims``,
  ``coordinates`` -> ``coords``, ``attributes`` -> ``attrs``.
- New ``xray.Dataset.load_data`` and ``xray.Dataset.close``
  methods for datasets facilitate lower level of control of data loaded from
  disk.

v0.1.1 (20 May 2014)
--------------------

xray 0.1.1 is a bug-fix release that includes changes that should be almost
entirely backwards compatible with v0.1:

- Python 3 support (:issue:`53`)
- Required numpy version relaxed to 1.7 (:issue:`129`)
- Return numpy.datetime64 arrays for non-standard calendars (:issue:`126`)
- Support for opening datasets associated with NetCDF4 groups (:issue:`127`)
- Bug-fixes for concatenating datetime arrays (:issue:`134`)

Special thanks to new contributors Thomas Kluyver, Joe Hamman and Alistair
Miles.

v0.1 (2 May 2014)
-----------------

Initial release.
.. _ecosystem:

Xarray related projects
-----------------------

Below is a list of existing open source projects that build
functionality upon xarray. See also section :ref:`internals` for more
details on how to build xarray extensions. We also maintain the
`xarray-contrib <https://github.com/xarray-contrib>`_ GitHub organization
as a place to curate projects that build upon xarray.

Geosciences
~~~~~~~~~~~

- `aospy <https://aospy.readthedocs.io>`_: Automated analysis and management of gridded climate data.
- `argopy <https://github.com/euroargodev/argopy>`_: xarray-based Argo data access, manipulation and visualisation for standard users as well as Argo experts.
- `climpred <https://climpred.readthedocs.io>`_: Analysis of ensemble forecast models for climate prediction.
- `geocube <https://corteva.github.io/geocube>`_: Tool to convert geopandas vector data into rasterized xarray data.
- `GeoWombat <https://github.com/jgrss/geowombat>`_: Utilities for analysis of remotely sensed and gridded raster data at scale (easily tame Landsat, Sentinel, Quickbird, and PlanetScope).
- `infinite-diff <https://github.com/spencerahill/infinite-diff>`_: xarray-based finite-differencing, focused on gridded climate/meteorology data
- `marc_analysis <https://github.com/darothen/marc_analysis>`_: Analysis package for CESM/MARC experiments and output.
- `MetPy <https://unidata.github.io/MetPy/dev/index.html>`_: A collection of tools in Python for reading, visualizing, and performing calculations with weather data.
- `MPAS-Analysis <http://mpas-analysis.readthedocs.io>`_: Analysis for simulations produced with Model for Prediction Across Scales (MPAS) components and the Accelerated Climate Model for Energy (ACME).
- `OGGM <http://oggm.org/>`_: Open Global Glacier Model
- `Oocgcm <https://oocgcm.readthedocs.io/>`_: Analysis of large gridded geophysical datasets
- `Open Data Cube <https://www.opendatacube.org/>`_: Analysis toolkit of continental scale Earth Observation data from satellites.
- `Pangaea: <https://pangaea.readthedocs.io/en/latest/>`_: xarray extension for gridded land surface & weather model output).
- `Pangeo <https://pangeo-data.github.io>`_: A community effort for big data geoscience in the cloud.
- `PyGDX <https://pygdx.readthedocs.io/en/latest/>`_: Python 3 package for
  accessing data stored in GAMS Data eXchange (GDX) files. Also uses a custom
  subclass.
- `pyinterp <https://pangeo-pyinterp.readthedocs.io/en/latest/>`_: Python 3 package for interpolating geo-referenced data used in the field of geosciences.
- `pyXpcm <https://pyxpcm.readthedocs.io>`_: xarray-based Profile Classification Modelling (PCM), mostly for ocean data.
- `Regionmask <https://regionmask.readthedocs.io/>`_: plotting and creation of masks of spatial regions
- `rioxarray <https://corteva.github.io/rioxarray>`_: geospatial xarray extension powered by rasterio
- `salem <https://salem.readthedocs.io>`_: Adds geolocalised subsetting, masking, and plotting operations to xarray's data structures via accessors.
- `SatPy <https://satpy.readthedocs.io/>`_ : Library for reading and manipulating meteorological remote sensing data and writing it to various image and data file formats.
- `Spyfit <https://spyfit.readthedocs.io/en/master/>`_: FTIR spectroscopy of the atmosphere
- `windspharm <https://ajdawson.github.io/windspharm/index.html>`_: Spherical
  harmonic wind analysis in Python.
- `wradlib <https://wradlib.org/>`_: An Open Source Library for Weather Radar Data Processing.
- `wrf-python <https://wrf-python.readthedocs.io/>`_: A collection of diagnostic and interpolation routines for use with output of the Weather Research and Forecasting (WRF-ARW) Model.
- `xarray-simlab <https://xarray-simlab.readthedocs.io>`_: xarray extension for computer model simulations.
- `xarray-spatial <https://makepath.github.io/xarray-spatial>`_: Numba-accelerated raster-based spatial processing tools (NDVI, curvature, zonal-statistics, proximity, hillshading, viewshed, etc.)
- `xarray-topo <https://gitext.gfz-potsdam.de/sec55-public/xarray-topo>`_: xarray extension for topographic analysis and modelling.
- `xbpch <https://github.com/darothen/xbpch>`_: xarray interface for bpch files.
- `xclim <https://xclim.readthedocs.io/>`_: A library for calculating climate science indices with unit handling built from xarray and dask.
- `xESMF <https://pangeo-xesmf.readthedocs.io/>`_: Universal regridder for geospatial data.
- `xgcm <https://xgcm.readthedocs.io/>`_: Extends the xarray data model to understand finite volume grid cells (common in General Circulation Models) and provides interpolation and difference operations for such grids.
- `xmitgcm <http://xgcm.readthedocs.io/>`_: a python package for reading `MITgcm <http://mitgcm.org/>`_ binary MDS files into xarray data structures.
- `xnemogcm <https://github.com/rcaneill/xnemogcm/>`_: a package to read `NEMO <https://nemo-ocean.eu/>`_ output files and add attributes to interface with xgcm.

Machine Learning
~~~~~~~~~~~~~~~~
- `ArviZ <https://arviz-devs.github.io/arviz/>`_: Exploratory analysis of Bayesian models, built on top of xarray.
- `Darts <https://github.com/unit8co/darts/>`_: User-friendly modern machine learning for time series in Python.
- `Elm <https://ensemble-learning-models.readthedocs.io>`_: Parallel machine learning on xarray data structures
- `sklearn-xarray (1) <https://phausamann.github.io/sklearn-xarray>`_: Combines scikit-learn and xarray (1).
- `sklearn-xarray (2) <https://sklearn-xarray.readthedocs.io/en/latest/>`_: Combines scikit-learn and xarray (2).

Other domains
~~~~~~~~~~~~~
- `ptsa <https://pennmem.github.io/ptsa/html/index.html>`_: EEG Time Series Analysis
- `pycalphad <https://pycalphad.org/docs/latest/>`_: Computational Thermodynamics in Python
- `pyomeca <https://pyomeca.github.io/>`_: Python framework for biomechanical analysis

Extend xarray capabilities
~~~~~~~~~~~~~~~~~~~~~~~~~~
- `Collocate <https://github.com/cistools/collocate>`_: Collocate xarray trajectories in arbitrary physical dimensions
- `eofs <https://ajdawson.github.io/eofs/>`_: EOF analysis in Python.
- `hypothesis-gufunc <https://hypothesis-gufunc.readthedocs.io/en/latest/>`_: Extension to hypothesis. Makes it easy to write unit tests with xarray objects as input.
- `nxarray <https://github.com/nxarray/nxarray>`_: NeXus input/output capability for xarray.
- `xarray-compare <https://github.com/astropenguin/xarray-compare>`_: xarray extension for data comparison.
- `xarray-dataclasses <https://github.com/astropenguin/xarray-dataclasses>`_: xarray extension for typed DataArray and Dataset creation.
- `xarray_extras <https://github.com/crusaderky/xarray_extras>`_: Advanced algorithms for xarray objects (e.g. integrations/interpolations).
- `xpublish <https://xpublish.readthedocs.io/>`_: Publish Xarray Datasets via a Zarr compatible REST API.
- `xrft <https://github.com/rabernat/xrft>`_: Fourier transforms for xarray data.
- `xr-scipy <https://xr-scipy.readthedocs.io>`_: A lightweight scipy wrapper for xarray.
- `X-regression <https://github.com/kuchaale/X-regression>`_: Multiple linear regression from Statsmodels library coupled with Xarray library.
- `xskillscore <https://github.com/xarray-contrib/xskillscore>`_: Metrics for verifying forecasts.
- `xyzpy <http://xyzpy.readthedocs.io>`_: Easily generate high dimensional data, including parallelization.

Visualization
~~~~~~~~~~~~~
- `datashader <https://datashader.org>`_, `geoviews <http://geoviews.org>`_, `holoviews <http://holoviews.org/>`_, : visualization packages for large data.
- `hvplot <https://hvplot.pyviz.org/>`_ : A high-level plotting API for the PyData ecosystem built on HoloViews.
- `psyplot <https://psyplot.readthedocs.io>`_: Interactive data visualization with python.
- `xarray-leaflet <https://github.com/davidbrochart/xarray_leaflet>`_: An xarray extension for tiled map plotting based on ipyleaflet.
- `xtrude <https://github.com/davidbrochart/xtrude>`_: An xarray extension for 3D terrain visualization based on pydeck.

Non-Python projects
~~~~~~~~~~~~~~~~~~~
- `xframe <https://github.com/QuantStack/xframe>`_: C++ data structures inspired by xarray.
- `AxisArrays <https://github.com/JuliaArrays/AxisArrays.jl>`_ and
  `NamedArrays <https://github.com/davidavdav/NamedArrays.jl>`_: similar data structures for Julia.

More projects can be found at the `"xarray" Github topic <https://github.com/topics/xarray>`_.
.. Generate API reference pages, but don't display these in tables.
.. This extra page is a work around for sphinx not having any support for
.. hiding an autosummary table.

:orphan:

.. currentmodule:: xarray

.. autosummary::
   :toctree: generated/

   core.coordinates.DatasetCoordinates.get
   core.coordinates.DatasetCoordinates.items
   core.coordinates.DatasetCoordinates.keys
   core.coordinates.DatasetCoordinates.merge
   core.coordinates.DatasetCoordinates.to_dataset
   core.coordinates.DatasetCoordinates.to_index
   core.coordinates.DatasetCoordinates.update
   core.coordinates.DatasetCoordinates.values
   core.coordinates.DatasetCoordinates.dims
   core.coordinates.DatasetCoordinates.indexes
   core.coordinates.DatasetCoordinates.variables

   core.rolling.DatasetCoarsen.boundary
   core.rolling.DatasetCoarsen.coord_func
   core.rolling.DatasetCoarsen.obj
   core.rolling.DatasetCoarsen.side
   core.rolling.DatasetCoarsen.trim_excess
   core.rolling.DatasetCoarsen.windows

   core.rolling.DatasetRolling.center
   core.rolling.DatasetRolling.dim
   core.rolling.DatasetRolling.min_periods
   core.rolling.DatasetRolling.obj
   core.rolling.DatasetRolling.rollings
   core.rolling.DatasetRolling.window

   core.weighted.DatasetWeighted.obj
   core.weighted.DatasetWeighted.weights

   Dataset.load_store
   Dataset.dump_to_store

   DataArray.astype
   DataArray.item

   core.coordinates.DataArrayCoordinates.get
   core.coordinates.DataArrayCoordinates.items
   core.coordinates.DataArrayCoordinates.keys
   core.coordinates.DataArrayCoordinates.merge
   core.coordinates.DataArrayCoordinates.to_dataset
   core.coordinates.DataArrayCoordinates.to_index
   core.coordinates.DataArrayCoordinates.update
   core.coordinates.DataArrayCoordinates.values
   core.coordinates.DataArrayCoordinates.dims
   core.coordinates.DataArrayCoordinates.indexes
   core.coordinates.DataArrayCoordinates.variables

   core.rolling.DataArrayCoarsen.boundary
   core.rolling.DataArrayCoarsen.coord_func
   core.rolling.DataArrayCoarsen.obj
   core.rolling.DataArrayCoarsen.side
   core.rolling.DataArrayCoarsen.trim_excess
   core.rolling.DataArrayCoarsen.windows

   core.rolling.DataArrayRolling.center
   core.rolling.DataArrayRolling.dim
   core.rolling.DataArrayRolling.min_periods
   core.rolling.DataArrayRolling.obj
   core.rolling.DataArrayRolling.window
   core.rolling.DataArrayRolling.window_labels

   core.weighted.DataArrayWeighted.obj
   core.weighted.DataArrayWeighted.weights

   core.accessor_dt.DatetimeAccessor.ceil
   core.accessor_dt.DatetimeAccessor.floor
   core.accessor_dt.DatetimeAccessor.round
   core.accessor_dt.DatetimeAccessor.strftime
   core.accessor_dt.DatetimeAccessor.calendar
   core.accessor_dt.DatetimeAccessor.date
   core.accessor_dt.DatetimeAccessor.day
   core.accessor_dt.DatetimeAccessor.dayofweek
   core.accessor_dt.DatetimeAccessor.dayofyear
   core.accessor_dt.DatetimeAccessor.days_in_month
   core.accessor_dt.DatetimeAccessor.daysinmonth
   core.accessor_dt.DatetimeAccessor.hour
   core.accessor_dt.DatetimeAccessor.is_leap_year
   core.accessor_dt.DatetimeAccessor.is_month_end
   core.accessor_dt.DatetimeAccessor.is_month_start
   core.accessor_dt.DatetimeAccessor.is_quarter_end
   core.accessor_dt.DatetimeAccessor.is_quarter_start
   core.accessor_dt.DatetimeAccessor.is_year_end
   core.accessor_dt.DatetimeAccessor.is_year_start
   core.accessor_dt.DatetimeAccessor.isocalendar
   core.accessor_dt.DatetimeAccessor.microsecond
   core.accessor_dt.DatetimeAccessor.minute
   core.accessor_dt.DatetimeAccessor.month
   core.accessor_dt.DatetimeAccessor.nanosecond
   core.accessor_dt.DatetimeAccessor.quarter
   core.accessor_dt.DatetimeAccessor.season
   core.accessor_dt.DatetimeAccessor.second
   core.accessor_dt.DatetimeAccessor.time
   core.accessor_dt.DatetimeAccessor.week
   core.accessor_dt.DatetimeAccessor.weekday
   core.accessor_dt.DatetimeAccessor.weekday_name
   core.accessor_dt.DatetimeAccessor.weekofyear
   core.accessor_dt.DatetimeAccessor.year

   core.accessor_dt.TimedeltaAccessor.ceil
   core.accessor_dt.TimedeltaAccessor.floor
   core.accessor_dt.TimedeltaAccessor.round
   core.accessor_dt.TimedeltaAccessor.days
   core.accessor_dt.TimedeltaAccessor.microseconds
   core.accessor_dt.TimedeltaAccessor.nanoseconds
   core.accessor_dt.TimedeltaAccessor.seconds

   core.accessor_str.StringAccessor.capitalize
   core.accessor_str.StringAccessor.casefold
   core.accessor_str.StringAccessor.cat
   core.accessor_str.StringAccessor.center
   core.accessor_str.StringAccessor.contains
   core.accessor_str.StringAccessor.count
   core.accessor_str.StringAccessor.decode
   core.accessor_str.StringAccessor.encode
   core.accessor_str.StringAccessor.endswith
   core.accessor_str.StringAccessor.extract
   core.accessor_str.StringAccessor.extractall
   core.accessor_str.StringAccessor.find
   core.accessor_str.StringAccessor.findall
   core.accessor_str.StringAccessor.format
   core.accessor_str.StringAccessor.get
   core.accessor_str.StringAccessor.get_dummies
   core.accessor_str.StringAccessor.index
   core.accessor_str.StringAccessor.isalnum
   core.accessor_str.StringAccessor.isalpha
   core.accessor_str.StringAccessor.isdecimal
   core.accessor_str.StringAccessor.isdigit
   core.accessor_str.StringAccessor.islower
   core.accessor_str.StringAccessor.isnumeric
   core.accessor_str.StringAccessor.isspace
   core.accessor_str.StringAccessor.istitle
   core.accessor_str.StringAccessor.isupper
   core.accessor_str.StringAccessor.join
   core.accessor_str.StringAccessor.len
   core.accessor_str.StringAccessor.ljust
   core.accessor_str.StringAccessor.lower
   core.accessor_str.StringAccessor.lstrip
   core.accessor_str.StringAccessor.match
   core.accessor_str.StringAccessor.normalize
   core.accessor_str.StringAccessor.pad
   core.accessor_str.StringAccessor.partition
   core.accessor_str.StringAccessor.repeat
   core.accessor_str.StringAccessor.replace
   core.accessor_str.StringAccessor.rfind
   core.accessor_str.StringAccessor.rindex
   core.accessor_str.StringAccessor.rjust
   core.accessor_str.StringAccessor.rpartition
   core.accessor_str.StringAccessor.rsplit
   core.accessor_str.StringAccessor.rstrip
   core.accessor_str.StringAccessor.slice
   core.accessor_str.StringAccessor.slice_replace
   core.accessor_str.StringAccessor.split
   core.accessor_str.StringAccessor.startswith
   core.accessor_str.StringAccessor.strip
   core.accessor_str.StringAccessor.swapcase
   core.accessor_str.StringAccessor.title
   core.accessor_str.StringAccessor.translate
   core.accessor_str.StringAccessor.upper
   core.accessor_str.StringAccessor.wrap
   core.accessor_str.StringAccessor.zfill

   Variable.all
   Variable.any
   Variable.argmax
   Variable.argmin
   Variable.argsort
   Variable.astype
   Variable.broadcast_equals
   Variable.chunk
   Variable.clip
   Variable.coarsen
   Variable.compute
   Variable.concat
   Variable.conj
   Variable.conjugate
   Variable.copy
   Variable.count
   Variable.cumprod
   Variable.cumsum
   Variable.equals
   Variable.fillna
   Variable.get_axis_num
   Variable.identical
   Variable.isel
   Variable.isnull
   Variable.item
   Variable.load
   Variable.max
   Variable.mean
   Variable.median
   Variable.min
   Variable.no_conflicts
   Variable.notnull
   Variable.pad
   Variable.prod
   Variable.quantile
   Variable.rank
   Variable.reduce
   Variable.roll
   Variable.rolling_window
   Variable.round
   Variable.searchsorted
   Variable.set_dims
   Variable.shift
   Variable.squeeze
   Variable.stack
   Variable.std
   Variable.sum
   Variable.to_base_variable
   Variable.to_coord
   Variable.to_dict
   Variable.to_index
   Variable.to_index_variable
   Variable.to_variable
   Variable.transpose
   Variable.unstack
   Variable.var
   Variable.where
   Variable.T
   Variable.attrs
   Variable.chunks
   Variable.data
   Variable.dims
   Variable.dtype
   Variable.encoding
   Variable.imag
   Variable.nbytes
   Variable.ndim
   Variable.real
   Variable.shape
   Variable.size
   Variable.sizes
   Variable.values

   IndexVariable.all
   IndexVariable.any
   IndexVariable.argmax
   IndexVariable.argmin
   IndexVariable.argsort
   IndexVariable.astype
   IndexVariable.broadcast_equals
   IndexVariable.chunk
   IndexVariable.clip
   IndexVariable.coarsen
   IndexVariable.compute
   IndexVariable.concat
   IndexVariable.conj
   IndexVariable.conjugate
   IndexVariable.copy
   IndexVariable.count
   IndexVariable.cumprod
   IndexVariable.cumsum
   IndexVariable.equals
   IndexVariable.fillna
   IndexVariable.get_axis_num
   IndexVariable.get_level_variable
   IndexVariable.identical
   IndexVariable.isel
   IndexVariable.isnull
   IndexVariable.item
   IndexVariable.load
   IndexVariable.max
   IndexVariable.mean
   IndexVariable.median
   IndexVariable.min
   IndexVariable.no_conflicts
   IndexVariable.notnull
   IndexVariable.pad
   IndexVariable.prod
   IndexVariable.quantile
   IndexVariable.rank
   IndexVariable.reduce
   IndexVariable.roll
   IndexVariable.rolling_window
   IndexVariable.round
   IndexVariable.searchsorted
   IndexVariable.set_dims
   IndexVariable.shift
   IndexVariable.squeeze
   IndexVariable.stack
   IndexVariable.std
   IndexVariable.sum
   IndexVariable.to_base_variable
   IndexVariable.to_coord
   IndexVariable.to_dict
   IndexVariable.to_index
   IndexVariable.to_index_variable
   IndexVariable.to_variable
   IndexVariable.transpose
   IndexVariable.unstack
   IndexVariable.var
   IndexVariable.where
   IndexVariable.T
   IndexVariable.attrs
   IndexVariable.chunks
   IndexVariable.data
   IndexVariable.dims
   IndexVariable.dtype
   IndexVariable.encoding
   IndexVariable.imag
   IndexVariable.level_names
   IndexVariable.name
   IndexVariable.nbytes
   IndexVariable.ndim
   IndexVariable.real
   IndexVariable.shape
   IndexVariable.size
   IndexVariable.sizes
   IndexVariable.values

   ufuncs.angle
   ufuncs.arccos
   ufuncs.arccosh
   ufuncs.arcsin
   ufuncs.arcsinh
   ufuncs.arctan
   ufuncs.arctan2
   ufuncs.arctanh
   ufuncs.ceil
   ufuncs.conj
   ufuncs.copysign
   ufuncs.cos
   ufuncs.cosh
   ufuncs.deg2rad
   ufuncs.degrees
   ufuncs.exp
   ufuncs.expm1
   ufuncs.fabs
   ufuncs.fix
   ufuncs.floor
   ufuncs.fmax
   ufuncs.fmin
   ufuncs.fmod
   ufuncs.fmod
   ufuncs.frexp
   ufuncs.hypot
   ufuncs.imag
   ufuncs.iscomplex
   ufuncs.isfinite
   ufuncs.isinf
   ufuncs.isnan
   ufuncs.isreal
   ufuncs.ldexp
   ufuncs.log
   ufuncs.log10
   ufuncs.log1p
   ufuncs.log2
   ufuncs.logaddexp
   ufuncs.logaddexp2
   ufuncs.logical_and
   ufuncs.logical_not
   ufuncs.logical_or
   ufuncs.logical_xor
   ufuncs.maximum
   ufuncs.minimum
   ufuncs.nextafter
   ufuncs.rad2deg
   ufuncs.radians
   ufuncs.real
   ufuncs.rint
   ufuncs.sign
   ufuncs.signbit
   ufuncs.sin
   ufuncs.sinh
   ufuncs.sqrt
   ufuncs.square
   ufuncs.tan
   ufuncs.tanh
   ufuncs.trunc

   plot.plot
   plot.line
   plot.step
   plot.hist
   plot.contour
   plot.contourf
   plot.imshow
   plot.pcolormesh
   plot.scatter
   plot.surface

   plot.FacetGrid.map_dataarray
   plot.FacetGrid.set_titles
   plot.FacetGrid.set_ticks
   plot.FacetGrid.map

   CFTimeIndex.all
   CFTimeIndex.any
   CFTimeIndex.append
   CFTimeIndex.argsort
   CFTimeIndex.argmax
   CFTimeIndex.argmin
   CFTimeIndex.asof
   CFTimeIndex.asof_locs
   CFTimeIndex.astype
   CFTimeIndex.calendar
   CFTimeIndex.ceil
   CFTimeIndex.contains
   CFTimeIndex.copy
   CFTimeIndex.days_in_month
   CFTimeIndex.delete
   CFTimeIndex.difference
   CFTimeIndex.drop
   CFTimeIndex.drop_duplicates
   CFTimeIndex.droplevel
   CFTimeIndex.dropna
   CFTimeIndex.duplicated
   CFTimeIndex.equals
   CFTimeIndex.factorize
   CFTimeIndex.fillna
   CFTimeIndex.floor
   CFTimeIndex.format
   CFTimeIndex.get_indexer
   CFTimeIndex.get_indexer_for
   CFTimeIndex.get_indexer_non_unique
   CFTimeIndex.get_level_values
   CFTimeIndex.get_loc
   CFTimeIndex.get_slice_bound
   CFTimeIndex.get_value
   CFTimeIndex.groupby
   CFTimeIndex.holds_integer
   CFTimeIndex.identical
   CFTimeIndex.insert
   CFTimeIndex.intersection
   CFTimeIndex.is_
   CFTimeIndex.is_boolean
   CFTimeIndex.is_categorical
   CFTimeIndex.is_floating
   CFTimeIndex.is_integer
   CFTimeIndex.is_interval
   CFTimeIndex.is_mixed
   CFTimeIndex.is_numeric
   CFTimeIndex.is_object
   CFTimeIndex.is_type_compatible
   CFTimeIndex.isin
   CFTimeIndex.isna
   CFTimeIndex.isnull
   CFTimeIndex.item
   CFTimeIndex.join
   CFTimeIndex.map
   CFTimeIndex.max
   CFTimeIndex.memory_usage
   CFTimeIndex.min
   CFTimeIndex.notna
   CFTimeIndex.notnull
   CFTimeIndex.nunique
   CFTimeIndex.putmask
   CFTimeIndex.ravel
   CFTimeIndex.reindex
   CFTimeIndex.rename
   CFTimeIndex.repeat
   CFTimeIndex.round
   CFTimeIndex.searchsorted
   CFTimeIndex.set_names
   CFTimeIndex.set_value
   CFTimeIndex.shift
   CFTimeIndex.slice_indexer
   CFTimeIndex.slice_locs
   CFTimeIndex.sort
   CFTimeIndex.sort_values
   CFTimeIndex.sortlevel
   CFTimeIndex.strftime
   CFTimeIndex.symmetric_difference
   CFTimeIndex.take
   CFTimeIndex.to_datetimeindex
   CFTimeIndex.to_flat_index
   CFTimeIndex.to_frame
   CFTimeIndex.to_list
   CFTimeIndex.to_native_types
   CFTimeIndex.to_numpy
   CFTimeIndex.to_series
   CFTimeIndex.tolist
   CFTimeIndex.transpose
   CFTimeIndex.union
   CFTimeIndex.unique
   CFTimeIndex.value_counts
   CFTimeIndex.view
   CFTimeIndex.where

   CFTimeIndex.T
   CFTimeIndex.array
   CFTimeIndex.asi8
   CFTimeIndex.date_type
   CFTimeIndex.day
   CFTimeIndex.dayofweek
   CFTimeIndex.dayofyear
   CFTimeIndex.dtype
   CFTimeIndex.empty
   CFTimeIndex.freq
   CFTimeIndex.has_duplicates
   CFTimeIndex.hasnans
   CFTimeIndex.hour
   CFTimeIndex.inferred_type
   CFTimeIndex.is_all_dates
   CFTimeIndex.is_monotonic
   CFTimeIndex.is_monotonic_increasing
   CFTimeIndex.is_monotonic_decreasing
   CFTimeIndex.is_unique
   CFTimeIndex.microsecond
   CFTimeIndex.minute
   CFTimeIndex.month
   CFTimeIndex.name
   CFTimeIndex.names
   CFTimeIndex.nbytes
   CFTimeIndex.ndim
   CFTimeIndex.nlevels
   CFTimeIndex.second
   CFTimeIndex.shape
   CFTimeIndex.size
   CFTimeIndex.values
   CFTimeIndex.year

   backends.NetCDF4DataStore.close
   backends.NetCDF4DataStore.encode
   backends.NetCDF4DataStore.encode_attribute
   backends.NetCDF4DataStore.encode_variable
   backends.NetCDF4DataStore.get_attrs
   backends.NetCDF4DataStore.get_dimensions
   backends.NetCDF4DataStore.get_encoding
   backends.NetCDF4DataStore.get_variables
   backends.NetCDF4DataStore.load
   backends.NetCDF4DataStore.open
   backends.NetCDF4DataStore.open_store_variable
   backends.NetCDF4DataStore.prepare_variable
   backends.NetCDF4DataStore.set_attribute
   backends.NetCDF4DataStore.set_attributes
   backends.NetCDF4DataStore.set_dimension
   backends.NetCDF4DataStore.set_dimensions
   backends.NetCDF4DataStore.set_variable
   backends.NetCDF4DataStore.set_variables
   backends.NetCDF4DataStore.store
   backends.NetCDF4DataStore.store_dataset
   backends.NetCDF4DataStore.sync
   backends.NetCDF4DataStore.autoclose
   backends.NetCDF4DataStore.ds
   backends.NetCDF4DataStore.format
   backends.NetCDF4DataStore.is_remote
   backends.NetCDF4DataStore.lock

   backends.H5NetCDFStore.autoclose
   backends.H5NetCDFStore.close
   backends.H5NetCDFStore.encode
   backends.H5NetCDFStore.encode_attribute
   backends.H5NetCDFStore.encode_variable
   backends.H5NetCDFStore.format
   backends.H5NetCDFStore.get_attrs
   backends.H5NetCDFStore.get_dimensions
   backends.H5NetCDFStore.get_encoding
   backends.H5NetCDFStore.get_variables
   backends.H5NetCDFStore.is_remote
   backends.H5NetCDFStore.load
   backends.H5NetCDFStore.lock
   backends.H5NetCDFStore.open
   backends.H5NetCDFStore.open_store_variable
   backends.H5NetCDFStore.prepare_variable
   backends.H5NetCDFStore.set_attribute
   backends.H5NetCDFStore.set_attributes
   backends.H5NetCDFStore.set_dimension
   backends.H5NetCDFStore.set_dimensions
   backends.H5NetCDFStore.set_variable
   backends.H5NetCDFStore.set_variables
   backends.H5NetCDFStore.store
   backends.H5NetCDFStore.store_dataset
   backends.H5NetCDFStore.sync
   backends.H5NetCDFStore.ds

   backends.PydapDataStore.close
   backends.PydapDataStore.get_attrs
   backends.PydapDataStore.get_dimensions
   backends.PydapDataStore.get_encoding
   backends.PydapDataStore.get_variables
   backends.PydapDataStore.load
   backends.PydapDataStore.open
   backends.PydapDataStore.open_store_variable

   backends.ScipyDataStore.close
   backends.ScipyDataStore.encode
   backends.ScipyDataStore.encode_attribute
   backends.ScipyDataStore.encode_variable
   backends.ScipyDataStore.get_attrs
   backends.ScipyDataStore.get_dimensions
   backends.ScipyDataStore.get_encoding
   backends.ScipyDataStore.get_variables
   backends.ScipyDataStore.load
   backends.ScipyDataStore.open_store_variable
   backends.ScipyDataStore.prepare_variable
   backends.ScipyDataStore.set_attribute
   backends.ScipyDataStore.set_attributes
   backends.ScipyDataStore.set_dimension
   backends.ScipyDataStore.set_dimensions
   backends.ScipyDataStore.set_variable
   backends.ScipyDataStore.set_variables
   backends.ScipyDataStore.store
   backends.ScipyDataStore.store_dataset
   backends.ScipyDataStore.sync
   backends.ScipyDataStore.ds

   backends.FileManager.acquire
   backends.FileManager.acquire_context
   backends.FileManager.close

   backends.CachingFileManager.acquire
   backends.CachingFileManager.acquire_context
   backends.CachingFileManager.close

   backends.DummyFileManager.acquire
   backends.DummyFileManager.acquire_context
   backends.DummyFileManager.close

   backends.BackendArray
   backends.BackendEntrypoint.guess_can_open
   backends.BackendEntrypoint.open_dataset

   core.indexing.IndexingSupport
   core.indexing.explicit_indexing_adapter
   core.indexing.BasicIndexer
   core.indexing.OuterIndexer
   core.indexing.VectorizedIndexer
   core.indexing.LazilyIndexedArray
   core.indexing.LazilyVectorizedIndexedArray

   conventions.decode_cf_variables

   coding.variables.UnsignedIntegerCoder
   coding.variables.CFMaskCoder
   coding.variables.CFScaleOffsetCoder

   coding.strings.CharacterArrayCoder
   coding.strings.EncodedStringCoder

   coding.times.CFTimedeltaCoder
   coding.times.CFDatetimeCoder
Gallery
=======

Here's a list of examples on how to use xarray. We will be adding more examples soon.
Contributions are highly welcomed and appreciated. So, if you are interested in contributing, please consult the
:doc:`contributing` guide.



Notebook Examples
-----------------

.. panels::
    :column: text-center col-lg-6 col-md-6 col-sm-12 col-xs-12 p-2
    :card: +my-2
    :img-top-cls: w-75 m-auto p-2
    :body: d-none

    ---
    :img-top: _static/thumbnails/toy-weather-data.png
    ++++
    .. link-button:: examples/weather-data
        :type: ref
        :text: Toy weather data
        :classes: btn-outline-dark btn-block stretched-link

    ---
    :img-top: _static/thumbnails/monthly-means.png
    ++++
    .. link-button:: examples/monthly-means
        :type: ref
        :text: Calculating Seasonal Averages from Timeseries of Monthly Means
        :classes: btn-outline-dark btn-block stretched-link

    ---
    :img-top: _static/thumbnails/area_weighted_temperature.png
    ++++
    .. link-button:: examples/area_weighted_temperature
        :type: ref
        :text: Compare weighted and unweighted mean temperature
        :classes: btn-outline-dark btn-block stretched-link

    ---
    :img-top: _static/thumbnails/multidimensional-coords.png
    ++++
    .. link-button:: examples/multidimensional-coords
        :type: ref
        :text: Working with Multidimensional Coordinates
        :classes: btn-outline-dark btn-block stretched-link

    ---
    :img-top: _static/thumbnails/visualization_gallery.png
    ++++
    .. link-button:: examples/visualization_gallery
        :type: ref
        :text: Visualization Gallery
        :classes: btn-outline-dark btn-block stretched-link

    ---
    :img-top: _static/thumbnails/ROMS_ocean_model.png
    ++++
    .. link-button:: examples/ROMS_ocean_model
        :type: ref
        :text: ROMS Ocean Model Example
        :classes: btn-outline-dark btn-block stretched-link

    ---
    :img-top: _static/thumbnails/ERA5-GRIB-example.png
    ++++
    .. link-button:: examples/ERA5-GRIB-example
        :type: ref
        :text: GRIB Data Example
        :classes: btn-outline-dark btn-block stretched-link

    ---
    :img-top: _static/dataset-diagram-square-logo.png
    ++++
    .. link-button:: examples/apply_ufunc_vectorize_1d
        :type: ref
        :text: Applying unvectorized functions with apply_ufunc
        :classes: btn-outline-dark btn-block stretched-link


.. toctree::
    :maxdepth: 1
    :hidden:

    examples/weather-data
    examples/monthly-means
    examples/area_weighted_temperature
    examples/multidimensional-coords
    examples/visualization_gallery
    examples/ROMS_ocean_model
    examples/ERA5-GRIB-example
    examples/apply_ufunc_vectorize_1d


External Examples
-----------------


.. panels::
    :column: text-center col-lg-6 col-md-6 col-sm-12 col-xs-12 p-2
    :card: +my-2
    :img-top-cls: w-75 m-auto p-2
    :body: d-none

    ---
    :img-top: _static/dataset-diagram-square-logo.png
    ++++
    .. link-button:: https://corteva.github.io/rioxarray/stable/examples/examples.html
        :type: url
        :text: Managing raster data with rioxarray
        :classes: btn-outline-dark btn-block stretched-link

    ---
    :img-top: https://avatars.githubusercontent.com/u/60833341?s=200&v=4
    ++++
    .. link-button:: http://gallery.pangeo.io/
        :type: url
        :text: Xarray and dask on the cloud with Pangeo
        :classes: btn-outline-dark btn-block stretched-link

    ---
    :img-top: _static/dataset-diagram-square-logo.png
    ++++
    .. link-button:: https://examples.dask.org/xarray.html
        :type: url
        :text: Xarray with Dask Arrays
        :classes: btn-outline-dark btn-block stretched-link

Tutorials and Videos
====================


Tutorials
----------

- `Xarray's Tutorials`_ repository
- The `UW eScience Institute's Geohackweek`_ tutorial on xarray for geospatial data scientists.
- `Nicolas Fauchereau's 2015 tutorial`_ on xarray for netCDF users.



Videos
-------

.. panels::
  :card: text-center

  ---
  Xdev Python Tutorial Seminar Series 2021 seminar introducing xarray (1 of 2) | Anderson Banihirwe
  ^^^
  .. raw:: html

    <iframe width="100%" src="https://www.youtube.com/embed/Ss4ryKukhi4" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

  ---
  Xarray's virtual tutorial | October 2020 | Anderson Banihirwe, Deepak Cherian, and Martin Durant
  ^^^
  .. raw:: html

    <iframe width="100%" src="https://www.youtube.com/embed/a339Q5F48UQ" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

  ---
  Xarray's Tutorial presented at the 2020 SciPy Conference | Joe Hamman, Ryan Abernathey,
  Deepak Cherian, and Stephan Hoyer
  ^^^
  .. raw:: html

    <iframe width="100%" src="https://www.youtube.com/embed/mecN-Ph_-78" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

  ---
  Scipy 2015 talk introducing xarray to a general audience | Stephan Hoyer
  ^^^
  .. raw:: html

    <iframe width="100%" src="https://www.youtube.com/embed/X0pAhJgySxk" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

  ---
  2015 Unidata Users Workshop talk and tutorial with (`with answers`_) introducing
  xarray to users familiar with netCDF | Stephan Hoyer
  ^^^
  .. raw:: html

    <iframe width="100%" src="https://www.youtube.com/embed/J9ypQOnt5l8" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

Books, Chapters and Articles
-----------------------------

- Stephan Hoyer and Joe Hamman's `Journal of Open Research Software paper`_ describing the xarray project.


.. _Xarray's Tutorials: https://xarray-contrib.github.io/xarray-tutorial/
.. _Journal of Open Research Software paper: http://doi.org/10.5334/jors.148
.. _UW eScience Institute's Geohackweek : https://geohackweek.github.io/nDarrays/
.. _tutorial: https://github.com/Unidata/unidata-users-workshop/blob/master/notebooks/xray-tutorial.ipynb
.. _with answers: https://github.com/Unidata/unidata-users-workshop/blob/master/notebooks/xray-tutorial-with-answers.ipynb
.. _Nicolas Fauchereau's 2015 tutorial: http://nbviewer.iPython.org/github/nicolasfauchereau/metocean/blob/master/notebooks/xray.ipynb
.. _roadmap:

Development roadmap
===================

Authors: Xarray developers

Date: September 7, 2021

Xarray is an open source Python library for labeled multidimensional
arrays and datasets.

Our philosophy
--------------

Why has xarray been successful? In our opinion:

-  Xarray does a great job of solving **specific use-cases** for
   multidimensional data analysis:

   -  The dominant use-case for xarray is for analysis of gridded
      dataset in the geosciences, e.g., as part of the
      `Pangeo <http://pangeo.io>`__ project.
   -  Xarray is also used more broadly in the physical sciences, where
      we've found the needs for analyzing multidimensional datasets are
      remarkably consistent (e.g., see
      `SunPy <https://github.com/sunpy/ndcube>`__ and
      `PlasmaPy <https://github.com/PlasmaPy/PlasmaPy/issues/59>`__).
   -  Finally, xarray is used in a variety of other domains, including
      finance, `probabilistic
      programming <https://arviz-devs.github.io/arviz/>`__ and
      genomics.

-  Xarray is also a **domain agnostic** solution:

   -  We focus on providing a flexible set of functionality related
      labeled multidimensional arrays, rather than solving particular
      problems.
   -  This facilitates collaboration between users with different needs,
      and helps us attract a broad community of contributors.
   -  Importantly, this retains flexibility, for use cases that don't
      fit particularly well into existing frameworks.

-  Xarray **integrates well** with other libraries in the scientific
   Python stack.

   -  We leverage first-class external libraries for core features of
      xarray (e.g., NumPy for ndarrays, pandas for indexing, dask for
      parallel computing)
   -  We expose our internal abstractions to users (e.g.,
      ``apply_ufunc()``), which facilitates extending xarray in various
      ways.

Together, these features have made xarray a first-class choice for
labeled multidimensional arrays in Python.

We want to double-down on xarray's strengths by making it an even more
flexible and powerful tool for multidimensional data analysis. We want
to continue to engage xarray's core geoscience users, and to also reach
out to new domains to learn from other successful data models like those
of `yt <https://yt-project.org>`__ or the `OLAP
cube <https://en.wikipedia.org/wiki/OLAP_cube>`__.

Specific needs
--------------

The user community has voiced a number specific needs related to how
xarray interfaces with domain specific problems. Xarray may not solve
all of these issues directly, but these areas provide opportunities for
xarray to provide better, more extensible, interfaces. Some examples of
these common needs are:

-  Non-regular grids (e.g., staggered and unstructured meshes).
-  Physical units.
-  Lazily computed arrays (e.g., for coordinate systems).
-  New file-formats.

Technical vision
----------------

We think the right approach to extending xarray's user community and the
usefulness of the project is to focus on improving key interfaces that
can be used externally to meet domain-specific needs.

We can generalize the community's needs into three main categories:

-  More flexible grids/indexing.
-  More flexible arrays/computing.
-  More flexible storage backends.
-  More flexible data structures.

Each of these are detailed further in the subsections below.

Flexible indexes
~~~~~~~~~~~~~~~~

.. note::
   Work on flexible grids and indexes is currently underway. See
   `GH Project #1 <https://github.com/pydata/xarray/projects/1>`__ for more detail.

Xarray currently keeps track of indexes associated with coordinates by
storing them in the form of a ``pandas.Index`` in special
``xarray.IndexVariable`` objects.

The limitations of this model became clear with the addition of
``pandas.MultiIndex`` support in xarray 0.9, where a single index
corresponds to multiple xarray variables. MultiIndex support is highly
useful, but xarray now has numerous special cases to check for
MultiIndex levels.

A cleaner model would be to elevate ``indexes`` to an explicit part of
xarray's data model, e.g., as attributes on the ``Dataset`` and
``DataArray`` classes. Indexes would need to be propagated along with
coordinates in xarray operations, but will no longer would need to have
a one-to-one correspondance with coordinate variables. Instead, an index
should be able to refer to multiple (possibly multidimensional)
coordinates that define it. See `GH
1603 <https://github.com/pydata/xarray/issues/1603>`__ for full details

Specific tasks:

-  Add an ``indexes`` attribute to ``xarray.Dataset`` and
   ``xarray.Dataset``, as dictionaries that map from coordinate names to
   xarray index objects.
-  Use the new index interface to write wrappers for ``pandas.Index``,
   ``pandas.MultiIndex`` and ``scipy.spatial.KDTree``.
-  Expose the interface externally to allow third-party libraries to
   implement custom indexing routines, e.g., for geospatial look-ups on
   the surface of the Earth.

In addition to the new features it directly enables, this clean up will
allow xarray to more easily implement some long-awaited features that
build upon indexing, such as groupby operations with multiple variables.

Flexible arrays
~~~~~~~~~~~~~~~

.. note::
   Work on flexible arrays is currently underway. See
   `GH Project #2 <https://github.com/pydata/xarray/projects/2>`__ for more detail.

Xarray currently supports wrapping multidimensional arrays defined by
NumPy, dask and to a limited-extent pandas. It would be nice to have
interfaces that allow xarray to wrap alternative N-D array
implementations, e.g.:

-  Arrays holding physical units.
-  Lazily computed arrays.
-  Other ndarray objects, e.g., sparse, xnd, xtensor.

Our strategy has been to pursue upstream improvements in NumPy (see
`NEP-22 <http://www.numpy.org/neps/nep-0022-ndarray-duck-typing-overview.html>`__)
for supporting a complete duck-typing interface using with NumPy's
higher level array API. Improvements in NumPy's support for custom data
types would also be highly useful for xarray users.

By pursuing these improvements in NumPy we hope to extend the benefits
to the full scientific Python community, and avoid tight coupling
between xarray and specific third-party libraries (e.g., for
implementing untis). This will allow xarray to maintain its domain
agnostic strengths.

We expect that we may eventually add some minimal interfaces in xarray
for features that we delegate to external array libraries (e.g., for
getting units and changing units). If we do add these features, we
expect them to be thin wrappers, with core functionality implemented by
third-party libraries.

Flexible storage
~~~~~~~~~~~~~~~~

.. note::
   Work on flexible storage backends is currently underway. See
   `GH Project #3 <https://github.com/pydata/xarray/projects/3>`__ for more detail.

The xarray backends module has grown in size and complexity. Much of
this growth has been "organic" and mostly to support incremental
additions to the supported backends. This has left us with a fragile
internal API that is difficult for even experienced xarray developers to
use. Moreover, the lack of a public facing API for building xarray
backends means that users can not easily build backend interface for
xarray in third-party libraries.

The idea of refactoring the backends API and exposing it to users was
originally proposed in `GH
1970 <https://github.com/pydata/xarray/issues/1970>`__. The idea would
be to develop a well tested and generic backend base class and
associated utilities for external use. Specific tasks for this
development would include:

-  Exposing an abstract backend for writing new storage systems.
-  Exposing utilities for features like automatic closing of files,
   LRU-caching and explicit/lazy indexing.
-  Possibly moving some infrequently used backends to third-party
   packages.

Flexible data structures
~~~~~~~~~~~~~~~~~~~~~~~~

Xarray provides two primary data structures, the ``xarray.DataArray`` and
the ``xarray.Dataset``. This section describes two possible data model
extensions.

Tree-like data structure
++++++++++++++++++++++++

.. note::
   Work on developing a hierarchical data structure in xarray is just
   beginning. See `Datatree <https://github.com/TomNicholas/datatree>`__
   for an early prototype.

Xarray’s highest-level object is currently an ``xarray.Dataset``, whose data
model echoes that of a single netCDF group. However real-world datasets are
often better represented by a collection of related Datasets. Particular common
examples include:

-  Multi-resolution datasets,
-  Collections of time series datasets with differing lengths,
-  Heterogeneous datasets comprising multiple different types of related
   observational or simulation data,
-  Bayesian workflows involving various statistical distributions over multiple
   variables,
-  Whole netCDF files containing multiple groups.
-  Comparison of output from many similar models (such as in the IPCC's Coupled Model Intercomparison Projects)

A new tree-like data structure which is essentially a structured hierarchical
collection of Datasets could represent these cases, and would instead map to
multiple netCDF groups (see `GH4118 <https://github.com/pydata/xarray/issues/4118>`__.).

Currently there are several libraries which have wrapped xarray in order to build
domain-specific data structures (e.g. `xarray-multiscale <https://github.com/JaneliaSciComp/xarray-multiscale>`__.),
but a general ``xarray.DataTree`` object would obviate the need for these and]
consolidate effort in a single domain-agnostic tool, much as xarray has already achieved.

Labeled array without coordinates
+++++++++++++++++++++++++++++++++

There is a need for a lightweight array structure with named dimensions for
convenient indexing and broadcasting. Xarray includes such a structure internally
(``xarray.Variable``). We want to factor out xarray's “Variable”  object into a
standalone package with minimal dependencies for integration with libraries that
don't want to inherit xarray's dependency on pandas (e.g. scikit-learn).
The new “Variable” class will follow established array protocols and the new
data-apis standard. It will be capable of wrapping multiple array-like objects
(e.g. NumPy, Dask, Sparse, Pint, CuPy, Pytorch). While “DataArray” fits some of
these requirements, it offers a more complex data model than is desired for
many applications and depends on pandas.

Engaging more users
-------------------

.. note::
   Work on improving xarray’s documentation and user engagement is
   currently underway. See `GH Project #4 <https://github.com/pydata/xarray/projects/4>`__
   for more detail.

Like many open-source projects, the documentation of xarray has grown
together with the library's features. While we think that the xarray
documentation is comprehensive already, we acknowledge that the adoption
of xarray might be slowed down because of the substantial time
investment required to learn its working principles. In particular,
non-computer scientists or users less familiar with the pydata ecosystem
might find it difficult to learn xarray and realize how xarray can help
them in their daily work.

In order to lower this adoption barrier, we propose to:

-  Develop entry-level tutorials for users with different backgrounds. For
   example, we would like to develop tutorials for users with or without
   previous knowledge of pandas, NumPy, netCDF, etc. These tutorials may be
   built as part of xarray's documentation or included in a separate repository
   to enable interactive use (e.g. mybinder.org).
-  Document typical user workflows in a dedicated website, following the example
   of `dask-stories
   <https://matthewrocklin.com/blog/work/2018/07/16/dask-stories>`__.
-  Write a basic glossary that defines terms that might not be familiar to all
   (e.g. "lazy", "labeled", "serialization", "indexing", "backend").


Administrative
--------------

NumFOCUS
~~~~~~~~

On July 16, 2018, Joe and Stephan submitted xarray's fiscal sponsorship
application to NumFOCUS.
.. currentmodule:: xarray

.. _howdoi:

How do I ...
============

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - How do I...
     - Solution
   * - add a DataArray to my dataset as a new variable
     - ``my_dataset[varname] = my_dataArray`` or :py:meth:`Dataset.assign` (see also :ref:`dictionary_like_methods`)
   * - add variables from other datasets to my dataset
     - :py:meth:`Dataset.merge`
   * - add a new dimension and/or coordinate
     - :py:meth:`DataArray.expand_dims`, :py:meth:`Dataset.expand_dims`
   * - add a new coordinate variable
     - :py:meth:`DataArray.assign_coords`
   * - change a data variable to a coordinate variable
     - :py:meth:`Dataset.set_coords`
   * - change the order of dimensions
     - :py:meth:`DataArray.transpose`, :py:meth:`Dataset.transpose`
   * - reshape dimensions
     - :py:meth:`DataArray.stack`, :py:meth:`Dataset.stack`, :py:meth:`Dataset.coarsen.construct`, :py:meth:`DataArray.coarsen.construct`
   * - remove a variable from my object
     - :py:meth:`Dataset.drop_vars`, :py:meth:`DataArray.drop_vars`
   * - remove dimensions of length 1 or 0
     - :py:meth:`DataArray.squeeze`, :py:meth:`Dataset.squeeze`
   * - remove all variables with a particular dimension
     - :py:meth:`Dataset.drop_dims`
   * - convert non-dimension coordinates to data variables or remove them
     - :py:meth:`DataArray.reset_coords`, :py:meth:`Dataset.reset_coords`
   * - rename a variable, dimension or coordinate
     - :py:meth:`Dataset.rename`, :py:meth:`DataArray.rename`, :py:meth:`Dataset.rename_vars`, :py:meth:`Dataset.rename_dims`,
   * - convert a DataArray to Dataset or vice versa
     - :py:meth:`DataArray.to_dataset`, :py:meth:`Dataset.to_array`, :py:meth:`Dataset.to_stacked_array`, :py:meth:`DataArray.to_unstacked_dataset`
   * - extract variables that have certain attributes
     - :py:meth:`Dataset.filter_by_attrs`
   * - extract the underlying array (e.g. NumPy or Dask arrays)
     - :py:attr:`DataArray.data`
   * - convert to and extract the underlying NumPy array
     - :py:attr:`DataArray.values`
   * - convert to a pandas DataFrame
     - :py:attr:`Dataset.to_dataframe`
   * - sort values
     - :py:attr:`Dataset.sortby`
   * - find out if my xarray object is wrapping a Dask Array
     - :py:func:`dask.is_dask_collection`
   * - know how much memory my object requires
     - :py:attr:`DataArray.nbytes`, :py:attr:`Dataset.nbytes`
   * - Get axis number for a dimension
     - :py:meth:`DataArray.get_axis_num`
   * - convert a possibly irregularly sampled timeseries to a regularly sampled timeseries
     - :py:meth:`DataArray.resample`, :py:meth:`Dataset.resample` (see :ref:`resampling` for more)
   * - apply a function on all data variables in a Dataset
     - :py:meth:`Dataset.map`
   * - write xarray objects with complex values to a netCDF file
     - :py:func:`Dataset.to_netcdf`, :py:func:`DataArray.to_netcdf` specifying ``engine="h5netcdf", invalid_netcdf=True``
   * - make xarray objects look like other xarray objects
     - :py:func:`~xarray.ones_like`, :py:func:`~xarray.zeros_like`, :py:func:`~xarray.full_like`, :py:meth:`Dataset.reindex_like`, :py:meth:`Dataset.interp_like`, :py:meth:`Dataset.broadcast_like`, :py:meth:`DataArray.reindex_like`, :py:meth:`DataArray.interp_like`, :py:meth:`DataArray.broadcast_like`
   * - Make sure my datasets have values at the same coordinate locations
     - ``xr.align(dataset_1, dataset_2, join="exact")``
   * - replace NaNs with other values
     - :py:meth:`Dataset.fillna`, :py:meth:`Dataset.ffill`, :py:meth:`Dataset.bfill`, :py:meth:`Dataset.interpolate_na`, :py:meth:`DataArray.fillna`, :py:meth:`DataArray.ffill`, :py:meth:`DataArray.bfill`, :py:meth:`DataArray.interpolate_na`
   * - extract the year, month, day or similar from a DataArray of time values
     - ``obj.dt.month`` for example where ``obj`` is a :py:class:`~xarray.DataArray` containing ``datetime64`` or ``cftime`` values. See :ref:`dt_accessor` for more.
   * - round off time values to a specified frequency
     - ``obj.dt.ceil``, ``obj.dt.floor``, ``obj.dt.round``. See :ref:`dt_accessor` for more.
   * - make a mask that is ``True`` where an object contains any of the values in a array
     - :py:meth:`Dataset.isin`, :py:meth:`DataArray.isin`
   * - Index using a boolean mask
     - :py:meth:`Dataset.query`, :py:meth:`DataArray.query`, :py:meth:`Dataset.where`, :py:meth:`DataArray.where`
   * - preserve ``attrs`` during (most) xarray operations
     - ``xr.set_options(keep_attrs=True)``
xarray: N-D labeled arrays and datasets in Python
=================================================

**xarray** (formerly **xray**) is an open source project and Python package
that makes working with labelled multi-dimensional arrays simple,
efficient, and fun!

Xarray introduces labels in the form of dimensions, coordinates and
attributes on top of raw NumPy_-like arrays, which allows for a more
intuitive, more concise, and less error-prone developer experience.
The package includes a large and growing library of domain-agnostic functions
for advanced analytics and visualization with these data structures.

Xarray is inspired by and borrows heavily from pandas_, the popular data
analysis package focused on labelled tabular data.
It is particularly tailored to working with netCDF_ files, which were the
source of xarray's data model, and integrates tightly with dask_ for parallel
computing.

.. _NumPy: http://www.numpy.org
.. _pandas: http://pandas.pydata.org
.. _dask: http://dask.org
.. _netCDF: http://www.unidata.ucar.edu/software/netcdf


.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: For users

   Getting Started <getting-started-guide/index>
   User Guide <user-guide/index>
   Gallery <gallery>
   Tutorials & Videos <tutorials-and-videos>
   API Reference <api>
   How do I ... <howdoi>
   Ecosystem <ecosystem>

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: For developers/contributors

   Contributing Guide <contributing>
   Xarray Internals <internals/index>
   Development Roadmap <roadmap>
   Team <team>
   Developers Meeting <developers-meeting>
   What’s New <whats-new>
   GitHub repository <https://github.com/pydata/xarray>

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Community

   GitHub discussions <https://github.com/pydata/xarray/discussions>
   StackOverflow <https://stackoverflow.com/questions/tagged/python-xarray>




Get in touch
------------

- If you have a question like "How do I concatenate a list of datasets?", ask on `GitHub discussions`_ or `StackOverflow`_.
  Please include a self-contained reproducible example if possible.
- Report bugs, suggest features or view the source code `on GitHub`_.
- For less well defined questions or ideas, or to announce other projects of
  interest to xarray users, use `GitHub discussions`_ or the `mailing list`_.

.. _StackOverFlow: https://stackoverflow.com/questions/tagged/python-xarray
.. _Github discussions: https://github.com/pydata/xarray/discussions
.. _mailing list: https://groups.google.com/forum/#!forum/xarray
.. _on GitHub: https://github.com/pydata/xarray

NumFOCUS
--------

.. image:: _static/numfocus_logo.png
   :scale: 50 %
   :target: https://numfocus.org/

Xarray is a fiscally sponsored project of NumFOCUS_, a nonprofit dedicated
to supporting the open source scientific computing community. If you like
Xarray and want to support our mission, please consider making a donation_
to support our efforts.

.. _donation: https://numfocus.salsalabs.org/donate-to-xarray/


History
-------

Xarray is an evolution of an internal tool developed at `The Climate
Corporation`__. It was originally written by Climate Corp researchers Stephan
Hoyer, Alex Kleeman and Eugene Brevdo and was released as open source in
May 2014. The project was renamed from "xray" in January 2016. Xarray became a
fiscally sponsored project of NumFOCUS_ in August 2018.

__ http://climate.com/
.. _NumFOCUS: https://numfocus.org

License
-------

Xarray is available under the open source `Apache License`__.

__ http://www.apache.org/licenses/LICENSE-2.0.html
Developers meeting
------------------

Xarray developers meet bi-weekly every other Wednesday.

The meeting occurs on `Zoom <https://us02web.zoom.us/j/88251613296?pwd=azZsSkU1UWJZTVFKNnhIUVdZcENUZz09>`__.

Notes for the meeting are kept `here <https://hackmd.io/@U4W-olO3TX-hc-cvbjNe4A/xarray-dev-meeting/edit>`__.

There is a `GitHub issue <https://github.com/pydata/xarray/issues/4001>`__ for changes to the meeting.

You can subscribe to this calendar to be notified of changes:

* `Google Calendar <https://calendar.google.com/calendar/embed?src=ucar.edu_2gjd5fuugcj4ol6ij7knj8krn8%40group.calendar.google.com&ctz=America%2FLos_Angeles>`__
* `iCal <https://calendar.google.com/calendar/ical/ucar.edu_2gjd5fuugcj4ol6ij7knj8krn8%40group.calendar.google.com/public/basic.ics>`__

.. raw:: html

   <iframe id="calendariframe" src="https://calendar.google.com/calendar/embed?ctz=local&amp;src=ucar.edu_2gjd5fuugcj4ol6ij7knj8krn8%40group.calendar.google.com" style="border: 0" width="800" height="600" frameborder="0" scrolling="no"></iframe>
   <script>document.getElementById("calendariframe").src = document.getElementById("calendariframe").src.replace("ctz=local", "ctz=" + Intl.DateTimeFormat().resolvedOptions().timeZone)</script>
:orphan:

xarray
------

You can find information about building the docs at our `Contributing page <http://xarray.pydata.org/en/latest/contributing.html#contributing-to-the-documentation>`_.
Team
-----

Current core developers
~~~~~~~~~~~~~~~~~~~~~~~

Xarray core developers are responsible for the ongoing organizational maintenance and technical direction of the xarray project.

The current core developers team comprises:

.. panels::
   :column: col-lg-4 col-md-4 col-sm-6 col-xs-12 p-2
   :card: text-center

   ---
   .. image:: https://avatars.githubusercontent.com/u/1217238?v=4
   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   :link-badge:`https://github.com/shoyer,"Stephan Hoyer",cls=btn badge-light`

   ---
   .. image:: https://avatars.githubusercontent.com/u/1197350?v=4
   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   :link-badge:`https://github.com/rabernat,"Ryan Abernathey",cls=btn badge-light`

   ---
   .. image:: https://avatars.githubusercontent.com/u/2443309?v=4
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   :link-badge:`https://github.com/jhamman,"Joe Hamman",cls=btn badge-light`

   ---
   .. image:: https://avatars.githubusercontent.com/u/4160723?v=4
   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   :link-badge:`https://github.com/benbovy,"Benoit Bovy",cls=btn badge-light`

   ---
   .. image:: https://avatars.githubusercontent.com/u/10050469?v=4
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   :link-badge:`https://github.com/fmaussion,"Fabien Maussion",cls=btn badge-light`

   ---
   .. image:: https://avatars.githubusercontent.com/u/6815844?v=4
   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   :link-badge:`https://github.com/fujiisoup,"Keisuke Fujii",cls=btn badge-light`

   ---
   .. image:: https://avatars.githubusercontent.com/u/5635139?v=4
   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   :link-badge:`https://github.com/max-sixty,"Maximilian Roos",cls=btn badge-light`

   ---
   .. image:: https://avatars.githubusercontent.com/u/2448579?v=4
   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   :link-badge:`https://github.com/dcherian,"Deepak Cherian",cls=btn badge-light`

   ---
   .. image:: https://avatars.githubusercontent.com/u/6628425?v=4
   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   :link-badge:`https://github.com/spencerkclark,"Spencer Clark",cls=btn badge-light`

   ---
   .. image:: https://avatars.githubusercontent.com/u/35968931?v=4
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   :link-badge:`https://github.com/TomNicholas,"Tom Nicholas",cls=btn badge-light`

   ---
   .. image:: https://avatars.githubusercontent.com/u/6213168?v=4
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   :link-badge:`https://github.com/crusaderky,"Guido Imperiale",cls=btn badge-light`

   ---
   .. image:: https://avatars.githubusercontent.com/u/14808389?v=4
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   :link-badge:`https://github.com/keewis,"Justus Magin",cls=btn badge-light`

   ---
   .. image:: https://avatars.githubusercontent.com/u/10194086?v=4
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   :link-badge:`https://github.com/mathause,"Mathias Hauser",cls=btn badge-light`

   ---
   .. image:: https://avatars.githubusercontent.com/u/13301940?v=4
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   :link-badge:`https://github.com/andersy005,"Anderson Banihirwe",cls=btn badge-light`

   ---
   .. image:: https://avatars.githubusercontent.com/u/14371165?v=4
   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   :link-badge:`https://github.com/Illviljan,"Jimmy Westling",cls=btn badge-light`


The full list of contributors is on our `GitHub Contributors Page <https://github.com/pydata/xarray/graphs/contributors>`__.
.. currentmodule:: xarray

.. _api:

#############
API reference
#############

This page provides an auto-generated summary of xarray's API. For more details
and examples, refer to the relevant chapters in the main part of the
documentation.

See also: :ref:`public api`

Top-level functions
===================

.. autosummary::
   :toctree: generated/

   apply_ufunc
   align
   broadcast
   concat
   merge
   combine_by_coords
   combine_nested
   where
   infer_freq
   full_like
   zeros_like
   ones_like
   cov
   corr
   cross
   dot
   polyval
   map_blocks
   show_versions
   set_options
   get_options
   unify_chunks

Dataset
=======

Creating a dataset
------------------

.. autosummary::
   :toctree: generated/

   Dataset
   decode_cf

Attributes
----------

.. autosummary::
   :toctree: generated/

   Dataset.dims
   Dataset.sizes
   Dataset.data_vars
   Dataset.coords
   Dataset.attrs
   Dataset.encoding
   Dataset.indexes
   Dataset.chunks
   Dataset.chunksizes
   Dataset.nbytes

Dictionary interface
--------------------

Datasets implement the mapping interface with keys given by variable names
and values given by ``DataArray`` objects.

.. autosummary::
   :toctree: generated/

   Dataset.__getitem__
   Dataset.__setitem__
   Dataset.__delitem__
   Dataset.update
   Dataset.get
   Dataset.items
   Dataset.keys
   Dataset.values

Dataset contents
----------------

.. autosummary::
   :toctree: generated/

   Dataset.copy
   Dataset.assign
   Dataset.assign_coords
   Dataset.assign_attrs
   Dataset.pipe
   Dataset.merge
   Dataset.rename
   Dataset.rename_vars
   Dataset.rename_dims
   Dataset.swap_dims
   Dataset.expand_dims
   Dataset.drop_vars
   Dataset.drop_dims
   Dataset.set_coords
   Dataset.reset_coords
   Dataset.convert_calendar
   Dataset.interp_calendar
   Dataset.get_index

Comparisons
-----------

.. autosummary::
   :toctree: generated/

   Dataset.equals
   Dataset.identical
   Dataset.broadcast_equals

Indexing
--------

.. autosummary::
   :toctree: generated/

   Dataset.loc
   Dataset.isel
   Dataset.sel
   Dataset.drop_sel
   Dataset.drop_isel
   Dataset.head
   Dataset.tail
   Dataset.thin
   Dataset.squeeze
   Dataset.interp
   Dataset.interp_like
   Dataset.reindex
   Dataset.reindex_like
   Dataset.set_index
   Dataset.reset_index
   Dataset.reorder_levels
   Dataset.query

Missing value handling
----------------------

.. autosummary::
   :toctree: generated/

   Dataset.isnull
   Dataset.notnull
   Dataset.combine_first
   Dataset.count
   Dataset.dropna
   Dataset.fillna
   Dataset.ffill
   Dataset.bfill
   Dataset.interpolate_na
   Dataset.where
   Dataset.isin

Computation
-----------

.. autosummary::
   :toctree: generated/

   Dataset.map
   Dataset.reduce
   Dataset.groupby
   Dataset.groupby_bins
   Dataset.rolling
   Dataset.rolling_exp
   Dataset.weighted
   Dataset.coarsen
   Dataset.resample
   Dataset.diff
   Dataset.quantile
   Dataset.differentiate
   Dataset.integrate
   Dataset.map_blocks
   Dataset.polyfit
   Dataset.curvefit

Aggregation
-----------

.. autosummary::
   :toctree: generated/

   Dataset.all
   Dataset.any
   Dataset.argmax
   Dataset.argmin
   Dataset.idxmax
   Dataset.idxmin
   Dataset.max
   Dataset.min
   Dataset.mean
   Dataset.median
   Dataset.prod
   Dataset.sum
   Dataset.std
   Dataset.var
   Dataset.cumsum
   Dataset.cumprod

ndarray methods
---------------

.. autosummary::
   :toctree: generated/

   Dataset.argsort
   Dataset.astype
   Dataset.clip
   Dataset.conj
   Dataset.conjugate
   Dataset.imag
   Dataset.round
   Dataset.real
   Dataset.rank

Reshaping and reorganizing
--------------------------

.. autosummary::
   :toctree: generated/

   Dataset.transpose
   Dataset.stack
   Dataset.unstack
   Dataset.to_stacked_array
   Dataset.shift
   Dataset.roll
   Dataset.pad
   Dataset.sortby
   Dataset.broadcast_like

Plotting
--------

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_method.rst

   Dataset.plot.scatter
   Dataset.plot.quiver
   Dataset.plot.streamplot

DataArray
=========

.. autosummary::
   :toctree: generated/

   DataArray

Attributes
----------

.. autosummary::
   :toctree: generated/

   DataArray.values
   DataArray.data
   DataArray.coords
   DataArray.dims
   DataArray.sizes
   DataArray.name
   DataArray.attrs
   DataArray.encoding
   DataArray.indexes
   DataArray.chunksizes

ndarray attributes
------------------

.. autosummary::
   :toctree: generated/

   DataArray.ndim
   DataArray.nbytes
   DataArray.shape
   DataArray.size
   DataArray.dtype
   DataArray.nbytes
   DataArray.chunks


DataArray contents
------------------

.. autosummary::
   :toctree: generated/

   DataArray.assign_coords
   DataArray.assign_attrs
   DataArray.pipe
   DataArray.rename
   DataArray.swap_dims
   DataArray.expand_dims
   DataArray.drop_vars
   DataArray.drop_duplicates
   DataArray.reset_coords
   DataArray.copy
   DataArray.convert_calendar
   DataArray.interp_calendar
   DataArray.get_index
   DataArray.astype
   DataArray.item

Indexing
--------

.. autosummary::
   :toctree: generated/

   DataArray.__getitem__
   DataArray.__setitem__
   DataArray.loc
   DataArray.isel
   DataArray.sel
   DataArray.drop_sel
   DataArray.drop_isel
   DataArray.head
   DataArray.tail
   DataArray.thin
   DataArray.squeeze
   DataArray.interp
   DataArray.interp_like
   DataArray.reindex
   DataArray.reindex_like
   DataArray.set_index
   DataArray.reset_index
   DataArray.reorder_levels
   DataArray.query

Missing value handling
----------------------

.. autosummary::
  :toctree: generated/

  DataArray.isnull
  DataArray.notnull
  DataArray.combine_first
  DataArray.count
  DataArray.dropna
  DataArray.fillna
  DataArray.ffill
  DataArray.bfill
  DataArray.interpolate_na
  DataArray.where
  DataArray.isin

Comparisons
-----------

.. autosummary::
   :toctree: generated/

   DataArray.equals
   DataArray.identical
   DataArray.broadcast_equals

Computation
-----------

.. autosummary::
   :toctree: generated/

   DataArray.reduce
   DataArray.groupby
   DataArray.groupby_bins
   DataArray.rolling
   DataArray.rolling_exp
   DataArray.weighted
   DataArray.coarsen
   DataArray.resample
   DataArray.get_axis_num
   DataArray.diff
   DataArray.dot
   DataArray.quantile
   DataArray.differentiate
   DataArray.integrate
   DataArray.polyfit
   DataArray.map_blocks
   DataArray.curvefit

Aggregation
-----------

.. autosummary::
   :toctree: generated/

   DataArray.all
   DataArray.any
   DataArray.argmax
   DataArray.argmin
   DataArray.idxmax
   DataArray.idxmin
   DataArray.max
   DataArray.min
   DataArray.mean
   DataArray.median
   DataArray.prod
   DataArray.sum
   DataArray.std
   DataArray.var
   DataArray.cumsum
   DataArray.cumprod

ndarray methods
---------------

.. autosummary::
   :toctree: generated/

   DataArray.argsort
   DataArray.clip
   DataArray.conj
   DataArray.conjugate
   DataArray.imag
   DataArray.searchsorted
   DataArray.round
   DataArray.real
   DataArray.T
   DataArray.rank


String manipulation
-------------------

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor.rst

   DataArray.str

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_method.rst

   DataArray.str.capitalize
   DataArray.str.casefold
   DataArray.str.cat
   DataArray.str.center
   DataArray.str.contains
   DataArray.str.count
   DataArray.str.decode
   DataArray.str.encode
   DataArray.str.endswith
   DataArray.str.extract
   DataArray.str.extractall
   DataArray.str.find
   DataArray.str.findall
   DataArray.str.format
   DataArray.str.get
   DataArray.str.get_dummies
   DataArray.str.index
   DataArray.str.isalnum
   DataArray.str.isalpha
   DataArray.str.isdecimal
   DataArray.str.isdigit
   DataArray.str.islower
   DataArray.str.isnumeric
   DataArray.str.isspace
   DataArray.str.istitle
   DataArray.str.isupper
   DataArray.str.join
   DataArray.str.len
   DataArray.str.ljust
   DataArray.str.lower
   DataArray.str.lstrip
   DataArray.str.match
   DataArray.str.normalize
   DataArray.str.pad
   DataArray.str.partition
   DataArray.str.repeat
   DataArray.str.replace
   DataArray.str.rfind
   DataArray.str.rindex
   DataArray.str.rjust
   DataArray.str.rpartition
   DataArray.str.rsplit
   DataArray.str.rstrip
   DataArray.str.slice
   DataArray.str.slice_replace
   DataArray.str.split
   DataArray.str.startswith
   DataArray.str.strip
   DataArray.str.swapcase
   DataArray.str.title
   DataArray.str.translate
   DataArray.str.upper
   DataArray.str.wrap
   DataArray.str.zfill

Datetimelike properties
-----------------------

**Datetime properties**:

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_attribute.rst

   DataArray.dt.year
   DataArray.dt.month
   DataArray.dt.day
   DataArray.dt.hour
   DataArray.dt.minute
   DataArray.dt.second
   DataArray.dt.microsecond
   DataArray.dt.nanosecond
   DataArray.dt.dayofweek
   DataArray.dt.weekday
   DataArray.dt.weekday_name
   DataArray.dt.dayofyear
   DataArray.dt.quarter
   DataArray.dt.days_in_month
   DataArray.dt.daysinmonth
   DataArray.dt.season
   DataArray.dt.time
   DataArray.dt.date
   DataArray.dt.calendar
   DataArray.dt.is_month_start
   DataArray.dt.is_month_end
   DataArray.dt.is_quarter_end
   DataArray.dt.is_year_start
   DataArray.dt.is_leap_year

**Datetime methods**:

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_method.rst

   DataArray.dt.floor
   DataArray.dt.ceil
   DataArray.dt.isocalendar
   DataArray.dt.round
   DataArray.dt.strftime

**Timedelta properties**:

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_attribute.rst

   DataArray.dt.days
   DataArray.dt.seconds
   DataArray.dt.microseconds
   DataArray.dt.nanoseconds

**Timedelta methods**:

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_method.rst

   DataArray.dt.floor
   DataArray.dt.ceil
   DataArray.dt.round


Reshaping and reorganizing
--------------------------

.. autosummary::
   :toctree: generated/

   DataArray.transpose
   DataArray.stack
   DataArray.unstack
   DataArray.to_unstacked_dataset
   DataArray.shift
   DataArray.roll
   DataArray.pad
   DataArray.sortby
   DataArray.broadcast_like

Plotting
--------

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_callable.rst

   DataArray.plot

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_method.rst

   DataArray.plot.contourf
   DataArray.plot.contour
   DataArray.plot.hist
   DataArray.plot.imshow
   DataArray.plot.line
   DataArray.plot.pcolormesh
   DataArray.plot.step
   DataArray.plot.surface

.. _api.ufuncs:

Universal functions
===================

.. warning::

   With recent versions of NumPy, Dask and xarray, NumPy ufuncs are now
   supported directly on all xarray and Dask objects. This obviates the need
   for the ``xarray.ufuncs`` module, which should not be used for new code
   unless compatibility with versions of NumPy prior to v1.13 is
   required. They will be removed once support for NumPy prior to
   v1.17 is dropped.

These functions are copied from NumPy, but extended to work on NumPy arrays,
dask arrays and all xarray objects. You can find them in the ``xarray.ufuncs``
module:

:py:attr:`~ufuncs.angle`
:py:attr:`~ufuncs.arccos`
:py:attr:`~ufuncs.arccosh`
:py:attr:`~ufuncs.arcsin`
:py:attr:`~ufuncs.arcsinh`
:py:attr:`~ufuncs.arctan`
:py:attr:`~ufuncs.arctan2`
:py:attr:`~ufuncs.arctanh`
:py:attr:`~ufuncs.ceil`
:py:attr:`~ufuncs.conj`
:py:attr:`~ufuncs.copysign`
:py:attr:`~ufuncs.cos`
:py:attr:`~ufuncs.cosh`
:py:attr:`~ufuncs.deg2rad`
:py:attr:`~ufuncs.degrees`
:py:attr:`~ufuncs.exp`
:py:attr:`~ufuncs.expm1`
:py:attr:`~ufuncs.fabs`
:py:attr:`~ufuncs.fix`
:py:attr:`~ufuncs.floor`
:py:attr:`~ufuncs.fmax`
:py:attr:`~ufuncs.fmin`
:py:attr:`~ufuncs.fmod`
:py:attr:`~ufuncs.fmod`
:py:attr:`~ufuncs.frexp`
:py:attr:`~ufuncs.hypot`
:py:attr:`~ufuncs.imag`
:py:attr:`~ufuncs.iscomplex`
:py:attr:`~ufuncs.isfinite`
:py:attr:`~ufuncs.isinf`
:py:attr:`~ufuncs.isnan`
:py:attr:`~ufuncs.isreal`
:py:attr:`~ufuncs.ldexp`
:py:attr:`~ufuncs.log`
:py:attr:`~ufuncs.log10`
:py:attr:`~ufuncs.log1p`
:py:attr:`~ufuncs.log2`
:py:attr:`~ufuncs.logaddexp`
:py:attr:`~ufuncs.logaddexp2`
:py:attr:`~ufuncs.logical_and`
:py:attr:`~ufuncs.logical_not`
:py:attr:`~ufuncs.logical_or`
:py:attr:`~ufuncs.logical_xor`
:py:attr:`~ufuncs.maximum`
:py:attr:`~ufuncs.minimum`
:py:attr:`~ufuncs.nextafter`
:py:attr:`~ufuncs.rad2deg`
:py:attr:`~ufuncs.radians`
:py:attr:`~ufuncs.real`
:py:attr:`~ufuncs.rint`
:py:attr:`~ufuncs.sign`
:py:attr:`~ufuncs.signbit`
:py:attr:`~ufuncs.sin`
:py:attr:`~ufuncs.sinh`
:py:attr:`~ufuncs.sqrt`
:py:attr:`~ufuncs.square`
:py:attr:`~ufuncs.tan`
:py:attr:`~ufuncs.tanh`
:py:attr:`~ufuncs.trunc`

IO / Conversion
===============

Dataset methods
---------------

.. autosummary::
   :toctree: generated/

   open_dataset
   load_dataset
   open_mfdataset
   open_rasterio
   open_zarr
   Dataset.to_netcdf
   Dataset.to_pandas
   Dataset.as_numpy
   Dataset.to_zarr
   save_mfdataset
   Dataset.to_array
   Dataset.to_dataframe
   Dataset.to_dask_dataframe
   Dataset.to_dict
   Dataset.from_dataframe
   Dataset.from_dict
   Dataset.close
   Dataset.compute
   Dataset.persist
   Dataset.load
   Dataset.chunk
   Dataset.unify_chunks
   Dataset.filter_by_attrs
   Dataset.info

DataArray methods
-----------------

.. autosummary::
   :toctree: generated/

   open_dataarray
   load_dataarray
   DataArray.to_dataset
   DataArray.to_netcdf
   DataArray.to_pandas
   DataArray.to_series
   DataArray.to_dataframe
   DataArray.to_numpy
   DataArray.as_numpy
   DataArray.to_index
   DataArray.to_masked_array
   DataArray.to_cdms2
   DataArray.to_iris
   DataArray.from_iris
   DataArray.to_dict
   DataArray.from_series
   DataArray.from_cdms2
   DataArray.from_dict
   DataArray.close
   DataArray.compute
   DataArray.persist
   DataArray.load
   DataArray.chunk
   DataArray.unify_chunks

Coordinates objects
===================

.. autosummary::
   :toctree: generated/

   core.coordinates.DataArrayCoordinates
   core.coordinates.DatasetCoordinates

GroupBy objects
===============

.. currentmodule:: xarray.core.groupby

Dataset
-------

.. autosummary::
   :toctree: generated/

   DatasetGroupBy
   DatasetGroupBy.map
   DatasetGroupBy.reduce
   DatasetGroupBy.assign
   DatasetGroupBy.assign_coords
   DatasetGroupBy.first
   DatasetGroupBy.last
   DatasetGroupBy.fillna
   DatasetGroupBy.quantile
   DatasetGroupBy.where
   DatasetGroupBy.all
   DatasetGroupBy.any
   DatasetGroupBy.count
   DatasetGroupBy.max
   DatasetGroupBy.mean
   DatasetGroupBy.median
   DatasetGroupBy.min
   DatasetGroupBy.prod
   DatasetGroupBy.std
   DatasetGroupBy.sum
   DatasetGroupBy.var
   DatasetGroupBy.dims
   DatasetGroupBy.groups

DataArray
---------

.. autosummary::
   :toctree: generated/

   DataArrayGroupBy
   DataArrayGroupBy.map
   DataArrayGroupBy.reduce
   DataArrayGroupBy.assign_coords
   DataArrayGroupBy.first
   DataArrayGroupBy.last
   DataArrayGroupBy.fillna
   DataArrayGroupBy.quantile
   DataArrayGroupBy.where
   DataArrayGroupBy.all
   DataArrayGroupBy.any
   DataArrayGroupBy.count
   DataArrayGroupBy.max
   DataArrayGroupBy.mean
   DataArrayGroupBy.median
   DataArrayGroupBy.min
   DataArrayGroupBy.prod
   DataArrayGroupBy.std
   DataArrayGroupBy.sum
   DataArrayGroupBy.var
   DataArrayGroupBy.dims
   DataArrayGroupBy.groups


Rolling objects
===============

.. currentmodule:: xarray.core.rolling

Dataset
-------

.. autosummary::
   :toctree: generated/

   DatasetRolling
   DatasetRolling.construct
   DatasetRolling.reduce
   DatasetRolling.argmax
   DatasetRolling.argmin
   DatasetRolling.count
   DatasetRolling.max
   DatasetRolling.mean
   DatasetRolling.median
   DatasetRolling.min
   DatasetRolling.prod
   DatasetRolling.std
   DatasetRolling.sum
   DatasetRolling.var

DataArray
---------

.. autosummary::
   :toctree: generated/

   DataArrayRolling
   DataArrayRolling.construct
   DataArrayRolling.reduce
   DataArrayRolling.argmax
   DataArrayRolling.argmin
   DataArrayRolling.count
   DataArrayRolling.max
   DataArrayRolling.mean
   DataArrayRolling.median
   DataArrayRolling.min
   DataArrayRolling.prod
   DataArrayRolling.std
   DataArrayRolling.sum
   DataArrayRolling.var

Coarsen objects
===============

Dataset
-------

.. autosummary::
   :toctree: generated/

   DatasetCoarsen
   DatasetCoarsen.all
   DatasetCoarsen.any
   DatasetCoarsen.construct
   DatasetCoarsen.count
   DatasetCoarsen.max
   DatasetCoarsen.mean
   DatasetCoarsen.median
   DatasetCoarsen.min
   DatasetCoarsen.prod
   DatasetCoarsen.reduce
   DatasetCoarsen.std
   DatasetCoarsen.sum
   DatasetCoarsen.var

DataArray
---------

.. autosummary::
   :toctree: generated/

   DataArrayCoarsen
   DataArrayCoarsen.all
   DataArrayCoarsen.any
   DataArrayCoarsen.construct
   DataArrayCoarsen.count
   DataArrayCoarsen.max
   DataArrayCoarsen.mean
   DataArrayCoarsen.median
   DataArrayCoarsen.min
   DataArrayCoarsen.prod
   DataArrayCoarsen.reduce
   DataArrayCoarsen.std
   DataArrayCoarsen.sum
   DataArrayCoarsen.var

Exponential rolling objects
===========================

.. currentmodule:: xarray.core.rolling_exp

.. autosummary::
   :toctree: generated/

   RollingExp
   RollingExp.mean
   RollingExp.sum

Weighted objects
================

.. currentmodule:: xarray.core.weighted

Dataset
-------

.. autosummary::
   :toctree: generated/

   DatasetWeighted
   DatasetWeighted.mean
   DatasetWeighted.sum
   DatasetWeighted.std
   DatasetWeighted.var
   DatasetWeighted.sum_of_weights
   DatasetWeighted.sum_of_squares

DataArray
---------

.. autosummary::
   :toctree: generated/

   DataArrayWeighted
   DataArrayWeighted.mean
   DataArrayWeighted.sum
   DataArrayWeighted.std
   DataArrayWeighted.var
   DataArrayWeighted.sum_of_weights
   DataArrayWeighted.sum_of_squares

Resample objects
================

.. currentmodule:: xarray.core.resample

Dataset
-------

.. autosummary::
   :toctree: generated/

   DatasetResample
   DatasetResample.asfreq
   DatasetResample.backfill
   DatasetResample.interpolate
   DatasetResample.nearest
   DatasetResample.pad
   DatasetResample.all
   DatasetResample.any
   DatasetResample.apply
   DatasetResample.assign
   DatasetResample.assign_coords
   DatasetResample.bfill
   DatasetResample.count
   DatasetResample.ffill
   DatasetResample.fillna
   DatasetResample.first
   DatasetResample.last
   DatasetResample.map
   DatasetResample.max
   DatasetResample.mean
   DatasetResample.median
   DatasetResample.min
   DatasetResample.prod
   DatasetResample.quantile
   DatasetResample.reduce
   DatasetResample.std
   DatasetResample.sum
   DatasetResample.var
   DatasetResample.where
   DatasetResample.dims
   DatasetResample.groups


DataArray
---------

.. autosummary::
   :toctree: generated/

   DataArrayResample
   DataArrayResample.asfreq
   DataArrayResample.backfill
   DataArrayResample.interpolate
   DataArrayResample.nearest
   DataArrayResample.pad
   DataArrayResample.all
   DataArrayResample.any
   DataArrayResample.apply
   DataArrayResample.assign_coords
   DataArrayResample.bfill
   DataArrayResample.count
   DataArrayResample.ffill
   DataArrayResample.fillna
   DataArrayResample.first
   DataArrayResample.last
   DataArrayResample.map
   DataArrayResample.max
   DataArrayResample.mean
   DataArrayResample.median
   DataArrayResample.min
   DataArrayResample.prod
   DataArrayResample.quantile
   DataArrayResample.reduce
   DataArrayResample.std
   DataArrayResample.sum
   DataArrayResample.var
   DataArrayResample.where
   DataArrayResample.dims
   DataArrayResample.groups

Accessors
=========

.. currentmodule:: xarray

.. autosummary::
   :toctree: generated/

   core.accessor_dt.DatetimeAccessor
   core.accessor_dt.TimedeltaAccessor
   core.accessor_str.StringAccessor

Custom Indexes
==============
.. autosummary::
   :toctree: generated/

   CFTimeIndex

Creating custom indexes
-----------------------
.. autosummary::
   :toctree: generated/

   cftime_range
   date_range
   date_range_like

Faceting
--------
.. autosummary::
   :toctree: generated/

   plot.FacetGrid
   plot.FacetGrid.add_colorbar
   plot.FacetGrid.add_legend
   plot.FacetGrid.add_quiverkey
   plot.FacetGrid.map
   plot.FacetGrid.map_dataarray
   plot.FacetGrid.map_dataarray_line
   plot.FacetGrid.map_dataset
   plot.FacetGrid.set_axis_labels
   plot.FacetGrid.set_ticks
   plot.FacetGrid.set_titles
   plot.FacetGrid.set_xlabels
   plot.FacetGrid.set_ylabels

Tutorial
========

.. autosummary::
   :toctree: generated/

   tutorial.open_dataset
   tutorial.open_rasterio
   tutorial.load_dataset

Testing
=======

.. autosummary::
   :toctree: generated/

   testing.assert_equal
   testing.assert_identical
   testing.assert_allclose
   testing.assert_chunks_equal

Exceptions
==========

.. autosummary::
   :toctree: generated/

   MergeError
   SerializationWarning

Advanced API
============

.. autosummary::
   :toctree: generated/

   Dataset.variables
   DataArray.variable
   Variable
   IndexVariable
   as_variable
   Context
   register_dataset_accessor
   register_dataarray_accessor
   Dataset.set_close
   backends.BackendArray
   backends.BackendEntrypoint

These backends provide a low-level interface for lazily loading data from
external file-formats or protocols, and can be manually invoked to create
arguments for the ``load_store`` and ``dump_to_store`` Dataset methods:

.. autosummary::
   :toctree: generated/

   backends.NetCDF4DataStore
   backends.H5NetCDFStore
   backends.PydapDataStore
   backends.ScipyDataStore
   backends.FileManager
   backends.CachingFileManager
   backends.DummyFileManager

Deprecated / Pending Deprecation
================================

.. autosummary::
   :toctree: generated/

   Dataset.drop
   DataArray.drop
   Dataset.apply
   core.groupby.DataArrayGroupBy.apply
   core.groupby.DatasetGroupBy.apply

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_attribute.rst

   DataArray.dt.weekofyear
   DataArray.dt.week
{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module.split('.')[0] }}

.. autoaccessorcallable:: {{ (module.split('.')[1:] + [objname]) | join('.') }}.__call__
{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module.split('.')[0] }}

.. autoaccessormethod:: {{ (module.split('.')[1:] + [objname]) | join('.') }}
{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module.split('.')[0] }}

.. autoaccessorattribute:: {{ (module.split('.')[1:] + [objname]) | join('.') }}
{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module.split('.')[0] }}

.. autoaccessor:: {{ (module.split('.')[1:] + [objname]) | join('.') }}
##############
Quick overview
##############

Here are some quick examples of what you can do with :py:class:`xarray.DataArray`
objects. Everything is explained in much more detail in the rest of the
documentation.

To begin, import numpy, pandas and xarray using their customary abbreviations:

.. ipython:: python

    import numpy as np
    import pandas as pd
    import xarray as xr

Create a DataArray
------------------

You can make a DataArray from scratch by supplying data in the form of a numpy
array or list, with optional *dimensions* and *coordinates*:

.. ipython:: python

    data = xr.DataArray(np.random.randn(2, 3), dims=("x", "y"), coords={"x": [10, 20]})
    data

In this case, we have generated a 2D array, assigned the names *x* and *y* to the two dimensions respectively and associated two *coordinate labels* '10' and '20' with the two locations along the x dimension. If you supply a pandas :py:class:`~pandas.Series` or :py:class:`~pandas.DataFrame`, metadata is copied directly:

.. ipython:: python

    xr.DataArray(pd.Series(range(3), index=list("abc"), name="foo"))

Here are the key properties for a ``DataArray``:

.. ipython:: python

    # like in pandas, values is a numpy array that you can modify in-place
    data.values
    data.dims
    data.coords
    # you can use this dictionary to store arbitrary metadata
    data.attrs


Indexing
--------

Xarray supports four kinds of indexing. Since we have assigned coordinate labels to the x dimension we can use label-based indexing along that dimension just like pandas. The four examples below all yield the same result (the value at `x=10`) but at varying levels of convenience and intuitiveness.

.. ipython:: python

    # positional and by integer label, like numpy
    data[0, :]

    # loc or "location": positional and coordinate label, like pandas
    data.loc[10]

    # isel or "integer select":  by dimension name and integer label
    data.isel(x=0)

    # sel or "select": by dimension name and coordinate label
    data.sel(x=10)


Unlike positional indexing, label-based indexing frees us from having to know how our array is organized. All we need to know are the dimension name and the label we wish to index i.e. ``data.sel(x=10)`` works regardless of whether ``x`` is the first or second dimension of the array and regardless of whether ``10`` is the first or second element of ``x``. We have already told xarray that x is the first dimension when we created ``data``: xarray keeps track of this so we don't have to. For more, see :ref:`indexing`.


Attributes
----------

While you're setting up your DataArray, it's often a good idea to set metadata attributes. A useful choice is to set ``data.attrs['long_name']`` and ``data.attrs['units']`` since xarray will use these, if present, to automatically label your plots. These special names were chosen following the `NetCDF Climate and Forecast (CF) Metadata Conventions <http://cfconventions.org/cf-conventions/cf-conventions.html>`_. ``attrs`` is just a Python dictionary, so you can assign anything you wish.

.. ipython:: python

    data.attrs["long_name"] = "random velocity"
    data.attrs["units"] = "metres/sec"
    data.attrs["description"] = "A random variable created as an example."
    data.attrs["random_attribute"] = 123
    data.attrs
    # you can add metadata to coordinates too
    data.x.attrs["units"] = "x units"


Computation
-----------

Data arrays work very similarly to numpy ndarrays:

.. ipython:: python

    data + 10
    np.sin(data)
    # transpose
    data.T
    data.sum()

However, aggregation operations can use dimension names instead of axis
numbers:

.. ipython:: python

    data.mean(dim="x")

Arithmetic operations broadcast based on dimension name. This means you don't
need to insert dummy dimensions for alignment:

.. ipython:: python

    a = xr.DataArray(np.random.randn(3), [data.coords["y"]])
    b = xr.DataArray(np.random.randn(4), dims="z")

    a
    b

    a + b

It also means that in most cases you do not need to worry about the order of
dimensions:

.. ipython:: python

    data - data.T

Operations also align based on index labels:

.. ipython:: python

    data[:-1] - data[:1]

For more, see :ref:`comput`.

GroupBy
-------

Xarray supports grouped operations using a very similar API to pandas (see :ref:`groupby`):

.. ipython:: python

    labels = xr.DataArray(["E", "F", "E"], [data.coords["y"]], name="labels")
    labels
    data.groupby(labels).mean("y")
    data.groupby(labels).map(lambda x: x - x.min())

Plotting
--------

Visualizing your datasets is quick and convenient:

.. ipython:: python

    @savefig plotting_quick_overview.png
    data.plot()

Note the automatic labeling with names and units. Our effort in adding metadata attributes has paid off! Many aspects of these figures are customizable: see :ref:`plotting`.

pandas
------

Xarray objects can be easily converted to and from pandas objects using the :py:meth:`~xarray.DataArray.to_series`, :py:meth:`~xarray.DataArray.to_dataframe` and :py:meth:`~pandas.DataFrame.to_xarray` methods:

.. ipython:: python

    series = data.to_series()
    series

    # convert back
    series.to_xarray()

Datasets
--------

:py:class:`xarray.Dataset` is a dict-like container of aligned ``DataArray``
objects. You can think of it as a multi-dimensional generalization of the
:py:class:`pandas.DataFrame`:

.. ipython:: python

    ds = xr.Dataset(dict(foo=data, bar=("x", [1, 2]), baz=np.pi))
    ds


This creates a dataset with three DataArrays named ``foo``, ``bar`` and ``baz``. Use dictionary or dot indexing to pull out ``Dataset`` variables as ``DataArray`` objects but note that assignment only works with dictionary indexing:

.. ipython:: python

    ds["foo"]
    ds.foo


When creating ``ds``, we specified that ``foo`` is identical to ``data`` created earlier, ``bar`` is one-dimensional with single dimension ``x`` and associated values '1' and '2', and ``baz`` is a scalar not associated with any dimension in ``ds``. Variables in datasets can have different ``dtype`` and even different dimensions, but all dimensions are assumed to refer to points in the same shared coordinate system i.e. if two variables have dimension ``x``, that dimension must be identical in both variables.

For example, when creating ``ds`` xarray automatically *aligns* ``bar`` with ``DataArray`` ``foo``, i.e., they share the same coordinate system so that ``ds.bar['x'] == ds.foo['x'] == ds['x']``. Consequently, the following works without explicitly specifying the coordinate ``x`` when creating ``ds['bar']``:

.. ipython:: python

    ds.bar.sel(x=10)



You can do almost everything you can do with ``DataArray`` objects with
``Dataset`` objects (including indexing and arithmetic) if you prefer to work
with multiple variables at once.

Read & write netCDF files
-------------------------

NetCDF is the recommended file format for xarray objects. Users
from the geosciences will recognize that the :py:class:`~xarray.Dataset` data
model looks very similar to a netCDF file (which, in fact, inspired it).

You can directly read and write xarray objects to disk using :py:meth:`~xarray.Dataset.to_netcdf`, :py:func:`~xarray.open_dataset` and
:py:func:`~xarray.open_dataarray`:

.. ipython:: python

    ds.to_netcdf("example.nc")
    xr.open_dataset("example.nc")

.. ipython:: python
    :suppress:

    import os

    os.remove("example.nc")


It is common for datasets to be distributed across multiple files (commonly one file per timestep). Xarray supports this use-case by providing the :py:meth:`~xarray.open_mfdataset` and the :py:meth:`~xarray.save_mfdataset` methods. For more, see :ref:`io`.
.. _installing:

Installation
============

Required dependencies
---------------------

- Python (3.8 or later)
- `numpy <https://www.numpy.org/>`__ (1.18 or later)
- `pandas <https://pandas.pydata.org/>`__ (1.1 or later)

.. _optional-dependencies:

Optional dependencies
---------------------

.. note::

  If you are using pip to install xarray, optional dependencies can be installed by
  specifying *extras*. :ref:`installation-instructions` for both pip and conda
  are given below.

For netCDF and IO
~~~~~~~~~~~~~~~~~

- `netCDF4 <https://github.com/Unidata/netcdf4-python>`__: recommended if you
  want to use xarray for reading or writing netCDF files
- `scipy <http://scipy.org/>`__: used as a fallback for reading/writing netCDF3
- `pydap <http://www.pydap.org/>`__: used as a fallback for accessing OPeNDAP
- `h5netcdf <https://github.com/shoyer/h5netcdf>`__: an alternative library for
  reading and writing netCDF4 files that does not use the netCDF-C libraries
- `PyNIO <https://www.pyngl.ucar.edu/Nio.shtml>`__: for reading GRIB and other
  geoscience specific file formats. Note that PyNIO is not available for Windows and
  that the PyNIO backend may be moved outside of xarray in the future.
- `zarr <http://zarr.readthedocs.io/>`__: for chunked, compressed, N-dimensional arrays.
- `cftime <https://unidata.github.io/cftime>`__: recommended if you
  want to encode/decode datetimes for non-standard calendars or dates before
  year 1678 or after year 2262.
- `PseudoNetCDF <http://github.com/barronh/pseudonetcdf/>`__: recommended
  for accessing CAMx, GEOS-Chem (bpch), NOAA ARL files, ICARTT files
  (ffi1001) and many other.
- `rasterio <https://github.com/mapbox/rasterio>`__: for reading GeoTiffs and
  other gridded raster datasets.
- `iris <https://github.com/scitools/iris>`__: for conversion to and from iris'
  Cube objects
- `cfgrib <https://github.com/ecmwf/cfgrib>`__: for reading GRIB files via the
  *ECMWF ecCodes* library.

For accelerating xarray
~~~~~~~~~~~~~~~~~~~~~~~

- `scipy <http://scipy.org/>`__: necessary to enable the interpolation features for
  xarray objects
- `bottleneck <https://github.com/pydata/bottleneck>`__: speeds up
  NaN-skipping and rolling window aggregations by a large factor
- `numbagg <https://github.com/shoyer/numbagg>`_: for exponential rolling
  window operations

For parallel computing
~~~~~~~~~~~~~~~~~~~~~~

- `dask.array <http://dask.pydata.org>`__: required for :ref:`dask`.

For plotting
~~~~~~~~~~~~

- `matplotlib <http://matplotlib.org/>`__: required for :ref:`plotting`
- `cartopy <http://scitools.org.uk/cartopy/>`__: recommended for :ref:`plot-maps`
- `seaborn <http://seaborn.pydata.org/>`__: for better
  color palettes
- `nc-time-axis <https://github.com/SciTools/nc-time-axis>`__: for plotting
  cftime.datetime objects

Alternative data containers
~~~~~~~~~~~~~~~~~~~~~~~~~~~
- `sparse <https://sparse.pydata.org/>`_: for sparse arrays
- `pint <https://pint.readthedocs.io/>`_: for units of measure
- Any numpy-like objects that support
  `NEP-18 <https://numpy.org/neps/nep-0018-array-function-protocol.html>`_.
  Note that while such libraries theoretically should work, they are untested.
  Integration tests are in the process of being written for individual libraries.


.. _mindeps_policy:

Minimum dependency versions
---------------------------
Xarray adopts a rolling policy regarding the minimum supported version of its
dependencies:

- **Python:** 24 months
  (`NEP-29 <https://numpy.org/neps/nep-0029-deprecation_policy.html>`_)
- **numpy:** 18 months
  (`NEP-29 <https://numpy.org/neps/nep-0029-deprecation_policy.html>`_)
- **all other libraries:** 12 months

This means the latest minor (X.Y) version from N months prior. Patch versions (x.y.Z)
are not pinned, and only the latest available at the moment of publishing the xarray
release is guaranteed to work.

You can see the actual minimum tested versions:

`<https://github.com/pydata/xarray/blob/main/ci/requirements/py38-min-all-deps.yml>`_

.. _installation-instructions:

Instructions
------------

Xarray itself is a pure Python package, but its dependencies are not. The
easiest way to get everything installed is to use conda_. To install xarray
with its recommended dependencies using the conda command line tool::

    $ conda install -c conda-forge xarray dask netCDF4 bottleneck

.. _conda: http://conda.io/

If you require other :ref:`optional-dependencies` add them to the line above.

We recommend using the community maintained `conda-forge <https://conda-forge.github.io/>`__ channel,
as some of the dependencies are difficult to build. New releases may also appear in conda-forge before
being updated in the default channel.

If you don't use conda, be sure you have the required dependencies (numpy and
pandas) installed first. Then, install xarray with pip::

    $ python -m pip install xarray

We also maintain other dependency sets for different subsets of functionality::

    $ python -m pip install "xarray[io]"        # Install optional dependencies for handling I/O
    $ python -m pip install "xarray[accel]"     # Install optional dependencies for accelerating xarray
    $ python -m pip install "xarray[parallel]"  # Install optional dependencies for dask arrays
    $ python -m pip install "xarray[viz]"       # Install optional dependencies for visualization
    $ python -m pip install "xarray[complete]"  # Install all the above

The above commands should install most of the `optional dependencies`_. However,
some packages which are either not listed on PyPI or require extra
installation steps are excluded. To know which dependencies would be
installed, take a look at the ``[options.extras_require]`` section in
``setup.cfg``:

.. literalinclude:: ../../setup.cfg
   :language: ini
   :start-at: [options.extras_require]
   :end-before: [options.package_data]


Testing
-------

To run the test suite after installing xarray, install (via pypi or conda) `py.test
<https://pytest.org>`__ and run ``pytest`` in the root directory of the xarray
repository.


Performance Monitoring
~~~~~~~~~~~~~~~~~~~~~~

..
   TODO: uncomment once we have a working setup
         see https://github.com/pydata/xarray/pull/5066

   A fixed-point performance monitoring of (a part of) our code can be seen on
   `this page <https://pandas.pydata.org/speed/xarray/>`__.

To run these benchmark tests in a local machine, first install

- `airspeed-velocity <https://asv.readthedocs.io/en/latest/>`__: a tool for benchmarking
  Python packages over their lifetime.

and run
``asv run  # this will install some conda environments in ./.asv/envs``
Overview: Why xarray?
=====================

Xarray introduces labels in the form of dimensions, coordinates and attributes on top of
raw NumPy-like multidimensional arrays, which allows for a more intuitive, more concise,
and less error-prone developer experience.

What labels enable
------------------

Multi-dimensional (a.k.a. N-dimensional, ND) arrays (sometimes called
"tensors") are an essential part of computational science.
They are encountered in a wide range of fields, including physics, astronomy,
geoscience, bioinformatics, engineering, finance, and deep learning.
In Python, NumPy_ provides the fundamental data structure and API for
working with raw ND arrays.
However, real-world datasets are usually more than just raw numbers;
they have labels which encode information about how the array values map
to locations in space, time, etc.

Xarray doesn't just keep track of labels on arrays -- it uses them to provide a
powerful and concise interface. For example:

-  Apply operations over dimensions by name: ``x.sum('time')``.
-  Select values by label (or logical location) instead of integer location:
   ``x.loc['2014-01-01']`` or ``x.sel(time='2014-01-01')``.
-  Mathematical operations (e.g., ``x - y``) vectorize across multiple
   dimensions (array broadcasting) based on dimension names, not shape.
-  Easily use the `split-apply-combine <https://vita.had.co.nz/papers/plyr.pdf>`_
   paradigm with ``groupby``:
   ``x.groupby('time.dayofyear').mean()``.
-  Database-like alignment based on coordinate labels that smoothly
   handles missing values: ``x, y = xr.align(x, y, join='outer')``.
-  Keep track of arbitrary metadata in the form of a Python dictionary:
   ``x.attrs``.

The N-dimensional nature of xarray's data structures makes it suitable for dealing
with multi-dimensional scientific data, and its use of dimension names
instead of axis labels (``dim='time'`` instead of ``axis=0``) makes such
arrays much more manageable than the raw numpy ndarray: with xarray, you don't
need to keep track of the order of an array's dimensions or insert dummy dimensions of
size 1 to align arrays (e.g., using ``np.newaxis``).

The immediate payoff of using xarray is that you'll write less code. The
long-term payoff is that you'll understand what you were thinking when you come
back to look at it weeks or months later.

Core data structures
--------------------

Xarray has two core data structures, which build upon and extend the core
strengths of NumPy_ and pandas_. Both data structures are fundamentally N-dimensional:

- :py:class:`~xarray.DataArray` is our implementation of a labeled, N-dimensional
  array. It is an N-D generalization of a :py:class:`pandas.Series`. The name
  ``DataArray`` itself is borrowed from Fernando Perez's datarray_ project,
  which prototyped a similar data structure.
- :py:class:`~xarray.Dataset` is a multi-dimensional, in-memory array database.
  It is a dict-like container of ``DataArray`` objects aligned along any number of
  shared dimensions, and serves a similar purpose in xarray to the
  :py:class:`pandas.DataFrame`.

The value of attaching labels to numpy's :py:class:`numpy.ndarray` may be
fairly obvious, but the dataset may need more motivation.

The power of the dataset over a plain dictionary is that, in addition to
pulling out arrays by name, it is possible to select or combine data along a
dimension across all arrays simultaneously. Like a
:py:class:`~pandas.DataFrame`, datasets facilitate array operations with
heterogeneous data -- the difference is that the arrays in a dataset can have
not only different data types, but also different numbers of dimensions.

This data model is borrowed from the netCDF_ file format, which also provides
xarray with a natural and portable serialization format. NetCDF is very popular
in the geosciences, and there are existing libraries for reading and writing
netCDF in many programming languages, including Python.

Xarray distinguishes itself from many tools for working with netCDF data
in-so-far as it provides data structures for in-memory analytics that both
utilize and preserve labels. You only need to do the tedious work of adding
metadata once, not every time you save a file.

Goals and aspirations
---------------------

Xarray contributes domain-agnostic data-structures and tools for labeled
multi-dimensional arrays to Python's SciPy_ ecosystem for numerical computing.
In particular, xarray builds upon and integrates with NumPy_ and pandas_:

- Our user-facing interfaces aim to be more explicit versions of those found in
  NumPy/pandas.
- Compatibility with the broader ecosystem is a major goal: it should be easy
  to get your data in and out.
- We try to keep a tight focus on functionality and interfaces related to
  labeled data, and leverage other Python libraries for everything else, e.g.,
  NumPy/pandas for fast arrays/indexing (xarray itself contains no compiled
  code), Dask_ for parallel computing, matplotlib_ for plotting, etc.

Xarray is a collaborative and community driven project, run entirely on
volunteer effort (see :ref:`contributing`).
Our target audience is anyone who needs N-dimensional labeled arrays in Python.
Originally, development was driven by the data analysis needs of physical
scientists (especially geoscientists who already know and love
netCDF_), but it has become a much more broadly useful tool, and is still
under active development.
See our technical :ref:`roadmap` for more details, and feel free to reach out
with questions about whether xarray is the right tool for your needs.

.. _datarray: https://github.com/fperez/datarray
.. _Dask: http://dask.org
.. _matplotlib: http://matplotlib.org
.. _netCDF: http://www.unidata.ucar.edu/software/netcdf
.. _NumPy: http://www.numpy.org
.. _pandas: http://pandas.pydata.org
.. _SciPy: http://www.scipy.org
.. _faq:

Frequently Asked Questions
==========================

.. ipython:: python
    :suppress:

    import numpy as np
    import pandas as pd
    import xarray as xr

    np.random.seed(123456)


Your documentation keeps mentioning pandas. What is pandas?
-----------------------------------------------------------

pandas_ is a very popular data analysis package in Python
with wide usage in many fields. Our API is heavily inspired by pandas —
this is why there are so many references to pandas.

.. _pandas: https://pandas.pydata.org


Do I need to know pandas to use xarray?
---------------------------------------

No! Our API is heavily inspired by pandas so while knowing pandas will let you
become productive more quickly, knowledge of pandas is not necessary to use xarray.


Should I use xarray instead of pandas?
--------------------------------------

It's not an either/or choice! xarray provides robust support for converting
back and forth between the tabular data-structures of pandas and its own
multi-dimensional data-structures.

That said, you should only bother with xarray if some aspect of data is
fundamentally multi-dimensional. If your data is unstructured or
one-dimensional, pandas is usually the right choice: it has better performance
for common operations such as ``groupby`` and you'll find far more usage
examples online.


Why is pandas not enough?
-------------------------

pandas is a fantastic library for analysis of low-dimensional labelled data -
if it can be sensibly described as "rows and columns", pandas is probably the
right choice.  However, sometimes we want to use higher dimensional arrays
(`ndim > 2`), or arrays for which the order of dimensions (e.g., columns vs
rows) shouldn't really matter. For example, the images of a movie can be
natively represented as an array with four dimensions: time, row, column and
color.

pandas has historically supported N-dimensional panels, but deprecated them in
version 0.20 in favor of xarray data structures. There are now built-in methods
on both sides to convert between pandas and xarray, allowing for more focused
development effort. Xarray objects have a much richer model of dimensionality -
if you were using Panels:

- You need to create a new factory type for each dimensionality.
- You can't do math between NDPanels with different dimensionality.
- Each dimension in a NDPanel has a name (e.g., 'labels', 'items',
  'major_axis', etc.) but the dimension names refer to order, not their
  meaning. You can't specify an operation as to be applied along the "time"
  axis.
- You often have to manually convert collections of pandas arrays
  (Series, DataFrames, etc) to have the same number of dimensions.
  In contrast, this sort of data structure fits very naturally in an
  xarray ``Dataset``.

You can :ref:`read about switching from Panels to xarray here <panel transition>`.
pandas gets a lot of things right, but many science, engineering and complex
analytics use cases need fully multi-dimensional data structures.

How do xarray data structures differ from those found in pandas?
----------------------------------------------------------------

The main distinguishing feature of xarray's ``DataArray`` over labeled arrays in
pandas is that dimensions can have names (e.g., "time", "latitude",
"longitude"). Names are much easier to keep track of than axis numbers, and
xarray uses dimension names for indexing, aggregation and broadcasting. Not only
can you write ``x.sel(time='2000-01-01')`` and  ``x.mean(dim='time')``, but
operations like ``x - x.mean(dim='time')`` always work, no matter the order
of the "time" dimension. You never need to reshape arrays (e.g., with
``np.newaxis``) to align them for arithmetic operations in xarray.


Why don't aggregations return Python scalars?
---------------------------------------------

Xarray tries hard to be self-consistent: operations on a ``DataArray`` (resp.
``Dataset``) return another ``DataArray`` (resp. ``Dataset``) object. In
particular, operations returning scalar values (e.g. indexing or aggregations
like ``mean`` or ``sum`` applied to all axes) will also return xarray objects.

Unfortunately, this means we sometimes have to explicitly cast our results from
xarray when using them in other libraries. As an illustration, the following
code fragment

.. ipython:: python

    arr = xr.DataArray([1, 2, 3])
    pd.Series({"x": arr[0], "mean": arr.mean(), "std": arr.std()})

does not yield the pandas DataFrame we expected. We need to specify the type
conversion ourselves:

.. ipython:: python

    pd.Series({"x": arr[0], "mean": arr.mean(), "std": arr.std()}, dtype=float)

Alternatively, we could use the ``item`` method or the ``float`` constructor to
convert values one at a time

.. ipython:: python

    pd.Series({"x": arr[0].item(), "mean": float(arr.mean())})


.. _approach to metadata:

What is your approach to metadata?
----------------------------------

We are firm believers in the power of labeled data! In addition to dimensions
and coordinates, xarray supports arbitrary metadata in the form of global
(Dataset) and variable specific (DataArray) attributes (``attrs``).

Automatic interpretation of labels is powerful but also reduces flexibility.
With xarray, we draw a firm line between labels that the library understands
(``dims`` and ``coords``) and labels for users and user code (``attrs``). For
example, we do not automatically interpret and enforce units or `CF
conventions`_. (An exception is serialization to and from netCDF files.)

.. _CF conventions: http://cfconventions.org/latest.html

An implication of this choice is that we do not propagate ``attrs`` through
most operations unless explicitly flagged (some methods have a ``keep_attrs``
option, and there is a global flag, accessible with :py:func:`xarray.set_options`,
for setting this to be always True or False). Similarly, xarray does not check
for conflicts between ``attrs`` when combining arrays and datasets, unless
explicitly requested with the option ``compat='identical'``. The guiding
principle is that metadata should not be allowed to get in the way.

What other netCDF related Python libraries should I know about?
---------------------------------------------------------------

`netCDF4-python`__ provides a lower level interface for working with
netCDF and OpenDAP datasets in Python. We use netCDF4-python internally in
xarray, and have contributed a number of improvements and fixes upstream. Xarray
does not yet support all of netCDF4-python's features, such as modifying files
on-disk.

__ https://github.com/Unidata/netcdf4-python

Iris_ (supported by the UK Met office) provides similar tools for in-
memory manipulation of labeled arrays, aimed specifically at weather and
climate data needs. Indeed, the Iris :py:class:`~iris.cube.Cube` was direct
inspiration for xarray's :py:class:`~xarray.DataArray`. Xarray and Iris take very
different approaches to handling metadata: Iris strictly interprets
`CF conventions`_. Iris particularly shines at mapping, thanks to its
integration with Cartopy_.

.. _Iris: https://scitools-iris.readthedocs.io/en/stable/
.. _Cartopy: http://scitools.org.uk/cartopy/docs/latest/

`UV-CDAT`__ is another Python library that implements in-memory netCDF-like
variables and `tools for working with climate data`__.

__ http://uvcdat.llnl.gov/
__ http://drclimate.wordpress.com/2014/01/02/a-beginners-guide-to-scripting-with-uv-cdat/

We think the design decisions we have made for xarray (namely, basing it on
pandas) make it a faster and more flexible data analysis tool. That said, Iris
and CDAT have some great domain specific functionality, and xarray includes
methods for converting back and forth between xarray and these libraries. See
:py:meth:`~xarray.DataArray.to_iris` and :py:meth:`~xarray.DataArray.to_cdms2`
for more details.

What other projects leverage xarray?
------------------------------------

See section :ref:`ecosystem`.

How should I cite xarray?
-------------------------

If you are using xarray and would like to cite it in academic publication, we
would certainly appreciate it. We recommend two citations.

  1. At a minimum, we recommend citing the xarray overview journal article,
     published in the Journal of Open Research Software.

     - Hoyer, S. & Hamman, J., (2017). xarray: N-D labeled Arrays and
       Datasets in Python. Journal of Open Research Software. 5(1), p.10.
       DOI: http://doi.org/10.5334/jors.148

       Here’s an example of a BibTeX entry::

           @article{hoyer2017xarray,
             title     = {xarray: {N-D} labeled arrays and datasets in {Python}},
             author    = {Hoyer, S. and J. Hamman},
             journal   = {Journal of Open Research Software},
             volume    = {5},
             number    = {1},
             year      = {2017},
             publisher = {Ubiquity Press},
             doi       = {10.5334/jors.148},
             url       = {http://doi.org/10.5334/jors.148}
           }

  2. You may also want to cite a specific version of the xarray package. We
     provide a `Zenodo citation and DOI <https://doi.org/10.5281/zenodo.598201>`_
     for this purpose:

        .. image:: https://zenodo.org/badge/doi/10.5281/zenodo.598201.svg
           :target: https://doi.org/10.5281/zenodo.598201

       An example BibTeX entry::

           @misc{xarray_v0_8_0,
                 author = {Stephan Hoyer and Clark Fitzgerald and Joe Hamman and others},
                 title  = {xarray: v0.8.0},
                 month  = aug,
                 year   = 2016,
                 doi    = {10.5281/zenodo.59499},
                 url    = {https://doi.org/10.5281/zenodo.59499}
                }

.. _public api:

What parts of xarray are considered public API?
-----------------------------------------------

As a rule, only functions/methods documented in our :ref:`api` are considered
part of xarray's public API. Everything else (in particular, everything in
``xarray.core`` that is not also exposed in the top level ``xarray`` namespace)
is considered a private implementation detail that may change at any time.

Objects that exist to facilitate xarray's fluent interface on ``DataArray`` and
``Dataset`` objects are a special case. For convenience, we document them in
the API docs, but only their methods and the ``DataArray``/``Dataset``
methods/properties to construct them (e.g., ``.plot()``, ``.groupby()``,
``.str``) are considered public API. Constructors and other details of the
internal classes used to implemented them (i.e.,
``xarray.plot.plotting._PlotMethods``, ``xarray.core.groupby.DataArrayGroupBy``,
``xarray.core.accessor_str.StringAccessor``) are not.
################
Getting Started
################

The getting started guide aims to get you using xarray productively as quickly as possible.
It is designed as an entry point for new users, and it provided an introduction to xarray's main concepts.

.. toctree::
   :maxdepth: 2
   :hidden:

   why-xarray
   installing
   quick-overview
   faq
.. currentmodule:: xarray

.. _comput:

###########
Computation
###########


The labels associated with :py:class:`~xarray.DataArray` and
:py:class:`~xarray.Dataset` objects enables some powerful shortcuts for
computation, notably including aggregation and broadcasting by dimension
names.

Basic array math
================

Arithmetic operations with a single DataArray automatically vectorize (like
numpy) over all array values:

.. ipython:: python
    :suppress:

    import numpy as np
    import pandas as pd
    import xarray as xr

    np.random.seed(123456)

.. ipython:: python

    arr = xr.DataArray(
        np.random.RandomState(0).randn(2, 3), [("x", ["a", "b"]), ("y", [10, 20, 30])]
    )
    arr - 3
    abs(arr)

You can also use any of numpy's or scipy's many `ufunc`__ functions directly on
a DataArray:

__ http://docs.scipy.org/doc/numpy/reference/ufuncs.html

.. ipython:: python

    np.sin(arr)

Use :py:func:`~xarray.where` to conditionally switch between values:

.. ipython:: python

    xr.where(arr > 0, "positive", "negative")

Use `@` to perform matrix multiplication:

.. ipython:: python

    arr @ arr

Data arrays also implement many :py:class:`numpy.ndarray` methods:

.. ipython:: python

    arr.round(2)
    arr.T

.. _missing_values:

Missing values
==============

Xarray objects borrow the :py:meth:`~xarray.DataArray.isnull`,
:py:meth:`~xarray.DataArray.notnull`, :py:meth:`~xarray.DataArray.count`,
:py:meth:`~xarray.DataArray.dropna`, :py:meth:`~xarray.DataArray.fillna`,
:py:meth:`~xarray.DataArray.ffill`, and :py:meth:`~xarray.DataArray.bfill`
methods for working with missing data from pandas:

.. ipython:: python

    x = xr.DataArray([0, 1, np.nan, np.nan, 2], dims=["x"])
    x.isnull()
    x.notnull()
    x.count()
    x.dropna(dim="x")
    x.fillna(-1)
    x.ffill("x")
    x.bfill("x")

Like pandas, xarray uses the float value ``np.nan`` (not-a-number) to represent
missing values.

Xarray objects also have an :py:meth:`~xarray.DataArray.interpolate_na` method
for filling missing values via 1D interpolation.

.. ipython:: python

    x = xr.DataArray(
        [0, 1, np.nan, np.nan, 2],
        dims=["x"],
        coords={"xx": xr.Variable("x", [0, 1, 1.1, 1.9, 3])},
    )
    x.interpolate_na(dim="x", method="linear", use_coordinate="xx")

Note that xarray slightly diverges from the pandas ``interpolate`` syntax by
providing the ``use_coordinate`` keyword which facilitates a clear specification
of which values to use as the index in the interpolation.
Xarray also provides the ``max_gap`` keyword argument to limit the interpolation to
data gaps of length ``max_gap`` or smaller. See :py:meth:`~xarray.DataArray.interpolate_na`
for more.

Aggregation
===========

Aggregation methods have been updated to take a `dim` argument instead of
`axis`. This allows for very intuitive syntax for aggregation methods that are
applied along particular dimension(s):

.. ipython:: python

    arr.sum(dim="x")
    arr.std(["x", "y"])
    arr.min()


If you need to figure out the axis number for a dimension yourself (say,
for wrapping code designed to work with numpy arrays), you can use the
:py:meth:`~xarray.DataArray.get_axis_num` method:

.. ipython:: python

    arr.get_axis_num("y")

These operations automatically skip missing values, like in pandas:

.. ipython:: python

    xr.DataArray([1, 2, np.nan, 3]).mean()

If desired, you can disable this behavior by invoking the aggregation method
with ``skipna=False``.

.. _comput.rolling:

Rolling window operations
=========================

``DataArray`` objects include a :py:meth:`~xarray.DataArray.rolling` method. This
method supports rolling window aggregation:

.. ipython:: python

    arr = xr.DataArray(np.arange(0, 7.5, 0.5).reshape(3, 5), dims=("x", "y"))
    arr

:py:meth:`~xarray.DataArray.rolling` is applied along one dimension using the
name of the dimension as a key (e.g. ``y``) and the window size as the value
(e.g. ``3``).  We get back a ``Rolling`` object:

.. ipython:: python

    arr.rolling(y=3)

Aggregation and summary methods can be applied directly to the ``Rolling``
object:

.. ipython:: python

    r = arr.rolling(y=3)
    r.reduce(np.std)
    r.mean()

Aggregation results are assigned the coordinate at the end of each window by
default, but can be centered by passing ``center=True`` when constructing the
``Rolling`` object:

.. ipython:: python

    r = arr.rolling(y=3, center=True)
    r.mean()

As can be seen above, aggregations of windows which overlap the border of the
array produce ``nan``\s.  Setting ``min_periods`` in the call to ``rolling``
changes the minimum number of observations within the window required to have
a value when aggregating:

.. ipython:: python

    r = arr.rolling(y=3, min_periods=2)
    r.mean()
    r = arr.rolling(y=3, center=True, min_periods=2)
    r.mean()

From version 0.17, xarray supports multidimensional rolling,

.. ipython:: python

    r = arr.rolling(x=2, y=3, min_periods=2)
    r.mean()

.. tip::

   Note that rolling window aggregations are faster and use less memory when bottleneck_ is installed. This only applies to numpy-backed xarray objects with 1d-rolling.

.. _bottleneck: https://github.com/pydata/bottleneck/

We can also manually iterate through ``Rolling`` objects:

.. code:: python

    for label, arr_window in r:
        # arr_window is a view of x
        ...

.. _comput.rolling_exp:

While ``rolling`` provides a simple moving average, ``DataArray`` also supports
an exponential moving average with :py:meth:`~xarray.DataArray.rolling_exp`.
This is similar to pandas' ``ewm`` method. numbagg_ is required.

.. _numbagg: https://github.com/shoyer/numbagg

.. code:: python

    arr.rolling_exp(y=3).mean()

The ``rolling_exp`` method takes a ``window_type`` kwarg, which can be ``'alpha'``,
``'com'`` (for ``center-of-mass``), ``'span'``, and ``'halflife'``. The default is
``span``.

Finally, the rolling object has a ``construct`` method which returns a
view of the original ``DataArray`` with the windowed dimension in
the last position.
You can use this for more advanced rolling operations such as strided rolling,
windowed rolling, convolution, short-time FFT etc.

.. ipython:: python

    # rolling with 2-point stride
    rolling_da = r.construct(x="x_win", y="y_win", stride=2)
    rolling_da
    rolling_da.mean(["x_win", "y_win"], skipna=False)

Because the ``DataArray`` given by ``r.construct('window_dim')`` is a view
of the original array, it is memory efficient.
You can also use ``construct`` to compute a weighted rolling sum:

.. ipython:: python

    weight = xr.DataArray([0.25, 0.5, 0.25], dims=["window"])
    arr.rolling(y=3).construct(y="window").dot(weight)

.. note::
  numpy's Nan-aggregation functions such as ``nansum`` copy the original array.
  In xarray, we internally use these functions in our aggregation methods
  (such as ``.sum()``) if ``skipna`` argument is not specified or set to True.
  This means ``rolling_da.mean('window_dim')`` is memory inefficient.
  To avoid this, use ``skipna=False`` as the above example.


.. _comput.weighted:

Weighted array reductions
=========================

:py:class:`DataArray` and :py:class:`Dataset` objects include :py:meth:`DataArray.weighted`
and :py:meth:`Dataset.weighted` array reduction methods. They currently
support weighted ``sum``, ``mean``, ``std`` and ``var``.

.. ipython:: python

    coords = dict(month=("month", [1, 2, 3]))

    prec = xr.DataArray([1.1, 1.0, 0.9], dims=("month",), coords=coords)
    weights = xr.DataArray([31, 28, 31], dims=("month",), coords=coords)

Create a weighted object:

.. ipython:: python

    weighted_prec = prec.weighted(weights)
    weighted_prec

Calculate the weighted sum:

.. ipython:: python

    weighted_prec.sum()

Calculate the weighted mean:

.. ipython:: python

    weighted_prec.mean(dim="month")

The weighted sum corresponds to:

.. ipython:: python

    weighted_sum = (prec * weights).sum()
    weighted_sum

the weighted mean to:

.. ipython:: python

    weighted_mean = weighted_sum / weights.sum()
    weighted_mean

the weighted variance to:

.. ipython:: python

    weighted_var = weighted_prec.sum_of_squares() / weights.sum()
    weighted_var

and the weighted standard deviation to:

.. ipython:: python

    weighted_std = np.sqrt(weighted_var)
    weighted_std

However, the functions also take missing values in the data into account:

.. ipython:: python

    data = xr.DataArray([np.NaN, 2, 4])
    weights = xr.DataArray([8, 1, 1])

    data.weighted(weights).mean()

Using ``(data * weights).sum() / weights.sum()`` would (incorrectly) result
in 0.6.


If the weights add up to to 0, ``sum`` returns 0:

.. ipython:: python

    data = xr.DataArray([1.0, 1.0])
    weights = xr.DataArray([-1.0, 1.0])

    data.weighted(weights).sum()

and ``mean``, ``std`` and ``var`` return ``NaN``:

.. ipython:: python

    data.weighted(weights).mean()


.. note::
  ``weights`` must be a :py:class:`DataArray` and cannot contain missing values.
  Missing values can be replaced manually by ``weights.fillna(0)``.

.. _comput.coarsen:

Coarsen large arrays
====================

:py:class:`DataArray` and :py:class:`Dataset` objects include a
:py:meth:`~xarray.DataArray.coarsen` and :py:meth:`~xarray.Dataset.coarsen`
methods. This supports the block aggregation along multiple dimensions,

.. ipython:: python

    x = np.linspace(0, 10, 300)
    t = pd.date_range("1999-12-15", periods=364)
    da = xr.DataArray(
        np.sin(x) * np.cos(np.linspace(0, 1, 364)[:, np.newaxis]),
        dims=["time", "x"],
        coords={"time": t, "x": x},
    )
    da

In order to take a block mean for every 7 days along ``time`` dimension and
every 2 points along ``x`` dimension,

.. ipython:: python

    da.coarsen(time=7, x=2).mean()

:py:meth:`~xarray.DataArray.coarsen` raises an ``ValueError`` if the data
length is not a multiple of the corresponding window size.
You can choose ``boundary='trim'`` or ``boundary='pad'`` options for trimming
the excess entries or padding ``nan`` to insufficient entries,

.. ipython:: python

    da.coarsen(time=30, x=2, boundary="trim").mean()

If you want to apply a specific function to coordinate, you can pass the
function or method name to ``coord_func`` option,

.. ipython:: python

    da.coarsen(time=7, x=2, coord_func={"time": "min"}).mean()


.. _compute.using_coordinates:

Computation using Coordinates
=============================

Xarray objects have some handy methods for the computation with their
coordinates. :py:meth:`~xarray.DataArray.differentiate` computes derivatives by
central finite differences using their coordinates,

.. ipython:: python

    a = xr.DataArray([0, 1, 2, 3], dims=["x"], coords=[[0.1, 0.11, 0.2, 0.3]])
    a
    a.differentiate("x")

This method can be used also for multidimensional arrays,

.. ipython:: python

    a = xr.DataArray(
        np.arange(8).reshape(4, 2), dims=["x", "y"], coords={"x": [0.1, 0.11, 0.2, 0.3]}
    )
    a.differentiate("x")

:py:meth:`~xarray.DataArray.integrate` computes integration based on
trapezoidal rule using their coordinates,

.. ipython:: python

    a.integrate("x")

.. note::
    These methods are limited to simple cartesian geometry. Differentiation
    and integration along multidimensional coordinate are not supported.


.. _compute.polyfit:

Fitting polynomials
===================

Xarray objects provide an interface for performing linear or polynomial regressions
using the least-squares method. :py:meth:`~xarray.DataArray.polyfit` computes the
best fitting coefficients along a given dimension and for a given order,

.. ipython:: python

    x = xr.DataArray(np.arange(10), dims=["x"], name="x")
    a = xr.DataArray(3 + 4 * x, dims=["x"], coords={"x": x})
    out = a.polyfit(dim="x", deg=1, full=True)
    out

The method outputs a dataset containing the coefficients (and more if `full=True`).
The inverse operation is done with :py:meth:`~xarray.polyval`,

.. ipython:: python

    xr.polyval(coord=x, coeffs=out.polyfit_coefficients)

.. note::
    These methods replicate the behaviour of :py:func:`numpy.polyfit` and :py:func:`numpy.polyval`.


.. _compute.curvefit:

Fitting arbitrary functions
===========================

Xarray objects also provide an interface for fitting more complex functions using
:py:func:`scipy.optimize.curve_fit`. :py:meth:`~xarray.DataArray.curvefit` accepts
user-defined functions and can fit along multiple coordinates.

For example, we can fit a relationship between two ``DataArray`` objects, maintaining
a unique fit at each spatial coordinate but aggregating over the time dimension:

.. ipython:: python

    def exponential(x, a, xc):
        return np.exp((x - xc) / a)


    x = np.arange(-5, 5, 0.1)
    t = np.arange(-5, 5, 0.1)
    X, T = np.meshgrid(x, t)
    Z1 = np.random.uniform(low=-5, high=5, size=X.shape)
    Z2 = exponential(Z1, 3, X)
    Z3 = exponential(Z1, 1, -X)

    ds = xr.Dataset(
        data_vars=dict(
            var1=(["t", "x"], Z1), var2=(["t", "x"], Z2), var3=(["t", "x"], Z3)
        ),
        coords={"t": t, "x": x},
    )
    ds[["var2", "var3"]].curvefit(
        coords=ds.var1,
        func=exponential,
        reduce_dims="t",
        bounds={"a": (0.5, 5), "xc": (-5, 5)},
    )

We can also fit multi-dimensional functions, and even use a wrapper function to
simultaneously fit a summation of several functions, such as this field containing
two gaussian peaks:

.. ipython:: python

    def gaussian_2d(coords, a, xc, yc, xalpha, yalpha):
        x, y = coords
        z = a * np.exp(
            -np.square(x - xc) / 2 / np.square(xalpha)
            - np.square(y - yc) / 2 / np.square(yalpha)
        )
        return z


    def multi_peak(coords, *args):
        z = np.zeros(coords[0].shape)
        for i in range(len(args) // 5):
            z += gaussian_2d(coords, *args[i * 5 : i * 5 + 5])
        return z


    x = np.arange(-5, 5, 0.1)
    y = np.arange(-5, 5, 0.1)
    X, Y = np.meshgrid(x, y)

    n_peaks = 2
    names = ["a", "xc", "yc", "xalpha", "yalpha"]
    names = [f"{name}{i}" for i in range(n_peaks) for name in names]
    Z = gaussian_2d((X, Y), 3, 1, 1, 2, 1) + gaussian_2d((X, Y), 2, -1, -2, 1, 1)
    Z += np.random.normal(scale=0.1, size=Z.shape)

    da = xr.DataArray(Z, dims=["y", "x"], coords={"y": y, "x": x})
    da.curvefit(
        coords=["x", "y"],
        func=multi_peak,
        param_names=names,
        kwargs={"maxfev": 10000},
    )

.. note::
    This method replicates the behavior of :py:func:`scipy.optimize.curve_fit`.


.. _compute.broadcasting:

Broadcasting by dimension name
==============================

``DataArray`` objects automatically align themselves ("broadcasting" in
the numpy parlance) by dimension name instead of axis order. With xarray, you
do not need to transpose arrays or insert dimensions of length 1 to get array
operations to work, as commonly done in numpy with :py:func:`numpy.reshape` or
:py:data:`numpy.newaxis`.

This is best illustrated by a few examples. Consider two one-dimensional
arrays with different sizes aligned along different dimensions:

.. ipython:: python

    a = xr.DataArray([1, 2], [("x", ["a", "b"])])
    a
    b = xr.DataArray([-1, -2, -3], [("y", [10, 20, 30])])
    b

With xarray, we can apply binary mathematical operations to these arrays, and
their dimensions are expanded automatically:

.. ipython:: python

    a * b

Moreover, dimensions are always reordered to the order in which they first
appeared:

.. ipython:: python

    c = xr.DataArray(np.arange(6).reshape(3, 2), [b["y"], a["x"]])
    c
    a + c

This means, for example, that you always subtract an array from its transpose:

.. ipython:: python

    c - c.T

You can explicitly broadcast xarray data structures by using the
:py:func:`~xarray.broadcast` function:

.. ipython:: python

    a2, b2 = xr.broadcast(a, b)
    a2
    b2

.. _math automatic alignment:

Automatic alignment
===================

Xarray enforces alignment between *index* :ref:`coordinates` (that is,
coordinates with the same name as a dimension, marked by ``*``) on objects used
in binary operations.

Similarly to pandas, this alignment is automatic for arithmetic on binary
operations. The default result of a binary operation is by the *intersection*
(not the union) of coordinate labels:

.. ipython:: python

    arr = xr.DataArray(np.arange(3), [("x", range(3))])
    arr + arr[:-1]

If coordinate values for a dimension are missing on either argument, all
matching dimensions must have the same size:

.. ipython::
    :verbatim:

    In [1]: arr + xr.DataArray([1, 2], dims="x")
    ValueError: arguments without labels along dimension 'x' cannot be aligned because they have different dimension size(s) {2} than the size of the aligned dimension labels: 3


However, one can explicitly change this default automatic alignment type ("inner")
via :py:func:`~xarray.set_options()` in context manager:

.. ipython:: python

    with xr.set_options(arithmetic_join="outer"):
        arr + arr[:1]
    arr + arr[:1]

Before loops or performance critical code, it's a good idea to align arrays
explicitly (e.g., by putting them in the same Dataset or using
:py:func:`~xarray.align`) to avoid the overhead of repeated alignment with each
operation. See :ref:`align and reindex` for more details.

.. note::

    There is no automatic alignment between arguments when performing in-place
    arithmetic operations such as ``+=``. You will need to use
    :ref:`manual alignment<align and reindex>`. This ensures in-place
    arithmetic never needs to modify data types.

.. _coordinates math:

Coordinates
===========

Although index coordinates are aligned, other coordinates are not, and if their
values conflict, they will be dropped. This is necessary, for example, because
indexing turns 1D coordinates into scalar coordinates:

.. ipython:: python

    arr[0]
    arr[1]
    # notice that the scalar coordinate 'x' is silently dropped
    arr[1] - arr[0]

Still, xarray will persist other coordinates in arithmetic, as long as there
are no conflicting values:

.. ipython:: python

    # only one argument has the 'x' coordinate
    arr[0] + 1
    # both arguments have the same 'x' coordinate
    arr[0] - arr[0]

Math with datasets
==================

Datasets support arithmetic operations by automatically looping over all data
variables:

.. ipython:: python

    ds = xr.Dataset(
        {
            "x_and_y": (("x", "y"), np.random.randn(3, 5)),
            "x_only": ("x", np.random.randn(3)),
        },
        coords=arr.coords,
    )
    ds > 0

Datasets support most of the same methods found on data arrays:

.. ipython:: python

    ds.mean(dim="x")
    abs(ds)

Datasets also support NumPy ufuncs (requires NumPy v1.13 or newer), or
alternatively you can use :py:meth:`~xarray.Dataset.map` to map a function
to each variable in a dataset:

.. ipython:: python

    np.sin(ds)
    ds.map(np.sin)

Datasets also use looping over variables for *broadcasting* in binary
arithmetic. You can do arithmetic between any ``DataArray`` and a dataset:

.. ipython:: python

    ds + arr

Arithmetic between two datasets matches data variables of the same name:

.. ipython:: python

    ds2 = xr.Dataset({"x_and_y": 0, "x_only": 100})
    ds - ds2

Similarly to index based alignment, the result has the intersection of all
matching data variables.

.. _comput.wrapping-custom:

Wrapping custom computation
===========================

It doesn't always make sense to do computation directly with xarray objects:

  - In the inner loop of performance limited code, using xarray can add
    considerable overhead compared to using NumPy or native Python types.
    This is particularly true when working with scalars or small arrays (less
    than ~1e6 elements). Keeping track of labels and ensuring their consistency
    adds overhead, and xarray's core itself is not especially fast, because it's
    written in Python rather than a compiled language like C. Also, xarray's
    high level label-based APIs removes low-level control over how operations
    are implemented.
  - Even if speed doesn't matter, it can be important to wrap existing code, or
    to support alternative interfaces that don't use xarray objects.

For these reasons, it is often well-advised to write low-level routines that
work with NumPy arrays, and to wrap these routines to work with xarray objects.
However, adding support for labels on both :py:class:`~xarray.Dataset` and
:py:class:`~xarray.DataArray` can be a bit of a chore.

To make this easier, xarray supplies the :py:func:`~xarray.apply_ufunc` helper
function, designed for wrapping functions that support broadcasting and
vectorization on unlabeled arrays in the style of a NumPy
`universal function <https://docs.scipy.org/doc/numpy-1.13.0/reference/ufuncs.html>`_ ("ufunc" for short).
``apply_ufunc`` takes care of everything needed for an idiomatic xarray wrapper,
including alignment, broadcasting, looping over ``Dataset`` variables (if
needed), and merging of coordinates. In fact, many internal xarray
functions/methods are written using ``apply_ufunc``.

Simple functions that act independently on each value should work without
any additional arguments:

.. ipython:: python

    squared_error = lambda x, y: (x - y) ** 2
    arr1 = xr.DataArray([0, 1, 2, 3], dims="x")
    xr.apply_ufunc(squared_error, arr1, 1)

For using more complex operations that consider some array values collectively,
it's important to understand the idea of "core dimensions" from NumPy's
`generalized ufuncs <http://docs.scipy.org/doc/numpy/reference/c-api/generalized-ufuncs.html>`_. Core dimensions are defined as dimensions
that should *not* be broadcast over. Usually, they correspond to the fundamental
dimensions over which an operation is defined, e.g., the summed axis in
``np.sum``. A good clue that core dimensions are needed is the presence of an
``axis`` argument on the corresponding NumPy function.

With ``apply_ufunc``, core dimensions are recognized by name, and then moved to
the last dimension of any input arguments before applying the given function.
This means that for functions that accept an ``axis`` argument, you usually need
to set ``axis=-1``. As an example, here is how we would wrap
:py:func:`numpy.linalg.norm` to calculate the vector norm:

.. code-block:: python

    def vector_norm(x, dim, ord=None):
        return xr.apply_ufunc(
            np.linalg.norm, x, input_core_dims=[[dim]], kwargs={"ord": ord, "axis": -1}
        )

.. ipython:: python
    :suppress:

    def vector_norm(x, dim, ord=None):
        return xr.apply_ufunc(
            np.linalg.norm, x, input_core_dims=[[dim]], kwargs={"ord": ord, "axis": -1}
        )

.. ipython:: python

    vector_norm(arr1, dim="x")

Because ``apply_ufunc`` follows a standard convention for ufuncs, it plays
nicely with tools for building vectorized functions, like
:py:func:`numpy.broadcast_arrays` and :py:class:`numpy.vectorize`. For high performance
needs, consider using Numba's :doc:`vectorize and guvectorize <numba:user/vectorize>`.

In addition to wrapping functions, ``apply_ufunc`` can automatically parallelize
many functions when using dask by setting ``dask='parallelized'``. See
:ref:`dask.automatic-parallelization` for details.

:py:func:`~xarray.apply_ufunc` also supports some advanced options for
controlling alignment of variables and the form of the result. See the
docstring for full details and more examples.
.. currentmodule:: xarray
.. _terminology:

Terminology
===========

*Xarray terminology differs slightly from CF, mathematical conventions, and
pandas; so we've put together a glossary of its terms. Here,* ``arr`` *
refers to an xarray* :py:class:`DataArray` *in the examples. For more
complete examples, please consult the relevant documentation.*

.. glossary::

    DataArray
        A multi-dimensional array with labeled or named
        dimensions. ``DataArray`` objects add metadata such as dimension names,
        coordinates, and attributes (defined below) to underlying "unlabeled"
        data structures such as numpy and Dask arrays. If its optional ``name``
        property is set, it is a *named DataArray*.

    Dataset
        A dict-like collection of ``DataArray`` objects with aligned
        dimensions. Thus, most operations that can be performed on the
        dimensions of a single ``DataArray`` can be performed on a
        dataset. Datasets have data variables (see **Variable** below),
        dimensions, coordinates, and attributes.

    Variable
        A `NetCDF-like variable
        <https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_data_set_components.html#variables>`_
        consisting of dimensions, data, and attributes which describe a single
        array. The main functional difference between variables and numpy arrays
        is that numerical operations on variables implement array broadcasting
        by dimension name. Each ``DataArray`` has an underlying variable that
        can be accessed via ``arr.variable``. However, a variable is not fully
        described outside of either a ``Dataset`` or a ``DataArray``.

        .. note::

            The :py:class:`Variable` class is low-level interface and can
            typically be ignored. However, the word "variable" appears often
            enough in the code and documentation that is useful to understand.

    Dimension
        In mathematics, the *dimension* of data is loosely the number of degrees
        of freedom for it. A *dimension axis* is a set of all points in which
        all but one of these degrees of freedom is fixed. We can think of each
        dimension axis as having a name, for example the "x dimension".  In
        xarray, a ``DataArray`` object's *dimensions* are its named dimension
        axes, and the name of the ``i``-th dimension is ``arr.dims[i]``. If an
        array is created without dimension names, the default dimension names are
        ``dim_0``, ``dim_1``, and so forth.

    Coordinate
        An array that labels a dimension or set of dimensions of another
        ``DataArray``. In the usual one-dimensional case, the coordinate array's
        values can loosely be thought of as tick labels along a dimension. There
        are two types of coordinate arrays: *dimension coordinates* and
        *non-dimension coordinates* (see below). A coordinate named ``x`` can be
        retrieved from ``arr.coords[x]``. A ``DataArray`` can have more
        coordinates than dimensions because a single dimension can be labeled by
        multiple coordinate arrays. However, only one coordinate array can be a
        assigned as a particular dimension's dimension coordinate array. As a
        consequence, ``len(arr.dims) <= len(arr.coords)`` in general.

    Dimension coordinate
        A one-dimensional coordinate array assigned to ``arr`` with both a name
        and dimension name in ``arr.dims``. Dimension coordinates are used for
        label-based indexing and alignment, like the index found on a
        :py:class:`pandas.DataFrame` or :py:class:`pandas.Series`. In fact,
        dimension coordinates use :py:class:`pandas.Index` objects under the
        hood for efficient computation. Dimension coordinates are marked by
        ``*`` when printing a ``DataArray`` or ``Dataset``.

    Non-dimension coordinate
        A coordinate array assigned to ``arr`` with a name in ``arr.coords`` but
        *not* in ``arr.dims``. These coordinates arrays can be one-dimensional
        or multidimensional, and they are useful for auxiliary labeling. As an
        example, multidimensional coordinates are often used in geoscience
        datasets when :doc:`the data's physical coordinates (such as latitude
        and longitude) differ from their logical coordinates
        <../examples/multidimensional-coords>`. However, non-dimension coordinates
        are not indexed, and any operation on non-dimension coordinates that
        leverages indexing will fail. Printing ``arr.coords`` will print all of
        ``arr``'s coordinate names, with the corresponding dimension(s) in
        parentheses. For example, ``coord_name (dim_name) 1 2 3 ...``.

    Index
        An *index* is a data structure optimized for efficient selecting and
        slicing of an associated array. Xarray creates indexes for dimension
        coordinates so that operations along dimensions are fast, while
        non-dimension coordinates are not indexed. Under the hood, indexes are
        implemented as :py:class:`pandas.Index` objects. The index associated
        with dimension name ``x`` can be retrieved by ``arr.indexes[x]``. By
        construction, ``len(arr.dims) == len(arr.indexes)``

    name
        The names of dimensions, coordinates, DataArray objects and data
        variables can be anything as long as they are :term:`hashable`. However,
        it is preferred to use :py:class:`str` typed names.

    scalar
        By definition, a scalar is not an :term:`array` and when converted to
        one, it has 0 dimensions. That means that, e.g., :py:class:`int`,
        :py:class:`float`, and :py:class:`str` objects are "scalar" while
        :py:class:`list` or :py:class:`tuple` are not.

    duck array
        `Duck arrays`__ are array implementations that behave
        like numpy arrays. They have to define the ``shape``, ``dtype`` and
        ``ndim`` properties. For integration with ``xarray``, the ``__array__``,
        ``__array_ufunc__`` and ``__array_function__`` protocols are also required.

        __ https://numpy.org/neps/nep-0022-ndarray-duck-typing-overview.html
.. _interp:

Interpolating data
==================

.. ipython:: python
    :suppress:

    import numpy as np
    import pandas as pd
    import xarray as xr

    np.random.seed(123456)

Xarray offers flexible interpolation routines, which have a similar interface
to our :ref:`indexing <indexing>`.

.. note::

  ``interp`` requires `scipy` installed.


Scalar and 1-dimensional interpolation
--------------------------------------

Interpolating a :py:class:`~xarray.DataArray` works mostly like labeled
indexing of a :py:class:`~xarray.DataArray`,

.. ipython:: python

    da = xr.DataArray(
        np.sin(0.3 * np.arange(12).reshape(4, 3)),
        [("time", np.arange(4)), ("space", [0.1, 0.2, 0.3])],
    )
    # label lookup
    da.sel(time=3)

    # interpolation
    da.interp(time=2.5)


Similar to the indexing, :py:meth:`~xarray.DataArray.interp` also accepts an
array-like, which gives the interpolated result as an array.

.. ipython:: python

    # label lookup
    da.sel(time=[2, 3])

    # interpolation
    da.interp(time=[2.5, 3.5])

To interpolate data with a :py:doc:`numpy.datetime64 <reference/arrays.datetime>` coordinate you can pass a string.

.. ipython:: python

    da_dt64 = xr.DataArray(
        [1, 3], [("time", pd.date_range("1/1/2000", "1/3/2000", periods=2))]
    )
    da_dt64.interp(time="2000-01-02")

The interpolated data can be merged into the original :py:class:`~xarray.DataArray`
by specifying the time periods required.

.. ipython:: python

    da_dt64.interp(time=pd.date_range("1/1/2000", "1/3/2000", periods=3))

Interpolation of data indexed by a :py:class:`~xarray.CFTimeIndex` is also
allowed.  See :ref:`CFTimeIndex` for examples.

.. note::

  Currently, our interpolation only works for regular grids.
  Therefore, similarly to :py:meth:`~xarray.DataArray.sel`,
  only 1D coordinates along a dimension can be used as the
  original coordinate to be interpolated.


Multi-dimensional Interpolation
-------------------------------

Like :py:meth:`~xarray.DataArray.sel`, :py:meth:`~xarray.DataArray.interp`
accepts multiple coordinates. In this case, multidimensional interpolation
is carried out.

.. ipython:: python

    # label lookup
    da.sel(time=2, space=0.1)

    # interpolation
    da.interp(time=2.5, space=0.15)

Array-like coordinates are also accepted:

.. ipython:: python

    # label lookup
    da.sel(time=[2, 3], space=[0.1, 0.2])

    # interpolation
    da.interp(time=[1.5, 2.5], space=[0.15, 0.25])


:py:meth:`~xarray.DataArray.interp_like` method is a useful shortcut. This
method interpolates an xarray object onto the coordinates of another xarray
object. For example, if we want to compute the difference between
two :py:class:`~xarray.DataArray` s (``da`` and ``other``) staying on slightly
different coordinates,

.. ipython:: python

    other = xr.DataArray(
        np.sin(0.4 * np.arange(9).reshape(3, 3)),
        [("time", [0.9, 1.9, 2.9]), ("space", [0.15, 0.25, 0.35])],
    )

it might be a good idea to first interpolate ``da`` so that it will stay on the
same coordinates of ``other``, and then subtract it.
:py:meth:`~xarray.DataArray.interp_like` can be used for such a case,

.. ipython:: python

    # interpolate da along other's coordinates
    interpolated = da.interp_like(other)
    interpolated

It is now possible to safely compute the difference ``other - interpolated``.


Interpolation methods
---------------------

We use :py:class:`scipy.interpolate.interp1d` for 1-dimensional interpolation and
:py:func:`scipy.interpolate.interpn` for multi-dimensional interpolation.

The interpolation method can be specified by the optional ``method`` argument.

.. ipython:: python

    da = xr.DataArray(
        np.sin(np.linspace(0, 2 * np.pi, 10)),
        dims="x",
        coords={"x": np.linspace(0, 1, 10)},
    )

    da.plot.line("o", label="original")
    da.interp(x=np.linspace(0, 1, 100)).plot.line(label="linear (default)")
    da.interp(x=np.linspace(0, 1, 100), method="cubic").plot.line(label="cubic")
    @savefig interpolation_sample1.png width=4in
    plt.legend()

Additional keyword arguments can be passed to scipy's functions.

.. ipython:: python

    # fill 0 for the outside of the original coordinates.
    da.interp(x=np.linspace(-0.5, 1.5, 10), kwargs={"fill_value": 0.0})
    # 1-dimensional extrapolation
    da.interp(x=np.linspace(-0.5, 1.5, 10), kwargs={"fill_value": "extrapolate"})
    # multi-dimensional extrapolation
    da = xr.DataArray(
        np.sin(0.3 * np.arange(12).reshape(4, 3)),
        [("time", np.arange(4)), ("space", [0.1, 0.2, 0.3])],
    )

    da.interp(time=4, space=np.linspace(-0.1, 0.5, 10), kwargs={"fill_value": None})


Advanced Interpolation
----------------------

:py:meth:`~xarray.DataArray.interp` accepts :py:class:`~xarray.DataArray`
as similar to :py:meth:`~xarray.DataArray.sel`, which enables us more advanced interpolation.
Based on the dimension of the new coordinate passed to :py:meth:`~xarray.DataArray.interp`, the dimension of the result are determined.

For example, if you want to interpolate a two dimensional array along a particular dimension, as illustrated below,
you can pass two 1-dimensional :py:class:`~xarray.DataArray` s with
a common dimension as new coordinate.

.. image:: ../_static/advanced_selection_interpolation.svg
    :height: 200px
    :width: 400 px
    :alt: advanced indexing and interpolation
    :align: center

For example:

.. ipython:: python

    da = xr.DataArray(
        np.sin(0.3 * np.arange(20).reshape(5, 4)),
        [("x", np.arange(5)), ("y", [0.1, 0.2, 0.3, 0.4])],
    )
    # advanced indexing
    x = xr.DataArray([0, 2, 4], dims="z")
    y = xr.DataArray([0.1, 0.2, 0.3], dims="z")
    da.sel(x=x, y=y)

    # advanced interpolation
    x = xr.DataArray([0.5, 1.5, 2.5], dims="z")
    y = xr.DataArray([0.15, 0.25, 0.35], dims="z")
    da.interp(x=x, y=y)

where values on the original coordinates
``(x, y) = ((0.5, 0.15), (1.5, 0.25), (2.5, 0.35))`` are obtained by the
2-dimensional interpolation and mapped along a new dimension ``z``.

If you want to add a coordinate to the new dimension ``z``, you can supply
:py:class:`~xarray.DataArray` s with a coordinate,

.. ipython:: python

    x = xr.DataArray([0.5, 1.5, 2.5], dims="z", coords={"z": ["a", "b", "c"]})
    y = xr.DataArray([0.15, 0.25, 0.35], dims="z", coords={"z": ["a", "b", "c"]})
    da.interp(x=x, y=y)

For the details of the advanced indexing,
see :ref:`more advanced indexing <more_advanced_indexing>`.


Interpolating arrays with NaN
-----------------------------

Our :py:meth:`~xarray.DataArray.interp` works with arrays with NaN
the same way that
`scipy.interpolate.interp1d <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html>`_ and
`scipy.interpolate.interpn <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interpn.html>`_ do.
``linear`` and ``nearest`` methods return arrays including NaN,
while other methods such as ``cubic`` or ``quadratic`` return all NaN arrays.

.. ipython:: python

    da = xr.DataArray([0, 2, np.nan, 3, 3.25], dims="x", coords={"x": range(5)})
    da.interp(x=[0.5, 1.5, 2.5])
    da.interp(x=[0.5, 1.5, 2.5], method="cubic")

To avoid this, you can drop NaN by :py:meth:`~xarray.DataArray.dropna`, and
then make the interpolation

.. ipython:: python

    dropped = da.dropna("x")
    dropped
    dropped.interp(x=[0.5, 1.5, 2.5], method="cubic")

If NaNs are distributed randomly in your multidimensional array,
dropping all the columns containing more than one NaNs by
:py:meth:`~xarray.DataArray.dropna` may lose a significant amount of information.
In such a case, you can fill NaN by :py:meth:`~xarray.DataArray.interpolate_na`,
which is similar to :py:meth:`pandas.Series.interpolate`.

.. ipython:: python

    filled = da.interpolate_na(dim="x")
    filled

This fills NaN by interpolating along the specified dimension.
After filling NaNs, you can interpolate:

.. ipython:: python

    filled.interp(x=[0.5, 1.5, 2.5], method="cubic")

For the details of :py:meth:`~xarray.DataArray.interpolate_na`,
see :ref:`Missing values <missing_values>`.


Example
-------

Let's see how :py:meth:`~xarray.DataArray.interp` works on real data.

.. ipython:: python

    # Raw data
    ds = xr.tutorial.open_dataset("air_temperature").isel(time=0)
    fig, axes = plt.subplots(ncols=2, figsize=(10, 4))
    ds.air.plot(ax=axes[0])
    axes[0].set_title("Raw data")

    # Interpolated data
    new_lon = np.linspace(ds.lon[0], ds.lon[-1], ds.dims["lon"] * 4)
    new_lat = np.linspace(ds.lat[0], ds.lat[-1], ds.dims["lat"] * 4)
    dsi = ds.interp(lat=new_lat, lon=new_lon)
    dsi.air.plot(ax=axes[1])
    @savefig interpolation_sample3.png width=8in
    axes[1].set_title("Interpolated data")

Our advanced interpolation can be used to remap the data to the new coordinate.
Consider the new coordinates x and z on the two dimensional plane.
The remapping can be done as follows

.. ipython:: python

    # new coordinate
    x = np.linspace(240, 300, 100)
    z = np.linspace(20, 70, 100)
    # relation between new and original coordinates
    lat = xr.DataArray(z, dims=["z"], coords={"z": z})
    lon = xr.DataArray(
        (x[:, np.newaxis] - 270) / np.cos(z * np.pi / 180) + 270,
        dims=["x", "z"],
        coords={"x": x, "z": z},
    )

    fig, axes = plt.subplots(ncols=2, figsize=(10, 4))
    ds.air.plot(ax=axes[0])
    # draw the new coordinate on the original coordinates.
    for idx in [0, 33, 66, 99]:
        axes[0].plot(lon.isel(x=idx), lat, "--k")
    for idx in [0, 33, 66, 99]:
        axes[0].plot(*xr.broadcast(lon.isel(z=idx), lat.isel(z=idx)), "--k")
    axes[0].set_title("Raw data")

    dsi = ds.interp(lon=lon, lat=lat)
    dsi.air.plot(ax=axes[1])
    @savefig interpolation_sample4.png width=8in
    axes[1].set_title("Remapped data")
.. currentmodule:: xarray

Working with numpy-like arrays
==============================

.. warning::

   This feature should be considered experimental. Please report any bug you may find on
   xarray’s github repository.

NumPy-like arrays (:term:`duck array`) extend the :py:class:`numpy.ndarray` with
additional features, like propagating physical units or a different layout in memory.

:py:class:`DataArray` and :py:class:`Dataset` objects can wrap these duck arrays, as
long as they satisfy certain conditions (see :ref:`internals.duck_arrays`).

.. note::

    For ``dask`` support see :ref:`dask`.


Missing features
----------------
Most of the API does support :term:`duck array` objects, but there are a few areas where
the code will still cast to ``numpy`` arrays:

- dimension coordinates, and thus all indexing operations:

  * :py:meth:`Dataset.sel` and :py:meth:`DataArray.sel`
  * :py:meth:`Dataset.loc` and :py:meth:`DataArray.loc`
  * :py:meth:`Dataset.drop_sel` and :py:meth:`DataArray.drop_sel`
  * :py:meth:`Dataset.reindex`, :py:meth:`Dataset.reindex_like`,
    :py:meth:`DataArray.reindex` and :py:meth:`DataArray.reindex_like`: duck arrays in
    data variables and non-dimension coordinates won't be casted

- functions and methods that depend on external libraries or features of ``numpy`` not
  covered by ``__array_function__`` / ``__array_ufunc__``:

  * :py:meth:`Dataset.ffill` and :py:meth:`DataArray.ffill` (uses ``bottleneck``)
  * :py:meth:`Dataset.bfill` and :py:meth:`DataArray.bfill` (uses ``bottleneck``)
  * :py:meth:`Dataset.interp`, :py:meth:`Dataset.interp_like`,
    :py:meth:`DataArray.interp` and :py:meth:`DataArray.interp_like` (uses ``scipy``):
    duck arrays in data variables and non-dimension coordinates will be casted in
    addition to not supporting duck arrays in dimension coordinates
  * :py:meth:`Dataset.rolling` and :py:meth:`DataArray.rolling` (requires ``numpy>=1.20``)
  * :py:meth:`Dataset.rolling_exp` and :py:meth:`DataArray.rolling_exp` (uses
    ``numbagg``)
  * :py:meth:`Dataset.interpolate_na` and :py:meth:`DataArray.interpolate_na` (uses
    :py:class:`numpy.vectorize`)
  * :py:func:`apply_ufunc` with ``vectorize=True`` (uses :py:class:`numpy.vectorize`)

- incompatibilities between different :term:`duck array` libraries:

  * :py:meth:`Dataset.chunk` and :py:meth:`DataArray.chunk`: this fails if the data was
    not already chunked and the :term:`duck array` (e.g. a ``pint`` quantity) should
    wrap the new ``dask`` array; changing the chunk sizes works.


Extensions using duck arrays
----------------------------
Here's a list of libraries extending ``xarray`` to make working with wrapped duck arrays
easier:

- `pint-xarray <https://github.com/xarray-contrib/pint-xarray>`_
.. currentmodule:: xarray

.. _dask:

Parallel computing with Dask
============================

Xarray integrates with `Dask <http://dask.pydata.org/>`__ to support parallel
computations and streaming computation on datasets that don't fit into memory.
Currently, Dask is an entirely optional feature for xarray. However, the
benefits of using Dask are sufficiently strong that Dask may become a required
dependency in a future version of xarray.

For a full example of how to use xarray's Dask integration, read the
`blog post introducing xarray and Dask`_. More up-to-date examples
may be found at the `Pangeo project's gallery <http://gallery.pangeo.io/>`_
and at the `Dask examples website <https://examples.dask.org/xarray.html>`_.

.. _blog post introducing xarray and Dask: http://stephanhoyer.com/2015/06/11/xray-dask-out-of-core-labeled-arrays/

What is a Dask array?
---------------------

.. image:: ../_static/dask_array.png
   :width: 40 %
   :align: right
   :alt: A Dask array

Dask divides arrays into many small pieces, called *chunks*, each of which is
presumed to be small enough to fit into memory.

Unlike NumPy, which has eager evaluation, operations on Dask arrays are lazy.
Operations queue up a series of tasks mapped over blocks, and no computation is
performed until you actually ask values to be computed (e.g., to print results
to your screen or write to disk). At that point, data is loaded into memory
and computation proceeds in a streaming fashion, block-by-block.

The actual computation is controlled by a multi-processing or thread pool,
which allows Dask to take full advantage of multiple processors available on
most modern computers.

For more details on Dask, read `its documentation <http://dask.pydata.org/>`__.
Note that xarray only makes use of ``dask.array`` and ``dask.delayed``.

.. _dask.io:

Reading and writing data
------------------------

The usual way to create a ``Dataset`` filled with Dask arrays is to load the
data from a netCDF file or files. You can do this by supplying a ``chunks``
argument to :py:func:`~xarray.open_dataset` or using the
:py:func:`~xarray.open_mfdataset` function.

.. ipython:: python
    :suppress:

    import numpy as np
    import pandas as pd
    import xarray as xr

    np.random.seed(123456)
    np.set_printoptions(precision=3, linewidth=100, threshold=100, edgeitems=3)

    ds = xr.Dataset(
        {
            "temperature": (
                ("time", "latitude", "longitude"),
                np.random.randn(30, 180, 180),
            ),
            "time": pd.date_range("2015-01-01", periods=30),
            "longitude": np.arange(180),
            "latitude": np.arange(89.5, -90.5, -1),
        }
    )
    ds.to_netcdf("example-data.nc")

.. ipython:: python

    ds = xr.open_dataset("example-data.nc", chunks={"time": 10})
    ds

In this example ``latitude`` and ``longitude`` do not appear in the ``chunks``
dict, so only one chunk will be used along those dimensions.  It is also
entirely equivalent to opening a dataset using :py:meth:`~xarray.open_dataset`
and then chunking the data using the ``chunk`` method, e.g.,
``xr.open_dataset('example-data.nc').chunk({'time': 10})``.

To open multiple files simultaneously in parallel using Dask delayed,
use :py:func:`~xarray.open_mfdataset`::

    xr.open_mfdataset('my/files/*.nc', parallel=True)

This function will automatically concatenate and merge datasets into one in
the simple cases that it understands (see :py:func:`~xarray.combine_by_coords`
for the full disclaimer). By default, :py:meth:`~xarray.open_mfdataset` will chunk each
netCDF file into a single Dask array; again, supply the ``chunks`` argument to
control the size of the resulting Dask arrays. In more complex cases, you can
open each file individually using :py:meth:`~xarray.open_dataset` and merge the result, as
described in :ref:`combining data`. Passing the keyword argument ``parallel=True`` to :py:meth:`~xarray.open_mfdataset` will speed up the reading of large multi-file datasets by
executing those read tasks in parallel using ``dask.delayed``.

You'll notice that printing a dataset still shows a preview of array values,
even if they are actually Dask arrays. We can do this quickly with Dask because
we only need to compute the first few values (typically from the first block).
To reveal the true nature of an array, print a DataArray:

.. ipython:: python

    ds.temperature

Once you've manipulated a Dask array, you can still write a dataset too big to
fit into memory back to disk by using :py:meth:`~xarray.Dataset.to_netcdf` in the
usual way.

.. ipython:: python

    ds.to_netcdf("manipulated-example-data.nc")

By setting the ``compute`` argument to ``False``, :py:meth:`~xarray.Dataset.to_netcdf`
will return a ``dask.delayed`` object that can be computed later.

.. ipython:: python

    from dask.diagnostics import ProgressBar

    # or distributed.progress when using the distributed scheduler
    delayed_obj = ds.to_netcdf("manipulated-example-data.nc", compute=False)
    with ProgressBar():
        results = delayed_obj.compute()

.. note::

    When using Dask's distributed scheduler to write NETCDF4 files,
    it may be necessary to set the environment variable `HDF5_USE_FILE_LOCKING=FALSE`
    to avoid competing locks within the HDF5 SWMR file locking scheme. Note that
    writing netCDF files with Dask's distributed scheduler is only supported for
    the `netcdf4` backend.

A dataset can also be converted to a Dask DataFrame using :py:meth:`~xarray.Dataset.to_dask_dataframe`.

.. ipython:: python
    :okwarning:

    df = ds.to_dask_dataframe()
    df

Dask DataFrames do not support multi-indexes so the coordinate variables from the dataset are included as columns in the Dask DataFrame.

.. ipython:: python
    :suppress:

    import os

    os.remove("example-data.nc")
    os.remove("manipulated-example-data.nc")

Using Dask with xarray
----------------------

Nearly all existing xarray methods (including those for indexing, computation,
concatenating and grouped operations) have been extended to work automatically
with Dask arrays. When you load data as a Dask array in an xarray data
structure, almost all xarray operations will keep it as a Dask array; when this
is not possible, they will raise an exception rather than unexpectedly loading
data into memory. Converting a Dask array into memory generally requires an
explicit conversion step. One notable exception is indexing operations: to
enable label based indexing, xarray will automatically load coordinate labels
into memory.

.. tip::

   By default, dask uses its multi-threaded scheduler, which distributes work across
   multiple cores and allows for processing some datasets that do not fit into memory.
   For running across a cluster, `setup the distributed scheduler <https://docs.dask.org/en/latest/setup.html>`_.

The easiest way to convert an xarray data structure from lazy Dask arrays into
*eager*, in-memory NumPy arrays is to use the :py:meth:`~xarray.Dataset.load` method:

.. ipython:: python

    ds.load()

You can also access :py:attr:`~xarray.DataArray.values`, which will always be a
NumPy array:

.. ipython::
    :verbatim:

    In [5]: ds.temperature.values
    Out[5]:
    array([[[  4.691e-01,  -2.829e-01, ...,  -5.577e-01,   3.814e-01],
            [  1.337e+00,  -1.531e+00, ...,   8.726e-01,  -1.538e+00],
            ...
    # truncated for brevity

Explicit conversion by wrapping a DataArray with ``np.asarray`` also works:

.. ipython::
    :verbatim:

    In [5]: np.asarray(ds.temperature)
    Out[5]:
    array([[[  4.691e-01,  -2.829e-01, ...,  -5.577e-01,   3.814e-01],
            [  1.337e+00,  -1.531e+00, ...,   8.726e-01,  -1.538e+00],
            ...

Alternatively you can load the data into memory but keep the arrays as
Dask arrays using the :py:meth:`~xarray.Dataset.persist` method:

.. ipython:: python

    ds = ds.persist()

:py:meth:`~xarray.Dataset.persist` is particularly useful when using a
distributed cluster because the data will be loaded into distributed memory
across your machines and be much faster to use than reading repeatedly from
disk.

.. warning::

   On a single machine :py:meth:`~xarray.Dataset.persist` will try to load all of
   your data into memory. You should make sure that your dataset is not larger than
   available memory.

.. note::
   For more on the differences between :py:meth:`~xarray.Dataset.persist` and
   :py:meth:`~xarray.Dataset.compute` see this `Stack Overflow answer <https://stackoverflow.com/questions/41806850/dask-difference-between-client-persist-and-client-compute>`_ and the `Dask documentation <https://distributed.readthedocs.io/en/latest/manage-computation.html#dask-collections-to-futures>`_.

For performance you may wish to consider chunk sizes.  The correct choice of
chunk size depends both on your data and on the operations you want to perform.
With xarray, both converting data to a Dask arrays and converting the chunk
sizes of Dask arrays is done with the :py:meth:`~xarray.Dataset.chunk` method:

.. ipython:: python
    :suppress:

    ds = ds.chunk({"time": 10})

.. ipython:: python

    rechunked = ds.chunk({"latitude": 100, "longitude": 100})

You can view the size of existing chunks on an array by viewing the
:py:attr:`~xarray.Dataset.chunks` attribute:

.. ipython:: python

    rechunked.chunks

If there are not consistent chunksizes between all the arrays in a dataset
along a particular dimension, an exception is raised when you try to access
``.chunks``.

.. note::

    In the future, we would like to enable automatic alignment of Dask
    chunksizes (but not the other way around). We might also require that all
    arrays in a dataset share the same chunking alignment. Neither of these
    are currently done.

NumPy ufuncs like ``np.sin`` transparently work on all xarray objects, including those
that store lazy Dask arrays:

.. ipython:: python

    import numpy as np

    np.sin(rechunked)

To access Dask arrays directly, use the
:py:attr:`DataArray.data <xarray.DataArray.data>` attribute. This attribute exposes
array data either as a Dask array or as a NumPy array, depending on whether it has been
loaded into Dask or not:

.. ipython:: python

    ds.temperature.data

.. note::

    ``.data`` is also used to expose other "computable" array backends beyond Dask and
    NumPy (e.g. sparse and pint arrays).

.. _dask.automatic-parallelization:

Automatic parallelization with ``apply_ufunc`` and ``map_blocks``
-----------------------------------------------------------------

Almost all of xarray's built-in operations work on Dask arrays. If you want to
use a function that isn't wrapped by xarray, and have it applied in parallel on
each block of your xarray object, you have three options:

1. Extract Dask arrays from xarray objects (``.data``) and use Dask directly.
2. Use :py:func:`~xarray.apply_ufunc` to apply functions that consume and return NumPy arrays.
3. Use :py:func:`~xarray.map_blocks`, :py:meth:`Dataset.map_blocks` or :py:meth:`DataArray.map_blocks`
   to apply functions that consume and return xarray objects.


``apply_ufunc``
~~~~~~~~~~~~~~~

Another option is to use xarray's :py:func:`~xarray.apply_ufunc`, which can
automate `embarrassingly parallel
<https://en.wikipedia.org/wiki/Embarrassingly_parallel>`__ "map" type operations
where a function written for processing NumPy arrays should be repeatedly
applied to xarray objects containing Dask arrays. It works similarly to
:py:func:`dask.array.map_blocks` and :py:func:`dask.array.blockwise`, but without
requiring an intermediate layer of abstraction.

For the best performance when using Dask's multi-threaded scheduler, wrap a
function that already releases the global interpreter lock, which fortunately
already includes most NumPy and Scipy functions. Here we show an example
using NumPy operations and a fast function from
`bottleneck <https://github.com/pydata/bottleneck>`__, which
we use to calculate `Spearman's rank-correlation coefficient <https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient>`__:

.. code-block:: python

    import numpy as np
    import xarray as xr
    import bottleneck


    def covariance_gufunc(x, y):
        return (
            (x - x.mean(axis=-1, keepdims=True)) * (y - y.mean(axis=-1, keepdims=True))
        ).mean(axis=-1)


    def pearson_correlation_gufunc(x, y):
        return covariance_gufunc(x, y) / (x.std(axis=-1) * y.std(axis=-1))


    def spearman_correlation_gufunc(x, y):
        x_ranks = bottleneck.rankdata(x, axis=-1)
        y_ranks = bottleneck.rankdata(y, axis=-1)
        return pearson_correlation_gufunc(x_ranks, y_ranks)


    def spearman_correlation(x, y, dim):
        return xr.apply_ufunc(
            spearman_correlation_gufunc,
            x,
            y,
            input_core_dims=[[dim], [dim]],
            dask="parallelized",
            output_dtypes=[float],
        )

The only aspect of this example that is different from standard usage of
``apply_ufunc()`` is that we needed to supply the ``output_dtypes`` arguments.
(Read up on :ref:`comput.wrapping-custom` for an explanation of the
"core dimensions" listed in ``input_core_dims``.)

Our new ``spearman_correlation()`` function achieves near linear speedup
when run on large arrays across the four cores on my laptop. It would also
work as a streaming operation, when run on arrays loaded from disk:

.. ipython::
    :verbatim:

    In [56]: rs = np.random.RandomState(0)

    In [57]: array1 = xr.DataArray(rs.randn(1000, 100000), dims=["place", "time"])  # 800MB

    In [58]: array2 = array1 + 0.5 * rs.randn(1000, 100000)

    # using one core, on NumPy arrays
    In [61]: %time _ = spearman_correlation(array1, array2, 'time')
    CPU times: user 21.6 s, sys: 2.84 s, total: 24.5 s
    Wall time: 24.9 s

    In [8]: chunked1 = array1.chunk({"place": 10})

    In [9]: chunked2 = array2.chunk({"place": 10})

    # using all my laptop's cores, with Dask
    In [63]: r = spearman_correlation(chunked1, chunked2, "time").compute()

    In [64]: %time _ = r.compute()
    CPU times: user 30.9 s, sys: 1.74 s, total: 32.6 s
    Wall time: 4.59 s

One limitation of ``apply_ufunc()`` is that it cannot be applied to arrays with
multiple chunks along a core dimension:

.. ipython::
    :verbatim:

    In [63]: spearman_correlation(chunked1, chunked2, "place")
    ValueError: dimension 'place' on 0th function argument to apply_ufunc with
    dask='parallelized' consists of multiple chunks, but is also a core
    dimension. To fix, rechunk into a single Dask array chunk along this
    dimension, i.e., ``.rechunk({'place': -1})``, but beware that this may
    significantly increase memory usage.

This reflects the nature of core dimensions, in contrast to broadcast (non-core)
dimensions that allow operations to be split into arbitrary chunks for
application.

.. tip::

    For the majority of NumPy functions that are already wrapped by Dask, it's
    usually a better idea to use the pre-existing ``dask.array`` function, by
    using either a pre-existing xarray methods or
    :py:func:`~xarray.apply_ufunc()` with ``dask='allowed'``. Dask can often
    have a more efficient implementation that makes use of the specialized
    structure of a problem, unlike the generic speedups offered by
    ``dask='parallelized'``.


``map_blocks``
~~~~~~~~~~~~~~

Functions that consume and return xarray objects can be easily applied in parallel using :py:func:`map_blocks`.
Your function will receive an xarray Dataset or DataArray subset to one chunk
along each chunked dimension.

.. ipython:: python

    ds.temperature

This DataArray has 3 chunks each with length 10 along the time dimension.
At compute time, a function applied with :py:func:`map_blocks` will receive a DataArray corresponding to a single block of shape 10x180x180
(time x latitude x longitude) with values loaded. The following snippet illustrates how to check the shape of the object
received by the applied function.

.. ipython:: python

    def func(da):
        print(da.sizes)
        return da.time


    mapped = xr.map_blocks(func, ds.temperature)
    mapped

Notice that the :py:meth:`map_blocks` call printed
``Frozen({'time': 0, 'latitude': 0, 'longitude': 0})`` to screen.
``func`` is received 0-sized blocks! :py:meth:`map_blocks` needs to know what the final result
looks like in terms of dimensions, shapes etc. It does so by running the provided function on 0-shaped
inputs (*automated inference*). This works in many cases, but not all. If automatic inference does not
work for your function, provide the ``template`` kwarg (see below).

In this case, automatic inference has worked so let's check that the result is as expected.

.. ipython:: python

    mapped.load(scheduler="single-threaded")
    mapped.identical(ds.time)

Note that we use ``.load(scheduler="single-threaded")`` to execute the computation.
This executes the Dask graph in `serial` using a for loop, but allows for printing to screen and other
debugging techniques. We can easily see that our function is receiving blocks of shape 10x180x180 and
the returned result is identical to ``ds.time`` as expected.


Here is a common example where automated inference will not work.

.. ipython:: python
    :okexcept:

    def func(da):
        print(da.sizes)
        return da.isel(time=[1])


    mapped = xr.map_blocks(func, ds.temperature)

``func`` cannot be run on 0-shaped inputs because it is not possible to extract element 1 along a
dimension of size 0. In this case we need to tell :py:func:`map_blocks` what the returned result looks
like using the ``template`` kwarg. ``template`` must be an xarray Dataset or DataArray (depending on
what the function returns) with dimensions, shapes, chunk sizes, attributes, coordinate variables *and* data
variables that look exactly like the expected result. The variables should be dask-backed and hence not
incur much memory cost.

.. note::

    Note that when ``template`` is provided, ``attrs`` from ``template`` are copied over to the result. Any
    ``attrs`` set in ``func`` will be ignored.


.. ipython:: python

    template = ds.temperature.isel(time=[1, 11, 21])
    mapped = xr.map_blocks(func, ds.temperature, template=template)


Notice that the 0-shaped sizes were not printed to screen. Since ``template`` has been provided
:py:func:`map_blocks` does not need to infer it by running ``func`` on 0-shaped inputs.

.. ipython:: python

    mapped.identical(template)


:py:func:`map_blocks` also allows passing ``args`` and ``kwargs`` down to the user function ``func``.
``func`` will be executed as ``func(block_xarray, *args, **kwargs)`` so ``args`` must be a list and ``kwargs`` must be a dictionary.

.. ipython:: python

    def func(obj, a, b=0):
        return obj + a + b


    mapped = ds.map_blocks(func, args=[10], kwargs={"b": 10})
    expected = ds + 10 + 10
    mapped.identical(expected)


.. tip::

   As :py:func:`map_blocks` loads each block into memory, reduce as much as possible objects consumed by user functions.
   For example, drop useless variables before calling ``func`` with :py:func:`map_blocks`.



Chunking and performance
------------------------

The ``chunks`` parameter has critical performance implications when using Dask
arrays. If your chunks are too small, queueing up operations will be extremely
slow, because Dask will translate each operation into a huge number of
operations mapped across chunks. Computation on Dask arrays with small chunks
can also be slow, because each operation on a chunk has some fixed overhead from
the Python interpreter and the Dask task executor.

Conversely, if your chunks are too big, some of your computation may be wasted,
because Dask only computes results one chunk at a time.

A good rule of thumb is to create arrays with a minimum chunksize of at least
one million elements (e.g., a 1000x1000 matrix). With large arrays (10+ GB), the
cost of queueing up Dask operations can be noticeable, and you may need even
larger chunksizes.

.. tip::

   Check out the dask documentation on `chunks <https://docs.dask.org/en/latest/array-chunks.html>`_.


Optimization Tips
-----------------

With analysis pipelines involving both spatial subsetting and temporal resampling, Dask performance can become very slow in certain cases. Here are some optimization tips we have found through experience:

1. Do your spatial and temporal indexing (e.g. ``.sel()`` or ``.isel()``) early in the pipeline, especially before calling ``resample()`` or ``groupby()``. Grouping and resampling triggers some computation on all the blocks, which in theory should commute with indexing, but this optimization hasn't been implemented in Dask yet. (See `Dask issue #746 <https://github.com/dask/dask/issues/746>`_).

2. Save intermediate results to disk as a netCDF files (using ``to_netcdf()``) and then load them again with ``open_dataset()`` for further computations. For example, if subtracting temporal mean from a dataset, save the temporal mean to disk before subtracting. Again, in theory, Dask should be able to do the computation in a streaming fashion, but in practice this is a fail case for the Dask scheduler, because it tries to keep every chunk of an array that it computes in memory. (See `Dask issue #874 <https://github.com/dask/dask/issues/874>`_)

3. Specify smaller chunks across space when using :py:meth:`~xarray.open_mfdataset` (e.g., ``chunks={'latitude': 10, 'longitude': 10}``). This makes spatial subsetting easier, because there's no risk you will load chunks of data referring to different chunks (probably not necessary if you follow suggestion 1).

4. Using the h5netcdf package by passing ``engine='h5netcdf'`` to :py:meth:`~xarray.open_mfdataset`
   can be quicker than the default ``engine='netcdf4'`` that uses the netCDF4 package.

5. Some dask-specific tips may be found `here <https://docs.dask.org/en/latest/array-best-practices.html>`_.

6. The dask `diagnostics <https://docs.dask.org/en/latest/understanding-performance.html>`_ can be
   useful in identifying performance bottlenecks.
.. _reshape:

###############################
Reshaping and reorganizing data
###############################

These methods allow you to reorganize

.. ipython:: python
    :suppress:

    import numpy as np
    import pandas as pd
    import xarray as xr

    np.random.seed(123456)

Reordering dimensions
---------------------

To reorder dimensions on a :py:class:`~xarray.DataArray` or across all variables
on a :py:class:`~xarray.Dataset`, use :py:meth:`~xarray.DataArray.transpose`. An
ellipsis (`...`) can be use to represent all other dimensions:

.. ipython:: python

    ds = xr.Dataset({"foo": (("x", "y", "z"), [[[42]]]), "bar": (("y", "z"), [[24]])})
    ds.transpose("y", "z", "x")
    ds.transpose(..., "x")  # equivalent
    ds.transpose()  # reverses all dimensions

Expand and squeeze dimensions
-----------------------------

To expand a :py:class:`~xarray.DataArray` or all
variables on a :py:class:`~xarray.Dataset` along a new dimension,
use :py:meth:`~xarray.DataArray.expand_dims`

.. ipython:: python

    expanded = ds.expand_dims("w")
    expanded

This method attaches a new dimension with size 1 to all data variables.

To remove such a size-1 dimension from the :py:class:`~xarray.DataArray`
or :py:class:`~xarray.Dataset`,
use :py:meth:`~xarray.DataArray.squeeze`

.. ipython:: python

    expanded.squeeze("w")

Converting between datasets and arrays
--------------------------------------

To convert from a Dataset to a DataArray, use :py:meth:`~xarray.Dataset.to_array`:

.. ipython:: python

    arr = ds.to_array()
    arr

This method broadcasts all data variables in the dataset against each other,
then concatenates them along a new dimension into a new array while preserving
coordinates.

To convert back from a DataArray to a Dataset, use
:py:meth:`~xarray.DataArray.to_dataset`:

.. ipython:: python

    arr.to_dataset(dim="variable")

The broadcasting behavior of ``to_array`` means that the resulting array
includes the union of data variable dimensions:

.. ipython:: python

    ds2 = xr.Dataset({"a": 0, "b": ("x", [3, 4, 5])})

    # the input dataset has 4 elements
    ds2

    # the resulting array has 6 elements
    ds2.to_array()

Otherwise, the result could not be represented as an orthogonal array.

If you use ``to_dataset`` without supplying the ``dim`` argument, the DataArray will be converted into a Dataset of one variable:

.. ipython:: python

    arr.to_dataset(name="combined")

.. _reshape.stack:

Stack and unstack
-----------------

As part of xarray's nascent support for :py:class:`pandas.MultiIndex`, we have
implemented :py:meth:`~xarray.DataArray.stack` and
:py:meth:`~xarray.DataArray.unstack` method, for combining or splitting dimensions:

.. ipython:: python

    array = xr.DataArray(
        np.random.randn(2, 3), coords=[("x", ["a", "b"]), ("y", [0, 1, 2])]
    )
    stacked = array.stack(z=("x", "y"))
    stacked
    stacked.unstack("z")

As elsewhere in xarray, an ellipsis (`...`) can be used to represent all unlisted dimensions:

.. ipython:: python

    stacked = array.stack(z=[..., "x"])
    stacked

These methods are modeled on the :py:class:`pandas.DataFrame` methods of the
same name, although in xarray they always create new dimensions rather than
adding to the existing index or columns.

Like :py:meth:`DataFrame.unstack<pandas.DataFrame.unstack>`, xarray's ``unstack``
always succeeds, even if the multi-index being unstacked does not contain all
possible levels. Missing levels are filled in with ``NaN`` in the resulting object:

.. ipython:: python

    stacked2 = stacked[::2]
    stacked2
    stacked2.unstack("z")

However, xarray's ``stack`` has an important difference from pandas: unlike
pandas, it does not automatically drop missing values. Compare:

.. ipython:: python

    array = xr.DataArray([[np.nan, 1], [2, 3]], dims=["x", "y"])
    array.stack(z=("x", "y"))
    array.to_pandas().stack()

We departed from pandas's behavior here because predictable shapes for new
array dimensions is necessary for :ref:`dask`.

.. _reshape.stacking_different:

Stacking different variables together
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These stacking and unstacking operations are particularly useful for reshaping
xarray objects for use in machine learning packages, such as `scikit-learn
<http://scikit-learn.org/stable/>`_, that usually require two-dimensional numpy
arrays as inputs. For datasets with only one variable, we only need ``stack``
and ``unstack``, but combining multiple variables in a
:py:class:`xarray.Dataset` is more complicated. If the variables in the dataset
have matching numbers of dimensions, we can call
:py:meth:`~xarray.Dataset.to_array` and then stack along the the new coordinate.
But :py:meth:`~xarray.Dataset.to_array` will broadcast the dataarrays together,
which will effectively tile the lower dimensional variable along the missing
dimensions. The method :py:meth:`xarray.Dataset.to_stacked_array` allows
combining variables of differing dimensions without this wasteful copying while
:py:meth:`xarray.DataArray.to_unstacked_dataset` reverses this operation.
Just as with :py:meth:`xarray.Dataset.stack` the stacked coordinate is
represented by a :py:class:`pandas.MultiIndex` object. These methods are used
like this:

.. ipython:: python

    data = xr.Dataset(
        data_vars={"a": (("x", "y"), [[0, 1, 2], [3, 4, 5]]), "b": ("x", [6, 7])},
        coords={"y": ["u", "v", "w"]},
    )
    data
    stacked = data.to_stacked_array("z", sample_dims=["x"])
    stacked
    unstacked = stacked.to_unstacked_dataset("z")
    unstacked

In this example, ``stacked`` is a two dimensional array that we can easily pass to a scikit-learn or another generic
numerical method.

.. note::

    Unlike with ``stack``,  in ``to_stacked_array``, the user specifies the dimensions they **do not** want stacked.
    For a machine learning task, these unstacked dimensions can be interpreted as the dimensions over which samples are
    drawn, whereas the stacked coordinates are the features. Naturally, all variables should possess these sampling
    dimensions.


.. _reshape.set_index:

Set and reset index
-------------------

Complementary to stack / unstack, xarray's ``.set_index``, ``.reset_index`` and
``.reorder_levels`` allow easy manipulation of ``DataArray`` or ``Dataset``
multi-indexes without modifying the data and its dimensions.

You can create a multi-index from several 1-dimensional variables and/or
coordinates using :py:meth:`~xarray.DataArray.set_index`:

.. ipython:: python

    da = xr.DataArray(
        np.random.rand(4),
        coords={
            "band": ("x", ["a", "a", "b", "b"]),
            "wavenumber": ("x", np.linspace(200, 400, 4)),
        },
        dims="x",
    )
    da
    mda = da.set_index(x=["band", "wavenumber"])
    mda

These coordinates can now be used for indexing, e.g.,

.. ipython:: python

    mda.sel(band="a")

Conversely, you can use :py:meth:`~xarray.DataArray.reset_index`
to extract multi-index levels as coordinates (this is mainly useful
for serialization):

.. ipython:: python

    mda.reset_index("x")

:py:meth:`~xarray.DataArray.reorder_levels` allows changing the order
of multi-index levels:

.. ipython:: python

    mda.reorder_levels(x=["wavenumber", "band"])

As of xarray v0.9 coordinate labels for each dimension are optional.
You can also use ``.set_index`` / ``.reset_index`` to add / remove
labels for one or several dimensions:

.. ipython:: python

    array = xr.DataArray([1, 2, 3], dims="x")
    array
    array["c"] = ("x", ["a", "b", "c"])
    array.set_index(x="c")
    array = array.set_index(x="c")
    array = array.reset_index("x", drop=True)

.. _reshape.shift_and_roll:

Shift and roll
--------------

To adjust coordinate labels, you can use the :py:meth:`~xarray.Dataset.shift` and
:py:meth:`~xarray.Dataset.roll` methods:

.. ipython:: python

    array = xr.DataArray([1, 2, 3, 4], dims="x")
    array.shift(x=2)
    array.roll(x=2, roll_coords=True)

.. _reshape.sort:

Sort
----

One may sort a DataArray/Dataset via :py:meth:`~xarray.DataArray.sortby` and
:py:meth:`~xarray.DataArray.sortby`.  The input can be an individual or list of
1D ``DataArray`` objects:

.. ipython:: python

    ds = xr.Dataset(
        {
            "A": (("x", "y"), [[1, 2], [3, 4]]),
            "B": (("x", "y"), [[5, 6], [7, 8]]),
        },
        coords={"x": ["b", "a"], "y": [1, 0]},
    )
    dax = xr.DataArray([100, 99], [("x", [0, 1])])
    day = xr.DataArray([90, 80], [("y", [0, 1])])
    ds.sortby([day, dax])

As a shortcut, you can refer to existing coordinates by name:

.. ipython:: python

    ds.sortby("x")
    ds.sortby(["y", "x"])
    ds.sortby(["y", "x"], ascending=False)
.. _indexing:

Indexing and selecting data
===========================

.. ipython:: python
    :suppress:

    import numpy as np
    import pandas as pd
    import xarray as xr

    np.random.seed(123456)

Xarray offers extremely flexible indexing routines that combine the best
features of NumPy and pandas for data selection.

The most basic way to access elements of a :py:class:`~xarray.DataArray`
object is to use Python's ``[]`` syntax, such as ``array[i, j]``, where
``i`` and ``j`` are both integers.
As xarray objects can store coordinates corresponding to each dimension of an
array, label-based indexing similar to ``pandas.DataFrame.loc`` is also possible.
In label-based indexing, the element position ``i`` is automatically
looked-up from the coordinate values.

Dimensions of xarray objects have names, so you can also lookup the dimensions
by name, instead of remembering their positional order.

Thus in total, xarray supports four different kinds of indexing, as described
below and summarized in this table:

.. |br| raw:: html

   <br />

+------------------+--------------+---------------------------------+--------------------------------+
| Dimension lookup | Index lookup | ``DataArray`` syntax            | ``Dataset`` syntax             |
+==================+==============+=================================+================================+
| Positional       | By integer   | ``da[:, 0]``                    | *not available*                |
+------------------+--------------+---------------------------------+--------------------------------+
| Positional       | By label     | ``da.loc[:, 'IA']``             | *not available*                |
+------------------+--------------+---------------------------------+--------------------------------+
| By name          | By integer   | ``da.isel(space=0)`` or |br|    | ``ds.isel(space=0)`` or |br|   |
|                  |              | ``da[dict(space=0)]``           | ``ds[dict(space=0)]``          |
+------------------+--------------+---------------------------------+--------------------------------+
| By name          | By label     | ``da.sel(space='IA')`` or |br|  | ``ds.sel(space='IA')`` or |br| |
|                  |              | ``da.loc[dict(space='IA')]``    | ``ds.loc[dict(space='IA')]``   |
+------------------+--------------+---------------------------------+--------------------------------+

More advanced indexing is also possible for all the methods by
supplying :py:class:`~xarray.DataArray` objects as indexer.
See :ref:`vectorized_indexing` for the details.


Positional indexing
-------------------

Indexing a :py:class:`~xarray.DataArray` directly works (mostly) just like it
does for numpy arrays, except that the returned object is always another
DataArray:

.. ipython:: python

    da = xr.DataArray(
        np.random.rand(4, 3),
        [
            ("time", pd.date_range("2000-01-01", periods=4)),
            ("space", ["IA", "IL", "IN"]),
        ],
    )
    da[:2]
    da[0, 0]
    da[:, [2, 1]]

Attributes are persisted in all indexing operations.

.. warning::

    Positional indexing deviates from the NumPy when indexing with multiple
    arrays like ``da[[0, 1], [0, 1]]``, as described in
    :ref:`vectorized_indexing`.

Xarray also supports label-based indexing, just like pandas. Because
we use a :py:class:`pandas.Index` under the hood, label based indexing is very
fast. To do label based indexing, use the :py:attr:`~xarray.DataArray.loc` attribute:

.. ipython:: python

    da.loc["2000-01-01":"2000-01-02", "IA"]

In this example, the selected is a subpart of the array
in the range '2000-01-01':'2000-01-02' along the first coordinate `time`
and with 'IA' value from the second coordinate `space`.

You can perform any of the label indexing operations `supported by pandas`__,
including indexing with individual, slices and arrays of labels, as well as
indexing with boolean arrays. Like pandas, label based indexing in xarray is
*inclusive* of both the start and stop bounds.

__ http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-label

Setting values with label based indexing is also supported:

.. ipython:: python

    da.loc["2000-01-01", ["IL", "IN"]] = -10
    da


Indexing with dimension names
-----------------------------

With the dimension names, we do not have to rely on dimension order and can
use them explicitly to slice data. There are two ways to do this:

1. Use a dictionary as the argument for array positional or label based array
   indexing:

    .. ipython:: python

        # index by integer array indices
        da[dict(space=0, time=slice(None, 2))]

        # index by dimension coordinate labels
        da.loc[dict(time=slice("2000-01-01", "2000-01-02"))]

2. Use the :py:meth:`~xarray.DataArray.sel` and :py:meth:`~xarray.DataArray.isel`
   convenience methods:

    .. ipython:: python

        # index by integer array indices
        da.isel(space=0, time=slice(None, 2))

        # index by dimension coordinate labels
        da.sel(time=slice("2000-01-01", "2000-01-02"))

The arguments to these methods can be any objects that could index the array
along the dimension given by the keyword, e.g., labels for an individual value,
Python :py:class:`slice` objects or 1-dimensional arrays.

.. note::

    We would love to be able to do indexing with labeled dimension names inside
    brackets, but unfortunately, Python `does yet not support`__ indexing with
    keyword arguments like ``da[space=0]``

__ http://legacy.python.org/dev/peps/pep-0472/


.. _nearest neighbor lookups:

Nearest neighbor lookups
------------------------

The label based selection methods :py:meth:`~xarray.Dataset.sel`,
:py:meth:`~xarray.Dataset.reindex` and :py:meth:`~xarray.Dataset.reindex_like` all
support ``method`` and ``tolerance`` keyword argument. The method parameter allows for
enabling nearest neighbor (inexact) lookups by use of the methods ``'pad'``,
``'backfill'`` or ``'nearest'``:

.. ipython:: python

    da = xr.DataArray([1, 2, 3], [("x", [0, 1, 2])])
    da.sel(x=[1.1, 1.9], method="nearest")
    da.sel(x=0.1, method="backfill")
    da.reindex(x=[0.5, 1, 1.5, 2, 2.5], method="pad")

Tolerance limits the maximum distance for valid matches with an inexact lookup:

.. ipython:: python

    da.reindex(x=[1.1, 1.5], method="nearest", tolerance=0.2)

The method parameter is not yet supported if any of the arguments
to ``.sel()`` is a ``slice`` object:

.. ipython::
   :verbatim:

   In [1]: da.sel(x=slice(1, 3), method="nearest")
   NotImplementedError

However, you don't need to use ``method`` to do inexact slicing. Slicing
already returns all values inside the range (inclusive), as long as the index
labels are monotonic increasing:

.. ipython:: python

    da.sel(x=slice(0.9, 3.1))

Indexing axes with monotonic decreasing labels also works, as long as the
``slice`` or ``.loc`` arguments are also decreasing:

.. ipython:: python

    reversed_da = da[::-1]
    reversed_da.loc[3.1:0.9]


.. note::

  If you want to interpolate along coordinates rather than looking up the
  nearest neighbors, use :py:meth:`~xarray.Dataset.interp` and
  :py:meth:`~xarray.Dataset.interp_like`.
  See :ref:`interpolation <interp>` for the details.


Dataset indexing
----------------

We can also use these methods to index all variables in a dataset
simultaneously, returning a new dataset:

.. ipython:: python

    da = xr.DataArray(
        np.random.rand(4, 3),
        [
            ("time", pd.date_range("2000-01-01", periods=4)),
            ("space", ["IA", "IL", "IN"]),
        ],
    )
    ds = da.to_dataset(name="foo")
    ds.isel(space=[0], time=[0])
    ds.sel(time="2000-01-01")

Positional indexing on a dataset is not supported because the ordering of
dimensions in a dataset is somewhat ambiguous (it can vary between different
arrays). However, you can do normal indexing with dimension names:

.. ipython:: python

    ds[dict(space=[0], time=[0])]
    ds.loc[dict(time="2000-01-01")]

Dropping labels and dimensions
------------------------------

The :py:meth:`~xarray.Dataset.drop_sel` method returns a new object with the listed
index labels along a dimension dropped:

.. ipython:: python

    ds.drop_sel(space=["IN", "IL"])

``drop_sel`` is both a ``Dataset`` and ``DataArray`` method.

Use :py:meth:`~xarray.Dataset.drop_dims` to drop a full dimension from a Dataset.
Any variables with these dimensions are also dropped:

.. ipython:: python

    ds.drop_dims("time")

.. _masking with where:

Masking with ``where``
----------------------

Indexing methods on xarray objects generally return a subset of the original data.
However, it is sometimes useful to select an object with the same shape as the
original data, but with some elements masked. To do this type of selection in
xarray, use :py:meth:`~xarray.DataArray.where`:

.. ipython:: python

    da = xr.DataArray(np.arange(16).reshape(4, 4), dims=["x", "y"])
    da.where(da.x + da.y < 4)

This is particularly useful for ragged indexing of multi-dimensional data,
e.g., to apply a 2D mask to an image. Note that ``where`` follows all the
usual xarray broadcasting and alignment rules for binary operations (e.g.,
``+``) between the object being indexed and the condition, as described in
:ref:`comput`:

.. ipython:: python

    da.where(da.y < 2)

By default ``where`` maintains the original size of the data.  For cases
where the selected data size is much smaller than the original data,
use of the option ``drop=True`` clips coordinate
elements that are fully masked:

.. ipython:: python

    da.where(da.y < 2, drop=True)

.. _selecting values with isin:

Selecting values with ``isin``
------------------------------

To check whether elements of an xarray object contain a single object, you can
compare with the equality operator ``==`` (e.g., ``arr == 3``). To check
multiple values, use :py:meth:`~xarray.DataArray.isin`:

.. ipython:: python

    da = xr.DataArray([1, 2, 3, 4, 5], dims=["x"])
    da.isin([2, 4])

:py:meth:`~xarray.DataArray.isin` works particularly well with
:py:meth:`~xarray.DataArray.where` to support indexing by arrays that are not
already labels of an array:

.. ipython:: python

    lookup = xr.DataArray([-1, -2, -3, -4, -5], dims=["x"])
    da.where(lookup.isin([-2, -4]), drop=True)

However, some caution is in order: when done repeatedly, this type of indexing
is significantly slower than using :py:meth:`~xarray.DataArray.sel`.

.. _vectorized_indexing:

Vectorized Indexing
-------------------

Like numpy and pandas, xarray supports indexing many array elements at once in a
`vectorized` manner.

If you only provide integers, slices, or unlabeled arrays (array without
dimension names, such as ``np.ndarray``, ``list``, but not
:py:meth:`~xarray.DataArray` or :py:meth:`~xarray.Variable`) indexing can be
understood as orthogonally. Each indexer component selects independently along
the corresponding dimension, similar to how vector indexing works in Fortran or
MATLAB, or after using the :py:func:`numpy.ix_` helper:

.. ipython:: python

    da = xr.DataArray(
        np.arange(12).reshape((3, 4)),
        dims=["x", "y"],
        coords={"x": [0, 1, 2], "y": ["a", "b", "c", "d"]},
    )
    da
    da[[0, 2, 2], [1, 3]]

For more flexibility, you can supply :py:meth:`~xarray.DataArray` objects
as indexers.
Dimensions on resultant arrays are given by the ordered union of the indexers'
dimensions:

.. ipython:: python

    ind_x = xr.DataArray([0, 1], dims=["x"])
    ind_y = xr.DataArray([0, 1], dims=["y"])
    da[ind_x, ind_y]  # orthogonal indexing
    da[ind_x, ind_x]  # vectorized indexing

Slices or sequences/arrays without named-dimensions are treated as if they have
the same dimension which is indexed along:

.. ipython:: python

    # Because [0, 1] is used to index along dimension 'x',
    # it is assumed to have dimension 'x'
    da[[0, 1], ind_x]

Furthermore, you can use multi-dimensional :py:meth:`~xarray.DataArray`
as indexers, where the resultant array dimension is also determined by
indexers' dimension:

.. ipython:: python

    ind = xr.DataArray([[0, 1], [0, 1]], dims=["a", "b"])
    da[ind]

Similar to how NumPy's `advanced indexing`_ works, vectorized
indexing for xarray is based on our
:ref:`broadcasting rules <compute.broadcasting>`.
See :ref:`indexing.rules` for the complete specification.

.. _advanced indexing: https://docs.scipy.org/doc/numpy-1.13.0/reference/arrays.indexing.html

Vectorized indexing also works with ``isel``, ``loc``, and ``sel``:

.. ipython:: python

    ind = xr.DataArray([[0, 1], [0, 1]], dims=["a", "b"])
    da.isel(y=ind)  # same as da[:, ind]

    ind = xr.DataArray([["a", "b"], ["b", "a"]], dims=["a", "b"])
    da.loc[:, ind]  # same as da.sel(y=ind)

These methods may also be applied to ``Dataset`` objects

.. ipython:: python

    ds = da.to_dataset(name="bar")
    ds.isel(x=xr.DataArray([0, 1, 2], dims=["points"]))

Vectorized indexing may be used to extract information from the nearest
grid cells of interest, for example, the nearest climate model grid cells
to a collection specified weather station latitudes and longitudes.

.. ipython:: python

    ds = xr.tutorial.open_dataset("air_temperature")

    # Define target latitude and longitude (where weather stations might be)
    target_lon = xr.DataArray([200, 201, 202, 205], dims="points")
    target_lat = xr.DataArray([31, 41, 42, 42], dims="points")

    # Retrieve data at the grid cells nearest to the target latitudes and longitudes
    da = ds["air"].sel(lon=target_lon, lat=target_lat, method="nearest")
    da

.. tip::

  If you are lazily loading your data from disk, not every form of vectorized
  indexing is supported (or if supported, may not be supported efficiently).
  You may find increased performance by loading your data into memory first,
  e.g., with :py:meth:`~xarray.Dataset.load`.

.. note::

  If an indexer is a :py:meth:`~xarray.DataArray`, its coordinates should not
  conflict with the selected subpart of the target array (except for the
  explicitly indexed dimensions with ``.loc``/``.sel``).
  Otherwise, ``IndexError`` will be raised.


.. _assigning_values:

Assigning values with indexing
------------------------------

To select and assign values to a portion of a :py:meth:`~xarray.DataArray` you
can use indexing with ``.loc`` :

.. ipython:: python

    ds = xr.tutorial.open_dataset("air_temperature")

    # add an empty 2D dataarray
    ds["empty"] = xr.full_like(ds.air.mean("time"), fill_value=0)

    # modify one grid point using loc()
    ds["empty"].loc[dict(lon=260, lat=30)] = 100

    # modify a 2D region using loc()
    lc = ds.coords["lon"]
    la = ds.coords["lat"]
    ds["empty"].loc[
        dict(lon=lc[(lc > 220) & (lc < 260)], lat=la[(la > 20) & (la < 60)])
    ] = 100

or :py:meth:`~xarray.where`:

.. ipython:: python

    # modify one grid point using xr.where()
    ds["empty"] = xr.where(
        (ds.coords["lat"] == 20) & (ds.coords["lon"] == 260), 100, ds["empty"]
    )

    # or modify a 2D region using xr.where()
    mask = (
        (ds.coords["lat"] > 20)
        & (ds.coords["lat"] < 60)
        & (ds.coords["lon"] > 220)
        & (ds.coords["lon"] < 260)
    )
    ds["empty"] = xr.where(mask, 100, ds["empty"])



Vectorized indexing can also be used to assign values to xarray object.

.. ipython:: python

    da = xr.DataArray(
        np.arange(12).reshape((3, 4)),
        dims=["x", "y"],
        coords={"x": [0, 1, 2], "y": ["a", "b", "c", "d"]},
    )
    da
    da[0] = -1  # assignment with broadcasting
    da

    ind_x = xr.DataArray([0, 1], dims=["x"])
    ind_y = xr.DataArray([0, 1], dims=["y"])
    da[ind_x, ind_y] = -2  # assign -2 to (ix, iy) = (0, 0) and (1, 1)
    da

    da[ind_x, ind_y] += 100  # increment is also possible
    da

Like ``numpy.ndarray``, value assignment sometimes works differently from what one may expect.

.. ipython:: python

    da = xr.DataArray([0, 1, 2, 3], dims=["x"])
    ind = xr.DataArray([0, 0, 0], dims=["x"])
    da[ind] -= 1
    da

Where the 0th element will be subtracted 1 only once.
This is because ``v[0] = v[0] - 1`` is called three times, rather than
``v[0] = v[0] - 1 - 1 - 1``.
See `Assigning values to indexed arrays`__ for the details.

__ https://docs.scipy.org/doc/numpy/user/basics.indexing.html#assigning-values-to-indexed-arrays


.. note::
  Dask array does not support value assignment
  (see :ref:`dask` for the details).

.. note::

  Coordinates in both the left- and right-hand-side arrays should not
  conflict with each other.
  Otherwise, ``IndexError`` will be raised.

.. warning::

  Do not try to assign values when using any of the indexing methods ``isel``
  or ``sel``::

    # DO NOT do this
    da.isel(space=0) = 0

  Assigning values with the chained indexing using ``.sel`` or ``.isel`` fails silently.

  .. ipython:: python

      da = xr.DataArray([0, 1, 2, 3], dims=["x"])
      # DO NOT do this
      da.isel(x=[0, 1, 2])[1] = -1
      da

You can also assign values to all variables of a :py:class:`Dataset` at once:

.. ipython:: python

    ds_org = xr.tutorial.open_dataset("eraint_uvz").isel(
        latitude=slice(56, 59), longitude=slice(255, 258), level=0
    )
    # set all values to 0
    ds = xr.zeros_like(ds_org)
    ds

    # by integer
    ds[dict(latitude=2, longitude=2)] = 1
    ds["u"]
    ds["v"]

    # by label
    ds.loc[dict(latitude=47.25, longitude=[11.25, 12])] = 100
    ds["u"]

    # dataset as new values
    new_dat = ds_org.loc[dict(latitude=48, longitude=[11.25, 12])]
    new_dat
    ds.loc[dict(latitude=47.25, longitude=[11.25, 12])] = new_dat
    ds["u"]

The dimensions can differ between the variables in the dataset, but all variables need to have at least the dimensions specified in the indexer dictionary.
The new values must be either a scalar, a :py:class:`DataArray` or a :py:class:`Dataset` itself that contains all variables that also appear in the dataset to be modified.

.. _more_advanced_indexing:

More advanced indexing
-----------------------

The use of :py:meth:`~xarray.DataArray` objects as indexers enables very
flexible indexing. The following is an example of the pointwise indexing:

.. ipython:: python

    da = xr.DataArray(np.arange(56).reshape((7, 8)), dims=["x", "y"])
    da
    da.isel(x=xr.DataArray([0, 1, 6], dims="z"), y=xr.DataArray([0, 1, 0], dims="z"))


where three elements at ``(ix, iy) = ((0, 0), (1, 1), (6, 0))`` are selected
and mapped along a new dimension ``z``.

If you want to add a coordinate to the new dimension ``z``,
you can supply a :py:class:`~xarray.DataArray` with a coordinate,

.. ipython:: python

    da.isel(
        x=xr.DataArray([0, 1, 6], dims="z", coords={"z": ["a", "b", "c"]}),
        y=xr.DataArray([0, 1, 0], dims="z"),
    )

Analogously, label-based pointwise-indexing is also possible by the ``.sel``
method:

.. ipython:: python

    da = xr.DataArray(
        np.random.rand(4, 3),
        [
            ("time", pd.date_range("2000-01-01", periods=4)),
            ("space", ["IA", "IL", "IN"]),
        ],
    )
    times = xr.DataArray(
        pd.to_datetime(["2000-01-03", "2000-01-02", "2000-01-01"]), dims="new_time"
    )
    da.sel(space=xr.DataArray(["IA", "IL", "IN"], dims=["new_time"]), time=times)

.. _align and reindex:

Align and reindex
-----------------

Xarray's ``reindex``, ``reindex_like`` and ``align`` impose a ``DataArray`` or
``Dataset`` onto a new set of coordinates corresponding to dimensions. The
original values are subset to the index labels still found in the new labels,
and values corresponding to new labels not found in the original object are
in-filled with `NaN`.

Xarray operations that combine multiple objects generally automatically align
their arguments to share the same indexes. However, manual alignment can be
useful for greater control and for increased performance.

To reindex a particular dimension, use :py:meth:`~xarray.DataArray.reindex`:

.. ipython:: python

    da.reindex(space=["IA", "CA"])

The :py:meth:`~xarray.DataArray.reindex_like` method is a useful shortcut.
To demonstrate, we will make a subset DataArray with new values:

.. ipython:: python

    foo = da.rename("foo")
    baz = (10 * da[:2, :2]).rename("baz")
    baz

Reindexing ``foo`` with ``baz`` selects out the first two values along each
dimension:

.. ipython:: python

    foo.reindex_like(baz)

The opposite operation asks us to reindex to a larger shape, so we fill in
the missing values with `NaN`:

.. ipython:: python

    baz.reindex_like(foo)

The :py:func:`~xarray.align` function lets us perform more flexible database-like
``'inner'``, ``'outer'``, ``'left'`` and ``'right'`` joins:

.. ipython:: python

    xr.align(foo, baz, join="inner")
    xr.align(foo, baz, join="outer")

Both ``reindex_like`` and ``align`` work interchangeably between
:py:class:`~xarray.DataArray` and :py:class:`~xarray.Dataset` objects, and with any number of matching dimension names:

.. ipython:: python

    ds
    ds.reindex_like(baz)
    other = xr.DataArray(["a", "b", "c"], dims="other")
    # this is a no-op, because there are no shared dimension names
    ds.reindex_like(other)

.. _indexing.missing_coordinates:

Missing coordinate labels
-------------------------

Coordinate labels for each dimension are optional (as of xarray v0.9). Label
based indexing with ``.sel`` and ``.loc`` uses standard positional,
integer-based indexing as a fallback for dimensions without a coordinate label:

.. ipython:: python

    da = xr.DataArray([1, 2, 3], dims="x")
    da.sel(x=[0, -1])

Alignment between xarray objects where one or both do not have coordinate labels
succeeds only if all dimensions of the same name have the same length.
Otherwise, it raises an informative error:

.. ipython::
    :verbatim:

    In [62]: xr.align(da, da[:2])
    ValueError: arguments without labels along dimension 'x' cannot be aligned because they have different dimension sizes: {2, 3}

Underlying Indexes
------------------

Xarray uses the :py:class:`pandas.Index` internally to perform indexing
operations.  If you need to access the underlying indexes, they are available
through the :py:attr:`~xarray.DataArray.indexes` attribute.

.. ipython:: python

    da = xr.DataArray(
        np.random.rand(4, 3),
        [
            ("time", pd.date_range("2000-01-01", periods=4)),
            ("space", ["IA", "IL", "IN"]),
        ],
    )
    da
    da.indexes
    da.indexes["time"]

Use :py:meth:`~xarray.DataArray.get_index` to get an index for a dimension,
falling back to a default :py:class:`pandas.RangeIndex` if it has no coordinate
labels:

.. ipython:: python

    da = xr.DataArray([1, 2, 3], dims="x")
    da
    da.get_index("x")


.. _copies_vs_views:

Copies vs. Views
----------------

Whether array indexing returns a view or a copy of the underlying
data depends on the nature of the labels.

For positional (integer)
indexing, xarray follows the same rules as NumPy:

* Positional indexing with only integers and slices returns a view.
* Positional indexing with arrays or lists returns a copy.

The rules for label based indexing are more complex:

* Label-based indexing with only slices returns a view.
* Label-based indexing with arrays returns a copy.
* Label-based indexing with scalars returns a view or a copy, depending
  upon if the corresponding positional indexer can be represented as an
  integer or a slice object. The exact rules are determined by pandas.

Whether data is a copy or a view is more predictable in xarray than in pandas, so
unlike pandas, xarray does not produce `SettingWithCopy warnings`_. However, you
should still avoid assignment with chained indexing.

.. _SettingWithCopy warnings: http://pandas.pydata.org/pandas-docs/stable/indexing.html#returning-a-view-versus-a-copy


.. _multi-level indexing:

Multi-level indexing
--------------------

Just like pandas, advanced indexing on multi-level indexes is possible with
``loc`` and ``sel``. You can slice a multi-index by providing multiple indexers,
i.e., a tuple of slices, labels, list of labels, or any selector allowed by
pandas:

.. ipython:: python

    midx = pd.MultiIndex.from_product([list("abc"), [0, 1]], names=("one", "two"))
    mda = xr.DataArray(np.random.rand(6, 3), [("x", midx), ("y", range(3))])
    mda
    mda.sel(x=(list("ab"), [0]))

You can also select multiple elements by providing a list of labels or tuples or
a slice of tuples:

.. ipython:: python

    mda.sel(x=[("a", 0), ("b", 1)])

Additionally, xarray supports dictionaries:

.. ipython:: python

    mda.sel(x={"one": "a", "two": 0})

For convenience, ``sel`` also accepts multi-index levels directly
as keyword arguments:

.. ipython:: python

    mda.sel(one="a", two=0)

Note that using ``sel`` it is not possible to mix a dimension
indexer with level indexers for that dimension
(e.g., ``mda.sel(x={'one': 'a'}, two=0)`` will raise a ``ValueError``).

Like pandas, xarray handles partial selection on multi-index (level drop).
As shown below, it also renames the dimension / coordinate when the
multi-index is reduced to a single index.

.. ipython:: python

    mda.loc[{"one": "a"}, ...]

Unlike pandas, xarray does not guess whether you provide index levels or
dimensions when using ``loc`` in some ambiguous cases. For example, for
``mda.loc[{'one': 'a', 'two': 0}]`` and ``mda.loc['a', 0]`` xarray
always interprets ('one', 'two') and ('a', 0) as the names and
labels of the 1st and 2nd dimension, respectively. You must specify all
dimensions or use the ellipsis in the ``loc`` specifier, e.g. in the example
above, ``mda.loc[{'one': 'a', 'two': 0}, :]`` or ``mda.loc[('a', 0), ...]``.


.. _indexing.rules:

Indexing rules
--------------

Here we describe the full rules xarray uses for vectorized indexing. Note that
this is for the purposes of explanation: for the sake of efficiency and to
support various backends, the actual implementation is different.

0. (Only for label based indexing.) Look up positional indexes along each
   dimension from the corresponding :py:class:`pandas.Index`.

1. A full slice object ``:`` is inserted for each dimension without an indexer.

2. ``slice`` objects are converted into arrays, given by
   ``np.arange(*slice.indices(...))``.

3. Assume dimension names for array indexers without dimensions, such as
   ``np.ndarray`` and ``list``, from the dimensions to be indexed along.
   For example, ``v.isel(x=[0, 1])`` is understood as
   ``v.isel(x=xr.DataArray([0, 1], dims=['x']))``.

4. For each variable in a ``Dataset`` or  ``DataArray`` (the array and its
   coordinates):

   a. Broadcast all relevant indexers based on their dimension names
      (see :ref:`compute.broadcasting` for full details).

   b. Index the underling array by the broadcast indexers, using NumPy's
      advanced indexing rules.

5. If any indexer DataArray has coordinates and no coordinate with the
   same name exists, attach them to the indexed object.

.. note::

  Only 1-dimensional boolean arrays can be used as indexers.
.. _combining data:

Combining data
--------------

.. ipython:: python
    :suppress:

    import numpy as np
    import pandas as pd
    import xarray as xr

    np.random.seed(123456)

* For combining datasets or data arrays along a single dimension, see concatenate_.
* For combining datasets with different variables, see merge_.
* For combining datasets or data arrays with different indexes or missing values, see combine_.
* For combining datasets or data arrays along multiple dimensions see combining.multi_.

.. _concatenate:

Concatenate
~~~~~~~~~~~

To combine arrays along existing or new dimension into a larger array, you
can use :py:func:`~xarray.concat`. ``concat`` takes an iterable of ``DataArray``
or ``Dataset`` objects, as well as a dimension name, and concatenates along
that dimension:

.. ipython:: python

    da = xr.DataArray(
        np.arange(6).reshape(2, 3), [("x", ["a", "b"]), ("y", [10, 20, 30])]
    )
    da.isel(y=slice(0, 1))  # same as da[:, :1]
    # This resembles how you would use np.concatenate:
    xr.concat([da[:, :1], da[:, 1:]], dim="y")
    # For more friendly pandas-like indexing you can use:
    xr.concat([da.isel(y=slice(0, 1)), da.isel(y=slice(1, None))], dim="y")

In addition to combining along an existing dimension, ``concat`` can create a
new dimension by stacking lower dimensional arrays together:

.. ipython:: python

    da.sel(x="a")
    xr.concat([da.isel(x=0), da.isel(x=1)], "x")

If the second argument to ``concat`` is a new dimension name, the arrays will
be concatenated along that new dimension, which is always inserted as the first
dimension:

.. ipython:: python

    xr.concat([da.isel(x=0), da.isel(x=1)], "new_dim")

The second argument to ``concat`` can also be an :py:class:`~pandas.Index` or
:py:class:`~xarray.DataArray` object as well as a string, in which case it is
used to label the values along the new dimension:

.. ipython:: python

    xr.concat([da.isel(x=0), da.isel(x=1)], pd.Index([-90, -100], name="new_dim"))

Of course, ``concat`` also works on ``Dataset`` objects:

.. ipython:: python

    ds = da.to_dataset(name="foo")
    xr.concat([ds.sel(x="a"), ds.sel(x="b")], "x")

:py:func:`~xarray.concat` has a number of options which provide deeper control
over which variables are concatenated and how it handles conflicting variables
between datasets. With the default parameters, xarray will load some coordinate
variables into memory to compare them between datasets. This may be prohibitively
expensive if you are manipulating your dataset lazily using :ref:`dask`.

.. _merge:

Merge
~~~~~

To combine variables and coordinates between multiple ``DataArray`` and/or
``Dataset`` objects, use :py:func:`~xarray.merge`. It can merge a list of
``Dataset``, ``DataArray`` or dictionaries of objects convertible to
``DataArray`` objects:

.. ipython:: python

    xr.merge([ds, ds.rename({"foo": "bar"})])
    xr.merge([xr.DataArray(n, name="var%d" % n) for n in range(5)])

If you merge another dataset (or a dictionary including data array objects), by
default the resulting dataset will be aligned on the **union** of all index
coordinates:

.. ipython:: python

    other = xr.Dataset({"bar": ("x", [1, 2, 3, 4]), "x": list("abcd")})
    xr.merge([ds, other])

This ensures that ``merge`` is non-destructive. ``xarray.MergeError`` is raised
if you attempt to merge two variables with the same name but different values:

.. ipython::

    @verbatim
    In [1]: xr.merge([ds, ds + 1])
    MergeError: conflicting values for variable 'foo' on objects to be combined:
    first value: <xarray.Variable (x: 2, y: 3)>
    array([[ 0.4691123 , -0.28286334, -1.5090585 ],
           [-1.13563237,  1.21211203, -0.17321465]])
    second value: <xarray.Variable (x: 2, y: 3)>
    array([[ 1.4691123 ,  0.71713666, -0.5090585 ],
           [-0.13563237,  2.21211203,  0.82678535]])

The same non-destructive merging between ``DataArray`` index coordinates is
used in the :py:class:`~xarray.Dataset` constructor:

.. ipython:: python

    xr.Dataset({"a": da.isel(x=slice(0, 1)), "b": da.isel(x=slice(1, 2))})

.. _combine:

Combine
~~~~~~~

The instance method :py:meth:`~xarray.DataArray.combine_first` combines two
datasets/data arrays and defaults to non-null values in the calling object,
using values from the called object to fill holes.  The resulting coordinates
are the union of coordinate labels. Vacant cells as a result of the outer-join
are filled with ``NaN``. For example:

.. ipython:: python

    ar0 = xr.DataArray([[0, 0], [0, 0]], [("x", ["a", "b"]), ("y", [-1, 0])])
    ar1 = xr.DataArray([[1, 1], [1, 1]], [("x", ["b", "c"]), ("y", [0, 1])])
    ar0.combine_first(ar1)
    ar1.combine_first(ar0)

For datasets, ``ds0.combine_first(ds1)`` works similarly to
``xr.merge([ds0, ds1])``, except that ``xr.merge`` raises ``MergeError`` when
there are conflicting values in variables to be merged, whereas
``.combine_first`` defaults to the calling object's values.

.. _update:

Update
~~~~~~

In contrast to ``merge``, :py:meth:`~xarray.Dataset.update` modifies a dataset
in-place without checking for conflicts, and will overwrite any existing
variables with new values:

.. ipython:: python

    ds.update({"space": ("space", [10.2, 9.4, 3.9])})

However, dimensions are still required to be consistent between different
Dataset variables, so you cannot change the size of a dimension unless you
replace all dataset variables that use it.

``update`` also performs automatic alignment if necessary. Unlike ``merge``, it
maintains the alignment of the original array instead of merging indexes:

.. ipython:: python

    ds.update(other)

The exact same alignment logic when setting a variable with ``__setitem__``
syntax:

.. ipython:: python

    ds["baz"] = xr.DataArray([9, 9, 9, 9, 9], coords=[("x", list("abcde"))])
    ds.baz

Equals and identical
~~~~~~~~~~~~~~~~~~~~

Xarray objects can be compared by using the :py:meth:`~xarray.Dataset.equals`,
:py:meth:`~xarray.Dataset.identical` and
:py:meth:`~xarray.Dataset.broadcast_equals` methods. These methods are used by
the optional ``compat`` argument on ``concat`` and ``merge``.

:py:attr:`~xarray.Dataset.equals` checks dimension names, indexes and array
values:

.. ipython:: python

    da.equals(da.copy())

:py:attr:`~xarray.Dataset.identical` also checks attributes, and the name of each
object:

.. ipython:: python

    da.identical(da.rename("bar"))

:py:attr:`~xarray.Dataset.broadcast_equals` does a more relaxed form of equality
check that allows variables to have different dimensions, as long as values
are constant along those new dimensions:

.. ipython:: python

    left = xr.Dataset(coords={"x": 0})
    right = xr.Dataset({"x": [0, 0, 0]})
    left.broadcast_equals(right)

Like pandas objects, two xarray objects are still equal or identical if they have
missing values marked by ``NaN`` in the same locations.

In contrast, the ``==`` operation performs element-wise comparison (like
numpy):

.. ipython:: python

    da == da.copy()

Note that ``NaN`` does not compare equal to ``NaN`` in element-wise comparison;
you may need to deal with missing values explicitly.

.. _combining.no_conflicts:

Merging with 'no_conflicts'
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``compat`` argument ``'no_conflicts'`` is only available when
combining xarray objects with ``merge``. In addition to the above comparison
methods it allows the merging of xarray objects with locations where *either*
have ``NaN`` values. This can be used to combine data with overlapping
coordinates as long as any non-missing values agree or are disjoint:

.. ipython:: python

    ds1 = xr.Dataset({"a": ("x", [10, 20, 30, np.nan])}, {"x": [1, 2, 3, 4]})
    ds2 = xr.Dataset({"a": ("x", [np.nan, 30, 40, 50])}, {"x": [2, 3, 4, 5]})
    xr.merge([ds1, ds2], compat="no_conflicts")

Note that due to the underlying representation of missing values as floating
point numbers (``NaN``), variable data type is not always preserved when merging
in this manner.

.. _combining.multi:

Combining along multiple dimensions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For combining many objects along multiple dimensions xarray provides
:py:func:`~xarray.combine_nested` and :py:func:`~xarray.combine_by_coords`. These
functions use a combination of ``concat`` and ``merge`` across different
variables to combine many objects into one.

:py:func:`~xarray.combine_nested` requires specifying the order in which the
objects should be combined, while :py:func:`~xarray.combine_by_coords` attempts to
infer this ordering automatically from the coordinates in the data.

:py:func:`~xarray.combine_nested` is useful when you know the spatial
relationship between each object in advance. The datasets must be provided in
the form of a nested list, which specifies their relative position and
ordering. A common task is collecting data from a parallelized simulation where
each processor wrote out data to a separate file. A domain which was decomposed
into 4 parts, 2 each along both the x and y axes, requires organising the
datasets into a doubly-nested list, e.g:

.. ipython:: python

    arr = xr.DataArray(
        name="temperature", data=np.random.randint(5, size=(2, 2)), dims=["x", "y"]
    )
    arr
    ds_grid = [[arr, arr], [arr, arr]]
    xr.combine_nested(ds_grid, concat_dim=["x", "y"])

:py:func:`~xarray.combine_nested` can also be used to explicitly merge datasets
with different variables. For example if we have 4 datasets, which are divided
along two times, and contain two different variables, we can pass ``None``
to ``'concat_dim'`` to specify the dimension of the nested list over which
we wish to use ``merge`` instead of ``concat``:

.. ipython:: python

    temp = xr.DataArray(name="temperature", data=np.random.randn(2), dims=["t"])
    precip = xr.DataArray(name="precipitation", data=np.random.randn(2), dims=["t"])
    ds_grid = [[temp, precip], [temp, precip]]
    xr.combine_nested(ds_grid, concat_dim=["t", None])

:py:func:`~xarray.combine_by_coords` is for combining objects which have dimension
coordinates which specify their relationship to and order relative to one
another, for example a linearly-increasing 'time' dimension coordinate.

Here we combine two datasets using their common dimension coordinates. Notice
they are concatenated in order based on the values in their dimension
coordinates, not on their position in the list passed to ``combine_by_coords``.

.. ipython:: python
    :okwarning:

    x1 = xr.DataArray(name="foo", data=np.random.randn(3), coords=[("x", [0, 1, 2])])
    x2 = xr.DataArray(name="foo", data=np.random.randn(3), coords=[("x", [3, 4, 5])])
    xr.combine_by_coords([x2, x1])

These functions can be used by :py:func:`~xarray.open_mfdataset` to open many
files as one dataset. The particular function used is specified by setting the
argument ``'combine'`` to ``'by_coords'`` or ``'nested'``. This is useful for
situations where your data is split across many files in multiple locations,
which have some known relationship between one another.
.. _groupby:

GroupBy: split-apply-combine
----------------------------

Xarray supports `"group by"`__ operations with the same API as pandas to
implement the `split-apply-combine`__ strategy:

__ http://pandas.pydata.org/pandas-docs/stable/groupby.html
__ http://www.jstatsoft.org/v40/i01/paper

- Split your data into multiple independent groups.
- Apply some function to each group.
- Combine your groups back into a single data object.

Group by operations work on both :py:class:`~xarray.Dataset` and
:py:class:`~xarray.DataArray` objects. Most of the examples focus on grouping by
a single one-dimensional variable, although support for grouping
over a multi-dimensional variable has recently been implemented. Note that for
one-dimensional data, it is usually faster to rely on pandas' implementation of
the same pipeline.

Split
~~~~~

Let's create a simple example dataset:

.. ipython:: python
    :suppress:

    import numpy as np
    import pandas as pd
    import xarray as xr

    np.random.seed(123456)

.. ipython:: python

    ds = xr.Dataset(
        {"foo": (("x", "y"), np.random.rand(4, 3))},
        coords={"x": [10, 20, 30, 40], "letters": ("x", list("abba"))},
    )
    arr = ds["foo"]
    ds

If we groupby the name of a variable or coordinate in a dataset (we can also
use a DataArray directly), we get back a ``GroupBy`` object:

.. ipython:: python

    ds.groupby("letters")

This object works very similarly to a pandas GroupBy object. You can view
the group indices with the ``groups`` attribute:

.. ipython:: python

    ds.groupby("letters").groups

You can also iterate over groups in ``(label, group)`` pairs:

.. ipython:: python

    list(ds.groupby("letters"))

You can index out a particular group:

.. ipython:: python

    ds.groupby("letters")["b"]

Just like in pandas, creating a GroupBy object is cheap: it does not actually
split the data until you access particular values.

Binning
~~~~~~~

Sometimes you don't want to use all the unique values to determine the groups
but instead want to "bin" the data into coarser groups. You could always create
a customized coordinate, but xarray facilitates this via the
:py:meth:`~xarray.Dataset.groupby_bins` method.

.. ipython:: python

    x_bins = [0, 25, 50]
    ds.groupby_bins("x", x_bins).groups

The binning is implemented via :func:`pandas.cut`, whose documentation details how
the bins are assigned. As seen in the example above, by default, the bins are
labeled with strings using set notation to precisely identify the bin limits. To
override this behavior, you can specify the bin labels explicitly. Here we
choose `float` labels which identify the bin centers:

.. ipython:: python

    x_bin_labels = [12.5, 37.5]
    ds.groupby_bins("x", x_bins, labels=x_bin_labels).groups


Apply
~~~~~

To apply a function to each group, you can use the flexible
:py:meth:`~xarray.core.groupby.DatasetGroupBy.map` method. The resulting objects are automatically
concatenated back together along the group axis:

.. ipython:: python

    def standardize(x):
        return (x - x.mean()) / x.std()


    arr.groupby("letters").map(standardize)

GroupBy objects also have a :py:meth:`~xarray.core.groupby.DatasetGroupBy.reduce` method and
methods like :py:meth:`~xarray.core.groupby.DatasetGroupBy.mean` as shortcuts for applying an
aggregation function:

.. ipython:: python

    arr.groupby("letters").mean(dim="x")

Using a groupby is thus also a convenient shortcut for aggregating over all
dimensions *other than* the provided one:

.. ipython:: python

    ds.groupby("x").std(...)

.. note::

    We use an ellipsis (`...`) here to indicate we want to reduce over all
    other dimensions


First and last
~~~~~~~~~~~~~~

There are two special aggregation operations that are currently only found on
groupby objects: first and last. These provide the first or last example of
values for group along the grouped dimension:

.. ipython:: python

    ds.groupby("letters").first(...)

By default, they skip missing values (control this with ``skipna``).

Grouped arithmetic
~~~~~~~~~~~~~~~~~~

GroupBy objects also support a limited set of binary arithmetic operations, as
a shortcut for mapping over all unique labels. Binary arithmetic is supported
for ``(GroupBy, Dataset)`` and ``(GroupBy, DataArray)`` pairs, as long as the
dataset or data array uses the unique grouped values as one of its index
coordinates. For example:

.. ipython:: python

    alt = arr.groupby("letters").mean(...)
    alt
    ds.groupby("letters") - alt

This last line is roughly equivalent to the following::

    results = []
    for label, group in ds.groupby('letters'):
        results.append(group - alt.sel(letters=label))
    xr.concat(results, dim='x')

Squeezing
~~~~~~~~~

When grouping over a dimension, you can control whether the dimension is
squeezed out or if it should remain with length one on each group by using
the ``squeeze`` parameter:

.. ipython:: python

    next(iter(arr.groupby("x")))

.. ipython:: python

    next(iter(arr.groupby("x", squeeze=False)))

Although xarray will attempt to automatically
:py:attr:`~xarray.DataArray.transpose` dimensions back into their original order
when you use apply, it is sometimes useful to set ``squeeze=False`` to
guarantee that all original dimensions remain unchanged.

You can always squeeze explicitly later with the Dataset or DataArray
:py:meth:`~xarray.DataArray.squeeze` methods.

.. _groupby.multidim:

Multidimensional Grouping
~~~~~~~~~~~~~~~~~~~~~~~~~

Many datasets have a multidimensional coordinate variable (e.g. longitude)
which is different from the logical grid dimensions (e.g. nx, ny). Such
variables are valid under the `CF conventions`__. Xarray supports groupby
operations over multidimensional coordinate variables:

__ http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#_two_dimensional_latitude_longitude_coordinate_variables

.. ipython:: python

    da = xr.DataArray(
        [[0, 1], [2, 3]],
        coords={
            "lon": (["ny", "nx"], [[30, 40], [40, 50]]),
            "lat": (["ny", "nx"], [[10, 10], [20, 20]]),
        },
        dims=["ny", "nx"],
    )
    da
    da.groupby("lon").sum(...)
    da.groupby("lon").map(lambda x: x - x.mean(), shortcut=False)

Because multidimensional groups have the ability to generate a very large
number of bins, coarse-binning via :py:meth:`~xarray.Dataset.groupby_bins`
may be desirable:

.. ipython:: python

    da.groupby_bins("lon", [0, 45, 50]).sum()

These methods group by `lon` values. It is also possible to groupby each
cell in a grid, regardless of value, by stacking multiple dimensions,
applying your function, and then unstacking the result:

.. ipython:: python

    stacked = da.stack(gridcell=["ny", "nx"])
    stacked.groupby("gridcell").sum(...).unstack("gridcell")
.. _data structures:

Data Structures
===============

.. ipython:: python
    :suppress:

    import numpy as np
    import pandas as pd
    import xarray as xr

    np.random.seed(123456)
    np.set_printoptions(threshold=10)

DataArray
---------

:py:class:`xarray.DataArray` is xarray's implementation of a labeled,
multi-dimensional array. It has several key properties:

- ``values``: a :py:class:`numpy.ndarray` holding the array's values
- ``dims``: dimension names for each axis (e.g., ``('x', 'y', 'z')``)
- ``coords``: a dict-like container of arrays (*coordinates*) that label each
  point (e.g., 1-dimensional arrays of numbers, datetime objects or
  strings)
- ``attrs``: :py:class:`dict` to hold arbitrary metadata (*attributes*)

Xarray uses ``dims`` and ``coords`` to enable its core metadata aware operations.
Dimensions provide names that xarray uses instead of the ``axis`` argument found
in many numpy functions. Coordinates enable fast label based indexing and
alignment, building on the functionality of the ``index`` found on a pandas
:py:class:`~pandas.DataFrame` or :py:class:`~pandas.Series`.

DataArray objects also can have a ``name`` and can hold arbitrary metadata in
the form of their ``attrs`` property. Names and attributes are strictly for
users and user-written code: xarray makes no attempt to interpret them, and
propagates them only in unambiguous cases
(see FAQ, :ref:`approach to metadata`).

.. _creating a dataarray:

Creating a DataArray
~~~~~~~~~~~~~~~~~~~~

The :py:class:`~xarray.DataArray` constructor takes:

- ``data``: a multi-dimensional array of values (e.g., a numpy ndarray,
  :py:class:`~pandas.Series`, :py:class:`~pandas.DataFrame` or ``pandas.Panel``)
- ``coords``: a list or dictionary of coordinates. If a list, it should be a
  list of tuples where the first element is the dimension name and the second
  element is the corresponding coordinate array_like object.
- ``dims``: a list of dimension names. If omitted and ``coords`` is a list of
  tuples, dimension names are taken from ``coords``.
- ``attrs``: a dictionary of attributes to add to the instance
- ``name``: a string that names the instance

.. ipython:: python

    data = np.random.rand(4, 3)
    locs = ["IA", "IL", "IN"]
    times = pd.date_range("2000-01-01", periods=4)
    foo = xr.DataArray(data, coords=[times, locs], dims=["time", "space"])
    foo

Only ``data`` is required; all of other arguments will be filled
in with default values:

.. ipython:: python

    xr.DataArray(data)

As you can see, dimension names are always present in the xarray data model: if
you do not provide them, defaults of the form ``dim_N`` will be created.
However, coordinates are always optional, and dimensions do not have automatic
coordinate labels.

.. note::

  This is different from pandas, where axes always have tick labels, which
  default to the integers ``[0, ..., n-1]``.

  Prior to xarray v0.9, xarray copied this behavior: default coordinates for
  each dimension would be created if coordinates were not supplied explicitly.
  This is no longer the case.

Coordinates can be specified in the following ways:

- A list of values with length equal to the number of dimensions, providing
  coordinate labels for each dimension. Each value must be of one of the
  following forms:

  * A :py:class:`~xarray.DataArray` or :py:class:`~xarray.Variable`
  * A tuple of the form ``(dims, data[, attrs])``, which is converted into
    arguments for :py:class:`~xarray.Variable`
  * A pandas object or scalar value, which is converted into a ``DataArray``
  * A 1D array or list, which is interpreted as values for a one dimensional
    coordinate variable along the same dimension as it's name

- A dictionary of ``{coord_name: coord}`` where values are of the same form
  as the list. Supplying coordinates as a dictionary allows other coordinates
  than those corresponding to dimensions (more on these later). If you supply
  ``coords`` as a dictionary, you must explicitly provide ``dims``.

As a list of tuples:

.. ipython:: python

    xr.DataArray(data, coords=[("time", times), ("space", locs)])

As a dictionary:

.. ipython:: python

    xr.DataArray(
        data,
        coords={
            "time": times,
            "space": locs,
            "const": 42,
            "ranking": ("space", [1, 2, 3]),
        },
        dims=["time", "space"],
    )

As a dictionary with coords across multiple dimensions:

.. ipython:: python

    xr.DataArray(
        data,
        coords={
            "time": times,
            "space": locs,
            "const": 42,
            "ranking": (("time", "space"), np.arange(12).reshape(4, 3)),
        },
        dims=["time", "space"],
    )

If you create a ``DataArray`` by supplying a pandas
:py:class:`~pandas.Series`, :py:class:`~pandas.DataFrame` or
``pandas.Panel``, any non-specified arguments in the
``DataArray`` constructor will be filled in from the pandas object:

.. ipython:: python

    df = pd.DataFrame({"x": [0, 1], "y": [2, 3]}, index=["a", "b"])
    df.index.name = "abc"
    df.columns.name = "xyz"
    df
    xr.DataArray(df)

DataArray properties
~~~~~~~~~~~~~~~~~~~~

Let's take a look at the important properties on our array:

.. ipython:: python

    foo.values
    foo.dims
    foo.coords
    foo.attrs
    print(foo.name)

You can modify ``values`` inplace:

.. ipython:: python

    foo.values = 1.0 * foo.values

.. note::

    The array values in a :py:class:`~xarray.DataArray` have a single
    (homogeneous) data type. To work with heterogeneous or structured data
    types in xarray, use coordinates, or put separate ``DataArray`` objects
    in a single :py:class:`~xarray.Dataset` (see below).

Now fill in some of that missing metadata:

.. ipython:: python

    foo.name = "foo"
    foo.attrs["units"] = "meters"
    foo

The :py:meth:`~xarray.DataArray.rename` method is another option, returning a
new data array:

.. ipython:: python

    foo.rename("bar")

DataArray Coordinates
~~~~~~~~~~~~~~~~~~~~~

The ``coords`` property is ``dict`` like. Individual coordinates can be
accessed from the coordinates by name, or even by indexing the data array
itself:

.. ipython:: python

    foo.coords["time"]
    foo["time"]

These are also :py:class:`~xarray.DataArray` objects, which contain tick-labels
for each dimension.

Coordinates can also be set or removed by using the dictionary like syntax:

.. ipython:: python

    foo["ranking"] = ("space", [1, 2, 3])
    foo.coords
    del foo["ranking"]
    foo.coords

For more details, see :ref:`coordinates` below.

Dataset
-------

:py:class:`xarray.Dataset` is xarray's multi-dimensional equivalent of a
:py:class:`~pandas.DataFrame`. It is a dict-like
container of labeled arrays (:py:class:`~xarray.DataArray` objects) with aligned
dimensions. It is designed as an in-memory representation of the data model
from the `netCDF`__ file format.

__ http://www.unidata.ucar.edu/software/netcdf/

In addition to the dict-like interface of the dataset itself, which can be used
to access any variable in a dataset, datasets have four key properties:

- ``dims``: a dictionary mapping from dimension names to the fixed length of
  each dimension (e.g., ``{'x': 6, 'y': 6, 'time': 8}``)
- ``data_vars``: a dict-like container of DataArrays corresponding to variables
- ``coords``: another dict-like container of DataArrays intended to label points
  used in ``data_vars`` (e.g., arrays of numbers, datetime objects or strings)
- ``attrs``: :py:class:`dict` to hold arbitrary metadata

The distinction between whether a variable falls in data or coordinates
(borrowed from `CF conventions`_) is mostly semantic, and you can probably get
away with ignoring it if you like: dictionary like access on a dataset will
supply variables found in either category. However, xarray does make use of the
distinction for indexing and computations. Coordinates indicate
constant/fixed/independent quantities, unlike the varying/measured/dependent
quantities that belong in data.

.. _CF conventions: http://cfconventions.org/

Here is an example of how we might structure a dataset for a weather forecast:

.. image:: ../_static/dataset-diagram.png

In this example, it would be natural to call ``temperature`` and
``precipitation`` "data variables" and all the other arrays "coordinate
variables" because they label the points along the dimensions. (see [1]_ for
more background on this example).

.. _dataarray constructor:

Creating a Dataset
~~~~~~~~~~~~~~~~~~

To make an :py:class:`~xarray.Dataset` from scratch, supply dictionaries for any
variables (``data_vars``), coordinates (``coords``) and attributes (``attrs``).

- ``data_vars`` should be a dictionary with each key as the name of the variable
  and each value as one of:

  * A :py:class:`~xarray.DataArray` or :py:class:`~xarray.Variable`
  * A tuple of the form ``(dims, data[, attrs])``, which is converted into
    arguments for :py:class:`~xarray.Variable`
  * A pandas object, which is converted into a ``DataArray``
  * A 1D array or list, which is interpreted as values for a one dimensional
    coordinate variable along the same dimension as it's name

- ``coords`` should be a dictionary of the same form as ``data_vars``.

- ``attrs`` should be a dictionary.

Let's create some fake data for the example we show above:

.. ipython:: python

    temp = 15 + 8 * np.random.randn(2, 2, 3)
    precip = 10 * np.random.rand(2, 2, 3)
    lon = [[-99.83, -99.32], [-99.79, -99.23]]
    lat = [[42.25, 42.21], [42.63, 42.59]]

    # for real use cases, its good practice to supply array attributes such as
    # units, but we won't bother here for the sake of brevity
    ds = xr.Dataset(
        {
            "temperature": (["x", "y", "time"], temp),
            "precipitation": (["x", "y", "time"], precip),
        },
        coords={
            "lon": (["x", "y"], lon),
            "lat": (["x", "y"], lat),
            "time": pd.date_range("2014-09-06", periods=3),
            "reference_time": pd.Timestamp("2014-09-05"),
        },
    )
    ds

Here we pass :py:class:`xarray.DataArray` objects or a pandas object as values
in the dictionary:

.. ipython:: python

    xr.Dataset(dict(bar=foo))


.. ipython:: python

    xr.Dataset(dict(bar=foo.to_pandas()))

Where a pandas object is supplied as a value, the names of its indexes are used as dimension
names, and its data is aligned to any existing dimensions.

You can also create an dataset from:

- A :py:class:`pandas.DataFrame` or ``pandas.Panel`` along its columns and items
  respectively, by passing it into the :py:class:`~xarray.Dataset` directly
- A :py:class:`pandas.DataFrame` with :py:meth:`Dataset.from_dataframe <xarray.Dataset.from_dataframe>`,
  which will additionally handle MultiIndexes See :ref:`pandas`
- A netCDF file on disk with :py:func:`~xarray.open_dataset`. See :ref:`io`.

Dataset contents
~~~~~~~~~~~~~~~~

:py:class:`~xarray.Dataset` implements the Python mapping interface, with
values given by :py:class:`xarray.DataArray` objects:

.. ipython:: python

    "temperature" in ds
    ds["temperature"]

Valid keys include each listed coordinate and data variable.

Data and coordinate variables are also contained separately in the
:py:attr:`~xarray.Dataset.data_vars` and :py:attr:`~xarray.Dataset.coords`
dictionary-like attributes:

.. ipython:: python

    ds.data_vars
    ds.coords

Finally, like data arrays, datasets also store arbitrary metadata in the form
of `attributes`:

.. ipython:: python

    ds.attrs

    ds.attrs["title"] = "example attribute"
    ds

Xarray does not enforce any restrictions on attributes, but serialization to
some file formats may fail if you use objects that are not strings, numbers
or :py:class:`numpy.ndarray` objects.

As a useful shortcut, you can use attribute style access for reading (but not
setting) variables and attributes:

.. ipython:: python

    ds.temperature

This is particularly useful in an exploratory context, because you can
tab-complete these variable names with tools like IPython.

.. _dictionary_like_methods:

Dictionary like methods
~~~~~~~~~~~~~~~~~~~~~~~

We can update a dataset in-place using Python's standard dictionary syntax. For
example, to create this example dataset from scratch, we could have written:

.. ipython:: python

    ds = xr.Dataset()
    ds["temperature"] = (("x", "y", "time"), temp)
    ds["temperature_double"] = (("x", "y", "time"), temp * 2)
    ds["precipitation"] = (("x", "y", "time"), precip)
    ds.coords["lat"] = (("x", "y"), lat)
    ds.coords["lon"] = (("x", "y"), lon)
    ds.coords["time"] = pd.date_range("2014-09-06", periods=3)
    ds.coords["reference_time"] = pd.Timestamp("2014-09-05")

To change the variables in a ``Dataset``, you can use all the standard dictionary
methods, including ``values``, ``items``, ``__delitem__``, ``get`` and
:py:meth:`~xarray.Dataset.update`. Note that assigning a ``DataArray`` or pandas
object to a ``Dataset`` variable using ``__setitem__`` or ``update`` will
:ref:`automatically align<update>` the array(s) to the original
dataset's indexes.

You can copy a ``Dataset`` by calling the :py:meth:`~xarray.Dataset.copy`
method. By default, the copy is shallow, so only the container will be copied:
the arrays in the ``Dataset`` will still be stored in the same underlying
:py:class:`numpy.ndarray` objects. You can copy all data by calling
``ds.copy(deep=True)``.

.. _transforming datasets:

Transforming datasets
~~~~~~~~~~~~~~~~~~~~~

In addition to dictionary-like methods (described above), xarray has additional
methods (like pandas) for transforming datasets into new objects.

For removing variables, you can select and drop an explicit list of
variables by indexing with a list of names or using the
:py:meth:`~xarray.Dataset.drop_vars` methods to return a new ``Dataset``. These
operations keep around coordinates:

.. ipython:: python

    ds[["temperature"]]
    ds[["temperature", "temperature_double"]]
    ds.drop_vars("temperature")

To remove a dimension, you can use :py:meth:`~xarray.Dataset.drop_dims` method.
Any variables using that dimension are dropped:

.. ipython:: python

    ds.drop_dims("time")

As an alternate to dictionary-like modifications, you can use
:py:meth:`~xarray.Dataset.assign` and :py:meth:`~xarray.Dataset.assign_coords`.
These methods return a new dataset with additional (or replaced) values:

.. ipython:: python

    ds.assign(temperature2=2 * ds.temperature)

There is also the :py:meth:`~xarray.Dataset.pipe` method that allows you to use
a method call with an external function (e.g., ``ds.pipe(func)``) instead of
simply calling it (e.g., ``func(ds)``). This allows you to write pipelines for
transforming your data (using "method chaining") instead of writing hard to
follow nested function calls:

.. ipython:: python

    # these lines are equivalent, but with pipe we can make the logic flow
    # entirely from left to right
    plt.plot((2 * ds.temperature.sel(x=0)).mean("y"))
    (ds.temperature.sel(x=0).pipe(lambda x: 2 * x).mean("y").pipe(plt.plot))

Both ``pipe`` and ``assign`` replicate the pandas methods of the same names
(:py:meth:`DataFrame.pipe <pandas.DataFrame.pipe>` and
:py:meth:`DataFrame.assign <pandas.DataFrame.assign>`).

With xarray, there is no performance penalty for creating new datasets, even if
variables are lazily loaded from a file on disk. Creating new objects instead
of mutating existing objects often results in easier to understand code, so we
encourage using this approach.

Renaming variables
~~~~~~~~~~~~~~~~~~

Another useful option is the :py:meth:`~xarray.Dataset.rename` method to rename
dataset variables:

.. ipython:: python

    ds.rename({"temperature": "temp", "precipitation": "precip"})

The related :py:meth:`~xarray.Dataset.swap_dims` method allows you do to swap
dimension and non-dimension variables:

.. ipython:: python

    ds.coords["day"] = ("time", [6, 7, 8])
    ds.swap_dims({"time": "day"})

.. _coordinates:

Coordinates
-----------

Coordinates are ancillary variables stored for ``DataArray`` and ``Dataset``
objects in the ``coords`` attribute:

.. ipython:: python

    ds.coords

Unlike attributes, xarray *does* interpret and persist coordinates in
operations that transform xarray objects. There are two types of coordinates
in xarray:

- **dimension coordinates** are one dimensional coordinates with a name equal
  to their sole dimension (marked by ``*`` when printing a dataset or data
  array). They are used for label based indexing and alignment,
  like the ``index`` found on a pandas :py:class:`~pandas.DataFrame` or
  :py:class:`~pandas.Series`. Indeed, these "dimension" coordinates use a
  :py:class:`pandas.Index` internally to store their values.

- **non-dimension coordinates** are variables that contain coordinate
  data, but are not a dimension coordinate. They can be multidimensional (see
  :ref:`/examples/multidimensional-coords.ipynb`), and there is no
  relationship between the name of a non-dimension coordinate and the
  name(s) of its dimension(s).  Non-dimension coordinates can be
  useful for indexing or plotting; otherwise, xarray does not make any
  direct use of the values associated with them.  They are not used
  for alignment or automatic indexing, nor are they required to match
  when doing arithmetic (see :ref:`coordinates math`).

.. note::

  Xarray's terminology differs from the `CF terminology`_, where the
  "dimension coordinates" are called "coordinate variables", and the
  "non-dimension coordinates" are called "auxiliary coordinate variables"
  (see :issue:`1295` for more details).

.. _CF terminology: http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#terminology


Modifying coordinates
~~~~~~~~~~~~~~~~~~~~~

To entirely add or remove coordinate arrays, you can use dictionary like
syntax, as shown above.

To convert back and forth between data and coordinates, you can use the
:py:meth:`~xarray.Dataset.set_coords` and
:py:meth:`~xarray.Dataset.reset_coords` methods:

.. ipython:: python

    ds.reset_coords()
    ds.set_coords(["temperature", "precipitation"])
    ds["temperature"].reset_coords(drop=True)

Notice that these operations skip coordinates with names given by dimensions,
as used for indexing. This mostly because we are not entirely sure how to
design the interface around the fact that xarray cannot store a coordinate and
variable with the name but different values in the same dictionary. But we do
recognize that supporting something like this would be useful.

Coordinates methods
~~~~~~~~~~~~~~~~~~~

``Coordinates`` objects also have a few useful methods, mostly for converting
them into dataset objects:

.. ipython:: python

    ds.coords.to_dataset()

The merge method is particularly interesting, because it implements the same
logic used for merging coordinates in arithmetic operations
(see :ref:`comput`):

.. ipython:: python

    alt = xr.Dataset(coords={"z": [10], "lat": 0, "lon": 0})
    ds.coords.merge(alt.coords)

The ``coords.merge`` method may be useful if you want to implement your own
binary operations that act on xarray objects. In the future, we hope to write
more helper functions so that you can easily make your functions act like
xarray's built-in arithmetic.

Indexes
~~~~~~~

To convert a coordinate (or any ``DataArray``) into an actual
:py:class:`pandas.Index`, use the :py:meth:`~xarray.DataArray.to_index` method:

.. ipython:: python

    ds["time"].to_index()

A useful shortcut is the ``indexes`` property (on both ``DataArray`` and
``Dataset``), which lazily constructs a dictionary whose keys are given by each
dimension and whose the values are ``Index`` objects:

.. ipython:: python

    ds.indexes

MultiIndex coordinates
~~~~~~~~~~~~~~~~~~~~~~

Xarray supports labeling coordinate values with a :py:class:`pandas.MultiIndex`:

.. ipython:: python

    midx = pd.MultiIndex.from_arrays(
        [["R", "R", "V", "V"], [0.1, 0.2, 0.7, 0.9]], names=("band", "wn")
    )
    mda = xr.DataArray(np.random.rand(4), coords={"spec": midx}, dims="spec")
    mda

For convenience multi-index levels are directly accessible as "virtual" or
"derived" coordinates (marked by ``-`` when printing a dataset or data array):

.. ipython:: python

    mda["band"]
    mda.wn

Indexing with multi-index levels is also possible using the ``sel`` method
(see :ref:`multi-level indexing`).

Unlike other coordinates, "virtual" level coordinates are not stored in
the ``coords`` attribute of ``DataArray`` and ``Dataset`` objects
(although they are shown when printing the ``coords`` attribute).
Consequently, most of the coordinates related methods don't apply for them.
It also can't be used to replace one particular level.

Because in a ``DataArray`` or ``Dataset`` object each multi-index level is
accessible as a "virtual" coordinate, its name must not conflict with the names
of the other levels, coordinates and data variables of the same object.
Even though xarray sets default names for multi-indexes with unnamed levels,
it is recommended that you explicitly set the names of the levels.

.. [1] Latitude and longitude are 2D arrays because the dataset uses
   `projected coordinates`__. ``reference_time`` refers to the reference time
   at which the forecast was made, rather than ``time`` which is the valid time
   for which the forecast applies.

__ http://en.wikipedia.org/wiki/Map_projection
.. currentmodule:: xarray
.. _pandas:

===================
Working with pandas
===================

One of the most important features of xarray is the ability to convert to and
from :py:mod:`pandas` objects to interact with the rest of the PyData
ecosystem. For example, for plotting labeled data, we highly recommend
using the visualization `built in to pandas itself`__ or provided by the pandas
aware libraries such as `Seaborn`__.

__ http://pandas.pydata.org/pandas-docs/stable/visualization.html
__ http://seaborn.pydata.org/

.. ipython:: python
    :suppress:

    import numpy as np
    import pandas as pd
    import xarray as xr

    np.random.seed(123456)

Hierarchical and tidy data
~~~~~~~~~~~~~~~~~~~~~~~~~~

Tabular data is easiest to work with when it meets the criteria for
`tidy data`__:

* Each column holds a different variable.
* Each rows holds a different observation.

__ http://www.jstatsoft.org/v59/i10/

In this "tidy data" format, we can represent any :py:class:`Dataset` and
:py:class:`DataArray` in terms of :py:class:`~pandas.DataFrame` and
:py:class:`~pandas.Series`, respectively (and vice-versa). The representation
works by flattening non-coordinates to 1D, and turning the tensor product of
coordinate indexes into a :py:class:`pandas.MultiIndex`.

Dataset and DataFrame
---------------------

To convert any dataset to a ``DataFrame`` in tidy form, use the
:py:meth:`Dataset.to_dataframe()` method:

.. ipython:: python

    ds = xr.Dataset(
        {"foo": (("x", "y"), np.random.randn(2, 3))},
        coords={
            "x": [10, 20],
            "y": ["a", "b", "c"],
            "along_x": ("x", np.random.randn(2)),
            "scalar": 123,
        },
    )
    ds
    df = ds.to_dataframe()
    df

We see that each variable and coordinate in the Dataset is now a column in the
DataFrame, with the exception of indexes which are in the index.
To convert the ``DataFrame`` to any other convenient representation,
use ``DataFrame`` methods like :py:meth:`~pandas.DataFrame.reset_index`,
:py:meth:`~pandas.DataFrame.stack` and :py:meth:`~pandas.DataFrame.unstack`.

For datasets containing dask arrays where the data should be lazily loaded, see the
:py:meth:`Dataset.to_dask_dataframe()` method.

To create a ``Dataset`` from a ``DataFrame``, use the
:py:meth:`Dataset.from_dataframe` class method or the equivalent
:py:meth:`pandas.DataFrame.to_xarray` method:

.. ipython:: python

    xr.Dataset.from_dataframe(df)

Notice that that dimensions of variables in the ``Dataset`` have now
expanded after the round-trip conversion to a ``DataFrame``. This is because
every object in a ``DataFrame`` must have the same indices, so we need to
broadcast the data of each array to the full size of the new ``MultiIndex``.

Likewise, all the coordinates (other than indexes) ended up as variables,
because pandas does not distinguish non-index coordinates.

DataArray and Series
--------------------

``DataArray`` objects have a complementary representation in terms of a
:py:class:`~pandas.Series`. Using a Series preserves the ``Dataset`` to
``DataArray`` relationship, because ``DataFrames`` are dict-like containers
of ``Series``. The methods are very similar to those for working with
DataFrames:

.. ipython:: python

    s = ds["foo"].to_series()
    s
    # or equivalently, with Series.to_xarray()
    xr.DataArray.from_series(s)

Both the ``from_series`` and ``from_dataframe`` methods use reindexing, so they
work even if not the hierarchical index is not a full tensor product:

.. ipython:: python

    s[::2]
    s[::2].to_xarray()

Multi-dimensional data
~~~~~~~~~~~~~~~~~~~~~~

Tidy data is great, but it sometimes you want to preserve dimensions instead of
automatically stacking them into a ``MultiIndex``.

:py:meth:`DataArray.to_pandas()` is a shortcut that lets you convert a
DataArray directly into a pandas object with the same dimensionality, if
available in pandas (i.e., a 1D array is converted to a
:py:class:`~pandas.Series` and 2D to :py:class:`~pandas.DataFrame`):

.. ipython:: python

    arr = xr.DataArray(
        np.random.randn(2, 3), coords=[("x", [10, 20]), ("y", ["a", "b", "c"])]
    )
    df = arr.to_pandas()
    df

To perform the inverse operation of converting any pandas objects into a data
array with the same shape, simply use the :py:class:`DataArray`
constructor:

.. ipython:: python

    xr.DataArray(df)

Both the ``DataArray`` and ``Dataset`` constructors directly convert pandas
objects into xarray objects with the same shape. This means that they
preserve all use of multi-indexes:

.. ipython:: python

    index = pd.MultiIndex.from_arrays(
        [["a", "a", "b"], [0, 1, 2]], names=["one", "two"]
    )
    df = pd.DataFrame({"x": 1, "y": 2}, index=index)
    ds = xr.Dataset(df)
    ds

However, you will need to set dimension names explicitly, either with the
``dims`` argument on in the ``DataArray`` constructor or by calling
:py:class:`~Dataset.rename` on the new object.

.. _panel transition:

Transitioning from pandas.Panel to xarray
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Panel``, pandas' data structure for 3D arrays, was always a second class
data structure compared to the Series and DataFrame. To allow pandas
developers to focus more on its core functionality built around the
DataFrame, pandas removed ``Panel`` in favor of directing users who use
multi-dimensional arrays to xarray.

Xarray has most of ``Panel``'s features, a more explicit API (particularly around
indexing), and the ability to scale to >3 dimensions with the same interface.

As discussed :ref:`elsewhere <data structures>` in the docs, there are two primary data structures in
xarray: ``DataArray`` and ``Dataset``. You can imagine a ``DataArray`` as a
n-dimensional pandas ``Series`` (i.e. a single typed array), and a ``Dataset``
as the ``DataFrame`` equivalent (i.e. a dict of aligned ``DataArray`` objects).

So you can represent a Panel, in two ways:

- As a 3-dimensional ``DataArray``,
- Or as a ``Dataset`` containing a number of 2-dimensional DataArray objects.

Let's take a look:

.. ipython:: python

    data = np.random.RandomState(0).rand(2, 3, 4)
    items = list("ab")
    major_axis = list("mno")
    minor_axis = pd.date_range(start="2000", periods=4, name="date")

With old versions of pandas (prior to 0.25), this could stored in a ``Panel``:

.. ipython::
    :verbatim:

    In [1]: pd.Panel(data, items, major_axis, minor_axis)
    Out[1]:
    <class 'pandas.core.panel.Panel'>
    Dimensions: 2 (items) x 3 (major_axis) x 4 (minor_axis)
    Items axis: a to b
    Major_axis axis: m to o
    Minor_axis axis: 2000-01-01 00:00:00 to 2000-01-04 00:00:00

To put this data in a ``DataArray``, write:

.. ipython:: python

    array = xr.DataArray(data, [items, major_axis, minor_axis])
    array

As you can see, there are three dimensions (each is also a coordinate). Two of
the axes of were unnamed, so have been assigned ``dim_0`` and ``dim_1``
respectively, while the third retains its name ``date``.

You can also easily convert this data into ``Dataset``:

.. ipython:: python

    array.to_dataset(dim="dim_0")

Here, there are two data variables, each representing a DataFrame on panel's
``items`` axis, and labeled as such. Each variable is a 2D array of the
respective values along the ``items`` dimension.

While the xarray docs are relatively complete, a few items stand out for Panel users:

- A DataArray's data is stored as a numpy array, and so can only contain a single
  type. As a result, a Panel that contains :py:class:`~pandas.DataFrame` objects
  with multiple types will be converted to ``dtype=object``. A ``Dataset`` of
  multiple ``DataArray`` objects each with its own dtype will allow original
  types to be preserved.
- :ref:`Indexing <indexing>` is similar to pandas, but more explicit and
  leverages xarray's naming of dimensions.
- Because of those features, making much higher dimensional data is very
  practical.
- Variables in ``Dataset`` objects can use a subset of its dimensions. For
  example, you can have one dataset with Person x Score x Time, and another with
  Person x Score.
- You can use coordinates are used for both dimensions and for variables which
  _label_ the data variables, so you could have a coordinate Age, that labelled
  the Person dimension of a Dataset of Person x Score x Time.

While xarray may take some getting used to, it's worth it! If anything is unclear,
please post an issue on `GitHub <https://github.com/pydata/xarray>`__ or
`StackOverflow <http://stackoverflow.com/questions/tagged/python-xarray>`__,
and we'll endeavor to respond to the specific case or improve the general docs.
.. currentmodule:: xarray

.. _weather-climate:

Weather and climate data
========================

.. ipython:: python
    :suppress:

    import xarray as xr

Xarray can leverage metadata that follows the `Climate and Forecast (CF) conventions`_ if present. Examples include automatic labelling of plots with descriptive names and units if proper metadata is present (see :ref:`plotting`) and support for non-standard calendars used in climate science through the ``cftime`` module (see :ref:`CFTimeIndex`). There are also a number of geosciences-focused projects that build on xarray (see :ref:`ecosystem`).

.. _Climate and Forecast (CF) conventions: http://cfconventions.org

.. _cf_variables:

Related Variables
-----------------

Several CF variable attributes contain lists of other variables
associated with the variable with the attribute.  A few of these are
now parsed by xarray, with the attribute value popped to encoding on
read and the variables in that value interpreted as non-dimension
coordinates:

- ``coordinates``
- ``bounds``
- ``grid_mapping``
- ``climatology``
- ``geometry``
- ``node_coordinates``
- ``node_count``
- ``part_node_count``
- ``interior_ring``
- ``cell_measures``
- ``formula_terms``

This decoding is controlled by the ``decode_coords`` kwarg to
:py:func:`open_dataset` and :py:func:`open_mfdataset`.

The CF attribute ``ancillary_variables`` was not included in the list
due to the variables listed there being associated primarily with the
variable with the attribute, rather than with the dimensions.

.. _metpy_accessor:

CF-compliant coordinate variables
---------------------------------

`MetPy`_ adds a	``metpy`` accessor that allows accessing coordinates with appropriate CF metadata using generic names ``x``, ``y``, ``vertical`` and ``time``. There is also a `cartopy_crs` attribute that provides projection information, parsed from the appropriate CF metadata, as a `Cartopy`_ projection object. See `their documentation`_ for more information.

.. _`MetPy`: https://unidata.github.io/MetPy/dev/index.html
.. _`their documentation`:	https://unidata.github.io/MetPy/dev/tutorials/xarray_tutorial.html#coordinates
.. _`Cartopy`: https://scitools.org.uk/cartopy/docs/latest/crs/projections.html

.. _CFTimeIndex:

Non-standard calendars and dates outside the Timestamp-valid range
------------------------------------------------------------------

Through the standalone ``cftime`` library and a custom subclass of
:py:class:`pandas.Index`, xarray supports a subset of the indexing
functionality enabled through the standard :py:class:`pandas.DatetimeIndex` for
dates from non-standard calendars commonly used in climate science or dates
using a standard calendar, but outside the `Timestamp-valid range`_
(approximately between years 1678 and 2262).

.. note::

   As of xarray version 0.11, by default, :py:class:`cftime.datetime` objects
   will be used to represent times (either in indexes, as a
   :py:class:`~xarray.CFTimeIndex`, or in data arrays with dtype object) if
   any of the following are true:

   - The dates are from a non-standard calendar
   - Any dates are outside the Timestamp-valid range.

   Otherwise pandas-compatible dates from a standard calendar will be
   represented with the ``np.datetime64[ns]`` data type, enabling the use of a
   :py:class:`pandas.DatetimeIndex` or arrays with dtype ``np.datetime64[ns]``
   and their full set of associated features.

For example, you can create a DataArray indexed by a time
coordinate with dates from a no-leap calendar and a
:py:class:`~xarray.CFTimeIndex` will automatically be used:

.. ipython:: python

    from itertools import product
    from cftime import DatetimeNoLeap

    dates = [
        DatetimeNoLeap(year, month, 1)
        for year, month in product(range(1, 3), range(1, 13))
    ]
    da = xr.DataArray(np.arange(24), coords=[dates], dims=["time"], name="foo")

Xarray also includes a :py:func:`~xarray.cftime_range` function, which enables
creating a :py:class:`~xarray.CFTimeIndex` with regularly-spaced dates.  For
instance, we can create the same dates and DataArray we created above using:

.. ipython:: python

    dates = xr.cftime_range(start="0001", periods=24, freq="MS", calendar="noleap")
    da = xr.DataArray(np.arange(24), coords=[dates], dims=["time"], name="foo")

Mirroring pandas' method with the same name, :py:meth:`~xarray.infer_freq` allows one to
infer the sampling frequency of a :py:class:`~xarray.CFTimeIndex` or a 1-D
:py:class:`~xarray.DataArray` containing cftime objects. It also works transparently with
``np.datetime64[ns]`` and ``np.timedelta64[ns]`` data.

.. ipython:: python

    xr.infer_freq(dates)

With :py:meth:`~xarray.CFTimeIndex.strftime` we can also easily generate formatted strings from
the datetime values of a :py:class:`~xarray.CFTimeIndex` directly or through the
``dt`` accessor for a :py:class:`~xarray.DataArray`
using the same formatting as the standard `datetime.strftime`_ convention .

.. _datetime.strftime: https://docs.python.org/3/library/datetime.html#strftime-strptime-behavior

.. ipython:: python

    dates.strftime("%c")
    da["time"].dt.strftime("%Y%m%d")

Conversion between non-standard calendar and to/from pandas DatetimeIndexes is
facilitated with the :py:meth:`xarray.Dataset.convert_calendar` method (also available as
:py:meth:`xarray.DataArray.convert_calendar`). Here, like elsewhere in xarray, the ``use_cftime``
argument controls which datetime backend is used in the output. The default (``None``) is to
use `pandas` when possible, i.e. when the calendar is standard and dates are within 1678 and 2262.

.. ipython:: python

    dates = xr.cftime_range(start="2001", periods=24, freq="MS", calendar="noleap")
    da_nl = xr.DataArray(np.arange(24), coords=[dates], dims=["time"], name="foo")
    da_std = da.convert_calendar("standard", use_cftime=True)

The data is unchanged, only the timestamps are modified. Further options are implemented
for the special ``"360_day"`` calendar and for handling missing dates. There is also
:py:meth:`xarray.Dataset.interp_calendar` (and :py:meth:`xarray.DataArray.interp_calendar`)
for `interpolating` data between calendars.

For data indexed by a :py:class:`~xarray.CFTimeIndex` xarray currently supports:

- `Partial datetime string indexing`_:

.. ipython:: python

    da.sel(time="0001")
    da.sel(time=slice("0001-05", "0002-02"))

.. note::


   For specifying full or partial datetime strings in cftime
   indexing, xarray supports two versions of the `ISO 8601 standard`_, the
   basic pattern (YYYYMMDDhhmmss) or the extended pattern
   (YYYY-MM-DDThh:mm:ss), as well as the default cftime string format
   (YYYY-MM-DD hh:mm:ss).  This is somewhat more restrictive than pandas;
   in other words, some datetime strings that would be valid for a
   :py:class:`pandas.DatetimeIndex` are not valid for an
   :py:class:`~xarray.CFTimeIndex`.

- Access of basic datetime components via the ``dt`` accessor (in this case
  just "year", "month", "day", "hour", "minute", "second", "microsecond",
  "season", "dayofyear", "dayofweek", and "days_in_month") with the addition
  of "calendar", absent from pandas:

.. ipython:: python

    da.time.dt.year
    da.time.dt.month
    da.time.dt.season
    da.time.dt.dayofyear
    da.time.dt.dayofweek
    da.time.dt.days_in_month
    da.time.dt.calendar

- Rounding of datetimes to fixed frequencies via the ``dt`` accessor:

.. ipython:: python

    da.time.dt.ceil("3D")
    da.time.dt.floor("5D")
    da.time.dt.round("2D")

- Group-by operations based on datetime accessor attributes (e.g. by month of
  the year):

.. ipython:: python

    da.groupby("time.month").sum()

- Interpolation using :py:class:`cftime.datetime` objects:

.. ipython:: python

    da.interp(time=[DatetimeNoLeap(1, 1, 15), DatetimeNoLeap(1, 2, 15)])

- Interpolation using datetime strings:

.. ipython:: python

    da.interp(time=["0001-01-15", "0001-02-15"])

- Differentiation:

.. ipython:: python

    da.differentiate("time")

- Serialization:

.. ipython:: python

    da.to_netcdf("example-no-leap.nc")
    xr.open_dataset("example-no-leap.nc")

.. ipython:: python
    :suppress:

    import os

    os.remove("example-no-leap.nc")

- And resampling along the time dimension for data indexed by a :py:class:`~xarray.CFTimeIndex`:

.. ipython:: python

    da.resample(time="81T", closed="right", label="right", base=3).mean()

.. _Timestamp-valid range: https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#timestamp-limitations
.. _ISO 8601 standard: https://en.wikipedia.org/wiki/ISO_8601
.. _partial datetime string indexing: https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#partial-string-indexing
###########
User Guide
###########

In this user guide, you will find detailed descriptions and
examples that describe many common tasks that you can accomplish with xarray.


.. toctree::
   :maxdepth: 2
   :hidden:

   terminology
   data-structures
   indexing
   interpolation
   computation
   groupby
   reshaping
   combining
   time-series
   weather-climate
   pandas
   io
   dask
   plotting
   duckarrays
.. currentmodule:: xarray
.. _io:

Reading and writing files
=========================

Xarray supports direct serialization and IO to several file formats, from
simple :ref:`io.pickle` files to the more flexible :ref:`io.netcdf`
format (recommended).

.. ipython:: python
    :suppress:

    import numpy as np
    import pandas as pd
    import xarray as xr

    np.random.seed(123456)

.. _io.netcdf:

netCDF
------

The recommended way to store xarray data structures is `netCDF`__, which
is a binary file format for self-described datasets that originated
in the geosciences. Xarray is based on the netCDF data model, so netCDF files
on disk directly correspond to :py:class:`Dataset` objects (more accurately,
a group in a netCDF file directly corresponds to a :py:class:`Dataset` object.
See :ref:`io.netcdf_groups` for more.)

NetCDF is supported on almost all platforms, and parsers exist
for the vast majority of scientific programming languages. Recent versions of
netCDF are based on the even more widely used HDF5 file-format.

__ http://www.unidata.ucar.edu/software/netcdf/

.. tip::

    If you aren't familiar with this data format, the `netCDF FAQ`_ is a good
    place to start.

.. _netCDF FAQ: http://www.unidata.ucar.edu/software/netcdf/docs/faq.html#What-Is-netCDF

Reading and writing netCDF files with xarray requires scipy or the
`netCDF4-Python`__ library to be installed (the latter is required to
read/write netCDF V4 files and use the compression options described below).

__ https://github.com/Unidata/netcdf4-python

We can save a Dataset to disk using the
:py:meth:`Dataset.to_netcdf` method:

.. ipython:: python

    ds = xr.Dataset(
        {"foo": (("x", "y"), np.random.rand(4, 5))},
        coords={
            "x": [10, 20, 30, 40],
            "y": pd.date_range("2000-01-01", periods=5),
            "z": ("x", list("abcd")),
        },
    )

    ds.to_netcdf("saved_on_disk.nc")

By default, the file is saved as netCDF4 (assuming netCDF4-Python is
installed). You can control the format and engine used to write the file with
the ``format`` and ``engine`` arguments.

.. tip::

   Using the `h5netcdf <https://github.com/shoyer/h5netcdf>`_  package
   by passing ``engine='h5netcdf'`` to :py:meth:`open_dataset` can
   sometimes be quicker than the default ``engine='netcdf4'`` that uses the
   `netCDF4 <https://github.com/Unidata/netcdf4-python>`_ package.


We can load netCDF files to create a new Dataset using
:py:func:`open_dataset`:

.. ipython:: python

    ds_disk = xr.open_dataset("saved_on_disk.nc")
    ds_disk

Similarly, a DataArray can be saved to disk using the
:py:meth:`DataArray.to_netcdf` method, and loaded
from disk using the :py:func:`open_dataarray` function. As netCDF files
correspond to :py:class:`Dataset` objects, these functions internally
convert the ``DataArray`` to a ``Dataset`` before saving, and then convert back
when loading, ensuring that the ``DataArray`` that is loaded is always exactly
the same as the one that was saved.

A dataset can also be loaded or written to a specific group within a netCDF
file. To load from a group, pass a ``group`` keyword argument to the
``open_dataset`` function. The group can be specified as a path-like
string, e.g., to access subgroup 'bar' within group 'foo' pass
'/foo/bar' as the ``group`` argument. When writing multiple groups in one file,
pass ``mode='a'`` to ``to_netcdf`` to ensure that each call does not delete the
file.

Data is *always* loaded lazily from netCDF files. You can manipulate, slice and subset
Dataset and DataArray objects, and no array values are loaded into memory until
you try to perform some sort of actual computation. For an example of how these
lazy arrays work, see the OPeNDAP section below.

There may be minor differences in the :py:class:`Dataset` object returned
when reading a NetCDF file with different engines. For example,
single-valued attributes are returned as scalars by the default
``engine=netcdf4``, but as arrays of size ``(1,)`` when reading with
``engine=h5netcdf``.

It is important to note that when you modify values of a Dataset, even one
linked to files on disk, only the in-memory copy you are manipulating in xarray
is modified: the original file on disk is never touched.

.. tip::

    Xarray's lazy loading of remote or on-disk datasets is often but not always
    desirable. Before performing computationally intense operations, it is
    often a good idea to load a Dataset (or DataArray) entirely into memory by
    invoking the :py:meth:`Dataset.load` method.

Datasets have a :py:meth:`Dataset.close` method to close the associated
netCDF file. However, it's often cleaner to use a ``with`` statement:

.. ipython:: python

    # this automatically closes the dataset after use
    with xr.open_dataset("saved_on_disk.nc") as ds:
        print(ds.keys())

Although xarray provides reasonable support for incremental reads of files on
disk, it does not support incremental writes, which can be a useful strategy
for dealing with datasets too big to fit into memory. Instead, xarray integrates
with dask.array (see :ref:`dask`), which provides a fully featured engine for
streaming computation.

It is possible to append or overwrite netCDF variables using the ``mode='a'``
argument. When using this option, all variables in the dataset will be written
to the original netCDF file, regardless if they exist in the original dataset.


.. _io.netcdf_groups:

Groups
~~~~~~

NetCDF groups are not supported as part of the :py:class:`Dataset` data model.
Instead, groups can be loaded individually as Dataset objects.
To do so, pass a ``group`` keyword argument to the
:py:func:`open_dataset` function. The group can be specified as a path-like
string, e.g., to access subgroup ``'bar'`` within group ``'foo'`` pass
``'/foo/bar'`` as the ``group`` argument.
In a similar way, the ``group`` keyword argument can be given to the
:py:meth:`Dataset.to_netcdf` method to write to a group
in a netCDF file.
When writing multiple groups in one file, pass ``mode='a'`` to
:py:meth:`Dataset.to_netcdf` to ensure that each call does not delete the file.

.. _io.encoding:

Reading encoded data
~~~~~~~~~~~~~~~~~~~~

NetCDF files follow some conventions for encoding datetime arrays (as numbers
with a "units" attribute) and for packing and unpacking data (as
described by the "scale_factor" and "add_offset" attributes). If the argument
``decode_cf=True`` (default) is given to :py:func:`open_dataset`, xarray will attempt
to automatically decode the values in the netCDF objects according to
`CF conventions`_. Sometimes this will fail, for example, if a variable
has an invalid "units" or "calendar" attribute. For these cases, you can
turn this decoding off manually.

.. _CF conventions: http://cfconventions.org/

You can view this encoding information (among others) in the
:py:attr:`DataArray.encoding` and
:py:attr:`DataArray.encoding` attributes:

.. ipython::
    :verbatim:

    In [1]: ds_disk["y"].encoding
    Out[1]:
    {'zlib': False,
     'shuffle': False,
     'complevel': 0,
     'fletcher32': False,
     'contiguous': True,
     'chunksizes': None,
     'source': 'saved_on_disk.nc',
     'original_shape': (5,),
     'dtype': dtype('int64'),
     'units': 'days since 2000-01-01 00:00:00',
     'calendar': 'proleptic_gregorian'}

    In [9]: ds_disk.encoding
    Out[9]:
    {'unlimited_dims': set(),
     'source': 'saved_on_disk.nc'}

Note that all operations that manipulate variables other than indexing
will remove encoding information.

.. ipython:: python
    :suppress:

    ds_disk.close()


.. _combining multiple files:

Reading multi-file datasets
...........................

NetCDF files are often encountered in collections, e.g., with different files
corresponding to different model runs or one file per timestamp.
Xarray can straightforwardly combine such files into a single Dataset by making use of
:py:func:`concat`, :py:func:`merge`, :py:func:`combine_nested` and
:py:func:`combine_by_coords`. For details on the difference between these
functions see :ref:`combining data`.

Xarray includes support for manipulating datasets that don't fit into memory
with dask_. If you have dask installed, you can open multiple files
simultaneously in parallel using :py:func:`open_mfdataset`::

    xr.open_mfdataset('my/files/*.nc', parallel=True)

This function automatically concatenates and merges multiple files into a
single xarray dataset.
It is the recommended way to open multiple files with xarray.
For more details on parallel reading, see :ref:`combining.multi`, :ref:`dask.io` and a
`blog post`_ by Stephan Hoyer.
:py:func:`open_mfdataset` takes many kwargs that allow you to
control its behaviour (for e.g. ``parallel``, ``combine``, ``compat``, ``join``, ``concat_dim``).
See its docstring for more details.


.. note::

    A common use-case involves a dataset distributed across a large number of files with
    each file containing a large number of variables. Commonly, a few of these variables
    need to be concatenated along a dimension (say ``"time"``), while the rest are equal
    across the datasets (ignoring floating point differences). The following command
    with suitable modifications (such as ``parallel=True``) works well with such datasets::

         xr.open_mfdataset('my/files/*.nc', concat_dim="time", combine="nested",
     	              	   data_vars='minimal', coords='minimal', compat='override')

    This command concatenates variables along the ``"time"`` dimension, but only those that
    already contain the ``"time"`` dimension (``data_vars='minimal', coords='minimal'``).
    Variables that lack the ``"time"`` dimension are taken from the first dataset
    (``compat='override'``).


.. _dask: http://dask.pydata.org
.. _blog post: http://stephanhoyer.com/2015/06/11/xray-dask-out-of-core-labeled-arrays/

Sometimes multi-file datasets are not conveniently organized for easy use of :py:func:`open_mfdataset`.
One can use the ``preprocess`` argument to provide a function that takes a dataset
and returns a modified Dataset.
:py:func:`open_mfdataset` will call ``preprocess`` on every dataset
(corresponding to each file) prior to combining them.


If :py:func:`open_mfdataset` does not meet your needs, other approaches are possible.
The general pattern for parallel reading of multiple files
using dask, modifying those datasets and then combining into a single ``Dataset`` is::

     def modify(ds):
         # modify ds here
         return ds


     # this is basically what open_mfdataset does
     open_kwargs = dict(decode_cf=True, decode_times=False)
     open_tasks = [dask.delayed(xr.open_dataset)(f, **open_kwargs) for f in file_names]
     tasks = [dask.delayed(modify)(task) for task in open_tasks]
     datasets = dask.compute(tasks)  # get a list of xarray.Datasets
     combined = xr.combine_nested(datasets)  # or some combination of concat, merge


As an example, here's how we could approximate ``MFDataset`` from the netCDF4
library::

    from glob import glob
    import xarray as xr

    def read_netcdfs(files, dim):
        # glob expands paths with * to a list of files, like the unix shell
        paths = sorted(glob(files))
        datasets = [xr.open_dataset(p) for p in paths]
        combined = xr.concat(datasets, dim)
        return combined

    combined = read_netcdfs('/all/my/files/*.nc', dim='time')

This function will work in many cases, but it's not very robust. First, it
never closes files, which means it will fail if you need to load more than
a few thousand files. Second, it assumes that you want all the data from each
file and that it can all fit into memory. In many situations, you only need
a small subset or an aggregated summary of the data from each file.

Here's a slightly more sophisticated example of how to remedy these
deficiencies::

    def read_netcdfs(files, dim, transform_func=None):
        def process_one_path(path):
            # use a context manager, to ensure the file gets closed after use
            with xr.open_dataset(path) as ds:
                # transform_func should do some sort of selection or
                # aggregation
                if transform_func is not None:
                    ds = transform_func(ds)
                # load all data from the transformed dataset, to ensure we can
                # use it after closing each original file
                ds.load()
                return ds

        paths = sorted(glob(files))
        datasets = [process_one_path(p) for p in paths]
        combined = xr.concat(datasets, dim)
        return combined

    # here we suppose we only care about the combined mean of each file;
    # you might also use indexing operations like .sel to subset datasets
    combined = read_netcdfs('/all/my/files/*.nc', dim='time',
                            transform_func=lambda ds: ds.mean())

This pattern works well and is very robust. We've used similar code to process
tens of thousands of files constituting 100s of GB of data.


.. _io.netcdf.writing_encoded:

Writing encoded data
~~~~~~~~~~~~~~~~~~~~

Conversely, you can customize how xarray writes netCDF files on disk by
providing explicit encodings for each dataset variable. The ``encoding``
argument takes a dictionary with variable names as keys and variable specific
encodings as values. These encodings are saved as attributes on the netCDF
variables on disk, which allows xarray to faithfully read encoded data back into
memory.

It is important to note that using encodings is entirely optional: if you do not
supply any of these encoding options, xarray will write data to disk using a
default encoding, or the options in the ``encoding`` attribute, if set.
This works perfectly fine in most cases, but encoding can be useful for
additional control, especially for enabling compression.

In the file on disk, these encodings are saved as attributes on each variable, which
allow xarray and other CF-compliant tools for working with netCDF files to correctly
read the data.

Scaling and type conversions
............................

These encoding options work on any version of the netCDF file format:

- ``dtype``: Any valid NumPy dtype or string convertible to a dtype, e.g., ``'int16'``
  or ``'float32'``. This controls the type of the data written on disk.
- ``_FillValue``:  Values of ``NaN`` in xarray variables are remapped to this value when
  saved on disk. This is important when converting floating point with missing values
  to integers on disk, because ``NaN`` is not a valid value for integer dtypes. By
  default, variables with float types are attributed a ``_FillValue`` of ``NaN`` in the
  output file, unless explicitly disabled with an encoding ``{'_FillValue': None}``.
- ``scale_factor`` and ``add_offset``: Used to convert from encoded data on disk to
  to the decoded data in memory, according to the formula
  ``decoded = scale_factor * encoded + add_offset``.

These parameters can be fruitfully combined to compress discretized data on disk. For
example, to save the variable ``foo`` with a precision of 0.1 in 16-bit integers while
converting ``NaN`` to ``-9999``, we would use
``encoding={'foo': {'dtype': 'int16', 'scale_factor': 0.1, '_FillValue': -9999}}``.
Compression and decompression with such discretization is extremely fast.

.. _io.string-encoding:

String encoding
...............

Xarray can write unicode strings to netCDF files in two ways:

- As variable length strings. This is only supported on netCDF4 (HDF5) files.
- By encoding strings into bytes, and writing encoded bytes as a character
  array. The default encoding is UTF-8.

By default, we use variable length strings for compatible files and fall-back
to using encoded character arrays. Character arrays can be selected even for
netCDF4 files by setting the ``dtype`` field in ``encoding`` to ``S1``
(corresponding to NumPy's single-character bytes dtype).

If character arrays are used:

- The string encoding that was used is stored on
  disk in the ``_Encoding`` attribute, which matches an ad-hoc convention
  `adopted by the netCDF4-Python library <https://github.com/Unidata/netcdf4-python/pull/665>`_.
  At the time of this writing (October 2017), a standard convention for indicating
  string encoding for character arrays in netCDF files was
  `still under discussion <https://github.com/Unidata/netcdf-c/issues/402>`_.
  Technically, you can use
  `any string encoding recognized by Python <https://docs.python.org/3/library/codecs.html#standard-encodings>`_ if you feel the need to deviate from UTF-8,
  by setting the ``_Encoding`` field in ``encoding``. But
  `we don't recommend it <http://utf8everywhere.org/>`_.
- The character dimension name can be specified by the ``char_dim_name`` field of a variable's
  ``encoding``. If the name of the character dimension is not specified, the default is
  ``f'string{data.shape[-1]}'``. When decoding character arrays from existing files, the
  ``char_dim_name`` is added to the variables ``encoding`` to preserve if encoding happens, but
  the field can be edited by the user.

.. warning::

  Missing values in bytes or unicode string arrays (represented by ``NaN`` in
  xarray) are currently written to disk as empty strings ``''``. This means
  missing values will not be restored when data is loaded from disk.
  This behavior is likely to change in the future (:issue:`1647`).
  Unfortunately, explicitly setting a ``_FillValue`` for string arrays to handle
  missing values doesn't work yet either, though we also hope to fix this in the
  future.

Chunk based compression
.......................

``zlib``, ``complevel``, ``fletcher32``, ``continguous`` and ``chunksizes``
can be used for enabling netCDF4/HDF5's chunk based compression, as described
in the `documentation for createVariable`_ for netCDF4-Python. This only works
for netCDF4 files and thus requires using ``format='netCDF4'`` and either
``engine='netcdf4'`` or ``engine='h5netcdf'``.

.. _documentation for createVariable: http://unidata.github.io/netcdf4-python/#netCDF4.Dataset.createVariable

Chunk based gzip compression can yield impressive space savings, especially
for sparse data, but it comes with significant performance overhead. HDF5
libraries can only read complete chunks back into memory, and maximum
decompression speed is in the range of 50-100 MB/s. Worse, HDF5's compression
and decompression currently cannot be parallelized with dask. For these reasons, we
recommend trying discretization based compression (described above) first.

Time units
..........

The ``units`` and ``calendar`` attributes control how xarray serializes ``datetime64`` and
``timedelta64`` arrays to datasets on disk as numeric values. The ``units`` encoding
should be a string like ``'days since 1900-01-01'`` for ``datetime64`` data or a string
like ``'days'`` for ``timedelta64`` data. ``calendar`` should be one of the calendar types
supported by netCDF4-python: 'standard', 'gregorian', 'proleptic_gregorian' 'noleap',
'365_day', '360_day', 'julian', 'all_leap', '366_day'.

By default, xarray uses the ``'proleptic_gregorian'`` calendar and units of the smallest time
difference between values, with a reference time of the first time value.


.. _io.coordinates:

Coordinates
...........

You can control the ``coordinates`` attribute written to disk by specifying ``DataArray.encoding["coordinates"]``.
If not specified, xarray automatically sets ``DataArray.encoding["coordinates"]`` to a space-delimited list
of names of coordinate variables that share dimensions with the ``DataArray`` being written.
This allows perfect roundtripping of xarray datasets but may not be desirable.
When an xarray ``Dataset`` contains non-dimensional coordinates that do not share dimensions with any of
the variables, these coordinate variable names are saved under a "global" ``"coordinates"`` attribute.
This is not CF-compliant but again facilitates roundtripping of xarray datasets.

Invalid netCDF files
~~~~~~~~~~~~~~~~~~~~

The library ``h5netcdf`` allows writing some dtypes (booleans, complex, ...) that aren't
allowed in netCDF4 (see
`h5netcdf documentation <https://github.com/shoyer/h5netcdf#invalid-netcdf-files>`_).
This feature is available through :py:meth:`DataArray.to_netcdf` and
:py:meth:`Dataset.to_netcdf` when used with ``engine="h5netcdf"``
and currently raises a warning unless ``invalid_netcdf=True`` is set:

.. ipython:: python
    :okwarning:

    # Writing complex valued data
    da = xr.DataArray([1.0 + 1.0j, 2.0 + 2.0j, 3.0 + 3.0j])
    da.to_netcdf("complex.nc", engine="h5netcdf", invalid_netcdf=True)

    # Reading it back
    xr.open_dataarray("complex.nc", engine="h5netcdf")

.. ipython:: python
    :suppress:

    import os

    os.remove("complex.nc")

.. warning::

  Note that this produces a file that is likely to be not readable by other netCDF
  libraries!

.. _io.iris:

Iris
----

The Iris_ tool allows easy reading of common meteorological and climate model formats
(including GRIB and UK MetOffice PP files) into ``Cube`` objects which are in many ways very
similar to ``DataArray`` objects, while enforcing a CF-compliant data model. If iris is
installed, xarray can convert a ``DataArray`` into a ``Cube`` using
:py:meth:`DataArray.to_iris`:

.. ipython:: python

    da = xr.DataArray(
        np.random.rand(4, 5),
        dims=["x", "y"],
        coords=dict(x=[10, 20, 30, 40], y=pd.date_range("2000-01-01", periods=5)),
    )

    cube = da.to_iris()
    cube

Conversely, we can create a new ``DataArray`` object from a ``Cube`` using
:py:meth:`DataArray.from_iris`:

.. ipython:: python

    da_cube = xr.DataArray.from_iris(cube)
    da_cube


.. _Iris: http://scitools.org.uk/iris


OPeNDAP
-------

Xarray includes support for `OPeNDAP`__ (via the netCDF4 library or Pydap), which
lets us access large datasets over HTTP.

__ http://www.opendap.org/

For example, we can open a connection to GBs of weather data produced by the
`PRISM`__ project, and hosted by `IRI`__ at Columbia:

__ http://www.prism.oregonstate.edu/
__ http://iri.columbia.edu/

.. ipython source code for this section
   we don't use this to avoid hitting the DAP server on every doc build.

   remote_data = xr.open_dataset(
       'http://iridl.ldeo.columbia.edu/SOURCES/.OSU/.PRISM/.monthly/dods',
       decode_times=False)
   tmax = remote_data.tmax[:500, ::3, ::3]
   tmax

   @savefig opendap-prism-tmax.png
   tmax[0].plot()

.. ipython::
    :verbatim:

    In [3]: remote_data = xr.open_dataset(
       ...:     "http://iridl.ldeo.columbia.edu/SOURCES/.OSU/.PRISM/.monthly/dods",
       ...:     decode_times=False,
       ...: )

    In [4]: remote_data
    Out[4]:
    <xarray.Dataset>
    Dimensions:  (T: 1422, X: 1405, Y: 621)
    Coordinates:
      * X        (X) float32 -125.0 -124.958 -124.917 -124.875 -124.833 -124.792 -124.75 ...
      * T        (T) float32 -779.5 -778.5 -777.5 -776.5 -775.5 -774.5 -773.5 -772.5 -771.5 ...
      * Y        (Y) float32 49.9167 49.875 49.8333 49.7917 49.75 49.7083 49.6667 49.625 ...
    Data variables:
        ppt      (T, Y, X) float64 ...
        tdmean   (T, Y, X) float64 ...
        tmax     (T, Y, X) float64 ...
        tmin     (T, Y, X) float64 ...
    Attributes:
        Conventions: IRIDL
        expires: 1375315200

.. TODO: update this example to show off decode_cf?

.. note::

    Like many real-world datasets, this dataset does not entirely follow
    `CF conventions`_. Unexpected formats will usually cause xarray's automatic
    decoding to fail. The way to work around this is to either set
    ``decode_cf=False`` in ``open_dataset`` to turn off all use of CF
    conventions, or by only disabling the troublesome parser.
    In this case, we set ``decode_times=False`` because the time axis here
    provides the calendar attribute in a format that xarray does not expect
    (the integer ``360`` instead of a string like ``'360_day'``).

We can select and slice this data any number of times, and nothing is loaded
over the network until we look at particular values:

.. ipython::
    :verbatim:

    In [4]: tmax = remote_data["tmax"][:500, ::3, ::3]

    In [5]: tmax
    Out[5]:
    <xarray.DataArray 'tmax' (T: 500, Y: 207, X: 469)>
    [48541500 values with dtype=float64]
    Coordinates:
      * Y        (Y) float32 49.9167 49.7917 49.6667 49.5417 49.4167 49.2917 ...
      * X        (X) float32 -125.0 -124.875 -124.75 -124.625 -124.5 -124.375 ...
      * T        (T) float32 -779.5 -778.5 -777.5 -776.5 -775.5 -774.5 -773.5 ...
    Attributes:
        pointwidth: 120
        standard_name: air_temperature
        units: Celsius_scale
        expires: 1443657600

    # the data is downloaded automatically when we make the plot
    In [6]: tmax[0].plot()

.. image:: ../_static/opendap-prism-tmax.png

Some servers require authentication before we can access the data. For this
purpose we can explicitly create a :py:class:`backends.PydapDataStore`
and pass in a `Requests`__ session object. For example for
HTTP Basic authentication::

    import xarray as xr
    import requests

    session = requests.Session()
    session.auth = ('username', 'password')

    store = xr.backends.PydapDataStore.open('http://example.com/data',
                                            session=session)
    ds = xr.open_dataset(store)

`Pydap's cas module`__ has functions that generate custom sessions for
servers that use CAS single sign-on. For example, to connect to servers
that require NASA's URS authentication::

  import xarray as xr
  from pydata.cas.urs import setup_session

  ds_url = 'https://gpm1.gesdisc.eosdis.nasa.gov/opendap/hyrax/example.nc'

  session = setup_session('username', 'password', check_url=ds_url)
  store = xr.backends.PydapDataStore.open(ds_url, session=session)

  ds = xr.open_dataset(store)

__ http://docs.python-requests.org
__ http://pydap.readthedocs.io/en/latest/client.html#authentication

.. _io.pickle:

Pickle
------

The simplest way to serialize an xarray object is to use Python's built-in pickle
module:

.. ipython:: python

    import pickle

    # use the highest protocol (-1) because it is way faster than the default
    # text based pickle format
    pkl = pickle.dumps(ds, protocol=-1)

    pickle.loads(pkl)

Pickling is important because it doesn't require any external libraries
and lets you use xarray objects with Python modules like
:py:mod:`multiprocessing` or :ref:`Dask <dask>`. However, pickling is
**not recommended for long-term storage**.

Restoring a pickle requires that the internal structure of the types for the
pickled data remain unchanged. Because the internal design of xarray is still
being refined, we make no guarantees (at this point) that objects pickled with
this version of xarray will work in future versions.

.. note::

  When pickling an object opened from a NetCDF file, the pickle file will
  contain a reference to the file on disk. If you want to store the actual
  array values, load it into memory first with :py:meth:`Dataset.load`
  or :py:meth:`Dataset.compute`.

.. _dictionary io:

Dictionary
----------

We can convert a ``Dataset`` (or a ``DataArray``) to a dict using
:py:meth:`Dataset.to_dict`:

.. ipython:: python

    d = ds.to_dict()
    d

We can create a new xarray object from a dict using
:py:meth:`Dataset.from_dict`:

.. ipython:: python

    ds_dict = xr.Dataset.from_dict(d)
    ds_dict

Dictionary support allows for flexible use of xarray objects. It doesn't
require external libraries and dicts can easily be pickled, or converted to
json, or geojson. All the values are converted to lists, so dicts might
be quite large.

To export just the dataset schema without the data itself, use the
``data=False`` option:

.. ipython:: python

    ds.to_dict(data=False)

This can be useful for generating indices of dataset contents to expose to
search indices or other automated data discovery tools.

.. ipython:: python
    :suppress:

    import os

    os.remove("saved_on_disk.nc")

.. _io.rasterio:

Rasterio
--------

GeoTIFFs and other gridded raster datasets can be opened using `rasterio`_, if
rasterio is installed. Here is an example of how to use
:py:func:`open_rasterio` to read one of rasterio's `test files`_:

.. deprecated:: 0.20.0

        Deprecated in favor of rioxarray.
        For information about transitioning, see:
        https://corteva.github.io/rioxarray/stable/getting_started/getting_started.html

.. ipython::
    :verbatim:

    In [7]: rio = xr.open_rasterio("RGB.byte.tif")

    In [8]: rio
    Out[8]:
    <xarray.DataArray (band: 3, y: 718, x: 791)>
    [1703814 values with dtype=uint8]
    Coordinates:
      * band     (band) int64 1 2 3
      * y        (y) float64 2.827e+06 2.826e+06 2.826e+06 2.826e+06 2.826e+06 ...
      * x        (x) float64 1.021e+05 1.024e+05 1.027e+05 1.03e+05 1.033e+05 ...
    Attributes:
        res:        (300.0379266750948, 300.041782729805)
        transform:  (300.0379266750948, 0.0, 101985.0, 0.0, -300.041782729805, 28...
        is_tiled:   0
        crs:        +init=epsg:32618


The ``x`` and ``y`` coordinates are generated out of the file's metadata
(``bounds``, ``width``, ``height``), and they can be understood as cartesian
coordinates defined in the file's projection provided by the ``crs`` attribute.
``crs`` is a PROJ4 string which can be parsed by e.g. `pyproj`_ or rasterio.
See :ref:`/examples/visualization_gallery.ipynb#Parsing-rasterio-geocoordinates`
for an example of how to convert these to longitudes and latitudes.


Additionally, you can use `rioxarray`_ for reading in GeoTiff, netCDF or other
GDAL readable raster data using `rasterio`_ as well as for exporting to a geoTIFF.
`rioxarray`_ can also handle geospatial related tasks such as re-projecting and clipping.

.. ipython::
    :verbatim:

    In [1]: import rioxarray

    In [2]: rds = rioxarray.open_rasterio("RGB.byte.tif")

    In [3]: rds
    Out[3]:
    <xarray.DataArray (band: 3, y: 718, x: 791)>
    [1703814 values with dtype=uint8]
    Coordinates:
      * band         (band) int64 1 2 3
      * y            (y) float64 2.827e+06 2.826e+06 ... 2.612e+06 2.612e+06
      * x            (x) float64 1.021e+05 1.024e+05 ... 3.389e+05 3.392e+05
        spatial_ref  int64 0
    Attributes:
        STATISTICS_MAXIMUM:  255
        STATISTICS_MEAN:     29.947726688477
        STATISTICS_MINIMUM:  0
        STATISTICS_STDDEV:   52.340921626611
        transform:           (300.0379266750948, 0.0, 101985.0, 0.0, -300.0417827...
        _FillValue:          0.0
        scale_factor:        1.0
        add_offset:          0.0
        grid_mapping:        spatial_ref

    In [4]: rds.rio.crs
    Out[4]: CRS.from_epsg(32618)

    In [5]: rds4326 = rds.rio.reproject("epsg:4326")

    In [6]: rds4326.rio.crs
    Out[6]: CRS.from_epsg(4326)

    In [7]: rds4326.rio.to_raster("RGB.byte.4326.tif")


.. _rasterio: https://rasterio.readthedocs.io/en/latest/
.. _rioxarray: https://corteva.github.io/rioxarray/stable/
.. _test files: https://github.com/mapbox/rasterio/blob/master/tests/data/RGB.byte.tif
.. _pyproj: https://github.com/pyproj4/pyproj

.. _io.zarr:

Zarr
----

`Zarr`_ is a Python package that provides an implementation of chunked, compressed,
N-dimensional arrays.
Zarr has the ability to store arrays in a range of ways, including in memory,
in files, and in cloud-based object storage such as `Amazon S3`_ and
`Google Cloud Storage`_.
Xarray's Zarr backend allows xarray to leverage these capabilities, including
the ability to store and analyze datasets far too large fit onto disk
(particularly :ref:`in combination with dask <dask>`).

Xarray can't open just any zarr dataset, because xarray requires special
metadata (attributes) describing the dataset dimensions and coordinates.
At this time, xarray can only open zarr datasets that have been written by
xarray. For implementation details, see :ref:`zarr_encoding`.

To write a dataset with zarr, we use the :py:meth:`Dataset.to_zarr` method.

To write to a local directory, we pass a path to a directory:

.. ipython:: python
    :suppress:

    ! rm -rf path/to/directory.zarr

.. ipython:: python

    ds = xr.Dataset(
        {"foo": (("x", "y"), np.random.rand(4, 5))},
        coords={
            "x": [10, 20, 30, 40],
            "y": pd.date_range("2000-01-01", periods=5),
            "z": ("x", list("abcd")),
        },
    )
    ds.to_zarr("path/to/directory.zarr")

(The suffix ``.zarr`` is optional--just a reminder that a zarr store lives
there.) If the directory does not exist, it will be created. If a zarr
store is already present at that path, an error will be raised, preventing it
from being overwritten. To override this behavior and overwrite an existing
store, add ``mode='w'`` when invoking :py:meth:`~Dataset.to_zarr`.

To store variable length strings, convert them to object arrays first with
``dtype=object``.

To read back a zarr dataset that has been created this way, we use the
:py:func:`open_zarr` method:

.. ipython:: python

    ds_zarr = xr.open_zarr("path/to/directory.zarr")
    ds_zarr

Cloud Storage Buckets
~~~~~~~~~~~~~~~~~~~~~

It is possible to read and write xarray datasets directly from / to cloud
storage buckets using zarr. This example uses the `gcsfs`_ package to provide
an interface to `Google Cloud Storage`_.

From v0.16.2: general `fsspec`_ URLs are parsed and the store set up for you
automatically when reading, such that you can open a dataset in a single
call. You should include any arguments to the storage backend as the
key ``storage_options``, part of ``backend_kwargs``.

.. code:: python

    ds_gcs = xr.open_dataset(
        "gcs://<bucket-name>/path.zarr",
        backend_kwargs={
            "storage_options": {"project": "<project-name>", "token": None}
        },
        engine="zarr",
    )


This also works with ``open_mfdataset``, allowing you to pass a list of paths or
a URL to be interpreted as a glob string.

For older versions, and for writing, you must explicitly set up a ``MutableMapping``
instance and pass this, as follows:

.. code:: python

    import gcsfs

    fs = gcsfs.GCSFileSystem(project="<project-name>", token=None)
    gcsmap = gcsfs.mapping.GCSMap("<bucket-name>", gcs=fs, check=True, create=False)
    # write to the bucket
    ds.to_zarr(store=gcsmap)
    # read it back
    ds_gcs = xr.open_zarr(gcsmap)

(or use the utility function ``fsspec.get_mapper()``).

.. _fsspec: https://filesystem-spec.readthedocs.io/en/latest/
.. _Zarr: http://zarr.readthedocs.io/
.. _Amazon S3: https://aws.amazon.com/s3/
.. _Google Cloud Storage: https://cloud.google.com/storage/
.. _gcsfs: https://github.com/dask/gcsfs

Zarr Compressors and Filters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are many different options for compression and filtering possible with
zarr. These are described in the
`zarr documentation <http://zarr.readthedocs.io/en/stable/tutorial.html#compressors>`_.
These options can be passed to the ``to_zarr`` method as variable encoding.
For example:

.. ipython:: python
    :suppress:

    ! rm -rf foo.zarr

.. ipython:: python

    import zarr

    compressor = zarr.Blosc(cname="zstd", clevel=3, shuffle=2)
    ds.to_zarr("foo.zarr", encoding={"foo": {"compressor": compressor}})

.. note::

    Not all native zarr compression and filtering options have been tested with
    xarray.

.. _io.zarr.consolidated_metadata:

Consolidated Metadata
~~~~~~~~~~~~~~~~~~~~~

Xarray needs to read all of the zarr metadata when it opens a dataset.
In some storage mediums, such as with cloud object storage (e.g. amazon S3),
this can introduce significant overhead, because two separate HTTP calls to the
object store must be made for each variable in the dataset.
As of xarray version 0.18, xarray by default uses a feature called
*consolidated metadata*, storing all metadata for the entire dataset with a
single key (by default called ``.zmetadata``). This typically drastically speeds
up opening the store. (For more information on this feature, consult the
`zarr docs <https://zarr.readthedocs.io/en/latest/tutorial.html#consolidating-metadata>`_.)

By default, xarray writes consolidated metadata and attempts to read stores
with consolidated metadata, falling back to use non-consolidated metadata for
reads. Because this fall-back option is so much slower, xarray issues a
``RuntimeWarning`` with guidance when reading with consolidated metadata fails:

    Failed to open Zarr store with consolidated metadata, falling back to try
    reading non-consolidated metadata. This is typically much slower for
    opening a dataset. To silence this warning, consider:

    1. Consolidating metadata in this existing store with
       :py:func:`zarr.consolidate_metadata`.
    2. Explicitly setting ``consolidated=False``, to avoid trying to read
       consolidate metadata.
    3. Explicitly setting ``consolidated=True``, to raise an error in this case
       instead of falling back to try reading non-consolidated metadata.

.. _io.zarr.appending:

Appending to existing Zarr stores
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Xarray supports several ways of incrementally writing variables to a Zarr
store. These options are useful for scenarios when it is infeasible or
undesirable to write your entire dataset at once.

.. tip::

    If you can load all of your data into a single ``Dataset`` using dask, a
    single call to ``to_zarr()`` will write all of your data in parallel.

.. warning::

    Alignment of coordinates is currently not checked when modifying an
    existing Zarr store. It is up to the user to ensure that coordinates are
    consistent.

To add or overwrite entire variables, simply call :py:meth:`~Dataset.to_zarr`
with ``mode='a'`` on a Dataset containing the new variables, passing in an
existing Zarr store or path to a Zarr store.

To resize and then append values along an existing dimension in a store, set
``append_dim``. This is a good option if data always arives in a particular
order, e.g., for time-stepping a simulation:

.. ipython:: python
    :suppress:

    ! rm -rf path/to/directory.zarr

.. ipython:: python

    ds1 = xr.Dataset(
        {"foo": (("x", "y", "t"), np.random.rand(4, 5, 2))},
        coords={
            "x": [10, 20, 30, 40],
            "y": [1, 2, 3, 4, 5],
            "t": pd.date_range("2001-01-01", periods=2),
        },
    )
    ds1.to_zarr("path/to/directory.zarr")
    ds2 = xr.Dataset(
        {"foo": (("x", "y", "t"), np.random.rand(4, 5, 2))},
        coords={
            "x": [10, 20, 30, 40],
            "y": [1, 2, 3, 4, 5],
            "t": pd.date_range("2001-01-03", periods=2),
        },
    )
    ds2.to_zarr("path/to/directory.zarr", append_dim="t")

Finally, you can use ``region`` to write to limited regions of existing arrays
in an existing Zarr store. This is a good option for writing data in parallel
from independent processes.

To scale this up to writing large datasets, the first step is creating an
initial Zarr store without writing all of its array data. This can be done by
first creating a ``Dataset`` with dummy values stored in :ref:`dask <dask>`,
and then calling ``to_zarr`` with ``compute=False`` to write only metadata
(including ``attrs``) to Zarr:

.. ipython:: python
    :suppress:

    ! rm -rf path/to/directory.zarr

.. ipython:: python

    import dask.array

    # The values of this dask array are entirely irrelevant; only the dtype,
    # shape and chunks are used
    dummies = dask.array.zeros(30, chunks=10)
    ds = xr.Dataset({"foo": ("x", dummies)})
    path = "path/to/directory.zarr"
    # Now we write the metadata without computing any array values
    ds.to_zarr(path, compute=False)

Now, a Zarr store with the correct variable shapes and attributes exists that
can be filled out by subsequent calls to ``to_zarr``. The ``region`` provides a
mapping from dimension names to Python ``slice`` objects indicating where the
data should be written (in index space, not coordinate space), e.g.,

.. ipython:: python

    # For convenience, we'll slice a single dataset, but in the real use-case
    # we would create them separately possibly even from separate processes.
    ds = xr.Dataset({"foo": ("x", np.arange(30))})
    ds.isel(x=slice(0, 10)).to_zarr(path, region={"x": slice(0, 10)})
    ds.isel(x=slice(10, 20)).to_zarr(path, region={"x": slice(10, 20)})
    ds.isel(x=slice(20, 30)).to_zarr(path, region={"x": slice(20, 30)})

Concurrent writes with ``region`` are safe as long as they modify distinct
chunks in the underlying Zarr arrays (or use an appropriate ``lock``).

As a safety check to make it harder to inadvertently override existing values,
if you set ``region`` then *all* variables included in a Dataset must have
dimensions included in ``region``. Other variables (typically coordinates)
need to be explicitly dropped and/or written in a separate calls to ``to_zarr``
with ``mode='a'``.

.. _io.cfgrib:

.. ipython:: python
    :suppress:

    import shutil

    shutil.rmtree("foo.zarr")
    shutil.rmtree("path/to/directory.zarr")

GRIB format via cfgrib
----------------------

Xarray supports reading GRIB files via ECMWF cfgrib_ python driver,
if it is installed. To open a GRIB file supply ``engine='cfgrib'``
to :py:func:`open_dataset`:

.. ipython::
    :verbatim:

    In [1]: ds_grib = xr.open_dataset("example.grib", engine="cfgrib")

We recommend installing cfgrib via conda::

    conda install -c conda-forge cfgrib

.. _cfgrib: https://github.com/ecmwf/cfgrib

.. _io.pynio:

Formats supported by PyNIO
--------------------------

Xarray can also read GRIB, HDF4 and other file formats supported by PyNIO_,
if PyNIO is installed. To use PyNIO to read such files, supply
``engine='pynio'`` to :py:func:`open_dataset`.

We recommend installing PyNIO via conda::

    conda install -c conda-forge pynio

.. warning::

    PyNIO is no longer actively maintained and conflicts with netcdf4 > 1.5.3.
    The PyNIO backend may be moved outside of xarray in the future.

.. _PyNIO: https://www.pyngl.ucar.edu/Nio.shtml

.. _io.PseudoNetCDF:

Formats supported by PseudoNetCDF
---------------------------------

Xarray can also read CAMx, BPCH, ARL PACKED BIT, and many other file
formats supported by PseudoNetCDF_, if PseudoNetCDF is installed.
PseudoNetCDF can also provide Climate Forecasting Conventions to
CMAQ files. In addition, PseudoNetCDF can automatically register custom
readers that subclass PseudoNetCDF.PseudoNetCDFFile. PseudoNetCDF can
identify readers either heuristically, or by a format specified via a key in
`backend_kwargs`.

To use PseudoNetCDF to read such files, supply
``engine='pseudonetcdf'`` to :py:func:`open_dataset`.

Add ``backend_kwargs={'format': '<format name>'}`` where `<format name>`
options are listed on the PseudoNetCDF page.

.. _PseudoNetCDF: http://github.com/barronh/PseudoNetCDF


CSV and other formats supported by pandas
-----------------------------------------

For more options (tabular formats and CSV files in particular), consider
exporting your objects to pandas and using its broad range of `IO tools`_.
For CSV files, one might also consider `xarray_extras`_.

.. _xarray_extras: https://xarray-extras.readthedocs.io/en/latest/api/csv.html

.. _IO tools: http://pandas.pydata.org/pandas-docs/stable/io.html


Third party libraries
---------------------

More formats are supported by extension libraries:

- `xarray-mongodb <https://xarray-mongodb.readthedocs.io/en/latest/>`_: Store xarray objects on MongoDB
.. currentmodule:: xarray
.. _plotting:

Plotting
========

Introduction
------------

Labeled data enables expressive computations. These same
labels can also be used to easily create informative plots.

Xarray's plotting capabilities are centered around
:py:class:`DataArray` objects.
To plot :py:class:`Dataset` objects
simply access the relevant DataArrays, i.e. ``dset['var1']``.
Dataset specific plotting routines are also available (see :ref:`plot-dataset`).
Here we focus mostly on arrays 2d or larger. If your data fits
nicely into a pandas DataFrame then you're better off using one of the more
developed tools there.

Xarray plotting functionality is a thin wrapper around the popular
`matplotlib <http://matplotlib.org/>`_ library.
Matplotlib syntax and function names were copied as much as possible, which
makes for an easy transition between the two.
Matplotlib must be installed before xarray can plot.

To use xarray's plotting capabilities with time coordinates containing
``cftime.datetime`` objects
`nc-time-axis <https://github.com/SciTools/nc-time-axis>`_ v1.2.0 or later
needs to be installed.

For more extensive plotting applications consider the following projects:

- `Seaborn <http://seaborn.pydata.org/>`_: "provides
  a high-level interface for drawing attractive statistical graphics."
  Integrates well with pandas.

- `HoloViews <http://holoviews.org/>`_
  and `GeoViews <https://geoviews.org/>`_: "Composable, declarative
  data structures for building even complex visualizations easily." Includes
  native support for xarray objects.

- `hvplot <https://hvplot.pyviz.org/>`_: ``hvplot`` makes it very easy to produce
  dynamic plots (backed by ``Holoviews`` or ``Geoviews``) by adding a ``hvplot``
  accessor to DataArrays.

- `Cartopy <http://scitools.org.uk/cartopy/>`_: Provides cartographic
  tools.

Imports
~~~~~~~

.. ipython:: python
    :suppress:

    # Use defaults so we don't get gridlines in generated docs
    import matplotlib as mpl

    mpl.rcdefaults()

The following imports are necessary for all of the examples.

.. ipython:: python

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import xarray as xr

For these examples we'll use the North American air temperature dataset.

.. ipython:: python

    airtemps = xr.tutorial.open_dataset("air_temperature")
    airtemps

    # Convert to celsius
    air = airtemps.air - 273.15

    # copy attributes to get nice figure labels and change Kelvin to Celsius
    air.attrs = airtemps.air.attrs
    air.attrs["units"] = "deg C"

.. note::
   Until :issue:`1614` is solved, you might need to copy over the metadata in ``attrs`` to get informative figure labels (as was done above).


DataArrays
----------

One Dimension
~~~~~~~~~~~~~

================
 Simple Example
================

The simplest way to make a plot is to call the :py:func:`DataArray.plot()` method.

.. ipython:: python
    :okwarning:

    air1d = air.isel(lat=10, lon=10)

    @savefig plotting_1d_simple.png width=4in
    air1d.plot()

Xarray uses the coordinate name along with metadata ``attrs.long_name``, ``attrs.standard_name``, ``DataArray.name`` and ``attrs.units`` (if available) to label the axes. The names ``long_name``, ``standard_name`` and ``units`` are copied from the `CF-conventions spec <http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/build/ch03s03.html>`_. When choosing names, the order of precedence is ``long_name``, ``standard_name`` and finally ``DataArray.name``. The y-axis label in the above plot was constructed from the ``long_name`` and ``units`` attributes of ``air1d``.

.. ipython:: python

    air1d.attrs

======================
 Additional Arguments
======================

Additional arguments are passed directly to the matplotlib function which
does the work.
For example, :py:func:`xarray.plot.line` calls
matplotlib.pyplot.plot_ passing in the index and the array values as x and y, respectively.
So to make a line plot with blue triangles a matplotlib format string
can be used:

.. _matplotlib.pyplot.plot: http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot

.. ipython:: python
    :okwarning:

    @savefig plotting_1d_additional_args.png width=4in
    air1d[:200].plot.line("b-^")

.. note::
    Not all xarray plotting methods support passing positional arguments
    to the wrapped matplotlib functions, but they do all
    support keyword arguments.

Keyword arguments work the same way, and are more explicit.

.. ipython:: python
    :okwarning:

    @savefig plotting_example_sin3.png width=4in
    air1d[:200].plot.line(color="purple", marker="o")

=========================
 Adding to Existing Axis
=========================

To add the plot to an existing axis pass in the axis as a keyword argument
``ax``. This works for all xarray plotting methods.
In this example ``axes`` is an array consisting of the left and right
axes created by ``plt.subplots``.

.. ipython:: python
    :okwarning:

    fig, axes = plt.subplots(ncols=2)

    axes

    air1d.plot(ax=axes[0])
    air1d.plot.hist(ax=axes[1])

    plt.tight_layout()

    @savefig plotting_example_existing_axes.png width=6in
    plt.draw()

On the right is a histogram created by :py:func:`xarray.plot.hist`.

.. _plotting.figsize:

=============================
 Controlling the figure size
=============================

You can pass a ``figsize`` argument to all xarray's plotting methods to
control the figure size. For convenience, xarray's plotting methods also
support the ``aspect`` and ``size`` arguments which control the size of the
resulting image via the formula ``figsize = (aspect * size, size)``:

.. ipython:: python
    :okwarning:

    air1d.plot(aspect=2, size=3)
    @savefig plotting_example_size_and_aspect.png
    plt.tight_layout()

.. ipython:: python
    :suppress:

    # create a dummy figure so sphinx plots everything below normally
    plt.figure()

This feature also works with :ref:`plotting.faceting`. For facet plots,
``size`` and ``aspect`` refer to a single panel (so that ``aspect * size``
gives the width of each facet in inches), while ``figsize`` refers to the
entire figure (as for matplotlib's ``figsize`` argument).

.. note::

    If ``figsize`` or ``size`` are used, a new figure is created,
    so this is mutually exclusive with the ``ax`` argument.

.. note::

    The convention used by xarray (``figsize = (aspect * size, size)``) is
    borrowed from seaborn: it is therefore `not equivalent to matplotlib's`_.

.. _not equivalent to matplotlib's: https://github.com/mwaskom/seaborn/issues/746


.. _plotting.multiplelines:

=========================
 Determine x-axis values
=========================

Per default dimension coordinates are used for the x-axis (here the time coordinates).
However, you can also use non-dimension coordinates, MultiIndex levels, and dimensions
without coordinates along the x-axis. To illustrate this, let's calculate a 'decimal day' (epoch)
from the time and assign it as a non-dimension coordinate:

.. ipython:: python
    :okwarning:

    decimal_day = (air1d.time - air1d.time[0]) / pd.Timedelta("1d")
    air1d_multi = air1d.assign_coords(decimal_day=("time", decimal_day.data))
    air1d_multi

To use ``'decimal_day'`` as x coordinate it must be explicitly specified:

.. ipython:: python
    :okwarning:

    air1d_multi.plot(x="decimal_day")

Creating a new MultiIndex named ``'date'`` from ``'time'`` and ``'decimal_day'``,
it is also possible to use a MultiIndex level as x-axis:

.. ipython:: python
    :okwarning:

    air1d_multi = air1d_multi.set_index(date=("time", "decimal_day"))
    air1d_multi.plot(x="decimal_day")

Finally, if a dataset does not have any coordinates it enumerates all data points:

.. ipython:: python
    :okwarning:

    air1d_multi = air1d_multi.drop("date")
    air1d_multi.plot()

The same applies to 2D plots below.

====================================================
 Multiple lines showing variation along a dimension
====================================================

It is possible to make line plots of two-dimensional data by calling :py:func:`xarray.plot.line`
with appropriate arguments. Consider the 3D variable ``air`` defined above. We can use line
plots to check the variation of air temperature at three different latitudes along a longitude line:

.. ipython:: python
    :okwarning:

    @savefig plotting_example_multiple_lines_x_kwarg.png
    air.isel(lon=10, lat=[19, 21, 22]).plot.line(x="time")

It is required to explicitly specify either

1. ``x``: the dimension to be used for the x-axis, or
2. ``hue``: the dimension you want to represent by multiple lines.

Thus, we could have made the previous plot by specifying ``hue='lat'`` instead of ``x='time'``.
If required, the automatic legend can be turned off using ``add_legend=False``. Alternatively,
``hue`` can be passed directly to :py:func:`xarray.plot.line` as `air.isel(lon=10, lat=[19,21,22]).plot.line(hue='lat')`.


========================
 Dimension along y-axis
========================

It is also possible to make line plots such that the data are on the x-axis and a dimension is on the y-axis. This can be done by specifying the appropriate ``y`` keyword argument.

.. ipython:: python
    :okwarning:

    @savefig plotting_example_xy_kwarg.png
    air.isel(time=10, lon=[10, 11]).plot(y="lat", hue="lon")

============
 Step plots
============

As an alternative, also a step plot similar to matplotlib's ``plt.step`` can be
made using 1D data.

.. ipython:: python
    :okwarning:

    @savefig plotting_example_step.png width=4in
    air1d[:20].plot.step(where="mid")

The argument ``where`` defines where the steps should be placed, options are
``'pre'`` (default), ``'post'``, and ``'mid'``. This is particularly handy
when plotting data grouped with :py:meth:`Dataset.groupby_bins`.

.. ipython:: python
    :okwarning:

    air_grp = air.mean(["time", "lon"]).groupby_bins("lat", [0, 23.5, 66.5, 90])
    air_mean = air_grp.mean()
    air_std = air_grp.std()
    air_mean.plot.step()
    (air_mean + air_std).plot.step(ls=":")
    (air_mean - air_std).plot.step(ls=":")
    plt.ylim(-20, 30)
    @savefig plotting_example_step_groupby.png width=4in
    plt.title("Zonal mean temperature")

In this case, the actual boundaries of the bins are used and the ``where`` argument
is ignored.


Other axes kwargs
~~~~~~~~~~~~~~~~~


The keyword arguments ``xincrease`` and ``yincrease`` let you control the axes direction.

.. ipython:: python
    :okwarning:

    @savefig plotting_example_xincrease_yincrease_kwarg.png
    air.isel(time=10, lon=[10, 11]).plot.line(
        y="lat", hue="lon", xincrease=False, yincrease=False
    )

In addition, one can use ``xscale, yscale`` to set axes scaling; ``xticks, yticks`` to set axes ticks and ``xlim, ylim`` to set axes limits. These accept the same values as the matplotlib methods ``Axes.set_(x,y)scale()``, ``Axes.set_(x,y)ticks()``, ``Axes.set_(x,y)lim()`` respectively.


Two Dimensions
~~~~~~~~~~~~~~

================
 Simple Example
================

The default method :py:meth:`DataArray.plot` calls :py:func:`xarray.plot.pcolormesh` by default when the data is two-dimensional.

.. ipython:: python
    :okwarning:

    air2d = air.isel(time=500)

    @savefig 2d_simple.png width=4in
    air2d.plot()

All 2d plots in xarray allow the use of the keyword arguments ``yincrease``
and ``xincrease``.

.. ipython:: python
    :okwarning:

    @savefig 2d_simple_yincrease.png width=4in
    air2d.plot(yincrease=False)

.. note::

    We use :py:func:`xarray.plot.pcolormesh` as the default two-dimensional plot
    method because it is more flexible than :py:func:`xarray.plot.imshow`.
    However, for large arrays, ``imshow`` can be much faster than ``pcolormesh``.
    If speed is important to you and you are plotting a regular mesh, consider
    using ``imshow``.

================
 Missing Values
================

Xarray plots data with :ref:`missing_values`.

.. ipython:: python
    :okwarning:

    bad_air2d = air2d.copy()

    bad_air2d[dict(lat=slice(0, 10), lon=slice(0, 25))] = np.nan

    @savefig plotting_missing_values.png width=4in
    bad_air2d.plot()

========================
 Nonuniform Coordinates
========================

It's not necessary for the coordinates to be evenly spaced. Both
:py:func:`xarray.plot.pcolormesh` (default) and :py:func:`xarray.plot.contourf` can
produce plots with nonuniform coordinates.

.. ipython:: python
    :okwarning:

    b = air2d.copy()
    # Apply a nonlinear transformation to one of the coords
    b.coords["lat"] = np.log(b.coords["lat"])

    @savefig plotting_nonuniform_coords.png width=4in
    b.plot()

====================
 Other types of plot
====================

There are several other options for plotting 2D data.

Contour plot using :py:meth:`DataArray.plot.contour()`

.. ipython:: python
    :okwarning:

    @savefig plotting_contour.png width=4in
    air2d.plot.contour()

Filled contour plot using :py:meth:`DataArray.plot.contourf()`

.. ipython:: python
    :okwarning:

    @savefig plotting_contourf.png width=4in
    air2d.plot.contourf()

Surface plot using :py:meth:`DataArray.plot.surface()`

.. ipython:: python
    :okwarning:

    @savefig plotting_surface.png width=4in
    # transpose just to make the example look a bit nicer
    air2d.T.plot.surface()

====================
 Calling Matplotlib
====================

Since this is a thin wrapper around matplotlib, all the functionality of
matplotlib is available.

.. ipython:: python
    :okwarning:

    air2d.plot(cmap=plt.cm.Blues)
    plt.title("These colors prove North America\nhas fallen in the ocean")
    plt.ylabel("latitude")
    plt.xlabel("longitude")
    plt.tight_layout()

    @savefig plotting_2d_call_matplotlib.png width=4in
    plt.draw()

.. note::

    Xarray methods update label information and generally play around with the
    axes. So any kind of updates to the plot
    should be done *after* the call to the xarray's plot.
    In the example below, ``plt.xlabel`` effectively does nothing, since
    ``d_ylog.plot()`` updates the xlabel.

    .. ipython:: python
        :okwarning:

        plt.xlabel("Never gonna see this.")
        air2d.plot()

        @savefig plotting_2d_call_matplotlib2.png width=4in
        plt.draw()

===========
 Colormaps
===========

Xarray borrows logic from Seaborn to infer what kind of color map to use. For
example, consider the original data in Kelvins rather than Celsius:

.. ipython:: python
    :okwarning:

    @savefig plotting_kelvin.png width=4in
    airtemps.air.isel(time=0).plot()

The Celsius data contain 0, so a diverging color map was used. The
Kelvins do not have 0, so the default color map was used.

.. _robust-plotting:

========
 Robust
========

Outliers often have an extreme effect on the output of the plot.
Here we add two bad data points. This affects the color scale,
washing out the plot.

.. ipython:: python
    :okwarning:

    air_outliers = airtemps.air.isel(time=0).copy()
    air_outliers[0, 0] = 100
    air_outliers[-1, -1] = 400

    @savefig plotting_robust1.png width=4in
    air_outliers.plot()

This plot shows that we have outliers. The easy way to visualize
the data without the outliers is to pass the parameter
``robust=True``.
This will use the 2nd and 98th
percentiles of the data to compute the color limits.

.. ipython:: python
    :okwarning:

    @savefig plotting_robust2.png width=4in
    air_outliers.plot(robust=True)

Observe that the ranges of the color bar have changed. The arrows on the
color bar indicate
that the colors include data points outside the bounds.

====================
 Discrete Colormaps
====================

It is often useful, when visualizing 2d data, to use a discrete colormap,
rather than the default continuous colormaps that matplotlib uses. The
``levels`` keyword argument can be used to generate plots with discrete
colormaps. For example, to make a plot with 8 discrete color intervals:

.. ipython:: python
    :okwarning:

    @savefig plotting_discrete_levels.png width=4in
    air2d.plot(levels=8)

It is also possible to use a list of levels to specify the boundaries of the
discrete colormap:

.. ipython:: python
    :okwarning:

    @savefig plotting_listed_levels.png width=4in
    air2d.plot(levels=[0, 12, 18, 30])

You can also specify a list of discrete colors through the ``colors`` argument:

.. ipython:: python
    :okwarning:

    flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
    @savefig plotting_custom_colors_levels.png width=4in
    air2d.plot(levels=[0, 12, 18, 30], colors=flatui)

Finally, if you have `Seaborn <http://seaborn.pydata.org/>`_
installed, you can also specify a seaborn color palette to the ``cmap``
argument. Note that ``levels`` *must* be specified with seaborn color palettes
if using ``imshow`` or ``pcolormesh`` (but not with ``contour`` or ``contourf``,
since levels are chosen automatically).

.. ipython:: python
    :okwarning:

    @savefig plotting_seaborn_palette.png width=4in
    air2d.plot(levels=10, cmap="husl")
    plt.draw()

.. _plotting.faceting:

Faceting
~~~~~~~~

Faceting here refers to splitting an array along one or two dimensions and
plotting each group.
Xarray's basic plotting is useful for plotting two dimensional arrays. What
about three or four dimensional arrays? That's where facets become helpful.
The general approach to plotting here is called “small multiples”, where the same kind of plot is repeated multiple times, and the specific use of small multiples to display the same relationship conditioned on one ore more other variables is often called a “trellis plot”.

Consider the temperature data set. There are 4 observations per day for two
years which makes for 2920 values along the time dimension.
One way to visualize this data is to make a
separate plot for each time period.

The faceted dimension should not have too many values;
faceting on the time dimension will produce 2920 plots. That's
too much to be helpful. To handle this situation try performing
an operation that reduces the size of the data in some way. For example, we
could compute the average air temperature for each month and reduce the
size of this dimension from 2920 -> 12. A simpler way is
to just take a slice on that dimension.
So let's use a slice to pick 6 times throughout the first year.

.. ipython:: python

    t = air.isel(time=slice(0, 365 * 4, 250))
    t.coords

================
 Simple Example
================

The easiest way to create faceted plots is to pass in ``row`` or ``col``
arguments to the xarray plotting methods/functions. This returns a
:py:class:`xarray.plot.FacetGrid` object.

.. ipython:: python
    :okwarning:

    @savefig plot_facet_dataarray.png
    g_simple = t.plot(x="lon", y="lat", col="time", col_wrap=3)

Faceting also works for line plots.

.. ipython:: python
    :okwarning:

    @savefig plot_facet_dataarray_line.png
    g_simple_line = t.isel(lat=slice(0, None, 4)).plot(
        x="lon", hue="lat", col="time", col_wrap=3
    )

===============
 4 dimensional
===============

For 4 dimensional arrays we can use the rows and columns of the grids.
Here we create a 4 dimensional array by taking the original data and adding
a fixed amount. Now we can see how the temperature maps would compare if
one were much hotter.

.. ipython:: python
    :okwarning:

    t2 = t.isel(time=slice(0, 2))
    t4d = xr.concat([t2, t2 + 40], pd.Index(["normal", "hot"], name="fourth_dim"))
    # This is a 4d array
    t4d.coords

    @savefig plot_facet_4d.png
    t4d.plot(x="lon", y="lat", col="time", row="fourth_dim")

================
 Other features
================

Faceted plotting supports other arguments common to xarray 2d plots.

.. ipython:: python
    :suppress:

    plt.close("all")

.. ipython:: python
    :okwarning:

    hasoutliers = t.isel(time=slice(0, 5)).copy()
    hasoutliers[0, 0, 0] = -100
    hasoutliers[-1, -1, -1] = 400

    @savefig plot_facet_robust.png
    g = hasoutliers.plot.pcolormesh(
        "lon",
        "lat",
        col="time",
        col_wrap=3,
        robust=True,
        cmap="viridis",
        cbar_kwargs={"label": "this has outliers"},
    )

===================
 FacetGrid Objects
===================

The object returned, ``g`` in the above examples, is a :py:class:`~xarray.plot.FacetGrid` object
that links a :py:class:`DataArray` to a matplotlib figure with a particular structure.
This object can be used to control the behavior of the multiple plots.
It borrows an API and code from `Seaborn's FacetGrid
<http://seaborn.pydata.org/tutorial/axis_grids.html>`_.
The structure is contained within the ``axes`` and ``name_dicts``
attributes, both 2d NumPy object arrays.

.. ipython:: python

    g.axes

    g.name_dicts

It's possible to select the :py:class:`xarray.DataArray` or
:py:class:`xarray.Dataset` corresponding to the FacetGrid through the
``name_dicts``.

.. ipython:: python

    g.data.loc[g.name_dicts[0, 0]]

Here is an example of using the lower level API and then modifying the axes after
they have been plotted.

.. ipython:: python
    :okwarning:

    g = t.plot.imshow("lon", "lat", col="time", col_wrap=3, robust=True)

    for i, ax in enumerate(g.axes.flat):
        ax.set_title("Air Temperature %d" % i)

    bottomright = g.axes[-1, -1]
    bottomright.annotate("bottom right", (240, 40))

    @savefig plot_facet_iterator.png
    plt.draw()


:py:class:`~xarray.plot.FacetGrid` objects have methods that let you customize the automatically generated
axis labels, axis ticks and plot titles. See :py:meth:`~xarray.plot.FacetGrid.set_titles`,
:py:meth:`~xarray.plot.FacetGrid.set_xlabels`, :py:meth:`~xarray.plot.FacetGrid.set_ylabels` and
:py:meth:`~xarray.plot.FacetGrid.set_ticks` for more information.
Plotting functions can be applied to each subset of the data by calling :py:meth:`~xarray.plot.FacetGrid.map_dataarray` or to each subplot by calling :py:meth:`~xarray.plot.FacetGrid.map`.

TODO: add an example of using the ``map`` method to plot dataset variables
(e.g., with ``plt.quiver``).

.. _plot-dataset:

Datasets
--------

Xarray has limited support for plotting Dataset variables against each other.
Consider this dataset

.. ipython:: python

    ds = xr.tutorial.scatter_example_dataset()
    ds


Scatter
~~~~~~~

Suppose we want to scatter ``A`` against ``B``

.. ipython:: python
    :okwarning:

    @savefig ds_simple_scatter.png
    ds.plot.scatter(x="A", y="B")

The ``hue`` kwarg lets you vary the color by variable value

.. ipython:: python
    :okwarning:

    @savefig ds_hue_scatter.png
    ds.plot.scatter(x="A", y="B", hue="w")

When ``hue`` is specified, a colorbar is added for numeric ``hue`` DataArrays by
default and a legend is added for non-numeric ``hue`` DataArrays (as above).
You can force a legend instead of a colorbar by setting ``hue_style='discrete'``.
Additionally, the boolean kwarg ``add_guide`` can be used to prevent the display of a legend or colorbar (as appropriate).

.. ipython:: python
    :okwarning:

    ds = ds.assign(w=[1, 2, 3, 5])
    @savefig ds_discrete_legend_hue_scatter.png
    ds.plot.scatter(x="A", y="B", hue="w", hue_style="discrete")

The ``markersize`` kwarg lets you vary the point's size by variable value. You can additionally pass ``size_norm`` to control how the variable's values are mapped to point sizes.

.. ipython:: python
    :okwarning:

    @savefig ds_hue_size_scatter.png
    ds.plot.scatter(x="A", y="B", hue="z", hue_style="discrete", markersize="z")

Faceting is also possible

.. ipython:: python
    :okwarning:

    @savefig ds_facet_scatter.png
    ds.plot.scatter(x="A", y="B", col="x", row="z", hue="w", hue_style="discrete")


For more advanced scatter plots, we recommend converting the relevant data variables to a pandas DataFrame and using the extensive plotting capabilities of ``seaborn``.

Quiver
~~~~~~

Visualizing vector fields is supported with quiver plots:

.. ipython:: python
    :okwarning:

    @savefig ds_simple_quiver.png
    ds.isel(w=1, z=1).plot.quiver(x="x", y="y", u="A", v="B")


where ``u`` and ``v`` denote the x and y direction components of the arrow vectors. Again, faceting is also possible:

.. ipython:: python
    :okwarning:

    @savefig ds_facet_quiver.png
    ds.plot.quiver(x="x", y="y", u="A", v="B", col="w", row="z", scale=4)

``scale`` is required for faceted quiver plots. The scale determines the number of data units per arrow length unit, i.e. a smaller scale parameter makes the arrow longer.

Streamplot
~~~~~~~~~~

Visualizing vector fields is also supported with streamline plots:

.. ipython:: python
    :okwarning:

    @savefig ds_simple_streamplot.png
    ds.isel(w=1, z=1).plot.streamplot(x="x", y="y", u="A", v="B")


where ``u`` and ``v`` denote the x and y direction components of the vectors tangent to the streamlines. Again, faceting is also possible:

.. ipython:: python
    :okwarning:

    @savefig ds_facet_streamplot.png
    ds.plot.streamplot(x="x", y="y", u="A", v="B", col="w", row="z")

.. _plot-maps:

Maps
----

To follow this section you'll need to have Cartopy installed and working.

This script will plot the air temperature on a map.

.. ipython:: python
    :okwarning:

    import cartopy.crs as ccrs

    air = xr.tutorial.open_dataset("air_temperature").air

    p = air.isel(time=0).plot(
        subplot_kws=dict(projection=ccrs.Orthographic(-80, 35), facecolor="gray"),
        transform=ccrs.PlateCarree(),
    )
    p.axes.set_global()

    @savefig plotting_maps_cartopy.png width=100%
    p.axes.coastlines()

When faceting on maps, the projection can be transferred to the ``plot``
function using the ``subplot_kws`` keyword. The axes for the subplots created
by faceting are accessible in the object returned by ``plot``:

.. ipython:: python
    :okwarning:

    p = air.isel(time=[0, 4]).plot(
        transform=ccrs.PlateCarree(),
        col="time",
        subplot_kws={"projection": ccrs.Orthographic(-80, 35)},
    )
    for ax in p.axes.flat:
        ax.coastlines()
        ax.gridlines()
    @savefig plotting_maps_cartopy_facetting.png width=100%
    plt.draw()


Details
-------

Ways to Use
~~~~~~~~~~~

There are three ways to use the xarray plotting functionality:

1. Use ``plot`` as a convenience method for a DataArray.

2. Access a specific plotting method from the ``plot`` attribute of a
   DataArray.

3. Directly from the xarray plot submodule.

These are provided for user convenience; they all call the same code.

.. ipython:: python
    :okwarning:

    import xarray.plot as xplt

    da = xr.DataArray(range(5))
    fig, axes = plt.subplots(ncols=2, nrows=2)
    da.plot(ax=axes[0, 0])
    da.plot.line(ax=axes[0, 1])
    xplt.plot(da, ax=axes[1, 0])
    xplt.line(da, ax=axes[1, 1])
    plt.tight_layout()
    @savefig plotting_ways_to_use.png width=6in
    plt.draw()

Here the output is the same. Since the data is 1 dimensional the line plot
was used.

The convenience method :py:meth:`xarray.DataArray.plot` dispatches to an appropriate
plotting function based on the dimensions of the ``DataArray`` and whether
the coordinates are sorted and uniformly spaced. This table
describes what gets plotted:

=============== ===========================
Dimensions      Plotting function
--------------- ---------------------------
1               :py:func:`xarray.plot.line`
2               :py:func:`xarray.plot.pcolormesh`
Anything else   :py:func:`xarray.plot.hist`
=============== ===========================

Coordinates
~~~~~~~~~~~

If you'd like to find out what's really going on in the coordinate system,
read on.

.. ipython:: python

    a0 = xr.DataArray(np.zeros((4, 3, 2)), dims=("y", "x", "z"), name="temperature")
    a0[0, 0, 0] = 1
    a = a0.isel(z=0)
    a

The plot will produce an image corresponding to the values of the array.
Hence the top left pixel will be a different color than the others.
Before reading on, you may want to look at the coordinates and
think carefully about what the limits, labels, and orientation for
each of the axes should be.

.. ipython:: python
    :okwarning:

    @savefig plotting_example_2d_simple.png width=4in
    a.plot()

It may seem strange that
the values on the y axis are decreasing with -0.5 on the top. This is because
the pixels are centered over their coordinates, and the
axis labels and ranges correspond to the values of the
coordinates.

Multidimensional coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See also: :ref:`/examples/multidimensional-coords.ipynb`.

You can plot irregular grids defined by multidimensional coordinates with
xarray, but you'll have to tell the plot function to use these coordinates
instead of the default ones:

.. ipython:: python
    :okwarning:

    lon, lat = np.meshgrid(np.linspace(-20, 20, 5), np.linspace(0, 30, 4))
    lon += lat / 10
    lat += lon / 10
    da = xr.DataArray(
        np.arange(20).reshape(4, 5),
        dims=["y", "x"],
        coords={"lat": (("y", "x"), lat), "lon": (("y", "x"), lon)},
    )

    @savefig plotting_example_2d_irreg.png width=4in
    da.plot.pcolormesh("lon", "lat")

Note that in this case, xarray still follows the pixel centered convention.
This might be undesirable in some cases, for example when your data is defined
on a polar projection (:issue:`781`). This is why the default is to not follow
this convention when plotting on a map:

.. ipython:: python
    :okwarning:

    import cartopy.crs as ccrs

    ax = plt.subplot(projection=ccrs.PlateCarree())
    da.plot.pcolormesh("lon", "lat", ax=ax)
    ax.scatter(lon, lat, transform=ccrs.PlateCarree())
    ax.coastlines()
    @savefig plotting_example_2d_irreg_map.png width=4in
    ax.gridlines(draw_labels=True)

You can however decide to infer the cell boundaries and use the
``infer_intervals`` keyword:

.. ipython:: python
    :okwarning:

    ax = plt.subplot(projection=ccrs.PlateCarree())
    da.plot.pcolormesh("lon", "lat", ax=ax, infer_intervals=True)
    ax.scatter(lon, lat, transform=ccrs.PlateCarree())
    ax.coastlines()
    @savefig plotting_example_2d_irreg_map_infer.png width=4in
    ax.gridlines(draw_labels=True)

.. note::
    The data model of xarray does not support datasets with `cell boundaries`_
    yet. If you want to use these coordinates, you'll have to make the plots
    outside the xarray framework.

.. _cell boundaries: http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#cell-boundaries

One can also make line plots with multidimensional coordinates. In this case, ``hue`` must be a dimension name, not a coordinate name.

.. ipython:: python
    :okwarning:

    f, ax = plt.subplots(2, 1)
    da.plot.line(x="lon", hue="y", ax=ax[0])
    @savefig plotting_example_2d_hue_xy.png
    da.plot.line(x="lon", hue="x", ax=ax[1])
.. _time-series:

================
Time series data
================

A major use case for xarray is multi-dimensional time-series data.
Accordingly, we've copied many of features that make working with time-series
data in pandas such a joy to xarray. In most cases, we rely on pandas for the
core functionality.

.. ipython:: python
    :suppress:

    import numpy as np
    import pandas as pd
    import xarray as xr

    np.random.seed(123456)

Creating datetime64 data
------------------------

Xarray uses the numpy dtypes ``datetime64[ns]`` and ``timedelta64[ns]`` to
represent datetime data, which offer vectorized (if sometimes buggy) operations
with numpy and smooth integration with pandas.

To convert to or create regular arrays of ``datetime64`` data, we recommend
using :py:func:`pandas.to_datetime` and :py:func:`pandas.date_range`:

.. ipython:: python

    pd.to_datetime(["2000-01-01", "2000-02-02"])
    pd.date_range("2000-01-01", periods=365)

Alternatively, you can supply arrays of Python ``datetime`` objects. These get
converted automatically when used as arguments in xarray objects:

.. ipython:: python

    import datetime

    xr.Dataset({"time": datetime.datetime(2000, 1, 1)})

When reading or writing netCDF files, xarray automatically decodes datetime and
timedelta arrays using `CF conventions`_ (that is, by using a ``units``
attribute like ``'days since 2000-01-01'``).

.. _CF conventions: http://cfconventions.org

.. note::

   When decoding/encoding datetimes for non-standard calendars or for dates
   before year 1678 or after year 2262, xarray uses the `cftime`_ library.
   It was previously packaged with the ``netcdf4-python`` package under the
   name ``netcdftime`` but is now distributed separately. ``cftime`` is an
   :ref:`optional dependency<installing>` of xarray.

.. _cftime: https://unidata.github.io/cftime


You can manual decode arrays in this form by passing a dataset to
:py:func:`~xarray.decode_cf`:

.. ipython:: python

    attrs = {"units": "hours since 2000-01-01"}
    ds = xr.Dataset({"time": ("time", [0, 1, 2, 3], attrs)})
    xr.decode_cf(ds)

One unfortunate limitation of using ``datetime64[ns]`` is that it limits the
native representation of dates to those that fall between the years 1678 and
2262. When a netCDF file contains dates outside of these bounds, dates will be
returned as arrays of :py:class:`cftime.datetime` objects and a :py:class:`~xarray.CFTimeIndex`
will be used for indexing.  :py:class:`~xarray.CFTimeIndex` enables a subset of
the indexing functionality of a :py:class:`pandas.DatetimeIndex` and is only
fully compatible with the standalone version of ``cftime`` (not the version
packaged with earlier versions ``netCDF4``).  See :ref:`CFTimeIndex` for more
information.

Datetime indexing
-----------------

Xarray borrows powerful indexing machinery from pandas (see :ref:`indexing`).

This allows for several useful and succinct forms of indexing, particularly for
`datetime64` data. For example, we support indexing with strings for single
items and with the `slice` object:

.. ipython:: python

    time = pd.date_range("2000-01-01", freq="H", periods=365 * 24)
    ds = xr.Dataset({"foo": ("time", np.arange(365 * 24)), "time": time})
    ds.sel(time="2000-01")
    ds.sel(time=slice("2000-06-01", "2000-06-10"))

You can also select a particular time by indexing with a
:py:class:`datetime.time` object:

.. ipython:: python

    ds.sel(time=datetime.time(12))

For more details, read the pandas documentation and the section on `Indexing Using Datetime Components <datetime_component_indexing>`_ (i.e. using the ``.dt`` acessor).

.. _dt_accessor:

Datetime components
-------------------

Similar `to pandas`_, the components of datetime objects contained in a
given ``DataArray`` can be quickly computed using a special ``.dt`` accessor.

.. _to pandas: http://pandas.pydata.org/pandas-docs/stable/basics.html#basics-dt-accessors

.. ipython:: python

    time = pd.date_range("2000-01-01", freq="6H", periods=365 * 4)
    ds = xr.Dataset({"foo": ("time", np.arange(365 * 4)), "time": time})
    ds.time.dt.hour
    ds.time.dt.dayofweek

The ``.dt`` accessor works on both coordinate dimensions as well as
multi-dimensional data.

Xarray also supports a notion of "virtual" or "derived" coordinates for
`datetime components`__ implemented by pandas, including "year", "month",
"day", "hour", "minute", "second", "dayofyear", "week", "dayofweek", "weekday"
and "quarter":

__ http://pandas.pydata.org/pandas-docs/stable/api.html#time-date-components

.. ipython:: python

    ds["time.month"]
    ds["time.dayofyear"]

For use as a derived coordinate, xarray adds ``'season'`` to the list of
datetime components supported by pandas:

.. ipython:: python

    ds["time.season"]
    ds["time"].dt.season

The set of valid seasons consists of 'DJF', 'MAM', 'JJA' and 'SON', labeled by
the first letters of the corresponding months.

You can use these shortcuts with both Datasets and DataArray coordinates.

In addition, xarray supports rounding operations ``floor``, ``ceil``, and ``round``. These operations require that you supply a `rounding frequency as a string argument.`__

__ http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases

.. ipython:: python

    ds["time"].dt.floor("D")

The ``.dt`` accessor can also be used to generate formatted datetime strings
for arrays utilising the same formatting as the standard `datetime.strftime`_.

.. _datetime.strftime: https://docs.python.org/3/library/datetime.html#strftime-strptime-behavior

.. ipython:: python

    ds["time"].dt.strftime("%a, %b %d %H:%M")

.. _datetime_component_indexing:

Indexing Using Datetime Components
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You can use use the ``.dt`` accessor when subsetting your data as well. For example, we can subset for the month of January using the following:

.. ipython:: python

    ds.isel(time=(ds.time.dt.month == 1))

You can also search for multiple months (in this case January through March), using ``isin``:

.. ipython:: python

    ds.isel(time=ds.time.dt.month.isin([1, 2, 3]))

.. _resampling:

Resampling and grouped operations
---------------------------------

Datetime components couple particularly well with grouped operations (see
:ref:`groupby`) for analyzing features that repeat over time. Here's how to
calculate the mean by time of day:

.. ipython:: python
    :okwarning:

    ds.groupby("time.hour").mean()

For upsampling or downsampling temporal resolutions, xarray offers a
:py:meth:`~xarray.Dataset.resample` method building on the core functionality
offered by the pandas method of the same name. Resample uses essentially the
same api as ``resample`` `in pandas`_.

.. _in pandas: http://pandas.pydata.org/pandas-docs/stable/timeseries.html#up-and-downsampling

For example, we can downsample our dataset from hourly to 6-hourly:

.. ipython:: python
    :okwarning:

    ds.resample(time="6H")

This will create a specialized ``Resample`` object which saves information
necessary for resampling. All of the reduction methods which work with
``Resample`` objects can also be used for resampling:

.. ipython:: python
    :okwarning:

    ds.resample(time="6H").mean()

You can also supply an arbitrary reduction function to aggregate over each
resampling group:

.. ipython:: python

    ds.resample(time="6H").reduce(np.mean)

For upsampling, xarray provides six methods: ``asfreq``, ``ffill``, ``bfill``, ``pad``,
``nearest`` and ``interpolate``. ``interpolate`` extends ``scipy.interpolate.interp1d``
and supports all of its schemes. All of these resampling operations work on both
Dataset and DataArray objects with an arbitrary number of dimensions.

In order to limit the scope of the methods ``ffill``, ``bfill``, ``pad`` and
``nearest`` the ``tolerance`` argument can be set in coordinate units.
Data that has indices outside of the given ``tolerance`` are set to ``NaN``.

.. ipython:: python

    ds.resample(time="1H").nearest(tolerance="1H")


For more examples of using grouped operations on a time dimension, see
:doc:`../examples/weather-data`.

.. _internals.duck_arrays:

Integrating with duck arrays
=============================

.. warning::

    This is a experimental feature.

Xarray can wrap custom :term:`duck array` objects as long as they define numpy's
``shape``, ``dtype`` and ``ndim`` properties and the ``__array__``,
``__array_ufunc__`` and ``__array_function__`` methods.

In certain situations (e.g. when printing the collapsed preview of
variables of a ``Dataset``), xarray will display the repr of a :term:`duck array`
in a single line, truncating it to a certain number of characters. If that
would drop too much information, the :term:`duck array` may define a
``_repr_inline_`` method that takes ``max_width`` (number of characters) as an
argument:

.. code:: python

    class MyDuckArray:
        ...

        def _repr_inline_(self, max_width):
            """format to a single line with at most max_width characters"""
            ...

        ...

To avoid duplicated information, this method must omit information about the shape and
:term:`dtype`. For example, the string representation of a ``dask`` array or a
``sparse`` matrix would be:

.. ipython:: python

    import dask.array as da
    import xarray as xr
    import sparse

    a = da.linspace(0, 1, 20, chunks=2)
    a

    b = np.eye(10)
    b[[5, 7, 3, 0], [6, 8, 2, 9]] = 2
    b = sparse.COO.from_numpy(b)
    b

    xr.Dataset(dict(a=("x", a), b=(("y", "z"), b)))

Extending xarray
================

.. ipython:: python
    :suppress:

    import xarray as xr


Xarray is designed as a general purpose library, and hence tries to avoid
including overly domain specific functionality. But inevitably, the need for more
domain specific logic arises.

One standard solution to this problem is to subclass Dataset and/or DataArray to
add domain specific functionality. However, inheritance is not very robust. It's
easy to inadvertently use internal APIs when subclassing, which means that your
code may break when xarray upgrades. Furthermore, many builtin methods will
only return native xarray objects.

The standard advice is to use `composition over inheritance`__, but
reimplementing an API as large as xarray's on your own objects can be an onerous
task, even if most methods are only forwarding to xarray implementations.

__ https://github.com/pydata/xarray/issues/706

If you simply want the ability to call a function with the syntax of a
method call, then the builtin :py:meth:`~xarray.DataArray.pipe` method (copied
from pandas) may suffice.

To resolve this issue for more complex cases, xarray has the
:py:func:`~xarray.register_dataset_accessor` and
:py:func:`~xarray.register_dataarray_accessor` decorators for adding custom
"accessors" on xarray objects. Here's how you might use these decorators to
write a custom "geo" accessor implementing a geography specific extension to
xarray:

.. literalinclude:: ../examples/_code/accessor_example.py

In general, the only restriction on the accessor class is that the ``__init__`` method
must have a single parameter: the ``Dataset`` or ``DataArray`` object it is supposed
to work on.

This achieves the same result as if the ``Dataset`` class had a cached property
defined that returns an instance of your class:

.. code-block:: python

    class Dataset:
        ...

        @property
        def geo(self):
            return GeoAccessor(self)

However, using the register accessor decorators is preferable to simply adding
your own ad-hoc property (i.e., ``Dataset.geo = property(...)``), for several
reasons:

1. It ensures that the name of your property does not accidentally conflict with
   any other attributes or methods (including other accessors).
2. Instances of accessor object will be cached on the xarray object that creates
   them. This means you can save state on them (e.g., to cache computed
   properties).
3. Using an accessor provides an implicit namespace for your custom
   functionality that clearly identifies it as separate from built-in xarray
   methods.

.. note::

   Accessors are created once per DataArray and Dataset instance. New
   instances, like those created from arithmetic operations or when accessing
   a DataArray from a Dataset (ex. ``ds[var_name]``), will have new
   accessors created.

Back in an interactive IPython session, we can use these properties:

.. ipython:: python
    :suppress:

    exec(open("examples/_code/accessor_example.py").read())

.. ipython:: python

    ds = xr.Dataset({"longitude": np.linspace(0, 10), "latitude": np.linspace(0, 20)})
    ds.geo.center
    ds.geo.plot()

The intent here is that libraries that extend xarray could add such an accessor
to implement subclass specific functionality rather than using actual subclasses
or patching in a large number of domain specific methods. For further reading
on ways to write new accessors and the philosophy behind the approach, see
:issue:`1080`.

To help users keep things straight, please `let us know
<https://github.com/pydata/xarray/issues>`_ if you plan to write a new accessor
for an open source library. In the future, we will maintain a list of accessors
and the libraries that implement them on this page.

To make documenting accessors with ``sphinx`` and ``sphinx.ext.autosummary``
easier, you can use `sphinx-autosummary-accessors`_.

.. _sphinx-autosummary-accessors: https://sphinx-autosummary-accessors.readthedocs.io/
.. _add_a_backend:

How to add a new backend
------------------------

Adding a new backend for read support to Xarray does not require
to integrate any code in Xarray; all you need to do is:

- Create a class that inherits from Xarray :py:class:`~xarray.backends.BackendEntrypoint`
  and implements the method ``open_dataset`` see :ref:`RST backend_entrypoint`

- Declare this class as an external plugin in your ``setup.py``, see :ref:`RST backend_registration`

If you also want to support lazy loading and dask see :ref:`RST lazy_loading`.

Note that the new interface for backends is available from Xarray
version >= 0.18 onwards.

.. _RST backend_entrypoint:

BackendEntrypoint subclassing
+++++++++++++++++++++++++++++

Your ``BackendEntrypoint`` sub-class is the primary interface with Xarray, and
it should implement the following attributes and methods:

- the ``open_dataset`` method (mandatory)
- the ``open_dataset_parameters`` attribute (optional)
- the ``guess_can_open`` method (optional).

This is what a ``BackendEntrypoint`` subclass should look like:

.. code-block:: python

    from xarray.backends import BackendEntrypoint


    class MyBackendEntrypoint(BackendEntrypoint):
        def open_dataset(
            self,
            filename_or_obj,
            *,
            drop_variables=None,
            # other backend specific keyword arguments
            # `chunks` and `cache` DO NOT go here, they are handled by xarray
        ):
            return my_open_dataset(filename_or_obj, drop_variables=drop_variables)

        open_dataset_parameters = ["filename_or_obj", "drop_variables"]

        def guess_can_open(self, filename_or_obj):
            try:
                _, ext = os.path.splitext(filename_or_obj)
            except TypeError:
                return False
            return ext in {".my_format", ".my_fmt"}

``BackendEntrypoint`` subclass methods and attributes are detailed in the following.

.. _RST open_dataset:

open_dataset
^^^^^^^^^^^^

The backend ``open_dataset`` shall implement reading from file, the variables
decoding and it shall instantiate the output Xarray class :py:class:`~xarray.Dataset`.

The following is an example of the high level processing steps:

.. code-block:: python

    def open_dataset(
        self,
        filename_or_obj,
        *,
        drop_variables=None,
        decode_times=True,
        decode_timedelta=True,
        decode_coords=True,
        my_backend_option=None,
    ):
        vars, attrs, coords = my_reader(
            filename_or_obj,
            drop_variables=drop_variables,
            my_backend_option=my_backend_option,
        )
        vars, attrs, coords = my_decode_variables(
            vars, attrs, decode_times, decode_timedelta, decode_coords
        )  #  see also conventions.decode_cf_variables

        ds = xr.Dataset(vars, attrs=attrs, coords=coords)
        ds.set_close(my_close_method)

        return ds


The output :py:class:`~xarray.Dataset` shall implement the additional custom method
``close``, used by Xarray to ensure the related files are eventually closed. This
method shall be set by using :py:meth:`~xarray.Dataset.set_close`.


The input of ``open_dataset`` method are one argument
(``filename_or_obj``) and one keyword argument (``drop_variables``):

- ``filename_or_obj``: can be any object but usually it is a string containing a path or an instance of
  :py:class:`pathlib.Path`.
- ``drop_variables``: can be `None` or an iterable containing the variable
  names to be dropped when reading the data.

If it makes sense for your backend, your ``open_dataset``  method
should implement in its interface the following boolean keyword arguments, called
**decoders**, which default to ``None``:

- ``mask_and_scale``
- ``decode_times``
- ``decode_timedelta``
- ``use_cftime``
- ``concat_characters``
- ``decode_coords``

Note: all the supported decoders shall be declared explicitly
in backend ``open_dataset`` signature and adding a ``**kargs`` is not allowed.

These keyword arguments are explicitly defined in Xarray
:py:func:`~xarray.open_dataset` signature. Xarray will pass them to the
backend only if the User explicitly sets a value different from ``None``.
For more details on decoders see :ref:`RST decoders`.

Your backend can also take as input a set of backend-specific keyword
arguments. All these keyword arguments can be passed to
:py:func:`~xarray.open_dataset` grouped either via the ``backend_kwargs``
parameter or explicitly using the syntax ``**kwargs``.


If you don't want to support the lazy loading, then the
:py:class:`~xarray.Dataset` shall contain values as a :py:class:`numpy.ndarray`
and your work is almost done.

.. _RST open_dataset_parameters:

open_dataset_parameters
^^^^^^^^^^^^^^^^^^^^^^^

``open_dataset_parameters`` is the list of backend ``open_dataset`` parameters.
It is not a mandatory parameter, and if the backend does not provide it
explicitly, Xarray creates a list of them automatically by inspecting the
backend signature.

If ``open_dataset_parameters`` is not defined, but ``**kwargs`` and ``*args``
are in the backend ``open_dataset`` signature, Xarray raises an error.
On the other hand, if the backend provides the ``open_dataset_parameters``,
then ``**kwargs`` and ``*args`` can be used in the signature.
However, this practice is discouraged unless there is a good reasons for using
``**kwargs`` or ``*args``.

.. _RST guess_can_open:

guess_can_open
^^^^^^^^^^^^^^

``guess_can_open`` is used to identify the proper engine to open your data
file automatically in case the engine is not specified explicitly. If you are
not interested in supporting this feature, you can skip this step since
:py:class:`~xarray.backends.BackendEntrypoint` already provides a
default :py:meth:`~xarray.backends.BackendEntrypoint.guess_can_open`
that always returns ``False``.

Backend ``guess_can_open`` takes as input the ``filename_or_obj`` parameter of
Xarray :py:meth:`~xarray.open_dataset`, and returns a boolean.

.. _RST decoders:

Decoders
^^^^^^^^
The decoders implement specific operations to transform data from on-disk
representation to Xarray representation.

A classic example is the “time” variable decoding operation. In NetCDF, the
elements of the “time” variable are stored as integers, and the unit contains
an origin (for example: "seconds since 1970-1-1"). In this case, Xarray
transforms the pair integer-unit in a :py:class:`numpy.datetime64`.

The standard coders implemented in Xarray are:

- :py:class:`xarray.coding.strings.CharacterArrayCoder()`
- :py:class:`xarray.coding.strings.EncodedStringCoder()`
- :py:class:`xarray.coding.variables.UnsignedIntegerCoder()`
- :py:class:`xarray.coding.variables.CFMaskCoder()`
- :py:class:`xarray.coding.variables.CFScaleOffsetCoder()`
- :py:class:`xarray.coding.times.CFTimedeltaCoder()`
- :py:class:`xarray.coding.times.CFDatetimeCoder()`

Xarray coders all have the same interface. They have two methods: ``decode``
and ``encode``. The method ``decode`` takes a ``Variable`` in on-disk
format and returns a ``Variable`` in Xarray format. Variable
attributes no more applicable after the decoding, are dropped and stored in the
``Variable.encoding`` to make them available to the ``encode`` method, which
performs the inverse transformation.

In the following an example on how to use the coders ``decode`` method:

.. ipython:: python

    var = xr.Variable(
        dims=("x",), data=np.arange(10.0), attrs={"scale_factor": 10, "add_offset": 2}
    )
    var

    coder = xr.coding.variables.CFScaleOffsetCoder()
    decoded_var = coder.decode(var)
    decoded_var
    decoded_var.encoding

Some of the transformations can be common to more backends, so before
implementing a new decoder, be sure Xarray does not already implement that one.

The backends can reuse Xarray’s decoders, either instantiating the coders
and using the method ``decode`` directly or using the higher-level function
:py:func:`~xarray.conventions.decode_cf_variables` that groups Xarray decoders.

In some cases, the transformation to apply strongly depends on the on-disk
data format. Therefore, you may need to implement your own decoder.

An example of such a case is when you have to deal with the time format of a
grib file. grib format is very different from the NetCDF one: in grib, the
time is stored in two attributes dataDate and dataTime as strings. Therefore,
it is not possible to reuse the Xarray time decoder, and implementing a new
one is mandatory.

Decoders can be activated or deactivated using the boolean keywords of
Xarray :py:meth:`~xarray.open_dataset` signature: ``mask_and_scale``,
``decode_times``, ``decode_timedelta``, ``use_cftime``,
``concat_characters``, ``decode_coords``.
Such keywords are passed to the backend only if the User sets a value
different from ``None``.  Note that the backend does not necessarily have to
implement all the decoders, but it shall declare in its ``open_dataset``
interface only the boolean keywords related to the supported decoders.

.. _RST backend_registration:

How to register a backend
+++++++++++++++++++++++++++

Define a new entrypoint in your ``setup.py`` (or ``setup.cfg``) with:

- group: ``xarray.backends``
- name: the name to be passed to :py:meth:`~xarray.open_dataset`  as ``engine``
- object reference: the reference of the class that you have implemented.

You can declare the entrypoint in ``setup.py`` using the following syntax:

.. code-block::

    setuptools.setup(
        entry_points={
            "xarray.backends": ["my_engine=my_package.my_module:MyBackendEntryClass"],
        },
    )

in ``setup.cfg``:

.. code-block:: cfg

    [options.entry_points]
    xarray.backends =
        my_engine = my_package.my_module:MyBackendEntryClass


See https://packaging.python.org/specifications/entry-points/#data-model
for more information

If you are using `Poetry <https://python-poetry.org/>`_ for your build system, you can accomplish the same thing using "plugins". In this case you would need to add the following to your ``pyproject.toml`` file:

.. code-block:: toml

    [tool.poetry.plugins."xarray_backends"]
    "my_engine" = "my_package.my_module:MyBackendEntryClass"

See https://python-poetry.org/docs/pyproject/#plugins for more information on Poetry plugins.

.. _RST lazy_loading:

How to support Lazy Loading
+++++++++++++++++++++++++++
If you want to make your backend effective with big datasets, then you should
support lazy loading.
Basically, you shall replace the :py:class:`numpy.ndarray` inside the
variables with a custom class that supports lazy loading indexing.
See the example below:

.. code-block:: python

    backend_array = MyBackendArray()
    data = indexing.LazilyIndexedArray(backend_array)
    var = xr.Variable(dims, data, attrs=attrs, encoding=encoding)

Where:

- :py:class:`~xarray.core.indexing.LazilyIndexedArray` is a class
  provided by Xarray that manages the lazy loading.
- ``MyBackendArray`` shall be implemented by the backend and shall inherit
  from :py:class:`~xarray.backends.BackendArray`.

BackendArray subclassing
^^^^^^^^^^^^^^^^^^^^^^^^

The BackendArray subclass shall implement the following method and attributes:

- the ``__getitem__`` method that takes in input an index and returns a
  `NumPy <https://numpy.org/>`__ array
- the ``shape`` attribute
- the ``dtype`` attribute.


Xarray supports different type of
`indexing <http://xarray.pydata.org/en/stable/indexing.html>`__, that can be
grouped in three types of indexes
:py:class:`~xarray.core.indexing.BasicIndexer`,
:py:class:`~xarray.core.indexing.OuterIndexer` and
:py:class:`~xarray.core.indexing.VectorizedIndexer`.
This implies that the implementation of the method ``__getitem__`` can be tricky.
In oder to simplify this task, Xarray provides a helper function,
:py:func:`~xarray.core.indexing.explicit_indexing_adapter`, that transforms
all the input  ``indexer`` types (`basic`, `outer`, `vectorized`) in a tuple
which is interpreted correctly by your backend.

This is an example ``BackendArray`` subclass implementation:

.. code-block:: python

    from xarray.backends import BackendArray


    class MyBackendArray(BackendArray):
        def __init__(
            self,
            shape,
            dtype,
            lock,
            # other backend specific keyword arguments
        ):
            self.shape = shape
            self.dtype = lock
            self.lock = dtype

        def __getitem__(
            self, key: xarray.core.indexing.ExplicitIndexer
        ) -> np.typing.ArrayLike:
            return indexing.explicit_indexing_adapter(
                key,
                self.shape,
                indexing.IndexingSupport.BASIC,
                self._raw_indexing_method,
            )

        def _raw_indexing_method(self, key: tuple) -> np.typing.ArrayLike:
            # thread safe method that access to data on disk
            with self.lock:
                ...
                return item

Note that ``BackendArray.__getitem__`` must be thread safe to support
multi-thread processing.

The :py:func:`~xarray.core.indexing.explicit_indexing_adapter` method takes in
input the ``key``, the array ``shape`` and the following parameters:

- ``indexing_support``: the type of index supported by ``raw_indexing_method``
- ``raw_indexing_method``: a method that shall take in input a key in the form
  of a tuple and return an indexed :py:class:`numpy.ndarray`.

For more details see
:py:class:`~xarray.core.indexing.IndexingSupport` and :ref:`RST indexing`.

In order to support `Dask <http://dask.pydata.org/>`__ distributed and
:py:mod:`multiprocessing`, ``BackendArray`` subclass should be serializable
either with :ref:`io.pickle` or
`cloudpickle <https://github.com/cloudpipe/cloudpickle>`__.
That implies that all the reference to open files should be dropped. For
opening files, we therefore suggest to use the helper class provided by Xarray
:py:class:`~xarray.backends.CachingFileManager`.

.. _RST indexing:

Indexing Examples
^^^^^^^^^^^^^^^^^
**BASIC**

In the ``BASIC`` indexing support, numbers and slices are supported.

Example:

.. ipython::
    :verbatim:

    In [1]: # () shall return the full array
       ...: backend_array._raw_indexing_method(())
    Out[1]: array([[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]])

    In [2]: # shall support integers
       ...: backend_array._raw_indexing_method(1, 1)
    Out[2]: 5

    In [3]: # shall support slices
       ...: backend_array._raw_indexing_method(slice(0, 3), slice(2, 4))
    Out[3]: array([[2, 3], [6, 7], [10, 11]])

**OUTER**

The ``OUTER`` indexing shall support number, slices and in addition it shall
support also lists of integers. The the outer indexing is equivalent to
combining multiple input list with ``itertools.product()``:

.. ipython::
    :verbatim:

    In [1]: backend_array._raw_indexing_method([0, 1], [0, 1, 2])
    Out[1]: array([[0, 1, 2], [4, 5, 6]])

    # shall support integers
    In [2]: backend_array._raw_indexing_method(1, 1)
    Out[2]: 5


**OUTER_1VECTOR**

The ``OUTER_1VECTOR`` indexing shall supports number, slices and at most one
list. The behaviour with the list shall be the same of ``OUTER`` indexing.

If you support more complex indexing as `explicit indexing` or
`numpy indexing`, you can have a look to the implemetation of Zarr backend and Scipy backend,
currently available in :py:mod:`~xarray.backends` module.

.. _RST preferred_chunks:

Backend preferred chunks
^^^^^^^^^^^^^^^^^^^^^^^^

The backend is not directly involved in `Dask <http://dask.pydata.org/>`__
chunking, since it is internally managed by Xarray. However, the backend can
define the preferred chunk size inside the variable’s encoding
``var.encoding["preferred_chunks"]``. The ``preferred_chunks`` may be useful
to improve performances with lazy loading. ``preferred_chunks`` shall be a
dictionary specifying chunk size per dimension like
``{“dim1”: 1000, “dim2”: 2000}``  or
``{“dim1”: [1000, 100], “dim2”: [2000, 2000, 2000]]}``.

The ``preferred_chunks`` is used by Xarray to define the chunk size in some
special cases:

- if ``chunks`` along a dimension is ``None`` or not defined
- if ``chunks`` is ``"auto"``.

In the first case Xarray uses the chunks size specified in
``preferred_chunks``.
In the second case Xarray accommodates ideal chunk sizes, preserving if
possible the "preferred_chunks". The ideal chunk size is computed using
:py:func:`dask.array.core.normalize_chunks`, setting
``previous_chunks = preferred_chunks``.
.. currentmodule:: xarray

.. _zarr_encoding:

Zarr Encoding Specification
============================

In implementing support for the `Zarr <https://zarr.readthedocs.io/>`_ storage
format, Xarray developers made some *ad hoc* choices about how to store
NetCDF data in Zarr.
Future versions of the Zarr spec will likely include a more formal convention
for the storage of the NetCDF data model in Zarr; see
`Zarr spec repo <https://github.com/zarr-developers/zarr-specs>`_ for ongoing
discussion.

First, Xarray can only read and write Zarr groups. There is currently no support
for reading / writting individual Zarr arrays. Zarr groups are mapped to
Xarray ``Dataset`` objects.

Second, from Xarray's point of view, the key difference between
NetCDF and Zarr is that all NetCDF arrays have *dimension names* while Zarr
arrays do not. Therefore, in order to store NetCDF data in Zarr, Xarray must
somehow encode and decode the name of each array's dimensions.

To accomplish this, Xarray developers decided to define a special Zarr array
attribute: ``_ARRAY_DIMENSIONS``. The value of this attribute is a list of
dimension names (strings), for example ``["time", "lon", "lat"]``. When writing
data to Zarr, Xarray sets this attribute on all variables based on the variable
dimensions. When reading a Zarr group, Xarray looks for this attribute on all
arrays, raising an error if it can't be found. The attribute is used to define
the variable dimension names and then removed from the attributes dictionary
returned to the user.

Because of these choices, Xarray cannot read arbitrary array data, but only
Zarr data with valid ``_ARRAY_DIMENSIONS`` attributes on each array.

After decoding the ``_ARRAY_DIMENSIONS`` attribute and assigning the variable
dimensions, Xarray proceeds to [optionally] decode each variable using its
standard CF decoding machinery used for NetCDF data (see :py:func:`decode_cf`).

Finally, it's worth noting that Xarray writes (and attempts to read)
"consolidated metadata" by default (the ``.zmetadata`` file), which is another
non-standard Zarr extension, albeit one implemented upstream in Zarr-Python.
You do not need to write consolidated metadata to make Zarr stores readable in
Xarray, but because Xarray can open these stores much faster, users will see a
warning about poor performance when reading non-consolidated stores unless they
explicitly set ``consolidated=False``. See :ref:`io.zarr.consolidated_metadata`
for more details.

As a concrete example, here we write a tutorial dataset to Zarr and then
re-open it directly with Zarr:

.. ipython:: python

    import os
    import xarray as xr
    import zarr

    ds = xr.tutorial.load_dataset("rasm")
    ds.to_zarr("rasm.zarr", mode="w")

    zgroup = zarr.open("rasm.zarr")
    print(os.listdir("rasm.zarr"))
    print(zgroup.tree())
    dict(zgroup["Tair"].attrs)
Variable objects
================

The core internal data structure in xarray is the :py:class:`~xarray.Variable`,
which is used as the basic building block behind xarray's
:py:class:`~xarray.Dataset` and :py:class:`~xarray.DataArray` types. A
``Variable`` consists of:

- ``dims``: A tuple of dimension names.
- ``data``: The N-dimensional array (typically, a NumPy or Dask array) storing
  the Variable's data. It must have the same number of dimensions as the length
  of ``dims``.
- ``attrs``: An ordered dictionary of metadata associated with this array. By
  convention, xarray's built-in operations never use this metadata.
- ``encoding``: Another ordered dictionary used to store information about how
  these variable's data is represented on disk. See :ref:`io.encoding` for more
  details.

``Variable`` has an interface similar to NumPy arrays, but extended to make use
of named dimensions. For example, it uses ``dim`` in preference to an ``axis``
argument for methods like ``mean``, and supports :ref:`compute.broadcasting`.

However, unlike ``Dataset`` and ``DataArray``, the basic ``Variable`` does not
include coordinate labels along each axis.

``Variable`` is public API, but because of its incomplete support for labeled
data, it is mostly intended for advanced uses, such as in xarray itself or for
writing new backends. You can access the variable objects that correspond to
xarray objects via the (readonly) :py:attr:`Dataset.variables
<xarray.Dataset.variables>` and
:py:attr:`DataArray.variable <xarray.DataArray.variable>` attributes.
.. _internals:

xarray Internals
================

Xarray builds upon two of the foundational libraries of the scientific Python
stack, NumPy and pandas. It is written in pure Python (no C or Cython
extensions), which makes it easy to develop and extend. Instead, we push
compiled code to :ref:`optional dependencies<installing>`.


.. toctree::
   :maxdepth: 2
   :hidden:

   variable-objects
   duck-arrays-integration
   extending-xarray
   zarr-encoding-spec
   how-to-add-new-backend
