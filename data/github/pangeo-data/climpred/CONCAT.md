Our Code of Conduct can be found on our documentation [here](https://climpred.readthedocs.io/en/stable/code_of_conduct.html).
A guide on how to contribute to `climpred` can be found on our documentation [here](https://climpred.readthedocs.io/en/stable/contributing.html).
# Description

Please include a summary of the change and which issue is fixed. Please also include relevant motivation and context. List any dependencies that are required for this change.

Closes #(issue)

## To-Do List

<!-- Feel free to add a checklist of steps to be performed while you are working through creating this PR. -->

## Type of change

<!-- Please delete options that are not relevant.-->

-   [ ]  Bug fix (non-breaking change which fixes an issue)
-   [ ]  New feature (non-breaking change which adds functionality)
-   [ ]  Breaking change (fix or feature that would cause existing functionality to not work as expected)
-   [ ]  This change requires a documentation update
-   [ ]  Performance (if you modified existing code run `asv` to detect performance changes)
-   [ ]  Refactoring
-   [ ]  Improved Documentation

# How Has This Been Tested?

<!-- Please describe the tests that you ran to verify your changes. This could point to a cell in the updated notebooks. Or a snippet of code with accompanying figures here. -->

-   [ ]  Tests added for `pytest`, if necessary.

## Checklist (while developing)

-   [ ]  I have added docstrings to all new functions.
-   [ ]  I have commented my code, particularly in hard-to-understand areas.
-   [ ]  I have updated the sphinx documentation, if necessary.
-   [ ]  Any new functions are added to the API. (See contribution guide)
-   [ ]  CHANGELOG is updated with reference to this PR.

## References

<!-- Please add any references to manuscripts, textbooks, etc. if applicable -->
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**Code Sample**

A "Minimal, Complete and Verifiable Example" will make it much easier for maintainers to help you:
<http://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports>. I.e., this shouldn't rely upon a file on your local machine. Can you recreate this bug with a small dummy dataset of randomly generated `numpy`/`xarray` data?

```python
# Your code here
```

**Expected behavior**
A clear and concise description of what you expected to happen.

**Output of `climpred.show_versions()`**

<details>
# Paste the output here climpred.show_versions() here
</details>

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: feature request
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
title: 'climpred: Verification of weather and climate forecasts'
tags:
  - python
  - climate
  - forecasting
  - geospatial
  - big data
authors:
  - name: Riley X. Brady
    orcid: 0000-0002-2309-8245
    affiliation: 1
  - name: Aaron Spring
    orcid: 0000-0003-0216-2241
    affiliation: "2, 3"
affiliations:
 - name: Department of Atmospheric and Oceanic Sciences and Institute of Arctic and Alpine Research, University of Colorado Boulder, Boulder, Colorado USA
   index: 1
 - name: Max Planck Institute for Meteorology, Hamburg, Germany
   index: 2
 - name: International Max Planck Research School on Earth System Modelling, Hamburg, Germany
   index: 3
date: 22 February 2021
bibliography: paper.bib
---

<!-- 991 words total -->

# Summary
<!-- 71 words -->
Predicting extreme events and variations in weather and climate provides crucial
information for economic, social, and environmental decision-making [@Merryfield:2020].
However, quantifying prediction skill for multi-dimensional geospatial model output is
computationally expensive and a difficult coding challenge. The large datasets
(order gigabytes to terabytes) require parallel and out-of-memory computing to be
analyzed efficiently. Further, aligning the many forecast initializations with differing
observational products is a straight-forward, but exhausting and error-prone exercise
for researchers.

<!-- 68 words -->
To simplify and standardize forecast verification across scales from hourly weather to
decadal climate forecasts, we built `climpred`: a community-driven python package for
computationally efficient and methodologically consistent verification of ensemble
prediction models. The code base is maintained through open-source development. It
leverages `xarray` [@Hoyer:2017] to anticipate core prediction ensemble dimensions
(ensemble `member`, `init`ialization date and `lead` time) and `dask`
[@dask; @Rocklin:2015] to perform out-of-memory and parallelized computations on
large datasets.

 <!-- 54 words -->
`climpred` aims to offer a comprehensive set of analysis tools for assessing the quality
of dynamical forecasts relative to verification products (_e.g._, observations,
reanalysis products, control simulations). The package includes a suite of deterministic
and probabilistic verification metrics that are constantly expanded by the community and
are generally organized in our companion package, `xskillscore`.

<!-- 173 words -->
# Statement of Need
While other climate verification packages exist (_e.g._, `s2dverification`
[@Manubens:2018] written in R and `MurCSS` [@Illing:2014] written with python-based
`CDO`-bindings [@CDO]), `climpred` is unique for many reasons.

1. `climpred` spans broad temporal scales of prediction, supporting the weather,
subseasonal-to-seasonal (S2S), and seasonal-to-decadal (S2D) communities.

2. `climpred` is highly modular and supports the research process from end-to-end,
from loading in model output, to interactive pre-processing and analysis, to
visualization.

3. `climpred` supports `dask` [@dask; @Rocklin:2015] and thus works across all
computational scales, from personal laptops to supercomputers (HPC).

4. Flexibility and scaling leads to verification of global 5° x 5° resolution climate
predictions in 8 seconds, compared to the 8 minutes required by `MurCSS`. However,
note that `climpred` modularizes its workflow such that the verification step is
performed on already pre-processed output, while `MurCSS` uses a more rigid framework
that always required pre-processing. This time scale of seconds allows for a truly
interactive analysis experience.

5. `climpred` is part of the wider scientific python community, `pangeo`
[@Abernathey:2017; @Eynard:2019]. A wide adoption of `climpred` could standardize prediction model
evaluation and make verification reproducible [@Irving:2015].

<!-- 207 words -->
# Prediction Simulation Types
Weather and climate modeling institutions typically run so-called “hindcasts," where
dynamical models are retrospectively initialized from many past observed climate states
[@Meehl:2009]. Initializations are then slightly perturbed to generate an ensemble of
forecasts that diverge solely due to their sensitive dependence on initial conditions
[@Lorenz:1963]. Hindcasts are evaluated by using some statistical metric to score their
performance against historical observations. “Skill" is established by comparing these
results to the performance of some “reference" forecast (@Jolliffe:2012;
_e.g._, a persistence forecast). The main assumption is that the skill established
relative to the past will propagate to forecasts of the future.

A more idealized approach is the so-called “perfect-model" framework, which is ideal for
investigating processes leading to potentially exploitable predictability
[@Griffies:1997; @Bushuk:2018; @Seferian:2018; @Spring:2020]. Ensemble members are spun
off an individual model (by slightly perturbing its state) to predict its own evolution.
This avoids initialization shocks [@Kroger:2017], since the framework is self-contained.
However, it cannot predict the real world. The perfect-model setup rather estimates the
theoretical upper limit timescale after which the value of dynamical initialization is
lost due to chaos in the Earth system, assuming that the model perfectly replicates the
dynamics of the real world. Skill quantification is accomplished by considering one
ensemble member as the verification data and the remaining members as the forecasts
[@Griffies:1997].

<!-- 360 words -->
# Climpred Classes and Object-Oriented Verification
`climpred` supports both prediction system formats, offering `HindcastEnsemble` and
`PerfectModelEnsemble` objects. `HindcastEnsemble` is instantiated with an `initialized`
hindcast ensemble dataset and requires an `observation`al dataset against which to
verify. `PerfectModelEnsemble` is instantiated with an `initialized` perfect-model
ensemble dataset and also accepts a `control` dataset against which to evaluate
forecasts. Both objects can also track an `uninitialized` dataset, which represents a
historical simulation that evolves solely due to random internal climate variability or
can be used to isolate the influence of external forcing [e.g., @Kay:2014].

Assessing skill for `PredictionEnsemble` objects (the parent class to `HindcastEnsemble`
and `PerfectModelEnsemble`) is standardized into a one-liner:

```python
PredictionEnsemble.verify(
    # Score forecast using the Anomaly Correlation Coefficient.
    metric='acc',
    # Compare the ensemble mean to observations.
    comparison='e2o',
    # Keep the same set of initializations at each lead time.
    alignment='same_inits',
    # Reduce the verification over the initialization dimension.
    dim='init',
    # Score performance of a persistence forecast as well.
    reference='persistence',
)
```

Each keyword argument allows flexibility from the user’s end—one can select from a
library of metrics, comparison types, alignment strategies, dimensional reductions, and
reference forecasts. The most unique feature to `climpred`, however, is the ability for
users to choose the alignment strategy to pair initialization dates with verification
dates over numerous lead times. In other words, initialization dates need to be
converted to target forecast dates by shifting them using the lead time coordinate. This
is tedious, since one must remedy disparities in calendar types between the model and
observations and account for the time span of or gaps in observations relative to the
time span of the model.

There is seemingly no unified approach to how hindcast initialization dates are aligned
with observational dates in the academic literature. The authors of `climpred` thus
identified three techniques, which can be selected by the user:

1. Maximize the degrees of freedom by selecting all initialization dates that verify with
the available observations at each lead. In turn, initializations and verification dates
are not held constant for each lead.

2. Use the identical set of initializations that can verify over the given observational
window at all leads. However, the verification dates change at each lead.

3. Use the identical verification window at each lead, while allowing the set of
initializations used at each lead to change.

These strategies are shown graphically and explained in more
detail in the documentation. Note that `climpred` offers extensive analysis functionality
in addition to forecast verification, such as spatiotemporal smoothing [@Goddard:2013],
bias removal [@Boer:2016], significance testing [@Goddard:2013; @Boer:2016; @DelSole:2016],
and a graphics library.

<!-- 58 words -->
# Use in Academic Literature
`climpred` has been used to drive analysis in three academic papers so far. @Brady:2020
used the `HindcastEnsemble` class to highlight multi-year predictability of ocean
acidification in the California Current; @Spring:2020 and @Spring:2021 used the
`PerfectModelEnsemble` class to highlight predictability horizons in the global carbon
cycle; and @Krumhardt:2020 used the `HindcastEnsemble` class to illuminate multi-year
predictability in marine Net Primary Productivity.

# Acknowledgements
We thank Andrew Huang for early stage refactoring and continued feedback on
`climpred`. We also thank Kathy Pegion for pioneering the seasonal, monthly,
and subseasonal time resolutions. Thanks in addition to Ray Bell for
initiating and maintaining `xskillscore`, which serves to host the majority of metrics
used in `climpred`.

# References
https://github.com/pydata/xarray/blob/main/asv_bench/benchmarks/README_CI.md
==========
What's New
==========

.. ipython:: python
    :suppress:

    import climpred
    from climpred import HindcastEnsemble
    import matplotlib as mpl

    mpl.rcdefaults()
    mpl.use("Agg")
    # cut border when saving (for maps)
    mpl.rcParams["savefig.bbox"] = "tight"


climpred unreleased
===================

Internals/Minor Fixes
---------------------
- Refactor ``asv`` benchmarking. Add ``run-benchmarks`` label to ``PR`` to run ``asv``
  via Github Actions. (:issue:`664`, :pr:`718`) `Aaron Spring`_.
- Remove ``ipython`` from ``requirements.txt``. (:pr:`720`) `Aaron Spring`_.


climpred v2.2.0 (2021-12-20)
============================

Bug Fixes
---------
- Fix when creating ``valid_time`` from ``lead.attrs["units"]`` in
  ``["seasons", "years"]`` with multi-month stride in ``init``.
  (:issue:`698`, :pr:`700`) `Aaron Spring`_.
- Fix ``seasonality="season"`` in ``reference="climatology"``.
  (:issue:`641`, :pr:`703`) `Aaron Spring`_.

New Features
------------
- Upon instantiation, :py:class:`.PredictionEnsemble` generates new
  2-dimensional coordinate ``valid_time`` for ``initialized`` from ``init`` and
  ``lead``, which is matched with ``time`` from ``verification`` during alignment.
  (:issue:`575`, :pr:`675`, :pr:`678`) `Aaron Spring`_.

.. :: python

>>> hind = climpred.tutorial.load_dataset("CESM-DP-SST")
>>> hind.lead.attrs["units"] = "years"
>>> climpred.HindcastEnsemble(hind).get_initialized()
<xarray.Dataset>
Dimensions:     (lead: 10, member: 10, init: 64)
Coordinates:
  * lead        (lead) int32 1 2 3 4 5 6 7 8 9 10
  * member      (member) int32 1 2 3 4 5 6 7 8 9 10
  * init        (init) object 1954-01-01 00:00:00 ... 2017-01-01 00:00:00
    valid_time  (lead, init) object 1955-01-01 00:00:00 ... 2027-01-01 00:00:00
Data variables:
    SST         (init, lead, member) float64 ...

- Allow ``lead`` as ``float`` also if ``calendar="360_day"`` or ``lead.attrs["units"]``
  not in ``["years","seasons","months"]``. (:issue:`564`, :pr:`675`) `Aaron Spring`_.
- Implement :py:meth:`.HindcastEnsemble.generate_uninitialized` resampling years
  without replacement from ``initialized``. (:issue:`589`, :pr:`591`) `Aaron Spring`_.
- Implement Logarithmic Ensemble Skill Score :py:func:`~climpred.metrics._less`.
  (:issue:`239`, :pr:`687`) `Aaron Spring`_.
- :py:meth:`.HindcastEnsemble.remove_seasonality` and
  :py:meth:`.PerfectModelEnsemble.remove_seasonality` remove the
  seasonality of all ``climpred`` datasets. (:issue:`530`, :pr:`688`) `Aaron Spring`_.
- Add keyword ``groupby`` in :py:meth:`.HindcastEnsemble.verify`,
  :py:meth:`.PerfectModelEnsemble.verify`, :py:meth:`.HindcastEnsemble.bootstrap` and
  :py:meth:`.PerfectModelEnsemble.bootstrap` to group skill by
  initializations seasonality. (:issue:`635`, :pr:`690`) `Aaron Spring`_.


.. :: python

>>> import climpred
>>> hind = climpred.tutorial.load_dataset("NMME_hindcast_Nino34_sst")
>>> obs = climpred.tutorial.load_dataset("NMME_OIv2_Nino34_sst")
>>> hindcast = climpred.HindcastEnsemble(hind).add_observations(obs)
>>> # skill for each init month separated
>>> skill = hindcast.verify(
...     metric="rmse",
...     dim="init",
...     comparison="e2o",
...     skipna=True,
...     alignment="maximize",
...     groupby="month",
... )
>>> skill
<xarray.Dataset>
Dimensions:  (month: 12, lead: 12, model: 12)
Coordinates:
  * lead     (lead) float64 0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0
  * model    (model) object 'NCEP-CFSv2' 'NCEP-CFSv1' ... 'GEM-NEMO'
    skill    <U11 'initialized'
  * month    (month) int64 1 2 3 4 5 6 7 8 9 10 11 12
Data variables:
    sst      (month, lead, model) float64 0.4127 0.3837 0.3915 ... 1.255 3.98
>>> skill.sst.plot(hue="model", col="month", col_wrap=3)

- :py:meth:`.HindcastEnsemble.plot_alignment` shows how forecast and
  observations are aligned based on the `alignment <alignment.html>`_ keyword.
  This may help understanding which dates are matched for the different ``alignment``
  approaches. (:issue:`701`, :pr:`702`) `Aaron Spring`_.

  .. ipython:: python
      :okwarning:
      :okexcept:

      from climpred.tutorial import load_dataset

      hindcast = climpred.HindcastEnsemble(
          load_dataset("CESM-DP-SST")
      ).add_observations(load_dataset("ERSST"))
      @savefig plot_alignment_example.png width=100%
      hindcast.plot_alignment(edgecolor="w")

- Add ``attrs`` to new ``coordinates`` created by ``climpred``.
  (:issue:`695`, :pr:`697`) `Aaron Spring`_.
- Add ``seasonality="weekofyear"`` in ``reference="climatology"``.
  (:pr:`703`) `Aaron Spring`_.
- Compute ``reference="persistence"`` in
  :py:class:`.PerfectModelEnsemble` from ``initialized`` first ``lead``
  if :py:class:`~climpred.options.set_options`
  ``(PerfectModel_persistence_from_initialized_lead_0=True)`` (``False`` by default)
  using :py:func:`~climpred.reference.compute_persistence_from_first_lead`.
  (:issue:`637`, :pr:`706`) `Aaron Spring`_.


Internals/Minor Fixes
---------------------
- Reduce dependencies. (:pr:`686`) `Aaron Spring`_.
- Add `typing <https://docs.python.org/3/library/typing.html>`_.
  (:issue:`685`, :pr:`692`) `Aaron Spring`_.
- refactor ``add_attrs`` into :py:meth:`.HindcastEnsemble.verify` and
  :py:meth:`.HindcastEnsemble.bootstrap`. Now all keywords are
  captured in the skill dataset attributes ``.attrs``.
  (:issue:`475`, :pr:`694`) `Aaron Spring`_.
- docstrings formatting with `blackdocs <https://github.com/keewis/blackdoc>`_.
  (:pr:`708`) `Aaron Spring`_.

Documentation
-------------
- Refresh all docs with ``sphinx_book_theme`` and ``myst_nb``.
  (:issue:`707`, :pr:`708`, :pr:`709`, :pr:`710`) `Aaron Spring`_.


climpred v2.1.6 (2021-08-31)
============================

Adding on to ``v2.1.5``, more bias reduction methods wrapped from
`xclim <https://xclim.readthedocs.io/en/latest/sdba.html>`__
are implemented.

Bug Fixes
---------
- Fix ``results="p"`` in :py:meth:`.HindcastEnsemble.bootstrap` and
  :py:meth:`.PerfectModelEnsemble.bootstrap` when
  ``reference='climatology'``.
  (:issue:`668`, :pr:`670`) `Aaron Spring`_.
- :py:meth:`.HindcastEnsemble.remove_bias` for ``how`` in
  ``["modified_quantile", "basic_quantile", "gamma_mapping", "normal_mapping"]``
  from `bias_correction <https://github.com/pankajkarman/bias_correction>`__
  takes all ``member`` to create model distribution. (:pr:`667`) `Aaron Spring`_.

New Features
------------
- allow more `bias reduction <bias_removal.html>`_ methods wrapped from
  `xclim <https://xclim.readthedocs.io/en/stable/sdba_api.html>`__ in
  :py:meth:`.HindcastEnsemble.remove_bias`:

    * ``how="EmpiricalQuantileMapping"``:
      :py:class:`xclim.sdba.adjustment.EmpiricalQuantileMapping`
    * ``how="DetrendedQuantileMapping"``:
      :py:class:`xclim.sdba.adjustment.DetrendedQuantileMapping`
    * ``how="PrincipalComponents"``:
      :py:class:`xclim.sdba.adjustment.PrincipalComponents`
    * ``how="QuantileDeltaMapping"``:
      :py:class:`xclim.sdba.adjustment.QuantileDeltaMapping`
    * ``how="Scaling"``: :py:class:`xclim.sdba.adjustment.Scaling`
    * ``how="LOCI"``: :py:class:`xclim.sdba.adjustment.LOCI`

  These methods do not respond to ``OPTIONS['seasonality']`` like the other methods.
  Provide ``group="init.month"`` to group by month or ``group='init'`` to skip grouping.
  Provide ``group=None`` or skip ``group`` to use ``init.{OPTIONS['seasonality']}``.
  (:issue:`525`, :pr:`662`, :pr:`666`, :pr:`671`) `Aaron Spring`_.


climpred v2.1.5 (2021-08-12)
============================

While ``climpred`` has used in the
`ASP summer colloquium 2021 <https://asp.ucar.edu/asp-colloquia>`_,
many new features in :py:meth:`.HindcastEnsemble.remove_bias` were
implemented.

Breaking changes
----------------
- renamed ``cross_validate`` to ``cv=False`` in
  :py:meth:`.HindcastEnsemble.remove_bias`.
  Only used when ``train_test_split='unfair-cv'``.
  (:issue:`648`, :pr:`655`). `Aaron Spring`_.

Bug Fixes
---------
- Shift back ``init`` by ``lead`` after
  :py:meth:`.HindcastEnsemble.verify`.
  (:issue:`644`, :pr:`645`) `Aaron Spring`_.

New Features
------------
- :py:meth:`.HindcastEnsemble.remove_bias` accepts new keyword
  ``train_test_split='fair/unfair/unfair-cv'`` (default ``unfair``) following
  `Risbey et al. 2021 <http://www.nature.com/articles/s41467-021-23771-z>`_.
  (:issue:`648`, :pr:`655`) `Aaron Spring`_.
- allow more `bias reduction <bias_removal.html>`_ methods in
  :py:meth:`.HindcastEnsemble.remove_bias`:

    * ``how="additive_mean"``: correcting the mean forecast additively
      (already implemented)
    * ``how="multiplicative_mean"``: correcting the mean forecast multiplicatively
    * ``how="multiplicative_std"``: correcting the standard deviation multiplicatively

  Wrapped from `bias_correction <https://github.com/pankajkarman/bias_correction/blob/master/bias_correction.py>`__:

    * ``how="modified_quantile"``: `Bai et al. 2016 <https://www.sciencedirect.com/science/article/abs/pii/S0034425716302000?via%3Dihub>`_
    * ``how="basic_quantile"``: `Themeßl et al. 2011 <https://rmets.onlinelibrary.wiley.com/doi/pdf/10.1002/joc.2168>`_
    * ``how="gamma_mapping"`` and ``how="normal_mapping"``: `Switanek et al. 2017 <https://www.hydrol-earth-syst-sci.net/21/2649/2017/>`_

- :py:meth:`.HindcastEnsemble.remove_bias` now does
  `leave-one-out cross validation <https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.LeaveOneOut.html>`_
  when passing ``cv='LOO'`` and ``train_test_split='unfair-cv'``.
  ``cv=True`` falls  back to ``cv='LOO'``. (:issue:`643`, :pr:`646`) `Aaron Spring`_.
- Add new metrics :py:func:`~climpred.metrics._spread` and
  :py:func:`~climpred.metrics._mul_bias` (:pr:`638`) `Aaron Spring`_.
- Add new tutorial datasets: (:pr:`651`) `Aaron Spring`_.

    * ``NMME_OIv2_Nino34_sst`` and ``NMME_hindcast_Nino34_sst`` with monthly leads
    * ``Observations_Germany`` and ``ECMWF_S2S_Germany`` with daily leads

- Metadata from `CF convenctions <http://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html>`_
  are automatically attached by
  `cf_xarray <https://cf-xarray.readthedocs.io/en/latest/generated/xarray.DataArray.cf.add_canonical_attributes.html#xarray.DataArray.cf.add_canonical_attributes>`_.
  (:issue:`639`, :pr:`656`) `Aaron Spring`_.
- Raise warning when dimensions ``time``, ``init`` or ``member`` are chunked to show
  user how to circumvent ``xskillscore`` chunking ``ValueError`` when passing these
  dimensions as ``dim`` in :py:meth:`.HindcastEnsemble.verify` or
  :py:meth:`.HindcastEnsemble.bootstrap`.
  (:issue:`509`, :pr:`658`) `Aaron Spring`_.
- Implement ``PredictionEnsemble.chunks``. (:pr:`658`) `Aaron Spring`_.


Documentation
-------------
- Speed up `ENSO monthly example <examples/monseas/monthly-enso-subx-example.ipynb>`_
  with IRIDL server-side preprocessing
  (see `context <https://twitter.com/realaaronspring/status/1406980080883150848?s=21>`_)
  (:issue:`594`, :pr:`633`) `Aaron Spring`_.
- Add `CITATION.cff <https://github.com/pangeo-data/climpred/blob/main/CITATION.cff>`_.
  Please cite
  `Brady and Spring, 2020 <https://joss.theoj.org/papers/10.21105/joss.02781>`_.
  (`GH <https://github.com/pangeo-data/climpred/commit/eceb3f46d78c7dd8eb25243b2e0b673ddd78a4b2>`_) `Aaron Spring`_.
- Use ``NMME_OIv2_Nino34_sst`` and ``NMME_hindcast_Nino34_sst`` with monthly leads for
  `bias reduction <bias_removal.html>`_ demonstrating
  :py:meth:`.HindcastEnsemble.remove_bias`.
  (:pr:`646`) `Aaron Spring`_.


climpred v2.1.4 (2021-06-28)
============================

New Features
------------
- Allow ``hours``, ``minutes`` and ``seconds`` as ``lead.attrs['units']``.
  (:issue:`404`, :pr:`603`) `Aaron Spring`_.
- Allow to set ``seasonality`` via :py:class:`~climpred.options.set_options` to specify
  how to group in ``verify(reference='climatology'`` or in
  :py:meth:`.HindcastEnsemble.remove_bias`.
  (:issue:`529`, :pr:`593`, :pr:`603`) `Aaron Spring`_.
- Allow ``weekofyear`` via ``datetime`` in
  :py:meth:`.HindcastEnsemble.remove_bias`, but not yet implemented in
  ``verify(reference='climatology')``. (:issue:`529`, :pr:`603`) `Aaron Spring`_.
- Allow more dimensions in ``initialized`` than in ``observations``. This is particular
  useful if you have forecasts from multiple models (in a ``model`` dimension) and want
  to verify against the same observations.
  (:issue:`129`, :issue:`528`, :pr:`619`) `Aaron Spring`_.
- Automatically rename dimensions to ``CLIMPRED_ENSEMBLE_DIMS``
  [``"init"``, ``"member"``, ``"lead"``] if CF standard_names in coordinate attributes
  match: (:issue:`613`, :pr:`622`) `Aaron Spring`_.

    * ``"init"``: ``"forecast_reference_time"``
    * ``"member"``: ``"realization"``
    * ``"lead"``: ``"forecast_period"``
- If ``lead`` coordinate is ``pd.Timedelta``,
  :py:class:`.PredictionEnsemble` converts ``lead`` coordinate upon
  instantiation to integer ``lead`` and corresponding ``lead.attrs["units"]``.
  (:issue:`606`, :pr:`627`) `Aaron Spring`_.
- Require ``xskillscore >= 0.0.20``.
  :py:func:`~climpred.metrics._rps` now works with different ``category_edges``
  for observations and forecasts, see
  `daily ECMWF example <examples/subseasonal/daily-S2S-ECMWF.html#biweekly-aggregates>`_.
  (:issue:`629`, :pr:`630`) `Aaron Spring`_.
- Set options ``warn_for_failed_PredictionEnsemble_xr_call``,
  ``warn_for_rename_to_climpred_dims``, ``warn_for_init_coords_int_to_annual``,
  ``climpred_warnings`` via :py:class:`~climpred.options.set_options`.
  (:issue:`628`, :pr:`631`) `Aaron Spring`_.
- :py:class:`.PredictionEnsemble` acts like
  :py:class:`xarray.Dataset` and understands ``data_vars``, ``dims``, ``sizes``,
  ``coords``, ``nbytes``, ``equals``, ``identical``, ``__iter__``, ``__len__``,
  ``__contains__``, ``__delitem__``. (:issue:`568`, :pr:`632`) `Aaron Spring`_.


Documentation
-------------
- Add `documentation page about publicly available initialized datasets and
  corresponding `climpred` examples <initialized-datasets.html>`_.
  (:issue:`510`, :issue:`561`, :pr:`600`) `Aaron Spring`_.
- Add `GEFS example <examples/NWP/NWP_GEFS_6h_forecasts.html>`_ for numerical weather
  prediction. (:issue:`602`, :pr:`603`) `Aaron Spring`_.
- Add subseasonal `daily ECMWF example <examples/subseasonal/daily-S2S-ECMWF.html>`__
  using `climetlab <https://github.com/ecmwf-lab/climetlab-s2s-ai-challenge>`_ to access
  hindcasts from ECMWF cloud.  (:issue:`587`, :pr:`603`) `Aaron Spring`_.
- Add subseasonal `daily S2S example <examples/subseasonal/daily-S2S-IRIDL.html>`_
  accessing `S2S <http://s2sprediction.net/>`_ output on
  `IRIDL <https://iridl.ldeo.columbia.edu/SOURCES/.ECMWF/.S2S/>`_ with a cookie and
  working with "on-the-fly" reforecasts with ``hdate`` dimension.
  (:issue:`588`, :pr:`593`) `Aaron Spring`_.
- Added example `climpred on GPU <examples/misc/climpred_gpu.ipynb>`_. Running
  :py:meth:`.PerfectModelEnsemble.verify` on GPU with `cupy-xarray
  <https://github.com/xarray-contrib/cupy-xarray>`_ finishes 10x faster.
  (:issue:`592`, :pr:`607`) `Aaron Spring`_.
- How to work with biweekly aggregates in ``climpred``, see
  `daily ECMWF example <examples/subseasonal/daily-S2S-ECMWF.html#biweekly-aggregates>`__.
  (:issue:`625`, :pr:`630`) `Aaron Spring`_.


Internals/Minor Fixes
---------------------
- Add weekly upstream CI, which raises issues for failures. Adapted from ``xarray``.
  Manually trigger by ``git commit -m '[test-upstream]'``. Skip climpred_testing CI by
  ``git commit -m '[skip-ci]'``
  (:issue:`518`, :pr:`596`) `Aaron Spring`_.


climpred v2.1.3 (2021-03-23)
============================

Breaking changes
----------------

New Features
------------
- :py:meth:`.HindcastEnsemble.verify`,
  :py:meth:`.PerfectModelEnsemble.verify`,
  :py:meth:`.HindcastEnsemble.bootstrap` and
  :py:meth:`.PerfectModelEnsemble.bootstrap`
  accept reference ``climatology``. Furthermore, reference ``persistence`` also allows
  probabilistic metrics (:issue:`202`, :issue:`565`, :pr:`566`) `Aaron Spring`_.
- Added new metric  :py:class:`~climpred.metrics._roc` Receiver Operating
  Characteristic as ``metric='roc'``. (:pr:`566`) `Aaron Spring`_.

Bug fixes
---------
- :py:meth:`.HindcastEnsemble.verify` and
  :py:meth:`.HindcastEnsemble.bootstrap` accept ``dim`` as ``list``,
  ``set``, ``tuple`` or ``str`` (:issue:`519`, :pr:`558`) `Aaron Spring`_.
- :py:meth:`.PredictionEnsemble.map` now does not fail silently when
  applying a function to all ``xr.Datasets`` of
  :py:class:`.PredictionEnsemble`. Instead, ``UserWarning``s are
  raised. Furthermore, ``PredictionEnsemble.map(func, *args, **kwargs)``
  applies only function to Datasets with matching dims if ``dim="dim0_or_dim1"`` is
  passed as ``**kwargs``. (:issue:`417`, :issue:`437`, :pr:`552`) `Aaron Spring`_.
- :py:class:`~climpred.metrics._rpc` was fixed in ``xskillscore>=0.0.19`` and hence is
  not falsely limited to 1 anymore (:issue:`562`, :pr:`566`) `Aaron Spring`_.

Internals/Minor Fixes
---------------------
- Docstrings are now tested in GitHub actions continuous integration.
  (:issue:`545`, :pr:`560`) `Aaron Spring`_.
- Github actions now cancels previous commits, instead of running the full
  testing suite on every single commit. (:pr:`560`) `Aaron Spring`_.
- :py:meth:`.PerfectModelEnsemble.verify` does not add
  climpred attributes to skill by default anymore.
  (:pr:`560`) `Aaron Spring`_.
- Drop ``python==3.6`` support. (:pr:`573`) `Aaron Spring`_.
- Notebooks are now linted with
  `nb_black <https://github.com/dnanhkhoa/nb_black>`_ using
  ``%load_ext nb_black`` or ``%load_ext lab_black`` for
  `Jupyter <https://jupyter.org>`_ notebooks and
  `Jupyter <https://jupyter.org>`_ lab.
  (:issue:`526`, :pr:`572`) `Aaron Spring`_.
- Reduce dependencies to install climpred.
  (:issue:`454`, :pr:`572`) `Aaron Spring`_.
- Examples from documentation available via `Binder <https://mybinder.org/v2/gh/pangeo-data/climpred/master?urlpath=lab%2Ftree%2Fdocs%2Fsource%2Fquick-start.ipynb>`_.
  Find further examples in the ``examples`` folder.
  (:issue:`549`, :pr:`578`) `Aaron Spring`_.
- Rename branch ``master`` to ``main``. (:pr:`579`) `Aaron Spring`_.


climpred v2.1.2 (2021-01-22)
============================

This release is the fixed version for our Journal of Open Source Software (JOSS)
article about ``climpred``, see `review
<https://github.com/openjournals/joss-reviews/issues/2781>`_.

New Features
------------
- Function to calculate predictability horizon
  :py:func:`~climpred.predictability_horizon.predictability_horizon` based on condition.
  (:issue:`46`, :pr:`521`) `Aaron Spring`_.

Bug fixes
---------
- :py:meth:`.PredictionEnsemble.smooth` now carries ``lead.attrs``
  (:issue:`527`, pr:`521`) `Aaron Spring`_.
- :py:meth:`.PerfectModelEnsemble.verify` now works with ``references``
  also for geospatial inputs, which returned ``NaN`` before.
  (:issue:`522`, pr:`521`) `Aaron Spring`_.
- :py:meth:`.PredictionEnsemble.plot` now shifts composite lead
  frequencies like ``days``, ``pentads``, ``seasons`` correctly.
  (:issue:`532`, :pr:`533`) `Aaron Spring`_.
- Adapt to ``xesmf>=0.5.2`` for spatial xesmf smoothing. (:issue:`543`, :pr:`548`)
  `Aaron Spring`_.
- :py:meth:`.HindcastEnsemble.remove_bias` now carries attributes.
  (:issue:`531`, :pr:`551`) `Aaron Spring`_.


climpred v2.1.1 (2020-10-13)
============================

Breaking changes
----------------

This version introduces a lot of breaking changes. We are trying to overhaul
``climpred`` to have an intuitive API that also forces users to think about methodology
choices when running functions. The main breaking changes we introduced are for
:py:meth:`.HindcastEnsemble.verify` and
:py:meth:`.PerfectModelEnsemble.verify`. Now, instead of assuming
defaults for most keywords, we require the user to define ``metric``, ``comparison``,
``dim``, and ``alignment`` (for hindcast systems). We also require users to designate
the number of ``iterations`` for bootstrapping.

- User now has to designate number of iterations with ``iterations=...`` in
  :py:meth:`.HindcastEnsemble.bootstrap` (:issue:`384`, :pr:`436`)
  `Aaron Spring`_ and `Riley X. Brady`_.
- Make ``metric``, ``comparison``, ``dim``, and ``alignment`` required (previous default
  ``None``) arguments for :py:meth:`.HindcastEnsemble.verify`
  (:issue:`384`, :pr:`436`) `Aaron Spring`_ and `Riley X. Brady`_.
- Metric :py:class:`~climpred.metrics._brier_score` and
  :py:func:`~climpred.metrics._threshold_brier_score` now requires callable keyword
  argument ``logical`` instead of ``func`` (:pr:`388`) `Aaron Spring`_.
- :py:meth:`.HindcastEnsemble.verify` does not correct ``dim``
  automatically to ``member`` for probabilistic metrics.
  (:issue:`282`, :pr:`407`) `Aaron Spring`_.
- Users can no longer add multiple observations to
  :py:class:`.HindcastEnsemble`. This will make current and future
  development much easier on maintainers (:issue:`429`, :pr:`453`) `Riley X. Brady`_.
- Standardize the names of the output coordinates for
  :py:meth:`.PredictionEnsemble.verify` and
  :py:meth:`.PredictionEnsemble.bootstrap` to ``initialized``,
  ``uninitialized``, and ``persistence``. ``initialized`` showcases the metric result
  after comparing the initialized ensemble to the verification data; ``uninitialized``
  when comparing the uninitialized (historical) ensemble to the verification data;
  ``persistence`` is the evaluation of the persistence forecast
  (:issue:`460`, :pr:`478`, :issue:`476`, :pr:`480`) `Aaron Spring`_.
- ``reference`` keyword in :py:meth:`.HindcastEnsemble.verify` should
  be choosen from [``uninitialized``, ``persistence``]. ``historical`` no longer works.
  (:issue:`460`, :pr:`478`, :issue:`476`, :pr:`480`) `Aaron Spring`_.
- :py:meth:`.HindcastEnsemble.verify` returns no ``skill`` dimension
  if ``reference=None``  (:pr:`480`) `Aaron Spring`_.
- ``comparison`` is not applied to uninitialized skill in
  :py:meth:`.HindcastEnsemble.bootstrap`.
  (:issue:`352`, :pr:`418`) `Aaron Spring`_.

New Features
------------

This release is accompanied by a bunch of new features. Math operations can now be used
with our :py:class:`.PredictionEnsemble` objects and their variables
can be sub-selected. Users can now quick plot time series forecasts with these objects.
Bootstrapping is available for :py:class:`.HindcastEnsemble`. Spatial
dimensions can be passed to metrics to do things like pattern correlation. New metrics
have been implemented based on Contingency tables. We now include an early version
of bias removal for :py:class:`.HindcastEnsemble`.

- Use math operations like ``+-*/`` with :py:class:`.HindcastEnsemble`
  and :py:class:`.PerfectModelEnsemble`. See
  `demo <prediction-ensemble-object.html>`_
  Arithmetic-Operations-with-PredictionEnsemble-Objects. (:pr:`377`) `Aaron Spring`_.
- Subselect data variables from :py:class:`.PerfectModelEnsemble` as
  from :py:class:`xarray.Dataset`:
  ``PredictionEnsemble[["var1", "var3"]]`` (:pr:`409`) `Aaron Spring`_.
- Plot all datasets in :py:class:`.HindcastEnsemble` or
  :py:class:`.PerfectModelEnsemble` by
  :py:meth:`.PredictionEnsemble.plot` if no other spatial dimensions
  are present. (:pr:`383`) `Aaron Spring`_.
- Bootstrapping now available for :py:class:`.HindcastEnsemble` as
  :py:meth:`.HindcastEnsemble.bootstrap`, which is analogous to
  the :py:class:`.PerfectModelEnsemble` method.
  (:issue:`257`, :pr:`418`) `Aaron Spring`_.
- :py:meth:`.HindcastEnsemble.verify` allows all dimensions from
  ``initialized`` ensemble as ``dim``. This allows e.g. spatial dimensions to be used
  for pattern correlation. Make sure to use ``skipna=True`` when using spatial
  dimensions and output has NaNs (in the case of land, for instance).
  (:issue:`282`, :pr:`407`) `Aaron Spring`_.
- Allow binary forecasts at when calling
  :py:meth:`.HindcastEnsemble.verify`,
  rather than needing to supply binary results beforehand. In other words,
  ``hindcast.verify(metric='bs', comparison='m2o', dim='member', logical=logical)``
  is now the same as
  ``hindcast.map(logical).verify(metric='brier_score', comparison='m2o', dim='member'``.
  (:pr:`431`) `Aaron Spring`_.
- Check ``calendar`` types when using
  :py:meth:`.HindcastEnsemble.add_observations`,
  :py:meth:`.HindcastEnsemble.add_uninitialized`,
  :py:meth:`.PerfectModelEnsemble.add_control` to ensure that the
  verification data calendars match that of the initialized ensemble.
  (:issue:`300`, :pr:`452`, :issue:`422`, :pr:`462`)
  `Riley X. Brady`_ and `Aaron Spring`_.
- Implement new metrics which have been ported over from
  https://github.com/csiro-dcfp/doppyo/ to ``xskillscore`` by `Dougie Squire`_.
  (:pr:`439`, :pr:`456`) `Aaron Spring`_

    * rank histogram :py:func:`~climpred.metrics._rank_histogram`
    * discrimination :py:func:`~climpred.metrics._discrimination`
    * reliability :py:func:`~climpred.metrics._reliability`
    * ranked probability score :py:func:`~climpred.metrics._rps`
    * contingency table and related scores :py:func:`~climpred.metrics._contingency`

- Perfect Model :py:meth:`.PerfectModelEnsemble.verify`
  no longer requires ``control`` in :py:class:`.PerfectModelEnsemble`.
  It is only required when ``reference=['persistence']``. (:pr:`461`) `Aaron Spring`_.
- Implemented bias removal
  :py:class:`~climpred.classes.HindcastEnsemble.remove_bias`.
  ``remove_bias(how='mean')`` removes the mean bias of initialized hindcasts with
  respect to observations. See `example <bias_removal.html>`__.
  (:pr:`389`, :pr:`443`, :pr:`459`) `Aaron Spring`_ and `Riley X. Brady`_.

Deprecated
----------

- ``spatial_smoothing_xrcoarsen`` no longer used for spatial smoothing.
  (:pr:`391`) `Aaron Spring`_.
- ``compute_metric``, ``compute_uninitialized`` and ``compute_persistence`` no longer
  in use for :py:class:`.PerfectModelEnsemble` in favor of
  :py:meth:`.PerfectModelEnsemble.verify` with the ``reference``
  keyword instead.
  (:pr:`436`, :issue:`468`, :pr:`472`) `Aaron Spring`_ and `Riley X. Brady`_.
- ``'historical'`` no longer a valid choice for ``reference``. Use ``'uninitialized'``
  instead. (:pr:`478`) `Aaron Spring`_.

Bug Fixes
---------

- :py:meth:`.PredictionEnsemble.verify` and
  :py:meth:`.PredictionEnsemble.bootstrap` now accept
  ``metric_kwargs``. (:pr:`387`) `Aaron Spring`_.
- :py:meth:`.PerfectModelEnsemble.verify` now accepts
  ``'uninitialized'`` as a reference. (:pr:`395`) `Riley X. Brady`_.
- Spatial and temporal smoothing :py:meth:`.PredictionEnsemble.smooth`
  now work as expected and rename time dimensions after
  :py:meth:`~climpred.classes.PredictionEnsembleEnsemble.verify`.
  (:pr:`391`) `Aaron Spring`_.
- ``PredictionEnsemble.verify(comparison='m2o', references=['uninitialized',
  'persistence']`` does not fail anymore. (:issue:`385`, :pr:`400`) `Aaron Spring`_.
- Remove bias using ``dayofyear`` in
  :py:meth:`.HindcastEnsemble.reduce_bias`.
  (:pr:`443`) `Aaron Spring`_.
- ``climpred`` works with ``dask=>2.28``. (:issue:`479`, :pr:`482`) `Aaron Spring`_.

Documentation
-------------
- Updates ``climpred`` tagline to "Verification of weather and climate forecasts."
  (:pr:`420`) `Riley X. Brady`_.
- Adds section on how to use arithmetic with
  :py:class:`.HindcastEnsemble`.
  (:pr:`378`) `Riley X. Brady`_.
- Add docs section for similar open-source forecasting packages.
  (:pr:`432`) `Riley X. Brady`_.
- Add all metrics to main API in addition to metrics page.
  (:pr:`438`) `Riley X. Brady`_.
- Add page on bias removal `Aaron Spring`_.

Internals/Minor Fixes
---------------------
- :py:meth:`.PredictionEnsemble.verify` replaces deprecated
  ``PerfectModelEnsemble.compute_metric()`` and accepts ``reference`` as keyword.
  (:pr:`387`) `Aaron Spring`_.
- Cleared out unnecessary statistics functions from ``climpred`` and migrated them to
  ``esmtools``. Add ``esmtools`` as a required package. (:pr:`395`) `Riley X. Brady`_.
- Remove fixed pandas dependency from ``pandas=0.25`` to stable ``pandas``.
  (:issue:`402`, :pr:`403`) `Aaron Spring`_.
- ``dim`` is expected to be a list of strings in
  :py:func:`~climpred.prediction.compute_perfect_model` and
  :py:func:`~climpred.prediction.compute_hindcast`.
  (:issue:`282`, :pr:`407`) `Aaron Spring`_.
- Update ``cartopy`` requirement to 0.0.18 or greater to release lock on
  ``matplotlib`` version. Update ``xskillscore`` requirement to 0.0.18 to
  cooperate with new ``xarray`` version. (:pr:`451`, :pr:`449`)
  `Riley X. Brady`_
- Switch from Travis CI and Coveralls to Github Actions and CodeCov.
  (:pr:`471`) `Riley X. Brady`_
- Assertion functions added for :py:class:`.PerfectModelEnsemble`:
  :py:func:`~climpred.testing.assert_PredictionEnsemble`. (:pr:`391`) `Aaron Spring`_.
- Test all metrics against synthetic data. (:pr:`388`) `Aaron Spring`_.


climpred v2.1.0 (2020-06-08)
============================

Breaking Changes
----------------

- Keyword ``bootstrap`` has been replaced with ``iterations``. We feel that this more
  accurately describes the argument, since "bootstrap" is really the process as a whole.
  (:pr:`354`) `Aaron Spring`_.

New Features
------------

- :py:class:`.HindcastEnsemble` and
  :py:class:`.PerfectModelEnsemble` now use an HTML representation,
  following the more recent versions of ``xarray``. (:pr:`371`) `Aaron Spring`_.
- ``HindcastEnsemble.verify()`` now takes ``reference=...`` keyword. Current options are
  ``'persistence'`` for a persistence forecast of the observations and
  ``'uninitialized'`` for an uninitialized/historical reference, such as an
  uninitialized/forced run. (:pr:`341`) `Riley X. Brady`_.
- We now only enforce a union of the initialization dates with observations if
  ``reference='persistence'`` for :py:class:`.HindcastEnsemble`.
  This is to ensure that the same set of initializations is used by the observations to
  construct a persistence forecast. (:pr:`341`) `Riley X. Brady`_.
- :py:func:`~climpred.prediction.compute_perfect_model` now accepts initialization
  (``init``) as ``cftime`` and ``int``. ``cftime`` is now implemented into the
  bootstrap uninitialized functions for the perfect model configuration.
  (:pr:`332`) `Aaron Spring`_.
- New explicit keywords in bootstrap functions for ``resampling_dim`` and
  ``reference_compute`` (:pr:`320`) `Aaron Spring`_.
- Logging now included for ``compute_hindcast`` which displays the ``inits`` and
  verification dates used at each lead (:pr:`324`) `Aaron Spring`_,
  (:pr:`338`) `Riley X. Brady`_. See (`logging <alignment.html#Logging>`__).
- New explicit keywords added for ``alignment`` of verification dates and
  initializations. (:pr:`324`) `Aaron Spring`_. See (`alignment <alignment.html>`__)

    * ``'maximize'``: Maximize the degrees of freedom by slicing ``hind`` and
      ``verif`` to a common time frame at each lead. (:pr:`338`) `Riley X. Brady`_.
    * ``'same_inits'``: slice to a common init frame prior to computing
      metric. This philosophy follows the thought that each lead should be
      based on the same set of initializations. (:pr:`328`) `Riley X. Brady`_.
    * ``'same_verifs'``: slice to a common/consistent verification time frame prior
      to computing metric. This philosophy follows the thought that each lead
      should be based on the same set of verification dates. (:pr:`331`)
      `Riley X. Brady`_.

Performance
-----------

The major change for this release is a dramatic speedup in bootstrapping functions, led
by `Aaron Spring`_. We focused on scalability with ``dask`` and found many places we
could compute skill simultaneously over all bootstrapped ensemble members rather than
at each iteration.

- Bootstrapping uninitialized skill in the perfect model framework is now sped up
  significantly for annual lead resolution. (:pr:`332`) `Aaron Spring`_.
- General speedup in :py:func:`~climpred.bootstrap.bootstrap_hindcast` and
  :py:func:`~climpred.bootstrap.bootstrap_perfect_model`: (:pr:`285`) `Aaron Spring`_.

    * Properly implemented handling for lazy results when inputs are chunked.

    * User gets warned when chunking potentially unnecessarily and/or inefficiently.

Bug Fixes
---------
- Alignment options now account for differences in the historical time series if
  ``reference='historical'``. (:pr:`341`) `Riley X. Brady`_.

Internals/Minor Fixes
---------------------
- Added a `Code of Conduct <code_of_conduct.html>`__ (:pr:`285`) `Aaron Spring`_.
- Gather ``pytest.fixture in ``conftest.py``. (:pr:`313`) `Aaron Spring`_.
- Move ``x_METRICS`` and ``COMPARISONS`` to ``metrics.py`` and ``comparisons.py`` in
  order to avoid circular import dependencies. (:pr:`315`) `Aaron Spring`_.
- ``asv`` benchmarks added for :py:class:`.HindcastEnsemble`
  (:pr:`285`) `Aaron Spring`_.
- Ignore irrelevant warnings in ``pytest`` and mark slow tests
  (:pr:`333`) `Aaron Spring`_.
- Default ``CONCAT_KWARGS`` now in all ``xr.concat`` to speed up bootstrapping.
  (:pr:`330`) `Aaron Spring`_.
- Remove ``member`` coords for ``m2c`` comparison for probabilistic metrics.
  (:pr:`330`) `Aaron Spring`_.
- Refactored :py:func:`~climpred.prediction.compute_hindcast` and
  :py:func:`~climpred.prediction.compute_perfect_model`. (:pr:`330`) `Aaron Spring`_.
- Changed lead0 coordinate modifications to be compliant with ``xarray=0.15.1`` in
  :py:func:`~climpred.reference.compute_persistence`. (:pr:`348`) `Aaron Spring`_.
- Exchanged ``my_quantile`` with ``xr.quantile(skipna=False)``.
  (:pr:`348`) `Aaron Spring`_.
- Remove ``sig`` from
  :py:func:`~climpred.graphics.plot_bootstrapped_skill_over_leadyear`.
  (:pr:`351`) `Aaron Spring`_.
- Require ``xskillscore v0.0.15`` and use their functions for effective sample
  size-based metrics. (:pr: `353`) `Riley X. Brady`_.
- Faster bootstrapping without replacement used in threshold functions of
  ``climpred.stats`` (:pr:`354`) `Aaron Spring`_.
- Require ``cftime v1.1.2``, which modifies their object handling to create 200-400x
  speedups in some basic operations. (:pr:`356`) `Riley X. Brady`_.
- Resample first and then calculate skill in
  :py:func:`~climpred.bootstrap.bootstrap_perfect_model` and
  :py:func:`~climpred.bootstrap.bootstrap_hindcast` (:pr:`355`) `Aaron Spring`_.

Documentation
-------------
- Added demo to setup your own raw model output compliant to ``climpred``
  (:pr:`296`) `Aaron Spring`_. See (`here <examples/misc/setup_your_own_data.html>`__).
- Added demo using ``intake-esm`` with ``climpred``.
  See `demo <examples/misc/setup_your_own_data.html#intake-esm-for-cmorized-output>`__.
  (:pr:`296`) `Aaron Spring`_.
- Added `Verification Alignment <alignment.html>`_ page explaining how initializations
  are selected and aligned with verification data. (:pr:`328`) `Riley X. Brady`_.
  See (`here <alignment.html>`__).


climpred v2.0.0 (2020-01-22)
============================

New Features
------------
- Add support for ``days``, ``pentads``, ``weeks``, ``months``, ``seasons`` for lead
  time resolution. ``climpred`` now requires a ``lead`` attribute "units" to decipher
  what resolution the predictions are at. (:pr:`294`) `Kathy Pegion`_ and
  `Riley X. Brady`_.

.. :: python

    >>> hind = climpred.tutorial.load_dataset("CESM-DP-SST")
    >>> hind.lead.attrs["units"] = "years"

- :py:class:`.HindcastEnsemble` now has
  :py:meth:`.HindcastEnsemble.add_observations` and
  :py:meth:`.HindcastEnsemble.get_observations`
  methods. These are the same as ``.add_reference()`` and ``.get_reference()``, which
  will be deprecated eventually. The name change clears up confusion, since "reference"
  is the appropriate name for a reference forecast, e.g. ``"persistence"``. (:pr:`310`)
  `Riley X. Brady`_.

- :py:class:`.HindcastEnsemble` now has ``.verify()`` function, which
  duplicates the ``.compute_metric()`` function. We feel that ``.verify()`` is more
  clear and easy to write, and follows the terminology of the field.
  (:pr:`310`) `Riley X. Brady`_.

- ``e2o`` and ``m2o`` are now the preferred keywords for comparing hindcast ensemble
  means and ensemble members to verification data, respectively. (:pr:`310`)
  `Riley X. Brady`_.

Documentation
-------------
- New example pages for subseasonal-to-seasonal prediction using ``climpred``.
  (:pr:`294`) `Kathy Pegion`_

    * Calculate the skill of the MJO index as a function of lead time
      (`link <examples/subseasonal/daily-subx-example.html>`__).

    * Calculate the skill of the MJO index as a function of lead time for weekly data
      (`link <examples/subseasonal/weekly-subx-example.html>`__).

    * Calculate ENSO skill as a function of initial month vs. lead time
      (`link <examples/monseas/monthly-enso-subx-example.html>`__).

    * Calculate Seasonal ENSO skill
      (`link <examples/monseas/seasonal-enso-subx-example.html>`__).

- `Comparisons <comparisons.html>`__ page rewritten for more clarity. (:pr:`310`)
  `Riley X. Brady`_.

Bug Fixes
---------
- Fixed `m2m` broken comparison issue and removed correction.
  (:pr:`290`) `Aaron Spring`_.

Internals/Minor Fixes
---------------------
- Updates to ``xskillscore`` v0.0.12 to get a 30-50% speedup in compute functions that
  rely on metrics from there. (:pr:`309`) `Riley X. Brady`_.
- Stacking dims is handled by ``comparisons``, no need for internal keyword
  ``stack_dims``. Therefore ``comparison`` now takes ``metric`` as argument instead.
  (:pr:`290`) `Aaron Spring`_.
- ``assign_attrs`` now carries `dim` (:pr:`290`) `Aaron Spring`_.
- ``reference`` changed to ``verif`` throughout hindcast compute functions. This is more
  clear, since ``reference`` usually refers to a type of forecast, such as persistence.
  (:pr:`310`) `Riley X. Brady`_.
- ``Comparison`` objects can now have aliases. (:pr:`310`) `Riley X. Brady`_.



climpred v1.2.1 (2020-01-07)
============================

Depreciated
-----------
- ``mad`` no longer a keyword for the median absolute error metric. Users should now
  use ``median_absolute_error``, which is identical to changes in ``xskillscore``
  version 0.0.10. (:pr:`283`) `Riley X. Brady`_
- ``pacc`` no longer a keyword for the p value associated with the Pearson
  product-moment correlation, since it is used by the correlation coefficient.
  (:pr:`283`) `Riley X. Brady`_
- ``msss`` no longer a keyword for the Murphy's MSSS, since it is reserved for the
  standard MSSS. (:pr:`283`) `Riley X. Brady`_

New Features
------------
- Metrics ``pearson_r_eff_p_value`` and ``spearman_r_eff_p_value`` account for
  autocorrelation in computing p values. (:pr:`283`) `Riley X. Brady`_
- Metric ``effective_sample_size`` computes number of independent samples between two
  time series being correlated. (:pr:`283`) `Riley X. Brady`_
- Added keywords for metrics: (:pr:`283`) `Riley X. Brady`_

    * ``'pval'`` for ``pearson_r_p_value``
    * ``['n_eff', 'eff_n']`` for ``effective_sample_size``
    * ``['p_pval_eff', 'pvalue_eff', 'pval_eff']`` for ``pearson_r_eff_p_value``
    * ``['spvalue', 'spval']`` for ``spearman_r_p_value``
    * ``['s_pval_eff', 'spvalue_eff', 'spval_eff']`` for ``spearman_r_eff_p_value``
    * ``'nev'`` for ``nmse``

Internals/Minor Fixes
---------------------
- ``climpred`` now requires ``xarray`` version 0.14.1 so that the ``drop_vars()``
  keyword used in our package does not throw an error. (:pr:`276`) `Riley X. Brady`_
- Update to ``xskillscore`` version 0.0.10 to fix errors in weighted metrics with
  pairwise NaNs. (:pr:`283`) `Riley X. Brady`_
- ``doc8`` added to ``pre-commit`` to have consistent formatting on ``.rst`` files.
  (:pr:`283`) `Riley X. Brady`_
- Remove ``proper`` attribute on ``Metric`` class since it isn't used anywhere.
  (:pr:`283`) `Riley X. Brady`_
- Add testing for effective p values. (:pr:`283`) `Riley X. Brady`_
- Add testing for whether metric aliases are repeated/overwrite each other.
  (:pr:`283`) `Riley X. Brady`_
- ``ppp`` changed to ``msess``, but keywords allow for ``ppp`` and ``msss`` still.
  (:pr:`283`) `Riley X. Brady`_

Documentation
-------------
- Expansion of `metrics documentation <metrics.html>`_ with much more
  detail on how metrics are computed, their keywords, references, min/max/perfect
  scores, etc. (:pr:`283`) `Riley X. Brady`_
- Update `terminology page <terminology.html>`_ with more information on metrics
  terminology. (:pr:`283`) `Riley X. Brady`_


climpred v1.2.0 (2019-12-17)
============================

Depreciated
-----------
- Abbreviation ``pval`` depreciated. Use ``p_pval`` for ``pearson_r_p_value`` instead.
  (:pr:`264`) `Aaron Spring`_.

New Features
------------
- Users can now pass a custom ``metric`` or ``comparison`` to compute functions.
  (:pr:`268`) `Aaron Spring`_.

    * See `user-defined-metrics <metrics.html#user-defined-metrics>`_ and
      `user-defined-comparisons <comparisons.html#user-defined-comparisons>`_.

- New deterministic metrics (see `metrics <metrics.html>`_). (:pr:`264`)
  `Aaron Spring`_.

    * Spearman ranked correlation (spearman_r_)
    * Spearman ranked correlation p-value (spearman_r_p_value_)
    * Mean Absolute Deviation (mad_)
    * Mean Absolute Percent Error (mape_)
    * Symmetric Mean Absolute Percent Error (smape_)

.. _spearman_r: metrics.html#spearman-anomaly-correlation-coefficient-sacc
.. _spearman_r_p_value: metrics.html#spearman-anomaly-correlation-coefficient-sacc
.. _mad: metrics.html#median-absolute-deviation-mad
.. _mape: metrics.html#mean-absolute-percentage-error-mape
.. _smape: metrics.html#symmetric-mean-absolute-percentage-error-smape

- Users can now apply arbitrary ``xarray`` methods to
  :py:class:`.HindcastEnsemble` and
  :py:class:`.PerfectModelEnsemble`. (:pr:`243`) `Riley X. Brady`_.

    * See the
      `Prediction Ensemble objects demo page <prediction-ensemble-object.html>`_.

- Add "getter" methods to :py:class:`.HindcastEnsemble` and
  :py:class:`.PerfectModelEnsemble` to retrieve ``xarray`` datasets
  from the objects. (:pr:`243`) `Riley X. Brady`_.

.. :: python

>>> hind = climpred.tutorial.load_dataset("CESM-DP-SST")
>>> ref = climpred.tutorial.load_dataset("ERSST")
>>> hindcast = climpred.HindcastEnsemble(hind)
>>> hindcast = hindcast.add_reference(ref, "ERSST")
>>> print(hindcast)
<climpred.HindcastEnsemble>
Initialized Ensemble:
    SST      (init, lead, member) float64 ...
ERSST:
    SST      (time) float32 ...
Uninitialized:
    None
>>> print(hindcast.get_initialized())
<xarray.Dataset>
Dimensions:  (init: 64, lead: 10, member: 10)
Coordinates:
* lead     (lead) int32 1 2 3 4 5 6 7 8 9 10
* member   (member) int32 1 2 3 4 5 6 7 8 9 10
* init     (init) float32 1954.0 1955.0 1956.0 1957.0 ... 2015.0 2016.0 2017.0
Data variables:
    SST      (init, lead, member) float64 ...
>>> print(hindcast.get_reference("ERSST"))
<xarray.Dataset>
Dimensions:  (time: 61)
Coordinates:
* time     (time) int64 1955 1956 1957 1958 1959 ... 2011 2012 2013 2014 2015
Data variables:
    SST      (time) float32 ...

- ``metric_kwargs`` can be passed to :py:class:`~climpred.metrics.Metric`.
  (:pr:`264`) `Aaron Spring`_.

    * See ``metric_kwargs`` under `metrics <metrics.html>`_.

Bug Fixes
---------
- :py:meth:`.HindcastEnsemble.compute_metric` doesn't drop coordinates
  from the initialized hindcast ensemble anymore. (:pr:`258`) `Aaron Spring`_.
- Metric ``uacc`` does not crash when ``ppp`` negative anymore. (:pr:`264`)
  `Aaron Spring`_.
- Update ``xskillscore`` to version 0.0.9 to fix all-NaN issue with ``pearson_r`` and
  ``pearson_r_p_value`` when there's missing data. (:pr:`269`) `Riley X. Brady`_.

Internals/Minor Fixes
---------------------
- Rewrote :py:func:`~climpred.stats.varweighted_mean_period` based on ``xrft``.
  Changed ``time_dim`` to ``dim``. Function no longer drops coordinates. (:pr:`258`)
  `Aaron Spring`_
- Add ``dim='time'`` in :py:func:`~climpred.stats.dpp`. (:pr:`258`) `Aaron Spring`_
- Comparisons ``m2m``, ``m2e`` rewritten to not stack dims into supervector because
  this is now done in ``xskillscore``. (:pr:`264`) `Aaron Spring`_
- Add ``tqdm`` progress bar to :py:func:`~climpred.bootstrap.bootstrap_compute`.
  (:pr:`244`) `Aaron Spring`_
- Remove inplace behavior for :py:class:`.HindcastEnsemble` and
  :py:class:`.PerfectModelEnsemble`. (:pr:`243`) `Riley X. Brady`_

    * See `demo page on prediction ensemble objects <prediction-ensemble-object.html>`_

- Added tests for chunking with ``dask``. (:pr:`258`) `Aaron Spring`_
- Fix test issues with esmpy 8.0 by forcing esmpy 7.1 (:pr:`269`). `Riley X. Brady`_
- Rewrote ``metrics`` and ``comparisons`` as classes to accomodate custom metrics and
  comparisons. (:pr:`268`) `Aaron Spring`_

    * See `user-defined-metrics <metrics.html#user-defined-metrics>`_ and
      `user-defined-comparisons <comparisons.html#user-defined-comparisons>`_.

Documentation
-------------
- Add examples notebook for
  `temporal and spatial smoothing <examples/smoothing.html>`_. (:pr:`244`)
  `Aaron Spring`_
- Add documentation for computing a metric over a
  `specified dimension <comparisons.html#compute-over-dimension>`_.
  (:pr:`244`) `Aaron Spring`_
- Update `API <api.html>`_ to be more organized with individual function/class pages.
  (:pr:`243`) `Riley X. Brady`_.
- Add `page <prediction-ensemble-object.html>`_ describing the
  :py:class:`.HindcastEnsemble` and
  :py:class:`.PerfectModelEnsemble` objects more clearly.
  (:pr:`243`) `Riley X. Brady`_
- Add page for `publications <publications.html>`_ and
  `helpful links <helpful-links.html>`_. (:pr:`270`) `Riley X. Brady`_.

climpred v1.1.0 (2019-09-23)
============================

Features
--------
- Write information about skill computation to netcdf attributes(:pr:`213`)
  `Aaron Spring`_
- Temporal and spatial smoothing module (:pr:`224`) `Aaron Spring`_
- Add metrics `brier_score`, `threshold_brier_score` and `crpss_es` (:pr:`232`)
  `Aaron Spring`_
- Allow `compute_hindcast` and `compute_perfect_model` to specify which dimension `dim`
  to calculate metric over (:pr:`232`) `Aaron Spring`_

Bug Fixes
---------
- Correct implementation of probabilistic metrics from `xskillscore` in
  `compute_perfect_model`, `bootstrap_perfect_model`, `compute_hindcast` and
  `bootstrap_hindcast`, now requires xskillscore>=0.05 (:pr:`232`) `Aaron Spring`_

Internals/Minor Fixes
---------------------
- Rename .stats.DPP to dpp (:pr:`232`) `Aaron Spring`_
- Add `matplotlib` as a main dependency so that a direct pip installation works
  (:pr:`211`) `Riley X. Brady`_.
- ``climpred`` is now installable from conda-forge (:pr:`212`) `Riley X. Brady`_.
- Fix erroneous descriptions of sample datasets (:pr:`226`) `Riley X. Brady`_.
- Benchmarking time and peak memory of compute functions with `asv` (:pr:`231`)
  `Aaron Spring`_

Documentation
-------------
- Add scope of package to docs for clarity for users and developers. (:pr:`235`)
  `Riley X. Brady`_.

climpred v1.0.1 (2019-07-04)
============================

Bug Fixes
---------
- Accomodate for lead-zero within the ``lead`` dimension (:pr:`196`) `Riley X. Brady`_.
- Fix issue with adding uninitialized ensemble to
  :py:class:`.HindcastEnsemble` object
  (:pr:`199`) `Riley X. Brady`_.
- Allow ``max_dof`` keyword to be passed to ``compute_metric`` and
  ``compute_persistence`` for :py:class:`.HindcastEnsemble`.
  (:pr:`199`) `Riley X. Brady`_.

Internals/Minor Fixes
---------------------
- Force ``xskillscore`` version 0.0.4 or higher to avoid ``ImportError``
  (:pr:`204`) `Riley X. Brady`_.
- Change ``max_dfs`` keyword to ``max_dof`` (:pr:`199`) `Riley X. Brady`_.
- Add tests for :py:class:`.HindcastEnsemble` and
  ``PerfectModelEnsemble``. (:pr:`199`) `Riley X. Brady`_

climpred v1.0.0 (2019-07-03)
============================
``climpred`` v1.0.0 represents the first stable release of the package. It includes
:py:class:`.HindcastEnsemble` and ``PerfectModelEnsemble`` objects to
perform analysis with.
It offers a suite of deterministic and probabilistic metrics that are optimized to be
run on single time series or grids of data (e.g., lat, lon, and depth). Currently,
``climpred`` only supports annual forecasts.

Features
--------
- Bootstrap prediction skill based on resampling with replacement consistently in
  ``ReferenceEnsemble`` and ``PerfectModelEnsemble``. (:pr:`128`) `Aaron Spring`_
- Consistent bootstrap function for ``climpred.stats`` functions via ``bootstrap_func``
  wrapper. (:pr:`167`) `Aaron Spring`_
- many more metrics: ``_msss_murphy``, ``_less`` and probabilistic ``_crps``,
  ``_crpss`` (:pr:`128`) `Aaron Spring`_

Bug Fixes
---------
- ``compute_uninitialized`` now trims input data to the same time window.
  (:pr:`193`) `Riley X. Brady`_
- ``rm_poly`` now properly interpolates/fills NaNs. (:pr:`192`) `Riley X. Brady`_

Internals/Minor Fixes
---------------------
- The ``climpred`` version can be printed. (:pr:`195`) `Riley X. Brady`_
- Constants are made elegant and pushed to a separate module. (:pr:`184`)
  `Andrew Huang`_
- Checks are consolidated to their own module. (:pr:`173`) `Andrew Huang`_

Documentation
-------------
- Documentation built extensively in multiple PRs.


climpred v0.3 (2019-04-27)
==========================

``climpred`` v0.3 really represents the entire development phase leading up to the
version 1 release. This was done in collaboration between `Riley X. Brady`_,
`Aaron Spring`_, and `Andrew Huang`_. Future releases will have less additions.

Features
--------
- Introduces object-oriented system to ``climpred``, with classes
  ``ReferenceEnsemble`` and ``PerfectModelEnsemble``. (:pr:`86`) `Riley X. Brady`_
- Expands bootstrapping module for perfect-module configurations. (:pr:`78`, :pr:`87`)
  `Aaron Spring`_
- Adds functions for computing Relative Entropy (:pr:`73`) `Aaron Spring`_
- Sets more intelligible dimension expectations for ``climpred``
  (:pr:`98`, :pr:`105`) `Riley X. Brady`_ and `Aaron Spring`_:

    -   ``init``:  initialization dates for the prediction ensemble
    -   ``lead``:  retrospective forecasts from prediction ensemble;
        returned dimension for prediction calculations
    -   ``time``:  time dimension for control runs, references, etc.
    -   ``member``:  ensemble member dimension.
- Updates ``open_dataset`` to display available dataset names when no argument is
  passed. (:pr:`123`) `Riley X. Brady`_
- Change ``ReferenceEnsemble`` to :py:class:`.HindcastEnsemble`.
  (:pr:`124`) `Riley X. Brady`_
- Add probabilistic metrics to ``climpred``. (:pr:`128`) `Aaron Spring`_
- Consolidate separate perfect-model and hindcast functions into singular functions
  (:pr:`128`) `Aaron Spring`_
- Add option to pass proxy through to ``open_dataset`` for firewalled networks.
  (:pr:`138`) `Riley X. Brady`_

Bug Fixes
---------
- ``xr_rm_poly`` can now operate on Datasets and with multiple variables.
  It also interpolates across NaNs in time series. (:pr:`94`) `Andrew Huang`_
- Travis CI, ``treon``, and ``pytest`` all run for automated testing of new features.
  (:pr:`98`, :pr:`105`, :pr:`106`) `Riley X. Brady`_ and `Aaron Spring`_
- Clean up ``check_xarray`` decorators and make sure that they work. (:pr:`142`)
  `Andrew Huang`_
- Ensures that ``help()`` returns proper docstring even with decorators.
  (:pr:`149`) `Andrew Huang`_
- Fixes bootstrap so p values are correct. (:pr:`170`) `Aaron Spring`_

Internals/Minor Fixes
---------------------
- Adds unit testing for all perfect-model comparisons. (:pr:`107`) `Aaron Spring`_
- Updates CESM-LE uninitialized ensemble sample data to have 34 members.
  (:pr:`113`) `Riley X. Brady`_
- Adds MPI-ESM hindcast, historical, and assimilation sample data.
  (:pr:`119`) `Aaron Spring`_
- Replaces ``check_xarray`` with a decorator for checking that input arguments are
  xarray objects. (:pr:`120`) `Andrew Huang`_
- Add custom exceptions for clearer error reporting. (:pr:`139`) `Riley X. Brady`_
- Remove "xr" prefix from stats module. (:pr:`144`) `Riley X. Brady`_
- Add codecoverage for testing. (:pr:`152`) `Riley X. Brady`_
- Update exception messages for more pretty error reporting. (:pr:`156`) `Andrew Huang`_
- Add ``pre-commit`` and ``flake8``/``black`` check in CI. (:pr:`163`) `Riley X. Brady`_
- Change ``loadutils`` module to ``tutorial`` and ``open_dataset`` to
  ``load_dataset``. (:pr:`164`) `Riley X. Brady`_
- Remove predictability horizon function to revisit for v2. (:pr:`165`)
  `Riley X. Brady`_
- Increase code coverage through more testing. (:pr:`167`) `Aaron Spring`_
- Consolidates checks and constants into modules. (:pr:`173`) `Andrew Huang`_

climpred v0.2 (2019-01-11)
==========================

Name changed to ``climpred``, developed enough for basic decadal prediction tasks on a
perfect-model ensemble and reference-based ensemble.

climpred v0.1 (2018-12-20)
==========================

Collaboration between Riley Brady and Aaron Spring begins.

.. _`Anderson Banihirwe`: https://github.com/andersy005
.. _`Ray Bell`: https://github.com/raybellwaves
.. _`Riley X. Brady`: https://github.com/bradyrx
.. _`Andrew Huang`: https://github.com/ahuang11
.. _`Kathy Pegion`: https://github.com/kpegion
.. _`Aaron Spring`: https://github.com/aaronspring
.. _`Dougie Squire`: https://github.com/dougiesquire
Release Procedure
-----------------

We follow semantic versioning, e.g., ``v1.0.0``. A major version causes incompatible API
changes, a minor version adds functionality, and a patch covers bug fixes.

#. Create a new branch ``release-v1.0.0`` with the version for the release.

 * Update `CHANGELOG.rst <CHANGELOG.html>`_.
 * Make sure all new changes and features are reflected in the documentation.

#. Open a new pull request for this branch targeting ``main``

#. After all tests pass and the PR has been approved, merge the PR into ``main``

#. Tag a release and push to github::

    $ git tag -a v1.0.0 -m "Version 1.0.0"
    $ git push upstream main --tags

#. We use Github Actions to automate the new release being published to PyPI.
   Simply confirm that the new release is reflected at
   https://pypi.org/project/climpred/. There is typically a delay, but check Github
   Actions if there's an issue, and otherwise you can manually do it with the
   following::

    $ git clean -xfd  # remove any files not checked into git
    $ python setup.py sdist bdist_wheel --universal  # build package
    $ twine upload dist/*  # register and push to pypi

#. Next, update the stable branch with ``main``. This will trigger a stable build
   for ReadTheDocs::

    $ git checkout stable
    $ git rebase main
    $ git push -f upstream stable
    $ git checkout main

#. Go to https://readthedocs.org and add the new version to ``"Active Versions"``
   under the version tab. Force-build ``"stable"`` if it isn't already building.

#. Update climpred conda-forge feedstock

 * Fork `climpred-feedstock repository <https://github.com/conda-forge/climpred-feedstock>`_
 * Clone this fork and edit recipe::

        $ git clone git@github.com:username/climpred-feedstock.git
        $ cd climpred-feedstock
        $ cd recipe
        $ # edit meta.yaml

 * Update version
 * Get ``sha256`` from pypi.org for `climpred <https://pypi.org/project/climpred/#files>`_
 * Check that ``requirements.txt`` from the main ``climpred`` repo is accounted for
   in ``meta.yaml`` from the feedstock.
 * Fill in the rest of information as described
   `here <https://github.com/conda-forge/climpred-feedstock#updating-climpred-feedstock>`_
 * Commit and submit a PR
.. image:: https://i.imgur.com/HPOdOsR.png

Verification of weather and climate forecasts.

..
    Table version of badges inspired by pySTEPS.

.. list-table::
    :stub-columns: 1
    :widths: 10 90

    * - docs
      - |docs| |joss| |doi|
    * - tests
      - |ci| |upstream| |codecov| |precommit|
    * - package
      - |conda| |conda downloads| |pypi| |pypi downloads|
    * - license
      - |license|
    * - community
      - |gitter| |contributors| |forks| |stars| |issues| |PRs|
    * - tutorials
      - |gallery| |workshop| |cloud|

.. |docs| image:: https://img.shields.io/readthedocs/climpred/stable.svg?style=flat
    :target: https://climpred.readthedocs.io/en/stable/?badge=stable
    :alt: Documentation Status

.. |joss| image:: https://joss.theoj.org/papers/246d440e3fcb19025a3b0e56e1af54ef/status.svg
    :target: https://joss.theoj.org/papers/246d440e3fcb19025a3b0e56e1af54ef
    :alt: JOSS paper

.. |doi| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4556085.svg
    :target: https://doi.org/10.5281/zenodo.4556085
    :alt: DOI

.. |ci| image:: https://github.com/pangeo-data/climpred/workflows/climpred%20testing/badge.svg
    :target: https://github.com/pangeo-data/climpred/actions/workflows/climpred_testing.yml
    :alt: CI

.. |upstream| image:: https://github.com/pangeo-data/climpred/actions/workflows/upstream-dev-ci.yml/badge.svg
    :target: https://github.com/pangeo-data/climpred/actions/workflows/upstream-dev-ci.yml
    :alt: CI upstream

.. |codecov| image:: https://codecov.io/gh/pangeo-data/climpred/branch/main/graph/badge.svg
      :target: https://codecov.io/gh/pangeo-data/climpred
      :alt: coverage

.. |precommit| image:: https://results.pre-commit.ci/badge/github/pangeo-data/climpred/main.svg
   :target: https://results.pre-commit.ci/latest/github/pangeo-data/climpred/main
   :alt: pre-commit.ci status

.. |conda| image:: https://img.shields.io/conda/vn/conda-forge/climpred.svg
    :target: https://anaconda.org/conda-forge/climpred
    :alt: Conda Version

.. |pypi| image:: https://img.shields.io/pypi/v/climpred.svg
   :target: https://pypi.python.org/pypi/climpred/
   :alt: pypi Version

.. |license| image:: https://img.shields.io/github/license/pangeo-data/climpred.svg
    :alt: license
    :target: LICENSE.txt

.. |gitter| image:: https://badges.gitter.im/Join%20Chat.svg
    :target: https://gitter.im/climpred
    :alt: gitter chat

.. |contributors| image:: https://img.shields.io/github/contributors/pangeo-data/climpred
    :alt: GitHub contributors
    :target: https://github.com/pangeo-data/climpred/graphs/contributors

.. |conda downloads| image:: https://img.shields.io/conda/dn/conda-forge/climpred
    :alt: Conda downloads
    :target: https://anaconda.org/conda-forge/climpred

.. |pypi downloads| image:: https://pepy.tech/badge/climpred
    :alt: pypi downloads
    :target: https://pepy.tech/project/climpred

.. |gallery| image:: https://img.shields.io/badge/climpred-examples-ed7b0e.svg
    :alt: climpred gallery
    :target: https://mybinder.org/v2/gh/pangeo-data/climpred/main?urlpath=lab%2Ftree%2Fdocs%2Fsource%2Fquick-start.ipynb

.. |workshop| image:: https://img.shields.io/badge/climpred-workshop-f5a252
    :alt: climpred workshop
    :target: https://mybinder.org/v2/gh/bradyrx/climpred_workshop/master

.. |cloud| image:: https://img.shields.io/badge/climpred-cloud_demo-f9c99a
    :alt: climpred cloud demo
    :target: https://github.com/aaronspring/climpred-cloud-demo

.. |forks| image:: https://img.shields.io/github/forks/pangeo-data/climpred
    :alt: GitHub forks
    :target: https://github.com/pangeo-data/climpred/network/members

.. |stars| image:: https://img.shields.io/github/stars/pangeo-data/climpred
    :alt: GitHub stars
    :target: https://github.com/pangeo-data/climpred/stargazers

.. |issues| image:: https://img.shields.io/github/issues/pangeo-data/climpred
    :alt: GitHub issues
    :target: https://github.com/pangeo-data/climpred/issues

.. |PRs| image:: https://img.shields.io/github/issues-pr/pangeo-data/climpred
    :alt: GitHub PRs
    :target: https://github.com/pangeo-data/climpred/pulls

..


.. note::

    We are actively looking for new contributors for climpred! Riley moved to McKinsey's
    Climate Analytics team. Aaron is finishing his PhD, but will stay in academia.
    We especially hope for python enthusiasts from seasonal, subseasonal or weather
    prediction community. In our past coding journey, collaborative coding, feedbacking
    issues and pull requests advanced our code and thinking about forecast verification
    more than we could have ever expected.
    `Aaron <https://github.com/aaronspring/>`_ can provide guidance on
    implementing new features into ``climpred``. Feel free to implement
    your own new feature or take a look at the
    `good first issue <https://github.com/pangeo-data/climpred/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22>`_
    tag in the issues. Please reach out to us via `gitter <https://gitter.im/climpred>`_.


Installation
============

You can install the latest release of ``climpred`` using ``pip`` or ``conda``:

.. code-block:: bash

    pip install climpred[complete]

.. code-block:: bash

    conda install -c conda-forge climpred

You can also install the bleeding edge (pre-release versions) by cloning this
repository or installing directly from GitHub:

.. code-block:: bash

    git clone https://github.com/pangeo-data/climpred.git
    cd climpred
    pip install . --upgrade

.. code-block:: bash

    pip install git+https://github.com/pangeo-data/climpred.git


Documentation
=============

Documentation is in development and can be found on readthedocs_.

.. _readthedocs: https://climpred.readthedocs.io/en/latest/
####################
Significance Testing
####################

Significance testing is important for assessing whether a given initialized prediction
system is skillful. Some questions that significance testing can answer are:

    - Is the correlation coefficient of a lead time series significantly different from
      zero?

    - What is the probability that the retrospective forecast is more valuable than a
      historical/uninitialized simulation?

    - Are correlation coefficients statistically significant despite temporal and
      spatial autocorrelation?

All of these questions deal with statistical significance. See below on how to use
``climpred`` to address these questions.
Please also have a look at the
`significance testing example <examples/decadal/significance.html>`__.

p value for temporal correlations
#################################

For the correlation `metrics <metrics.html>`__, like
:py:func:`~climpred.metrics._pearson_r` and :py:func:`~climpred.metrics._spearman_r`,
``climpred`` also hosts the associated p-value, like
:py:func:`~climpred.metrics._pearson_r_p_value`,
that this correlation is significantly different from zero.
:py:func:`~climpred.metrics._pearson_r_eff_p_value` also incorporates the reduced
degrees of freedom due to temporal autocorrelation. See
`example <examples/decadal/significance.html#p-value-for-temporal-correlations>`__.

Bootstrapping with replacement
##############################

Testing statistical significance through bootstrapping is commonly used in the field of
climate prediction. Bootstrapping relies on
resampling the underlying data with replacement for a large number of ``iterations``, as
proposed by the decadal prediction framework :cite:p:`Goddard2013,Boer2016`.
This means that the ``initialized`` ensemble is resampled with replacement along a
dimension (``init`` or ``member``) and then that resampled ensemble is verified against
the observations. This leads to a distribution of ``initialized`` skill. Further, a
``reference`` forecast uses the resampled ``initialized`` ensemble, which creates a
``reference`` skill distribution. Lastly, an ``uninitialized`` skill distribution is
created from the underlying historical members or the control simulation.

The probability or p value is the fraction of these resampled ``initialized`` metrics
beaten by the ``uninitialized`` or resampled reference metrics calculated from their
respective distributions. Confidence intervals using these distributions are also
calculated.

This behavior is incorporated by :py:meth:`.HindcastEnsemble.bootstrap` and
:py:meth:`.PerfectModelEnsemble.bootstrap`, see
`example <examples/decadal/significance.html#Bootstrapping-with-replacement>`__.


Field significance
##################

Please use :py:func:`esmtools.testing.multipletests` to control the false discovery
rate (FDR) in geospatial data from the above obtained p-values :cite:p:`Wilks2016`.
See the `FDR example <examples/decadal/significance.html#Field-significance>`__.


Sign test
#########

Use DelSole's sign test relying on the statistics of a random walk to decide whether
one forecast is significantly better than another forecast
:cite:p:`Benjamini1994,DelSole2016`, see :py:func:`xskillscore.sign_test` and
`sign test example <examples/decadal/significance.html#sign-test>`__.

References
##########

.. bibliography::
  :filter: docname in docnames
**********************
Prediction Terminology
**********************

Terminology is often confusing and highly variable amongst those that make predictions
in the geoscience community. Here we define some common terms in climate prediction and
how we use them in ``climpred``.

Simulation Design
#################

*Hindcast Ensemble* (:py:class:`.HindcastEnsemble`):
Ensemble members are initialized from a simulation (generally a reconstruction from
reanalysis) or an analysis (representing the current state of the atmosphere, land, and
ocean by assimilation of observations) at initialization dates and integrated for some
lead years :cite:p:`Boer2016`.

*Perfect Model Experiment* (:py:class:`.PerfectModelEnsemble`):
Ensemble members are initialized from a control simulation
(:py:meth:`.PerfectModelEnsemble.add_control`) at randomly chosen
initialization dates and integrated for some lead years :cite:p:`Griffies1997`.

*Reconstruction/Assimilation*: (:py:meth:`.HindcastEnsemble.add_observations`)
A "reconstruction" is a model solution that uses
observations in some capacity to approximate historical or current conditions of the
atmosphere, ocean, sea ice, and/or land. This could be done via a forced simulation,
such as an OMIP run that uses a dynamical ocean/sea ice core with reanalysis forcing
from atmospheric winds. This could also be a fully data assimilative model, which
assimilates observations into the model solution.  For weather, subseasonal, and
seasonal predictions, the terms re-analysis and analysis are the terms typically used,
while reconstruction is more commonly used for decadal predictions.

*Uninitialized Ensemble*: (:py:meth:`.HindcastEnsemble.add_uninitialized`)
In this framework, an *uninitialized ensemble* is one that
is generated by perturbing initial conditions only at one point in the historical run.
These are generated via micro (round-off error perturbations) or macro (starting from
completely different restart files) methods. Uninitialized ensembles are used to
approximate the magnitude of internal climate variability and to confidently extract
the forced response (ensemble mean) in the climate system. In ``climpred``, we use
uninitialized ensembles as a baseline for how important (reoccurring) initializations
are for lending predictability to the system. Some modeling centers (such as NCAR)
provide a dynamical uninitialized ensemble (the CESM Large Ensemble) along with their
initialized prediction system (the CESM Decadal Prediction Large Ensemble). If this
isn't available, one can approximate the unintiailized response by bootstrapping a
control simulation.

Forecast Assessment
###################

*Accuracy*: The average degree of correspondence between individual pairs of forecasts
and observations :cite:p:`Murphy1988,Jolliffe2011`. Examples include Mean Absolute Error
(MAE) :py:func:`~.climpred.metrics._mae` and Mean Square Error (MSE)
:py:func:`~.climpred.metrics._mse`. See `metrics <metrics.html>`_.

*Association*: The overall strength of the relationship between individual pairs of
forecasts and observations :cite:p:`Jolliffe2011`. The primary measure of association
is the Anomaly Correlation Coefficient (ACC), which can be measured using the Pearson
product-moment correlation :py:func:`~.climpred.metrics._pearson_r` or
Spearman's Rank correlation :py:func:`~.climpred.metrics._spearman_r`. See
`metrics <metrics.html>`_.

*(Potential) Predictability*: This characterizes the "ability to be predicted"
rather than the current "capability to predict." One estimates this by computing a
metric (like the anomaly correlation coefficient (ACC)) between the prediction
ensemble and a member (or collection of members) selected as the verification member(s)
(in a perfect-model setup) or the reconstruction that initialized it
(in a hindcast setup) :cite:p:`Meehl2013,Pegion2019`.

*(Prediction) Skill*: (:py:meth:`.HindcastEnsemble.verify`)
This characterizes the current ability of the ensemble
forecasting system to predict the real world. This is derived by computing a metric
between the prediction ensemble and observations, reanalysis, or analysis of the real
world :cite:p:`Meehl2013,Pegion2019`.

*Skill Score*: The most generic skill score can be defined as the following
:cite:t:`Murphy1988`:

.. math::
    S = \frac{A_{f} - A_{r}}{A_{p} - A_{r}},

where :math:`A_{f}`, :math:`A_{p}`, and :math:`A_{r}` represent the accuracy of the
forecast being assessed, the accuracy of a perfect forecast, and the accuracy of the
reference forecast (e.g. persistence), respectively :cite:p:`Murphy1985`. Here,
:math:`S` represents the improvement in accuracy of the forecasts over the reference
forecasts relative to the total possible improvement in accuracy. They are typically
designed to take a value of 1 for a perfect forecast and 0 for equivalent to the
reference forecast :cite:p:`Jolliffe2011`.

Forecasting
###########

*Hindcast*: Retrospective forecasts of the past initialized from a reconstruction
integrated forward in time, also called re-forcasts.  Depending on the length of time
of the integration, external forcings may or may not be included.  The longer the
integration (e.g. decadal vs. daily), the more important it is to include external
forcing :cite:p:`Boer2016`.  Because they represent so-called forecasts over periods
that already occurred, their prediction skill can be evaluated.

*Prediction*: Forecasts initialized from a reconstruction integrated into the future.
Depending on the length of time of the integration, external forcings may or may not
be included.  The longer the integration (e.g. decadal vs. daily), the more important
it is to include external forcing :cite:p:`Boer2016`. Because predictions are made into
the future, it is necessary to wait until the forecast occurs before one can quantify
the skill of the forecast.

*Projection* An estimate of the future climate that is dependent on the externally
forced climate response, such as anthropogenic greenhouse gases, aerosols, and
volcanic eruptions :cite:p:`Meehl2013`.

References
##########

.. bibliography::
  :filter: docname in docnames
=====================
Contribution Guide
=====================

Contributions are highly welcomed and appreciated.  Every little help counts,
so do not hesitate! You can make a high impact on ``climpred`` just by using
it and reporting `issues <https://github.com/pangeo-data/climpred/issues>`__.

The following sections cover some general guidelines
regarding development in ``climpred`` for maintainers and contributors.

Please also review our `Code of Conduct <code_of_conduct.html>`__.

Nothing here is set in stone and can't be changed.
Feel free to suggest improvements or changes in the workflow.


.. _submitfeedback:

Feature requests and feedback
-----------------------------

We are eager to hear about your requests for new features and any suggestions
about the API, infrastructure, and so on. Feel free to submit these as
`issues <https://github.com/pangeo-data/climpred/issues/new>`__ with the label
``"feature request"``.

Please make sure to explain in detail how the feature should work and keep the
scope as narrow as possible. This will make it easier to implement in small
PRs.


.. _reportbugs:

Report bugs
-----------

Report bugs for ``climpred`` in the
`issue tracker <https://github.com/pangeo-data/climpred/issues>`_ with the
label "bug".

If you are reporting a bug, please include:

* Any details about your local setup that might be helpful in troubleshooting,
  specifically the Python interpreter version, installed libraries, and
  ``climpred`` version.
* Detailed steps to reproduce the bug.

If you can write a demonstration test that currently fails but should pass,
that is a very useful commit to make as well, even if you cannot fix the bug
itself.


.. _fixbugs:

Bug Fix
-------

Look through the
`GitHub issues for bugs <https://github.com/pangeo-data/climpred/labels/bug>`_.

Talk to developers to find out how you can fix specific bugs.


Write documentation
-------------------

``climpred`` could always use more documentation.  What exactly is needed?

* More complementary documentation.  Have you perhaps found something unclear?
* Example notebooks with different Earth System Models, lead times, etc. --
  they're all very appreciated.

You can also edit documentation files directly in the GitHub web interface,
without using a local copy.  This can be convenient for small fixes.

Our documentation is written in reStructuredText. You can follow our
conventions in already written documents. Some helpful guides are located
`rst-quickref <http://docutils.sourceforge.net/docs/user/rst/quickref.html>`__
and
`rst-cheatsheet <https://github.com/ralsina/rst-cheatsheet/blob/master/rst-cheatsheet.rst>`__.

.. note::
    Build the documentation locally with the following command:

    .. code:: bash

        $ conda env update -f ci/requirements/climpred-dev.yml
        $ cd docs
        $ make html

    The built documentation should be available in the ``docs/build/``.

If you need to add new functions to the API, run
``sphinx-autogen -o api api.rst`` from the ``docs/source`` directory after
adding functions to ``api.rst``.

 .. _`pull requests`:
 .. _pull-requests:

Preparing Pull Requests
-----------------------

#. Fork the `climpred GitHub repository <https://github.com/pangeo-data/climpred>`__.
   It's fine to use ``climpred`` as your fork repository name because it will
   live under your user.

#. Clone your fork locally using `git <https://git-scm.com/>`_, connect your
   repository to the upstream (main project), and create a branch::

    $ git clone git@github.com:YOUR_GITHUB_USERNAME/climpred.git
    $ cd climpred
    $ git remote add upstream git@github.com:pangeo-data/climpred.git

    # now, to fix a bug or add feature create your own branch off "main":

    $ git checkout -b your-bugfix-feature-branch-name main

   If you need some help with Git, follow this quick start
   `guide <https://git.wiki.kernel.org/index.php/QuickStart>`_.

#. Install dependencies into a new
   `conda <https://conda.io/projects/conda/en/latest/user-guide/getting-started.html>`_
   environment::

    $ conda env create -f ci/requirements/climpred-dev.yml
    $ conda activate climpred-dev

#. Make an editable install of ``climpred`` by running::

    $ pip install -e .

#. Install `pre-commit <https://pre-commit.com>`_ and its hook on the
   ``climpred`` repo::

     $ pip install --user pre-commit
     $ pre-commit install

   ``pre-commit`` automatically beautifies the code, makes it more
   maintainable and catches syntax errors. Afterwards ``pre-commit`` will run
   whenever you commit.

   Now you have an environment called ``climpred-dev`` that you can work in.
   You’ll need to make sure to activate that environment next time you want
   to use it after closing the terminal or your system.

   You can now edit your local working copy and run/add tests as necessary.
   Please try to follow
   `PEP-8 <https://www.python.org/dev/peps/pep-0008/#naming-conventions>`_ for
   naming. When committing, ``pre-commit`` will modify the files as
   needed, or will generally be quite clear about what you need to do to pass
   the commit test.

   ``pre-commit`` also runs::

    * `mypy <http://mypy-lang.org/>`_ for static type checking on
      `type hints <https://docs.python.org/3/library/typing.html>`_.
    * `isort <https://pycqa.github.io/isort/>`_ sorting imports
    * `black <https://black.readthedocs.io/en/stable/>`_ code formatting
    * `flake8 <https://flake8.pycqa.org/en/latest/>`_ code linting
    * `blackdoc <https://blackdoc.readthedocs.io/en/latest/>`_ docstring code
      formatter


#. Break your edits up into reasonably sized commits::

    $ git commit -a -m "<commit message>"
    $ git push -u

#. Run all tests

   Once commits are pushed to ``origin``, GitHub Actions runs continuous
   integration of all tests on all new commits. However, you are already
   run tests locally::

    $ pytest climpred

   Check that `doctests <https://docs.pytest.org/en/stable/doctest.html>`_ are
   passing::

    $ pytest --doctest-modules climpred --ignore climpred/tests

   Check that your contribution is covered by tests and therefore increases
   the overall test coverage::

    $ coverage run --source climpred -m py.test
    $ coverage report
    $ coveralls

   Please stick to
   `xarray <http://xarray.pydata.org/en/stable/contributing.html>`_'s testing
   recommendations.

#. Running the performance test suite

   If you considerably changed to core of code of ``climpred``, it is worth
   considering whether your code has introduced performance regressions.
   ``climpred`` has a suite of benchmarking tests using
   `asv <https://asv.readthedocs.io/en/stable/>`_
   to enable easy monitoring of the performance of critical ``climpred``
   operations. These benchmarks are all found in the ``asv_bench`` directory.

   If you need to run a benchmark, change your directory to ``asv_bench/`` and
   run::

      $ asv continuous -f 1.1 upstream/main HEAD

   You can replace ``HEAD`` with the name of the branch you are working on,
   and report benchmarks that changed by more than 10%.
   The command uses ``conda`` by default for creating the benchmark
   environments.

   Running the full benchmark suite can take up to half an hour and use up a
   few GBs of RAM. Usually it is sufficient to paste only a subset of the
   results into the pull request to show that the committed changes do not
   cause unexpected performance regressions.
   If you want to only run a specific group of tests from a file, you can do it
   using ``.`` as a separator. For example::

      $ asv continuous -f 1.1 upstream/main HEAD -b benchmarks_PredictionEnsemble.GenerateHindcastEnsembleSmall.time_bootstrap

   will only run the ``time_bootstrap`` benchmark of class
   ``GenerateHindcastEnsembleSmall`` defined in ``benchmarks_PredictionEnsemble.py``.


#. Create a new changelog entry in `CHANGELOG.rst <CHANGELOG.html>`_:

   The entry should be entered as:

   ``<description>`` (``:pr:`#<pull request number>```) ```<author's names>`_``

   where ``<description>`` is the description of the PR related to the change
   and ``<pull request number>`` is the pull request number and
   ``<author's names>`` are your first and last names.

   Add yourself to list of authors at the end of `CHANGELOG.rst <CHANGELOG.html>`_ file if
   not there yet, in alphabetical order.

#. Add yourself to the `contributors <https://climpred.readthedocs.io/en/latest/contributors.html>`_ list via ``docs/source/contributors.rst``.

#. Finally, submit a `Pull Request <https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests>`_ through the GitHub website using this data::

    head-fork: YOUR_GITHUB_USERNAME/climpred
    compare: your-branch-name

    base-fork: pangeo-data/climpred
    base: main

Note that you can create the ``Pull Request`` while you're working on this.
The PR will update as you add more commits. ``climpred`` developers and
contributors can then review your code and offer suggestions.
************
Contributors
************

.. note::
  We are actively looking for new contributors for climpred! Riley moved to McKinsey's
  Climate Analytics team. Aaron is finishing his PhD, but will stay in academia.
  We especially hope for python enthusiasts from seasonal, subseasonal or weather
  prediction community. In our past coding journey, collaborative coding, feedbacking
  issues and pull requests advanced our code and thinking about forecast verification
  more than we could have ever expected.
  `Aaron <https://github.com/aaronspring/>`_ can provide guidance on
  implementing new features into climpred. Feel free to implement
  your own new feature or take a look at the
  `good first issue <https://github.com/pangeo-data/climpred/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22>`_
  tag in the issues. Please reach out to us via `gitter <https://gitter.im/climpred>`_.

Core Developers
===============
* Aaron Spring (`github <https://github.com/aaronspring/>`__)
* Riley X. Brady (`github <https://github.com/bradyrx/>`__)

Contributors
============
* Andrew Huang (`github <https://github.com/ahuang11/>`__)
* Kathy Pegion (`github <https://github.com/kpegion/>`__)
* Anderson Banihirwe (`github <https://github.com/andersy005/>`__)
* Ray Bell (`github <https://github.com/raybellwaves/>`__)

For a list of all the contributions, see the github
`contribution graph <https://github.com/pangeo-data/climpred/graphs/contributors>`_.
Examples
########

.. note::

    Please use the ``climpred-dev`` environment

    .. code-block:: bash

       conda env create -f ci/requirements/climpred-dev.yml

    to ensure that all dependencies are installed to complete all example
    notebooks listed here.


Numerical Weather Prediction
============================
.. toctree::
  :maxdepth: 1

  examples/NWP/NWP_GEFS_6h_forecasts.ipynb


Subseasonal
===========
.. toctree::
  :maxdepth: 1

  examples/subseasonal/daily-subx-example.ipynb
  examples/subseasonal/daily-S2S-IRIDL.ipynb
  examples/subseasonal/weekly-subx-example.ipynb
  examples/subseasonal/daily-S2S-ECMWF.ipynb


Monthly and Seasonal
====================
.. toctree::
  :maxdepth: 1

  examples/monseas/monthly-enso-subx-example.ipynb
  examples/monseas/seasonal-enso-subx-example.ipynb


Decadal
=======
.. toctree::
  :maxdepth: 1

  examples/decadal/perfect-model-predictability-demo.ipynb
  examples/decadal/tropical-pacific-ssts.ipynb
  examples/decadal/diagnose-potential-predictability.ipynb
  examples/decadal/Significance.ipynb


Misc
====
.. toctree::
  :maxdepth: 1

  examples/misc/efficient_dask.ipynb
  examples/misc/climpred_gpu.ipynb
  examples/misc/setup_your_own_data.ipynb
*************
Helpful Links
*************

We hope to curate in the ``climpred`` documentation a comprehensive report of
terminology, best practices, analysis methods, etc. in the prediction community.
Here we suggest other resources for initialized prediction of the Earth system to round
out the information provided in our documentation.

Forecast Verification
#####################

* `CAWCR Forecast Verification Overview <https://www.cawcr.gov.au/projects/verification/>`_:
  A nice overview of forecast verification, including a suite of metrics and their
  derivation.
********************
Initialized Datasets
********************

Probably the hardest part in working with ``climpred`` is getting the ``initialized``
dataset complying to the expectations and data model of ``climpred``.
For names, data types and conventions of :py:class:`xarray.Dataset` dimensions and
coordinates, please refer to `Setting up your Dataset <setting-up-data.html>`_.

Here, we list publicly available initialized datasets and corresponding ``climpred``
examples:

.. list-table:: List of initialized Datasets
   :widths: 25 15 40 40 25 25
   :header-rows: 1

   * - Short Name
     - Community
     - Description
     - Data Source
     - Reference Paper
     - Example
   * - DCPP
     - decadal
     - `Decadal Climate Prediction Project (DCPP) contribution to CMIP6 <https://www.wcrp-climate.org/dcp-overview>`_
     - `ESGF <https://esgf-data.dkrz.de/search/cmip6-dkrz/>`_, `pangeo <https://pangeo-data.github.io/pangeo-cmip6-cloud/accessing_data.html#loading-an-esm-collection>`_
     - :cite:t:`Boer2016`
     - `with intake-esm <examples/misc/setup_your_own_data.html#intake-esm-for-cmorized-output>`_, `Anderson <https://github.com/andersy005>`_ at NOAA's 45th CDP Workshop: `slides <https://talks.andersonbanihirwe.dev/climpred-cdpw-2020.html>`_, `Notebook <https://nbviewer.jupyter.org/github/andersy005/talks/blob/gh-pages/notebooks/climpred-demo.ipynb>`_
   * - CESM-DPLE
     - decadal
     - `Decadal Prediction Large Ensemble Project <http://www.cesm.ucar.edu/projects/community-projects/DPLE/>`_
     - `Data <https://www.earthsystemgrid.org/dataset/ucar.cgd.ccsm4.CESM1-CAM5-DP.html>`_
     - :cite:t:`Yeager2018`
     - many standard climpred `examples <quick-start.html>`_
   * - NMME
     - seasonal
     - `The North American Multimodel Ensemble: Phase-1 Seasonal-to-Interannual Prediction <https://www.cpc.ncep.noaa.gov/products/NMME/>`_
     - `IRIDL <http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/>`__
     - :cite:t:`Kirtman2014`
     - `seasonal SubX <examples.html#monthly-and-seasonal>`_
   * - SubX
     - subseasonal
     - `A Multimodel Subseasonal Prediction Experiment <http://cola.gmu.edu/subx/>`_
     - `IRIDL <http://iridl.ldeo.columbia.edu/SOURCES/.Models/.SubX/>`__
     - :cite:t:`Pegion2019`
     - `subseasonal SubX <examples.html#subseasonal>`_
   * - S2S
     - subseasonal
     - `The Subseasonal to Seasonal (S2S) Prediction Project Database <http://wwww.s2sprediction.net/>`_
     - `IRIDL <https://iridl.ldeo.columbia.edu/SOURCES/.ECMWF/.S2S/>`__, `climetlab <https://github.com/ecmwf-lab/climetlab-s2s-ai-challenge>`_
     - :cite:t:`Vitart2017`
     - `IRIDL <examples/subseasonal/daily-S2S-IRIDL.html>`_, `EWC Cloud/climetlab <examples/subseasonal/daily-S2S-ECMWF.html>`_
   * - GEFS
     - weather
     - `Global Ensemble Forecast System (GEFS) <https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-ensemble-forecast-system-gefs>`_
     - `NOAA THREDDS <https://www.ncei.noaa.gov/thredds/catalog/model-gefs-003/catalog.html>`_
     - :cite:t:`Toth1993`
     - `GEFS NWP <examples/NWP/NWP_GEFS_6h_forecasts.html>`_
   * - name
     - weather
     - please add a `Pull Request <contributing.html>`_ for numerical weather prediction
     - dataset
     - appreciated
     - `examples to add <https://github.com/pangeo-data/climpred/issues/602>`_

If you find or use another publicly available initialized datasets, please consider
adding a `Pull Request <contributing.html>`_.

References
##########

.. bibliography::
  :filter: docname in docnames
Scope of ``climpred``
=====================

``climpred`` aims to be the primary package used to analyze output from initialized
dynamical forecast models, ranging from short-term weather forecasts to decadal climate
forecasts. The code base is driven by the geoscientific prediction community through
open source development. It leverages `xarray <http://xarray.pydata.org/en/stable/>`_
to keep track of core prediction ensemble dimensions (e.g., ensemble member,
initialization date, and lead time) and `dask <https://dask.org/>`_ to perform
out-of-memory computations on large datasets.

The primary goal of ``climpred`` is to offer a comprehensive set of analysis tools for
assessing the forecasts relative to a validation product (e.g., observations,
reanalysis products, control simulations, baseline forecasts). This ranges from simple
deterministic and probabilistic verification `metrics <metrics.html>`_ — such as, e.g.
mean absolute error or rank histogram — to more advanced contingency table-derived
metrics. ``climpred`` expects users to handle their domain-specific post-processing of
model output, so that the package can focus on the actual analysis of forecasts.

Finally, the ``climpred`` documentation will serve as a repository of unified analysis
methods through `jupyter <https://jupyter.org/>`_ notebook `examples <examples.html>`_,
and collects relevant references and literature.
.. include:: ../../HOWTORELEASE.rst
**********
Literature
**********

References used in the documentation and used for studying predictability.

.. bibliography::
  :all:
.. currentmodule:: climpred.metrics

.. ipython:: python
    :suppress:

    from climpred.metrics import __ALL_METRICS__ as all_metrics

    metric_aliases = {}
    for m in all_metrics:
        if m.aliases is not None:
            metric_list = [m.name] + m.aliases
        else:
            metric_list = [m.name]
        metric_aliases[m.name] = metric_list

#######
Metrics
#######

All high-level functions like :py:meth:`.HindcastEnsemble.verify`,
:py:meth:`.HindcastEnsemble.bootstrap`, :py:meth:`.PerfectModelEnsemble.verify` and
:py:meth:`.PerfectModelEnsemble.bootstrap` have a ``metric`` argument
that has to be called to determine which metric is used in computing predictability.

.. note::

    We use the term 'observations' ``o`` here to refer to the 'truth' data to which
    we compare the forecast ``f``. These metrics can also be applied relative
    to a control simulation, reconstruction, observations, etc. This would just
    change the resulting score from quantifying skill to quantifying potential
    predictability.

Internally, all metric functions require ``forecast`` and ``observations`` as inputs.
The dimension ``dim`` has to be set to specify over which dimensions
the ``metric`` is applied and are hence reduced.
See :ref:`comparisons` for more on the ``dim`` argument.

*************
Deterministic
*************

Deterministic metrics assess the forecast as a definite prediction of the future, rather
than in terms of probabilities. Another way to look at deterministic metrics is that
they are a special case of probabilistic metrics where a value of one is assigned to
one category and zero to all others :cite:p:`Jolliffe2011`.

Correlation Metrics
===================

The below metrics rely fundamentally on correlations in their computation. In the
literature, correlation metrics are typically referred to as the Anomaly Correlation
Coefficient (ACC). This implies that anomalies in the forecast and observations
are being correlated. Typically, this is computed using the linear
`Pearson Product-Moment Correlation <#pearson-product-moment-correlation-coefficient>`_.
However, ``climpred`` also offers the
`Spearman's Rank Correlation <#spearman-s-rank-correlation-coefficient>`_.

Note that the p value associated with these correlations is computed via a separate
metric. Use :py:func:`~climpred.metrics._pearson_r_p_value` or
:py:func:`~climpred.metrics._spearman_r_p_value` to compute p values assuming
that all samples in the correlated time series are independent. Use
:py:func:`~climpred.metrics._pearson_r_eff_p_value` or
:py:func:`~climpred.metrics._spearman_r_eff_p_value` to account for autocorrelation
in the time series by calculating the
:py:func:`~climpred.metrics._effective_sample_size`.

Pearson Product-Moment Correlation Coefficient
----------------------------------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['pearson_r']}")

.. autofunction:: _pearson_r

Pearson Correlation p value
---------------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['pearson_r_p_value']}")

.. autofunction:: _pearson_r_p_value

Effective Sample Size
---------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['effective_sample_size']}")

.. autofunction:: _effective_sample_size

Pearson Correlation Effective p value
-------------------------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['pearson_r_eff_p_value']}")

.. autofunction:: _pearson_r_eff_p_value


Spearman's Rank Correlation Coefficient
---------------------------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['spearman_r']}")

.. autofunction:: _spearman_r


Spearman's Rank Correlation Coefficient p value
-----------------------------------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['spearman_r_p_value']}")

.. autofunction:: _spearman_r_p_value

Spearman's Rank Correlation Effective p value
---------------------------------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['spearman_r_eff_p_value']}")

.. autofunction:: _spearman_r_eff_p_value

Distance Metrics
================

This class of metrics simply measures the distance (or difference) between forecasted
values and observed values.

Mean Squared Error (MSE)
------------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['mse']}")

.. autofunction:: _mse


Root Mean Square Error (RMSE)
-----------------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['rmse']}")

.. autofunction:: _rmse


Mean Absolute Error (MAE)
-------------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['mae']}")

.. autofunction:: _mae


Median Absolute Error
---------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['median_absolute_error']}")

.. autofunction:: _median_absolute_error


Spread
------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['spread']}")

.. autofunction:: _spread


Multiplicative bias
-------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['mul_bias']}")

.. autofunction:: _mul_bias



Normalized Distance Metrics
===========================

Distance metrics like ``mse`` can be normalized to 1. The normalization factor
depends on the comparison type choosen. For example, the distance between an ensemble
member and the ensemble mean is half the distance of an ensemble member with other
ensemble members. See :py:func:`~climpred.metrics._get_norm_factor`.

Normalized Mean Square Error (NMSE)
-----------------------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['nmse']}")

.. autofunction:: _nmse


Normalized Mean Absolute Error (NMAE)
-------------------------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['nmae']}")

.. autofunction:: _nmae


Normalized Root Mean Square Error (NRMSE)
-----------------------------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['nrmse']}")

.. autofunction:: _nrmse


Mean Square Error Skill Score (MSESS)
-------------------------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['msess']}")

.. autofunction:: _msess


Mean Absolute Percentage Error (MAPE)
-------------------------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['mape']}")

.. autofunction:: _mape

Symmetric Mean Absolute Percentage Error (sMAPE)
------------------------------------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['smape']}")

.. autofunction:: _smape


Unbiased Anomaly Correlation Coefficient (uACC)
-----------------------------------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['uacc']}")

.. autofunction:: _uacc


Murphy Decomposition Metrics
============================

Metrics derived in :cite:p:`Murphy1988` which decompose the ``MSESS`` into a correlation term,
a conditional bias term, and an unconditional bias term. See
https://www-miklip.dkrz.de/about/murcss/ for a walk through of the decomposition.

Standard Ratio
--------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['std_ratio']}")

.. autofunction:: _std_ratio

Conditional Bias
----------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['conditional_bias']}")

.. autofunction:: _conditional_bias

Unconditional Bias
------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['unconditional_bias']}")

Simple bias of the forecast minus the observations.

.. autofunction:: _unconditional_bias

Bias Slope
----------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['bias_slope']}")

.. autofunction:: _bias_slope

Murphy's Mean Square Error Skill Score
--------------------------------------

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['msess_murphy']}")

.. autofunction:: _msess_murphy

*************
Probabilistic
*************

Probabilistic metrics include the spread of the ensemble simulations in their
calculations and assign a probability value between 0 and 1 to their forecasts
:cite:p:`Jolliffe2011`.

Continuous Ranked Probability Score (CRPS)
==========================================

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['crps']}")

.. autofunction:: _crps

Continuous Ranked Probability Skill Score (CRPSS)
=================================================

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['crpss']}")

.. autofunction:: _crpss

Continuous Ranked Probability Skill Score Ensemble Spread
=========================================================

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['crpss_es']}")

.. autofunction:: _crpss_es

Brier Score
===========

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['brier_score']}")

.. autofunction:: _brier_score

Threshold Brier Score
=====================

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['threshold_brier_score']}")

.. autofunction:: _threshold_brier_score

Ranked Probability Score
========================

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['rps']}")

.. autofunction:: _rps

Reliability
===========

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['reliability']}")

.. autofunction:: _reliability

Discrimination
==============

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['discrimination']}")

.. autofunction:: _discrimination

Rank Histogram
==============

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['rank_histogram']}")

.. autofunction:: _rank_histogram

Logarithmic Ensemble Spread Score
=================================

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['less']}")

.. autofunction:: _less

*************************
Contingency-based metrics
*************************

Contingency
===========

A number of metrics can be derived from a `contingency table <https://www.cawcr.gov.au/projects/verification/#Contingency_table>`_. To use this in ``climpred``, run ``.verify(metric='contingency', score=...)`` where score can be chosen from `xskillscore <https://xskillscore.readthedocs.io/en/stable/api.html#contingency-based-metrics>`_.

.. autofunction:: _contingency

Receiver Operating Characteristic
=================================

.. ipython:: python

    # Enter any of the below keywords in ``metric=...`` for the compute functions.
    print(f"Keywords: {metric_aliases['roc']}")

.. autofunction:: _roc


********************
User-defined metrics
********************

You can also construct your own metrics via the :py:class:`climpred.metrics.Metric`
class.

.. autosummary:: Metric

First, write your own metric function, similar to the existing ones with required
arguments ``forecast``, ``observations``, ``dim=None``, and ``**metric_kwargs``::

  from climpred.metrics import Metric

  def _my_msle(forecast, observations, dim=None, **metric_kwargs):
      """Mean squared logarithmic error (MSLE).
      https://peltarion.com/knowledge-center/documentation/modeling-view/build-an-ai-model/loss-functions/mean-squared-logarithmic-error."""
      # function
      return ( (np.log(forecast + 1) + np.log(observations + 1) ) ** 2).mean(dim)

Then initialize this metric function with :py:class:`climpred.metrics.Metric`::

  _my_msle = Metric(
      name='my_msle',
      function=_my_msle,
      probabilistic=False,
      positive=False,
      unit_power=0,
      )

Finally, compute skill based on your own metric::

  skill = hindcast.verify(metric=_my_msle, comparison='e2o', alignment='same_verif', dim='init')

Once you come up with an useful metric for your problem, consider contributing
this metric to `climpred`, so all users can benefit from your metric, see
`contributing <contributing.html>`_.

**********
References
**********

.. bibliography::
  :filter: docname in docnames
*******************************
Publications Using ``climpred``
*******************************

Below is a list of publications that have made use of ``climpred`` in their analysis.
We appreciate a reference to ``climpred``, e.g., in your acknowledgements section to
help build the community. Please cite :cite:p:`Brady2021`.

Feel free to open a `Pull Request <contributing.html>`_ to add your publication to the
list!

:cite:empty:`Brady2020,Krumhardt2020,Spring2020,Brady2021,Spring2021a,Spring2021b`

2021
####

.. bibliography::
  :filter: cited and year and (year == "2021") and docname in docnames

2020
####

.. bibliography::
  :filter: cited and year and (year == "2020") and docname in docnames
*******************
Reference Forecasts
*******************

To quantify the quality of an initialized forecast, it is useful to judge it against
some simple reference forecast. ``climpred`` currently supports a several reference
forecasts, and we are open to adding other reference forecasts. Consider opening a
`Pull Request <contributing.html>`_ for additional references.

**Persistence Forecast**: Whatever is observed at the time of initialization is
forecasted to persist into the forecast period :cite:p:`Jolliffe2011`.
You can compute this by passing ``reference="persistence"`` into
:py:meth:`.HindcastEnsemble.verify`, :py:meth:`.HindcastEnsemble.bootstrap`,
:py:meth:`.PerfectModelEnsemble.verify` and :py:meth:`.PerfectModelEnsemble.bootstrap`.

**Damped Persistence Forecast**: (*Not Implemented*) The amplitudes of the anomalies
reduce in time exponentially at a time scale of the local autocorrelation :cite:p:`Yuan2016`.

.. math::

    v_{dp}(t) = v(0)e^{-\alpha t}

**Climatology**: The average values at the temporal forecast resolution (e.g., annual,
monthly, daily) over some long period, which is usually 30 years :cite:p:`Jolliffe2011`.
You can compute this by passing ``reference="climatology"`` into
:py:meth:`.HindcastEnsemble.verify`, :py:meth:`.HindcastEnsemble.bootstrap`,
:py:meth:`.PerfectModelEnsemble.verify` and :py:meth:`.PerfectModelEnsemble.bootstrap`.

**Uninitialized**: Uninitialized ensembles are generated by perturbing initial
conditions only at one point in the historical run.
These are generated via micro (round-off error perturbations) or macro (starting from
completely different restart files) methods. Uninitialized ensembles are used to
approximate the magnitude of internal climate variability and to confidently extract
the forced response (ensemble mean) in the climate system. In ``climpred``, we use
uninitialized ensembles as a baseline for how important (reoccurring) initializations
are for lending predictability to the system.
You can compute this by passing ``reference="uninitialized"`` into
:py:meth:`.HindcastEnsemble.verify`, :py:meth:`.HindcastEnsemble.bootstrap`,
:py:meth:`.PerfectModelEnsemble.verify` and :py:meth:`.PerfectModelEnsemble.bootstrap`.
Some modeling centers (such as NCAR)
provide a dynamical uninitialized ensemble (the CESM Large Ensemble) along with their
initialized prediction system (the CESM Decadal Prediction Large Ensemble).
Use :py:meth:`.HindcastEnsemble.add_uninitialized` or
:py:meth:`.PerfectModelEnsemble.add_uninitialized`.
This could be, for example, output from an ``uninitialized`` Large Ensemble.
If ``uninitialzed`` isn't available, one can run
:py:meth:`.HindcastEnsemble.generate_uninitialized` or
:py:meth:`.PerfectModelEnsemble.generate_uninitialized`, which
resamples the ``initialized`` from :py:class:`.HindcastEnsemble` or
``control`` from :py:class:`.PerfectModelEnsemble` to an
``uninitialized`` forecast.

**Random Mechanism**: (*Not Implemented*) A probability distribution is assigned to the
possible range of the variable being forecasted, and a sequence of forecasts is
produced by taking a sequence of independent values from that distribution
:cite:p:`Jolliffe2011`. This would be similar to computing an ``uninitialized``
forecast.

References
##########

.. bibliography::
  :filter: docname in docnames
Overview: Why climpred?
=======================

There are many packages out there related to computing metrics on initialized
geoscience predictions. However, we didn't find any one package that unified all our
needs.

Output from earth system prediction hindcast (also called re-forecast) experiments is
difficult to work with. A typical output file could contain the dimensions
``initialization``, ``lead time``, ``ensemble member``, ``latitude``, ``longitude``,
``depth``. ``climpred`` leverages the labeled dimensions of ``xarray`` to handle the
headache of bookkeeping for you. We offer :py:class:`.HindcastEnsemble` and
:py:class:`.PerfectModelEnsemble` objects that carry products to verify against (e.g.,
control runs, reconstructions, uninitialized ensembles) along with your initialized
prediction output.

When computing lead-dependent skill scores, ``climpred`` handles all of the
``init+lead-valid_time``-matching for you, properly aligning the multiple
``time`` dimensions between the hindcast and verification datasets.
We offer a suite of vectorized `deterministic <metrics.html#deterministic>`_
and `probabilistic <metrics.html#probabilistic>`_ metrics that can be applied to time
series and grids. It's as easy as concatenating your initialized prediction output into
one :py:class:`xarray.Dataset` and running the :py:meth:`.HindcastEnsemble.verify`
command:

.. :: python

>>> HindcastEnsemble.verify(
...     metric="rmse", comparison="e2o", dim="init", alignment="maximize"
... )
****************
Related Packages
****************

We're big fans of open-source software at ``climpred`` and want to support and
collaborate with any other forecasting-related packages.
Below is a list of packages that have some place in the Earth System forecasting
world. Please reach out to us if you are aware of any open-source packages in this
domain that are not on the list.

* `Microsoft forecasting <https://microsoft.github.io/forecasting/>`_:
  A collection of best practices for time series forecasting.
* `MurCSS <https://github.com/illing2005/murcss>`_:
  A tool for standardized evaluation of decadal hindcast systems. See also
  `this <https://www-miklip.dkrz.de/about/problems/>`__ and
  `this <https://www-miklip.dkrz.de/about/murcss/>`__ page for additional resources
  on scoring from MurCSS.
* `properscoring <https://github.com/TheClimateCorporation/properscoring>`_:
  Probabilistic forecast metrics in python.
  (We've since wrapped these functions in
  `xskillscore <https://xskillscore.readthedocs.io>`_.)
* `pySTEPS <https://pysteps.github.io/>`_:
  Framework for short-term ensemble forecasting of precipitation.
* `S2D Verification <https://cran.r-project.org/web/packages/s2dverification/s2dverification.pdf>`_:
  An R package for a common set of tools for forecast verification.
* `SwirlsPy <https://docs.com-swirls.org/latest/>`_:
  Analysis and prediction of nowcasts for precipitation and weather phenomena.
* `xskillscore <https://xskillscore.readthedocs.io>`_:
  Metrics for verifying forecasts (a key dependency to ``climpred``).
* `doppyo <https://github.com/csiro-dcfp/doppyo>`_ with many metrics transferred to
  `xskillscore <https://xskillscore.readthedocs.io>`_
climpred: verification of weather and climate forecasts
=======================================================

..
    Table version of badges inspired by pySTEPS.

.. list-table::
    :stub-columns: 1
    :widths: 10 90

    * - docs
      - |docs| |joss| |doi|
    * - tests
      - |ci| |upstream| |codecov| |precommit|
    * - package
      - |conda| |conda downloads| |pypi| |pypi downloads|
    * - license
      - |license|
    * - community
      - |gitter| |contributors| |forks| |stars| |issues| |PRs|
    * - tutorials
      - |gallery| |workshop| |cloud|

.. |docs| image:: https://img.shields.io/readthedocs/climpred/stable.svg?style=flat
    :target: https://climpred.readthedocs.io/en/stable/?badge=stable
    :alt: Documentation Status

.. |joss| image:: https://joss.theoj.org/papers/246d440e3fcb19025a3b0e56e1af54ef/status.svg
    :target: https://joss.theoj.org/papers/246d440e3fcb19025a3b0e56e1af54ef
    :alt: JOSS paper

.. |doi| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4556085.svg
    :target: https://doi.org/10.5281/zenodo.4556085
    :alt: DOI

.. |ci|  image:: https://github.com/pangeo-data/climpred/workflows/climpred%20testing/badge.svg
    :target: https://github.com/pangeo-data/climpred/actions/workflows/climpred_testing.yml
    :alt: CI

.. |upstream| image:: https://github.com/pangeo-data/climpred/actions/workflows/upstream-dev-ci.yml/badge.svg
    :target: https://github.com/pangeo-data/climpred/actions/workflows/upstream-dev-ci.yml
    :alt: CI upstream

.. |codecov| image:: https://codecov.io/gh/pangeo-data/climpred/branch/main/graph/badge.svg
    :target: https://codecov.io/gh/pangeo-data/climpred
    :alt: codecov

.. |precommit| image:: https://results.pre-commit.ci/badge/github/pangeo-data/climpred/main.svg
   :target: https://results.pre-commit.ci/latest/github/pangeo-data/climpred/main
   :alt: pre-commit.ci status

.. |conda| image:: https://img.shields.io/conda/vn/conda-forge/climpred.svg
    :target: https://anaconda.org/conda-forge/climpred
    :alt: Conda Version

.. |pypi| image:: https://img.shields.io/pypi/v/climpred.svg
    :target: https://pypi.python.org/pypi/climpred/
    :alt: pypi Version

.. |license| image:: https://img.shields.io/github/license/pangeo-data/climpred.svg
    :target: ../../LICENSE.txt
    :alt: license

.. |gitter| image:: https://badges.gitter.im/Join%20Chat.svg
    :target: https://gitter.im/climpred
    :alt: gitter chat

.. |contributors| image:: https://img.shields.io/github/contributors/pangeo-data/climpred
    :alt: GitHub contributors
    :target: https://github.com/pangeo-data/climpred/graphs/contributors

.. |conda downloads| image:: https://img.shields.io/conda/dn/conda-forge/climpred
    :alt: Conda downloads
    :target: https://anaconda.org/conda-forge/climpred

.. |pypi downloads| image:: https://pepy.tech/badge/climpred
    :alt: pypi downloads
    :target: https://pepy.tech/project/climpred

.. |gallery| image:: https://img.shields.io/badge/climpred-examples-ed7b0e.svg
    :alt: climpred gallery
    :target: https://mybinder.org/v2/gh/pangeo-data/climpred/main?urlpath=lab%2Ftree%2Fdocs%2Fsource%2Fquick-start.ipynb

.. |workshop| image:: https://img.shields.io/badge/climpred-workshop-f5a252
    :alt: climpred workshop
    :target: https://mybinder.org/v2/gh/bradyrx/climpred_workshop/master

.. |cloud| image:: https://img.shields.io/badge/climpred-cloud_demo-f9c99a
    :alt: climpred cloud demo
    :target: https://github.com/aaronspring/climpred-cloud-demo

.. |forks| image:: https://img.shields.io/github/forks/pangeo-data/climpred
    :alt: GitHub forks
    :target: https://github.com/pangeo-data/climpred/network/members

.. |stars| image:: https://img.shields.io/github/stars/pangeo-data/climpred
    :alt: GitHub stars
    :target: https://github.com/pangeo-data/climpred/stargazers

.. |issues| image:: https://img.shields.io/github/issues/pangeo-data/climpred
    :alt: GitHub issues
    :target: https://github.com/pangeo-data/climpred/issues

.. |PRs| image:: https://img.shields.io/github/issues-pr/pangeo-data/climpred
    :alt: GitHub PRs
    :target: https://github.com/pangeo-data/climpred/pulls


.. note::
    We are actively looking for new contributors for climpred! Riley moved to McKinsey's
    Climate Analytics team. Aaron is finishing his PhD, but will stay in academia.
    We especially hope for python enthusiasts from seasonal, subseasonal or weather
    prediction community. In our past coding journey, collaborative coding, feedbacking
    issues and pull requests advanced our code and thinking about forecast verification
    more than we could have ever expected.
    `Aaron <https://github.com/aaronspring/>`_ can provide guidance on
    implementing new features into climpred. Feel free to implement
    your own new feature or take a look at the
    `good first issue <https://github.com/pangeo-data/climpred/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22>`_
    tag in the issues. Please reach out to us via `gitter <https://gitter.im/climpred>`_.



Installation
============

You can install the latest release of ``climpred`` using ``pip`` or ``conda``:

.. code-block:: bash

    pip install climpred[complete]

.. code-block:: bash

    conda install -c conda-forge climpred

You can also install the bleeding edge (pre-release versions) by cloning this
repository or installing directly from GitHub:

.. code-block:: bash

    git clone https://github.com/pangeo-data/climpred.git
    cd climpred
    pip install . --upgrade

.. code-block:: bash

    pip install git+https://github.com/pangeo-data/climpred.git



**Getting Started**


* :doc:`why-climpred`
* :doc:`scope`
* :doc:`quick-start`
* :doc:`examples`

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Getting Started

    why-climpred
    scope
    quick-start.ipynb
    examples

**User Guide**

* :doc:`setting-up-data`
* :doc:`initialized-datasets`
* :doc:`prediction-ensemble-object`
* :doc:`alignment`
* :doc:`metrics`
* :doc:`comparisons`
* :doc:`significance`
* :doc:`bias_removal`
* :doc:`smoothing`
* :doc:`terminology`
* :doc:`reference_forecast`

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: User Guide

    setting-up-data
    initialized-datasets
    prediction-ensemble-object.ipynb
    alignment.ipynb
    metrics
    comparisons
    significance
    bias_removal.ipynb
    smoothing
    terminology
    reference_forecast

**Help & Reference**

* :doc:`api`
* :doc:`changelog`
* :doc:`code_of_conduct`
* :doc:`contributing`
* :doc:`contributors`
* :doc:`helpful-links`
* :doc:`literature`
* :doc:`publications`
* :doc:`related-packages`
* :doc:`release_procedure`

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Help & Reference

    api
    changelog
    code_of_conduct
    contributing
    contributors
    helpful-links
    literature
    publications
    related-packages
    release_procedure
===============
Code of Conduct
===============

Our Pledge
----------

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

Our Standards
-------------

Examples of behavior that contributes to creating a positive environment
include:

-   Using welcoming and inclusive language
-   Being respectful of differing viewpoints and experiences
-   Gracefully accepting constructive criticism
-   Focusing on what is best for the community
-   Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

-   The use of sexualized language or imagery and unwelcome sexual attention or
    advances
-   Trolling, insulting/derogatory comments, and personal or political attacks
-   Public or private harassment
-   Publishing others' private information, such as a physical or electronic
    address, without explicit permission
-   Other conduct which could reasonably be considered inappropriate in a
    professional setting

Our Responsibilities
--------------------

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

Scope
-----

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

Enforcement
-----------

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team.
All complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

Attribution
-----------

This Code of Conduct is adapted from the `Contributor Covenant homepage <https://www.contributor-covenant.org/version/1/4/code-of-conduct.html>`__, and `xgcm <https://github.com/xgcm/xgcm>`__. For answers to common questions about this code of conduct, see `Contributor Covenant <https://www.contributor-covenant.org/faq>`__.
***********************
Setting Up Your Dataset
***********************

``climpred`` relies on a consistent naming system for
`xarray <https://xarray.pydata.org/en/stable/>`_ dimensions.
This allows things to run more easily under-the-hood.

:py:class:`.PredictionEnsemble` expects at the minimum to contain dimensions
``init`` and ``lead``.

``init`` is the initialization dimension, that relays the time
steps at which the ensemble was initialized.
``init`` is known as ``forecast_reference_time`` in the `CF convention <http://cfconventions.org/Data/cf-standard-names/77/build/cf-standard-name-table.html>`_.
``init`` must be of type :py:class:`pandas.DatetimeIndex`, or
:py:class:`xarray.CFTimeIndex`.
If ``init`` is of type ``int``, it is assumed to be annual data starting Jan 1st.
A UserWarning is issues when this assumption is made.

``lead`` is the lead time of the forecasts from initialization.
``lead`` is known as ``forecast_period`` in the `CF convention <http://cfconventions.org/Data/cf-standard-names/77/build/cf-standard-name-table.html>`_.
``lead`` must be ``int`` or ``float``.
The units for the ``lead`` dimension must be specified in as an attribute.
Valid options are ``["years", "seasons", "months"]`` and
``["weeks", "pentads", "days", "hours", "minutes", "seconds"]``.
If ``lead`` is provided as :py:class:`pandas.Timedelta` up to ``"weeks"``, ``lead``
is converted to ``int`` and a corresponding ``lead.attrs["units"]``.
For larger ``lead`` as :py:class:`pandas.Timedelta`
``["months", "seasons" or "years"]``, no conversion is possible.

``valid_time=init+lead`` will be calculated in :py:class:`.PredictionEnsemble` upon
instantiation.

Another crucial dimension is ``member``, which holds the various ensemble members,
which is only required for probabilistic metrics. ``member`` is known as
``realization`` in the `CF convention <http://cfconventions.org/Data/cf-standard-names/77/build/cf-standard-name-table.html>`_

Any additional dimensions will be broadcasted: these could be dimensions like ``lat``,
``lon``, ``depth``, etc.

If the expected dimensions are not found, but the matching `CF convention <http://cfconventions.org/Data/cf-standard-names/77/build/cf-standard-name-table.html>`_
``standard_name`` in a coordinate attribute, the dimension is renamed to the
corresponding ``climpred`` ensemble dimension.

Check out the demo to setup a ``climpred``-ready prediction ensemble
`from your own data <examples/misc/setup_your_own_data.html>`_ or via
`intake-esm <https://intake-esm.readthedocs.io/>`_ from `CMIP DCPP <examples/misc/setup_your_own_data.html#intake-esm-for-cmorized-output>`_.

**Verification products** are expected to contain the ``time`` dimension at the minimum.
For best use of ``climpred``, their ``time`` dimension should cover the full length of
``init`` and be the same calendar type as the accompanying prediction ensemble.
The ``time`` dimension must be :py:class:`pandas.DatetimeIndex`, or
:py:class:`xarray.CFTimeIndex`.
``time`` dimension of type ``int`` is assumed to be annual data starting Jan 1st.
A UserWarning is issued when this assumption is made.
These products can also include additional dimensions, such as ``lat``, ``lon``,
``depth``, etc.

See the below table for a summary of dimensions used in ``climpred``, and data types
that ``climpred`` supports for them.

.. list-table:: List of ``climpred`` dimension and coordinates
   :widths: 25 25 25 25 25
   :header-rows: 1

   * - Short Name
     - Types
     - Long name
     - `CF convention <http://cfconventions.org/Data/cf-standard-names/77/build/cf-standard-name-table.html>`_
     - Attribute(s)
   * - ``lead``
     - ``int``, ``float`` or :py:class:`pandas.Timedelta` up to ``weeks``
     - lead timestep after initialization ``init``
     - ``forecast_period``
     - units (str) [``years``, ``seasons``, ``months``, ``weeks``, ``pentads``, ``days``, ``hours``, ``minutes``, ``seconds``] or :py:class:`pandas.Timedelta`
   * - ``init``
     -  :py:class:`pandas.DatetimeIndex` or :py:class:`xarray.CFTimeIndex`.
     - initialization as start date of experiment
     - ``forecast_reference_time``
     - None
   * - ``member``
     - ``int``, ``str``
     - ensemble member
     - ``realization``
     - None

Probably the most challenging part is concatenating
(:py:func:`xarray.concat`) raw model output with dimension ``time`` of
multiple simulations to a multi-dimensional :py:class:`xarray.Dataset` containing
dimensions ``init``, (``member``) and ``lead``, where ``time`` becomes
``valid_time=init+lead``. One way of doing it is
:py:func:`climpred.preprocessing.shared.load_hindcast`.
***********
Comparisons
***********

Forecasts have to be verified against some product to evaluate their performance.
However, when verifying against a product, there are many different ways one can
compare the ensemble of forecasts. Here, we cover the comparison options for both
:py:class:`.HindcastEnsemble` and
:py:class:`.PerfectModelEnsemble`.
See `terminology <terminology.html>`__ for clarification on the differences between
these two experimental setups.

All high-level functions like :py:meth:`.HindcastEnsemble.verify`,
:py:meth:`.HindcastEnsemble.bootstrap`, :py:meth:`.PerfectModelEnsemble.verify` and
:py:meth:`.PerfectModelEnsemble.bootstrap` take a ``comparison`` keyword to select the
comparison style. See below for a detailed description on the differences between these
comparisons.

Hindcast Ensembles
##################

In :py:class:`.HindcastEnsemble`, the ensemble mean forecast
(``comparison="e2o"``) is expected to perform better than individual ensemble members
(``comparison="m2o"``) as the chaotic component of forecasts is expected to be
suppressed by this averaging, while the memory of the system sustains. :cite:p:`Boer2016`

.. currentmodule:: climpred.comparisons

``keyword: "e2o", "e2r"``

.. autosummary:: _e2o

``keyword: "m2o", "m2r"``

.. autosummary:: _m2o


Perfect Model Ensembles
#######################

In :py:class:`.PerfectModelEnsemble`, there are many more ways of
verifying forecasts. :cite:t:`Seferian2018` uses a comparison of all ensemble members against
the control run (``comparison="m2c"``) and all ensemble members against all other
ensemble members (``comparison="m2m"``). Furthermore, the ensemble mean forecast can
be verified against one control member (``comparison="e2c"``) or all members
(``comparison="m2e"``) as done in :cite:t:`Griffies1997`.

``keyword: "m2e"``

.. autosummary:: _m2e

``keyword: "m2c"``

.. autosummary:: _m2c

``keyword: "m2m"``

.. autosummary:: _m2m

``keyword: "e2c"``

.. autosummary:: _e2c


Normalization
#############

The goal of a normalized distance metric is to get a constant or comparable value of
typically ``1`` (or ``0`` for metrics defined as ``1 - metric``) when the metric
saturates and the predictability horizon is reached (see `metrics <metrics.html>`__).

A factor is added in the normalized metric formula :cite:p:`Seferian2018` to accomodate
different comparison styles. For example, ``metric="nrmse"`` gets smaller in
comparison ``"m2e"``.
than ``"m2m"`` by design, since the ensembe mean is always closer to individual members
than the ensemble members to each other. In turn, the normalization factor is ``2`` for
comparisons ``"m2c"``, ``"m2m"``, and ``"m2o"``. It is 1 for ``"m2e"``, ``"e2c"``, and
``"e2o"``.

Interpretation of Results
#########################

When :py:class:`.HindcastEnsemble` skill is computed over all
initializations ``dim="init"`` of the hindcast, the resulting skill is a mean forecast
skill over all initializations.

:py:class:`.PerfectModelEnsemble` skill is computed over a
supervector comprised of all
initializations and members, which allows the computation of the ACC-based skill
:cite:p:`Bushuk2018`, but also returns a mean forecast skill over all initializations.

The supervector approach shown in :cite:t:`Bushuk2018` and just calculating a distance-based
metric like ``rmse`` over the member dimension as in :cite:t:`Griffies1997` yield very similar
results.

Compute over dimension
######################

The argument ``dim`` defines over which dimension a metric is computed. We can
apply a metric over all dimensions from the ``initialized`` dataset expect ``lead``.
The resulting skill is then
reduced by this ``dim``. Therefore, applying a metric over ``dim="member"`` or
``dim=[]`` creates a skill for all initializations individually.
This can show the initial conditions dependence of skill.
Likewise when computing skill over ``"init"``, we get skill for each member.
This ``dim`` argument is different from the ``comparison`` argument which
just specifies how ``forecast`` and ``observations`` are defined.

However, this above logic applies to deterministic metrics. Probabilistic metrics need
to be applied to the ``member`` dimension and ``comparison`` from
``["m2c", "m2m"]`` in :py:meth:`.PerfectModelEnsemble.verify` and ``"m2o"`` comparison
in :py:meth:`.HindcastEnsemble.verify`.

``dim`` should not contain ``member`` when the comparison already computes ensemble
means as in ``["e2o", "e2c"]``.


User-defined comparisons
########################

You can also construct your own comparisons via the
:py:class:`climpred.comparisons.Comparison` class.

.. autosummary:: Comparison

First, write your own comparison function, similar to the existing ones. If a
comparison should also be used for probabilistic metrics, make sure that probabilistic
metrics returns ``forecast`` with ``member`` dimension and
``observations`` without. For deterministic metrics, return ``forecast`` and
``observations`` with identical dimensions but without an identical comparison::

  from climpred.comparisons import Comparison, M2M_MEMBER_DIM

  def _my_m2median_comparison(initialized, metric=None):
      """Identical to m2e but median."""
      observations_list = []
      forecast_list = []
      supervector_dim = "member"
      for m in initialized.member.values:
          forecast = initialized.drop_sel(member=m).median("member")
          observations = initialized.sel(member=m).squeeze()
          forecast_list.append(forecast)
          observations_list.append(observations)
      observations = xr.concat(observations_list, M2M_MEMBER_DIM)
      forecast = xr.concat(forecast_list, M2M_MEMBER_DIM)
      forecast[M2M_MEMBER_DIM] = np.arange(forecast[M2M_MEMBER_DIM].size)
      observations[M2M_MEMBER_DIM] = np.arange(observations[M2M_MEMBER_DIM].size)
      return forecast, observations

Then initialize this comparison function with
:py:class:`climpred.comparisons.Comparison`::

  __my_m2median_comparison = Comparison(
      name="m2me",
      function=_my_m2median_comparison,
      probabilistic=False,
      hindcast=False)

Finally, compute skill based on your own comparison::

  PerfectModelEnsemble.verify(
    metric="rmse",
    comparison=__my_m2median_comparison,
    dim=[],
  )

Once you come up with an useful comparison for your problem, consider contributing this
comparison to ``climpred``, so all users can benefit from your comparison, see
`contributing <contributing.html>`_.


References
##########

.. bibliography::
  :filter: docname in docnames
.. include:: ../../CHANGELOG.rst
API Reference
=============

This page provides an auto-generated summary of ``climpred``'s API.
For more details and examples, refer to the relevant chapters in the main part of the
documentation.

High-Level Classes
------------------

.. currentmodule:: climpred.classes


A primary feature of ``climpred`` is our prediction ensemble objects,
:py:class:`.HindcastEnsemble` and :py:class:`.PerfectModelEnsemble`.
Users can add their initialized ensemble to these classes, as well as verification
products (assimilations, reconstructions, observations), control runs, and uninitialized
ensembles.

PredictionEnsemble
~~~~~~~~~~~~~~~~~~

:py:class:`.PredictionEnsemble` is the base class for :py:class:`.HindcastEnsemble` and
:py:class:`.PerfectModelEnsemble`. :py:class:`.PredictionEnsemble` cannot be called
directly, but :py:class:`.HindcastEnsemble` and :py:class:`.PerfectModelEnsemble`
inherit the common base functionality.

.. autosummary::
    :toctree: api/

    PredictionEnsemble
    PredictionEnsemble.__init__

-------
Builtin
-------

.. autosummary::
    :toctree: api/

    PredictionEnsemble.__len__
    PredictionEnsemble.__iter__
    PredictionEnsemble.__delitem__
    PredictionEnsemble.__contains__
    PredictionEnsemble.__add__
    PredictionEnsemble.__sub__
    PredictionEnsemble.__mul__
    PredictionEnsemble.__truediv__
    PredictionEnsemble.__getitem__
    PredictionEnsemble.__getattr__

----------
Properties
----------

.. autosummary::
    :toctree: api/

    PredictionEnsemble.coords
    PredictionEnsemble.nbytes
    PredictionEnsemble.sizes
    PredictionEnsemble.dims
    PredictionEnsemble.chunks
    PredictionEnsemble.chunksizes
    PredictionEnsemble.data_vars
    PredictionEnsemble.equals
    PredictionEnsemble.identical


HindcastEnsemble
~~~~~~~~~~~~~~~~

A :py:class:`.HindcastEnsemble` is a prediction ensemble that is
initialized off of some form of observations (an assimilation, reanalysis, etc.). Thus,
it is anticipated that forecasts are verified against observation-like products. Read
more about the terminology `here <terminology.html>`_.

.. autosummary::
    :toctree: api/

    HindcastEnsemble
    HindcastEnsemble.__init__

-------------------------
Add and Retrieve Datasets
-------------------------

.. autosummary::
    :toctree: api/

    HindcastEnsemble.add_observations
    HindcastEnsemble.add_uninitialized
    HindcastEnsemble.get_initialized
    HindcastEnsemble.get_observations
    HindcastEnsemble.get_uninitialized

------------------
Analysis Functions
------------------

.. autosummary::
    :toctree: api/

    HindcastEnsemble.verify
    HindcastEnsemble.bootstrap

-------------
Generate Data
-------------

.. autosummary::
    :toctree: api/

    HindcastEnsemble.generate_uninitialized

--------------
Pre-Processing
--------------

.. autosummary::
    :toctree: api/

    HindcastEnsemble.smooth
    HindcastEnsemble.remove_bias
    HindcastEnsemble.remove_seasonality

-------------
Visualization
-------------

.. autosummary::
    :toctree: api/

    HindcastEnsemble.plot
    HindcastEnsemble.plot_alignment


PerfectModelEnsemble
~~~~~~~~~~~~~~~~~~~~

A ``PerfectModelEnsemble`` is a prediction ensemble that is initialized off of a
control simulation for a number of randomly chosen initialization dates. Thus,
forecasts cannot be verified against real-world observations.
Instead, they are `compared <comparisons.html>`_ to one another and to the
original control run. Read more about the terminology `here <terminology.html>`_.

.. autosummary::
    :toctree: api/

    PerfectModelEnsemble
    PerfectModelEnsemble.__init__

-------------------------
Add and Retrieve Datasets
-------------------------

.. autosummary::
    :toctree: api/

    PerfectModelEnsemble.add_control
    PerfectModelEnsemble.get_initialized
    PerfectModelEnsemble.get_control
    PerfectModelEnsemble.get_uninitialized

------------------
Analysis Functions
------------------

.. autosummary::
    :toctree: api/

    PerfectModelEnsemble.verify
    PerfectModelEnsemble.bootstrap

-------------
Generate Data
-------------

.. autosummary::
    :toctree: api/

    PerfectModelEnsemble.generate_uninitialized

--------------
Pre-Processing
--------------

.. autosummary::
    :toctree: api/

    PerfectModelEnsemble.smooth
    PerfectModelEnsemble.remove_seasonality

-------------
Visualization
-------------

.. autosummary::
    :toctree: api/

    PerfectModelEnsemble.plot


Direct Function Calls
---------------------

While not encouraged anymore, a user can directly call functions in ``climpred``.
This requires entering more arguments, e.g. the initialized ensemble directly as
well as a verification product. Our object
:py:class:`.HindcastEnsemble` and
:py:class:`.PerfectModelEnsemble` wrap most of these functions, making
the analysis process much simpler.

Bootstrap
~~~~~~~~~
.. currentmodule:: climpred.bootstrap

.. autosummary::
    :toctree: api/

    bootstrap_compute
    bootstrap_hindcast
    bootstrap_perfect_model
    bootstrap_uninit_pm_ensemble_from_control_cftime
    bootstrap_uninitialized_ensemble
    dpp_threshold
    varweighted_mean_period_threshold

Prediction
~~~~~~~~~~
.. currentmodule:: climpred.prediction

.. autosummary::
    :toctree: api/

    compute_hindcast
    compute_perfect_model


Reference
~~~~~~~~~
.. currentmodule:: climpred.reference

.. autosummary::
    :toctree: api/

    compute_persistence
    compute_persistence_from_first_lead
    compute_uninitialized
    compute_climatology

Horizon
~~~~~~~
.. currentmodule:: climpred.horizon

.. autosummary::
    :toctree: api/

    horizon

Statistics
~~~~~~~~~~
.. currentmodule:: climpred.stats

.. autosummary::
    :toctree: api/

    decorrelation_time
    dpp
    varweighted_mean_period
    rm_poly
    rm_trend

Tutorial
~~~~~~~~
.. currentmodule:: climpred.tutorial

.. autosummary::
    :toctree: api/

    load_dataset

Preprocessing
~~~~~~~~~~~~~
.. currentmodule:: climpred.preprocessing.shared

.. autosummary::
    :toctree: api/

    load_hindcast
    rename_to_climpred_dims
    rename_SLM_to_climpred_dims
    set_integer_time_axis

.. currentmodule:: climpred.preprocessing.mpi

.. autosummary::
    :toctree: api/

    get_path


Smoothing
~~~~~~~~~
.. currentmodule:: climpred.smoothing

.. autosummary::
    :toctree: api/

    temporal_smoothing
    spatial_smoothing_xesmf

Visualization
~~~~~~~~~~~~~
.. currentmodule:: climpred.graphics

.. autosummary::
    :toctree: api/

    plot_bootstrapped_skill_over_leadyear
    plot_ensemble_perfect_model
    plot_lead_timeseries_hindcast


Metrics
-------

For a thorough look at our metrics library, please see the
`metrics </metrics.html>`_ page.

.. currentmodule:: climpred.metrics

.. autosummary::
    :toctree: api/

    Metric
    Metric.__init__
    Metric.__repr__
    _get_norm_factor
    _pearson_r
    _pearson_r_p_value
    _effective_sample_size
    _pearson_r_eff_p_value
    _spearman_r
    _spearman_r_p_value
    _spearman_r_eff_p_value
    _mse
    _rmse
    _mae
    _median_absolute_error
    _nmse
    _nmae
    _nrmse
    _msess
    _mape
    _smape
    _uacc
    _std_ratio
    _conditional_bias
    _unconditional_bias
    _bias_slope
    _msess_murphy
    _crps
    _crpss
    _crpss_es
    _brier_score
    _threshold_brier_score
    _rps
    _discrimination
    _reliability
    _rank_histogram
    _contingency
    _roc
    _spread
    _mul_bias
    _less

Comparisons
-----------

For a thorough look at our metrics library, please see the
`comparisons </comparisons.html>`_ page.

.. currentmodule:: climpred.comparisons

.. autosummary::
    :toctree: api/

    Comparison
    Comparison.__init__
    Comparison.__repr__
    _e2o
    _m2o
    _m2m
    _m2e
    _m2c
    _e2c

Config
------
Set options analogous to
`xarray <http://xarray.pydata.org/en/stable/generated/xarray.set_options.html>`_.

.. currentmodule:: climpred.options

.. autosummary::
    :toctree: api/

    set_options
﻿climpred.classes.PredictionEnsemble.\_\_iter\_\_
================================================

.. currentmodule:: climpred.classes

.. automethod:: PredictionEnsemble.__iter__
climpred.metrics.\_mse
======================

.. currentmodule:: climpred.metrics

.. autofunction:: _mse
climpred.metrics.\_rps
======================

.. currentmodule:: climpred.metrics

.. autofunction:: _rps
climpred.classes.PerfectModelEnsemble.get\_initialized
======================================================

.. currentmodule:: climpred.classes

.. automethod:: PerfectModelEnsemble.get_initialized
climpred.metrics.\_spearman\_r\_eff\_p\_value
=============================================

.. currentmodule:: climpred.metrics

.. autofunction:: _spearman_r_eff_p_value
﻿climpred.comparisons.Comparison.\_\_init\_\_
=============================================

.. currentmodule:: climpred.comparisons

.. automethod:: Comparison.__init__
climpred.classes.PerfectModelEnsemble.smooth
============================================

.. currentmodule:: climpred.classes

.. automethod:: PerfectModelEnsemble.smooth
climpred.classes.PredictionEnsemble.dims
========================================

.. currentmodule:: climpred.classes

.. autoproperty:: PredictionEnsemble.dims
climpred.classes.HindcastEnsemble.add\_observations
===================================================

.. currentmodule:: climpred.classes

.. automethod:: HindcastEnsemble.add_observations
climpred.classes.PredictionEnsemble.sizes
=========================================

.. currentmodule:: climpred.classes

.. autoproperty:: PredictionEnsemble.sizes
climpred.metrics.\_spearman\_r
==============================

.. currentmodule:: climpred.metrics

.. autofunction:: _spearman_r
climpred.classes.HindcastEnsemble.\_\_init\_\_
==============================================

.. currentmodule:: climpred.classes

.. automethod:: HindcastEnsemble.__init__
climpred.preprocessing.shared.rename\_to\_climpred\_dims
========================================================

.. currentmodule:: climpred.preprocessing.shared

.. autofunction:: rename_to_climpred_dims
climpred.classes.PerfectModelEnsemble.\_\_init\_\_
==================================================

.. currentmodule:: climpred.classes

.. automethod:: PerfectModelEnsemble.__init__
climpred.metrics.\_spread
=========================

.. currentmodule:: climpred.metrics

.. autofunction:: _spread
climpred.classes.PredictionEnsemble.chunks
==========================================

.. currentmodule:: climpred.classes

.. autoproperty:: PredictionEnsemble.chunks
climpred.classes.PredictionEnsemble.identical
=============================================

.. currentmodule:: climpred.classes

.. automethod:: PredictionEnsemble.identical
climpred.classes.HindcastEnsemble.generate\_uninitialized
=========================================================

.. currentmodule:: climpred.classes

.. automethod:: HindcastEnsemble.generate_uninitialized
climpred.classes.PredictionEnsemble
===================================

.. currentmodule:: climpred.classes

.. autoclass:: PredictionEnsemble


   .. automethod:: __init__


   .. rubric:: Methods

   .. autosummary::

      ~PredictionEnsemble.__init__
      ~PredictionEnsemble.equals
      ~PredictionEnsemble.get_initialized
      ~PredictionEnsemble.get_uninitialized
      ~PredictionEnsemble.identical
      ~PredictionEnsemble.plot
      ~PredictionEnsemble.remove_seasonality
      ~PredictionEnsemble.smooth





   .. rubric:: Attributes

   .. autosummary::

      ~PredictionEnsemble.chunks
      ~PredictionEnsemble.chunksizes
      ~PredictionEnsemble.coords
      ~PredictionEnsemble.data_vars
      ~PredictionEnsemble.dims
      ~PredictionEnsemble.mathType
      ~PredictionEnsemble.nbytes
      ~PredictionEnsemble.sizes
climpred.metrics.\_msess
========================

.. currentmodule:: climpred.metrics

.. autofunction:: _msess
climpred.metrics.\_roc
======================

.. currentmodule:: climpred.metrics

.. autofunction:: _roc
climpred.classes.PerfectModelEnsemble
=====================================

.. currentmodule:: climpred.classes

.. autoclass:: PerfectModelEnsemble


   .. automethod:: __init__


   .. rubric:: Methods

   .. autosummary::

      ~PerfectModelEnsemble.__init__
      ~PerfectModelEnsemble.add_control
      ~PerfectModelEnsemble.bootstrap
      ~PerfectModelEnsemble.generate_uninitialized
      ~PerfectModelEnsemble.get_control
      ~PerfectModelEnsemble.get_initialized
      ~PerfectModelEnsemble.get_uninitialized
      ~PerfectModelEnsemble.plot
      ~PerfectModelEnsemble.smooth
      ~PerfectModelEnsemble.verify
﻿climpred.comparisons.\_e2o
===========================

.. currentmodule:: climpred.comparisons

.. autofunction:: _e2o
climpred.preprocessing.shared.rename\_SLM\_to\_climpred\_dims
=============================================================

.. currentmodule:: climpred.preprocessing.shared

.. autofunction:: rename_SLM_to_climpred_dims
climpred.metrics.\_median\_absolute\_error
==========================================

.. currentmodule:: climpred.metrics

.. autofunction:: _median_absolute_error
climpred.metrics.\_conditional\_bias
====================================

.. currentmodule:: climpred.metrics

.. autofunction:: _conditional_bias
climpred.classes.PredictionEnsemble.\_\_getattr\_\_
===================================================

.. currentmodule:: climpred.classes

.. automethod:: PredictionEnsemble.__getattr__
climpred.classes.HindcastEnsemble.smooth
========================================

.. currentmodule:: climpred.classes

.. automethod:: HindcastEnsemble.smooth
climpred.preprocessing.shared.set\_integer\_time\_axis
======================================================

.. currentmodule:: climpred.preprocessing.shared

.. autofunction:: set_integer_time_axis
climpred.classes.PerfectModelEnsemble.bootstrap
===============================================

.. currentmodule:: climpred.classes

.. automethod:: PerfectModelEnsemble.bootstrap
climpred.options.set\_options
=============================

.. currentmodule:: climpred.options

.. autoclass:: set_options


   .. automethod:: __init__


   .. rubric:: Methods

   .. autosummary::

      ~set_options.__init__
climpred.classes.HindcastEnsemble.get\_uninitialized
====================================================

.. currentmodule:: climpred.classes

.. automethod:: HindcastEnsemble.get_uninitialized
climpred.metrics.\_reliability
==============================

.. currentmodule:: climpred.metrics

.. autofunction:: _reliability
climpred.metrics.\_crpss
========================

.. currentmodule:: climpred.metrics

.. autofunction:: _crpss
climpred.classes.HindcastEnsemble.plot
======================================

.. currentmodule:: climpred.classes

.. automethod:: HindcastEnsemble.plot
climpred.metrics.\_std\_ratio
=============================

.. currentmodule:: climpred.metrics

.. autofunction:: _std_ratio
climpred.bootstrap.bootstrap\_uninitialized\_ensemble
=====================================================

.. currentmodule:: climpred.bootstrap

.. autofunction:: bootstrap_uninitialized_ensemble
climpred.horizon.horizon
========================

.. currentmodule:: climpred.horizon

.. autofunction:: horizon
climpred.metrics.\_brier\_score
===============================

.. currentmodule:: climpred.metrics

.. autofunction:: _brier_score
climpred.prediction.compute\_hindcast
=====================================

.. currentmodule:: climpred.prediction

.. autofunction:: compute_hindcast
climpred.classes.PredictionEnsemble.data\_vars
==============================================

.. currentmodule:: climpred.classes

.. autoproperty:: PredictionEnsemble.data_vars
climpred.stats.rm\_trend
========================

.. currentmodule:: climpred.stats

.. autofunction:: rm_trend
﻿climpred.comparisons.\_e2c
===========================

.. currentmodule:: climpred.comparisons

.. autofunction:: _e2c
climpred.metrics.\_rank\_histogram
==================================

.. currentmodule:: climpred.metrics

.. autofunction:: _rank_histogram
climpred.preprocessing.shared.load\_hindcast
============================================

.. currentmodule:: climpred.preprocessing.shared

.. autofunction:: load_hindcast
climpred.prediction.compute\_perfect\_model
===========================================

.. currentmodule:: climpred.prediction

.. autofunction:: compute_perfect_model
climpred.classes.HindcastEnsemble.remove\_seasonality
=====================================================

.. currentmodule:: climpred.classes

.. automethod:: HindcastEnsemble.remove_seasonality
climpred.stats.rm\_poly
=======================

.. currentmodule:: climpred.stats

.. autofunction:: rm_poly
﻿climpred.comparisons.Comparison
================================

.. currentmodule:: climpred.comparisons

.. autoclass:: Comparison


   .. automethod:: __init__


   .. rubric:: Methods

   .. autosummary::

      ~Comparison.__init__
climpred.classes.PerfectModelEnsemble.get\_uninitialized
========================================================

.. currentmodule:: climpred.classes

.. automethod:: PerfectModelEnsemble.get_uninitialized
climpred.stats.varweighted\_mean\_period
========================================

.. currentmodule:: climpred.stats

.. autofunction:: varweighted_mean_period
﻿climpred.metrics.Metric.\_\_repr\_\_
=====================================

.. currentmodule:: climpred.metrics

.. automethod:: Metric.__repr__
climpred.classes.PredictionEnsemble.chunksizes
==============================================

.. currentmodule:: climpred.classes

.. autoproperty:: PredictionEnsemble.chunksizes
climpred.classes.PredictionEnsemble.\_\_init\_\_
================================================

.. currentmodule:: climpred.classes

.. automethod:: PredictionEnsemble.__init__
climpred.preprocessing.mpi.get\_path
====================================

.. currentmodule:: climpred.preprocessing.mpi

.. autofunction:: get_path
climpred.metrics.\_unconditional\_bias
======================================

.. currentmodule:: climpred.metrics

.. autofunction:: _unconditional_bias
﻿climpred.metrics.Metric
========================

.. currentmodule:: climpred.metrics

.. autoclass:: Metric


   .. automethod:: __init__


   .. rubric:: Methods

   .. autosummary::

      ~Metric.__init__
climpred.bootstrap.bootstrap\_uninit\_pm\_ensemble\_from\_control\_cftime
=========================================================================

.. currentmodule:: climpred.bootstrap

.. autofunction:: bootstrap_uninit_pm_ensemble_from_control_cftime
climpred.reference.compute\_persistence
=======================================

.. currentmodule:: climpred.reference

.. autofunction:: compute_persistence
climpred.metrics.\_rmse
=======================

.. currentmodule:: climpred.metrics

.. autofunction:: _rmse
climpred.metrics.\_threshold\_brier\_score
==========================================

.. currentmodule:: climpred.metrics

.. autofunction:: _threshold_brier_score
climpred.metrics.\_mape
=======================

.. currentmodule:: climpred.metrics

.. autofunction:: _mape
climpred.stats.decorrelation\_time
==================================

.. currentmodule:: climpred.stats

.. autofunction:: decorrelation_time
climpred.classes.PerfectModelEnsemble.get\_control
==================================================

.. currentmodule:: climpred.classes

.. automethod:: PerfectModelEnsemble.get_control
climpred.classes.PerfectModelEnsemble.generate\_uninitialized
=============================================================

.. currentmodule:: climpred.classes

.. automethod:: PerfectModelEnsemble.generate_uninitialized
climpred.classes.HindcastEnsemble.bootstrap
===========================================

.. currentmodule:: climpred.classes

.. automethod:: HindcastEnsemble.bootstrap
climpred.stats.dpp
==================

.. currentmodule:: climpred.stats

.. autofunction:: dpp
climpred.metrics.\_crpss\_es
============================

.. currentmodule:: climpred.metrics

.. autofunction:: _crpss_es
climpred.classes.PerfectModelEnsemble.remove\_seasonality
=========================================================

.. currentmodule:: climpred.classes

.. automethod:: PerfectModelEnsemble.remove_seasonality
﻿climpred.comparisons.Comparison.\_\_repr\_\_
=============================================

.. currentmodule:: climpred.comparisons

.. automethod:: Comparison.__repr__
climpred.classes.HindcastEnsemble
=================================

.. currentmodule:: climpred.classes

.. autoclass:: HindcastEnsemble


   .. automethod:: __init__


   .. rubric:: Methods

   .. autosummary::

      ~HindcastEnsemble.__init__
      ~HindcastEnsemble.add_observations
      ~HindcastEnsemble.add_uninitialized
      ~HindcastEnsemble.bootstrap
      ~HindcastEnsemble.get_initialized
      ~HindcastEnsemble.get_observations
      ~HindcastEnsemble.get_uninitialized
      ~HindcastEnsemble.plot
      ~HindcastEnsemble.remove_bias
      ~HindcastEnsemble.smooth
      ~HindcastEnsemble.verify
climpred.metrics.\_bias\_slope
==============================

.. currentmodule:: climpred.metrics

.. autofunction:: _bias_slope
climpred.metrics.\_mae
======================

.. currentmodule:: climpred.metrics

.. autofunction:: _mae
climpred.reference.compute\_climatology
=======================================

.. currentmodule:: climpred.reference

.. autofunction:: compute_climatology
climpred.tutorial.load\_dataset
===============================

.. currentmodule:: climpred.tutorial

.. autofunction:: load_dataset
climpred.metrics.\_mul\_bias
============================

.. currentmodule:: climpred.metrics

.. autofunction:: _mul_bias
﻿climpred.comparisons.\_m2m
===========================

.. currentmodule:: climpred.comparisons

.. autofunction:: _m2m
climpred.classes.PredictionEnsemble.\_\_len\_\_
===============================================

.. currentmodule:: climpred.classes

.. automethod:: PredictionEnsemble.__len__
﻿climpred.comparisons.\_m2o
===========================

.. currentmodule:: climpred.comparisons

.. autofunction:: _m2o
climpred.classes.HindcastEnsemble.get\_observations
===================================================

.. currentmodule:: climpred.classes

.. automethod:: HindcastEnsemble.get_observations
climpred.metrics.\_spearman\_r\_p\_value
========================================

.. currentmodule:: climpred.metrics

.. autofunction:: _spearman_r_p_value
climpred.classes.PerfectModelEnsemble.add\_control
==================================================

.. currentmodule:: climpred.classes

.. automethod:: PerfectModelEnsemble.add_control
climpred.graphics.plot\_bootstrapped\_skill\_over\_leadyear
===========================================================

.. currentmodule:: climpred.graphics

.. autofunction:: plot_bootstrapped_skill_over_leadyear
﻿climpred.comparisons.\_m2c
===========================

.. currentmodule:: climpred.comparisons

.. autofunction:: _m2c
climpred.bootstrap.bootstrap\_compute
=====================================

.. currentmodule:: climpred.bootstrap

.. autofunction:: bootstrap_compute
climpred.bootstrap.bootstrap\_perfect\_model
============================================

.. currentmodule:: climpred.bootstrap

.. autofunction:: bootstrap_perfect_model
﻿climpred.comparisons.\_m2e
===========================

.. currentmodule:: climpred.comparisons

.. autofunction:: _m2e
climpred.metrics.\_smape
========================

.. currentmodule:: climpred.metrics

.. autofunction:: _smape
climpred.graphics.plot\_ensemble\_perfect\_model
================================================

.. currentmodule:: climpred.graphics

.. autofunction:: plot_ensemble_perfect_model
climpred.bootstrap.varweighted\_mean\_period\_threshold
=======================================================

.. currentmodule:: climpred.bootstrap

.. autofunction:: varweighted_mean_period_threshold
climpred.metrics.\_nmae
=======================

.. currentmodule:: climpred.metrics

.. autofunction:: _nmae
climpred.metrics.\_contingency
==============================

.. currentmodule:: climpred.metrics

.. autofunction:: _contingency
climpred.bootstrap.dpp\_threshold
=================================

.. currentmodule:: climpred.bootstrap

.. autofunction:: dpp_threshold
﻿climpred.metrics.Metric.\_\_init\_\_
=====================================

.. currentmodule:: climpred.metrics

.. automethod:: Metric.__init__
climpred.classes.PredictionEnsemble.\_\_add\_\_
===============================================

.. currentmodule:: climpred.classes

.. automethod:: PredictionEnsemble.__add__
climpred.metrics.\_pearson\_r\_p\_value
=======================================

.. currentmodule:: climpred.metrics

.. autofunction:: _pearson_r_p_value
climpred.classes.HindcastEnsemble.get\_initialized
==================================================

.. currentmodule:: climpred.classes

.. automethod:: HindcastEnsemble.get_initialized
climpred.smoothing.spatial\_smoothing\_xesmf
============================================

.. currentmodule:: climpred.smoothing

.. autofunction:: spatial_smoothing_xesmf
climpred.metrics.\_nrmse
========================

.. currentmodule:: climpred.metrics

.. autofunction:: _nrmse
climpred.metrics.\_discrimination
=================================

.. currentmodule:: climpred.metrics

.. autofunction:: _discrimination
climpred.classes.PredictionEnsemble.nbytes
==========================================

.. currentmodule:: climpred.classes

.. autoproperty:: PredictionEnsemble.nbytes
climpred.classes.HindcastEnsemble.remove\_bias
==============================================

.. currentmodule:: climpred.classes

.. automethod:: HindcastEnsemble.remove_bias
climpred.classes.PredictionEnsemble.\_\_sub\_\_
===============================================

.. currentmodule:: climpred.classes

.. automethod:: PredictionEnsemble.__sub__
climpred.metrics.\_less
=======================

.. currentmodule:: climpred.metrics

.. autofunction:: _less
climpred.metrics.\_uacc
=======================

.. currentmodule:: climpred.metrics

.. autofunction:: _uacc
climpred.classes.PredictionEnsemble.\_\_contains\_\_
====================================================

.. currentmodule:: climpred.classes

.. automethod:: PredictionEnsemble.__contains__
climpred.metrics.\_nmse
=======================

.. currentmodule:: climpred.metrics

.. autofunction:: _nmse
climpred.classes.HindcastEnsemble.plot\_alignment
=================================================

.. currentmodule:: climpred.classes

.. automethod:: HindcastEnsemble.plot_alignment
climpred.metrics.\_msess\_murphy
================================

.. currentmodule:: climpred.metrics

.. autofunction:: _msess_murphy
climpred.graphics.plot\_lead\_timeseries\_hindcast
==================================================

.. currentmodule:: climpred.graphics

.. autofunction:: plot_lead_timeseries_hindcast
climpred.smoothing.temporal\_smoothing
======================================

.. currentmodule:: climpred.smoothing

.. autofunction:: temporal_smoothing
climpred.classes.PredictionEnsemble.\_\_delitem\_\_
===================================================

.. currentmodule:: climpred.classes

.. automethod:: PredictionEnsemble.__delitem__
climpred.classes.HindcastEnsemble.add\_uninitialized
====================================================

.. currentmodule:: climpred.classes

.. automethod:: HindcastEnsemble.add_uninitialized
climpred.metrics.\_get\_norm\_factor
====================================

.. currentmodule:: climpred.metrics

.. autofunction:: _get_norm_factor
climpred.classes.PerfectModelEnsemble.verify
============================================

.. currentmodule:: climpred.classes

.. automethod:: PerfectModelEnsemble.verify
climpred.metrics.\_crps
=======================

.. currentmodule:: climpred.metrics

.. autofunction:: _crps
climpred.metrics.\_pearson\_r
=============================

.. currentmodule:: climpred.metrics

.. autofunction:: _pearson_r
climpred.reference.compute\_uninitialized
=========================================

.. currentmodule:: climpred.reference

.. autofunction:: compute_uninitialized
climpred.classes.PredictionEnsemble.equals
==========================================

.. currentmodule:: climpred.classes

.. automethod:: PredictionEnsemble.equals
climpred.metrics.\_effective\_sample\_size
==========================================

.. currentmodule:: climpred.metrics

.. autofunction:: _effective_sample_size
climpred.classes.HindcastEnsemble.verify
========================================

.. currentmodule:: climpred.classes

.. automethod:: HindcastEnsemble.verify
climpred.classes.PredictionEnsemble.plot
========================================

.. currentmodule:: climpred.classes

.. automethod:: PredictionEnsemble.plot
climpred.classes.PerfectModelEnsemble.plot
==========================================

.. currentmodule:: climpred.classes

.. automethod:: PerfectModelEnsemble.plot
climpred.bootstrap.bootstrap\_hindcast
======================================

.. currentmodule:: climpred.bootstrap

.. autofunction:: bootstrap_hindcast
climpred.reference.compute\_persistence\_from\_first\_lead
==========================================================

.. currentmodule:: climpred.reference

.. autofunction:: compute_persistence_from_first_lead
climpred.classes.PredictionEnsemble.\_\_truediv\_\_
===================================================

.. currentmodule:: climpred.classes

.. automethod:: PredictionEnsemble.__truediv__
climpred.classes.PredictionEnsemble.coords
==========================================

.. currentmodule:: climpred.classes

.. autoproperty:: PredictionEnsemble.coords
climpred.classes.PredictionEnsemble.\_\_mul\_\_
===============================================

.. currentmodule:: climpred.classes

.. automethod:: PredictionEnsemble.__mul__
climpred.metrics.\_pearson\_r\_eff\_p\_value
============================================

.. currentmodule:: climpred.metrics

.. autofunction:: _pearson_r_eff_p_value
climpred.classes.PredictionEnsemble.\_\_getitem\_\_
===================================================

.. currentmodule:: climpred.classes

.. automethod:: PredictionEnsemble.__getitem__
