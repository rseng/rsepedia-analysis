[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5905676.svg)](https://doi.org/10.5281/zenodo.5905676)
[![Research Software Directory Badge](https://img.shields.io/badge/rsd-storyboards-00a3e3.svg)](https://www.research-software.nl/software/storyboards)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=eucp-project_storyboards&metric=alert_status)](https://sonarcloud.io/summary/new_code?id=eucp-project_storyboards)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu)

# Storyboards

This repository contains the source code for the
[EUCP](https://www.eucp-project.eu/) [Storyboard
website](https://eucp-project.github.io/storyboards/).

The website has been built with [Nuxt](https://nuxtjs.org), using
[nuxt-content](https://content.nuxtjs.org/) for authoring stories and
[tailwindcss](https://tailwindcss.com/docs/installation) + [tailwind
typography](https://tailwindcss.com/docs/typography-plugin) for styling. It is
hosted on [GitHub pages](https://nuxtjs.org/deployments/github-pages/).

## Writing a story

All stories are stored in the `static/stories` folder. Each story consists of a
markdown file and a folder with some assets (images, etc.) belonging to the
story. In addition to standard markdown, the `:::Chapter{}` directive is used to
break the story into parts that can be displayed individually.

An example story might look like this:

`static/stories/example-story.md`
```markdown
---
title: Example story
author: Peter Kalverla et al., Netherlands eScience Center
thumbnail: "intro.png"
category: EUCP data and products
trl: high
id: 13
---
:::Chapter{headline="Introduction" image="intro.png"}
## This is the first Chapter
You can format text using markdown.

The headline property will be used for the chapter navigation blocks.

The image property will be used for the main display image of this chapter.

Even though it's called 'image', you can also add standalone HTML pages, such as
an exported mapbox file.

The three colons below mark the end of the first chapter.
:::

:::Chapter{headline="Methods" image="concept.png"}
## This is the second Chapter

and so on...
:::
```

This produces the following layout:

![Screenshot of example-story](example-story.png)

The images should be stored in a directory with the same name as the story, but
with a leading underscore, like so:

```bash
- static/
  - stories/
    - example-story.md
    - _example-story/
      - intro.png
      - concept.png
```

The frontmatter (title, author, etc.) will be used to show the story on the
stories overview page. The ID and TRL (technical readiness level) properties are
currently not used, but they are still here for legacy reasons.

## Adding your story to our collection

If you want your story to be included on
[eucp-project.github.io/storyboards](https://eucp-project.github.io/storyboards),
you can make a pull request to this repository. We will review it and if
everything is okay, we'll merge the story into the main branch. Continous
deployment will then automatically update the site.

If you want to add a story but are unsure about the github workflow, please
don't hesitate to get in touch. We are happy to help.

## Serving the site locally

You can also make a local build of the site, if you want to check that your
story is formatted correctly before making a pull request. The following
instructions are the default instructions from a new nuxt project. After cloning
the repository:

```bash
# install dependencies
$ npm install

# serve with hot reload at localhost:3000
$ npm run dev

# build for production and launch server
$ npm run build
$ npm run start

# generate static project
$ npm run generate
```

For detailed explanation on how things work, check out the [documentation](https://nuxtjs.org).

# Reusing the storyboards format for a different project

The source code (excluding the stories content) is licenced under Apache 2. You
can fork this repo and add your own content, modify the styling, and do whatever
you want. We'd appreciate it if you inform us about your re-using the software.
We're also happy to help setting it up for you.

# Reusing the storyboard materials
The content of the storyboards is licenced under CC-BY 4.0. Please don't
hesitate to contact the storyboard authors if you're interested in their work.
---
id: 10
trl: medium
category: EUCP science and methodologies
title: Multiple lines of evidence
author: T. Crocker et al., MetOffice
thumbnail: "intro_lines.png"
---

:::Chapter{headline="Introduction and aim" image="intro_lines.png"}
## Multiple lines of evidence
Users of climate information now have many sources of climate data to choose
from.

The different sources have different characteristics such as resolution which
means that some sources are better suited to particular applications.

However, the different datasets sample different sources of uncertainty and
therefore will offer different uncertainty ranges.

Although high resolution data (e.g. from convection permitting models) is often
preferred by users, the number of model simulations available is often limited
compared to lower resolution simulations.

Can information from other data sources, or lines of evidence, be used to better
inform the interpretations of results from impact studies?
:::

:::Chapter{headline="Introduction and aim" image="intro_users.png"}
## Multiple lines of evidence
**EUCP is trying to help users by assembling data from across different datasets
for common regions and providing guidance on interpreting those differences.**

This is important because looking across the different datasets provide really
important context for interpreting results based on a particularly dataset.

**The following case study shows how this might look in practice.**
:::

:::Chapter{headline="Example case study" image="casestudy.png"}
## Example case study: Landslide risk in Romania
Landslides are a common hazard in Romania with impacts on property and
infrastructure. Rainfall is a key triggering factor in landslide events and
there is interest around how landslide risk in the future will change under
different climate projections.

Reference:

[Balteanu, D., Chendes, V., & Sima, M. (2009). GIS landslide hazard map of
Romania. GIM International, 23(4),
13-15.](https://www.gim-international.com/content/article/gis-landslide-hazard-map-of-romania)
:::

:::Chapter{headline="Recent research" image="recent_research.png"}
## Recent research into future climate changes and their impact on landslide risk.

[Niculiţă, M. (2020). Landslide Hazard Induced by Climate Changes in
North-Eastern Romania. Climate Change Management, (May),
245-265.](https://doi.org/10.1007/978-3-030-37425-9_13)

This study finds that climate projections imply a return to higher rainfall
intensities that were last seen in the 1970s  and 1980s and which coincided with
increased landslide events.

However, the study is based on a subset of the CORDEX multi-model ensemble,
consisting of just 10 models.

Can further information that puts these models in the context of the wider
CORDEX ensemble, as well as other climate model ensembles, help with
interpreting the results of this study?
:::

:::Chapter{headline="Examining the information" image="examining_info.png"}
## Examining the information from multiple lines of evidence
By visualising data from other sources of climate projections, it is possible to
interpret the results of the study in a wider context and gain greater
information about the uncertainties involved in interpreting the results.

Each point in the figure represents the % change in the multi annual mean of JJA
(June, July, August) of precipitation in Romania from the period 1996-2005 to
2041-2050.

Separating out the models used as drivers for downscaling helps with
understanding how representative the driving models are of the wider ensemble,
but further analysis is useful for understanding the impact downscaling itself…
:::

:::Chapter{headline="Impact of dynamical downscaling" image="impact.png"}
## Impact of dynamical downscaling
The impact of downscaling the CMIP5 models is mixed, causing both increases and
decreases in the precipitation change signal relative to the driving GCM
depending on the RCM used. The variance is more positively skewed though, with
larger changes on the wetter side. Choosing other CMIP5 models to downscale
could conceivably produce wetter or drier results than the existing CORDEX
range.

For the two EUCP CPM models, the impact of downscaling the driving CORDEX models
gives either a similar or a wetter result compared to the driving model.
:::

:::Chapter{headline="Summary" image="summary.png"}
## Summary
Examining results from multiple lines of evidence makes it easier to put the
results of a study in context.

In this example, information suggests that the models used in the study are a
reasonable representation of the most likely outcomes indicated by most other
lines of evidence.

However, there is a significant tail risk of a drier future outcome than
suggested by the models used in the study. There is also a less likely tail risk
of wetter outcomes.
:::
---
id: 9
trl: high
category: Application of EUCP innovations
title: Estimating regionalized hydrological impacts of climate change over Europe
author: F. Sperna Weiland et al., Deltares
thumbnail: panel_1.png
---
:::Chapter{headline="Constructing a multi-model ensemble" image="panel_1.png"}

The EURO-CORDEX initiative developed an ensemble of Regional Climate Model (RCM)
simulations at a horizontal resolution of ~11 km - in total 20 RCMs were used
here.

This study builds upon a set of readily available European scale hydrological
simulations. The variation between the models provides insights on the influence
of bias-correction, calibration and hydrological model selected.

Three hydrological models were selected:

1. **LISFLOOD** from JRC, the model is calibrated and the CORDEX data is
   bias-corrected
2. **CWatM** from IIASA, the model is calibrated on a regional scale, no
   bias-correction is applied
3. **wflow_sbm**, the model is parameterized based on soil, land use and
   vegetation characteristics, no bias-correction is applied

Changes in discharge are assessed between the period 1981-2010 and 2031-2060 for
RCP8p5.
:::

:::Chapter{headline="9 river basins with different characteristics" image="panel_2.png"}

- Past studies indicated future wetter conditions for Northern Europe and drier
  conditions for Southern Europe.
- Discharge changes will likely depend on climatological conditions and
  catchment characteristics
- The discharge response to rainfall changes is non linear and thus requires
- modelling
:::


:::Chapter{headline="Historical discharge regimes" image="panel_3_rearranged.png"}
- Each column represents simulations run by a different model
- The LISFLOOD simulations match the observations (in black) best. This is a
  result of bias-correction and calibration
- The similarity between CWatM and Wflow is larger. There is a wide spread
  between the different ensemble members because no bias-correction was applied.
:::

:::Chapter{headline="Projected changes – minimum flows" image="panel_4.png"}
- There are clear indications for decreases in minimum flow in Southern European
  basins and increases in Northern basins.
- The consensus between the different hydrological models and climate models is
  large for the Northern basins (Glomma, Angermanalven) and the Southern basins
  (Ebro and Tanaro)
- For the Central European basins the signal is mixed and the choice of
  hydrological model and RCM does influence the impact estimation.
:::

:::Chapter{headline="Performance-based weighting" image="panel_5.png"}
Performance based weighting is used to obtain a more robust change signal from
the ensemble of future projections.

Projected changes are displayed for weighted change in annual mean discharge
(mean Q), annual maximum discharge (maxQ) and 7 day minimum discharge (NM7Q)

Two weighting methods have been applied:

  * **ClimWIP**: Climate model Weighting by Independence and Performance. This
    method evaluates the historical performance of GCMs and their independence.
    Focus of the weighting was on meteorological variables
  * **REA**: Reliability Ensemble Averaging. This method evaluates the
    historical performance of the models and the distance of its change
    projections to the future projected ensemble weighted change. Focus of the
    weighting was on river discharges near the outlet of the basin.

For part of the basins the influence of the weighting was limited. There is
little difference between the weighted and non-weigthed ensemble mean change.

For other basins the REA method deviates from the other change projections. We
assume there is a strong influence of the relatively large historical bias of
the RCM based simulations from observations. A single RCM with the smallest
historical bias is highly favored by the method.
:::

:::Chapter{headline="More information" image="panel_6.png"}
# Further reading

[Sperna Weiland F.C., R.F. Visser, P. Greve, B. Bisselin, L. Brunner, A.H.
Weerts, Frontiers in Water, 3, 143, doi:
10.3389/frwa.2021.713537](https://doi.org/10.3389/frwa.2021.713537)
:::
---
id: 5
trl: high
category: Application of EUCP innovations
title: Sandy beach erosion induced by sea level rise
author: P. Athanasiou, et al. Deltares.
thumbnail: overview.png
---

:::Chapter{headline="Introduction" image="overview.png"}
## Sandy beach erosion induced by sea level rise

- Almost 41% of European Union’s population lives near the coast
- Big part of European coastline comprises erodible sandy stretches
- Sea level rise (SLR) will cause shoreline retreat
- Need to identify vulnerable coastal areas to focus adaptation measures
- Large scale coastal erosion assessments can aid in identifying hotspots
:::

:::Chapter{headline="Framework" image="mechanisms.png"}
## Framework for coastal land loss projections

Global probabilistic sea level rise projections were used under the RCP 4.5 and
RCP 8.5 scenarios. For RCP 8.5, a high-end RCP 8.5 scenario based on IPCC AR5
but with a higher contribution from Antarctic and Greenland ice sheets was used.

Two different maps of the location of sandy beaches and two different approaches
for their nearshore slopes are used to account for epistemic input data
uncertainty.

Shoreline retreat is projected in the feature at sandy beach locations on 1 km
alongshore grid around the European coastline.

Sandy beach land loss is calculated by aggregating the total shoreline retreat
per NUTS3 region.
:::

:::Chapter{headline="Shoreline retreat" image="rcpSenarios.png"}
## Shoreline retreat in Europe

Using satellite derived sandy beach locations and spatially varying nearshore
slopes data, projections of shoreline retreat are made.

A European averaged median shoreline retreat of 31 m (22 m) is projected under
RCP 8.5 (4.5) by year 2050, relative to the baseline year 2010.

These values change to 97 m (54 m) under RCP 8.5 (4.5) by year 2100.

Spatial variability along the European coastline is pronounced depending on the
nearshore slopes variability.

This retreat would translate to 2,500 km<sup>2</sup> (1,400 km<sup>2</sup>) of
coastal land loss by 2100.
:::

:::Chapter{headline="Land-loss at regions scale" image="regions.png"}
## Coastal land-loss at European regions scale

Aggregating the projected shoreline retreat per NUTS3 region and normalizing by
the region’s coastline length, helps to identify erosion hotspots around the
European coastline

Highly vulnerable regions are identified on the
- Italian Adriatic coast
- the French Atlantic coast
- Belgium
- The Netherlands
- Denmark
- Lithuania
- Latvia
:::

:::Chapter{headline="Uncertainty analysis" image="uncertainty.png"}
## Uncertainty analysis

A variance-based sensitivity analysis was used to quantify the uncertainty of
the coastal land loss estimates. This method measures the contribution of each
of the uncertain input parameters; 1) nearshore slope data, 2) sandy beach
location data, 3) RCP scenario, 4) Sea level rise (SLR).

This allows for quantifying and comparing the uncertainties related to the
choice of input geophysical datasets (epistemic) versus the inherent uncertainty
associated with RCPs and SLR projections.

On average epistemic uncertainty accounts for 45% (26%) by 2050 (2100), with
high variability between countries

This provides insights on areas that geophysical datasets should be updated with
monitoring and observations.
:::


:::Chapter{headline="Reference" image="logos.png"}
## More information

Athanasiou, P., van Dongeren, A., Giardino, A., Vousdoukas, M. I., Ranasinghe,
R., and Kwadijk, J. (2020). Uncertainties in projections of sandy beach erosion
due to sea level rise: an analysis at the European scale. Sci. Rep. 10, 1–14.
doi: [10.1038/s41598-020-68576-0](https://doi.org/10.1038/s41598-020-68576-0).

The dataset corresponding to this story is:

Projections of sea level rise induced shoreline retreat in Europe
(https://doi.org/10.4121/uuid:8e73cab0-960b-46a8-bf67-ee0eadcc1e7d).
:::
---
id: 13
trl: medium
category: EUCP data and products
title: Atlas of (un)constrained climate projections
author: B. Booth et al., UK Met Office
thumbnail: "comparing.png"
---

:::Chapter{headline=Introduction image=ipcc.png}
## Observationally constrained projections gaining greater visibility
Observations have been used to constrain climate projection spread in National Climate Scenarios (e.g. Swiss, UK, Australian) but have not yet been commonly applied elsewhere.  They represent a way to assess which existing climate projections are more plausible given how well they reproduce historical climate.  They are starting to see wider adoption and can be an important tool in the climate projection toolbox.

Observational constraints used in IPCC to provide context, given many more high
temperature response simulations that are not thought to be more likely.

> For the first time in an IPCC report, assessed future changes in global
> surface temperature, ocean warming and sea level are constructed by combining
> multi-model projections with observational constraints based on past simulated
> warming, as well as the AR6 assessment of climate sensitivity.
> -- <cite>Summary for Policy Makers, AR6, IPCC</cite>
:::

:::Chapter{headline="Underlying understanding" image=methods.png}
## EUCP has led underlying work to understand user value in these approaches
The work within EUCP has been on understanding different potential quantitative approaches to using historical validation of current projections. We provided the first common evaluation of their projected changes, across a range of currently available methodologies.

We have assessed the reliability of these approaches using blind out-of-sample evaluation, using new CMIP6 simulations as a proxy for real world historical and future responses. This has enabled us to look at the skill of each method across a range of metrics, which provides an objective assessment of a particular methodological approach.
:::

:::Chapter{headline="Assessing the value added" image=summary.png}
## Assessing the value added

An important output from our EUCP work has been to develop an objective assessment of the value added by using a methodology that uses observations to validate projections. We coordinated a blind common out-of-sample assessment of each method's skill, using new CMIP6 simulations as pseudo-observations and then validated the constrained projections using the future CMIP6 responses.

This is a tough test of the methods as CMIP6 includes process not represented in earlier climate projection and includes some projected changes that lie outside this range. This enabled us to assess to identified where these methods provided more reliable estimates than just taking the raw climate model range.

We have been able to summarise the implied skill of each method. Based on this
out-of-sample analysis, we find that:

* Use of observations leads to improved temperature projections for all methods
* However, little evidence of consistent skill for rainfall

These results also point to where using a particular methodology can actually make the results worse than just using the CMIP spread. Rainfall projections are worse using ASK, whereas CALL showed poor results in the Mediterranean for temperature. These issues relate to the particular properties of the methodologies that make them more or less suitable for some applications.

The summary of the out-of-sample assessment, presented here, provides a guide to where applying one of these methodologies is likely to add value.
:::

:::Chapter{headline=Atlas image=comparing.png}
## Atlas of (un)constrained projections
Given the added skill in these methods, we have been working to make this data
more widely available.

The [EUCP Atlas](https://eucp-project.github.io/) provides a portal for
preprocessed projections for Europe. It allows exploration of the different
methods, comparison betweem multiple views, and answers to frequently asked
questions.

The Atlas has its own DOI:
[10.5281/zenodo.5654741](https://doi.org/10.5281/zenodo.5654741). Source code is
[available on GitHub](https://github.com/eucp-project/).
:::


:::Chapter{headline=Guidance image=example_user.png}
### Example use cases
To provide guidance on how to use the constrained projections, we have worked
out 4 hypothetical use cases.

These use cases are also available in the atlas.
:::

:::Chapter{headline="Data availability" image=zenodo.png}
### Data availability

We have gone to great lengths to provide the data from different methods in a
standardized format. Where possible, we followed CMOR standards and
CF-conventions.

The data is available through Zenodo and can be cited using the information
provided there. Its DOI is
[10.5281/zenodo.5645153](https://doi.org/10.5281/zenodo.5645153)
:::
---
id: 6
trl: high
category: EUCP science and methodologies
title: Skillful decadal prediction of southern European summer temperature
authors: L. F. Borchert, et al. Institut Pierre Simon Laplace (IPSL).
thumbnail: "correlation.png"
---

:::Chapter{headline="Problem" image="correlation.png"}
## Outline of the problem

Predictions of European summer temperature (EUST; defined as the box in the
right figure) 2 to 9 years ahead (decadal predictions) can provide useful
information for decision makers. A promising tool for performing such
predictions are global climate models, which simulate interactions between
ocean, atmosphere, ice and vegetation, and the human impact on the climate
system. In recent climate models, decadal predictions agree well with observed
EUST during the period 1970-2014.

Much of this observed skill in decadal EUST predictions stems from
anthropogenic, volcanic and solar forcing to the climate system. But naturally
occurring unforced internal climate variations can also have a pronounced impact
on changes of EUST. Unfortunately, unforced variability of EUST is not directly
predictable in current climate models. Our aim is to suggest a statistical fix
to the climate models' inability to predict unforced EUST on the decadal time
scale.
:::

:::Chapter{headline="Idea and approach" image="steps.png"}
## The idea behind this study

The struggle of climate models to predict unforced decadal variations of
temperature over continents is well-documented. Skilful dynamical model
predictions for unforced variations are mostly reported for oceanic properties
such as heat content or sea surface temperature (SST). Our idea: identify an
observed (statistical) connection between the dynamically predictable ocean and
dynamically unpredictable EUST to provide skilful dynamical-statistical decadal
prediction of unforced EUST. Such a dyn-stat prediction would rely on skilful
decadal unforced SST prediction to infer statistically EUST predictions using
SST-EUST pathways.

## The residuals approach
These investigations require an estimate of observed unforced variability. This
variability is difficult to estimate. One approach is to construct the mean of
historical model simulations with global climate models. These are different
except in their external forcing. Their mean can separate the response of the
climate system to forcing, so the forced response can be subtracted from the
observed climate signal. This residuals method was introduced by Smith et al.
(2019). Here, we use it for the first time to also examine unforced interactions
in the climate system. We construct the residuals from an average over 28
different climate models with a total of more than 250 simulations from the
Coupled Model Intercomparison Project Phase 6 (CMIP 6).

Reference

[Smith, D.M., Eade, R., Scaife, A.A. et al. (2019) Robust skill of decadal
climate predictions. npj Clim Atmos Sci 2,
13.](https://doi.org/10.1038/s41612-019-0071-y)
:::

:::Chapter{headline="An unforced observed link" image="linkSST.png"}
## An unforced observed link between North Atlantic SST and EUST

We examine residuals for the period 1900-1969 in HadISST and HadCRU
observations. This period was chosen to avoid overlap with the prediction test
period 1970-2014. We find an unforced link between spring North Atlantic SST and
EUST.

Particularly the eastern North Atlantic and Mediterranean region stands out in
terms of unforced spring SST influence on EUST. However, another important
criterion for the construction of our dyn-stat model is that the SST region is
predictable. Comparing a prediction skill map for unforced spring SST with the
areas that significantly impact EUST identifies Eastern North Atlantic -
Mediterranean (ENAMED; blue outlines in the figure) SST as a promising index for
dyn-stat prediction of unforced EUST.
:::

:::Chapter{headline="Physics behind the link" image="linkSLP.png"}
## Short excursion: the physics behind the link

Analysing sea level pressure (SLP) to get an indication of the atmospheric
dynamics involved in transporting the unforced SST-EUST signal indicates little
dynamical influence beyond local heating from the ocean. It is therefore likely
that local ocean heat release together with thermodynamic processes accomplish
the unforced seasonal ENAMED SST - EUST link.
:::

:::Chapter{headline="Dynamical-statistical predictions " image="dynamicalPred.png"}
## Dynamical-statistical predictions of unforced EUST

In our framework, a prediction of unforced EUST is constructed from a dynamical
prediction of unforced spring ENAMED SST, which is then rescaled based on the
local historically observed summer temperature variance at each grid point. Such
predictions outperform purely dynamical predictions over most of southern
Europe, leading to significant prediction skill. In some areas, the skill
increase is even significant.

Our findings highlight the value of using the skill of dynamical decadal ocean
predictions to harvest skill for the unforced signal of EUST variations 2-9
years ahead.
:::

:::Chapter{headline="More information" image="info.png"}
## More information

This work was published open access in Environmental Research Letters:
[Borchert, L.F., V. Koul, M.B. Menary, D.J. Befort, D. Swingedouw, G. Sgubin, J.
Mignot (2021) Skillful decadal prediction of unforced southern European summer
temperature variations. Environ. Res. Lett. 16
104017](https://doi.org/10.1088/1748-9326/ac20f5)
:::
---
id: 3
trl: low
category: EUCP science and methodologies
title: Comparing methods to constrain future climate projections
author: L. Brunner, et al. ETH Zurich.
thumbnail: constraining_papers_collage.png
---

:::Chapter{headline="Introduction" image="constraining_papers_collage.png"}
## Background

Many different methods exist to quantify and constrain projections of future
climate, differing in both their evaluation against observed changes and their
underlying methodological assumptions.

In order to investigate if different methods agree with each other, a
coordinated testing framework is needed. Otherwise, method results reported in
the scientific literature are often not comparable because they differ not only
in their approach (which we want to analyse) but also in their setup such as:

* variable (e.g., temperature vs precipitation)
* region (e.g., global vs Europe)
* season (e.g., annual vs summer)
* time period (e.g., 2041-2060 vs 2080-2099)
* type of change (e.g., absolute vs relative change)
* models used (including ensemble members used)
* uncertainties included (e.g., internal variability vs forced response only)
* reported results (e.g., likely vs interquartile range)
:::

:::Chapter{headline="Framework" image="Brunner2020_fig1.png"}
## A common framework for method comparison
Within EUCP we brought together 8 groups working on quantifying uncertainty in
future climate projections and developed a common framework to enable a fair
comparison:

* temperature and precipitation
* Europe, SREX regions, 4 grid points
* summer (June, July, August)
* change in 2041-60 relative to 1995-2014
* CMIP5 models under RCP8.5
* same model pool (if possible)
* including 20-year internal variability
* median, 50%, and 80% range
* 2.5x2.5 horizontal resolution
* ocean masked using grid cell center
:::

:::Chapter{headline="Methods and institutions" image="Brunner2020_tab1.png"}
## Participating methods and Institutions

Here is a list of all the participating methods and institutions.
:::

:::Chapter{headline="CEU temperature change" image="Brunner2020_fig2c.gif"}
## Central European summer temperature change

* Most methods show a slightly lower median warming than in the unconstrained
  case
* Most methods show a reduction in spread
* Some divergence in the constrained distributions from the different methods
  (see next panel)
* Also the unconstrained distributions differ  (common model pool **not** used
  here)
:::

:::Chapter{headline=""CEU temperature change: same model pool"" image="EUCP_newsletter.png"}
## Central European summer temperature change: using the same model pool

* Using the same model pool reduces differences in the unconstrained
  distribution to a minimum. Remaining difference for HistC:
    1. calculation of percentiles (empirically vs Gaussian)
    2. slightly different handling of internal variability
* Methods consistently narrow the uncertainty range and agree on slightly less warming
:::

:::Chapter{headline="Remaining differences and Conclusions" image="Brunner2020c_tab3.png"}
## Remaining differences in the underlying method properties
There are remaining differences in the underlying method properties. See the
table in the left.

## Conclusions and outlook
* A common testing framework allows to consistently compare constrained
  projections from several different methods
* Methods agree on several aspects but disagree on others
* Agreement depends on the target (variable, region, mean/extreme changes)
* Some divergence in the constrained projections from different methods is
  expected due to differences in their underlying properties

- How can we choose between methods (if they disagree)?
    1. Considering the decision context (mean or extreme changes important?)
    2. Combining different methods (e.g., Hegerl et al. 2021)
    3. Selecting methods based on an objective skill measure (e.g., O’Reilly et
    al. in preparation)
:::
---
id: 8
trl: high
category: Application of EUCP innovations
title: Alpine flash floods
author: M. Zander, et al. Deltares.
thumbnail: intro_eucp.png
---

:::Chapter{headline="Introduction-EUCP data" image="intro_eucp.png"}
## Introduction and Aim: EUCP data

### WP3 simulations with high-resolution convection-permitting regional climate models

- Better representation of extreme rainfall (Ban et al., 2021)
- Heavy to severe rainfall would happen more often in the future (Pichelli et
  al., 2021)

### So what would this mean?

Translate high-resolution convection-permitting climate models to hydrological
impacts
:::

:::Chapter{headline="Introduction-Flash floods" image="intro_flash.png"}
## Introduction and Aim: Flash floods in the Alps

- Dangerous fast sudden floods
- Triggered by extreme local rainstorms
- Very deadly!

### Question
- If extreme rainfall changes, how will this affect flash flood risks?

Photo's: damage after storm hits Southerm Tyrol, Italy causing flash flooding in
the Adige catchment 4 August 2012. [Photo credits: Firefighter Department,
published by Alto
Adige](https://www.altoadige.it/foto/locale/nubifragi-in-alto-adige-allagamenti-frane-e-molti-danni-ecco-le-immagini-1.2938053#21).
:::

:::Chapter{headline="Methodology" image="method.png"}
## Methodology
Simulate high-resolution distributed hydrological model for 4 periods:

- validation: ERA5
- historical climate: climate model ERA-Interim 2000-2012
- current climate: climate model 1998-2007
- future climate: climate model 2095-2105
:::

:::Chapter{headline="Results-simulation" image="simulation.png"}
## Results

- Hydrological model plus climate model perform reasonable.
- Able to simulate recorded historical flash floods in the modelling domain.

Graph: observed discharge and simulated discharge for the ERA5 and ERA-Interim
driven simulations at Bellinzona station on the Ticino.
:::

:::Chapter{headline="Results-performance" image="erainterim_ukmo_kge_13.html"}
## Results
- Kling-Gupta efficiency for 130 gauging stations in the modelled domain for the
  ERA-Interim driven simulation. Assessment period 2002-2012 for daily
  discharges.

Rhine at Basel: KGE 0.62

Model chain of ERA-Interim driven convection-permitting climate model and
distributed hydrological model gives reasonable performance. Recorded flash
floods from databases are reproduced in the modelling results (not shown).
:::

:::Chapter{headline="Results-frequency" image="frequency.png"}
## Results - Contrasting response summer and autumn - frequency

- Adapted definition of [Amponsah et al
  (2018)](https://doi.org/10.5194/essd-10-1783-2018) for flash floods:
  - Maximal size of affected 3000 km2
  - Peak specific discharge is at least 0.5 m3s-1km-2 discharge divided by
    upstream area
:::

:::Chapter{headline="Results-intensity" image="intensity.png"}
## Results - Contrasting response summer and autumn - intensity

Magnitudes of peak specific discharge for the catchments when threshold for
flash flood is reached.

- Boxplot: 25th, 50th and 75th percentile
- Whiskers: 2.5th and 97.5th percentile
- Diamonds: maximum

Adapted definition of [Amponsah et al
(2018)](https://doi.org/10.5194/essd-10-1783-2018) for flash floods:

- Maximal size of affected 3000 km2
- Peak specific discharge is at least 0.5 m3s-1km-2 discharge divided by
  upstream area
:::

:::Chapter{headline="Conclusions" image="references.png"}
## Conclusions and Future work
- Frequency decrease in summer, but more severe extremes
- Frequency increase in autumn with more severe extremes
- Ensemble study
---
id: 4
trl: high
category: EUCP data and products
title: Infrastructure in support of EUropean Climate Prediction
author: P. Kalverla et al. Netherlands eScience Center.
thumbnail: "steroids.png"
---

:::Chapter{headline="Introduction" image="steroids.png"}
## EUCP infrastructure

EUCP aims to provide digital infrastructure in support of its overall objective
of delivering accessible, authoritative and actionable climate information.

The infrastructure should facilitate aggregation and blending of 'raw' climate
data from various sources into high-level data products that are useful for
end-users. To this end, EUCP develops methods for merging and constraining
climate model output, and we strive to make these methods easily reusable.

Storage and compute resources from [SURF](https://www.surf.nl/en) are used for
data exchange during the project, and to develop an initial prototype of the
envisioned infrastructure. For long-term sustainability and public access, we
are leveraging the European climate infrastructure developed by
[IS-ENES](https://portal.enes.org/services).

ESMValTool is seen as a key component of this infrastructure, providing a
protocol for sharing climate analytics workflows, simplifying data access, and
facilitating the development of new methods.
:::

:::Chapter{headline="Shared infrastructure" image="shared_infra.png"}
## IS-ENES

IS-ENES is the InfraStructure for the European Network for Earth System
modelling. Currently, this infrastructure is maintained and further developed
under the IS-ENES3 project.

IS-ENES provides access to climate data through ESGF. It also provides compute
resources, with many commonly used software libraries pre-installed.

In IS-ENES3, data access options are greatly improved, with a complete redesign
of the Climate4Impact portal, as well as access through JupyterHub, among
others.

## IS-ENES & EUCP

EUCP leverages the IS-ENES infrastructure to provide access to its data and
methods. DCPP and CPM data produced within EUCP have been or will be uploaded to
ESGF. Some of our methods have already been integrated into the ESMValTool. We
have also extended the ESMValTool with a Python API so as to make it usable from
within the IS-ENES Jupyter environments.

## More information

- [ENES portal](https://portal.enes.org/)
- [IS-ENES3](https://is.enes.org/)
- [JupyterHub by DKRZ](https://jupyterhub.dkrz.de)
- [ESGF](https://esgf.llnl.gov/mission.html)
:::

:::Chapter{headline="ESMValTool" image="esmvaltool.png"}
## ESMValTool

EMSValTool is *A community diagnostic and performance metrics tool for routine
evaluation of earth system models in CMIP.* It comes with a large collection of
reusable "recipes" that, by their nature, document a climate analytics workflow.
These recipes are quality controlled with both a technical and a scientific
review, and well documented.

By providing common pre-processing functions, ESMValTool helps to standardize
scientific workflows. The tool also tracks provenance information. As such,
ESMValTool facilitates a FAIR and open research process.

## EUCP contributions

New recipes:

* [ClimWIP](https://docs.esmvaltool.org/en/latest/recipes/recipe_climwip.html): Climate model weighting by independence and performance
* [KCS](https://docs.esmvaltool.org/en/latest/recipes/recipe_kcs.html): KNMI climate scenarios
* [...](https://github.com/MetOffice/EUCP_WP5_Lines_of_Evidence): Evaluation of ensemble projections accross projects/scales (work in progress)

Added functionality:

* Added a Python API (see next tab)
* Better support for scientific regions
* Improved multimodel statistics and added ensemble statistics (WIP)
* Improved support for CORDEX data and added support for FPS-convection data (WIP)
* Contributed to, and taught, the [ESMValTool tutorial](https://esmvalgroup.github.io/ESMValTool_Tutorial/)

## More information

* [ESMValTool documentation](https://docs.esmvaltool.org/)
* [ESMValTool on GitHub](https://github.com/ESMValGroup/ESMValTool)
* [ESMValTool on Zenodo](https://zenodo.org/record/5140083)
* [ESMValTool tutorial](https://esmvalgroup.github.io/ESMValTool_Tutorial/)
:::

:::Chapter{headline="EMSValTool on steroids" image="api.png"}
## Bringing ESMValTool to the JupyterLab

By adding a Python API, we have made it possible to use ESMValTool in an
interactive Python session. The user can browse and search through
the available recipes, inspect their documentation, execute the recipes, and
immediately work with the generated output.

This enables a much more intuitive research and development process, where
ESMValTool preprocessor functions can be used from the start, while
usecase-specific code can be written directly in the notebook, until it is ready
to be refactored into an ESMValTool diagnostic script.

The Python API also makes it possible to use ESMValTool in a Jupyter session.
The JupyterHub access developed in IS-ENES3 provides a very user-friendly
interface to the ESGF data and compute resources. ESMValTool now extends this
research infrastructure with access to existing analytics workflows and
standardized pre-processing pipelines.

## More information

* [Documentation](https://docs.esmvaltool.org/projects/esmvalcore/en/latest/api/esmvalcore.api.html)
* [Example notebooks and setup instructions](https://github.com/ESMValGroup/ESMValTool-JupyterLab)
* [Walkthrough video of example noteboook](https://www.youtube.com/watch?v=pY9gWckRQYs)
:::

:::Chapter{headline="Output web pages" image="html_output.png"}
## Output web pages

The new Python API paved the way for automatic generation of HTML reports,
including:

* Description
* Authors/maintainers
* Projects
* References
* Output figures
* Output files
* Provenance information

These reports have been developed under IS-ENES3 funding. Development of a
revamped ESMValTool result browser is also foreseen.
:::

:::Chapter{headline="Data products" image="data_products.png"}
## Custom data products

While the Jupyter interface is very suitable for researchers and developers, end
users typically have different requirements. Therefore, in EUCP, we are
developing a number of prototype web applications that can inform future service
products.

For example, we are working on an [atlas](https://eucp-project.github.io/atlas/)
of constrained vs unconstrained climate projections derived from different EUCP
methods. We will add guidance and further documentation, and iterate on the
design in a process of co-development with users.

While the content of this atlas is currently not coupled to the shared
infrastructure, it would be a natural extension. In fact, multiple groups are
currently experimenting with adding custom *front ends* to ESMValTool output, to
enable interactive output products.

:::Chapter{headline="Storyboards" image="storyboards2.png"}
## Project overview at a glance

Research output comes in many forms. Datasets, scientific publications, improved
models, software implementations, figures, video's, increased inderstanding,
example (web) applications, (technical) reports, et cetera.

Even if all information is publicly available, it is not always easy to find.
This inhibits obtaining a quick yet comprehensive overview of project outputs.

With these storyboards, we attempt to provide such an overview of important
project outputs. The storyboards focus on high-level information, linking to
external sources for more information where appropriate.

We hope they will prove useful.
:::

:::Chapter{headline="Impact" image="impact.png"}
## A long road to impact

A shared infrastructure provides access to both data and compute, and with
libraries like ESMValTool we see the dawn of a new age, where reusability and
reproducibility are axiomatic. This provides a new and (to some extent)
standardized way of sharing codified methods, enabling an efficient scientific
exchange mechanism.

The ESMValTool community has come a long way in establishing a robust
development process which includes automated testing and peer-review of both
technical and scientific aspects.

We believe this is an opportunity for climate services. The developments
described in this storyboard show a pathway from exploratory analysis to routine
evaluation, and from routine evaluation to actionable climate information. In a
sense, it looks similar to a continous integration framework. It is interesting
to consider how this framework can best be exploited in an operational setting.

## Feedback

If you have any additions, corrections, suggestions, thoughts, ideas, or other
kinds of input or feedback, please contact p.kalverla@esciencecenter.nl.
:::
---
id: 14
trl: low
category: EUCP data and products
title: Benefits and added value of convection-permitting climate modeling over Fenno-Scandinavia
author: Petter Lind et al., SMHI
thumbnail: modelling.png
---

:::Chapter{headline="Convection-permitting modeling in Fenno-Scandinavia" image="modelling.png"}
## Benefits and added value of convection-permitting climate modeling over Fenno-Scandinavia
Snow, ice, the long coastline, and the Scandes mountains, they all make the
Fenno-Scandinavia a complex region. The box plots show observed distribution of
daily mean 2 meter temperature (T2m, top panel) and precipitation (Pr, bottom)
for winter (DJF, left half) and summer (JJA, right half) seasons.

## A mini-ensemble
A mini-ensemble is used to explore the climate in this region. The modelling
region is illustrated in the left figure and the configuration is listed below:
- HCLIM-ALADIN @ 12km resolution, HCLIM-AROME @ 3km
- ERA-I driven runs; 1998–2018
- GCM driven runs; 1985–2005 / 2040–2060 / 2080–2100
- 2 forcing GCMs: EC-Earth & GFDL-CM3 (CMIP5)
- 2 emission scenarios: RCP8.5 and RCP4.5
:::

:::Chapter{headline="Added value of CPM model" image="precipitation.png"}
## Why the use of CPM model?
- For precipitation, there are obvious arguments in favor of using CPM models
  over Fenno-Scandinavia (and elsewhere).
- Left panel: Probability density function (PDF) of summer (JJA) hourly
  precipitation over Sweden comparing HCLIM3 CPM and HCLIM12 regional climate
  model (RCM) with radar observations.
- Right panel: Maps of 99.9th percentile of summer (JJA) hourly precipitation
  over Sweden in HCLIM12, HCLIM3 and radar observations.
:::

:::Chapter{headline="More precipitation extremes in scenarios" image="extreme.png"}
## More precipitation extremes in scenarios
- End-of-C (2080-2100), RCP8.5 changes in precipitation distributions;
- Map plots: CPM change signal in average (left) and heavy (right) summer (JJA)
  precipitation,
- Line plots: ASoP analysis for present climate (top) and changes (scn - ctrl)
  (bottom) of summer (JJA) hourly precipitation (Pr).
- Distinct shifts in precipitation distributions - reductions in low
  intensities, increases in higher intensities.
- More emphasized change in CPM compared to RCM.
:::

:::Chapter{headline="Seasonally distinct precipitation responses" image="seasonal.png"}
## More precipitation extremes in scenarios
- Changes (in % per degree warming) of high precipitation (Pr) percentiles
  (hourly data) in winter and summer over Fenno-Scandinavia.
- Stronger GCM influence in winter when RCM/CPM combinations are relatively
  close.
- Differences in RCM and CPM responses are much more pronounced in summer when
  forcing from boundaries is not as constraining (compared to DJF) and local
  convection becomes more prominent.
:::

:::Chapter{headline="Running to the hills" image="hills.png"}
## Running to the hills...
- Ambition to set focus on the Scandes mountains; if and how CPM may provide
  added value.
- These figures show changes in snow climate variables (x axes) as a function of
  altitude (y axes), by end-of-C in RCP4.5 and 8.5.
- Solid Precipitation (top left), solid Precipitation/total Precipitation ratio
  (top right), snow cover (bottom left), # days with large snow depth (bottom
  right).
:::

:::Chapter{headline="Orographic precipitation" image="orographic.png"}
## Orographic precipitation in winter
- Separating windward/leeward model responses in the Scandes indicate systematic
  differences between CPM and RCM.
- If a consistent response it will be an important aspect of expected regional
  climate change in Scandes.
- Further investigation needed.
:::

:::Chapter{headline="Other variables" image="winds.png"}
## Not just precipitation?
- Added value from CPM is evident for precipitation - what about other
  variables?
- We could expect improved representations of winds (e.g. convective wind gusts,
  orographic down-slope wind storms, thermal circulations).

Figures show ongoing work:
- Top panels: PDFs of winter near surface wind speed, comparing CPM (red) and
  RCM (blue) with station values in Sweden (black). Left panel is for low-lying
  stations and right panel for stations in mountains.

In mountains there is clear indication of added value in CPM model, especially
in the tail distribution.

- Bottom figure shows the wind speed as a function of altitude in present
  climate in RCM (blue/green) and CPM (red/purple).
:::

:::Chapter{headline="Summary" image="scandinavia.png"}
## How do we benefit from using CPM models in Fenno-Scandinavia?
- More realistic representation of precipitation, especially in the summer
  season (improvements in frequency, duration and intensity)
- Potential added value through systematically different climate change response
  in heavy summer precipitation.
- Some suggestions of different climate change responses between CPM and RCM in
  solid versus liquid precipitation in the Scandes mountains.
- Near-surface wind speeds are also better represented in the Scandes. Ongoing
  work to investigate potential added value in climate scenarios as well as more
  process based perspectives (e.g. orographic down-slope wind-storms).
:::
---
id: 7
trl: medium
category: EUCP science and methodologies
title: Quasi-stationary intense rainstorms spread across Europe under climate change
author: A. Kahraman, et al. Newcastle University.
thumbnail: background.png
---

:::Chapter{headline="Background" image="background.png"}
## The background and method

The ingredients-based approach for forecasting heavy precipitation (Doswell et
al. 1986) suggests that the amount of precipitation in a location is equal to
the product of average rainfall rate and duration. For a convective storm, the
average rainfall rate is a function of specific humidity (available moisture),
vertical velocity (condensation rate of the moisture), and precipitation
efficiency (fraction of condensed water falling to ground). Local duration
depends on many factors including organization, size and movement of the storm
system, but we focus on the (slower) storm motion only, to maximize local
exposure.

Using 2.2km pan-European simulations for the present (1998-2007) and RCP8.5
future (~2100) simulations, we calculate the co-existence of high moisture
(specific humidity ≥ 10 g kg−1 at 850 hPa) and ascent (vertical velocity ≥ 2 m
s−1 ascent at 700 hPa) with 3-hourly intervals, namely "Extreme Precipitation
Potential" (EPP). Neglecting the effects of highly uncertain precipitation
efficiency, this approximates to hourly precipitation extremes. Adding another
condition for storm motion (Corfidi Vector ≤3 m s−1), we look for "Slow-moving
Extreme Precipitation Potential" (SEPP).
:::

:::Chapter{headline="Metrics: Present and future" image="epp_sepp.png"}
## Extreme Precipitation Potential and Slow-moving Extreme Precipitation Potential in the present and end-of-century (RCP8.5) climates

All of Europe is prone to intense rainstorms as measured by EPP, but the central
Mediterranean experiences the highest frequency of cases, both currently and in
the future. In contrast, All of Europe is prone to intense rainstorms as
measured by EPP, and these are 7.4x more frequent by the end of century (from
∼24 per 100×100km area in the present climate to ~175 in the future). On the
other hand, the slow-moving ones, i.e. SEPPs are relatively rare in the current
climate but become widespread across the continent by 2100, with a 10.6x figure
(changing from ~0.7 to 7.2 per 100×100km area).

In concert with these increases, the number of events with precipitation ≥100 mm
h−1 (about 20% of EPPs in the current climate) increases threefold, while ≥150
mm h−1 is experienced 4× more frequently, and 200 mm h−1 5.2× more frequently by
2100.
:::

:::Chapter{headline="Monthly distribution" image="distribution.png"}
## Monthly distribution

Both EPP and SEPP peaks in the late summer and early autumn, with land areas
more towards the summer and sea areas to autumn. Changes over the land are more
pronounced compared to the sea. Currently, 52% of EPP cases occur over land; by
2100 this jumps to 61%, with increases in land (sea) EPP frequencies of 8.6×
(6×) respectively. SEPPs occur relatively equally over land (48%) and sea (52%)
in the current climate but an enormous increase in land SEPPs (14.3×), with a
smaller, but still significant, increase in sea SEPPs (7.3×) increases the land
fraction to 65% by 2100.
:::

:::Chapter{headline="Tracking results" image="tracking.png"}
## Tracking results from model data

Running a precipitation-tracking algorithm on re-gridded hourly precipitation
output from the model, precipitation areas with at least one 12 km grid point
with ≥20 mmh−1 are detected and the movement speed distribution of such storms
are analyzed. The frequency of such storms robustly increases in the future.

Similar to SEPPs, slow-moving storm systems (≤3 ms−1) analyzed with this
approach are most frequent in autumn in the current climate, but become much
more frequent, with a frequency increase higher than that of faster-moving
systems, in the future. Furthermore, storm speed distributions in all seasons
become skewed to the left; and are thus slower on average year-round.
:::

:::Chapter{headline="Change in components" image="change.png"}
## Change in components

The large increase in EPPs stems mainly from the moisture ingredient, with a
dramatic increase of cases with q ≥ 10 gkg−1 in all months, while projected w ≥
2 ms−1 cases are higher only from June to September. By 2100, the CPM projects
29× more cases with very moist environments (exceeding the q threshold) than for
the current climate. The average q increases by ∼35%. The moistest environment
is found in August (5.42 g kg−1 for current climate, increasing to 7.26 g kg−1
by 2100), but change in q peaks in November with a 51% increase. By 2100, the
CPM projects 29× more cases with very moist environments (exceeding the q
threshold) than for the current climate.

By 2100, there is higher frequency of slow-storm environments, except during
February-April, with the annual number of slow-moving storms projected to
increase by 20%. We find the largest increases in August-November, ranging from
31% to 65%, and peaking in September. Since this coincides with the peak
exceedance months of thermodynamic thresholds, the kinematic environment of
these high accumulation storms also contributes to, and enhances, the extreme
precipitation rate for a given locality, resulting in an almost 11-fold increase
in SEPPs compared to a 7.4-fold increase in EPPs by 2100.
:::

:::Chapter{headline="Storm motion" image="precipitation_storm.png"}
## Hourly precipitation extremes and storm motion

A comparison of the maximum hourly precipitation around EPP cases versus the
Corfidi Vector magnitude shows a negative relationship for both the current and
future climate simulations. Clearly, the extreme values of maximum accumulated
hourly precipitation occur with slower storm motions. The PDF distribution of
the Corfidi Vector magnitude for all EPPs vs EPPs with 150 mm or higher hourly
precipitation around them also agrees with this statement. Here, the shift of
the peak towards slower storm motion is also evident for the future climate.
:::

:::Chapter{headline="Relative humidity" image="humidity.png"}
## Relative humidity in the mid-troposphere

Although we exclude the precipitation efficiency factor, we have processed an
analysis of mid-tropospheric relative humidity to investigate one of the main
contributing parameters to precipitation efficiency. The average relative
humidity for the whole domain at the 500 hPa level is seasonally variable, with
lower values in the warm season for both current and future climates. The future
simulation shows decreases in relative humidity, with differences compared to
the current climate simulation highest in the spring. The decrease is more
limited in August, September and October, which corresponds to our high season
for EPP and SEPP cases.

For both current and future simulations, extremely high hourly precipitation
amounts occur exclusively within high relative humidity environments; which
validates the role of environmental relative humidity in terms of precipitation
efficiency.
:::

:::Chapter{headline="Conclusion" image="conclusion.png"}
## Conclusion

This study uses a novel method to examine changes to the ingredients of heavy
precipitation. Our results suggest that storms will have higher peak intensity,
longer duration and will be more frequent across the whole of Europe. Current
storms already produce a large number of flash floods, with their potential
impact depending on land use, terrain slope, drainage, and other factors. SEPP
increases would significantly increase this flash flood potential, as an MCS
would be more likely to "stagnate" on a locality, exposing it to extreme
precipitation of longer duration.

Additionally, storm system movement is associated with upper level winds and,
hence, large-scale dynamics. This may increase fluvial flood risk through
consecutive events in one favorable synoptic setting spanning days or more;
blocking is also favored by large-scale meandering patterns associated with
slower flows.

Understanding the underlying ingredients for heavy precipitation change is
crucial from an impacts perspective, helping to discriminate controlling
factors, which have wider applicability beyond those for a single accumulation
period, and to identify the reliability of projected changes. This suggests that
future studies should focus on precipitation accumulations over space and time.

### Reference

Original paper: [Quasi-Stationary Intense Rainstorms Spread Across Europe Under
Climate Change](https://doi.org/10.1029/2020GL092361).
:::
---
id: 12
trl: medium
category: EUCP data and products
title: Representation and identification of 3 historic "Heavy Precipitation Events"
author: S. Muller, et al., International Centre for Theoretical Physics (ICTP)
thumbnail: "intro.png"
---

:::Chapter{headline="Introduction" image="intro.png"}
## Introduction
Thunderstorms are among the most fascinating weather phenomena. They are
characterized by lightning and thunder, cause strong wind gusts and typically
they bring intense rain. In this way the most severe thunderstorms can be
perceived as heavy precipitation events (HPEs). We here will present a
methodology that, upon objective criteria, identifies HPEs in
convection-permitting observations and numerical models.

Flash floods rank among the most severe hazards caused by weather extremes.
Typically the precipitation events have a strong convective character, are
relatively small in spatial scale and short-lived. Thus it is challenging to
accurately forecast them and to represent them in climate models with coarse
grid spacing.
:::

:::Chapter{headline="Motivation" image="motivation.png"}
## Motivation
Regional Climate Models integrated on convection permitting grids (dx<=4km;
cpRCMs) improve the representation of precipitation extremes, as shown by
[Pichelli et al.
(2021)](https://link.springer.com/article/10.1007/s00382-021-05657-4), e.g. in
their Figure 3. However, information about the properties of the events, causing
these extremes, is not yielded.
:::

:::Chapter{headline="Tracking algorithms" image="algorithms.png"}
## Tracking algorithms
Transform the problem from the Eulerian into the Lagrangian frame of reference
and provide valuable information about the heavy precipitation events
themselves. The schematic shows the Intensity in precipitation on the y-axis and
the scale of events on the x-axis. In this way we identify here convective
phenomena sorted by their typical scales and intensity.

Trackers in principle work the following: first a threshold is imposed to
precipition rate (prTH), second grid cells are clustered in space and time to
form objects, and third all objects that are smaller than a minimal volume
(minvol) are dismissed.

We here apply a well-established tracker (MET MODE-TD, see
https://dtcenter.org/met-online-tutorial-metv8-0/mode-time-domain) on a
composite of high-resolution observations and two cpRCMs. With this we both
demonstrate the capability of convection-permitting resolutions and the
usefulness of the tracker.
:::

:::Chapter{headline="Composite" image="composite.html"}
## Interactive composit map

We identify from a composite of convection-permitting observational datasets,
covering France, Italy, Switzerland and Germany and the years 2001 to 2009,
heavy precipitation events using the MODE-MTD tracker.

The setup is given by a precipitation threshold of 5 mm/hour, a minimum volume
of 100 grid cells and a smoothing in space by 3*3 grid cells.

In this interactive map one can zoom and find historic events, and look into
their properties, which show up when hovering over a circle, which represents a
HPE.

Properties given are: mean area, duration, distance traveled, propagation speed,
total precipitation, mean and maximum precipitation rate, and several
percentiles of precipitation rate. The diameter of each circle is proportional
to the mean area and the color represents the intensity, by the 90th percentile
of precipitation rate.
:::

:::Chapter{headline="Carrara, Italy" image="carrara.png"}
## Carrara, Italy, in September 2003
This HP and flash flood event is summarized on Italian Wikipedia:
https://it.wikipedia.org/wiki/Alluvione_di_Carrara_del_2003. The region was
affected by torrential rain: "in about two and a half hours 200 mm of rain was
falling". It may as well be described as a landfalling event, hot and humid air
from the Mediterranean sea, the tense south-western flow and the orography of
upper Tuscany, were among the triggers that played a major role.

The tracker identifies the event in the observations and attributes fitting
properties: It is short-lived, of small scale and intense.

Also in a cpRCM such events may be well-represented, as indicated.
:::

:::Chapter{headline="Gard, France" image="gard.png"}
## Gard, France, in September 2002
This event is documented in detail by [Delrieu et al.,2005](https://doi.org/10.1175/JHM-400.1).
A landfalling mesoscale convective system remained stationary for almost two
days, trapped by the orography. It brough horrendous amounts of rainfall and
caused wide-spread flooding.

It is well-represented in the observations on the left. Our tracker identifies
it and describes it correctly with a duration of 32 hours, a high intensity and
a large scale and total precipitation amount. Also in a cpRCM such events are
presented and can be identified.
:::

:::Chapter{headline="Western Germany" image="germany.png"}
## Western Germany, 2016/05/26-2016/06/09
A stationary synoptic weather pattern, characterized by a ridge stretching from
Great Britain to Scandinavia, caused a blocking situation. The pressure gradient
and the resulting wind speed in the lower troposphere were very weak.
Consequently, thunderstorms were almost stationary, resulting in large
precipitation accumulations in local areas.

The tracker identifies a multitude of precipitation events in the observations.
:::

:::Chapter{headline="Conclusions" image="conclusion.png"}
## Conclusions
Through the application of trackers valuable information about the events
causing hazards like flash floods can be obtained. Here statistics obtained from
9 year-periods of an observational dataset, and two cpRCM ensembles, one driven
by reanalysis ("evaluation") and one driven by CMIP5-GCMs are presented. On left
we find properties describing the scale of the HPEs and on the right properties
that describe their intensity. In this way models can be evluated against
observations.

Eventually the same can be done for climate projections, promising valuable
information about how HPEs are going to change in a future warmer climate.

Models and Observations provided through the project CORDEX-FPSCONV:
https://www.hymex.org/cordexfps-convection/wiki/doku.php?id=Home

References:

[Pichelli, E., Coppola, E., Sobolowski, S., Ban, N., Giorgi, F., Stocchi, P.,
... & Vergara-Temprado, J. (2021). The first multi-model ensemble of regional
climate simulations at kilometer-scale resolution part 2: historical and future
simulations of precipitation. Climate Dynamics, 56(11),
3581-3602](https://link.springer.com/article/10.1007/s00382-021-05657-4).

[Delrieu, G., Nicol, J., Yates, E., Kirstetter, P. E., Creutin, J. D., Anquetin,
S., ... & Wobrock, W. (2005). The catastrophic flash-flood event of 8–9
September 2002 in the Gard Region, France: a first case study for the
Cévennes–Vivarais Mediterranean Hydrometeorological Observatory. Journal of
hydrometeorology, 6(1), 34-52](https://doi.org/10.1175/JHM-400.1).

[Piper, D., Kunz, M., Ehmele, F., Mohr, S., Mühr, B., Kron, A., & Daniell, J.
(2016). Exceptional sequence of severe thunderstorms and related flash floods in
May and June 2016 in Germany–Part 1: Meteorological background. Natural Hazards
and Earth System Sciences, 16(12),
2835-2850](https://doi.org/10.5194/nhess-16-2835-2016).

Images of flash floods:

- [Gard](https://www.taimsalu.com/uzes/p-uzes/uzes-out/slides/pdg_flood-714.jpg)
- [Carrara](https://it.wikipedia.org/wiki/Alluvione_di_Carrara_del_2003#/media/File:Carrione.jpg)
:::
---
id: 11
trl: medium
category: Application of EUCP innovations
title: Attribution of a small scale, heavy flash flood event to climate change
author: D. Matte et al., University of Copenhagen
thumbnail: radar.gif
---

:::Chapter{headline="The Copenhagen flood of 2011" image="radar.gif"}
## The Copenhagen flood of 2011
In the evening on July 2, 2011 a severe cloudburst occurred over Copenhagen,
Denmark. Between 90 and 135 mm of precipitation was recorded in less than 2
hours, which flooded the city and caused hundreds of millions of Euros in
insured damages.

The animation depicts the evolution of the radar reflectivity of the event.

![picture of flash flood event](stories/flood/Istedgade_skybrud_2011-07-02.jpg)

The image above is taken from
[Wikimedia](https://commons.m.wikimedia.org/wiki/File:Istedgade_skybrud_2011-07-02.jpg#mw-jump-to-license)
where it is shared under a
[CC-BY-SA](https://creativecommons.org/licenses/by-sa/2.0/deed.en).
:::

:::Chapter{headline="Simulating intense convective storms" image="precipitation.png"}
## Simulating intense convective storms
Simulating an intense convective storm is very challenging due, among other
things, to a very high sensitivity to initial conditions and a general
difficulty of the models to resolve small-scale interactions. Here, we use a
sophisticated ensemble approach with the HARMONIE-AROME limited-area numerical
weather prediction model to represent the event.

The figure shows the hourly precipitation at 16 UTC for all ensemble 13 members
(a,d-o), a calculated risk of exceeding a threshold of 60 mm/h (b) and the
associated ensemble mean (c). The red star indicates the location of Copenhagen.
:::

:::Chapter{headline="Risk assessment" image="probability.png"}
## Risk assessment
Our approach indicates that we have been able to catch the main characteristics
of the storm and hence assessing a related risk of exceeding certain
precipitation intensity thresholds. The next objective is to better understand
the influence on the event from due to overall global warming  and projected
further future warming levels.

Using an adapted pseudo-global warming approach, we investigate the related risk
of exceedance of the same event using conditions mimicking pre-industrial
conditions, But also to analyse the related risk of the event using conditions
mimicking different future warming levels.

The figure shows maps of the calculated risk of exceeding a threshold of 60 mm/h
between 15 UTC and 20 UTC for all warming levels (-1°C, control, +1°C, +2°C and
+3°C respectively).
:::

:::Chapter{headline="Take home messages" image="closure.png"}
## Take home messages
- The Harmonie CP NWP model can capture the nature of the extreme Copenhagen
  event of July 2, 2011
- While it is not necessarily possible to simulate the Copenhagen extreme event
  with perfection, it is possible to capture the possible severity and a
  reasonably accurate location when taking an ensemble approach.
- The risk of achieving precipitation amounts such as those observed has
  increased due to global warming.
- The risk of achieving even more intense precipitation rates over the affected
  region increases with further warming.
:::
---
id: 0
trl: high
category: Application of EUCP innovations
title: Multi-year prediction of drought and heat stress in the wheat sector
author: B. Solaraju-Murali, et al. Barcelona Supercomputing Center.
thumbnail: "img_1.png"
---

:::Chapter{headline="Introduction and Aim" image="img_1.png"}
## Introduction and Aim

Wheat sector is heavily influenced by recent alterations in the frequency and
severity of extreme climate events. There is a growing need of climate
information for effective planning and adaptive actions to the near-term climate
variability and change.

The aim is to assess the forecast quality in predicting the evolution of
multi-year drought and heat stress conditions by using user-relevant
agro-climatic index over global wheat growing region. The considered
agro-climatic index includes Standardized Potential Evapotranspiration Index
(SPEI)  and Heat Magnitude Day Index (HMDI).

## Wheat harvest months
The harvested month is retrieved from the MIRCA2000 dataset. For each grid box,
the forecast quality evaluation of initialized decadal prediction is carried-out
for the months preceding the harvest of wheat on a global spatial scale and
therefore the presented results in this evaluation will correspond to different
times of the year for different regions.
:::

:::Chapter{headline="Data and Methods" image="img_2.png"}
## Data and Methods

Datasets used in this project are:

- NCAR Decadal Prediction Large Ensemble (DPLE): Forecasts that were run by
  explicitly prescribing the contemporaneous state of the climate system at the
  start of the simulation annually.
- Observational reference: JRA-55 for monthly mean temperature and GPCCv2018 for
  precipitation.

Methodology is described in the figure.
:::

:::Chapter{headline="Forecast quality assessment: Prediction skill" image="img_3.png"}
## Forecast quality assessment of SPEI6 and HMDI3: Prediction skill

The probabilistic skill of initialized decadal forecast in predicting multi-year
averaged SPEI6 and HMDI3 is presented for the months preceding the wheat harvest
on a global spatial scale.

Fair Ranked probability skill score (FRPSS) for tercile events of calibrated
SPEI6 (left) and HMDI3  (right) averaged over forecast years 1-5 using decadal
prediction, with respect to the observational climatology during the winter
wheat harvest months for the period 1961-2014.
:::

:::Chapter{headline="Forecast quality assessment: Reliability" image="img_4.png"}
## Forecast quality assessment of SPEI6 and HMDI3: Reliability

The reliability of initialized decadal forecast in predicting multi-year
averaged SPEI6 and HMDI3 over global wheat harvesting region.

This figure shows reliability diagrams (lines) for probabilistic categorical
forecasts (tercile events) of unadjusted and calibrated SPEI6 (top) and HMDI3
(bottom) estimate over Global wheat growing region. The frequency of occurrence
of each bin for each category is shown to the right of the panels.
:::

:::Chapter{headline="Decadal prediction: Map" image="img_5.png"}
## Decadal prediction for wheat sector

An illustration of a potential climate service product for the wheat sector
based on categorical events.

This figure shows high probability of dry conditions (below-normal category of
SPEI6) and heat stress events (above-normal category of HMDI3) over most of the
wheat growing regions for the forecast years 2015-2019.

Multi-year probabilistic calibrated forecast (a, b) and observed (c, d) most
likely tercile category of SPEI6 (left) and HMDI3 (right) for the year 2014-2018
over the wheat harvesting regions that presents positive FRPSS values.
:::

:::Chapter{headline="Decadal prediction: Timeseries" image="img_6.png"}
## Decadal prediction for wheat sector

Here is another illustration of a potential climate service product for the
wheat sector based on categorical events.

Jesi is one of the identified regions by the stakeholders in the EU
H2020-MED-GOLD project to be among the leading regions used for the production
of durum wheat in Italy.

This figure shows time series of multi-year averaged SPEI6 (top) and HMDI3
(bottom) over Jesi, Italy (43.82° N, 13.75° E). The small (large) gray dots
correspond to the ensemble members of the hindcasts (ensemble members mean) for
each start date. The black dots correspond to the observed values. The brown
(blue) and green (red) horizontal lines show its lower and upper observed
terciles of SPEI6 (HMDI3), respectively.
:::

:::Chapter{headline="General conclusions" image="img_7.png"}
## General conclusions

Our findings lead to the general conclusions:

- Initialized decadal predictions hold promise in predicting drought and heat
  stress conditions on a multi-annual timescale.
- Positive probabilistic skill is found for SPEI6 and HMDI3 index, particularly
  over Europe and the Middle Eastern countries during the wheat harvest months.
- Calibration of the indices improves the skill over several regions, which
  presents high negative values.
- Using Attributes (or Reliability) Diagrams, it is evident that Calibration (or
  bias adjustment) of the index is essential for providing trustworthy and
  robust forecasts to the interested end users in the agricultural sector.

## More information

Detailed information is available in the [reference
paper](https://www.nature.com/articles/s41612-021-00189-4).
:::
---
id: 2
trl: medium
category: EUCP science and methodologies
title: Physical storylines of future European drought events
author: K. van der Wiel, et al. Royal Netherlands Meteorological Institute (KNMI).
thumbnail: img_1.gif
---

:::Chapter{headline="Advantages of large ensemble simulations" image="img_1.gif"}
## (Some) advantages of large ensemble simulations

- **Large ensembles** are many realisations of the same climate model
  experiment, the many realisations capture **internal** variability within the
  climate system (Deser et al., 2020, NCC).
- Large ensembles thus contain a large diversity of weather events, we can use
  this to **select events** and **investigate drivers or changes** (Van der Wiel
  et al., 2020, ERL).
:::

:::Chapter{headline="Advantages of storylines" image="img_2.jpg"}
## (Some) advantages of storylines

- People are very well trained in listening to and understanding stories. Just
  think: "Once upon a time…".
- **Storylines** improve risk awareness and decision-making processes. This is
  due to their **connection with observed events** and user-focus (Hazeleger et
  al., 2015, NCC; Shepherd et al., 2018, CC).
- In this work we aim to combine the positives of large ensembles and storylines
  to answer societally relevant questions.
:::

:::Chapter{headline="Proposed methodology" image="img_3.jpg"}
## Proposed methodology

Two step process:

1. Select simulated events like the observed event: 'analogues'
1. Compare the analogues in present and future climates

To verify the quality of the results, we have to consider:

- How to define the observed event, how to select analogues?
- Are there analogues that resemble the observed event?
- Does this provide robust and reliable climate change information?

In the remainder of this storyboard we provide a **proof-of-concept**, creating
storylines of **future European drought events like 2018**.
:::

:::Chapter{headline="Selection of simulated analogues" image="img_4.jpg"}
## Selection of simulated analogues

- We define metrics that describes (part of) the event of interest, this is done
  by choosing a variable and its spatial and temporal properties.
- Using this event definition, we select simulated events from the large
  ensemble dataset by means of maximum similarity.
- In the case study we used three event metrics, metric 1 is shown here. Twenty
  simulated events were selected, each of these with very high August-October
  precipitation deficit as the observed event in 2018. We took the mean over
  these twenty events, and used this composite 'analogue' in the remainder of
  the study.
:::

:::Chapter{headline="Investigation of climate change and analogues" image="img_5.jpg"}
## Investigation of climate change and analogues

- We compute analogues using three large ensemble datasets, describing
  present-day climate and two future warmer climates (pre-industrial + 2C
  warming and +3C warming).
- These analogues can then be compared to identify the influence of climate
  change on the event of interest. The analogues show what the event might look
  like in a warmer future climate.
- In our case study, we find that **future droughts are more severe**. Lower
  precipitation and higher evapotranspiration lead to **higher future peak
  precipitation deficits**. The exact timing and magnitude of changes are metric
  dependent.
:::

:::Chapter{headline="Conclusions" image="img_6.jpg"}
## Conclusions

- The large ensemble storyline method works. It allows to create analogues of
  observed events in future warmer climates.
- The ensemble data set has to be quite bit, this is necessary to find the right
  balance between analogue composite size (to decrease influence of natural
  variability, noise) and analogue return period (to sample an appropriate
  extreme event).
- More information in published paper [Van der Wiel, K., G. Lenderink, H. de
  Vries: Physical storylines of future European drought events like 2018 based
  on ensemble climate modelling, Weather and Climate Extremes, 33, pp.
  100350](http://doi.org/10.1016/j.wace.2021.100350).
- The dataset used in this study is available through the following link:
  [https://doi.org/10.5281/zenodo.5083159](https://zenodo.org/record/5083160#.YZ-1Db3MKcY)
:::
---
id: 1
trl: low
category: EUCP data and products
title: The EUCP Caribbean runs
author: H. de Vries, et al. Royal Netherlands Meteorological Institute (KNMI).
thumbnail: carribean.png
---

:::chapter{headline=Introduction image=carribean.png}
## Caribbean

Many of the regional climate simulations carried out as part of EUCP were
situated over Europe. On explicit request by the EU, EUCP also conducted a
limited number of simulations with very high resolution regional climate models
for a number outer-European domains. One of the domains is the Caribbean, known
for its sunshine and sandy beaches, but also for devastating tropical cyclones,
monster waves and coastal erosion.

## Convection permitting (CP) simulations

The research question addressed here is what very high resolution CPM
simulations could add to the story of climate change in the region, with a
special interest in the wind and precipitation characteristics of tropical
cyclones. Additionally, it could perhaps be used to apply methodology derived
primarily for the European domains, to these other domains. Two other sets of
outer-European CPM simulations were created within EUCP: one for Madeira, and
one for the island of La Réunion. These are not discussed in this presentation.
:::

:::Chapter{headline="Tropical cyclones variability", image="TCvariability.png"}
## Tropical cyclones variability is huge

As the previous slide already showed, there is huge spatial and temporal
variability in the tracks of tropical cyclones. There is a large interannual
variability with some years being relatively quiet and others very active.

Because of finite resources, only a few hurricane seasons could be simulated.
Therefore it made sense to attempt a pseudo-global warming approach, in which
observed hurricane seasons are simulated twice: once with observed (ERA5)
forcing at lateral boundaries, and a second time in which the same seasons are
'futurized' by adding seasonally varying delta-fields obtained from
climate-change experiments (subset of CMIP5).
:::

:::Chapter{headline="Literature" image="literature.png"}
## Existing literature

Existing literature on the subject (CPM simulations of the domain) is not very
abundant yet. The paper by Gutmann et al (2018) found that a Pseudo-Global
Warming (PGW) approach gave satisfying results, with TC having faster maximum
winds, slower storm translation, lower central pressures and higher
precipitation rates. An important remark was that despite the PGW approach,
there was surprisingly large variability between the storms, necessitating the
creation of a large ensemble.

**Gutmann et al. (2018)**: Changes in Hurricanes from a 13-Yr
Convection-Permitting Pseudo-Global Warming Simulation
- **ERA-Interim**. Period 2001-2013 (32 storms).
- **PGW**: Changes from CMIP5
- CPM: WRF@4km. 5440 km (east-west) x 4064km (north-south)

**Key results**: on average faster maximum winds, slower storm translation
speeds, lower central pressures, higher precipitation rates.

Important remark: changes varied substantially between individual storms, and
not all storms were simulated in both simulations!
:::

:::Chapter{headline="Models and strategies" image="strategy.png"}
## Contributing groups/models
Several groups signed up to contribute with CPM (and RCM) simulations.

**CPM**
- KNMI - HCLIM38-AROME (nested in RACMO)
- SMHI - HCLIM38-AROME (direct)
- ICTP - RegCM4-7 (direct)
- GERICS - REMO-NH (direct)

**RCM**
- KNMI - RACMO (1979-2020)
- DMI - HCLIM38-AROME 12km(!)

## Simulation strategy
Simulation strategy was to come up with a top-10 list of hurricane seasons and
simulate these using ERA5 as boundary and initial conditions (REF). A second set
of simulations (PGW) was then created, by adding CMIP5-derived delta fields to
the ERA5 boundary conditions. In this way, entire seasons were 'futurized'.

Note that within the domain, the models generate their own internal dynamics so
it is not guaranteed that each REF and PGW run gives the same number of tropical
cyclones. KNMI ran a third set of simulations using a uniform temperature
increase.
:::

:::Chapter{headline="Top-10 list of seasons" image="top10.png"}
## Top-10 list
The top-10 list of seasons was created on the basis of the IBTrACS database and
a couple of criteria (see the title of the plot). The approximate internal
domain is shown, as well as the named hurricanes occurring in the top-10
seasons.
:::

:::Chapter{headline="First results of Irma" image="irma.gif"}
## First results
First result is of Irma, the big tropical cyclone of 2017. Shown is the
simulation of the hourly windfield. The black line is the observed track from
the IBTrACS database. As you can see the TC is well resolved (clear eye
formation) and follows the observed track rather accurately, making it quite
suitable for impact studies in e.g. WP4. At the end of the simulation period,
the next tropical cyclone (Jose) is just appearing.
:::

:::Chapter{headline="Poor man's tracking" image="tracking.png"}
## Poor man's tracking
As a second example we count the number of hours with hurricane-force windspeeds
in each of the models in the 2017 hurricane season. This could be considered as
a poor man's tracking approach. Although most models resolve the TC from 2017,
there is still some model spread (remember the lateral forcing is all
identical).

Note that the models labelled SMHI and KNMI are identical (HCLIM-AROME), but
that KNMI uses RACMO as the intermediate model, whereas SMHI uses a direct
approach. [SMHI and DMI are also the same model, but their resolutions differ].
Track coherence seems weakest in GERICS (REMO-NH model).
:::

:::Chapter{headline="PGW perturbations" image="trackingPGW.png"}
## Poor man's tracking with PGW perturbations added
In the second experiment, the same seasons were simulated with PGW perturbations
added. Again we see considerable differences, although some of the major
hurricanes still feature on the map. But clearly it is not guaranteed that each
TC get more intense under PGW. To the contrary, most of the TC shown here seem
to weaken (except in the simulations by GERICS (REMO-NH). Also the tracks are
not identical, implying that studying impact changes is not straightforward,
unless we aggregate / pool results. A similar pattern is shown when we use all
years. Many storms weaken, some seem to intensify, and most display small
changes in track location.
:::

:::Chapter{headline="Precipitation" image="precipitation.png"}
## Precipitation
A first glance at the changes in precipitation can be obtained from histograms
of the hourly spatial maximum of the precipitation field, aggregated over all
hours and all years. By looking at the field-maximum we thus discard information
on the size of events. Furthermore, only sea area is considered here.

There are substantial differences between the models yet all show unmistakably a
tendency towards higher maxima in future. A more complete analysis of the
distribution, that includes the size of events is in preparation.
:::

:::Chapter{headline="Conclusions" image="summary.png"}
## Summary and conclusions

EUCP Caribbean CPM and RCM simulations: ERA5 + PGW, focus on Hurricane season.
Top-10 seasons according to criteria. Includes 2017 (Irma). Data will be used by
WP4 for eg., storm surge modelling.

Preliminary results: wind speed and precipitation. Robust increase in
precipitation maxima, wind speed change less obvious. Poor man’s tracking: Not
all tropical cyclones reproduced in REF and PGW: Optimal conditions in
present-day climate not necessarily more favorable for development in future
climate. Further analysis necessary.

Still to do:
- More extended analysis, including pooling of results.
- Comparison to observations
- Impact modelling by WP4: storm surges
- Discussion of relevance of PGW versus GCM-forced.

## Reference

[Gutmann, E. D. et al. Changes in hurricanes from a 13-Yr convection-permitting
pseudo- global warming simulation. J. Clim. 31, 3643–3657
(2018).](https://journals.ametsoc.org/view/journals/clim/31/9/jcli-d-17-0391.1.xml?tab_body=pdf)

[Lenderink, G. et al. Systematic increases in the thermodynamic response of
hourly precipitation extremes in an idealized warming experiment with a
convection-permitting climate model. Environ. Res. Lett. 14,
(2019).](https://iopscience.iop.org/article/10.1088/1748-9326/ab214a/meta)

[Liu, C. et al. Continental-scale convection-permitting modeling of the current
and future climate of North America. Clim. Dyn. 49, 71–95
(2017).](https://link.springer.com/article/10.1007/s00382-016-3327-9)
:::
---
id: 15
trl: high
category: Application of EUCP innovations
title: Assessment and attribution of the changes in wind energy in Europe
author: Shuang Yu et al., IPSL
thumbnail: "intro.png"
---
:::Chapter{headline="Introduction" image="intro.png"}
## Introduction

The decline of wind speed in some regions of the world is expected to decrease
wind power production.

As the occurrence of wind drought receives increasing attention, two questions
should be discussed.

1. Will the decrease in wind energy resources continue in the future?
2. Are the current "wind energy drought" caused by climate change?

**Aim**: Systematic attribution work of wind power production in Europe.
:::

:::Chapter{headline="Data and Methods" image="data_methods.png"}
## Data and Methods

Datasets used in this project are:

- ERA5 forecast dataset:
  - 10m wind speed (w10), w100, 2m temperature (t2m), 850hPa temperature (850t)
- Model simulation:
  - 10 models from EURO-CORDEX ensemble, variables including: w10, t2m, 850t.
    Machine Learning method was used to get the w100 from 10 model ensembles.
- Mast stations


**Load factor** - the ratio of average electrical power to the rated electrical
power, is used to characterize wind energy.
:::

:::Chapter{headline="Wind drought in the Whole Europe" image="europe.png"}
## Wind drought in the Whole Europe

The figure shows the changes in return periods and return values for 10 model
ensembles in three different periods.

-> Previous:1971-2000

-> Current: 2000-2029

-> Future: 2041-2069

**For Whole Europe**, the low wind energy winter becomes more frequent in the
recent year and the increase will be continue in the future. During 2041-2070,
the risk ratio is about 1.52, which implies that the probability that the load
factor is lower than the threshold would be increased in the future climate, and
it would be more than 1.52 times as much as that in the current climate.

**In Britanny**, changes of wind energy drought are not significant during three
different periods.
:::

:::Chapter{headline="Wind energy drought in different countries" image="othercountries.png"}
## Wind energy drought in different countries

The figure shows the risk ratios in each country during different climate
periods.

When the risk ratio >1, the low wind energy events are more frequent compared
with the earlier climate period.

Most countries are experiencing more frequent low wind energy events compared
with the previous period. In the future, low wind energy events will become more
common in most countries, especially in northern countries where wind energy
drought will occur more frequently.
:::

:::Chapter{headline="Wind energy drought in a warming world" image="warmingworld.png"}
## Wind energy drought in a warming world

The figure shows the return periods for model ensemble under 1.5, 2, and 3℃
global warming scenarios. The 5% and 95% of the medians show the uncertainty
interval.

The changes of low wind energy in Europe are particularly vulnerable to climate
change, experience adverse effect of global warming in the whole Europe. With
global warming, low wind energy winter in Europe occur more frequently.
:::

:::Chapter{headline="General conclusion" image="conclusion.png"}
## General conclusion

- Whole Europe:

  Low wind energy winter in Europe is expected to become more common in a
  warming world. The return period of 11 years return period events have been
  shortened to about 10 years in the whole Europe and it would be shortened to
  about 7 years during 2041-2070 under RCP 8.5 scenarios.

- Regional differences:

  Low wind energy winter are becoming more common in most countries from the
  past to the current, and it will become more frequent in the future,
  particularly in northern countries, which would have a negative impact on the
  development of wind power generation in these regions.

:::
# STORE

**This directory is not required, you can delete it if you don't want to use it.**

This directory contains your Vuex Store files.
Vuex Store option is implemented in the Nuxt.js framework.

Creating a file in this directory automatically activates the option in the framework.

More information about the usage of this directory in [the documentation](https://nuxtjs.org/guide/vuex-store).
