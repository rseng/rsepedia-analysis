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
reported to the community leaders responsible for enforcement:
Clemens Schmid (schmid@shh.mpg.de).
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

[![Project Status: Inactive – The project has reached a stable, usable
state but is no longer being actively developed; support/maintenance
will be provided as time
allows.](https://www.repostatus.org/badges/latest/inactive.svg)](https://www.repostatus.org/#inactive)
![GitHub R package
version](https://img.shields.io/github/r-package/v/nevrome/bleiglas)
[![R-CMD-check](https://github.com/nevrome/bleiglas/actions/workflows/check-release.yaml/badge.svg)](https://github.com/nevrome/bleiglas/actions/workflows/check-release.yaml)
[![Coverage
Status](https://img.shields.io/codecov/c/github/nevrome/bleiglas/master.svg)](https://codecov.io/github/nevrome/bleiglas?branch=master)
[![license](https://img.shields.io/github/license/nevrome/bleiglas)](https://www.r-project.org/Licenses/MIT)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03092/status.svg)](https://doi.org/10.21105/joss.03092)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# bleiglas

bleiglas is an R package that employs
[Voro++](http://math.lbl.gov/voro++/) for the calculation of three
dimensional Voronoi diagrams from input point clouds. This is a special
form of tessellation where each polygon is defined as the area closest
to one particular seed point. Voronoi diagrams have useful applications
in - among others - astronomy, material science or geography and
bleiglas provides functions to make 3D tessellation more readily
available as a mean for data visualisation and interpolation. It can be
used for any 3D point cloud, but the output is optimized for
spatiotemporal applications in archaeology.

1.  This README (see Quickstart guide below) describes a basic workflow
    with code and explains some of my thought process when writing this
    package.
2.  A [JOSS paper](https://doi.org/10.21105/joss.03092) gives some
    background, introduces the core functions from a more technical
    point of view and presents an example application.
3.  A (rather technical) vignette presents all the code necessary to
    reproduce the “real world” example application in said JOSS paper.
    When bleiglas is installed you can open the vignette in R with
    `vignette("bleiglas_case_study")`.

If you have questions beyond this documentation feel free to open an
[issue](https://github.com/nevrome/bleiglas/issues) here on Github.
Please also see our [contributing guide](CONTRIBUTING.md).

## Installation

You can install bleiglas from github

``` r
if(!require('remotes')) install.packages('remotes')
remotes::install_github("nevrome/bleiglas", build_vignettes = TRUE)
```

For the main function `tessellate` you also have to [install the Voro++
software](http://math.lbl.gov/voro++/download/). The package is already
available in all major Linux software repositories (on Debian/Ubuntu you
can simply run `sudo apt-get install voro++`.). MacOS users should be
able to install it via homebrew (`brew install voro++`).

## Quickstart

For this quickstart, we assume you have packages `tidyverse`, `sf`,
`rgeos` (which in turn requires the Unix package `geos`) and `c14bazAAR`
installed.

#### Getting some data

I decided to use Dirk Seidenstickers [*Archives des datations
radiocarbone d’Afrique
centrale*](https://github.com/dirkseidensticker/aDRAC) dataset for this
purpose. It includes radiocarbon datings from Central Africa that
combine spatial (x & y) and temporal (z) position with some meta
information.

<details>
<summary>
Click here for the data preparation steps
</summary>
<p>

I selected dates from Cameroon between 1000 and 3000 uncalibrated BP and
projected them into a worldwide cylindrical reference system (epsg
[4088](https://epsg.io/4088)). As Cameroon is close to the equator this
projection should represent distances, angles and areas sufficiently
correct for this example exercise. As a minor pre-processing step, I
here also remove samples with equal position in all three dimensions for
the tessellation.

``` r
# download raw data with the data access package c14bazAAR
# c14bazAAR can be installed with
# install.packages("c14bazAAR", repos = c(ropensci = "https://ropensci.r-universe.dev"))
c14_cmr <- c14bazAAR::get_c14data("adrac") %>% 
  # filter data
  dplyr::filter(!is.na(lat) & !is.na(lon), c14age > 1000, c14age < 3000, country == "CMR") 
```

    ##   |                                                          |                                                  |   0%  |                                                          |++++++++++++++++++++++++++++++++++++++++++++++++++|  99%  |                                                          |++++++++++++++++++++++++++++++++++++++++++++++++++| 100%

``` r
# remove doubles
c14_cmr_unique <- c14_cmr %>%
  dplyr::mutate(
    rounded_coords_lat = round(lat, 3),
    rounded_coords_lon = round(lon, 3)
  ) %>%
  dplyr::group_by(rounded_coords_lat, rounded_coords_lon, c14age) %>%
  dplyr::filter(dplyr::row_number() == 1) %>%
  dplyr::ungroup()

# transform coordinates
coords <- data.frame(c14_cmr_unique$lon, c14_cmr_unique$lat) %>% 
  sf::st_as_sf(coords = c(1, 2), crs = 4326) %>% 
  sf::st_transform(crs = 4088) %>% 
  sf::st_coordinates()

# create active dataset
c14 <- c14_cmr_unique %>% 
  dplyr::transmute(
    id = seq_len(nrow(.)),
    x = coords[,1], 
    y = coords[,2], 
    z = c14age,
    period = period
)
```

</p>
</details>
<details>
<summary>
Data: <b>c14</b>
</summary>
<p>

``` r
c14 
```

    ## # A tibble: 393 × 5
    ##       id        x       y     z period
    ##    <int>    <dbl>   <dbl> <int> <chr> 
    ##  1     1 1284303. 450340.  1920 EIA   
    ##  2     2 1101276. 321798.  2340 EIA   
    ##  3     3 1101276. 321798.  2520 LSA   
    ##  4     4 1093159. 264311.  2000 <NA>  
    ##  5     5 1132077. 340034.  1670 <NA>  
    ##  6     6 1101276. 321798.  2200 <NA>  
    ##  7     7 1101276. 321798.  2030 <NA>  
    ##  8     8 1101276. 321798.  1760 EIA   
    ##  9     9 1093159. 264311.  1710 <NA>  
    ## 10    10 1093159. 264311.  1940 <NA>  
    ## # … with 383 more rows

</p>
</details>

#### 3D tessellation

[Tessellation](https://en.wikipedia.org/wiki/Tessellation) means filling
space with polygons so that neither gaps nor overlaps occur. This is an
exciting application for art (e.g. textile art or architecture) and an
interesting challenge for mathematics. As a computational archaeologist
I was already aware of one particular tessellation algorithm that has
quite some relevance for geostatistical analysis like spatial
interpolation: Voronoi tilings that are created with [Delaunay
triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation).
These are tessellations where each polygon covers the space closest to
one of a set of sample points.

<table style="width:100%">
<tr>
<th>
<figure>
<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/6/66/Ceramic_Tile_Tessellations_in_Marrakech.jpg/320px-Ceramic_Tile_Tessellations_in_Marrakech.jpg" height="150" />
<br>
<figcaption>
Islamic mosaic with tile tessellations in Marrakech, Morocco.
<a href="https://en.wikipedia.org/wiki/File:Ceramic_Tile_Tessellations_in_Marrakech.jpg">wiki</a>
</figcaption>
</figure>
</th>
<th>
<figure>
<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/5/56/Delaunay_Voronoi.svg/441px-Delaunay_Voronoi.svg.png" height="150" />
<br>
<figcaption>
Delaunay triangulation and its Voronoi diagram.
<a href="https://commons.wikimedia.org/wiki/File:Delaunay_Voronoi.svg">wiki</a>
</figcaption>
</figure>
</th>
<th>
<figure>
<img src="https://aip.scitation.org/na101/home/literatum/publisher/aip/journals/content/cha/2009/cha.2009.19.issue-4/1.3215722/production/images/medium/1.3215722.figures.f4.gif" height="150" />
<br>
<figcaption>
Output example of Voro++ rendered with POV-Ray.
<a href="http://math.lbl.gov/voro++">math.lbl.gov</a>
</figcaption>
</figure>
</th>
<tr>
</table>

It turns out that Voronoi tessellation can be calculated not just for 2D
surfaces, but also for higher dimensions. The
[Voro++](http://math.lbl.gov/voro++/) software library does exactly this
for 3 dimensions. This makes it useful for spatio-temporal applications.

`bleiglas::tessellate()` is a minimal wrapper function that calls the
Voro++ command line interface (therefore you have to install Voro++ to
use it) for datasets like the one introduced above. We can apply it like
this:

``` r
raw_voro_output <- bleiglas::tessellate(
  c14[, c("id", "x", "y", "z")],
  x_min = min(c14$x) - 150000, x_max = max(c14$x) + 150000, 
  y_min = min(c14$y) - 150000, y_max = max(c14$y) + 150000,
  unit_scaling = c(0.001, 0.001, 1)
)
```

A critical step when using tessellation for spatio-temporal data is a
suitable conversion scale between time- and spatial units. Since 3D
tessellation crucially depends on the concept of a 3D-distance, we need
to make a decision how to combine length- and time-units. Here, for the
purpose of this example, we have 1 kilometre correspond to 1 year. Since
after the coordinate conversion our spatial units are given in meters,
we divide all spatial distances by a factor 1000 to achieve this
correspondence: `unit_scaling = c(0.001, 0.001, 1)`.

I decided to increase the size of the tessellation box by 150 kilometres
to each (spatial) direction to cover the area of Cameroon. Mind that the
scaling factors in `unit_scaling` are also applied to the box size
parameters `x_min`, `x_max`, ….

The output of Voro++ is highly customizable, and structurally complex.
With the `-v` flag, the voro++ CLI interface prints some config info,
which is also the output of `bleiglas::tesselate`:

    Container geometry        : [937.154:1936.57] [63.1609:1506.58] [1010:2990]
    Computational grid size   : 3 by 5 by 6 (estimated from file)
    Filename                  : /tmp/RtmpVZjBW3/file3aeb5f400f38
    Output string             : %i*%P*%t
    Total imported particles  : 392 (4.4 per grid block)
    Total V. cells computed   : 392
    Total container volume    : 2.8563e+09
    Total V. cell volume      : 2.8563e+09

It then produces an output file (`*.vol`) that contains all sorts of
geometry information for the calculated 3D polygons. `tesselate` returns
the content of this file as a character vector with the additionally
attached attribute `unit_scaling`
(`attributes(raw_voro_output)$unit_scaling`), which is just the scaling
vector we put in above.

I focussed on the edges of the polygons and wrote a parser function
`bleiglas::read_polygon_edges()` that can transform the complex Voro++
output for this specific output case to a tidy data.table with six
columns: the coordinates (x, y, z) of the start (a) and end point (b) of
each polygon edge. A data.table is a tabular R data structure very
similar to the standard data.frame. Read more about it
[here](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html).

``` r
polygon_edges <- bleiglas::read_polygon_edges(raw_voro_output)
```

`read_polygon_edges` automatically reverses the rescaling introduced in
`tesselate` with the `unit_scaling` attribute.

<details>
<summary>
Data: <b>polygon\_edges</b>
</summary>
<p>

    ##            x.a     y.a     z.a     x.b    y.b     z.b polygon_id
    ##     1:  937154  374130 1307.99 1201480 392161 1299.80         25
    ##     2: 1289460  241706 1324.42 1201480 392161 1299.80         25
    ##     3: 1212280  387619 1290.18 1201480 392161 1299.80         25
    ##     4: 1190480  335990 1202.59 1233970 377206 1268.57         25
    ##     5: 1352310  233958 1240.81 1233970 377206 1268.57         25
    ##    ---                                                          
    ## 24916: 1683200  887252 2655.00 1645270 892489 2655.00        290
    ## 24917: 1622180  900165 2682.50 1645270 892489 2655.00        290
    ## 24918: 1393030 1012170 2682.50 1622180 900165 2682.50        290
    ## 24919: 1596200  911750 2731.50 1622180 900165 2682.50        290
    ## 24920: 1645270  892489 2655.00 1622180 900165 2682.50        290

</p>
</details>
<details>
<summary>
We can plot these polygon edges (black) together with the input sample
points (red) in 3D.
</summary>
<p>

``` r
rgl::axes3d()
rgl::points3d(c14$x, c14$y, c14$z, color = "red")
rgl::aspect3d(1, 1, 1)
rgl::segments3d(
  x = as.vector(t(polygon_edges[,c(1,4)])),
  y = as.vector(t(polygon_edges[,c(2,5)])),
  z = as.vector(t(polygon_edges[,c(3,6)]))
)
rgl::view3d(userMatrix = view_matrix, zoom = 0.9)
```

</p>
</details>

<img src="README_files/figure-gfm/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

#### Cutting the polygons

This 3D plot, even if rotatable using mouse input, is of rather limited
value since it’s very hard to read. I therefore wrote
`bleiglas::cut_polygons()` that can cut the 3D polygons at different
levels of the z-axis. As the function assumes that x and y represent
geographic coordinates, the cuts produce sets of spatial 2D polygons for
different values of z – in our example different points in time. The
parameter `cuts` takes a numeric vector of cutting points on the z axis.
`bleiglas::cut_polygons()` yields a rather raw format for specifying
polygons. Another function, `bleiglas::cut_polygons_to_sf()`, transforms
it to `sf`. Here `crs` defines the spatial coordinate reference system
of x and y to project the resulting 2D polygons correctly.

``` r
cut_surfaces <- bleiglas::cut_polygons(
  polygon_edges, 
  cuts = c(2500, 2000, 1500)
) %>%
  bleiglas::cut_polygons_to_sf(crs = 4088)
```

<details>
<summary>
Data: <b>cut\_surfaces</b>
</summary>
<p>

    ## Simple feature collection with 76 features and 2 fields
    ## Geometry type: POLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: 937154 ymin: 63160.9 xmax: 1936570 ymax: 1506580
    ## Projected CRS: World Equidistant Cylindrical (Sphere)
    ## First 10 features:
    ##                                 x    z  id
    ## 1  POLYGON ((1195386 319810.5,... 2500   3
    ## 2  POLYGON ((1936570 809055.4,... 2500  31
    ## 3  POLYGON ((1146675 374628.2,... 2500  38
    ## 4  POLYGON ((1215947 365177.1,... 2500  40
    ## 5  POLYGON ((1416056 455852, 1... 2500  69
    ## 6  POLYGON ((1083149 968036.7,... 2500 103
    ## 7  POLYGON ((1936570 315020.3,... 2500 105
    ## 8  POLYGON ((1386575 333838.1,... 2500 135
    ## 9  POLYGON ((1116416 63160.9, ... 2500 144
    ## 10 POLYGON ((1377347 63160.9, ... 2500 185

</p>
</details>
<details>
<summary>
With this data we can plot a matrix of maps that show the cut surfaces.
</summary>
<p>

``` r
cut_surfaces %>%
  ggplot() +
  geom_sf(
    aes(fill = z), 
    color = "white",
    lwd = 0.2
  ) +
  geom_sf_text(aes(label = id)) +
  facet_wrap(~z) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
```

</p>
</details>

<img src="README_files/figure-gfm/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

<details>
<summary>
As all input dates come from Cameroon it makes sense to cut the polygon
surfaces to the outline of this administrative unit.
</summary>
<p>

``` r
cameroon_border <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>% 
  dplyr::filter(name == "Cameroon") %>% 
  sf::st_transform(4088)

cut_surfaces_cropped <- cut_surfaces %>% sf::st_intersection(cameroon_border)
```

``` r
cut_surfaces_cropped %>%
  ggplot() +
  geom_sf(
    aes(fill = z), 
    color = "white",
    lwd = 0.2
  ) +
  facet_wrap(~z) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
```

<p>
</details>

<img src="README_files/figure-gfm/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

<details>
<summary>
Finally, we can also visualise any point-wise information in our input
data as a feature of the tessellation polygons.
</summary>
<p>

``` r
cut_surfaces_material <- cut_surfaces_cropped %>%
  dplyr::left_join(
    c14, by = "id"
  )
```

``` r
cut_surfaces_material %>%
  ggplot() +
  geom_sf(
    aes(fill = period), 
    color = "white",
    lwd = 0.2
  ) +
  facet_wrap(~z.x) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
```

</p>
</details>

<img src="README_files/figure-gfm/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />

This quickstart was a simple primer on how to use this package. If you
think the final use case wasn’t too impressive, take a look at this
analysis of Bronze Age burial types through time, as performed in our
[JOSS
paper](https://github.com/nevrome/bleiglas/blob/master/paper/paper.md)
and the
[vignette](https://github.com/nevrome/bleiglas/blob/master/vignettes/complete_example.Rmd).

<!-- Add JOSS paper figure here? Just a suggestion as a further teaser. It's just beautiful -->

## Citation


    To cite bleiglas in publications use:

      Schmid and Schiffels (2021). bleiglas: An R package for interpolation
      and visualisation of spatiotemporal data with 3D tessellation.
      Journal of Open Source Software, 6(60), 3092,
      https://doi.org/10.21105/joss.03092

    A BibTeX entry for LaTeX users is

      @Article{,
        title = {{bleiglas}: An {R} package for interpolation and visualisation of spatiotemporal data with 3D tessellation},
        author = {Clemens Schmid and Stephan Schiffels},
        journal = {Journal of Open Source Software},
        volume = {6},
        number = {60},
        pages = {3092},
        year = {2021},
        doi = {10.21105/joss.03092},
        url = {https://doi.org/10.21105/joss.03092},
      }
- V 1.0.1: Fixed a small memory allocation bug
- V 1.0.0: Release version# Contributing

We love pull requests from everyone. By participating in this project, you
agree to abide by our [code of conduct](CONDUCT.md).

## Getting Started

* Make sure you have a [GitHub account](https://github.com/signup/free). If you are not familiar with git and GitHub, take a look at <http://happygitwithr.com/> to get started.
* [Submit a post for your issue](https://github.com/nevrome/bleiglas/issues/), assuming one does not already exist.
  * Clearly describe your issue, including steps to reproduce when it is a bug, or some justification for a proposed improvement.
* [Fork](https://github.com/nevrome/bleiglas/#fork-destination-box) the repository on GitHub to make a copy of the repository on your account. Or use this line in your shell terminal:

    `git clone git@github.com:your-username/bleiglas.git`
    
## Making changes

* Edit the files, save often, and make commits of logical units, where each commit indicates one concept
* Follow a good [style guide](http://adv-r.had.co.nz/Style.html).
* Make sure you write [good commit messages](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html).
* Make sure you have added the necessary tests for your code changes.
* Run _all_ the tests using `devtools::check()` to assure nothing else was accidentally broken.
* If you need help or unsure about anything, post an update to [your issue](https://github.com/nevrome/bleiglas/issues/).

## Submitting your changes

Push to your fork and [submit a pull request](https://github.com/nevrome/bleiglas/pulls).

At this point you're waiting on us. We like to at least comment on pull requests
within a few days. We may suggest some changes or improvements or alternatives.
---
title: 'bleiglas: An R package for interpolation and visualisation of spatiotemporal data with 3D tessellation'
tags:
  - R
  - 3D data analysis
  - Tessellation
  - Voronoi diagrams
  - Spatiotemporal analysis
  - Archaeology
authors:
  - name: Clemens Schmid
    orcid: 0000-0003-3448-5715
    affiliation: 1
  - name: Stephan Schiffels
    orcid: 0000-0002-1017-9150
    affiliation: 1
affiliations:
  - name: Department of Archaeogenetics, Max Planck Institute for the Science of Human History, Kahlaische Strasse 10, 07745 Jena, Germany
    index: 1
date: 11 February 2020
bibliography: paper.bib
---

# Background

The open source software library [Voro++](http://math.lbl.gov/voro++) [@Rycroft2009-rp] allows fast calculation of Voronoi diagrams in three dimensions. Voronoi diagrams are a special form of tessellation (i.e., filling space with geometric shapes without gaps or overlaps) where each polygon is defined as the area closest to one particular seed point. Imagine a volume in three dimensional space and an arbitrary distribution of unique points within this volume. Voro++ creates a polygon around each point so that everything within this polygon is closest to the corresponding point and farther away from every other point.

Voronoi tessellation has useful applications in all kinds of scientific contexts spanning astronomy (e.g., @Paranjape2020-sg), material science (e.g., @Tsuru2020-ep), and geography (e.g., @Liu2019-fw). In computational and landscape archaeology, Delaunay triangulation and Voronoi diagrams have also been applied [@Nakoinz2016-bq], but to our knowledge, mostly limited to an entirely spatial 2D perspective. 3D tessellation could be employed here to add a third dimension, most intriguingly a temporal one. This could allow for new methods of spatiotemporal data interpolation, analysis, and visualisation.

# Statement of need

The ``bleiglas`` R package serves as an R interface to the Voro++ command line tool. It adds a number of utility functions for particular data manipulation applications, including but not limited to automatic 2D cutting of the 3D Voro++ output for subsequent mapping and grid sampling for position and value uncertainty mitigation. The relevant workflows are explained below. Although we wrote this package for our own needs in archaeology and archaeogenetics, the code is by no means restricted to data from these fields, just as Voronoi tessellation is a generic, subject-agnostic method with a huge range of use cases.

Voronoi tessellation is implemented in many R packages, perhaps most prominently in the deldir (*Delaunay Triangulation and Dirichlet (Voronoi) Tessellation*) [@Turner2021] and tripack (*Triangulation of Irregularly Spaced Data*) [@Gebhardt2020] packages, which were specifically designed for this application. ggvoronoi (*Voronoi Diagrams and Heatmaps with 'ggplot2'*) [@Garrett2018] and `dismo::voronoi()` [@Hijmans2020] build on deldir. Other implementations such as `sf::st_voronoi()` [@Pebesma2018] and `geos::geos_voronoi_polygons()` [@Dunnington2021] rely on the [GEOS library](https://trac.osgeo.org/geos) [@geos2021]. All of these packages and functions focus on 2D data, and to our knowledge, none offer a toolset to handle 3D Voronoi diagrams comparable to that introduced with ``bleiglas``.

# Core functionality

``bleiglas`` provides the `bleiglas::tessellate()` function, which is a command line utility wrapper for Voro++. It requires Voro++ to be installed locally. `tessellate()` takes input points in the form of a `data.frame` with an integer ID and three numeric coordinate columns. Additional Voro++ [options](http://math.lbl.gov/voro++/doc/cmd.html) can be set with a character argument `options` and only the [output format definition](http://math.lbl.gov/voro++/doc/custom.html) (`-c`) is lifted to an extra character argument `output_definition`. `tessellate()` returns a character vector containing the raw output of Voro++ with one vector element corresponding to one row. Depending on the structure of this raw output, different parsing functions are required to transform it to a useful R object. At the moment, ``bleiglas`` provides one such function: `bleiglas::read_polygon_edges()`. It is configured to read data produced with the Voro++ output format string `%i*%P*%t`, which returns polygon edge coordinates. These are necessary for the default ``bleiglas`` workflow illustrated in the example below. Future versions of the package may include other parsing functions for different pipelines.

The output of `read_polygon_edges()` is (for performance reasons) a `data.table` [@Dowle2019] object that can be used with `bleiglas::cut_polygons()`. This function now shoulders the core task of cutting the 3D, polgyon-filled Voro++ output box into 2D slices. 3D data is notoriously difficult to plot and read. Extracting and visualizing slices is therefore indispensable. The necessary algorithm for each 3D polygon can be summarised as finding the cutting point of each polygon edge line with the requested cutting surface and then defining the convex hull of the cutting points as a result 2D polygon. We implemented the line-segment-plane-intersection operation via Rcpp [@Eddelbuettel2017] in C++ for better performance and used `grDevices::chull()` for the convex hull search. The output of `cut_polygons()` is a list (for each cut surface) of lists (for each polygon) of `data.table`s (3D coordinates for each 2D cutting point). Optionally, and in case of horizontal (z-axis) cuts with spatial coordinates on the x- and y-axis, this output can be transformed to an `sf` [@Pebesma2018] object via `bleiglas::cut_polygons_to_sf()`. This significantly simplifies subsequent map plotting.

The final core function of ``bleiglas`` is `bleiglas::predict_grid()`, which paves the way for more complex applications and data subject to a higher degree of positional uncertainty. It employs the tessellation output to predict values at arbitrary positions by determining in which 3D polygons they are located. Therefore `bleiglas::predict_grid()` should theoretically mimic the outcome of a nearest neighbor search, but is fully implemented with the above mentioned tessellation workflow. The core algorithm `bleiglas::attribute_grid_points_to_polygons()` uses `cut_polygons()` to cut the tessellation volume at the z-axis level for each prediction point. It then checks in which 2D polygon the point is located to attribute it accordingly. This is done with custom C++ code initially developed for our recexcavAAR package [@Schmid2017]. `bleiglas::predict_grid()` can be used to automatically rerun tessellation multiple times for data with uncertain position in one or multiple of the three dimensions. In this case the resulting wiggling of input points and therefore output polygons is recorded with the static prediction grid. Different per-prediction-point observations in the grid can eventually be summarised to calculate mean outcomes and deviation. One concrete archaeological application of this feature is temporal resampling from post-calibration radiocarbon age probability distributions.

A prerequisite for performing tessellation in three dimensions is the normalisation or mapping of length units across three dimensions. If all three dimensions have the same units (as is the case for 3D spatial data), this is not an issue, and tessellation works as expected. However, if dimensions have different units, the outcome and meaning of the tessellation depends crucially on how these units are mapped to each other. This is the case for spatiotemporal data, in which one axis denotes time and the other two axes denote a 2D spatial position. In such cases it is critical to use external information to inform on an appropriate scaling. For example, one might set 1 km to correspond to 1 year, in which case two contemporaneous points 100 km apart are considered "as close as" two points 100 years apart. What scaling to use clearly depends on the dataset and how to query it.

# Example: Burial rite distributions in Bronze Age Europe

One strength of ``bleiglas`` is visualisation of spatiotemporal data. Here we show an example of Bronze Age burial rites as measured on radiocarbon dates from burials in Central, Northern, and Northwestern Europe between 2200 and 800 calBC. Information about source data (taken from the [RADON-B database](https://radon-b.ufg.uni-kiel.de) [@kneiselRadonB2013]), data preparation, and meaning are presented in @Schmid2019. A vignette in ``bleiglas`` (`vignette("bleiglas_case_study")`) contains the complete code to reproduce the following figures.

Bronze Age burials can be classified by two main aspects: inhumation vs. cremation (*burial type*) and flat grave vs. burial mound (*burial construction*). \autoref{fig:plot_map} is a map of burials through time for which we have some information about these variables. Each grave has a position in space (coordinates) and in time (median calibrated radiocarbon age). For \autoref{fig:plot_3D} and \autoref{fig:plot_bleiglas}, we only look at the *burial type* aspect. The burials are distributed in a three dimensional, spatiotemporal space and therefore can be subjected to Voronoi tessellation with Voro++. As detailed above, the outcome depends on the relative scaling of the input dimensions - for this example we choose $1\text{kilometer}=1\text{year}$, informed by some intuition about the range of human movement through time.

For \autoref{fig:plot_bleiglas} we cut these polygons into 2D time slices that can be visualized in a map matrix (*bleiglas plot*). We believe this matrix is a visually appealing and highly informative way to communicate processes derived from 3D point patterns. It conveys both the main trends (here, the general switch from inhumation to cremation from the Middle Bronze Age onwards) as well as how much data is available in certain areas and periods. The latter is especially relevant regarding the derived question which resolution could be expected from a model based on this input data.

For the final \autoref{fig:plot_prediction_grid}, we applied temporal resampling with `bleiglas::predict_grid()` to record observational error caused by radiocarbon dating uncertainty. For reasons of computational performance we kept the number of resampling runs and the spatial resolution in this example low. Nevertheless the advantages of the resampling approach are easily seen: areas and periods with uncertain attribution to one burial type or the other are indicated as such.

![Graves in the research area (rectangular frame) dating between 2200 and 800 calBC as extracted from the \protect\href{radon-b.ufg.uni-kiel.de}{Radon-B database}. The classes of the variable burial type are distinguished by colour; the ones of burial construction by shape. The map projection is EPSG:102013 (Europe Albers Equal Area Conic) and the base layer data is taken from the \protect\href{https://www.naturalearthdata.com}{Natural Earth project}.\label{fig:plot_map}](03_map_plot.jpeg)

![Graves in 3D space defined by two spatial (x and y in km) and one temporal (z in years calBC) dimensions with Voronoi polygons constructed by Voro++. Each red dot represents one grave with known burial type, the fine black lines the edges of the result polygons, and the rectangular wireframe box the research area now in space and time.\label{fig:plot_3D}](05_3D_plot.jpeg)

![*bleiglas plot*. Map matrix of 2D cuts through 3D Voronoi polygons as presented in \autoref{fig:plot_3D}. Each subplot shows one timeslice between 2200 and 800 calBC. As each 2D polygon belongs to one input grave and data density in some areas and time periods is very low, some graves are represented in multiple subplots. Color coding and map background is as in \autoref{fig:plot_map}.\label{fig:plot_bleiglas}](06_bleiglas_plot.jpeg)

![Temporal resampling version of the bleiglas map matrix. 30 resampling runs, and a spatial resolution of 100*100 cells in the research area shape bounding box. Color coding and map background again as in \autoref{fig:plot_map}.\label{fig:plot_prediction_grid}](07_prediction_grid_plot.jpeg)

# Acknowledgements

The package benefitted from valuable comments by Joscha Gretzinger, who also suggested the name *bleiglas* (German *Bleiglasfenster* for English *Leadlight*) inspired by the appearance of the cut surface plots.

# References
---
output: github_document
editor_options: 
  chunk_output_type: console
always_allow_html: true
---

[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](https://www.repostatus.org/badges/latest/inactive.svg)](https://www.repostatus.org/#inactive)
![GitHub R package version](https://img.shields.io/github/r-package/v/nevrome/bleiglas)
[![R-CMD-check](https://github.com/nevrome/bleiglas/actions/workflows/check-release.yaml/badge.svg)](https://github.com/nevrome/bleiglas/actions/workflows/check-release.yaml)
[![Coverage Status](https://img.shields.io/codecov/c/github/nevrome/bleiglas/master.svg)](https://codecov.io/github/nevrome/bleiglas?branch=master)
[![license](https://img.shields.io/github/license/nevrome/bleiglas)](https://www.r-project.org/Licenses/MIT)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03092/status.svg)](https://doi.org/10.21105/joss.03092)

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, echo = FALSE}
library(magrittr)
library(knitr)
library(rgl)
library(ggplot2)
knit_hooks$set(webgl = hook_rgl)
view_matrix <- structure(c(0.586383819580078, 0.356217533349991, -0.727502763271332, 
0, -0.810031354427338, 0.257360488176346, -0.526888787746429, 
0, -0.000456457957625389, 0.898260772228241, 0.439460128545761, 
0, 0, 0, 0, 1), .Dim = c(4L, 4L))
```

# bleiglas

bleiglas is an R package that employs [Voro++](http://math.lbl.gov/voro++/) for the calculation of three dimensional Voronoi diagrams from input point clouds. This is a special form of tessellation where each polygon is defined as the area closest to one particular seed point. Voronoi diagrams have useful applications in - among others - astronomy, material science or geography and bleiglas provides functions to make 3D tessellation more readily available as a mean for data visualisation and interpolation. It can be used for any 3D point cloud, but the output is optimized for spatiotemporal applications in archaeology.

1. This README (see Quickstart guide below) describes a basic workflow with code and explains some of my thought process when writing this package.
2. A [JOSS paper](https://doi.org/10.21105/joss.03092) gives some background, introduces the core functions from a more technical point of view and presents an example application.
3. A (rather technical) vignette presents all the code necessary to reproduce the "real world" example application in said JOSS paper. When bleiglas is installed you can open the vignette in R with `vignette("bleiglas_case_study")`.

If you have questions beyond this documentation feel free to open an [issue](https://github.com/nevrome/bleiglas/issues) here on Github. Please also see our [contributing guide](CONTRIBUTING.md).

## Installation 

You can install bleiglas from github

```{r, eval=FALSE}
if(!require('remotes')) install.packages('remotes')
remotes::install_github("nevrome/bleiglas", build_vignettes = TRUE)
```

For the main function `tessellate` you also have to [install the Voro++ software](http://math.lbl.gov/voro++/download/). The package is already available in all major Linux software repositories (on Debian/Ubuntu you can simply run `sudo apt-get install voro++`.). MacOS users should be able to install it via homebrew (`brew install voro++`).

## Quickstart

For this quickstart, we assume you have packages `tidyverse`, `sf`, `rgeos` (which in turn requires the Unix package `geos`) and `c14bazAAR` installed. 

#### Getting some data

I decided to use Dirk Seidenstickers [*Archives des datations radiocarbone d'Afrique centrale*](https://github.com/dirkseidensticker/aDRAC) dataset for this purpose. It includes radiocarbon datings from Central Africa that combine spatial (x & y) and temporal (z) position with some meta information.

<details><summary>Click here for the data preparation steps</summary>
<p>

I selected dates from Cameroon between 1000 and 3000 uncalibrated BP and projected them into a worldwide cylindrical reference system (epsg [4088](https://epsg.io/4088)). As Cameroon is close to the equator this projection should represent distances, angles and areas sufficiently correct for this example exercise. As a minor pre-processing step, I here also remove samples with equal position in all three dimensions for the tessellation.

```{r, message=FALSE}
# download raw data with the data access package c14bazAAR
# c14bazAAR can be installed with
# install.packages("c14bazAAR", repos = c(ropensci = "https://ropensci.r-universe.dev"))
c14_cmr <- c14bazAAR::get_c14data("adrac") %>% 
  # filter data
  dplyr::filter(!is.na(lat) & !is.na(lon), c14age > 1000, c14age < 3000, country == "CMR") 

# remove doubles
c14_cmr_unique <- c14_cmr %>%
  dplyr::mutate(
    rounded_coords_lat = round(lat, 3),
    rounded_coords_lon = round(lon, 3)
  ) %>%
  dplyr::group_by(rounded_coords_lat, rounded_coords_lon, c14age) %>%
  dplyr::filter(dplyr::row_number() == 1) %>%
  dplyr::ungroup()

# transform coordinates
coords <- data.frame(c14_cmr_unique$lon, c14_cmr_unique$lat) %>% 
  sf::st_as_sf(coords = c(1, 2), crs = 4326) %>% 
  sf::st_transform(crs = 4088) %>% 
  sf::st_coordinates()

# create active dataset
c14 <- c14_cmr_unique %>% 
  dplyr::transmute(
    id = seq_len(nrow(.)),
    x = coords[,1], 
    y = coords[,2], 
    z = c14age,
    period = period
)
```

</p>
</details>

<details><summary>Data: <b>c14</b></summary>
<p>

```{r}
c14 
```

</p>
</details>

#### 3D tessellation

[Tessellation](https://en.wikipedia.org/wiki/Tessellation) means filling space with polygons so that neither gaps nor overlaps occur. This is an exciting application for art (e.g. textile art or architecture) and an interesting challenge for mathematics. As a computational archaeologist I was already aware of one particular tessellation algorithm that has quite some relevance for geostatistical analysis like spatial interpolation: Voronoi tilings that are created with [Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation). These are tessellations where each polygon covers the space closest to one of a set of sample points.

<table style="width:100%">
  <tr>
    <th>
      <figure><img src="https://upload.wikimedia.org/wikipedia/commons/thumb/6/66/Ceramic_Tile_Tessellations_in_Marrakech.jpg/320px-Ceramic_Tile_Tessellations_in_Marrakech.jpg" height="150" />
      <br>
      <figcaption>Islamic mosaic with tile tessellations in Marrakech, Morocco. <a href="https://en.wikipedia.org/wiki/File:Ceramic_Tile_Tessellations_in_Marrakech.jpg">wiki</a></figcaption></figure>
    </th>
    <th>
      <figure><img src="https://upload.wikimedia.org/wikipedia/commons/thumb/5/56/Delaunay_Voronoi.svg/441px-Delaunay_Voronoi.svg.png" height="150" />
      <br>
      <figcaption>Delaunay triangulation and its Voronoi diagram. <a href="https://commons.wikimedia.org/wiki/File:Delaunay_Voronoi.svg">wiki</a></figcaption></figure>
    </th>
    <th>
      <figure><img src="https://aip.scitation.org/na101/home/literatum/publisher/aip/journals/content/cha/2009/cha.2009.19.issue-4/1.3215722/production/images/medium/1.3215722.figures.f4.gif" height="150" />
      <br>
      <figcaption>Output example of Voro++ rendered with POV-Ray. <a href="http://math.lbl.gov/voro++">math.lbl.gov</a></figcaption></figure>
    </th>
  <tr>
</table>

It turns out that Voronoi tessellation can be calculated not just for 2D surfaces, but also for higher dimensions. The [Voro++](http://math.lbl.gov/voro++/) software library does exactly this for 3 dimensions. This makes it useful for spatio-temporal applications.

`bleiglas::tessellate()` is a minimal wrapper function that calls the Voro++ command line interface (therefore you have to install Voro++ to use it) for datasets like the one introduced above. We can apply it like this:

```{r}
raw_voro_output <- bleiglas::tessellate(
  c14[, c("id", "x", "y", "z")],
  x_min = min(c14$x) - 150000, x_max = max(c14$x) + 150000, 
  y_min = min(c14$y) - 150000, y_max = max(c14$y) + 150000,
  unit_scaling = c(0.001, 0.001, 1)
)
```

A critical step when using tessellation for spatio-temporal data is a suitable conversion scale between time- and spatial units. Since 3D tessellation crucially depends on the concept of a 3D-distance, we need to make a decision how to combine length- and time-units. Here, for the purpose of this example, we have 1 kilometre correspond to 1 year. Since after the coordinate conversion our spatial units are given in meters, we divide all spatial distances by a factor 1000 to achieve this correspondence: `unit_scaling = c(0.001, 0.001, 1)`.

I decided to increase the size of the tessellation box by 150 kilometres to each (spatial) direction to cover the area of Cameroon. Mind that the scaling factors in `unit_scaling` are also applied to the box size parameters `x_min`, `x_max`, ....

The output of Voro++ is highly customizable, and structurally complex. With the `-v` flag, the voro++ CLI interface prints some config info, which is also the output of `bleiglas::tesselate`:

```
Container geometry        : [937.154:1936.57] [63.1609:1506.58] [1010:2990]
Computational grid size   : 3 by 5 by 6 (estimated from file)
Filename                  : /tmp/RtmpVZjBW3/file3aeb5f400f38
Output string             : %i*%P*%t
Total imported particles  : 392 (4.4 per grid block)
Total V. cells computed   : 392
Total container volume    : 2.8563e+09
Total V. cell volume      : 2.8563e+09
```

It then produces an output file (`*.vol`) that contains all sorts of geometry information for the calculated 3D polygons. `tesselate` returns the content of this file as a character vector with the additionally attached attribute `unit_scaling` (`attributes(raw_voro_output)$unit_scaling`), which is just the scaling vector we put in above. 

I focussed on the edges of the polygons and wrote a parser function `bleiglas::read_polygon_edges()` that can transform the complex Voro++ output for this specific output case to a tidy data.table with six columns: the coordinates (x, y, z) of the start (a) and end point (b) of each polygon edge. A data.table is a tabular R data structure very similar to the standard data.frame. Read more about it [here](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html).

```{r}
polygon_edges <- bleiglas::read_polygon_edges(raw_voro_output)
```

`read_polygon_edges` automatically reverses the rescaling introduced in `tesselate` with the `unit_scaling` attribute.

<details><summary>Data: <b>polygon_edges</b></summary>
<p>

```{r, echo=FALSE}
polygon_edges
```

</p>
</details>

<details><summary>We can plot these polygon edges (black) together with the input sample points (red) in 3D.</summary>
<p>

```{r, webgl=TRUE, fig.width=10, fig.align="center", eval=FALSE}
rgl::axes3d()
rgl::points3d(c14$x, c14$y, c14$z, color = "red")
rgl::aspect3d(1, 1, 1)
rgl::segments3d(
  x = as.vector(t(polygon_edges[,c(1,4)])),
  y = as.vector(t(polygon_edges[,c(2,5)])),
  z = as.vector(t(polygon_edges[,c(3,6)]))
)
rgl::view3d(userMatrix = view_matrix, zoom = 0.9)
```

</p>
</details>

```{r, webgl=TRUE, fig.width=10, fig.align="center", echo=FALSE}
rgl::axes3d()
rgl::points3d(c14$x, c14$y, c14$z, color = "red")
rgl::aspect3d(1, 1, 1)
rgl::segments3d(
  x = as.vector(t(polygon_edges[,c(1,4)])),
  y = as.vector(t(polygon_edges[,c(2,5)])),
  z = as.vector(t(polygon_edges[,c(3,6)]))
)
rgl::view3d(userMatrix = view_matrix, zoom = 0.9)
```

#### Cutting the polygons

This 3D plot, even if rotatable using mouse input, is of rather limited value since it's very hard to read. I therefore wrote `bleiglas::cut_polygons()` that can cut the 3D polygons at different levels of the z-axis. As the function assumes that x and y represent geographic coordinates, the cuts produce sets of spatial 2D polygons for different values of z -- in our example different points in time. The parameter `cuts` takes a numeric vector of cutting points on the z axis. `bleiglas::cut_polygons()` yields a rather raw format for specifying polygons. Another function, `bleiglas::cut_polygons_to_sf()`, transforms it to `sf`. Here `crs` defines the spatial coordinate reference system of x and y to project the resulting 2D polygons correctly.

```{r}
cut_surfaces <- bleiglas::cut_polygons(
  polygon_edges, 
  cuts = c(2500, 2000, 1500)
) %>%
  bleiglas::cut_polygons_to_sf(crs = 4088)
```

<details><summary>Data: <b>cut_surfaces</b></summary>
<p>

```{r, echo=FALSE}
cut_surfaces
```

</p>
</details>

<details><summary>With this data we can plot a matrix of maps that show the cut surfaces.</summary>
<p>

```{r, fig.width=8, fig.align="center", eval=FALSE}
cut_surfaces %>%
  ggplot() +
  geom_sf(
    aes(fill = z), 
    color = "white",
    lwd = 0.2
  ) +
  geom_sf_text(aes(label = id)) +
  facet_wrap(~z) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
```

</p>
</details>

```{r, fig.width=8, fig.align="center", echo=FALSE}
cut_surfaces %>%
  ggplot() +
  geom_sf(
    aes(fill = z), 
    color = "white",
    lwd = 0.2
  ) +
  facet_wrap(~z) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
```

<details><summary>As all input dates come from Cameroon it makes sense to cut the polygon surfaces to the outline of this administrative unit.</summary>
<p>

```{r, warning=FALSE}
cameroon_border <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>% 
  dplyr::filter(name == "Cameroon") %>% 
  sf::st_transform(4088)

cut_surfaces_cropped <- cut_surfaces %>% sf::st_intersection(cameroon_border)
```

```{r, fig.width=8, fig.align="center", eval=FALSE}
cut_surfaces_cropped %>%
  ggplot() +
  geom_sf(
    aes(fill = z), 
    color = "white",
    lwd = 0.2
  ) +
  facet_wrap(~z) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
```

<p>
</details>

```{r, fig.width=8, fig.align="center", echo=FALSE}
cut_surfaces_cropped %>%
  ggplot() +
  geom_sf(
    aes(fill = z), 
    color = "white",
    lwd = 0.2
  ) +
  facet_wrap(~z) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
```


<details><summary>Finally, we can also visualise any point-wise information in our input data as a feature of the tessellation polygons.</summary>
<p>

```{r, warning=FALSE}
cut_surfaces_material <- cut_surfaces_cropped %>%
  dplyr::left_join(
    c14, by = "id"
  )
```

```{r, fig.width=8, fig.align="center", eval=FALSE}
cut_surfaces_material %>%
  ggplot() +
  geom_sf(
    aes(fill = period), 
    color = "white",
    lwd = 0.2
  ) +
  facet_wrap(~z.x) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
```

</p>
</details>

```{r, fig.width=8, fig.align="center", echo=FALSE}
cut_surfaces_material %>%
  ggplot() +
  geom_sf(
    aes(fill = period), 
    color = "white",
    lwd = 0.2
  ) +
  facet_wrap(~z.x) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
```

This quickstart was a simple primer on how to use this package. If you think the final use case wasn't too impressive, take a look at this analysis of Bronze Age burial types through time, as performed in our [JOSS paper](https://github.com/nevrome/bleiglas/blob/master/paper/paper.md) and the [vignette](https://github.com/nevrome/bleiglas/blob/master/vignettes/complete_example.Rmd).

<!-- Add JOSS paper figure here? Just a suggestion as a further teaser. It's just beautiful --> 

## Citation

```{r, echo=F,comment=""}
citation("bleiglas")
```
---
title: "Bleiglas Bronze Age burial rite distribution case study"
output: pdf_document
vignette: >
  %\VignetteIndexEntry{Bleiglas Bronze Age burial rite distribution case study}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

This vignette contains all the code behind the example in the Journal of Open Source Software paper. This can be considered a real life application initially introduced in Schmid 2019 (https://doi.org/10.1177/1059712319860842), so it contains different files in `inst/workflow` in the package directory and does not omit any code for data preparation, manipulation or plotting. See the package README for a minimal Quickstart guide to bleiglas instead.

```{r, echo = FALSE}
code_files_full_path <- list.files(
  system.file("workflow_example", package = "bleiglas", mustWork = T),
  pattern = ".R$", full.names = T
)
code_files <- basename(code_files_full_path)
sourcecode <- lapply(code_files_full_path, function(x) { paste(readLines(x), collapse="\n") })
```

Running the code in this vignette requires some additional packages. You can install them with

``` {r eval=FALSE}
install.packages(
  c("Bchron", "bleiglas", "c14bazAAR", "data.table", "dplyr", 
    "ggplot2", "magrittr", "pbapply", "purrr", "raster", 
    "rnaturalearth", "scatterplot3d", "sf", "tibble"),
  repos = c(
    CRAN = "https://cloud.r-project.org", 
    ropensci = "https://ropensci.r-universe.dev"
  )
)
```

## `r code_files[1]`

This first script downloads and prepares a set of spatial data objects which will later be used for plotting. The research area for this example was arbitrarily defined as a rectangle covering the most dense data accumulations in the relevant subset of the Radon-B database.

```{r eval=FALSE, code=sourcecode[[1]]}
```

## `r code_files[2]`

This script contains the code to filter and prepare radiocarbon dates on graves from Radon-B. For the purpose of tessellation we need the dates to be transformed to a simple table with columns for the spatiotemporal position as well as for the burial type context.

```{r eval=FALSE, code=sourcecode[[2]]}
```

### `dates_prepared`

```{r, echo=FALSE} 
load(file.path(system.file("workflow_example", package = "bleiglas", mustWork = T), "dates_prepared.RData"))
data.table::as.data.table(dates_prepared)
```

## `r code_files[3]`

The script for the first plot in the JOSS paper: A simple map showing the spatial distribution of the observed dates or the graves they represent.

```{r eval=FALSE, code=sourcecode[[3]]}
```

```{r, echo=FALSE, out.width = '100%'}
knitr::include_graphics(
  file.path(system.file("workflow_example", package = "bleiglas", mustWork = T), "03_map_plot.jpeg")
)
```

## `r code_files[4]`

In this script the tessellation is finally performed. The output is a `sf` object with 2D spatial polygons for different time cuts through the 3D tessellated cube. The time slices are already cut to the boundaries of the research area and the European land outline within the latter.

```{r eval=FALSE, code=sourcecode[[4]]}
```

### `vertices`

```{r, echo=FALSE} 
load(file.path(system.file("workflow_example", package = "bleiglas", mustWork = T), "tesselation_calage_center_burial_type.RData"))
```

```{r, echo=FALSE} 
data.table::as.data.table(vertices)
```

### `polygon_edges`

```{r, echo=FALSE} 
polygon_edges
```

### `cut_surfaces[[1]][1:3]`

```{r, echo=FALSE} 
cut_surfaces[[1]][1:3]
```

## `r code_files[5]`

This script illustrates one way to create a (printable) 3D plot of the 3D tessellation output -- so again a plot script behind the second figure in the JOSS article. The input radiocarbon dates are plotted as red dots surrounded by the edges of the 3D polygons voro++ constructs. The spatial research area becomes a spatiotemporal cube.

```{r eval=FALSE, code=sourcecode[[5]]}
```

```{r, echo=FALSE, out.width = '100%'}
knitr::include_graphics(
  file.path(system.file("workflow_example", package = "bleiglas", mustWork = T), "05_3D_plot.jpeg")
)
```

## `r code_files[6]`

Yet another plot script for what we call the "bleiglas" plot (Figure 3 in the JOSS article). The 2D cut polygons are projected onto a map in a diachronic plot matrix. It allows to inspect the 3D tessellation in a easily digestible way, much more human readable than the afore produced 3D plot.

```{r eval=FALSE, code=sourcecode[[6]]}
```

```{r, echo=FALSE, out.width = '100%'}
knitr::include_graphics(
  file.path(system.file("workflow_example", package = "bleiglas", mustWork = T), "06_bleiglas_plot.jpeg")
)
```

## `r code_files[7]`

The final script and cradle of the last JOSS figure applies the bleiglas grid prediction method to account for the temporal uncertainty of the C14 dates. The result is a plot with less aesthetic, but more scientific value, as the known input errors are indicated.

```{r eval=FALSE, code=sourcecode[[7]]}
```

```{r, echo=FALSE, out.width = '100%'}
knitr::include_graphics(
  file.path(system.file("workflow_example", package = "bleiglas", mustWork = T), "07_prediction_grid_plot.jpeg")
)
```

### `prediction (first 10 rows)`

```{r, echo=FALSE} 
load(file.path(system.file("workflow_example", package = "bleiglas", mustWork = T), "prediction_grid_example.RData"))
prediction_grid_example
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bleiglas.R
\docType{package}
\name{bleiglas-package}
\alias{bleiglas}
\alias{bleiglas-package}
\title{bleiglas: Spatiotemporal Data Interpolation and Visualisation based on 3D Tessellation}
\description{
Employs Voro++ for the calculation of three dimensional Voronoi diagrams from
    input point clouds. This is a special form of tessellation where each polygon is defined 
    as the area closest to one particular seed point. Voronoi diagrams have useful applications
    in - among others - astronomy, material science or geography and bleiglas provides functions 
    to make 3D tessellation more readily available as a mean for data visualisation and interpolation.
    It can be used for any 3D point cloud, but the output is optimized for spatiotemporal applications 
    in archaeology.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/nevrome/bleiglas}
  \item Report bugs at \url{https://github.com/nevrome/bleiglas/issues}
}

}
\author{
\strong{Maintainer}: Clemens Schmid \email{clemens@nevrome.de} (\href{https://orcid.org/0000-0003-3448-5715}{ORCID}) [copyright holder]

Other contributors:
\itemize{
  \item Stephan Schiffels (\href{https://orcid.org/0000-0002-1017-9150}{ORCID}) [contributor]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_grid.R
\name{predict_grid}
\alias{predict_grid}
\alias{attribute_grid_points_to_polygons}
\title{predict_grid}
\usage{
predict_grid(x, prediction_grid, unit_scaling = c(1, 1, 1), ...)

attribute_grid_points_to_polygons(prediction_grid, polygon_edges)
}
\arguments{
\item{x}{List of data.tables/data.frames with the input points that define
the tessellation model:
\itemize{
  \item id: id number that is passed to the output polygon (integer)
  \item x: x-axis coordinate (numeric)
  \item y: y-axis coordinate (numeric)
  \item z: z-axis coordinate (numeric)
  \item ...: arbitrary variables
}}

\item{prediction_grid}{data.table/data.frame with the points that should be
predicted by the tessellation model:
\itemize{
  \item x: x-axis coordinate (numeric)
  \item y: y-axis coordinate (numeric)
  \item z: z-axis coordinate (numeric)
}}

\item{unit_scaling}{passed to \link{tessellate} - see the documentation there}

\item{...}{Further variables passed to \code{pbapply::pblapply} (e.g. \code{cl})}

\item{polygon_edges}{polygon points as returned by \code{bleiglas::read_polygon_edges}}
}
\value{
list of data.tables with polygon attribution and predictions
}
\description{
\code{predict_grid} allows to conveniently use the tessellation output as a
model to predict values for arbitrary points. See the bleiglass JOSS paper and
\code{vignette("complete_example", "bleiglas")} for an example application.
\code{attribute_grid_points_to_polygons} is a helper function that does the
important step of point-to-polygon attribution, which might be useful by
itself.
}
\examples{
x <- lapply(1:5, function(i) {
  current_iteration <- data.table::data.table(
    id = 1:5,
    x = c(1, 2, 3, 2, 1) + rnorm(5, 0, 0.3),
    y = c(3, 1, 4, 4, 3) + rnorm(5, 0, 0.3),
    z = c(1, 3, 4, 2, 5) + rnorm(5, 0, 0.3),
    value1 = c("Brot", "Kaese", "Wurst", "Gurke", "Brot"),
    value2 = c(5.3, 5.1, 5.8, 1.0, 1.2)
  )
  data.table::setkey(current_iteration, "x", "y", "z")
  unique(current_iteration)
})

all_iterations <- data.table::rbindlist(x)

prediction_grid <- expand.grid(
  x = seq(min(all_iterations$x), max(all_iterations$x), length.out = 10),
  y = seq(min(all_iterations$y), max(all_iterations$y), length.out = 10),
  z = seq(min(all_iterations$z), max(all_iterations$z), length.out = 5)
)

bleiglas::predict_grid(x, prediction_grid, cl = 1)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_polygon_edges.R
\name{read_polygon_edges}
\alias{read_polygon_edges}
\title{read_polygon_edges}
\usage{
read_polygon_edges(x, rescale = TRUE)
}
\arguments{
\item{x}{character vector with raw, linewise output of voro++ as produced with
\link{tessellate} when \code{output_definition = "\%i*\%P*\%t"}}

\item{rescale}{Should the output of \link{tessellate} be back-rescaled according to its 
\code{unit_scaling} attribute? Ignored if \code{x} does not have this attribute}
}
\value{
\link[data.table]{data.table} with columns for the coordinates x, y and z of the starting and
end point of each polygon edge
}
\description{
Special reader function for polygon edge output of voro++.
}
\examples{
random_unique_points <- unique(data.table::data.table(
  id = NA,
  x = runif(10, 0, 100000),
  y = runif(10, 0, 100000),
  z = runif(10, 0, 100)
))
random_unique_points$id <- seq_len(nrow(random_unique_points))

voro_output <- tessellate(random_unique_points, unit_scaling = c(0.001, 0.001, 1))

polygon_points <- read_polygon_edges(voro_output)

cut_surfaces <- cut_polygons(polygon_points, c(20, 40, 60))

cut_surfaces_sf <- cut_polygons_to_sf(cut_surfaces, crs = 25832)
\donttest{
polygons_z_20 <- sf::st_geometry(cut_surfaces_sf[cut_surfaces_sf$z == 20, ])
plot(polygons_z_20, col = sf::sf.colors(10, categorical = TRUE))
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tessellate.R
\name{tessellate}
\alias{tessellate}
\title{tessellate}
\usage{
tessellate(
  x,
  x_min = NA,
  x_max = NA,
  y_min = NA,
  y_max = NA,
  z_min = NA,
  z_max = NA,
  unit_scaling = c(1, 1, 1),
  output_definition = "\%i*\%P*\%t",
  options = "-v",
  voro_path = "voro++"
)
}
\arguments{
\item{x}{data.table/data.frame with the input points described by four variables (named columns):
\itemize{
  \item id: id number that is passed to the output polygon (integer)
  \item x: x-axis coordinate (numeric)
  \item y: y-axis coordinate (numeric)
  \item z: z-axis coordinate (numeric)
}}

\item{x_min}{minimum x-axis coordinate of the tessellation box. Default: min(x).
These values are automatically multiplied by the scaling factor in \code{unit_scaling}!}

\item{x_max}{maximum x-axis coordinate of the tessellation box. Default: max(x)}

\item{y_min}{minimum y-axis coordinate of the tessellation box. Default: min(y)}

\item{y_max}{maximum y-axis coordinate of the tessellation box. Default: max(y)}

\item{z_min}{minimum z-axis coordinate of the tessellation box. Default: min(z)}

\item{z_max}{maximum z-axis coordinate of the tessellation box. Default: max(z)}

\item{unit_scaling}{numeric vector with 3 scaling factors for x, y and z axis values.
As a default setting (c(1,1,1)) tesselate assumes that the values given as x, y and z are comparable in units.
If you input spatio-temporal data, make sure that you have units that determine your 3D distance
metric the way you intend it to be. For example, if you need 1km=1year, use those units in the input. 
Otherwise, rescale appropriately. Mind that the values of *_min and *_max are adjusted 
as well by these factors. The unit_scaling parameter is stored as an attribute of the output
to scale the output back automatically in \link{read_polygon_edges}.}

\item{output_definition}{string that describes how the output file of voro++ should be structured.
This is passed to the -c option of the command line interface. All possible customization options
are documented \href{http://math.lbl.gov/voro++/doc/custom.html}{here}. Default: "\%i*\%P*\%t"}

\item{options}{string with additional options passed to voro++. All options are documented
\href{http://math.lbl.gov/voro++/doc/cmd.html}{here}. Default: "-v"}

\item{voro_path}{system path to the voro++ executable. Default: "voro++"}
}
\value{
raw, linewise output of voro++ in a character vector with an attribute "unit scaling" (see above)
}
\description{
Command line utility wrapper for the \href{http://math.lbl.gov/voro++}{voro++} software library.
voro++ must be installed on your system to use this function.
}
\examples{
random_unique_points <- unique(data.table::data.table(
  id = NA,
  x = runif(10, 0, 100000),
  y = runif(10, 0, 100000),
  z = runif(10, 0, 100)
))
random_unique_points$id <- seq_len(nrow(random_unique_points))

voro_output <- tessellate(random_unique_points, unit_scaling = c(0.001, 0.001, 1))

polygon_points <- read_polygon_edges(voro_output)

cut_surfaces <- cut_polygons(polygon_points, c(20, 40, 60))

cut_surfaces_sf <- cut_polygons_to_sf(cut_surfaces, crs = 25832)
\donttest{
polygons_z_20 <- sf::st_geometry(cut_surfaces_sf[cut_surfaces_sf$z == 20, ])
plot(polygons_z_20, col = sf::sf.colors(10, categorical = TRUE))
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cut_polygons.R
\name{cut_polygons}
\alias{cut_polygons}
\title{cut_polygons}
\usage{
cut_polygons(x, cuts)
}
\arguments{
\item{x}{\link[data.table]{data.table} with output of voro++ as produced with
\link{tessellate} and then \link{read_polygon_edges}}

\item{cuts}{numeric vector with z-axis coordinates where cuts should be applied}
}
\value{
list of lists. One list element for each cutting surface and within these
data.tables for each 2D polygon that resulted from the cutting operation.
Each data.table holds the corner coordinates for one 2D polygon.
}
\description{
Figuratively cut horizontal slices of a 3D, tessellated cube.
}
\examples{
random_unique_points <- unique(data.table::data.table(
  id = NA,
  x = runif(10, 0, 100000),
  y = runif(10, 0, 100000),
  z = runif(10, 0, 100)
))
random_unique_points$id <- seq_len(nrow(random_unique_points))

voro_output <- tessellate(random_unique_points, unit_scaling = c(0.001, 0.001, 1))

polygon_points <- read_polygon_edges(voro_output)

cut_surfaces <- cut_polygons(polygon_points, c(20, 40, 60))

cut_surfaces_sf <- cut_polygons_to_sf(cut_surfaces, crs = 25832)
\donttest{
polygons_z_20 <- sf::st_geometry(cut_surfaces_sf[cut_surfaces_sf$z == 20, ])
plot(polygons_z_20, col = sf::sf.colors(10, categorical = TRUE))
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cut_polygons_to_sf.R
\name{cut_polygons_to_sf}
\alias{cut_polygons_to_sf}
\title{cut_polygons_to_sf}
\usage{
cut_polygons_to_sf(x, crs)
}
\arguments{
\item{x}{list of lists of data.tables. Output of cut_polygons}

\item{crs}{coordinate reference system of the resulting 2D polygons.
Integer with the EPSG code, or character with proj4string}
}
\value{
sf object with one row for each 2D polygon
}
\description{
Transform the polygon cut slices to the sf format. This only makes sense
if the x and y coordinate of your input dataset are spatial coordinates.
}
\examples{
random_unique_points <- unique(data.table::data.table(
  id = NA,
  x = runif(10, 0, 100000),
  y = runif(10, 0, 100000),
  z = runif(10, 0, 100)
))
random_unique_points$id <- seq_len(nrow(random_unique_points))

voro_output <- tessellate(random_unique_points, unit_scaling = c(0.001, 0.001, 1))

polygon_points <- read_polygon_edges(voro_output)

cut_surfaces <- cut_polygons(polygon_points, c(20, 40, 60))

cut_surfaces_sf <- cut_polygons_to_sf(cut_surfaces, crs = 25832)
\donttest{
polygons_z_20 <- sf::st_geometry(cut_surfaces_sf[cut_surfaces_sf$z == 20, ])
plot(polygons_z_20, col = sf::sf.colors(10, categorical = TRUE))
}

}
