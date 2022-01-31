# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who 
contribute through reporting issues, posting feature requests, updating documentation,
submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for
everyone, regardless of level of experience, gender, gender identity and expression,
sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or
imagery, derogatory comments or personal attacks, trolling, public or private harassment,
insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments,
commits, code, wiki edits, issues, and other contributions that are not aligned to this 
Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed 
from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by 
opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the Contributor Covenant 
(http:contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
# ecoengine

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
![CRAN/GitHub 1.10.0_/1.10.99](https://img.shields.io/badge/CRAN/GitHub-1.10.0_/0.1.9999-blue.svg)

## R interface to the Berkeley Ecoinformatics Engine


**Providing access to UC Berkeley's Natural History Data**


This package provides a R wrapper for the newly available [ecoinformatics engine from UC Berkeley](http://ecoengine.berkeley.edu/). The API is very new and currently provides access to two types of data.

* Georeferenced data from the Wieslander project
* Data on > 2 million georeferenced Berkley museum specimens.

## Package Status and Installation

[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/ecoengine?branch=master&svg=true)](https://ci.appveyor.com/project/ropensci/ecoengine)
[![Travis-CI Build Status](https://travis-ci.org/ropensci/ecoengine.svg?branch=master)](https://travis-ci.org/)
 [![codecov](https://codecov.io/gh/ropensci/ecoengine/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/ecoengine)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/ecoengine?color=blue)](https://github.com/metacran/cranlogs.app)

__Installation instructions__

__Stable Version__

A stable version can be installed from the central CRAN repository

```coffee
install.packages("ecoengine", dependencies = TRUE)
```

__Development Version__
Install the package with `devtools` to obtain the latest development version.

```coffee
# If you don't already have the devtools package installed, run
# install.packages("devtools")
# unlike most packages, devtools requires additional non-R dependencies depending on your OS.
# See → https://github.com/karthik/dlab-advanced-r/blob/master/installation.md#installation
library(devtools)
install_github('ropensci/ecoengine')
```

## Usage 
### Documentation

A quick start guide is available both as [markdown](https://github.com/ropensci/ecoengine/blob/master/inst/Using_ecoengine.md) and a [nicely formatted PDF](https://github.com/ropensci/ecoengine/blob/master/inst/Using_ecoengine.pdf?raw=true) or you can go through a [set of slides](http://karthik.github.io/eeguide) from a recent talk.

## Citation

```coffee
To cite package ‘ecoengine’ in publications use:

  Karthik Ram (2014). ecoengine: Programmatic interface to the API
  serving UC Berkeley's Natural History Data. R package version 1.9.
  https://github.com/ropensci/ecoengine

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {ecoengine: Programmatic interface to the API serving UC Berkeley's Natural History
Data},
    author = {Karthik Ram},
    year = {2014},
    note = {R package version 1.9},
    url = {https://github.com/ropensci/ecoengine},
  }
```

## Resources

* [Ecoengine API documentation](http://ecoengine.berkeley.edu/developers/)
* [Berkeley Natural History Museums](http://bnhm.berkeley.edu/)
* [Guide to using the ecoengine (slides)](http://karthik.github.io/eeguide)
* [How to make large requests to the ecoengine](https://gist.github.com/9360037)


[Bug reports](https://github.com/ropensci/ecoengine/issues/new), feature requests and suggestions (especially as pull requests) are most welcome.


---
  
This package is part of a richer suite called [SPOCC Species Occurrence Data](https://github.com/ropensci/spocc), along with several other packages, that provide access to occurrence records from multiple databases. We recommend using SPOCC as the primary R interface to ecoengine unless your needs are limited to this single source.    


---

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md).
By participating in this project you agree to abide by its terms.



[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
Uwe Ligges writes:

Thanks, we see:

   Found the following (possibly) invalid file URI:
     URI: CONDUCT.md
       From: README.md

Please include the file in the package, link to it via a fully specified
URL, or omit the link in case you do not need it.


I have ignored both files from the build.Guide to using the ecoengine R package
======================================

The Berkeley Ecoengine (<http://ecoengine.berkeley.edu>) provides an
open API to a wealth of museum data contained in the [Berkeley natural
history museums](https://bnhm.berkeley.edu/). This R package provides a
programmatic interface to this rich repository of data allowing for the
data to be easily analyzed and visualized or brought to bear in other
contexts. This vignette provides a brief overview of the package's
capabilities.

The API documentation is available at
<http://ecoengine.berkeley.edu/developers/>. As with most APIs it is
possible to query all the available endpoints that are accessible
through the API itself. Ecoengine has something similar.

    library(ecoengine)
    ee_about()

<table style="width:50%;">
<caption>Table continues below</caption>
<colgroup>
<col width="50%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">type</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">wieslander_vegetation_type_mapping</td>
</tr>
<tr class="even">
<td align="left">wieslander_vegetation_type_mapping</td>
</tr>
<tr class="odd">
<td align="left">wieslander_vegetation_type_mapping</td>
</tr>
<tr class="even">
<td align="left">wieslander_vegetation_type_mapping</td>
</tr>
<tr class="odd">
<td align="left">data</td>
</tr>
<tr class="even">
<td align="left">data</td>
</tr>
<tr class="odd">
<td align="left">data</td>
</tr>
<tr class="even">
<td align="left">data</td>
</tr>
<tr class="odd">
<td align="left">actions</td>
</tr>
<tr class="even">
<td align="left">meta-data</td>
</tr>
<tr class="odd">
<td align="left">meta-data</td>
</tr>
<tr class="even">
<td align="left">meta-data</td>
</tr>
</tbody>
</table>

<table style="width:75%;">
<colgroup>
<col width="75%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">endpoint</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/vtmplots_trees/" class="uri">https://ecoengine.berkeley.edu/api/vtmplots_trees/</a></td>
</tr>
<tr class="even">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/vtmplots/" class="uri">https://ecoengine.berkeley.edu/api/vtmplots/</a></td>
</tr>
<tr class="odd">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/vtmplots_brushes/" class="uri">https://ecoengine.berkeley.edu/api/vtmplots_brushes/</a></td>
</tr>
<tr class="even">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/vtmveg/" class="uri">https://ecoengine.berkeley.edu/api/vtmveg/</a></td>
</tr>
<tr class="odd">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/checklists/" class="uri">https://ecoengine.berkeley.edu/api/checklists/</a></td>
</tr>
<tr class="even">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/sensors/" class="uri">https://ecoengine.berkeley.edu/api/sensors/</a></td>
</tr>
<tr class="odd">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/observations/" class="uri">https://ecoengine.berkeley.edu/api/observations/</a></td>
</tr>
<tr class="even">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/photos/" class="uri">https://ecoengine.berkeley.edu/api/photos/</a></td>
</tr>
<tr class="odd">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/search/" class="uri">https://ecoengine.berkeley.edu/api/search/</a></td>
</tr>
<tr class="even">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/layers/" class="uri">https://ecoengine.berkeley.edu/api/layers/</a></td>
</tr>
<tr class="odd">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/series/" class="uri">https://ecoengine.berkeley.edu/api/series/</a></td>
</tr>
<tr class="even">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/sources/" class="uri">https://ecoengine.berkeley.edu/api/sources/</a></td>
</tr>
</tbody>
</table>

The ecoengine class
-------------------

The data functions in the package include ones that query obervations,
checklists, photos, and vegetation records. These data are all formatted
as a common `S3` class called `ecoengine`. The class includes 4 slots.

-   \[`Total results on server`\] A total result count (not necessarily
    the results in this particular object but the total number available
    for a particlar query)
-   \[`Args`\] The arguments (So a reader can replicate the results or
    rerun the query using other tools.)  
-   \[`Type`\] The type (`photos`, `observation`, or `checklist`)  
-   \[`Number of results retrieved`\] The data. Data are most often
    coerced into a `data.frame`. To access the data simply use
    `result_object$data`.

The default `print` method for the class will summarize the object.

Notes on downloading large data requests
----------------------------------------

For the sake of speed, results are paginated at `1000` results per page.
It is possible to request all pages for any query by specifying
`page = all` in any function that retrieves data. However, this option
should be used if the request is reasonably sized. With larger requests,
there is a chance that the query might become interrupted and you could
lose any data that may have been partially downloaded. In such cases the
recommended practice is to use the returned observations to split the
request. You can always check the number of requests you'll need to
retreive data for any query by running `ee_pages(obj)` where `obj` is an
object of class `ecoengine`.

    request <- ee_photos(county = "Santa Clara County", quiet = TRUE, progress = FALSE)
    # Use quiet to suppress messages. Use progress = FALSE to suppress progress
    # bars which can clutter up documents.
    ee_pages(request)

    #>  [1] 1

    # Now it's simple to parallelize this request You can parallelize across
    # number of cores by passing a vector of pages from 1 through the total
    # available.

### Specimen Observations

The database contains over 2 million records (3427932 total). Many of
these have already been georeferenced. There are two ways to obtain
observations. One is to query the database directly based on a partial
or exact taxonomic match. For example

    pinus_observations <- ee_observations(scientific_name = "Pinus", page = 1, quiet = TRUE, 
        progress = FALSE)
    pinus_observations

    #>  [Total results on the server]: 59543 
    #>  [Args]: 
    #>  country = United States 
    #>  scientific_name = Pinus 
    #>  extra = last_modified 
    #>  georeferenced = FALSE 
    #>  page_size = 1000 
    #>  page = 1 
    #>  [Type]: FeatureCollection 
    #>  [Number of results retrieved]: 1000

For additional fields upon which to query, simply look through the help
for `?ee_observations`. In addition to narrowing data by taxonomic
group, it's also possible to add a bounding box (add argument `bbox`) or
request only data that have been georeferenced (set
`georeferenced = TRUE`).

    lynx_data <- ee_observations(genus = "Lynx", georeferenced = TRUE, quiet = TRUE, 
        progress = FALSE)
    lynx_data

    #>  [Total results on the server]: 725 
    #>  [Args]: 
    #>  country = United States 
    #>  genus = Lynx 
    #>  extra = last_modified 
    #>  georeferenced = True 
    #>  page_size = 1000 
    #>  page = 1 
    #>  [Type]: FeatureCollection 
    #>  [Number of results retrieved]: 725

    # Notice that we only for the first 1000 rows.  But since 795 is not a big
    # request, we can obtain this all in one go.
    lynx_data <- ee_observations(genus = "Lynx", georeferenced = TRUE, page = "all", 
        progress = FALSE)

    #>  Search contains 725 observations (downloading 1 of 1 pages)

    lynx_data

    #>  [Total results on the server]: 725 
    #>  [Args]: 
    #>  country = United States 
    #>  genus = Lynx 
    #>  extra = last_modified 
    #>  georeferenced = True 
    #>  page_size = 1000 
    #>  page = all 
    #>  [Type]: FeatureCollection 
    #>  [Number of results retrieved]: 725

**Other search examples**

    animalia <- ee_observations(kingdom = "Animalia")
    Artemisia <- ee_observations(scientific_name = "Artemisia douglasiana")
    asteraceae <- ee_observationss(family = "asteraceae")
    vulpes <- ee_observations(genus = "vulpes")
    Anas <- ee_observations(scientific_name = "Anas cyanoptera", page = "all")
    loons <- ee_observations(scientific_name = "Gavia immer", page = "all")
    plantae <- ee_observations(kingdom = "plantae")
    # grab first 10 pages (250 results)
    plantae <- ee_observations(kingdom = "plantae", page = 1:10)
    chordata <- ee_observations(phylum = "chordata")
    # Class is clss since the former is a reserved keyword in SQL.
    aves <- ee_observations(clss = "aves")

**Additional Features**

As of July 2014, the API now allows you exclude or request additional
fields from the database, even if they are not directly exposed by the
API. The list of fields are:

`id`, `record`, `source`, `remote_resource`, `begin_date`, `end_date`,
`collection_code`, `institution_code`, `state_province`, `county`,
`last_modified`, `original_id`, `geometry`,
`coordinate_uncertainty_in_meters`, `md5`, `scientific_name`,
`observation_type`, `date_precision`, `locality`,
`earliest_period_or_lowest_system`, `latest_period_or_highest_system`,
`kingdom`, `phylum`, `clss`, `order`, `family`, `genus`,
`specific_epithet`, `infraspecific_epithet`, `minimum_depth_in_meters`,
`maximum_depth_in_meters`, `maximum_elevation_in_meters`,
`minimum_elevation_in_meters`, `catalog_number`, `preparations`, `sex`,
`life_stage`, `water_body`, `country`, `individual_count`,
`associated_resources`

*To request additional fields*

Just pass then in the `extra` field with multiple ones separated by
commas.

    aves <- ee_observations(clss = "aves", extra = "kingdom,genus")

    #>  Search contains 237673 observations (downloading 1 of 238 pages)

    names(aves$data)

    #>   [1] "longitude"                        "latitude"                        
    #>   [3] "type"                             "url"                             
    #>   [5] "record"                           "observation_type"                
    #>   [7] "scientific_name"                  "country"                         
    #>   [9] "state_province"                   "begin_date"                      
    #>  [11] "end_date"                         "source"                          
    #>  [13] "remote_resource"                  "locality"                        
    #>  [15] "coordinate_uncertainty_in_meters" "recorded_by"                     
    #>  [17] "kingdom"                          "genus"                           
    #>  [19] "last_modified"

Similarly use `exclude` to exclude any fields that might be returned by
default.

    aves <- ee_observations(clss = "aves", exclude = "source,remote_resource")

    #>  Search contains 237673 observations (downloading 1 of 238 pages)

    names(aves$data)

    #>   [1] "longitude"                        "latitude"                        
    #>   [3] "type"                             "url"                             
    #>   [5] "record"                           "observation_type"                
    #>   [7] "scientific_name"                  "country"                         
    #>   [9] "state_province"                   "begin_date"                      
    #>  [11] "end_date"                         "locality"                        
    #>  [13] "coordinate_uncertainty_in_meters" "recorded_by"                     
    #>  [15] "last_modified"

**Mapping observations**

The development version of the package includes a new function
`ee_map()` that allows users to generate interactive maps from
observation queries using Leaflet.js.

    lynx_data <- ee_observations(genus = "Lynx", georeferenced = TRUE, page = "all", 
        quiet = TRUE)
    ee_map(lynx_data)

![Map of Lynx observations across North America](map.png)

### Photos

The ecoengine also contains a large number of photos from various
sources. It's easy to query the photo database using similar arguments
as above. One can search by taxa, location, source, collection and much
more.

    photos <- ee_photos(quiet = TRUE, progress = FALSE)
    photos

    #>  [Total results on the server]: 74468 
    #>  [Args]: 
    #>  page_size = 1000 
    #>  georeferenced = 0 
    #>  page = 1 
    #>  [Type]: photos 
    #>  [Number of results retrieved]: 1000

The database currently holds 74468 photos. Photos can be searched by
state province, county, genus, scientific name, authors along with date
bounds. For additional options see `?ee_photos`.

#### Searching photos by author

    charles_results <- ee_photos(authors = "Charles Webber", quiet = TRUE, progress = FALSE)
    charles_results

    #>  [Total results on the server]: 4907 
    #>  [Args]: 
    #>  page_size = 1000 
    #>  authors = Charles Webber 
    #>  georeferenced = FALSE 
    #>  page = 1 
    #>  [Type]: photos 
    #>  [Number of results retrieved]: 1000

    # Let's examine a couple of rows of the data
    charles_results$data[1:2, ]

    #>  # A tibble: 2 x 18
    #>                                                                            url
    #>                                                                          <chr>
    #>  1 https://ecoengine.berkeley.edu/api/photos/CalPhotos%3A0024%2B3291%2B2018%2B
    #>  2 https://ecoengine.berkeley.edu/api/photos/CalPhotos%3A0024%2B3291%2B1998%2B
    #>  # ... with 17 more variables: record <chr>, authors <chr>, locality <chr>,
    #>  #   county <chr>, photog_notes <chr>, begin_date <dttm>, end_date <dttm>,
    #>  #   collection_code <chr>, scientific_name <chr>, url <chr>,
    #>  #   license <chr>, media_url <chr>, remote_resource <chr>, source <chr>,
    #>  #   geometry.type <chr>, longitude <chr>, latitude <chr>

------------------------------------------------------------------------

#### Browsing these photos

    view_photos(charles_results)

This will launch your default browser and render a page with thumbnails
of all images returned by the search query. You can do this with any
`ecoengine` object of type `photos`. Suggestions for improving the photo
browser are welcome.

![](browse_photos.png)

Other photo search examples

    # All the photos in the CDGA collection
    all_cdfa <- ee_photos(collection_code = "CDFA", page = "all", progress = FALSE)
    # All Racoon pictures
    racoons <- ee_photos(scientific_name = "Procyon lotor", quiet = TRUE, progress = FALSE)

------------------------------------------------------------------------

### Species checklists

There is a wealth of checklists from all the source locations. To get
all available checklists from the engine, run:

    all_lists <- ee_checklists()

    #>  Returning 52 checklists

    head(all_lists[, c("footprint", "subject")])

    #>                                                          footprint
    #>  1   https://ecoengine.berkeley.edu/api/footprints/angelo-reserve/
    #>  2   https://ecoengine.berkeley.edu/api/footprints/angelo-reserve/
    #>  3   https://ecoengine.berkeley.edu/api/footprints/angelo-reserve/
    #>  4 https://ecoengine.berkeley.edu/api/footprints/hastings-reserve/
    #>  5   https://ecoengine.berkeley.edu/api/footprints/angelo-reserve/
    #>  6 https://ecoengine.berkeley.edu/api/footprints/hastings-reserve/
    #>       subject
    #>  1    Mammals
    #>  2     Mosses
    #>  3    Beetles
    #>  4    Spiders
    #>  5 Amphibians
    #>  6       Ants

Currently there are 52 lists available. We can drill deeper into any
list to get all the available data. We can also narrow our checklist
search to groups of interest (see `unique(all_lists$subject)`). For
example, to get the list of Spiders:

    spiders <- ee_checklists(subject = "Spiders")

    #>  Found 1 checklists

    spiders

    #>                  record
    #>  4 bigcb:specieslist:15
    #>                                                          footprint
    #>  4 https://ecoengine.berkeley.edu/api/footprints/hastings-reserve/
    #>                                                                        url
    #>  4 https://ecoengine.berkeley.edu/api/checklists/bigcb%3Aspecieslist%3A15/
    #>                                            source subject
    #>  4 https://ecoengine.berkeley.edu/api/sources/18/ Spiders

Now we can drill deep into each list. For this tutorial I'll just
retrieve data from the the two lists returned above.

    library(plyr)

    #>  -------------------------------------------------------------------------

    #>  You have loaded plyr after dplyr - this is likely to cause problems.
    #>  If you need functions from both plyr and dplyr, please load plyr first, then dplyr:
    #>  library(plyr); library(dplyr)

    #>  -------------------------------------------------------------------------

    #>  
    #>  Attaching package: 'plyr'

    #>  The following objects are masked from 'package:dplyr':
    #>  
    #>      arrange, count, desc, failwith, id, mutate, rename, summarise,
    #>      summarize

    #>  The following object is masked from 'package:purrr':
    #>  
    #>      compact

    spider_details <- ldply(spiders$url, checklist_details)
    names(spider_details)

    #>   [1] "url"                              "observation_type"                
    #>   [3] "scientific_name"                  "collection_code"                 
    #>   [5] "institution_code"                 "country"                         
    #>   [7] "state_province"                   "county"                          
    #>   [9] "locality"                         "begin_date"                      
    #>  [11] "end_date"                         "kingdom"                         
    #>  [13] "phylum"                           "clss"                            
    #>  [15] "order"                            "family"                          
    #>  [17] "genus"                            "specific_epithet"                
    #>  [19] "infraspecific_epithet"            "source"                          
    #>  [21] "remote_resource"                  "earliest_period_or_lowest_system"
    #>  [23] "latest_period_or_highest_system"

    unique(spider_details$scientific_name)

    #>   [1] "Holocnemus pluchei"        "Oecobius navus"           
    #>   [3] "Uloborus diversus"         "Neriene litigiosa"        
    #>   [5] "Theridion "                "Tidarren "                
    #>   [7] "Dictyna "                  "Mallos "                  
    #>   [9] "Yorima "                   "Hahnia sanjuanensis"      
    #>  [11] "Cybaeus "                  "Zanomys "                 
    #>  [13] "Anachemmis "               "Titiotus "                
    #>  [15] "Oxyopes scalaris"          "Zora hespera"             
    #>  [17] "Drassinella "              "Phrurotimpus mateonus"    
    #>  [19] "Scotinella "               "Castianeira luctifera"    
    #>  [21] "Meriola californica"       "Drassyllus insularis"     
    #>  [23] "Herpyllus propinquus"      "Micaria utahna"           
    #>  [25] "Trachyzelotes lyonneti"    "Ebo evansae"              
    #>  [27] "Habronattus oregonensis"   "Metaphidippus "           
    #>  [29] "Platycryptus californicus" "Calymmaria "              
    #>  [31] "Frontinella communis"      "Undetermined "            
    #>  [33] "Latrodectus hesperus"

Our resulting dataset now contains 33 unique spider species.

### Searching the engine

The search is elastic by default. One can search for any field in
`ee_observations()` across all available resources. For example,

    # The search function runs an automatic elastic search across all resources
    # available through the engine.
    lynx_results <- ee_search(query = "genus:Lynx")
    lynx_results[, -3]
    # This gives you a breakdown of what's available allowing you dig deeper.

<table style="width:43%;">
<caption>Table continues below</caption>
<colgroup>
<col width="30%" />
<col width="12%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">field</th>
<th align="left">results</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">California</td>
<td align="left">470</td>
</tr>
<tr class="even">
<td align="left">Nevada</td>
<td align="left">105</td>
</tr>
<tr class="odd">
<td align="left">Alaska</td>
<td align="left">82</td>
</tr>
<tr class="even">
<td align="left">British Columbia</td>
<td align="left">47</td>
</tr>
<tr class="odd">
<td align="left">Arizona</td>
<td align="left">36</td>
</tr>
<tr class="even">
<td align="left">Baja California Sur</td>
<td align="left">25</td>
</tr>
<tr class="odd">
<td align="left">Montana</td>
<td align="left">19</td>
</tr>
<tr class="even">
<td align="left">Baja California</td>
<td align="left">16</td>
</tr>
<tr class="odd">
<td align="left">New Mexico</td>
<td align="left">14</td>
</tr>
<tr class="even">
<td align="left">Oregon</td>
<td align="left">13</td>
</tr>
</tbody>
</table>

<table>
<colgroup>
<col width="100%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">search_url</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22California%22">https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22California%22</a></td>
</tr>
<tr class="even">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22Nevada%22">https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22Nevada%22</a></td>
</tr>
<tr class="odd">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22Alaska%22">https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22Alaska%22</a></td>
</tr>
<tr class="even">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22British+Columbia%22">https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22British+Columbia%22</a></td>
</tr>
<tr class="odd">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22Arizona%22">https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22Arizona%22</a></td>
</tr>
<tr class="even">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22Baja+California+Sur%22">https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22Baja+California+Sur%22</a></td>
</tr>
<tr class="odd">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22Montana%22">https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22Montana%22</a></td>
</tr>
<tr class="even">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22Baja+California%22">https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22Baja+California%22</a></td>
</tr>
<tr class="odd">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22New+Mexico%22">https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22New+Mexico%22</a></td>
</tr>
<tr class="even">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22Oregon%22">https://ecoengine.berkeley.edu/api/search/?q=genus%3ALynx&amp;selected_facets=state_province_exact%3A%22Oregon%22</a></td>
</tr>
</tbody>
</table>

Similarly it's possible to search through the observations in a detailed
manner as well.

    all_lynx_data <- ee_search_obs(query = "Lynx", page = "all", progress = FALSE)

    #>  Search contains 1033 observations (downloading 2 of 2 pages)

    all_lynx_data

    #>  [Total results on the server]: 1033 
    #>  [Args]: 
    #>  q = Lynx 
    #>  page_size = 1000 
    #>  page = all 
    #>  [Type]: observations 
    #>  [Number of results retrieved]: 1000

------------------------------------------------------------------------

### Miscellaneous functions

**Footprints**

`ee_footprints()` provides a list of all the footprints.

    footprints <- ee_footprints()
    footprints[, -3]  # To keep the table from spilling over

<table style="width:33%;">
<caption>Table continues below</caption>
<colgroup>
<col width="33%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">name</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">Angelo Reserve</td>
</tr>
<tr class="even">
<td align="left">Sagehen Reserve</td>
</tr>
<tr class="odd">
<td align="left">Hastings Reserve</td>
</tr>
<tr class="even">
<td align="left">Blue Oak Ranch Reserve</td>
</tr>
</tbody>
</table>

<table style="width:99%;">
<colgroup>
<col width="98%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">url</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/footprints/angelo-reserve/" class="uri">https://ecoengine.berkeley.edu/api/footprints/angelo-reserve/</a></td>
</tr>
<tr class="even">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/footprints/sagehen-reserve/" class="uri">https://ecoengine.berkeley.edu/api/footprints/sagehen-reserve/</a></td>
</tr>
<tr class="odd">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/footprints/hastings-reserve/" class="uri">https://ecoengine.berkeley.edu/api/footprints/hastings-reserve/</a></td>
</tr>
<tr class="even">
<td align="left"><a href="https://ecoengine.berkeley.edu/api/footprints/blue-oak-ranch-reserve/" class="uri">https://ecoengine.berkeley.edu/api/footprints/blue-oak-ranch-reserve/</a></td>
</tr>
</tbody>
</table>

**Data sources**

`ee_sources()` provides a list of data sources for the specimens
contained in the museum.

    source_list <- ee_sources()
    unique(source_list$name)

<table style="width:39%;">
<colgroup>
<col width="38%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">name</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">MVZ Birds Observations</td>
</tr>
<tr class="even">
<td align="left">MVZ Birds Eggs and Nests</td>
</tr>
<tr class="odd">
<td align="left">MVZ Herp Collection</td>
</tr>
<tr class="even">
<td align="left">MVZ Herp Observations</td>
</tr>
<tr class="odd">
<td align="left">MVZ Hildebrand Collection</td>
</tr>
<tr class="even">
<td align="left">MVZ Mammals</td>
</tr>
<tr class="odd">
<td align="left">BIGCB Species Checklists</td>
</tr>
<tr class="even">
<td align="left">MVZ Mammals Observations</td>
</tr>
<tr class="odd">
<td align="left">Essig Museum of Entymology</td>
</tr>
<tr class="even">
<td align="left">UC Museum for Paleontology</td>
</tr>
</tbody>
</table>

    devtools::session_info()

    #>  Session info -------------------------------------------------------------

    #>   setting  value                       
    #>   version  R version 3.3.3 (2017-03-06)
    #>   system   x86_64, darwin13.4.0        
    #>   ui       X11                         
    #>   language (EN)                        
    #>   collate  en_US.UTF-8                 
    #>   tz       America/Los_Angeles         
    #>   date     2017-08-15

    #>  Packages -----------------------------------------------------------------

    #>   package    * version    date       source                           
    #>   assertthat   0.2.0      2017-04-11 CRAN (R 3.3.2)                   
    #>   backports    1.1.0      2017-05-22 CRAN (R 3.3.2)                   
    #>   base       * 3.3.3      2017-03-07 local                            
    #>   bindr        0.1        2016-11-13 cran (@0.1)                      
    #>   bindrcpp     0.2        2017-06-17 CRAN (R 3.3.2)                   
    #>   broom        0.4.2      2017-02-13 CRAN (R 3.3.2)                   
    #>   cellranger   1.1.0      2016-07-27 CRAN (R 3.3.0)                   
    #>   codetools    0.2-15     2016-10-05 CRAN (R 3.3.3)                   
    #>   colorspace   1.3-2      2016-12-14 CRAN (R 3.3.2)                   
    #>   coyote       0.3        2017-01-06 Github (karthik/coyote@acad466)  
    #>   curl         2.8.1      2017-07-21 CRAN (R 3.3.2)                   
    #>   data.table   1.10.4     2017-02-01 CRAN (R 3.3.2)                   
    #>   datasets   * 3.3.3      2017-03-07 local                            
    #>   devtools   * 1.13.2     2017-06-02 CRAN (R 3.3.3)                   
    #>   digest       0.6.12     2017-01-27 CRAN (R 3.3.2)                   
    #>   dplyr      * 0.7.2      2017-07-20 CRAN (R 3.3.2)                   
    #>   ecoengine  * 1.11       2017-07-01 local                            
    #>   evaluate     0.10.1     2017-06-24 CRAN (R 3.3.3)                   
    #>   forcats      0.2.0      2017-01-23 CRAN (R 3.3.2)                   
    #>   foreign      0.8-69     2017-06-21 CRAN (R 3.3.2)                   
    #>   formatR      1.5        2017-04-25 CRAN (R 3.3.2)                   
    #>   ggplot2    * 2.2.1.9000 2017-05-04 Github (hadley/ggplot2@f4398b6)  
    #>   glue         1.1.1      2017-06-21 CRAN (R 3.3.2)                   
    #>   graphics   * 3.3.3      2017-03-07 local                            
    #>   grDevices  * 3.3.3      2017-03-07 local                            
    #>   grid         3.3.3      2017-03-07 local                            
    #>   gtable       0.2.0      2016-02-26 CRAN (R 3.3.0)                   
    #>   haven        1.1.0      2017-07-09 CRAN (R 3.3.2)                   
    #>   hms          0.3        2016-11-22 CRAN (R 3.3.2)                   
    #>   htmltools    0.3.6      2017-04-28 CRAN (R 3.3.3)                   
    #>   httr         1.2.1      2016-07-03 CRAN (R 3.3.0)                   
    #>   jsonlite     1.5        2017-06-01 CRAN (R 3.3.2)                   
    #>   knitr        1.16       2017-05-18 CRAN (R 3.3.3)                   
    #>   lattice      0.20-35    2017-03-25 CRAN (R 3.3.2)                   
    #>   lazyeval     0.2.0      2016-06-12 CRAN (R 3.3.0)                   
    #>   lubridate    1.6.0      2016-09-13 CRAN (R 3.3.0)                   
    #>   magrittr     1.5        2014-11-22 CRAN (R 3.3.0)                   
    #>   memoise      1.1.0      2017-04-21 CRAN (R 3.3.2)                   
    #>   methods    * 3.3.3      2017-03-07 local                            
    #>   mnormt       1.5-5      2016-10-15 CRAN (R 3.3.0)                   
    #>   modelr       0.1.1      2017-07-24 CRAN (R 3.3.2)                   
    #>   munsell      0.4.3      2016-02-13 CRAN (R 3.3.0)                   
    #>   nlme         3.1-131    2017-02-06 CRAN (R 3.3.2)                   
    #>   pander     * 0.6.1      2017-08-06 CRAN (R 3.3.2)                   
    #>   parallel     3.3.3      2017-03-07 local                            
    #>   pkgconfig    2.0.1      2017-03-21 cran (@2.0.1)                    
    #>   plyr       * 1.8.4      2016-06-08 CRAN (R 3.3.0)                   
    #>   psych        1.7.5      2017-05-03 CRAN (R 3.3.3)                   
    #>   purrr      * 0.2.3      2017-08-02 CRAN (R 3.3.2)                   
    #>   R6           2.2.2      2017-06-17 CRAN (R 3.3.2)                   
    #>   Rcpp         0.12.12    2017-07-15 CRAN (R 3.3.2)                   
    #>   readr      * 1.1.1      2017-05-16 CRAN (R 3.3.3)                   
    #>   readxl       1.0.0      2017-04-18 CRAN (R 3.3.3)                   
    #>   reshape2     1.4.2      2016-10-22 CRAN (R 3.3.0)                   
    #>   rlang        0.1.2      2017-08-09 CRAN (R 3.3.2)                   
    #>   rmarkdown  * 1.6        2017-06-15 CRAN (R 3.3.2)                   
    #>   rprojroot    1.2        2017-01-16 cran (@1.2)                      
    #>   rvest        0.3.2      2016-06-17 CRAN (R 3.3.0)                   
    #>   scales       0.4.1      2016-11-09 CRAN (R 3.3.2)                   
    #>   stats      * 3.3.3      2017-03-07 local                            
    #>   stringi      1.1.5      2017-04-07 CRAN (R 3.3.2)                   
    #>   stringr      1.2.0      2017-02-18 cran (@1.2.0)                    
    #>   tibble     * 1.3.3      2017-05-30 Github (tidyverse/tibble@b2275d5)
    #>   tidyr      * 0.6.3      2017-05-15 CRAN (R 3.3.2)                   
    #>   tidyverse  * 1.1.1      2017-01-27 CRAN (R 3.3.2)                   
    #>   tools        3.3.3      2017-03-07 local                            
    #>   utils      * 3.3.3      2017-03-07 local                            
    #>   whisker      0.3-2      2013-04-28 CRAN (R 3.3.0)                   
    #>   withr        2.0.0      2017-07-28 CRAN (R 3.3.2)                   
    #>   xml2         1.1.1      2017-01-24 CRAN (R 3.3.2)

Please send any comments, questions, or ideas for new functionality or
improvements to
&lt;[karthik.ram@berkeley.edu](karthik.ram@berkeley.edu)&gt;. The code
lives on GitHub [under the rOpenSci
account](https://github.com/ropensci/ecoengine). Pull requests and [bug
reports](https://github.com/ropensci/ecoengine/issues?state=open) are
most welcome.

Karthik Ram  
Aug, 2017  
*Santa Clara, California*


```{r, echo = FALSE, message = FALSE}
knitr::opts_chunk$set(
  comment = "#> ",
  error = FALSE,
  cache = TRUE,
  warning = FALSE, 
  tidy = TRUE
)
```

# Guide to using the ecoengine R package

The Berkeley Ecoengine ([http://ecoengine.berkeley.edu](http://ecoengine.berkeley.edu)) provides an open API to a wealth of museum data contained in the [Berkeley natural history museums](https://bnhm.berkeley.edu/). This R package provides a programmatic interface to this rich repository of data allowing for the data to be easily analyzed and visualized or brought to bear in other contexts. This vignette provides a brief overview of the package's capabilities. 

The API documentation is available at [http://ecoengine.berkeley.edu/developers/](http://ecoengine.berkeley.edu/developers/). As with most APIs it is possible to query all the available endpoints that are accessible through the API itself. Ecoengine has something similar.

```{r, about, echo = FALSE}
suppressPackageStartupMessages(library(ecoengine))
suppressPackageStartupMessages(library(pander))
```

```{r, about_ee_dont, eval = FALSE, tidy = TRUE}
library(ecoengine)
ee_about()
```

```{r, about_ee, results = "asis", echo = FALSE}
pandoc.table(ee_about(), justify = "left")
```

## The ecoengine class

The data functions in the package include ones that query obervations, checklists, photos, and vegetation records. These data are all formatted as a common `S3` class called `ecoengine`. The class includes 4 slots.

- [`Total results on server`] A total result count (not necessarily the results in this particular object but the total number available for a particlar query)
- [`Args`] The arguments  (So a reader can replicate the results or rerun the query using other tools.)  
- [`Type`] The type (`photos`, `observation`, or `checklist`)  
- [`Number of results retrieved`] The data. Data are most often coerced into a `data.frame`. To access the data simply use `result_object$data`.  

The default `print` method for the class will summarize the object.

## Notes on downloading large data requests

For the sake of speed, results are paginated at `1000` results per page. It is possible to request all pages for any query by specifying `page = all` in any function that retrieves data. However, this option should be used if the request is reasonably sized. With larger requests, there is a chance that the query might become interrupted and you could lose any data that may have been partially downloaded. In such cases the recommended practice is to use the returned observations to split the request. You can always check the number of requests you'll need to retreive data for any query by running `ee_pages(obj)` where `obj` is an object of class `ecoengine`.

```{r, pagination, eval = TRUE}
request <- ee_photos(county = "Santa Clara County", quiet = TRUE, progress = FALSE)
# Use quiet to suppress messages. Use progress = FALSE to suppress progress bars which can clutter up documents. 
ee_pages(request)
# Now it's simple to parallelize this request
# You can parallelize across number of cores by passing a vector of pages from 1 through the total available.
```


### Specimen Observations

```{r, obs_counts, echo = FALSE, message = FALSE}
x <- ee_observations(quiet = TRUE, progress = FALSE)
```

The database contains over 2 million records (`r format(x$results, nsmall = 0)` total). Many of these have already been georeferenced. There are two ways to obtain observations. One is to query the database directly based on a partial or exact taxonomic match. For example

```{r, observations_1}
pinus_observations <- ee_observations(scientific_name = "Pinus", page = 1, quiet = TRUE, progress = FALSE)
pinus_observations
```

For additional fields upon which to query, simply look through the help for `?ee_observations`. In addition to narrowing data by taxonomic group, it's also possible to add a bounding box (add argument `bbox`) or request only data that have been georeferenced (set `georeferenced = TRUE`). 

```{r, lynx_data, cache = TRUE}
lynx_data <- ee_observations(genus = "Lynx",georeferenced = TRUE, quiet = TRUE, progress = FALSE)
lynx_data
# Notice that we only for the first 1000 rows.
# But since 795 is not a big request, we can obtain this all in one go.
lynx_data <- ee_observations(genus = "Lynx", georeferenced = TRUE, page = "all", progress = FALSE)
lynx_data
```

__Other search examples__

```{r, other_obs, eval = FALSE}
animalia <- ee_observations(kingdom = "Animalia")
Artemisia <- ee_observations(scientific_name = "Artemisia douglasiana")
asteraceae <- ee_observationss(family = "asteraceae")
vulpes <- ee_observations(genus = "vulpes")
Anas <- ee_observations(scientific_name = "Anas cyanoptera", page = "all")
loons <- ee_observations(scientific_name = "Gavia immer", page = "all")
plantae <- ee_observations(kingdom = "plantae")
# grab first 10 pages (250 results)
plantae <- ee_observations(kingdom = "plantae", page = 1:10)
chordata <- ee_observations(phylum = "chordata")
# Class is clss since the former is a reserved keyword in SQL.
aves <- ee_observations(clss = "aves")
```

__Additional Features__

As of July 2014, the API now allows you exclude or request additional fields from the database, even if they are not directly exposed by the API. The list of fields are:


`id`, `record`, `source`, `remote_resource`, `begin_date`, `end_date`, `collection_code`, `institution_code`, `state_province`, `county`, `last_modified`, `original_id`, `geometry`, `coordinate_uncertainty_in_meters`, `md5`, `scientific_name`, `observation_type`, `date_precision`, `locality`, `earliest_period_or_lowest_system`, `latest_period_or_highest_system`, `kingdom`, `phylum`, `clss`, `order`, `family`, `genus`, `specific_epithet`, 
`infraspecific_epithet`, `minimum_depth_in_meters`, `maximum_depth_in_meters`, `maximum_elevation_in_meters`, `minimum_elevation_in_meters`, `catalog_number`, `preparations`, `sex`, `life_stage`, `water_body`, `country`, `individual_count`, `associated_resources`


_To request additional fields_

Just pass then in the `extra` field with multiple ones separated by commas.

```{r request_fields}
aves <- ee_observations(clss = "aves", extra = "kingdom,genus")
names(aves$data)
```
Similarly use `exclude` to exclude any fields that might be returned by default.

```{r exclude_fields}
aves <- ee_observations(clss = "aves", exclude = "source,remote_resource")
names(aves$data)
```


__Mapping observations__

The development version of the package includes a new function `ee_map()` that allows users to generate interactive maps from observation queries using Leaflet.js. 

```{r, eval = FALSE}
lynx_data <- ee_observations(genus = "Lynx", georeferenced = TRUE, page = "all", quiet = TRUE)
ee_map(lynx_data)
```

![Map of Lynx observations across North America](map.png)

### Photos  

The ecoengine also contains a large number of photos from various sources. It's easy to query the photo database using similar arguments as above. One can search by taxa, location, source, collection and much more.


```{r, photo_count}
photos <- ee_photos(quiet = TRUE, progress = FALSE)
photos
```
The database currently holds `r format(photos$results, nsmall = 0)` photos. Photos can be searched by state province, county, genus, scientific name, authors along with date bounds. For additional options see `?ee_photos`.


#### Searching photos by author

```{r, photos_by_author, tidy = TRUE, width.cutoff = 60, background = '#F7F7F7'}
charles_results <- ee_photos(authors = "Charles Webber", quiet = TRUE, progress = FALSE)
charles_results
# Let's examine a couple of rows of the data
charles_results$data[1:2, ]
```
---  

#### Browsing these photos

```{r, browsing_photos, eval = FALSE}
view_photos(charles_results)
```
This will launch your default browser and render a page with thumbnails of all images returned by the search query. You can do this with any `ecoengine` object of type `photos`. Suggestions for improving the photo browser are welcome.


![](browse_photos.png)

Other photo search examples


```{r, photo_examples, eval = FALSE}
# All the photos in the CDGA collection
all_cdfa <- ee_photos(collection_code = "CDFA", page = "all", progress = FALSE)
# All Racoon pictures
racoons <- ee_photos(scientific_name = "Procyon lotor", quiet = TRUE, progress = FALSE)
```

---  


### Species checklists

There is a wealth of checklists from all the source locations. To get all available checklists from the engine, run: 
  
```{r, checklists}
all_lists  <- ee_checklists()
head(all_lists[, c("footprint", "subject")])
```  
Currently there are `r nrow(all_lists)` lists available. We can drill deeper into any list to get all the available data. We can also narrow our checklist search to groups of interest (see `unique(all_lists$subject)`). For example, to get the list of Spiders:

```{r, checklist_spiders, cache = TRUE}
spiders  <- ee_checklists(subject = "Spiders")
spiders
```  

Now we can drill deep into each list. For this tutorial I'll just retrieve data from the the two lists returned above.

```{r, checklist_details}
library(plyr)
spider_details <- ldply(spiders$url, checklist_details)
names(spider_details)
unique(spider_details$scientific_name)
```

Our resulting dataset now contains `r length(unique(spider_details$scientific_name))` unique spider species. 


### Searching the engine  

The search is elastic by default. One can search for any field in `ee_observations()` across all available resources. For example, 


```{r, search, eval = FALSE}
# The search function runs an automatic elastic search across all resources available through the engine. 
lynx_results <- ee_search(query = "genus:Lynx")
lynx_results[, -3]
# This gives you a breakdown of what's available allowing you dig deeper.
```

```{r, search_print, eval = TRUE, results = "asis", echo = FALSE}
lynx_results <- ee_search(query = "genus:Lynx")
pander::pandoc.table(lynx_results[[1]], justify = "left")
```

Similarly it's possible to search through the observations in a detailed manner as well.

```{r, ee_obs_search}
all_lynx_data <- ee_search_obs(query  = "Lynx", page = "all", progress = FALSE)
all_lynx_data
```

---


### Miscellaneous functions

__Footprints__

`ee_footprints()` provides a list of all the footprints.

```{r, footprints_notrun, results = "asis", eval = FALSE, echo = TRUE}
footprints <- ee_footprints()
footprints[, -3] # To keep the table from spilling over
```

```{r, footprints, results = "asis", echo = FALSE}
footprints <- ee_footprints()
pandoc.table(footprints[, -3], justify = "left")
```


__Data sources__

`ee_sources()` provides a list of data sources for the specimens contained in the museum.

```{r, results = "asis", eval = FALSE}
source_list <- ee_sources()
unique(source_list$name)
```

```{r, results = "asis", echo = FALSE}
source_list <- ee_sources()
pandoc.table(data.frame(name = unique(source_list$name)), justify = "left")
```

```{r, version}
devtools::session_info()
```

Please send any comments, questions, or ideas for new functionality or improvements to <[karthik.ram@berkeley.edu](karthik.ram@berkeley.edu)>. The code lives on GitHub [under the rOpenSci account](https://github.com/ropensci/ecoengine). Pull requests and [bug reports](https://github.com/ropensci/ecoengine/issues?state=open) are most welcome.

```{r location, eval = TRUE, echo = FALSE}
library(httr)
x <- content(GET("http://ipinfo.io/"), as = "parsed")
```

 Karthik Ram  
 `r library(lubridate); as.character(lubridate::month(lubridate::now(), label = TRUE))`, `r lubridate::year(lubridate::now())`   
 _`r x$city`, `r x$region`_
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_footprints.R
\name{ee_footprints}
\alias{ee_footprints}
\title{ee_footprints}
\usage{
ee_footprints(foptions = list())
}
\arguments{
\item{foptions}{A list of additional arguments for \code{httr}. There is no reason to use this argument except for debugging purposes.}
}
\value{
data.frame
}
\description{
List of ecoengine footprints.
}
\examples{
# footprints <- ee_footprints()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ecoengine_package.R
\docType{package}
\name{ecoengine}
\alias{ecoengine}
\title{ecoengine}
\description{
\href{http://ecoengine.berkeley.edu}{The Berkeley Ecoengine} provides an open API to a wealth of museum data contained in the  \href{http://bnhm.berkeley.edu/}{Berkeley natural history museums}. This R package provides a programmatic interface to this rich repository of data allowing for the data to be easily analyzed and visualized or brought to bear in other contexts. This vignette provides a brief overview of the package's capabilities. The current version of the package works with API version 0.0.3
}
\details{
\href{http://ecoengine.berkeley.edu}{The Berkeley Ecoengine} provides an open API to a wealth of museum data contained in the  \href{http://bnhm.berkeley.edu/}{Berkeley natural history museums}. This R package provides a programmatic interface to this rich repository of data allowing for the data to be easily analyzed and visualized or brought to bear in other contexts. This vignette provides a brief overview of the package's capabilities. The current version of the package works with API version 0.0.3
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_checklists.R
\name{checklist_details}
\alias{checklist_details}
\title{Checklist details}
\usage{
checklist_details(list_name, ...)
}
\arguments{
\item{list_name}{URL of a checklist}

\item{...}{Additional arguments (currently not implemented)}
}
\value{
\code{data.frame}
}
\description{
Will return details on any checklist
}
\examples{
\dontrun{
spiders  <- ee_checklists(subject = "Spiders")
# Now retrieve all the details for each species on both lists
library(plyr)
spider_details <- ldply(spiders$url, checklist_details)
}
}
\seealso{
\code{\link{ee_checklists}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_search.R
\name{ee_search}
\alias{ee_search}
\title{Ecoengine search}
\usage{
ee_search(query = NULL, foptions = list())
}
\arguments{
\item{query}{search term}

\item{foptions}{A list of additional arguments for \code{httr}. There is no reason to use this argument except for debugging purposes.}
}
\description{
Search across the entire ecoengine database.
}
\examples{
\dontrun{
# lynx_results <- ee_search(query = "genus:Lynx")
}
}
\keyword{search}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ecoengine_package.R
\docType{data}
\name{full_sensor_list}
\alias{full_sensor_list}
\title{full_sensor_list}
\format{A data table with 40 rows and 13 variables. Contains units, variable, data_url, source, record, site_code,   method_name, begin_date, end_date, station_name, geojson.type, geojson.coordinates1, geojson.coordinates2}
\usage{
data(full_sensor_list)
}
\description{
A data.frame containing all the sensors currently available through the engine. This list can be updated anytime by calling the ee_sensors() function.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ecoengine_package.R
\docType{data}
\name{california_counties}
\alias{california_counties}
\title{california_counties}
\format{A data table with 58 rows and 3 variables. Contains county name, FIPS state code (06 in this case) and FIPS county code. See \href{http://en.wikipedia.org/wiki/FIPS_county_code}{Wikipedia} for more information.}
\usage{
data(california_counties)
}
\description{
A data.frame containing list of California counties. Useful for most functions since the data is primarily state based.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_sources.R
\name{ee_sources}
\alias{ee_sources}
\title{Ecoengine data sources}
\usage{
ee_sources(foptions = list())
}
\arguments{
\item{foptions}{A list of additional arguments for \code{httr}. There is no reason to use this argument except for debugging purposes.}
}
\value{
\code{data.frame}
}
\description{
Returns a full list of data sources supported by the ecoengine
}
\examples{
# source_list <- ee_sources()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_about.R
\name{ee_about}
\alias{ee_about}
\title{About the Berkeley Ecoinformatics Engine}
\usage{
ee_about(as.df = TRUE, type = NA)
}
\arguments{
\item{as.df}{\code{FALSE} Returns a list unless this set to \code{TRUE}}

\item{type}{The type of end point. Options include \code{data}, \code{meta-data}, and \code{actions}}
}
\value{
\code{list}
}
\description{
Function returns the current status of fast-evolving API. Returns endpoints and category. Default return is a \code{list} but one can also request a nicely formatted \code{data.frame} by setting the \code{as.df} argument to \code{TRUE}.
}
\examples{
  
\dontrun{
ee_about()
# set as.df = FALSE to return a list rather than a data.frame
ee_about(as.df = FALSE)
# You can also filter by methods by data, meta-data, and actions.
# ee_about(type = "data")
# ee_about(type = "meta-data")
# ee_about(type = "actions")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_utils.R
\name{ee_pages}
\alias{ee_pages}
\title{ee_pages - Returns total number of pages for any ecoengine request}
\usage{
ee_pages(ee, page_size = 1000)
}
\arguments{
\item{ee}{Object of class \code{ecoengine}}

\item{page_size}{Default page size. Currently set to \code{1000} package wide.}
}
\value{
integer
}
\description{
ee_pages - Returns total number of pages for any ecoengine request
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_browse_photos.R
\name{view_photos}
\alias{view_photos}
\title{Browse photo queries in your default browser.}
\usage{
view_photos(input = NULL, output = NULL, browse = TRUE)
}
\arguments{
\item{input}{Input, usually output from a call to \code{\link[ecoengine]{ee_photos}}}

\item{output}{Path and file name for output file. If NULL, a temp file is used.}

\item{browse}{Browse file in your default browse immediately after file creation.
If \code{FALSE}, the file is written, but not opened.}
}
\description{
Browse photo queries in your default browser.
}
\examples{
# view_photos(ee_photos())
# Pictures of racoons
# view_photos(ee_photos(scientific_name = "Procyon lotor", quiet = TRUE))
# or the California Condor
# view_photos(ee_photos(scientific_name = "Gymnogyps californianus", quiet = TRUE))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_utils.R
\name{ee_paginator}
\alias{ee_paginator}
\title{ecoengine paginator}
\usage{
ee_paginator(page, total_obs, page_size = 1000)
}
\arguments{
\item{page}{requested page number or page range. Can also be "all"}

\item{total_obs}{Total number of records available for any search query}

\item{page_size}{Default is \code{25}. Set higher if needed.}
}
\description{
Takes a page range and total number of observations to return the right sequence of pages that need to be crawled.
}
\examples{
\dontrun{
ee_paginator(1, 100)
ee_paginator("all", 100)
ee_paginator(1:2, 1000)
ee_paginator(1:4, 3800)
# This will return an error since there are only 4 pages per 100 observations
ee_paginator(1:5, 10000)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_search.R
\name{ee_search_obs}
\alias{ee_search_obs}
\title{Search observations}
\usage{
ee_search_obs(
  query = NULL,
  page = NULL,
  page_size = 1000,
  quiet = FALSE,
  progress = TRUE,
  foptions = list()
)
}
\arguments{
\item{query}{= The search term}

\item{page}{Page number}

\item{page_size}{The number of observations per page returned by the query. The package default is \code{250} observations. The API default, however, is \code{10}. Setting this number higher will require fewer calls (i.e. pages) but a large size may also fail or timeout.}

\item{quiet}{Default is \code{FALSE}. Set to \code{TRUE} to suppress messages.}

\item{progress}{Progress bar is set to \code{TRUE} by default. Use \code{FALSE} especially when calling functions inside documents.}

\item{foptions}{A list of additional arguments for \code{httr}. There is no reason to use this argument except for debugging purposes.}
}
\value{
data.frame
}
\description{
A powerful way to search through the observations.
}
\examples{
# general_lynx_query <- ee_search_obs(query  = "Lynx")
 \dontrun{
# lynx_data <- ee_search_obs(query  = "genus:Lynx")
# all_lynx_data <- ee_search_obs(query  = "Lynx", page = "all")
}
}
\seealso{
\code{\link{ee_search})}
}
\keyword{search}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_photos.R
\name{ee_photos}
\alias{ee_photos}
\title{ee_photos}
\usage{
ee_photos(
  page = NULL,
  state_province = NULL,
  county = NULL,
  genus = NULL,
  scientific_name = NULL,
  authors = NULL,
  remote_id = NULL,
  collection_code = NULL,
  source = NULL,
  min_date = NULL,
  max_date = NULL,
  related_type = NULL,
  related = NULL,
  page_size = 1000,
  quiet = FALSE,
  georeferenced = FALSE,
  progress = FALSE,
  other_catalog_numbers = NULL,
  foptions = list()
)
}
\arguments{
\item{page}{Page number}

\item{state_province}{Need to describe these parameters}

\item{county}{California counties. Package include a full list of counties. To load dataset \code{data(california_counties)}}

\item{genus}{genus name}

\item{scientific_name}{scientific name}

\item{authors}{author name}

\item{remote_id}{remote id}

\item{collection_code}{Type of collection. Can be \code{CalAcademy}, \code{Private}, \code{VTM}, \code{CDFA}. \code{CalFlora} Others TBA}

\item{source}{data source. See \code{\link{ee_sources}}}

\item{min_date}{Lower date bound. Most be in format Y-m-D h:m:s. This package uses lubridate's \code{ymd_hms} function to format dates.}

\item{max_date}{upper date bound. Most be in format Y-m-D h:m:s. This package uses lubridate's \code{ymd_hms} function to format dates.}

\item{related_type}{Need to describe these parameters}

\item{related}{Need to describe these parameters}

\item{page_size}{The number of observations per page returned by the query. The package default is \code{250} observations. The API default, however, is \code{10}. Setting this number higher will require fewer calls (i.e. pages) but a large size may also fail or timeout.}

\item{quiet}{Default is \code{FALSE}. Set to \code{TRUE} to suppress messages.}

\item{georeferenced}{Default is \code{FALSE}. Set to \code{TRUE} to filter by photos that have geo data.}

\item{progress}{Progress bar is set to \code{TRUE} by default. Use \code{FALSE} especially when calling functions inside documents.}

\item{other_catalog_numbers}{Need to describe these parameters}

\item{foptions}{A list of additional arguments for \code{httr}. There is no reason to use this argument except for debugging purposes.}
}
\description{
Search the photos methods in the ecoengine API.
}
\examples{
# Request all photos. This request will paginate.
# merced <- ee_photos(county = "Merced County")
 ee_photos(page_size = 10)
# Search by collection code. See notes above on options
# ee_photos(collection_code = "CalAcademy")
# ee_photos(collection_code = "VTM")
# ee_photos(collection_code = "CalFlora")
# ee_photos(collection_code = "CDFA")
# Search by county.
# sc_county <- ee_photos(county = "Santa Clara County")
# merced <- ee_photos(county = "Merced County")
# merced <- ee_photos(county = "Merced County", page = "all")
# The package also contains a full list of counties
data(california_counties)
# alameda <- ee_photos(county = california_counties[1, 1])
# alameda$data
# You can also get all the data for Alameda county with one request
# alameda <- ee_photos(county = "Alameda county", page = "all")
# Spidering through the rest of the counties can easily be automated.
# Or by author
# charles_results <- ee_photos(author = "Charles Webber", page = 1:2)
# You can also request all pages in a single call by using ee_photos()
# In this example below, there are 6 pages of results (52 result items).
# Function will return all at once.
# racoons <- ee_photos(scientific_name = "Procyon lotor", page = "all")
}
\seealso{
related: \code{\link{ee_photos}} \code{\link{california_counties}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_checklists.R
\name{ee_checklists}
\alias{ee_checklists}
\title{Ecoengine checklists}
\usage{
ee_checklists(subject = NULL, foptions = list())
}
\arguments{
\item{subject}{Enter one of the following subjects: Mammals, Mosses, Beetles, Spiders, Amphibians, Ants, Fungi, Lichen, Plants.}

\item{foptions}{A list of additional arguments for \code{httr}. There is no reason to use this argument except for debugging purposes.}
}
\value{
data.frame
}
\description{
Retrieves existing checklists from the ecoengine database
}
\examples{
# all_lists  <- ee_checklists()
# mammals_list  <- ee_checklists(subject = "Mammals")
# spiders  <- ee_checklists(subject = "Spiders")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_utils.R
\name{ee_cbind}
\alias{ee_cbind}
\title{ee_cbind}
\usage{
ee_cbind(results)
}
\arguments{
\item{results}{A list of objects of class \code{ecoengine}}
}
\description{
Allows for combining split ecoengine calls (e.g. paginated calls) back into one single result object
}
\examples{
\dontrun{
x1 <- ee_observations(genus = "Lynx", page = 1)
x2 <- ee_observations(genus = "Lynx", page = 2)
x12 <- ee_cbind(list(x1, x2))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_observations.R
\name{ee_observations}
\alias{ee_observations}
\title{Observations List}
\usage{
ee_observations(
  page = NULL,
  page_size = 1000,
  country = "United States",
  state_province = NULL,
  county = NULL,
  kingdom = NULL,
  phylum = NULL,
  order = NULL,
  clss = NULL,
  family = NULL,
  genus = NULL,
  scientific_name = NULL,
  kingdom__exact = NULL,
  phylum__exact = NULL,
  order__exact = NULL,
  clss__exact = NULL,
  family__exact = NULL,
  genus__exact = NULL,
  scientific_name__exact = NULL,
  remote_id = NULL,
  collection_code = NULL,
  source = NULL,
  min_date = NULL,
  max_date = NULL,
  georeferenced = FALSE,
  bbox = NULL,
  exclude = NULL,
  extra = NULL,
  quiet = FALSE,
  progress = FALSE,
  foptions = list()
)
}
\arguments{
\item{page}{Page number}

\item{page_size}{The number of observations per page returned by the query. The package default is \code{250} observations. The API default, however, is \code{10}. Setting this number higher will require fewer calls (i.e. pages) but a large size may also fail or timeout.}

\item{country}{country name}

\item{state_province}{description needed.}

\item{county}{California county. See \code{data(california_counties)}}

\item{kingdom}{kingdom name}

\item{phylum}{phylum name}

\item{order}{order name}

\item{clss}{class name}

\item{family}{family name}

\item{genus}{genus name.}

\item{scientific_name}{A full scientific name}

\item{kingdom__exact}{exact kingdom name}

\item{phylum__exact}{exact phylum name}

\item{order__exact}{exact order name}

\item{clss__exact}{class name}

\item{family__exact}{exact family name}

\item{genus__exact}{exact genus name}

\item{scientific_name__exact}{exact scientific name}

\item{remote_id}{remote ID}

\item{collection_code}{collections code}

\item{source}{data source. See \code{\link{ee_sources}}}

\item{min_date}{Lower date bound. Most be in format Y-m-D h:m:s. This package uses lubridate's \code{ymd_hms} function to format dates.}

\item{max_date}{upper date bound. Most be in format Y-m-D h:m:s. This package uses lubridate's \code{ymd_hms} function to format dates.}

\item{georeferenced}{Default is \code{FALSE}. Set to \code{TRUE} to return only georeferenced records.}

\item{bbox}{Set a bounding box for your search. Use format \code{bbox=-124,32,-114,42}. Order is min Longitude , min Latitude , max Longitude , max Latitude. Use \code{http://boundingbox.klokantech.com/} this website to quickly grab a bounding box (set format to csv on lower right)}

\item{exclude}{Default is \code{NULL}. Pass a list of fields to exclude.}

\item{extra}{Default is \code{NULL}. Pass a list of extra fields to be returned. Additional fields include:
# "id", "record", "source", "remote_resource", "begin_date", "end_date",
# "collection_code", "institution_code", "state_province", "county",
# "last_modified", "original_id", "geometry", "coordinate_uncertainty_in_meters"
# "md5", "scientific_name", "observation_type", "date_precision", "locality"
# "earliest_period_or_lowest_system", "latest_period_or_highest_system", "kingdom"
# "phylum", "clss", "order", "family", "genus", "specific_epithet",
# "infraspecific_epithet", "minimum_depth_in_meters", "maximum_depth_in_meters",
# "maximum_elevation_in_meters", "minimum_elevation_in_meters", "catalog_number"
# "preparations", "sex", "life_stage", "water_body", "country", "individual_count",
"associated_resources"}

\item{quiet}{Default is \code{FALSE}. Set to \code{TRUE} to suppress messages.}

\item{progress}{Progress bar is set to \code{TRUE} by default. Use \code{FALSE} especially when calling functions inside documents.}

\item{foptions}{A list of additional arguments for \code{httr}. There is no reason to use this argument except for debugging purposes.}
}
\value{
\code{data.frame}
}
\description{
API endpoint that represents a list of observations.
}
\examples{
# vulpes <- ee_observations(genus = "vulpes")
\dontrun{
# pinus <- ee_observations(scientific_name = "Pinus", page_size = 100)
# lynx_data <- ee_observations(genus = "Lynx")
# Georeferenced data only
# lynx_data <- ee_observations(genus = "Lynx", georeferenced = TRUE)
# animalia <- ee_observations(kingdom = "Animalia")
# Artemisia <- ee_observations(scientific_name = "Artemisia douglasiana")
# asteraceae <- ee_observationss(family = "asteraceae")
# vulpes <- ee_observations(genus = "vulpes")
# Anas <- ee_observations(scientific_name = "Anas cyanoptera", page = "all")
# loons <- ee_observations(scientific_name = "Gavia immer", page = "all")
# plantae <- ee_observations(kingdom = "plantae")
# chordata <- ee_observations(phylum = "chordata")
# Class is clss since the former is a reserved keyword in SQL.
# aves <- ee_observations(clss = "aves")
# You can also exclude fields using the exclude field
#  aves <- ee_observations(clss = "aves", exclude = "remote_resource")
#  aves <- ee_observations(clss = "aves", exclude = "source,remote_resource")
#  or request additional fields using extra
#  aves <- ee_observations(clss = "aves", extra = "kingdom,genus")
# aves <- ee_observations(clss = "aves", bbox = '-124,32,-114,42')
# aves <- ee_observations(clss = "aves", county = "Alameda county")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ee_utils.R
\name{print.ecoengine}
\alias{print.ecoengine}
\title{Print a summary for an ecoengine object}
\usage{
\method{print}{ecoengine}(x, ...)
}
\arguments{
\item{x}{An object of class \code{ecoengine}}

\item{...}{additional arguments}
}
\description{
Print a summary for an ecoengine object
}
