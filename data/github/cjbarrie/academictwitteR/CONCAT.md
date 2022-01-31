---
title: 'academictwitteR: an R package to access the Twitter Academic Research Product Track v2 API endpoint'
tags:
  - R
  - twitter
  - social media
  - API
authors:
  - name: Christopher Barrie
    orcid: 0000-0002-9156-990X
    affiliation: 1
  - name: Justin Chun-ting Ho
    orcid: 0000-0002-7884-1059
    affiliation: 2
affiliations:
 - name: School of Social and Political Sciences, University of Edinburgh, Scotland, UK.
   index: 1
 - name: Centre for European Studies and Comparative Politics, Sciences Po, France.
   index: 2
date: 23 April 2021
bibliography: paper.bib
---


# Statement of need

In January, 2021, Twitter announced the "Academic Research Product Track." This provides academic researchers with greatly expanded access to Twitter data. Existing R packages for querying the Twitter API, such as the popular ``rtweet`` package [@rtweet], are yet to introduce functionality to allow users to connect to the new v2 API endpoints with Academic Research Product Track credentials. The ``academictwitteR`` package [@academictwitteR] is built with academic research in mind. It encourages efficient and responsible storage of data, given the likely large amounts of data being collected, as well as a number of shortcut and query building functions to access new v2 API endpoints.

# Summary

The Twitter Application Programming Interface, or API, was first introduced in 2006. It was designed principally with commercial objectives in mind. Over time, however, researchers began to repurpose the Twitter API for academic ends. In January, 2021, [Twitter announced the "Academic Research Product Track"](https://blog.twitter.com/developer/en_us/topics/tools/2021/enabling-the-future-of-academic-research-with-the-twitter-api.html), noting that "[t]oday, academic researchers are one of the largest groups of people using the Twitter API."

Authorization for the Academic Research Product Track provides access to the Twitter v2 API endpoints, introduced in 2020, as well as much improved data access. In summary the Academic Research product track allows the authorized user:

1. Access to the full archive of (as-yet-undeleted) tweets published on Twitter;
2. A higher monthly tweet cap (10m---or 20x what was previously possible with the standard v1.1 API);
3. Ability to access these data with more precise filters permitted by the v2 API.

The ``academictwitteR`` package was designed: 1) to make the Academic Research Product Track easily accessible for R users by providing dedicated functions to query the the v2 API endpoints; 2) to encourage academic researchers efficiently and safely to store their data.

<<<<<<< HEAD
The functions allow the user to collect tweets from (or to) specified users and to collect tweets containing specified words or sets of words. In particular, queries that include so-called "conjunction-required"" operators can also be accessed via a set of shortcut functions for accessing e.g. tweets containing media content, tweets containing geographic location information, or tweets containing urls. Additionally, separate query builder functions allow the user to specify complex queries to incorporate into the API call. 
=======
The functions allow the user to collect tweets from (or to) specified users and to collect tweets containing specified words or sets of words. In particular, queries that include so-called "conjunction-required" operators can also be accessed via a set of shortcut functions for accessing e.g. tweets containing media content, tweets containing geographic location information, or tweets containing urls. Additionally, separate query builder functions allow the user to specify complex queries to incorporate into the API call. 
>>>>>>> a55475f5d8ea852e5173b04c3eec68a35762d0b8

Data is stored in serialized form as RDS files or as separate JSON files. The former represents the most efficient storage solution for native R data-file formats; the latter helps mitigate loss by storing data as separate JSONs for each pagination token (or up to 500 tweets). Convenience functions are also included to bind tweet- and user-level information stored as JSON files, and to pick up data collection where it left off in the case of unplanned interruption.

# References


# academictwitteR <img src="man/figures/academictwitteRhex.png" width="160px" align="right" />

<!-- badges: start -->
[![v2](https://img.shields.io/endpoint?url=https%3A%2F%2Ftwbadges.glitch.me%2Fbadges%2Fv2)](https://developer.twitter.com/en/docs/twitter-api)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03272/status.svg)](https://doi.org/10.21105/joss.03272) 
[![](https://www.r-pkg.org/badges/version/academictwitteR)](https://cran.r-project.org/package=academictwitteR)
![Downloads](https://cranlogs.r-pkg.org/badges/academictwitteR)
[![](http://cranlogs.r-pkg.org/badges/grand-total/academictwitteR)](https://cran.r-project.org/package=academictwitteR)
[![Codecov test coverage](https://codecov.io/gh/cjbarrie/academictwitteR/branch/master/graph/badge.svg)]( https://app.codecov.io/gh/cjbarrie/academictwitteR?branch=master)
<!-- badges: end -->


[![Twitter](https://img.shields.io/twitter/url/https/twitter.com/cbarrie.svg?style=social&label=Follow%20%40cbarrie)](https://twitter.com/cbarrie)
[![Twitter](https://img.shields.io/twitter/url/https/twitter.com/justin_ct_ho.svg?style=social&label=Follow%20%40justin_ct_ho)](https://twitter.com/justin_ct_ho)

Repo containing code to for R package <tt>academictwitteR</tt> to collect tweets from v2 API endpoint for the Academic Research Product Track.

To cite package ‘academictwitteR’ in publications use:

  - Barrie, Christopher and Ho, Justin Chun-ting. (2021). academictwitteR: an R package to access the Twitter Academic Research Product Track v2 API endpoint. *Journal of Open Source Software*, 6(62), 3272, https://doi.org/10.21105/joss.03272

A BibTeX entry for LaTeX users is:

```
@article{BarrieHo2021,
  doi = {10.21105/joss.03272},
  url = {https://doi.org/10.21105/joss.03272},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {62},
  pages = {3272},
  author = {Christopher Barrie and Justin Chun-ting Ho},
  title = {academictwitteR: an R package to access the Twitter Academic Research Product Track v2 API endpoint},
  journal = {Journal of Open Source Software}
}

  
```

## Installation

You can install the package with:
``` r
install.packages("academictwitteR")
```

Alternatively, you can install the development version with:
``` r
devtools::install_github("cjbarrie/academictwitteR", build_vignettes = TRUE)
```

Get started by reading `vignette("academictwitteR-intro")`.

To use the package, it first needs to be loaded with:

```r

library(academictwitteR)

```

The <tt>academictwitteR</tt> package has been designed with the efficient storage of data in mind. Queries to the API include arguments to specify whether tweets be stored as a .rds file using the `file` argument or as separate JSON files for tweet- and user-level information separately with argument `data_path`.

Tweets are returned as a data.frame object and, when a `file` argument has been included, will also be saved as a .rds file.

When collecting large amounts of data, we recommend the workflow described below, which allows the user : 1) to efficiently store authorization credentials; 2) to efficiently store returned data; 3) bind the data into a data.frame object or tibble ;4) resume collection in case of interruption; and 5) update collection in case of need.

## Authorization

The first task is set authorization credentials with the `set_bearer()` function, which allows the user to store their bearer token in the .Renviron file.

To do so, use:

```r
set_bearer()
```

and enter authorization credentials as below:

![](vignettes/files/TWITTER_BEARER.gif)

This will mean that the bearer token is automatically called during API calls. It also avoids the inadvisable practice of hard-coding authorization credentials into scripts. 

See the vignette documentation `vignette("academictwitteR-auth")` for further information on obtaining a bearer token.

## Collection

The workhorse function is `get_all_tweets()`, which is able to collect tweets matching a specific search query or all tweets by a specific set of users.

```r

tweets <-
  get_all_tweets(
    query = "#BlackLivesMatter",
    start_tweets = "2020-01-01T00:00:00Z",
    end_tweets = "2020-01-05T00:00:00Z",
    file = "blmtweets",
    data_path = "data/",
    n = 1000000,
  )
  
```

Here, we are collecting tweets containing a hashtag related to the Black Lives Matter movement over the period January 1, 2020 to January 5, 2020. 

We have also set an upper limit of one million tweets. When collecting large amounts of Twitter data we recommend including a `data_path` and setting `bind_tweets = FALSE` such that data is stored as JSON files and can be bound at a later stage upon completion of the API query.

```r

tweets <-
  get_all_tweets(
    users = c("jack", "cbarrie"),
    start_tweets = "2020-01-01T00:00:00Z",
    end_tweets = "2020-01-05T00:00:00Z",
    file = "blmtweets",
    n = 1000
  )
  
```

Whereas here we are not specifying a search query and instead are requesting all tweets by users @jack and @cbarrie over the period January 1, 2020 to January 5, 2020. Here, we set an upper limit of 1000 tweets.

The search query and user query arguments can be combined in a single API call as so:

```r

get_all_tweets(
  query = "twitter",
  users = c("cbarrie", "jack"),
  start_tweets = "2020-01-01T00:00:00Z",
  end_tweets = "2020-05-01T00:00:00Z",
  n = 1000
)

```

Where here we would be collecting tweets by users @jack and @cbarrie over the period January 1, 2020 to January 5, 2020 containing the word "twitter."

```r

get_all_tweets(
  query = c("twitter", "social"),
  users = c("cbarrie", "jack"),
  start_tweets = "2020-01-01T00:00:00Z",
  end_tweets = "2020-05-01T00:00:00Z",
  n = 1000
)

```

While here we are collecting tweets by users @jack and @cbarrie over the period January 1, 2020 to January 5, 2020 containing the words "twitter" or "social."

Note that the "AND" operator is implicit when specifying more than one character string in the query. See [here](https://developer.twitter.com/en/docs/twitter-api/tweets/search/integrate/build-a-query) for information on building queries for search tweets. Thus, when searching for all elements of a character string, a call may look like:

```r

get_all_tweets(
  query = c("twitter social"),
  users = c("cbarrie", "jack"),
  start_tweets = "2020-01-01T00:00:00Z",
  end_tweets = "2020-05-01T00:00:00Z",
  n = 1000
)

```

, which will capture tweets containing *both* the words "twitter" and "social." The same logics apply for hashtag queries.

Whereas if we specify our query as separate elements of a character vector like this:

```r

get_all_tweets(
  query = c("twitter", "social"),
  users = c("cbarrie", "jack"),
  start_tweets = "2020-01-01T00:00:00Z",
  end_tweets = "2020-05-01T00:00:00Z",
  n = 1000
)

```
, this will be capturing tweets by users @cbarrie or @jack containing the words "twitter" *or* social. 

Finally, we may wish to query an exact phrase. To do so, we can either input the phrase in escape quotes, e.g., `query ="\"Black Lives Matter\""` or we can use the optional parameter `exact_phrase = T` (in devt. version) to search for tweets containing the exact phrase string:

```r

tweets <-
  get_all_tweets(
    query = "Black Lives Matter",
    exact_phrase = T,
    start_tweets = "2021-01-04T00:00:00Z",
    end_tweets = "2021-01-04T00:45:00Z",
    n = Inf
  )

```

See the vignette documentation `vignette("academictwitteR-build")` for further information on building more complex API calls.

## Data storage

Files are stores as JSON files in specified directory when a `data_path` is specified. Tweet-level data is stored in files beginning "data_"; user-level data is stored in files beginning "users_".

If a filename is supplied, the functions will save the resulting tweet-level information as a .rds file.

Functions always return a data.frame object unless a `data_path` is specified and `bind_tweets` is set to `FALSE`. When collecting large amounts of data, we recommend using the `data_path` option with `bind_tweets = FALSE`. This mitigates potential data loss in case the query is interrupted. 

See the vignette documentation `vignette("academictwitteR-intro")` for further information on data storage conventions.

## Reformatting

Users can then use the `bind_tweets` convenience function to bundle the JSONs into a data.frame object for analysis in R as such:

```r
tweets <- bind_tweets(data_path = "data/")
users <- bind_tweets(data_path = "data/", user = TRUE)
```

To bind JSONs into tidy format, users can also specify a tidy output format. 

```r
bind_tweets(data_path = "tweetdata", output_format = "tidy")
```

See the vignette documentation `vignette("academictwitteR-tidy")` for further information on alternative output formats.

## Interruption and Continuation

The package offers two functions to deal with interruption and continue previous data collection session. If you have set a data_path and export_query was set to "TRUE" during the original collection, you can use `resume_collection()` to resume a previous interrupted collection session. An example would be:

```r
resume_collection(data_path = "data")
```

If a previous data collection session is completed, you can use `update_collection()` to continue data collection with a new end date. This function is particularly useful for getting data for ongoing events. An example would be:

```r
update_collection(data_path = "data", end_tweets = "2020-05-10T00:00:00Z")
```

## Note on v2 Twitter API

For more information on the parameters and fields available from the v2 Twitter API endpoint see: [https://developer.twitter.com/en/docs/twitter-api/tweets/search/api-reference/get-tweets-search-all](https://developer.twitter.com/en/docs/twitter-api/tweets/search/api-reference/get-tweets-search-all).

## Arguments

`get_all_tweets()` accepts a range of arguments, which can be combined to generate a more precise query.

| Arguments   |     Description      |
|----------|:-------------:|
|query | Search query or queries e.g. "cat"
|exact_phrase | If `TRUE`, only tweets will be returned matching the exact phrase
|users | string or character vector, user handles to collect tweets from the specified users
|reply_to | string or character vector, user handles to collect replies to the specified users
|retweets_of| string or character vector, user handles to collects retweets of tweets by the specified users
|exclude | string or character vector, tweets containing the keyword(s) will be excluded
|is_retweet | If `TRUE`, only retweets will be returned; if `FALSE`, retweets will not be returned, only tweets will be returned; if `NULL`, both retweets and tweets will be returned.
|is_reply | If `TRUE`, only reply tweets will be returned
|is_quote | If `TRUE`, only quote tweets will be returned
|is_verified |If `TRUE`, only tweets whose authors are verified by Twitter will be returned
|remove_promoted | If `TRUE`, tweets created for promotion only on ads.twitter.com are removed
|has_hashtags | If `TRUE`, only tweets containing hashtags will be returned
|has_cashtags | If `TRUE`, only tweets containing cashtags will be returned
|has_links | If `TRUE`, only tweets containing links and media will be returned
|has_mentions |If `TRUE`, only tweets containing mentions will be returned
|has_media |If `TRUE`, only tweets containing a recognized media object, such as a photo, GIF, or video, as determined by Twitter will be returned
|has_images |If `TRUE`, only tweets containing a recognized URL to an image will be returned
|has_videos |If `TRUE`, only tweets containing contain native Twitter videos, uploaded directly to Twitter will be returned
|has_geo |If `TRUE`, only tweets containing Tweet-specific geolocation data provided by the Twitter user will be returned
|place | Name of place e.g. "London"
|country | Name of country as ISO alpha-2 code e.g. "GB"
|point_radius | A vector of two point coordinates latitude, longitude, and point radius distance (in miles)
|bbox | A vector of four bounding box coordinates from west longitude to north latitude
|lang | A single BCP 47 language identifier e.g. "fr"
|url | string, return tweets containing specified url
|conversation_id| string, return tweets that share the specified conversation ID

## Batch Compliance

There are three functions to work with Twitter's Batch Compliance endpoints: `create_compliance_job()` creates a new compliance job and upload the dataset; `list_compliance_jobs` lists all created jobs and their job status; `get_compliance_result()` downloads the result.

## Acknowledgements

Function originally inspired by [Gist](https://gist.github.com/schochastics/1ff42c0211916d73fc98ba8ad0dcb261#file-get_tweets-r-L14) from [https://github.com/schochastics](https://github.com/schochastics).


## Code of Conduct

Please note that the academictwitteR project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.
# academictwitteR 0.3.0
* Added support for batch compliance
* Added functions `get_user_id`, `get_retweeted_by`, and `convert_json`.
* Added parameter `exact_phrase` for `build_query` (also for the downstream function `get_all_tweets`).

# academictwitteR 0.2.1
* Fixed error 400 when fetching tweets with the context annotation field

# academictwitteR 0.2.0

* Support Likes lookup, followers, following, liked tweets, and liking user endpoints
* A function to counts tweets by query string: `count_all_tweets`
* Added the `n` argument
* Autosleep when hitting rate limit
* `bind_tweets` allows transforming the collected tweets to various formats, e.g. tidy data frame
* `set_bearer` and `get_bearer` for managing the bearer token
* Many wrappers to `get_all_tweets` are deprecated
* Bug fixes and efficiency improvements

# academictwitteR 0.1.0
* Initial CRAN version.
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

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
* Focusing on what is best not just for us as individuals, but for the overall
community

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

Community leaders are responsible for clarifying and enforcing our standards
of acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies
when an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at christopher.barrie@ed.ac.uk or justin.chunting.ho@sciencespo.fr. All complaints will be reviewed and investigated promptly and fairly.

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

**Community Impact**: A violation through a single incident or series of
actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or permanent
ban.

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
standards, including sustained inappropriate behavior, harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within the
community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0,
available at https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at https://
www.contributor-covenant.org/translations.
## Test environments
* MACOS 10.15.6, R 4.0.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs. # visual output

    Code
      capture_warnings(x <- get_all_tweets(query = "#commtwitter", start_tweets = "2021-06-01T00:00:00Z",
        end_tweets = "2021-06-05T00:00:00Z", context_annotations = FALSE, verbose = TRUE))
    Output
      query:  #commtwitter 
      Total pages queried: 1 (tweets captured this page: 5).
      This is the last page for #commtwitter : finishing collection.
      [1] "Recommended to specify a data path in order to mitigate data loss when ingesting large amounts of data."                  
      [2] "Tweets will not be stored as JSONs or as a .rds file and will only be available in local memory if assigned to an object."
    Code
      capture_warnings(y <- get_all_tweets(query = "#commtwitter", start_tweets = "2021-06-01T00:00:00Z",
        end_tweets = "2021-06-05T00:00:00Z", context_annotations = TRUE, verbose = TRUE))
    Output
      page_n is limited to 100 due to the restriction imposed by Twitter API
      query:  #commtwitter 
      Total pages queried: 1 (tweets captured this page: 5).
      This is the last page for #commtwitter : finishing collection.
      [1] "Recommended to specify a data path in order to mitigate data loss when ingesting large amounts of data."                  
      [2] "Tweets will not be stored as JSONs or as a .rds file and will only be available in local memory if assigned to an object."
    Code
      capture_warnings(z <- get_all_tweets(query = "#commtwitter", start_tweets = "2021-06-01T00:00:00Z",
        end_tweets = "2021-06-05T00:00:00Z", context_annotations = TRUE, page_n = 99,
        verbose = TRUE))
    Output
      query:  #commtwitter 
      Total pages queried: 1 (tweets captured this page: 5).
      This is the last page for #commtwitter : finishing collection.
      [1] "Recommended to specify a data path in order to mitigate data loss when ingesting large amounts of data."                  
      [2] "Tweets will not be stored as JSONs or as a .rds file and will only be available in local memory if assigned to an object."
    Code
      capture_warnings(x1 <- get_all_tweets(query = "#commtwitter", start_tweets = "2021-06-01T00:00:00Z",
        end_tweets = "2021-06-05T00:00:00Z", context_annotations = FALSE, verbose = FALSE))
    Output
      character(0)
    Code
      capture_warnings(y1 <- get_all_tweets(query = "#commtwitter", start_tweets = "2021-06-01T00:00:00Z",
        end_tweets = "2021-06-05T00:00:00Z", context_annotations = TRUE, verbose = FALSE))
    Output
      character(0)
    Code
      capture_warnings(z1 <- get_all_tweets(query = "#commtwitter", start_tweets = "2021-06-01T00:00:00Z",
        end_tweets = "2021-06-05T00:00:00Z", context_annotations = TRUE, page_n = 99,
        verbose = FALSE))
    Output
      character(0)

# running two make_query in rapid succession will not trigger HTTP 429

    Code
      academictwitteR:::make_query(url = endpoint_url, params = params, bearer_token = get_bearer())
    Output
      $meta
      $meta$result_count
      [1] 0
      
      

---

    Code
      academictwitteR:::make_query(url = endpoint_url, params = params, bearer_token = get_bearer())
    Output
      $meta
      $meta$result_count
      [1] 0
      
      

# test .trigger_sleep

    Code
      academictwitteR:::.trigger_sleep(r, ref_time = as.POSIXlt("2021-06-25 11:50:30",
        tz = "UTC"), verbose = TRUE, really_sleep = FALSE, tzone = "UTC")
    Output
      Rate limit reached. Rate limit will reset at 2021-06-25 11:50:37 
      Sleeping for 8 seconds. 
      ================================================================================

---

    Code
      academictwitteR:::.trigger_sleep(r, ref_time = as.POSIXlt("2021-06-25 11:50:30",
        tz = "UTC"), verbose = FALSE, really_sleep = FALSE, tzone = "UTC")

---

    Code
      academictwitteR:::.trigger_sleep(r, ref_time = as.POSIXlt("2021-06-25 11:50:33",
        tz = "UTC"), verbose = TRUE, really_sleep = TRUE, tzone = "UTC")
    Output
      Rate limit reached. Rate limit will reset at 2021-06-25 11:50:37 
      Sleeping for 5 seconds. 
      ================================================================================

---

    Code
      academictwitteR:::.trigger_sleep(r, ref_time = as.POSIXlt("2021-06-25 11:50:36",
        tz = "UTC"), verbose = FALSE, really_sleep = TRUE, tzone = "UTC")

# disable adaptive sleeping if ref_time is later than reset_time in r, #213

    Code
      academictwitteR:::.trigger_sleep(r, ref_time = as.POSIXlt("2021-06-25 11:50:40",
        tz = "UTC"), verbose = TRUE, really_sleep = FALSE, tzone = "UTC")
    Output
      Rate limit reached. Cannot estimate adaptive sleep time. Sleeping for 900 seconds. 
      ================================================================================

---

    Code
      academictwitteR:::.trigger_sleep(r, ref_time = as.POSIXlt("2021-06-25 11:50:40",
        tz = "UTC"), verbose = FALSE, really_sleep = FALSE, tzone = "UTC")

---
title: "Intro. to academictwitteR"
author: Christopher Barrie and Justin Ho
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Intro. to academictwitteR}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Introduction

This vignette provides an introduction to the R package `academictwitteR`. The package is useful solely for querying the Twitter Academic Research Product Track v2. API endpoint. 

This version of the Twitter API allows researchers to access larger volumes of Twitter data. For more information on the the Twitter API, including how to apply for access to the Academic Research Product Track, see the Twitter Developer platform.

The following vignette will guide you through how to use the package.

We will begin by describing the thinking behind the development of this package and, specifically, the data storage conventions we have established when querying the API.

## The Twitter Academic Research Product Track

The Academic Research Product Track permits the user to access larger volumes of data, over a far longer time range, than was previously possible. From the Twitter release for the new track:

> "The Academic Research product track includes full-archive search, as well as increased access and other v2 endpoints and functionality designed to get more precise and complete data for analyzing the public conversation, at no cost for qualifying researchers. Since the Academic Research track includes specialized, greater levels of access, it is reserved solely for non-commercial use".

The new "v2 endpoints" refer to the v2 API, introduced around the same time as the new Academic Research Product Track. Full details of the v2 endpoints are available on the Twitter Developer platform.

In summary the Academic Research product track allows the authorized user:

1. Access to the full archive of (as-yet-undeleted) tweets published on Twitter
2. A higher monthly tweet cap (10m--or 20x what was previously possible with the standard v1.1 API)
3. Ability to access these data with more precise filters permitted by the v2 API

## Setting up your bearer token

Please refer to [this vignette](academictwitteR-auth.html) on how to obtain your own bearer token. You can supply this bearer token in every request. The more advisable and secure approach is to set up your bearer token in your .Renviron file.

We begin by loading the package with:

```{r setup}
library(academictwitteR)
```

And then launch `set_bearer`. This will open your .Renviron file in the home directory. Enter your bearer token as below (the bearer token used below is not real).

![](files/TWITTER_BEARER.gif){width=75%}

**For this environment variable to be recognized, you first have to restart R** 

You can then obtain your bearer token with `get_bearer`. This is also the default for all data collection functions.

You can check that this works with: 

```{r, eval = FALSE}
get_bearer()
```

```{r, echo=F}
dummy_bearer <- "AAAAAAAAAAAAAAAAAAAAAPwXWFFlLLDVC6G0PFo4shkDVg02DwVxGQIVKvhPVE3vdV"
dummy_bearer

```

## Querying the Twitter API with `academictwitteR`

The workhorse function of `academictwitteR` for collecting tweets is `get_all_tweets()`.

```{r, eval=F}

tweets <-
  get_all_tweets(
    query = "#BlackLivesMatter",
    start_tweets = "2020-01-01T00:00:00Z",
    end_tweets = "2020-01-05T00:00:00Z",
    file = "blmtweets"
  )
  
```

Here, we are collecting tweets containing a hashtag related to the Black Lives Matter movement over the period January 1, 2020 to January 5, 2020. 

Note that once we have stored our bearer token with set_bearer, it will be called within the function automatically.

This query will only capture a maximum of 100 tweets as we have not changed the default 

If you have not set your bearer token, you can do also do so within the function as follows:

```{r, eval=F}

tweets <-
  get_all_tweets(
    query = "#BlackLivesMatter",
    start_tweets = "2020-01-01T00:00:00Z",
    end_tweets = "2020-01-05T00:00:00Z",
    bearer_token = "AAAAAAAAAAAAAAAAAAAAAPwXWFFlLLDVC6G0Pg02DwVxGQIVKTHISISNOTAREALTOKEN",
    file = "blmtweets"
  )
  
```

This is **not** recommended as by convention it is not advisable to keep API authorization tokens within your scripts. 

## Storage conventions in `academictwitteR`

Given the sizeable increase in the volume of data potentially retrievable with the Academic Research Product Track, it is advisable that researchers establish clear storage conventions to mitigate data loss caused by e.g. the unplanned interruption of an API query.

We first draw your attention first to the `file` argument in the code for the API query above.

In the file path, the user can specify the name of a file to be stored with a ".rds" extension, which includes all of the tweet-level information collected for a given query.

Alternatively, the user can specify a `data_path` as follows:

```{r, eval=F}

tweets <-
  get_all_tweets(
    query = "#BlackLivesMatter",
    start_tweets = "2015-01-01T00:00:00Z",
    end_tweets = "2020-01-05T00:00:00Z",
    data_path = "data/",
    bind_tweets = FALSE,
    n= 1000000
  )
  
```

In the data path, the user can either specify a directory that already exists or name a new directory.

The data is stored in this folder as a series of JSONs. Tweet-level data is stored as a series of JSONs beginning "data_"; User-level data is stored as a series of JSONs beginning "users_".

Note that the `get_all_tweets()` function always returns a data.frame object unless `data_path` is specified and `bind_tweets` is set to `FALSE`. 

When collecting large amounts of data, we recommend using the `data_path` option with `bind_tweets = FALSE`. This mitigates potential data loss in case the query is interrupted, and avoids system memory usage errors.

Note finally that here we are setting an upper limit of tweets of one million. The default limit is set to 100. **For almost all applications, users will wish to change this**. We can also set `n = Inf` if we do not require any upper limit. This will collect all available tweets matching the query. 

## Binding JSON files into data.frame objects

When `bind_tweets` is `FALSE`, no data.frame object will be returned. In order to get the tweets into a data.frame, you can then use the `bind_tweets()` helper function to bundle the JSONs into a data.frame object for analysis in R as such:

```{r, eval=FALSE}
tweets <- bind_tweets(data_path = "data/")
```

If you want to bundle together the user-level data, you can achieve this with the same helper function. The only change is that `user` is now set to `TRUE`, meaning we want to bundle user-level data:

```{r, eval=FALSE}
users <- bind_tweets(data_path = "data/", user = TRUE)
```

Note: v0.2 of the package incorporates functionality to convert JSONs into multiple data frame formats. Most usefully, these additions permit the incorporation of user-level and tweet-level data into a single tibble.---
title: "Understanding API errors"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Understanding API errors}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Introduction

The table below gives an overview of the common errors you may encounter when using the Twitter Academic Research Product Track API. The meaning and solutions for these errors are taken from the Twitter API Response codes and errors help page [here](https://developer.twitter.com/en/support/twitter-api/error-troubleshooting).

The most common error that is fixable is a 400 status code error. This means that the query has been misspecified. In these cases, return to the query and consult the documentation to ensure that it has been appropriately specified.


| Error | Text                  | Meaning                                                                                                                                                                                                                          | Solution                                                                                                                                                                                                                                                                                                            |
| ----- | --------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 400   | Bad request           | The request was invalid or cannot be otherwise served. An accompanying error message will explain further. Requests without authentication or with invalid query parameters are considered invalid and will yield this response. | Double check that your query is valid.                                                                                                                                                                                                                                                                              |
| 401   | Unauthorized          | There was a problem authenticating your request. This could be due to missing or incorrect authentication credentials. This may also be returned in other undefined circumstances.                                               | Check that you are using the correct authentication method and that your credentials are correct. Consult the authorization vignette [here](academictwitteR-auth.html) on how to get authorization for using the Academic Research Product Track and make sure you have correctly specified your bearer token. See \`?get\_bearer\` for more details. |
| 403   | Forbidden             | The request is understood, but it has been refused or access is not allowed. An accompanying error message will explain why.                                                                                                     | Check that your data plan includes access to the endpoint you’re trying to use. You may also need to get your App allowlisted.                                                                                                                                                                                      |
| 404   | Not Found             | The URI requested is invalid or the resource requested, such as a user, does not exist.                                                                                                                                          | Check that you are using valid parameters and the correct URI for the endpoint you’re using.                                                                                                                                                                                                                        |
| 429   | Too Many Requests     | Returned when a request cannot be served due to the App's rate limit having been exhausted for the resource.                                                                                                                     | Check the number of requests per timeframe allowed with the endpoint you’re using. Wait for the timeframe to reset.                                                                                                                                                                                                 |
| 500   | Internal Server Error | Something is broken. This is usually a temporary error, for example in a high load situation or if an endpoint is temporarily having issues.                                                                                     | Check the Twitter API status page or the developer community forum in case others are having similar issues, or simply wait and try again later.                                                                                                                                                                    |
| 503   | Service Unavailable   | The Twitter servers are up, but overloaded with requests. Try again later.                                                                                                                                                       | Check the Twitter API status page or the developer community forum in case others are having similar issues, or simply wait and try again later.                                                                                                                                                                    |                              |---
title: "Batch Compliance"
author: Christopher Barrie, Justin Ho, Chung-hong Chan
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Batch Compliance}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Introduction

It is crucial for Twitter data holders (e.g. users of academictwitteR) to ensure that the data reflects user intent and the current state of content on Twitter. Specifically, if a tweet is deleted or modified on Twitter, data holders must delete or modify any content stored offline accordingly. The batch compliance endpoints and corresponding functions offer researchers an easy way to help maintain Twitter data in compliance with the [Twitter Developer Agreement and Policy](https://developer.twitter.com/en/developer-terms/policy).

## Upload Dataset

To create a compliance job and upload a dataset, you can use the `create_compliance_job()` function:

There are two ways to use `create_compliance_jobs`. You can supply a vector of tweet/user ids.

```r
jobid <- create_compliance_job(x = c("1233905174730682369", "1233905174659256320",
                                     "1233905174650986497"),
                               type = "tweets")
jobid
```

Another way is to provide a text file, in which each line contains a Tweet ID or user ID. The content of a text file looks like this:

```
1233905174730682369
1233905174659256320
1233905174650986497
1233905174655139841
1233905174663458816
1233905174835449856
1233905174789443584
1233905174533545987
1233905174936080384
```

A sample can be found [here](https://raw.githubusercontent.com/echen102/COVID-19-TweetIDs/master/2020-03/coronavirus-tweet-id-2020-03-01-00.txt). And then run `create_compliance_job` with the file name (e.g. "tweet_ids.txt").

```{r, eval = FALSE}
jobid <- create_compliance_job(x = "tweet_ids.txt",
                               type = "tweets")
jobid
```

The function will return a Job ID, which will be used to download the results.

## Check Job Status

Optionally, you can check your job status by running:

```{r, eval = FALSE}
list_compliance_jobs()
```

## Download Results
When a job is completed, you can download the compliance results using:

```{r, eval = FALSE}
get_compliance_result(id = "1460077048991555585")
```

The function will automatically check if the job has been completed.

If you forgot your job id, you can retrieve it using `list_compliance_jobs()`.
---
title: "Building a query in academictwitteR"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Building a query in academictwitteR}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Introduction

The v2 Twitter API allows for greater precision when making queries. A query might just be a single string like "happy new year" if you're interested on how people are celebrating on the night of December 31. Alternatively, the query might involve several additional operators that filter tweets with greater precision to return specific tweet content. 

This vignette guides you through the logics underpinning queries to the Twitter API. For full information on these logics you may additionally wish to consult the Twitter API documentation on how to build a query [here](https://developer.twitter.com/en/docs/twitter-api/tweets/search/integrate/build-a-query).

## Query strings

We first load our package into memory with:

```{r, eval =FALSE}
library(academictwitteR)
```

```{r, echo =FALSE}
library(devtools)
load_all()
```

We then make sure we have set our bearer token appropriately by calling:

```{r, eval =FALSE}
get_bearer()
```

```{r, echo =FALSE}
bearer_example <- "AAAAAAAAAAAAAAAAAAAAAPw%2BJQEAAAAAq5Ot8BBYyYlAqT9nLMuVuR1jI5fA%3DqG9HTHISISNOTAREALTOKEN"
bearer_example
```

If your bearer token is not set appropriately, consult the guide in [the authorization vignette](academictwitteR-auth.html).

Let's say we were interested in what people were talking about on New Year's Eve. We might do something like this:

```{r, eval = F}

tweets <-
  get_all_tweets(
    query = "happy",
    start_tweets = "2019-12-31T10:00:00Z",
    end_tweets = "2020-01-01T10:00:00Z",
    n = 10000
  )

```

Note here that we have also specified an upper limit of 10,000 tweets. The default is 100. For most applications, the user will need to specify a higher n than the default. 

The default upper limit is set to 100 in order to prevent unnecessary ingests of data when e.g. trialling an API call.

As an alternative to this, the user might also wish to use the `count_all_tweets()` function in order to get an idea of how many tweets match the specified API query.

## Additional parameters

In the above we search for all tweets between two dates that contain the string "happy." But what if we were only interested in a particular region or written in a particular language?

Let's say we were only interested in tweets written in English and originating from the US. We would add several operators to our query to filter by these characteristics:

```{r, eval = F}

tweets <-
  get_all_tweets(
    query = "happy",
    start_tweets = "2019-12-31T10:00:00Z",
    end_tweets = "2020-01-01T10:00:00Z",
    country = "US", 
    lang = "en"
  )

```

In fact, the `get_all_tweets()` function can be combined with multiple additional filtering parameters. The example below includes numerous additional filters, keeping only tweets with images, hashtags, and mentions:

```{r, eval = F}

tweets <-
  get_all_tweets(
    query = "happy",
    start_tweets = "2019-12-31T10:00:00Z",
    end_tweets = "2020-01-01T10:00:00Z",
    country = "US", 
    lang = "en",
    has_images = TRUE,
    has_hashtags = TRUE,
    has_mentions = TRUE
  )

```

We might then decide that our geo filter is not accurate enough. We don't just want tweets originating from the US but we want tweets from Seattle in particular. This would mean adding more operators to our query:

```{r, eval=F}
tweets <-
  get_all_tweets(
    query = "happy",
    start_tweets = "2019-12-31T10:00:00Z",
    end_tweets = "2020-01-01T10:00:00Z",
    country = "US", 
    place = "seattle",
    lang = "en",
    has_images = TRUE,
    has_hashtags = TRUE,
    has_mentions = TRUE
  )
  
```

What if we were unsatisfied with the accuracy of our geo parameters and we wanted to be sure that our tweets were actually coming from a particular place? Let's say we are interested in central Seattle, as shown in the map below.

![](files/seattle.png){width=70%}

Twitter also allows us to query tweets originating from within a particular geographical buffer too. Here, we simply specify the longitude and latitude of the southwest and then the northeast corners of this bounding box. Note, this image is taken from a screenshot of the website [http://bboxfinder.com](http://bboxfinder.com). 

Many such websites exist that allow you to find the bounding box coordinates of a place of interest, including [https://www.openstreetmap.org](https://www.openstreetmap.org) and [https://boundingbox.klokantech.com/](https://boundingbox.klokantech.com/).

We can then input this information with the Twitter "bounding_box" operator using the `bbox` argument as below:


```{r, eval=F}

tweets <-
  get_all_tweets(
    query = "happy",
    start_tweets = "2019-12-31T10:00:00Z",
    end_tweets = "2020-01-01T10:00:00Z",
    country = "US", 
    place = "seattle",
    lang = "en",
    has_images = TRUE,
    has_hashtags = TRUE,
    has_mentions = TRUE,
    bbox = c(-122.375679, 47.563554, -122.266159, 47.643417)
  )
  
  
```

The alternative `point_radius` argument requires three pieces of information: the longitude and latitude of a target coordinate, and the buffer size around that coordinate.

```{r, eval=F}

tweets <-
  get_all_tweets(
    query = "happy",
    start_tweets = "2019-12-31T10:00:00Z",
    end_tweets = "2020-01-01T10:00:00Z",
    country = "US", 
    place = "seattle",
    lang = "en",
    point_radius = c(-122.33795253639994, 47.60900846404393, 25)
  )
  
```

Note that the maximum radius for the buffer is 25 miles. Similarly, the maximum height and width of any bounding box is 25 miles. Inputting coordinate information that exceeds these bounds will result in a 400 status code error.

## Multiple query strings

In the above we were specifying just one query string. If we were interested in multiple query strings, however, we could easily do this by just defining a character vector and passing this to the `get_all_tweets()` function as our query.

So let's say we wanted to search for tweets with the words "happy," "new" or "year" or any combination of these words within the same tweet. 

```{r, eval = F}

tweets <-
  get_all_tweets(
    query = c("happy", "new", "year"),
    start_tweets = "2019-12-31T10:00:00Z",
    end_tweets = "2020-01-01T10:00:00Z"
  )

```

The same is true for hashtags. We would just need to add the # before the string. So if we were interested in tweets containing #BLM or #BlackLivesMatter then we would do the following:

```{r, eval = F}

tweets <-
  get_all_tweets(
    query = c("#BLM", "#BlackLivesMatter"),
    start_tweets = "2019-12-31T10:00:00Z",
    end_tweets = "2020-01-01T10:00:00Z"
  )

```

This is the equivalent of the OR operator logics detailed by Twitter [here](https://developer.twitter.com/en/docs/twitter-api/tweets/search/integrate/build-a-query).

Note that the "AND" operator is implicit when specifying more than one character string in the query. Thus, when searching for all elements of a character string, a call may look like:

```{r, eval=F}

tweets <- get_all_tweets(
  query = c("twitter social"),
  users = c("cbarrie", "jack"),
  start_tweets = "2020-01-01T00:00:00Z",
  end_tweets = "2020-05-01T00:00:00Z",
  n = 1000
)

```

, which will capture tweets containing *both* the words "twitter" and "social." The same logics apply for hashtag queries.

Finally, we can search for *exact* phrases by using an additional optional parameter `exact_phrase`. So, if we wanted to search tweets containing, for example, the exact phrase "Black Lives Matter," we could do the following:

```{r, eval = F}

tweets <-
  get_all_tweets(
    query = "Black Lives Matter",
    exact_phrase = T,
    start_tweets = "2021-01-04T00:00:00Z",
    end_tweets = "2021-01-04T00:45:00Z",
    n = Inf
  )

```

## Checking your query

When building your query you can also check the query you are building separately with the `build_query()`, which is called inside the `get_all_tweets()` function. 

So, the above query focusing on Seattle is querying the following, which we can build separately as so:

```{r}

build_query(
  query = "happy",
  country = "US",
  place = "seattle",
  lang = "en",
  point_radius = c(-122.33795253639994, 47.60900846404393, 25)
)

```

Here we can see and check the format of the query sent to the Twitter API. 

## Getting user tweets

Finally, we can combine all of the above search functionality when searching by a set of users. To do this we simply need to specify the user or set of users as follows:

```{r, eval = F}

tweets <-
  get_all_tweets(
    users = c("cbarrie", "jack"),
    start_tweets = "2021-01-01T00:00:00Z",
    end_tweets = "2021-06-01T00:00:00Z",
    n = 1000
  )

```

Notice that here we are not specifying any string query. Instead, we are getting all the tweets by users @cbarrie and @jack over a given time period, collecting up to a limit of 1000 tweets.

However, we could also search user tweets and filter by what they're saying at the same time. Here, we would be getting tweets by users @cbarrie and @jack when they were talking about a shared love, "twitter":

```{r, eval = F}

get_all_tweets(
  query = "twitter",
  users = c("cbarrie", "jack"),
  start_tweets = "2021-01-01T00:00:00Z",
  end_tweets = "2021-06-01T00:00:00Z",
  n = 1000
)

```---
title: "Authorization for Twitter Academic Research Product Track"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Authorization for Twitter Academic Research Product Track}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Introduction

In order to use the Twitter Academic Research Product Track you will first need to obtain an authorization token. You will find additional details about the process of obtaining authorization on the Twitter Developer platform.

**In order to gain authorization you first need a Twitter account.**

First, Twitter will ask for details about your academic profile. Per the documentation linked above, they will ask for the following:

> Your full name as it is appears on your institution’s documentation
> 
>   Links to webpages that help establish your identity; provide one or more of the following:
> 
>   - A link to your profile in your institution’s faculty or student directory
>   - A link to your Google Scholar profile
>   - A link to your research group, lab or departmental website where you are listed
> 
>   Information about your academic institution: its name, country, state, and city
> 
>   Your department, school, or lab name
> 
>   Your academic field of study or discipline at this institution
> 
>   Your current role as an academic (whether you are a graduate student, doctoral candidate,       post-doc, professor, research scientist, or other faculty member)

Twitter will then ask for details of the proposed research project. Here, questions include:

> 1. What is the name of your research project?
>
> 2. Does this project receive funding from outside your academic institution? If yes, please list all your sources of funding.
>
> 3. In English, describe your research project. Minimum 200 characters.
>
> 4. In English, describe how Twitter data via the Twitter API will be used in your research project. Minimum 200 characters.
>
> 5. In English, describe your methodology for analyzing Twitter data, Tweets, and/or Twitter users. Minimum 200 characters.
>
> 6. Will your research present Twitter data individually or in aggregate?
>
> 7. In English, describe how you will share the outcomes of your research (include tools, data, and/or other resources you hope to build and share). Minimum 200 characters.
>
> 8. Will your analysis make Twitter content or derived information available to a government entity?

Once you have gained authorization for your project you will be able to see the new project on your Twitter developer portal. First click on the developer portal as below. 


<center>
![](files/twitterdev2.png){width=80%}
</center>


Here you will see your new project, and the name you gave it, appear on the left hand side. Once you have associated an App with this project, it will also appear below the name of the project. Here, I have several Apps authorized to query the basic API. I have one App, named "gencap", that is associated with my Academic Research Product Track project. 

<center>
![](files/twitterdev3.png){width=80%}
</center>

When you click on the project, you will first see how much of your monthly cap of 10m tweets you have spent. You will also see the App associated with your project below the monthly tweet cap usage information.

<center>
![](files/twitterdev4.png){width=80%}
</center>

By clicking on the Settings icons for the App, you will be taken through to the information about the App associated with the project. Here, you will see two options listed, for "Settings" and "Keys and Tokens."

<center>
![](files/twitterdev5.png){width=80%}
</center>

Beside the panel for Bearer Token, you will see an option to Regenerate the token. You can do this if you have not stored the information about the token and no longer have access to it. It is important to store information on the Bearer Token to avoid having to continually regenerate the Bearer Token information.

<center>
![](files/twitterdev6.png){width=80%}
</center>

Once you have the Bearer Token, you are ready to use `academictwitteR`.

It is not advisable to save your bearer token within your scripts. Instead, we advise saving it as an environment variable. Please consult `?get_bearer` for more details.

---
title: "Building a tidy data frame"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Building a tidy data frame}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Introduction

In v0.2 of the package, we include functionality to convert JSON files to various data frame formats. In order to use these features, we recommend the following workflow.

First, you should build your query using the `build_query` function.

```{r}
require(academictwitteR)
require(tibble)
my_query <- build_query(c("#ichbinhanna", "#ichwarhanna"), place = "Berlin")
my_query
```

Then, use the `get_all_tweets` to collect data. Make sure to specify `data_path` and set `bind_tweets` to FALSE.

```{r, eval = F}
get_all_tweets(
  query = my_query,
  start_tweets = "2021-06-01T00:00:00Z",
  end_tweets = "2021-06-20T00:00:00Z",
  n = Inf,
  data_path = "tweetdata",
  bind_tweets = FALSE
)
```

The first format is the so-called "vanilla" format. This vanilla format is the direct output from `jsonlite::read_json`. It can display columns such as `text` just fine. But some columns such as `retweet_count` are nested in list-columns. 

In order to extract user information, it is additionally necessary to set `user = TRUE`. Please also note that the data frame returned in this format is not a tibble. As such, we first need to convert it to a tibble.

```{r, eval = FALSE}
bind_tweets(data_path = "tweetdata") %>% as_tibble
```

```{r, echo = FALSE}
bind_tweets(system.file("extdata", "tweetdata", package = "academictwitteR")) %>% as_tibble
```

The second format is the "raw" format. It is a list of data frames containing all of the data extracted in the API call. Please note that not all data frames are in Boyce-Codd 3rd Normal form, i.e. some columns are still list-column.

```{r, eval = FALSE}
bind_tweets(data_path = "tweetdata", output_format = "raw") %>% names
```

```{r, echo = FALSE}
bind_tweets(system.file("extdata", "tweetdata", package = "academictwitteR"), output_format = "raw") %>% names
```

The third format is the "tidy" format. It is an "opinionated" format, which we believe to contain all essential columns for social media research. By default, it is a tibble.

```{r, eval = FALSE}
bind_tweets(data_path = "tweetdata", output_format = "tidy")
```

```{r, echo = FALSE}
bind_tweets(system.file("extdata", "tweetdata", package = "academictwitteR"), output_format = "tidy")
```

It has the following features / caveats:

1. It has both the data about tweets, their authors, and "source tweets", a.k.a. referenced tweets. Columns are named according to these three sources. The primary keys of these three sources are named `tweet_id`, `author_id` and `sourcetweet_id` respectively.
2. By default, the `text` field of a retweet is truncated. However, the full-text original tweet is located in `sourcetweet_text`.
3. The replied tweets of a reply is not counted as `sourcetweet_text`. If you need that data, please follow the clue using the `conversation_id`.
4. Many data extracted from `text` by Twitter are not available in the tidy format, e.g. list of hashtags, cashtags, urls, entities, context annotations etc. If you need those columns, please consider using the "raw" format above.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compliance.R
\name{get_compliance_result}
\alias{get_compliance_result}
\title{Get Compliance Result}
\usage{
get_compliance_result(id, bearer_token = get_bearer(), verbose = TRUE)
}
\arguments{
\item{id}{string, the job id}

\item{bearer_token}{string, bearer token}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}
}
\value{
a data frame
}
\description{
This function retrieves the information for a single compliance job.
}
\examples{
\dontrun{
get_compliance_result("1460077048991555585")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_functions.R
\name{get_geo_tweets}
\alias{get_geo_tweets}
\title{Get tweets for query containing geo information
`r lifecycle::badge("deprecated")}
\usage{
get_geo_tweets(
  query,
  start_tweets,
  end_tweets,
  bearer_token = get_bearer(),
  n = 100,
  file = NULL,
  data_path = NULL,
  bind_tweets = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{query}{string or character vector, search query or queries}

\item{start_tweets}{string, starting date}

\item{end_tweets}{string, ending date}

\item{bearer_token}{string, bearer token}

\item{n}{integer, upper limit of tweets to be fetched}

\item{file}{string, name of the resulting RDS file}

\item{data_path}{string, if supplied, fetched data can be saved to the designated path as jsons}

\item{bind_tweets}{If \code{TRUE}, tweets captured are bound into a data.frame for assignment}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{...}{arguments will be passed to \code{\link[=build_query]{build_query()}} function. See \code{?build_query()} for further information.}
}
\value{
a data.frame
}
\description{
This function collects tweets containing strings or
hashtags between specified date ranges that also contain Tweet-specific geolocation data provided by the
Twitter user. This can be either a location in the form of a Twitter place, with the corresponding display
name, geo polygon, and other fields, or in rare cases, a geo lat-long coordinate. Note: Operators matching
on place (Tweet geo) will only include matches from original tweets. Tweet-level data is stored in a data/
path as a series of JSONs beginning "data_"; User-level data is stored as a series of JSONs beginning "users_".
If a filename is supplied, this function will save the result as a RDS file, otherwise it will return the results as a dataframe.
}
\examples{
\dontrun{
bearer_token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
get_geo_tweets("protest", "2020-01-01T00:00:00Z", "2020-01-05T00:00:00Z", 
               bearer_token, data_path = "data/")
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_user_profile.R
\name{get_user_profile}
\alias{get_user_profile}
\title{Get user profile}
\usage{
get_user_profile(x, bearer_token = get_bearer())
}
\arguments{
\item{x}{string containing one user id or a vector of user ids}

\item{bearer_token}{string, bearer token}
}
\value{
a data frame
}
\description{
This function fetches user-level information for a vector of user IDs.
}
\examples{
\dontrun{
bearer_token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
users <- c("2244994945", "6253282")
get_user_profile(users, bearer_token)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_functions.R
\name{get_video_tweets}
\alias{get_video_tweets}
\title{Get tweets for query containing videos
`r lifecycle::badge("deprecated")}
\usage{
get_video_tweets(
  query,
  start_tweets,
  end_tweets,
  bearer_token = get_bearer(),
  n = 100,
  file = NULL,
  data_path = NULL,
  bind_tweets = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{query}{string or character vector, search query or queries}

\item{start_tweets}{string, starting date}

\item{end_tweets}{string, ending date}

\item{bearer_token}{string, bearer token}

\item{n}{integer, upper limit of tweets to be fetched}

\item{file}{string, name of the resulting RDS file}

\item{data_path}{string, if supplied, fetched data can be saved to the designated path as jsons}

\item{bind_tweets}{If \code{TRUE}, tweets captured are bound into a data.frame for assignment}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{...}{arguments will be passed to \code{\link[=build_query]{build_query()}} function. See \code{?build_query()} for further information.}
}
\value{
a data frame
}
\description{
This function collects tweets containing strings or hashtags between specified date ranges
that also contain native Twitter videos, uploaded directly to Twitter. This will not match
on videos created with Periscope, or Tweets with links to other video hosting sites. Tweet-level
data is stored in a data/ path as a series of JSONs beginning "data_"; User-level data is stored
as a series of JSONs beginning "users_". If a filename is supplied, this function will save the
result as a RDS file, otherwise it will return the results as a dataframe.
}
\examples{
\dontrun{
bearer_token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
get_video_tweets("#BLM", "2020-01-01T00:00:00Z", "2020-01-05T00:00:00Z",
                  bearer_token, data_path = "data/")
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_all_tweets.R
\name{count_all_tweets}
\alias{count_all_tweets}
\title{Count tweets from full archive search}
\usage{
count_all_tweets(
  query = NULL,
  start_tweets,
  end_tweets,
  bearer_token = get_bearer(),
  n = 100,
  file = NULL,
  data_path = NULL,
  export_query = TRUE,
  bind_tweets = TRUE,
  granularity = "day",
  verbose = TRUE,
  ...
)
}
\arguments{
\item{query}{string or character vector, search query or queries}

\item{start_tweets}{string, starting date}

\item{end_tweets}{string, ending date}

\item{bearer_token}{string, bearer token}

\item{n}{integer, upper limit of tweet counts to be fetched (i.e., for 365 days n must be at least 365). Default is 100.}

\item{file}{string, name of the resulting RDS file}

\item{data_path}{string, if supplied, fetched data can be saved to the designated path as jsons}

\item{export_query}{If \code{TRUE}, queries are exported to data_path}

\item{bind_tweets}{If \code{TRUE}, tweets captured are bound into a data.frame for assignment}

\item{granularity}{string, the granularity for the search counts results. Options are "day"; "hour"; "minute". Default is day.}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{...}{arguments will be passed to \code{build_query()} function. See \code{?build_query()} for further information.}
}
\value{
a data.frame
}
\description{
This function returns aggregate counts of tweets by query string or strings
between specified date ranges.
}
\examples{
\dontrun{

count_all_tweets(query = "Hogmanay", 
          start_tweets = "2019-12-2700:00:00Z", 
          end_tweets = "2020-01-05T00:00:00Z", 
          bearer_token = get_bearer())
          
count_all_tweets(query = "Hogmanay", 
          start_tweets = "2019-12-27T00:00:00Z", 
          end_tweets = "2020-01-05T00:00:00Z", 
          bearer_token = get_bearer(),
          granularity = "hour",
          n = 500)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_collection.R
\name{update_collection}
\alias{update_collection}
\title{Update previous collection session}
\usage{
update_collection(
  data_path,
  end_tweets,
  bearer_token = get_bearer(),
  verbose = TRUE,
  ...
)
}
\arguments{
\item{data_path}{string, name of an existing data_path}

\item{end_tweets}{string, ending date}

\item{bearer_token}{string, bearer token}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{...}{arguments will be passed to \code{get_all_tweets()} function. See \code{?get_all_tweets()} for further information.}
}
\value{
a data.frame
}
\description{
This function continues a previous collection session with a new end date.
For this function to work, export_query must be set to "TRUE" during the original collection.
}
\examples{
\dontrun{
update_collection(data_path = "data", "2020-01-03T00:00:00Z", bearer_token = get_bearer())
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_functions.R
\name{get_media_tweets}
\alias{get_media_tweets}
\title{Get tweets for query containing media
`r lifecycle::badge("deprecated")}
\usage{
get_media_tweets(
  query,
  start_tweets,
  end_tweets,
  bearer_token = get_bearer(),
  n = 100,
  file = NULL,
  data_path = NULL,
  bind_tweets = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{query}{string or character vector, search query or queries}

\item{start_tweets}{string, starting date}

\item{end_tweets}{string, ending date}

\item{bearer_token}{string, bearer token}

\item{n}{integer, upper limit of tweets to be fetched}

\item{file}{string, name of the resulting RDS file}

\item{data_path}{string, if supplied, fetched data can be saved to the designated path as jsons}

\item{bind_tweets}{If \code{TRUE}, tweets captured are bound into a data.frame for assignment}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{...}{arguments will be passed to \code{\link[=build_query]{build_query()}} function. See \code{?build_query()} for further information.}
}
\value{
a data.frame
}
\description{
This function collects tweets containing the strings or hashtags between specified date ranges
that also contain a media object, such as a photo, GIF, or video, as determined by Twitter. Tweet-level
data is stored in a data/ path as a series of JSONs beginning "data_"; User-level data is stored as a
series of JSONs beginning "users_". If a filename is supplied, this function will save the result
as a RDS file, otherwise it will return the results as a dataframe.
}
\examples{
\dontrun{
bearer_token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
get_media_tweets("#BLM", "2020-01-01T00:00:00Z", "2020-01-05T00:00:00Z",
                 bearer_token, data_path = "data/")
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_retweeted_by.R
\name{get_retweeted_by}
\alias{get_retweeted_by}
\title{Get users who has retweeted a tweet}
\usage{
get_retweeted_by(
  x,
  bearer_token = get_bearer(),
  data_path = NULL,
  verbose = TRUE
)
}
\arguments{
\item{x}{string containing one tweet id or a vector of tweet ids}

\item{bearer_token}{string, bearer token}

\item{data_path}{string, if supplied, fetched data can be saved to the designated path as jsons}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}
}
\value{
a data frame
}
\description{
This function fetches users who retweeted a tweet
}
\examples{
\dontrun{
tweets <- c("1392887366507970561","1409931481552543749")
get_retweeted_by(tweets, bearer_token = get_bearer())
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_functions.R
\name{bind_user_jsons}
\alias{bind_user_jsons}
\title{Bind user information stored as JSON files
`r lifecycle::badge("deprecated")}
\usage{
bind_user_jsons(data_path, verbose = TRUE)
}
\arguments{
\item{data_path}{string, file path to directory of stored tweets data saved as users_\emph{id}.json}

\item{verbose}{If \code{FALSE}, messages are suppressed}
}
\value{
a data.frame
}
\description{
Bind user information stored as JSON files
`r lifecycle::badge("deprecated")
}
\examples{
\dontrun{
bind_user_jsons("data_path = "data/"")
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_functions.R
\name{get_country_tweets}
\alias{get_country_tweets}
\title{Get tweets with country parameter
`r lifecycle::badge("deprecated")}
\usage{
get_country_tweets(
  query,
  country,
  start_tweets,
  end_tweets,
  bearer_token = get_bearer(),
  n = 100,
  file = NULL,
  data_path = NULL,
  bind_tweets = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{query}{string or character vector, search query or queries}

\item{country, }{string, name of country as ISO alpha-2 code e.g. "GB"}

\item{start_tweets}{string, starting date}

\item{end_tweets}{string, ending date}

\item{bearer_token}{string, bearer token}

\item{n}{integer, upper limit of tweets to be fetched}

\item{file}{string, name of the resulting RDS file}

\item{data_path}{string, if supplied, fetched data can be saved to the designated path as jsons}

\item{bind_tweets}{If \code{TRUE}, tweets captured are bound into a data.frame for assignment}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{...}{arguments will be passed to \code{\link[=build_query]{build_query()}} function. See \code{?build_query()} for further information.}
}
\value{
a data.frame
}
\description{
This function collects tweets containing strings or hashtags
between specified date ranges filtering by country. Tweet-level data is stored in a data/
path as a series of JSONs beginning "data_"; User-level data is stored as a series of
JSONs beginning "users_". If a filename is supplied, this function will
save the result as a RDS file, otherwise it will return the results as a data.frame.
}
\examples{
\dontrun{
bearer_token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
get_country_tweets("happy", country = "GB",
                   "2021-01-01T00:00:00Z", "2021-01-01T00:10:00Z",
                   bearer_token, data_path = "data/")
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_functions.R
\name{get_image_tweets}
\alias{get_image_tweets}
\title{Get tweets containing images
`r lifecycle::badge("deprecated")}
\usage{
get_image_tweets(
  query,
  start_tweets,
  end_tweets,
  bearer_token = get_bearer(),
  n = 100,
  file = NULL,
  data_path = NULL,
  bind_tweets = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{query}{string or character vector, search query or queries}

\item{start_tweets}{string, starting date}

\item{end_tweets}{string, ending date}

\item{bearer_token}{string, bearer token}

\item{n}{integer, upper limit of tweets to be fetched}

\item{file}{string, name of the resulting RDS file}

\item{data_path}{string, if supplied, fetched data can be saved to the designated path as jsons}

\item{bind_tweets}{If \code{TRUE}, tweets captured are bound into a data.frame for assignment}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{...}{arguments will be passed to \code{\link[=build_query]{build_query()}} function. See \code{?build_query()} for further information.}
}
\value{
a data.frame
}
\description{
This function collects tweets containing strings or
hashtags between specified date ranges that also contain (a recognized URL to) an image. Tweet-level data
is stored in a data/ path as a series of JSONs beginning "data_"; User-level data is stored as a series
of JSONs beginning "users_". If a filename is supplied, this function will save the result as a RDS file,
otherwise it will return the results as a dataframe.
}
\examples{
\dontrun{
bearer_token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
get_image_tweets("#BLM", "2020-01-01T00:00:00Z", "2020-01-05T00:00:00Z",
                 bearer_token, data_path = "data/")
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_functions.R
\name{get_to_tweets}
\alias{get_to_tweets}
\title{Get tweets to users
`r lifecycle::badge("deprecated")}
\usage{
get_to_tweets(
  users,
  start_tweets,
  end_tweets,
  bearer_token = get_bearer(),
  n = 100,
  file = NULL,
  data_path = NULL,
  bind_tweets = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{users}{character vector, user handles from which to collect data}

\item{start_tweets}{string, starting date}

\item{end_tweets}{string, ending date}

\item{bearer_token}{string, bearer token}

\item{n}{integer, upper limit of tweets to be fetched}

\item{file}{string, name of the resulting RDS file}

\item{data_path}{string, if supplied, fetched data can be saved to the designated path as jsons}

\item{bind_tweets}{If \code{TRUE}, tweets captured are bound into a data.frame for assignment}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{...}{arguments will be passed to \code{\link[=build_query]{build_query()}} function. See \code{?build_query()} for further information.}
}
\value{
a data frame
}
\description{
This function collects tweets between specified date ranges that are
in reply to the specified user(s). Tweet-level data is stored in a data/ path as a series of JSONs beginning
"data_"; User-level data is stored as a series of JSONs beginning "users_". If a filename is supplied,
this function will save the result as a RDS file, otherwise it will return the results as a dataframe.
}
\examples{
\dontrun{
bearer_token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
users <- c("uoessps", "spsgradschool")
get_to_tweets(users, "2020-01-01T00:00:00Z", "2020-01-05T00:00:00Z",
             bearer_token, data_path = "data/")
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_functions.R
\name{get_retweets_of_user}
\alias{get_retweets_of_user}
\title{Get retweets of user
`r lifecycle::badge("deprecated")}
\usage{
get_retweets_of_user(
  users,
  start_tweets,
  end_tweets,
  bearer_token = get_bearer(),
  n = 100,
  file = NULL,
  data_path = NULL,
  bind_tweets = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{users}{character vector, user handles from which to collect data}

\item{start_tweets}{string, starting date}

\item{end_tweets}{string, ending date}

\item{bearer_token}{string, bearer token}

\item{n}{integer, upper limit of tweets to be fetched}

\item{file}{string, name of the resulting RDS file}

\item{data_path}{string, if supplied, fetched data can be saved to the designated path as jsons}

\item{bind_tweets}{If \code{TRUE}, tweets captured are bound into a data.frame for assignment}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{...}{arguments will be passed to \code{\link[=build_query]{build_query()}} function. See \code{?build_query()} for further information.}
}
\value{
a data frame
}
\description{
This function collects retweets of tweets by a user or set of users between specified date ranges.
Tweet-level data is stored in a data/ path as a series of JSONs beginning "data_"; User-level data is stored as
a series of JSONs beginning "users_". If a filename is supplied, this function will save the result as a RDS file,
otherwise it will return the results as a dataframe.
}
\examples{
\dontrun{
bearer_token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
users <- c("cbarrie", "justin_ct_ho")
get_retweets_of_user(users, "2020-01-01T00:00:00Z", "2020-04-05T00:00:00Z",
                     bearer_token, data_path = "data/")
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_user_edges.R
\name{get_user_followers}
\alias{get_user_followers}
\title{Get user followers}
\usage{
get_user_followers(x, bearer_token = get_bearer(), ...)
}
\arguments{
\item{x}{string containing one user id or a vector of user ids}

\item{bearer_token}{string, bearer token}

\item{...}{arguments passed to other backend functions}
}
\value{
a data frame
}
\description{
This function fetches users who are followers of the specified user ID.
}
\examples{
\dontrun{
bearer_token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
users <- "2244994945"
get_user_followers(users, bearer_token = get_bearer())
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_user_timeline.R
\name{get_user_timeline}
\alias{get_user_timeline}
\title{Get tweets by a single user}
\usage{
get_user_timeline(
  x,
  start_tweets,
  end_tweets,
  bearer_token = get_bearer(),
  n = 100,
  file = NULL,
  data_path = NULL,
  export_query = TRUE,
  bind_tweets = TRUE,
  page_n = 100,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{x}{string containing one user id or a vector of user ids}

\item{start_tweets}{string, starting date}

\item{end_tweets}{string, ending date}

\item{bearer_token}{string, bearer token}

\item{n}{integer, upper limit of tweets to be fetched}

\item{file}{string, name of the resulting RDS file}

\item{data_path}{string, if supplied, fetched data can be saved to the designated path as jsons}

\item{export_query}{If \code{TRUE}, queries are exported to data_path}

\item{bind_tweets}{If \code{TRUE}, tweets captured are bound into a data.frame for assignment}

\item{page_n}{integer, amount of tweets to be returned by per page}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{...}{arguments will be passed to \code{build_query()} function. See \code{?build_query()} for further information.}
}
\value{
a data.frame
}
\description{
This function collects tweets by an user ID from the users endpoint.
}
\details{
Only the most recent 3,200 Tweets can be retrieved.

If a filename is supplied, the function will
save the result as a RDS file.

If a data path is supplied, the function will also return
tweet-level data in a data/ path as a series of JSONs beginning "data_";
while user-level data will be returned as a series of JSONs beginning "users_".

When bind_tweets is \code{TRUE}, the function returns a data frame.
}
\examples{
\dontrun{

get_user_timeline("2244994945",
                  start_tweets = "2020-01-01T00:00:00Z", 
                  end_tweets = "2021-05-14T00:00:00Z",
                  bearer_token = get_bearer(),
                  n = 200)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_functions.R
\name{get_mentions_tweets}
\alias{get_mentions_tweets}
\title{Get tweets for query containing mentions of another user
`r lifecycle::badge("deprecated")}
\usage{
get_mentions_tweets(
  query,
  start_tweets,
  end_tweets,
  bearer_token = get_bearer(),
  n = 100,
  file = NULL,
  data_path = NULL,
  bind_tweets = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{query}{string or character vector, search query or queries}

\item{start_tweets}{string, starting date}

\item{end_tweets}{string, ending date}

\item{bearer_token}{string, bearer token}

\item{n}{integer, upper limit of tweets to be fetched}

\item{file}{string, name of the resulting RDS file}

\item{data_path}{string, if supplied, fetched data can be saved to the designated path as jsons}

\item{bind_tweets}{If \code{TRUE}, tweets captured are bound into a data.frame for assignment}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{...}{arguments will be passed to \code{\link[=build_query]{build_query()}} function. See \code{?build_query()} for further information.}
}
\value{
a data frame
}
\description{
This function collects tweets containing strings or
hashtags between specified date ranges that also contain mentions of another Twitter user. Tweet-level data
is stored in a data/ path as a series of JSONs beginning "data_"; User-level data is stored as a series of
JSONs beginning "users_". If a filename is supplied, this function will save the result as a RDS file, otherwise,
it will return the results as a dataframe.
}
\examples{
\dontrun{
bearer_token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
get_mentions_tweets("#nowplaying", "2020-01-01T00:00:00Z", "2020-01-05T00:00:00Z",
                    bearer_token, data_path = "data/")
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compliance.R
\name{create_compliance_job}
\alias{create_compliance_job}
\title{Create Compliance Job}
\usage{
create_compliance_job(
  x,
  type = "tweets",
  bearer_token = get_bearer(),
  force_ids = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{x}{either a character vector of Tweet IDs or user IDs; or a plain text file that each line contains a Tweet ID or user ID.}

\item{type}{the type of the job, whether "tweets" or "users".}

\item{bearer_token}{string, bearer token}

\item{force_ids}{logical, make sure \code{x} is treated as a character vector of Tweet IDs or user IDs.}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}
}
\value{
the job ID (invisibly)
}
\description{
This function creates a new compliance job and upload the Tweet IDs or user IDs. By default, the parameter \code{x} with the length of one is assumed to be a text file containing either Tweet IDs or iser IDs. This default behavior can be bypassed using \code{force_ids} For example, if you want to check for just a single Tweet ID.
}
\examples{
\dontrun{
create_compliance_job(x = "tweetids.txt", type = "tweets")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bind_tweets.R
\name{bind_tweets}
\alias{bind_tweets}
\alias{convert_json}
\title{Bind information stored as JSON files}
\usage{
bind_tweets(data_path, user = FALSE, verbose = TRUE, output_format = NA)

convert_json(data_file, output_format = "tidy")
}
\arguments{
\item{data_path}{string, file path to directory of stored tweets data saved as data_\emph{id}.json and users_\emph{id}.json}

\item{user}{If \code{FALSE}, this function binds JSON files into a data frame containing tweets; data frame containing user information otherwise. Ignore if \code{output_format} is not NA}

\item{verbose}{If \code{FALSE}, messages are suppressed}

\item{output_format}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}} string, if it is not NA, this function return an unprocessed data.frame containing either tweets or user information. Currently, this function supports the following format(s)
\itemize{
\item{"raw"}{List of data frames; Note: not all data frames are in Boyce-Codd 3rd Normal Form}
\item{"tidy"}{Tidy format; all essential columns are available}
}}

\item{data_file}{string, a single file path to a JSON file; or a vector of file paths to JSON files of stored tweets data saved as data_\emph{id}.json}
}
\value{
a data.frame containing either tweets or user information
}
\description{
This function binds information stored as JSON files. The experimental function \code{convert_json} converts individual JSON files into either "raw" or "tidy" format.
}
\details{
By default, \code{bind_tweets} binds into a data frame containing tweets (from data_\emph{id}.json files).

If users is TRUE, it binds into a data frame containing user information (from users_\emph{id}.json).
}
\examples{
\dontrun{
# bind json files in the directory "data" into a data frame containing tweets
bind_tweets(data_path = "data/")

# bind json files in the directory "data" into a data frame containing user information
bind_tweets(data_path = "data/", user = TRUE)

# bind json files in the directory "data" into a "tidy" data frame / tibble
bind_tweets(data_path = "data/", user = TRUE, output_format = "tidy")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compliance.R
\name{list_compliance_jobs}
\alias{list_compliance_jobs}
\title{List Compliance Jobs}
\usage{
list_compliance_jobs(type = "tweets", bearer_token = get_bearer())
}
\arguments{
\item{type}{the type of the job, whether "tweets" or "users".}

\item{bearer_token}{string, bearer token}
}
\value{
a data frame
}
\description{
This function lists all compliance jobs.
}
\examples{
\dontrun{
list_compliance_jobs()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_liking_users.R
\name{get_liking_users}
\alias{get_liking_users}
\title{Get liking users}
\usage{
get_liking_users(x, bearer_token = get_bearer(), verbose = TRUE)
}
\arguments{
\item{x}{string containing one tweet id or a vector of tweet ids}

\item{bearer_token}{string, bearer token}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}
}
\value{
a data frame
}
\description{
This function fetches users who liked a tweet or tweets.
}
\examples{
\dontrun{
tweet <- "1387744422729748486"
get_liking_users(tweet, bearer_token = get_bearer())
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_all_tweets.R
\name{get_all_tweets}
\alias{get_all_tweets}
\title{Get tweets from full archive search}
\usage{
get_all_tweets(
  query = NULL,
  start_tweets,
  end_tweets,
  bearer_token = get_bearer(),
  n = 100,
  file = NULL,
  data_path = NULL,
  export_query = TRUE,
  bind_tweets = TRUE,
  page_n = 500,
  context_annotations = FALSE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{query}{string or character vector, search query or queries}

\item{start_tweets}{string, starting date}

\item{end_tweets}{string, ending date}

\item{bearer_token}{string, bearer token}

\item{n}{integer, upper limit of tweets to be fetched}

\item{file}{string, name of the resulting RDS file}

\item{data_path}{string, if supplied, fetched data can be saved to the designated path as jsons}

\item{export_query}{If \code{TRUE}, queries are exported to data_path}

\item{bind_tweets}{If \code{TRUE}, tweets captured are bound into a data.frame for assignment}

\item{page_n}{integer, amount of tweets to be returned by per page}

\item{context_annotations}{If \code{TRUE}, context_annotations will be fetched. Note it will limit the page_n to 100 due restrictions of Twitter API.}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{...}{arguments will be passed to \code{\link[=build_query]{build_query()}} function. See \code{?build_query()} for further information.}
}
\value{
When bind_tweets is \code{TRUE} (default), the function returns a data frame. Nothing otherwise.
}
\description{
This function collects tweets by query string or strings
between specified date ranges.
}
\details{
The function can also collect tweets by users. These may be specified alongside
a query string or without. When no query string is supplied, the function collects
all tweets by that user.

If a filename is supplied, the function will
save the result as a RDS file.

If a data path is supplied, the function will also return
tweet-level data in a data/ path as a series of JSONs beginning "data_";
while user-level data will be returned as a series of JSONs beginning "users_".
}
\examples{
\dontrun{
bearer_token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"

get_all_tweets(query = "BLM", 
               start_tweets = "2020-01-01T00:00:00Z", 
               end_tweets = "2020-01-05T00:00:00Z", 
               bearer_token = get_bearer(), 
               data_path = "data",
               n = 500)
  
get_all_tweets(users = c("cbarrie", "jack"),
               start_tweets = "2021-01-01T00:00:00Z", 
               end_tweets = "2021-06-01T00:00:00Z",
               bearer_token = get_bearer(), 
               n = 1000)
                            
get_all_tweets(start_tweets = "2021-01-01T00:00:00Z", 
               end_tweets = "2021-06-01T00:00:00Z",
               bearer_token = get_bearer(), 
               n = 1500, 
               conversation_id = "1392887366507970561")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_user_edges.R
\name{get_user_following}
\alias{get_user_following}
\title{Get user following}
\usage{
get_user_following(x, bearer_token = get_bearer(), ...)
}
\arguments{
\item{x}{string containing one user id or a vector of user ids}

\item{bearer_token}{string, bearer token}

\item{...}{arguments passed to other backend functions}
}
\value{
a data frame
}
\description{
This function fetches a list of users the specified user ID is following.
}
\examples{
\dontrun{
bearer_token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
users <- "2244994945"
get_user_following(users, bearer_token)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_functions.R
\name{get_lang_tweets}
\alias{get_lang_tweets}
\title{Get tweets in particular language
`r lifecycle::badge("deprecated")}
\usage{
get_lang_tweets(
  query,
  lang,
  start_tweets,
  end_tweets,
  bearer_token = get_bearer(),
  n = 100,
  file = NULL,
  data_path = NULL,
  bind_tweets = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{query}{string or character vector, search query or queries}

\item{lang, }{string, a single BCP 47 language identifier e.g. "fr"}

\item{start_tweets}{string, starting date}

\item{end_tweets}{string, ending date}

\item{bearer_token}{string, bearer token}

\item{n}{integer, upper limit of tweets to be fetched}

\item{file}{string, name of the resulting RDS file}

\item{data_path}{string, if supplied, fetched data can be saved to the designated path as jsons}

\item{bind_tweets}{If \code{TRUE}, tweets captured are bound into a data.frame for assignment}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{...}{arguments will be passed to \code{\link[=build_query]{build_query()}} function. See \code{?build_query()} for further information.}
}
\value{
a data.frame
}
\description{
This function collects tweets containing strings or hashtags
between specified date ranges filtering by language. Tweet-level data is stored in a data/
path as a series of JSONs beginning "data_"; User-level data is stored as a series of
JSONs beginning "users_". If a filename is supplied, this function will
save the result as a RDS file, otherwise it will return the results as a dataframe.
}
\examples{
\dontrun{
bearer_token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
get_lang_tweets("bonne", lang= "fr",
                "2021-01-01T00:00:00Z", "2021-01-01T00:10:00Z",
                bearer_token, data_path = "data/")
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\arguments{
\item{lhs}{A value or the magrittr placeholder.}

\item{rhs}{A function call using the magrittr semantics.}
}
\value{
The result of calling \code{rhs(lhs)}.
}
\description{
See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_functions.R
\name{get_place_tweets}
\alias{get_place_tweets}
\title{Get tweets with place parameter
`r lifecycle::badge("deprecated")}
\usage{
get_place_tweets(
  query,
  place,
  start_tweets,
  end_tweets,
  bearer_token = get_bearer(),
  n = 100,
  file = NULL,
  data_path = NULL,
  bind_tweets = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{query}{string or character vector, search query or queries}

\item{place, }{string, name of place e.g. "new york city"}

\item{start_tweets}{string, starting date}

\item{end_tweets}{string, ending date}

\item{bearer_token}{string, bearer token}

\item{n}{integer, upper limit of tweets to be fetched}

\item{file}{string, name of the resulting RDS file}

\item{data_path}{string, if supplied, fetched data can be saved to the designated path as jsons}

\item{bind_tweets}{If \code{TRUE}, tweets captured are bound into a data.frame for assignment}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{...}{arguments will be passed to \code{\link[=build_query]{build_query()}} function. See \code{?build_query()} for further information.}
}
\value{
a data.frame
}
\description{
This function collects tweets containing strings or hashtags
between specified date ranges filtering by place. Tweet-level data is stored in a data/
path as a series of JSONs beginning "data_"; User-level data is stored as a series of
JSONs beginning "users_". If a filename is supplied, this function will
save the result as a RDS file, otherwise, it will return the results as a data.frame.
}
\examples{
\dontrun{
bearer_token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
get_place_tweets("happy", place = "London",
                 "2021-01-01T00:00:00Z", "2021-01-01T00:10:00Z",
                 bearer_token, data_path = "data/")
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hydrate_tweets.R
\name{hydrate_tweets}
\alias{hydrate_tweets}
\title{Hydrate Tweets Based On Tweet IDs}
\usage{
hydrate_tweets(
  ids,
  bearer_token = get_bearer(),
  data_path = NULL,
  context_annotations = FALSE,
  bind_tweets = TRUE,
  verbose = TRUE,
  errors = FALSE
)
}
\arguments{
\item{ids}{a character vector of Tweet IDs}

\item{bearer_token}{string, bearer token}

\item{data_path}{string, if supplied, fetched data can be saved to the designated path as jsons}

\item{context_annotations}{If \code{TRUE}, context_annotations will be fetched.}

\item{bind_tweets}{If \code{TRUE}, tweets captured are bound into a data.frame for assignment}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{errors}{logical, if \code{TRUE}, the error capturing mechanism is enabled. See details below.}
}
\value{
When bind_tweets is \code{TRUE}, the function returns a data frame. The \code{data_path} (invisibly) if \code{bind_tweets} is \code{FALSE}
}
\description{
This function is helpful for hydrating Tweet IDs (i.e. getting the full content of tweets from a list of Tweet IDs).
}
\details{
When the error capturing mechanism is enabled, Tweets IDs that cannot be queried (e.g. with error) are stored as \code{errors_*.json} files. If \code{bind_tweets} is TRUE, those error Tweets IDs are retained in the returned data.frame with the column \code{error} indicating the error.
}
\examples{
\dontrun{
hydrate_tweets(c("1266876474440761346", "1266868259925737474", "1266867327079002121",
"1266866660713127936", "1266864490446012418", "1266860737244336129",
"1266859737615826944", "1266859455586676736", "1266858090143588352",
"1266857669157097473"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_bearer.R
\name{set_bearer}
\alias{set_bearer}
\title{Set bearer token}
\usage{
set_bearer()
}
\description{
This function lets the user add their bearer token to the \code{.Renviron} file.
}
\details{
It is in general not safe to 1) hard code your bearer token in your R script or
2) have your bearer token in your command history.

\code{set_bearer} opens the .Renviron file
for the user and provides instructions on how to add the bearer token, which requires the
addition of just one line in the \code{.Renviron} file, following the format TWITTER_BEARER=YOURTOKENHERE.

Replace YOURTOKENHERE with your own token.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_user_id.R
\name{get_user_id}
\alias{get_user_id}
\title{Get user id}
\usage{
get_user_id(
  usernames,
  bearer_token = get_bearer(),
  all = FALSE,
  keep_na = TRUE
)
}
\arguments{
\item{usernames}{character vector containing screen names to be queried}

\item{bearer_token}{string, bearer token}

\item{all}{logical, default FALSE to get a character vector of user IDs. Set it to TRUE to get a data frame, see below}

\item{keep_na}{logical, default TRUE to keep usernames that cannot be queried. Set it to TRUE to exclude those usernames. Only useful when all is FALSE}
}
\value{
a string vector with the id of each of the users unless all = TRUE. If all = TRUE, a data.frame with ids, names (showed on the screen) and usernames is returned.
}
\description{
This function get the user IDs (e.g. 1349149096909668363) of given usernames, e.g. "potus".
}
\examples{
\dontrun{
bearer_token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
users <- c("Twitter", "TwitterDev")
get_user_id(users, bearer_token)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_functions.R
\name{get_url_tweets}
\alias{get_url_tweets}
\title{Get tweets containing URL
`r lifecycle::badge("deprecated")}
\usage{
get_url_tweets(
  query,
  start_tweets,
  end_tweets,
  bearer_token = get_bearer(),
  n = 100,
  file = NULL,
  data_path = NULL,
  bind_tweets = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{query}{string, url}

\item{start_tweets}{string, starting date}

\item{end_tweets}{string, ending date}

\item{bearer_token}{string, bearer token}

\item{n}{integer, upper limit of tweets to be fetched}

\item{file}{string, name of the resulting RDS file}

\item{data_path}{string, if supplied, fetched data can be saved to the designated path as jsons}

\item{bind_tweets}{If \code{TRUE}, tweets captured are bound into a data.frame for assignment}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{...}{arguments will be passed to \code{\link[=build_query]{build_query()}} function. See \code{?build_query()} for further information.}
}
\value{
a data.frame
}
\description{
This function collects tweets containing a given url between specified date ranges.
Tweet-level data is stored in a data/ path as a series of JSONs beginning "data_"; User-level data is stored as a series of
JSONs beginning "users_". If a filename is supplied, this function will save the result as a RDS file, otherwise
it will return the results as a dataframe.
}
\examples{
\dontrun{
bearer_token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
tweets <- get_url_tweets("https://www.theguardian.com/",
"2020-01-01T00:00:00Z", "2020-04-04T00:00:00Z", bearer_token, data_path = "data/")
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_functions.R
\name{get_bbox_tweets}
\alias{get_bbox_tweets}
\title{Get tweets within bounding box
`r lifecycle::badge("deprecated")}
\usage{
get_bbox_tweets(
  query,
  bbox,
  start_tweets,
  end_tweets,
  bearer_token = get_bearer(),
  n = 100,
  file = NULL,
  data_path = NULL,
  bind_tweets = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{query}{string or character vector, search query or queries}

\item{bbox}{numeric, a vector of four bounding box coordinates from west longitude to north latitude}

\item{start_tweets}{string, starting date}

\item{end_tweets}{string, ending date}

\item{bearer_token}{string, bearer token}

\item{n}{integer, upper limit of tweets to be fetched}

\item{file}{string, name of the resulting RDS file}

\item{data_path}{string, if supplied, fetched data can be saved to the designated path as jsons}

\item{bind_tweets}{If \code{TRUE}, tweets captured are bound into a data.frame for assignment}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{...}{arguments will be passed to \code{\link[=build_query]{build_query()}} function. See \code{?build_query()} for further information.}
}
\value{
a data.frame
}
\description{
This function collects tweets containing strings or hashtags
between specified date ranges filtering by bounding box. Tweet-level data is stored in a data/
path as a series of JSONs beginning "data_"; User-level data is stored as a series of
JSONs beginning "users_". If a filename is supplied, this function will
save the result as a RDS file, otherwise it will return the results as a dataframe.
Note: width and height of the bounding box must be less than 25mi.
}
\examples{
\dontrun{
bearer_token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
tweets <- get_bbox_tweets("happy", bbox= c(-0.222473,51.442453,0.072784,51.568534),
                           "2021-01-01T00:00:00Z", "2021-02-01T10:00:00Z",
                           bearer_token, data_path = "data/")
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resume_collection.R
\name{resume_collection}
\alias{resume_collection}
\title{Resume previous collection}
\usage{
resume_collection(data_path, bearer_token = get_bearer(), verbose = TRUE, ...)
}
\arguments{
\item{data_path}{string, name of an existing data_path}

\item{bearer_token}{string, bearer token}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{...}{arguments will be passed to \code{get_all_tweets()} function. See \code{?get_all_tweets()} for further information.}
}
\value{
a data.frame
}
\description{
This function resumes a previous interrupted collection session.
}
\details{
For this function to work, export_query must be set to "TRUE" during the original collection.
}
\examples{
\dontrun{
resume_collection(data_path = "data", bearer_token = get_bearer())
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_bearer.R
\name{get_bearer}
\alias{get_bearer}
\title{Manage bearer token}
\usage{
get_bearer()
}
\value{
string represents your bearer token, if it the environmental variable "TWITTER_BEARER" has been preset.
}
\description{
This function attempts to retrieve your bearer token from the environmental variable "TWITTER_BEARER".
The easiest way to setup this environmental variable is to use \code{set_bearer()}
and insert your bearer token to \code{.Renviron} file following the format: TWITTER_BEARER=YOURTOKENHERE.
Replace YOURTOKENHERE with your own token.
}
\details{
Note: for \code{get_bearer()} to retrieve your bearer token you will need to restart the
R session after storing in \code{.Renviron}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_functions.R
\name{get_radius_tweets}
\alias{get_radius_tweets}
\title{Get tweets within radius buffer
`r lifecycle::badge("deprecated")}
\usage{
get_radius_tweets(
  query,
  radius,
  start_tweets,
  end_tweets,
  bearer_token = get_bearer(),
  n = 100,
  file = NULL,
  data_path = NULL,
  bind_tweets = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{query}{string or character vector, search query or queries}

\item{radius}{numeric, a vector of two point coordinates latitude, longitude, and point radius distance (in miles)}

\item{start_tweets}{string, starting date}

\item{end_tweets}{string, ending date}

\item{bearer_token}{string, bearer token}

\item{n}{integer, upper limit of tweets to be fetched}

\item{file}{string, name of the resulting RDS file}

\item{data_path}{string, if supplied, fetched data can be saved to the designated path as jsons}

\item{bind_tweets}{If \code{TRUE}, tweets captured are bound into a data.frame for assignment}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{...}{arguments will be passed to \code{\link[=build_query]{build_query()}} function. See \code{?build_query()} for further information.}
}
\value{
a data.frame
}
\description{
This function collects tweets containing strings or hashtags
between specified date ranges filtering by radius buffer. Tweet-level data is stored in a data/
path as a series of JSONs beginning "data_"; User-level data is stored as a series of
JSONs beginning "users_". If a filename is supplied, this function will
save the result as a RDS file, otherwise it will return the results as a data.frame.
Note: radius must be less than 25mi.
}
\examples{
\dontrun{
bearer_token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
tweets <- get_radius_tweets("happy", radius = c(-0.131969125179604,51.50847878040284, 25),
                           start_tweets = "2021-01-01T00:00:00Z",
                           end_tweets = "2021-01-01T10:00:00Z",
                           bearer_token = bearer_token, data_path = "data/")
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_queryv2.R
\name{build_query}
\alias{build_query}
\title{Build tweet query}
\usage{
build_query(
  query = NULL,
  exact_phrase = NULL,
  users = NULL,
  reply_to = NULL,
  retweets_of = NULL,
  exclude = NULL,
  is_retweet = NULL,
  is_reply = NULL,
  is_quote = NULL,
  is_verified = NULL,
  remove_promoted = FALSE,
  has_hashtags = NULL,
  has_cashtags = NULL,
  has_links = NULL,
  has_mentions = NULL,
  has_media = NULL,
  has_images = NULL,
  has_videos = NULL,
  has_geo = NULL,
  place = NULL,
  country = NULL,
  point_radius = NULL,
  bbox = NULL,
  lang = NULL,
  conversation_id = NULL,
  url = NULL
)
}
\arguments{
\item{query}{string or character vector, search query or queries}

\item{exact_phrase}{If \code{TRUE}, only tweets will be returned matching the exact phrase}

\item{users}{string or character vector, user handles to collect tweets from the specified users}

\item{reply_to}{string or character vector, user handles to collect replies to the specified users}

\item{retweets_of}{string or character vector, user handles to collects retweets of tweets by the specified users}

\item{exclude}{string or character vector, tweets containing the keyword(s) will be excluded}

\item{is_retweet}{If \code{TRUE}, only retweets will be returned; if \code{FALSE}, retweets will be excluded; if \code{NULL}, both retweets and other tweet types will be returned.}

\item{is_reply}{If \code{TRUE}, only replies will be returned; if \code{FALSE}, replies will be excluded; if \code{NULL}, both replies and other tweet types will be returned.}

\item{is_quote}{If \code{TRUE}, only quote tweets will be returned; if \code{FALSE}, quote tweets will be excluded; if \code{NULL}, both quote tweets and other tweet types will be returned.}

\item{is_verified}{If \code{TRUE}, only tweets from verified accounts will be returned; if \code{FALSE}, tweets from verified accounts will be excluded; if \code{NULL}, both verified account tweets and tweets from non-verified accounts will be returned.}

\item{remove_promoted}{If \code{TRUE}, tweets created for promotion only on ads.twitter.com are removed}

\item{has_hashtags}{If \code{TRUE}, only tweets containing hashtags will be returned; if \code{FALSE}, tweets containing hashtags will be excluded; if \code{NULL}, both tweets containing hashtags and tweets without hashtags will be returned.}

\item{has_cashtags}{If \code{TRUE}, only tweets containing cashtags will be returned; if \code{FALSE}, tweets containing cashtags will be excluded; if \code{NULL}, both tweets containing cashtags and tweets without cashtags will be returned.}

\item{has_links}{If \code{TRUE}, only tweets containing links (and media) will be returned; if \code{FALSE}, tweets containing links (and media) will be excluded; if \code{NULL}, both tweets containing links (and media) and tweets without links (and media) will be returned.}

\item{has_mentions}{If \code{TRUE}, only tweets containing mentions will be returned; if \code{FALSE}, tweets containing mentions will be excluded; if \code{NULL}, both tweets containing mentions and tweets without mentions will be returned.}

\item{has_media}{If \code{TRUE}, only tweets containing media such as a photo, GIF, or video (as determined by Twitter) will be returned will be returned; if \code{FALSE}, tweets containing media will be excluded; if \code{NULL}, both tweets containing media and tweets without media will be returned.}

\item{has_images}{If \code{TRUE}, only tweets containing (recognized URLs to) images will be returned will be returned will be returned; if \code{FALSE}, tweets containing images will be excluded; if \code{NULL}, both tweets containing images and tweets without images will be returned.}

\item{has_videos}{If \code{TRUE},  only tweets containing contain videos (recognized as native videos uploaded directly to Twitter) will be returned will be returned; if \code{FALSE}, tweets containing videos will be excluded; if \code{NULL}, both tweets containing videos and tweets without videos will be returned.}

\item{has_geo}{If \code{TRUE}, only tweets containing geo information (Tweet-specific geolocation data provided by the Twitter user) will be returned; if \code{FALSE}, tweets containing geo information will be excluded; if \code{NULL}, both tweets containing geo information and tweets without geo information will be returned.}

\item{place}{string, name of place e.g. "London"}

\item{country}{string, name of country as ISO alpha-2 code e.g. "GB"}

\item{point_radius}{numeric, a vector of two point coordinates latitude, longitude, and point radius distance (in miles)}

\item{bbox}{numeric, a vector of four bounding box coordinates from west longitude to north latitude}

\item{lang}{string, a single BCP 47 language identifier e.g. "fr"}

\item{conversation_id}{string, return tweets that share the specified conversation ID}

\item{url}{string, url}
}
\value{
a query string
}
\description{
Build tweet query according to targeted parameters.
}
\details{
This function is already called within the main
\code{\link{get_all_tweets}} function.

It may also be called separately and the output saved as
a character object query string to be input as query parameter to \code{\link{get_all_tweets}}.
}
\examples{
\dontrun{
query <- build_query(query = "happy", is_retweet = FALSE,
                     country = "US",
                     place = "seattle",
                     point_radius = c(-122.33795253639994, 47.60900846404393, 25),
                     lang = "en")
                     
query <- build_query(query = "twitter",
                     point_radius = c(-122.33795253639994, 47.60900846404393, 25),
                     lang = "en")
                     
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_liked_tweets.R
\name{get_liked_tweets}
\alias{get_liked_tweets}
\title{Get liked tweets}
\usage{
get_liked_tweets(x, bearer_token = get_bearer(), ...)
}
\arguments{
\item{x}{string containing one user id or a vector of user ids}

\item{bearer_token}{string, bearer token}

\item{...}{arguments passed to other backend functions}
}
\value{
a data frame
}
\description{
This function fetches returns tweets liked by a user or users.
}
\examples{
\dontrun{
users <- c("2244994945", "95226101")
get_liked_tweets(users, bearer_token = get_bearer())
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_functions.R
\name{get_user_tweets}
\alias{get_user_tweets}
\title{Get tweets from user
`r lifecycle::badge("deprecated")

This function collects tweets of a user or set of users between specified date ranges.
Tweet-level data is stored in a data/ path as a series of JSONs beginning "data_"; User-level
data is stored as a series of JSONs beginning "users_". If a filename is supplied, this
function will save the result as a RDS file, otherwise it will return the results as a dataframe.}
\usage{
get_user_tweets(
  users,
  start_tweets,
  end_tweets,
  bearer_token = get_bearer(),
  n = 100,
  file = NULL,
  data_path = NULL,
  bind_tweets = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{users}{character vector, user handles from which to collect data}

\item{start_tweets}{string, starting date}

\item{end_tweets}{string, ending date}

\item{bearer_token}{string, bearer token}

\item{n}{integer, upper limit of tweets to be fetched}

\item{file}{string, name of the resulting RDS file}

\item{data_path}{string, if supplied, fetched data can be saved to the designated path as jsons}

\item{bind_tweets}{If \code{TRUE}, tweets captured are bound into a data.frame for assignment}

\item{verbose}{If \code{FALSE}, query progress messages are suppressed}

\item{...}{arguments will be passed to \code{\link[=build_query]{build_query()}} function. See \code{?build_query()} for further information.}
}
\value{
a data frame
}
\description{
Get tweets from user
`r lifecycle::badge("deprecated")

This function collects tweets of a user or set of users between specified date ranges.
Tweet-level data is stored in a data/ path as a series of JSONs beginning "data_"; User-level
data is stored as a series of JSONs beginning "users_". If a filename is supplied, this
function will save the result as a RDS file, otherwise it will return the results as a dataframe.
}
\examples{
\dontrun{
bearer_token <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
users <- c("uoessps", "spsgradschool")
get_user_tweets(users, "2020-01-01T00:00:00Z", "2020-01-05T00:00:00Z",
                bearer_token, data_path = "data/")
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_functions.R
\name{bind_tweet_jsons}
\alias{bind_tweet_jsons}
\title{Bind tweets stored as JSON files
`r lifecycle::badge("deprecated")}
\usage{
bind_tweet_jsons(data_path, verbose = TRUE)
}
\arguments{
\item{data_path}{string, file path to directory of stored tweets data saved as data_\emph{id}.json}

\item{verbose}{If \code{FALSE}, messages are suppressed}
}
\value{
a data.frame
}
\description{
Bind tweets stored as JSON files
`r lifecycle::badge("deprecated")
}
\examples{
\dontrun{
bind_tweet_jsons(data_path = "data/")
}
}
\keyword{internal}
