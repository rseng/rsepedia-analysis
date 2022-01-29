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

