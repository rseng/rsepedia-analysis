---
title: "rtweet: Collecting and analyzing Twitter data"
authors:
  - name: Michael W. Kearney
    orcid: 0000-0002-0730-4694
    affiliation: '1'
affiliations:
  - name: School of Journalism, Informatics Institute, University of Missouri
    index: '1'
date: 13 May 2019
bibliography: paper.bib
tags:
  - R
  - twitter
  - social media
  - API
---

# Statement of need

Following the [announced (2016) deprecation of the ``twitteR``
package](https://github.com/ropensci/rtweet/issues/1#issuecomment-492753003)
[@twitteR], R users seeking to interact with Twitter APIs have been encouraged
to use the ``rtweet`` package. Use of the up-to-date and actively-maintained
``rtweet`` package is especially important in light of changes to Twitter's APIs
since 2016. Most notably is the increased character limit (from 140 to 280) for
Twitter statuses [@tweetmodeextended]. In addition to providing an updated
interface with similar functionality to ``twitteR``, allowing R users to
communicate with various endpoints from Twitter's REST API, the ``rtweet``
package also provides support for communicating with Twitter's stream API.

# Summary

Interest in Twitter data continues to grow, but for many the task of actually
collecting and analyzing data via Twitter APIs remains daunting. For example, in
order to interact with Twitter's APIs users must, in addition to identifying and
digesting the relevant information from [Twitter's developer
documentation](https://developer.twitter.com), build/send/receive requests,
manage rate limits, and wrangle nested and real-time response objects into
analysis-friendly data structures. Fortunately, the ``rtweet`` R package [@rtweet]
is designed to simplify these processes, making interacting with Twitter's APIs
more accessible to a wider range of users.

The main goals of the ``rtweet`` package are two-fold. The first goal is to make
interacting with Twitter's APIs more approachable and streamlined for less
computationally-experienced users. The second goal is to assist in the analysis
of Twitter data via converting information returned by Twitter's APIs into
tabular data structures and providing several convenience functions for common
analytical techniques such as examining Twitter networks or the frequency of
tweets over time. In short, although it is certainly possible for users to write
their own Twitter API wrapper functions, the heavy-lifting done by ``rtweet`` to
(a) streamline the building, authorizing, and sending of API requests, (b)
wrangle deeply nested JSON data into tabular structures, and (c) provide
convenience functions for relevant and popular analytical techniques, make
it a valuable contribution in the area of collecting and analyzing Twitter data.

Although ``rtweet`` provides some coverage to user context-behaviors (e.g.,
posting statuses, liking tweets, following users, etc.), the primary audience
for the package to date has been researchers. Accordingly, ``rtweet`` has been
featured in numerous popular press [e.g.,
@bajak2019democrats; @machlis2019r; @riley2019twitter] and academic publications
[e.g.,
@bossetta2018simulated; @bradley2019major; @buscema2018media;
@erlandsen2018twitter; @gitto2019brand; @kearney2019analyzing;
@kearney2018analyzing; @li2019sentiment; @lutkenhaus2019tailoring;
@lutkenhaus2019mapping; @molyneux2018media; @tsoi2018can; @unsihuay2018topic;
@valls2017urban; @wu2018finding].

# References



<!-- README.md is generated from README.Rmd. Please edit that file -->

# rtweet <img src="man/figures/logo.png" width="160px" align="right" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/ropensci/rtweet/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rtweet/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/rtweet)](https://cran.r-project.org/package=rtweet)
[![Coverage
Status](https://codecov.io/gh/ropensci/rtweet/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/rtweet?branch=master)
![Downloads](https://cranlogs.r-pkg.org/badges/rtweet)
[![ZENODO](https://zenodo.org/badge/64161359.svg)](https://zenodo.org/badge/latestdoi/64161359)
[![rOpenSci](https://badges.ropensci.org/302_status.svg)](https://github.com/ropensci/software-review/issues/302)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.01829/status.svg)](https://doi.org/10.21105/joss.01829)
<!-- badges: end -->

Use twitter from R. Get started by reading `vignette("rtweet")`.

## Installation

To get the current released version from CRAN:

``` r
install.packages("rtweet")
```

## Usage

All users must be authenticated to interact with Twitter’s APIs. The
easiest way to authenticate is to use your personal twitter account -
this will happen automatically (via a browser popup) the first time you
use an rtweet function. See `auth_setup_default()` for details. Using
your personal account is fine for casual use, but if you are trying to
collect a lot of data it’s a good idea to authentication with your own
Twitter “app”. See `vignette("auth", package = "rtweet")` for details.

``` r
library(rtweet)
```

rtweet should be used in strict accordance with Twitter’s [developer
terms](https://developer.twitter.com/en/developer-terms/more-on-restricted-use-cases).

### Search tweets or users

Search for up to 10,000 tweets containing #rstats, the common hashtag
used to refer to the R language, excluding retweets:

``` r
rt <- search_tweets("#rstats", n = 10000, include_rts = FALSE)
```

Twitter rate limits cap the number of search results returned to 18,000
every 15 minutes. To request more than that, set
`retryonratelimit = TRUE` and rtweet will wait for rate limit resets for
you.

Search for 1,000 users with the #rstats in their profile:

``` r
usrs <- search_users("#rstats", n = 1000)
```

### Stream tweets

Randomly sample (approximately 1%) from the live stream of all tweets:

``` r
rt <- stream_tweets("")
```

Stream all geo-located tweets from London for 60 seconds:

``` r
rt <- stream_tweets(location = lookup_coords("london"), timeout = 60)
```

### Get friends and followers

Get all accounts followed by a user:

``` r
## get user IDs of accounts followed by R Foundation
R_Foundation_fds <- get_friends("_R_Foundation")

## lookup data on those accounts
R_Foundation_fds_data <- lookup_users(R_Foundation_fds$user_id)
```

Get all accounts following a user:

``` r
## get user IDs of accounts following R Foundation
R_Foundation_flw <- get_followers("_R_Foundation", n = 10000)
R_Foundation_flw_data <- lookup_users(R_Foundation_flw$user_id)
```

If you want *all* followers, you’ll need you’ll need to set `n = Inf`
and `retryonratelimit = TRUE` but be warned that this might take a
*long* time.

### Get timelines

Get the most recent 3,200 tweets from R Foundation:

``` r
## get user IDs of accounts followed by R Foundation
tmls <- get_timelines("_R_Foundation", n = 3200)
```

### Get favorites

Get the 3,000 most recently favorited statuses by R Foundation:

``` r
jkr <- get_favorites("_R_Foundation", n = 3000)
```

## Contact

Communicating with Twitter’s APIs relies on an internet connection,
which can sometimes be inconsistent. With that said, if you encounter an
obvious bug for which there is not already an active
[issue](https://github.com/ropensci/rtweet/issues), please [create a new
issue](https://github.com/ropensci/rtweet/issues/new) with all code used
(preferably a reproducible example) on Github.

# Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.
# rtweet (development version)

- The default value of `retryonratelimit` comes from the option
  `rtweet.retryonratelimit` so you can globally set it to `TRUE` if desired
  (#173).
  
- `auth_as()` accepts path to an authentication to make it easier to use authentication outside a user account (#602)

- All paginated functions that don't return tweets now use a consistent 
  pagination interface. They all store the "next cursor" in an `rtweet_cursor`
  attribute, which will be automatically retrieved when you use the `cursor`
  argument.
  
- Message are now properly capitalized (#565, @jsta)  

- Banned or protected accounts now trigger a warning instead of an error (#590)

- `get_friends()` and `get_followers()` return similar formatted output with 
  two columns "from_id" and "to_id" (#308)

- `lookup_users()` and `search_users()` now returns a data frame containing
  all information about each user (not their latest tweet). If you want to get 
  that data you can use `tweets_data()`.

- `parse = FALSE` always means return the raw "JSON". Previously some functions
  (e.g. `my_friendships()`) would return the raw HTTP response instead (#504).

- Functions that return tweets (e.g. `get_favorites()`, `get_my_timeline()`, 
  `get_timeline()`, `get_mentions()`, `lists_statuses()` and `search_tweets()`)
  now expose a consistent pagination interface. They all support `max_id` and 
  `since_id` to find earlier and later tweets respectively, as well as
  `retryonratelimit` to wait as long as needed when rate limited (#510).

- `lookup_collections()` and `get_collections()` has been hard deprecated 
  because the underlying Twitter API has been deprecated.

- `previous_cursor()` has been hard deprecated. It could only be used with 
  `lists_memberships()` and it has been dropped in favour of making regular
  pagination better.

- The `home` argument to `get_timeline()` has been deprecated. You can only
  retrieve the home timeline for the logged in user, and that's the job of
  `get_my_timeline()` (#550).

- All functions that perform multiple requests on your behalf now display
  a progress bar so you know what's happening. If you don't want it, you can 
  turn it off with `verbose = FALSE` (#518). 

- `lookup_statuses()` has been deprecated in favour of `lookup_tweets()`.

- `as_userid()` has been deprecated since in case of ambiguity the default is
  to assume a numeric string is a user id (#520). All functions now use a 
  single `user_type()` function so behaviour is identical for all rtweet 
  functions.

- `get_timelines()` has been deprecated since it does that same thing as
  `get_timeline()` (#509).

- rtweet no longer re-exports the pipe; if you want to continue using it, you'll
  need to `library(magrittr)` or `library(dplyr)` (#522).

- `stream_tweets()` has been overhauled to only write valid data. This obsoletes
  all previous strategy to cleen up bad data after the fact (#350, #356).

- `stream_tweets2()` and `parse_stream()` have been deprecated in favour of
  `stream_tweets()` and `jsonlite::stream_in()`.

- rtweet 1.0.0 implements a consistent strategy for handling rate limits. 
  By default, if a paginated function (i.e. a rtweet function that performs 
  multiple calls to the twitter API) is rate-limited it will return all results 
  received up to that point, along with a warning telling you how to get more
  results. Alternatively, if you want to automatically wait until the 
  rate-limit is reset, you can set `retryratelimit = TRUE`.

- The `rate_limit()` interface has been drastically simplified.

- `suggested_slugs()`, `suggested_users()`, `suggested_users_all()` have been
  removed as they stopped working when Twitter remove the suggested users
  endpoint in June 2019 (https://twittercommunity.com/t/124732).

- Added support for posting alt-text metadata with images tweeted with status 
  updated via `post_tweet()`. (#425, @hrbrmstr)
  
- Update to new rOpenSci Code of Conduct: https://ropensci.org/code-of-conduct/

## Authentication

rtweet's authentication system has been completely written. It is now based around three authentication options: `rtweet_user()`, `rtweet_app()`, and `rtweet_bot()`. Authentication no longer touches `~/.Renviron` file; instead `auth_save()` and `auth_as()` allow you to explicitly save and load authentication mechanisms from a system config directory. See `vignette("auth")` for more details.

- The httpuv package is now only suggested, since it's only needed for 
  interactive auth, and you'll be prompted to install it when needed.

- `bearer_token()` has been deprecated in favour of `rtweet_app()`, which takes 
  the bearer token found in your Twitter developer portal. `invalidate_bearer()`
  has been deprecated since this is something you should do yourself in the
  Twitter developer portal.

- `create_token()` has been deprecated in favour of the combination of
  `rtweet_user()`/`rtweet_bot()`/`rtweet_app()` + `auth_as()` + `auth_save()`.

- `get_token()` and `get_tokens()` have been deprecated in favour of 
  `auth_get()` and `auth_list()`.

# rtweet 0.7.0

- Added paper.md as part of ROpenSci submission.
- Added contributing template.
- Added explanation of requirements and usage to `bearer_token()` docs.
- Transferred repo to ropensci
- Fixed numerous typos and grammatical mistakes (thank you several pull requests)
- More robust testing setup with encrypted keys entered via travis-ci web UI
- Data parsing: various bug fixes and stability improvements
- Added extended tweet mode to `list_statuses()` endpoint

# rtweet 0.6.9
- Better tweet-validating in streaming data–interrupted statuses/broken lines 
  are now returned
- Added network-graph convenience functions `network_data()` and `network_graph()`
- Added experimental support for premium APIs

# rtweet 0.6.8
- Users can now create read-only using the built-in rtweet client!

# rtweet 0.6.7
- `lookup_coords()` now requires a Google Maps API key. It will be stored for 
  easy future use once supplied.
- Improved documentation for authentication/token creation.
- Various bug fixes and improvements.

# rtweet 0.6.6
- Added `bearer_token()` option for access to more generous rate limits.
- Fixed issues with `create_token()` when using browse-based authentication 
  method.

# rtweet 0.6.5
- Added list management functionality via `post_list()`, which now allows users
  to create and populate lists as well as delete lists on behalf of one's own
  Twitter account.
- `lists_memberships()` and now scrolls through multiple pages of results to 
  automate collection of larger numbers of lists.
- Various bug fixes and improvements.

# rtweet 0.6.4
- Added new oauth method to `create_token()` which allows for creation of token
  non-interactive sessions via accepting inputs for consumer key, consumer 
  secret (always required), oauth key, and oauth secret (optional, if supplied
  then non-browser sign method is used).
- `ts_*()` functions now offer a `tz` (timezone) argument, allowing users to 
  more easily print and plot in non-UTC time.
- Users can now delete tweets by passing the status ID (of the desired tweet to
  be deleted) to the `destroy_id` argument in `post_tweet()`
- Various bug fixes and stability improvements.

# rtweet 0.6.3
- Fixed bug in `join_rtweet()`, which omitted users who didn't have 
  available tweets.
- Various bug fixes and stability improvements.

# rtweet 0.6.2
- Added `all_suggested_users()`, which automates the collection of Twitter's
  suggested users data.
- Various bug fixes and stability improvements.
- Significant upgrades to `save_as_csv()`, including addition of new 
  `prep_as_csv()` as convenience function for flattening Twitter data frames.
- Tokens have been retooled. For at least the time being, users must 
  create a Twitter app in order to be authorized to interact with the 
  REST and stream APIs.
- Joined data: instead of returning users/tweets data with its
  complementary tweets/users data stored as an attribute, functions now
  return a joined data frame, consisting of the tweets-level data 
  joined with the newest (most recent) observation for each user 
  This means functions now return a more consistent and intuitive 
  data object where one row is always equal to one tweet. 
- Overhauled `save_as_csv()` with improved flattening and ID-preserving 
  saving methods. The function now saves a single [joined] data set as 
  well.
- Fixed major bugs in `get_favorites()` and in several `lists_*()` 
  functions.
- Tweaked date-time aggregator internals to make time-rounding more 
  precise.

# rtweet 0.6.0
- Introduced new API authorization method, which leverages an embedded
  rtweet Twitter app that is authorized locally by the user. Creating
  Twitter apps is non longer necessary. Users need only click "okay"
  to create and store their API authorization token.
- Improved parsing and line-reading internals for `stream_tweets()`
- Added `stream_tweets2()` function for more robust streaming
  method. Streams JSON files to directory and reconnects following
  premature disruptions.
- Various bug fixes and numerous documentation improvements.

# rtweet 0.5.0
- Added access to direct messages, mentions, list subscriptions, list
  users, list members, and list memberships
- Various fixes to parsing, integrating tibble for output, and
  streaming geolocation-related functions and data.
- Fixed issues with streaming and parsing streamed data.

# rtweet 0.4.9
- Functions `get_timeline()`, `get_favorites()`, `get_friends()`, and
 `get_followers()` now accept vectors of length > 1.
- Fixed bugs related to users data and its extracter, `users_data()`
- New stream parser, `stream_data()`, designed to parse files that cannot
  wholely fit into memory. `stream_data()` can now work in parallel as well.

# rtweet 0.4.8
- Support for additional APIs has been added--including APIs designed
  to return information related to lists and retweets.
- The `post_status()` function has been fixed and can now be used to
  upload media.
- Several adjustments have been made in response to various changes in
  Twitter's APIs.
- Thanks to all the great feedback on Github, numerous bug fixes
  and improvements have been included as well. In general, things
  should become a lot more stable across functions and data
  structures.

# rtweet 0.4.7
- The relatively lightweight tibble package is now a package dependency.
- Speed boosts to parsing process. It's possible to convert from json to
  data frames in parallel, but I'm not sure minimal gains are worth the
  headache. Regardless, the current version should return more data,
  more reliably, and faster.
- By default, functions now return data frames (tibbles) with recursive
  lists (e.g., the 3rd observation of `mentions_screen_name` may consist of
  4 screen names).
- To revert back to the flattened/delim object, use the `flatten()` function.
  Exporting functions such as `save_as_csv` will apply flatten by default.
- Three different sets of coordinate variables are now returned: `coords_coords`,
  `geo_coords`, and `bbox_coords` bounding box. The first two come in
  pairs of coords (a list column) and bbox_coords comes with 8
  values (longX4 latX4). This should allow users to maximize returns
  on
  geo-location data.

# rtweet 0.4.6
- More efficient iterations through pages of results.
- Added to documentation, including new package documentation domain:
  http://rtweet.info.
- Improvements made in collecting and using geo data.

# rtweet 0.4.5
- Convenience function `plain_tweets()` added for textual analysis.
- Overhaul of `ts_plot()` with improved time-aggregating method. Now a
  wrapper around `ts_data()`, deprecating `ts_filter`.

# rtweet 0.4.4
- Lots of query-building features added to search tweets, including
  ability to search by geolocation.
- Post actions now include replying to status ID.
- Other various bug fixes and speed improvements.

# rtweet 0.4.3
- Now returns tibbles (tibble is a recommended dependency)
- Various bug fixes and code improvements.

# rtweet 0.4.2
* Various bug fixes
* Integration with ggplot2 as a suggested dependency

# rtweet 0.4.1
* Fixed bugs with `mutate_coords()` and `retryonratelimit`.
* Now returns full text of tweets exceeding 140 characters. This
  change was necessary due to recent changes in Twitter's API.

# rtweet 0.4.0
* CRAN release featuring major additions to documentation and support in
  addition to new and improved functions like `ts_plot()`, `ts_filter()`
  and more!

# rtweet 0.3.96
* For dev: added package builder for better versioning and
  more frequent updates to NEWS.md file.
* Added new live streaming vignette as well as updated
  and improved tokens vignette
* Various bug fixes and improvements to tokens, parse, and plot functions.

# rtweet 0.3.93
* All interactive/posting functions have been modified with the prefix
  `post_`. This was done to clearly distinguish write functions from
  retrieval functions.
* More bug fixes and various improvements.
* The `ts_plot()` function is now more robust with more adaptive
  characteristics for variations in the number of filters, the method
  of distinguishing lines, the position of the legend, and the
  aesthetics of the themes.
* Added `ts_filter()` function which allows users to convert Twitter
  data into a time series-like data frame. Users may also provide
  filtering rules with which `ts_filter()` will subset the data as it
  converts it to multiple time series, which it then outputs as a
  long-form (tidy) data frame.

# rtweet 0.3.92

* `search_tweets` now includes `retryonratelimit` argument to
allow for searches requesting more than 18,000 tweets. This
automates what was previously possible through use of `max_id`.
* Various bug fixes and improvements to parsing and pagination-
assisting functions.
* Fixed bug in encoding with `stream_tweets`.

# rtweet 0.3.91

* Major improvements to ts_plot including SIX different
themes from which users may choose
* More parsing fixes and misc stability improvements
* Minor renaming of variables along with returning more
variables overall

# rtweet 0.3.9

* Fixes minor problems with `parse.piper` function
* More additions to plotting and data wrangling for the
purpose of plotting

# rtweet 0.3.8

* Functions by default use a new faster parser that returns more
variables
* Text analysis functions provided for convenience
* Plotting with maps
* Tidyverse consistencies

# rtweet 0.3.8

* Fixed issue with geo tracking in stream_tweets
* Various bug fixes and stability improvements

# rtweet 0.3.7

* Reworked `ts_plot` to enable different filtered time series and
an aesthetic overhaul of the plot function as well.

# rtweet 0.3.6

* Added `as_double` argument to provide flexibility in handling
id variables (as_double provides performance boost but can create
problems when printing and saving, depending on format). By default
functions will return IDs as character vectors.
* Numerous improvements made to parsing and bug fixes to lookup
and search functions.

# rtweet 0.3.5
* `clean_tweets` argument provided to allow user more control over
encoding and handling of non-ascii characters.
* Fixed issue with `search_users` and implemented several
improvements to `stream_tweets` and `plot_ts`.

# rtweet 0.3.4

* Implemented robust methods to fetch tokens (whether set as
environment variable, .httr-oauth file, or if the tokens exist
in the global environment). Functions now search for variations
in the labeling of tokens---i.e., if your token(s) are saved as
`twitter_tokens`, `twitter_token`, `tokens`, or `token`, rtweet
will find it.
* Fixed issues with parsing tweets and users data.
* Stability improvements to `search_tweets` and `stream_tweeets`

# rtweet 0.3.3

* Flattened recursive columns for more reliable parsing and various
speed enhancements

# rtweet 0.3.2

* Added built-in, encrypted tokens
* Fixed issues with tweets parsing and reading streams
* Numerous speed improvements

# rtweet 0.3.1

* `include_retweets` arg added to `search_tweets()` function.
* `user_id` class changed to double when parsed. double is significantly
faster and consumes less space. it's also capable of handling the length of
id scalars, so the only downside is truncated printing.

# rtweet 0.3.0

* New CRAN version!
* Lots of improvements to stability and entirely new functions to
play around with (see previous news updates for more info).
* Added more documentation all round, including help features, examples, and
vignette infrastructure.

# rtweet 0.2.92

* Added gzip option for `stream_tweets()`

# rtweet 0.2.91

* Added sample method for `stream_tweets()` function. By default,
the streaming query argument, `q`, is now set to an empty string,
`q = ""`, which returns a random sample of all Tweets
(pretty cool, right?).

# rtweet 0.2.9

* Added `post_tweet()` function. Users can now post tweets from their R console.

# rtweet 0.2.8

* Added `get_favorites()` function
* Update tests
* Exports tweets and users classes with show and plot methods

# rtweet 0.2.7

* Added screen_name variable for user mentions (in addition to user_id).

# rtweet 0.2.6

* Added `lookup_statuses()` function, which is the counterpart to
`lookup_users()`. Supply a vector of status IDs and return tweet data
for each status. `lookup_statuses()` is particularly powerful when
combined with other methods designed to collect older Tweets. Early
experiments with doing this all through R have turned out surprisingly
well, but packaging it in a way that makes it easy to do on other
machines is unlikely to happen in the short term.

* Removed dplyr dependencies. Everyone should install and use `dplyr`,
but for sake of parsimony, it's been removed from rtweet.

* Continued development of S4 classes and methods. Given removal of
dplyr dependencies, I've started to integrate print/show methods that
will limit the number of rows (and width of columns) when printed.
Given the amount of data returned in a relatively short period of time,
printing entire data frames quickly becomes headache-inducing.

# rtweet 0.2.5

* S4 class and methods integration

# rtweet 0.2.4

* Added new trends functions. Find what trending locations are
available with `trends_available()` and/or search for trends
worldwide or by geographical location using `get_trends()`.

* Stability improvements including integration with Travis CI and
code analysis via codecov. Token encryption method also means API
testing conducted on multiple machines and systems.

# rtweet 0.2.3

* Added new `search_users()` function! Search for users by keyword,
name, or interest and return data on the first 1000 hits.

# rtweet 0.2.2

* Output for `search_tweets()`, `stream_tweets()`, and
`get_timeline()` now consists of tweets data and contains users data
attribute.

* Output for `lookup_users()` now consists of users data and contains
tweets data attribute.

* To access users data from a tweets object or vice-versa, use
`users_data()` and `tweets_data()` functions on objects output by major 
rtweet retrieval functions.

* Updated testthat tests

# rtweet 0.2.1

* Output for `get_friends()` and `get_followers()` is now a tibble
of "ids". To retrieve next cursor value, use new `next_cursor()`
function.

* Major stability improvements via testthat tests for every major
function.

# rtweet 0.2.0

* Since previous CRAN release, numerous new features and improvements
to functions returning tweets, user data, and ids.

* Search function now optimized to return more tweets per search.

* Numerous improvements to stability, error checks, and namespace
management.

# rtweet 0.1.91

* Improvements to `get_friends` and `get_followers`. Returns list
with value (`next_cursor`) used for next page of results. When
this value is 0, all results have been returned.

* Functions `get_friends` and `get_followers` now return the list
of user ids as a tibble data table, which makes the print out much
cleaner.

# rtweet 0.1.9

* Improved scrolling methods such that `search_tweets` and
`get_timeline` should return a lot more now

* Added `parser` function to return status (tweets) AND user (users)
data frames when available. As a result, the parsed output for some
functions now comes as a list containing two data frames.

# rtweet 0.1.8

* Added `get_timeline` function that returns tweets from selected user

* Added vignettes covering tokens and search tweets

* Fixed issue with `count` argument in search and user functions

# rtweet 0.1.7

* Fixed parsing issue for return objects with omitted variables

* Added `clean_tweets` convenience function for text analysis

* More examples included in documentation.

# rtweet 0.1.6

* Added `recode_error` argument to `get_friends` function. This is
especially useful for tracking networks over time.

* Further integrated `ROAuth` methods/objects to increase
compatibility with `twitteR` authorization procedures.

* Improved token checking procedures.

# rtweet 0.1.4

* Added `NEWS.md` file

* Added `key features` and more descriptions to `README.md`.

# rtweet 0.1.3

* There are now two stable parse (convert json obj to data frame)
types. For user objects (e.g., output of `lookup_users`), there
is `parse_user`. For tweet objects (e.g., output of `search_tweets`
or `stream_tweets`), there is `parse_tweets`.

* New parse functions are now exported, so they should available
for use with compatible Twitter packages or user-defined API
request operations.

# rtweet 0.1.2

* More parsing improvements

* Added `format_date` function

* Various stability improvements

# rtweet 0.1.1

* Improvements to parse functions

# rtweet 0.1.0

* Initial release
## Test environments
* local R installation, R 3.6.1
* ubuntu 16.04 (on travis-ci), R 3.6.1
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings
# screen_name has print and [ methods

    Code
      x
    Output
      <rwteet_screen_name>
      [1] "123456"

# gives useful errors

    Code
      search_tweets(c(1:10), verbose = FALSE)
    Error <simpleError>
      length(q) == 1L is not TRUE

---

    Code
      search_tweets("stats", type = "all")
    Error <rlang_error>
      `type` must be one of "mixed", "recent", or "popular".

# get_token() and get_tokens() are deprecated

    Code
      . <- get_token()
    Warning <lifecycle_warning_deprecated>
      `get_token()` was deprecated in rtweet 1.0.0.
      Please use `auth_get()` instead.
    Code
      . <- get_tokens()
    Warning <lifecycle_warning_deprecated>
      `get_tokens()` was deprecated in rtweet 1.0.0.
      Please use `auth_get()` instead.

# create_token is deprecated

    Code
      token <- suppressMessages(create_token("my-app", "x", "x", "y", "y"))
    Warning <lifecycle_warning_deprecated>
      `create_token()` was deprecated in rtweet 1.0.0.
      See vignette('auth') for details

# Check geo-related inputs for post_tweet

    Code
      msg <- paste("test geolocated error", Sys.time())
      post_tweet(msg, lat = "x", long = 0)
    Error <simpleError>
      `lat` must be numeric.
    Code
      post_tweet(msg, lat = 0, long = "x")
    Error <simpleError>
      `long` must be numeric.
    Code
      post_tweet(msg, lat = 91, long = 0)
    Error <simpleError>
      `lat` must be between -90 and 90 degrees.
    Code
      post_tweet(msg, lat = 0, long = 181)
    Error <simpleError>
      `long` must be between -180 and 180 degrees.
    Code
      post_tweet(msg, lat = 0, long = 0, display_coordinates = "error")
    Error <simpleError>
      `display_coordinates` must be TRUE/FALSE.

# bearer token doesn't accidentally expose secrets

    Code
      rtweet_app("abc")
    Output
      <Twitter bearer token>

# find auth errors politely

    Code
      find_auth(1:10)
    Error <rlang_error>
      Unrecognised input to `auth`
    Code
      find_auth("not-present")
    Error <rlang_error>
      Can't find saved auth with name 'not-present'

# default_cached_auth() handles 0, 1, and n saved

    Code
      default_cached_auth()
    Error <rlang_error>
      No default authentication found. Please call `auth_setup_default()`

---

    Code
      default_cached_auth()
    Error <rlang_error>
      No default authentication found. Pick existing auth with:
      * auth_as('test1')
      * auth_as('test2')

# collections API has been deprecated

    Code
      lookup_collections("custom-539487832448843776")
    Error <lifecycle_error_deprecated>
      `lookup_collections()` was deprecated in rtweet 1.0.0 and is now defunct.
    Code
      get_collections(status_id = "925172982313570306")
    Error <lifecycle_error_deprecated>
      `get_collections()` was deprecated in rtweet 1.0.0 and is now defunct.

# next_cursor generates informative errors

    Code
      next_cursor(letters)
    Error <rlang_error>
      `cursor` must be a string or data frame
    Code
      next_cursor(mtcars)
    Error <rlang_error>
      `cursor` must have a `rtweet_cursor` attribute

# max_id and since_id generate informative erorrs

    Code
      max_id(10)
    Error <rlang_error>
      `max_id` must be a character vector or data frame
    Code
      max_id(mtcars)
    Error <rlang_error>
      `max_id` must contain a `id` column
    Code
      since_id(10)
    Error <rlang_error>
      `since_id` must be a character vector or data frame
    Code
      since_id(mtcars)
    Error <rlang_error>
      `since_id` must contain a `id` column

# get_timelines() is deprecated

    Code
      x <- get_timelines("cnn", n = 10)
    Warning <lifecycle_warning_deprecated>
      `get_timelines()` was deprecated in rtweet 1.0.0.
      Please use `get_timeline()` instead.

# Contributing to rtweet

This outlines how to propose a change to rtweet. For more detailed
info about contributing to this, and other tidyverse packages, please see the
[**development contributing guide**](https://rstd.io/tidy-contrib).

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the Travis and AppVeyor build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.  
*  New code should follow the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with
[Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html), 
for documentation.  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the rtweet project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See tidyverse [development contributing guide](https://rstd.io/tidy-contrib)
for further details.
<!-- This is an issue template for bugs and requests for R pkg rtweet -->

<!-- If you've encountered a likely bug in rtweet, please take a few seconds to 
look through existing issues for a similar issue. If you don't see a related 
issue, please complete the prompts below to make it easier to replicate and 
[hopefully] resolve your issue.  -->

<!-- If you have some question about how to use this package, please post it on
https://discuss.ropensci.org/tag/rtweet -->

### Problem

<!-- Succinctly describe the problem (be as specific as you think necessary) -->

### Expected behavior

<!-- Describe the behavior/result you expected -->

### Reproduce the problem

<!-- Describe and provide relevant code to reproduce the problem -->
<!-- If code doesn't always produce error, provide approximate code anyway -->

``` r
## insert code here

```

### rtweet version

<!-- run the code below and copy/paste the output -->

``` r
## copy/paste output
packageVersion("rtweet")
```


### Session info

<!-- run the code below and copy/paste the output -->

``` r
## copy/paste output
sessionInfo()
```

<!-- If you think the problem may be related to features/limitations of 
Twitter's API, you can find more information about Twitter's APIs here: 
https://developer.twitter.com/en/docs.html -->

<!-- Thank you for using and improving rtweet!  -->
# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 3.6.1 (2019-07-05) |
|os       |Ubuntu 18.04.3 LTS           |
|system   |x86_64, linux-gnu            |
|ui       |RStudio                      |
|language |(EN)                         |
|collate  |en_US.UTF-8                  |
|ctype    |en_US.UTF-8                  |
|tz       |America/Chicago              |
|date     |2019-12-19                   |

# Dependencies

|package     |old      |new      |Δ  |
|:-----------|:--------|:--------|:--|
|rtweet      |0.6.9    |0.7.0    |*  |
|askpass     |1.1      |1.1      |   |
|assertthat  |0.2.1    |0.2.1    |   |
|backports   |1.1.5    |1.1.5    |   |
|BH          |1.72.0-2 |1.72.0-2 |   |
|cli         |2.0.0    |2.0.0    |   |
|crayon      |1.3.4    |1.3.4    |   |
|curl        |4.3      |4.3      |   |
|digest      |0.6.23   |0.6.23   |   |
|ellipsis    |0.3.0    |0.3.0    |   |
|fansi       |0.4.0    |0.4.0    |   |
|glue        |1.3.1    |1.3.1    |   |
|hms         |0.5.2    |0.5.2    |   |
|httpuv      |1.5.2    |1.5.2    |   |
|httr        |1.4.1    |1.4.1    |   |
|jsonlite    |1.6      |1.6      |   |
|later       |1.0.0    |1.0.0    |   |
|magrittr    |1.5      |1.5      |   |
|mime        |0.8      |0.8      |   |
|openssl     |1.4.1    |1.4.1    |   |
|pillar      |1.4.2    |1.4.2    |   |
|pkgconfig   |2.0.3    |2.0.3    |   |
|prettyunits |1.0.2    |1.0.2    |   |
|progress    |1.2.2    |1.2.2    |   |
|promises    |1.1.0    |1.1.0    |   |
|R6          |2.4.1    |2.4.1    |   |
|Rcpp        |1.0.3    |1.0.3    |   |
|rlang       |0.4.2    |0.4.2    |   |
|sys         |3.3      |3.3      |   |
|tibble      |2.1.3    |2.1.3    |   |
|utf8        |1.1.4    |1.1.4    |   |
|vctrs       |0.2.1    |0.2.1    |   |
|zeallot     |0.1.0    |0.1.0    |   |

# Revdeps

## Failed to check (2)

|package   |version |error |warning |note |
|:---------|:-------|:-----|:-------|:----|
|bdpar     |?       |      |        |     |
|VOSONDash |?       |      |        |     |

*Wow, no problems at all. :)*# bdpar

<details>

* Version: 
* Source code: ???
* URL: https://CRAN.R-project.org/package=rtweet
* BugReports: https://github.com/ropensci/rtweet/issues
* Number of recursive dependencies: 0

Run `revdep_details(,"")` for more info

</details>

## Error before installation

### Devel

```






```
### CRAN

```






```
# VOSONDash

<details>

* Version: 
* Source code: ???
* URL: https://CRAN.R-project.org/package=rtweet
* BugReports: https://github.com/ropensci/rtweet/issues
* Number of recursive dependencies: 0

Run `revdep_details(,"")` for more info

</details>

## Error before installation

### Devel

```






```
### CRAN

```






```
# rtweet 0.6.9
- Better tweet-validating in streaming data–interrupted statuses/broken lines 
  are now returned
- Added network-graph convenience functions `network_data()` and `network_graph()`
- Added experimental support for premium APIs

# rtweet 0.6.8
- Users can now create read-only using the built-in rtweet client!

# rtweet 0.6.7
- `lookup_coords()` now requires a Google Maps API key. It will be stored for 
  easy future use once supplied.
- Improved documentation for authentication/token creation.
- Various bug fixes and improvements.

# rtweet 0.6.6
- Added `bearer_token()` option for access to more generous rate limits.
- Fixed issues with `create_token()` when using browse-based authentication 
  method.

# rtweet 0.6.5
- Added list management functionality via `post_list()`, which now allows users
  to create and populate lists as well as delete lists on behalf of one's own
  Twitter account.
- `lists_memberships()` and now scrolls through multiple pages of results to 
  automate collection of larger numbers of lists.
- Various bug fixes and improvements.

# rtweet 0.6.4
- Added new oauth method to `create_token()` which allows for creation of token
  non-interactive sessions via accepting inputs for consumer key, consumer 
  secret (always required), oauth key, and oauth secret (optional, if supplied
  then non-browser sign method is used).
- `ts_*()` functions now offer a `tz` (timezone) argument, allowing users to 
  more easily print and plot in non-UTC time.
- Users can now delete tweets by passing the status ID (of the desired tweet to
  be deleted) to the `destroy_id` argument in `post_tweet()`
- Various bug fixes and stability improvements.

# rtweet 0.6.3
- Fixed bug in `join_rtweet()`, which omitted users who didn't have 
  available tweets.
- Various bug fixes and stability improvements.

# rtweet 0.6.2
- Added `all_suggested_users()`, which automates the collection of Twitter's
  suggested users data.
- Various bug fixes and stability improvements.
- Significant upgrades to `save_as_csv()`, including addition of new 
  `prep_as_csv()` as convenience function for flattening Twitter data frames.
- Tokens have been retooled. For at least the time being, users must 
  create a Twitter app in order to be authorized to interact with the 
  REST and stream APIs.
- Joined data: instead of returning users/tweets data with its
  complementary tweets/users data stored as an attribute, functions now
  return a joined data frame, consisting of the tweets-level data 
  joined with the newest (most recent) observation for each user 
  This means functions now return a more consistent and intuitive 
  data object where one row is always equal to one tweet. 
- Overhauled `save_as_csv()` with improved flattening and ID-preserving 
  saving methods. The function now saves a single [joined] data set as 
  well.
- Fixed major bugs in `get_favorites()` and in several `lists_*()` 
  functions.
- Tweaked date-time aggregator internals to make time-rounding more 
  precise.

# rtweet 0.6.0
- Introduced new API authorization method, which leverages an embedded
  rtweet Twitter app that is authorized locally by the user. Creating
  Twitter apps is non longer necessary. Users need only click "okay"
  to create and store their API authorization token.
- Improved parsing and line-reading internals for `stream_tweets()`
- Added `stream_tweets2()` function for more robust streaming
  method. Streams JSON files to directory and reconnects following
  premature disruptions.
- Various bug fixes and numerous documentation improvements.

# rtweet 0.5.0
- Added access to direct messages, mentions, list subscriptions, list
  users, list members, and list memberships
- Various fixes to parsing, integrating tibble for output, and
  streamling geolocation-related functions and data.
- Fixed issues with streaming and parsing streamed data.

# rtweet 0.4.9
- Functions `get_timeline()`, `get_favorites()`, `get_friends()`, and
 `get_followers()` now accept vectors of length > 1.
- Fixed bugs related to users data and its extracter, `users_data()`
- New stream parser, `stream_data()`, designed to parse files that cannot
  wholely fit into memory. `stream_data()` can now work in parallel as well.

# rtweet 0.4.8
- Support for additional APIs has been added--including APIs designed
  to return information related to lists and retweets.
- The `post_status()` function has been fixed and can now be used to
  upload media.
- Several adjustments have been made in response to various changes in
  Twitter's APIs.
- Thanks to all the great feedback on Github, numerous bug fixes
  and improvements have been included as well. In general, things
  should become a lot more stable across functions and data
  structures.

# rtweet 0.4.7
- The relatively lightweight tibble package is now a package dependency.
- Speed boosts to parsing process. It's possible to convert from json to
  data frames in parallel, but I'm not sure minimal gains are worth the
  headache. Regardless, the current version should return more data,
  more reliably, and faster.
- By default, functions now return data frames (tibbles) with recursive
  lists (e.g., the 3rd observation of `mentions_screen_name` may consist of
  4 screen names).
- To revert back to the flattened/delim object, use the `flatten()` function.
  Exporting functions such as `save_as_csv` will apply flatten by default.
- Three different sets of coordinate variables are now returned: `coords_coords`,
  `geo_coords`, and `bbox_coords` bounding box. The first two come in
  pairs of coords (a list column) and bbox_coords comes with 8
  values (longX4 latX4). This should allow users to maximize returns
  on
  geo-location data.

# rtweet 0.4.6
- More efficient iterations through pages of results.
- Added to documentation, including new package documentation domain:
  http://rtweet.info.
- Improvements made in collecting and using geo data.

# rtweet 0.4.5
- Convenience function `plain_tweets()` added for textual analysis.
- Overhaul of `ts_plot()` with improved time-aggregating method. Now a
  wrapper around `ts_data()`, deprecating `ts_filter`.

# rtweet 0.4.4
- Lots of query-building features added to search tweets, including
  ability to search by geolocation.
- Post actions now include replying to status ID.
- Other various bug fixes and speed improvements.

# rtweet 0.4.3
- Now returns tibbles (tibble is a recommended dependency)
- Various bug fixes and code improvements.

# rtweet 0.4.2
* Various bug fixes
* Integration with ggplot2 as a suggested dependency

# rtweet 0.4.1
* Fixed bugs with `mutate_coords()` and `retryonratelimit`.
* Now returns full text of tweets exceeding 140 characters. This
  change was necessary due to recent changes in Twitter's API.

# rtweet 0.4.0
* CRAN release featuring major additions to documentation and support in
  addition to new and improved functions like `ts_plot()`, `ts_filter()`
  and more!

# rtweet 0.3.96
* For dev: added package builder for better versioning and
  more frequent updates to NEWS.md file.
* Added new live streaming vignette as well as updated
  and improved tokens vignette
* Various bug fixes and improvements to tokens, parse, and plot functions.

# rtweet 0.3.93
* All interactive/posting functions have been modified with the prefix
  `post_`. This was done to clearly distinguish write functions from
  retrieval functions.
* More bug fixes and various improvements.
* The `ts_plot()` function is now more robust with more adaptive
  characteristics for variations in the number of filters, the method
  of distinguishing lines, the position of the legend, and the
  aesthetics of the themes.
* Added `ts_filter()` function which allows users to convert Twitter
  data into a time series-like data frame. Users may also provide
  filtering rules with which `ts_filter()` will subset the data as it
  converts it to multiple time series, which it then outputs as a
  long-form (tidy) data frame.

# rtweet 0.3.92

* `search_tweets` now includes `retryonratelimit` argument to
allow for searches requesting more than 18,000 tweets. This
automates what was previously possible through use of `max_id`.
* Various bug fixes and improvements to parsing and pagination-
assisting functions.
* Fixed bug in encoding with `stream_tweets`.

# rtweet 0.3.91

* Major improvements to ts_plot including SIX different
themes from which users may choose
* More parsing fixes and misc stability improvements
* Minor renaming of variables along with returning more
variables overall

# rtweet 0.3.9

* Fixes minor problems with `parse.piper` function
* More additions to plotting and data wrangling for the
purpose of plotting

# rtweet 0.3.8

* Functions by default use a new faster parser that returns more
variables
* Text analysis functions provided for convenience
* Plotting with maps
* Tidyverse consistencies

# rtweet 0.3.8

* Fixed issue with geo tracking in stream_tweets
* Various bug fixes and stability improvements

# rtweet 0.3.7

* Reworked `ts_plot` to enable different filtered time series and
an aesthetic overhaul of the plot function as well.

# rtweet 0.3.6

* Added `as_double` argument to provide flexibility in handling
id variables (as_double provides performance boost but can create
problems when printing and saving, depending on format). By default
functions will return IDs as character vectors.
* Numerous improvements made to parsing and bug fixes to lookup
and search functions.

# rtweet 0.3.5
* `clean_tweets` argument provided to allow user more control over
encoding and handling of non-ascii characters.
* Fixed issue with `search_users` and implemented several
improvements to `stream_tweets` and `plot_ts`.

# rtweet 0.3.4

* Implemented robust methods to fetch tokens (whether set as
environment variable, .httr-oauth file, or if the tokens exist
in the global environment). Functions now search for variations
in the labeling of tokens---i.e., if your token(s) are saved as
`twitter_tokens`, `twitter_token`, `tokens`, or `token`, rtweet
will find it.
* Fixed issues with parsing tweets and users data.
* Stability improvements to `search_tweets` and `stream_tweeets`

# rtweet 0.3.3

* Flattened recursive columns for more reliable parsing and various
speed enhancements

# rtweet 0.3.2

* Added built-in, encrypted tokens
* Fixed issues with tweets parsing and reading streams
* Numerous speed improvements

# rtweet 0.3.1

* `include_retweets` arg added to `search_tweets()` function.
* `user_id` class changed to double when parsed. double is significantly
faster and consumes less space. it's also capable of handling the length of
id scalars, so the only downside is truncated printing.

# rtweet 0.3.0

* New CRAN version!
* Lots of improvements to stability and entirely new functions to
play around with (see previous news updates for more info).
* Added more documentation all round, including help features, examples, and
vignette infrastructure.

# rtweet 0.2.92

* Added gzip option for `stream_tweets()`

# rtweet 0.2.91

* Added sample method for `stream_tweets()` function. By default,
the streaming query argument, `q`, is now set to an empty string,
`q = ""`, which returns a random sample of all Tweets
(pretty cool, right?).

# rtweet 0.2.9

* Added `post_tweet()` function. Users can now post tweets from their R console.

# rtweet 0.2.8

* Added `get_favorites()` function
* Update tests
* Exports tweets and users classes with show and plot methods

# rtweet 0.2.7

* Added screen_name variable for user mentions (in addition to user_id).

# rtweet 0.2.6

* Added `lookup_statuses()` function, which is the counterpart to
`lookup_users()`. Supply a vector of status IDs and return tweet data
for each status. `lookup_statuses()` is particularly powerful when
combined with other methods designed to collect older Tweets. Early
experiments with doing this all through R have turned out surprisingly
well, but packaging it in a way that makes it easy to do on other
machines is unlikely to happen in the short term.

* Removed dplyr dependencies. Everyone should install and use `dplyr`,
but for sake of parsimony, it's been removed from rtweet.

* Continued development of S4 classes and methods. Given removal of
dplyr dependencies, I've started to integrate print/show methods that
will limit the number of rows (and width of columns) when printed.
Given the amount of data returned in a relatively short period of time,
printing entire data frames quickly becomes headache-inducing.

# rtweet 0.2.5

* S4 class and methods integration

# rtweet 0.2.4

* Added new trends functions. Find what trending locations are
available with `trends_available()` and/or search for trends
worldwide or by geographical location using `get_trends()`.

* Stability improvements including integration with Travis CI and
code analysis via codecov. Token encryption method also means API
testing conducted on multiple machines and systems.

# rtweet 0.2.3

* Added new `search_users()` function! Search for users by keyword,
name, or interest and return data on the first 1000 hits.

# rtweet 0.2.2

* Output for `search_tweets()`, `stream_tweets()`, and
`get_timeline()` now consists of tweets data and contains users data
attribute.

* Output for `lookup_users()` now consists of users data and contains
tweets data attribute.

* To access users data from a tweets object or vice-versa, use
`users_data()` and `tweets_data()` functions on objects output by major 
rtweet retrieval functions.

* Updated testthat tests

# rtweet 0.2.1

* Output for `get_friends()` and `get_followers()` is now a tibble
of "ids". To retrieve next cursor value, use new `next_cursor()`
function.

* Major stability improvements via testthat tests for every major
function.

# rtweet 0.2.0

* Since previous CRAN release, numerous new features and improvements
to functions returning tweets, user data, and ids.

* Search function now optimized to return more tweets per search.

* Numerous improvements to stability, error checks, and namespace
management.

# rtweet 0.1.91

* Improvements to `get_friends` and `get_followers`. Returns list
with value (`next_cursor`) used for next page of results. When
this value is 0, all results have been returned.

* Functions `get_friends` and `get_followers` now return the list
of user ids as a tibble data table, which makes the print out much
cleaner.

# rtweet 0.1.9

* Improved scrolling methods such that `search_tweets` and
`get_timeline` should return a lot more now

* Added `parser` function to return status (tweets) AND user (users)
data frames when available. As a result, the parsed output for some
functions now comes as a list containing two data frames.

# rtweet 0.1.8

* Added `get_timeline` function that returns tweets from selected user

* Added vignettes covering tokens and search tweets

* Fixed issue with `count` argument in search and user functions

# rtweet 0.1.7

* Fixed parsing issue for return objects with omitted variables

* Added `clean_tweets` convenience function for text analysis

* More examples included in documentation.

# rtweet 0.1.6

* Added `recode_error` argument to `get_friends` function. This is
especially useful for tracking networks over time.

* Further integrated `ROAuth` methods/objects to increase
compatibility with `twitteR` authorization procedures.

* Improved token checking procedures.

# rtweet 0.1.4

* Added `NEWS.md` file

* Added `key features` and more descriptions to `README.md`.

# rtweet 0.1.3

* There are now two stable parse (convert json obj to data frame)
types. For user objects (e.g., output of `lookup_users`), there
is `parse_user`. For tweet objects (e.g., output of `search_tweets`
or `stream_tweets`), there is `parse_tweets`.

* New parse functions are now exported, so they should available
for use with compatible Twitter packages or user-defined API
request operations.

# rtweet 0.1.2

* More parsing improvements

* Added `format_date` function

* Various stability improvements

# rtweet 0.1.1

* Improvements to parse functions

# rtweet 0.1.0

* Initial release
# rtweet 0.6.8
- Users can now create read-only using the built-in rtweet client!

# rtweet 0.6.7
- `lookup_coords()` now requires a Google Maps API key. It will be stored for 
  easy future use once supplied.
- Improved documentation for authentication/token creation.
- Various bug fixes and improvements.

# rtweet 0.6.6
- Added `bearer_token()` option for access to more generous rate limits.
- Fixed issues with `create_token()` when using browse-based authentication 
  method.

# rtweet 0.6.5
- Added list management functionality via `post_list()`, which now allows users
  to create and populate lists as well as delete lists on behalf of one's own
  Twitter account.
- `lists_memberships()` and now scrolls through multiple pages of results to 
  automate collection of larger numbers of lists.
- Various bug fixes and improvements.

# rtweet 0.6.4
- Added new oauth method to `create_token()` which allows for creation of token
  non-interactive sessions via accepting inputs for consumer key, consumer 
  secret (always required), oauth key, and oauth secret (optional, if supplied
  then non-browser sign method is used).
- `ts_*()` functions now offer a `tz` (timezone) argument, allowing users to 
  more easily print and plot in non-UTC time.
- Users can now delete tweets by passing the status ID (of the desired tweet to
  be deleted) to the `destroy_id` argument in `post_tweet()`
- Various bug fixes and stability improvements.

# rtweet 0.6.3
- Fixed bug in `join_rtweet()`, which omitted users who didn't have 
  available tweets.
- Various bug fixes and stability improvements.

# rtweet 0.6.2
- Added `all_suggested_users()`, which automates the collection of Twitter's
  suggested users data.
- Various bug fixes and stability improvements.
- Significant upgrades to `save_as_csv()`, including addition of new 
  `prep_as_csv()` as convience function for flattening Twitter data frames.
- Tokens have been retooled. For at least the time being, users must 
  create a Twitter app in order to be authorized to interact with the 
  REST and stream APIs.
- Joined data: instead of returning users/tweets data with its
  complementary tweets/users data stored as an attribute, functions now
  return a joined data frame, consisting of the tweets-level data 
  joined with the newest (most recent) observation for each user 
  This means functions now return a more consistent and intuitive 
  data object where one row is always equal to one tweet. 
- Overhauled `save_as_csv()` with improved flattening and ID-preserving 
  saving methods. THe function now saves a single [joined] data set as 
  well.
- Fixed major bugs in `get_favorites()` and in several `lists_*()` 
  functions.
- Tweaked date-time aggregator internals to make time-rounding more 
  precise.

# rtweet 0.6.0
- Introduced new API authorization method, which leverages an embedded
  rtweet Twitter app that is authorized locally by the user. Creating
  Twitter apps is non longer necessary. Users need only click "okay"
  to create and store their API authorization token.
- Improved parsing and line-reading internals for `stream_tweets()`
- Added `stream_tweets2()` function for more robust streaming
  method. Streams JSON files to directory and reconnects following
  premature disruptions.
- Various bug fixes nad numerous documentation improvements.

# rtweet 0.5.0
- Added access to direct messages, mentions, list subscriptions, list
  users, list members, and list memberships
- Various fixes to parsing, integrating tibble for output, and
  streamling geolocation-related functions and data.
- Fixed issues with streaming and parsing streamed data.

# rtweet 0.4.9
- Functions `get_timeline()`, `get_favorites()`, `get_friends()`, and
 `get_followers()` now accept vectors of length > 1.
- Fixed bugs related to users data and its extracter, `users_data()`
- New stream parser, `stream_data()`, designed to parse files that cannot
  wholely fit into memory. `stream_data()` can now work in parallel as well.

# rtweet 0.4.8
- Support for additional APIs has been added--including APIs designed
  to return information related to lists and retweets.
- The `post_status()` function has been fixed and can now be used to
  upload media.
- Several adjustments have been made in response to various changes in
  Twitter's APIs.
- Thanks to all the great feedback on Github, numerous bug fixes
  and improvements have been included as well. In general, things
  should become a lot more stable across functions and data
  structures.

# rtweet 0.4.7
- The relatively lightweight tibble package is now a package dependency.
- Speed boosts to parsing process. It's possible to convert from json to
  data frames in parallel, but I'm not sure minimal gains are worth the
  headache. Regardless, the current version should return more data,
  more reliably, and faster.
- By default, functions now return data frames (tibbles) with recursive
  lists (e.g., the 3rd observation of `mentions_screen_name` may consist of
  4 screen names).
- To revert back to the flattened/delim object, use the `flatten()` function.
  Exporting functions such as `save_as_csv` will apply flatten by default.
- Three different sets of coordinate variables are now returned: `coords_coords`,
  `geo_coords`, and `bbox_coords` bounding box. The first two come in
  pairs of coords (a list column) and bbox_coords comes with 8
  values (longX4 latX4). This should allow users to maximize returns
  on
  geo-location data.

# rtweet 0.4.6
- More efficient iterations through pages of results.
- Added to documentation, including new package documentation domain:
  http://rtweet.info.
- Improvements made in collecting and using geo data.

# rtweet 0.4.5
- Convenience function `plain_tweets()` added for textual analysis.
- Overhaul of `ts_plot()` with improved time-aggregating method. Now a
  wrapper around `ts_data()`, deprecating `ts_filter`.

# rtweet 0.4.4
- Lots of query-building features added to search tweets, including
  ability to search by geolocation.
- Post actions now include replying to status ID.
- Other various bug fixes and speed improvements.

# rtweet 0.4.3
- Now returns tibbles (tibble is a recommended dependency)
- Various bug fixes and code improvements.

# rtweet 0.4.2
* Various bug fixes
* Integration with ggplot2 as a suggested dependency

# rtweet 0.4.1
* Fixed bugs with `mutate_coords()` and `retryonratelimit`.
* Now returns full text of tweets exceeding 140 characters. This
  change was necessary due to recent changes in Twitter's API.

# rtweet 0.4.0
* CRAN release featuring major additions to documentation and support in
  addition to new and improved functions like `ts_plot()`, `ts_filter()`
  and more!

# rtweet 0.3.96
* For dev: added package builder for better versioning and
  more frequent updates to NEWS.md file.
* Added new live streaming vignette as well as updated
  and improved tokens vignette
* Various bug fixes and improvements to tokens, parse, and plot functions.

# rtweet 0.3.93
* All interactive/posting functions have been modified with the prefix
  `post_`. This was done to clearly distinguish write functions from
  retrieval functions.
* More bug fixes and various improvements.
* The `ts_plot()` function is now more robust with more adaptive
  characteristics for variations in the number of filters, the method
  of distiguishing lines, the position of the legend, and the
  aesthetics of the themes.
* Added `ts_filter()` function which allows users to convert Twitter
  data into a time series-like data frame. Users may also provide
  filtering rules with which `ts_filter()` will subset the data as it
  converts it to multiple time series, which it then outputs as a
  long-form (tidy) data frame.

# rtweet 0.3.92

* `search_tweets` now includes `retryonratelimit` argument to
allow for searches requesting more than 18,000 tweets. This
automates what was previously possible through use of `max_id`.
* Various bug fixes and improvements to parsing and pagination-
assisting functions.
* Fixed bug in encoding with `stream_tweets`.

# rtweet 0.3.91

* Major improvements to ts_plot including SIX different
themes from which users may choose
* More parsing fixes and misc stability improvements
* Minor renamig of variables along with returning more
variables overall

# rtweet 0.3.9

* Fixes minor problems with `parse.piper` function
* More additions to plotting and data wrangling for the
purpose of plotting

# rtweet 0.3.8

* Functions by default use a new faster parser that returns more
variables
* Text analysis functions provided for convenience
* Plotting with maps
* Tidyverse consistencies

# rtweet 0.3.8

* Fixed issue with geo tracking in stream_tweets
* Various bug fixes and stability improvements

# rtweet 0.3.7

* Reworked `ts_plot` to enable different filtered time series and
an aesthetic overhaul of the plot function as well.

# rtweet 0.3.6

* Added `as_double` argument to provide flexibility in handling
id variables (as_double provides performance boost but can create
problems when printing and saving, depending on format). By default
functions will return IDs as character vectors.
* Numerous improvements made to parsing and bug fixes to lookup
and search functions.

# rtweet 0.3.5
* `clean_tweets` argument provided to allow user more control over
encoding and handling of non-ascii characters.
* Fixed issue with `search_users` and implemented several
improvements to `stream_tweets` and `plot_ts`.

# rtweet 0.3.4

* Implemented robust methods to fetch tokens (whether set as
environment variable, .httr-oauth file, or if the tokens exist
in the global environment). Functions now search for variations
in the labeling of tokens---i.e., if your token(s) are saved as
`twitter_tokens`, `twitter_token`, `tokens`, or `token`, rtweet
will find it.
* Fixed issues with parsing tweets and users data.
* Stability improvements to `search_tweets` and `stream_tweeets`

# rtweet 0.3.3

* Flattened recursive columns for more reliable parsing and various
speed enhancements

# rtweet 0.3.2

* Added built-in, encrypted tokens
* Fixed issues with tweets parsing and reading streams
* Numerous speed improvements

# rtweet 0.3.1

* `include_retweets` arg added to `search_tweets()` function.
* `user_id` class changed to double when parsed. double is significantly
faster and consumes less space. it's also capable of handling the length of
id scalars, so the only downside is truncated printing.

# rtweet 0.3.0

* New CRAN version!
* Lots of improvements to stability and entirely new functions to
play around with (see previous news updates for more info).
* Added more documentation all round, including help features, examples, and
vignette infrastructure.

# rtweet 0.2.92

* Added gzip option for `stream_tweets()`

# rtweet 0.2.91

* Added sample method for `stream_tweets()` function. By default,
the streaming query argument, `q`, is now set to an empty string,
`q = ""`, which returns a random sample of all Tweets
(pretty cool, right?).

# rtweet 0.2.9

* Added `post_tweet()` function. Users can now post tweets from their R console.

# rtweet 0.2.8

* Added `get_favorites()` function
* Update tests
* Exports tweets and users classes with show and plot methods

# rtweet 0.2.7

* Added screen_name variable for user mentions (in addition to user_id).

# rtweet 0.2.6

* Added `lookup_statuses()` function, which is the counterpart to
`lookup_users()`. Supply a vector of status IDs and return tweet data
for each status. `lookup_statuses()` is particularly powerful when
combined with other methods designed to collect older Tweets. Early
experiments with doing this all through R have turned out surprisingly
well, but packaging it in a way that makes it easy to do on other
machines is unlikely to happen in the short term.

* Removed dplyr dependencies. Everyone should install and use `dplyr`,
but for sake of parsimony, it's been removed from rtweet.

* Continued development of S4 classes and methods. Given removal of
dplyr dependencies, I've started to integrate print/show methods that
will limit the number of rows (and width of columns) when printed.
Given the amount of data returned in a relatively short period of time,
printing entire data frames quickly becomes headache-inducing.

# rtweet 0.2.5

* S4 class and methods integration

# rtweet 0.2.4

* Added new trends functions. Find what trending locations are
available with `trends_available()` and/or search for trends
worldwide or by geogaphical location using `get_trends()`.

* Stability improvements including integration with Travis CI and
code analysis via codecov. Token encryption method also means API
testing conducted on multiple machines and systems.

# rtweet 0.2.3

* Added new `search_users()` function! Search for users by keyword,
name, or interest and return data on the first 1000 hits.

# rtweet 0.2.2

* Output for `search_tweets()`, `stream_tweets()`, and
`get_timeline()` now consists of tweets data and contains users data
attribute.

* Output for `lookup_users()` now consists of users data and contains
tweets data attribute.

* To access users data from a tweets object or vice-versa, use
`users_data()` and `tweets_data()` functions on objects outputed
by major rtweet retrieval functions.

* Updated testthat tests

# rtweet 0.2.1

* Output for `get_friends()` and `get_followers()` is now a tibble
of "ids". To retrieve next cursor value, use new `next_cursor()`
function.

* Major stability improvements via testthat tests for every major
function.

# rtweet 0.2.0

* Since previous CRAN release, numerous new features and improvements
to functions returning tweets, user data, and ids.

* Search function now optimized to return more tweets per search.

* Numerous improvements to stability, error checks, and namespace
management.

# rtweet 0.1.91

* Improvements to `get_friends` and `get_followers`. Returns list
with value (`next_cursor`) used for next page of results. When
this value is 0, all results have been returned.

* Functions `get_friends` and `get_followers` now return the list
of user ids as a tibble data table, which makes the print out much
cleaner.

# rtweet 0.1.9

* Improved scrolling methods such that `search_tweets` and
`get_timeline` should return a lot more now

* Added `parser` function to return status (tweets) AND user (users)
data frames when available. As a result, the parsed output for some
functions now comes as a list containing two data frames.

# rtweet 0.1.8

* Added `get_timeline` function that returns tweets from selected user

* Added vignettes covering tokens and search tweets

* Fixed issue with `count` argument in search and user functions

# rtweet 0.1.7

* Fixed parsing issue for return objects with omitted variables

* Added `clean_tweets` convenience function for text analysis

* More examples included in documentation.

# rtweet 0.1.6

* Added `recode_error` argument to `get_friends` function. This is
especially useful for tracking networks over time.

* Further integrated `ROAuth` methods/objects to increase
compatibility with `twitteR` authorization procedures.

* Improved token checking procedures.

# rtweet 0.1.4

* Added `NEWS.md` file

* Added `key features` and more descriptions to `README.md`.

# rtweet 0.1.3

* There are now two stable parse (convert json obj to data frame)
types. For user objects (e.g., output of `lookup_users`), there
is `parse_user`. For tweet objects (e.g., output of `search_tweets`
or `stream_tweets`), there is `parse_tweets`.

* New parse functions are now exported, so they should available
for use with compatible Twitter packages or user-defined API
request operations.

# rtweet 0.1.2

* More parsing improvements

* Added `format_date` function

* Various stability improvements

# rtweet 0.1.1

* Improvements to parse functions

# rtweet 0.1.0

* Initial release
Tests and Coverage
================
17 January, 2019 12:18:36

This output is created by
[covrpage](https://github.com/metrumresearchgroup/covrpage).

## Coverage

Coverage summary is created using the
[covr](https://github.com/r-lib/covr) package.

| Object                                               | Coverage (%) |
| :--------------------------------------------------- | :----------: |
| carbonate                                            |    48.58     |
| [R/carbonate.R](../R/carbonate.R)                    |     0.00     |
| [R/selenium\_functions.R](../R/selenium_functions.R) |     0.00     |
| [R/uri\_functions.R](../R/uri_functions.R)           |    52.63     |
| [R/carbon.R](../R/carbon.R)                          |    64.29     |
| [R/helpers.R](../R/helpers.R)                        |    70.09     |
| [R/set\_get\_functions.R](../R/set_get_functions.R)  |    100.00    |

<br>

## Unit Tests

Unit Test summary is created using the
[testthat](https://github.com/r-lib/testthat)
package.

| file                              | n |  time | error | failed | skipped | warning | icon |
| :-------------------------------- | -: | ----: | ----: | -----: | ------: | ------: | :--- |
| [test-set.R](testthat/test-set.R) | 3 | 0.013 |     0 |      0 |       0 |       0 |      |
| [test-uri.R](testthat/test-uri.R) | 8 | 1.352 |     0 |      0 |       1 |       0 | 🔶    |
| [test-yml.R](testthat/test-yml.R) | 6 | 0.017 |     0 |      0 |       0 |       0 |      |

<details open>

<summary> Show Detailed Test Results
</summary>

| file                                  | context | test                                        | status  | n |  time | icon |
| :------------------------------------ | :------ | :------------------------------------------ | :------ | -: | ----: | :--- |
| [test-set.R](testthat/test-set.R#L8)  | set\_   | set functions: set\_template                | PASS    | 1 | 0.011 |      |
| [test-set.R](testthat/test-set.R#L13) | set\_   | set functions: set\_font\_family            | PASS    | 1 | 0.001 |      |
| [test-set.R](testthat/test-set.R#L18) | set\_   | set functions: set\_windows\_control\_theme | PASS    | 1 | 0.001 |      |
| [test-uri.R](testthat/test-uri.R#L9)  | uri     | options: benchmark                          | PASS    | 1 | 0.002 |      |
| [test-uri.R](testthat/test-uri.R#L17) | uri     | uri: benchmark                              | PASS    | 1 | 0.002 |      |
| [test-uri.R](testthat/test-uri.R#L21) | uri     | uri: 200                                    | PASS    | 1 | 0.721 |      |
| [test-uri.R](testthat/test-uri.R#L27) | uri     | encode: encode character                    | PASS    | 1 | 0.002 |      |
| [test-uri.R](testthat/test-uri.R#L31) | uri     | encode: no encode character                 | PASS    | 1 | 0.001 |      |
| [test-uri.R](testthat/test-uri.R#L37) | uri     | tiny: valid tiny                            | PASS    | 1 | 0.620 |      |
| [test-uri.R](testthat/test-uri.R#L42) | uri     | tiny: clipboard                             | SKIPPED | 1 | 0.001 | 🔶    |
| [test-uri.R](testthat/test-uri.R#)    | uri     | bad template: error uri                     | PASS    | 1 | 0.003 |      |
| [test-yml.R](testthat/test-yml.R#L24) | yml     | yaml fields: rgba                           | PASS    | 1 | 0.006 |      |
| [test-yml.R](testthat/test-yml.R#L29) | yml     | yaml fields: template                       | PASS    | 1 | 0.003 |      |
| [test-yml.R](testthat/test-yml.R#L34) | yml     | yaml fields: bad font family                | PASS    | 1 | 0.002 |      |
| [test-yml.R](testthat/test-yml.R#L39) | yml     | yaml fields: pv                             | PASS    | 1 | 0.002 |      |
| [test-yml.R](testthat/test-yml.R#L44) | yml     | yaml fields: ph                             | PASS    | 1 | 0.002 |      |
| [test-yml.R](testthat/test-yml.R#L59) | yml     | namesless palette: fill in palette          | PASS    | 1 | 0.002 |      |

| Failed | Warning | Skipped |
| :----- | :------ | :------ |
| 🛑      | ⚠️      | 🔶       |

</details>

<details>

<summary> Session Info </summary>

| Field    | Value                               |
| :------- | :---------------------------------- |
| Version  | R version 3.5.1 (2018-07-02)        |
| Platform | x86\_64-apple-darwin15.6.0 (64-bit) |
| Running  | macOS 10.14.2                       |
| Language | en\_US                              |
| Timezone | America/Chicago                     |

| Package  | Version    |
| :------- | :--------- |
| testthat | 2.0.0.9000 |
| covr     | 3.2.0      |
| covrpage | 0.0.69     |

</details>

<!--- Final Status : skipped/warning --->
Tests and Coverage
================
17 January, 2019 12:18:36

This output is created by
[covrpage](https://github.com/metrumresearchgroup/covrpage).

## Coverage

Coverage summary is created using the
[covr](https://github.com/r-lib/covr) package.

| Object                                               | Coverage (%) |
| :--------------------------------------------------- | :----------: |
| carbonate                                            |    48.58     |
| [R/carbonate.R](../R/carbonate.R)                    |     0.00     |
| [R/selenium\_functions.R](../R/selenium_functions.R) |     0.00     |
| [R/uri\_functions.R](../R/uri_functions.R)           |    52.63     |
| [R/carbon.R](../R/carbon.R)                          |    64.29     |
| [R/helpers.R](../R/helpers.R)                        |    70.09     |
| [R/set\_get\_functions.R](../R/set_get_functions.R)  |    100.00    |

<br>

## Unit Tests

Unit Test summary is created using the
[testthat](https://github.com/r-lib/testthat)
package.

| file                              | n |  time | error | failed | skipped | warning | icon |
| :-------------------------------- | -: | ----: | ----: | -----: | ------: | ------: | :--- |
| [test-set.R](testthat/test-set.R) | 3 | 0.013 |     0 |      0 |       0 |       0 |      |
| [test-uri.R](testthat/test-uri.R) | 8 | 1.352 |     0 |      0 |       1 |       0 | 🔶    |
| [test-yml.R](testthat/test-yml.R) | 6 | 0.017 |     0 |      0 |       0 |       0 |      |

<details open>

<summary> Show Detailed Test Results
</summary>

| file                                  | context | test                                        | status  | n |  time | icon |
| :------------------------------------ | :------ | :------------------------------------------ | :------ | -: | ----: | :--- |
| [test-set.R](testthat/test-set.R#L8)  | set\_   | set functions: set\_template                | PASS    | 1 | 0.011 |      |
| [test-set.R](testthat/test-set.R#L13) | set\_   | set functions: set\_font\_family            | PASS    | 1 | 0.001 |      |
| [test-set.R](testthat/test-set.R#L18) | set\_   | set functions: set\_windows\_control\_theme | PASS    | 1 | 0.001 |      |
| [test-uri.R](testthat/test-uri.R#L9)  | uri     | options: benchmark                          | PASS    | 1 | 0.002 |      |
| [test-uri.R](testthat/test-uri.R#L17) | uri     | uri: benchmark                              | PASS    | 1 | 0.002 |      |
| [test-uri.R](testthat/test-uri.R#L21) | uri     | uri: 200                                    | PASS    | 1 | 0.721 |      |
| [test-uri.R](testthat/test-uri.R#L27) | uri     | encode: encode character                    | PASS    | 1 | 0.002 |      |
| [test-uri.R](testthat/test-uri.R#L31) | uri     | encode: no encode character                 | PASS    | 1 | 0.001 |      |
| [test-uri.R](testthat/test-uri.R#L37) | uri     | tiny: valid tiny                            | PASS    | 1 | 0.620 |      |
| [test-uri.R](testthat/test-uri.R#L42) | uri     | tiny: clipboard                             | SKIPPED | 1 | 0.001 | 🔶    |
| [test-uri.R](testthat/test-uri.R#)    | uri     | bad template: error uri                     | PASS    | 1 | 0.003 |      |
| [test-yml.R](testthat/test-yml.R#L24) | yml     | yaml fields: rgba                           | PASS    | 1 | 0.006 |      |
| [test-yml.R](testthat/test-yml.R#L29) | yml     | yaml fields: template                       | PASS    | 1 | 0.003 |      |
| [test-yml.R](testthat/test-yml.R#L34) | yml     | yaml fields: bad font family                | PASS    | 1 | 0.002 |      |
| [test-yml.R](testthat/test-yml.R#L39) | yml     | yaml fields: pv                             | PASS    | 1 | 0.002 |      |
| [test-yml.R](testthat/test-yml.R#L44) | yml     | yaml fields: ph                             | PASS    | 1 | 0.002 |      |
| [test-yml.R](testthat/test-yml.R#L59) | yml     | namesless palette: fill in palette          | PASS    | 1 | 0.002 |      |

| Failed | Warning | Skipped |
| :----- | :------ | :------ |
| 🛑      | ⚠️      | 🔶       |

</details>

<details>

<summary> Session Info </summary>

| Field    | Value                               |
| :------- | :---------------------------------- |
| Version  | R version 3.5.1 (2018-07-02)        |
| Platform | x86\_64-apple-darwin15.6.0 (64-bit) |
| Running  | macOS 10.14.2                       |
| Language | en\_US                              |
| Timezone | America/Chicago                     |

| Package  | Version    |
| :------- | :--------- |
| testthat | 2.0.0.9000 |
| covr     | 3.2.0      |
| covrpage | 0.0.69     |

</details>

<!--- Final Status : skipped/warning --->
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  eval = FALSE,
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# rtweet <img src="man/figures/logo.png" width="160px" align="right" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/ropensci/rtweet/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rtweet/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/rtweet)](https://cran.r-project.org/package=rtweet)
[![Coverage Status](https://codecov.io/gh/ropensci/rtweet/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/rtweet?branch=master)
![Downloads](https://cranlogs.r-pkg.org/badges/rtweet)
[![ZENODO](https://zenodo.org/badge/64161359.svg)](https://zenodo.org/badge/latestdoi/64161359)
[![rOpenSci](https://badges.ropensci.org/302_status.svg)](https://github.com/ropensci/software-review/issues/302)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.01829/status.svg)](https://doi.org/10.21105/joss.01829)
<!-- badges: end -->

Use twitter from R. Get started by reading `vignette("rtweet")`.

## Installation

To get the current released version from CRAN:

```{r}
install.packages("rtweet")
```

## Usage

All users must be authenticated to interact with Twitter's APIs. The easiest way to authenticate is to use your personal twitter account - this will happen automatically (via a browser popup) the first time you use an rtweet function. See `auth_setup_default()` for details. Using your personal account is fine for casual use, but if you are trying to collect a lot of data it's a good idea to authentication with your own Twitter "app". See `vignette("auth", package = "rtweet")` for details.

```{r}
library(rtweet)
```

rtweet should be used in strict accordance with Twitter's [developer  terms](https://developer.twitter.com/en/developer-terms/more-on-restricted-use-cases).

### Search tweets or users

Search for up to 10,000 tweets containing #rstats, the common hashtag used to refer to the R language, excluding retweets:

```{r}
rt <- search_tweets("#rstats", n = 10000, include_rts = FALSE)
```

Twitter rate limits cap the number of search results returned to 18,000 every 15 minutes. To request more than that, set `retryonratelimit = TRUE` and rtweet will wait for rate limit
resets for you.

Search for 1,000 users with the #rstats in their profile:

```{r}
usrs <- search_users("#rstats", n = 1000)
```

### Stream tweets

Randomly sample (approximately 1%) from the live stream of all tweets:

```{r}
rt <- stream_tweets("")
```

Stream all geo-located tweets from London for 60 seconds:

```{r}
rt <- stream_tweets(location = lookup_coords("london"), timeout = 60)
```

### Get friends and followers

Get all accounts followed by a user:

```{r}
## get user IDs of accounts followed by R Foundation
R_Foundation_fds <- get_friends("_R_Foundation")

## lookup data on those accounts
R_Foundation_fds_data <- lookup_users(R_Foundation_fds$user_id)
```

Get all accounts following a user:

```{r}
## get user IDs of accounts following R Foundation
R_Foundation_flw <- get_followers("_R_Foundation", n = 10000)
R_Foundation_flw_data <- lookup_users(R_Foundation_flw$user_id)
```

If you want *all* followers, you'll need you'll need to set `n = Inf` and `retryonratelimit = TRUE` but be warned that this might take a *long* time.

### Get timelines

Get the most recent 3,200 tweets from R Foundation:

```{r}
## get user IDs of accounts followed by R Foundation
tmls <- get_timelines("_R_Foundation", n = 3200)
```

### Get favorites

Get the 3,000 most recently favorited statuses by R Foundation:

```{r}
jkr <- get_favorites("_R_Foundation", n = 3000)
```

## Contact

Communicating with Twitter's APIs relies on an internet connection, which can sometimes be inconsistent. With that said, if you encounter an obvious bug for which there is not already an active [issue](https://github.com/ropensci/rtweet/issues), please [create a new issue](https://github.com/ropensci/rtweet/issues/new) with all code used (preferably a reproducible example) on Github.

# Code of Conduct

Please note that this package is released with a [Contributor
Code of Conduct](https://ropensci.org/code-of-conduct/). 
By contributing to this project, you agree to abide by its terms.
---
title: "Intro to rtweet: Collecting Twitter Data"
author: "Michael W. Kearney"
output:
  rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Intro to rtweet}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE, eval = FALSE, comment = "#>", collapse = TRUE)
```

This vignette provides a quick tour of the R package `rtweet: Collecting Twitter Data`.

### Search tweets

Search for up to 18,000 (non-retweeted) tweets containing the rstats hashtag.
```{r, eval=FALSE}
## search for 18000 tweets using the rstats hashtag
rt <- search_tweets(
  "#rstats", n = 18000, include_rts = FALSE
)

## preview tweets data
rt

## preview users data
users_data(rt)

## plot time series (if ggplot2 is installed)
ts_plot(rt)
```


Quickly visualize frequency of tweets over time using `ts_plot()`.
```{r, eval=FALSE}
## plot time series of tweets
ts_plot(rt, "3 hours") +
  ggplot2::theme_minimal() +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")) +
  ggplot2::labs(
    x = NULL, y = NULL,
    title = "Frequency of #rstats Twitter statuses from past 9 days",
    subtitle = "Twitter status (tweet) counts aggregated using three-hour intervals",
    caption = "\nSource: Data collected from Twitter's REST API via rtweet"
  )
```

Twitter rate limits cap the number of search results returned to
18,000 every 15 minutes. To request more than that, simply set
`retryonratelimit = TRUE` and rtweet will wait for rate limit
resets for you.
```{r, eval=FALSE}
## search for 250,000 tweets containing the word data
rt <- search_tweets(
  "data", n = 250000, retryonratelimit = TRUE
)
```

Search by geo-location---for example, find 10,000 tweets in the English
language sent from the United States.
```{r, eval=FALSE}
## search for 10,000 tweets sent from the US
rt <- search_tweets(
  "lang:en", geocode = lookup_coords("usa"), n = 10000
)

## create lat/lng variables using all available tweet and profile geo-location data
rt <- lat_lng(rt)

## plot state boundaries
par(mar = c(0, 0, 0, 0))
maps::map("state", lwd = .25)

## plot lat and lng points onto state map
with(rt, points(lng, lat, pch = 20, cex = .75, col = rgb(0, .3, .7, .75)))
```

### Stream tweets

Randomly sample (approximately 1%) from the live stream of all tweets.
```{r, eval=FALSE}
## random sample for 30 seconds (default)
rt <- stream_tweets("")
```

Stream all geo enabled tweets from London for 60 seconds.
```{r, eval=FALSE}
## stream tweets from london for 60 seconds
rt <- stream_tweets(lookup_coords("london, uk"), timeout = 60)
```

Stream all tweets mentioning #rstats for a week.
```{r, eval=FALSE}
## stream london tweets for a week (60 secs x 60 mins * 24 hours *  7 days)
stream_tweets(
  "#rstats",
  timeout = 60 * 60 * 24 * 7,
  file_name = "tweets_about_R.json",
  parse = FALSE
)

## read in the data as a tidy tbl data frame
djt <- parse_stream("tweets_about_R.json.json")
```

### Get friends

Retrieve a list of all the accounts a **user follows**.
```{r, eval=FALSE}
## get user IDs of accounts followed by R Foundation
R_foundation_fds <- get_friends("_R_Foundation")

## lookup data on those accounts
R_foundation_fds_data <- lookup_users(R_foundation$user_id)
```

### Get followers

Retrieve a list of the **accounts following** a user.
```{r, eval=FALSE}
## get user IDs of accounts following R Foundation
R_foundation_flw <- get_followers("_R_Foundation", n = 75000)

## lookup data on those accounts
R_foundation_flw_data <- lookup_users(R_foundation_flw$user_id)
```

Or if you really want ALL of their followers:
```{r, eval=FALSE}
## how many total follows does R Foundation have?
R_foundation <- lookup_users("_R_Foundation")

## get them all (this would take a little over 5 days)
R_foundation_flw <- get_followers(
  "_R_Foundation", n = R_foundation$followers_count, retryonratelimit = TRUE
)
```

### Get timelines

Get the most recent 3,200 tweets from R Foundation.
```{r, eval=FALSE}
## get user IDs of accounts followed by R Foundation
tmls <- get_timelines("_R_Foundation", n = 3200)

## plot the frequency of tweets for each user over time
tmls %>%
  dplyr::filter(created_at > "2017-10-29") %>%
  dplyr::group_by(screen_name) %>%
  ts_plot("days", trim = 1L) +
  ggplot2::geom_point() +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.title = ggplot2::element_blank(),
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold")) +
  ggplot2::labs(
    x = NULL, y = NULL,
    title = "Frequency of Twitter statuses posted by news organization",
    subtitle = "Twitter status (tweet) counts aggregated by day from October/November 2017",
    caption = "\nSource: Data collected from Twitter's REST API via rtweet"
  )
```

### Get favorites

Get the 3,000 most recently favorited statuses by R Foundation.
```{r, eval=FALSE}
jkr <- get_favorites("_R_Foundation", n = 3000)
```

### Search users

Search for 1,000 users with the rstats hashtag in their profile bios.
```{r, eval=FALSE}
## search for users with #rstats in their profiles
usrs <- search_users("#rstats", n = 1000)
```

### Get trends

Discover what's currently trending in San Francisco.
```{r, eval=FALSE}
sf <- get_trends("san francisco")
```

### Lookup users

```{r, eval=FALSE}
## lookup users by screen_name or user_id
users <- c("KimKardashian", "justinbieber", "taylorswift13",
           "espn", "JoelEmbiid", "cstonehoops", "KUHoops",
           "upshotnyt", "fivethirtyeight", "hadleywickham",
           "cnn", "foxnews", "msnbc", "maddow", "seanhannity",
           "potus", "epa", "hillaryclinton", "realdonaldtrump",
           "natesilver538", "ezraklein", "annecoulter")
famous_tweeters <- lookup_users(users)

## preview users data
famous_tweeters

# extract most recent tweets data from the famous tweeters
tweets_data(famous_tweeters)
```


### Posting statuses

```{r, eval=FALSE}
post_tweet("my first rtweet #rstats")
```

### Following users

```{r, eval=FALSE}
## ty for the follow ;)
post_follow("kearneymw")
```
---
title: "Authentication with rtweet"
output:
  rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Authentication with rtweet}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, eval = FALSE, comment = "#>", collapse = TRUE)
```

rtweet's default authentication mechanism (`auth_setup_default()`) allows you to act on behalf of your personal Twitter account, as if you were performing actions on twitter.com.
This mechanism is very simple to set up, but is best suited for casual use.

If you want to collect a lot of data or implement a bot, you should instead use one of rtweet's two other authentication mechanisms:

-   **App authentication** allows you to act as if you were a Twitter app.
    You can't perform operations that a user can (like posting a tweet or reading a DM), but you get higher rate limits on data collection operations.

-   **Bot authentication** allows you to create a fully automated Twitter bot that performs actions on its own behalf rather than on behalf of a human.

In either case, you'll need to create your own Twitter app, so we'll start by discussing what an app is and how to create one on the Twitter website.
Next, you'll learn how to use the `rtweet_app()` and `rtweet_bot()` functions to tell rtweet about your app config.
You'll then learn how to set the default authentication mechanism for the current R session, and how to save it so you can use it in a future session.

```{r}
library(rtweet)
```

## Creating a Twitter app

You're already familiar with using twitter, either through [the website](https://twitter.com) or an app that you installed on your phone or computer.
To use twitter from R, you'll need to learn a little more about what's going on behind the scenes.
The first important concept to grasp is that every request to the Twitter API has to go through an "app".
Normally, someone else has created the app for you, but now that you're using twitter programmatically, you need to create your own app.
(It's still called an app even though you'll be using it through an R package).

To create a Twitter app, you need to first apply for a developer account by following the instructions at <https://developer.twitter.com>.
Once you have been approved (which may take several hours), navigate to the [developer portal](https://developer.twitter.com/en/portal/projects-and-apps) and click the "Create App" button at the bottom of the page.
You'll need to name your app: the name is unimportant for our purposes, but needs to be unique across all twitter apps.

After you've created your app, you'll see a screen that gives you some important information.
You'll only see this once, so make sure to record it in a secure location.

![](app-info.png){width="548"}

(Don't worry if you forget to save this data: you can always regenerate new values by clicking the "regenerate" button on the "keys and tokens" page.)

## Setup

Now that you have an app, you have to tell rtweet about it.
You'll use either `rtweet_app()` or `rtweet_bot()` depending on whether you want app-style authentication or bot-style authentication as described above.

### rtweet_app()

To use app based authentication, run this code:

```{r, eval = FALSE}
auth <- rtweet_app()
```

This will prompt you to enter the bearer token that you recorded earlier.

It's good practice to only provide secrets interactively, because that makes it harder to accidentally share them in either your `.Rhistory` or an `.R` file.

### `rtweet_bot()`

Bot based authentication works similarly:

```{r, eval = FALSE}
auth <- rtweet_bot()
```

But you'll need more data --- as well as the API key and secret you recorded earlier, you'll also need to generate a "Access token and secret" which you can get by clicking the "Generate" button on the Keys and Tokens page:

![](keys-tokens.png){width="362"}

Again, you'll want to record this data in a secure place because you only get to see it once.

## Default auth

Now you have an auth object that you can provide to the `token` argument of any rtweet function:

```{r, eval = FALSE}
df <- search_tweets("#rstats", token = auth)
```

It's a good idea to do this once to check that you've entered all the app data correctly, but it'd be annoying if you had to pass this object around every single time.
Instead, you can call `auth_as()` to set this as the default for the remainder of the session:

```{r, eval = FALSE}
auth_as(auth)
```

## Saving and loading

`auth_as()` only lasts for a single session; if you close and re-open R, you'll need to repeat the whole process.
This would be annoying (!) so rtweet also provides a way to save and reload auth across sessions:

```{r, eval = FALSE}
auth_save(auth, "some-name")
```

The second argument to `auth_save()` can be any string.
It just needs to be meaningful to you so that you remember exactly what you're loading when you use it a future session:

```{r, eval = FALSE}
auth_as("some-name")
```

You can see all the authentication options you have saved with `auth_list()`.
---
title: "Live streaming tweets"
subtitle: "rtweet: Collecting Twitter Data"
output: 
  rmarkdown::html_vignette:
    self_contained: false
vignette: >
  %\VignetteIndexEntry{Live streaming tweets}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, eval = FALSE, comment = "#>", collapse = TRUE)
```

## Installing and loading package

Prior to streaming, make sure to install and load rtweet. This
vignette assumes users have already setup user access tokens (see: the "auth" vignette).

```{r}
## Install rtweet
install.packages("rtweet")
## Load rtweet
library(rtweet)
```

## Specifying parameters

In addition to accessing Twitter's REST API (e.g., `search_tweets`,
`get_timeline`), rtweet makes it possible to capture live streams of
Twitter data using the `stream_tweets()` function. By default,
`stream_tweets` will stream for 30 seconds and return a random sample
of tweets. To modify the default settings, `stream_tweets` accepts
several parameters, including `q` (query used to filter tweets),
`timeout` (duration or time of stream), and `file_name` (path name for
saving raw json data).

```{r}
## Stream keywords used to filter tweets
q <- "hillaryclinton,imwithher,realdonaldtrump,maga,electionday"

## Stream time in seconds so for one minute set timeout = 60
## For larger chunks of time, I recommend multiplying 60 by the number
## of desired minutes. This method scales up to hours as well
## (x * 60 = x mins, x * 60 * 60 = x hours)
## Stream for 30 minutes
streamtime <- 30 * 60

## Filename to save json data (backup)
filename <- "rtelect.json"
```

## stream_tweets()

Once these parameters are specified, initiate the stream. Note:
Barring any disconnection or disruption of the API, streaming will
occupy your current instance of R until the specified time has
elapsed. It is possible to start a new instance or R---streaming
itself usually isn't very memory intensive---but operations may drag a
bit during the parsing process which takes place immediately after
streaming ends.

```{r}
## Stream election tweets
rt <- stream_tweets(q = q, timeout = streamtime, file_name = filename)
```

Parsing larger streams can take quite a bit of time (in addition to
time spent streaming) due to a somewhat time-consuming simplifying
process used to convert a json file into an R object. For example, the
stream above yielded a little over 140,000 tweets and took my Macbook
Air, which has 4gb of RAM, about 10 minutes to process.

## Saving files

Given a lengthy parsing process, users may want to stream tweets into
json files upfront and parse those files later on. To do this, simply
add `parse = FALSE` and make sure you provide a path (file name) to a
location you can find later. To ensure the stream automatically reconnects following
any interruption prior to the specified stream time, use
`stream_tweets2()`.

Regardless of whether you decide to setup an organizational system for
streaming data, the process of streaming a file to disk and parsing it
at a later point in space-time is the same, as illustrated in the example
below.

```{r}
## No upfront-parse save as json file instead method
stream_tweets(
  q = q,
  parse = FALSE,
  timeout = streamtime,
  file_name = filename
)
## Parse from json file
rt <- parse_stream(filename)

## stream_tweets2 method
twoweeks <- 60L * 60L * 24L * 7L * 2L
congress <- "congress,senate,house of representatives,representatives,senators,legislative"
stream_tweets2(
  q = congress,
  parse = FALSE,
  timeout = twoweeks,
  dir = "congress-stream"
)

## Parse from json file
rt <- parse_stream("congress-stream.json")
```

## Returned data object

The parsed object should be the same whether a user parses up-front or
from a json file in a later session. The returned object should be a
data frame consisting of tweets data.

```{r}
## Preview tweets data
rt
```

The returned object should also include a data frame of users data,
which Twitter's stream API automatically returns along with tweets
data. To access users data, use the `users_data` function.

```{r}
## Preview users data
users_data(rt)
```

## Plotting

Once parsed, `ts_plot()` provides a quick visual of the frequency of
tweets. By default, `ts_plot()` will try to aggregate time by the
day. Because I know the stream only lasted 30 minutes, I've opted to
aggregate tweets by the second. It'd also be possible to aggregate by
the minute, i.e., `by = "mins"`, or by some value of seconds, e.g.,`by
= "15 secs"`. I usually fiddle around with this a bit until the plot
looks good.

```{r}
## Plot time series of all tweets aggregated by second
ts_plot(rt, by = "secs")
```

<p align="center">
<img src="files/stream-ts.png" alt="stream-ts">
</p>

## Plotting with filters

The `ts_plot()` function can also generate multiple time series for
grouped data frames.

```{r}
## plot multiple time series by first grouping the data by screen name
rt %>%
  dplyr::group_by(screen_name) %>%
  ts_plot() +
  ggplot2::labs(
    title = "Tweets during election day for the 2016 U.S. election",
    subtitle = "Tweets collected, parsed, and plotted using `rtweet`"
  )
```

Often times these plots kinda resemble a frowny face with the
first and last points appearing significantly lower than the
rest. This is because the first and last intervals of time are
artificially shrunken by connection and disconnection processes. To
remedy this, users can specify `trim = 1` to tell R to drop the
first and last observation for each time series. This usually yields a
much more attractive looking plot.

<p align="center">
<img src="files/stream-filter.png" alt="stream-filter">
</p>
---
title: "Intro to rtweet: Collecting Twitter Data"
author: "Michael W. Kearney"
output:
  rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Intro to rtweet}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE, eval = FALSE, comment = "#>", collapse = TRUE)
```

This vignette provides a quick tour of the R package `rtweet: Collecting Twitter Data`.

### Search tweets

Search for up to 18,000 (non-retweeted) tweets containing the rstats hashtag.
```{r, eval=FALSE}
## search for 18000 tweets using the rstats hashtag
rt <- search_tweets(
  "#rstats", n = 18000, include_rts = FALSE
)

## preview tweets data
rt

## preview users data
users_data(rt)

## plot time series (if ggplot2 is installed)
ts_plot(rt)
```


Quickly visualize frequency of tweets over time using `ts_plot()`.
```{r, eval=FALSE}
## plot time series of tweets
ts_plot(rt, "3 hours") +
  ggplot2::theme_minimal() +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")) +
  ggplot2::labs(
    x = NULL, y = NULL,
    title = "Frequency of #rstats Twitter statuses from past 9 days",
    subtitle = "Twitter status (tweet) counts aggregated using three-hour intervals",
    caption = "\nSource: Data collected from Twitter's REST API via rtweet"
  )
```

Twitter rate limits cap the number of search results returned to
18,000 every 15 minutes. To request more than that, simply set
`retryonratelimit = TRUE` and rtweet will wait for rate limit
resets for you.
```{r, eval=FALSE}
## search for 250,000 tweets containing the word data
rt <- search_tweets(
  "data", n = 250000, retryonratelimit = TRUE
)
```

Search by geo-location---for example, find 10,000 tweets in the English
language sent from the United States.
```{r, eval=FALSE}
## search for 10,000 tweets sent from the US
rt <- search_tweets(
  "lang:en", geocode = lookup_coords("usa"), n = 10000
)

## create lat/lng variables using all available tweet and profile geo-location data
rt <- lat_lng(rt)

## plot state boundaries
par(mar = c(0, 0, 0, 0))
maps::map("state", lwd = .25)

## plot lat and lng points onto state map
with(rt, points(lng, lat, pch = 20, cex = .75, col = rgb(0, .3, .7, .75)))
```

### Stream tweets

Randomly sample (approximately 1%) from the live stream of all tweets.
```{r, eval=FALSE}
## random sample for 30 seconds (default)
rt <- stream_tweets("")
```

Stream all geo enabled tweets from London for 60 seconds.
```{r, eval=FALSE}
## stream tweets from london for 60 seconds
rt <- stream_tweets(lookup_coords("london, uk"), timeout = 60)
```

Stream all tweets mentioning realDonaldTrump or Trump for a week.
```{r, eval=FALSE}
## stream london tweets for a week (60 secs x 60 mins * 24 hours *  7 days)
stream_tweets(
  "realdonaldtrump,trump",
  timeout = 60 * 60 * 24 * 7,
  file_name = "tweetsabouttrump.json",
  parse = FALSE
)

## read in the data as a tidy tbl data frame
djt <- parse_stream("tweetsabouttrump.json")
```

### Get friends

Retrieve a list of all the accounts a **user follows**.
```{r, eval=FALSE}
## get user IDs of accounts followed by CNN
cnn_fds <- get_friends("cnn")

## lookup data on those accounts
cnn_fds_data <- lookup_users(cnn_fds$user_id)
```

### Get followers

Retrieve a list of the **accounts following** a user.
```{r, eval=FALSE}
## get user IDs of accounts following CNN
cnn_flw <- get_followers("cnn", n = 75000)

## lookup data on those accounts
cnn_flw_data <- lookup_users(cnn_flw$user_id)
```

Or if you really want ALL of their followers:
```{r, eval=FALSE}
## how many total follows does cnn have?
cnn <- lookup_users("cnn")

## get them all (this would take a little over 5 days)
cnn_flw <- get_followers(
  "cnn", n = cnn$followers_count, retryonratelimit = TRUE
)
```

### Get timelines

Get the most recent 3,200 tweets from cnn, BBCWorld, and foxnews.
```{r, eval=FALSE}
## get user IDs of accounts followed by CNN
tmls <- get_timelines(c("cnn", "BBCWorld", "foxnews"), n = 3200)

## plot the frequency of tweets for each user over time
tmls %>%
  dplyr::filter(created_at > "2017-10-29") %>%
  dplyr::group_by(screen_name) %>%
  ts_plot("days", trim = 1L) +
  ggplot2::geom_point() +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.title = ggplot2::element_blank(),
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold")) +
  ggplot2::labs(
    x = NULL, y = NULL,
    title = "Frequency of Twitter statuses posted by news organization",
    subtitle = "Twitter status (tweet) counts aggregated by day from October/November 2017",
    caption = "\nSource: Data collected from Twitter's REST API via rtweet"
  )
```

### Get favorites

Get the 3,000 most recently favorited statuses by JK Rowling.
```{r, eval=FALSE}
jkr <- get_favorites("jk_rowling", n = 3000)
```

### Search users

Search for 1,000 users with the rstats hashtag in their profile bios.
```{r, eval=FALSE}
## search for users with #rstats in their profiles
usrs <- search_users("#rstats", n = 1000)
```

### Get trends

Discover what's currently trending in San Francisco.
```{r, eval=FALSE}
sf <- get_trends("san francisco")
```

### Lookup users

```{r, eval=FALSE}
## lookup users by screen_name or user_id
users <- c("KimKardashian", "justinbieber", "taylorswift13",
           "espn", "JoelEmbiid", "cstonehoops", "KUHoops",
           "upshotnyt", "fivethirtyeight", "hadleywickham",
           "cnn", "foxnews", "msnbc", "maddow", "seanhannity",
           "potus", "epa", "hillaryclinton", "realdonaldtrump",
           "natesilver538", "ezraklein", "annecoulter")
famous_tweeters <- lookup_users(users)

## preview users data
famous_tweeters

# extract most recent tweets data from the famous tweeters
tweets_data(famous_tweeters)
```


### Posting statuses

```{r, eval=FALSE}
post_tweet("my first rtweet #rstats")
```

### Following users

```{r, eval=FALSE}
## ty for the follow ;)
post_follow("kearneymw")
```
---
title: "Obtaining and using access tokens"
subtitle: "rtweet: Collecting Twitter Data"
output:
  rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Obtaining and using access tokens}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, eval = FALSE, comment = "#>", collapse = TRUE)
```

## rtweet

This vignette covers how to obtain and use Twitter API access tokens for use in the `rtweet` package.

## Creating a Twitter App

- To create a Twitter app, navigate to [apps.twitter.com](https://apps.twitter.com/) and create a new app by providing a `Name`, `Description`, and `Website` of your choosing (example screenshot provided below).

- **Important** In the `Callback URL` field, make sure to enter the following: `http://127.0.0.1:1410`

- Check yes if you agree and then click "Create your Twitter application".

<p align="center">
<img src="files/creating.png" alt="creating">
</p>

## Authorization methods

Users can create their personal access token in two different ways. Each method is outlined below.

### 1. Browser-based authentication

- Authentication via web browser requires the `httpuv` package to be installed.

```{r}
## install httpuv if not already
if (!requireNamespace("httpuv", quietly = TRUE)) {
  install.packages("httpuv")
}
```

- Click the tab labeled `Keys and Access Tokens` to retrieve your keys.

<p align="center">
<img src="files/created.png" alt="created">
</p>

- In the `Keys and Access Tokens` tab, locate the values `Consumer Key` (aka "API Key") and `Consumer Secret` (aka "API Secret").

<p align="center">
<img src="files/keys.png" alt="keys">
</p>

- Copy and paste the two keys (along with the name of your app) into an R script file and pass them along to `create_token()`.

```{r}
## autheticate via web browser
token <- create_token(
  app = "rtweet_token",
  consumer_key = "XYznzPFOFZR2a39FwWKN1Jp41",
  consumer_secret = "CtkGEWmSevZqJuKl6HHrBxbCybxI1xGLqrD5ynPd9jG0SoHZbD")
```

- A browser window should pop up. Click to approve (must be signed into twitter.com) and return to R.

- The `create_token()` function should automatically save your token as an environment variable for you. To make sure it worked, compare the created token object to the object returned by `get_token()`

```{r}
## check to see if the token is loaded
identical(twitter_token, get_token())
```


### 2. Access token/secret method

- Click the tab labeled `Keys and Access Tokens` to retrieve your keys.

<p align="center">
<img src="files/created.png" alt="created">
</p>

- In the `Keys and Access Tokens` tab, locate and copy/paste values `Consumer Key` (aka "API Key") and `Consumer Secret` (aka "API Secret") into an R script.

<p align="center">
<img src="files/keys.png" alt="keys">
</p>

- In the `Keys and Access Tokens` tab, scroll down to `Token Actions` and click `Create my access token`.

<p align="center">
<img src="files/gen_token.png" alt="gen_token">
</p>

- That should generate two access keys `Access Token` and `Access Token Secret`

<p align="center">
<img src="files/accesskeys.png" alt="acesskeys">
</p>

- Locate and copy/paste the `Consumer Key` (aka "API Key"), `Consumer Secret` (aka "API Secret"), `Access Token`, and `Access Token Secret` values and pass them along to `create_token()`, storing the output as a `token` object.

```{r}
## authenticate via access token
token <- create_token(
  app = "my_twitter_research_app",
  consumer_key = "XYznzPFOFZR2a39FwWKN1Jp41",
  consumer_secret = "CtkGEWmSevZqJuKl6HHrBxbCybxI1xGLqrD5ynPd9jG0SoHZbD",
  acess_token = "9551451262-wK2EmA942kxZYIwa5LMKZoQA4Xc2uyIiEwu2YXL",
  access_secret = "9vpiSGKg1fIPQtxc5d5ESiFlZQpfbknEN1f1m2xe5byw7")
```

- The `create_token()` function should automatically save your token as an environment variable for you. To make sure it worked, compare the created token object to the object returned by `get_token()`

```{r}
## check to see if the token is loaded
identical(twitter_token, get_token())
```


That's it!
---
title: "FAQ"
subtitle: "rtweet: Collecting Twitter Data"
output:
  rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{FAQ}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, eval = FALSE, comment = "#>", collapse = TRUE)
```

## `Error in init_oauth1.0(...`

#### Context

Occurs when attempting to create a token

```{r}
## these are fake keys
#> create_token(
#>   app = "rtweet_token",
#>   consumer_key = "XYznzPFOFZR2a39FwWKN1Jp41",
#>   consumer_secret = "CtkGEWmSevZqJuKl6HHrBxbCybxI1xGLqrD5ynPd9jG0SoHZbD")
# `Error in init_oauth1.0(endpoint, app, permission = params$permission) :
#  client error: (401) Unauthorized`
```

#### Solutions

1. Make sure you have **at least** rtweet version `0.6.6`
1. Check callback URL
1. Make sure `Callback URL` option in "Settings" tab at https://apps.twitter.com **match exactly** the following: `http://127.0.0.1:1410`
1. Make sure API keys **match exactly** the values for your Twitter app (found under the `Keys and Access Tokens` tab at https://apps.twitter.com)
1. In your app page at https://apps.twitter.com/ under the `Keys and Access Tokens` tab, click `Regenerate Consumer Key and Secret`. Create a new token using the new key and secret.
1. Update the **{httr}** package
1. Update R


## `Error in oauth_listener(...`

#### Context

Occurs when attempting to create a token (using in-browser authorization method)

```{r}
## these are fake keys
#> create_token(
#>   app = "rtweet_token",
#>   consumer_key = "XYznzPFOFZR2a39FwWKN1Jp41",
#>   consumer_secret = "CtkGEWmSevZqJuKl6HHrBxbCybxI1xGLqrD5ynPd9jG0SoHZbD")
#Error in oauth_listener(authorize_url, is_interactive) :
#  httpuv package required to capture OAuth credentials.
```

#### Solutions
1. Install the **{httpuv}** package

```{r}
install.packages("httpuv")
```

## `Warning: 89 - Invalid or expired token.`

#### Context

Occurs when sending request to Twitter API. It means you have invalid or expired keys stored in access token.

```{r}
#> search_tweets("lang:en")
#Warning: 89 - Invalid or expired token.
```

#### Solutions
1. Create new token using the keys from your previously created Twitter application found at https://apps.twitter.com
1. In your app page at https://apps.twitter.com/ under the `Keys and Access Tokens` tab, click `Regenerate Consumer Key and Secret`. Create a new token using the new key and secret.

```{r}
## these are fake keys
#> create_token(
#>   app = "rtweet_token",
#>   consumer_key = "XYznzPFOFZR2a39FwWKN1Jp41",
#>   consumer_secret = "CtkGEWmSevZqJuKl6HHrBxbCybxI1xGLqrD5ynPd9jG0SoHZbD")
```
---
title: "Live streaming tweets"
subtitle: "rtweet: Collecting Twitter Data"
output:
  rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Live streaming tweets}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE, eval = FALSE, comment = "#>", collapse = TRUE)
```

## Installing and loading package

Prior to streaming, make sure to install and load rtweet. This
vignette assumes users have already setup user access tokens (see: [obtaining and using access tokens](https://mkearney.github.io/rtweet/articles/auth.html)).

```{r}
## Install rtweet
install.packages("rtweet")
## Load rtweet
library(rtweet)
```

## Specifying parameters

In addition to accessing Twitter's REST API (e.g., `search_tweets`,
`get_timeline`), rtweet makes it possible to capture live streams of
Twitter data using the `stream_tweets()` function. By default,
`stream_tweets` will stream for 30 seconds and return a random sample
of tweets. To modify the default settings, `stream_tweets` accepts
several parameters, including `q` (query used to filter tweets),
`timeout` (duration or time of stream), and `file_name` (path name for
saving raw json data).

```{r}
## Stream keywords used to filter tweets
q <- "hillaryclinton,imwithher,realdonaldtrump,maga,electionday"

## Stream time in seconds so for one minute set timeout = 60
## For larger chunks of time, I recommend multiplying 60 by the number
## of desired minutes. This method scales up to hours as well
## (x * 60 = x mins, x * 60 * 60 = x hours)
## Stream for 30 minutes
streamtime <- 30 * 60

## Filename to save json data (backup)
filename <- "rtelect.json"
```

## stream_tweets()

Once these parameters are specified, initiate the stream. Note:
Barring any disconnection or disruption of the API, streaming will
occupy your current instance of R until the specified time has
elapsed. It is possible to start a new instance or R---streaming
itself usually isn't very memory intensive---but operations may drag a
bit during the parsing process which takes place immediately after
streaming ends.

```{r}
## Stream election tweets
rt <- stream_tweets(q = q, timeout = streamtime, file_name = filename)
```

Parsing larger streams can take quite a bit of time (in addition to
time spent streaming) due to a somewhat time-consuming simplifying
process used to convert a json file into an R object. For example, the
stream above yielded a little over 140,000 tweets and took my Macbook
Air, which has 4gb of RAM, about 10 minutes to process.

## Saving files

Given a lengthy parsing process, users may want to stream tweets into
json files upfront and parse those files later on. To do this, simply
add `parse = FALSE` and make sure you provide a path (file name) to a
location you can find later. To ensure the stream automatically reconnects following
any interruption prior to the specified stream time, use
`stream_tweets2()`.

Regardless of whether you decide to setup an organizational system for
streaming data, the process of streaming a file to disk and parsing it
at a later point in space-time is the same, as illustrated in the example
below.

```{r}
## No upfront-parse save as json file instead method
stream_tweets(
  q = q,
  parse = FALSE,
  timeout = streamtime,
  file_name = filename
)
## Parse from json file
rt <- parse_stream(filename)

## stream_tweets2 method
twoweeks <- 60L * 60L * 24L * 7L * 2L
congress <- "congress,senate,house of representatives,representatives,senators,legislative"
stream_tweets2(
  q = congress,
  parse = FALSE,
  timeout = twoweeks,
  dir = "congress-stream"
)

## Parse from json file
rt <- parse_stream("congress-stream.json")
```

## Returned data object

The parsed object should be the same whether a user parses up-front or
from a json file in a later session. The returned object should be a
data frame consisting of tweets data.

```{r}
## Preview tweets data
rt
```

The returned object should also include a data frame of users data,
which Twitter's stream API automatically returns along with tweets
data. To access users data, use the `users_data` function.

```{r}
## Preview users data
users_data(rt)
```

## Plotting

Once parsed, `ts_plot()` provides a quick visual of the frequency of
tweets. By default, `ts_plot()` will try to aggregate time by the
day. Because I know the stream only lasted 30 minutes, I've opted to
aggregate tweets by the second. It'd also be possible to aggregate by
the minute, i.e., `by = "mins"`, or by some value of seconds, e.g.,`by
= "15 secs"`. I usually fiddle around with this a bit until the plot
looks good.

```{r}
## Plot time series of all tweets aggregated by second
ts_plot(rt, by = "secs")
```

<p align="center">
<img src="files/stream-ts.png" alt="stream-ts">
</p>

## Plotting with filters

The `ts_plot()` function can also generate multiple time series for
grouped data frames.

```{r}
## plot multiple time series by first grouping the data by screen name
rt %>%
  dplyr::group_by(screen_name) %>%
  ts_plot() +
  ggplot2::labs(
    title = "Tweets during election day for the 2016 U.S. election",
    subtitle = "Tweets collected, parsed, and plotted using `rtweet`"
  )
```

Often times these plots kinda resemble a frowny face with the
first and last points appearing significantly lower than the
rest. This is because the first and last intervals of time are
artificially shrunken by connection and disconnection processes. To
remedy this, users can specify `trim = 1` to tell R to drop the
first and last observation for each time series. This usually yields a
much more attractive looking plot.

<p align="center">
<img src="files/stream-filter.png" alt="stream-filter">
</p>
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trends.R
\name{get_trends}
\alias{get_trends}
\title{Get Twitter trends data.}
\usage{
get_trends(
  woeid = 1,
  lat = NULL,
  lng = NULL,
  exclude_hashtags = FALSE,
  token = NULL,
  parse = TRUE
)
}
\arguments{
\item{woeid}{Numeric, WOEID (Yahoo! Where On Earth ID) or character
string of desired town or country. Users may also supply latitude
and longitude coordinates to fetch the closest available trends
data given the provided location. Latitude/longitude coordinates
should be provided as WOEID value consisting of 2 numeric values
or via one latitude value and one longitude value (to the
appropriately named parameters).  To browse all available trend
places, see \code{\link[=trends_available]{trends_available()}}}

\item{lat}{Optional alternative to WOEID. Numeric, latitude in
degrees.  If two coordinates are provided for WOEID, this
function will coerce the first value to latitude.}

\item{lng}{Optional alternative to WOEID. Numeric, longitude in
degrees.  If two coordinates are provided for WOEID, this
function will coerce the second value to longitude.}

\item{exclude_hashtags}{Logical, indicating whether or not to
exclude hashtags. Defaults to FALSE--meaning, hashtags are
included in returned trends.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}
}
\value{
Tibble data frame of trends data for a given geographical area.
}
\description{
Get Twitter trends data.
}
\examples{

\dontrun{

## Retrieve available trends
trends <- trends_available()
trends

## Store WOEID for Worldwide trends
worldwide <- trends$woeid[grep("world", trends$name, ignore.case = TRUE)[1]]

## Retrieve worldwide trends datadata
ww_trends <- get_trends(worldwide)

## Preview trends data
ww_trends

## Retrieve trends data using latitude, longitude near New York City
nyc_trends <- get_trends_closest(lat = 40.7, lng = -74.0)

## should be same result if lat/long supplied as first argument
nyc_trends <- get_trends_closest(c(40.7, -74.0))

## Preview trends data
nyc_trends

## Provide a city or location name using a regular expression string to
## have the function internals do the WOEID lookup/matching for you
(luk <- get_trends("london"))

}

}
\seealso{
Other trends: 
\code{\link{trends_available}()}
}
\concept{trends}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lists_subscriptions.R
\name{lists_subscriptions}
\alias{lists_subscriptions}
\title{Get list subscriptions of a given user.}
\usage{
lists_subscriptions(
  user,
  n = 20,
  cursor = "-1",
  parse = TRUE,
  retryonratelimit = NULL,
  verbose = TRUE,
  token = NULL
)
}
\arguments{
\item{user}{Character vector of screen names or user ids.
See \code{\link[=as_screenname]{as_screenname()}} for more details.}

\item{n}{Desired number of results to return. Results are downloaded
in pages when \code{n} is large; the default value will download a single
page. Set \code{n = Inf} to download as many results as possible.

The Twitter API rate limits the number of requests you can perform
in each 15 minute period. The easiest way to download more than that is
to use \code{retryonratelimit = TRUE}.

You are not guaranteed to get exactly \code{n} results back. You will get
fewer results when tweets have been deleted or if you hit a rate limit.
You will get more results if you ask for a number of tweets that's not
a multiple of page size, e.g. if you request \code{n = 150} and the page
size is 200, you'll get 200 results back.}

\item{cursor}{Which page of results to return. The default will return
the first page; you can supply the result from a previous call to
continue pagination from where it left off.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}

\item{retryonratelimit}{If \code{TRUE}, and a rate limit is exhausted, will wait
until it refreshes. Most Twitter rate limits refresh every 15 minutes.
If \code{FALSE}, and the rate limit is exceeded, the function will terminate
early with a warning; you'll still get back all results received up to
that point. The default value, \code{NULL}, consults the option
\code{rtweet.retryonratelimit} so that you can globally set it to \code{TRUE},
if desired.

If you expect a query to take hours or days to perform, you should not
rely soley on \code{retryonratelimit} because it does not handle other common
failure modes like temporarily losing your internet connection.}

\item{verbose}{Show progress bars and other messages indicating current
progress?}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}
}
\description{
Get list subscriptions of a given user.
}
\examples{

\dontrun{

## get kearneymw subscriptions
rstats <- lists_subscriptions(
  user = "kearneymw",
  n = 1000
)

}

}
\references{
\url{https://developer.twitter.com/en/docs/twitter-api/v1/accounts-and-users/create-manage-lists/api-reference/get-lists-subscriptions}
}
\seealso{
Other lists: 
\code{\link{lists_members}()},
\code{\link{lists_statuses}()},
\code{\link{lists_subscribers}()},
\code{\link{lists_users}()}
}
\concept{lists}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/next_cursor.R
\name{previous_cursor}
\alias{previous_cursor}
\title{Previous cursor}
\usage{
previous_cursor(x)
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
Reverse pagination is no longer supported.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post-tweet.R
\name{post_tweet}
\alias{post_tweet}
\alias{post_status}
\title{Posts status update to user's Twitter account}
\usage{
post_tweet(
  status = "my first rtweet #rstats",
  media = NULL,
  token = NULL,
  in_reply_to_status_id = NULL,
  destroy_id = NULL,
  retweet_id = NULL,
  auto_populate_reply_metadata = FALSE,
  media_alt_text = NULL,
  lat = NULL,
  long = NULL,
  display_coordinates = FALSE
)
}
\arguments{
\item{status}{Character, tweet status. Must be 280 characters or less.}

\item{media}{Length 1 character vector with a file path to video media \strong{OR}
up-to length 4 character vector with file paths to static images to be included in tweet.
\strong{The caller is responsible for managing this.}}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}

\item{in_reply_to_status_id}{Status ID of tweet to which you'd like to reply.
Note: in line with the Twitter API, this parameter is ignored unless the
author of the tweet this parameter references is mentioned within the
status text.}

\item{destroy_id}{To delete a status, supply the single status ID here. If a
character string is supplied, overriding the default (NULL), then a destroy
request is made (and the status text and media attachments) are irrelevant.}

\item{retweet_id}{To retweet a status, supply the single status ID here. If a
character string is supplied, overriding the default (NULL), then a retweet
request is made (and the status text and media attachments) are irrelevant.}

\item{auto_populate_reply_metadata}{If set to TRUE and used with
in_reply_to_status_id, leading @mentions will be looked up from the
original Tweet, and added to the new Tweet from there. Defaults to FALSE.}

\item{media_alt_text}{attach additional \href{https://en.wikipedia.org/wiki/Alt_attribute}{alt text}
metadata to the \code{media} you are uploading. Should be same length as
\code{media} (i.e. as many alt text entries as there are \code{media} entries). See
\href{https://developer.twitter.com/en/docs/media/upload-media/api-reference/post-media-metadata-create}{the official API documentation}
for more information.}

\item{lat}{A numeric value representing the latitude of the location the
tweet refers to. Range should be between -90 and 90 (north). Note that you
should enable the "Precise location" option in your account via \emph{Settings
and privacy > Privacy and Safety > Location}. See
\href{https://help.twitter.com/en/safety-and-security/twitter-location-services-for-mobile}{the official Help Center section}.}

\item{long}{A numeric value representing the longitude of the location the
tweet refers to. Range should be between -180 and 180 (west). See
\code{lat} parameter.}

\item{display_coordinates}{Put a pin on the exact coordinates a tweet has
been sent from. Value should be TRUE or FALSE. This parameter would apply
only if you have provided a valid \code{lat/long} pair of valid values.}
}
\description{
Posts status update to user's Twitter account
}
\examples{
\dontrun{
## generate data to make/save plot (as a .png file)
x <- rnorm(300)
y <- x + rnorm(300, 0, .75)
col <- c(rep("#002244aa", 50), rep("#440000aa", 50))
bg <- c(rep("#6699ffaa", 50), rep("#dd6666aa", 50))

## create temporary file name
tmp <- tempfile(fileext = ".png")

## save as png
png(tmp, 6, 6, "in", res = 127.5)
par(tcl = -.15, family = "Inconsolata",
    font.main = 2, bty = "n", xaxt = "l", yaxt = "l",
    bg = "#f0f0f0", mar = c(3, 3, 2, 1.5))
plot(x, y, xlab = NULL, ylab = NULL, pch = 21, cex = 1,
     bg = bg, col = col,
     main = "This image was uploaded by rtweet")
grid(8, lwd = .15, lty = 2, col = "#00000088")
dev.off()

## post tweet with media attachment
post_tweet("a tweet with media attachment", media = tmp)

# example of replying within a thread
## first post
post_tweet(status="first in a thread")

## lookup status_id
my_timeline <- get_my_timeline()

## ID for reply
reply_id <- my_timeline$status_id[1]

## post reply
post_tweet("second in the thread",
  in_reply_to_status_id = reply_id)
}
}
\references{
Tweet: \url{https://developer.twitter.com/en/docs/twitter-api/v1/tweets/post-and-engage/api-reference/post-statuses-update}
Retweet: \url{https://developer.twitter.com/en/docs/twitter-api/v1/tweets/post-and-engage/api-reference/post-statuses-retweet-id}
Media: \url{https://developer.twitter.com/en/docs/twitter-api/v1/media/upload-media/api-reference/post-media-upload}
Alt-text: \url{https://developer.twitter.com/en/docs/twitter-api/v1/media/upload-media/api-reference/post-media-metadata-create}
}
\seealso{
Other post: 
\code{\link{post_favorite}()},
\code{\link{post_follow}()},
\code{\link{post_friendship}()}
}
\concept{post}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mentions.R
\name{get_mentions}
\alias{get_mentions}
\title{Get mentions for the authenticating user.}
\usage{
get_mentions(
  n = 200,
  since_id = NULL,
  max_id = NULL,
  parse = TRUE,
  retryonratelimit = NULL,
  verbose = TRUE,
  token = NULL,
  ...
)
}
\arguments{
\item{n}{Desired number of results to return. Results are downloaded
in pages when \code{n} is large; the default value will download a single
page. Set \code{n = Inf} to download as many results as possible.

The Twitter API rate limits the number of requests you can perform
in each 15 minute period. The easiest way to download more than that is
to use \code{retryonratelimit = TRUE}.

You are not guaranteed to get exactly \code{n} results back. You will get
fewer results when tweets have been deleted or if you hit a rate limit.
You will get more results if you ask for a number of tweets that's not
a multiple of page size, e.g. if you request \code{n = 150} and the page
size is 200, you'll get 200 results back.}

\item{since_id}{Supply a vector of ids or a data frame of previous results to
find tweets \strong{newer} than \code{since_id}.}

\item{max_id}{Supply a vector of ids or a data frame of previous results to
find tweets \strong{older} than \code{max_id}.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}

\item{retryonratelimit}{If \code{TRUE}, and a rate limit is exhausted, will wait
until it refreshes. Most Twitter rate limits refresh every 15 minutes.
If \code{FALSE}, and the rate limit is exceeded, the function will terminate
early with a warning; you'll still get back all results received up to
that point. The default value, \code{NULL}, consults the option
\code{rtweet.retryonratelimit} so that you can globally set it to \code{TRUE},
if desired.

If you expect a query to take hours or days to perform, you should not
rely soley on \code{retryonratelimit} because it does not handle other common
failure modes like temporarily losing your internet connection.}

\item{verbose}{Show progress bars and other messages indicating current
progress?}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}

\item{...}{Other arguments passed as parameters in composed API
query.}
}
\value{
Tibble of mentions data.
}
\description{
Returns data on up to 200 of the most recent mentions (Tweets
containing a users's screen_name) of the authenticating user.
The timeline returned is the equivalent of the one seen when you view
your mentions on twitter.com.
}
\examples{

\dontrun{
tw <- get_mentions()
tw

# newer mentions
get_mentions(since_id = tw)
}
}
\references{
\url{https://developer.twitter.com/en/docs/tweets/timelines/api-reference/get-statuses-mentions_timeline}
}
\seealso{
Other tweets: 
\code{\link{get_favorites}()},
\code{\link{get_timeline}()},
\code{\link{lists_statuses}()},
\code{\link{lookup_tweets}()},
\code{\link{search_tweets}()}
}
\concept{tweets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/friends.R
\name{lookup_friendships}
\alias{lookup_friendships}
\title{Lookup friendship information between two specified users.}
\usage{
lookup_friendships(source, target, parse = TRUE, token = NULL)
}
\arguments{
\item{source}{Screen name or user id of source user.}

\item{target}{Screen name or user id of target user.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}
}
\description{
Gets information on friendship between two Twitter users.
}
\references{
\url{https://developer.twitter.com/en/docs/twitter-api/v1/accounts-and-users/follow-search-get-users/api-reference/get-friendships-show}
}
\seealso{
Other friends: 
\code{\link{my_friendships}()}
}
\concept{friends}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{langs}
\alias{langs}
\title{Language codes recognized by Twitter data.}
\format{
A tibble with five variables and 486 observations.
}
\usage{
langs
}
\description{
This data comes from the Library of Congress,
\url{http://www.loc.gov/standards/iso639-2/ISO-639-2_utf-8.txt}. The data are
descriptions and codes associated with internationally recognized languages.
Variables include translations for each language represented as
bibliographic, terminologic, alpha, english, and french.
}
\examples{
head(langs)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtweet-package.R
\docType{package}
\name{rtweet-package}
\alias{rtweet}
\alias{rtweet-package}
\title{rtweet: Collect Twitter data from R}
\description{
rtweet provides users a range of functions designed to extract data
from Twitter's REST and streaming APIs. It has three main goals:
\itemize{
\item Formulate and send requests to Twitter's REST and stream APIs.
\item Retrieve and iterate over returned data.
\item Wrangling data into tidy structures.
}

Get started by reading \code{vignette("rtweet")}.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/rtweet}
  \item \url{https://github.com/ropensci/rtweet}
  \item Report bugs at \url{https://github.com/ropensci/rtweet/issues}
}

}
\author{
\strong{Maintainer}: Michael W. Kearney \email{kearneymw@missouri.edu} (\href{https://orcid.org/0000-0002-0730-4694}{ORCID})

Authors:
\itemize{
  \item Lluís Revilla Sancho (\href{https://orcid.org/0000-0001-9747-2570}{ORCID})
  \item Hadley Wickham (\href{https://orcid.org/0000-0003-4757-117X}{ORCID})
}

Other contributors:
\itemize{
  \item Andrew Heiss (\href{https://orcid.org/0000-0002-3948-3914}{ORCID}) [reviewer]
  \item Francois Briatte [reviewer]
  \item Jonathan Sidi \email{yonicd@gmail.com} (\href{https://orcid.org/0000-0002-4222-1819}{ORCID}) [contributor]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/direct_messages.R
\name{direct_messages_received}
\alias{direct_messages_received}
\alias{direct_messages_sent}
\title{(DEPRECATED) Get the most recent direct messages sent to the authenticating user.}
\usage{
direct_messages_received(
  since_id = NULL,
  max_id = NULL,
  n = 200,
  parse = TRUE,
  token = NULL
)

direct_messages_sent(
  since_id = NULL,
  max_id = NULL,
  n = 200,
  parse = TRUE,
  token = NULL
)
}
\value{
Return object converted to nested list. If status code of
response object is not 200, the response object is returned
directly.
}
\description{
Retrieves up to 200 of the most recently received direct messages
by the authenticating (home) user. This function requires access
token with read, write, and direct messages access.
}
\details{
Includes detailed information about the sender and
recipient user. You can request up to 200 direct messages per
call, and only the most recent 200 direct messages will be available using
this endpoint.
}
\examples{

\dontrun{

## get my direct messages
dms <- direct_messages_received()

## inspect data structure
str(dms)

## get direct messages I've sent
sdms <- direct_messages_sent()

## inspect data structure
str(dms)

}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ts_plot.R
\name{ts_plot}
\alias{ts_plot}
\title{Plots tweets data as a time series-like data object.}
\usage{
ts_plot(data, by = "days", trim = 0L, tz = "UTC", ...)
}
\arguments{
\item{data}{Data frame or grouped data frame.}

\item{by}{Desired interval of time expressed as numeral plus one of
"secs", "mins", "hours", "days", "weeks", "months", or
"years". If a numeric is provided, the value is assumed to be in
seconds.}

\item{trim}{The number of observations to drop off the beginning
and end of the time series.}

\item{tz}{Time zone to be used, defaults to "UTC" (Twitter default)}

\item{...}{Other arguments passed to
\code{\link[ggplot2:geom_path]{ggplot2::geom_line()}}.}
}
\value{
If
\href{https://cran.r-project.org/package=ggplot2}{ggplot2} is
installed then a \code{\link[ggplot2:ggplot]{ggplot2::ggplot()}} plot object.
}
\description{
Creates a ggplot2 plot of the frequency of tweets over a specified
interval of time.
}
\examples{

\dontrun{
## search for tweets containing "rstats"
rt <- search_tweets("rstats", n = 10000)

## plot frequency in 1 min intervals
ts_plot(rt, "mins")

## plot multiple time series--retweets vs non-retweets
ts_plot(dplyr::group_by(tmls, is_retweet), "hours")

## compare account activity for some important US political figures
tmls <- get_timeline(
  c("SenSchumer", "SenGillibrand", "realDonaldTrump"),
  n = 3000
)

## examine all Twitter activity using weekly intervals
ts_plot(tmls, "weeks")

## group by screen name and plot each time series
ts_plot(dplyr::group_by(tmls, screen_name), "weeks")

## group by screen name and is_retweet
ts_plot(dplyr::group_by(tmls, screen_name, is_retweet), "months")

}
}
\concept{ts_data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{stopwordslangs}
\alias{stopwordslangs}
\title{Twitter stop words in multiple languages data.}
\format{
A tibble with three variables and 24,000 observations
}
\usage{
stopwordslangs
}
\description{
This data comes form a group of Twitter searches conducted at
several times during the calendar year of 2017. The data are
commonly observed words associated with 10 different languages,
including \code{c("ar", "en", "es", "fr", "in", "ja", "pt", "ru", "tr", "und")}. Variables include \code{"word"} (potential stop
words), \code{"lang"} (two or three word code), and \code{"p"}
(probability value associated with frequency position along a
normal distribution with higher values meaning the word occurs more
frequently and lower values meaning the words occur less
frequently).
}
\examples{
head(stopwordslangs)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/favorites.R
\name{get_favorites}
\alias{get_favorites}
\title{Get tweets favorited by one or more users}
\usage{
get_favorites(
  user,
  n = 200,
  since_id = NULL,
  max_id = NULL,
  parse = TRUE,
  retryonratelimit = NULL,
  verbose = TRUE,
  token = NULL
)
}
\arguments{
\item{user}{Character vector of screen names or user ids.
See \code{\link[=as_screenname]{as_screenname()}} for more details.}

\item{n}{Desired number of results to return. Results are downloaded
in pages when \code{n} is large; the default value will download a single
page. Set \code{n = Inf} to download as many results as possible.

The Twitter API rate limits the number of requests you can perform
in each 15 minute period. The easiest way to download more than that is
to use \code{retryonratelimit = TRUE}.

You are not guaranteed to get exactly \code{n} results back. You will get
fewer results when tweets have been deleted or if you hit a rate limit.
You will get more results if you ask for a number of tweets that's not
a multiple of page size, e.g. if you request \code{n = 150} and the page
size is 200, you'll get 200 results back.}

\item{since_id}{Supply a vector of ids or a data frame of previous results to
find tweets \strong{newer} than \code{since_id}.}

\item{max_id}{Supply a vector of ids or a data frame of previous results to
find tweets \strong{older} than \code{max_id}.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}

\item{retryonratelimit}{If \code{TRUE}, and a rate limit is exhausted, will wait
until it refreshes. Most Twitter rate limits refresh every 15 minutes.
If \code{FALSE}, and the rate limit is exceeded, the function will terminate
early with a warning; you'll still get back all results received up to
that point. The default value, \code{NULL}, consults the option
\code{rtweet.retryonratelimit} so that you can globally set it to \code{TRUE},
if desired.

If you expect a query to take hours or days to perform, you should not
rely soley on \code{retryonratelimit} because it does not handle other common
failure modes like temporarily losing your internet connection.}

\item{verbose}{Show progress bars and other messages indicating current
progress?}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}
}
\value{
A tibble with one row for each tweet.
}
\description{
Returns up to 3,000 tweets favorited for each \code{user}.
}
\examples{
\dontrun{
# get likes for a single user
kfc <- get_favorites("KFC")
kfc
# get newer likes since last request 
newer <- get_favorites("KFC", since_id = kfc)

# get likes from multiple users
favs <- get_favorites(c("Lesdoggg", "pattonoswalt", "meganamram"))
favs
}
}
\references{
\url{https://developer.twitter.com/en/docs/tweets/post-and-engage/api-reference/get-favorites-list}
}
\seealso{
Other tweets: 
\code{\link{get_mentions}()},
\code{\link{get_timeline}()},
\code{\link{lists_statuses}()},
\code{\link{lookup_tweets}()},
\code{\link{search_tweets}()}
}
\concept{tweets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/http.R
\name{TWIT_paginate_max_id}
\alias{TWIT_paginate_max_id}
\alias{TWIT_paginate_cursor}
\alias{TWIT_paginate_chunked}
\title{Pagination}
\usage{
TWIT_paginate_max_id(
  token,
  api,
  params,
  get_id = function(x) x$id_str,
  n = 1000,
  page_size = 200,
  since_id = NULL,
  max_id = NULL,
  count_param = "count",
  retryonratelimit = NULL,
  verbose = TRUE
)

TWIT_paginate_cursor(
  token,
  api,
  params,
  n = 5000,
  page_size = 5000,
  cursor = "-1",
  get_id = function(x) x$ids,
  retryonratelimit = NULL,
  verbose = TRUE
)

TWIT_paginate_chunked(
  token,
  api,
  params_list,
  retryonratelimit = NULL,
  verbose = TRUE
)
}
\arguments{
\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}

\item{get_id}{A single argument function that returns a vector of ids given
the JSON response. The defaults are chosen to cover the most common cases,
but you'll need to double check whenever implementing pagination for
a new endpoint.}

\item{n}{Desired number of results to return. Results are downloaded
in pages when \code{n} is large; the default value will download a single
page. Set \code{n = Inf} to download as many results as possible.

The Twitter API rate limits the number of requests you can perform
in each 15 minute period. The easiest way to download more than that is
to use \code{retryonratelimit = TRUE}.

You are not guaranteed to get exactly \code{n} results back. You will get
fewer results when tweets have been deleted or if you hit a rate limit.
You will get more results if you ask for a number of tweets that's not
a multiple of page size, e.g. if you request \code{n = 150} and the page
size is 200, you'll get 200 results back.}

\item{since_id}{Supply a vector of ids or a data frame of previous results to
find tweets \strong{newer} than \code{since_id}.}

\item{max_id}{Supply a vector of ids or a data frame of previous results to
find tweets \strong{older} than \code{max_id}.}

\item{retryonratelimit}{If \code{TRUE}, and a rate limit is exhausted, will wait
until it refreshes. Most Twitter rate limits refresh every 15 minutes.
If \code{FALSE}, and the rate limit is exceeded, the function will terminate
early with a warning; you'll still get back all results received up to
that point. The default value, \code{NULL}, consults the option
\code{rtweet.retryonratelimit} so that you can globally set it to \code{TRUE},
if desired.

If you expect a query to take hours or days to perform, you should not
rely soley on \code{retryonratelimit} because it does not handle other common
failure modes like temporarily losing your internet connection.}

\item{verbose}{Show progress bars and other messages indicating current
progress?}

\item{cursor}{Which page of results to return. The default will return
the first page; you can supply the result from a previous call to
continue pagination from where it left off.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}
These are internal functions used for pagination inside of rtweet.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{emojis}
\alias{emojis}
\title{Emojis codes and descriptions data.}
\format{
A tibble with two variables and 2,623 observations.
}
\usage{
emojis
}
\description{
This data comes from "Unicode.org",
\url{http://unicode.org/emoji/charts/full-emoji-list.html}. The data are
codes and descriptions of Emojis.
}
\examples{
head(emojis)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/friends.R
\name{get_friends}
\alias{get_friends}
\title{Get user IDs of accounts followed by target user(s).}
\usage{
get_friends(
  users,
  n = 5000,
  retryonratelimit = NULL,
  cursor = "-1",
  parse = TRUE,
  verbose = TRUE,
  token = NULL,
  page = lifecycle::deprecated()
)
}
\arguments{
\item{users}{Screen name or user ID of target user from which the
user IDs of friends (accounts followed BY target user) will be
retrieved.}

\item{n}{Desired number of results to return. Results are downloaded
in pages when \code{n} is large; the default value will download a single
page. Set \code{n = Inf} to download as many results as possible.

The Twitter API rate limits the number of requests you can perform
in each 15 minute period. The easiest way to download more than that is
to use \code{retryonratelimit = TRUE}.

You are not guaranteed to get exactly \code{n} results back. You will get
fewer results when tweets have been deleted or if you hit a rate limit.
You will get more results if you ask for a number of tweets that's not
a multiple of page size, e.g. if you request \code{n = 150} and the page
size is 200, you'll get 200 results back.}

\item{retryonratelimit}{If \code{TRUE}, and a rate limit is exhausted, will wait
until it refreshes. Most Twitter rate limits refresh every 15 minutes.
If \code{FALSE}, and the rate limit is exceeded, the function will terminate
early with a warning; you'll still get back all results received up to
that point. The default value, \code{NULL}, consults the option
\code{rtweet.retryonratelimit} so that you can globally set it to \code{TRUE},
if desired.

If you expect a query to take hours or days to perform, you should not
rely soley on \code{retryonratelimit} because it does not handle other common
failure modes like temporarily losing your internet connection.}

\item{cursor}{Which page of results to return. The default will return
the first page; you can supply the result from a previous call to
continue pagination from where it left off.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}

\item{verbose}{Show progress bars and other messages indicating current
progress?}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}

\item{page}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Please use \code{cursor} instead.}
}
\value{
A tibble data frame with two columns, "from_id" for name or ID of target
user and "to_id" for accounts ID they follow.
}
\description{
Returns a list of user IDs for the accounts following BY one or
more specified users.
}
\details{
Generally, you should not need to set \code{n} to more than 5,000 since Twitter
limits the number of people that you can follow (i.e. to follow more than
5,000 people at least 5,000 people need to follow you).
}
\examples{

\dontrun{
users <- get_friends("ropensci")
users
}
}
\references{
\url{https://developer.twitter.com/en/docs/accounts-and-users/follow-search-get-users/api-reference/get-friends-ids}

\url{https://developer.twitter.com/en/docs/twitter-api/v1/accounts-and-users/follow-search-get-users/api-reference/get-friends-ids}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post-user.R
\name{post_friendship}
\alias{post_friendship}
\alias{friendship_update}
\title{Updates friendship notifications and retweet abilities.}
\usage{
post_friendship(user, device = FALSE, retweets = FALSE, token = NULL)
}
\arguments{
\item{user}{Character vector of screen names or user ids.
See \code{\link[=as_screenname]{as_screenname()}} for more details.}

\item{device}{Logical indicating whether to enable or disable
device notifications from target user behaviors. Defaults
to false.}

\item{retweets}{Logical indicating whether to enable or disable
retweets from target user behaviors. Defaults to false.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}
}
\description{
Updates friendship notifications and retweet abilities.
}
\references{
\url{https://developer.twitter.com/en/docs/twitter-api/v1/accounts-and-users/follow-search-get-users/api-reference/post-friendships-update}
}
\seealso{
Other post: 
\code{\link{post_favorite}()},
\code{\link{post_follow}()},
\code{\link{post_tweet}()}
}
\concept{post}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/users.R
\name{lookup_users}
\alias{lookup_users}
\title{Get Twitter users data for given users (user IDs or screen names).}
\usage{
lookup_users(
  users,
  parse = TRUE,
  token = NULL,
  retryonratelimit = NULL,
  verbose = TRUE
)
}
\arguments{
\item{users}{User id or screen name of target user.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}

\item{retryonratelimit}{If \code{TRUE}, and a rate limit is exhausted, will wait
until it refreshes. Most Twitter rate limits refresh every 15 minutes.
If \code{FALSE}, and the rate limit is exceeded, the function will terminate
early with a warning; you'll still get back all results received up to
that point. The default value, \code{NULL}, consults the option
\code{rtweet.retryonratelimit} so that you can globally set it to \code{TRUE},
if desired.

If you expect a query to take hours or days to perform, you should not
rely soley on \code{retryonratelimit} because it does not handle other common
failure modes like temporarily losing your internet connection.}

\item{verbose}{Show progress bars and other messages indicating current
progress?}
}
\value{
A tibble of users data.
}
\description{
Get Twitter users data for given users (user IDs or screen names).
}
\examples{

\dontrun{
users <- c(
  "potus", "hillaryclinton", "realdonaldtrump",
  "fivethirtyeight", "cnn", "espn", "twitter"
)
users <- lookup_users(users)
users

# latest tweet from each user
tweets_data(users)
}

}
\references{
\url{https://developer.twitter.com/en/docs/accounts-and-users/follow-search-get-users/api-reference/get-users-lookup}
}
\seealso{
Other users: 
\code{\link{as_screenname}()},
\code{\link{lists_subscribers}()},
\code{\link{search_users}()}
}
\concept{users}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lists_memberships.R
\name{lists_memberships}
\alias{lists_memberships}
\title{Get Twitter list memberships (lists containing a given user)}
\usage{
lists_memberships(
  user = NULL,
  n = 200,
  cursor = "-1",
  filter_to_owned_lists = FALSE,
  token = NULL,
  parse = TRUE,
  retryonratelimit = NULL,
  verbose = TRUE,
  previous_cursor = NULL
)
}
\arguments{
\item{user}{Character vector of screen names or user ids.
See \code{\link[=as_screenname]{as_screenname()}} for more details.}

\item{n}{Desired number of results to return. Results are downloaded
in pages when \code{n} is large; the default value will download a single
page. Set \code{n = Inf} to download as many results as possible.

The Twitter API rate limits the number of requests you can perform
in each 15 minute period. The easiest way to download more than that is
to use \code{retryonratelimit = TRUE}.

You are not guaranteed to get exactly \code{n} results back. You will get
fewer results when tweets have been deleted or if you hit a rate limit.
You will get more results if you ask for a number of tweets that's not
a multiple of page size, e.g. if you request \code{n = 150} and the page
size is 200, you'll get 200 results back.}

\item{cursor}{Which page of results to return. The default will return
the first page; you can supply the result from a previous call to
continue pagination from where it left off.}

\item{filter_to_owned_lists}{When \code{TRUE}, will return only lists that
authenticating user owns.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}

\item{retryonratelimit}{If \code{TRUE}, and a rate limit is exhausted, will wait
until it refreshes. Most Twitter rate limits refresh every 15 minutes.
If \code{FALSE}, and the rate limit is exceeded, the function will terminate
early with a warning; you'll still get back all results received up to
that point. The default value, \code{NULL}, consults the option
\code{rtweet.retryonratelimit} so that you can globally set it to \code{TRUE},
if desired.

If you expect a query to take hours or days to perform, you should not
rely soley on \code{retryonratelimit} because it does not handle other common
failure modes like temporarily losing your internet connection.}

\item{verbose}{Show progress bars and other messages indicating current
progress?}

\item{previous_cursor}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Please use
\code{cursor} instead.}
}
\description{
Due to deleted or removed lists, the returned number of memberships
is often less than the provided n value. This is a reflection of the API and
not a unique quirk of rtweet.
}
\examples{
\dontrun{

## get up to 1000 Twitter lists that include Nate Silver
ns538 <- lists_memberships("NateSilver538", n = 1000)

## view data
ns538

}

}
\references{
\url{https://developer.twitter.com/en/docs/twitter-api/v1/accounts-and-users/create-manage-lists/api-reference/get-lists-memberships}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post-message.R
\name{post_message}
\alias{post_message}
\title{Posts direct message from user's Twitter account}
\usage{
post_message(text, user, media = NULL, token = NULL)
}
\arguments{
\item{text}{Character, text of message.}

\item{user}{Character vector of screen names or user ids.
See \code{\link[=as_screenname]{as_screenname()}} for more details.}

\item{media}{File path to image or video media to be
included in tweet.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}
}
\description{
Posts direct message from user's Twitter account
}
\references{
\url{https://developer.twitter.com/en/docs/twitter-api/v1/direct-messages/sending-and-receiving/api-reference/new-event}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auth.R
\name{auth_as}
\alias{auth_as}
\title{Set default authentication for the current session}
\usage{
auth_as(auth = NULL)
}
\arguments{
\item{auth}{One of the following options:
\itemize{
\item \code{NULL}, the default, will look for rtweet's "default" authentication
which uses your personal Twitter account. If it's not found, it will
call \code{\link[=auth_setup_default]{auth_setup_default()}} to set it up.
\item A string giving the name of a saved auth file made by \code{\link[=auth_save]{auth_save()}}.
\item An auth object created by \code{\link[=rtweet_app]{rtweet_app()}}, \code{\link[=rtweet_bot]{rtweet_bot()}}, or
\code{\link[=rtweet_user]{rtweet_user()}}.
}}
}
\value{
Invisibly returns the previous authentication mechanism.
}
\description{
\code{auth_as()} sets up the default authentication mechanism used by all
rtweet API calls. See \code{\link[=rtweet_user]{rtweet_user()}} to learn more about the three
available authentication options.
}
\examples{
\dontrun{
# Use app auth for the remainder of this session:
my_app <- rtweet_app()
auth_as(my_app)

# Switch back to the default user based auth
auth_as()

# Load auth saved by auth_save()
auth_as("my-saved-app")
}
}
\seealso{
Other authentication: 
\code{\link{auth_get}()},
\code{\link{auth_save}()},
\code{\link{auth_setup_default}()},
\code{\link{rtweet_user}()}
}
\concept{authentication}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/save_as_csv.R
\name{write_as_csv}
\alias{write_as_csv}
\alias{save_as_csv}
\title{Save Twitter data as a comma separated value file.}
\usage{
write_as_csv(x, file_name, prepend_ids = TRUE, na = "", fileEncoding = "UTF-8")

save_as_csv(x, file_name, prepend_ids = TRUE, na = "", fileEncoding = "UTF-8")
}
\arguments{
\item{x}{Data frame returned by an rtweet function.}

\item{file_name}{Desired name to save file as. If \code{file_name} does not
include the extension ".csv" it will be added automatically.}

\item{prepend_ids}{Logical indicating whether to prepend an "x"
before all Twitter IDs (for users, statuses, lists, etc.). It's
recommended when saving to CSV as these values otherwise get
treated as numeric and as a result the values are often less
precise due to rounding or other class-related quirks. Defaults
to true.}

\item{na}{Value to be used for missing (NA)s. Defaults to empty
character, "".}

\item{fileEncoding}{Encoding to be used when saving to
CSV. defaults to "UTF-8".}
}
\value{
Saved CSV files in current working directory.
}
\description{
Saves as flattened CSV file of Twitter data.
}
\seealso{
Other datafiles: 
\code{\link{flatten}()},
\code{\link{read_twitter_csv}()}

Other datafiles: 
\code{\link{flatten}()},
\code{\link{read_twitter_csv}()}
}
\concept{datafiles}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bearer_token.R
\name{invalidate_bearer}
\alias{invalidate_bearer}
\title{Invalidate bearer token}
\usage{
invalidate_bearer(token = NULL)
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
Invalidation of the bearer token is no longer the responsibility of rtweet.
This is something you should instead perform in the \href{https://developer.twitter.com/en/portal/projects-and-apps}{Twitter developer portal}.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auth.R
\name{auth_save}
\alias{auth_save}
\alias{auth_list}
\title{Save an authentication mechanism for use in a future session}
\usage{
auth_save(auth, name)

auth_list()
}
\arguments{
\item{auth}{One of \code{\link[=rtweet_app]{rtweet_app()}}, \code{\link[=rtweet_bot]{rtweet_bot()}}, or \code{\link[=rtweet_user]{rtweet_user()}}.}

\item{name}{Cache name to use.}
}
\description{
Use \code{auth_save()} with \code{\link[=auth_as]{auth_as()}} to avoid repeatedly entering app
credentials, making it easier to share auth between projects.
Use \code{auth_list()} to list all saved credentials.
}
\examples{
\dontrun{
# save app auth for use in other sessions
auth <- rtweet_app()
auth_save(auth, "my-app")

# later, in a different session...
auth_as("my-app")
}
}
\seealso{
Other authentication: 
\code{\link{auth_as}()},
\code{\link{auth_get}()},
\code{\link{auth_setup_default}()},
\code{\link{rtweet_user}()}
}
\concept{authentication}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tokens.R
\name{create_token}
\alias{create_token}
\title{Create custom Twitter OAuth token}
\usage{
create_token(
  app = "mytwitterapp",
  consumer_key = NULL,
  consumer_secret = NULL,
  access_token = NULL,
  access_secret = NULL,
  set_renv = TRUE
)
}
\arguments{
\item{app}{Name of user created Twitter application}

\item{consumer_key, consumer_secret}{App API key and secret.}

\item{access_token, access_secret}{Access token and secret.}

\item{set_renv}{Should the token be cached?}
}
\value{
Twitter OAuth token(s) (Token1.0).
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
By default, \code{create_token()} does three things: it creates an authentication
"token", sets it as the default token for the current session, and save it
to disk to make it the default for future sessions.

These three components have now been split up into three separate pieces:
use \code{\link[=rtweet_user]{rtweet_user()}}/\code{\link[=rtweet_app]{rtweet_app()}}/\code{\link[=rtweet_bot]{rtweet_bot()}} to create the token,
\code{\link[=auth_as]{auth_as()}} to make it the default for this session, and \code{\link[=auth_save]{auth_save()}} to
use it in future sessions. See \code{vignette("auth")} for full details.
}
\seealso{
Other tokens: 
\code{\link{get_token}()},
\code{\link{rate_limit}()}
}
\concept{tokens}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rate_limit.R
\name{rate_limit}
\alias{rate_limit}
\alias{rate_limit_reset}
\alias{rate_limit_wait}
\title{Rate limit helpers}
\usage{
rate_limit(resource_match = NULL, token = NULL)

rate_limit_reset(endpoint, token = NULL)

rate_limit_wait(endpoint, token = NULL)
}
\arguments{
\item{resource_match}{An optional regular expression used to filter the
resources listed in returned rate limit data.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}

\item{endpoint}{Name of Twitter endpoint like \code{"lookup/users"},
\code{"/media/upload"}, or \code{"/feedback/show/:id"}.
#' @references \url{https://developer.twitter.com/en/docs/developer-utilities/rate-limit-status/api-reference/get-application-rate_limit_status}}
}
\description{
\itemize{
\item \code{rate_limit()} returns a tibble of info about all rate limits
\item \code{rate_limit_reset()} returns the next reset for a endpoint
\item \code{rate_limit_wait()} waits for the next reset for an endpoint
}

You should not need to use these function in the usual operation of rtweet
because all paginated functions will wait on your behalf if you set
\code{retryonratelimit = TRUE}.
}
\examples{

\dontrun{
rate_limit()
}
}
\references{
\url{https://developer.twitter.com/en/docs/twitter-api/v1/developer-utilities/rate-limit-status/api-reference/get-application-rate_limit_status}
}
\seealso{
Other tokens: 
\code{\link{create_token}()},
\code{\link{get_token}()}
}
\concept{tokens}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auth.R
\name{auth_get}
\alias{auth_get}
\title{Get the current authentication mechanism}
\usage{
auth_get()
}
\description{
If no authentication has been set up for this session, \code{auth_get()} will
call \code{\link[=auth_as]{auth_as()}} to set it up.
}
\seealso{
Other authentication: 
\code{\link{auth_as}()},
\code{\link{auth_save}()},
\code{\link{auth_setup_default}()},
\code{\link{rtweet_user}()}
}
\concept{authentication}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/save_as_csv.R
\name{flatten}
\alias{flatten}
\alias{unflatten}
\title{flatten/unflatten data frame}
\usage{
flatten(x)

unflatten(x)
}
\arguments{
\item{x}{Data frame with list columns or converted-to-character (flattened)
columns.}
}
\value{
If flattened, then data frame where non-recursive list
columns---that is, list columns that contain only atomic, or non-list,
elements---have been converted to character vectors. If unflattened,
this function splits on spaces columns originally returned as lists
by functions in rtweet package. See details for more information.
}
\description{
Converts list columns that containing all atomic elements into
character vectors and vice versa (for appropriate named variables
according to the rtweet package)
}
\details{
If recursive list columns are contained within the data frame,
relevant columns will still be converted to atomic types but output
will also be accompanied with a warning message.

\code{flatten} flattens list columns by pasting them into a single string for
each observations. For example, a tweet that mentions four other users,
for the mentions_user_id variable, it will include the four user IDs
separated by a space.

`unflatten`` splits on spaces to convert into list columns any
columns with the following names: hashtags, symbols, urls_url,
urls_t.co, urls_expanded_url, media_url, media_t.co,
media_expanded_url, media_type, ext_media_url, ext_media_t.co,
ext_media_expanded_url, mentions_user_id, mentions_screen_name,
geo_coords, coords_coords, bbox_coords, mentions_screen_name
}
\seealso{
Other datafiles: 
\code{\link{read_twitter_csv}()},
\code{\link{write_as_csv}()}

Other datafiles: 
\code{\link{read_twitter_csv}()},
\code{\link{write_as_csv}()}
}
\concept{datafiles}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stream.R
\name{parse_stream}
\alias{parse_stream}
\title{Converts Twitter stream data (JSON file) into parsed data frame.}
\usage{
parse_stream(path, ...)
}
\arguments{
\item{path}{Character, name of JSON file with data collected by
\code{\link[=stream_tweets]{stream_tweets()}}.}

\item{...}{Other arguments passed on to internal data_from_stream
function.}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
Please use \code{jsonlite::stream_in()} instead.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ts_plot.R
\name{ts_data}
\alias{ts_data}
\title{Converts tweets data into time series-like data object.}
\usage{
ts_data(data, by = "days", trim = 0L, tz = "UTC")
}
\arguments{
\item{data}{Data frame or grouped data frame.}

\item{by}{Desired interval of time expressed as numeral plus one of
"secs", "mins", "hours", "days", "weeks", "months", or
"years". If a numeric is provided, the value is assumed to be in
seconds.}

\item{trim}{Number of observations to trim off the front and end of
each time series}

\item{tz}{Time zone to be used, defaults to "UTC" (Twitter default)}
}
\value{
Data frame with time, n, and grouping column if applicable.
}
\description{
Returns data containing the frequency of tweets over a specified
interval of time.
}
\examples{

\dontrun{

## handles of women senators
sens <- c("SenatorBaldwin", "SenGillibrand", "PattyMurray", "SenatorHeitkamp")

## get timelines for each
sens <- get_timeline(sens, n = 3200)

## get single time series for tweets
ts_data(sens)

## using weekly intervals
ts_data(sens, "weeks")

## group by screen name and then use weekly intervals
ts_plot(dplyr::group_by(sens, screen_name), "weeks")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trends.R
\name{trends_available}
\alias{trends_available}
\title{Available Twitter trends along with associated WOEID.}
\usage{
trends_available(token = NULL, parse = TRUE)
}
\arguments{
\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}
}
\value{
Data frame with WOEID column. WOEID is a Yahoo! Where On
Earth ID.
}
\description{
Available Twitter trends along with associated WOEID.
}
\examples{
\dontrun{
## Retrieve available trends
trends <- trends_available()
trends

}

}
\references{
\url{https://developer.twitter.com/en/docs/twitter-api/v1/trends/locations-with-trending-topics/api-reference/get-trends-available}
}
\seealso{
Other trends: 
\code{\link{get_trends}()}
}
\concept{trends}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/next_cursor.R
\name{max_id}
\alias{max_id}
\title{Extract min/max id (for id based pagination)}
\usage{
max_id(x)
}
\arguments{
\item{x}{Either a data frame of tweets or a character vector of status ids.}
}
\description{
These internal helpers extract the ids passed on to the \code{max_id}
and \code{since_id} arguments to functions that use \code{\link[=TWIT_paginate_max_id]{TWIT_paginate_max_id()}}.
}
\examples{
\dontrun{
tw <- search_tweets("#rstats")

# retrieve older tweets
older <- search_tweets("#rstats", max_id = tw)
even_older <- search_tweets("#rstats", max_id = older)

# retrieve newer tweets
newer <- search_tweets("#rstats", since_id = tw)
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/followers.R
\name{get_followers}
\alias{get_followers}
\title{Get user IDs for accounts following target user.}
\usage{
get_followers(
  user,
  n = 5000,
  cursor = "-1",
  retryonratelimit = NULL,
  parse = TRUE,
  verbose = TRUE,
  token = NULL,
  page = lifecycle::deprecated()
)
}
\arguments{
\item{user}{Character vector of screen names or user ids.
See \code{\link[=as_screenname]{as_screenname()}} for more details.}

\item{n}{Desired number of results to return. Results are downloaded
in pages when \code{n} is large; the default value will download a single
page. Set \code{n = Inf} to download as many results as possible.

The Twitter API rate limits the number of requests you can perform
in each 15 minute period. The easiest way to download more than that is
to use \code{retryonratelimit = TRUE}.

You are not guaranteed to get exactly \code{n} results back. You will get
fewer results when tweets have been deleted or if you hit a rate limit.
You will get more results if you ask for a number of tweets that's not
a multiple of page size, e.g. if you request \code{n = 150} and the page
size is 200, you'll get 200 results back.}

\item{cursor}{Which page of results to return. The default will return
the first page; you can supply the result from a previous call to
continue pagination from where it left off.}

\item{retryonratelimit}{If \code{TRUE}, and a rate limit is exhausted, will wait
until it refreshes. Most Twitter rate limits refresh every 15 minutes.
If \code{FALSE}, and the rate limit is exceeded, the function will terminate
early with a warning; you'll still get back all results received up to
that point. The default value, \code{NULL}, consults the option
\code{rtweet.retryonratelimit} so that you can globally set it to \code{TRUE},
if desired.

If you expect a query to take hours or days to perform, you should not
rely soley on \code{retryonratelimit} because it does not handle other common
failure modes like temporarily losing your internet connection.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}

\item{verbose}{Show progress bars and other messages indicating current
progress?}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}

\item{page}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Please use \code{cursor} instead.}
}
\value{
A tibble data frame with one column named "from_id" with the
followers and another one "to_id" with the user used as input.
}
\description{
Returns a list of user IDs for the accounts following specified
user.
}
\examples{

\dontrun{
users <- get_followers("KFC")
users

# use `cursor` to find the next "page" of results
more_users <- get_followers("KFC", cursor = users)

}
}
\references{
\url{https://developer.twitter.com/en/docs/accounts-and-users/follow-search-get-users/api-reference/get-followers-ids}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tweet_embed.R
\name{tweet_embed}
\alias{tweet_embed}
\title{Create a Tweet Embed}
\usage{
tweet_embed(screen_name, status_id, ...)
}
\arguments{
\item{screen_name}{character, screen name of the user}

\item{status_id}{character, status id}

\item{...}{parameters to pass to the GET call. See
\url{https://developer.twitter.com/en/docs/tweets/post-and-engage/api-reference/get-statuses-oembed}
for details.}
}
\value{
character
}
\description{
Twitter API GET call to retieve the tweet in embedded form.
}
\examples{
name   <- 'kearneymw'
status <- '1087047171306856451'

tweet_embed(screen_name = name, status_id = status)

tweet_embed(
 screen_name = name,
 status_id = status,
 hide_thread = TRUE, 
 hide_media = FALSE, 
 align = 'center'
)

}
\seealso{
\code{\link[httr:GET]{httr::GET()}},\code{\link[httr:content]{httr::content()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/friends.R
\name{my_friendships}
\alias{my_friendships}
\title{Lookup friendship information between users.}
\usage{
my_friendships(user, parse = FALSE, token = NULL)
}
\arguments{
\item{user}{Character vector of screen names or user ids.
See \code{\link[=as_screenname]{as_screenname()}} for more details.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}
}
\description{
Gets information on friendship between authenticated user and up
to 100 other users.
}
\references{
\url{https://developer.twitter.com/en/docs/twitter-api/v1/accounts-and-users/follow-search-get-users/api-reference/get-friendships-lookup}
}
\seealso{
Other friends: 
\code{\link{lookup_friendships}()}
}
\concept{friends}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stream.R
\name{stream_tweets}
\alias{stream_tweets}
\title{Collect a live stream of Twitter data}
\usage{
stream_tweets(
  q = "",
  timeout = 30,
  parse = TRUE,
  token = NULL,
  file_name = NULL,
  verbose = TRUE,
  append = TRUE,
  ...
)
}
\arguments{
\item{q}{Query used to select and customize streaming collection
method.  There are four possible methods:
\enumerate{
\item The default, \code{q = ""}, returns a small random sample of all
publicly available Twitter statuses.
\item To filter by keyword, provide a comma separated character string with
the desired phrase(s) and keyword(s).
\item Track users by providing a comma separated list of user IDs or
screen names.
\item Use four latitude/longitude bounding box points to stream by geo
location. This must be provided via a vector of length 4, e.g.,
\code{c(-125, 26, -65, 49)}.
}}

\item{timeout}{Integer specifying number of seconds to stream tweets for.
Stream indefinitely with \code{timeout = Inf}.

The stream can be interrupted at any time, and \code{file_name} will still be
valid file.}

\item{parse}{Use \code{FALSE} to opt-out of parsing the tweets.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}

\item{file_name}{Character with name of file. If not specified,
will write to \code{stream_tweets.json} in the current working directory.}

\item{verbose}{If \code{TRUE}, display a progress bar.}

\item{append}{If \code{TRUE}, will append to the end of \code{file_name}; if
\code{FALSE}, will overwrite.}

\item{...}{Other arguments passed in to query parameters.}
}
\value{
A tibble with one row per tweet
}
\description{
Streams public statuses to a file via one of the following four methods:
\enumerate{
\item Sampling a small random sample of all publicly available tweets
\item Filtering via a search-like query (up to 400 keywords)
\item Tracking via vector of user ids (up to 5000 user_ids)
\item Location via geo coordinates (1-360 degree location boxes)
}

Learn more in \code{vignette("stream", package = "rtweet")}
}
\examples{
\dontrun{
# stream tweets mentioning "election" for 10 seconds
e <- stream_tweets("election", timeout = 10)
e

# Download another 10s worth of data to the same file
e <- stream_tweets("election", timeout = 10, append = TRUE)

# stream tweets about continential USA for 5 minutes
usa <- stream_tweets(location = lookup_coords("usa"), file_name = "usa.json", timeout = 300)

}
}
\references{
\url{https://developer.twitter.com/en/docs/twitter-api/v1/tweets/sample-realtime/api-reference/get-statuses-sample},
\url{https://developer.twitter.com/en/docs/twitter-api/v1/tweets/filter-realtime/overview}

Stream: \url{https://developer.twitter.com/en/docs/twitter-api/v1/tweets/sample-realtime/api-reference/get-statuses-sample}
Filter: \url{https://developer.twitter.com/en/docs/twitter-api/v1/tweets/filter-realtime/api-reference/post-statuses-filter}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post-favorite.R
\name{post_favorite}
\alias{post_favorite}
\alias{post_favourite}
\alias{favorite_tweet}
\title{Favorites target status id.}
\usage{
post_favorite(
  status_id,
  destroy = FALSE,
  include_entities = FALSE,
  token = NULL
)
}
\arguments{
\item{status_id}{Status id of target tweet.}

\item{destroy}{Logical indicating whether to post (add) or
remove (delete) target tweet as favorite.}

\item{include_entities}{Logical indicating whether to
include entities object in return.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}
}
\description{
Favorites target status id.
}
\examples{
\dontrun{
rt <- search_tweets("#rstats", n = 1)
post_favorite(rt$status_id)
}
}
\references{
Create: \url{https://developer.twitter.com/en/docs/twitter-api/v1/tweets/post-and-engage/api-reference/post-favorites-create}
Destroy: \url{https://developer.twitter.com/en/docs/twitter-api/v1/tweets/post-and-engage/api-reference/post-favorites-destroy}
}
\seealso{
Other post: 
\code{\link{post_follow}()},
\code{\link{post_friendship}()},
\code{\link{post_tweet}()}
}
\concept{post}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/timeline.R
\name{get_timeline}
\alias{get_timeline}
\alias{get_my_timeline}
\alias{get_timelines}
\title{Get one or more user timelines}
\usage{
get_timeline(
  user = NULL,
  n = 100,
  since_id = NULL,
  max_id = NULL,
  home = FALSE,
  parse = TRUE,
  check = TRUE,
  retryonratelimit = NULL,
  verbose = TRUE,
  token = NULL,
  ...
)

get_my_timeline(
  n = 100,
  since_id = NULL,
  max_id = NULL,
  parse = TRUE,
  check = TRUE,
  retryonratelimit = NULL,
  verbose = TRUE,
  token = NULL,
  ...
)
}
\arguments{
\item{user}{Character vector of screen names or user ids.
See \code{\link[=as_screenname]{as_screenname()}} for more details.}

\item{n}{Desired number of results to return. Results are downloaded
in pages when \code{n} is large; the default value will download a single
page. Set \code{n = Inf} to download as many results as possible.

The Twitter API rate limits the number of requests you can perform
in each 15 minute period. The easiest way to download more than that is
to use \code{retryonratelimit = TRUE}.

You are not guaranteed to get exactly \code{n} results back. You will get
fewer results when tweets have been deleted or if you hit a rate limit.
You will get more results if you ask for a number of tweets that's not
a multiple of page size, e.g. if you request \code{n = 150} and the page
size is 200, you'll get 200 results back.}

\item{since_id}{Supply a vector of ids or a data frame of previous results to
find tweets \strong{newer} than \code{since_id}.}

\item{max_id}{Supply a vector of ids or a data frame of previous results to
find tweets \strong{older} than \code{max_id}.}

\item{home}{Logical, indicating whether to return a "user" timeline
(the default, what a user has tweeted/retweeted) or a "home" timeline
(what the user would see if they logged into twitter).}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}

\item{check}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}}

\item{retryonratelimit}{If \code{TRUE}, and a rate limit is exhausted, will wait
until it refreshes. Most Twitter rate limits refresh every 15 minutes.
If \code{FALSE}, and the rate limit is exceeded, the function will terminate
early with a warning; you'll still get back all results received up to
that point. The default value, \code{NULL}, consults the option
\code{rtweet.retryonratelimit} so that you can globally set it to \code{TRUE},
if desired.

If you expect a query to take hours or days to perform, you should not
rely soley on \code{retryonratelimit} because it does not handle other common
failure modes like temporarily losing your internet connection.}

\item{verbose}{Show progress bars and other messages indicating current
progress?}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}

\item{...}{Further arguments passed on as parameters in API query.}
}
\value{
A tbl data frame of tweets data with users data attribute.
}
\description{
\code{get_timeline()} returns the timeline of any Twitter user (i.e. what they
have tweeted). \code{get_my_timeline()} returns the home timeline for the
authenticated user (i.e. the tweets you see when you log into Twitter).
}
\examples{

\dontrun{
tw <- get_timeline("JustinBieber")
tw

# get tweets that arrived since the last request
get_timeline("JustinBieber", since_id = tw)
# get earlier tweets
get_timeline("JustinBieber", max_id = tw)

# get timelines for multiple users
tw <- get_timeline(c("KFC", "PizzaHut", "McDonalds"))
tw
}

}
\references{
\url{https://developer.twitter.com/en/docs/tweets/timelines/api-reference/get-statuses-user_timeline}
}
\seealso{
Other tweets: 
\code{\link{get_favorites}()},
\code{\link{get_mentions}()},
\code{\link{lists_statuses}()},
\code{\link{lookup_tweets}()},
\code{\link{search_tweets}()}
}
\concept{tweets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_tweets.R
\name{search_tweets}
\alias{search_tweets}
\alias{search_tweets2}
\title{Get tweets data on statuses identified via search query.}
\usage{
search_tweets(
  q,
  n = 100,
  type = c("mixed", "recent", "popular"),
  include_rts = TRUE,
  geocode = NULL,
  since_id = NULL,
  max_id = NULL,
  parse = TRUE,
  token = NULL,
  retryonratelimit = NULL,
  verbose = TRUE,
  ...
)

search_tweets2(...)
}
\arguments{
\item{q}{Query to be searched, used to filter and select tweets to
return from Twitter's REST API. Must be a character string not to
exceed maximum of 500 characters. Spaces behave like boolean
"AND" operator. To search for tweets containing at least one of
multiple possible terms, separate each search term with spaces
and "OR" (in caps). For example, the search \code{q = "data science"} looks for tweets containing both "data" and
"science" located anywhere in the tweets and in any order.
When "OR" is entered between search terms, \code{query = "data OR science"}, Twitter's REST API should return any tweet
that contains either "data" or "science." It is also possible to
search for exact phrases using double quotes. To do this, either
wrap single quotes around a search query using double quotes,
e.g., \code{q = '"data science"'} or escape each internal double
quote with a single backslash, e.g., \verb{q = "\\"data science\\""}.

Some other useful query tips:

\itemize{
\item Exclude retweets via \code{"-filter:retweets"}
\item Exclude quotes via \code{"-filter:quote"}
\item Exclude replies via \code{"-filter:replies"}
\item Filter (return only) verified via \code{"filter:verified"}
\item Exclude verified via \code{"-filter:verified"}
\item Get everything (firehose for free) via \code{"-filter:verified OR filter:verified"}
\item Filter (return only) tweets with links to news articles via \code{"filter:news"}
\item Filter (return only) tweets with media \code{"filter:media"}
}}

\item{n}{Desired number of results to return. Results are downloaded
in pages when \code{n} is large; the default value will download a single
page. Set \code{n = Inf} to download as many results as possible.

The Twitter API rate limits the number of requests you can perform
in each 15 minute period. The easiest way to download more than that is
to use \code{retryonratelimit = TRUE}.

You are not guaranteed to get exactly \code{n} results back. You will get
fewer results when tweets have been deleted or if you hit a rate limit.
You will get more results if you ask for a number of tweets that's not
a multiple of page size, e.g. if you request \code{n = 150} and the page
size is 200, you'll get 200 results back.}

\item{type}{Character string specifying which type of search
results to return from Twitter's REST API. The current default is
\code{type = "recent"}, other valid types include \code{type = "mixed"} and \code{type = "popular"}.}

\item{include_rts}{Logical, indicating whether to include retweets
in search results. Retweets are classified as any tweet generated
by Twitter's built-in "retweet" (recycle arrows) function. These
are distinct from quotes (retweets with additional text provided
from sender) or manual retweets (old school method of manually
entering "RT" into the text of one's tweets).}

\item{geocode}{Geographical limiter of the template
"latitude,longitude,radius" e.g., \code{geocode = "37.78,-122.40,1mi"}.}

\item{since_id}{Supply a vector of ids or a data frame of previous results to
find tweets \strong{newer} than \code{since_id}.}

\item{max_id}{Supply a vector of ids or a data frame of previous results to
find tweets \strong{older} than \code{max_id}.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}

\item{retryonratelimit}{If \code{TRUE}, and a rate limit is exhausted, will wait
until it refreshes. Most Twitter rate limits refresh every 15 minutes.
If \code{FALSE}, and the rate limit is exceeded, the function will terminate
early with a warning; you'll still get back all results received up to
that point. The default value, \code{NULL}, consults the option
\code{rtweet.retryonratelimit} so that you can globally set it to \code{TRUE},
if desired.

If you expect a query to take hours or days to perform, you should not
rely soley on \code{retryonratelimit} because it does not handle other common
failure modes like temporarily losing your internet connection.}

\item{verbose}{Show progress bars and other messages indicating current
progress?}

\item{...}{Further arguments passed as query parameters in request
sent to Twitter's REST API. To return only English language
tweets, for example, use \code{lang = "en"}. For more options see
Twitter's API documentation.}
}
\value{
List object with tweets and users each returned as a
data frame.

A tbl data frame with additional "query" column.
}
\description{
Returns Twitter statuses matching a user provided search
query. ONLY RETURNS DATA FROM THE PAST 6-9 DAYS.

search_tweets2 Passes all arguments to search_tweets. Returns data from
one OR MORE search queries.
}
\details{
Twitter API documentation recommends limiting searches to
10 keywords and operators. Complex queries may also produce API
errors preventing recovery of information related to the query.
It should also be noted Twitter's search API does not consist of
an index of all Tweets. At the time of searching, the search API
index includes between only 6-9 days of Tweets.
}
\examples{

\dontrun{
tweets <- search_tweets("weather")
tweets

# data about the users who made those tweets
users_data(tweets)

# Retrieve all the tweets made since the previous request
# (there might not be any if people aren't tweeting about the weather)
newer <- search_tweets("weather", since_id = tweets)
# Retrieve tweets made before the previous request
older <- search_tweets("weather", max_id = tweets)

# Restrict to English only, and ignore retweets
tweets2 <- search_tweets("weather", lang = "en", include_rts = FALSE)
}

\dontrun{

## search using multilple queries
st2 <- search_tweets2(
  c("\"data science\"", "rstats OR python"),
  n = 500
)

## preview tweets data
st2

## preview users data
users_data(st2)

## check breakdown of results by search query
table(st2$query)

}

}
\references{
\url{https://developer.twitter.com/en/docs/tweets/search/api-reference/get-search-tweets}

\url{https://developer.twitter.com/en/docs/twitter-api/v1/tweets/search/api-reference/get-search-tweets}
}
\seealso{
Other tweets: 
\code{\link{get_favorites}()},
\code{\link{get_mentions}()},
\code{\link{get_timeline}()},
\code{\link{lists_statuses}()},
\code{\link{lookup_tweets}()}
}
\concept{tweets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lists_subscribers.R
\name{lists_subscribers}
\alias{lists_subscribers}
\title{Get subscribers of a specified list.}
\usage{
lists_subscribers(
  list_id = NULL,
  slug = NULL,
  owner_user = NULL,
  n = 20,
  cursor = "-1",
  parse = TRUE,
  retryonratelimit = NULL,
  verbose = TRUE,
  token = NULL
)
}
\arguments{
\item{list_id}{required The numerical id of the list.}

\item{slug, owner_user}{The list name (slug) and owner.}

\item{n}{Desired number of results to return. Results are downloaded
in pages when \code{n} is large; the default value will download a single
page. Set \code{n = Inf} to download as many results as possible.

The Twitter API rate limits the number of requests you can perform
in each 15 minute period. The easiest way to download more than that is
to use \code{retryonratelimit = TRUE}.

You are not guaranteed to get exactly \code{n} results back. You will get
fewer results when tweets have been deleted or if you hit a rate limit.
You will get more results if you ask for a number of tweets that's not
a multiple of page size, e.g. if you request \code{n = 150} and the page
size is 200, you'll get 200 results back.}

\item{cursor}{Which page of results to return. The default will return
the first page; you can supply the result from a previous call to
continue pagination from where it left off.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}

\item{retryonratelimit}{If \code{TRUE}, and a rate limit is exhausted, will wait
until it refreshes. Most Twitter rate limits refresh every 15 minutes.
If \code{FALSE}, and the rate limit is exceeded, the function will terminate
early with a warning; you'll still get back all results received up to
that point. The default value, \code{NULL}, consults the option
\code{rtweet.retryonratelimit} so that you can globally set it to \code{TRUE},
if desired.

If you expect a query to take hours or days to perform, you should not
rely soley on \code{retryonratelimit} because it does not handle other common
failure modes like temporarily losing your internet connection.}

\item{verbose}{Show progress bars and other messages indicating current
progress?}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}
}
\description{
Get subscribers of a specified list.
}
\examples{

\dontrun{

## get subscribers of new york times politics list
rstats <- lists_subscribers(
  slug = "new-york-times-politics",
  owner_user = "nytpolitics",
  n = 1000
)

}

}
\references{
\url{https://developer.twitter.com/en/docs/twitter-api/v1/accounts-and-users/create-manage-lists/api-reference/get-lists-subscribers}
}
\seealso{
Other lists: 
\code{\link{lists_members}()},
\code{\link{lists_statuses}()},
\code{\link{lists_subscriptions}()},
\code{\link{lists_users}()}

Other users: 
\code{\link{as_screenname}()},
\code{\link{lookup_users}()},
\code{\link{search_users}()}
}
\concept{lists}
\concept{users}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lists_users.R
\name{lists_users}
\alias{lists_users}
\title{Get all lists a specified user subscribes to, including their own.}
\usage{
lists_users(user = NULL, reverse = FALSE, token = NULL, parse = TRUE)
}
\arguments{
\item{user}{Character vector of screen names or user ids.
See \code{\link[=as_screenname]{as_screenname()}} for more details.}

\item{reverse}{optional Set this to true if you would like owned lists
to be returned first. See description above for information on
how this parameter works.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}
}
\value{
data
}
\description{
Get all lists a specified user subscribes to, including their own.
}
\examples{
\dontrun{

## get lists subsribed to by Nate Silver
lists_users("NateSilver538")

}

}
\seealso{
Other lists: 
\code{\link{lists_members}()},
\code{\link{lists_statuses}()},
\code{\link{lists_subscribers}()},
\code{\link{lists_subscriptions}()}
}
\concept{lists}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lists_members.R
\name{lists_members}
\alias{lists_members}
\title{Get Twitter list members (users on a given list).}
\usage{
lists_members(
  list_id = NULL,
  slug = NULL,
  owner_user = NULL,
  n = 5000,
  cursor = "-1",
  token = NULL,
  retryonratelimit = NULL,
  verbose = TRUE,
  parse = TRUE,
  ...
)
}
\arguments{
\item{list_id}{required The numerical id of the list.}

\item{slug}{required You can identify a list by its slug instead of
its numerical id. If you decide to do so, note that you'll also
have to specify the list owner using the owner_id or
owner_user parameters.}

\item{owner_user}{optional The screen name or user ID of the user}

\item{n}{Desired number of results to return. Results are downloaded
in pages when \code{n} is large; the default value will download a single
page. Set \code{n = Inf} to download as many results as possible.

The Twitter API rate limits the number of requests you can perform
in each 15 minute period. The easiest way to download more than that is
to use \code{retryonratelimit = TRUE}.

You are not guaranteed to get exactly \code{n} results back. You will get
fewer results when tweets have been deleted or if you hit a rate limit.
You will get more results if you ask for a number of tweets that's not
a multiple of page size, e.g. if you request \code{n = 150} and the page
size is 200, you'll get 200 results back.}

\item{cursor}{Which page of results to return. The default will return
the first page; you can supply the result from a previous call to
continue pagination from where it left off.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}

\item{retryonratelimit}{If \code{TRUE}, and a rate limit is exhausted, will wait
until it refreshes. Most Twitter rate limits refresh every 15 minutes.
If \code{FALSE}, and the rate limit is exceeded, the function will terminate
early with a warning; you'll still get back all results received up to
that point. The default value, \code{NULL}, consults the option
\code{rtweet.retryonratelimit} so that you can globally set it to \code{TRUE},
if desired.

If you expect a query to take hours or days to perform, you should not
rely soley on \code{retryonratelimit} because it does not handle other common
failure modes like temporarily losing your internet connection.}

\item{verbose}{Show progress bars and other messages indicating current
progress?}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}

\item{...}{Other arguments used as parameters in query composition.}
}
\description{
Get Twitter list members (users on a given list).
}
\examples{
\dontrun{

## get list members for a list of polling experts using list_id
(pollsters <- lists_members("105140588"))

## get list members of cspan's senators list
sens <- lists_members(slug = "senators", owner_user = "cspan")
sens

## get list members for an rstats list using list topic slug
## list owner's screen name
rstats <- lists_members(slug = "rstats", owner_user = "scultrera")
rstats

}

}
\references{
\url{https://developer.twitter.com/en/docs/twitter-api/v1/accounts-and-users/create-manage-lists/api-reference/get-lists-members}
}
\seealso{
Other lists: 
\code{\link{lists_statuses}()},
\code{\link{lists_subscribers}()},
\code{\link{lists_subscriptions}()},
\code{\link{lists_users}()}
}
\concept{lists}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lat_lng.R
\name{lat_lng}
\alias{lat_lng}
\title{Adds single-point latitude and longitude variables to tweets data.}
\usage{
lat_lng(x, coords = c("coords_coords", "bbox_coords", "geo_coords"))
}
\arguments{
\item{x}{Parsed Twitter data as returned by various rtweet
functions. This should be a data frame with variables such as
"bbox_coords", "coords_coords", and "geo_coords" (among
other non-geolocation Twitter variables).}

\item{coords}{Names of variables containing latitude and longitude
coordinates.  Priority is given to bounding box coordinates (each
obs consists of eight entries) followed by the supplied order of
variable names. Defaults to "bbox_coords",
"coords_coords", and "geo_coords") (which are the default column
names of data returned by most status-oriented rtweet functions).}
}
\value{
Returns updated data object with full information latitude
and longitude vars.
}
\description{
Appends parsed Twitter data with latitude and longitude variables
using all available geolocation information.
}
\details{
On occasion values may appear to be outliers given a
previously used query filter (e.g., when searching for tweets
sent from the continental US).  This is typically because those
tweets returned a large bounding box that overlapped with the
area of interest. This function converts boxes into their
geographical midpoints, which works well in the vast majority of
cases, but sometimes includes an otherwise puzzling result.
}
\examples{

\dontrun{

## stream tweets sent from the US
rt <- stream_tweets(lookup_coords("usa"), timeout = 10)

## use lat_lng to recover full information geolocation data
rtll <- lat_lng(rt)

## plot points
with(rtll, plot(lng, lat))

}

}
\seealso{
Other geo: 
\code{\link{lookup_coords}()}
}
\concept{geo}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tweet_shot.R
\name{tweet_shot}
\alias{tweet_shot}
\title{Capture an image of a tweet/thread}
\usage{
tweet_shot(statusid_or_url, zoom = 3, scale = TRUE)
}
\arguments{
\item{statusid_or_url}{a valid Twitter status id (e.g. "\code{947082036019388416}") or
a valid Twitter status URL (e.g. "\verb{https://twitter.com/jhollist/status/947082036019388416}").}

\item{zoom}{a positive number >= 1. See the help for \verb{[webshot::webshot()]} for more information.}

\item{scale}{auto-scale the image back to 1:1? Default it \code{TRUE}, which means \code{magick}
will be used to return a "normal" sized tweet. Set it to \code{FALSE} to perform your
own image manipulation.}
}
\value{
\code{magick} object
}
\description{
Provide a status id or a full Twitter link to a tweet and this function
will capture an image of the tweet --- or tweet + thread (if there are
Twitter-linked replies) --- from the mobile version of said tweet/thread.
}
\details{
For this to work, you will need to ensure the packages in \verb{Suggests:} are
installed as they will be loaded upon the first invocation of this function.

Use the \code{zoom} factor to get more pixels which may improve the text rendering
of the tweet/thread.
}
\examples{
\dontrun{
tweet_shot("947082036019388416")
tweet_shot("https://twitter.com/jhollist/status/947082036019388416")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/retweets.R
\name{get_retweets}
\alias{get_retweets}
\alias{get_retweeters}
\title{Get the most recent retweets/retweeters}
\usage{
get_retweets(status_id, n = 100, parse = TRUE, token = NULL, ...)

get_retweeters(status_id, n = 100, parse = TRUE, token = NULL)
}
\arguments{
\item{status_id}{Tweet/status id.}

\item{n}{Number of results to retrieve. Must be <= 100.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}

\item{...}{Other arguments used as parameters in the query sent to
Twitter's rest API, for example, \code{trim_user = TRUE}.}
}
\value{
Tweets data of the most recent retweets/retweeters of a given status
}
\description{
\code{get_retweets()} returns the 100 most recent retweets of a tweet;
\code{get_retweeters()} retursn the 100 most recent users who retweeted them.
}
\references{
\url{https://developer.twitter.com/en/docs/twitter-api/v1/tweets/post-and-engage/api-reference/get-statuses-retweets-id}

\url{https://developer.twitter.com/en/docs/twitter-api/v1/tweets/post-and-engage/api-reference/get-statuses-retweeters-ids}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/direct_messages.R
\name{direct_messages}
\alias{direct_messages}
\title{Get direct messages sent to and received by the authenticating user from the
past 30 days}
\usage{
direct_messages(
  n = 50,
  cursor = NULL,
  next_cursor = NULL,
  parse = TRUE,
  token = NULL,
  retryonratelimit = NULL,
  verbose = TRUE
)
}
\arguments{
\item{n}{Desired number of results to return. Results are downloaded
in pages when \code{n} is large; the default value will download a single
page. Set \code{n = Inf} to download as many results as possible.

The Twitter API rate limits the number of requests you can perform
in each 15 minute period. The easiest way to download more than that is
to use \code{retryonratelimit = TRUE}.

You are not guaranteed to get exactly \code{n} results back. You will get
fewer results when tweets have been deleted or if you hit a rate limit.
You will get more results if you ask for a number of tweets that's not
a multiple of page size, e.g. if you request \code{n = 150} and the page
size is 200, you'll get 200 results back.}

\item{cursor}{Which page of results to return. The default will return
the first page; you can supply the result from a previous call to
continue pagination from where it left off.}

\item{next_cursor}{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} Use \code{cursor} instead.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}

\item{retryonratelimit}{If \code{TRUE}, and a rate limit is exhausted, will wait
until it refreshes. Most Twitter rate limits refresh every 15 minutes.
If \code{FALSE}, and the rate limit is exceeded, the function will terminate
early with a warning; you'll still get back all results received up to
that point. The default value, \code{NULL}, consults the option
\code{rtweet.retryonratelimit} so that you can globally set it to \code{TRUE},
if desired.

If you expect a query to take hours or days to perform, you should not
rely soley on \code{retryonratelimit} because it does not handle other common
failure modes like temporarily losing your internet connection.}

\item{verbose}{Show progress bars and other messages indicating current
progress?}
}
\value{
A list with one element for each page of results.
}
\description{
Returns all Direct Message events (both sent and received) within the last 30
days. Sorted in reverse-chronological order. Includes detailed information
about the sender and recipient.
}
\examples{

\dontrun{

## get my direct messages
dms <- direct_messages()

## inspect data structure
str(dms)

}
}
\references{
\url{https://developer.twitter.com/en/docs/twitter-api/v1/direct-messages/sending-and-receiving/api-reference/list-events}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post-list.R
\name{post_list}
\alias{post_list}
\title{Manage Twitter lists}
\usage{
post_list(
  users = NULL,
  name = NULL,
  description = NULL,
  private = FALSE,
  destroy = FALSE,
  list_id = NULL,
  slug = NULL,
  token = NULL
)
}
\arguments{
\item{users}{Character vectors of users to be added to list.}

\item{name}{Name of new list to create.}

\item{description}{Optional, description of list (single character string).}

\item{private}{Logical indicating whether created list should be private.
Defaults to false, meaning the list would be public. Not applicable if list
already exists.}

\item{destroy}{Logical indicating whether to delete a list. Either \code{list_id} or
\code{slug} must be provided if \code{destroy = TRUE}.}

\item{list_id}{Optional, numeric ID of list.}

\item{slug}{Optional, list slug.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}
}
\value{
Response object from HTTP request.
}
\description{
Create, add users, and destroy Twitter lists
}
\examples{
\dontrun{

## R related Twitter accounts
users <- c("_R_Foundation", "R_dev_news", "rweekly_live", "RConsortium", "rstats4ds",
  "icymi_r", "rstatstweet", "RLadiesGlobal")

## create r-accounts list with 8 total users
(r_lst <- post_list(users,
  "r-accounts", description = "R related accounts"))

## view list in browser
browseURL(sprintf("https://twitter.com/\%s/lists/r-accounts",
  rtweet:::api_screen_name()))

## search for more rstats users
r_users <- search_users("rstats", n = 200)

## filter and select more users to add to list
more_users <- r_users$screen_name[r_users$verified]

## add more users to list- note: can only add up to 100 at a time
post_list(users = more_users, slug = "r-accounts")

## view updated list in browser (should be around 100 users)
browseURL(sprintf("https://twitter.com/\%s/lists/r-accounts",
  rtweet:::api_screen_name()))

drop_users <- "icymi_r"

## drop these users from the R list
post_list(users = drop_users, slug = "r-accounts",
  destroy = TRUE)

## view updated list in browser (should be around 100 users)
browseURL(sprintf("https://twitter.com/\%s/lists/r-accounts",
  rtweet:::api_screen_name()))

## delete list entirely
post_list(slug = "r-accounts", destroy = TRUE)

}
}
\references{
Create: \url{https://developer.twitter.com/en/docs/twitter-api/v1/accounts-and-users/create-manage-lists/api-reference/post-lists-create}
Destroy: \url{https://developer.twitter.com/en/docs/twitter-api/v1/accounts-and-users/create-manage-lists/api-reference/post-lists-destroy}
Add users: \url{https://developer.twitter.com/en/docs/twitter-api/v1/accounts-and-users/create-manage-lists/api-reference/post-lists-members-create},
\url{https://developer.twitter.com/en/docs/twitter-api/v1/accounts-and-users/create-manage-lists/api-reference/post-lists-members-create_all}
Remove users: \url{https://developer.twitter.com/en/docs/twitter-api/v1/accounts-and-users/create-manage-lists/api-reference/post-lists-members-destroy},
\url{https://developer.twitter.com/en/docs/twitter-api/v1/accounts-and-users/create-manage-lists/api-reference/post-lists-members-destroy_all}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post_destroy.R
\name{post_destroy}
\alias{post_destroy}
\title{Delete status of user's Twitter account}
\usage{
post_destroy(destroy_id, token = NULL)
}
\arguments{
\item{destroy_id}{To delete a status, supply the single status ID here. If a
character string is supplied, overriding the default (NULL), then a destroy
request is made (and the status text and media attachments) are irrelevant.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}
}
\description{
Deletes a status of user's profile.
}
\examples{
\dontrun{
pt <- post_tweet()
crt <- httr::content(pt)
post_destroy(crt$id_str)
}
}
\references{
\url{https://developer.twitter.com/en/docs/twitter-api/v1/tweets/post-and-engage/api-reference/post-statuses-destroy-id}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_users.R
\name{search_users}
\alias{search_users}
\title{Search for users}
\usage{
search_users(q, n = 100, parse = TRUE, token = NULL, verbose = TRUE)
}
\arguments{
\item{q}{As string providing the search query. Try searching by interest,
full name, company name, or location. Exact match searches are not
supported.}

\item{n}{Desired number of results to return. Results are downloaded
in pages when \code{n} is large; the default value will download a single
page. Set \code{n = Inf} to download as many results as possible.

The Twitter API rate limits the number of requests you can perform
in each 15 minute period. The easiest way to download more than that is
to use \code{retryonratelimit = TRUE}.

You are not guaranteed to get exactly \code{n} results back. You will get
fewer results when tweets have been deleted or if you hit a rate limit.
You will get more results if you ask for a number of tweets that's not
a multiple of page size, e.g. if you request \code{n = 150} and the page
size is 200, you'll get 200 results back.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}

\item{verbose}{Show progress bars and other messages indicating current
progress?}
}
\value{
Data frame with one row for each matching user.
}
\description{
Search for Twitter users. The Twitter API limits the results to at most
1,000 users.
}
\examples{

\dontrun{
users <- search_users("#rstats", n = 300)
users

# latest tweet from each user
tweets_data(users)
}

}
\references{
\url{https://developer.twitter.com/en/docs/twitter-api/v1/accounts-and-users/follow-search-get-users/api-reference/get-users-search}
}
\seealso{
Other users: 
\code{\link{as_screenname}()},
\code{\link{lists_subscribers}()},
\code{\link{lookup_users}()}
}
\concept{users}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph-network.R
\name{network_data}
\alias{network_data}
\alias{network_graph}
\title{Network data}
\usage{
network_data(x, e = c("mention", "retweet", "reply", "quote"))

network_graph(x, e = c("mention", "retweet", "reply", "quote"))
}
\arguments{
\item{x}{Data frame returned by rtweet function}

\item{e}{Type of edge/link–i.e., "mention", "retweet", "quote", "reply".
This must be a character vector of length one or more. This value will be
split on punctuation and space (so you can include multiple types in the
same string separated by a comma or space). The values "all" and
"semantic" are assumed to mean all edge types, which is equivalent to the
default value of \code{c("mention", "retweet", "reply", "quote")}}
}
\value{
A from/to data edge data frame

An igraph object
}
\description{
\itemize{
\item \code{network_data()} returns a data frame that can easily be converted to
various network classes.
\item \code{network_graph()} returns a igraph object
}
}
\details{
Retrieve data to know which users are connected to which users.
}
\examples{

\dontrun{
  ## search for #rstats tweets
  rstats <- search_tweets("#rstats", n = 200)

  ## create from-to data frame representing retweet/mention/reply connections
  rstats_net <- network_data(rstats, c("retweet","mention","reply"))

  ## view edge data frame
  rstats_net

  ## view user_id->screen_name index
  attr(rstats_net, "idsn")

  ## if igraph is installed...
  if (requireNamespace("igraph", quietly = TRUE)) {

    ## (1) convert directly to graph object representing semantic network
    rstats_net <- network_graph(rstats)

    ## (2) plot graph via igraph.plotting
    plot(rstats_net)
  }
}
}
\seealso{
network_graph
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bearer_token.R
\name{bearer_token}
\alias{bearer_token}
\title{Bearer token}
\usage{
bearer_token(token = NULL)
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}

\code{bearer_token()} has been deprecated because it used a rather indirect
method of generating a bearer token. Instead, use \code{\link[=rtweet_app]{rtweet_app()}}, copying
in the bearer token directly from the
\href{https://developer.twitter.com/en/portal/projects-and-apps}{Twitter developer portal}.
See \code{vignette("auth")} for full details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/save_as_csv.R
\name{read_twitter_csv}
\alias{read_twitter_csv}
\title{Read comma separated value Twitter data.}
\usage{
read_twitter_csv(file, unflatten = FALSE)
}
\arguments{
\item{file}{Name of CSV file.}

\item{unflatten}{Logical indicating whether to unflatten (separate hasthags
and mentions columns on space, converting characters to lists), defaults
to FALSE.}
}
\value{
A tbl data frame of Twitter data
}
\description{
Reads Twitter data that was previously saved as a CSV file.
}
\examples{

\dontrun{

## read in data.csv
rt <- read_twitter_csv("data.csv")

}
}
\seealso{
Other datafiles: 
\code{\link{flatten}()},
\code{\link{write_as_csv}()}
}
\concept{datafiles}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractors.R
\name{users_data}
\alias{users_data}
\alias{tweets_data}
\title{Get tweets from users, or users from tweets}
\usage{
users_data(tweets)

tweets_data(users)
}
\arguments{
\item{tweets}{A data frame of tweets.}

\item{users}{A data frame of users.}
}
\value{
\code{user_data()} returns a data frame of users; \code{tweets_data()}
returns a data frame of tweets.
}
\description{
Twitter API endpoints that return tweets also return data about the users who
tweeted, and most endpoints that return users also return their last tweet.
Showing these additional columns would clutter the default display, so
rtweet instead stores in special attributes and allows you to retrieve them
with the \code{user_data()} and \code{tweets_data()} helpers.
}
\examples{
\dontrun{
# find users from tweets
tweets <- search_tweets("r")
users_data(tweets)

# from tweets from users
users <- search_users("r")
tweets_data(users)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auth.R
\name{auth_setup_default}
\alias{auth_setup_default}
\title{Set up default authentication}
\usage{
auth_setup_default()
}
\description{
You'll need to run this function once per computer so that rtweet can use
your personal Twitter account. See \code{\link[=rtweet_app]{rtweet_app()}}/\link{rtweet_bot} and
\code{\link[=auth_save]{auth_save()}} for other authentication options.
}
\seealso{
Other authentication: 
\code{\link{auth_as}()},
\code{\link{auth_get}()},
\code{\link{auth_save}()},
\code{\link{rtweet_user}()}
}
\concept{authentication}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coords.R
\name{lookup_coords}
\alias{lookup_coords}
\title{Get coordinates of specified location.}
\usage{
lookup_coords(address, components = NULL, apikey = NULL, ...)
}
\arguments{
\item{address}{Desired location typically in the form of place
name, subregion, e.g., address = "lawrence, KS". Also accepts the
name of countries, e.g., address = "usa", address = "brazil" or
states, e.g., address = "missouri" or cities, e.g., address =
"chicago". In most cases using only address should be sufficient.}

\item{components}{Unit of analysis for address e.g., components =
"country:US". Potential components include postal_code, country,
administrative_area, locality, route.}

\item{apikey}{A valid Google Maps API key. If NULL, \code{lookup_coords()} will
look for a relevant API key stored as an environment variable (e.g.,
\code{GOOGLE_MAPS_KEY}).}

\item{...}{Additional arguments passed as parameters in the HTTP
request}
}
\value{
Object of class coords.
}
\description{
Convenience function for looking up latitude/longitude coordinate
information for a given location. Returns data as a special
"coords" object, which is specifically designed to interact
smoothly with other relevant package functions. NOTE: USE OF THIS FUNCTION
REQUIRES A VALID GOOGLE MAPS API KEY.
}
\details{
Since Google Maps implemented stricter API requirements, sending
requests to Google's API isn't very convenient. To enable basic uses
without requiring a Google Maps API key, a number of the major cities
throughout the word and the following two larger locations are
baked into this function: 'world' and 'usa'. If 'world' is supplied then
a bounding box of maximum latitutde/longitude values, i.e.,
\code{c(-180, -90, 180, 90)}, and a center point \code{c(0, 0)} are
returned. If 'usa' is supplied then estimates of the United States'
bounding box and mid-point are returned. To specify a city, provide the
city name followed by a space and then the US state abbreviation or
country name. To see a list of all included cities, enter
\code{rtweet:::citycoords} in the R console to see coordinates data.
}
\examples{

\dontrun{

## get coordinates associated with the following addresses/components
sf <- lookup_coords("san francisco, CA", "country:US")
usa <- lookup_coords("usa")
lnd <- lookup_coords("london")
bz <- lookup_coords("brazil")

## pass a returned coords object to search_tweets
bztw <- search_tweets(geocode = bz)

## or stream tweets
ustw <- stream_tweets(usa, timeout = 10)

}

}
\seealso{
Other geo: 
\code{\link{lat_lng}()}
}
\concept{geo}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/do_call_rbind.R
\name{do_call_rbind}
\alias{do_call_rbind}
\title{Binds list of data frames while preserving attribute (tweets or users) data.}
\usage{
do_call_rbind(x)
}
\arguments{
\item{x}{List of parsed tweets data or users data, each of which
presumably contains an attribute of the other (i.e., users data
contains tweets attribute; tweets data contains users attribute).}
}
\value{
A single merged (by row) data frame (tbl) of tweets or
users data that also contains as an attribute a merged (by row)
data frame (tbl) of its counterpart, making it accessible via the
\code{\link[=users_data]{users_data()}} or \code{\link[=tweets_data]{tweets_data()}} extractor
functions.
}
\description{
Row bind lists of tweets/users data whilst also preserving and binding
users/tweets attribute data.
}
\examples{

\dontrun{

## lapply through three different search queries
lrt <- lapply(
  c("rstats OR tidyverse", "data science", "python"),
  search_tweets,
  n = 5000
)

## convert list object into single parsed data rame
rt <- do_call_rbind(lrt)

## preview tweets data
rt

## preview users data
users_data(rt)

}

}
\concept{parsing}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user_id.R
\name{as_screenname}
\alias{as_screenname}
\alias{as_userid}
\title{Mark a user id as a screen name}
\usage{
as_screenname(x)
}
\arguments{
\item{x}{A character vector of Twitter screen names.}
}
\description{
There are two ways to identify a Twitter user: a screen name (e.g.
"justinbieber") or a user identifier (e.g. "27260086"). User identifiers
look like regular numbers, but are actually 64-bit integers which can't be
stored in R's numeric vectors. For this reason, rtweet always returns ids as
strings.

Unfortunately this introduces an ambiguity, because user names can
also consist solely of numbers (e.g. "123456") so it's not obvious whether
a string consisting only of numbers is a screen name or a user id. By
default, rtweet will assume its a user id, so if you have a screen name
that consists only of numbers, you'll need to use \code{as_screenname()} to
tell rtweet that it's actually a screen name.

Note that in general, you are best off using user ids; screen names are
not static and may change over longer periods of time.
}
\examples{
\dontrun{
# Look up user with id 123456 (screen name harperreed)
lookup_users("123456") 

# Look up user with name 123456
lookup_users(as_screenname("123456"))
}
}
\seealso{
Other users: 
\code{\link{lists_subscribers}()},
\code{\link{lookup_users}()},
\code{\link{search_users}()}
}
\concept{users}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auth.R
\name{rtweet_user}
\alias{rtweet_user}
\alias{rtweet_bot}
\alias{rtweet_app}
\title{Authentication options}
\usage{
rtweet_user(api_key = NULL, api_secret = NULL)

rtweet_bot(
  api_key = ask_pass("API key"),
  api_secret = ask_pass("API secret"),
  access_token = ask_pass("access token"),
  access_secret = ask_pass("access token")
)

rtweet_app(bearer_token = ask_pass("bearer token"))
}
\arguments{
\item{api_key, api_secret}{Application API key and secret. These are
generally not required for \code{tweet_user()} since the defaults will use
the built-in rtweet app.}

\item{access_token, access_secret}{Access token and secret.}

\item{bearer_token}{App bearer token.}
}
\description{
There are three ways that you can authenticate with the Twitter API:
\itemize{
\item \code{rtweet_user()} interactively authenticates an existing Twitter user.
This form is most appropriate if you want rtweet to control your
Twitter account.
\item \code{rtweet_app()} authenticates as a Twitter application. An application can't
perform actions (i.e. it can't tweet) but otherwise has generally higher
rate limits (i.e. you can do more searches). See details
at \url{https://developer.twitter.com/en/docs/basics/rate-limits.html}.
This form is most appropriate if you are collecting data.
\item \code{rtweet_bot()} authenticates as bot that takes actions on behalf of an app.
This form is most appropriate if you want to create a Twitter account that
is run by a computer, rather than a human.
}

To use \code{rtweet_app()} or \code{rtweet_bot()} you will need to create your own
Twitter app following the instructions in \code{vignette("auth.Rmd")}.
\code{rtweet_user()} \emph{can be} used with your own app, but generally there is
no need to because it uses the Twitter app provided by rtweet.

Use \code{\link[=auth_as]{auth_as()}} to set the default auth mechanism for the current session,
and \code{\link[=auth_save]{auth_save()}} to save an auth mechanism for use in future sessions.
}
\section{Security}{
All of the arguments to these functions are roughly equivalent to
passwords so should generally not be typed into the console (where they
the will be recorded in \code{.Rhistory}) or recorded in a script (which is
easy to accidentally share). Instead, call these functions without arguments
since the default behaviour is to use \code{\link[askpass:askpass]{askpass::askpass()}} to interactively
prompt you for the values.
}

\seealso{
Other authentication: 
\code{\link{auth_as}()},
\code{\link{auth_get}()},
\code{\link{auth_save}()},
\code{\link{auth_setup_default}()}
}
\concept{authentication}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lists_statuses.R
\name{lists_statuses}
\alias{lists_statuses}
\title{Get a timeline of tweets authored by members of a specified list.}
\usage{
lists_statuses(
  list_id = NULL,
  slug = NULL,
  owner_user = NULL,
  since_id = NULL,
  max_id = NULL,
  n = 200,
  include_rts = TRUE,
  parse = TRUE,
  retryonratelimit = NULL,
  verbose = TRUE,
  token = NULL
)
}
\arguments{
\item{list_id}{required The numerical id of the list.}

\item{slug}{required You can identify a list by its slug instead of
its numerical id. If you decide to do so, note that you'll also have
to specify the list owner using the owner_id or owner_screen_name
parameters.}

\item{owner_user}{optional The screen name or user ID of the user
who owns the list being requested by a slug.}

\item{since_id}{Supply a vector of ids or a data frame of previous results to
find tweets \strong{newer} than \code{since_id}.}

\item{max_id}{Supply a vector of ids or a data frame of previous results to
find tweets \strong{older} than \code{max_id}.}

\item{n}{Desired number of results to return. Results are downloaded
in pages when \code{n} is large; the default value will download a single
page. Set \code{n = Inf} to download as many results as possible.

The Twitter API rate limits the number of requests you can perform
in each 15 minute period. The easiest way to download more than that is
to use \code{retryonratelimit = TRUE}.

You are not guaranteed to get exactly \code{n} results back. You will get
fewer results when tweets have been deleted or if you hit a rate limit.
You will get more results if you ask for a number of tweets that's not
a multiple of page size, e.g. if you request \code{n = 150} and the page
size is 200, you'll get 200 results back.}

\item{include_rts}{optional When set to either true, t or 1,
the list timeline will contain native retweets (if they exist) in
addition to the standard stream of tweets. The output format of
retweeted tweets is identical to the representation you see in
home_timeline.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}

\item{retryonratelimit}{If \code{TRUE}, and a rate limit is exhausted, will wait
until it refreshes. Most Twitter rate limits refresh every 15 minutes.
If \code{FALSE}, and the rate limit is exceeded, the function will terminate
early with a warning; you'll still get back all results received up to
that point. The default value, \code{NULL}, consults the option
\code{rtweet.retryonratelimit} so that you can globally set it to \code{TRUE},
if desired.

If you expect a query to take hours or days to perform, you should not
rely soley on \code{retryonratelimit} because it does not handle other common
failure modes like temporarily losing your internet connection.}

\item{verbose}{Show progress bars and other messages indicating current
progress?}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}
}
\value{
data
}
\description{
Get a timeline of tweets authored by members of a specified list.
}
\examples{
\dontrun{
(ls_senators <- lists_statuses(slug = "senators", owner_user = "cspan"))
}
}
\seealso{
Other lists: 
\code{\link{lists_members}()},
\code{\link{lists_subscribers}()},
\code{\link{lists_subscriptions}()},
\code{\link{lists_users}()}

Other tweets: 
\code{\link{get_favorites}()},
\code{\link{get_mentions}()},
\code{\link{get_timeline}()},
\code{\link{lookup_tweets}()},
\code{\link{search_tweets}()}
}
\concept{lists}
\concept{tweets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stream.R
\name{stream_tweets2}
\alias{stream_tweets2}
\title{A more robust version of stream_tweets}
\usage{
stream_tweets2(..., dir = NULL, append = FALSE)
}
\arguments{
\item{dir}{Name of directory in which json files should be written.
The default, NULL, will create a timestamped "stream" folder in the
current working directory. If a dir name is provided that does not
already exist, one will be created.}

\item{append}{Logical indicating whether to append or overwrite
file_name if the file already exists. Defaults to FALSE, meaning
this function will overwrite the preexisting file_name (in other
words, it will delete any old file with the same name as
file_name) meaning the data will be added as new lines to file if
pre-existing.}
}
\value{
Returns data as expected using original search_tweets
function.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
Please use \code{\link[=stream_tweets]{stream_tweets()}} instead.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ts_plot.R
\name{round_time}
\alias{round_time}
\title{A generic function for rounding date and time values}
\usage{
round_time(x, n, tz)
}
\arguments{
\item{x}{A vector of class POSIX or Date.}

\item{n}{Unit to round to. Defaults to mins. Numeric values treated
as seconds. Otherwise this should be one of "mins", "hours", "days",
"weeks", "months", "years" (plural optional).}

\item{tz}{Time zone to be used, defaults to "UTC" (Twitter default)}
}
\value{
If POSIXct then POSIX. If date then Date.
}
\description{
A generic function for rounding date and time values
}
\examples{

## class posixct
round_time(Sys.time(), "12 hours")

## class date
unique(round_time(seq(Sys.Date(), Sys.Date() + 100, "1 day"), "weeks"))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tweet_threading.R
\name{tweet_threading}
\alias{tweet_threading}
\title{Collect statuses contained in a thread}
\usage{
tweet_threading(
  tw,
  traverse = c("backwards", "forwards"),
  n = 10,
  verbose = FALSE
)
}
\arguments{
\item{tw}{\code{\link[=lookup_tweets]{lookup_tweets()}} output containing
at least the last status in the thread or an id of a tweet.}

\item{traverse}{character, direction to traverse from origin status in tw,
Default: c('backwards','forwards')}

\item{n}{numeric, timeline to fetch to start forwards traversing, Default: 10}

\item{verbose}{logical, Output to console status of traverse, Default: FALSE}
}
\value{
\code{\link[=lookup_tweets]{lookup_tweets()}} tibble
}
\description{
Return all statuses that are part of a thread (Replies from a user to their
own tweets). By default the function traverses first backwards from the
origin status_id of the thread up to the root, then checks if there are any
child statuses that were posted after the origin status.
}
\examples{
\dontrun{
tw_thread <- tweet_threading("1461776330584956929")
tw_thread
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collections.R
\name{lookup_collections}
\alias{lookup_collections}
\alias{get_collections}
\title{Collections API}
\usage{
lookup_collections(id, n = 200, parse = TRUE, token = NULL, ...)

get_collections(
  user = NULL,
  status_id = NULL,
  n = 200,
  cursor = NULL,
  parse = TRUE,
  token = NULL
)
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}

\code{get_collections()} and \code{lookup_collections()} have been deprecated
since the underlying Twitter API has been deprecated:
\url{https://developer.twitter.com/en/docs/twitter-for-websites/timelines/guides/collection}.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statuses.R
\name{lookup_tweets}
\alias{lookup_tweets}
\alias{lookup_statuses}
\title{Get tweets data for given statuses (status IDs).}
\usage{
lookup_tweets(
  statuses,
  parse = TRUE,
  token = NULL,
  retryonratelimit = NULL,
  verbose = TRUE
)
}
\arguments{
\item{statuses}{User id or screen name of target user.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}

\item{retryonratelimit}{If \code{TRUE}, and a rate limit is exhausted, will wait
until it refreshes. Most Twitter rate limits refresh every 15 minutes.
If \code{FALSE}, and the rate limit is exceeded, the function will terminate
early with a warning; you'll still get back all results received up to
that point. The default value, \code{NULL}, consults the option
\code{rtweet.retryonratelimit} so that you can globally set it to \code{TRUE},
if desired.

If you expect a query to take hours or days to perform, you should not
rely soley on \code{retryonratelimit} because it does not handle other common
failure modes like temporarily losing your internet connection.}

\item{verbose}{Show progress bars and other messages indicating current
progress?}
}
\value{
A tibble of tweets data.
}
\description{
Get tweets data for given statuses (status IDs).
}
\examples{

\dontrun{
statuses <- c(
  "567053242429734913",
  "266031293945503744",
  "440322224407314432"
)

## lookup tweets data for given statuses
tw <- lookup_tweets(statuses)
tw
}
}
\references{
\url{https://developer.twitter.com/en/docs/tweets/post-and-engage/api-reference/get-statuses-lookup}
}
\seealso{
Other tweets: 
\code{\link{get_favorites}()},
\code{\link{get_mentions}()},
\code{\link{get_timeline}()},
\code{\link{lists_statuses}()},
\code{\link{search_tweets}()}
}
\concept{tweets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/premium.R
\name{search_fullarchive}
\alias{search_fullarchive}
\alias{search_30day}
\title{Premium Twitter searches}
\usage{
search_fullarchive(
  q,
  n = 100,
  fromDate = NULL,
  toDate = NULL,
  env_name = NULL,
  safedir = NULL,
  parse = TRUE,
  token = NULL
)

search_30day(
  q,
  n = 100,
  fromDate = NULL,
  toDate = NULL,
  env_name = NULL,
  safedir = NULL,
  parse = TRUE,
  token = NULL
)
}
\arguments{
\item{q}{Search query on which to match/filter tweets. See details for
information about available search operators.}

\item{n}{Desired number of results to return. Results are downloaded
in pages when \code{n} is large; the default value will download a single
page. Set \code{n = Inf} to download as many results as possible.

The Twitter API rate limits the number of requests you can perform
in each 15 minute period. The easiest way to download more than that is
to use \code{retryonratelimit = TRUE}.

You are not guaranteed to get exactly \code{n} results back. You will get
fewer results when tweets have been deleted or if you hit a rate limit.
You will get more results if you ask for a number of tweets that's not
a multiple of page size, e.g. if you request \code{n = 150} and the page
size is 200, you'll get 200 results back.}

\item{fromDate}{Oldest date-time (YYYYMMDDHHMM) from which tweets should be
searched for.}

\item{toDate}{Newest date-time (YYYYMMDDHHMM) from which tweets should be
searched for.}

\item{env_name}{Name/label of developer environment to use for the search.}

\item{safedir}{Name of directory to which each response object should be
saved. If the directory doesn't exist, it will be created. If NULL (the
default) then a dir will be created in the current working directory. To
override/deactivate safedir set this to FALSE.}

\item{parse}{If \code{TRUE}, the default, returns a tidy data frame. Use \code{FALSE}
to return the "raw" list corresponding to the JSON returned from the
Twitter API.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}
}
\value{
A tibble data frame of Twitter data
}
\description{
Search 30day or fullarchive products
}
\section{Developer Account}{

Users must have an approved developer account and an active/labeled
environment to access Twitter's premium APIs. For more information, to check
your current Subscriptions and Dev Environments, or to apply for a developer
account visit \url{https://developer.twitter.com}.
}

\section{Search operators}{

\emph{Note: Bolded operators ending with a colon should be immediately
followed by a word or quoted phrase (if appropriate)–e.g.,} \code{lang:en}
}

\section{Keyword}{

\itemize{
\item \strong{""}           ~~ match exact phrase
\item \strong{#}               ~~ hashtag
\item \strong{@}               ~~ at mentions)
\item \strong{url:}            ~~ found in URL
\item \strong{lang:}           ~~ language of tweet
}
}

\section{Accounts of interest}{

\itemize{
\item \strong{from:}           ~~ authored by
\item \strong{to:}             ~~ sent to
\item \strong{retweets_of:}    ~~ retweet author
}
}

\section{Tweet attributes}{

\itemize{
\item \strong{is:retweet}      ~~ only retweets
\item \strong{has:mentions}    ~~ uses mention(s)
\item \strong{has:hashtags}    ~~ uses hashtags(s)
\item \strong{has:media}       ~~ includes media(s)
\item \strong{has:videos}      ~~ includes video(s)
\item \strong{has:images}      ~~ includes image(s)
\item \strong{has:links}       ~~ includes URL(s)
\item \strong{is:verified}     ~~ from verified accounts
}
}

\section{Geospatial}{

\itemize{
\item \strong{bounding_box:[west_long south_lat east_long north_lat]} ~~ lat/long coordinates box
\item \strong{point_radius:[lon lat radius]} ~~ center of search radius
\item \strong{has:geo}           ~~ uses geotagging
\item \strong{place:}            ~~ by place
\item \strong{place_country:}    ~~ by country
\item \strong{has:profile_geo}   ~~ geo associated with profile
\item \strong{profile_country:}  ~~ country associated with profile
\item \strong{profile_region:}   ~~ region associated with profile
\item \strong{profile_locality:} ~~ locality associated with profile
}
}

\examples{

\dontrun{
## search fullarchive for up to 300 rstats tweets sent in Jan 2014
rt <- search_fullarchive("#rstats", n = 300, env_name = "research",
  fromDate = "201401010000", toDate = "201401312359")

toDate <- format(Sys.time() - 60 * 60 * 24 * 7, "\%Y\%m\%d\%H\%M")

## search 30day for up to 300 rstats tweets sent before the last week
rt <- search_30day("#rstats", n = 300,
  env_name = "research", toDate = toDate)
}

}
\references{
\url{https://developer.twitter.com/en/docs/twitter-api/premium/search-api/api-reference/premium-search}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tokens.R
\name{get_token}
\alias{get_token}
\alias{get_tokens}
\title{Fetch Twitter OAuth token}
\usage{
get_token()

get_tokens()
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}
Please use \code{\link[=auth_get]{auth_get()}} instead.
}
\seealso{
Other tokens: 
\code{\link{create_token}()},
\code{\link{rate_limit}()}
}
\concept{tokens}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plain_tweets.R
\name{plain_tweets}
\alias{plain_tweets}
\title{Clean up character vector (tweets) to more of a plain text.}
\usage{
plain_tweets(x)
}
\arguments{
\item{x}{The desired character vector or data frame/list with named column/element
"text" to be cleaned and processed.}
}
\value{
Data reformatted with ascii encoding and normal ampersands and
without URL links, line breaks, fancy spaces/tabs, fancy apostrophes,
}
\description{
Removes links, linebreaks, fancy spaces and apostrophes and convert everything to ASCII text.
Deprecated to be defunct for next release as there are better text processing tools.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tweets_and_users.R
\name{tweets_with_users}
\alias{tweets_with_users}
\alias{users_with_tweets}
\title{Parsing data into tweets/users data tibbles}
\usage{
tweets_with_users(x)

users_with_tweets(x)
}
\arguments{
\item{x}{A list of responses, with one element for each page.}
}
\value{
A tweets/users tibble with users/tweets attribute.
}
\description{
For internal use only
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post-user.R
\name{post_follow}
\alias{post_follow}
\alias{follow_user}
\alias{post_unfollow_user}
\alias{unfollow_user}
\alias{post_mute}
\alias{mute_user}
\title{Follows target Twitter user.}
\usage{
post_follow(
  user,
  destroy = FALSE,
  mute = FALSE,
  notify = FALSE,
  retweets = TRUE,
  token = NULL
)

post_unfollow_user(user, token = NULL)

post_mute(user, token = NULL)
}
\arguments{
\item{user}{Character vector of screen names or user ids.
See \code{\link[=as_screenname]{as_screenname()}} for more details.}

\item{destroy}{Logical indicating whether to post (add) or
remove (delete) target tweet as favorite.}

\item{mute}{Logical indicating whether to mute the intended
friend (you must already be following this account prior
to muting them)}

\item{notify}{Logical indicating whether to enable notifications
for target user. Defaults to false.}

\item{retweets}{Logical indicating whether to enable retweets
for target user. Defaults to true.}

\item{token}{Expert use only. Use this to override authentication for
a single API call. In most cases you are better off changing the
default for all calls. See \code{\link[=auth_as]{auth_as()}} for details.}
}
\description{
Follows target Twitter user.
}
\examples{
\dontrun{
post_follow("BarackObama")
}
}
\references{
Update: \url{https://developer.twitter.com/en/docs/twitter-api/v1/accounts-and-users/follow-search-get-users/api-reference/post-friendships-update}
Create: \url{https://developer.twitter.com/en/docs/twitter-api/v1/accounts-and-users/follow-search-get-users/api-reference/post-friendships-create}
Destroy: \url{https://developer.twitter.com/en/docs/twitter-api/v1/accounts-and-users/follow-search-get-users/api-reference/post-friendships-destroy}
Mute: \url{https://developer.twitter.com/en/docs/twitter-api/v1/accounts-and-users/mute-block-report-users/api-reference/post-mutes-users-create}
}
\seealso{
Other post: 
\code{\link{post_favorite}()},
\code{\link{post_friendship}()},
\code{\link{post_tweet}()}
}
\concept{post}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/next_cursor.R
\name{next_cursor}
\alias{next_cursor}
\alias{since_id}
\title{Extract cursor (for cursor based pagination)}
\usage{
next_cursor(x)

since_id(x)
}
\description{
This internal helper extracts the cursor from the object passed to the
\code{cursor} argument of the functions that use \code{\link[=TWIT_paginate_cursor]{TWIT_paginate_cursor()}}.
}
\examples{
\dontrun{
page1 <- get_followers("potus")
page2 <- get_followers("potus", cursor = page1)
}
}
\keyword{internal}
