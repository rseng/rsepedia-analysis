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
