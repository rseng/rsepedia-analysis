# roreviewapi

<!-- badges: start -->

[![R build
status](https://github.com/ropensci-review-tools/roreviewapi/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci-review-tools/roreviewapi/actions?query=workflow%3AR-CMD-check)
[![Project Status:
Concept](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
<!-- badges: end -->

Plumber API to generate reports on package structure and function for
the [`@ropensci-review-bot`](https://github.com/ropensci-review-bot).
The package is not intended for general use, and these documents are
primarily intended for maintainers of this package, although they may
serve as useful templates for similar endeavours. Please contact us if
you have any questions.

Uses functionality provided by the
[`pkgcheck`](https://github.com/ropensci-review-tools/pkgcheck) and
[`pkgstats`](https://github.com/ropensci-review-tools/pkgstats)
packages. This suite of three packages requires a few system installs,
two for [`pkgstats`](https://github.com/ropensci-review-tools/pkgstats)
of [`ctags`](https://ctags.io) and [GNU
`global`](https://www.gnu.org/software/global/). Procedures to install
these libraries on various operating systems are described in the
[`pkgstats` package](https://docs.ropensci.org/pkgstats). This package
itself also requires both the [GitHub command-line-interface (cli),
`gh`](https://cli.github.com/) and
[`dos2unix`](https://sourceforge.net/projects/dos2unix/).

A local GitHub token also needs to be stored as an environment variable
named `GITHUB_TOKEN`, and not `GITHUB_PAT` or anything else; the `gh`
cli only recognises the former name.

The package also works by locally caching previously analysed packages,
in a `pkgcheck` subdirectory of the location determined by

``` r
rappdirs::user_cache_dir()
```

You may manually erase the contents of this subdirectory at any time at
no risk.

## Dockerfile

The server associated with this package can be built by cloning this
repository, and modifying the associated
[`Dockerfile`](https://github.com/ropensci-review-tools/roreviewapi/blob/master/Dockerfile)
by inserting a GitHub token, and associated `git config` options.

## Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.
---
title: "roreviewapi"
output:
  md_document:
    variant: gfm

  rmarkdown::html_vignette:
    self_contained: no
---

# roreviewapi

<!-- badges: start -->
[![R build status](https://github.com/ropensci-review-tools/roreviewapi/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci-review-tools/roreviewapi/actions?query=workflow%3AR-CMD-check)
[![Project Status: Concept](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
<!-- badges: end -->

Plumber API to generate reports on package structure and function for the
[`@ropensci-review-bot`](https://github.com/ropensci-review-bot). The package
is not intended for general use, and these documents are primarily intended for
maintainers of this package, although they may serve as useful templates for
similar endeavours. Please contact us if you have any questions.

Uses functionality provided by the
[`pkgcheck`](https://github.com/ropensci-review-tools/pkgcheck) and
[`pkgstats`](https://github.com/ropensci-review-tools/pkgstats) packages. This suite of
three packages requires a few system installs, two for 
[`pkgstats`](https://github.com/ropensci-review-tools/pkgstats) of
[`ctags`](https://ctags.io) and [GNU
`global`](https://www.gnu.org/software/global/). Procedures to install these
libraries on various operating systems are described in the [`pkgstats`
package](https://docs.ropensci.org/pkgstats). This package itself also requires
both the [GitHub command-line-interface (cli), `gh`](https://cli.github.com/)
and [`dos2unix`](https://sourceforge.net/projects/dos2unix/).

A local GitHub token also needs to be stored as an environment
variable named `GITHUB_TOKEN`, and not `GITHUB_PAT` or anything else; the `gh`
cli only recognises the former name.

The package also works by locally caching previously analysed packages, in
a `pkgcheck` subdirectory of the location determined by
```{r, rappdirs, eval = FALSE}
rappdirs::user_cache_dir()
```
You may manually erase the contents of this subdirectory at any time at no
risk.


## Dockerfile

The server associated with this package can be built by cloning this
repository, and modifying the associated
[`Dockerfile`](https://github.com/ropensci-review-tools/roreviewapi/blob/master/Dockerfile)
by inserting a GitHub token, and associated `git config` options.

## Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.
---
title: "Debugging the package check API"
author: 
  - "Mark Padgham"
date: "`r Sys.Date()`"
vignette: >
  %\VignetteIndexEntry{How to debug the package check API}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = TRUE,
  message = TRUE,
  width = 120,
  comment = "#>",
  fig.retina = 2,
  fig.path = "README-"
)
```

The main API endpoint for the Digital Ocean server is the
[`editorcheck`](https://docs.ropensci.org/roreviewapi/articles/endpoints.html#editorcheck),
called by the `ropensci-review-bot` on every submission, and also manually
(re-)triggered by the command `@ropensci-review-bot check package`. This
vignette describes procedures which may be used to debug any instances in which
package checks fail.

Most failures at the server end will deliver HTTP error codes of "500" [like
this
example](https://github.com/ropensci/software-review/issues/463#issuecomment-920941332).
This code *may* indicate that the server is currently unavailable, which will
arise once per week during rebuild and redeploy processes, scheduled every
Monday at 23:59 UTC. If a submission happens to be roughly at this time, the
recommended procedure is to wait at least an hour before manually trying
`@ropensci-review-bot check package`. Other 500 codes should be reported
straight to maintainers of the check system, currently to
[\@mpadge](https://github.com/mpadge) and
[\@noamross](https://github.com/noamross) as respective first and second points
of contact. This vignette describes the procedures they would follow to debug
any problems. 

Debugging generally requires stepping through the code as called by the main
[`editorcheck`
endpoint](https://github.com/ropensci-review-tools/roreviewapi/blob/main/R/plumber.R),
and into the [`editor_check()`
function](https://github.com/ropensci-review-tools/roreviewapi/blob/main/R/editor-check.R)
called by the endpoint as a background process.

# Debugging Procedure

## Check that the API is online

The following code will confirm that the API is online, by returning a single numeric value:

```{r do-check, eval = FALSE}
httr::GET ("http://<ip-address>:8000/mean") |>
    httr::content ()
```

## Check log of recent requests

The log entries are illustrated in the [endpoints
vignette](https://docs.ropensci.org/roreviewapi/articles/endpoints.html#5-log),
and can be used to ensure that all variables were successfully delivered in the
request. Any missing or malformed variables most likely indicate problems with
the submission template. These should be caught by the initial call make by the
[`editorcheck` function to the `check_issue_template()`
function](https://github.com/ropensci-review-tools/roreviewapi/blob/82c9724a712094e4ccabb3974cb952f68ad5180f/R/plumber.R#L23). Potential issues can be debugged by calling that function locally:

```{r check-template, eval = FALSE}
roreviewapi::check_issue_template (orgrepo = "ropensci/software-review",
                                   issue_num = <issue_num>)
```

That should return an empty string with an additional attribute,
`"proceed_with_checks"`, which should be `TRUE`. Any other return value
indicates an issue with the submission template, for which the return value
should contain text describing the problem.

## Check installation of system dependencies

This and the following checks are components of the [`editor_check()`
function](https://github.com/ropensci-review-tools/roreviewapi/blob/main/R/editor-check.R)
called by the main endpoint as a background process, the first step of which is
to identify and install system dependencies. This step is error prone, as is
also the case on all CRAN machines. A final step of the check is to identify
any packages which were unable to be installed, and to post a list of these
directly in the submission issue. Such instances should generally be temporary,
and fixed by forthcoming CRAN updates. These happen because of temporary
breakages in one package which lead to other packages becoming unable to be
installed from CRAN, with these temporary breakages in turn generally arising
because of changes to external (non-R) system libraries. If any notified
inabilities to install packages adversely affect final check results, debugging
may require manually running checks on the Digital Ocean server, as described
in the final section below.



## Check system output and error logs

The system output and error logs from checks for a specified package are
accessible from [the `stdlogs`
endpoint](https://docs.ropensci.org/roreviewapi/articles/endpoints.html#7-stdlogs),
which requires the single parameter of the `repourl` of the repository being
checked. Logs are kept for the latest git head, and are accessible if and only
if the latest GitHub commit matches the points at which the logs were
generated. In those cases, the return value may offer useful diagnostic
insights into potential problems either with submitted packages, or the check
system itself.

# Manually running checks

If all else fails, checks can be manually run directly from the Digital Ocean
server, and sent straight to the relevant issue. The following procedure
describe how, generally following the main [`editor_check()`
function](https://github.com/ropensci-review-tools/roreviewapi/blob/main/R/editor-check.R).

1. Enter the `roreviewapi` Docker image via
   `docker -it --rm --entrypoint /bin/bash roreviewapi` and start R.
2. Set `repourl <- <url>` and run `path <- roreviewapi::dl_gh_repo (repourl)`
   to download a clone of the repository.
3. Type `os <- "ubuntu"` and `os_release <- "20.04"`.
4. Run `p <- roreviewapi::pkgrep_install_deps (path, os, os_release)`. The
   value, `p`, will list any packages which were unable to be installed.
   These will then need to be manually installed, generally through finding the
   remote/dev URLs for the packages, and running `remotes::install_github()` or
   similar. Note that successful installation may only be possible in a
   particular order, and in the worst cases may be a process of trial and
   error.
5. Finally generate the main checks by running
   `checks <- pkgcheck::pkgcheck(path)`, during which diagnostic output will be
   dumped directly to the screen.
6. Convert checks to markdown format by running
   `out <- roreviewapi::collate_editor_check (checks)`.
7. Finally, post the checks to the desired issue with
   `out <- roreviewapi::post_to_issue (out, orgrepo, issue_num)`. Alternative
   values for `orgrepo` and `issue_num` can be used to first confirm that the
   checks look okay before posting them to `ropensci/software-review`.
---
title: "API Endpoints"
author: 
  - "Mark Padgham"
date: "`r Sys.Date()`"
vignette: >
  %\VignetteIndexEntry{How to extend checks}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = TRUE,
  message = TRUE,
  width = 120,
  comment = "#>",
  fig.retina = 2,
  fig.path = "README-"
)
```

This vignette provides brief summaries of the endpoints of the `roreviewapi`
package. These are encoded within the identical
[`R/plumber.R`](https://github.com/ropensci-review-tools/roreviewapi/blob/main/R/plumber.R)
and
[`inst/plumber.R`](https://github.com/ropensci-review-tools/roreviewapi/blob/main/inst/plumber.R)
files. All endpoints respond only to `GET` calls.

---

## 1. editorcheck

This is the main endpoint called by the `ropensci-review-bot` in response to
package submission. The call itself is configured as part of an [external
service call in
`ropensci-org/buffy](https://github.com/ropensci-org/buffy/blob/82dd29bae4aeaa6bf5ca77b27be82cacd3a1ba04/config/settings-production.yml#L18-L32), which passes the parameters specified there of:

1. `repo` The GitHub repository from where the call originates, generally
   `ropensci/software-review`;
2. `issue_id` as the number of the issue in `repo` describing the software
   submission; and
3. `repourl` as specified in the submission template, and specifying the GitHub
   repository of the software being submitted, also in the format
   `<org>/<repo>`.

This endpoint implements the following steps:

1. Call
   [`roreviewapi::check_issue_template()`](https://docs.ropensci.org/roreviewapi/reference/check_issue_template.html)
   to check the existence and format of HTML variables included within the
   submission template. This function returns an empty string if the template
   is okay; otherwise a descriptive error message. The return value also
   includes a binary attribute, `"proceed_with_checks"`, which is set to
   `FALSE` only if `repourl` is improperly specified. In this case the function
   returns immediately with a text string describing the error. Otherwise the string is carried through to the next step:
2. The [`pkgcheck::pkgcheck()`
   function](https://github.com/ropensci-review-tools/pkgcheck) is started as
   a background process, dumping both `stdout` and `stderr` messages to
   specified logfiles (see `stdlogs` endpoint, below).
3. Any messages generated above are prepended to a return message that the
   package checks have started, that message delivered back to the bot, and
   ultimately dumped in the issue thread.

All messages, and the results of the [`pkgcheck::pkgcheck()`
process](https://github.com/ropensci-review-tools/pkgcheck), are dumped to the
specified `issue_id` in the specified `repo`.

---

## 2. editorcheck_contents

The `editorcheck_contents` endpoint implements the main check functionality of
the `editorcheck` endpoint without dumping any results to the specified issue.
It is primarily intended to aid debugging any issues arising within checks,
through the use of the `stdlogs` endpoint described below. This endpoint
accepts the single argument of `repourl` only.


---

## 3. mean

A simple `mean` endpoint can be used to confirm that the server is running. It
accepts a single integer value of `n`, and returns the value of
`mean(rnorm(n))`.

---

## 4. stats_badge

This endpoint is used by the bot to extract the stats badge from those issues
which have one, in the form `"6\approved-bronze-v0.0.1"`. This is used in turn
by the bot to respond to `mint` commands used to change badge grades.

---

## 5. log

The `log` endpoint accepts a single parameter, `n`, specifying the number of
latest log entries to retrieve. An example of the log entry for [this
submission](https://github.com/ropensci/software-review/issues/470) follows:


```{r log, echo = FALSE}
"INFO [2021-10-07 16:48:14] 3.236.83.25 \"Faraday v1.7.1\" <ip>:8000 GET /editorcheck ?bot_name=ropensci-review-bot&issue_author=ewallace&issue_id=470&repo=ropensci%2Fsoftware-review&repourl=https%3A%2F%2Fgithub.com%2Fewallace%2Ftidyqpcr&sender=ewallace 200 1.964"
```

Each entry contains the following information:

1. Date and time at which call was made;
2. IP address and machine from which call was sent;
3. Method used to send call;
4. IP address to which call was delivered (always the address hosting the
   `roreviewapi` instance);
5. `http` method for the call (always `GET` for all endpoints encoded here);
6. The endpoint called (one of the methods listed above);
7. The parameters submitted along with the call;
8. The HTTP status of the call (hopefully 200); and
9. The total duration of the call response.

The 7th item of parameters submitted along with the call is particularly useful
for debugging purposes; and is specified in [this line of the `R/api.R`
file](https://github.com/ropensci-review-tools/roreviewapi/blob/e912885f516198efac885f9923318467c304df5a/R/api.R#L86).

---

## 6. clear_cache

This endpoint can be used to clear the server's cache whenever desired or
required. This cache is mainly used to store the results of calls to
[`pkgcheck::pkgcheck()`](https://github.com/ropensci-review-tools/pkgcheck).
The only effect of clearing the cache will be extra time taken to regenerate
any calls which were previously cached.

---

## 7. stdlogs

This is the most important endpoint for debugging problems within the 
[`pkgcheck`](https://github.com/ropensci-review-tools/pkgcheck) process itself.
The endpoint accepts the single parameter of `repourl`, and will return the
results of both `stdout` and `stderr` connections produced during 
[`pkgcheck`](https://github.com/ropensci-review-tools/pkgcheck). These checks
are hashed with the latest git head, ensuring that the endpoint returns checks
for the latest commit.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/editor-check.R
\name{collate_editor_check}
\alias{collate_editor_check}
\title{Collate list of checks to single concatenated character string}
\usage{
collate_editor_check(checks)
}
\arguments{
\item{checks}{Result of \code{pkgcheck::pgkcheck} function}
}
\value{
Single character
}
\description{
Collate list of checks to single concatenated character string
}
\note{
Exported only for access by plumber; not intended for general external
usage.
}
\seealso{
Other main: 
\code{\link{editor_check}()},
\code{\link{serve_api}()}
}
\concept{main}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{check_cache}
\alias{check_cache}
\title{check_cache}
\usage{
check_cache(org, repo, cache_dir = tempdir())
}
\arguments{
\item{org}{Github organization}

\item{repo}{Github repository}

\item{cache_dir}{Directory in which packages are cached}
}
\value{
FALSE If a package has previously been cached, and the github repo
has not changed; otherwise TRUE.
}
\description{
Check whether a package has been cached, and if so, whether commits have been
added to the github repo since cached version.
}
\note{
This function is not intended to be called directly, and is only
exported to enable it to be used within the \pkg{plumber} API.
}
\seealso{
Other utils: 
\code{\link{pkgrep_install_deps}()},
\code{\link{stdout_stderr_cache}()},
\code{\link{symbol_crs}()},
\code{\link{symbol_tck}()}
}
\concept{utils}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/srr.R
\name{stats_badge}
\alias{stats_badge}
\title{Get stats badge grade and standards version for a submission}
\usage{
stats_badge(repo = "ropensci/software-review", issue_num = 258)
}
\arguments{
\item{repo}{The submission repo}

\item{issue_num}{GitHub issue number of submission}
}
\value{
A single character containing the label used directly for the issue
badge
}
\description{
Get stats badge grade and standards version for a submission
}
\seealso{
Other ropensci: 
\code{\link{check_issue_template}()},
\code{\link{push_to_gh_pages}()}
}
\concept{ropensci}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dl_repo.R
\name{dl_gh_repo}
\alias{dl_gh_repo}
\title{Download a GitHub repo to local cache}
\usage{
dl_gh_repo(u)
}
\arguments{
\item{u}{URL of GitHub repository}
}
\value{
Path to locally cached '.zip' version of repository
}
\description{
Download a GitHub repo to local cache
}
\seealso{
Other github: 
\code{\link{file_pkgcheck_issue}()},
\code{\link{post_to_issue}()}
}
\concept{github}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gh-issue.R
\name{check_issue_template}
\alias{check_issue_template}
\title{Check template variables in GitHub issue}
\usage{
check_issue_template(orgrepo, issue_num)
}
\arguments{
\item{orgrepo}{GitHub organization and repo as single string separated by
forward slash (\code{org/repo}).}

\item{issue_num}{Number of issue from which to extract opening comment}
}
\value{
Comment as character string
}
\description{
Check template variables in GitHub issue
}
\seealso{
Other ropensci: 
\code{\link{push_to_gh_pages}()},
\code{\link{stats_badge}()}
}
\concept{ropensci}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/editor-check.R
\name{file_pkgcheck_issue}
\alias{file_pkgcheck_issue}
\title{File an issue in the GitHub \code{pkgcheck} repo on any packages which fail checks}
\usage{
file_pkgcheck_issue(
  repourl = NULL,
  repo = "ropensci-review-tools/pkgcheck",
  issue_id = NULL
)
}
\arguments{
\item{repourl}{The URL for the repo being checked.}

\item{repo}{The 'context.repo' parameter defining the repository from which
the command was invoked, passed in 'org/repo' format.}

\item{issue_id}{The id (number) of the issue from which the command was
invoked.}
}
\description{
File an issue in the GitHub \code{pkgcheck} repo on any packages which fail checks
}
\note{
Exported only for access by plumber; not intended for general external
usage.
}
\seealso{
Other github: 
\code{\link{dl_gh_repo}()},
\code{\link{post_to_issue}()}
}
\concept{github}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gh-pages.R
\name{push_to_gh_pages}
\alias{push_to_gh_pages}
\title{Push static \code{html} files to \code{gh-pages} branch of this repo to serve via
GitHub pages.}
\usage{
push_to_gh_pages(check)
}
\arguments{
\item{check}{Return result of \link{editor_check} function.}
}
\value{
Vector of two paths containing URLs of the \code{srr} and \code{network} files.
}
\description{
Push static \code{html} files to \code{gh-pages} branch of this repo to serve via
GitHub pages.
}
\seealso{
Other ropensci: 
\code{\link{check_issue_template}()},
\code{\link{stats_badge}()}
}
\concept{ropensci}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{symbol_crs}
\alias{symbol_crs}
\title{Cross symbol, exported for direct use in plumber API}
\usage{
symbol_crs()
}
\description{
Cross symbol, exported for direct use in plumber API
}
\seealso{
Other utils: 
\code{\link{check_cache}()},
\code{\link{pkgrep_install_deps}()},
\code{\link{stdout_stderr_cache}()},
\code{\link{symbol_tck}()}
}
\concept{utils}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roreviewapi-package.R
\docType{package}
\name{roreviewapi-package}
\alias{roreviewapi}
\alias{roreviewapi-package}
\title{roreviewapi: Plumber API to report package structure and function}
\description{
Plumber API to report package structure and function.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/roreviewapi}
  \item \url{https://github.com/ropensci-review-tools/roreviewapi}
  \item Report bugs at \url{https://github.com/ropensci-review-tools/roreviewapi/issues}
}

}
\author{
\strong{Maintainer}: Mark Padgham \email{mark.padgham@email.com} (\href{https://orcid.org/0000-0003-2172-5265}{ORCID})

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api.R
\name{serve_api}
\alias{serve_api}
\title{serve plumber API to report on packages}
\usage{
serve_api(port = 8000L, cache_dir = NULL, os = "", os_release = "")
}
\arguments{
\item{port}{Port for API to be exposed on}

\item{cache_dir}{Directory where previously downloaded repositories are
cached}

\item{os}{Name of operating system, passed to \pkg{remotes} package to
install system dependencies.}

\item{os_release}{Name of operating system release, passed to \pkg{remotes}
package to install system dependencies.}
}
\value{
Nothing; calling this starts a blocking process.
}
\description{
The API exposes the single POST points of \code{report} to download software from
the given URL and return a textual analysis of its structure and
functionality.
}
\seealso{
Other main: 
\code{\link{collate_editor_check}()},
\code{\link{editor_check}()}
}
\concept{main}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/editor-check.R
\name{editor_check}
\alias{editor_check}
\title{Body of main 'editorcheck' response}
\usage{
editor_check(repourl, repo, issue_id, post_to_issue = TRUE)
}
\arguments{
\item{repourl}{The URL for the repo being checked.}

\item{repo}{The 'context.repo' parameter defining the repository from which
the command was invoked, passed in 'org/repo' format.}

\item{issue_id}{The id (number) of the issue from which the command was
invoked.}

\item{post_to_issue}{Integer value > 0 will post results back to issue (via
'gh' cli); otherwise just return character string with result.}
}
\value{
If \code{!post_to_issue}, a markdown-formatted response body from static
package checks, otherwise URL of the issue comment to which response body has
been posted.
}
\description{
Body of main 'editorcheck' response
}
\seealso{
Other main: 
\code{\link{collate_editor_check}()},
\code{\link{serve_api}()}
}
\concept{main}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{stdout_stderr_cache}
\alias{stdout_stderr_cache}
\title{Set up stdout & stderr cache files for \code{r_bg} process}
\usage{
stdout_stderr_cache(repourl)
}
\arguments{
\item{repourl}{The URL of the repo being checked}
}
\value{
Vector of two strings holding respective local paths to \code{stdout} and
\code{stderr} files for \code{r_bg} process controlling the main \link{editor_check}
function.
}
\description{
Set up stdout & stderr cache files for \code{r_bg} process
}
\note{
These files are needed for the \pkg{callr} \code{r_bg} process which
controls the main \link{editor_check} call. See
\url{https://github.com/r-lib/callr/issues/204}. The \code{stdout} and \code{stderr}
pipes from the process are stored in the cache directory so they can be
inspected via their own distinct endpoint calls.
}
\seealso{
Other utils: 
\code{\link{check_cache}()},
\code{\link{pkgrep_install_deps}()},
\code{\link{symbol_crs}()},
\code{\link{symbol_tck}()}
}
\concept{utils}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dependencies.R
\name{pkgrep_install_deps}
\alias{pkgrep_install_deps}
\title{Install all system and package dependencies of an R package}
\usage{
pkgrep_install_deps(path, os, os_release)
}
\arguments{
\item{path}{Path to local package}

\item{os}{Name of operating system, passed to \pkg{remotes} package to
install system dependencies.}

\item{os_release}{Name of operating system release, passed to \pkg{remotes}
package to install system dependencies.}
}
\value{
Hopefully a character vector of length zero, otherwise a list of any
R packages unable to be installed.
}
\description{
Install all system and package dependencies of an R package
}
\seealso{
Other utils: 
\code{\link{check_cache}()},
\code{\link{stdout_stderr_cache}()},
\code{\link{symbol_crs}()},
\code{\link{symbol_tck}()}
}
\concept{utils}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post-to-issue.R
\name{post_to_issue}
\alias{post_to_issue}
\title{Post review checks to GitHub issue}
\usage{
post_to_issue(cmt, repo, issue_id)
}
\arguments{
\item{cmt}{Single character string with comment to post.}

\item{repo}{The repository to post to (obtained directly from bot).}

\item{issue_id}{The number of the issue to post to.}
}
\value{
URL of the comment within the nominated issue
}
\description{
Post review checks to GitHub issue
}
\seealso{
Other github: 
\code{\link{dl_gh_repo}()},
\code{\link{file_pkgcheck_issue}()}
}
\concept{github}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{is_user_authorized}
\alias{is_user_authorized}
\title{Check whether a user, identified from GitHub API token, is authorized to call
endpoints.}
\usage{
is_user_authorized()
}
\value{
Logical value indicating whether or not a user is authorized.
}
\description{
This function is used only in the \pkg{plumber} endpoints, to prevent them
being called by unauthorized users.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{symbol_tck}
\alias{symbol_tck}
\title{Tick symbol, exported for direct use in plumber API}
\usage{
symbol_tck()
}
\description{
Tick symbol, exported for direct use in plumber API
}
\seealso{
Other utils: 
\code{\link{check_cache}()},
\code{\link{pkgrep_install_deps}()},
\code{\link{stdout_stderr_cache}()},
\code{\link{symbol_crs}()}
}
\concept{utils}
