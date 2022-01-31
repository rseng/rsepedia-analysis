# rrlite

[![Build Status](https://travis-ci.org/ropensci/rrlite.png?branch=master)](https://travis-ci.org/ropensci/rrlite)
[![Coverage Status](https://coveralls.io/repos/ropensci/rrlite/badge.svg?branch=master)](https://coveralls.io/r/ropensci/rrlite?branch=master)



R interface to [rlite](https://github.com/seppo0010/rlite).  rlite is a "self-contained, serverless, zero-configuration, transactional redis-compatible database engine. rlite is to Redis what SQLite is to SQL.  And `Redis` is a *data structures* server; at the simplest level it can be used as a key-value store, but it can store other data types (hashes, lists, sets and more).

This package is designed to follow exactly the same interface as [redux](https://github.com/richfitz/redux).

# Usage

See [`redux`](https://github.com/richfitz/redux) for more details.

The main function here is `rrlite::hirlite` that creates a `redis_api` object that exposes the full Redis API.


```r
con <- rrlite::hirlite()
```


```r
con
```

```
## <redis_api>
##   Redis commands:
##     APPEND: function
##     AUTH: function
##     BGREWRITEAOF: function
##     BGSAVE: function
##     ...
##     ZSCORE: function
##     ZUNIONSTORE: function
##   Other public methods:
##     clone: function
##     command: function
##     config: function
##     initialize: function
##     pipeline: function
##     reconnect: function
##     subscribe: function
##     type: function
```

This object has all the same methods as the corresponding object created by `redux::hiredis()` but operating on a rlite database.  The default database uses the magic path `:memory:` but persistent on-disk storage is possible (see `?rlite_config`).

All the usual Redis-type things work:


```r
con$SET("mykey", "mydata")
```

```
## [Redis: OK]
```

```r
con$GET("mykey")
```

```
## [1] "mydata"
```

As with redux, commands are vectorised:


```r
con$MSET(c("a", "b", "c"), c(1, 2, 3))
```

```
## [Redis: OK]
```

```r
con$MGET(c("a", "b", "c"))
```

```
## [[1]]
## [1] "1"
##
## [[2]]
## [1] "2"
##
## [[3]]
## [1] "3"
```

# Approach

This package aims to be a drop-in self-contained replacement for `redux` without requiring  `Redis` server.  Therefore almost the entire package (and tests) is automaticaly generated from `redux`.  The only installed files not generated are:

* R/hirlite.R (because documentation)
* src/subscribe.c (just a stub)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rrlite/issues).
* License: GPL
* Get citation information for `rrlite` in R by doing `citation(package = 'rrlite')`

[![rofooter](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
# rrlite

[![Build Status](https://travis-ci.org/ropensci/rrlite.png?branch=master)](https://travis-ci.org/ropensci/rrlite)
[![Coverage Status](https://coveralls.io/repos/ropensci/rrlite/badge.svg?branch=master)](https://coveralls.io/r/ropensci/rrlite?branch=master)
[![](https://badges.ropensci.org/6_status.svg)](https://github.com/ropensci/onboarding/issues/6)

```{r, echo=FALSE, results="hide"}
knitr::opts_chunk$set(error=FALSE)
```

R interface to [rlite](https://github.com/seppo0010/rlite).  rlite is a "self-contained, serverless, zero-configuration, transactional redis-compatible database engine. rlite is to Redis what SQLite is to SQL.  And `Redis` is a *data structures* server; at the simplest level it can be used as a key-value store, but it can store other data types (hashes, lists, sets and more).

This package is designed to follow exactly the same interface as [redux](https://github.com/richfitz/redux).

# Usage

See [`redux`](https://github.com/richfitz/redux) for more details.

The main function here is `rrlite::hirlite` that creates a `redis_api` object that exposes the full Redis API.

```{r}
con <- rrlite::hirlite()
```

```{r, eval=FALSE}
con
```
```{r, echo=FALSE}
res <- capture.output(print(con))
res <- c(res[1:6], "    ...",
         res[(max(grep("\\s+[A-Z]", res)) - 2):length(res)])
writeLines(res)
```

This object has all the same methods as the corresponding object created by `redux::hiredis()` but operating on a rlite database.  The default database uses the magic path `:memory:` but persistent on-disk storage is possible (see `?rlite_config`).

All the usual Redis-type things work:

```{r}
con$SET("mykey", "mydata")
con$GET("mykey")
```

As with redux, commands are vectorised:

```{r}
con$MSET(c("a", "b", "c"), c(1, 2, 3))
con$MGET(c("a", "b", "c"))
```

# Approach

This package aims to be a drop-in self-contained replacement for `redux` without requiring  `Redis` server.  Therefore almost the entire package (and tests) is automaticaly generated from `redux`.  The only installed files not generated are:

* R/hirlite.R (because documentation)
* src/subscribe.c (just a stub)

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/rrlite/issues).
* License: GPL
* Get citation information for `rrlite` in R by doing `citation(package = 'rrlite')`

[![rofooter](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connection.R
\name{rlite_connection}
\alias{rlite_connection}
\title{Create a rlite connection}
\usage{
rlite_connection(config = rlite_config())
}
\arguments{
\item{config}{Configuration parameters as generated by
\code{\link{rlite_config}}}
}
\description{
Create a rlite connection.  This function is designed to be used
in other packages, and not directly by end-users.  However, it is
possible and safe to use.  See the \code{\link{hirlite}} package
for the user friendly interface.
}
\details{
This function creates a list of functions, appropriately bound to
a pointer to a rlite connection.  This is designed for package
authors to use so without having to ever deal with the actual
pointer itself (which cannot be directly manipulated from R
anyway).

The returned list has elements, all of which are functions:

\describe{
\item{\code{config()}}{The configuration information}

\item{\code{reconnect()}}{Attempt reconnection of a connection
that has been closed, through serialisation/deserialiation or
through loss of internet connection.}

\item{command(cmd)}{Run a Redis command.  The format of this
command will be documented elsewhere.}

\item{pipeline(cmds)}{Run a pipeline of Redis commands.}

\item{subscribe(channel, pattern, callback, envir)}{Subscribe to a
channel or pattern specifying channels.  Here, \code{channel} must
be a character vector, \code{pattern} a logical indicating if
\code{channel} should be interpreted as a pattern, \code{callback}
is a function to apply to each recieved message, returning
\code{TRUE} when subscription should stop, and \code{envir} is the
environment in which to evaluate \code{callback}.  See below.}

}
}
\section{Subscriptions}{


  The callback function must take a single argument; this will be
  the recieved message with named elements \code{type} (which will
  be message), \code{channel} (the name of the channel) and
  \code{value} (the message contents).  If \code{pattern} was
  \code{TRUE}, then an additional element \code{pattern} will be
  present (see the Redis docs).  The callback must return
  \code{TRUE} or \code{FALSE}; this indicates if the client should
  continue quit (i.e., \code{TRUE} means return control to R,
  \code{FALSE} means keep going).

  Because the \code{subscribe} function is blocking and returns
  nothing, so all data collection needs to happen as a side-effect
  of the callback function.

  There is currently no way of interrupting the client while it is
  waiting for a message.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hirlite.R
\name{hirlite}
\alias{hirlite}
\alias{rlite_available}
\title{Interface to rlite}
\usage{
hirlite(...)

rlite_available(...)
}
\arguments{
\item{...}{Named configuration options passed to
\code{\link{redis_config}}, used to create the environment
(notable keys include \code{host}, \code{port}, and the
environment variable \code{REDIS_URL}).  In addition to the
\code{Redux} treatment of the configuration, \code{RLITE_URL}
takes precendence over \code{REDIS_URL}, and a host of
\code{localhost} or \code{127.0.0.1} will be treated as an
in-memory database (\code{:memory:}).}
}
\description{
Create an interface to rlite, with a generated interface to all
rlite commands (using \code{Redux}).
}
\examples{
r <- hirlite()
r$PING()
r$SET("foo", "bar")
r$GET("foo")
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rlite_config.R
\name{rlite_config}
\alias{rlite_config}
\title{rlite configuration}
\usage{
rlite_config(...)
}
\arguments{
\item{...}{Arguments passed to \code{\link{redis_config}}; see
that file for more information.}
}
\description{
rlite configuration settings.  Based on the \code{redis_config}
function but with additional tweaks for rlite.  The differences
between this configuration and \code{\link{redis_config}} is that:
}
\details{
\itemize{

\item{\code{RLITE_URL} takes precendence over \code{REDIS_URL} if
both are present (otherwise \code{REDIS_URL} will still be used).}

\item{A host of \code{localhost} or \code{127.0.0.1}, which is
\code{redis_config}'s default, will map to a filename of
\code{:memory:} for a transient in-memory store.}
}

The \code{port} entry will be ignored, but the \code{password} and
\code{db} entries will be used if present.  \code{path} is
equivalent to \code{host}.
}

