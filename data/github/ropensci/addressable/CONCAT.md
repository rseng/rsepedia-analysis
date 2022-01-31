addressable
============

![R-CMD-check](https://github.com/ropensci/addressable/workflows/R-CMD-check/badge.svg)
[![codecov](https://codecov.io/gh/ropensci/addressable/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/addressable)

Email Address Validation

## Install


```r
remotes::install_github("ropensci/addressable@main")
```


```r
library("addressable")
```

## Address


```r
x <- Address$new("User+tag@example.com")
x$host$host_name
#> [1] "example.com"
x$local$local
#> [1] "user+tag"
x$valid()
#> [1] TRUE
x$fail()
#> NULL
```


```r
x <- Address$new("user1")
x$valid()
#> [1] FALSE
x$fail()
#> [1] "Invalid Domain Name"
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/addressable/issues).
* License: MIT
* Get citation information for `addressable` in R doing `citation(package = 'addressable')`
* Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[coc]: https://github.com/ropensci/addressable/blob/maddressable/CODE_OF_CONDUCT.md
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
(https://contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/
<!--- Provide a general summary of your changes in the Title above -->

## Description
<!--- Describe your changes in detail -->

## Related Issue
<!--- if this closes an issue make sure include e.g., "fix #4"
or similar - or if just relates to an issue make sure to mention
it like "#4" -->

## Example
<!--- if introducing a new feature or changing behavior of existing
methods/functions, include an example if possible to do in brief form -->

<!--- Did you remember to include tests? Unless you're just changing
grammar, please include new tests for your change -->
# CONTRIBUTING #

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/addressable/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/addressable.git`
* Make sure to track progress upstream (i.e., on our version of `addressable` at `ropensci/addressable`) by doing `git remote add upstream https://github.com/ropensci/addressable.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/addressable`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
addressable
============

```{r echo=FALSE}
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

![R-CMD-check](https://github.com/ropensci/addressable/workflows/R-CMD-check/badge.svg)
[![codecov](https://codecov.io/gh/ropensci/addressable/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/addressable)


Email Address Validation

## Install

```{r eval=FALSE}
remotes::install_github("ropensci/addressable@main")
```

```{r}
library("addressable")
```

## Address

```{r}
x <- Address$new("User+tag@example.com")
x$host$host_name
x$local$local
x$valid()
x$fail()
```

```{r}
x <- Address$new("user1")
x$valid()
x$fail()
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/addressable/issues).
* License: MIT
* Get citation information for `addressable` in R doing `citation(package = 'addressable')`
* Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[coc]: https://github.com/ropensci/addressable/blob/maddressable/CODE_OF_CONDUCT.md
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Address.R
\name{Address}
\alias{Address}
\title{Address}
\description{
Address Class
}
\examples{
\dontrun{
x <- Address$new("User+tag@example.com")
x
x$local
x$local$valid()
x$host
x$host$valid()
x$valid()
x$fail()
x$munge()
x$normal()

x <- Address$new("User+tag@example.com",
  config = list(munge_string = "*^*^*^*"))
x$munge()
}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{method}}{(character) one or more URLs}

\item{\code{url}}{(character) one or more URLs}
}
\if{html}{\out{</div>}}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{method}}{(character) one or more URLs}

\item{\code{url}}{(character) one or more URLs}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Address$new()}}
\item \href{#method-split_local_host}{\code{Address$split_local_host()}}
\item \href{#method-normal}{\code{Address$normal()}}
\item \href{#method-to_s}{\code{Address$to_s()}}
\item \href{#method-munge}{\code{Address$munge()}}
\item \href{#method-canonical}{\code{Address$canonical()}}
\item \href{#method-base}{\code{Address$base()}}
\item \href{#method-reference}{\code{Address$reference()}}
\item \href{#method-valid}{\code{Address$valid()}}
\item \href{#method-fail}{\code{Address$fail()}}
\item \href{#method-clone}{\code{Address$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new Address object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Address$new(email_address, config = list())}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{email_address}}{(character) an email address}

\item{\code{config}}{(list) list of config options}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-split_local_host"></a>}}
\if{latex}{\out{\hypertarget{method-split_local_host}{}}}
\subsection{Method \code{split_local_host()}}{
split local host
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Address$split_local_host(email)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{email}}{(character) email address}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
character string
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-normal"></a>}}
\if{latex}{\out{\hypertarget{method-normal}{}}}
\subsection{Method \code{normal()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Address$normal()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-to_s"></a>}}
\if{latex}{\out{\hypertarget{method-to_s}{}}}
\subsection{Method \code{to_s()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Address$to_s()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-munge"></a>}}
\if{latex}{\out{\hypertarget{method-munge}{}}}
\subsection{Method \code{munge()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Address$munge()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-canonical"></a>}}
\if{latex}{\out{\hypertarget{method-canonical}{}}}
\subsection{Method \code{canonical()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Address$canonical()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-base"></a>}}
\if{latex}{\out{\hypertarget{method-base}{}}}
\subsection{Method \code{base()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Address$base()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-reference"></a>}}
\if{latex}{\out{\hypertarget{method-reference}{}}}
\subsection{Method \code{reference()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Address$reference()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-valid"></a>}}
\if{latex}{\out{\hypertarget{method-valid}{}}}
\subsection{Method \code{valid()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Address$valid(options = list())}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fail"></a>}}
\if{latex}{\out{\hypertarget{method-fail}{}}}
\subsection{Method \code{fail()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Address$fail()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Address$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Host.R
\name{Host}
\alias{Host}
\title{Host}
\description{
Host Class
}
\examples{
\dontrun{
x <- Host$new("example.com")
x
x$host_name
x$domain_name
x$registration_name
x$tld
x$valid()
# x$munge()
x$error
x$error_message
x$ip_address
x$config$config

x <- Host$new("gmail.com", config = list(dns_lookup = "off"))
x$config$config$dns_lookup
x$valid()

x <- Host$new("gmail.com")
x
x$valid()
x$fail()

x <- Host$new("gm ail.com")
x
x$valid()
x$fail()

x <- Host$new("user1")
x$domain_name
x$dns_name
x$registration_name
}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Host$new()}}
\item \href{#method-name}{\code{Host$name()}}
\item \href{#method-canonical}{\code{Host$canonical()}}
\item \href{#method-munge}{\code{Host$munge()}}
\item \href{#method-parse}{\code{Host$parse()}}
\item \href{#method-parse_comment}{\code{Host$parse_comment()}}
\item \href{#method-make_host_name}{\code{Host$make_host_name()}}
\item \href{#method-fully_qualified_domain_name}{\code{Host$fully_qualified_domain_name()}}
\item \href{#method-hosted_service}{\code{Host$hosted_service()}}
\item \href{#method-find_provider}{\code{Host$find_provider()}}
\item \href{#method-set_provider}{\code{Host$set_provider()}}
\item \href{#method-ip}{\code{Host$ip()}}
\item \href{#method-ipv4}{\code{Host$ipv4()}}
\item \href{#method-ipv6}{\code{Host$ipv6()}}
\item \href{#method-matches}{\code{Host$matches()}}
\item \href{#method-registration_name_matches}{\code{Host$registration_name_matches()}}
\item \href{#method-tld_matches}{\code{Host$tld_matches()}}
\item \href{#method-provider_matches}{\code{Host$provider_matches()}}
\item \href{#method-domain_matches}{\code{Host$domain_matches()}}
\item \href{#method-ip_matches}{\code{Host$ip_matches()}}
\item \href{#method-dns_enabled}{\code{Host$dns_enabled()}}
\item \href{#method-dns_a_record}{\code{Host$dns_a_record()}}
\item \href{#method-exchangers}{\code{Host$exchangers()}}
\item \href{#method-valid}{\code{Host$valid()}}
\item \href{#method-fail}{\code{Host$fail()}}
\item \href{#method-valid_dns}{\code{Host$valid_dns()}}
\item \href{#method-valid_mx}{\code{Host$valid_mx()}}
\item \href{#method-valid_format}{\code{Host$valid_format()}}
\item \href{#method-valid_ip}{\code{Host$valid_ip()}}
\item \href{#method-localhost}{\code{Host$localhost()}}
\item \href{#method-connect}{\code{Host$connect()}}
\item \href{#method-clone}{\code{Host$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$new(host_name, config = list())}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-name"></a>}}
\if{latex}{\out{\hypertarget{method-name}{}}}
\subsection{Method \code{name()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$name()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-canonical"></a>}}
\if{latex}{\out{\hypertarget{method-canonical}{}}}
\subsection{Method \code{canonical()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$canonical()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-munge"></a>}}
\if{latex}{\out{\hypertarget{method-munge}{}}}
\subsection{Method \code{munge()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$munge()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-parse"></a>}}
\if{latex}{\out{\hypertarget{method-parse}{}}}
\subsection{Method \code{parse()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$parse(host)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-parse_comment"></a>}}
\if{latex}{\out{\hypertarget{method-parse_comment}{}}}
\subsection{Method \code{parse_comment()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$parse_comment(host)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-make_host_name"></a>}}
\if{latex}{\out{\hypertarget{method-make_host_name}{}}}
\subsection{Method \code{make_host_name()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$make_host_name(name)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fully_qualified_domain_name"></a>}}
\if{latex}{\out{\hypertarget{method-fully_qualified_domain_name}{}}}
\subsection{Method \code{fully_qualified_domain_name()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$fully_qualified_domain_name(host_part)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-hosted_service"></a>}}
\if{latex}{\out{\hypertarget{method-hosted_service}{}}}
\subsection{Method \code{hosted_service()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$hosted_service()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-find_provider"></a>}}
\if{latex}{\out{\hypertarget{method-find_provider}{}}}
\subsection{Method \code{find_provider()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$find_provider()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_provider"></a>}}
\if{latex}{\out{\hypertarget{method-set_provider}{}}}
\subsection{Method \code{set_provider()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$set_provider(name, provider_config = list())}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ip"></a>}}
\if{latex}{\out{\hypertarget{method-ip}{}}}
\subsection{Method \code{ip()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$ip()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ipv4"></a>}}
\if{latex}{\out{\hypertarget{method-ipv4}{}}}
\subsection{Method \code{ipv4()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$ipv4()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ipv6"></a>}}
\if{latex}{\out{\hypertarget{method-ipv6}{}}}
\subsection{Method \code{ipv6()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$ipv6()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-matches"></a>}}
\if{latex}{\out{\hypertarget{method-matches}{}}}
\subsection{Method \code{matches()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$matches(rules)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-registration_name_matches"></a>}}
\if{latex}{\out{\hypertarget{method-registration_name_matches}{}}}
\subsection{Method \code{registration_name_matches()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$registration_name_matches(rule)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-tld_matches"></a>}}
\if{latex}{\out{\hypertarget{method-tld_matches}{}}}
\subsection{Method \code{tld_matches()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$tld_matches(rule)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-provider_matches"></a>}}
\if{latex}{\out{\hypertarget{method-provider_matches}{}}}
\subsection{Method \code{provider_matches()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$provider_matches(rule)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-domain_matches"></a>}}
\if{latex}{\out{\hypertarget{method-domain_matches}{}}}
\subsection{Method \code{domain_matches()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$domain_matches(rule)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ip_matches"></a>}}
\if{latex}{\out{\hypertarget{method-ip_matches}{}}}
\subsection{Method \code{ip_matches()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$ip_matches(cidr)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-dns_enabled"></a>}}
\if{latex}{\out{\hypertarget{method-dns_enabled}{}}}
\subsection{Method \code{dns_enabled()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$dns_enabled()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-dns_a_record"></a>}}
\if{latex}{\out{\hypertarget{method-dns_a_record}{}}}
\subsection{Method \code{dns_a_record()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$dns_a_record()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-exchangers"></a>}}
\if{latex}{\out{\hypertarget{method-exchangers}{}}}
\subsection{Method \code{exchangers()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$exchangers()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-valid"></a>}}
\if{latex}{\out{\hypertarget{method-valid}{}}}
\subsection{Method \code{valid()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$valid(rules = list())}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fail"></a>}}
\if{latex}{\out{\hypertarget{method-fail}{}}}
\subsection{Method \code{fail()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$fail()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-valid_dns"></a>}}
\if{latex}{\out{\hypertarget{method-valid_dns}{}}}
\subsection{Method \code{valid_dns()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$valid_dns()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-valid_mx"></a>}}
\if{latex}{\out{\hypertarget{method-valid_mx}{}}}
\subsection{Method \code{valid_mx()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$valid_mx()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-valid_format"></a>}}
\if{latex}{\out{\hypertarget{method-valid_format}{}}}
\subsection{Method \code{valid_format()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$valid_format()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-valid_ip"></a>}}
\if{latex}{\out{\hypertarget{method-valid_ip}{}}}
\subsection{Method \code{valid_ip()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$valid_ip()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-localhost"></a>}}
\if{latex}{\out{\hypertarget{method-localhost}{}}}
\subsection{Method \code{localhost()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$localhost(x)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-connect"></a>}}
\if{latex}{\out{\hypertarget{method-connect}{}}}
\subsection{Method \code{connect()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$connect()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Host$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Local.R
\name{Local}
\alias{Local}
\title{Local}
\description{
Local Class
}
\examples{
\dontrun{
host <- Host$new("gmail.com")
x <- Local$new("very.common", list(), host)
x
x$valid()
x$local
x$format()
x$munge()
x$error
x$error_message
x$config$config
}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Local$new()}}
\item \href{#method-set_local}{\code{Local$set_local()}}
\item \href{#method-parse}{\code{Local$parse()}}
\item \href{#method-parse_comment}{\code{Local$parse_comment()}}
\item \href{#method-parse_tag}{\code{Local$parse_tag()}}
\item \href{#method-is_ascii}{\code{Local$is_ascii()}}
\item \href{#method-is_unicode}{\code{Local$is_unicode()}}
\item \href{#method-is_redacted}{\code{Local$is_redacted()}}
\item \href{#method-format}{\code{Local$format()}}
\item \href{#method-to_s}{\code{Local$to_s()}}
\item \href{#method-conventional}{\code{Local$conventional()}}
\item \href{#method-canonical}{\code{Local$canonical()}}
\item \href{#method-relax}{\code{Local$relax()}}
\item \href{#method-standard}{\code{Local$standard()}}
\item \href{#method-munge}{\code{Local$munge()}}
\item \href{#method-valid}{\code{Local$valid()}}
\item \href{#method-valid_size}{\code{Local$valid_size()}}
\item \href{#method-valid_size_checks}{\code{Local$valid_size_checks()}}
\item \href{#method-valid_encoding}{\code{Local$valid_encoding()}}
\item \href{#method-is_conventional}{\code{Local$is_conventional()}}
\item \href{#method-is_relaxed}{\code{Local$is_relaxed()}}
\item \href{#method-is_standard}{\code{Local$is_standard()}}
\item \href{#method-clone}{\code{Local$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$new(local, config = list(), host = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_local"></a>}}
\if{latex}{\out{\hypertarget{method-set_local}{}}}
\subsection{Method \code{set_local()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$set_local(raw)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-parse"></a>}}
\if{latex}{\out{\hypertarget{method-parse}{}}}
\subsection{Method \code{parse()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$parse(raw)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-parse_comment"></a>}}
\if{latex}{\out{\hypertarget{method-parse_comment}{}}}
\subsection{Method \code{parse_comment()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$parse_comment(raw)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-parse_tag"></a>}}
\if{latex}{\out{\hypertarget{method-parse_tag}{}}}
\subsection{Method \code{parse_tag()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$parse_tag(raw)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-is_ascii"></a>}}
\if{latex}{\out{\hypertarget{method-is_ascii}{}}}
\subsection{Method \code{is_ascii()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$is_ascii()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-is_unicode"></a>}}
\if{latex}{\out{\hypertarget{method-is_unicode}{}}}
\subsection{Method \code{is_unicode()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$is_unicode()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-is_redacted"></a>}}
\if{latex}{\out{\hypertarget{method-is_redacted}{}}}
\subsection{Method \code{is_redacted()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$is_redacted()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-format"></a>}}
\if{latex}{\out{\hypertarget{method-format}{}}}
\subsection{Method \code{format()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$format(form = self$config$config[["local_format"]] \%||\% "conventional")}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-to_s"></a>}}
\if{latex}{\out{\hypertarget{method-to_s}{}}}
\subsection{Method \code{to_s()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$to_s()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-conventional"></a>}}
\if{latex}{\out{\hypertarget{method-conventional}{}}}
\subsection{Method \code{conventional()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$conventional()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-canonical"></a>}}
\if{latex}{\out{\hypertarget{method-canonical}{}}}
\subsection{Method \code{canonical()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$canonical()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-relax"></a>}}
\if{latex}{\out{\hypertarget{method-relax}{}}}
\subsection{Method \code{relax()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$relax()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-standard"></a>}}
\if{latex}{\out{\hypertarget{method-standard}{}}}
\subsection{Method \code{standard()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$standard()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-munge"></a>}}
\if{latex}{\out{\hypertarget{method-munge}{}}}
\subsection{Method \code{munge()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$munge()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-valid"></a>}}
\if{latex}{\out{\hypertarget{method-valid}{}}}
\subsection{Method \code{valid()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$valid(
  format = self$config$config[["local_format"]] \%||\% "conventional"
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-valid_size"></a>}}
\if{latex}{\out{\hypertarget{method-valid_size}{}}}
\subsection{Method \code{valid_size()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$valid_size()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-valid_size_checks"></a>}}
\if{latex}{\out{\hypertarget{method-valid_size_checks}{}}}
\subsection{Method \code{valid_size_checks()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$valid_size_checks(range)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-valid_encoding"></a>}}
\if{latex}{\out{\hypertarget{method-valid_encoding}{}}}
\subsection{Method \code{valid_encoding()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$valid_encoding(
  enc = self$config$config[["local_encoding"]] \%||\% "ascii"
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-is_conventional"></a>}}
\if{latex}{\out{\hypertarget{method-is_conventional}{}}}
\subsection{Method \code{is_conventional()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$is_conventional()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-is_relaxed"></a>}}
\if{latex}{\out{\hypertarget{method-is_relaxed}{}}}
\subsection{Method \code{is_relaxed()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$is_relaxed()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-is_standard"></a>}}
\if{latex}{\out{\hypertarget{method-is_standard}{}}}
\subsection{Method \code{is_standard()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$is_standard()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Local$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addressable-package.R
\docType{package}
\name{addressable-package}
\alias{addressable-package}
\alias{addressable}
\title{addressable}
\description{
Email Address Validation
}
\author{
Scott Chamberlain \email{sckott@protonmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Config.R
\name{Config}
\alias{Config}
\title{Config}
\description{
Config Class
}
\examples{
\dontrun{
x <- Config$new()
x
x$config
x$config$dns_timeout

z <- Config$new(list(dns_timeout = 4))
z
z$config
z$config$dns_timeout
}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-all_settings}{\code{Config$all_settings()}}
\item \href{#method-new}{\code{Config$new()}}
\item \href{#method-provider}{\code{Config$provider()}}
\item \href{#method-error_message}{\code{Config$error_message()}}
\item \href{#method-error_messages}{\code{Config$error_messages()}}
\item \href{#method-configure}{\code{Config$configure()}}
\item \href{#method-clone}{\code{Config$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-all_settings"></a>}}
\if{latex}{\out{\hypertarget{method-all_settings}{}}}
\subsection{Method \code{all_settings()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Config$all_settings(x)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Config$new(config = list())}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-provider"></a>}}
\if{latex}{\out{\hypertarget{method-provider}{}}}
\subsection{Method \code{provider()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Config$provider(name, config = list())}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-error_message"></a>}}
\if{latex}{\out{\hypertarget{method-error_message}{}}}
\subsection{Method \code{error_message()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Config$error_message(name, locale = "en")}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-error_messages"></a>}}
\if{latex}{\out{\hypertarget{method-error_messages}{}}}
\subsection{Method \code{error_messages()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Config$error_messages(x = list(), locale = "en")}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-configure"></a>}}
\if{latex}{\out{\hypertarget{method-configure}{}}}
\subsection{Method \code{configure()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Config$configure(settings)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Config$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Exchanger.R
\name{Exchanger}
\alias{Exchanger}
\title{Exchanger}
\description{
Exchanger Class
}
\examples{
\dontrun{
x <- Exchanger$new("gmail.com")
x
x$cached("gmail.com")
x
x$mxers_
x$mxers()
}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-cached}{\code{Exchanger$cached()}}
\item \href{#method-new}{\code{Exchanger$new()}}
\item \href{#method-provider}{\code{Exchanger$provider()}}
\item \href{#method-mxers}{\code{Exchanger$mxers()}}
\item \href{#method-mx_ips}{\code{Exchanger$mx_ips()}}
\item \href{#method-matches}{\code{Exchanger$matches()}}
\item \href{#method-in_cidr}{\code{Exchanger$in_cidr()}}
\item \href{#method-clone}{\code{Exchanger$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cached"></a>}}
\if{latex}{\out{\hypertarget{method-cached}{}}}
\subsection{Method \code{cached()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Exchanger$cached(host, config = list())}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Exchanger$new(host, config = list())}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-provider"></a>}}
\if{latex}{\out{\hypertarget{method-provider}{}}}
\subsection{Method \code{provider()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Exchanger$provider()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mxers"></a>}}
\if{latex}{\out{\hypertarget{method-mxers}{}}}
\subsection{Method \code{mxers()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Exchanger$mxers()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-mx_ips"></a>}}
\if{latex}{\out{\hypertarget{method-mx_ips}{}}}
\subsection{Method \code{mx_ips()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Exchanger$mx_ips()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-matches"></a>}}
\if{latex}{\out{\hypertarget{method-matches}{}}}
\subsection{Method \code{matches()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Exchanger$matches(rules)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-in_cidr"></a>}}
\if{latex}{\out{\hypertarget{method-in_cidr}{}}}
\subsection{Method \code{in_cidr()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Exchanger$in_cidr(cidr)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Exchanger$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
