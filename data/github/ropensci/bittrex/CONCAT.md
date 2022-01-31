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
<!-- README.md is generated from README.Rmd. Please edit that file -->

bittrex: An R client for the [Bittrex Crypto-Currency Exchange](https://bittrex.com)
=====================================================

**Authors:** Michael J. Kane<br/>
**License:** [LGPL-2](https://opensource.org/licenses/LGPL-2.1)

[![Build Status](https://travis-ci.org/ropensci/bittrex.svg?branch=master)](https://travis-ci.org/ropensci/bittrex)
[![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/bittrex/master.svg)](https://codecov.io/github/ropensci/bittrex?branch=master)
[![](https://badges.ropensci.org/120_status.svg)](https://github.com/ropensci/onboarding/issues/120)
[![DOI](https://zenodo.org/badge/88363176.svg)](https://zenodo.org/badge/latestdoi/88363176)

Disclaimer
===

This software is in no way affiliated, endorsed, or approved by the
[Bittrex crypto-currency exchange](https://bittrex.com/) or any of its affiliates. 
It comes with absolutely no warranty and should not be used in actual trading 
unless the user can read and understand the source and know what you are doing.

Overview
===

Package 'bittrex' is an R implementation of the REST interface used by the [Bittrex
crypto-currency exchange](https://bittrex.com/). It provides functions 
for endpoints supported by the exchange. This includes the ability 
to retrieve price, volume, and order book information as well as the ability
to trade crypto-currencies.

Calls to the exchange are categorized as either public, which includes 
requests for price, volume, and order book information, and private, which 
includes all requests requiring an account including placing buy or sell 
orders. Public calls can be used directly by installing the package. 
Private calls require that you 
[create an account](https://https://bittrex.com/account/Register) and create an API and secret 
key with appropriate permissions.

Private calls retrieve the API and secret key using the BITTREX_API_KEY and 
BITTREX_SECRET_KEY environment variables. These may be set by the user before 
opening the R session or, they can be set using the 'bittrex_authenticate' 
function.

Quickstart
===

Install
---

The package is available from GitHub and will be uploaded to CRAN
shortly. If you wish to install the development version then install the 
[devtools package](https://CRAN.R-project.org/package=devtools), available 
from CRAN. 


```r
#install.packages("devtools")
devtools::install_github("ropensci/bittrex")
```

Using the Package
---

After installation, you may query the exchange with any of the public
calls. For example, if we want to see the spread of the cost of doge 
coins in bitcoins, we can use the following code.


```r
library(bittrex)
library(scales)
library(ggplot2)

# The price of doge coins in bitcoins.
doge_btc = bt_getmarkethistory(market='btc-doge')$result

ggplot(doge_btc, aes(x=time_stamp, y=price, group=order_type, 
  color=order_type)) + geom_line() + 
  scale_x_datetime(breaks=date_breaks("hours"), 
    labels=date_format("%m-%d %H:%M")) + xlab("Date and Time") +
  ylab("Price") + scale_colour_discrete(name="Order Type")
```

![](inst/README_files/readme-viz.png)

Contributing
---

If you would like to contribute to the project please contact
the maintainer directly. Please note that this project is
released with a [Contributor Code of Conduct](CONDUCT.md).
By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# 0.1.8

## BUG FIXES
* Fixes to documentation.

## MINOR IMPROVEMENTS
* Documentation enhancements.


# bittrex 0.1.0

* First version of the package including:
    * Market summaries, price, order book, and volume for Bittrex crypto-currencies.
    * The ability to trade crypto-currencies with a user account.
---
title: 'bittrex: An R client for the Bittrex Crypto-Currency Exchange'
tags:
  - finance
  - crypto-currency
  - open exchanges
authors:
 - name: Michael Kane
   orcid: 0000-0003-1899-6662
   affiliation: 1
affiliations:
 - name: Yale University
   index: 1
date: 11 July 2017
bibliography: paper.bib
---

# Summary

Package ```bittrex``` [@bittrex] is a RESTful R [@R] client for the
Bittrex crypto-currency exchange [@bittrex_site]. The package provides 
functions for all endpoints supported by the exchange including the
ability to retrieve price, volume, and order book information as well as
the ability to trade crypto-currencies.

Calls to the exchange are categorized as either public, which includes
requests for price, volume, and order book information, and private,
which includes all requests requiring an account including placing buy
or sell orders. Public calls can be used directly by installing the
package. Private calls require that you [create an
account](https://https://bittrex.com/account/Register) and create an API 
and secret key with appropriate permissions.

Private calls retrieve the API and secret key using the BITTREX\_API\_KEY
and BITTREX\_SECRET\_KEY environment variables. These may be set by the
user before opening the R session or, they can be set using the
'bittrex\_authenticate' function.

While this package can be used to trade currencies, its primary goal is to 
allow finance and economics researchers 
to gather exchange data, including transaction prices,
orderbooks, and volumes for individual crypto-currencies. A few research
avenues currently being explored include:
- Comparisons between semi-regulated markets (bittrex) to traditional, regulated markets.
- Microfinancing for start-ups using initial coin offerings.
- Quantifying the value of non-traditional products and services.
- Systemic risk of multiple crypto-currencies on the same block chain.
- Crypto-coins as a proxy for sports betting where sports betting is illegal.


# References
  
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/private-interface.r
\name{bt_getorder}
\alias{bt_getorder}
\title{Order Data for a Specified Order}
\usage{
bt_getorder(uuid)
}
\arguments{
\item{uuid}{the uuid of the order.}
}
\value{
A named list with the following elements:
\itemize{
\item{success: }{a boolean indicating whether the request was successful.}
\item{message: }{a string describing the error if the request was not
successful, otherwise an empty string.}
\item{result:  }{a \code{data.frame} providing information about the
open order including (but not limited to) the market, quantity remaining
in the order, the type of order, and when the order was opened.
}
}
}
\description{
The \code{bt_getorder()} function retrieves open order data
on \url{https://bittrex.com}. This function can only be used after you provide
information for authentication.
}
\examples{
\dontrun{
# Note you must authenticate and define a uuid first.
bt_getorder(uuid)
# $success
# [1] TRUE
#
# $message
# [1] ""
#
# $result
#   account_id                           order_uuid exchange      type quantity
# 1         NA 63181c27-dd14-476e-960c-1bd8366bb312  ETH-LTC LIMIT_BUY        1
#   quantity_remaining limit reserved reserve_remaining commission_reserved
# 1                  1 0.001    0.001             0.001             2.5e-06
#   commission_reserve_remaining commission_paid price price_per_unit
# 1                      2.5e-06               0     0             NA
#                    opened closed is_open                             sentinel
# 1 2017-07-11T17:05:55.583     NA    TRUE 78b210fa-a156-4f4d-8043-51163740056f
#   cancel_initiated immediate_or_cancel is_conditional condition
# 1            FALSE               FALSE          FALSE      NONE
#   condition_target
# 1               NA
}
}
\references{
https://bittrex.com/api/v1.1/account/getorder
}
\seealso{
\code{\link[=bt_authenticate]{bt_authenticate()}}, \code{\link[=bt_getopenorders]{bt_getopenorders()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public-interface.r
\name{bt_api_check}
\alias{bt_api_check}
\title{Check the Connection to the Bittrex Exchange}
\usage{
bt_api_check(warn = TRUE)
}
\arguments{
\item{warn}{if the request is not successful, should a warning be provided
with the status code.}
}
\value{
A named logical indicating if you can connect to the exchange
through the public interface.
}
\description{
The \code{bt_api_check()} function checks to see
if you can successfully connect to \url{https://bittrex.com}.
}
\examples{
\dontrun{
bt_api_check()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package-doc.r
\docType{package}
\name{bittrex-package}
\alias{bittrex-package}
\alias{bittrex}
\title{Client for the Bittrex Crypto-Currency Exchange}
\description{
This software is in no way affiliated, endorsed, or approved by
the Bittrex crypto-currency exchange or any of its affiliates. It comes with
absolutely no warranty and should not be used in actual trading
unless the user can read and understand the source and knows what they are
doing.

The \code{bittrex} package is an R implementation of the REST interface used
by the Bittrex crypto-currency exchange (\url{https://bittrex.com}). It
provides functions for all endpoints currently (as of May 16, 2017)
supported by the exchange. This includes the ability to retrieve price,
volume, and orderbook information as well as the ability to trade
crypto-currencies.

Calls to the exchange are categorized as either public, which includes
requests for price, volume, and order book information, and private, which
includes all requests requiring an account including placing buy or sell
orders. Public calls can be used immediately after installing the package.
Private calls require creating an account at \url{https://bittrex.com} and
creating API and secret keys with appropriate permissions.

Private calls retrieve the API and secret key using the BITTREX_API_KEY and
BITTREX_SECRET_KEY environment variables. These may be set by the user
before opening the R session or, they can be set using the
\code{\link[=bt_authenticate]{bt_authenticate()}} function.
}
\section{Public Function Calls}{

\itemize{
\item{\code{\link[=bt_api_check]{bt_api_check()}}  }{check if the bittrex REST API is working}
\item{\code{\link[=bt_getcurrencies]{bt_getcurrencies()}}  }{all supported currencies at Bittrex along with other meta data}
\item{\code{\link[=bt_getmarkethistory]{bt_getmarkethistory()}}  }{the latest trades that have occurred for a specified market}
\item{\code{\link[=bt_getmarkets]{bt_getmarkets()}}  }{the open and available trading markets at Bittrex along with other meta data}
\item{\code{\link[=bt_getmarketsummary]{bt_getmarketsummary()}}  }{the last 24 hours' summary of a specific market}
\item{\code{\link[=bt_getmarketsummaries]{bt_getmarketsummaries()}}  }{the last 24 hours' summary of all active markets}
\item{\code{\link[=bt_getorderbook]{bt_getorderbook()}}  }{the orderbook for a given market}
\item{\code{\link[=bt_getticker]{bt_getticker()}}  }{the current tick values for a market}
}
}

\section{Private Function Calls}{

\itemize{
\item{\code{\link[=bt_authenticate]{bt_authenticate()}}  }{provide user authentication data}
\item{\code{\link[=bt_buy]{bt_buy()}}  }{place a buy limit order}
\item{\code{\link[=bt_cancel]{bt_cancel()}}  }{cancel a buy or sell order}
\item{\code{\link[=bt_getbalances]{bt_getbalances()}}  }{account balances for all currencies}
\item{\code{\link[=bt_getbalance]{bt_getbalance()}}  }{account balance for a specified currency}
\item{\code{\link[=bt_getdepositaddress]{bt_getdepositaddress()}}  }{retrieve or generate an address for a specified
currency}
\item{\code{\link[=bt_getdeposithistory]{bt_getdeposithistory()}}  }{retrieve your deposit history}
\item{\code{\link[=bt_getopenorders]{bt_getopenorders()}}  }{retrieve data for all open orders}
\item{\code{\link[=bt_getorder]{bt_getorder()}}  }{retrieve a single order by uuid}
\item{\code{\link[=bt_getorderhistory]{bt_getorderhistory()}}  }{recent order history for an account }
\item{\code{\link[=bt_getwithdrawalhistory]{bt_getwithdrawalhistory()}}  }{retrieve your withdrawal history}
\item{\code{\link[=bt_sell]{bt_sell()}}  }{place a sell limit order}
\item{\code{\link[=bt_withdraw]{bt_withdraw()}}  }{withdraw funds from your account}
}
}

\references{
\url{https://bittrex.com},
\url{https://github.com/kaneplusplus/bittrex}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/private-interface.r
\name{bt_getdepositaddress}
\alias{bt_getdepositaddress}
\title{Retrieve the Address for a Specified Currency}
\usage{
bt_getdepositaddress(currency)
}
\arguments{
\item{currency}{currency to get the deposit address}
}
\value{
A named list with the following elements:
\itemize{
\item{success: }{a boolean indicating whether the request was successful.}
\item{message: }{a string describing the error if the request was not
successful, otherwise an empty string. If an
address has not been generated this field will have value
"ADDRESS_GENERATING" until it is available.}
\item{result:  }{a \code{data.frame} with the one row and columns
providing the currency and the address.}
}
}
\description{
The \code{bt_getdepositaddress()} retrieves or creates the account
deposit address for a specified currency.
}
\examples{
\dontrun{
# Note you must authenticate first.
bt_getdepositaddress("btc")
# $success
# [1] TRUE
#
# $message
# [1] ""
#
# $result
#   currency                            address
# 1      BTC 1Q6WissSMNF7NCNw3sDXQ2F7AbrSCYouj2
}
}
\references{
https://bittrex.com/api/v1.1/account/getdepositaddress
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/private-interface.r
\name{bt_authenticate}
\alias{bt_authenticate}
\title{Provide User Authentication Data}
\usage{
bt_authenticate(api_key, secret_key)
}
\arguments{
\item{api_key}{the api key provided by the exchange}

\item{secret_key}{the secret key provided by the exchange}
}
\description{
The \code{bt_authenicate()} function sets the
BITTREX_API_KEY and BITTREX_SECRET_KEY environment variables in your current
session to access your account on \url{https://bittrex.com}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public-interface.r
\name{bt_getmarkets}
\alias{bt_getmarkets}
\title{Available Markets and Other Meta Data}
\usage{
bt_getmarkets()
}
\value{
A named list with the following elements:
\itemize{
\item{success: }{a boolean indicating whether the request was successful.}
\item{message: }{a string describing the error if the request was not
successful, otherwise and empty string.}
\item{result:  }{A \code{data.frame} with the market currencies,
base currencies, base currency long name, minimum trade size, market name,
if the market is active, when the market was created, market notices,
if the market is sponsored, and the location of the market logo.}
}
}
\description{
The \code{bt_getmarkets()} function returns all of the currently
available markets on \url{https://bittrex.com} along with other information
including when the exchange was created and the minimum order size.
}
\examples{
\dontrun{
markets <- bt_getmarkets()$result
head(markets)
}
}
\references{
https://bittrex.com/api/v1.1/public/getmarkets
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/private-interface.r
\name{bt_sell}
\alias{bt_sell}
\title{Place a Sell Order}
\usage{
bt_sell(market, quantity, rate, type = c("limit", "market"))
}
\arguments{
\item{market}{the market to place the buy limit order on.}

\item{quantity}{how much of the quote currency to sell.}

\item{rate}{the price per unit of the quote currency that you wish to sell
at.}

\item{type}{either \code{"market"} or \code{"limit"}. Note that market orders are
disabled as of July 7, 2017 (default is limit).}
}
\value{
A named list with the following elements:
\itemize{
\item{success: }{a boolean indicating whether the request is successful.}
\item{message: }{a string describing the error if the request was not
successful, otherwise an empty string.}
\item{result: }{a named list, called "uuid" whose element is an integer
identifying the order. This value is used to query the status of of the
order with either the \code{\link[=bt_getorder]{bt_getorder()}} or \code{\link[=bt_getopenorders]{bt_getopenorders()}} function.
When the order is fulfilled it appears in the order history \code{data.frame}
returned by the \code{\link[=bt_getorderhistory]{bt_getorderhistory()}} function. }
}
}
\description{
The \code{bt_sell()} function places a sell order on
\url{https://bittrex.com}. This function only works if you have set up
authentication.
}
\examples{
\dontrun{
# Note you must authenticate first.
# Sell half a bitcoin at 17000 USDT per bitcoin.
bt_sell("usdt-btc", 0.5, 17000)
# $success
# [1] TRUE
#
# $message
# [1] ""
#
# $result
# $result$uuid
# [1] "2d6l69e9-17fb-4f2a-8aff-37418b515624"

# Sell 3.4 ETH at 0.063 BTC per ETH
bt_sell("btc-eth", 3.4, 0.063)

# Sell 2 LTC at 255 USDT per LTC
bt_sell("usdt-ltc", 2, 255)
}
}
\references{
https://bittrex.com/api/v1.1/market/selllimit
}
\seealso{
\code{\link[=bt_authenticate]{bt_authenticate()}}, \code{\link[=bt_buy]{bt_buy()}}, \code{\link[=bt_getopenorders]{bt_getopenorders()}},
\code{\link[=bt_getorderhistory]{bt_getorderhistory()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/private-interface.r
\name{bt_getopenorders}
\alias{bt_getopenorders}
\title{Order Data for all Open Orders}
\usage{
bt_getopenorders(market)
}
\arguments{
\item{market}{(optional) the market on which you would like to see all
open orders. If not specified, then all open orders
for all markets are returned.}
}
\value{
A named list with the following elements:
\itemize{
\item{success: }{a boolean indicating whether the request was successful.}
\item{message: }{a string describing the error if the request was not
successful, otherwise an empty string.}
\item{result:  }{a \code{data.frame} providing information about the
open orders including (but not limited to) the market, quantity remaining
in the order, the type of order, and when the order was opened.
}
}
}
\description{
The \code{bt_getopenorders()} function retrieves all the user's open
orders on \url{https://bittrex.com}. This function can only be used after
you provide information for authentication.
}
\examples{
\dontrun{
# Note you must authenticate first.
bt_getopenorders()
# $success
# [1] TRUE
#
# $message
# [1] ""
#
# $result
#   uuid                           order_uuid exchange order_type quantity
# 1   NA 2d6169e9-17fb-4f2a-8aff-37418b515624  ETH-LTC  LIMIT_BUY        1
#   quantity_remaining limit commission_paid price price_per_unit
# 1                  1 0.001               0     0             NA
#                   opened closed cancel_initiated immediate_or_cancel
# 1 2017-07-11T18:53:30.07     NA            FALSE               FALSE
#   is_conditional condition condition_target
# 1          FALSE      NONE               NA
}
}
\references{
https://bittrex.com/api/v1.1/market/getopenorders
}
\seealso{
\code{\link[=bt_authenticate]{bt_authenticate()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public-interface.r
\name{bt_getticker}
\alias{bt_getticker}
\title{Get the Ticker Values for a Market}
\usage{
bt_getticker(market)
}
\arguments{
\item{market}{the market to get the ticker for.}
}
\value{
A named list with the following elements:
\itemize{
\item{success: }{a boolean indicating whether the request was successful.}
\item{message: }{a string describing the error if the request was not
successful, otherwise and empty string.}
\item{result:  }{A \code{data.frame} with the bid, ask, and last
transaction prices.}
}
}
\description{
The \code{bt_getticker()} function returns the bid, ask, and last
transaction price for a specified market on \url{https://bittrex.com}.
The complete list of markets is available via \code{\link[=bt_getmarkets]{bt_getmarkets()}}.
}
\examples{
\dontrun{
getticker("btc-ltc")
}
}
\references{
https://bittrex.com/api/v1.1/public/getticker
}
\seealso{
\code{\link[=bt_getmarkets]{bt_getmarkets()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/private-interface.r
\name{bt_getbalance}
\alias{bt_getbalance}
\title{Account Balance for a Specified Currency}
\usage{
bt_getbalance(currency)
}
\arguments{
\item{currency}{currency to retrieve the account balance for.}
}
\value{
A named list with the following elements:
\itemize{
\item{success: }{a boolean indicating whether the request was successful.}
\item{message: }{a string describing the error if the request was not
successful, otherwise an empty string.}
\item{result:  }{a \code{data.frame} with the currency, balance,
available funds, the amount of any pending transactions, and
cryptographic addresses that can be used to receive funding.
}
}
}
\description{
The \code{bt_getbalance()} function retrieves the account balance
for a specified currency on \url{https://bittrex.com}. This function
can only be used after you provide information for authentication.
}
\examples{
\dontrun{
Note you must authenticate first.
bt_getbalance("btc")$result
#   currency balance available pending                     crypto_address
# 1      BTC       0         0       0 1Q6WissSMNF7NCNw3sDXQ2F7AbrSCYouj2
}
}
\references{
https://bittrex.com/api/v1.1/account/getbalances
}
\seealso{
\code{\link[=bt_authenticate]{bt_authenticate()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/private-interface.r
\name{bt_cancel}
\alias{bt_cancel}
\title{Cancel an Open Order}
\usage{
bt_cancel(uuid)
}
\arguments{
\item{uuid}{the uuid of the order you would like to cancel.}
}
\value{
A named list with the following elements:
\itemize{
\item{success: }{a boolean indicating whether the request was successful.}
\item{message: }{a string describing the error if the request was not
successful, otherwise an empty string.}
\item{result:  }{always NULL}
}
}
\description{
The \code{bt_cancel()} function cancels an open order on
\url{https://bittrex.com}. This function
is called after using either \code{\link[=bt_buy]{bt_buy()}} or \code{\link[=bt_sell]{bt_sell()}}.
}
\examples{
\dontrun{
# Note you must authenticate and define a uuid first.
bt_cancel(uuid)
# $success
# [1] TRUE
#
# $message
# [1] ""
#
# $result
}
}
\references{
https://bittrex.com/api/v1.1/account/cancel
}
\seealso{
\code{\link[=bt_authenticate]{bt_authenticate()}}, \code{\link[=bt_getopenorders]{bt_getopenorders()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public-interface.r
\name{bt_getorderbook}
\alias{bt_getorderbook}
\title{Order Book for a Market}
\usage{
bt_getorderbook(market, type = c("both", "buy", "sell"), depth = 50)
}
\arguments{
\item{market}{the market from which the order book will be retrieved.}

\item{type}{type of orders to retrieve (default is "both")}

\item{depth}{how deep the returned order book should be (default and
maximum are 50). This is the size and price of bids whose price is
lower than the highest bid and higher than the lowest ask.}
}
\value{
A named list with the following elements:
\itemize{
\item{success: }{a boolean indicating whether the request was successful.}
\item{message: }{a string describing the error if the request was not
successful, otherwise and empty string.}
\item{result:  }{A named list with the buy and sell orders (depending
on the specified \code{type} parameter). If \code{type} is "buy" or
"both" then the list will contain a element named "buy" with
a \code{data.frame} of the buy orders.}
}
}
\description{
The \code{bt_getorderbook()} function returns the order book
for a specified market on \url{https://bittrex.com}.
}
\examples{
\dontrun{
head(bt_getorderbook("btc-eth")$result)
}
}
\references{
https://bittrex.com/api/v1.1/public/getorderbook
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/private-interface.r
\name{bt_withdraw}
\alias{bt_withdraw}
\title{Withdraw Funds from an Account}
\usage{
bt_withdraw(currency, quantity, address, paymentid)
}
\arguments{
\item{currency}{the currency to withdraw.}

\item{quantity}{the quantity of the currency to withdraw.}

\item{address}{where to send the funds.}

\item{paymentid}{CryptoNotes/BitShareX/Nxt optional field (memo/paymentid).}
}
\value{
A named list with the following elements:
\itemize{
\item{success: }{a boolean indicating whether the request was successful.}
\item{message: }{a string describing the error if the request was not
successful, otherwise an empty string.}
\item{result:  }{a named list, with element "uuid" whose element is an
integer identifying the order.}
}
}
\description{
The \code{bt_withdraw()} function moves funds from a
\url{https://bittrex.com} account to a specified address. It does not
include the transaction fee.
}
\examples{
\dontrun{
# Note you must authenticate first.
# Send the author your bitcoins.
bt_widthdraw("btc", 10, "1Q6WissSMNF7NCNw3sDXQ2F7AbrSCYouj2")
}
}
\references{
https://bittrex.com/api/v1.1/account/withdraw
}
\seealso{
\code{\link[=bt_getbalance]{bt_getbalance()}}, \code{\link[=bt_getbalances]{bt_getbalances()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/private-interface.r
\name{bt_getbalances}
\alias{bt_getbalances}
\title{Account Balances for All Currencies}
\usage{
bt_getbalances()
}
\value{
A named list with the following elements:
\itemize{
\item{success: }{a boolean indicating whether the request was successful.}
\item{message: }{a string describing the error if the request was not
successful, otherwise an empty string.}
\item{result:  }{a \code{data.frame} with the currencies, balances,
available funds, the amount of any pending transactions, and
cryptographic addresses that can be used to receive funding.
}
}
}
\description{
The \code{bt_getbalances()} function retrieves the account balance
for all currencies on \url{https://bittrex.com}. This function
can only be used after you provide information for authentication.
}
\examples{
\dontrun{
# Note you must authenticate first.
balances <- bt_getbalances()$result
# $success
# [1] TRUE
#
# $message
# [1] ""
#
# $result
#   currency   balance available pending
# 1      BTC 0.0000000 0.0000000       0
# 2      ETH 0.2187638 0.2187638       0
# 3      LTC 0.0000000 0.0000000       0
#                               crypto_address
# 1         1Q6WissSMNF7NCNw3sDXQ2F7AbrSCYouj2
# 2 0x0ceac821a72037b07df691a53e201d797252b5a6
# 3         Li71CUBjxFH6PfEZn2phqfPhoasydfNfqF
}
}
\references{
https://bittrex.com/api/v1.1/account/getbalances
}
\seealso{
\code{\link[=bt_authenticate]{bt_authenticate()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public-interface.r
\name{bt_getcurrencies}
\alias{bt_getcurrencies}
\title{Retrieve all Available Currencies on the Exchange}
\usage{
bt_getcurrencies()
}
\value{
A named list with the following elements:
\itemize{
\item{success: }{a boolean indicating whether the request was successful.}
\item{message: }{a string describing the error if the request was not
successful, otherwise and empty string.}
\item{result:  }{A \code{data.frame} with the currency ticker,
currency, a minimum confirmation number, the transaction fee, if the
currency is active, the coin type, the base address, and currency
notices.}
}
}
\description{
The \code{bt_getcurrencies()} function returns the available
currencies on \url{https://bittrex.com}.
}
\examples{
\dontrun{
currencies <- bt_getcurrencies()$result
}
}
\references{
https://bittrex.com/api/v1.1/public/getcurrencies
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public-interface.r
\name{bt_getmarketsummaries}
\alias{bt_getmarketsummaries}
\title{Summary of All Active Markets}
\usage{
bt_getmarketsummaries()
}
\value{
A named list with the following elements:
\itemize{
\item{success: }{a boolean indicating whether the request was successful.}
\item{message: }{a string describing the error if the request was not
successful, otherwise and empty string.}
\item{result:  }{A \code{data.frame} with one row per market and, for
each market: the market name, the high, the low, the
volume, the last trade, the last trade price, the
base currency volume, a time stamp for the last
transaction, the current bid, the current ask, the number
of open buy orders, the number of open sell orders, the
the previous day close, and when the market was created.
}
}
}
\description{
the \code{bt_getmarketsummaries()} retrieves a summary of all
active markets on \url{https://bittrex.com} for the last 24 hours.
}
\examples{
\dontrun{
ms <- bt_getmarketsummaries()$result
head(ms)
}
}
\references{
https://bittrex.com/api/v1.1/public/getmarketsummaries
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/private-interface.r
\name{bt_getdeposithistory}
\alias{bt_getdeposithistory}
\title{Retrieve Deposit History}
\usage{
bt_getdeposithistory(currency)
}
\arguments{
\item{currency}{(optional) the currency to retrieve the deposits for. If
this is not specified then deposit history for all currencies is retrieved.}
}
\value{
A named list with the following elements:
\itemize{
\item{success: }{a boolean indicating whether the request was successful.}
\item{message: }{a string describing the error if the request was not
successful, otherwise an empty string.}
\item{result:  }{a \code{data.frame} providing data about
previously completed orders including the order uuid, the currency,
the time of the withdraw, the quantity, etc.
}
}
}
\description{
The \code{bt_getdeposithistory()} function retrieves the
deposit history for an account on \url{https://bittrex.com}.
This function can only be used after you provide authentication information.
}
\examples{
\dontrun{
bt_getdeposithistory()
# $success
# [1] TRUE
#
# $message
# [1] ""
#
# $result
#         id     amount currency confirmations        last_updated
# 1 20774372 0.39125728      ETH            49 2017-06-22 16:05:50
# 2 18255803 0.05936286      BTC             6 2017-05-19 16:28:36
#                                                                tx_id
# 1 0xbecc44384d8b94f1d03834ffb9324e97e4fa2a8161e17e61116aaabd5fb35050
# 2   7084cad99373475d8137547ce947b1472bfcb2d23b5160b05010f1f15e3c6287
#                               crypto_address
# 1 0x0ceac821a72037b07df691a53e201d797252b5a6
# 2         1Q6WissSMNF7NCNw3sDXQ2F7AbrSCYouj2
}
}
\references{
https://bittrex.com/api/v1.1/account/getdeposithistory
}
\seealso{
\code{\link[=bt_authenticate]{bt_authenticate()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/private-interface.r
\name{bt_getorderhistory}
\alias{bt_getorderhistory}
\title{Order History for an Account}
\usage{
bt_getorderhistory(market)
}
\arguments{
\item{market}{(optional) the market on which you would like to see all
open orders. If not specified, then completed orders for all markets are
returned.}
}
\value{
A named list with the following elements:
\itemize{
\item{success: }{a boolean indicating whether the request was successful.}
\item{message: }{a string describing the error if the request was not
successful, otherwise an empty string.}
\item{result:  }{a \code{data.frame} providing data about
previously completed orders including the order uuid, the exchange
the time of the order, the order type, the limit, the quantity, the
quantity remaining, the commission, the price, the price per unit,
and whether or not it was a conditional trade.
}
}
}
\description{
The \code{bt_getorderhistory()} function retrieves order history
data on \url{https://bittrex.com}.
This function can only be used after you provide authentication information.
}
\examples{
\dontrun{
bt_getorderhistory()
# $success
# [1] TRUE
#
# $message
# [1] ""
#
# $result
#                             order_uuid exchange time_stamp order_type   limit
# 1 c04bc07b-e6a9-4f47-a2c8-f9eb3c9a7fa1  BTC-ETH       <NA> LIMIT_SELL 0.11771
# 2 319fcc92-0b0d-43f1-a538-18f790d85ffa  BTC-ETH       <NA> LIMIT_SELL 0.14170
# 3 cd052594-d655-4e79-bec6-383ba7be8302  BTC-ETH       <NA> LIMIT_SELL 0.14020
#   quantity quantity_remaining commission      price price_per_unit
# 1    0.400                  0 0.00011772 0.04708801      0.1177200
# 2    0.100                  0 0.00003542 0.01417000      0.1417000
# 3    0.175                  0 0.00006133 0.02453499      0.1401999
#   is_conditional condition condition_target immediate_or_cancel
# 1          FALSE      NONE               NA               FALSE
# 2          FALSE      NONE               NA               FALSE
# 3          FALSE      NONE               NA               FALSE
#                    closed
# 1 2017-06-22T20:06:23.973
# 2   2017-06-13T15:33:20.4
# 3 2017-06-13T14:59:13.923
}
}
\references{
https://bittrex.com/api/v1.1/account/getorderhistory
}
\seealso{
\code{\link[=bt_authenticate]{bt_authenticate()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public-interface.r
\name{bt_getmarkethistory}
\alias{bt_getmarkethistory}
\title{Recent History for a Market}
\usage{
bt_getmarkethistory(market)
}
\arguments{
\item{market}{the market from which history data will be retrieved.}
}
\value{
A named list with the following elements:
\itemize{
\item{success: }{a boolean indicating whether the request was successful.}
\item{message: }{a string describing the error if the request was not
successful, otherwise and empty string.}
\item{result:  }{A \code{data.frame} containing recent trade information
including the order type, time, quantity, price, and fill type.}
}
}
\description{
the \code{bt_getmarkethistory()} function retrieves recent trade
information for a specified market on \url{https://bittrex.com}.
}
\examples{
\dontrun{
mh <- bt_getmarkethistory("btc-eth")$result
head(mh)
}
}
\references{
https://bittrex.com/api/v1.1/public/getmarkethistory
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/private-interface.r
\name{bt_buy}
\alias{bt_buy}
\title{Place a Buy Limit Order}
\usage{
bt_buy(market, quantity, rate, type = c("limit", "market"))
}
\arguments{
\item{market}{the market to place the buy limit order on.}

\item{quantity}{how much of the quote currency you want to buy.}

\item{rate}{the price per unit of the quote currency that you wish to buy at.}

\item{type}{either "market" or "limit". Note that market orders are
disabled as of July 7, 2017. (default is "limit")}
}
\value{
A named list with the following elements:
\itemize{
\item{success: }{a boolean indicating whether the request is successful.}
\item{message: }{a string describing the error if the request was not
successful, otherwise an empty string.}
\item{result:  }{a named list, with element "uuid" whose element is an
integer identifying the order. This value is used to query the status of
of the order with either the \code{\link[=bt_getorder]{bt_getorder()}} or \code{\link[=bt_getopenorders]{bt_getopenorders()}}
function. When the order is fulfilled it appears in the order history
\code{data.frame} returned by the \code{\link[=bt_getorderhistory]{bt_getorderhistory()}} function.}
}
}
\description{
The \code{bt_buy()} function places a buy order on
\url{https://bittrex.com}. This function only works after you have set up
authentication.
}
\examples{
\dontrun{
# Note you must authenticate first.
# Buy half a bitcoin at 17000 USDT per bitcoin.
bt_buy("usdt-btc", 0.5, 17000)
# $success
# [1] TRUE
#
# $message
# [1] ""
#
# $result
# $result$uuid
# [1] "2d6169e9-17fb-4f2a-8aff-37418b515624"

# Buy 3.4 ETH at 0.063 BTC per ETH
bt_buy("btc-eth", 3.4, 0.063)

# Buy 2 LTC at 255 USDT per LTC
bt_buy("usdt-ltc", 2, 255)
}
}
\references{
https://bittrex.com/api/v1.1/market/buylimit
}
\seealso{
\code{\link[=bt_authenticate]{bt_authenticate()}}, \code{\link[=bt_sell]{bt_sell()}}, \code{\link[=bt_getorder]{bt_getorder()}},
\code{\link[=bt_getopenorders]{bt_getopenorders()}}, \code{\link[=bt_getorderhistory]{bt_getorderhistory()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/private-interface.r
\name{bt_getwithdrawalhistory}
\alias{bt_getwithdrawalhistory}
\title{Retrieve Withdrawal History}
\usage{
bt_getwithdrawalhistory(currency)
}
\arguments{
\item{currency}{(optional) the currency to retrieve the withdraw for. If
this is not specified then withdraw history for all currencies is
retrieved.}
}
\value{
A named list with the following elements:
\itemize{
\item{success: }{a boolean indicating whether the request was successful.}
\item{message: }{a string describing the error if the request was not
successful, otherwise an empty string.}
\item{result:  }{a \code{data.frame} providing data about
previously completed orders including the order uuid, the currency,
the time of the withdraw, the quantity, etc.
}
}
}
\description{
The \code{bt_getwithdrawalhistory()} function retrieves the
withdraw history for an account on \url{https://bittrex.com}. This function
can only be used after you provide authentication information.
}
\examples{
\dontrun{
# Note you must authenticate first.
bt_getwithdrawalhistory()
# $success
# [1] TRUE
#
# $message
# [1] ""
#
# $result
#                           payment_uuid currency     amount
# 1 ba0c85de-1fd8-423e-939d-e34d2aad34fd      BTC 0.04597029
# 2 f10c3536-fcf2-48eb-9ce4-253271d2c8e8      BTC 0.01313458
# 3 5dae07ad-7a8f-40a0-ac6e-225a3d0d6d8a      BTC 0.02347405
#                              address              opened authorized
# 1 1C31WQL12CDZqnidra9tMfs4DLkebAuNgc 2017-06-22 20:08:26       TRUE
# 2 1C31WQL12CDZqnidra9tMfs4DLkebAuNgc 2017-06-13 15:34:47       TRUE
# 3 1C31WQL12CDZqnidra9tMfs4DLkebAuNgc 2017-06-13 15:03:23       TRUE
#   pending_payment tx_cost
# 1           FALSE   0.001
# 2           FALSE   0.001
# 3           FALSE   0.001
#                                                              tx_id canceled
# 1 e628848ed92be4baee877f97e3a48b22f5ee2f7ca35c2908282b8c9ee2f4b94a    FALSE
# 2 fbafe847d02761d089b19a2cafecff561030219ded1eb03cc796c8c2eac0dd5c    FALSE
# 3 c981f7dc569188db16753cff4ab24aef148039964b68428603e2bfd18c754df4    FALSE
#   invalid_address
# 1           FALSE
# 2           FALSE
# 3           FALSE
}
}
\references{
https://bittrex.com/api/v1.1/account/getwithdrawalhistory
}
\seealso{
\code{\link[=bt_authenticate]{bt_authenticate()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public-interface.r
\name{bt_getmarketsummary}
\alias{bt_getmarketsummary}
\title{Summary of a Markets}
\usage{
bt_getmarketsummary(market)
}
\arguments{
\item{market}{the market to retrieve the summary for.}
}
\value{
A named list with the following elements:
\itemize{
\item{success: }{a boolean indicating whether the request was successful.}
\item{message: }{a string describing the error if the request was not
successful, otherwise and empty string.}
\item{result:  }{A \code{data.frame} with one row and columns corresponding
to: the market name, the high, the low, the
volume, the last trade, the last trade price, the
base currency volume, a time stamp for the last
transaction, the current bid, the current ask, the number
of open buy orders, the number of open sell orders, the
the previous day close, and when the market was created.
}
}
}
\description{
the \code{bt_getmarketsummary()} retrieves a summary for a specified
market on \url{https://bittrex.com}.
}
\examples{
\dontrun{
ms <- bt_getmarketsummary("btc-eth")$result
head(ms)
}
}
\references{
https://bittrex.com/api/v1.1/public/getmarketsummary
}
