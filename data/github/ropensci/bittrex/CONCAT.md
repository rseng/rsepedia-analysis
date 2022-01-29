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
  
