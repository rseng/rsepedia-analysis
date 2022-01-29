---
title: 'Finnish Media Scrapers'
tags:
  - Web scraping
  - Media research
  - Python
authors:
  - name: Eetu Mäkelä^[corresponding author]
    orcid: 0000-0002-8366-8414
    affiliation: 1
  - name: Pihla Toivanen^[co-first author]
    orcid: 0000-0003-0872-7098
    affiliation: 1
affiliations:
 - name: University of Helsinki
   index: 1
date: 21 June 2021
bibliography: paper.bib

---

# Summary

Finnish Media Scrapers is a package for extracting articles from Finnish journalistic media websites.
Included are scrapers for the four biggest Finnish journalistic media: [YLE](https://www.yle.fi/uutiset/), [Helsingin Sanomat](https://www.hs.fi/), [Iltalehti](https://www.iltalehti.fi/) and [Iltasanomat](https://www.is.fi/). The scrapers have been designed for researchers needing a local corpus of news article texts matching a specified set of query keywords as well as temporal limitations. As a design principle, these scrapers have been designed to extract the articles in as trustworthy a manner as possible, as required for content-focused research targetting the text of those articles. Further, the process is split into distinct parts that 1) query, 2) fetch, 3) convert to text and 4) post-filter the articles separately. Each of these steps also records its output as separate files. This way, the tools can be used in a versatile manner. Further, a good record is maintained of the querying and filtering process for reproducibility as well as error analysis.

# Statement of need

There is an increasing need for user-friendly computational tools in the humanities and social sciences. For example, a common workflow in media research is to collect a large amount of data and combine quantitative and qualitative methods in the analysis phase (@Koivunen, @Weber2011). This package responds to the research needs by providing easy-to-use tools for scraping Finnish media articles and extracting the article texts from the scraped HTML files. At the same time, the functionality has also been packaged as a Python module for the benefit of more computationally-savvy users.

The scripts were originally developed for a data journalism article (@SKarticle) analyzing how Finnish members of parliament were represented in the media in 2020. Further developing and packaging the scripts into a reusable package was based on an expressed interest from the Finnish computational science community. Since initial beta release a couple of months ago, the package is now known to be already used in at least two research projects targeting Finnish media analysis.

## Related work

For a more general library for crawling media articles, have a look at [newspaper3k](https://newspaper.readthedocs.io/en/latest/index.html) as well as [news-please](https://github.com/fhamborg/news-please), which has been built on top of it. Do note however that at the time of writing this, it is [unclear](https://github.com/codelucas/newspaper/issues/878) whether newspaper3k is being maintained any more. More importantly for content research purposes, note that 1) newspaper3k does not handle the Finnish news sources targeted by this crawler very well and 2) it is based more on a best-effort principle (suitable for extracting masses of data for e.g. NLP training) as opposed to completeness and verisimilitude (required for trustworthy content-focused research targetting a particular set of news). Thus, given an article URL, newspaper3k will happily try to return something from it, but not guarantee completeness. This crawler on the other hand has been designed to be conservative, and to complain loudly through logging whenever it encounters problems that may hinder extracting the actual text of the article, such as article layouts that haven't been yet handled and verified to extract correctly.

# General workflow

The general workflow for using the scrapers is as follows:
1. The scrapers support specifying a keyword as well as a timespan for extraction, and output a CSV of all matching articles with links.
2. A second set of scripts then allows downloading the matched articles in HTML format.
3. Third, there are further scripts for extracting plain text versions of the article texts out of the HTML.
4. Finally, a script exists to post-filter the resulting plain texts again with keywords.

Important to know when applying the workflow is that due to the fact that all the sources use some kind of stemming for their search, they can often return also spurious hits. Further, if searching for multiple words, the engines often perform a search for either word instead of the complete phrase. The post-filtering script above exists to counteract this by allowing the refiltering of the results more rigorously and uniformly locally.

At the same time and equally importantly, the stemming for a particular media may not cover e.g. all inflectional forms of words. Thus, it often makes sense to query for at least all common inflected variants and merge the results. For a complete worked up example of this kind of use, see the [members_of_parliament](https://github.com/hsci-r/finnish-media-scraper/tree/master/members_of_parliament) folder, which demonstrates how one can collect and count how many articles in each media mention the members of the Finnish Parliament.

# Acknowledgements

We acknowledge contributions from the Suomen Kuvalehti team (Samuel Nyroos, Salla Vuorikoski and Leena Sharma) during the testing phase of the scrapers.

# References
# Finnish Media Scrapers

[![PyPI version](https://badge.fury.io/py/finnish-media-scrapers.svg)](https://badge.fury.io/py/finnish-media-scrapers) [![DOI](https://zenodo.org/badge/335605978.svg)](https://zenodo.org/badge/latestdoi/335605978) [![Documentation Status](https://readthedocs.org/projects/finnish-media-scrapers/badge/?version=latest)](https://finnish-media-scrapers.readthedocs.io/en/latest/?badge=latest)

Scrapers for extracting articles from Finnish journalistic media websites by the [University of Helsinki](https://www.helsinki.fi/) [Human Sciences – Computing Interaction research group](https://heldig.fi/hsci/). Included are scrapers for [YLE](https://www.yle.fi/uutiset/), [Helsingin Sanomat](https://www.hs.fi/), [Iltalehti](https://www.iltalehti.fi/) and [Iltasanomat](https://www.is.fi/).

The scrapers have been designed for researchers needing a local corpus of news article texts matching a specified set of query keywords as well as temporal limitations. As a design principle, these scrapers have been designed to extract the articles in as trustworthy a manner as possible, as required for content-focused research targetting the text of those articles (for an example of such research, see e.g. [here](https://researchportal.helsinki.fi/en/publications/a-year-in-the-spotlight-who-got-the-attention-of-the-media-who-wa)). Thus, the scrapers will complain loudly for example if your search query matches more articles than the APIs are willing to return, or if the plain text extractors encounter new article layouts that have not yet been verified to extract correctly. Further, the process is split into distinct parts that 1) query, 2) fetch, 3) convert to text and 4) post-filter the articles separately. Each of these steps also records its output as separate files. Each of these steps also records its output as separate files. This way, the tools can be used in a versatile manner. Further, a good record is maintained of the querying and filtering process for reproducibility as well as error analysis.

## Installation

Install the scripts (and Python module) using `pip install finnish-media-scrapers`. After this, the scripts should be useable from the command line, and the functionality importable from Python. Or, if you have [pipx](https://pypa.github.io/pipx/) and just want the command line scripts, use `pipx install finnish-media-scrapers` instead.

## General workflow

![Data collection workflow](https://github.com/hsci-r/finnish_media_scrapers/raw/master/images/fms_datacollection_50border.png)

The general workflow for using the scrapers is as follows:

1. Query YLE/HS/IL/IS APIs for matching articles using the scripts `fms-query-{yle|hs|il|is}`, which output all matching articles with links into CSVs.
2. Fetch the matching articles using `fms-fetch-{hs|open}`. These save the articles as HTML files in a specified directory.
3. Extract the plain text from the article HMTL using `fms-html-to-text-{yle|svyle|hs|il|is}`.
4. Optionally refilter the results using `fms-post-filter`.

Important to know when applying the workflow is that due to the fact that all the sources use some kind of stemming for their search, they can often return also spurious hits. Further, if searching for multiple words, the engines often perform a search for either word instead of the complete phrase. The post-filtering script above exists to counteract this by allowing the refiltering of the results more rigorously and uniformly locally.

At the same time and equally importantly, the stemming for a particular media may not cover e.g. all inflectional forms of words. Thus, it often makes sense to query for at least all common inflected variants and merge the results. For a complete worked up example of this kind of use, see the [members_of_parliament](https://github.com/hsci-r/finnish-media-scraper/tree/master/members_of_parliament) folder, which demonstrates how one can collect and count how many articles in each media mention chairperson of National Coalition Party (Petteri Orpo) or alternatively all members of the Finnish Parliament.

To be a good netizen, when using the scripts, by default there is a one second delay between each web request to the media websites to ensure that scraping will not cause undue load on their servers. This is however configurable using command line parameters.

Apart from using the scripts, the functionality of the package is also provided as a python module that you may use programmatically from within Python. For the functionalities thus provided, see the [module documentation](https://finnish-media-scrapers.readthedocs.io/en/latest/).

## Media-specific instructions and caveats

### Helsingin Sanomat

First, query the articles you want using `fms-query-hs`. For example, `fms-query-hs -f 2020-02-16 -t 2020-02-18 -o hs-sdp.csv -q SDP`.

For downloading articles, use `fms-fetch-hs` with adding credentials. For example `fms-fetch-hs -i hs-sdp.csv -o hs-sdp -u username -p password`. This scraper requires paid Helsingin Sanomat credentials (user id and password). You can create them in [https://www.hs.fi/](https://www.hs.fi/) with clicking "Kirjaudu" button and following the instructions for a news subscription.

Technically, the scraper uses [pyppeteer](https://pypi.org/project/pyppeteer/) to control a headless Chromium browser to log in and ensure the dynamically rendered content in HS articles is captured. To ensure a compatible Chromium, when first running the tool, pyppeteer will download an isolated version of Chromium for itself, causing some ~150MB of network traffic and disk space usage.

After fetching the articles, extract texts with e.g. `fms-html-to-text-hs -o hs-sdp-output hs-sdp`.

Known special considerations:

- The search engine used seems to be employing some sort of stemming/lemmatization, so e.g. the query string `kok` seems to match `kokki`, `koko` and `koki`.
- A single query can return at most 9,950 hits. This can be sidestepped by invoking the script multiple times with smaller query time spans.

### Yle

example: `fms-query-yle -f 2020-02-16 -t 2020-02-18 -o yle-sdp.csv -q SDP` + `fms-fetch-open -i yle-sdp.csv -o yle-sdp` + `fms-html-to-text-yle -o yle-sdp-output yle-sdp` (or `fms-html-to-text-svyle -o svyle-sdp-output svyle-sdp` if articles come from Svenska YLE)

Known special considerations:

- A single query can return at most 10,000 hits. This can be sidestepped by invoking the script multiple times with smaller query time spans.

### Iltalehti

example: `fms-query-il -f 2020-02-16 -t 2020-02-18 -o il-sdp.csv -q SDP` + `fms-fetch-open -i il-sdp.csv -o il-sdp` + `fms-html-to-text-il -o il-sdp-output il-sdp`

### Iltasanomat

example: `fms-query-is -f 2020-02-16 -t 2020-02-18 -o is-sdp.csv -q SDP` + `fms-fetch-open -i is-sdp.csv -o is-sdp` + `fms-html-to-text-is -o is-sdp-output is-sdp`

Known special considerations:

- The search engine used seems to be employing some sort of stemming/lemmatization, so e.g. the query string `kok` seems to match `kokki`, `koko` and `koki`.
- A single query can return at most 9,950 hits. This can be sidestepped by invoking the script multiple times with smaller query time spans.

### Using the fms-post-filter script

For example, after collecting texts from Helsingin Sanomat with the example above, run:
`fms-post-filter -i hs-sdp.csv -t hs-sdp-output/ -o hs-sdp-filtered.csv -q SDP`

where `-i` parameter specifies the query output file, `-t` the folder name to search extracted texts, `-o` the output filename and `-q` search word to filter.

There is also an option `-ci` for configuring the case-insensitiveness (default false).

## Contact

For more information on the scrapers, please contact associate professor [Eetu Mäkelä](http://iki.fi/eetu.makela). For support on using them or for reporting problems or issues, we suggest you to use the facilities provided by GitHub.

## Development

Pull requests welcome! To set up a development environment, you need [poetry](https://python-poetry.org/). Then, use poetry to install and manage the dependencies and build process (`poetry install`).

## Citation 

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03504/status.svg)](https://doi.org/10.21105/joss.03504)

```
@article{Mäkelä2021,
  doi = {10.21105/joss.03504},
  url = {https://doi.org/10.21105/joss.03504},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {68},
  pages = {3504},
  author = {Eetu Mäkelä and Pihla Toivanen},
  title = {Finnish Media Scrapers},
  journal = {Journal of Open Source Software}
}
```

## Related work

For a more general library for crawling media articles, have a look at [newspaper3k](https://newspaper.readthedocs.io/en/latest/index.html) as well as [news-please](https://github.com/fhamborg/news-please), which has been built on top of it. Do note however that at the time of writing this, it is [unclear](https://github.com/codelucas/newspaper/issues/878) whether newspaper3k is being maintained any more. More importantly for content research purposes, note that 1) newspaper3k does not handle the Finnish news sources targeted by this crawler very well and 2) it is based more on a best-effort principle (suitable for extracting masses of data for e.g. NLP training) as opposed to completeness and verisimilitude (required for trustworthy content-focused research targetting a particular set of news). Thus, given an article URL, newspaper3k will happily try to return something from it, but not guarantee completeness. This crawler on the other hand has been designed to be conservative, and to complain loudly through logging whenever it encounters problems that may hinder extracting the actual text of the article, such as article layouts that haven't been yet handled and verified to extract correctly.
# Example: Petteri Orpo (or all members of parliament)

These example scripts collect articles about Petteri Orpo, chairperson of National Coalition Party (mp-expanded-orpo.csv) during a specific timeframe (from 2020-01-01 to 2020-01-07).

If you want to collect articles about every member of parliament in Finland during 2020, change mp-expanded-orpo.csv to mp-expanded.csv in files fetch-mp-inflections.py and query-mp-inflections.py

`mkdir ../queries`

`python3 query-mp-inflections.py`

`python3 fetch-mp-inflections.py -u <insert your user> -p  <insert your password>`
