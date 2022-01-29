# Lysander: Finding political journalists on Twitter
[![DOI](https://zenodo.org/badge/87090173.svg)](https://zenodo.org/badge/latestdoi/87090173)
[![Research Software Directory](https://img.shields.io/badge/rsd-Research%20Software%20Directory-00a3e3.svg)](https://research-software.nl/software/finding-journalists)

This directory contains software developed in the pilot of
the project [Automated Analysis of Online Behaviour on
Social
Media](https://www.esciencecenter.nl/project/automated-analysis-of-online-behaviour-on-social-media),
a cooperation of the University of Groningen and the
Netherlands eScience Center. The main project software
repository is called [machine
learning](https://github.com/online-behaviour/machine-learning).

The goal of the project is to analyze tweets of politicians
and political journalists. Names of relevant politicians can
be collected from documents like parliament member lists
and ballots but it is much harder to find the names of
relevant journalists. The software in this directory aims to
find such journalists by examining the follower links
between politicians and other users on Twitter.

## Usage

Run like:

```
python getFollowers.py markrutte sybrandbuma apechtold > getFollowers.out 
```

to collect the people that are followed by the users in
your seed list. The script needs your Twitter account data
to be stored in a file definitions.py in the format:

```
# twitter.com authentication keys
token = "???"
token_secret = "???"
consumer_key = "???"
consumer_secret = "???"
```

Replace the strings "???" with the key information from
https://apps.twitter.com , see
https://www.slickremix.com/docs/how-to-get-api-keys-and-tokens-for-twitter/
for instructions

In order to find more relevant users, you can run this command 
after getFollowers.py is finished:

```
python makevec.py markrutte sybrandbuma apechtold < getFollowers.out > makevec.out 2> makevec.err
```

It generates a selection of relevant users (makevec.err) and 
a vector representations for these users (makevec.out)

Input data files for Dutch politicians can be requested
from Erik Tjong Kim Sang e.tjong.kim.sang(at)esciencenter.nl

