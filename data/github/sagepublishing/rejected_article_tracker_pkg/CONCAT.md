
# SAGE Rejected Article Tracker: Check if an article has been published elsewhere.

## What is it?

The SAGE rejected article tracker
This package will take metadata(title, authors etc...) for research articles ("query articles") and will find the published versions in CrossRef ("match articles"). 

The published versions are selected using very simple machine-learning. 
1. Numerical values are calculated for the difference in the titles of each query article and each candidate match article returned by CrossRef.
2. ArXiv preprints with known DOIs are used to train a simple logistic regression to find the correct candidate for each query. 

This process is approximately 95% accurate. However, the package includes code to recreate the training dataset using the ArXiv and CrossRef APIs. This training dataset can be used to test other approaches to this problem and similar problems (such as detecting duplicate articles).

Presently, the package takes JSON input, but using a library such as [pandas](https://github.com/pandas-dev/pandas) to convert it, spreadsheet data can be used.  

If you are a ScholarOne user, you can produce this easily. More detailed instructions are below. 

## Where to get it?

Although the source code is on Github, the latest release is always published to the [Python package index](https://pypi.org/project/rejected-article-tracker).

```
pip install rejected_article_tracker
```

## Dependencies
- [pandas](https://pypi.org/project/pandas)
- [openpyxl](https://pypi.org/project/openpyxl)
- [xlsxwriter](https://pypi.org/project/xlsxwriter)
- [fuzzywuzzy](https://pypi.org/project/fuzzywuzzy)
- [requests](https://pypi.org/project/requests)
- [sklearn](https://pypi.org/project/sklearn)
- [numpy](https://pypi.org/project/numpy)


## Typical Usage

```python
# Import package
from rejected_article_tracker import RejectedArticlesMatch

# Create a list of article data dicts.
# PLEASE NOTE: Due to reliance on 3rd party APIs, the more articles the longer it takes.
 
articles = [ # Fake values fabricated for the example
{
      "journal_name": "FakeJ2",
      "manuscript_title": "Multiplex Genome Engineering Using CRISPR/Cas Systems",
      "authors": "Cong, Le; Ran, F. Ann; Cox, David",
      "final_decision": "",
      "decision_date": "2012-02-15T13:29:58.999Z", 
      "submission_date": "2012-02-14T13:29:58.999Z",
      "manuscript_id": "ABC-18-070",
    }
]

# @see below for configuration details.
config = {
    "threshold": 70, # Filters out matches which have a fuzz.ratio below this value (fuzz.ratio is a normalised form of Levenshtein distance)
}


# The CrossRef API requires an email address for lookups.    
email = "someome@example.com"


# Define a 'results' list. 
# This is a little unconventional and not great practice. 
# However, injecting a results variable, allows us grab already processed articles
# in case of an issue with 3rd party APIs. 
# This would then mean only the remaining articles would need to be rerun.
results = []

# Run match
RejectedArticlesMatch(
            articles=articles,
            config=config,
            email=email,
            results=results
        ).match()

print(results)
```
---
### Article format 
This is how each article should be structured:

| Property | Contents |
| --- | --- |
| `journal_name` | Simply the name of the journal e.g. ‘Proceedings of mysociety’ |
| `manuscript_title` | The text of the title of the manuscript you are looking for. The article will not be found without this. |
| `authors` | This is a critical part of the input data and MUST be formatted ‘LastName1, FirstName1; LastName2, FirstName2 …’ So note that authors names appear with their lastname first and then a comma, followed by the first name. Then, when we have multiple authors on a paper, we separate them with semicolons. |
| `final_decision` | E.g. ‘Accept’ or ‘Reject’ (the application will simply ignore everything that is listed as ‘Accept’) |
| `decision_date` | The date of final rejection from your journal. Should be a date string format, like: “YYYY-MM-DD”, “YYYY-MM-DD hh:mm:ss” Either format is fine. |
| `submission_date` | Same format as decision date. |
| `manuscript_id` | This can be any unique id that you use to identify the article. If it is formatted in the way that ScholarOne format their ids, then revision numbers will be ignored so that each article is only searched for once. So ‘ABC-20-123’, ‘ABC-20-123.R1’ and ‘ABC-20-123.R2’ all become ‘ABC-20-123’ |

---
<br>
<br>

## Usage with directly downloaded Scholar One data
**Preparing the input:**

If using Scholar One, under `‘Peer Review Details Reports’` select `‘Build Your Own Reports’`. 
The report should have the following columns:

```
‘Journal Name', 'Manuscript ID', 'Manuscript Title', 'Author Names', 'Submission Date', 'Decision Date', 'Accept or Reject Final Decision'
```

It’s ok to include other columns in the Excel file, but they are not needed. 
_IMPORTANTLY:_ remember to download your report as Excel 2007 Data Format.

If not using ScholarOne, then simply prepare your submission data with the above column headings ensuring that
- Author lists are correctly formatted `last_name1, first_name1; last_name2, first_name2;...`
- Date columns (`Submission Date` and `Decision Date`) are correctly formatted, either as datetime strings (as below) or as simple date strings, e.g. `2010-12-31` for 31 December 2010.

**Usage:**
```python

from rejected_article_tracker import ScholarOneRejectedArticlesMatch
import pandas as pd

df = pd.read_excel("/path/to/file")
allowed_cols = [
    'Journal Name',
    'Manuscript ID',
    'Manuscript Title',
    'Author Names',
    'Submission Date',
    'Decision Date',
    'Accept or Reject Final Decision'
]
articles = df[allowed_cols].to_dict('records')

# Which might look like:
"""  
articles = [ # Fake values fabricated for the example
{
      "Journal Name": "FakeJ1",
      "Manuscript Title": "Cosmological Consequences of a Rolling Homogeneous Scalar Field",
      "Author Names": "Ratra, Bharat; Peebles, P J E",
      "Accept or Reject Final Decision": "",
      "Decision Date": "1987-06-14T13:29:58.999Z", 
      "Submission Date": "1987-06-14T13:29:58.999Z",
      "Manuscript ID": "ABC-18-070",
    }
]
"""

# @see below for configuration details.
config = {
    "threshold": 70, # Filters out matches which are less than this nubmer  
}

# The CrossRef API requires an email address for lookups.    
email = "someome@example.com"

# Define a 'results' list.
results = []

# Run match
ScholarOneRejectedArticlesMatch(
    articles=articles,
    config=config,
    email=email,
    results=results
    ).match()

print(results)
```


## Example output

Example output when match found:
```json
[
  {
  "manuscript_id": "ABC-18-070", 
  "raw_manuscript_id": "ABC-18-070", 
  "journal_name": "FakeJ1", 
  "manuscript_title": "Cosmological Consequences of a Rolling Homogeneous Scalar Field", 
  "submission_date": "1987-06-14", 
  "decision_date": "1987-06-14", 
  "authors": "Bharat+Ratra, P J E+Peebles", 
  "text_sub_date": "1987-06-14", 
  "final_decision": "", 
  "match_doi": "10.1103/physrevd.37.3406", 
  "match_type": "journal-article", 
  "match_title": "Cosmological consequences of a rolling homogeneous scalar field", 
  "match_authors": "Bharat+Ratra, P. J. E.+Peebles", 
  "match_publisher": "American Physical Society (APS)", 
  "match_journal": "Physical Review D", 
  "match_pub_date": "1988-6-15", 
  "match_earliest_date": "2002-07-27", 
  "match_similarity": 89, 
  "match_one": 1, 
  "match_all": 1, 
  "match_crossref_score": 84.5133, 
  "match_crossref_cites": 2881, 
  "match_rank": 1, 
  "match_total_decision_days": 5521
  }
  ]
```

Example out when NO match found:

```json
[
  {
    "manuscript_id": "ABC-18-070",
    "raw_manuscript_id": "ABC-18-070",
    "journal_name": "The International Journal of Robotics Research",
    "manuscript_title": "Learning hand-eye coordination for robotic grasping with deep learning and large-scale data collection",
    "submission_date": "2018-10-01",
    "decision_date": "2019-01-01",
    "authors": "Levine, Sergey; Pastor, Peter; Krizhevsky, Alex; Ibarz, Julian; Quillen, Deirdre",
    "text_sub_date": "2018-07-20",
    "final_decision": "",
    "match_doi": "No Match",
    "match_type": "No Match",
    "match_title": "No Match",
    "match_authors": "No Match",
    "match_publisher": "No Match",
    "match_journal": "No Match",
    "match_pub_date": "No Match",
    "match_earliest_date": "No Match",
    "match_similarity": "No Match",
    "match_one": "No Match",
    "match_all": "No Match",
    "match_crossref_score": "No Match",
    "match_crossref_cites": "No Match",
    "match_rank": "No Match",
    "match_total_decision_days": "No Match"
  }
]
``` 

**To rebuild the training dataset and train a new model**

Note that rebuilding the training dataset relies on external APIs and can be a very slow process (a few days depending on response times). However, once acquired, model training and testing takes seconds.

```python
from rejected_article_tracker.src.ML.Train import LogReg

LogReg().best_model_to_file()
```


---
## Configuration
Configuration is set using a dictionary. The following values can be set: 

| Name | Description | Example
| --- | --- | --- |  
| `threshold` | An integer value which determines the minimum "cut off" for scoring matching articles. Any matching articles below this score will not be considered. | `70` |   
---


## License
[MIT](LICENSE.md)

## Contributing
All contributions, bug reports, bug fixes, documentation improvements, enhancements, and ideas are welcome.

Please see the [contributing](CONTRIBUTING.md) for instructions 
## Code of Conduct

### Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

### Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualised language or imagery and unwelcome sexual attention or
advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

### Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

### Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

### Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at Adam,Day@sagepub.co.uk or Andy.Hails@sagepub.co.uk. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

### Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/# Contributing
When contributing to this repository, please first discuss the change you wish to make via issue, email, or any other method with the owners of this repository before making a change.

Please note we have a [Code of Conduct](CODE_OF_CONDUCT.md), please follow it in all your interactions with the project.

## Pull Request Process

1. Ensure any install or build dependencies are removed before the end of the layer when doing a 
   build.
2. Update the README.md with details of changes to the interface, this includes new environment 
   variables, exposed ports, useful file locations and container parameters.
3. Increase the version numbers in any examples files and the README.md to the new version that this
   Pull Request would represent. The versioning scheme we use is [SemVer](http://semver.org/).
4. You may merge the Pull Request in once you have the sign-off of two other developers, or if you 
   do not have permission to do that, you may request the second reviewer to merge it for you.


## Notes for deployment 

Prerequisite - tools to be installed 

```
pip install --upgrade pip setuptools wheel
pip install tqdm
pip install --user --upgrade twine
```

**To compile for local development:**
```
scripts/compile.sh
```

**To push to pypi:**
The GitHub CI takes care of compilation and deployment to pypi.
 
You just need to the following to deploy:
1. Update version in [`setup.py`](setup.py) 
2. Commit your code.
2. Run the release script from the root of this project:

```
$ scripts/release.sh
```---
title: 'The SAGE Rejected Article Tracker'
tags:
  - Python
  - peer-review
  - academic publishing
  - crossref
authors:
  - name: Andrew Hails
    orcid: 0000-0002-1014-621X
    affiliation: 1
  - name: Adam Day
    orcid: 0000-0002-8529-9990
    affiliation: 1
affiliations:
  - name: "SAGE Publishing, 1 Oliver's Yard, London, EC1Y 1SP"
  - index: 1
date: 4 November 2020
bibliography: paper.bib

---
# Summary
Over 3m peer-reviewed research papers are published in academic journals each year [@STM]. An unknown number of research papers are also rejected by peer-review. There is little understanding of what happens to those rejected papers.

[CrossRef](https://www.crossref.org/about/) is a not-for-profit organisation which maintains a large set of metadata describing the majority of published peer-reviewed research papers. The SAGE rejected article tracker extracts knowledge from that dataset by analysing data from the [CrossRef REST API](https://github.com/CrossRef/rest-api-doc). 

Given metadata for a rejected article, the rejected article tracker will:

* search CrossRef's API to retrieve a list of possible matches and 

* select the most likely correct result from that list using simple machine learning.

The target audience for the tracker is researchers studying rejected articles. The task performed by the tracker is record-linkage, i.e., finding the correct CrossRef metadata record given incomplete data about a paper. So, while the intended use of the tracker is to track rejected articles, it can also be used by researchers performing record-linkage for other reasons, such as connecting preprints to their published versions, e.g., work by @Cabanac2021. This is a particularly topical application at the current time due to the rapid growth of preprint servers in recent years [@Hoy2020].

The tracker is available as [a Python package](https://github.com/ad48/rejected_article_tracker_pkg) with [a temporary live demonstration](https://rejectedarticlestorage.z6.web.core.windows.net/) scheduled to run until mid-2021.

# Statement of need

As the rate of creation of research manuscripts continues to grow at a rapid pace, the need to understand the peer-review process, improve efficiencies and tackle abuse becomes all the more pressing. 

Rejected article tracking has been performed in a number of research settings [@Wijnhoven2010; @Docherty2017; @Citerio2018; @Chung2020]. Typically, this is done by manually searching for rejected articles over a small dataset. However, commercial rejected article trackers are available [@HighWire; @Dimensions2021]. To date, the lack of open source tools has prevented easy acquisition of data on rejected articles for analysis.

Data acquired by rejected article tracking makes various insights into the peer-review and publication processes possible. For example: 

- It is possible to measure the rate at which rejected articles are published and cited. This provides evidence for the effectiveness of journal peer-review in identifying (or failing to identify) flaws in research.

Rejected article tracking is also valuable to the study of scientific misconduct (examples: [@Hesselmann2017; @Ding2019; @Bozzo2017]). 

Common forms of author-misconduct can be identified and studied. 

- Dual submission (where an author has submitted the same article to multiple journals simultaneously) can be detected retrospectively with a high-degree of confidence. 

- In a similar way, self-plagiarism can potentially be detected quickly and cheaply by checking new-submissions against CrossRef with the tracker. However, the well-established [CrossCheck](https://www.crossref.org/services/similarity-check/) service based on [iThenticate](https://www.ithenticate.com/) should yield superior results. 

- When a rejected article has been later published _and then retracted_ due to fraud or other misconduct, this can allow the publisher who rejected the paper to identify that case of misconduct in their own part of the peer-review system.

Finally, the rejected article tracker can also be used to link preprints with their published versions. Due to the rapid recent growth in preprint servers [@Fraser2020], there is a growing need to improve the data-quality surrounding preprints. 

The rejected article tracker is set up, by default, to accept data in the format exported by [ScholarOne](https://clarivate.com/webofsciencegroup/solutions/scholarone/), a popular system for managing peer-review. However, the input data required is minimal, so data from any peer-review management system should be easily adapted for the tracker. Instructions are given in the `Readme.md` file of [the GitHub repository](https://github.com/ad48/rejected_article_tracker_pkg).

## How the matching algorithm works

The CrossRef API is often used to perform record-linkage. A typical use-case is adding metadata to incomplete references in the reference-list of a research paper [@Tkaczyk2018]. Under this typical use-case, there are often other data available, such as journal name and publication date as well as issue, volume, or page numbers. However, if we wish to track rejected articles, it is likely that we only have the title and author names for an article and there is a lower chance that it exists in CrossRef's data (since not all rejected articles are published). So, searching the API for just these 2 things often results in incorrect results being retrieved.

We begin with a dataset of ArXiv preprint metadata retrieved from the [ArXiv OAI-PMH API](https://arxiv.org/help/oa/index). This dataset resembles journal submission data in that it includes the titles and author names of preprints. In many cases, this data also includes the DOI of the same article when it was published. The published version of the article is known as the "version of record" (VOR). This means that we know the correct result of a record-linkage process for this article. We find that the title and author lists are not always identical. Titles often undergo minor (and occasionally major) changes in the time between appearing as a preprint and publication. Author lists can also change in a number of ways (full names might be used instead of initials, or perhaps new authors are added to an author list at some point in the process).


We search the CrossRef API for each preprint's DOI as well as the best incorrect search result. This means that we can fill out 2 rows of a table of data for each preprint.

| ArXiv <br/>title | ArXiv <br/>authors | VOR <br/>title | VOR <br/>authors | correct/ <br/>incorrect |
|-|-|-|-|-|
| title1| author_list1 | title2 | author_list2 | correct |
| title1| author_list1 | title3 | author_list3 | incorrect |

We then:

- Calculate the Levenshtein distance between the titles in each row. This is normalised to a number between 0 and 100 using the `fuzz.ratio` method from the [Python `fuzzywuzzy` package](https://pypi.org/project/fuzzywuzzy/).

- Normalise all author names to a single string of `first_initial+last_name` in lower case. Then calculate 2 boolean values: one showing if there is a 100% match in author lists and one showing if there is at least 1 author name matching in the 2 lists. 

This gives us a table of numerical data:

| levenshtein_similarity | authors_match_one | authors_match_all | label |
|-|-|-|-|
| 98 | 1 | 1 | 1 |
| 70 | 1 | 0 | 0 |

Then this data can be used to create a Logistic Regression classifier model using [Scikit-Learn](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html) with `label` as our target variable. The model essentially learns what the typical difference is between a preprint title and its published version - a bespoke form of fuzzy matching. 

This model can then be used to classify search results from the CrossRef API in order to find the correct DOI, and other metadata, for a rejected article.

The complete code required to build and customise the training dataset is included in [the SAGE Rejected Article Tracker](https://github.com/ad48/rejected_article_tracker_pkg).

## The dataset

 The training dataset is also useful for other tasks such as identifying duplicate submissions. For example, if an author submits a paper to two or more journals at once, fuzzy matching on titles and author lists is an effective way to identify this behaviour. 

 A dataset similar to the one used to train the SAGE Rejected Article Tracker is available to download from [Zenodo](http://doi.org/10.5281/zenodo.5122848) [@SAGERATData].


# Acknowledgements

We thank Helen King and Martha Sedgwick for support and advice in the development of this application. We also thank the community at [_Journal of Open Source Software_ (JOSS)](https://joss.theoj.org/) for knowledgeable, patient and highly constructive feedback, particularly: [Martin Fenner](https://github.com/mfenner) (Reviewer), [Daniel Himmelstein](https://github.com/dhimmel) (Reviewer) and [Daniel S. Katz](https://github.com/danielskatz) (Editor).

# References
