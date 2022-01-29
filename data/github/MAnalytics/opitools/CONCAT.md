---
title: "R-Opitools – An Opinion Analytical Tool for Big Digital Text Document (DTD)"
date: "30th July 2021"
bibliography: paper.bib
output: pdf_document
affiliations:
- name: Crime and Well-being Big Data Centre, Manchester Metropolitan University, United Kingdom
  index: 1
tags:
- digital text document
- sentiment analysis
- opinion mining
- randomization testing
authors:
- name: Monsuru Adepeju
  orcid: 0000-0002-9006-4934
  affiliation: 1
---

# Statement of Need

Since the year 2000, various computational intelligence techniques have been developed for analyzing sentiments of users in the field of natural language processing (NLP). To date, the majority of the techniques as deployed across various fields, including social sciences [@Somasundaran:2010; @Nikolovska:2020; @Ansari:2020] and market research [@Feldman:2011; @Otaibi:2018], have focused largely on detecting subjectivity, and/or extracting and classifying sentiments and opinions in a text document. Building on this existing work, the current paper advances an opinion impact analytical tool, named `Opitools`, that not only extracts inherent themes from within a digital text document (DTD), but also evaluates the extent to which a specified theme may have contributed to the overall opinions expressed by the document. Based on this advancement, `Opitools` has wider applications in the aforementioned application fields. For example, in law enforcement, the package can be deployed to understand factors (themes) that drive public perception of police services [@Adepeju:2021]; and in product marketing, to identify factors that underlie customers satisfaction in a product.

# Implementation

Having extracted a set of thematic keywords from a digital text document, the goal is to computationally classify the sentiments expressed in each text record into positive, negative or a neutral sentiment, using a lexicon-based classification approach [@Nielsen:2011; @Adepeju:2021]. The resulting sentiment scores are combined in order to estimate the overall opinion score of the document. To assess the impacts of a selected theme (or a subject) on the estimated opinion score, we simply ask the question; *If expected opinion scores were generated under the null hypothesis, how likely would we be to find a score higher than the estimated score?*. The question is answered by employing a non-parametric randomization testing strategy [@Fisher:1935; @Good:2006] which involves random re-assignment of sentiment labels of the original text document to derive the expectation distribution, which is then compared with the observed score to obtain the statistical significance of the impacts.


# Key Functionalities

The package contains text exploratory functions for extracting themes from a digital text document. In order to conduct impact analysis, a user can draw on a number of interrelated functions to compute the required measures, such as the observed opinion score, the expectation distribution, and the statistical significance of impacts. Whilst different types of opinion score functions are embedded in the package, there is also a provision that allows a user to integrate his/her own pre-defined user score function. This feature is to further facilitate the uptake of the package in more application fields.

# Acknowledgment

We gratefully acknowledge the Economic and Social Research Council (ESRC), who funded the Understanding Inequalities project (Grant Reference ES/P009301/1) through which this research was conducted.

# References

# "Opitools"

An R-package for analyzing Opinions in Big Digital Text Document (DTD)

## Description

The `opitools` is an R package for exploring a digital text document (DTD) as well as performing the impact analysis of the opinions expressed by the DTD. The package is particularly suited for opinion-based DTD, such as individual-level posts (e.g. commentaries on Twitter or Facebook) and comments (as in an online product reviews). First, an `exploratory` function `word_imp` can be used to identify words relating to different themes (or subjects) that are referenced in a DTD. The `opi_impact` function can then be utilized to investigate whether an identified/specified theme (or subject) impacts the overall opinion expressed by the document. The potentials of `opitools` for application across a wide range of domains is described briefly here (see the `vignette` for details). 

[Click here](https://manalytics.github.io/opitools/index.html) to visit the built website for the package.

## Installation from `CRAN`

From an R console, type:

```{r}
#install.packages("opitools")
library(opitools)

```

To install the developmental version of the package, type:
`remotes::install_github("MAnalytics/opitools")`. (Note: `remotes` is an extra package that needed to be installed prior to the running of this code). 

Please, report any installation problems in the issues.

## Example usage

Below is an example usage of how `opitools` can be employed to identify themes (or subjects) from a DTD and then deployed to perform opinion impact analysis. 


### Importing the dataset

The `policing_dtd` - a DTD comprising Twitter posts about police/policing in a neighbourhood, will be used in this demonstration.

```r
> data(policing_dtd)

```

### Identify themes (subjects)

Utilize `word_imp` function to highlights terms or words in accordance to their importance in the DTD. Through visual inspection, collate related terms or words that denote specific theme or subject from the generated wordcloud. Several words relating to the COVId-19 pandemic, including 'infect', 'pandemic', and 'lockdown' can be identified as been important in the document. These words are provided as `covid_theme` in the package. The function can be ran as follows (see the documentation for details): 

```r
> p1a <- word_imp(textdoc = policing_dtd, metric= "tf", 
                           words_to_filter=c("police","policing"))
                           
#Note: `policing_dtd` is a dataframe
```

### Impact analysis

The impact analysis can be conducted as follows:

```r
#Running the analysis

results <- opi_impact(textdoc = policing_dtd, theme_keys=covid_theme, metric = 1,
                       fun = NULL, nsim = 99, alternative="two.sided", quiet=FALSE)
                       
#Note: `policing_dtd` is a dataframe    

print(results)

$test
[1] "Test of significance (Randomization testing)"

$criterion
[1] "two.sided"

$exp_summary
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 -8.240  -5.880  -5.880  -5.028  -3.530  -1.180 

$p_table


|observed_score |S_beat |nsim |pvalue  |signif |
|:--------------|:------|:----|:-------|:------|
|-5.88          |56     |99   |0.52    |'      |

$p_key
[1] "0.99'"   "0.05*"   "0.025**" "0.01***"

$p_formula
[1] "(S_beat + 1)/(nsim + 1)"

#......

```

The research question of the analysis above can be stated as follows:

***RQ1***: "Does COVID-19 pandemic influence public opinion on neighourhood policing?"

The output shows an overall negative opinion (`-5.88`) of the public on the neighbourhood policing, and that the pandemic has not had a significant impacts (`pvalue` = 0.52) on the opinion expressed by the public. (More detailed explanation can be found in the study [adepeju, M. and Jimoh, F. (2021)](https://www.scirp.org/journal/paperinformation.aspx?paperid=107836)). 

### Other applications

Table 1 summarizes the analysis using different example datasets provided in the package. The output from the law enforcement application (as in above) is entered in the first row of the table. Other research questions investigated are as follow:

***RQ2***: "How does the democratic candidate (Hillary Clinton) affects viewers’ opinion of the presidential debate?"

***RQ3a***: "Do the refreshment outlets/items impact customers’ opinion of the Piccadilly train services?"

***RQ4b***: "Do the signages influence customers’ opinion of the Piccadilly train services?"


***Table 1. Impact analysis results***

```r
|  RQs   | Primary data | Theme_keys        | Score function   | Observed Score (S) | P-value    |
|:-----: | :----------: | :---------------: | :---------------:| :-----------------:| :---------:| 
|  RQ1   | policing_dtd | covid_theme       | 'Polarity score' | -5.88              | 0.52       |
|  RQ2   | debate_dtd   | direct input      | 'Polarity score' | -0.33              | 0.93       |
|  RQ3a  | reviews_dtd  | refreshment_theme | 'Polarity score' | 67.92              | 0.01       |
|  RQ3b  | reviews_dtd  | signage_theme     | 'Polarity score' | 67.92              | 0.1        |

```

In each example, the same opinion score function is employed (`metric = 1`, i.e. the `polarity score = (P - N)/(P + N)*100`, where `P` and `N` represent `positive` and `negative` sentiments, respectively). See the documentation for details. Employing a threshold of `p=0.05`, any p-values less or equal to the threshold (e.g. RQ3a) represent a significant impact of the specified theme (i.e. `refreshment_theme`) on the overall opinion score computed based on the `reviews_dtd`.


### References
1. Adepeju, M. and Jimoh, F. (2021). An Analytical Framework for Measuring Inequality in the Public Opinions on Policing – Assessing the impacts of COVID-19 Pandemic using Twitter Data. [click here:](https://www.scirp.org/journal/paperinformation.aspx?paperid=107836)


## Code of Conduct

Please note that the `opitools` package is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.
---
title: "NEWS.md"
authors: "Monsuru Adepeju"
date: "29 July 2021"
output: html_document
---



##
Date: 29th April 2021
##

Updates: 

1. Added new function (`word_imp`)
2. Added new example datasets (i.e. `reviews_dtd.rda` and `debate_dtd.rda`)
3. Added more sample analyses in the vignette.

Thanks,

Monsuru
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
and learning from the experience
* Focusing on what is best not just for us as individuals, but for the overall
community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards
of acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies
when an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at 
m.adepeju@mmu.ac.uk. All complaints will be reviewed and investigated promptly 
and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series of
actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or permanent
ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior, harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within the
community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0,
available at https://www.contributor-covenant.org/version/2/0/
code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at https://
www.contributor-covenant.org/translations.
---
name: Feature request
about: Kindly suggest an idea (an argument, a function, etc.) for this project
title: "[Add feature]:"
labels: documentation
assignees: MAnalytics

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Are alternatives you've considered? Please briefly describe**
A clear and concise description of any alternative solutions or features you've considered.

**Additional information**
Add any other information or screenshots about the feature request here.
---
name: Bug report
about: Let us know what is wrong
title: "[BUG]: "
labels: bug
assignees: MAnalytics

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
If possible, state the steps to reproduce the behavior:

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Additional information**
Add any other information about the problem here.
