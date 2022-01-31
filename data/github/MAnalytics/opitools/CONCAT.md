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
---
title: "An Opinion Analytical Tool for Big Digital Text Document - A User Guide"

author: |
  | `Author:`
  | `Adepeju, M.`
  | `Big Data Centre, Manchester Metropolitan University, Manchester, M15 6BH, UK`
  
date: |
  | `Date:`
  | ``r Sys.Date()``

output:
  html_vignette
  
#output:
  #pdf_document: default
  
#dev: png
#output:
  #word_document: default
  #always_allow_html: yes
#  pdf_document: default
always_allow_html: yes
#fig_caption: yes
bibliography: references.bib

abstract: The development of `'opitools'` is instigated by the lack of tools for performing opinion impact analysis of a digital text document (DTD). The package includes a number of interrelated functions for exploring and analyzing text records in order to complete an opinion impact analysis of a digital text document (DTD). The utility of `Opitools` is demonstrated with examples from the law enforcement, electoral politics and product marketing.

vignette: >
  %\VignetteIndexEntry{An Opinion Analytical Tool for Big Digital Text Document - A User Guide}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

<style type="text/css">

h1.title {
  font-size: 26px;
  line-height: 130%;
  color: Black;
  text-align: center;
}

h2.subtitle {
  font-size: 13px;
  line-height: 120%;
  color: Black;
  text-align: center;
}

h4.author { /* Header 4 - and the author and data headers use this too  */
  font-size: 17px;
  font-family: "Arial";
  color: Black;
  text-align: center;
}
h4.date { /* Header 4 - and the author and data headers use this too  */
  font-size: 17px;
  font-family: "Arial", Times, serif;
  color: Black;
  text-align: center;
}

h4.abstract { /* Header 4 - and the author and data headers use this too  */
  font-size: 10px;
  font-family: "Arial", Times, serif;
  color: black;
  text-align: center;
}

h4.institute{ /* Header 4 - and the author and data headers use this too  */
  font-size: 10px;
  font-family: "Arial", Times, serif;
  color: black;
  text-align: center;
}

body, td {
   font-size: 14px;
}
code.r{
  font-size: 13px;
}
pre {
  font-size: 13px
}
h1 { /* Header 1 */
  font-size: 16px;
  color: DarkBlue;
}
h2 { /* Header 2 */
    font-size: 16px;
  color: DarkBlue;
}
h3 { /* Header 3 */
  font-size: 15px;
  font-family: "Times New Roman", Times, serif;
  color: DarkBlue;

</style>

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r functions, include=FALSE}
# A function for captioning and referencing images
fig <- local({
    i <- 0
    ref <- list()
    list(
        cap=function(refName, text) {
            i <<- i + 1
            ref[[refName]] <<- i
            paste("Figure ", i, ": ", text, sep="")
        },
        ref=function(refName) {
            ref[[refName]]
        })
})
```




# Introduction

The `'opitools'` provides the mechanism for carrying out opinion impact analysis of a digital text document (DTD). The package functions can be categorized into two groups, namely (a) exploratory - for exploring terms in a document, e.g. importance of words and their statistical distribution, and (b) analytical - for computing metrics on the impacts of a theme or subject on the opinion expressed in a document. The potentials of `opitools` for application across a wide variety of domains is demonstrated with examples from the `law enforcement` to examine the impacts of covid-19 on the citizens' opinion of their neighbourhood policing; from `electoral politics` to examine the extent to which a political candidate drives viewers opinion of a political debate; and from `product marketing` to examine what features in and around a train station influence customers reviews of the train services;  


### Set the Working directory

```{r, message=FALSE, eval=FALSE}
WORKING_DIR <- 'C:/R/Github/JGIS_Policing_COVID-19'

#setting working directory
setwd(WORKING_DIR)
```


### Instal libraries
```{r, include=TRUE, message=FALSE, eval=TRUE}
library(opitools) #for impact analysis
require(knitr) #for rendering the vignette
library(rvest)
library(kableExtra) #for designing tables
library(dplyr) #for data analysis
library(cowplot) #for plot design


```

# Example datasets

The following datasets will be employed in our demonstration. The datasets are automatically installed upon the installation of `Opitools`. 

```{r, echo=FALSE, include=FALSE}
col1 <- c("1", "2", "3")
col2 <- c("`policing_dtd`", "`debate_dtd`", "`reviews_dtd`")
col3 <- c("`Law Enforcement`", "`Electoral Politics`", "`Product Marketing`")
col4 <- c("A digital text document (DTD) containing twitter posts, within a geographical neighbourhood, on police/policing during the 2020 COVID-19 pandemic", "A DTD containing individual comments on the video showing the debate between two United States presidential candidates (Donald Trump and Hillary Clinton) in September of 2016. (Credit: NBC News).", "A DTD containing customers reviews of the Piccadilly train station (Manchester, UK). The records cover from July 2016 to March 2021.")
col5 <- c("www.twitter.com", "www.youtube.com", "www.tripadvisor.co.uk")
tble1 <- data.frame(col1, col2, col3, col4, col5)
tble1 <- tble1
```



```{r table1, results='asis', echo=FALSE, tidy.opts=list(width.cutoff=50)}
knitr::kable(tble1, caption = "Table 1. `Example datasets`", col.names = c("SN","Data","Application","Details", "Data Source")) %>%
  kable_styling(full_width = F) %>%
  column_spec(1, bold = T, border_right = T) %>%
  column_spec(2, width = "8em", background = "white") %>%
  column_spec(3, width = "12em", background = "white") %>%
  column_spec(4, width = "16em", background = "white")#%>%
  #row_spec(3:5, bold = T, color = "white", background = "#D7261E")
```
 

# Functions

The function in `'opitools'` can be categorized into two groups, namely (1) exploratory - for exploring terms or words in a text document, and (2) analytical - for computing metrics for the impact analysis. 

### Exploratory function

Table 2 shows two key exploratory functions embedded in `'opitools'`, namely `'word_distrib'` and `'word_imp'`. Details as follow:

```{r, echo=FALSE, include=FALSE}
col1 <- c("1", "2")
col2 <- c("`word_distrib`","`word_imp`")
col3 <- c("`Words Distribution`","`Importance of words (terms) embedded in a text document`")
col4 <- c("Examines the extent to which the terms (words) in a DTD follow the Zipf's distribution (Zipf 1934) - the ideal natural language model", "Produces a table or graphic that highlights the importance of individual terms (or words) in a DTD.")
tble2 <- data.frame(col1, col2, col3, col4)
tble2 <- tble2
```

```{r table2, results='asis', echo=FALSE, tidy.opts=list(width.cutoff=50)}
knitr::kable(tble2, caption = "Table 2. `Exploratory` functions", col.names = c("SN","Function","Title","Description")) %>%
  kable_styling(full_width = F) %>%
  column_spec(1, bold = T, border_right = T) %>%
  column_spec(2, width = "8em", background = "white") %>%
  column_spec(3, width = "12em", background = "white") %>%
  column_spec(4, width = "16em", background = "white")#%>%
  #row_spec(3:5, bold = T, color = "white", background = "#D7261E")
```


### Example data exploration

Below, the use of `'word_distrib'` and `'word_imp'` functions is demonstrated with the example datasets `'tweets'` (see documentations). The `'word_distrib'` can be used to answer a RQ such as; "How similar is a given DTD to the ideal natural language text?". In other words, does the word usage in a DTD follow the ideal natural language usage represeted by the Zipf's distribution [@Zipf1936]? Basically, the function examines the log rank-frequency of the terms across the document. The level of similarity is examined by the extent to which the frequency of each term is inversely proportional to their rank in a frequency table. Using a randomized Twitter data (embedded in the package), 

```{r, message=FALSE, include = TRUE, eval=FALSE}

# Load data
data(tweets)

# Get the texts
tweets_dat <- as.data.frame(tweets[,1])

# Run function
plt = word_distrib(textdoc = tweets_dat)

#Note: `tweets_dat` is a dataframe 

#Show Zipf's distribution:
plt$plot

```

```{r figs1, echo=FALSE, fig.width=5,fig.height=6,fig.align="center", fig.cap=fig$cap("figs1", "Data freq. plot vs. Zipf's distribution")}
knitr::include_graphics("zipf.png")
```

For a natural language text, the Zipf's distribution plot has a negative slope with each term falling on a straight line. Any variance from this (ideal) trend can be attributed to imperfections in the term usage. For example, the presence of strange terms or 'made-up' words can cause a deviation from the ideal trend. The result of our example dataaset is shown in Figure `r fig$ref("figs1")`. The graph can be divided into the three sections: the upper, the middle and the lower sections. By fitting a regression line (an ideal Zipf's distribution), we can see what the slope of the upper section is quite different from the middle and the lower sections of the graph. The variance at the high rank terms (usually, common terms, such as 'the', 'of', and 'at',) indicates an imperfection because of lack of perfect alignment with the straight line. For Twitter data, this variance suggests a significant use of abbreviated terms, such as using "&" or "nd" instead of the term "`and`". Apart from the small variance at the upper section of the graph, we can state that the law holds within most parts of the document.

The second exploratory function `'word_imp'` can be used to highlight the importance of each terms in a DTD. The level of importance can then be used to identify different themes (or subjects) from the DTD. Two measures, namely (i) `'term frequency (tf)'` and `term frequency inverse document frequency (tf-idf)` [@Silge2016] can be used to highlight the importance of terms in a DTD. In Figure 2, results 1, 2 and 3 are generated from the datasets; `policing_dtd`, `debate_dtd` `reviews_dtd`, respectively. The results for each data sample based on the `'tf'` and `'tf_idf'` measures are labeled as A and B, respectively. The function draws from `wordcloud2` package in R [@Dawei2018], in order to highlight words importance as shown in Figure 2.

```{r, message=FALSE, include = TRUE, eval=FALSE}

#Load datasets

data("policing_dtd")
data("debate_dtd")
data("reviews_dtd")


#Words importance of public tweets on neighbourhood policing 
#based on  (a) ‘tf’ and (b) 'tf-idf' metrics
p1a <- word_imp(textdoc = policing_dtd, metric= "tf", 
                           words_to_filter=c("police","policing"))

#Note: `policing_dtd` is a dataframe 

p1b <- word_imp(textdoc = policing_dtd, metric= "tf-idf", 
                           words_to_filter=c("police","policing"))

#Note: 'words_to_filter' parameter is used to eliminate non-necessary words that 
#may be too dominant in the DTD.

#Words importance  of comments on a video of a political debate 
#based on (a) ‘tf’ and (b) 'tf-idf' metrics
p2a <- word_imp(textdoc = debate_dtd, metric= "tf", 
                words_to_filter=c("trump","hillary")) 

p2b <- word_imp(textdoc = debate_dtd, metric= "tf-idf", 
                words_to_filter=c("trump","hillary")) 

#Words importance of customer reviews of a transport service 
#based on (a) ‘tf’ and (b) 'tf-idf' metrics 
p3a <- word_imp(textdoc = reviews_dtd, metric= "tf", 
                words_to_filter=c("station")) 

p3b <- word_imp(textdoc = reviews_dtd, metric= "tf-idf",  
                words_to_filter=c("station")) 

#outputs
p1a$plot; p1b$plot; p2a$plot; p2b$plot; p3a$plot; p3b$plot

```

```{r figs2, echo=FALSE, fig.width=3,fig.height=4,fig.align="center", fig.cap=fig$cap("figs2", "Highlighting words importance from a DTD")}
knitr::include_graphics("wordcloud.png")
```

In Figure `r fig$ref("figs2")`, the size of a word represents its level of importance according to a selected metric (i.e. 'tf' or 'tf-idf'). From each representation, a user can easily identify (related) words that denote certain themes or subjects within the document. In Figure 1A, words such as 'lockdown', 'infect' and 'covid', refer to the sentiments relating to the pandemic under which the police were operating as at the time tweets were posted. "Does the pandemic has any significant impacts on the public opinion of the police activities during the period?". This is an example of the RQs from social science research that can be investigated using `Opitools`. On the other hand, Figure 1B shows that the word importance may vary significantly by the choice of the measure between 'tf' and 'tf-idf'. Figure 1B shows that a significant amount of the most important words using `tf-idf` are distinct from those ones resulting from using `tf`. 

In marketing research (Figure 3A), we see that related words, such as 'restaurant', 'shops', 'food' and 'coffee', denoting `refreshment items`, are highlighted as been the most important to the customers. "`Do these refreshment items impact the overall customers' opinion of the train services?`". Similarly, related words such as 'board', 'map' and 'display', denoting `signages` around the station may have influenced customers opinion of the station.


### Analytical functions

Table 3 shows the list of analytical functions of the `'Opitools'` package.

```{r, echo=FALSE, include=FALSE}
col1 <- c("3", "4", "5")
col2 <- c("`opi_score`","`opi_sim`", "`opi_impact`")
col3 <- c("`Computes the overall opinion score of a text document`",
          "`Simulates the opinion expectation distribution of a text document`",
          "`Computes the statistical significance of impacts of a specified/selected theme (or subject) on the overall opinion score of a document`")
col4 <- c("Given a text document, this function computes the overall opinion score based on the proportion of text records classified as expressing positive, negative or a neutral sentiment about the subject.",  "This function simulates the expectation distribution of the observed opinion score (computed using the `opi_score` function).",  "This function assesses the impacts of a theme (or subject) on the overall opinion computed for a text document. The text records relating to the theme in question should be identified and provided as input to this function")
tble3 <- data.frame(col1, col2, col3, col4)
tble3<- tble3
```

```{r table3, results='asis', echo=FALSE, tidy.opts=list(width.cutoff=50)}
knitr::kable(tble3, caption = "Table 3. `Impact Analytical` function", col.names = c("SN","Function","Title","Description")) %>%
  kable_styling(full_width = F) %>%
  column_spec(1, bold = T, border_right = T) %>%
  column_spec(2, width = "8em", background = "white") %>%
  column_spec(3, width = "12em", background = "white") %>%
  column_spec(4, width = "16em", background = "white")#%>%
  #row_spec(3:5, bold = T, color = "white", background = "#D7261E")
```

The key analytical function is the `opi_impact()` which draws from `opi_score()` and `opi_sim()` to compute the observed opinion score and the expectations, respectively. The observed opinion score is compared with the expectation in order to derive the statistical significance of impacts of a specified/selected theme. Using the example dataset `policing_dtd`:  

The impact analysis can be performed as follows: 

```{r, message=FALSE, include = TRUE, eval=FALSE}

#Application: Law enforcement

#Load DTD
data(policing_dtd)

#Load theme keywords
data(covid_theme)

# Run the analysis
output1 <- opi_impact(policing_dtd, theme_keys=covid_theme, metric = 1,
                       fun = NULL, nsim = 99, alternative="two.sided",
                       quiet=TRUE)

#Note: `policing_dtd` is a dataframe 

print(output1)

```

To print results: 

```{r, echo=TRUE, message=FALSE, eval=FALSE}

> output1

$test
[1] "Test of significance (Randomization testing)"

$criterion
[1] "two.sided"

$exp_summary
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 -8.240  -5.880  -5.880  -4.837  -3.530  -1.180 

$p_table

|observed_score |S_beat |nsim |pvalue |signif |
|:--------------|:------|:----|:------|:------|
|-5.88          |51     |99   |0.52   |'      |

$p_key
[1] "0.99'"   "0.05*"   "0.025**" "0.01***"

$p_formula
[1] "(S_beat + 1)/(nsim + 1)"

$plot
```

* The descriptions of output variables are as follow:

  + `test` - title of the analysis

  + `criterion` - criterion for determining the significance value

  + `exp_summary` - summary of expected opinion scores
  
  + `p_table` - details of Statistical Significance

  + `p_key` - keys for interpreting the statistical significance value

  + `p_formula` - function of opinion score employed
  
  + `plot` - plot showing Percentage proportion of classes


The output shows an overall negative opinion (`-5.88`) of the public on the neighbourhood policing, and that the pandemic has not had a significant impacts (`pvalue` = 0.52) on the opinion expressed by the public. 

To display the graphics showing the proportion of various sentiment classes (as in Figure `r fig$ref("figs3")`), type `output$plot` in the console.

```{r figs3, echo=FALSE, fig.width=5,fig.height=6,fig.align="center", fig.cap=fig$cap("figs3", "Percentage proportion of classes")}
knitr::include_graphics("likert.png")
```

Table 4 summarizes the results of impact analysis using all the example datasets in Table 1 with their respective research questions (RQs). The outputs from the law enforcement application (as in above) is entered as the first record of the table. In each example, we employed the same score function (`metric = 1`, i.e. `polarity score = (P - N)/(P + N)*100`, where `P` and `N` represent `positive` and `negative` sentiments, respectively). See the documentation for details.

```{r, echo=FALSE, include=FALSE}
col1 <- c("1", "2", "3a","3b")
col2 <- c("Does COVID-19 pandemic influence public opinion on neighourhood policing?", "How does the democratic candidate (Hillary Clinton) affects viewers' opinion of the presidential debate?", "Do the refreshment outlets/items impact customers' opinion of the Piccadilly train services?", "Do the signages influence customers' opinion of the Piccadilly train services?")
col3 <- c("`policing_dtd`","`debate_dtd`","`reviews_dtd`", "`reviews_dtd`")
col4 <- c("`covid_theme`","direct input","`refreshment_theme`", "`signage_theme`")
col5 <- c("two.sided", "two.sided", "two.sided", "two.sided")
col6 <- c("Polarity score", "Polarity score", "Polarity score", "Polarity score")
col7 <- c("-5.88", "-0.33", "67.92", "67.92")
col8 <- c("0.52", "0.93", "0.01", "0.1")
tble4 <- data.frame(col1, col2, col3, col4, col5, col6, col7, col8)
tble4 <- tble4
```


```{r table4, results='asis', echo=FALSE, tidy.opts=list(width.cutoff=50)}
knitr::kable(tble4, caption = "Table 4. `Impact analysis`", col.names = c("SN.","RQs","Primary data","theme_keys","criterion", "Score function", "Observed score (S)", "p-value")) %>%
  kable_styling(full_width = F) %>%
  column_spec(1, bold = T, border_right = T) %>%
  column_spec(2, width = "26em", background = "white") %>%
  column_spec(3, width = "8em", background = "white") %>%
  column_spec(4, width = "8em", background = "white") %>%
  column_spec(5, width = "8em", background = "white") %>%
  column_spec(6, width = "8em", background = "white") %>%
  column_spec(7, width = "8em", background = "white")%>%
  column_spec(8, width = "8em", background = "white")#%>%
  #row_spec(3:5, bold = T, color = "white", background = "#D7261E")
```
 
In the application to electoral politics, we attempt to identify factors that may have contributed to viewers opinion of a presidential debate. The statistics (`S=-0.33, pvalue = 0.93`) show that viewers opinion of the debate is negative. However, viewers opinion of the democratic candidate (Hillary Clinton) has no impacts of their overall negative opinion of the debate.

Lastly in the marketing application (3a & 3b), we tested the impacts of two themes, namely; `refreshment items/outlets` and the `signages`, on the customers' opinions of the Piccadilly train services. The statistics (`S=67.92, pvalue = 0.01`) show an overall positive customers' opinions of the train station or services, with the refreshments items/outlets (in and around the station) contributing significantly towards the positive opinions. On the other hand, the the signages (in and around the station) have no significant influence (p=0.1) on the positive review of the train service. 


### Using a user-defined opinion score function

An opinion score function may be defined in a variety of ways depending on the application domain. For instance, [@Razorfish2019] defines customers' opinion score of a product brand as `score = (P + O - N)/(P + O + N)`, where `P`, `O`, and `N`, represent the amount/proportion of positive, neutral and negative, sentiments, respectively. In order for a user to be able to plug-in their own opinion score function, we provided the parameter `fun` in the `opi_impact` function. This allow users to integrate new functions and test new theories. Employing the example function, we can re-run the analysis as follows: 

(a) Define the user-defined function: 

```{r, echo=TRUE, message=FALSE, eval=FALSE}

#The equation
Score = (P + O - N)/(P + O + N)

#Corresponding function
myfun <- function(P, N, O){
   score <- (P + O - N)/(P + O + N)
   return(score)
}

```

(b) Integrate the function and run analysis (for example)


```{r, echo=TRUE, message=FALSE, eval=FALSE}

output <- opi_impact(debate_dtd, theme_keys=keys, metric = 5,
                       fun = myfun, nsim = 99, alternative="two.sided",
                       quiet=TRUE)
```


# Conclusion

The `opitools` package has been developed in order to address the lack of mechanism for carrying out impact analysis of opinions for digital text document (DTD). The utility of the package was originally demonstrated in [@Adepeju2021] in its application in law enforcement domain. In addition, this vignette demonstrates the potentials of the package for application across a wide variety of areas, with examples from electoral politics and product marketing. This package is being updated on a regular basis to add more functionalities. 

We encourage users to report any bugs encountered while using the package so that they can be fixed immediately. Welcome contributions to this package which will be acknowledged accordingly. 

# References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opi_impact.R
\name{opi_impact}
\alias{opi_impact}
\title{Statistical assessment of
impacts of a specified theme from a DTD.}
\usage{
opi_impact(textdoc, theme_keys=NULL, metric = 1,
fun = NULL, nsim = 99, alternative="two.sided",
quiet=TRUE)
}
\arguments{
\item{textdoc}{An \code{n} x \code{1} list (dataframe) of
individual text records, where \code{n} is the total
number of individual records.}

\item{theme_keys}{(a list) A one-column dataframe (of any
number of length) containing a list of keywords relating
to the theme or secondary subject to be investigated.
The keywords can also be defined as a vector of characters.}

\item{metric}{(an integer) Specify the metric to utilize
for the calculation of opinion score. Default: \code{1}.
See detailed documentation
in the \code{opi_score} function.}

\item{fun}{A user-defined function given that parameter
\code{metric} (above) is set equal to \code{5}.
See detailed documentation
in the \code{opi_score} function.}

\item{nsim}{(an integer) Number of replicas (ESD) to generate.
See detailed documentation in the \code{opi_sim} function.
Default: \code{99}.}

\item{alternative}{(a character) Default: \code{"two.sided"},
indicating a two-tailed test. A user can override
this default value by specifying \code{“less”} or \code{“greater”}
to run the analysis as one-tailed test when the observed score
is located at the lower or upper regions of the expectation
distribution, respectively. Note: for \code{metric=1},
the \code{alternative} parameter should be
set equal to \code{"two.sided"} because the opinion score is
bounded by both positive and negative values. For an opinion
score bounded by positive values, such as when
\code{metric = 2, 3 or 4}, the \code{alternative} parameter
should be set as "greater", and set as "less" otherwise.
If metric parameter is set equal to \code{5}, with a user-defined
opinion score function (i.e. \code{fun} not NULL ), the user is required
to determine the limits of the opinion scores, and set the
\code{alternative} argument appropriately.}

\item{quiet}{(TRUE or FALSE) To suppress processing
messages. Default: \code{TRUE}.}
}
\value{
Details of statistical significance of impacts
of a secondary subject \code{B} on the opinion concerning the
primary subject \code{A}.
}
\description{
This function assesses the impacts of a theme
(or subject) on the overall opinion computed for a DTD
Different themes in a DTD can be identified by the keywords
used in the DTD. These keywords (or words) can be extracted by
any analytical means available to the users, e.g.
\code{word_imp} function. The keywords must be collated and
supplied this function through the \code{theme_keys} argument
(see below).
}
\details{
This function calculates the statistical
significance value (\code{p-value}) of an opinion score
by comparing the observed score (from the \code{opi_score}
function) with the expected scores (distribution) (from the
\code{opi_sim} function). The formula is given as
\code{p = (S.beat+1)/(S.total+1)}, where \code{S_total} is the total
number of replicas (\code{nsim}) specified, \code{S.beat} is number of replicas
in which their expected scores are than the observed score (See
further details in Adepeju and Jimoh, 2021).
}
\examples{

# Application in marketing research:

# data -> 'reviews_dtd'
# theme_keys -> 'refreshment_theme'

#RQ2a: "Do the refreshment outlets impact customers'
#opinion of the services at the Piccadilly train station?"

##execute function
output <- opi_impact(textdoc = reviews_dtd,
          theme_keys=refreshment_theme, metric = 1,
          fun = NULL, nsim = 99, alternative="two.sided",
          quiet=TRUE)

#To print results
print(output)

#extracting the pvalue in order to answer RQ2a
output$pvalue


}
\references{
(1) Adepeju, M. and Jimoh, F. (2021). An Analytical
Framework for Measuring Inequality in the Public Opinions on
Policing – Assessing the impacts of COVID-19 Pandemic using
Twitter Data. https://doi.org/10.31235/osf.io/c32qh
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{osd_data}
\alias{osd_data}
\title{Observed sentiment document (OSD).}
\format{
A dataframe with the following variables:
\itemize{
\item ID: numeric id of text record with valid
resultant sentiments score and classification.
\item sentiment: Containing the sentiment classes.
\item keywords: Indicator to show whether a secondary
keyword in present or absent in a text record.
}
}
\usage{
osd_data
}
\description{
A tidy-format list (dataframe) showing the resulting
classification of each text record into positive, negative
or neutral sentiment. The second column of the dataframe consists of
labels variables \code{present} and \code{absent} to indicate whether any of
the secondary keywords exist in a text record.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{tweets}
\alias{tweets}
\title{Fake Twitter posts on police/policing 2}
\format{
A dataframe with the following variables:
\itemize{
\item text: individual text records
\item group: real/arbitrary groups of text records
}
}
\usage{
tweets
}
\description{
A text document (an DTD) containing twitter posts
(for an anonymous geographical location 2) on police/policing
(primary subject A). The DTD includes
posts that express sentiments on policing in relation to
the COVID-19 pandemic (Secondary subject B)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/word_distrib.R
\name{word_distrib}
\alias{word_distrib}
\title{Words Distribution}
\usage{
word_distrib(textdoc)
}
\arguments{
\item{textdoc}{\code{n} x \code{1} list (dataframe)
of individual text records, where \code{n} is the number
of individual records.}
}
\value{
A list of word ranks and their respective
frequencies, and a plot showing the relationship between
the two variables.
}
\description{
This function examines whether the distribution
of word frequencies in a text document follows the Zipf distribution
(Zipf 1934). The Zipf's distribution is considered the ideal
distribution of a perfect natural language text.
}
\details{
The Zipf's distribution is most easily observed by
plotting the data on a log-log graph, with the axes being
log(word rank order) and log(word frequency). For a perfect
natural language text, the relationship between the word rank
and the word frequency should have a negative slope with all points
falling on a straight line. Any deviation from the straight
line can be considered an imperfection attributable to the
texts within the document.
}
\examples{

#Get an \code{n} x 1 text document
tweets_dat <- data.frame(text=tweets[,1])
plt = word_distrib(textdoc = tweets_dat)

plt

}
\references{
Zipf G (1936). The Psychobiology of Language.
London: Routledge; 1936.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{signage_theme}
\alias{signage_theme}
\title{Keywords relating to signages at train stations}
\format{
A dataframe containing one variable:
\itemize{
\item keys: list of keywords
}
}
\usage{
signage_theme
}
\description{
List of signages at the Piccadilly Train
Station (Manchester)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{refreshment_theme}
\alias{refreshment_theme}
\title{Keywords relating to facilities at train stations}
\format{
A dataframe containing one variable:
\itemize{
\item keys: list of keywords
}
}
\usage{
refreshment_theme
}
\description{
List of words relating to refreshments that can
be found at the Piccadilly Train
Station (Manchester)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opi_sim.R
\name{opi_sim}
\alias{opi_sim}
\title{Simulates the opinion expectation distribution
of a digital text document.}
\usage{
opi_sim(osd_data, nsim=99, metric = 1, fun = NULL, quiet=TRUE)
}
\arguments{
\item{osd_data}{A list (dataframe). An \code{n} x \code{3}
OSD, in which \code{n} represents the length of the
text records that have been successfully classified as
expressing positive, negative or a neutral sentiment.
Column \code{1} of the OSD is the text record ID,
column \code{2} shows the sentiment classes (i.e. positive,
negative, or neutral), while column \code{3} contains two
variables: \code{present} and \code{absent} indicating records that
include and records that do not include any of the specified
theme keywords, respectively.}

\item{nsim}{(an integer) Number of replicas (ESD) to simulate.
Recommended values are: 99, 999, 9999, and so on. Since the run time
is proportional to the number of replicas, a moderate number of
simulation, such as 999, is recommended. Default: \code{99}.}

\item{metric}{(an integer) Specify the metric to utilize for the
calculation of the opinion score. Default: \code{1}. See
details in the documentation of \code{opi_score} function.
The input argument here must correspond to that of \code{opi_score}
function in order to compute a statistical significance value (p-value).}

\item{fun}{A user-defined function given that parameter
\code{metric} is set equal to \code{5}. See details in the
documentation of the \code{opi_score} function.}

\item{quiet}{(TRUE or FALSE) To suppress processing
messages. Default: \code{TRUE}.}
}
\value{
Returns a list of expected opinion scores with length equal
to the number of simulation (\code{nsim}) specified.
}
\description{
This function simulates the expectation distribution of the
observed opinion score (computed using the \code{opi_score} function).
The resulting tidy-format dataframe can be described as the
\verb{expected sentiment document (ESD)} (Adepeju and Jimoh, 2021).
}
\details{
Employs non-parametric randomization testing approach in
order to generate the expectation distribution of the observed
opinion scores (see details in Adepeju and Jimoh 2021).
}
\examples{

#Prepare an osd data from the output
#of opi_score function.

score <- opi_score(textdoc = policing_dtd,
                     metric = 1, fun = NULL)
#extract OSD
OSD <- score$OSD
#note that OSD is shorter in length
#than policing_dtd, meaning that some
#text records were not classified

#Bind a fictitious indicator column
osd_data2 <- data.frame(cbind(OSD,
           keywords = sample(c("present","absent"), nrow(OSD),
           replace=TRUE, c(0.35, 0.65))))

#generate expected distribution
exp_score <- opi_sim(osd_data2, nsim=99, metric = 1,
                                 fun = NULL, quiet=TRUE)
#preview the distribution
hist(exp_score)

}
\references{
(1) Adepeju, M. and Jimoh, F. (2021). An Analytical
Framework for Measuring Inequality in the Public Opinions on
Policing – Assessing the impacts of COVID-19 Pandemic using
Twitter Data. https://doi.org/10.31235/osf.io/c32qh
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{reviews_dtd}
\alias{reviews_dtd}
\title{Customer reviews from tripadvisor website}
\format{
A dataframe containing one variable
\itemize{
\item text: individual text records
}
}
\usage{
reviews_dtd
}
\description{
A text document (an DTD) containing the  customer reviews
of the Piccadilly train station (Manchester) downloaded
from the www.tripadvisor.co.uk'. The reviews cover from
July 2016 to March 2021.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{covid_theme}
\alias{covid_theme}
\title{keywords relating to COVID-19 pandemics}
\format{
A dataframe containing one variable:
\itemize{
\item keys: list of keywords
}
}
\usage{
covid_theme
}
\description{
A list of keywords relating to the COVID-19 pandemic
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{debate_dtd}
\alias{debate_dtd}
\title{Comments on a video of a political debate.}
\format{
A dataframe containing one variable
\itemize{
\item text: individual text records
}
}
\usage{
debate_dtd
}
\description{
A DTD containing individual comments on a video
showing the first debate between two US presidential
nominees (Donald Trump and Hillary Clinton)
in Sept. 2016. (Credit: NBC News).
}
\details{
The DTD only include the comments within the first 24hrs
in which the video was posted. All individual comments
in which the names of both candidates are mentioned
are filtered out.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/word_imp.R
\name{word_imp}
\alias{word_imp}
\title{Importance of words (terms) embedded
in a text document}
\usage{
word_imp(textdoc, metric= "tf",
words_to_filter=NULL)
}
\arguments{
\item{textdoc}{An \code{n} x \code{1} list (dataframe) of
individual text records, where \code{n} is the total
number of individual records. An \code{n} x code{2} dataframe can
also be supplied, in which the second column represents the
labels of the pre-defined groupings of the text records,
e.g. labels of geographical areas where each text record
originates.
For an \code{n} x \code{1} dataframe, an arbitrary grouping is
automatically imposed.}

\item{metric}{(character) The measure for determining the level of
importance of each word within the text document. Options
include \code{'tf'}
representing \verb{term frequency} and \code{'tf-idf'}
representing \verb{term frequency inverse document frequency}
(Silge & Robinson, 2016).}

\item{words_to_filter}{A pre-defined vector of words (terms) to
filter out from the DTD prior to highlighting words importance.
default: \code{NULL}. This parameter helps to eliminate
non-necessary words that may be too dominant in the results.}
}
\value{
Graphical representation of words importance
according to a specified metric. A wordcloud is used
to represent words importance if \code{tf} is specified, while
facet wrapped histogram is used if \code{tf-idf} is specified.
A wordcloud is represents each word with a size corresponding
to its level of importance. In the facet wrapped histograms
words are ranked in each group (histogram) in their order
of importance.
}
\description{
Produces a wordcloud which represents the
level of importance of each word (across different text groups)
within a text document, according to a specified measure.
}
\details{
The function determines the most important words
across various grouping of a text document. The measure
options include the \code{tf} and \code{tf-idf}. The idea of \code{tf}
is to rank words in the order of their number of occurrences
across the text document, whereas \code{tf-idf} finds words that
are not used very much, but appear across
many groups in the document.
}
\examples{
#words to filter out
wf <- c("police","policing")
output <- word_imp(textdoc = policing_dtd, metric= "tf",
words_to_filter= wf)
}
\references{
Silge, J. and Robinson, D. (2016) tidytext:
Text mining and analysis using tidy data principles in R.
Journal of Open Source Software, 1, 37.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opi_score.R
\name{opi_score}
\alias{opi_score}
\title{Opinion score of a digital text document (DTD)}
\usage{
opi_score(textdoc, metric = 1, fun = NULL)
}
\arguments{
\item{textdoc}{An \code{n} x \code{1} list (dataframe) of
individual text records, where \code{n} is the total
number of individual records.}

\item{metric}{(an integer) Specify the metric to utilize for
the calculation of opinion score. Valid values include
\code{1, 2, ...,5}.
Assuming \code{P}, \code{N} and \code{O} represent positive,
negative, and neutral record sentiments, respectively,
the followings are the details of the opinion score function
represented by the numerical arguments above:
\code{1}: Polarity (percentage difference)
\code{((P - N)/(P + N))*100}, (Bound: -100\%, +100\%);
\code{2}: Polarity (proportional difference)
\code{((abs(P - N) / (P + N + O))*100},
(Bound: 0, +100\%);
\code{3}: Positivity \code{(P/ (P + N + O))*100},
(Bound: 0, +100\%); \code{4}: Negativity \code{(N / (P + N + O))*100},
(Bound: 0, +100\%) (Malshe, A. 2019;
Lowe et al. 2011). \code{5}: To pass a
user-defined opinion score function (also see the \code{fun}
parameter below.}

\item{fun}{A user-defined function given that \code{metric}
parameter (above) is set equal to \code{5}.
For example, given a defined opinion score function
\code{myfun} <- \verb{function(P, N, O)\{}
\code{("some tasks to do")}; \verb{return("a value")\}}, the input
argument of \code{fun} parameter then becomes \code{fun = myfun}.
Default: \code{NULL}.}
}
\value{
Returns an \code{opi_object} containing details of the
opinion measures from the text document.
}
\description{
Given a DTD,
this function computes the overall opinion score based on the
proportion of text records classified as expressing positive,
negative or a neutral sentiment.
The function first transforms
the text document into a tidy-format dataframe, described as the
\verb{observed sentiment document (OSD)} (Adepeju and Jimoh, 2021),
in which each text record is assigned a sentiment class based
on the summation of all sentiment scores expressed by the words in
the text record.
}
\details{
An opinion score is derived from all the sentiments
(i.e. positive, negative (and neutral) expressed within a
text document. We deploy a lexicon-based approach
(Taboada et al. 2011) using the \code{AFINN} lexicon
(Nielsen, 2011).
}
\examples{
# Use police/pandemic posts on Twitter
# Experiment with a standard metric (e.g. metric 1)
score <- opi_score(textdoc = policing_dtd, metric = 1, fun = NULL)
#print result
print(score)

#Example using a user-defined opinion score -
#a demonstration with a component of SIM opinion
#Score function (by Razorfish, 2009). The opinion
#function can be expressed as:

myfun <- function(P, N, O){
  score <- (P + O - N)/(P + O + N)
return(score)
}

#Run analysis
score <- opi_score(textdoc = policing_dtd, metric = 5, fun = myfun)
#print results
print(score)


}
\references{
(1) Adepeju, M. and Jimoh, F. (2021). An
Analytical Framework for Measuring Inequality in the
Public Opinions on Policing – Assessing the impacts
of COVID-19 Pandemic using Twitter Data.
https://doi.org/10.31235/osf.io/c32qh
(2) Malshe, A. (2019) Data Analytics Applications.
Online book available at:
https://ashgreat.github.io/analyticsAppBook/index.html.
Date accessed: 15th December 2020.
(3) Taboada, M.et al. (2011).
Lexicon-based methods for sentiment analysis. Computational
linguistics, 37(2), pp.267-307.
(4) Lowe, W. et al. (2011).
Scaling policy preferences from coded political texts.
Legislative studies quarterly, 36(1), pp.123-155.
(5) Razorfish (2009) Fluent: The Razorfish Social Influence
Marketing Report. Accessed: 24th February, 2021.
(6) Nielsen, F. A. (2011), “A new ANEW: Evaluation of a word
list for sentiment analysis in microblogs”, Proceedings of the
ESWC2011 Workshop on 'Making Sense of Microposts': Big things
come in small packages (2011) 93-98.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{policing_dtd}
\alias{policing_dtd}
\title{Twitter posts on police/policing}
\format{
A dataframe containing one variable
\itemize{
\item text: individual text records
}
}
\usage{
policing_dtd
}
\description{
A text document (an DTD) containing twitter posts
(for an anonymous geographical location 'A') on police/policing.
The DTD also includes
posts that express sentiments on policing in relation to
the COVID-19 pandemic (Secondary subject B)
}
\keyword{datasets}
