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
monkeylearn
===========

  [![Project Status: Abandoned – Initial development has started, but there has not yet been a stable, usable release; the project has been abandoned and the author(s) do not intend on continuing development.](https://www.repostatus.org/badges/latest/abandoned.svg)](https://www.repostatus.org/#abandoned)
[![](https://badges.ropensci.org/45_status.svg)](https://github.com/ropensci/onboarding/issues/45)

This package has been archived.

This R package is an interface to the [MonkeyLearn API](http://docs.monkeylearn.com/article/api-reference/). MonkeyLearn is a Machine Learning platform on the cloud that allows software companies and developers to easily extract actionable data from text. :monkey:
# monkeylearn 0.2.0

* New functions `monkey_classify()` and `monkey_extract()` that:
    * Accept as input both a vector and a dataframe and named column
    * Always return a tibble explicitly relating each input to its classification, allowing for the removal of the MD5 hash
    * Have an `unnest` flag to unnest the output (turn 1 row per input into 1 row per output)
    * Have a `.keep_all` flag to retain other columns if input is a dataframe
    * Coerce `NULL` values and empty vectors returned from MonkeyLearn to `NA`s
    * Include inputs that could not be processed as `NA`s in the output
    * Message the first 20 indices of inputs that are not sent to the API (these now include `NA` and `NULL` values as well as empty strings)
    * Message the currently processing batch

* Bug fixes and improvements to `monkeylearn_classify()` and `monkeylearn_extract()`
    * `monkeylearn_classify()` can now accept `params`
    * Fix to messaging when unable to connect to MonkeyLearn API
    * Default texts per request is set to 200 now (the recommended number), rather than 20
    * Addition of message suggesting that users switch to newer functions

* Implementation of `ratelimitr`. Creation and documentation of two environment variables allowing smarter rate handling when querying the MonkeyLearn API.

* Creation of `pkgdown` website

* Programmatic test coverage to re-use common tests for multiple circumstances.

* Use of a `cowsay` monkey when verbose=TRUE.


# monkeylearn 0.1.3

* Better states the dependency on tibble, it is tibble >= 1.2.

* Better handles blank text in input, outputs an empty tibble and a warning if the request is only blank, and a message if only parts of the request are blank.


# monkeylearn 0.1.2

* Disables HTTP2 for now because of a bug for Windows users. Fix by Jeroen Ooms.

# monkeylearn 0.1.1

* Added a `NEWS.md` file to track changes to the package.



## Test environments
* local x86_64-w64-mingw32/x64 install, R 3.3.1
* Ubuntu 12.04 (on Travis CI), R devel, release and oldrel
* Windows on Appveyor CI (stable, patched and oldrel)

## R CMD check results

0 errors | 0 warnings | 0 note

## Release summary


* New functions `monkey_classify()` and `monkey_extract()` that:
    * Accept as input both a vector and a dataframe and named column
    * Always return a tibble explicitly relating each input to its classification, allowing for the removal of the MD5 hash
    * Have an `unnest` flag to unnest the output (turn 1 row per input into 1 row per output)
    * Have a `.keep_all` flag to retain other columns if input is a dataframe
    * Coerce `NULL` values and empty vectors returned from MonkeyLearn to `NA`s
    * Include inputs that could not be processed as `NA`s in the output
    * Message the first 20 indices of inputs that are not sent to the API (these now include `NA` and `NULL` values as well as empty strings)
    * Message the currently processing batch

* Bug fixes and improvements to `monkeylearn_classify()` and `monkeylearn_extract()`
    * `monkeylearn_classify()` can now accept `params`
    * Fix to messaging when unable to connect to MonkeyLearn API
    * Default texts per request is set to 200 now (the recommended number), rather than 20
    * Addition of message suggesting that users switch to newer functions

* Implementation of `ratelimitr`. Creation and documentation of two environment variables allowing smarter rate handling when querying the MonkeyLearn API.

* Creation of `pkgdown` website

* Programmatic test coverage to re-use common tests for multiple circumstances.

* Use of a `cowsay` monkey when verbose=TRUE.


---
---
title: "monkeylearn, a R Package for Natural Language Processing Using Monkeylearn Existing Modules"
author: "M. Salmon, A. Dobbyn"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{intro}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, echo = FALSE, warning=FALSE, message=FALSE}
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
pat <- Sys.getenv("MONKEYLEARN_KEY")
IS_THERE_KEY <- (pat != "")
NOT_CRAN <- ifelse(IS_THERE_KEY, NOT_CRAN, FALSE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  purl = NOT_CRAN,
  eval = NOT_CRAN
)
```

# Intro

This package is an interface to the [MonkeyLearn API](http://docs.monkeylearn.com/article/api-reference/). MonkeyLearn is a Machine Learning platform on the cloud that allows software companies and developers to easily extract actionable data from text.

The goal of the package is not to support machine learning algorithms development with R or the API, but only to *reap the benefits of the existing modules on Monkeylearn*. Therefore, there are only two functions, one for using *extractors*, and one for using *classifiers*. The difference between extractors and classifiers is that extractors output information about words, whereas classifiers output information about each text as a whole. Named entity recognition is an extraction task, whereas assigning a topic to a text is a classification task.

## Setup

To get an API key for MonkeyLearn, register at http://monkeylearn.com/. Note that MonkeyLearn supports registration through GitHub, which makes the registration process really easy. For ease of use, save your API key as an environment variable as described at http://stat545.com/bit003_api-key-env-var.html. You might also want to use the `usethis::edit_r_environ()` function to modify .Renviron.

All functions of the package will conveniently look for your API key using `Sys.getenv("MONKEYLEARN_KEY")` so if your API key is an environment variable called "MONKEYLEARN\_KEY" you don't need to input it manually.

Please also create a "MONKEYLEARN\_PLAN" environment variable indicating whether your [Monkeylearn plan](https://app.monkeylearn.com/main/my-account/tab/change-plan/) is "free", "team", "business" or "custom". If you do not indicate it by default it will be "free" with a message. If your plan is "custom" you'll need a third environment variable "MONKEYLEARN\_RATE" indicating the maximum amount of requests per minute that you can make to the API. If you do not indicate it, by default it will be 120 with a message.

## So many monkeys/functions

The packages exports `monkeylearn_classify`, `monkey_classify`, `monkeylearn_extract`, `monkey_extract`. The `monkey_` functions are the newer and better ones, so if you don't have legacy code, just start using those!

For inspiration beyond this vignette, you can see external examples of the package in action [on this page](http://ropensci.github.io/monkeylearn/). In particular you'll find examples using the older set of functions but we now recommend using `monkey_extract` and `monkey_classify`, see more later in the vignette. 

# Extract

## A first example 

```{r, message = FALSE}
library(monkeylearn)
library(magrittr)

text <- "In the 19th century, the major European powers had gone to great lengths to maintain a balance of power throughout Europe, resulting in the existence of a complex network of political and military alliances throughout the continent by 1900.[7] These had started in 1815, with the Holy Alliance between Prussia, Russia, and Austria. Then, in October 1873, German Chancellor Otto von Bismarck negotiated the League of the Three Emperors (German: Dreikaiserbund) between the monarchs of Austria-Hungary, Russia and Germany."
output <- monkey_extract(input = text,
                         extractor_id = "ex_isnnZRbS")
output
attr(output, "headers")
```

## Parameters

If the documentation of the extractor you use states it has parameters, you can pass them as a named list, see below.

```{r}
text <- "A panel of Goldman Sachs employees spent a recent Tuesday night at the
Columbia University faculty club trying to convince a packed room of potential
recruits that Wall Street, not Silicon Valley, was the place to be for computer
scientists.\n\n The Goldman employees knew they had an uphill battle. They were
fighting against perceptions of Wall Street as boring and regulation-bound and
Silicon Valley as the promised land of flip-flops, beanbag chairs and million-dollar
stock options.\n\n Their argument to the room of technologically inclined students
was that Wall Street was where they could find far more challenging, diverse and,
yes, lucrative jobs working on some of the worlds most difficult technical problems."

output <- monkey_extract(text,
                        extractor_id = "ex_y7BPYzNG",
                        params = list(max_keywords = 3))
output
output2 <- monkey_extract(text,
                          extractor_id = "ex_y7BPYzNG",
                          params = list(max_keywords = 1))
output2
attr(output2, "headers")
```

## How to find extractors?

You can find extractors and their IDs, including extractors for text in Spanish, at https://app.monkeylearn.com/main/explore 

There is no endpoint for automatically finding all extractors, but if you find one in the website you particularly like and use a lot in your language and application, you could choose to save its id as an environment variable as explained [here]( http://stat545.com/bit003_api-key-env-var.html). Reading about extractors on the website will give you a good overview of their characteristics and original application.

Here are a few ones for text in English:

* [Entity extractor](https://app.monkeylearn.com/extraction/extractors/ex_isnnZRbS/tab/description-tab), `extractor_id = "ex_isnnZRbS"` (used in the first example). Extract Entities from text using Named Entity Recognition (NER). NER labels sequences of words in a text which are the names of things, such as person and company names. This implementation labels 3 classes: PERSON, ORGANIZATION and LOCATION. This NER tagger is implemented using Conditional Random Field (CRF) sequence models.

* [Keyword extractor](https://app.monkeylearn.com/extraction/extractors/ex_y7BPYzNG/tab/description-tab), `extractor_id = "ex_y7BPYzNG"`. Extract keywords from text in English. Keywords can be compounded by one or more words and are defined as the important topics in your content and can be used to index data, generate tag clouds or for searching. This keyword extraction algorithm employs statistical algorithms and natural language processing technology to analyze your content and identify the relevant keywords.

```{r, message = FALSE}
text <- "A panel of Goldman Sachs employees spent a recent Tuesday night at the Columbia University faculty club trying to convince a packed room of potential recruits that Wall Street, not Silicon Valley, was the place to be for computer scientists.

The Goldman employees knew they had an uphill battle. They were fighting against perceptions of Wall Street as boring and regulation-bound and Silicon Valley as the promised land of flip-flops, beanbag chairs and million-dollar stock options.

Their argument to the room of technologically inclined students was that Wall Street was where they could find far more challenging, diverse and, yes, lucrative jobs working on some of the world’s most difficult technical problems.

“Whereas in other opportunities you might be considering, it is working one type of data or one type of application, we deal in hundreds of products in hundreds of markets, with thousands or tens of thousands of clients, every day, millions of times of day worldwide,” Afsheen Afshar, a managing director at Goldman Sachs, told the students."

monkey_extract(text, extractor_id = "ex_y7BPYzNG")
```

* [Useful data extractor](https://app.monkeylearn.com/extraction/extractors/ex_dqRio5sG/tab/description-tab), `extractor_id = "ex_dqRio5sG"`. Extract useful data from text. This algorithm can be used to detect many different useful data: links, phones, ips, prices, times, emails, bitcoin addresses, dates, ipv6s, hex colors and credit cards.

When using this extractor, the format of the API output is a bit different than for other extractors, see below how the output looks like.

```{r, message = FALSE}
text <- "Hi, my email is john@example.com and my credit card is 4242-4242-4242-4242 so you can charge me with $10. My phone number is 15555 9876. We can get in touch on April 16, at 10:00am"
text2 <- "Hi, my email is mary@example.com and my credit card is 4242-4232-4242-4242. My phone number is 16655 9876. We can get in touch on April 16, at 10:00am"

monkey_extract(c(text, text2), extractor_id = "ex_dqRio5sG", unnest = TRUE)
```


# Classify

## A first example

```{r, message = FALSE}
text1 <- "my dog is an avid rice eater"
text2 <- "i want to buy an iphone"
request <- c(text1, text2)

monkey_classify(request, classifier_id = "cl_oFKL5wft")
```
## How to find classifiers?

You can find classifiers and their IDs at https://app.monkeylearn.com/main/explore or you can use the `monkeylearn_classifiers` function, choosing to show all classifiers or only the private ones with `private = TRUE`. The first column of the resulting data.frame is the `classifier_id` to be used in `monkeylearn_classify`.

```{r}
monkeylearn_classifiers(private = FALSE)
```

Here are a few other examples:

* [Language detection](https://app.monkeylearn.com/categorizer/projects/cl_oJNMkt2V/tab/main-tab), `classifier_id = "cl_oJNMkt2V"`. Detect language in text. New languages were added for a total of 48 different languages arranged in language families.

```{r, message = FALSE}
text1 <- "Hauràs de dirigir-te al punt de trobada del grup al que et vulguis unir."
text2 <- "i want to buy an iphone"
text3 <- "Je déteste ne plus avoir de dentifrice."
request <- c(text1, text2, text3)

monkey_classify(request, classifier_id = "cl_oJNMkt2V")
```

* [Profanity and abuse detection](https://app.monkeylearn.com/categorizer/projects/cl_KFXhoTdt/tab/main-tab), `classifier_id = "cl_KFXhoTdt"`.

```{r, message = FALSE}
text1 <- "I think this is awesome."
text2 <- "Holy shit! You did great!"
request <- c(text1, text2)

monkey_classify(request, classifier_id = "cl_KFXhoTdt")
```

* [General topic classifier](https://app.monkeylearn.com/categorizer/projects/cl_5icAVzKR/tab/), `classifier_id = "cl_5icAVzKR"`.

```{r, message = FALSE}
text1 <- "Let me tell you about my dog and my cat. They are really friendly and like going on walks. They both like chasing mice."
text2 <- "My first R package was probably a disaster but I keep learning how to program."
request <- c(text1, text2)
monkey_classify(request, classifier_id = "cl_5icAVzKR")

```


# Get what you paid for

Monkeylearn offers a different service based on your current plan, that is, "free", "team" or "business". These plans will both influence your _rate limiting_ (how fast?) and your _query limiting_ (how many queries?). See https://monkeylearn.com/pricing/. Thanks to your MONKEYLEARN_PLAN environment variable, the rate will be handled automatically thanks to [`ratelimitr`](https://github.com/tarakc02/ratelimitr).

## Check the number of remaining calls

After each call to a function you can check how many calls to the API you can still make  using `attr(output, "headers")$x.query.limit.remaining` and `attr(output, "headers")$x.query.limit.limit`. The period after which `attr(output, "headers")$x.query.limit.remaining` depends on your subscription and is not included in the output.



# Fit `monkeylearn` into your pipeline!

You can:

* Send a vector of texts *or* a dataframe and a named column (unquoted)
* Output either a nested or unnested dataframe
    * Nested = 1 row per input; unnested = 1 row per output
* This output
    * Relates each input text to its (usually) multiple classifications/extractions
    * Retains a record of inputs that could not be classified/extracted (e.g., empty strings)
* Batch requests


## In a bit more detail

You can classify or extract a vector or dataframe of texts while relating the original input text to its classifications. This is important, because the input:output relationship may not always (and in fact, is not usually) 1:1. These functions retain the tie between each `input`[^1] element and all of its output elements.

```{r monkey_input}
input <- c("Emma Woodhouse, handsome, clever, and rich, with a comfortable home",     
 "and happy disposition, seemed to unite some of the best blessings of",  
 "existence; and had lived nearly twenty-one years in the world with very", 
 "little to distress or vex her.",                                          
 "",                   # <--- note the empty string!                                                   
 "She was the youngest of the two daughters of a most affectionate,",       
 "indulgent father; and had, in consequence of her sister's marriage, been",
 "mistress of his house from a very early period. Her mother had died",     
 "too long ago for her to have more than an indistinct remembrance of",     
 "her caresses; and her place had been supplied by an excellent woman as",  
 "governess, who had fallen little short of a mother in affection.")
```

That is true even if you have inputs that cannot be processed. For instance, empty string and `NA` input elements are not sent to the API for classification/extraction. (You'll get a warning of this if `verbose = TRUE`.) We've got one above to illustrate and elements that returned no classifications/extractions are included in the resulting dataframe. This way you'll know which inputs could not be processed.

```{r monkey_output}
(output <- monkey_classify(input, unnest = FALSE))
```

<br>

If there are more than 20 empty inputs, we save your console by messaging only the first 20 indices.

```{r very_empty_input}
(very_empty_input <- rep("", 25) %>% c(input) %>% sample())
```


Since the entire original input is represented in the output, if you need to find all of the empty inputs you can easily filter the output to all of the rows containing empty strings.
```{r}
monkey_classify(very_empty_input, unnest = FALSE)
```


### Configuring the Output

The default output is a nested dataframe with the same number of rows as your input dataframe or the same length as your input vector, depending on which one you sent in. 

Let's take a look at the `res` output column. 
```{r}
output$res
```

You can easily choose an unnested output by setting the **unnest flag** to TRUE (which it is by default) to get one row per classification/extraction. 

```{r unnest_true}
(output_unnested <- monkey_classify(input, verbose = FALSE, unnest = TRUE))
```

We could have gotten the same result by sending in a dataframe and a named column. If a dataframe is supplied input column is not renamed to `req` as it is when input is a vector; the original column name is retained.

```{r compare_df}
input_df <- tibble::tibble(text = input) 
output_df_unnested <- monkey_classify(input_df, text, unnest = TRUE, verbose = FALSE) %>% 
    dplyr::rename(req = text)

testthat::expect_equal(output_unnested, output_df_unnested)
```

<br>

If the input is a dataframe, setting the `.keep_all` option to TRUE allows you to retain all input columns. If FALSE, only the column you specify for classification will be retained. 

```{r keep_all}
sw <- dplyr::starwars %>% 
  dplyr::select(name, height) %>% 
  dplyr::sample_n(nrow(input_df))

sw_input_df <- input_df %>% 
  dplyr::bind_cols(sw)

sw_input_df %>% monkey_classify(text, unnest = FALSE, verbose = FALSE)
```


### Batching

Retaining the relationship between input and output doesn't mean you'll need to send requests one-by-one. **Batch requests** by setting the `texts_per_req` value which governs the number of texts that are sent per request. Per the [MonkeyLearn documentation](http://help.monkeylearn.com/frequently-asked-questions/queries/can-i-classify-or-extract-more-than-one-text-with-one-api-request), the maximum we recommend sending at once is 200 requests. 

If `texts_per_req` is NULL, the default, we try to optimize the response time from the API by setting `texts_per_req` to 200 when your input has more than 200 texts or to the length of the `input` if you've got fewer. You'll see a significant speedup by batching your requests this way. However, batching doesn't save you on queries; a batch of 150 texts still uses up 150 queries. 

These functions also include some more verbose **progress reporting**, letting you know what batch you're on out of the total, and which texts are set to be processed in that batch.

```{r one_by_one, warning=FALSE}
one_by_one <- system.time(output <- monkey_classify(input, texts_per_req = 1))
```

```{r batch_of_five, warning=FALSE}
batch_of_five <- system.time(output <- monkey_classify(input, texts_per_req = 5))
```

How much does sending 5 texts in a batch vs. 1 text improve our processing time?
```{r speedup}
(speedup <- one_by_one[1] / batch_of_five[1])
```


A 3-4x speedup isn't so bad! Worth keeping in mind that if you need the blazing fast speeds you might consider upgrading to a higher MonkeyLearn price tier. 

<br>





***

<br>


# Meta

* Please [report any issues or bugs](https://github.com/ropensci/monkeylearn/issues).
* License: GPL
* Get citation information for `monkeylearn` in R doing `citation(package = 'monkeylearn')`
* Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
* This package is part of the [rOpenSci project](https://ropensci.org/).


[^1]: Thanks to [Julia Silge](https://juliasilge.com/)'s fantastic [`janeaustenr`](https://github.com/juliasilge/janeaustenr) package for this text!

---
title: "monkeylearn, a R Package for Natural Language Processing Using Monkeylearn Existing Modules"
author: "M. Salmon, A. Dobbyn"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{intro}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, echo = FALSE, warning=FALSE, message=FALSE}
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
pat <- Sys.getenv("MONKEYLEARN_KEY")
IS_THERE_KEY <- (pat != "")
NOT_CRAN <- ifelse(IS_THERE_KEY, NOT_CRAN, FALSE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  purl = NOT_CRAN,
  eval = NOT_CRAN
)
```

# Intro

This package is an interface to the [MonkeyLearn API](http://docs.monkeylearn.com/article/api-reference/). MonkeyLearn is a Machine Learning platform on the cloud that allows software companies and developers to easily extract actionable data from text.

The goal of the package is not to support machine learning algorithms development with R or the API, but only to *reap the benefits of the existing modules on Monkeylearn*. Therefore, there are only two functions, one for using *extractors*, and one for using *classifiers*. The difference between extractors and classifiers is that extractors output information about words, whereas classifiers output information about each text as a whole. Named entity recognition is an extraction task, whereas assigning a topic to a text is a classification task.

## Setup

To get an API key for MonkeyLearn, register at http://monkeylearn.com/. Note that MonkeyLearn supports registration through GitHub, which makes the registration process really easy. For ease of use, save your API key as an environment variable as described at http://stat545.com/bit003_api-key-env-var.html. You might also want to use the `usethis::edit_r_environ()` function to modify .Renviron.

All functions of the package will conveniently look for your API key using `Sys.getenv("MONKEYLEARN_KEY")` so if your API key is an environment variable called "MONKEYLEARN\_KEY" you don't need to input it manually.

Please also create a "MONKEYLEARN\_PLAN" environment variable indicating whether your [Monkeylearn plan](https://app.monkeylearn.com/main/my-account/tab/change-plan/) is "free", "team", "business" or "custom". If you do not indicate it by default it will be "free" with a message. If your plan is "custom" you'll need a third environment variable "MONKEYLEARN\_RATE" indicating the maximum amount of requests per minute that you can make to the API. If you do not indicate it, by default it will be 120 with a message.

## So many monkeys/functions

The packages exports `monkeylearn_classify`, `monkey_classify`, `monkeylearn_extract`, `monkey_extract`. The `monkey_` functions are the newer and better ones, so if you don't have legacy code, just start using those!

For inspiration beyond this vignette, you can see external examples of the package in action [on this page](http://docs.ropensci.org/monkeylearn/). In particular you'll find examples using the older set of functions but we now recommend using `monkey_extract` and `monkey_classify`, see more later in the vignette. 

# Extract

## A first example 

```{r, message = FALSE}
library(monkeylearn)
library(magrittr)

text <- "In the 19th century, the major European powers had gone to great lengths to maintain a balance of power throughout Europe, resulting in the existence of a complex network of political and military alliances throughout the continent by 1900.[7] These had started in 1815, with the Holy Alliance between Prussia, Russia, and Austria. Then, in October 1873, German Chancellor Otto von Bismarck negotiated the League of the Three Emperors (German: Dreikaiserbund) between the monarchs of Austria-Hungary, Russia and Germany."
output <- monkey_extract(input = text,
                         extractor_id = "ex_isnnZRbS")
output
attr(output, "headers")
```

## Parameters

If the documentation of the extractor you use states it has parameters, you can pass them as a named list, see below.

```{r}
text <- "A panel of Goldman Sachs employees spent a recent Tuesday night at the
Columbia University faculty club trying to convince a packed room of potential
recruits that Wall Street, not Silicon Valley, was the place to be for computer
scientists.\n\n The Goldman employees knew they had an uphill battle. They were
fighting against perceptions of Wall Street as boring and regulation-bound and
Silicon Valley as the promised land of flip-flops, beanbag chairs and million-dollar
stock options.\n\n Their argument to the room of technologically inclined students
was that Wall Street was where they could find far more challenging, diverse and,
yes, lucrative jobs working on some of the worlds most difficult technical problems."

output <- monkey_extract(text,
                        extractor_id = "ex_y7BPYzNG",
                        params = list(max_keywords = 3))
output
output2 <- monkey_extract(text,
                          extractor_id = "ex_y7BPYzNG",
                          params = list(max_keywords = 1))
output2
attr(output2, "headers")
```

## How to find extractors?

You can find extractors and their IDs, including extractors for text in Spanish, at https://app.monkeylearn.com/main/explore 

There is no endpoint for automatically finding all extractors, but if you find one in the website you particularly like and use a lot in your language and application, you could choose to save its id as an environment variable as explained [here]( http://stat545.com/bit003_api-key-env-var.html). Reading about extractors on the website will give you a good overview of their characteristics and original application.

Here are a few ones for text in English:

* [Entity extractor](https://app.monkeylearn.com/extraction/extractors/ex_isnnZRbS/tab/description-tab), `extractor_id = "ex_isnnZRbS"` (used in the first example). Extract Entities from text using Named Entity Recognition (NER). NER labels sequences of words in a text which are the names of things, such as person and company names. This implementation labels 3 classes: PERSON, ORGANIZATION and LOCATION. This NER tagger is implemented using Conditional Random Field (CRF) sequence models.

* [Keyword extractor](https://app.monkeylearn.com/extraction/extractors/ex_y7BPYzNG/tab/description-tab), `extractor_id = "ex_y7BPYzNG"`. Extract keywords from text in English. Keywords can be compounded by one or more words and are defined as the important topics in your content and can be used to index data, generate tag clouds or for searching. This keyword extraction algorithm employs statistical algorithms and natural language processing technology to analyze your content and identify the relevant keywords.

```{r, message = FALSE}
text <- "A panel of Goldman Sachs employees spent a recent Tuesday night at the Columbia University faculty club trying to convince a packed room of potential recruits that Wall Street, not Silicon Valley, was the place to be for computer scientists.

The Goldman employees knew they had an uphill battle. They were fighting against perceptions of Wall Street as boring and regulation-bound and Silicon Valley as the promised land of flip-flops, beanbag chairs and million-dollar stock options.

Their argument to the room of technologically inclined students was that Wall Street was where they could find far more challenging, diverse and, yes, lucrative jobs working on some of the world’s most difficult technical problems.

“Whereas in other opportunities you might be considering, it is working one type of data or one type of application, we deal in hundreds of products in hundreds of markets, with thousands or tens of thousands of clients, every day, millions of times of day worldwide,” Afsheen Afshar, a managing director at Goldman Sachs, told the students."

monkey_extract(text, extractor_id = "ex_y7BPYzNG")
```

* [Useful data extractor](https://app.monkeylearn.com/extraction/extractors/ex_dqRio5sG/tab/description-tab), `extractor_id = "ex_dqRio5sG"`. Extract useful data from text. This algorithm can be used to detect many different useful data: links, phones, ips, prices, times, emails, bitcoin addresses, dates, ipv6s, hex colors and credit cards.

When using this extractor, the format of the API output is a bit different than for other extractors, see below how the output looks like.

```{r, message = FALSE}
text <- "Hi, my email is john@example.com and my credit card is 4242-4242-4242-4242 so you can charge me with $10. My phone number is 15555 9876. We can get in touch on April 16, at 10:00am"
text2 <- "Hi, my email is mary@example.com and my credit card is 4242-4232-4242-4242. My phone number is 16655 9876. We can get in touch on April 16, at 10:00am"

monkey_extract(c(text, text2), extractor_id = "ex_dqRio5sG", unnest = TRUE)
```


# Classify

## A first example

```{r, message = FALSE}
text1 <- "my dog is an avid rice eater"
text2 <- "i want to buy an iphone"
request <- c(text1, text2)

monkey_classify(request, classifier_id = "cl_oFKL5wft")
```
## How to find classifiers?

You can find classifiers and their IDs at https://app.monkeylearn.com/main/explore or you can use the `monkeylearn_classifiers` function, choosing to show all classifiers or only the private ones with `private = TRUE`. The first column of the resulting data.frame is the `classifier_id` to be used in `monkeylearn_classify`.

```{r}
monkeylearn_classifiers(private = FALSE)
```

Here are a few other examples:

* [Language detection](https://app.monkeylearn.com/categorizer/projects/cl_oJNMkt2V/tab/main-tab), `classifier_id = "cl_oJNMkt2V"`. Detect language in text. New languages were added for a total of 48 different languages arranged in language families.

```{r, message = FALSE}
text1 <- "Hauràs de dirigir-te al punt de trobada del grup al que et vulguis unir."
text2 <- "i want to buy an iphone"
text3 <- "Je déteste ne plus avoir de dentifrice."
request <- c(text1, text2, text3)

monkey_classify(request, classifier_id = "cl_oJNMkt2V")
```

* [Profanity and abuse detection](https://app.monkeylearn.com/categorizer/projects/cl_KFXhoTdt/tab/main-tab), `classifier_id = "cl_KFXhoTdt"`.

```{r, message = FALSE}
text1 <- "I think this is awesome."
text2 <- "Holy shit! You did great!"
request <- c(text1, text2)

monkey_classify(request, classifier_id = "cl_KFXhoTdt")
```

* [General topic classifier](https://app.monkeylearn.com/categorizer/projects/cl_5icAVzKR/tab/), `classifier_id = "cl_5icAVzKR"`.

```{r, message = FALSE}
text1 <- "Let me tell you about my dog and my cat. They are really friendly and like going on walks. They both like chasing mice."
text2 <- "My first R package was probably a disaster but I keep learning how to program."
request <- c(text1, text2)
monkey_classify(request, classifier_id = "cl_5icAVzKR")

```


# Get what you paid for

Monkeylearn offers a different service based on your current plan, that is, "free", "team" or "business". These plans will both influence your _rate limiting_ (how fast?) and your _query limiting_ (how many queries?). See https://monkeylearn.com/pricing/. Thanks to your MONKEYLEARN_PLAN environment variable, the rate will be handled automatically thanks to [`ratelimitr`](https://github.com/tarakc02/ratelimitr).

## Check the number of remaining calls

After each call to a function you can check how many calls to the API you can still make  using `attr(output, "headers")$x.query.limit.remaining` and `attr(output, "headers")$x.query.limit.limit`. The period after which `attr(output, "headers")$x.query.limit.remaining` depends on your subscription and is not included in the output.



# Fit `monkeylearn` into your pipeline!

You can:

* Send a vector of texts *or* a dataframe and a named column (unquoted)
* Output either a nested or unnested dataframe
    * Nested = 1 row per input; unnested = 1 row per output
* This output
    * Relates each input text to its (usually) multiple classifications/extractions
    * Retains a record of inputs that could not be classified/extracted (e.g., empty strings)
* Batch requests


## In a bit more detail

You can classify or extract a vector or dataframe of texts while relating the original input text to its classifications. This is important, because the input:output relationship may not always (and in fact, is not usually) 1:1. These functions retain the tie between each `input`[^1] element and all of its output elements.

```{r monkey_input}
input <- c("Emma Woodhouse, handsome, clever, and rich, with a comfortable home",     
 "and happy disposition, seemed to unite some of the best blessings of",  
 "existence; and had lived nearly twenty-one years in the world with very", 
 "little to distress or vex her.",                                          
 "",                   # <--- note the empty string!                                                   
 "She was the youngest of the two daughters of a most affectionate,",       
 "indulgent father; and had, in consequence of her sister's marriage, been",
 "mistress of his house from a very early period. Her mother had died",     
 "too long ago for her to have more than an indistinct remembrance of",     
 "her caresses; and her place had been supplied by an excellent woman as",  
 "governess, who had fallen little short of a mother in affection.")
```

That is true even if you have inputs that cannot be processed. For instance, empty string and `NA` input elements are not sent to the API for classification/extraction. (You'll get a warning of this if `verbose = TRUE`.) We've got one above to illustrate and elements that returned no classifications/extractions are included in the resulting dataframe. This way you'll know which inputs could not be processed.

```{r monkey_output}
(output <- monkey_classify(input, unnest = FALSE))
```

<br>

If there are more than 20 empty inputs, we save your console by messaging only the first 20 indices.

```{r very_empty_input}
(very_empty_input <- rep("", 25) %>% c(input) %>% sample())
```


Since the entire original input is represented in the output, if you need to find all of the empty inputs you can easily filter the output to all of the rows containing empty strings.
```{r}
monkey_classify(very_empty_input, unnest = FALSE)
```


### Configuring the Output

The default output is a nested dataframe with the same number of rows as your input dataframe or the same length as your input vector, depending on which one you sent in. 

Let's take a look at the `res` output column. 
```{r}
output$res
```

You can easily choose an unnested output by setting the **unnest flag** to TRUE (which it is by default) to get one row per classification/extraction. 

```{r unnest_true}
(output_unnested <- monkey_classify(input, verbose = FALSE, unnest = TRUE))
```

We could have gotten the same result by sending in a dataframe and a named column. If a dataframe is supplied input column is not renamed to `req` as it is when input is a vector; the original column name is retained.

```{r compare_df}
input_df <- tibble::tibble(text = input) 
output_df_unnested <- monkey_classify(input_df, text, unnest = TRUE, verbose = FALSE) %>% 
    dplyr::rename(req = text)

testthat::expect_equal(output_unnested, output_df_unnested)
```

<br>

If the input is a dataframe, setting the `.keep_all` option to TRUE allows you to retain all input columns. If FALSE, only the column you specify for classification will be retained. 

```{r keep_all}
sw <- dplyr::starwars %>% 
  dplyr::select(name, height) %>% 
  dplyr::sample_n(nrow(input_df))

sw_input_df <- input_df %>% 
  dplyr::bind_cols(sw)

sw_input_df %>% monkey_classify(text, unnest = FALSE, verbose = FALSE)
```


### Batching

Retaining the relationship between input and output doesn't mean you'll need to send requests one-by-one. **Batch requests** by setting the `texts_per_req` value which governs the number of texts that are sent per request. Per the [MonkeyLearn documentation](http://help.monkeylearn.com/frequently-asked-questions/queries/can-i-classify-or-extract-more-than-one-text-with-one-api-request), the maximum we recommend sending at once is 200 requests. 

If `texts_per_req` is NULL, the default, we try to optimize the response time from the API by setting `texts_per_req` to 200 when your input has more than 200 texts or to the length of the `input` if you've got fewer. You'll see a significant speedup by batching your requests this way. However, batching doesn't save you on queries; a batch of 150 texts still uses up 150 queries. 

These functions also include some more verbose **progress reporting**, letting you know what batch you're on out of the total, and which texts are set to be processed in that batch.

```{r one_by_one, warning=FALSE}
one_by_one <- system.time(output <- monkey_classify(input, texts_per_req = 1))
```

```{r batch_of_five, warning=FALSE}
batch_of_five <- system.time(output <- monkey_classify(input, texts_per_req = 5))
```

How much does sending 5 texts in a batch vs. 1 text improve our processing time?
```{r speedup}
(speedup <- one_by_one[1] / batch_of_five[1])
```


A 3-4x speedup isn't so bad! Worth keeping in mind that if you need the blazing fast speeds you might consider upgrading to a higher MonkeyLearn price tier. 

<br>





***

<br>


# Meta

* Please [report any issues or bugs](https://github.com/ropensci/monkeylearn/issues).
* License: GPL
* Get citation information for `monkeylearn` in R doing `citation(package = 'monkeylearn')`
* Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
* This package is part of the [rOpenSci project](https://ropensci.org/).


[^1]: Thanks to [Julia Silge](https://juliasilge.com/)'s fantastic [`janeaustenr`](https://github.com/juliasilge/janeaustenr) package for this text!

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/monkey_extract.R
\name{monkey_extract}
\alias{monkey_extract}
\title{Monkeylearn extract from a dataframe column or vector of texts}
\usage{
monkey_extract(input, col = NULL, key = monkeylearn_key(quiet = TRUE),
  extractor_id = "ex_isnnZRbS", params = NULL, texts_per_req = NULL,
  unnest = TRUE, .keep_all = TRUE, verbose = TRUE, ...)
}
\arguments{
\item{input}{A dataframe or vector of texts (each text smaller than 50kB)}

\item{col}{If input is a dataframe, the unquoted name of the character column containing text to extract from}

\item{key}{The API key}

\item{extractor_id}{The ID of the extractor}

\item{params}{Parameters for the module as a named list.}

\item{texts_per_req}{Number of texts to be processed per requests. Minimum value is the number of texts in input; max is 200, as per
[Monkeylearn documentation](docs.monkeylearn.com/article/api-reference/). If NULL, we default to 200, or, if there are fewer than 200 texts, the length of the input.}

\item{unnest}{Should the output column be unnested?}

\item{.keep_all}{If \code{input} is a dataframe, should non-\code{col} columns be retained in the output?}

\item{verbose}{Whether to output messages about batch requests and progress of processing.}

\item{...}{Other arguments}
}
\value{
A data.frame (tibble) with the cleaned input (empty strings removed) and a new column, nested by default, containing the extraction for that particular row.
Attribute is a data.frame (tibble) "headers" including the number of remaining queries as "x.query.limit.remaining".
}
\description{
Independent extractions for each row of a dataframe using the Monkeylearn extractor modules
}
\details{
Find IDs of extractors using \url{https://app.monkeylearn.com/main/explore}.

This function relates the rows in your original dataframe or elements in your vector to an extraction particular to that row.
This allows you to know which row of your original dataframe is associated with which extraction.
Each row of the dataframe is extracted separately from all of the others, but the number of extractions a particular input row
is assigned may vary (unless you specify a fixed number of outputs in \code{params}).

The \code{texts_per_req} parameter simply specifies the number of rows to feed the API at a time; it does not lump these together
for extraction as a group. Varying this parameter does not affect the final output, but does affect speed: one batched request of
x texts is faster than x single-text requests:
\url{http://help.monkeylearn.com/frequently-asked-questions/queries/can-i-classify-or-extract-more-than-one-text-with-one-api-request}.
Even if batched, each text still counts as one query, so batching does not save you on hits to the API.
See the [Monkeylearn API docs](docs.monkeylearn.com/article/api-reference/) for more details.

You can check the number of calls you can still make in the API using \code{attr(output, "headers")$x.query.limit.remaining}
and \code{attr(output, "headers")$x.query.limit.limit}.

Find IDs of extractors using \url{https://app.monkeylearn.com/main/explore}.
Within the free plan, you can make up to 20 requests per minute.

You can use batch to send up to 200 texts to be analyzed within the API
(classification or extraction) with each request.
So for example, if you need to analyze 6000 tweets,
instead of doing 6000 requests to the API, you can use batch to send 30 requests,
each request with 200 tweets.
The function automatically makes these batch calls and waits if there is a throttle limit error,
but you might want to control the process yourself using several calls to the function.

You can check the number of calls you can still make in the API using \code{attr(output, "headers")$x.query.limit.remaining}
and \code{attr(output, "headers")$x.query.limit.limit}.
}
\examples{
\dontrun{
text <- "In the 19th century, the major European powers had gone to great lengths
to maintain a balance of power throughout Europe, resulting in the existence of
a complex network of political and military alliances throughout the continent by 1900.[7]
These had started in 1815, with the Holy Alliance between Prussia, Russia, and Austria.
Then, in October 1873, German Chancellor Otto von Bismarck negotiated the League of
the Three Emperors (German: Dreikaiserbund) between the monarchs of Austria-Hungary,
Russia and Germany."
output <- monkeylearn_extract(request = text)
output


# Example with parameters
text <- "A panel of Goldman Sachs employees spent a recent Tuesday night at the
Columbia University faculty club trying to convince a packed room of potential
recruits that Wall Street, not Silicon Valley, was the place to be for computer
scientists.\\n\\n The Goldman employees knew they had an uphill battle. They were
fighting against perceptions of Wall Street as boring and regulation-bound and
Silicon Valley as the promised land of flip-flops, beanbag chairs and million-dollar
stock options.\\n\\n Their argument to the room of technologically inclined students
was that Wall Street was where they could find far more challenging, diverse and,
yes, lucrative jobs working on some of the worlds most difficult technical problems."

output <- monkey_extract(text,
                            extractor_id = "ex_y7BPYzNG",
                            params = list(max_keywords = 3,
                            use_company_names = 1))
attr(output, "headers")}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractor.R
\name{monkeylearn_extract}
\alias{monkeylearn_extract}
\title{monkeylearn_extract}
\usage{
monkeylearn_extract(request, key = monkeylearn_key(quiet = TRUE),
  extractor_id = "ex_isnnZRbS", texts_per_req = 200, verbose = TRUE,
  params = NULL)
}
\arguments{
\item{request}{A vector of characters (each text smaller than 50kB)}

\item{key}{The API key}

\item{extractor_id}{The ID of the extractor}

\item{texts_per_req}{Number of texts to be fed through per request (max 200). Does not affect output, but may affect speed of processing.}

\item{verbose}{Whether to output messages about batch requests}

\item{params}{Parameters for the module as a named list. See the second example.}
}
\value{
A data.frame with the results whose attribute is a data.frame (tibble) "headers" including the number of remaining queries as "x.query.limit.remaining".
Both data.frames include a column with the (list of) md5 checksum(s) of the corresponding text(s) computed using the \code{digest digest} function.
}
\description{
Access to Monkeylearn extractors modules
}
\details{
Find IDs of extractors using \url{https://app.monkeylearn.com/main/explore}.
Within the free plan, you can make up to 20 requests per minute.

You can use batch to send up to 200 texts to be analyzed within the API
(classification or extraction) with each request.
So for example, if you need to analyze 6000 tweets,
instead of doing 6000 requests to the API, you can use batch to send 30 requests,
each request with 200 tweets.
The function automatically makes these batch calls and waits if there is a throttle limit error,
but you might want to control the process yourself using several calls to the function.

You can check the number of calls you can still make in the API using \code{attr(output, "headers")$x.query.limit.remaining}
and \code{attr(output, "headers")$x.query.limit.limit}.
}
\examples{
\dontrun{
text <- "In the 19th century, the major European powers had gone to great lengths
to maintain a balance of power throughout Europe, resulting in the existence of
 a complex network of political and military alliances throughout the continent by 1900.[7]
  These had started in 1815, with the Holy Alliance between Prussia, Russia, and Austria.
  Then, in October 1873, German Chancellor Otto von Bismarck negotiated the League of
   the Three Emperors (German: Dreikaiserbund) between the monarchs of Austria-Hungary,
    Russia and Germany."
output <- monkeylearn_extract(request = text)
output
# example with parameters
text <- "A panel of Goldman Sachs employees spent a recent Tuesday night at the
Columbia University faculty club trying to convince a packed room of potential
recruits that Wall Street, not Silicon Valley, was the place to be for computer
scientists.\\n\\n The Goldman employees knew they had an uphill battle. They were
fighting against perceptions of Wall Street as boring and regulation-bound and
Silicon Valley as the promised land of flip-flops, beanbag chairs and million-dollar
stock options.\\n\\n Their argument to the room of technologically inclined students
was that Wall Street was where they could find far more challenging, diverse and,
yes, lucrative jobs working on some of the worlds most difficult technical problems."

output <- monkeylearn_extract(text,
                              extractor_id = "ex_y7BPYzNG",
                              params = list(max_keywords = 3,
                                            use_company_names = 1))
attr(output, "headers")}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/monkeylearn-package.R
\docType{package}
\name{monkeylearn-package}
\alias{monkeylearn}
\alias{monkeylearn-package}
\title{monkeylearn: Accesses the Monkeylearn API for Text Classifiers and Extractors}
\description{
Allows using some services of Monkeylearn <http://monkeylearn.com/> which is
a Machine Learning platform on the cloud for text analysis (classification and extraction).
}
\seealso{
Useful links:
\itemize{
  \item \url{http://github.com/ropensci/monkeylearn}
  \item \url{http://docs.ropensci.org/monkeylearn/}
  \item Report bugs at \url{http://github.com/ropensci/monkeylearn/issues}
}

}
\author{
\strong{Maintainer}: Maëlle Salmon \email{maelle.salmon@yahoo.se} (0000-0002-2815-0399)

Authors:
\itemize{
  \item Amanda Dobbyn \email{amanda.e.dobbyn@gmail.com}
}

Other contributors:
\itemize{
  \item Thomas Leeper (Thomas Leeper reviewed the package for rOpenSci, see https://github.com/ropensci/onboarding/issues/45) [reviewer]
  \item Jeroen Ooms [contributor]
  \item rOpenSci (https://ropensci.org/) [funder]
  \item Earlybird Software (http://earlybird.co/) [funder]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify.R
\name{monkeylearn_classify}
\alias{monkeylearn_classify}
\title{monkeylearn_classify}
\usage{
monkeylearn_classify(request, key = monkeylearn_key(quiet = TRUE),
  classifier_id = "cl_oFKL5wft", texts_per_req = 200, verbose = TRUE,
  params = NULL)
}
\arguments{
\item{request}{A vector of characters (each text smaller than 50kB)}

\item{key}{The API key}

\item{classifier_id}{The ID of the classifier}

\item{texts_per_req}{Number of texts to be fed through per request (max 200). Does not affect output, but may affect speed of processing.}

\item{verbose}{Whether to output messages about batch requests}

\item{params}{Parameters for the module as a named list. See the second example.}
}
\value{
A data.frame (tibble) with the results whose attribute is a data.frame (tibble) "headers" including the number of remaining queries as "x.query.limit.remaining".
Both data.frames include a column with the (list of) md5 checksum(s) of the corresponding text(s) computed using the \code{digest digest} function.
}
\description{
Access to Monkeylearn classifiers modules
}
\details{
Find IDs of classifiers using \url{https://app.monkeylearn.com/main/explore}.

You can use batch to send up to 200 texts to be analyzed within the API
(classification or extraction) with each request.
So for example, if you need to analyze 6000 tweets,
instead of doing 6000 requests to the API, you can use batch to send 30 requests,
each request with 200 tweets.
The function automatically makes these batch calls and waits if there is a throttle limit error,
but you might want to control the process yourself using several calls to the function.

You can check the number of calls you can still make in the API using \code{attr(output, "headers")$x.query.limit.remaining}
and \code{attr(output, "headers")$x.query.limit.limit}.
}
\examples{
\dontrun{
text1 <- "my dog is an avid rice eater"
text2 <- "i want to buy an iphone"
request <- c(text1, text2)
output <- monkeylearn_classify(request)
output
attr(output, "headers")}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classifiers.R
\name{monkeylearn_classifiers}
\alias{monkeylearn_classifiers}
\title{monkeylearn_classifiers}
\usage{
monkeylearn_classifiers(private = FALSE, key = monkeylearn_key(quiet =
  TRUE))
}
\arguments{
\item{private}{default is FALSE, whether to show private modules only instead of private and public modules}

\item{key}{The API key}
}
\value{
A data.frame (tibble) with details about the
classifiers including their classifier_id which should be used in
\code{monkeylearn_classify}.
}
\description{
List of Monkeylearn classifiers modules
}
\details{
If you don't have any private modules,
\code{monkeylearn_classifiers(private = TRUE)} returns an empty data.frame.
}
\examples{
\dontrun{
monkeylearn_classifiers(private = FALSE)
monkeylearn_classifiers(private = TRUE)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{monkeylearn_key}
\alias{monkeylearn_key}
\title{Retrieve Monkeylearn API key}
\usage{
monkeylearn_key(quiet = TRUE)
}
\value{
An Monkeylearn API Key
}
\description{
Retrieve Monkeylearn API key
}
\details{
Looks in env var \code{MONKEYLEARN_KEY}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/monkey_classify.R
\name{monkey_classify}
\alias{monkey_classify}
\title{Monkeylearn classify from a dataframe column or vector of texts}
\usage{
monkey_classify(input, col = NULL, key = monkeylearn_key(quiet = TRUE),
  classifier_id = "cl_oFKL5wft", params = NULL, texts_per_req = NULL,
  unnest = TRUE, .keep_all = TRUE, verbose = TRUE, ...)
}
\arguments{
\item{input}{A dataframe or vector of texts (each text smaller than 50kB)}

\item{col}{If input is a dataframe, the unquoted name of the character column containing text to classify}

\item{key}{The API key}

\item{classifier_id}{The ID of the classifier}

\item{params}{Parameters for the module as a named list.}

\item{texts_per_req}{Number of texts to be processed per requests. Minimum value is the number of texts in input; max is 200, as per
[Monkeylearn documentation](docs.monkeylearn.com/article/api-reference/). If NULL, we default to 200, or, if there are fewer than 200 texts, the length of the input.}

\item{unnest}{Should the output column be unnested?}

\item{.keep_all}{If \code{input} is a dataframe, should non-\code{col} columns be retained in the output?}

\item{verbose}{Whether to output messages about batch requests and progress of processing.}

\item{...}{Other arguments}
}
\value{
A data.frame (tibble) with the cleaned input (empty strings removed) and a new column, nested by default, containing the classification for that particular row.
Attribute is a data.frame (tibble) "headers" including the number of remaining queries as "x.query.limit.remaining".
}
\description{
Independent classifications for each row of a dataframe using the Monkeylearn classifiers modules
}
\details{
Find IDs of classifiers using \url{https://app.monkeylearn.com/main/explore}.

This function relates the rows in your original dataframe or elements in your vector to a classification particular to that row.
This allows you to know which row of your original dataframe is associated with which classification.
Each row of the dataframe is classified separately from all of the others, but the number of classifications a particular input row
is assigned may vary (unless you specify a fixed number of outputs in \code{params}).

The \code{texts_per_req} parameter simply specifies the number of rows to feed the API at a time; it does not lump these together
for classification as a group. Varying this parameter does not affect the final output, but does affect speed: one batched request of
x texts is faster than x single-text requests:
\url{http://help.monkeylearn.com/frequently-asked-questions/queries/can-i-classify-or-extract-more-than-one-text-with-one-api-request}.
Even if batched, each text still counts as one query, so batching does not save you on hits to the API.
See the [Monkeylearn API docs](docs.monkeylearn.com/article/api-reference/) for more details.

You can check the number of calls you can still make in the API using \code{attr(output, "headers")$x.query.limit.remaining}
and \code{attr(output, "headers")$x.query.limit.limit}.
}
\examples{
\dontrun{
text1 <- "Hauràs de dirigir-te al punt de trobada del grup al que et vulguis unir."
text2 <- "i want to buy an iphone"
text3 <- "Je déteste ne plus avoir de dentifrice."
text_4 <- "I hate not having any toothpaste."
request_df <- tibble::as_tibble(list(txt = c(text1, text2, text3, text_4)))
monkey_classify(request_df, txt, texts_per_req = 2, unnest = TRUE)
attr(output, "headers")}

}
