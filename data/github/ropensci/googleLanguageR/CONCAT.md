# googleLanguageR - R client for the Google Translation API, Natural Language API, Speech-to-Text API and Text-to-Speech API

[![CRAN](https://www.r-pkg.org/badges/version/googleLanguageR)](https://cran.r-project.org/package=googleLanguageR)
[![Build
Status](https://travis-ci.org/ropensci/googleLanguageR.png?branch=master)](https://travis-ci.org/ropensci/googleLanguageR)
[![codecov.io](https://codecov.io/github/ropensci/googleLanguageR/coverage.svg?branch=master)](http://codecov.io/github/ropensci/googleLanguageR?branch=master)
[![](https://badges.ropensci.org/127_status.svg)](https://github.com/ropensci/onboarding/issues/127)

## Language tools for R via Google Machine Learning APIs

Read the [introduction blogpost on rOpenSci's blog](https://ropensci.org/blog/2017/10/03/googlelanguager/)

This package contains functions for analysing language through the
[Google Cloud Machine Learning
APIs](https://cloud.google.com/products/machine-learning/)

Note all are paid services, you will need to provide your credit card
details for your own Google Project to use them.

The package can be used by any user who is looking to take advantage of
Google’s massive dataset to train these machine learning models. Some
applications include:

  - Translation of speech into another language text, via speech-to-text
    then translation and having the results spoen back to you
  - Talking Shiny apps
  - Identification of sentiment within text, such as from Twitter feeds
  - Pulling out the objects of a sentence, to help classify texts and
    get metadata links from Wikipedia about them.

The applications of the API results could be relevant to business or
researchers looking to scale text analysis.

## Google Natural Language API

> Google Natural Language API reveals the structure and meaning of text
> by offering powerful machine learning models in an easy to use REST
> API. You can use it to extract information about people, places,
> events and much more, mentioned in text documents, news articles or
> blog posts. You can also use it to understand sentiment about your
> product on social media or parse intent from customer conversations
> happening in a call center or a messaging app.

Read more [on the Google Natural Language
API](https://cloud.google.com/natural-language/)

## Google Cloud Translation API

> Google Cloud Translation API provides a simple programmatic interface
> for translating an arbitrary string into any supported language.
> Translation API is highly responsive, so websites and applications can
> integrate with Translation API for fast, dynamic translation of source
> text from the source language to a target language (e.g. French to
> English).

Read more [on the Google Cloud Translation
Website](https://cloud.google.com/translate/)

## Google Cloud Speech-to-Text API

> Google Cloud Speech-to-Text API enables you to convert audio to text
> by applying neural network models in an easy to use API. The API
> recognizes over 80 languages and variants, to support your global user
> base. You can transcribe the text of users dictating to an
> application’s microphone or enable command-and-control through voice
> among many other use cases.

Read more [on the Google Cloud Speech
Website](https://cloud.google.com/speech/)

## Google Cloud Text-to-Speech API

> Google Cloud Text-to-Speech enables developers to synthesize
> natural-sounding speech with 30 voices, available in multiple
> languages and variants. It applies DeepMind’s groundbreaking research
> in WaveNet and Google’s powerful neural networks to deliver the
> highest fidelity possible. With this easy-to-use API, you can create
> lifelike interactions with your users, across many applications and
> devices.

Read more [on the Google Cloud Text-to-Speech
Website](https://cloud.google.com/text-to-speech/)

## Installation

1.  Create a [Google API Console
    Project](https://cloud.google.com/resource-manager/docs/creating-managing-projects)
2.  Within your project, add a [payment method to the
    project](https://support.google.com/cloud/answer/6293589)
3.  Within your project, check the relevant APIs are activated

<!-- end list -->

  - [Google Natural Language
    API](https://console.cloud.google.com/apis/api/language.googleapis.com/overview)
  - [Google Cloud Translation
    API](https://console.cloud.google.com/apis/api/translate.googleapis.com/overview)
  - [Google Cloud Speech-to-Text
    API](https://console.cloud.google.com/apis/api/speech.googleapis.com/overview)
  - [Google Cloud Text-to-Speech API](https://console.cloud.google.com/apis/library/texttospeech.googleapis.com)

<!-- end list -->

4.  [Generate a service account
    credential](https://cloud.google.com/storage/docs/authentication#generating-a-private-key)
    as a JSON file
5.  Return to R, and install the official release via
    `install.packages("googleLanguageR")`, or the development version
    with `remotes::install_github("ropensci/googleLanguageR")`
    
### Docker image

Some Docker images are publicly available.  In general `gcr.io/gcer-public/googleLanguageR:$BRANCH_NAME` carries that GitHub branch's version.

* `gcr.io/gcer-public/googleLanguageR:CRAN` - the latest CRAN version  [![CRAN](https://www.r-pkg.org/badges/version/googleLanguageR)](https://cran.r-project.org/package=googleLanguageR)
* `gcr.io/gcer-public/googleLanguageR:master` - latest GitHub master version [![Build
Status](https://travis-ci.org/ropensci/googleLanguageR.png?branch=master)](https://travis-ci.org/ropensci/googleLanguageR)
* `gcr.io/gcer-public/googleLanguageR:feature` - a feature branch from GitHub

## Usage

### Authentication

The best way to authenticate is to use an environment file. See
`?Startup`. I usually place this in my home directory. (e.g. if using
RStudio, click on `Home` in the file explorer, create a new `TEXT` file
and call it `.Renviron`)

Set the file location of your download Google Project JSON file in a
`GL_AUTH` argument:

    #.Renviron
    GL_AUTH=location_of_json_file.json

Then, when you load the library you should auto-authenticate:

``` r
library(googleLanguageR)
```

You can also authenticate directly using the `gl_auth` function pointing
at your JSON auth file:

``` r
library(googleLanguageR)
gl_auth("location_of_json_file.json")
```

You can then call the APIs via the functions:

  - `gl_nlp()` - Natural Langage API
  - `gl_speech()` - Cloud Speech-to-Text API
  - `gl_translate()` - Cloud Translation API
  - `gl_talk()` - Cloud Text-to-Speech API

## Natural Language API

The Natural Language API returns natural language understanding
technolgies. You can call them individually, or the default is to return
them all. The available returns are:

  - *Entity analysis* - Finds named entities (currently proper names and
    common nouns) in the text along with entity types, salience,
    mentions for each entity, and other properties. If possible, will
    also return metadata about that entity such as a Wikipedia URL. If
    using the **v1beta2** endpoint this also includes sentiment for each
    entity.
  - *Syntax* - Analyzes the syntax of the text and provides sentence
    boundaries and tokenization along with part of speech tags,
    dependency trees, and other properties.
  - *Sentiment* - The overall sentiment of the text, represented by a
    magnitude `[0, +inf]` and score between `-1.0` (negative sentiment)
    and `1.0` (positive sentiment).

### Demo for Entity Analysis

You can pass a vector of text which will call the API for each element.
The return is a list of responses, each response being a list of tibbles
holding the different types of
analysis.

``` r
texts <- c("to administer medicince to animals is frequently a very difficult matter, and yet sometimes it's necessary to do so", 
         "I don't know how to make a text demo that is sensible")
nlp_result <- gl_nlp(texts)

# two results of lists of tibbles
str(nlp_result, max.level = 2)
```

See more examples and details [on the
website](http://code.markedmondson.me/googleLanguageR/articles/nlp.html)
or via `vignette("nlp", package = "googleLanguageR")`

## Google Translation API

You can detect the language via `gl_translate_detect`, or translate and
detect language via `gl_translate`

Note this is a lot more refined than the free version on Google’s
translation
website.

``` r
text <- "to administer medicine to animals is frequently a very difficult matter, and yet sometimes it's necessary to do so"
## translate British into Danish
gl_translate(text, target = "da")$translatedText
```

See more examples and details [on the
website](http://code.markedmondson.me/googleLanguageR/articles/translation.html)
or via `vignette("translate", package = "googleLanguageR")`

## Google Cloud Speech-to-Text API

The Cloud Speech-to-Text API provides audio transcription. Its
accessible via the `gl_speech` function.

A test audio file is installed with the package which reads:

> “To administer medicine to animals is frequently a very difficult
> matter, and yet sometimes it’s necessary to do so”

The file is sourced from the University of Southampton’s speech
detection (`http://www-mobile.ecs.soton.ac.uk/newcomms/`) group and is
fairly difficult for computers to parse, as we see below:

``` r
## get the sample source file
test_audio <- system.file("woman1_wb.wav", package = "googleLanguageR")

## its not perfect but...:)
gl_speech(test_audio)$transcript


    ## # A tibble: 1 x 2
    ##   transcript                                                    confidence
    ##   <chr>                                                         <chr>     
    ## 1 to administer medicine to animals is frequency of very diffi… 0.9180294
```

See more examples and details [on the
website](http://code.markedmondson.me/googleLanguageR/articles/speech.html)
or via `vignette("speech", package = "googleLanguageR")`

## Google Cloud Text-to-Speech API

The Cloud Text-to-Speech API turns text into talk audio files. Its
accessible via the `gl_talk` function.

To use, supply your text to the function:

``` r
gl_talk("This is a talking computer.  Hello Dave.")
```

See more examples and details [on the
website](http://code.markedmondson.me/googleLanguageR/articles/text-to-speech.html)
or via `vignette("text-to-speech", package =
"googleLanguageR")`

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# 0.3.0.9000

* ...

# 0.3.0

* Improved error handling for vectorised `gl_nlp()` (#55)
* `gl_nlp()`'s classifyText returns list of data.frames, not data.frame
* Fix `gl_nlp` when `nlp_type='classifyText'`
* `customConfig` available for `gl_speech`
* Add support for SSML for `gl_talk()` (#66)
* Add support for device profiles for `gl_talk()` (#67)
* Add support for tuneR wave objects in `gl_speech()` - (#62 thanks @muschellij2)
* Add check for file size for audio source - (#62 thanks @muschellij2)

# 0.2.0

* Added an example Shiny app that calls the Speech API
* Fixed bug where cbind of any missing API content raised an error (#28)
* Add Google text to speech via `gl_talk()` (#39)
* Add classify text endpoint for `gl_nlp()` (#20)

# 0.1.1

* Fix bug where `gl_speech()` only returned first few seconds of translation when asynch (#23)
* CRAN version carries stable API, GitHub version for beta features

# 0.1.0


* Natural language API via `gl_nlp`
* Speech annotation via `gl_speech_recognise`
* Translation detection and performance via `gl_translate_detect` and `gl_translate`
* Vectorised support for inputs
* Translating HTML support
* Tibble outputs
# Contributing to googleLanguageR

Thank you for your interest in contributing to this project!

To run unit tests, saved API calls are cached in the `tests/testthat/mock` folder.  These substitute for an API call to Google and avoid authentication.

The API calls to the Cloud Speech API and translation varies slightly, so test success is judged if the string is within 10 characters of the test string. For this test then, the `stringdist` package is needed (under the package `Suggests`)

To run integration tests that hit the API, you will need to add your own authentication service JSON file from Google Cloud projects.  Save this file to your computer and then set an environment variable `GL_AUTH` pointing to the file location. If not present, (such as on CRAN or Travis) the integration tests will be skipped.

You will need to enable the following APIs:

* [Google Cloud Speech API](https://console.developers.google.com/apis/api/speech.googleapis.com/overview)
* [Google Cloud Natural Language API](https://console.developers.google.com/apis/api/language.googleapis.com/overview)
* [Google Cloud Translation API](https://console.developers.google.com/apis/api/translate.googleapis.com/overview)

To create new mock files, it needs to fully load the package so do:

```
remotes::install_github("ropensci/googleLanguageR")
setwd("tests/testthat")
source("test_unit.R")
```

## Contributor Covenant Code of Conduct.

* Any contributors should be doing so for the joy of creating and sharing and advancing knowledge.  Treat each other with the respect this deserves.  
* The main language is English. 
* The community will tolerate everything but intolerance. It should be assumed everyone is trying to be tolerant until they repeatedly prove otherwise. 
* Don't break any laws or copyrights during contributions, credit sources where it is due. 
## Test environments
* local OS X install, R 3.6.3
* ubuntu 14.04.5 LTS R 3.6.3
* rhub - Windows Server 2008 R2 SP1, R-devel, 32/64 bit, R 3.6.3

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.

# Google Cloud Speech API Shiny app

This is a demo on using the [Cloud Speech API](https://cloud.google.com/speech/) with Shiny. 

It uses `library(tuneR)` to process the audio file, and a JavaScript audio library from [Web Audio Demos](https://webaudiodemos.appspot.com/AudioRecorder/index.html) to capture the audio in your browser.

You can also optionally send your transcription to the [Cloud Translation API](https://cloud.google.com/translate/)

The results are then spoken back to you using the `gl_talk()` functions.

## Screenshot

![](babelfish.png)

