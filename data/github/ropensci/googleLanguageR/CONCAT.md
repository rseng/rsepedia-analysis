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

---
title: "Google Cloud Text-to-Speech API"
author: "Mark Edmondson"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Google Cloud Text-to-Speech API}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

Google Cloud Text-to-Speech enables developers to synthesize natural-sounding speech with 30 voices, available in multiple languages and variants. It applies DeepMind’s groundbreaking research in WaveNet and Google’s powerful neural networks to deliver the highest fidelity possible. With this easy-to-use API, you can create lifelike interactions with your users, across many applications and devices.

Read more [on the Google Cloud Text-to-Speech Website](https://cloud.google.com/text-to-speech/)

The Cloud Text-to-Speech API turns text into sound files of the spoken words.  Its accessible via the `gl_talk` function.

Arguments include:

* `input` - The text to turn into speech
* `output` Where to save the speech audio file
* `languageCode` The language of the voice as a [`BCP-47` language tag](https://tools.ietf.org/html/bcp47)
* `name` Name of the voice, see list via `gl_talk_languages()` or [online](https://cloud.google.com/text-to-speech/docs/voices) for supported voices.  If not set, then the service will choose a voice based on `languageCode` and `gender`.
* `gender` The gender of the voice, if available
* `audioEncoding` Format of the requested audio stream - can be a choice of `.wav`, `.mp3` or `.ogg`
* `speakingRate` Speaking rate/speed
* `pitch` Speaking pitch
* `volumeGainDb` Volumne gain in dB
* `sampleRateHertz` Sample rate for returned audio

## Returned structure

The API returns an audio file which is saved to the location specified in `output` - by default this is `output.wav` - if you don't rename this file it will be overwritten by the next API call.  

It is advised to set the appropriate file extension if you change the audio encoding (e.g. to one of `.wav`, `.mp3` or `.ogg`) so audio payers recognise the file format. 

## Talk Languages

The API can talk several different languages, with more being added over time.  You can get a current list via the function `gl_talk_languages()` or [online](https://cloud.google.com/text-to-speech/docs/voices)

```r
gl_talk_languages()
# A tibble: 32 x 4
   languageCodes name             ssmlGender naturalSampleRateHertz
   <chr>         <chr>            <chr>                       <int>
 1 es-ES         es-ES-Standard-A FEMALE                      24000
 2 ja-JP         ja-JP-Standard-A FEMALE                      22050
 3 pt-BR         pt-BR-Standard-A FEMALE                      24000
 4 tr-TR         tr-TR-Standard-A FEMALE                      22050
 5 sv-SE         sv-SE-Standard-A FEMALE                      22050
 6 nl-NL         nl-NL-Standard-A FEMALE                      24000
 7 en-US         en-US-Wavenet-A  MALE                        24000
 8 en-US         en-US-Wavenet-B  MALE                        24000
 9 en-US         en-US-Wavenet-C  FEMALE                      24000
10 en-US         en-US-Wavenet-D  MALE                        24000
```

If you are looking a specific language, specify that in the function call e.g. to see only Spanish (`es`)
voices issue:

```r
gl_talk_languages(languageCode = "es")
# A tibble: 1 x 4
  languageCodes name             ssmlGender naturalSampleRateHertz
  <chr>         <chr>            <chr>                       <int>
1 es-ES         es-ES-Standard-A FEMALE                      24000
```

You can then specify that voice when calling the API via the `name` argument, which overrides the `gender` and `languageCode` argument:

```r
gl_talk("Hasta la vista", name = "es-ES-Standard-A")
```

Otherwise, specify your own `gender` and `languageCode` and the voice will be picked for you:

```r
gl_talk("Would you like a cup of tea?", gender = "FEMALE", languageCode = "en-GB")
```

Some languages are not yet supported, such as Danish.  The API will return an error in those cases. 

## Support for SSML

Support is also included for Speech Synthesis Markup Language (SSML) - more details on using this to insert pauses, sounds and breaks in your audio can be found here: `https://cloud.google.com/text-to-speech/docs/ssml`

To use, send in your SSML markup around the text you want to talk and set `inputType= "ssml"`:

```r
# using SSML
gl_talk('<speak>The <say-as interpret-as=\"characters\">SSML</say-as>
  standard <break time=\"1s\"/>is defined by the
  <sub alias=\"World Wide Web Consortium\">W3C</sub>.</speak>',
  inputType =  "ssml")
```

## Effect Profiles

You can output audio files that are optimised for playing on various devices. 

To use audio profiles, supply a character vector of the available audio profiles listed here: `https://cloud.google.com/text-to-speech/docs/audio-profiles` - the audio profiles are applied in the order given. 

For instance `effectsProfileIds="wearable-class-device"` will optimise output for smart watches, `effectsProfileIds=c("wearable-class-device","telephony-class-application")` will apply sound filters optimised for smart watches, then telephonic devices.

```r
# using effects profiles
gl_talk("This sounds great on headphones",
        effectsProfileIds = "headphone-class-device")
```

## Browser Speech player

Creating and clicking on the audio file to play it can be a bit of a drag, so you also have a function that will play the audio file for you, launching via the browser.  This can be piped via the tidyverse's `%>%`

```r
library(magrittr)
gl_talk("This is my audio player") %>% gl_talk_player()

## non-piped equivalent
gl_talk_player(gl_talk("This is my audio player"))
```

The `gl_talk_player()` creates a HTML file called `player.html` in your working directory by default.

### Using with Shiny

You can do this in Shiny too, which is demonstrated in the [example Shiny app](https://github.com/ropensci/googleLanguageR/tree/master/inst/shiny/capture_speech) included with the package.

Click the link for a video tutorial on how to [integrate text-to-speech into a Shiny app](https://www.youtube.com/watch?v=Ny0e7vHFu6o&t=8s) - the demo uses text-to-speech to [talk through a user's Google Analytics statistics](https://github.com/MarkEdmondson1234/verbal_ga_shiny).

<iframe width="560" height="315" src="https://www.youtube.com/embed/Ny0e7vHFu6o" frameborder="0" allow="autoplay; encrypted-media" allowfullscreen></iframe>

A shiny module has been created to help integrate text-to-speech into your Shiny apps, demo in the video above and below:

```r
library(shiny)
library(googleLanguageR)  # assume auto auth setup

ui <- fluidPage(
  gl_talk_shinyUI("talk")
)

server <- function(input, output, session){

     transcript <- reactive({
        paste("This is a demo talking Shiny app!")
     })

     callModule(gl_talk_shiny, "talk", transcript = transcript)
}


shinyApp(ui = ui, server = server)
```


---
title: "Google Cloud Translation API"
author: "Mark Edmondson"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Google Cloud Translation API}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

The Google Cloud Translation API provides a simple programmatic interface for translating an arbitrary string into any supported language. Translation API is highly responsive, so websites and applications can integrate with Translation API for fast, dynamic translation of source text from the source language to a target language (e.g. French to English). 

Read more [on the Google Cloud Translation Website](https://cloud.google.com/translate/)

You can detect the language via `gl_translate_detect`, or translate and detect language via `gl_translate`

### Language Translation

Translate text via `gl_translate`.  Note this is a lot more refined than the free version on Google's translation website.

```r
library(googleLanguageR)

text <- "to administer medicince to animals is frequently a very difficult matter, and yet sometimes it's necessary to do so"
## translate British into Danish
gl_translate(text, target = "da")$translatedText
```

You can choose the target language via the argument `target`.  The function will automatically detect the language if you do not define an argument `source`.  This function which will also detect the langauge. As it costs the same as `gl_translate_detect`, its usually cheaper to detect and translate in one step.

You can pass a vector of text which will first be attempted to translate in one API call - if that fails due to being greater than the API limits, it will attempt again but vectorising the API calls.  This will result in more calls and be slower, but cost the same as you are charged per character translated, not per API call. 

#### HTML support

You can also supply web HTML and select the `format='html'` which will handle HTML tags to give you a cleaner translation. 

Consider removing anything not needed to be translated first, such as JavaScript and CSS scripts using the tools of `rvest` - an example is shown below:

```r
# translate webpages
library(rvest)
library(googleLanguageR)

my_url <- "http://www.dr.dk/nyheder/indland/greenpeace-facebook-og-google-boer-foelge-apples-groenne-planer"

## in this case the content to translate is in css select .wcms-article-content
read_html(my_url) %>% # read html
  html_node(css = ".wcms-article-content") %>%   # select article content
  html_text %>% # extract text
  gl_translate(format = "html") %>% # translate with html flag
  dplyr::select(translatedText) # show translatedText column of output tibble

```


### Language Detection

This function only detects the language:

```r
## which language is this?
gl_translate_detect("katten sidder på måtten")
```

The more text it has, the better.  And it helps if its not Danish...

It may be better to use [`cld2`](https://github.com/ropensci/cld2) to translate offline first, to avoid charges if the translation is unnecessary (e.g. already in English).  You could then verify online for more uncertain cases.

```r
cld2::detect_language("katten sidder på måtten")
```

### Translation API limits

The API limits in three ways: characters per day, characters per 100 seconds, and API requests per 100 seconds. All can be set in the API manager in Google Cloud console: `https://console.developers.google.com/apis/api/translate.googleapis.com/quotas`

The library will limit the API calls for the characters and API requests per 100 seconds. The API will automatically retry if you are making requests too quickly, and also pause to make sure you only send `100000` characters per 100 seconds. 
---
title: "Google Natural Language API"
author: "Mark Edmondson"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Google Natural Language API}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

The Google Natural Language API reveals the structure and meaning of text by offering powerful machine learning models in an easy to use REST API. You can use it to extract information about people, places, events and much more, mentioned in text documents, news articles or blog posts. You can also use it to understand sentiment about your product on social media or parse intent from customer conversations happening in a call center or a messaging app. 

Read more [on the Google Natural Language API](https://cloud.google.com/natural-language/)

The Natural Language API returns natural language understanding technologies.  You can call them individually, or the default is to return them all.  The available returns are:

* *Entity analysis* - Finds named entities (currently proper names and common nouns) in the text along with entity types, salience, mentions for each entity, and other properties.  If possible, will also return metadata about that entity such as a Wikipedia URL. 
* *Syntax* - Analyzes the syntax of the text and provides sentence boundaries and tokenization along with part of speech tags, dependency trees, and other properties.
* *Sentiment* - The overall sentiment of the text, represented by a magnitude `[0, +inf]` and score between `-1.0` (negative sentiment) and `1.0` (positive sentiment).
* *Content Classification* - Analyzes a document and returns a list of content categories that apply to the text found in the document. A complete list of content categories can be found [here](https://cloud.google.com/natural-language/docs/categories).


### Demo for Entity Analysis

You can pass a vector of text which will call the API for each element.  The return is a list of responses, each response being a list of tibbles holding the different types of analysis.

```r
library(googleLanguageR)

# random text form wikipedia
texts <- c("Norma is a small constellation in the Southern Celestial Hemisphere between Ara and Lupus, one of twelve drawn up in the 18th century by French astronomer Nicolas Louis de Lacaille and one of several depicting scientific instruments. Its name refers to a right angle in Latin, and is variously considered to represent a rule, a carpenter's square, a set square or a level. It remains one of the 88 modern constellations. Four of Norma's brighter stars make up a square in the field of faint stars. Gamma2 Normae is the brightest star with an apparent magnitude of 4.0. Mu Normae is one of the most luminous stars known, but is partially obscured by distance and cosmic dust. Four star systems are known to harbour planets. ", 
         "Solomon Wariso (born 11 November 1966 in Portsmouth) is a retired English sprinter who competed primarily in the 200 and 400 metres.[1] He represented his country at two outdoor and three indoor World Championships and is the British record holder in the indoor 4 × 400 metres relay.")
nlp_result <- gl_nlp(texts)
```

Each text has its own entry in returned tibbles

```r
str(nlp_result, max.level = 2)
List of 7
 $ sentences        :List of 2
  ..$ :'data.frame':	7 obs. of  4 variables:
  ..$ :'data.frame':	1 obs. of  4 variables:
 $ tokens           :List of 2
  ..$ :'data.frame':	139 obs. of  17 variables:
  ..$ :'data.frame':	54 obs. of  17 variables:
 $ entities         :List of 2
  ..$ :Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	52 obs. of  9 variables:
  ..$ :Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	8 obs. of  9 variables:
 $ language         : chr [1:2] "en" "en"
 $ text             : chr [1:2] "Norma is a small constellation in the Southern Celestial Hemisphere between Ara and Lupus, one of twelve drawn "| __truncated__ "Solomon Wariso (born 11 November 1966 in Portsmouth) is a retired English sprinter who competed primarily in th"| __truncated__
 $ documentSentiment:Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	2 obs. of  2 variables:
  ..$ magnitude: num [1:2] 2.4 0.1
  ..$ score    : num [1:2] 0.3 0.1
 $ classifyText     :Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	1 obs. of  2 variables:
  ..$ name      : chr "/Science/Astronomy"
  ..$ confidence: num 0.93
```

Sentence structure and sentiment:

```r
## sentences structure
nlp_result$sentences[[2]]

content
1 Solomon Wariso (born 11 November 1966 in Portsmouth) is a retired English sprinter who competed primarily in the 200 and 400 metres.[1] He represented his country at two outdoor and three indoor World Championships and is the British record holder in the indoor 4 × 400 metres relay.
  beginOffset magnitude score
1           0       0.1   0.1
```

Information on what words (tokens) are within each text:

```r
# word tokens data
str(nlp_result$tokens[[1]])
'data.frame':	139 obs. of  17 variables:
 $ content       : chr  "Norma" "is" "a" "small" ...
 $ beginOffset   : int  0 6 9 11 17 31 34 38 47 57 ...
 $ tag           : chr  "NOUN" "VERB" "DET" "ADJ" ...
 $ aspect        : chr  "ASPECT_UNKNOWN" "ASPECT_UNKNOWN" "ASPECT_UNKNOWN" "ASPECT_UNKNOWN" ...
 $ case          : chr  "CASE_UNKNOWN" "CASE_UNKNOWN" "CASE_UNKNOWN" "CASE_UNKNOWN" ...
 $ form          : chr  "FORM_UNKNOWN" "FORM_UNKNOWN" "FORM_UNKNOWN" "FORM_UNKNOWN" ...
 $ gender        : chr  "GENDER_UNKNOWN" "GENDER_UNKNOWN" "GENDER_UNKNOWN" "GENDER_UNKNOWN" ...
 $ mood          : chr  "MOOD_UNKNOWN" "INDICATIVE" "MOOD_UNKNOWN" "MOOD_UNKNOWN" ...
 $ number        : chr  "SINGULAR" "SINGULAR" "NUMBER_UNKNOWN" "NUMBER_UNKNOWN" ...
 $ person        : chr  "PERSON_UNKNOWN" "THIRD" "PERSON_UNKNOWN" "PERSON_UNKNOWN" ...
 $ proper        : chr  "PROPER" "PROPER_UNKNOWN" "PROPER_UNKNOWN" "PROPER_UNKNOWN" ...
 $ reciprocity   : chr  "RECIPROCITY_UNKNOWN" "RECIPROCITY_UNKNOWN" "RECIPROCITY_UNKNOWN" "RECIPROCITY_UNKNOWN" ...
 $ tense         : chr  "TENSE_UNKNOWN" "PRESENT" "TENSE_UNKNOWN" "TENSE_UNKNOWN" ...
 $ voice         : chr  "VOICE_UNKNOWN" "VOICE_UNKNOWN" "VOICE_UNKNOWN" "VOICE_UNKNOWN" ...
 $ headTokenIndex: int  1 1 4 4 1 4 9 9 9 5 ...
 $ label         : chr  "NSUBJ" "ROOT" "DET" "AMOD" ...
 $ value         : chr  "Norma" "be" "a" "small" ...
```

What entities within text have been identified, with optional wikipedia URL if its available.

```r
nlp_result$entities
[[1]]
# A tibble: 52 x 9
   name           type         salience mid   wikipedia_url magnitude score beginOffset mention_type
   <chr>          <chr>           <dbl> <chr> <chr>             <dbl> <dbl>       <int> <chr>       
 1 angle          OTHER         0.0133  NA    NA                  0     0           261 COMMON      
 2 Ara            ORGANIZATION  0.0631  NA    NA                  0     0            76 PROPER      
 3 astronomer     NA           NA       NA    NA                 NA    NA           144 COMMON      
 4 carpenter      PERSON        0.0135  NA    NA                  0     0           328 COMMON      
 5 constellation  OTHER         0.150   NA    NA                  0     0            17 COMMON      
 6 constellations OTHER         0.0140  NA    NA                  0.9   0.9         405 COMMON      
 7 distance       OTHER         0.00645 NA    NA                  0     0           649 COMMON      
 8 dust           OTHER         0.00645 NA    NA                  0.3  -0.3         669 COMMON      
 9 field          LOCATION      0.00407 NA    NA                  0.6  -0.6         476 COMMON      
10 French         LOCATION      0.0242  NA    NA                  0     0           137 PROPER      
# ... with 42 more rows

[[2]]
# A tibble: 8 x 9
  name                type         salience mid         wikipedia_url    magnitude score beginOffset mention_type
  <chr>               <chr>           <dbl> <chr>       <chr>                <dbl> <dbl>       <int> <chr>       
1 British             LOCATION       0.0255 NA          NA                     0     0           226 PROPER      
2 country             LOCATION       0.0475 NA          NA                     0     0           155 COMMON      
3 English             OTHER          0.0530 NA          NA                     0     0            66 PROPER      
4 Portsmouth          LOCATION       0.0530 /m/0619_    https://en.wiki…       0     0            41 PROPER      
5 record holder       PERSON         0.0541 NA          NA                     0     0           234 COMMON      
6 Solomon Wariso      ORGANIZATION   0.156  /g/120x5nf6 https://en.wiki…       0     0             0 PROPER      
7 sprinter            PERSON         0.600  NA          NA                     0     0            74 COMMON      
8 World Championships EVENT          0.0113 NA          NA                     0.1   0.1         195 PROPER      

```

Sentiment of the entire text:

```r
nlp_result$documentSentiment
# A tibble: 2 x 2
  magnitude score
      <dbl> <dbl>
1       2.4   0.3
2       0.1   0.1
```

The category for the text as defined by the list [here](https://cloud.google.com/natural-language/docs/categories).

```r
nlp_result$classifyText
# A tibble: 1 x 2
  name               confidence
  <chr>                   <dbl>
1 /Science/Astronomy       0.93
```

The language for the text:

```r
nlp_result$language
# [1] "en" "en"
```

The original passed in text, to aid with working with the output:

```r
nlp_result$text
[1] "Norma is a small constellation in the Southern Celestial Hemisphere between Ara and Lupus, one of twelve drawn up in the 18th century by French astronomer Nicolas Louis de Lacaille and one of several depicting scientific instruments. Its name refers to a right angle in Latin, and is variously considered to represent a rule, a carpenter's square, a set square or a level. It remains one of the 88 modern constellations. Four of Norma's brighter stars make up a square in the field of faint stars. Gamma2 Normae is the brightest star with an apparent magnitude of 4.0. Mu Normae is one of the most luminous stars known, but is partially obscured by distance and cosmic dust. Four star systems are known to harbour planets."
[2] "Solomon Wariso (born 11 November 1966 in Portsmouth) is a retired English sprinter who competed primarily in the 200 and 400 metres.[1] He represented his country at two outdoor and three indoor World Championships and is the British record holder in the indoor 4 × 400 metres relay."
```
---
title: "Introduction to googleLanguageR"
author: "Mark Edmondson"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to googleLanguageR}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

`googleLanguageR` contains functions for analysing language through the [Google Cloud Machine Learning APIs](https://cloud.google.com/products/machine-learning/)

Note all are paid services, you will need to provide your credit card details for your own Google Project to use them.

The package can be used by any user who is looking to take advantage of Google's massive dataset to train these machine learning models.  Some applications include:

* Translation of speech into another language text, via speech-to-text then translation
* Identification of sentiment within text, such as from Twitter feeds
* Pulling out the objects of a sentence, to help classify texts and get metadata links from Wikipedia about them.

The applications of the API results could be relevant to business or researchers looking to scale text analysis.

## Google Natural Language API

> Google Natural Language API reveals the structure and meaning of text by offering powerful machine learning models in an easy to use REST API. You can use it to extract information about people, places, events and much more, mentioned in text documents, news articles or blog posts. You can also use it to understand sentiment about your product on social media or parse intent from customer conversations happening in a call center or a messaging app. 

Read more [on the Google Natural Language API](https://cloud.google.com/natural-language/)

## Google Cloud Translation API

> Google Cloud Translation API provides a simple programmatic interface for translating an arbitrary string into any supported language. Translation API is highly responsive, so websites and applications can integrate with Translation API for fast, dynamic translation of source text from the source language to a target language (e.g. French to English). 

Read more [on the Google Cloud Translation Website](https://cloud.google.com/translate/)

## Google Cloud Speech API

> Google Cloud Speech API enables you to convert audio to text by applying neural network models in an easy to use API. The API recognizes over 80 languages and variants, to support your global user base. You can transcribe the text of users dictating to an application’s microphone or enable command-and-control through voice among many other use cases. 

Read more [on the Google Cloud Speech Website](https://cloud.google.com/speech/)

## Installation

1. Create a [Google API Console Project](https://cloud.google.com/resource-manager/docs/creating-managing-projects)
2. Within your project, add a [payment method to the project](https://support.google.com/cloud/answer/6293589)
3. Within your project, check the relevant APIs are activated
  - [Google Natural Language API](https://console.cloud.google.com/apis/api/language.googleapis.com/overview)
  - [Google Cloud Translation API](https://console.cloud.google.com/apis/api/translate.googleapis.com/overview)
  - [Google Cloud Speech API](https://console.cloud.google.com/apis/api/speech.googleapis.com/overview)
4. [Generate a service account credential](https://cloud.google.com/storage/docs/authentication#generating-a-private-key) as a JSON file
5. Return to R, and install this library via `devtools::install_github("MarkEdmondson1234/googleLanguageR")`

## Usage

### Authentication

The best way to authenticate is to use an environment file.  See `?Startup`.  I usually place this in my home directory. (e.g. if using RStudio, click on `Home` in the file explorer, create a new `TEXT` file and call it `.Renviron`)  

Set the file location of your download Google Project JSON file in a `GL_AUTH` argument:

```
#.Renviron
GL_AUTH=location_of_json_file.json
```

Then, when you load the library you should auto-authenticate:

```r
library(googleLanguageR)
# Setting scopes to https://www.googleapis.com/auth/cloud-platform
# Set any additional scopes via options(googleAuthR.scopes.selected = c('scope1', 'scope2')) before loading library.
# Successfully authenticated via location_of_json_file.json
```

You can also authenticate directly using the `gl_auth` function pointing at your JSON auth file:

```r
library(googleLanguageR)
gl_auth("location_of_json_file.json")
```

You can then call the APIs via the functions:

* `gl_nlp()` - Natural Langage API
* `gl_speech()` - Cloud Speech API
* `gl_translate()` - Cloud Translation API

---
title: "Google Cloud Speech-to-Text API"
author: "Mark Edmondson"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Google Cloud Speech-to-Text API}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

The Google Cloud Speech-to-Text API enables you to convert audio to text by applying neural network models in an easy to use API. The API recognizes over 80 languages and variants, to support your global user base. You can transcribe the text of users dictating to an application’s microphone or enable command-and-control through voice among many other use cases. 

Read more [on the Google Cloud Speech-to-Text Website](https://cloud.google.com/speech/)

The Cloud Speech API provides audio transcription.  Its accessible via the `gl_speech` function.

Arguments include:

* `audio_source` - this is a local file in the correct format, or a Google Cloud Storage URI. This can also be a `Wave` class object from the package `tuneR`
* `encoding` - the format of the sound file - `LINEAR16` is the common `.wav` format, other formats include `FLAC` and `OGG_OPUS`
* `sampleRate` - this needs to be set to what your file is recorded at.  
* `languageCode` - specify the language spoken as a [`BCP-47` language tag](https://tools.ietf.org/html/bcp47)
* `speechContexts` - you can supply keywords to help the translation with some context. 

### Returned structure

The API returns a list of two data.frame tibbles - `transcript` and `timings`.

Access them via the returned object and `$transcript` and `$timings`

```r
return <- gl_speech(test_audio, languageCode = "en-GB")

return$transcript
# A tibble: 1 x 2
#                                                                                                         transcript confidence
#                                                                                                              <chr>      <chr>
#1 to administer medicine to animals is frequently a very difficult matter and yet sometimes it's necessary to do so  0.9711006

return$timings
#   startTime endTime       word
#1         0s  0.100s         to
#2     0.100s  0.700s administer
#3     0.700s  0.700s   medicine
#4     0.700s  1.200s         to
# etc...
```

### Demo for Google Cloud Speech-to-Text API


A test audio file is installed with the package which reads:

> "To administer medicine to animals is frequently a very difficult matter, and yet sometimes it's necessary to do so"

The file is sourced from the University of Southampton's speech detection (`http://www-mobile.ecs.soton.ac.uk/`) group and is fairly difficult for computers to parse, as we see below:

```r
library(googleLanguageR)
## get the sample source file
test_audio <- system.file("woman1_wb.wav", package = "googleLanguageR")

## its not perfect but...:)
gl_speech(test_audio)$transcript

## get alternative transcriptions
gl_speech(test_audio, maxAlternatives = 2L)$transcript

gl_speech(test_audio, languageCode = "en-GB")$transcript

## help it out with context for "frequently"
gl_speech(test_audio, 
            languageCode = "en-GB", 
            speechContexts = list(phrases = list("is frequently a very difficult")))$transcript
```

### Word transcripts

The API [supports timestamps](https://cloud.google.com/speech/reference/rest/v1/speech/recognize#WordInfo) on when words are recognised. These are outputted into a second data.frame that holds three entries: `startTime`, `endTime` and the `word`.


```r
str(result$timings)
#'data.frame':	152 obs. of  3 variables:
# $ startTime: chr  "0s" "0.100s" "0.500s" "0.700s" ...
# $ endTime  : chr  "0.100s" "0.500s" "0.700s" "0.900s" ...
# $ word     : chr  "a" "Dream" "Within" "A" ...

result$timings
#     startTime endTime       word
#1          0s  0.100s          a
#2      0.100s  0.500s      Dream
#3      0.500s  0.700s     Within
#4      0.700s  0.900s          A
#5      0.900s      1s      Dream
```

## Custom configurations

You can also send in other arguments which can help shape the output, such as speaker diagrization (labelling different speakers) - to use such custom configurations create a [`RecognitionConfig`](https://cloud.google.com/speech-to-text/docs/reference/rest/v1p1beta1/RecognitionConfig) object.  This can be done via R lists which are converted to JSON via `library(jsonlite)` and an example is shown below:

```r
## Use a custom configuration
my_config <- list(encoding = "LINEAR16",
                  diarizationConfig = list(
                    enableSpeakerDiarization = TRUE,
                    minSpeakerCount = 2,
                    maxSpeakCount = 3
                  ))

# languageCode is required, so will be added if not in your custom config
gl_speech(my_audio, languageCode = "en-US", customConfig = my_config)
```

## Asynchronous calls

For speech files greater than 60 seconds of if you don't want your results straight away, set `asynch = TRUE` in the call to the API.

This will return an object of class `"gl_speech_op"` which should be used within the `gl_speech_op()` function to check the status of the task.  If the task is finished, then it will return an object the same form as the non-asynchronous case. 

```r
async <- gl_speech(test_audio, asynch = TRUE)
async
## Send to gl_speech_op() for status
## 4625920921526393240

result <- gl_speech_op(async)
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/translate.R
\name{gl_translate_languages}
\alias{gl_translate_languages}
\title{Lists languages from Google Translate API}
\usage{
gl_translate_languages(target = "en")
}
\arguments{
\item{target}{If specified, language names are localized in target language}
}
\value{
A tibble of supported languages
}
\description{
Returns a list of supported languages for translation.
}
\details{
Supported language codes, generally consisting of its ISO 639-1 identifier. (E.g. \code{'en', 'ja'}).
In certain cases, BCP-47 codes including language + region identifiers are returned (e.g. \code{'zh-TW', 'zh-CH'})
}
\examples{

\dontrun{

# default english names of languages supported
gl_translate_languages()

# specify a language code to get other names, such as Danish
gl_translate_languages("da")

}
}
\seealso{
\url{https://cloud.google.com/translate/docs/reference/languages}

Other translations: 
\code{\link{gl_translate_detect}()},
\code{\link{gl_translate}()}
}
\concept{translations}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/text-to-speech.R
\name{gl_talk_shiny}
\alias{gl_talk_shiny}
\title{Speak in Shiny module (server)}
\usage{
gl_talk_shiny(
  input,
  output,
  session,
  transcript,
  ...,
  autoplay = TRUE,
  controls = TRUE,
  loop = FALSE,
  keep_wav = FALSE
)
}
\arguments{
\item{input}{shiny input}

\item{output}{shiny output}

\item{session}{shiny session}

\item{transcript}{The (reactive) text to talk}

\item{...}{
  Arguments passed on to \code{\link[=gl_talk]{gl_talk}}
  \describe{
    \item{\code{languageCode}}{The language of the voice as a \code{BCP-47} language code}
    \item{\code{name}}{Name of the voice, see list via \link{gl_talk_languages} for supported voices.  Set to \code{NULL} to make the service choose a voice based on \code{languageCode} and \code{gender}.}
    \item{\code{gender}}{The gender of the voice, if available}
    \item{\code{audioEncoding}}{Format of the requested audio stream}
    \item{\code{speakingRate}}{Speaking rate/speed between \code{0.25} and \code{4.0}}
    \item{\code{pitch}}{Speaking pitch between \code{-20.0} and \code{20.0} in semitones.}
    \item{\code{volumeGainDb}}{Volumne gain in dB}
    \item{\code{sampleRateHertz}}{Sample rate for returned audio}
    \item{\code{inputType}}{Choose between \code{text} (the default) or SSML markup. The \code{input} text must be SSML markup if you choose \code{ssml}}
    \item{\code{effectsProfileIds}}{Optional. An identifier which selects 'audio effects' profiles that are applied on (post synthesized) text to speech. Effects are applied on top of each other in the order they are given}
  }}

\item{autoplay}{passed to the HTML audio player - default \code{TRUE} plays on load}

\item{controls}{passed to the HTML audio player - default \code{TRUE} shows controls}

\item{loop}{passed to the HTML audio player - default \code{FALSE} does not loop}

\item{keep_wav}{keep the generated wav files if TRUE.}
}
\description{
Call via \code{shiny::callModule(gl_talk_shiny, "your_id")}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/speech-to-text.R
\name{gl_speech_op}
\alias{gl_speech_op}
\title{Get a speech operation}
\usage{
gl_speech_op(operation = .Last.value)
}
\arguments{
\item{operation}{A speech operation object from \link{gl_speech} when \code{asynch = TRUE}}
}
\value{
If the operation is still running, another operation object.  If done, the result as per \link{gl_speech}
}
\description{
For asynchronous calls of audio over 60 seconds, this returns the finished job
}
\examples{

\dontrun{

test_audio <- system.file("woman1_wb.wav", package = "googleLanguageR")

## make an asynchronous API request (mandatory for sound files over 60 seconds)
asynch <- gl_speech(test_audio, asynch = TRUE)

## Send to gl_speech_op() for status or finished result
gl_speech_op(asynch)

}

}
\seealso{
\link{gl_speech}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/text-to-speech.R
\name{gl_talk_shinyUI}
\alias{gl_talk_shinyUI}
\title{Speak in Shiny module (ui)}
\usage{
gl_talk_shinyUI(id)
}
\arguments{
\item{id}{The Shiny id}
}
\description{
Speak in Shiny module (ui)
}
\details{
Shiny Module for use with \link{gl_talk_shiny}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{is.NullOb}
\alias{is.NullOb}
\title{A helper function that tests whether an object is either NULL _or_
a list of NULLs}
\usage{
is.NullOb(x)
}
\description{
A helper function that tests whether an object is either NULL _or_
a list of NULLs
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/text-to-speech.R
\name{gl_talk}
\alias{gl_talk}
\title{Perform text to speech}
\usage{
gl_talk(
  input,
  output = "output.wav",
  languageCode = "en",
  gender = c("SSML_VOICE_GENDER_UNSPECIFIED", "MALE", "FEMALE", "NEUTRAL"),
  name = NULL,
  audioEncoding = c("LINEAR16", "MP3", "OGG_OPUS"),
  speakingRate = 1,
  pitch = 0,
  volumeGainDb = 0,
  sampleRateHertz = NULL,
  inputType = c("text", "ssml"),
  effectsProfileIds = NULL
)
}
\arguments{
\item{input}{The text to turn into speech}

\item{output}{Where to save the speech audio file}

\item{languageCode}{The language of the voice as a \code{BCP-47} language code}

\item{gender}{The gender of the voice, if available}

\item{name}{Name of the voice, see list via \link{gl_talk_languages} for supported voices.  Set to \code{NULL} to make the service choose a voice based on \code{languageCode} and \code{gender}.}

\item{audioEncoding}{Format of the requested audio stream}

\item{speakingRate}{Speaking rate/speed between \code{0.25} and \code{4.0}}

\item{pitch}{Speaking pitch between \code{-20.0} and \code{20.0} in semitones.}

\item{volumeGainDb}{Volumne gain in dB}

\item{sampleRateHertz}{Sample rate for returned audio}

\item{inputType}{Choose between \code{text} (the default) or SSML markup. The \code{input} text must be SSML markup if you choose \code{ssml}}

\item{effectsProfileIds}{Optional. An identifier which selects 'audio effects' profiles that are applied on (post synthesized) text to speech. Effects are applied on top of each other in the order they are given}
}
\value{
The file output name you supplied as \code{output}
}
\description{
Synthesizes speech synchronously: receive results after all text input has been processed.
}
\details{
Requires the Cloud Text-To-Speech API to be activated for your Google Cloud project.

Supported voices are here \url{https://cloud.google.com/text-to-speech/docs/voices} and can be imported into R via \link{gl_talk_languages}

To play the audio in code via a browser see \link{gl_talk_player}

To use Speech Synthesis Markup Language (SSML) select \code{inputType=ssml} - more details on using this to insert pauses, sounds and breaks in your audio can be found here: \url{https://cloud.google.com/text-to-speech/docs/ssml}

To use audio profiles, supply a character vector of the available audio profiles listed here: \url{https://cloud.google.com/text-to-speech/docs/audio-profiles} - the audio profiles are applied in the order given.  For instance \code{effectsProfileIds="wearable-class-device"} will optimise output for smart watches, \code{effectsProfileIds=c("wearable-class-device","telephony-class-application")} will apply sound filters optimised for smart watches, then telephonic devices.
}
\examples{

\dontrun{
library(magrittr)
gl_talk("The rain in spain falls mainly in the plain",
        output = "output.wav")

gl_talk("Testing my new audio player") \%>\% gl_talk_player()

# using SSML
gl_talk('<speak>The <say-as interpret-as=\"characters\">SSML</say-as>
  standard <break time=\"1s\"/>is defined by the
  <sub alias=\"World Wide Web Consortium\">W3C</sub>.</speak>',
  inputType =  "ssml")

# using effects profiles
gl_talk("This sounds great on headphones",
        effectsProfileIds = "headphone-class-device")

}

}
\seealso{
\url{https://cloud.google.com/text-to-speech/docs/}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/googleLanguageR.R
\docType{package}
\name{googleLanguageR}
\alias{googleLanguageR}
\title{googleLanguageR}
\description{
This package contains functions for analysing language through the
  Google Cloud Machine Learning APIs
}
\details{
For examples and documentation see the vignettes and the website:

\url{http://code.markedmondson.me/googleLanguageR/}
}
\seealso{
\url{https://cloud.google.com/products/machine-learning/}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/speech-to-text.R
\name{gl_speech}
\alias{gl_speech}
\title{Call Google Speech API}
\usage{
gl_speech(
  audio_source,
  encoding = c("LINEAR16", "FLAC", "MULAW", "AMR", "AMR_WB", "OGG_OPUS",
    "SPEEX_WITH_HEADER_BYTE"),
  sampleRateHertz = NULL,
  languageCode = "en-US",
  maxAlternatives = 1L,
  profanityFilter = FALSE,
  speechContexts = NULL,
  asynch = FALSE,
  customConfig = NULL
)
}
\arguments{
\item{audio_source}{File location of audio data, or Google Cloud Storage URI}

\item{encoding}{Encoding of audio data sent}

\item{sampleRateHertz}{Sample rate in Hertz of audio data. Valid values \code{8000-48000}. Optimal and default if left \code{NULL} is \code{16000}}

\item{languageCode}{Language of the supplied audio as a \code{BCP-47} language tag}

\item{maxAlternatives}{Maximum number of recognition hypotheses to be returned. \code{0-30}}

\item{profanityFilter}{If \code{TRUE} will attempt to filter out profanities}

\item{speechContexts}{An optional character vector of context to assist the speech recognition}

\item{asynch}{If your \code{audio_source} is greater than 60 seconds, set this to TRUE to return an asynchronous call}

\item{customConfig}{[optional] A \code{RecognitionConfig} object that will be converted from a list to JSON via \code{\link[jsonlite]{toJSON}} - see \href{https://cloud.google.com/speech-to-text/docs/reference/rest/v1p1beta1/RecognitionConfig}{RecognitionConfig documentation}. The \code{languageCode} will be taken from this functions arguments if not present since it is required.}
}
\value{
A list of two tibbles:  \code{$transcript}, a tibble of the \code{transcript} with a \code{confidence}; \code{$timings}, a tibble that contains \code{startTime}, \code{endTime} per \code{word}.  If maxAlternatives is greater than 1, then the transcript will return near-duplicate rows with other interpretations of the text.
 If \code{asynch} is TRUE, then an operation you will need to pass to \link{gl_speech_op} to get the finished result.
}
\description{
Turn audio into text
}
\details{
Google Cloud Speech API enables developers to convert audio to text by applying powerful
neural network models in an easy to use API.
The API recognizes over 80 languages and variants, to support your global user base.
You can transcribe the text of users dictating to an application’s microphone,
enable command-and-control through voice, or transcribe audio files, among many other use cases.
Recognize audio uploaded in the request, and integrate with your audio storage on Google Cloud Storage,
by using the same technology Google uses to power its own products.
}
\section{AudioEncoding}{


Audio encoding of the data sent in the audio message. All encodings support only 1 channel (mono) audio.
Only FLAC and WAV include a header that describes the bytes of audio that follow the header.
The other encodings are raw audio bytes with no header.
For best results, the audio source should be captured and transmitted using a
lossless encoding (FLAC or LINEAR16).
Recognition accuracy may be reduced if lossy codecs, which include the other codecs listed in this section,
are used to capture or transmit the audio, particularly if background noise is present.

Read more on audio encodings here \url{https://cloud.google.com/speech/docs/encoding}
}

\section{WordInfo}{



\code{startTime} - Time offset relative to the beginning of the audio, and corresponding to the start of the spoken word.

\code{endTime} - Time offset relative to the beginning of the audio, and corresponding to the end of the spoken word.

\code{word} - The word corresponding to this set of information.
}

\examples{

\dontrun{

test_audio <- system.file("woman1_wb.wav", package = "googleLanguageR")
result <- gl_speech(test_audio)

result$transcript
result$timings

result2 <- gl_speech(test_audio, maxAlternatives = 2L)
result2$transcript

result_brit <- gl_speech(test_audio, languageCode = "en-GB")


## make an asynchronous API request (mandatory for sound files over 60 seconds)
asynch <- gl_speech(test_audio, asynch = TRUE)

## Send to gl_speech_op() for status or finished result
gl_speech_op(asynch)

## Upload to GCS bucket for long files > 60 seconds
test_gcs <- "gs://mark-edmondson-public-files/googleLanguageR/a-dream-mono.wav"
gcs <- gl_speech(test_gcs, sampleRateHertz = 44100L, asynch = TRUE)
gl_speech_op(gcs)

## Use a custom configuration
my_config <- list(encoding = "LINEAR16",
                  diarizationConfig = list(
                    enableSpeakerDiarization = TRUE,
                    minSpeakerCount = 2,
                    maxSpeakCount = 3
                    ))

# languageCode is required, so will be added if not in your custom config
gl_speech(my_audio, languageCode = "en-US", customConfig = my_config)

}



}
\seealso{
\url{https://cloud.google.com/speech/reference/rest/v1/speech/recognize}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/text-to-speech.R
\name{gl_talk_languages}
\alias{gl_talk_languages}
\title{Get a list of voices available for text to speech}
\usage{
gl_talk_languages(languageCode = NULL)
}
\arguments{
\item{languageCode}{A \code{BCP-47} language tag.  If specified, will only return voices that can be used to synthesize this languageCode}
}
\description{
Returns a list of voices supported for synthesis.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{rmNullObs}
\alias{rmNullObs}
\title{Recursively step down into list, removing all such objects}
\usage{
rmNullObs(x)
}
\description{
Recursively step down into list, removing all such objects
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/translate.R
\name{gl_translate}
\alias{gl_translate}
\title{Translate the language of text within a request}
\usage{
gl_translate(
  t_string,
  target = "en",
  format = c("text", "html"),
  source = "",
  model = c("nmt", "base")
)
}
\arguments{
\item{t_string}{A character vector of text to detect language for}

\item{target}{The target language}

\item{format}{Whether the text is plain or HTML}

\item{source}{Specify the language to translate from. Will detect it if left default}

\item{model}{What translation model to use}
}
\value{
A tibble of \code{translatedText} and \code{detectedSourceLanguage}
  and \code{text} of length equal to the vector of text you passed in.
}
\description{
Translate character vectors via the Google Translate API
}
\details{
You can translate a vector of strings, although if too many for one call then it will be
  broken up into one API call per element.
  This is the same cost as charging is per character translated, but will take longer.

If translating HTML set the \code{format = "html"}.
Consider removing anything not needed to be translated first,
  such as JavaScript and CSS scripts. See example on how to do this with \code{rvest}

The API limits in three ways: characters per day, characters per 100 seconds,
  and API requests per 100 seconds.
All can be set in the API manager
  \url{https://console.developers.google.com/apis/api/translate.googleapis.com/quotas}
}
\examples{

\dontrun{

text <- "to administer medicine to animals is frequently a very difficult matter,
  and yet sometimes it's necessary to do so"

gl_translate(text, target = "ja")

# translate webpages using rvest to process beforehand
library(rvest)
library(googleLanguageR)

# translate webpages

# dr.dk article
my_url <- "http://bit.ly/2yhrmrH"

## in this case the content to translate is in css selector '.wcms-article-content'
read_html(my_url) \%>\%
  html_node(css = ".wcms-article-content") \%>\%
  html_text \%>\%
  gl_translate(format = "html")

}

}
\seealso{
\url{https://cloud.google.com/translate/docs/reference/translate}

Other translations: 
\code{\link{gl_translate_detect}()},
\code{\link{gl_translate_languages}()}
}
\concept{translations}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auth.R
\name{gl_auth}
\alias{gl_auth}
\alias{gl_auto_auth}
\title{Authenticate with Google language API services}
\usage{
gl_auth(json_file)

gl_auto_auth(...)
}
\arguments{
\item{json_file}{Authentication json file you have downloaded from your Google Project}

\item{...}{additional argument to
pass to \code{\link{gar_attach_auto_auth}}.}
}
\description{
Authenticate with Google language API services
}
\details{
The best way to authenticate is to use an environment argument pointing at your authentication file.

Set the file location of your download Google Project JSON file in a \code{GL_AUTH} argument

Then, when you load the library you should auto-authenticate

However, you can authenticate directly using this function pointing at your JSON auth file.
}
\examples{

\dontrun{
library(googleLanguageR)
gl_auth("location_of_json_file.json")
}

\dontrun{
library(googleLanguageR)
gl_auto_auth()
gl_auto_auth(environment_var = "GAR_AUTH_FILE")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/translate.R
\name{gl_translate_detect}
\alias{gl_translate_detect}
\title{Detect the language of text within a request}
\usage{
gl_translate_detect(string)
}
\arguments{
\item{string}{A character vector of text to detect language for}
}
\value{
A tibble of the detected languages with columns \code{confidence}, \code{isReliable}, \code{language}, and \code{text} of length equal to the vector of text you passed in.
}
\description{
Detect the language of text within a request
}
\details{
Consider using \code{library(cld2)} and \code{cld2::detect_language} instead offline,
since that is free and local without needing a paid API call.

\link{gl_translate} also returns a detection of the language,
so you could also wish to do it in one step via that function.
}
\examples{

\dontrun{

gl_translate_detect("katten sidder på måtten")
# Detecting language: 39 characters - katten sidder på måtten...
# confidence isReliable language                    text
# 1   0.536223      FALSE       da katten sidder på måtten


}

}
\seealso{
\url{https://cloud.google.com/translate/docs/reference/detect}

Other translations: 
\code{\link{gl_translate_languages}()},
\code{\link{gl_translate}()}
}
\concept{translations}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/text-to-speech.R
\name{gl_talk_player}
\alias{gl_talk_player}
\title{Play audio in a browser}
\usage{
gl_talk_player(audio = "output.wav", html = "player.html")
}
\arguments{
\item{audio}{The file location of the audio file.  Must be supported by HTML5}

\item{html}{The html file location that will be created host the audio}
}
\description{
This uses HTML5 audio tags to play audio in your browser
}
\details{
A platform neutral way to play audio is not easy, so this uses your browser to play it instead.
}
\examples{

\dontrun{

gl_talk("Testing my new audio player") \%>\% gl_talk_player()

}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/natural-language.R
\name{gl_nlp}
\alias{gl_nlp}
\title{Perform Natural Language Analysis}
\usage{
gl_nlp(
  string,
  nlp_type = c("annotateText", "analyzeEntities", "analyzeSentiment", "analyzeSyntax",
    "analyzeEntitySentiment", "classifyText"),
  type = c("PLAIN_TEXT", "HTML"),
  language = c("en", "zh", "zh-Hant", "fr", "de", "it", "ja", "ko", "pt", "es"),
  encodingType = c("UTF8", "UTF16", "UTF32", "NONE")
)
}
\arguments{
\item{string}{A vector of text to detect language for, or Google Cloud Storage URI(s)}

\item{nlp_type}{The type of Natural Language Analysis to perform.  The default \code{annotateText} will perform all features in one call.}

\item{type}{Whether input text is plain text or a HTML page}

\item{language}{Language of source, must be supported by API.}

\item{encodingType}{Text encoding that the caller uses to process the output}
}
\value{
A list of the following objects, if those fields are asked for via \code{nlp_type}:

\itemize{
 \item{sentences - }{\href{https://cloud.google.com/natural-language/docs/reference/rest/v1/Sentence}{Sentences in the input document}}
 \item{tokens - }{\href{https://cloud.google.com/natural-language/docs/reference/rest/v1/Token}{Tokens, along with their syntactic information, in the input document}}
 \item{entities - }{\href{https://cloud.google.com/natural-language/docs/reference/rest/v1/Entity}{Entities, along with their semantic information, in the input document}}
 \item{documentSentiment - }{\href{https://cloud.google.com/natural-language/docs/reference/rest/v1/Sentiment}{The overall sentiment for the document}}
 \item{classifyText -}{\href{https://cloud.google.com/natural-language/docs/classifying-text}{Classification of the document}}
 \item{language - }{The language of the text, which will be the same as the language specified in the request or, if not specified, the automatically-detected language}
 \item{text - }{The original text passed into the API. \code{NA} if not passed due to being zero-length etc. }
}
}
\description{
Analyse text entities, sentiment, syntax and categorisation using the Google Natural Language API
}
\details{
\code{string} can be a character vector, or a location of a file content on Google cloud Storage.
  This URI must be of the form \code{gs://bucket_name/object_name}

Encoding type can usually be left at default \code{UTF8}.
  \href{https://cloud.google.com/natural-language/docs/reference/rest/v1/EncodingType}{Read more here}

The current language support is available \href{https://cloud.google.com/natural-language/docs/languages}{here}
}
\examples{

\dontrun{

text <- "to administer medicince to animals is frequently a very difficult matter,
  and yet sometimes it's necessary to do so"
nlp <- gl_nlp(text)

nlp$sentences

nlp$tokens

nlp$entities

nlp$documentSentiment

## vectorised input
texts <- c("The cat sat one the mat", "oh no it didn't you fool")
nlp_results <- gl_nlp(texts)



}

}
\seealso{
\url{https://cloud.google.com/natural-language/docs/reference/rest/v1/documents}
}
