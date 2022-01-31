<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build
Status](https://travis-ci.com/ropenscilabs/qcoder.svg?branch=master)](https://travis-ci.com/ropenscilabs/qcoder)

QCoder
======

Lightweight package to conduct qualitative coding.

<img src="images/hex/imgHex.png" width="150" />

tl;dr
-----

To test with some sample data:

    install.packages("remotes")
    remotes::install_github("ropenscilabs/qcoder")
    library(qcoder)
    create_qcoder_project("my_qcoder_project", sample = TRUE)
    import_project_data(project = "my_qcoder_project")
    qcode()

Click "Select project folder" and "my_qcoder_project."  
There are two ways to add codes. To use an existing code, highlight the
text to be coded, select the code, click "Add selected code" and then
"Save changes." Text to be assigned a new (or existing) code should be
surrounded by (QCODER) (/QCODER) tags. The closing tag is followed
immediately by the code enclosed in curly brackets and prefixed with a
\# for example {\#samplecode}

Installation
------------

To install the latest development version, run

    install.packages("remotes")
    remotes::install_github("ropenscilabs/qcoder")
    library(qcoder)

Please note that this is not a release-ready version and should be
considered experimental and subject to changes. Still, we encourage you
to install and send us feedback on our issue tracker.

Motivation
----------

The motivation stems from the need for a free, open source option for
analyzing textual qualitative data. Textual qualitative data refers to
text from interview transcripts, observation notes, memos, jottings and
primary source/archival documents. A detailed discussion of the
motivation and other software can be found in our [motivation
document](https://github.com/ropenscilabs/qcoder/blob/master/motivation.Rmd).

Using QCoder
------------

QCoder is designed to be easy to use and to require minimal knowledge of
computer systems and code. Like all software, including other
applications for QDA there will be a learning period, but as we develop
Qcoder our goal will be to keep the interface simple and steadily
improve it. Currently we have a very minimal prototype.

Once you have installed QCoder, load it with the library command.

    library(qcoder)

This readme file is going to use sample data to illustrate basic QCoder
functionality. We will be using the simplest approach which is to use
the QCoder defaults for file names and folders. If you follow those same
patterns and conventions with your data you can use QCoder in the same
way. A full vignette will explain how to use non standard names and file
locations.

To begin we will create a QCoder project with sample data. (To create an
empty project leave out the sample option.)

    create_qcoder_project("my_qcoder_project", sample = TRUE)

This will create one main folder and four subfolders. Unless you
specified otherwise it will be in your current working directory (you
can find this with the `getwd()` command at the console). If you have a
specific location where you want to put the folder change your working
directory.

These will hold the documents to be coded, information about the codes,
unit information and the r data frames that will be the core of the
analysis. For this example the folder and file structures for the sample
data will look similar to this.

![](images/folderstructure.png)

### Documents

In our example we've already placed our documents into the "documents"
folder. At this point we only have tested support for txt files. If you
have documents in other formats you can use "Save As" to convert to txt.
If you have doc, docx, html, pdf, rtf or some other formats these can be
processed if you install the `textreadr` package. For many users this
will simply require

    install.packages("textreadr")

However for other users, particularly those on linux systems, additional
steps are required. Please follow the [installation instructions for
pdftools](https://github.com/ropensci/pdftools).

### Codes

QCoder has the option to import a list of predefined codes from a CSV
file (if you have this in a spreadsheet you can "Save As" csv). This
file should have exactly 3 columns with headings:

-   code\_id (A unique number for each code)
-   code (One word description, can use underscores or hyphens)
-   code.description (Longer description of the code, must be enclosed
    in quotation marks.)

To use project defaults, this file should be called *codes.csv*. Here
are the contents of the sample data csv file that comes with QCoder.

     code_id,code,code.description
    1,"harassment","define or describes harassing behavior"
    2,"person_talk","naming a specific person to talk to if there are violations"
    3,"gender","mentions gender"
    4,"gender_id","mentions gender identity or expression"
    5,"consequences","Detailing what happens if someone violates the code of conduct"
    6,"license","The license for this code of conduct"

You are not restricted to using the listed codes in the csv file, but
this file allows you to produce a detailed codebook including
descriptions. (Creating a user interface for adding new codes is high
priority item on the project road map.)

### Units

Units represent the unit of analysis data are about. Often this is
individual people, but it may also be organizations, events or
locations. Units may be associated with multiple documents. In the
sample data a minimum units file is used, but additional columns can be
used to assign attribute data.

The default file name is units.csv; if stored in a spreadsheet this can
be created by using "Save As" csv.

(Treatment of units is a work in progress and subject to change.)

    Filename,unit_id,Name
    1,"rOpenSci"
    2,"LIBD Rstats Club"
    3,"Carpentries"
    4,"Rladies

A second file (and data frame once imported) connects units to
documents.  
Our framework allows each unit to be associated with multiple documents
and each document with multiple units. (Note that the sample data is
designed to allow you to add more unit-document links and hence does not
link each unit to a document.)

    doc_path,unit_id
    CoC_Example1_mod_MU.txt,1
    CoC_Example1_MU.txt,2
    CoC_Example3_MU.txt,3
    CoC_Example4_MU.txt,4

### Importing the data

To import this data into Qcode user the `import_project_data()`
function. 

    import_project_data(project = "my_qcoder_project")

Now the data_frames folder will contain the imported files.
![](images/data_frames_folder.png)

Now it's time to start coding.

Coding uses a "Shiny App" to provide a user interface to the data. To
launch the app use the function `qcode()`.

    qcode()

Which will launch this application. If your current working directory is not the 
location of your project, use the `use_wd = FALSE` option.
However, on Windows this will not work unless you have set a HOME or R_USER.

Coding
------

Once you have selected your project there will be a drop down menu on
the "Add codes to text" tab to allow you to pick a specific document to
code. This will pull a document into the editor.

![](images/coding_step1.png)

Select your project folder.

![](images/coding_step2_folder_select.png)

Once you have a project, use the drop down menu to select a particular
document to code. This will open in an editor. When done coding
(instructions below), click Save changes.

Select your project folder.

![](images/coding_step3_document_select.png)

Switching to the "Codes" tab a list of codes from the codes file is
displayed.

![](images/codestab.png)

Our sample data already has some coding done, and the code-text data is
displayed on the "Coded data" tab.

![](images/codeddata.png)

### Coding the data

To add codes to the documents uses a tagging system. Text to be assigned
a code should be surrounded by (QCODER) (/QCODER) tags. The closing tag
is followed immediately by the code enclosed in curly brackets and
prefixed with a \# for example {\#samplecode}

(QCODE)This is the text that is being assigned a
code.(/QCODE){\#instructions}

One pair of {} can contain multiple codes, each with at \# and separated
by commas.

Alternatively, to use an existing code, highlight the text to be coded,
select the code or codes, click "Add selected code."

When you have finished coding a document press the "Save changes"
button.

### Cautions and known issues

Each time you save, Qcoder makes a backup copy of your documents data
frame. This is for safety and reproducability. This can end up with a
lot of files if you save often. You may want to periodically delete some
backups to save storage space.  An important goal is to move to using 
git for this purpose.

Currently when you create a new code while coding, this code will be
displayed on the Coded data tab, but not on the Codes or Summary tabs.
You must go to the first tab of the the qcode application to update those displays. This is
a high priority development item.

### Road map

QCoder can be used right now for coding. However, we are not yet ready
for release.

Our immediate goal is to create a somewhat more advanced minimum viable
product. Please see the issue tracker for a list of short-term and longer-term goals. 
These goals include interoperability with other QDA packages.

The most important thing is to have more people try qcoder and give us feedback! We
do not want to release and then discover that our testing has missed problems that
are obvious to our intended user base.

Contributors
------------

-   [Elin Waring](https://github.com/elinw)
-   [Dan Sholler](https://github.com/dsholler)
-   [Jenny Draper](https://github.com/learithe)
-   [Beth Duckles](https://github.com/bduckles)
---
title: "Contributing"
output: md_document
---

# Contributing

Contributions to `qcoder` whether in the form of bug fixes, issue reports, new
code or documentation improvement are welcome. Please use the github issue
tracker. For any pull request please link to or open a corresponding issue in the
issue tracker. Please ensure that you have notifications turned on and respond to
questions, comments or needed changes promptly.

In making contributions and suggestions please keep in mind the general approach
and road map outlined below as well as the motivation for qcoder. The intent
of this project is not to replace a full featured proprietary QDA application.
In our view, the advantage of open source is that, instead, qcoder can offer
the distinctive features of qualitative coding while providing ways to 
interact with other R packages that implement specific analytical techniques.
In deciding whether to implement specific features we will always first
ask whether it would be possible to interact with (or possibly help to build)
a separate R package. Our main focus is on coding and providing ways to
interact with other packages, not analysis.

Our approach is currently highly focused on ease of use for 
people who have not used R previously and who have possibly 
not used QDA software at all.  Therefore
we have made the decision to first focus on a simple, tightly structured 
project model that does not require the installation of Java or of a data base.

Once we have a satisfactory interface there are a number of extensions that
we will make that will allow more flexibility and power and hence support 
larger and more complex projects and users with more experience.  These will
include more flexiblity on the location of files (including the ability to
use text accessible on the Web) and storage in a data base. Therefore, all
code contributions should be made with this future in mind.  Functions should
work from the command line, which also ensures that they can be efficiently
tested. We welcome interested people, both users and developers, to discuss
these plans in our issue tracker and help make them happen.

If you would like to make a code or documentation contribution, please open
an issue in the tracker that you link to a pull request. 

---
title: "Motivation"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Motivation

The motivation stems from the need for a free, open source option for analyzing
textual qualitative data. Textual qualitative data refers to text from 
interview transcripts, observation notes, memos, jottings and primary 
source/archival documents. 

Qualitative data analysis (QDA) processes, particularly those developed by 
Corbin and Strauss  (2014), Miles, Huberman, and Saldana (2013), and 
Glaser and Strauss (2017), can be thought of as layering interpretation onto 
the text. The researcher starts with open coding, meaning that she is free to 
tag snippets of text with whatever descriptions she deems appropriate. For 
instance, if the researcher is coding observation notes and senses that 
conversations between two individuals will be relevant to the research 
questions, she might tag the instances in which the two individuals speak with 
the code "conversation." In the next round of coding, she might classify what 
the participants discuss with a finer tag, like "conversation_package" if they 
were talking about creating packages in R.  In another round, she might get 
even more specific with codes such as "conversation_package_nomoney" if a 
participant discussed not having money to create a package in R. Later rounds
involve conflating codes that might mean the same thing, relating codes to one
another (often by documenting their meanings in a similar way to software 
documentation), and eliminating codes that no longer make sense.

Sometimes the QDA process involves looking for instances that demonstrate some 
concept, mechanism, or theory from the academic literature on the subject. This 
approach is common toward the end of the process in fields such as organization 
studies, management, and other disciplines and is prominent as an approach from
the beginning in fields with experimental or clinical roots (e.g., psychology).

Using software for QDA allows researchers to nest codes, then begin to see the 
number of instances in which a particular code has been applied. Perhaps more 
importantly, software allows easy visualization and analysis of where codes 
co-occur (i.e., where multiple codes have been applied to the same snippets of 
text), and other linking activities that help researchers identify and specify 
themes in the data. Some software packages allow the researcher to visualize 
codes and the relationships between them in new and innovative ways; however, a 
cursory review of the literature suggests that qualitative researchers often 
use very basic features; some even use software such as Endnote on Microsoft 
Word or even analog systems such as paper or sticky notes.

QDA is an iterative process. Researchers will often change, lump together, split,
or re-organize codes as they analyze their data. Depending on the coding approach, researchers might also create a codebook or research notes for each code that 
defines the code and specifies instances in which the code should be applied. 

## Current QDA software

To date, researchers who conduct QDA largely rely upon proprietary software 
such as those listed below. Free and open source options are limited, not 
easy to integrate with R and   
(credit to [bduckles](https://github.com/bduckles) for the descriptions): 

* [NVivo – QSR International](http://www.qsrinternational.com/product)
Desktop based software with both Windows and Mac support. Functionality 
includes video and images as data and auto-coding processes. Tech support 
and training is available and a relatively large user base. 

* [Atlas.TI](http://atlasti.com/) Desktop based software for Mac and Windows,
also has mobile apps for android and iOS. Similar functionality as NVivo 
including using video and visuals as data. Has training and support. 

* [MAXQDA](http://www.maxqda.com/) Desktop based software that works for 
both Mac and Windows. Training and several tiers of licenses available.

* [QDA Miner – Provalis Research]
(https://provalisresearch.com/products/qualitative-data-analysis-software/) 
A Windows application to analyze qualitative text, can be used with some 
visuals as well. They also have a “lite” package which is free. 

* [HyperResearch – Researchware]
(http://www.researchware.com/products/hyperresearch.html) Desktop based 
software that works on both Mac and Windows. 

* [Dedoose](http://www.dedoose.com/) Web based software, usable on any 
platform. Data stored in the web instead of on your device. Includes text, 
photos, audio, video and spreadsheet information. 

* [Annotations](http://www.annotationsapp.com/) Mac only app that allows 
you to highlight, keyword and create notes for text documents. An inexpensive
way to do basic qualitative data analysis. Last updated in 2014. 

* [RQDA](http://rqda.r-forge.r-project.org/) A QDA package for R that is free 
and open source. Bugs in the program make it challenging to use. Last updated 
in 2012. 

* [Coding Analysis Toolkit – CAT](http://cat.texifter.com/) A free, web based, 
open source service that was created by University of Pittsburgh. Can be used 
on any platform. Not supported. 

* [Weft QDA](http://www.pressure.to/qda/) A free and open source Windows 
program that has no support. The website says that bugs in the program can 
result in the loss of data. 

## Limitations of Existing Software

Each extant software package has its limitations. The foremost limitation 
is cost, which can prohibit students and underfunded qualitative researchers 
from conducting analyses systematically, and efficiently. Furthermore, the 
mature software packages (e.g., Atlas.TI, NVIVO) offer features that exceed 
the needs of many users and, as a result, suffer speed issues (particularly 
for those researchers who may not benefit from advanced hardware). The sharing
process for proprietary QDA outputs is equally unwieldy, relying on 
non-intuitive bundling and unbundling processes, steep learning curves, and 
non-transferable skill development. 

Open source languages such as R offer the opportunity to involve qualitative
researchers in open source software development. Greater involvement of 
qualitative researchers serves to expand the scope of R users and could 
create inroads to connect qualitative and quantitative R packages. For 
instance, better integration of qualitative research packages into R would 
make it possible for existing text analysis programs to work alongside 
qualitative coding.  

## User Needs

We began this project at rOpenSci’s 
[runconf18](https://github.com/ropensci/unconf18) based on
[bduckles](https://github.com/bduckles) issue. We began by discussing
the common challenges we and our peers face in conducting QDA and considering
what we might learn from existing quantitative and text analysis open source
packages. We identified a number of unique challenges qualitative researchers
encounter when analyzing data. 

* *GUI for broader adoption.* Quantitative researchers are often familiar with
conducting analyses from the command line. Qualitative researchers, on the other
hand, often rely upon Graphical User Interfaces (GUI) to perform analysis tasks.
A proposed, in-development solution is integrating Shiny apps to permit a 
personalized GUI for QCoder. 

* *Collaboration.* Several aspects of collaboration in qualitative research 
present challenges for development and use of QDA software. The first challenge
relates to version control and tracking changes. Whereas contributors to the
open source software community, and perhaps the software development community more
generally, have developed processes for sharing code, data, and other work products,
qualitative researchers often rely on email, DropBox/box, and other sharing 
technologies that lack appropriate version control features. One proposed solution
is to implement a system of record-keeping that maintains all files akin to
Atlas.TI’s bundling and unbundling process, but that permits simultaneous work 
and comparison/merging of changes. We might also consider urging users to learn 
git and facilitate skill development by offering tutorials. 
  The second collaboration challenge relates to the need to establish 
inter-rater or inter-coder reliability. Existing methods of assessing 
inter-rater reliability for QDA include Cronbach’s Alpha and Krippendorff's
Alpha, among others. Such evaluations, though, do not account for 
complexities such as decisions not to code a segment of text. In other words, 
by deciding not to apply a code to a snippet of text might be considered an 
equivalent level of agreement as a decision to apply the same code to the same 
snippet of text. QCoder aims to make inter-rater reliability assessments easier 
and more effective than traditional approaches by allowing direct comparisons 
of differences in the codebooks it generates. 

* *Privacy concerns.* Qualitative data often includes personally-identifying 
information. Quantitative approaches to deidentification are not always 
applicable to qualitative datasets because text-based data includes rich 
descriptions that may elevate the risk of re-identification. Furthermore, 
qualitative data tends to be unstructured in ways that make systematic 
approaches to de-identification difficult to implement. QDA software, then, 
needs to account for these complexities.

* *Reproducibility.* Qualitative research communities, particularly those 
communities employing methods to probe human and social behavior (such as 
ethnography), fiercely debate whether or not qualitative research is meant 
to be reproducible. In other words, some researchers argue that human and 
social behavior is contextually bound by moments in time, physical settings, 
and the generally fluid nature of the phenomena under study. Accepting this 
dispute as valid, we are striving to develop a package that can facilitate 
some level of reproduction of data analysis. 
  Our approach centers on the production of codebooks. Qualitative researchers 
create codebooks while coding their data. These codebooks match the code with 
a description of how and when the code is used. The creation of the codebook 
is an iterative process that is concurrent with the coding process. Codes are 
routinely combined, split and shifted as the researcher does their analysis.
Codebooks are a way for the researcher to indicate to other members of their 
research team how to apply codes to the data. Historically, these code books 
have not been standard across different QDA packages. 
  There could be a benefit to having a standard codebook format which could be
compared across projects and even as a project iterates. Using Git or other 
version control as a codebook is created could prevent mistakes and be a 
eaching and learning tool for qualitative research. Additionally, most 
codebooks do not contain sensitive or private data and they could be shared 
and made publicly available. One could envision the capacity for researchers 
to share their codebooks so that research that follows could draw on existing 
codes for future research. A key characteristic of QCoder is that it exposes 
QDA processes to critique in ways that proprietary software excludes. For 
example, it is unreasonable to expect that researchers interested in 
reviewing or reproducing the analysis done for a manuscript pay for a license
to perform the task at hand. QCoder can be easily installed and potentially 
maintained so that critical mistakes in analysis might be caught and corrected
in review processes.

* *Money.* Some of the disciplines in which qualitative work occurs do not 
enjoy the same level of funding as, say, natural sciences and engineering.
Expensive software packages, then, may not be available to researchers. The 
case is equally dire for students and researchers in under-resourced communities 
and countries. An R-based QDA package enables broader participation in 
qualitative research and, by extension, broader perspectives on the nature of
human and social behavior.  

## Conceptualizing a Minimum Viable Product: A vignette

Questa the graduate student is working on one chapter of her dissertation research 
where she is trying to understand how different codes of conduct in the 
software community welcome underrepresented communities. She wants to understand 
language use and to trace how these codes of conduct are created from differing 
communities. She also is planning to do interviews with people who have created 
and used these codes of conduct to understand how their adoption influences 
inclusion.  
 
None of the grants that Questa sent out have come back with any money to do 
this research. She’s discouraged but passionate about her work and she knows 
she needs to finish this chapter of her dissertation so she can graduate 
already. She’s done some work with a nonprofit and they support her work. 
She’s pretty sure she has a postdoc with them when she finishes her degree. 
While they should be able to hire her and do support her research, it’s 
unlikely they’d have enough funds for her to get a license for one of the 
QDA packages. She has been able to try out some of the larger QDA packages 
and they work ok with her student license, but her computer is on its last 
legs and the program keeps crashing her computer. 

She’s planning to start up her interviews and overall she thinks she’ll have 
maybe around 50-75 text files that include 1) codes of conduct 2) interview 
transcripts 3) field research at a conference. 

She’s not sure how many codes she’d need, but she’s guessing that the first 
round of the coding would be a lot of different codes - maybe 150 - 200 
codes to start. Then she’d likely change them and boil it down into groups 
of codes and categories of those codes. So ideally the program would make it
easy for her to edit, change and move around the codes. She also needs to 
write up a simple codebook as she does the coding. 


## References

Corbin, J., & Strauss, A. L. (2014). Basics of qualitative research. 
Thousand Oaks, CA: Sage.

Glaser, B. G., & Strauss, A. L. (2017). Discovery of grounded theory: 
Strategies for qualitative research. London, UK: Routledge.

Miles, M. B., & Huberman, A. M. (1994). Qualitative data analysis: 
An expanded sourcebook. Thousand Oaks, CA: Sage.

Miles, M. B., Huberman, A. M., & Saldana, J. (2013). Qualitative data 
analysis. Thousand Oaks, CA: Sage.

Saldaña, J. (2015). The coding manual for qualitative researchers. 
Thousand Oaks, CA: Sage.
---
title: "Contributing"
output: md_document
---

# Contributing

Contributions to `qcoder` whether in the form of bug fixes, issue reports, new
code or documentation improvement are welcome. Please use the github issue
tracker. For any pull request please link to or open a corresponding issue in the
issue tracker. Please ensure that you have notifications turned on and respond to
questions, comments or needed changes promptly.

In making contributions and suggestions please keep in mind the general approach
and road map outlined below as well as the motivation for qcoder. The intent
of this project is not to replace a full featured proprietary QDA application.
In our view, the advantage of open source is that, instead, qcoder can offer
the distinctive features of qualitative coding while providing ways to 
interact with other R packages that implement specific analytical techniques.
In deciding whether to implement specific features we will always first
ask whether it would be possible to interact with (or possibly help to build)
a separate R package. Our main focus is on coding and providing ways to
interact with other packages, not analysis.

Our approach is currently highly focused on ease of use for 
people who have not used R previously and who have possibly 
not used QDA software at all.  Therefore
we have made the decision to first focus on a simple, tightly structured 
project model that does not require the installation of Java or of a data base.

Once we have a satisfactory interface there are a number of extensions that
we will make that will allow more flexibility and power and hence support 
larger and more complex projects and users with more experience.  These will
include more flexiblity on the location of files (including the ability to
use text accessible on the Web) and storage in a data base. Therefore, all
code contributions should be made with this future in mind.  Functions should
work from the command line, which also ensures that they can be efficiently
tested. We welcome interested people, both users and developers, to discuss
these plans in our issue tracker and help make them happen.

If you would like to make a code or documentation contribution, please open
an issue in the tracker that you link to a pull request. 

---
title: ""
output: html_document
---


 
```{r, echo=FALSE, message=FALSE}
# Load your libraries here

library('qcoder')
library('dplyr')
library('pander')
```

```{r, echo=FALSE, message=FALSE}
# Put the exact code that this report should be about inside of quotation marks.
# Example "gender"
code_to_report_on <- CODE

# This can be copied from the qcode() application. Include quotation marks.
project_path <- PATH_TO_PROJECT_HERE


```

```{r, echo=FALSE, message=FALSE}
# Do the code work
    project_name <- basename(project_path)
    docs_df_path <- paste0(project_path,  
                           "/data_frames/qcoder_documents_", 
                          project_name, ".rds")
    codes_df_path <- paste0(project_path,  
                            "/data_frames/qcoder_codes_", 
                            project_name, ".rds")


        text_df <- readRDS(docs_df_path)
       parsed <-  qcoder::parse_qcodes(text_df)
       parsed <- parsed %>% dplyr::filter(qcode == code_to_report_on) %>%
                        dplyr::select(doc, text)
       
```

## Project `r project_name`
## Report on `r code_to_report_on` 

```{r, echo=FALSE, message=FALSE}
pander::pander(parsed)
```



Type your text here.
---
title: "Project Summary"
output: html_document
---


 
```{r, echo=FALSE, message=FALSE}
# Load your libraries here

library('qcoder')
library('dplyr')
library('knitr')
```

```{r, echo=FALSE, message=FALSE}

# This can be copied from the qcode() application. Include quotation marks.
project_path <- PATH_TO_PROJECT_HERE

```

```{r, echo=FALSE, message=FALSE}
# Do the code work
    project_name <- basename(project_path)
    docs_df_path <- paste0(project_path,  
                           "/data_frames/qcoder_documents_", 
                          project_name, ".rds")
    codes_df_path <- paste0(project_path,  
                            "/data_frames/qcoder_codes_", 
                            project_name, ".rds")
    units_df_path <- paste0(project_path,  
                            "/data_frames/qcoder_units_", 
                            basename(project_path), ".rds")


      text_df <- readRDS(docs_df_path)
      parsed <-  qcoder::parse_qcodes(text_df)
      codes_df <- readRDS(codes_df_path)
      units_df <- readRDS(units_df_path)
       
```

## Project `r project_name`
## Summary

Number of documents: `r nrow(text_df)`  
Number of codes: `r nrow(codes_df)`  
Number of segments coded: `r nrow(parsed)`
Number of units:  `r nrow(units_df)`

Frequency of codes
```{r}
parsed %>% group_by(qcode) %>% summarise(n = n())
```




---
title: "Project Summary"
output: html_document
---


 
```{r, echo=FALSE, message=FALSE}
# Load your libraries here

library('qcoder')
library('knitr')
```

```{r, echo=FALSE, message=FALSE}

# This can be copied from the qcode() application. Include quotation marks.
project_path <- PATH_TO_PROJECT_HERE

```

```{r, echo=FALSE, message=FALSE}
# Do the code work
    project_name <- basename(project_path)
    codes_df_path <- paste0(project_path,  
                            "/data_frames/qcoder_codes_", 
                            project_name, ".rds")

    codes_df <- readRDS(codes_df_path)
```

## Project `r project_name`
## Codebook


```{r}
knitr::kable(codes_df)
```




---
title: "Motivation"
output:
  pdf_document: default
  html_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Motivation

The motivation stems from the need for a free, open source option for analyzing textual qualitative data. Textual qualitative data refers to text from interview transcripts, observation notes, memos, jottings and primary source/archival documents. 

Qualitative data analysis (QDA) processes, particularly those developed by Corbin and Strauss  (2014), Miles, Huberman, and Saldana (2013), and Glaser and Strauss (2017), can be thought of as layering interpretation onto the text. The researcher starts with open coding, meaning that she is free to tag snippets of text with whatever descriptions she deems appropriate. For instance, if the researcher is coding observation notes and senses that conversations between two individuals will be relevant to the research questions, she might tag the instances in which the two individuals speak with the code "conversation." In the next round of coding, she might classify what the participants discuss with a finer tag, like "conversation_package" if they were talking about creating packages in R.  In another round, she might get even more specific with codes such as "conversation_package_nomoney" if a participant discussed not having money to create a package in R. Later rounds involve conflating codes that might mean the same thing, relating codes to one another (often by documenting their meanings in a similar way to software documentation), and eliminating codes that no longer make sense.

Sometimes the QDA process involves looking for instances that demonstrate some concept, mechanism, or theory from the academic literature on the subject. This approach is common toward the end of the process in fields such as organization studies, management, and other disciplines and is prominent as an approach from the beginning in fields with experimental or clinical roots (e.g., psychology).

Using software for QDA allows researchers to nest codes, then begin to see the number of instances in which a particular code has been applied. Perhaps more importantly, software allows easy visualization and analysis of where codes co-occur (i.e., where multiple codes have been applied to the same snippets of text), and other linking activities that help researchers identify and specify themes in the data. Some software packages allow the researcher to visualize codes and the relationships between them in new and innovative ways; however, a cursory review of the literature suggests that qualitative researchers often use very basic features; some even use software such as Endnote on Microsoft Word or even analog systems such as paper or sticky notes.

QDA is an iterative process. Researchers will often change, lump together, split, or re-organize codes as they analyze their data. Depending on the coding approach, researchers might also create a codebook or research notes for each code that defines the code and specifies instances in which the code should be applied. 

## Current QDA software

To date, researchers who conduct QDA largely rely upon proprietary software such as those listed below. Free and open source options are limited, not easy to integrate with R and   (credit to [bduckles](https://github.com/bduckles) for the descriptions): 

* [NVivo - QSR International](http://www.qsrinternational.com/product)
Desktop based software with both Windows and Mac support. Functionality includes video and images as data and auto-coding processes. Tech support and training is available and a relatively large user base. 

* [Atlas.TI](http://atlasti.com/) Desktop based software for Mac and Windows, also has mobile apps for android and iOS. Similar functionality as NVivo including using video and visuals as data. Has training and support. 

* [MAXQDA](http://www.maxqda.com/) Desktop based software that works for both Mac and Windows. Training and several tiers of licenses available.

* [QDA Miner - Provalis Research](https://provalisresearch.com/products/qualitative-data-analysis-software/) A Windows application to analyze qualitative text, can be used with some visuals as well. They also have a "lite" package which is free. 

* [HyperResearch - Researchware](http://www.researchware.com/products/hyperresearch.html) Desktop based software that works on both Mac and Windows. 

* [Dedoose](http://www.dedoose.com/) Web based software, usable on any platform. Data stored in the web instead of on your device. Includes text, photos, audio, video and spreadsheet information. 

* [Annotations](http://www.annotationsapp.com/) Mac only app that allows you to highlight, keyword and create notes for text documents. An inexpensive way to do basic qualitative data analysis. Last updated in 2014. 

* [RQDA](http://rqda.r-forge.r-project.org/) A QDA package for R that is free and open source. Bugs in the program make it challenging to use. Last updated in 2012. 

* [Coding Analysis Toolkit - CAT](http://cat.texifter.com/) A free, web based, open source service that was created by University of Pittsburgh. Can be used on any platform. Not supported. 

* [Weft QDA](http://www.pressure.to/qda/) A free and open source Windows program that has no support. The website says that bugs in the program can result in the loss of data. 

## Limitations of Existing Software

Each extant software package has its limitations. The foremost limitation is cost, which can prohibit students and underfunded qualitative researchers from conducting analyses systematically, and efficiently. Furthermore, the mature software packages (e.g., Atlas.TI, NVIVO) offer features that exceed the needs of many users and, as a result, suffer speed issues (particularly for those researchers who may not benefit from advanced hardware). The sharing process for proprietary QDA outputs is equally unwieldy, relying on non-intuitive bundling and unbundling processes, steep learning curves, and non-transferable skill development. 

Open source languages such as R offer the opportunity to involve qualitative researchers in open source software development. Greater involvement of qualitative researchers serves to expand the scope of R users and could create inroads to connect qualitative and quantitative R packages. For instance, better integration of qualitative research packages into R would make it possible for existing text analysis programs to work alongside qualitative coding.  

## User Needs

We began this project at rOpenSci's [runconf18](https://github.com/ropensci/unconf18) based on [bduckles](https://github.com/bduckles) issue. We began by discussing the common challenges we and our peers face in conducting QDA and considering what we might learn from existing quantitative and text analysis open source packages. We identified a number of unique challenges qualitative researchers encounter when analyzing data. 

* *GUI for broader adoption.* Quantitative researchers are often familiar with conducting analyses from the command line. Qualitative researchers, on the other hand, often rely upon Graphical User Interfaces (GUI) to perform analysis tasks. A proposed, in-development solution is integrating Shiny apps to permit a personalized GUI for QCoder. 

* *Collaboration.* Several aspects of collaboration in qualitative research present challenges for development and use of QDA software. The first challenge relates to version control and tracking changes. Whereas contributors to the open source software community, and perhaps the software development community more generally, have developed processes for sharing code, data, and other work products, qualitative researchers often rely on email, DropBox/box, and other sharing technologies that lack appropriate version control features. One proposed solution is to implement a system of record-keeping that maintains all files akin to Atlas.TI's bundling and unbundling process, but that permits simultaneous work and comparison/merging of changes. We might also consider urging users to learn git and facilitate skill development by offering tutorials. 
  The second collaboration challenge relates to the need to establish inter-rater or inter-coder reliability. Existing methods of assessing inter-rater reliability for QDA include Cronbach's Alpha and Krippendorff's Alpha, among others. Such evaluations, though, do not account for complexities such as decisions not to code a segment of text. In other words, by deciding not to apply a code to a snippet of text might be considered an equivalent level of agreement as a decision to apply the same code to the same snippet of text. QCoder aims to make inter-rater reliability assessments easier and more effective than traditional approaches by allowing direct comparisons of differences in the codebooks it generates. 

* *Privacy concerns.* Qualitative data often includes personally-identifying information. Quantitative approaches to deidentification are not always applicable to qualitative datasets because text-based data includes rich descriptions that may elevate the risk of re-identification. Furthermore, qualitative data tends to be unstructured in ways that make systematic approaches to de-identification difficult to implement. QDA software, then, needs to account for these complexities.

* *Reproducibility.* Qualitative research communities, particularly those communities employing methods to probe human and social behavior (such as ethnography), fiercely debate whether or not qualitative research is meant to be reproducible. In other words, some researchers argue that human and social behavior is contextually bound by moments in time, physical settings, and the generally fluid nature of the phenomena under study. Accepting this dispute as valid, we are striving to develop a package that can facilitate some level of reproduction of data analysis. 
  Our approach centers on the production of codebooks. Qualitative researchers create codebooks while coding their data. These codebooks match the code with a description of how and when the code is used. The creation of the codebook is an iterative process that is concurrent with the coding process. Codes are routinely combined, split and shifted as the researcher does their analysis. Codebooks are a way for the researcher to indicate to other members of their research team how to apply codes to the data. Historically, these code books have not been standard across different QDA packages. 
  There could be a benefit to having a standard codebook format which could be compared across projects and even as a project iterates. Using Git or other version control as a codebook is created could prevent mistakes and be a teaching and learning tool for qualitative research. Additionally, most codebooks do not contain sensitive or private data and they could be shared and made publicly available. One could envision the capacity for researchers to share their codebooks so that research that follows could draw on existing codes for future research. A key characteristic of QCoder is that it exposes QDA processes to critique in ways that proprietary software excludes. For example, it is unreasonable to expect that researchers interested in reviewing or reproducing the analysis done for a manuscript pay for a license to perform the task at hand. QCoder can be easily installed and potentially maintained so that critical mistakes in analysis might be caught and corrected in review processes.

* *Money.* Some of the disciplines in which qualitative work occurs do not enjoy the same level of funding as, say, natural sciences and engineering. Expensive software packages, then, may not be available to researchers. The case is equally dire for students and researchers in under-resourced communities and countries. An R-based QDA package enables broader participation in qualitative research and, by extension, broader perspectives on the nature of human and social behavior.  

## Conceptualizing a Minimum Viable Product: A vignette

Questa the graduate student is working on one chapter of her dissertation research where she is trying to understand how different codes of conduct in the software community welcome underrepresented communities. She wants to understand language use and to trace how these codes of conduct are created from differing communities. She also is planning to do interviews with people who have created and used these codes of conduct to understand how their adoption influences inclusion.  
 
None of the grants that Questa sent out have come back with any money to do this research. She's discouraged but passionate about her work and she knows she needs to finish this chapter of her dissertation so she can graduate already. She's done some work with a nonprofit and they support her work. She's pretty sure she has a postdoc with them when she finishes her degree. While they should be able to hire her and do
support her research, it's unlikely they'd have enough funds for her to get a license for one of the QDA packages. She has been able to try out some of the larger QDA packages and they work ok with her student license, but her computer is on its last legs and the program keeps crashing her computer. 

She's planning to start up her interviews and overall she thinks she'll have maybe around 50-75 text files that include 1) codes of conduct 2) interview transcripts 3) field research at a conference. 

She's not sure how many codes she'd need, but she's guessing that the first round of the coding would be a lot of different codes -- maybe 150 -- 200 codes to start. Then she'd likely change them and boil it down into groups of codes and categories of those codes. So ideally the program would make it easy for her to edit, change and move around the codes. She also needs to write up a simple codebook as she does the coding. 


## References

Corbin, J., & Strauss, A. L. (2014). Basics of qualitative research. Thousand Oaks, CA: Sage.

Glaser, B. G., & Strauss, A. L. (2017). Discovery of grounded theory: Strategies for qualitative research. London, UK: Routledge.

Miles, M. B., & Huberman, A. M. (1994). Qualitative data analysis: An expanded sourcebook. Thousand Oaks, CA: Sage.

Miles, M. B., Huberman, A. M., & Saldana, J. (2013). Qualitative data analysis. Thousand Oaks, CA: Sage.

Saldana, J. (2015). The coding manual for qualitative researchers. Thousand Oaks, CA: Sage.
---
title: "Contributing"
output:
  pdf_document: default
  html_document: default
---

## Contributing

Contributions to `qcoder` whether in the form of bug fixes, issue reports, new
code or documentation improvement are welcome. Please use the [GitHub issue
tracker](https://github.com/ropenscilabs/qcoder/issues). For any pull request please link to or open a corresponding issue in the
issue tracker. Please ensure that you have notifications turned on and respond to
questions, comments or needed changes promptly.

In making contributions and suggestions please keep in mind the general approach
and road map outlined below as well as the motivation for qcoder. The intent
of this project is not to replace a full featured proprietary QDA application.
In our view, the advantage of open source is that, instead, qcoder can offer
the distinctive features of qualitative coding while providing ways to 
interact with other R packages that implement specific analytical techniques.
In deciding whether to implement specific features we will always first
ask whether it would be possible to interact with (or possibly help to build)
a separate R package. Our main focus is on coding and providing ways to
interact with other packages, not analysis.

Our approach is currently highly focused on ease of use for 
people who have not used R previously and who have possibly 
not used QDA software at all. Therefore
we have made the decision to first focus on a simple, tightly structured 
project model that does not require the installation of Java or of a data base.

Once we have a satisfactory interface there are a number of extensions that
we will make that will allow more flexibility and power and hence support 
larger and more complex projects and users with more experience. These will
include more flexibility on the location of files (including the ability to
use text accessible on the Web) and storage in a data base. Therefore, all
code contributions should be made with this future in mind. Functions should
work from the command line, which also ensures that they can be efficiently
tested. We welcome interested people, both users and developers, to discuss
these plans in our issue tracker and help make them happen.

If you would like to make a code or documentation contribution, please open
an issue in the tracker that you link to a pull request. 
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{create_qcoder_project}
\alias{create_qcoder_project}
\title{Create a standard set of folders for a QCoder project}
\usage{
create_qcoder_project(project_name, sample = FALSE)
}
\arguments{
\item{project_name}{A string project name to be located in the
current working directory or a path to a project folder.}

\item{sample}{Logical that indicates that the sample data should be copied to
the project.}
}
\description{
Create a standard set of folders for a QCoder project
}
\examples{
create_qcoder_project(project_name = "my_qcoder_project")
unlink("./my_qcoder_project", recursive=TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{validate_project}
\alias{validate_project}
\title{Check for a valid qcoder project}
\usage{
validate_project(path_to_test)
}
\arguments{
\item{path_to_test}{Path to possible project folder}
}
\value{
NULL for valid project, Error otherwise.
}
\description{
Check for a valid qcoder project
}
\examples{
create_qcoder_project(project_name = "_my_qcoder_project")
validate_project("_my_qcoder_project")
unlink("./_my_qcoder_project", recursive=TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_in_files.R
\name{add_new_documents}
\alias{add_new_documents}
\title{Add new documents
Adds new document or documents to an existing documents data frame.}
\usage{
add_new_documents(files, docs_df_path = "", file_path = "")
}
\arguments{
\item{files}{Vector of new files to be added}

\item{docs_df_path}{Path to existing data frame of text documents}

\item{file_path}{Full path to the data set of documents including
trailing slash}
}
\description{
Add new documents
Adds new document or documents to an existing documents data frame.
}
\examples{
create_qcoder_project(project_name = "my_qcoder_project", sample = TRUE)

unlink("./my_qcoder_project", recursive=TRUE)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_in_files.R
\name{create_empty_docs_file}
\alias{create_empty_docs_file}
\title{Create an empty documents data set}
\usage{
create_empty_docs_file(path)
}
\arguments{
\item{path}{Full path to data frame to be created.}
}
\description{
Used to create a codes data frame with no data but that can
have data added. File is placed in the data_frames folder.
}
\examples{
create_qcoder_project(project_name = "_my_qcoder_project")
path <- file.path(getwd(),
  "_my_qcoder_project/data_frames/qcoder_docs__my_qcoder_project")
create_empty_docs_file(path)
unlink("./_my_qcoder_project", recursive=TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{update_links}
\alias{update_links}
\title{Update document to unit links
Saves or updates the links between observation units and documents}
\usage{
update_links(
  checked = "",
  docs_df_path = "",
  this_doc_path = "",
  units_docs_path = ""
)
}
\arguments{
\item{checked}{vector of new or updated links}

\item{docs_df_path}{full path to document dataset}

\item{this_doc_path}{value of doc_path for the document}

\item{units_docs_path}{full path of the data frame of unit to docs links}
}
\description{
Update document to unit links
Saves or updates the links between observation units and documents
}
\examples{

unlink("./_my_qcoder_project", recursive=TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_in_files.R
\name{read_documents_data}
\alias{read_documents_data}
\title{Create a data frame of documents}
\usage{
read_documents_data(
  project_name,
  data_path = "documents/",
  df_path = "data_frames",
  data_frame_name = "qcoder_documents",
  project_path = ""
)
}
\arguments{
\item{project_name}{Name of the Qcoder project}

\item{data_path}{path to a folder contain text files to be analyzed.}

\item{df_path}{Full path to the docs data frame.}

\item{data_frame_name}{The name of the RDS file that the data frame
will be stored in.}

\item{project_path}{Full path to the project folder.}
}
\description{
Create a data frame of documents
}
\examples{
 \dontrun{
read_documents_data("_my_qcoder_project")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selection.R
\name{add_codes_to_selection}
\alias{add_codes_to_selection}
\title{Adds codes surrounding the selected text}
\usage{
add_codes_to_selection(selection, codes)
}
\arguments{
\item{selection}{The selection of text to be coded}

\item{codes}{The code or codes to be added to the document}
}
\description{
Adds codes surrounding the selected text
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{txt2html}
\alias{txt2html}
\title{Format text as HTML
Minimal conversion of a text to html}
\usage{
txt2html(text)
}
\arguments{
\item{text}{text to be converted}
}
\description{
Format text as HTML
Minimal conversion of a text to html
}
\examples{
txt2html("The quick brown (QCODE)fox(/QCODE){#animal} jumped over ")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{qcode_custom}
\alias{qcode_custom}
\title{This launches the coder custom Shiny app}
\usage{
qcode_custom()
}
\description{
This launches the coder custom Shiny app
}
\examples{
if (interactive()) {
  qcode_custom()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_qcodes.R
\name{add_discovered_code}
\alias{add_discovered_code}
\title{Update codes data frame
Add discovered codes to the codes data frame}
\usage{
add_discovered_code(
  codes_list = "",
  code_data_frame = NULL,
  codes_df_path = ""
)
}
\arguments{
\item{codes_list}{A list of codes (usually from a coded document)}

\item{code_data_frame}{Existing data frame of QCODE codes}

\item{codes_df_path}{The path where the updated code data frame should be saved}
}
\description{
Update codes data frame
Add discovered codes to the codes data frame
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{add_code}
\alias{add_code}
\title{Add code
Append a new unit record to the existing data frame}
\usage{
add_code(codes_df, new_code, new_code_desc, codes_df_path)
}
\arguments{
\item{codes_df}{Existing codes data frame}

\item{new_code}{text name of a new code (single name only)}

\item{new_code_desc}{text description of the code}

\item{codes_df_path}{full path to the codes data frame}
}
\description{
Add code
Append a new unit record to the existing data frame
}
\examples{
unlink("./_my_qcoder_project", recursive=TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_in_files.R
\name{read_unit_document_map_data}
\alias{read_unit_document_map_data}
\title{Create a data frame of unit to document links from csv file
Use this is you have a spreadsheet already created.}
\usage{
read_unit_document_map_data(
  project_name,
  data_path = "units/unit_document_map.csv",
  data_frame_name = "qcoder_unit_document_map",
  project_path = "",
  df_path = "data_frames"
)
}
\arguments{
\item{project_name}{Name of project}

\item{data_path}{path to a file containing  unit document map data in csv.}

\item{data_frame_name}{The name of the RDS file that the data frame
will be stored in.}

\item{project_path}{Full path to the project folder}

\item{df_path}{Full path to the documents data frame.}
}
\description{
Create a data frame of unit to document links from csv file
Use this is you have a spreadsheet already created.
}
\examples{
create_qcoder_project(project_name = "_my_qcoder_project", sample = TRUE)
project_name = "_my_qcoder_project"
read_unit_document_map_data( project_name = "_my_qcoder_project")
unlink("./_my_qcoder_project", recursive=TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_in_files.R
\name{create_empty_code_file}
\alias{create_empty_code_file}
\title{Create an empty codes data set}
\usage{
create_empty_code_file(path)
}
\arguments{
\item{path}{Full path to data frame to be created.}
}
\description{
Used to create a codes data frame with no data but that can
have data added. File is placed in the data_frames folder.
}
\examples{
create_qcoder_project(project_name = "_my_qcoder_project")
path <- file.path(getwd(),
  "_my_qcoder_project/data_frames/qcoder_codes__my_qcoder_project")
create_empty_docs_file(path)
unlink("./_my_qcoder_project", recursive=TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qcoder-package.R
\docType{package}
\name{qcoder-package}
\alias{qcoder-package}
\alias{qcoder}
\title{Code Qualitative Data
A light weight approach to qualitative coding and analysis}
\description{
Light weight coding
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{qcode}
\alias{qcode}
\title{This launches the coder Shiny app}
\usage{
qcode(use_wd = TRUE)
}
\arguments{
\item{use_wd}{Whether or not the current working directory when launching
qcoder should be used as the base from which the project file is selected.}
}
\description{
This launches the coder Shiny app
}
\examples{
if (interactive()) {
 qcode()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_in_files.R
\name{build_paths}
\alias{build_paths}
\title{Build the paths for file creation}
\usage{
build_paths(
  project_name,
  data_path = "",
  data_frame_name = "",
  df_path = "data_frames",
  project_path = ""
)
}
\arguments{
\item{project_name}{Name of the project}

\item{data_path}{Path segment to the data. The format for this may
depend on the function using the paths.}

\item{data_frame_name}{Name of the data frame that will contain the data}

\item{df_path}{path segment(s) to the created data frame file
from the project path}

\item{project_path}{Path to the project (not including project_name). This
will be set to getwd() if a value of "" is passed in.}
}
\value{
A named list of paths. "data_frame_path" is the path to the data frame
and "data" is the path to the data.
}
\description{
Builds the paths to the data to be imported and to the data
frame where the imported data is to be stored.
The project name is required, all other parameters may be set.
These each represent a segment of the path to a file.
If a project path is not set or set to "" it will be set to
the current working directory via getwd().
The Shiny qcode application assumes that the data frame folder will
be "data_frames" and that any new documents to be imported will
be in a folder called "documents".
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_in_files.R
\name{create_empty_units_file}
\alias{create_empty_units_file}
\title{Define an empty  units data frame}
\usage{
create_empty_units_file(path)
}
\arguments{
\item{path}{Full path to data frame to be created.}
}
\description{
Define an empty  units data frame
}
\examples{
create_qcoder_project(project_name = "_my_qcoder_project")
path <- file.path(getwd(),
  "_my_qcoder_project/data_frames/qcoder_units__my_qcoder_project")
create_empty_docs_file(path)
unlink("./_my_qcoder_project", recursive=TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_in_files.R
\name{read_code_data}
\alias{read_code_data}
\title{Create a file of codes from csv file
Use this if you have a spreadsheet of codes already created.}
\usage{
read_code_data(
  project_name,
  data_path = "codes/codes.csv",
  df_path = "data_frames",
  data_frame_name = "qcoder_codes",
  project_path = ""
)
}
\arguments{
\item{project_name}{Name of the project, which matches folder name}

\item{data_path}{Path to a file containing  code data in csv.}

\item{df_path}{Full path to the codes data frame.}

\item{data_frame_name}{The name of the RDS file that the data frame
will be stored in.}

\item{project_path}{Full path to the project folder}
}
\description{
Create a file of codes from csv file
Use this if you have a spreadsheet of codes already created.
}
\examples{
create_qcoder_project(project_name = "_my_qcoder_project", sample = TRUE)
read_code_data(project_name = "_my_qcoder_project")
unlink("./_my_qcoder_project", recursive=TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{add_unit}
\alias{add_unit}
\title{Add unit
Append a new unit record to the existing data frame}
\usage{
add_unit(units_df, new_unit, units_df_path)
}
\arguments{
\item{units_df}{Existing units data frame}

\item{new_unit}{text name of a new unit (single name only)}

\item{units_df_path}{full path to the units data frame}
}
\description{
Add unit
Append a new unit record to the existing data frame
}
\examples{
unlink("./_my_qcoder_project", recursive=TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_qcodes.R
\name{parse_splititem}
\alias{parse_splititem}
\title{Parse a single item within a document}
\usage{
parse_splititem(splititem, df, doc_id, dots)
}
\arguments{
\item{splititem}{String usually generated by parse_one_document()}

\item{df}{Data frame for storing parsed data}

\item{doc_id}{The doc_id for the document this string is part of.}

\item{dots}{List of additional options passed in.}
}
\description{
Parse a single item within a document
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_qcodes.R
\name{get_codes}
\alias{get_codes}
\title{Extract codes from text
Take coded text and extract the codes, assuming they are correctly formatted.}
\usage{
get_codes(doc_text)
}
\arguments{
\item{doc_text}{The text data for a single document}
}
\description{
Extract codes from text
Take coded text and extract the codes, assuming they are correctly formatted.
}
\examples{
unlink("./my_qcoder_project", recursive=TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qcoder-package.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{\%>\%}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{magrittr}{\code{\link[magrittr:pipe]{\%>\%}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_in_files.R
\name{read_unit_data}
\alias{read_unit_data}
\title{Create a data frame of units from csv file
Use this is you have a spreadsheet of units already created.}
\usage{
read_unit_data(
  data_path = "units/units.csv",
  data_frame_name = "qcoder_units",
  project_name,
  project_path = "",
  df_path = "data_frames"
)
}
\arguments{
\item{data_path}{path to a file containing  unit data in csv.}

\item{data_frame_name}{The name of the RDS file that the data frame
will be stored in.}

\item{project_name}{Name of project if available}

\item{project_path}{Full path to the project folder.}

\item{df_path}{Full path to the units data frame.}
}
\description{
Create a data frame of units from csv file
Use this is you have a spreadsheet of units already created.
}
\examples{
create_qcoder_project(project_name = "_my_qcoder_project", sample = TRUE)
read_unit_data(project_name = "_my_qcoder_project")
unlink("./_my_qcoder_project", recursive=TRUE)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_qcodes.R
\name{parse_qcodes}
\alias{parse_qcodes}
\title{Parse coded text}
\usage{
parse_qcodes(x, ...)
}
\arguments{
\item{x}{A data frame containing the text to be coded; requires columns
"doc_id" and "document_text"}

\item{...}{Other parameters optionally passed in}
}
\value{
If the data frame contains coded text in the \code{document_text}
column, output will be a data frame with three columns: "doc",
"qcode", and "text".\preformatted{    The \code{doc} is the \code{doc_id} from the input data frame.

    \code{qcode} is the code that the captured text was marked up with.

    \code{text} is the text that was captured.
}
}
\description{
Take a data frame of coded text documents and return a data frame of the
codes captured within.
}
\details{
This function takes a text document containing coded text of the form:
\preformatted{"stuff to ignore (QCODE) coded text we care about (/QCODE){#my_code}
more stuff to ignore"} and turns it into a data frame with one row per coded
item, of the form: \code{docid,qcode,text}

\code{parse_qcodes} assumes that it is being passed a data frame, the
\code{\link{parse_one_document}} function is called to do the heavy lifting
extracting the coded text from the \code{document_text} column.

Newline characters are replaced with an HTML \code{<br>} in the captured text.

If no valid qcodes are found, \code{parse_qcodes} returns an empty data frame
(no rows).
}
\examples{
parse_qcodes(my_documents)

# Data frames can be piped into this function
my_documents \%>\%
  parse_qcodes()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_qcodes.R
\name{error_check}
\alias{error_check}
\title{Check for coding errors}
\usage{
error_check(document)
}
\arguments{
\item{document}{A string to be scanned for errors.}
}
\value{
A \code{warning} message as a character string.
}
\description{
Checks the current document for coding errors.
}
\details{
This function takes a string (such as the contents of a document), and conducts some basic linting. It returns a warning if there aren't a matching number of \code{(QCODE)} tags, or if text has been marked to be captured but the capture is missing a tag (missing \code{{#my_tag}}).
}
\examples{
error_check("An (QCODE)unmatched set of (QCODE) gives (/QCODE){#tag} a warning.")
error_check("A (QCODE) qcode with a missing tag gives a warning.")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{read_data}
\alias{read_data}
\title{This launches the data-reader Shiny app}
\usage{
read_data()
}
\description{
This launches the data-reader Shiny app
}
\examples{
 \dontrun{
 read_data()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{validate_project_files}
\alias{validate_project_files}
\title{Check for required imported data frames.}
\usage{
validate_project_files(path_to_test)
}
\arguments{
\item{path_to_test}{Path to possible project folder}
}
\value{
NULL for valid project, Error otherwise.
}
\description{
Check for required imported data frames.
}
\examples{
create_qcoder_project(project_name = "_my_qcoder_project", sample = TRUE)
import_project_data("_my_qcoder_project")
validate_project_files("_my_qcoder_project")
unlink("./_my_qcoder_project", recursive=TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{do_update_document}
\alias{do_update_document}
\title{Update document
Updates the text field of the documents data frame, typically after
pressing Save button in the Shiny App.  May
also be used in the console.}
\usage{
do_update_document(updated, docs_df_path, this_doc_path)
}
\arguments{
\item{updated}{The updated text as a character string}

\item{docs_df_path}{Location of the documents rds file.}

\item{this_doc_path}{Name of record to be updated, as recorded in "doc_path"
field of data frame.}
}
\description{
Update document
Updates the text field of the documents data frame, typically after
pressing Save button in the Shiny App.  May
also be used in the console.
}
\examples{
unlink("./_my_qcoder_project", recursive=TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_qcodes.R
\name{parse_one_document}
\alias{parse_one_document}
\title{Parse one document}
\usage{
parse_one_document(doc, df, qcoder_documents, dots = NULL)
}
\arguments{
\item{doc}{A single document from qcoder_data}

\item{df}{The data frame that will contain the parsed data}

\item{qcoder_documents}{The full documents data frame}

\item{dots}{Other parameters that may be passed in.}
}
\description{
Parse one document
}
\examples{
unlink("./my_qcoder_project", recursive=TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_in_files.R
\name{import_project_data}
\alias{import_project_data}
\title{Read data into a project
Convenience method to read raw data from standard locations and using
standard names in a project folder structure.}
\usage{
import_project_data(project_name)
}
\arguments{
\item{project_name}{The project name. This should represent the folder
holding the project.}
}
\description{
Read data into a project
Convenience method to read raw data from standard locations and using
standard names in a project folder structure.
}
\examples{
create_qcoder_project(project_name = "_my_qcoder_project", sample = TRUE)
import_project_data("_my_qcoder_project")
unlink("./_my_qcoder_project", recursive=TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_in_files.R
\name{create_empty_unit_doc_file}
\alias{create_empty_unit_doc_file}
\title{Define an empty many to many unit to document map}
\usage{
create_empty_unit_doc_file(path)
}
\arguments{
\item{path}{Full path to data frame to be created.}
}
\description{
Define an empty many to many unit to document map
}
\examples{
create_qcoder_project(project_name = "_my_qcoder_project")
path <- file.path(getwd(),
  "_my_qcoder_project/data_frames/qcoder_units_document_map__my_qcoder_project")
create_empty_docs_file(path)
unlink("./_my_qcoder_project", recursive=TRUE)
}
