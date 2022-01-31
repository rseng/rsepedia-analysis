
<!-- README.md is generated from README.Rmd. Please edit that file -->

# reviewer

[![Project Status: Abandoned – Initial development has started, but there has not yet been a stable, usable release; the project has been abandoned and the author(s) do not intend on continuing development.](https://www.repostatus.org/badges/latest/abandoned.svg)](https://www.repostatus.org/#abandoned)


Improving the track changes and reviewing experience in R markdown.
`reviewer` provides two main functions:

  - an RStudio addin that adds the required JavaScript code to an
    rmarkdown document, so that when rendered to HTML it can be
    annotated using the Hypothes.is service
  - the capability to compare two versions of an rmarkdown document and
    display their differences in a nicely-formatted manner.

## Annotating web pages

### Important note

In order to use the annotation functionality it is needed to sign-up at
[the Hypothes.is website](https://hypothes.is/signup)

-----

## Differences between rmarkdown files

The `diff_rmd` function can be used to produce a nicely-formatted
document showing the differences between two rmarkdown files. This
function can be used to compare two files, or a file with previous
versions of itself (within a git repository).

See the [package
vignette](https://ropenscilabs.github.io/reviewer/articles/reviewer.html)
for a demonstration.

## Related packages

  - [trackmd](https://github.com/ropenscilabs/trackmd) is similar to
    `reviewer`, but:
    
      - is an RStudio-specific addin, whereas `reviewer` can be used
        outside of the RStudio environment (e.g. with your preferred
        text editor)
      - shows changes only in the *rendered* rmarkdown file (i.e. once
        it has been converted to its HTML document format). `reviewer`
        can show changes in either the raw rmarkdown document or its
        rendered output.

  - [latexdiffr](https://github.com/hughjonesd/latexdiffr) similarly
    shows differences in the *rendered* document, but uses the
    `latexdiff` utility to do so (you need `latexdiff` installed on your
    system to use it). It can also be used outside of RStudio.

  - [diffobj](https://github.com/brodieG/diffobj) provides a colourized
    depiction of the differences between arbitrary R objects. This could
    be used to compare two rmarkdown documents by e.g. reading their
    contents into character vectors and applying the `diffChr` function.

  - [rmdrive](https://github.com/ekothe/rmdrive) allows easy
    round-tripping of an rmarkdown document to Google Drive, where it
    can be edited by non-R-using collaborators, and back again. The
    edited changes could then be viewed using `reviewer`.

  - [markdrive](https://github.com/MilesMcBain/markdrive) is similar to
    `rmdrive`, but pushes the *rendered* rmarkdown document to and from
    Google.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```
# reviewer
Improving the track changes and reviewing experience in R markdown. `reviewer` provides two main functions:

- an RStudio addin that adds the required JavaScript code to an rmarkdown document, so that when rendered to HTML it can be annotated using the Hypothes.is service
- the capability to compare two versions of an rmarkdown document and display their differences in a nicely-formatted manner.

## Installation

You can install the development version of `reviewer` from [GitHub](https://github.com/ropenscilabs/reviewer) with:

```{r eval = FALSE}
remotes::install_github("ropenscilabs/reviewer")
```

## Annotating web pages

### Important note
In order to use the annotation functionality it is needed to sign-up at [the Hypothes.is website](https://hypothes.is/signup)

----

## Differences between rmarkdown files

The `diff_rmd` function can be used to produce a nicely-formatted document showing the differences between two rmarkdown files. This function can be used to compare two files, or a file with previous versions of itself (within a git repository).

See the [package vignette](https://ropenscilabs.github.io/reviewer/articles/reviewer.html) for a demonstration.

## Related packages

- [trackmd](https://github.com/ropenscilabs/trackmd) is similar to `reviewer`, but:

  - is an RStudio-specific addin, whereas `reviewer` can be used outside of the RStudio environment (e.g. with your preferred text editor)
  - shows changes only in the *rendered* rmarkdown file (i.e. once it has been converted to its HTML document format). `reviewer` can show changes in either the raw rmarkdown document or its rendered output.

- [latexdiffr](https://github.com/hughjonesd/latexdiffr) similarly shows differences in the *rendered* document, but uses the `latexdiff` utility to do so (you need `latexdiff` installed on your system to use it). It can also be used outside of RStudio.

- [diffobj](https://github.com/brodieG/diffobj) provides a colourized depiction of the differences between arbitrary R objects. This could be used to compare two rmarkdown documents by e.g. reading their contents into character vectors and applying the `diffChr` function.

- [rmdrive](https://github.com/ekothe/rmdrive) allows easy round-tripping of an rmarkdown document to Google Drive, where it can be edited by non-R-using collaborators, and back again. The edited changes could then be viewed using `reviewer`.

- [markdrive](https://github.com/MilesMcBain/markdrive) is similar to `rmdrive`, but pushes the *rendered* rmarkdown document to and from Google.
---
title: "Git Cheat Sheet"
author: "Amy Stringer"
date: "06/01/2018"
output: pdf_document
---

<!-- this file has been edited just for demonstration purposes within the reviewer package -->

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE)
```

Git is a type of version control software that stores the history of changes made to files in a particular repository. This document contains a brief rundown of the main commands used within the terminal to run git from your personal computer. Towards the end there will be details on how to use git to collaborate with other people using GitHub.

# Some lingo 

The **working directory** is the folder on your local machine that you are making changes in. 

The **staging area** is the middle ground between the working directory and the git repository. You *stage* changes to a file if you plan to add them to the repo.

The **repository** is where all files already committed go. You can take files from the repository and edit them, but the respository will only hold the unedited version until your changes are staged and then commited. 

**Untracked files** are files that exist in your working directory that have not yet been added to the repository and hence any changes made to this will not be tracked and you will not have access to previous versions. 

# Using git on your local computer only 

All commands listed below can be run through terminal provided git is installed on your local computer. If git is not installed, you can download it from [the git website](https://git-scm.com/downloads). If you have a windows computer, you should make sure to select that gitbash is included in your download, this will act as your terminal. 

## git Functions 

`git init` initialises a git repository. That is, if you run `git init` it will add a file called `.git` to your directory which contains all the back of house information, and sets up your repository. You will not need to access this file, but if you want to check that your initialisation worked, type `ls -a` into terminal and check that it is there. 

`git status` shows you the status of your repository. If you run this immediately after initiallising it will tell you that you have not yet made any commits to you repository, that is, your repository is empty. While in this state a repository cannot be cloned (git reccomends adding a README file to each newly initialised repo to describe the project, though it is not necessary if you already have files to commit to your repo.) `git status` will also display any untracked files that currently exist in your working directory. 

`git clone` clones an already existing repository, generally through a link from github or a file path to another folder on your local computer. This copies everything within that repository including the commit history so you can see any changes made to the repo by other people. 

`git log` prints a list of the previous commits to the repository. The latest comit is the top entry, and each commit should have a unique commit ID and a commit message detailing the change made 

`git diff` allows you to look at differences between different versions of a file. Though it can be used in different ways to investigate different things 

- `git diff commit1 commit2` looks at the differences between two separate commits in the history 
- `git diff --staged` looks at the differences between the a file in the staging area and the last known commit of that file
- `git diff` with no arguments looks at differences between the working directory and the staging area 

`git checkout` allows you to check out a previous commit, which means if any files in the folder are now opened,   

`git branch add branchname` creates a new branch called branchname. 

`git branch` with no arguments shows a list of the known branches, with the currently checked out branch marked with an asterisk 

`git merge branch1 branch2` merges two branches into each other for the most updated version of the code. For example, the master branch on a game may have continued receiving updates while another collaborator worked on another branch `addcolor`. Once the collaborator has the colour code working, they will want to merge this back with the master branch without erasing any of the changes made to the master since the `addcolor` branch was checked out. This is what git merge does. There was be conflicts, and when this happens, you will need to correct them manually. Git will show you where the conflicts are. This is a tricky part of git, and likely won't be needed much, but contact Amy is there are any issues. 
 
`git add filename` will add a file to the staging area once you have made edits. This is typically done once you intend to commit changes. 

`git commit` commits all files that are currently in the staging area to the repository. You must provide a message explaining what is contained in that particular commit. That is, what changed. These messages should be simple and descriptive, for example 'change distance calculation' and should be made after each new addition to the code in order for us to effectively use version control. Be careful here, this is equivalent to saving your work, only manually. Saving too little could results in an error being introduced that is quite hard to find in amongst a large number of changes, but committing too much could make it hard to pinpoint which commit introduced an error. 

Here is a random additional sentence.

# Using Git for Collaboration 

Using github allows you to set up a remote repository online that you can push data to or pull data from. If you clone a repository from github, that means that the repository you cloned is automatically set up as a remote from your local device. Otherwise, in order to access a new repo from github, you can set the remote manually and then pull from the repo, and this will have the same effect as cloning 

## Some functions 

`git remote` with no arguments lists any remotes you currently have by name 

`git remote add https://link` will add a remote to your local device where that remote repository can be found at http://link. This means you can save to an online repo that others can access once you have added files or made changes to files on you local computer. Typically the computer you are working from will automatically create the branch `master`, and automatically set your remote to be called `origin`.

`git pull origin master` means pull all data from origin (your remote) and store it in master (you master branch on you local computer). This downloads the most recent copy of all files within your repository, including the commits made by other collaborators.  

`git push origin master` means push to origin from the master branch. This has the opposite effect of the previous command, and this should be used each time you make a new commit, so that anyone else using the repo can see your changes right away. 






---
title: "Git Cheat Sheet"
author: "Amy Stringer"
date: "06/01/2018"
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

Git is a type of version control software that stores the history of changes made to files in a particular repository. Contained within is a brief rundown of the main commands used within the terminal to run git from your personal computer. Towards the end there will be details on how to use git to collaborate with others on files by using GitHub.

# Some lingo 

The **working directory** is the folder on your local machine that you are making changes in. 

The **staging area** is the middle ground between the working directory and the git repository. You *stage* changes to a file if you plan to add them to the repo.

The **repository** is where all files already committed go. You can take files from the repository and edit them, but the respository will only hold the unedited version until your changes are staged and then commited. 

**Untracked files** are files that exist in your working directory that have not yet been added to the repository and hence any changes made to this will not be tracked and you will not have access to previous versions. 

# Using git on your local computer only 

All commands listed below can be run through terminal provided git is installed on your local computer. If git is not installed, you can download it from [the git website](https://git-scm.com/downloads). If you have a windows computer, you should make sure to select that gitbash is included in your download, this will act as your terminal. 

## git Functions 

`git init` initialises a git repository. That is, if you run `git init` it will add a file called `.git` to your directory which contains all the back of house information, and sets up your repository. You will not need to access this file, but if you want to check that your initialisation worked, type `ls -a` into terminal and check that it is there. 

`git status` shows you the status of your repository. If you run this immediately after initiallising it will tell you that you have not yet made any commits to you repository, that is, your repository is empty. While in this state a repository cannot be cloned (git reccomends adding a README file to each newly initialised repo to describe the project, though it is not necessary if you already have files to commit to your repo.) `git status` will also display any untracked files that currently exist in your working directory. 

`git clone` clones an already existing repository, generally through a link from github or a file path to another folder on your local computer. This copies everything within that repository including the commit history so you can see any changes made to the repo by other people. 

`git log` prints a list of the previous commits to the repository. The latest comit is the top entry, and each commit should have a unique commit ID and a commit message detailing the change made 

`git diff` allows you to look at differences between different versions of a file. Though it can be used in different ways to investigate different things 

- `git diff commit1 commit2` looks at the differences between two separate commits in the history 
- `git diff --staged` looks at the differences between the a file in the staging area and the last known commit of that file
- `git diff` with no arguments looks at differences between the working directory and the staging area 

`git checkout` allows you to check out a previous commit, which means if any files in the folder are now opened,   

`git branch add branchname` creates a new branch called branchname. 

`git branch` with no arguments shows a list of the known branches, with the currently checked out branch marked with an asterisk 

`git merge branch1 branch2` merges two branches into each other for the most updated version of the code. E.g. Master branch on a game may have continued receiveing updates while another collaborator worked on another branch `addcolor`. Once the collaborator has the colour code working, they will want to merge this back with the master branch without erasing any of the changes made to the master since the `addcolor` branch was checked out. This is what git merge does. There was be conflicts, and when this happens, you will need to correct them manually. Git will show you where the conflicts are. This is a tricky part of git, and likely won't be needed much, but contact Amy is there are any issues. 
 
`git add filename` will add a file to the staging area once you have made edits. This is typically done once you intend to commit changes. 

`git commit` commits all files currently in the staging area to the repository. You must provide a message explaining what is contained in that particular commit. That is, what changed. These messages should be simple and descriptive, for example 'change distance calculation' and should be made after each new addition to the code in order for us to effectively use version control. Be careful here, this is equivalent to saving your work, only manually. Saving too little could results in an error being introduced that is quite hard to find in amongst a large number of changes, but committing too much could make it hard to pinpoint which commit introduced an error. 


# Using Git for Collaboration 

Using github allows you to set up a remote repository online that you can push data to or pull data from. If you clone a repository from github, that means that the repository you cloned is automatically set up as a remote from your local device. Otherwise, in order to access a new repo from github, you can set the remote manually and then pull from the repo, and this will have the same effect as cloning 

## Some functions 

`git remote` with no arguments lists any remotes you currently have by name 

`git remote add https://link` will add a remote to your local device where that remote repository can be found at http://link. This means you can save to an online repo that others can access once you have added files or made changes to files on you local computer. Typically the computer you are working from will automatically create the branch `master`, and automatically set your remote to be called `origin`.

`git pull origin master` means pull all data from origin (your remote) and store it in master (you master branch on you local computer). This downloads the most recent copy of all files within your repository, including the commits made by other collaborators.  

`git push origin master` means push to origin from the master branch. This has the opposite effect of the previous command, and this should be used each time you make a new commit, so that anyone else using the repo can see your changes right away. 






---
title: "Reviewer"
author: "Amy Stringer, Ben Raymond, Mark Dulhunty, Laura de Jong"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{reviewer}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Differences between rmarkdown files

The `diff_rmd` function can be used to produce a nicely-formatted document showing the differences between two rmarkdown files. For the purposes of demonstration, we'll use the example files bundled with the package. We'll compare `modified_file` to its earlier version `reference_file`:

```{r exfile1}
modified_file <- system.file("extdata/CheatSheet-modified.Rmd", package = "reviewer")
reference_file <- system.file("extdata/CheatSheet.Rmd", package = "reviewer")
```

```{r exfile2, echo = FALSE}
## for demo purposes, just take the first few lines of these files
tempfile1 <- tempfile()
writeLines(head(readLines(modified_file), 14), con = tempfile1)
modified_file <- tempfile1
tempfile2 <- tempfile()
writeLines(head(readLines(reference_file), 12), con = tempfile2)
reference_file <- tempfile2
```

Compare:

```{r dodiff}
library(reviewer)
result <- diff_rmd(modified_file, reference_file)
```

And our output:

```{r dodiff-hidden, echo = FALSE}
## some minor trickery to make the rendered output appear properly in the rendered file
htmltools::includeHTML(result$raw)
```

We can also compare the current version of a document to a previous version in stored in a git repository. (These examples are not run here).

If a `reference_file` argument is not provided, by default the `modified_file` will be compared to the most recent copy in the git repo:

```{r diff-ex-2, eval = FALSE}
result <- diff_rmd(modified_file)
```

Or we can compare it to how it appeared in the git repository after a particular commit (here, the commit with reference 750ab4):

```{r diff-ex-3, eval = FALSE}
result <- diff_rmd(reference_file, "750ab4")
```

## Annotating web pages

Anyone can view annotations, but to fully use the annotation functionality and
add your own annotations an account is needed at [the Hypothes.is
website](https://hypothes.is/signup). Annotations can be enabled by using the
RStudio addin `Insert html annotation js snippet`.  Annotations have even been
enabled for this document!  It works by adding the following javascript snippet
to the bottom of your RMarkdown document:

`<script src="https://hypothes.is/embed.js" async></script>`

<script src="https://hypothes.is/embed.js" async></script>
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enable_html_annotation.R
\name{enable_html_annotation}
\alias{enable_html_annotation}
\title{Enable html annotations using the Hypothes.is web annontation client
(https://web.hypothes.is/)}
\usage{
enable_html_annotation()
}
\value{
Hypothes.is web client javascript text snippet
}
\description{
Enable html annotations using the Hypothes.is web annontation client
(https://web.hypothes.is/)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reviewer-package.R
\docType{package}
\name{reviewer-package}
\alias{reviewer}
\alias{reviewer-package}
\title{reviewer: Improving the Track Changes and Reviewing Experience in R Markdown}
\description{
Provides functionality to compare two versions of an rmarkdown document and display their differences in a nicely-formatted manner, along with an RStudio addin that adds the required JavaScript code to an rmarkdown document, so that when rendered to HTML it can be annotated using the Hypothes.is service.
}
\author{
\strong{Maintainer}: Amy Stringer \email{amy.stringer@me.com}

Authors:
\itemize{
  \item Ben Raymond
  \item Mark Dulhunty
  \item Laura de Jong
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff_rmd_addin.R
\name{diff_rmd_addin}
\alias{diff_rmd_addin}
\title{Rstudio addin to display the raw (unrendered) differences between two
rmarkdown files}
\usage{
diff_rmd_addin()
}
\value{
Displays viewable html of the diff in the RStudio Viewer pane. If
  file is identical to previous version a message is provided. This may occur
  if changes since the last commit haven't been saved.
}
\description{
This addin is a small wrapper for \code{diff_rmd()}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmd_diffs.R
\name{diff_rmd}
\alias{diff_rmd}
\alias{diff_rmd_css}
\title{Render the differences between two rmarkdown files}
\usage{
diff_rmd(current_file, reference_file = "HEAD", show = "raw",
  output_format = "html_document", keep_intermediate = FALSE,
  quiet = TRUE, css = diff_rmd_css())

diff_rmd_css()
}
\arguments{
\item{current_file}{string: path to file after changes}

\item{reference_file}{string: path to file before changes}

\item{show}{string: \code{"raw"} (show differences in the raw rmarkdown file) or \code{"rendered"} (show differences in the rendered output). Note that the \code{"rendered"} output is unlikely to work well if a change has been made inside a code block in the rmarkdown document. The primary use for this package is anticipated to be the \code{"raw"} output}

\item{output_format}{string: format of the output file (currently only \code{"html_document"})}

\item{keep_intermediate}{logical: keep the intermediate rmarkdown file?}

\item{quiet}{logical: if \code{TRUE}, suppress pandoc output (Only applicable if \code{show="rendered"})}

\item{css}{character vector: css specification to apply to changed sections. Defaults to \code{diff_rmd_css()}; specify \code{NULL} to not include a \code{<style>} section in the output}
}
\value{
A list containing one or more elements \code{rendered} (the path to the rendered diff file, if \code{show="rendered"}), \code{intermediate} (the path to the intermediate file, if \code{keep_intermediate = TRUE}), and \code{raw} (if \code{show="raw"})

The path to the rendered file showing the differences
}
\description{
Render the differences between two rmarkdown files
}
\examples{
\dontrun{
  result <- diff_rmd(my_current_file, my_reference_file)
  browseURL(result)
}

}
