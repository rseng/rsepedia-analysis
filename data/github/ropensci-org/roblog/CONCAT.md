
<!-- README.md is generated from README.Rmd. Please edit that file -->

# roblog

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
<!-- badges: end -->

Are you preparing a blog post or tech note for rOpenSci, after getting
the go-ahead and a tentative publication date from our community manager
[Stefanie Butland](https://ropensci.org/authors/stefanie-butland/)? Stay
here\! The goal of roblog is to provide practical helpers for authors of
blog posts and tech notes on rOpenSci website.

Best paired with [rOpenSci blog
guidance](https://blogguide.ropensci.org/).

## Installation

``` r
remotes::install_github("ropenscilabs/roblog")
```
---
slug: "post-template"
title: Wonderful Title
package_version: 0.1.0
authors:
  - Author Name
date: 2019-06-04
categories: blog
topicid:
tags:
  - Software Peer Review
  - R
  - community
# delete the line below
# if you have no preferred image
# for Twitter cards
twitterImg: img/blog-images/2019-06-04-post-template/name-of-image.png
---

Save this file under /content/blog/YEAR-MONTH-DAY-slug.md in the local copy of your roweb2 fork.
{{< figure src="orange-mug-near-macbook-3219546.jpg" width="300" link="https://www.pexels.com/photo/orange-mug-near-macbook-3219546/" alt="Laptop keyboard with a tree leaf beside it" class="center" caption="Another type of leaf. [Engin Akyurt on Pexels](https://www.pexels.com/photo/orange-mug-near-macbook-3219546/)." >}}
s
{{< figure src="orange-mug-near-macbook-3219546.jpg" width="300" link="https://www.pexels.com/photo/orange-mug-near-macbook-3219546/" alt="Laptop keyboard with a tree leaf beside it" >}}
### Heading in sentence case

#### Another sentence as heading
---
slug: "post-template"
title: Wonderful title
package_version: 0.1.0
authors:
  - Author Name
date: 2019-06-04
categories: blog
topicid:
tags:
  - Software Peer Review
  - R
  - community
# delete the line below
# if you have no preferred image
# for Twitter cards
twitterImg: img/blog-images/2019-06-04-post-template/name-of-image.png
---

Save this file under /content/blog/YEAR-MONTH-DAY-slug.md in the local copy of your roweb2 fork.

[Cool blog](/blog/)

[Broken blog](https://masalmon.eu/40004)

[Broken blog again](https://masalmon.eu/400040)

Beware! If you want to generate this post from R Markdown, use the R Markdown template instead!

  Everywhere in this template (YAML, paths to images), you should change "post-template" to the slug of your post, and "2019-06-04" to the publication date.

Introduction including outline of the post.

### First awesome section

I like Hugo[^1]. Yes, that is how you add a footnote.

#### First awesome subsection of the first awesome section

Here's how to use a Hugo shortcode to add an image.

{{< figure src = "/img/blog-images/2019-06-04-post-template/name-of-image.png" width = "200">}}

{{< figure src = "/img/blog-images/2019-06-04-post-template/name-of-image.png" width = "200" alt = "too short">}}

<blockquote class="twitter-tweet"><p lang="en" dir="ltr">Finally... hello subtools 1.0! 🥳 Read, write and manipulate subtitles in <a href="https://twitter.com/hashtag/rstats?src=hash&amp;ref_src=twsrc%5Etfw">#rstats</a>. Substantially re-written to integrate with tidytext by <a href="https://twitter.com/juliasilge?ref_src=twsrc%5Etfw">@juliasilge</a> and <a href="https://twitter.com/drob?ref_src=twsrc%5Etfw">@drob</a> <a href="https://t.co/QmCWGk9NOX">https://t.co/QmCWGk9NOX</a> cc <a href="https://twitter.com/ma_salmon?ref_src=twsrc%5Etfw">@ma_salmon</a> <a href="https://t.co/7576oktL7k">pic.twitter.com/7576oktL7k</a></p>&mdash; Francois Keck (@FrancoisKeck) <a href="https://twitter.com/FrancoisKeck/status/1200040510540386304?ref_src=twsrc%5Etfw">November 28, 2019</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script> 


<blockquote class="twitter-tweet"><p lang="en" dir="ltr">When I try to become acquainted with a new (to me) <a href="https://twitter.com/hashtag/rstats?src=hash&amp;ref_src=twsrc%5Etfw">#rstats</a> package, I prefer to read ___________</p>&mdash; Jonathan Carroll (@carroll_jono) <a href="https://twitter.com/carroll_jono/status/969442252610191361?ref_src=twsrc%5Etfw">March 2, 2018</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script> 

<!--html_preserve--> 
{{% figure src = "/img/blog-images/2019-06-04-post-template/name-of-image.png" width = "200" % alt = "an actually good description!"%}}
<!--/html_preserve-->

Ropensci was very nice, cool ropensci.

#### Second awesome subsection

Here's how to use a Hugo shortcode to embed a tweet. We recommend the use of [Hugo shortcodes](https://gohugo.io/content-management/shortcodes/) to include tweets, Vimeo or Youtube videos, gists, etc.

{{< tweet 1138216112808529920 >}}

### Conclusion

Have fun writing your blog post!

 [rOpenSci blog](https://ropensci.org/blog)

 Here's how to add the footnote text for your reference above.

[^1]: Hugo! https://gohugo.io/
---
slug: "post-template"
title: Wonderful title that is not in title case
package_version: 0.1.0
authors:
  - Author Name
date: 2019-06-04
categories: blog
topicid:
tags:
  - Software Peer Review
  - R
  - community
# delete the line below
# if you have no preferred image
# for Twitter cards
twitterImg: img/blog-images/2019-06-04-post-template/name-of-image.png
---

Save this file under /content/blog/YEAR-MONTH-DAY-slug.md in the local copy of your roweb2 fork.
---
slug: "post-template"
title: Wonderful Title
package_version: 0.1.0
authors:
  - Author Name
date: 2019-06-04
categories: blog
topicid:
tags:
  - Software Peer Review
  - R
  - community
# delete the line below
# if you have no preferred image
# for Twitter cards
twitterImg: img/blog-images/2019-06-04-post-template/name-of-image.png
---

Save this file under /content/blog/YEAR-MONTH-DAY-slug.md in the local copy of your roweb2 fork.

### Heading in Title Case

#### Another sentence as heading

#### Another Title as Heading
---
slug: "post-template"
title: Wonderful Title
package_version: 0.1.0
authors:
  - Author Name
date: 2019-06-04
categories: blog
topicid:
tags:
  - Software Peer Review
  - R
  - community
# delete the line below
# if you have no preferred image
# for Twitter cards
twitterImg: img/blog-images/2019-06-04-post-template/name-of-image.png
---

Save this file under /content/blog/YEAR-MONTH-DAY-slug.md in the local copy of your roweb2 fork.

### Heading in sentence case

#### Another sentence as heading

[Good link](https://twitter.com), [Click here](https://twitter.com), [here](https://twitter.com).
---
slug: "post-template"
title: Wonderful Title
package_version: 0.1.0
authors:
  - Author Name
date: 2019-06-04
categories: blog
topicid:
tags:
  - Software Peer Review
  - R
  - community
# delete the line below
# if you have no preferred image
# for Twitter cards
twitterImg: img/blog-images/2019-06-04-post-template/name-of-image.png
---

Save this file under /content/blog/YEAR-MONTH-DAY-slug.md in the local copy of your roweb2 fork.

### Heading in sentence case

#### Another sentence as heading

![altbabla](imagepath/on/laptop.png)
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>",
fig.path = "man/figures/README-",
out.width = "100%"
)
```

# roblog

<!-- badges: start -->
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
<!-- badges: end -->

Are you preparing a blog post or tech note for rOpenSci, after getting the go-ahead and a tentative publication date from our community manager [Stefanie Butland](https://ropensci.org/authors/stefanie-butland/)? Stay here! The goal of roblog is to provide practical helpers for authors of blog posts and tech notes on rOpenSci website. 

Best paired with [rOpenSci blog guidance](https://blogguide.ropensci.org/).

## Installation

``` r
remotes::install_github("ropenscilabs/roblog")
```
---
title: "Check your post"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{checks-guidance}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

Functions to be run on the path to your blog post (rendered, not the Rmd).

## Best practice

`ro_lint_md()` should identify some potential problems and enforce: 

* the use of complete alternative descriptions for image;

* the use of Title Case for the title;

* the use of Sentence case for other headings;

* the absence of "click here" as text for links;

* the proper case (lowerCamelCase) for rOpenSci name;

* the use of Hugo shortcodes for figures;

* the use of Hugo shortcodes for tweets;

* the use of relative links for links to rOpenSci website.

You need to run `render_one` on the path to the Markdown file. 
Some Markdown examples and the corresponding messages below.

---

```{r lint, results = "asis", echo = FALSE}
render_one <- function(path) {
  result <- as.character(suppressMessages(roblog::ro_lint_md(path)))
  glue::glue_collapse(c(
    details::details(path, lang = "Markdown", summary = path),"",
  result), 
  sep = "\n")
}

glue::glue_collapse(
  purrr::map_chr(dir(system.file("examples", package = "roblog"), full.names = TRUE), 
               render_one),
  sep = "\n\n----\n\n")

```

## URL validity

`ro_check_urls()` will identify possibly broken URLs.


```{r urls}
path1 <- system.file(file.path("examples", "bad-no-alt.md"),
                                         package = "roblog")

roblog::ro_check_urls(path1)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ro_check_urls.R
\name{ro_check_urls}
\alias{ro_check_urls}
\title{Check URLs in Markdown post}
\usage{
ro_check_urls(path)
}
\arguments{
\item{path}{Path to the Markdown post (not source Rmd!)}
}
\description{
Check URLs in Markdown post
}
\examples{
\dontrun{
path <- system.file(file.path("examples", "bad-no-alt.md"),
                    package = "roblog")
ro_check_urls(path)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lint-posts.R
\name{ro_lint_md}
\alias{ro_lint_md}
\title{Lint Markdown post for rOpenSci blog}
\usage{
ro_lint_md(path)
}
\arguments{
\item{path}{Path to the Markdown post (not source Rmd!)}
}
\description{
Lint Markdown post for rOpenSci blog
}
\examples{
\dontrun{
path <- system.file(file.path("examples", "bad-no-alt.md"),
                    package = "roblog")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ro_blog_post.R
\name{ro_blog_post_author}
\alias{ro_blog_post_author}
\title{Create a new rOpenSci author file}
\usage{
ro_blog_post_author()
}
\description{
Create a new author file, in RStudio.
}
\details{
Call them via the add-in or directly or get the \href{https://ropensci-org.github.io/blog-guidance/templatemd.html}{templates online}.

In any case, an internet connection is required as templates are downloaded
fresh from their source

\code{ro_blog_post_author()} creates Markdown files,
RStudio might warn you against saving them as ".md" but ignore that.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
}
\keyword{internal}
