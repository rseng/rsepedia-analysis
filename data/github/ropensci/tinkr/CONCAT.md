
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tinkr

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R build
status](https://github.com/ropensci/tinkr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/tinkr/actions)
[![Coverage
status](https://codecov.io/gh/ropensci/tinkr/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/tinkr?branch=master)
<!-- badges: end -->

The goal of tinkr is to convert (R)Markdown files to XML and back to
allow their editing with xml2 (XPath!) instead of numerous complicated
regular expressions. If these words mean nothing to you, see our list of
[resources to get started](#background--pre-requisites).

## Use Cases

Possible applications are R scripts using tinkr, and XPath via xml2 to:

-   change levels of headers, cf [this
    script](inst/scripts/roweb2_headers.R) and [this pull request to
    roweb2](https://github.com/ropensci/roweb2/pull/279);
-   change chunk labels and options;
-   extract all runnable code (including inline code);
-   insert arbitrary Markdown elements;
-   modify link URLs;
-   your idea, please [report use
    cases](https://discuss.ropensci.org/c/usecases/10)!

## Workflow

Only the body of the (R) Markdown file is cast to XML, using the
Commonmark specification via the [`commonmark`
package](https://github.com/jeroen/commonmark). YAML metadata could be
edited using the [`yaml` package](https://github.com/viking/r-yaml),
which is not the goal of this package.

We have created an [R6 class](https://r6.r-lib.org/) object called
**yarn** to store the representation of both the YAML and the XML data,
both of which are accessible through the `$body` and `$yaml` elements.
In addition, the namespace prefix is set to “md” in the `$ns` element.

You can perform XPath queries using the `$body` and `$ns` elements:

``` r
library("tinkr")
library("xml2")
path <- system.file("extdata", "example1.md", package = "tinkr")
head(readLines(path))
#| [1] "---"                                                                               
#| [2] "title: \"What have these birds been studied for? Querying science outputs with R\""
#| [3] "slug: birds-science"                                                               
#| [4] "authors:"                                                                          
#| [5] "  - name: Maëlle Salmon"                                                           
#| [6] "    url: https://masalmon.eu/"
ex1 <- tinkr::yarn$new(path)
# find all ropensci.org blog links
xml_find_all(
  x = ex1$body, 
  xpath = ".//md:link[contains(@destination,'ropensci.org/blog')]", 
  ns = ex1$ns
)
#| {xml_nodeset (7)}
#| [1] <link destination="https://ropensci.org/blog/2018/08/21/birds-radolfzell/ ...
#| [2] <link destination="https://ropensci.org/blog/2018/09/04/birds-taxo-traits ...
#| [3] <link destination="https://ropensci.org/blog/2018/08/21/birds-radolfzell/ ...
#| [4] <link destination="https://ropensci.org/blog/2018/08/14/where-to-bird/" t ...
#| [5] <link destination="https://ropensci.org/blog/2018/08/21/birds-radolfzell/ ...
#| [6] <link destination="https://ropensci.org/blog/2018/08/28/birds-ocr/" title ...
#| [7] <link destination="https://ropensci.org/blog/2018/09/04/birds-taxo-traits ...
```

## Installation

Wanna try the package and tell us what doesn’t work yet?

``` r
install.packages("tinkr", repos = "https://ropensci.r-universe.dev")
```

## Examples

### Markdown

This is a basic example. We read “example1.md”, change all headers 3 to
headers 1, and save it back to md. Because the xml2 objects are [passed
by
reference](https://blog.penjee.com/wp-content/uploads/2015/02/pass-by-reference-vs-pass-by-value-animation.gif),
manipulating them does not require reassignment.

``` r
library("magrittr")
library("tinkr")
# From Markdown to XML
path <- system.file("extdata", "example1.md", package = "tinkr")
# Level 3 header example:
cat(tail(readLines(path, 40)), sep = "\n")
#| ### Getting a list of 50 species from occurrence data
#| 
#| For more details about the following code, refer to the [previous post
#| of the series](https://ropensci.org/blog/2018/08/21/birds-radolfzell/).
#| The single difference is our adding a step to keep only data for the
#| most recent years.
ex1  <- tinkr::yarn$new(path)
# transform level 3 headers into level 1 headers
ex1$body %>%
  xml2::xml_find_all(xpath = ".//md:heading[@level='3']", ex1$ns) %>% 
  xml2::xml_set_attr("level", 1)

# Back to Markdown
tmp <- tempfile(fileext = "md")
ex1$write(tmp)
# Level three headers are now Level one:
cat(tail(readLines(tmp, 40)), sep = "\n")
#| # Getting a list of 50 species from occurrence data
#| 
#| For more details about the following code, refer to the [previous post
#| of the series](https://ropensci.org/blog/2018/08/21/birds-radolfzell/).
#| The single difference is our adding a step to keep only data for the
#| most recent years.
unlink(tmp)
```

### R Markdown

For R Markdown files, to ease editing of chunk label and options,
`to_xml` munges the chunk info into different attributes. E.g. below you
see that `code_blocks` can have a `language`, `name`, `echo` attributes.

``` r
path <- system.file("extdata", "example2.Rmd", package = "tinkr")
rmd <- tinkr::yarn$new(path)
rmd$body
#| {xml_document}
#| <document xmlns="http://commonmark.org/xml/1.0">
#|  [1] <code_block xml:space="preserve" language="r" name="setup" include="FALS ...
#|  [2] <heading level="2">\n  <text xml:space="preserve">R Markdown</text>\n</h ...
#|  [3] <paragraph>\n  <text xml:space="preserve">This is an </text>\n  <striket ...
#|  [4] <paragraph>\n  <text xml:space="preserve">When you click the </text>\n   ...
#|  [5] <code_block xml:space="preserve" language="r" name="" eval="TRUE" echo=" ...
#|  [6] <heading level="2">\n  <text xml:space="preserve">Including Plots</text> ...
#|  [7] <paragraph>\n  <text xml:space="preserve">You can also embed plots, for  ...
#|  [8] <code_block xml:space="preserve" language="python" name="" fig.cap="&quo ...
#|  [9] <code_block xml:space="preserve" language="python" name="">plot(pressure ...
#| [10] <paragraph>\n  <text xml:space="preserve">Non-RMarkdown blocks are also  ...
#| [11] <code_block info="bash" xml:space="preserve" name="">echo "this is an un ...
#| [12] <code_block xml:space="preserve" name="">This is an ambiguous code block ...
#| [13] <paragraph>\n  <text xml:space="preserve">Note that the </text>\n  <code ...
#| [14] <table>\n  <table_header>\n    <table_cell align="left">\n      <text xm ...
#| [15] <paragraph>\n  <text xml:space="preserve">blabla</text>\n</paragraph>
```

Note that all of the features in tinkr work for both Markdown and R
Markdown.

### Inserting new Markdown elements

Inserting new nodes into the AST is surprisingly difficult if there is a
default namespace, so we have provided a method in the **yarn** object
that will take plain Markdown and translate it to XML nodes and insert
them into the document for you. For example, you can add a new code
block:

``` r
path <- system.file("extdata", "example2.Rmd", package = "tinkr")
rmd <- tinkr::yarn$new(path)
xml2::xml_find_first(rmd$body, ".//md:code_block", rmd$ns)
#| {xml_node}
#| <code_block space="preserve" language="r" name="setup" include="FALSE" eval="TRUE">
new_code <- c(
  "```{r xml-block, message = TRUE}",
  "message(\"this is a new chunk from {tinkr}\")",
  "```")
new_table <- data.frame(
  package = c("xml2", "xslt", "commonmark", "tinkr"),
  cool = TRUE
)
# Add chunk into document after the first chunk
rmd$add_md(new_code, where = 1L)
# Add a table after the second chunk:
rmd$add_md(knitr::kable(new_table), where = 2L)
# show the first 21 lines of modified document
rmd$head(21)
#| ---
#| title: "Untitled"
#| author: "M. Salmon"
#| date: "September 6, 2018"
#| output: html_document
#| ---
#| 
#| ```{r setup, include=FALSE, eval=TRUE}
#| knitr::opts_chunk$set(echo = TRUE)
#| ```
#| 
#| ```{r xml-block, message=TRUE}
#| message("this is a new chunk from {tinkr}")
#| ```
#| 
#| | package                    | cool                | 
#| | :------------------------- | :------------------ |
#| | xml2                       | TRUE                | 
#| | xslt                       | TRUE                | 
#| | commonmark                 | TRUE                | 
#| | tinkr                      | TRUE                |
```

## Background / pre-requisites

If you are not closely following one of the examples provided, what
background knowledge do you need before using tinkr?

-   That XPath, a language for querying XML & HTML, exists, and [some
    basics](https://www.w3schools.com/xml/xpath_intro.asp).
-   Basics of how [xml2
    works](https://blog.r-hub.io/2020/01/22/mutable-api/#exposing-the-c-api-in-xml2):
    how to find, replace, remove nodes etc.
-   How to use R6 classes… although reading the examples should help you
    get the gist.
-   If you are not happy with [our default
    stylesheet](#general-principles-and-solution), then understanding
    [XSLT](https://ropensci.org/blog/2017/01/10/xslt-release/) will help
    you create your own. Refer to this good resource on [XSLT for XML
    transformations](https://www.w3schools.com/xml/xsl_intro.asp).

## Loss of Markdown style

### General principles and solution

The (R)md to XML to (R)md loop on which `tinkr` is based is slightly
lossy because of Markdown syntax redundancy, so the loop from (R)md to
R(md) via `to_xml` and `to_md` will be a bit lossy. For instance

-   lists can be created with either “+”, “-” or "\*“. When using
    `tinkr`, the (R)md after editing will only use”-" for lists.

-   Links built like `[word][smallref]` with a bottom anchor
    `[smallref]: URL` will have the anchor moved to the bottom of the
    document.

-   Characters are escaped (e.g. “\[” when not for a link).

-   [x] GitHub tickboxes are preserved (only for `yarn` objects)

-   Block quotes lines all get “&gt;” whereas in the input only the
    first could have a “&gt;” at the beginning of the first line.

-   For tables see the next subsection.

Such losses make your (R)md different, and the git diff a bit harder to
parse, but should *not* change the documents your (R)md is rendered to.
If it does, report a bug in the issue tracker!

A solution to not loose your Markdown style, e.g. your preferring "\*"
over “-” for lists is to tweak [our XSL
stylesheet](inst/extdata/xml2md_gfm.xsl) and provide its filepath as
`stylesheet_path` argument to `to_md`.

### The special case of tables

-   Tables are supposed to remain/become pretty after a full loop
    `to_xml` + `to_md`. If you notice something amiss, e.g. too much
    space compared to what you were expecting, please open an issue.

### LaTeX equations

While Markdown parsers like pandoc know what LaTeX is, commonmark does
not, and that means LaTeX equations will end up with extra markup due to
commonmark’s desire to escape characters.

However, if you have LaTeX equations that use either `$` or `$$` to
delimit them, you can protect them from formatting changes with the
`$protect_math()` method (for users of the `yarn` object) or the
`protect_math()` funciton (for those using the output of `to_xml()`).
Below is a demonstration using the `yarn` object:

``` r
path <- system.file("extdata", "math-example.md", package = "tinkr")
math <- tinkr::yarn$new(path)
math$tail() # malformed
#| 
#| $$
#| Q\_{N(norm)}=\\frac{C\_N +C\_{N-1}}2\\times
#| \\frac{\\sum *{i=N-n}^{N}Q\_i} {\\sum*{j=N-n}^{N}{(\\frac{C\_j+C\_{j-1}}2)}}
#| $$
math$protect_math()$tail() # success!
#| 
#| $$
#| Q_{N(norm)}=\frac{C_N +C_{N-1}}2\times
#| \frac{\sum _{i=N-n}^{N}Q_i} {\sum_{j=N-n}^{N}{(\frac{C_j+C_{j-1}}2)}}
#| $$
```

Note, however, that there are a few caveats for this:

1.  The dollar notation for inline math must be adjacent to the text.
    E.G. `$\alpha$` is valid, but `$ \alpha$` and `$\alpha $` are not
    valid.

2.  We do not currently have support for bracket notation

3.  If you use a postfix dollar sign in your prose (e.g. BASIC commands
    or a Burroughs-Wheeler Transformation demonstration), you must be
    sure to either use punctuation after the trailing dollar sign OR
    format the text as code. (i.e. `` `INKEY$` `` is good, but `INKEY$`
    by itself is not good and will be interepreted as LaTeX code,
    throwing an error:

    ``` r
    path <- system.file("extdata", "basic-math.md", package = "tinkr")
    math <- tinkr::yarn$new(path)
    math$head(15) # malformed
    #| ---
    #| title: basic math
    #| ---
    #| 
    #| BASIC programming can make things weird:
    #| 
    #| - Give you $2 to tell me what INKEY$ means.
    #| - Give you $2 to *show* me what INKEY$ means.
    #| - Give you $2 to *show* me what `INKEY$` means.
    #| 
    #| Postfix dollars mixed with prefixed dollars can make things weird:
    #| 
    #| - We write $2 but say 2$ verbally.
    #| - We write $2 but *say* 2$ verbally.
    math$protect_math() #error
    #| Error: Inline math delimiters are not balanced.
    #| 
    #| HINT: If you are writing BASIC code, make sure you wrap variable
    #|       names and code in backtics like so: `INKEY$`. 
    #| 
    #| Below are the pairs that were found:
    #|            start...end
    #|            -----...---
    #|  Give you $2 to ... me what INKEY$ means.
    #|  Give you $2 to ... 2$ verbally.
    #| We write $2 but ...
    ```

## Meta

Please note that the ‘tinkr’ project is released with a [Contributor
Code of Conduct](https://ropensci.org/code-of-conduct). By contributing
to this project, you agree to abide by its terms.
# tinkr dev

* xml and yaml objects are now stored in an R6 class called `yarn`.
* testthat edition 3 is now being used with snapshot testing.
* Tables are now pretty after a full loop `to_xml()` + `to_md()` (@pdaengeli, #9)
* 2021-05-04: yarn objects remember the `sourcepos` and `encoding` options 
  when using the `$reset()` method.
* 2021-05-06: `protect_math()` function and method protects LaTeX math (dollar 
  notation) from escaping by commonmark (@zkamvar, #39).
* 2021-05-06: GitHub-flavored markdown ticks/checkboxes are now protected by
  default (@zkamvar, #39).
* 2021-05-11: `md_ns()` is a new convenience function to provide the `md` 
  namespace prefix for commonmark xml documents (@zkamvar, #39).
* 2021-05-11: `stylesheet()` returns the path to the internal {tinkr} stylesheet
  so that it can easily be discovered by other packages
* 2021-05-11: yarn methods `show()`, `head()`, and `tail()` all gain 
  `stylesheet_path` arguments so the modified stylesheets can be used.
* 2021-05-24: reference style links (i.e. `[text][link-ref]` with `[link-ref]: 
  <link>` on another place in the document will be preserved and the anchor will
  sink to the bottom of the document.
* 2021-09-14: numeric options fig.width and fig.height will no longer be quoted;
  `transform_params()` is simplified and no longer requires glue.
* 2021-10-15: math with embedded code and punctuation following are now allowed
  (@zkamvar #56)
* 2021-10-18: links and asis nodes that are at the beginning of paragraphs are
  no longer escaped (@zkamvar, #58)

# tinkr 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
# mal-formed inline math throws an informative error

    Inline math delimiters are not balanced.
    
    HINT: If you are writing BASIC code, make sure you wrap variable
          names and code in backtics like so: `INKEY$`. 
    
    Below are the pairs that were found:
               start...end
               -----...---
     Give you $2 to ... me what INKEY$ means.
     Give you $2 to ... 2$ verbally.
    We write $2 but ...

# block math can be protected

    Code
      show_user(m$protect_math()$tail(48), force = TRUE)
    Output
      $$
      \begin{align} % This mode aligns the equations to the '&=' signs
      \begin{split} % This mode groups the equations into one.
      \bar{r}_d &= \frac{\sum\sum{cov_{j,k}}}{
                         \sum\sum{\sqrt{var_{j} \cdot var_{k}}}} \\
                &= \frac{V_O - V_E}{2\sum\sum{\sqrt{var_{j} \cdot var_{k}}}}
      \end{split}
      \end{align}
      $$
      ```
      
      $$
      \begin{align} % This mode aligns the equations to the '&=' signs
      \begin{split} % This mode groups the equations into one.
      \bar{r}_d &= \frac{\sum\sum{cov_{j,k}}}{
      \sum\sum{\sqrt{var_{j} \cdot var_{k}}}} \
      &= \frac{V_O - V_E}{2\sum\sum{\sqrt{var_{j} \cdot var_{k}}}}
      \end{split}
      \end{align}
      $$
      
      When $a \ne 0$, there are two solutions to $ax^2 + bx + c = 0$ and they are
      
      ```latex
      $$
      x = {-b \pm \sqrt{b^2-4ac} \over 2a}
      $$
      ```
      
      $$
      x = {-b \pm \sqrt{b^2-4ac} \over 2a}
      $$
      
      Below is an example from [https://github.com/ropensci/tinkr/issues/38](https://github.com/ropensci/tinkr/issues/38)
      $\frac{\sum _{i=N-n}^{N}Q_i} {\sum_{j=N-n}^{N}{(\frac{C_j+C_{j-1}}2)}}$
      
      ```latex
      $$
      Q_{N(norm)}=\frac{C_N +C_{N-1}}2\times 
      \frac{\sum _{i=N-n}^{N}Q_i} {\sum_{j=N-n}^{N}{(\frac{C_j+C_{j-1}}2)}}
      $$
      ```
      
      $$
      Q_{N(norm)}=\frac{C_N +C_{N-1}}2\times
      \frac{\sum _{i=N-n}^{N}Q_i} {\sum_{j=N-n}^{N}{(\frac{C_j+C_{j-1}}2)}}
      $$
      

# tick boxes are protected by default

    Code
      show_user(m$head(15), force = TRUE)
    Output
      ---
      title: An example with math elements
      ---
      
      This is cheap, it only costs 10$!
      
      This example has $\LaTeX$ elements embedded in the
      text. It is intended to demonstrate that m $\alpha_\tau$ h
      mode can work with tinkr. $y =
      mx + b$
      
      - [ ] This is an empty checkbox
      - [x] This is a checked checkbox
      - [This is a link](https://ropensci.org)
      - \[this is an example\]

# anchored links are processed by default

    Code
      show_user(m$show(), force = TRUE)
    Output
      ---
      title: this tests links
      ---
      
      ## These are some links that are valid in basic markdown
      
      This is some text [that contains links][this fun link1] which
      [can be `inline`](https://example.com/2) or [can be spread across multiple lines
      because the link text is JUST TOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      LONG, y'all][link3].
      
      Mainly, we want to see if [link text
      by reference][link4] and if links [can be referenced multiple times][this fun link1]
      
      This should also [include non-reference links](https://example.com/5)
      
      If you write \[some link text\]\[link2\], that will appear as [some link text][link2]
      but you can also [test][racehorse] [sub][sub-link1] [links][sub-link2]
      
      [pizza \& icecream][pizzaicecream] = fun
      
      ```markdown
      you can write links like [a link](https://example.com/racehorse) or using
      [reference style][racehorce]
      
      [racehorse]: https://example.com/racehorse/   
      ```
      
      ## This is some extended markdown content {#extended .callout}
      
      This should also include references that use [standalone][standalone] links and
      footnotes should not be properly parsed and will be considered 'asis' nodes,
      at least that's what I *believe*\[^footy\]. Maybe this might not pan out \[^but who
      knows? footnotes are **WEIRD**, man\].
      
      <!-- links go here! -->
      
      \[^footy\]: this is a footnote that
      should be preserved
      
      [this fun link1]: https://example.com/1
      [link3]: https://example.com/3
      [link4]: https://example.com/4
      [link2]: https://example.com/2 "link with title!"
      [racehorse]: https://example.com/racehorse/
      [sub-link1]: https://example.com/racehorse/1/1 "One One Won One"
      [sub-link2]: https://example.com/racehorse/2/2/ "Two Two Won One Two"
      [pizzaicecream]: https://example.com/pizza&icecream
      [standalone]: https://example.com/standalone
      
      

# users can turn off anchor links

    Code
      show_user(m$show(), force = TRUE)
    Output
      ---
      title: this tests links
      ---
      
      ## These are some links that are valid in basic markdown
      
      This is some text [that contains links](https://example.com/1) which
      [can be `inline`](https://example.com/2) or [can be spread across multiple lines
      because the link text is JUST TOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      LONG, y'all](https://example.com/3).
      
      Mainly, we want to see if [link text
      by reference](https://example.com/4) and if links [can be referenced multiple times](https://example.com/1)
      
      This should also [include non-reference links](https://example.com/5)
      
      If you write \[some link text\]\[link2\], that will appear as [some link text](https://example.com/2 "link with title!")
      but you can also [test](https://example.com/racehorse/) [sub](https://example.com/racehorse/1/1 "One One Won One") [links](https://example.com/racehorse/2/2/ "Two Two Won One Two")
      
      [pizza \& icecream](https://example.com/pizza&icecream) = fun
      
      ```markdown
      you can write links like [a link](https://example.com/racehorse) or using
      [reference style][racehorce]
      
      [racehorse]: https://example.com/racehorse/   
      ```
      
      ## This is some extended markdown content {#extended .callout}
      
      This should also include references that use [standalone](https://example.com/standalone) links and
      footnotes should not be properly parsed and will be considered 'asis' nodes,
      at least that's what I *believe*\[^footy\]. Maybe this might not pan out \[^but who
      knows? footnotes are **WEIRD**, man\].
      
      <!-- links go here! -->
      
      \[^footy\]: this is a footnote that
      should be preserved
      

# yarn show, head, and tail methods work

    Code
      show_user(res <- y1$show(), TRUE)
    Output
      ---
      title: "Untitled"
      author: "M. Salmon"
      date: "September 6, 2018"
      output: html_document
      ---
      
      ```{r setup, include=FALSE, eval=TRUE}
      knitr::opts_chunk$set(echo = TRUE)
      ```
      
      ## R Markdown
      
      This is an ~~R Markdown document~~. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see [http://rmarkdown.rstudio.com](http://rmarkdown.rstudio.com).
      
      When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:
      
      ```{r, eval=TRUE, echo=TRUE}
      summary(cars)
      ```
      
      ## Including Plots
      
      You can also embed plots, for example:
      
      ```{python, fig.cap="pretty plot", echo=-c(1, 2), eval=TRUE}
      plot(pressure)
      ```
      
      ```{python}
      plot(pressure)
      ```
      
      Non-RMarkdown blocks are also considered
      
      ```bash
      echo "this is an unevaluted bash block"
      ```
      
      ```
      This is an ambiguous code block
      ```
      
      Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
      
      | scientific\_name            | common\_name         | n   | 
      | :------------------------- | :------------------ | --: |
      | Corvus corone              | Carrion Crow        | 288 | 
      | Turdus merula              | Eurasian Blackbird  | 285 | 
      | Anas platyrhynchos         | Mallard             | 273 | 
      | Fulica atra                | Eurasian Coot       | 268 | 
      | Parus major                | Great Tit           | 266 | 
      | Podiceps cristatus         | Great Crested Grebe | 254 | 
      | Ardea cinerea              | Gray Heron          | 236 | 
      | Cygnus olor                | Mute Swan           | 234 | 
      | Cyanistes caeruleus        | Eurasian Blue Tit   | 233 | 
      | Chroicocephalus ridibundus | Black-headed Gull   | 223 | 
      
      blabla
      

---

    Code
      show_user(res <- y1$head(10), TRUE)
    Output
      ---
      title: "Untitled"
      author: "M. Salmon"
      date: "September 6, 2018"
      output: html_document
      ---
      
      ```{r setup, include=FALSE, eval=TRUE}
      knitr::opts_chunk$set(echo = TRUE)
      ```

---

    Code
      show_user(res <- y1$tail(11), TRUE)
    Output
      | Anas platyrhynchos         | Mallard             | 273 | 
      | Fulica atra                | Eurasian Coot       | 268 | 
      | Parus major                | Great Tit           | 266 | 
      | Podiceps cristatus         | Great Crested Grebe | 254 | 
      | Ardea cinerea              | Gray Heron          | 236 | 
      | Cygnus olor                | Mute Swan           | 234 | 
      | Cyanistes caeruleus        | Eurasian Blue Tit   | 233 | 
      | Chroicocephalus ridibundus | Black-headed Gull   | 223 | 
      
      blabla
      

---
title: "What have these birds been studied for? Querying science outputs with R"
slug: birds-science
authors:
  - name: Maëlle Salmon
    url: https://masalmon.eu/
date: 2018-09-11
topicid: 1347
preface: The blog post series corresponds to the material for a talk Maëlle will give at the [Animal Movement Analysis summer school in Radolfzell, Germany on September the 12th](http://animove.org/animove-2019-evening-keynotes/), in a Max Planck Institute of Ornithology.
tags:
- rebird
- birder
- fulltext
- dataone
- EML
- literature
output:
  md_document:
    variant: markdown_github
    preserve_yaml: true
---

In the [second post of the series where we obtained data from
eBird](https://ropensci.org/blog/2018/08/21/birds-radolfzell/) we
determined what birds were observed in the county of Constance, and we
complemented this knowledge with some taxonomic and trait information in
[the fourth post of the
series](https://ropensci.org/blog/2018/09/04/birds-taxo-traits/). Now,
we could be curious about the occurrence of these birds in *scientific
work*. In this post, we will query the scientific literature and an open
scientific data repository for species names: what have these birds been
studied for? Read on if you want to learn how to use R packages allowing
to do so!

### Getting a list of 50 species from occurrence data

For more details about the following code, refer to the [previous post
of the series](https://ropensci.org/blog/2018/08/21/birds-radolfzell/).
The single difference is our adding a step to keep only data for the
most recent years.

```r
# polygon for filtering
landkreis_konstanz <- osmdata::getbb("Landkreis Konstanz",
                             format_out = "sf_polygon")
crs <- sf::st_crs(landkreis_konstanz)

# get and filter data
f_out_ebd <- "ebird/ebd_lk_konstanz.txt"

library("magrittr")

ebd <- auk::read_ebd(f_out_ebd) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"),
                crs = crs)

in_indices <- sf::st_within(ebd, landkreis_konstanz)

ebd <- dplyr::filter(ebd, lengths(in_indices) > 0)

ebd <- as.data.frame(ebd)

ebd <- dplyr::filter(ebd, approved, lubridate::year(observation_date) > 2010)
```

For the sake of simplicity, we shall only use the 50 species observed
the most often.

```r
species <- ebd %>%
  dplyr::count(common_name, sort = TRUE) %>%
  head(n = 50) %>%
  dplyr::pull(common_name)
```

The species are Carrion Crow, Eurasian Blackbird, Mallard, Eurasian
Coot, Great Tit, Great Crested Grebe, Mute Swan, Great Cormorant,
Eurasian Blue Tit, Gray Heron, Black-headed Gull, Common Chaffinch,
Common Chiffchaff, Tufted Duck, European Starling, White Wagtail,
European Robin, Little Grebe, Common Wood-Pigeon, Red-crested Pochard,
Ruddy Shelduck, Graylag Goose, Red Kite, Common Buzzard, Eurasian
Blackcap, Great Spotted Woodpecker, Eurasian Magpie, Gadwall, Common
Pochard, Eurasian Nuthatch, Green-winged Teal, House Sparrow, Eurasian
Jay, Yellow-legged Gull, Yellowhammer, Eurasian Green Woodpecker, Eared
Grebe, Eurasian Reed Warbler, Barn Swallow, Northern Shoveler, Eurasian
Moorhen, Black Redstart, Great Egret, White Stork, Eurasian Wren,
Long-tailed Tit, Common House-Martin, Eurasian Kestrel, European
Goldfinch and European Greenfinch
[(`glue::glue_collapse(species, sep = ", ", last = " and ")`)](https://twitter.com/LucyStats/status/1031938964796657665?s=19).

### Querying the scientific literature

Just like rOpenSci has a taxonomic toolbelt
([`taxize`](https://github.com/ropensci/taxize)) and a species
occurrence data toolbelt ([`spocc`](https://github.com/ropensci/spocc)),
it has a super package for querying the scientific literature:
[`fulltext`](https://github.com/ropensci/fulltext)! This package
supports search for "PLOS via the rplos package, Crossref via the
rcrossref package, Entrez via the rentrez package, arXiv via the aRxiv
package, and BMC, Biorxiv, EuroPMC, and Scopus via internal helper
functions".

We shall use `fulltext` to retrieve the titles and abstracts of
scientific articles mentioning each species, and will use `tidytext` to
compute the most prevalent words in these works.

We first define a function retrieving the titles and abstracts of works
obtained as result when querying one species name.

We use `dplyr::bind_rows` because we want all results for one species at
once, while `fulltext` returns a list of data.frames with one data.frame
by data source.

```r
.get_papers <- function(species){
  species %>%
    tolower() %>%
    fulltext::ft_search() %>%
    fulltext::ft_get() %>%
    fulltext::ft_collect() %>%
    fulltext::ft_chunks(c("title", "abstract")) %>%
    fulltext::ft_tabularize() %>%
    dplyr::bind_rows()
}

.get_papers(species[1]) %>%
  dplyr::pull(title)
```

```
##  [1] "Great spotted cuckoo nestlings have no antipredatory effect on magpie or carrion crow host nests in southern Spain"
##  [2] "Donor-Control of Scavenging Food Webs at the Land-Ocean Interface"
##  [3] "Formal comment to Soler et al.: Great spotted cuckoo nestlings have no antipredatory effect on magpie or carrion crow host nests in southern Spain"
##  [4] "Socially Driven Consistent Behavioural Differences during Development in Common Ravens and Carrion Crows"
##  [5] "Behavioral Responses to Inequity in Reward Distribution and Working Effort in Crows and Ravens"
##  [6] "Early Duplication of a Single MHC IIB Locus Prior to the Passerine Radiations"
##  [7] "Investigating the impact of media on demand for wildlife: A case study of Harry Potter and the UK trade in owls"
##  [8] "New Caledonian Crows Rapidly Solve a Collaborative Problem without Cooperative Cognition"
##  [9] "Nest Predation Deviates from Nest Predator Abundance in an Ecologically Trapped Bird"
## [10] "Dietary Compositions and Their Seasonal Shifts in Japanese Resident Birds, Estimated from the Analysis of Volunteer Monitoring Data"
```

If we were working on a scientific study, we'd add a few more filters,
e.g. having the species mentioned in the abstract, and not only
somewhere in the paper which is probably the way the different
literature search providers define a match. But we're not, so we can
keep our query quite free! My favourite paper involving the Carrion Crow
is ["Investigating the impact of media on demand for wildlife: A case
study of Harry Potter and the UK trade in
owls"](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0182368)
because it's a fun and important scientific question, and is supported
by open data (by the way you can access CITES trade data (international
trade in endangered species) in R using
[`cites`](https://github.com/ecohealthalliance/cites/) and CITES
Speciesplus database using
[`rcites`](https://ibartomeus.github.io/rcites/)).

We then apply this function to all 50 species and keep each article only
once.

```r
get_papers <- ratelimitr::limit_rate(.get_papers,
                                     rate = ratelimitr::rate(1, 2))

all_papers <- purrr::map_df(species, get_papers)

nrow(all_papers)
```

```
## [1] 522
```

```r
all_papers <- unique(all_papers)

nrow(all_papers)
```

```
## [1] 378
```

Now, we get the most common words from titles and abstracts. For that we
first append the title to the abstract which is a quick hack.

```r
library("tidytext")
library("rcorpora")

stopwords <- corpora("words/stopwords/en")$stopWords

all_papers %>%
  dplyr::group_by(title, abstract) %>%
  dplyr::summarise(text = paste(title, abstract)) %>%
  dplyr::ungroup() %>%
  unnest_tokens(word, text) %>%
  dplyr::filter(!word %in% stopwords) %>%
  dplyr::count(word, sort = TRUE) -> words
```

So, what are the most common words in these papers?

```r
head(words, n = 10)
```

```
##           word   n
## 1      species 754
## 2        birds 514
## 3        virus 270
## 4        avian 268
## 5         bird 262
## 6        study 243
## 7     breeding 231
## 8         wild 227
## 9  populations 217
## 10  population 213
```

Not too surprising, and obviously less entertaining than looking at
individual species' results. Maybe a wordcloud can give us a better idea
of the wide area of topics of studies involving our 50 most frequent
bird species. We use the [`wordcloud`
package](https://cran.r-project.org/web/packages/wordcloud/index.html).

```r
library("wordcloud")

with(words, wordcloud(word, n, max.words = 100))
```

![wordcloud of titles and abstracts of scientific
papers](/img/blog-images/2018-09-11-birds-science/wordcloud-1.png)

We see that topics include ecological words such as "foraging" but also
epidemiological questions since "influenza" and "h5n1" come up. Now, how
informative as this wordcloud can be, it's a bit ugly, so we'll prettify
it using the [`wordcloud2`
package](https://github.com/Lchiffon/wordcloud2) instead, and the
silhouette of a bird [from
Phylopic](http://phylopic.org/image/6209c9be-060e-4d7f-bc74-a75f3ccf4629/).

```r
bird <- words %>%
  head(n = 100) %>%
  wordcloud2::wordcloud2(figPath = "bird.png",
                       color = "black", size = 1.5)
# https://www.r-graph-gallery.com/196-the-wordcloud2-library/
htmlwidgets::saveWidget(bird,
                        "tmp.html",
                        selfcontained = F)
```

I wasn't able to `webshot` the resulting html despite increasing the
`delay` parameter so I screenshot it by hand!

```r
magick::image_read("screenshot.png")
```

<img src="/img/blog-images/2018-09-11-birds-science/wordcloud2-1.png" alt="wordcloud shaped as a bird" width="1366" />
<p class="caption">
wordcloud shaped as a bird
</p>

The result is a bit kitsch, doesn't include the word "species", one
needs to know it's the silhouette of a bird to recognize it, and we'd
need to work a bit on not reshaping the silhouette, but it's fun as it
is.

### Querying scientific open data

There are quite a few scientific open data repositories out there, among
which the giant [DataONE](https://www.dataone.org/) that has an API
interfaced with an R package. We shall use it to perform a search
similar to the previous section, but looking at the data indexed on
DataONE. Since DataONE specializes in ecological and environmental data,
we expect to find rather ecological data.

We first define a function to retrieve metadata of datasets for one
species. It looks the species names in the abstract.

```r
.get_meta <- function(species){

  cn <- dataone::CNode("PROD")
  search <- list(q = glue::glue("abstract:{species}"),
                        fl = "id,title,abstract",
                        sort = "dateUploaded+desc")

  result <- dataone::query(cn, solrQuery = search,
                           as="data.frame")

  if(nrow(result) == 0){
    NULL
  }else{
    # otherwise one line by version
  result <- unique(result)

  tibble::tibble(species = species,
                 title = result$title,
                 abstract = result$abstract)
  }
}
```

Note that DataONE searching could be more precise: one can choose to
search from a given data source only for instance. See the [searching
DataONE
vignette](https://github.com/DataONEorg/rdataone/blob/master/vignettes/searching-dataone.Rmd).

```r
get_meta <- ratelimitr::limit_rate(.get_meta,
                                     rate = ratelimitr::rate(1, 2))

all_meta <- purrr::map_df(species, get_meta)

nrow(all_meta)
```

```
## [1] 266
```

```r
length(unique(all_meta$species))
```

```
## [1] 35
```

35 species are represented.

```r
all_meta <- unique(all_meta[,c("title", "abstract")])

nrow(all_meta)
```

```
## [1] 104
```

We then extract the most common words.

```r
all_meta %>%
  dplyr::group_by(title, abstract) %>%
  dplyr::summarise(text = paste(title, abstract)) %>%
  dplyr::ungroup() %>%
  unnest_tokens(word, text) %>%
  dplyr::filter(!word %in% stopwords) %>%
  dplyr::count(word, sort = TRUE) -> data_words

head(data_words, n = 10)
```

```
## # A tibble: 10 x 2
##    word           n
##    <chr>      <int>
##  1 data         153
##  2 species      120
##  3 birds         94
##  4 breeding      87
##  5 feeding       75
##  6 population    65
##  7 bird          60
##  8 genetic       58
##  9 study         56
## 10 effects       54
```

Data is the most common word which is quite logical for metadata of
actual datasets. Let's also have a look at a regular wordcloud.

```r
with(data_words, wordcloud(word, n, max.words = 100))
```

![wordcloud of titles and abstracts of scientific
metadata](/img/blog-images/2018-09-11-birds-science/wordcloud3-1.png)

As expected, the words seem more focused on ecology than when looking at
scientific papers. DataONE is a gigantic data catalogue, where one could

- study the results of such queries (e.g. meta studies of number of,
  say, versions by datasets)

- or find data to integrate to a new study. If you want to *download*
  data from DataONE, refer to the [download data
  vignette](https://github.com/DataONEorg/rdataone/blob/master/vignettes/download-data.Rmd).

### Conclusion

In this post, we used the rOpenSci `fulltext` package, and the DataONE
`dataone` package, to search for bird species names in scientific papers
and scientific open datasets. We were able to draw wordclouds
representing the diversity of topics of studies in which the birds had
been mentioned or studied. Such a search could be fun to do for your
favourite bird(s)! And in general, following the same approach you could
answer your own specific research question.

#### Scientific literature access

As a reminder, the pipeline to retrieve abstracts and titles of works
mentioning a bird species was quite smooth:

```r
species %>%
    tolower() %>%
    fulltext::ft_search() %>%
    fulltext::ft_get() %>%
    fulltext::ft_collect() %>%
    fulltext::ft_chunks(c("title", "abstract")) %>%
    fulltext::ft_tabularize() %>%
    dplyr::bind_rows()
```

`fulltext` gives you a lot of power! Other rOpenSci accessing literature
data include [`europepmc`](https://github.com/ropensci/europepmc), R
Interface to Europe PMC RESTful Web Service;
[`jstor`](https://github.com/ropensci/jstor);
[`suppdata`](https://github.com/ropensci/suppdata) for extracting
supplemental information, and [much
more](https://ropensci.org/packages/).

#### Scientific data access… and publication with R

In this post we used the [`dataone`
package](https://github.com/DataONEorg/rdataone) to access data from
DataONE. That same package allows uploading data to DataONE. The
rOpenSci suite features the
[`rfigshare`](https://github.com/ropensci/rfigshare) package for getting
data from, and publishing data to, [Figshare](https://figshare.com/).
For preparing your own data and its documentation for publication, check
out the [`EML` package](https://github.com/ropensci/EML) for writing
metadata respecting the Ecological Metadata Standard, and the [unconf
`dataspice` project](https://github.com/ropenscilabs/dataspice) for
simpler metadata entry.

Explore more of our packages suite, including and beyond access to
scientific literature \&data and data publication,
[here](https://ropensci.org/packages/).

#### No more birding? No, your turn!

This was the last post of this series, that hopefully provided an
overview of how rOpenSci packages can help you learn more about birds,
and can support your workflow. As a reminder, in this series we saw

- [How to identify spots for birding using open geographical
  data](https://ropensci.org/blog/2018/08/14/where-to-bird/).
  Featuring `opencage` for geocoding, `bbox` for bounding box
  creation, `osmdata` for OpenStreetMap's Overpass API querying,
  `osmplotr` for map drawing using OpenStreetMap's data.

- [How to obtain bird occurrence data in
  R](https://ropensci.org/blog/2018/08/21/birds-radolfzell/).
  Featuring `rebird` for interaction with the eBird's API, and `auk`
  for munging of the whole eBird dataset.

- [How to extract text from old natural history
  drawings](https://ropensci.org/blog/2018/08/28/birds-ocr/).
  Featuring `magick` for image manipulation, `tesseract` for Optical
  Character Recognition, `cld2` and `cld3` for language detection, and
  `taxize::gnr_resolve` for taxonomic name resolution.

- [How to complement an occurrence dataset with taxonomy and trait
  information](https://ropensci.org/blog/2018/09/04/birds-taxo-traits/).
  Featuring `taxize`, taxonomic toolbelt for R, and `traits`,
  providing access to species traits data.

- How to query the scientific literature and scientific open data
  repositories. This is the post you've just read!

That's a wrap! But now, don't *you* hesitate to explore our packages
suite for your own needs, and to share about your use cases of rOpenSci
packages as a birder or not via [our friendly discussion
forum](https://discuss.ropensci.org/c/usecases)! Happy birding!


---
title: "What have these birds been studied for? Querying science outputs with R"
---

# TABLE HERE

[KILROY](https://en.wikipedia.org/wiki/Kilroy_was_here) WAS **HERE**

stop copying me!

STOP COPYING ME!

| A column | Another column | A third column |  |  | 
| -------- | -------------- | -------------- |  |  |
| a        | a              | aaaaaaaaaaaaaa |  |  | 
| a        | aaaaaaaaaaaaaa | a              |  |  | 
|          |                |                |  |  | 
| aaaaaaaa | a              | a              |  |  | 
|          |                |                |  |  | 

| a | b | c | d |
| : | - | :: | : |
| l | n | c | r |


---
title: "What have these birds been studied for? Querying science outputs with R"
---

| A column | Another column | A third column |  |  | 
| -------- | -------------- | -------------- |  |  |
| a        | a              | aaaaaaaaaaaaaa |  |  | 
| a        | aaaaaaaaaaaaaa | a              |  |  | 
|          |                |                |  |  | 
| aaaaaaaa | a              | a              |  |  | 
|          |                |                |  |  | 

| a | b | c | d |
| : | - | :: | : |
| l | n | c | r |


---
title: "What have these birds been studied for? Querying science outputs with R"
slug: birds-science
authors:
  - name: Maëlle Salmon
    url: https://masalmon.eu/
date: 2018-09-11
topicid: 1347
preface: The blog post series corresponds to the material for a talk Maëlle will give at the [Animal Movement Analysis summer school in Radolfzell, Germany on September the 12th](http://animove.org/animove-2019-evening-keynotes/), in a Max Planck Institute of Ornithology.
tags:
- rebird
- birder
- fulltext
- dataone
- EML
- literature
output:
  md_document:
    variant: markdown_github
    preserve_yaml: true
---

In the [second post of the series where we obtained data from
eBird](https://ropensci.org/blog/2018/08/21/birds-radolfzell/) we
determined what birds were observed in the county of Constance, and we
complemented this knowledge with some taxonomic and trait information in
[the fourth post of the
series](https://ropensci.org/blog/2018/09/04/birds-taxo-traits/). Now,
we could be curious about the occurrence of these birds in *scientific
work*. In this post, we will query the scientific literature and an open
scientific data repository for species names: what have these birds been
studied for? Read on if you want to learn how to use R packages allowing
to do so!

# Getting a list of 50 species from occurrence data

For more details about the following code, refer to the [previous post
of the series](https://ropensci.org/blog/2018/08/21/birds-radolfzell/).
The single difference is our adding a step to keep only data for the
most recent years.

```r
# polygon for filtering
landkreis_konstanz <- osmdata::getbb("Landkreis Konstanz",
                             format_out = "sf_polygon")
crs <- sf::st_crs(landkreis_konstanz)

# get and filter data
f_out_ebd <- "ebird/ebd_lk_konstanz.txt"

library("magrittr")

ebd <- auk::read_ebd(f_out_ebd) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"),
                crs = crs)

in_indices <- sf::st_within(ebd, landkreis_konstanz)

ebd <- dplyr::filter(ebd, lengths(in_indices) > 0)

ebd <- as.data.frame(ebd)

ebd <- dplyr::filter(ebd, approved, lubridate::year(observation_date) > 2010)
```

For the sake of simplicity, we shall only use the 50 species observed
the most often.

```r
species <- ebd %>%
  dplyr::count(common_name, sort = TRUE) %>%
  head(n = 50) %>%
  dplyr::pull(common_name)
```

The species are Carrion Crow, Eurasian Blackbird, Mallard, Eurasian
Coot, Great Tit, Great Crested Grebe, Mute Swan, Great Cormorant,
Eurasian Blue Tit, Gray Heron, Black-headed Gull, Common Chaffinch,
Common Chiffchaff, Tufted Duck, European Starling, White Wagtail,
European Robin, Little Grebe, Common Wood-Pigeon, Red-crested Pochard,
Ruddy Shelduck, Graylag Goose, Red Kite, Common Buzzard, Eurasian
Blackcap, Great Spotted Woodpecker, Eurasian Magpie, Gadwall, Common
Pochard, Eurasian Nuthatch, Green-winged Teal, House Sparrow, Eurasian
Jay, Yellow-legged Gull, Yellowhammer, Eurasian Green Woodpecker, Eared
Grebe, Eurasian Reed Warbler, Barn Swallow, Northern Shoveler, Eurasian
Moorhen, Black Redstart, Great Egret, White Stork, Eurasian Wren,
Long-tailed Tit, Common House-Martin, Eurasian Kestrel, European
Goldfinch and European Greenfinch
[(`glue::glue_collapse(species, sep = ", ", last = " and ")`)](https://twitter.com/LucyStats/status/1031938964796657665?s=19).

# Querying the scientific literature

Just like rOpenSci has a taxonomic toolbelt
([`taxize`](https://github.com/ropensci/taxize)) and a species
occurrence data toolbelt ([`spocc`](https://github.com/ropensci/spocc)),
it has a super package for querying the scientific literature:
[`fulltext`](https://github.com/ropensci/fulltext)! This package
supports search for "PLOS via the rplos package, Crossref via the
rcrossref package, Entrez via the rentrez package, arXiv via the aRxiv
package, and BMC, Biorxiv, EuroPMC, and Scopus via internal helper
functions".

We shall use `fulltext` to retrieve the titles and abstracts of
scientific articles mentioning each species, and will use `tidytext` to
compute the most prevalent words in these works.

We first define a function retrieving the titles and abstracts of works
obtained as result when querying one species name.

We use `dplyr::bind_rows` because we want all results for one species at
once, while `fulltext` returns a list of data.frames with one data.frame
by data source.

```r
.get_papers <- function(species){
  species %>%
    tolower() %>%
    fulltext::ft_search() %>%
    fulltext::ft_get() %>%
    fulltext::ft_collect() %>%
    fulltext::ft_chunks(c("title", "abstract")) %>%
    fulltext::ft_tabularize() %>%
    dplyr::bind_rows()
}

.get_papers(species[1]) %>%
  dplyr::pull(title)
```

```
##  [1] "Great spotted cuckoo nestlings have no antipredatory effect on magpie or carrion crow host nests in southern Spain"
##  [2] "Donor-Control of Scavenging Food Webs at the Land-Ocean Interface"
##  [3] "Formal comment to Soler et al.: Great spotted cuckoo nestlings have no antipredatory effect on magpie or carrion crow host nests in southern Spain"
##  [4] "Socially Driven Consistent Behavioural Differences during Development in Common Ravens and Carrion Crows"
##  [5] "Behavioral Responses to Inequity in Reward Distribution and Working Effort in Crows and Ravens"
##  [6] "Early Duplication of a Single MHC IIB Locus Prior to the Passerine Radiations"
##  [7] "Investigating the impact of media on demand for wildlife: A case study of Harry Potter and the UK trade in owls"
##  [8] "New Caledonian Crows Rapidly Solve a Collaborative Problem without Cooperative Cognition"
##  [9] "Nest Predation Deviates from Nest Predator Abundance in an Ecologically Trapped Bird"
## [10] "Dietary Compositions and Their Seasonal Shifts in Japanese Resident Birds, Estimated from the Analysis of Volunteer Monitoring Data"
```

If we were working on a scientific study, we'd add a few more filters,
e.g. having the species mentioned in the abstract, and not only
somewhere in the paper which is probably the way the different
literature search providers define a match. But we're not, so we can
keep our query quite free! My favourite paper involving the Carrion Crow
is ["Investigating the impact of media on demand for wildlife: A case
study of Harry Potter and the UK trade in
owls"](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0182368)
because it's a fun and important scientific question, and is supported
by open data (by the way you can access CITES trade data (international
trade in endangered species) in R using
[`cites`](https://github.com/ecohealthalliance/cites/) and CITES
Speciesplus database using
[`rcites`](https://ibartomeus.github.io/rcites/)).

We then apply this function to all 50 species and keep each article only
once.

```r
get_papers <- ratelimitr::limit_rate(.get_papers,
                                     rate = ratelimitr::rate(1, 2))

all_papers <- purrr::map_df(species, get_papers)

nrow(all_papers)
```

```
## [1] 522
```

```r
all_papers <- unique(all_papers)

nrow(all_papers)
```

```
## [1] 378
```

Now, we get the most common words from titles and abstracts. For that we
first append the title to the abstract which is a quick hack.

```r
library("tidytext")
library("rcorpora")

stopwords <- corpora("words/stopwords/en")$stopWords

all_papers %>%
  dplyr::group_by(title, abstract) %>%
  dplyr::summarise(text = paste(title, abstract)) %>%
  dplyr::ungroup() %>%
  unnest_tokens(word, text) %>%
  dplyr::filter(!word %in% stopwords) %>%
  dplyr::count(word, sort = TRUE) -> words
```

So, what are the most common words in these papers?

```r
head(words, n = 10)
```

```
##           word   n
## 1      species 754
## 2        birds 514
## 3        virus 270
## 4        avian 268
## 5         bird 262
## 6        study 243
## 7     breeding 231
## 8         wild 227
## 9  populations 217
## 10  population 213
```

Not too surprising, and obviously less entertaining than looking at
individual species' results. Maybe a wordcloud can give us a better idea
of the wide area of topics of studies involving our 50 most frequent
bird species. We use the [`wordcloud`
package](https://cran.r-project.org/web/packages/wordcloud/index.html).

```r
library("wordcloud")

with(words, wordcloud(word, n, max.words = 100))
```

![wordcloud of titles and abstracts of scientific
papers](/img/blog-images/2018-09-11-birds-science/wordcloud-1.png)

We see that topics include ecological words such as "foraging" but also
epidemiological questions since "influenza" and "h5n1" come up. Now, how
informative as this wordcloud can be, it's a bit ugly, so we'll prettify
it using the [`wordcloud2`
package](https://github.com/Lchiffon/wordcloud2) instead, and the
silhouette of a bird [from
Phylopic](http://phylopic.org/image/6209c9be-060e-4d7f-bc74-a75f3ccf4629/).

```r
bird <- words %>%
  head(n = 100) %>%
  wordcloud2::wordcloud2(figPath = "bird.png",
                       color = "black", size = 1.5)
# https://www.r-graph-gallery.com/196-the-wordcloud2-library/
htmlwidgets::saveWidget(bird,
                        "tmp.html",
                        selfcontained = F)
```

I wasn't able to `webshot` the resulting html despite increasing the
`delay` parameter so I screenshot it by hand!

```r
magick::image_read("screenshot.png")
```

<img src="/img/blog-images/2018-09-11-birds-science/wordcloud2-1.png" alt="wordcloud shaped as a bird" width="1366" />
<p class="caption">
wordcloud shaped as a bird
</p>

The result is a bit kitsch, doesn't include the word "species", one
needs to know it's the silhouette of a bird to recognize it, and we'd
need to work a bit on not reshaping the silhouette, but it's fun as it
is.

# Querying scientific open data

There are quite a few scientific open data repositories out there, among
which the giant [DataONE](https://www.dataone.org/) that has an API
interfaced with an R package. We shall use it to perform a search
similar to the previous section, but looking at the data indexed on
DataONE. Since DataONE specializes in ecological and environmental data,
we expect to find rather ecological data.

We first define a function to retrieve metadata of datasets for one
species. It looks the species names in the abstract.

```r
.get_meta <- function(species){

  cn <- dataone::CNode("PROD")
  search <- list(q = glue::glue("abstract:{species}"),
                        fl = "id,title,abstract",
                        sort = "dateUploaded+desc")

  result <- dataone::query(cn, solrQuery = search,
                           as="data.frame")

  if(nrow(result) == 0){
    NULL
  }else{
    # otherwise one line by version
  result <- unique(result)

  tibble::tibble(species = species,
                 title = result$title,
                 abstract = result$abstract)
  }
}
```

Note that DataONE searching could be more precise: one can choose to
search from a given data source only for instance. See the [searching
DataONE
vignette](https://github.com/DataONEorg/rdataone/blob/master/vignettes/searching-dataone.Rmd).

```r
get_meta <- ratelimitr::limit_rate(.get_meta,
                                     rate = ratelimitr::rate(1, 2))

all_meta <- purrr::map_df(species, get_meta)

nrow(all_meta)
```

```
## [1] 266
```

```r
length(unique(all_meta$species))
```

```
## [1] 35
```

35 species are represented.

```r
all_meta <- unique(all_meta[,c("title", "abstract")])

nrow(all_meta)
```

```
## [1] 104
```

We then extract the most common words.

```r
all_meta %>%
  dplyr::group_by(title, abstract) %>%
  dplyr::summarise(text = paste(title, abstract)) %>%
  dplyr::ungroup() %>%
  unnest_tokens(word, text) %>%
  dplyr::filter(!word %in% stopwords) %>%
  dplyr::count(word, sort = TRUE) -> data_words

head(data_words, n = 10)
```

```
## # A tibble: 10 x 2
##    word           n
##    <chr>      <int>
##  1 data         153
##  2 species      120
##  3 birds         94
##  4 breeding      87
##  5 feeding       75
##  6 population    65
##  7 bird          60
##  8 genetic       58
##  9 study         56
## 10 effects       54
```

Data is the most common word which is quite logical for metadata of
actual datasets. Let's also have a look at a regular wordcloud.

```r
with(data_words, wordcloud(word, n, max.words = 100))
```

![wordcloud of titles and abstracts of scientific
metadata](/img/blog-images/2018-09-11-birds-science/wordcloud3-1.png)

As expected, the words seem more focused on ecology than when looking at
scientific papers. DataONE is a gigantic data catalogue, where one could

- study the results of such queries (e.g. meta studies of number of,
  say, versions by datasets)

- or find data to integrate to a new study. If you want to *download*
  data from DataONE, refer to the [download data
  vignette](https://github.com/DataONEorg/rdataone/blob/master/vignettes/download-data.Rmd).

# Conclusion

In this post, we used the rOpenSci `fulltext` package, and the DataONE
`dataone` package, to search for bird species names in scientific papers
and scientific open datasets. We were able to draw wordclouds
representing the diversity of topics of studies in which the birds had
been mentioned or studied. Such a search could be fun to do for your
favourite bird(s)! And in general, following the same approach you could
answer your own specific research question.

#### Scientific literature access

As a reminder, the pipeline to retrieve abstracts and titles of works
mentioning a bird species was quite smooth:

```r
species %>%
    tolower() %>%
    fulltext::ft_search() %>%
    fulltext::ft_get() %>%
    fulltext::ft_collect() %>%
    fulltext::ft_chunks(c("title", "abstract")) %>%
    fulltext::ft_tabularize() %>%
    dplyr::bind_rows()
```

`fulltext` gives you a lot of power! Other rOpenSci accessing literature
data include [`europepmc`](https://github.com/ropensci/europepmc), R
Interface to Europe PMC RESTful Web Service;
[`jstor`](https://github.com/ropensci/jstor);
[`suppdata`](https://github.com/ropensci/suppdata) for extracting
supplemental information, and [much
more](https://ropensci.org/packages/).

#### Scientific data access… and publication with R

In this post we used the [`dataone`
package](https://github.com/DataONEorg/rdataone) to access data from
DataONE. That same package allows uploading data to DataONE. The
rOpenSci suite features the
[`rfigshare`](https://github.com/ropensci/rfigshare) package for getting
data from, and publishing data to, [Figshare](https://figshare.com/).
For preparing your own data and its documentation for publication, check
out the [`EML` package](https://github.com/ropensci/EML) for writing
metadata respecting the Ecological Metadata Standard, and the [unconf
`dataspice` project](https://github.com/ropenscilabs/dataspice) for
simpler metadata entry.

Explore more of our packages suite, including and beyond access to
scientific literature \&data and data publication,
[here](https://ropensci.org/packages/).

#### No more birding? No, your turn!

This was the last post of this series, that hopefully provided an
overview of how rOpenSci packages can help you learn more about birds,
and can support your workflow. As a reminder, in this series we saw

- [How to identify spots for birding using open geographical
  data](https://ropensci.org/blog/2018/08/14/where-to-bird/).
  Featuring `opencage` for geocoding, `bbox` for bounding box
  creation, `osmdata` for OpenStreetMap's Overpass API querying,
  `osmplotr` for map drawing using OpenStreetMap's data.

- [How to obtain bird occurrence data in
  R](https://ropensci.org/blog/2018/08/21/birds-radolfzell/).
  Featuring `rebird` for interaction with the eBird's API, and `auk`
  for munging of the whole eBird dataset.

- [How to extract text from old natural history
  drawings](https://ropensci.org/blog/2018/08/28/birds-ocr/).
  Featuring `magick` for image manipulation, `tesseract` for Optical
  Character Recognition, `cld2` and `cld3` for language detection, and
  `taxize::gnr_resolve` for taxonomic name resolution.

- [How to complement an occurrence dataset with taxonomy and trait
  information](https://ropensci.org/blog/2018/09/04/birds-taxo-traits/).
  Featuring `taxize`, taxonomic toolbelt for R, and `traits`,
  providing access to species traits data.

- How to query the scientific literature and scientific open data
  repositories. This is the post you've just read!

That's a wrap! But now, don't *you* hesitate to explore our packages
suite for your own needs, and to share about your use cases of rOpenSci
packages as a birder or not via [our friendly discussion
forum](https://discuss.ropensci.org/c/usecases)! Happy birding!


<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If you've updated a file in the man-roxygen directory, make sure to update the man/ files by running devtools::document() or similar as .Rd files should be affected by your change -->

<!--- Provide a general summary of your changes in the Title above -->

## Description
<!--- Describe your changes in detail -->

## Related Issue
<!--- if this closes an issue make sure include e.g., "fix #4"
or similar - or if just relates to an issue make sure to mention
it like "#4" -->

## Example
<!--- if introducing a new feature or changing behavior of existing
methods/functions, include an example if possible to do in brief form -->

<!--- Did you remember to include tests? Unless you're just changing
grammar, please include new tests for your change -->
# CONTRIBUTING #

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the Travis and AppVeyor build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.  
*  We recommend the tidyverse [style guide](http://style.tidyverse.org).
You can use the [styler](https://CRAN.R-project.org/package=styler) package to
apply these styles, but please don't restyle code that has nothing to do with 
your PR.  
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the tinkr project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://ropensci.github.io/dev_guide/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
Submitting Author: Maëlle Salmon (@maelle)  
Other Package Authors: Zhian N. Kamvar (@zkamvar)
Repository:  <!--repourl-->https://github.com/ropensci/tinkr/<!--end-repourl-->
Submission type: <!--submission-type-->Pre-submission<!--end-submission-type-->

---

-   Paste the full DESCRIPTION file inside a code block below:

```
Package: tinkr
Title: Casts (R)Markdown Files to XML and Back Again
Version: 0.0.0.9000
Authors@R: 
    c(person(given = "Maëlle",
             family = "Salmon",
             role = c("aut", "cre"),
             email = "msmaellesalmon@gmail.com",
             comment = c(ORCID = "0000-0002-2815-0399")),
      person(given = "Zhian N.",
             family = "Kamvar",
             role = "aut",
             email = "zkamvar@gmail.com",
             comment = c(ORCID = "0000-0003-1458-7108")),
      person(given = "Jeroen",
             family = "Ooms",
             role = "aut"),
      person(given = "Nick",
             family = "Wellnhofer",
             role = "cph",
             comment = "Nick Wellnhofer wrote the XSLT stylesheet."),
      person(given = "rOpenSci",
             role = "fnd",
             comment = "https://ropensci.org/"),
      person(given = "Peter",
             family = "Daengeli",
             role = "ctb"))
Description: Casts (R)Markdown files to XML and back to allow their
    editing via XPath.
License: GPL-3
URL: https://docs.ropensci.org/tinkr/, https://github.com/ropensci/tinkr
BugReports: https://github.com/ropensci/tinkr/issues
Imports: 
    commonmark (>= 1.6),
    fs,
    glue,
    knitr,
    magrittr,
    purrr,
    R6,
    stringr,
    xml2,
    xslt,
    yaml
Suggests: 
    covr,
    testthat (>= 3.0.0),
    withr
Config/testthat/edition: 3
Encoding: UTF-8
LazyData: true
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.1.1.9001
```


## Scope 

- Please indicate which category or categories from our [package fit policies](https://ropensci.github.io/dev_guide/policies.html#package-categories) this package falls under: (Please check an appropriate box below.:

	- [ ] data retrieval
	- [ ] data extraction
	- [ ] database access
	- [x] data munging
	- [ ] data deposition
	- [ ] workflow automation
	- [ ] version control
	- [ ] citation management and bibliometrics
	- [ ] scientific software wrappers
	- [ ] database software bindings
	- [ ] geospatial data
	- [ ] text analysis
	

- Explain how and why the package falls under these categories (briefly, 1-2 sentences).  Please note any areas you are unsure of:

With tinkr one can extract structured data out of R Markdown or Markdown files with XPath, instead of regular expressions.
A further application, that is however not in scope as it might be viewed as "general tools for literate programming", is _modifying_ such files, e.g. adding a code chunk to a bunch of R Markdown files at once.

-   Who is the target audience and what are scientific applications of this package?  

The target audience would be any scientist using R Markdown or Markdown files as data source. Salient examples of applications are extracting/modifying links from markdown documents (e.g. for CRAN checks), analyzing patterns of markdown features in documents across repositories (e.g. https://carpentries.github.io/pegboard/articles/swc-survey.html#summary-of-solutions-1), and transformation of markdown documents in a systematic way (e.g. https://carpentries.github.io/pegboard/#manipulation).

-   Are there other R packages that accomplish the same thing? If so, how does yours differ or meet [our criteria for best-in-category](https://ropensci.github.io/dev_guide/policies.html#overlap)?

The [parsermd](https://rundel.github.io/parsermd/articles/parsermd.html) package by Colin Rundel aims at "extracting the content of an R Markdown file to allow for programmatic interactions with the document’s contents (i.e. code chunks and markdown text)".
However, parsermd is focused on R Markdown documents, and as written in its docs "The goal is to capture the fundamental structure of the document and as such we do not attempt to parse every detail of the Rmd" whereas tinkr parses everything into XML according to the commonmark style.

-   (If applicable) Does your package comply with our [guidance around _Ethics, Data Privacy and Human Subjects Research_](https://devguide.ropensci.org/policies.html#ethics-data-privacy-and-human-subjects-research)?

Not applicable.

-  Any other questions or issues we should be aware of?:

No.
> the simple example of a blockquote
> the simple example of a blockquote
> the simple example of a blockquote
> the simple example of a blockquote
> ... continuation
> ... continuation
> ... continuation
> ... continuation

empty blockquote:

> 
> > > > > > deeply nested blockquote
> > > > > > deeply nested blockquote
> > > > > > deeply nested blockquote
> > > > > > deeply nested blockquote
> > > > > > deeply nested blockquote
> > > > > > deeply nested blockquote

> deeply nested blockquote
> 
> > deeply nested blockquote
> > 
> > > deeply nested blockquote
> > > 
> > > > deeply nested blockquote
> > > > 
> > > > > deeply nested blockquote
> > > > > 
> > > > > > deeply nested blockquote

```
    an
    example

    of



    a code
    block
```

````text
an
example
```
of


a fenced
```
code
block
````

# heading

### heading

##### heading

# heading

### heading

##### heading ##########

\############ not a heading

***

***

***

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\* text

<div class="this is an html block">

blah blah

</div>

<table>
  <tr>
    <td>
      **test**
    </td>
  </tr>
</table>

<table>

  <tr>

```
<td>

  test

</td>
```

  </tr>

</table>

<![CDATA[
  [[[[[[[[[[[... *cdata section - this should not be parsed* ...]]]]]]]]]]]
]]>

## heading

# heading

not a heading
\----------------------------------- text

- tidy

- bullet

- list

- loose

- bullet

- list

0. ordered
1. list
2. example

- - - - 
1. 2. 3. 
- an example
  of a list item
  with a continuation
  
  this part is inside the list

this part is just a paragraph

1. test

- test

1. test

- test

111111111111111111111111111111111111111111. is this a valid bullet?

- ***

- this

- is
  
  a
  
  long

- loose

- list

- with

- some
  
  tidy

- list

- items

- in

- between

- ***

- this
  
  - is
    - a
      - deeply
        - nested
          - bullet
            - list

1. this
   2\. is
   3\. a
   4\. deeply
   5\. nested
   6\. unordered
   7\. list

- 1

- 2

- 3

- 4

- 5

- 6

- 7

- 6

- 5

- 4

- 3

- 2

- 1

- - - - - - - - - deeply-nested one-element item

[1](http://something.example.com/foo/bar) [2](http://something.example.com/foo/bar "test") [3](http://foo/bar) [1](http://something.example.com/foo/bar) [2](http://something.example.com/foo/bar "test") [3](http://foo/bar)

[looooooooooooooooooooooooooooooooooooooooooooooooooong label](111 "test")

\[\[\[\[\[\[\[\[\[\[\[\[\[\[\[\[\[\[\[\[ this should not slow down anything \]\]\]\]\]\]\]\]\]\]\]\]\]\]\]\]\]\]\]\]: q
(as long as it is not referenced anywhere)

\[\[\[\[\[\[\[\[\[\[\[\[\[\[\[\[\[\[\[\[\]: this is not a valid reference
\[\[\[\[\[\[\[foo\]\]\]\]\]\]\]

\[\[\[\[\[\[\[foo\]\]\]\]\]\]\]: bar
\[\[\[\[\[\[foo\]\]\]\]\]\]: bar
\[\[\[\[\[foo\]\]\]\]\]: bar
\[\[\[\[foo\]\]\]\]: bar
\[\[\[foo\]\]\]: bar
\[\[foo\]\]: bar
\[foo\]: bar

\[*\[*\[*\[*\[foo\]*\]*\]*\]*\]

\[*\[*\[*\[*\[foo\]*\]*\]*\]*\]: bar
\[*\[*\[*\[foo\]*\]*\]*\]: bar
\[*\[*\[foo\]*\]*\]: bar
\[*\[foo\]*\]: bar
\[foo\]: bar
closed (valid) autolinks:

[ftp://1.2.3.4:21/path/foo](ftp://1.2.3.4:21/path/foo)
[http://foo.bar.baz?q=hello\&id=22\&boolean](http://foo.bar.baz?q=hello&id=22&boolean)
[http://veeeeeeeeeeeeeeeeeeery.loooooooooooooooooooooooooooooooong.autolink/](http://veeeeeeeeeeeeeeeeeeery.loooooooooooooooooooooooooooooooong.autolink/)
[teeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeest@gmail.com](mailto:teeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeest@gmail.com)

these are not autolinks:

\<[ftp://1.2.3.4:21/path/foo](ftp://1.2.3.4:21/path/foo)
\<[http://foo.bar.baz?q=hello\&id=22\&boolean](http://foo.bar.baz?q=hello&id=22&boolean)
\<[http://veeeeeeeeeeeeeeeeeeery.loooooooooooooooooooooooooooooooong.autolink](http://veeeeeeeeeeeeeeeeeeery.loooooooooooooooooooooooooooooooong.autolink)
\<[teeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeest@gmail.com](mailto:teeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeest@gmail.com)\
\< [http://foo.bar.baz?q=hello\&id=22\&boolean](http://foo.bar.baz?q=hello&id=22&boolean) >
`lots`of`backticks`

`i`wonder`how`this`will`be`parsed`
*this* *is* *your* *basic* *boring* *emphasis*

*this* *is* *your* *basic* *boring* *emphasis*

**this** **is** **your** **basic** **boring** **emphasis**
*this *is *a *bunch* of* nested* emphases*

**this **is **a **bunch** of** nested** emphases**

***this ***is ***a ***bunch*** of*** nested*** emphases***\
\*this \*is \*a \*worst \*case \*for \*em \*backtracking

\_\_this \_\_is \_\_a \_\_worst \_\_case \_\_for \_\_em \_\_backtracking

\*\*\*this \*\*\*is \*\*\*a \*\*\*worst \*\*\*case \*\*\*for \*\*\*em \*\*\*backtracking
entities:

  \& © Æ Ď ¾ ℋ ⅆ ∲

\# Ӓ Ϡ �

non-entities:

\&18900987654321234567890; \&1234567890098765432123456789009876543212345678987654;

\&qwertyuioppoiuytrewqwer; \&oiuytrewqwertyuioiuytrewqwertyuioytrewqwertyuiiuytri;

\\t\\e\\s\\t\\i\\n\\g \\e\\s\\c\\a\\p\\e \\s\\e\\q\\u\\e\\n\\c\\e\\s

!\\"#$%\&'()\*+,./:;\<=>?

@ \[ \] ^ \_ \` { | } ~ - '

  
\\
\\  
\\\\
\\\\\\

\<this> \<is> \<not> \<html>

Taking commonmark tests from the spec for benchmarking here:

<a><bab><c2c>

<a/><b2/>

<a  /><b2
data="foo" >

<a foo="bar" bam = 'baz <em>"</em>'
_boolean zoop:33=zoop:33 />

\<33> \<\_\_>

\<a h\*#ref="hi">

\<a href="hi'> \<a href=hi'>

\< a>\<
foo>\<bar/ >

\<a href='bar'title=title>

</a>
</foo >

\</a href="foo">

foo <!-- this is a
comment - with hyphen -->

foo \<!-- not a comment -- two hyphens -->

foo <?php echo $a; ?>

foo <!ELEMENT br EMPTY>

foo <![CDATA[>&<]]>

<a href="&ouml;">

<a href="\*">

\<a href=""">
Valid links:

[this is a link]()
[this is a link](http://something.example.com/foo/bar)
[this is a link](http://something.example.com/foo/bar "test")
![this is an image]()
![this is an image](http://something.example.com/foo/bar)
![this is an image](http://something.example.com/foo/bar "test")

[escape test](>>>>>>>>>>>>>> "''''''''''''''")
[escape test \]\]\]\]\]\]\]\]\]\]\]\]\]\]\]\]](\)\)\)\)\)\)\)\)\)\)\)\)\)\))

Invalid links:

\[this is not a link

\[this is not a link\](

\[this is not a link\]([http://something.example.com/foo/bar](http://something.example.com/foo/bar) 'test'

\[this is not a link\](((((((((((((((((((((((((((((((((((((((((((((((

[this is not a link](\(\(\(\(\(\(\(\(\(\(\)\)\)\)\)\)\)\)\)\) "((((((((("))))))))))
Valid links:

\[\[\[\[\[\[\[[](test)\](test)\](test)\](test)\](test)\](test)\](test)\]

\[ \[\[\[\[\[\[\[\[\[\[\[\[\[\[\[\[\[\[ [](test) \]\]\]\]\]\]\]\]\]\]\]\]\]\]\]\]\]\] \](test)

Invalid links:

\[\[\[\[\[\[\[\[\[

\[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[ \[

!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[!\[

this  
should  
be  
separated  
by  
newlines

this  
should  
be  
separated  
by  
newlines  
too

this
should
not
be
separated
by
newlines

Lorem ipsum dolor sit amet, **consectetur** adipiscing elit. Cras imperdiet nec erat ac condimentum. Nulla vel rutrum ligula. Sed hendrerit interdum orci a posuere. Vivamus ut velit aliquet, mollis purus eget, iaculis nisl. Proin posuere malesuada ante. Proin auctor orci eros, ac molestie lorem dictum nec. Vestibulum sit amet erat est. Morbi luctus sed elit ac luctus. Proin blandit, enim vitae egestas posuere, neque elit ultricies dui, vel mattis nibh enim ac lorem. Maecenas molestie nisl sit amet velit dictum lobortis. Aliquam erat volutpat.

Vivamus sagittis, diam in [vehicula](https://github.com/markdown-it/markdown-it) lobortis, sapien arcu mattis erat, vel aliquet sem urna et risus. Ut feugiat sapien vitae mi elementum laoreet. Suspendisse potenti. Aliquam erat nisl, aliquam pretium libero aliquet, sagittis eleifend nunc. In hac habitasse platea dictumst. Integer turpis augue, tincidunt dignissim mauris id, rhoncus dapibus purus. Maecenas et enim odio. Nullam massa metus, varius quis vehicula sed, pharetra mollis erat. In quis viverra velit. Vivamus placerat, est nec hendrerit varius, enim dui hendrerit magna, ut pulvinar nibh lorem vel lacus. Mauris a orci iaculis, hendrerit eros sed, gravida leo. In dictum mauris vel augue varius, ac ullamcorper nisl ornare. In eu posuere velit, ac fermentum arcu. Interdum et malesuada fames ac ante ipsum primis in faucibus. Nullam sed malesuada leo, at interdum elit.

Nullam ut tincidunt nunc. [Pellentesque](http://something.example.com/foo/bar) metus lacus, commodo eget justo ut, rutrum varius nunc. Sed non rhoncus risus. Morbi sodales gravida pulvinar. Duis malesuada, odio volutpat elementum vulputate, massa magna scelerisque ante, et accumsan tellus nunc in sem. Donec mattis arcu et velit aliquet, non sagittis justo vestibulum. Suspendisse volutpat felis lectus, nec consequat ipsum mattis id. Donec dapibus vehicula facilisis. In tincidunt mi nisi, nec faucibus tortor euismod nec. Suspendisse ante ligula, aliquet vitae libero eu, vulputate dapibus libero. Sed bibendum, sapien at posuere interdum, libero est sollicitudin magna, ac gravida tellus purus eu ipsum. Proin ut quam arcu.

Suspendisse potenti. Donec ante velit, ornare at augue quis, tristique laoreet sem. Etiam in ipsum elit. Nullam cursus dolor sit amet nulla feugiat tristique. Phasellus ac tellus tincidunt, imperdiet purus eget, ullamcorper ipsum. Cras eu tincidunt sem. Nullam sed dapibus magna. Lorem ipsum dolor sit amet, consectetur adipiscing elit. In id venenatis tortor. In consectetur sollicitudin pharetra. Etiam convallis nisi nunc, et aliquam turpis viverra sit amet. Maecenas faucibus sodales tortor. Suspendisse lobortis mi eu leo viverra volutpat. Pellentesque velit ante, vehicula sodales congue ut, elementum a urna. Cras tempor, ipsum eget luctus rhoncus, arcu ligula fermentum urna, vulputate pharetra enim enim non libero.

Proin diam quam, elementum in eleifend id, elementum et metus. Cras in justo consequat justo semper ultrices. Sed dignissim lectus a ante mollis, nec vulputate ante molestie. Proin in porta nunc. Etiam pulvinar turpis sed velit porttitor, vel adipiscing velit fringilla. Cras ac tellus vitae purus pharetra tincidunt. Sed cursus aliquet aliquet. Cras eleifend commodo malesuada. In turpis turpis, ullamcorper ut tincidunt a, ullamcorper a nunc. Etiam luctus tellus ac dapibus gravida. Ut nec lacus laoreet neque ullamcorper volutpat.

Nunc et leo erat. Aenean mattis ultrices lorem, eget adipiscing dolor ultricies eu. In hac habitasse platea dictumst. Vivamus cursus feugiat sapien quis aliquam. Mauris quam libero, porta vel volutpat ut, blandit a purus. Vivamus vestibulum dui vel tortor molestie, sit amet feugiat sem commodo. Nulla facilisi. Sed molestie arcu eget tellus vestibulum tristique.

this is a test for tab expansion, be careful not to replace them with spaces

1	4444
22	333
333	22
4444	1

```
tab-indented line
space-indented line
tab-indented line
```

a lot of                                                spaces in between here

a lot of												tabs in between here

| scientific\_name| common\_name| n| 
 |  --- | --- | --- |
| Corvus corone| Carrion Crow| 288| 
| Turdus merula| Eurasian Blackbird| 285| 
| Anas platyrhynchos| Mallard| 273| 
| Fulica atra| Eurasian Coot| 268| 
| Parus major| Great Tit| 266| 
| Podiceps cristatus| Great Crested Grebe| 254| 
| Ardea cinerea| Gray Heron| 236| 
| Cygnus olor| Mute Swan| 234| 
| Cyanistes caeruleus| Eurasian Blue Tit| 233| 
| Chroicocephalus ridibundus| Black-headed Gull| 223| 

~~blabla~~


---
title: "What have these birds been studied for? Querying science outputs with R"
---



| A column | Another column | A third column |   |   |
|----------|----------------|----------------|---|---|
| a        | a              | aaaaaaaaaaaaaa |   |   |
| a        | aaaaaaaaaaaaaa | a              |   |   |
|          |                |                |   |   |
| aaaaaaaa | a              | a              |   |   |
|          |                |                |   |   |


| a | b | c | d | 
| : | - | :: | : |
| l | n | c | r | 
---
title: basic math
---

BASIC programming can make things weird:

 - Give you $2 to tell me what INKEY$ means.
 - Give you $2 to _show_ me what INKEY$ means.
 - Give you $2 to _show_ me what `INKEY$` means.

Postfix dollars mixed with prefixed dollars can make things weird:

 - We write $2 but say 2$ verbally.
 - We write $2 but _say_ 2$ verbally.
---
title: this tests links
---


## These are some links that are valid in basic markdown

This is some text [that contains links][this fun link1] which 
[can be `inline`](https://example.com/2) or [can be spread across multiple lines
because the link text is JUST TOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
LONG, y'all](https://example.com/3).

Mainly, we want to see if [link text
by reference][link4] and if links [can be referenced multiple times][this fun link1]

This should also [include non-reference links](https://example.com/5)

If you write \[some link text\]\[link2\], that will appear as [some link text][link2]
but you can also [test][racehorse] [sub][sub-link1] [links][sub-link2]

[pizza & icecream][pizzaicecream] = fun

```markdown
you can write links like [a link](https://example.com/racehorse) or using
[reference style][racehorce]

[racehorse]: https://example.com/racehorse/   
```


[this fun link1]: https://example.com/1
[link2]: https://example.com/2 "link with title!"
[link3]: https://example.com/3
[link4]: https://example.com/4
[racehorse]: https://example.com/racehorse/   
[sub-link1]: https://example.com/racehorse/1/1 "One One Won One"
[sub-link2]: https://example.com/racehorse/2/2/ "Two Two Won One Two"
[pizzaicecream]: https://example.com/pizza&icecream

## This is some extended markdown content {#extended .callout}

This should also include references that use [standalone] links and 
footnotes should not be properly parsed and will be considered 'asis' nodes,
at least that's what I *believe*[^footy]. Maybe this might not pan out [^but who
knows? footnotes are **WEIRD**, man].

<!-- links go here! -->

[standalone]: https://example.com/standalone
[^footy]: this is a footnote that
  should be preserved
> the simple example of a blockquote 
> the simple example of a blockquote
> the simple example of a blockquote
> the simple example of a blockquote
... continuation
... continuation
... continuation
... continuation

empty blockquote:

>
>
>
>

>>>>>> deeply nested blockquote
>>>>> deeply nested blockquote
>>>> deeply nested blockquote
>>> deeply nested blockquote
>> deeply nested blockquote
> deeply nested blockquote

> deeply nested blockquote
>> deeply nested blockquote
>>> deeply nested blockquote
>>>> deeply nested blockquote
>>>>> deeply nested blockquote
>>>>>> deeply nested blockquote

        an
        example

        of



        a code
        block


``````````text
an
example
```
of


a fenced
```
code
block
``````````

# heading
### heading
##### heading

# heading #
### heading ###
##### heading \#\#\#\#\######

############ not a heading

 * * * * *

 -  -  -  -  -

 ________


 ************************* text

<div class="this is an html block">

blah blah

</div>

<table>
  <tr>
    <td>
      **test**
    </td>
  </tr>
</table>

<table>

  <tr>

    <td>

      test

    </td>

  </tr>

</table>

<![CDATA[
  [[[[[[[[[[[... *cdata section - this should not be parsed* ...]]]]]]]]]]]
]]>

heading
---

heading
===================================

not a heading
----------------------------------- text
 - tidy
 - bullet
 - list


 - loose

 - bullet

 - list


 0. ordered
 1. list
 2. example


 -
 -
 -
 -


 1.
 2.
 3.


 -  an example
of a list item
       with a continuation

    this part is inside the list

   this part is just a paragraph  


 1. test
 -  test
 1. test
 -  test


111111111111111111111111111111111111111111. is this a valid bullet?

 - _________________________

 - this
 - is

   a

   long
 - loose
 - list

 - with
 - some

   tidy

 - list
 - items
 - in

 - between
 - _________________________

 - this
   - is
     - a
       - deeply
         - nested
           - bullet
             - list
   

 1. this
    2. is
       3. a
          4. deeply
             5. nested
                6. unordered
                   7. list


 - 1
  - 2
   - 3
    - 4
     - 5
      - 6
       - 7
      - 6
     - 5
    - 4
   - 3
  - 2
 - 1


 - - - - - - - - - deeply-nested one-element item

[1] [2] [3] [1] [2] [3]

[looooooooooooooooooooooooooooooooooooooooooooooooooong label]

 [1]: <http://something.example.com/foo/bar>
 [2]: http://something.example.com/foo/bar 'test'
 [3]:
 http://foo/bar
 [    looooooooooooooooooooooooooooooooooooooooooooooooooong   label    ]:
 111
 'test'
 [[[[[[[[[[[[[[[[[[[[ this should not slow down anything ]]]]]]]]]]]]]]]]]]]]: q
 (as long as it is not referenced anywhere)

 [[[[[[[[[[[[[[[[[[[[]: this is not a valid reference
[[[[[[[foo]]]]]]]

[[[[[[[foo]]]]]]]: bar
[[[[[[foo]]]]]]: bar
[[[[[foo]]]]]: bar
[[[[foo]]]]: bar
[[[foo]]]: bar
[[foo]]: bar
[foo]: bar

[*[*[*[*[foo]*]*]*]*]

[*[*[*[*[foo]*]*]*]*]: bar
[*[*[*[foo]*]*]*]: bar
[*[*[foo]*]*]: bar
[*[foo]*]: bar
[foo]: bar
closed (valid) autolinks:

 <ftp://1.2.3.4:21/path/foo>
 <http://foo.bar.baz?q=hello&id=22&boolean>
 <http://veeeeeeeeeeeeeeeeeeery.loooooooooooooooooooooooooooooooong.autolink/>
 <teeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeest@gmail.com>

these are not autolinks:

 <ftp://1.2.3.4:21/path/foo
 <http://foo.bar.baz?q=hello&id=22&boolean
 <http://veeeeeeeeeeeeeeeeeeery.loooooooooooooooooooooooooooooooong.autolink
 <teeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeest@gmail.com
 < http://foo.bar.baz?q=hello&id=22&boolean >
`lots`of`backticks`

``i``wonder``how``this``will``be``parsed``
*this* *is* *your* *basic* *boring* *emphasis*

_this_ _is_ _your_ _basic_ _boring_ _emphasis_

**this** **is** **your** **basic** **boring** **emphasis**
*this *is *a *bunch* of* nested* emphases* 

__this __is __a __bunch__ of__ nested__ emphases__ 

***this ***is ***a ***bunch*** of*** nested*** emphases*** 
*this *is *a *worst *case *for *em *backtracking

__this __is __a __worst __case __for __em __backtracking

***this ***is ***a ***worst ***case ***for ***em ***backtracking
entities:

&nbsp; &amp; &copy; &AElig; &Dcaron; &frac34; &HilbertSpace; &DifferentialD; &ClockwiseContourIntegral;

&#35; &#1234; &#992; &#98765432;

non-entities:

&18900987654321234567890; &1234567890098765432123456789009876543212345678987654;

&qwertyuioppoiuytrewqwer; &oiuytrewqwertyuioiuytrewqwertyuioytrewqwertyuiiuytri;

\t\e\s\t\i\n\g \e\s\c\a\p\e \s\e\q\u\e\n\c\e\s

\!\\\"\#\$\%\&\'\(\)\*\+\,\.\/\:\;\<\=\>\?

\@ \[ \] \^ \_ \` \{ \| \} \~ \- \'

\
\\
\\\
\\\\
\\\\\

\<this\> \<is\> \<not\> \<html\>

Taking commonmark tests from the spec for benchmarking here:

<a><bab><c2c>

<a/><b2/>

<a  /><b2
data="foo" >

<a foo="bar" bam = 'baz <em>"</em>'
_boolean zoop:33=zoop:33 />

<33> <__>

<a h*#ref="hi">

<a href="hi'> <a href=hi'>

< a><
foo><bar/ >

<a href='bar'title=title>

</a>
</foo >

</a href="foo">

foo <!-- this is a
comment - with hyphen -->

foo <!-- not a comment -- two hyphens -->

foo <?php echo $a; ?>

foo <!ELEMENT br EMPTY>

foo <![CDATA[>&<]]>

<a href="&ouml;">

<a href="\*">

<a href="\"">
Valid links:

 [this is a link]()
 [this is a link](<http://something.example.com/foo/bar>)
 [this is a link](http://something.example.com/foo/bar 'test')
 ![this is an image]()
 ![this is an image](<http://something.example.com/foo/bar>)
 ![this is an image](http://something.example.com/foo/bar 'test')
 
 [escape test](<\>\>\>\>\>\>\>\>\>\>\>\>\>\>> '\'\'\'\'\'\'\'\'\'\'\'\'\'\'')
 [escape test \]\]\]\]\]\]\]\]\]\]\]\]\]\]\]\]](\)\)\)\)\)\)\)\)\)\)\)\)\)\))

Invalid links:

 [this is not a link

 [this is not a link](

 [this is not a link](http://something.example.com/foo/bar 'test'
 
 [this is not a link](((((((((((((((((((((((((((((((((((((((((((((((
 
 [this is not a link]((((((((((()))))))))) (((((((((()))))))))))
Valid links:

[[[[[[[[](test)](test)](test)](test)](test)](test)](test)]

[ [[[[[[[[[[[[[[[[[[ [](test) ]]]]]]]]]]]]]]]]]] ](test)

Invalid links:

[[[[[[[[[

[ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [ [

![![![![![![![![![![![![![![![![![![![![![![![![![![![![![![![![![![![![![![

this\
should\
be\
separated\
by\
newlines

this  
should  
be  
separated  
by  
newlines  
too

this
should
not
be
separated
by
newlines

Lorem ipsum dolor sit amet, __consectetur__ adipiscing elit. Cras imperdiet nec erat ac condimentum. Nulla vel rutrum ligula. Sed hendrerit interdum orci a posuere. Vivamus ut velit aliquet, mollis purus eget, iaculis nisl. Proin posuere malesuada ante. Proin auctor orci eros, ac molestie lorem dictum nec. Vestibulum sit amet erat est. Morbi luctus sed elit ac luctus. Proin blandit, enim vitae egestas posuere, neque elit ultricies dui, vel mattis nibh enim ac lorem. Maecenas molestie nisl sit amet velit dictum lobortis. Aliquam erat volutpat.

Vivamus sagittis, diam in [vehicula](https://github.com/markdown-it/markdown-it) lobortis, sapien arcu mattis erat, vel aliquet sem urna et risus. Ut feugiat sapien vitae mi elementum laoreet. Suspendisse potenti. Aliquam erat nisl, aliquam pretium libero aliquet, sagittis eleifend nunc. In hac habitasse platea dictumst. Integer turpis augue, tincidunt dignissim mauris id, rhoncus dapibus purus. Maecenas et enim odio. Nullam massa metus, varius quis vehicula sed, pharetra mollis erat. In quis viverra velit. Vivamus placerat, est nec hendrerit varius, enim dui hendrerit magna, ut pulvinar nibh lorem vel lacus. Mauris a orci iaculis, hendrerit eros sed, gravida leo. In dictum mauris vel augue varius, ac ullamcorper nisl ornare. In eu posuere velit, ac fermentum arcu. Interdum et malesuada fames ac ante ipsum primis in faucibus. Nullam sed malesuada leo, at interdum elit.

Nullam ut tincidunt nunc. [Pellentesque][1] metus lacus, commodo eget justo ut, rutrum varius nunc. Sed non rhoncus risus. Morbi sodales gravida pulvinar. Duis malesuada, odio volutpat elementum vulputate, massa magna scelerisque ante, et accumsan tellus nunc in sem. Donec mattis arcu et velit aliquet, non sagittis justo vestibulum. Suspendisse volutpat felis lectus, nec consequat ipsum mattis id. Donec dapibus vehicula facilisis. In tincidunt mi nisi, nec faucibus tortor euismod nec. Suspendisse ante ligula, aliquet vitae libero eu, vulputate dapibus libero. Sed bibendum, sapien at posuere interdum, libero est sollicitudin magna, ac gravida tellus purus eu ipsum. Proin ut quam arcu.

Suspendisse potenti. Donec ante velit, ornare at augue quis, tristique laoreet sem. Etiam in ipsum elit. Nullam cursus dolor sit amet nulla feugiat tristique. Phasellus ac tellus tincidunt, imperdiet purus eget, ullamcorper ipsum. Cras eu tincidunt sem. Nullam sed dapibus magna. Lorem ipsum dolor sit amet, consectetur adipiscing elit. In id venenatis tortor. In consectetur sollicitudin pharetra. Etiam convallis nisi nunc, et aliquam turpis viverra sit amet. Maecenas faucibus sodales tortor. Suspendisse lobortis mi eu leo viverra volutpat. Pellentesque velit ante, vehicula sodales congue ut, elementum a urna. Cras tempor, ipsum eget luctus rhoncus, arcu ligula fermentum urna, vulputate pharetra enim enim non libero.

Proin diam quam, elementum in eleifend id, elementum et metus. Cras in justo consequat justo semper ultrices. Sed dignissim lectus a ante mollis, nec vulputate ante molestie. Proin in porta nunc. Etiam pulvinar turpis sed velit porttitor, vel adipiscing velit fringilla. Cras ac tellus vitae purus pharetra tincidunt. Sed cursus aliquet aliquet. Cras eleifend commodo malesuada. In turpis turpis, ullamcorper ut tincidunt a, ullamcorper a nunc. Etiam luctus tellus ac dapibus gravida. Ut nec lacus laoreet neque ullamcorper volutpat.

Nunc et leo erat. Aenean mattis ultrices lorem, eget adipiscing dolor ultricies eu. In hac habitasse platea dictumst. Vivamus cursus feugiat sapien quis aliquam. Mauris quam libero, porta vel volutpat ut, blandit a purus. Vivamus vestibulum dui vel tortor molestie, sit amet feugiat sem commodo. Nulla facilisi. Sed molestie arcu eget tellus vestibulum tristique.

[1]: https://github.com/markdown-it

this is a test for tab expansion, be careful not to replace them with spaces

1	4444
22	333
333	22
4444	1


	tab-indented line
    space-indented line
	tab-indented line


a lot of                                                spaces in between here

a lot of												tabs in between here

| scientific\_name           | common\_name        |    n|
|:---------------------------|:--------------------|----:|
| Corvus corone              | Carrion Crow        |  288|
| Turdus merula              | Eurasian Blackbird  |  285|
| Anas platyrhynchos         | Mallard             |  273|
| Fulica atra                | Eurasian Coot       |  268|
| Parus major                | Great Tit           |  266|
| Podiceps cristatus         | Great Crested Grebe |  254|
| Ardea cinerea              | Gray Heron          |  236|
| Cygnus olor                | Mute Swan           |  234|
| Cyanistes caeruleus        | Eurasian Blue Tit   |  233|
| Chroicocephalus ridibundus | Black-headed Gull   |  223|

~~blabla~~
---
title: An example with math elements
---

This is cheap, it only costs 10$!

This example has $\LaTeX$ elements embedded in the
text. It is intended to demonstrate that m $\alpha_\tau$ h
mode can work with tinkr. $y = 
mx + b$

 - [ ] This is an empty checkbox
 - [x] This is a checked checkbox
 - [This is a link](https://ropensci.org)
 - \[this is an example\]

Here is an example from the mathjax website:

When $a \ne 0$, there are two solutions to \(ax^2 + bx + c = 0\) and they are
\[x = {-b \pm \sqrt{b^2-4ac} \over 2a}.\]

```latex
$$
\begin{align} % This mode aligns the equations to the '&=' signs
\begin{split} % This mode groups the equations into one.
\bar{r}_d &= \frac{\sum\sum{cov_{j,k}}}{
                   \sum\sum{\sqrt{var_{j} \cdot var_{k}}}} \\
          &= \frac{V_O - V_E}{2\sum\sum{\sqrt{var_{j} \cdot var_{k}}}}
\end{split}
\end{align}
$$
```

$$
\begin{align} % This mode aligns the equations to the '&=' signs
\begin{split} % This mode groups the equations into one.
\bar{r}_d &= \frac{\sum\sum{cov_{j,k}}}{
                   \sum\sum{\sqrt{var_{j} \cdot var_{k}}}} \\
          &= \frac{V_O - V_E}{2\sum\sum{\sqrt{var_{j} \cdot var_{k}}}}
\end{split}
\end{align}
$$

When $a \ne 0$, there are two solutions to $ax^2 + bx + c = 0$ and they are

```latex
$$
x = {-b \pm \sqrt{b^2-4ac} \over 2a}
$$
```


$$
x = {-b \pm \sqrt{b^2-4ac} \over 2a}
$$

Below is an example from https://github.com/ropensci/tinkr/issues/38
$\frac{\sum _{i=N-n}^{N}Q_i} {\sum_{j=N-n}^{N}{(\frac{C_j+C_{j-1}}2)}}$

```latex
$$
Q_{N(norm)}=\frac{C_N +C_{N-1}}2\times 
\frac{\sum _{i=N-n}^{N}Q_i} {\sum_{j=N-n}^{N}{(\frac{C_j+C_{j-1}}2)}}
$$
```

$$
Q_{N(norm)}=\frac{C_N +C_{N-1}}2\times 
\frac{\sum _{i=N-n}^{N}Q_i} {\sum_{j=N-n}^{N}{(\frac{C_j+C_{j-1}}2)}}
$$
---
title: "What have these birds been studied for? Querying science outputs with R"
slug: birds-science
authors:
  - name: Maëlle Salmon
    url: https://masalmon.eu/
date: 2018-09-11
topicid: 1347
preface: The blog post series corresponds to the material for a talk Maëlle will give at the [Animal Movement Analysis summer school in Radolfzell, Germany on September the 12th](http://animove.org/animove-2019-evening-keynotes/), in a Max Planck Institute of Ornithology.
tags:
- rebird
- birder
- fulltext
- dataone
- EML
- literature
output:
  md_document:
    variant: markdown_github
    preserve_yaml: true
---

In the [second post of the series where we obtained data from
eBird](https://ropensci.org/blog/2018/08/21/birds-radolfzell/) we
determined what birds were observed in the county of Constance, and we
complemented this knowledge with some taxonomic and trait information in
[the fourth post of the
series](https://ropensci.org/blog/2018/09/04/birds-taxo-traits/). Now,
we could be curious about the occurrence of these birds in *scientific
work*. In this post, we will query the scientific literature and an open
scientific data repository for species names: what have these birds been
studied for? Read on if you want to learn how to use R packages allowing
to do so!

### Getting a list of 50 species from occurrence data

For more details about the following code, refer to the [previous post
of the series](https://ropensci.org/blog/2018/08/21/birds-radolfzell/).
The single difference is our adding a step to keep only data for the
most recent years.

``` r
# polygon for filtering
landkreis_konstanz <- osmdata::getbb("Landkreis Konstanz",
                             format_out = "sf_polygon")
crs <- sf::st_crs(landkreis_konstanz)

# get and filter data
f_out_ebd <- "ebird/ebd_lk_konstanz.txt"

library("magrittr")

ebd <- auk::read_ebd(f_out_ebd) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"),
                crs = crs)

in_indices <- sf::st_within(ebd, landkreis_konstanz)

ebd <- dplyr::filter(ebd, lengths(in_indices) > 0)

ebd <- as.data.frame(ebd)

ebd <- dplyr::filter(ebd, approved, lubridate::year(observation_date) > 2010)
```

For the sake of simplicity, we shall only use the 50 species observed
the most often.

``` r
species <- ebd %>%
  dplyr::count(common_name, sort = TRUE) %>%
  head(n = 50) %>%
  dplyr::pull(common_name)
```

The species are Carrion Crow, Eurasian Blackbird, Mallard, Eurasian
Coot, Great Tit, Great Crested Grebe, Mute Swan, Great Cormorant,
Eurasian Blue Tit, Gray Heron, Black-headed Gull, Common Chaffinch,
Common Chiffchaff, Tufted Duck, European Starling, White Wagtail,
European Robin, Little Grebe, Common Wood-Pigeon, Red-crested Pochard,
Ruddy Shelduck, Graylag Goose, Red Kite, Common Buzzard, Eurasian
Blackcap, Great Spotted Woodpecker, Eurasian Magpie, Gadwall, Common
Pochard, Eurasian Nuthatch, Green-winged Teal, House Sparrow, Eurasian
Jay, Yellow-legged Gull, Yellowhammer, Eurasian Green Woodpecker, Eared
Grebe, Eurasian Reed Warbler, Barn Swallow, Northern Shoveler, Eurasian
Moorhen, Black Redstart, Great Egret, White Stork, Eurasian Wren,
Long-tailed Tit, Common House-Martin, Eurasian Kestrel, European
Goldfinch and European Greenfinch
[(`glue::glue_collapse(species, sep = ", ", last = " and ")`)](https://twitter.com/LucyStats/status/1031938964796657665?s=19).

### Querying the scientific literature

Just like rOpenSci has a taxonomic toolbelt
([`taxize`](https://github.com/ropensci/taxize)) and a species
occurrence data toolbelt ([`spocc`](https://github.com/ropensci/spocc)),
it has a super package for querying the scientific literature:
[`fulltext`](https://github.com/ropensci/fulltext)! This package
supports search for “PLOS via the rplos package, Crossref via the
rcrossref package, Entrez via the rentrez package, arXiv via the aRxiv
package, and BMC, Biorxiv, EuroPMC, and Scopus via internal helper
functions”.

We shall use `fulltext` to retrieve the titles and abstracts of
scientific articles mentioning each species, and will use `tidytext` to
compute the most prevalent words in these works.

We first define a function retrieving the titles and abstracts of works
obtained as result when querying one species name.

We use `dplyr::bind_rows` because we want all results for one species at
once, while `fulltext` returns a list of data.frames with one data.frame
by data source.

``` r
.get_papers <- function(species){
  species %>%
    tolower() %>%
    fulltext::ft_search() %>%
    fulltext::ft_get() %>%
    fulltext::ft_collect() %>%
    fulltext::ft_chunks(c("title", "abstract")) %>%
    fulltext::ft_tabularize() %>%
    dplyr::bind_rows()
}

.get_papers(species[1]) %>%
  dplyr::pull(title)
```

    ##  [1] "Great spotted cuckoo nestlings have no antipredatory effect on magpie or carrion crow host nests in southern Spain"
    ##  [2] "Donor-Control of Scavenging Food Webs at the Land-Ocean Interface"
    ##  [3] "Formal comment to Soler et al.: Great spotted cuckoo nestlings have no antipredatory effect on magpie or carrion crow host nests in southern Spain"
    ##  [4] "Socially Driven Consistent Behavioural Differences during Development in Common Ravens and Carrion Crows"
    ##  [5] "Behavioral Responses to Inequity in Reward Distribution and Working Effort in Crows and Ravens"
    ##  [6] "Early Duplication of a Single MHC IIB Locus Prior to the Passerine Radiations"
    ##  [7] "Investigating the impact of media on demand for wildlife: A case study of Harry Potter and the UK trade in owls"
    ##  [8] "New Caledonian Crows Rapidly Solve a Collaborative Problem without Cooperative Cognition"
    ##  [9] "Nest Predation Deviates from Nest Predator Abundance in an Ecologically Trapped Bird"
    ## [10] "Dietary Compositions and Their Seasonal Shifts in Japanese Resident Birds, Estimated from the Analysis of Volunteer Monitoring Data"

If we were working on a scientific study, we’d add a few more filters,
e.g. having the species mentioned in the abstract, and not only
somewhere in the paper which is probably the way the different
literature search providers define a match. But we’re not, so we can
keep our query quite free! My favourite paper involving the Carrion Crow
is [“Investigating the impact of media on demand for wildlife: A case
study of Harry Potter and the UK trade in
owls”](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0182368)
because it’s a fun and important scientific question, and is supported
by open data (by the way you can access CITES trade data (international
trade in endangered species) in R using
[`cites`](https://github.com/ecohealthalliance/cites/) and CITES
Speciesplus database using
[`rcites`](https://ibartomeus.github.io/rcites/)).

We then apply this function to all 50 species and keep each article only
once.

``` r
get_papers <- ratelimitr::limit_rate(.get_papers,
                                     rate = ratelimitr::rate(1, 2))

all_papers <- purrr::map_df(species, get_papers)

nrow(all_papers)
```

    ## [1] 522

``` r
all_papers <- unique(all_papers)

nrow(all_papers)
```

    ## [1] 378

Now, we get the most common words from titles and abstracts. For that we
first append the title to the abstract which is a quick hack.

``` r
library("tidytext")
library("rcorpora")

stopwords <- corpora("words/stopwords/en")$stopWords

all_papers %>%
  dplyr::group_by(title, abstract) %>%
  dplyr::summarise(text = paste(title, abstract)) %>%
  dplyr::ungroup() %>%
  unnest_tokens(word, text) %>%
  dplyr::filter(!word %in% stopwords) %>%
  dplyr::count(word, sort = TRUE) -> words
```

So, what are the most common words in these papers?

``` r
head(words, n = 10)
```

    ##           word   n
    ## 1      species 754
    ## 2        birds 514
    ## 3        virus 270
    ## 4        avian 268
    ## 5         bird 262
    ## 6        study 243
    ## 7     breeding 231
    ## 8         wild 227
    ## 9  populations 217
    ## 10  population 213

Not too surprising, and obviously less entertaining than looking at
individual species’ results. Maybe a wordcloud can give us a better idea
of the wide area of topics of studies involving our 50 most frequent
bird species. We use the [`wordcloud`
package](https://cran.r-project.org/web/packages/wordcloud/index.html).

``` r
library("wordcloud")

with(words, wordcloud(word, n, max.words = 100))
```

![wordcloud of titles and abstracts of scientific
papers](/img/blog-images/2018-09-11-birds-science/wordcloud-1.png)

We see that topics include ecological words such as “foraging” but also
epidemiological questions since “influenza” and “h5n1” come up. Now, how
informative as this wordcloud can be, it’s a bit ugly, so we’ll prettify
it using the [`wordcloud2`
package](https://github.com/Lchiffon/wordcloud2) instead, and the
silhouette of a bird [from
Phylopic](http://phylopic.org/image/6209c9be-060e-4d7f-bc74-a75f3ccf4629/).

``` r
bird <- words %>%
  head(n = 100) %>%
  wordcloud2::wordcloud2(figPath = "bird.png",
                       color = "black", size = 1.5)
# https://www.r-graph-gallery.com/196-the-wordcloud2-library/
htmlwidgets::saveWidget(bird,
                        "tmp.html",
                        selfcontained = F)
```

I wasn’t able to `webshot` the resulting html despite increasing the
`delay` parameter so I screenshot it by hand!

``` r
magick::image_read("screenshot.png")
```

<img src="/img/blog-images/2018-09-11-birds-science/wordcloud2-1.png" alt="wordcloud shaped as a bird" width="1366" />
<p class="caption">
wordcloud shaped as a bird
</p>

The result is a bit kitsch, doesn’t include the word “species”, one
needs to know it’s the silhouette of a bird to recognize it, and we’d
need to work a bit on not reshaping the silhouette, but it’s fun as it
is.

### Querying scientific open data

There are quite a few scientific open data repositories out there, among
which the giant [DataONE](https://www.dataone.org/) that has an API
interfaced with an R package. We shall use it to perform a search
similar to the previous section, but looking at the data indexed on
DataONE. Since DataONE specializes in ecological and environmental data,
we expect to find rather ecological data.

We first define a function to retrieve metadata of datasets for one
species. It looks the species names in the abstract.

``` r
.get_meta <- function(species){

  cn <- dataone::CNode("PROD")
  search <- list(q = glue::glue("abstract:{species}"),
                        fl = "id,title,abstract",
                        sort = "dateUploaded+desc")

  result <- dataone::query(cn, solrQuery = search,
                           as="data.frame")

  if(nrow(result) == 0){
    NULL
  }else{
    # otherwise one line by version
  result <- unique(result)

  tibble::tibble(species = species,
                 title = result$title,
                 abstract = result$abstract)
  }
}
```

Note that DataONE searching could be more precise: one can choose to
search from a given data source only for instance. See the [searching
DataONE
vignette](https://github.com/DataONEorg/rdataone/blob/master/vignettes/searching-dataone.Rmd).

``` r
get_meta <- ratelimitr::limit_rate(.get_meta,
                                     rate = ratelimitr::rate(1, 2))

all_meta <- purrr::map_df(species, get_meta)

nrow(all_meta)
```

    ## [1] 266

``` r
length(unique(all_meta$species))
```

    ## [1] 35

35 species are represented.

``` r
all_meta <- unique(all_meta[,c("title", "abstract")])

nrow(all_meta)
```

    ## [1] 104

We then extract the most common words.

``` r
all_meta %>%
  dplyr::group_by(title, abstract) %>%
  dplyr::summarise(text = paste(title, abstract)) %>%
  dplyr::ungroup() %>%
  unnest_tokens(word, text) %>%
  dplyr::filter(!word %in% stopwords) %>%
  dplyr::count(word, sort = TRUE) -> data_words

head(data_words, n = 10)
```

    ## # A tibble: 10 x 2
    ##    word           n
    ##    <chr>      <int>
    ##  1 data         153
    ##  2 species      120
    ##  3 birds         94
    ##  4 breeding      87
    ##  5 feeding       75
    ##  6 population    65
    ##  7 bird          60
    ##  8 genetic       58
    ##  9 study         56
    ## 10 effects       54

Data is the most common word which is quite logical for metadata of
actual datasets. Let’s also have a look at a regular wordcloud.

``` r
with(data_words, wordcloud(word, n, max.words = 100))
```

![wordcloud of titles and abstracts of scientific
metadata](/img/blog-images/2018-09-11-birds-science/wordcloud3-1.png)

As expected, the words seem more focused on ecology than when looking at
scientific papers. DataONE is a gigantic data catalogue, where one could

-   study the results of such queries (e.g. meta studies of number of,
    say, versions by datasets)

-   or find data to integrate to a new study. If you want to *download*
    data from DataONE, refer to the [download data
    vignette](https://github.com/DataONEorg/rdataone/blob/master/vignettes/download-data.Rmd).

### Conclusion

In this post, we used the rOpenSci `fulltext` package, and the DataONE
`dataone` package, to search for bird species names in scientific papers
and scientific open datasets. We were able to draw wordclouds
representing the diversity of topics of studies in which the birds had
been mentioned or studied. Such a search could be fun to do for your
favourite bird(s)! And in general, following the same approach you could
answer your own specific research question.

#### Scientific literature access

As a reminder, the pipeline to retrieve abstracts and titles of works
mentioning a bird species was quite smooth:

``` r
species %>%
    tolower() %>%
    fulltext::ft_search() %>%
    fulltext::ft_get() %>%
    fulltext::ft_collect() %>%
    fulltext::ft_chunks(c("title", "abstract")) %>%
    fulltext::ft_tabularize() %>%
    dplyr::bind_rows()
```

`fulltext` gives you a lot of power! Other rOpenSci accessing literature
data include [`europepmc`](https://github.com/ropensci/europepmc), R
Interface to Europe PMC RESTful Web Service;
[`jstor`](https://github.com/ropensci/jstor);
[`suppdata`](https://github.com/ropensci/suppdata) for extracting
supplemental information, and [much
more](https://ropensci.org/packages/).

#### Scientific data access… and publication with R

In this post we used the [`dataone`
package](https://github.com/DataONEorg/rdataone) to access data from
DataONE. That same package allows uploading data to DataONE. The
rOpenSci suite features the
[`rfigshare`](https://github.com/ropensci/rfigshare) package for getting
data from, and publishing data to, [Figshare](https://figshare.com/).
For preparing your own data and its documentation for publication, check
out the [`EML` package](https://github.com/ropensci/EML) for writing
metadata respecting the Ecological Metadata Standard, and the [unconf
`dataspice` project](https://github.com/ropenscilabs/dataspice) for
simpler metadata entry.

Explore more of our packages suite, including and beyond access to
scientific literature &data and data publication,
[here](https://ropensci.org/packages/).

#### No more birding? No, your turn!

This was the last post of this series, that hopefully provided an
overview of how rOpenSci packages can help you learn more about birds,
and can support your workflow. As a reminder, in this series we saw

-   [How to identify spots for birding using open geographical
    data](https://ropensci.org/blog/2018/08/14/where-to-bird/).
    Featuring `opencage` for geocoding, `bbox` for bounding box
    creation, `osmdata` for OpenStreetMap’s Overpass API querying,
    `osmplotr` for map drawing using OpenStreetMap’s data.

-   [How to obtain bird occurrence data in
    R](https://ropensci.org/blog/2018/08/21/birds-radolfzell/).
    Featuring `rebird` for interaction with the eBird’s API, and `auk`
    for munging of the whole eBird dataset.

-   [How to extract text from old natural history
    drawings](https://ropensci.org/blog/2018/08/28/birds-ocr/).
    Featuring `magick` for image manipulation, `tesseract` for Optical
    Character Recognition, `cld2` and `cld3` for language detection, and
    `taxize::gnr_resolve` for taxonomic name resolution.

-   [How to complement an occurrence dataset with taxonomy and trait
    information](https://ropensci.org/blog/2018/09/04/birds-taxo-traits/).
    Featuring `taxize`, taxonomic toolbelt for R, and `traits`,
    providing access to species traits data.

-   How to query the scientific literature and scientific open data
    repositories. This is the post you’ve just read!

That’s a wrap! But now, don’t *you* hesitate to explore our packages
suite for your own needs, and to share about your use cases of rOpenSci
packages as a birder or not via [our friendly discussion
forum](https://discuss.ropensci.org/c/usecases)! Happy birding!
---
title: "Untitled"
author: "M. Salmon"
date: "September 6, 2018"
output: html_document
---

```{r setup, include="FALSE", eval="TRUE"}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see [http://rmarkdown.rstudio.com](http://rmarkdown.rstudio.com).

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r, eval="TRUE", echo="TRUE"}
summary(cars)
```

## Including Plots

You can also embed plots, for example:

```{python, fig.cap="pretty plot", echo=-c(1, 2), eval=TRUE}
plot(pressure)
```

```{python}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.

| scientific\_name| common\_name| n| 
 |  --- | --- | --- |
| Corvus corone| Carrion Crow| 288| 
| Turdus merula| Eurasian Blackbird| 285| 
| Anas platyrhynchos| Mallard| 273| 
| Fulica atra| Eurasian Coot| 268| 
| Parus major| Great Tit| 266| 
| Podiceps cristatus| Great Crested Grebe| 254| 
| Ardea cinerea| Gray Heron| 236| 
| Cygnus olor| Mute Swan| 234| 
| Cyanistes caeruleus| Eurasian Blue Tit| 233| 
| Chroicocephalus ridibundus| Black-headed Gull| 223| 

blabla


