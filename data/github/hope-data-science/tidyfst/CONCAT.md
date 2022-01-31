# tidyfst: Tidy Verbs for Fast Data Manipulation<img src="man/figures/hex-tidyfst_url.png" align="right" alt="" width="120" />

 [![](https://www.r-pkg.org/badges/version/tidyfst?color=orange)](https://cran.r-project.org/package=tidyfst) [![](https://img.shields.io/badge/devel%20version-2.0.0-purple.svg)](https://github.com/hope-data-science/tidyfst) ![](https://img.shields.io/badge/lifecycle-stable-green.svg)  [![downloads](http://cranlogs.r-pkg.org/badges/grand-total/tidyfst?color=yellow)](https://r-pkg.org/pkg/tidyfst)

 [![download](https://cranlogs.r-pkg.org/badges/tidyfst?color=red)](https://rdrr.io/cran/tidyfst/) [![downloads](https://cranlogs.r-pkg.org/badges/last-week/tidyfst?color=ff69b4)](https://crantastic.org/packages/tidyfst) [![downloads](https://cranlogs.r-pkg.org/badges/last-day/tidyfst?color=9cf)](https://crantastic.org/packages/tidyfst)

 [![ZENODO DOI](https://zenodo.org/badge/240626994.svg)](https://zenodo.org/badge/latestdoi/240626994) [![JOSS DOI](http://joss.theoj.org/papers/10.21105/joss.02388/status.svg)](https://joss.theoj.org/papers/10.21105/joss.02388)



## Overview

*tidyfst* is a toolkit of tidy data manipulation verbs with *data.table* as the backend . Combining the merits of syntax elegance from *dplyr* and computing performance from *data.table*,  *tidyfst* intends to provide users with state-of-the-art data manipulation tools with least pain. This package is an extension of *data.table*, while enjoying a tidy syntax, it also wraps combinations of efficient functions to facilitate frequently-used data operations.  Also, *tidyfst* would introduce more tidy data verbs from other packages, including but not limited to *tidyverse* and *data.table*. If you are a *dplyr* user but have to use *data.table* for speedy computation,  or *data.table* user looking for readable coding syntax, *tidyfst* is designed for you (and me of course). For further details and tutorials, see [vignettes](https://hope-data-science.github.io/tidyfst/). Both [Chinese](https://hope-data-science.github.io/tidyfst/articles/chinese_tutorial.html) and [English](https://hope-data-science.github.io/tidyfst/articles/english_tutorial.html) tutorials could be found there.

Till now, *tidyfst* has an API that might even transcend its predecessors (e.g. [`select_dt`](https://hope-data-science.github.io/tidyfst/reference/select.html) could accept nearly anything for super column selection). Enjoy the efficient data operations in *tidyfst* !

PS: For extreme performance in tidy syntax, try *tidyfst*'s mirror package [tidyft](https://github.com/hope-data-science/tidyft). 



## Features

- Receives any data.frame (tibble/data.table/data.frame) and returns a data.table.
- Show the variable class of data.table as default.
- Never use in place replacement (also known as modification by reference, which means the original variable would not be modified without notification).
- Use suffix ("_dt") rather than prefix to increase the efficiency (especially when you have IDE with automatic code completion).
- More flexible verbs (e.g. [pairwise_count_dt](https://hope-data-science.github.io/tidyfst/reference/pairwise.html)) for big data manipulation.
- Supporting data importing and parsing with *fst*, which saves both time and memory. Details see [parse_fst/select_fst/filter_fst](https://hope-data-science.github.io/tidyfst/reference/fst.html) and [import_fst/export_fst](https://hope-data-science.github.io/tidyfst/reference/fst_io.html).
- Low and stable dependency on mature packages (data.table, fst, stringr)



## Installation

```R
install.packages("tidyfst")
```



## Example

```R
library(tidyfst)

iris %>%
  mutate_dt(group = Species,sl = Sepal.Length,sw = Sepal.Width) %>%
  select_dt(group,sl,sw) %>%
  filter_dt(sl > 5) %>%
  arrange_dt(group,sl) %>%
  distinct_dt(sl,.keep_all = T) %>%
  summarise_dt(sw = max(sw),by = group)
#>         group  sw
#>        <fctr> <num>
#> 1:     setosa 4.4
#> 2: versicolor 3.4
#> 3:  virginica 3.8

iris %>%
  count_dt(Species) %>%
  add_prop()
#>       Species     n      prop prop_label
#>        <fctr> <int>     <num>     <char>
#> 1:     setosa    50 0.3333333      33.3%
#> 2: versicolor    50 0.3333333      33.3%
#> 3:  virginica    50 0.3333333      33.3%

iris[3:8,] %>%
  mutate_when(Petal.Width == .2,
              one = 1,Sepal.Length=2)
#>    Sepal.Length Sepal.Width Petal.Length Petal.Width Species one
#>          <num>       <num>        <num>       <num>  <fctr> <num>
#> 1:          2.0         3.2          1.3         0.2  setosa   1
#> 2:          2.0         3.1          1.5         0.2  setosa   1
#> 3:          2.0         3.6          1.4         0.2  setosa   1
#> 4:          5.4         3.9          1.7         0.4  setosa  NA
#> 5:          4.6         3.4          1.4         0.3  setosa  NA
#> 6:          2.0         3.4          1.5         0.2  setosa   1


```



## Future plans

*tidyfst* will keep up with the [updates](https://github.com/Rdatatable/data.table/blob/master/NEWS.md) of *data.table* , in the next step would introduce more new features to improve the performance and flexibility to facilitate fast data manipulation in tidy syntax. 



## Vignettes
- [Example 1: Basic usage](https://hope-data-science.github.io/tidyfst/articles/example1_intro.html)
- [Example 2: Join tables](https://hope-data-science.github.io/tidyfst/articles/example2_join.html)
- [Example 3: Reshape](https://hope-data-science.github.io/tidyfst/articles/example3_reshape.html)
- [Example 4: Nest](https://hope-data-science.github.io/tidyfst/articles/example4_nest.html)
- [Example 5: Fst](https://hope-data-science.github.io/tidyfst/articles/example5_fst.html) 
- [Example 6: Dt](https://hope-data-science.github.io/tidyfst/articles/example6_dt.html) 

## Cheat sheet

<a href="docs/tidyfst_cheatsheet.pdf"><img src="tidyfst_cheatsheet.png"/></a>



## Related work

- [data.table](https://github.com/Rdatatable/data.table)
- [fst](https://github.com/fstpackage/fst)
- [tidyr](https://github.com/tidyverse/tidyr)
- [dplyr](https://github.com/tidyverse/dplyr)
- [dtplyr](https://github.com/tidyverse/dtplyr)



## Acknowledgement

The author of [maditr](https://github.com/gdemin/maditr), [Gregory Demin](https://github.com/gdemin) and the author of [fst](https://github.com/fstpackage/fst), [Marcus Klik](https://github.com/MarcusKlik) have helped me a lot in the development of this work. It is so lucky to have them (and many other selfless contributors) in the same open source community of R.


## 1.5.0
1.Set `options("datatable.print.trunc.cols" = TRUE)`, so as to let the printing work like tibbles in dplyr.
2.Make functions in tidyfst could be used in other functions. Details see <https://stackoverflow.com/questions/69098157/how-to-past-parameters-in-r-functions-using-substitute-and-eval-to-make-data>. Some functions have replaced the previous `eval` to `eval.parent`.
3. Export `%like%` from data.table.
4. Add function `sql_join_dt` to implement case insensitive joining for data.frame.
5. Add function `percent` and `add_prop` to calculate percentage conveniently.
6. Add function `pairwise_count_dt` to count pairs of items within a group.

## 1.0.0
Date:20210908
1. Add "fromLast" parameter to `distinct_dt`
2. Add a new function named `col_max` and `col_min` to get the max/min column name
3. Upgrade `dummy_dt` to be faster

## 0.9.9
Date:20200901
1. Do not truncate the columns by default.
2. Add `print_options` to control global printing od data.table.
3. Add citation in the package, linking to the JOSS paper(<https://doi.org/10.21105/joss.02388>)
4. Add `rec_num` and `rec_char` function for variable recoding.
5. Get a cheat sheet for tidyfst.
6. Export `between` from data.table.
7. Support summarisation of multiple functions on multiple columns in `summarise_vars`.

## 0.9.8
Date:20200801
1. Add `rename_with_dt` like dplyr's `rename_with`
2. Update `slice_dt` to support `.N`
3. Update vignette "english_turoial" to remove the outdated codes
4. Improve `count_dt` by using `select_dt` inside
5. Correct error in example of `impute_dt` for user defined functions
6. Export `rleid` and `rleidv` from data.table
7. Add ".name" paramter to `nest_dt` and `squeeze_dt`
8. Debug `slice_max_dt` and `slice_min_dt`
9. Give the slice* family a "by" parameter to slice by group
10. Debug `select_dt`
11. Update the vignette of English tutorial
12. Update `filter_dt` and do not support comma as "&" any more
13. Use testthat package to implement unit test for tidyfst
14. Give sample functions a "by" parameter to sample by group
15. Correct errors in the English tutorial
16. Import data.table v1.13.0 and use its new features


## 0.9.7
Date:20200528
1. Update `separate_dt` to accepte `NA` in parameter "into".
2. Add a new collection of `slice*` function to match dplyr 1.0.0.
3. Simplify the joining functions.
4. Debug `complete_dt` to suppress unnecessary warning in special cases.
5. Debug `nest_dt` to use full join to unnest multiple columns.
6. Debug the joining functions to make it robust for non-data.table data frames.

## 0.9.6
Date:20200502
1. Update Chinense tutorial.
2. Add `impute_dt` to impute missing values using mean, mode and median.
3. Improve `t_dt` to be faster.
4. Add set operations including `union_dt`,etc. This could be used on non-data.table data.frames, which is considered to be convenient.
5. Update "Example 2" vignette.


## 0.9.5
Date: 20200410
0. Reason for update: The update of `as_dt` is very important(see point 5), becasue it is used everywhere in tidyfst. This update might be minor inside the function, but it can improve the performance by large, especially for extremly large data sets (this means in version before 0.9.5[<=0.9.4], operation on large data frames could be quite slow because copies are made in every movement). 
1. Improve `distinct_dt` to receive variables more flexibly.
2. Add `summary_fst` to get info of the fst table.
3. Upgrade "mcols" in `nest_dt` to accept more flexibly by using `select_dt`.
4. Debug `anti_join` and `semi_join` to become more efficient and robust.
5. Update `as_dt` and many functions, which make it faster by reducing data copying when possible, but still stick to principals that never modify by reference. Suppressing the copy when possible, but copies are still made when necessary(using `as.data.table`). 
6. Improve `separate_dt` and `unite_dt`.
7. Improve `replace_dt`.
8. For every `summarise_` and `mutate_`, give a "by" parameter.
9. Add `summarise_when`.

## 0.9.4
Date: 20200402
0. Reason for update: The former introduction of modification by reference is violating the principals of the package, remove them. Modification by reference might be good, I build another package named 'tidyft' to realize it.
1. Add `mat_df` and `df_mat` to covert between named matrix and tidy data.frame, using base-r only.
2. Add `rn_col` and `col_rn`.
3. Add "by" parameter for `summarise_vars` and `mutate_vars`.
4. Make `filter_fst` more robust.
5. Update the vignette of `fst`.
6. Add a new set of join functions with another syntax.
7. Improve `select_fst` with `select_dt`
8. Remove facilities of modification by reference in tidyfst, including `set*` family and "inplace" parameter in `group_by_dt`

## 0.9.3
Date: 20200324
0. Reason for update: The rmarkdown has a poor support of Chinese, which makes the vignette name messy on the CRAN page (see the vignette part of <https://CRAN.R-project.org/package=tidyfst>). Therefore, have to change it to an English name. Also, as many new adjustments coming in, there are some substantial changes for tidyfst to be safer (robust), faster, simpler and feature richer.
1. Improve `group_by_dt` to let it be more flexiable. Now it can receive what `select_dt` receives.
2. Improve `select_fst`, can select one single column by number now.
3. Improve `fill_na_dt` to make it faster with `setnafill`, `shift` and `fcoalesce`.
4. Change the parameter `data` to `.data`. This change of API would be applied to all functions and some other parameters too (start with dot).
5. Remove `drop_all_na_cols` and `drop_all_na_rows`, use `delete_na_cols` and `delete_na_rows` instead to remove columns or rows with NAs larger than a threshold in proportion or number.
6. Rewrite `rename_dt` to be safer.
7. Improve `relocate_dt` to make it faster, by moving names but not data.frame itself, only move at the final step.
8. Remove `mutate_ref`. Design a new family for `set_` to modify by reference. Details see `?set_in_dt`.
9. Add `as_fst` to save a data.frame  as "fst" in tempfile and parse it back in fst_table. 
10. Improve `longer_dt` and `wider_dt` by using `select_mix` to select unchanged columns. Also, change the parameter API to make it more concise. Now it should be easier to use. The vignette of reshape(example 3) is updated too.
11. Make `separate_dt` to be more robust by receiving non-character as column. This means you can use `df %>% separate_dt(x, c("A", "B"))` now. See examples in `?separate_dt`.
12. Give a "by" parameter to `mutate_dt` and `transmute_dt` to mutate by group.
13. Fix a bug in `select_dt`.
14. Remove`all-at-if` collection, use `mutate_vars` and `summarise_vars` instead.
15. Add `replace_dt` to replace any value(s) in data.table.
16. Add an english tutorial and test many basic and complicated examples.
17. Debug `wider_dt` and add a new functionality to take `list` as aggregated function and unchop automatically.
18. Improve `mutate_vars` with raw data.table codes, which is faster.

## 0.8.8
Date: 20200315
0. Reason for update: Check every function in `data.table`, `dplyr` and `tidyr`, optimize and add functionalities when possible, and keep up with the updates of `dplyr` (the upcoming v1.0.0). There are so many substantial updates, so I think an upgrade of version should be proposed. This package is driving to a stable stage later (if no fatal bugs coming after weeks), and the next minor updates will only come after the major updates of data.table (waiting for the release of v1.12.9) and the potential new bugs reported by users. 
1. Get better understanding on non-standard evaluation, update functions that could be optimized. The updated functions include: `mutate_dt`, `transmute_dt`,`arrange_dt`,`distinct_dt`,`slice_dt`,`top_n_dt`,`top_frac_dt`,`mutate_when`. Therefore, now these functions should be faster than before.
2. Add `nth` to extract element of vector via position, useful when we want a single element from the bottom.
3. The API of `longer_dt` has been changed to be more powerful, and update the examples in `wider_dt`. Update the `Example 3: Reshape` vignette.
4. Rewrite the nest part, `nest_by` and `unnest_col` are deprecated, switch to `nest_dt` and `unnest_dt` for new APIs and features. 
5. Design `squeeze_dt` and add `chop_dt`/`unchop_dt` for new usage of nesting.
6. Exporting `frollapply` from data.table, this is a powerful function for aggregation on sliding window.
7. Enhances `select_dt` once more, does not export `select_if_dt` now, merges this functionality directly into `select_dt`. Also, we could now use `-` or `!` to select the negative columns for regular expressions.
8. Optimize `top_n` using `frank` (faster with less memory).
9. Add `sys_time_print` to get the running time more intuitively.
10. Add `uncount_dt`, works just like `tidyr::uncount`.
11. Add `rowwise_dt`, could carry out analysis like `dplyr::rowwise`.
12. Add `relocate_dt` to rearrange columns in data.table.
13. Add `top_dt` and `sample_dt` for convenience.
14. Add `mutate_vars` to complement `all_dt`/`if_dt`/`at_dt`.
15. Add `set_dt` and `mutate_ref` for fast operation by reference of data.table.
16. Add "fun" paramter to `wider_dt` for multiple aggregation.
17. Debug `separate_dt`.
18. Add a Chinese vignette for folks in China (titled as "tidyfst包实例分析").
19. Shorten the description file to be more specific.
20. Add `group_by_dt` and `group_exe_dt` to perform more convenient and efficient group operation.
21. Add `select_mix` for super selection of columns.
22. Fix typos in description.

## 0.7.7
Date: 20200305
0. Reason for update: I've been using `tidyfst` on my daily work by adding `_dt` to many past and current tasks. In these experience, I debug some important functions (they run well on simple tasks, but not on complicated ones), and add more functions. These features are so many that I think an update is necessary for users to get a better tookit earlier. If the update is too frequent, please accept my apology.
1. Optimize `group_dt`. First, it is faster than before because I use `[][]` instead of ` %>% `. (Using `%>%` for `.SD` is slow) Second, I design an alternative to use `.SD` directly in `group_dt`, which might improve the efficiency further.
2. Debug `filter_dt`.
3. Add `fill_na_dt` to fill NAs in data.table. Debug all missing functions. Examples are refreshed.
4. Debug `mutate_when`.
5. Add `complete_dt` to complete a data.frame like `tidyr::complete`.
6. Add `dummy_dt` to get dummy variables from columns.
7. Add `t_dt` to transpose data frame efficiently.
8. Two functions:`as_dt` and `in_dt` to create a short cut to data.table facilities. Add vignette as tutorial in this feature.
9. Add `unite_dt` and `separate_dt` for simple usage.
10. Debug `mutate_dt`.

## 0.6.9
Date: 20200227
0. Reason for urgent update: The use of `show_tibble` violates the principals of programming. I hope this idea would not spread in the vignette. See changes in 4.
1. Improve `select_dt` to let it accept `a:c`-like inputs. Add example `iris %>% select_dt(Sepal.Length:Petal.Length)`. Moreover, now `select_dt` supports delete columns with `-` symbol.
2. Improve `group_dt` to let "by" parameter also accept list of variables, which means we could not use `mtcars %>% group_dt(by =list(vs,am),summarise_dt(avg = mean(mpg)))`.
3. Fix a few typos in description and vignettes.
4. Show the class of variables by default, using `options("datatable.print.class" = TRUE)`, and remove the inappropriate use of `show_tibble`. Details see <https://github.com/tidyverse/tibble/issues/716>.
5. Add `select_if_dt` function. Moreover, support negative conditional selection in `if_dt`.
6. Delete the vignette entitled "Example 5: Tibble", as this feature is not used any more.
7. Add vignette "Example 5:Fst" for better introduction of the feature.
8. Update vignette "Example 1:Basic usage".
 
## 0.6.6 
Date:20200224
1. Change all `print` and `cat` function to `message`.
2. Use `tempdir()` to write file and read it back in the example of `parse_fst`.
3. Fix the bug in `count_dt` and `add_count_dt` and add examples in the function.
4. Add `show_tibble` function, and now the package can use the printing form of tibble to get better information of the data.table. This is not used by default, but might be preferred for tidyverse users.
5. Remove all the unnecessary `\donttest` and use `\dontrun` when have to write files to directory, only to make an example of how to use it(refer to `utils::write.table` document). This should make the best example for real usage.
6. Add URL to Description file.
7. More vignettes added.

## 0.5.7
0. Major updates:(1) Change package name to `tidyfst` (according to the suggestions from CRAN);(2) Do not use `maditr` codes any more (change the description), based on `stringr` and `data.table` only; (3) Support `fst` package with tidy syntax; (4) Add 4 vignettes
1. Support 'fst' package in various ways (see functions end with "_fst")
2. Test the functions and get three vignettes for comparison
3. Totally support group computing with `group_dt` function
4. Correct various typos in the document
5. Rewrite `nest_by` and `unnest_col`. Did not use "_dt" name because they are different from the `tidyverse` API. They might be even more efficient and simple to use.
6. Add "negate" parameter to `select_dt` function.
7. Add `all_dt`,`at_dt` and `if_dt` functions for flexible mutate and summarise.


## 0.3.1
Fix some bugs and add a vignette.

## 0.3.0
Rewrite all functions and use only `data.table` and `stringr` as imported packages.
Have changed the license to MIT.
This time, `tidydt` is lightweight,efficient and powerful. It is totally different from the previous version in many ways.
The previous version would be archived in <https://github.com/hope-data-science/tidydt0>.

## 0.2.1
Some issue seems to happen, check <https://github.com/hope-data-science/tidydt0/issues/1>. Hope to get an offical answer from CRAN. 
Done in the mailing list, keep moving. [20200129]

## v0.2.0 20200123
1. Use new API for `rename_dt`, more like the `rename` in `dplyr`.
2. Change some API name, e.g. `topn_dt` to `top_n_dt`.
3. Add functions to deal with missing values(`replace_na_dt`,`drop_na_dt`).
4. Change the `on_attach.R` file to change the hints.
5. Add `pull_dt`, which I use a lot and so may many others.
6. Add `mutate_when` for another advanced `case_when` utility.
7. Fix according to CRAN suggestions.


---
title: 'tidyfst: Tidy Verbs for Fast Data Manipulation'
tags:

  - R

  - data.table

  - data aggregation

  - data manipulation

  - dplyr
 
  - tidyfst

authors:

  - name: Tian-Yuan Huang
    orcid: 0000-0002-4151-3764
    affiliation: "1" 
  - name: Bin Zhao
    orcid: 0000-0002-3530-2469
    affiliation: "1" 
affiliations:
    
 - name: School of Life Science, Fudan University
   index: 1
date: 3 June 2020
bibliography: paper.bib
---




# Summary

The tidyfst package [@Huang-488] is an R package [@R-Core-Team-479] for fast data manipulation in tidy syntax. The top-level design is inherited from tidy data structure proposed by Hadley Wickham [@Wickham-458], which are: (1) each variable is a column; (2) each observation is a row; (3) each type of observational unit is a table. Moreover, the function names as well as parameter settings are very much borrowed from dplyr [@WickhamFrancois-481] and tidyr [@WickhamHenry-490], reducing the learning cost for tidyverse [@WickhamAverick-493] users. At the bottom, tidyfst is backed by the high performance package data.table [@DowleSrinivasan-480], which is speedy, stable (with little dependency), memory efficient and feature rich.

Sharing similar goals, both data.table and dplyr have gained much popularity among R users. Their features have been compared widely in the community, with pros and cons suggested in ideas and tested in examples (e.g. <https://stackoverflow.com/questions/21435339/data-table-vs-dplyr-can-one-do-something-well-the-other-cant-or-does-poorly/27840349#27840349>). While some opinions might be subjective, consensus could be made on at least two points: (1) data.table could handle data manipulation in less time than dplyr; (2) dplyr has a possibly more user-friendly syntax for learning and communication than data.table. The tidyfst package is designed to combine the merits of dplyr and data.table, so as to provide a suite of tidy verbs for fast data manipulation.

Note that tidyfst is neither the only nor the first package to make trade-offs between data.table and dplyr. Many similar works have been published on CRAN, including dtplyr [@Henry-482], maditr [@Demin-484], table.express [@Sarda-Espinosa-486], tidyfast [@Barrett-487], tidytable [@Fairbanks-491], etc. Nevertheless, tidyfst holds its unique features that no alternative compares so far. One important feature is the support of data manipulation on fst file supported by fst package [@Klik-483]. It means the users could parse the data frames stored in disk first and load the minimum needed subsets to compute on. Other features include convenient column selection in various forms (regular expression, index, etc.), more concise parameter settings and new verbs for frequently-used data operations. 

Furthermore, to save memory and lift speed to a higher level, tidyft [@Huang-489]  has been designed, which utilizes modification by reference feature from data.table whenever possible. The tidyfst and tidyft share similar parameter settings, but function names of tidyft are even simpler (functions in tidyfst usually ends with “_dt”, which tidyft does not). Though tidyft has better performance than tidyfst, it is less robust and demands the users to have deeper understanding on the concepts of modification by reference in data.table. 

Hopefully, tidyfst could provide some reference for the design of dplyr and bring convenience to even data.table users by wrapping some complicated operations in concise steps.



# Acknowledgement

The author of [maditr](https://github.com/gdemin/maditr), [Gregory Demin](https://github.com/gdemin) and the author of [fst](https://github.com/fstpackage/fst), [Marcus Klik](https://github.com/MarcusKlik) have helped us a lot in the development of this work. It is so lucky to have them (and many other selfless contributors) in the same open source community of R.

# References
---
title: "Example 2: Join tables"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{example2_join}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```
  This post has referred to a vignette from `dplyr`, you can find it in <https://dplyr.tidyverse.org/articles/two-table.html>. We'll try to display how to join data tables in this vignette. First, load the packages we need and get some data.
```{r setup}
library(tidyfst)
library(nycflights13)

flights2 <- flights %>% 
  select_dt(year,month,day, hour, origin, dest, tailnum, carrier)

```
  Do a left join with a simple:
```{r}
flights2 %>% 
  left_join_dt(airlines)
```
  
## Controlling how the tables are matched
  Join works the same as `dplyr`:
```{r}
flights2 %>% left_join_dt(weather)
flights2 %>% left_join_dt(planes, by = "tailnum")
flights2 %>% left_join_dt(airports, c("dest" = "faa"))
flights2 %>% left_join_dt(airports, c("origin" = "faa"))
```

## Types of join

```{r}
df1 <- data.table(x = c(1, 2), y = 2:1)
df2 <- data.table(x = c(1, 3), a = 10, b = "a")

df1 %>% inner_join_dt(df2) 
df1 %>% left_join_dt(df2)
df1 %>% right_join_dt(df2)
df1 %>% full_join_dt(df2)

```
  If all you have is a data.frame or tibble, you have no need to change the format. Feed the data directly:
```{r}
df1 <- data.frame(x = c(1, 1, 2), y = 1:3)
df2 <- data.frame(x = c(1, 1, 2), z = c("a", "b", "a"))

df1 %>% left_join_dt(df2)
```
  The "_dt" suffix should remind you that this is backed up by `data.table` and will always return a data.table in the end.

## Filtering joins
  Filtering joins have also been supported in `tidyfst`.
```{r}
flights %>% 
  anti_join_dt(planes, by = "tailnum") %>% 
  count_dt(tailnum, sort = TRUE)
```
  Other examples (`semi_join_dt()` and `anti_join_dt()` never duplicate; they only ever remove observations.):
```{r}
df1 <- data.frame(x = c(1, 1, 3, 4), y = 1:4)
df2 <- data.frame(x = c(1, 1, 2), z = c("a", "b", "a"))

# Four rows to start with:
df1 %>% nrow()

# And we get four rows after the join
df1 %>% inner_join_dt(df2, by = "x") %>% nrow()

# But only two rows actually match
df1 %>% semi_join_dt(df2, by = "x") %>% nrow()
```
  
## Set operations
  For set operations, wrap `data.table`'s function directly, but the functions will automatically turn any data.frame into data.table. Examples are listed as below:
  
```{r}
x = iris[c(2,3,3,4),]
x2 = iris[2:4,]
y = iris[c(3:5),]

intersect_dt(x, y)            # intersect
intersect_dt(x, y, all=TRUE)  # intersect all
setdiff_dt(x, y)              # except
setdiff_dt(x, y, all=TRUE)    # except all
union_dt(x, y)                # union
union_dt(x, y, all=TRUE)      # union all
setequal_dt(x, x2, all=FALSE) # setequal
setequal_dt(x, x2)     
```
  
For more details, just find the help from `data.table` using `?setops`.

---
title: "Use data.table the tidy way: An ultimate tutorial of tidyfst"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{english_tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  #eval = FALSE,
  comment = "#>"
)
```

  I love the tidy syntax of dplyr and the ultimate speed of data.table. Why not take the both? That is why I have started the work of tidyfst, bridging the tidy syntax and computation performance via translation. This tool is especially friendly for dplyr users who want to learn some data.table, but data.table could also benefit from it (more or less).<br><br>
  A great comparison of data.table and dplyr was displayed at <https://atrebas.github.io/post/2019-03-03-datatable-dplyr/> (thanks to Atrebas). I love this tutorial very much because it dig rather deep into many features from both packages. Here I'll try to implement all operations from that tutorial, and the potential users could find why they would prefer tidyfst for some (if not most) tasks.<br><br>
  The below examples have all been checked with tidyfst. Now let's begin.

## Create example data
```{r setup}
library(tidyfst)

set.seed(1L)

## Create a data table
DF <- data.table(V1 = rep(c(1L, 2L), 5)[-10],
                V2 = 1:9,
                V3 = c(0.5, 1.0, 1.5),
                V4 = rep(LETTERS[1:3], 3))
copy(DF) -> DT

class(DF)
DF
```

# Basic operations
## Filter rows
```{r}
### Filter rows using indices
slice_dt(DF, 3:4)

### Discard rows using negative indices
slice_dt(DF, -(3:7))

### Filter rows using a logical expression
filter_dt(DF, V2 > 5)
filter_dt(DF, V4 %in% c("A", "C"))
filter_dt(DF, V4 %chin% c("A", "C")) # fast %in% for character

### Filter rows using multiple conditions
filter_dt(DF, V1 == 1, V4 == "A")
# equals to
filter_dt(DF, V1 == 1 & V4 == "A")

### Filter unique rows
distinct_dt(DF) # unique(DF)
distinct_dt(DF, V1,V4)

### Discard rows with missing values
drop_na_dt(DF) # na.omit(DF)

### Other filters
sample_n_dt(DF, 3)      # n random rows
sample_frac_dt(DF, 0.5) # fraction of random rows
top_n_dt(DF, 1, V1)     # top n entries (includes equals)

filter_dt(DT,V4 %like% "^B")
filter_dt(DT,V2 %between% c(3, 5))
filter_dt(DT,between(V2, 3, 5, incbounds = FALSE))
filter_dt(DT,V2 %inrange% list(-1:1, 1:3))  # see also ?inrange
```

## Sort rows
```{r}
### Sort rows by column
arrange_dt(DF, V3)

### Sort rows in decreasing order
arrange_dt(DF, -V3)

### Sort rows based on several columns
arrange_dt(DF, V1, -V2)
```

## Select columns
```{r}
### Select one column using an index (not recommended)
pull_dt(DT,3) # returns a vector
select_dt(DT,3) # returns a data.table

### Select one column using column name
select_dt(DF, V2) # returns a data.table
pull_dt(DF, V2)   # returns a vector

### Select several columns
select_dt(DF, V2, V3, V4)
select_dt(DF, V2:V4) # select columns between V2 and V4

### Exclude columns
select_dt(DF, -V2, -V3)

### Select/Exclude columns using a character vector
cols <- c("V2", "V3")
select_dt(DF,cols = cols)
select_dt(DF,cols = cols,negate = TRUE)

### Other selections
select_dt(DF, cols = paste0("V", 1:2))
relocate_dt(DF, V4) # reorder columns
select_dt(DF, "V")
select_dt(DF, "3$")
select_dt(DF, ".2")
select_dt(DF, "V1")
select_dt(DF, -"^V2")
# remove variables using "-" prior to function
```

## Summarise data
```{r}
### Summarise one column
summarise_dt(DF, sum(V1)) # returns a data.table
summarise_dt(DF, sumV1 = sum(V1)) # returns a data.table

### Summarise several columns
summarise_dt(DF, sum(V1), sd(V3))

### Summarise several columns and assign column names
DF %>%
  summarise_dt(sumv1 = sum(V1),
              sdv3  = sd(V3))

### Summarise a subset of rows
DT[1:4, sum(V1)]
DF %>%
  slice_dt(1:4) %>%
  summarise_dt(sum(V1))

### Misc
summarise_dt(DF, nth(V3,1))
summarise_dt(DF, nth(V3,-1))
summarise_dt(DF, nth(V3, 5))
summarise_dt(DF, uniqueN(V4))
uniqueN(DF)
```


## group computation (by)
```{r}
### By group
# not recommended
DF %>%
  group_dt(
    by = V4,
    summarise_dt(sumV2 = sum(V2))
  )
# recommended
DF %>% 
  summarise_dt(sumV2 = sum(V2),by = V4)


### By several groups
DF %>% 
  summarise_dt(sumV2 = sum(V2),by = .(V1,V4))

### Calling function in by
DF %>% 
  summarise_dt(sumV2 = sum(V2),by = tolower(V4))


### Assigning column name in by
DF %>% 
  summarise_dt(sumV2 = sum(V2),by = .(abc = tolower(V4)))


### Using a condition in by
DF %>% 
  summarise_dt(sumV2 = sum(V2),by = V4 == "A")



### By on a subset of rows
DF %>% 
  slice_dt(1:5) %>% 
  summarise_dt(sumV1 = sum(V1),by = V4)

### Count number of observations for each group
count_dt(DF, V4)

### Add a column with number of observations for each group
add_count_dt(DF, V1)

### Retrieve the first/last/nth observation for each group
DF %>% summarise_dt(by = V4,nth(V2,1))
DF %>% summarise_dt(by = V4,nth(V2,-1))
DF %>% summarise_dt(by = V4,nth(V2,2))

```

# Going further

## Advanced columns manipulation

```{r}
### Summarise all the columns
summarise_vars(DT,.func = max)

### Summarise several columns
summarise_vars(DT,c("V1", "V2"),mean)

### Summarise several columns by group
DT %>% 
    summarise_vars(c("V1", "V2"),mean,by = V4)

## using patterns (regex)
DT %>% 
    summarise_vars("V1|V2",mean,by = V4)

## Summarise with more than one function by group
# when you can't find a way, you can always use `in_dt` to use data.table
DT %>% 
  in_dt(, by = V4, 
     c(lapply(.SD, sum),
       lapply(.SD, mean)))

### Summarise using a condition
summarise_vars(DF, is.numeric, mean)

### Modify all the columns
mutate_vars(DF, .func = rev)

### Modify several columns (dropping the others)
DF %>% 
  select_dt(cols = c("V1", "V2")) %>% 
  mutate_vars(.func = sqrt)
DF %>% 
  select_dt(-V4) %>% 
  mutate_vars(.func = exp)

### Modify several columns (keeping the others)
DF %>% 
  mutate_vars(c("V1", "V2"), sqrt)
DF %>% 
  mutate_vars(-"V4",  exp)

### Modify columns using a condition (dropping the others)
select_dt(DT,is.numeric)

### Modify columns using a condition (keeping the others)
mutate_vars(DT,is.numeric,as.integer)

### Use a complex expression
DF %>% 
  group_dt(
    by = V4,
    slice_dt(1:2) %>% 
      transmute_dt(V1 = V1,
            V2 = "X")
  )

### Use multiple expressions (with DT[,{j}])
DT %>% 
  in_dt(,{
      print(V1) #  comments here!
      print(summary(V1))
      x <- V1 + sum(V2)
     .(A = 1:.N, B = x) # last list returned as a data.table
     }
  )

```

## Advanced use of by

```{r}
### Select first/last/… row by group
DT %>% 
  group_dt(
    by = V4,
    head(1)
  )

DT %>% 
  group_dt(
    by = V4,
    tail(2)
  )

DT %>% 
  group_dt(
    by = V4,
    slice_dt(1,.N)
  )

### Select rows using a nested query
DF %>% 
  group_dt(
    by = V4,
    arrange_dt(V2) %>% 
      slice_dt(1)
  )

### Add a group counter column
DT %>% 
  mutate_dt(Grp = .GRP,by = .(V4, V1))

### Get row number of first (and last) observation by group
DT %>% summarise_dt(I = .I,by = V4)
DT %>% summarise_dt(I = .I[1],by = V4)
DT %>% summarise_dt(I = .I[c(1,.N)],by = V4)

### Handle list-columns by group

DT %>% 
  select_dt(V1,V4) %>% 
  chop_dt(V1) # return V1 as a list
DT %>% nest_dt(V4) # subsets of the data

### Grouping sets (multiple by at once)
# use data.table directly, tidyfst does not provide new methods for it yet
data.table::rollup(DT,
       .(SumV2 = sum(V2)),
       by = c("V1", "V4"))

data.table::rollup(DT,
       .(SumV2 = sum(V2), .N),
       by = c("V1", "V4"),
       id = TRUE)

data.table::cube(DT,
     .(SumV2 = sum(V2), .N),
     by = c("V1", "V4"),
     id = TRUE)

data.table::groupingsets(DT,
             .(SumV2 = sum(V2), .N),
             by   = c("V1", "V4"),
             sets = list("V1", c("V1", "V4")),
             id   = TRUE)

```

# Miscellaneous

## Read / Write data
tidyfst exports `data.table::fread` and `data.table::fwrite` directly.

```{r,eval=FALSE}
### Write data to a csv file
fwrite(DT, "DT.csv")

### Write data to a tab-delimited file
fwrite(DT, "DT.txt", sep = "\t")

### Write list-column data to a csv file
fwrite(setDT(list(0, list(1:5))), "DT2.csv")
#
### Read a csv / tab-delimited file
fread("DT.csv")
# fread("DT.csv", verbose = TRUE) # full details
fread("DT.txt", sep = "\t")

### Read a csv file selecting / droping columns
fread("DT.csv", select = c("V1", "V4"))
fread("DT.csv", drop = "V4")
# NA
### Read and rbind several files
rbindlist(lapply(c("DT.csv", "DT.csv"), fread))
# c("DT.csv", "DT.csv") %>% lapply(fread) %>% rbindlist

```

## Reshape data

```{r}
### Melt data (from wide to long)
fsetequal(DT,DF)

mDT = DT %>% longer_dt(V3,V4)
mDF = DF %>% longer_dt(-"V1|V2")

fsetequal(mDT,mDF)
mDT

### Cast data (from long to wide)
mDT %>% 
  wider_dt(V4,name = "name",value = "value")
# below is a special case and could only be done in tidyfst
mDT %>% 
  wider_dt(V4,name = "name",value = "value",fun = list)

mDT %>% 
  wider_dt(V4,name = "name",value = "value",fun = sum)

### Split
split(DT, by = "V4")

```

## Other
```{r}
### Lead/Lag
lag_dt(1:10,n = 1)
lag_dt(1:10,n = 1:2)
lead_dt(1:10,n = 1)

```

# Join/Bind data sets

## Join
```{r}
x <- data.table(Id  = c("A", "B", "C", "C"),
                X1  = c(1L, 3L, 5L, 7L),
                XY  = c("x2", "x4", "x6", "x8"),
                key = "Id")

y <- data.table(Id  = c("A", "B", "B", "D"),
                Y1  = c(1L, 3L, 5L, 7L),
                XY  = c("y1", "y3", "y5", "y7"),
                key = "Id")

### left join
left_join_dt(x, y, by = "Id")

### right join
right_join_dt(x, y, by = "Id")

### inner join
inner_join_dt(x, y, by = "Id")

### full join
full_join_dt(x, y, by = "Id")

### semi join
semi_join_dt(x, y, by = "Id")

### anti join
anti_join_dt(x, y, by = "Id")

```

## Bind
```{r}
x <- data.table(1:3)
y <- data.table(4:6)
z <- data.table(7:9, 0L)
### Bind rows
rbind(x, y)
rbind(x, z, fill = TRUE)

### Bind rows using a list
rbindlist(list(x, y), idcol = TRUE)

### Bind columns
cbind(x, y)

```

## Set operations
```{r}
x <- data.table(c(1, 2, 2, 3, 3))
y <- data.table(c(2, 2, 3, 4, 4))
### Intersection
fintersect(x, y)
fintersect(x, y, all = TRUE)

### Difference
fsetdiff(x, y)
fsetdiff(x, y, all = TRUE)

### Union
funion(x, y)
funion(x, y, all = TRUE)

### Equality
fsetequal(x, x[order(-V1),])
all.equal(x, x) # S3 method
setequal(x, x[order(-V1),])

```

# Summary
To break all these codes through, tidyfst has improved bit by bit. If you are using tidyfst frequently, you'll find that while it enjoys a tidy syntax, it is more like you are using data.table in another style. Compared with many other packages with similar goals, tidyfst sticks to many principles of data.table (and is more like data.table in many ways). However, the ultimate goal is still clear: providing users with state-of-the-art data manipulation tools with least pain. Therefore, keep it simple and make it fast. Enjoy tidyfst~



---
title: "Performance"
output: rmarkdown::html_vignette
author: Tian-Yuan Huang (huang.tian-yuan@qq.com)
vignette: >
  %\VignetteIndexEntry{Performance}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  comment = "#>"
)
```
  One may wonder how fast is *tidyfst*. Well, it depends. Generally, it is as fast as *data.table* because it is backed by it, but it would spend extra time on the generation of data.table codes. This extra time is marginal on large (and even small) data sets.<br><br>
  Now let's do a test to compare the performance of *tidyfst*, *data.table* and *dplyr*. In the vignette we'll use a small data set. The example was provided by the *data.table* package (<https://h2oai.github.io/db-benchmark/>) and tweaked here. These tests are based on computation by groups.<br><br>
  First let's load the package and generate some data.<br><br>
```{r setup}
# load packages
library(tidyfst)
library(data.table)
library(dplyr)
library(bench)

# generate the data
# if you have a HPC and want to try larger data sets, increase N
N = 1e4 
K = 1e2

set.seed(2020)

cat(sprintf("Producing data of %s rows and %s K groups factors\n", N, K))

DT = data.table(

  id1 = sample(sprintf("id%03d",1:K), N, TRUE),      # large groups (char)

  id2 = sample(sprintf("id%03d",1:K), N, TRUE),      # large groups (char)

  id3 = sample(sprintf("id%010d",1:(N/K)), N, TRUE), # small groups (char)

  id4 = sample(K, N, TRUE),                          # large groups (int)

  id5 = sample(K, N, TRUE),                          # large groups (int)

  id6 = sample(N/K, N, TRUE),                        # small groups (int)

  v1 =  sample(5, N, TRUE),                          # int in range [1,5]

  v2 =  sample(5, N, TRUE),                          # int in range [1,5]

  v3 =  round(runif(N,max=100),4)                    # numeric e.g. 23.5749

)

object_size(DT)
```
  This data is rather small, the size is around 527 Kb. However, with the *bench* package, we could detect the difference by increasing iteration times. In this way, examples listed here could be implemented even on relatively low performance computers. 
  
## Q1 
  Here, we try to get median and standard deviation by groups.After dplyr v1.0.0, the regrouping feature could be confusing sometimes (comes with warning message). If you are using it, make sure they are in the right groups before grouped computation. In tidyfst and data.table, we have "by" parameter to specify the groups. Here we would not check if the results are equal, because dplyr will return a tibble class even when we input a data.table in the first place. The iteration time is 10 for each of the test below.
```{r}
bench::mark(
  data.table = DT[,.(median_v3 = median(v3),
                     sd_v3 = sd(v3)),
                  by = .(id4,id5)],
  tidyfst = DT %>%
    summarise_dt(
      by = "id4,id5",
      median_v3 = median(v3),
      sd_v3 = sd(v3)
    ),
  dplyr = DT %>%
    group_by(id4,id5,.drop = TRUE) %>%
    summarise(median_v3 = median(v3),sd_v3 = sd(v3)),
  check = FALSE,iterations = 10
) -> q1

q1
```
  We could find that spent time of tidyfst and data.table are quite similar, but much less than dplyr.

## Q2 
  This example performs quite similar to the above one. *tidyfst* might spend a tiny little more time and space on code translation than *data.table*, but still performs much better than *dplyr*.
```{r}
bench::mark(
  data.table =DT[,.(range_v1_v2 = max(v1) - min(v2)),by = id3],
  tidyfst = DT %>% summarise_dt(
    by = id3,
    range_v1_v2 = max(v1) - min(v2)
  ),
  dplyr = DT %>%
    group_by(id3,.drop = TRUE) %>%
    summarise(range_v1_v2 = max(v1) - min(v2)),
  check = FALSE,iterations = 10
) -> q2

q2
```

## Q3
  Here we'll display a rather different test to show the flexibly in *tidyfst*. In *tidyfst*, if your code writes more like *data.table*, the codes could speed up. If you write it more like *dplyr*, the codes might be more readable but slows down. In *tidyfst*, there is `in_dt` function for you to write data.table codes to gain speed when you meet a bottomneck.<br><br>
  In the following example, we use the exact same syntax of *data.table* in `tidyfst::in_dt`. 
```{r}
bench::mark(
  data.table =DT[order(-v3),.(largest2_v3 = head(v3,2L)),by = id6],
  tidyfst = DT %>%
    in_dt(order(-v3),.(largest2_v3 = head(v3,2L)),by = id6),
  dplyr = DT %>%
    select(id6,largest2_v3 = v3) %>%
    group_by(id6) %>%
    slice_max(largest2_v3,n = 2,with_ties = FALSE),
  check = FALSE,iterations = 10
) -> q3

q3
```
## Q4
  To summarise multiple columns by group, *tidyfst* has designed a function named `summarise_vars`, which is even more convenient than the `across` function in *dplyr*. It first choose the columns, then tell it what to do, and you can provide the "by" parameter to operate by groups (optional).
```{r}
bench::mark(
  data.table =DT[,lapply(.SD,mean),by = id4,.SDcols = v1:v3],
  tidyfst = DT %>%
    summarise_vars(
      v1:v3,
      mean,
      by = id4
    ),
  dplyr = DT %>%
    group_by(id4) %>%
    summarise(across(v1:v3,mean)),
  check = FALSE,iterations = 10
) -> q4

q4
```
   Take a look at the performance, *tidyfst* still lies between *data.table* and *dplyr*.
   
## Q5
  Now let's try more groups, here we use all the id (id1~id6) as group, and get the sum and count. Note that *tidyfst* is written in *data.table*, so it do not use `n()` in *dplyr* but `.N` in *data.table* to get counts by group.

```{r}
bench::mark(
  data.table =DT[,.(v3 = sum(v3),count = .N),by = id1:id6],
  tidyfst = DT %>%
    summarise_dt(
      by = id1:id6,
      v3 = sum(v3),
      count = .N
    ),
  dplyr = DT %>%
    group_by(id1,id2,id3,id4,id5,id6) %>%
    summarise(v3 = sum(v3),count = n()),
  check = FALSE,iterations = 10
) -> q5

q5
```

## Last words
  While in a data set of ~0.5 Mb we find that the performance of *tidyfst* lies between *data.table* and *dplyr*, we could discover that the speed is much closer to *data.table*. In fact, if you try a much larger data set in a computer with large RAM and multiple cores, you'll find that the performance of *tidyfst* sticks close to *data.table*. If you are interested and has a high-performance computer, try to generate a larger data set and test out.  Moreover, while the dplyr user might find these data manipulation verbs friendly, the innate syntax of *tidyfst* is more like *data.table*, and could be a good companion of *data.table* for some frequently used complex tasks.

## Session information
```{r}
sessionInfo()
```







---
title: "Example 1: Basic usage"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{example1_intro}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)
```

# Use tidyfst just like dplyr
This part of vignette has referred to `dplyr`'s vignette in <https://dplyr.tidyverse.org/articles/dplyr.html>. We'll try to reproduce all the results. First load the needed packages.

```{r}
library(tidyfst)
library(nycflights13)
library(data.table)

data.table(flights)
```

## Filter rows with `filter_dt()`
```{r}
filter_dt(flights, month == 1 & day == 1)
```
Note that comma could not be used in the expressions. Which means `filter_dt(flights, month == 1,day == 1)` would return error.
## Arrange rows with `arrange_dt()`
```{r}
arrange_dt(flights, year, month, day)
```

  Use `-` (minus symbol) to order a column in descending order:
```{r}
arrange_dt(flights, -arr_delay)
```

## Select columns with `select_dt()`

```{r}
select_dt(flights, year, month, day)
```

  `select_dt(flights, year:day)` and `select_dt(flights, -(year:day))` are not supported. But I have added a feature to help select with regular expression, which means you can:
```{r}
select_dt(flights, "^dep")
```
  The rename process is almost the same as that in `dplyr`:
```{r}
select_dt(flights, tail_num = tailnum)
rename_dt(flights, tail_num = tailnum)
```
  
## Add new columns with `mutate_dt()`
```{r}
mutate_dt(flights,
  gain = arr_delay - dep_delay,
  speed = distance / air_time * 60
)
```
  
  However, if you just create the column, please split them. The following codes would not work:
```{r,eval=FALSE}
mutate_dt(flights,
  gain = arr_delay - dep_delay,
  gain_per_hour = gain / (air_time / 60)
)
```
  Instead, use:
```{r}
mutate_dt(flights,gain = arr_delay - dep_delay) %>%
  mutate_dt(gain_per_hour = gain / (air_time / 60))
```
  
  If you only want to keep the new variables, use `transmute_dt()`:
```{r}
transmute_dt(flights,
  gain = arr_delay - dep_delay
)
```
  
## Summarise values with `summarise_dt()`
```{r}
summarise_dt(flights,
  delay = mean(dep_delay, na.rm = TRUE)
)
```

## Randomly sample rows with `sample_n_dt()` and `sample_frac_dt()`
```{r}
sample_n_dt(flights, 10)
sample_frac_dt(flights, 0.01)
```

## Grouped operations
  For the below `dplyr` codes:
```{r,eval=FALSE}
by_tailnum <- group_by(flights, tailnum)
delay <- summarise(by_tailnum,
  count = n(),
  dist = mean(distance, na.rm = TRUE),
  delay = mean(arr_delay, na.rm = TRUE))
delay <- filter(delay, count > 20, dist < 2000)
```
  We could get it via:
```{r}
flights %>% 
  summarise_dt( count = .N,
  dist = mean(distance, na.rm = TRUE),
  delay = mean(arr_delay, na.rm = TRUE),by = tailnum)
```
  `summarise_dt` (or `summarize_dt`) has a parameter "by", you can specify the group.
  We could find the number of planes and the number of flights that go to each possible destination:
```{r}
# the dplyr syntax:
# destinations <- group_by(flights, dest)
# summarise(destinations,
#   planes = n_distinct(tailnum),
#   flights = n()
# )

summarise_dt(flights,planes = uniqueN(tailnum),flights = .N,by = dest) %>% 
  arrange_dt(dest)

```
  If you need to group by many variables, use:
```{r}
# the dplyr syntax:
# daily <- group_by(flights, year, month, day)
# (per_day   <- summarise(daily, flights = n()))

flights %>% 
  summarise_dt(by = .(year,month,day),flights = .N)

# (per_month <- summarise(per_day, flights = sum(flights)))
flights %>% 
  summarise_dt(by = .(year,month,day),flights = .N) %>% 
  summarise_dt(by = .(year,month),flights = sum(flights))

# (per_year  <- summarise(per_month, flights = sum(flights)))
flights %>% 
  summarise_dt(by = .(year,month,day),flights = .N) %>% 
  summarise_dt(by = .(year,month),flights = sum(flights)) %>% 
  summarise_dt(by = .(year),flights = sum(flights))
```
  
# Comparison with data.table syntax
  *tidyfst* provides a tidy syntax for *data.table*. For such design, *tidyfst* never runs faster than the analogous *data.table* codes. Nevertheless, it facilitate the dplyr-users to gain the computation performance in no time and guide them to learn more about data.table for speed.
  Below, we'll compare the syntax of `tidyfst` and `data.table` (referring to  [Introduction to data.table](https://rdatatable.gitlab.io/data.table/articles/datatable-intro.html)). This could let you know how they are different, and let users to choose their preference. Ideally, *tidyfst* will lead even more users to learn more about *data.table* and its wonderful features, so as to design more extentions for *tidyfst* in the future.
  
## Data
Because we want a more stable data source, here we'll use the flight data from the above `nycflights13` package.
```{r}
library(tidyfst)
library(data.table)
library(nycflights13)

flights = data.table(flights) %>% na.omit()
```

## Subset rows
```{r}
# data.table
head(flights[origin == "JFK" & month == 6L])
flights[1:2]
flights[order(origin, -dest)] 

# tidyfst
flights %>% 
  filter_dt(origin == "JFK" & month == 6L) %>% 
  head()
flights %>% slice_dt(1:2)
flights %>% arrange_dt(origin,-dest)
```

## Select column(s)
```{r}
# data.table
flights[, list(arr_delay)]
flights[, .(arr_delay, dep_delay)]
flights[, .(delay_arr = arr_delay, delay_dep = dep_delay)]

# tidyfst
flights %>% select_dt(arr_delay)
flights %>% select_dt(arr_delay, dep_delay)
flights %>% transmute_dt(delay_arr = arr_delay, delay_dep = dep_delay)
```

## Mixed computation
```{r}
# data.table
flights[, sum( (arr_delay + dep_delay) < 0)]
flights[origin == "JFK" & month == 6L,
               .(m_arr = mean(arr_delay), m_dep = mean(dep_delay))]
flights[origin == "JFK" & month == 6L, length(dest)]
flights[origin == "JFK" & month == 6L, .N]

# tidyfst
flights %>% summarise_dt(sum( (arr_delay + dep_delay) < 0))
flights %>% 
  filter_dt(origin == "JFK" & month == 6L) %>% 
  summarise_dt(m_arr = mean(arr_delay), m_dep = mean(dep_delay))
flights %>% 
  filter_dt(origin == "JFK" & month == 6L) %>% 
  nrow()
flights %>% 
  filter_dt(origin == "JFK" & month == 6L) %>% 
  count_dt()
flights %>% 
  filter_dt(origin == "JFK" & month == 6L) %>% 
  summarise_dt(.N)
```
  In the above examples, we could learn that in *tidyfst*, you could still use the methods in data.table, such as `.N`.
  
## Refer to columns by names
```{r}
# data.table
flights[, c("arr_delay", "dep_delay")]

select_cols = c("arr_delay", "dep_delay")
flights[ , ..select_cols]
flights[ , select_cols, with = FALSE]

flights[, !c("arr_delay", "dep_delay")]
flights[, -c("arr_delay", "dep_delay")]

# returns year,month and day
flights[, year:day]
# returns day, month and year
flights[, day:year]
# returns all columns except year, month and day
flights[, -(year:day)]
flights[, !(year:day)]

# tidyfst
flights %>% select_dt(c("arr_delay", "dep_delay"))

select_cols = c("arr_delay", "dep_delay")
flights %>% select_dt(cols = select_cols)

flights %>% select_dt(-arr_delay,-dep_delay)

flights %>% select_dt(year:day)
flights %>% select_dt(day:year)
flights %>% select_dt(-(year:day))
flights %>% select_dt(!(year:day))
```

## Aggregations
```{r}
# data.table
flights[, .N, by = .(origin)]
flights[carrier == "AA", .N, by = origin]
flights[carrier == "AA", .N, by = .(origin, dest)]
flights[carrier == "AA",
        .(mean(arr_delay), mean(dep_delay)),
        by = .(origin, dest, month)]

# tidyfst
flights %>% count_dt(origin) # sort by default
flights %>% filter_dt(carrier == "AA") %>% count_dt(origin)
flights %>% filter_dt(carrier == "AA") %>% count_dt(origin,dest)
flights %>% filter_dt(carrier == "AA") %>% 
  summarise_dt(mean(arr_delay), mean(dep_delay),
               by = .(origin, dest, month))
```
  Note that currently `keyby` is not used in *tidyfst*. This featuer might be included in the future for better performance in order-independent tasks. Moreover, `count_dt` is sorted automatically by the counted number, this could be controlled by the parameter "sort".
```{r}
# data.table
flights[carrier == "AA", .N, by = .(origin, dest)][order(origin, -dest)]
flights[, .N, .(dep_delay>0, arr_delay>0)]

# tidyfst
flights %>% 
  filter_dt(carrier == "AA") %>% 
  count_dt(origin,dest,sort = FALSE) %>% 
  arrange_dt(origin,-dest)
flights %>% 
  summarise_dt(.N,by = .(dep_delay>0, arr_delay>0))
```
  Now let's try a more complex example:
```{r}
# data.table
flights[carrier == "AA", 
        lapply(.SD, mean), 
        by = .(origin, dest, month), 
        .SDcols = c("arr_delay", "dep_delay")] 

# tidyfst
flights %>% 
  filter_dt(carrier == "AA") %>% 
  group_dt(
    by = .(origin, dest, month),
    at_dt("_delay",summarise_dt,mean)
           )
```
  Let me explain what happens here, especially in `group_dt`. First filter by condition `carrier == "AA"`, then group by three variables, which are `origin, dest, month`. Last, summarise by columns with "_delay" in the column names and get the mean value of all such variables(with "_delay" in their column names). This is a very creative design, utilizing `.SD` in *data.table* and upgrade the `group_by` function in *dplyr* (because you never need to `ungroup` now, just put the group operations in the `group_dt`). And **you can pipe in the group_dt function**. Let's play with it a little bit further:
```{r}
flights %>% 
  filter_dt(carrier == "AA") %>% 
  group_dt(
    by = .(origin, dest, month),
    at_dt("_delay",summarise_dt,mean) %>% 
      mutate_dt(sum = dep_delay + arr_delay)
           )
```
  However, I don't recommend using it if you don't acutually need it for group computation (just start another pipe follows `group_dt`).
  Now let's end with some easy examples:
```{r}
# data.table
flights[, head(.SD, 2), by = month]

# tidyfst
flights %>% 
  group_dt(by = month,head(2))
```
  Deep inside, *tidyfst* is born from *dplyr* and *data.table*, and use *stringr* to make flexible APIs, so as to bring their superiority into full play.
  
  
  
  


---
title: "Example 6: Dt"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{example6_dt}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```
  For absolute physical speed, use *data.table* directly. While the learning curve might be longer, the improvement of computation performance pays off if you are dealing with large datasets frequently. There are several ways to cut into *data.table* syntax to gain higher performance in *tidyfst*. A convenient way is to use the `DT[I,J,BY]` syntax after the pipe(`%>%`).
```{r setup}
library(tidyfst)
iris %>% 
  as_dt()%>%  #coerce a data.frame to data.table
  .[,.SD[1],by = Species]

```
  This syntax is not so consistent with the tidy syntax, therefore `in_dt` is also designed for the short cut to *data.table* method, which could be used as:
```{r}
iris %>% 
  in_dt(,.SD[1],by = Species)
```
  
  `in_dt` follows the basic principals of *tidyfst*, which include: (1) Never use in place replacement. Therefore, the in place functions like `:=` will still return the results. (2) Always recieves a data frame (data.frame/tibble/data.table) and returns a data.table. This means you don't have to write `as.data.table` or `as_dt` all the time as long as you are working on data frames in R.
 
  
---
title: "Example 5: Fst"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{example5_fst}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)
```
 
 The [fst](https://github.com/fstpackage/fst) package for R provides a fast, easy and flexible way to serialize data frames. It has very amazing features, such as fast read and write of R data frames, super file compression and parse data frames without reading it. Considering all these features, now [tidyfst](https://github.com/hope-data-science/tidyfst) could provide a new workflow to manipulate data more efficiently. The core idea is: We never need the whole data all at once, we only need the things we want and aggregate them to get the summary to provide target information. 
 
 *tidyfst* have provided the following functions to facilitate the workflow:
 
- parse_fst: Get information of the data.frame without reading it
- slice_fst: Select the target rows by number
- select_fst: Select the target columns for the task
- filter_fst: Conditional selection of rows
- import_fst: Read a fst file like `fst::read_fst` but always return a data.table
- export_fst: Write a fst file like `fst::write_fst` but always use largest compress factor (which yields smallest file)  

In such a workflow, you never need to read the whole data.frame into your RAM, you just select the target data, process them instantly and get the results all at once. You do not have to read the data to know the structure of data.frame, because we have `parse_fst`(a wrapper for `fst` in *fst* package). Now let's give it a try.
```{r setup}
library(tidyfst)

# Generate some random data frame with 10 million rows and various column types
nr_of_rows <- 1e7

df <- data.frame(
    Logical = sample(c(TRUE, FALSE, NA), prob = c(0.85, 0.1, 0.05), nr_of_rows, replace = TRUE),
    Integer = sample(1L:100L, nr_of_rows, replace = TRUE),
    Real = sample(sample(1:10000, 20) / 100, nr_of_rows, replace = TRUE),
    Factor = as.factor(sample(labels(UScitiesD), nr_of_rows, replace = TRUE))
  )

# write the fst file, make sure you do not have the file with same name in the directory
export_fst(df,"fst_test.fst")

# remove all variables in the environment
rm(list = ls())
```
 Now, we want to know the information in the data frame.
```{r}
parse_fst("fst_test.fst") -> ft
ft
```
 If we want to get the information in the `Factor` column, use:
```{r}
ft %>% 
  select_fst(Factor) %>% 
  count_dt(Factor) -> factor_info

factor_info
```
 If we want to calculate the mean of `Integer` by the group of `Factor`, use:
```{r}
ft %>% 
  select_fst(Integer,Factor) %>% 
  summarise_dt(avg = mean(Integer),by = Factor) -> avg_info

avg_info

```
 In this workflow, we only select/filter/slice the data we need, and get the results directly from the pipeline. Therefore, we read the minimum needed data into RAM and release it and save only the results we want. This workflow could save memory for many exploratory big data analysis.
 Last, let's delete the output file:
```{r}
# delete the output file
unlink("fst_test.fst")
```
 After (>=) version 0.9.3, tidyfst has also added a function `as_fst()`, which can turn any data.frame into a fst table and saved the data in the temporary file. This means that we might never have to save the object in the RAM ever (as long as it is a data.frame)! A small example:
```{r}
iris %>% as_fst() -> iris_fst
mtcars %>% as_fst() -> mtcars_fst

iris_fst
mtcars_fst
```
 So when you have generated a pretty large data.frame and do not want it to consume the cache in your computer, just save it and read it when needed using `as_fst`.
 
 
 

---
title: "tidyfst包实例分析"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{chinese_tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```


我的R语言小伙伴最近分享了自己使用R来做工业级数据清洗的[经验](https://mp.weixin.qq.com/s/NVlCPss32j6Ohdrc9Edx-A)，最近我自己在不断测试我的新包tidyfst，因此就拿这个data.table的案例来尝试一下。

## 测试数据构造
  本次测试，将不会加载data.table包，但是其实tidyfst里面无处不是data.table的元素，而且也导出了很多内置的data.table函数，比如as.data.table和data.table。所以这些代码在tidyfst中就可以自如地使用。
```{r}
library(tidyfst)
diamonds <- ggplot2::diamonds
n = 1e5  #如果想做工业级测试，可以继续增加数量
set.seed(2020)
dtranges <- seq.Date(from = as.Date("2011-01-01"),
                     to = as.Date("2020-01-01"),
                     by = 1)
n1 <- sample(nrow(diamonds), n, replace = TRUE)
dat1 <- as.data.table(diamonds[n1, ])
dat1[, "dt"] <- sample(dtranges, n, replace = TRUE)  # 增加dt列
n2 <- sample(nrow(dat1), nrow(dat1)/1000)
dat1[n2, "price"] <- NA # price列构造千分之一缺失值
dat2 <- data.table(dt = sample(dtranges, min(n/1000, length(dtranges))),
                   price1 = sample(1000, min(n/1000, length(dtranges)), replace = TRUE))

dat3 <- data.table(dt = sample(dtranges, min(n/1000, length(dtranges))),
                   price2 = sample(1000, min(n/1000, length(dtranges)), replace = TRUE))

print(dat1)
```

## 基础
### 小技巧
后面的分析，经常要根据日期进行计算。所以，先对日期进行排序，就能够提高运行速度。在tidyfst中，可以使用`arrange_dt`函数来对数据进行原位的各种操作，其中就包括排序。
```{r}
dat1 = arrange_dt(dat1,dt)
dat1
```
那么，现在dat1的数据就按照日期排好序了。

### 聚合
#### 1.求每种切割类型、每种颜色钻石的平均价格、中位数价格与最高价格
在tidyfst中，我设置了一个`sys_time_print`函数，可以方便地输出`system.time()`函数返回的结果。
```{r}
sys_time_print({
  r1_1 <- dat1 %>% 
    summarise_dt(
      by = .(cut,color),
      mean_price = mean(price, na.rm = TRUE),
      median_price = median(price, na.rm = TRUE),
      max_price = max(price, na.rm = TRUE)
    )
})
r1_1
```

tidyfst是永远不可能比data.table快的，但是如果你觉得上面的代码更容易掌握、更容易读懂，而在日常工作中多花零点几秒的运行时间没有太大问题（实际上节省了大家的交流时间，甚至就是节省将来自己再次读懂自己代码的时间），tidyfst就值得拥有。

#### 2.求每天最高出售价格对应的那笔订单

```{r}
sys_time_print({
  r1_2 <- dat1 %>% 
    arrange_dt(dt,-price) %>% 
    drop_na_dt(price) %>% 
    group_dt(
      by = dt,
      head(1)
    )
})
r1_2
```
### join
#### 1.dat1与dat2以dt列左连接
实质上，merge函数已经优化得很好。tidyfst设计`*_join`系列函数的时候，只是为了一种不一样的语法结构来帮助实现不同的连接，因为它确实更加直观一些。但是实质上它还是merge.data.table函数的包装版本。
```{r}
sys_time_print({
  r2_1 <- dat1 %>% 
    left_join_dt(dat2,by = "dt")
})
r2_1
```
#### 2.多重join
  
```{r}
sys_time_print({
  mymerge <- function(x, y) left_join_dt(x, y, by = "dt")
  r2_2 <- Reduce(mymerge, list(dat1, dat2, dat3))
})
r2_2
```
### 长宽表转换
#### 1.长表转宽表

```{r}
sys_time_print({
  mean1 <- function(x) mean(x, na.rm = TRUE)
  max1 <- function(x) max(x, na.rm = TRUE)
  r3_1 <-dat1 %>% 
    wider_dt(cut,
             value = c("depth", "price"),
             name = "color",
             fun = list(mean1,max1))
})
r3_1
```

#### 2.宽表转长表

```{r}
sys_time_print({
  r3_2 <-dat1 %>% 
    select_dt(cut,color,x,y,z) %>% 
    longer_dt(cut,color,
              name = "xyz",
              value = "xyzvalue")
})

r3_2 
```

## 高阶
### 向上/下填充空值
对于填充空值来说，可以这样操作：
```{r}
sys_time_print({
  dat1 %>% fill_na_dt(price) -> dat1
})
dat1

```


### 添加子维度聚合结果为新列

#### 1.以dat1为例，添加两列，一列为以cut、color聚合求price的均值，另一列是求标准差

```{r}

sys_time_print({
  mutate_dt(dat1,
           mean_price = mean(price, na.rm = TRUE),
           sd_price = sd(price, na.rm = TRUE),
           by = .(cut, color))
})

dat1
```

#### 2.以dat1为例，以dt分组添加一列序号id
```{r}

sys_time_print({
  dat1 %>% 
  group_dt(
    by = dt,
    mutate_dt(id = seq(.N))
  ) -> dat1
})
dat1
```

### 移动函数

```{r}

sys_time_print({
  dat1 %>% 
    group_dt(
      by = color,
      mutate_dt(
        MA10_price = frollmean(price, 10),
        MSD10_price = frollapply(price, 10, FUN = sd)
      )
    ) -> dat1
})

dat1
```


## 系统参数
```{r}
sessionInfo()
```

---
title: "Example 3: Reshape"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{example3_reshape}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```
  When I design `longer_dt` and `wider_dt`, I could find the `pivot_longer` and `pivot_wider` in `tidyr` and `melt` and `dcast` in `data.table`. Still, designing this API is not easy, as my goal is to let users use it with least pain. Here we would try to reproduce the results in the vignette of `tidyr`(<https://tidyr.tidyverse.org/articles/pivot.html>). First load the packages:
```{r setup}
library(tidyfst)
library(tidyr)
```

## Longer
  First inspect the data:
```{r}
relig_income
```
  In `tidyr`, to get the longer format you need:
```{r,eval=FALSE}
relig_income %>% 
  pivot_longer(-religion, names_to = "income", values_to = "count")
```
  In `tidyfst`, we have:
```{r,warning=FALSE}
relig_income %>% 
  longer_dt("religion",name = "income",value = "count")
```
  Another example from `tidyr`:
```{r,eval=FALSE}
billboard

# tidyr way:
 billboard %>%
   pivot_longer(
     cols = starts_with("wk"),
     names_to = "week",
     values_to = "rank",
     values_drop_na = TRUE
   )

# tidyfst way:
billboard %>% 
  longer_dt(-"wk",
            name = "week",
            value = "rank",
            na.rm = TRUE
            )
# regex could select the groups to keep, and minus could select the reverse
```
A warning would could come out because the merging column has different data types and do the coercion automatically.
    
## Wider

```{r}
## data
fish_encounters

## tidyr way:
fish_encounters %>% 
  pivot_wider(names_from = station, values_from = seen)

## tidyfst way:
fish_encounters %>% 
  wider_dt(name = "station",value = "seen")

# if no keeped groups are selected, use all except for name and value columns
```
    
  If you want to fill with 0s, use:
```{r}
fish_encounters %>% 
  wider_dt(name = "station",value = "seen",fill = 0)
```
  Note that the parameter of `name` and `value` should always be provided and should be explicit called (with the parameter names attached).

## More complicated example
This example comes from data.table (<https://rdatatable.gitlab.io/data.table/articles/datatable-reshape.html>), and has been used in tidyr too. We'll try to do it in tidyfst in this example.
If we have a data.frame as below:
```{r}
family <- fread("family_id age_mother dob_child1 dob_child2 dob_child3 gender_child1 gender_child2 gender_child3
1         30 1998-11-26 2000-01-29         NA             1             2            NA
2         27 1996-06-22         NA         NA             2            NA            NA
3         26 2002-07-11 2004-04-05 2007-09-02             2             2             1
4         32 2004-10-10 2009-08-27 2012-07-21             1             1             1
5         29 2000-12-05 2005-02-28         NA             2             1            NA")

family
```
And want to reshape the data.table to be like this:
```{r}
#     family_id age_mother  child        dob gender
#         <int>      <int> <char>     <char> <char>
#  1:         1         30 child1 1998-11-26      1
#  2:         1         30 child2 2000-01-29      2
#  3:         1         30 child3       <NA>   <NA>
#  4:         2         27 child1 1996-06-22      2
#  5:         2         27 child2       <NA>   <NA>
#  6:         2         27 child3       <NA>   <NA>
#  7:         3         26 child1 2002-07-11      2
#  8:         3         26 child2 2004-04-05      2
#  9:         3         26 child3 2007-09-02      1
# 10:         4         32 child1 2004-10-10      1
# 11:         4         32 child2 2009-08-27      1
# 12:         4         32 child3 2012-07-21      1
# 13:         5         29 child1 2000-12-05      2
# 14:         5         29 child2 2005-02-28      1
# 15:         5         29 child3       <NA>   <NA>
```

The `data.table::dcast` and `tidyr::pivot_longer` could transfer it in one step, however, not so easy to understand. Here we'll do it step by step to see what actually happens in this transfer.

```{r,warning=FALSE}
family %>% 
  longer_dt(1:2) %>% 
  separate_dt("name",into = c("class","child")) %>% 
  wider_dt(-"class|value",
           name = "class",
           value = "value")
```

In such a process, we could find that we actually get a longer table, then separate it, and wider it later. *tidyfst* is not going to support the complicated transfer in one step, because it might be easier to implement, but much harder to understand 3 procedures in 1 step. If you still prefer that way, use `data.table::dcast` and `tidyr::pivot_longer` instead.




---
title: "Example 4: Nest"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{example4_nest}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```
  This vignette displays how to use nesting in `tidyfst`. It has referred to `tidyr`s vignette in <https://tidyr.tidyverse.org/articles/nest.html>. Now fist, we nest the "mtcars" data.frame by "cyl" column.
```{r setup}
library(tidyfst)

# nest by "cyl" column
mtcars_nested <- mtcars %>% 
  nest_dt(cyl) # you can use "cyl" too, very flexible

# inspect the output data.table
mtcars_nested

```

  Now, we want to do a regression within the nested group "cyl". We'll use the famous `lapply` to complete this:

```{r}
mtcars_nested2 <- mtcars_nested %>% 
  mutate_dt(model = lapply(ndt,function(df) lm(mpg ~ wt, data = df)))

mtcars_nested2
```
  We could see that the model is stored in the column "model".
  Now, we try to get the fitted value in the model.
```{r}
mtcars_nested3 <- mtcars_nested2 %>% 
  mutate_dt(model_predict = lapply(model, predict))
mtcars_nested3$model_predict
```
  We could find that the "model_predict" is a list of numeric vectors. Let's try to unnest the target column "model_predict". 
```{r}
mtcars_nested3 %>% unnest_dt(model_predict)
```
  This process would remove all the other list column automatically. For instance, in our case, the column "ndt" is removed.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/complete.R
\name{complete_dt}
\alias{complete_dt}
\title{Complete a data frame with missing combinations of data}
\usage{
complete_dt(.data, ..., fill = NA)
}
\arguments{
\item{.data}{data.frame}

\item{...}{Specification of columns to expand.The selection of columns is
supported by the flexible \code{\link[tidyfst]{select_dt}}.
To find all unique combinations of provided columns, including those not found in the data,
supply each variable as a separate argument. But the two modes (select the
needed columns and fill outside values) could not be mixed,
find more details in examples.}

\item{fill}{Atomic value to fill into the missing cell, default uses \code{NA}.}
}
\value{
data.table
}
\description{
Turns implicit missing values into explicit missing values.
 All the combinations of column values (should be unique) will be constructed.
 Other columns will be filled with NAs or constant value.
}
\details{
When the provided columns with addtion data are of different length,
all the unique combinations would be returned. This operation should be used
only on unique entries, and it will always returned the unique entries.

If you supply fill parameter, these values will also replace existing explicit missing values in the data set.
}
\examples{
df <- data.table(
  group = c(1:2, 1),
  item_id = c(1:2, 2),
  item_name = c("a", "b", "b"),
  value1 = 1:3,
  value2 = 4:6
)

df \%>\% complete_dt(item_id,item_name)
df \%>\% complete_dt(item_id,item_name,fill = 0)
df \%>\% complete_dt("item")
df \%>\% complete_dt(item_id=1:3)
df \%>\% complete_dt(item_id=1:3,group=1:2)
df \%>\% complete_dt(item_id=1:3,group=1:3,item_name=c("a","b","c"))

}
\seealso{
\code{\link[tidyr]{complete}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/topn_dt.R
\name{top_dt}
\alias{top_dt}
\alias{top_n_dt}
\alias{top_frac_dt}
\title{Select top (or bottom) n rows (by value)}
\usage{
top_dt(.data, wt = NULL, n = NULL, prop = NULL)

top_n_dt(.data, n, wt = NULL)

top_frac_dt(.data, prop, wt = NULL)
}
\arguments{
\item{.data}{data.frame}

\item{wt}{(Optional). The variable to use for ordering.
If not specified, defaults to the last variable in the data.frame.}

\item{n}{Number of rows to return.
Will include more rows if there are ties.
If \code{n} is positive, selects the top rows.
If negative, select the bottom rows.}

\item{prop}{Fraction of rows to return.
Will include more rows if there are ties.
If \code{prop} is positive, selects the top rows.
If negative, select the bottom rows.}
}
\value{
data.table
}
\description{
Get the top entries (rows) according to the values of specified columns.
 One can get the top or bottom ones according to number or proportion.

In \code{top_dt}, you can use an API for both
functionalities in `top_n_dt()` and `top_frac_dt()`.
}
\examples{
iris \%>\% top_n_dt(10,Sepal.Length)
iris \%>\% top_n_dt(-10,Sepal.Length)
iris \%>\% top_frac_dt(.1,Sepal.Length)
iris \%>\% top_frac_dt(-.1,Sepal.Length)

# For `top_dt`, you can use both modes above
iris \%>\% top_dt(Sepal.Length,n = 10)
iris \%>\% top_dt(Sepal.Length,prop = .1)
}
\seealso{
\code{\link[dplyr]{top_n}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pairwise.R
\name{pairwise_count_dt}
\alias{pairwise_count_dt}
\title{Count pairs of items within a group}
\usage{
pairwise_count_dt(
  .data,
  .group,
  .value,
  upper = FALSE,
  diag = FALSE,
  sort = TRUE
)
}
\arguments{
\item{.data}{A data.frame.}

\item{.group}{Column name of counting group.}

\item{.value}{Item to count pairs, will end up in \code{V1} and \code{V2} columns.}

\item{upper}{When \code{FALSE}(Default), duplicated combinations would be removed.}

\item{diag}{Whether to include diagonal (V1==V2) in output. Default uses \code{FALSE}.}

\item{sort}{Whether to sort rows by counts. Default uses \code{TRUE}.}
}
\value{
A data.table with 3 columns (named as "V1","V2" and "n"), containing combinations
 in "V1" and "V2", and counts in "n".
}
\description{
Count the number of times each pair
of items appear together within a group.
For example, this could count the number of times two words appear within documents.
This function has referred to \code{pairwise_count} in \strong{widyr} package,
but with very different defaults on several parameters.
}
\examples{

dat <- data.table(group = rep(1:5, each = 2),
              letter = c("a", "b",
                         "a", "c",
                         "a", "c",
                         "b", "e",
                         "b", "f"))
pairwise_count_dt(dat,group,letter)
pairwise_count_dt(dat,group,letter,sort = FALSE)
pairwise_count_dt(dat,group,letter,diag = TRUE)
pairwise_count_dt(dat,group,letter,diag = TRUE,upper = TRUE)

# The column name could be specified using character.
pairwise_count_dt(dat,"group","letter")
}
\seealso{
\code{\link[widyr]{pairwise_count}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lag_lead.R
\name{lead_dt}
\alias{lead_dt}
\alias{lag_dt}
\title{Fast lead/lag for vectors}
\usage{
lead_dt(x, n = 1L, fill = NA)

lag_dt(x, n = 1L, fill = NA)
}
\arguments{
\item{x}{A vector}

\item{n}{a positive integer of length 1,
giving the number of positions to lead or lag by. Default uses 1}

\item{fill}{Value to use for padding when the window goes beyond the input length.
Default uses \code{NA}}
}
\value{
A vector
}
\description{
Find the "next" or "previous" values in a vector.
It has wrapped \pkg{data.table}'s \code{shift} function.
}
\examples{
lead_dt(1:5)
lag_dt(1:5)
lead_dt(1:5,2)
lead_dt(1:5,n = 2,fill = 0)
}
\seealso{
\code{\link[dplyr]{lead}},\code{\link[data.table]{shift}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_dt.R
\name{filter_dt}
\alias{filter_dt}
\title{Filter entries in data.frame}
\usage{
filter_dt(.data, ...)
}
\arguments{
\item{.data}{data.frame}

\item{...}{List of variables or name-value pairs of summary/modifications
functions.}
}
\value{
data.table
}
\description{
Choose rows where conditions are true.
}
\examples{
iris \%>\% filter_dt(Sepal.Length > 7)
iris \%>\% filter_dt(Sepal.Length == max(Sepal.Length))

# comma is not supported in tidyfst after v0.9.8
# which means you can't use:
\dontrun{
 iris \%>\% filter_dt(Sepal.Length > 7, Sepal.Width > 3)
}
# use following code instead
iris \%>\% filter_dt(Sepal.Length > 7 & Sepal.Width > 3)

}
\seealso{
\code{\link[dplyr]{filter}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_by.R
\name{group_by_dt}
\alias{group_by_dt}
\alias{group_exe_dt}
\title{Group by variable(s) and implement operations}
\usage{
group_by_dt(.data, ..., cols = NULL)

group_exe_dt(.data, ...)
}
\arguments{
\item{.data}{A data frame}

\item{...}{Variables to group by for \code{group_by_dt},
namely the columns to sort by. Do not quote the column names.
Any data manipulation arguments that could be
implemented on a data.frame for \code{group_exe_dt}.
It can receive what \code{select_dt} receives.}

\item{cols}{A character vector of column names to group by.}
}
\value{
A data.table with keys
}
\description{
Carry out data manipulation within specified groups. Different from \code{group_dt},
the implementation is split into two operations, namely grouping and implementation.

Using \code{setkey} and \code{setkeyv} in \pkg{data.table}
to carry out \code{group_by}-like functionalities in \pkg{dplyr}. This is
not only convenient but also efficient in computation.
}
\details{
\code{group_by_dt} and \code{group_exe_dt} are a pair of functions
to be used in combination. It utilizes the feature of key setting in data.table,
which provides high performance for group operations, especially when you have
to operate by specific groups frequently.
}
\examples{

# aggregation after grouping using group_exe_dt
as.data.table(iris) -> a
a \%>\%
  group_by_dt(Species) \%>\%
  group_exe_dt(head(1))

a \%>\%
  group_by_dt(Species) \%>\%
  group_exe_dt(
    head(3) \%>\%
      summarise_dt(sum = sum(Sepal.Length))
  )

mtcars \%>\%
  group_by_dt("cyl|am") \%>\%
  group_exe_dt(
    summarise_dt(mpg_sum = sum(mpg))
  )
# equals to
mtcars \%>\%
  group_by_dt(cols = c("cyl","am")) \%>\%
  group_exe_dt(
    summarise_dt(mpg_sum = sum(mpg))
  )
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distinct_dt.R
\name{distinct_dt}
\alias{distinct_dt}
\title{Select distinct/unique rows in data.frame}
\usage{
distinct_dt(.data, ..., .keep_all = FALSE)
}
\arguments{
\item{.data}{data.frame}

\item{...}{Optional variables to use when determining uniqueness.
If there are multiple rows for a given combination of inputs,
only the first row will be preserved.
If omitted, will use all variables.}

\item{.keep_all}{If \code{TRUE}, keep all variables in data.frame. If a combination of ... is not distinct,
this keeps the first row of values.}
}
\value{
data.table
}
\description{
Select only unique/distinct rows from a data frame.
}
\examples{
iris \%>\% distinct_dt()
iris \%>\% distinct_dt(Species)
iris \%>\% distinct_dt(Species,.keep_all = TRUE)
mtcars \%>\% distinct_dt(cyl,vs)
mtcars \%>\% distinct_dt(cyl,vs,.keep_all = TRUE)

}
\seealso{
\code{\link[dplyr]{distinct}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_prop.R
\name{percent}
\alias{percent}
\alias{add_prop}
\title{Add percentage to counts in data.frame}
\usage{
percent(x, digits = 1)

add_prop(.data, count_name = last(names(.data)), digits = 1)
}
\arguments{
\item{x}{A number (numeric).}

\item{digits}{How many digits to keep in the percentage. Default uses 1.}

\item{.data}{A data frame.}

\item{count_name}{Column name of counts (Character).
Default uses the last column of data.frame.}
}
\description{
Add percentage for counts in the data.frame, both numeric and
 character with `%` would be provided.
}
\examples{

 percent(0.9057)
 percent(0.9057,3)

 iris \%>\%
   count_dt(Species) \%>\%
   add_prop()

 iris \%>\%
   count_dt(Species) \%>\%
   add_prop(count_name = "n",digits = 2)

}
\references{
https://stackoverflow.com/questions/7145826/how-to-format-a-number-as-percentage-in-r
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rownames.R
\name{rn_col}
\alias{rn_col}
\alias{col_rn}
\title{Tools for working with row names}
\usage{
rn_col(.data, var = "rowname")

col_rn(.data, var = "rowname")
}
\arguments{
\item{.data}{A data.frame.}

\item{var}{Name of column to use for rownames.}
}
\value{
\code{rn_col} returns a data.table,
\code{col_rn} returns a data frame.
}
\description{
The enhanced data.frame, including tibble and data.table, do not
support row names. To link to some base r facilities, there should be functions
to save information in row names. These functions are analogous to
\code{rownames_to_column} and \code{column_to_rownames} in \pkg{tibble}.
}
\examples{

 mtcars \%>\% rn_col()
 mtcars \%>\% rn_col("rn")

 mtcars \%>\% rn_col() -> new_mtcars

 new_mtcars \%>\% col_rn() -> old_mtcars
 old_mtcars
 setequal(mtcars,old_mtcars)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sql_join.R
\name{sql_join}
\alias{sql_join}
\alias{sql_join_dt}
\title{Case insensitive table joining like SQL}
\usage{
sql_join_dt(x, y, by = NULL, type = "inner", suffix = c(".x", ".y"))
}
\arguments{
\item{x}{A data.table}

\item{y}{A data.table}

\item{by}{(Optional) A character vector of variables to join by.

  If `NULL`, the default, `*_join_dt()` will perform a natural join, using all
  variables in common across `x` and `y`. A message lists the variables so that you
  can check they're correct; suppress the message by supplying `by` explicitly.

  To join by different variables on `x` and `y`, use a named vector.
  For example, `by = c("a" = "b")` will match `x$a` to `y$b`.

  To join by multiple variables, use a vector with length > 1.
  For example, `by = c("a", "b")` will match `x$a` to `y$a` and `x$b` to
  `y$b`. Use a named vector to match different variables in `x` and `y`.
  For example, `by = c("a" = "b", "c" = "d")` will match `x$a` to `y$b` and
  `x$c` to `y$d`.

  Notice that in `sql_join`, the joining variables would turn to upper case
  in the output table.}

\item{type}{Which type of join would you like to use?
Default uses "inner", other options include
"left", "right", "full", "anti", "semi".}

\item{suffix}{If there are non-joined duplicate variables in x and y, these
suffixes will be added to the output to disambiguate them. Should be a
character vector of length 2.}
}
\value{
A data.table
}
\description{
Work like the `*_join_dt` series functions, joining
tables with common or customized keys  in various ways. The only
difference is the joining is case insensitive like SQL.
}
\examples{
dt1 = data.table(x = c("A","b"),y = 1:2)
dt2 = data.table(x = c("a","B"),z = 4:5)
sql_join_dt(dt1,dt2)
}
\seealso{
\code{\link[tidyfst]{join}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impute_dt.R
\name{impute_dt}
\alias{impute_dt}
\title{Impute missing values with mean, median or mode}
\usage{
impute_dt(.data, ..., .func = "mode")
}
\arguments{
\item{.data}{A data.frame}

\item{...}{Columns to select}

\item{.func}{Character, "mode" (default), "mean" or "median".
Could also define it by oneself.}
}
\value{
A data.table
}
\description{
Impute the columns of data.frame with its mean, median or mode.
}
\examples{

Pclass <- c(3, 1, 3, 1, 3, 2, 2, 3, NA, NA)
Sex <- c('male', 'male', 'female', 'female', 'female',
         'female', NA, 'male', 'female', NA)
Age <- c(22, 38, 26, 35, NA,
         45, 25, 39, 28, 40)
SibSp <- c(0, 1, 3, 1, 2, 3, 2, 2, NA, 0)
Fare <- c(7.25, 71.3, 7.92, NA, 8.05, 8.46, 51.9, 60, 32, 15)
Embarked <- c('S', NA, 'S', 'Q', 'Q', 'S', 'C', 'S', 'C', 'S')
data <- data.frame('Pclass' = Pclass,
 'Sex' = Sex, 'Age' = Age, 'SibSp' = SibSp,
 'Fare' = Fare, 'Embarked' = Embarked)

data
data \%>\% impute_dt() # defalut uses "mode" as `.func`
data \%>\% impute_dt(is.numeric,.func = "mean")
data \%>\% impute_dt(is.numeric,.func = "median")

my_fun = function(x){
  x[is.na(x)] = (max(x,na.rm = TRUE) - min(x,na.rm = TRUE))/2
  x
}
data \%>\% impute_dt(is.numeric,.func = my_fun)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_dt.R
\name{count_dt}
\alias{count_dt}
\alias{add_count_dt}
\title{Count observations by group}
\usage{
count_dt(.data, ..., sort = TRUE, .name = "n")

add_count_dt(.data, ..., .name = "n")
}
\arguments{
\item{.data}{data.table/data.frame data.frame will be automatically converted
to data.table.}

\item{...}{Variables to group by, could receive what `select_dt` receives.}

\item{sort}{logical. If TRUE result will be sorted in desending order by resulting variable.}

\item{.name}{character. Name of resulting variable. Default uses "n".}
}
\value{
data.table
}
\description{
Count the unique values of one or more variables.
}
\examples{
iris \%>\% count_dt(Species)
iris \%>\% count_dt(Species,.name = "count")
iris \%>\% add_count_dt(Species)
iris \%>\% add_count_dt(Species,.name = "N")

mtcars \%>\% count_dt(cyl,vs)
mtcars \%>\% count_dt("cyl|vs")
mtcars \%>\% count_dt(cyl,vs,.name = "N",sort = FALSE)
mtcars \%>\% add_count_dt(cyl,vs)
mtcars \%>\% add_count_dt("cyl|vs")

}
\seealso{
\code{\link[dplyr]{count}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slice.R
\name{slice_dt}
\alias{slice_dt}
\alias{slice_head_dt}
\alias{slice_tail_dt}
\alias{slice_max_dt}
\alias{slice_min_dt}
\alias{slice_sample_dt}
\title{Subset rows using their positions}
\usage{
slice_dt(.data, ..., by = NULL)

slice_head_dt(.data, n, by = NULL)

slice_tail_dt(.data, n, by = NULL)

slice_max_dt(.data, order_by, n, by = NULL, with_ties = TRUE)

slice_min_dt(.data, order_by, n, by = NULL, with_ties = TRUE)

slice_sample_dt(.data, n, replace = FALSE, by = NULL)
}
\arguments{
\item{.data}{A data.table}

\item{...}{Provide either positive values to keep, or negative values to drop.
The values provided must be either all positive or all negative.}

\item{by}{Slice by which group(s)?}

\item{n}{When larger than or equal to 1, the number of rows.
When between 0 and 1, the proportion of rows to select.}

\item{order_by}{Variable or function of variables to order by.}

\item{with_ties}{Should ties be kept together? The default, `TRUE`,
may return more rows than you request. Use `FALSE` to ignore ties,
and return the first `n` rows.}

\item{replace}{Should sampling be performed with (`TRUE`) or without
(`FALSE`, the default) replacement.}
}
\value{
A data.table
}
\description{
`slice_dt()` lets you index rows by their (integer) locations. It allows you
to select, remove, and duplicate rows. It is accompanied by a number of
helpers for common use cases:

* `slice_head_dt()` and `slice_tail_dt()` select the first or last rows.
* `slice_sample_dt()` randomly selects rows.
* `slice_min_dt()` and `slice_max_dt()` select rows with highest or lowest values
  of a variable.
}
\examples{

a = iris
slice_dt(a,1,2)
slice_dt(a,2:3)
slice_dt(a,141:.N)
slice_dt(a,1,.N)
slice_head_dt(a,5)
slice_head_dt(a,0.1)
slice_tail_dt(a,5)
slice_tail_dt(a,0.1)
slice_max_dt(a,Sepal.Length,10)
slice_max_dt(a,Sepal.Length,10,with_ties = FALSE)
slice_min_dt(a,Sepal.Length,10)
slice_min_dt(a,Sepal.Length,10,with_ties = FALSE)
slice_sample_dt(a,10)
slice_sample_dt(a,0.1)


# use by to slice by group

## following codes get the same results
slice_dt(a,1:3,by = "Species")
slice_dt(a,1:3,by = Species)
slice_dt(a,1:3,by = .(Species))

slice_head_dt(a,2,by = Species)
slice_tail_dt(a,2,by = Species)

slice_max_dt(a,Sepal.Length,3,by = Species)
slice_max_dt(a,Sepal.Length,3,by = Species,with_ties = FALSE)
slice_min_dt(a,Sepal.Length,3,by = Species)
slice_min_dt(a,Sepal.Length,3,by = Species,with_ties = FALSE)

# in `slice_sample_dt`, "by" could only take character class
slice_sample_dt(a,.1,by = "Species")
slice_sample_dt(a,3,by = "Species")
slice_sample_dt(a,51,replace = TRUE,by = "Species")

}
\seealso{
\code{\link[dplyr]{slice}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cummean.R
\name{cummean}
\alias{cummean}
\title{Cumulative mean}
\usage{
cummean(x)
}
\arguments{
\item{x}{a numeric or complex object,
or an object that can be coerced to one of these.}
}
\description{
Returns a vector whose elements are the cumulative mean of the elements of the argument.
}
\examples{
cummean(1:10)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mutate_dt.R
\name{mutate_dt}
\alias{mutate_dt}
\alias{transmute_dt}
\title{Mutate columns in data.frame}
\usage{
mutate_dt(.data, ..., by)

transmute_dt(.data, ..., by)
}
\arguments{
\item{.data}{data.frame}

\item{...}{List of variables or name-value pairs of summary/modifications
functions.}

\item{by}{(Optional) Mutate by what group?}
}
\value{
data.table
}
\description{
Adds or updates columns in data.frame.
}
\examples{

iris \%>\% mutate_dt(one = 1,Sepal.Length = Sepal.Length + 1)
iris \%>\% transmute_dt(one = 1,Sepal.Length = Sepal.Length + 1)
# add group number with symbol `.GRP`
iris \%>\% mutate_dt(id = 1:.N,grp = .GRP,by = Species)

}
\seealso{
\code{\link[dplyr]{mutate}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fst.R
\name{fst}
\alias{fst}
\alias{parse_fst}
\alias{slice_fst}
\alias{select_fst}
\alias{filter_fst}
\alias{summary_fst}
\title{Parse,inspect and extract data.table from fst file}
\usage{
parse_fst(path)

slice_fst(ft, row_no)

select_fst(ft, ...)

filter_fst(ft, ...)

summary_fst(ft)
}
\arguments{
\item{path}{path to fst file}

\item{ft}{An object of class fst_table, returned by \code{parse_fst}}

\item{row_no}{An integer vector (Positive)}

\item{...}{The filter conditions}
}
\value{
\code{parse_fst} returns a fst_table class.

\code{select_fst} and \code{filter_fst} returns a data.table.
}
\description{
A tookit of APIs for reading fst file as data.table, could select by column, row and conditional filtering.
}
\details{
\code{summary_fst} could provide some basic information about
the fst table.
}
\examples{

\dontrun{
  fst::write_fst(iris,"iris_test.fst")
  # parse the file but not reading it
  parse_fst("iris_test.fst") -> ft
  ft

  class(ft)
  lapply(ft,class)
  names(ft)
  dim(ft)
  summary_fst(ft)

  # get the data by query
  ft \%>\% slice_fst(1:3)
  ft \%>\% slice_fst(c(1,3))

  ft \%>\% select_fst(Sepal.Length)
  ft \%>\% select_fst(Sepal.Length,Sepal.Width)
  ft \%>\% select_fst("Sepal.Length")
  ft \%>\% select_fst(1:3)
  ft \%>\% select_fst(1,3)
  ft \%>\% select_fst("Se")
  ft \%>\% select_fst("nothing")
  ft \%>\% select_fst("Se|Sp")
  ft \%>\% select_fst(cols = names(iris)[2:3])

  ft \%>\% filter_fst(Sepal.Width > 3)
  ft \%>\% filter_fst(Sepal.Length > 6 , Species == "virginica")
  ft \%>\% filter_fst(Sepal.Length > 6 & Species == "virginica" & Sepal.Width < 3)

  unlink("iris_test.fst")
}
}
\seealso{
\code{\link[fst]{fst}}, \code{\link[fst]{metadata_fst}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dt.R
\name{in_dt}
\alias{in_dt}
\alias{as_dt}
\title{Short cut to data.table}
\usage{
in_dt(.data, ...)

as_dt(.data)
}
\arguments{
\item{.data}{A data.frame}

\item{...}{Recieve \code{B} in data.table's \code{A[B]} syntax.}
}
\description{
To use facilities provided by \pkg{data.table}, but do not have to
load \pkg{data.table} package.
}
\details{
The \code{as_dt} could turn any data frame to data.table class. If the data is
not a data frame, return error.

The \code{in_dt} function creates a virtual environment in data.table, it could be
piped well because it still follows the principals of \pkg{tidyfst}, which are: (1) Never
use in place replacement and (2) Always recieves a data frame (data.frame/tibble/data.table)
and returns a data.table. Therefore, the in place functions like \code{:=} will still
return the results.
}
\examples{
iris \%>\% as_dt()
iris \%>\% in_dt(order(-Sepal.Length),.SD[.N],by=Species)
}
\seealso{
\code{\link[data.table]{data.table}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transpose.R
\name{t_dt}
\alias{t_dt}
\title{Efficient transpose of data.frame}
\usage{
t_dt(.data)
}
\arguments{
\item{.data}{A data.frame/data.table/tibble}
}
\value{
A transposed data.frame
}
\description{
An efficient way to transpose data frames(data.frame/data.table/tibble).
}
\details{
This function would return the original data.frame structure,
keeping all the row names and column names. If the row names are not
available or, "V1,V2..." will be provided.
}
\examples{

t_dt(iris)
t_dt(mtcars)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nest_dt.R
\name{nest_dt}
\alias{nest_dt}
\alias{unnest_dt}
\alias{squeeze_dt}
\alias{chop_dt}
\alias{unchop_dt}
\title{Nest and unnest}
\usage{
nest_dt(.data, ..., mcols = NULL, .name = "ndt")

unnest_dt(.data, ...)

squeeze_dt(.data, ..., .name = "ndt")

chop_dt(.data, ...)

unchop_dt(.data, ...)
}
\arguments{
\item{.data}{data.table, nested or unnested}

\item{...}{The variables for nest group(for \code{nest_dt}),
columns to be nested(for \code{squeeze_dt} and \code{chop_dt}),
or column(s) to be unnested(for \code{unnest_dt}).
Could recieve anything that \code{\link[tidyfst]{select_dt}} could receive.}

\item{mcols}{Name-variable pairs in the list, form like}

\item{.name}{Character. The nested column name. Defaults to "ndt".
\code{list(petal="^Pe",sepal="^Se")}, see example.}
}
\value{
data.table, nested or unnested
}
\description{
Create or melt list columns in data.frame.

Analogous function for \code{nest} and \code{unnest} in \pkg{tidyr}.
\code{unnest_dt} will automatically remove other list-columns except for the
target list-columns (which would be unnested later). Also, \code{squeeze_dt} is
designed to merge multiple columns into list column.
}
\details{
In the \code{nest_dt}, the data would be nested to a column named `ndt`,
 which is short for nested data.table.

The \code{squeeze_dt} would not remove the originial columns.

The \code{unchop_dt} is the reverse operation of \code{chop_dt}.

These functions are experiencing the experimental stage, especially
the \code{unnest_dt}. If they don't work on some circumtances, try \pkg{tidyr}
package.
}
\examples{

# examples for nest_dt
# nest by which columns?
 mtcars \%>\% nest_dt(cyl)
 mtcars \%>\% nest_dt("cyl")
 mtcars \%>\% nest_dt(cyl,vs)
 mtcars \%>\% nest_dt(vs:am)
 mtcars \%>\% nest_dt("cyl|vs")
 mtcars \%>\% nest_dt(c("cyl","vs"))

 # change the nested column name
 mtcars \%>\% nest_dt(cyl,.name = "data")

# nest two columns directly
iris \%>\% nest_dt(mcols = list(petal="^Pe",sepal="^Se"))

# nest more flexibly
iris \%>\% nest_dt(mcols = list(ndt1 = 1:3,
  ndt2 = "Pe",
  ndt3 = Sepal.Length:Sepal.Width))

# examples for unnest_dt
# unnest which column?
 mtcars \%>\% nest_dt("cyl|vs") \%>\%
   unnest_dt(ndt)
 mtcars \%>\% nest_dt("cyl|vs") \%>\%
   unnest_dt("ndt")

df <- data.table(
  a = list(c("a", "b"), "c"),
  b = list(c(TRUE,TRUE),FALSE),
  c = list(3,c(1,2)),
  d = c(11, 22)
)

df
df \%>\% unnest_dt(a)
df \%>\% unnest_dt(2)
df \%>\% unnest_dt("c")
df \%>\% unnest_dt(cols = names(df)[3])

# You can unnest multiple columns simultaneously
df \%>\% unnest_dt(1:3)
df \%>\% unnest_dt(a,b,c)
df \%>\% unnest_dt("a|b|c")

# examples for squeeze_dt
# nest which columns?
iris \%>\% squeeze_dt(1:2)
iris \%>\% squeeze_dt("Se")
iris \%>\% squeeze_dt(Sepal.Length:Petal.Width)
iris \%>\% squeeze_dt(1:2,.name = "data")

# examples for chop_dt
df <- data.table(x = c(1, 1, 1, 2, 2, 3), y = 1:6, z = 6:1)
df \%>\% chop_dt(y,z)
df \%>\% chop_dt(y,z) \%>\% unchop_dt(y,z)
}
\references{
https://www.r-bloggers.com/much-faster-unnesting-with-data-table/

https://stackoverflow.com/questions/25430986/create-nested-data-tables-by-collapsing-rows-into-new-data-tables
}
\seealso{
\code{\link[tidyr]{nest}}, \code{\link[tidyr]{chop}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sys_time_print.R
\name{sys_time_print}
\alias{sys_time_print}
\title{Convenient print of time taken}
\usage{
sys_time_print(expr)
}
\arguments{
\item{expr}{Valid R expression to be timed.}
}
\value{
A character vector of the form HH:MM:SS,
or SS.MMMsec if under 60 seconds. See examples.
}
\description{
Convenient printing of time elapsed. A wrapper of
\code{data.table::timetaken}, but showing the results more directly.
}
\examples{

sys_time_print(Sys.sleep(1))

a = iris
sys_time_print({
  res = iris \%>\%
    mutate_dt(one = 1)
})
res
}
\seealso{
\code{\link[data.table]{timetaken}}, \code{\link[base]{system.time}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_dt.R
\name{group_dt}
\alias{group_dt}
\alias{rowwise_dt}
\title{Data manipulation within groups}
\usage{
group_dt(.data, by = NULL, ...)

rowwise_dt(.data, ...)
}
\arguments{
\item{.data}{A data.frame}

\item{by}{Variables to group by,unquoted name of grouping variable of list of unquoted names of grouping variables.}

\item{...}{Any data manipulation arguments that could be implemented on a data.frame.}
}
\value{
data.table
}
\description{
Carry out data manipulation within specified groups.
}
\details{
If you want to use \code{summarise_dt} and \code{mutate_dt} in
\code{group_dt}, it is better to use the "by" parameter in those functions,
that would be much faster because you don't have to use \code{.SD} (which takes
extra time to copy).
}
\examples{
iris \%>\% group_dt(by = Species,slice_dt(1:2))
iris \%>\% group_dt(Species,filter_dt(Sepal.Length == max(Sepal.Length)))
iris \%>\% group_dt(Species,summarise_dt(new = max(Sepal.Length)))

# you can pipe in the `group_dt`
iris \%>\% group_dt(Species,
                  mutate_dt(max= max(Sepal.Length)) \%>\%
                    summarise_dt(sum=sum(Sepal.Length)))

# for users familiar with data.table, you can work on .SD directly
# following codes get the first and last row from each group
iris \%>\%
  group_dt(
    by = Species,
    rbind(.SD[1],.SD[.N])
  )

#' # for summarise_dt, you can use "by" to calculate within the group
mtcars \%>\%
  summarise_dt(
   disp = mean(disp),
   hp = mean(hp),
   by = cyl
)

  # but you could also, of course, use group_dt
 mtcars \%>\%
   group_dt(by =.(vs,am),
     summarise_dt(avg = mean(mpg)))

  # and list of variables could also be used
 mtcars \%>\%
   group_dt(by =list(vs,am),
            summarise_dt(avg = mean(mpg)))

# examples for `rowwise_dt`
df <- data.table(x = 1:2, y = 3:4, z = 4:5)

df \%>\% mutate_dt(m = mean(c(x, y, z)))

df \%>\% rowwise_dt(
  mutate_dt(m = mean(c(x, y, z)))
)
}
\references{
https://stackoverflow.com/questions/36802385/use-by-each-row-for-data-table
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/separate.R
\name{separate_dt}
\alias{separate_dt}
\title{Separate a character column into two columns using
a regular expression separator}
\usage{
separate_dt(
  .data,
  separated_colname,
  into,
  sep = "[^[:alnum:]]+",
  remove = TRUE
)
}
\arguments{
\item{.data}{A data frame.}

\item{separated_colname}{Column to be separated, can be a character or alias.}

\item{into}{Character vector of length 2.}

\item{sep}{Separator between columns.}

\item{remove}{If \code{TRUE}, remove input column from output data frame.}
}
\description{
Given either regular expression,
\code{separate_dt()} turns a single character column into two columns.
}
\examples{
df <- data.frame(x = c(NA, "a.b", "a.d", "b.c"))
df \%>\% separate_dt(x, c("A", "B"))
# equals to
df \%>\% separate_dt("x", c("A", "B"))

# If you just want the second variable:
df \%>\% separate_dt(x,into = c(NA,"B"))
}
\seealso{
\code{\link[tidyr]{separate}}, \code{\link[tidyfst]{unite_dt}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/replace_dt.R
\name{replace_dt}
\alias{replace_dt}
\title{Fast value replacement in data frame}
\usage{
replace_dt(.data, ..., from = is.nan, to = NA)
}
\arguments{
\item{.data}{A data.frame}

\item{...}{Colunms to be replaced. If not specified, use all columns.}

\item{from}{A value, a vector of values or a function returns a logical value.
Defaults to \code{is.nan}.}

\item{to}{A value. Defaults to \code{NA}.}
}
\value{
A data.table.
}
\description{
While \code{replace_na_dt} could replace all NAs to another
value, \code{replace_dt} could replace any value(s) to another specific
value.
}
\examples{
iris \%>\% mutate_vars(is.factor,as.character) -> new_iris

new_iris \%>\%
  replace_dt(Species, from = "setosa",to = "SS")
new_iris \%>\%
  replace_dt(Species,from = c("setosa","virginica"),to = "sv")
new_iris \%>\%
  replace_dt(Petal.Width, from = .2,to = 2)
new_iris \%>\%
  replace_dt(from = .2,to = NA)
new_iris \%>\%
  replace_dt(is.numeric, from = function(x) x > 3, to = 9999 )
}
\seealso{
\code{\link[tidyfst]{replace_na_dt}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wider_dt.R
\name{wider_dt}
\alias{wider_dt}
\title{Pivot data from long to wide}
\usage{
wider_dt(.data, ..., name, value = NULL, fun = NULL, fill = NA)
}
\arguments{
\item{.data}{A data.frame}

\item{...}{Optional. The unchanged group in the transformation.
Could use integer vector, could receive what \code{select_dt} receives.}

\item{name}{Chracter.One column name of class to spread}

\item{value}{Chracter.One column name of value to spread.
If \code{NULL}, use all other variables.}

\item{fun}{Should the data be aggregated before casting?
Defaults to \code{NULL}, which uses \code{length} for aggregation.
If a function is provided, with aggregated by this function.}

\item{fill}{Value with which to fill missing cells. Default uses \code{NA}.}
}
\value{
data.table
}
\description{
Transform a data frame from long format to wide by increasing the number of columns and decreasing the number of rows.
}
\details{
The parameter of `name` and `value` should always
be provided and should be explicit called (with the parameter names attached).
}
\examples{
 stocks = data.frame(
   time = as.Date('2009-01-01') + 0:9,
   X = rnorm(10, 0, 1),
   Y = rnorm(10, 0, 2),
   Z = rnorm(10, 0, 4)
 ) \%>\%
   longer_dt(time) -> longer_stocks

 longer_stocks

 longer_stocks \%>\%
   wider_dt("time",
            name = "name",
            value = "value")

 longer_stocks \%>\%
   mutate_dt(one = 1) \%>\%
   wider_dt("time",
            name = "name",
            value = "one")

## using "fun" parameter for aggregation
DT <- data.table(v1 = rep(1:2, each = 6),
                 v2 = rep(rep(1:3, 2), each = 2),
                 v3 = rep(1:2, 6),
                 v4 = rnorm(6))
## for each combination of (v1, v2), add up all values of v4
DT \%>\%
  wider_dt(v1,v2,
           value = "v4",
           name = ".",
           fun = sum)
}
\seealso{
\code{\link[tidyfst]{longer_dt}},
 \code{\link[data.table]{dcast}},
 \code{\link[tidyr]{pivot_wider}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/missing.R
\name{drop_na_dt}
\alias{drop_na_dt}
\alias{replace_na_dt}
\alias{delete_na_cols}
\alias{delete_na_rows}
\alias{fill_na_dt}
\alias{shift_fill}
\title{Dump, replace and fill missing values in data.frame}
\usage{
drop_na_dt(.data, ...)

replace_na_dt(.data, ..., to)

delete_na_cols(.data, prop = NULL, n = NULL)

delete_na_rows(.data, prop = NULL, n = NULL)

fill_na_dt(.data, ..., direction = "down")

shift_fill(x, direction = "down")
}
\arguments{
\item{.data}{data.frame}

\item{...}{Colunms to be replaced or filled. If not specified, use all columns.}

\item{to}{What value should NA replace by?}

\item{prop}{If proportion of NAs is larger than or equal to "prop", would be deleted.}

\item{n}{If number of NAs is larger than or equal to "n", would be deleted.}

\item{direction}{Direction in which to fill missing values.
Currently either "down" (the default) or "up".}

\item{x}{A vector with missing values to be filled.}
}
\value{
data.table
}
\description{
A set of tools to deal with missing values in data.frames. 
It can dump, replace, fill (with next or previous observation) or delete entries according to their missing values.
}
\details{
\code{drop_na_dt} drops the entries with NAs in specific columns.
\code{fill_na_dt} fill NAs with observations ahead ("down") or below ("up"),
which is also known as last observation carried forward (LOCF) and
next observation carried backward(NOCB).

\code{delete_na_cols} could drop the columns with NA proportion larger
than or equal to "prop" or NA number larger than or equal to "n",
\code{delete_na_rows} works alike but deals with rows.

\code{shift_fill} could fill a vector with missing values.
}
\examples{
df <- data.table(x = c(1, 2, NA), y = c("a", NA, "b"))
 df \%>\% drop_na_dt()
 df \%>\% drop_na_dt(x)
 df \%>\% drop_na_dt(y)
 df \%>\% drop_na_dt(x,y)

 df \%>\% replace_na_dt(to = 0)
 df \%>\% replace_na_dt(x,to = 0)
 df \%>\% replace_na_dt(y,to = 0)
 df \%>\% replace_na_dt(x,y,to = 0)

 df \%>\% fill_na_dt(x)
 df \%>\% fill_na_dt() # not specified, fill all columns
 df \%>\% fill_na_dt(y,direction = "up")

x = data.frame(x = c(1, 2, NA, 3), y = c(NA, NA, 4, 5),z = rep(NA,4))
x
x \%>\% delete_na_cols()
x \%>\% delete_na_cols(prop = 0.75)
x \%>\% delete_na_cols(prop = 0.5)
x \%>\% delete_na_cols(prop = 0.24)
x \%>\% delete_na_cols(n = 2)

x \%>\% delete_na_rows(prop = 0.6)
x \%>\% delete_na_rows(n = 2)

# shift_fill
y = c("a",NA,"b",NA,"c")

shift_fill(y) # equals to
shift_fill(y,"down")

shift_fill(y,"up")
}
\references{
https://stackoverflow.com/questions/23597140/how-to-find-the-percentage-of-nas-in-a-data-frame

https://stackoverflow.com/questions/2643939/remove-columns-from-dataframe-where-all-values-are-na

https://stackoverflow.com/questions/7235657/fastest-way-to-replace-nas-in-a-large-data-table
}
\seealso{
\code{\link[tidyr]{drop_na}},\code{\link[tidyr]{replace_na}},
\code{\link[tidyr]{fill}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/longer_dt.R
\name{longer_dt}
\alias{longer_dt}
\title{Pivot data from wide to long}
\usage{
longer_dt(.data, ..., name = "name", value = "value", na.rm = FALSE)
}
\arguments{
\item{.data}{A data.frame}

\item{...}{Pattern for unchanged group or unquoted names. Pattern can accept
regular expression to match column names. It can recieve what \code{select_dt}
recieves.}

\item{name}{Name for the measured variable names column.
The default name is 'name'.}

\item{value}{Name for the molten data values column(s).
The default name is 'value'.}

\item{na.rm}{If \code{TRUE}, \code{NA} values will be removed from the molten data.}
}
\value{
A data.table
}
\description{
Turning a wide table to its longer form. It takes multiple columns and collapses into key-value pairs.
}
\examples{

## Example 1:
stocks = data.frame(
  time = as.Date('2009-01-01') + 0:9,
  X = rnorm(10, 0, 1),
  Y = rnorm(10, 0, 2),
  Z = rnorm(10, 0, 4)
)

stocks

stocks \%>\%
  longer_dt(time)

stocks \%>\%
  longer_dt("ti")

# Example 2:

\donttest{
  library(tidyr)

  billboard \%>\%
    longer_dt(
      -"wk",
      name = "week",
      value = "rank",
      na.rm = TRUE
    )

  # or use:
  billboard \%>\%
    longer_dt(
      artist,track,date.entered,
      name = "week",
      value = "rank",
      na.rm = TRUE
    )

  # or use:
  billboard \%>\%
    longer_dt(
      1:3,
      name = "week",
      value = "rank",
      na.rm = TRUE
    )
}
}
\seealso{
\code{\link[tidyfst]{wider_dt}},
  \code{\link[data.table]{melt}},
  \code{\link[tidyr]{pivot_longer}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setops.R
\name{intersect_dt}
\alias{intersect_dt}
\alias{union_dt}
\alias{setdiff_dt}
\alias{setequal_dt}
\title{Set operations for data frames}
\usage{
intersect_dt(x, y, all = FALSE)

union_dt(x, y, all = FALSE)

setdiff_dt(x, y, all = FALSE)

setequal_dt(x, y, all = TRUE)
}
\arguments{
\item{x}{A data.frame}

\item{y}{A data.frame}

\item{all}{Logical. When \code{FALSE} (default),
removes duplicate rows on the result.}
}
\value{
A data.table
}
\description{
Wrappers of set operations in \pkg{data.table}.
Only difference is it could be applied to non-data.table data frames by
recognizing and coercing them to data.table automatically.
}
\examples{

x = iris[c(2,3,3,4),]
x2 = iris[2:4,]
y = iris[c(3:5),]

intersect_dt(x, y)            # intersect
intersect_dt(x, y, all=TRUE)  # intersect all
setdiff_dt(x, y)              # except
setdiff_dt(x, y, all=TRUE)    # except all
union_dt(x, y)                # union
union_dt(x, y, all=TRUE)      # union all
setequal_dt(x, x2, all=FALSE) # setequal
setequal_dt(x, x2)            # setequal all

}
\seealso{
\code{\link[data.table]{setops}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dummy_dt.R
\name{dummy_dt}
\alias{dummy_dt}
\title{Fast creation of dummy variables}
\usage{
dummy_dt(.data, ..., longname = TRUE)
}
\arguments{
\item{.data}{data.frame}

\item{...}{Columns you want to create dummy variables from.
Very flexible, find in the examples.}

\item{longname}{logical. Should the output column labeled with the
original column name? Default uses \code{TRUE}.}
}
\value{
data.table
}
\description{
Quickly create dummy (binary) columns from character and factor type columns in the inputted data (and numeric columns if specified.)
This function is useful for statistical analysis when you want binary columns rather than character columns.
}
\details{
If no columns provided, will return the original data frame.

This function is inspired by \pkg{fastDummies} package, but provides
simple and precise usage, whereas \code{fastDummies::dummy_cols} provides more
features for statistical usage.
}
\examples{
iris \%>\% dummy_dt(Species)
iris \%>\% dummy_dt(Species,longname = FALSE)

mtcars \%>\% head() \%>\% dummy_dt(vs,am)
mtcars \%>\% head() \%>\% dummy_dt("cyl|gear")
}
\seealso{
\code{\link[fastDummies]{dummy_cols}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unite.R
\name{unite_dt}
\alias{unite_dt}
\title{Unite multiple columns into one by pasting strings together}
\usage{
unite_dt(
  .data,
  united_colname,
  ...,
  sep = "_",
  remove = FALSE,
  na2char = FALSE
)
}
\arguments{
\item{.data}{A data frame.}

\item{united_colname}{The name of the new column, string only.}

\item{...}{A selection of columns. If want to select all columns,
pass "" to the parameter. See example.}

\item{sep}{Separator to use between values.}

\item{remove}{If \code{TRUE}, remove input columns from output data frame.}

\item{na2char}{If \code{FALSE}, missing values would be merged into \code{NA},
otherwise \code{NA} is treated as character "NA". This is different from
\pkg{tidyr}.}
}
\description{
Convenience function to paste together multiple columns into one.
}
\examples{
df <- expand.grid(x = c("a", NA), y = c("b", NA))
df

# Treat missing value as NA, default
df \%>\% unite_dt("z", x:y, remove = FALSE)
# Treat missing value as character "NA"
df \%>\% unite_dt("z", x:y, na2char = TRUE, remove = FALSE)
df \%>\%
  unite_dt("xy", x:y)

# Select all columns
iris \%>\% unite_dt("merged_name","")
}
\seealso{
\code{\link[tidyr]{unite}},\code{\link[tidyfst]{separate_dt}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rename_dt.R
\name{rename_dt}
\alias{rename_dt}
\alias{rename_with_dt}
\title{Rename column in data.frame}
\usage{
rename_dt(.data, ...)

rename_with_dt(.data, .fn, ...)
}
\arguments{
\item{.data}{data.frame}

\item{...}{statements of rename, e.g. `sl = Sepal.Length` means the column named
as "Sepal.Length" would be renamed to "sl"}

\item{.fn}{A function used to transform the selected columns.
Should return a character vector the same length as the input.}
}
\value{
data.table
}
\description{
Rename one or more columns in the data.frame.
}
\examples{
iris \%>\%
  rename_dt(sl = Sepal.Length,sw = Sepal.Width) \%>\%
  head()
iris \%>\% rename_with_dt(toupper)
iris \%>\% rename_with_dt(toupper,"^Pe")

}
\seealso{
\code{\link[dplyr]{rename}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select_dt.R
\name{select_dt}
\alias{select_dt}
\alias{select_mix}
\title{Select column from data.frame}
\usage{
select_dt(.data, ..., cols = NULL, negate = FALSE)

select_mix(.data, ..., rm.dup = TRUE)
}
\arguments{
\item{.data}{data.frame}

\item{...}{List of variables or name-value pairs of summary/modifications
functions. It can also recieve conditional function to select columns.
When starts with `-`(minus symbol) or `!`, return the negative columns.}

\item{cols}{(Optional)A numeric or character vector.}

\item{negate}{Applicable when regular expression and "cols" is used.
If \code{TRUE}, return the non-matched pattern. Default uses \code{FALSE}.}

\item{rm.dup}{Should duplicated columns be removed? Defaults to \code{TRUE}.}
}
\value{
data.table
}
\description{
Select specific column(s) via various ways. One can select columns by their column names, indexes or regular expression recognizing the column name(s).
}
\examples{
iris \%>\% select_dt(Species)
iris \%>\% select_dt(Sepal.Length,Sepal.Width)
iris \%>\% select_dt(Sepal.Length:Petal.Length)
iris \%>\% select_dt(-Sepal.Length)
iris \%>\% select_dt(-Sepal.Length,-Petal.Length)
iris \%>\% select_dt(-(Sepal.Length:Petal.Length))
iris \%>\% select_dt(c("Sepal.Length","Sepal.Width"))
iris \%>\% select_dt(-c("Sepal.Length","Sepal.Width"))
iris \%>\% select_dt(1)
iris \%>\% select_dt(-1)
iris \%>\% select_dt(1:3)
iris \%>\% select_dt(-(1:3))
iris \%>\% select_dt(1,3)
iris \%>\% select_dt("Pe")
iris \%>\% select_dt(-"Se")
iris \%>\% select_dt(!"Se")
iris \%>\% select_dt("Pe",negate = TRUE)
iris \%>\% select_dt("Pe|Sp")
iris \%>\% select_dt(cols = 2:3)
iris \%>\% select_dt(cols = 2:3,negate = TRUE)
iris \%>\% select_dt(cols = c("Sepal.Length","Sepal.Width"))
iris \%>\% select_dt(cols = names(iris)[2:3])

iris \%>\% select_dt(is.factor)
iris \%>\% select_dt(-is.factor)
iris \%>\% select_dt(!is.factor)

# select_mix could provide flexible mix selection
select_mix(iris, Species,"Sepal.Length")
select_mix(iris,1:2,is.factor)

select_mix(iris,Sepal.Length,is.numeric)
# set rm.dup to FALSE could save the duplicated column names
select_mix(iris,Sepal.Length,is.numeric,rm.dup = FALSE)

}
\seealso{
\code{\link[dplyr]{select}}, \code{\link[dplyr]{select_if}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/_global_setting.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{\%>\%}
\alias{\%like\%}
\alias{data.table}
\alias{as.data.table}
\alias{fread}
\alias{fwrite}
\alias{fintersect}
\alias{fsetdiff}
\alias{funion}
\alias{fsetequal}
\alias{frollapply}
\alias{fcoalesce}
\alias{uniqueN}
\alias{rbindlist}
\alias{tables}
\alias{like}
\alias{copy}
\alias{key}
\alias{CJ}
\alias{rleid}
\alias{rleidv}
\alias{fcase}
\alias{between}
\alias{set}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{data.table}{\code{\link[data.table:like]{\%like\%}}, \code{\link[data.table:J]{CJ}}, \code{\link[data.table]{as.data.table}}, \code{\link[data.table]{between}}, \code{\link[data.table]{copy}}, \code{\link[data.table]{data.table}}, \code{\link[data.table]{fcase}}, \code{\link[data.table:coalesce]{fcoalesce}}, \code{\link[data.table:setops]{fintersect}}, \code{\link[data.table]{fread}}, \code{\link[data.table:froll]{frollapply}}, \code{\link[data.table:setops]{fsetdiff}}, \code{\link[data.table:setops]{fsetequal}}, \code{\link[data.table:setops]{funion}}, \code{\link[data.table]{fwrite}}, \code{\link[data.table:setkey]{key}}, \code{\link[data.table]{like}}, \code{\link[data.table]{rbindlist}}, \code{\link[data.table]{rleid}}, \code{\link[data.table:rleid]{rleidv}}, \code{\link[data.table:assign]{set}}, \code{\link[data.table]{tables}}, \code{\link[data.table:duplicated]{uniqueN}}}

  \item{stringr}{\code{\link[stringr:pipe]{\%>\%}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nth.R
\name{nth}
\alias{nth}
\title{Extract the nth value from a vector}
\usage{
nth(v, n = 1)
}
\arguments{
\item{v}{A vector}

\item{n}{A single integer specifying the position. Default uses \code{1}.
Negative integers index from the end
 (i.e. -1L will return the last value in the vector).
 If a double is supplied, it will be silently truncated.}
}
\value{
A single value.
}
\description{
Get the value from a vector with its position.
}
\examples{

x = 1:10
nth(x, 1)
nth(x, 5)
nth(x, -2)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pull_dt.R
\name{pull_dt}
\alias{pull_dt}
\title{Pull out a single variable}
\usage{
pull_dt(.data, col)
}
\arguments{
\item{.data}{data.frame}

\item{col}{A name of column or index (should be positive).}
}
\value{
vector
}
\description{
Extract vector from data.frame, works likt `[[`. Analogous function for \code{pull} in \pkg{dplyr}
}
\examples{
mtcars \%>\% pull_dt(2)
mtcars \%>\% pull_dt(cyl)
mtcars \%>\% pull_dt("cyl")
}
\seealso{
\code{\link[dplyr]{pull}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mutate_when.R
\name{mutate_when}
\alias{mutate_when}
\alias{mutate_vars}
\title{Conditional update of columns in data.table}
\usage{
mutate_when(.data, when, ..., by)

mutate_vars(.data, .cols = NULL, .func, ..., by)
}
\arguments{
\item{.data}{data.frame}

\item{when}{An object which can be coerced to logical mode}

\item{...}{Name-value pairs of expressions for \code{mutate_when}.
Additional parameters to be passed to parameter '.func' in \code{mutate_vars}.}

\item{by}{(Optional) Mutate by what group?}

\item{.cols}{Any types that can be accepted by \code{\link[tidyfst]{select_dt}}.}

\item{.func}{Function to be run within each column, should return a value or
vectors with same length.}
}
\value{
data.table
}
\description{
Update or add columns when the given condition is met.

\code{mutate_when} integrates \code{mutate} and \code{case_when}
in \pkg{dplyr} and make a new tidy verb for data.table. \code{mutate_vars} is
 a super function to do updates in specific columns according to conditions.
}
\examples{
iris[3:8,]
iris[3:8,] \%>\%
  mutate_when(Petal.Width == .2,
              one = 1,Sepal.Length=2)

iris \%>\% mutate_vars("Pe",scale)
iris \%>\% mutate_vars(is.numeric,scale)
iris \%>\% mutate_vars(-is.factor,scale)
iris \%>\% mutate_vars(1:2,scale)
iris \%>\% mutate_vars(.func = as.character)
}
\seealso{
\code{\link[tidyfst]{select_dt}}, \code{\link[dplyr]{case_when}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/col_max.R
\name{col_max}
\alias{col_max}
\alias{col_min}
\title{Get the column name of the max/min number each row}
\usage{
col_max(.data, .name = "max_col")

col_min(.data, .name = "min_col")
}
\arguments{
\item{.data}{A data.frame with numeric column(s)}

\item{.name}{The column name of the new added column}
}
\value{
A data.table
}
\description{
For a data.frame with numeric values, add a new column
specifying the column name of the first max/min value each row.
}
\examples{
set.seed(199057)
DT <- data.table(matrix(sample(10, 100, TRUE), ncol=10))
DT
col_max(DT)
col_max(DT,.name = "max_col_name")
col_min(DT)

col_max(iris)
}
\references{
https://stackoverflow.com/questions/17735859/for-each-row-return-the-column-name-of-the-largest-value
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{print_options}
\alias{print_options}
\title{Set global printing method for data.table}
\usage{
print_options(
  topn = 5,
  nrows = 100,
  class = TRUE,
  row.names = TRUE,
  col.names = "auto",
  print.keys = TRUE,
  trunc.cols = FALSE
)
}
\arguments{
\item{topn}{The number of rows to be printed from the beginning and
end of tables with more than \code{nrow} rows.}

\item{nrows}{The number of rows which will be printed before truncation is enforced.}

\item{class}{If \code{TRUE}, the resulting output will include above each column its storage class (or a self-evident abbreviation thereof).}

\item{row.names}{If \code{TRUE}, row indices will be printed.}

\item{col.names}{One of three flavours for controlling the display of column names in output. \code{"auto"} includes column names above the data, as well as below the table if \code{nrow(x) > 20}. \code{"top"} excludes this lower register when applicable, and \code{"none"} suppresses column names altogether (as well as column classes if \code{class = TRUE}.}

\item{print.keys}{If \code{TRUE}, any \code{\link{key}} and/or \code{\link[=indices]{index}} currently assigned to \code{x} will be printed prior to the preview of the data.}

\item{trunc.cols}{If \code{TRUE}, only the columns that can be printed in the console without wrapping the columns to new lines will be printed (similar to \code{tibbles}).}
}
\value{
None. This function is used for its side effect of changing options.
}
\description{
This function allow user to define how data.table is printed.
}
\details{
Notice that \pkg{tidyfst} has a slightly different printing default for data.table,
 which is it always prints the keys and variable class (not like \pkg{data.table}).
}
\examples{

iris \%>\% as.data.table()
print_options(topn = 3,trunc.cols = TRUE)
iris \%>\% as.data.table()

# set all settings to default in tidyfst
print_options()
iris \%>\% as.data.table()

}
\seealso{
\code{\link[data.table]{print.data.table}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rec.R
\name{rec}
\alias{rec}
\alias{rec_num}
\alias{rec_char}
\title{Recode number or strings}
\usage{
rec_num(x, rec, keep = TRUE)

rec_char(x, rec, keep = TRUE)
}
\arguments{
\item{x}{A numeric or character vector.}

\item{rec}{String with recode pairs of old and new values.
Find the usage in examples.}

\item{keep}{Logical. Decide whether to keep the original values if not recoded.
Defaults to \code{TRUE}.}
}
\value{
A vector.
}
\description{
Recode discrete variables, including numerice and character
variable.
}
\examples{

x = 1:10
x
rec_num(x, rec = "1=10; 4=2")
rec_num(x, rec = "1:3=1; 4:6=2")
rec_num(x, rec = "1:3=1; 4:6=2",keep = FALSE)

y = letters[1:5]
y
rec_char(y,rec = "a=A;b=B")
rec_char(y,rec = "a,b=A;c,d=B")
rec_char(y,rec = "a,b=A;c,d=B",keep = FALSE)

}
\seealso{
\code{\link[sjmisc]{rec}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utf8_encoding.R
\name{utf8_encoding}
\alias{utf8_encoding}
\title{Use UTF-8 for character encoding in a data frame}
\usage{
utf8_encoding(.data)
}
\arguments{
\item{.data}{A data.frame.}
}
\value{
A data.table with characters in UTF-8 encoding
}
\description{
\code{fread} from \pkg{data.table} could not recognize the encoding
and return the correct form, this could be unconvenient for text mining tasks. The
\code{utf8-encoding} could use "UTF-8" as the encoding to override the current
encoding of characters in a data frame.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uncount.R
\name{uncount_dt}
\alias{uncount_dt}
\title{"Uncount" a data frame}
\usage{
uncount_dt(.data, wt, .remove = TRUE)
}
\arguments{
\item{.data}{A data.frame}

\item{wt}{A vector of weights.}

\item{.remove}{Should the column for \code{weights} be removed?
Default uses \code{TRUE}.}
}
\description{
Duplicating rows according to a weighting variable.
 This is the opposite operation of `count_dt`.
 Analogous to `tidyr::uncount`.
}
\examples{

df <- data.table(x = c("a", "b"), n = c(1, 2))
uncount_dt(df, n)
uncount_dt(df,n,FALSE)
}
\seealso{
\code{\link[dplyr]{count}}, \code{\link[tidyr]{uncount}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample.R
\name{sample_dt}
\alias{sample_dt}
\alias{sample_n_dt}
\alias{sample_frac_dt}
\title{Sample rows randomly from a table}
\usage{
sample_dt(.data, n = NULL, prop = NULL, replace = FALSE, by = NULL)

sample_n_dt(.data, size, replace = FALSE, by = NULL)

sample_frac_dt(.data, size, replace = FALSE, by = NULL)
}
\arguments{
\item{.data}{A data.frame}

\item{n}{Number of rows to select}

\item{prop}{Fraction of rows to select}

\item{replace}{Sample with or without replacement? Default uses \code{FALSE}.}

\item{by}{(Optional) Character. Specify if you want to sample by group.}

\item{size}{For \code{sample_n_dt}, the number of rows to select.
For \code{sample_frac_dt}, the fraction of rows to select.}
}
\value{
data.table
}
\description{
Select a number or proportion of rows randomly from the data frame

\code{sample_dt} is a merged version of \code{sample_n_dt} and
\code{sample_frac_dt}, this could be convenient.
}
\examples{
sample_n_dt(mtcars, 10)
sample_n_dt(mtcars, 50, replace = TRUE)
sample_frac_dt(mtcars, 0.1)
sample_frac_dt(mtcars, 1.5, replace = TRUE)


sample_dt(mtcars,n=10)
sample_dt(mtcars,prop = 0.1)


# sample by group(s)
iris \%>\% sample_n_dt(2,by = "Species")
iris \%>\% sample_frac_dt(.1,by = "Species")

mtcars \%>\% sample_n_dt(1,by = "cyl,vs")
# equals to
mtcars \%>\% sample_n_dt(1,by = c("cyl","vs"))
}
\seealso{
\code{\link[dplyr]{sample_n}},\code{\link[dplyr]{sample_frac}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/join.R
\name{join}
\alias{join}
\alias{inner_join_dt}
\alias{left_join_dt}
\alias{right_join_dt}
\alias{full_join_dt}
\alias{anti_join_dt}
\alias{semi_join_dt}
\title{Join tables}
\usage{
inner_join_dt(x, y, by = NULL, on = NULL, suffix = c(".x", ".y"))

left_join_dt(x, y, by = NULL, on = NULL, suffix = c(".x", ".y"))

right_join_dt(x, y, by = NULL, on = NULL, suffix = c(".x", ".y"))

full_join_dt(x, y, by = NULL, on = NULL, suffix = c(".x", ".y"))

anti_join_dt(x, y, by = NULL, on = NULL)

semi_join_dt(x, y, by = NULL, on = NULL)
}
\arguments{
\item{x}{A data.table}

\item{y}{A data.table}

\item{by}{(Optional) A character vector of variables to join by.

  If `NULL`, the default, `*_join_dt()` will perform a natural join, using all
  variables in common across `x` and `y`. A message lists the variables so that you
  can check they're correct; suppress the message by supplying `by` explicitly.

  To join by different variables on `x` and `y`, use a named vector.
  For example, `by = c("a" = "b")` will match `x$a` to `y$b`.

  To join by multiple variables, use a vector with length > 1.
  For example, `by = c("a", "b")` will match `x$a` to `y$a` and `x$b` to
  `y$b`. Use a named vector to match different variables in `x` and `y`.
  For example, `by = c("a" = "b", "c" = "d")` will match `x$a` to `y$b` and
  `x$c` to `y$d`.}

\item{on}{(Optional)
Indicate which columns in x should be joined with which columns in y.
Examples included:
  1.\code{.by = c("a","b")} (this is a must for \code{set_full_join_dt});
  2.\code{.by = c(x1="y1", x2="y2")};
  3.\code{.by = c("x1==y1", "x2==y2")};
  4.\code{.by = c("a", V2="b")};
  5.\code{.by = .(a, b)};
  6.\code{.by = c("x>=a", "y<=b")} or \code{.by = .(x>=a, y<=b)}.}

\item{suffix}{If there are non-joined duplicate variables in x and y, these
suffixes will be added to the output to disambiguate them. Should be a
character vector of length 2.}
}
\value{
A data.table
}
\description{
The mutating joins add columns from `y` to `x`,
matching rows based on the keys:

* `inner_join_dt()`: includes all rows in `x` and `y`.
* `left_join_dt()`: includes all rows in `x`.
* `right_join_dt()`: includes all rows in `y`.
* `full_join_dt()`: includes all rows in `x` or `y`.

Filtering joins filter rows from `x` based on the presence or absence
of matches in `y`:

* `semi_join_dt()` return all rows from `x` with a match in `y`.
* `anti_join_dt()` return all rows from `x` without a match in `y`.
}
\examples{

workers = fread("
    name company
    Nick Acme
    John Ajax
    Daniela Ajax
")

positions = fread("
    name position
    John designer
    Daniela engineer
    Cathie manager
")

workers \%>\% inner_join_dt(positions)
workers \%>\% left_join_dt(positions)
workers \%>\% right_join_dt(positions)
workers \%>\% full_join_dt(positions)

# filtering joins
workers \%>\% anti_join_dt(positions)
workers \%>\% semi_join_dt(positions)

# To suppress the message, supply 'by' argument
workers \%>\% left_join_dt(positions, by = "name")

# Use a named 'by' if the join variables have different names
positions2 = setNames(positions, c("worker", "position")) # rename first column in 'positions'
workers \%>\% inner_join_dt(positions2, by = c("name" = "worker"))

# the syntax of 'on' could be a bit different
workers \%>\% inner_join_dt(positions2,on = "name==worker")


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_fst.R
\name{as_fst}
\alias{as_fst}
\title{Save a data.frame as a fst table}
\usage{
as_fst(.data)
}
\arguments{
\item{.data}{A data.frame}
}
\value{
An object of class \code{fst_table}
}
\description{
This function first export the data.frame to a temporal file,
and then parse it back as a fst table (class name is "fst_table").
}
\examples{

\dontrun{
  iris \%>\%
    as_fst() -> iris_fst
  iris_fst
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fst_io.R
\name{export_fst}
\alias{export_fst}
\alias{import_fst}
\title{Read and write fst files}
\usage{
export_fst(x, path, compress = 100, uniform_encoding = TRUE)

import_fst(
  path,
  columns = NULL,
  from = 1,
  to = NULL,
  as.data.table = TRUE,
  old_format = FALSE
)
}
\arguments{
\item{x}{a data frame to write to disk}

\item{path}{path to fst file}

\item{compress}{value in the range 0 to 100, indicating the amount of compression to use.
Lower values mean larger file sizes. The default compression is set to 50.}

\item{uniform_encoding}{If `TRUE`, all character vectors will be assumed to have elements with equal encoding.
The encoding (latin1, UTF8 or native) of the first non-NA element will used as encoding for the whole column.
This will be a correct assumption for most use cases.
If `uniform.encoding` is set to `FALSE`, no such assumption will be made and all elements will be converted
to the same encoding. The latter is a relatively expensive operation and will reduce write performance for
character columns.}

\item{columns}{Column names to read. The default is to read all columns.}

\item{from}{Read data starting from this row number.}

\item{to}{Read data up until this row number. The default is to read to the last row of the stored dataset.}

\item{as.data.table}{If TRUE, the result will be returned as a \code{data.table} object. Any keys set on
dataset \code{x} before writing will be retained. This allows for storage of sorted datasets. This option
requires \code{data.table} package to be installed.}

\item{old_format}{must be FALSE, the old fst file format is deprecated and can only be read and
converted with fst package versions 0.8.0 to 0.8.10.}
}
\value{
`import_fst` returns a data.table with the selected columns and rows. `export_fst`
writes `x` to a `fst` file and invisibly returns `x` (so you can use this function in a pipeline).
}
\description{
Wrapper for \code{\link[fst]{read_fst}} and \code{\link[fst]{write_fst}}
from \pkg{fst}, but use a different default. For data import, always return a data.table.
For data export, always compress the data to the smallest size.
}
\examples{
\dontrun{
export_fst(iris,"iris_fst_test.fst")
iris_dt = import_fst("iris_fst_test.fst")
iris_dt
unlink("iris_fst_test.fst")
}
}
\seealso{
\code{\link[fst]{read_fst}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arrange_dt.R
\name{arrange_dt}
\alias{arrange_dt}
\title{Arrange entries in data.frame}
\usage{
arrange_dt(.data, ...)
}
\arguments{
\item{.data}{data.frame}

\item{...}{Arrange by what group? Minus symbol means arrange by
descending order.}
}
\value{
data.table
}
\description{
Order the rows of a data frame rows by the values of selected columns.
}
\examples{

iris \%>\% arrange_dt(Sepal.Length)

# minus for decreasing order
iris \%>\% arrange_dt(-Sepal.Length)

# arrange by multiple variables
iris \%>\% arrange_dt(Sepal.Length,Petal.Length)

}
\seealso{
\code{\link[dplyr]{arrange}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidymat.R
\name{mat_df}
\alias{mat_df}
\alias{df_mat}
\title{Conversion between tidy table and named matrix}
\usage{
mat_df(m)

df_mat(df, row, col, value)
}
\arguments{
\item{m}{A matrix}

\item{df}{A data.frame with at least 3 columns, one for row name,
one for column name, and one for values. The names for column and
row should be unique.}

\item{row}{Unquoted expression of column name for row}

\item{col}{Unquoted expression of column name for column}

\item{value}{Unquoted expression of column name for values}
}
\value{
For \code{mat_df}, a data.frame.
For \code{df_mat}, a named matrix.
}
\description{
Convenient fucntions to implement conversion between
 tidy table and named matrix.
}
\examples{

mm = matrix(c(1:8,NA),ncol = 3,dimnames = list(letters[1:3],LETTERS[1:3]))
mm
tdf = mat_df(mm)
tdf
mat = df_mat(tdf,row,col,value)
setequal(mm,mat)

tdf \%>\%
  setNames(c("A","B","C")) \%>\%
  df_mat(A,B,C)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarise_dt.R
\name{summarise_dt}
\alias{summarise_dt}
\alias{summarize_dt}
\alias{summarise_when}
\alias{summarize_when}
\alias{summarise_vars}
\alias{summarize_vars}
\title{Summarise columns to single values}
\usage{
summarise_dt(.data, ..., by = NULL)

summarize_dt(.data, ..., by = NULL)

summarise_when(.data, when, ..., by = NULL)

summarize_when(.data, when, ..., by = NULL)

summarise_vars(.data, .cols = NULL, .func, ..., by)

summarize_vars(.data, .cols = NULL, .func, ..., by)
}
\arguments{
\item{.data}{data.frame}

\item{...}{List of variables or name-value pairs of summary/modifications
functions for \code{summarise_dt}.Additional parameters to be passed to
 parameter '.func' in \code{summarise_vars}.}

\item{by}{unquoted name of grouping variable of list of unquoted names of
grouping variables. For details see \link[data.table]{data.table}}

\item{when}{An object which can be coerced to logical mode}

\item{.cols}{Columns to be summarised.}

\item{.func}{Function to be run within each column, should return a value or vectors with same length.}
}
\value{
data.table
}
\description{
Summarise group of values into one value for each group. If there is only one group, then only one value would be returned.
 The summarise function should always return a single value.
}
\details{
\code{summarise_vars} could complete summarise on specific columns.
}
\examples{
iris \%>\% summarise_dt(avg = mean(Sepal.Length))
iris \%>\% summarise_dt(avg = mean(Sepal.Length),by = Species)
mtcars \%>\% summarise_dt(avg = mean(hp),by = .(cyl,vs))

# the data.table way
mtcars \%>\% summarise_dt(cyl_n = .N, by = .(cyl, vs)) # `.` is short for list

iris \%>\% summarise_vars(is.numeric,min)
iris \%>\% summarise_vars(-is.factor,min)
iris \%>\% summarise_vars(1:4,min)

iris \%>\% summarise_vars(is.numeric,min,by ="Species")
mtcars \%>\% summarise_vars(is.numeric,mean,by = "vs,am")

# use multiple functions on multiple columns
iris \%>\%
  summarise_vars(is.numeric,.func = list(mean,sd,median))
iris \%>\%
  summarise_vars(is.numeric,.func = list(mean,sd,median),by = Species)

}
\seealso{
\code{\link[dplyr]{summarise}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/object_size.R
\name{object_size}
\alias{object_size}
\title{Nice printing of report the Space Allocated for an Object}
\usage{
object_size(object)
}
\arguments{
\item{object}{an R object.}
}
\value{
An object of class "object_size"
}
\description{
Provides an estimate of the memory that is being used to store an R object.
A wrapper of `object.size`, but use a nicer printing unit.
}
\examples{

iris \%>\% object_size()

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/relocate_dt.R
\name{relocate_dt}
\alias{relocate_dt}
\title{Change column order}
\usage{
relocate_dt(.data, ..., how = "first", where = NULL)
}
\arguments{
\item{.data}{A data.frame}

\item{...}{Columns to move}

\item{how}{The mode of movement, including "first","last","after","before".
Default uses "first".}

\item{where}{Destination of columns selected by \code{...}.
Applicable for "after" and "before" mode.}
}
\value{
A data.table with rearranged columns.
}
\description{
Change the position of columns,
 using the same syntax as `select_dt()`. Check similar function
 as `relocate` in \pkg{dplyr}.
}
\examples{
df <- data.table(a = 1, b = 1, c = 1, d = "a", e = "a", f = "a")
df
df \%>\% relocate_dt(f)
df \%>\% relocate_dt(a,how = "last")

df \%>\% relocate_dt(is.character)
df \%>\% relocate_dt(is.numeric, how = "last")
df \%>\% relocate_dt("[aeiou]")

df \%>\% relocate_dt(a, how = "after",where = f)
df \%>\% relocate_dt(f, how = "before",where = a)
df \%>\% relocate_dt(f, how = "before",where = c)
df \%>\% relocate_dt(f, how = "after",where = c)

df2 <- data.table(a = 1, b = "a", c = 1, d = "a")
df2 \%>\% relocate_dt(is.numeric,
                    how = "after",
                    where = is.character)
df2 \%>\% relocate_dt(is.numeric,
                    how="before",
                    where = is.character)
}
\seealso{
\code{\link[dplyr]{relocate}}
}
