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
