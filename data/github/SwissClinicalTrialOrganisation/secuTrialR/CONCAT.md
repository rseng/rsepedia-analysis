<!-- README.md is generated from README.Rmd. Please edit that file -->

secuTrialR
==========

[![](https://img.shields.io/badge/dev%20version-1.0.9-blue.svg)](https://github.com/SwissClinicalTrialOrganisation/secuTrialR)
[![](https://www.r-pkg.org/badges/version/secuTrialR?color=green)](https://cran.r-project.org/package=secuTrialR)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/SwissClinicalTrialOrganisation/secuTrialR?branch=master&svg=true)](https://ci.appveyor.com/project/SwissClinicalTrialOrganisation/secuTrialR)
[![travis](https://api.travis-ci.com/SwissClinicalTrialOrganisation/secuTrialR.svg?branch=master)](https://travis-ci.com/github/SwissClinicalTrialOrganisation/secuTrialR)
[![Actions
Status](https://github.com/SwissClinicalTrialOrganisation/secuTrialR/workflows/R-CMD-check/badge.svg)](https://github.com/SwissClinicalTrialOrganisation/secuTrialR/actions)
[![codecov](https://codecov.io/github/SwissClinicalTrialOrganisation/secuTrialR/branch/master/graphs/badge.svg)](https://codecov.io/github/SwissClinicalTrialOrganisation/secuTrialR)

An R package to handle data from the clinical data management system
(CDMS) [secuTrial](https://www.secutrial.com/en/).

Installing from GitHub with devtools
------------------------------------

Please note that `R versions >= 3.5` should be used to run `secuTrialR`.

    devtools::install_github("SwissClinicalTrialOrganisation/secuTrialR")

Recommended export options
--------------------------

While the package strives to allow loading of as many types of secuTrial
data exports as possible, there are certain export options which are
less likely to cause issues. If possible it is suggested to export data
which adheres to a suggested option set. Thus, we suggest to work with
exports which: \* are **zipped** \* are **English** \* have **reference
values** stored **in a separate table** \* contain **Add-IDs**, **centre
information**, **structure information**, **form status**, **project
setup** \* do **NOT** have the **meta data duplicated** into all tables
\* are **UTF-8** encoded \* are **“CSV format”** or **“CSV format for MS
Excel”** \* do **NOT** contain form **data of hidden fields**

If you use `read_secuTrial()` to read your export then it will inform
you regarding deviations.

Basic usage
-----------

An extensive applied manual/vignette is available
[here](https://github.com/SwissClinicalTrialOrganisation/secuTrialR/blob/master/vignettes/secuTrialR-package-vignette.pdf)
and probably the best place to get started.

Load the package

    library(secuTrialR)

Load a dataset

    export_location <- system.file("extdata", "sT_exports", "lnames",
                                   "s_export_CSV-xls_CTU05_long_ref_miss_en_utf8.zip",
                                   package = "secuTrialR")
    ctu05 <- read_secuTrial(export_location)

    ## Read export successfully.

    ## The following export options deviate from the suggested specifications:
    ## Data from hidden fields is part of the export.

This will load all sheets from the export into an object of class
`secuTrialdata`, which is basically a list. It will always contain
`export_details` (which are parsed from the HTML ExportOptions file that
secuTrial generates). By default, it will also contain all other files
in the dataset. secuTrialR automatically strips file names of dates. The
new file names can be seen via `ctu05$export_options$data_names`. The
function also adds [labels to variables](#variable-labels) and
data.frames, converts [categorical variables to
`factor`s](#prepare-factors) and ensures that [dates are `Date`s and
date-times are `POSIXct`](#prepare-dates). `read_secuTrial` is a wrapper
for the functions described below, so it is possible to achieve more
flexibility by using the individual functions (if necessary). Individual
tables can be extracted from the `ctu05` object via `tab <- ctu05$tab`,
where `tab` is the table of interest.

<details>

<summary>Wrapped functions</summary>

#### Load the dataset

    # prepare path to example export
    export_location <- system.file("extdata", "sT_exports", "BMD",
                                   "s_export_CSV-xls_BMD_short_en_utf8.zip",
                                   package = "secuTrialR")
    # load all export data
    bmd_export <- read_secuTrial_raw(data_dir = export_location)

    # load a second dataset
    export_location <- system.file("extdata", "sT_exports", "lnames",
                                   "s_export_CSV-xls_CTU05_long_ref_miss_en_utf8.zip",
                                   package = "secuTrialR")
    ctu05_raw <- read_secuTrial_raw(export_location)

    # View names of the bmd_export object
    names(bmd_export)

    ##  [1] "export_options" "fs"             "cn"             "ctr"           
    ##  [5] "is"             "qs"             "qac"            "vp"            
    ##  [9] "vpfs"           "atcn"           "atcvp"          "cts"           
    ## [13] "bmd"            "atbmd"

`read_secuTrial_raw` returns an object of class `secuTrialdata`, which
is basically a list. It will always contain `export_details` (which are
parsed from the HTML ExportOptions file that secuTrial generates). By
default, it will also contain all other files in the dataset. secuTrialR
automatically strips file names of dates. The new file names can be seen
via `bmd_export$export_options$data_names`.
<!-- DEDICATED ACCESSOR FUNCTION FOR DATA_NAMES? might already be implemented in the print method -->

`bmd_export` is a list, with class `secuTrialdata`. To prevent it from
printing all data to the console, a special print method returns some
useful information about the objects within `bmd_export` instead. The
information returned includes the original file name in the datafile,
it’s name in the `secuTrialdata` object, together with the number of
rows and columns and a column indicating whether the object is metadata
or not:

    bmd_export

    ## secuTrial data imported from:
    ## /Users/runner/work/_temp/Library/secuTrialR/extdata/sT_exports/BMD/s_export_CSV-
    ## xls_BMD_short_en_utf8.zip 
    ##  table nrow ncol  meta original_name
    ##     vp    1   10  TRUE        vp.xls
    ##   vpfs    1    2  TRUE      vpfs.xls
    ##     fs    1    7  TRUE        fs.xls
    ##     qs    1    7  TRUE        qs.xls
    ##     is    3    8  TRUE        is.xls
    ##    ctr    1    3  TRUE       ctr.xls
    ##     cn  113   13  TRUE        cn.xls
    ##   atcn    0    6  TRUE      atcn.xls
    ##  atcvp    0   11  TRUE     atcvp.xls
    ##    qac    0   10  TRUE       qac.xls
    ##    cts    0    8  TRUE       cts.xls
    ##    bmd  504   27 FALSE       bmd.xls
    ##  atbmd    0   28 FALSE     atbmd.xls

Individual tables can be extracted from the `bmd_export` object via
`tab <- bmd_export$tab`, where `tab` is the table of interest.
<!-- accessor function? -->

#### Variable labels

For creating tables, it is often useful to have access to variable
labels. secuTrialR supports two main methods for handling them - a named
list, or via variable attributes. The list approach works as follows.

    labs <- labels_secuTrial(bmd_export)
    # query the list with the variable name of interest
    labs[["age"]]

    ## [1] "Age"

The attribute based approach adds labels as an attribute to a variable,
which can then be accessed via `label(var)`.

    labelled <- label_secuTrial(bmd_export)
    label(labelled$bmd$age)

    ## [1] "Age"

Labels can be added to new variables or changed via

    label(labelled$bmd$age) <- "Age (years)"
    label(labelled$bmd$age)

    ## [1] "Age (years)"

Where units have been defined in the SecuTrial database, they can be
accessed or changed analogously (here, age had no unit assigned, but we
can add one).

    units(labelled$bmd$age)

    ## NULL

    units(labelled$bmd$age) <- "years"
    units(labelled$bmd$age)

    ## [1] "years"

There is a drawback to the attribute based approach - labels will not be
propagated if variables are derived and may be lost if variables are
edited.

Currently, `label_secuTrial` should be used prior to `dates_secuTrial`
or `factorize_secuTrial` so that labels and units are propagated to
factor and date variables.

#### Prepare factors

It is often useful to have categorical variables as factors (R knows how
to handle factors). secuTrialR can prepare factors easily.

    factors <- factorize_secuTrial(ctu05_raw)

This functions loops through each table of the dataset, creating new
factor variables where necessary. The new variables are the same as the
original but with `.factor` appended (i.e. a new variable called
`sex.factor` would be added to the relevant form).

    # original variable
    str(factors$ctu05baseline$gender)

    ##  int [1:17] 1 NA NA 2 1 2 1 NA NA 1 ...

    # factor
    str(factors$ctu05baseline$gender.factor)

    ##  Factor w/ 2 levels "male","female": 1 NA NA 2 1 2 1 NA NA 1 ...

    # cross tabulation
    table(original = factors$ctu05baseline$gender, factor = factors$ctu05baseline$gender.factor)

    ##         factor
    ## original male female
    ##        1    5      0
    ##        2    0      5

#### Prepare dates

Date(time)s are a very common data type. They cannot be easily used
though in their export format. This is also easily rectified in
secuTrialR:

    dates <- dates_secuTrial(ctu05_raw)

Date variables are converted to `Date` class, and datetimes are
converted to `POSIXct` class. Rather than overwriting the original
variable, new variables are added with the new class. This is a safetly
mechanism in case `NA`s are accidentally created.

    dates$ctu05baseline[c(1, 7), c("aspirin_start", "aspirin_start.date",
                                  "hiv_date", "hiv_date.datetime")]

    ##   aspirin_start aspirin_start.date     hiv_date   hiv_date.datetime
    ## 1            NA               <NA> 201903052356 2019-03-05 23:56:00
    ## 7      20060301         2006-03-01           NA                <NA>

secuTrial exports containing date variables sometimes include incomplete
dates. e.g. the day or the month may be missing. During date conversion
(i.e. `dates_secuTrial()`) `secuTrialR` currently creates `NA`s from
such incomplete date entries.

Incomplete dates are not approximated to exact dates, since this can
lead to false conclusions and biases. Users are, however, informed about
this behaviour with a `warning()`. Subsequent approximation of
incomplete dates can be manually performed.

Recommended literature on incomplete dates/date imputation:  
[Dubois and Hebert
2001](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/F50311F9FFAB56176CDDC9FFBF66F655/S1041610202008025a.pdf/imputation_of_missing_dates_of_death_or_institutionalization_for_timetoevent_analyses_in_the_canadian_study_of_health_and_aging.pdf)  
[Bowman 2006](https://www.lexjansen.com/phuse/2006/po/PO11.pdf)  

#### Recommended approach if not using `read_secuTrial`

    f <- "PATH_TO_FILE"
    d <- read_secuTrial_raw(f)
    l <- label_secuTrial(d)
    fa <- factorize_secuTrial(l)
    dat <- dates_secuTrial(fa)

    # or, if you like pipes
    library(magrittr)
    f <- "PATH_TO_FILE"
    d <- read_secuTrial_raw(f)
    dat <- d %>% 
      label_secuTrial() %>%
      factorize_secuTrial() %>%
      dates_secuTrial()

</details>

### Exploratory helpers

`secuTrialR` has a couple of functions to help get to grips with a
secuTrial data export. They are intended to be used in an exploratory
manner only.

#### as.data.frame

Working with a list can be tiresome so `secuTrialR` provides a
`as.data.frame` method to save the `data.frames` in the list to an
environment of your choice. As a demonstration, we’ll create a new
environment (`env`) and create the `data.frame`s in there. In practice,
using `.GlobalEnv` would probably be more useful.

    env <- new.env()
    ls(env)

    ## character(0)

    names(ctu05)

    ##  [1] "export_options"          "forms"                  
    ##  [3] "casenodes"               "centres"                
    ##  [5] "items"                   "questions"              
    ##  [7] "queries"                 "visitplan"              
    ##  [9] "visitplanforms"          "atcasenodes"            
    ## [11] "atcasevisitplans"        "comments"               
    ## [13] "miv"                     "cl"                     
    ## [15] "atmiv"                   "ctu05baseline"          
    ## [17] "atmnpctu05baseline"      "ctu05outcome"           
    ## [19] "atmnpctu05outcome"       "ctu05treatment"         
    ## [21] "atmnpctu05treatment"     "ctu05allmedi"           
    ## [23] "atmnpctu05allmedi"       "ctu05studyterminat"     
    ## [25] "atmnpctu05studyterminat" "ctu05ae"                
    ## [27] "atmnpctu05ae"            "ctu05sae"               
    ## [29] "atmnpctu05sae"           "emnpctu05surgeries"     
    ## [31] "atemnpctu05surgeries"    "atadverseevents"

    as.data.frame(ctu05, envir = env)
    ls(env)

    ##  [1] "atadverseevents"         "atemnpctu05surgeries"   
    ##  [3] "atmiv"                   "atmnpctu05ae"           
    ##  [5] "atmnpctu05allmedi"       "atmnpctu05baseline"     
    ##  [7] "atmnpctu05outcome"       "atmnpctu05sae"          
    ##  [9] "atmnpctu05studyterminat" "atmnpctu05treatment"    
    ## [11] "ctu05ae"                 "ctu05allmedi"           
    ## [13] "ctu05baseline"           "ctu05outcome"           
    ## [15] "ctu05sae"                "ctu05studyterminat"     
    ## [17] "ctu05treatment"          "emnpctu05surgeries"

There are also options for selecting specific forms (option
`data.frames`), changing names based on regex (options `regex` and
`rep`) and specifying whether metadata objects should be returned
(option `meta`).

#### Recruitment over time

Recruitment is an important cornerstone for every clinical trial.
`secuTrialR` allows for straigt forward visualizion of recuitment over
time for a given export file.

    # show plot
    # note that there is no line for Universitätsspital 
    # Basel because only one participant is registered for this centre
    plot_recruitment(ctu05, cex = 1.5, rm_regex = "\\(.*\\)$")

![](README_files/figure-markdown_strict/unnamed-chunk-16-1.png)

    # return the plot data
    plot_recruitment(ctu05, return_data = TRUE)

    ## [[1]]
    ##          date centre_id pat_count                      centre_name
    ## 11 2018-05-01       441         1 Universitätsspital Basel (RPACK)
    ## 1  2019-04-01       462         2           Charité Berlin (RPACK)
    ## 2  2019-04-02       462         3           Charité Berlin (RPACK)
    ## 3  2019-04-03       462         4           Charité Berlin (RPACK)
    ## 4  2019-04-04       462         5           Charité Berlin (RPACK)
    ## 5  2019-04-05       462         6           Charité Berlin (RPACK)
    ## 6  2019-04-11       461         7         Inselspital Bern (RPACK)
    ## 7  2019-04-12       461         8         Inselspital Bern (RPACK)
    ## 8  2019-04-13       461         9         Inselspital Bern (RPACK)
    ## 9  2019-04-14       461        10         Inselspital Bern (RPACK)
    ## 10 2019-04-15       461        11         Inselspital Bern (RPACK)
    ## 
    ## [[2]]
    ##         date centre_id pat_count            centre_name
    ## 1 2019-04-01       462         1 Charité Berlin (RPACK)
    ## 2 2019-04-02       462         2 Charité Berlin (RPACK)
    ## 3 2019-04-03       462         3 Charité Berlin (RPACK)
    ## 4 2019-04-04       462         4 Charité Berlin (RPACK)
    ## 5 2019-04-05       462         5 Charité Berlin (RPACK)
    ## 
    ## [[3]]
    ##          date centre_id pat_count              centre_name
    ## 6  2019-04-11       461         1 Inselspital Bern (RPACK)
    ## 7  2019-04-12       461         2 Inselspital Bern (RPACK)
    ## 8  2019-04-13       461         3 Inselspital Bern (RPACK)
    ## 9  2019-04-14       461         4 Inselspital Bern (RPACK)
    ## 10 2019-04-15       461         5 Inselspital Bern (RPACK)
    ## 
    ## [[4]]
    ##          date centre_id pat_count                      centre_name
    ## 11 2018-05-01       441         1 Universitätsspital Basel (RPACK)

Furthermore, recruitment per year and center can be returned.

    annual_recruitment(ctu05, rm_regex = "\\(.*\\)$")

    ##                     Center Total 2018 2019
    ## 1                      All    11    1   10
    ## 2           Charité Berlin     5    0    5
    ## 3         Inselspital Bern     5    0    5
    ## 4 Universitätsspital Basel     1    1    0

#### Form status summary statistics

If you are not sure about how complete the data in you export is, it may
be useful to get a quick overview of how well the forms have been
filled.

    count_summary <- form_status_summary(ctu05)
    tail(count_summary)

    ##             form_name partly_filled completely_filled empty with_warnings
    ## 5        ctu05allmedi             1                16     0             0
    ## 6       ctu05baseline             3                14     0             0
    ## 7        ctu05outcome             1                12     0             0
    ## 8            ctu05sae             0                 2     0             0
    ## 9  ctu05studyterminat             0                10     0             0
    ## 10     ctu05treatment             0                11     0             0
    ##    with_errors partly_filled.percent completely_filled.percent empty.percent
    ## 5            0            0.05882353                 0.9411765             0
    ## 6            0            0.17647059                 0.8235294             0
    ## 7            0            0.07692308                 0.9230769             0
    ## 8            0            0.00000000                 1.0000000             0
    ## 9            0            0.00000000                 1.0000000             0
    ## 10           0            0.00000000                 1.0000000             0
    ##    with_warnings.percent with_errors.percent form_count
    ## 5                      0                   0         17
    ## 6                      0                   0         17
    ## 7                      0                   0         13
    ## 8                      0                   0          2
    ## 9                      0                   0         10
    ## 10                     0                   0         11

As you can see, the majority of forms has been completeley filled. None
of the forms were saved empty, with warnings or with errors. For a more
participant id centered statistic you can perform the following.

    form_status_counts(ctu05)

This will give you a count based overview per participant id and form.
Please note that both `form_status_summary` and `form_status_counts`
only work with saved forms since unsaved form data is not available in
secuTrial exports.

#### Visit plan

secuTrialR can provide a depiction of the visit structure, although only
where the visit plan is fixed:

    vs <- visit_structure(ctu05)
    plot(vs)

<!-- PLOT METHOD DIRECTLY FOR secuTrialdata objects? -->

#### Linking different forms

Linkages amongst forms can be explored with the `links_secuTrial`
function. This relies on the `igraph` package to create a network. It is
possible to interact with the network, e.g. move nodes around in order
to read the labels better. The device ID is returned to the console, but
can be ignored. Forms are plotted in deep yellow, variables in light
blue.

    links_secuTrial(bmd_export)

![](inst/extdata/graphics/map.png)
<!-- Figure has to be generated outside of the Rmd file - resize the window and select view/"fit to screen", export it to a PDF and then convert it to a PNG -->

#### Sampling random participants

During study monitoring it is common practice to check random
participants from a study database. These participants should be
retrieved in a reproducible fashion. The below function allows this for
a loaded secuTrial data export.

    # retrieve at least 25 percent of participants recorded after March 18th 2019 
    # from the centres "Inselspital Bern" and "Charité Berlin"
    return_random_participants(ctu05, percent = 0.25, seed = 1337, date = "2019-03-18",
                               centres = c("Inselspital Bern (RPACK)", "Charité Berlin (RPACK)"))

    ## $participants
    ##          mnpaid                   centre mnpvisstartdate
    ## 2 RPACK-INS-012 Inselspital Bern (RPACK)      2019-04-12
    ## 4 RPACK-INS-014 Inselspital Bern (RPACK)      2019-04-14
    ## 5 RPACK-CBE-005   Charité Berlin (RPACK)      2019-04-05
    ## 3 RPACK-CBE-003   Charité Berlin (RPACK)      2019-04-03
    ## 
    ## $rng_config
    ## [1] "Mersenne-Twister" "Inversion"        "Rejection"

For contributors
----------------

### Testing with devtools

    # run tests
    devtools::test("secuTrialR")
    # spell check -> will contain some technical terms beyond the below list which is fine
    ignore_words <- c("AdminTool", "allforms", "casenodes", "CDMS", "codebook",
                      "codebooks", "datetime" ,"dir" ,"Hmisc" ,"igraph",
                      "labelled", "mnp", "savedforms", "secutrial", "secuTrial", 
                      "secuTrialdata", "tcltk", "tibble")
    devtools::spell_check("secuTrialR", ignore = ignore_words)

### Linting with lintr

    # lint the package -> should be clean
    library(lintr)
    lint_package("secuTrialR", linters = with_defaults(camel_case_linter = NULL,
                                                       object_usage_linter = NULL,
                                                       line_length_linter(125)))

### Building the vignette

    library(rmarkdown)
    render("vignettes/secuTrialR-package-vignette.Rmd",
           output_format=c("pdf_document"))

### Generating the README file

The README file is automatically generated on GitHub via a GitHub
action.

### Handling dependencies

Dependencies to other R packages are to be declared in the `DESCRIPTION`
file under `Imports:` and in the specific `roxygen2` documentation of
the functions relying on the dependency. It is suggested to be as
explicit as possible. i.e. Just import functions that are needed and not
entire packages.

Example to import `str_match` `str_length` `str_wrap` from the `stringr`
package (see [read\_secuTrial\_raw.R](R/read_secuTrial_raw.R)):

    #' @importFrom stringr str_match str_length str_wrap

### Preparing a release on CRAN

    # build the package archive
    R CMD build secuTrialR
    # check the archive (should return "Status: OK", no WARNINGs, no NOTEs)
    # in this example for version 0.9.0
    R CMD check secuTrialR_0.9.0.tar.gz

### Versioning and releases

The version number is made up of three digits. The first digit is
reserved for major releases which may break backwards compatibility. The
second and third digits are used for medium and minor changes
respectively. Versions released on CRAN will be tagged and saved as
releases on GitHub. The version released on CRAN is regarded as the
stable version while the master branch on GitHub is regarded as the
current development version.

#### Release checklist

Compile/Update: \* README.Rmd \* vignette \* pkgdown page \* NEWS.md

### Guidelines for contributors

Requests for new features and bug fixes should first be documented as an
[Issue](https://github.com/SwissClinicalTrialOrganisation/secuTrialR/issues)
on GitHub. Subsequently, in order to contribute to this R package you
should fork the main repository. After you have made your changes please
run the [tests](README.md#testing-with-devtools) and
[lint](README.md#linting-with-lintr) your code as indicated above.
Please also increment the version number and recompile the `README.md`
to increment the dev-version badge (requires installing the package
after editing the `DESCRIPTION` file). If all tests pass and linting
confirms that your coding style conforms you can send a pull request
(PR). Changes should also be mentioned in the `NEWS` file. The PR should
have a description to help the reviewer understand what has been
added/changed. New functionalities must be thoroughly documented, have
examples and should be accompanied by at least one
[test](tests/testthat/) to ensure long term robustness. The PR will only
be reviewed if all travis checks are successful. The person sending the
PR should not be the one merging it.

A depiction of the core functionalities for loading can be found
[here](inst/extdata/graphics/secuTrialR.png).

### Citation [![DOI](https://joss.theoj.org/papers/10.21105/joss.02816/status.svg)](https://doi.org/10.21105/joss.02816)

If you use and benefit from `secuTrialR` in your work please cite it
as:  
Wright et al., (2020). secuTrialR: Seamless interaction with clinical
trial databases in R. Journal of Open Source Software, 5(55), 2816,
<a href="https://doi.org/10.21105/joss.02816" class="uri">https://doi.org/10.21105/joss.02816</a>
# secuTrialR 1.0.9
* transferred the `Maintainer` tag to `Alan Haynes`

# secuTrialR 1.0.8
* reverted `factorize_secuTrial()` back to 1.0.3 version due to problems with lookup table factorization (#224)
* added citation (doi: 10.21105/joss.02816)

# secuTrialR 1.0.7
* moved `tcltk` and `igraph` dependency to suggested (#223)

# secuTrialR 1.0.6
* added `skip` parameter to `read_validation_overview()` (#212)

# secuTrialR 1.0.5
* fixed a bug where `subset_secuTrial()` would drop labels during the subsetting process (#203)

# secuTrialR 1.0.4
* improved import speed, specifically through changes in `dates_secutrial()` and `factorize_secutrial()` (#204)

# secuTrialR 1.0.1, 1.0.2, 1.0.3
* adjustments to handle review feedback from CRAN (#190)

# secuTrialR 1.0.0
* clarify correct options in `read_secuTrial()` failure message (#187)
* check for `project_setup` in `visit_structure()` (#181)

# secuTrialR 0.9.1
* added "Form meta data: Structure" export option information to `export_options` (#182)
* added error handling for missing structure data when running `annual_recruitment()` and `return_random_participants()` (#182)

# secuTrialR 0.9.0
* restructuring in preparation for a release on CRAN

# secuTrialR 0.8.9
* added suggestion to *NOT* export form data of hidden fields (#177)

# secuTrialR 0.8.8
* added check to make sure that specified centres are part of the export in `return_random_participants()` (#151)

# secuTrialR 0.8.7
* extended failure comment in `read_secuTrial()` to indicate that the problem could be a rectangular export file (#168)
* added "Form data of hidden fields" export option information to `export_options` (#171)
* added `return_hidden_items()` function (#172)

# secuTrialR 0.8.6
* bug fix: presence of the audit trail was incorrectly identified due to a comment in the source file of the export options (see #155, comments from @suvi-subra and @OliviaEbnerIAS)

# secuTrialR 0.8.5
* Added sorting option to visit_structure. (#152)

# secuTrialR 0.8.4
* adjusted warning message in `label_secuTrial()`
* only allow unique labels in `label_secuTrial()`
* added "Frequent warning messages" paragraph to the vignette (#156)

# secuTrialR 0.8.3
* added up-to-date vignette (#99)
* path in `print.secuTrialdata` is now wrapped at 80 characters

# secuTrialR 0.8.2
* `secutrialoptions` class is now `secuTrialoptions`.

# secuTrialR 0.8.1
* add appveyor testing, pkgdown site
* fix possible bug on windows due to regex in .prep_line_items (used in plot_recruitment) (#147)

# secuTrialR 0.8.0
* Changed license for the package from GPL-2 to MIT.

# secuTrialR 0.7.9
* The general nomenclature for a study subject will from now on be participant (pat). All variations of this
(e.g. case, patient) have been adjusted in the code and the documentation.

# secuTrialR 0.7.8
* Removed generic `plot()` function for `secuTrialdata` objects. (#139)

# secuTrialR 0.7.7
* `read_secuTrial()` and `read_secuTrial_raw()` now check if the input file exists. (#137)

# secuTrialR 0.7.6
* `factorize_secuTrial()` warning messages have been adjusted to improve trouble shooting experience. (#134, #135)

# secuTrialR 0.7.5
* `dates_secuTrial()` incomplete date warnings are now concatenated and returned as one warning per form instead of many. (#124)

# secuTrialR 0.7.4
* Fixed issue #121 on GitHub. `factorize_secuTrial()` can now handle exports which have the reset option
enabled in radio buttons.

# secuTrialR 0.7.3
* `write_secuTrial()` now allows xpt version 8 files to be written. (closes #57)

# secuTrialR 0.7.2
* `check_export_options()` function was added. It informs on deviations from suggested export options. (closes #17)
* Removed tracking of obsolete export options (`partial_date_string`, `partial_date_handling`, `unknown_date_string`).
* Added `format_info` (e.g. "CSV format for MS Excel") to `export_options`.

# secuTrialR 0.7.1
* Fixed issue #116 on GitHub.

# secuTrialR 0.7.0
* `subset_secuTrial()` function was added. It allows subsetting of secuTrialdata based on patient ID and/or study centre name.
* `get_participants()` function was added. It allows easy extractions of participant info from a secuTrialdata object.

# secuTrialR 0.6.5
* `return_random_cases()` now returns a list. The first element are the cases and the second element is the output of `RNGkind()`.

# secuTrialR 0.6.4
* New function `diff_secuTrial()` added to allow light weight comparison of the setup of two secuTrial exports.

# secuTrialR 0.6.3
* Metadata variables are now also transformed to date and datetime formats, whenever appropriate.

# secuTrialR 0.6.2
* `factorize_secuTrial()` now no longer triggers an unexpected warning when the name of a secuTrial lookuptable is equal to the name of the variable it is being used in. (PR #108)

# secuTrialR 0.6.1
* `return_random_cases()` has been added to the package. It allows to sample a random subset of cases from a secuTrial export in a reproducible fashion.

# secuTrialR 0.6.0
* `read_secuTrial_raw()` and `read_secuTrial()` no longer fail due to missing Add-ID, centre information or project setup in export data. Instead, adding of `pat_id` (no Add-ID), `centre` (no centre information) and `visit_name` (no project setup) to the data tables is now omitted if the relevant data for the operation is not available.

# secuTrialR 0.5.5
* `read_secuTrial_raw()` and `read_secuTrial()` no longer fail due to missing "Description" in export options.

# secuTrialR 0.5.4
* `dates_secuTrial()` now warns if not all dates were parsed (expected if there are incomplete dates).
* `factorize_secuTrial()` now warns if there are issues with the factorization (not expected to trigger).

# secuTrialR 0.5.2
* New function `build_secuTrial_url()` has been added. It allows users to easily compose URLs to specific secuTrial forms.

# secuTrialR 0.5.0
* The function name of `read_secuTrial_export()` has been changed to `read_secuTrial_raw()`
  to avoid confusion with `read_secuTrial()`.
 
# secuTrialR 0.4.16
* As of version 0.4.17, changes will be recorded in the NEWS file.
---
title: 'secuTrialR: Seamless interaction with clinical trial databases in R'
tags:
  - R software
  - clinical trials
  - data management
  - descriptive statistics
  - secuTrial
authors:
 - name: Patrick R. Wright
   orcid: 0000-0002-1153-0846 
   affiliation: "1, 3"
 - name: Alan G. Haynes
   orcid: 0000-0003-1374-081X
   affiliation: "2, 4"
 - name: Milica Markovic
   orcid: 0000-0002-6973-6679
   affiliation: "1, 3"
affiliations:
 - name: University Hospital Basel, Clinical Trial Unit, Basel, Switzerland
   index: 1
 - name: CTU Bern, University of Bern
   index: 2
 - name: Data Management Platform of the Swiss Clinical Trial Organisation (SCTO)
   index: 3
 - name: Statistics and Methodology Platform of the Swiss Clinical Trial Organisation (SCTO)
   index: 4
date: 14 April 2020
bibliography: paper.bib
---

# Summary

Elementary clinical trials have been conducted for hundreds of years [@meinert1986clinical]. The most famous early example
is the proof that sailors' scurvy can be cured by the consumption of citrus fruit [@lind_2014] performed by James Lind
in the 18th century. Since those initial days of clinical research, trials have significantly evolved methodically, ethically,
and technologically. While it was viable and legitimate to collect clinical trials data in unversioned
spreadsheets in the past, this is no longer true and digital clinical data management systems (CDMS) have taken over.
CDMS allow constraint-based and version-controlled data entry into a clinical trial database, which ensures traceability, 
integrity, and quality of study data.  

There is a vast market of heterogeneous CDMS solutions, each of which has individual advantages and limitations [@kuchinke_etal_2010].
One limitation can be the interaction with the data after it has been collected. Specifically, a CDMS may be
tailored for optimal data capture while, at least to some extent, disregarding ease-of-use of study data after
the conclusion of data entry. It is, however, vital that the interaction between data sources and data analysts is
fast and seamless in order to avoid loss of valuable time due to technical overhead. This point has been prominently
highlighted by the currently ongoing coronavirus pandemic [@callaway_etal_2020] in which issues have been reported
regarding the timely and complete transfer of information for the preparation of up-to-date infection
counts [@spiegel_meldeluecke; @bbc_excel]. These issues led to confusion and may have ultimately
delayed important actions. While this is a stark example, it still serves to show how severe the influence
of technical friction between digital systems can be.  

To this end we have developed the open source R statistics [@r_citation] software package `secuTrialR`, which enables
seamless interaction with data collected in the commercially available CDMS
[secuTrial](https://www.secutrial.com) (vendor [interActive Systems Berlin](https://interactive-systems.de/)).
In addition to parsing and reading the data, it performs data transformation for dates, date times, and categorical data
to reduce the data preparation overhead and to allow a swift transition into the analytical phase.
Furthermore, `secuTrialR` includes standard functionalities to
show descriptive statistics such as study recruitment or completeness of entered data per case report form
for secuTrial data exports.

# Statement of need

Due to the size and complexity of clinical trial and registry databases, technical friction during the initial interaction
with data exported from secuTrial can be expected. Our own first-hand experience revealed that this overhead can sometimes
significantly redirect scarce time and energy away from analysis and towards data management. The amount of time
spent on data management should be as small as possible. The use of `secuTrialR` leads to a pronounced reduction of time
necessary for data management, enables swift quantitative analyses through preimplemented functionalities, and most importantly,
standardizes the interaction with data exports from secuTrial, thus allowing robust and reproducible science.

While some CDMS provide APIs (e.g., REDCap [@Harris2009; @Harris2019]) or Open Database Connectivity (ODBC) connections (e.g., 
[2mt's WebSpirit](http://www.2mt-software.de)) to download data easily, using secuTrial's SOAP API involves querying
individual datapoints. This results in an extraordinarily high number of 
queries even to download a relatively small database, and high demand on servers. As such, approaches such as those 
for REDCap (e.g., the [REDCapR](https://CRAN.R-project.org/package=REDCapR) package, which can interface to REDCap's REST 
API and download all data in a single query, but does no data preparation) are not suitable for secuTrial. 
Another approach is to parse data exported manually from websites (e.g., the [ox](https://github.com/acobos/ox) package for importing [OpenClinica](https://www.openclinica.com) exports into R). This approach is used in `secuTrialR`.

# Design

All secuTrial data exports share a certain common technical structure independent of the specific database at hand.
In `secuTrialR` we make use of this information to build an S3 object of class `secuTrialdata`, which is a list, while the
data is being read into R. All downstream functions implemented in `secuTrialR` expect a `secuTrialdata` object as input
but custom analyses with other compenents of R statistics are also an option (see Figure 1).
While editing the `secuTrialdata` object is technically possible, this is not advisable.
Instead, it should be treated as raw data archive from which data can be extracted for analysis. However, if necessary,
it is possible to extract subsets of `secuTrialdata` objects with the `subset_secuTrial()` function and return
intact `secuTrialdata` objects. The individual elements of the secuTrialdata object can be accessed via regular list 
access operations or the `as.data.frame()` method, which assigns all objects to an environment of choice.

![secuTrialR information flow](secuTrialR_information_flow.png)
Figure 1: Information flow from secuTrial to R statistics and within R. Arrows indicate the 
direction from gray towards black. "..." indicates further functions working with `secuTrialdata`
objects.

# Availability

`secuTrialR` is available on [GitHub](https://github.com/SwissClinicalTrialOrganisation/secuTrialR),
[CRAN](https://cran.r-project.org/package=secuTrialR), [Anaconda Cloud](https://anaconda.org/conda-forge/r-secutrialr), and
should be functional on all major operating systems.

# Dependencies

The development of `secuTrialR` made extensive use of the `tidyverse` [@tidyverse_cit] and greatly benefited from
the `devtools` package [@devtools_cit] and `RStudio` [@rstudio_cit]. Furthermore, `tcltk` and `igraph` [@igraph_cit]
are incorporated.

# interActive Systems statement

InterActive Systems (iAS) has given permission for the open source development of this software
package but accepts no responsibility for the correctness of any functionalities within.

iAS has read and approved this manuscript.

# Acknowledgements

The authors thank Pascal Benkert, Nicole Bruni, Gilles Dutilh, Olivia Ebner, Stefanie von Felten, 
Thomas Fabbro, Inessa Kraft, Arnaud Künzi, Daniel Lengwiler, Armando Lenz, Pia Neuschwander, Henry Owusu, Hans Rock, Claudia Rokitta,
Marie Roumet, Constantin Sluka, Klaus Steigmiller, Suvitha Subramaniam, Miriam Wegmann, Laura Werlen, and Thomas Zumbrunn for ideas,
testing, and constructive feedback on the `secuTrialR` package. We also thank [Michael Sachs](https://github.com/sachsmc) 
and [Francisco Estupiñán-Romero](https://github.com/pacoramon) for kindly reviewing this manuscript and the R package and making additional 
recommendations, and [Charlotte Soneson](https://github.com/csoneson) for acting as editor.
Furthermore, the authors thank the State Secretariat of Education, Research and Innovation and the Swiss National
Science Foundation for the funding of this project and the Swiss Clinical Trial Organisation for its ongoing support.

# Conflict of interest

The authors are not employees but customers of interActive Systems (iAS). The authors therefore declare
no conflict of interest.

# References
---
name: Custom issue template
about: Describe this issue template's purpose here.
title: ''
labels: ''
assignees: ''

---


---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
# Readme for extdata

This readme file contains information about file structure and naming conventions used in the directory inst/extdata/sT_exports.

## File structure

- **sT_exports** - contains all secuTrial exports delivered with this package, including all exports used for tests and examples
  - see extdata/sT_exports/README.md for detailed description of directory contents and naming conventions
- **dictionaries** - contains .csv files with dictionary tables used for example for internationalization of the package
  - *dict_items_table.csv* - contains a table translating the secuTrial items table into different languages
  - *dict_form_status_mnpfc.csv* - contains secuTrial form completion status, validation status and data entry status codes translated into different languages
  - *dict_export_options_settings* - contains a translation table of secuTrial export option settings into different languages
  - *dict_export_options_keys* - contains a translation table of secuTrial export option keys into different languages
- **graphics** - contains all graphics delivered with this package that are used in the vignette or in the README.md
# Readme for sT_exports
This readme file contains information about file structure and naming conventions used in the directory inst/extdata/sT_exports.

## Directory structure

- **BMD** - contains secuTrial exports associated with the Bone Mineral Density (BMD) dataset
- **encodings** - contains secuTrial exports relevant for encoding
- **lnames** - contains a collection of secuTrial exports with long table names
- **snames** - contains a collection of secuTrial exports with short table names
- **exp_opt** - contains a collection of secuTrial exports with differing export options
- **change_tracking** - contains exports with differing project setups and without any case data
- **subset** - contains partial secuTrial exports based on study centre. This is needed for testing the secuTrialdata subsetting functionality.

## File naming

Following naming convention is used for all secuTrial exports contained within inst/extdata/sT_exports:

- **s_export_CSV-xls** - all exports start with this string, followed by underscore separated tags listed below
- **project tag** - short alpha numeric tag in capital letters that stands for the secuTrial database the data was extracted from
- **rt** - export in rectangular table format
- **short** - export with shortened table names
- **long** - export with full table names
- **meta** - export with duplicated form metadata in all tables
- **ref** - export with reference values in a separate table
- **no-** - export without Add-ID and without Pat-ID
- **no** - export specifically omitting certain export options
- **miss** - export with missing values
- **language tag** - export in a language. possible tags: en / de / fr / es / it / pl / unsupported
- **encoding tag** - export encoding settings e.g. "utf16", "utf8", etc.

<!-- README.md is generated from README.Rmd. Please edit that file -->



# secuTrialR

`r badger::badge_custom("dev version", as.character(packageVersion("secuTrialR")), "blue", "https://github.com/SwissClinicalTrialOrganisation/secuTrialR")` [![](https://www.r-pkg.org/badges/version/secuTrialR?color=green)](https://cran.r-project.org/package=secuTrialR)   [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/SwissClinicalTrialOrganisation/secuTrialR?branch=master&svg=true)](https://ci.appveyor.com/project/SwissClinicalTrialOrganisation/secuTrialR) [![travis](https://api.travis-ci.com/SwissClinicalTrialOrganisation/secuTrialR.svg?branch=master)](https://travis-ci.com/github/SwissClinicalTrialOrganisation/secuTrialR) [![Actions Status](https://github.com/SwissClinicalTrialOrganisation/secuTrialR/workflows/R-CMD-check/badge.svg)](https://github.com/SwissClinicalTrialOrganisation/secuTrialR/actions) [![codecov](https://codecov.io/github/SwissClinicalTrialOrganisation/secuTrialR/branch/master/graphs/badge.svg)](https://codecov.io/github/SwissClinicalTrialOrganisation/secuTrialR)

An R package to handle data from the clinical data management system (CDMS) [secuTrial](https://www.secutrial.com/en/).

## Installing from GitHub with devtools

Please note that `R versions >= 3.5` should be used to run `secuTrialR`.

```{r, eval = FALSE}
devtools::install_github("SwissClinicalTrialOrganisation/secuTrialR")
```

## Recommended export options

While the package strives to allow loading of as many types of secuTrial data exports
as possible, there are certain export options which are less likely to cause issues.
If possible it is suggested to export data which adheres to a suggested option set.
Thus, we suggest to work with exports which:
* are **zipped**
* are **English**
* have **reference values** stored **in a separate table**
* contain **Add-IDs**, **centre information**, **structure information**,  **form status**, **project setup**
* do **NOT** have the **meta data duplicated** into all tables
* are **UTF-8** encoded
* are **"CSV format"** or **"CSV format for MS Excel"**
* do **NOT** contain form **data of hidden fields**

If you use `read_secuTrial()` to read your export then it will inform you regarding deviations.

## Basic usage

An extensive applied manual/vignette is available
[here](https://github.com/SwissClinicalTrialOrganisation/secuTrialR/blob/master/vignettes/secuTrialR-package-vignette.pdf)
and probably the best place to get started.

Load the package
```{r, echo = TRUE, warning=FALSE, message=FALSE}
library(secuTrialR)
```
Load a dataset 
```{r}
export_location <- system.file("extdata", "sT_exports", "lnames",
                               "s_export_CSV-xls_CTU05_long_ref_miss_en_utf8.zip",
                               package = "secuTrialR")
ctu05 <- read_secuTrial(export_location)
```
This will load all sheets from the export into an object of class `secuTrialdata`, which is basically a list. It will always contain `export_details` (which are parsed from the HTML ExportOptions file that secuTrial generates). By default, it will also contain all other files in the dataset. secuTrialR automatically strips file names of dates. The new file names can be seen via `ctu05$export_options$data_names`. The function also adds [labels to variables](#variable-labels) and data.frames, converts [categorical variables to `factor`s](#prepare-factors) and ensures that [dates are `Date`s and date-times are `POSIXct`](#prepare-dates).
`read_secuTrial` is a wrapper for the functions described below, so it is possible to achieve more flexibility by using the individual functions (if necessary).
Individual tables can be extracted from the `ctu05` object via `tab <- ctu05$tab`, where `tab` is the table of interest.

<details><summary>Wrapped functions</summary>


#### Load the dataset
```{r}
# prepare path to example export
export_location <- system.file("extdata", "sT_exports", "BMD",
                               "s_export_CSV-xls_BMD_short_en_utf8.zip",
                               package = "secuTrialR")
# load all export data
bmd_export <- read_secuTrial_raw(data_dir = export_location)

# load a second dataset
export_location <- system.file("extdata", "sT_exports", "lnames",
                               "s_export_CSV-xls_CTU05_long_ref_miss_en_utf8.zip",
                               package = "secuTrialR")
ctu05_raw <- read_secuTrial_raw(export_location)

# View names of the bmd_export object
names(bmd_export)
```

`read_secuTrial_raw` returns an object of class `secuTrialdata`, which is basically a list. It will always contain `export_details` (which are parsed from the HTML ExportOptions file that secuTrial generates). By default, it will also contain all other files in the dataset. secuTrialR automatically strips file names of dates. The new file names can be seen via `bmd_export$export_options$data_names`.
<!-- DEDICATED ACCESSOR FUNCTION FOR DATA_NAMES? might already be implemented in the print method -->

`bmd_export` is a list, with class `secuTrialdata`. To prevent it from printing all data to the console, a special print method returns some useful information about the objects within `bmd_export` instead. The information returned includes the original file name in the datafile, it's name in the `secuTrialdata` object, together with the number of rows and columns and a column indicating whether the object is metadata or not:
```{r}
bmd_export
```

Individual tables can be extracted from the `bmd_export` object via `tab <- bmd_export$tab`, where `tab` is the table of interest.
<!-- accessor function? -->


#### Variable labels
For creating tables, it is often useful to have access to variable labels. secuTrialR supports two main methods for handling them - a named list, or via variable attributes. The list approach works as follows.
```{r}
labs <- labels_secuTrial(bmd_export)
# query the list with the variable name of interest
labs[["age"]]

```

The attribute based approach adds labels as an attribute to a variable, which can then be accessed via `label(var)`.
```{r}
labelled <- label_secuTrial(bmd_export)
label(labelled$bmd$age)
```
Labels can be added to new variables or changed via 
```{r}
label(labelled$bmd$age) <- "Age (years)"
label(labelled$bmd$age)
```
Where units have been defined in the SecuTrial database, they can be accessed or changed analogously (here, age had no unit assigned, but we can add one).
```{r}
units(labelled$bmd$age)
units(labelled$bmd$age) <- "years"
units(labelled$bmd$age)
```
There is a drawback to the attribute based approach - labels will not be propagated if variables are derived and may be lost if variables are edited.

Currently, `label_secuTrial` should be used prior to `dates_secuTrial` or `factorize_secuTrial` so that labels and units are propagated to factor and date variables.


 
#### Prepare factors
It is often useful to have categorical variables as factors (R knows how to handle factors). secuTrialR can prepare factors easily.
```{r, error=TRUE}
factors <- factorize_secuTrial(ctu05_raw)
```
This functions loops through each table of the dataset, creating new factor variables where necessary. The new variables are the same as the original but with `.factor` appended (i.e. a new variable called `sex.factor` would be added to the relevant form).

```{r}
# original variable
str(factors$ctu05baseline$gender)
# factor
str(factors$ctu05baseline$gender.factor)
# cross tabulation
table(original = factors$ctu05baseline$gender, factor = factors$ctu05baseline$gender.factor)
```


#### Prepare dates
Date(time)s are a very common data type. They cannot be easily used though in their export format. This is also easily rectified in secuTrialR:


```{r}
dates <- dates_secuTrial(ctu05_raw)
```

Date variables are converted to `Date` class, and datetimes are converted to `POSIXct` class. Rather than overwriting the original variable, new variables are added with the new class. This is a safetly mechanism in case `NA`s are accidentally created.

```{r}
dates$ctu05baseline[c(1, 7), c("aspirin_start", "aspirin_start.date",
                              "hiv_date", "hiv_date.datetime")]
```

secuTrial exports containing date variables sometimes include incomplete dates. e.g. the day or the month may be missing.
During date conversion (i.e. `dates_secuTrial()`) `secuTrialR` currently creates `NA`s from such incomplete date entries.

Incomplete dates are not approximated to exact dates, since this can lead to false conclusions and biases.
Users are, however, informed about this behaviour with a `warning()`. Subsequent approximation of incomplete dates can be manually performed.

Recommended literature on incomplete dates/date imputation:\
[Dubois and Hebert 2001](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/F50311F9FFAB56176CDDC9FFBF66F655/S1041610202008025a.pdf/imputation_of_missing_dates_of_death_or_institutionalization_for_timetoevent_analyses_in_the_canadian_study_of_health_and_aging.pdf) \
[Bowman 2006](https://www.lexjansen.com/phuse/2006/po/PO11.pdf) \


#### Recommended approach if not using `read_secuTrial`

```{r, eval=FALSE}
f <- "PATH_TO_FILE"
d <- read_secuTrial_raw(f)
l <- label_secuTrial(d)
fa <- factorize_secuTrial(l)
dat <- dates_secuTrial(fa)

# or, if you like pipes
library(magrittr)
f <- "PATH_TO_FILE"
d <- read_secuTrial_raw(f)
dat <- d %>% 
  label_secuTrial() %>%
  factorize_secuTrial() %>%
  dates_secuTrial()
```

</details>

### Exploratory helpers
`secuTrialR` has a couple of functions to help get to grips with a secuTrial data export. They are intended to be used in an exploratory manner only.

#### as.data.frame
Working with a list can be tiresome so `secuTrialR` provides a `as.data.frame` method to save the `data.frames` in the list to an environment of your choice. 
As a demonstration, we'll create a new environment (`env`) and create the `data.frame`s in there. In practice, using `.GlobalEnv` would probably be more useful.

```{r}
env <- new.env()
ls(env)
names(ctu05)
as.data.frame(ctu05, envir = env)
ls(env)
```

There are also options for selecting specific forms (option `data.frames`), changing names based on regex (options `regex` and `rep`) and specifying whether metadata objects should be returned (option `meta`).


#### Recruitment over time
Recruitment is an important cornerstone for every clinical trial. `secuTrialR` allows for straigt forward visualizion of recuitment
over time for a given export file.

```{r, eval = TRUE}
# show plot
# note that there is no line for Universitätsspital 
# Basel because only one participant is registered for this centre
plot_recruitment(ctu05, cex = 1.5, rm_regex = "\\(.*\\)$")
# return the plot data
plot_recruitment(ctu05, return_data = TRUE)
```

Furthermore, recruitment per year and center can be returned.

```{r, eval = TRUE}
annual_recruitment(ctu05, rm_regex = "\\(.*\\)$")
```


#### Form status summary statistics
If you are not sure about how complete the data in you export is, it may be useful to get a quick overview of how well the forms
have been filled.

```{r, eval = TRUE}
count_summary <- form_status_summary(ctu05)
tail(count_summary)
```

As you can see, the majority of forms has been completeley filled. None of the forms were saved empty, with warnings or with errors.
For a more participant id centered statistic you can perform the following.

```{r, eval = FALSE}
form_status_counts(ctu05)
```

This will give you a count based overview per participant id and form. Please note that both `form_status_summary` 
and `form_status_counts` only work with saved forms since unsaved form data is not available in secuTrial exports.

#### Visit plan
secuTrialR can provide a depiction of the visit structure, although only where the visit plan is fixed:
```{r, eval = FALSE}
vs <- visit_structure(ctu05)
plot(vs)
```
<!-- PLOT METHOD DIRECTLY FOR secuTrialdata objects? -->


#### Linking different forms

Linkages amongst forms can be explored with the `links_secuTrial` function. This relies on the `igraph` package to create a network. It is possible to interact with the network, e.g. move nodes around in order to read the labels better. The device ID is returned to the console, but can be ignored. Forms are plotted in deep yellow, variables in light blue.

```{r, eval=FALSE}
links_secuTrial(bmd_export)
```
![](inst/extdata/graphics/map.png)
<!-- Figure has to be generated outside of the Rmd file - resize the window and select view/"fit to screen", export it to a PDF and then convert it to a PNG -->

#### Sampling random participants

During study monitoring it is common practice to check random participants from a study database. These
participants should be retrieved in a reproducible fashion. The below function allows this for a loaded 
secuTrial data export.

```{r}
# retrieve at least 25 percent of participants recorded after March 18th 2019 
# from the centres "Inselspital Bern" and "Charité Berlin"
return_random_participants(ctu05, percent = 0.25, seed = 1337, date = "2019-03-18",
                           centres = c("Inselspital Bern (RPACK)", "Charité Berlin (RPACK)"))
```

## For contributors
### Testing with devtools

```{r, eval = FALSE}
# run tests
devtools::test("secuTrialR")
# spell check -> will contain some technical terms beyond the below list which is fine
ignore_words <- c("AdminTool", "allforms", "casenodes", "CDMS", "codebook",
                  "codebooks", "datetime" ,"dir" ,"Hmisc" ,"igraph",
                  "labelled", "mnp", "savedforms", "secutrial", "secuTrial", 
                  "secuTrialdata", "tcltk", "tibble")
devtools::spell_check("secuTrialR", ignore = ignore_words)
```

### Linting with lintr


```{r, eval = FALSE}
# lint the package -> should be clean
library(lintr)
lint_package("secuTrialR", linters = with_defaults(camel_case_linter = NULL,
                                                   object_usage_linter = NULL,
                                                   line_length_linter(125)))
```

### Building the vignette
```{r, eval = FALSE}
library(rmarkdown)
render("vignettes/secuTrialR-package-vignette.Rmd",
       output_format=c("pdf_document"))
```

### Generating the README file

The README file is automatically generated on GitHub via a GitHub action.

### Handling dependencies

Dependencies to other R packages are to be declared in the `DESCRIPTION` file under `Imports:` and in
the specific `roxygen2` documentation of the functions relying on the dependency. It is suggested to
be as explicit as possible. i.e. Just import functions that are needed and not entire packages.

Example to import `str_match` `str_length` `str_wrap` from the `stringr` package (see [read_secuTrial_raw.R](R/read_secuTrial_raw.R)):
```{r, eval = FALSE}
#' @importFrom stringr str_match str_length str_wrap
```

### Preparing a release on CRAN

```bash
# build the package archive
R CMD build secuTrialR
# check the archive (should return "Status: OK", no WARNINGs, no NOTEs)
# in this example for version 0.9.0
R CMD check secuTrialR_0.9.0.tar.gz
```

### Versioning and releases

The version number is made up of three digits. The first digit
is reserved for major releases which may break backwards compatibility.
The second and third digits are used for medium and minor changes respectively.
Versions released on CRAN will be tagged and saved as releases on GitHub.
The version released on CRAN is regarded as the stable version while
the master branch on GitHub is regarded as the current development version.

#### Release checklist

Compile/Update:
* README.Rmd
* vignette
* pkgdown page
* NEWS.md

### Guidelines for contributors

Requests for new features and bug fixes should first be documented as an [Issue](https://github.com/SwissClinicalTrialOrganisation/secuTrialR/issues) on GitHub.
Subsequently, in order to contribute to this R package you should fork the main repository.
After you have made your changes please run the 
[tests](README.md#testing-with-devtools)
and 
[lint](README.md#linting-with-lintr) your code as 
indicated above. Please also increment the version number and recompile the `README.md` to increment the dev-version badge (requires installing the package after editing the `DESCRIPTION` file). If all tests pass and linting confirms that your 
coding style conforms you can send a pull request (PR). Changes should also be mentioned in the `NEWS` file.
The PR should have a description to help the reviewer understand what has been 
added/changed. New functionalities must be thoroughly documented, have examples 
and should be accompanied by at least one [test](tests/testthat/) to ensure long term 
robustness. The PR will only be reviewed if all travis checks are successful. 
The person sending the PR should not be the one merging it.

A depiction of the core functionalities for loading can be found [here](inst/extdata/graphics/secuTrialR.png).

### Citation  [![DOI](https://joss.theoj.org/papers/10.21105/joss.02816/status.svg)](https://doi.org/10.21105/joss.02816)

If you use and benefit from `secuTrialR` in your work please cite it as:  
Wright et al., (2020). secuTrialR: Seamless interaction with clinical trial databases in R.
Journal of Open Source Software, 5(55), 2816, https://doi.org/10.21105/joss.02816
---
title: "secuTrialR - a walkthrough"
author: "Patrick R. Wright, Milica Markovic, Alan G. Haynes"
date: "`r Sys.Date()`"
output: rmarkdown::pdf_document
toc: true
vignette: >
  %\VignetteIndexEntry{secuTrialR-package-vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

\newpage

# Introduction

> “If I had just five minutes to chop down a tree I would spend the first two and a half minutes sharpening my axe.”

> `r tufte::quote_footer('--- An anonymous woodsman')`

\vspace{13pt}

This R package provides functions for handling data from the clinical data management system (CDMS) 
[\textcolor{blue}{secuTrial}](https://www.secutrial.com/en/).
The most important components are related to reading data exports from
secuTrial into R. In brief, the package aims to enable swift execution of
repetitive tasks in order to allow spending more time on the unique aspects
of a dataset. It is developed and maintained by the Swiss Clinical Trial Organisation
([\textcolor{blue}{SCTO}](https://www.scto.ch/en/news.html)).

If you are still challenged by more basic operations in R we suggest reading
[\textcolor{blue}{Hands-On Programming with R}](https://rstudio-education.github.io/hopr/),
which serves as an excellent introduction to the basic concepts of R.

This vignette will teach you how to use the `secuTrialR` package and you will likely
learn quite a bit about secuTrial exports in general along the way.
Throughout the `secuTrialR` package and within this vignette we refer to
patients, cases, subjects etc. enrolled in a secuTrial database as participants.

```{r, include = FALSE}
# needed so that the as.data.frame part of the vignette
# does not need a restart of the session everytime the
# vignette is built
#rm(list = ls()) # removed this at the request of the CRAN submission
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# Install

Please note that `R versions >= 3.5` should be used to run `secuTrialR`.
\vspace{5pt}

**Stable release from CRAN**
\vspace{5pt}
```{r, eval = FALSE}
install.packages("secuTrialR", dependencies = TRUE)
```

**Development release from GitHub with `devtools`**

For this you will need to have `devtools` installed.
If you are working on Windows and would like to install with `devtools` you
will likely need to install
[\textcolor{blue}{Rtools}](https://cran.r-project.org/bin/windows/Rtools/).
Installing everything, including the dependencies, from scratch may take a while (20-30 minutes).

\vspace{5pt}

```{r, eval = FALSE}
# install
devtools::install_github("SwissClinicalTrialOrganisation/secuTrialR")
```

```{r}
# load silently
suppressMessages(library(secuTrialR))
# show secuTrialR version
packageVersion("secuTrialR")
```

# The CTU05 dataset

Before we continue with the functionalities let's briefly talk about the test data which is
delivered as a part of the package. We refer to it as the `CTU05` (clinical trial unit project 05)
data. This dataset has been fabricated for demonstration purposes only and is not
real clinical data. Principally it is made up of eight forms. These are called "surgeries",
"baseline", "outcome", "treatment", "allmedi", "studyterminat", "ae" and "sae".
You will see these names again later when the data has been read into R.
The project setup includes most data types implementable in secuTrial. It is, however,
not exhaustive. Since the data is delivered with the installation of the `secuTrialR`
package we can point to it via the `system.file()` function.
\vspace{5pt}

```{r}
ctu05_data_location <- system.file("extdata", "sT_exports", "exp_opt",
                                   "s_export_CSV-xls_CTU05_all_info.zip",
                                   package = "secuTrialR")
```

If you work on your own datasets you can specify a path as a regular character string without
using `system.file()`.

# secuTrial export options

Prior to reading your data into R you need to export it with the secuTrial ExportSearchTool.
We suggest exporting **non-rectangular**, **zipped**, **English** data with 
**reference values** stored **in a separate table** including 
**Add-IDs**, **centre information**, **structure information**,  **form status**,
**project setup**, **without duplicated meta data** and **without form data of hidden fields**.
Furthermore, it is important to use
**"CSV format"**/**"CSV format for MS Excel"** and suggested to select **UTF-8** encoding.
Most of these options are truly optional and reading your data should
work even with differences from the above specifications.

A description of how data can be exported from secuTrial can be found
[\textcolor{blue}{here}](https://swissclinicaltrialorganisation.github.io/secuTrial_recipes/export_data/).
This description includes screenshots of the export options configuration interface.

# Reading a secuTrial data export into R

There is one principle function to read your data (i.e. `read_secuTrial()`).
Below you can see it in action with the `CTU05` dataset.
\vspace{5pt}

```{r}
ctu05_data <- read_secuTrial(data_dir = ctu05_data_location)
```

If the "Read export successfully." message appears your data was correctly read.
In this example you are also warned that hidden data fields are in the export which is a deviation
from the suggested export option configuraion.

# The `secuTrialdata` object

If you inspect the `class()` of `ctu05_data` you will find that it is a `secuTrialdata` object.
\vspace{5pt}

```{r}
class(ctu05_data)
```

Really this is only a `list` containing all the information from your secuTrial export.
\vspace{5pt}

```{r}
typeof(ctu05_data)
```

\newpage

## The data tables in the `secuTrialdata` object

We have implemented a custom variation of the `print()` function for `secuTrialdata` objects.

```{r}
print(ctu05_data)
```

It shows you where the export archive of your `secuTrialdata` object is located, tells you which
data tables (i.e. `table`) it contains, what the source files (i.e. `original_name`) are and
specifies each table's dimensions (i.e. `ncol`, `nrow`).

By now you have possibly realized that all the forms specified earlier
(i.e. "surgeries", "baseline", "outcome", "treatment", "allmedi", "studyterminat", "ae"
and "sae") are present, but also that there are many tables that do not correspond to the
previously introduced forms.

The majority of the unfamiliar tables are tagged as `TRUE` in the `meta` column.
This means that they are metadata tables. Their names and data structures
are fixed in secuTrial exports. In the following we will briefly explain
which information the most relevant meta tables contain.

* `vp` - visitplan definition
* `vpfs` - visitplan form linkage
* `fs` - forms information
* `qs` - questions
* `is` - items i.e. variable definitions
* `ctr` - centre information
* `cn` - casenodes i.e. table of entered study participants
* `cl` - information how the data in the variables is coded

Furthermore, there is a set of tables whose names start with "at". These are
audit trail tables. They are only relevant if you need to investigate changes
in the data over time. For example certain values may be corrected (i.e. changed)
due to findings during monitoring visits at study centres.
Last but not least you may have also realized that the "surgeries" table is
called `esurgeries`. This is because it is a so-called repetition form.
Repetition forms are labelled with a leading "e" and are implemented
as subforms in other forms. In this case, `esurgeries` is a subform in
`baseline` and the linkage is defined by the `mnpdocid` column in both tables.
If this sounds cryptic to you we suggest you talk so someone who has implemented
a database in secuTrial and let them explain it with a specific example.
It is pretty straight forward when you look at a concrete implementation.

## Accessing the tables and values

Since the `secuTrialdata` object is a `list` and the data tables within this `list` are
`data.frame`s you can simply access the tables using `$`. Let's say you would like to have
a look at the placebo to verum ratio in your `treatment` data or what types of other
medication were entered in `allmedi`.
\vspace{5pt}

```{r}
table(ctu05_data$treatment$rando_treatment)
table(ctu05_data$allmedi$med_product)
```

## Data transformations

During the loading process, coded categorical data is transformed. For example
the `gender` variable in the `baseline` form is categorical. The raw data is
accessible via `gender` and the transformed version of the data is added during
the reading process and becomes accessible via `gender.factor`. Thus, data is not overwritten
but added with the `.factor` extension. If there are issues during factorization
a `warning()` will inform you of this.
\vspace{5pt}

```{r}
# raw gender data
ctu05_data$baseline$gender

# transformed gender data
ctu05_data$baseline$gender.factor

# raw more meds
ctu05_data$allmedi$no_more_meds

# transformed more meds
ctu05_data$allmedi$no_more_meds.factor
```

Note that descriptive labels have also been automatically added to the data.
\vspace{5pt}

```{r}
label(ctu05_data$allmedi$no_more_meds.factor)
label(ctu05_data$baseline$gender.factor)
label(ctu05_data$esurgeries$surgery_organ.factor)
```

Datetime data is also transformed and similarly to the factorization
process the names are concatenated with
`.date` or `.datetime`.
\vspace{5pt}

```{r}
# raw
ctu05_data$baseline$visit_date

# processed
ctu05_data$baseline$visit_date.date

# raw only head
head(ctu05_data$baseline$hiv_date)

# processed only head
head(ctu05_data$baseline$hiv_date.datetime)

# classes
class(ctu05_data$baseline$visit_date)
class(ctu05_data$baseline$visit_date.date)
class(ctu05_data$baseline$hiv_date)
class(ctu05_data$baseline$hiv_date.datetime)

```

Depending on the setup, incomplete dates can be valid entries in a secuTrial database.
Thus they will also occasionally appear in your exports. The datetime conversion
does not work in these cases and `NA`s are created. If this happens, `secuTrialR` will
warn you accordingly and you should have a closer look into the affected
datetime variables and whether you would like to perform so-called date imputation.

## Export options

The `secuTrialdata` object also contains information on the export options.
\vspace{5pt}

```{r}
ctu05_data$export_options
```

`export_options` itself is a `list`. If you are interested in more information
than is printed you can also access it. Let's assume you would like to know
the `project_name` and `encoding`.
\vspace{5pt}

```{r}
ctu05_data$export_options$project_name
ctu05_data$export_options$encoding
```

Much more information is stored in the elements of `export_options`.
The names of the elements should be descriptive enough to infer the contents.
\vspace{5pt}

```{r}
names(ctu05_data$export_options)
```

# Generic functions for `secuTrialdata` objects

Now that you understand the `secuTrialdata` object we will show you some generic functions
you can use on objects of this class.

## Show the study participants

First off you may be interested in a table of participants.
\vspace{5pt}

```{r}
get_participants(ctu05_data)
```

Please note that the `mnpaid` column in this table corresponds to the `pat_id`
column in other tables.

\newpage

## Recruitment over time

You can extract information about participant recruitment per centre and year by applying
`annual_recruitment()` on a `secuTrialdata` object.
\vspace{5pt}

```{r}
annual_recruitment(ctu05_data)
```

Since the centre names often have a systematic addition (e.g. (RPACK)) we have
enabled the option to remove certain parts of the centre descriptions via
regular expressions (i.e. `rm_regex` argument). In this case the regular expression
removes trailing parentheses and everything they enclose.
\vspace{5pt}

```{r}
annual_recruitment(ctu05_data, rm_regex = "\\(.*\\)$")
```

It is also possible to plot the recruitment over time.
\vspace{5pt}

```{r, fig.height = 3.6, fig.width = 8}
plot_recruitment(ctu05_data, cex = 1.2, rm_regex = "\\(.*\\)$")
```

\newpage

## Visit plan visualization

`secuTrialR` can provide a depiction of the visit structure, although only where
the visit plan is fixed. Black rectangles in the grid represent a form to be filled (x)
during one of the visits (y).
\vspace{5pt}

```{r, fig.height = 3.9, fig.width = 3.9}
vs <- visit_structure(ctu05_data)
plot(vs)
```

## Completeness of forms

If you are not sure about how complete the data in your export is, it may be useful to
get a quick overview of how well the forms have been filled. The below table shows
both absolute and relative numbers for a few forms.
\vspace{5pt}

```{r}
fss <- form_status_summary(ctu05_data)
tail(fss, n = 5)
```

Please note that a form is only complete if all required fields have been filled.
Thus, a whole study may have 99% completeness on variable basis while showing
0% completeness on form basis. It is currently not technically possible to assess
completeness on variable basis in a generic way. Hence, high completeness on form
basis implies high completeness on variable basis but **NOT** vice versa.

If you would rather retrieve information on form completeness for each participant individually
you can perform the following.
\vspace{5pt}

```{r}
fsc <- form_status_counts(ctu05_data)
# show the top
head(fsc)
```

## Form linkage

Linkages amongst forms can be explored with the `links_secuTrial()` function.
This relies on the `igraph` package to create a network. It is possible to
interact with the network, e.g. move nodes around in order to read the labels better.
The R graphics device ID is returned to the console, but can be ignored. Forms are plotted in
deep yellow, variables in light blue.
\vspace{5pt}

```{r, eval = FALSE}
links_secuTrial(ctu05_data)
```

The output can not be shown within this vignette but you should give it a try.
Please note that the linkage plot is likely most useful **without** the audit trail
data in the export.

\newpage

## Sampling random participants

During study monitoring it is common practice to check random participants
from a study database. These participants should be retrieved in a reproducible
fashion, which can be achieved by setting a so-called seed.
The below function allows reproducible retrieval for a loaded secuTrial data export.
\vspace{5pt}

```{r}
# randomly retrieve at least 25 percent of participants recorded after March 18th 2019
# from the centres "Inselspital Bern" and "Charité Berlin"
return_random_participants(ctu05_data,
                           percent = 0.25,
                           seed = 1337,
                           date = "2019-03-18",
                           centres = c("Inselspital Bern (RPACK)",
                                       "Charité Berlin (RPACK)"))
```

Please note that earlier R versions may return different results because
there is a different `rng_config` (i.e. `RNGkind()`).
For this reason we have added the `rng_config` to the output.

## Retrieve score variables

secuTrial allows implementing calculated variables (i.e. scores). Data is
not directly entered into these variables but rather calculated automatically.
Scores are defined by a set of rules and use the data in other variables
as basis. For example the age of a study participant at data entry can be
calculated as the difference between the participant's birthday and the day
of data entry.

It is advisable to recalculate or validate score variable data before data analysis.
A rule of thumb: The more complex a score is and the more data from different 
forms is necessary for its calculation the more likely its value should be recalculated.
The below function will allow you to detect which variables this concerns.
\vspace{5pt}

```{r}
return_scores(ctu05_data)
```

## Retrieve hidden variables

Sometimes, during a study, certain fields may be hidden because data should
no longer be entered into them. If this is the case and the data of these
fields is part of your export is likely good to know about it. In this case
nothing is hidden.
\vspace{5pt}

```{r}
return_hidden_items(ctu05_data)
```


## Finding changes/differences in project setup implementations

In ongoing studies it is possible that changes to the secuTrial data entry
interface (i.e. the electronic case report forms) are made. Sometimes these
changes may call for adjustments in analysis code. It is considered good practice
to run `diff_secuTrial()` on the last export and the current export of a project
to at least make yourself aware of potential changes in the setup. If there are
differences, the results of this function should be interpreted as a first indicator
since they may not cover all alterations. Information is returned on forms and
variables. A detailed list of changes can be produced in the secuTrial
FormBuilder with "Compare project setup".

For the below `diff_secuTrial()` showcase we emulated a changed setup of CTU05
by copying the setup and importing it in the FormBuilder as a new secuTrial
project (CTU06). From this, we created a data export (v1) and then made a few
minor changes and exported again (v2). If this sounds confusing, never mind.
CTU06 v1 is simply a copy of CTU05. CTU06 v2 is a slighly altered version of CTU06 v1.
\vspace{5pt}

``` {r}
ctu06_v1 <- read_secuTrial(system.file("extdata", "sT_exports", "change_tracking",
                                       "s_export_CSV-xls_CTU06_version1.zip",
                                       package = "secuTrialR"))

ctu06_v2 <- read_secuTrial(system.file("extdata", "sT_exports", "change_tracking",
                                       "s_export_CSV-xls_CTU06_version2.zip",
                                       package = "secuTrialR"))

diff_secuTrial(ctu06_v1, ctu06_v2)

```

As you can see `ctu06_v2` contains the two additional forms `mnpctu06anewform` and
`mnpctu06anothernewform` and the two additional variables `new_item_in_fu` and
`new_item_in_new_form`.

## Conversion to SPSS, STATA, SAS

Given that you are working with R it is unlikely that you need such conversions for
yourself. However, collaborators may ask for data which is readily importable into SPSS,
STATA or SAS. For this you can use `write_secuTrial()`.  

Since this has not been heavily tested or used there may be issues and
you might prefer doing this manually with the `haven` package. One particular 
sticking point is the length of variable names - R is not restrictive in this 
respect, but other software can be. `secuTrialR` does not truncate names, prefering 
to leave this to the user, which can cause `write_secuTrial()` to fail with an error.
\vspace{5pt}

```{r, eval = FALSE}
# retrieve path to a temporary directory
tdir <- tempdir()
# write spss
write_secuTrial(ctu05_data, format = "sav", path = tdir)
```

## Subsetting `secuTrialdata`

In some cases it may be useful to subset your `secuTrialdata` object.
For example if you have cohort data and would like to supply a subset of
the data for a retrospective study. We have implemented this
option with `subset_secuTrial()`. It will truncate your `secuTrialdata` object
and return a new `secuTrialdata` object which is a subset of the original data.
It is possible to subset by including or excluding specific participant
ids or centres.
\vspace{5pt}

```{r}
# initialize some subset identifiers
participants <- c("RPACK-INS-011", "RPACK-INS-014", "RPACK-INS-015")
centres <- c("Inselspital Bern (RPACK)", "Universitätsspital Basel (RPACK)")

# exclude Bern and Basel
ctu05_data_berlin <- subset_secuTrial(ctu05_data, centre = centres, exclude = TRUE)
# exclude Berlin
ctu05_data_bern_basel <- subset_secuTrial(ctu05_data, centre = centres)
# keep only subset of participants
ctu05_data_pids <- subset_secuTrial(ctu05_data, participant = participants)

class(ctu05_data_berlin)
class(ctu05_data_bern_basel)
class(ctu05_data_pids)
```

If you subset based on centres all traces of deleted centres will be removed.
If you remove based on participant ids all traces of deleted participants will be removed.
\vspace{5pt}

```{r}
# only Berlin remains
ctu05_data_berlin$ctr

# all centres remain even though all three participant ids are from Bern
ctu05_data_pids$ctr
```

\newpage

Since the truncated object's class remains unchanged (i.e. `secuTrialdata`) you can
still use the generic functions on it.
Let's say you would only like to look at the recruitment plot for Bern alone.
\vspace{5pt}

```{r, fig.height = 3.8, fig.width = 8}
# keep only Bern
ctu05_data_bern <- subset_secuTrial(ctu05_data, centre = "Inselspital Bern (RPACK)")
# plot
plot_recruitment(ctu05_data_bern)
```

... or Bern and Berlin.
\vspace{5pt}

```{r, fig.height = 3.8, fig.width = 8}
# keep only Bern and Berlin
ctu05_data_bern_berlin <- subset_secuTrial(ctu05_data,
                                           centre = c("Inselspital Bern (RPACK)",
                                                      "Charité Berlin (RPACK)"))
# plot
plot_recruitment(ctu05_data_bern_berlin)
```

# Building URLs to your secuTrial server

If you are creating reports in which you would like to directly
link to specific pages of your secuTrial DataCapture you can
use `build_secuTrial_url`. If you are no expert regarding the
secuTrial server architecture you would like to build links for, you should
talk to the server admin or consult the `build_secuTrial_url` help page.
They will be able to guide you regarding the information for the `server`,
`instance`, `customer` and `project` parameters. The `docid`, however, is
included in your export data in the non-meta data tables of the
`secuTrialdata` object and can be found in the `mnpdocid` columns.

```{r}
head(ctu05_data$treatment$mnpdocid)
head(ctu05_data$baseline$mnpdocid)
```

To demsonstrate `build_secuTrial_url` we will use imaginary data
for the `server`, `instance`, `customer` and `project` parameters.
The real counterparts on your server will likely look structurally
similar.

```{r}
server <- "server.secutrial.com"
instance <- "ST21-setup-DataCapture"
customer <- "TES"
project <- "7036"

# make three links with the first three baseline docids
bl_docids <- head(ctu05_data$baseline$mnpdocid, n = 3)
links <- build_secuTrial_url(server, instance, customer,
                             project, bl_docids)
```
**These are the links**:  
`r links[1]`  

`r links[2]`  

`r links[3]`  

Of course they are dead ends but maybe you can use them to make
out the arguments for your server.

# The `as.data.frame` function

This vignette has been working with the `secuTrialdata` object, which is of type `list`.
For some users, working with a `list` can be tiresome so `secuTrialR` provides an
`as.data.frame()` method to save the `data.frame`s in the `secuTrialdata` object to
an environment of your choice.

As an example, we will create an environment called `env` and check that it's empty
before running `as.data.frame()`...
\vspace{5pt}

```{r}
env <- new.env()
ls(env)
```

```{r}
# add files to env
as.data.frame(ctu05_data, envir = env)
```

... and afterwards.
\vspace{5pt}
```{r}
ls(env)
```

Substituting `env` with `.GlobalEnv` instead would also be an option and would make the data.frames immediately accessible without having to refer to an environment.

# Frequent warning messages

Certain `warning` messages can occur quite frequently when running `read_secuTrial()`.
Some of them may call for deliberate action and thus it is important to understand them.
We briefly mentioned some of them earlier in this document but will now more closely
explain how they can be interpreted.

Please note that `warning` messages may "pile up" depending on the export you are reading.
For example this may happen if there are many date variables with incomplete data. This is
no reason for concern. We suggest that you read them and interpret them based on the
explanations below. We use `a_form_name` and `a_variable_name` as place holders in the examples.
If in doubt you can always work with the raw data because it is never overwritten.

## Dates

The below warning tells you that some data in a date variable could not be converted
during the process of date conversion (i.e. `dates_secuTrial()`). This ususally occurs
if incomplete date entries are present. Since the raw data is not overwritten but rather a
`variable_name.date` or `variable_name.datetime` column are added to the dataset you can specifically see which
values could not be converted because the raw data will contain data while the corresponding `.date`/`.datetime`
entires will be `NA`. The `warning` also indicates where to look. The dummy example below indicates to look
at the variable `a_variable_name` in form `a_form_name`.

```{r, echo = FALSE, results = TRUE}

# incomplete dates
warning(
"In dates_secuTrial.data.frame(tmp, datevars, timevars, dateformat,  :
Not all dates were converted for
variable: 'a_variable_name'
in form: 'a_form_name'
This is likely due to incomplete date entries."
)

```

## Factors

In some cases secuTrial allows differently coded data to be decoded to the same target value for the
same variable. For instance this can happen if hierarchical lookuptables have been implemented
in the database. Because this interferes with the factorization (i.e. `factorize_secuTrial()`)
we add the code to the duplicate decoded value and return the below message to make you aware.

If you run into this `warning` message we suggest running the `table()` function on the
variable in question. This will likely clarify the above explanation.

```{r, echo = FALSE, results = TRUE}
# duplicate factors
warning(
"In factorize_secuTrial.data.frame(curr_form_data, cl = object$cl,  :
Duplicate values found during factorization of a_variable_name")
```

## Labels

Sometimes the labels of variables in a secuTrial database implementation may be changed
after release of the database. In these cases all labels (current and previous versions) are
added to the `label` attribute during labelling (i.e. `label_secuTrial()`) and the below `warning`
is triggered. It indicates which variables in which forms are affected.

```{r, echo = FALSE, results = TRUE}
# duplicate labels
warning(
"In label_secuTrial.secuTrialdata(d) :
The labels attribute may be longer than 1 for the following variables and forms.
Likely the label was changed from its original state in the secuTrial project setup.
variables: a_variable_name
forms: a_form_name"
)
```

# Merging forms in the `secuTrialdata` object

Naturally, you will sometimes need to merge/join some of the data from the individual form data stored
in your `secuTrialdata` object. To achieve this you can use the base R `merge()` function.
For our dataset we might be interested in merging the `baseline` form data with that of the `treatment`
form. For this we can use the `mnpcvpid` which uniquely identifies each participant visit. Since we
are only interested in the `rando_treatment` variable we will shrink the data in the `treatment` form prior to merging.

\vspace{5pt}
```{r}
treatment_shrink <- ctu05_data$treatment[, c("mnpcvpid", "rando_treatment")]
```

Because we do not want to drop non-matching rows from `baseline` we set `all.x = TRUE`.
As you can see from the `dim()` calls, one column has been added after the merge.
This corresponds to the `rando_treatment` variable.

```{r}
bl_treat <- merge(x = ctu05_data$baseline, y = treatment_shrink,
                  by = "mnpcvpid", all.x = TRUE)
# check dimensions
dim(ctu05_data$baseline)
dim(bl_treat)
```

Another common task may be to merge repetition form data to its parent form.
In our case `esurgeries` can be naively merged with `baseline` via the
`mnpdocid` (from the secuTrial manual: "Each eCRF document record has a unique document identifier."):

\vspace{5pt}
```{r}
bl_surg <- merge(x = ctu05_data$baseline, y = ctu05_data$esurgeries, by = "mnpdocid")
```

Please note, that such naive merging can cause duplication of data if the ids that the merge is directed by are not unique.
This also happened during the production of `bl_surg` in the code above. Participant "RPACK-INS-012" exhibits the
`mnpdocid` 234 twice in the `esurgeries` repetition form which causes a duplication of the `baseline` data matching
`mnpdocid` 234.

\vspace{5pt}
```{r}
table(ctu05_data$esurgeries$mnpdocid)
```

Lets briefly illustrate the consequences by looking at a `table()` of the `height` variable from the `baseline`
form before and after merging.

\vspace{5pt}
```{r}
# before merge
table(ctu05_data$baseline$height)

# after merge
table(bl_surg$height)
```

A closer look reveals that `180` now appears four times instead of three, which can be attributed to
the duplication. This is not a favourable outcome because it can cause confusion and misinterpretation.
A better approach is to change the structure of your repetition form before merging to make the ids you merge by unique.
For this you need to investigate which data you would like to merge and design an appropriate stategy.
In our example case we are interested in the `surgery_organ`. Of course we also need to drag the `mnpdocid` along
to perform the actual merge.

\vspace{5pt}
```{r}
# write a temporary object
surg <- ctu05_data$esurgeries[, c("mnpdocid", "surgery_organ.factor")]
# only retain non NA rows
surg <- surg[which(! is.na(surg$surgery_organ.factor)), ]
# show it
surg
```

\newpage

In order to prevent duplication we can restructure the data before merging.

\vspace{5pt}
```{r}
library(tidyr) # pivot_wider
# add a count
surg$count <- 1
# show the data
surg

# make it wide
surg_wide <- pivot_wider(surg, names_from = surgery_organ.factor, values_from = count)
# show the wide data
surg_wide

```

Checking the dimensions before and after merging reveals that the structure, especially the line count, remains
the same except for the data added from the `esurgeries` repetition form (i.e. `Stomach` and `Other`). Also,
as expected, the `table()` of the height variable returns the expected result (i.e. 180 is present three not four times).

\vspace{5pt}
```{r}
# merge
bl_surg_no_dup <- merge(x = ctu05_data$baseline, y = surg_wide,
                        by = "mnpdocid", all.x = TRUE)

# compare dimensions
dim(bl_surg_no_dup)
dim(ctu05_data$baseline)

# check the height variable
table(bl_surg_no_dup$height)
```

The above description only provides a very brief and simplified example. Merging
strategies need to be individually tailored and require a good understanding of the data at hand.
The `links_secuTrial()` function may be helpful to understand which variables will allow you to
merge forms.

# A note on `mnp*` variables

There is a plethora of variables in the tables of secuTrial exports whose names start
with `mnp`. These are metadata variables which are e.g. important to logically link the
different tables. Explaining them all is beyond the scope of this vignette.
For detailed explanations, please refer to the secuTrial "Export Formats" user manual.

\newpage

```{r}
sessionInfo()
```

\vspace{185pt}

**Disclaimer**

The descriptions of the secuTrial exports used in this vignette and other `secuTrialR`
documentation correspond to our understanding of them and come with no warranty.
For in depth details please refer to the original secuTrial manuals.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/return_scores.R
\name{return_scores}
\alias{return_scores}
\title{Returns the score (calculated) items from \code{secuTrialdata} objects}
\usage{
return_scores(x)
}
\arguments{
\item{x}{a \code{secuTrialdata} object}
}
\value{
a data.frame (columns: name, itemtype, label) that pinpoints which items are scores/calculated.
}
\description{
secuTrial allows to set up calculated fields (i.e. scores) that depend on other items. It is not
             suggested to use the scores calculated by secuTrial to perform reliable analyses. To this end,
             calling \code{return_scores} will return all items in the secuTrial export which are scores and should
             be manually recalculated before data analysis.
}
\examples{
# export location
expot_loc <- system.file("extdata", "sT_exports", "lnames",
                         "s_export_CSV-xls_CTU05_long_ref_miss_en_utf8.zip",
                         package = "secuTrialR")
# read export
sT_export <- read_secuTrial(expot_loc)

# return scores
return_scores(sT_export)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/return_hidden_items.R
\name{return_hidden_items}
\alias{return_hidden_items}
\title{Returns hidden items (variables) from \code{secuTrialdata} objects}
\usage{
return_hidden_items(x)
}
\arguments{
\item{x}{a \code{secuTrialdata} object}
}
\value{
a data.frame (columns: name, itemtype, label) that pinpoints which items are hidden
}
\description{
Sometimes, during a study, certain fields may be hidden because data should
             no longer be entered into them. If this is the case and the data of these
             fields is part of your export is likely good to know about it.
}
\examples{
# export location
expot_loc <- system.file("extdata", "sT_exports", "lnames",
                         "s_export_CSV-xls_CTU05_long_ref_miss_en_utf8.zip",
                         package = "secuTrialR")
# read export
sT_export <- read_secuTrial(expot_loc)

# return scores
return_hidden_items(sT_export)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff_secuTrial.R
\name{diff_secuTrial}
\alias{diff_secuTrial}
\title{Returns differences in the setup of two \code{secuTrialdata} objects}
\usage{
diff_secuTrial(x, y)
}
\arguments{
\item{x}{a \code{secuTrialdata} object (the older export)}

\item{y}{a \code{secuTrialdata} object (the newer export)}
}
\value{
If there are differences, \code{diff_secuTrial()} will produce a list of vectors.
        The fist vector informs about new forms and the second vector informs about
        new variables.
}
\description{
During ongoing studies it is possible that changes to the DataCapture interface
             are made. Sometimes these changes may call for adjustments in analysis code.
             It is considered good practice to run \code{diff_secuTrial()} on the last export
             and the current export of a project to at least make yourself aware of
             potential changes to the setup. If there are differences, the results of this function should
             be interpreted as a first indicator since they may not cover all alterations.
             Information is returned on new forms and variables.
             A detailed list of changes can be produced in the FormBuilder with
             "Compare project setup".
}
\examples{
# read exports

# v1 is essentially a clone of the CTU05 setup
ctu06_v1 <- read_secuTrial(system.file("extdata", "sT_exports", "change_tracking",
                                       "s_export_CSV-xls_CTU06_version1.zip",
                                       package = "secuTrialR"))
# v2 contains 2 additional forms (mnpctu06anewform, mnpctu06anothernewform) and
# 2 additional variables (new_item_in_fu, new_item_in_new_form)
ctu06_v2 <- read_secuTrial(system.file("extdata", "sT_exports", "change_tracking",
                                       "s_export_CSV-xls_CTU06_version2.zip",
                                       package = "secuTrialR"))
# return diff
diff_secuTrial(ctu06_v1, ctu06_v2)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_recruitment.R
\name{plot_recruitment}
\alias{plot_recruitment}
\title{Plots the recruitment over time for \code{secuTrialdata} objects}
\usage{
plot_recruitment(
  x,
  return_data = FALSE,
  show_centres = TRUE,
  cex = 1,
  rm_regex = ""
)
}
\arguments{
\item{x}{a \code{secuTrialdata} object}

\item{return_data}{logical - return the data used to produce the plot instead of the plot}

\item{show_centres}{logical - subset the data into centres}

\item{cex}{double - specifies font size in legend}

\item{rm_regex}{character - specifies a regular expression to be removed from the centre names in the legend.
e.g. rm_regex = "\\\(.*\\\)$" will remove trailing brackets and their contents.}
}
\value{
a simple line plot showing recruitment over time
        or a list of data.frames if return_data is set to TRUE
}
\description{
secuTrial exports inherently contain the information on which participant was
             registered at which point in time. This function makes use of this property
             to plot recruitment over time. Centers indicated with a black line in the
             legend did not recruit at all.
}
\examples{
# export location
expot_loc <- system.file("extdata", "sT_exports", "lnames",
                         "s_export_CSV-xls_CTU05_long_ref_miss_en_utf8.zip",
                         package = "secuTrialR")
# read export
sT_export <- read_secuTrial(expot_loc)

# plot recruitment
plot_recruitment(sT_export)

# plot without trailing bracket
plot_recruitment(sT_export, rm_regex = "\\\\(.*\\\\)$")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dates_secuTrial.R
\name{dates_secuTrial}
\alias{dates_secuTrial}
\alias{dates_secuTrial.secuTrialdata}
\title{Methods to handle date(times)s in secuTrial exports}
\usage{
dates_secuTrial(object, ...)

\method{dates_secuTrial}{secuTrialdata}(object, ...)
}
\arguments{
\item{object}{\code{secuTrialdata} object}

\item{...}{further parameters}
}
\value{
same as the original object with date variables converted to \code{Date}s.
}
\description{
Converts dates and datetime variables to \code{Date} or \code{POSIXct} class, as appropriate.
}
\details{
New variables are created appended with \code{.date} or \code{.datetime}.
         This is a safety mechanism in case NAs are inadvertently introduced.
}
\examples{
# prepare path to example export
export_location <- system.file("extdata", "sT_exports", "lnames",
                               "s_export_CSV-xls_CTU05_long_ref_miss_en_utf8.zip",
                               package = "secuTrialR")
# load all export data
sT_export <- read_secuTrial_raw(data_dir = export_location)
# prepare dates
sT_export_dates <- dates_secuTrial(sT_export)

# show parsed datetime example
sT_export_dates$ctu05baseline$hiv_date.datetime[1]
# [1] "2019-03-05 23:56:00 CET"
# show parsed date example
sT_export_dates$ctu05baseline$paracetamol_start.date[1]
# [1] "2019-03-05"
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_secuTrial.R
\name{write_secuTrial}
\alias{write_secuTrial}
\alias{write_secuTrial.secuTrialdata}
\alias{write_secuTrial.data.frame}
\title{Write secuTrial exports to other formats}
\usage{
write_secuTrial(object, ...)

\method{write_secuTrial}{secuTrialdata}(object, format = "dta", metadata = FALSE, ...)

\method{write_secuTrial}{data.frame}(df, filename, path = "", format = "dta", ...)
}
\arguments{
\item{object}{\code{secuTrialdata} object}

\item{...}{further parameters}

\item{format}{format in which to save the export (one of "dta", "sas", "sav", "xpt")}

\item{metadata}{if TRUE then metadate files will also be written}

\item{df}{a data.frame}

\item{filename}{file name}

\item{path}{directory where the files should be saved}
}
\value{
a list of filenames
}
\description{
Convert the export prepared in R and export it to
             SPSS (sav), Stata (dta) or SAS (sas, xpt version 8)
             using the haven package.
}
\details{
Due to variable naming limitations in other packages, date variables are
         appended with _d (rather than _date), datetime/POSIX variables are appended
         with _dt (rather than _datetime) and factors with _f (rather than _factor).
         Further variable names may be altered in the conversion process.
         For details please refer to the \code{haven} documentation.
}
\examples{
# prepare path to example export
export_location <- system.file("extdata", "sT_exports", "lnames",
                               "s_export_CSV-xls_CTU05_long_ref_miss_en_utf8.zip",
                               package = "secuTrialR")
# load all export data
sT_export <- read_secuTrial(data_dir = export_location)
tdir <- tempdir()
write_secuTrial(sT_export, format = "dta", path = tdir)
list.files(tdir)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/labels_secuTrial.R
\name{labels_secuTrial}
\alias{labels_secuTrial}
\alias{label_secuTrial}
\alias{label}
\alias{units}
\alias{label<-}
\alias{units<-}
\title{Get variable labels for secuTrialdata objects}
\usage{
labels_secuTrial(object, form = NULL)

label_secuTrial(object, ...)

label(x)

units(x)

label(x) <- value

units(x) <- value
}
\arguments{
\item{object}{a \code{secuTrialdata} object}

\item{form}{which form (string)}

\item{...}{further parameters}

\item{x}{any object}

\item{value}{any object}
}
\value{
\code{labels_secuTrial} returns a named vector
\code{label_secuTrial} returns the same object as \code{object}, but with labels added to variables
and data.frames
\code{label} and \code{units} return strings with the appropriate labels

\code{secuTrialdata} object with labels applied to each variable
}
\description{
Variable labels are important for understanding the contents of a variable.
             \code{secuTrialR} offers two main methods to get those labels. \code{labels_secuTrial}
             returns a named list of labels. \code{label_secuTrial} adds labels and units to
             variables (and data.frames) which can then be queried via \code{label} or \code{units}.
}
\details{
For \code{labels_secuTrial}, regular expressions are used with \code{form}
         (specifically, it is inserted between \code{(} and \code{)$} to identify the form).
         Consequently, if \code{form} matches multiple forms (because the beginning is different),
         multiple forms may be returned. You could be more specific with the regular expression,
         remembering that it is inserted between \code{(} and \code{)$}.
}
\note{
The \code{label_secuTrial}/\code{label} syntax is similar to that used in Hmisc, with the
      advantage that it does not change data types (Hmisc coerces everything to labelled integer).
      Similar to Hmisc, however, most operations will remove the labels.
}
\examples{
# APPROACH 1: labels_secuTrial
# ex. 1
# prepare path to example export
export_location <- system.file("extdata", "sT_exports", "BMD",
                               "s_export_CSV-xls_BMD_short_en_utf8.zip",
                               package = "secuTrialR")
# load all export data
sT_export <- read_secuTrial_raw(data_dir = export_location)
# get all labels
labels <- labels_secuTrial(sT_export)
labels[["age"]]

# ex. 2
# load export
sT_export <- read_secuTrial_raw(system.file("extdata", "sT_exports", "lnames",
                                            "s_export_CSV-xls_CTU05_long_miss_en_utf8.zip",
                                            package = "secuTrialR"))

# get labels for sae, treatment and surgeries forms
labels <- labels_secuTrial(sT_export, form = c("sae", "treatment", "surgeries"))


# APPROACH 2: label_secuTrial
# load secuTrial export with separate reference table
sT_export <- read_secuTrial_raw(system.file("extdata", "sT_exports", "lnames",
                                            "s_export_CSV-xls_CTU05_long_ref_miss_en_utf8.zip",
                                            package = "secuTrialR"))
# label the secuTrialdata object
sT_export_labelled <- label_secuTrial(sT_export)
# form label
label(sT_export_labelled$ctu05baseline)
# variable label
label(sT_export_labelled$ctu05baseline$visit_date)
# sampling units
units(sT_export_labelled$ctu05baseline$height)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_secuTrial.R
\name{read_secuTrial}
\alias{read_secuTrial}
\title{Read secuTrial export}
\usage{
read_secuTrial(data_dir, labels = TRUE, factor = TRUE, dates = TRUE)
}
\arguments{
\item{data_dir}{string - location of the export}

\item{labels}{logical - add labels to variables and table}

\item{factor}{logical - convert categorical variables to factor variables
(ignored when reference values are not in a separate table)}

\item{dates}{logical - convert date variables}
}
\value{
\code{secuTrialdata} object - a list with one data.frame for each file on the export
        and a list containing the export options
}
\description{
Convenience wrapper for \code{read_secuTrial_raw}, \code{label_secuTrial},
             \code{factorize_secuTrial} and \code{dates_secuTrial}.
}
\examples{
export_location <- system.file("extdata", "sT_exports", "lnames",
                               "s_export_CSV-xls_CTU05_long_ref_miss_en_utf8.zip",
                               package = "secuTrialR")
# read all export data
sT_export <- read_secuTrial(data_dir = export_location)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annual_recruitment.R
\name{annual_recruitment}
\alias{annual_recruitment}
\title{Returns a data.frame showing the annual recruitment per center}
\usage{
annual_recruitment(x, rm_regex = "")
}
\arguments{
\item{x}{a \code{secuTrialdata} object}

\item{rm_regex}{character - specifies a regular expression to be removed from the centre names in the legend.
e.g. rm_regex = "\\\(.*\\\)$" will remove trailing brackets and their contents.}
}
\value{
a data.frame showing the annual recruitment counts per center
}
\description{
secuTrial exports inherently contain the information on which participant was
             registered at which point in time. This function makes use of this property
             to show annual recruitment.
}
\note{
This function wraps plot_recruitment to retrieve the data.
}
\examples{
# export location
expot_loc <- system.file("extdata", "sT_exports", "lnames",
                         "s_export_CSV-xls_CTU05_long_ref_miss_en_utf8.zip",
                         package = "secuTrialR")
# read export
sT_export <- read_secuTrial(expot_loc)

# plot recruitment
annual_recruitment(sT_export)

# show without trailing bracket
annual_recruitment(sT_export, rm_regex = "\\\\(.*\\\\)$")

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.R
\name{get_participants}
\alias{get_participants}
\title{Retrieves participants present in \code{secuTrialdata}}
\usage{
get_participants(dat)
}
\arguments{
\item{dat}{\code{secuTrialdata} object containing participant IDs and centre information}
}
\value{
data.frame containing participants present in dat
}
\description{
Given a \code{secuTrialdata} object, this function simply returns a list of participants.
Information included are participant IDs and corresponding study centre information,
if available.
}
\examples{

path <- system.file("extdata", "sT_exports", "exp_opt",
                    "s_export_CSV-xls_CTU05_all_info.zip",
                    package = "secuTrialR")
sT_export <- read_secuTrial(path)

# show participants
participants <- get_participants(sT_export)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_secuTrial_raw.R
\name{secuTrialdata}
\alias{secuTrialdata}
\alias{read_secuTrial_raw}
\alias{print.secuTrialdata}
\title{This function loads a secuTrial export}
\usage{
read_secuTrial_raw(data_dir)

\method{print}{secuTrialdata}(x, ...)
}
\arguments{
\item{data_dir}{string The data_dir specifies the path to the secuTrial data export.}

\item{x}{secuTrialdata object as returned by \code{read_secuTrial_raw}}

\item{...}{further parameters}
}
\value{
\code{secuTrialdata} object containing a list of
        export options and data.frames with all the data loaded from
        the secuTrial export.
        The list will contain at least the metadata data.frames and
        export_options list.

data.frame with a row for each table in the export. For each table it
        contains the name, number of rows and columns, an indicator for
        whether the table is a metadata table and the files original name.
}
\description{
This function will always load the full set of meta and data tables.
The export options are also loaded and written into export_options.
}
\examples{
# prepare path to example export
export_location <- system.file("extdata", "sT_exports", "BMD",
                               "s_export_CSV-xls_BMD_short_en_utf8.zip",
                               package = "secuTrialR")
# read all export data
sT_export <- read_secuTrial_raw(data_dir = export_location)

# Print method
print(sT_export)
# or
sT_export
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_validation_overview.R
\name{read_validation_overview}
\alias{read_validation_overview}
\title{This function loads a multi-page secuTrial 'Validation Overview' report into an R tibble.}
\usage{
read_validation_overview(data_dir, skip = 7)
}
\arguments{
\item{data_dir}{Path to the Validation Overview (must be an *.xlsx file).}

\item{skip}{Equivalent parameter in read_excel().
The validation overview xlsx files contain some information in the first
few lines of each sheet which need to be skipped in order to produce the
correct header in R. Prior to reading the validation overview with read_validation_overview()
it is likely a good idea to check how many lines need to be skipped. This
parameter has been added because the amount of lines can differ between different
versions of secuTrial.}
}
\value{
tibble with the 'Validation Overview' data
}
\description{
This function loads a multi-page secuTrial 'Validation Overview' report into an R tibble.
}
\examples{
val_ovv_location <- system.file("extdata", "sT_exports", "BMD",
                                "bmd_validation_overview.xlsx",
                                package = "secuTrialR")
val_ovv <- read_validation_overview(data_dir = val_ovv_location)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visit_structure.R
\name{visit_structure}
\alias{visit_structure}
\alias{plot.secuTrialvisit}
\title{Get the visit structure of \code{secuTrialdata} objects}
\usage{
visit_structure(x, sorted = TRUE)

\method{plot}{secuTrialvisit}(x, ...)
}
\arguments{
\item{x}{a \code{secuTrialdata} object}

\item{sorted}{logical if TRUE sorted by first visit}

\item{...}{further parameters}
}
\value{
data.frame with 1 for whether a form (rows) was collected during a particular visit (columns)

plot of the visit plan
}
\description{
Get the visit structure of \code{secuTrialdata} objects
}
\note{
Requires a fixed visit structure - an error will be returned for projects without
      a visit structure or one with flexible visits
}
\examples{
export_location <- system.file("extdata", "sT_exports", "lnames",
                               "s_export_CSV-xls_CTU05_long_ref_miss_en_utf8.zip",
                               package = "secuTrialR")
# read all export data
sT_export <- read_secuTrial(data_dir = export_location)
# get visit structure
vs <- visit_structure(sT_export)
# plot
plot(vs)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asdataframe.R
\name{as.data.frame.secuTrialdata}
\alias{as.data.frame.secuTrialdata}
\title{as.data.frame method for secuTrialdata objects
Make the data from the exports more easily accessible by placing them in
another environment (e.g. place them in the global environment
(\code{.GlobalEnv}) and you can reference them without referring to the
\code{secuTrialdata} object anymore. Ie. they become regular \code{data.frame}s).}
\usage{
\method{as.data.frame}{secuTrialdata}(
  x,
  ...,
  envir,
  data.frames = NULL,
  meta = FALSE,
  regex = NULL,
  rep = ""
)
}
\arguments{
\item{x}{\code{secuTrialdata} object}

\item{...}{further parameters}

\item{envir}{environment in which to put the data (e.g. \code{.GlobalEnv})}

\item{data.frames}{character vector of data.frame names to turn into data.frames}

\item{meta}{logical should metadata be returned}

\item{regex}{regex syntax to remove from names}

\item{rep}{replacement for regex}
}
\value{
each \code{data.frame} in the \code{secuTrialdata} object is saved to it's
own \code{data.frame} in the designated environment
}
\description{
as.data.frame method for secuTrialdata objects
Make the data from the exports more easily accessible by placing them in
another environment (e.g. place them in the global environment
(\code{.GlobalEnv}) and you can reference them without referring to the
\code{secuTrialdata} object anymore. Ie. they become regular \code{data.frame}s).
}
\details{
\code{envir} must be specifically defined. For simplicity,
\code{.GlobalEnv} would probably be the easiest (assigning it to another
environment would still entail referring to that environment).
}
\examples{
# prepare path to example export
export_location <- system.file("extdata", "sT_exports", "lnames",
                               "s_export_CSV-xls_CTU05_long_ref_miss_en_utf8.zip",
                               package = "secuTrialR")
# load all export data
sT_export <- read_secuTrial_raw(data_dir = export_location)
# add files to a new environment called env1
env1 <- new.env()
as.data.frame(sT_export, envir = env1)
# add files to a new environment called env2, removing the project name from
# the file names
env2 <- new.env()
as.data.frame(sT_export, regex = "ctu05", envir = env2)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/links_secuTrial.R
\name{links_secuTrial}
\alias{links_secuTrial}
\title{Show links between forms}
\usage{
links_secuTrial(
  object,
  forms = NULL,
  formcol = "#d8b365",
  varcol = "#e5f5f9",
  plot = TRUE
)
}
\arguments{
\item{object}{\code{secuTrialdata} object}

\item{forms}{a regular expression for which forms should be included}

\item{formcol}{color for form name circles}

\item{varcol}{color for variable name circles}

\item{plot}{boolean specifies if the plot should be shown}
}
\value{
a tcltk plot window.
}
\description{
secuTrial creates a large number of files and identifiers with which to link files together.
Understanding the links can be difficult. This function produces a map linking the forms
with common variables.
}
\details{
We recommend to resize the tcltk window and and click view/"fit to screen" to improve readability.
         Forms are colored dull orange, variables are colored light blue.
}
\note{
Note that where a form name is also a variable name, it is appended by \code{_form}
      (igraph requires uniquely named nodes).
}
\examples{
\donttest{
# ex. 1
# prepare path to example export
export_location <- system.file("extdata", "sT_exports", "BMD",
                               "s_export_CSV-xls_BMD_short_en_utf8.zip",
                               package = "secuTrialR")
# load all export data
sT_export <- read_secuTrial_raw(data_dir = export_location)
# plot links
links_secuTrial(sT_export)

# ex. 2
# prepare path to example export
export_location <- system.file("extdata", "sT_exports", "lnames",
                               "s_export_CSV-xls_CTU05_long_ref_miss_en_utf8.zip",
                               package = "secuTrialR")
# load all export data
sT_export <- read_secuTrial_raw(data_dir = export_location)
# plot links for form names starting with "ctu05"
links_secuTrial(sT_export, forms = "^ctu05")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/return_random_participants.R
\name{return_random_participants}
\alias{return_random_participants}
\title{Returns the random participants from a secuTrial export}
\usage{
return_random_participants(
  x,
  centres = "all",
  percent = 0.1,
  date = "1900-01-01",
  seed = 1
)
}
\arguments{
\item{x}{a \code{secuTrialdata} object}

\item{centres}{A character vector of centres for which participants should be returned. If left
unspecified it will return participants for every study centre.}

\item{percent}{A number greater than 0 and smaller than 1 specifying the approximate percentage
of participants to be returned per centre.}

\item{date}{If only participants after a specific date should be considered this can be entered here.
Format should be "YYYY-MM-DD" (e.g. "2011-03-26" for March 26th 2011).
This date is checked against mnpvisstartdate in the casenodes table.}

\item{seed}{Allows to configure a seed for id sampling. Every centre will use a small variation of this seed.}
}
\value{
list of two elements. First, a data.frame that contains the
         random participants from each specified centre. This is performed based on a specified seed to retain
         reproducibilty. Second, the configuration of the randomization (i.e. result of \code{RNGkind()}).
         If the percentage does not yield an integer for a centre the number is tranformed into an integer
         under application of the \code{ceiling()} function (i.e. it is rounded up).
}
\description{
There are situations (e.g. randomized monitoring) in which you may want to
             return a list of random participants from a secuTrial export.
}
\examples{
# export location
expot_loc <- system.file("extdata", "sT_exports", "lnames",
                         "s_export_CSV-xls_CTU05_long_ref_miss_en_utf8.zip",
                         package = "secuTrialR")
# read export
sT_export <- read_secuTrial(expot_loc)

# return random participants
return_random_participants(sT_export, percent = 0.25, seed = 1337, date = "2019-03-18",
                           centres = c("Inselspital Bern (RPACK)", "Charité Berlin (RPACK)"))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_secuTrial_url.R
\name{build_secuTrial_url}
\alias{build_secuTrial_url}
\title{Compose a secuTrial URL}
\usage{
build_secuTrial_url(
  server,
  instance = NA,
  customer = NA,
  projid = NA,
  docid = NA
)
}
\arguments{
\item{server}{string containing a server URL}

\item{instance}{(optional) string containing secuTrial instance name}

\item{customer}{(optional) string containing secuTrial customer label}

\item{projid}{(optional) string containing secuTrial project identifier}

\item{docid}{(optional) secuTrial document/form identifer}
}
\value{
string containing a URL to desired secuTrial page. Currently we provide no
        guarantee that the returned URL is valid.
}
\description{
Given a secuTrial server URL, and optionally instance, customer, project id and document id,
this function composes a URL to a specific secuTrial instance, customer or form.
}
\details{
To find the server and instance of a secuTrial database, simply extract the information from the URL that
you usually use to log in. For example in:

\emph{https://server.secutrial.com/apps/WebObjects/ST21-setup-DataCapture.woa/wa/choose?customer=TES}

\itemize{
\item server id is: \strong{server.secutrial.com}
\item instance id is: \strong{ST21-setup-DataCapture}
\item you can find the customer id at the end of the link i.e \strong{TES}

      Alternatively it can be found in the
      CustomerAdminTool -> click "Kunden" -> Value in DB column of table on login page.
\item you can find the project id in the AdminTool -> Projects -> left side ID column in the table.
\item you can find docids in secuTrial exports as values in the "mnpdocid" column.

      Alternatively they are noted in the form footers of the DataCapture.
}

Also note that only the server address has to be provided, the other arguments are optional.
Thus, there are different scenarios:
\itemize{
\item if only server address is provided, the output will point to the secuTrial server
\item if secuTrial server and instance are provided, the output will point to the secuTrial instance.
\item if secutrial server, instance and customer id are provided, the output will point to the customer page.
\item if secuTrial server, instance, customer, project and document id are provided,
      the output will point to a specific secuTrial form.
}
}
\examples{

# This example, builds pseudo-urls that do not point to an active secuTrial instance.

server <- "server.secutrial.com"
instance <- "ST21-setup-DataCapture"
customer <- "TES"
project <- "7036"
docid <- "181"

build_secuTrial_url(server)
build_secuTrial_url(server, instance)
build_secuTrial_url(server, instance, customer)
build_secuTrial_url(server, instance, customer, project)
build_secuTrial_url(server, instance, customer, project, docid)

# examples of docids (mnpdocid)
path <- system.file("extdata", "sT_exports", "lnames",
                    "s_export_CSV-xls_CTU05_long_ref_miss_en_utf8.zip",
                    package = "secuTrialR")
sT_export <- read_secuTrial(path)

# return docids
docids <- sT_export$ctu05baseline$mnpdocid

# make several links with all docids
build_secuTrial_url(server, instance, customer, project, docids)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/completeness.R
\name{form_status_counts}
\alias{form_status_counts}
\title{A function to assess the status of forms}
\usage{
form_status_counts(object, ...)
}
\arguments{
\item{object}{\code{secuTrialdata} object}

\item{...}{further parameters}
}
\value{
data.frame informing on the status of saved forms per participant
}
\description{
This function returns a data.frame informing on the status of saved
forms per participant. There is a line for every combination of
form type and participant id. The numbers are occurrence counts.
}
\examples{
# prepare path to example export
export_location <- system.file("extdata", "sT_exports", "snames",
                               "s_export_CSV-xls_CTU05_short_ref_miss_en_utf8.zip",
                               package = "secuTrialR")
# load all export data
sT_export <- read_secuTrial(data_dir = export_location)

# get form status
form_status_counts(sT_export)

}
\keyword{completeness}
\keyword{form}
\keyword{status}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/factorize.R
\name{factorize}
\alias{factorize}
\alias{factorize_secuTrial}
\alias{factorize_secuTrial.secuTrialdata}
\title{Add factors to \code{secuTrialdata} objects}
\usage{
factorize_secuTrial(object, ...)

\method{factorize_secuTrial}{secuTrialdata}(object, ...)
}
\arguments{
\item{object}{\code{secuTrialdata} object with additional factor variables in study forms containing categorical data}

\item{...}{further parameters}
}
\value{
factorized \code{secuTrialdata} object
}
\description{
secuTrial can return a codebook of codes and labels for categorical variables, including lookup
             type variables, if this option is selected in the export tool ('reference values as separate table').
             This allows factors to be easily created. Factorize methods exist for \code{secuTrialdata} objects,
             \code{data.frames}, \code{integer}s and \code{logical}s, but the intent is that only the former be
             used by users. The other methods could be used with customized codebooks.
}
\details{
factorize_secuTrial will return an error if the appropriate codebook is not available.
}
\examples{
# load secuTrial export with separate reference table
sT_export <- read_secuTrial_raw(system.file("extdata", "sT_exports", "lnames",
                                            "s_export_CSV-xls_CTU05_long_ref_miss_en_utf8.zip",
                                            package = "secuTrialR"))
# factorize the secuTrialdata object
sT_export_factorized <- factorize_secuTrial(sT_export)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_export_options.R
\name{check_export_options}
\alias{check_export_options}
\title{Returns deviations from suggested export options}
\usage{
check_export_options(dat)
}
\arguments{
\item{dat}{\code{secuTrialdata} object}
}
\description{
Given a \code{secuTrialdata} object, this function returns information on deviations
from suggested export options.
}
\details{
While the package strives to allow loading of as many types of secuTrial data exports
         as possible, there are certain export options which are less likely to cause issues.
         If possible it is suggested to export data which adheres to a suggested option set.
         This function points out deviations from the suggested set of options which are: \cr
         is_zip == TRUE \cr
         refvals_separate == TRUE \cr
         add_id == TRUE \cr
         duplicate_meta == FALSE \cr
         encoding == "UTF-8" \cr
         form_status == TRUE \cr
         centre_info == TRUE \cr
         proj_setup == TRUE \cr
         dict_items$lang == "en" \cr
         hidden_fields == FALSE \cr
         structure == TRUE
}
\examples{
path <- system.file("extdata", "sT_exports", "exp_opt",
                    "s_export_CSV-xls_CTU05_only_column_names.zip",
                    package = "secuTrialR")
sT_export <- read_secuTrial_raw(path)

secuTrialR:::check_export_options(sT_export)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assess_form_variable_completeness.R
\name{assess_form_variable_completeness}
\alias{assess_form_variable_completeness}
\title{Asses completeness of data for each variable in a secuTrial export}
\usage{
assess_form_variable_completeness(
  form,
  casenodes_table,
  validation_overview,
  completeness = "allforms",
  occ_in_vp = 1,
  omit_mnp = TRUE
)
}
\arguments{
\item{form}{data.frame Form for which to assess variable completeness
(i.e. a list element returned by read_secuTrial_raw).}

\item{casenodes_table}{data.frame The central casenodes record file
(i.e. 'casenodes' from the list returned by read_secuTrial_raw).}

\item{validation_overview}{tibble returned by read_validation_overview.}

\item{completeness}{string Specifies if completeness is assessed for all forms ('allforms')
or only for saved forms ('savedforms').}

\item{occ_in_vp}{integer Specifies how often the form occurs in the visit plan
(only relevant if completeness = 'allforms' is set).}

\item{omit_mnp}{boolean Removes variable names from the result table that start with mnp.}
}
\value{
data.frame showing percent completeness per variable.
}
\description{
NOTE: This is not exported currently since it is not generic enough.
      It can still be used but in depth knowledge of the function is required
      to evaluate if the result is correct.
      As a rule of thumb: rigid CDMA setups with low complexity will likely
                          work as expected.

Variable completeness is defined as percentage of data entered for a variable
when considering all occurrences of the variable throughout the visit plan and
all participants registered in the secuTrial data base. The completeness can
be assessed for 'allforms', which will also take unsaved forms into account, and
for 'savedforms' which will only consider forms that have been actively saved.
In order to asses variable completeness the function requires a loaded secuTrial
validation overview and a loaded secuTrial standard export. Both the validation
overview and the data export should be from the same time (not more than minutes
apart if possible) to ensure data integrity.
Variable completeness is considered based on rules in the secuTrial setup. Thus,
if a variable is not marked as necessary for completeness, this variable will
always be considered as 100 percent complete if only 'savedforms' are assessed.
Please note that variable completeness for repetition forms should only be
assessed for saved forms.
}
\examples{
# prepare path to example export
export_location <- system.file("extdata", "sT_exports", "BMD",
                               "s_export_CSV-xls_BMD_short_en_utf8.zip",
                               package = "secuTrialR")
# read all export data
sT_export <- read_secuTrial_raw(data_dir = export_location)

# read validation overview
val_ovv_location <- system.file("extdata", "sT_exports", "BMD",
                                "bmd_validation_overview.xlsx",
                                package = "secuTrialR")
val_ovv <- read_validation_overview(data_dir = val_ovv_location)

secuTrialR:::assess_form_variable_completeness(form = sT_export$bmd,
                                               casenodes_table = sT_export$cn,
                                               validation_overview = val_ovv,
                                               completeness = "allforms",
                                               occ_in_vp = 5)

}
\seealso{
read_validation_overview, read_secuTrial_raw
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subset.R
\name{subset_secuTrial}
\alias{subset_secuTrial}
\title{Subsets a \code{secuTrialdata} object}
\usage{
subset_secuTrial(dat, participant = NULL, centre = NULL, exclude = FALSE)
}
\arguments{
\item{dat}{\code{secuTrialdata} object containing participant IDs and centre information}

\item{participant}{character vector with a selection of participant IDs (mnpaid) used for subsetting}

\item{centre}{character vector with a selection of centre names (mnpctrname) used for subsetting}

\item{exclude}{boolean which if true excludes participants and centres from dat}
}
\value{
\code{secuTrialdata} object containing only those participants that meet the selection criteria.
}
\description{
Given a \code{secuTrialdata} object, and subsetting parameters,
this function filters the data object to include only the desired study participants.
Subsetting is possible based on participants and based on centres. In order to subset
based on participants, participant IDs (mnpaid) musst be present in the export.
In order to subset based on centres, centre information must be included in the export.
}
\details{
Subsetting based on participants only, centers only, or based on both is possible. The value of parameter exclude
determines whether the output will include participants that meet selection criteria (when exclude = FALSE),
or exclude them (when exclude = TRUE). When selecting based on both participants and centres,
exclude = FALSE will include the intersection of participants meeting the selection criteria.
If exclude = TRUE, a complement of union of participant and centre sets is returned.
}
\examples{

path <- system.file("extdata", "sT_exports", "exp_opt",
                    "s_export_CSV-xls_CTU05_all_info.zip",
                    package = "secuTrialR")
sT <- read_secuTrial(path)
participants <- c("RPACK-INS-011", "RPACK-INS-014", "RPACK-INS-015")
centres <- c("Inselspital Bern (RPACK)", "Universitätsspital Basel (RPACK)")

# show all participants
get_participants(sT)

# subset sT_export
sT_subset1 <- subset_secuTrial(dat = sT, participant = participants)
get_participants(sT_subset1)
sT_subset2 <- subset_secuTrial(dat = sT, participant = participants, exclude = TRUE)
get_participants(sT_subset2)
sT_subset3 <- subset_secuTrial(dat = sT, centre = centres, exclude = TRUE)
get_participants(sT_subset3)
sT_subset4 <- subset_secuTrial(dat = sT, participant = participants,
                               centre = centres, exclude = FALSE)
get_participants(sT_subset4)
sT_subset5 <- subset_secuTrial(dat = sT, participant = participants,
                               centre = centres[2], exclude = FALSE)
get_participants(sT_subset5)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/completeness.R
\name{form_status_summary}
\alias{form_status_summary}
\title{A function to show summary statistics for form statuses}
\usage{
form_status_summary(object, ...)
}
\arguments{
\item{object}{\code{secuTrialdata} object}

\item{...}{further parameters}
}
\value{
data.frame summarizing the statuses for each form
}
\description{
This function warps form_status_counts and returns a
data.frame summarizing the statuses for each form.
Only saved forms are considered for the statistic.
}
\examples{
# prepare path to example export
export_location <- system.file("extdata","sT_exports", "snames",
                               "s_export_CSV-xls_CTU05_short_ref_miss_en_utf8.zip",
                               package = "secuTrialR")
# load all export data
sT_export <- read_secuTrial(data_dir = export_location)

# get form status
form_status_summary(sT_export)

}
\keyword{completeness}
\keyword{form}
\keyword{status}
