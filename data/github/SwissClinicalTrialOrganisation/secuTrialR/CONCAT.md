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
