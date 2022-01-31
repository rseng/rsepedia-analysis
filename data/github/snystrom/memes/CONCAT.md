
<!-- README.md is generated from README.Rmd. Please edit that file -->

# memes

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Codecov test
coverage](https://codecov.io/gh/snystrom/memes/branch/master/graph/badge.svg)](https://codecov.io/gh/snystrom/memes?branch=master)
![R-CMD-check-bioc](https://github.com/snystrom/memes/workflows/R-CMD-check-bioc/badge.svg)
![Bioconductor Build
Status](https://bioconductor.org/shields/build/devel/bioc/memes.svg)
![Bioconductor
Lifetime](https://bioconductor.org/shields/years-in-bioc/memes.svg)
<!-- badges: end -->

An R interface to the [MEME Suite](http://meme-suite.org/) family of
tools, which provides several utilities for performing motif analysis on
DNA, RNA, and protein sequences. memes works by detecting a local
install of the MEME suite, running the commands, then importing the
results directly into R.

## Installation

### Bioconductor

memes is currently available on the Bioconductor `devel` branch:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("memes")
```

### Development Version (Github)

You can install the development version of memes from
[GitHub](https://github.com/snystrom/memes) with:

``` r
if (!requireNamespace("remotes", quietly=TRUE))
  install.packages("remotes")
remotes::install_github("snystrom/memes")

# To temporarily bypass the R version 4.1 requirement, you can pull from the following branch:
remotes::install_github("snystrom/memes", ref = "no-r-4")
```

### Docker Container

``` shell
# Get development version from dockerhub
docker pull snystrom/memes_docker:devel
# the -v flag is used to mount an analysis directory, 
# it can be excluded for demo purposes
docker run -e PASSWORD=<password> -p 8787:8787 -v <path>/<to>/<project>:/mnt/<project> snystrom/memes_docker:devel
```

## Detecting the MEME Suite

memes relies on a local install of the [MEME
Suite](http://meme-suite.org/). For installation instructions for the
MEME suite, see the [MEME Suite Installation
Guide](http://meme-suite.org/doc/install.html?man_type=web).

memes needs to know the location of the `meme/bin/` directory on your
local machine. You can tell memes the location of your MEME suite
install in 4 ways. memes will always prefer the more specific definition
if it is a valid path. Here they are ranked from most- to
least-specific:

1.  Manually passing the install path to the `meme_path` argument of all
    memes functions
2.  Setting the path using `options(meme_bin = "/path/to/meme/bin/")`
    inside your R script
3.  Setting `MEME_BIN=/path/to/meme/bin/` in your `.Renviron` file
4.  memes will try the default MEME install location `~/meme/bin/`

If memes fails to detect your install at the specified location, it will
fall back to the next option.

To verify memes can detect your MEME install, use `check_meme_install()`
which uses the search herirarchy above to find a valid MEME install. It
will report whether any tools are missing, and print the path to MEME
that it sees. This can be useful for troubleshooting issues with your
install.

``` r
library(memes)

# Verify that memes detects your meme install
# (returns all green checks if so)
check_meme_install()
#> checking main install
#> ✓ /opt/meme/bin
#> checking util installs
#> ✓ /opt/meme/bin/dreme
#> ✓ /opt/meme/bin/ame
#> ✓ /opt/meme/bin/fimo
#> ✓ /opt/meme/bin/tomtom
#> ✓ /opt/meme/bin/meme
#> x /opt/meme/bin/streme
```

``` r
# You can manually input a path to meme_path
# If no meme/bin is detected, will return a red X
check_meme_install(meme_path = 'bad/path')
#> checking main install
#> x bad/path
```

## The Core Tools

| Function Name |              Use               | Sequence Input | Motif Input | Output                                                           |
|:-------------:|:------------------------------:|:--------------:|:-----------:|:-----------------------------------------------------------------|
| `runStreme()` | Motif Discovery (short motifs) |      Yes       |     No      | `universalmotif_df`                                              |
| `runDreme()`  | Motif Discovery (short motifs) |      Yes       |     No      | `universalmotif_df`                                              |
|  `runAme()`   |        Motif Enrichment        |      Yes       |     Yes     | data.frame (optional: `sequences` column)                        |
|  `runFimo()`  |         Motif Scanning         |      Yes       |     Yes     | GRanges of motif positions                                       |
| `runTomTom()` |        Motif Comparison        |       No       |     Yes     | `universalmotif_df` w/ `best_match_motif` and `tomtom` columns\* |
|  `runMeme()`  | Motif Discovery (long motifs)  |      Yes       |     No      | `universalmotif_df`                                              |

\* **Note:** if `runTomTom()` is run using a `universalmotif_df` the
results will be joined with the `universalmotif_df` results as extra
columns. This allows easy comparison of *de-novo* discovered motifs with
their matches.

**Sequence Inputs** can be any of:

1.  Path to a .fasta formatted file
2.  `Biostrings::XStringSet` (can be generated from GRanges using
    `get_sequence()` helper function)
3.  A named list of `Biostrings::XStringSet` objects (generated by
    `get_sequence()`)

**Motif Inputs** can be any of:

1.  A path to a .meme formatted file of motifs to scan against
2.  A `universalmotif` object, or list of `universalmotif` objects
3.  A `runDreme()` results object (this allows the results of
    `runDreme()` to pass directly to `runTomTom()`)
4.  A combination of all of the above passed as a `list()`
    (e.g. `list("path/to/database.meme", "dreme_results" = dreme_res)`)

**Output Types**:

`runDreme()`, `runStreme()`, `runMeme()` and `runTomTom()` return
`universalmotif_df` objects which are data.frames with special columns.
The `motif` column contains a `universalmotif` object, with 1 entry per
row. The remaining columns describe the properties of each returned
motif. The following column names are special in that their values are
used when running `update_motifs()` and `to_list()` to alter the
properties of the motifs stored in the `motif` column. Be careful about
changing these values as these changes will propagate to the `motif`
column when calling `update_motifs()` or `to_list()`.

-   name
-   altname
-   family
-   organism
-   strand
-   nsites
-   bkgsites
-   pval
-   qval
-   eval

memes is built around the [universalmotif
package](https://www.bioconductor.org/packages/release/bioc/html/universalmotif.html)
which provides a framework for manipulating motifs in R.
`universalmotif_df` objects can interconvert between data.frame and
`universalmotif` list format using the `to_df()` and `to_list()`
functions, respectively. This allows use of `memes` results with all
other Bioconductor motif packages, as `universalmotif` objects can
convert to any other motif type using `convert_motifs()`.

`runTomTom()` returns a special column: `tomtom` which is a `data.frame`
of all match data for each input motif. This can be expanded out using
`tidyr::unnest(tomtom_results, "tomtom")`, and renested with
`nest_tomtom()`. The `best_match_` prefixed columns returned by
`runTomTom()` indicate values for the motif which was the best match to
the input motif.

## Quick Examples

### Motif Discovery with DREME

``` r
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(GenomicRanges))

# Example transcription factor peaks as GRanges
data("example_peaks", package = "memes")

# Genome object
dm.genome <- BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6
```

The `get_sequence` function takes a `GRanges` or `GRangesList` as input
and returns the sequences as a `BioStrings::XStringSet`, or list of
`XStringSet` objects, respectively. `get_sequence` will name each fasta
entry by the genomic coordinates each sequence is from.

``` r
# Generate sequences from 200bp about the center of my peaks of interest
sequences <- example_peaks %>% 
  resize(200, "center") %>% 
  get_sequence(dm.genome)
```

`runDreme()` accepts XStringSet or a path to a fasta file as input. You
can use other sequences or shuffled input sequences as the control
dataset.

``` r
# runDreme accepts all arguments that the commandline version of dreme accepts
# here I set e = 50 to detect motifs in the limited example peak list
# In a real analysis, e should typically be < 1
dreme_results <- runDreme(sequences, control = "shuffle", e = 50)
```

memes is built around the
[universalmotif](https://www.bioconductor.org/packages/release/bioc/html/universalmotif.html)
package. The results are returned in `universalmotif_df` format, which
is an R data.frame that can seamlessly interconvert between data.frame
and `universalmotif` format using `to_list()` to convert to
`universalmotif` list format, and `to_df()` to convert back to
data.frame format. Using `to_list()` allows using `memes` results with
all `universalmotif` functions like so:

``` r
library(universalmotif)

dreme_results %>% 
  to_list() %>% 
  view_motifs()
```

<img src="man/figures/README-unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

### Matching motifs using TOMTOM

Discovered motifs can be matched to known TF motifs using `runTomTom()`,
which can accept as input a path to a .meme formatted file, a
`universalmotif` list, or the results of `runDreme()`.

TomTom uses a database of known motifs which can be passed to the
`database` parameter as a path to a .meme format file, or a
`universalmotif` object.

Optionally, you can set the environment variable `MEME_DB` in
`.Renviron` to a file on disk, or the `meme_db` value in `options` to a
valid .meme format file and memes will use that file as the database.
memes will always prefer user input to the function call over a global
variable setting.

``` r
options(meme_db = system.file("extdata/flyFactorSurvey_cleaned.meme", package = "memes"))
m <- create_motif("CMATTACN", altname = "testMotif")
tomtom_results <- runTomTom(m)
```

``` r
tomtom_results
#>         motif  name   altname consensus alphabet strand icscore type
#> 1 <mot:motif> motif testMotif  CMATTACN      DNA     +-      13  PPM
#>                      bkg best_match_name best_match_altname
#> 1 0.25, 0.25, 0.25, 0.25      prd_FlyReg                prd
#>              best_db_name best_match_offset best_match_pval best_match_eval
#> 1 flyFactorSurvey_cleaned                 0        9.36e-05           0.052
#>   best_match_qval best_match_strand
#> 1          0.0353                 +
#>                                                       best_match_motif
#> 1 <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  tomtom
#> 1 prd_FlyReg, tup_SOLEXA_10, CG13424_Cell, CG11085_Cell, BH2_Cell, CG13424_SOLEXA_2, Tup_Cell, Tup_SOLEXA, Bsh_Cell, Exex_SOLEXA, Odsh_SOLEXA, Unc4_Cell, Ubx_FlyReg, Unc4_SOLEXA, E5_Cell, inv_SOLEXA_5, BH2_SOLEXA, Zen_SOLEXA, CG33980_SOLEXA_2_10, BH1_SOLEXA, CG33980_SOLEXA_2_0, Hgtx_Cell, NK7.1_Cell, Slou_Cell, CG13424_SOLEXA, Zen2_Cell, AbdA_SOLEXA, Antp_SOLEXA, Btn_Cell, Dfd_SOLEXA, Eve_SOLEXA, Ftz_Cell, Hmx_SOLEXA, Hmx_Cell, CG34031_Cell, zen2_SOLEXA_2, En_Cell, Pb_SOLEXA, Slou_SOLEXA, Unpg_Cell, inv_SOLEXA_2, ovo_FlyReg, lim_SOLEXA_2, C15_SOLEXA, Ems_Cell, Btn_SOLEXA, Unpg_SOLEXA, Pb_Cell, Bsh_SOLEXA, Scr_SOLEXA, Zen2_SOLEXA, CG34031_SOLEXA, Eve_Cell, Pph13_Cell, BH1_Cell, CG11085_SOLEXA, CG32532_Cell, en_FlyReg, Dll_SOLEXA, Dfd_Cell, Dr_SOLEXA, Ap_Cell, Ro_Cell, CG4136_SOLEXA, CG33980_SOLEXA, Hbn_SOLEXA, Lbl_Cell, Otp_Cell, Rx_Cell, CG32532_SOLEXA, NK7.1_SOLEXA, Dr_Cell, Odsh_Cell, Al_SOLEXA, Antp_Cell, Hgtx_SOLEXA, Ftz_SOLEXA, Lab_SOLEXA, Dfd_FlyReg, Ap_SOLEXA, Awh_SOLEXA, CG11294_SOLEXA, CG4136_Cell, E5_SOLEXA, Ro_SOLEXA, PhdP_SOLEXA, CG12361_SOLEXA_2, Ind_Cell, Scr_Cell, CG9876_Cell, CG18599_Cell, CG9876_SOLEXA, Otp_SOLEXA, Lbl_SOLEXA, Ubx_Cell, Ubx_SOLEXA, en_SOLEXA_2, Pph13_SOLEXA, Rx_SOLEXA, CG15696_SOLEXA, CG18599_SOLEXA, Ems_SOLEXA, Repo_Cell, Dll_Cell, C15_Cell, CG12361_SOLEXA, Abd-A_FlyReg, Repo_SOLEXA, Zen_Cell, Inv_Cell, En_SOLEXA, Lim3_Cell, Lim1_SOLEXA, CG15696_Cell, Crc_CG6272_SANGER_5, Lab_Cell, CG32105_SOLEXA, Bap_SOLEXA, CG9437_SANGER_5, AbdA_Cell, pho_FlyReg, CG33980_Cell, Cad_SOLEXA, CG4328_SOLEXA, CG4328_Cell, Gsc_Cell, vri_SANGER_5, AbdB_SOLEXA, Xrp1_CG6272_SANGER_5, Al_Cell, Exex_Cell, br-Z4_FlyReg, CG11294_Cell, Aef1_FlyReg, CG7745_SANGER_5, PhdP_Cell, Awh_Cell, prd, tup, lms, CG11085, B-H2, lms, tup, tup, bsh, exex, OdsH, unc-4, Ubx, unc-4, E5, inv, B-H2, zen, CG33980, B-H1, CG33980, HGTX, NK7.1, slou, lms, zen2, abd-A, Antp, btn, Dfd, eve, ftz, Hmx, Hmx, CG34031, zen2, en, pb, slou, unpg, inv, ovo, Lim1, C15, ems, btn, unpg, pb, bsh, Scr, zen2, CG34031, eve, Pph13, B-H1, CG11085, CG32532, en, Dll, Dfd, Dr, ap, ro, CG4136, CG33980, hbn, lbl, otp, Rx, CG32532, NK7.1, Dr, OdsH, al, Antp, HGTX, ftz, lab, Dfd, ap, Awh, CG11294, CG4136, E5, ro, PHDP, CG12361, ind, Scr, CG9876, CG18599, CG9876, otp, lbl, Ubx, Ubx, en, Pph13, Rx, CG15696, CG18599, ems, repo, Dll, C15, CG12361, abd-A, repo, zen, inv, en, Lim3, Lim1, CG15696, crc, lab, CG32105, bap, CG9437, abd-A, pho, CG33980, cad, CG4328, CG4328, Gsc, vri, Abd-B, Xrp1, al, exex, br, CG11294, Aef1, CG7745, PHDP, Awh, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, <S4 class 'universalmotif' [package "universalmotif"] with 20 slots>, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, flyFactorSurvey_cleaned, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, -3, 0, 4, 1, 0, 0, 0, 0, -2, 0, -2, 1, 0, -2, 1, 0, 2, -1, 0, 9.36e-05, 0.000129, 0.00013, 0.000194, 0.000246, 0.000247, 0.000319, 0.000364, 0.000427, 0.000435, 0.000473, 0.000473, 0.000556, 0.000643, 0.000694, 0.000728, 0.000748, 0.000748, 0.000785, 0.000828, 0.000849, 0.000868, 0.000868, 0.000868, 0.00102, 0.00113, 0.00116, 0.00116, 0.00116, 0.00116, 0.00116, 0.00116, 0.00118, 0.00126, 0.00134, 0.00141, 0.00144, 0.00145, 0.00145, 0.00154, 0.00161, 0.00165, 0.00173, 0.00177, 0.00177, 0.00179, 0.00179, 0.0019, 0.00191, 0.00191, 0.00191, 0.00203, 0.00203, 0.00203, 0.00212, 0.00214, 0.00219, 0.00232, 0.00242, 0.00247, 0.0026, 0.00264, 0.00264, 0.00268, 0.00282, 0.00282, 0.00282, 0.00282, 0.00282, 0.00286, 0.00286, 0.00296, 0.00305, 0.00316, 0.0032, 0.0032, 0.00326, 0.00326, 0.00341, 0.00348, 0.00348, 0.00348, 0.00348, 0.00348, 0.00348, 0.0035, 0.00359, 0.00359, 0.00359, 0.00364, 0.00372, 0.00372, 0.00372, 0.00387, 0.00387, 0.00412, 0.00415, 0.00424, 0.00424, 0.00439, 0.00452, 0.00452, 0.00467, 0.00483, 0.00497, 0.00497, 0.00528, 0.00549, 0.0059, 0.00597, 0.00624, 0.00635, 0.00709, 0.00761, 0.0079, 0.00856, 0.00859, 0.00912, 0.00935, 0.0097, 0.00974, 0.0101, 0.0109, 0.0116, 0.0123, 0.0123, 0.0127, 0.0138, 0.0138, 0.014, 0.014, 0.0147, 0.0148, 0.0155, 0.0165, 0.0166, 0.0177, 0.052, 0.0718, 0.0725, 0.108, 0.137, 0.137, 0.177, 0.202, 0.238, 0.242, 0.263, 0.263, 0.309, 0.357, 0.386, 0.405, 0.416, 0.416, 0.436, 0.46, 0.472, 0.483, 0.483, 0.483, 0.566, 0.629, 0.646, 0.646, 0.646, 0.646, 0.646, 0.646, 0.654, 0.702, 0.745, 0.785, 0.8, 0.808, 0.808, 0.858, 0.897, 0.92, 0.961, 0.985, 0.985, 0.993, 0.993, 1.05, 1.06, 1.06, 1.06, 1.13, 1.13, 1.13, 1.18, 1.19, 1.22, 1.29, 1.34, 1.38, 1.45, 1.47, 1.47, 1.49, 1.57, 1.57, 1.57, 1.57, 1.57, 1.59, 1.59, 1.65, 1.7, 1.76, 1.78, 1.78, 1.81, 1.81, 1.9, 1.94, 1.94, 1.94, 1.94, 1.94, 1.94, 1.95, 2, 2, 2, 2.02, 2.07, 2.07, 2.07, 2.15, 2.15, 2.29, 2.31, 2.36, 2.36, 2.44, 2.52, 2.52, 2.6, 2.68, 2.76, 2.76, 2.94, 3.05, 3.28, 3.32, 3.47, 3.53, 3.94, 4.23, 4.39, 4.76, 4.77, 5.07, 5.2, 5.39, 5.41, 5.61, 6.06, 6.42, 6.81, 6.81, 7.03, 7.66, 7.65, 7.78, 7.78, 8.15, 8.26, 8.59, 9.18, 9.2, 9.85, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0353, 0.0368, 0.0369, 0.0369, 0.0369, 0.0369, 0.0369, 0.0372, 0.0372, 0.0372, 0.0372, 0.0372, 0.0372, 0.0372, 0.0372, 0.0372, 0.0372, 0.0372, 0.0372, 0.0372, 0.0372, 0.0372, 0.0379, 0.0379, 0.0381, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0396, 0.0404, 0.0404, 0.0424, 0.0424, 0.0424, 0.0424, 0.0435, 0.0439, 0.0439, 0.0449, 0.046, 0.0464, 0.0464, 0.0489, 0.0504, 0.0536, 0.0538, 0.0557, 0.0562, 0.0622, 0.0662, 0.0681, 0.0727, 0.0727, 0.0765, 0.0778, 0.0797, 0.0797, 0.0819, 0.0877, 0.0923, 0.0964, 0.0964, 0.0987, 0.106, 0.106, 0.106, 0.106, 0.11, 0.111, 0.114, 0.121, 0.121, 0.128, +, -, -, -, -, -, -, -, -, -, -, -, +, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, +, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, +, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, -, +, -, -, -, -, -, -, -, +, -, -, -, +, -, -, -, -, -, -, -, +, -, +, -, -, -, +, +, +, -, -
#> 
#> [Hidden empty columns: family, organism, nsites, bkgsites, pval, qval,
#>   eval.]
```

### Using runDreme results as TOMTOM input

`runTomTom()` will add its results as columns to a `runDreme()` results
data.frame.

``` r
full_results <- dreme_results %>% 
  runTomTom()
```

### Motif Enrichment using AME

AME is used to test for enrichment of known motifs in target sequences.
`runAme()` will use the `MEME_DB` entry in `.Renviron` or
`options(meme_db = "path/to/database.meme")` as the motif database.
Alternately, it will accept all valid inputs similar to `runTomTom()`.

``` r
# here I set the evalue_report_threshold = 30 to detect motifs in the limited example sequences
# In a real analysis, evalue_report_threshold should be carefully selected
ame_results <- runAme(sequences, control = "shuffle", evalue_report_threshold = 30)
ame_results
#> # A tibble: 2 x 17
#>    rank motif_db motif_id motif_alt_id consensus  pvalue adj.pvalue evalue tests
#>   <int> <chr>    <chr>    <chr>        <chr>       <dbl>      <dbl>  <dbl> <int>
#> 1     1 /usr/lo… Eip93F_… Eip93F       ACWSCCRA… 5.14e-4     0.0339   18.8    67
#> 2     2 /usr/lo… Cf2-PB_… Cf2          CSSHNKDT… 1.57e-3     0.0415   23.1    27
#> # … with 8 more variables: fasta_max <dbl>, pos <int>, neg <int>,
#> #   pwm_min <dbl>, tp <int>, tp_percent <dbl>, fp <int>, fp_percent <dbl>
```

## Visualizing Results

`view_tomtom_hits` allows comparing the input motifs to the top hits
from TomTom. Manual inspection of these matches is important, as
sometimes the top match is not always the correct assignment. Altering
`top_n` allows you to show additional matches in descending order of
their rank.

``` r
full_results %>% 
  view_tomtom_hits(top_n = 1)
#> $m01_AGAGC
```

<img src="man/figures/README-unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

It can be useful to view the results from `runAme()` as a heatmap.
`plot_ame_heatmap()` can create complex visualizations for analysis of
enrichment between different region types (see vignettes for details).
Here is a simple example heatmap.

``` r
ame_results %>% 
  plot_ame_heatmap()
```

<img src="man/figures/README-unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

# Scanning for motif occurances using FIMO

The FIMO tool is used to identify matches to known motifs. `runFimo`
will return these hits as a `GRanges` object containing the genomic
coordinates of the motif match.

``` r
# Query MotifDb for a motif
e93_motif <- MotifDb::query(MotifDb::MotifDb, "Eip93F") %>% 
  universalmotif::convert_motifs()
#> See system.file("LICENSE", package="MotifDb") for use restrictions.

# Scan for the E93 motif within given sequences
fimo_results <- runFimo(sequences, e93_motif, thresh = 1e-3)

# Visualize the sequences matching the E93 motif
plot_sequence_heatmap(fimo_results$matched_sequence)  
```

<img src="man/figures/README-unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

## Importing Data from previous runs

memes also supports importing results generated using the MEME suite
outside of R (for example, running jobs on
[meme-suite.org](meme-suite.org), or running on the commandline). This
enables use of preexisting MEME suite results with downstream memes
functions.

| MEME Tool |    Function Name    | File Type  |
|:---------:|:-------------------:|:----------:|
|  Streme   | `importStremeXML()` | streme.xml |
|   Dreme   | `importDremeXML()`  | dreme.xml  |
|  TomTom   | `importTomTomXML()` | tomtom.xml |
|    AME    |    `importAme()`    | ame.tsv\*  |
|   FIMO    |   `importFimo()`    |  fimo.tsv  |
|   Meme    |   `importMeme()`    |  meme.txt  |

\* `importAME()` can also use the “sequences.tsv” output when AME used
`method = "fisher"`, this is optional.

# FAQs

### How do I use memes/MEME on Windows?

The MEME Suite does not currently support Windows, although it can be
installed under [Cygwin](https://www.cygwin.com/) or the [Windows Linux
Subsytem](https://docs.microsoft.com/en-us/windows/wsl/install-win10)
(WSL). Please note that if MEME is installed on Cygwin or WSL, you must
also run R inside Cygwin or WSL to use memes.

An alternative solution is to use
[Docker](https://www.docker.com/get-started) to run a virtual
environment with the MEME Suite installed. We provide a [memes docker
container](https://github.com/snystrom/memes_docker)  
that ships with the MEME Suite, R studio, and all `memes` dependencies
pre-installed.

# Citation

memes is a wrapper for a select few tools from the MEME Suite, which
were developed by another group. In addition to citing memes, please
cite the MEME Suite tools corresponding to the tools you use.

If you use `runDreme()` in your analysis, please cite:

Timothy L. Bailey, “DREME: Motif discovery in transcription factor
ChIP-seq data”, Bioinformatics, 27(12):1653-1659, 2011. [full
text](https://academic.oup.com/bioinformatics/article/27/12/1653/257754)

If you use `runTomTom()` in your analysis, please cite:

Shobhit Gupta, JA Stamatoyannopolous, Timothy Bailey and William
Stafford Noble, “Quantifying similarity between motifs”, Genome Biology,
8(2):R24, 2007. [full text](http://genomebiology.com/2007/8/2/R24)

If you use `runAme()` in your analysis, please cite:

Robert McLeay and Timothy L. Bailey, “Motif Enrichment Analysis: A
unified framework and method evaluation”, BMC Bioinformatics, 11:165,
2010, <doi:10.1186/1471-2105-11-165>. [full
text](http://www.biomedcentral.com/1471-2105/11/165)

If you use `runFimo()` in your analysis, please cite:

Charles E. Grant, Timothy L. Bailey, and William Stafford Noble, “FIMO:
Scanning for occurrences of a given motif”, Bioinformatics,
27(7):1017-1018, 2011. [full
text](http://bioinformatics.oxfordjournals.org/content/early/2011/02/16/bioinformatics.btr064.full)

## Licensing Restrictions

The MEME Suite is free for non-profit use, but for-profit users should
purchase a license. See the [MEME Suite Copyright
Page](http://meme-suite.org/doc/copyright.html) for details.
# memes 1.2.4
* updated NAMESPACE to fix R CMD CHECK note.
# memes 1.2.3
* fixed a bug in `importTomTomXML` where `tomtom` list column would contain missing data if `tomtom` was run using multiple database sources as input.
# memes 1.2.2
* fixed a bug in `runStreme` causing failures on data import for STREME version >= 5.4.1
# memes 1.2.1
* fixed a bug in `runStreme` causing failures when using BStringSetLists as input
# memes 1.1.3
* fixed a bug in `runTomTom` where setting `norc = TRUE` failed on data import
# memes 1.1.1
* `runFimo` now returns `NULL` and prints a message when `text = FALSE` and FIMO detects no matches instead of throwing a cryptic error message
# memes 0.99.11
* Add support for STREME with `runStreme()`. STREME will supercede DREME in a future MEME Suite release.
# memes 0.99.10
* Fixed a bug where paths weren't correctly expanded when used as `database` entry under certain conditions
# memes 0.99.8
* Removed inline `r` call in integrative_analysis vignette to fix issue on bioc build machine
# memes 0.99.7
* Version bump to force pkg rebuild
# memes 0.99.6
* added list S3 method for `plot_sequence_heatmap` so now named lists of sequences can be passed natively to this function.
  * updated ChIP-seq vignette to demonstrate this

# memes 0.99.5
* added `plot_sequence_heatmap` for making heatmaps of sequence lists
* Added significantly more explanation to the ChIP-seq vignette
* renamed `ame_plot_heatmap` -> `plot_ame_heatmap` for consistency 

# memes 0.1.2
* `runFimo()` `skip_matched_sequence` default is now `FALSE`. Set this to `TRUE` if fimo takes a long time to run, then use `add_sequence()` to add it back if needed.
* `runTomTom()` `dist` default is now `ed` (changed from `pearson`).

# memes 0.1.0
* Removed `as_universalmotif_df()`, `as_universalmotif()`, and `update_motifs()`.
  * These functions are replaced by `universalmotif::to_df()`, `universalmotif::to_list()`, and `universalmotif::update_motifs()`
* `runDreme` and `runTomTom` results are now returned in `universalmotif_df` format (behaves just like a data.frame)
  * The `motif` column of `universalmotif_df` objects can no longer be called directly in `universalmotif` operations like `view_motifs(df$motif)`. Use `to_list()` for this behavior instead.
  * To support this change, the `pvalue`, `evalue`, and `qvalue` columns are renamed `pval`, `eval`, and `qval`. The same is true for tomtom output columns `match_pvalue` -> `match_pval`, `best_match_pvalue` -> `best_match_pval`, etc.
  * Updated example datasets to use `unviversalmotif_df` type
* `ame_plot_heatmap` ranking issue is resolved, plots now sort correctly
* Added `remove_duplicate_motifs` and `has_duplicate_motifs` for detecting and removing duplicated matrices in a universalmotif object or data.frame
* Overhauled the Tidying Motifs vignette for more extensive EDA and a demo of deduplication
  * Updated the `flyFactorSurvey_cleaned.meme` example database to reflect new changes to the vignette

---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  fig.align = "center"
  #out.width = "100%"
)
```

# memes

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Codecov test coverage](https://codecov.io/gh/snystrom/memes/branch/master/graph/badge.svg)](https://codecov.io/gh/snystrom/memes?branch=master)
![R-CMD-check-bioc](https://github.com/snystrom/memes/workflows/R-CMD-check-bioc/badge.svg)
![Bioconductor Build Status](https://bioconductor.org/shields/build/devel/bioc/memes.svg)
![Bioconductor Lifetime](https://bioconductor.org/shields/years-in-bioc/memes.svg)
<!-- badges: end -->

An R interface to the [MEME Suite](http://meme-suite.org/) family of tools,
which provides several utilities for performing motif analysis on DNA, RNA, and
protein sequences. memes works by detecting a local install of the MEME suite,
running the commands, then importing the results directly into R.

## Installation

### Bioconductor

memes is currently available on the Bioconductor `devel` branch:

```{r, eval=F}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("memes")
```

### Development Version (Github)
You can install the development version of memes from [GitHub](https://github.com/snystrom/memes) with:

```{r, eval=F}
if (!requireNamespace("remotes", quietly=TRUE))
  install.packages("remotes")
remotes::install_github("snystrom/memes")

# To temporarily bypass the R version 4.1 requirement, you can pull from the following branch:
remotes::install_github("snystrom/memes", ref = "no-r-4")
```

### Docker Container
```{shell}
# Get development version from dockerhub
docker pull snystrom/memes_docker:devel
# the -v flag is used to mount an analysis directory, 
# it can be excluded for demo purposes
docker run -e PASSWORD=<password> -p 8787:8787 -v <path>/<to>/<project>:/mnt/<project> snystrom/memes_docker:devel
```


## Detecting the MEME Suite

memes relies on a local install of the [MEME Suite](http://meme-suite.org/).
For installation instructions for the MEME suite, see the [MEME Suite
Installation Guide](http://meme-suite.org/doc/install.html?man_type=web).

memes needs to know the location of the `meme/bin/` directory on your local machine.
You can tell memes the location of your MEME suite install in 4 ways. memes
will always prefer the more specific definition if it is a valid path. Here they
are ranked from most- to least-specific:

1. Manually passing the install path to the `meme_path` argument of all memes functions
2. Setting the path using `options(meme_bin = "/path/to/meme/bin/")` inside your R script
3. Setting `MEME_BIN=/path/to/meme/bin/` in your `.Renviron` file
4. memes will try the default MEME install location `~/meme/bin/`

If memes fails to detect your install at the specified location, it will fall
back to the next option.

To verify memes can detect your MEME install, use `check_meme_install()` which
uses the search herirarchy above to find a valid MEME install. It will report
whether any tools are missing, and print the path to MEME that it sees. This can
be useful for troubleshooting issues with your install.
```{r check_install_works}
library(memes)

# Verify that memes detects your meme install
# (returns all green checks if so)
check_meme_install()
```

```{r check_install_fails}
# You can manually input a path to meme_path
# If no meme/bin is detected, will return a red X
check_meme_install(meme_path = 'bad/path')
```

## The Core Tools

| Function Name | Use              | Sequence Input | Motif Input | Output |
|:-------------:|:----------------:|:--------------:|:-----------:|:-------------------------------------------------------|
| `runStreme()` | Motif Discovery (short motifs)  | Yes | No      | `universalmotif_df`                                    |
| `runDreme()`  | Motif Discovery (short motifs)  | Yes | No      | `universalmotif_df`                                    |
| `runAme()`    | Motif Enrichment                | Yes | Yes     | data.frame (optional: `sequences` column)              |
| `runFimo()`   | Motif Scanning                  | Yes | Yes     | GRanges of motif positions                             |
| `runTomTom()` | Motif Comparison                | No  | Yes     | `universalmotif_df` w/ `best_match_motif` and `tomtom` columns* |
| `runMeme()`   | Motif Discovery (long motifs)   | Yes | No      | `universalmotif_df`                                    |

\* **Note:** if `runTomTom()` is run using a `universalmotif_df`
the results will be joined with the `universalmotif_df` results as extra
columns. This allows easy comparison of *de-novo* discovered motifs with their
matches.

**Sequence Inputs** can be any of:

1. Path to a .fasta formatted file
2. `Biostrings::XStringSet` (can be generated from GRanges using `get_sequence()` helper function)
3. A named list of `Biostrings::XStringSet` objects (generated by `get_sequence()`)

**Motif Inputs** can be any of:

1. A path to a .meme formatted file of motifs to scan against
2. A `universalmotif` object, or list of `universalmotif` objects
3. A `runDreme()` results object (this allows the results of `runDreme()` to pass directly to `runTomTom()`)
4. A combination of all of the above passed as a `list()` (e.g. `list("path/to/database.meme", "dreme_results" = dreme_res)`)

**Output Types**:

`runDreme()`, `runStreme()`, `runMeme()` and `runTomTom()` return
`universalmotif_df` objects which are data.frames with special columns. The
`motif` column contains a `universalmotif` object, with 1 entry per row. The
remaining columns describe the properties of each returned motif. The following
column names are special in that their values are used when running
`update_motifs()` and `to_list()` to alter the properties of the motifs stored
in the `motif` column. Be careful about changing these values as these changes
will propagate to the `motif` column when calling `update_motifs()` or
`to_list()`.

 - name
 - altname
 - family
 - organism
 - strand
 - nsites
 - bkgsites
 - pval
 - qval
 - eval

memes is built around the [universalmotif package](https://www.bioconductor.org/packages/release/bioc/html/universalmotif.html)
which provides a framework for manipulating motifs in R. `universalmotif_df`
objects can interconvert between data.frame and `universalmotif` list format
using the `to_df()` and `to_list()` functions, respectively. This allows use of
`memes` results with all other Bioconductor motif packages, as `universalmotif`
objects can convert to any other motif type using `convert_motifs()`.

`runTomTom()` returns a special column: `tomtom` which is a `data.frame` of all
match data for each input motif. This can be expanded out using
`tidyr::unnest(tomtom_results, "tomtom")`, and renested with `nest_tomtom()`.
The `best_match_` prefixed columns returned by `runTomTom()` indicate values for
the motif which was the best match to the input motif.

## Quick Examples
### Motif Discovery with DREME

```{r}
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(GenomicRanges))

# Example transcription factor peaks as GRanges
data("example_peaks", package = "memes")

# Genome object
dm.genome <- BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6
```
The `get_sequence` function takes a `GRanges` or `GRangesList` as input and
returns the sequences as a `BioStrings::XStringSet`, or list of `XStringSet`
objects, respectively. `get_sequence` will name each fasta entry by the genomic
coordinates each sequence is from.
```{r}
# Generate sequences from 200bp about the center of my peaks of interest
sequences <- example_peaks %>% 
  resize(200, "center") %>% 
  get_sequence(dm.genome)
```

`runDreme()` accepts XStringSet or a path to a fasta file as input. You can use
other sequences or shuffled input sequences as the control dataset.
```{r}
# runDreme accepts all arguments that the commandline version of dreme accepts
# here I set e = 50 to detect motifs in the limited example peak list
# In a real analysis, e should typically be < 1
dreme_results <- runDreme(sequences, control = "shuffle", e = 50)
```
memes is built around the
[universalmotif](https://www.bioconductor.org/packages/release/bioc/html/universalmotif.html)
package. The results are returned in `universalmotif_df` format, which is an R data.frame that can seamlessly interconvert between data.frame and `universalmotif` format using `to_list()` to convert to `universalmotif` list format, and `to_df()` to convert back to data.frame format. Using `to_list()` allows using `memes` results with all `universalmotif` functions like so:

```{r}
library(universalmotif)

dreme_results %>% 
  to_list() %>% 
  view_motifs()
```

### Matching motifs using TOMTOM
Discovered motifs can be matched to known TF motifs using `runTomTom()`, which can accept as input a path to a .meme formatted file, a `universalmotif` list, or the results of `runDreme()`.

TomTom uses a database of known motifs which can be passed to the `database` parameter as a path to a .meme format file, or a `universalmotif` object.

Optionally, you can set the environment variable `MEME_DB` in `.Renviron` to a file on disk, or
the `meme_db` value in `options` to a valid .meme format file and memes will
use that file as the database. memes will always prefer user input to the
function call over a global variable setting.

```{r}
options(meme_db = system.file("extdata/flyFactorSurvey_cleaned.meme", package = "memes"))
m <- create_motif("CMATTACN", altname = "testMotif")
tomtom_results <- runTomTom(m)
```
```{r}
tomtom_results
```

### Using runDreme results as TOMTOM input
`runTomTom()` will add its results as columns to a `runDreme()` results data.frame.
```{r}
full_results <- dreme_results %>% 
  runTomTom()
```


### Motif Enrichment using AME

AME is used to test for enrichment of known motifs in target sequences. `runAme()`
will use the `MEME_DB` entry in `.Renviron` or `options(meme_db =
"path/to/database.meme")` as the motif database. Alternately, it will accept all
valid inputs similar to `runTomTom()`.
```{r}
# here I set the evalue_report_threshold = 30 to detect motifs in the limited example sequences
# In a real analysis, evalue_report_threshold should be carefully selected
ame_results <- runAme(sequences, control = "shuffle", evalue_report_threshold = 30)
ame_results
```


## Visualizing Results

`view_tomtom_hits` allows comparing the input motifs to the top hits from
TomTom. Manual inspection of these matches is important, as sometimes the top
match is not always the correct assignment. Altering `top_n` allows you to show
additional matches in descending order of their rank.

```{r}
full_results %>% 
  view_tomtom_hits(top_n = 1)
```

It can be useful to view the results from `runAme()` as a heatmap. 
`plot_ame_heatmap()` can create complex visualizations for analysis of enrichment
between different region types (see vignettes for details). Here is a simple
example heatmap.
```{r, fig.height=3, fig.width=5}
ame_results %>% 
  plot_ame_heatmap()
```

# Scanning for motif occurances using FIMO

The FIMO tool is used to identify matches to known motifs. `runFimo` will return
these hits as a `GRanges` object containing the genomic coordinates of the motif
match.
```{r, fig.height=4, fig.width=3}
# Query MotifDb for a motif
e93_motif <- MotifDb::query(MotifDb::MotifDb, "Eip93F") %>% 
  universalmotif::convert_motifs()

# Scan for the E93 motif within given sequences
fimo_results <- runFimo(sequences, e93_motif, thresh = 1e-3)

# Visualize the sequences matching the E93 motif
plot_sequence_heatmap(fimo_results$matched_sequence)  
```

## Importing Data from previous runs

memes also supports importing results generated using the MEME suite outside of
R (for example, running jobs on [meme-suite.org](meme-suite.org), or running on
the commandline). This enables use of preexisting MEME suite results with
downstream memes functions.

| MEME Tool | Function Name       | File Type        |
|:---------:|:-------------------:|:----------------:|
| Streme    | `importStremeXML()` | streme.xml       |
| Dreme     | `importDremeXML()`  | dreme.xml        |
| TomTom    | `importTomTomXML()` | tomtom.xml       |
| AME       | `importAme()`       | ame.tsv*         |
| FIMO      | `importFimo()`      | fimo.tsv         | 
| Meme      | `importMeme()`      | meme.txt         | 

\* `importAME()` can also use the "sequences.tsv" output when AME used `method = "fisher"`, this is optional.


# FAQs
### How do I use memes/MEME on Windows?
The MEME Suite does not currently support Windows, although it can be
installed under [Cygwin](https://www.cygwin.com/) or the [Windows Linux Subsytem](https://docs.microsoft.com/en-us/windows/wsl/install-win10) (WSL).
Please note that if MEME is installed on Cygwin or WSL, you must also run R
inside Cygwin or WSL to use memes.

An alternative solution is to use [Docker](https://www.docker.com/get-started)
to run a virtual environment with the MEME Suite installed. We provide a [memes docker container](https://github.com/snystrom/memes_docker)  
that ships with the MEME Suite, R studio, and all `memes` dependencies
pre-installed. 

# Citation

memes is a wrapper for a select few tools from the MEME Suite, which were
developed by another group. In addition to citing memes, please cite the MEME
Suite tools corresponding to the tools you use.

If you use `runDreme()` in your analysis, please cite:

Timothy L. Bailey, "DREME: Motif discovery in transcription factor ChIP-seq data", Bioinformatics, 27(12):1653-1659, 2011. [full text](https://academic.oup.com/bioinformatics/article/27/12/1653/257754)

If you use `runTomTom()` in your analysis, please cite:

Shobhit Gupta, JA Stamatoyannopolous, Timothy Bailey and William Stafford Noble, "Quantifying similarity between motifs", Genome Biology, 8(2):R24, 2007. [full text](http://genomebiology.com/2007/8/2/R24)

If you use `runAme()` in your analysis, please cite:

Robert McLeay and Timothy L. Bailey, "Motif Enrichment Analysis: A unified framework and method evaluation", BMC Bioinformatics, 11:165, 2010, doi:10.1186/1471-2105-11-165. [full text](http://www.biomedcentral.com/1471-2105/11/165)

If you use `runFimo()` in your analysis, please cite:

Charles E. Grant, Timothy L. Bailey, and William Stafford Noble, "FIMO: Scanning for occurrences of a given motif", Bioinformatics, 27(7):1017-1018, 2011. [full text](http://bioinformatics.oxfordjournals.org/content/early/2011/02/16/bioinformatics.btr064.full)

## Licensing Restrictions
The MEME Suite is free for non-profit use, but for-profit users should purchase a
license. See the [MEME Suite Copyright Page](http://meme-suite.org/doc/copyright.html) for details.
---
title: "Motif Scanning using FIMO"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Motif Scanning using FIMO}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include=FALSE}
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")

eval_vignette <- NOT_CRAN & memes::meme_is_installed()

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  purl = eval_vignette,
  eval = eval_vignette
)
```

# See package website for full vignette

The Bioconductor build system does not have the MEME Suite installed, therefore
these vignettes will not contain any R output. To view the full vignette, visit
this article page on the memes website [at this link](https://snystrom.github.io/memes-manual/articles/core_fimo.html)

```{r setup}
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(universalmotif))
library(memes)
```

## Inputs
FIMO searches input sequences for occurrances of a motif. `runFimo()` has two
required inputs: fasta-format sequences, with optional genomic coordinate
headers, and a set of motifs to detect within the input sequences.

### Sequence Inputs:
Sequence input to `runFimo()` can be as a path to a .fasta formatted file, or as
a `Biostrings::XStringSet` object. Unlike other memes functions, `runFimo()`
**does not** accept a `Biostrings::BStringSetList` as input. This is to simplify
ranged join operations (see [joins](#joins)) on output data.

By default, `runFimo()` will parse genomic coordinates from sequence entries
from the fasta headers. These are generated automatically if using
`get_sequences()` to generate sequences for input from a `GRanges` object.

```{r}
data("example_chip_summits", package = "memes")

dm.genome <- BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3

# Take 100bp windows around ChIP-seq summits
summit_flank <- example_chip_summits %>% 
  plyranges::anchor_center() %>% 
  plyranges::mutate(width = 100) 

# Get sequences in peaks as Biostring::BStringSet
sequences <- summit_flank %>% 
  get_sequence(dm.genome)

# get_sequence includes genomic coordinate as the fasta header name
names(sequences)[1:2]
```


### Motif Inputs:
Motif input to `runFimo()` can be as a path to a .meme formatted file, a
list of `universalmotif` objects, or a singular `universalmotif` object.
`runFimo()` will not use any of the default search path behavior for a motif
database as in `runAme()` or `runTomTom()`.

```{r}
e93_motif <- MotifDb::MotifDb %>% 
  # Query the database for the E93 motif using it's gene name
  MotifDb::query("Eip93F") %>% 
  # Convert from motifdb format to universalmotif format
  universalmotif::convert_motifs() %>% 
  # The result is a list, to simplify the object, return it as a single universalmotif
  .[[1]]

# Rename the motif from it's flybase gene number to a more user-friendly name
e93_motif["name"] <- "E93_FlyFactor"
```

## Note about default settings 
`runFimo()` is configured to use different default behavior relative to the
commandline and MEME-Suite Server versions. By default, `runFimo()` runs using
`text` mode, which greatly increases speed and allows returning all detected
matches to the input motifs. By default `runFimo()` will add the sequence of the
matched region to the output data; however, this operation can be very slow for
large sets of regions and can drastically increase the size of the output data.
To speed up the operation and decrease data size, set `skip_matched_sequence =
TRUE`. Sequence information can be added back later using `add_sequence()`. 

```{r}
fimo_results <- runFimo(sequences, e93_motif)
```


## Data integration with join operations {#joins}
The `plyranges` package provides an extended framework for performing
range-based operations in R. While several of its utilities are useful for
range-based analyses, the `join_` functions are particularly useful for
integrating FIMO results with input peak information. A few common examples are
briefly highlighted below:

`plyranges::join_overlap_left()` can be used to add peak-level metadata to motif position information:
```{r}
fimo_results_with_peak_info <- fimo_results %>% 
  plyranges::join_overlap_left(summit_flank)

fimo_results_with_peak_info[1:5]
```

`plyranges::intersect_()` can be used to simultaneously subset input peaks to
the ranges overlapping motif hits while appending motif-level metadata to each
overlap.
```{r}
input_intersect_hits <- summit_flank %>% 
  plyranges::join_overlap_intersect(fimo_results)

input_intersect_hits[1:5]
```

## Identifying matched sequence {#matched-sequence}

When setting `skip_match_sequence = TRUE`, FIMO does not automatically return
the matched sequence within each hit. These sequences can be easily recovered in
R using `add_sequence()` on the FIMO results `GRanges` object.
```{r}
fimo_results_with_seq <- fimo_results %>% 
  plyranges::join_overlap_left(summit_flank) %>% 
  add_sequence(dm.genome)
```

Returning the sequence of the matched regions can be used to re-derive PWMs from different match categories as follows (here done for different binding categories):
```{r}
motifs_by_binding <- fimo_results_with_seq %>% 
  # Split on parameter of interest
  split(mcols(.)$peak_binding_description) %>% 
  # Convert GRangesList to regular list() to use `purrr`
  as.list() %>% 
  # imap passes the list entry as .x and the name of that object to .y
  purrr::imap(~{
    # Pass the sequence column to create_motif to generate a PCM
    create_motif(.x$sequence, 
                 # Append the binding description to the motif name
                 name = paste0("E93_", .y))
    })
```

Motifs from each category can be visualized with `universalmotif::view_motifs()`
```{r}
motifs_by_binding %>% 
  view_motifs()
```

To allow better comparison to the reference motif, we can append it to the list as follows:
```{r}
motifs_by_binding <- c(
  # Add the E93 FlyFactor motif to the list as a reference
  list("E93_FlyFactor" = e93_motif),
  motifs_by_binding
)
```

Visualizing the motifs as ICMs reveals subtle differences in E93 motif sequence
between each category.
```{r}
motifs_by_binding %>% 
  view_motifs()
```

Visualizing the results as a position-probability matrix (PPM) does a better job
of demonstrating that the primary differences between each category are coming
from positions 1-4 in the matched sequences.
```{r}
motifs_by_binding %>% 
  view_motifs(use.type = "PPM")
```

Finally, the sequence-level information can be used to visualize all sequences
and their contribution to the final PWM using `plot_sequence_heatmap`.

```{r, fig.height=5, fig.width=3}
plot_sequence_heatmap(fimo_results_with_seq$sequence)
```


## Importing Data from previous FIMO Runs

`importFimo()` can be used to import an `fimo.tsv` file from a previous run
on the MEME server or on the commandline. Details for how to save data from the
FIMO webserver are below.

### Saving data from FIMO Web Server

To download TSV data from the FIMO Server, right-click the FIMO TSV output link
and "Save Target As" or "Save Link As" (see example image below), and save as
`<filename>.tsv`. This file can be read using `importFimo()`. 

![](save_fimo.png)

# Citation

memes is a wrapper for a select few tools from the MEME Suite, which were
developed by another group. In addition to citing memes, please cite the MEME
Suite tools corresponding to the tools you use.

If you use `runFimo()` in your analysis, please cite:

Charles E. Grant, Timothy L. Bailey, and William Stafford Noble, "FIMO: Scanning for occurrences of a given motif", Bioinformatics, 27(7):1017-1018, 2011. [full text](http://bioinformatics.oxfordjournals.org/content/early/2011/02/16/bioinformatics.btr064.full)

## Licensing Restrictions
The MEME Suite is free for non-profit use, but for-profit users should purchase a
license. See the [MEME Suite Copyright Page](http://meme-suite.org/doc/copyright.html) for details.

# Session Info
```{r}
sessionInfo()
```
---
title: "Install MEME"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Install MEME}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")

eval_vignette <- NOT_CRAN & memes::meme_is_installed()

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  purl = eval_vignette,
  eval = eval_vignette
)
```

# See package website for full vignette

The Bioconductor build system does not have the MEME Suite installed, therefore
these vignettes will not contain any R output. To view the full vignette, visit
this article page on the memes website [at this link](https://snystrom.github.io/memes-manual/articles/install_guide.html)

# Introduction

memes is an R interface to the [MEME Suite](http://meme-suite.org/) family of tools,
which provides several utilities for performing motif analysis on DNA, RNA, and
protein sequences. It works by detecting a local install of the MEME suite,
running the commands, then importing the results directly into R.

## Installing the MEME Suite

memes relies on a local install of the [MEME Suite](http://meme-suite.org/).
For installation instructions for the MEME suite, see the [MEME Suite Installation Guide](http://meme-suite.org/doc/install.html?man_type=web).

Briefly, the MEME suite can be installed to a default location (`~/meme/`) on
Linux, MacOS, Cygwin, and Windows Linux Subsystem using the following shell
commands:

```{bash, eval=F}
# As of December 2021, version 5.4.1 is the most recent MEME-Suite version
# Please check the install guide (linked above) for more recent information
version=5.4.1
wget http://meme-suite.org/meme-software/$version/meme-$version.tar.gz
tar zxf meme-$version.tar.gz
cd meme-$version
./configure --prefix=$HOME/meme --with-url=http://meme-suite.org/ --enable-build-libxml2 --enable-build-libxslt
make
make test
make install
```

For additional troubleshooting or to learn more about install configuration, please see the [Installation Guide](http://meme-suite.org/doc/install.html?man_type=web).

## Detecting the MEME Suite

memes needs to know the location of the `meme/bin/` directory on your local machine.
You can tell memes the location of your MEME suite install in 4 ways. memes
will always prefer the more specific definition if it is a valid path. Here they
are ranked from most- to least-specific:

1. Manually passing the install path to the `meme_path` argument of all memes functions
2. Setting the path using `options(meme_bin = "/path/to/meme/bin/")` inside your R script
3. Setting `MEME_BIN=/path/to/meme/bin/` in your `.Renviron` file, or `export MEME_BIN=/path/to/meme/bin` in your `~/.bashrc`
4. memes will try the default MEME install location `~/meme/bin/`

If memes fails to detect your install at the specified location, it will fall
back to the next option.

To verify memes can detect your MEME install, use `check_meme_install()` which
uses the search herirarchy above to find a valid MEME install. It will report
whether any tools are missing, and print the path to MEME that it sees. This can
be useful for troubleshooting issues with your install.

```{r check_install_works}
library(memes)

# Verify that memes detects your meme install
# (returns all green checks if so)
check_meme_install()
```

```{r check_install_fails}
# You can manually input a path to meme_path
# If no meme/bin is detected, will return a red X
check_meme_install(meme_path = "bad/path")
```


## FAQS

I get the following error: installation of package 'R.oo' had non-zero exit status

 - Problem: Your R installation likely lacks the `R.css` file
 - Solution: when installing the package, set `remotes::install_github("snystrom/memes", INSTALL_opts = c("--no-html"))`
 - NOTE: all help documents for memes will be parsed as plain-text so will lack links or other formatting.

# Session Info
```{r}
sessionInfo()
```
---
title: "Motif Comparison using TomTom"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Motif Comparison using TomTom}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# See package website for full vignette

The Bioconductor build system does not have the MEME Suite installed, therefore
these vignettes will not contain any R output. To view the full vignette, visit
this article page on the memes website [at this link](https://snystrom.github.io/memes-manual/articles/core_tomtom.html)

# Introduction
TomTom is a tool for comparing motifs to a known set of motifs. It takes as
input a set of motifs and a database of known motifs to return a ranked list of
the significance of the match between the input and known motifs. TomTom can be
run using the `runTomTom()` function.

```{r setup}
library(memes)
library(magrittr)
```

## Accepted database formats

`runTomTom()` can accept a variety of inputs to use as the "known" motif database. The formats are as follows:
 - a path to a .meme format file (eg `"fly_factor_survey.meme"`)
 - a list of universalmotifs
 - the output object from `runDreme()`
 - a `list()` of all the above. If entries are named, `runTomTom()` will use those names as the database identifier

### Setting a default database 

memes can be configured to use a default .meme format file as the query
database, which it will use if the user does not provide a value to `database`
when calling `runTomTom()`. The following locations will be searched in order:

1. The `meme_db` option, defined using `options(meme_db = "path/to/database.meme")`
 - The `meme_db` option can also be set to an R object, like a universalmotif list.
2. The `MEME_DB` environment variable defined in `.Renviron`
 - The `MEME_DB` variable will only accept a path to a .meme file

**NOTE:** if an invalid location is found at one option, `runTomTom()` will fall
back to the next location if valid (eg if the `meme_db` option is set to an
invalid file, but the `MEME_DB` environment variable is a valid file, the
`MEME_DB` path will be used.

```{r}
options(meme_db = system.file("extdata/flyFactorSurvey_cleaned.meme", package = "memes", mustWork = TRUE))
```
 
 
## Input types

To use TomTom on existing motifs, `runTomTom()` will accept any motifs in
`universalmotif` format. The `universalmotif` package provides several utilities
for importing data from various sources.
```{r, eval=F}
library(universalmotif)
example_motif <- create_motif("CCRAAAW")
runTomTom(example_motif)
```

`runTomTom()` can also take the output of `runDreme` as input. This allows users
to easily discover denovo motifs, then match them to as set of known motifs.
When run on the output of `runDreme`, all `runTomTom()` output columns will be
appended to the `runDreme()` output data.frame, so no information will be lost. 
```{r, eval=F}
data("dreme_example")
runTomTom(dreme_example) 
```

## Output data
```{r}
# This is a pre-build dataset packaged with memes 
# that mirrors running:
# options(meme_db = system.file("inst/extdata/db/fly_factor_survey_id.meme", package = "memes"))
# example_motif <- create_motif("CCRAAAW")
# example_tomtom <- runTomTom(example_motif)
data("example_tomtom")
```

When run using a `universalmotif` object as input, `runTomTom` returns the following columns:
```{r}
names(example_tomtom)
```

Columns preappended with `best_` indicate the data corresponding to the best match to the motif listed in `name`.

The `tomtom` column is a special column which contains a nested `data.frame` of
the rank-order list of TomTom hits for the motif listed in `name`. 
```{r}
names(example_tomtom$tomtom[[1]])
```

The `best_match_motif` column contains the universalmotif representation of the best match motif.
```{r, fig.height=2, fig.width=5}
library(universalmotif)
view_motifs(example_tomtom$best_match_motif)
```

The `match_motif` column of `tomtom` contains the universalmotif format motif
from the database corresponding to each match in descending order.
```{r, fig.height=3, fig.width=5}
example_tomtom$tomtom[[1]]$match_motif[1:2] %>% 
  view_motifs()
```


The `drop_best_match()` function drops all the `best_match_*` columns from the `runTomTom()` output.
```{r}
example_tomtom %>% 
  drop_best_match() %>% 
  names
```

To unnest the `tomtom` data.frame column, use `tidyr::unnest()`. The
`drop_best_match()` function can be useful when doing this to clean up the
unnested data.frame.
```{r}
unnested <- example_tomtom %>% 
  drop_best_match() %>% 
  tidyr::unnest(tomtom) 
names(unnested)
```

To re-nest the tomtom results, use `nest_tomtom()` (Note: that `best_match_`
columns will be automatically updated based on the rank-order of the `tomtom`
data.frame)

```{r}
unnested %>% 
  nest_tomtom() %>% 
  names
```

### Manipulating the assigned best match

While TomTom can be useful for limiting the search-space for potential true
motif matches, often times the default "best match" is not the correct
assignment. Users should use their domain-specific knowledge in conjunction with
the data returned by TomTom to make this judgement (see below for more details).
memes provides a few convenience functions for reassigning these values.

First, the `update_best_match()` function will update the values of the
`best_match*` columns to reflect the values stored in the first row of the
`tomtom` data.frame entry. This means that **the rank of the `tomtom` data is
meaningful**, and users should only manipulate it if intending to create
side-effects.

If the user can force motifs to contain a certain motif as their best match
using the `force_best_match()` function. `force_best_match()` takes a named
vector as input, where the name corresponds to the input motif `name`, and the
value corresponds to a `match_name` found in the `tomtom` list data (**NOTE:**
this means that users cannot force the best match to be a motif that TomTom did
not return as a potential match).

For example, below the example motif could match either "Eip93F_SANGER_10", or "Lag1_Cell".
```{r}
example_tomtom$tomtom[[1]] %>% head(3)
```

The current best match is listed as "Eip93F_SANGER_10".
```{r}
example_tomtom %>% 
  dplyr::select(name, best_match_name)
```

To force "example_motif" to have the best match as "Lag1_Cell", do the following:
```{r}
new_tomtom <- example_tomtom %>% 
  # multiple motifs can be updated at a time by passing additional name-value pairs.
  force_best_match(c("example_motif" = "Lag1_Cell"))
```

The `best_match_*` columns will be updated to reflect the modifications.
```{r}
# original best match:
example_tomtom$best_match_name
# new best match after forcing:
new_tomtom$best_match_name
```
 
## Visualize data

`view_tomtom_hits()` can be used to compare the hits from tomtom to each input
motif. Hits are shown in descending order by rank. By default, all hits are
shown, or the user can pass an integer to `top_n` to view the top number of
motifs. This can be a useful plot for determining which of the matches appear to
be the "best" hit. 

For example, it appears that indeed "Eip93F_SANGER_10" is the best of the top 3
hits, as most of the matching sequences in the "Lag1_Cell" and "pho_SOLEXA_5"
motifs correspond to low information-content regions of the matched motifs.
```{r, fig.height=4, fig.width=5}
example_tomtom %>% 
  view_tomtom_hits(top_n = 3)
```

## Importing previous data

`importTomTomXML()` can be used to import a `tomtom.xml` file from a previous
run on the MEME server or on the commandline. Details for how to save data from
the TomTom webserver are below.

### Saving data from TomTom Web Server
To download XML data from the MEME Server, right-click the TomTom XML output link
and "Save Target As" or "Save Link As" (see example image below), and save as
`<filename>.xml`. This file can be read using `importTomTomXML()`

![](save_tomtom.png)

# Citation

memes is a wrapper for a select few tools from the MEME Suite, which were
developed by another group. In addition to citing memes, please cite the MEME
Suite tools corresponding to the tools you use.

If you use `runTomTom()` in your analysis, please cite:

Shobhit Gupta, JA Stamatoyannopolous, Timothy Bailey and William Stafford Noble, "Quantifying similarity between motifs", Genome Biology, 8(2):R24, 2007. [full text](http://genomebiology.com/2007/8/2/R24)

## Licensing Restrictions
The MEME Suite is free for non-profit use, but for-profit users should purchase a
license. See the [MEME Suite Copyright Page](http://meme-suite.org/doc/copyright.html) for details.

# Session Info
```{r}
sessionInfo()
```
---
title: "Tidying Motif Metadata"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Tidying Motif Metadata}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

memes uses the `universalmotif` package to simplify working with motif metadata.
`universalmotif` objects can be represented in an alternative form, the
`unviversalmotif_df` which allows users to manipulate motif metadata just as
they would a normal R data.frame (in fact, they are just R data.frames).

These objects are useful for tidying motif metadata to prepare a motif database
for use with memes, or performing any other data-driven tasks involving motifs.
Here I describe one way these data structures can be used to construct a motif
database for use with memes.

```{r setup}
library(memes)
library(magrittr)
library(universalmotif)
```

The `MotifDb` package makes it easy to query thousands of motifs from public
databases. Here, I will describe how to use one of these queries as input to
memes functions, and how to manipulate the resulting motifs to prepare them for
MEME suite tools.

I will use the motifs from the [FlyFactorSurvey](https://mccb.umassmed.edu/ffs/)
as an example. They can be accessed from `MotifDb` using the following query.
```{r}
flyFactorDb <- MotifDb::MotifDb %>% 
  MotifDb::query("FlyFactorSurvey")
```

Use `universalmotif::convert_motifs()` to convert a `MotifDb` query into motif
objects. In many cases, the resulting list can be used directly as input to
memes functions, like `runTomTom`, or `runAme`.
```{r}
flyFactorMotifs <- flyFactorDb %>% 
  convert_motifs()
```

But there are some issues with this database. For example, the following motif
name is the FlyBase gene number, and the alternate name is the actual
informative name of the PWM. The MEME Suite relies more heavily on the primary
name, so it would be nice if the database used interpretable names.
```{r}
flyFactorMotifs %>% 
  head(1)
```

The universalmotif function `to_df()` converts universalmotif lists into
`universalmotif_df` format which can be used to update motif entries. This is
particularly useful when dealing with several motifs at once.
```{r}
flyFactor_data <- flyFactorMotifs %>% 
  to_df()
```

The columns of the `universalmotif_df` can be freely changed to edit the
properties of the motifs stored in the `motif` column. Just like standard
data.frames, additional columns can be added to store additional metadata. For
more details on these objects see the help page: `?universalmotif::to_df`.

```{r}
# The following columns can be changed to update motif metadata
flyFactor_data %>% 
  names
```


using the `universalmotif_df`, we can quickly see that the issue with FBgn
numbers only applies to certain entries. And TFs which are represented by
multiple motifs in the database are assigned the same name. The MEME Suite tools
which use a motif database (like TomTom and AME) require that the entries have
unique primary identifiers, therefore the default names will not be appropriate.

```{r}
flyFactor_data %>% 
  head(5)
```
However, the `altname` slots from the `motifDb` query are already unique, so we can make them the primary name.
```{r}
length(flyFactor_data$altname) == length(unique(flyFactor_data$altname))
```

An easy way is to use `dplyr::rename` to swap the columns.
```{r}
flyFactor_data %<>% 
  dplyr::rename("altname" = "name", 
                "name" = "altname")
```

The `name` column now contains the full motif name. 
```{r}
flyFactor_data %>% 
  head(3)
```






Next to solve the issue with the FBgn's. FBgn numbers are unique identifiers for
a gene within a given FlyBase reference assembly. However, FBgn numbers are not
stable over time (i.e. the same gene may have a different FBgn number between
reference assemblies), therefore they are unreliable values to determine the
correct gene symbol. [FlyBase.org](https://flybase.org) has a nice [conversion tool](https://flybase.org/convert/id) 
which can be used to update FBgn numbers. 

As of this writing in March 2021, the FBgn entries provided by the Fly Factor
Survey database are out of date. In order to demonstrate an example of methods
for tidying motif metadata, I won't use the FlyBase conversion tool, but will
instead highlight some approaches which may be more generally useful when
working with motif databases from disparate sources.

For this example, we will try to grab the correct gene name from the motif name,
which is stored in the first field of the name, formatted as follows:
"<gene>_<sequencing-platform>_<FBgn>".

We use `tidyr::separate` to split out the first entry to the `tifd` column, then
only use this value if the altname contains an FBgn.
```{r}
flyFactor_data %<>% 
  # Critical to set remove = FALSE to keep the `name` column
  tidyr::separate(name, c("tfid"), remove = FALSE, extra = "drop") %>% 
  # Only use the tfid if the altname contains an FBgn
  dplyr::mutate(altname = ifelse(grepl("^FBgn", altname), tfid, altname))
```

Now, the first two entries are listed as "ab" instead of "FBgn0259750".
```{r}
flyFactor_data %>% 
  head(3)
```

Next, because the FBgn's are out of date, we will remove them from the "names"
to shorten up the motif names. This also makes the motif name more comparable to
the original motif names from the [FlyFactor Survey](https://mccb.umassmed.edu/ffs/).
```{r}
flyFactor_data %<>% 
  dplyr::mutate(name = gsub("_FBgn\\d+", "", name))
```

## A few reality checks

It's worth taking a look at the instances where the `altname` and our parsed
`tfid` do not match. This is a good way to ensure we haven't missed any
important edge cases in the data. As new edge cases are encountered, we can develop new rules for tidying the data to ensure a high quality set of motifs.

Start by simply filtering for all instance where there is a mismatch between
`altname` and `tfid`.

Carefully compare the `altname`, `name`, and `tfid` columns. Why might the
values differ? Are there instances that make you question the data?
```{r}
flyFactor_data %>% 
  dplyr::filter(altname != tfid) %>% 
  # I'm only showing the first 5 rows for brevity, but take a look at the full
  # data and see what patterns you notice
  head(5)
```

One thing that becomes obvious is that many motifs have mismatched
`altname`/`tfid` values because of capitalization or hyphenation differences.
You can use domain-specific knowledge to assess which one is correct. For
*Drosophila*, "abd-A" is correct over "AbdA", for example.

After manually inspecting these rows, I determined that instances of different
capitalization, hyphenation, or names that contain "." or "()" can be ignored.
To further investigate the data, I will ignore capitalization and special character
differences as follows:

```{r}
flyFactor_data %>% 
  # calling tolower() on both columns removes capitalization as a difference
  dplyr::filter(tolower(altname) != tolower(tfid),
                # Select all altnames that do not contain "-", "." or "("
                !grepl("-|\\.|\\(", altname),
                ) %>% 
  # I'll visalize only these columns for brevity
  dplyr::select(altname, tfid, name, consensus) %>% 
  head(10)
 
```
Next, what is obvious is that several `altnames` set to "da" have a high number
of mismatched `tfid`s. For instance, `amos_da_SANGER_10`. When checking the
[FlyFactorSurvey page for
da](https://mccb.umassmed.edu/ffs/TFdetails.php?FlybaseID=FBgn0000413), it
reveals only 1 motif corresponds to this factor. Checking the [page for amos](https://mccb.umassmed.edu/ffs/TFdetails.php?FlybaseID=FBgn0003270) shows a
match to `amos_da_SANGER_10`. Therefore, we can conclude that factors assigned
the name of `da` are incorrectly assigned, and we should prefer our parsed
`tfid`.

```{r}
flyFactor_data %<>% 
  # rename all "da" instances using their tfid value instead
  dplyr::mutate(altname = ifelse(altname == "da", tfid, altname))
```

Now we've handled the "da" mismatches, we filter them out to identify new special cases.
```{r}
flyFactor_data %>% 
  dplyr::filter(tolower(altname) != tolower(tfid),
                !grepl("-|\\.|\\(", altname)) %>% 
  dplyr::select(altname, tfid, name, consensus) %>% 
  head(10)
```
The next thing to notice about these data is that entries with "CG" prefixed
tfids are often mismatched. This is because when the FlyFactor survey was
conducted, many genes were unnamed, and thus assigned a CG from FlyBase. As time
has gone on, some CG's have been named. Checking the [FlyBase page for
CG10267](http://flybase.org/reports/FBgn0037446) reveals that it has been
renamed "Zif". This matches with the `altname`, so we conclude that rows with a
"CG" `tfid` can be safely skipped as their `altname` contains the new gene symbol.

```{r}
flyFactor_data %>% 
  dplyr::filter(tolower(altname) != tolower(tfid),
                !grepl("-|\\.|\\(", altname),
                # Remove CG genes from consideration
                !grepl("CG\\d+", tfid)
                ) %>% 
  dplyr::select(altname, tfid, name, consensus)
```
The remaining rows (only 20 values) can be manually inspected for any
discrepancies. I went through each entry by hand, looking up their motifs on
[FlyFactor](https://mccb.umassmed.edu/ffs/), and their gene names on
[FlyBase](flybase.org) to determine the best way to handle these motifs.
Sometimes the best way to be sure your data are high quality is to carefully inspect it!

I determined from this that a few `altnames` need swapping, and one motif I will
remove because it is unusual
([Bgb](https://mccb.umassmed.edu/ffs/TFdetails.php?FlybaseID=FBgn0013753) has an
identical motif to
[run](https://mccb.umassmed.edu/ffs/TFdetails.php?FlybaseID=FBgn0003300), but
the motif is marked "run" on the FlyFactor website).

I'll make those changes to the data:
```{r}
swap_alt_id <- c("CG6272", "Clk", "Max", "Mnt", "Jra")
remove <- "Bgb"

flyFactor_data %<>% 
  dplyr::mutate(altname = ifelse(altname %in% swap_alt_id, tfid, altname)) %>% 
  dplyr::filter(!(altname %in% remove))
```

Finally, the remaining motif metadata is also OK based on my manual inspection.
```{r}
flyFactor_data %>% 
  dplyr::filter(tolower(altname) != tolower(tfid),
                !grepl("-|\\.|\\(", altname),
                # Remove CG genes from consideration
                !grepl("CG\\d+", tfid)
                ) %>% 
  dplyr::select(altname, tfid, name, consensus)
```

## Removing duplicate motif matrices

Just because the metadata for each entry is unique, this does not mean that the
motif matrix for each entry is unique. There are many reasons why two different
factors could have identical motifs: some biological, others technical. In the
case of the FlyFactorSurvey, some entries are duplicated in MotifDb which should not be. 

For instance, the following motif is a duplicate where the tidied metadata
matches: 

```{r}
flyFactor_data %>% 
  dplyr::filter(consensus == "MMCACCTGYYV")
```
It is difficult to determine in a high-throughput way whether any matrix entries
are identical in a large database, and it is not always possible to rely on
metadata to determine matrix duplication.

In order to identify and remove duplicate motif matrices, memes provides
`remove_duplicate_motifs()`, which can be used to deduplicate a list of motifs
based solely on their motif matrices (i.e. it ignores motif name & other
metadata). We will use this strategy to deduplicate the flyFactor data.

(NOTE: When working with other motif databases, it is critical to understand the data
source to determine appropriate measures for handling duplicated entries.)
```{r}
# This operation takes a while to run on large motif lists
flyFactor_dedup <- remove_duplicate_motifs(flyFactor_data)
```

Duplicate removal identifies and removes `r nrow(flyFactor_data) - nrow(flyFactor_dedup)` identical matrices.
```{r}
# Rows before cleanup
nrow(flyFactor_data)
# Rows after cleanup
nrow(flyFactor_dedup)
```

Using the example from before now shows only 1 motif corresponding to this sequence.
```{r}
flyFactor_dedup %>% 
  dplyr::filter(consensus == "MMCACCTGYYV")
```
Finally, now that the database has been tidied and deduplicated, the resulting
data.frame can be converted back into a universalmotif list using
`to_list()`. To discard the additional columns we created so they are not passed on to the `universalmotif`, set `extrainfo = FALSE`.
```{r}
# extrainfo = FALSE drops the extra columns we added during data cleaning which are now unneeded
flyFactorMotifs_final <- to_list(flyFactor_dedup, extrainfo = FALSE)
```

The resulting universalmotif list object now reflects the changes we made to the
`data.frame` and can now be exported as a .meme format file using
`universalmotif::write_meme` or can be used directly as input to tools like
`runTomTom` or `runAme`.
```{r}
flyFactorMotifs_final %>% 
  head(1)
```
This cleaned-up version of the FlyFactorSurvey data is packaged with memes in
`system.file("extdata/flyFactorSurvey_cleaned.meme", package = "memes")`.
```{r,eval=F,include=T}
write_meme(flyFactorMotifs_final, "flyFactorSurvey_cleaned.meme")
```

# Session Info
```{r}
sessionInfo()
```---
title: "Motif Enrichment Testing using AME"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Motif Enrichment Testing using AME}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include=FALSE}
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")

eval_vignette <- NOT_CRAN & memes::meme_is_installed()

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  purl = eval_vignette,
  eval = eval_vignette
)
```

```{r setup}
library(memes)
suppressPackageStartupMessages(library(GenomicRanges))
library(magrittr)
```

## Sequence Input
AME requires a series of input sequences to scan for motif enrichment.
`runAme()` accepts sequence input in the following formats:

 - a `Biostrings::XStringSet` object
 - a named list of `Biostrings::XStringSet` objects
 - a path to a .fasta file
 
**NOTE** `XStringSet` inputs can be easily generated for DNA sequences from a
GRanges object using the `get_sequence()` function

```{r}
data("example_peaks", package = "memes")

dm.genome <- BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3

sequence <- example_peaks %>% 
  get_sequence(dm.genome)
```

## Database Input
AME scans input sequences against a database of known motifs and tests for
enrichment of each motif in the database. `runAme()` can accept a database in
the following formats:

 - a list of universalmotif objects
 - a single universalmotif object
 - the results object from `runDreme`
 - a path to a .meme format file

### Setting a default database 

memes can be configured to use a default .meme format file as the query
database, which it will use if the user does not provide a value to `database`
when calling `runAme()`. The following locations will be searched in order.

1. The `meme_db` option, defined using `options(meme_db = "path/to/database.meme")`
 - The `meme_db` option can also be set to an R object, like a universalmotif list.
2. The `MEME_DB` environment variable defined in `.Renviron`
 - The `MEME_DB` variable will only accept a path to a .meme file

**NOTE:** if an invalid location is found at one option, `runAme()` will fall
back to the next location if valid (eg if the `meme_db` option is set to an
invalid file, but the `MEME_DB` environment variable is a valid file, the
`MEME_DB` path will be used.

```{r}
options(meme_db = system.file("extdata/flyFactorSurvey_cleaned.meme", package = "memes", mustWork = TRUE))
```

## Running AME

`runAme()` supports running AME using three modes:

| AME Mode          | Description                | Command                                         |
|:-----------------:|:--------------------------:|:------------------------------------------------|
| Vs Shuffled       | Input vs Shuffled Sequence | `runAme(input = sequence, control = "shuffle")` |
| Discriminative    | Input vs Control Sequence  | `runAme(input = sequence, control = control)`   |
| Partitioning      | Rank Input by fasta score  | `runAme(input = sequence, control = NA)`        |

```{r, eval=F}
ame_vs_shuffle <- runAme(sequence)
```

```{r, eval=F}
ame_vs_control <- runAme(sequence[1:5], sequence[6:10])
```


To run AME using partitioning mode, the fasta header must contain a score value
for each entry in the form: ">entry_name score". The `get_sequences()` `score`
argument allows users to set the score value to a column value from input
regions.

```{r}
sequence_scored <- example_peaks %>% 
  plyranges::mutate(score = seq_along(.)) %>% 
  get_sequence(dm.genome, score = "score")

names(sequence_scored)[1]
```
```{r, eval=F}
ame_partition <- runAme(sequence_scored, control = NA)
```

### Running AME on multiple groups

If using a list input to `runAme()`, it will dispatch multiple AME runs for each
object in the list.

```{r}
data("example_chip_summits", package = "memes")

seq_by_behavior <- example_chip_summits %>% 
  plyranges::mutate(width = 100) %>% 
  split(mcols(.)$e93_sensitive_behavior) %>% 
  get_sequence(dm.genome)
```
```{r, eval=F}
ame_by_behavior <- runAme(seq_by_behavior)
```

#### Discriminative analysis using list input

If the input to `runAme()` is a named list of `XStringSet` objects, `control`
can be set to one or more values from `names(input)` to use those regions as
background. It will skip running those regions as the input. The following code
will result in these comparisons:

1. Increasing vs Static
2. Decreasing vs Static
```{r}
ame_by_behavior_vs_static <- runAme(seq_by_behavior, control = "Static")
```

If multiple names are used in the `control` section, they will be combined
together to make a single control set which will be used for all comparisons.
Here, we use "Static" and "Decreasing" sites as the control, which will result
in only running 1 comparison: Increasing vs Static+Decreasing.
```{r, eval=F}
runAme(seq_by_behavior, control = c("Static", "Decreasing"))
```
```{r, echo=F}
data("example_ame", package = "memes")
ame_by_behavior_vs_static <- example_ame
```

## Output Format

AME will return different output formats depending on the `method` used. For
detailed information about these values see the [AME Output description
webpage](http://meme-suite.org/doc/ame-output-format.html). As a general rule of
thumb, `runAme()` will return the same column names described in the webpage,
except dashes are removed and all column names are lowercase.

```{r}
ame_by_behavior_vs_static$Decreasing %>% 
  names
```

If `runAme()` is run with `method = "fisher"`, the sequences output can be added
to the results by setting `sequences = TRUE`. This will be added as a list
column named `sequences` that can be unnested using `tidyr::unnest()`.

## Visualizing Results as Heatmap

The `plot_ame_heatmap()` function provides a method to easily generate
visualizations of AME results.

To plot results from multiple runs together, they must first be joined into 1
data frame. The `ame_by_behavior_vs_static` object is a list whose names correspond to the
E93 response (Increasing or Decreasing). The list can be combined into a data.frame using
`dplyr::bind_rows`. Setting `.id = "behavior` creates a new column
`behavior` that contains the names from the `ame_by_behavior_vs_static` list. In this
way, the resulting data.frame contains all AME results for each run, which can
be distinguished by the `behavior` column.

```{r}
ame_by_behavior_vs_static %>% 
  # AME results in list format are easily combined using dplyr::bind_rows
  # .id will specify a column to hold the list object names
  dplyr::bind_rows(.id = "behavior") %>% 
  # setting group to a column name will split the results on the y-axis
  plot_ame_heatmap(group = behavior)
```

### Complex Heatmap Example

There are several nuances when making heatmap visualizations of these data. The
following examples highlight some of these issues and provide alternative
approaches and solutions.

We start by using different binding site categories as input.
```{r}
seq_by_binding <- example_chip_summits %>% 
  plyranges::mutate(width = 100) %>% 
  split(mcols(.)$peak_binding_description) %>% 
  get_sequence(dm.genome)
```

```{r, eval=F}
ame_by_binding <- seq_by_binding %>% 
  runAme
```

```{r, echo=F}
# Allows vignette to build on systems w/o functioning MEME suite
data("example_ame_large", package = "memes")
ame_by_binding <- example_ame_large
```

```{r}
ame_res <- ame_by_binding %>% 
  dplyr::bind_rows(.id = "binding_type")
```

It is possible to aggregate results from multiple runs into a heatmap by setting
the `group` parameter in `plot_ame_heatmap()`.

This is too many hits to properly view in this vignette, but you can see that the heatmap
will plot motifs by their overlap across groups, where unique motifs are on the
left, and shared motifs are on the right.
```{r, fig.height=5,fig.width=15}
ame_res %>% 
  plot_ame_heatmap(group = binding_type)
```

## Issues with Heatmap Visualization

The dynamic range of p-values in these data varies between groups. For this
reason, a simple heatmap scaled using all data values will make it more
difficult to interpret within groups with a lower dynamic range of values. In
other words, because the dynamic range of values are different between
experiments, placing them on the default scale for comparison may not always be
the most optimal visualization.

We can partially overcome this limitation by filling the heatmap with the
normalized rank value for each TF, which accounts for differences in total number of
discovered motifs between AME runs. Although it does not completely abrogate
differences, the signal values for high-ranked motifs within groups will be more
comparable. However, **the normalized rank visualization eliminates all real values related to statistical significance!** 
Instead, this visualization represents the relative ranks of hits within an AME
run, which already pass a significance threshold set during `runAME()`. This
means that even if several motifs have similar or even identical p-values, their
heatmap representation will be a different color value based on their ranked
order in the results list. This tends to only be useful when there are a large
number of hits (>=100). Both visualizations can be useful and reveal different
properties of the data to the user. **If in doubt**, prefer the
`-log10(adj.pvalue)` representation.

Below is a comparison of the distribution of values when using
`-log10(adj.pvalue)` (A) vs normalized ranks (B). Because orphan sites tend to
have smaller p-values overall, the heatmap scale will be skewed towards the high
values in the orphan data, making ectopic and entopic heat values lighter by
comparison.
```{r, fig.height=3, fig.width=7.5}
ame_res %>% 
  ame_compare_heatmap_methods(group = binding_type)
```

To use the normalized rank value, set `value = "normalize"` in `plot_ame_heatmap()`.

This plot reveals that the motifs which tend to be shared across all 3
categories tend to be higher ranks in the output than the motifs unique to the
different categories, which tend to come from lower ranks. This *suggests* that
although there are differences in motif content across the three categories,
they may be largely similar in motif makeup. We will investigate this question
in more detail in the "Denovo motif similarity" section.
```{r, fig.height=3, fig.width=25}
library(ggplot2)

(normalize_heatmap <- ame_res %>% 
  dplyr::group_by(binding_type, motif_alt_id) %>% 
  dplyr::filter(adj.pvalue == min(adj.pvalue)) %>% 
  plot_ame_heatmap(group = binding_type, id = motif_alt_id, value = "normalize") +
    # All ggplot functions can be used to extend or edit the heatmap plots
    ggtitle("value = \"normalize\""))
```

An additional third option exists to rescale the `-log10(adj.pvalue)` heatmap to
change the heatmap's maxiumum color value. This allows the user to maintain
values which represent significance, but rescale the data to capture the lower
end of the dynamic range. Using the cumulative distribution plot above, a
reasonable cutoff is anywhere between 7 & 10, which captures > 90% of the data
for ectopic and entopic sites. 

A comparison of all three methods can be seen below. 

```{r, fig.height=3, fig.width=25}
pval_heatmap <- ame_res %>% 
  dplyr::group_by(binding_type, motif_alt_id) %>% 
  dplyr::filter(adj.pvalue == min(adj.pvalue)) %>% 
  plot_ame_heatmap(group = binding_type, id = motif_alt_id) +
    ggtitle("value = -log10(adj.pvalue)")

scale_heatmap <- ame_res %>% 
  dplyr::group_by(binding_type, motif_alt_id) %>% 
  dplyr::filter(adj.pvalue == min(adj.pvalue)) %>% 
  plot_ame_heatmap(group = binding_type, id = motif_alt_id, scale_max = 7.5) +
    ggtitle("value = -log10(adj.pvalue) (scale capped at 7.5)")
```

Below is a comparison using the `-log10(adj.pvalue)` vs `normalize` methods for
plotting the heatmap. Note how the different plots highlight different data
properties. The `-log10(adj.pvalue)` plot shows overall significance of each
hit, while `normalize` method shows the relative rank of each hit within a
`binding_type`. Lowering the maximum scale value in C) does a better job than A)
at visualizing differences in significance along the ectopic and entopic rows at
the cost of decreasing the dynamic range of the orphan row. Selecting a
visualization for publication will depend heavily on context, but if in doubt,
prefer one which includes information of statistical significance as in A) or C).
```{r, fig.height=9, fig.width=25}
cowplot::plot_grid(pval_heatmap, 
                   normalize_heatmap, 
                   scale_heatmap, 
                   ncol = 1, labels = "AUTO")
```

## Importing Previous Data

`importAme()` can be used to import an `ame.tsv` file from a previous run
on the MEME server or on the commandline. Details for how to save data from the
AME webserver are below.

Optionally, if AME was run on the commandline with `--method fisher`, the user
can pass a path to the `sequences.tsv` file to the `sequences` argument of
`importAme()` to append the sequence information to the AME results.

### Saving data from AME Web Server
To download TSV data from the MEME Server, right-click the AME TSV output link
and "Save Target As" or "Save Link As" (see example image below), and save as
`<filename>.tsv`. This file can be read using `importAme()`. 

![saving AME tsv results](save_ame.png)

# Citation

memes is a wrapper for a select few tools from the MEME Suite, which were
developed by another group. In addition to citing memes, please cite the MEME
Suite tools corresponding to the tools you use.

If you use `runAme()` in your analysis, please cite:

Robert McLeay and Timothy L. Bailey, "Motif Enrichment Analysis: A unified framework and method evaluation", BMC Bioinformatics, 11:165, 2010, doi:10.1186/1471-2105-11-165. [full text](http://www.biomedcentral.com/1471-2105/11/165)

## Licensing Restrictions
The MEME Suite is free for non-profit use, but for-profit users should purchase a
license. See the [MEME Suite Copyright Page](http://meme-suite.org/doc/copyright.html) for details.

# Session Info
```{r}
sessionInfo()
```---
title: "Denovo Motif Discovery Using DREME"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Denovo Motif Discovery Using DREME}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include=FALSE}
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")

eval_vignette <- NOT_CRAN & memes::meme_is_installed()

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  purl = eval_vignette,
  eval = eval_vignette
)
```

# See package website for full vignette

The Bioconductor build system does not have the MEME Suite installed, therefore
these vignettes will not contain any R output. To view the full vignette, visit
this article page on the memes website [at this link](https://snystrom.github.io/memes-manual/articles/core_dreme.html)

```{r setup}
library(memes)
```

```{r}
# Verify that memes detects your meme install
# should return all green checks if so.
check_meme_install()
```

```{r}
fa <- system.file("extdata/fasta_ex/fa1.fa", package = "memes")
```

```{r}
# NOTE: setting e > 1 is usually not recomended. 
# the example fasta file only has 1 sequence in it
# to keep the file size low and let the example run quickly.
# I set evalue = 39 because dreme cannot detect high confidence motifs from only 1 sequence.
dreme_out <- runDreme(fa, "shuffle", evalue = 39, outdir = tempdir())
```

[DREME Commandline Documentation](http://meme-suite.org/doc/dreme.html)

### Aliased flags
| memes alias | DREME Flag | description                               |
|:------------:|:----------:|:------------------------------------------|
| nmotifs      | m*         | max number of motifs to discover          |
| sec          | t          | max number of seconds to run              |
| evalue       | e          | max E-value cutoff                        |
| seed         | s*         | random seed if using "shuffle" as control | 
| ngen         | g          | number of REs to generalize               |

\* flags marked with \* must be assigned using their alias
```{r, eval=F}
# equivalent to above
runDreme(fa, "shuffle", evalue = 39, outdir = tempdir())
runDreme(fa, "shuffle", e = 39, outdir = tempdir(), nmotifs = 2)
```

dreme results are a `data.frame`. The `motif` column contains a `universalmotif`
object with the PCM information for each *de-novo* discovered motif. This is so
that any filtering of the results object also simply filter the available
motifs. For more details about each column see the "Value" section of `?runDreme`.

```{r}
dreme_out
```

The results can be converted back to `universalmotif` format using `to_list()`.
The `view_motifs()` function accepts a `universalmotif` list and can be used to
visualize the motifs.
```{r}
library(universalmotif)
library(dplyr)

dreme_out %>% 
  to_list() %>% 
  view_motifs()
```

The primary advantage of using the `data.frame` output allows simple integration
with base subsetting, piping, and the `tidyverse`.
```{r}
dreme_out %>% 
  # after filtering with dplyr, only motifs with length 3 will be plotted
  filter(length == 3) %>% 
  to_list() %>% 
  universalmotif::view_motifs()
```

`universalmotif` manipulations can easily be executed on the motifs as well. For example:
```{r, fig.height=1.5}
dreme_out %>% 
  to_list() %>% 
  merge_motifs() %>% 
  view_motifs()
```

#### Updating motif information
Occasionally, it can be useful to update the metadata associated with a
dicovered motif (for example, to assign a new name to a denovo motif). memes
provides a few utilities to accomplish this.

`update_motifs()` will search for specific column names which describe
properties of the `motif` column and update the metadata in the `motif` column
to reflect those values. See `?update_motifs` for details.

`as_universalmotif()` will convert one of these special universalmotif
data.frames into a universalmotif list after updating the metadata to reflect
values as in `update_motifs()`.
```{r}
# update_motifs will update the values in the motif column
# to values in the data.frame
dreme_edit <- dreme_out %>% 
  dplyr::mutate(name = c("one", "two", "three", "four", "five")) %>% 
  update_motifs()

# to_list() will first update motif information
# before returning only the motif column
edit_motifs <- dreme_out %>% 
  dplyr::mutate(name = c("one", "two", "three", "four", "five")) %>% 
  to_list()

# The following outputs are identical
# where edit_motifs is a list of motifs
# and dreme_edit is a data.frame with a motif list column
identical(edit_motifs$one, dreme_edit$motif$one)
```


### Notes about shuffled control sequences

Setting `control = "shuffle"` will use dreme's random number generator to
shuffle the input sequences. By default, dreme will use `1` as the random seed,
so repeat runs of the same shuffle command will produce the same output. To
change the random seed, pass `seed = [your random seed]` to `runDreme()`.
**NOTE:** beware system-specific differences. As of MEME v5, dreme will compile
using the default python installation on a system (either python2.7 or python3).
The random number generator changed between python2.7 and python3, so results
will not be reproducible between systems using python2.7 vs 3 
**even if setting the same random seed**.

One way to overcome this is to manually shuffle the sequences within R. This can
be done easily using `universalmotif::shuffle_sequences()`. Set `k = 2` to
preserve dinucleotide frequency (similar to dreme's built-in shuffle), and set
`rng.seed` to any number to create a reproducible shuffle. The output of this
function can be used directly as the control sequences.

```{r}
# Create random sequences to use for this example
seq <- create_sequences(rng.seed = 100)
# Shuffle sequences preserving dinucleotide frequency
shuffle <- shuffle_sequences(seq, k = 2, rng.seed = 100)
```

### Analysis on Multiple Groups and Differential Analysis

Often, users want to perform motif analysis on many groups of sequences.
For example, here we have ChIP-seq peaks for a transcription factor, E93.
Analysis of chromatin accessibility in E93 peaks revealed sites that Increase
accessibility, Decrease accessibility, or remain Static following E93 binding.
```{r}
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(plyranges))

data("example_chip_summits", package = "memes")
peaks <- example_chip_summits
```

To examine whether there are differences in motif content between increasing,
decreasing, and static sites, we split the peaks into a list by their response
to E93.
```{r}
by_behavior <- peaks %>% 
  anchor_center() %>% 
  mutate(width = 100) %>% 
  split(mcols(.)$e93_sensitive_behavior)
```

Next, this list can be used directly in `get_sequences()` to generate a list of sequences for each set of peaks.
```{r}
dm.genome <- BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3

seq_by_behavior <- by_behavior %>% 
  get_sequence(dm.genome)
```

```{r}
names(seq_by_behavior)
```

To run DREME on each set using shuffled input sequence as background, run:
```{r, eval = F}
runDreme(seq_by_behavior, control = "shuffle")
```

#### Discriminative analysis using list input
For this analysis, however, we are most interested in identifying motifs
associated with increasing and decreasing that do not involve E93 binding.
Therefore, a more appropriate control is to use the Static sites as
background. 

As always, an `XStringSet` object can be used as the control regions. However, running dreme in this way will run 3 jobs:

1. Increasing vs Static
2. Decreasing vs Static
3. Static vs Static

This will waste time, as job #3 will detect no motifs (since input & control are
identical), but will still take a long time to run. `runDreme()` has additional
functionality to help avoid these issues, and to facilitate more complicated
analysis designs.
```{r, eval=F}
runDreme(seq_by_behavior, control = seq_by_behavior$Static)
```

If the input to `runDreme` is a named list of `XStringSet` objects, `control`
can be set to one or more values from `names(input)` to use those regions as
background. It will skip running those regions as the input. The following code
will result in these comparisons:

1. Increasing vs Static
2. Decreasing vs Static
```{r, eval=F}
runDreme(seq_by_behavior, control = "Static")
```

If multiple names are used in the `control` section, they will be combined
together to make a single control set which will be used for all comparisons.
Here, we use "Static" and "Decreasing" sites as the control, which will result
in only running 1 comparison: Increasing vs Static+Decreasing.
```{r, eval=F}
runDreme(seq_by_behavior, control = c("Static", "Decreasing"))
```

## Importing previous data

`importDremeXML()` can be used to import a `dreme.xml` file from a previous run
on the MEME server or on the commandline. Details for how to save data from the
DREME webserver are below.

### Saving data from DREME Web Server
To download XML data from the MEME Server, right-click the DREME XML output link
and "Save Target As" or "Save Link As" (see example image below), and save as
`<filename>.xml`. This file can be read using `importDremeXML()`

![](save_dreme.png)

# Citation

memes is a wrapper for a select few tools from the MEME Suite, which were
developed by another group. In addition to citing memes, please cite the MEME
Suite tools corresponding to the tools you use.

If you use `runDreme()` in your analysis, please cite:

Timothy L. Bailey, "DREME: Motif discovery in transcription factor ChIP-seq data", Bioinformatics, 27(12):1653-1659, 2011. [full text](https://academic.oup.com/bioinformatics/article/27/12/1653/257754)

## Licensing Restrictions
The MEME Suite is free for non-profit use, but for-profit users should purchase a
license. See the [MEME Suite Copyright Page](http://meme-suite.org/doc/copyright.html) for details.

# Session Info
```{r}
sessionInfo()
```
---
title: "ChIP-seq Analysis"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{ChIP-seq Analysis}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include=FALSE}
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")

eval_vignette <- NOT_CRAN & memes::meme_is_installed()

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  purl = eval_vignette,
  eval = eval_vignette
)
```

```{r setup}
library(memes)
library(magrittr)
library(ggplot2)
suppressPackageStartupMessages(library(GenomicRanges))
```

# See package website for full vignette

The Bioconductor build system does not have the MEME Suite installed, therefore
these vignettes will not contain any R output. To view the full vignette, visit
this article page on the memes website [at this link](https://snystrom.github.io/memes-manual/articles/integrative_analysis.html)

# Introduction
In this vignette, we'll explore using `memes` to deeply analyze a set of
ChIP-seq peaks to identify motifs to explain differences in transcription factor
binding, and consequences to chromatin accessibility at these ChIP peaks.

We will use a dataset from (Nystrom, 2020) which investigated the binding of the
transcription factor E93 during *Drosophila* wing development. These data are
useful because they contain multiple groupings with interesting motif properties
between groups.

The ChIP experiments were performed under two conditions: ectopic expression of
E93 during early development, and under wild-type conditions in a late stage of
development. Peaks are annotated based on whether they are found in both
conditions ("Endogenous" peaks), only in the ectopic expression condition
("Ectopic" peaks), or only found under wild-type conditions ("Orphan" peaks).

A key question from this observation is whether there are motifs which can
distinguish these three categories of binding from eachother. To investigate
this question, we will first look for enrichment of known transcription factor
motifs within each binding category using the AME tool.

To start, ensure memes can detect your local install of the MEME Suite (see
`vignette("install_guide")` for help).

```{r}
check_meme_install()
```

### Prepare peaks for analysis

A subset of annotated E93 ChIP summits are shipped with the memes package and can be loaded as follows:
```{r, message=F}
data("example_chip_summits", package = "memes")
head(example_chip_summits, 5)
```

The ChIP summits have 3 additional annotation columns: `id`,
`peak_binding_description`, and `e93_sensitive_behavior`. `id` is a unique
identifier for each peak, `peak_binding_description` describes the condition in
which E93 is observed to bind to this region (described in detail below), and
`e93_sensitive_behavior` describes the change in chromatin accessibility
observed at this site following E93 expression. These annotations will be used
throughout this vignette to identify motifs that are associated with differences in E93
binding, and differences in chromatin accessibility at E93 bound sites.

The summit of a ChIP peak represents the single base-pair position with the
highest amplitude ChIP signal within a ChIP peak. Typically, motif analysis
gives the best results when performed in a small window (100 - 200 bp) around
the summit. We can resize our summits to use a 100bp window using `plyranges`:

```{r}
# These data use the dm3 reference genome
dm.genome <- BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3

# Get sequences in a 100bp window around the peak summit
summit_flank <- example_chip_summits %>%
  # this ensures we take 50bp up & downstream of the summit for a total width of 100bp
  plyranges::anchor_center() %>% 
  plyranges::mutate(width = 100)
```

## Determinants of ectopic and orphan binding

One approach to identify enriched motifs is to search for enrichment of known
transcription factor (TF) motifs within sequences of interest. The Fly Factor
Survey is a database of known *Drosophila* TF motifs that we will use throught
this vignette. These data are also packaged with `memes`:

```{r}
meme_db_path <- system.file("extdata/flyFactorSurvey_cleaned.meme", package = "memes", mustWork = TRUE)
```

### Pre-filtering database for expressed transcription factors

A key consideration when using a database of known transcription factors is whether to include motifs for TFs which are not expressed. On one hand, including all motifs allows maximum discovery of potential regulators, but on the other hand may detect factors which are not relevant under your experimental conditions. Another consideration for many tools is that increasing the number of tested motifs increases the multiple testing penaly, so limiting the motif database to a smaller number of relevant candidates can be a way to improve statistical power. There is no one right answer to this question. For this vignette, I will give an example for how one could go about filtering a database to consider only expressed factors. For other analyses, this question will need to be carefully considered.

First, we load the motif database into R using `read_meme`. To more easily manipulate the motif metadata, we can convert this to a data.frame using `to_df`.

```{r}
library(universalmotif)
meme_db <- read_meme(meme_db_path) %>% 
  to_df()
```

To simplify the vignette, we've preproccessed RNAseq data from our tissue of interest and shipped them with `memes`. These are FPKM normalized counts for each transcription factor gene at early and late developmental times. 
```{r}
data("example_rnaseq", package = "memes")
head(example_rnaseq, 4)
```

```{r, fig.height=5, fig.width=5, eval=F, include=F, echo=F}
library(ggplot2)

example_rnaseq %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::filter(fpkm == max(fpkm)) %>% 
  ggplot(aes(fpkm)) +
    geom_histogram(binwidth = 1) +
    coord_cartesian(xlim = c(0, 25)) +
    theme_bw() +
    labs(title = "Largest FPKM value in RNAseq dataset for each gene",
         y = "Number of genes",
         x = "Largest FPKM Value")
```

Choosing a method to decide whether a gene is expressed or not is a non-trivial topic. **DO NOT BLINDLY DO WHAT I DO HERE**. To keep this vignette on topic, I will use an overly simplistic cutoff of FPKM > 5 to define expressed genes. In this example, I will consider any gene that exceeds the cutoff at either timepoint as expressed.
```{r}
expressed_genes <- example_rnaseq %>% 
  # For each gene, keep only those with max(FPKM) greater or equal to 5.
  dplyr::group_by(symbol) %>% 
  dplyr::filter(max(fpkm) >= 5)
```

Finally, I will filter the full motif database to select only those motifs corresponding to expressed gene.
```{r}
meme_db_expressed <- meme_db %>% 
  # the altname slot of meme_db contains the gene symbol
  # (this is database-specific)
  dplyr::filter(altname %in% expressed_genes$symbol)
```

Below is a comparison of the number of motifs in the database before and after filtering.
```{r}
#Number of motifs pre-filtering: 
nrow(meme_db)
#Number of motifs post-filter: 
nrow(meme_db_expressed)
```

Using this filtering strategy, we've decreased our multiple-testing penalty by ~54%.

As a convenience function, we can set a default motif database to use for each `memes` function by setting the R option `meme_db` to a `universalmotif` object. This can be done as follows:

```{r}
# to_list() converts the database back from data.frame format to a standard `universalmotif` object.
options(meme_db = to_list(meme_db_expressed, extrainfo = FALSE))
```

Alternatively, you can pass a universalmotif object to the `database` parameter for any `memes` function and this will override the global option.

### Examination of binding categories with AME

Next, to test for enrichment of known motifs between E93 binding categories, we
first collect the sequences of each category like so:

```{r}
by_binding <- summit_flank %>%
  # Get a list of chip peaks belonging to each set
  split(mcols(.)$peak_binding_description) %>%
  # look up the DNA sequence of each peak within each group
  get_sequence(dm.genome)
```

The data are returned as a Biostrings List, where each list entry represents the sequences of each E93 binding category.

```{r}
head(by_binding)
```

Finally, we test each set of sequencing using AME with the `runAme` function. This will run AME on each set of input sequences in the Biostrings List.

```{r}
ame_by_binding <- by_binding %>% 
  runAme
```

This returns a list object where each entry is the AME results for each
category. For example, here is what the results for the ectopic peaks look like:

```{r}
head(ame_by_binding$ectopic, 5)
```

### Visualizing AME results 

The `plot_ame_heatmap()` function provides a quick way to visualize AME results.
It is built on top of `ggplot2`, so all `ggplot2` functions can be used to
further modify the plot.

By default, it uses the -log10(adjusted p-value) as the heat values. See the
documentation (`?plot_ame_heatmap`) for additional notes on customization.

We can visualize the top 10 hits from the ectopic binding results as follows:

```{r, fig.height=3.5, fig.width=7}
ame_by_binding$ectopic %>% 
  dplyr::filter(rank %in% 1:10) %>% 
  plot_ame_heatmap(group_name = "Ectopic Sites")  +
    ggtitle("Top 10 AME Hits in Ectopic Sites")
```

To plot results from multiple runs together, they must first be joined into a single 
data frame. The `ame_by_binding` object is a list whose names correspond to the
E93 binding category. The list can be combined into a data.frame using
`dplyr::bind_rows`. Setting `.id = "binding_type"` creates a new column
`binding_type` that contains the names from the `ame_by_binding` list. In this
way, the `ame_res` data.frame contains all AME results for each run, which can
be distinguished by the `binding_type` column.

```{r}
ame_res <- ame_by_binding %>% 
  dplyr::bind_rows(.id = "binding_type")
  
```

It is possible to aggregate results from multiple runs into a heatmap by setting
the `group` parameter in `plot_ame_heatmap()`. This will stratify the y-axis by
the entries in the column passed to `group`. This results in clustering motifs
with shared hits across categories along the x-axis for quickly visualizing
motifs that may be shared or distinct between multiple groups.

```{r, fig.height=5,fig.width=15}
ame_res %>% 
  plot_ame_heatmap(group = binding_type)
```

### Reducing redundant motif hits

Another key consideration for the above visualization is that
in the FlyFactorSurvey database we used, different TFs can have multiple motif
entries in the database which are all detected separately by AME. Here, when
returning the top 5 hits from each group, you can see, for example, that a motif
matching "Aef1" is reported 2 times within the top 5 hits of the ectopic sites.
In this situation, it makes sense to summarize the data at the TF level, instead
of the motif level. **Note:** There may be exceptions to this if, for example, a
TF has multiple DNA binding sequences it can recognize, in which case having
multiple hits may reflect a biological property of your sequences. You will have
to handle this on a case-by-case basis for interesting hits and different motif
databases. Here, we can see that at least for Aef1, the consensus
sequences are very similar.

```{r}
ame_res %>% 
  dplyr::filter(binding_type == "ectopic", rank %in% 1:5) %>% 
  head(5) %>% 
  dplyr::select(binding_type, rank, motif_id, motif_alt_id, consensus)
```

How to solve this problem will vary with different motif databases (For details
on how to pre-process a motif database, see `vignette("tidy_motifs")`). In the
tidy version of the FlyFactorSurvey database, the `altname`s of the motifs are
set to the transcription factor gene symbol. This information is included in the
AME results as the `motif_alt_id` column. In order to reduce the occurrance of
redundant hits in the heatmap, we select the hit for each TF with the lowest
p-value (i.e. the most significant hit) as follows:

```{r, fig.height=3, fig.width=12}
ame_res %>%
  # perform the next dplyr operation on each TF within each binding type
  dplyr::group_by(binding_type, motif_alt_id) %>%
  # within each binding type, select the TF hit with the lowest adjusted p-value
  dplyr::filter(adj.pvalue == min(adj.pvalue)) %>% 
  plot_ame_heatmap(group = binding_type, id = motif_alt_id) +
    labs(y = "Binding Category",
         x = "Transcription Factor Motif")
```


### AME Heatmap Visualization

If you have read `vignette("core_ame")`, the "normalized rank" heatmap was
introduced as an alternative visualization method useful when AME produces very
large numbers of hits, or when p-values are on very different scales between groups. 
The `ame_compare_heatmap_methods()` function can be used to visualize why.

Below is a comparison of the distribution of values when using
`-log10(adj.pvalue)` (A) vs normalized ranks (B). Because there are relatively
few hits in the results (~30), and the number of hits between groups varies more
than the `-log10(p-value)` distributions, the "normalize" method will produce a
misleading heatmap relative to the `-log10(p-value)` map.

```{r, fig.height=3, fig.width=7.5}
ame_res %>% 
  dplyr::group_by(binding_type, motif_alt_id) %>% 
  dplyr::filter(adj.pvalue == min(adj.pvalue)) %>% 
  ame_compare_heatmap_methods(group = binding_type)
```

A third option for adjusting this visualization is to cap the heatmap scale
values at a certain value. A good rule of thumb for selecting this cutoff is to
view the cumulative distribution plot above in (A), and select a value on the
x-axis which captures a majority of the data. Here, we see that
`-log10(adj.pvalue) < 10` captures >90% of the entopic & ectopic hits, but only
~50% of the orphan hits. By selecting 10 as the heatmap scale cap, we can
improve the dynamic range of the signal for values <10, while preserving the
information that many orphan sites generally have higher scores than the other
categories. The heatmap scale can be capped at 10 by setting `scale_max = 10`. Below is a
comparison of the 3 different visualizations. Notice how the intensity of the
orphan vs other categories hits changes. 


```{r, fig.height=3, fig.width=12}
best_ame_hits <- ame_res %>% 
  dplyr::group_by(binding_type, motif_alt_id) %>% 
  dplyr::filter(adj.pvalue == min(adj.pvalue))
  
pval_heatmap <- best_ame_hits %>% 
  plot_ame_heatmap(group = binding_type, id = motif_alt_id) +
    labs(x = NULL,
         title = "-log10(adj.pvalue)")
  

norm_heatmap <- best_ame_hits %>% 
  plot_ame_heatmap(group = binding_type, id = motif_alt_id, value = "normalize") +
    labs(x = NULL,
         title = "normalize")

pval_scaled_heatmap <- best_ame_hits %>% 
  plot_ame_heatmap(group = binding_type, id = motif_alt_id, scale_max = 10) +
    labs(x = NULL,
         title = "-log10(adj.pvalue) (scale capped at 10)")
```

```{r, fig.height=9, fig.width=12}
cowplot::plot_grid(pval_heatmap,
                   norm_heatmap,
                   pval_scaled_heatmap,
                   ncol = 1,
                   labels = "AUTO")
```

## *De-novo* motif similarity by binding

Discovery of *de-novo* motifs within binding categories can be another way to
identify meaningful motifs that distinguish between ectopic, entopic, and orphan
sites which does not rely on known motif information. The MEME tool, Dreme,
discovers short *de-novo* motifs in input sequences. We begin the analysis by
searching for 5 *de-novo* motifs per binding category and shuffling the input
sequences to use as the control set.

We pass the same `by_binding` list as the input to `runDreme` which will run
Dreme on each set of sequences independently. The results are returned in list
format, which we combine into a data.frame using `dplyr::bind_rows` as above.

```{r dreme_by_binding, eval=F}
dreme_by_binding <- by_binding %>% 
  runDreme("shuffle", nmotifs = 5) %>% 
  dplyr::bind_rows(.id = "binding_type")
```

```{r}
# The above code chunk takes a long time to run.
# memes is packaged with the results of this run in the "example_dreme_by_binding" dataset
# which can be loaded as follows:
data("example_dreme_by_binding", package = "memes")
dreme_by_binding <- example_dreme_by_binding %>% 
  dplyr::bind_rows(.id = "binding_type")
```

Next, we want to compare the de-novo motifs between each binding category to
identify motifs which could distinguish the three groups. One way to interrogate
this is to examine the correlation score between each motif.

Rename the motifs to indicate the binding category they were discovered in (for visualization purposes).
```{r}
dreme_by_binding_renamed <- dreme_by_binding %>% 
  dplyr::mutate(name = paste(binding_type, seq, sep = "_")) %>% 
  # update_motifs updates the information in the special `motif` column
  update_motifs()
```

Create the correlation heatmap:
```{r, fig.height=5, fig.width=8}
# Set the color values for the heatmap scale
cols <- colorRampPalette(c("white", "dodgerblue4"))(255)

# This is for adding the colored annotation blocks indicating group membership
# to the heatmap
anno.df <- dreme_by_binding_renamed %>% 
  dplyr::select(name, binding_type) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("name")

dreme_by_binding_renamed %>%
  # Convert to universalmotif format 
  to_list() %>%
  # Compute the pearson correlation for each motif with all other motifs
  universalmotif::compare_motifs(method = "PCC") %>%
  # Plot the correlation matrix along with the annotations
  pheatmap::pheatmap(color = cols,
                     # This sets the heatmap range to be from 0-1
                     breaks = seq(0, 1, by = 1/255),
                     annotation_col = anno.df,
                     # the cutree options are just cosmetic to add some spacing
                     # between some of the clusters
                     cutree_rows = 6,
                     cutree_cols = 6,
                     show_colnames = FALSE) 
```

This heatmap reveals that most motifs discovered in each group are highly
similar to the motifs found in other groups. 


### Test *de-novo* motif enrichment using AME

The above analysis suggests the motif content of the different binding
categories are highly similar in sequence composition. To extend these analyses,
we can use AME to test for motif enrichment of the *de-novo* discovered motifs
within each binding category, and determine whether the motifs detected in one
category are indeed enriched in another. To do this, we can provide the
*de-novo* motifs as the AME database to test for their enrichment in each sequence category.

`runAme()` allows using a `runDreme()` results object as the `database` input by
passing it within a `list()`. Naming the `list()` entry produces an informative
`motif_db` name in the results data.frame.

```{r ame_by_binding}
ame_denovo_by_binding <- by_binding %>% 
  runAme(database = list("denovo_binding_motifs" = dreme_by_binding_renamed)) %>% 
  dplyr::bind_rows(.id = "binding_type") 
```

Plotting the heatmap of results reveals that indeed a majority of the *de-novo*
motifs discovered within a single category are detected in all 3 categories,
supporting the conclusion that orphan, ectopic, and entopic sites are highly
similar in sequence content.
```{r, fig.height=4, fig.width=10}
ame_denovo_by_binding %>% 
  plot_ame_heatmap(group = binding_type, scale_max = 10)
```

However, there are 2 interesting motifs which distinguish orphan and ectopic
sites from entopic sites. To help identify which TFs these motifs might belong
to, we can use TomTom to match them to known TF motifs.

First, we select the motif id's that are not found in entopic sites.
```{r}
entopic_motifs <- ame_denovo_by_binding %>% 
  dplyr::filter(binding_type == "entopic") %>% 
  dplyr::pull(motif_id)

ame_denovo_binding_unique <- ame_denovo_by_binding %>% 
  dplyr::filter(!(motif_id %in% entopic_motifs))
```

Next, we use the motif id's from the unique AME results to select those entries in the Dreme
results object, and run TomTom on that subset.
```{r}
dreme_by_binding_unique <- dreme_by_binding_renamed %>% 
  dplyr::filter(name %in% ame_denovo_binding_unique$motif_id) %>% 
  runTomTom(dist = "ed")
```

Finally, we visualize the TomTom results to identify candidate TFs driving the presence of ectopic and orphan sites.
```{r, fig.height=4, fig.width=12}
dreme_by_binding_unique %>% 
  view_tomtom_hits(3) %>%
  cowplot::plot_grid(plotlist = ., nrow = 1, labels = "AUTO", byrow = TRUE)
```

## Motifs in opening vs closing sites

Another category of ChIP peaks in this dataset are ones that overlap regions
with dynamic changes to chromatin accessibility in response to E93 expression.
These peaks are annotated in the `e93_sensitive_behavior` column which indicates
whether sites have higher accessibility after ectopic E93 expression ("Increasing"),
lower accessibility after ectopic E93 expression ("Decreasing"), or if the accessibility
does not change ("Static").

To identify motifs associated with Increasing or Decreasing accessibility in
response to ectopic E93 expression, we first grab the DNA sequences of the Increasing,
Decreasing, and Static sites by splitting the ChIP peaks on the
`e93_sensitive_behavior` column.

```{r}
# split by response to E93 binding
by_sens <- summit_flank %>% 
  split(mcols(.)$e93_sensitive_behavior) %>% 
  get_sequence(dm.genome)
```

Proper selection of the control set of sequences is critical to identifying
biologically relevant motifs. Because we are interested in motifs associated
with Increasing or Decreasing behavior, we want to exclude detection of motifs
that are simply associated with E93 binding. To do this, we can use the "Static"
category of E93 ChIP peaks, because these theoretically represent E93 binding
sites without motifs that influence chromatin accessibility. To start, we will
search for *de-novo* motifs within dynamic E93 binding sites using Dreme.

The `by_sens` list has three entries: `Increasing`, `Decreasing`, and `Static`.
When calling `runDreme` on a sequence list, you can pass the name of a list
entry as the `control` argument to use that set as the background sequences for
the remaining sets of sequences. Below we search for *de-novo* motifs in
Increasing and Decreasing peaks using the Static peaks as the background set. 


```{r dreme_by_sens_vs_static}
# By setting "Static" as the control, runDreme will run:
# - Increasing vs Static
# - Decreasing vs Static

dreme_by_sens_vs_static <- runDreme(by_sens, "Static")
```

The results are returned as a list with two entries: `Increasing`, and
`Decreasing`. These can be joined into a single data.frame using
`dplyr::bind_rows`.

```{r}
dreme_results <- dreme_by_sens_vs_static %>% 
  dplyr::bind_rows(.id = "e93_response")
```

Next, we match the *de-novo* motifs with known TF motifs using TomTom. The TomTom
results will be added as new columns to the input data. TomTom data columns
begin with `best_match_` summarizing the top hit, and the `tomtom` column is a
nested data.frame storing the full results for each input motif.

```{r}
denovo_dynamic_motifs <- dreme_results %>% 
  runTomTom()
```

We can visualize some of the top TomTom hits for each *de-novo* motif using
`view_tomtom_hits()`. Passing a number to `view_tomtom_hits()` shows that many
matches ranked in descending order. Passing no value will show all hits. It's a
good idea to inspect all hits visually. This can be an extremely revealing step,
as it will give you insight into which motifs are generally similar to your
input motif. You may also discover proteins with identical motifs, in which case
you may need to decide which TF is the better candidate. For brevity, I'll only
show the top 3 hits for each motif.

```{r, fig.height=8, fig.width=8.5}
denovo_dynamic_motifs %>% 
  # Plot the top 3 matches
  view_tomtom_hits(3) %>% 
  # Arrange into figure panels
  # You don't need to do this if you're just exploring your data, I do it here to make a pretty figure
  cowplot::plot_grid(plotlist = ., labels = "AUTO")
```

After manually inspecting the above hits, we may decide that a lower-ranked
motif is a better match by eye, or because of domain-specific knowledge that
makes a different TF a better candidate. We can swap the best matched motif
using `force_best_match()`. This will update the `best_match_` columns and
re-rank the `tomtom` nested data column to set a lower ranked match to the top
hit. The updates `best_match_` columns will reflect the values corresponding to
the new match.

As an example, I'll set the top hit for `m03_AKGG` to the lower ranked
`pho_SANGER_10` since it has weaker consensus sequence in the regions that don't
overlap the *de-novo* motif.

```{r}
# use force_best_match to use update the best match info
# %<>% is a shortcut from magrittr for updating the data in-place
denovo_dynamic_motifs %<>% 
  force_best_match(c("m03_AKGG" = "pho_SANGER_10"))
```

Replotting the top hits, you can see that "pho_SANGER_10" is now listed as the top hit for "m03_AKGG".

```{r, fig.height=4.5, fig.width=8.5}
denovo_dynamic_motifs %>% 
  view_tomtom_hits(1) %>% 
  cowplot::plot_grid(plotlist = ., labels = "AUTO")
```

Finally, we can make publication quality figures using cowplot. This code is a
little complicated, but hopefully gives an idea for how you can use R graphics
to generate nice motif figures!

A striking result from this is that the *de-novo* motifs associated with
Decreasing sites match the E93 motif itself.

```{r, fig.height=3, fig.width=8.5}
denovo_dynamic_motifs %>%
  # Create a new label that indicates whether the motif is found in Increasing or Decreasing sites
  dplyr::mutate(label = paste0(e93_response, " in response to E93")) %>%
  # Sort the top hits first
  dplyr::arrange(rank) %>%
  # Run the plot command separately for Increasing and Decreasing sets
  split(.$label) %>% 
  # For each set:
  purrr::imap(~{
    # Plot the top hit for each motif in the set as a panel figure
    top_hits <- view_tomtom_hits(.x, 1) %>% 
      cowplot::plot_grid(plotlist = ., nrow = 1, labels = "AUTO")

    # Create a figure title from the label we created above
    title <- cowplot::ggdraw() +
        cowplot::draw_text(.y)

    # Combine the title & the motif hits figure
    cowplot::plot_grid(plotlist = list(title, top_hits), 
                       ncol = 1,
                       rel_heights = c(0.1, 1)
    )
  })
```


## Scanning for motif matches using FIMO

Because we discover a *de-novo* motif from DREME that matches to E93 in
Decreasing sites, this suggests that the E93 motifs in Decreasing sites have
different properties relative to E93 motifs in Increasing or Static sites. For
example, the E93 motif may have higher similarity to the canonical sequence, may
be present in greater number, or may have different positioning in Decreasing
sites relative to Increasing or Static sites. To test these questions directly,
we can identify individual E93 motif matches within each sequence and
investigate their properties.  The FIMO tool is used to identify the positions
of motif matches in a set of sequences. It will return the coordinates,
sequence, and a score (higher values are more similar to the consensus) for each
motif found in the input sequences.

In order to scan for the E93 motif within target sequences, we first need a copy
of the E93 motif. Although we could import the motifs from our local meme
database, it is also possible to use motifs pulled from a
[`MotifDb`](https://bioconductor.org/packages/release/bioc/html/MotifDb.html)
query as follows. I'll use that approach here as an example, but if this were a
real analysis, I'd take the result from our local database.

```{r}
e93_flyfactor <- MotifDb::MotifDb %>% 
  # Query the database for the E93 motif using it's gene name
  MotifDb::query("Eip93F") %>% 
  # Convert from motifdb format to universalmotif format
  universalmotif::convert_motifs() %>% 
  # The result is a list, to simplify the object, return it as a single universalmotif
  .[[1]]

# Rename the motif from it's flybase gene number to a more user-friendly name
e93_flyfactor["name"] <- "E93_FlyFactor"

# Display the cleaned motif
e93_flyfactor
```

Just like `runDreme()` and `runAme()`, `runFimo()` takes a
`Biostrings::XStringSet` as input. Using `get_sequence()` to generate these for
DNA sequences is preferred because it creates sequence inputs that can be parsed
into genomic coordinates by Fimo.

```{r}
fimo_res <- summit_flank %>% 
  get_sequence(dm.genome) %>% 
  runFimo(motifs = e93_flyfactor, thresh = 1e-3)
```

The results of `runFimo()` are returned as a `GRanges` object containing the
positions of each motif discovered in the input sequences. The best way to
integrate these data with our input peaks is to use the `plyranges` suite of
tools for performing overlaps and joins between `GRanges` objects.

```{r}
fimo_res
```

### Counting the number of motifs per peak

`plyranges` extends `dplyr`-like syntax to range objects. Here we add a count
for each motif per peak in the `n_motifs` column. We also add a column
`has_motif` which will be a binary indicator of whether a peak contains any
motifs.

```{r}
summit_flank %<>% 
  # Note, if running FIMO with multiple motifs, this solution will not work
  # as it will count all motifs within the fimo-results without splitting by motif_id
  plyranges::mutate(n_motifs = plyranges::count_overlaps(., fimo_res), 
                    has_motif = n_motifs > 0)
```

First, we want to determine whether E93 sensitive sites are more likely to
have E93 motifs in certain response types. Here we can see that sensitive
decreasing sites are more likely to have E93 motifs than sensitive increasing or
insensitive static sites.

```{r, fig.height=4, fig.width=5}
summit_flank %>% 
  data.frame %>% 
  dplyr::mutate(has_match = dplyr::if_else(has_motif, "Match", "No Match")) %>% 
  ggplot(aes(e93_sensitive_behavior)) +
    geom_bar(aes(fill = forcats::fct_rev(has_match)), position = "fill") +
    scale_fill_manual(values = c("Match" = "firebrick",
                                 "No Match" = "Black")) +
    labs(fill = "E93 Motif Match",
         y = "Fraction of Sites",
         x = "Response to E93 binding")
```

To investigate whether E93-sensitive sites that increase, decrease, or remain static have
different numbers of E93 motifs, we plot the fraction of sites with each number
of motifs per group. Here we can see that in addition to being more likely to contain an
E93 motif, sensitive decreasing sites are more likely to to contain 2 or more
matches, where 10% contain at least 2 motifs.

```{r, fig.height=5, fig.width=7}
summit_flank %>% 
  # currently, group operations are faster as data.frames, so we convert to data.frame
  data.frame %>%
  dplyr::group_by(e93_sensitive_behavior, n_motifs) %>% 
  dplyr::count() %>% 
  dplyr::group_by(e93_sensitive_behavior) %>% 
  dplyr::mutate(frac = n/sum(n)) %>% 
  ggplot(aes(n_motifs, frac)) +
    geom_line(aes(color = e93_sensitive_behavior), size = 1) +
    labs(y = "Fraction of Sites",
         x = "Number of E93 Motifs",
         color = "Response to E93 Binding") +
    theme_bw()
```

Finally, we want to assess whether the quality of E93 motifs is different
between sensitivity categories. To examine this, we need to determine which
motifs are found in which peaks. We use `plyranges::join_overlap_intersect`
to return motif entries appended with peak-level metadata, like the peak id each
motif is found within.

```{r}
# return position of each motif match w/ peak metadata
intersect <- fimo_res %>% 
  plyranges::join_overlap_intersect(summit_flank)
```

We use the FIMO `score` as a proxy for quality, where higher scores are better
matches to the motif. Here we examine only the best match (highest score) motif
within each peak.

```{r}
best_motifs <- intersect %>%
  # group by ChIP peak id
  plyranges::group_by(id) %>% 
  # Keep only motifs w/ the highest score within a peak
  plyranges::filter(score == max(score)) %>% 
  plyranges::ungroup()
```

Finally, visualize the FIMO scores for the best E93 motif within each ChIP peak
across categories. This reveals that Decreasing sites appear to have higher FIMO
scores, suggesting the E93 motifs within Decreasing sites are higher quality E93
motifs than those in Increasing or Static sites.

```{r, fig.height=4, fig.width=3}
best_motifs %>% 
  data.frame %>% 
  ggplot(aes(e93_sensitive_behavior, score)) +
    geom_boxplot(aes(fill = e93_sensitive_behavior), notch = TRUE, width = 0.5) +
    guides(fill = "none") +
    labs(x = "Response to E93",
         y = "FIMO Score") +
    theme_bw()
```

We can test whether the distribution of FIMO scores between multiple groups is
statistically significantly different using
an [Anderson-Darling test](https://en.wikipedia.org/wiki/Anderson%E2%80%93Darling_test). For a
pairwise comparison (in case you only have 2 groups to test), a [Kolmogorov-Smirnov test](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test) is commonly
used instead. The `PMCMRplus::adAllPairsTest()` performs the Anderson-Darling
test and corrects for multiple testing. This reveals that Decreasing FIMO scores
are significantly different from Increasing and Static scores (p < 0.05), but Increasing
and Static scores are not significantly different from eachother (p > 0.05), supporting the
conclusion that E93 motifs in Decreasing sites are higher quality than those in
Increasing or Static sites.

```{r}
best_motifs %>% 
  data.frame %>% 
  dplyr::mutate(behavior = factor(e93_sensitive_behavior)) %>% 
  PMCMRplus::adAllPairsTest(score ~ behavior, ., p.adjust.method = "fdr")
```

We can visualize the differences in E93 motif quality between groups by
reconstructing PCM's from the matched motif sequences. Here, it is clear that
although each category contains an E93 motif match, the sequence content of
those motifs differs between categories, providing a nice visual representation
of the statistical analysis above.

```{r, fig.height=6, fig.width=5}
best_motifs %>% 
  # Get sequences of each matched motif
  add_sequence(dm.genome) %>% 
  # Split by E93 response category
  split(mcols(.)$e93_sensitive_behavior) %>% 
  # Convert from GRangesList to list() to use `purrr::imap`
  as.list() %>% 
  # imap passes the list entry as .x and the name of that object to .y
  # here, .y holds the e93_sensitive_behavior name, and .x holds the motif positions
  purrr::imap(~{
    # This builds a new motif using the matched sequences from each group
    create_motif(.x$sequence, 
                 # Append the response category to the motif name
                 name = paste0("E93_", .y))
    }) %>% 
  # Prepend E93 Fly Factor motif to beginning of motif list so it is plotted at the top
  c(list("E93_FlyFactor" = e93_flyfactor), .) %>% 
  # Plot each motif
  view_motifs()
```

We can also build a more complex representation of these data by making a
heatmap of the sequence variability for each matched sequence using
`plot_sequence_heatmap()`. This function will plot each sequence as a row of a
heatmap where each column represents the position relative to the motif start
coordinate. The cells are colored by the DNA sequence found at that position for
each sequence. Using this type of plot it becomes even more clear how
variability at position 6 differs between Decreasing and Increasing sites, where
Decreasing sites have a very consistent "C" at position 6, while Increasing
sites are more likely to have non-C bases at this position.

```{r, fig.height=5.25, fig.width=8}
best_motifs %>% 
  # Get sequences of each matched motif
  add_sequence(dm.genome) %>% 
  data.frame %>% 
  # Split by E93 response category
  split(.$e93_sensitive_behavior) %>% 
  # Extract the "sequence" column from each category
  purrr::map("sequence") %>% 
  # Plot a heatmap for each category
  plot_sequence_heatmap(title_hjust = 0.5) %>% 
  # Arrange into a grid
  cowplot::plot_grid(plotlist = ., nrow = 1)
```


### Centrality of E93 motif

Next we want to visualize whether the E93 motif has altered positioning across
the response categories. To do this we add the metadata for each nearest motif
to our peak summits. Setting `distance = TRUE` in `plyranges::join_nearest` adds
a `distance` column indicating the distance to the nearest joined region.

```{r}
summit_nearest_e93_motif <- summit_flank %>% 
  plyranges::anchor_center() %>% 
  plyranges::mutate(width = 1) %>% 
  plyranges::join_nearest(fimo_res, distance = TRUE)
```

We can visualize these data as a cumuluative distribution plot. In this example,
it does appear that decreasing sites tend to be closer to the peak summit than
increasing or static sites.

```{r}
summit_nearest_e93_motif %>% 
  data.frame %>%
  # Limit our analysis to peaks that contain motifs
  dplyr::filter(has_motif == TRUE) %>% 
  ggplot(aes(distance)) +
    stat_ecdf(aes(color = e93_sensitive_behavior), size = 1, pad = FALSE) +
    labs(x = "Distance to nearest E93 Motif",
         y = "Fraction of Sites",
         color = "E93 Response") +
    theme_linedraw() 
```

```{r, include=F, eval=F, echo=F}
summit_nearest_e93_motif %>% 
  data.frame %>% 
  dplyr::filter(has_motif == TRUE) %>% 
  dplyr::mutate(behavior = factor(e93_sensitive_behavior)) %>% 
  # Anderson-darling test for multiple distributions
  PMCMRplus::adAllPairsTest(distance ~ behavior, ., p.adjust.method = "fdr") 

```

# Conclusion

This is just one example of how `memes` can integrate multiple data sources to
perform deep analysis of genomics data. Here we went from a set of annotated
ChIP peaks do multifaceted interrogation of motif content. Each approach has
advantages and disadvantages. The most important thing is to critically evaluate
the data at every step: visualize *de-novo* motif hits with `view_motifs`,
compare matches using `view_tomtom_hits`, and think carefully about which
sequences to use as control regions.

# References
Nystrom SL, Niederhuber MJ, McKay DJ. Expression of E93 provides an instructive cue to control dynamic enhancer activity and chromatin accessibility during development. Development 2020 Mar 16;147(6).

# Session Info
```{r}
sessionInfo()
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example_dreme_by_binding}
\alias{example_dreme_by_binding}
\title{runDreme() output for example_chip_summits split by binding description}
\format{
a runDreme results data.frame
}
\usage{
example_dreme_by_binding
}
\description{
See vignette("integrative_analysis", package = "memes") for details
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{add_sequence}
\alias{add_sequence}
\title{Add nucleic acid sequence of regions to metadata column}
\usage{
add_sequence(ranges, genome, name = "sequence")
}
\arguments{
\item{ranges}{GRanges object}

\item{genome}{BSgenome object or any other valid input to
`Biostrings::getSeq()` (Do `showMethods(Biostrings::getSeq)` to show valid
types)}

\item{name}{name of metadata column to hold sequence information (default:
"sequence"). Note, this will overwrite existing columns without warning if
the name already exists.}
}
\value{
`ranges` with new metadata column named "sequence" (or another value
  passed to `name`) holding the DNA or RNA sequence from `genome`
}
\description{
Add nucleic acid sequence of regions to metadata column
}
\examples{
data(example_peaks, package = "memes")
dm.genome <- BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3
add_sequence(example_peaks, dm.genome)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_duplicate_motifs.R
\name{has_duplicate_motifs}
\alias{has_duplicate_motifs}
\title{Check for duplicated motif matrices}
\usage{
has_duplicate_motifs(x)
}
\arguments{
\item{x}{a universalmotif list or universalmotif_df}
}
\value{
logical value indicating presence or absence of duplicated motif matrices
}
\description{
This function identifies whether any motif matrices in the input
universalmotif list or universalmotif_df are identical to each other. Note:
this operation is slow on large motif lists
}
\examples{
motif <- universalmotif::create_motif()
duplicated <- c(motif, motif)
has_duplicate_motifs(duplicated)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fimo.R
\name{importFimo}
\alias{importFimo}
\title{Import FIMO results}
\usage{
importFimo(fimo_tsv)
}
\arguments{
\item{fimo_tsv}{path to fimo.tsv output file}
}
\value{
GenomicRanges object for each match position. Note unless coordinates
  are genomic positions, each `seqnames` entry will be the fasta header, and
  start/end will be the position within that sequence of the match.
}
\description{
Import FIMO results
}
\examples{
fimo_tsv <- system.file("extdata", "fimo.tsv", package = "memes")
importFimo(fimo_tsv)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{write_fasta}
\alias{write_fasta}
\title{Write fasta file from stringset}
\usage{
write_fasta(seq, path = tempfile(fileext = ".fa"))
}
\arguments{
\item{seq}{a `Biostrings::XStringSet`}

\item{path}{path of fasta file to write (default: temporary file)}
}
\value{
path to created fasta file
}
\description{
Write fasta file from stringset
}
\examples{
seq <- universalmotif::create_sequences()
\donttest{
write_fasta(seq)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tomtom.R
\name{runTomTom}
\alias{runTomTom}
\title{Run TomTom on target motifs}
\usage{
runTomTom(
  input,
  database = NULL,
  outdir = "auto",
  thresh = 10,
  min_overlap = 5,
  dist = "ed",
  evalue = TRUE,
  silent = TRUE,
  meme_path = NULL,
  ...
)
}
\arguments{
\item{input}{path to .meme format file of motifs, a list of universalmotifs,
or a universalmotif data.frame object (such as the output of \code{runDreme()})}

\item{database}{path to .meme format file to use as reference database (or
list of universalmotifs). \strong{NOTE:} p-value estimates are inaccurate when
the database has fewer than 50 entries.}

\item{outdir}{directory to store tomtom results (will be overwritten if
exists). Default: location of input fasta file, or temporary location if using universalmotif input.}

\item{thresh}{report matches less than or equal to this value. If evalue =
TRUE (default), set an e-value threshold (default = 10). If evalue = FALSE,
set a value between 0-1 (default = 0.5).}

\item{min_overlap}{only report matches that overlap by this value or more,
unless input motif is shorter, in which case the shorter length is used as
the minimum value}

\item{dist}{distance metric. Valid arguments: \code{allr | ed | kullback | pearson | sandelin | blic1 | blic5 | llr1 | llr5}.
Default: \code{ed} (euclidean distance).}

\item{evalue}{whether to use E-value as significance threshold (default:
\code{TRUE}). If evalue = FALSE, uses \emph{q-value} instead.}

\item{silent}{suppress printing stderr to console (default: TRUE).}

\item{meme_path}{path to "meme/bin/" (optional). If unset, will check R
environment variable "MEME_DB (set in \code{.Renviron}), or option
"meme_db" (set with \code{option(meme_db = "path/to/meme/bin")})}

\item{...}{additional flags passed to tomtom using {cmdfun} formating (see table below for details)}
}
\value{
data.frame of match results. Contains \code{best_match_motif} column of
\code{universalmotif} objects with the matched PWM from the database, a series
of \verb{best_match_*} columns describing the TomTom results of the match, and a
\code{tomtom} list column storing the ranked list of possible matches to each
motif. If a universalmotif data.frame is used as input, these columns are
appended to the data.frame. If no matches are returned, \code{tomtom} and
\code{best_match_motif} columns will be set to \code{NA} and a message indicating
this will print.
}
\description{
TomTom compares input motifs to a database of known, user-provided motifs to
identify matches.
}
\details{
runTomTom will rank matches by significance and return a
best match motif for each input (whose properties are stored in the \verb{best_match_*}
columns) as well as a ranked list of all possible matches stored in the
\code{tomtom} list column.

Additional arguments

runTomTom() can accept all valid tomtom arguments passed to \code{...} as described in the
\href{http://meme-suite.org/doc/tomtom.html?man_type=web}{tomtom commandline reference}. For
convenience, below is a table of valid arguments, their default values, and
their description.\tabular{cccl}{
   TomTom Flag \tab allowed values \tab default \tab description \cr
   bfile \tab file path \tab \code{NULL} \tab path to background model for converting frequency matrix to log-odds score (not used when \code{dist} is set to "ed", "kullback", "pearson", or "sandelin" \cr
   motif_pseudo \tab \code{numeric} \tab 0.1 \tab pseudocount to add to motifs \cr
   xalph \tab \code{logical} \tab FALSE \tab convert alphabet of target database to alphabet of query database \cr
   norc \tab \code{logical} \tab FALSE \tab Do not score reverse complements of motifs \cr
   incomplete_scores \tab \code{logical} \tab FALSE \tab Compute scores using only aligned columns \cr
   thresh \tab \code{numeric} \tab 0.5 \tab only report matches with significance values <= this value. Unless \code{evalue = TRUE}, this value must be < 1. \cr
   internal \tab \code{logical} \tab FALSE \tab forces the shorter motif to be completely contained in the longer motif \cr
   min_overlap \tab \code{integer} \tab 1 \tab only report matches that overlap by this number of positions or more. If query motif is smaller than this value, its width is used as the min overlap for that query \cr
   time \tab \code{integer} \tab \code{NULL} \tab Maximum runtime in CPU seconds (default: no limit) \cr
}
}
\section{Citation}{
If you use \code{runTomTom()} in your analysis, please cite:

Shobhit Gupta, JA Stamatoyannopolous, Timothy Bailey and William Stafford
Noble, "Quantifying similarity between motifs", Genome Biology, 8(2):R24,
2007. \href{http://genomebiology.com/2007/8/2/R24}{full text}
\subsection{Licensing}{

The MEME Suite is free for non-profit use, but for-profit users should purchase a
license. See the \href{http://meme-suite.org/doc/copyright.html}{MEME Suite Copyright Page} for details.
}
}

\examples{
if (meme_is_installed()) {
motif <- universalmotif::create_motif("CCRAAAW")
database <- system.file("extdata", "flyFactorSurvey_cleaned.meme", package = "memes")

runTomTom(motif, database)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ame_methods.R
\name{plot_ame_heatmap}
\alias{plot_ame_heatmap}
\title{Plot AME heatmap clustered by similarity in detected motifs}
\usage{
plot_ame_heatmap(
  ame,
  id = motif_id,
  group = NULL,
  value = -log10(adj.pvalue),
  group_name = NULL,
  scale_max = NA
)
}
\arguments{
\item{ame}{ame results data.frame}

\item{id}{column of motif ids to use (default: motif_id).}

\item{group}{grouping column if comparing across multiple ame runs (optional,
default: NULL).}

\item{value}{value to display as heatmap intensity. Default:
-log10(adj.pvalue). Takes function or column name as input. If set to
"normalize", will use normalized rank within `group` as the heatmap values.
**If in doubt**, prefer the -log10(adj.pvalue) plot potentially in
conjunction with adjusting `scale_max`. (See "Normalized rank
visualization" section below for more details on how to interpret these
data)}

\item{group_name}{when group = NULL, name to use for input regions. Ignored
if group is set.}

\item{scale_max}{max heatmap value to limit upper-value of scale. Useful if
distribution of `value`s vary greatly between groups. Usually a better to
tweak this option than to use value = "normalize". The cumulative
distribution plot generated by `ame_compare_heatmap_methods()` can be
useful for selecting this value, try to pick a value which captures the
largest fraction of hits across all groups while excluding outliers.}
}
\value{
`ggplot` object
}
\description{
Plot AME heatmap clustered by similarity in detected motifs
}
\details{
Normalized rank visualization
**NOTE:** The normalized rank visualization eliminates all real values
related to statistical significance! Instead, this visualization represents
the relative ranks of hits within an AME run, which already pass a
significance threshold set during `runAME()`. This means that even if several
motifs have similar or even identical pvalues, their heatmap representation
will be a different color value based on their ranked order in the results
list. This also means that using the normalized rank visualization will
be misleading if there are only a few AME hits; it is only worth using if the
number of hits is very large (>100). Both visualizations can be useful and
reveal different properties of the data to the user during data
exploration. Use `ame_compare_heatmap_methods()` to help assess
differences in the two visualizations. **If in doubt**, prefer the
`-log10(adj.pvalue)` representation.

Common mistake: if `value` is set to a string that is not "normalize", it
will return: "Error: Discrete value supplied to continuous scale". To use a
column by name, do not quote the column name.
}
\examples{
data("example_ame", package = "memes")

# Plot a single category heatmap
plot_ame_heatmap(example_ame$Decreasing)

# Plot a multi category heatmap
grouped_ame <- dplyr::bind_rows(example_ame, .id = "category")
plot_ame_heatmap(grouped_ame, group = category)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_xml.R
\name{importDremeXML}
\alias{importDremeXML}
\title{Import Dreme output from previous run}
\usage{
importDremeXML(dreme_xml_path)
}
\arguments{
\item{dreme_xml_path}{path to dreme.xml file}
}
\value{
data.frame with statistics for each discovered motif. The `motifs`
  column contains a universalmotif object representation in PCM format of
  each DREME motif. If no motifs are discovered, returns NULL.
}
\description{
Import Dreme output from previous run
}
\examples{
dreme_xml <- system.file("extdata", "dreme.xml", package = "memes")
importDremeXML(dreme_xml)
}
\seealso{
[runDreme()]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match_utils.R
\name{nest_tomtom}
\alias{nest_tomtom}
\title{Nest TomTom results columns into a data.frame column named "tomtom"}
\usage{
nest_tomtom(data)
}
\arguments{
\item{data}{tomtom results data.frame after unnesting the `tomtom` column}
}
\value{
the input data.frame with the match_* columns nested into a column named `tomtom`
}
\description{
This is a convienience function for re-nesting the `tomtom` list column if
the user unnests it. Additionally, it will update the best_match information
based on the ranking of the resulting `tomtom` data.frame. This avoids having
out-of-date best_match information after manipulating the `tomtom` entries.
}
\details{
**NOTE:** that the resulting columns may not be in the same order, so
operations like `identical()` before & after a nest/renest operation may fail
even though the column values are unchanged.
}
\examples{
if (meme_is_installed()){
motif <- universalmotif::create_motif("CCRAAAW")
db <- system.file("extdata/flyFactorSurvey_cleaned.meme", package = "memes")
res <- runTomTom(motif, database = db)
data <- tidyr::unnest(res, "tomtom")
identical(nest_tomtom(data), res)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{get_sequence}
\alias{get_sequence}
\title{Get sequence from GRanges}
\usage{
get_sequence(regions, genome, score_column, ...)
}
\arguments{
\item{regions}{GRanges, or GRangesList object. Will also accept a data.frame
as long as it can be coerced to a GRanges object, or a string in the
format: "chr:start-end" (NOTE: use 1-based closed
intervals, not BED format 0-based half-open intervals).}

\item{genome}{object of any valid type in `showMethods(Biostrings::getSeq)`.
Commonly a BSgenome object, or fasta file. Used to look up sequences in regions.}

\item{score_column}{optional name of column (in mcols() of `regions`)
containing a fasta score that is added to the fasta header of each entry.
Used when using [runAme()] in partitioning mode. (default: `NULL`)}

\item{...}{additional arguments passed to Biostrings::getSeq.}
}
\value{
`Biostrings::DNAStringSet` object with names corresponding to genomic
  coordinates. If input is a list object, output will be a
  `Biostrings::BStringSetList` with list names corresponding to input list
  names.
}
\description{
A light wrapper around Biostrings::getSeq to return named DNAStringSets, from
input genomic coordinates.
}
\examples{
# using character string as coordinates
# using BSgenome object for genome
drosophila.genome <- BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6
get_sequence("chr2L:100-200", drosophila.genome)

# using GRanges object for coordinates
data(example_peaks, package = "memes")
get_sequence(example_peaks, drosophila.genome)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{view_tomtom_hits}
\alias{view_tomtom_hits}
\title{Compare top tomtom hits to original motif}
\usage{
view_tomtom_hits(results, top_n = "all")
}
\arguments{
\item{results}{results data.frame from [runTomTom()]}

\item{top_n}{number of matched motifs to return in plot (default: "all")}
}
\value{
plot of input motif vs the top n number of tomtom matched motifs. If
  no match found, will plot "No Match". Note: the "No Match" plots are not
  amenable to ggplot theme() manipulations, while all others are.
}
\description{
Although TomTom does a good job of matching unknown motifs to known motifs,
sometimes the top hit is not the correct assignment. It can be useful to
manually inspect the hits. This function provides a quick utility to compare
matches.
}
\details{
This is intended to be a function used interactively and may not always be
the best tool for creating publication-quality figures. Results with matches
return ggseqlogo outputs which can be further manipulated using
[ggplot2::theme()] calls, but results containing no matches are static plots.
}
\examples{
results <- importTomTomXML(system.file("extdata", "tomtom.xml", package = "memes"))
# show top 3 hits
view_tomtom_hits(results, top_n = 3)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match_utils.R
\name{drop_best_match}
\alias{drop_best_match}
\title{Drop best match columns from tomtom results}
\usage{
drop_best_match(res)
}
\arguments{
\item{res}{results of runTomTom}
}
\value{
`res` without the tomtom best_match_ columns
}
\description{
Convenience function for dropping all columns created by runTomTom prefixed
by "best_match_" and the "best_db_name" column. Keeps the "tomtom" data.frame
column. Can be useful if you want to unnest the `tomtom` data without
propagating these columns.
}
\examples{
data("example_dreme_tomtom")
names(example_dreme_tomtom)
names(drop_best_match(example_dreme_tomtom))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_xml.R
\name{importTomTomXML}
\alias{importTomTomXML}
\title{Import tomtom data from previous run}
\usage{
importTomTomXML(tomtom_xml_path)
}
\arguments{
\item{tomtom_xml_path}{path to tomtom.xml}
}
\value{
will return data.frame with input motifs & results for best match.
  `tomtom` list column contains full tomtom data for each input motif.
  NOTE: if tomtom detects no matches for any input motif, currently will
  print a message & return NA values for `tomtom`, `best_match_name`, and
  `best_match_motif`.
}
\description{
Import tomtom data from previous run
}
\details{
tomtom list column format
the `tomtom` list column contains data.frames with the following format:
    - name: name of query PWM
    - altname: alternate name of query PWM
    - match_name: name of matched PWM
    - match_altname: alt name of matched PWM
    - match_pval: p-value of match
    - match_eval: E-value of match
    - match_qval: q-value of match
    - match_offset: number of letters the query was offset from the target match
    - match_strand: whether the motif was found on input strand (+) or as reverse-complement (-)
    - db_name: database source of matched motif
    - match_motif: universalmotif object containing the PWM that was matched
}
\examples{
tomtom_xml <- system.file("extdata", "tomtom.xml", package = "memes")
importTomTomXML(tomtom_xml)
}
\seealso{
[runTomTom()]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example_dreme}
\alias{example_dreme}
\title{Example runDreme() output}
\format{
a runDreme results data.frame
}
\usage{
example_dreme
}
\description{
Result when running dreme using 100bp window around example_chip_summits
using "Decreasing" sites as foreground, and "Static" sites as background.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{check_meme_install}
\alias{check_meme_install}
\title{Check user's MEME install}
\usage{
check_meme_install(meme_path = NULL)
}
\arguments{
\item{meme_path}{path to "meme/bin" (if unset will search \code{MEME_BIN}
environment variable or \code{meme_bin} option)}
}
\value{
message indicating which MEME utilities are installed and their
location on disk
}
\description{
In order to use the run* family of functions, memes must detect a local
install of the MEME Suite. MEME is installed in a directory named meme/bin/
which can be located anywhere on the filesystem, but is typically found in \verb{~/meme/bin}.
If the MEME Suite is installed at \verb{~/meme/bin}, memes can autodetect the install. However,
in the case that the MEME Suite is found at a nonstandard location, the user
may specify the location of their meme/bin in three ways:
}
\details{
\enumerate{
\item provide the full path to \code{meme/bin} to the \code{meme_path} argument to each \verb{run*} function.
\item set the \code{meme_bin} option using \code{options(meme_bin = "path/to/meme/bin")} once per R session.
\item set the \code{MEME_BIN} environment variable either in \code{.Renviron} or \verb{~/.bashrc} with the path to \code{meme/bin}
}

To aid the user in determining if memes can detect their \code{meme/bin} install,
\code{check_meme_install()} will search the aforementioned locations for a valid
\code{meme/bin}, returning green checks for each detected tool, or red X's for
undetected tools. Alternatively, users can run \code{meme_is_installed()} to get a
boolean value indicating whether their MEME Suite can be detected.

\code{check_meme_install()} searches using the following heirarchy. This heirarchy
mimics how all \verb{run*} functions search for \code{meme/bin}, thus the paths printed
from \code{check_meme_install()} will indicate the paths used by each \verb{run*}
function. The heirarchy is as follows:
\enumerate{
\item the \code{meme_path} function argument if set
\item the \code{meme_bin} option
\item the \code{MEME_BIN} environment variable
\item the default location at \verb{~/meme/bin}
}
}
\examples{
check_meme_install()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example_tomtom}
\alias{example_tomtom}
\title{Example runTomTom() output}
\format{
a data.frame
}
\usage{
example_tomtom
}
\description{
Result when running `runTomTom(example_dreme$motif)` using FlyFactorSurvey as database
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example_ame}
\alias{example_ame}
\title{Example runAme() output}
\format{
A list object of AME results data.frames
\describe{
 \item{Increasing}{`runAme()` Results object for Increasing sites vs Static sites}
 \item{Decreasing}{`runAme()` Results object for Decreasing sites vs Static sites}
}
}
\usage{
example_ame
}
\description{
Result when running AME using 100bp window around `example_chip_summits` for
"Increasing" and "Decreasing" sites, using "Static" as background.
}
\examples{
# Data can be combined into 1 large data.frame using:
# where the "behavior" column will hold the "Increasing"/"Decreasing" information
dplyr::bind_rows(example_ame, .id = "behavior")
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meme.R
\name{runMeme}
\alias{runMeme}
\alias{runMeme.list}
\alias{runMeme.BStringSetList}
\alias{runMeme.default}
\title{Identify motifs with MEME}
\usage{
runMeme(
  input,
  control = NA,
  outdir = "auto",
  alph = "dna",
  parse_genomic_coord = TRUE,
  combined_sites = FALSE,
  silent = TRUE,
  meme_path = NULL,
  ...
)

\method{runMeme}{list}(
  input,
  control = NA,
  outdir = "auto",
  alph = "dna",
  parse_genomic_coord = TRUE,
  combined_sites = FALSE,
  silent = TRUE,
  meme_path = NULL,
  ...
)

\method{runMeme}{BStringSetList}(
  input,
  control = NA,
  outdir = "auto",
  alph = "dna",
  parse_genomic_coord = TRUE,
  combined_sites = FALSE,
  silent = TRUE,
  meme_path = NULL,
  ...
)

\method{runMeme}{default}(
  input,
  control = NA,
  outdir = "auto",
  alph = "dna",
  parse_genomic_coord = TRUE,
  combined_sites = FALSE,
  silent = TRUE,
  meme_path = NULL,
  ...
)
}
\arguments{
\item{input}{path to fasta, Biostrings::BStringSet list, or list of
Biostrings::BStringSet (can generate using \code{get_sequence()})}

\item{control}{any data type as in \code{input}, or a character vector of
\code{names(input)} to use those regions as control sequences. Using sequences
as background requires an alternative objective function. Users must pass a non-default value of
\code{objfun} to \code{...} if using a non-NA control set (default: NA)}

\item{outdir}{(default: "auto") Directory where output data will be stored.}

\item{alph}{one of c("dna", "rna", "protein") or path to alphabet file (default: "dna").}

\item{parse_genomic_coord}{\code{logical(1)} whether to parse genomic coordinates
from fasta headers. Requires headers are in the form: "chr:start-end", or
will result in an error. Automatically set to \code{FALSE} if \code{alph = "protein"}. This setting only needs to be changed if using a custom-built
fasta file without genomic coordinates in the header.}

\item{combined_sites}{\code{logical(1)} whether to return combined sites
information (coerces output to list) (default: FALSE)}

\item{silent}{Whether to suppress printing stdout to terminal (default: TRUE)}

\item{meme_path}{path to "meme/bin/". If unset, will use default search
behavior:
\enumerate{
\item \code{meme_path} setting in \code{options()}
\item \code{MEME_PATH} setting in \code{.Renviron} or \code{.bashrc}
}}

\item{...}{additional arguments passed to MEME (see below)}
}
\value{
MEME results in universalmotif_df format (see:
\code{\link[universalmotif:tidy-motifs]{universalmotif::to_df()}}). \code{sites_hits} is a nested data.frame
column containing the position within each input sequence of matches to the
identified motif.
}
\description{
MEME performs \emph{de-novo} discovery of ungapped motifs present in the input
sequences. It can be used in both discriminative and non-discriminative
modes.
}
\details{
Note that MEME can take a long time to run. The more input sequences used,
the wider the motifs searched for, and the more motifs MEME is asked to
discover will drastically affect runtime. For this reason, MEME usually
performs best on a few (<50) short (100-200 bp) sequences, although this is
not a requirement. Additional details on how data size affects runtime can be
found \href{https://groups.google.com/g/meme-suite/c/7b7PBr6RzJk}{here}.

MEME works best when specifically tuned to the analysis question. The default
settings are unlikely to be ideal. It has several complex arguments
\href{http://meme-suite.org/doc/meme.html}{documented here}, which \code{runMeme()}
accepts as R function arguments (see details below).

If discovering motifs within ChIP-seq, ATAC-seq, or similar peaks, MEME may perform
best if using sequences flaking the summit (the site of maximum signal) of
each peak rather than the center. ChIP-seq or similar data can also benefit
from setting \verb{revcomp = TRUE, minw = 5, maxw = 20}. For more tips on using
MEME to analyze ChIP-seq data, see the following
\href{https://groups.google.com/forum/#\%21topic/meme-suite/rIbjIHbcpAE}{tips page}.
\subsection{Additional arguments}{

\code{\link[=runMeme]{runMeme()}} accepts all valid arguments to meme as arguments passed to \code{...}.
For flags without values, pass them as \code{flag = TRUE}. The \code{dna}, \code{rna}, and
\code{protein} flags should instead be passed to the \code{alph} argument of
\code{\link[=runMeme]{runMeme()}}.  The arguments passed to MEME often have many interactions
with each other, for a detailed description of each argument see
\href{meme-suite.org/doc/meme.html}{MEME Commandline Documentation}.
}
}
\section{Citation}{
If you use \code{runMeme()} in your analysis, please cite:

Timothy L. Bailey and Charles Elkan, "Fitting a mixture model by expectation
maximization to discover motifs in biopolymers", Proceedings of the Second
International Conference on Intelligent Systems for Molecular Biology, pp.
28-36, AAAI Press, Menlo Park, California, 1994.
\href{https://tlbailey.bitbucket.io/papers/ismb94.pdf}{pdf}
}

\section{Licensing}{
The MEME Suite is free for non-profit use, but for-profit users should purchase a
license. See the \href{http://meme-suite.org/doc/copyright.html}{MEME Suite Copyright Page} for details.
}

\examples{
if (meme_is_installed()) {
seqs <- universalmotif::create_sequences("CCRAAAW", seqnum = 4)
names(seqs) <- 1:length(seqs)
runMeme(seqs, parse_genomic_coord = FALSE)

}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_xml.R
\name{importStremeXML}
\alias{importStremeXML}
\title{Import Streme output from previous run}
\usage{
importStremeXML(streme_xml_path)
}
\arguments{
\item{streme_xml_path}{path to streme.xml file}
}
\value{
data.frame with statistics for each discovered motif. The `motifs`
  column contains a universalmotif object representation in PCM format of
  each DREME motif. If no motifs are discovered, returns NULL.
}
\description{
Import Streme output from previous run
}
\examples{
streme_xml <- system.file("extdata", "streme.xml", package = "memes")
importStremeXML(streme_xml)
}
\seealso{
[runStreme()]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{meme_is_installed}
\alias{meme_is_installed}
\title{Returns logical vector indicating valid MEME-Suite install status}
\usage{
meme_is_installed(path = NULL)
}
\arguments{
\item{path}{optional path to "meme/bin/". If unset, will follow the search
heirarchy listed above.}
}
\value{
\code{logical(1)} indicating whether meme is installed with all supported utilities
}
\description{
Checks for a valid meme install using same heirarchy as \code{check_meme_install()}.
Returns \code{TRUE} if all supported utilities are found in the meme install
location, \code{FALSE} if not.
}
\details{
The search heirarchy is as follows:
\enumerate{
\item the \code{meme_path} function argument if set
\item the \code{meme_bin} option
\item the \code{MEME_BIN} environment variable
\item the default location at \verb{~/meme/bin}
}
}
\examples{
meme_is_installed()
}
\seealso{
\code{\link[=check_meme_install]{check_meme_install()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ame.R, R/generics.R
\name{runAme.list}
\alias{runAme.list}
\alias{runAme.BStringSetList}
\alias{runAme.default}
\alias{runAme}
\title{Motif enrichment using AME}
\usage{
\method{runAme}{list}(
  input,
  control = "shuffle",
  outdir = "auto",
  method = "fisher",
  database = NULL,
  meme_path = NULL,
  sequences = FALSE,
  silent = TRUE,
  ...
)

\method{runAme}{BStringSetList}(
  input,
  control = "shuffle",
  outdir = "auto",
  method = "fisher",
  database = NULL,
  meme_path = NULL,
  sequences = FALSE,
  silent = TRUE,
  ...
)

\method{runAme}{default}(
  input,
  control = "shuffle",
  outdir = "auto",
  method = "fisher",
  database = NULL,
  meme_path = NULL,
  sequences = FALSE,
  silent = TRUE,
  ...
)

runAme(
  input,
  control = "shuffle",
  outdir = "auto",
  method = "fisher",
  database = NULL,
  meme_path = NULL,
  sequences = FALSE,
  silent = TRUE,
  ...
)
}
\arguments{
\item{input}{path to fasta, or DNAstringset (optional: DNAStringSet object
names contain fasta score, required for partitioning mode)}

\item{control}{default: "shuffle", or set to
DNAstringset or path to fasta file to use those sequences in discriminative
mode. If \code{input} is a list of DNAStringSet objects, and \code{control} is set to
a character vector of names in \code{input}, those regions will be used as
background in discriminitive mode and AME will skip running on any
\code{control} entries (NOTE: if \code{input} contains an entry named "shuffle" and
control is set to "shuffle", it will use the \code{input} entry, not the AME
shuffle algorithm). If \code{control} is a Biostrings::BStringSetList (generated
using \code{get_sequence}), all sequences in the list will be combined as the
control set. Set to \code{NA} for partitioning based on input fasta score (see
\code{get_sequence()} for assigning fasta score). If input sequences are not
assigned fasta scores but AME is run in partitioning mode, the sequence
order will be used as the score.}

\item{outdir}{Path to output directory location to save data files. If set to "auto",
will use location of input files if passing file paths, otherwise will
write to a temporary directory. default: "auto"}

\item{method}{default: fisher (allowed values: fisher, ranksum, pearson, spearman, 3dmhg, 4dmhg)}

\item{database}{path to .meme format file, universalmotif list object, dreme
results data.frame, or list() of multiple of these. If objects are assigned names in the list,
that name will be used as the database id in the output. It is highly
recommended you set a name if not using a file path so the database name
will be informative, otherwise the position in the list will be used as the
database id. For example, if the input is: list("motifs.meme",
list_of_motifs), the database id's will be: "motifs.meme" and "2". If the
input is list("motifs.meme", "customMotifs" = list_of_motifs), the database
id's will be "motifs.meme" and "customMotifs".}

\item{meme_path}{path to "meme/bin/" (default: \code{NULL}). Will use default
search behavior as described in \code{check_meme_install()} if unset.}

\item{sequences}{\code{logical(1)} add results from \code{sequences.tsv} to \code{sequences}
list column to returned data.frame. Valid only if method = "fisher". See
\href{http://alternate.meme-suite.org/doc/ame-output-format.html}{AME outputs}
webpage for more information (Default: FALSE).}

\item{silent}{whether to suppress stdout (default: TRUE), useful for debugging.}

\item{...}{additional arguments passed to AME (see AME Flag table below)}
}
\value{
a data.frame of AME enrichment results. If run using a BStringsSetList
input, returns a list of data.frames.
}
\description{
AME identifies known motifs (provided by the user) that are enriched in your input sequences.
}
\details{
AME can be run using several modes:
\enumerate{
\item Discriminative mode: to discover motifs enriched relative to shuffled
input, or a set of provided control sequences
\item Partitioning mode: in which AME uses some biological signal to identify
the association between the signal and motif presence.
}

To read more about how AME works, see the \href{http://meme-suite.org/doc/ame-tutorial.html}{AME Tutorial}

Additional AME arguments

memes allows passing any valid flag to it's target programs via \code{...}. For
additional details for all valid AME arguments, see the \href{http://meme-suite.org/doc/ame.html}{AME Manual} webpage. For convenience, a table
of valid parameters, and brief descriptions of their function are provided
below:\tabular{cccl}{
   AME Flag \tab allowed values \tab default \tab description \cr
   kmer \tab \code{integer(1)} \tab 2 \tab kmer frequency to preserve when shuffling control sequences \cr
   seed \tab \code{integer(1)} \tab 1 \tab seed for random number generator when shuffling control sequences \cr
   scoring \tab "avg", "max", "sum", "totalhits" \tab "avg" \tab Method for scoring a sequence for matches to a PWM (avg, max, sum, totalhits) \cr
   hit_lo_fraction \tab \code{numeric(1)} \tab 0.25 \tab fraction of hit log odds score to exceed to be considered a "hit" \cr
   evalue_report_threshold \tab \code{numeric(1)} \tab 10 \tab E-value threshold for reporting a motif as significantly enriched \cr
   fasta_threshold \tab \code{numeric(1)} \tab 0.001 \tab AME will classify sequences with FASTA scores below this value as positives. Only valid when \verb{method = "fisher", poslist = "pwm", control = NA, fix_partition = NULL}. \cr
   fix_partition \tab \code{numeric(1)} \tab \code{NULL} \tab AME evaluates only the partition of the first N sequences. Only works when \code{control = NA} and \code{poslist = "fasta"} \cr
   poslist \tab "pwm", "fasta" \tab "fasta" \tab When using paritioning mode (\code{control = NA}), test thresholds on either PWM or Fasta score \cr
   log_fscores \tab \code{logical(1)} \tab FALSE \tab Convert FASTA scores into log-space (only used when \code{method = "pearson"}) \cr
   log_pwmscores \tab \code{logical(1)} \tab FALSE \tab Convert PWM scores into log-space (only used for \code{method = "pearson"} or \verb{method = "spearman}) \cr
   lingreg_switchxy \tab \code{logical(1)} \tab FALSE \tab Make the x-points FASTA scores and y-points PWM scores (only used for \code{method = "pearson"} or \verb{method = "spearman}) \cr
   xalph \tab file path \tab \code{NULL(1)} \tab alphabet file to use if input motifs are in different alphabet than input sequences \cr
   bfile \tab "motif", "motif-file", "uniform", path to file \tab \code{NULL} \tab source of 0-order background model. If "motif" or "motif-file" 0-order letter frequencies in the first motif file are used. If "uniform" uses uniform letter frequencies. \cr
   motif_pseudo \tab \code{numeric(1)} \tab 0.1 \tab Addd this pseudocount when converting from frequency matrix to log-odds matrix \cr
   inc \tab \code{character(1)} \tab \code{NULL} \tab use only motifs with names matching this regex \cr
   exc \tab \code{character(1)} \tab \code{NULL} \tab exclude motifs with names matching this regex \cr
}
}
\section{Citation}{
If you use \code{runAme()} in your analysis, please cite:

Robert McLeay and Timothy L. Bailey, "Motif Enrichment Analysis: A unified
framework and method evaluation", BMC Bioinformatics, 11:165, 2010,
doi:10.1186/1471-2105-11-165. \href{http://www.biomedcentral.com/1471-2105/11/165}{full text}
\subsection{Licensing}{

The MEME Suite is free for non-profit use, but for-profit users should purchase a
license. See the \href{http://meme-suite.org/doc/copyright.html}{MEME Suite Copyright Page} for details.
}
}

\examples{
if (meme_is_installed()) {
# Create random named sequences as input for example
seqs <- universalmotif::create_sequences(rng.seed = 123)
names(seqs) <- seq_along(seqs)

# An example path to a motif database file in .meme format
motif_file <- system.file("extdata", "flyFactorSurvey_cleaned.meme", package = "memes")

runAme(seqs, database = motif_file)

# Dreme results dataset for example
dreme_xml <- system.file("extdata", "dreme.xml", package = "memes")
dreme_results <- importDremeXML(dreme_xml)

# database can be set to multiple values like so: 
runAme(seqs, database = list(motif_file, "my_dreme_motifs" = dreme_results))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example_chip_summits}
\alias{example_chip_summits}
\title{Annotated Transcription Factor ChIP-seq summits}
\format{
A GRanges object of ChIP summit position with 2 metadata columns
 \describe{
\item{peak_binding_description}{Binding profiles between Early and Late E93
were compared. Peaks are annotated as whether they are bound in Early wings
only ("ectopic"), both Early and Late wings ("entopic"), or only bound in
Late wings ("orphan").}
  \item{e93_sensitive_behavior}{change in chromatin accessibility in response to E93 binding: Increasing, Decreasing, or Static}
 }
}
\source{
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE141738&format=file&file=GSE141738%5Funion%5Fchip%5Fpeaks%5Fannotated%2Ecsv%2Egz
}
\usage{
example_chip_summits
}
\description{
ChIP-seq summit positions on Drosophila melanogaster chromosome 3 for the
transcription factor E93 using a union set of peaks from third-instar larval
wings ("Early") and 24 hour Pupal ("Late") wings.
}
\details{
E93 is a transcription factor normally present only in Late wings. An
experimental perturbation precociously expressed E93 during Early stages.
Binding profiles between Early and Late E93 were compared. Peaks are
annotated as whether they are bound in Early wings only ("ectopic"), both
Early and Late wings ("entopic"), or only bound in Late
wings ("orphan").

DNA elements can be made "open" or "closed" in response to binding of
transcription factors like E93. Accessibility of E93 binding sites before and
after E93 expression was measured using FAIRE-seq. ChIP peaks are annotated
by how their accessibility changes in response to E93 binding . Peaks can
become more open ("Increasing"), more closed ("Decreasing"), or unchanged in
accessibility ("Static"). These experiments demonstrate a causal relationship
between E93 binding and both opening and closing of DNA elements.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match_utils.R
\name{force_best_match}
\alias{force_best_match}
\title{Force best tomtom match by id}
\usage{
force_best_match(res, matches)
}
\arguments{
\item{res}{results from runTomTom}

\item{matches}{named vector where name is the input motif name, and value is
the match_name to use as the new best match}
}
\value{
`res` with new best_* columns and re-ranked tomtom data in the
  `tomtom` list column for the updated entries.
}
\description{
Although TomTom assigns a best match, this is not always the most
biologically relevant match. In these cases, it is useful to "force" the best
match to another lower ranked, but still significant TomTom match. This
function allows users to select a new best match motif from the set of
lower-ranked matches in the `tomtom` list column. This function also reorders
the `tomtom` data.frame such that the forced match is the first row of the
`tomtom` entry.
}
\examples{
if (meme_is_installed()){
motif <- universalmotif::create_motif("CCRAAAW", name = "example_motif")
db <- system.file("extdata", "flyFactorSurvey_cleaned.meme", package = "memes")
res <- runTomTom(motif, database = db)
res$best_match_name
res2 <- force_best_match(res, c("example_motif" = "Eip93F_SANGER_10"))
res2$best_match_name
}
}
\seealso{
[update_best_match()]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fimo.R
\name{runFimo}
\alias{runFimo}
\title{Find instances of motifs using FIMO}
\usage{
runFimo(
  sequences,
  motifs,
  bfile = "motif",
  outdir = "auto",
  parse_genomic_coord = TRUE,
  skip_matched_sequence = FALSE,
  max_strand = TRUE,
  text = TRUE,
  meme_path = NULL,
  silent = TRUE,
  ...
)
}
\arguments{
\item{sequences}{path to fasta file, or stringset input.}

\item{motifs}{path to .meme format file, or universalmotif/universalmotif list input.}

\item{bfile}{path to background file, or special values: "motif" to use
0-order frequencies contained in the motif, or "uniform" to use a uniform
letter distribution. (default: "motif")}

\item{outdir}{output directory location. Only used if text = FALSE. Default:
"auto" to autogenerate directory name. Note: if not using a fasta path as
input, this will be a temporary location unless explicity set.}

\item{parse_genomic_coord}{\code{logical(1)} whether to parse genomic position
from fasta headers. Fasta headers must be UCSC format positions (ie
"chr:start-end"), but base 1 indexed (GRanges format). If names of fasta
entries are genomic coordinates and parse_genomic_coord == TRUE, results
will contain genomic coordinates of motif matches, otherwise FIMO will return
relative coordinates (i.e. positions from 1 to length of the fasta entry).}

\item{skip_matched_sequence}{\code{logical(1)} whether or not to include the DNA
sequence of the match. Default: \code{FALSE}. Note: jobs will complete faster if
set to \code{TRUE}. \code{add_sequence()} can be used to lookup the sequence after data import if
\code{parse_genomic_coord} is \code{TRUE}, so setting this flag is not strictly needed.}

\item{max_strand}{if match is found on both strands, only report strand with
best match (default: TRUE).}

\item{text}{\code{logical(1)} (default: \code{TRUE}). No output files will be created
on the filesystem. The results are unsorted and no q-values are computed.
This setting allows fast searches on very large inputs. When set to \code{FALSE}
FIMO will discard 50\% of the lower significance matches if >100,000 matches are
detected. \code{text = FALSE} will also incur a performance penalty because it
must first read a file to disk, then read it into memory. For these reasons,
I suggest keeping \code{text = TRUE}.}

\item{meme_path}{path to \verb{meme/bin/} (optional). Defaut: \code{NULL}, searches
"MEME_PATH" environment variable or "meme_path" option for path to "meme/bin/".}

\item{silent}{\code{logical(1)} whether to suppress stdout/stderr printing to
console (default: TRUE). If the command is failing or giving unexpected
output, setting \code{silent = FALSE} can aid troubleshooting.}

\item{...}{additional commandline arguments to fimo. See the FIMO Flag table below.}
}
\value{
GRanges object containing positions of each match. Note: if
\code{parse_genomic_coords = FALSE}, each \code{seqnames} entry will be the full fasta
header, and start/end will be the relative position within that sequence of the
match. The GRanges object has the following additional \code{mcols}:
* motif_id = primary name of matched motif
* motif_alt_id = alternate name of matched motif
* score = score of match (higher score is a better match)
* pvalue = pvalue of the match
* qvalue = qvalue of the match
* matched_sequence = sequence that matches the motif
}
\description{
FIMO scans input sequences to identify the positions of matches to each input
motif. FIMO has no sequence length or motif number restrictions.
}
\details{
Additional arguments passed to \code{...}. See: \href{http://meme-suite.org/doc/fimo.html?man_type=web}{Fimo web manual}
for a complete description of FIMO flags.\tabular{cccl}{
   FIMO Flag \tab allowed values \tab default \tab description \cr
   alpha \tab \code{numeric(1)} \tab 1 \tab alpha for calculating position-specific priors. Represents fraction of sites that are binding sites of TF of interest. Used in conjunction with \code{psp} \cr
   bfile \tab "motif", "motif-file", "uniform", file path, \tab "motif" \tab If "motif" or "motif-file", use 0-order letter frequencies from motif. "uniform" sets uniform letter frequencies. \cr
   max_stored_scores \tab \code{integer(1)} \tab NULL \tab maximum number of scores to be stored for computing q-values. used when \code{text = FALSE}, see FIMO webpage for details \cr
   motif_pseudo \tab \code{numeric(1)} \tab 0.1 \tab pseudocount added to motif matrix \cr
   no_qvalue \tab \code{logical(1)} \tab FALSE \tab only needed when \code{text = FALSE}, do not compute q-value for each p-value \cr
   norc \tab \code{logical(1)} \tab FALSE \tab Do not score reverse complement strand \cr
   prior_dist \tab file path \tab NULL \tab file containing binned distribution of priors \cr
   psp \tab file path \tab NULL \tab file containing position specific priors. Requires \code{prior_dist} \cr
   qv_thresh \tab \code{logical(1)} \tab FALSE \tab use q-values for the output threshold \cr
   thresh \tab \code{numeric(1)} \tab \code{1e-4} \tab output threshold for returning a match, only matches with values less than \code{thresh} are returned. \cr
}

\subsection{Licensing}{

The MEME Suite is free for non-profit use, but for-profit users should purchase a
license. See the \href{http://meme-suite.org/doc/copyright.html}{MEME Suite Copyright Page} for details.
}
}
\section{Citation}{
If you use \code{runFimo()} in your analysis, please cite:

Charles E. Grant, Timothy L. Bailey, and William Stafford Noble, "FIMO:
Scanning for occurrences of a given motif", Bioinformatics, 27(7):1017-1018,
2011. \href{http://bioinformatics.oxfordjournals.org/content/early/2011/02/16/bioinformatics.btr064.full}{full text}
}

\examples{
if (meme_is_installed()){
# Generate some example input sequences
seq <- universalmotif::create_sequences()
# sequences must have values in their fasta headers
names(seq) <- seq_along(seq)
# Create random example motif to search for
motif <- universalmotif::create_motif()

# Search for motif in sequences
# parse_genomic_coord set to FALSE since fasta headers aren't in "chr:start-end" format.
runFimo(seq, motif, parse_genomic_coord = FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example_ame_large}
\alias{example_ame_large}
\title{runAme() output for example_chip_summits split by binding description}
\format{
a list of runAme() results data.frames
}
\usage{
example_ame_large
}
\description{
AME was run for "ectopic", "entopic", and "orphan" sites using shuffled background.
}
\details{
see `vignette("integrative_analysis", package = "memes")` for details.
}
\examples{
# Data can be combined into 1 large data.frame using:
dplyr::bind_rows(example_ame_large, .id = "binding_type")
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_duplicate_motifs.R
\name{remove_duplicate_motifs}
\alias{remove_duplicate_motifs}
\title{Remove duplicated motif entries}
\usage{
remove_duplicate_motifs(x)
}
\arguments{
\item{x}{a universalmotif list or universalmotif_df}
}
\value{
A deduplicated list or universalmotif_df
}
\description{
This function identifies motif matrices which are duplicated in a
universalmotif list or universalmotif_df and removes them. This operation
ignores motif metadata and operates by removing all entries whose motif
matrices are identical. The first instance of a duplicated motif in the input
list is the one returned.
}
\examples{
motif <- universalmotif::create_motif()
duplicated <- c(motif, motif)
remove_duplicate_motifs(duplicated)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example_peaks}
\alias{example_peaks}
\title{Example ChIP-seq peaks}
\format{
An object of class \code{GRanges} of length 10.
}
\source{
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE141738&format=file&file=GSE141738%5Funion%5Fchip%5Fpeaks%5Fannotated%2Ecsv%2Egz
}
\usage{
example_peaks
}
\description{
10 ChIP-seq peaks from GSE141738
}
\details{
A small number of transcription factor ChIP-seq peaks as a GRanges object,
taken from [GSE141738](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141738)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/memes-package.R
\docType{package}
\name{memes-package}
\alias{memes}
\alias{memes-package}
\title{memes: motif matching, comparison, and de novo discovery using the MEME Suite}
\description{
A seamless interface to the MEME Suite family of tools for motif
    analysis. 'memes' provides data aware utilities for using GRanges objects as
    entrypoints to motif analysis, data structures for examining & editing motif
    lists, and novel data visualizations. 'memes' functions and data structures are
    amenable to both base R and tidyverse workflows.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://snystrom.github.io/memes/}
  \item \url{https://github.com/snystrom/memes}
  \item Report bugs at \url{https://github.com/snystrom/memes/issues}
}

}
\author{
\strong{Maintainer}: Spencer Nystrom \email{nystromdev@gmail.com} (\href{https://orcid.org/0000-0003-1000-1579}{ORCID}) [copyright holder]

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example_fimo}
\alias{example_fimo}
\title{Example runFimo() output}
\format{
A GRanges object of E93 motif positions within 100bp windows of `example_chip_summits`
}
\usage{
example_fimo
}
\description{
Run using 100bp windows around `example_chip_summits`, using E93 motif as database.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example_dreme_tomtom}
\alias{example_dreme_tomtom}
\title{Example runDreme() output after passing to runTomTom()}
\format{
a runDreme results data.frame with runTomTom results columns
}
\usage{
example_dreme_tomtom
}
\description{
Result when running `runTomTom(example_dreme)` using FlyFactorSurvey as database.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ame_methods.R
\name{ame_compare_heatmap_methods}
\alias{ame_compare_heatmap_methods}
\title{Compare AME heatmap methods}
\usage{
ame_compare_heatmap_methods(ame, group, value = -log10(adj.pvalue))
}
\arguments{
\item{ame}{ame results data.frame}

\item{group}{optional name of group to split results by}

\item{value}{value to compare to "normalize" method (default:
-log10(adj.pvalue))}
}
\value{
a cowplot 2 panel plot comparing the distribution of `value` to
  normalized rank values
}
\description{
This helper function allows the user to visualize the distribution of their
AME results data on different scales to help understand the implications of
using different values in `plot_ame_heatmap()`
}
\examples{
data("example_ame", package = "memes")
ame_compare_heatmap_methods(example_ame$Decreasing)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/streme.R
\name{runStreme}
\alias{runStreme}
\title{Denovo motif discovery of target regions using STREME}
\usage{
runStreme(
  input,
  control,
  outdir = "auto",
  objfun = "de",
  alph = "dna",
  meme_path = NULL,
  silent = TRUE,
  ...
)
}
\arguments{
\item{input}{regions to scan for motifs. If using `objfun = "cd"` to test for
centrally enriched motifs, be sure to include sufficient flanking sequence
(e.g. +/- 500bp) for an accurate estimate. Can be any of:
- path to fasta file
- DNAStringSet object (can be generated from GRanges using `get_sequence()`)
- List of DNAStringSet objects (generated from `get_sequence()`)
- *NOTE:* if using StringSet inputs, each entry must be named (set with `names()`).
- *NOTE:* If you want to retain the raw streme output files, you must use a
path to fasta file as input, or specify an "outdir"}

\item{control}{regions to use as background for motif search. These should
have a similar length distribution as the input sequences. Can be any of:
- path to fasta file
- DNAStringSet object (can be generated from GRanges using get_sequence)
- A Biostrings::BStringSetList (generated using `get_sequence`), in which
case all sequences in the list will be combined as the control set.
- if `input` is a list of DNAStringSet objects, a character vector of names
in `input` will use those sequences as background. runstreme will not scan
those regions as input.
- "shuffle" to use streme's built-in dinucleotide shuffle feature (NOTE: if
`input` is a list object with an entry named "shuffle", the list entry will
be used instead).
Optionally can also pass `seed = <any number>` to `...` to use as the random
seed during shuffling. If no seed is passed, streme will use 0 as the random
seed, so results will be reproducible if rerunning.}

\item{outdir}{path to output directory of streme files, or "auto" to
autogenerate path. Default: location of input fasta in dir named
"\<input\>_vs_\<control\>". If input is DNAstringset, will be temporary
path. This means that if you want to save the raw output files, you must
use fasta files as input or use an informative (and unique) outdir name.
memes will **not check** if it overwrites files in a directory. Directories
will be recursively created if needed. (default: "auto")}

\item{objfun}{one of c("de", "cd"). Default: "de" for differential
enrichment. "cd" for central distance (control must be set to NA for "cd").}

\item{alph}{one of c("dna", "rna", "protein") or a path to a MEME format alph file. (default: "dna")}

\item{meme_path}{path to "meme/bin"}

\item{silent}{Whether to suppress printing stdout & stderr to console
(default: TRUE). Warnings are always printed regardless of this setting.}

\item{...}{pass any commandline options as R function arguments. For a
complete list of STREME options, see 
[the STREME manual](https://meme-suite.org/meme/doc/streme.html).}
}
\value{
a `universalmotif_df` of STREME Motifs
}
\description{
STREME discovers short, ungapped, *de-novo* motifs that are enriched or
relatively enriched relative to a control set of sequences. STREME can be run
to discover motifs relative to a shuffled set of input sequences, against a
separately provided set of "control" sequences, or to determine whether
motifs are centrally enriched within input sequences.
}
\details{
Properly setting the `control` parameter is key to discovering biologically
relevant motifs. Often, using `control = "shuffle"` will produce a suboptimal
set of motifs; however, some discriminative analysis designs don't have
proper "control" regions other than to shuffle.


If you have fewer than 50 sequences, consider using [runMeme()] instead.



 # Citation

 If you use `runStreme()` in your analysis, please cite:

 Timothy L. Bailey, "STREME: Accurate and versatile sequence motif
 discovery", Bioinformatics, 2021.
 https://doi.org/10.1093/bioinformatics/btab203
 
 # Licensing
 The MEME Suite is free for non-profit use, but for-profit users should purchase a
 license. See the [MEME Suite Copyright Page](http://meme-suite.org/doc/copyright.html) for details.
}
\seealso{
`?universalmotif::tidy-motifs`
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sequence_heatmap.R
\name{plot_sequence_heatmap}
\alias{plot_sequence_heatmap}
\title{Visualize a heatmap of aligned sequences}
\usage{
plot_sequence_heatmap(
  sequence,
  title = NULL,
  logo = TRUE,
  alph = c("DNA", "RNA", "AA"),
  title_hjust = 0,
  heights = c(1, 5),
  legend = "none"
)
}
\arguments{
\item{sequence}{character vector of sequences, plot will be ranked in order
  of the sequences. Each sequence must be equal length. Alternately, sequence
can be a named list in which case each plot will be titled by the names of
the list entries.}

\item{title}{title of the plot. Default: NULL. If sequence is a named list of
sequences, title defaults to the list entry names. Set to NULL to override
this behavior. To use a different title than the list entry name, pass a
vector of names to `title`.}

\item{logo}{whether to include a sequence logo above the heatmap}

\item{alph}{alphabet colorscheme to use. One of: DNA, RNA, AA.}

\item{title_hjust}{value from 0 to 1 determining the horizontal justification
of the title. Default: 0.}

\item{heights}{ratio of logo:heatmap heights. Given as: c(logo_height,
heatmap_height). Values are not absolute. Ignored when logo = FALSE.}

\item{legend}{passed to ggplot2::theme(legend.position). Default: "none".
Values can be: "none", "left", "right", "top", "bottom", or coordinates in
c(x,y) format.}
}
\value{
a ggplot object of the sequence heatmap ranked by the order of
  sequences
}
\description{
Sometimes it is useful to visualize individual motif matches in aggregate to
understand how sequence variability contributes to motif matches. This
function creates a heatmap where each row represents a single sequence and
each column represents a position. Cells are colored by the sequence at that
position. Sequences are optionally aggregated into a sequence logo aligned in
register with the heatmap to visualize how sequence variability contributes
to motif makeup.
}
\examples{
data(example_fimo, package = "memes")
genome <- BSgenome.Dmelanogaster.UCSC.dm3::BSgenome.Dmelanogaster.UCSC.dm3
motifs <- add_sequence(example_fimo, genome)
plot_sequence_heatmap(motifs$sequence)

# Use on named list
sequences <- list("set 1" = motifs$sequence[1:100], 
                  "set 2" = motifs$sequence[101:200])
plot_sequence_heatmap(sequences)

# Use different titles for list input
plot_sequence_heatmap(sequences, title = c("A", "B"))
}
\seealso{
runFimo
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meme.R
\name{importMeme}
\alias{importMeme}
\title{Import MEME results}
\usage{
importMeme(meme_txt, parse_genomic_coord = FALSE, combined_sites = FALSE)
}
\arguments{
\item{meme_txt}{path to "meme.txt" output}

\item{parse_genomic_coord}{whether to parse sequence headers into genomic
coordinates for motif position information, only works if fasta files were
written such that the sequence headers are in the form: "chr:start-end", or
some variation of this form (delimiters can be any of: "[^[:alnum:]]+" (ie
non-alphanumeric characters)), (default = FALSE).}

\item{combined_sites}{whether to add `combined_sites` output which contains
coordinates of each sequence, the motif sequence (if `parse_genomic_coord =
TRUE`), and the `diagram` column raw output from MEME indicating the
relative locations of motifs in the sequence.}
}
\value{
MEME results in universalmotif data.frame format (see:
  [as_universalmotif_dataframe()]). `sites_hits` is a nested data.frame
  column containing the position within each input sequence of matches to the
  identified motif.
}
\description{
This is a light wrapper around [universalmotif::read_meme()] that imports
MEME results as universalmotif data.frame. If MEME is run with genomic
coordinates in the fasta header, in "chr:start-end" format (base 1 indexed),
the genomic coordinates of the motif match from input sequences can be
parsed from the header.
}
\examples{
example_meme_txt <- system.file("extdata", "meme_full.txt", package = "universalmotif")
importMeme(example_meme_txt)
}
\seealso{
[runMeme()] [universalmotif::read_meme()]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ame.R
\name{importAme}
\alias{importAme}
\title{Parse AME output}
\usage{
importAme(path, method = "fisher", sequences = FALSE)
}
\arguments{
\item{path}{path to ame results file ("ame.tsv")}

\item{method}{ame run method used (one of: c("fisher", "ranksum", "dmhg3",
"dmhg4", "pearson", "spearman")). Default: "fisher".}

\item{sequences}{FALSE or path to sequences file (only valid for method = "fisher")}
}
\value{
data.frame with method-specific results. See [AME
  results](http://meme-suite.org/doc/ame-output-format.html) webpage for more
  information. If sequences is set to a path to the sequences.tsv and method
  = "fisher", the list-column `sequences` will be appended to resulting
  data.frame.
}
\description{
This imports AME results using the "ame.tsv" output, and optionally the
"sequences.tsv" output if run with "method = fisher". AME results differ
based on the method used, thus this must be set on import or the column
names will be incorrect.
}
\examples{
ame_tsv <- system.file("extdata", "ame.tsv", package = "memes", mustWork = TRUE)
importAme(ame_tsv)
}
\seealso{
[runAme()]
}
\concept{import}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match_utils.R
\name{update_best_match}
\alias{update_best_match}
\title{Update best match info by ranking of tomtom data}
\usage{
update_best_match(res)
}
\arguments{
\item{res}{results from runTomTom}
}
\value{
`res` with updated best_* columns
}
\description{
This function updates the best_match columns based on the rankings on the
tomtom list data. By re-ordering the entries of a `tomtom` object, the
best_match columns can be updated to reflect the new rankings using
[update_best_match()], where the first row of the `tomtom` data.frame is
selected as the best match.
}
\examples{
data("example_dreme_tomtom")
# best match is "CG2052_SANGER_2.5"
example_dreme_tomtom$best_match_name[1]
# reorder the `tomtom` data.frame
example_dreme_tomtom$tomtom[[1]] <- dplyr::arrange(example_dreme_tomtom$tomtom[[1]],
                                                   dplyr::desc(match_eval))
# update_best_match will use the new order of rows, taking the top row as the new best match
new_res <- update_best_match(example_dreme_tomtom)
# best match is now altered:
new_res$best_match_name[1]
}
\seealso{
[force_best_match()]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example_rnaseq}
\alias{example_rnaseq}
\title{RNAseq data from Early and Late Drosophila wings}
\format{
A data.frame of RNAseq FPKM data
\describe{
 \item{symbol}{The FlyBase gene symbol}
 \item{time}{Developmental stage of RNAseq experiment}
 \item{fpkm}{RNAseq count in Fragments per Kilobase Million (FPKM)}
}
}
\source{
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97956&format=file&file=GSE97956%5FgeneFPKMs%5FL3%5F24hr%5F44hr%2Exlsx"
}
\usage{
example_rnaseq
}
\description{
These data are a subset of RNAseq counts from the full FPKM table in
GSE97956. Includes counts for all Drosophila transcription factors and ~200
other random genes.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{runDreme}
\alias{runDreme}
\title{Denovo motif discovery of target regions using DREME}
\usage{
runDreme(input, control, outdir = "auto", meme_path = NULL, silent = TRUE, ...)
}
\arguments{
\item{input}{regions to scan for motifs. Can be any of:
\itemize{
\item path to fasta file
\item DNAStringSet object (can be generated from GRanges using \code{get_sequence()})
\item List of DNAStringSet objects (generated from \code{get_sequence()})
\item \emph{NOTE:} if using StringSet inputs, each entry must be named (set with \code{names()}).
\item \emph{NOTE:} If you want to retain the raw dreme output files, you must use a
path to fasta file as input, or specify an "outdir"
}}

\item{control}{regions to use as background for motif search. Can be any of:
\itemize{
\item path to fasta file
\item DNAStringSet object (can be generated from GRanges using get_sequence)
\item A Biostrings::BStringSetList (generated using \code{get_sequence}), in which
case all sequences in the list will be combined as the control set.
\item if \code{input} is a list of DNAStringSet objects, a character vector of names
in \code{input} will use those sequences as background. runDreme will not scan
those regions as input.
\item "shuffle" to use dreme's built-in dinucleotide shuffle feature (NOTE: if
\code{input} is a list object with an entry named "shuffle", the list entry will
be used instead).
Optionally can also pass \verb{seed = <any number>} to \code{...} to use as the random
seed during shuffling. If no seed is passed, dreme will use 1 as the random
seed, so results will be reproducible if rerunning. \strong{NOTE:} beware
system-specific differences. As of v5, dreme will compile using the default
python installation on a system (either python2.7 or python3). The random
number generator changed between python2.7 and python3, so results will not
be reproducible between systems using python2.7 vs 3.
}}

\item{outdir}{path to output directory of dreme files, or "auto" to autogenerate path. Default: location of
input fasta in dir named "\<input\>\emph{vs}\<control\>". If input is
DNAstringset, will be temporary path. This means that if you want to save
the raw output files, you must use fasta files as input or use an
informative (and unique) outdir name. memes will \strong{not check} if it
overwrites files in a directory. Directories will be recursively created if needed.}

\item{meme_path}{optional, path to "meme/bin/" on your local machine.
runDreme will search 3 places in order for meme if this flag is unset:
\enumerate{
\item the option "meme_bin" (set with options(meme_bin = "path/to/meme/bin"))
\item the environment variable "MEME_PATH" (set in .Renviron)
\item "~/meme/bin/" as the default location
}
\itemize{
\item If the user sets meme_path in the function call, this value will always be used
}}

\item{silent}{whether to suppress printing dreme stdout as a message when
finishing with no errors. Can be useful for troubleshooting in situations
where no motifs are discovered, but command completes successfully.
(default: TRUE)}

\item{...}{dreme flags can be passed as R function arguments to use
non-default behavior. For a full list of valid arguments, run your local
install of dreme -h, or visit the dreme documentation
\href{http://meme-suite.org/doc/dreme.html}{website}. See list below for aliases
of common flags. To set flags with no values (ex. \code{-dna}), pass the
argument as a boolean value (ex. \code{dna = TRUE}).}
}
\value{
\code{universalmotif_df} with statistics for each discovered motif. The \code{motif}
column contains a universalmotif object representation in PCM format of
each DREME motif. If no motifs are discovered, returns NULL. The column
values are as follows:
\itemize{
\item rank = ranked order of discovered motif
\item name = primary name of motif
\item altname = alternative name of motif
\item seq = consensus sequence of the motif
\item length = length of discovered motif
\item nsites = number of times the motif is found in input sequences
\item positive_hits = number of sequences in input containing at least 1 of the motif
\item negative_hits = number of sequences in control containing at least 1 of the motif
\item pval = p-value
\item eval = E-value
\item unerased_eval = Unerased E-Value
\item positive_total = number of sequences in input
\item negative_total = number of sequences in control
\item pos_frac = fraction of positive sequences with a hit
\item neg_frac = fraction of negative sequences with a hit
\item motif = a universalmotif object of the discovered motif
}
}
\description{
DREME discovers short, ungapped, \emph{de-novo} motifs that are relatively
enriched relative to a control set of sequences. DREME can be run to discover
motifs relative to a shuffled set of input sequences, or against a separately
provided set of "control" sequences.
}
\details{
Properly setting the \code{control} parameter is key to discovering biologically
relevant motifs. Often, using \code{control = "shuffle"} will produce a suboptimal
set of motifs; however, some discriminative analysis designs don't have
proper "control" regions other than to shuffle.

As of MEME version 5.2.0, DREME is deprecated. Consider \code{\link[=runStreme]{runStreme()}} instead.

In addition to allowing any valid flag of dreme to be passed to \code{...}, we
provide a few user-friendly aliases for common flags which are more readable (see list below).
For example, e = 1 will use a max evalue cutoff of 1. This is equivalent to
setting evalue = 1. For additional details about each DREME flag, see the
\href{http://meme-suite.org/doc/dreme.html}{DREME Manual Webpage}.

List of values which can be passed to \code{...}:
\strong{NOTE:} values must be referred to using their name in the "memes alias"
column, not the "DREME Flag" column.\tabular{cclc}{
   memes alias \tab DREME Flag \tab description \tab default \cr
   nmotifs \tab m \tab max number of motifs to discover \tab NULL \cr
   sec \tab t \tab max number of seconds to run \tab NULL \cr
   evalue \tab e \tab max E-value cutoff \tab 0.05 \cr
   seed \tab s \tab random seed if using "shuffle" as control \tab 1 \cr
   ngen \tab g \tab nuber of REs to generalize \tab 100 \cr
   mink \tab mink \tab minimum motif width to search \tab 3 \cr
   maxk \tab maxk \tab maximum motif width to search \tab 7 \cr
   k \tab k \tab set mink and maxk to this value \tab NULL \cr
   norc \tab norc \tab search only the input strand for sequences \tab FALSE \cr
   dna \tab dna \tab use DNA alphabet \tab TRUE \cr
   rna \tab rna \tab use RNA alphabet \tab FALSE \cr
   protein \tab protein \tab use protein alphabet (NOT RECCOMENDED) \tab FALSE \cr
}
}
\section{Citation}{
If you use \code{runDreme()} in your analysis, please cite:

Timothy L. Bailey, "DREME: Motif discovery in transcription factor ChIP-seq
data", Bioinformatics, 27(12):1653-1659, 2011.
\href{https://academic.oup.com/bioinformatics/article/27/12/1653/257754}{full text}
\subsection{Licensing}{

The MEME Suite is free for non-profit use, but for-profit users should purchase a
license. See the \href{http://meme-suite.org/doc/copyright.html}{MEME Suite Copyright Page} for details.
}
}

\examples{
if (meme_is_installed()) {
# Create random named sequences as input for example
seqs <- universalmotif::create_sequences(rng.seed = 123)
names(seqs) <- seq_along(seqs)

# Runs dreme with default settings, shuffles input as background
runDreme(seqs, "shuffle")

# Runs searching for max 2 motifs, e-value cutoff = 0.1, explicitly using the DNA alphabet
runDreme(seqs, "shuffle", nmotifs = 2, e = 0.1, dna = TRUE)
}
}
