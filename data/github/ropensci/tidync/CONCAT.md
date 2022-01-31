# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who 
contribute through reporting issues, posting feature requests, updating documentation,
submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for
everyone, regardless of level of experience, gender, gender identity and expression,
sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or
imagery, derogatory comments or personal attacks, trolling, public or private harassment,
insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments,
commits, code, wiki edits, issues, and other contributions that are not aligned to this 
Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed 
from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by 
opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the Contributor Covenant 
(http:contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/

<!-- README.md is generated from README.Rmd. Please edit that file -->

# tidync

<!-- badges: start -->

[![Travis-CI Build
Status](http://badges.herokuapp.com/travis/ropensci/tidync?branch=master&env=BUILD_NAME=trusty_release&label=linux)](https://travis-ci.org/ropensci/tidync)
[![Build
Status](http://badges.herokuapp.com/travis/ropensci/tidync?branch=master&env=BUILD_NAME=osx_release&label=osx)](https://travis-ci.org/ropensci/tidync)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/ropensci/tidync?branch=master&svg=true)](https://ci.appveyor.com/project/mdsumner/tidync-nvwwt)
[![Coverage
status](https://codecov.io/gh/ropensci/tidync/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/tidync?branch=master)
[![](https://badges.ropensci.org/174_status.svg)](https://github.com/ropensci/onboarding/issues/174)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/tidync)](https://cran.r-project.org/package=tidync)
[![CRAN\_Download\_Badge](http://cranlogs.r-pkg.org/badges/tidync)](https://cran.r-project.org/package=tidync)
<!-- badges: end -->

The goal of tidync is to ease exploring the contents of a NetCDF source
and to simplify the process of data extraction.

When extracting, data can be accessed as array/s, or in long-form as a
data frame. In contrast to other packages tidync helps reduce the volume
of code required to discover and read the contents of NetCDF, with
simple steps:

  - Connect and summarize `tidync()`.
  - (optionally) Specify source variables `activate()`.
  - (optionally) Specify array sub-setting (slicing) `hyper_filter()`.
  - Read array data in native form `hyper_array()` or long-form
    `hyper_tibble()` or bespoke form `hyper_tbl_cube()`.

NetCDF is **Network Common Data Form** a very common, and very general
way to store and work with scientific array-based data. NetCDF is
defined and provided by
[Unidata](https://www.unidata.ucar.edu/software/netcdf/). R has
(independent) support for NetCDF via the
[ncdf4](https://CRAN.R-project.org/package=ncdf4),
[rhdf5](https://bioconductor.org/packages/release/bioc/html/rhdf5.html),
[RNetCDF](https://CRAN.R-project.org/package=RNetCDF),
[rgdal](https://CRAN.R-project.org/package=rgdal),
[sf](https://CRAN.R-project.org/package=sf) and
[vapour](https://CRAN.R-project.org/package=vapour) packages.

This project uses RNetCDF for the primary access to the NetCDF library,
as well as the ncdf4 package in some cases. The wrapper provided by
[ncmeta](https://CRAN.R-project.org/package=ncmeta) over RNetCDF is used
to obtain information about data sources.

## Installation

Install tidync from CRAN.

``` r
install.packages("tidync")
```

You can install the development version from github with:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/tidync", dependencies = TRUE)
```

The package packages ncdf4 and RNetCDF are required, so first make sure
you can install and use these if it doesn’t work the first time.

``` r
install.packages("ncdf4")
install.packages("RNetCDF")
```

If you have problems, please see the [INSTALL instructions for
RNetCDF](https://CRAN.R-project.org/package=RNetCDF/INSTALL), these
should work as well for ncdf4. Below I note specifics for different
operating systems, notably Ubuntu/Debian where I work the most - these
aren’t comprehensive details but might be helpful.

### Windows

On Windows, everything should be easy as ncdf4 and RNetCDF are supported
by CRAN. The RNetCDF package now includes OpenDAP/Thredds for 64-bit
Windows (not 32-bit), and so tidync will work for those sources too.

### MacOS

On MacOS, it should also be easy as there are binaries for ncdf4 and
RNetCDF available on CRAN. As far as I know, only RNetCDF will support
Thredds.

### Ubuntu/Debian

On Linux you will need at least the following installed by an
administrator, here tested on Ubuntu Xenial 16.04.

``` bash
apt update 
apt upgrade --assume-yes

## Install 3rd parties for NetCDF
apt install libnetcdf-dev libudunits2-dev

## install 3rd parties needed for devtools + openssl git2r httr
apt install libssl-dev
```

Then in R

``` r
install.packages("remotes")
remotes::install_github("ropensci/tidync")
```

More general information about system dependencies libnetcdf-dev and
libudunits2-dev is available from
[Unidata](https://www.unidata.ucar.edu).

## Usage

This is a basic example which shows how to connect to a file.

``` r
file <- system.file("extdata", "oceandata", "S20080012008031.L3m_MO_CHL_chlor_a_9km.nc", package = "tidync")
library(tidync)
tidync(file) 
#> 
#> Data Source (1): S20080012008031.L3m_MO_CHL_chlor_a_9km.nc ...
#> 
#> Grids (4) <dimension family> : <associated variables> 
#> 
#> [1]   D1,D0 : chlor_a    **ACTIVE GRID** ( 9331200  values per variable)
#> [2]   D3,D2 : palette
#> [3]   D0    : lat
#> [4]   D1    : lon
#> 
#> Dimensions 4 (2 active): 
#>   
#>   dim   name  length    min   max start count   dmin  dmax unlim coord_dim 
#>   <chr> <chr>  <dbl>  <dbl> <dbl> <int> <int>  <dbl> <dbl> <lgl> <lgl>     
#> 1 D0    lat     2160  -90.0  90.0     1  2160  -90.0  90.0 FALSE TRUE      
#> 2 D1    lon     4320 -180.  180.      1  4320 -180.  180.  FALSE TRUE      
#>   
#> Inactive dimensions:
#>   
#>   dim   name          length   min   max unlim coord_dim 
#>   <chr> <chr>          <dbl> <dbl> <dbl> <lgl> <lgl>     
#> 1 D2    rgb                3     1     3 FALSE FALSE     
#> 2 D3    eightbitcolor    256     1   256 FALSE FALSE
```

There are two main ways of using tidync, interactively to explore what
is there, and for extraction. The functions `tidync` and `activate` and
`hyper_filter` allow us to hone in on the part/s of the data we want,
and functions `hyper_array`, `hyper_tibble` and `hyper_tbl_cube` give
raw-array or data frames.

### Interactive

Use `tidync()` and `hyper_filter()` to discern what variables and
dimensions are available, and to craft axis-filtering expressions by
value or by index. (Use the name of the variable on the LHS to target
it, use its name to filter by value and the special name `index` to
filter it by its index).

``` r
filename <- system.file("extdata/argo/MD5903593_001.nc", package = "tidync")
## discover the available entities, and the active grid's dimensions and variables
tidync(filename)
#> 
#> Data Source (1): MD5903593_001.nc ...
#> 
#> Grids (16) <dimension family> : <associated variables> 
#> 
#> [1]   D0,D9,D11,D8 : SCIENTIFIC_CALIB_DATE
#> [2]   D6,D9,D11,D8 : PARAMETER
#> [3]   D7,D9,D11,D8 : SCIENTIFIC_CALIB_EQUATION, SCIENTIFIC_CALIB_COEFFICIENT, SCIENTIFIC_CALIB_COMMENT
#> [4]   D6,D9,D8     : STATION_PARAMETERS
#> [5]   D10,D8       : PRES, PRES_QC, PRES_ADJUSTED, PRES_ADJUSTED_QC, PRES_ADJUSTED_ERROR, TEMP, TEMP_QC, TEMP_ADJUSTED, TEMP_ADJUSTED_QC, TEMP_ADJUSTED_ERROR, PSAL, PSAL_QC, PSAL_ADJUSTED, PSAL_ADJUSTED_QC, PSAL_ADJUSTED_ERROR, DOXY, DOXY_QC, DOXY_ADJUSTED, DOXY_ADJUSTED_QC, DOXY_ADJUSTED_ERROR, CHLA, CHLA_QC, CHLA_ADJUSTED, CHLA_ADJUSTED_QC, CHLA_ADJUSTED_ERROR, BBP700, BBP700_QC, BBP700_ADJUSTED, BBP700_ADJUSTED_QC, BBP700_ADJUSTED_ERROR, NITRATE, NITRATE_QC, NITRATE_ADJUSTED, NITRATE_ADJUSTED_QC, NITRATE_ADJUSTED_ERROR    **ACTIVE GRID** ( 986  values per variable)
#> [6]   D1,D8        : DATA_CENTRE
#> [7]   D2,D8        : DATA_STATE_INDICATOR, WMO_INST_TYPE
#> [8]   D3,D8        : PLATFORM_NUMBER, POSITIONING_SYSTEM
#> [9]   D5,D8        : DC_REFERENCE, PLATFORM_TYPE, FLOAT_SERIAL_NO, FIRMWARE_VERSION
#> [10]   D6,D8        : PROJECT_NAME, PI_NAME
#> [11]   D7,D8        : VERTICAL_SAMPLING_SCHEME
#> [12]   D9,D8        : PARAMETER_DATA_MODE
#> [13]   D0           : REFERENCE_DATE_TIME, DATE_CREATION, DATE_UPDATE
#> [14]   D2           : FORMAT_VERSION, HANDBOOK_VERSION
#> [15]   D5           : DATA_TYPE
#> [16]   D8           : CYCLE_NUMBER, DIRECTION, DATA_MODE, JULD, JULD_QC, JULD_LOCATION, LATITUDE, LONGITUDE, POSITION_QC, CONFIG_MISSION_NUMBER, PROFILE_PRES_QC, PROFILE_TEMP_QC, PROFILE_PSAL_QC, PROFILE_DOXY_QC, PROFILE_CHLA_QC, PROFILE_BBP700_QC, PROFILE_NITRATE_QC
#> 
#> Dimensions 14 (2 active): 
#>   
#>   dim   name     length   min   max start count  dmin  dmax unlim coord_dim 
#>   <chr> <chr>     <dbl> <dbl> <dbl> <int> <int> <dbl> <dbl> <lgl> <lgl>     
#> 1 D8    N_PROF        2     1     2     1     2     1     2 FALSE FALSE     
#> 2 D10   N_LEVELS    493     1   493     1   493     1   493 FALSE FALSE     
#>   
#> Inactive dimensions:
#>   
#>    dim   name       length   min   max unlim coord_dim 
#>    <chr> <chr>       <dbl> <dbl> <dbl> <lgl> <lgl>     
#>  1 D0    DATE_TIME      14     1    14 FALSE FALSE     
#>  2 D1    STRING2         2     1     2 FALSE FALSE     
#>  3 D2    STRING4         4     1     4 FALSE FALSE     
#>  4 D3    STRING8         8     1     8 FALSE FALSE     
#>  5 D4    STRING16       16    NA    NA FALSE FALSE     
#>  6 D5    STRING32       32     1    32 FALSE FALSE     
#>  7 D6    STRING64       64     1    64 FALSE FALSE     
#>  8 D7    STRING256     256     1   256 FALSE FALSE     
#>  9 D9    N_PARAM         7     1     7 FALSE FALSE     
#> 10 D11   N_CALIB         1     1     1 FALSE FALSE     
#> 11 D12   N_HISTORY       0    NA    NA TRUE  FALSE     
#> 12 D13   N_VALUES41     41    NA    NA FALSE FALSE

## activate a different grid
grid_identifier <- "D7,D9,D11,D8"
tidync(filename) %>% activate(grid_identifier)
#> 
#> Data Source (1): MD5903593_001.nc ...
#> 
#> Grids (16) <dimension family> : <associated variables> 
#> 
#> [1]   D0,D9,D11,D8 : SCIENTIFIC_CALIB_DATE
#> [2]   D6,D9,D11,D8 : PARAMETER
#> [3]   D7,D9,D11,D8 : SCIENTIFIC_CALIB_EQUATION, SCIENTIFIC_CALIB_COEFFICIENT, SCIENTIFIC_CALIB_COMMENT    **ACTIVE GRID** ( 3584  values per variable)
#> [4]   D6,D9,D8     : STATION_PARAMETERS
#> [5]   D10,D8       : PRES, PRES_QC, PRES_ADJUSTED, PRES_ADJUSTED_QC, PRES_ADJUSTED_ERROR, TEMP, TEMP_QC, TEMP_ADJUSTED, TEMP_ADJUSTED_QC, TEMP_ADJUSTED_ERROR, PSAL, PSAL_QC, PSAL_ADJUSTED, PSAL_ADJUSTED_QC, PSAL_ADJUSTED_ERROR, DOXY, DOXY_QC, DOXY_ADJUSTED, DOXY_ADJUSTED_QC, DOXY_ADJUSTED_ERROR, CHLA, CHLA_QC, CHLA_ADJUSTED, CHLA_ADJUSTED_QC, CHLA_ADJUSTED_ERROR, BBP700, BBP700_QC, BBP700_ADJUSTED, BBP700_ADJUSTED_QC, BBP700_ADJUSTED_ERROR, NITRATE, NITRATE_QC, NITRATE_ADJUSTED, NITRATE_ADJUSTED_QC, NITRATE_ADJUSTED_ERROR
#> [6]   D1,D8        : DATA_CENTRE
#> [7]   D2,D8        : DATA_STATE_INDICATOR, WMO_INST_TYPE
#> [8]   D3,D8        : PLATFORM_NUMBER, POSITIONING_SYSTEM
#> [9]   D5,D8        : DC_REFERENCE, PLATFORM_TYPE, FLOAT_SERIAL_NO, FIRMWARE_VERSION
#> [10]   D6,D8        : PROJECT_NAME, PI_NAME
#> [11]   D7,D8        : VERTICAL_SAMPLING_SCHEME
#> [12]   D9,D8        : PARAMETER_DATA_MODE
#> [13]   D0           : REFERENCE_DATE_TIME, DATE_CREATION, DATE_UPDATE
#> [14]   D2           : FORMAT_VERSION, HANDBOOK_VERSION
#> [15]   D5           : DATA_TYPE
#> [16]   D8           : CYCLE_NUMBER, DIRECTION, DATA_MODE, JULD, JULD_QC, JULD_LOCATION, LATITUDE, LONGITUDE, POSITION_QC, CONFIG_MISSION_NUMBER, PROFILE_PRES_QC, PROFILE_TEMP_QC, PROFILE_PSAL_QC, PROFILE_DOXY_QC, PROFILE_CHLA_QC, PROFILE_BBP700_QC, PROFILE_NITRATE_QC
#> 
#> Dimensions 14 (4 active): 
#>   
#>   dim   name      length   min   max start count  dmin  dmax unlim coord_dim 
#>   <chr> <chr>      <dbl> <dbl> <dbl> <int> <int> <dbl> <dbl> <lgl> <lgl>     
#> 1 D7    STRING256    256     1   256     1   256     1   256 FALSE FALSE     
#> 2 D8    N_PROF         2     1     2     1     2     1     2 FALSE FALSE     
#> 3 D9    N_PARAM        7     1     7     1     7     1     7 FALSE FALSE     
#> 4 D11   N_CALIB        1     1     1     1     1     1     1 FALSE FALSE     
#>   
#> Inactive dimensions:
#>   
#>    dim   name       length   min   max unlim coord_dim 
#>    <chr> <chr>       <dbl> <dbl> <dbl> <lgl> <lgl>     
#>  1 D0    DATE_TIME      14     1    14 FALSE FALSE     
#>  2 D1    STRING2         2     1     2 FALSE FALSE     
#>  3 D2    STRING4         4     1     4 FALSE FALSE     
#>  4 D3    STRING8         8     1     8 FALSE FALSE     
#>  5 D4    STRING16       16    NA    NA FALSE FALSE     
#>  6 D5    STRING32       32     1    32 FALSE FALSE     
#>  7 D6    STRING64       64     1    64 FALSE FALSE     
#>  8 D10   N_LEVELS      493     1   493 FALSE FALSE     
#>  9 D12   N_HISTORY       0    NA    NA TRUE  FALSE     
#> 10 D13   N_VALUES41     41    NA    NA FALSE FALSE

## pass named expressions to subset dimension by value or index (step)
(subs <- tidync(filename) %>% hyper_filter(N_PROF = N_PROF > 1, STRING256 = index > 10))
#> Warning in hyper_filter.tidync(., N_PROF = N_PROF > 1, STRING256 = index > :
#> 'STRING256' not found in active grid, ignoring
#> 
#> Data Source (1): MD5903593_001.nc ...
#> 
#> Grids (16) <dimension family> : <associated variables> 
#> 
#> [1]   D0,D9,D11,D8 : SCIENTIFIC_CALIB_DATE
#> [2]   D6,D9,D11,D8 : PARAMETER
#> [3]   D7,D9,D11,D8 : SCIENTIFIC_CALIB_EQUATION, SCIENTIFIC_CALIB_COEFFICIENT, SCIENTIFIC_CALIB_COMMENT
#> [4]   D6,D9,D8     : STATION_PARAMETERS
#> [5]   D10,D8       : PRES, PRES_QC, PRES_ADJUSTED, PRES_ADJUSTED_QC, PRES_ADJUSTED_ERROR, TEMP, TEMP_QC, TEMP_ADJUSTED, TEMP_ADJUSTED_QC, TEMP_ADJUSTED_ERROR, PSAL, PSAL_QC, PSAL_ADJUSTED, PSAL_ADJUSTED_QC, PSAL_ADJUSTED_ERROR, DOXY, DOXY_QC, DOXY_ADJUSTED, DOXY_ADJUSTED_QC, DOXY_ADJUSTED_ERROR, CHLA, CHLA_QC, CHLA_ADJUSTED, CHLA_ADJUSTED_QC, CHLA_ADJUSTED_ERROR, BBP700, BBP700_QC, BBP700_ADJUSTED, BBP700_ADJUSTED_QC, BBP700_ADJUSTED_ERROR, NITRATE, NITRATE_QC, NITRATE_ADJUSTED, NITRATE_ADJUSTED_QC, NITRATE_ADJUSTED_ERROR    **ACTIVE GRID** ( 986  values per variable)
#> [6]   D1,D8        : DATA_CENTRE
#> [7]   D2,D8        : DATA_STATE_INDICATOR, WMO_INST_TYPE
#> [8]   D3,D8        : PLATFORM_NUMBER, POSITIONING_SYSTEM
#> [9]   D5,D8        : DC_REFERENCE, PLATFORM_TYPE, FLOAT_SERIAL_NO, FIRMWARE_VERSION
#> [10]   D6,D8        : PROJECT_NAME, PI_NAME
#> [11]   D7,D8        : VERTICAL_SAMPLING_SCHEME
#> [12]   D9,D8        : PARAMETER_DATA_MODE
#> [13]   D0           : REFERENCE_DATE_TIME, DATE_CREATION, DATE_UPDATE
#> [14]   D2           : FORMAT_VERSION, HANDBOOK_VERSION
#> [15]   D5           : DATA_TYPE
#> [16]   D8           : CYCLE_NUMBER, DIRECTION, DATA_MODE, JULD, JULD_QC, JULD_LOCATION, LATITUDE, LONGITUDE, POSITION_QC, CONFIG_MISSION_NUMBER, PROFILE_PRES_QC, PROFILE_TEMP_QC, PROFILE_PSAL_QC, PROFILE_DOXY_QC, PROFILE_CHLA_QC, PROFILE_BBP700_QC, PROFILE_NITRATE_QC
#> 
#> Dimensions 14 (2 active): 
#>   
#>   dim   name     length   min   max start count  dmin  dmax unlim coord_dim 
#>   <chr> <chr>     <dbl> <dbl> <dbl> <int> <int> <dbl> <dbl> <lgl> <lgl>     
#> 1 D8    N_PROF        2     1     2     2     1     2     2 FALSE FALSE     
#> 2 D10   N_LEVELS    493     1   493     1   493     1   493 FALSE FALSE     
#>   
#> Inactive dimensions:
#>   
#>    dim   name       length   min   max unlim coord_dim 
#>    <chr> <chr>       <dbl> <dbl> <dbl> <lgl> <lgl>     
#>  1 D0    DATE_TIME      14     1    14 FALSE FALSE     
#>  2 D1    STRING2         2     1     2 FALSE FALSE     
#>  3 D2    STRING4         4     1     4 FALSE FALSE     
#>  4 D3    STRING8         8     1     8 FALSE FALSE     
#>  5 D4    STRING16       16    NA    NA FALSE FALSE     
#>  6 D5    STRING32       32     1    32 FALSE FALSE     
#>  7 D6    STRING64       64     1    64 FALSE FALSE     
#>  8 D7    STRING256     256     1   256 FALSE FALSE     
#>  9 D9    N_PARAM         7     1     7 FALSE FALSE     
#> 10 D11   N_CALIB         1     1     1 FALSE FALSE     
#> 11 D12   N_HISTORY       0    NA    NA TRUE  FALSE     
#> 12 D13   N_VALUES41     41    NA    NA FALSE FALSE

## with the saved filtering from above, choose data frame or tbl_cube output
## optionally with only selected variables
subs %>% hyper_tibble()
#> # A tibble: 493 x 37
#>     PRES PRES_QC PRES_ADJUSTED PRES_ADJUSTED_QC PRES_ADJUSTED_E…  TEMP TEMP_QC
#>    <dbl> <chr>           <dbl> <chr>                       <dbl> <dbl> <chr>  
#>  1  7.70 1                7.79 1                            2.40  13.2 1      
#>  2 11.8  1               11.9  1                            2.40  13.2 1      
#>  3 16.3  1               16.4  1                            2.40  13.2 1      
#>  4 21.6  1               21.7  1                            2.40  13.2 1      
#>  5 26.7  1               26.8  1                            2.40  13.2 1      
#>  6 31.7  1               31.8  1                            2.40  13.2 1      
#>  7 36.6  1               36.7  1                            2.40  13.2 1      
#>  8 41.4  1               41.5  1                            2.40  13.2 1      
#>  9 46.5  1               46.6  1                            2.40  13.2 1      
#> 10 51.8  1               51.9  1                            2.40  13.2 1      
#> # … with 483 more rows, and 30 more variables: TEMP_ADJUSTED <dbl>,
#> #   TEMP_ADJUSTED_QC <chr>, TEMP_ADJUSTED_ERROR <dbl>, PSAL <dbl>,
#> #   PSAL_QC <chr>, PSAL_ADJUSTED <dbl>, PSAL_ADJUSTED_QC <chr>,
#> #   PSAL_ADJUSTED_ERROR <dbl>, DOXY <dbl>, DOXY_QC <chr>, DOXY_ADJUSTED <dbl>,
#> #   DOXY_ADJUSTED_QC <chr>, DOXY_ADJUSTED_ERROR <dbl>, CHLA <dbl>,
#> #   CHLA_QC <chr>, CHLA_ADJUSTED <dbl>, CHLA_ADJUSTED_QC <chr>,
#> #   CHLA_ADJUSTED_ERROR <dbl>, BBP700 <dbl>, BBP700_QC <chr>,
#> #   BBP700_ADJUSTED <dbl>, BBP700_ADJUSTED_QC <chr>,
#> #   BBP700_ADJUSTED_ERROR <dbl>, NITRATE <dbl>, NITRATE_QC <chr>,
#> #   NITRATE_ADJUSTED <dbl>, NITRATE_ADJUSTED_QC <chr>,
#> #   NITRATE_ADJUSTED_ERROR <dbl>, N_LEVELS <int>, N_PROF <int>
subs %>% hyper_tbl_cube(select_var = c("PRES", "PRES_QC", "PSAL_ADJUSTED"))
#> $mets
#> Class: tidync_data (list of tidync data arrays)
#> Variables (3): 'PRES', 'PRES_QC', 'PSAL_ADJUSTED'
#> Dimension (1): N_LEVELS,N_PROF (493)
#> Source: /perm_storage/home/mdsumner/R/x86_64-pc-linux-gnu-library/3.6/tidync/extdata/argo/MD5903593_001.nc
#> 
#> $dims
#> $dims$N_LEVELS
#>   [1]   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
#>  [19]  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36
#>  [37]  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54
#>  [55]  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72
#>  [73]  73  74  75  76  77  78  79  80  81  82  83  84  85  86  87  88  89  90
#>  [91]  91  92  93  94  95  96  97  98  99 100 101 102 103 104 105 106 107 108
#> [109] 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126
#> [127] 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144
#> [145] 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162
#> [163] 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180
#> [181] 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198
#> [199] 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216
#> [217] 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234
#> [235] 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252
#> [253] 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270
#> [271] 271 272 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288
#> [289] 289 290 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 306
#> [307] 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324
#> [325] 325 326 327 328 329 330 331 332 333 334 335 336 337 338 339 340 341 342
#> [343] 343 344 345 346 347 348 349 350 351 352 353 354 355 356 357 358 359 360
#> [361] 361 362 363 364 365 366 367 368 369 370 371 372 373 374 375 376 377 378
#> [379] 379 380 381 382 383 384 385 386 387 388 389 390 391 392 393 394 395 396
#> [397] 397 398 399 400 401 402 403 404 405 406 407 408 409 410 411 412 413 414
#> [415] 415 416 417 418 419 420 421 422 423 424 425 426 427 428 429 430 431 432
#> [433] 433 434 435 436 437 438 439 440 441 442 443 444 445 446 447 448 449 450
#> [451] 451 452 453 454 455 456 457 458 459 460 461 462 463 464 465 466 467 468
#> [469] 469 470 471 472 473 474 475 476 477 478 479 480 481 482 483 484 485 486
#> [487] 487 488 489 490 491 492 493
#> 
#> $dims$N_PROF
#> [1] 2
#> 
#> 
#> attr(,"class")
#> [1] "tbl_cube"
```

A grid is a “virtual table” in the sense of a database source. It’s
possible to activate a grid via a variable within it, so all variables
are available by default. Grids have identifiers based on which
dimensions they are defined with, so use i.e. “D1,D0” and can otherwise
be activated by their count identifier (starting at 1). The “D0” is an
identifier, it matches the internal 0-based indexing and identity used
by NetCDF itself.

Please note that `hyper_filter()` expressions must be unique, unlike
with `dplyr::filter()` we cannot load multiple comparisons into one.

While dplyr filter can load up multiple comparisons:

``` r
df %>% dplyr::filter(longitude > 100, longitude < 150)
```

in hyper\_filter we must load them into one named expression.

``` r
tidync(filename) %>% hyper_filter(longitude = longitude > 100 & longitude < 150)
```

### Extractive

Use what we learned interactively to extract the data, either in data
frame or raw-array (hyper slice) form.

``` r
## we'll see a column for the variable activated, and whatever other 
## variables the grid has
tidync(filename) %>% activate("JULD") %>% 
  hyper_filter(N_PROF = N_PROF == 1) %>% 
  hyper_tibble()
#> # A tibble: 98 x 5
#>    SCIENTIFIC_CALIB_DATE DATE_TIME N_PARAM N_CALIB N_PROF
#>    <chr>                     <int>   <int>   <int>  <int>
#>  1 2                             1       1       1      1
#>  2 0                             2       1       1      1
#>  3 1                             3       1       1      1
#>  4 7                             4       1       1      1
#>  5 0                             5       1       1      1
#>  6 4                             6       1       1      1
#>  7 1                             7       1       1      1
#>  8 0                             8       1       1      1
#>  9 1                             9       1       1      1
#> 10 4                            10       1       1      1
#> # … with 88 more rows


## native array form, we'll see a (list of) R arrays with a dimension for 
## each seen by tidync(filename) %>% activate("JULD")
tidync(filename) %>% activate("JULD") %>% 
  hyper_filter(N_PROF = N_PROF == 1) %>% 
  hyper_array()
#> Class: tidync_data (list of tidync data arrays)
#> Variables (1): 'SCIENTIFIC_CALIB_DATE'
#> Dimension (4): DATE_TIME,N_PARAM,N_CALIB,N_PROF (14, 7, 1, 1)
#> Source: /perm_storage/home/mdsumner/R/x86_64-pc-linux-gnu-library/3.6/tidync/extdata/argo/MD5903593_001.nc
```

It’s important to not actual request the data extraction until the
expressions above would result in an efficient size (don’t try a data
frame version of a 20Gb ROMs variable …). Use the interactive modes to
determine the likely size of the output you will receive.

Functions seamlessly build the actual index values required by the
NetCDF library. This can be used to debug the process or to define your
own tools for the extraction. Currently each `hyper_*` function can take
the filtering expressions, but it’s not obvious if this is a good idea
or not.

See the vignettes for more:

``` r
browseVignettes(package = "tidync")
```

## Limitations

Please get in touch if you have specific workflows that `tidync` is not
providing. There’s a lot of room for improvement\!

  - we can’t do “grouped filters”" (i.e. polygon-overlay extraction),
    but it’s in the works
  - compound types are not supported, though see the “rhdf5” branch on
    Github
  - NetCDF groups are not exposed (groups are like a “files within a
    file”, analogous to a file system directory)

I’m interested in lighter and rawer access to the NetCDF library, I’ve
explored that here and it may or may not be a good idea:

<https://github.com/hypertidy/ncapi>

## Terminology

  - **slab**, **hyperslab** - array variable that may be read from a
    NetCDF
  - **shape**, **grid** - set of dimensions that define variables in
    NetCDF
  - **activation** - choice of a given grid to apply subsetting and read
    operations to

-----

## Code of conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/ropensci/tidync/blob/master/CONDUCT.md). By
participating in this project you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# tidync 0.2.5

* Improved documentation for large data check in `hyper_array()` and `hyper_tibble()` and `hyper_tbl_cube()`, 
and allowed user-controlled option to avoid this check. Thanks to Alessandro Bigi for feedback. 

# tidync 0.2.4

* Now inports cubelyr rather than dplyr for `tbl_cube()` thanks to Hadley Wickham PR #102. 

# tidync 0.2.3

* Remove test burden. 

# tidync 0.2.2

* Fix broken tests due to changes in RNetCDF, thanks to CRAN for notification. 

* Added example to readme and docs, and checks for non-unique `hyper_filter()` names, as per #93. 
 Thanks to @everydayduffy for the suggestion. 

# tidync 0.2.1

* Now depends on R 3.5.0. 

* Removed obscure numeric 'what' logic for `activate()`. 

* Fixed tidync method for `hyper_array()` output. 

* Found a really bad bug in `hyper_array()` (#92), now fixed. Axis order in 
  transforms was sometimes reversed, which caused garbled results from `hyper_tibble()`. 
  The effect is to ruin any ggplot2 figures for some source files. It's very likely
  that no other system yet uses the `hyper_array()` format so impact is low.  

# tidync 0.2.0

* tidync is now part of the rOpenSci project. 

* Fixes in tests and examples to avoid version-incapable NetCDF problems on Solaris. 
 
* A number of improved tests and documentation fixes. 

* Deal with warnings from tidyr version > 0.8.3. 

* A huge thank you to all rOpenSci contributions, especially very helpful reviews from   
  [&#x0040;Nowosad](https://github.com/Nowosad) and [&#x0040;timcdlucas](https://github.com/timcdlucas) as well as organizers [&#x0040;karthik](https://github.com/karthik), [&#x0040;stefaniebutland](https://github.com/stefaniebutland)  and [&#x0040;sckott](https://github.com/sckott). 
  
  Helpful input was also provided via issues  from [&#x0040;adrfantini](https://github.com/adrfantini), [&#x0040;JustBerkhout](https://github.com/JustBerkhout), [&#x0040;matteodefelice](https://github.com/matteodefelice), [&#x0040;rensa](https://github.com/rensa), 
  [&#x0040;rmendels](https://github.com/rmendels), [&#x0040;sw-rifai](https://github.com/sw-rifai), and [&#x0040;tremenyi](https://github.com/tremenyi). 
  
  Pull requests were contributed by [&#x0040;dlebauer](https://github.com/dlebauer), [&#x0040;edzer](https://github.com/edzer). 

# tidync 0.1.1

* Package improvements thanks to CRAN feedback, clarified Description and added
  more examples. Replaced `cat()` and `print()` calls with `message()` and `warning()`. 

* New class 'tidync_data' for output of `hyper_array()`, no underlying change to 
  the object which is simply a list of arrays from each variable, and axis transforms 
  stored in an attribute. 

* Old deprecated function `axis_transforms()` now Defunct. 


# tidync 0.1.0

* FIRST RELEASE, tidync was greatly improved via help from the rOpenSci review 
  process. 

* New function `hyper_grids()` to report available grid names. 

* A printing error of dimension value ranges is now fixed, thanks to James Goldie 
  (#84). 

* Now supports 'NC_CHAR' type, by exploding these into the array size expected. 

* Breaking change: when using `tidync$grid`it's now expected that this must be  
  `tidyr::unnest()`ed in order to expand out the grid list per variable, in line 
  with https://github.com/hypertidy/ncmeta/issues/26. 
  
* The `hyper_array` function now stores the relevant transforms table as an attribute 
  `transforms` so that objects can be constructed directly from the native array 
  output. 

# tidync 0.0.3

* Function rename `hyper_array()` now matches `hyper_tibble()` indicating the form 
  of the output (rather than the action used, was `hyper_slice()`). 

* New functions `hyper_vars()` and `hyper_dims()` for reporting on the currently active 
  variables and dimensions and their status. 

* Now dependent on ncmeta >= 0.0.2, partly to avoid crashing on invalid source 
  or file strings

* Removed `hyper_index()` and incorporated that into `hyper_filter()`, there's 
  now only one delay-capable class which is 'tidync'


# tidync 0.0.2

* Function `hyper_filter()` now uses a selection idiom, to record the state of 
  the axis rather than explicitly filter it. This means we can have more flexibility 
  on what the axis transform tables can be used for, and removes some unwieldy 
  handling code. All the available axes are on the object from first contact, which 
  means we can program against the entire space in the source which will help for 
  complex mapping scenarios. 
 
* Function `hyper filter()` print now handles the case of char-type coordinate 
  values by setting the min and max to NA_real_

* Various improvements and fixes for the print method for tidync

* Support coordinate-less dimensions has been added, there is new information in 
  the print summary about which dimensions are a "coord_dim" and this results in 
  the axis transform tables using the index as the coordinate value.  

* Functions `hyper_slice()` and `hyper_tibble()` now return all variables that 
  exist within a grid

* This version sees a new model where activation is on 'grids', effectively a 
  space composed of dimensions. In addition to the variables, dimension, and attributes 
  entities we add 'grid' defined by a set of dimensions, and 'axis' which is an 
  instance of a particular dimension as used by a variable. 

* Sources without recognizable variables now gracefully handled, with help from 
  ncmeta. 

# tidync 0.0.1

* Now imports ncdump > 0.0.3. 

* Installed external example data from [Unidata website](https://www.unidata.ucar.edu/software/netcdf/examples/files.html)

* First working version now has `tidync()`, and `hyper_*()` family of functions. 

* Migrated from ncdump. 


## tidync 0.2.4

Hello, update to link to cubelyr rather than dplyr, for dplyr 1.0.0. 

Thank you. 

## Test environments
* local Ubuntu 18.04 install, R 3.6.3
* win-builder (release and devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

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

Please note that the tidync project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://ropensci.github.io/dev_guide/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if

* you have a question, an use case, or otherwise not a bug or feature request for the software itself.
* you think your issue requires a longer form discussion.

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
---
output: github_document
editor_options: 
  chunk_output_type: console
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

# tidync

<!-- badges: start -->
[![Travis-CI Build Status](http://badges.herokuapp.com/travis/ropensci/tidync?branch=master&env=BUILD_NAME=trusty_release&label=linux)](https://travis-ci.org/ropensci/tidync)
[![Build Status](http://badges.herokuapp.com/travis/ropensci/tidync?branch=master&env=BUILD_NAME=osx_release&label=osx)](https://travis-ci.org/ropensci/tidync)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/ropensci/tidync?branch=master&svg=true)](https://ci.appveyor.com/project/mdsumner/tidync-nvwwt)
[![Coverage status](https://codecov.io/gh/ropensci/tidync/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/tidync?branch=master)
[![](https://badges.ropensci.org/174_status.svg)](https://github.com/ropensci/onboarding/issues/174)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing) [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/tidync)](https://cran.r-project.org/package=tidync)
[![CRAN_Download_Badge](http://cranlogs.r-pkg.org/badges/tidync)](https://cran.r-project.org/package=tidync)
  <!-- badges: end -->


The goal of tidync is to ease exploring the contents of a NetCDF source and to simplify the process of data extraction. 

When extracting, data can be accessed as array/s, or in long-form as a data frame.
In contrast to other packages tidync helps reduce the volume of code required to discover and read the contents of NetCDF, with simple steps:

* Connect and summarize `tidync()`.
* (optionally) Specify source variables `activate()`.
* (optionally) Specify array sub-setting (slicing) `hyper_filter()`.
* Read array data in native form `hyper_array()` or long-form `hyper_tibble()` or bespoke form `hyper_tbl_cube()`. 

NetCDF is **Network Common Data Form** a very common, and very general way to store and work with 
scientific array-based data. NetCDF is defined and provided by [Unidata](https://www.unidata.ucar.edu/software/netcdf/). R has (independent) support for NetCDF via the [ncdf4](https://CRAN.R-project.org/package=ncdf4), [rhdf5](https://bioconductor.org/packages/release/bioc/html/rhdf5.html), [RNetCDF](https://CRAN.R-project.org/package=RNetCDF), [rgdal](https://CRAN.R-project.org/package=rgdal), [sf](https://CRAN.R-project.org/package=sf) and [vapour](https://CRAN.R-project.org/package=vapour) packages. 

This project uses RNetCDF for the primary access to the NetCDF library, as well as the ncdf4 package in some cases. The wrapper provided by [ncmeta](https://CRAN.R-project.org/package=ncmeta) over RNetCDF is used to obtain information about data sources. 

## Installation

Install tidync from CRAN. 

```{r cran, eval = FALSE}
install.packages("tidync")
```

You can install the development version from github with:

```{r gh-installation, eval = FALSE}
# install.packages("remotes")
remotes::install_github("ropensci/tidync", dependencies = TRUE)
```

The package packages ncdf4 and RNetCDF are required, so first make sure you can install and use these if it doesn't work the first time. 

```{r,eval=FALSE}
install.packages("ncdf4")
install.packages("RNetCDF")
```

If you have problems, please see the [INSTALL instructions for RNetCDF](https://CRAN.R-project.org/package=RNetCDF/INSTALL), these should work as well for ncdf4. Below I note specifics for different operating systems, notably Ubuntu/Debian where I work the most - these aren't comprehensive details but might be helpful.

### Windows

On Windows, everything should be easy as ncdf4 and RNetCDF are supported by CRAN. The RNetCDF package now includes OpenDAP/Thredds for 64-bit Windows (not 32-bit), and so tidync will work for those sources too. 

### MacOS

On MacOS, it should also be easy as there are binaries for ncdf4 and RNetCDF available on CRAN. As far as I know, only RNetCDF will support Thredds. 

### Ubuntu/Debian

On Linux you will need at least the following installed by an administrator,
here tested on Ubuntu Xenial 16.04. 

```bash
apt update 
apt upgrade --assume-yes

## Install 3rd parties for NetCDF
apt install libnetcdf-dev libudunits2-dev

## install 3rd parties needed for devtools + openssl git2r httr
apt install libssl-dev
```

Then in R

```R
install.packages("remotes")
remotes::install_github("ropensci/tidync")
```



More general information about system dependencies libnetcdf-dev and libudunits2-dev is available from [Unidata](https://www.unidata.ucar.edu). 

## Usage

This is a basic example which shows how to connect to a file. 

```{r start}
file <- system.file("extdata", "oceandata", "S20080012008031.L3m_MO_CHL_chlor_a_9km.nc", package = "tidync")
library(tidync)
tidync(file) 
```

There are two main ways of using tidync, interactively to explore what is there, and for extraction. The functions
`tidync` and `activate` and `hyper_filter` allow us to hone in on the part/s of the data we want, and functions `hyper_array`, `hyper_tibble` and `hyper_tbl_cube` give raw-array or data frames. 


### Interactive

Use `tidync()` and `hyper_filter()` to discern what variables and dimensions are available, and to craft axis-filtering
expressions by value or by index. (Use the name of the variable on the LHS to target it, use its name to filter by value and the special name `index` to filter it by its index). 

```{r interactive,eval=TRUE}
filename <- system.file("extdata/argo/MD5903593_001.nc", package = "tidync")
## discover the available entities, and the active grid's dimensions and variables
tidync(filename)

## activate a different grid
grid_identifier <- "D7,D9,D11,D8"
tidync(filename) %>% activate(grid_identifier)

## pass named expressions to subset dimension by value or index (step)
(subs <- tidync(filename) %>% hyper_filter(N_PROF = N_PROF > 1, STRING256 = index > 10))

## with the saved filtering from above, choose data frame or tbl_cube output
## optionally with only selected variables
subs %>% hyper_tibble()
subs %>% hyper_tbl_cube(select_var = c("PRES", "PRES_QC", "PSAL_ADJUSTED"))
```


A grid is a "virtual table" in the sense of a database source. It's possible to activate a grid via a variable within it, so all variables are available by default. Grids have identifiers based on which dimensions they are defined with, so use i.e. "D1,D0" and can otherwise be activated by their count identifier (starting at 1). The "D0" is an identifier, it matches the internal 0-based indexing and identity used by NetCDF itself. 

Please note that `hyper_filter()` expressions must be unique, unlike with `dplyr::filter()` we cannot load multiple comparisons into one. 

While dplyr filter can load up multiple comparisons: 

```{r pseudo1, eval=FALSE}
df %>% dplyr::filter(longitude > 100, longitude < 150)
```

in hyper_filter we must load them into one named expression.  

```{r pseudo2, eval=FALSE}
tidync(filename) %>% hyper_filter(longitude = longitude > 100 & longitude < 150)
```


### Extractive

Use what we learned interactively to extract the data, either in data frame or raw-array (hyper slice) form. 

```{r extraction, eval=TRUE}
## we'll see a column for the variable activated, and whatever other 
## variables the grid has
tidync(filename) %>% activate("JULD") %>% 
  hyper_filter(N_PROF = N_PROF == 1) %>% 
  hyper_tibble()


## native array form, we'll see a (list of) R arrays with a dimension for 
## each seen by tidync(filename) %>% activate("JULD")
tidync(filename) %>% activate("JULD") %>% 
  hyper_filter(N_PROF = N_PROF == 1) %>% 
  hyper_array()

```

It's important to not actual request the data extraction until the expressions above would result in an efficient size (don't try a data frame version of a 20Gb ROMs variable ...). Use the interactive modes to determine the likely size of the output you will receive. 

Functions seamlessly build the actual index values required by the NetCDF library. This can be used to debug the process or to define your own tools for the extraction. Currently each `hyper_*` function can take the filtering expressions, but it's not obvious if this is a good idea or not. 


See the vignettes for more: 

```R
browseVignettes(package = "tidync")
```

## Limitations

Please get in touch if you have specific workflows that `tidync` is not providing. There's a lot of room for improvement!

* we can't do "grouped filters"" (i.e. polygon-overlay extraction), but it's in the works
* compound types are not supported, though see the "rhdf5" branch on Github
* NetCDF groups are not exposed (groups are like a "files within a file", analogous to a file system directory)

I'm interested in lighter and rawer access to the NetCDF library, I've explored that here and it may or may not be a good idea:

https://github.com/hypertidy/ncapi

## Terminology

* **slab**, **hyperslab** - array variable that may be read from a NetCDF 
* **shape**, **grid** - set of dimensions that define variables in NetCDF
* **activation** - choice of a given grid to apply subsetting and read operations to


---

## Code of conduct

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/ropensci/tidync/blob/master/CONDUCT.md). By participating in this project you agree to abide by its terms.



[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
<!--
We must produce these tables on demand in whatever most efficient process is required. There are also types of forms that do not need to be stored in a single table, drawing a matrix of pixel values on screen being one clear example. When it comes to *scientific array data* in R, if we ignore underlying patterns in *structured-observations* such as coordinates or use the long-form single table form naively then it is very inefficient. 

The first significant tidyverse package was ggplot2 and this had an enormouse influence on graphics in R, but from the beginning required an inefficient intermediate form for both raster and vector spatial graphics. This was pilloried by powerful voices in R and the community has never bridged the divide, with most users repelled by the toxic communications on display. 

Similar tensions exist between imagery used in a GIS (Geographical Information System) context and scientific array data used in environmental observations and physical modelling. In general, GIS contexts use an *affine transformation* and modelling contexts will use a rectilinear or curvilinear transform. There's no clearer illustration of this divide than that the *NetCDF format* does not support the affine transformation method, even when the use of rectilinear or curvilinear coordinates is degenerate or redundant [^1]. 

# Array spaces

Array data live inside a *space* and if that is more than the rows and columns (and slices ...) of the array then there's three ways to specify them. 

* affine transform
* rectilinear transform
* curvilinear transform. 

Consider a simple matrix of 12 values. 

```{r matrix}
m <- matrix(1:12, nrow = 3L, ncol 4L)
```


Let's assume that these define data points that exist in the following intervals. 

```{r intervals}
xi <- c(1, 4)
yi <- c(1, 3)
```

We have everything needed to plot these data as if it were a spatial image. 

```{r image}
image(x = seq(xi[1L], xi[2L], length.out = nrow(m) + 1L), 
      y = seq(yi[1L], yi[2L]), length.out = ncol(m) + 1), z = m)
```
If we store data in *long form* as per the tidyverse, it really expects us to not only store every variable but also *every coordinate* of the array elements. This is what leads to (sometimes angry) complaints about the efficiency and scalability of parts of the tidyverse for scientific array data. 

There are several approaches to this, and tidync tries to insert itself into the middle. 

*NB: This post is not a comprehensive introduction to NetCDF, but does explore some of the problems encountered in some detail. The hope is that tidync will provide easier ways to explore your data, and allow array extraction or expansion into a data frame that is easier than coding from scratch.* 


[^1]: NetCDF *can be used to store an affine transformation*. It can be used to store anything as you will often hear - it's just not terribly well suited to some tasks, and this usage is rare and usually bumps against the usual assumptions made by users of the format.  The GMT variant of NetCDF does use a form of affine transformation, but also stores arrays in long form 1D vectors with the dimensional metadata stored seaprately. 

-->---
slug: tidync
title: 'tidync: scientific array data from NetCDF in R'
package_version: 0.2.2
authors:
  - Michael Sumner
date: 2019-11-05
categories: blog
topicid:
tags:
- Software Peer Review
- R
- community
- software
- packages
- tidync
- NetCDF
- array
- tidyverse
- data
twitterImg: img/blog-images/2019-11-05-tidync/oisst-data-single-line-1.png
output: 
  html_document:
    keep_md: true
---

```{r setup, include=FALSE}
## build the blog, Knit, then 
## remove ../../../roweb2/themes/ropensci/static from the .md
## remove local path from raster
knitr::opts_chunk$set(echo = TRUE, fig.path = "../../../roweb2/themes/ropensci/static/img/blog-images/2019-11-05-tidync/")
knitr::opts_chunk$set(root.dir=normalizePath(here::here("../../themes/ropensci/static")))
library(tidync)
library(dplyr)
library(tidyr)
library(ncmeta)
library(stars)
library(raster)
```

In May 2019 version 0.2.0 of [tidync](https://docs.ropensci.org/tidync/) was approved by rOpenSci and accepted to CRAN. Here we provide a [quick overview](#overview) of the typical workflow with some pseudo-code for the main functions in tidync. This overview is enough to read if you just want to try out the package on your own data.  The tidync package is focussed on *efficient data extraction* for developing your own software, and this somewhat long post takes the time to explain the concepts in detail.  

There is a section about the [NetCDF](#netcdf)  data model itself. Then there is  a [detailed illustration](#raster-data-in-netcdf) of a raster data set in R including some of the challenges faced by R users. This is followed by sections on how tidync sees [metadata](#metadata) and [coordinates](#transforms) in NetCDF, how we can [slice a dataset](#slicing) and control the format of the output. We then discuss [some limitations](#limitations) and [future work](#future) , and then (most importantly) reflect on the [rOpenSci process](#review) of package review. 

<br>

### NetCDF in R {#netcdf-in-R}

[NetCDF](https://www.unidata.ucar.edu/software/netcdf/) is a very widely used system for storing and distributing scientific array data. A NetCDF data source typically stores one or more arrays of data, along with metadata that describe the data array space (grid), and any metadata describing array coordinates, units, and interpretation. A NetCDF source may be a file or an online URL. If you want to automate your own workflow around a series of NetCDF data sources then [tidync](https://docs.ropensci.org/tidync/) provides all the flexibility and power required with as little pain as possible.

The [tidyverse](https://www.tidyverse.org/) has had an enormous impact on the use of R with a strict approach to [*variables* and *observations*](https://r4ds.had.co.nz/tidy-data.html) (in short, tidy data are tabular, with each variable having its own column and each observation having its own row). This tidy-data-frame form can be used for a wide range of data, but it does have some shortcomings. It can be inefficient in terms of storage, which may be problematic with large data. If the data are accompanied by additional metadata (as is generally the case with NetCDF data) there is often no neat way to store this information in the same table, and these inherent properties of the original data can be lost. 

There is a tension between the **tidyverse** and **scientific array data** that comes down to this difference in data storage, and the intermediate forms used to get between one form and another. We also think this tension has been exaggerated in unproductive ways. 

<br>

### tidync {#overview}

The official website for tidync is https://docs.ropensci.org/tidync/ and the latest release can be found on [CRAN](https://CRAN.R-project.org/package=tidync). 

The tidync package provides a compromise position, allowing efficient interaction with NetCDF files.  It will produce native-array *or* tidy-data-frame output as desired. It delays any data-reading activity until after the output format is chosen. In particular, tidync exists in order to reduce the amount of plumbing code required to get to the data. It allows an interactive way to convert between the different spaces (coordinates and indices) in which the data can be referenced.  

In pseudo-code, there are only a few simple steps, at each step we can save the result and explore a summary. 

1. Connect to a data source and retrieve metadata, and read a summary: 

```R
src <- tidync(<netcdf-source>)
print(src)
```

2. By default the largest array-space (grid) is *activated* and usually this will be the right choice - if required we can nominate a different grid using `activate()`. 

```R
src <- src %>% activate(<a different grid>)
```

3. Apply subsetting to *slice arrays* by coordinate or index, this step is optional but very important for large and complicated data sources. 

```R
## lazy subsetting by value or index
src_slc <- src %>% hyper_filter(<filter expressions on dimensions>)
```

4. Finally, choose an output format - list of arrays, a data frame, or a `tbl_cube`. 

```R
src_slc %>% hyper_array()

src_slc %>% hyper_tibble()

src_slc %>% hyper_tbl_cube()
```

There are various other packages for NetCDF in R, the main ones being [RNetCDF](https://CRAN.r-project.org/package=RNetCDF) and [ncdf4](https://CRAN.r-project.org/package=ncdf4). These are both *lower-level* tools than tidync - they are interfaces to the [underlying NetCDF library](https://github.com/Unidata/netcdf-c), and tidync uses both to read information and data. The [raster](https://CRAN.r-project.org/package=raster) and [stars](https://CRAN.r-project.org/package=stars) packages provide quite different approaches and stars is more general than raster, but is similarly *higher-level* than tidync. 

To follow along with the code below requires all of the following packages, and we assume that recent versions are in use, particularly `ncmeta (>= 0.2.0)`, 
`tidync (>= 0.2.2)`, and `tidyr (>= 1.0.0)`. 

```R
install.packages(c("ncmeta", "tidync", "maps", "stars", "ggplot2", "devtools", 
                   "stars", "RNetCDF", "raster", "dplyr", "tidyr"))
```


<br>

### NetCDF {#netcdf}

NetCDF is a very widely used file format for storing array-based data as
**variables** with **dimensions** and **attributes**. 

The *space* (or ***grid***) occupied by a **variable** is defined by its **dimensions** and their **attributes**. Dimensions are by definition
*one-dimensional arrays* (i.e. an atomic vector in R of length 1 or more). An array can include coordinate metadata, units, type, and interpretation; **attributes** define all of this extra information for dimensions and variables.  The *space* of a variable (the ***grid*** it lives in) is defined by one or more of the dimensions in the file. 

A given variable won't necessarily use all the available dimensions and no dimensions are mandatory or particularly special. We consider the existence of a dimension within a grid to be an *instance of that dimension* and call that an ***axis***, subtly different to the dimension on its own.

NetCDF is very general and used in many different ways. It is quite common to see subcultures that rally around the way their particular domain's data are used and stored while ignoring many other valid ways of using NetCDF. The tidync approach is to be as general as possible, sacrificing high level interpretations for lower-level control if that generality is at risk. 

<br>

### Raster data in NetCDF {#raster-data-in-netcdf}

NetCDF can be used to store *raster data*, and very commonly data are provided as a global grid of scientific data. Here we use a snapshot of global ocean surface temperature for a single day. The file used is called `reduced.nc` in the stars package, derived from the daily [OISSTV2 product](https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html). 

We will explore this data set in detail to give an introduction to the tidync summary and functions. The data set is very commonly used in marine research as it includes a very long time series of daily global maps of sea surface temperatures. The spatial resolution is 0.25 degree (1/4°) and the coverage is complete for the entire world's oceans because it is a blend of remote sensing and direct observations. We call it **OISST** (or *oisst*) after its [official description](https://www.ncdc.noaa.gov/oisst). 

> The NOAA 1/4° daily *Optimum Interpolation Sea Surface Temperature* (or daily OISST) is an analysis constructed by combining observations from different platforms (satellites, ships, buoys) on a regular global grid. A spatially complete SST map is produced by interpolating to fill in gaps. 

The file we use is a simplified version from 1981-12-31 (the series started in 1981-09-01). This has been reduced in resolution so that it can be stored in an R package. There are four variables in the data **sst** (sea surface temperature, in Celsius), **anom** (anomaly of sst of this day from the long term mean), **err** (estimated error standard deviation of sst), and **ice** (sea ice concentration). The ice concentration acts as a mask, if there is sea ice present in the pixel then the temperature is undefined. 

```{r basic-raster}
oisstfile <- system.file("nc/reduced.nc", package = "stars")
```

To connect to this file use `tidync()`. 

```{r basic-raster-connect}
library(tidync)
oisst <- tidync(oisstfile)
```

*Note: this is not a [file connection](https://stat.ethz.ch/R-manual/R-devel/library/base/html/connections.html), like that used by ncdf4 or RNetCDF; tidync functions always open the file in read-only mode, extract information and/or data, and then close the open file connection.*

To see the available data in the file print a summmary of the source. 

```{r basic-raster-print}
print(oisst)
```

There are three kinds of information 

*  one (1) Data Source, our one file
* five (5) Grids, available *spaces* in the source 
* four (4) Dimensions, orthogonal axes from which Grids are composed

There is only one grid available for multidimensional data in this file, the first one "D0,D1,D2,D3" - all other grids are one-dimensional. This 4D grid has four variables `sst`, `anom`, `ice`, and `err` and each 1D grid has a single variable. 

*Note: it's only a coincidence that this 4D grid also has 4 variables*. 

The 1D grids have a corresponding dimension `dim` and variable `name`, making these *coordinate dimensions* (see `coord_dim` in the dimensions table). It's not necessarily true that a 1D grid will have a single 1D variable, it may have more than one variable, and it may only have an *index variable*, i.e. only the position values `1:length(dimension)`. 

Each dimension's name, length, valid minimum and maximum values are seen in the Dimensions table and these values can never change, also see flags  `coord_dim` and `unlim`. This refers to an unlimited dimension, used when a data time series is spread across multiple files. 

The other Dimensions columns `start`, `count`, `dmin`, `dmax` apply when we slice into data variables with `hyper_filter()`. 

<br>

#### Metadata in tidync: ncmeta {#metadata}

NetCDF is such a general form for storing data that there are many ways to approach its use. We wanted to focus on a *tidy approach to NetCDF* and so built on existing packages to do the lower level tasks. tidync relies on the package [ncmeta](https://CRAN.r-project.org/package=ncmeta) to extract information about NetCDF sources. There are functions to find available variables, dimensions, attributes, grids, and axes in ncmeta. 

We also want functions that return information about each kind of entity in a straightforward way. Since there is a complex relationship between variables, dimensions and grids we cannot store this information well in a single structure.  

See that there are 5 grids and 8 variables, with a row for each. 

```{r ncmeta-grids}
ncmeta::nc_grids(oisstfile)

ncmeta::nc_vars(oisstfile)
```

Each grid has a name, dimensionality (`ndims`), and set of variables. Each grid is listed only once, an important pattern when we are programming, and the same applies to variables. The relationship to the tidyverse starts here with the metadata; there are five grids observed and we have four columns of information for each grid. These are the grid's name, number of dimensions, the NetCDF-variables defined with it, and their number. When dealing with metadata we can also use tidy principles as we do with the data itself. 

Some grids have more than one variable, so they are nested in the grid rows - use `tidyr::unnest()` to see all variables with their parent grid. 

```{r ncemeta-grids-expandvars}
ncmeta::nc_grids(oisstfile) %>% tidyr::unnest(cols = c(variables))
```

Similar functions exist for dimensions and variables. 

```{r ncmeta-dimensions-variables}
ncmeta::nc_dims(oisstfile)

ncmeta::nc_vars(oisstfile)
```

There are corresponding functions to find out more about individual variables, dimensions and attributes by name or by index. 


Note that we can use the *internal index* (a zero-based count) of a variable as well as its name (`anom` is
the variable at the 5-index  as shown by `nc_vars()` above). 

```{r ncmeta-variable1}
ncmeta::nc_var(oisstfile, "anom")
ncmeta::nc_var(oisstfile, 5)

```

Similarly we can use name or index for dimensions and attributes, but attributes for a variable can only be found by name. 

```{r ncmeta-dimension1}
ncmeta::nc_dim(oisstfile, "lon")
ncmeta::nc_dim(oisstfile, 0)

ncmeta::nc_atts(oisstfile)
ncmeta::nc_atts(oisstfile, "zlev")
```

We can find the internal metadata for each variable by expanding the value. 

```{r ncmeta-time-attributes}
ncmeta::nc_atts(oisstfile, "time") %>% tidyr::unnest(cols = c(value))
```

With this information we may now apply the right interpretation to the time values. In the print of the tidync object above we see the value `1460`, which is given without context in the dimensions table. 

We can get that value by activating the right grid and extracting, the `time` is
a single integer value. 

```{r activate-time}
oisst <- tidync(oisstfile)
time_ex <- oisst %>% activate("D3") %>% hyper_array()
time_ex$time
```

As mentioned, tidync considers time and its metadata a bit dangerous. A record about these can often be wrong, inconsistent, include time zone issues, sometimes *extra seconds* to account for ... and so we prefer to leave these interpretations to be validated manually,  before automation. 

Obtain the time units information and then use it to convert the raw value (`1460`) into a date-time understood by R. 

```{r ncmeta-time-atts}
tunit <- ncmeta::nc_atts(oisstfile, "time") %>% tidyr::unnest(cols = c(value)) %>% dplyr::filter(name == "units")

print(tunit)

time_parts <- RNetCDF::utcal.nc(tunit$value, time_ex$time)

## convert to date-time
ISOdatetime(time_parts[,"year"], 
            time_parts[,"month"], 
            time_parts[,"day"], 
            time_parts[,"hour"], 
            time_parts[,"minute"], 
            time_parts[,"second"])
```


Alternatively we can do this by hard-coding. We find that different cases are best handled in different ways, and especially after some careful checking. 

```{r ncmeta-time-conversion}
as.POSIXct("1978-01-01 00:00:00", tz = "UTC") + time_ex$time * 24 * 3600
```

Finally, we can check that other independent systems provide the same information. 

```{r raster-stars}
raster::brick(oisstfile, varname = "anom")
stars::read_stars(oisstfile)
```

In terms of *interpreting the meaning of stored metadata*, tidync shies away from doing this automatically. There are simply too many ways for automatic tools to get intentions wrong. So, used in combination the ncmeta and tidync packages provide the tools to program around the vagaries presented by NetCDF sources. If your data is pretty clean and standardized there is higher-level software that can easily interpret these things automatically. Some examples are [stars](https://CRAN.r-project.org/package=stars) the R package, and outside of R itself there is also  [xarray](http://xarray.pydata.org), [GDAL](https://gdal.org/), [ferret](https://ferret.pmel.noaa.gov/Ferret/), and [Panoply](https://www.giss.nasa.gov/tools/panoply/). 

<br>

#### Axes versus dimensions

Previously we mentioned the concept of an ***axis*** as an instance of a dimension. This distinction arose early on, sometimes the dimension on its own is important, at other times we want to know where it occurs. The functions `nc_axes()` and `nc_dims()` make this clear, every instance of a dimension across variables is listed as an axis, but they are derived from only four dimensions. 

```{r axis-vs-dimension}
ncmeta::nc_axes(oisstfile)

ncmeta::nc_dims(oisstfile)
```

<br>

#### Degenerate dimensions

See that both `zlev` and `time` are listed as dimensions but have length 1, and also their min and max values are constants. The `zlev` tells us that this grid exists at elevation = 0 (the sea surface) and `time` that the data applies to `time = 1460`. The time is not expressed as a duration even though it presumably applies to the entire day. These are *degenerate dimensions*, i.e. the data are really 2D but we have a record of a 4D space from which they are expressed as a slice. For time, we know that the neighbouring days exist in other OISST files, but for `zlev` it only records sea-level. This can cause problems as we would usually treat this data as a matrix in R, and so the ncdf4 and RNetCDF package read functions have arguments that are analogous to R's array indexing argument `drop = TRUE`. If a dimension of length 1 is encountered the 'to drop' means to ignore it. tidync will also drop dimensions by default when reading data, see the `drop` argument in `?hyper_array`. 

<br>

#### Reading the OISST data

At this point only metadata has been read, so let's read some sea surface temperatures!

The fastest way to get all the data is to call the function `hyper_array`, this is the lowest level and is very close to using the ncdf4 or RNetCDF package directly. 

```{r read-data}
(oisst_data <- oisst %>% hyper_array())
```

What happened there? We got a classed object `tidync_data`; this is a list with arrays. 

```{r oisst-data}
length(oisst_data)
names(oisst_data)
dim(oisst_data[[1]])
image(oisst_data[[1]])
```

This is exactly the data provided by `ncdf4::ncvar_get()` or `RNetCDF::var.get.nc()` but we can do it in a single line of code. Without tidync we must find the variable names and loop over them. We automatically get variables from the largest grid that is available, which was `activate()`-d by default. 

```{r oisst-data-single-line}
oisst_data <- tidync(oisstfile) %>% hyper_array()
op <- par(mfrow = n2mfrow(length(oisst_data)))
pals <- c("YlOrRd", "viridis", "Grays", "Blues")
for (i in seq_along(oisst_data)) {
  image(oisst_data[[i]], main = names(oisst_data)[i], col = hcl.colors(20, pals[i], rev = i ==1))
}
par(op)
```

<br>

### Transforms {#transforms}

In this context ***transform*** means the conversion between index and geographic coordinate for grid cells, and in this data set this means the longitude and latitude assigned to the centre of each cell.

We have done nothing with the spatial side of these data, ignoring the lon and lat values completely. 

```{r oisst-data-dims}
oisst_data

lapply(oisst_data, dim)
```

The print summary of the `oisst_data` object shows that it knows there are four variables and that they each have two dimensions (`zlev` and `time` were *dropped*). This is now stored as a list of native R arrays, but there is also the transforms attribute available with `hyper_transforms()`. 

The values on each transform table may be used directly. 

```{r oisst-data-transforms}
(trans <- attr(oisst_data, "transforms"))

image(trans$lon$lon, trans$lat$lat,  oisst_data[[1]])
maps::map("world2", add = TRUE)
```


In this case these *transforms* are somewhat redundant, there is a value stored for every step in `lon` and every step in `lat`. They are completely regular series whereas the usual approach in graphics is to store an *offset and scale* rather than each step's coordinate. Sometimes these coordinate values are not reducible this way and we would call them *rectilinear*, we would have to store the sequence of each 1D coordinate step. 

<br>

### Slicing {#slicing}

We can slice into these dimensions using a tidyverse approach. For example, to slice out only the data for the waters of the Pacific Ocean, we need a range in longitude and in latitude. 

<br>

#### Old style slicing

This section illustrates the old laborious way to access a subset of data from NetCDF, a subset shown in this plot. 

```{r slicing-long-lat}
lonrange <- c(144, 247)
latrange <- c(-46, 47)

image(trans$lon$lon, trans$lat$lat,  oisst_data[[1]])
rect(lonrange[1], latrange[1], lonrange[2], latrange[2])
```


It's common on the internet to see posts that explain how to drive the NetCDF library with *start* and *count* indices, to do that we need to compare our ranges with the *transforms* of each dimension. 

```{r start-count}
xs <- findInterval(lonrange, trans$lon$lon)
ys <- findInterval(latrange, trans$lat$lat)
print(xs)
print(ys)
start <- c(xs[1], ys[1])
count <- c(diff(xs), diff(ys))

print(start)
print(count)


```

The idea here is that `xs` and `ys` tell us the columns and rows of interest, based on our geographic input in longitude latitude values that we understand. 

Let's try to read with NetCDF.  Hmmm .... what goes wrong. 

```{r read-RNetCDF-fail}
con <- RNetCDF::open.nc(oisstfile)
try(sst_matrix <- RNetCDF::var.get.nc(con, "sst", start = start, count = count))
```

We have been bitten by thinking that this source data are 2D!  So we just add start and count of 1 for each extra dimension. (Consider that it could 3D, or 5D, and maybe with different dimension order; all of these things complicate the general case for these otherwise simple solutions). 

```{r read-RNetCDF-succeed}
start <- c(start, 1, 1)
count <- c(count, 1, 1)
sst_matrix <- RNetCDF::var.get.nc(con, "sst", start = start, count = count)

```

And we're good! Except, we now don't have the coordinates for the mapping. We have to slice the lon and lat values as well, but let's cut to the chase and go back to tidync. 

<br>

#### tidync style slicing

Rather than slice the arrays read into memory, we can *filter* the object that understands the source and it does *not do any data slicing at all*, but records slices *to be done in future*.  This is the lazy beauty of the tidyverse, applied to NetCDF. 

Here we use standard R inequality syntax for `lon` and `lat`.  *We don't have to specify the redundant slice into zlev or time*. 

```{r tidync-slice}
library(dplyr)
oisst_slice <- oisst %>% hyper_filter(lon = lon > lonrange[1] & lon <= lonrange[2], 
                       lat = lat > latrange[1] & lat <= latrange[2])

oisst_slice
```


The print summary has updated the `start` and `count` columns now to match our laboriously acquired versions above. 

The `dmin` and `dmax` (data-min, data-max) columns are also updated, reporting the coordinate value at the start and end of the slice we have specified. 

Now we can break the lazy chain and call for the data. 

```{r hyper-array}
oisst_slice_data <- oisst_slice %>% hyper_array()
trans <- attr(oisst_slice_data, "transforms")
```

One unfortunate issue here is that we cannot use the transforms directly, they *have* been updated by changing the value of the `selected` column from `TRUE` to `FALSE`. Then we have to be aware of using only the values that remain *selected* (i.e. not filtered out).  

First filter the lon and lat transforms based on the `selected` column. 

```{r hyper-array-slice}
lon <- trans$lon %>% dplyr::filter(selected)
lat <- trans$lat %>% dplyr::filter(selected)

image(lon$lon, lat$lat, oisst_slice_data[[1]])
maps::map("world2", add = TRUE)
```

We do have to do extra work with `hyper_array()` but it gives total control over what we get. 

It's much easier to use other output types. 


```{r tbl-cube}
tcube <- tidync(oisstfile) %>% 
  hyper_filter(lon = between(lon, lonrange[1], lonrange[2]), 
                       lat = lat > latrange[1] & lat <= latrange[2]) %>% 
  hyper_tbl_cube()

tcube
```

We can also read our slice in directly as a tibble data frame, and plot with `geom_raster()`. 

```{r geom_raster}
tdata <- tidync(oisstfile) %>% 
  hyper_filter(lon = between(lon, lonrange[1], lonrange[2]), 
                       lat = lat > latrange[1] & lat <= latrange[2]) %>% 
  hyper_tibble()

library(ggplot2)
ggplot(tdata, aes(lon, lat, fill = sst)) + geom_raster()
```

By default, all variables are available but we can limit with `select_var`. 


```{r select-var}
tidync(oisstfile) %>% 
  hyper_filter(lon = between(lon, lonrange[1], lonrange[2]), 
                       lat = lat > latrange[1] & lat <= latrange[2]) %>% 
  hyper_tibble(select_var = c("err", "ice"))

```

<br>

#### slicing into multidimensional time series

As a further example, now open a *time-series* NetCDF file. We apply a spatial subset on the `lon` and `lat` dimensions, convert to tidy data frame and plot the `tos` variable over time. 

```{r time-series}
tos <- tidync(system.file("nc/tos_O1_2001-2002.nc", package = "stars"))
library(dplyr)
stos <- tos %>% hyper_filter(lon = between(lon, 140, 220), 
                     lat = between(lat, -60, 0)) %>% hyper_tibble()

library(ggplot2)
ggplot(stos, aes(lon, lat, fill = tos)) + geom_raster() + facet_wrap(~time)
```

We can alternatively choose the middle value of longitude (it lies at index = 90) and plot the `tos` variable as a function of latitude over time. We can easily re-orient our approach to this data set and it works as well with more complicated multi-dimensional sources as well. 


```{r time-lat-series}
lon180 <- tos %>%  hyper_filter(lon = index == 90, 
                     lat = between(lat, -60, 0)) %>% hyper_tibble()
ggplot(lon180, aes(time, lat, fill = tos)) + geom_raster() 

```

<br>

### Limitations {#limitations}

There are some limitations, specific to the tidync R package that are unrelated to the capabilities of the latest NetCDF library. 

* No groups, a group can be specified by providing the group-within-a-source *as a source*. 
* No compound types. 
* No attribute metadata, coordinates of 1D axes are stored as *transform tables*, but coordinates of pairs (or higher sets) of axes are not explicitly linked to their array data.  
* Curvilinear coordinates are not automatically expanded, this is because they exist (usually) on a different grid to the active one. 
* Unknowns about what is supported on what platforms. This is surprisingly tricky and unstable, there are a lot of things that are possible on one operating system at a given time, but not on others. The situation changes fairly slowly but is always changing due to library versions and releases, package and tooling support on CRAN, and operating system details. 

If you have problems with a given source please get in touch ([open an issue on Github issues](https://github.com/ropensci/tidync/), [chat on twitter](https://twitter.com/mdsumner/)) so we can learn more about the overall landscape. 

<br>

### Future helpers  {#future}

<br>

#### Coordinate expansion

A feature being considered for an upcoming version is to expand out all available linked coordinates. This occurs when an array has a dimension but only stores its index. When a dimension stores values directly this is known as a *dim-coord*, and usually occurs for time values. One way to expand this out would be to include an `expand_coords` argument to `hyper_tibble()` and have it run the following code: 

```{r internal-expand}
#' Expand coordinates stored against dimensions
#'
#' @param x tidync object
#' @param ... ignored
#'
#' @return data frame of all variables and any linked-coordinates 
#' @noRd
#'
#' @examples

full_expand <- function(x, ...) {
  ad <- active(x)
  spl <- strsplit(ad, ",")[[1L]]
  out <- hyper_tibble(x)
  
  for (i in seq_along(spl)) {
    out <- dplyr::inner_join(out, activate(x, spl[i]) %>% hyper_tibble())
  } 
  out
}
```

It's not clear how consistently this fits in the wider variants found in the NetCDF world, so any feedback is welcome. 

A real world example is available in the `ncdfgeom` package. This package provides much more in terms of storing geometry within a NetCDF file, but here we only extract the lon, lat and station name that `hyper_tibble()` isn't seeing by default. 

```{r ncdfgeom-example}
huc <- system.file('extdata','example_huc_eta.nc', package = 'ncdfgeom')

full_expand(tidync(huc))

hyper_tibble(tidync(huc))
```

<br>

#### Tidy approaches to other data sources

This approach could be applied to other array-based data systems, such as the [ff package](https://CRAN.r-project.org/package=ff), the [matter package](https://CRAN.r-project.org/package=matter) GDAL [raster](https://gdal.org/tutorials/index.html#raster) or [multi-dimensional](https://gdal.org/tutorials/index.html#multidimensional-raster) data sources, and [HDF5](https://www.hdfgroup.org/solutions/hdf5/) or [GRIB](https://en.wikipedia.org/wiki/GRIB) sources. 

We have experimented with this for non-NetCDF formats, please get in touch ([open an issue on Github issues](https://github.com/ropensci/tidync/), [chat on twitter](https://twitter.com/mdsumner/)) if you are interested.  

The [stars project](https://github.com/r-spatial/stars/) takes another perspective on a tidy approach to scientific array data. It is very high-level and may be a drop-in solution for well-behaved data so it's recommended to try that as well. 

<br>

### rOpenSci package review {#review}

The `tidync` package made it to CRAN after a fairly long review process on [rOpenSci](https://github.com/ropensci/software-review/issues/174). The package itself was inspired by many years of experience and discussions with [Tom Remenyi](https://github.com/tremenyi/), [Simon Wotherspoon](https://github.com/SWotherspoon/), [Sophie Bestley](https://github.com/snowpeaSoho/), and [Ben Raymond](https://github.com/raymondben/). In early 2018 I really wasn't sure if it could be finished at all in a neat way and was a bit overwhelmed, but thanks to very helpful reviewers and also some key insights about [obscure types](https://github.com/ropensci/tidync/issues/75#issuecomment-468064627) it was done. The package benefitted greatly from the review feedback provided by [Jakub Nowosad](https://github.com/Nowosad) and [Tim Lucas](https://github.com/timcdlucas). I really appreciated the clarity provided by these reviews, it helped to finalize some design decisions on the naming of functions and their intended use. There are various aspects that I thought were obstacles in completing the project, and having reviews that did not share my concerns and also gave positive feedback and suggestions for more relevant changes was extremely helpful. 

Thanks also to [rOpenSci](https://ropensci.org/) community members for encouragement and support!



---
title: "NetCDF with tidync"
author: "Michael D. Sumner"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    fig_width: 10
    fig_height: 10
vignette: >
  %\VignetteIndexEntry{NetCDF with tidync}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(pillar.subtle = FALSE, pillar.sigfig = 4)
```

The goal of tidync is to ease exploring the contents of a NetCDF source and to simplify the process of data extraction.  

Array dimension expressions virtually *slice* (i.e. *index*) into data arrays with immediate feedback, and ease programming for automating data extraction. No data is read until explicitly requested, and can be pulled directly as an array, or in long-form as a data frame. 

NetCDF is **Network Common Data Form**, a [data model][not-a-format] and an
[API](https://en.wikipedia.org/wiki/Application_programming_interface).

This is a very common, and *very general* way to store and work with scientific
array-based data. NetCDF is defined and provided by
[Unidata](https://www.unidata.ucar.edu/software/netcdf/).

Here we introduce traditional concepts of NetCDF, and show examples with
built-in files and online sources to demonstrate tidync functionality.

# NetCDF

NetCDF is a very widely used file format for storing array-based data as
*variables*. The **space** occupied by a **variable** is defined by its
**dimensions** and their metadata. Dimensions are by definition
*one-dimensional* (i.e. an atomic vector in R of length 1 or more), an array with coordinate metadata, units,
type and interpretation. The **space** of a variable is defined as one or more
of the dimensions in the file. A given variable won't necessarily use all the
available dimensions and no dimensions are mandatory or particularly special.

Some conventions exist to define usage and minimal standards for metadata for particular file schemas, but these are many and varied, and not always adhered to. Tidync does not provide a data interpretation as it is intended for general use, including by tools that create data formats. 

The R community is not particularly strong with use of NetCDF, though it
is common and widely used it pales compared to use in general climate science work, and there the most used tool is the [CDO Climate Data Operators](https://code.mpimet.mpg.de/projects/cdo).  In R the most common
tools used are ncdf4 and raster (which uses ncdf4). 

Both the RNetCDF and ncdf4 packages provide a traditional summary format,
familiar to many NetCDF users as the output of the command line program
[`ncdump`](https://www.unidata.ucar.edu/software/netcdf/workshops/most-recent/utilities/Ncdump.html).


```{r}
ice_file <- system.file("extdata", "ifremer", "20171002.nc", package = "tidync", mustWork = TRUE)
library(RNetCDF)
print.nc(open.nc(ice_file))
```

Using `ncdump` at the command line on a suitable system would yield very similar
output to the print above.

```bash 
ncdump -h /path/to/extdata/ifremer/20171002.nc
```

With the ncdf4 package it's a slightly different approach, but gives the same
result.

```R
print(ncdf4::nc_open(ice_file))
```

Notice how the listing above is organized by *dimension* and then by *variable*.
It's not particularly obvious that some variables are defined within the same
set of dimensions as others.

A NetCDF file is a container for simple array based data structures. There is
limited capacity ^[I.e. it's not possible to query a file for an arbitrary sparse set of values, without constructing a degenerate hyperslab query for each point or reading a hyperslab containing cells not in the set.  Do you know different? Please let me know!] in the formal API for accessing data randomly
within a variable, the primary mechanism is to define offset and stride (start
and count) hyperslab indexes.



## tidync

Tidync aims to ease exploration of the contents of a NetCDF file and provides methods extract arbitrary hyperslabs. These can be used directly in array contexts, or in "long form" database contexts. 

On first contact with the file, the available variables are classified by grid and
dimension.  The "active" grid is the one that queries may be made against, and may be changed with the `activate` function. 


```{r}
library(tidync)
tidync(ice_file)
```

Here we see variables are clearly grouped by the *grid* they exist in, where
grid is a specific (and ordered!) set of dimensions. This allows us to see the
set of variables that implicitly co-exist, they have the same *shape*.  The
first grid "D0,D1,D2" has two variables, *concentration* and *quality_flag*, and
the second "D2" has only one variable *time*. There are no general rules here, a
file might have any number of dimensions and variables, and any variable might
be defined by one or more dimensions.

In this case the D2 grid has only one variable in its single dimension, and it
happens to be a special kind of variable - a "coordinate dimension", as
indicated by the `coord_dim` flag. In the traditional `ncdump` summary above
it's easy to see there's only really one data grid, in `ni,nj,time` that it
holds two variables, and that time is a special coordinate dimension - in
contrast neither `ni` or `nj` have an explicit 1-dimension variable. When there
are many dimensions and/or many variables those patterns are *not* easy to see.

A particular grid was chosen by default, this is the "D0,D1,D2" grid composed of
3 dimensions, generally the largest grid will be chosen as that is *usually* the
target we would be after. To choose a different grid we may nominate it by name,
or by member variable.

By name we choose a grid composed of only one dimension, and the summary print
makes a distinction based on which dimensions are *active*.

```{r activate}
tidync(ice_file) %>% activate("D2")
```

It is also possible to choose a grid by *variable name*, and in tidyverse style
here we do not need to quote the name. If we choose a variable in the default
grid, then it's no different to the first case.

```{r  NSE-activate}
tidync(ice_file) %>% activate(time)

## choose grid by variable name, which happens to be the default grid here
tidync(ice_file) %>% activate(quality_flag)

## same as the default
tidync(ice_file) %>% activate("D0,D1,D2")
```

The data variables available on a grid can be expanded out as a single data
frame, which all the coordinates copied out - this is not efficient(!) but if we
craft our queries sensibly to read only what we need, it's a very easy way to
explore the data in a file.

The 'hyper_filter' function allows specification of expressions to subset a variable based on each dimension's coordinate values.  

If no expressions are included we are presented with a table containing a row
for each dimension, its extent in coordinates and its length. For convenience we
also assign the activate form to an R variable, though we could just chain the
entire operation without this.

```{r}
concentration <- tidync(ice_file) 
concentration %>% hyper_filter() 
```


By specifying inequality expressions we see an *implicit* subsetting of the
array. Everything so far is implicit to delay any file-based computation
required to actually interact with the file and read from it.

Notice that these are "name = expr" paired expressions, because the right hand
side may be quite general we need the left hand side name to be assured of the
name of the dimension referred to.

```{r}
concentration %>% hyper_filter(nj = nj < 20)
```

We can also use the special internal variable 'index', which will test against
position in the dimension elements '1:length' rather than the values. It's not
different in this case because ni and nj are just position dimensions anyway.
The special 'dplyr' adverbs like 'between' will work.

```{r}
concentration %>% 
  hyper_filter(ni = index < 50, 
               nj = dplyr::between(index, 30, 100))
```

## Data extraction

How to use these idioms to extract actual data?

We can now exercise these variable choice and dimension filters to return actual
data, either in by slicing out a  "slab" in array-form, or as a data frame.

```{r}
hf <- concentration %>% 
  hyper_filter(ni = index > 150, 
               nj = dplyr::between(index, 30, 100))

## as an array
arr <- hf %>% hyper_array()
str(arr)

## as a data frame
hf %>% 
  hyper_tibble() %>% 
  dplyr::filter(!is.na(concentration)) %>% dplyr::distinct(concentration, quality_flag)
```

## Interrogating data by dimensions

The connection object 'hf' is available for determining what is available and how we might cut into it. 'time' interestingly is of length 1, so perhaps adds nothing to the information about this otherwise 2D data set. If we were to open this file in `ncdf4` or `RNetCDF` and wanted to take a subset of the file out, we would have to specify and `start` of 1 and a `count` of ` even though it's completely redundant. 


```{r sea-ice-example}
hf
```


The 'start' and 'count' values reported are directly
useable by the traditional API tools, and in particular by the functions
`ncdf4::ncvar_get()` (`varid`, `start`, `count`), its counterpart in
`RNetCDF::var.get.nc()` and command line tools like CDO.


But, it's a pain to have to know the dimension of the variable and specify every slot. Code in the traditional API that looks like this

```{r eval= FALSE}
## WARNING, pseudocode
var_get(con, variable, start = c(1, 1, 1), count = c(10, 5, 1))
```

Obviously, it should be reasonable to specify only the count of any dimensions that *we don't want the entirety of*. This problem manifests exactly in R arrays generally, we can't provide information only about the dimensions we want, they have to be specified explicitly - even if we mean *all of it*. 

It does not matter what we include in the filter query, it can be all, some or none of the available dimensions, in any order. 

```{r dimension-index,eval=TRUE}
hf %>% hyper_filter(nj = index < 20, ni = ni > 20)

hf %>% hyper_filter(nj = index < 20)
```

The other requirement we have is the ability to automate these task, and so far we have only interacted with the dimensionality information from the print out. For programming, there are functions `hyper_vars()`, `hyper_dims()` and `hyper_grids()` to report on elements in the source. The value of `hyper_vars()` and `hyper_dims()` is restricted to the *active* grid. The function `hyper_grids()` reports all available grids, with the currently active one indicated. The name of the current grid is also available via `active()`. 


```{r}
hyper_vars(hf)
hyper_dims(hf)

## change the active grid
hf %>% activate("D2") %>% 
  hyper_vars()

active(hf)

hf %>% activate("D2") %>%
  active()
```

## Transforms

Under the hood tidync manages the relationship between dimensions and coordinates via *transforms* tables. This is a 1-dimensional mapping between dimension index and its coordinate value (if there's no explicit value it is assumed to be equal to the 1-based index). These tables are available via the `hyper_transforms()` function. 

Each table has columns `<dimension-name>`, `index`, `id`, `name`, `coord_dim` (whether the dimension has explicit coordinates, in the `<dimension-name>` column), and `selected`. The `selected` column records which of the dimension elements is currently requested by a `hyper_filter` query, and by default is set to `TRUE`. Expressions in the filter function work by updating this column. 


```{r}
hyper_transforms(hf)
```

## Future improvements

* support multiple sources at once for lazy read of a virtual grid
* support groups
* support compound types
* allow a group-by function for a polygon layer, against a pair of dimensions to
classify cells
* allow a truly DBI and dplyr level of lazy read, with more filter, select,
mutate and collect idioms
* provide converters to raster format, stars format


[not-a-format]: https://twitter.com/TedHabermann/status/958034585002041344
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidync_data.R
\name{print.tidync_data}
\alias{print.tidync_data}
\title{Print tidync data}
\usage{
\method{print}{tidync_data}(x, ...)
}
\arguments{
\item{x}{'tidync_data' object (from \code{\link[=hyper_array]{hyper_array()}})}

\item{...}{reserved args}
}
\value{
the input object invisibly
}
\description{
Print method for the 'tidync_data' list of arrays returned by \code{\link[=hyper_array]{hyper_array()}}.
}
\details{
The output lists the variables and their dimensions of an object from a
previous call to \code{\link[=tidync]{tidync()}}, and possibly \code{\link[=hyper_filter]{hyper_filter()}}. The available
data will differ from the source in terms of variables (via \code{select_var} in
\link{hyper_array}) and the lengths of each dimension (via named expressions in
\code{\link[=hyper_filter]{hyper_filter()}}).
}
\examples{
argofile <- system.file("extdata/argo/MD5903593_001.nc", package = "tidync")
argodata <- tidync(argofile) \%>\% hyper_filter(N_LEVELS = index < 5) \%>\% 
              hyper_array(select_var = c("TEMP_ADJUSTED", "PRES"))
print(argodata)
}
\seealso{
tidync_data
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidync-package.R
\docType{package}
\name{tidync-package}
\alias{tidync-package}
\title{Tidy tools for NetCDF data.}
\description{
Provides easy to use idioms for working with NetCDF data for extraction,
manipulation and visualization. NetCDF is Network Common Data Form
\url{https://www.unidata.ucar.edu/software/netcdf/}.
}
\details{
See \code{\link[=print.tidync]{print.tidync()}} for details on the printed version of a tidync object.

There is a family of functions "hyper_verb" around exploring and
extracting data.
\tabular{ll}{
\code{\link{active}} \tab report the currently active grid \cr
\code{\link{activate}} \tab active a grid \cr
\code{\link{tidync}} \tab core NetCDF source object for tidync functions\cr
\code{\link{hyper_filter}} \tab apply dimension expressions to specify array slices\cr
\code{\link{hyper_array}} \tab extracts a raw data array based on a NetCDF index \cr
\code{\link{hyper_tbl_cube}} \tab extracts data as a dplyr tbl_cube \cr
\code{\link{hyper_tibble}} \tab extracts data as a data frame with all dimension values\cr
\code{\link{hyper_transforms}} \tab extract the active (or all) dimension transforms\cr
\code{\link{hyper_vars}} \tab information on active variables \cr
\code{\link{hyper_dims}} \tab information on active dimensions \cr
\code{\link{hyper_grids}} \tab information on grids \cr
}
The scheme generally processes dimension filters into NetCDF extraction
indexes and these are always available to each function, and are expressed
in printed output.
}
\examples{
argofile <- system.file("extdata/argo/MD5903593_001.nc", package = "tidync")
argo <- tidync(argofile)
argo \%>\% active()
argo \%>\% activate("D3,D8") \%>\% hyper_array()
argo \%>\% hyper_filter(N_LEVELS = index < 4)
argo \%>\% hyper_tbl_cube()
argo \%>\% hyper_tibble(select_var = c("TEMP_QC"))
argo \%>\% hyper_transforms()
argo \%>\% hyper_vars()
argo \%>\% hyper_dims()
argo \%>\% hyper_grids()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/activate.R
\name{activate}
\alias{activate}
\alias{active}
\alias{active<-}
\alias{activate.tidync}
\alias{active.tidync}
\alias{active.default}
\alias{active<-.default}
\title{Activate a NetCDF grid}
\usage{
activate(.data, what, ..., select_var = NULL)

\method{activate}{tidync}(.data, what, ..., select_var = NULL)

\method{active}{tidync}(x)

active(x)

\method{active}{default}(x)

active(x) <- value

\method{active}{default}(x) <- value
}
\arguments{
\item{.data}{NetCDF object}

\item{what}{name of a grid or variable}

\item{...}{reserved, currently ignored}

\item{select_var}{optional argument to set selected state of variable/s by
name}

\item{x}{NetCDF object}

\item{value}{name of grid or variable to be active}
}
\value{
NetCDF object
}
\description{
A grid in NetCDF is a particular shape and size available for array
variables, and consists of sets of dimensions. To activate a grid is to set
the context for downstream operations, for querying, summarizing and reading
data. There's no sense in performing these operations on more than one grid
at a time, but multiple variables may exist in a single grid. There may be
only one significant grid in a source or many, individual dimensions are
themselves grids.
}
\details{
There may be more than one grid and one is always activated by default. A
grid may be activated by name in the form of 'D1,D0' where one or more
numbered dimensions indicates the grid. The grid definition names are printed
as part of the summary of in the tidync object and may be obtained directly
with \code{\link[=hyper_grids]{hyper_grids()}} on the tidync object.

Activation of a grid sets the context for downstream operations (slicing and
reading data) from NetCDF, and as there may be several grids in a single
source activation allows a different choice of available variables.  By
default the largest grid is activated. Once activated, all downstream tasks
apply to the set of variables that exist on that grid.

If \code{\link[=activate]{activate()}} is called with a variable name, it puts the variable first.
The function \code{\link[=active]{active()}} gets and sets the active grid. To restrict ultimate
read to particular variables use the \code{select_var} argument to
\code{\link[=hyper_filter]{hyper_filter()}}, \code{\link[=hyper_tibble]{hyper_tibble()}} and \code{\link[=hyper_tbl_cube]{hyper_tbl_cube()}}.

Scalar variables are not currently available to tidync, and it's not obvious
how activation would occur for scalars, but in future perhaps \code{activate("S")}
could be the right way forward.
}
\examples{
if (!tolower(Sys.info()[["sysname"]]) == "sunos") {
 l3file <- "S20080012008031.L3m_MO_CHL_chlor_a_9km.nc"
 rnc <- tidync(system.file("extdata", "oceandata", l3file,
 package = "tidync"))
 activate(rnc, "palette")

 ## extract available grid names
 hyper_grids(rnc)
}
}
\seealso{
hyper_filter hyper_tibble hyper_tbl_cube
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hyper_filter.R
\name{hyper_filter}
\alias{hyper_filter}
\alias{hyper_filter.tidync}
\title{Subset NetCDF variable by expression}
\usage{
hyper_filter(.x, ...)

\method{hyper_filter}{tidync}(.x, ...)
}
\arguments{
\item{.x}{NetCDF file, connection object, or \code{tidync} object}

\item{...}{currently ignored}
}
\value{
data frame
}
\description{
The \code{\link[=hyper_filter]{hyper_filter()}} acts on a \link{tidync} object by matching one or more
filtering expressions like with \code{dplyr::filter}. This allows us to lazily
specify a subset from a NetCDF array without pulling  any data. The modified
object may be printed to see the effects of subsetting, or saved for further
use.
}
\details{
The function \code{\link[=hyper_filter]{hyper_filter()}} will act on an existing tidync object or a
source string.

Filter arguments must be named as per the dimensions in the variable in form
\code{dimname = dimname < 10}. This is a restrictive variant of \code{\link[dplyr:filter]{dplyr::filter()}},
with a syntax more like \code{\link[dplyr:mutate]{dplyr::mutate()}}. This ensures that each element is
named, so we know which dimension to apply this to, but also that the
expression evaluated against can do some extra work for a nuanced test.

There are special columns provided with each axis, one is 'index' so that
exact matching can be done by position, or to ignore the actual value of the
coordinate. That means we can use a form like \code{dimname = index < 10} to
subset by position in the array index, without necessarily knowing the
values along that dimension.
}
\examples{
f <- "S20080012008031.L3m_MO_CHL_chlor_a_9km.nc"
l3file <- system.file("extdata/oceandata", f, package= "tidync")
## filter by value
tidync(l3file) \%>\% hyper_filter(lon = lon < 100)
## filter by index
tidync(l3file) \%>\% hyper_filter(lon = index < 100)

## be careful that multiple comparisons must occur in one expression
 tidync(l3file) \%>\% hyper_filter(lon = lon < 100 & lon > 50)

## filter in combination/s
tidync(l3file) \%>\% hyper_filter(lat = abs(lat) < 10, lon = index < 100)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hyper_array.R
\name{hyper_array}
\alias{hyper_array}
\alias{tidync_data}
\alias{hyper_slice}
\alias{hyper_array.tidync}
\alias{hyper_array.character}
\title{Extract NetCDF data as an array}
\usage{
hyper_array(
  x,
  select_var = NULL,
  ...,
  raw_datavals = FALSE,
  force = FALSE,
  drop = TRUE
)

hyper_slice(
  x,
  select_var = NULL,
  ...,
  raw_datavals = FALSE,
  force = FALSE,
  drop = TRUE
)

\method{hyper_array}{tidync}(
  x,
  select_var = NULL,
  ...,
  raw_datavals = FALSE,
  force = FALSE,
  drop = TRUE
)

\method{hyper_array}{character}(x, select_var = NULL, ..., raw_datavals = FALSE, drop = TRUE)
}
\arguments{
\item{x}{NetCDF file, connection object, or \link{tidync} object}

\item{select_var}{optional vector of variable names to select}

\item{...}{passed to \code{\link[=hyper_filter]{hyper_filter()}}}

\item{raw_datavals}{logical to control whether scaling in the NetCDF is
applied or not}

\item{force}{ignore caveats about large extraction and just do it}

\item{drop}{collapse degenerate dimensions, defaults to \code{TRUE}}
}
\description{
Extract the raw array data as a list of  one or more arrays. This can be the
entire variable/s or after dimension-slicing using \code{\link[=hyper_filter]{hyper_filter()}}
expressions. This is a delay-breaking function and causes data to be read
from the source into R native arrays. This list of arrays is
lightly classed as \link{tidync_data}, with methods for \code{\link[=print]{print()}} and \code{\link[=tidync]{tidync()}}.
}
\details{
The function \code{\link[=hyper_array]{hyper_array()}} is used by \code{\link[=hyper_tibble]{hyper_tibble()}} and \code{\link[=hyper_tbl_cube]{hyper_tbl_cube()}}
to actually extract data arrays from NetCDF, if a result would be particularly large
there is a check made and user-opportunity to cancel. This is controllable as an
option \code{getOption('tidync.large.data.check')}, and can be set to never check with
\code{options(tidync.large.data.check = FALSE)}.

The function \code{\link[=hyper_array]{hyper_array()}} will act on an existing tidync object or a source
string.

By default all variables in the active grid are returned, use \code{select_var} to
specify one or more desired variables.

The transforms are stored as a list of tables in an attribute `transforms``,
access these with \code{\link[=hyper_transforms]{hyper_transforms()}}.
}
\examples{
f <- "S20080012008031.L3m_MO_CHL_chlor_a_9km.nc"
l3file <- system.file("extdata/oceandata", f, package= "tidync")

## extract a raw list by filtered dimension
library(dplyr)
araw1 <- tidync(l3file) \%>\%
 hyper_filter(lat = between(lat, -78, -75.8), 
              lon = between(lon, 165, 171)) \%>\%
 hyper_array()

araw <- tidync(l3file) \%>\% 
         hyper_filter(lat = abs(lat) < 10, 
                     lon = index < 100) \%>\%
  hyper_array()

## hyper_array will pass the expressions to hyper_filter
braw <- tidync(l3file) \%>\% 
  hyper_array(lat = abs(lat) < 10, lon = index < 100)

## get the transforms tables (the axis coordinates)
lapply(attr(braw, "transforms"), 
   function(x) nrow(dplyr::filter(x, selected)))
## the selected axis coordinates should match in order and in size
lapply(braw, dim)
}
\seealso{
\link{print.tidync_data} for a description of the print summary,
\code{\link[=hyper_tbl_cube]{hyper_tbl_cube()}} and \code{\link[=hyper_tibble]{hyper_tibble()}} which are also delay-breaking
functions that cause data to be read
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidync.R
\name{print.tidync}
\alias{print.tidync}
\title{Print tidync object}
\usage{
\method{print}{tidync}(x, ...)
}
\arguments{
\item{x}{NetCDF object}

\item{...}{reserved}
}
\description{
Provide a summary of variables and dimensions, organized by their 'grid' (or
'shape') and with a summary of any slicing operations provided as 'start' and
'count' summaries for each dimension in the active grid.
}
\details{
See \link[=tidync]{tidync} for detail about the object, and
\link[=hyper_vars]{hyper_vars} for programmatic access to the active grid's
variable and dimension information.

The print summary is organized in two sections, the first is available grids
(sets of dimensions) and their associated variables, the second is the
dimensions, separated into active and inactive. All dimensions may be active
in some NetCDF sources.

Individual \emph{active} dimensions include the following components: * 'dim'    -
dimension label, D0, D1, D2, ... * 'name'   - dimension name * 'length' -
size of the dimension * 'min'    - minimum value of the dimension * 'max' -
maximum value of the dimension * 'start'  - start index of subsetting *
'count'  - length of subsetting index * 'dmin'   - minimum value of the
subset dimension * 'dmax'   - maximum value of the subset dimension * 'unlim'
\itemize{
\item indicates whether dimension is unlimited (spread across other files,
usually the time-step) * 'coord_dim' - indicates whether dimension is a
coordinate-dimension (i.e. listed as a 1-D grid)
}

The \emph{inactive} dimension summary does not include 'start', 'count', 'dmin',
'dmax' as these are identical to the values of 1, 'length', 'min', 'max' when
no array subsetting has been applied.
}
\examples{
argofile <- system.file("extdata/argo/MD5903593_001.nc", package = "tidync")
argo <- tidync(argofile)
print(argo)

## the print is modified by choosing a new grid or running filters
argo \%>\% activate("D7,D9,D11,D8")

argo \%>\% hyper_filter(N_LEVELS = index > 300)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hyper_vars_dims.R
\name{hyper_vars}
\alias{hyper_vars}
\alias{hyper_dims}
\alias{hyper_grids}
\title{Grid status}
\usage{
hyper_vars(x, ...)

hyper_dims(x, ...)

hyper_grids(x, ...)
}
\arguments{
\item{x}{tidync object}

\item{...}{ignored}
}
\value{
data frame
}
\description{
Functions to report on the current status of the \code{active} grid. Information
on the active dimensions and variables are listed in a data frame with
multiple columns.
}
\details{
The dimensions and variables of the active grid are identified in the
\link[=print.tidync]{print} method of the tidync object, these functions exist to
provide that information directly.

\code{hyper_vars()} will list the ids, data type, name, dimension number, number
of attributes and and coordinate status of the variables on the currently
active grid.

\code{hyper_dims()} will list the names, lengths, start/count index, ids, and
status of dimensions on the currently active grid. records on the currently
active dimensions.

\code{hyper_grids()} will list the names, number of dimension, and number of
variables and active status of each grid in the source.
}
\examples{
f <- "S20080012008031.L3m_MO_CHL_chlor_a_9km.nc"
l3file <- system.file("extdata/oceandata", f, package= "tidync")
tnc <- tidync(l3file)
hyper_vars(tnc)
hyper_dims(tnc)
hyper_dims(tnc \%>\% hyper_filter(lat = lat < 20))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hyper_tbl_cube.R
\name{hyper_tbl_cube}
\alias{hyper_tbl_cube}
\alias{hyper_tbl_cube.tidync}
\alias{hyper_tbl_cube.character}
\title{A dplyr cube tbl}
\usage{
hyper_tbl_cube(x, ..., force = FALSE)

\method{hyper_tbl_cube}{tidync}(x, ..., force = FALSE)

\method{hyper_tbl_cube}{character}(x, ..., force = FALSE)
}
\arguments{
\item{x}{tidync object}

\item{...}{arguments for \code{\link[=hyper_filter]{hyper_filter()}}}

\item{force}{ignore caveats about large extraction and just do it}
}
\value{
tbl_cube

\code{dplyr::tbl_cube}
}
\description{
Produce a \link[cubelyr:tbl_cube]{tbl_cube} from NetCDF. This is a
delay-breaking function and causes data to be read from the source
into the tbl cube format defined in the \link[cubelyr:tbl_cube]{dplyr}
package.
}
\details{
The size of an extraction is checked and if \emph{quite large} there is an a user-controlled
prompt to proceed or cancel. This can be disabled with \code{options(tidync.large.data.check = FALSE)}
\itemize{
\item please see \code{\link[=hyper_array]{hyper_array()}} for more details.
}

The tbl cube is a very general and arbitrarily-sized array that
can be used with tidyverse functionality. Dimension coordinates are
stored with the tbl cube, derived from the grid
\link[=hyper_transforms]{transforms}.
}
\examples{
f <- "S20080012008031.L3m_MO_CHL_chlor_a_9km.nc"
l3file <- system.file("extdata/oceandata", f, package= "tidync")
(cube <- hyper_tbl_cube(tidync(l3file) \%>\%
activate(chlor_a), lon = lon > 107, lat = abs(lat) < 30))
ufile <- system.file("extdata", "unidata", "test_hgroups.nc", 
 package = "tidync", mustWork = TRUE)
 
## some versions of NetCDF don't support this file
## (4.1.3 tidync/issues/82)
group_nc <- try(tidync(ufile), silent = TRUE)
if (!inherits(group_nc, "try-error")) {
 res <-  hyper_tbl_cube(tidync(ufile))
 print(res)
} else {
 ## the error was
 writeLines(c(group_nc))
}
}
\seealso{
\code{\link[=hyper_array]{hyper_array()}} and \code{\link[=hyper_tibble]{hyper_tibble()}} which are also delay-breaking
functions that cause data to be read
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidync.R, R/tidync_data.R
\name{tidync}
\alias{tidync}
\alias{tidync.character}
\alias{tidync.tidync_data}
\title{Tidy NetCDF}
\usage{
tidync(x, what, ...)

\method{tidync}{character}(x, what, ...)

\method{tidync}{tidync_data}(x, what, ...)
}
\arguments{
\item{x}{path to a NetCDF file}

\item{what}{(optional) character name of grid (see \code{ncmeta::nc_grids}) or
(bare) name of variable (see \code{ncmeta::nc_vars}) or index of grid to
\code{activate}}

\item{...}{reserved for arguments to methods, currently ignored}
}
\description{
Connect to a NetCDF source and allow use of \verb{hyper_*} verbs for slicing with
\code{\link[=hyper_filter]{hyper_filter()}}, extracting data with \code{\link[=hyper_array]{hyper_array()}} and  [hyper_tibble()
from an activated grid. By default the largest \emph{grid} encountered is
activated, see\code{\link[=activate]{activate()}}.
}
\details{
The print method for tidync includes a lot of information about which
variables exist on which dimensions, and if any slicing (\code{\link[=hyper_filter]{hyper_filter()}})
operations have occurred these are summarized as 'start' and 'count'
modifications relative to the dimension lengths. See \link[=print.tidync]{print}
for these details, and \link[=hyper_vars]{hyper_vars} for programmatic access to
this information

Many NetCDF forms are supported and tidync tries to reduce the interpretation
applied to a given source. The NetCDF system defines a 'grid' for storing
array data, where 'grid' is the array 'shape', or 'set of dimensions'). There
may be several grids in a single source and so we introduce the concept of
grid 'activation'. Once activated, all downstream tasks apply to the set of
variables that exist on that grid.

NetCDF sources with numeric types are chosen by default, even if existing
'NC_CHAR' type variables are on the largest grid. When read any 'NC_CHAR'
type variables are exploded into single character elements so that dimensions
match the source.
}
\section{Grids}{
 A grid is an instance of a particular set of dimensions,
which can be shared by more than one variable. This is not the 'rank' of a
variable (the number of dimensions) since a single data set may have many
3D variables composed of different sets of axes/dimensions. There's no
formality around the concept of 'shape', as far as we know.

A dimension may have length zero, but this is a special case for a
"measure" dimension, we think. (It doesn't mean the product of the
dimensions is zero, for example).
}

\section{Limitations}{
 Files with compound types are not yet supported and
should fail gracefully. Groups are not yet supported.

We haven't yet explored 'HDF5' in detail, so any feedback is appreciated.
Major use of compound types is made by \url{https://github.com/sosoc/croc}.
}

\examples{
## a SeaWiFS (S) Level-3 Mapped (L3m) monthly (MO) chlorophyll-a (CHL)
## remote sensing product at 9km resolution (at the equator)
## from the NASA ocean colour group in NetCDF4 format (.nc)
## for 31 day period January 2008 (S20080012008031) 
f <- "S20080012008031.L3m_MO_CHL_chlor_a_9km.nc"
l3file <- system.file("extdata/oceandata", f, package= "tidync")
## skip on Solaris
if (!tolower(Sys.info()[["sysname"]]) == "sunos") {
tnc <- tidync(l3file)
print(tnc)
}

## very simple Unidata example file, with one dimension
\dontrun{
uf <- system.file("extdata/unidata", "test_hgroups.nc", package = "tidync")
recNum <- tidync(uf) \%>\% hyper_tibble()
print(recNum)
}
## a raw grid of Southern Ocean sea ice concentration from IFREMER
## it is 12.5km resolution passive microwave concentration values
## on a polar stereographic grid, on 2 October 2017, displaying the 
## "hole in the ice" made famous here:
## https://tinyurl.com/ycbchcgn
ifr <- system.file("extdata/ifremer", "20171002.nc", package = "tidync")
ifrnc <- tidync(ifr)
ifrnc \%>\% hyper_tibble(select_var = "concentration")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/magrittr-pipe.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{\%>\%}
\title{Objects exported from other packages}
\examples{
system.file("extdata/argo/MD5903593_001.nc", package = "tidync") \%>\% 
     tidync()
}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{magrittr}{\code{\link[magrittr]{\%>\%}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hyper_tibble.R
\name{hyper_tibble}
\alias{hyper_tibble}
\alias{hyper_tibble.character}
\alias{hyper_tibble.tidync}
\title{Extract NetCDF data as an expanded table.}
\usage{
hyper_tibble(x, ..., na.rm = TRUE, force = FALSE)

\method{hyper_tibble}{character}(x, ..., na.rm = TRUE, force = FALSE)

\method{hyper_tibble}{tidync}(x, ..., na.rm = TRUE, force = FALSE)
}
\arguments{
\item{x}{NetCDF file, connection object, or \code{tidync} object}

\item{...}{arguments to `hyper_filter``}

\item{na.rm}{if \code{TRUE} these rows are not included in the output when all
variables are \code{NA}}

\item{force}{ignore caveats about large extraction and just do it}
}
\value{
a \code{tbl_df}
}
\description{
Extract the raw array data as an expanded data frame. This can be the entire
variable/s or after dimension-slicing using \code{\link[=hyper_filter]{hyper_filter()}} expressions with
dimension values expanded appropriately for each element in the arrays (one
row per element).
}
\details{
The size of an extraction is checked and if \emph{quite large} there is an a user-controlled
prompt to proceed or cancel. This can be disabled with \code{options(tidync.large.data.check = FALSE)}
\itemize{
\item please see \code{\link[=hyper_array]{hyper_array()}} for more details.
}

The function \code{\link[=hyper_tibble]{hyper_tibble()}} will act on an existing tidync object or a source
string.

By default all variables in the active grid are returned, use \code{select_var} to
limit.
}
\examples{
l3file <- "S20080012008031.L3m_MO_CHL_chlor_a_9km.nc"
f <- system.file("extdata", "oceandata", l3file, package= "tidync")
rnc <- tidync(f)
hyper_filter(rnc)
library(dplyr)
lapply(hyper_array(f, lat = lat > 0, lon = index > 3000), dim)

 ht <- hyper_tibble(rnc) \%>\%
 filter(!is.na(chlor_a))
ht
library(ggplot2)
ggplot(ht \%>\% filter(!is.na(chlor_a)),
aes(x = lon, y = lat, fill = chlor_a)) + geom_tile()
}
\seealso{
\code{\link[=hyper_array]{hyper_array()}} and \code{\link[=hyper_tbl_cube]{hyper_tbl_cube()}}  which are also delay-breaking
functions that cause data to be read
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hyper_transforms.R
\name{hyper_transforms}
\alias{hyper_transforms}
\alias{hyper_transforms.default}
\title{Axis transforms}
\usage{
hyper_transforms(x, all = FALSE, ...)

\method{hyper_transforms}{default}(x, all = FALSE, ...)
}
\arguments{
\item{x}{tidync object}

\item{all}{set to \code{TRUE} to return all transforms, not only active ones}

\item{...}{ignored}
}
\value{
list of axis transforms
}
\description{
Axis 'transforms' are data frames of each dimension in a NetCDF source.
\code{hyper_transforms} returns a list of the active transforms by default,
}
\details{
Each transform is available by name from a list, and each data frame has the
coordinate of the dimension, its index, and a 'selected' variable set by the
filtering expressions in \code{hyper_filter} and used by the read-functions
\code{hyper_array} and \code{hyper_tibble}.

Use \code{hyper_transforms} to interrogate and explore the available dimension
manually, or for development of custom functions.
}
\examples{
l3file <- "S20080012008031.L3m_MO_CHL_chlor_a_9km.nc"
f <- system.file("extdata", "oceandata", l3file, package = "tidync")
ax <- tidync(f) \%>\% hyper_transforms()
names(ax)
lapply(ax, dim)

## this function returns the transforms tidync knows about for this source
str(tidync(f)$transforms)
names(hyper_transforms(tidync(f), all = TRUE))
}
