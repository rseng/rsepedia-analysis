<!-- README.md is generated from README.Rmd. Please edit that file -->
RGPDD
=====

***DEPRECATED*** 

This package is no longer necessary, as the GPDD can now be downloaded directly as `.csv` files from the KNB repostitory: <https://doi.org/10.5063/F1BZ63Z8>


This package originally provided an R interface the [Global Population Dynamics Database](http://www3.imperial.ac.uk/cpb/databases/gpdd), where data was served from a no-longer-maintained API and provided only in the Microsoft Access DataBase format, both of which made access more difficult for R users.  

## Test environments
* local OS X install, R 3.1.2
* ubuntu 12.04 (on travis-ci), R 3.1.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking dependencies in R code ... NOTE
  Namespace in Imports field not imported from: 'R6'

  R6 is a build-time dependency.

## Downstream dependencies
I have also run R CMD check on downstream dependencies of httr 
(https://github.com/wch/checkresults/blob/master/httr/r-release). All packages 
that I could install passed except:

* XYZ:...
Data downloaded as MDB file.  
