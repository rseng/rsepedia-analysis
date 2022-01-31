
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Webchem

<!-- badges: start -->

[![status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN](https://www.r-pkg.org/badges/version/webchem)](https://CRAN.R-project.org/package=webchem)
[![R build
status](https://github.com/ropensci/webchem/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/webchem/actions)
[![Coverage](https://codecov.io/github/ropensci/webchem/coverage.svg?branch=master)](https://codecov.io/gh/ropensci/webchem/branch/master)
[![Downloads](https://cranlogs.r-pkg.org/badges/webchem)](https://cran.r-project.org/package=webchem)
[![Total
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/webchem?color=blue)](https://cran.r-project.org/package=webchem)

<!-- badges: end -->

`webchem` is a R package to retrieve chemical information from the web.
This package interacts with a suite of web APIs to retrieve chemical
information.

The functions in the package that hit a specific API have a prefix and
suffix separated by an underscore (`prefix_suffix()`). They follow the
format of `source_functionality`, with the exception of functions that
retrieve database identifiers which follow the format of
`get_identifier`. e.g.`cs_compinfo` uses ChemSpider to retrieve compound
informations and `get_csid()` retrieves ChemSpider IDs.

## Chemical databases currently accessed by webchem

At least some of the data in the following sources is accesible through
`webchem` functions. To learn more about what is available, browse the
documentation
[here](https://docs.ropensci.org/webchem/reference/index.html).

-   [BCPC Compendium of Pesticide Common
    Names](https://pesticidecompendium.bcpc.org) (formerly Alan Wood’s
    Compendium of Pesticide Common Names)
-   [ChEBI](https://www.ebi.ac.uk/chebi/)
-   [Chemical Identifier Resolver
    (CIR)](https://cactus.nci.nih.gov/chemical/structure)
-   [Chemical Translation Service
    (CTS)](http://cts.fiehnlab.ucdavis.edu/)
-   [ChemIDplus](https://chem.nlm.nih.gov/chemidplus/)
-   [ChemSpider](http://www.chemspider.com/) (requires an [API
    token](https://developer.rsc.org/))
-   [ETOX](http://webetox.uba.de/webETOX/index.do)
-   [Flavornet](http://www.flavornet.org)
-   [NIST](https://webbook.nist.gov) (currently gas chromatography
    retention indices only)
-   [OPSIN](http://opsin.ch.cam.ac.uk/instructions.html)
-   [PAN Pesticide Database](https://www.pesticideinfo.org/)
-   [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
-   [U.S. EPA Substance Registry Service
    (SRS)](https://cdxnodengn.epa.gov/cdx-srs-rest/)
-   [Wikidata](https://www.wikidata.org/wiki/Wikidata:WikiProject_Chemistry)

#### API keys

Some ChemSpider functions require an API key. Please register at RSC
(<https://developer.rsc.org/>) to retrieve an API key.

## Installation

#### Install from CRAN (stable version)

``` r
install.packages("webchem")
```

#### Install from Github (development version)

``` r
install.packages("devtools")
library("devtools")
install_github("ropensci/webchem")
```

### Use Cases

Have you used `webchem` in your work? Please let us know by opening an
issue or making a pull request to edit this section!

-   Allaway RJ, La Rosa S, Guinney J, Gosline SJC (2018) Probing the
    chemical–biological relationship space with the Drug Target
    Explorer. Journal of Cheminformatics 10:41.
    <https://doi.org/10.1186/s13321-018-0297-4>
-   Bergmann AJ, Points GL, Scott RP, et al (2018) Development of
    quantitative screen for 1550 chemicals with GC-MS. Anal Bioanal Chem
    410:3101–3110. <https://doi.org/10.1007/s00216-018-0997-7>
-   Brokl M, Morales V, Bishop L, et al (2019) Comparison of Mainstream
    Smoke Composition from CR20 Resin Filter and Empty-Cavity Filter
    Cigarettes by Headspace SPME Coupled with GC×GC TOFMS and
    Chemometric Analysis. Beiträge zur Tabakforschung
    International/Contributions to Tobacco Research 28:231–249.
    <https://doi.org/10.2478/cttr-2019-0004>
-   Münch D, Galizia CG (2016) DoOR 2.0 - Comprehensive Mapping of
    Drosophila melanogaster Odorant Responses. Scientific Reports
    6:21841. <https://doi.org/10.1038/srep21841>

### Citation

If you use `webchem` in a publication, please cite our paper:

-   Szöcs E, Stirling T, Scott ER, et al (2020) webchem: An R Package to
    Retrieve Chemical Information from the Web. J Stat Soft 93:.
    <https://doi.org/10.18637/jss.v093.i13>

### Acknowledgements

Without the fantastic web services `webchem` wouldn’t be here.
Therefore, kudos to the web service providers and developers! Please
remember to acknowledge these data resources in your work using
`webchem`.

### Want to contribute?

Check out our [contribution guide
here](https://github.com/ropensci/webchem/blob/master/CONTRIBUTING.md).

### Meta

-   Please [report any issues, bugs or feature
    requests](https://github.com/ropensci/webchem/issues).
-   License: MIT
-   Get citation information for `webchem` in R with
    `citation("webchem")`
-   Please note that this package is released with a [Contributor Code
    of Conduct](https://ropensci.org/code-of-conduct/). By contributing
    to this project, you agree to abide by its terms.

[![ropensci](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
# webchem 1.1.2.9001

## BUG FIXES

* get_cid() became more robust to smiles queries with special characters.

# webchem 1.1.2

## NEW FEATURES

* Export chemical structures in Mol format with write_mol().

## BUG FIXES

* ci_query() can no longer query chemicals by name.
* Non-exported function ping_pubchem_pw() was incorrectly reporting that PUG VIEW was down.  This has been fixed.
* is.cas() now catches whitespaces correctly.
* aw_query() was renamed and adapted to bcpc_query, as the alanwood site has moved

## MINOR IMPROVEMENTS

* webchem functions now default to global options regarding verbose messages.

# webchem 1.1.1.

## NEW FEATURES

* Fetch LIPID MAPS and SwissLipids identifiers from Wikidata.

## BUG FIXES

* Fix get_csid() so it doesn't break when a query is invalid.

# webchem 1.1.0.

## NEW FEATURES

* Download images of substances from Chemical Identifier Resolver (CIR) with `cir_img()`.
* Download images of substances from ChemSpider with `cs_img()`.
* `find_db()` checks if a query gets a hit in most databases integrated in webchem. Useful for deciding which of several databases to focus on given a set of chemicals.

## MINOR IMPROVEMENTS

* Most functions now use httr::RETRY() to access webservices.
* Verbose messages are now harmonized.
* The `"type"` argument in `ci_query()` and `aw_query()` has been changed to `"from"` for consistency with other functions.
* `fn_percept()` and `cts_compinfo()` now have `"query"` and `"from"` arguments for consistency with other functions.
* Possible values for `"from"` have been made more consistent across functions.
* `pc_synonyms()`, `cts_convert()`, `cir_query()` have been changed to use the `match` argument instead of `choices` for consistency with other functions.
* `get_etoxid()` output changed slightly so that the matched chemical name string no longer includes the etoxid in parentheses.
* `is.cas()` is now vectorized.

## BUG FIXES

* Fix URL encoding so SMILES queries don't fail on some special characters.

# webchem 1.0.0

## NEW FEATURES

* get_cid() now can search by registry IDs (e.g. CAS RN), and can handle more complex requests like searching for similar compounds.
* Retrieve chemical data from PubChem content pages with pc_sect().
* get_etoxid() now can search by CAS, EC, GSBL and RTECS numbers. Added `from = ` argument.
* nist_ri() now can search by name, InChI, InChIKey, or CAS.  The `cas` argument is deprecated.  Use `query` instead with `from = "cas"`.

## MINOR IMPROVEMENTS

* All `get_*()` functions now output tibbles with a column for the query and a column for the retrieved ID.
* Changes to arguments in `get_*()` functions to make them more consistent.
* aw_idx.rda is no longer included in the package as a data set. Instead, it is built by build_aw_idx() to tempdir().


## BUG FIXES

* nist_ri() returned malformed tables or errored if there was only one entry for a query.
* get_csid() now returns all csids when queried from formula.
* get_csid() returned an error when query was NA.
* get_chebiid() and chebi_comp_entity() fixed for invalid queries.
* get_cid() returned the PubChem ID of sodium when the query was NA.
* aw_query() returned a list for successful queries, NA for unsuccessful queries.

## DEPRECATED FUNCTIONS

## DEFUNCT FUNCTIONS

# webchem 0.5.0


## NEW FEATURES

* Retrieve data from ChEBI (https://www.ebi.ac.uk/chebi/) webservice with get_chebiid() and chebi_comp_entity(). ChEBI comprises a rich data base on chemicals with bilogical interest.
* Retrieve retention indices from NIST (https://webbook.nist.gov) with nist_ri().
* Get record details from US EPA Substance Registry Services (https://cdxnodengn.epa.gov/cdx-srs-rest/) with srs_query().
* "first" argument in cts_convert() and cir_query() and "interactive" argument in pc_synonyms() deprecated.  Use "choices" instead to return either a list of all results, only the first result, or an interactive menu to choose a result to return.
* ChemSpider functions now look for an API token stored in .Renviron or .Rprofile by default so you can keep them hidden more easily.

## MINOR IMPROVEMENTS

* as.cas() added.
* Removed documentation files for non-exported functions that were only used internally.

## BUG FIXES

* cs_prop() failed with duplicated return values.
* pp_query() failed when compound present, but no properties.
* ci_query() failed when missing table.
* get_csid() failed because of a major change in the ChemSpider API.
* multiple functions failed because of a major change in the ChemSpider API.
* cir_query() mistook NA for sodium.
* fixed functions that communicate with the ChemSpider API.
* get_etoxid() printed incorrect results for certain match types.

## DEPRECATED FUNCTIONS

* cs_extcompinfo() cannot be fixed as there is no equivalent in the new ChemSpider API yet.

## DEFUNCT FUNCTIONS

* ppdb_parse() has been removed. webchem no longer offers any support for PPDB.
* pp_query() has been removed. Physprop API is no longer active.
* cs_prop() has been removed.

# webchem 0.4.0

## NEW FEATURES

## MINOR IMPROVEMENTS

## BUG FIXES

* extr_num() did not work properly with decimal numbers.
* cs_prop() failed when epi-suite data was not available.
* cs_prop() failed with invalid html.
* cs_prop() gave incorrect answer, if entries were not available.
* cs_prop() did not parse scientific number correctly.
* is.smiles() failed because of changes in rcdk.
* cir_query() failed with identifiers containing spaces (e.g. 'acetic acid').
* Aeveral other functions failed with identifiers containing spaces & returned wrong distance.

## DEPRECATED FUNCTIONS

## DEFUNCT FUNCTIONS


# webchem 0.3.0


## NEW FEATURES

## MINOR IMPROVEMENTS

* cs_prop() now also return experimental data for Boiling and Melting Points.
* pc_synonyms gained an argument 'interactive' to enter an interactive mode for selecting synonyms.
* cts_convert now returns NA if no matches are found.

## BUG FIXES

* cs_prop() failed with some CSIDs.
* wd_ident() failed if multiple entries where found. Now returns the first hit only.
* ci_query() did not return fully cleaned smiles and inchi.

## DEPRECATED FUNCTIONS

## DEFUNCT FUNCTIONS


# webchem 0.2.0


## NEW FEATURES

* fn_percept() extracts flavor percepts using CAS numbers from www.flavornet.org. Flavornet is a database of 738 compounds with human-detectible odors.

## MINOR IMPROVEMENTS

## BUG FIXES

## DEPRECATED FUNCTIONS

## DEFUNCT FUNCTIONS



# webchem 0.1.1


## NEW FEATURES
* Added ping_pubchem() to check whether pubchem is up & running.
* Added cs_web_ping () to check whether the chemspider webpage is functional.

## MINOR IMPROVEMENTS
* Updated allan wood index.

## BUG FIXES
* pc_prop() returned to many rows if last cid supplied was NA.
* Switched to https for NCBI, chemspider & chemid.
* get_wdid() failed if non-ascii characters where returned by wikipedia.
* rcdk:parse.smiles() now returns NA if a SMILES string could not be parsed.
   => broke is.smiles

## DEPRECATED FUNCTIONS

## DEFUNCT FUNCTIONS



# webchem 0.1.0


## NEW FEATURES
* Added cts_to() and cts_from() to retrieve possible ids that can be queried.
* cts_*(), pp_query(), cir_query(), get_cid(), get_etoxid(), etox_*(), pan_query() get_wdid(), aw_query(), get_csid(), cs_prop(), cs_compinfo() and ci_query() can handle multiple inputs.
* pc_prop() queries properties and pc_synonmy() synonyms from PUG-REST.
* Added extractors for webchem objects: cas(), inchikey() and smiles().


## MINOR IMPROVEMENTS
* Rewrite of pubchem functions using PUG-REST.
* ChemSpider: better use of NA in input (=return NA).
* More robust matching in get_etoxid.

## BUG FIXES

* pan_query() did not return numeric values.
* get_cid() failed with multiple results.

## DEPRECATED FUNCTIONS


## DEFUNCT FUNCTIONS

* ppdb_query() has been removed due to copyright issues.
The new ppdb_parse() parses only a html, but does not interact with the database.
* pan()
* alanwood()
* get_cid()
* cid_compinfo()
* chemid()
* physprop()



# webchem 0.0.5

## NEW FEATURES

* is.smiles() checks SMILES strings, by parsing via (R)CDK.
* get_wdid() and wd_indent() to retrieve information from wikidata.
* get_etoxid() can handle multi inputs (interactive mode, best match, first match, NA and all matches).
* ci_query() can handle multi inputs (interactive mode, best match, first match and NA).
* cs_prop() queries predictions (ACD and EPiSuite) from ChemSpider.

## MINOR IMPROVEMENTS

* webchem uses exclusively xml2 (instead of XML).
* All function return source_url for (micro-)attribution of sources.
* cs_compinfo(): names of returned list changed.
* cs_extcompinfo():
  - names of returned list changed.
  - result is numeric where appropriate.
* cir(): result is numeric where appropriate.
* Unified naming scheme of functions.
* is.inchikey_cs() has been integrated into is.inchikey().
* aw_query() returns multiple inchikey if found.
* pan() now returns chemical name and matched synonym.

## BUG FIXES

* Utility functions are not vectorized and throw an error.
* chemid() did mot work with inchikey as input.
* ppdb_idx returned duplicated CAS values, which caused ppdb() to fail.
* ppdb() failed in some cases because of false encoding.
* etox_*() functions are more robust.
* ci_query() failed if multi hits were found. Now returns first hit.
* aw_fuery() failed if inchikey was not found.

## DEPRECATED FUNCTIONS

* pan_query() replaces pan().
* aw_query() replaces alanwood().
* get_pcid() replaces get_cid().
* pc_compinfo() replaces cid_compinfo().
* ci_query() replaces chemid().
* pp_query() replaces physprop().

## DEFUNCT FUNCTIONS

* csid_compinfo()
* csid_extcompinfo()




# webchem 0.0.4


## NEW FEATURES

* chemid() to query ChemIDplus http://chem.sis.nlm.nih.gov/chemidplus/.
* is.inchikey() and is.cas() to check if a string is valid inchikey or CAS registry number.
* parse_mol(): A simple molfile parser.
* Functions to work with ChemSpider InChI API:
  + cs_csid_mol() : convert csid to mol.
  + cs_inchikey_csid() : convert inchikey to csid.
  + cs_inchikey_inchi() : convert inchikey to inchi.
  + cs_inchikey_mol() : convert inchikey to Molfile.
  + cs_inchi_csid() : convert inchi to csid.
  + cs_inchi_inchikey : convert inchi to inchikey.
  + cs_inchi_mol() : convert inchi to molfile.
  + cs_inchi_smiles() : convert inchi to smiles.
  + cs_smiles_inchi() : convert smiles to inchi.
  + These are all wrapped into cs_convert().
  + is.inchikey_cs() : Check via ChemSpider if inchikey is valid.
* webchem has now a zenodo doi, please cite if you use it.


## MINOR IMPROVEMENTS

* cts_compinfo() checks if input is a inchikey (via exported function is.inchikey()).
* cts_compinfo() is now more robust and verbose, if problems are encountered.
* alanwood() returns separate inchi and ichikeys in case of isomers.
* alanwood() returns also subactvity (e.g. $Fluazinam$activity [1] "fungicides" and $Fluazinam$subactivity [1] "pyridine fungicides").
* physprop() also returns boiling and melting points. Moreover, values are now numeric.


## BUG FIXES

* alanwood() returns only results for first match in case of multiple links found.
* physprop() stopped working after change of SRC to https, fixed now.
* Changed etox_* functions to https.


## DEPRECATED FUNCTIONS

* ppdb() replaces ppdb_query() and accepts individual index as created by ppdb_buildidx().
* cir() replaces cir_query().
* cs_compinfo() replaces csid_compinfo().
* cs_extcompinfo() replaces csid_extcompinfo().


## DEFUNCT FUNCTIONS

* allanwood()




# webchem 0.0.3


## NEW FEATURES

* Query SRC PHYSPROP Database with physprop().
* Query the ETOX ID with get_etoxid(); query basic information with etox_basic();
  quality targets with etox_targets() and test results with etox_tests().
* Query PPDB with ppdb_query().

## MINOR IMPROVEMENTS

* Added exceptions/checks to tests.
* Improved robustness of cir_query().

## BUG FIXES

* Correct the spelling of Alan Wood and rename function allanwood() to alanwood().



# webchem 0.0.2


## NEW FEATURES

* Query the PAN Pesticides Database with pan().
* Query Allan Woods Compendium of Pesticide Common Names with allanwood().

## MINOR IMPROVEMENTS

* Added checks for user input.
* Fixed documentation, added example for bulk processing.
* cts_convert() returns NA if no result was found.
* Set 'verbose = TRUE' as default for all functions.
* Added unit tests.
* All functions return silently NA, if API is not reachable.

## BUG FIXES

* cts_convert() does not ignore 'first' argument.
* get_csid() did not return NA, if there was a problem with the API.
* Many functions returned 'NA2+' if NA was given - now return NA by default.
* Many fixes in NA handling, e.g. when no hit was found.
# CONTRIBUTING #

Thank you for your interest in contributing to the project! The goal of this guide is to help you with your endevour. There are many ways to contribute and we have outlined some opportunities which might be interesting to you. If you have any questions or suggestions, feel free to contact us at <webchem@ropensci.org>.

### Fill out the survey

You can fill out our survey at https://forms.gle/V7dfGGn73dkesn5L6.

The `webchem` survey allows us to learn which databases you use and how you interact with chemical data. This is extremely valuable information for us and guides our development efforts. The survey takes about 5 minutes to fill out.

### Share a use case

Write us an e-mail and show us a full example of how you use or how you would like to use `webchem` in your data analysis! This would give us ideas about new features and also help us create better vignettes that help others get started. Please send your e-mails to <webchem@ropensci.org>.

### Raise a new issue or join discussion on an existing issue

If you found a bug either in the code or the documentation, or a data source you would like us to integrate into ```webchem```, maybe you dreamed up a new functionality that would be nice to implement, raise an issue and let's discuss it! Even if you don't have the time or the coding background to resolve the issue yourself, maybe others do, and so just by giving a good problem you might help others who are looking for interesting problems to solve. You can raise an issue [here](https://github.com/ropensci/webchem/issues). Feel free to join discussions on existing issues as well!

### Code contributions

If you know some coding, you can also add code contributions.

1. **Fork** this repo to your Github account.
2. **Clone** your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/webchem.git`.
3. Make sure to **track upstream** progress (i.e., on our version of `webchem` at `ropensci/webchem`) by doing `git remote add upstream https://github.com/ropensci/webchem.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
4. Make your **changes**. Bonus points for making changes on a new branch.

Creating new branches is good practice. This is because if you finish with a topic, open a pull request and start working on another topic without starting a new branch, any further commits you push to your account will be automatically added to your pull request as well, making it much harder for us to evaluate your request. To aboid this, open a new branch for each new topic.

5. **Push** up to your account.
6. Submit a **pull request** to home base at `ropensci/webchem`.

### Guidelines for code contributions

You can find the rOpenSci developer guide at https://devguide.ropensci.org/

We are happy to help at any point in your work.

1. We follow the [tidyverse](https://tidyverse.org) style. You can find the style guide [here](https://style.tidyverse.org/). Before committing your code, we encourage you to use ```lintr::lint_file()``` to check for nonconformances.

2. We use [`roxygen2`](https://cran.r-project.org/web/packages/roxygen2/index.html) for documentation. Please make sure you update the package to the latest version before you update the documentation with `devtools::document()`. Use `@noRd` for non exported functions.

3. Please use the [`xml2`](https://cran.r-project.org/web/packages/xml2/index.html) package instead of the `XML` package. The maintainance of xml2 is much better.

4. Please use the lightweight [`jsonlite`](https://cran.r-project.org/web/packages/jsonlite/index.html) package for handling JSON.

5. Use utilities in `webchem::utils.R` when possible to keep function style consistent across the package.

6. Be nice to the resources! Minimise interaction with the servers. Use appropriate timeouts.

7. Within test files always include a check whether the webservice is running and skip all tests when it is not. See `R/ping.R` for more details.  

Some consistency guidelines:

8. Functions that query a database for one or more database specific identifiers should follow the naming convention `get_*`, e.g. the function that queries ChEBI IDs is called `get_chebiid()`. These functions should take a vector of queries and return a single [tibble](https://cran.r-project.org/web/packages/tibble/index.html). Whenever possible these functions should have arguments `query`, `from`, `match`, `verbose` and `...`. The first column of the tibble should contain the ID-s and the last should contain the queries. Invalid queries should return a row of NA-s (apart from the last element of the row which should be the query itself).

9. The naming of functions that query a database for chemical information should start with the name of the database, followed by the functionality, e.g. `pc_synonyms()` searches for synonyms in PubChem. These functions should take a vector of queries and return a list of responses. Invalid queries should return `NA`.

10. Functions should always validate their input when appropriate. Use `match.arg()` for input validation.

11. Make sure `NA` is not confused with sodium.

12. Functions that retrieve images should follow the naming convention `*_img`, e.g. the function that retrieves images from ChemSpider is called `cs_img()`. These functions should take a vector of arguments and download images into a user defined directory. They should not keep images in memory, should not implement image processing functionality, and should not return anything to the console. Functions should include arguments `dir`, `overwrite = TRUE` and `verbose = TRUE`. `dir` should not have a default value.

13. SMILES strings may use special characters like "#". `URLencode()` does not encode this as "%23" by default, so use `URLencode(query, reserved = TRUE)` instead. It's important to note that it's the query that has to be encoded like this, not the full url.

14. Print verbose messages. Use `httr::message_for_status()` and `webchem_message()` functions to generate standard messages when possible.

### Data Sources

You might think all webscraping is perfectly legal but it is unfortunately not that simple.

Some services allow you to browse their website but do not allow you programmable access, for various reasons. Therefore, we always have to check the Terms & Conditons and any other legal documents that might restrict programmable access. `webchem` only provides access to databases where programmable access is clearly approved by the database provider. A provider might create a publicly accessible API, and if they do not have a restrictive T&C, this indicates their implicit approval for programmatically accessing their data. In all other cases explicit approval is required, i.e. either the T&C has to state that scraping is allowed, or we have to acquire written consent from the database provider before developing functions that scrape their website.

And there is a big difference between scraping and crawling. `webchem` does provide some scraping functionality but it does not provide crawling functionality.

### Thanks for contributing!# webchem 1.1.2

* https://developer.rsc.org/, SSL certificate problem: unable to get local issuer certificate: The URLs work and including these URLs in the package is important.## Pull Request

**Please do not include API keys in your pull request**. You should store any API keys or tokens in `.Renviron` or `.Rprofile` and be sure to have these files added to `.gitignore`. Do not include API keys in function examples or tests. See `cs_check_key()` for an example.

If this pull request adds a new database to the list of databases `webchem` can access, please provide some evidence that this database is OK with being accessed by a third-party---for example, a link to an About or FAQ page.

Before you submit a pull request, please do the following:

* Add an entry to NEWS concisely describing what you changed.
* If appropriate, add unit tests in the tests/testthat directory.
* Run Build->Check Package in the RStudio IDE, or `devtools::check()`, to make sure your change did not add any messages, warnings, or errors.

Doing these things will make it easier to evaluate your pull request.

We will try to be responsive and provide feedback.


## Testing

You can test the package using `devtools::test()`. 
Tests are disable in the `master` branch (commented out in `~tests\testthat.R`).
For local testing you can uncomment, or run `test_check("webchem")`.
Please do not include changes to `~tests\testthat.R` in your PR.


Delete these instructions once you have read them.

---

Brief description of the PR


PR task list:
- [ ] Update NEWS
- [ ] Add tests (if appropriate)
- [ ] Update documentation with `devtools::document()`
- [ ] Check package passed---
name: Database suggestion
about: Suggest a databaset to be integrated into webchem
---

Please answer the following questions to the best of your ability so contributors can get a sense of whether they might be excited about working on adding this functionality to `webchem`.  If you're intersted in contributing code to access this database yourself, please indicate this and also see CONTRIBUTING.md

### Suggested database: <name and link to home page>

1. **Feasibility**.  Does the database provide an API?  If not, does the database allow scraping?  This information might be in an FAQ or About page.  If you can't find it, don't worry, we'll look into it.

2. **Scope**. What kind of data does the database contain? It should be primarily chemical properties, but databases of chemical identifiers or synonyms may also be considered.

3. **Overlap**.  How much does the database contents overlap with current databases that webchem can access? Does it provide unique properties and/or data on a unique set of chemicals? 


---
name: Bug report
about: Describe a bug you've seen or problem you've encountered
---
Please briefly describe your problem and what output you expect. 
If you have a question, please try using stackoverflow <http://stackoverflow.com> first.

## Minimal reproducible example

Please include a minimal reproducible example (reprex). 
The goal of a reprex is to make it as easy as possible for me to recreate your problem so that we can fix it. 
If you've never heard of a reprex before, start by reading <https://github.com/jennybc/reprex#what-is-a-reprex> or <http://tinyurl.com/reproducible-000>. 


Delete these instructions once you have read them.

---

Brief description of the problem

```r
# insert reprex here
```---
output: github_document
editor_options: 
  chunk_output_type: console
---
<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r echo=FALSE}
knitr::opts_chunk$set(
  warning = FALSE, 
  message = FALSE,
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%",
  cache = TRUE
)
```
# Webchem

<!-- badges: start -->

[![status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN](https://www.r-pkg.org/badges/version/webchem)](https://CRAN.R-project.org/package=webchem) 
[![R build status](https://github.com/ropensci/webchem/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/webchem/actions)
[![Coverage](https://codecov.io/github/ropensci/webchem/coverage.svg?branch=master)](https://codecov.io/gh/ropensci/webchem/branch/master) 
[![Downloads](https://cranlogs.r-pkg.org/badges/webchem)](https://cran.r-project.org/package=webchem)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/webchem?color=blue)](https://cran.r-project.org/package=webchem)

<!-- badges: end -->

`webchem` is a R package to retrieve chemical information from  the web. 
This package interacts with a suite of web APIs to retrieve chemical information.

The functions in the package that hit a specific API have a prefix and suffix separated by an underscore (`prefix_suffix()`).
They follow the format of `source_functionality`, with the exception of functions that retrieve database identifiers which follow the format of `get_identifier`. e.g.`cs_compinfo` uses ChemSpider to retrieve compound informations and `get_csid()` retrieves ChemSpider IDs.

## Chemical databases currently accessed by webchem

At least some of the data in the following sources is accesible through `webchem` functions.  To learn more about what is available, browse the documentation [here](https://docs.ropensci.org/webchem/reference/index.html).

- [BCPC Compendium of Pesticide Common Names](https://pesticidecompendium.bcpc.org) (formerly Alan Wood's Compendium of Pesticide Common Names)
- [ChEBI](https://www.ebi.ac.uk/chebi/)
- [Chemical Identifier Resolver (CIR)](http://cactus.nci.nih.gov/chemical/structure)
- [Chemical Translation Service (CTS)](http://cts.fiehnlab.ucdavis.edu/)
- [ChemIDplus](https://chem.nlm.nih.gov/chemidplus/) 
- [ChemSpider](http://www.chemspider.com/) (requires an [API token]((https://developer.rsc.org/)))
- [ETOX](http://webetox.uba.de/webETOX/index.do)
- [Flavornet](http://www.flavornet.org) 
- [NIST](https://webbook.nist.gov) (currently gas chromatography retention indices only)
- [OPSIN](http://opsin.ch.cam.ac.uk/instructions.html)
- [PAN Pesticide Database](http://www.pesticideinfo.org/)
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
- [U.S. EPA Substance Registry Service (SRS)](https://cdxnodengn.epa.gov/cdx-srs-rest/)
- [Wikidata](https://www.wikidata.org/wiki/Wikidata:WikiProject_Chemistry)

#### API keys

Some ChemSpider functions require an API key. 
Please register at RSC (https://developer.rsc.org/) to retrieve an API key.

## Installation
#### Install from CRAN (stable version)

```{r install_cran, eval=FALSE}
install.packages("webchem")
```


#### Install from Github (development version)

```{r install_github, eval=FALSE}
install.packages("devtools")
library("devtools")
install_github("ropensci/webchem")
```

### Use Cases

Have you used `webchem` in your work?  Please let us know by opening an issue or making a pull request to edit this section!

- Allaway RJ, La Rosa S, Guinney J, Gosline SJC (2018) Probing the chemical–biological relationship space with the Drug Target Explorer. Journal of Cheminformatics 10:41. https://doi.org/10.1186/s13321-018-0297-4
- Bergmann AJ, Points GL, Scott RP, et al (2018) Development of quantitative screen for 1550 chemicals with GC-MS. Anal Bioanal Chem 410:3101–3110. https://doi.org/10.1007/s00216-018-0997-7
- Brokl M, Morales V, Bishop L, et al (2019) Comparison of Mainstream Smoke Composition from CR20 Resin Filter and Empty-Cavity Filter Cigarettes by Headspace SPME Coupled with GC×GC TOFMS and Chemometric Analysis. Beiträge zur Tabakforschung International/Contributions to Tobacco Research 28:231–249. https://doi.org/10.2478/cttr-2019-0004
- Münch D, Galizia CG (2016) DoOR 2.0 - Comprehensive Mapping of Drosophila melanogaster Odorant Responses. Scientific Reports 6:21841. https://doi.org/10.1038/srep21841

### Citation

If you use `webchem` in a publication, please cite our paper:

- Szöcs E, Stirling T, Scott ER, et al (2020) webchem: An R Package to Retrieve Chemical Information from the Web. J Stat Soft 93:. https://doi.org/10.18637/jss.v093.i13


### Acknowledgements
Without the fantastic web services `webchem` wouldn't be here. Therefore, kudos to the web service providers and developers! Please remember to acknowledge these data resources in your work using `webchem`.

### Want to contribute?

Check out our [contribution guide here](https://github.com/ropensci/webchem/blob/master/CONTRIBUTING.md).

### Meta

- Please [report any issues, bugs or feature requests](https://github.com/ropensci/webchem/issues).
- License: MIT
- Get citation information for `webchem` in R with `citation("webchem")`
- Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
---
title: "Getting started with webchem"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Getting started with webchem}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
<!-- Only edit webchem.Rmd.orig and compile using the script in vignettes/precompile.R -->




```r
library(webchem)
library(dplyr)
```

The `lc50` dataset provided with `webchem` contains acute ecotoxicity of 124 insecticides.  We'll work with a subset of these to obtain chemical names and octanal/water partitioning coefficients from PubChem, and gas chromatography retention indices from the NIST Web Book.


```r
head(lc50)
#>        cas        value
#> 4  50-29-3    12.415277
#> 12 52-68-6     1.282980
#> 15 55-38-9    12.168138
#> 18 56-23-5 35000.000000
#> 21 56-38-2     1.539119
#> 36 57-74-9    98.400000

lc50_sub <- lc50[1:15, ]
```

## Getting Identifiers

Usually a `webchem` workflow starts with translating and retrieving chemical identifiers since most chemical information databases use their own internal identifiers.

First, we will covert CAS numbers to InChIKey identifiers using the Chemical Translation Service.  Then, we'll use these InChiKeys to get Pubchem CompoundID numbers, to use for retrieving chemical properties from PubChem.


```r
lc50_sub$inchikey <- cts_convert(lc50_sub$cas, from = "CAS", to = "InChIKey", choices = 1, verbose = FALSE)
head(lc50_sub)
#>        cas        value                    inchikey
#> 4  50-29-3    12.415277 YVGGHNCTFXOJCH-UHFFFAOYSA-N
#> 12 52-68-6     1.282980 NFACJZMKEDPNKN-UHFFFAOYSA-N
#> 15 55-38-9    12.168138 PNVJTZOFSHSLTO-UHFFFAOYSA-N
#> 18 56-23-5 35000.000000 VZGDMQKNWNREIO-UHFFFAOYSA-N
#> 21 56-38-2     1.539119 LCCNCVORNKJIRZ-UHFFFAOYSA-N
#> 36 57-74-9    98.400000 BIWJNBZANLAXMG-YQELWRJZSA-N
any(is.na(lc50_sub$inchikey))
#> [1] FALSE
```

Great, now we can retrieve PubChem CIDs.  All `get_*()` functions return a data frame containing the query and the retrieved identifier.  We can merge this with our dataset with `dplyr::full_join()`


```r
x <- get_cid(lc50_sub$inchikey, from = "inchikey", match = "first", verbose = FALSE)
library(dplyr)
lc50_sub2 <- full_join(lc50_sub, x, by = c("inchikey" = "query"))
head(lc50_sub2)
#>       cas        value                    inchikey      cid
#> 1 50-29-3    12.415277 YVGGHNCTFXOJCH-UHFFFAOYSA-N     3036
#> 2 52-68-6     1.282980 NFACJZMKEDPNKN-UHFFFAOYSA-N     5853
#> 3 55-38-9    12.168138 PNVJTZOFSHSLTO-UHFFFAOYSA-N     3346
#> 4 56-23-5 35000.000000 VZGDMQKNWNREIO-UHFFFAOYSA-N     5943
#> 5 56-38-2     1.539119 LCCNCVORNKJIRZ-UHFFFAOYSA-N      991
#> 6 57-74-9    98.400000 BIWJNBZANLAXMG-YQELWRJZSA-N 11954021
```

## Retrieving Chemical Properties

Functions that query chemical information databases begin with a prefix that matches the database.  For example, functions to query PubChem begin with `pc_` and functions to query ChemSpider begin with `cs_`. In this example, we'll get the names and log octanal/water partitioning coefficients for each compound using PubChem, and the WHO acute toxicity rating from the PAN Pesticide database.


```r
y <- pc_prop(lc50_sub2$cid, properties = c("IUPACName", "XLogP"))
#> https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/IUPACName,XLogP/JSON
y$CID <- as.character(y$CID)
lc50_sub3 <- full_join(lc50_sub2, y, by = c("cid" = "CID"))
head(lc50_sub3)
#>       cas        value                    inchikey      cid
#> 1 50-29-3    12.415277 YVGGHNCTFXOJCH-UHFFFAOYSA-N     3036
#> 2 52-68-6     1.282980 NFACJZMKEDPNKN-UHFFFAOYSA-N     5853
#> 3 55-38-9    12.168138 PNVJTZOFSHSLTO-UHFFFAOYSA-N     3346
#> 4 56-23-5 35000.000000 VZGDMQKNWNREIO-UHFFFAOYSA-N     5943
#> 5 56-38-2     1.539119 LCCNCVORNKJIRZ-UHFFFAOYSA-N      991
#> 6 57-74-9    98.400000 BIWJNBZANLAXMG-YQELWRJZSA-N 11954021
#>                                                                      IUPACName XLogP
#> 1                  1-chloro-4-[2,2,2-trichloro-1-(4-chlorophenyl)ethyl]benzene   6.9
#> 2                                 2,2,2-trichloro-1-dimethoxyphosphorylethanol   0.5
#> 3 dimethoxy-(3-methyl-4-methylsulfanylphenoxy)-sulfanylidene-lambda5-phosphane   4.1
#> 4                                                           tetrachloromethane   2.8
#> 5                    diethoxy-(4-nitrophenoxy)-sulfanylidene-lambda5-phosphane   3.8
#> 6            (1R,7S)-1,3,4,7,8,9,10,10-octachlorotricyclo[5.2.1.02,6]dec-8-ene   4.9
```

The IUPAC names are long and unwieldy, and one could use `pc_synonyms()` to choose better names. Several other functions return synonyms as well, even though they are not explicitly translator type functions.  We'll see an example of that next.

Many of the chemical databases `webchem` can query contain vast amounts of information in a variety of structures.  Therefore, some `webchem` functions return nested lists rather than data frames.  `pan_query()` is one such function.


```r
out <- pan_query(lc50_sub3$cas, verbose = FALSE)
#> Warning in lapply(out[tonum], as.numeric): NAs introduced by coercion

#> Warning in lapply(out[tonum], as.numeric): NAs introduced by coercion

#> Warning in lapply(out[tonum], as.numeric): NAs introduced by coercion

#> Warning in lapply(out[tonum], as.numeric): NAs introduced by coercion

#> Warning in lapply(out[tonum], as.numeric): NAs introduced by coercion

#> Warning in lapply(out[tonum], as.numeric): NAs introduced by coercion

#> Warning in lapply(out[tonum], as.numeric): NAs introduced by coercion
```

`out` is a nested list which you can inspect with `View()`.  It has an element for each query, and within each query, many elements corresponding to different properties in the database.  To extract a single property from all queries, we need to use a mapping function such as `sapply()` or one of the `map_*()` functions from the `purrr` package.


```r
lc50_sub3$who_tox <- sapply(out, function(y) y$`WHO Acute Toxicity`)
lc50_sub3$common_name <- sapply(out, function(y) y$`Chemical name`)

# #equivalent with purrr package:
# lc50_sub3$who_tox <- map_chr(out, pluck, "WHO Acute Toxicity")
# lc50_sub3$common_name <- map_chr(out, pluck, "Chemical name")
```


```r
#tidy up columns
lc50_done <- dplyr::select(lc50_sub3, common_name, cas, inchikey, XLogP, who_tox)
head(lc50_done)
#>            common_name     cas                    inchikey XLogP                  who_tox
#> 1            DDT, p,p' 50-29-3 YVGGHNCTFXOJCH-UHFFFAOYSA-N   6.9 II, Moderately Hazardous
#> 2          Trichlorfon 52-68-6 NFACJZMKEDPNKN-UHFFFAOYSA-N   0.5 II, Moderately Hazardous
#> 3             Fenthion 55-38-9 PNVJTZOFSHSLTO-UHFFFAOYSA-N   4.1 II, Moderately Hazardous
#> 4 Carbon tetrachloride 56-23-5 VZGDMQKNWNREIO-UHFFFAOYSA-N   2.8               Not Listed
#> 5            Parathion 56-38-2 LCCNCVORNKJIRZ-UHFFFAOYSA-N   3.8  Ia, Extremely Hazardous
#> 6            Chlordane 57-74-9 BIWJNBZANLAXMG-YQELWRJZSA-N   4.9 II, Moderately Hazardous
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chemspider.R
\name{cs_control}
\alias{cs_control}
\title{Control ChemSpider API requests}
\usage{
cs_control(
  datasources = vector(),
  order_by = "default",
  order_direction = "default",
  include_all = FALSE,
  complexity = "any",
  isotopic = "any"
)
}
\arguments{
\item{datasources}{character; specifies the databases to query. Use
\code{cs_datasources()} to retrieve available ChemSpider data sources.}

\item{order_by}{character; specifies the sort order for the results.
Valid values are \code{"default"}, \code{"recordId"}, \code{"massDefect"},
\code{"molecularWeight"}, \code{"referenceCount"}, \code{"dataSourceCount"},
\code{"pubMedCount"}, \code{"rscCount"}.}

\item{order_direction}{character; specifies the sort order for the results.
Valid values are \code{"default"}, \code{"ascending"}, \code{"descending"}.}

\item{include_all}{logical; see details.}

\item{complexity}{character; see details.
Valid values are \code{"any"} \code{"single"}, \code{"multiple"}.}

\item{isotopic}{character; see details.
Valid values are \code{"any"}, \code{"labeled"}, \code{"unlabeled"}.}
}
\value{
Returns a list of specified control options.
}
\description{
For some ChemSpider API requests, you can also specify various control
options. This function is used to set these control options.
}
\details{
The only function that currently uses \code{databases} is
\code{get_csid()} and only when you query a CSID from a formula. This
parameter is disregarded in all other queries.

Setting \code{include_all} to \code{TRUE} will consider records
which contain all of the filter criteria specified in the request. Setting
it to \code{FALSE} will consider records which contain any of the filter
criteria.

A compound with a  \code{complexity} of \code{"multiple"} has more
than one disconnected system in it or a metal atom or ion.
}
\note{
This is a full list of all API control options.
However, not all of these options are used in all functions.
Each API uses a subset of these controls.
The controls that are available for a given function are indicated within the
documentation of the function.
}
\examples{
cs_control()
cs_control(order_direction = "descending")
}
\references{
\url{https://developer.rsc.org/docs/compounds-v1-trial/1/overview}
}
\seealso{
\code{\link{get_csid}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{is.cas}
\alias{is.cas}
\title{Check if input is a valid CAS}
\usage{
is.cas(x, verbose = getOption("verbose"))
}
\arguments{
\item{x}{character; input CAS}

\item{verbose}{logical; print messages during processing to console?}
}
\value{
a logical
}
\description{
This function checks if a string is a valid CAS registry number.
A valid CAS is 1) separated by two hyphes into three parts; 2) the first part
consists from two up to seven digits; 3) the second of two digits; 4) the
third of one digit (check digit); 5) the check digits corresponds the
checksum. The checksum is found by taking the last digit (excluding the check
digit) multiplyingit with 1, the second last multiplied with 2, the
third-last multiplied with 3 etc. The modulo 10 of the sum of these is the
checksum.
}
\note{
This function can only handle one CAS string
}
\examples{
is.cas('64-17-5')
is.cas('64175')
is.cas('4-17-5')
is.cas('64-177-6')
is.cas('64-17-55')
is.cas('64-17-6')
}
\references{
Eduard Szöcs, Tamás Stirling, Eric R. Scott, Andreas Scharmüller,
Ralf B. Schäfer (2020). webchem: An R Package to Retrieve Chemical
Information from the Web. Journal of Statistical Software, 93(13).
\doi{10.18637/jss.v093.i13}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chemspider.R
\name{cs_datasources}
\alias{cs_datasources}
\title{Retrieve ChemSpider data sources}
\usage{
cs_datasources(apikey = NULL, verbose = getOption("verbose"))
}
\arguments{
\item{apikey}{character; your API key. If NULL (default),
\code{cs_check_key()} will look for it in .Renviron or .Rprofile.}

\item{verbose}{should a verbose output be printed on the console?}
}
\value{
Returns a character vector.
}
\description{
The function returns a vector of available data sources used by ChemSpider.
Some ChemSpider functions allow you to restrict which sources are used to
lookup the requested query. Restricting the sources makes these queries
faster.
}
\note{
An API key is needed. Register at \url{https://developer.rsc.org/}
for an API key. Please respect the Terms & Conditions. The Terms & Conditions
can be found at \url{https://developer.rsc.org/terms}.
}
\examples{
\dontrun{
cs_datasources()
}
}
\references{
\url{https://developer.rsc.org/docs/compounds-v1-trial/1/overview}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opsin.R
\name{opsin_query}
\alias{opsin_query}
\title{OPSIN web interface}
\usage{
opsin_query(query, verbose = getOption("verbose"), ...)
}
\arguments{
\item{query}{character;  chemical name that should be queryed.}

\item{verbose}{logical; should a verbose output be printed on the console?}

\item{...}{currently not used.}
}
\value{
a tibble with six columnns: "query", inchi", "stdinchi", "stdinchikey", "smiles", "message", and "status"
}
\description{
Query the OPSIN  (Open Parser for Systematic IUPAC nomenclature) web service
\url{https://opsin.ch.cam.ac.uk/instructions.html}.
}
\examples{
\donttest{
opsin_query('Cyclopropane')
opsin_query(c('Cyclopropane', 'Octane'))
opsin_query(c('Cyclopropane', 'Octane', 'xxxxx'))
}
}
\references{
Lowe, D. M., Corbett, P. T., Murray-Rust, P., & Glen, R. C. (2011).
Chemical Name to Structure: OPSIN, an Open Source Solution. Journal of Chemical Information and Modeling,
51(3), 739–753. \doi{10.1021/ci100384d}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chebi.R
\name{chebi_comp_entity}
\alias{chebi_comp_entity}
\title{Retrieve Complete Entity from ChEBI}
\usage{
chebi_comp_entity(chebiid, verbose = getOption("verbose"), ...)
}
\arguments{
\item{chebiid}{character; search term (i.e. chebiid).}

\item{verbose}{logical; should a verbose output be printed on the console?}

\item{...}{optional arguments}
}
\value{
returns a list of data.frames or lists containing a complete ChEBI
entity
}
\description{
Returns a list of Complete ChEBI entities.
ChEBI data are parsed as data.frames ("properties", "chebiid_snd",
"synonyms", "iupacnames", "formulae", "regnumbers", "citations", "dblinks",
"parents", "children", "comments", "origins") or
as a list ("chem_structure") in the list.
The SOAP protocol is used \url{https://www.ebi.ac.uk/chebi/webServices.do}.
}
\examples{
\donttest{
# might fail if API is not available
chebi_comp_entity('CHEBI:27744')

# multiple inputs
comp <- c('CHEBI:27744', 'CHEBI:27744')
chebi_comp_entity(comp)

}
}
\references{
Hastings J, Owen G, Dekker A, Ennis M, Kale N, Muthukrishnan V,
  Turner S, Swainston N, Mendes P, Steinbeck C. (2016). ChEBI in 2016:
  Improved services and an expanding collection of metabolites. Nucleic Acids
  Res.

  Hastings, J., de Matos, P., Dekker, A., Ennis, M., Harsha, B., Kale, N.,
  Muthukrishnan, V., Owen, G., Turner, S., Williams, M., and Steinbeck, C.
  (2013) The ChEBI reference database and ontology for biologically relevant
  chemistry: enhancements for 2013. Nucleic Acids Res.

  de Matos, P., Alcantara, R., Dekker, A., Ennis, M., Hastings, J., Haug, K.,
  Spiteri, I., Turner, S., and Steinbeck, C. (2010) Chemical entities of
  biological interest: an update. Nucleic Acids Res. Degtyarenko, K.,
  Hastings, J., de Matos, P., and Ennis, M. (2009). ChEBI: an open
  bioinformatics and cheminformatics resource. Current protocols in
  bioinformatics / editoral board, Andreas D. Baxevanis et al., Chapter 14.

  Degtyarenko, K., de Matos, P., Ennis, M., Hastings, J., Zbinden, M.,
  McNaught, A., Alcántara, R., Darsow, M., Guedj, M. and Ashburner, M. (2008)
  ChEBI: a database and ontology for chemical entities of biological
  interest. Nucleic Acids Res. 36, D344–D350.

Eduard Szöcs, Tamás Stirling, Eric R. Scott, Andreas Scharmüller,
Ralf B. Schäfer (2020). webchem: An R Package to Retrieve Chemical
Information from the Web. Journal of Statistical Software, 93(13).
\doi{10.18637/jss.v093.i13}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chemid.R
\name{ci_query}
\alias{ci_query}
\title{Retrieve information from ChemIDPlus}
\usage{
ci_query(query, from = c("rn", "inchikey"), verbose = getOption("verbose"))
}
\arguments{
\item{query}{character; query string}

\item{from}{character; type of query string, can be one of \code{"rn"} for
CAS registry numbers or \code{"inchikey"}.}

\item{verbose}{logical; should a verbose output be printed on the console?}
}
\value{
A list of 8 entries: name (vector), synonyms (vector), cas (vector),
inchi (vector), inchikey (vector), smiles(vector), toxicity (data.frame),
physprop (data.frame) and source_url.
}
\description{
Retrieve information from ChemIDPlus
\url{https://chem.nlm.nih.gov/chemidplus}
}
\note{
Please respect the Terms and Conditions of the National Library of
Medicine, \url{https://www.nlm.nih.gov/databases/download.html}.
}
\examples{
\dontrun{
# might fail if API is not available
y1 <- ci_query('50-00-0', from = 'rn')
y1[['50-00-0']]$inchikey

# query by inchikey
y2 <- ci_query('WSFSSNUMVMOOMR-UHFFFAOYSA-N', from = 'inchikey')
y2[[1]]$name

# query multiple compounds
comps <- c("50-00-0", "64-17-5")
y3 <- ci_query(comps, from = "rn")

# extract log-P
sapply(y3, function(y){
 if (length(y) == 1 && is.na(y))
   return(NA)
 y$physprop$Value[y$physprop$`Physical Property` == 'log P (octanol-water)']
 })
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chebi.R
\name{get_chebiid}
\alias{get_chebiid}
\title{Retrieve Lite Entity (identifiers) from ChEBI}
\usage{
get_chebiid(
  query,
  from = c("all", "chebi id", "chebi name", "definition", "name", "iupac name",
    "citations", "registry numbers", "manual xrefs", "automatic xrefs", "formula",
    "mass", "monoisotopic mass", "charge", "inchi", "inchikey", "smiles", "species"),
  match = c("all", "best", "first", "ask", "na"),
  max_res = 200,
  stars = c("all", "two only", "three only"),
  verbose = getOption("verbose"),
  ...
)
}
\arguments{
\item{query}{character; search term.}

\item{from}{character; type of input.  \code{"all"} searches all types and
\code{"name"} searches all names. Other options include \code{'chebi id'},
\code{'chebi name'}, \code{'definition'}, \code{'iupac name'},
\code{'citations'}, \code{'registry numbers'}, \code{'manual xrefs'},
\code{'automatic xrefs'}, \code{'formula'}, \code{'mass'},
\code{'monoisotopic mass'},\code{'charge'}, \code{'inchi'},
\code{'inchikey'}, \code{'smiles'}, and \code{'species'}}

\item{match}{character; How should multiple hits be handled?, \code{"all"}
all matches are returned, \code{"best"} the best matching (by the ChEBI
searchscore) is returned, \code{"ask"} enters an interactive mode and the
user is asked for input, \code{"na"} returns NA if multiple hits are found.}

\item{max_res}{integer; maximum number of results to be retrieved from the
web service}

\item{stars}{character; "three only" restricts results to those manualy
annotated by the ChEBI team.}

\item{verbose}{logical; should a verbose output be printed on the console?}

\item{...}{currently unused}
}
\value{
returns a list of data.frames containing a chebiid, a chebiasciiname,
  a searchscore and stars if matches were found. If not, data.frame(NA) is
  returned
}
\description{
Returns a data.frame with a ChEBI entity ID (chebiid),
a ChEBI entity name (chebiasciiname), a search score (searchscore) and
stars (stars) using the SOAP protocol:
\url{https://www.ebi.ac.uk/chebi/webServices.do}
}
\examples{
\donttest{
# might fail if API is not available
get_chebiid('Glyphosate')
get_chebiid('BPGDAMSIGCZZLK-UHFFFAOYSA-N')

# multiple inputs
comp <- c('Iron', 'Aspirin', 'BPGDAMSIGCZZLK-UHFFFAOYSA-N')
get_chebiid(comp)

}
}
\references{
Hastings J, Owen G, Dekker A, Ennis M, Kale N, Muthukrishnan V,
  Turner S, Swainston N, Mendes P, Steinbeck C. (2016). ChEBI in 2016:
  Improved services and an expanding collection of metabfolites. Nucleic
  Acids Res.

  Hastings, J., de Matos, P., Dekker, A., Ennis, M., Harsha, B., Kale, N.,
  Muthukrishnan, V., Owen, G., Turner, S., Williams, M., and Steinbeck, C.
  (2013) The ChEBI reference database and ontology for biologically relevant
  chemistry: enhancements for 2013. Nucleic Acids Res.

  de Matos, P., Alcantara, R., Dekker, A., Ennis, M., Hastings, J., Haug, K.,
  Spiteri, I., Turner, S., and Steinbeck, C. (2010) Chemical entities of
  biological interest: an update. Nucleic Acids Res. Degtyarenko, K.,
  Hastings, J., de Matos, P., and Ennis, M. (2009). ChEBI: an open
  bioinformatics and cheminformatics resource. Current protocols in
  bioinformatics / editoral board, Andreas D. Baxevanis et al., Chapter 14.

  Degtyarenko, K., de Matos, P., Ennis, M., Hastings, J., Zbinden, M.,
  McNaught, A., Alcántara, R., Darsow, M., Guedj, M. and Ashburner, M. (2008)
  ChEBI: a database and ontology for chemical entities of biological
  interest. Nucleic Acids Res. 36, D344–D350.

Eduard Szöcs, Tamás Stirling, Eric R. Scott, Andreas Scharmüller,
  Ralf B. Schäfer (2020). webchem: An R Package to Retrieve Chemical
  Information from the Web. Journal of Statistical Software, 93(13).
  \doi{10.18637/jss.v093.i13}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integration.R
\name{find_db}
\alias{find_db}
\title{Check data source coverage of compounds}
\usage{
find_db(
  query,
  from,
  sources = c("etox", "pc", "chebi", "cs", "bcpc", "fn", "pan", "srs"),
  plot = FALSE
)
}
\arguments{
\item{query}{character; the search term}

\item{from}{character; the format or type of query.  Commonly accepted values
are "name", "cas", "inchi", and "inchikey"}

\item{sources}{character; which data sources to check.  Data sources are
identified by the prefix associated with webchem functions that query those
databases.  If not specified, all data sources listed will be checked.}

\item{plot}{logical; plot a graphical representation of results.}
}
\value{
a tibble of logical values where \code{TRUE} indicates that a data
  source contains a record for the query
}
\description{
Checks if entries are found in (most) data sources included in webchem
}
\examples{
\dontrun{
find_db("hexane", from = "name")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chemspider.R
\name{cs_img}
\alias{cs_img}
\title{Download images from ChemSpider}
\usage{
cs_img(
  csid,
  dir,
  overwrite = TRUE,
  apikey = NULL,
  verbose = getOption("verbose")
)
}
\arguments{
\item{csid}{numeric; the ChemSpider ID (CSID) of the substance. This will
also be the name of the image file.}

\item{dir}{character; the download directory. \code{dir} accepts both
absolute and relative paths.}

\item{overwrite}{logical; should existing files in the directory with the
same name be overwritten?}

\item{apikey}{character; your API key. If NULL (default),
\code{cs_check_key()} will look for it in .Renviron or .Rprofile.}

\item{verbose}{logical; should a verbose output be printed on the console?}
}
\description{
Retrieve images of substances from ChemSpider and export them
in PNG format.
}
\note{
An API key is needed. Register at \url{https://developer.rsc.org/}
for an API key. Please respect the Terms & Conditions. The Terms & Conditions
can be found at \url{https://developer.rsc.org/terms}.
}
\examples{
\dontrun{
cs_img(c(582, 682), dir = tempdir())
}
}
\references{
\url{https://developer.rsc.org/docs/compounds-v1-trial/1/overview}
}
\seealso{
\code{\link{get_csid}}, \code{\link{cs_check_key}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chemspider.R
\name{get_csid}
\alias{get_csid}
\title{ChemSpider ID from compound name, formula, SMILES, InChI or InChIKey}
\usage{
get_csid(
  query,
  from = c("name", "formula", "inchi", "inchikey", "smiles"),
  match = c("all", "first", "ask", "na"),
  verbose = getOption("verbose"),
  apikey = NULL,
  ...
)
}
\arguments{
\item{query}{character; search term.}

\item{from}{character; the type of the identifier to convert from. Valid
values are \code{"name"}, \code{"formula"}, \code{"smiles"},
\code{"inchi"}, \code{"inchikey"}. The default value is \code{"name"}.}

\item{match}{character; How should multiple hits be handled?, "all" all
matches are returned, "best" the best matching is returned, "ask" enters an
interactive mode and the user is asked for input, "na" returns NA if
multiple hits are found.}

\item{verbose}{logical; should a verbose output be printed on the console?}

\item{apikey}{character; your API key. If NULL (default),
\code{cs_check_key()} will look for it in .Renviron or .Rprofile.}

\item{...}{furthrer arguments passed to \code{\link{cs_control}}}
}
\value{
Returns a tibble.
}
\description{
Query one or more compunds by name, formula, SMILES, InChI or InChIKey and
return a vector of ChemSpider IDs.
}
\details{
Queries by SMILES, InChI or InChiKey do not use \code{cs_control}
  options. Queries by name use \code{order_by} and \code{order_direction}.
  Queries by formula also use \code{datasources}. See \code{cs_control()} for
  a full list of valid values for these control options.

\code{formula} can be expressed with and without LaTeX syntax.
}
\note{
An API key is needed. Register at \url{https://developer.rsc.org/} for
  an API key. Please respect the Terms & conditions:
  \url{https://developer.rsc.org/terms}.
}
\examples{
\dontrun{
get_csid("triclosan")
get_csid(c("carbamazepine", "naproxene","oxygen"))
get_csid("C2H6O", from = "formula")
get_csid("C_{2}H_{6}O", from = "formula")
get_csid("CC(O)=O", from = "smiles")
get_csid("InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)", from = "inchi")
get_csid("QTBSBXVTEAMEQO-UHFFFAOYAR", from = "inchikey")
}
}
\references{
\url{https://developer.rsc.org/docs/compounds-v1-trial/1/overview}

Eduard Szöcs, Tamás Stirling, Eric R. Scott, Andreas Scharmüller,
Ralf B. Schäfer (2020). webchem: An R Package to Retrieve Chemical
Information from the Web. Journal of Statistical Software, 93(13).
\doi{10.18637/jss.v093.i13}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cts.R
\name{cts_to}
\alias{cts_to}
\title{Return a list of all possible ids}
\usage{
cts_to(verbose = getOption("verbose"))
}
\arguments{
\item{verbose}{logical; should a verbose output be printed on the console?}
}
\value{
a character vector.
}
\description{
Return a list of all possible ids that can be used in the 'to' argument
}
\details{
See also \url{http://cts.fiehnlab.ucdavis.edu/services}
}
\examples{
\donttest{
cts_from()
}
}
\references{
Wohlgemuth, G., P. K. Haldiya, E. Willighagen, T. Kind, and O.
Fiehn 2010The Chemical Translation Service -- a Web-Based Tool to Improve
Standardization of Metabolomic Reports. Bioinformatics 26(20): 2647–2648.
}
\seealso{
\code{\link{cts_convert}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chemspider.R
\name{cs_check_key}
\alias{cs_check_key}
\title{Retrieve ChemSpider API key}
\usage{
cs_check_key()
}
\value{
an API key
}
\description{
Look for and retrieve ChemSpider API key stored in .Renviron or .Rprofile.
}
\details{
To use the any of the functions in \code{webchem} that access the
ChemSpider database, you'll need to obtain an API key. Register at
\url{https://developer.rsc.org/} for an API key. Please respect the Terms &
Conditions \url{https://developer.rsc.org/terms}.

You can store your API key as \code{CHEMSPIDER_KEY = <your key>} in
.Renviron or as \code{options(chemspider_key = <your key>)} in .Rprofile.
This will allow you to use ChemSpider without adding your API key in the
beginning of each session, and will also allow you to share your analysis
without sharing your API key. Keeping your API key hidden is good practice.
}
\examples{
\dontrun{
cs_check_key()
}
}
\seealso{
\code{\link[usethis]{edit_r_environ}}
\code{\link[usethis]{edit_r_profile}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/etox.R
\name{etox_basic}
\alias{etox_basic}
\title{Get basic information from a ETOX ID}
\usage{
etox_basic(id, verbose = getOption("verbose"))
}
\arguments{
\item{id}{character; ETOX ID}

\item{verbose}{logical; print message during processing to console?}
}
\value{
a list with lists of four entries: cas (the CAS numbers), ec (the EC
  number), gsbl (the gsbl number), a data.frame synonys with synonyms and the
  source url.
}
\description{
Query ETOX: Information System Ecotoxicology and Environmental Quality
Targets \url{https://webetox.uba.de/webETOX/index.do} for basic information
}
\note{
Before using this function, please read the disclaimer
  \url{https://webetox.uba.de/webETOX/disclaimer.do}.
}
\examples{
\dontrun{
id <- get_etoxid('Triclosan', match = 'best')
etox_basic(id$etoxid)

# Retrieve data for multiple inputs
ids <- c("20179", "9051")
out <- etox_basic(ids)
out

# extract cas numbers
sapply(out, function(y) y$cas)
}
}
\references{
Eduard Szöcs, Tamás Stirling, Eric R. Scott, Andreas Scharmüller,
Ralf B. Schäfer (2020). webchem: An R Package to Retrieve Chemical
Information from the Web. Journal of Statistical Software, 93(13).
\doi{10.18637/jss.v093.i13}.
}
\seealso{
\code{\link{get_etoxid}} to retrieve ETOX IDs,
  \code{\link{etox_basic}} for basic information, \code{\link{etox_targets}}
  for quality targets and \code{\link{etox_tests}} for test results
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/etox.R
\name{get_etoxid}
\alias{get_etoxid}
\title{Get ETOX ID}
\usage{
get_etoxid(
  query,
  from = c("name", "cas", "ec", "gsbl", "rtecs"),
  match = c("all", "best", "first", "ask", "na"),
  verbose = getOption("verbose")
)
}
\arguments{
\item{query}{character; The searchterm}

\item{from}{character; Type of input, can be one of "name" (chemical name),
"cas" (CAS Number), "ec" (European Community number for regulatory purposes),
"gsbl" (Identifier used by \url{https://www.gsbl.de}) and "rtecs" (Identifier used
by the Registry of Toxic Effects of Chemical Substances database).}

\item{match}{character; How should multiple hits be handeled? "all" returns
all matched IDs, "first" only the first match, "best" the best matching (by
name) ID, "ask" is a interactive mode and the user is asked for input, "na"
returns \code{NA} if multiple hits are found.}

\item{verbose}{logical; print message during processing to console?}
}
\value{
a tibble with 3 columns: the query, the match, and the etoxID
}
\description{
Query ETOX: Information System Ecotoxicology and Environmental Quality
Targets \url{https://webetox.uba.de/webETOX/index.do} for their substance ID
}
\note{
Before using this function, please read the disclaimer
\url{https://webetox.uba.de/webETOX/disclaimer.do}.
}
\examples{
\dontrun{
# might fail if API is not available
get_etoxid("Triclosan")
# multiple inputs
comps <- c("Triclosan", "Glyphosate")
get_etoxid(comps)
get_etoxid(comps, match = "all")
get_etoxid("34123-59-6", from = "cas") # Isoproturon
get_etoxid("133483", from = "gsbl") # 3-Butin-1-ol
get_etoxid("203-157-5", from = "ec") # Paracetamol
}
}
\references{
Eduard Szöcs, Tamás Stirling, Eric R. Scott, Andreas Scharmüller,
Ralf B. Schäfer (2020). webchem: An R Package to Retrieve Chemical
Information from the Web. Journal of Statistical Software, 93(13).
\doi{10.18637/jss.v093.i13}.
}
\seealso{
\code{\link{etox_basic}} for basic information,
\code{\link{etox_targets}} for quality targets and
\code{\link{etox_tests}} for test results.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pubchem.R
\name{pc_sect}
\alias{pc_sect}
\title{Retrieve data from PubChem content pages}
\usage{
pc_sect(
  id,
  section,
  domain = c("compound", "substance", "assay", "gene", "protein", "patent"),
  verbose = getOption("verbose")
)
}
\arguments{
\item{id}{numeric or character; a vector of PubChem identifiers to search
for.}

\item{section}{character; the section of the content page to be imported.}

\item{domain}{character; the query domain. Can be one of \code{"compound"},
\code{"substance"}, \code{"assay"}, \code{"gene"}, \code{"protein"} or
\code{"patent"}.}

\item{verbose}{logical; should a verbose output be printed on the console?}
}
\value{
Returns a tibble of query results. In the returned tibble,
\code{SourceName} is the name of the depositor, and \code{SourceID} is the
ID of the search term within the depositor's database. You can browse
\url{https://pubchem.ncbi.nlm.nih.gov/sources/} for more information about
the depositors.
}
\description{
When you search for an entity at \url{https://pubchem.ncbi.nlm.nih.gov/},
e.g. a compound or a substance, and select the record you are interested in,
you will be forwarded to a PubChem content page. When you look at a PubChem
content page, you can see that chemical information is organised into
sections, subsections, etc. The chemical data live at the lowest levels of
these sections. Use this function to retrieve the lowest level information
from PubChem content pages.
}
\details{
\code{section} is not case sensitive but it is sensitive to typing
errors and it requires the full name of the section as it is printed on the
content page. The PubChem Table of Contents Tree can also be found at
\url{https://pubchem.ncbi.nlm.nih.gov/classification/#hid=72}.
}
\note{
Please respect the Terms and Conditions of the National Library of
Medicine, \url{https://www.nlm.nih.gov/databases/download.html} the data
usage policies of National Center for Biotechnology Information,
\url{https://www.ncbi.nlm.nih.gov/home/about/policies/},
\url{https://pubchemdocs.ncbi.nlm.nih.gov/programmatic-access}, and the data
usage policies of the individual data sources
\url{https://pubchem.ncbi.nlm.nih.gov/sources/}.
}
\examples{
# might fail if API is not available
\donttest{
pc_sect(176, "Dissociation Constants")
pc_sect(c(176, 311), "density")
pc_sect(2231, "depositor-supplied synonyms", "substance")
pc_sect(780286, "modify date", "assay")
pc_sect(9023, "Ensembl ID", "gene")
pc_sect("1ZHY_A", "Sequence", "protein")
}
}
\references{
Kim, S., Thiessen, P.A., Cheng, T. et al. PUG-View: programmatic
access to chemical annotations integrated in PubChem. J Cheminform 11, 56
(2019). \doi{10.1186/s13321-019-0375-2}.
}
\seealso{
\code{\link{get_cid}}, \code{\link{pc_prop}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chemspider.R
\name{cs_extcompinfo}
\alias{cs_extcompinfo}
\title{Get extended record details by ChemSpider ID}
\usage{
cs_extcompinfo(csid, token, verbose = getOption("verbose"), ...)
}
\arguments{
\item{csid}{character,  ChemSpider ID.}

\item{token}{character; security token.}

\item{verbose}{logical; should a verbose output be printed on the console?}

\item{...}{currently not used.}
}
\value{
a data.frame with entries: 'csid', 'mf' (molecular formula),
'smiles', 'inchi' (non-standard), 'inchikey' (non-standard), 'average_mass',
'mw' (Molecular weight), 'monoiso_mass' (MonoisotopicMass), nominal_mass',
'alogp', 'xlogp', 'common_name' and 'source_url'
}
\description{
Get extended info from ChemSpider, see \url{https://www.chemspider.com/}
}
\note{
A security token is needed. Please register at RSC
\url{https://www.rsc.org/rsc-id/register}
for a security token.
Please respect the Terms & conditions
\url{https://www.rsc.org/help-legal/legal/terms-conditions/}.

use \code{\link{cs_compinfo}} to retrieve standard inchikey.
}
\examples{
\dontrun{
token <- "<redacted>"
csid <- get_csid("Triclosan")
cs_extcompinfo(csid, token)

csids <- get_csid(c('Aspirin', 'Triclosan'))
cs_compinfo(csids)
}
}
\seealso{
\code{\link{get_csid}} to retrieve ChemSpider IDs,
\code{\link{cs_compinfo}} for extended compound information.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wikidata.R
\name{wd_ident}
\alias{wd_ident}
\title{Retrieve identifiers from Wikidata}
\usage{
wd_ident(id, verbose = getOption("verbose"))
}
\arguments{
\item{id}{character; identifier, as returned by \code{\link{get_wdid}}}

\item{verbose}{logical; print message during processing to console?}
}
\value{
A data.frame of identifiers. Currently these are 'smiles', 'cas', 'cid', 'einecs', 'csid', 'inchi', 'inchikey',
'drugbank', 'zvg', 'chebi', 'chembl', 'unii', 'lipidmaps', 'swisslipids' and source_url.
}
\description{
Retrieve identifiers from Wikidata
}
\note{
Only matches in labels are returned. If more than one unique hit is found,
only the first is returned.
}
\examples{
\dontrun{
 id <- c("Q408646", "Q18216")
 wd_ident(id)
}
}
\references{
Willighagen, E., 2015. Getting CAS registry numbers out of WikiData. The Winnower.
\doi{10.15200/winn.142867.72538}

Mitraka, Elvira, Andra Waagmeester, Sebastian Burgstaller-Muehlbacher, et al. 2015
Wikidata: A Platform for Data Integration and Dissemination for the Life Sciences and beyond. bioRxiv: 031971.
}
\seealso{
\code{\link{get_wdid}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pan.R
\name{pan_query}
\alias{pan_query}
\title{Query the PAN Pesticide database}
\usage{
pan_query(
  query,
  from = c("name", "cas"),
  match = c("best", "all", "first", "na"),
  verbose = TRUE,
  ...
)
}
\arguments{
\item{query}{character; searchterm, e.g. chemical name or CAS.}

\item{from}{character; one of "name" or "cas".}

\item{match}{character; \code{match="all"} returns all matches,
\code{match="first"} the first one and \code{match="best"} (recommended) the hit with the lowest
 Levenshtein distance between query and matching synonym.}

\item{verbose}{logical; should a verbose output be printed on the console?}

\item{...}{currently not used.}
}
\value{
a named list of 73 entries,
  see \url{https://www.pesticideinfo.org/about/overview.html} for more information.
  If \code{match="best"} an additional entry \code{match_score} with the normalized
  Levenshtein distance (0 = perfect match, 1 = worst match).

CAS Number; U.S. EPAPC Code; CA ChemCode;
Use Type; Chemical Class; Molecular Weight; U.S. EPARegistered ; CA Reg Status;
PIC; POPs; WHO Obsolete; EPA HAP; CA TAC; Ground Water Contaminant;
CA Grnd Water Contam.; Acute Aquatic Toxcity; Chronic Aquatic Toxicity;
PAN BadActor Chem; Dirty Dozen; Acute Toxicity Summary; Cholinesterase Inhibitor;
Acute rating from U.S. EPA product label; U.S. NTP Acute Toxicity Studies;
Material Safety Data Sheets; TRI Acute Hazard; WHO Acute Toxicity; Cancer Rating;
U.S. EPA Carcinogens; IARC Carcinogens; U.S. NTP Carcinogens;
California Prop 65 Known Carcinogens; TRI Carcinogen;
Developmental or Reproductive Toxicity; CA Prop 65 Developmental Toxin;
U.S. TRI Developmental Toxin; CA Prop 65 Female Reproductive Toxin;
CA Prop 65 Male Reproductive Toxin ; U.S. TRI Reproductive Toxin;
Endocrine Disruption; E.U. ED Rating; Benbrook list; Denmark Inert list;
Colborn list; Illinois EPA list; Keith list; Water Solubility (Avg, mg/L);
Adsorption Coefficient (Koc); Hydrolysis Half-life (Avg, Days);
Aerobic Soil Half-life (Avg, Days); Anaerobic Soil Half-life (Avg, Days);
Maximum Contaminant Level (MCL) (ug/L); Maximum Contaminant Level Goal (MCLG) (ug/L);
One Day Exposure Health Advisory Level (ug/L); Ten Day Exposure Health Advisory Level (ug/L);
Reference Dose (ug/kg/day); U.S. Drinking Water Equivalent Level (ug/L);
Lifetime Exposure Health Advisory Level (ug/L);
Lifetime Estimated Cancer Risk (cases per 1,000,000);
Maximum Acceptable Concentration (MAC) (ug/L);
Interim Maximum Acceptable Concentration (IMAC) (ug/L);
Aesthetic Objectives (ug/L); Fresh Water Quality Criteria Continuous Exposure (ug/L);
Fresh Water Quality Criteria Maximum Peak (ug/L); Salt Water Quality Criteria Continuous Exposure (ug/L);
Salt Water Quality Criteria Max (ug/L); Human Consumption of Organisms from Water Source (ug/L);
Human Consumption of Water and Organisms from Water Source (ug/L);
Taste and Odor Criteria (ug/L);
Fresh Water Guidelines (ug/L); Salt Water Guidelines (ug/L);
Irrigation Water Guidelines (ug/L); Livestock Water Guidelines (ug/L);
Chemical Name; matching synonym; source URL
}
\description{
Retrieve information from the PAN database (\url{https://www.pesticideinfo.org/}).
This function is currently broken.
}
\examples{
\dontrun{
 # might fail if API is not available

 # return all hits
 pan_query('2,4-dichlorophenol')[[1]][c(1, 2, 5, 74)]
 # return only first hit
 pan_query('2,4-dichlorophenol', match = 'first')[[1]][c(1, 2, 5, 74)]
 # return only best hit
 pan_query('2,4-dichlorophenol', match = 'best')[[1]][c(1, 2, 5, 74)]

 out <- pan_query(c('Glyphosate', 'Rotenone'), from = "name", match = 'best')
 out

 # extract Acute Toxicity Summary
 # sapply(out, function(y) y$`Acute Toxicity Summary`)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{is.smiles}
\alias{is.smiles}
\title{Check if input is a SMILES string}
\usage{
is.smiles(x, verbose = getOption("verbose"))
}
\arguments{
\item{x}{character; input SMILES.}

\item{verbose}{logical; print messages during processing to console?}
}
\value{
a logical
}
\description{
This function checks if a string is a valid SMILES by checking
if (R)CDK can parse it. If it cannot be parsed by rcdk FALSE is returned,
else TRUE.
}
\note{
This function can handle only one SMILES string.
}
\examples{
\dontrun{
# might fail if rcdk is not working properly
is.smiles('Clc(c(Cl)c(Cl)c1C(=O)O)c(Cl)c1Cl')
is.smiles('Clc(c(Cl)c(Cl)c1C(=O)O)c(Cl)c1ClJ')
}
}
\references{
Egon Willighagen (2015). How to test SMILES strings in
Supplementary Information.
\url{https://chem-bla-ics.blogspot.nl/2015/10/how-to-test-smiles-strings-in.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chemspider.R
\name{cs_compinfo}
\alias{cs_compinfo}
\title{Retrieve record details by ChemSpider ID}
\usage{
cs_compinfo(csid, fields, verbose = getOption("verbose"), apikey = NULL)
}
\arguments{
\item{csid}{numeric; can be obtained using \code{\link{get_csid}}}

\item{fields}{character; see details.}

\item{verbose}{logical; should a verbose output be printed on the console?}

\item{apikey}{character; your API key. If NULL (default),
\code{cs_check_key()} will look for it in .Renviron or .Rprofile.}
}
\value{
Returns a data frame.
}
\description{
Submit a ChemSpider ID (CSID) and the fields you are interested in, and
retrieve the record details for your query.
}
\details{
Valid values for \code{fields} are \code{"SMILES"},
\code{"Formula"}, \code{"InChI"}, \code{"InChIKey"}, \code{"StdInChI"},
\code{"StdInChIKey"}, \code{"AverageMass"}, \code{"MolecularWeight"},
\code{"MonoisotopicMass"}, \code{"NominalMass"}, \code{"CommonName"},
\code{"ReferenceCount"}, \code{"DataSourceCount"}, \code{"PubMedCount"},
\code{"RSCCount"}, \code{"Mol2D"}, \code{"Mol3D"}. You can specify any
number of fields.
}
\note{
An API key is needed. Register at \url{https://developer.rsc.org/}
for an API key. Please respect the Terms & Conditions. The Terms & Conditions
can be found at \url{https://developer.rsc.org/terms}.
}
\examples{
\dontrun{
cs_compinfo(171, c("SMILES", "CommonName"))
cs_compinfo(171:182, "SMILES")
}
}
\references{
\url{https://developer.rsc.org/docs/compounds-v1-trial/1/overview}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{as.cas}
\alias{as.cas}
\title{Format numbers as CAS numbers}
\usage{
as.cas(x)
}
\arguments{
\item{x}{numeric vector, or character vector of CAS numbers missing the
hyphens}
}
\value{
character vector of valid CAS numbers
}
\description{
This function attempts to format numeric (or character) vectors
as character vectors of CAS numbers.  If they cannot be converted to CAS
format or don't pass \code{\link{is.cas}}, \code{NA} is returned
}
\examples{
x = c(58082, 123456, "hexenol")
as.cas(x)

}
\seealso{
\code{\link{is.cas}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nist.R
\name{nist_ri}
\alias{nist_ri}
\title{Retrieve retention indices from NIST}
\usage{
nist_ri(
  query,
  from = c("cas", "inchi", "inchikey", "name"),
  type = c("kovats", "linear", "alkane", "lee"),
  polarity = c("polar", "non-polar"),
  temp_prog = c("isothermal", "ramp", "custom"),
  cas = NULL,
  verbose = getOption("verbose")
)
}
\arguments{
\item{query}{character; the search term}

\item{from}{character; type of search term. can be one of \code{"name"},
\code{"inchi"}, \code{"inchikey"}, or \code{"cas"}. Using an identifier is
preferred to \code{"name"} since \code{NA} is returned in the event of
multiple matches to a query. Using an identifier other than a CAS number
will cause this function to run slower as CAS numbers are used as internal
identifiers by NIST.}

\item{type}{Retention index type. One of \code{"kovats"}, \code{"linear"},
\code{"alkane"}, or \code{"lee"}. See details for more.}

\item{polarity}{Column polarity. One of \code{"polar"} or \code{"non-polar"}
to get RIs calculated for polar or non-polar columns.}

\item{temp_prog}{Temperature program. One of \code{"isothermal"},
\code{"ramp"}, or \code{"custom"}.}

\item{cas}{deprecated.  Use \code{query} instead.}

\item{verbose}{logical; should a verbose output be printed on the console?}
}
\value{
returns a tibble of literature RIs with the following columns:
\itemize{
\item{\code{CAS} is the CAS number}
\item{\code{type} is the column type, e.g. "capillary"}
\item{\code{phase} is the stationary phase (column phase)}
\item{\code{RI} is retention index}
\item{\code{length} is column length in meters}
\item{\code{gas} is the carrier gas used}
\item{\code{substrate}}
\item{\code{diameter} is the column diameter in mm}
\item{\code{thickness} is the phase thickness in µm}
\item{\code{program}. various columns depending on the value of
\code{temp_prog}}
\item{\code{reference} is where this retention index was published}
\item{\code{comment}. I believe this denotes the database these data
      were aggregated from}
}
}
\description{
This function scrapes NIST for literature retention indices
  given CAS numbers as an input.
}
\details{
The types of retention indices included in NIST include Kovats
  (\code{"kovats"}), Van den Dool and Kratz (\code{"linear"}), normal alkane
  (\code{"alkane"}), and Lee (\code{"lee"}). Details about how these are
  calculated are available on the NIST website:
  \url{https://webbook.nist.gov/chemistry/gc-ri/}
}
\note{
Copyright for NIST Standard Reference Data is governed by the Standard
Reference Data Act, \url{https://www.nist.gov/srd/public-law}.
}
\examples{
\dontrun{
myRIs <- nist_ri(c("78-70-6", "13474-59-4"), from = "cas", "linear",
"non-polar", "ramp")
}
}
\references{
NIST Mass Spectrometry Data Center, William E. Wallace, director,
  "Retention Indices" in NIST Chemistry WebBook, NIST Standard Reference
  Database Number 69, Eds. P.J. Linstrom and W.G. Mallard,
  National Institute of Standards and Technology, Gaithersburg MD, 20899,
  \doi{10.18434/T4D303}.
}
\seealso{
\code{\link{is.cas}} \code{\link{as.cas}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{is.inchikey_cs}
\alias{is.inchikey_cs}
\title{Check if input is a valid inchikey using ChemSpider API}
\usage{
is.inchikey_cs(x, verbose = getOption("verbose"))
}
\arguments{
\item{x}{character; input string}

\item{verbose}{logical; print messages during processing to console?}
}
\value{
a logical
}
\description{
Check if input is a valid inchikey using ChemSpider API
}
\examples{
\donttest{
# might fail if API is not available
is.inchikey_cs('BQJCRHHNABKAKU-KBQPJGBKSA-N')
is.inchikey_cs('BQJCRHHNABKAKU-KBQPJGBKSA')
is.inchikey_cs('BQJCRHHNABKAKU-KBQPJGBKSA-5')
is.inchikey_cs('BQJCRHHNABKAKU-KBQPJGBKSA-n')
is.inchikey_cs('BQJCRHHNABKAKU/KBQPJGBKSA/N')
is.inchikey_cs('BQJCRHHNABKAKU-KBQPJGBKXA-N')
is.inchikey_cs('BQJCRHHNABKAKU-KBQPJGBKSB-N')
}
}
\seealso{
\code{\link{is.inchikey}} for a pure-R implementation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pubchem.R
\name{pc_prop}
\alias{pc_prop}
\title{Retrieve compound properties from a pubchem CID}
\usage{
pc_prop(cid, properties = NULL, verbose = getOption("verbose"), ...)
}
\arguments{
\item{cid}{character; Pubchem ID (CID).}

\item{properties}{character vector; properties to retrieve, e.g.
c("MolecularFormula", "MolecularWeight"). If NULL (default) all available
properties are retrieved. See
\url{https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest}
for a list of all available properties.}

\item{verbose}{logical; should a verbose output be printed to the console?}

\item{...}{currently not used.}
}
\value{
a data.frame
}
\description{
Retrieve compound information from pubchem CID, see
\url{https://pubchem.ncbi.nlm.nih.gov/}
}
\note{
Please respect the Terms and Conditions of the National Library of
Medicine, \url{https://www.nlm.nih.gov/databases/download.html} the data
usage policies of National Center for Biotechnology Information,
\url{https://www.ncbi.nlm.nih.gov/home/about/policies/},
\url{https://pubchemdocs.ncbi.nlm.nih.gov/programmatic-access}, and the data
usage policies of the indicidual data sources
\url{https://pubchem.ncbi.nlm.nih.gov/sources/}.
}
\examples{
\donttest{
# might fail if API is not available
pc_prop(5564)

###
# multiple CIDS
comp <- c("Triclosan", "Aspirin")
cids <- get_cid(comp)
pc_prop(cids$cid, properties = c("MolecularFormula", "MolecularWeight",
"CanonicalSMILES"))
}
}
\references{
Wang, Y., J. Xiao, T. O. Suzek, et al. 2009 PubChem: A Public
Information System for
Analyzing Bioactivities of Small Molecules. Nucleic Acids Research 37:
623–633.

Kim, Sunghwan, Paul A. Thiessen, Evan E. Bolton, et al. 2016
PubChem Substance and Compound Databases. Nucleic Acids Research 44(D1):
D1202–D1213.

Kim, S., Thiessen, P. A., Bolton, E. E., & Bryant, S. H. (2015).
PUG-SOAP and PUG-REST: web services for programmatic access to chemical
information in PubChem. Nucleic acids research, gkv396.

Eduard Szöcs, Tamás Stirling, Eric R. Scott, Andreas Scharmüller,
Ralf B. Schäfer (2020). webchem: An R Package to Retrieve Chemical
Information from the Web. Journal of Statistical Software, 93(13).
\doi{10.18637/jss.v093.i13}.
}
\seealso{
\code{\link{get_cid}}, \code{\link{pc_sect}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{parse_mol}
\alias{parse_mol}
\title{Parse Molfile (as returned by ChemSpider) into a R-object.}
\usage{
parse_mol(string)
}
\arguments{
\item{string}{molfile as one string}
}
\value{
A list with of four entries: header (eh), counts line (cl), atom
block (ab) and bond block (bb).

header: a = number of atoms, b = number of bonds, l = number of atom lists,
f = obsolete, c = chiral flag (0=not chiral, 1 = chiral), s = number of stext
entries, x, r, p, i = obsolete, m = 999, v0 version

atom block: x, y, z = atom coordinates, a = mass difference, c= charge,
s= stereo parity, h = hydrogen count 1, b = stereo care box, v = valence,
h = h0 designator, r, i = not used, m = atom-atom mapping number,
n = inversion/retention flag, e = exact change flag

bond block:
1 = first atom, 2 = second atom, t = bond type, s = stereo type, x = not
used, r = bond typology, c = reacting center status.
}
\description{
Parse Molfile (as returned by ChemSpider) into a R-object.
}
\references{
Grabner, M., Varmuza, K., & Dehmer, M. (2012). RMol:
a toolset for transforming SD/Molfile structure information into R objects.
Source Code for Biology and Medicine, 7, 12.
\doi{10.1186/1751-0473-7-12}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cts.R
\name{cts_compinfo}
\alias{cts_compinfo}
\title{Get record details from Chemical Translation Service (CTS)}
\usage{
cts_compinfo(
  query,
  from = "inchikey",
  verbose = getOption("verbose"),
  inchikey
)
}
\arguments{
\item{query}{character; InChIkey.}

\item{from}{character; currently only accepts "inchikey".}

\item{verbose}{logical; should a verbose output be printed on the console?}

\item{inchikey}{deprecated}
}
\value{
a list of lists (for each supplied inchikey):
a list of 7. inchikey, inchicode, molweight, exactmass, formula, synonyms and
externalIds
}
\description{
Get record details from CTS, see \url{http://cts.fiehnlab.ucdavis.edu/}
}
\examples{
\donttest{
# might fail if API is not available
out <- cts_compinfo("XEFQLINVKFYRCS-UHFFFAOYSA-N")
# = Triclosan
str(out)
out[[1]][1:5]

### multiple inputs
inchikeys <- c("XEFQLINVKFYRCS-UHFFFAOYSA-N","BSYNRYMUTXBXSQ-UHFFFAOYSA-N" )
out2 <- cts_compinfo(inchikeys)
str(out2)
# a list of two
# extract molecular weight
sapply(out2, function(y) y$molweight)
}
}
\references{
Wohlgemuth, G., P. K. Haldiya, E. Willighagen, T. Kind, and O.
Fiehn 2010The Chemical Translation Service -- a Web-Based Tool to Improve
Standardization of Metabolomic Reports. Bioinformatics 26(20): 2647–2648.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{is.inchikey}
\alias{is.inchikey}
\title{Check if input is a valid inchikey}
\usage{
is.inchikey(
  x,
  type = c("format", "chemspider"),
  verbose = getOption("verbose")
)
}
\arguments{
\item{x}{character; input InChIKey}

\item{type}{character; How should be checked? Either, by format (see above)
('format') or by ChemSpider ('chemspider').}

\item{verbose}{logical; print messages during processing to console?}
}
\value{
a logical
}
\description{
This function checks if a string is a valid inchikey.
Inchikey must fulfill the following criteria:
1) consist of 27 characters;
2) be all uppercase, all letters (no numbers);
3) contain two hyphens at positions 15 and 26;
4) 24th character (flag character) be 'S' (Standard InChI) or 'N' (non-standard)
5) 25th character (version character) must be 'A' (currently).
}
\note{
This function can handle only one inchikey string.
}
\examples{
is.inchikey('BQJCRHHNABKAKU-KBQPJGBKSA-N')
is.inchikey('BQJCRHHNABKAKU-KBQPJGBKSA')
is.inchikey('BQJCRHHNABKAKU-KBQPJGBKSA-5')
is.inchikey('BQJCRHHNABKAKU-KBQPJGBKSA-n')
is.inchikey('BQJCRHHNABKAKU/KBQPJGBKSA/N')
is.inchikey('BQJCRHHNABKAKU-KBQPJGBKXA-N')
is.inchikey('BQJCRHHNABKAKU-KBQPJGBKSB-N')
}
\references{
Heller, Stephen R., et al. "InChI, the IUPAC International
Chemical Identifier." Journal of Cheminformatics 7.1 (2015): 23.

Eduard Szöcs, Tamás Stirling, Eric R. Scott, Andreas Scharmüller,
Ralf B. Schäfer (2020). webchem: An R Package to Retrieve Chemical
Information from the Web. Journal of Statistical Software, 93(13).
\doi{10.18637/jss.v093.i13}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integration.R
\name{with_cts}
\alias{with_cts}
\title{Auto-translate identifiers and search databases}
\usage{
with_cts(query, from, .f, .verbose = getOption("verbose"), ...)
}
\arguments{
\item{query}{character; the search term}

\item{from}{character; the format or type of query.  Commonly accepted values
are "name", "cas", "inchi", and "inchikey"}

\item{.f}{character; the (quoted) name of a webchem function}

\item{.verbose}{logical; print a message when translating query?}

\item{...}{other arguments passed to the function specified with \code{.f}}
}
\value{
returns results from \code{.f}
}
\description{
Supply a query of any type (e.g. SMILES, CAS, name, InChI, etc.) along with
any webchem function that has \code{query} and \code{from} arguments.  If the
function doesn't accept the type of query you've supplied, this will try to
automatically translate it using CTS and run the query.
}
\note{
During the translation step, only the first hit from CTS is used.
  Therefore, using this function to translate on the fly is not foolproof and
  care should be taken to verify the results.
}
\examples{
\dontrun{
with_cts("XDDAORKBJWWYJS-UHFFFAOYSA-N", from = "inchikey", .f = "get_etoxid")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wikidata.R
\name{get_wdid}
\alias{get_wdid}
\title{Get Wikidata Item ID}
\usage{
get_wdid(
  query,
  match = c("best", "first", "all", "ask", "na"),
  verbose = getOption("verbose"),
  language = "en"
)
}
\arguments{
\item{query}{character; The searchterm}

\item{match}{character; How should multiple hits be handeled? 'all' returns
all matched IDs, 'first' only the first match, 'best' the best matching (by
name) ID, 'ask' is a interactive mode and the user is asked for input, na'
returns NA if multiple hits are found.}

\item{verbose}{logical; print message during processing to console?}

\item{language}{character; the language to search in}
}
\value{
if match = 'all' a list with ids, otherwise a dataframe with 4 columns:
id, matched text, string distance to match and the queried string
}
\description{
Search www.wikidata.org for wikidata item identifiers.  Note that this search
is currently not limited to chemical substances, so be sure to check your
results.
}
\note{
Only matches in labels are returned.
}
\examples{
\dontrun{
get_wdid('Triclosan', language = 'de')
get_wdid('DDT')
get_wdid('DDT', match = 'all')

# multiple inputs
comps <- c('Triclosan', 'Glyphosate')
get_wdid(comps)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pubchem.R
\name{pc_synonyms}
\alias{pc_synonyms}
\title{Search synonyms in pubchem}
\usage{
pc_synonyms(
  query,
  from = c("name", "cid", "sid", "aid", "smiles", "inchi", "inchikey"),
  match = c("all", "first", "ask", "na"),
  verbose = getOption("verbose"),
  arg = NULL,
  choices = NULL,
  ...
)
}
\arguments{
\item{query}{character; search term.}

\item{from}{character; type of input, can be one of "name" (default), "cid",
"sid", "aid", "smiles", "inchi", "inchikey"}

\item{match}{character; How should multiple hits be handled? \code{"all"}
returns all matches, \code{"first"} returns only the first result,
\code{"ask"} enters an interactive mode and the user is asked for input,
\code{"na"} returns \code{NA} if multiple hits are found.}

\item{verbose}{logical; should a verbose output be printed on the console?}

\item{arg}{character; optional arguments like "name_type=word" to match
individual words.}

\item{choices}{deprecated.  Use the \code{match} argument instead.}

\item{...}{currently unused}
}
\value{
a named list.
}
\description{
Search synonyms using PUG-REST,
see \url{https://pubchem.ncbi.nlm.nih.gov/}.
}
\note{
Please respect the Terms and Conditions of the National Library of
Medicine, \url{https://www.nlm.nih.gov/databases/download.html} the data
usage policies of National Center for Biotechnology Information,
\url{https://www.ncbi.nlm.nih.gov/home/about/policies/},
\url{https://pubchemdocs.ncbi.nlm.nih.gov/programmatic-access}, and the data
usage policies of the indicidual data sources
\url{https://pubchem.ncbi.nlm.nih.gov/sources/}.
}
\examples{
\donttest{
pc_synonyms("Aspirin")
pc_synonyms(c("Aspirin", "Triclosan"))
pc_synonyms(5564, from = "cid")
pc_synonyms(c("Aspirin", "Triclosan"), match = "ask")
}
}
\references{
Wang, Y., J. Xiao, T. O. Suzek, et al. 2009 PubChem: A Public
Information System for
Analyzing Bioactivities of Small Molecules. Nucleic Acids Research 37:
623–633.

Kim, Sunghwan, Paul A. Thiessen, Evan E. Bolton, et al. 2016
PubChem Substance and Compound Databases. Nucleic Acids Research 44(D1):
D1202–D1213.

Kim, S., Thiessen, P. A., Bolton, E. E., & Bryant, S. H. (2015).
PUG-SOAP and PUG-REST: web services for programmatic access to chemical
information in PubChem. Nucleic acids research, gkv396.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/webchem-defunct.R
\name{webchem-defunct}
\alias{webchem-defunct}
\alias{ppdb_query}
\alias{ppdb_parse}
\alias{ppdb}
\alias{cir}
\alias{pp_query}
\alias{cs_prop}
\title{Defunct function(s) in the webchem package}
\usage{
ppdb_query()

ppdb_parse()

ppdb()

cir()

pp_query()

cs_prop()
}
\description{
These functions are defunct and no longer available.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/webchem-package.R
\docType{data}
\name{lc50}
\alias{lc50}
\title{Acute toxicity data from U.S. EPA ECOTOX}
\format{
A data frame with 124 rows and 2 variables:
\describe{
  \item{cas}{CAS registry number}
  \item{value}{LC50value}
}
}
\source{
\url{https://cfpub.epa.gov/ecotox/}
}
\usage{
lc50
}
\description{
This dataset comprises acute ecotoxicity data of 124 insecticides.
The data is publicly available and can be retrieved from the EPA ECOTOX database
(\url{https://cfpub.epa.gov/ecotox/})
It comprises acute toxicity data (D. magna, 48h, Laboratory, 48h) and has been
preprocessed (remove non-insecticides, aggregate multiple value, keep only numeric data etc).
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cts.R
\name{cts_from}
\alias{cts_from}
\title{Return a list of all possible ids}
\usage{
cts_from(verbose = getOption("verbose"))
}
\arguments{
\item{verbose}{logical; should a verbose output be printed on the console?}
}
\value{
a character vector.
}
\description{
Return a list of all possible ids that can be used in the 'from' argument
}
\details{
See also \url{http://cts.fiehnlab.ucdavis.edu/services}
}
\examples{
\donttest{
cts_from()
}
}
\references{
Wohlgemuth, G., P. K. Haldiya, E. Willighagen, T. Kind, and O.
Fiehn 2010The Chemical Translation Service -- a Web-Based Tool to Improve
Standardization of Metabolomic Reports. Bioinformatics 26(20): 2647–2648.
}
\seealso{
\code{\link{cts_convert}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/etox.R
\name{etox_targets}
\alias{etox_targets}
\title{Get Quality Targets from a ETOX ID}
\usage{
etox_targets(id, verbose = getOption("verbose"))
}
\arguments{
\item{id}{character; ETOX ID}

\item{verbose}{logical; print message during processing to console?}
}
\value{
A list of lists of two: \code{res} a data.frame with quality targets
  from the ETOX database, and source_url.
}
\description{
Query ETOX: Information System Ecotoxicology and Environmental Quality
Targets \url{https://webetox.uba.de/webETOX/index.do} for quality targets
}
\note{
Before using this function, please read the disclaimer
  \url{https://webetox.uba.de/webETOX/disclaimer.do}.
}
\examples{
\dontrun{
id <- get_etoxid('Triclosan', match = 'best')
out <- etox_targets(id$etoxid)
out[ , c('Substance', 'CAS_NO', 'Country_or_Region', 'Designation',
'Value_Target_LR', 'Unit')]
etox_targets( c("20179", "9051"))

}
}
\references{
Eduard Szöcs, Tamás Stirling, Eric R. Scott, Andreas Scharmüller,
Ralf B. Schäfer (2020). webchem: An R Package to Retrieve Chemical
Information from the Web. Journal of Statistical Software, 93(13).
\doi{10.18637/jss.v093.i13}.
}
\seealso{
\code{\link{get_etoxid}} to retrieve ETOX IDs,
  \code{\link{etox_basic}} for basic information, \code{\link{etox_targets}}
  for quality targets and \code{\link{etox_tests}} for test results
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/srs.R
\name{srs_query}
\alias{srs_query}
\title{Get record details from U.S. EPA Substance Registry Servives (SRS)}
\usage{
srs_query(
  query,
  from = c("itn", "cas", "epaid", "tsn", "name"),
  verbose = getOption("verbose"),
  ...
)
}
\arguments{
\item{query}{character; query ID.}

\item{from}{character; type of query ID, e.g. \code{'itn'} , \code{'cas'},
\code{'epaid'}, \code{'tsn'}, \code{'name'}.}

\item{verbose}{logical; should a verbose output be printed on the console?}

\item{...}{not currently used.}
}
\value{
a list of lists (for each supplied query): a list of 22. subsKey,
 internalTrackingNumber, systematicName, epaIdentificationNumber,
 currentCasNumber, currentTaxonomicSerialNumber, epaName, substanceType,
 categoryClass, kingdomCode, iupacName, pubChemId, molecularWeight,
 molecularFormula, inchiNotation, smilesNotation, classifications,
 characteristics, synonyms, casNumbers, taxonomicSerialNumbers, relationships
}
\description{
Get record details from SRS, see \url{https://cdxnodengn.epa.gov/cdx-srs-rest/}
}
\examples{
\donttest{
# might fail if API is not available
srs_query(query = '50-00-0', from = 'cas')

### multiple inputs
casrn <- c('50-00-0', '67-64-1')
srs_query(query = casrn, from = 'cas')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ping.R
\name{ping_service}
\alias{ping_service}
\title{Ping an API used in webchem to see if it's working.}
\usage{
ping_service(
  service = c("bcpc", "chebi", "ci", "cs", "cs_web", "cir", "cts", "etox", "fn",
    "nist", "opsin", "pan", "pc", "srs", "wd")
)
}
\arguments{
\item{service}{character; the same abbreviations used as prefixes in \code{webchem} functions, with the exception of \code{"cs_web"}, which only checks if the ChemSpider website is up, and thus doesn't require an API key.}
}
\value{
A logical, TRUE if the service is available or FALSE if it isn't
}
\description{
Ping an API used in webchem to see if it's working.
}
\examples{
\dontrun{
ping_service("pan")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/webchem-package.R
\docType{data}
\name{jagst}
\alias{jagst}
\title{Organic plant protection products in the river Jagst / Germany in 2013}
\format{
A data frame with 442 rows and 4 variables:
\describe{
  \item{date}{sampling data}
  \item{substance}{substance names}
  \item{value}{concentration in ug/L}
  \item{qual}{qualifier, indicating values < LOQ}
}
}
\source{
\url{https://udo.lubw.baden-wuerttemberg.de/?highlightglobalid=gewaesserguetedaten}
}
\usage{
jagst
}
\description{
This dataset comprises environmental monitoring data of organic plant protection products
in the year 2013 in the river Jagst, Germany.
The data is publicly available and can be retrieved from the
LUBW Landesanstalt für Umwelt, Messungen und Naturschutz Baden-Württemberg.
It has been preprocessed and comprises measurements of 34 substances.
Substances without detects have been removed.
on 13 sampling occasions.
Values are given in ug/L.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/webchem-package.R
\docType{package}
\name{webchem}
\alias{webchem}
\title{webchem: An R package to retrieve chemical information from the web.}
\description{
Chemical information from around the web. This package interacts with a suite
of web APIs for chemical information.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{write_mol}
\alias{write_mol}
\title{Export a Chemical Structure in .mol Format.}
\usage{
write_mol(x, file = "")
}
\arguments{
\item{x}{a character string of a chemical structure in mol format.}

\item{file}{a character vector of file names}
}
\description{
Some webchem functions return character strings that contain a chemical
structure in Mol format. This function exports a character string as a .mol
file so it can be imported with other chemistry software.
}
\examples{
\dontrun{
# export Mol file
csid <- get_csid("bergapten")
mol3d <- cs_compinfo(csid$csid, field = "Mol3D")
write_mol(mol3d$mol3D, file = mol3d$id)

# export multiple Mol files
csids <- get_csid(c("bergapten", "xanthotoxin"))
mol3ds <- cs_compinfo(csids$csid, field = "Mol3D")
mapply(function(x, y) write_mol(x, y), x = mol3ds$mol3D, y = mol3ds$id)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flavornet.R
\name{fn_percept}
\alias{fn_percept}
\title{Retrieve flavor percepts from www.flavornet.org}
\usage{
fn_percept(query, from = "cas", verbose = getOption("verbose"), CAS, ...)
}
\arguments{
\item{query}{character; CAS number to search by. See \code{\link{is.cas}} for correct formatting}

\item{from}{character; currently only CAS numbers are accepted.}

\item{verbose}{logical; should a verbose output be printed on the console?}

\item{CAS}{deprecated}

\item{...}{currently unused}
}
\value{
A named character vector containing flavor percepts or NA's in the case of CAS numbers that are not found
}
\description{
Retreive flavor percepts from \url{http://www.flavornet.org}.  Flavornet is a database of 738 compounds with odors
perceptible to humans detected using gas chromatography olfactometry (GCO).
}
\examples{
\dontrun{
# might fail if website is not available
fn_percept("123-32-0")

CASs <- c("75-07-0",  "64-17-5",  "109-66-0", "78-94-4",  "78-93-3")
fn_percept(CASs)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcpc.R
\name{bcpc_query}
\alias{bcpc_query}
\title{Query https://pesticidecompendium.bcpc.org}
\usage{
bcpc_query(
  query,
  from = c("name", "cas"),
  verbose = getOption("verbose"),
  type,
  ...
)
}
\arguments{
\item{query}{character; search string}

\item{from}{character; type of input ('cas' or 'name')}

\item{verbose}{logical; print message during processing to console?}

\item{type}{deprecated}

\item{...}{additional arguments to internal utility functions}
}
\value{
A list of eight entries: common-name, status, preferred IUPAC Name,
IUPAC Name, cas, formula, activity, subactivity, inchikey, inchi and source
url.
}
\description{
Query the BCPC Compendium of Pesticide Common Names
\url{https://pesticidecompendium.bcpc.org}
formerly known as Alan Woods Compendium of Pesticide Common Names
}
\note{
for from = 'cas' only the first matched link is returned.
Please respect Copyright, Terms and Conditions
\url{https://pesticidecompendium.bcpc.org/legal.html}!
}
\examples{
\dontrun{
bcpc_query('Fluazinam', from = 'name')
out <- bcpc_query(c('Fluazinam', 'Diclofop'), from = 'name')
out
# extract subactivity from object
sapply(out, function(y) y$subactivity[1])

# use CAS-numbers
bcpc_query("79622-59-6", from = 'cas')
}
}
\references{
Eduard Szöcs, Tamás Stirling, Eric R. Scott, Andreas Scharmüller,
Ralf B. Schäfer (2020). webchem: An R Package to Retrieve Chemical
Information from the Web. Journal of Statistical Software, 93(13).
\doi{10.18637/jss.v093.i13}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cir.R
\name{cir_query}
\alias{cir_query}
\title{Query Chemical Identifier Resolver}
\usage{
cir_query(
  identifier,
  representation = "smiles",
  resolver = NULL,
  match = c("all", "first", "ask", "na"),
  verbose = getOption("verbose"),
  choices = NULL,
  ...
)
}
\arguments{
\item{identifier}{character; chemical identifier.}

\item{representation}{character; what representation of the identifier should
be returned. See details for possible representations.}

\item{resolver}{character; what resolver should be used? If NULL (default)
the identifier type is detected and the different resolvers are used in turn.
See details for possible resolvers.}

\item{match}{character; How should multiple hits be handled? \code{"all"}
returns all matches, \code{"first"} returns only the first result,
\code{"ask"} enters an interactive mode and the user is asked for input,
\code{"na"} returns \code{NA} if multiple hits are found.}

\item{verbose}{logical; should a verbose output be printed on the console?}

\item{choices}{deprecated.  Use the \code{match} argument instead.}

\item{...}{currently not used.}
}
\value{
A list of character vectors.
}
\description{
A interface to the Chemical Identifier Resolver (CIR).
 (\url{https://cactus.nci.nih.gov/chemical/structure_documentation}).
}
\details{
CIR can resolve can be of the following \code{identifier}: Chemical Names,
 IUPAC names,
 CAS Numbers, SMILES strings, IUPAC InChI/InChIKeys, NCI/CADD Identifiers,
 CACTVS HASHISY, NSC number, PubChem SID, ZINC Code, ChemSpider ID,
 ChemNavigator SID, eMolecule VID.

 \code{cir_query()} can handle only a part of all possible conversions of CIR.
 Possible \code{representations} are:
 \itemize{
     \item \code{'smiles'}(SMILES strings),
     \item \code{'names'} (Names),
     \item \code{'cas'} (CAS numbers),
     \item \code{'stdinchikey'} (Standard InChIKey),
     \item \code{'stdinchi'} (Standard InChI),
     \item \code{'ficts'} (FICTS Identifier),
     \item \code{'ficus'} (FICuS Indetifier),
     \item \code{'uuuuu'} (uuuuu Identifier),
     \item \code{'mw'} (Molecular weight),
     \item \code{'monoisotopic_mass'} (Monoisotopic Mass),
     \item \code{'formula'} (Chemical Formula),
     \item \code{'chemspider_id'} (ChemSpider ID),
     \item \code{'pubchem_sid'} (PubChem SID),
     \item \code{'chemnavigator_sid'} (ChemNavigator SID),
     \item \code{'h_bond_donor_count'} (Number of Hydrogen Bond Donors),
     \item \code{'h_bond_acceptor_count'} (Number of Hydrogen Bond Acceptors),
     \item \code{'h_bond_center_count'} (Number of Hydrogen Bond Centers),
     \item \code{'rule_of_5_violation_count'} (Number of Rule of 5 Violations),
     \item \code{'rotor_count'} (Number of Freely Rotatable Bonds),
     \item \code{'effective_rotor_count'} (Number of Effectively Rotatable Bonds),
     \item \code{'ring_count'} (Number of Rings),
     \item \code{'ringsys_count'} (Number of Ring Systems),
     \item \code{'xlogp2'} (octanol-water partition coefficient),
     \item \code{'aromatic'} (is the compound aromatic),
     \item \code{'macrocyclic'} (is the compound macrocyclic),
     \item \code{'heteroatom_count'} (heteroatom count),
     \item \code{'hydrogen_atom_count'} (H atom count),
     \item \code{'heavy_atom_count'} ( Heavy atom count),
     \item \code{'deprotonable_group_count'} (Number of deprotonable groups),
     \item \code{'protonable_group_count'} (Number of protonable groups).
 }

 CIR first tries to determine the identifier type submitted and then
 uses 'resolvers' to look up the data.
 If no \code{resolver} is supplied, CIR tries different resolvers in
 turn till a hit is found.
 E.g. for names CIR tries first to look up in OPSIN and if this fails
 the local name index of CIR.
 However, it can be also specified which resolvers to use
 (if you know e.g. know your identifier type)
 Possible \code{resolvers} are:
 \itemize{
   \item \code{'name_by_cir'} (Lookup in name index of CIR),
   \item \code{'name_by_opsin'} (Lookup in OPSIN),
   \item \code{'name_by_chemspider'} (Lookup in ChemSpider,
   \url{https://cactus.nci.nih.gov/blog/?p=1386}),
   \item \code{'smiles'} (Lookup SMILES),
   \item \code{'stdinchikey'}, \code{'stdinchi'} (InChI),
   \item \code{'cas_number'} (CAS Number),
   \item \code{'name_pattern'} (Google-like pattern search
   (\url{https://cactus.nci.nih.gov/blog/?p=1456})
   Note, that the pattern search can be combined with other resolvers,
   e.g. \code{resolver = 'name_by_chemspider,name_pattern'}.

 }
}
\note{
You can only make 1 request per second (this is a hard-coded feature).
}
\examples{
\donttest{
# might fail if API is not available
cir_query("Triclosan", "cas")
cir_query("3380-34-5", "cas", match = "first")
cir_query("3380-34-5", "cas", resolver = "cas_number")
cir_query("3380-34-5", "smiles")
cir_query("Triclosan", "mw")

# multiple inputs
comp <- c("Triclosan", "Aspirin")
cir_query(comp, "cas", match = "first")

}
}
\references{
\code{cir} relies on the great CIR web service created by the CADD
Group at NCI/NIH! \cr
\url{https://cactus.nci.nih.gov/chemical/structure_documentation}, \cr
\url{https://cactus.nci.nih.gov/blog/?cat=10}, \cr
\url{https://cactus.nci.nih.gov/blog/?p=1386}, \cr
\url{https://cactus.nci.nih.gov/blog/?p=1456}, \cr
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{is.inchikey_format}
\alias{is.inchikey_format}
\title{Check if input is a valid inchikey using format}
\usage{
is.inchikey_format(x, verbose = getOption("verbose"))
}
\arguments{
\item{x}{character; input string}

\item{verbose}{logical; print messages during processing to console?}
}
\value{
a logical
}
\description{
Inchikey must fulfill the following criteria:
1) consist of 27 characters;
2) be all uppercase, all letters (no numbers);
3) contain two hyphens at positions 15 and 26;
4) 24th character (flag character) be 'S' (Standard InChI) or 'N'
(non-standard)
5) 25th character (version character) must be 'A' (currently).
}
\examples{
\donttest{
# might fail if API is not available
is.inchikey_format('BQJCRHHNABKAKU-KBQPJGBKSA-N')
is.inchikey_format('BQJCRHHNABKAKU-KBQPJGBKSA')
is.inchikey_format('BQJCRHHNABKAKU-KBQPJGBKSA-5')
is.inchikey_format('BQJCRHHNABKAKU-KBQPJGBKSA-n')
is.inchikey_format('BQJCRHHNABKAKU/KBQPJGBKSA/N')
is.inchikey_format('BQJCRHHNABKAKU-KBQPJGBKXA-N')
is.inchikey_format('BQJCRHHNABKAKU-KBQPJGBKSB-N')
}
}
\seealso{
\code{\link{is.inchikey}} for a pure-R implementation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/webchem-deprecated.R
\name{webchem-deprecated}
\alias{webchem-deprecated}
\alias{cid_compinfo}
\alias{aw_query}
\title{Deprecated function(s) in the webchem package}
\usage{
cid_compinfo(...)

aw_query(...)
}
\arguments{
\item{...}{Parameters to be passed to the modern version of the function}
}
\description{
These functions are provided for compatibility with older version of
the webchem package.  They may eventually be completely
removed.
}
\details{
Deprecated functions are:
\tabular{rl}{
  \code{pc_prop} \tab was formerly \code{\link{cid_compinfo}}\cr
  \code{bcpc_query} \tab was formerly \code{\link{aw_query}}\cr
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/etox.R
\name{etox_tests}
\alias{etox_tests}
\title{Get Tests from a ETOX ID}
\usage{
etox_tests(id, verbose = getOption("verbose"))
}
\arguments{
\item{id}{character; ETOX ID}

\item{verbose}{logical; print message during processing to console?}
}
\value{
A list of lists of two: A data.frame with test results from the ETOX database and the source_url.
}
\description{
Query ETOX: Information System Ecotoxicology and Environmental Quality Targets
\url{https://webetox.uba.de/webETOX/index.do} for tests
}
\note{
Before using this function, please read the disclaimer
\url{https://webetox.uba.de/webETOX/disclaimer.do}.
}
\examples{
\dontrun{
id <- get_etoxid('Triclosan', match = 'best')
out <- etox_tests(id$etoxid)
out[ , c('Organism', 'Effect', 'Duration', 'Time_Unit',
'Endpoint', 'Value', 'Unit')]
etox_tests( c("20179", "9051"))
}
}
\seealso{
\code{\link{get_etoxid}} to retrieve ETOX IDs, \code{\link{etox_basic}} for basic information,
\code{\link{etox_targets}} for quality targets and \code{\link{etox_tests}} for test results
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pubchem.R
\name{get_cid}
\alias{get_cid}
\title{Retrieve Pubchem Compound ID (CID)}
\usage{
get_cid(
  query,
  from = "name",
  domain = c("compound", "substance", "assay"),
  match = c("all", "first", "ask", "na"),
  verbose = getOption("verbose"),
  arg = NULL,
  first = NULL,
  ...
)
}
\arguments{
\item{query}{character; search term, one or more compounds.}

\item{from}{character; type of input. See details for more information.}

\item{domain}{character; query domain, can be one of \code{"compound"},
\code{"substance"}, \code{"assay"}.}

\item{match}{character; How should multiple hits be handled?, \code{"all"}
all matches are returned, \code{"first"} the first matching is returned,
\code{"ask"} enters an interactive mode and the user is asked for input,
\code{"na"} returns NA if multiple hits are found.}

\item{verbose}{logical; should a verbose output be printed on the console?}

\item{arg}{character; optinal arguments like "name_type=word" to match
individual words.}

\item{first}{deprecated. Use `match` instead.}

\item{...}{currently unused.}
}
\value{
a tibble.
}
\description{
Retrieve compound IDs (CIDs) from PubChem.
}
\details{
Valid values for the \code{from} argument depend on the
\code{domain}:
\itemize{
\item{\code{compound}: \code{"name"}, \code{"smiles"}, \code{"inchi"},
\code{"inchikey"}, \code{"formula"}, \code{"sdf"}, <xref>,
<structure search>, <fast search>.}
\item{\code{substance}: \code{"name"}, \code{"sid"},
\code{<xref>}, \code{"sourceid/<source id>"} or \code{"sourceall"}.}
\item{\code{assay}: \code{"aid"}, \code{<assay target>}.}
}

<structure search> is assembled as "{\code{substructure} |
\code{superstructure} | \code{similarity} | \code{identity}} / {\code{smiles}
 | \code{inchi} | \code{sdf} | \code{cid}}", e.g.
 \code{from = "substructure/smiles"}.

\code{<xref>} is assembled as "\code{xref}/\{\code{RegistryID} |
\code{RN} | \code{PubMedID} | \code{MMDBID} | \code{ProteinGI},
\code{NucleotideGI} | \code{TaxonomyID} | \code{MIMID} | \code{GeneID} |
\code{ProbeID} | \code{PatentID}\}", e.g. \code{from = "xref/RN"} will query
by CAS RN.

<fast search> is either \code{fastformula} or it is assembled as
"{\code{fastidentity} | \code{fastsimilarity_2d} | \code{fastsimilarity_3d} |
\code{fastsubstructure} | \code{fastsuperstructure}}/{\code{smiles} |
\code{smarts} | \code{inchi} | \code{sdf} | \code{cid}}", e.g.
\code{from = "fastidentity/smiles"}.

\code{<source id>} is any valid PubChem Data Source ID. When
\code{from = "sourceid/<source id>"}, the query is the ID of the substance in
the depositor's database.

If \code{from = "sourceall"} the query is one or more valid Pubchem
depositor names. Depositor names are not case sensitive.

Depositor names and Data Source IDs can be found at
\url{https://pubchem.ncbi.nlm.nih.gov/sources/}.

\code{<assay target>} is assembled as "\code{target}/\{\code{gi} |
\code{proteinname} | \code{geneid} | \code{genesymbol} | \code{accession}\}",
e.g. \code{from = "target/geneid"} will query by GeneID.
}
\note{
Please respect the Terms and Conditions of the National Library of
Medicine, \url{https://www.nlm.nih.gov/databases/download.html} the data
usage policies of National Center for Biotechnology Information,
\url{https://www.ncbi.nlm.nih.gov/home/about/policies/},
\url{https://pubchemdocs.ncbi.nlm.nih.gov/programmatic-access}, and the data
usage policies of the indicidual data sources
\url{https://pubchem.ncbi.nlm.nih.gov/sources/}.
}
\examples{
\donttest{
# might fail if API is not available
get_cid("Triclosan")
get_cid("Triclosan", arg = "name_type=word")
# from SMILES
get_cid("CCCC", from = "smiles")
# from InChI
get_cid("InChI=1S/CH5N/c1-2/h2H2,1H3", from = "inchi")
# from InChIKey
get_cid("BPGDAMSIGCZZLK-UHFFFAOYSA-N", from = "inchikey")
# from formula
get_cid("C26H52NO6P", from = "formula")
# from CAS RN
get_cid("56-40-6", from = "xref/rn")
# similarity
get_cid(5564, from = "similarity/cid")
get_cid("CCO", from = "similarity/smiles")
# from SID
get_cid("126534046", from = "sid", domain = "substance")
# sourceid
get_cid("VCC957895", from = "sourceid/23706", domain = "substance")
# sourceall
get_cid("Optopharma Ltd", from = "sourceall", domain = "substance")
# from AID (CIDs of substances tested in the assay)
get_cid(170004, from = "aid", domain = "assay")
# from GeneID (CIDs of substances tested on the gene)
get_cid(25086, from = "target/geneid", domain = "assay")

# multiple inputs
get_cid(c("Triclosan", "Aspirin"))

}
}
\references{
Wang, Y., J. Xiao, T. O. Suzek, et al. 2009 PubChem: A Public
Information System for
Analyzing Bioactivities of Small Molecules. Nucleic Acids Research 37:
623–633.

Kim, Sunghwan, Paul A. Thiessen, Evan E. Bolton, et al. 2016
PubChem Substance and Compound Databases. Nucleic Acids Research 44(D1):
D1202–D1213.

Kim, S., Thiessen, P. A., Bolton, E. E., & Bryant, S. H. (2015).
PUG-SOAP and PUG-REST: web services for programmatic access to chemical
information in PubChem. Nucleic acids research, gkv396.

Eduard Szöcs, Tamás Stirling, Eric R. Scott, Andreas Scharmüller,
Ralf B. Schäfer (2020). webchem: An R Package to Retrieve Chemical
Information from the Web. Journal of Statistical Software, 93(13).
\doi{10.18637/jss.v093.i13}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cir.R
\name{cir_img}
\alias{cir_img}
\title{Query Chemical Identifier Resolver Images}
\usage{
cir_img(
  query,
  dir,
  format = c("png", "gif"),
  width = 500,
  height = 500,
  linewidth = 2,
  symbolfontsize = 16,
  bgcolor = NULL,
  antialiasing = TRUE,
  atomcolor = NULL,
  bondcolor = NULL,
  csymbol = c("special", "all"),
  hsymbol = c("special", "all"),
  hcolor = NULL,
  header = NULL,
  footer = NULL,
  frame = NULL,
  verbose = getOption("verbose"),
  ...
)
}
\arguments{
\item{query}{character; Search term. Can be any common chemical identifier
(e.g. CAS, INCHI(KEY), SMILES etc.)}

\item{dir}{character; Directory to save the image.}

\item{format}{character; Output format of the image. Can be one of "png",
"gif".}

\item{width}{integer; Width of the image.}

\item{height}{integer; Height of the image.}

\item{linewidth}{integer; Width of lines.}

\item{symbolfontsize}{integer; Fontsize of atoms in the image.}

\item{bgcolor}{character; E.g. transparent, white, \%23AADDEE}

\item{antialiasing}{logical; Should antialiasing be used?}

\item{atomcolor}{character; Color of the atoms in the image.}

\item{bondcolor}{character; Color of the atom bond lines.}

\item{csymbol}{character; Can be one of "special" (default - i.e. only
hydrogen atoms in functional groups or defining stereochemistry) or "all".}

\item{hsymbol}{character; Can be one of "special" (default - i.e. none are
shown) or "all" (all are printed).}

\item{hcolor}{character; Color of the hydrogen atoms.}

\item{header}{character; Should a header text be added to the image? Can be
any string.}

\item{footer}{character; Should a footer text be added to the image? Can be
any string.}

\item{frame}{integer; Should a frame be plotted? Can be on of NULL (default)
or 1.}

\item{verbose}{logical; Should a verbose output be printed on the console?}

\item{...}{currently not used.}
}
\value{
image written to disk
}
\description{
A interface to the Chemical Identifier Resolver (CIR).
 (\url{https://cactus.nci.nih.gov/chemical/structure_documentation}).
}
\details{
CIR can resolve can be of the following \code{identifier}: Chemical Names,
 IUPAC names,
 CAS Numbers, SMILES strings, IUPAC InChI/InChIKeys, NCI/CADD Identifiers,
 CACTVS HASHISY, NSC number, PubChem SID, ZINC Code, ChemSpider ID,
 ChemNavigator SID, eMolecule VID.

 For an image with transparent background use ‘transparent’ as color name and
 switch off antialiasing (i.e. antialiasing = 0).
}
\note{
You can only make 1 request per second (this is a hard-coded feature).
}
\examples{
\donttest{
# might fail if API is not available
cir_img("CCO", dir = tempdir()) # SMILES

# multiple query strings and different formats
query = c("Glyphosate", "Isoproturon", "BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
cir_img(query, dir = tempdir(), bgcolor = "transparent", antialising = 0)

# all parameters
query  = "Triclosan"
cir_img(query,
        dir = tempdir(),
        format = "png",
        width = 600,
        height = 600,
        linewidth = 5,
        symbolfontsize = 30,
        bgcolor = "red",
        antialiasing = FALSE,
        atomcolor = "green",
        bondcolor = "yellow",
        csymbol = "all",
        hsymbol = "all",
        hcolor = "purple",
        header = "My funky chemical structure..",
        footer = "..is just so awesome!",
        frame = 1,
        verbose = getOption("verbose"))
}
}
\references{
\code{cir} relies on the great CIR web service created by the CADD
Group at NCI/NIH! \cr
\url{https://cactus.nci.nih.gov/chemical/structure_documentation}, \cr
\url{https://cactus.nci.nih.gov/blog/?cat=10}, \cr
\url{https://cactus.nci.nih.gov/blog/?p=1386}, \cr
\url{https://cactus.nci.nih.gov/blog/?p=1456}, \cr
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cts.R
\name{cts_convert}
\alias{cts_convert}
\title{Convert Ids using Chemical Translation Service (CTS)}
\usage{
cts_convert(
  query,
  from,
  to,
  match = c("all", "first", "ask", "na"),
  verbose = getOption("verbose"),
  choices = NULL,
  ...
)
}
\arguments{
\item{query}{character; query ID.}

\item{from}{character; type of query ID, e.g. \code{'Chemical Name'} ,
\code{'InChIKey'}, \code{'PubChem CID'}, \code{'ChemSpider'}, \code{'CAS'}.}

\item{to}{character; type to convert to.}

\item{match}{character; How should multiple hits be handled? \code{"all"}
returns all matches, \code{"first"} returns only the first result,
\code{"ask"} enters an interactive mode and the user is asked for input,
\code{"na"} returns \code{NA} if multiple hits are found.}

\item{verbose}{logical; should a verbose output be printed on the console?}

\item{choices}{deprecated.  Use the \code{match} argument instead.}

\item{...}{currently not used.}
}
\value{
a list of character vectors or if \code{choices} is used, then a
single named vector.
}
\description{
Convert Ids using Chemical Translation Service (CTS), see
\url{http://cts.fiehnlab.ucdavis.edu/}
}
\details{
See also \url{http://cts.fiehnlab.ucdavis.edu/}
for possible values of from and to.
}
\note{
When this version of webchem was released, CTS was temporarily unable
to convert chemical names to IDs.
}
\examples{
\donttest{
# might fail if API is not available
cts_convert("XEFQLINVKFYRCS-UHFFFAOYSA-N", "inchikey", "Chemical Name")

### multiple inputs
keys <- c("XEFQLINVKFYRCS-UHFFFAOYSA-N", "VLKZOEOYAKHREP-UHFFFAOYSA-N")
cts_convert(keys, "inchikey", "cas")
}
}
\references{
Wohlgemuth, G., P. K. Haldiya, E. Willighagen, T. Kind, and O.
Fiehn 2010The Chemical Translation Service -- a Web-Based Tool to Improve
Standardization of Metabolomic Reports. Bioinformatics 26(20): 2647–2648.
}
\seealso{
\code{\link{cts_from}} for possible values in the 'from' argument
and \code{\link{cts_to}} for possible values in the 'to' argument.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chemspider.R
\name{cs_convert}
\alias{cs_convert}
\title{Convert identifiers using ChemSpider}
\usage{
cs_convert(query, from, to, verbose = getOption("verbose"), apikey = NULL)
}
\arguments{
\item{query}{character; query ID.}

\item{from}{character; type of query ID.}

\item{to}{character; type to convert to.}

\item{verbose}{logical; should a verbose output be printed on the console?}

\item{apikey}{character; your API key. If NULL (default),
\code{cs_check_key()} will look for it in .Renviron or .Rprofile.}
}
\value{
Returns a vector containing the converted identifier(s).
}
\description{
Submit one or more identifiers (CSID, SMILES, InChI, InChIKey or Mol) and
return one or more identifiers in another format (CSID, SMILES, InChI,
InChIKey or Mol).
}
\details{
Not all conversions are supported. Allowed conversions:
\itemize{
\item CSID <-> InChI
\item CSID <-> InChIKey
\item CSID <-> SMILES
\item CSID -> Mol file
\item InChI <-> InChIKey
\item InChI <-> SMILES
\item InChI -> Mol file
\item InChIKey <-> Mol file
}
}
\note{
An API key is needed. Register at \url{https://developer.rsc.org/}
for an API key. Please respect the Terms & Conditions. The Terms & Conditions
can be found at \url{https://developer.rsc.org/terms}.
}
\examples{
\dontrun{
cs_convert("BQJCRHHNABKAKU-KBQPJGBKSA-N",
  from = "inchikey", to = "csid"
)
cs_convert("BQJCRHHNABKAKU-KBQPJGBKSA-N",
  from = "inchikey", to = "inchi"
)
cs_convert("BQJCRHHNABKAKU-KBQPJGBKSA-N",
  from = "inchikey", to = "mol"
)
cs_convert(160, from = "csid", to = "smiles")
}
}
\references{
\url{https://developer.rsc.org/docs/compounds-v1-trial/1/overview}

Eduard Szöcs, Tamás Stirling, Eric R. Scott, Andreas Scharmüller,
Ralf B. Schäfer (2020). webchem: An R Package to Retrieve Chemical
Information from the Web. Journal of Statistical Software, 93(13).
\doi{10.18637/jss.v093.i13}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractors.R
\name{extractors}
\alias{extractors}
\alias{cas}
\alias{inchikey}
\alias{smiles}
\title{Extract parts from webchem objects}
\usage{
cas(x, ...)

inchikey(x, ...)

smiles(x, ...)
}
\arguments{
\item{x}{object}

\item{...}{currently not used.}
}
\value{
a vector.
}
\description{
Extract parts from webchem objects
}
\references{
Eduard Szöcs, Tamás Stirling, Eric R. Scott, Andreas Scharmüller,
Ralf B. Schäfer (2020). webchem: An R Package to Retrieve Chemical
Information from the Web. Journal of Statistical Software, 93(13).
\doi{10.18637/jss.v093.i13}.
}
