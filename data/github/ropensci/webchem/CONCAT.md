
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
```