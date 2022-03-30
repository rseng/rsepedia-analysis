# patRoon 2.0.1

- Fixed: `newProject()` would not add suspect annotation to the output script if the example suspect list or sets workflows were chosen.
- Fixed: Default optimization range for KPIC2 `min_width` was incorrect (PR #31, thanks to @@coltonlloyd)
- `installPatRoon()` improvements in determining what is already installed
- Fixed: Group qualities/scores were not transferred to new `featureGroups` objects after calling `screenSuspects()` or `unset()`
- Checking of MS file extensions (e.g. for `generateAnalysisInfo()`) is now case insensitive (see issue #34 and #43)
- Fixed: `newProject()`: `reportCSV()` call in generated script included non-existing `MSPeakLists` argument.
- Fixed: Suspect screening would result in an error if the `adduct` argument was specified (Corey Griffith)
- Refactoring of the reference documentation pages (issue #35)
    - Workflow data generation functions and their algorithms specific counterparts (e.g. `findFeatures`, `findFeaturesKPIC2`) are now documented on separate pages.
    - The plotting functions for `featureGroups` are now documented in a separate page (?`feature-plotting`)
    - Many small textual improvements were made in the process
- Fixed: the `findFeaturesKPIC2()` and `importFeaturesKPIC2()` now have correct casing (was lower case 'f')
- Fixed: some function arguments for `convertMSFiles()` were not properly verified
- More checks to verify if input mzML/mzXML data actually exists
- Improved reference documentation for `analysis-information` (issue #33)
- Handbook: detailed overview of all workflow functions and classes (issue #41, special thanks to @@hechth)
- Fixed: annotations slot for `featureGroups` was not updated when removing groups with `delete()` (except sets workflows)
- Fixed: better handle missing spectrum data with spectrumSimilarity()
- Made `nontarget` an optional dependency and install it from GitHub with CI and in the installation docs (see issue #48)
- Fixed: MS files were not always correctly found
- Fixed: `newProject()` ignored group/blank input for sets workflows
- Fixed: better error handling for suspect lists with only one valid column
- Improve suspect list handling when input `mz`/`rt` columns are not numeric
- `checkFeatures()`: don't show multiple rows if a suspect was matched with multiple feature groups. This change removed the option to show specific suspect columns.
- Fixed: `checkFeatures()` errored if Plot mode was 'Top most replicates' or 'All'
- Fixed: several issues with `topMost` plotting of EICs for sets data
- Fixed: `plotChroms()`: peak area filling (`showPeakAreas=TRUE`) didn't work if the peak height exceeded `ylim`
- `screenSuspects()` with sets workflows: don't warn about set specific suspect data if all data is NA
- `checkFeatures()`/`checkComponents()` now cleanup unavailable selections when saving the session
- Fixed: `reportHTML()`: Don't try to report TP components if no data is available
- `formulasSet` method for `plotSpectrum()`: don't try to plot a comparison plot for candidates without MS/MS data
- `reportHTML()` don't try to plot a comparison plot for formula candidates without MS/MS data


# patRoon 2.0

This release adds a significant amount of new functionality and changes. Please see the updated Handbook and sections below for more information.

Users of previous `patRoon` versions should inform themselves with the important changes highlighted in the next section. Furthermore, it is highly recommended to remove any cached data, i.e. by running `clearCache("all")` or manually removing the `cache.sqlite` file from your project directory.


## Important changes

- Features
    - XMCS(3): Renamed argument `exportedData` to `loadRawData`.
    - The `mzWindow` and `EICMzWindow` arguments were renamed to `mzExpWindow` / `EICMzExpWindow` and are now with slightly different meaning (please see the reference manual).
    - OpenMS: `minFWHM`/`maxFWHM` defaults lowered for `findFeatures` and feature optimization.
- Annotation
    - ggplot2 support for several plotting functions (i.e. `useGGPlot2` argument) is removed (not often used and a maintenance burden).
    - the `precursor` argument to the `plotSpectrum()`, `annotatedSpectrum()` and `plotScores()` methods for `formulas` now expects the neutral formula instead of the ionized formula. This change was done for consistency with compound annotations and sets workflows.
    - The way of obtaining candidate formulae from different features in the same group (i.e. _feature formula_ consensus) was changed.
        - Fixes were applied to improve thresholding with `featThreshold`.
        - A second and new threshold, `featThresholdAnn`, only takes annotated features into account.
        - The default of `featThreshold` is now `0`, for `featThresholdAnn` it is the same as the previous default for `featThreshold`.
        - Candidate results: renamed `analysis` column to `analysis_from` and added `analyses` column that lists all analyses from which the consensus was made.
        - if multiple annotations are available for a single MS/MS peak (eg due to differences between feature annotations) then only the annotation with lowest m/z deviation is kept (and a warning is emitted).
        - Scores of annotated fragments from different features are now averaged.
    - The storage classes and interface for formula and compound annotation was harmonized
        - The `formulas` and `compounds` classes now derive from the new `featureAnnotations` class. Most of the functionality common to formulas/compounds are defined for this class.
        - Storage of formula annotation results mostly follow the format that was already used for compound annotations.
        - The `maxFormulas`/`maxFragFormulas` argument for `as.data.table()` were removed, as these don't make much sense with the new format.
        - The `elements` filter now applies to neutral formulae for both formula and compound annotations (`fragElements` still applies to the ionized fragment formula).
    - Formula annotation with Bruker now require `MSPeakLists`. Since all algorithms now require peak lists, `generateFormulas` now has a mandatory `MSPeakLists` argument (similar to `generateCompounds`).
    - Formula candidates (formula and compound annotations) are now reported in the `ion_formula` (ionized) and `neutral_formula` (neutral) columns. Similarly, the `formula_mz` column was renamed to `ion_formula_mz`.
- Suspect screening
    - The methodology to match m/z values of suspects and features was changed. This was mainly for consistency and compatibility with sets workflows. Please see the updated section on suspect screening in the Handbook.


## Major new functionality

### Transformation product screening

The most important new functionality in `patRoon 2.0` are transformation product (TP) screening workflows. This release adds functionality to predict TPs (with `BioTransformer` or metabolic logic) or search TPs in `PubChem` or custom databases. Furthermore,  other data such as MS/MS similarity or feature classification data can be used to relate parent/TP features. Other TP screening functionality includes TP prioritization and automatic generation of TP compound database for `MetFrag` annotation. The workflows follow the classical design of `patRoon`, where flexible workflows can be executed with a combination of established algorithms and new functionality. For more information, please see the dedicated chapter about TP screening in the Handbook.

### Sets workflows

Another major change in this release is the addition of sets workflows. These workflows are typically used to simultaneously process positive and negative ionization data. Advantages of sets workflows include simplification of data processing, combining positive and negative data to improve e.g. feature annotations and easily comparing features across polarities. A sets workflow requires minimal changes to a 'classical workflow', and most of the additional work needed to process both polarities is done automatically behind the scenes. For more information, please see the dedicated chapter about sets workflows in the Handbook.

### Features

The following new feature detection/grouping algorithms were integrated: `SIRIUS`, `KPIC2` and `SAFD`. Furthermore, integration with `MetaClean` was implemented for the calculation of peak qualities and machine learning based classification of pass/fail peaks. In addition, the peak qualities are used to calculate peak scores, which can be used for quick assessment and prioritization.

### Data curation

Interactive curation of feature data with `checkChromatograms()` was replaced with `checkFeatures()`, which is much faster, is better suitable for larger datasets, customizable and has an improved user interface. Furthermore, this tool can be used for training/assessing `MetaClean` models. Similarly, `checkComponents()` is a function that allows interactive curation of component data.

The `delete()` generic function allows straightforward deletion of parts of the workflow data, such as features, components and annotations. Furthermore, this function makes it easy to implement customized filters.

### Adduct annotation

The algorithms of `OpenMS` (`MetaboliteAdductDecharger`) and `cliqueMS` were integrated for additional ways to detect adducts/isotopes through componentization. Furthermore, the new `selectIons()` method uses these annotations to prioritize features (e.g. by only retaining those with preferable adducts). In addition, this function stores the adduct annotations for the remaining feature groups, which can then be automatically used for e.g. formula and compound annotation.

## Other new functionality

- `newProject`
    - Updated for new functionality such as sets and TP workflows and adduct annotation.
    - Completely re-designed code generation to improve extensibility. The generated code will have a slightly different layout and some parameter defaults were changed.
- Features
    - `as.data.table()`
        - intensity normalization (`normFunc` argument)
        - customized averaging (`averageFunc` argument)
        - calculation of Fold-changes (`FCParams` argument)
        - report peak qualities/scores (`qualities` argument)
    - New `plotVolcano()` method function to plot fold changes.
    - `topMostByRGroup/EICTopMostByRGroup` arguments for plotting/reporting EIC data of only the top most intense replicate(s).
    - `reportHTML()` now only plots the EIC of the most intense feature of each replicate group (i.e. `EICTopMostByRGroup=TRUE` and `EICTopMost=1`).
    - `XCMS3`
        - `loadRawData` argument for feature grouping and `comparison()`
        - `...` argument for `findFeaturesXCMS3`
        - `preGroupParam` to specify grouping parameters used prior to RT alignment (suggested by Ricardo Cunha)
    - The internal `XCMS` feature (group) objects are synchronized as much as possible when feature data is changed.
    - Feature groups: print feature counts with `show()` and `filter()` methods.
    - `OpenMS` feature finding: `useFFMIntensities` argument to speed up intensity loading (experimental).
    - `reportHTML()` now reports general feature information in a separate tab.
    - Feature groups: `results` argument to `[` (subset) and `filter()` to quickly synchronize feature groups between objects (e.g. to quickly remove feature groups without annotation results).
- Annotation
    - The methodology of `plotSpectrum()` to automatically calculate the space necessary for formula annotation texts and candidate structures was improved. Annotation texts are now automatically resized if there is insufficient space, and the maximum size and resolution for candidate structures can be controlled with the `maxMolSize`/`molRes` parameters.
    - `filter()` method for `MSPeakLists`: `minMSMSPeaks` filter to only retain MSMS peak lists with a minimum number of peaks.
    - `filter()` method for `MSPeakLists`: `annotatedBy` filter to only keep peaks with formula/compound annotations.
- Suspect screening
    - The `screenSuspects()` method now supports the `amend` argument, which allows combining results of different `screenSuspects()` calls (see the Handbook for details).
- Components
    - A new algorithm, `specclust`, which generates components based on hierarchically clustering feature groups with high MS/MS similarities.


## Minor changes

- Features
    - `groupFeatures`: renamed the `feat` argument to `obj`.
    - Improved performance for some feature group filters.
    - `reportHTML()`: EICs are shared between tabs to avoid duplicated plotting
    - The `features` object embedded in `featureGroups` objects is now synchronized, and any features not present in any group are removed accordingly. This reduces memory usage and indirectly causes `reportCSV()` to only report features still present.
    - `plotInt()`: now has `xnames` and `showLegend` arguments to adjust plotting.
- Annotation
    - The `[` (subset) and `filter()` methods for `MSPeakLists` now only re-average peak lists if the new `reAverage` argument is set to `TRUE` (default `FALSE`). This change was mainly done as (1) the effects are usually minor and (2) re-averaging invalidates any formula/compound annotations done prior to filtering.
    - `filter()` method for `MSPeakLists`: the `withMSMS` filter is now applied after all other filters.
    - `MetFrag`: the raw unprocessed annotation formulas are now additionally stored in the `fragInfo` tables (`formula_MF` column).
    - `MetFrag`: the precursor ion m/z is now taken from peak list data instead of the feature group to improve annotation.
    - `MetFrag`: the `useSmiles` parameter is now set to `true` as this seems to improve results sometimes.
    - `as.data.table()` method for `formulas`: if `average=TRUE` then all column data that cannot be reasonably averaged are excluded.
    - `annotatedPeakList()`: also add annotation columns for missing results (for consistency).
    - Compound MS/MS annotations do not include peak intensities anymore (already stored in MS peak lists).
    - The `minMaxNormalization` argument to the `consensus()` method for `compounds` was removed (unused).
    - `filter()` for `formulas`/`compounds`: if algorithm consensus results are filtered with `scoreLimits`, and a score term exists multiple times for a candidate, only one of the terms needs to fall within the specified limits for the candidate to be kept (was all).
    - `plotSpectrum()` for `compounds`: `plotStruct` is now defaulted to `FALSE`.
    - `MSPeakLists` data now store an unique identifier for each mass peak in the `ID` column. These IDs are used by e.g. formula/compound annotations, and stored in the `PLID` column in their `fragInfo` data. This replaces the `PLIndex` column in `fragInfo` data, which was only row based, and therefore invalidated in case peak lists were filtered afterwards.
- Adducts
    - `GenFormAdducts()` and `MetFragAdducts()` now additionally return adducts in generic format and use cached data for efficiency. 
    - `err` argument to `as.character()` to control if an error or `NA` should returned if conversion fails.
    - `as.adduct()` now removes any whitespace and performs stricter format checks to make conversion more robust.
    - Standardized GenForm/MetFrag element addition/subtraction data to improve consistency for conversions (eg NH4 --> H4N).
    - Conversion from/to adduct formats of OpenMS (`MetaboliteAdductDecharger`) and `cliqueMS`.
    - `calculateIonFormula()` and `calculateNeutralFormula()` now Hill sort their result
    - The embedded GenForm code was updated to the latest version.
- Suspect screening
    - `as.data.table()`: Suspect screening specific columns are now prefixed with `susp_`.
    - The `suspFormRank` and `suspCompRank` suspect annotation data columns were renamed to `formRank`/`suspCompRank` (the previous change made prefixing unnecessary).
    - Several updates for Bruker TASQ import.
    - `logPath` argument for `annotateSuspects()` to specify the file path for log files are disable logging completely.
    - suspect names are now trimmed to 150 characters to avoid logging issues on e.g. Windows
- Components
    - Intensity clusters now use `fastcluster` for hierarchical clustering.
    - Changed column `rt` to `ret` for consistency.
    - `show()`: show unique feature group counts.
    - `filter()`: allow negative `rtIncrement` values.
    - `nontarget`: replaced `extraOpts` argument with `...`.
    - `nontarget`: store links as character string indices instead of numeric indices.
    - `RAMClustR`: moved position of `ionization` argument to improve consistency.
    - The 'reduced components' mechanism, where a components object was returned without any algorithm specific data (using the `componentsReduced` class) when filtering/subsetting components, was removed. This system was quite unintuitive and imposed unnecessary limitations. Instead, functions that cannot work after component data is changed (e.g. those specific to intensity clustering) will throw an error if needed.
    - The objects returned from `intclust` components are now derived from a general `componentsClust` class, which is shared with `specclust` components. The common functionality for both algorithms is implemented for this class.
- Misc
    - `show()` methods now print class inheritance tree
    - The `progressr` package is not used anymore, thus, it is not necessary to set up progress bars with `future` based multiprocessing.
    - `newProject`: Moved order of componentization step (now before annotation & suspect screening).
    - Plots of chromatograms, spectra etc that are without data now reflect this in the generated plot.


## Fixes

- Features
    - Blank filter: don't subtract blanks from each other
    - Fixed: when `xlim`/`ylim` was used with `plotChroms` then peaks were not always correctly filled
    - `retMin` argument to `plot()` method for `featureGroupsComparison` wasn't properly used/defaulted.
- Annotations
    - Fixed: `plotSpectrum()` if `xlim` is set and this yields no data then an empty plot is shown.
    - Fixed: `plotSpectrum()` automatic `ylim` determination was incorrect if only one peak is shown.
    - Fixed: consensus from feature formulas possibly could have fragment m/zs not in group MS/MS peak lists.
    - Fixed: consensus from feature formulas possibly could have fragment m/zs that deviated from those in in group MS/MS peaklists.
    - Fixed: formula algorithm consensus wrongly ranked candidates not ubiquitously present in all algorithms.
    - Fixed: the `scoreLimits` filter for formulas could ignore results not obtained with MS/MS data.
    - Fixed: MetFrag was using a wrong/inconsistent cache name.
    - Fixed: `as.data.table(compounds, fragments=TRUE)` returned empty results for candidates without fragment annotations.
    - Fixed: `topX` arguments for the `MSPeakLists` method for `filter()` would re-order peak lists, thereby invaliding any annotations.
    - Fixed: conversion of adducts with multiple additions/subtractions to GenForm/MetFrag format failed.
    - Fixed: Hill ordering: H wasn't sorted alphabetically if no C is present.
    - Several fixes were applied to improve handling of `SIRIUS` 'adduct fragments'.
    - formula/compound annotation consensus ranking is now properly scaled.
    - Fixed: `generateMSPeakListsDAFMF()` potentially used wrong DA compound data in case features were filtered.
- Suspect screening
    - Fixed: `numericIDLevel()` now properly handles `NA` values.
    - `importFeatureGroupsBrukerTASQ()`: Improved handling of absent analyses in imported results files.
    - Fixed: Automatic _m/z_ calculation for suspects
        - Improperly handled electron masses for adducts involving element subtract (e.g. `[M-H]-`), resulting in ~1.5 mDa deviations
        - Adduct conversion didn't handle multiple molecules (e.g. `[2M+H]+`) and multiple charges (e.g. `[M+2H]2+`)
- Components
    - `RAMClustR`: ensure that columns are the right type if all values are NA.
    - `CAMERA`: correctly handle cases when `minSize` filter results in zero components.
    - `plotGraph()`: improve error handling with empty objects.
- Misc
    - Future multiprocessing: make sure that logs are created even when an error occurs.
    - Classic multiprocessing: intermediate results are cached again.
    - Fixed: parallelization issues with cached data (thanks to https://blog.r-hub.io/2021/03/13/rsqlite-parallel/)
    - `newProject()`: correctly handle DIA with Bruker MS peak lists.


# patRoon 1.2.1

- Fixed: XCMS feature grouping didn't work when the `peakgroups` alignment method was used (fixes issue #22)
- Fixed: (harmless) `mapply` warning was shown with `newProject()`
- `newProject()`: don't show _Remove_ button in analyses select screen when the script option is selected, as this will not work properly.
- `IPO`: add default limits for OpenMS `traceTermOutliers`
- `IPO` optimization fix: integer parameters are properly rounded
- Fixed: `generateFeatureOptPSet("xcms3", method="matchedFilter")` would return a parameter set with `step` instead of `binSize` (issue #23)
- Fixed: `newProject()` would generate an ID levels configuration file even when no suspect list was selected.
- MCS calculation: handle `NULL` values that may occasionally be returned by `rcdk::get.mcs`
- Fixed: intensity filter failed if previous filters lead to zero feature groups.
- Fixed: `reportHTML()` annotation table was paged.
- Fixed: Check final path lengths of log files and truncate where necessary (reported by Corey Griffith)
- Fixed: in some cases the checking of `analysisInfo` validity may result in an error (reported by Tiago Sobreira)
- Fixed: `convertMSFiles()` error with `dirs=TRUE` (reported by Tiago Sobreira)
- Small updates/fixes for `installPatRoon()`
- Fixed: `screenSuspects()` did not take `onlyHits` into account for caching
- `screenSuspects()`: The original suspect name is stored in the `name_orig` column
- `XCMS3` feature group optimization: `binSize` and `minFraction` values were rounded while they shouldn't (issue #27)


# patRoon 1.2.0

This releases focuses on a significantly changed suspect screening interface, which brings several utilities to assist suspect annotation, prioritization, mixing of suspect and full NTA workflows and other general improvements.

**IMPORTANT**: The suspect screening interface has changed significantly. Please read the documentation (`?screenSuspects` and the handbook) for more details. If you want to quickly update your code without using any new functionality:

Change your existing code, e.g.

```r
scr <- screenSuspects(fGroups, suspectList, ...)
fGroupsScr <- groupFeaturesScreening(fGroups, scr)
```

to

```r
fGroupsScr <- screenSuspects(fGroups, suspectList, ..., onlyHits = TRUE)
```

**Major changes**

* New suspect screening interface
    * By default, feature groups without suspect hit are _not_ removed (unless `onlyHits=TRUE`). This allows straightforward mixing of suspect and full non-target workflows.
    * The feature groups are _not_ renamed tot the suspect name anymore. If you quickly want to assess which suspects were found, use the `screenInfo()` or `as.data.table()` methods.
    * Subsetting of suspsect screening results can be done with the `suspects` argument to `[`, e.g. `fGroupsScr[, suspects = "carbamazepine"]`
    * A new method, `annotateSuspects()`, allows combining the annotation workflow data (peak lists, formulas, compounds) to perform a detailed annotation for the suspects found during the workflow. This method calculates properties such as
        * Rankings: where is the suspect formula/compound ranked in the candidates from the workflow data
        * Annotation similarity: how well does the MS/MS spectrum of the suspect matches with formula/compound annotations.
        * An _estimation_ of identification levels to assist in quickly evaluating how well the suspect was annotated. The rules for identification levels are fully configurable.
    * A dedicated `filter()` method for suspect screening results, which allows you to easily prioritize data, for instance, by selecting minimum annotation ranks and similarities, identification levels and automatically choosing the best match in case multiple suspects are assigned to one feature (and vice versa).
    * A dedicated `as.data.table()` method and reporting functionality for suspect screening results to quickly inspect their annotation data.
    * Please refer to the updated suspect screening sections in the handbook and `?screenSuspects` and `?annotateSuspects` for  more information.
* Changes to suspect lists
    * Whenever possible, suspect information such as formulae, neutral masses, InChIKeys etc will be calculated for the input suspect list (obtainable afterwards with `screenInfo()`).
    * The suspect names will be checked to be file compatible, and automatically adjusted if necessary.
    * If MS/MS fragments are known for a suspect (formula or `m/z`), these can be included in the suspect list to improve suspect annotation.
    * The old suspect screening support for `features` objects was removed. The same and much more functionality can be obtained by the workflow for feature groups.
* The `reportCSV()` function was simplified and uses `as.data.table()` to generate the CSV data. This should give more consistent results.
* The `individualMoNAScore` MetFrag scoring is now enabled by default.

Other changes

* `reportHTML()` now allows toggling visibility for the columns shown in the feature annotation table.
* The `plotVenn()` method for `featureGroups` now allows to compare combinations of multiple replicate groups with each other. See `?plotVenn` for more information.
* Fix: locating `SIRIUS` binary on `macOS` did not work properly
* Fix: timeout warning for `GenForm` resulted in an error (https://github.com/rickhelmus/patRoon/issues/18)
* Fix: plotting structures resulted in a Java error on the RStudio Docker image (https://github.com/rickhelmus/patRoon/issues/18)

            
# patRoon 1.1

* **IMPORTANT**: The `plotEIC()`, `groups()` and `plotSpec()` methods were renamed to `plotChroms()`, `groupTable()` and `plotSpectrum()`. This was done to avoid name clashes with `XCMS` and `CAMERA`. The old functions still work (with a warning), but please update your scripts as these will be removed in the future.
* **IMPORTANT**: Major changes to the parallelization functionality
    * `patRoon` now supports an additional method to perform parallelization for tools such as `MetFrag`, `SIRIUS` etc. The main purpose of this method is to allow you to perform such calculations on external computer clusters. Please see the updated parallelization section in the handbook for more details.
    * The `logPath` and `maxProcAmount` arguments to functions such `generateFormulas`, `generateCompounds` etc were removed. These should now solely be configured through package options (see `?patRoon`).
    * The `patRoon.maxProcAmount` package option was renamed to `patRoon.MP.maxProcs`.
* Changes related to `SIRIUS`
    * **IMPORTANT:** Support for SIRIUS 4.5.0. Please update to this version since these changes break support for older versions.
    * Fix: SIRIUS formula calculation with `calculateFeatures=TRUE` would try to calculate formulas for features even if not present (eg after being removed by subsetting or filtering steps).
    * The `SIRBatchSize` argument to formula and compound generation functions was renamed to `splitBatches` and its meaning has slightly changed. Please see the reference manual (e.g. `?generateFormulas`) for more details.
* Changes related to MetFrag
    * Paths to local database files for MetFrag are now normalized, which makes handling of relative paths more reliable.
    * Changes in the specified local MetFrag database files are now inspected to improve caching.
    * Consistency: `generateCompoundsMetfrag` was renamed to `generateCompoundsMetFrag`.
* Optimized loading of spectra and EIC data.
* New utility functions
    * `withOpt()` to temporarily change (`patRoon`) package options.
    * `printPackageOpts()`: display current package options of `patRoon`.
* Finding features with `OpenMS`: potentially large temporary files are removed when possible to avoid clogging up disk space (especially relevant on some Linux systems where `/tmp` is small).
* Several packages such as `XCMS` are not attached by default, which significantly speeds up loading `patRoon` (e.g. with `library()`).
* The `compoundViewer()` function was marked as defunct, as it hasn;t been working for some time and its functionality is largely replaced by `reportHTML()`.
* `generateComponentsNontarget()`: update homolog statistics for merged series.
* `checkChromatograms()`: fix error when `fGroups` has only one replicate group
* `convertMSFiles()`: If `algorithm="pwiz"` and vendor centroiding is used then any extra filters are now correctly put after the `peakPicking` filter.
* `getXCMSnExp()` is now properly exported and documented.

# patRoon 1.0.4
* Small compatibility fixes for macOS
* Updated support for latest PubChemLite
* The `annoTypeCount` score for annotated compounds with PubChemLite is now not normalized by default anymore when reporting results.
* `reportHTML()` now correctly handles relative paths while opening the final report in a browser.


# patRoon 1.0.3
* `componentsNT`: include algorithm data returned by `nontarget::homol.search` in `homol` slot (suggested by Vittorio Albergamo)
* several `convertMSFiles()` fixes (issue #14)
    * prevent error when no input files are found
    * only allow one input/output format (didn't properly work before)
    * recognize that Waters files are directories
    * `cwt` option is now available for conversion with ProteoWizard
* minor fixes for subsetting XCMS `features` objects
* Fixed: `generateCompoundsMetFrag()`: compound names could be sometimes be interpreted as dates (reported by Corey Griffith)
* Fixed: on very rare cases empty peaklists could be present after averaging.
* Fixed: `SIRIUS` annotation didn't use set adduct but used default instead
* `SIRIUS` results are better handled if choosen adduct is not `[M+H]+` or `[M+H]+`
* More fixes for loading `data.table` objects properly from cache.
* RStudio Docker image: see the updated installation instructions in the handbook (thanks to Thanh Wang for some last minute fixes!)


# patRoon 1.0.2

* Fixed: avoid errors when SIRIUS returns zero results (reported by Vittorio Albergamo)
* Fixed: `plotGraph()` didn't properly handle components without linked series (reported by Vittorio Albergamo)
* Keep signal to noise data when importing/exporting XCMS data (`sn` column) (suggested by Ricardo Cunha)
* Reversed argument order of `exportedData`/`verbose` to `getXCMSSet()` functions to avoid ambiguities
* Automated tests for importing/exporting XCMS(3) data + small fixes for surfaced bugs
* `generateComponentsNontarget()`: allow wider _m/z_ deviation for proper linkage of different series (controlled by `absMzDevLink` argument).
* Fixed: `addAllDAEICs()` sometimes used wrong names for EICs
* Improved handling of empty feature groups objects when reporting
* Fixed: `reportPDF()` may report formula annotated spectra of results not present in input `featureGroups`
* Fixed: Loading `data.table` data from cache now calls `data.table::setalloccol()` to ensure proper behavior if `data.table::set()` is called on cached data.
* Fixed: plotSpec() for `compounds` with `useGGPlot2=TRUE` would try to plot formulas for non-annotated peaks (resulting in many extra diagonal lines)  
* Fixed: some functions involved in caching plot data for HTML reports sometimes returned invalid data.
* Fixed: EICs plotted by `reportPDF()` where not properly placed in a grid (as specified by `EICGrid` argument)
* Small tweaks/fixes for `reportHTML()`
    * now displays subscripted neutral formulae
    * Fixed: x axis title of EICs in annotation tab was cut-off
    * Fixed: The rt vs mz plot in the summary page now uses minutes for retention times if `retMin=TRUE`
* Updates for SIRIUS 4.4.29


# patRoon 1.0.1

* Perform neutral mass calculation for suspect screening with OpenBabel to avoid some possible Java/RCDK bugs on Linux.


# patRoon 1.0

## June 2020

* Fixed: `newProject()` didn't show polarity selection if only a compound identification algorithm was selected.
* Updated external dependency versions in installer script. 
* Fixed: `groupFeaturesXCMS3()` didn't properly cache results.
* `MSPeakLists`: results for averaged peak lists are now the same order as the input feature groups
* Fixed: XCMS(3) feature group import used wrong variable name internally (reported by Ricardo Cunha)

## May 2020

* **IMPORTANT** Major changes were made related to `SIRIUS` support 
    * Multiple features can now be annotated at once by `SIRIUS` (configurable with new `SIRBatchSize` function argument). This dramatically improves overal calculation times (thanks to Markus Fleischauer for pointing out this possibility!).
    * `generateFormulasSirius()` and `generateCompoundsSirius()` are now properly capitalized to `generateFormulasSIRIUS()` and `generateCompoundsSIRIUS()`
    * Support for `SIRIUS` 4.4.
    * If all features are annotated at once then `SIRIUS` output is directly shown on the console.
    * The amount of cores used by `SIRIUS` can be specified with the `cores` function arguments.
    * More extra commandline options can be given to `SIRIUS`
* Fixed: `groupNames()`, `analyses()` and similar methods sometimes returned `NULL` instead of an empty `character` vector for empty objects.
* `plotHeatMap()` with `interactive=TRUE`: switch from now removed `d3heatmap` package to `heatmaply`
* Fixed: `reportHTML()` didn't split PubChem URLs when multiple identifiers were reported.
* `PWizBatchSize` argument for `convertMSFiles()`


## April 2020

* `extraOptsRT`/`extraOptsGroup` arguments for OpenMS feature grouping to allow custom command line options.
* `importFeatureGroupsBrukerTASQ`
    * now correctly takes retention times of suspects into account when creating feature groups.
    * retention times / _m/z_ values are now averaged over grouped suspects.
* The `plot()` method for `featureGroups` now allows drawing legends when `colourBy="fGroups"` and sets `colourBy="none"` by default, both for consistency with `plotEIC()`.
* All documentation is now available as PDF files on the website (https://rickhelmus.github.io/patRoon/)
* `newProject()` now uses XCMS3 algorithms instead of the older XCMS interface.
* Fixed: features in objects generated by `xcms` (not `xcms3`) could not be subset with zero analyses (which resulted in errors by e.g. `unique()` and `reportHTML()`). Reported by Corey Griffith.


## March 2020

* Fixed: Normalization of scorings for formulae/compounds potentially uses wrong data after subsetting/filtering of `formulas`/`compounds` objects
* Suspect screening
    * Fixed: Errors/Warnings of missing data in suspect list were not shown if using cached data
    * If a value is missing in any of the columns from the suspect list used for automatic ion mass calculation (e.g. SMILES, formula, ...) then data from another suitable column is tried.
    * Fixed: invalid neutral mass calculation for suspects with charged SMILES/InChIs
    * Default adduct can be specified in `newProject()` dialog
* Small compatibility fix for feature finding with OpenMS 2.5 (reported by Thanh Wang)
* RAMClustR is now supported from CRAN and no need to install github package anymore
* pubchemlite identifiers are now URL linked in HTML reports
* related CIDs are now reported for PubChemLite results.
* MetFrag compound generation: removed `addTrivialNames` option as it never worked very well.
* `reportHTML()`: only components with reported feature groups are now reported.


## February 2020

* Several small improvements and fixes for TASQ import
* Suspect screening:
    * now also support chemical formula for automatic m/z calculation
    * more robust loading of suspect lists (e.g. skip suspects with missing/invalid data)


## January 2020

* Ignore user specified scorings for local databases such as CompTox that are not actually present in the DB. This makes it easier to use e.g. different DB versions with differing scorings.
* Add scorings from wastewater and smoking metadata comptox MetFrag databases
* Windows install script now install latest (March2019) CompTox
* Updates for latest PubChemLite relaease (Jan2020)
* Suspect screening now doesn't require pre-calculated ion `m/z` values. Instead, suspect lists can contain SMILES, InChI or neutral mass values which are used for automatic ion `m/z` calculation. See `?screenSuspects` for more details.


## December 2019

* Added missing score terms for latest CompTox MetFrag database
* labels parameter for formulas/compounds methods of `consensus()`
* Fixed: colour assignment for scores plotting of merged formulae/compound results might be incorrect (reported by Emma Schymanski)
* Fixed: analysis table in `newProject()` UI only showed partial amount of rows.
* Fixed: don't print normalized instead of absolute design parameters when only one parameter is optimized in DoEs (fixes issue #10)


## November 2019

* **IMPORTANT** The `addFormulaScoring()` function now uses a different algorithm to calculate formula scores for compound candidates. The score is now based on the actual formula ranking in the provided `formulas` object, and is fixed between _zero_ (no match) and _one_ (best match).
* Formula feature consensus:
    * All scorings are now averaged, including those that are not fragment specific (e.g. precursor m/z error)
    * This also improves ranking in certain specific cases
* Vectorized plotting of MS spectra to make it potentially faster
* Added PubChemLite support for MetFrag


## October 2019

* Fixed: `convertMSFiles` correctly checks if input exists
* Specific optimizations after benchmarking results:
    * `maxProcAmount` (i.e. number of parallel processes) now defaults to amount of physical cores instead of total number of CPU threads.
    * Decreased `batchSize` to `8` for GenForm formula calculation.
* `plot()` for `featureGroups` can now highlight unique/shared features across replicates (suggested by V Albergamo)
* Linking of homologous series:
    * Improved info descriptions for `plotGraph()`
    * Series are now properly unlinked when merging (was too greedy)
    * Better algorithm to detect conflicting series
    * Fixed bug when updating removed links
* `concs` option for `generateAnalysisInfo()` to set concentration data


## September 2019

* Labels for objects in a `featureGroupsComparison` can be customized (useful for e.g. plotting)
* Caching and progress bar support for suspect screening
* Updated/Fixed JDK installation for installation script
* Fixed missing pipe operator import (`%>%`)


## August 2019

* `topMost` argument for GenForm formula calculation.
* Added XCMS3 support for finding and grouping features, importing/exisiting data and parameter optimization (i.e. mostly on-par with classic XCMS support).
* Changed compound result column name from InChi to InChI


## June 2019

* **IMPORTANT** Several things are renamed for clarity/consistency
    * The column to specify replicate groups for blank subtraction in the analysis information is re-named from `ref` to `blank`. Similarly, the `refs` argument to `generateAnalysisInfo()` is now called `blanks`.
    * `reportMD()` is renamed to `reportHTML()`
    *  `filter()` method for `formulas`: `minExplainedFragPeaks` is now called `minExplainedPeaks`
    * `screenTargets` and its `targets` parameter have been renamed to `screenSuspects()` / `suspects`
* Fixed incorrect selection after feature table (or other interactive tables) have been manually re-ordered (reported by Thanh Wang)
* `groups()` and `as.data.table()` methods for `featureGroups`: optionally consider feature areas instead of peak intensities.
* `plotSilhouettes()` method for `compoundsCluster`
* Added `rGroups` argument to subset operator for `featureGroups` to subset by replicate groups (equivalent to `rGroups` argument to `filter()`).
* Improved logging of output from CLI tools (e.g. OpenMS, MetFrag, SIRUS, ...) 

## May 2019

* Formula updates
    * `GenForm` formula calculation with `MSMode="both"` (the default): instead of repeating calculations with and without MS/MS data and combining the data, it now simply does either of the two depending on MS/MS data availability. The old behavior turned out to be redundant, hence, calculation is now a bit faster.
    * `GenForm` now perform _precursor isolation_ to cleanup MS1 data prior to formula calculation. During this step any mass peaks that are unlikely part of the isotopic pattern of the feature are removed, which otherwise would penalize the isotopic scoring. The result is that isotopic scoring is dramatically improved now. This filter step is part of new filter functionality for `MSPeakLists`, see `?MSPeakLists` and `?generateFormulas` for more information.
    * When formula consensus are made from multiple features the scorings and mass errors are now averaged (instead of taking the values from the best ranked feature).
    * Improved ranking of candidates from a consensus of multiple formula objects (see `?formulas`).
* Consensus for compounds are now similarly ranked as formulas.
* More consistent minimum abundance arguments for `consensus()` (`absMinAbundance` and `relMinAbundance`)
* `MetFrag`: for-ident database and new statistical scores are now supported
* `as.data.table()` / `as.data.frame()` for `featureGroups` now optionally reports regression information, which may be useful for quantitative purposes. This replaces the (defunct) `regression()` method and limited support from `screenTargets()`.
* `plotGraph()` method to visually inspect linked homologous series.

## April 2019

* Misc small tweaks and fixes for `newProject()` (e.g. loading of example data).
* Improved graphical output of various common plotting functions.
* Updated tutorial vignette and added handbook


## March 2019
* `reportMD()`: most time consuming plots are now cached. Hence, re-reporting should be signficiantly faster now.
* Updates to MS data file conversion:
    * `convertMSFiles()` now (optionally) takes analysis information (`anaInfo`) for file input.
    * `convertMSFiles()` now supports Bruker DataAnalysis as conversion algorithm (replaces now deprecated `exportDAFiles()` function).
    * `MSFileFormats()` function to list supported input conversion formats.
    * `generateAnalysisInfo()` now recognizes more file formats. This is mainly useful so its output can be used with `convertMSFiles()`.
    * `convertMSFiles()` now has the `centroid` argument to more easily perform centroiding.
* Updates to `newProject()`:
    * The analyses selector recognizes more data file formats. This way you can select analyses that have not been converted yet.
    * Data pre-treatment options now include more sophisticated file conversion options (_e.g._ using ProteoWizard). This and the new analysis selector functionality ensures that data files in all major vendor formats do not have to be converted prior to generating a script.
    * Re-organized tabs to mirror non-target workflow.
    * Suspect screening support.
    * Improved layout of output script.
* `withMSMS` filter for MS peak lists.
* Timeout for formula calculation with GenForm to avoid excessive calculation times.
* `importFeatures()` generic function
* Reporting functions renamed arguments related to compounds reporting (e.g. compoundTopMost to compound**s**TopMost)


## February 2019
* Compound scorings are now normalized towards their original min/max values. This ensures that the `score` column of MetFrag results stays correct.
* plotScores(): Option to only report scorings that have been used for ranking
* as.data.table()/as.data.frame() method for compounds: optionally normalize scores.
* `reportPDF()`/`reportMD()` now report only 5 top most candidate compounds by default (controlled by `compoundsTopMost` argument).
* metadata for MS peak lists
* `plotSpec()` now displays subscripted formulae
* **IMPORTANT** Several major changes were made to the `filter()` methods for `features` and `featureGroups`. Please carefully read the updated documentation for these methods! (i.e. `` ?`filter,features-method` `` and `` ?`filter,featureGroups-method` ``).
    * Most argument have been renamed for consistency, simplicity and clarity.
    * The order when multiple filters are specified to the `featureGroups` method was adjusted, notably to improve reliability of blank filtration. Again, please see `` ?`filter,featureGroups-method` ``.
    * The following new filters were added:
        * mass defect range (`mzDefectRange` argument)
        * maximum relative standard deviation (RSD) of intensities between replicates (`maxReplicateIntRSD` argument)
        * minimum number of features within analyses (`absMinFeatures` and `relMinFeatures` arguments).
        * pre-intensity filters (`preAbsMinIntensity` and `preRelMinIntensity` arguments)
        * most existing filters now accept both relative and absolute values.
    * The script generation functionality of `newScript()` has been updated and supports more filter types.
    * The `repetitions` argument is not needed anymore for the new algorithm and has been removed.
    * `Inf` values now should be used to specify no maximum for range filters (was `-1`).
* Fixed: GenForm now always uses Hill sorting.
* `annotatedPeakList()` method for `formulas` and `compounds`. Also used by `reportMD` for improved annotation peak tables.
* Tweaked default mzR MS peak lists settings (halved `maxRtMSWidth` and `precursorMzWindow`)
* Fixed: Make sure that MetFrag web doesn't try to set unsupported database
* **IMPORTANT** Several changes were made to improve clarity and consensistency for arguments that specify retention/mz windows or allowable deviations.
    * Functions with changed argument names: `generateComponentsNontarget`, `generateComponentsRAMClustR`, `generateCompoundsSirius`, `generateFormulasGenForm`, `generateFormulasSirius`, `generateMSPeakListsDA`, `generateMSPeakListsMzR`, `importFeatureGroupsBrukerPA`
    * The `maxRtMSWidth` argument to `generateMSPeakListsDA`, `generateMSPeakListsMzR` (now `maxMSRtWindow`) now specifies a retention time window (\emph{i.e.} +/- retention time feature) instead of total retention width around a feature. Hence, _current input values should be halved_.
* CAMERA and RAMClustR components: both now have `minSize` and `relMinReplicates` (replaces `ubiquitous` for CAMERA) arguments. Note that their defaults may filter out (feature groups from) components. See their documentation for more info.
* Changed capitalisation of MetFrag CL location option from `patRoon.path.metFragCL` to `patRoon.path.MetFragCL`. The old name still works for backward compatability.
* Documented usage of the CompTox database with MetFrag. See `?generateCompounds`.
* Default normalization of MetFrag scorings now follows MetFrag web behaviour.
* `topMostFormulas` argument for SIRIUS compound generation.
* Fixed GenForm ranking in case both MS and MS/MS formulae are present.
* `reportPDF()`/`reportMD()` now report only 5 top most candidate formulae by default (controlled by `formulasTopMost` argument).
* Added `verifyDependencies()` function to let the user verify if external tools can be found.
* The meaning of the `dirs` argument to `convertMSFiles()` was slightly changed: if `TRUE` (the default) the input can either be paths to analyses files or to directories containing the analyses files.
* More effective locating ProteoWizard binaries by using the Windows registry.
* Nicer default graphics for `featureGroups` method for `plot()`.
* `reportMD()`: Don't plot Chord if <3 (non-empty) replicate groups are available. 
* All `filter()` methods now support negation by `negate` argument.


## January 2019
* minSize and ubiquitous arguments for CAMERA component generation. See ?generateComponentsCamera.
* Various tweaks for plotEIC() and plotSpec() methods
* Various small additions to newProject()
* `reportMD()`: added table with annotated fragments for compounds/formulas
* `consensus()` updates
    * `consensus()` methods now support extracting unique data. This also replaces the `unique()` method that was defined for `featureGroupsComparison`.
    * `comparison()` now automatically determines object names from algorithm (consistency with `consensus()` method for other objects).
    * Fixed: coverage calculation for consensus formulas now correctly based on precursor overlap (was overlap of precursor+fragment).    
* `plotVenn()` and `plotUpSet()` methods to compare different compounds or formulas objects.
* `filter()` method for components.
* DataAnalysis formula generation: fixed neutral formula calculation if `MSMode="msms"`, now needs `adduct` argument.
* Neutral loss filter for compounds.
* **IMPORTANT** Adduct specification is now simplified and more generic by usage of a new `adduct` class. This means that `generateCompounds()` and `generateFormulas()` now expect slightly differing arguments. Please see their manual pages.
* Workaround for homologous series generation with nontarget (see https://github.com/blosloos/nontarget/issues/6)
* Improvements to terminate background commandline processes when e.g. R is terminated. 
* `clearCache()` now supports removal of caches via regular expressions.
* Added/Improved `topMost` and `extraOpts` arguments for SIRIUS formula/compound generation.
* Annotated fragments from SIRIUS compounds now correctly contain charged molecular form.
* `filter()` method for compounds now support generic scoring filtering and on elements of precursor and fragment formulae.
* **IMPORTANT** Several changes were made to the MetFrag compound generation interface in order to simplify it and make it more generic. See `?generateCompounds` for more details (notably the Scorings section).
* More MS peak list updates
    * Precursor peaks are now flagged in MS peak list data and `plotSpec()`
    * Prune MS peak lists (not MS/MS) if no precursor could be determined (enabled by default, see `pruneMissingPrecursorMS` option in `?generateMSPeakLists`).
    * Better retain precursor peaks after filtering steps: only intensity thresholds may remove precursors (always for MS data, optional for MS/MS with `retainPrecursorMSMS` function arguments, see `?MSPeakLists` and `?generateMSPeakLists`).
* All major workflow classes now have `algorithm()` and `as.data.table()/as.data.frame()` methods. The latter replaces and enhances the `makeTable()` (`formulas` class) and `groupTable()` (`featureGroups` class) methods.


## December 2018
* Moved OpenMS XML writing code from `R` to `C++`: significantly reduces time required for grouping large amount of features.
* Several updates for functionality that uses Bruker DataAnalyses
    * Improved verification and consistency for handling processed data from DataAnalysis
    * Automatic saving & closing of analyses processed with DataAnalysis. Files are now generally closed by default to limit the resources needed by DataAnalysis. 
    * `revertDAAnalyses()` function: brings back set of Bruker analyses to their unprocessed state.
    * Minimum intensity arguments for Bruker DataAnalysis MS peak lists.
    * Slightly different `doFMF` behaviour for DataAnalysis feature finding.
* Several important updates were made to fomula calculation functionality.
    * The interface has been simplified as the functionality from the `formula` and `formulaConsensus` classes are now merged: there is no need to call `consensus()` anymore after `generateFormulas()`.
    * Formulae can now directly be calculated for feature groups by using group averaged MS peak lists (by setting `calculateFeatures=FALSE`). This can greatly speed up calulcation, especially with many analyses.
    * The new `filter()` and `as.data.table()`/`as.data.frame` methods bring new functionalities related to filtering, extracting data and performing several processing steps commonly performed for organic matter (OM) characterization.
    * Other updates on formulas
        * length now returns number of unique precursor formulas (was total number of results)
        * Fixed: Reported fragment formulas from SIRIUS were incorrectly assumed to be charged. Charged fragment formulas are now calculated manually (the neutral form is stored in the `frag_neutral_formula` column). This ensures correct comparison when a consensus is made.
        * `reportCSV()` now splits formulas for each feature group in separate CSV files (similar to `compounds` reporting).
        * Fixed: `reportPDF()` now actually includes formula annotations in annotated compound spectra when formulas are specified.
        * New oc argument when using GenForm: if enabled only organic formulae are accepted (i.e. with at least one carbon atom). Also incorporated a small fix for the actual GenForm binary to make this option work (https://sourceforge.net/p/genform/tickets/1/).
        * Fixed: coverage calculation of formulae across features treated formulae calculated only from MS data separately.
        * GenForm now also includes precursor annotation from MS/MS data.
* `file` argument for `clearCache()`
* Updates on MS peak lists
    * More consistent naming for algorithm specific MS peak list generators (i.e. `generateMSPeakListsX` where X is the algo).
    * Additional MS peak lists are generated by averaging the lists of features within a feature group.
    * `generateCompounds()` and plotting functionality now uses averaged group peak lists instead of peak list of most intense analysis.
    * `plotSpec()` method for MSPeakLists: plot (non-annotated) MS and MS/MS spectra.
    * Minimum intensity filter option that is applied after averaging.
    * Now uses "hclust" method for averaging by default, which now uses the [fastcluster] package.


## November 2018
* Default value for `maxRtMSWidth` argument used for peak list generation.
* Fixed: `maxRtMSWidth` argument for mzR peak list generation had no effect.
* Preliminary EPA DSSTox support (via LocalCSV).
* Added `addAllDAEICs()` function.
* Renamed `mzWidth` argument of `addDAEIC()` to `mzWindow`.
* Normalization of compound scores: normalization method can now be set and specified scorings can be excluded.
* Store/report IUPACName (as compoundName) from MetFrag PubChem data.
* Renamed trivialName to compoundName for compound tables.
* `convertMSFiles`: changed interface with more options, parallelization and ProteoWizard support.
* Automatic optimization of parameters necessary for feature finding and grouping. Heavily based on the IPO R package. See the 'feature-optimization' manual page.
* **IMPORTANT** `getXcmsSet()` is renamed to `getXCMSSet()`
* verbose option for `findFeatures()` / `groupFeatures()`
* Changed `nintersects` default for plotUpSet so that all intersections are plotted by default.
* plotChord() now properly stops if nothing overlaps.
* replicateGroupSubtract() now removes replicate groups that were subtracted.
* Fixed: replicateGroupSubtract() now correctly takes maximum mean intensity for threshold determination when multiple rGroups are specified.
* Fixed: Wrong compound clusters plotted in reportMD().
* Fixed: Added timeout after restarting failed command (e.g. MetFrag CL) to prevent rare error "The requested operation cannot be performed on a file with a user-mapped section open".


## October 2018
* OpenMS `features` class objects now store number of isotopes found for each feature.
* **IMPORTANT** Added all relevant options of FeatureFinderMetabo as function arguments to findFeaturesOpenMS() and renamed/reordered current options for more conistent style. Please check ?findFeatures for the updated function arguments!
* openReport option for reportMD(). If TRUE the generated report will be opened with a web browser.
* reportPlots option for reportMD() which collapses reportFGroups, reportChord and reportFormulaSpectra and adds control to plot Venn and UpSet diagrams.
* plotUpSet() methods to compare feature groups by UpSet plots. See e.g. http://caleydo.org/tools/upset/
* filter() method for features.
* EICs now loaded via faster C++ code thats uses mzR instead of XCMS
* Moved feature intensity loading code for OpenMS features to C++. This results in much faster feature finding.


## September 2018
* Removed filterBy methods: these are now deprecated with new subset operators and groupNames()/analyses() methods. Example: `fGroups <- fGroups[, groupNames(compounds)]`
* subset/extraction operators ("[", "[[" and "$") for features, featureGroups, MSPeakLists, formulas, formulaConsensus, compounds, compoundsCluster and components classes.
* analyses() and groupNames() generics to get analyses and feature group names of the data within an object.
* "[" method for featureGroups: empty feature groups now always dropped, drop argument now ignored.
* reportMD(): The layout to show compounds, formulas and components is now done with DataTables (DT package). This change allows faster initial loading of results. Furthermore, several small tweaks were done to improve general design.
* plotSpec() (compounds method): remove unused normalizeScores flag
* plotSpec() (compounds method): plotting of embedded structure now optional (plotStruct argument)
* plotSpec() (compounds method): automatic calculation of necessary extra height to plot labels/structure


## Augustus 2018
* The XML code required to load feature (group) data generated by OpenMS is now moved to a C++ interface that uses [Rcpp] and [pugixml]. This results in a significant reduction of required processing time. In addition, files are now processed in chunks, allowing even very large feature sets (>10000) without clogging up system memory.
* Improved general numeric comparisons, resulting in e.g. improved EIC generation.
* Tweaked OpenMS feature intensity loading: now takes intensity from data point closest to retention time instead of max intensity from datapoints in the search window. Furthermore, the search window for datapoints was reduced and made configurable.


## July 2018
* getMCS() method for compounds
* plotStructure() method for compounds will draw MCS when mutiple indices are specified


## June 2018
* Added removeRefAnalyses argument to filter() (featureGroups method) to easily remove e.g. analyses that are used as blanks after blank subtraction.
* Added filterBy() method which removes any feature groups from a featureGroups object of which no results are present in a specified object. Methods are defined for MSPeakLists, formulaConsenus, compounds and components. This method replaces some of the functionality of the filter() method for featureGroups (formConsensus and compounds arguments).
* Added mz and chromatographic peak width range options to filter() method for feature groups.
* Moved intensity clustering code (makeHCluster) to new component type (see componentsIntClust class documentation).


## May 2018
* Added compound clustering (see makeHCluster method for compounds). This is an useful tool to get an overview of all the candidate chemical structures after compound identification. The clustering will reduce data complexity. Furthermore, maximum common sucstructures (MCS) can be calculated and plotted for each cluster to get a quick impression of the different structures of candidates.
* Added function arguments checks using [checkmate]. This guards all exported functions and methods from wrong user input arguments. If any problems are found (e.g. a wrong data type or range was specified) then the user is informed about this and what kind of input is to be expected.
* Added workaround (removed latex dependencies added automatically by `kableExtra` package) that may cause memory leakage when `reportMD()` is called repeatedly.


## April 2018
* Added unit tests (using [testthat]) and fixed several small bugs that were revealed in the process.
* Continuous integration (CI) with running tests on [CircleCI] (Linux builds) and [AppVeyor] (Windows builds) and testing coverage on [Codecov]. Docker images with patRoon and all its dependencies are automatically pushed on [Docker Hub][DH].
* Many small bug fixes.




[Rcpp]: http://www.rcpp.org/
[pugixml]: https://pugixml.org/
[checkmate]: https://github.com/mllg/checkmate
[testthat]: https://github.com/r-lib/testthat
[CircleCI]: https://circleci.com/gh/rickhelmus/patRoon
[AppVeyor]: https://ci.appveyor.com/project/rickhelmus/patroon/branch/master
[Codecov]: https://codecov.io/gh/rickhelmus/patRoon
[DH]: https://hub.docker.com/r/patroonorg/patroon/
[fastcluster]: https://cran.r-project.org/web/packages/fastcluster/index.html
# Priority

# Lower priority


## General

- add 'keep.rownames = FALSE' to all as.data.table methods (or see if there is a work-around)
- remove mz column from patRoonData suspects?
- convertMSFiles()
    - Agilent .d is also a directory?
    - Remove necessity to have different input/output formats? (at least OK for pwiz)
- allow specifying average function in other places where as.data.table() is used (eg clustering, plotting etc)
- delete() for other classes
- generalize makeLegend() to new plot util


## Features

- misc
    - OpenMS: alignment may yield negative RTs...
    - OpenMS MapAligner exception
        - seems to be related when little overlap between sets --> add note in doc?
- adduct annotations
    - selectIons(): prefer adducts based on MS/MS? eg handy for Na/K adducts
    - what to do with unsupported adducts for annotation?
	    - default selectIons() to only consider 'common' adducts? or change default adducts for componentization algos?
	    - check better for what is supported by SIRIUS?
- import XCMS features: verify anaInfo (or remove necessity, eg with importAnaInfo func)
- getEICsForFeatures method for kpic2?
- optimize hashing? Or further avoid hashing big objects/objects with lists?
- load OpenMS intensities in parallel
    - either with futures or with MP and cache intensities afterwards
- XCMS: multiple features grouped in same analysis?
    - can be, but now handled by default method="medret" param. Make this configurable?
- updatePICSet(): also sync peaks list? otherwise doc
- Somehow integrate XCMS::fillChromPeaks
- Export generic EIC generation, i.e. without the need of feature data


## Annotation

- SusDat MF support
- parallel MSPeakLists generation?
- somehow handle different fragment formula annotations when making a consensus between formula/compounds objects
- DA formulas: also rank formula results like GF/SIRIUS?
- plotSpectrum/spectrumSimilarity: allow separate MSLevel for comparisons
- Support multiple MS/MS formula annotation candidates (ie same MS/MS peak annotated with different formulas)
    - mainly relevant for GenForm
    

## Components

- feature components
    - cliqueMS
        - change checkPackage GH link once PRs are merged
        - maxCharge --> chargeMax (same as OpenMS)? update docs
        - apply sort fix? https://github.com/osenan/cliqueMS/issues/8
- RC: check spearmans correlation
- NT: minimum size argument, combine rows for multiple rGroups?
- int and others: also use calculateComponentIntensities() for intensities?
- plot doesnt work for componentsReduced that originates from cluster components
    - maybe drop reduced mechanism?
- intclust
    - optionally take areas instead of intensities
    - cache results
- import check sessions?
    - needs way to match component names
    

## TPs

- predictTPsBioTransformer()
    - do we still need to check for non-calculated formulae?
- spectrumSimilarity
    - defaults OK for sim params?
        - precursor FALSE?
        - thresholds not really handy for formulas/compounds
            - at least doc that annotation results may disappear


## tests

- checkComponents() / checkFeatures()
    - server tests?
    - import
- MSPeakLists and others?: also test object that is fully empty (now still has analyses)
- ensure peaklists are sorted
- features
    - new multiple blank filtering
    - syncing of XCMS/KPIC2 objects
    - check if featindex and groups slots are in sync with features
    - subsetting and groupScores
    - as.data.table: normalization?
    - plotInt sets arg?
- components: somehow verify adductConflictsUsePref


## docs

- ref docs
    - delete()
        - mention j=DT for fGroups method?
    - order of sets sections
- handbook
    - TPs
        - add MSPL filtering of annotated peaks (in examples)?
    - update
        - introduction
            - mention changes of patRoon 2.0?
- tutorial
    - check if all is OK
- TPs
    - doc that merging TPs (same fGroup/TP) could be done with suspect screening
    - logic: mention results can be filtered with TP components?
    - mention how sets filter work for componentsTPs?


# Future


## General

- test negative subset indices
- convertMSFiles(): Support OpenMS vendor conversion? (eg thermo)
- newProject()
    - also allow suspect annotation with only peak lists? currently only selectable if formulas/compounds selected
- Reduce non-exported class only methods
- future MP
    - delayBetweenProc?
    - batch mode
- msPurity integration
- algorithmObject() generic: for xset, xsa, rc, ...
- more withr wrapping? (par)
- newProject()
    - concentration column for anaInfo
    - generate more detailed script with e.g. commented examples of subsetting, extraction etc
	- import Bruker seq file?
    - fix multi line delete (when possible)


## Features

- makeSet(): also support fGroups method via comparison?
- feature optim:
    - keep retcor_done?
    - get rid of getXCMSSet() calls?
- filter()
    - document which filters work on feature level (e.g. chromWidth)
    - remove zero values for maxReplicateIntRSD?
- integrate OpenMS feature scoring and isotopes and PPS in general (also include filters?)
- OpenMS: Support KD grouper?
- Integration of mzMine features (package pending...), MS-DIAL and peakonly?
- suspect screening
    - automatic suspect list name assignment if that's lacking? might be handy for some NORMAN lists
- topMost filter that accepts rGroups, either as AND or OR

## Annotation

- MSPeakLists
    - isotope tagging is lost after averaging
    - test avg params
    - metadata() generic?
    - DA
        - generateMSPeakListsDA: find precursor masses with larger window
        - tests
            - utils? EICs with export/vdiffr?
            - test MS peak lists deisotoping?
- metadata for Bruker peaklists?
- SIRIUS: use --auto-charge instead of manually fixing charge of fragments (or not? conflicting docs on what it does)
- test score normalization?
- timeouts for SIRIUS?
- do something about negative H explained fragments by MF?
- MetFrag: auto-include suspect results if suspectListScore is selected?
- do something with sirius fingerprints? --> comparison?
- fix compoundViewer
- add new MF HD scorings and make sure default normalization equals that of MF web
- CFM-ID and MS-FINDER integration
- utility functions to make custom DBs for MetFrag and SIRIUS and support to use them with the latter
- DBE calculation for SIRIUS?
- OM reporting
- as.data.table: option to average per replicate group?
- ID levels for non-suspects
    - function to calculate ID levels from suspect list (to take RTs/MSMS if available), formulas, compounds
    - store in compounds?
    - does it make sense for formula candidates?
    - add into reporting
        - also mark if in suspect list

## Suspects

- ID level rules: add scorings for SIRIUS/DA
- interface
    - also convert TASQ?
    - annotateSuspects()
        - check why it's is sometimes slow
            - seems to be logging, disable by default? --> only slow with testthat?
    - don't assign level <1 if suspect is a target? or give the choice (or make filter?)
- misc
    - prepareSuspectList(): export?
        - mainly to allow merging of lists, perhaps make util for that instead? Would also be handy for MF databases
            - could also fix column names, replace "-" with NAs etc
        - if yes, mention in ref docs for screenSuspects()
- expand reporting
    - eg include suspect name in EICs
        - already now in featInfo
    - mention suspect similarities/ranks etc for candidates (or somehow in compounds?)
    - optionally report with collapsed suspects



## components
- mass defect components
- split peak correlation and adduct etc annotation? would allow better non-target integration
- fillPeaks for CAMERA (and RAMClustR?)
- feature components
    - cliqueMS
        - current adduct conversion to this format doesn't mimic Cat and 2H/2Na etc
            - Perhaps just document limitation?
    - minimal annotation abundance across analyses (eg adduct must be annotated in >=X analyses)?
    - prefAdducts: also include eg Na by default?


## Sets
- compound/formula set consensus
    - weights for ranking (like compound consensus)?


## TPs

- filter on stability/persistence/toxicity of TP?


## Reporting

- add more options to reportPlots argument of reportHTML()?
- onlyAnnotated argument

# patRoon

[![CircleCI](https://circleci.com/gh/rickhelmus/patRoon.svg?style=svg)](https://circleci.com/gh/rickhelmus/patRoon)
[![Build status](https://ci.appveyor.com/api/projects/status/52nnpq8kqpkjqc92/branch/master?svg=true)](https://ci.appveyor.com/project/rickhelmus/patroon/branch/master)
[![codecov](https://codecov.io/gh/rickhelmus/patRoon/branch/master/graph/badge.svg)](https://codecov.io/gh/rickhelmus/patRoon)
[![Docker image](https://img.shields.io/docker/image-size/patroonorg/patroonrs/latest)][DockerImg]

`patRoon` aims to provide comprehensive mass spectrometry based non-target analysis (NTA) workflows for environmental
analysis. The name is derived from a Dutch word that means _pattern_ and may also be an acronym for _hyPhenated mAss
specTROmetry nOn-target aNalysis_.

## Project news

**December 2021** `patRoon 2.0` is now available. This major new release adds functionality to automatically screen and
identify transformation products, process positive and negative ionization MS data simultaneously and combine the
results, new algorithms for feature and adduct detection, interactive data curation and more. Please see the [Project
NEWS][NEWS] for details.

## Introduction

Mass spectrometry based non-target analysis is used to screen large numbers of chemicals simultaneously. For this
purpose, high resolution mass spectrometry instruments are used which are typically coupled (or _hyphenated_) with
chromatography (_e.g._ LC or GC). The size and complexity of resulting data makes manual processing impractical. Many
software tools were/are developed to facilitate a more automated approach. However, these tools are generally not
optimized for environmental workflows and/or only implement parts of the functionality required.

`patRoon` combines established software tools with novel functionality in order to provide comprehensive NTA workflows.
The different algorithms are provided through a consistent interface, which removes the need to know all the details of
each individual software tool and performing tedious data conversions during the workflow. The table below outlines the
major functionality of `patRoon`.

Functionality | Description | Algorithms
---------------------- | ------------------------------------------------------------------------ | -----------------------
Raw data pre-treatment | MS format conversion (e.g. vendor to `mzML`) and calibration.            | [ProteoWizard], [OpenMS], [DataAnalysis]
Feature extraction     | Finding features and grouping them across analyses.                      | [XCMS], [OpenMS], [enviPick], [DataAnalysis], [KPIC2], [SIRIUS], [SAFD]
Suspect screening      | Finding features with suspected presence by MS and chromatographic data. Estimation of identification confidence levels. | Native
MS data extraction     | Automatic extraction and averaging of feature MS(/MS) peak lists.        | Native, [mzR], [DataAnalysis]
Formula annotation     | Automatic calculation of formula candidates for features.                | [GenForm], [SIRIUS], [DataAnalysis]
Compound annotation    | Automatic (_in silico_) compound annotation of features.                 | [MetFrag], [SIRIUS], Native
Componentization & adduct annotation | Grouping of related features based on chemistry (e.g. isotopes, adducts and homologs), hierarchical clustering or MS/MS similarity into components. Using adduct and isotope annotations for prioritizing features and improving formula/compound annotations. | [RAMClustR], [CAMERA], [nontarget R package][nontarget], [OpenMS], [cliqueMS], Native
Combining algorithms   | Combine data from different algorithms (e.g. features, annotations) and generate a consensus. | Native
_Sets workflows_       | Simultaneous processing and combining +/- MS ionization data             | Native
Transformation product (TP) screening | Automatic screening of TPs using library/_in-silico_ data, MS similarities and classifications. Tools to improve compound TP annotation. | [BioTransformer], [PubChemLite][PubChemLiteTR], Native
Reporting              | Automatic reporting in _CSV_, _PDF_ and (interactive) _HTML_ formats. An example HTML report can be viewed [here][example]. | Native
Data clean-up & prioritization | Filters for blanks, replicates, intensity thresholds, neutral losses, annotation scores, identification levels and many more. | Native
Data curation          | Several graphical interactive tools and functions to inspect and remove unwanted data. | Native

The workflow of non-target analysis typically depends on the aims and requirements of the study and the instrumentation
and methodology used for sample analysis. For this reason, `patRoon` does not enforce a certain workflow. Instead, most
workflow steps are optional, fully configurable and algorithms can easily be mixed or even combined.

## Implementation details

* `patRoon` is implemented as an [R] package, which allows easy interfacing with the many other `R` based MS tools and other data processing functionality from `R`.
* Fully open-source (GPLv3).
* Developed on Windows, Linux and macOS
* S4 classes and generics are used to implement a consistent interface to all supported algorithms.
* Continuous integration is used for automated unit testing, automatically updating the [Website][ghweb] and documentation and maintaining a [miniCRAN] [repository][patRoonDeps] and [Docker image][DockerImg] to simplify installation (see [the handbook][handbook-inst] for more details).
* Supports all major instrument vendor input formats (through usage of [ProteoWizard] and [DataAnalysis]).
* Optimizations
    * `data.table` is used internally as a generally much more efficient alternative to `data.frame`.
    * The [processx] and [future] `R` packages are used for parallelization.
    * Results from workflow steps are cached within a [SQLite] database to avoid repeated computations.
    * Code for loading MS and EIC data, MS similarity calculations and others were implemented in `C++` to reduce computational times.
* The [RDCOMClient] `R` package is used to interface with Bruker DataAnalysis algorithms.
* The [Shiny] `R` package was used to implement several GUI tools.


## Installation

`patRoon` itself can be installed as any other `R` package, however, some additional installation steps are needed to
install its dependencies. Alternatively, [R Studio][RStudio] based Docker images are available to easily deploy a
complete `patRoon` environment. Please see the [installation section in the handbook][handbook-inst] for more
information.


## Getting started

For a very quick start:

``` r
library(patRoon)
newProject()
```

The `newProject()` function will pop-up a dialog screen (requires [R Studio][RStudio]), which will allow you to quickly
select the analyses and common workflow options to subsequently generate a template `R` processing script.

However, for a better guide to get started it is recommended to read the [tutorial]. Afterwards the [handbook] is a
recommended read if you want to know more about advanced usage of `patRoon`. Finally, the [reference] outlines all the
details of the `patRoon` package.


## Citing

When you use `patRoon` please cite its publications:

Rick Helmus, Thomas L. ter Laak, Annemarie P. van Wezel, Pim de Voogt and Emma L. Schymanski. [patRoon: open source
software platform for environmental mass spectrometry based non-target
screening](https://doi.org/10.1186/s13321-020-00477-w). _Journal of Cheminformatics_ **13**, 1 (2021)

Rick Helmus, Bas van de Velde, Andrea M. Brunner, Thomas L. ter Laak, Annemarie P. van Wezel and Emma L. Schymanski.
**patRoon 2.0: Improved non-target analysis workflows including automated transformation product screening**. _In
preparation_

`patRoon` builds on many open-source software tools and open data sources. Therefore, it is important to also cite their
work when using these algorithms via `patRoon`.

## Contributing

For bug reports, code contributions (pull requests), questions, suggestions and general feedback please use the [GitHub page](https://github.com/rickhelmus/patRoon).


[R]: https://www.r-project.org/
[NEWS]: https://github.com/rickhelmus/patRoon/blob/master/NEWS.md
[XCMS]: https://github.com/sneumann/xcms
[OpenMS]: http://openms.de/
[enviPick]: https://cran.r-project.org/web/packages/enviPick/index.html
[KPIC2]: https://github.com/hcji/KPIC2
[SAFD]: https://bitbucket.org/SSamanipour/safd.jl/src/master/
[DataAnalysis]: https://www.bruker.com/
[ProfileAnalysis]: https://www.bruker.com/
[mzR]: https://github.com/sneumann/mzR/
[GenForm]: https://sourceforge.net/projects/genform
[SIRIUS]: https://bio.informatik.uni-jena.de/software/sirius/
[MetFrag]: http://c-ruttkies.github.io/MetFrag/
[RAMClustR]: https://github.com/sneumann/RAMClustR
[CAMERA]: http://msbi.ipb-halle.de/msbi/CAMERA/
[nontarget]: https://cran.r-project.org/web/packages/nontarget/index.html
[cliqueMS]: https://github.com/osenan/cliqueMS
[BioTransformer]: https://bitbucket.org/djoumbou/biotransformer/src/master/
[PubChemLiteTR]: https://doi.org/10.5281/zenodo.5644560
[future]: https://github.com/HenrikBengtsson/future
[pngquant]: https://pngquant.org/
[tutorial]: https://rickhelmus.github.io/patRoon/articles/tutorial.html
[handbook]: https://rickhelmus.github.io/patRoon/handbook_bd/index.html
[handbook-inst]: https://rickhelmus.github.io/patRoon/handbook_bd/installation.html
[reference]: https://rickhelmus.github.io/patRoon/reference/index.html
[RStudio]: https://www.rstudio.com/
[processx]: https://github.com/r-lib/processx
[SQLite]: https://www.sqlite.org/index.html
[RDCOMClient]: http://www.omegahat.net/RDCOMClient/
[Shiny]: https://shiny.rstudio.com/
[example]: https://rickhelmus.github.io/patRoon/examples/report.html
[ProteoWizard]: http://proteowizard.sourceforge.net/
[ghweb]: https://rickhelmus.github.io/patRoon/
[patRoonDeps]: https://github.com/rickhelmus/patRoonDeps
[miniCRAN]: https://cran.r-project.org/web/packages/miniCRAN/index.html
[DockerImg]: https://hub.docker.com/r/patroonorg/patroonrs
---
title: "patRoon tutorial"
author: "Rick Helmus"
date: "`r Sys.Date()`"
header-includes:
    - \usepackage{fvextra}
    - \DefineVerbatimEnvironment{Highlighting}{Verbatim}{breaklines,commandchars=\\\{\}}
vignette: >
    %\VignetteIndexEntry{patRoon tutorial}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
    error = FALSE
)

source(file.path("shared", "init.R"))
doCompounds <- TRUE # curl::has_internet()
```

```{css code=readLines("styles.css"),echo=FALSE}
```


# Introduction

In this tutorial you will learn how to perform a simple non-target analysis with `patRoon`. This tutorial is not meant to give a detailed overview of `patRoon`. Instead, it serves as a quick introduction on how to use `patRoon` to setup and perform a full non-target analysis workflow.

The workflow in this tutorial consists of the following steps:

```{r workflow,echo=FALSE,out.width="75%"}
plotGV("
digraph rmarkdown {
graph [ rankdir = LR ]
node [ shape = box,
       fixedsize = true,
       width = 2.2,
       height = 1,
       fontsize = 18,
       fillcolor = darkseagreen1,
       style = filled ]
'New project' -> Features -> Annotation -> Reporting [ arrowhead = vee ]
}
", height = 100, width = 750)
```


# Data

In this tutorial we will use example data provided within the [patRoonData] package. Please make sure this package is installed (see the [Handbook] for brief installation instructions). The example dataset contains LC-MS data for a standard mixture with known composition ('standard-X') and a blank solvent ('solvent-X'), both in triplicate and measured with positive and negative ionization. While this may not seem like the most exciting data, it does allow to demonstrate the most important functionality of `patRoon`.

The provided analyses already have been exported to an open format (_.mzML_) and are ready to use. For your own data it may be necessary to first export your data to _mzXML_ or _mzML_ and perform other data pre-treatment steps such as mass re-calibration. This can be done using the tools from [ProteoWizard] or software from your mass spectrometer vendor. Alternatively, `patRoon` can do this automatically for your analyses with the `convertMSFiles()` function. Please see the handbook and reference manual for its usage.

# New project

Whenever you start a new non-target analysis it is highly recommended to start this from a fresh _project directory_. This directory will contain your `R` processing script(s) and any other output generated during the workflow. Note that this directory does not have to contain the raw MS data files. In fact, keeping these files separate may be handy, for instance, if you want to run multiple non-target analyses on these files or store the analysis files on a shared location.

Starting a new project typically consists of

1. Creating a new directory (unsurprisingly!)
2. Changing the active working directory to the project directory (e.g. with `setwd()`).
3. Create (or copy) an `R` processing script.

Note that step 2 is important as any output files (e.g. reports and cached results) are stored to the current working directory by default. Consequently, always take care to ensure that this directory is active, for instance, after restarting `R`.

Steps 1-3 can be easily performed with the `newProject()` function. Alternatively, you can of course also perform these steps yourself. Both approaches will be discussed in the next sections.

## Automatic project creation

Ensure that RStudio is active and start the new project utility:

```{r eval=FALSE}
patRoon::newProject()
```

> **_NOTE_** Currently `newProject()` _only_ works when using RStudio.

A dialog should pop-up (see screenshot below) where you can specify where and how the new project will be generated, which analyses you want to include and define a basic workflow. Based on this input a new project with a template script will be automatically generated.

![](newp.png){width=450px}

For this tutorial make the following selections

* **Destination tab** Select your desired location of the new project. Leave other settings as they are.
* **Analyses tab** Here you normally select your analyses. However, for this tutorial simply select the _Example data_ option.
* **Data pre-treatment tab** Since the example data is already ready to use you can simply skip this tab.
* **Features tab** Leave the default OpenMS algorithm for feature finding and grouping.
* **Annotation tab** Select _GenForm_, _MetFrag_ and _mzR_ for the formula generation, compound identification and peak list generator options, respectively (note that the last will become visible when selecting either of the other options).
* **TP screening tab** No need to do anything here for this tutorial.
* **Reporting tab** Make sure to enable HTML reporting.

## Manual project creation

For RStudio users it is easiest to simply create a new RStudio project (e.g. _File_ --> _New Project_). This will create a new directory and ensure that the working directory is set whenever you re-open it. Alternatively, you can do this manually, for instance:

```{r eval=FALSE}
projDir <- "~/myProjectDir"
dir.create(projDir)
setwd(projDir)
```

The next step is to create a new `R` script. For this tutorial simply copy the script that is shown in the next section to a new `.R` file.

## Template R script

After you ran `newProject()` the file below will be created. Before running this script, however, we still have to add and modify some of its code. In the next sections you will learn more about each part of the script, make the necessary changes and run its code.

```{r scriptPre,code=readLines("script-pre.R"),eval=FALSE}
```


# Workflow

Now that you have generated a new project with a template script it is time to make some minor modifications and run it afterwards. In the next sections each major part of the script (initialization, finding and grouping features, annotation and reporting) will be discussed separately. Each section will briefly discuss the code, what needs to be modified and finally you will run the code. In addition, several functions will be demonstrated that you can use to inspect generated data. 

## Initialization

The first part of the script loads `patRoon`, makes sure the current working directory is set correctly and loads information about the analyses. This part in your script looks more or less like this:

```{r init,eval=FALSE}
library(patRoon)

workPath <- "C:/my_project"
setwd(workPath)


# Example data from patRoonData package (triplicate solvent blank + triplicate standard)
anaInfo <- patRoonData::exampleAnalysisInfo("positive")
```

```{r doInit,include=FALSE}
# negative numeric eval option to disable lines is useless nowadays as it
# visually comments out the code --> guess we need two blocks
# https://github.com/yihui/knitr/issues/558

library(patRoon)
anaInfo <- patRoonData::exampleAnalysisInfo("positive")
```

After you ran this part the analysis information should be stored in the `anaInfo` variable. This information is important as it will be required for subsequent steps in the workflow. Lets peek at its contents:

```{r anaInfo}
anaInfo
```

As you can see the generated `data.frame` consists of four columns:

* *path*: the directory path of the file containing the analysis data
* *analysis*: the name of the analysis. This should be the file name _without_ file extension.
* *group*: to which _replicate group_ the analysis belongs. All analysis which are replicates of each other get the same name.
* *blank*: which replicate group should be used for blank subtraction.

The latter two columns are especially important for [data cleanup](#data-cleanup), which will be discussed later. For now keep in mind that the analyses for the solvents and standards each belong to a different replicate group (`"solvent"` and `"standard"`) and that the solvents should be used for blank subtraction.

In this tutorial the analysis information was just copied directly from [patRoonData]. The `generateAnalysisInfo()` function can be used to generate such a table for your own sample analyses. This function scans a given directory for MS data files and automatically fills in the `path` and `analysis` columns from this information. In addition, you can pass replicate group and blank information to this function. Example:

```{r genAnaInfo}
generateAnalysisInfo(patRoonData::exampleDataPath(), groups = c(rep("solvent-pos", 3), rep("standard-pos", 3)),
                     blanks = "solvent")
```

> **_NOTE_** Of course nothing stops you from creating a `data.frame` with analysis information manually within `R` or load the information from a _csv_ file. In fact, when you create a new project with `newProject()` you can select to generate a separate _csv_ file with analysis information (i.e. by filling in the right information in the analysis tab).

> **_NOTE_** The blanks for the solvent analyses are set to themselves. This will remove any features from the solvents later in the workflow, which is generally fine as we are usually not interested in the blanks anyway.


# Find and group features

The first step of a LC-MS non-target analysis workflow is typically the extraction of so called 'features'. While sometimes slightly different definitions are used, a feature can be seen as a single peak within an extracted ion chromatogram. For a complex sample it is not uncommon that hundreds to thousands of features can extracted. Because these large numbers this process is typically automatized nowadays.

To obtain all the features within your dataset the `findFeatures` function is used. This function requires data on the analysis information (`anaInfo` variable created earlier) and the desired algorithm that should be used. On top of that there are many more options that can significantly influence the feature finding process, hence, it is important to evaluate results afterwards.

In this tutorial we will use the [OpenMS] software to find features and stick with default parameters:

```{r features}
fList <- findFeatures(anaInfo, "openms", noiseThrInt = 1000, chromSNR = 3, chromFWHM = 5, minFWHM = 1, maxFWHM = 30)
```

After some processing time (especially for larger datasets), the next step is to _group features_. During this step, features from different analysis are grouped, optionally after alignment of their retention times. This grouping is necessary because it is common that instrumental errors will result in (slight) variations in both retention time and _m/z_ values which may complicate comparison of features between analyses. The resulting groups are referred to as **feature groups** and are crucial input for subsequent workflow steps.

To group features the `groupFeatures()` function is used, which has similar argument requirements as `findFeatures()` and likewise has many more options to tune the process. 

```{r fGroups}
fGroups <- groupFeatures(fList, "openms", rtalign = TRUE)
```

## Data clean-up {#data-cleanup}

The next step is to perform some basic rule based filtering with the filter() function. As its name suggests this function has several ways to filter data. It is a so called generic function and methods exists for various data types, such as the feature groups object that was made in the previous section (stored in the the `fGroups` variable). Note that in this tutorial the `absMinIntensity` was increased to `1E5` to simplify the results.

```{r filtering}
fGroups <- filter(fGroups, preAbsMinIntensity = 100, absMinIntensity = 1E5,
                  relMinReplicateAbundance = 1, maxReplicateIntRSD = 0.75,
                  blankThreshold = 5, removeBlanks = TRUE,
                  retentionRange = NULL, mzRange = NULL)
```


The following filterings steps will be performed:

* Features are removed if their intensity is below a defined intensity threshold (set by `absMinIntensity`). This filter is an effective way to  not only remove 'noisy' data, but, for instance, can also be used to remove any low intensity features which likely miss MS/MS data.
* If a feature is not ubiquitously present in (part of) replicate analyses it will be filtered out from that replicate group. This is controlled by setting `relMinReplicateAbundance`. The value is relative, for instance, a value of `0.5` would mean that a feature needs to be present in half of the replicates. In this tutorial we use a value of `1` which means that a feature should be present in all replicate samples. This is a _very_ effective filter in removing any outliers, for instance, caused by features which don't actually represent a well defined chromatographic peak.
* Similarly, features with within a replicate group are removed if the relative standard deviation (RSD) of their intensities exceeds that of the value set by the `maxReplicateIntRSD` argument.
* Features are filtered out that do not have a significantly higher intensity than the blank intensity. This is controlled by `blankThreshold`: the given value of `5` means that the intensity of a feature needs to be at least five times higher compared to the (average) blank signal.

The `removeBlanks` argument tells will remove all blank analyses after filtering. The `retentionRange` and `mzRange` arguments are not used here, but could be used to filter out any features outside a give retention or _m/z_ range. There are many more filters: see `?filter()` for more information.

As you may have noticed quite a large part of the features are removed as a result of the filtering step. However, using the right settings is a very effective way to separate interesting data from the rest.

<!-- UNDONE: add link below to other vignette when it contains prioritization -->

The next logical step in a non-target workflow is often to perform further prioritization of data. However, this will not be necessary in this tutorial as our samples are just known standard mixtures.

To simplify processing, we only continue with the first 25 feature groups:

```{r subSetFG}
fGroups <- fGroups[, 1:25]
```

## Inspecting results

In order to have a quick peek at the results we can use the default printing method:

```{r showFG}
fGroups
```

Furthermore, the `as.data.table()` function can be used to have a look at generated feature groups and their intensities (_i.e._ peak heights) across all analyses:

```{r gTable}
head(as.data.table(fGroups))
```

An overview of group properties is returned by the `groupInfo()` method:
```{r gInfo}
head(groupInfo(fGroups))
```

Finally, we can have a quick look at our data by plotting some nice extracted ion chromatograms (EICs) for all feature groups:

```{r plotChroms,fig.width=7}
plotChroms(fGroups, colourBy = "fGroups", showFGroupRect = FALSE, showPeakArea = TRUE,
           topMost = 1, showLegend = FALSE)
```

Note that we only plot the most intense feature of a feature group here (as set by `topMost=1`). See the reference docs for many more parameters to these functions (e.g. `?plotChroms`).

# Annotation

## MS peak lists

After obtaining a good dataset with features of interest we can start moving to find their chemical identity. Before doing so, however, the first step is to extract all relevant MS data that will be used for annotation. The tutorial data was obtained with data-dependent MS/MS, so in the ideal case we can obtain both MS and MS/MS data for each feature group.

The `generateMSPeakLists()` function will perform this action for us and will generate so called _MS peak lists_ in the process. These lists are basically (averaged) spectra in a tabular form. We will use algorithms from the [mzR] package to do so:

```{r MSPeakLists}
avgPListParams <- getDefAvgPListParams(clusterMzWindow = 0.002)
mslists <- generateMSPeakLists(fGroups, "mzr", maxMSRtWindow = 5, precursorMzWindow = 4,
                              avgFeatParams = avgPListParams, avgFGroupParams = avgPListParams)
```

Note that we lowered the `clusterMzWindow` value to _0.002_. This window is used during averaging to cluster similar _m/z_ values together. In general the better the resolution of your MS instrument, the lower the value can be set.

Similar to feature groups the `filter()` generic function can be used to clean up the peak lists afterwards:

```{r MSPLF}
mslists <- filter(mslists, relMSMSIntThr = 0.02, topMSMSPeaks = 10)
```

Here, all MS/MS mass peaks with intensities below 2% are removed and from the remaining peaks no more than the ten most intense are retained.

## Formula calculation

Using the data from the MS peak lists generated during the previous step we can generate a list of formula candidates for each feature group which is based on measured _m/z_ values, isotopic patterns and presence of MS/MS fragments. In this tutorial we will use this data as an extra hint to score candidate chemical structures generated during the next step. The command below will use [GenForm] to perform this step. Again running this code may take some time.

```{r formulas}
formulas <- generateFormulas(fGroups, mslists, "genform", relMzDev = 5, adduct = "[M+H]+", elements = "CHNOPSCl",
                             oc = FALSE, calculateFeatures = TRUE, featThresholdAnn = 0.75)
```

Note that you need to change the elements parameter to this function to make sure that formulae with sulfur and chloride (S/Cl) are also accepted. It is highly recommended to limit the elements (by default it is just C, H, N, O and P) as this can significantly reduce processing time and improbable formula candidates. In this tutorial we already knew which compounds to expect so the choice was easy, but often a good guess can be made in advance.

> **_NOTE_** The `generateFormulas()` function returns an object that contains formula candidates assigned for each feature group. In the above call the `calculateFeatures` argument is set to `TRUE`: by doing so formulae are first calculated for individual features within a feature group. These results are then used to generate a consensus candidate formula list for the complete feature group. During this process any outliers (defined by `featThresholdAnn`) are automatically removed. In contrast, setting `calculateFeatures` to `FALSE` will calculate formulae directly for feature groups (by using MS peak lists that have been averaged for the whole group). This will be significantly faster, but might produce (slightly) less accurate results.

## Compound identification

<!-- UNDONE Move MF installation to front? -->

Now it is time to actually see what compounds we may be dealing with. In this tutorial we will use [MetFrag] to come up with a list of possible candidates structures for each feature group. Before we can start you have to make sure that MetFrag and the PubChemLite library can be found by `patRoon`. Please see the [Handbook] for installation instructions. 

Then `generateCompounds()` is used to execute MetFrag and generate the `compounds`.

```{r compounds,eval=doCompounds}
compounds <- generateCompounds(fGroups, mslists, "metfrag", method = "CL",
                               dbRelMzDev = 5, fragRelMzDev = 5, fragAbsMzDev = 0.002,
                               adduct = "[M+H]+", database = "pubchemlite", maxCandidatesToStop = 2500)
```

While `generateCompounds()` is running a list of candidate compound structures will be downloaded for every feature group and ranked according to various scoring parameters.

See `?generateCompounds()` for more information on possible databases and many other parameters that can be set.

> **_NOTE_** This is often one of the most time consuming steps during the workflow. For this reason you should always take care to prioritize your data before running this function!

Finally we use the `addFormulaScoring()` function to improve ranking of candidates by incorporating the formula calculation data from the previous step.

```{r fscoring,eval=doCompounds}
compounds <- addFormulaScoring(compounds, formulas, updateScore = TRUE)
```

```{r noint,include=FALSE,eval=!doCompounds}
compounds <- NULL
```

## Inspecting results

Similar as feature groups we can quickly peek at some results:

```{r annRes,eval=doCompounds}
mslists
formulas
compounds

as.data.table(mslists)
as.data.table(formulas)[, 1:7] # only show first columns for clarity
as.data.table(compounds)[, 1:5] # only show first columns for clarity
```

```{r plotSpectrum,fig.width=6,fig.height=3.5,eval=doCompounds}
plotSpectrum(mslists, "M215_R333_2275", MSLevel = 2)
plotSpectrum(formulas, 1, "M109_R192_158", MSPeakLists = mslists)
plotSpectrum(compounds, 1, "M120_R268_286", mslists, plotStruct = TRUE)
```

# Reporting

The last step of the workflow is typically reporting data: during this step all the collected data is transformed to graphical plots (`reportPDF()` and `reportHTML()`) or tabular csv data (`reportCSV()`).

```{r eval=FALSE}
reportCSV(fGroups, path = "report", formulas = formulas, compounds = compounds, MSPeakLists = mslists,
          components = NULL)
reportHTML(fGroups, path = "report", formulas = formulas, compounds = compounds, MSPeakLists = mslists,
           components = NULL, reportPlots = c("chord", "venn", "upset", "eics", "formulas"),
           selfContained = FALSE, openReport = TRUE)
```

The output of `reportHTML()` can be viewed [here](../examples/report.html).

Note that these functions can be called at any time during the workflow. This may be especially useful if you want evaluate results during optimization or exploring the various algorithms and their parameters. In this case you can simply cherry pick the data that you want to report, for instance:

```{r eval=FALSE}
# only report feature groups (i.e. the bare minimum)
reportCSV(fGroups, path = "report", reportFeatures = FALSE)

# report formulas. Note that MSPeakLists (mslists variable) are required for formula/compound reporting
reportHTML(fGroups, path = "report", formulas = formulas, MSPeakLists = mslists)
```

```{r echo=FALSE,eval=pkgdown::in_pkgdown()}
# do the actual reporting here

# ugly work around for nested rmarkdown call made by reportHTML; based on https://gist.github.com/jennybc/1f747c5bb84aa9be9c3c
tempF <- tempfile(); tempScript <- tempfile(fileext = ".R")
save(fGroups, formulas, compounds, mslists, file = tempF)
writeLines(sprintf('
library(patRoon)
setwd("%s")
load("%s")
options("patRoon.cache.fileName" = "%s")
reportHTML(fGroups, path = "../docs/examples", formulas = formulas, compounds = compounds,
           components = NULL, MSPeakLists = mslists, optimizePng = TRUE, openReport = FALSE)
', gsub("\\", "/", getwd(), fixed = TRUE), gsub("\\", "/", tempF, fixed = TRUE), getOption("patRoon.cache.fileName")), con = tempScript)
# devtools::clean_source(tempScript, quiet = TRUE)
callr::rscript(tempScript, show = FALSE)
```

# Final script

In the previous sections the different parts of the processing script were discussed and where necessary modified. As a reference, the final script look similar ot this:

```{r scriptPost,code=readLines("script-post.R"),eval=FALSE}
```


<!-- UNDONE: next steps section -->


```{r child="shared/_refs.Rmd"}
```
# Installation

`patRoon` depends on various other software tools to perform the non-target analysis workflow steps and to implement various other functionality. Most of these dependencies are automatically installed when you install the `patRoon` R package, however, some may need to be manually installed and/or configured.

> **_NOTE_**  It is highly recommended to perform installation steps in a 'clean' `R` session to avoid errors when installing or upgrading packages. As such it is recommended to close all open (R Studio) sessions and open a plain R console to perform the installation. 
 
## Automatic installation (Windows only)

An installation script is provided that automatically installs and configures all dependencies and finally installs `patRoon` itself. At this moment, this script only works with Microsoft Windows. You don't have to install anything else to use it, simply open `R` and execute these commands:

```{r eval=FALSE}
source("https://raw.githubusercontent.com/rickhelmus/patRoon/master/install_patRoon.R")
installPatRoon()
```

A simple text based wizard will start and asks you what to install and how to do it. You can re-run this installer at any time, for instance, if something went wrong or you want to install additional dependencies.

## Docker image

Docker images are provided to easily install a reproducible environment with `R`, `patRoon` and nearly all of its dependencies. This section assumes you have a basic understanding of [Docker] and have it installed. If not, please refer to the many guides available on the Internet. The Docker images of `patRoon` were originally only used for automated testing, however, since these contain a complete working environment of `patRoon` they are also suitable for using the software. They come with all external dependencies (except ProteoWizard), `R` dependencies and `MetFrag` libraries. Furthermore, the Docker image also contains [RStudio] server, which makes using `patRoon` even easier.

Below are some example shell commands on how to run the image.

```{bash, eval=FALSE}
# run an interactive R console session
docker run --rm -it patroonorg/patroonrs

# run a linux shell, from which R can be launched
docker run --rm -it patroonorg/patroonrs bash

# run rstudio server, accessible from localhost:8787
# login with rstudio/yourpasswordhere
docker run --rm -p 8787:8787 -u 0 -e PASSWORD=yourpasswordhere patroonorg/patroonrs /init

# same as above, but mount a local directory (~/myvolume) as local volume so it can be used for persistent storage
# please ensure that ~/myvolume exists!
docker run --rm -p 8787:8787 -u 0 -e PASSWORD=yourpasswordhere -v ~/myvolume:/home/rstudio/myvolume patroonorg/patroonrs /init
```

Note that the first two commands run as the default user `rstudio`, while the last two as `root`. The last commands launch [RStudio] server. You can access it by browsing to `localhost:8787` and logging in with user `rstudio` and the password defined by the `PASSWORD` variable from the command (`yourpasswordhere` in the above example). The last command also links a local volume in order to obtain persistence of files in the container's home directory. The Docker image is based on the excellent work from the [rocker project](https://www.rocker-project.org/). For more information on RStudio related options see their documentation for the [RStudio image].


## Manual installation

The manual installation is for users who don't use Windows or Docker, prefer to do a manual installation or simply want to know what happens behind the scenes. The manual installation consists of three phases:

1. Installing some prerequisite `R` packages
2. Install and configure (non-`R`) dependencies
3. Install `patRoon`

### R prerequisites

When installing `patRoon` Windows users have the option to install from a customized ([miniCRAN]) repository (`patRoonDeps`). This repository provides a central repository for `patRoon` and all its R packages. An advantage is that installation will be faster and you will not need to install [Rtools]. Note that you will need to have the latest `R` version installed in order to use this repository.

When you decide to use the `patRoonDeps` repository you can simply _skip_ this step. **Otherwise** (i.e. you will use regular repositories instead), the following will install the `R` dependencies.

```{r eval=FALSE}
install.packages(c("BiocManager", "remotes"))
BiocManager::install("CAMERA")

# needed for finding/grouping features with KPIC2
BiocManager::install("ropls")
remotes::install_github("rickhelmus/KPIC2") 

# needed for RAMClustR for componentization
remotes::install_github("cbroeckl/RAMClustR")

# needed for nontarget for componentization
remotes::install_github("blosloos/nontargetData")
remotes::install_github("blosloos/nontarget")

# needed for cliqueMS for componentization
remotes::install_github("rickhelmus/cliqueMS") # note: repository is a fork with minor bug fixes

# only needed for Bruker DataAnalysis integration
# NOTE: a fork is installed that fixes some issues with the latest R versions
remotes::install_github("BSchamberger/RDCOMClient")

# only when using the R interface (not recommended by default)
remotes::install_github("c-ruttkies/MetFragR/metfRag")

# needed for feature quality calculations with MetaClean
BiocManager::install(c("BiocStyle", "Rgraphviz")) 
remotes::install_github("KelseyChetnik/MetaClean")
```

Note that only the `CAMERA` installation is mandatory, the rest involves installation of _optional_ packages. If you are unsure which you need then you can always install the packages at a later stage.

### Other dependencies

Depending on which functionality is used, the following external dependencies may need to be installed:

Software                            | Remarks
----------------------------------- | -----------------------------------------------------
[Java JDK][JavaJDK]                 | **Mandatory** for e.g. plotting structures and using MetFrag.
[Rtools]                            | Necessary on Window and when `patRoon` is _not_ installed from `patRoonDeps`.
[ProteoWizard]                      | Needed for automatic data-pretreatment (e.g. data file conversion and centroiding, Bruker users may use DataAnalysis integration instead).
[OpenMS]                            | Recommended. Used for e.g. finding and grouping features.
[MetFrag CL][MetFragCL]             | Recommended. Used for annotation with MetFrag.
[MetFrag CompTox DB][CompTox-dl]    | Database files necessary for usage of the [CompTox] database with MetFrag. Note that a recent version  of MetFrag (>=2.4.5) is required. Note that the lists with additions for [smoking metadata][CompTox-smoke] and [wastewater metadata][CompTox-WW] are also supported.
[MetFrag PubChemLite DB][PCLite-dl] | Database file needed to use [PubChemLite][PCLite-paper] with MetFrag.
[SIRIUS]                            | For obtaining feature data and formula and/or compound annotation.
[BioTransformer]                    | For prediction of transformation products. See the [BioTransformer] page for installation details. If you have trouble compiling the jar file you can download it from [here](https://github.com/rickhelmus/patRoonDeps/raw/master/ext/biotransformer-3.0.0.jar).
[SAFD]                              | For finding features with [SAFD]. Please follow all the installation on the [SAFD webpage][SAFD].
[OpenBabel]                         | Used in some cases for suspect screening (e.g. to calculate molecular masses for suspects with only InChI information). Otherwise optional.
[pngquant]                          | Used to reduce size of HTML reports, definitely optional.

Most of these packages are optional and only needed if their algorithms are used during the workflow. If you are unsure which are needed, simply skip them for now and install them at a later stage.

After installation of these dependencies you may need to configure their filepaths (OpenMS and ProteoWizard are usually found automatically). To configure the file locations you should set some global package options with the `options()` R function, for instance:

```{r, eval=FALSE}
options(patRoon.path.pwiz = "C:/ProteoWizard") # location of ProteoWizard installation folder
options(patRoon.path.SIRIUS = "C:/sirius-win64-3.5.1") # location where SIRIUS was extracted
options(patRoon.path.OpenMS = "/usr/local/bin") # directory with the OpenMS binaries
options(patRoon.path.pngquant = "~/pngquant") # directory containing pngquant binary
options(patRoon.path.MetFragCL = "~/MetFrag2.4.5-CL.jar") # full location to the jar file
options(patRoon.path.MetFragCompTox = "C:/CompTox_17March2019_SelectMetaData.csv") # full location to desired CompTox CSV file
options(patRoon.path.MetFragPubChemLite = "~/PubChemLite_14Jan2020_tier0.csv") # full location to desired PubChemLite CSV file
options(patRoon.path.BioTransformer = "~/biotransformer/biotransformer-3.0.1.jar")
options(patRoon.path.obabel = "C:/Program Files/OpenBabel-3.0.0") # directory with OpenBabel binaries
```

These commands have to be executed everytime you start a new R session (e.g. as part of your script). However, it is probably easier to add them to your `~/.Rprofile` file so that they are executed automatically when you start R. If you don't have this file yet you can simply create it yourself (for more information see e.g. [this SO answer](https://stackoverflow.com/a/46819910)).

### patRoon installation

Finally, it is time to install `patRoon` itself. As mentioned before, Windows users (who have the latest `R` version) can install `patRoon` and all its package dependencies from the `patRoonDeps` repository:

```{r eval=FALSE}
# install from patRoonDeps (only Windows with latest R version)
install.packages("patRoon", repos = "https://rickhelmus.github.io/patRoonDeps/", type = "binary")
```

**Otherwise**, installation occurs directly from GitHub:

```{r eval=FALSE}
remotes::install_github("rickhelmus/patRoon")
```

Optional example data is installed via GitHub:

```{r eval=FALSE}
# optional example data
remotes::install_github("rickhelmus/patRoonData")
```


Afterwards, you can run the `verifyDependencies()` function to see if `patRoon` can find all its dependencies (you may need to restart R beforehand):

```{r eval=FALSE}
patRoon::verifyDependencies()
```

```{r child=file.path(vignDir, "shared", "_refs.Rmd")}
```
---
title: "patRoon handbook"
author: "Rick Helmus"
date: "`r Sys.Date()`"
bibliography: ["refs.bib"]
biblio-style: "springer-basic-brackets-no-et-al"
link-citations: true
header-includes:
    - \usepackage{fvextra}
    - \DefineVerbatimEnvironment{Highlighting}{Verbatim}{breaklines,commandchars=\\\{\}}
---

```{r include=FALSE}
vignDir <- normalizePath("..", winslash = "/")
```

```{r child="intro.Rmd"}
```

```{r child="setup.Rmd"}
```
# Generating workflow data

## Introduction

### Workflow functions

Each step in the non-target workflow is performed by a function that performs the heavy lifting of a workflow step behind the scenes and finally return the results. An important goal of `patRoon` is to support multiple algorithms for each workflow step, hence, when such a function is called you have to specify which algorithm you want to use. The available algorithms and their characteristics will be discussed in the next sections. An overview of all functions involved in generating workflow data is shown in the table below.

Workflow step         | Function                                          | Output S4 class
--------------------- | ------------------------------------------------- | ------------
Data pre-treatment    | `convertMSFiles()`, `recalibrarateDAFiles()`      | -
Finding features      | `findFeatures()`                                  | `features`
Grouping features     | `groupFeatures()`                                 | `featureGroups`
Suspect screening     | `screenSuspects()`                                | `featureGroupsScreening` 
Componentization      | `generateComponents()`                            | `components`
MS peak lists         | `generateMSPeakLists()`                           | `MSPeakLists`
Formula annotation    | `generateFormulas()`                              | `formulas`
Compound annotation   | `generateCompounds()`                             | `compounds`

### Workflow output

The output of each workflow step is stored in objects derived from so called S4 classes. Knowing the details about the S4 class system of `R` is generally not important when using `patRoon` (and well written resources are available if you want to know more). In brief, usage of this class system allows a general data format that is used irrespective of the algorithm that was used to generate the data. For instance, when features have been found by [OpenMS] or [XCMS] they both return the same data format.

Another advantage of the S4 class system is the usage of so called _generic functions_. To put simply: a generic function performs a certain task for different types of data objects. A good example is the `plotSpectrum()` function which plots an (annotated) spectrum from data of MS peak lists or from formula or compound annotation:

```{r plotSpectrum,eval=FALSE}
# mslists, formulas, compounds contain results for MS peak lists and
# formula/compound annotations, respectively.

plotSpectrum(mslists, ...) # plot raw MS spectrum
plotSpectrum(formulas, ...) # plot annotated spectrum from formula annotation data
plotSpectrum(compounds, ...) # likewise but for compound annotation.
```

### Overview of all functions and their output

```{r workflowFuncsClasses,echo=FALSE,out.width="100%",fig.cap="**Workflow functions and output classes.**"}

makeRefURL <- function(node, page, label, tt) glue::glue("'{ node }' [URL=\"https://rickhelmus.github.io/patRoon/reference/{ page }.html\" label=<<U>{ label }</U>> tooltip=\"{ tt }\"]", node = node, page = page, label = label, tt = tt)
makeRefURLClass <- function(class) makeRefURL(class, paste0(class, "-class"), class, paste(class, "class"))
makeRefURLFunc <- function(func, page = func) makeRefURL(func, page, paste0(func, "()"), paste0(func, "()"))
makeRefURLPage <- function(node, page) makeRefURL(node, page, node, node)

URLs <- paste0(makeRefURLClass(c("transformationProducts", "featureGroupsScreening", "features", "featureGroups",
                                 "MSPeakLists", "formulas", "compounds", "components")),
               makeRefURLFunc(c("findFeatures", "groupFeatures", "generateTPs", "generateMSPeakLists",
                                "generateFormulas", "generateCompounds", "generateComponents")),
               makeRefURLFunc(c("screenSuspects", "annotateSuspects", "generateAnalysisInfo", "addFormulaScoring",
                                "selectIons", "convertToSuspects", "convertToMFDB"),
                              c("suspect-screening", "featureGroupsScreening-class", "analysis-information",
                                "compounds-class", "featureGroups-class", "transformationProducts-class",
                                "generics")),
               makeRefURLPage(c("Analysis information", "Suspect list", "database"),
                              c("analysis-information", "suspect-screening", "generateCompoundsMetFrag")), collapse = "\n")

plotGV(sprintf("
digraph workflow {
  newrank = true
  graph [ rankdir = TB, compound = true, style = invis, splines=ortho ]
  node [ fixedsize = true,
         width = 3.4,
         height = 1,
         fontsize = 24,
         style = filled ]

  subgraph cluster_suspAnn {
    graph [style = invis; margin=20 ]
    'Suspect list' -> screenSuspects -> featureGroupsScreening -> annotateSuspects
    annotateSuspects -> featureGroupsScreening
  }

  subgraph cluster_fGroups {
    graph [style = invis; margin=25 ]
    'Raw MS data' -> 'Analysis information' -> 'findFeatures' -> 'features' -> 'groupFeatures' -> 'featureGroups'
    generateAnalysisInfo -> 'Analysis information'
  }

  subgraph cluster_tps {
    graph [style = invis ]
    Parents -> generateTPs -> transformationProducts
    transformationProducts -> convertToSuspects
    convertToSuspects -> 'Suspect list' [constraint=none]
    transformationProducts -> convertToMFDB -> database
  }

  subgraph cluster_MSPL {
    graph [style = invis ]
    generateMSPeakLists -> MSPeakLists
  }

  subgraph cluster_form {
    graph [style = invis; margin=30 ]
    generateFormulas -> formulas
  }
  
  subgraph cluster_comp {
    graph [style = invis; margin=30 ]
    generateCompounds -> compounds
  }

  subgraph cluster_compon {
    graph [style = invis margin=25 ]
    generateComponents -> components
    components -> selectIons
    selectIons -> featureGroups
  }

  featureGroups -> generateMSPeakLists
  featureGroups -> generateFormulas
  featureGroups -> generateCompounds
  featureGroups -> generateComponents
  featureGroups -> screenSuspects [constraint=none]
  MSPeakLists -> generateFormulas
  MSPeakLists -> generateCompounds
  MSPeakLists -> annotateSuspects [style=dashed, constraint=none]
  formulas -> annotateSuspects [style=dashed, constraint=none]
  compounds -> annotateSuspects [style=dashed, constraint=none]
  formulas -> addFormulaScoring
  compounds -> addFormulaScoring
  addFormulaScoring -> compounds
  transformationProducts -> generateComponents [constraint=none, style=dashed]
  database -> generateCompounds [constraint=none, style=dashed]

  { rank=same; 'findFeatures'; Parents; 'Suspect list' }
  { rank=same; generateMSPeakLists; generateComponents }
  { rank=same; addFormulaScoring; compounds }
  { rank=same; 'Analysis information'; generateAnalysisInfo }
                  
  Parents, 'Suspect list', 'Raw MS data', database [shape=cylinder; fillcolor=cadetblue1]
  generateAnalysisInfo, findFeatures, groupFeatures, generateTPs, convertToSuspects, convertToMFDB, screenSuspects,
    annotateSuspects, generateMSPeakLists, generateFormulas, generateCompounds, generateComponents,
    selectIons, addFormulaScoring [ shape=Mrecord; fillcolor=\"#FFE6CC\" ]
  transformationProducts, featureGroupsScreening, 'Analysis information', features, featureGroups,
    MSPeakLists, formulas, compounds, components [ shape=octagon; fillcolor=\"#D5E8D4\" ]
  
  %s
}", URLs), width = 920, height = 800)
```

The next sections in this chapter will further detail on how to actually perform the non-target workflow steps to generate data. The transformation product screening workflows are discussed in [a separate chapter](#TPs).

## Preparations

### Data pre-treatment

Prior to performing the actual non-target data processing workflow some preparations often need to be made. Often data has to be pre-treated, for instance, by converting it to an open format that is usable for subsequent workflow steps or to perform mass re-calibration. Some common functions are listed below.

Task                                |  Function                | Algorithms | Supported file formats
------------------------------------| ------------------------ | -------------------------------------- | -------
Conversion                          | `convertMSFiles()`       | [OpenMS], [ProteoWizard], DataAnalysis | All common (algorithm dependent)
Advanced (e.g. spectral filtering)  | `convertMSFiles()`       | [ProteoWizard]                         | All common
Mass re-calibration                 | `recalibrarateDAFiles()` | DataAnalysis                           | Bruker

The `convertMSFiles()` function supports conversion between many different file formats typically used in non-target analysis. Furthermore, other pre-treatment steps are available (e.g. centroiding, filtering) when the [ProteoWizard] algorithm is used. For an overview of these functionalities see the [MsConvert documentation](http://proteowizard.sourceforge.net/tools/msconvert.html). Some examples:

```{r convert,eval=FALSE}
# Converts a single mzXML file to mzML format
convertMSFiles("standard-1.mzXML", to = "mzML", algorithm = "openms")

# Converts all Thermo files with ProteoWizard (the default) in the analyses/raw
# directory and stores the mzML files in analyses/raw. Afterwards, only MS1
# spectra are retained.
convertMSFiles("analyses/raw", "analyses/mzml", from = "thermo",
               centroid = "vendor", filters = "msLevel 1")
```

> **_NOTE_** Most algorithms further down the workflow require the _mzML_ or _mzXML_ file format and additionally require that mass peaks have been centroided. When using the ProteoWizard algorithm (the default), centroiding by vendor algorithms is generally recommended (i.e. by setting `centroid="vendor"` as shown in the above example).

When Bruker MS data is used it can be automatically re-calibrated to improve its mass accuracy. Often this is preceeded by calling the `setDAMethod()` function to set a DataAnalysis method to all files in order to configure automatic re-calibration. The `recalibrarateDAFiles()` function performs the actual re-calibration. The `getDACalibrationError()` function can be used at anytime to request the current calibration error of each analysis. An example of these functions is shown below.

```{r brukerCalib,eval=FALSE}
# anaInfo is a data.frame with information on analyses (see next section)
setDAMethod(anaInfo, "path/to/DAMethod.m") # configure Bruker files with given method that has automatic calibration setup
recalibrarateDAFiles(anaInfo) # trigger re-calibration for each analysis
getDACalibrationError(anaInfo) # get calibration error for each analysis (NOTE: also shown when previous function is finished)
```

### Analysis information {#anaInfo}

The final bits of preparation is constructing the information for the analyses that need to be processed. In `patRoon` this is referred to as the _analysis information_ and often stored in a variable `anaInfo` (of course you are free to choose a different name!). The analysis information should be a `data.frame` with the following columns:

* **path**: the directory path of the file containing the analysis data
* **analysis**: the name of the analysis. This should be the file name _without_ file extension.
* **group**: to which _replicate group_ the analysis belongs. All analysis which are replicates of each other get the same name.
* **blank**: which replicate group should be used for blank subtraction.
* **conc** (optional, advanced) A numeric value describing the concentration or any other value for which the intensity in this sample may correlate, for instance, dilution factor, sampling time etc. This column is only required when you want to obtain quantitative information (e.g. concentrations) using the `as.data.table()` method function (see `?featureGroups` for more information).

The `generateAnalysisInfo()` function can be used to (semi-)automatically generate a suitable `data.frame` that contains all the required information for a set of analysis. For, instance:

```{r genAnaInfo,eval=TRUE}
# Take example data from patRoonData package (triplicate solvent blank + triplicate standard)
generateAnalysisInfo(paths = patRoonData::exampleDataPath(),
                     groups = c(rep("solvent-pos", 3), rep("standard-pos", 3)),
                     blanks = "solvent-pos")
```

(Note that for the example data the `patRoonData::exampleAnalysisInfo()` function can also be used.)

Alternatively, the `newProject()` function discussed in the next section can be used to interactively construct this information.

### Automatic project generation with newProject() {#newProject}

The previous sections already highlighted some steps that have to be performed prior to starting a new non-target analysis workflow: data pre-treatment and gathering information on the analysis. Most of the times you will put this and other `R` code a script file so you can recall what you have done before (i.e. reproducible research).

The `newProject()` function can be used to setup a new project. When you run this function it will launch a small tool (see screenshot below) where you can select your analyses and configure the various workflow steps which you want to execute (e.g. data pre-treatment, finding features, annotation etc). After setting everything up the function will generate a template script which can easily be edited afterwards. In addition, you have the option to create a new RStudio project, which is advantegeous as it neatly seperates your data processing work from the rest.

```{r newp,echo=FALSE,out.width="450px"}
if (knitr::is_html_output()) {
    knitr::include_graphics(knitr::image_uri(file.path(vignDir, "newp.png")), error = FALSE)
} else
    knitr::include_graphics(file.path(vignDir, "newp.png"))
```

> **_NOTE_** At the moment `newProject()` _only_ works with [RStudio].


## Features

Collecting features from the analyses consists of finding all features, grouping them across analyses (optionally after retention time alignment), and if desired suspect screening: 

```{r featWorkflow,echo=FALSE,out.width="75%"}
plotGV("
digraph Features {
  graph [ rankdir = LR, compound = true ]
  node [ shape = box,
         fixedsize = true,
         width = 2.2,
         height = 1,
         fontsize = 18,
         fillcolor = darkseagreen1,
         style = filled ]

    'Find features' -> 'Group features'
    'Suspect screening'
    'Group features' -> 'Suspect screening' [style=dashed]
}", height = 70, width = 500)
```

### Finding and grouping features

Several algorithms are available for finding features. These are listed in the table below alongside their usage and general remarks.

Algorithm        | Usage                                       | Remarks
---------------- | ------------------------------------------- | --------------
[OpenMS]         | `findFeatures(algorithm = "openms", ...)`   | Uses the [FeatureFinderMetabo] algorithm
[XCMS]           | `findFeatures(algorithm = "xcms", ...)`     | Uses `xcms::xcmsSet()` function
[XCMS] (import)  | `importFeatures(algorithm = "xcms", ...)`   | Imports an existing `xcmsSet` object
[XCMS3]          | `findFeatures(algorithm = "xcms3", ...)`    | Uses `xcms::findChromPeaks()` from the new XCMS3 interface
[XCMS3] (import) | `importFeatures(algorithm = "xcms3", ...)`  | Imports an existing `XCMSnExp` object
[enviPick]       | `findFeatures(algorithm = "envipick", ...)` | Uses `enviPick::enviPickwrap()`
[KPIC2]          | `findFeatures(algorithm = "kpic2", ...)`    | Uses the [KPIC2] `R` package
[KPIC2] (import) | `importFeatures(algorithm = "kpic2", ...)`  | Imports features from [KPIC2]
[SIRIUS]         | `findFeatures(algorithm = "sirius", ...)`   | Uses [SIRIUS] to find features
[SAFD]           | `findFeatures(algorithm = "safd", ...)`     | Uses the [SAFD] algorithm (experimental)
DataAnalysis     | `findFeatures(algorithm = "bruker", ...)`   | Uses Find Molecular Features from DataAnalysis (Bruker only)

Most often the performance of these algorithms heavily depend on the data and parameter settings that are used. Since obtaining a good feature dataset is crucial for the rest of the workflow, it is highly recommended to experiment with different settings (this process can also be automated, see [the feature optimization section](#fOpt) for more details). Some common parameters to look at are listed in the table below. However, there are many more parameters that can be set, please see the reference documentation for these (e.g. `?findFeatures`).

Algorithm        | Common parameters
---------------- | ---------------------------------------------------------------------------------
[OpenMS]         | `noiseThrInt`, `chromSNR`, `chromFWHM`, `mzPPM`, `minFWHM`, `maxFWHM` (see `?findFeatures`)
[XCMS] / [XCMS3] | `peakwidth`, `mzdiff`, `prefilter`, `noise` (assuming default `centWave` algorithm, see `?findPeaks.centWave` / `?CentWaveParam`)
[enviPick]       | `dmzgap`, `dmzdens`, `drtgap`, `drtsmall`, `drtdens`, `drtfill`, `drttotal`, `minpeak`, `minint`, `maxint` (see `?enviPickwrap`)
[KPIC2]          | `kmeans`, `level`, `min_snr` (see `?findFeatures` and `?getPIC` / `?getPIC.kmeans`)
[SIRIUS]         | The `sirius` algorithm is currently parameterless
[SAFD]           | `mzRange`, `maxNumbIter`, `resolution`, `minInt` (see `?findFeatures`)
DataAnalysis     | See _Find_ -> _Parameters..._ -> _Molecular Features_ in DataAnalysis.

> **_NOTE_** Support for SAFD is still experimental and some extra work is required to set everything up. Please see the reference documentation for this algorithm (`?findFeatures`).

> **_NOTE_** DataAnalysis feature settings have to be configured in DataAnalysis prior to calling `findFeatures()`.

Similarly, for grouping features across analyses several algorithms are supported.

Algorithm        | Usage                                                    | Remarks
---------------- | -------------------------------------------------------- | -------------------------------------
[OpenMS]         | `groupFeatures(algorithm = "openms", ...)`               | Uses the [FeatureLinkerUnlabeled] algorithm (and [MapAlignerPoseClustering] for retention alignment)
[XCMS]           | `groupFeatures(algorithm = "xcms", ...)`                 | Uses `xcms::group()` `xcms::retcor()` functions
[XCMS] (import)  | `importFeatureGroupsXCMS(...)`                           | Imports an existing `xcmsSet` object.
[XCMS3]          | `groupFeatures(algorithm = "xcms3", ...)`                | Uses `xcms::groupChromPeaks()` and `xcms::adjustRtime()` functions
[XCMS3] (import) | `importFeatureGroupsXCMS3(...)`                          | Imports an existing `XCMSnExp` object.
[KPIC2]          | `groupFeatures(algorithm = "kpic2", ...)`                | Uses the [KPIC2] package
[KPIC2] (import) | `importFeatureGroupsKPIC2(...)`                          | Imports a `PIC set` object
[SIRIUS]         | `groupFeatures(anaInfo, algorithm = "sirius")`           | Finds _and_ groups features with [SIRIUS]
ProfileAnalysis  | `importFeatureGroups(algorithm = "brukerpa", ...)`       | Import `.csv` file exported from Bruker ProfileAnalysis
TASQ             | `importFeatureGroups(algorithm = "brukertasq", ...)`     | Imports a _Global result table_ (exported to Excel file and then saved as `.csv` file)

> **_NOTE_**: Grouping features with the `sirius` algorithm will perform both finding and grouping features with [SIRIUS]. This algorithm cannot work with features from another algorithm.

Just like finding features, each algorithm has their own set of parameters. Often the defaults are a good start but it is recommended to have look at them. See `?groupFeatures` for more details.

When using the [XCMS] algorithms both the 'classical' interface and latest `XCMS3` interfaces are supported. Currently, both interfaces are mostly the same regarding functionalities and implementation. However, since future developments of XCMS are primarily focused the latter this interface is recommended.

Some examples of finding and grouping features are shown below.

```{r feat,eval=FALSE}
# The anaInfo variable contains analysis information, see the previous section

# Finding features
fListOMS <- findFeatures(anaInfo, "openms") # OpenMS, with default settings
fListOMS2 <- findFeatures(anaInfo, "openms", noiseThrInt = 500, chromSNR = 10) # OpenMS, adjusted minimum intensity and S/N
fListXCMS <- findFeatures(anaInfo, "xcms", ppm = 10) # XCMS
fListXCMSImp <- importFeatures(anaInfo, "xcms", xset) # import XCMS xcmsSet object
fListXCMS3 <- findFeatures(anaInfo, "xcms3", CentWaveParam(peakwidth = c(5, 15))) # XCMS3
fListEP <- findFeatures(anaInfo, "envipick", minint = 1E3) # enviPick
fListKPIC2 <- findFeatures(anaInfo, "kpic2", kmeans = TRUE, level = 1E4) # KPIC2
fListSIRIUS <- findFeatures(anaInfo, "sirius") # SIRIUS

# Grouping features
fGroupsOMS <- groupFeatures(fListOMS, "openms") # OpenMS grouping, default settings
fGroupsOMS2 <- groupFeatures(fListOMS2, "openms", rtalign = FALSE) # OpenMS grouping, no RT alignment
fGroupsOMS3 <- groupFeatures(fListXCMS, "openms", maxGroupRT = 6) # group XCMS features with OpenMS, adjusted grouping parameter
# group enviPick features with XCMS3, disable minFraction
fGroupsXCMS <- groupFeatures(fListEP, "xcms3",
                             xcms::PeakDensityParam(sampleGroups = analInfo$group,
                                                    minFraction = 0))
# group with KPIC2 and set some custom grouping/aligning parameters
fGroupsKPIC2 <- groupFeatures(fListKPIC2, "kpic2", groupArgs = list(tolerance = c(0.002, 18)),
                              alignArgs = list(move = "loess"))
fGroupsSIRIUS <- groupFeatures(anaInfo, "sirius") # find/group features with SIRIUS
```

### Suspect screening {#suspscr}

> **_NOTE_**: you may need to install [OpenBabel], for instance, when only InChI data is available for mass calculation.

After features have been grouped a so called suspect screening step may be performed to find features that may correspond to suspects within a given suspect list. The `screenSuspects()` function is used for this purpose, for instance:


```{r suspmz,eval=FALSE}
suspects <- data.frame(name = c("1H-benzotriazole", "N-Phenyl urea", "2-Hydroxyquinoline"),
                       mz = c(120.0556, 137.0709, 146.0600))
fGroupsSusp <- screenSuspects(fGroups, suspects)
```

#### Suspect list format

The example above has a very simple suspect list with just three compounds. The format of the suspect list is quite flexible, and can contain the following columns:

- `name`: The name of the suspect. Mandatory and should be unique and file-name compatible (if not, the name will be automatically re-named to make it compatible).
- `rt`: The retention time in seconds. Optional. If specified any feature groups with a different retention time will not be considered to match suspects.
- `mz`, `SMILES`, `InChI`, `formula`, `neutralMass`: _at least_ one of these columns must hold data for each suspect row. The `mz` column specifies the ionized mass of the suspect. If this is not available then data from any of the other columns is used to determine the suspect mass.
- `adduct`: The adduct of the suspect. Optional. Set this if you are sure that a suspect should be matched by a particular adduct ion and no data in the `mz` column is available.
- `fragments_mz` and `fragments_formula`: optional columns that may assist [suspect annotation](#suspAnn).

In most cases a suspect list is best made as a `csv` file which can then be imported with e.g. the `read.csv()` function. This is exactly what happen when you specify a suspect list when using the `newProject()` function.

Quite often, the ionized masses are not readily available and these have to be calculated. In this case, data in any of the `SMILES`/`InChI`/`formula`/`neutralMass` columns should be provided. Whenever possible, it is _strongly_ recommended to fill in `SMILES` column (or `InChI`), as this will assist [annotation](#suspAnn). Applying this to the above example:

```{r suspSMI,eval=runData}
suspects <- data.frame(name = c("1H-benzotriazole", "N-Phenyl urea", "2-Hydroxyquinoline"),
                       SMILES = c("[nH]1nnc2ccccc12", "NC(=O)Nc1ccccc1", "Oc1ccc2ccccc2n1"))
fGroupsSusp <- screenSuspects(fGroups, suspects, adduct = "[M+H]+")
```

Since suspect matching now occurs by the neutral mass it is required that the adduct information for the feature groups are set. This is done either by setting the `adduct` function argument to `screenSuspects` or by [feature group adduct annotations](#incorpAdductIso).

Finally, when the adduct is known for a suspect it can be specified in the suspect list:

```{r susp,eval=FALSE}
# Aldicarb is measured with a sodium adduct.
suspects <- data.frame(name = c("1H-benzotriazole", "N-Phenyl urea", "Aldicarb"),
                       SMILES = c("[nH]1nnc2ccccc12", "NC(=O)Nc1ccccc1", "CC(C)(C=NOC(=O)NC)SC"),
                       adduct = c("[M+H]+", "[M+H]+", "[M+Na]+"))
fGroupsSusp <- screenSuspects(fGroups, suspects)
```

To summarize:

* If a suspect has data in the `mz` column it will be directly matched with the _m/z_ value of a feature group.
* Otherwise, if the suspect has data in the `adduct` column, the `m/z` value for the suspect is calculated from its neutral mass and the adduct and then matched with the `m/z` of a feature group.
* Otherwise, suspects and feature groups are matched by their the neutral mass.

The `fragments_mz` and `fragments_formula` columns in the suspect list can be used to specify known fragments for a suspect, which can help [suspect annotation](#suspAnn). The former specifies the ionized _m/z_ of known MS/MS peaks, whereas the second specifies known formulas. Multiple values can be given by separating them with a semicolon:

```{r suspFrags,eval=FALSE}
suspects <- data.frame(name = c("1H-benzotriazole", "N-Phenyl urea", "2-Hydroxyquinoline"),
                       SMILES = c("[nH]1nnc2ccccc12", "NC(=O)Nc1ccccc1", "Oc1ccc2ccccc2n1"),
                       fragments_formula = c("C6H6N", "C6H8N;C7H6NO", ""),
                       fragments_mz = c("", "", "118.0652"))
```

#### Removing feature groups without hits

Note that any feature groups that were not matched to a suspect are _not_ removed by default. If you want to remove these, you can use the `onlyHits` parameter:

```{r suspOH,eval=FALSE}
fGroupsSusp <- screenSuspects(fGroups, suspects, onlyHits = TRUE) # remove any non-hits immediately
```

The advantage of removing non-hits is that it may significantly reduce the complexity of your dataset. On the other hand, retaining all features allows you to mix a full non-target analysis with a suspect screening workflow. The `filter()` function (discussed [here](#filtering)) can also be used to remove feature groups without a hit at a later stage.

#### Combining screening results

The `amend` function argument to `screenSuspects` can be used to combine screening results from different suspect lists. 

```{r suspAmend,eval=FALSE}
fGroupsSusp <- screenSuspects(fGroups, suspects)
fGroupsSusp <- screenSuspects(fGroupsSusp, suspects2, onlyHits = TRUE, amend = TRUE)
```

In this example the suspect lists defined in `suspects` and `suspects2` are both used for screening. By setting `amend=TRUE` the original screening results (i.e. from `suspects`) are preserved. Note that `onlyHits` should only be set in the final call to `screenSuspects` to ensure that all feature groups are screened.


## Componentization {#componentization}

In `patRoon` _componentization_ refers to grouping related feature groups together in components. There are different methodologies to generate components:

* Similarity on chromatographic elution profiles: feature groups with similar chromatographic behaviour which are assuming to be the same chemical compound (e.g. adducts or isotopologues).
* Homologous series: features with increasing _m/z_ and retention time.
* Intensity profiles: features that follow a similar intensity profile in the analyses.
* MS/MS similarity: feature groups with similar MS/MS spectra are clustered.
* Transformation products: Components are formed by grouping feature groups that have a parent/transformation product relationship. This is further discussed in [its own chapter](#TPs).

The following algorithms are currently supported:

Algorithm               | Usage                                              | Remarks
----------------------- | -------------------------------------------------- | --------------------------------------------
[CAMERA]                | `generateComponents(algorithm = "camera", ...)`    | Clusters feature groups with similar chromatographic elution profiles and annotate by known chemical rules (adducts, isotopologues, in-source fragments).
[RAMClustR]             | `generateComponents(algorithm = "ramclustr", ...)` | As above.
[cliqueMS]              | `generateComponents(algorithm = "cliquems", ...)`  | As above, but using _feature components_.
[OpenMS]                | `generateComponents(algorithm = "openms", ...)`    | As above. Uses [MetaboliteAdductDecharger].
[nontarget]             | `generateComponents(algorithm = "nontarget", ...)` | Uses the [nontarget] R package to perform unsupervised homologous series detection.
Intensity clustering    | `generateComponents(algorithm = "intclust", ...)`  | Groups features with similar intensity profiles across analyses by hierarchical clustering.
MS/MS clustering        | `generateComponents(algorithm = "specclust", ...)` | Clusters feature groups with similar MS/MS spectra.
Transformation products | `generateComponents(algorithm = "tp", ...)`        | Discussed in [its own chapter](#TPs).


### Features with similar chromatographic behaviour {#componentsChrom}

Isotopes, adducts and in-source fragments typically result in detection of multiple mass peaks by the mass spectrometer for a single chemical compound. While some feature finding algorithms already try to collapse (some of) these in to a single feature, this process is often incomplete (if performed at all) and it is not uncommon that multiple features will describe the same compound. To overcome this complexity several algorithms can be used to group features that undergo highly similar chromatographic behavior but have different _m/z_ values. Basic chemical rules are then applied to the resulting components to annotate adducts, in-source fragments and isotopologues, which may be highly useful for general identification purposes.

Note that some algorithms were primarily designed for datasets where features are generally present in the majority of the analyses (as is relatively common in metabolomics). For environmental analyses, however, this is often not the case. For instance, consider the following situation with three feature groups that chromatographically overlap and therefore could be considered a component:

Feature group | _m/z_      | analysis 1 | analysis 2 | analysis 3
------------- | ---------- | ---------- | ---------- | ----------
\#1           | 100.08827  | Present    | Present    | Absent
\#2           | 122.07021  | Present    | Present    | Absent
\#3           | 138.04415  | Absent     | Absent     | Present

Based on the mass differences from this example a cluster of `[M+H]+`, `[M+Na]+` and `[M+K]+` could be assumed. However, no features of the first two feature groups were detected in the third sample analysis, whereas the third feature group wasn't detected in the first two sample analysis. Based on this it seems unlikely that feature group _#3_ should be part of the component.

For the algorithms that operate on a 'feature group level' ([CAMERA] and [RAMClustR]), the `relMinReplicates` argument can be used to remove feature groups from a component that are not abundant. For instance, when this value is _0.5_ (the default), and all the features of a component were detected in four different replicate groups in total, then only those feature groups are kept for which its features were detected in at least two different replicate groups (_i.e._ half of four).

Another approach to reduce unlikely adduct annotations is to use algorithms that operate on a 'feature level' ([cliqueMS] and [OpenMS]). These algorithms generate components for each sample analysis individually. The 'feature components' are then merged by a consensus approach where unlikely annotations are removed (the algorithm is described further in the reference manual, `?generateComponents`).

Each algorithm supports many different parameters that may significantly influence the (quality of the) output. For instance, care has to be taken to avoid 'over-clustering' of feature groups which do not belong in the same component. This is often easily visible since the chromatographic peaks poorly overlap or are shaped differently. The `checkComponents` function (discussed [here](#intReview)) can be used to quickly verify componentization results. For a complete listing all arguments see the reference manual (e.g. `?generateComponents`).

Once the components with adduct and isotopes annotations are generated this data can be [used to prioritize and improve the workflow](#incorpAdductIso).

Some example usage is shown below.

```{r chromComps,eval=FALSE}
# Use CAMERA with defaults
componCAM <- generateComponents(fGroups, "camera", ionization = "positive")

# CAMERA with customized settings
componCAM2 <- generateComponents(fGroups, "camera", ionization = "positive",
                                 extraOpts = list(mzabs = 0.001, sigma = 5))

# Use RAMClustR with customized parameters
componRC <- generateComponents(fGroups, "ramclustr", ionization = "positive", hmax = 0.4,
                               extraOptsRC = list(cor.method = "spearman"),
                               extraOptsFM = list(ppm.error = 5))

# OpenMS with customized parameters
componOpenMS <- generateComponents(fGroups, "openms", ionization = "positive", chargeMax = 2,
                                   absMzDev = 0.002)

# cliqueMS with default parameters
componCliqueMS <- generateComponents(fGroups, "cliquems", ionization = "negative")
```


### Homologues series

Homologues series can be automatically detected by interfacing with the [nontarget] R package. Components are made from feature groups that show increasing _m/z_ and retention time values. Series are first detected within each replicate group. Afterwards, series from all replicates are linked in case (partial) overlap occurs and this overlap consists of the _same_ feature groups (see figure below). Linked series are then finally merged if this will not cause any conflicts with other series: such a conflict typically occurs when two series are not only linked to each other.

```{r HS,echo=FALSE,fig.cap="**Linking of homologues series** top: partial overlap and will be linked; bottom: no linkage due to different feature in overlapping series.",out.width="100%"}
plotGV("
digraph Homologs {
  graph [ rankdir = LR, compound = true, style = invis ]
  node [ shape = oval,
         fixedsize = true,
         width = 2.3,
         height = 0.8,
         fontsize = 25,
         fillcolor = darkseagreen1,
         style = filled ]

  subgraph cluster3 {
    I [fillcolor=skyblue]; J [fillcolor=skyblue]; K [fillcolor=skyblue]
    G -> H -> I -> J -> K [style=invis]
  }

  subgraph cluster4 {
    _I [label=I, fillcolor=skyblue]; X [fillcolor=indianred1]; _K [label=K, fillcolor=skyblue]
    _I -> X -> _K -> L [style=invis]
  }

  H -> _I [ltail=cluster3, lhead=cluster4, style=invis]
  
  I -> _I [constraint=false, style=dashed, arrowhead=none]
  K -> _K [constraint=false, style=dashed, arrowhead=none]
  
  subgraph cluster1 {
    C [fillcolor=skyblue]; D [fillcolor=skyblue]
    A -> B -> C -> D [style=invis]
  }

  subgraph cluster2 {
    _C [label=C, fillcolor=skyblue]; _D [label=D, fillcolor=skyblue]
    _C -> _D -> E -> F [style=invis]
  }
  
  B -> _C [ltail=cluster1, lhead=cluster2, style=invis]

  C -> _C [constraint=false, style=dashed, arrowhead=none]
  D -> _D [constraint=false, style=dashed, arrowhead=none]
}", height = 250, width = 600)
```

The series that are linked can be interactively explored with the `plotGraph()` function (discussed [here](#vis_feat_ann)).

Common function arguments to `generateComponents()` are listed below.

Argument             | Remarks
-------------------- | -----------------------------------------------------------------------
`ionization`         | Ionization mode: `"positive"` or `"negative"`. Not needed if [adduct annotations](#incorpAdductIso) are available.
`rtRange`, `mzRange` | Retention and _m/z_ increment range. Retention times can be negative to allow series with increasing _m/z_ values and decreasing retention times.
`elements`           | Vector with elements to consider.
`rtDev`, `absMzDev`  | Maximum retention time and _m/z_ deviation.
`...`                | Further arguments passed to the `homol.search()` function.

```{r compsNT,eval=FALSE}
# default settings
componNT <- generateComponents(fGroups, "nontarget", ionization = "positive")

# customized settings
componNT2 <- generateComponents(fGroups, "nontarget", ionization = "positive",
                                elements = c("C", "H"), rtRange = c(-60, 60))
```

### Intensity and MS/MS similarity {#compClust}

<!-- UNDONE: refs? -->

The previous componentization methods utilized chemical properties to relate features. The two componentization algorithms described in this section use a statistical approach based on hierarchical clustering. The first algorithm normalizes all feature intensities and then clusters features with similar intensity profiles across sample analyses together. The second algorithm compares all MS/MS spectra from all feature groups, and then uses hierarchical clustering to generate components from feature groups that have a high MS/MS spectrum similarity.

Some common arguments to `generateComponents()` are listed below. It is recommended to test various settings (especially for `method`) to optimize the clustering results.

Argument                                      | Algorithm   | Default          | Remarks
--------------------------------------------- | ----------- | ---------------- | ----------------------------------------------
`method`                                      | All         | `"complete"`     | Clustering method. See `?hclust`
`metric`                                      | `intclust`  | `"euclidean"`    | Metric used to calculate the distance matrix. See `?daisy`.
`normFunc`                                    | `intclust`  | `max`            | Function used to normalize data. Feature intensities within a feature group are divided by the result of when this function is called with their intensity values.
`average`                                     | `intclust`  | `TRUE`           | Whether intensities of replicates should first be averaged.
`MSPeakLists`                                 | `specclust` | -                | The [MS peak lists] object used for spectral similarity calculations
`specSimParams`                               | `specclust` | `getDefSpecSimParams()` | Parameters used for [spectral similarity calculation](#specSim).
`maxTreeHeight`, `deepSplit`, `minModuleSize` | All         | `1`, `TRUE`, `1` | Used for dynamic cluster assignment. See `?cutreeDynamicTree`.

The components are generated by automatically assigning clusters using the [dynamicTreeCut] `R` package. However, the cluster assignment can be performed manually or with different parameters, as is demonstrated below.

The resulting components are stored in an object from the `componentsIntClust` or `componentsSpecClust` S4 class, which are both derived from the `componentsClust` class (which in turn is derived from the `components` class). Several methods are defined that can be used on these objects to re-assign clusters, perform plotting operations and so on. Below are some examples. For plotting see the relevant [visualization section](#plotClust). More info can be found in the reference manual (e.g. `?componentsIntClust`, `?componentsSpecClust` and `?componentsClust`).

```{r intClust,eval=FALSE}
# generate intensity profile components with default settings
componInt <- generateComponents(fGroups, "intclust")

# manually re-assign clusters
componInt <- treeCut(componInt, k = 10)

# automatic re-assignment of clusters (adjusted max tree height)
componInt <- treeCutDynamic(componInt, maxTreeHeight = 0.7)

# MS/MS similarity components
componMSMS <- generateComponents(fGroups, "specclust", MSPeakLists = mslists)
```


```{r child=file.path(vignDir, "shared", "_refs.Rmd")}
```

## Incorporating adduct and isotopic data {#incorpAdductIso}

With mass spectrometry it is common that multiple _m/z_ values are detected for a single compound. These may be different adducts (e.g. `[M+H]+`, `[M+Na]+`, `[M-H]-`), the different isotopes of the molecule or a combination thereof. When multiple _m/z_ values are measured for the same compound, the feature finding algorithm may yield a distinct feature for each, which adds complexity to the data. In [the previous section](#componentsChrom) it was discussed how componentization can help to find feature groups that belong to the same adduct and/or isotope clusters. This section explains how this data can be used to simplify the feature dataset. Furthermore, this section also covers adduct annotations for feature groups which may improve and simplify the general workflow.

### Selecting features with preferential adducts/isotopes

The `selectIons` function forms the bridge between feature group and componentization data. This function uses the adduct and isotope annotations to select _preferential_ feature groups. For adduct clusters this means that only the feature group that has a preferential adduct (e.g. `[M+H]+`) is kept while others (e.g. `[M+Na]+`) are removed. If none of the adduct annotations are considered preferential, the most intense feature group is kept instead. For isotopic clusters typically only the feature group with the monoisotopic mass (i.e. _M0_) is kept.

The behavior of `selectIons` is configurable with the following parameters:

Argument         | Remarks
---------------- | -------------------------------------------------------
`prefAdduct`     | The _preferential adduct_. Usually `"[M+H]+"` or `"[M-H]-"`.
`onlyMonoIso`    | If `TRUE` and a feature group is with isotopic annotations then it is only kept if it is monoisotopic.
`chargeMismatch` | How charge mismatches between adduct and isotope annotations are dealt with. Valid options are `"isotope"`, `"adduct"`, `"none"` or `"ignore"`. See the reference manual for `selectIons` for more details.

In case componentization did not lead to an adduct annotation for a feature group it will never be removed and simply be annotated with the preferential adduct. Similarly, when no isotope annotations are available and `onlyMonoIso=TRUE`, the feature group will not be removed.

Although `selectIons` operates fairly conservative, it is still recommended to verify the componentization results in advance, for instance with the `checkComponents` function [discussed here](#intReview). Furthermore, the next subsection explains how adduct annotations can be corrected manually if needed.

An example usage is shown below.

```{r selectIons,eval=runData}
fGroupsSel <- selectIons(fGroups, componCAM, "[M+H]+")
```

### Setting adduct annotations for feature groups

The `adducts()` function can be used to obtain a character vector with adduct annotations for each feature group. When no adduct annotations are available it will simply return an empty character vector.

When the `selectIons` function is used it will automatically add adduct annotations based on the componentization data. In addition, the `adducts()<-` function can be used to manually add or change adduct annotations.

```{r featAdductAnn,eval=runData}
adducts(fGroups) # no adduct annotations
adducts(fGroupsSel)[1:5] # adduct annotations set by selectIons()

adducts(fGroupsSel)[3] <- "[M+Na]+" # modify annotation
adducts(fGroupsSel)[1:5] # verify
```

> **_NOTE_** Adduct annotations are always available with [sets workflows](#setsWorkflow).

### Using adduct annotations in the workflow {#useAdductAnn}

When feature groups have adduct annotations available this may simplify and improve the workflow. The `adduct` and `ionization` arguments used for suspect screening, formula/compound annotation and some componentization algorithms do not have to be set anymore, since this data can be obtained from the adduct annotations. Furthermore, these algorithms may improve their results, since the algorithms are now able to use adduct information for each feature group individually, instead of assuming that all feature groups have the same adduct.

## Annotation

The annotation consists of collecting MS peak lists and then formula and/or compound annotation:

```{r annWorkflow,echo=FALSE,out.width="50%"}
plotGV("
digraph Annotation {
  graph [ rankdir = LR ]
  node [ shape = box,
         fixedsize = true,
         width = 2.3,
         height = 1,
         fontsize = 18,
         fillcolor = darkseagreen1,
         style = filled ]

    'MS peak lists' -> 'Formula annotation'
    'MS peak lists' -> 'Compound annotation'
    'Formula annotation':n -> 'Suspect annotation':n [style = dotted]
    'Compound annotation':s -> 'Suspect annotation':s [style = dotted]
    'Formula annotation':e -> 'Compound annotation':e [style = dotted, constraint = false]
}", height = 120, width = 500)
```

Note that compound annotation is normally not dependent upon formula annotation. However, formula data can be used to improve ranking of candidates afterwards by the `addFormulaScoring()` function, which will be discussed later in this section. Furthermore, suspect annotation is not mandatory, and may use data from peak lists, formulae and/or comounds.

### MS peak lists

Algorithm    | Usage                                                   | Remarks
------------ | ------------------------------------------------------- | -----------------------------------------------------
[mzR]        | `generateMSPeakLists(algorithm = "mzr", ...)`           | Uses [mzR] for spectra retrieval. Recommended default.
DataAnalysis | `generateMSPeakLists(algorithm = "bruker", ...)`        | Loads data after automatically generating MS and MS/MS spectra in DataAnalysis
DataAnalysis FMF | `generateMSPeakLists(algorithm = "brukerfmf", ...)` | Uses spectra from the  _find molecular features_ algorithm.

The recommended default algorithm is `mzr`: this algorithm is generally faster and is not limited to a vendor data format as it will read the open `mzML` and `mzXML` file formats. On the other hand, when DataAnalysis is used with Bruker data the spectra can be automatically background subtracted and there is no need for file conversion. Note that the `brukerfmf` algorithm only works when `findFeatures()` was called with the `bruker` algorithm.

When `generateMSPeakists()` is called it will

1. Find all MS and MS/MS spectra that 'belong' to a feature. For MS spectra this means that all spectra close to the retention time of a feature will be collected. In addition, for MS/MS normally only spectra will be considered that have a precursor mass close to that of the feature (however, this can be disabled for data that was recorded with data independent acquisition (DIA, MS^E, bbCID, ...)).
2. Average all MS and MS/MS spectra to produce peak lists for each feature.
3. Average all peak lists for features within the same group.

Data from either (2) or (3) is used for subsequent annotation steps. Formula calculation can use either (as a trade-off between possibly more accurate results by outlier removal _vs_ speed), whereas compound annotation will always use data from (3) since annotating single features (as opposed to their groups) would take a very long time.

There are several common function arguments to `generateMSPeakLists()` that can be used to optimize its behaviour:

<!-- UNDONE: add topMost for Bruker when it's there -->

Argument                             | Algorithm(s)          | Remarks
------------------------------------ | --------------------- | ----------------------------------------------------------------------
`maxMSRtWindow`                      | `mzr`, `bruker`       | Maximum time window +/- the feature retention time (in seconds) to collect spectra for averaging. Higher values may significantly increase processing times.
`precursorMzWindow`                  | `mzr`                 | Maximum precursor _m/z_ search window to find MS/MS spectra. Set to `NULL` to disable (i.e. for DIA experiments).
`topMost`                            | `mzr`                 | Only retain feature data for no more than this amount analyses with highest intensity. For instance, a value of _1_ will only keep peak lists for the feature with highest intensity in a feature group.  
`bgsubtr`                            | `bruker`              | Perform background subtraction (if the spectra type supports this, e.g. MS and bbCID)
`minMSIntensity`, `minMSMSIntensity` | `bruker`, `brukerfmf` | Minimum MS and MS/MS intensity. Note that DataAnalysis reports many zero intensity peaks so a value of at least _1_ is recommended.
`MSMSType`                           | `bruker`              | The type of spectra that should be used for MSMS: `"BBCID"` for bbCID experiments, otherwise `"MSMS"` (the default).

In addition, several parameters can be set that affect spectral averaging. These parameters are passed as a `list` to the `avgFeatParams` (`mzr` algorithm only) and `avgFGroupParams` arguments, which affect averaging of feature and feature group data, respectively. Some typical parameters include:

* `clusterMzWindow`: Maximum _m/z_ window used to cluster mass peaks when averaging. The better the MS resolution, the lower this value should be.
* `topMost`: Retain no more than this amount of most intense mass peaks. Useful to filter out 'noisy' peaks.
* `minIntensityPre` / `minIntensityPost`: Mass peaks below this intensity will be removed before/after averaging.

See `?generateMSPeakLists` for all possible parameters.

A suitable list object to set averaging parameters can be obtained with the `getDefAvgPListParams()` function.

```{r avgMSPLParamas,eval=FALSE}
# lower default clustering window, other settings remain default
avgPListParams <- getDefAvgPListParams(clusterMzWindow = 0.001)

# Apply to both feature and feature group averaging
plists <- generateMSPeakLists(fGroups, "mzr", avgFeatParams = avgPListParams, avgFGroupParams = avgPListParams)
```

### Formulae

Formulae can be automatically calculated for all features using the `generateFormulas()` function. The following algorithms are currently supported:

Algorithm    | Usage                                          | Remarks
------------ | ---------------------------------------------- | ------------------------------------------------
[GenForm]    | `generateFormulas(algorithm = "genform", ...)` | Bundled with `patRoon`. Reasonable default.
[SIRIUS]     | `generateFormulas(algorithm = "sirius", ...)`  | Requires MS/MS data.
DataAnalysis | `generateFormulas(algorithm = "bruker", ...)`  | Requires FMF features (i.e. `findFeatures(algorithm = "bruker", ...)`). Uses _SmartFormula_ algorithms.

Calculation with [GenForm] is often a good default. It is fast and basic rules can be applied to filter out obvious non-existing formulae. A possible drawback of GenForm, however, is that may become slow when many candidates are calculated, for instance, due to a relative high feature _m/z_ (e.g. >600) or loose elemental restricitions. More thorough calculation is performed with [SIRIUS]: this algorithm often yields fewer and often more plausible results. However, [SIRIUS] requires MS/MS data (hence features without will not have results) and formula prediction may not work well for compounds that structurally deviate from the training sets used by [SIRIUS]. Calculation with DataAnalysis is only possible when features are obtained with DataAnalysis as well. An advantage is that analysis files do not have to be converted, however, compared to other algorithms calculation is often relative slow.

There are two methods for formula assignment:

1. Formulae are first calculated for each individual feature within a feature group. These results are then pooled, outliers are removed and remaining formulae are assigned to the feature group (i.e. `calculateFeatures = TRUE`).
2. Formulae are directly calculated for each feature group by using group averaged peak lists (see previous section) (i.e. `calculateFeatures = FALSE`).

The first method is more thorough and the possibility to remove outliers may sometimes result in better formula assignment. However, the second method is much faster and generally recommended for large number of analyses.

By default, formulae are either calculated by _only_ MS/MS data (SIRIUS) or with both MS _and_ MS/MS data (GenForm/Bruker). The latter also allows formula calculation when no MS/MS data is present. Furthermore, with Bruker algorithms, data from both MS and MS/MS formula data can be combined to allow inclusion of candidates that would otherwise be excluded by e.g. poor MS/MS data. However, a disadvantage is that formulae needs to be calculated twice. The `MSMode` argument (listed below) can be used to customize this behaviour.

An overview of common parameters that are typically set to customize formula calculation is listed below.

Argument          | Algorithm(s)        | Remarks
----------------- | ------------------- | -------------------------------------------------------------------------------
relMzDev          | `genform`, `sirius` | The maximum relative _m/z_ deviation for a formula to be considered (in _ppm_).
elements          | `genform`, `sirius` | Which elements to consider. By default `"CHNOP"`. Try to limit possible elements as much as possible.
calculateFeatures | `genform`, `sirius` | Whether formulae should be calculated first for all features (see discussion above) (always `TRUE` with DataAnalysis).
featThresholdAnn  | All                 | Minimum relative amount (_0-1_) that a candidate formula for a feature group should be found among all annotated features (e.g. _1_ means that a candidate is only considered if it was assigned to all annotated features). 
adduct            | All                 | The adduct to consider for calculation (e.g. `"[M+H]+"`, `"[M-H]-"`, more details in the [adduct section](#adducts)). Don't set this when [adduct annotations](#incorpAdductIso) are available.
MSMode            | `genform`, `bruker` | Whether formulae should be generated only from MS data (`"ms"`), MS/MS data (`"msms"`) or both (`"both"`). The latter is default, see discussion above.
profile           | `sirius`            | Instrument profile, e.g. `"qtof"`, `"orbitrap"`, `"fticr"`.

Some typical examples:

```{r forms,eval=FALSE}
formulasGF <- generateFormulas(fGroups, mslists, "genform") # GenForm, default settings
formulasGF2 <- generateFormulas(fGroups, mslists, "genform", calculateFeatures = FALSE) # direct feature group assignment (faster)
formulasSIR <- generateFormulas(fGroups, mslists, "sirius", elements = "CHNOPSClBr") # SIRIUS, common elements for pollutant
formulasSIR2 <- generateFormulas(fGroups, mslists, "sirius", adduct = "[M-H]-") # SIRIUS, negative ionization
formulasBr <- generateFormulas(fGroups, mslists, "bruker", MSMode = "MSMS") # Only consider MSMS data (SmartFormula3D)
```

### Compounds {#compounds}

An important step in a typical non-target workflow is structural identification for features of interest. Afterall, this information may finally reveal _what_ a feature is. The first step is to find all possible structures in a database that may be assigned to the feature (based on e.g. monoisotopic mass or formula). These candidates are then scored to rank likely candidates, for instance, on correspondence with in-silico or library MS/MS spectra and environmental relevance. 

Structure assignment in `patRoon` is performed automatically for all feature groups with the `generateCompounds()` function. Currently, this function supports two algorithms:

Algorithm                    | Usage                                           | Remarks
---------------------------- | ----------------------------------------------- | ------------------------------------------------
[MetFrag]                    | `generateCompounds(algorithm = "metfrag", ...)` | Supports many databases (including custom) and scorings for candidate ranking.
[SIRIUS] with [CSI:FingerID] | `generateCompounds(algorithm = "sirius", ...)`  | Incorporates prior comprohensive formula calculations.

Compound annotation is often a relative time and resource intensive procedure. For this reason, features are not annotated individually, but instead a feature group as a whole is annotated, which generally saves significant amounts of computational requirements. Nevertheless, it is not uncommon that this is the most time consuming step in the workflow. For this reason, prioritization of features is highly important, even more so to avoid 'abusing' servers when an online database is used for compound retrieval.

Selecting the right database is important for proper candidate assignment. Afterall, if the 'right' chemical compound is not present in the used database, it is impossible to assign the correct structure. Luckily, however, several large databases such as [PubChem] and [ChemSpider] are openly available which contain tens of millions of compounds. On the other hand, these databases may also lead to many unlikely candidates and therefore more specialized (or custom databases) may be preferred. Which database will be used is dictated by the `database` argument to `generateCompounds()`, currently the following options exist:

<!-- UNDONE: add missing databases when support is added -->

Database            | Algorithm(s)            | Remarks
------------------- | ----------------------- | -----------------------
`pubchem`           | `"metfrag"`, `"sirius"` | [PubChem] is currently the largest compound database and is used by default.
`chemspider`        | `"metfrag"`             | [ChemSpider] is another large database. Requires security token from [here](http://www.chemspider.com/AboutServices.aspx) (see next section).
`comptox`           | `"metfrag"`             | The EPA [CompTox] contains many compounds and scorings relevant to environmental studies. Needs manual download (see next section).
`pubchemlite`       | `"metfrag"`             | A specialized subset of the PubChem database. Needs manual download (see next section).
`for-ident`         | `"metfrag"`             | The [FOR-IDENT] (STOFF-IDENT) database for water related substances.
`kegg`              | `"metfrag"`, `"sirius"` | The [KEGG] database for biological compounds
`hmdb`              | `"metfrag"`, `"sirius"` | The [HMDB] contains many human metabolites.
`bio`               | `"sirius"`              | Selects all supports biological databases.
`csv`, `psv`, `sdf` | `"metfrag"`             | Custom database (see next section). [CSV example][csvDB-ex].

#### Configuring MetFrag databases and scoring

Some extra configuration may be necessary when using certain databases with MetFrag. In order to use the ChemSpider database a [security token](http://www.chemspider.com/AboutServices.aspx) should be requested and set with the `chemSpiderToken` argument to `generateCompounds()`. The CompTox and PubChemLite databases need to be manually downloaded from [CompTox][CompTox-dl] (or variations with [smoking][CompTox-smoke] or [wastewater][CompTox-WW] metadata) and [PubChemLite][PCLite-dl]. The file location of this and other local databases (`csv`, `psv`, `sdf`) needs to be manually configured, see the examples below and/or `?generateCompounds` for more information on how to do this.

```{r compDB,eval=FALSE}
# PubChem: the default
compsMF <- generateCompounds(fGroups, mslists, "metfrag", adduct = "[M+H]+")

# ChemSpider: needs security token
compsMF2 <- generateCompounds(fGroups, mslists, "metfrag", database = "chemspider",
                              chemSpiderToken = "MY_TOKEN_HERE", adduct = "[M+H]+")

# CompTox: set global option to database path
options(patRoon.path.MetFragCompTox = "~/CompTox_17March2019_SelectMetaData.csv")
compsMF3 <- generateCompounds(fGroups, mslists, "metfrag", database = "comptox", adduct = "[M+H]+")

# CompTox: set database location without global option
compsMF4 <- generateCompounds(fGroups, mslists, "metfrag", database = "comptox", adduct = "[M+H]+",
                              extraOpts = list(LocalDatabasePath = "~/CompTox_17March2019_SelectMetaData.csv"))

# Same, but for custom database
compsMF5 <- generateCompounds(fGroups, mslists, "metfrag", database = "csv", adduct = "[M+H]+",
                              extraOpts = list(LocalDatabasePath = "~/mydb.csv"))
```

An example of a custom _.csv_ database can be found [here][csvDB-ex].

With MetFrag compound databases are not only used to retrieve candidate structures but are also used to obtain metadata for further ranking. Each database has its own scorings, a table with currently supported scorings can be obtained with the `compoundScorings()` function (some columns omitted):

```{r compSc,echo=FALSE}
dt <- compoundScorings()[, c("name", "metfrag", "database", "default")]
if (knitr::is_html_output())
{
    k <- knitr::kable(dt, format = "html")
    k <- kableExtra::kable_styling(k, font_size = 11)
    k <- kableExtra::scroll_box(k, extra_css = "overflow-y: auto; height: 350px;")
} else
{
    k <- knitr::kable(dt, format = "latex", longtable = TRUE, booktabs = TRUE)
    # k <- kableExtra::kable_styling(k, latex_options = "repeat_header", font_size = 7)
    # UNDONE: repeat_header doesn't work anymore --> make bug report
    k <- kableExtra::kable_styling(k, font_size = 7)
}

k
```

The first two columns contain the generic and original MetFrag naming schemes for each scoring type. While both naming schemes can be used, the generic is often shorter and harmonized with other algorithms (e.g. SIRIUS). The _database_ column specifies for which databases a particular scoring is available (empty if not database specific). Most scorings are selected by default (as specified by the _default_ column), however, this behaviour can be customized by using the `scoreTypes` argument:

```{r compCustSc,eval=FALSE}
# Only in-silico and PubChem number of patents scorings
compsMF1 <- generateCompounds(fGroups, mslists, "metfrag", adduct = "[M+H]+",
                              scoreTypes = c("fragScore" "numberPatents"))

# Custom scoring in custom database
compsMF2 <- generateCompounds(fGroups, mslists, "metfrag", adduct = "[M+H]+",
                              database = "csv",
                              extraOpts = list(LocalDatabasePath = "~/mydb.csv"),
                              scoreTypes = c("fragScore", "myScore", "myScore2"))
```

By default ranking is performed with equal weight (i.e. _1_) for all scorings. This can be changed by the `scoreWeights` argument, which should be a `vector` containing the weights for all scorings following the order of `scoreTypes`, for instance:

```{r compCustScW,eval=FALSE}
compsMF <- generateCompounds(fGroups, mslists, "metfrag", adduct = "[M+H]+",
                             scoreTypes = c("fragScore" "numberPatents"),
                             scoreWeights = c(1, 2))
```

Sometimes thousands or more structural candidates are found when annotating a feature group. In this situation processing all these candidates will too involving (especially when external databases are used). To avoid this a default cut-off is set: when the number of candidates exceed a certain amount the search will be aborted and no results will be reported for that feature group. The maximum number of candidates can be set with the `maxCandidatesToStop` argument. The default value is relative conservative, especially for local databases it may be useful to increase this number.

#### Timeout and error handling

The use of online databases has the drawback that an error may occur, for instance, as a result of a connection error. Furthermore, MetFrag typically returns an error when too many candidates are found (as set by the `maxCandidatesToStop` argument). By default processing is restarted if an error has occurred (configured by the `errorRetries` argument). Similarly, the `timeoutRetries` and `timeout` arguments can be used to avoid being 'stuck' on obtaining results, for instance, due to an unstable internet connection.

If no compounds could be assigned due to an error a warning will be issued. In this case it is best to see what went wrong by manually checking the log files, which by default are stored in the _log/metfrag_ folder.

#### Formula scoring

Ranking of candidate structures may further be improved by incorporating formula information by using the `addFormulaScoring()` function:

```{r addFormSc,eval=FALSE}
comps <- addFormulaScoring(coms, formulas, updateScore = TRUE)
```

Here, corresponding formula and explained fragments will be used to calculate a _formulaScore_ for each candidate. Note that SIRIUS candidates are already based on calculated formulae, hence, running this function on SIRIUS results is less sensable unless scoring from another formula calculation algorithm is desired.

#### Further options and parameters

There are _many_ more options and parameters that affect compound annotation. For a full overview please have a look at the reference manual (e.g. by running `?generateCompounds`).

### Suspect annotation {#suspAnn}

The data obtained during the previously described annotation steps can be used to improve a suspect screening workflow. The `annotateSuspects()` method uses the annotation data to calculate various annotation properties for each suspect, such as their rank in formula/compound candidates, which fragments from the suspect list were matched, and a _rough_ indication of the identification level according to @Schymanski2014

```{r annSusps,eval=runData}
fGroupsSusp <- annotateSuspects(fGroupsSusp, MSPeakLists = mslists,
                                formulas = formulas, compounds = compounds)
```

The calculation of identification levels is performed by a set of pre-defined rules. The `genIDLevelRulesFile()` can be used to inspect the default rules or to create your own rules file, which can subsequently passed to `annotateSuspects()` with the `IDFile` argument. See `?annotateSuspects` for more details on the file format and options. The default identification levels can be summarized as follows:

Level | Description | Rules
----- | ----------- | --------------------------------
1     | Target match               | Retention time deviates <12 seconds from suspect list. At least `3` (or all if the suspect list contains less) fragments from the suspect list must match.
2a    | Good MS/MS library match   | Suspect is top ranked in the `compounds` results. The `individualMoNAScore` is at least `0.9` and all other candidates have no MoNA library score.
3a    | Fair library match | The `individualMoNAScore` is at least 0.4.
3b    | Known MS/MS match          | At least `3` (or all if the suspect list contains less) fragments from the suspect list must match.
3c    | Good in-silico MS/MS match | The annotation MS/MS similarity (`annSimComp` column) is at least `0.7`.
4a    | Good formula MS/MS match   | Suspect is top ranked formula candidate, annotation MS/MS similarity (`annSimForm` column) is at least `0.7` and isotopic match (`isoScore`) of at least 0.5. The latter two scores are at least `0.2` higher than next best ranked candidate.
4b    | Good formula isotopic pattern match | Suspect is top ranked formula candidate and isotopic match (`isoScore`) of at least `0.9` and at least `0.2` higher than next best ranked candidate.
5     | Unknown                    | All else.

In general, the more data provided by the suspect list and to `annotateSuspects()`, the better identification level estimation works. For instance, when considering the default rules, either the `fragments_mz` or `fragments_formula` column is necessary to be able assign a `level 3b`. Similarly, the suspect list needs retention times (as well as fragment data) to be able to assign `level 1`. As you can imagine, providing the annotation workflow objects (i.e. `MSPeakLists`, `formulas`, `compounds`) to `annotateSuspects()` is necessary for calculation of most levels.

The `annotateSuspects()` function will log decisions for identification level assignments to the `log/` sub-directory in the current working directory. This is useful to inspect level assignments and especially useful when you customized any rules.

> **_NOTE_**: The current identification level rules are _only_ optimized for GenForm and MetFrag annotation algorithms.
# References

<!-- Refs are always added at the end -->
# Advanced usage {#advanced_usage}

## Adducts {#adducts}

When generating formulae and compound annotations and some other functionalities it is required to specify the adduct species. Behind the scenes, different algorithms typically use different formats. For instance, in order to specify a protonated species...

* `GenForm` either accepts `"M+H"` and `"+H"`
* `MetFrag` either accepts the numeric code `1` or `"[M+H]+"`
* `SIRIUS` accepts `"[M+H]+"`

In addition, most algorithms only accept a limited set of possible adducts, which do not necessarily all overlap with each other. The `GenFormAdducts()` and `MetFragAdducts()` functions list the possible adducts for `GenForm` and `MetFrag`, respectively.

In order to simplify the situation `patRoon` internally uses its own format and converts it automatically to the algorithm specific format when necessary. Furthermore, during conversion it checks if the specified adduct format is actually allowed by the algorithm. Adducts in `patRoon` are stored in the `adduct` S4 class. Objects from this class specify which elements are added and/or subtracted, the final charge and the number of molecules present in the adduct (e.g. _2_ for a dimer).

```{r adduct1,eval=FALSE}
adduct(add = "H") # [M+H]+
adduct(sub = "H", charge = -1) # [M-H]-
adduct(add = "K", sub = "H2", charge = -1) # [M+K-H2]-
adduct(add = "H3", charge = 3) # [M+H3]3+
adduct(add = "H", molMult = 2) # [2M+H]+
```

A more easy way to generate adduct objects is by using the `as.adduct()` function:

```{r adduct2,eval=FALSE}
as.adduct("[M+H]+")
as.adduct("[M+H2]2+")
as.adduct("[2M+H]+")
as.adduct("[M-H]-")
as.adduct("+H", format = "genform")
as.adduct(1, isPositive = TRUE, format = "metfrag")
```

In fact, the `adduct` argument to workflow functions such as `generateFormulas()` and `generateCompounds()` is automatically converted to an `adduct` class with the `as.adduct()` function if necessary:

```{r adduct3,eval=FALSE}
formulas <- generateFormulas(..., adduct = adduct(sub = "H", charge = -1))
formulas <- generateFormulas(..., adduct = "[M-H]-") # same as above
```

More details can be found in the reference manual (`?adduct` and ``?`adduct-utils` ``).

## Feature parameter optimization {#fOpt}

Many different parameters exist that may affect the output quality of feature finding and grouping. To avoid time consuming manual experimentation, functionality is provided to largely automate the optimization process. The methodology, which uses design of experiments (DoE), is based on the excellent [Isotopologue Parameter Optimization (IPO)](IPO) `R` package. The functionality of this package is directly integrated in patRoon. Some functionality was added or changed, the most important being support for other feature finding and grouping algorithms besides [XCMS] and basic optimization support for qualitative parameters. Nevertheless, the core optimization algorithms are largely untouched.

This section will introduce the most important concepts and functionality. Please see the reference manual for more information (e.g. `` ?`feature-optimization` ``).

> **_NOTE_** The [SIRIUS] and [SAFD] algorithms are currently not (yet) supported.

### Parameter sets

Before starting an optimization experiment we have to define _parameter sets_. These sets contain the parameters and (initial) numeric ranges that should be tested. A parameter set is defined as a regular `list`, and can be easily constructed with the `generateFeatureOptPSet()` and  `generateFGroupsOptPSet()` functions (for feature finding and feature grouping, respectively).

```{r fOptPSet,eval=FALSE}
pSet <- generateFeatureOptPSet("openms") # default test set for OpenMS
pSet <- generateFeatureOptPSet("openms", chromSNR = c(5, 10)) # add parameter
# of course manually making a list is also possible (e.g. if you don't want to test the default parameters)
pSet <- list(noiseThrInt = c(1000, 5000))
```

When optimizing with [XCMS] or [KPIC2] a few things have to be considered. First of all, when using the XCMS3 interface (i.e. `algorithm="xcms3"`) the underlying method that should be used for finding and grouping features and retention alignment should be set. In case these are not set default methods will be used.

```{r fOptPSetXCMS1,eval=FALSE}
pSet <- list(method = "centWave", ppm = c(2, 8))
pSet <- list(ppm = c(2, 8)) # same: centWave is default

# get defaults, but for different grouping/alignment methods
pSetFG <- generateFGroupsOptPSet("xcms3", groupMethod = "nearest", retAlignMethod = "peakgroups")
```

In addition, when optimizing feature grouping (both `XCMS` interfaces and `KPIC2`) we need to set the grouping and retention alignment parameters in two different nested lists: these are `groupArgs`/`retcorArgs` (`algorithm="xcms"`), `groupParams`/`retAlignParams` (`algorithm="xcms3"`) or `groupArgs`/`alignArgs` (`algorithm="kpic2"`).

```{r fOptPSetXCMS2,eval=FALSE}
pSetFG <- list(groupParams = list(bw = c(20, 30))) # xcms3
pSetFG <- list(retcorArgs = list(gapInit = c(0, 7))) # xcms
pSetFG <- list(groupArgs = list(mz_weight = c(0.3, 0.9))) # kpic2
```

When a parameter set has been defined it should be used as input for the `optimizeFeatureFinding()` or `optimizeFeatureGrouping()` functions.

```{r fDoOpt,eval=FALSE}
ftOpt <- optimizeFeatureFinding(anaInfo, "openms", pSet)
fgOpt <- optimizeFeatureGrouping(fList, "openms", pSetFG) # fList is an existing features object
```

Similar to `findFeatures()`, the first argument to `optimizeFeatureFinding()` should be the [analysis information](#anaInfo). Note that it is not uncommon to perform the optimization with only a subset of (representative) analyses (i.e. to reduce processing time).

```{r fDoOpt2,eval=FALSE}
ftOpt <- optimizeFeatureFinding(anaInfo[1:2, ], "openms", pSet) # only use first two analyses
```

From the parameter set a design of experiment will be automatically created. Obviously, the more parameters are specified, the longer such an experiment will take. After an experiment has finished, the optimization algorithm will start a new experiment where numeric ranges for each parameter are increased or decreased in order to more accurately find optimum values. Hence, the numeric ranges specified in the parameter set are only _initial_ ranges, and will be changed in subsequent experiments. After each experiment iteration the results will be evaluated and a new experiment will be started as long as better results were obtained during the last experiment (although there is a hard limit defined by the `maxIterations` argument).

For some parameters it is recommended or even necessary to set hard limits on the numeric ranges that are allowed to be tested. For instance, setting a minimum feature intensity threshold is highly recommended to avoid excessive processing time and potentially suboptimal results due to excessive amounts of resulting features. Configuring absolute parameter ranges is done by setting the `paramRanges` argument.

```{r fOptPRanges,eval=FALSE}
# set minimum intensity threshold (but no max)
ftOpt <- optimizeFeatureFinding(anaInfo, "openms",
                                list(noiseThrInt = c(1000, 5000), # initial testing range
                                paramRanges = list(noiseThrInt = c(500, Inf))) # never test below 500
```

Depending on the used algorithm, several default absolute limits are imposed. These may be obtained with the `getDefFeaturesOptParamRanges()` and `getDefFGroupsOptParamRanges()` functions. 

The common operation is to optimize numeric parameters. However, parameters that are not numeric (i.e. _qualitative_) need a different approach. In this case you will need to define multiple parameter sets, where each set defines a different qualitative value.

```{r fOptQual,eval=FALSE}
ftOpt <- optimizeFeatureFinding(anaInfo, "openms",
                                list(chromFWHM = c(4, 8), isotopeFilteringModel = "metabolites (5% RMS)"),
                                list(chromFWHM = c(4, 8), isotopeFilteringModel = "metabolites (2% RMS)"))
```

In the above example there are two parameter sets: both define the numeric `chromFWHM` parameter, whereas the qualitative `isotopeFilteringModel` parameter has a different value for each. Note that we had to specify the `chromFWHM` twice, this can be remediated by using the `templateParams` argument:

```{r fOptQualT,eval=FALSE}
ftOpt <- optimizeFeatureFinding(anaInfo, "openms",
                                list(isotopeFilteringModel = "metabolites (5% RMS)"),
                                list(isotopeFilteringModel = "metabolites (2% RMS)"),
                                templateParams = list(chromFWHM = c(4, 8)))
```

As its name suggests, the `templateParams` argument serves as a template parameter set, and its values are essentially combined with each given parameter set. Note that current support for optimizing qualitative parameters is relatively basic and time consuming. This is because tests are essentially repeated for each parameter set (e.g. in the above example the `chromFWHM` parameter is optimized twice, each time with a different value for `isotopeFilteringModel`).

### Processing optmization results

The results of an optimization process are stored in objects from the S4 `optimizationResult` class. Several methods are defined to inspect and visualize these results.

The `optimizedParameters()` function is used to inspect the best parameter settings. Similarly, the `optimizedObject()` function can be used to obtain the object that was created with these settings (i.e. a `features` or `featureGroups` object).

```{r fOptProc1,eval=doOpt}
optimizedParameters(ftOpt) # best settings for whole experiment
optimizedObject(ftOpt) # features object with best settings for whole experiment
```

Results can also be obtained for specific parameter sets/iterations.

```{r fOptProc,eval=doOpt}
optimizedParameters(ftOpt, 1) # best settings for first parameter set
optimizedParameters(ftOpt, 1, 1) # best settings for first parameter set and experiment iteration
optimizedObject(ftOpt, 1) # features object with best settings for first parameter set
```

The `plot()` function can be used to visualize optimization results. This function will plot graphs for results of all tested parameter pairs. The graphs can be contour, image or perspective plots (as specified by the `type` argument).

```{r fOptPlot,eval=doOpt,out.width="70%"}
plot(ftOpt, paramSet = 1, DoEIteration = 1) # contour plots for first param set/experiment
plot(ftOpt, paramSet = 1, DoEIteration = 1, type = "persp") # pretty perspective plots
```

Please refer to the reference manual for more methods to inspect optimization results (e.g. `?optimizationResult`).

## Chromatographic peak qualities {#peakQual}

The algorithms used by `findFeatures` detect chromatographic peaks automatically to find the features. However, it is common that not all detected features have 'proper' chromatographic peaks, and some features could be just noise. The [MetaClean] `R` package supports various quality measures for chromatographic peaks. The quality measures include Gaussian fit, symmetry, sharpness and others. In addition, [MetaClean] averages all feature data for each feature group and adds a few additional group specific quality measures (_e.g._ retention time consistency). Please see @Chetnik2020 for more details. The calculations are integrated into `patRoon`, and are easily performed with the `calculatePeakQualities()` generic function.

```{r calcQ,eval=runData}
fList <- calculatePeakQualities(fList) # calculate for all features
fGroups <- calculatePeakQualities(fGroups) # calculate for all features and groups
```

Most often the `featureGroups` method is only used, unless you want to filter features (discussed below) prior to grouping.

An extension in `patRoon` is that the qualities are used to calculate _peak scores_. The score for each quality measure is calculated by normalizing and scaling the values into a `0-1` range, where zero is the worst and one the best. Note that most scores are relative, hence, the values should only be used to compare features among each other. Finally, a `totalScore` is calculated which sums all individual scores and serves as a rough overall score indicator for a feature (group).

The qualities and scores are easily obtained with the `as.data.table()` function.

```{r tabQ,eval=runData}
# (limit rows/columns for clarity)
as.data.table(fList)[1:5, 26:30]
# the qualities argument is necessary to include the scores.
# valid values are: "quality", "score" or "both"
as.data.table(fGroups, qualities = "both")[1:5, 25:29]
```

The feature quality values can also be reviewed interactively with reports generated with `reportHTML` (see [Reporting]) and with `checkFeatures` (see [here](#intReview)). The `filter` function can be used filter out low scoring features and feature groups:

```{r filterQ,eval=FALSE}
# only keep features with at least 0.3 Modality score and 0.5 symmetry score
fList <- filter(fList, qualityRange = list(ModalityScore = c(0.3, Inf),
                                           SymmetryScore = c(0.5, Inf)))

# same as above
fGroups <- filter(fGroups, featQualityRange = list(ModalityScore = c(0.3, Inf),
                                                   SymmetryScore = c(0.5, Inf)))

# filter group averaged data
fGroups <- filter(fGroups, groupQualityRange = list(totalScore = c(0.5, Inf)))
```

### Applying machine learning with MetaClean

An important feature of [MetaClean] is to use the quality measures to train a machine learning model to automatically recognize 'good' and 'bad' features. `patRoon` provides a few extensions to simplify training and using a model. Furthermore, while `MetaClean` was primarily designed to work with [XCMS], the extensions of `patRoon` allow the usage of data from all the algorithms supported by `patRoon`.

The `getMCTrainData` function can be used to convert data from a [feature check session](#intReview) to training data that can be used by [MetaClean]. This allows you to use interactively select good/bad peaks. The workflow looks like this:

```{r trainMC,eval=FALSE}
# untick the 'keep' checkbox for all 'bad' feature groups
checkFeatures(fGroupsTrain, "train_session.yml")

# get train data. This gives comparable data as MetaClean::getPeakQualityMetrics()
trainData <- getMCTrainData(fGroupsTrain, "train_session.yml")

# use train data with MetaClean with MetaClean::runCrossValidation(),
# MetaClean::getEvaluationMeasures(), MetaClean::trainClassifier() etc
# --> see the MetaClean vignette for details
```

Once you have created a model with [MetaClean] it can be used with the `predictCheckFeaturesSession()` function:

```{r useModelMC,eval=FALSE}
predictCheckFeaturesSession(fGroups, "model_session.yml", model)
```

This will generate another _check session file_: all the feature groups that are considered good will be with a 'keep' state, the others without. As [described elsewhere](#intReview), the `checkFeatures` function is used to review the results from a session and the `filter` function can be used to remove unwanted feature groups. Note that `calculatePeakQualitites()` must be called before `getMCTrainData`/`predictCheckFeaturesSession` can be used.

> **_NOTE_** `MetaClean` only predicts at the feature group level. Thus, only the kept feature groups from a _feature check session_ will be used for training, and any indivual features that were marked as removed will be ignored.


## Exporting and converting feature data

The feature group data obtained during the workflow can be exported to various formats with the `export()` generic function. There are currently three formats supported: `"brukerpa"` (Bruker ProfileAnalysis), `"brukertasq"` (Bruker TASQ) and `"mzmine"` (mzMine). The former exports a 'bucket table' which can be loaded in ProfileAnalysis, the second and third export a target list that can be processed with TASQ and mzMine, respectively.

The `getXCMSSet()` function converts a `features` or `featureGroups` object to an `xcmsSet` object which can be used for further processing with [xcms]. Similarly, the `getXCMSnExp()` function can be used for conversion to an XCMS3 style `XCMSnExp` object, and the `getPICSet()` function can be used to convert `features` to [KPIC2] data.

Some examples for these functions are shown below.

```{r fExpConv,eval=FALSE}
export(fGroups, "brukertasq", out = "my_targets.csv")

# convert features to xcmsSet.
# NOTE: loadRawData should only be FALSE when the analysis data files cannot be
# loaded by the algorithm (e.g. when features were obtained with DataAnalysis and data was not exported to mz(X)ML)
xset <- getXCMSSet(fList, loadRawData = TRUE)
xsetg <- getXCMSSet(fGroups, loadRawData = TRUE) # get grouped xcmsSet

# using the new XCMS3 interface
xdata <- getXCMSnExp(fList)
xdata <- getXCMSnExp(fGroups)

# KPIC2 conversion. Like XCMS it optionally loads the raw data.
picSet <- getPICSet(fList, loadRawData = TRUE)
```

## Algorithm consensus {#consensus}

With `patRoon` you have the option to choose between several algorithms for most workflow steps. Each algorithm is typically characterized by its efficiency, robustness, and may be optimized towards certain data properties. Comparing their output is therefore advantegeuous in order to design an optimum workflow. The `consensus()` generic function will compare different results from different algorithms and returns a _consensus_, which may be based on minimal overlap, uniqueness or simply a combination of all results from involved objects. The output from the `consensus()` function is of similar type as the input types and is therefore compatible to any 'regular' further data processing operations (e.g. input for other workflow steps or plotting). Note that a consensus can also be made from objects generated by the same algorithm, for instance, to compare or combine results obtained with different parameters (e.g. different databases used for compound annotation).

The `consensus()` generic is defined for most workflow objects. Some of its common function arguments are listed below.

Argument      | Classes                                            | Remarks
------------- | -------------------------------------------------- | -------------------------------------------------------------
`obj`, `...`  | All                                                | Two or more objects (of the same type) that should be compared to generate the consensus.
`compThreshold`, `relAbundance`, `absAbundance`, `formThreshold`   | `compounds`, `formulas`, `featureGroupsComparison` | The minimum overlap (relative/absolute) for a result (feature, candidate) to be kept.
`uniqueFrom`  | `compounds`, `formulas`, `featureGroupsComparison` | Only keep _unique_ results from specified objects.
`uniqueOuter` | `compounds`, `formulas`, `featureGroupsComparison` | Should be combined with `uniqueFrom`. If `TRUE` then only results are kept which are _also_ unique between the objects specified with `uniqueFrom`.

Note that current support for generating a consensus between `components` objects is very simplistic; here results are not compared, but the consensus simply consists a combination of all the components from each object.

Generating a consensus for feature groups involves first generating a `featureGroupsComparison` object. This step is necessary since (small) deviations between retention times and/or mass values reported by different feature finding/grouping algorithms complicates a direct comparison. The comparison objects are made by the `comparison()` function, and its results can be visualized by the plotting functions discussed [in the previous chapter](#visComp).

Some examples are shown below

```{r consensus,eval=FALSE}
compoundsCons <- consensus(compoundsMF, compoundsSIR) # combine MetFrag/SIRIUS results
compoundsCons <- consensus(compoundsMF, compoundsSIR,
                           compThreshold = 1) # only keep results that overlap

fGroupComp <- comparison(fGroupsXCMS, fGroupsOpenMS, fGroupsEnviPick,
                         groupAlgo = "openms")
plotVenn(fGroupComp) # visualize overlap/uniqueness
fGroupsCons <- consensus(fGroupComp,
                         uniqueFrom = 1:2) # only keep results unique in OpenMS+XCMS
fGroupsCons <- consensus(fGroupComp,
                         uniqueFrom = 1:2,
                         uniqueOuter = TRUE) # as above, but also exclude any overlap between OpenMS/XCMS
```

## Compound clustering {#compclust}

When large databases such as [PubChem] or [ChemSpider] are used for compound annotation, it is common to find _many_ candidate structures for even a single feature. While choosing the right scoring settings can significantly improve their ranking, it is still very much possible that many candidates of potential interest remain. In this situation it might help to perform _compound clustering_. During this process, all candidates for a feature are clustered hierarchically on basis of similar chemical structure. From the resulting cluster the _maximum common substructure_ (MCS) can be derived, which represents the largest possible substructure that 'fit' in all candidates. By visual inspection of the MCS it may be possible to identify likely common structural properties of a feature.

In order to perform compound clustering the `makeHCluster()` generic function should be used. This function heavily relies on chemical fingerprinting functionality provided by [rcdk].

```{r cClust,eval=FALSE}
compounds <- generateCompounds(...) # get our compounds
compsClust <- makeHCluster(compounds)
```

This function accepts several arguments to fine tune the clustering process:

* `method`: the clustering method (e.g. `"complete"` (default), `"ward.D2"`), see `?hclust` for options
* `fpType`: finger printing type (`"extended"` by default), see `?get.fingerprint`
* `fpSimMethod`: similarity method for generating the distance method (`"tanimoto"` by default), see `?fp.sim.matrix`

For all arguments see the reference manual (`?makeHClust`).

The resulting object is of type `compoundsCluster`. Several methods are defined to modify and inspect these results:

```{r cClustObj,eval=FALSE}
# plot MCS of first cluster from candidates of M109_R116_61
plotStructure(compsClust, groupName = "M109_R116_61", 1)

# plot dendrogram
plot(compsClust, groupName = "M109_R116_61")

# re-assign clusters for a feature group
compsClust <- treeCut(compsClust, k = 5, groupName = "M109_R116_61")
# ditto, but automatic cluster determination
compsClust <- treeCutDynamic(compsClust, minModuleSize = 3, groupName = "M109_R116_61")
```

For a complete overview see the reference manual (`?compoundsCluster`).

## Basic quantitative and regression analysis

While `patRoon` is currently primarily focused on qualitative analyses, some _basic_ quantitative analysis can also be performed, for instance, to estimate concentrations of features. In fact, other types of data that may be useful for regression analysis can be set such as sample dilution factor or sampling time. The latter may, for instance, be used to isolate features with increasing or decreasing intensity. Regardless of what kind of regression analysis is performed, here we simply refer the values to be calculated as _concentrations_. In order to use this functionality, an extra column (`conc`) should be added to the [analysis information](#anaInfo), for instance:

```{r quantAnaInfo,eval=FALSE}
# obtain analysis information as usual, but add some concentrations.
# The blanks are set to NA, whereas the standards are set to increasing levels.
anaInfo <- generateAnalysisInfo(paths = patRoonData::exampleDataPath(),
                                groups = c(rep("solvent", 3), rep("standard", 3)),
                                blanks = "solvent",
                                concs = c(NA, NA, NA, 1, 2, 3))
```

For analyses with known concentrations (e.g. standards) the `conc` column should be set; for all others the value should be set to `NA`.

The `as.data.table()` function (or `as.data.frame()`) can then be used to calculate regression data and estimate concentrations:

```{r quantDT,eval=FALSE}
# use areas for quantitation and make sure that feature data is reported
# (otherwise no concentrations are calculated)
# (only relevant columns are shown for clarity)
as.data.table(fGroups, areas = TRUE, features = TRUE, regression = TRUE)
```

```{r quantDTDo,echo=FALSE,eval=runData}
as.data.table(fGroupsConc, areas = TRUE, features = TRUE, regression = TRUE)[, c("group", "conc", "RSQ", "intercept", "slope", "conc_reg")]
```

Calculated concentrations are stored in the `conc_reg` column, alongside while other regression data (i.e. `RSQ`, `intercept`, `slope` columns). To perform basic trend analysis the `RSQ` (i.e. R squared) can be used:

```{r quantRSQ,eval=FALSE}
fGroupsTab <- as.data.table(fGroups, areas = TRUE, features = FALSE, regression = TRUE)
# subset fGroups with reasonable correlation
increasingFGroups <- fGroups[, fGroupsTab[RSQ >= 0.8, group]]
```

## Fold changes {#FCCalc}

A specific statistical way to prioritize feature data is by Fold changes (FC). This is a relative simple method to quickly identify (significant) changes between two sample groups. A typical use case is to compare the feature intensities before and after an experiment.

To perform FC calculations we first need to specify its parameters. This is best achieved with the `getFCParams()` function:

```{r FCParams,eval=TRUE}
getFCParams(c("before", "after"))
```

In this example we generate a `list` with parameters in order to make a comparison between two replicate groups: `before` and `after`. Several advanced parameters are available to tweak the calculation process. These are explained in the reference manual (`?featureGroups`).

The `as.data.table` function for feature groups is used to perform the FC calculations.

```{r include=FALSE,eval=runData}
fGroupsO <- fGroups
fGroups <- filter(fGroupsUF, absMinIntensity = 1E4)
```

```{r getFC,eval=runData}
myFCParams <- getFCParams(c("solvent-pos", "standard-pos")) # compare solvent/standard
as.data.table(fGroups, FCParams = myFCParams)[, c("group", "FC", "FC_log", "PV", "PV_log", "classification")]
```

The `classification` column allows you to easily identify if and how a feature changes between the two sample groups. This can also be used to prioritize feature groups:

```{r eval=FALSE}
tab <- as.data.table(fGroups, FCParams = myFCParams)
# only keep feature groups that significantly increase or decrease
fGroupsChanged <- fGroups[, tab[classification %in% c("increase", "decrease")]$group]
```

The `plotVolcano` function can be used to visually the FC data:

```{r plotVolcano,eval=runData}
plotVolcano(fGroups, myFCParams)
```

```{r include=FALSE,eval=runData}
fGroups <- fGroupsO
```

## Caching

In `patRoon` lengthy processing operations such as finding features and generating annotation data is _cached_. This means that when you run such a calculation again (without changing any parameters), the data is simply loaded from the cache data instead of re-generating it. This in turn is very useful, for instance, if you have closed your R session and want to continue with data processing at a later stage.

The cache data is stored in a [sqlite] database file. This file is stored by default under the name `cache.sqlite` in the current working directory (for this reason it is very important to always restore your working directory!). However, the name and location can be changed by setting a global package option:

```{r cacheFile,eval=FALSE}
options(patRoon.cache.fileName = "~/myCacheFile.sqlite")
```

For instance, this might be useful if you want to use a shared cache file between projects.

After a while you may see that your cache file can get quite large. This is especially true when testing different parameters to optimize your workflow. Furthermore, you may want to clear the cache after you have updated `patRoon` and want to make sure that the latest code is used to generate the data. At any point you can simply remove the cache file. A more fine tuned approach which doesn't wipe all your cached data is by using the `clearCache()` function. With this function you can selectively remove parts of the cache file. The function has two arguments: `what`, which specifies what should be removed, and `file` which specifies the path to the cache file. The latter only needs to be specified if you want to manage a different cache file.

In order to figure what is in the cache you can run `clearCache()` without any arguments:

```{r clearCache}
clearCache()
```

Using this output you can re-run the function again, for instance:

```{r clearCache2,eval=FALSE}
clearCache("featuresOpenMS")
clearCache(c("featureGroupsOpenMS", "formulasGenForm")) # clear multiple
clearCache("OpenMS") # clear all with OpenMS in name (ie partial matched)
clearCache("all") # same as simply removing the file
```

## Parallelization

Some steps in the non-target screening workflow are inherently computationally intensive. To reduce computational times `patRoon` is able to perform _parallelization_ for most of the important functionality. This is especially useful if you have a modern system with multiple CPU cores and sufficient RAM.

For various technical reasons several parallelization techniques are used, these can be categorized as _parallelization of `R` functions_ and _multiprocessing_. The next sections describe both parallelization approaches in order to let you optimize the workflow.

### Parellization of R functions

Several functions of `patRoon` support parallelization.

Function                 | Purpose                             | Remarks
------------------------ | ----------------------------------- | ------------------------
`findFeatures`           | Obtain feature data                 | Only `envipick` and `kpic2` algorithms
`generateComponents`     | Generate components                 | Only `cliquems` algorithm.
`optimizeFeatureFinding`, `optimizeFeatureGrouping` | Optimize feature finding/grouping parameters | Discussed [here](#fOpt).
`calculatePeakQualities` | Calculate feature (group) qualities | Discussed [here](#peakQual).

The parallelization is achieved with the [future] and [future.apply] `R` packages. To enable parallelization of these functions the `parallel` argument must be set to `TRUE` and the future framework must be properly configured in advance. For example:

```{r eval=FALSE}
# setup three workers to run in parallel
future::plan("multisession", workers = 3)

# find features with enviPick in parallel
fList <- findFeatures(anaInfo, "envipick", parallel = TRUE)
```

It is important to properly configure the right future plan. Please see the documentation of the [future] package for more details.

### Multiprocessing

`patRoon` relies on several external (command-line) tools to generate workflow data. These commands may be executed in _parallel_ to reduce computational times ('multiprocessing'). The table below outlines the tools that are executed in parallel.

Tool                  | Used by                                      | Notes
--------------------- | -------------------------------------------- | ---------------------------------
`msConvert`           | `convertMSFiles(algorithm="pwiz", ...)`      |
`FileConverter`       | `convertMSFiles(algorithm="openms", ...)`    |
`FeatureFinderMetabo` | `findFeatures(algorithm="openms", ...)`      |
`julia`               | `findFeatures(algorithm="safd", ...)`        |
`SIRIUS`              | `findFeatures(algorithm="sirius", ...)`      |
`MetaboliteAdductDecharger` | `generateComponents(algorithm="openms", ...)` |
`GenForm`             | `generateFormulas(agorithm="genform", ...)`  |
`SIRIUS`              | `generateFormulas(agorithm="sirius", ...)`, `generateCompounds(agorithm="sirius", ...)`  | Only if `splitBatches=TRUE`
`MetFrag`             | `generateCompounds(agorithm="metfrag", ...)` |
`pngquant`            | `reportHTML(...)`                            | Only if `optimizePng=TRUE`
`BioTransformer`      | `generateTPs(algorithm = "biotransformer")`  | Disabled by default (see `?generateTPs` for details).

Multiprocessing is either performed by executing processes in the background with the [processx] `R` package (_classic interface_) or by futures, which were introduced in the previous section. An overview of the characteristics of both parallelization techniques is shown below.

`classic`                                           |  `future`
--------------------------------------------------- | -------------------------------------------------------
requires little or no configuration                 | configuration needed to setup
works with all tools                                | doesn't work with `pngquant` and slower with `GenForm`
only supports parallelization on the local computer | allows both local and cluster computing

Which method is used is controlled by the `patRoon.MP.method` package option. Note that `reportHTML()` will always use the classic method for `pngquant`.

#### Classic multiprocessing interface

The classic interface is the 'original' method implemented in `patRoon`, and is therefore well tested and optimized. It is easier to setup, works well with all tools, and is therefore the default method. It is enabled as follows:

```{r eval=FALSE}
options(patRoon.MP.method = "classic")
```

The number of parallel processes is configured through the `patRoon.MP.maxProcs` option. By default it is set to the number of available CPU cores, which results usually in the best performance. However, you may want to lower this, for instance, to keep your computer more responsive while processing or limit the RAM used by the data processing workflow.

```{r eval=FALSE}
options(patRoon.MP.maxProcs = 2) # do not execute more than two tools in parallel. 
```

This will change the parallelization for the complete workflow. However, it may be desirable to change this for only a part the workflow. This is easily achieved with the `withOpt()` function.

```{r eval=FALSE}
# do not execute more than two tools in parallel.
options(patRoon.MP.maxProcs = 2)

# ... but execute up to four GenForm processes
withOpt(MP.maxProcs = 4, {
    formulas <- generateFormulas(fGroups, "genform", ...)
})
```

The `withOpt` function will temporarily change the given option(s) while executing a given code block and restore it afterwards (it is very similar to the `with_options()` function from the [withr] `R` package). Furthermore, notice how `withOpt()` does not require you to prefix the option names with `patRoon.`. 

#### Multiprocessing with futures

The primary goal of the "future" method is to allow parallel processing on one or more external computers. Since it uses the [future] `R` package, many approaches are supported, such as local parallelization (similar to the `classic` method), cluster computing via multiple networked computers and more advanced HPC approaches such as `slurm` via the [future.batchtools] `R` package. This parallelization method can be activated as follows:

```{r eval=FALSE}
options(patRoon.MP.method = "future")

# set a future plan

# example 1: start a local cluster with four nodes
future::plan("cluster", workers = 4)

# example 2: start a networked cluster with four nodes on PC with hostname "otherpc"
future::plan("cluster", workers = rep("otherpc", 4)) 
```

Please see the documentation of the respective packages (_e.g._ [future] and [future.batchtools]) for more details on how to configure the workers.

The `withOpt()` function introduced in the previous subsection can also be used to temporarily switch between parallelization approaches, for instance:

```{r eval=FALSE}
# default to future parallelization
options(patRoon.MP.method = "future")
future::plan("cluster", workers = 4)

# ... do workflow

# do classic parallelization for GenForm
withOpt(MP.method = "classic", {
    formulas <- generateFormulas(fGroups, "genform", ...)
})

# .. do more workflow
```

#### Logging

Most tools that are executed in parallel will log their output to text files. These files may contain valuable information, for instance, when an error occurred. By default, the logfiles are stored in the `log` directory placed in the current working directory. However, you can change this location by setting the `patRoon.MP.logPath` option. If you set this option to `FALSE` then no logging occurs.

### Notes when using parallelization with futures

Some important notes when using the `future` parallelization method:

* `GenForm` currently performs less optimal with future multiprocessing to the `classic` approach. Nevertheless, it may still be interesting to use the `future` method to move the computations to another system to free up resources on your local system.
* Behind the scenes the [future.apply] package is used to schedule the tools to be executed. The `patRoon.MP.futureSched` option sets the value for the `future.scheduling` argument to the `future_lapply()` function, and therefore allows you to tweak the scheduling.
* Make sure that `patRoon` and in case of multiprocessing that the tool to be executed (`MetFrag`, `SIRIUS` etc.) are exactly the _same_ version on all computing hosts.
* Make sure that `patRoon` is properly configured on all hosts, _e.g._ set the `patRoon.path.XXX` options to ensure all tools can be found.
* For `MetFrag` annotation: if a local database such as `PubChemLite` is used, it must be present on each computing node as well. Furthermore, the local computer (even if not used for the computations) _also_ must have this file present. Like the previous point, make sure that the `patRoon.path.XXX` options are set properly.
* If you encounter errors then it may be handy to switch to `future::plan("sequential")` and see if it works or you get more descriptive error messages.
* In order to restart the nodes, for instance after re-configuring `patRoon`, updating `R` packages etc, simply re-execute `future::plan(...)`.
* Setting the `future.debug` package option to `TRUE` may give you more insight what is happening to find problems.


```{r child=file.path(vignDir, "shared", "_refs.Rmd")}
```
# Sets workflows {#setsWorkflow}

In LC-HRMS screening workflows it is typically desired to be able to detect a broad range of chemicals. For this reason, the samples are often measured twice: with positive and negative ionization. Most data processing steps are only suitable for data with the same polarity, for instance, due to the fact that the _m/z_ values in mass spectra are inherently different (e.g. `[M+H]+` vs `[M-H]-`) and MS/MS fragmentation occurs differently. As a result, the screening workflow has to be done twice, which generally requires more time and complicates comparing and interpretation of the complete (positive and negative) dataset.

In `patRoon` version 2.0 the _sets workflow_ is introduced. This allows you to perform a single non-target screening workflow from different _sets_ of analyses files. Most commonly, each set represents a polarity, hence, there is a positive and negative set. However, more than two sets are supported, and other distinctions between sets are also possible, for instance, samples that were measured with different MS/MS techniques. Another important advantage of the sets workflow is that MS/MS data from different sets can be combined to provide more comprehensive annotations of features. The most important limitation is that (currently) the chromatographic method that was used when analyzing the samples from each set needs to be equal, since retention times are used to group features among the sets.

Performing a sets workflow usually only requires small modifications compared to a 'regular' `patRoon` workflow. This chapter outlines how to perform such workflows and how to use its unique functionality for data processing. It is assumed that the reader is already familiar with performing 'regular' workflows, which were discussed in the previous chapters.

## Initiating a sets workflow

A sets workflow is not much different than a 'regular' (or non-sets) workflow. For instance, consider the following workflow:

```{r initSets-reg,eval=FALSE}
anaInfo <- patRoonData::exampleAnalysisInfo("positive")
fList <- findFeatures(anaInfo, "openms")
fGroups <- groupFeatures(fList, "openms")
fGroups <- filter(fGroups, absMinIntensity = 10000, relMinReplicateAbundance = 1, maxReplicateIntRSD = 0.75,
                  blankThreshold = 5, removeBlanks = TRUE)

mslists <- generateMSPeakLists(fGroups, "mzr")
formulas <- generateFormulas(fGroups, mslists, "genform", adduct = "[M+H]+")
compounds <- generateCompounds(fGroups, mslists, "metfrag", adduct = "[M+H]+")

reportHTML(fGroups, MSPeakLists = mslists, formulas = formulas, compounds = compounds)
```

This example uses the example data from [patRoonData] to obtain a feature group dataset, which is cleaned-up afterwards. Then, feature groups are annotated and all the results are reported.

Converting this to a _sets workflow_:

```{r initSets,eval=FALSE}
anaInfoPos <- patRoonData::exampleAnalysisInfo("positive")
anaInfoNeg <- patRoonData::exampleAnalysisInfo("negative")
fListPos <- findFeatures(anaInfoPos, "openms")
fListNeg <- findFeatures(anaInfoNeg, "openms")
fList <- makeSet(fListPos, fListNeg, adducts = c("[M+H]+", "[M-H]-"))

fGroups <- groupFeatures(fList, "openms")
fGroups <- filter(fGroups, absMinIntensity = 10000, relMinReplicateAbundance = 1, maxReplicateIntRSD = 0.75,
                  blankThreshold = 5, removeBlanks = TRUE)

mslists <- generateMSPeakLists(fGroups, "mzr")
formulas <- generateFormulas(fGroups, mslists, "genform")
compounds <- generateCompounds(fGroups, mslists, "metfrag")

reportHTML(fGroups, MSPeakLists = mslists, formulas = formulas, compounds = compounds)
```

This workflow will do all the steps for positive _and_ negative data.

```{r setsWorkflow,echo=FALSE,out.width="75%"}
plotGV("
digraph Workflow {
  graph [ rankdir = TB, compound = true ]
  node [ shape = box,
         fixedsize = true,
         width = 2.2,
         height = 0.6,
         fontsize = 16,
         fillcolor = darkseagreen1,
         style = filled ]

    'Pre-treatment (+)' -> 'Find features (+)' -> 'makeSet'
    'Pre-treatment (-)' -> 'Find features (-)' -> 'makeSet'
    'makeSet' -> 'Group features' -> 'Annotation, ...'
}", height = 300, width = 250)
```

Only a few modifications were necessary:

* The [analysis information](#anaInfo) is obtained for positive and negative data (i.e. per set)
* Features are found for each set separately.
* `makeSet` is used to combine the feature data
* There is no need to specify the adduct anymore in the annotation steps.

> **_NOTE_** The `analysis` names for the [analysis information](#anaInfo) must be _unique_ for each row, even among sets. Furthermore, replicate groups should not contain analyses from different sets.

The key principle to make sets workflows work is performed by `makeSet`. This method function takes different `features` objects (or `featureGroups`, discussed later) to combine the feature data across sets. During this step features are _neutralized_: the feature _m/z_ data is converted to neutral feature masses. This step ensures that when features are grouped with `groupFeatures`, its algorithms are able to find the same feature among different sets, even when different MS ionization modes were used during acquisition. However, please note that (currently) no additional chromatographic alignment steps between sets are performed. For this reason, the chromatographic methodology that is used to acquire the data must be the same for all sets.

The feature neutralization step relies on adduct data. In the example above, it is simply assumed that all features measured with positive mode are protonated (M+H) species, and all negative features are deprotonated (M-H). It is also possible to use [adduct annotations](#incorpAdductIso) for neutralization; this is discussed later.

> **_NOTE_** The [newProject tool](#newProject) can be used to easily generate a sets workflow. Simply select "both" for the _Ionization_ option.

## Generating sets workflow data

As was shown in the previous section, the generation of workflow data with a sets workflow largely follows that as what was discussed in the previous chapters. The same generator functions are used:

Workflow step         | Function                                          | Output S4 class              
--------------------- | ------------------------------------------------- | ----------------------------
Grouping features     | `groupFeatures()`                                 | `featureGroupsSet`
Suspect screening     | `screenSuspects()`                                | `featureGroupsScreeningSet` 
MS peak lists         | `generateMSPeakLists()`                           | `MSPeakListsSet`
Formula annotation    | `generateFormulas()`                              | `formulasSet`
Compound annotation   | `generateCompounds()`                             | `compoundsSet`
Componentization      | `generateComponents()`                            | algorithm dependent

(the data pre-treatment and feature finding steps have been omitted as they are not specific to sets workflows).

While the same function generics are used to generate data, the class of the output objects differ (e.g. `formulasSet` instead of `formulas`). However, since all these classes _inherit_ from their non-sets workflow counterparts, using the workflow data in a sets workflow is nearly identical to what was discussed in the previous chapters (further discussed in the next section).

As discussed before, an important step is the neutralization of features. Other workflow steps also have internal mechanics to deal with data from different sets:

Workflow step               | Handling of set data
--------------------------- | ------------------------------------------------------
Finding/Grouping features   | Neutralization of _m/z_ values
Suspect screening           | Merging results from screening performed for each set
Componentization            | Algorithm dependent (discussed below)
MS peak lists               | MS data is obtained and stored per set. The final peak lists are combined (_not_ averaged)
Formula/Compound annotation | Annotation is performed for each set separately and used to generate a final consensus

In most cases the algorithms of the workflow steps are first performed for each set, and this data is then merged. To illustrate the importance of this, consider these examples

* A suspect screening with a suspect list that contains known MS/MS fragments
* Annotation where MS/MS fragments are used to predict the chemical formula
* Componentization in order to establish adduct assignments for the features

In all cases data is used that is highly dependent on the MS method (eg polarity) that was used to acquire the sample data. Nevertheless, all the steps needed to obtain and combine set data are performed automatically in the background, and are therefore largely invisible.

> **_NOTE_** Because feature groups in sets workflows always have [adduct annotations](#useAdductAnn), it is never required to specify the adduct or ionization mode when generating annotations, components or do suspect screening (_i.e._ the `adduct`/`ionization` arguments should not be specified).

### Componentization

When the componentization algorithms related to adduct/isotope annotations (e.g. [CAMERA], [RAMClustR] and [cliqueMS]) and [nontarget] are used, then componentization occurs per set and the final object (a `componentsSet` or `componentsNTSet`) contains all the components together. Since these algorithms are highly dependent upon MS data polarity, no attempt is made to merge components from different sets.

The other componentization algorithms work on the complete data. For more details, see the reference manual (`?generateComponents`).

### Formula and compound annotation

For formula and compound annotation, the data generated for each set is combined to generate a _set consensus_. The annotation tables are merged, scores are averaged and candidates are re-ranked. More details can be found in the reference manual (e.g. `?generateCompounds`). In addition, it possible to only keep candidates that exist in a minimum number of sets. For this, the `setThreshold` and `setThresholdAnn` argument can be used:

```{r setThreshold,eval=FALSE}
# candidate must be present in all sets
formulas <- generateFormulas(fGroups, mslists, "genform", setThreshold = 1)
# candidate must be present in all sets with annotation data
compounds <- generateCompounds(fGroups, mslists, "metfrag", setThresholdAnn = 1)
```

In the first example, a formula candidate for a feature group is only kept if it was found for all of the sets. In the second example, a compound candidate is only kept if it was present in all of the sets with annotation data available. The following examples of a common positive/negative sets workflow illustrate the differences:

Candidate | annotations | candidate present | `setThreshold=1` | `setThresholdAnn=1`
--------- | ----------- | ----------------- | ---------------- | ---------------------
\#1       | `+`, `-`    | `+`, `-`          | Keep             | Keep
\#2       | `+`, `-`    | `+`               | Remove           | Remove
\#3       | `+`         | `+`               | Remove           | Keep

For more information refer to the reference manual (e.g. `?generateCompounds`).

## Selecting adducts to improve grouping {#setsAdducts}

The `selectIons()` and `adduct()` functions [discussed before](#incorpAdductIso) can also improve sets workflows. This is because the adduct annotations can be used to improve feature neutralization, which in turn will improve grouping features between positive and negative ionization data. Once adduct annotations are set the features will be re-neutralized and re-grouped.

A typical workflow with  `selectIons` looks like this:

```{r setsSelectIonsWorkflow,echo=FALSE,out.width="75%"}
plotGV("
digraph Workflow {
  graph [ rankdir = TB, compound = true ]
  node [ shape = box,
         fixedsize = true,
         width = 2.2,
         height = 0.6,
         fontsize = 16,
         fillcolor = darkseagreen1,
         style = filled ]

    'Pre-treatment (+)' -> 'Find features (+)' -> 'makeSet'
    'Pre-treatment (-)' -> 'Find features (-)' -> 'makeSet'
    'makeSet' -> 'Group features' -> Componentization -> selectIons -> 'Annotation, ...'
}", height = 350, width = 250)
```

```{r setsSelectIons,eval=FALSE}
# as before ...
anaInfoPos <- patRoonData::exampleAnalysisInfo("positive")
anaInfoNeg <- patRoonData::exampleAnalysisInfo("negative")
fListPos <- findFeatures(anaInfoPos, "openms")
fListNeg <- findFeatures(anaInfoNeg, "openms")

fGroupsPos <- groupFeatures(fListPos, "openms")
fGroupsNeg <- groupFeatures(fListNeg, "openms")
fList <- makeSet(fListPos, fListNeg, adducts = c("[M+H]+", "[M-H]-"))

fGroups <- groupFeatures(fList, "openms")
fGroups <- filter(fGroups, absMinIntensity = 10000, relMinReplicateAbundance = 1, maxReplicateIntRSD = 0.75,
                  blankThreshold = 5, removeBlanks = TRUE)

components <- generateComponents(fGroups, "openms")
fGroups <- selectIons(fGroups, components, c("[M+H]+", "[M-H]-"))

# do rest of the workflow...
```

The first part of the workflow is exactly the same as was introduced in the beginning of this chapter. Furthermore, note that for sets workflows, `selectIons` needs a preferential adduct for each set.

The `adducts` function can also be used to obtain and modify adduct annotations. For sets workflows, these functions operate _per set_:

```{r setsAdduct,eval=FALSE}
adducts(fGroups, set = "positive")[1:5]
adducts(fGroups, set = "positive")[4] <- "[M+K]+"
```

If you want to modify annotations for multiple sets, it is best to delay the re-gouping step:

```{r setsAdductNoRegroup,eval=FALSE}
adducts(fGroups, set = "positive", reGroup = FALSE)[4] <- "[M+K]+"
adducts(fGroups, set = "negative", reGroup = TRUE)[10] <- "[M-H2O]-"
```

Setting `reGroup=FALSE` will not perform any re-neutralization and re-grouping, which preserves feature group names and safes processing time. However, it is **crucial** that the re-grouping step is eventually performed at the end.

## Processing data

All data objects that are generated during a sets workflow _inherit_ from the classes from a 'regular' workflow. This means that, with some minor exceptions, _all_ of the data processing functionality discussed in the [previous chapter](#processing) (e.g. subsetting, inspection, filtering, plotting, reporting) is also applicable to a sets workflow. For instance, the `as.data.table()` method can be used for general inspection:

```{r include=FALSE,eval=runData}
fGroupsO <- fGroups
fGroups <- fGroupsSets
mslistsO <- mslists
mslists <- mslistsSets
compoundsO <- compounds
compounds <- compoundsSets
```

```{r setsProcData,eval=runData}
as.data.table(compounds)[1:5, c("group", "score", "compoundName", "set")]
```

In addition, some the data processing functionality contains additional functionality for a sets workflow:

```{r setsProcDataExtra,eval=FALSE}

# only keep feature groups that have positive data
fGroupsPos <- fGroups[, sets = "positive"]
# only keep feature groups that have feature data for all sets
fGroupsF <- filter(fGroups, relMinSets = 1)

# only keep feature groups with features present in both polarities
fGroupsPosNeg <- overlap(fGroups, which = c("positive", "negative"), sets = TRUE)
# only keep feature groups with features that are present only in positive mode
fGroupsOnlyPos <- unique(fGroups, which = "positive", sets = TRUE)
```

And plotting:

```{r setsProcDataExtraPlotting,eval=runData,fig.show="hold"}
plotVenn(fGroups, sets = TRUE) # compare positive/negative features
plotSpectrum(compounds, index = 1, groupName = "M198_R317_272", MSPeakLists = mslists,
             plotStruct = TRUE)
```

```{r include=FALSE,eval=runData}
fGroups <- fGroupsO
mslists <- mslistsO
compounds <- compoundsO
```

The reference manual for the workflow objects contains specific notes applicable to sets workflows (`?featureGroups`, `?compounds` etc).

## Advanced

### Initiating a sets workflow from feature groups

The `makeSet` function can also be used to initiate a sets workflow from feature groups:

```{r initSetFG,eval=FALSE}
# as before ...
anaInfoPos <- patRoonData::exampleAnalysisInfo("positive")
anaInfoNeg <- patRoonData::exampleAnalysisInfo("negative")
fListPos <- findFeatures(anaInfoPos, "openms")
fListNeg <- findFeatures(anaInfoNeg, "openms")

fGroupsPos <- groupFeatures(fListPos, "openms")
fGroupsNeg <- groupFeatures(fListNeg, "openms")

fGroups <- makeSet(fGroupsPos, fGroupsNeg, groupAlgo = "openms",
                   adducts = c("[M+H]+", "[M-H]-"))

# do rest of the workflow...

```

In this case `makeSet` combines the positive and negative (un-grouped) features, neutralizes them and re-groups them all together (with the algorithm specified by `groupAlgo`).

While this option involves some extra steps, an advantage is that allows processing the feature data before they are combined, e.g.:

```{r initSetFGFilt,eval=FALSE}
fGroupsPos <- groupFeatures(fListPos, "openms")
fGroupsNeg <- groupFeatures(fListNeg, "openms")

# apply intensity theshold filters. Lower threshold for negative.
fGroupsPos <- filter(fGroupsPos, absMinIntensity = 1E4)
fGroupsNeg <- filter(fGroupsNeg, absMinIntensity = 1E3)

fGroups <- makeSet(fGroupsPos, fGroupsNeg, groupAlgo = "openms",
                   adducts = c("[M+H]+", "[M-H]-"))

```

Visually, this workflow looks like this:

```{r setsWorkflowG2,echo=FALSE,out.width="75%"}
plotGV("
digraph Workflow {
  graph [ rankdir = TB, compound = true ]
  node [ shape = box,
         fixedsize = true,
         width = 2.2,
         height = 0.6,
         fontsize = 16,
         fillcolor = darkseagreen1,
         style = filled ]

    'Find features (+)' -> 'Group features (+)' -> 'filter (+)' -> 'makeSet'
    'Find features (-)' -> 'Group features (-)' -> 'filter (-)' -> 'makeSet'
    'makeSet' -> '...'
}", height = 300, width = 250)
```

Of course, any other processing steps on the feature groups data such as subsetting and [visually checking features](#intReview) are also possible before the sets workflow is initiated. Furthermore, it is also possible to perform [adduct annotations](#incorpAdductIso) prior to grouping, which is an alternative way to improve neutralization to what [was discussed before](#setsAdducts). 

### Inspecting and converting set objects

The following generic functions may be used to inspect or convert data from sets workflows:

Generic      | Purpose                               | Notes
------------ | ------------------------------------- | ---------------------------------------
`sets`       | Return the names of the sets in this object.
`setObjects` | Obtain the raw data objects that were used to construct this object. | Not available for features and feature groups.
`unset`      | Converts this object to a regular workflow object. | The `set` argument must be given to specify which of the set data is to be converted. This function will restore the original _m/z_ values of features.

These methods are heavily used internally, but rarely needed otherwise. More details can be found in the reference manual.
```{r setup, include = FALSE}

# otherwise Linux will get into memory troubles...
knitr::knit_meta("latex_dependency", clean = TRUE)

knitr::opts_chunk$set(
    fig.width = 6, fig.height = 4, out.width = "50%"
)

source(file.path(vignDir, "shared", "init.R"))

runData <- T
doOpt <- runData
if (runData)
{
    # try to sync with tutorial so cache can be re-used
    anaInfo <- patRoonData::exampleAnalysisInfo()
    anaInfoRG <- anaInfo
    anaInfoRG$group <- c(rep("repl1", 2),
                         rep("repl2", 2),
                         rep("repl3", 1),
                         rep("repl4", 1))

    # set max proc to 1 to limit FFM memory usage a bit on CI
    getFeats <- function(ai) withr::with_options(list(patRoon.multiproc.max = 1), findFeatures(ai, "openms", noiseThrInt = 2E4))
    
    fList <- getFeats(anaInfo)
    fGroups <- fGroupsUF <- groupFeatures(fList, "openms")
    fGroups <- filter(fGroups, preAbsMinIntensity = 100, absMinIntensity = 10000,
                      relMinReplicateAbundance = 1, maxReplicateIntRSD = 0.75,
                      blankThreshold = 5, removeBlanks = TRUE,
                      retentionRange = c(120, Inf), mzRange = NULL)

    fListRG <- getFeats(anaInfoRG)
    fGroupsRG <- groupFeatures(fListRG, "openms")
    fGroupsRG <- filter(fGroupsRG, preAbsMinIntensity = 100, absMinIntensity = 10000)
    
    anaInfoConc <- generateAnalysisInfo(paths = patRoonData::exampleDataPath(),
                                        groups = c(rep("solvent", 3), rep("standard", 3)),
                                        blanks = "solvent",
                                        concs = c(NA, NA, NA, 1, 2, 3))
    fListConc <- getFeats(anaInfoConc)
    fGroupsConc <- groupFeatures(fListConc, "openms")
    fGroupsConc <- filter(fGroupsConc, preAbsMinIntensity = 100, absMinIntensity = 10000,
                          relMinReplicateAbundance = 1, maxReplicateIntRSD = 0.75,
                          blankThreshold = 5, removeBlanks = TRUE,
                          retentionRange = c(120, Inf), mzRange = NULL)
    
    fGroupsAnn <- screenSuspects(fGroups, patRoonData::suspectsPos[patRoonData::suspectsPos$name %in%
                                                                       c("1H-benzotriazole", "N-Phenyl urea",
                                                                         "2-Hydroxyquinoline", "DEET"), ],
                                 onlyHits = TRUE)
    avgPListParams <- getDefAvgPListParams(clusterMzWindow = 0.002)
    mslists <- generateMSPeakLists(fGroupsAnn, "mzr", maxMSRtWindow = 5, precursorMzWindow = 4,
                                   avgFeatParams = avgPListParams, avgFGroupParams = avgPListParams)
    mslists <- filter(mslists, relMSMSIntThr = 0.02, topMSMSPeaks = 10)
    formulas <- formsGF <- generateFormulas(fGroupsAnn, mslists, "genform", relMzDev = 5,
                                            adduct = "[M+H]+", elements = "CHNOPSCl",
                                            calculateFeatures = TRUE, featThresholdAnn = 0.75)
    formsSIR <- generateFormulas(fGroupsAnn, mslists, "sirius", elements = "CHNOPSCl", adduct = "[M+H]+",
                                 calculateFeatures = FALSE)

    compsMF <- compounds <-
        generateCompounds(fGroupsAnn, mslists, "metfrag", method = "CL",
                          dbRelMzDev = 5, fragRelMzDev = 5, fragAbsMzDev = 0.002,
                          adduct = "[M+H]+", database = "pubchemlite", maxCandidatesToStop = 5000)
    
    componCAM <- components <- generateComponents(fGroups, "camera", ionization = "positive")
    componInt <- generateComponents(fGroupsRG, "intclust")
    componNT <- generateComponents(fGroupsUF[, 1:200], "nontarget", minlength = 3, ionization = "positive")
    
    compsClust <- makeHCluster(compsMF)
    
    if (doOpt)
    {
        pSet <- generateFeatureOptPSet("openms")
        ftOpt <- optimizeFeatureFinding(anaInfo[1, ], "openms", pSet, maxIterations = 2,
                                        paramRanges = list(noiseThrInt = c(1500, Inf)))
    }
    
    fGroupsSets <- groupFeatures(makeSet(fList, getFeats(patRoonData::exampleAnalysisInfo("negative")),
                                         adducts = c("[M+H]+", "[M-H]-")), "openms")
    fGroupsSetsAnn <- screenSuspects(fGroupsSets, patRoonData::suspectsPos[patRoonData::suspectsPos$name
                                                                           == "Monuron", -2],
                                     onlyHits = TRUE)
    mslistsSets <- generateMSPeakLists(fGroupsSetsAnn, "mzr", maxMSRtWindow = 5, precursorMzWindow = 4,
                                       avgFeatParams = avgPListParams, avgFGroupParams = avgPListParams)
    mslistsSets <- filter(mslistsSets, topMSMSPeaks = 10)
    compoundsSets <- generateCompounds(fGroupsSetsAnn, mslistsSets, "metfrag", database = "pubchemlite")
    
    saveRDS(list(fList = fList, fGroups = fGroups, fGroupsUF = fGroupsUF, fListRG = fListRG,
                 fGroupsRG = fGroupsRG, mslists = mslists, formulas = formulas, compsMF = compsMF,
                 componCAM = componCAM, componInt = componInt, componNT = componNT,
                 compsClust = compsClust,
                 ftOpt = if (doOpt) ftOpt else NULL, fGroupsSets = fGroupsSets, fGroupsSetsAnn = fGroupsSetsAnn,
                 mslistsSets = mslistsSets, compoundsSets = compoundsSets),
            "~/handbook-obj.Rds")
}
```

```{css code=readLines(file.path(vignDir, "styles.css")),echo=FALSE,eval=knitr::is_html_output()}
```
# Processing workflow data {#processing}

The previous chapter mainly discussed how to create workflow data. This chapter will discuss how to _use_ the data.

## Inspecting results

Several generic functions exist that can be used to inspect data that is stored in a particular object (e.g. features, compounds etc):

Generic                                             | Classes                       | Remarks
--------------------------------------------------- | ----------------------------- | ------------------------------------------
`length()`                                          | All                           | Returns the length of the object (e.g. number of features, compounds etc)
`algorithm()`                                       | All                           | Returns the name of the algorithm used to generate the object. 
`groupNames()`                                      | All                           | Returns all the unique identitifiers (or names) of the feature groups for which this object contains results.
`names()`                                           | `featureGroups`, `components` | Returns names of the feature groups (similar to `groupNames()`) or components
`show()`                                            | All                           | Prints general information.
`"[["` / `"$"` operators                            | All                           | Extract general information, see below.
`as.data.table()` / `as.data.frame()`               | All                           | Convert data to a `data.table` or `data.frame`, see below.
`analysisInfo()`, `analyses()`, `replicateGroups()` | `features`, `featureGroups`   | Returns the [analysis information](#anaInfo), analyses or replicate groups for which this object contains data.
`groupInfo()`                                       | `featureGroups`               | Returns feature group information (_m/z_ and retention time values).
`screenInfo()`                                      | `featureGroupsScreening`      | Returns information on hits from [suspect screening](#suspscr).
`componentInfo()`                                   | `components`                  | Returns information for all components.
`annotatedPeakList()`                               | `formulas`, `compounds`       | Returns a table with annotated mass peaks (see below).

The common `R` extraction operators `"[["`, `"$"` can be used to obtain data for a particular feature groups, analysis etc:

```{r extrOp,eval=runData}
# Feature table (only first columns for readability)
fList[["standard-1"]][, 1:6]

# Feature group intensities
fGroups$M120_R268_30
fGroups[[1, "M120_R268_30"]] # only first analysis

# obtains MS/MS peak list  (feature group averaged data)
mslists[["M120_R268_30"]]$MSMS

# get all formula candidates for a feature group
formulas[["M120_R268_30"]][, 1:7]

# get all compound candidates for a feature group
compounds[["M120_R268_30"]][, 1:4]

# get a table with information of a component
components[["CMP7"]][, 1:6]
```

A more sophisticated way to obtain data from a workflow object is to use `as.data.table()` or `as.data.frame()`. These functions will convert _all_ information within the object to a table (`data.table` or `data.frame`) and allow various options to add extra information. An advantage is that this common data format can be used with many other functions within `R`. The output is in a [tidy format](https://r4ds.had.co.nz/tidy-data.html).

> **_NOTE_** If you are not familiar with `data.table` and want to know more see [data.table]. Briefly, this is a more efficient and largely compatible alternative to the regular `data.frame`.

> **_NOTE_** The `as.data.frame()` methods defined in `patRoon` simply convert the results from `as.data.table()`, hence, both functions are equal in their usage and are defined for the same object classes.

Some typical examples are shown below.

```{r asDT,eval=runData}
# obtain table with all features (only first columns for readability)
as.data.table(fList)[, 1:6]

# Returns group info and intensity values for each feature group
as.data.table(fGroups, average = TRUE) # average intensities for replicates

# As above, but with extra suspect screening information
# (select some columns to simplify the output below)
as.data.table(fGroupsSusp, average = TRUE, collapseSuspects = NULL,
              onlyHits = TRUE)[1:3, c("group", "susp_name", "susp_compRank", "susp_annSimBoth", "susp_estIDLevel")]

# Returns all peak lists for each feature group
as.data.table(mslists)

# Returns all formula candidates for each feature group with scoring
# information, neutral loss etc
as.data.table(formulas)[, 1:6]

# Returns all compound candidates for each feature group with scoring and other metadata
as.data.table(compounds)[, 1:4]

# Returns table with all components (including feature group info, annotations etc)
as.data.table(components)[, 1:6]
```

Finally, the `annotatedPeakList()` function is useful to inspect annotation results for a formula or compound candidate:

```{r annPList,eval=runData}
# formula annotations for the first formula candidate of feature group M137_R249_53
annotatedPeakList(formulas, index = 1, groupName = "M137_R249_53",
                  MSPeakLists = mslists)

# compound annotation for first candidate of feature group M137_R249_53
annotatedPeakList(compounds, index = 1, groupName = "M137_R249_53",
                  MSPeakLists = mslists)
```


More advanced examples for these functions are shown below.


```{r inspEx,eval=FALSE}
# Feature table, can also be accessed by numeric index
fList[[1]]
mslists[["standard-1", "M120_R268_30"]] # feature data (instead of feature group averaged)
formulas[[1, "M120_R268_30"]] # feature data (if available, i.e. calculateFeatures=TRUE)
components[["CMP1", 1]] # only for first feature group in component

as.data.frame(fList) # classic data.frame format, works for all objects
as.data.table(fGroups) # return non-averaged intensities (default)
as.data.table(fGroups, features = TRUE) # include feature information
as.data.table(mslists, averaged = FALSE) # peak lists for each feature
as.data.table(mslists, fGroups = fGroups) # add feature group information

as.data.table(formulas, countElements = c("C", "H")) # include C/H counts (e.g. for van Krevelen plots)
# add various information for organic matter characterization (common elemental
# counts/ratios, classifications etc)
as.data.table(formulas, OM = TRUE)

as.data.table(compounds, fGroups = fGroups) # add feature group information
as.data.table(compounds, fragments = TRUE) # include information of all annotated fragments

annotatedPeakList(formulas, index = 1, groupName = "M120_R268_30",
                  MSPeakLists = mslists, onlyAnnotated = TRUE) # only include annotated peaks
annotatedPeakList(compounds, index = 1, groupName = "M120_R268_30",
                  MSPeakLists = mslists, formulas = formulas) # include formula annotations
```

## Filtering {#filtering}

During a non-target workflow it is not uncommon that some kind of data-cleanup is necessary. Datasets are often highly complex, which makes separating data of interest from the rest highly important. Furthermore, general cleanup typically improves the quality of the dataset, for instance by removing low scoring annotation results or features that are unlikely to be 'correct' (e.g. noise or present in blanks). For this reason `patRoon` supports _many_ different filters that easily clean data produced during the workflow in a highly customizable way.

All major workflow objects (e.g. `featureGroups`, `compounds`, `components` etc.) support filtering operations by the `filter()` generic. This function takes the object to be filtered as first argument and any remaining arguments describe the desired filter options. The `filter()` generic function then returns the modified object back. Some examples are shown below.

```{r filtGen,eval=FALSE}
# remove low intensity (<500) features
features <- filter(features, absMinIntensity = 500)

# remove features with intensities lower than 5 times the blank
fGroups <- filter(fGroups, blankThreshold = 5)

# only retain compounds with >1 explained MS/MS peaks
compounds <- filter(compounds, minExplainedPeaks = 1)
```

The following sections will provide a more detailed overview of available data filters.

> **_NOTE_**  Some other `R` packages (notably `dplyr`) also provide a `filter()` generic function. To use the `filter()` function from different packages you may need to explicitly specify which one to use in your script. This can be done by prefixing it with the package name, e.g. `patRoon::filter(...)`, `dplyr::filter(...)` etc.


### Features

There are many filters available for feature data:

Filter                                     | Classes                     | Remarks
------------------------------------------ | --------------------------- | ---------------------------------------------------------
`absMinIntensity`, `relMinIntensity`       | `features`, `featureGroups` | Minimum intensity
`preAbsMinIntensity`, `preRelMinIntensity` | `featureGroups`             | Minimum intensity prior to other filtering (see below)
`retentionRange`, `mzRange`, `mzDefectRange`, `chromWidthRange` | `features`, `featureGroups` | Filter by feature properties
`absMinAnalyses`, `relMinAnalyses`         | `featureGroups`             | Minimum feature abundance in all analyses
`absMinReplicates`, `relMinReplicates`     | `featureGroups`             | Minimum feature abundance in different replicates
`absMinFeatures`, `relMinFeatures`         | `featureGroups`             | Only keep analyses with at least this amount of features 
`absMinReplicateAbundance`, `relMinReplicateAbundance` | `featureGroups` | Minimum feature abundance in a replicate group
`maxReplicateIntRSD`                       | `featureGroups`             | Maximum relative standard deviation of feature intensities in a replicate group.
`blankThreshold`                           | `featureGroups`             | Minimum intensity factor above blank intensity
`rGroups`                                  | `featureGroups`             | Only keep (features of) these replicate groups
`results`                                  | `featureGroups`             | Only keep feature groups with formula/compound annotations or componentization results

Application of filters to feature data is important for (environmental) non-target analysis. Especially blank and replicate filters (i.e. `blankThreshold` and `absMinReplicateAbundance`/`relMinReplicateAbundance`) are important filters and are highly recommended to always apply for cleaning up your dataset.

All filters are available for feature group data, whereas only a subset is available for feature objects. The main reason is that other filters need grouping of features between analyses. Regardless, in `patRoon` filtering feature data is less important, and typically only needed when the number of features are extremely large and direct grouping is undesired.

From the table above you can notice that many filters concern both _absolute_ and _relative_ data (i.e. as prefixed with `abs` and `rel`). When a relative filter is used the value is scaled between _0_ and _1_. For instance:

```{r filtFeatRel,eval=FALSE}
# remove features not present in at least half of the analyses within a replicate group
fGroups <- filter(fGroups, relMinReplicateAbundance = 0.5)
```

An advantage of relative filters is that you will not have to worry about the data size involved. For instance, in the above example the filter always takes half of the number of analyses within a replicate group, even when replicate groups have different number of analyses.

Note that multiple filters can be specified at once. Especially for feature group data the order of filtering may impact the final results, this is explained further in the reference manual (i.e. ``?`feature-filtering` ``).

Some examples are shown below.

```{r filtFeat,eval=FALSE}
# filter features prior to grouping: remove any features eluting before first 2 minutes
fList <- filter(fList, retentionRange = c(120, Inf))

# common filters for feature groups
fGroups <- filter(fGroups,
                  absMinIntensity = 500, # remove features <500 intensity
                  relMinReplicateAbundance = 1, # features should be in all analysis of replicate groups
                  maxReplicateIntRSD = 0.75, # remove features with intensity RSD in replicates >75%
                  blankThreshold = 5, # remove features <5x intensity of (average) blank intensity
                  removeBlanks = TRUE) # remove blank analyses from object afterwards

# filter by feature properties
fGroups <- filter(mzDefectRange = c(0.8, 0.9),
                  chromWidthRange = c(6, 120))

# remove features not present in at least 3 analyses
fGroups <- filter(fGroups, absMinAnalyses = 3)

# remove features not present in at least 20% of all replicate groups
fGroups <- filter(fGroups, relMinReplicates = 0.2)

# only keep data present in replicate groups "repl1" and "repl2"
# all other features and analyses will be removed
fGroups <- filter(fGroups, rGroups = c("repl1", "repl2"))

# only keep feature groups with compound annotations
fGroups <- filter(fGroups, results = compounds)
# only keep feature groups with formula or compound annotations
fGroups <- filter(fGroups, results = list(formulas, compounds))
```

### Suspect screening

Several additional filters are available for feature groups obtained with `screenSuspects()`:


Filter                                    | Classes                  | Remarks
----------------------------------------- | ------------------------ | ---------------------------------------------------------
`onlyHits`                                | `featureGroupsScreening` | Only retain feature groups assigned to one or more suspects.
`selectHitsBy`                            | `featureGroupsScreening` | Select the feature group that matches best with a suspect (in case there are multiple).
`selectBestFGroups`                       | `featureGroupsScreening` | Select the suspect that matches best with a feature group (in case there are multiple).
`maxLevel`, `maxFormRank`, `maxCompRank`  | `featureGroupsScreening` | Only retain suspect hits with identification/annotation ranks below a threshold.
`minAnnSimForm`, `minAnnSimComp`, `minAnnSimBoth` | `featureGroupsScreening` | Remove suspect hits with annotation similarity scores below this value.
`absMinFragMatches`, `relMinFragMatches`  | `featureGroupsScreening` | Only keep suspect hits with a minimum (relative) number of fragment matches from the suspect list.

> **_NOTE_**: most filters only remove suspect hit results. Set `onlyHits=TRUE` to also remove any feature groups that end up without suspect hits.

The `selectHitsBy` and `selectBestFGroups` filters are useful to remove duplicate hits (one suspect assigned to multiple feature groups or multiple feature groups assigned to the same suspect, respectively). The former selects based on either best identification level (`selectHitsBy="level"`) or highest mean intensity (`selectHitsBy="intensity"`). The `selectBestFGroups` can only be `TRUE`/`FALSE` and always selects by best identification level.

Some examples are shown below.

```{r filtSusp,eval=FALSE}
# only keep feature groups assigned to at least one suspect
fGroupsSusp <- filter(fGroupsSusp, onlyHits = TRUE)
# remove duplicate suspect to feature group matches and keep the best
fGroupsSusp <- filter(fGroupsSusp, selectHitsBy = "level")
# remove suspect hits with ID levels >3 and make sure no feature groups
# are present without suspect hits afterwards
fGroupsSusp <- filter(fGroupsSusp, maxLevel = 3, onlyHits = TRUE)
```

### Annotation

There are various filters available for handling annotation data:

Filter                                        | Classes                        | Remarks
--------------------------------------------- | ------------------------------ | -------------------------------------------------
`absMSIntThr`, `absMSMSIntThr`, `relMSIntThr`, `relMSMSIntThr` | `MSPeakLists` | Minimum intensity of mass peaks
`topMSPeaks`, `topMSMSPeaks`                  | `MSPeakLists`                  | Only keep most intense mass peaks
`withMSMS`                                    | `MSPeakLists`                  | Only keep results with MS/MS data
`minMSMSPeaks`                                | `MSPeakLists`                  | Only keep an MS/MS peak list if it contains a minimum number of peaks (excluding the precursor peak)
`annotatedBy`                                 | `MSPeakLists`                  | Only keep MS/MS peaks that have formula or compound annotations
`minExplainedPeaks`                           | `formulas`, `compounds`        | Minimum number of annotated mass peaks
`elements`, `fragElements`, `lossElements`    | `formulas`, `compounds`        | Restrain elemental composition
`topMost`                                     | `formulas`, `compounds`        | Only keep highest ranked candidates 
`minScore`, `minFragScore`, `minFormulaScore` | `compounds`                    | Minimum compound scorings
`scoreLimits`                                 | `formulas`, `compounds`        | Minimum/Maximum scorings
`OM`                                          | `formulas`, `compounds`        | Only keep candidates with likely elemental composition found in organic matter

Several intensity related filters are available to clean-up MS peak list data. For instance, the `topMSPeaks`/`topMSMSPeaks` filters provide a simple way to remove noisy data by only retaining a defined number of most intense mass peaks. Note that none of these filters will remove the precursor mass peak of the feature itself.

The filters applicable to formula and compound annotation generally concern minimal scoring or chemical properties. The former is useful to remove unlikely candidates, whereas the second is useful to focus on certain study specific chemical properties (e.g. known neutral losses).

Common examples are shown below.

```{r filtAnn,eval=FALSE}
# intensity filtering
mslists <- filter(mslists,
                  absMSIntThr = 500, # minimum MS mass peak intensity of 500
                  relMSMSIntThr = 0.1) # minimum MS/MS mass peak intensity of 10%

# only retain 10 most intens mass peaks
# (feature mass is always retained)
mslists <- filter(mslists, topMSPeaks = 10)

# remove MS/MS peaks without compound annotations
mslists <- filter(mslists, annotatedBy = compounds)

# remove MS/MS peaks not annotated by either a formula or compound candidate
mslists <- filter(mslists, annotatedBy = list(formulas, compounds))

# only keep formulae with 1-10 sulphur or phosphorus elements
formulas <- filter(formulas, elements = c("S1-10", "P1-10"))

# only keep candidates with MS/MS fragments that contain 1-10 carbons and 0-2 oxygens
formulas <- filter(formulas, fragElements = "C1-10O0-2")

# only keep candidates with CO2 neutral loss
formulas <- filter(formulas, lossElements = "CO2")

# only keep the 15 highest ranked candidates with at least 1 annotated MS/MS peak
compounds <- filter(compounds, minExplainedPeaks = 1, topMost = 15)

# minimum in-silico score
compounds <- filter(compounds, minFragScore = 10)

# candidate should be referenced in at least 1 patent
# (only works if database lists number of patents, e.g. PubChem)
compounds <- filter(compounds,
                    scoreLimits = list(numberPatents = c(1, Inf))
```

> **_NOTE_** As of `patRoon 2.0` MS peak lists are **not** re-generated after a filtering operation (unless the `reAverage` parameter is explicity set to `TRUE`). The reason for this change is that re-averaging invalidates any formula/compound annotation data (e.g. used for plotting and reporting) that were generated prior to the filter operation.


### Components

Finally several filters are available for components:

Filter                       | Remarks
---------------------------- | -----------------
`size`                       | Minimum component size
`adducts`, `isotopes`        | Filter features by adduct/istopes annotation
`rtIncrement`, `mzIncrement` | Filter homologs by retention/mz increment range

Note that these filters are only applied if the components contain the data the filter works on. For instance, filtering by adducts will _not_ affect components obtained from homologous series.

As before, some typical examples are shown below.

```{r filtComp,eval=FALSE}
# only keep components with at least 4 features
componInt <- filter(componInt, minSize = 4)

# remove all features from components are not annotated as an adduct
componRC <- filter(componRC, adducts = TRUE)

# only keep protonated and sodium adducts
componRC <- filter(componRC, adducts = c("[M+H]+", "[M+Na]+"))

# remove all features not recognized as isotopes
componRC <- filter(componRC, isotopes = FALSE)

# only keep monoisotopic mass
componRC <- filter(componRC, isotopes = 0)

# min/max rt/mz increments for homologs
componNT <- filter(componNT, rtIncrement = c(10, 30),
                   mzIncrement = c(16, 50))
```

> **_NOTE_** As mentioned before, components are still in a relative young development phase and results should always be verified!

### Negation

All filters support _negation_: if enabled all specified filters will be executed in an opposite manner. Negation may not be so commonly used, but allows greater flexibility which is sometimes needed for advanced filtering steps. Furthermore, it is also useful to specifically isolate the data that otherwise would have been removed. Some examples are shown below.

```{r filtNeg,eval=FALSE}
# keep all features/analyses _not_ present from replicate groups "repl1" and "repl2"
fGroups <- filter(fGroups, rGroups = c("repl1", "repl2"), negate = TRUE)

# only retain features with a mass defect outside 0.8-0.9
fGroups <- filter(mzDefectRange = c(0.8, 0.9), negate = TRUE)

# remove duplicate suspect hits and only keep the _worst_ hit
fGroupsSusp <- filter(fGroupsSusp, selectHitsBy = "level", negate = TRUE)

# remove candidates with CO2 neutral loss
formulas <- filter(formulas, lossElements = "CO2", negate = TRUE)

# select 15 worst ranked candidates
compounds <- filter(compounds, topMost = 15, negate = TRUE)

# only keep components with <5 features
componInt <- filter(componInt, minSize = 5, negate = TRUE)
```

## Subsetting {#subset}

The previous section discussed the `filter()` generic function to perform various data cleaning operations. A more generic way to select data is by _subsetting_: here you can manually specify which parts of an object should be retained. Subsetting is supported for all workflow objects and is performed by the R subset operator (`"["`). This operator either subsets by one or two arguments, which are referred to as the `i` and `j` arguments.

Class           | Argument `i`   | Argument `j`   | Remarks
--------------- | -------------- | -------------- | ------------------------------------------------
`features`      | analyses       |                |
`featureGroups` | analyses       | feature groups |
`MSPeakLists`   | analyses       | feature groups | peak lists for feature groups will be re-averaged when subset on analyses (by default)
`formulas`      | feature groups |                | 
`compounds`     | feature groups |                |
`components`    | components     | feature groups |

For objects that support two-dimensional subsetting (e.g. `featureGroups`, `MSPeakLists`), either the `i` or `j` argument is optional. Furthermore, unlike subsetting a `data.frame`, the position of `i` and `j` does not change when only one argument is specified:

```{r subsetArgs,eval=FALSE}
df[1, 1] # subset data.frame by first row/column
df[1] # subset by first column
df[1, ] # subset by first row

fGroups[1, 1] # subset by first analysis/feature group
fGroups[, 1] # subset by first feature group (i.e. column)
fGroups[1] # subset by first analysis (i.e. row)
```

The subset operator allows three types of input:

* A logical vector: elements are selected if corresponding values are `TRUE`.
* A numeric vector: select elements by numeric index.
* A character vector: select elements by their name.

When a logical vector is used as input it will be re-cycled if necessary. For instance, the following will select by the first, third, fifth, etc. analysis.

```{r subsetCyc,eval=FALSE}
fGroups[c(TRUE, FALSE)]
```

In order to select by a `character` you will need to know the names for each element. These can, for instance, be obtained by the `groupNames()` (feature group names), `analyses()` (analysis names) and `names()` (names for components or feature groups for `featureGroups` objects) generic functions.

Some more examples of common subsetting operations are shown below.

```{r subsetting,eval=FALSE}
# select first three analyses
fList[1:3]

# select first three analyses and first 500 feature groups
fGroups[1:3, 1:500]

# select all feature groups from first component
fGroupsNT <- fGroups[, componNT[[1]]$group]

# only keep feature groups with formula annotation results
fGroupsForms <- fGroups[, groupNames(formulas)]

# only keep feature groups with either formula or compound annotation results
fGroupsAnn <- fGroups[, union(groupNames(formulas), groupNames(compounds))]

# select first 15 components
components[1:15]

# select by name
components[c("CMP1", "CMP5")]

# only retain feature groups in components for which compound annotations are
# available
components[, groupNames(compounds)]
```

In addition, feature groups can also be subset by given replicate groups or annotation/componentization results (similar to `filter()`). Similarly, suspect screening results can also be subset by given suspect names.

```{r subsetRGSusp,eval=FALSE}
# equal as filter(fGroups, rGroups = ...)
fGroups[rGroups = c("repl1", "repl2")]
# equal as filter(fGroups, results = ...)
fGroups[results = compounds]
# only keep feature groups assigned to given suspects
fGroupsSusp[suspects = c("1H-benzotriazole", "2-Hydroxyquinoline")]
```

> **_NOTE_** As of `patRoon 2.0` MS peak lists are **not** re-generated after a subsetting operation (unless the `reAverage` parameter is explicity set to `TRUE`). The reason for this change is that re-averaging invalidates any formula/compound annotation data (e.g. used for plotting and reporting) that were generated prior to the subset operation.


### Prioritization workflow

An important use case of subsetting is prioritization of data. For instance, after statistical analysis only certain feature groups are deemed relevant for the rest of the workflow. A common prioritization workflow is illustrated below:

```{r prioriWorkflow,echo=FALSE,out.width="75%"}
plotGV("
digraph Prioritization {
  graph [ rankdir = LR ]
  node [ shape = box,
         fixedsize = true,
         width = 2.3,
         height = 1,
         fontsize = 18,
         fillcolor = darkseagreen1,
         style = filled ]

    'Object conversion' -> 'Prioritization' -> Subsetting
}", height = 120, width = 500)
```

During the first step the workflow object is converted to a suitable format, most often using the `as.data.frame()` function. The converted data is then used as input for the prioritization strategy. Finally, these results are then used to select the data of interest in the original object.

A very simplified example of such a process is shown below.

```{r prioriEx,eval=FALSE}
featTab <- as.data.frame(fGroups, average = TRUE)

# prioritization: sort by (averaged) intensity of the "sample" replicate group
# (from high to low) and then obtain the feature group identifiers of the top 5.
featTab <- featTab[order(featTab$standard, decreasing = TRUE), ]
groupsOfInterest <- featTab$group[1:5]

# subset the original data
fGroups <- fGroups[, groupsOfInterest]

# fGroups now only contains the feature groups for which intensity values in the
# "sample" replicate group were in the top 5
```

## Deleting data

The `delete()` generic function can be used to manually delete workflow data. This function is used internally within `patRoon` to implement filtering and subsetting operations, but may also be useful for advanced data processing.

Like the subset operator this function accepts a `i` and `j` parameter to specify which data should be operated on:

Class                   | Argument `i`  | Argument `j`
----------------------- | ------------- | -------------
`features`              | analysis      | feature index
`featureGroups`         | analysis      | feature group
`formulas`, `compounds` | feature group | candidate index
`components`            | component     | feature group

If `i` or `j` is not specified (`NULL`) then data is removed for the complete selection. Some examples are shown below:

```{r del,eval=FALSE}
# delete 2nd feature in analysis-1
fList <- delete(fList, i = "analysis-1", j = 2)
# delete first ten features in all analyses
fList <- delete(fList, i = NULL, j = 1:10)

# completely remove third/fourth analyses from feature groups
fGroups <- delete(fGroups, i = 3:4)
# delete specific feature group
fGroups <- delete(fGroups, j = "M120_R268_30")
# delete range of feature groups
fGroups <- delete(fGroups, j = 500:750)

# remove all results for a feature group
formulas <- delete(formulas, i = "M120_R268_30")

# remove top candidate for all feature groups
compounds <- delete(compounds, j = 1)

# remove a component
components <- delete(components, i = "CMP1")
# remove specific feature group from a component
components <- delete(components, i = "CMP1", j = "M120_R268_30")
# remove specific feature group from all components
components <- delete(components, j = "M120_R268_30")
```

The `j` parameter can also be a function: in this case it is called repeatedly on parts of the data to select what should be deleted. How the function is called and what it should return depends on the workflow data class:

Class                   | Called on every | First argument                   | Second argument    | Return value
----------------------- | --------------- | -------------------------------- | ------------------ | -------------------------------------------
`features`              | analysis        | `data.table` with features       | analysis name      | Features indices (as `integer` or `logical`)
`featureGroups`         | feature group   | vector with group intensities    | feature group name | The analyses of the features to remove (as `character`, `integer`, `logical`)
`formulas`, `compounds` | feature group   | `data.table` with annotations    | feature group name | Candidate indices (rows)
`components`            | component       | `data.table` with the component  | component name     | The feature groups (as `character`, `integer`)

Some examples for this:

```{r delF,eval=FALSE}
# remove features with intensities below 5000
fList <- delete(fList, j = function(f, ...) f$intensity <= 5E3)

# same, but for features in all feature groups from specific analyses
fGroups <- delete(i = 1:3, j = function(g, ...) g <= 5E3)

# remove formula candidates with high relative mass deviation
formulas <- delete(formulas, j = function(ft, ...) ft$error > 5)
```

## Unique and overlapping features {#unOv}

Often an analysis batch is composed of different sample groups, such as different treatments, influent/effluent etc. In such scenarios it may be highly interesting to evaluate uniqueness or overlap between these samples. Furthermore, extracting overlapping or unique features is a simple but effective prioritization strategy.

The `overlap()` and `unique()` functions can be used to extract overlapping and unique features between replicate groups, respectively. Both functions return a subset of the given `featureGroups` object. An overview of their arguments is given below.

Argument     | Function(s)             | Remarks
------------ | ----------------------- | ----------------------------------------------------------------
`which`      | `unique()`, `overlap()` | The replicate groups to compare.
`relativeTo` | `unique()`              | Only return unique features compared to these replicate groups (`NULL` for all). Replicate groups in `which` are ignored.
`outer`      | `unique()`              | If `TRUE` then only return features which are _also_ unique among the compared replicates groups.
`exclusive`  | `overlap`               | Only keep features that _only_ overlap between the compared replicate groups.

Some examples:

```{r compUn,eval=FALSE}
# only keep features uniquely present in replicate group "repl1"
fGroupsUn1 <- unique(fGroups, which = "repl1")
# only keep features in repl1/repl2 which are not in repl3
fGroupsUn2 <- unique(fGroups, which = c("repl1", "repl2"),
                     relativeTo = "repl3")
# only keep features that are only present in repl1 OR repl2
fGroupsUn3 <- unique(fGroups, which = c("repl1", "repl2"),
                     outer = TRUE)

# only keep features overlapping in repl1/repl2
fGroupsOv1 <- overlap(fGroups, which = c("repl1", "repl2"))
# only keep features overlapping in repl1/repl2 AND are not present in any other
# replicate group
fGroupsOv2 <- overlap(fGroups, which = c("repl1", "repl2"),
                      exclusive = TRUE)
```

In addition, several plotting functions are [discussed in the visualization section](#visComp) that visualize overlap and uniqueness of features.

## MS similarity {#specSim}

The _spectral similarity_ is used to compare spectra from different features. For this purpose the `spectrumSimilarity` function can be used. This function operates on [MS peak lists], and accepts the following function arguments:

Argument                   | Remarks
-------------------------- | ---------------------------------------------------------------------------------------
`MSPeakLists`              | The [MS peak lists] object from which peak lists data should be taken.
`groupName1`, `groupName2` | The name(s) of the first and second feature group(s) to compare
`analysis1`, `analysis2`   | The analysis names of the data to be compared. Set this when feature data (instead of feature group data) should be compared.
`MSLevel`                  | The MS level: `1` or `2` for MS and MS/MS, respectively.
`specSimParams`            | Parameters that define how similarities are calculated.
`NAToZero`                 | If `TRUE` then `NA` values are converted to zeros. `NA` values are reported if a comparison cannot be made because of missing peak list data.

The `specSimParams` argument defines the parameters for similarity calculations. It is a `list`, and the default values are obtained with the `getDefSpecSimParams()` function:

```{r defSpecSimParams,eval=TRUE}
getDefSpecSimParams()
```

The `method` field describes the calculation measure: this is either `"cosine"` or `"jaccard"`.

The `shift` field is primarily useful when comparing MS/MS data and defines if and how a spectral shift should be performed prior to similarity calculation:

* `"none"`: The default, no shifting is performed.
* `"precursor"` The mass difference between the precursor mass of both spectra (_i.e._ the feature mass) is first calculated. This difference is then subtracted from each of the mass peaks of the second spectrum. This shifting increases similarity if the MS fragmentation process itself occurs similarly (_i.e._ if both features show similar neutral losses).
* `"both`" This combines both shifting methods: first peaks are aligned that have the same mass, then the `precursor` strategy is applied for the remaining mass peaks. This shifting method yields higher similarities if either fragment masses or neutral losses are similar.

To override a default setting, simply pass it as an argument to `getDefSpecSimParams`:

```{r defSpecSimParamsOv,eval=FALSE}
getDefSpecSimParams(shift = "both")
```

For more details on the various similarity calculation parameters see the reference manual (`?getDefSpecSimParams`).

Some examples are shown below:

```{r specSim,eval=runData}
# similarity between MS spectra with default parameters
spectrumSimilarity(mslists, groupName1 = "M120_R268_30", groupName2 = "M137_R249_53")

# similarity between MS/MS spectra with default parameters
spectrumSimilarity(mslists, groupName1 = "M120_R268_30", groupName2 = "M192_R355_191",
                   MSLevel = 2)

# As above, with jaccard calculation
spectrumSimilarity(mslists, groupName1 = "M120_R268_30", groupName2 = "M192_R355_191",
                   MSLevel = 2, specSimParams = getDefSpecSimParams(method = "jaccard"))

# With shifting
spectrumSimilarity(mslists, groupName1 = "M120_R268_30", groupName2 = "M192_R355_191",
                   MSLevel = 2, specSimParams = getDefSpecSimParams(shift = "both"))
```

The `spectrumSimilarity` function can also be used to calculate _multiple_  similarities. Simply specify multiple feature group names for the `groupNameX` parameters. Alternatively, if you want to compare the same set of feature groups with each other pass their names only as the `groupName1` parameter:

```{r specSimMult,eval=runData}
# compare two pairs
spectrumSimilarity(mslists,
                   groupName1 = c("M120_R268_30", "M137_R249_53"),
                   groupName2 = c("M146_R309_68", "M192_R355_191"),
                   MSLevel = 2, specSimParams = getDefSpecSimParams(shift = "both"))

# compare all
spectrumSimilarity(mslists, groupName1 = groupNames(mslists),
                   MSLevel = 2, specSimParams = getDefSpecSimParams(shift = "both"))
```


## Visualization

### Features and annatation data {#vis_feat_ann}

Several generic functions are available to visualize feature and annotation data:

Generic           | Classes                                              | Remarks
----------------- | ---------------------------------------------------- | ---------------------------------------------------------------
`plot()`          | `featureGroups`, `featureGroupsComparison`           | Scatter plot for retention and _m/z_ values
`plotInt()`       | `featureGroups`                                      | Intensity profiles across analyses
`plotChroms()`    | `featureGroups`, `components`                        | Plot extracted ion chromatograms (EICs)
`plotSpectrum()`  | `MSPeakLists`, `formulas`, `compounds`, `components` | Plots (annotated) spectra
`plotStructure()` | `compounds`                                          | Draws candidate structures
`plotScores()`    | `formulas`, `compounds`                              | Barplot for candidate scoring
`plotGraph()`     | `componentsNT`                                       | Draws interactive graphs of linked homologous series

The most common plotting functions are `plotChroms()`, which plots chromatographic data for features, and `plotSpectrum()`, which will plot (annotated) spectra. An overview of their most important function arguments are shown below.

Argument                                      | Generic                   | Remarks
--------------------------------------------- | ------------------------- | -------------------------------------------------------------------
`rtWindow`                                    | `plotChroms()`                   | Extra time (in s) +/- retention limits of plotted features (useful to zoom out)
`retMin`                                      | `plotChroms()`                   | If `TRUE` plot retention times in minutes
`topMost`                                     | `plotChroms()`                   | Only draw this amount of highest intensity features in each group.
`topMostByRGroup`                             | `plotChroms()`                   | If `TRUE` then the `topMost` parameter specifies the top most intense features in each replicate group to draw (e.g. `topMost=1` would draw the most intense feature for each replicate group). 
`showPeakArea`, `showFGroupRect`              | `plotChroms()`                   | Fill peak areas / draw rectangles around feature groups?
`title`                                       | `plotChroms()`, `plotSpectrum()` | Override plot title
`colourBy`                                    | `plotChroms()`                   | Colour individual feature groups (`"fGroups"`) or replicate groups (`"rGroups"`). By default nothing is coloured (`"none"`)
`showLegend`                                  | `plotChroms()`                   | Display a legend? (only if `colourBy!="none"`)
`onlyPresent`                                 | `plotChroms()`                   | Only plot EICs for analyses where a feature was detected? Setting to `FALSE` is useful to inspect if a feature was 'missed'.
`xlim`, `ylim`                                | `plotChroms()`, `plotSpectrum()` | Override x/y axis ranges, i.e. to manually set plotting range.
`groupName`, `analysis`, `precursor`, `index` | `plotSpectrum()`                 | What to plot. See examples below.
`MSLevel`                                     | `plotSpectrum()`                 | Whether to plot an MS or MS/MS spectrum (only `MSPeakLists`)
`formulas`                                    | `plotSpectrum()`                 | Whether `formula` annotation should be added (only `compounds`)
`plotStruct`                                  | `plotSpectrum()`                 | Whether the structure should be added to the plot (only `compounds`)
`mincex`                                      | `plotSpectrum()`                 | Minimum annotation font size (only `formulas`/`compounds`)

Note that we can use [subsetting](#subset) to select which feature data we want to plot, e.g.

```{r plotChromsSub,eval=runData}
plotChroms(fGroups[1:2]) # only plot EICs from first and second analyses.
plotChroms(fGroups[, 1]) # only plot all features of first group
```

The `plotStructure()` function will draw a chemical structure for a compound candidate. In addition, this function can draw the maximum common substructure (MCS) of multiple candidates in order to assess common structural features.

```{r plotStruct,eval=runData,out.width="25%",fig.show="hold"}
# structure for first candidate
plotStructure(compounds, index = 1, groupName = "M120_R268_30")
# MCS for first three candidates
plotStructure(compounds, index = 1:3, groupName = "M120_R268_30")
```

Some other common and less common plotting operations are shown below.

```{r plotFeatAnn,eval=runData}
plot(fGroups) # simple scatter plot of retention and m/z values

plotChroms(fGroups) # plot EICs for all features
# get overview of all feature groups
plotChroms(fGroups,
           colourBy = "fGroup", # unique colour for each group
           topMost = 1, # only most intense feature in each group
           showPeakArea = TRUE, # show integrated areas
           showFGroupRect = FALSE,
           showLegend = FALSE) # no legend (too busy for many feature groups)
plotChroms(fGroups[, 1], # only plot all features of first group
           colourBy = "rGroup") # and mark them individually per replicate group
plotChroms(components, index = 7, fGroups = fGroups) # EICs from a component


plotSpectrum(mslists, "M120_R268_30") # non-annotated MS spectrum
plotSpectrum(mslists, "M120_R268_30", MSLevel = 2) # non-annotated MS/MS spectrum
# formula annotated spectrum
plotSpectrum(formulas, index = 1, groupName = "M120_R268_30",
             MSPeakLists = mslists)
# compound annotated spectrum, with added formula annotations
plotSpectrum(compounds, index = 1, groupName = "M120_R268_30", MSPeakLists = mslists,
             formulas = formulas, plotStruct = TRUE)
# custom intensity range (e.g. to zoom in)
plotSpectrum(compounds, index = 1, groupName = "M120_R268_30", MSPeakLists = mslists,
             ylim = c(0, 5000), plotStruct = FALSE)
plotSpectrum(components, index = 7) # component spectrum
```

```{r plotHomol,eval=runData}
# Inspect homologous series
plotGraph(componNT)
```


### Overlapping and unique data {#visComp}

There are three functions that can be used to visualize overlap and uniqueness between data:

Generic     | Classes
----------- | ------------------------------------------------------------------- 
`plotVenn`  | `featureGroups`, `featureGroupsComparison`, `formulas`, `compounds` 
`plotUpSet` | `featureGroups`, `featureGroupsComparison`, `formulas`, `compounds` 
`plotChord` | `featureGroups`, `featureGroupsComparison`

The most simple comparison plot is a Venn diagram (i.e. `plotVenn()`). This function is especially useful for two or three-way comparisons. More complex comparisons are better visualized with [UpSet] diagrams (i.e. `plotUpSet()`). Finally, chord diagrams (i.e. `plotChord()`) provide visually pleasing diagrams to assess overlap between data.

These functions can either be used to compare feature data or different objects of the same type. The former is typically used to compare overlap or uniqueness between features in different replicate groups, whereas comparison between objects is useful to visualize differences in algorithmic output. Besides visualization, note that both operations can also be performed to modify or combine objects (see [unique and overlapping features](#unOv) and [algorithm consensus](#consensus)).

As usual, some examples are shown below.

```{r prePlotComp,include=FALSE,eval=runData}
fGroupsO <- fGroups
fGroups <- fGroupsRG
```

```{r plotComp,eval=runData}
plotUpSet(fGroups) # compare replicate groups
plotVenn(fGroups, which = c("repl1", "repl2")) # compare some replicate groups
plotChord(fGroups, average = TRUE) # overlap between replicate groups
# compare with custom made groups
plotChord(fGroups, average = TRUE,
          outer = c(repl1 = "grp1", repl2 = "grp1", repl3 = "grp2", repl4 = "grp3"))

# compare GenForm and SIRIUS results
plotVenn(formsGF, formsSIR,
         labels = c("GF", "SIR")) # manual labeling
```

```{r postPlotComp,include=FALSE,eval=runData}
fGroups <- fGroupsO
```

### MS similarity

The `plotSpectrum` function is also useful to visually compare (annotated) spectra. This works for `MSPeakLists`, `formulas` and `compounds` object data.

```{r plotSim,eval=runData}
plotSpectrum(mslists, groupName = c("M120_R268_30", "M137_R249_53"), MSLevel = 2)
plotSpectrum(compounds, groupName = c("M120_R268_30", "M146_R309_68"), index = c(1, 1),
             MSPeakLists = mslists)
```

The `specSimParams` argument, which was discussed in [MS similarity](#specSim), can be used to configure the similarity calculation:

```{r plotSimParam,eval=runData}
plotSpectrum(mslists, groupName = c("M120_R268_30", "M137_R249_53"), MSLevel = 2,
             specSimParams = getDefSpecSimParams(shift = "both"))
```

### Hierarchical clustering results {#plotClust}

In `patRoon` hierarchical clustering is used [for some componentization algorithms](#compClust) and to cluster candidate compounds with similar chemical structure (see [compound clustering](#compclust)). The functions below can be used to visualize their results.

Generic             | Classes                                  | Remarks
------------------- | ---------------------------------------- | --------------------------------------------------
`plot()`            | All                                      | Plots a dendrogram
`plotInt()`         | `componentsIntClust`                     | Plots normalized intensity profiles in a cluster
`plotHeatMap()`     | `componentsIntClust`                     | Plots an heatmap
`plotSilhouettes()` | `componentsClust`                        | Plot silhouette information to determine the cluster amount
`plotStructure()`   | `compoundsCluster`                       | Plots the maximum common substructure (MCS) of a cluster


```{r plotClust,eval=runData}
plot(componInt) # dendrogram
plot(compsClust, groupName = "M120_R268_30") # dendrogram for clustered compounds
plotInt(componInt, index = 4) # intensities of 4th cluster
```

```{r plotClust2,eval=runData,fig.height=6}
plotHeatMap(componInt) # plot heatmap
```

```{r plotClust3,eval=runData}
plotHeatMap(componInt, interactive = TRUE) # interactive heatmap (with zoom-in!)
plotSilhouettes(componInt, 5:20) # plot silhouettes (e.g. to obtain ideal cluster amount)
```

### Generating EICs in DataAnalysis

If you have Bruker data and the DataAnalysis software installed, you can automatically add EIC data in a DataAnalysis session. The `addDAEIC()` will do this for a single _m/z_ in one analysis, whereas the `addAllDAEICs()` function adds EICs for all features in a `featureGroups` object.

```{r addDAEIC,eval=FALSE}
# add a single EIC with background subtraction
addDAEIC("mysample", "~/path/to/sample", mz = 120.1234, bgsubtr = TRUE)
# add TIC for MS/MS signal of precursor 120.1234 (value of mz is ignored for TICs)
addDAEIC("mysample", "~/path/to/sample", mz = 100, ctype = "TIC",
         mtype = "MSMS", fragpath = "120.1234", name = "MSMS 120")

addAllDAEICs(fGroups) # add EICs for all features
addAllDAEICs(fGroups[, 1:50]) # as usual, subsetting can be used for partial data
```

## Interactively explore and review data {#intReview}

The `checkFeatures` and `checkComponents` functions start a graphical user interface (GUI) which allows you to interactively explore and review feature and components data, respectively.

```{r eval=FALSE}
checkFeatures(fGroups) # inspect features and feature groups
checkComponents(componCAM, fGroups) # inspect components
```

Both functions allow you to easily explore the data in an interactive way. Furthermore, these functions allow you to remove unwanted data. This is useful to remove for example features that are actually noise and feature groups that shouldn't be in the same component. To remove an unwanted feature, feature group or components, simply uncheck its 'keep' checkbox. The next step is to save the selections you made. A _check session_ is a file that stores which data should be removed. Once the session file is saved the `filter` function can be used to actually remove the data:

```{r eval=FALSE}
fGroupsF <- filter(fGroups, checkFeaturesSession = TRUE)
componCAMF <- filter(componCAM, checkComponentsSession = TRUE)
```

If you saved the session and you re-launch the GUI it will restore the selections made earlier. The `clearSession` argument can be used to fully clear a session before starting the GUI, hence, all the data will be restored to their 'keep state'.

```{r eval=FALSE}
checkFeatures(fGroups, clearSession = TRUE) # start GUI with fresh session
```

It is also possible to use multiple different sessions. This is especially useful if you do not want to overwrite previous session data or want to inspect different objects. In this case the session file name should be specified:

```{r eval=FALSE}
checkFeatures(fGroups, "mysession.yml")
fGroupsF <- filter(fGroups, checkFeaturesSession = "mysession.yml")
```

The default session names are `"checked-features.yml"` and `"checked-components.yml"` for feature and component data, respectively.

The extension of session file names is `.yml` since the [YAML] file format is used. An advantage of this format is that it is easily readable and editable with a text editor.

Note that the session data is tied to the feature group names of your data. This means that, for instance, when you re-group your feature data after changing some parameters, the session data you prepared earlier cannot be used anymore. Since probably quite some manual work went into creating the session file, a special function is available to import a session that was made for previous data. This function tries its best to guess the new feature group name based on similarity of their retention times and _m/z_ values.

```{r eval=FALSE}
checkFeatures(fGroups) # do manual inspection

fGroups <- groupFeatures(fList, ...) # re-group with different parameters

importCheckFeaturesSession("checked-features.yml", "checked-features-new.yml", fGroups)

checkFeatures(fGroups, session = "checked-features-new.yml") # inspect new data
```

Take care to monitor the messages that `importCheckFeaturesSession` may output, as it may be possible that some 'old' feature groups are not found or are matched by multiple candidates of the new dataset.

Some additional parameters exist to the functions described in this section. As usualy check the reference manual for more details (_e.g._ `?checkFeatures`).

> **_NOTE_** Although the GUI tools described here allow you to easily filter out results, it is highly recommended to first prioritize your data to avoid doing a lot of unneeded manual work.

## Reporting {#report}

The previous sections showed various functionalies to inspect and visualize results. An easy and automated way to do this automatically is by using the _reporting_ functionality of `patRoon`. The following three reporting functions are available:

* `reportCSV()`: exports workflow data to comma-separated value (csv) files
* `reportPDF()`: generates simple reports by plotting workflow data in portable document files (PDFs)
* `reportHTML()`: generates interactive and easily explorable reports

There are many different arguments available to configure the reporting process. Some common arguments are listed below; for a complete listing see the reference manual (e.g. `?reporting`).

Argument                              | Functions                     | Remarks
------------------------------------- | ----------------------------- | ----------------------------------------------------------
`fGroups`, `formulas`, `compounds`, `formulas`, `components`, `compsCluster` | All | Objects to plot. Only `fGroups` is mandatory.
`MSPeakLists`                         | `reportPDF()`, `reportHTML()` | The `MSPeakLists` object that was used to generate annotation data. Only needs to be specified if `formulas` or `compounds` are reported.
`path`                                | All                           | Directory path where report files will be stored (`"report"` by default).
`formulasTopMost`, `compoundsTopMost` | `reportPDF()`, `reportHTML()` | Report no more than this amount of highest ranked candidates.
`EICOnlyPresent`                      | `reportPDF()`, `reportHTML()` | Only plot an EIC for an analysis if a feature was detected.
`selfContained`                       | `reportHTML()`                | Outputs to a single and self contained `.html` file. Handy to share reports, but not recommended for large amounts of data.

Which data will be reported is fully configurable. The only workflow object that must be specified are the feature groups (i.e. with the `fGroups` argument), all other data (e.g. `compounds`, `components`) are optional. This means that reporting can be performed at every stage during the workflow, which, for instance, can be useful to quickly inspect results when testing out various settings to generate workflow data.

When formula or compound results are reported with `reportPDF()` or `reportHTML()` then only the top ranked candidates are considered. This limtation is often necessary as reporting many candidates will take considerable time. By default the top 5 for each feature group are reported, however, this number can be changed with the `formulasTopMost` and `compoundsTopMost` arguments.

Some typical examples:

```{r reporting,eval=FALSE}
reportHTML(fGroups) # simple interactive report with feature data
# generate PDFs with feature and compound annotation data
reportPDF(fGroups, compounds = compounds, MSPeakLists = mslists)
reportCSV(fGroups, path = "myReport") # change destination path

# generate report with all workflow types and increase maximum number of
# compound candidates to top 10
reportHTML(fGroups, formulas = formulas, compounds = compounds,
           components = components, MSPeakLists = mslists,
           compsCluster = compsClust,
           compoundsTopMost = 10)
```

```{r child=file.path(vignDir, "shared", "_refs.Rmd")}
```
# Workflow concepts

<!-- UNDONE: -->
<!-- - usage of algorithms: algorithm arg vs specific fun -->
<!-- - S4 classes?: common output, hierarchy, generics, ... -->
<!-- - show images of feature definition etc? -->
    
In a non-target workflow both chromatographic and mass spectral data is automatically processed in order to provide a comprehensive chemical characterization of your samples. While the exact workflow is typically dependent on the type of study, it generally involves of the following steps:
    
<!-- UNDONE: include other data processing steps here? -->
    
```{r workflow,echo=FALSE,out.width="100%"}
plotGV("
digraph workflow {
  graph [ rankdir = LR, compound = true, style = invis ]
  node [ shape = box,
         fixedsize = true,
         width = 2.8,
         height = 1,
         fontsize = 20,
         fillcolor = darkseagreen1,
         style = filled ]

  subgraph cluster1 {
    'Data pre-treatment' -> 'Find features'
    'Find features' -> 'Group features'
  }

  subgraph cluster2 {
    'Suspect screening' -> 'Group features' [minlen=3, style=dashed]
    'Find features' -> 'Suspect screening' [style=invis]
    'Group features' -> 'Suspect screening' [style=dashed]
  }

  subgraph cluster3 {
    Componentization
  }
  
  subgraph cluster4 {
    graph [style = dashed ]
    color=blue
    'MS peak lists' 'Formula annotation'  'Compound Annotation'
  }

  'Group features' -> 'Componentization' [ltail=cluster1, lhead=cluster3, style=dashed, dir=both, minlen=3]
  'Group features' -> 'Formula annotation' [ltail=cluster1, lhead=cluster4, style=dashed, minlen=3]
  
}", height = 250, width = 750)
```

Note that `patRoon` supports flexible composition of workflows. In the scheme above you can recognize optional steps by a _dashed line_. The inclusion of each step is only necessary if a further steps depends on its data. For instance, annotation and componentization do not depend on each other and can therefore be executed in any order or simply be omitted. A brief description of all steps is given below.

During **data pre-treatment** raw MS data is prepared for further analysis. A common need for this step is to convert the data to an open format so that other tools are able to process it. Other pre-treatment steps may involve re-calibration of _m/z_ data or performing advanced filtering operations.

The next step is to extract **features** from the data. While different terminologies are used, a feature in `patRoon` refers to a single chromatographic peak in an extracted ion chromatogram for a single _m/z_ value (within a defined tolerance). Hence, a feature contains both chromatographic data (e.g. retention time and peak height) and mass spectral data (e.g. the accurate _m/z_). Note that with mass spectrometry multiple _m/z_ values may be detected for a single compound as a result of adduct formation, natural isotopes and/or in-source fragments. Some algorithms may try to combine these different masses in a single feature. However, in `patRoon` we generally assume this is not the case (and may optionally be done afterwards during the componentization step described below). Features are sometimes simply referred to as 'peaks'.

Features are found per analysis. Hence, in order to compare a feature across analyses, the next step is to group them. This step is essential as it finds equal features even if their retention time or _m/z_ values slightly differ due to analytical variability. The resulting **feature groups** are crucial input for subsequent workflow steps. Prior to grouping, _retention time alignment_ between analyses may be performed to improve grouping of features, especially when processing multiple analysis batches at once. Outside `patRoon` feature groups may also be defined as _profiles_, _aligned_ or _grouped features_ or _buckets_.

Depending on the study type, **suspect screening** is then performed to limit the features that should be considered for further processing. As its name suggests, with suspect screening only those features which are suspected to be present are considered for further processing. These suspects are retrieved from a suspect list which contains the _m/z_ and (optionally) retention times for each suspect. Typical suspect lists may be composed from databases with known pollutants or from predicted transformation products. Note that for a 'full' non-target analysis no suspect screening is performed, hence, this step is simply omitted and all features are to be considered.

The feature group data may then be subjected to **componentization**. A **component** is defined as a collection of multiple feature groups that are somehow related to each other. Typical examples are features that belong to the same chemical compound (i.e. with different _m/z_ values but equal retention time), such as adducts, isotopes and in-source fragments. Other examples are homologous series and features that display a similar intensity trend across samples. If adducts or isotopes were annotated during componentization then this data may be used to prioritize the feature groups.

The last step in the workflow commonly involves **annotation**. During this step MS and MS/MS data are collected in so called **MS peak lists**, which are then used as input for formula and compound annotation. Formula annotation involves automatic calculation of possible formulae for each feature based on its _m/z_, isotopic pattern and MS/MS fragments, whereas compound annotation (or identification) involves the assignment of actual chemical structures to each feature. Note that during formula and compound annotation typically multiple candidates are assigned to a single feature. To assist interpretation of this data each candidate is therefore ranked on characteristics such as isotopic fit, number of explained MS/MS fragments and metadata from an online database such as number of scientific references or presence in common suspect lists.

To summarize:
    
* **Data-pretreatment** involves preparing raw MS data for further processing (e.g. conversion to an open format)
* **Features** describe chromatographic and _m/z_ information (or 'peaks') in all analyses.
* A **feature group** consists of equal features across analyses.
* With **suspect screening** only features that are considered to be on a suspect list are considered further in the workflow.
* **Componentization** involves consolidating different feature groups that have a relationship to each other in to a single component.
* **MS peak lists** Summarizes all MS and MS/MS data that will be used for subsequent annotation.
* During **formula** and **compound annotation** candidate formulae/structures will be assigned and ranked for each feature.

The next chapters will discuss how to generate this data and process it. Afterwards, several advanced topics are discussed such as [combining positive and negative ionization data](#setsWorkflow), [screening for transformation products](#TPs) and [other advanced functionality](#advanced_usage).

```{r child=file.path(vignDir, "shared", "_refs.Rmd")}
```
# Transformation product screening {#TPs}

This chapter describes the various functionality for screening of _transformation products_ (TPs), which are introduced since `patRoon` 2.0. Screening for TPs, i.e. chemicals that are formed from a _parent_ chemical by e.g. chemical or biological processes, has broad applications. For this reason, the TP screening related functionality is designed to be flexible, thus allowing one to use a workflow that is best suited for a particular study.

Regardless, the TP screening workflow in `patRoon` can be roughly summarized as follows:

```{r TPWorkflow,echo=FALSE,out.width="100%"}
plotGV("
digraph Workflow {
  graph [ rankdir = LR ]
  node [ shape = box,
         fixedsize = true,
         width = 2.2,
         height = 1,
         fontsize = 18,
         fillcolor = darkseagreen1,
         style = filled ]

    'Parent screening' -> 'Obtaining TPs' -> 'TP screening' -> 'Linking parent/TPs'
}", height = 90, width = 500)
```

* **Parent screening** During this step a common `patRoon` workflow is used to screen for the parent chemicals of interest. This could be a full non-target analysis with compound annotation or a relative simple suspect or target screening.
* **Obtaining TPs** Data is obtained of potential TPs for the parents of interest. The TPs may originate from a library or predicted _in-silico_. Note that in some workflows this step is omitted (discussed later).
* **TP screening** A suspect screening is performed to find the TPs in the analysis data. 
* **Linking parents and TPs** In the step the parent features are linked with the TP features. Several post-processing functionality exists to improve and prioritize the data.

The next sections will outline more details on these steps are performed and configured. The [last section](#TPsExamples) in this chapter outlines several example workflows.

> **_NOTE_** The [newProject tool](#newProject) can be used to easily generate a workflow with transformation product screening.

## Obtaining transformation product data {#genTPs}

The `generateTPs` function is used to obtain TPs for a particular set of parents. Like other workflow generator functions (`findFeatures`, `generateCompounds`), several algorithms are available that do the actual work.

Algorithm        | Usage                                            | Remarks
---------------- | ------------------------------------------------ | ---------------------------------
[BioTransformer] | `generateTPs(algorithm = "biotransformer", ...)` | Predicts TPs with full structural information
Library          | `generateTPs(algorithm = "library", ...)`        | Obtains transformation products from a library ([PubChem transformations][PubChemLiteTR] or custom)
Metabolic logic  | `generateTPs(algorithm = "logic", ...)`          | Uses pre-defined logic to predict TPs based on common elemental differences (e.g. hydroxylation, demethylation). Based on @Scholle2015.

The `biotransformer` and `library` algorithms provide full structural information of the TPs (e.g. formula, SMILES, predicted Log P). However, these algorithms also depend on the full chemical structure of the parent compound. Hence, these algorithms are typically suitable when parents are known in advance or were found by a suspect screening. On the other hand, metabolic logic only requires the feature mass, and this simplicity allows it to predict TPs for all features. This algorithm is most suitable for full non-target analysis, however, extra care must be taken to rule out false positives.

An overview of common arguments for TP generation is listed below.

Argument                      | Algorithm(s)      | Remarks
----------------------------- | ----------------- | --------------------------------------------------------
`parents`                     | `biotransformer`, `library` | The input parents. See section below.
`fGroups`                     | `logic`           | The input feature groups to calculate TPs for.
`type`                        | `biotransformer`  | The prediction type: `"env"`, `"ecbased"`, `"cyp450"`, `"phaseII"`, `"hgut"`, `"superbio"`, `"allHuman"`. See [BioTransformer] for more details.
`TPLibrary`/`transformations` | `library`/`logic` | [Custom TP library/transformation rules](#TPsCustom).
`adduct`                      | `logic`           | The assumed adduct of the parents (e.g. `"[M+H]+"`). Not needed when [adduct annotations](#incorpAdductIso) are available.

### Parent input

The input parent structures for the `biotransformer` and `library` algorithms must be one the following:

* A suspect list (follows the same format as [suspect screening](#suspscr))
* A feature groups object with screening results (e.g. obtained with `screenSuspects`, see [suspect screening](#suspscr))
* A `compounds` object obtained with [compound annotation](#compounds)

In the former two cases the parent information is taken from the suspect list or from the hits in a suspect screening worklow, respectively. The last case is more suitable for when the parents are not completely known. In this case, the candidate structures from a [compound annotation](#compounds) are used as input to obtain TPs. Since _all_ the candidates are used, it is highly recommend to filter the object in advance, for instance, with the `topMost` filter. For `library`, the parent input is optional: if no parents are specified then TP data for _all_ parents in the database is used.

For the `logic` algorithm TPs are predicted directly for feature groups. Since this algorithm can only perform very basic validity checks, it is strongly recommended to first prioritize the feature group data.

Some typical examples:

```{r eval=FALSE}
# predict environmental TPs with BioTransformer for all parents in a suspect list
TPsBT <- generateTPs("biotransformer", parents = patRoonData::suspectsPos,
                     type = "env")
# obtain all TPs from the default library
TPsLib <- generateTPs("library")
# get TPs for the parents from a suspect screening
TPsLib <- generateTPs("library", parents = fGroupsScr)
# calculate TPs for all feature groups
TPsLogic <- generateTPs("logic", fGroups, adduct = "[M+H]+")
```

### Processing data

Similar to other workflow data, several generic functions are available to inspect the data:

Generic                            | Classes | Remarks                                                             
-----------------------------------|---------|---------------------------------------------------------------------
`length()`                         | All     | Returns the total number of transformation products
`names()`                          | All     | Returns the names of the parents
`parents()`                        | All     | Returns a table with information about the parents
`products()`                       | All     | Returns a `list` with for each parent a table with TPs
`as.data.table()`, `as.data.frame` | All     | Convert all the object information into a `data.table`/`data.frame`
`"[["` / `"$"` operators           | All     | Extract TP information for a specified parent

Some examples:

```{r include=FALSE,eval=TRUE}
# NOTE: this is always evaluated as it takes very little time...
TPs <- generateTPs("library")
```

```{r TPsProcInsp,eval=TRUE}
# just show a few columns in this example, there are many more!
# note: the double dot syntax (..cols) is necessary since the data is stored as data.tables
cols <- c("name", "formula", "InChIKey")
parents(TPs)[1:5, ..cols]
TPs[["DEET"]][, ..cols]
TPs[[2]][, ..cols]
as.data.table(TPs)[1:5, 1:3]
```

In addition, the following generic functions are available to modify or convert the object data:

Generic             | Classes                    | Remarks
------------------- | -------------------------- | --------------------------------------------------------
`"["` operator      | All                        | Subset this object on given parents
`filter`            | `transformationProductsBT` | Filters this object
`convertToSuspects` | All                        | Generates a suspect list of all TPs (and optionally parents) that is suitable for `screenSuspects`
`convertToMFDB`     | `transformationProductsBT`, `transformationProductsLibrary` | Generates a [MetFrag] database for all TPs (and optionally parents)

```{r TPsProcMod,eval=FALSE}
TPs2 <- TPs[1:10] # only keep results for first ten parents

# remove transformation products that are isomers to their parent or sibling TPs
# may simplify data as these are often difficult to identify
TPsF <- filter(TPs, removeParentIsomers = TRUE, removeTPIsomers = TRUE)

# remove duplicate transformation products from each parent
# these can occur if different pathways yield the same TPs
TPsF <- filter(TPs, removeDuplicates = TRUE)

# only keep TPs that have a structural similarity to their parent of >= 0.5
TPsF <- filter(TPs, minSimilarity = 0.5)

# do a suspect screening for all TPs and their parents
suspects <- convertToSuspects(TPs, includeParents = TRUE)
fGroupsScr <- screenSuspects(fGroups, suspects, onlyHits = TRUE)

# use the TP data for a specialized MetFrag database
convertToMFDB(TPs, "TP-database.csv", includeParents = FALSE)
compoundsTPs <- generateCompounds(fGroups, mslists, "metfrag", database = "csv",
                                  extraOpts = list(LocalDatabasePath = "TP-database.csv"))
```

The `convertToSuspects` function is always part of a workflow with `biotransformer` or `library` TPs. This is discussed further in the next section. The `convertToMFDB` function is especially handy with `biotransformer` workflows, as it allows generating a compound database for TPs that may not be available in other databases. This is further demonstrated in the [first example](#TPsEx1).

### Custom libraries and transformations {#TPsCustom}

By default the `library` and `logic` algorithms use data that is installed with `patRoon` (based on [PubChem transformations][PubChemLiteTR] and @Scholle2015, respectively). However, it is also possible to use custom data.

To use a custom TP library a simple `data.frame` is needed with the names, SMILES and optionally `log P` values for the parents and TPs. The `log P` values are used for prediction of the retention time direction of a TP compared to its parent, as is discussed further in the next section. The following small library has two TPs for benzotriazole and one for DEET:

```{r TPsCustomDB,eval=TRUE}
myTPLib <- data.frame(parent_name = c("1H-Benzotriazole", "1H-Benzotriazole", "DEET"),
                      parent_SMILES = c("C1=CC2=NNN=C2C=C1", "C1=CC2=NNN=C2C=C1", "CCN(CC)C(=O)C1=CC=CC(=C1)C"),
                      TP_name = c("1-Methylbenzotriazole", "1-Hydroxybenzotriazole", "N-ethyl-m-toluamide"),
                      TP_SMILES = c("CN1C2=CC=CC=C2N=N1", "C1=CC=C2C(=C1)N=NN2O", "CCNC(=O)C1=CC=CC(=C1)C"))
myTPLib
```

To use this library, simply pass it to the `TPLibrary` argument:

```{r eval=FALSE}
TPs <- generateTPs("library", TPLibrary = myTPLib)
```

Similarly, for `logic` a table with custom transformation rules can be specified for TP calculations:

```{r TPCustomLogic,eval=TRUE}
myTrans <- data.frame(transformation = c("hydroxylation", "demethylation"),
                      add = c("O", ""),
                      sub = c("", "CH2"),
                      retDir = c(-1, -1))
myTrans
```

The `add` and `sub` columns are used to denote the elements that are added or subtracted by the reaction. These are used to calculate mass differences between parents and TPs. The `retDir` column is used to indicate the retention time direction of the parent compared to the TP: `-1` (elutes before parent), `1` (elutes after parent) or `0` (similar or unknown). The next section describes how this data can be used to filter TPs. The custom rules can be used by passing them to the `transformations` argument:

```{r eval=FALSE}
TPs <- generateTPs("logic", fGroups, adduct = "[M+H]+", transformations = myTrans)
```

## Linking parent and transformation product features

This section discusses one of the most important steps in a TP screening workflow, which is to link feature groups of parents with those of candidate transformation products. During this step, _components_ are made, where each component consist of one or more feature groups of detected TPs for a particular parent. Note that componentization was [already introduced before](#componentization), but for very different algorithms. However, the data format for TP componentization is highly similar. After componentization, several filters are available to clean and prioritize the data. These can even allow workflows without obtaining potential TPs in advance, which is discussed in the last subsection.

### Componentization

Like [other algorithms](#components), the `generateComponents` generic function is used to generate TP components, by setting the `algorithm` parameter to `"tp"`.

The following arguments are of importance:

Argument        | Remarks
--------------- | --------------------------------------------------------------
`fGroups`       | The input feature groups for the _parents_
`fGroupsTPs`    | The input feature groups for the _TPs_
`ignoreParents` | Set to `TRUE` to ignore feature groups in `fGroupsTPs` that also occur in `fGroups`
`TPs`           | The input transformation products, ie as generated by `generateTPs()`
`MSPeakLists`, `formulas`, `compounds` | Annotation objects used for similarity calculation between the parent and its TPs
`minRTDiff`     | The minimum retention time difference (seconds) of a TP for it to be considered to elute differently than its parent.

#### Feature group input {#TPsFGroups}

The `fGroups`, `fGroupsTPs` and `ignoreParents` arguments are used by the componentization algorithm to identify which feature groups can be considered as parents and which as TPs. Three scenarios are possible:

1. `fGroups=fGroupsTPs` and `ignoreParents=FALSE`: in this case no distinction is made, and all feature groups are considered a parent or TP (default if `fGroupsTPs` is not specified).
2. `fGroups` and `fGroupsTPs` contain different subsets of the _same_ `featureGroups` object and `ignoreParents=FALSE`:  only the feature groups in `fGroups`/`fGroupsTPs` are considered as parents/TPs.
3. As above, but with `ignoreParents=TRUE`: the same distinction is made as above, but any feature groups in `fGroupsTPs` are ignored if also present in `fGroups`.

The first scenario is often used if it is unknown which feature groups may be parents or which are TPs. Furthermore, this scenario may also be used if the dataset is sufficiently simple, for instance, because a suspect screening with the results from `convertToSuspects` (discussed in the previous section) would reliably discriminate between parents and TPs. A workflow with the first scenario is demonstrated in the [second example](#TPsEx2).

In all other cases it is recommended to use either the second or third scenario, since making a prior distinction between parent and TP feature groups greatly simplifies the dataset and reduces false positives. A relative simple example where this can be used is when there are two sample groups: before and after treatment.

```{r eval=FALSE}
componTP <- generateComponents(algorithm = "tp",
                               fGroups = fGroups[rGroups = "before"],
                               fGroupsTPs = fGroups[rGroups = "after"])
```

In this example, only those feature groups present in the "before" replicate group are considered as parents, and those in "after" may be considered as a TP. Since it is likely that there will be some overlap in feature groups between both sample groups, the `ignoreParents` flag can be used to not consider any of the overlap for TP assignments:

```{r eval=FALSE}
componTP <- generateComponents(algorithm = "tp",
                               fGroups = fGroups[rGroups = "before"],
                               fGroupsTPs = fGroups[rGroups = "after"],
                               ignoreParents = TRUE)
```

More sophisticates ways are of course possible to provide an upfront distinction between parent/TP feature groups. In the [fourth example](#TPsEx4) a workflow is demonstrated where fold changes are used.

> **_NOTE_** The feature groups specified for `fGroups`/`fGroupsTPs` _must_ always originate from the same `featureGroups` object.

For the `library` and `biotransformer` algorithms it is mandatory that a suspect screening of parents and TPs is performed prior to componentization. This is necessary for the componentization algorithm to map the feature groups that belong to a particular parent or TP. To do so, the `convertToSuspects` function is used to prepare the suspect list:

```{r eval=FALSE}
# set includeParents to TRUE since both the parents and TPs are needed
suspects <- convertToSuspects(TPs, includeParents = TRUE)
fGroupsScr <- screenSuspects(fGroups, suspects, onlyHits = TRUE)

# do the componentization
# a similar distinction between fGroups/fGroupsScr as discussed above can of course also be done
componTP <- generateComponents(fGroups = fGroupsScr, ...)
```

If a parent screening was already performed in advance, for instance when the input parents to `generateTPs` are screening results, the screening results for parents and TPs can also be combined. The [second example](#TPsEx2) demonstrates this.

Note that in the case a parent suspect is matched to multiple feature groups, a component is made for each match. Similarly, if multiple feature groups match to a TP suspect, all of them will be incorporated in the component.

When TPs were generated with the `logic` algorithm a suspect screening must also be carried out in advance. However, in this case it is not necessary to include the parents (since each parent equals a feature group no mapping is necessary). The `onlyHits` variable to `screenSuspects` must not be set in order to keep the parents.

```{r eval=FALSE}
# only screen for TPs
suspects <- convertToSuspects(TPs, includeParents = FALSE)
# but keep all other feature groups as these may be parents
fGroupsScr <- screenSuspects(fGroups, suspects, onlyHits = FALSE)

# do the componentization...
```

#### Annotation similarity calculation

If additional annotation data for parents and TPs is given to the componentization algorithm, it will be used to calculate various similarity properties. Often, the chemical structure for a transformation product is similar to that of its parent. Hence, there is a good chance that a parent and its TPs also share similar MS/MS data.

Firstly, if MS peak lists are provided, then the [spectrum similarity](#specSim) is calculated between each parent and its potential TP candidates. This is performed with all the three different alignment shifts (see the [spectrum similarity section](#specSim) for more details).

In case `formulas` and/or `compounds` objects are specified, then a parent/TP comparison is made by counting the number of fragments and neutral losses that they share (by using the formula annotations). This property is mainly used for non-target workflows where the identity for a parent and TP is not yet well established. For this reason, fragments and neutral losses reported for _all_ candidates for the parent/TP feature group are considered. Hence, it is highly recommend to pre-treat the annotation objects, for instance, with the `topMost` filter. If both `formulas` and `compounds` are given the results are pooled. Note that each unique fragment/neutral loss is only counted once, thus multiple formula/compound candidates with the same annotations will not skew the results.

### Processing data {#TPsProc}

The output of TP componentization is an object of the `componentsTPs` class. This _derives_ from the 'regular' `components` class, therefore, all the data processing functionality described [before](#processing) (extraction, subsetting, filtering etc) are also valid for TP components.

Several additional filters are available to prioritize the data:

Filter        | Remarks
------------- | -----------------------------------
`retDirMatch` | If `TRUE` only keep TPs with an expected chromatographic retention direction compared to the parent.
`minSpecSim`, `minSpecPrec`, `minSpecSimBoth` | The minimum spectrum similarity between the parent and TP. Calculated with no, `"precursor"` and `"both"` alignment shifting (see [spectrum similarity](#specSim)).
`minFragMatches`, `minNLMatches` | Minimum number of formula fragment/neutral loss matches between parent and TP (discussed in previous section).
`formulas`    | A `formulas` object used to further verify candidate TPs that were generated by the `logic` algorithm.

The `retDirMatch` filter compares the expected and observed _retention time direction_ of a TP in order to decide if it should be kept. The direction is a value of either `-1` (TP elutes before parent), `+1` (TP elutes after parent) or `0` (TP elutes very close to the parent or its direction is unknown). The directions are taken from the [generated transformation products](#genTPs). For the `library` and `biotransformer` algorithms the log P values are compared of a TP and its parent. Here, it is assumed that lower log P values result in earlier elution (i.e. typical with reversed phase LC). For the `logic` algorithm the retention time direction is taken from the transformation rules table. Note that specifying a large enough value for the `minRTDiff` argument to `generateComponents` is important to ensure that some tolerance exists while comparing retention time directions of parent and TPs. This filter does nothing if either the observed or expected direction is zero.

When TPs data was generated with the `logic` algorithm it is recommended to use the `formulas` filter. This filter uses formula annotations to verify that (1) a parent feature group contains the elements that are subtracted during the transformation and (2) the TP feature group contains the elements that were added during the transformation. Since the 'right' candidate formula is most likely not yet known, this filter looks at _all_ candidates. Therefore, it is recommended to filter the `formulas` object, for instance, with the `topMost` filter.

### Omitting transformation product input

The `TPs` argument to `generateComponents` can also be omitted. In this case every feature group of `fGroupTPs` is considered to be a potential TP for the potential parents specified for `fGroups`. An advantage is that the screening workflow is not limited to any known TPs or transformations. However, such a workflow has high demands on prioritiation steps before and after the componentization to rule out the many false positives that may occur.

When no transformation data is supplied it is crucial to make [a prior distinction](#TPsFGroups) between parent and TP feature groups. Afterwards, the MS/MS spectral and other annotation similarity filters mentioned in the previous section may be a powerful way to further prioritize data.

The [fourth example](#TPsEx4) demonstrates such a workflow.

### Reporting TP components

The TP components can be reported with the `reportHTML` function. This is done by setting the `components` function argument (i.e. equally to all other component types). The results will be displayed with a customized format that allows easy exploring of each parent with its TPs.

```{r eval=FALSE}
reportHTML(fGroups, components = componTP)
```

## Example workflows {#TPsExamples}

The next subsections demonstrate several approaches to perform a TP screening workflow with `patRoon`. In all examples it is assumed that feature groups were already obtained (with the `findFeatures` and `groupFeatures` functions) and stored in the `fGroups` variable.

The workflows with `patRoon` are designed to be flexible, and the examples here are primarily meant to implement your own workflow. Furthermore, some of the techniques used in the examples can also be combined. For instance, the [Fold change](#FCCalc) classification and MS/MS similarity filters applied in the [fourth example](#TPsEx4) could also be applied to any of the other examples.

### Screen predicted TPs for targets {#TPsEx1}

The first example is a simple workflow where TPs are predicted for a set of given parents with [BioTransformer] and subsequently screened. A [MetFrag] compound database is generated and used for annotation.

```{r eval=FALSE}
# predict TPs for a fixed list of parents
TPs <- generateTPs("biotransformer", parents = patRoonData::suspectsPos)

# screen for the TPs
suspectsTPs <- convertToSuspects(TPs, includeParents = FALSE)
fGroupsTPs <- screenSuspects(fGroups, suspectsTPs, adduct = "[M+H]+", onlyHits = TRUE)

# perform annotation of TPs
mslistsTPs <- generateMSPeakLists(fGroupsTPs, "mzr")
convertToMFDB(TPs, "TP-database.csv", includeParents = FALSE) # generate MetFrag database
compoundsTPs <- generateCompounds(fGroupsTPs, mslistsTPs, "metfrag", adduct = "[M+H]+", database = "csv",
                                  extraOpts = list(LocalDatabasePath = "TP-database.csv"))
```

### Screening TPs from a library for suspects {#TPsEx2}

In this example TPs of interest are obtained for the parents that surfaced from of a suspect screening. The steps of this workflow are:

1. Suspect screening parents.
2. Obtain TPs for the suspect hits from a library.
3. A second suspect screening is performed for TPs and the original parent screening results are amended. Note that the parent data is needed for componentization.
4. Both parents and TPs are annotated using a database generated from their chemical structures.
5. Some prioritization is performed by
    a. Only keeping candidate structures for which _in-silico_ fragmentation resulted in at least one annotated MS/MS peak.
    b. Only keeping suspect hits with an estimated identification level of 3 or better.
6. The TP components are made and only feature groups with parent/TP assignments are kept.
7. All results are reported.

```{r eval=FALSE}
# step 1
fGroupsScr <- screenSuspects(fGroups, patRoonData::suspectsPos, adduct = "[M+H]+")
# step 2
TPs <- generateTPs("library", parents = fGroupsScr)

# step 3
suspects <- convertToSuspects(TPs)
fGroupsScr <- screenSuspects(fGroupsScr, suspects, adduct = "[M+H]+", onlyHits = TRUE, amend = TRUE)

# step 4
mslistsScr <- generateMSPeakLists(fGroupsScr, "mzr")
convertToMFDB(TPs, "TP-database.csv", includeParents = TRUE)
compoundsScr <- generateCompounds(fGroupsScr, mslistsScr, "metfrag", adduct = "[M+H]+", database = "csv",
                                  extraOpts = list(LocalDatabasePath = "TP-database.csv"))

# step 5a
compoundsScr <- filter(compoundsScr, minExplainedPeaks = 1)

# step 5b
fGroupsScrAnn <- annotateSuspects(fGroupsScr, MSPeakLists = mslistsScr,
                                  compounds = compoundsScr)
fGroupsScrAnn <- filter(fGroupsScrAnn, maxLevel = 3, onlyHits = TRUE)

# step 6
componTP <- generateComponents(fGroupsScrAnn, "tp", TPs = TPs, MSPeakLists = mslistsScr,
                               compounds = compoundsScr)
fGroupsScrAnn <- fGroupsScrAnn[results = componTP]

# step 7
reportHTML(fGroupsScrAnn, MSPeakLists = mslistsScr, compounds = compoundsScr,
           components = componTP)
```

### Non-target screening of predicted TPs {#TPsEx3}

This example uses metabolic logic to calculate possible TPs for all feature groups from a complete non-target screening. This example demonstrates how a workflow can be performed when little is known about the identity of the parents. The steps of this workflow are:

1. Formula annotations are performed for all feature groups.
2. These results are then limited to the top 5 candidates, and only feature groups with annotations are kept.
3. The TPs are calculated for all remaining feature groups.
4. A suspect screening is performed to find the TPs. Unlike the previous example feature groups without hits are kept ([discussed here](#TPsFGroups)).
5. The components are generated
6. The components are filtered:
    a. The TPs must follow an expected [retention time direction](#TPsProc)
    b. The parent/TPs should have at least one candidate formula that fits with the transformation.
7. Only feature groups are kept with parent/TP assignments and all results are reported.

```{r eval=FALSE}
# steps 1-2
mslists <- generateMSPeakLists(fGroups, "mzr")
formulas <- generateFormulas(fGroups, mslists, "genform", adduct = "[M+H]+")
formulas <- filter(formulas, topMost = 5)
fGroups <- fGroups[results = formulas]

# step 3
TPs <- generateTPs("logic", fGroups = fGroups, adduct = "[M+H]+")

# step 4
suspects <- convertToSuspects(TPs)
fGroupsScr <- screenSuspects(fGroups, suspects, adduct = "[M+H]+", onlyHits = FALSE)

# step 5
componTP <- generateComponents(fGroupsScr, "tp", TPs = TPs, MSPeakLists = mslists, formulas = formulas)

# step 6
componTP <- filter(componTP, retDirMatch = TRUE, formulas = formulas)

# step 7
fGroupsScr <- fGroupsScr[results = componTP]
reportHTML(fGroupsScr, MSPeakLists = mslists, formulas = formulas, components = componTP)
```

### Non-target screening of TPs by annotation similarities {#TPsEx4}

This example shows a workflow where no TP data from a prediction or library is used. Instead, this workflow relies on statistics and MS/MS data to find feature groups which may potentially have a parent - TP relationship. The workflow is similar to that of the previous example. The steps of this workflow are:

1. [Fold changes](#FCCalc) (FC) between two sample groups are calculated to classify which feature groups are decreasing (i.e. parents) or increasing (i.e. TPs).
2. Feature groups without classification are removed.
3. Formula annotations are performed like the previous example.
4. The componentization is performed and the FC classifications are used to specify which feature groups are to be considered parents or TPs.
5. Only TPs are kept that show a high MS/MS spectral similarity and share at least one fragment with their parent.
6. Only feature groups are kept with parent/TP assignments and all results are reported.

```{r eval=FALSE}
# step 1
tab <- as.data.table(fGroups, FCParams = getFCParams(c("before", "after")))
groupsParents <- tab[classification == "decrease"]$group
groupsTPs <- tab[classification == "increase"]$group

# step 2
fGroups <- fGroups[, union(groupsParents, groupsTPs)]

# step 3
mslists <- generateMSPeakLists(fGroups, "mzr")
formulas <- generateFormulas(fGroups, mslists, "genform", adduct = "[M+H]+")
formulas <- filter(formulas, topMost = 5)
fGroups <- fGroups[results = formulas]

# step 4
componTP <- generateComponents(algorithm = "tp",
                               fGroups = fGroups[, groupsParents],
                               fGroupsTPs = fGroups[, groupsTPs],
                               MSPeakLists = mslists, formulas = formulas)

# step 5
componTP <- filter(componTP, minSpecSimBoth = 0.75, minFragMatches = 1)

# step 6
fGroups <- fGroups[results = componTP]
reportHTML(fGroups, MSPeakLists = mslists, formulas = formulas, components = componTP)
```
# Introduction

Nowadays there are various software tools available to process data from non-target analysis (NTA) experiments. Individual tools such as [ProteoWizard], [XCMS], [OpenMS], [MetFrag] and mass spectrometry vendor tools are often combined to perform a complete data processing workflow. During this workflow, raw data files may undergo pre-treatment (e.g. conversion), chromatographic and mass spectral data are combined to extract so called _features_ (or 'peaks') and finally annotation is performed to elucidate chemical identities. The aim of `patRoon` is to harmonize the many available tools in order to provide a consistent user interface without the need to know all the details of each individual software tool and remove the need for tedious conversion of data when multiple tools are used. The name is derived from a Dutch word that means _pattern_ and may also be an acronym for _hyPhenated mAss specTROmetry nOn-target aNalysis_. The workflow of non-target analysis is typically highly dependent on several factors such as the analytical instrumentation used and requirements of the study. For this reason, `patRoon` does not enforce a certain workflow. Instead, most workflow steps are optional, are highly configurable and algorithms can easily be mixed or even combined. Furthermore, `patRoon` supplies a straightforward interface to easily inspect, select, visualize and report all data that is generated during the workflow.

The documentation of `patRoon` consists of three parts:

1. A tutorial (accessible at [here][tutorial])
2. This handbook
3. The reference manual (accessible in `R` with ``?`patRoon-package` `` or online [here][reference])

New users are highly recommended to start with the tutorial: this document provides an interactive introduction in performing a basic NTA processing workflow with `patRoon`. The handbook provides a more thorough overview of all concepts, functionalities and provides instructions and many examples on working with `patRoon`. Finally, the reference manual provides all the gritty details for all functionalities, and is meant if you want to know more details or need a quick reminder how a function should be used. 



```{r child=file.path(vignDir, "shared", "_refs.Rmd")}
```
[patroonData]: https://github.com/rickhelmus/patRoonData
[readme]: https://rickhelmus.github.io/patRoon/
[ProteoWizard]: http://proteowizard.sourceforge.net/
[OpenMS]: http://openms.de/
[mzR]: https://bioconductor.org/packages/release/bioc/html/mzR.html
[GenForm]: https://sourceforge.net/projects/genform/
[MetFrag]: http://ipb-halle.github.io/MetFrag/
[MetFragCL]: http://ipb-halle.github.io/MetFrag/projects/metfragcl/
[PubChem]: https://pubchem.ncbi.nlm.nih.gov/
[ChemSpider]: http://www.chemspider.com/
[CompTox]: https://comptox.epa.gov/dashboard
[FOR-IDENT]: https://water.for-ident.org
[XCMS]: https://www.bioconductor.org/packages/release/bioc/html/xcms.html
[XCMS3]: https://www.bioconductor.org/packages/release/bioc/html/xcms.html
[enviPick]: https://cran.r-project.org/web/packages/enviPick/index.html
[CompTox-dl]: ftp://newftp.epa.gov/COMPTOX/Sustainable_Chemistry_Data/Chemistry_Dashboard/MetFrag_metadata_files
[CompTox-smoke]: https://zenodo.org/record/3364464#.XnjM-XLvKUk
[CompTox-WW]: https://zenodo.org/record/3472781#.XnjMAHLvKUk
[PCLite-dl]: https://zenodo.org/record/4432124/
[PCLite-paper]: https://doi.org/10.1186/s13321-021-00489-0
[csvDB-ex]: https://raw.githubusercontent.com/rickhelmus/patRoon/master/tests/testthat/test_data/test-mf-db-isomers.csv
[SIRIUS]: https://bio.informatik.uni-jena.de/software/sirius/
[CSI:FingerID]: https://www.csi-fingerid.uni-jena.de/
[HMDB]: http://www.hmdb.ca/
[KEGG]: https://www.genome.jp/kegg/
[RAMClustR]: https://github.com/sneumann/RAMClustR
[CAMERA]: http://msbi.ipb-halle.de/msbi/CAMERA/
[nontarget]: https://cran.r-project.org/web/packages/nontarget/index.html
[dynamicTreeCut]: https://cran.r-project.org/package=dynamicTreeCut
[rcdk]: https://github.com/CDK-R/cdkr
[data.table]: https://github.com/Rdatatable/data.table/wiki
[sqlite]: https://www.sqlite.org/index.html
[msconvert]: http://proteowizard.sourceforge.net/tools/msconvert.html
[FeatureFinderMetabo]: https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureFinderMetabo.html
[FeatureLinkerUnlabeled]: https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureLinkerUnlabeled.html
[MapAlignerPoseClustering]: https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_MapAlignerPoseClustering.html
[UpSet]: https://caleydo.org/tools/upset/
[tutorial]: https://rickhelmus.github.io/patRoon/articles/tutorial.html
[reference]: https://rickhelmus.github.io/patRoon/reference/index.html
[miniCRAN]: https://github.com/andrie/miniCRAN
[Rtools]: https://cran.r-project.org/bin/windows/Rtools/
[JavaJDK]: https://www.oracle.com/technetwork/java/javase/downloads/index.html
[pngquant]: https://pngquant.org/
[OpenBabel]: http://openbabel.org/wiki/Main_Page
[Docker]: https://www.docker.com/
[RStudio]: https://rstudio.com/
[issues]: https://github.com/rickhelmus/patRoon/issues
[RStudio image]: https://hub.docker.com/r/rocker/rstudio
[future]: https://github.com/HenrikBengtsson/future
[future.apply]: https://github.com/HenrikBengtsson/future.batchtools
[future.batchtools]: https://github.com/HenrikBengtsson/future.batchtools
[withr]: https://cran.r-project.org/web/packages/withr/index.html
[processx]: https://github.com/r-lib/processx
[progressr]: https://github.com/HenrikBengtsson/progressr
[cliqueMS]: https://github.com/osenan/cliqueMS
[MetaboliteAdductDecharger]: https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/UTILS_MetaboliteAdductDecharger.html
[BioTransformer]: https://bitbucket.org/djoumbou/biotransformer/src/master/
[KPIC2]: https://github.com/hcji/KPIC2
[SAFD]: https://bitbucket.org/SSamanipour/safd.jl/src/master/
[YAML]: https://en.wikipedia.org/wiki/YAML
[MetaClean]: https://github.com/KelseyChetnik/MetaClean/
[PubChemLiteTR]: https://doi.org/10.5281/zenodo.5644560
[Handbook]: https://rickhelmus.github.io/patRoon/handbook_bd/index.html]
## {{ .annotationClass .compscl_ann-{ grpi } }}

### Dendrogram { grp } {{ .annotationClass .compscl_ann-{ grpi } }}

![]({ dendro })

## {{ .annotationClass .compscl_ann-{ grpi } }}

### {{ .annotationClass .compscl_ann-{ grpi } }}

{ mcs }

```{r echo=FALSE}
cTable <- componentTable(rmdVars$components)
cInfo <- componentInfo(rmdVars$components)
cNames <- names(rmdVars$components)

# the given fGroups may be a subset: make sure to only report components with
# given fGroups.
# NOTE: we cannot report a subset of the components object as it removes
# necessary metadata.
subComps <- rmdVars$components[, names(rmdVars$fGroups)]
indsWithFGroups <- which(names(rmdVars$components) %in% names(subComps))

message("Plotting components...")
prog <- openProgBar(0, length(indsWithFGroups))
allPlots <- vector("character", length(rmdVars$components) * 4)
curPlotInd <- 0
plotPathFull <- getPlotPath(FALSE)
plotPathLink <- getPlotPath(TRUE)

# HACK: this should be replaced some proper inheritance/methods at some point
isHClust <- inherits(rmdVars$components, "componentsIntClust")
isHomol <- inherits(rmdVars$components, "componentsNT")
isSet <- inherits(rmdVars$components, "componentsSet")

if (isHClust)
    clProps <- clusterProperties(rmdVars$components)

for (ci in indsWithFGroups)
{
    curPlotInd <- curPlotInd + 1
    allPlots[curPlotInd] <- file.path(plotPathFull, sprintf("component_spec_%d.png", ci))
    makeCachedPlot(allPlots[curPlotInd], "plotSpectrum",
                   list(rmdVars$components, ci,
                        main = sprintf("ret: %.1f; m/z: %.4f - %.4f", cInfo$ret[ci], min(cTable[[ci]]$mz), max(cTable[[ci]]$mz))),
                        7, 4.5, bg = NA, cacheDB = rmdVars$cacheDB)

    curPlotInd <- curPlotInd + 1
    allPlots[curPlotInd] <- file.path(plotPathFull, sprintf("component_eic_%d.png", ci))
    makeCachedPlot(allPlots[curPlotInd], "plotChroms",
                   list(rmdVars$components, ci, rmdVars$fGroups, rtWindow = rmdVars$EICRtWindow,
                        mzExpWindow = rmdVars$EICMzExpWindow, retMin = rmdVars$retMin, EICs = rmdVars$EICs),
                   7, 4.5, bg = NA, cacheDB = rmdVars$cacheDB)
    
    if (isHClust)
    {
        curPlotInd <- curPlotInd + 1
        allPlots[curPlotInd] <- file.path(plotPathFull, sprintf("component_int_norm_%d.png", ci))
        makeCachedPlot(allPlots[curPlotInd], "plotInt",
                       list(rmdVars$components, index = ci, main = "normalized"),
                       3.3, 3.3, bg = NA, cacheDB = rmdVars$cacheDB)
        
        curPlotInd <- curPlotInd + 1
        allPlots[curPlotInd] <- file.path(plotPathFull, sprintf("component_int_abs_%d.png", ci))
        fg <- fGroups[, unique(cTable[[ci]]$group)]
        makeCachedPlot(allPlots[curPlotInd], "plotInt", list(fg, average = clProps$average, main = "absolute"),
                       3.3, 3.3, bg = NA, cacheDB = rmdVars$cacheDB)
    }
    
    setTxtProgressBar(prog, ci)
}

setTxtProgressBar(prog, length(indsWithFGroups))
close(prog)

if (rmdVars$optimizePng && curPlotInd > 0)
    optimizePngPlots(allPlots[seq_len(curPlotInd)])
```


Components {data-orientation=rows}
===

```{r echo=FALSE,eval=isHClust}
rmdText <- knitr::knit(text = glue::glue("
## {{data-height=350}}

### heatmap

{ ticks } {{r fig.width=6, fig.height=5}}
plotHeatMap(rmdVars$components, interactive = { inter })
{ ticks }

### dendrogram

{ ticks } {{r fig.width=6, fig.height=5}}
plot(rmdVars$components)
{ ticks }

", ticks = "```", inter = as.character(rmdVars$interactiveHeat)))
```

```{r echo=FALSE,eval=isHomol}
rmdText <- knitr::knit(text = glue::glue("
##

### Linked series

{ ticks } {{r}}
plotGraph(rmdVars$components, onlyLinked = TRUE)
{ ticks }

", ticks = "```"))
```

`r if (isHClust || isHomol) rmdText`

<style> .components { overflow-x: auto; } </style>

```{r echo=FALSE}
makeCompDT <- function(s, scrollY)
{
    if (is.null(s))
        compInds <- indsWithFGroups
    else
        compInds <- intersect(cInfo[set == s, which = TRUE], indsWithFGroups)
    
    sppaths <- file.path(plotPathLink, sprintf("component_spec_%d.png", compInds))
    eicpaths <- file.path(plotPathLink, sprintf("component_eic_%d.png", compInds))

    if (rmdVars$selfContained)
    {
        sppaths <- sapply(sppaths, knitr::image_uri)
        eicpaths <- sapply(eicpaths, knitr::image_uri)
    }
    
    # clearout useless columns with only NA in them
    cTable <- sapply(cTable, function(ct)
    {
        ct[, sapply(ct, function(x) !all(is.na(x))), with = FALSE]
    }, simplify = FALSE)
    infoTables <- sapply(compInds, function(compi) knitr::kable(cTable[[compi]], "html") %>%
                             kableExtra::kable_styling(font_size = 11) %>%
                             kableExtra::scroll_box(extra_css = "overflow: auto; width: 350px; height: 300px;"))
    
    compTable <- data.table(component = names(rmdVars$components)[compInds],
                            info = infoTables,
                            EIC = imgTags(eicpaths))
    
    if (isHClust)
    {
        intnpaths <- file.path(plotPathLink, sprintf("component_int_norm_%d.png", compInds))
        intapaths <- file.path(plotPathLink, sprintf("component_int_abs_%d.png", compInds))
        
        if (rmdVars$selfContained)
        {
            intnpaths <- sapply(intnpaths, knitr::image_uri)
            intapaths <- sapply(intapaths, knitr::image_uri)
        }
        
        compTable[, intensities := paste0(imgTags(intnpaths), "<br>", imgTags(intapaths))]
    } else
        compTable[, spectrum := imgTags(sppaths)]
    
    eId <- if (!is.null(s)) paste0("componentsTable_", s) else "componentsTable"
    initDT <- DT::JS("function(settings, json) {",
                     "setTimeout(function() {",
                     sprintf("$(\"#%s .dataTable\").DataTable().columns.adjust().draw(); }", eId),
                     ", 25); }")
    DT::datatable(compTable, options = list(scrollX = TRUE, scrollY = paste0(scrollY, "px"), deferRender = TRUE,
                                            dom = "lrtip", pageLength = 25, autoWidth = FALSE,
                                            initComplete = initDT, ordering = FALSE),
                  class = "striped row-border", elementId = eId,
                  fillContainer = TRUE, rownames = FALSE, escape = FALSE)
}

compTableHeight <- if (isHClust || isHomol) 800 else 1200
hd <- function(cl) sprintf("## { data-height=%d .components %s }\n\n", compTableHeight, cl)
subhd <- function(t) paste0(sprintf("### %s { .components }\n\n", t),
                            "NOTE: only components with feature group data are shown here.\n\n")
body <- function(s) glue::glue("\n{ ticks } {{r}}\nmakeCompDT({ s }, { h })\n{ ticks }\n\n",
                               ticks = "```", s = s, h = compTableHeight - 200)

if (isSet)
{
    rmdText <- hd(".tabset")
    for (s in sets(rmdVars$components))
        rmdText <- c(rmdText, paste0(subhd(s), body(paste0("\"", s, "\""))))
    
} else
    rmdText <- paste0(hd(""), subhd("Components"), body("NULL"))

rmdText <- if (length(rmdText) > 0) paste0(knitr::knit(text = rmdText), collapse = "\n")
```

`r if (nzchar(rmdText)) rmdText`
```{r echo=FALSE, results="hide"}
formulas <- rmdVars$formulas; compounds <- rmdVars$compounds
if (!is.null(formulas) && !is.null(rmdVars$formulasTopMost))
    formulas <- filter(formulas, topMost = rmdVars$formulasTopMost)
if (!is.null(compounds) && !is.null(rmdVars$compoundsTopMost))
    compounds <- filter(compounds, topMost = rmdVars$compoundsTopMost)
```

<script>`r readAllFile(system.file("js", "utils-report.js", package = "patRoon"))`</script>

Annotation {data-orientation=rows}
===

##

### feature groups { data-width=300 .fGroups }

<style> .fGroups { overflow-x: auto; } </style>

```{r echo=FALSE}
compFGroups <- compsClustFGroups <- formFGroups <- componentFGroups <- mfWebFGroups <- character()
if (!is.null(compounds))
    compFGroups <- intersect(groupNames(compounds), rmdVars$groupNames)
if (!is.null(rmdVars$compsCluster))
{
    cl <- clusters(rmdVars$compsCluster)
    compsClustFGroups <- rmdVars$groupNames[rmdVars$groupNames %in% names(cl)]
}
if (rmdVars$includeMFWebLinks != "none")
{
    if (rmdVars$includeMFWebLinks == "compounds")
        mfWebFGroups <- compFGroups
    else if (rmdVars$includeMFWebLinks == "MSMS")
        mfWebFGroups <- rmdVars$groupNames[sapply(rmdVars$groupNames,
                                                  function(grp) any(sapply(peakLists(rmdVars$MSPeakLists),
                                                                           function(pa) !is.null(pa[[grp]]) &&
                                                                               !is.null(pa[[grp]][["MSMS"]]))),
                                                  USE.NAMES = FALSE)]
    mfWebLinks <- sapply(mfWebFGroups, function(grp)
    {
        if (grp %in% compFGroups && is(compounds, "compoundsMF"))
        {
            # make link with used MF settings
            set <- settings(compounds)
        }
        else
            set <- NULL
        
        return(buildMFLandingURL(set, rmdVars$MSPeakLists[[grp]][["MSMS"]],
                                 rmdVars$gInfo[grp, "mzs"]))
    })
}
if (!is.null(formulas) && "formulas" %in% rmdVars$reportPlots)
    formFGroups <- intersect(groupNames(formulas), rmdVars$groupNames)
if (!is.null(rmdVars$components))
{
    cTable <- componentTable(rmdVars$components)
    componentFGroups <- unique(unlist(sapply(cTable, "[[", "group")))
    componentFGroups <- componentFGroups[componentFGroups %in% rmdVars$groupNames]
}
plotGroups <- unique(c(compFGroups, compsClustFGroups, formFGroups))
allGroups <- unique(c(plotGroups, componentFGroups, mfWebFGroups))

# UNDONE: replace by proper inheritance
isSusp <- isScreening(rmdVars$fGroups)
keepSuspCols <- character()
colsRound5 <- "mz"
if (isSusp)
{
    fGroupsSpecDT <- as.data.table(rmdVars$fGroups, collapseSuspects = NULL)
    
    isSet <- isFGSet(rmdVars$fGroups)
    keepSuspCols <- c("susp_name", "susp_d_rt", "susp_d_mz")

    clR2 <- c("susp_d_rt", "susp_annSimForm", "susp_annSimComp", "susp_annSimBoth", "susp_maxFragMatchesRel")
    clR2 <- getAllSuspCols(clR2, names(fGroupsSpecDT), mergedConsensusNames(rmdVars$fGroups))
    
    if (length(clR2) > 0)
        fGroupsSpecDT[, (clR2) := lapply(mget(clR2), round, 2)]

    colsRound5 <- c(colsRound5, "susp_d_mz")

    # NOTE: all cols must be same type
    mergeCols <- function(fg, curColNames, parColNames, mergedName, is = isSet)
    {
        if (is)
        {
            for (s in sets(rmdVars$fGroups))
                fg <- mergeCols(fg, paste0(curColNames, "-", s), parColNames, paste0(mergedName, "-", s), FALSE)
            return(fg)
        }
        
        present <- which(curColNames %in% names(fg))
        if (length(present) == 0)
            return(fg)
        
        newColName <- paste(mergedName, paste0("(", paste0(parColNames[present], collapse = "/"), ")"))
        keepSuspCols <<- c(keepSuspCols, newColName)
        
        fg[, (newColName) := do.call(paste, c(mget(curColNames[present]), list(sep = " / ")))]
        return(fg)
    }

    fGroupsSpecDT <- mergeCols(fGroupsSpecDT, c("susp_formRank", "susp_compRank"), c("form", "comp"), "rank")
    fGroupsSpecDT <- mergeCols(fGroupsSpecDT, c("susp_annSimForm", "susp_annSimComp", "susp_annSimBoth"),
                               c("form", "comp", "both"), "annotated sim")
    fGroupsSpecDT <- mergeCols(fGroupsSpecDT, c("susp_maxFrags", "susp_maxFragMatches", "susp_maxFragMatchesRel"),
                               c("suspect", "max matches", "max matches rel"), "fragments")
    
    keepSuspCols <- c(keepSuspCols, getAllSuspCols("susp_estIDLevel", names(fGroupsSpecDT),
                                                   mergedConsensusNames(rmdVars$fGroups)))
    
    if (isSet)
        keepSuspCols <- c(keepSuspCols, "susp_sets")
} else
    fGroupsSpecDT <- as.data.table(rmdVars$fGroups)

fGroupsSpecDT <- fGroupsSpecDT[group %in% allGroups]
if (rmdVars$retMin)
    fGroupsSpecDT[, ret := ret / 60]

if (!is.null(fGroupsSpecDT[["neutralMass"]]))
    colsRound5 <- c(colsRound5, "neutralMass")

fGroupsSpecDT[, ret := round(ret, 2)]
fGroupsSpecDT[, (colsRound5) := lapply(mget(colsRound5), round, 5)]

fGroupsSpecDT <- fGroupsSpecDT[, intersect(c("group", "ret", "mz", "adduct", "neutralMass", keepSuspCols),
                                           names(fGroupsSpecDT)), with = FALSE]
fGroupsSpecDT[, groupInd := match(group, names(rmdVars$fGroups))]

showButton <- function(title, jsFunc, ...)
{
    args <- paste0("'", unlist(list(...)), "'", collapse = ", ")
    sprintf("<button onclick=\"%s(%s);\" style=\"padding: 0px 3px 0px 3px\">%s</button>", jsFunc, args, title)
}
maybeAddButton <- function(g, subGroups, ...) if (g %in% subGroups) showButton(...) else "^"

compButtons <- sapply(fGroupsSpecDT$group, function(g) maybeAddButton(g, compFGroups, "compounds", "showAnnotation",
                                                                      match(g, names(rmdVars$fGroups)), "compounds"))
compCLButtons <- sapply(fGroupsSpecDT$group, function(g) maybeAddButton(g, compsClustFGroups,
                                                                        "compounds clust", "showCompoundsCluster",
                                                                        match(g, names(rmdVars$fGroups))))
formButtons <- sapply(fGroupsSpecDT$group, function(g) maybeAddButton(g, formFGroups, "formulas",
                                                                      "showAnnotation", match(g, names(rmdVars$fGroups)),
                                                                      "formulas"))
mfWebButtons <- sapply(fGroupsSpecDT$group, function(g) maybeAddButton(g, mfWebFGroups, "MetFrag web",
                                                                       "window.open", mfWebLinks[g],
                                                                       "_blank"))

sp <- paste0(rep("&nbsp;", 4), collapse = "")
buttons <- paste(compButtons, compCLButtons, formButtons, mfWebButtons, sep = sp)
buttons <- gsub("\\^(&nbsp;)*", "", buttons) # remove placeholder (^) with accompanying spaces
fGroupsSpecDT[, show := buttons]
setcolorder(fGroupsSpecDT, c("groupInd", "group", "ret", "mz", "show"))

if (!is.null(rmdVars$components))
{
    annCols <- c("isogroup", "isonr", "charge", "ppm", # RAMClustR
                 "isotopes", "adnr", "adduct_rule", "adduct_charge", "adduct_nmol", "M_adduct", # CAMERA
                 "adduct_ion", # RC/CAMERA
                 "hsnr", "rGroup", # nontarget
                 "TP_name") # TPs
    
    fGroupsSpecDT[, components := sapply(group, function(grp)
    {
        cmps <- findFGroup(rmdVars$components, grp)
        if (length(cmps) == 0)
            return("")
        
        return(wrapStr(paste0(sapply(cmps, function(cmpi)
        {
            cline <- cTable[[cmpi]][group == grp]
            if (nrow(cline) > 1) # some components like NT/TP may have >1 row per fGroup
            {
                if (!is.null(cline[["rGroup"]])) # NT
                    cline[, rGroup := paste0(rGroup, collapse = "/")]
                cline <- cline[1]
            }
            cline <- cline[, sapply(cline, function(x) !is.na(x) && nzchar(x),
                                    USE.NAMES = FALSE), with = FALSE]
            annColsPresent <- annCols[annCols %in% names(cline)]
            cname <- names(rmdVars$components)[cmpi]
            
            if (length(annColsPresent) > 0)
            {
                cline <- cline[, annColsPresent, with = FALSE]
                for (j in seq_along(cline))
                {
                    if (is.numeric(cline[[j]]))
                        set(cline, 1L, j, round(cline[[j]], 5))
                }
                ann <- paste0(sprintf("%s: %s", names(cline), cline), collapse = ", ")
                return(sprintf("%s (%s)", cname, ann))
            }
            return(cname)
        }), collapse = ", "), width = 50, sep = "<br>"))
    }, USE.NAMES = FALSE)]
}

dtOpts <- list(paging = FALSE, pageLength = -1, scrollX = TRUE, scrollY = "200px",
               dom = "frtip",
               initComplete = DT::JS("function(settings, json)",
                                     "{ setTimeout(initAnnotation, 25); }"),
               order = list(list(1, "asc")),
               columnDefs = list(list(visible = FALSE, targets = 0),
                                 list(className = "dt-nowrap", targets = 4),
                                 list(className = "dt-center",
                                      targets = (seq_len(ncol(fGroupsSpecDT))[-5])-1)))

selCols <- ncol(fGroupsSpecDT) > 5 # selectable columns if >5 columns
if (selCols)
{
    dtOpts <- c(dtOpts, list(buttons = list(list(extend = "colvis", background = FALSE,
                                                 columns = seq(5, ncol(fGroupsSpecDT)-1),
                                                 collectionLayout = "three-column"))))
    dtOpts$dom <- paste0("B", dtOpts$dom)
}

dtArgs <- list(fGroupsSpecDT, escape = FALSE, rownames = FALSE, elementId = "fGroupsTable",
               options = dtOpts)
if (selCols)
    dtArgs <- c(dtArgs, list(extensions = "Buttons"))

do.call(DT::datatable, dtArgs)
```

### EIC { data-width=100 }

<img id=EICAnn style="display:none;"></img>
<div id=noAnnotationSelected>Push a **show** button to view annotation data for a feature group.</div>


```{r allAnnPlots, fig.keep='none', eval=length(plotGroups) > 0}
# Generate all plots in advance, since having many code chunks will cause a lot of overhead.

message("Generating spectra...")
prog <- openProgBar(0, length(plotGroups))

plotPathFull <- getPlotPath(FALSE)
plotPathLink <- getPlotPath(TRUE)

allPlots <- setNames(lapply(seq_along(plotGroups), function(i)
{
    grp <- plotGroups[i]
    grpi <- match(grp, names(rmdVars$fGroups))
    grpPlots <- list()
    
    if (grp %in% compFGroups)
    {
        cTable <- compounds[[grp]]
        compsSeq <- seq_len(nrow(cTable))
        grpPlots[["compoundScores"]] <- sapply(compsSeq, function(compi)
        {
            ret <- file.path(plotPathFull, sprintf("compscore_%d_%d.png", grpi, compi))
            makeCachedPlot(ret, "plotScores", list(compounds, compi, grp, rmdVars$compoundsNormalizeScores,
                                                   rmdVars$compoundsExclNormScores, rmdVars$compoundsOnlyUsedScorings),
                           4.5, 4.5, bg = NA, cacheDB = rmdVars$cacheDB)
            return(ret)
        })
        
        grpPlots[["compoundSpectra"]] <- sapply(compsSeq, function(compi)
        {
            ret <- file.path(plotPathFull, sprintf("compspec_%d_%d.png", grpi, compi))
            makeCachedPlot(ret, "plotSpectrum", list(compounds, compi,  grp, rmdVars$MSPeakLists,
                                                     formulas, FALSE),
                           7, 4.5, bg = NA, cacheDB = rmdVars$cacheDB)
            return(ret)
        })
        
        grpPlots[["compoundStructs"]] <- sapply(compsSeq, function(compi)
        {
            ret <- file.path(plotPathFull, sprintf("compstruct_%d_%d.png", grpi, compi))
            makeCachedPlot(ret, "plotStructure", list(compounds, compi, grp, width = 150, height = 150),
                           3, 3, bg = NA, cacheDB = rmdVars$cacheDB)
            return(ret)
        })
    }
    
    if (grp %in% compsClustFGroups)
    {
        plotf <- file.path(plotPathFull, sprintf("dendro_%d.png", grpi))
        makeCachedPlot(plotf, "plot", list(rmdVars$compsCluster, groupName = grp), 8, 4.5, cacheDB = rmdVars$cacheDB)
        grpPlots[["compClustDendro"]] <- plotf
        
        ct <- cutClusters(rmdVars$compsCluster)[[grp]]
        grpPlots[["compClustMCS"]] <- sapply(seq_along(unique(ct)), function(cli)
        {
            ret <- file.path(plotPathFull, sprintf("mcs_%d_%d.png", grpi, cli))
            makeCachedPlot(ret, "plotStructure", list(rmdVars$compsCluster, grp, cli, 100, 100),
                           3, 3, cacheDB = rmdVars$cacheDB)
            return(ret)
        })
    }
    
    if (grp %in% formFGroups)
    {
        fTable <- formulas[[grp]]
        formsSeq <- seq_len(nrow(fTable))

        grpPlots[["formulaScores"]] <- sapply(formsSeq, function(formi)
        {
            ret <- file.path(plotPathFull, sprintf("formscore_%d_%i.png", grpi, formi))
            makeCachedPlot(ret, "plotScores", list(formulas, formi, grp,
                                                   normalizeScores = rmdVars$formulasNormalizeScores,
                                                   excludeNormScores = rmdVars$formulasExclNormScores),
                           4.5, 4.5, bg = NA, cacheDB = rmdVars$cacheDB)
            return(ret)
        })
        
        grpPlots[["formulaSpecs"]] <- sapply(formsSeq, function(formi)
        {
            anPList <- annotatedPeakList(formulas, formi, grp,
                                         MSPeakLists = rmdVars$MSPeakLists, onlyAnnotated = TRUE)
            if (is.null(anPList))
                return("") # No MS/MS data available
            
            ret <- file.path(plotPathFull, sprintf("formspec_%d_%d.png", grpi, formi))
            makeCachedPlot(ret, "plotSpectrum", list(formulas, formi, grp, MSPeakLists = rmdVars$MSPeakLists),
                           6, 4.5, cacheDB = rmdVars$cacheDB)
            return(ret)
        })
    }
    
    setTxtProgressBar(prog, i)
    
    return(grpPlots)
}), plotGroups)
close(prog)

ap <- unlist(allPlots); ap <- ap[nzchar(ap)]
if (rmdVars$optimizePng && length(ap > 0))
    optimizePngPlots(ap)

if (rmdVars$selfContained)
    allPlots <- rapply(allPlots, function(ap) sapply(ap, function(p) if (nzchar(p)) knitr::image_uri(p) else ""), how = "replace")
```


## { .annotationClass .compounds }

### { .annotationClass .compounds }

<style> .compounds { overflow-x: auto; } </style>

```{r echo=FALSE, eval=length(compFGroups) > 0}
compoundsDT <- rbindlist(lapply(compFGroups, function(grp)
{
    ct <- compounds[[grp]]
    
    infoTexts <- sapply(seq_len(nrow(ct)), function(compi)
    {
        it <- paste0(getCompInfoList(ct, compi, mergedConsensusNames(compounds), TRUE), collapse = "<br>")
        if (isSusp)
        {
            # insert suspect names (if any)
            tbl <- as.data.table(rmdVars$fGroups, collapseSuspects = NULL)[group == grp]
            if (!is.null(tbl[["susp_compRank"]]) && any(tbl$susp_compRank == compi, na.rm = TRUE))
                it <- paste(paste("<strong>Suspect(s):</strong>", paste0(tbl[susp_compRank == compi]$susp_name, collapse = ", ")),
                             it, sep = "<br>")
        }
        return(it)
    })
    infoTexts <- makeInfoBox(infoTexts)
    
    fiTables <- sapply(seq_len(nrow(ct)), function(compi)
    {
        apl <- annotatedPeakList(compounds, index = compi, groupName = grp,
                                 MSPeakLists = rmdVars$MSPeakLists, formulas = formulas,
                                 onlyAnnotated = TRUE)
        
        if (is.null(apl) || nrow(apl) == 0)
            return("<div align=\"center\">No annotation available.</div>")
        
        apl[, ion_formula := subscriptFormulaHTML(ion_formula)]
        apl[, neutral_loss := subscriptFormulaHTML(neutral_loss)]
        
        knitr::kable(apl, "html", escape = FALSE) %>%
                           kableExtra::kable_styling(font_size = 11) %>%
                           kableExtra::scroll_box(extra_css = "overflow: auto; height: 125px;")
    })
    
    return(data.table(group = match(grp, names(rmdVars$fGroups)),
                      "#" = seq_len(nrow(ct)),
                      compound = paste0(imgTags(allPlots[[grp]]$compoundStructs), "<br>", infoTexts),
                      spectrum = paste0(imgTags(allPlots[[grp]]$compoundSpectra), "<br>", fiTables),
                      scores = imgTags(allPlots[[grp]]$compoundScores)))
}))

DT::datatable(compoundsDT, options = list(scrollX = TRUE, scrollY = "600px", deferRender = TRUE,
                                          dom = "lrtp", pageLength = 25, autoWidth = TRUE,
                                          ordering = FALSE,
                                          columnDefs = list(list(visible = FALSE, targets = 0))),
              rownames = FALSE, escape = FALSE, elementId = "compoundsTable")
```


## { .annotationClass .formulas }

### { .annotationClass .formulas }

<style> .formulas { overflow-x: auto; } </style>

```{r echo=FALSE, eval=length(formFGroups) > 0}
formulasDT <- rbindlist(lapply(formFGroups, function(grp)
{
    ft <- formulas[[grp]]
    
    infoTexts <- sapply(seq_len(nrow(ft)), function(formi)
    {
        it <- paste0(getFormInfoList(ft, formi, mergedConsensusNames(formulas), TRUE), collapse = "<br>")
        if (isSusp)
        {
            # insert suspect names (if any)
            tbl <- as.data.table(rmdVars$fGroups, collapseSuspects = NULL)[group == grp]
            if (!is.null(tbl[["susp_formRank"]]) && any(tbl$susp_formRank == formi, na.rm = TRUE))
                it <- paste(paste("<strong>Suspect(s):</strong>", paste0(tbl[susp_formRank == formi]$susp_name, collapse = ", ")),
                             it, sep = "<br>")
        }
        return(it)
    })
    
    infoTexts <- makeInfoBox(infoTexts)
    
    fiTables <- sapply(seq_len(nrow(ft)), function(formi)
    {
        apl <- annotatedPeakList(formulas, index = formi, groupName = grp, MSPeakLists = rmdVars$MSPeakLists,
                                 onlyAnnotated = TRUE)
        if (is.null(apl) || nrow(apl) == 0)
            return("<div align=\"center\">No annotation available.</div>")
        
        apl[, ion_formula := subscriptFormulaHTML(ion_formula)]
        apl[, neutral_loss := subscriptFormulaHTML(neutral_loss)]
        
        knitr::kable(apl, "html", escape = FALSE) %>%
            kableExtra::kable_styling(font_size = 11) %>%
            kableExtra::scroll_box(extra_css = "overflow: auto; height: 125px;")
    })
    
    ret <- data.table(group = match(grp, names(rmdVars$fGroups)),
                      neutral_formula = subscriptFormulaHTML(ft$neutral_formula),
                      spectrum = paste0(imgTags(allPlots[[grp]]$formulaSpecs), "<br>", fiTables),
                      scores = paste0(imgTags(allPlots[[grp]]$formulaScores), "<br>", infoTexts))

    return(ret)
}))

DT::datatable(formulasDT, options = list(scrollX = TRUE, scrollY = "600px", deferRender = TRUE,
                                         dom = "lrtp", pageLength = 25, autoWidth = TRUE,
                                         ordering = FALSE,
                                         columnDefs = list(list(visible = FALSE, targets = 0))),
              rownames = FALSE, escape = FALSE, elementId = "formulasTable")

```



```{r echo=FALSE, eval=length(compsClustFGroups) > 0}
rmdTexts <- vector("character", length = length(compsClustFGroups))

message("Generating compounds cluster layout... ", appendLF = FALSE)
# prog <- openProgBar(0, length(plotGroups))
compClustTempl <- readAllFile(system.file("templates", "comp-cluster.Rmd", package = "patRoon"))

cutcl <- cutClusters(rmdVars$compsCluster)
for (i in seq_along(compsClustFGroups))
{
    grp <- compsClustFGroups[i]
    grpi <- match(grp, names(rmdVars$fGroups))
    
    ct <- cutcl[[grp]]
    rmdTexts[i] <-
        paste0(glue::glue(compClustTempl,
                          grpi = grpi,
                          grp = grp,
                          dendro = allPlots[[grp]]$compClustDendro,
                          mcs = paste0(sprintf("![](%s)", allPlots[[grp]]$compClustMCS),
                                       collapse = "\n")),
               collapse = "\n")
}

rmdText <- paste0(rmdTexts, collapse = "\n")
message("Done!")
```

`r if (length(compsClustFGroups) > 0) rmdText`
Feature info {data-orientation=rows}
===

## { .ftInfo }

### Feature info { .ftInfo }

<style> .ftInfo { overflow-x: auto; } </style>

```{r echo=FALSE}
table <- as.data.table(rmdVars$fGroups, qualities = "score", average = TRUE)
table[, EIC := imgTags(chromPaths)]
setcolorder(table, c("group", "EIC"))
if (rmdVars$retMin)
    table[, ret := ret / 60]
for (col in names(table)[(sapply(table, is.numeric))])
    set(table, j = col, value = round(table[[col]], if (col == "mz") 5 else 2))

initDT <- DT::JS("function(settings, json) {",
                 "setTimeout(function() {",
                 "$(\".dataTable\").DataTable().columns.adjust().draw(); }",
                 ", 100); }")
buttonFunc <- function(cols)
{
    cols <- match(cols, names(table))
    DT::JS(paste("function (e, dt, node, config) {",
                 paste0(sprintf("dt.column(%d).visible(!dt.column(%d).visible());", cols, cols), collapse = "\n"),
                 "}"))
}

bts <- list(list(text = "Intensities", action = buttonFunc(replicateGroups(rmdVars$fGroups))))

hiddenCols <- integer()
if (hasFGroupScores(rmdVars$fGroups))
{
    setnames(table, "totalScore", "Total score")
    otherSc <- featureQualityNames(scores = TRUE, totScore = FALSE)
    bts <- c(bts, list(list(text = "Total score", action = buttonFunc("Total score")),
                       list(text = "Other scores", action = buttonFunc(otherSc))))
    hiddenCols <- match(otherSc, names(table))
}

if (isScreening(rmdVars$fGroups))
    bts <- c(bts, list(list(text = "Suspect name", action = buttonFunc("susp_name"))))

tabOpts <- list(dom = "Bfrtip", scrollX = TRUE, scrollY = "600px", deferRender = TRUE,
                paging = FALSE, autoWidth = FALSE,
                initComplete = initDT, order = list(list(0, "asc")),
                buttons = list(list(extend = "colvis", text = "Columns",
                                    background = TRUE, buttons = bts)))
if (length(hiddenCols))
    tabOpts <- c(tabOpts, list(columnDefs = list(list(visible = FALSE, targets = hiddenCols))))

DT::datatable(table, extensions = "Buttons", options = tabOpts, elementId = "ftInfoTable", escape = FALSE)
```
---
title: "report"
author: "`r getPackageName(create = FALSE)`"
date: "`r if (!rmdVars$noDate) date() else '' `"
output:
    flexdashboard::flex_dashboard:
        vertical_layout: scroll
        mathjax: null
        fig_mobile: false
        dev: png
---

```{r setup, include=FALSE}
# knitr::knit_hooks$set(optipng = knitr::hook_pngquant)
knitr::knit_hooks$set(pngquant = function(before, options, envir) suppressMessages(knitr::hook_pngquant(before, options, envir)))
knitr::opts_chunk$set(echo = FALSE, fig.keep = "all", fig.retina = 1, dpi = 72)
if (rmdVars$optimizePng)
    knitr::opts_chunk$set(pngquant = "")

# utility funcs for plotting
getPlotPath <- function(link)
{
    if (rmdVars$selfContained)
        ret <- "."
    else if (link)
        ret <- file.path("report_files", "plots")
    else
        ret <- file.path(rmdVars$outPath, "report_files", "plots")
    mkdirp(ret)
    return(ret)
}

imgTags <- function(img, style = "") 
{
    if (rmdVars$selfContained)
        ret <- sprintf("<img src=%s style='%s'></img>", img, style)
    else
    {
        # return(sprintf("<img src=file://%s></img>", img))
        ret <- sprintf("<img src='%s' style='%s'></img>", img, style)
    }
    return(ifelse(nzchar(img), ret, ""))
}

makeInfoBox <- function(txt)
{
    sprintf("<div style='max-width: 300px; max-height: 432px; border: 1px solid black; border-style: dotted; margin: 1px; padding: 1px; overflow: auto; white-space: nowrap; text-align: left;'>%s</div>", txt)
}

rGroupLenNonEmpty <- length(unique(analysisInfo(removeEmptyAnalyses(rmdVars$fGroups))$group))
rGroupLen <- length(replicateGroups(rmdVars$fGroups))
anyOverlap <- rGroupLen > 1 &&
    length(unique(rmdVars$fGroups, which = replicateGroups(rmdVars$fGroups), outer = TRUE)) < length(rmdVars$fGroups)
if (length(rmdVars$fGroups) > 0 && anyOverlap && rGroupLenNonEmpty > 1)
{
    doPlotChord <- "chord" %in% rmdVars$reportPlots && rGroupLenNonEmpty > 2
    doPlotVenn <- "venn" %in% rmdVars$reportPlots && rGroupLen < 6
    doPlotUpSet <- "upset" %in% rmdVars$reportPlots
} else
    doPlotChord <- doPlotVenn <- doPlotUpSet <- FALSE
doAnnotation <- !is.null(rmdVars$compounds) || !is.null(rmdVars$compsCluster) || !is.null(rmdVars$formulas) ||
    !is.null(rmdVars$components) || inherits(rmdVars$fGroups, "featureGroupsScreening")
doEICs <- length(rmdVars$fGroups) > 0 && "eics" %in% rmdVars$reportPlots
isComponentsTP <- !is.null(rmdVars$components) && inherits(rmdVars$components, "componentsTPs")
rmdText <- NULL
```

<style>
pre, code {
    white-space: pre !important;
    overflow-x: auto !important;
    max-height: 275px;
    overflow-y: auto;
}
</style>

```{r echo=FALSE,eval=doEICs || doAnnotation}
message("Generating chromatograms...")

plotPathFull <- getPlotPath(FALSE)

prog <- openProgBar(0, length(rmdVars$fGroups))
allPlots <- sapply(seq_len(length(rmdVars$fGroups)), function(grpi)
{
    path <- file.path(plotPathFull, sprintf("chrom_%d.png", grpi))
    makeCachedPlot(path, "plotChroms",
                   list(rmdVars$fGroups[, grpi], rmdVars$EICRtWindow,
                        rmdVars$EICMzExpWindow, rmdVars$retMin,
                        rmdVars$EICTopMost, rmdVars$EICTopMostByRGroup, rmdVars$EICs, colourBy = "rGroup"),
                   7, 4.5, bg = NA, cacheDB = rmdVars$cacheDB)
    setTxtProgressBar(prog, grpi)
    return(path)
})
close(prog)

if (rmdVars$optimizePng && length(allPlots) > 0)
    optimizePngPlots(allPlots)

chromPaths <- file.path(getPlotPath(TRUE), sprintf("chrom_%d.png", seq_len(length(rmdVars$fGroups))))
chromPathsFull <- file.path(plotPathFull, sprintf("chrom_%d.png", seq_len(length(rmdVars$fGroups))))
if (rmdVars$selfContained)
    chromPaths <- sapply(chromPaths, knitr::image_uri)

# stuff everything together: https://stackoverflow.com/a/21730473
rmdText <- sprintf("<script>var chromPaths = [ %s ];</script>",
                   paste0("'", chromPaths, "'", collapse = ", "))
```
`r if (!is.null(rmdText)) rmdText`

Summary {data-orientation=rows}
================

## { data-height=350 }

### EICs

```{r obj-plot, fig.width=10, fig.height=4}
par(mai = c(0.9, 0.8, 0.6, 0.1))
plotChroms(rmdVars$fGroups, rmdVars$EICRtWindow, rmdVars$EICMzExpWindow,
           rmdVars$retMin, 1, FALSE, rmdVars$EICs, TRUE, FALSE, colourBy = "fGroups", showLegend = FALSE,
           onlyPresent = rmdVars$EICOnlyPresent)
```


## { data-height=300 }

### Objects

```{r obj-show}
objToShow <- list(rmdVars$fGroups, rmdVars$MSPeakLists, rmdVars$formulas,
                  rmdVars$compounds, rmdVars$components)
objToShow <- objToShow[!sapply(objToShow, is.null)]
for (obji in seq_along(objToShow))
{
    show(objToShow[[obji]])
    cat("\n")
}
```


### Retention vs m/z
```{r fig.height=4}
par(mai = c(0.9, 0.8, 0.1, 0.1))
plot(rmdVars$fGroups, colourBy = "fGroups", showLegend = FALSE, retMin = rmdVars$retMin)
```


`r if (doPlotChord || doPlotVenn || doPlotUpSet) "## { data-height=425 } \n"`

`r if (doPlotChord) "### Chord diagram\n"`
```{r fig.height=5.5, eval=doPlotChord,out.height="400px"}
message("Creating chord diagram... ", appendLF = FALSE)
plotChord(rmdVars$fGroups, average = TRUE)
message("Done!")
```

`r if (doPlotVenn) "### Venn diagram\n"`
```{r fig.height=5.5, eval=doPlotVenn}
plotVenn(rmdVars$fGroups)
```

`r if (doPlotUpSet) "### UpSet diagram\n"`
```{r fig.height=5.5, eval=doPlotUpSet}
plotUpSet(rmdVars$fGroups)
```

`r if (doEICs) "EICs\n===\n"`
```{r results='asis', eval=doEICs}
cat(sprintf("![%s](%s)", names(rmdVars$fGroups), chromPaths), sep = "\n")
```


```{r child="featinfo.Rmd", eval=doEICs}
```


```{r child="components.Rmd", eval=!is.null(rmdVars$components) && length(rmdVars$components) > 0 && !isComponentsTP}
```


```{r child="annotation.Rmd", eval=doAnnotation }
```

```{r child="TPs.Rmd", eval=!is.null(rmdVars$components) && length(rmdVars$components) > 0 && isComponentsTP}
<script>`r readAllFile(system.file("js", "utils-report.js", package = "patRoon"))`</script>

Transformation Products {data-orientation=rows}
===

```{r TPPlots, fig.keep='none'}

# sync objects
components <- rmdVars$components[, intersect(names(rmdVars$fGroups), groupNames(rmdVars$components))]
components <- delete(components, i = !componentInfo(components)$parent_group %chin% names(rmdVars$fGroups))

cInfo <- componentInfo(components)

if (isScreening(rmdVars$fGroups))
{
    suspTbl <- as.data.table(rmdVars$fGroups, collapseSuspects = NULL)
    suspTbl <- suspTbl[, grepl("^susp_|group", names(suspTbl)), with = FALSE]
    setnames(suspTbl, sub("^susp_", "", names(suspTbl)))
}

plotPathFull <- getPlotPath(FALSE)

if (length(components) > 0)
{
    message("Generating parent plots...")
    prog <- openProgBar(0, nrow(cInfo))
    
    makeStructPlot <- function(SMI, out)
    {
        mol <- getMoleculesFromSMILES(SMI, emptyIfFails = TRUE)[[1]]
        withr::with_png(out, width = 1.5, height = 1.5, units = "in", res = 72, bg = NA, code = {
            withr::with_par(list(mar = rep(0, 4)), plot(getRCDKStructurePlot(mol, 150, 150)))
        })
    }
    
    plotIntMainArgs <- list(average = TRUE, normFunc = max)
    if (isFGSet(rmdVars$fGroups))
    {
        plotIntMainArgs <- c(plotIntMainArgs, list(sets = TRUE))
    } else
    {
        plotIntMainArgs <- c(plotIntMainArgs, list(col = "black"))
    }
    
    parentPlots <- setNames(Map(split(cInfo, seq_len(nrow(cInfo))), seq_len(nrow(cInfo)), f = function(parentRow, i)
    {
        grpi <- match(parentRow$parent_group, names(rmdVars$fGroups))
        plots <- list()
        fg <- rmdVars$fGroups[, parentRow$parent_group]
        
        plots[["int"]] <- file.path(plotPathFull, sprintf("int-parent_%d.png", grpi))
        makeCachedPlot(plots[["int"]], "plotInt", c(list(fg), plotIntMainArgs), 3.3, 3.3,
                       bg = NA, cacheDB = rmdVars$cacheDB)
        
        if (!is.null(parentRow[["parent_SMILES"]]))
        {
            plots[["struct"]] <- file.path(plotPathFull, sprintf("struct-parent_%d.png", grpi))
            makeStructPlot(parentRow$parent_SMILES, plots[["struct"]])
        }
        
        setTxtProgressBar(prog, i)
        
        return(plots)
    }), cInfo$name)
    
    message("Generating TP plots...")
    prog <- openProgBar(0, length(components))
    cmpTab <- as.data.table(components)
    # isFGScrAnnotated <- isScreening(rmdVars$fGroups) && screenInfo(rmdVars$fGroups)[[""]]
    
    TPPlotName <- function(cmp, grp) paste0(cmp, "_", grp)
    
    TPPlots <- setNames(Map(split(cmpTab, seq_len(nrow(cmpTab))), seq_len(nrow(cmpTab)), f = function(ctRow, i)
    {
        grpi <- match(ctRow$group, names(rmdVars$fGroups))
        plots <- list()
        fg <- rmdVars$fGroups[, ctRow$group]
        
        plots[["int"]] <- file.path(plotPathFull, sprintf("int-TP_%d.png", grpi))
        makeCachedPlot(plots[["int"]], "plotInt", c(list(fg), plotIntMainArgs), 3.3, 3.3, bg = NA,
                       cacheDB = rmdVars$cacheDB)
        
        SMI <- ctRow[["SMILES"]]
        if (!is.null(SMI))
        {
            plots[["struct"]] <- file.path(plotPathFull, sprintf("struct-parent_%d.png", grpi))
            makeStructPlot(SMI, plots[["struct"]])
        }
        
        if (!is.null(rmdVars[["MSPeakLists"]]))
        {
            # try to plot a mirror spectrum: use compounds if possible, otherwise try formulas or finally peak lists
            plSpecArgs <- list()
            
            if (isScreening(rmdVars$fGroups))
            {
                suspParRow <- suspTbl[name == ctRow$parent_name & group == ctRow$parent_group]
                suspTPRow <- suspTbl[name == ctRow$TP_name & group == ctRow$group]
                if (!is.null(rmdVars[["compounds"]]) && !is.null(suspTbl[["compRank"]]) &&
                    all(c(ctRow$parent_group, ctRow$group) %chin% groupNames(rmdVars$compounds)) &&
                    nrow(suspTPRow) == 1 && !is.na(suspParRow$compRank) && !is.na(suspTPRow$compRank))
                {
                    plSpecArgs <- list(obj = rmdVars$compounds, formulas = rmdVars[["formulas"]],
                                       index = c(suspParRow$compRank, suspTPRow$compRank),
                                       MSPeakLists = rmdVars$MSPeakLists, plotStruct = FALSE)
                }
                else if (!is.null(rmdVars[["formulas"]]) && !is.null(suspTbl[["formRank"]]) &&
                         all(c(ctRow$parent_group, ctRow$group) %chin% groupNames(rmdVars$formulas)) &&
                         nrow(suspTPRow) == 1 && !is.na(suspParRow$formRank) && !is.na(suspTPRow$formRank) &&
                         !is.null(rmdVars$MSPeakLists[[ctRow$parent_group]][["MSMS"]]) &&
                         !is.null(rmdVars$MSPeakLists[[ctRow$group]][["MSMS"]]))
                {
                    plSpecArgs <- list(obj = rmdVars$formulas,
                                       index = c(suspParRow$formRank, suspTPRow$formRank),
                                       MSPeakLists = rmdVars$MSPeakLists)
                }
            }
            
            if (length(plSpecArgs) == 0 && !is.null(rmdVars$MSPeakLists[[ctRow$parent_group]][["MSMS"]]) &&
                !is.null(rmdVars$MSPeakLists[[ctRow$group]][["MSMS"]]))
            {
                # no formulas/compounds, try peak lists
                plSpecArgs <- list(obj = rmdVars$MSPeakLists, MSLevel = 2)
            }
            
            if (length(plSpecArgs) > 0)
            {
                plots[["spec"]] <- file.path(plotPathFull, sprintf("spec-sim_%d.png", grpi))
                makeCachedPlot(plots[["spec"]], "plotSpectrum", c(plSpecArgs, list(groupName = c(ctRow$parent_group, ctRow$group),
                                                                                   specSimParams = rmdVars$specSimParams, title = "")),
                               5, 4.5, bg = NA, cacheDB = rmdVars$cacheDB)
            }
        }
        
        setTxtProgressBar(prog, i)
        
        return(plots)
    }), TPPlotName(cmpTab$name, cmpTab$group))
    
    prepPlots <- function(pl)
    {
        ap <- unlist(allPlots); ap <- ap[nzchar(ap)]
        if (rmdVars$optimizePng && length(ap > 0))
            optimizePngPlots(ap)
        
        if (rmdVars$selfContained)
            pl <- rapply(pl, function(ap) sapply(ap, function(p) if (nzchar(p)) knitr::image_uri(p) else ""), how = "replace")
        return(pl)
    }
    
    parentPlots <- prepPlots(parentPlots); TPPlots <- prepPlots(TPPlots)
}
```


##
    
### Parents { data-width=400 .parents }
    
<style> .parents { overflow-x: auto; } </style>
    
```{r echo=FALSE,eval=length(components)>0}
makeTPTab <- function(cr, cols)
{
    allCols <- unique(unlist(lapply(cols, function(cl) grep(paste0("^", cl), names(cr), value = TRUE))))
    
    roundCols <- function(t) t[, (names(t)) := lapply(.SD, function(x) if (is.double(x)) round(x, 2) else x)]
    
    if (!isTRUE(all.equal(cols, allCols, check.attributes = FALSE)))
    {
        sets <- unique(sub(".+\\-(.+)$", "\\1", allCols[grepl("-", allCols, fixed = TRUE)]))
        
        # NOTE: " " (a space) results in an unnamed column
        
        ret <- rbindlist(lapply(c(" ", sets), function(s)
        {
            takeCols <- if (s == " ") cols else paste0(cols, "-", s)
            whCols <- which(takeCols %in% names(cr))
            if (length(whCols) == 0)
                return(data.table())
            
            t <- cr[, takeCols[whCols], with = FALSE]
            setnames(t, cols[whCols])
            
            t <- roundCols(t)
            
            t[, set := s]
            setcolorder(t, "set")
            
            return(t)
        }), fill = TRUE)
        ret <- transpose(ret, keep.names = " ", make.names = "set")
    }
    else
    {
        ret <- roundCols(cr[, cols, with = FALSE])
        ret <- setnames(transpose(ret, keep.names = " "), 2, "value")
    }
    
    return(knitr::kable(ret, "html", escape = FALSE) %>%
               kableExtra::kable_styling(font_size = 11) %>%
               kableExtra::scroll_box(extra_css = "overflow-x: auto;"))

}

chromPlotStyle <- "width: auto; height: 250px;"

parentsDT <- data.table(compInd = seq_len(nrow(cInfo)))

splitCI <- split(cInfo, seq_len(nrow(cInfo)))

parentInfoCols <- intersect(c("name", "parent_name", "parent_group", "parent_formula", "parent_CID"), names(cInfo))
parentsDT[, parent := sapply(splitCI, function(cir)
{
    parInfo <- cir[, parentInfoCols, with = FALSE]
    setnames(parInfo, sub("^parent_", "", names(parInfo)))
    
    if (!is.null(parInfo[["formula"]]))
        parInfo[, formula := subscriptFormulaHTML(formula)]
    if (!is.null(parInfo[["CID"]]))
        parInfo[, CID := makeDBIdentLink("pubchem", CID)]

    par <- makeInfoBox(paste0(names(parInfo), ": ", parInfo, collapse = "<br>"))
    
    if (!is.null(parentPlots[[cir$name]][["struct"]]))
        par <- paste0(imgTags(parentPlots[[cir$name]][["struct"]]), "<br>", par)
    
    return(par)
})]

parentsDT[, EIC := imgTags(chromPaths[match(cInfo$parent_group, names(rmdVars$fGroups))], style = chromPlotStyle)]

if (isScreening(rmdVars$fGroups))
{
    suspCols <- c("formRank", "compRank", "annSimBoth", "estIDLevel")
    if (any(sapply(suspCols, grepl, names(suspTbl))))
    {
        parentsDT[, screening := sapply(splitCI, function(cir)
        {
            sr <- suspTbl[name == cir$parent_name & group == cir$parent_group]
            if (nrow(sr) > 0)
                makeTPTab(sr, suspCols)
            else
                ""
        })]
    }
}

parentsDT[, "intensity profile" := imgTags(sapply(parentPlots[cInfo$name], "[[", "int"))]

parentsDT[, show := { sprintf("<button onclick=\"showTPs('%s', %d);\" style=\"padding: 0px 3px 0px 3px\">Show</button>",
                              compInd, match(cInfo$parent_group, names(rmdVars$fGroups))) }]
setcolorder(parentsDT, c("compInd", "show"))

DT::datatable(parentsDT,
              extensions = "Buttons",
              options = list(paging = FALSE, pageLength = -1, scrollX = TRUE, scrollY = "300px",
                             dom = "tip",
                             initComplete = DT::JS("function(settings, json)",
                                                   "{ setTimeout(initTPs, 25); }"),
                             order = list(list(0, "asc")),
                             columnDefs = list(list(visible = FALSE, targets = 0),
                                               list(className = "dt-center",
                                                    targets = (seq_len(ncol(parentsDT)))-1))),
              escape = FALSE, rownames = FALSE, elementId = "parentsTable")
```


##

### Transformation Products { .TPsClass }

<style> .TPsClass { overflow-x: auto; } </style>
    
```{r echo=FALSE,eval=length(components)>0}

cTable <- componentTable(components)

TPsDT <- rbindlist(Map(cTable, seq_along(cTable), names(cTable), f = function(cmp, cInd, cName)
{
    ret <- data.table(compInd = cInd, "#" = seq_len(nrow(cmp)))
    
    splitCmp <- split(cmp, seq_len(nrow(cmp)))
    
    trCols <- intersect(c("formula", "retDir", "TP_retDir", "retDiff", "mzDiff", "formulaDiff", "set", "CID",
                          "SMILES"),
                        names(cmp))
    ret[, TP := sapply(splitCmp, function(cr)
    {
        tpInfo <- cr[, c("TP_name", "group", trCols), with = FALSE]
        if (rmdVars$retMin)
            tpInfo[, retDiff := retDiff / 60]
        tpInfo[, c("retDiff", "mzDiff") := .(round(retDiff, 2), round(mzDiff, 5))]
        setnames(tpInfo, "TP_name", "name")
        
        if (!is.null(tpInfo[["retDir"]]) && !is.null(tpInfo[["TP_retDir"]]))
        {
            tpInfo[, "retDir (predicted/actual)" := paste0(TP_retDir, "/", retDir)]
            tpInfo[, c("retDir", "TP_retDir") := NULL]
        }
        if (!is.null(tpInfo[["formula"]]))
            tpInfo[, formula := subscriptFormulaHTML(formula)]
        if (!is.null(tpInfo[["formulaDiff"]]))
            tpInfo[, formulaDiff := subscriptFormulaHTML(formulaDiff)]
        if (!is.null(tpInfo[["CID"]]))
            tpInfo[, CID := makeDBIdentLink("pubchem", CID)]
        
        tp <- makeInfoBox(paste0(names(tpInfo), ": ", tpInfo, collapse = "<br>"))
        
        if (!is.null(TPPlots[[TPPlotName(cName, cr$group)]][["struct"]]))
            tp <- paste0(imgTags(TPPlots[[TPPlotName(cName, cr$group)]][["struct"]]), "<br>", tp)
        
        return(tp)
    })]
    
    ret[, EIC := imgTags(chromPaths[match(cmp$group, names(rmdVars$fGroups))], style = chromPlotStyle)]
    
    if (isScreening(rmdVars$fGroups))
    {
        suspCols <- c("formRank", "compRank", "annSimBoth", "estIDLevel")
        if (any(sapply(suspCols, grepl, names(suspTbl))))
        {
            ret[, screening := sapply(splitCmp, function(cr)
            {
                sr <- suspTbl[name == cr$TP_name & group == cr$group]
                if (nrow(sr) > 0)
                    makeTPTab(sr, suspCols)
                else
                    ""
            })]
        }
    }
    
    if (any(grepl("^(specSimilarity|fragmentMatches|neutralLossMatches)", names(cmp))))
    {
        ret[, similarity := sapply(splitCmp, function(cr)
        {
            simt <- makeTPTab(cr, c("specSimilarity", "specSimilarityPrec", "specSimilarityBoth",
                                    "fragmentMatches", "neutralLossMatches"))
            return(simt)
        })]
    }
    
    ret[, spectrum := sapply(splitCmp, function(cr)
    {
        if (!is.null(TPPlots[[TPPlotName(cName, cr$group)]][["spec"]]))
            return(paste0(imgTags(TPPlots[[TPPlotName(cName, cr$group)]][["spec"]])))
        return("")
    })]
    if (!any(nzchar(ret$spectrum)))
        set(ret, j = "spectrum", value = NULL)

    ret[, "intensity profile" := imgTags(sapply(TPPlots[TPPlotName(cName, cmp$group)], "[[", "int"))]
    return(ret)
}), fill = TRUE) # fill: spectrum can be absent depending on candidate

DT::datatable(TPsDT, options = list(scrollX = TRUE, scrollY = "600px", deferRender = TRUE,
                                    dom = "Blrtp", pageLength = 25, autoWidth = FALSE,
                                    ordering = FALSE,
                                    columnDefs = list(list(visible = FALSE, targets = 0),
                                                      list(className = "dt-center",
                                                           targets = (seq_len(ncol(TPsDT)))-1)),
                                    buttons = list(list(extend = "colvis", background = FALSE,
                                                        columns = seq(3, ncol(TPsDT)-1)))),
              rownames = FALSE, escape = FALSE, elementId = "TPsTable")
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-exported.R
\name{COMStop}
\alias{COMStop}
\title{Internal fix for \pkg{RDCOMClient}, ignore.}
\usage{
COMStop(...)
}
\arguments{
\item{...}{ignore}
}
\description{
Internal fix for \pkg{RDCOMClient}, ignore.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mspeaklists-bruker.R
\name{generateMSPeakListsDAFMF}
\alias{generateMSPeakListsDAFMF}
\alias{generateMSPeakListsDAFMF,featureGroups-method}
\alias{generateMSPeakListsDAFMF,featureGroupsSet-method}
\title{Generate peak lists with Bruker DataAnalysis from bruker features}
\usage{
\S4method{generateMSPeakListsDAFMF}{featureGroups}(
  fGroups,
  minMSIntensity = 500,
  minMSMSIntensity = 500,
  close = TRUE,
  save = close,
  avgFGroupParams = getDefAvgPListParams()
)

\S4method{generateMSPeakListsDAFMF}{featureGroupsSet}(fGroups, ...)
}
\arguments{
\item{fGroups}{The \code{\link{featureGroups}} object from which MS peak lists should be extracted.}

\item{minMSIntensity}{Minimum intensity for peak lists obtained with DataAnalysis. Highly
recommended to set \samp{>0} as DA tends to report many very low intensity peaks.}

\item{minMSMSIntensity}{Minimum intensity for peak lists obtained with DataAnalysis. Highly
recommended to set \samp{>0} as DA tends to report many very low intensity peaks.}

\item{close}{If \code{TRUE} then Bruker files are closed and saved after
processing with DataAnalysis, respectively. Setting \code{close=TRUE}
prevents that many analyses might be opened simultaneously in DataAnalysis,
which otherwise may use excessive memory or become slow. By default
\code{save} is \code{TRUE} when \code{close} is \code{TRUE}, which is
likely what you want as otherwise any processed data is lost.}

\item{save}{If \code{TRUE} then Bruker files are closed and saved after
processing with DataAnalysis, respectively. Setting \code{close=TRUE}
prevents that many analyses might be opened simultaneously in DataAnalysis,
which otherwise may use excessive memory or become slow. By default
\code{save} is \code{TRUE} when \code{close} is \code{TRUE}, which is
likely what you want as otherwise any processed data is lost.}

\item{avgFGroupParams}{A \code{list} with parameters used for averaging of peak lists for feature groups. See
\code{\link{getDefAvgPListParams}} for more details.}

\item{...}{\setsWF Further arguments passed to the non-sets workflow method.}
}
\value{
A \code{\link{MSPeakLists}} object.
}
\description{
Uses 'compounds' that were generated by the Find Molecular Features (FMF) algorithm of Bruker DataAnalysis to extract
MS peak lists.
}
\details{
This function uses Bruker DataAnalysis with FMF to generate MS peak lists. This function is called when calling \code{generateMSPeakLists} with
  \code{algorithm="brukerfmf"}.

This function is similar to \code{\link{generateMSPeakListsDA}}, but uses 'compounds' that were generated by
  the Find Molecular Features (FMF) algorithm to extract MS peak lists. This is generally much faster , however, it
  only works when features were obtained with the \code{\link{findFeaturesBruker}} function. Since all MS spectra
  are generated in advance by Bruker DataAnalysis, only few parameters exist to customize its operation.
}
\note{
If any errors related to \command{DCOM} appear it might be necessary to
  terminate DataAnalysis (note that DataAnalysis might still be running as a
  background process). The \command{ProcessCleaner} application installed
  with DataAnalayis can be used for this.
}
\seealso{
\code{\link{generateMSPeakLists}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{sets-workflow}
\alias{sets-workflow}
\title{Sets workflows}
\description{
With sets workflows in \pkg{patRoon} a complete non-target (or suspect) screening workflow is performed with sample
analyses that were measured with different MS methods (typically positive and negative ionization).
}
\details{
The analyses files that were measured with a different method are grouped in \emph{sets}. In the most typical case,
there is a \code{"positive"} and \code{"negative"} set, for the positively/negatively ionized data, respectively.
However, other distinctions than polarity are also possible (although currently the chromatographic method should be
the same between sets). A sets workflow is typically initiated with the \code{\link{makeSet}} method. The handbook
contains much more details about sets workflows.
}
\seealso{
\code{\link{makeSet}} to initiate sets workflows, \code{\link{workflowStepSet}}, the \verb{Sets workflows}
  sections in other documentation pages and the \pkg{patRoon} handbook.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups-filter.R
\name{feature-filtering}
\alias{feature-filtering}
\alias{filter,featureGroups-method}
\alias{filter,featureGroupsSet-method}
\alias{replicateGroupSubtract,featureGroups-method}
\alias{replicateGroupSubtract}
\title{Filtering of grouped features}
\usage{
\S4method{filter}{featureGroups}(
  obj,
  absMinIntensity = NULL,
  relMinIntensity = NULL,
  preAbsMinIntensity = NULL,
  preRelMinIntensity = NULL,
  absMinAnalyses = NULL,
  relMinAnalyses = NULL,
  absMinReplicates = NULL,
  relMinReplicates = NULL,
  absMinFeatures = NULL,
  relMinFeatures = NULL,
  absMinReplicateAbundance = NULL,
  relMinReplicateAbundance = NULL,
  maxReplicateIntRSD = NULL,
  blankThreshold = NULL,
  retentionRange = NULL,
  mzRange = NULL,
  mzDefectRange = NULL,
  chromWidthRange = NULL,
  featQualityRange = NULL,
  groupQualityRange = NULL,
  rGroups = NULL,
  results = NULL,
  removeBlanks = FALSE,
  checkFeaturesSession = NULL,
  negate = FALSE
)

\S4method{filter}{featureGroupsSet}(
  obj,
  ...,
  negate = FALSE,
  sets = NULL,
  absMinSets = NULL,
  relMinSets = NULL
)

\S4method{replicateGroupSubtract}{featureGroups}(fGroups, rGroups, threshold = 0)
}
\arguments{
\item{absMinIntensity, relMinIntensity}{Minimum absolute/relative intensity for features to be kept. The relative
intensity is determined from the feature with highest intensity (of
all features from all groups). Set to \samp{0} or \code{NULL} to skip this step.}

\item{preAbsMinIntensity, preRelMinIntensity}{As \code{absMinIntensity}/\code{relMinIntensity}, but applied
\emph{before} any other filters. This is typically used to speed-up subsequent filter steps. However, care must be
taken that a sufficiently low value is chosen that is not expected to affect subsequent filtering steps. See below
why this may be important.}

\item{absMinAnalyses, relMinAnalyses}{Feature groups are only kept when they contain data for at least this (absolute
or relative) amount of analyses. Set to \code{NULL} to ignore.}

\item{absMinReplicates, relMinReplicates}{Feature groups are only kept when they contain data for at least this
(absolute or relative) amount of replicates. Set to \code{NULL} to ignore.}

\item{absMinFeatures, relMinFeatures}{Analyses are only kept when they contain at least this (absolute or relative)
amount of features. Set to \code{NULL} to ignore.}

\item{absMinReplicateAbundance, relMinReplicateAbundance}{Minimum absolute/relative abundance that a grouped feature
should be present within a replicate group. If this minimum is not met all features within the replicate group are
removed. Set to \code{NULL} to skip this step.}

\item{maxReplicateIntRSD}{Maximum relative standard deviation (RSD) of intensity values for features within a
replicate group. If the RSD is above this value all features within the replicate group are removed. Set to
\code{NULL} to ignore.}

\item{blankThreshold}{Feature groups that are also present in blank analyses (see
\link[=analysis-information]{analysis info}) are filtered out unless their relative intensity is above this
threshold. For instance, a value of \samp{5} means that only features with an intensity five times higher than that
of the blank are kept. The relative intensity values between blanks and non-blanks are determined from the mean of
all non-zero blank intensities. Set to \code{NULL} to skip this step.}

\item{retentionRange, mzRange, mzDefectRange, chromWidthRange}{Range of retention time (in seconds), \emph{m/z}, mass
defect (defined as the decimal part of \emph{m/z} values) or chromatographic peak width (in seconds), respectively.
Features outside this range will be removed. Should be a numeric vector with length of two containing the min/max
values. The maximum can be \code{Inf} to specify no maximum range. Set to \code{NULL} to skip this step.}

\item{featQualityRange}{Used to filter features by their peak qualities/scores
(see \code{\link{calculatePeakQualities}}). Should be a named \code{list} with min/max ranges for each
quality/score to be filtered (the \code{\link{featureQualityNames}} function can be used to obtain valid names).
Example: \code{featQualityRange=list(ModalityScore=c(0.3, Inf),
SymmetryScore=c(0.5, Inf))}. Set to \code{NULL} to ignore.}

\item{groupQualityRange}{Like \code{featQualityRange}, but filters on group specific or averaged qualities/scores.}

\item{rGroups}{A character vector of replicate groups that should be kept (\code{filter}) or subtracted from
(\code{replicateGroupSubtract}).}

\item{results}{Only keep feature groups that have results in the object specified by \code{results}. Valid classes
are \code{\link{featureAnnotations}} (\emph{e.g.} formula/compound annotations) and \code{\link{components}}. Can
also be a \code{list} with multiple objects: in this case a feature group is kept if it has a result in \emph{any}
of the objects. Set to \code{NULL} to ignore.}

\item{removeBlanks}{Set to \code{TRUE} to remove all analyses that belong to replicate groups that are specified as a
blank in the \link{analysis-information}. This is useful to simplify the analyses in the specified
\code{\link{featureGroups}} object after blank subtraction. When both \code{blankThreshold} and this argument are
set, blank subtraction is performed prior to removing any analyses.}

\item{checkFeaturesSession}{If set then features and/or feature groups are removed that were selected for removal
(see \link{check-GUI}). The session files are typically generated with the \code{\link{checkFeatures}} and
\code{\link{predictCheckFeaturesSession}} functions. The value of \code{checkFeaturesSession} should either by a
path to the session file or \code{TRUE}, in which case the default session file name is used. If \code{negate=TRUE}
then all non-selected features/feature groups are removed instead.}

\item{negate}{If set to \code{TRUE} then filtering operations are performed in opposite manner.}

\item{\dots}{\setsPassedArgs1{featureGroups}}

\item{sets}{\setsWF A \code{character} with name(s) of the sets to keep (or remove if \code{negate=TRUE}).}

\item{absMinSets, relMinSets}{\setsWF Feature groups are only kept when they contain data for at least this (absolute
or relative) amount of sets. Set to \code{NULL} to ignore.}

\item{fGroups, obj}{\code{\link{featureGroups}} object to which the filter is applied.}

\item{threshold}{Minimum relative threshold (compared to mean intensity of
replicate group being subtracted) for a feature group to be \emph{not}
removed. When \samp{0} a feature group is always removed when present in
the given replicate groups.}
}
\value{
A filtered \code{\link{featureGroups}} object. Feature groups that are filtered away have their intensity set
  to zero. In case a feature group is not present in any of the analyses anymore it will be removed completely.
}
\description{
Basic rule based filtering of feature groups.
}
\details{
\code{filter} performs common rule based filtering of feature groups such as blank subtraction, minimum
  intensity and minimum replicate abundance. Removing of features occurs by zeroing their intensity values.
  Furthermore, feature groups that are left completely empty (\emph{i.e.} all intensities are zero) will be
  automatically removed.

\code{replicateGroupSubtract} removes feature groups present in a
  given set of replicate groups (unless intensities are above a given
  threshold). The replicate groups that are subtracted will be removed.
}
\section{Sets workflows}{
 \setsWFChangedMethods{

  \item \code{filter} has specific arguments to filter by (feature presence in) sets. See the argument descriptions.

  }
}

\section{Filter order}{
 When multiple arguments are specified to \code{filter}, multiple filters are applied in
  sequence. Since some of these filters may affect each other, choosing their order correctly may be important for
  effective data filtering. For instance, when an intensity filter removes features from blank analyses, a subsequent
  blank filter may not adequately perform blank subtraction. Similarly, when intensity and blank filters are executed
  after the replicate abundance filter it may be necessary to ensure minimum replicate abundance again as the
  intensity and blank filters may have removed some features within a replicate group.

  With this in mind, filters (if specified) occur in the following order:

  \enumerate{

  \item Features/feature groups selected for removal by the session specified by \code{checkFeaturesSession}.

  \item Pre-Intensity filters (\emph{i.e.} \code{preAbsMinIntensity} and \code{preRelMinIntensity}).

  \item Chromatography and mass filters (\emph{i.e} \code{retentionRange}, \code{mzRange}, \code{mzDefectRange},
  \code{chromWidthRange}, \code{featQualityRange} and \code{groupQualityRange}).

  \item Replicate abundance filters (\emph{i.e.} \code{absMinReplicateAbundance}, \code{relMinReplicateAbundance} and
  \code{maxReplicateIntRSD}).

  \item Blank filter (\emph{i.e.} blankThreshold).

  \item Intensity filters (\emph{i.e.} \code{absMinIntensity} and \code{relMinIntensity}).

  \item Replicate abundance filters (2nd time, only if previous filters affected results).

  \item General abundance filters (\emph{i.e.} \code{absMinAnalyses}, \code{relMinAnalyses}, \code{absMinReplicates},
  \code{relMinReplicates}, \code{absMinFeatures} and \code{relMinFeatures}).

  \item Replicate group filter (\emph{i.e.} \code{rGroups}), results filter (\emph{i.e.} \code{results}) and blank
  analyses removal (\emph{i.e.} if \code{removeBlanks=TRUE}).

  }

  If another filtering order is desired then \code{filter} should be called multiple times with only one filter
  argument at a time.
}

\seealso{
\code{\link{featureGroups-class}} and \code{\link{groupFeatures}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/report.R
\name{reporting}
\alias{reporting}
\alias{reportCSV,featureGroups-method}
\alias{reportCSV}
\alias{reportPDF,featureGroups-method}
\alias{reportPDF}
\alias{reportHTML,featureGroups-method}
\alias{reportHTML}
\title{Report feature group data}
\usage{
\S4method{reportCSV}{featureGroups}(
  fGroups,
  path = "report",
  reportFeatures = FALSE,
  formulas = NULL,
  formulasNormalizeScores = "max",
  formulasExclNormScores = NULL,
  compounds = NULL,
  compoundsNormalizeScores = "max",
  compoundsExclNormScores = c("score", "individualMoNAScore", "annoTypeCount"),
  compsCluster = NULL,
  components = NULL,
  retMin = TRUE,
  clearPath = FALSE
)

\S4method{reportPDF}{featureGroups}(
  fGroups,
  path = "report",
  reportFGroups = TRUE,
  formulas = NULL,
  formulasTopMost = 5,
  formulasNormalizeScores = "max",
  formulasExclNormScores = NULL,
  reportFormulaSpectra = TRUE,
  compounds = NULL,
  compoundsNormalizeScores = "max",
  compoundsExclNormScores = c("score", "individualMoNAScore", "annoTypeCount"),
  compoundsOnlyUsedScorings = TRUE,
  compoundsTopMost = 5,
  compsCluster = NULL,
  components = NULL,
  MSPeakLists = NULL,
  retMin = TRUE,
  EICGrid = c(2, 1),
  EICRtWindow = 20,
  EICMzExpWindow = 0.001,
  EICTopMost = 1,
  EICTopMostByRGroup = TRUE,
  EICOnlyPresent = TRUE,
  clearPath = FALSE
)

\S4method{reportHTML}{featureGroups}(
  fGroups,
  path = "report",
  reportPlots = c("chord", "venn", "upset", "eics", "formulas"),
  formulas = NULL,
  formulasTopMost = 5,
  formulasNormalizeScores = "max",
  formulasExclNormScores = NULL,
  compounds = NULL,
  compoundsNormalizeScores = "max",
  compoundsExclNormScores = c("score", "individualMoNAScore", "annoTypeCount"),
  compoundsOnlyUsedScorings = TRUE,
  compoundsTopMost = 5,
  compsCluster = NULL,
  includeMFWebLinks = "compounds",
  components = NULL,
  interactiveHeat = FALSE,
  MSPeakLists = NULL,
  specSimParams = getDefSpecSimParams(),
  retMin = TRUE,
  EICRtWindow = 20,
  EICMzExpWindow = 0.001,
  EICTopMost = 1,
  EICTopMostByRGroup = TRUE,
  EICOnlyPresent = TRUE,
  selfContained = TRUE,
  optimizePng = FALSE,
  clearPath = FALSE,
  openReport = TRUE,
  noDate = FALSE
)
}
\arguments{
\item{fGroups}{The \code{\link{featureGroups}} object that should be used for
reporting data.}

\item{path}{The destination file path for files generated during reporting.
Will be generated if needed.}

\item{reportFeatures}{If set to \code{TRUE} then for each analysis a
\file{.csv} file will be generated with information about its detected
features.}

\item{formulas, compounds, compsCluster, components}{Further objects
(\code{\link{formulas}}, \code{\link{compounds}},
\code{\link{compoundsCluster}}, \code{\link{components}}) that should be
reported. Specify \code{NULL} to skip reporting a particular object.}

\item{compoundsNormalizeScores, formulasNormalizeScores}{A \code{character} that specifies how normalization of
annotation scorings occurs. Either
\code{"none"} (no normalization),
\code{"max"} (normalize to max value) or \code{"minmax"} (perform min-max
normalization). Note that normalization of negative scores (e.g. output by
\command{SIRIUS}) is always performed as min-max. Furthermore, currently
normalization for \code{compounds} takes the original min/max scoring
values into account when candidates were generated. Thus, for
\code{compounds} scoring, normalization is not affected when candidate
results were removed after they were generated (\emph{e.g.} by use of
\code{filter}).}

\item{compoundsExclNormScores, formulasExclNormScores}{A
  \code{character} vector specifying any compound scoring names that
  should \emph{not} be normalized. Set to \code{NULL} to normalize all
  scorings. Note that whether any normalization occurs is set by the
  \code{compoundsExclNormScores,formulasExclNormScores} argument.

  For \code{compounds}: By default \code{score} and
  \code{individualMoNAScore} are set to mimic the behavior of the
  \command{MetFrag} web interface.}

\item{retMin}{If \code{TRUE} then report retention times in minutes
(otherwise seconds).}

\item{clearPath}{If \code{TRUE} then the destination path will be
(recursively) removed prior to reporting.}

\item{reportFGroups}{If \code{TRUE} then feature group data will be reported.}

\item{formulasTopMost, compoundsTopMost}{Only this amount of top ranked
candidate formulae/compounds are reported. Lower values may significantly
speed up reporting. Set to \code{NULL} to ignore.}

\item{reportFormulaSpectra}{If \code{TRUE} then explained MS/MS spectra (if
available) for candidate formulae will be reported. Specifying
\code{formulas} and setting this argument to \code{FALSE} still allows
further annotation of compound MS/MS spectra.}

\item{compoundsOnlyUsedScorings}{If \code{TRUE} then only scorings are plotted
that actually have been used to rank data (see the \code{scoreTypes}
argument to \code{\link{generateCompoundsMetFrag}} for more details).}

\item{MSPeakLists}{A \code{\link{MSPeakLists}} object that is
\emph{mandatory} when spectra for formulae and/or compounds will be
reported.}

\item{EICGrid}{An integer vector in the form \code{c(columns, rows)} that is
used to determine the plotting grid when reporting EICs in PDF files.}

\item{EICRtWindow, EICMzExpWindow, EICTopMost, EICTopMostByRGroup, EICOnlyPresent}{Plotting parameters passed to \code{\link{plotChroms}} (\emph{i.e.}
\code{rtWindow}, \code{mzExpWindow}, \code{topMost}, \code{topMostByRGroup}
and \code{onlyPresent} arguments).}

\item{reportPlots}{A character vector specifying what should be plotted.
Valid options are: \code{"chord"}, \code{"venn"}, \code{"upset"} (plot a
chord, Venn and UpSet diagram, respectively), \code{"eics"} (plot EICs for
individual feature groups) and \code{"formulas"} (plot annotated formula
spectra). Set to \code{"none"} to plot none of these.}

\item{includeMFWebLinks}{A \code{character} specifying to which feature
groups a web link should be added in the annotation page to
\href{https://msbi.ipb-halle.de/MetFragBeta/index.xhtml}{MetFragWeb}.
Options are: \code{"compounds"} (only to those with compounds results),
\code{"MSMS"} (only to those with MSMS peak lists) or \code{"none"}.}

\item{interactiveHeat}{If \code{TRUE} an interactive heatmap HTML widget will
be generated to display hierarchical clustering results. Set to
\code{FALSE} for a 'regular' static plot.}

\item{specSimParams}{A named \code{list} with parameters that influence the calculation of MS spectra similarities.
See the \link[=specSimParams]{spectral similarity parameters} documentation for more details.}

\item{selfContained}{If \code{TRUE} the output will be a standalone HTML file
which contains all graphics and script dependencies. When \code{FALSE}, the
latter will be placed in an additional directory (\file{report_files})
which should remain present when viewing the output file. Especially on
Windows, a non-self contained output might be desirable when reporting
large amounts of data to prevent \command{pandoc} from running out of
memory.}

\item{optimizePng}{If \code{TRUE} then \command{pngquant} is used to reduce
the size of generated graphics. A significant reduction in disk space usage
may be seen, however, at the cost additional processing time. Multiple
\command{pngquant} processes will be executed in parallel, which can be
configured with \option{patRoon.MP.maxProcs} (parallelization will always
happen with the \code{"classic"} method, see
\link[=patRoon-package]{patRoon options}).}

\item{openReport}{If set to \code{TRUE} then the output report file will be
opened with the system browser.}

\item{noDate}{If \code{TRUE} then the current date is not added to the
report. This is mainly used for testing and its main purpose is to
guarentees equal report files when `reportHTML()` is called multiple times
with equal arguments.}
}
\description{
Functionality to report data produced by most workflow steps such as
features, feature groups, calculated chemical formulae and tentatively
identified compounds.
}
\details{
These functions are usually called at the very end of the workflow. It is
used to report various data on features and feature groups. In addition,
these functions may be used for reporting formulae and/or compounds that were
generated for the specified feature groups. Data can be reported in tabular
form (\emph{i.e.} \file{.csv} files) by \code{reportCSV} or graphically by
\code{reportPDF} and \code{reportHTML}. The latter functions will plot for
instance chromatograms and annotated mass spectra, which are useful to get a
graphical overview of results.

All functions have a wide variety of arguments that influence the reporting
process. Nevertheless, most parameters are optional and only required to be
given for fine tuning. In addition, only those objects (\emph{e.g.} formulae,
compounds, clustering) that are desired to be reported need to be specified.

\code{reportCSV} generates tabular data (\emph{i.e.} \file{.csv}
  files) for given data to be reported. This may also be useful to allow
  import by other tools for post processing.

\code{reportPDF} will report graphical data (\emph{e.g.}
  chromatograms and mass spectra) within PDF files. Compared
  to \code{reportHTML} this function may be faster and yield smaller report
  files, however, its functionality is a bit more basic and generated data is
  more 'scattered' around.

\code{reportHTML} will report graphical data (\emph{e.g.}
  chromatograms and mass spectra) and summary information in an easy
  browsable \code{HTML} file using \link{rmarkdown}, \link{flexdashboard} and
  \link{knitr}.
}
\note{
Any formulae and compounds for feature groups which are not present
  within \code{fGroups} (\emph{i.e.} because it has been subset afterwards)
  will not be reported.
}
\section{Parallelization}{
 \code{reportHTML} uses multiprocessing to parallelize
  computations. Please see the parallelization section in the handbook for
  more details and \link[=patRoon-package]{patRoon options} for configuration
  options.

 Currently, \code{reportHTML} only uses
  \code{"classic"} multiprocessing, regardless of the
  \option{patRoon.MP.method} option.
}

\references{
Creating MetFrag landing page URLs based on code from
  \href{https://github.com/Treutler/MetFamily}{MetFamily} R package. \cr\cr
  \addCitations{knitr}{2} \cr\cr \addCitations{knitr}{3}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features.R
\name{importFeatures}
\alias{importFeatures}
\title{Import features}
\usage{
importFeatures(analysisInfo, type, ...)
}
\arguments{
\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}

\item{type}{What type of data should be imported: \code{"xcms"}, \code{"xcms3"}, \code{"kpic2"} or \code{"envimass"}.}

\item{\dots}{Further arguments passed to the selected import algorithm function.}
}
\value{
An object of a class which is derived from \code{\link{features}}.
}
\description{
Generic function to import features produced by other software.
}
\details{
\code{importFeatures} is a generic function that will import features by one of the supported algorithms. The actual
  functionality is provided by algorithm specific functions such as \code{importFeaturesXCMS3} and \code{importFeaturesKPIC2}. While these
  functions may be called directly, \code{importFeatures} provides a generic interface and is therefore usually preferred.
}
\seealso{
The \code{\link{features}} output class and its methods and the algorithm specific functions:
  \code{\link{importFeaturesXCMS}}, \code{\link{importFeaturesXCMS3}}, \code{\link{importFeaturesKPIC2}}, \code{\link{importFeaturesEnviMass}}

\code{\link{findFeatures}} to find new features.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/components-intclust.R
\name{generateComponentsIntClust}
\alias{generateComponentsIntClust}
\alias{generateComponentsIntClust,featureGroups-method}
\title{Generate components based on intensity profiles}
\usage{
\S4method{generateComponentsIntClust}{featureGroups}(
  fGroups,
  method = "complete",
  metric = "euclidean",
  normFunc = max,
  average = TRUE,
  maxTreeHeight = 1,
  deepSplit = TRUE,
  minModuleSize = 1
)
}
\arguments{
\item{fGroups}{\code{\link{featureGroups}} object for which components should be generated.}

\item{method}{Clustering method that should be applied (passed to
\code{\link[fastcluster:hclust]{fastcluster::hclust}}).}

\item{metric}{Distance metric used to calculate the distance matrix (passed to \code{\link{daisy}}).}

\item{normFunc, average}{Passed to \code{\link[=as.data.table,featureGroups-method]{as.data.table}} to perform
normalization and averaging of data.}

\item{maxTreeHeight, deepSplit, minModuleSize}{Arguments used by
\code{\link{cutreeDynamicTree}}.}
}
\value{
The components are stored in objects derived from \code{\link{componentsIntClust}}.
}
\description{
Generates components based on intensity profiles of feature groups.
}
\details{
This function uses hierarchical clustering of intensity profiles to generate components. This function is called when calling \code{generateComponents} with
  \code{algorithm="intclust"}.

Hierarchical clustering is performed on normalized (and optionally replicate averaged) intensity data and
  the resulting dendrogram is automatically cut with \code{\link{cutreeDynamicTree}}. The distance matrix is
  calculated with \code{\link{daisy}} and clustering is performed with
  \code{\link[fastcluster:hclust]{fastcluster::hclust}}. The clustering of the resulting components can be further
  visualized and modified using the methods defined for \code{\link{componentsIntClust}}.
}
\section{Sets workflows}{
 In a \link[=sets-workflow]{sets workflow} normalization of feature intensities occur per
  set.
}

\references{
\addCitations{fastcluster}{1}

\insertRef{Scholle2018}{patRoon}
}
\seealso{
\code{\link{generateComponents}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups-openms.R
\name{groupFeaturesOpenMS}
\alias{groupFeaturesOpenMS}
\alias{groupFeaturesOpenMS,features-method}
\title{Group features using OpenMS}
\usage{
\S4method{groupFeaturesOpenMS}{features}(
  feat,
  rtalign = TRUE,
  QT = FALSE,
  maxAlignRT = 30,
  maxAlignMZ = 0.005,
  maxGroupRT = 12,
  maxGroupMZ = 0.005,
  extraOptsRT = NULL,
  extraOptsGroup = NULL,
  verbose = TRUE
)
}
\arguments{
\item{feat}{The \code{\link{features}} object with the features to be grouped.}

\item{rtalign}{Set to \code{TRUE} to enable retention time alignment.}

\item{QT}{If enabled, use \command{FeatureLinkerUnlabeledQT} instead of \command{FeatureLinkerUnlabeled} for feature
grouping.}

\item{maxAlignRT, maxAlignMZ}{Used for retention alignment. Maximum retention time or m/z difference (seconds/Dalton)
for feature pairing. Sets \code{-algorithm:pairfinder:distance_RT:max_difference} and
\code{-algorithm:pairfinder:distance_MZ:max_difference} otpions, respectively.}

\item{maxGroupRT, maxGroupMZ}{as \code{maxAlignRT} and \code{maxAlignMZ}, but for grouping of features. Sets
\code{-algorithm:distance_RT:max_difference} and \code{-algorithm:distance_MZ:max_difference} options,
respectively.}

\item{extraOptsRT, extraOptsGroup}{Named \code{list} containing extra options that will be passed to
\command{MapAlignerPoseClustering} or \command{FeatureLinkerUnlabeledQT/FeatureLinkerUnlabeled}, respectively. Any
options specified here will override any of the above. Example:
\code{extraOptsGroup=list("-algorithm:distance_RT:max_difference"=12)} (corresponds to setting
\code{maxGroupRT=12}). Set to \code{NULL} to ignore.}

\item{verbose}{if \code{FALSE} then no text output will be shown.}
}
\value{
An object of a class which is derived from \code{\link{featureGroups}}.

The \code{featuresSet} method (for \link[=sets-workflow]{sets workflows}) returns a
  \code{\link{featureGroupsSet}} object.
}
\description{
Group and align features with OpenMS tools
}
\details{
This function uses OpenMS to group features. This function is called when calling \code{groupFeatures} with
  \code{algorithm="openms"}.

Retention times may be aligned by the
  \href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_MapAlignerPoseClustering.html}{MapAlignerPoseClustering}
   TOPP tool. Grouping is achieved by either the
  \href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureLinkerUnlabeled.html}{FeatureLinkerUnlabeled}
   or
  \href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureLinkerUnlabeledQT.html}{FeatureLinkerUnlabeledQT}
   TOPP tools.
}
\references{
\insertRef{Rst2016}{patRoon} \cr\cr
  \href{https://pugixml.org/}{pugixml} (via \href{http://www.rcpp.org/}{Rcpp}) is used to process OpenMS XML output. \cr\cr
  \addCitations{Rcpp}{1} \cr\cr
  \addCitations{Rcpp}{2} \cr\cr
  \addCitations{Rcpp}{3}
}
\seealso{
\code{\link{groupFeatures}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/components-tps.R
\docType{class}
\name{componentsTPs-class}
\alias{componentsTPs-class}
\alias{componentsTPs}
\alias{as.data.table,componentsTPs-method}
\alias{filter,componentsTPs-method}
\alias{plotGraph,componentsTPs-method}
\title{Components based on parent and transformation product (TP) linkage.}
\usage{
\S4method{as.data.table}{componentsTPs}(x)

\S4method{filter}{componentsTPs}(
  obj,
  ...,
  retDirMatch = FALSE,
  minSpecSim = NULL,
  minSpecSimPrec = NULL,
  minSpecSimBoth = NULL,
  minFragMatches = NULL,
  minNLMatches = NULL,
  formulas = NULL,
  verbose = TRUE,
  negate = FALSE
)

\S4method{plotGraph}{componentsTPs}(obj, onlyLinked)
}
\arguments{
\item{x, obj}{A \code{componentsTPs} object.}

\item{\dots, verbose}{Further arguments passed to the base \code{\link[=filter,components-method]{filter method}}.}

\item{retDirMatch}{If set to \code{TRUE}, only keep TPs for which the retention time direction (\code{retDir}, see
Details in \link{componentsTPs}) matches with the observed direction. TPs will never be removed if the
expected/observed direction is \samp{0} (\emph{i.e.} unknown or not significantly different than the parent).}

\item{minSpecSim, minSpecSimPrec, minSpecSimBoth}{The minimum spectral similarity of a TP compared to its parent
(\samp{0-1}). The \code{minSpecSimPrec} and \code{minSpecSimBoth} apply to binned data that is shifted with the
\code{"precursor"} and \code{"both"} method, respectively (see \link[=specSimParams]{MS spectral similarity
parameters} for more details). Set to \code{NULL} to ignore.}

\item{minFragMatches, minNLMatches}{Minimum number of parent/TP fragment and neutral loss matches, respectively. Set
to \code{NULL} to ignore. See the \verb{Linking parents and transformation products} section in
\code{\link{generateComponentsTPs}} for more details.}

\item{formulas}{A \code{\link{formulas}} object. The formula annotation data in this object is to verify if elemental
additions/subtractions from metabolic logic reactions are possible (hence, it only works with data from
\code{\link{generateTPsLogic}}). To verify elemental additions, only TPs with at least one candidate formula that
has these elements are kept. Similarly, for elemental subtractions, any of the parent candidate formulae must
contain the subtraction elements. Note that TPs are currently not filtered if either the parent or the TP has no
formula annotations. Set to \code{NULL} to ignore.}

\item{negate}{If \code{TRUE} then filters are applied in opposite manner.}

\item{onlyLinked}{If \code{TRUE} then only components with links are plotted.}
}
\value{
\code{filter} returns a filtered \code{componentsTPs} object.

\code{plotGraph} returns the result of \code{\link{visNetwork}}.
}
\description{
This class is derived from \code{\link{components}} and is used to store components that result from linking feature
groups that are (predicted to be) parents with feature groups that (are predicted to be) transformation products. For
more details, see \code{\link{generateComponentsTPs}}.
}
\section{Methods (by generic)}{
\itemize{
\item \code{as.data.table}: Returns all component data as a \code{\link{data.table}}.

\item \code{filter}: Provides various rule based filtering options to clean and prioritize TP data.

\item \code{plotGraph}: Plots an interactive network graph for linked components. Components are linked with each
other if one or more transformation products overlap. The graph is constructed with the \pkg{\link{igraph}} package
and rendered with \pkg{\link{visNetwork}}.
}}

\note{
The intensity values for components (used by \code{plotSpectrum}) are set
  to a dummy value (1) as no single intensity value exists for this kind of
  components.
}
\section{S4 class hierarchy}{
 \itemize{   \item{\code{\link{components}}}   \itemize{     \item{\strong{\code{\link{componentsTPs}}}}   } }
}

\references{
\addCitations{igraph}{1} \cr \cr \addCitations{visNetwork}{1}
}
\seealso{
\code{\link{components}} for other relevant methods and \code{\link{generateComponents}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R, R/feature_groups-optimize.R,
%   R/features-optimize.R
\name{feature-optimization}
\alias{feature-optimization}
\alias{optimizeFeatureGrouping}
\alias{generateFGroupsOptPSet}
\alias{getDefFGroupsOptParamRanges}
\alias{optimizeFeatureFinding}
\alias{generateFeatureOptPSet}
\alias{getDefFeaturesOptParamRanges}
\title{Optimization of feature finding and grouping parameters}
\usage{
optimizeFeatureGrouping(
  features,
  algorithm,
  ...,
  templateParams = list(),
  paramRanges = list(),
  maxIterations = 50,
  maxModelDeviation = 0.1,
  parallel = TRUE
)

generateFGroupsOptPSet(algorithm, ...)

getDefFGroupsOptParamRanges(algorithm)

optimizeFeatureFinding(
  anaInfo,
  algorithm,
  ...,
  templateParams = list(),
  paramRanges = list(),
  isoIdent = if (algorithm == "openms") "OpenMS" else "IPO",
  checkPeakShape = "none",
  CAMERAOpts = list(),
  maxIterations = 50,
  maxModelDeviation = 0.1,
  parallel = TRUE
)

generateFeatureOptPSet(algorithm, ...)

getDefFeaturesOptParamRanges(algorithm, method = "centWave")
}
\arguments{
\item{features}{A \code{\link{features}} object with the features that should
be used to optimize grouping.}

\item{algorithm}{The algorithm used for finding or grouping features (see \code{\link{findFeatures}} and
\code{\link{groupFeatures}}).}

\item{\dots}{One or more lists with parameter sets (see below) (for \code{optimizeFeatureFinding} and
\code{optimizeFeatureGrouping}). Alternatively, named arguments that set (and possibly override) the parameters
that should be returned from \code{generateFeatureOptPSet} or \code{generateFGroupsOptPSet}.}

\item{templateParams}{Template parameter set (see below).}

\item{paramRanges}{A list with vectors containing absolute parameter ranges (minimum/maximum) that constrain numeric
parameters choosen during experiments. See the \code{\link{getDefFeaturesOptParamRanges}} and
\code{\link{getDefFGroupsOptParamRanges}} functions for defaults. Values should be \code{Inf} when no limit should
be used.}

\item{maxIterations}{Maximum number of iterations that may be performed to find optimimum values. Used to restrict
neededless long optimization procedures. In IPO this was fixed to \samp{50}.}

\item{maxModelDeviation}{See the \verb{Potential suboptimal results by optimization model} section below.}

\item{parallel}{If set to \code{TRUE} then code is executed in parallel through the \CRANpkg{futures} package. Please
see the parallelization section in the handbook for more details.}

\item{anaInfo}{\link[=analysis-information]{Analysis info table} (passed to \code{\link{findFeatures}}).}

\item{isoIdent}{Sets the algorithm used to identify isotopes. Valid values
are: \code{"IPO"}, \code{"CAMERA"} and \code{"OpenMS"}. The latter can only
be used when OpenMS is used to find features, and is highly recommended in
this situation.}

\item{checkPeakShape}{Additional peak shape checking of isotopes. Only used
if \code{isoIdent="IPO"}. Valid values: \code{"none"},
\code{"borderIntensity"}, \code{"sinusCurve"} or \code{"normalDistr"}.}

\item{CAMERAOpts}{A \code{list} with additional arguments passed to
\code{\link[CAMERA:findIsotopes-methods]{CAMERA::findIsotopes}} when \code{isoIdent="CAMERA"}.}

\item{method}{Method used by XCMS to find features (only if \code{algorithm="xcms"}).}
}
\value{
The \code{optimizeFeatureFinding} and \code{optimizeFeatureGrouping} return their results in a
  \code{\link{optimizationResult}} object.
}
\description{
Automatic optimization of feature finding and grouping parameters through Design of Experiments (DoE).
}
\details{
Many different parameters exist that may affect the output quality of feature finding and grouping. To avoid time
consuming manual experimentation, functionality is provided to largely automate the optimization process. The
methodology, which uses design of experiments (DoE), is based on the excellent
\href{https://github.com/rietho/IPO}{Isotopologue Parameter Optimization (IPO) R package}. The functionality of this
package is directly integrated in patRoon. Some functionality was added or changed, however, the principle algorithm
workings are nearly identical.

Compared to IPO, the following functionality was added or changed: \itemize{

\item The code was made more generic in order to include support for other feature finding/grouping algorithms
(\emph{e.g.} OpenMS, enviPick, XCMS3). \item The methodology of \command{FeatureFinderMetabo} (OpenMS) may be used to
find isotopes.

\item The \code{maxModelDeviation} parameter was added to potentially avoid suboptimal results
(\href{https://github.com/rietho/IPO/issues/61}{issue discussed here}).

\item The use of multiple 'parameter sets' (discussed below) which, for instance, allow optimizing qualitative
paremeters more easily (see \verb{examples}).

\item More consistent optimization code for feature finding/grouping.

\item More consistent output using S4 classes (\emph{i.e.} \code{\link{optimizationResult}} class).

\item Parallelization is performed via the \CRANpkg{future} package instead of \pkg{BiocParallel}. If this is enabled
(\code{parallel=TRUE}) then any parallelization supported by the feature finding or grouping algorithm is disabled.

}
}
\section{Parameter sets}{
 Which parameters should be optimized is determined by a \emph{parameter set}. A set is
  defined by a named \code{list} containing the minimum and maximum starting range for each parameter that should be
  tested. For instance, the set \code{list(chromFWHM = c(5, 10), mzPPM = c(5, 15))} specifies that the
  \code{chromFWHM} and \code{mzPPM} parameters (used by OpenMS feature finding) should be optimized within a range of
  \samp{5}-\samp{10} and \samp{5}-\samp{15}, respectively. Note that this range may be increased or decreased after a
  DoE iteration in order to find a better optimum. The absolute limits are controlled by the \code{paramRanges}
  function argument.

  Multiple parameter sets may be specified (\emph{i.e.} through the \dots function argument). In this situation, the
  optimization algorithm is repeated for each set, and the final optimum is determined from the parameter set with
  the best response. The \code{templateParams} function argument may be useful in this case to define a template for
  each parameter set. Actual parameter sets are then constructed by joining each parameter set with the set specified
  for \code{templateParams}. When a parameter is defined in both a regular and template set, the parameter in the
  regular set takes precedence.

  Parameters that should not be optimized but still need to be set for the feature finding/grouping functions should
  also be defined in a (template) parameter set. Which parameters should be optimized is determined whether its value
  is specified as a vector range or a single fixed value. For instance, when a set is defined as \code{list(chromFWHM
  = c(5, 10), mzPPM = 5)}, only the \code{chromFWHM} parameter is optimized, whereas \code{mzPPM} is kept constant at
  \samp{5}.

  Using multiple parameter sets with differing fixed values allows optimization of qualitative values (see examples
  below).

  The parameters specified in parameter sets are directly passed through the \code{\link{findFeatures}} or
  \code{\link{groupFeatures}} functions. Hence, grouping and retention time alignment parameters used by XCMS should
  (still) be set through the \code{groupArgs} and \code{retcorArgs} parameters.

  \strong{NOTE:} For XCMS3, which normally uses parameter classes for settings its options, the parameters must be
  defined in a named list like any other algorithm. The set parameters are then used passed to the constructor of the
  right parameter class object (e.g. \code{\link{CentWaveParam}}, \code{\link{ObiwarpParam}}). For grouping/alignment
  sets, these parameters need to be specified in nested lists called \code{groupParams} and \code{retAlignParams},
  respectively (similar to \code{groupArgs}/\code{retcorArgs} for \code{algorithm="xcms"}). Finally, the underlying
  XCMS method to be used should be defined in the parameter set (\emph{i.e.} by setting the \code{method} field for
  feature parameter sets and the \code{groupMethod} and \code{retAlignMethod} for grouping/aligning parameter sets).
  See the examples below for more details.

  \strong{NOTE:} Similar to IPO, the \code{peakwidth} and \code{prefilter} parameters for XCMS feature finding should
  be split in two different values: \itemize{

  \item The minimum and maximum ranges for \code{peakwidth} are optimized by setting \code{min_peakwidth} and
  \code{max_peakwidth}, respectively.

  \item The \code{k} and \code{I} parameters contained in \code{prefilter} are split in \code{prefilter} and
  \code{value_of_prefilter}, respectively.

  }

  \emph{Similary}, for KPIC2, the following parameters should be split: \itemize{

  \item the \code{width} parameter (feature optimization) is optimized by specifying the \code{min_width} and
  \code{max_width} parameters.

  \item the \code{tolerance} and \code{weight} parameters (feature grouping optimization) are optimized by setting
  \code{mz_tolerance}/\code{rt_tolerance} and \code{mz_weight}/\code{rt_weight} parameters, respectively.

  }
}

\section{Functions}{
 The \code{optimizeFeatureFinding} and \code{optimizeFeatureGrouping} are the functions to be used
  to optimize parameters for feature finding and grouping, respectively. These functions are analogous to
  \code{\link[IPO]{optimizeXcmsSet}} and \code{\link[IPO]{optimizeRetGroup}} from \pkg{IPO}.

  The \code{generateFeatureOptPSet} and \code{generateFGroupsOptPSet} functions may be used to generate a parameter
  set for feature finding and grouping, respectively. Some algorithm dependent default parameter optimization ranges
  will be returned. These functions are analogous to \code{\link[IPO]{getDefaultXcmsSetStartingParams}} and
  \code{\link[IPO]{getDefaultRetGroupStartingParams}} from \pkg{IPO}. However, unlike their IPO counterparts, these
  functions will not output default fixed values. The \code{generateFGroupsOptPSet} will only generate defaults for
  density grouping if \code{algorithm="xcms"}.

  The \code{getDefFeaturesOptParamRanges} and \code{getDefFGroupsOptParamRanges} return the default absolute
  optimization parameter ranges for feature finding and grouping, respectively. These functions are useful if you
  want to set the \code{paramRanges} function argument.
}

\section{Potential suboptimal results by optimization model}{
 After each experiment iteration an optimimum parameter
  set is found by generating a model containing the tested parameters and their responses. Sometimes the actual
  response from the parameters derived from the model is actually signficantly lower than expected. When the response
  is lower than the maximum reponse found during the experiment, the parameters belonging to this experimental
  maximum may be choosen instead. The \code{maxModelDeviation} argument sets the maximum deviation in response
  between the modelled and experimental maxima. The value is relative: \samp{0} means that experimental values will
  always be favored when leading to improved responses, whereas \code{1} will effectively disable this procedure (and
  return to 'regular' IPO behaviour).
}

\section{Source}{
 The code and methodology is a direct adaptation from the \href{https://github.com/rietho/IPO}{IPO R
  package}.
}

\examples{
\donttest{# example data from patRoonData package
dataDir <- patRoonData::exampleDataPath()
anaInfo <- generateAnalysisInfo(dataDir)
anaInfo <- anaInfo[1:2, ] # only focus on first two analyses (e.g. training set)

# optimize mzPPM and chromFWHM parameters
ftOpt <- optimizeFeatureFinding(anaInfo, "openms", list(mzPPM = c(5, 10), chromFWHM = c(4, 8)))

# optimize chromFWHM and isotopeFilteringModel (a qualitative parameter)
ftOpt2 <- optimizeFeatureFinding(anaInfo, "openms",
                                 list(isotopeFilteringModel = "metabolites (5\% RMS)"),
                                 list(isotopeFilteringModel = "metabolites (2\% RMS)"),
                                 templateParams = list(chromFWHM = c(4, 8)))

# perform grouping optimization with optimized features object
fgOpt <- optimizeFeatureGrouping(optimizedObject(ftOpt), "xcms",
                                 list(groupArgs = list(bw = c(22, 28)),
                                      retcorArgs = list(method = "obiwarp")))

# same, but using the XCMS3 interface
fgOpt2 <- optimizeFeatureGrouping(optimizedObject(ftOpt), "xcms3",
                                  list(groupMethod = "density", groupParams = list(bw = c(22, 28)),
                                       retAlignMethod = "obiwarp"))


# plot contour of first parameter set/DoE iteration
plot(ftOpt, paramSet = 1, DoEIteration = 1, type = "contour")

# generate parameter set with some predefined and custom parameters to be
# optimized.
pSet <- generateFeatureOptPSet("openms", chromSNR = c(3, 9),
                               useSmoothedInts = FALSE)
}
}
\references{
\insertRef{Libiseller2015}{patRoon}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features-safd.R
\name{findFeaturesSAFD}
\alias{findFeaturesSAFD}
\title{Find features using SAFD}
\usage{
findFeaturesSAFD(
  analysisInfo,
  profPath = NULL,
  mzRange = c(0, 400),
  maxNumbIter = 1000,
  maxTPeakW = 300,
  resolution = 30000,
  minMSW = 0.02,
  RThreshold = 0.75,
  minInt = 2000,
  sigIncThreshold = 5,
  S2N = 2,
  minPeakWS = 3,
  verbose = TRUE
)
}
\arguments{
\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}

\item{profPath}{A \code{character} vector with paths to the profile MS data for each analysis (will be re-cycled if
necessary). See the \verb{Using SAFD} section for more details.}

\item{mzRange}{The \emph{m/z} window to be imported (passed to the \code{import_files_MS1} function).}

\item{maxNumbIter, maxTPeakW, resolution, minMSW, RThreshold, minInt, sigIncThreshold, S2N, minPeakWS}{Parameters directly
passed to the \code{safd_s3D} function.}

\item{verbose}{If set to \code{FALSE} then no text output is shown.}
}
\value{
An object of a class which is derived from \code{\link{features}}.
}
\description{
Uses \href{https://bitbucket.org/SSamanipour/safd.jl/src/master/}{SAFD} to obtain features. This functionality is
still experimental. Please see the details below.
}
\details{
This function uses SAFD to automatically find features. This function is called when calling \code{findFeatures} with
  \code{algorithm="safd"}.

The support for SAFD is still experimental, and its interface might change in the future.

  In order to use SAFD, please make sure that its \code{julia} packages are installed and you have verified that
  everything works, \emph{e.g.} by running the test data.

  This algorithm supports profile and centroided MS data. If the use of profile data is desired, centroided data
  must still be available for other functionality of \code{patRoon}. The centroided data is specified through the
  'regular' \link[=analysis-information]{analysis info} mechanism. The location to any profile data is specified
  through the \code{profPath} argument (\code{NULL} for no profile data). The base file names (\emph{i.e.} the file
  name without path and extension) of both centroid and profile data must be the same. Furthermore, the format of the
  profile data must be \file{mzXML}.
}
\section{Parallelization}{
 \code{findFeaturesSAFD} uses multiprocessing to parallelize
  computations. Please see the parallelization section in the handbook for
  more details and \link[=patRoon-package]{patRoon options} for configuration
  options.

 Note that for caching purposes, the analyses files must always exist on the local host
  computer, even if it is not participating in computations.
}

\references{
\insertRef{Samanipour2019}{patRoon}
}
\seealso{
\code{\link{findFeatures}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups-screening.R,
%   R/feature_groups-screening-set.R
\docType{class}
\name{featureGroupsScreening-class}
\alias{featureGroupsScreening-class}
\alias{featureGroupsScreening}
\alias{screenInfo,featureGroupsScreening-method}
\alias{screenInfo}
\alias{show,featureGroupsScreening-method}
\alias{[,featureGroupsScreening,ANY,ANY,missing-method}
\alias{delete,featureGroupsScreening-method}
\alias{as.data.table,featureGroupsScreening-method}
\alias{annotateSuspects,featureGroupsScreening-method}
\alias{annotateSuspects}
\alias{filter,featureGroupsScreening-method}
\alias{featureGroupsScreeningSet-class}
\alias{featureGroupsScreeningSet}
\alias{screenInfo,featureGroupsScreeningSet-method}
\alias{show,featureGroupsScreeningSet-method}
\alias{[,featureGroupsScreeningSet,ANY,ANY,missing-method}
\alias{delete,featureGroupsScreeningSet-method}
\alias{as.data.table,featureGroupsScreeningSet-method}
\alias{annotateSuspects,featureGroupsScreeningSet-method}
\alias{filter,featureGroupsScreeningSet-method}
\alias{featureGroupsSetScreeningUnset-class}
\alias{featureGroupsSetScreeningUnset}
\alias{unset,featureGroupsScreeningSet-method}
\title{Class for suspect screened feature groups.}
\usage{
\S4method{screenInfo}{featureGroupsScreening}(obj)

\S4method{show}{featureGroupsScreening}(object)

\S4method{[}{featureGroupsScreening,ANY,ANY,missing}(x, i, j, ..., rGroups, suspects = NULL, drop = TRUE)

\S4method{delete}{featureGroupsScreening}(obj, i = NULL, j = NULL, ...)

\S4method{as.data.table}{featureGroupsScreening}(x, ..., collapseSuspects = ",", onlyHits = FALSE)

\S4method{annotateSuspects}{featureGroupsScreening}(
  fGroups,
  MSPeakLists,
  formulas,
  compounds,
  absMzDev = 0.005,
  relMinMSMSIntensity = 0.05,
  simMSMSMethod = "cosine",
  checkFragments = c("mz", "formula", "compound"),
  formulasNormalizeScores = "max",
  compoundsNormalizeScores = "max",
  IDFile = system.file("misc", "IDLevelRules.yml", package = "patRoon"),
  logPath = file.path("log", "ident")
)

\S4method{filter}{featureGroupsScreening}(
  obj,
  ...,
  onlyHits = NULL,
  selectHitsBy = NULL,
  selectBestFGroups = FALSE,
  maxLevel = NULL,
  maxFormRank = NULL,
  maxCompRank = NULL,
  minAnnSimForm = NULL,
  minAnnSimComp = NULL,
  minAnnSimBoth = NULL,
  absMinFragMatches = NULL,
  relMinFragMatches = NULL,
  negate = FALSE
)

\S4method{screenInfo}{featureGroupsScreeningSet}(obj)

\S4method{show}{featureGroupsScreeningSet}(object)

\S4method{[}{featureGroupsScreeningSet,ANY,ANY,missing}(x, i, j, ..., rGroups, suspects = NULL, sets = NULL, drop = TRUE)

\S4method{delete}{featureGroupsScreeningSet}(obj, i = NULL, j = NULL, ...)

\S4method{as.data.table}{featureGroupsScreeningSet}(x, ..., collapseSuspects = ",", onlyHits = FALSE)

\S4method{annotateSuspects}{featureGroupsScreeningSet}(
  fGroups,
  MSPeakLists = NULL,
  formulas = NULL,
  compounds = NULL,
  ...
)

\S4method{filter}{featureGroupsScreeningSet}(
  obj,
  ...,
  onlyHits = NULL,
  selectHitsBy = NULL,
  selectBestFGroups = FALSE,
  maxLevel = NULL,
  maxFormRank = NULL,
  maxCompRank = NULL,
  minAnnSimForm = NULL,
  minAnnSimComp = NULL,
  minAnnSimBoth = NULL,
  absMinFragMatches = NULL,
  relMinFragMatches = NULL,
  negate = FALSE
)

\S4method{unset}{featureGroupsScreeningSet}(obj, set)
}
\arguments{
\item{obj, object, x, fGroups}{The \code{featureGroupsScreening} object.}

\item{i, j, rGroups}{Used for subsetting data analyses, feature groups and
replicate groups, see \code{\link{featureGroups}}.}

\item{\dots}{Further arguments passed to the base method.}

\item{suspects}{An optional \code{character} vector with suspect names. If
specified, only \code{featureGroups} will be kept that are assigned to
these suspects.}

\item{drop}{Ignored.}

\item{collapseSuspects}{If a \code{character} then any suspects that were
matched to the same feature group are collapsed to a single row and suspect
names are separated by the value of \code{collapseSuspects}. If \code{NULL}
then no collapsing occurs, and each suspect match is reported on a single
row. Note that some columns will not be reported when collapsing is
enabled.}

\item{onlyHits}{For \code{as.data.table}: if \code{TRUE} then only feature groups with suspect hits are reported.

  For \code{filter} \itemize{

  \item if \code{negate=FALSE} and \code{onlyHits=TRUE} then all feature groups without suspect hits will be removed.
  Otherwise nothing will be done.

  \item if \code{negate=TRUE} then \code{onlyHits=TRUE} will select feature groups without suspect hits,
  \code{onlyHits=FALSE} will only retain feature groups with suspect matches and this filter is ignored if
  \code{onlyHits=NULL}.

  }}

\item{MSPeakLists, formulas, compounds}{Annotation data (\code{\link{MSPeakLists}}, \code{\link{formulas}} and
\code{\link{compounds}}) obtained for this \code{featureGroupsScreening} object. All arguments can be \code{NULL}
to exclude it from the annotation.}

\item{absMzDev}{Maximum absolute \emph{m/z} deviation.}

\item{relMinMSMSIntensity}{Minimum relative intensity (\samp{0-1}) threshold applied when calculating annotation
similarities.}

\item{simMSMSMethod}{Either \code{"cosine"} or \code{"jaccard"}: used to compare MS/MS peak lists for annotation
similarity calculation.}

\item{checkFragments}{Which type(s) of MS/MS fragments from workflow data should be checked to evaluate the number of
suspect fragment matches (\emph{i.e.} from the \code{fragments_mz}/\code{fragments_formula} columns in the suspect
list). Valid values are: \code{"mz"}, \code{"formula"}, \code{"compounds"}. The former uses \emph{m/z} values in
the specified \code{MSPeakLists} object, whereas the others use the formulae that were annotated to MS/MS peaks in
the given \code{formulas} or \code{compounds} objects. Multiple values are possible: in this case the maximum
number of fragment matches will be reported.}

\item{compoundsNormalizeScores, formulasNormalizeScores}{A \code{character} that specifies how normalization of
annotation scorings occurs. Either

\code{"max"} (normalize to max value) or \code{"minmax"} (perform min-max
normalization). Note that normalization of negative scores (e.g. output by
\command{SIRIUS}) is always performed as min-max. Furthermore, currently
normalization for \code{compounds} takes the original min/max scoring
values into account when candidates were generated. Thus, for
\code{compounds} scoring, normalization is not affected when candidate
results were removed after they were generated (\emph{e.g.} by use of
\code{filter}).}

\item{IDFile}{A file path to a YAML file with rules used for estimation of identification levels. See the
\verb{Suspect annotation} section for more details. If not specified then a default rules file will be used.}

\item{logPath}{A directory path to store logging information. If \code{NULL} then logging is disabled.}

\item{selectHitsBy}{Should be \code{"intensity"} or \code{"level"}. For cases
where the same suspect is matched to multiple feature groups, only the
suspect to the feature group with highest mean intensity
(\code{selectHitsBy="intensity"}) or best identification level
(\code{selectHitsBy="level"}) is kept. In case of ties only the first hit
is kept. Set to \code{NULL} to ignore this filter. If \code{negate=TRUE}
then only those hits with lowest mean intensity/poorest identification
level are kept.}

\item{selectBestFGroups}{If \code{TRUE} then for any cases where a single
feature group is matched to several suspects only the suspect assigned to
the feature group with best identification score is kept. In case of ties
only the first is kept.}

\item{maxLevel, maxFormRank, maxCompRank, minAnnSimForm, minAnnSimComp, minAnnSimBoth}{Filter suspects by maximum identification level (\emph{e.g.} \code{"3a"}),
formula/compound rank or with minimum formula/compound/combined annotation
similarity. Set to \code{NULL} to ignore.}

\item{absMinFragMatches, relMinFragMatches}{Only retain suspects with this
minimum number MS/MS matches with the fragments specified in the suspect
list (\emph{i.e.} \code{fragments_mz}/\code{fragments_formula}).
\code{relMinFragMatches} sets the minimum that is relative (\samp{0-1}) to
the maximum number of MS/MS fragments specified in the \code{fragments_*}
columns of the suspect list. Set to \code{NULL} to ignore.}

\item{negate}{If set to \code{TRUE} then filtering operations are performed
in opposite manner.}

\item{sets}{\setsWF A \code{character} with name(s) of the sets to keep (or remove if \code{negate=TRUE}).}

\item{set}{\setsWF The name of the set.}
}
\value{
\code{annotateSuspects} returns a \code{featureGroupsScreening} object, which is a
  \code{\link{featureGroups}} object amended with annotation data.

\code{filter} returns a filtered \code{featureGroupsScreening}
  object.
}
\description{
This class derives from \code{\link{featureGroups}} and adds suspect screening information.
}
\section{Methods (by generic)}{
\itemize{
\item \code{screenInfo}: Returns a table with screening information
(see \code{screenInfo} slot).

\item \code{show}: Shows summary information for this object.

\item \code{[}: Subset on analyses, feature groups and/or
suspects.

\item \code{as.data.table}: Obtain a summary table (a
\code{\link{data.table}}) with retention, \emph{m/z}, intensity and
optionally other feature data. Furthermore, the output table will be merged
with information from \code{screenInfo}, such as suspect names and other
properties and annotation data.

\item \code{annotateSuspects}: Incorporates annotation data obtained during the workflow to annotate suspects
with matched known MS/MS fragments, formula/candidate ranks and automatic estimation of identification levels. See
the \verb{Suspect annotation} section for more details. The estimation of identification levels for each suspect is
logged in the \code{log/ident} directory.

\item \code{filter}: Performs rule based filtering. This method
builds on the comprehensive filter functionality from the base
\code{\link{filter,featureGroups-method}}. It adds several filters to
select \emph{e.g.} the best ranked suspects or those with a minimum
estimated identification level. \strong{NOTE}: most filters \emph{only}
affect suspect hits, not feature groups. Set \code{onlyHits=TRUE} to
subsequently remove any feature groups that lost any suspect matches due to
other filter steps.
}}

\section{Slots}{

\describe{
\item{\code{screenInfo}}{A (\code{\link{data.table}}) with results from suspect screening. This table will be amended with
annotation data when \code{annotateSuspects} is run.}
}}

\note{
The \code{relMinMSMSIntensity} filter argument to \code{annotateSuspects} is applied \emph{after} removing the
  precursor ion from the peak lists (if present). Thus, intensity scales may be different when this filter is applied
  when the most abundant peak resulted from the precursor ion.

\code{filter} removes suspect hits with \code{NA} values when any of
  the filters related to minimum or maximum values are applied (unless
  \code{negate=TRUE}).
}
\section{Suspect annotation}{
 The \code{annotateSuspects} method is used to annotate suspects after
  \code{\link{screenSuspects}} was used to collect suspect screening results and other workflow steps such as formula
  and compound annotation steps have been completed. The annotation results, which can be acquired with the
  \code{as.data.table} and \code{screenInfo} methods, amends the current screening data with the following columns:

  \itemize{

  \item \code{formRank},\code{compRank} The rank of the suspect within the formula/compound annotation results.

  \item \code{annSimForm},\code{annSimComp},\code{annSimBoth} A similarity measure between measured and annotated
  MS/MS peaks from annotation of formulae, compounds or both. The similarity is calculated as the spectral similarity
  between a peaklist with (a) all MS/MS peaks and (b) only annotated peaks. Thus, a value of one means that all MS/MS
  peaks were annotated. If both formula and compound annotations are available then \code{annSimBoth} is calculated
  after combining all the annotated peaks, otherwise \code{annSimBoth} equals the available value for
  \code{annSimForm} or \code{annSimComp}. The similarity calculation can be configured with the
  \code{relMinMSMSIntensity} and \code{simMSMSMethod} arguments to \code{annotateSuspects}.

  \item \code{maxFrags} The maximum number of MS/MS fragments that can be matched for this suspect (based on the
  \code{fragments_*} columns from the suspect list).

  \item \code{maxFragMatches},\code{maxFragMatchesRel} The absolute and relative amount of experimental MS/MS peaks
  that were matched from the fragments specified in the suspect list. The value for \code{maxFragMatchesRel} is
  relative to the value for \code{maxFrags}. The calculation of this column is influenced by the
  \code{checkFragments} argument to \code{annotateSuspects}.

  \item \code{estIDLevel} Provides an \emph{estimation} of the identification level, roughly following that of
  \insertCite{Schymanski2014}{patRoon}. However, please note that this value is only an estimation, and manual
  interpretation is still necessary to assign final identification levels. The estimation is done through a set of
  rules, see the \verb{Identification level rules} section below.

  }

  Note that only columns are present is sufficient data is available for their calculation.
}

\section{Identification level rules}{
 The estimation of identification levels is configured through a YAML file which
  specifies the rules for each level. The default file is shown below.

 \preformatted{1:
    suspectFragments: 3
    retention: 12
2a:
    individualMoNAScore:
        min: 0.9
        higherThanNext: .inf
    rank:
        max: 1
        type: compound
3a:
    individualMoNAScore: 0.4
3b:
    suspectFragments: 3
3c:
    annMSMSSim:
        type: compound
        min: 0.7
4a:
    annMSMSSim:
        type: formula
        min: 0.7
        higherThanNext: 0.2
    isoScore:
        min: 0.5
        higherThanNext: 0.2
    rank:
        max: 1
        type: formula
4b:
    isoScore:
        min: 0.9
        higherThanNext: 0.2
    rank:
        max: 1
        type: formula
5:
    all: yes
}

 Most of the file should be self-explanatory. Some notes:

  \itemize{

  \item Each rule is either a field of \code{suspectFragments} (minimum number of MS/MS fragments matched from
  suspect list), \code{retention} (maximum retention deviation from suspect list), \code{rank} (the maximum
  annotation rank from formula or compound annotations), \code{all} (this level is always matched) or any of the
  scorings available from the formula or compound annotations.

  \item In case any of the rules could be applied to either formula or compound annotations, the annotation type must
  be specified with the \code{type} field (\code{formula} or \code{compound}).

  \item Identification levels should start with a number and may optionally be followed by a alphabetic character.
  The lowest levels are checked first.

  \item If \code{relative=yes} then the relative scoring will be used for testing.

  \item For \code{suspectFragments}: if the number of fragments from the suspect list (\code{maxFrags} column) is
  less then the minimum rule value, the minimum is adjusted to the number of available fragments.

  }

  A template rules file can be generated with the \code{\link{genIDLevelRulesFile}} function, and this file can
  subsequently passed to \code{annotateSuspects}. The file format is highly flexible and (sub)levels can be added or
  removed if desired. Note that the default file is currently only suitable when annotation is performed with GenForm
  and MetFrag, for other algorithms it is crucial to modify the rules.
}

\section{S4 class hierarchy}{
 \itemize{   \item{\code{\link{featureGroups}}}   \itemize{     \item{\strong{\code{\link{featureGroupsScreening}}}}     \itemize{       \item{\code{\link{featureGroupsSetScreeningUnset}}}     }   } }
}

\section{Source}{
 Cosine spectral similarity calculation was based on the code from \code{SpectrumSimilarity()}
  function of \pkg{OrgMassSpecR}.
}

\section{Sets workflows}{
 \setsWFClass{featureGroupsScreeningSet}{featureGroupsScreening}

  \setsWFNewMethodsSO{featureGroupsScreeningUnset}{Only the screening results present in the specified set are kept.}

  \setsWFChangedMethods{

  \item \code{annotateSuspects} Suspect annotation is performed per set. Thus, formula/compound ranks, estimated
  identification levels etc are calculated for each set. Subsequently, these results are merged in the final
  \code{screenInfo}. In addition, an overall \code{formRank} and \code{compRank} column is created based on the
  rankings of the suspect candidate in the set consensus data. Furthermore, an overall \code{estIDLevel} is generated
  that is based on the 'best' estimated identification level among the sets data (\emph{i.e.} the lowest). In case
  there is a tie between sub-levels (\emph{e.g.} \samp{3a} and \samp{3b}), then the sub-level is stripped
  (\emph{e.g.} \samp{3}).

  \item \code{filter} All filters releated to estimated identification levels and formula/compound rankings  are
  applied to the overall set data (see above). All others are applied to set specific data: in this case candidates
  are only removed if none of the set data confirms to the filter.

  }

  This class derives also from \code{\link{featureGroupsSet}}. Please see its documentation for more relevant details
  with sets workflows.

  Note that the \code{formRank} and \code{compRank} columns are \emph{not} updated when the data is subset.
}

\references{
\insertAllCited{} \cr \cr \insertRef{Stein1994}{patRoon}
}
\seealso{
\code{\link{featureGroups}}
}
\author{
Rick Helmus <\email{r.helmus@uva.nl}>, Emma Schymanski <\email{emma.schymanski@uni.lu}> (contributions to
  identification level rules), Bas van de Velde (contributions to spectral similarity calculation).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups-xcms.R
\name{groupFeaturesXCMS}
\alias{groupFeaturesXCMS}
\alias{groupFeaturesXCMS,features-method}
\alias{groupFeaturesXCMS,featuresSet-method}
\title{Group features using XCMS (old interface)}
\usage{
\S4method{groupFeaturesXCMS}{features}(
  feat,
  rtalign = TRUE,
  loadRawData = TRUE,
  groupArgs = list(mzwid = 0.015),
  retcorArgs = list(method = "obiwarp"),
  verbose = TRUE
)

\S4method{groupFeaturesXCMS}{featuresSet}(feat, groupArgs = list(mzwid = 0.015), verbose = TRUE)
}
\arguments{
\item{feat}{The \code{\link{features}} object with the features to be grouped.}

\item{rtalign}{Set to \code{TRUE} to enable retention time alignment.}

\item{loadRawData}{Set to \code{TRUE} if analyses are available as \code{mzXML} or \code{mzML} files. Otherwise MS
data is not loaded, and some dummy data (\emph{e.g.} file paths) is used in the returned object.}

\item{groupArgs}{named \code{character vector} that may contain extra grouping parameters to be used by
\code{\link[xcms:group-methods]{xcms::group}}}

\item{retcorArgs}{named \code{character vector} that may contain extra parameters to be used by
\code{\link[xcms:retcor-methods]{xcms::retcor}}.}

\item{verbose}{if \code{FALSE} then no text output will be shown.}
}
\value{
An object of a class which is derived from \code{\link{featureGroups}}.

The \code{featuresSet} method (for \link[=sets-workflow]{sets workflows}) returns a
  \code{\link{featureGroupsSet}} object.
}
\description{
Group and align features with the legacy \code{\link[xcms]{xcmsSet}} function from the \pkg{xcms} package.
}
\details{
This function uses XCMS to group features. This function is called when calling \code{groupFeatures} with
  \code{algorithm="xcms"}.

Grouping of features and
  alignment of their retention times are performed with the \code{\link[xcms:group-methods]{xcms::group}} and
  \code{\link[xcms:retcor-methods]{xcms::retcor}} functions, respectively. Both functions have an extensive list of
  parameters to modify their behavior and may therefore be used to potentially optimize results.
}
\section{Sets workflows}{
 \code{loadRawData} and arguments related to retention time alignment are currently not
  supported for \link[=sets-workflow]{sets workflows}.
}

\references{
\addCitations{xcms}{1} \cr\cr \addCitations{xcms}{2} \cr\cr \addCitations{xcms}{3}
}
\seealso{
\code{\link{groupFeatures}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/components-ramclustr.R
\name{generateComponentsRAMClustR}
\alias{generateComponentsRAMClustR}
\alias{generateComponentsRAMClustR,featureGroups-method}
\alias{generateComponentsRAMClustR,featureGroupsSet-method}
\title{Componentization of adducts, isotopes etc. with RAMClustR}
\usage{
\S4method{generateComponentsRAMClustR}{featureGroups}(
  fGroups,
  ionization = NULL,
  st = NULL,
  sr = NULL,
  maxt = 12,
  hmax = 0.3,
  normalize = "TIC",
  absMzDev = 0.002,
  relMzDev = 5,
  minSize = 2,
  relMinReplicates = 0.5,
  RCExperimentVals = list(design = list(platform = "LC-MS"), instrument =
    list(ionization = ionization, MSlevs = 1)),
  extraOptsRC = NULL,
  extraOptsFM = NULL
)

\S4method{generateComponentsRAMClustR}{featureGroupsSet}(fGroups, ionization = NULL, ...)
}
\arguments{
\item{fGroups}{\code{\link{featureGroups}} object for which components should be generated.}

\item{ionization}{Which ionization polarity was used to generate the data: should be \code{"positive"}
  or \code{"negative"}. If the \code{featureGroups} object has adduct annotations, and \code{ionization=NULL}, the
  ionization will be detected automatically.

  \setsWF This parameter is not supported for sets workflows, as the ionization will always be detected
  automatically.}

\item{st, sr, maxt, hmax, normalize}{Arguments to tune the behaviour of feature group clustering. See their documentation
from \code{\link[RAMClustR]{ramclustR}}. When \code{st} is \code{NULL} it will be automatically calculated as the
half of the median for all chromatographic peak widths.}

\item{absMzDev}{Maximum absolute \emph{m/z} deviation. Sets the \code{mzabs.error} argument to \code{\link[RAMClustR]{do.findmain}}}

\item{relMzDev}{Maximum relative mass deviation (\acronym{ppm}). Sets the \code{ppm.error} argument to
\code{\link[RAMClustR]{do.findmain}}.}

\item{minSize}{The minimum size of a component. Smaller components than this size will be removed. See note below. Sets the \code{minModuleSize} argument to \code{\link[RAMClustR]{ramclustR}}.}

\item{relMinReplicates}{Feature groups within a component are only kept when they contain data for at least this
(relative) amount of replicate analyses. For instance, \samp{0.5} means that at least half of the replicates should
contain data for a particular feature group in a component. In this calculation replicates that are fully absent
within a component are not taken in to account. See note below.}

\item{RCExperimentVals}{A named \code{list} containing two more \code{list}s: \code{design} and \code{instrument}.
These are used to construct the \code{ExpDes} argument passed to \code{\link[RAMClustR]{ramclustR}}.}

\item{extraOptsRC, extraOptsFM}{Named \code{list} with further arguments to be passed to
\code{\link[RAMClustR]{ramclustR}} and \code{\link[RAMClustR]{do.findmain}}. Set to \code{NULL} to ignore.}

\item{\dots}{\setsWF Further arguments passed to the non-sets workflow method.}
}
\value{
A \code{\link{components}} (derived) object containing all generated components.
}
\description{
Uses \href{https://github.com/cbroeckl/RAMClustR}{RAMClustR} to generate components from feature groups which follow
similar chromatographic retention profiles and annotate their relationships (\emph{e.g.} adducts and isotopes).
}
\details{
This function uses RAMClustR to generate components. This function is called when calling \code{generateComponents} with
  \code{algorithm="ramclustr"}.

This method uses the \code{\link[RAMClustR]{ramclustR}} functions for generating the components, whereas
  \code{\link[RAMClustR]{do.findmain}} is used for annotation.
}
\note{
The default value for  \code{relMinReplicates} results in
  extra filtering, hence, the final results may be different than what the algorithm normally would return.
}
\section{Sets workflows}{
 In a \link[=sets-workflow]{sets workflow} the componentization is first performed for each
  set independently. The resulting components are then all combined in a \code{\link{componentsSet}} object. Note that
  the components themselves are never merged. The components are renamed to include the set name from which they were
  generated (\emph{e.g.} \code{"CMP1"} becomes \code{"CMP1-positive"}).
}

\references{
\insertRef{Broeckling2013}{patRoon} \cr\cr \insertRef{Broeckling2014}{patRoon}
}
\seealso{
\code{\link{generateComponents}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups-comparison.R
\docType{class}
\name{featureGroups-compare}
\alias{featureGroups-compare}
\alias{featuresFromFeatGroups-class}
\alias{featuresFromFeatGroups}
\alias{featuresConsensus-class}
\alias{featuresConsensus}
\alias{featureGroupsConsensus-class}
\alias{featureGroupsConsensus}
\alias{comparison,featureGroups-method}
\alias{comparison}
\alias{plot,featureGroupsComparison,missing-method}
\alias{plotVenn,featureGroupsComparison-method}
\alias{plotUpSet,featureGroupsComparison-method}
\alias{plotChord,featureGroupsComparison-method}
\alias{consensus,featureGroupsComparison-method}
\alias{comparison,featureGroupsSet-method}
\alias{consensus,featureGroupsComparisonSet-method}
\title{Comparing feature groups}
\usage{
\S4method{comparison}{featureGroups}(..., groupAlgo, groupArgs = list(rtalign = FALSE))

\S4method{plot}{featureGroupsComparison,missing}(x, retMin = FALSE, ...)

\S4method{plotVenn}{featureGroupsComparison}(obj, which = NULL, ...)

\S4method{plotUpSet}{featureGroupsComparison}(obj, which = NULL, ...)

\S4method{plotChord}{featureGroupsComparison}(obj, addSelfLinks = FALSE, addRetMzPlots = TRUE, ...)

\S4method{consensus}{featureGroupsComparison}(
  obj,
  absMinAbundance = NULL,
  relMinAbundance = NULL,
  uniqueFrom = NULL,
  uniqueOuter = FALSE
)

\S4method{comparison}{featureGroupsSet}(..., groupAlgo, groupArgs = list(rtalign = FALSE))

\S4method{consensus}{featureGroupsComparisonSet}(obj, ...)
}
\arguments{
\item{\dots}{For \code{comparison}: \code{featureGroups} objects that should
  be compared. If the arguments are named (\emph{e.g.} \code{myGroups =
  fGroups}) then these are used for labelling, otherwise objects are
  automatically labelled by their \code{\link{algorithm}}.

  For \code{plot}, \code{plotVenn}, \code{plotChord}: further options passed
  to \code{plot}, \pkg{\link{VennDiagram}} plotting functions (\emph{e.g.}
  \code{\link{draw.pairwise.venn}}) and \code{\link{chordDiagram}}
  respectively.

  For \code{plotUpSet}: any further arguments passed to the \code{plotUpSet}
  method defined for \code{\link{featureGroups}}.}

\item{groupAlgo}{The \code{\link[=groupFeatures]{feature grouping
algorithm}} that should be used for grouping \emph{pseudo} features (see
details). Valid values are: \code{"xcms"}, \code{xcms3}, \code{kpic2} or \code{"openms"}.}

\item{groupArgs}{A \code{list} containing further parameters for
\code{\link[=groupFeatures]{feature grouping}}.}

\item{x, obj}{The \code{featureGroupsComparison} object.}

\item{retMin}{If \code{TRUE} retention times are plotted as minutes (seconds otherwise).}

\item{which}{A character vector specifying one or more labels of compared
feature groups. For \code{plotVenn}: if \code{NULL} then all compared
groups are used.}

\item{addSelfLinks}{If \code{TRUE} then 'self-links' are added which
represent non-shared data.}

\item{addRetMzPlots}{Set to \code{TRUE} to enable \emph{m/z} \emph{vs}
retention time scatter plots.}

\item{absMinAbundance, relMinAbundance}{Minimum absolute or relative
(\samp{0-1}) abundance across objects for a result to be kept. For
instance, \code{relMinAbundance=0.5} means that a result should be present
in at least half of the number of compared objects. Set to \samp{NULL} to
ignore and keep all results. Limits cannot be set when \code{uniqueFrom} is
not \code{NULL}.}

\item{uniqueFrom}{Set this argument to only retain feature groups that are unique
within one or more of the objects for which the consensus is made.
Selection is done by setting the value of \code{uniqueFrom} to a
\code{logical} (values are recycled), \code{numeric} (select by index) or a
\code{character} (as obtained with \code{algorithm(obj)}). For
\code{logical} and \code{numeric} values the order corresponds to the order
of the objects given for the consensus. Set to \code{NULL} to ignore.}

\item{uniqueOuter}{If \code{uniqueFrom} is not \code{NULL} and if
\code{uniqueOuter=TRUE}: only retain data that are also unique between
objects specified in \code{uniqueFrom}.}
}
\value{
\code{comparison} returns a \code{\link{featureGroupsComparison}}
  object.

\code{plotVenn} (invisibly) returns a list with the following fields: \itemize{
\item \code{gList} the \code{gList} object that was returned by
  the utilized \pkg{\link{VennDiagram}} plotting function.
\item \code{areas} The total area for each plotted group.
\item \code{intersectionCounts} The number of intersections between groups.
}

The order for the \code{areas} and \code{intersectionCounts} fields is the same as the parameter order
from the used plotting function (see \emph{e.g.} \code{\link{draw.pairwise.venn}} and
\code{\link{draw.triple.venn}}).

\code{consensus} returns a \code{\link{featureGroups}} object with a consensus from the compared feature
  groups.
}
\description{
Functionality to compare feature groups and make a consensus.
}
\details{
Feature groups objects originating from differing feature finding and/or
grouping algorithms (or their parameters) may be compared to assess their
output and generate a consensus.

The \code{comparison} method generates a
  \code{\link{featureGroupsComparison}} object from given feature groups
  objects, which in turn may be used for (visually) comparing presence of
  feature groups and generating a consensus. Internally, this function will
  collapse each feature groups object to \emph{pseudo} features objects by
  averaging their retention times, \emph{m/z} values and intensities, where
  each original feature groups object becomes an 'analysis'. All
  \emph{pseudo} features are then grouped using
  \link[=groupFeatures]{regular feature grouping algorithms} so that a
  comparison can be made.

\code{plot} generates an \emph{m/z} \emph{vs} retention time plot.

\code{plotVenn} plots a Venn diagram outlining unique and shared
  feature groups between up to five compared feature groups.

\code{plotUpSet} plots an UpSet diagram outlining unique and shared
  feature groups.

\code{plotChord} plots a chord diagram to visualize the distribution
  of feature groups.

\code{consensus} combines all compared feature groups and averages their retention, \emph{m/z} and intensity
  data. Not yet supported for \link[=sets-workflow]{sets workflows}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-mspeaklists.R
\name{getDefAvgPListParams}
\alias{getDefAvgPListParams}
\title{Parameters for averaging MS peak list data}
\usage{
getDefAvgPListParams(...)
}
\arguments{
\item{\dots}{Optional named arguments that override defaults.}
}
\value{
\code{getDefAvgPListParams} returns a \code{list} with the peak list averaging parameters.
}
\description{
Create parameter lists for averaging MS peak list data.
}
\details{
The parameters set used for averaging peak lists are set by the \code{avgFeatParams} and \code{avgFGroupParams}
arguments to \code{\link{generateMSPeakLists}} and its related algorithm specific functions. The parameters are
specified as a named \code{list} with the following values: \itemize{

\item \code{clusterMzWindow} \emph{m/z} window (in Da) used for clustering \emph{m/z} values when spectra are
averaged. For \code{method="hclust"} this corresponds to the cluster height, while for \code{method="distance"} this
value is used to find nearby masses (+/- window).  Too small windows will prevent clustering \emph{m/z} values (thus
erroneously treating equal masses along spectra as different), whereas too big windows may cluster unrelated
\emph{m/z} values from different or even the same spectrum together.

\item \code{topMost} Only retain this maximum number of MS peaks when generating averaged spectra. Lowering this
number may exclude more irrelevant (noisy) MS peaks and decrease processing time, whereas higher values may avoid
excluding lower intense MS peaks that may still be of interest.

\item \code{minIntensityPre} MS peaks with intensities below this value will be removed (applied prior to selection
by \code{topMost}) before averaging.

\item \code{minIntensityPost} MS peaks with intensities below this value will be removed after averaging.

\item \code{avgFun} Function that is used to calculate average \emph{m/z} values.

\item \code{method} Method used for producing averaged MS spectra. Valid values are \code{"hclust"}, used for
hierarchical clustering (using the \pkg{\link{fastcluster}} package), and \code{"distance"}, to use the between peak
distance. The latter method may reduces processing time and memory requirements, at the potential cost of reduced
accuracy.

\item \code{pruneMissingPrecursorMS} For MS data only: if \code{TRUE} then peak lists without a precursor peak are
removed. Note that even when this is set to \code{FALSE}, functionality that relies on MS (not MS/MS) peak lists
(\emph{e.g.} formulae calulcation) will still skip calculation if a precursor is not found.

\item \code{retainPrecursorMSMS} For MS/MS data only: if \code{TRUE} then always retain the precursor mass peak even
if is not amongst the \code{topMost} peaks. Note that MS precursor mass peaks are always kept. Furthermore, note that
precursor peaks in both MS and MS/MS data may still be removed by intensity thresholds (this is unlike the
\code{\link[=filter,MSPeakLists-method]{filter}} method function).

}

The \code{getDefAvgPListParams} function can be used to generate a default parameter list. The current defaults are:

\code{clusterMzWindow=0.005}; \code{topMost=50}; \code{minIntensityPre=500}; \code{minIntensityPost=500}; \code{avgFun=mean}; \code{method="hclust"}; \code{pruneMissingPrecursorMS=TRUE}; \code{retainPrecursorMSMS=TRUE}
}
\note{
With Bruker algorithms these parameters only control generation of feature groups averaged peak lists: how peak
  lists for features are generated is controlled by DataAnalysis.
}
\section{Source}{
 Averaging of mass spectra algorithms used by are based on the
  \href{https://github.com/zeehio/msProcess}{msProcess} R package (now archived on CRAN).
}

\references{
\addCitations{fastcluster}{1}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R, R/utils-exported.R
\name{analysis-information}
\alias{analysis-information}
\alias{generateAnalysisInfo}
\alias{generateAnalysisInfoFromEnviMass}
\title{Properties of sample analyses}
\usage{
generateAnalysisInfo(
  paths,
  groups = "",
  blanks = "",
  concs = NULL,
  formats = MSFileFormats()
)

generateAnalysisInfoFromEnviMass(path)
}
\arguments{
\item{paths}{A character vector containing one or more file paths that should be used for finding the analyses.}

\item{groups, blanks}{An (optional) character vector containing replicate groups and blanks, respectively (will be
recycled). If \code{groups} is an empty character string (\code{""}) the analysis name will be set as replicate
group.}

\item{concs}{An optional numeric vector containing concentration values for each analysis. Can be \code{NA} if
unknown. If the length of \code{concs} is less than the number of analyses the remainders will be set to \code{NA}.
Set to \code{NULL} to not include concentration data.}

\item{formats}{A character vector of analyses file types to consider. Analyses not present in these formats will be
ignored. For valid values see \code{\link{MSFileFormats}}.}

\item{path}{The path of the enviMass project.}
}
\description{
Properties for the sample analyses used in the workflow and utilities to automatically generate this information.
}
\details{
In \pkg{patRoon} a \emph{sample analysis}, or simply \emph{analysis}, refers to a single MS analysis file (sometimes
also called \emph{sample} or \emph{file}). The \emph{analysis information} summarizes several properties for the
analyses, and is used in various steps throughout the workflow, such as \code{\link{findFeatures}}, averaging
intensities of feature groups and blank subtraction. This information should be in a \code{data.frame}, with the
following columns:

\itemize{

\item \code{path} the full path to the directory of the analysis.

\item \code{analysis} the file name \strong{without} extension. Must be \strong{unique}, even if the \code{path} is
different.

\item \code{group} name of \emph{replicate group}. A replicate group is used to group analyses together that are
replicates of each other. Thus, the \code{group} column for all  analyses considered to be belonging to the same
replicate group should have an equal (but unique) value. Used for \emph{e.g.} averaging and
\code{\link[=filter,featureGroups-method]{filter}}.

\item \code{blank}: all analyses within this replicate group are used by the \code{featureGroups} method of
\code{\link[=filter,featureGroups-method]{filter}} for blank subtraction. Multiple entries can be entered by
separation with a comma.

\item \code{conc} a numeric value specifying the 'concentration' of the analysis. This can be actually any kind of
numeric value such as exposure time, dilution factor or anything else which may be used to form a linear
relationship.

}

Most workflows steps work with \file{mzXML} and \file{mzML} file formats. However, some algorithms only support
support one format (\emph{e.g.} \code{\link{findFeaturesOpenMS}}, \code{\link{findFeaturesEnviPick}}) or a
proprietary format (\code{\link{findFeaturesBruker}}). To mix such algorithms in the same workflow, the analyses
should be present in all required formats within the \emph{same} directory as specified by the \code{path} column.

Each analysis should only be specified \emph{once} in the analysis information, even if multiple file formats are
available. The \code{path} and \code{analysis} columns are internally used by \pkg{patRoon} to automatically find the
path of analysis files with the required format.

The \code{group} column is \emph{mandatory} and needs to be non-empty for each analysis. The \code{blank} column
should also be present, however, this may be empty (\code{""}) for analyses where no blank subtraction should occur.
The \code{conc} column is only required when obtaining regression information is desired with the
\code{\link[=as.data.table,featureGroups-method]{as.data.table}} method.

\code{generateAnalysisInfo} is an utility function that automatically generates a \code{data.frame} with
  analysis information. It scans the directories specified from the \code{paths} argument for analyses, and uses this
  to automatically fill in the \code{analysis} and \code{path} columns. Furthermore, this function also correctly
  handles analyses which are available in multiple formats.

\code{generateAnalysisInfoFromEnviMass} loads analysis information
  from an \pkg{enviMass} project. Note: this funtionality has only been
  tested with older versions of \pkg{enviMass}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_annotations.R
\docType{class}
\name{featureAnnotations-class}
\alias{featureAnnotations-class}
\alias{featureAnnotations}
\alias{annotations,featureAnnotations-method}
\alias{groupNames,featureAnnotations-method}
\alias{length,featureAnnotations-method}
\alias{[,featureAnnotations,ANY,missing,missing-method}
\alias{[[,featureAnnotations,ANY,missing-method}
\alias{$,featureAnnotations-method}
\alias{as.data.table,featureAnnotations-method}
\alias{delete,featureAnnotations-method}
\alias{filter,featureAnnotations-method}
\alias{plotVenn,featureAnnotations-method}
\alias{plotUpSet,featureAnnotations-method}
\title{Base feature annotations class}
\usage{
\S4method{annotations}{featureAnnotations}(obj)

\S4method{groupNames}{featureAnnotations}(obj)

\S4method{length}{featureAnnotations}(x)

\S4method{[}{featureAnnotations,ANY,missing,missing}(x, i, j, ..., drop = TRUE)

\S4method{[[}{featureAnnotations,ANY,missing}(x, i, j)

\S4method{$}{featureAnnotations}(x, name)

\S4method{as.data.table}{featureAnnotations}(
  x,
  fGroups = NULL,
  fragments = FALSE,
  countElements = NULL,
  countFragElements = NULL,
  OM = FALSE,
  normalizeScores = "none",
  excludeNormScores = defaultExclNormScores(x)
)

\S4method{delete}{featureAnnotations}(obj, i = NULL, j = NULL, ...)

\S4method{filter}{featureAnnotations}(
  obj,
  minExplainedPeaks = NULL,
  scoreLimits = NULL,
  elements = NULL,
  fragElements = NULL,
  lossElements = NULL,
  topMost = NULL,
  OM = FALSE,
  negate = FALSE
)

\S4method{plotVenn}{featureAnnotations}(obj, ..., labels = NULL, vennArgs = NULL)

\S4method{plotUpSet}{featureAnnotations}(
  obj,
  ...,
  labels = NULL,
  nsets = length(list(...)) + 1,
  nintersects = NA,
  upsetArgs = NULL
)
}
\arguments{
\item{obj, x}{\code{featureAnnotations} object to be accessed}

\item{i, j}{For \code{[}/\code{[[}: A numeric or character value which is used to select feature groups by
their index or name, respectively (for the order/names see \code{groupNames()}).\cr\cr For \code{[}: Can also be logical to perform logical selection
(similar to regular vectors). If missing all feature groups are selected.\cr\cr For \code{[[}: should be a scalar value.\cr\cr For \code{delete}: The data to remove from. \code{i} are the
feature groups as numeric index, logical or character, \code{j} the candidates as numeric indices (rows). If either is
\code{NULL} then data for all is removed. \code{j} may also be a function: it will be called for each 
feature group, with the annotation table (a \code{data.table}) as first argument, the feature group name as second argument, and any other arguments passed as
\code{\dots} to \code{delete}. The return value of this function specifies the candidate indices (rows) to be removed (specified as an \code{integer} or \code{logical} vector).}

\item{\dots}{For the \code{"["} operator: ignored.

  For \code{delete}: passed to the function specified as \code{j}.
  
  Others: Any further (and unique) \code{featureAnnotations} objects.}

\item{drop}{ignored.}

\item{name}{The feature group name (partially matched).}

\item{fGroups}{The \code{\link{featureGroups}} object that was used to
generate this object. If not \code{NULL} it is used to add feature group
information (retention and \emph{m/z} values).}

\item{fragments}{If \code{TRUE} then information on annotated fragments will be included. Automatically set to
\code{TRUE} if \code{countFragElements} is set.}

\item{countElements, countFragElements}{A \code{character} vector with elements that should be counted for each
candidate's formula. For instance, \code{c("C", "H")} adds columns for both carbon and hydrogen amounts of each
formula. Note that the neutral formula (\code{neutral_formula} column) is used to count elements of non-fragmented
formulae, whereas the charged formula of fragments (\code{ion_formula} column in \code{fragInfo} data) is used for
fragments. Set to \code{NULL} to not count any elements.}

\item{OM}{For \code{as.data.table}: if set to \code{TRUE} several columns with information relevant for organic
  matter (OM) characterization will be added (e.g. elemental ratios, classification). This will also make sure that
  \code{countElements} contains at least C, H, N, O, P and S.

  For \code{filter}: If \code{TRUE} then several filters are applied to exclude unlikely formula candidates present
  in organic matter (OM). See Source section for details.}

\item{normalizeScores}{A \code{character} that specifies how normalization of
annotation scorings occurs. Either
\code{"none"} (no normalization),
\code{"max"} (normalize to max value) or \code{"minmax"} (perform min-max
normalization). Note that normalization of negative scores (e.g. output by
\command{SIRIUS}) is always performed as min-max. Furthermore, currently
normalization for \code{compounds} takes the original min/max scoring
values into account when candidates were generated. Thus, for
\code{compounds} scoring, normalization is not affected when candidate
results were removed after they were generated (\emph{e.g.} by use of
\code{filter}).}

\item{excludeNormScores}{A
  \code{character} vector specifying any compound scoring names that
  should \emph{not} be normalized. Set to \code{NULL} to normalize all
  scorings. Note that whether any normalization occurs is set by the
  \code{excludeNormScores} argument.

  For \code{compounds}: By default \code{score} and
  \code{individualMoNAScore} are set to mimic the behavior of the
  \command{MetFrag} web interface.}

\item{minExplainedPeaks}{Minimum number of explained peaks. Set to \code{NULL} to ignore.}

\item{scoreLimits}{Filter results by their scores. Should be a named \code{list} that contains two-sized numeric
vectors with the minimum/maximum value of a score (use \code{-Inf}/\code{Inf} for no limits). The names of each
element should follow the name column of the table returned by \code{\link{formulaScorings}$name} and
\code{\link{compoundScorings}()$name}. For instance, \code{scoreLimits=list(numberPatents=c(10, Inf))} specifies
that \code{numberPatents} should be at least \samp{10}. Note that a result without a specified scoring is never
removed. If a score term exists multiple times, \emph{i.e.} due to a consensus, then a candidate is kept if at
least one of the terms falls within the range. Set to \code{NULL} to skip this filter.}

\item{elements}{Only retain candidate formulae (neutral form) that match a
given elemental restriction. The format of \code{elements} is a
\code{character} string with elements that should be present where each
element is followed by a valid amount or a range thereof. If no number is
specified then \samp{1} is assumed. For instance,
\code{elements="C1-10H2-20O0-2P"}, specifies that \samp{1-10}, \samp{2-20},
\samp{0-2} and \samp{1} carbon, hydrogen, oxygen and phosphorus atoms
should be present, respectively. When \code{length(elements)>1} formulas
are tested to follow at least one of the given elemental restrictions. For
instance, \code{elements=c("P", "S")} specifies that either one phosphorus
or one sulfur atom should be present. Set to \code{NULL} to ignore this
filter.}

\item{fragElements, lossElements}{Specifies elemental restrictions for
fragment or neutral loss formulae (charged form). Candidates are retained
if at least one of the fragment formulae follow (or not follow if
\code{negate=TRUE}) the given restrictions. See \code{elements} for the
used format.}

\item{topMost}{Only keep a maximum of \code{topMost} candidates with highest score (or least highest if
\code{negate=TRUE}). Set to \code{NULL} to ignore.}

\item{negate}{If \code{TRUE} then filters are applied in opposite manner.}

\item{labels}{A \code{character} with names to use for labelling. If \code{NULL} labels are automatically generated.}

\item{vennArgs}{A \code{list} with further arguments passed to \pkg{VennDiagram} plotting functions. Set to
\code{NULL} to ignore.}

\item{nsets, nintersects}{See \code{\link[UpSetR]{upset}}.}

\item{upsetArgs}{A list with any further arguments to be passed to \code{\link[UpSetR]{upset}}. Set to \code{NULL} to
ignore.}
}
\value{
\code{as.data.table} returns a \code{\link{data.table}}.

\code{delete} returns the object for which the specified data was removed.

\code{filter} returns a filtered \code{featureAnnotations} object.

\code{plotVenn} (invisibly) returns a list with the following fields: \itemize{
\item \code{gList} the \code{gList} object that was returned by
  the utilized \pkg{\link{VennDiagram}} plotting function.
\item \code{areas} The total area for each plotted group.
\item \code{intersectionCounts} The number of intersections between groups.
}

The order for the \code{areas} and \code{intersectionCounts} fields is the same as the parameter order
from the used plotting function (see \emph{e.g.} \code{\link{draw.pairwise.venn}} and
\code{\link{draw.triple.venn}}).
}
\description{
Holds information for all feature group annotations.
}
\details{
This class stores annotation data for feature groups, such as molecular formulae, SMILES identifiers, compound names
etc. The class of objects that are generated by formula and compound annotation (\code{\link{generateFormulas}} and
\code{\link{generateCompounds}}) are based on this class.
}
\section{Methods (by generic)}{
\itemize{
\item \code{annotations}: Accessor for the \code{groupAnnotations} slot.

\item \code{groupNames}: returns a \code{character} vector with the names of the
feature groups for which data is present in this object.

\item \code{length}: Obtain total number of candidates.

\item \code{[}: Subset on feature groups.

\item \code{[[}: Extracts annotation data for a feature group.

\item \code{$}: Extracts annotation data for a feature group.

\item \code{as.data.table}: Generates a table with all annotation data for each feature group and other
information such as element counts.

\item \code{delete}: Completely deletes specified annotations.

\item \code{filter}: Provides rule based filtering for feature group annotations. Useful to eliminate
unlikely candidates and speed up further processing.

\item \code{plotVenn}: plots a Venn diagram (using \pkg{\link{VennDiagram}}) outlining unique and shared
candidates of up to five different \code{featureAnnotations} objects.

\item \code{plotUpSet}: plots an UpSet diagram (using the \code{\link[UpSetR]{upset}} function) outlining
unique and shared formula candidates between different \code{featureAnnotations} objects.
}}

\section{Slots}{

\describe{
\item{\code{groupAnnotations}}{A \code{list} with for each annotated feature group a \code{data.table} with annotation data.
Use the \code{annotations} method for access.}

\item{\code{scoreTypes}}{A \code{character} with all the score types present in this object.}

\item{\code{scoreRanges}}{The minimum and maximum score values of all candidates for each feature group. Used for
normalization.}
}}

\section{Source}{
 Calculation of the aromaticity index (AI) and related double bond equivalents (DBE_AI) is performed
  as described in Koch 2015. Formula classification is performed by the rules described in Abdulla 2013. Filtering of
  OM related molecules is performed as described in Koch 2006 and Kujawinski 2006. (see references).
}

\section{S4 class hierarchy}{
 \itemize{   \item{\code{\link{workflowStep}}}   \itemize{     \item{\strong{\code{\link{featureAnnotations}}}}     \itemize{       \item{\code{\link{formulas}}}       \itemize{         \item{\code{\link{formulasConsensus}}}         \item{\code{\link{formulasSet}}}         \item{\code{\link{formulasUnset}}}       }       \item{\code{\link{compounds}}}       \itemize{         \item{\code{\link{compoundsConsensus}}}         \item{\code{\link{compoundsMF}}}         \item{\code{\link{compoundsSet}}}         \item{\code{\link{compoundsUnset}}}       }     }   } }
}

\references{
\insertRef{Koch2015}{patRoon} \cr\cr \insertRef{Abdulla2013}{patRoon} \cr\cr
  \insertRef{Koch2006}{patRoon} \cr\cr \insertRef{Kujawinski2006}{patRoon}

\insertRef{Conway2017}{patRoon} \cr\cr \insertRef{Lex2014}{patRoon}
}
\seealso{
\code{\link{formulas-class}} and \code{\link{compounds-class}}

The derived \code{\link{formulas}} and \code{\link{compounds}} classes.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TP-logic.R
\name{generateTPsLogic}
\alias{generateTPsLogic}
\alias{generateTPsLogic,featureGroups-method}
\alias{generateTPsLogic,featureGroupsSet-method}
\title{Obtain transformation products (TPs) with metabolic logic}
\usage{
\S4method{generateTPsLogic}{featureGroups}(fGroups, minMass = 40, adduct = NULL, transformations = NULL)

\S4method{generateTPsLogic}{featureGroupsSet}(fGroups, minMass = 40, transformations = NULL)
}
\arguments{
\item{fGroups}{A \code{\link{featureGroups}} object for which TPs should be calculated.}

\item{minMass}{A \code{numeric} that specifies the minimum mass of calculated TPs. If below this mass it will be
removed.}

\item{adduct}{An \code{\link{adduct}} object (or something that can be converted to it with \code{\link{as.adduct}}).
  Examples: \code{"[M-H]-"}, \code{"[M+Na]+"}. If the \code{featureGroups} object has
  adduct annotations then these are used if \code{adducts=NULL}.

  \setsWF The \code{adduct} argument is not supported for sets workflows, since the
  adduct annotations will then always be used.}

\item{transformations}{A \code{data.frame} with transformation reactions to be used for calculating the TPs (see
details below). If \code{NULL}, a default table from Schollee \emph{et al.} is used (see references).}
}
\value{
A \code{\link{transformationProducts}} (derived) object containing all generated TPs.
}
\description{
Automatically calculate potential transformation products with \emph{metabolic logic}.
}
\details{
This function uses metabolic logic to obtain transformation products. This function is called when calling \code{generateTPs} with
  \code{algorithm="logic"}.

With this algorithm TPs are predicted from common (environmental) chemical reactions, such as hydroxylation,
  demethylation etc. The generated TPs result from calculating the mass differences between a parent feature after it
  underwent the reaction. While this only results in little information on chemical properties of the TP, an
  advantage of this method is that it does not rely on structural information of the parent, which may be unknown in
  a full non-target analysis.
}
\section{Custom transformations}{
 The \code{transformations} argument to \code{generateTPsLogic} is used to specify
  custom rules to calculate transformation products. This should be a \code{data.frame} with the following columns:
  \itemize{

  \item \code{transformation} The name of the chemical transformation

  \item \code{add} The elements that are added by this reaction (\emph{e.g.} \code{"O"}).

  \item \code{sub} The elements that are removed by this reaction (\emph{e.g.} \code{"H2O"}).

  \item \code{retDir} The expected retention time direction relative to the parent (assuming a reversed phase like LC
  separation). Valid values are: \samp{-1} (elutes before the parent), \samp{1} (elutes after the parent) or \samp{0}
  (no significant change or unknown).

  }
}

\section{Source}{
 The algorithm of \code{generateTPsLogic} is directly based on the work done by Schollee \emph{et
  al.} (see references).
}

\references{
\insertRef{Scholle2015}{patRoon}
}
\seealso{
\code{\link{generateTPs}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features-kpic2.R
\name{getPICSet,features-method}
\alias{getPICSet,features-method}
\alias{getPICSet}
\alias{getPICSet,featuresKPIC2-method}
\title{Conversion to KPIC2 objects}
\usage{
\S4method{getPICSet}{features}(obj, loadRawData = TRUE)

\S4method{getPICSet}{featuresKPIC2}(obj, ...)
}
\arguments{
\item{obj}{The \code{features} object that should be converted.}

\item{loadRawData}{Set to \code{TRUE} if analyses are available as \code{mzXML} or \code{mzML} files. Otherwise MS
data is not loaded, and some dummy data (\emph{e.g.} file paths) is used in the returned object.}

\item{\dots}{Ignored}
}
\description{
Converts a \code{\link{features}} object to an \pkg{KPIC} object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/workflow-step-set.R
\docType{class}
\name{workflowStepSet-class}
\alias{workflowStepSet-class}
\alias{workflowStepSet}
\alias{setObjects,workflowStepSet-method}
\alias{sets,workflowStepSet-method}
\alias{show,workflowStepSet-method}
\title{(Virtual) base class for sets related workflow objects}
\usage{
\S4method{setObjects}{workflowStepSet}(obj)

\S4method{sets}{workflowStepSet}(obj)

\S4method{show}{workflowStepSet}(object)
}
\arguments{
\item{obj, object}{An object that is derived from \code{workflowStepSet}.}
}
\description{
This class is the base for many \link[=sets-workflow]{sets workflows} related classes. This class is virtual, and
therefore never created directly.
}
\details{
The most important purpose of this class is to hold data that is specific for a set. These \emph{set objects} are
typically objects with classes from a regular non-sets workflow (\emph{e.g.} \code{\link{components}},
\code{\link{compounds}}), and are used by the sets workflow object to \emph{e.g.} form a consensus. Since the set
objects may contain additional data, such as algorithm specific slots, it may in some cases be of interest to access
them directly with the \code{setObjects} method (described below).
}
\section{Methods (by generic)}{
\itemize{
\item \code{setObjects}: Accessor for the \code{setObjects} slot.

\item \code{sets}: Returns the names for each set in this object.

\item \code{show}: Shows summary information for this object.
}}

\section{Slots}{

\describe{
\item{\code{setObjects}}{A \code{list} with the \emph{set objects} (see the \verb{Details} section). The \code{list} is named
with the set names.}
}}

\section{S4 class hierarchy}{
 \itemize{   \item{\strong{\code{\link{workflowStepSet}}}}   \itemize{     \item{\code{\link{componentsSet}}}     \itemize{       \item{\code{\link{componentsNTSet}}}     }     \item{\code{\link{featureGroupsScreeningSet}}}     \item{\code{\link{compoundsSet}}}     \itemize{       \item{\code{\link{compoundsConsensusSet}}}     }     \item{\code{\link{formulasSet}}}     \itemize{       \item{\code{\link{formulasConsensusSet}}}     }     \item{\code{\link{MSPeakListsSet}}}   } }
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi-process.R
\name{executeMultiProcess}
\alias{executeMultiProcess}
\title{Simultaneous execution of system commands.}
\usage{
executeMultiProcess(
  commandQueue,
  finishHandler,
  timeoutHandler = function(...) TRUE,
  errorHandler = defMultiProcErrorHandler,
  prepareHandler = NULL,
  cacheName = NULL,
  setHash = NULL,
  procTimeout = NULL,
  printOutput = FALSE,
  printError = FALSE,
  logSubDir = NULL,
  showProgress = TRUE,
  waitTimeout = 50,
  batchSize = 1,
  delayBetweenProc = 0,
  method = NULL
)
}
\arguments{
\item{commandQueue}{A list with commands. Should contain \code{command}
(scalar string) and \code{args} (\code{character} vector). More user
defineds fields are allowed and useful to attach command information that
can be used in the finish, timeout and error handlers.}

\item{finishHandler}{A function that is called when a command has finished.
This function is typically used to process any results generated by the
command. The function is called right after spawning a new process, hence
processing results can occur while the next command is running in the
background. The function signature should be \code{function(cmd)} where
\code{cmd} is the queue data (from \code{commandQueue}) of the command that
has finished.}

\item{timeoutHandler}{A function that is called whenever a timeout for a
command occurs. Should return \code{TRUE} if execution of the command
should be retried. The function signature should be \code{function(cmd,
retries)} where \code{cmd} is the queue data for that command and
\code{retries} the number of times the command has been retried.}

\item{errorHandler}{Similar to \code{timeoutHandler}, but called whenever a
command has failed. The signature should be \code{function(cmd, exitStatus,
retries)}. The \code{exitStatus} argument is the exit code of the command
(may be \code{NA} in rare cases this is unknown). Other arguments are as
\code{timeoutHandler}. The return value should be as \code{timeoutHandler}
or a \code{character} with an error message which will be thrown with
\code{\link{stop}}.}

\item{prepareHandler}{A function that is called prior to execution of the
command. The function signature should be \code{function(cmd)} where
\code{cmd} is the queue data (from \code{commandQueue}) of the command to
be started. The return value must be (an updated) \code{cmd}.}

\item{cacheName, setHash}{Used for caching results. Set to \code{NULL} to
disable caching.}

\item{procTimeout}{The maximum time a process may consume before a timeout
occurs (in seconds). Set to \code{NULL} to disable timeouts. Ignored if
\code{patRoon.MP.method="future"}.}

\item{printOutput, printError}{Set to \code{TRUE} to print stdout/stderr
output to the console. Ignored if \option{patRoon.MP.method="future"}.}

\item{logSubDir}{The sub-directory used for log files. The final log file
path is constructed from \option{patRoon.MP.logPath}, \code{logSubDir} and
\code{logFile} set in the \code{commandQueue}.}

\item{showProgress}{Set to \code{TRUE} to display a progress bar. Ignored if
\option{patRoon.MP.method="future"}.}

\item{waitTimeout}{Number of milliseconds to wait before checking if a new
process should be spawned. Ignored if \option{patRoon.MP.method="future"}.}

\item{batchSize}{Number of commands that should be executed in sequence per
processes. See details. Ignored if \option{patRoon.MP.method="future"}.}

\item{delayBetweenProc}{Minimum number of milliseconds to wait before
spawning a new process. Might be needed to workaround errors. Ignored if
\option{patRoon.MP.method="future"}.}

\item{method}{Overrides \option{patRoon.MP.method} if not \code{NULL}.}
}
\description{
Execute a queue of system commands in parallel.
}
\details{
This function executes a given queue with system commands in parallel to
speed up computation. Commands are executed in the background using the
\pkg{processx} package. A configurable maximum amount of processes are
created to execute multiple commands in parallel.

Multiple commands may be executed in sequence that are launched from a single
parent process (as part of a batch script on Windows or combined with the
shell AND operator otherwise). Note that in this scenario still multiple
processes are spawned. Each of these processes will manage a chunk of the
command queue (size defined by \code{batchSize} argument). This approach is
typically suitable for fast running commands: the overhead of spawning a new
process for each command from R would in this case be significant enough to
loose most of the speedup otherwise gained with parallel execution. Note that
the actual batch size may be adjusted to ensure that a maximum number of
processes are running simultaneously.

Other functionalities of this function include timeout and error handling.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\docType{package}
\name{patRoon-package}
\alias{patRoon}
\alias{patRoon-package}
\title{Workflow solutions for mass-spectrometry based non-target analysis.}
\description{
\Sexpr[results=text,echo=FALSE]{packageDescription("patRoon", fields = "Description")}
}
\section{Package options}{


  The following package options (see \code{\link{options}}) can be set:

  \itemize{

  \item \code{patRoon.cache.mode}: A \code{character} setting the current caching mode: \code{"save"} and
  \code{"load"} will only save/load results to/from the cache, \code{"both"} (default) will do both and \code{"none"}
  to completely disable caching. This option can be changed anytime, which might be useful, for instance, to
  temporarily disable cached results before running a function.

  \item \code{patRoon.cache.fileName}: a \code{character} specifying the name of the cache file (default is
  \file{cache.sqlite}).

  \item \code{patRoon.MP.maxProcs}: The maximum number of processes that should be initiated in parallel. A good
  starting point is the number of physical cores, which is the default as detected by
  \code{\link[parallel]{detectCores}}. This option is only used when \option{patRoon.MP.method="classic"}.

  \item \code{patRoon.MP.method}: Either \code{"classic"} or \code{"future"}. The former is the default and uses
  \CRANpkg{processx} to execute multiple commands in parallel. When \code{"future"} the \code{\link{future.apply}}
  package is used for parallelization, which is especially useful for \emph{e.g.} cluster computing.

  \item \code{patRoon.MP.futureSched}: Sets the \code{future.scheduling} function argument for
  \code{\link{future_lapply}}. Only used if \option{patRoon.MP.method="future"}.

  \item \code{patRoon.MP.logPath}: The path used for logging of output from commands executed by multiprocess. Set to
  \code{FALSE} to disable logging.

  \item \code{patRoon.path.pwiz}: The path in which the \command{ProteoWizard} binaries are installed. If unset an
  attempt is made to find this directory from the Windows registry and \option{PATH} environment variable.

  \item \code{patRoon.path.GenForm}: The path to the \command{GenForm} executable. If not set (the default) the
  internal \code{GenForm} binary is used. Only set if you want to override the executable.

  \item \code{patRoon.path.MetFragCL}: The complete file path to the MetFrag CL \file{jar} file that \emph{must} be
  set when using \code{\link{generateCompoundsMetFrag}}. Example: \code{"C:/MetFrag2.4.2-CL.jar"}.

  \item \code{patRoon.path.MetFragCompTox}: The complete file path to the CompTox database \file{csv} file. See
  \code{\link{generateCompounds}} for more details.

  \item \code{patRoon.path.MetFragPubChemLite}: The complete file path to the PubChem database \file{csv} file. See
  \code{\link{generateCompounds}} for more details.

  \item \code{patRoon.path.SIRIUS}: The directory in which SIRIUS is installed. Unless the binaries can be located
  via the \option{PATH} environment variable, this \emph{must} be set when using \code{\link{generateFormulasSIRIUS}}
  or \code{\link{generateCompoundsSIRIUS}}. Example: \code{"C:/sirius-win64-3.5.1"}.

  \item \code{patRoon.path.OpenMS}: The path in which the \command{OpenMS} binaries are installed. Usually the
  location is added to the \option{PATH} environment variable when OpenMS is installed, in which case this option can
  be left empty.

  \item \code{patRoon.path.pngquant}: The path of the \command{pngquant} binary that is used when optimizing
  \file{.png} plots generated by \code{\link{reportHTML}} (with \code{optimizePng} set to \code{TRUE}). If the binary
  can be located through the \option{PATH} environment variable this option can remain empty. Note that some of the
  functionality of \code{reportHTML} only locates the binary through the \option{PATH} environment variable, hence,
  it is recommended to set up \option{PATH} instead.

  \item \code{patRoon.path.obabel}: The path in which the \command{OpenBabel} binaries are installed. Usually the
  location is added to the \option{PATH} environment variable when OpenBabel is installed, in which case this option
  can be left empty.

  \item \code{patRoon.path.BiotransFormer} The full file path to the \command{biotransformer} \file{.jar} command
  line utility. This needs to be set when \code{\link{generateTPsBioTransformer}} is used. For more details see
  \url{https://bitbucket.org/djoumbou/biotransformer/src/master}.

  }
}

\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/rickhelmus/patRoon}
  \item Report bugs at \url{https://github.com/rickhelmus/patRoon/issues}
}

}
\author{
\strong{Maintainer}: Rick Helmus \email{r.helmus@uva.nl} (\href{https://orcid.org/0000-0001-9401-3133}{ORCID})

Other contributors:
\itemize{
  \item Olaf Brock (\href{https://orcid.org/0000-0003-4727-8459}{ORCID}) [contributor]
  \item Vittorio Albergamo (\href{https://orcid.org/0000-0002-5347-1362}{ORCID}) [contributor]
  \item Andrea Brunner (\href{https://orcid.org/0000-0002-2801-1751}{ORCID}) [contributor]
  \item Emma Schymanski (\href{https://orcid.org/0000-0001-6868-8145}{ORCID}) [contributor]
  \item Bas van de Velde (\href{https://orcid.org/0000-0003-1292-3251}{ORCID}) [contributor]
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-exported.R
\name{featureQualityNames}
\alias{featureQualityNames}
\title{Returns chromatographic peak quality and score names for features and/or feature groups.}
\usage{
featureQualityNames(feat = TRUE, group = TRUE, scores = FALSE, totScore = TRUE)
}
\arguments{
\item{feat}{If \code{TRUE} then names specific to features are returned.}

\item{group}{If \code{TRUE} then names specific to groups are returned.}

\item{scores}{If \code{TRUE} the score names are returned, otherwise the quality names.}

\item{totScore}{If \code{TRUE} (and \code{scores=TRUE}) then the name of the total score is included.}
}
\description{
Returns chromatographic peak quality and score names for features and/or feature groups.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-compounds.R
\name{compoundScorings}
\alias{compoundScorings}
\title{Scorings terms for compound candidates}
\usage{
compoundScorings(
  algorithm = NULL,
  database = NULL,
  includeSuspectLists = TRUE,
  onlyDefault = FALSE,
  includeNoDB = TRUE
)
}
\arguments{
\item{algorithm}{The algorithm: \code{"metfrag"} or \code{"sirius"}. Set to \code{NULL} to return all scorings.}

\item{database}{The database for which results should be returned (\emph{e.g.} \code{"pubchem"}). Set to \code{NULL}
to return all scorings.}

\item{includeSuspectLists, onlyDefault, includeNoDB}{A logical specifying whether scoring terms related to suspect
lists, default scoring terms and non-database specific scoring terms should be included in the output,
respectively.}
}
\value{
A \code{data.frame} with information on which scoring terms are used, what their algorithm specific name is
  and other information such as to which database they apply and short remarks.
}
\description{
Returns an overview of scorings may be applied to rank candidate compounds.
}
\seealso{
generateCompounds
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert.R
\name{convertMSFiles}
\alias{convertMSFiles}
\alias{MSFileFormats}
\title{MS data conversion}
\usage{
MSFileFormats(algorithm = "pwiz", vendor = FALSE)

convertMSFiles(
  files = NULL,
  outPath = NULL,
  dirs = TRUE,
  anaInfo = NULL,
  from = NULL,
  to = "mzML",
  overWrite = FALSE,
  algorithm = "pwiz",
  centroid = algorithm != "openms",
  filters = NULL,
  extraOpts = NULL,
  PWizBatchSize = 1
)
}
\arguments{
\item{algorithm}{Either \code{"pwiz"} (implemented by \command{msConvert} of
ProteoWizard), \code{"openms"} (implemented by \command{FileConverter} of
OpenMS) or \code{"bruker"} (implemented by DataAnalysis).}

\item{vendor}{If \code{TRUE} only vendor formats are returned.}

\item{files, dirs}{The \code{files} argument should be a \code{character}
vector with input files. If \code{files} contains directories and
\code{dirs=TRUE} then files from these directories are also considered. An
alternative method to specify input files is by the \code{anaInfo}
argument. If the latter is specified \code{files} may be \code{NULL}.}

\item{outPath}{A character vector specifying directories that should be used
for the output. Will be re-cycled if necessary. If \code{NULL}, output
directories will be kept the same as the input directories.}

\item{anaInfo}{An \link[=analysis-information]{analysis info table} used to
retrieve input files. Either this argument or \code{files} (or both) should
be set (\emph{i.e.} not \code{NULL}).}

\item{from}{Input format (see below). These are used to find analyses when
\code{dirs=TRUE} or \code{anaInfo} is set.}

\item{to}{Output format: \code{"mzXML"} or \code{"mzML"}.}

\item{overWrite}{Should existing destination file be overwritten
(\code{TRUE}) or not (\code{FALSE})?}

\item{centroid}{Set to \code{TRUE} to enable centroiding (not supported if
\code{algorithm="openms"}). In addition, when \code{algorithm="pwiz"} the
value may be \code{"vendor"} to perform centroiding with the vendor
algorithm or \code{"cwt"} to use ProteoWizard's wavelet algorithm.}

\item{filters}{When \code{algorithm="pwiz"}: a \code{character} vector
specifying one or more filters. The elements of the specified vector are
directly passed to the \code{--filter} option (see
\href{http://proteowizard.sourceforge.net/tools/filters.html}{here})}

\item{extraOpts}{A \code{character} vector specifying any extra commandline
parameters passed to \command{msConvert} or \command{FileConverter}. Set to
\code{NULL} to ignore. For options: see
\href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FileConverter.html}{FileConverter}
 and
\href{http://proteowizard.sourceforge.net/tools/msconvert.html}{msConvert}.}

\item{PWizBatchSize}{When \code{algorithm="pwiz"}: the number of analyses to
process by a single call to \command{msConvert}. Usually a value of one is
most efficient. Set to zero to run all analyses all at once from a single
call.}
}
\description{
Conversion of MS analysis files between several open and closed data formats.
}
\details{
\code{MSFileFormats} returns a \code{character} with all supported
  input formats (see below).

\code{convertMSFiles} converts the data format of an analysis to
  another. It uses tools from
  \href{http://proteowizard.sourceforge.net/}{ProteoWizard}
  (\command{msConvert} command), \href{http://www.openms.de/}{OpenMS}
  (\command{FileConverter} command) or Bruker DataAnalysis to perform the
  conversion. Supported input and output formats include \file{mzXML},
  \file{.mzML} and several vendor formats, depending on which algorithm is
  used.
}
\section{Parallelization}{
 \code{convertMSFiles} (except if \code{algorithm="bruker"}) uses multiprocessing to parallelize
  computations. Please see the parallelization section in the handbook for
  more details and \link[=patRoon-package]{patRoon options} for configuration
  options.
}

\section{Conversion formats}{
 Possible output formats (\code{to} argument) are
  \code{mzXML} and \code{mzML}.

  Possible input formats (\code{from} argument) depend on the algorithm that
  was chosen and may include:

  \itemize{

  \item \code{thermo}: Thermo \file{.RAW} files (only
  \code{algorithm="pwiz"}).

  \item \code{bruker}: Bruker \file{.d}, \file{.yep}, \file{.baf} and
  \file{.fid} files (only \code{algorithm="pwiz"} or
  \code{algorithm="bruker"}).

  \item \code{agilent}: Agilent \file{.d} files (only
  \code{algorithm="pwiz"}).

  \item \code{ab}: AB Sciex \file{.wiff} files (only
  \code{algorithm="pwiz"}).

  \item \code{waters} Waters \file{.RAW} files (only
  \code{algorithm="pwiz"}).

  \item \code{mzXML}/\code{mzML}: Open format \file{.mzXML}/\file{.mzML}
  files (only \code{algorithm="pwiz"} or \code{algorithm="openms"}).

  }

  Note that the actual supported file formats of ProteoWizard depend on how
  it was installed (see
  \href{http://proteowizard.sourceforge.net/formats/index.html}{here}).
}

\examples{
\dontrun{
# Use FileConverter of OpenMS to convert between open mzXML/mzML format
convertMSFiles("standard-1.mzXML", to = "mzML", algorithm = "openms")

# Convert all Thermo .RAW files in the analyses/raw directory to mzML and
# store the files in analyses/mzml. During conversion files are centroided by
# the peakPicking filter and only MS 1 data is kept.
convertMSFiles("analyses/raw", "analyses/mzml", dirs = TRUE, from = "thermo",
               centroid = "vendor", filters = "msLevel 1")
}

}
\references{
\insertRef{Rst2016}{patRoon} \cr\cr
  \insertRef{Chambers2012}{patRoon}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-exported.R
\name{getFCParams}
\alias{getFCParams}
\title{Fold change calculation}
\usage{
getFCParams(rGroups, ...)
}
\arguments{
\item{rGroups}{A \code{character} vector with the names of the two replicate groups to be compared.}

\item{\dots}{Optional named arguments that override defaults.}
}
\description{
Fold change calculation
}
\details{
Fold change calculation can be used to easily identify significant changes between replicate groups. The
  calculation process is configured through a paramater list, which can be constructed with the \code{getFCParams}
  function. The parameter list has the following entries: \itemize{

  \item \code{rGroups} the name of the two replicate groups to compare (taken from the \code{rGroups} argument to
  \code{getFCParams}).

  \item \code{thresholdFC}: the threshold log FC for a feature group to be classified as increasing/decreasing.

  \item \code{thresholdPV}: the threshold log P for a feature group to be significantly different.

  \item \code{zeroMethod},\code{zeroValue}: how to handle zero values when calculating the FC: \code{add} adds an
  offset to zero values, \code{"fixed"} sets zero values to a fixed number and \code{"omit"} removes zero data. The
  number that is added/set by the former two options is defined by \code{zeroValue}.

  \item \code{PVTestFunc}: a function that is used to calculate P values (usually using \code{\link{t.test}}).

  \item \code{PVAdjFunc}: a function that is used to adjust P values (usually using \code{\link{p.adjust}})

  }
}
\seealso{
\code{\link{featureGroups-class}} and \link{feature-plotting}
}
\author{
The code to calculate and plot Fold change data was created by Bas van de Velde.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/workflow-step.R
\docType{class}
\name{workflowStep-class}
\alias{workflowStep-class}
\alias{workflowStep}
\alias{algorithm,workflowStep-method}
\alias{as.data.table,workflowStep-method}
\alias{as.data.frame,workflowStep-method}
\alias{show,workflowStep-method}
\title{(Virtual) Base class for all workflow objects.}
\usage{
\S4method{algorithm}{workflowStep}(obj)

\S4method{as.data.table}{workflowStep}(x, keep.rownames = FALSE, ...)

\S4method{as.data.frame}{workflowStep}(x, row.names = NULL, optional = FALSE, ...)

\S4method{show}{workflowStep}(object)
}
\arguments{
\item{obj, x, object}{An object (derived from) this class.}

\item{keep.rownames}{Ignored.}

\item{\dots}{Method specific arguments. Please see the documentation of the
derived classes.}

\item{row.names, optional}{Ignored.}
}
\description{
All workflow objects (\emph{e.g.} \code{\link{featureGroups}},
\code{\link{compounds}}, etc) are derived from this class. Objects from this
class are never created directly.
}
\section{Methods (by generic)}{
\itemize{
\item \code{algorithm}: Returns the algorithm that was used to generate an
object.

\item \code{as.data.table}: Summarizes the data in this object and returns this
as a \code{\link{data.table}}.

\item \code{as.data.frame}: This method simply calls \code{as.data.table} and
converts the result to a classic a \code{data.frame}.

\item \code{show}: Shows summary information for this object.
}}

\section{Slots}{

\describe{
\item{\code{algorithm}}{The algorithm that was used to generate this object. Use the
\code{algorithm} method for access.}
}}

\section{S4 class hierarchy}{
 \itemize{   \item{\strong{\code{\link{workflowStep}}}}   \itemize{     \item{\code{\link{transformationProducts}}}     \itemize{       \item{\code{\link{transformationProductsBT}}}       \item{\code{\link{transformationProductsLibrary}}}       \item{\code{\link{transformationProductsLogic}}}     }     \item{\code{\link{features}}}     \itemize{       \item{\code{\link{featuresSet}}}       \item{\code{\link{featuresUnset}}}       \item{\code{\link{featuresFromFeatGroups}}}       \item{\code{\link{featuresConsensus}}}       \item{\code{\link{featuresBruker}}}       \item{\code{\link{featuresEnviPick}}}       \item{\code{\link{featuresKPIC2}}}       \item{\code{\link{featuresOpenMS}}}       \item{\code{\link{featuresSAFD}}}       \item{\code{\link{featuresSIRIUS}}}       \item{\code{\link{featuresBrukerTASQ}}}       \item{\code{\link{featuresXCMS}}}       \item{\code{\link{featuresXCMS3}}}     }     \item{\code{\link{featureGroups}}}     \itemize{       \item{\code{\link{featureGroupsSet}}}       \itemize{         \item{\code{\link{featureGroupsScreeningSet}}}       }       \item{\code{\link{featureGroupsUnset}}}       \item{\code{\link{featureGroupsScreening}}}       \itemize{         \item{\code{\link{featureGroupsSetScreeningUnset}}}       }       \item{\code{\link{featureGroupsBruker}}}       \item{\code{\link{featureGroupsConsensus}}}       \item{\code{\link{featureGroupsEnviMass}}}       \item{\code{\link{featureGroupsKPIC2}}}       \item{\code{\link{featureGroupsOpenMS}}}       \item{\code{\link{featureGroupsSIRIUS}}}       \item{\code{\link{featureGroupsBrukerTASQ}}}       \item{\code{\link{featureGroupsXCMS}}}       \item{\code{\link{featureGroupsXCMS3}}}     }     \item{\code{\link{components}}}     \itemize{       \item{\code{\link{componentsCamera}}}       \item{\code{\link{componentsFeatures}}}       \itemize{         \item{\code{\link{componentsCliqueMS}}}         \item{\code{\link{componentsOpenMS}}}       }       \item{\code{\link{componentsClust}}}       \itemize{         \item{\code{\link{componentsIntClust}}}         \item{\code{\link{componentsSpecClust}}}       }       \item{\code{\link{componentsSet}}}       \itemize{         \item{\code{\link{componentsNTSet}}}       }       \item{\code{\link{componentsUnset}}}       \item{\code{\link{componentsNT}}}       \itemize{         \item{\code{\link{componentsNTUnset}}}       }       \item{\code{\link{componentsRC}}}       \item{\code{\link{componentsTPs}}}     }     \item{\code{\link{featureAnnotations}}}     \itemize{       \item{\code{\link{formulas}}}       \itemize{         \item{\code{\link{formulasConsensus}}}         \item{\code{\link{formulasSet}}}         \item{\code{\link{formulasUnset}}}       }       \item{\code{\link{compounds}}}       \itemize{         \item{\code{\link{compoundsConsensus}}}         \item{\code{\link{compoundsMF}}}         \item{\code{\link{compoundsSet}}}         \item{\code{\link{compoundsUnset}}}       }     }     \item{\code{\link{MSPeakLists}}}     \itemize{       \item{\code{\link{MSPeakListsSet}}}       \item{\code{\link{MSPeakListsUnset}}}     }   } }
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features-xcms3.R
\name{importFeaturesXCMS3}
\alias{importFeaturesXCMS3}
\title{Imports features from XCMS (new interface)}
\usage{
importFeaturesXCMS3(xdata, analysisInfo)
}
\arguments{
\item{xdata}{An \code{\link{XCMSnExp}} object.}

\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}
}
\value{
An object of a class which is derived from \code{\link{features}}.
}
\description{
Imports feature data generated from an existing \code{\link{XCMSnExp}} object generated by the \pkg{xcms} package.
}
\details{
This function imports data from XCMS3. This function is called when calling \code{importFeatures} with
  \code{type="xcms3"}.
}
\references{
\addCitations{xcms}{1} \cr\cr \addCitations{xcms}{2} \cr\cr \addCitations{xcms}{3}
}
\seealso{
\code{\link{importFeatures}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R, R/utils-bruker.R
\name{bruker-utils}
\alias{bruker-utils}
\alias{showDataAnalysis}
\alias{setDAMethod}
\alias{revertDAAnalyses}
\alias{recalibrarateDAFiles}
\alias{getDACalibrationError}
\alias{addDAEIC}
\alias{addAllDAEICs}
\title{Bruker DataAnalysis utilities}
\usage{
showDataAnalysis()

setDAMethod(anaInfo, method, close = TRUE)

revertDAAnalyses(anaInfo, close = TRUE, save = close)

recalibrarateDAFiles(anaInfo, close = TRUE, save = close)

getDACalibrationError(anaInfo)

addDAEIC(
  analysis,
  path,
  mz,
  mzWindow = 0.005,
  ctype = "EIC",
  mtype = "MS",
  polarity = "both",
  bgsubtr = FALSE,
  fragpath = "",
  name = NULL,
  hideDA = TRUE,
  close = FALSE,
  save = close
)

addAllDAEICs(
  fGroups,
  mzWindow = 0.005,
  ctype = "EIC",
  bgsubtr = FALSE,
  name = TRUE,
  onlyPresent = TRUE,
  hideDA = TRUE,
  close = FALSE,
  save = close
)
}
\arguments{
\item{anaInfo}{\link[=analysis-information]{Analysis info table}}

\item{method}{The full path of the DataAnalysis method.}

\item{close, save}{If \code{TRUE} then Bruker files are closed and saved after
processing with DataAnalysis, respectively. Setting \code{close=TRUE}
prevents that many analyses might be opened simultaneously in DataAnalysis,
which otherwise may use excessive memory or become slow. By default
\code{save} is \code{TRUE} when \code{close} is \code{TRUE}, which is
likely what you want as otherwise any processed data is lost.}

\item{analysis}{Analysis name (without file extension).}

\item{path}{path of the analysis.}

\item{mz}{\emph{m/z} (Da) value used for the chromatographic trace (if
applicable).}

\item{mzWindow}{\emph{m/z} window (in Da) used for the chromatographic trace
(if applicable).}

\item{ctype}{Type of the chromatographic trace. Valid options are:
\code{"EIC"} (extracted ion chromatogram), \code{"TIC"} (total ion
chromatogram, only for \code{addDAEIC}) and \code{"BPC"} (Base Peak
Chromatogram).}

\item{mtype}{MS filter for chromatographic trace. Valid values are:
\code{"all"}, \code{"MS"}, \code{"MSMS"}, \code{"allMSMS"} and
\code{"BBCID"}.}

\item{polarity}{Polarity filter for chromatographic trace. Valid values:
\code{"both"}, \code{"positive"} and \code{"negative"}.}

\item{bgsubtr}{If \code{TRUE} then background subtraction ('Spectral'
algorithm) will be performed.}

\item{fragpath}{Precursor \emph{m/z} used for MS/MS traces (\code{""} for
none).}

\item{name}{For \code{addDAEIC}: the name for the chromatographic trace. For
\code{addAllEICs}: \code{TRUE} to automatically set EIC names. Set to
\code{NULL} for none.}

\item{hideDA}{Hides DataAnalysis while adding the chromatographic trace
(faster).}

\item{fGroups}{The \code{\link{featureGroups}} object for which EICs should be made.}

\item{onlyPresent}{If \code{TRUE} then EICs are only generated for analyses
where the feature was detected.}
}
\value{
\code{getDACalibrationError} returns a \code{data.frame} with a
  column of all analyses (named \code{analysis}) and their mass error (named
  \code{error}).
}
\description{
Miscellaneous utility functions which interface with Bruker DataAnalysis
}
\details{
These functions communicate directly with Bruker DataAnalysis to provide
various functionality, such as calibrating and exporting data and adding
chromatographic traces. For this the \pkg{RDCOMClient} package is
required to be installed.

\code{showDataAnalysis} makes a hidden DataAnalysis window visible
  again. Most functions using DataAnalysis will hide the window during
  processing for efficiency reasons. If the window remains hidden
  (\emph{e.g.} because there was an error) this function can be used to make
  it visible again. This function can also be used to start DataAnalysis if
  it is not running yet.

\code{setDAMethod} Sets a given DataAnalysis method (\file{.m} file)
  to a set of analyses. \strong{NOTE}: as a workaround for a bug in
  DataAnalysis, this function will save(!), close and re-open any analyses
  that are already open prior to setting the new method. The \code{close}
  argument only controls whether the file should be closed after setting the
  method (files are always saved).

\code{revertDAAnalyses} Reverts a given set of analyses to their
  unprocessed raw state.

\code{recalibrarateDAFiles} Performs automatic mass recalibration of
  a given set of analyses. The current method settings for each analyses will
  be used.

\code{getDACalibrationError} is used to obtain the standard
  deviation of the current mass calibration (in ppm).

\code{addDAEIC} adds an Extracted Ion Chromatogram (EIC) or other
  chromatographic trace to a given analysis which can be used directly with
  DataAnalysis.

\code{addAllDAEICs} adds Extracted Ion Chromatograms (EICs) for all
  features within a \code{\link{featureGroups}} object.
}
\seealso{
\code{\link{analysis-information}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features-openms.R
\name{findFeaturesOpenMS}
\alias{findFeaturesOpenMS}
\title{Find features using OpenMS}
\usage{
findFeaturesOpenMS(
  analysisInfo,
  noiseThrInt = 1000,
  chromSNR = 3,
  chromFWHM = 5,
  mzPPM = 10,
  reEstimateMTSD = TRUE,
  traceTermCriterion = "sample_rate",
  traceTermOutliers = 5,
  minSampleRate = 0.5,
  minTraceLength = 3,
  maxTraceLength = -1,
  widthFiltering = "fixed",
  minFWHM = 1,
  maxFWHM = 30,
  traceSNRFiltering = FALSE,
  localRTRange = 10,
  localMZRange = 6.5,
  isotopeFilteringModel = "metabolites (5\% RMS)",
  MZScoring13C = FALSE,
  useSmoothedInts = TRUE,
  extraOpts = NULL,
  intSearchRTWindow = 3,
  useFFMIntensities = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}

\item{noiseThrInt}{Noise intensity threshold. Sets \code{algorithm:common:noise_threshold_int} option.}

\item{chromSNR}{Minimum S/N of a mass trace. Sets \code{algorithm:common:chrom_peak_snr} option.}

\item{chromFWHM}{Expected chromatographic peak width (in seconds). Sets \code{algorithm:common:chrom_fwhm} option.}

\item{mzPPM}{Allowed mass deviation (ppm) for trace detection. Sets \code{algorithm:mtd:mass_error_ppm}.}

\item{reEstimateMTSD}{If \code{TRUE} then enables dynamic re-estimation of m/z variance during mass trace collection
stage. Sets \code{algorithm:mtd:reestimate_mt_sd}.}

\item{traceTermCriterion, traceTermOutliers, minSampleRate}{Termination criterion for the extension of mass traces. See
\href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureFinderMetabo.html}{FeatureFinderMetabo}.
 Sets the \code{algorithm:mtd:trace_termination_criterion}, \code{algorithm:mtd:trace_termination_outliers} and
\code{algorithm:mtd:min_sample_rate} options, respectively.}

\item{minTraceLength, maxTraceLength}{Minimum/Maximum length of mass trace (seconds). Set negative value for maxlength
to disable maximum. Sets \code{algorithm:mtd:min_trace_length} and \code{algorithm:mtd:min_trace_length},
respectively.}

\item{widthFiltering, minFWHM, maxFWHM}{Enable filtering of unlikely peak widths. See
\href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureFinderMetabo.html}{FeatureFinderMetabo}.
 Sets \code{algorithm:epd:width_filtering}, \code{algorithm:epd:min_fwhm} and \code{algorithm:epd:max_fwhm},
respectively.}

\item{traceSNRFiltering}{If \code{TRUE} then apply post-filtering by signal-to-noise ratio after smoothing. Sets the
\code{algorithm:epd:masstrace_snr_filtering} option.}

\item{localRTRange, localMZRange}{Retention/MZ range where to look for coeluting/isotopic mass traces. Sets the
\code{algorithm:ffm:local_rt_range} and \code{algorithm:ffm:local_mz_range} options, respectively.}

\item{isotopeFilteringModel}{Remove/score candidate assemblies based on isotope intensities. See
\href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureFinderMetabo.html}{FeatureFinderMetabo}.
 Sets the \code{algorithm:ffm:isotope_filtering_model} option.}

\item{MZScoring13C}{Use the 13C isotope as the expected shift for isotope mass traces. See
\href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureFinderMetabo.html}{FeatureFinderMetabo}.
 Sets \code{algorithm:ffm:mz_scoring_13C}.}

\item{useSmoothedInts}{If \code{TRUE} then use LOWESS intensities instead of raw intensities. Sets the
\code{algorithm:ffm:use_smoothed_intensities} option.}

\item{extraOpts}{Named \code{list} containing extra options that will be passed to \command{FeatureFinderMetabo}. Any
options specified here will override any of the above. Example:
\code{extraOpts=list("-algorithm:common:noise_threshold_int"=1000)} (corresponds to setting
\code{noiseThrInt=1000}). Set to \code{NULL} to ignore.}

\item{intSearchRTWindow}{Retention time window (in seconds, +/- feature retention time) that is used to find the
closest data point to the retention time to obtain the intensity of a feature (this is needed since OpenMS does not
provide this data).}

\item{useFFMIntensities}{If \code{TRUE} then peak intensities are directly loaded from \command{FeatureFinderMetabo}
output. Otherwise, intensities are loaded afterwards from the input \file{mzML} files, which is potentially much
slower, especially with many analyses files. However, \code{useFFMIntensities=TRUE} is still somewhat experimental,
may be less accurate and requires a recent version of \command{OpenMS} (>=2.7).}

\item{verbose}{If set to \code{FALSE} then no text output is shown.}
}
\value{
An object of a class which is derived from \code{\link{features}}.
}
\description{
uses the
\href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_FeatureFinderMetabo.html}{FeatureFinderMetabo}
TOPP tool (see \url{http://www.openms.de}) to find features.
}
\details{
This function uses OpenMS to automatically find features. This function is called when calling \code{findFeatures} with
  \code{algorithm="openms"}.

This functionality has been tested with OpenMS version >= 2.0. Please make sure it is installed and its
  binaries are added to the PATH environment variable or the \code{patRoon.path.OpenMS} option is set.

  The file format of analyses must be \file{mzML}.

The input MS data files need to be centroided. The \code{\link{convertMSFiles}} function can be used to
  centroid data.
}
\section{Parallelization}{
 \code{findFeaturesOpenMS} uses multiprocessing to parallelize
  computations. Please see the parallelization section in the handbook for
  more details and \link[=patRoon-package]{patRoon options} for configuration
  options.

 Note that for caching purposes, the analyses files must always exist on the local host
  computer, even if it is not participating in computations.
}

\references{
\insertRef{Rst2016}{patRoon} \cr\cr
  \href{https://pugixml.org/}{pugixml} (via \href{http://www.rcpp.org/}{Rcpp}) is used to process OpenMS XML output. \cr\cr
  \addCitations{Rcpp}{1} \cr\cr
  \addCitations{Rcpp}{2} \cr\cr
  \addCitations{Rcpp}{3}
}
\seealso{
\code{\link{findFeatures}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups-envimass.R
\name{importFeatureGroupsEnviMass}
\alias{importFeatureGroupsEnviMass}
\title{Imports feature groups from enviMass}
\usage{
importFeatureGroupsEnviMass(path, feat, positive)
}
\arguments{
\item{path}{The path of the enviMass project.}

\item{feat}{The \code{\link{features}} object obtained with \code{\link{importFeaturesEnviMass}}.}

\item{positive}{Whether data from positive (\code{TRUE}) or negative (\code{FALSE}) should be loaded.}
}
\value{
An object of a class which is derived from \code{\link{featureGroups}}.

The \code{featuresSet} method (for \link[=sets-workflow]{sets workflows}) returns a
  \code{\link{featureGroupsSet}} object.
}
\description{
Imports a 'profiles' produced by \pkg{enviMass}.
}
\details{
This function imports data from enviMass. This function is called when calling \code{importFeatureGroups} with
  \code{type="envimass"}.

This function \emph{only} imports 'raw' profiles, \emph{not} any results from further componentization steps
  performed in \pkg{enviMass}. Furthermore, this functionality has only been tested with older versions of
  \pkg{enviMass}. Finally, please note that this function only supports features imported by
  \code{\link{importFeaturesEnviMass}} (obviously, the same project should be used for both importing functions).
}
\seealso{
\code{\link{importFeatureGroups}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mspeaklists.R
\name{generateMSPeakLists}
\alias{generateMSPeakLists}
\alias{generateMSPeakLists,featureGroups-method}
\title{Generation of MS Peak Lists}
\usage{
\S4method{generateMSPeakLists}{featureGroups}(fGroups, algorithm, ...)
}
\arguments{
\item{fGroups}{The \code{\link{featureGroups}} object from which MS peak lists should be extracted.}

\item{algorithm}{A character string describing the algorithm that should be
used: \code{"bruker"}, \code{"brukerfmf"}, \code{"mzr"}}

\item{\dots}{Any parameters to be passed to the selected MS peak lists generation algorithm.}
}
\value{
A \code{\link{MSPeakLists}} object.
}
\description{
Functionality to convert MS and MS/MS data into MS peak lists.
}
\details{
Formula calculation and identification tools rely on mass spectra that belong to features of interest. For
processing, MS (and MS/MS) spectra are typically reduced to a table with a column containing measured \emph{m/z}
values and a column containing their intensities. These 'MS peak lists' can then be used for
\link[=generateFormulas]{formula generation} and \link[=generateCompounds]{compound generation}.

MS and MS/MS peak lists are first generated for all features (or a subset, if the \code{topMost} argument is set).
During this step multiple spectra over the feature elution profile are averaged. Subsequently, peak lists will be
generated for each feature group by averaging peak lists of the features within the group. Functionality that uses
peak lists will either use data from individual features or from group averaged peak lists. For instance, the former
may be used by formulae calculation, while compound identification and plotting functionality typically uses group
averaged peak lists.

\code{generateMSPeakLists} is a generic function that will generateMSPeakLists by one of the supported algorithms. The actual
  functionality is provided by algorithm specific functions such as \code{generateMSPeakListsMzR} and \code{generateMSPeakListsDA}. While these
  functions may be called directly, \code{generateMSPeakLists} provides a generic interface and is therefore usually preferred.
}
\note{
In most cases it will be necessary to centroid your MS input files. The only exception is \command{Bruker},
  however, you will still need centroided \file{mzXML}/\file{mzML} files for \emph{e.g.} plotting chromatograms. In
  this case the centroided MS files should be stored in the same directory as the raw \command{Bruker} \file{.d}
  files. The \code{\link{convertMSFiles}} function can be used to centroid data.
}
\section{Sets workflows}{
 With a \link[=sets-workflow]{sets workflow}, the feature group averaged peak lists are made
  per set. This is important, because for averaging peak lists cannot be mixed, for instance, when different
  ionization modes were used to generate the sets. The group averaged peaklists are then simply combined and labelled
  in the final peak lists. However, please note that annotation and other functionality typically uses only the set
  specific peak lists, as this functionality cannot work with mixed peak lists.
}

\seealso{
The \code{\link{MSPeakLists}} output class and its methods and the algorithm specific functions:
  \code{\link{generateMSPeakListsDA}}, \code{\link{generateMSPeakListsDAFMF}}, \code{\link{generateMSPeakListsMzR}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TP.R, R/TP-logic.R
\docType{class}
\name{transformationProducts-class}
\alias{transformationProducts-class}
\alias{transformationProducts}
\alias{parents,transformationProducts-method}
\alias{parents}
\alias{products,transformationProducts-method}
\alias{products}
\alias{length,transformationProducts-method}
\alias{names,transformationProducts-method}
\alias{show,transformationProducts-method}
\alias{[,transformationProducts,ANY,missing,missing-method}
\alias{[[,transformationProducts,ANY,missing-method}
\alias{$,transformationProducts-method}
\alias{as.data.table,transformationProducts-method}
\alias{convertToSuspects,transformationProducts-method}
\alias{transformationProductsLogic-class}
\alias{transformationProductsLogic}
\title{Base transformation products (TP) class}
\usage{
\S4method{parents}{transformationProducts}(TPs)

\S4method{products}{transformationProducts}(TPs)

\S4method{length}{transformationProducts}(x)

\S4method{names}{transformationProducts}(x)

\S4method{show}{transformationProducts}(object)

\S4method{[}{transformationProducts,ANY,missing,missing}(x, i, j, ..., drop = TRUE)

\S4method{[[}{transformationProducts,ANY,missing}(x, i, j)

\S4method{$}{transformationProducts}(x, name)

\S4method{as.data.table}{transformationProducts}(x)

\S4method{convertToSuspects}{transformationProducts}(TPs, includeParents = FALSE)
}
\arguments{
\item{TPs, x, object}{\code{transformationProducts} object to be accessed}

\item{i}{For \code{[}/\code{[[}: A numeric or character value which is used to select parents by
their index or name, respectively (for the order/names see \code{names()}).\cr\cr For \code{[}: Can also be logical to perform logical selection
(similar to regular vectors). If missing all parents are selected.\cr\cr For \code{[[}: should be a scalar value.}

\item{\dots}{Unused.}

\item{drop, j}{ignored.}

\item{name}{The parent name (partially matched).}

\item{includeParents}{If \code{TRUE} then parents are also included in the returned suspect list.}
}
\description{
Holds information for all TPs for a set of parents.
}
\details{
This class holds all generated data for transformation products for a set of parents. The class is \code{virtual} and
derived objects are created by \link[=generateTPs]{TP generators}.

The TP data in objects from this class include a \code{retDir} column. These are \code{numeric} values that hint what
the the chromatographic retention order of a TP might be compared to its parent: a value of \samp{-1} means it will
elute earlier, \samp{1} it will elute later and \samp{0} that there is no significant difference or the direction is
unknown. These values are based on a typical reversed phase separation. When \code{\link{generateTPsBioTransformer}}
or \code{\link{generateTPsLibrary}} was used to generate the data, the \code{retDir} values are based on calculated
\code{log P} values of the parent and its TPs.
}
\section{Methods (by generic)}{
\itemize{
\item \code{parents}: Accessor method for the \code{parents} slot of a
\code{transformationProducts} class.

\item \code{products}: Accessor method for the \code{products} slot.

\item \code{length}: Obtain total number of transformation products.

\item \code{names}: Obtain the names of all parents in this object.

\item \code{show}: Show summary information for this object.

\item \code{[}: Subset on parents.

\item \code{[[}: Extracts a table with TPs for a parent.

\item \code{$}: Extracts a table with TPs for a parent.

\item \code{as.data.table}: Returns all TP data in a table.

\item \code{convertToSuspects}: Converts this object to a suspect list that can be used as input for
\code{\link{screenSuspects}}.
}}

\section{Slots}{

\describe{
\item{\code{parents}}{A \code{\link{data.table}} with metadata for all parents that have TPs in this object. Use the
\code{parents} method for access.}

\item{\code{products}}{A \code{list} with \code{\link{data.table}} entries with TP information for each parent. Use the
\code{products} method for access.}
}}

\section{S4 class hierarchy}{
 \itemize{   \item{\code{\link{workflowStep}}}   \itemize{     \item{\strong{\code{\link{transformationProducts}}}}     \itemize{       \item{\code{\link{transformationProductsBT}}}       \item{\code{\link{transformationProductsLibrary}}}       \item{\code{\link{transformationProductsLogic}}}     }   } }
}

\seealso{
Derived classes \code{\link{transformationProductsBT}} and \code{\link{transformationProductsLibrary}} for
  specific algorithm methods and \code{\link{generateTPs}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups-xcms3.R
\name{importFeatureGroupsXCMS3}
\alias{importFeatureGroupsXCMS3}
\title{Imports feature groups from XCMS (new interface)}
\usage{
importFeatureGroupsXCMS3(xdata, analysisInfo)
}
\arguments{
\item{xdata}{An \code{\link{XCMSnExp}} object.}

\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}
}
\value{
An object of a class which is derived from \code{\link{featureGroups}}.

The \code{featuresSet} method (for \link[=sets-workflow]{sets workflows}) returns a
  \code{\link{featureGroupsSet}} object.
}
\description{
Imports grouped features from a \code{\link{XCMSnExp}} object from the \pkg{xcms} package.
}
\references{
\addCitations{xcms}{1} \cr\cr \addCitations{xcms}{2} \cr\cr \addCitations{xcms}{3}
}
\seealso{
\code{\link{groupFeatures}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/utils-xcms.R
\name{getXCMSSet}
\alias{getXCMSSet}
\alias{getXCMSnExp}
\alias{getXCMSSet,features-method}
\alias{getXCMSSet,featuresXCMS-method}
\alias{getXCMSSet,featureGroups-method}
\alias{getXCMSSet,featureGroupsXCMS-method}
\alias{getXCMSSet,featuresSet-method}
\alias{getXCMSSet,featureGroupsSet-method}
\alias{getXCMSnExp,features-method}
\alias{getXCMSnExp,featuresXCMS3-method}
\alias{getXCMSnExp,featureGroups-method}
\alias{getXCMSnExp,featureGroupsXCMS3-method}
\alias{getXCMSnExp,featuresSet-method}
\alias{getXCMSnExp,featureGroupsSet-method}
\title{Conversion to XCMS objects}
\usage{
getXCMSSet(obj, verbose = TRUE, ...)

getXCMSnExp(obj, verbose = TRUE, ...)

\S4method{getXCMSSet}{features}(obj, verbose, loadRawData)

\S4method{getXCMSSet}{featuresXCMS}(obj, verbose = TRUE, ...)

\S4method{getXCMSSet}{featureGroups}(obj, verbose, loadRawData)

\S4method{getXCMSSet}{featureGroupsXCMS}(obj, verbose, loadRawData)

\S4method{getXCMSSet}{featuresSet}(obj, ..., set)

\S4method{getXCMSSet}{featureGroupsSet}(obj, ..., set)

\S4method{getXCMSnExp}{features}(obj, verbose, loadRawData)

\S4method{getXCMSnExp}{featuresXCMS3}(obj, verbose = TRUE, ...)

\S4method{getXCMSnExp}{featureGroups}(obj, verbose, loadRawData)

\S4method{getXCMSnExp}{featureGroupsXCMS3}(obj, verbose, loadRawData)

\S4method{getXCMSnExp}{featuresSet}(obj, ..., set)

\S4method{getXCMSnExp}{featureGroupsSet}(obj, ..., set)
}
\arguments{
\item{obj}{The object that should be converted.}

\item{verbose}{If \code{FALSE} then no text output is shown.}

\item{\dots}{\setsWF Further arguments passed to non-sets method.

  Otherwise ignored.}

\item{loadRawData}{Set to \code{TRUE} if analyses are available as \code{mzXML} or \code{mzML} files. Otherwise MS
data is not loaded, and some dummy data (\emph{e.g.} file paths) is used in the returned object.}

\item{set}{\setsWF The name of the set to be exported.}
}
\description{
Converts a \code{\link{features}} or \code{\link{featureGroups}} object to an \code{\link{xcmsSet}} or
\code{\link{XCMSnExp}} object.
}
\section{Sets workflows}{
 In a \link[=sets-workflow]{sets workflow}, \code{\link{unset}} is used to convert the
  feature (group) data before the object is exported.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups.R
\name{importFeatureGroups}
\alias{importFeatureGroups}
\title{Import feature groups from files}
\usage{
importFeatureGroups(path, type, ...)
}
\arguments{
\item{path}{The path that should be used for importing. See the algorithm specific functions for more details.}

\item{type}{Which file type should be imported: \code{"brukerpa"} (Bruker ProfileAnalysis), \code{"brukertasq"}
(Bruker TASQ) or \code{"envimass"} (\pkg{enviMass}).}

\item{\dots}{Further arguments passed to the selected import algorithm function.}
}
\value{
An object of a class which is derived from \code{\link{featureGroups}}.

The \code{featuresSet} method (for \link[=sets-workflow]{sets workflows}) returns a
  \code{\link{featureGroupsSet}} object.
}
\description{
Generic function to import feature groups produced by other software from files.
}
\details{
\code{importFeatureGroups} is a generic function that will import feature groups from files by one of the supported algorithms. The actual
  functionality is provided by algorithm specific functions such as \code{importFeatureGroupsBrukerTASQ} and \code{importFeatureGroupsBrukerPA}. While these
  functions may be called directly, \code{importFeatureGroups} provides a generic interface and is therefore usually preferred.
}
\seealso{
The \code{\link{featureGroups}} output class and its methods and the algorithm specific functions:
  \code{\link{importFeatureGroupsBrukerPA}}, \code{\link{importFeatureGroupsBrukerTASQ}}, \code{\link{importFeatureGroupsEnviMass}}

\code{\link{groupFeatures}} to group features. Other import functions:
  \code{\link{importFeatureGroupsXCMS}}, \code{\link{importFeatureGroupsXCMS3}} and
  \code{\link{importFeatureGroupsKPIC2}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features-kpic2.R
\name{findFeaturesKPIC2}
\alias{findFeaturesKPIC2}
\title{Find features using KPIC2}
\usage{
findFeaturesKPIC2(
  analysisInfo,
  kmeans = TRUE,
  level = 1000,
  ...,
  parallel = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}

\item{kmeans}{If \code{TRUE} then \code{\link[KPIC]{getPIC.kmeans}} is used to obtain PICs, otherwise it is
\code{\link[KPIC]{getPIC}}.}

\item{level}{Passed to \code{\link[KPIC]{getPIC}} or \code{\link[KPIC]{getPIC.kmeans}}}

\item{\dots}{Further parameters passed to \code{\link[KPIC]{getPIC}}/\code{\link[KPIC]{getPIC.kmeans}}}

\item{parallel}{If set to \code{TRUE} then code is executed in parallel through the \CRANpkg{futures} package. Please
see the parallelization section in the handbook for more details.}

\item{verbose}{If set to \code{FALSE} then no text output is shown.}
}
\value{
An object of a class which is derived from \code{\link{features}}.
}
\description{
Uses the \href{https://github.com/hcji/KPIC2}{KPIC2} \R package to extract features.
}
\details{
This function uses KPIC2 to automatically find features. This function is called when calling \code{findFeatures} with
  \code{algorithm="kpic2"}.

The MS files should be in the \code{mzML} or \code{mzXML} format.

The input MS data files need to be centroided. The \code{\link{convertMSFiles}} function can be used to
  centroid data.
}
\references{
\insertRef{Ji2017}{patRoon}
}
\seealso{
\code{\link{findFeatures}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features.R, R/features-set.R,
%   R/features-bruker.R, R/features-envipick.R, R/features-kpic2.R,
%   R/features-openms.R, R/features-safd.R, R/features-sirius.R,
%   R/features-xcms.R, R/features-xcms3.R
\docType{class}
\name{features-class}
\alias{features-class}
\alias{features}
\alias{length,features-method}
\alias{show,features-method}
\alias{featureTable,features-method}
\alias{analysisInfo,features-method}
\alias{analyses,features-method}
\alias{replicateGroups,features-method}
\alias{as.data.table,features-method}
\alias{filter,features-method}
\alias{[,features,ANY,missing,missing-method}
\alias{[[,features,ANY,missing-method}
\alias{$,features-method}
\alias{delete,features-method}
\alias{calculatePeakQualities,features-method}
\alias{featuresSet-class}
\alias{featuresSet}
\alias{sets,featuresSet-method}
\alias{show,featuresSet-method}
\alias{as.data.table,featuresSet-method}
\alias{[,featuresSet,ANY,missing,missing-method}
\alias{filter,featuresSet-method}
\alias{featuresUnset-class}
\alias{featuresUnset}
\alias{unset,featuresSet-method}
\alias{featuresBruker-class}
\alias{featuresBruker}
\alias{featuresEnviPick-class}
\alias{featuresEnviPick}
\alias{featuresKPIC2-class}
\alias{featuresKPIC2}
\alias{delete,featuresKPIC2-method}
\alias{featuresOpenMS-class}
\alias{featuresOpenMS}
\alias{featuresSAFD-class}
\alias{featuresSAFD}
\alias{featuresSIRIUS-class}
\alias{featuresSIRIUS}
\alias{featuresXCMS-class}
\alias{featuresXCMS}
\alias{delete,featuresXCMS-method}
\alias{featuresXCMS3-class}
\alias{featuresXCMS3}
\alias{delete,featuresXCMS3-method}
\title{Base features class}
\usage{
\S4method{length}{features}(x)

\S4method{show}{features}(object)

\S4method{featureTable}{features}(obj)

\S4method{analysisInfo}{features}(obj)

\S4method{analyses}{features}(obj)

\S4method{replicateGroups}{features}(obj)

\S4method{as.data.table}{features}(x)

\S4method{filter}{features}(
  obj,
  absMinIntensity = NULL,
  relMinIntensity = NULL,
  retentionRange = NULL,
  mzRange = NULL,
  mzDefectRange = NULL,
  chromWidthRange = NULL,
  qualityRange = NULL,
  negate = FALSE
)

\S4method{[}{features,ANY,missing,missing}(x, i, j, ..., drop = TRUE)

\S4method{[[}{features,ANY,missing}(x, i)

\S4method{$}{features}(x, name)

\S4method{delete}{features}(obj, i = NULL, j = NULL, ...)

\S4method{calculatePeakQualities}{features}(obj, weights, flatnessFactor, parallel = TRUE)

\S4method{sets}{featuresSet}(obj)

\S4method{show}{featuresSet}(object)

\S4method{as.data.table}{featuresSet}(x)

\S4method{[}{featuresSet,ANY,missing,missing}(x, i, ..., sets = NULL, drop = TRUE)

\S4method{filter}{featuresSet}(obj, ..., negate = FALSE, sets = NULL)

\S4method{unset}{featuresSet}(obj, set)

\S4method{delete}{featuresKPIC2}(obj, i = NULL, j = NULL, ...)

\S4method{delete}{featuresXCMS}(obj, i = NULL, j = NULL, ...)

\S4method{delete}{featuresXCMS3}(obj, i = NULL, j = NULL, ...)
}
\arguments{
\item{obj, x, object}{\code{features} object to be accessed}

\item{absMinIntensity, relMinIntensity}{Minimum absolute/relative intensity for features to be kept. The relative
intensity is determined from the feature with highest intensity (within the same analysis). Set to \samp{0} or \code{NULL} to skip this step.}

\item{retentionRange, mzRange, mzDefectRange, chromWidthRange}{Range of retention time (in seconds), \emph{m/z}, mass
defect (defined as the decimal part of \emph{m/z} values) or chromatographic peak width (in seconds), respectively.
Features outside this range will be removed. Should be a numeric vector with length of two containing the min/max
values. The maximum can be \code{Inf} to specify no maximum range. Set to \code{NULL} to skip this step.}

\item{qualityRange}{Used to filter features by their peak qualities/scores
(see \code{\link{calculatePeakQualities}}). Should be a named \code{list} with min/max ranges for each
quality/score to be filtered (the \code{\link{featureQualityNames}} function can be used to obtain valid names).
Example: \code{qualityRange=list(ModalityScore=c(0.3, Inf),
SymmetryScore=c(0.5, Inf))}. Set to \code{NULL} to ignore.}

\item{negate}{If set to \code{TRUE} then filtering operations are performed in opposite manner.}

\item{i, j}{For \code{[}/\code{[[}: A numeric or character value which is used to select analyses by
their index or name, respectively (for the order/names see \code{analyses()}).\cr\cr For \code{[}: Can also be logical to perform logical selection
(similar to regular vectors). If missing all analyses are selected.\cr\cr For \code{[[}: should be a scalar value.\cr\cr For \code{delete}: The data to remove from. \code{i} are the
analyses as numeric index, logical or character, \code{j} the features as numeric index (row) of the feature. If either is
\code{NULL} then data for all is removed. \code{j} may also be a function: it will be called for each 
analysis, with the feature table (a \code{data.table}) as first argument, the analysis name as second argument, and any other arguments passed as
\code{\dots} to \code{delete}. The return value of this function specifies the feature indices (rows) to be removed (specified as an \code{integer} or \code{logical} vector).}

\item{\dots}{For \code{delete}: passed to the function specified as \code{j}.

\setsPassedArgs1{features}}

\item{drop}{ignored.}

\item{name}{The analysis name (partially matched).}

\item{weights}{A named \code{numeric} vector that defines the weight for each score to calculate the
\verb{totalScore}. The names of the vector follow the score names. Unspecified weights are defaulted to \samp{1}.
Example: \code{weights=c(ApexBoundaryRatioScore=0.5, GaussianSimilarityScore=2)}.}

\item{flatnessFactor}{Passed to \pkg{MetaClean} as the \code{flatness.factor} argument to
\code{\link[MetaClean]{calculateJaggedness}} and \code{\link[MetaClean]{calculateModality}}.}

\item{parallel}{If set to \code{TRUE} then code is executed in parallel through the \CRANpkg{futures} package. Please
see the parallelization section in the handbook for more details.}

\item{sets}{\setsWF For \code{[} and \code{filter}: a \code{character} with name(s) of the sets to keep (or remove if
\code{negate=TRUE}).}

\item{set}{\setsWF The name of the set.}
}
\value{
\code{featureTable}: A \code{list} containing a
  \code{\link{data.table}} for each analysis with feature data

\code{analysisInfo}: A \code{data.frame} containing a column with
  analysis name (\code{analysis}), its path (\code{path}), and other columns
  such as replicate group name (\code{group}) and blank reference
  (\code{blank}).

\code{delete} returns the object for which the specified data was removed.

\code{calculatePeakQualities} returns a modified object amended with peak qualities and scores.
}
\description{
Holds information for all features present within a set of analysis.
}
\details{
This class provides a way to store intensity, retention times, \emph{m/z} and
other data for all features in a set of analyses. The class is \code{virtual}
and derived objects are created by 'feature finders' such as
\code{findFeaturesOpenMS}, \code{findFeaturesXCMS} and
\code{findFeaturesBruker}.
}
\section{Methods (by generic)}{
\itemize{
\item \code{length}: Obtain total number of features.

\item \code{show}: Shows summary information for this object.

\item \code{featureTable}: Get table with feature information

\item \code{analysisInfo}: Get analysis information

\item \code{analyses}: returns a \code{character} vector with the names of the
analyses for which data is present in this object.

\item \code{replicateGroups}: returns a \code{character} vector with the names of the
replicate groups for which data is present in this object.

\item \code{as.data.table}: Returns all feature data in a table.

\item \code{filter}: Performs common rule based filtering of features. Note
that this (and much more) functionality is also provided by the
\code{filter} method defined for \code{\link{featureGroups}}. However,
filtering a \code{features} object may be useful to avoid grouping large
amounts of features.

\item \code{[}: Subset on analyses.

\item \code{[[}: Extract a feature table for an analysis.

\item \code{$}: Extract a feature table for an analysis.

\item \code{delete}: Completely deletes specified features.

\item \code{calculatePeakQualities}: Calculates peak qualities for each feature. This uses
\href{https://github.com/KelseyChetnik/MetaClean/}{MetaClean} \R package to calculate the following metrics:
\verb{Apex-Boundary Ratio}, \verb{FWHM2Base}, \verb{Jaggedness}, \verb{Modality}, \verb{Symmetry}, \verb{Gaussian
Similarity}, \verb{Sharpness}, \verb{Triangle Peak Area Similarity Ratio} and \verb{Zig-Zag index}. Please see the
\pkg{MetaClean} publication (referenced below) for more details. For each metric, an additional score is calculated
by normalizing all feature values (unless the quality metric definition has a fixed range) and scale from \samp{0}
(worst) to \samp{1} (best). Then, a \verb{totalScore} for each feature is calculated by the (weighted) sum of all
score values.
}}

\section{Slots}{

\describe{
\item{\code{features}}{List of features per analysis file. Use the
\code{featureTable} method for access.}

\item{\code{analysisInfo}}{Analysis group information. Use the \code{analysisInfo} method
for access.}
}}

\note{
For \code{calculatePeakQualities}: sometimes \pkg{MetaClean} may return \code{NA} for the \verb{Gaussian
  Similarity} metric, in which case it will be set to \samp{0}.
}
\section{S4 class hierarchy}{
 \itemize{   \item{\code{\link{workflowStep}}}   \itemize{     \item{\strong{\code{\link{features}}}}     \itemize{       \item{\code{\link{featuresSet}}}       \item{\code{\link{featuresUnset}}}       \item{\code{\link{featuresFromFeatGroups}}}       \item{\code{\link{featuresConsensus}}}       \item{\code{\link{featuresBruker}}}       \item{\code{\link{featuresEnviPick}}}       \item{\code{\link{featuresKPIC2}}}       \item{\code{\link{featuresOpenMS}}}       \item{\code{\link{featuresSAFD}}}       \item{\code{\link{featuresSIRIUS}}}       \item{\code{\link{featuresBrukerTASQ}}}       \item{\code{\link{featuresXCMS}}}       \item{\code{\link{featuresXCMS3}}}     }   } }
}

\section{Sets workflows}{
 \setsWFClass{featuresSet}{features}

  \setsWFNewMethodsFeat{featuresUnset}{The adduct annotations for the selected set (\emph{e.g.} as passed to
  \code{makeSet}) are used to convert all feature masses to ionic \emph{m/z} values. }

  \setsWFChangedMethods{

  \item \code{filter} and the subset operator (\code{[}) have specific arguments to choose/filter by (feature
  presence in) sets. See the \code{sets} argument description.

  }
}

\references{
\insertRef{Chetnik2020}{patRoon}
}
\seealso{
\code{\link{findFeatures}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/components-specclust.R
\docType{class}
\name{componentsSpecClust-class}
\alias{componentsSpecClust-class}
\alias{componentsSpecClust}
\title{Components based on MS/MS similarity.}
\description{
This class is derived from \code{\link{componentsClust}} and is used to store components from feature groups that
were clustered based on their MS/MS similarities.
}
\details{
Objects from this class are generated by \code{\link{generateComponentsSpecClust}}
}
\note{
When the object is altered (\emph{e.g.} by filtering or subsetting it), methods that need the original
  clustered data such as plotting methods do not work anymore and stop with an error.
}
\section{S4 class hierarchy}{
 \itemize{   \item{\code{\link{componentsClust}}}   \itemize{     \item{\strong{\code{\link{componentsSpecClust}}}}   } }
}

\seealso{
\code{\link{componentsClust}} for other relevant methods and \code{\link{generateComponents}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-exported.R
\name{printPackageOpts}
\alias{printPackageOpts}
\title{Prints all the package options of \code{patRoon} and their currently set values.}
\usage{
printPackageOpts()
}
\description{
Prints all the package options of \code{patRoon} and their currently set values.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compounds-metfrag.R
\docType{class}
\name{compoundsMF-class}
\alias{compoundsMF-class}
\alias{compoundsMF}
\alias{settings,compoundsMF-method}
\alias{settings}
\title{Compounds list class for MetFrag results.}
\usage{
\S4method{settings}{compoundsMF}(compoundsMF)
}
\arguments{
\item{compoundsMF}{A \code{compoundsMF} object.}
}
\description{
This class is derived from \code{\link{compounds}} and contains additional
specific MetFrag data.
}
\details{
Objects from this class are generated by
\code{\link{generateCompoundsMetFrag}}
}
\section{Methods (by generic)}{
\itemize{
\item \code{settings}: Accessor method for the \code{settings} slot.
}}

\section{Slots}{

\describe{
\item{\code{settings}}{A list with all general configuration settings passed to
MetFrag. Feature specific items (\emph{e.g.} spectra and precursor masses)
are not contained in this list.}
}}

\section{S4 class hierarchy}{
 \itemize{   \item{\code{\link{compounds}}}   \itemize{     \item{\strong{\code{\link{compoundsMF}}}}   } }
}

\references{
\insertRef{Ruttkies2016}{patRoon}
}
\seealso{
\code{\link{compounds}} and \code{\link{generateCompoundsMetFrag}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/components-clust.R
\docType{class}
\name{componentsClust-class}
\alias{componentsClust-class}
\alias{componentsClust}
\alias{delete,componentsClust-method}
\alias{clusters,componentsClust-method}
\alias{cutClusters,componentsClust-method}
\alias{clusterProperties,componentsClust-method}
\alias{treeCut,componentsClust-method}
\alias{treeCutDynamic,componentsClust-method}
\alias{plot,componentsClust,missing-method}
\alias{plotSilhouettes,componentsClust-method}
\title{Base class for components that are based on hierarchical clustered data.}
\usage{
\S4method{delete}{componentsClust}(obj, ...)

\S4method{clusters}{componentsClust}(obj)

\S4method{cutClusters}{componentsClust}(obj)

\S4method{clusterProperties}{componentsClust}(obj)

\S4method{treeCut}{componentsClust}(obj, k = NULL, h = NULL)

\S4method{treeCutDynamic}{componentsClust}(obj, maxTreeHeight, deepSplit, minModuleSize)

\S4method{plot}{componentsClust,missing}(
  x,
  pal = "Paired",
  numericLabels = TRUE,
  colourBranches = length(x) < 50,
  showLegend = length(x) < 20,
  ...
)

\S4method{plotSilhouettes}{componentsClust}(obj, kSeq, pch = 16, type = "b", ...)
}
\arguments{
\item{\dots}{Further options passed to \code{\link{plot.dendrogram}} (\code{plot}) or \code{\link[graphics]{plot}}
(\code{plotSilhouettes}).}

\item{k, h}{Desired number of clusters or tree height to be used for cutting the dendrogram, respectively. One or the
other must be specified. Analogous to \code{\link{cutree}}.}

\item{maxTreeHeight, deepSplit, minModuleSize}{Arguments used by
\code{\link{cutreeDynamicTree}}.}

\item{x, obj}{A \code{componentsClust} (derived) object.}

\item{pal}{Colour palette to be used from \pkg{\link{RColorBrewer}}.}

\item{numericLabels}{Set to \code{TRUE} to label with numeric indices instead of (long) feature group names.}

\item{colourBranches}{Whether branches from cut clusters (and their labels)
should be coloured. Might be slow with large numbers of clusters, hence,
the default is only \code{TRUE} when this is not the case.}

\item{showLegend}{If \code{TRUE} and \code{colourBranches} is also
\code{TRUE} then a legend will be shown which outlines cluster numbers and
their colours. By default \code{TRUE} for small amount of clusters to avoid
overflowing the plot.}

\item{kSeq}{An integer vector containing the sequence that should be used for
average silhouette width calculation.}

\item{pch, type}{Passed to \code{\link[graphics]{plot}}.}
}
\description{
This base class is derived from \code{\link{components}} and is used to store components resulting from hierarchical
clustering information, for instance, generated by \code{\link{generateComponentsIntClust}} and
\code{\link{generateComponentsSpecClust}}.
}
\section{Methods (by generic)}{
\itemize{
\item \code{clusters}: Accessor method to the \code{clust} slot, which was generated by \code{\link{hclust}}.

\item \code{cutClusters}: Accessor method to the \code{cutClusters} slot. Returns a vector with cluster membership
for each candidate (format as \code{\link{cutree}}).

\item \code{clusterProperties}: Returns a list with properties on how the
clustering was performed.

\item \code{treeCut}: Manually (re-)cut the dendrogram.

\item \code{treeCutDynamic}: Automatically (re-)cut the dendrogram using the \code{\link{cutreeDynamicTree}} function
from \pkg{\link{dynamicTreeCut}}.

\item \code{plot}: generates a dendrogram from a given cluster object and optionally highlights resulting
branches when the cluster is cut.

\item \code{plotSilhouettes}: Plots the average silhouette width when the
clusters are cut by a sequence of k numbers. The k value with the highest
value (marked in the plot) may be considered as the optimal number of
clusters.
}}

\section{Slots}{

\describe{
\item{\code{distm}}{Distance matrix that was used for clustering (obtained with \code{\link{daisy}}).}

\item{\code{clust}}{Object returned by \code{\link{hclust}}.}

\item{\code{cutClusters}}{A \code{list} with assigned clusters (same format as what \code{\link{cutree}} returns).}

\item{\code{gInfo}}{The \code{\link{groupInfo}} of the feature groups object that was used.}

\item{\code{properties}}{A list containing general properties and parameters used for clustering.}

\item{\code{altered}}{Set to \code{TRUE} if the object was altered (\emph{e.g.} filtered) after its creation.}
}}

\note{
The intensity values for components (used by \code{plotSpectrum}) are set
  to a dummy value (1) as no single intensity value exists for this kind of
  components.

When the object is altered (\emph{e.g.} by filtering or subsetting it), methods that need the original
  clustered data such as plotting methods do not work anymore and stop with an error.
}
\section{S4 class hierarchy}{
 \itemize{   \item{\code{\link{components}}}   \itemize{     \item{\strong{\code{\link{componentsClust}}}}     \itemize{       \item{\code{\link{componentsIntClust}}}       \item{\code{\link{componentsSpecClust}}}     }   } }
}

\references{
\insertRef{Scholle2018}{patRoon}
}
\seealso{
\code{\link{components}} and \code{\link{generateComponents}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-mspeaklists.R
\name{specSimParams}
\alias{specSimParams}
\alias{getDefSpecSimParams}
\title{MS spectral similarity calculation parameters}
\usage{
getDefSpecSimParams(...)
}
\arguments{
\item{\dots}{optional named arguments that override defaults.}
}
\description{
Parameters relevant for calculation of similarities between mass spectra.
}
\details{
For the calculation of spectral similarities the following parameters exist:

\itemize{

\item \code{method} The similarity method: either \code{"cosine"} or \code{"jaccard"}.

\item \code{removePrecursor} If \code{TRUE} then precursor peaks (\emph{i.e.} the mass peak corresponding to the
feature) are removed prior to similarity calculation.

\item \code{mzWeight},\code{intWeight} Mass and intensity weights used for cosine calculation.

\item \code{absMzDev} Maximum absolute \emph{m/z} deviation between mass peaks, used for binning spectra.

\item \code{relMinIntensity} The minimum relative intensity for mass peaks (\samp{0-1}). Peaks with lower intensities
are not considered for similarity calculation. The relative intensities are called after the precursor peak is
removed when \code{removePrecursor=TRUE}.

\item \code{minPeaks} Only consider spectra that have at least this amount of peaks (\emph{after} the spectrum is
filtered).

\item \code{shift} If and how shifting is applied prior to similarity calculation. Valid options are: \code{"none"}
(no shifting), \code{"precursor"} (all mass peaks of the second spectrum are shifted by the mass difference between
the precursors of both spectra) or \code{"both"} (the spectra are first binned without shifting, and peaks still
unaligned are then shifted as is done when \code{shift="precursor"}).

\item \code{setCombinedMethod} \setsWF Determines how spectral similarities from different sets are combined.
Possible values are \code{"mean"}, \code{"min"} or \code{"max"}, which calculates the combined value as the mean,
minimum or maximum value, respectively. \code{NA} values (\emph{e.g.} if a set does not have peak list data to
combine) are removed in advance.

}

These parameters are typically passed as a named \code{list} as the \code{specSimParams} argument to functions that
do spectral similarity calculations. The \code{getDefSpecSimParams} function can be used to generate such parameter
list with defaults.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups-comparison.R
\docType{class}
\name{featureGroupsComparison-class}
\alias{featureGroupsComparison-class}
\alias{featureGroupsComparison}
\alias{names,featureGroupsComparison-method}
\alias{length,featureGroupsComparison-method}
\alias{[,featureGroupsComparison,ANY,missing,missing-method}
\alias{[[,featureGroupsComparison,ANY,missing-method}
\alias{$,featureGroupsComparison-method}
\alias{featureGroupsComparisonSet-class}
\alias{featureGroupsComparisonSet}
\title{Feature groups comparison class}
\usage{
\S4method{names}{featureGroupsComparison}(x)

\S4method{length}{featureGroupsComparison}(x)

\S4method{[}{featureGroupsComparison,ANY,missing,missing}(x, i, j, ..., drop = TRUE)

\S4method{[[}{featureGroupsComparison,ANY,missing}(x, i, j)

\S4method{$}{featureGroupsComparison}(x, name)
}
\arguments{
\item{x}{A \code{featureGroupsComparison} object.}

\item{i}{For \code{[}/\code{[[}: A numeric or character value which is used to select labels by
their index or name, respectively (for the order/names see \code{names()}).\cr\cr For \code{[}: Can also be logical to perform logical selection
(similar to regular vectors). If missing all labels are selected.\cr\cr For \code{[[}: should be a scalar value.}

\item{\dots}{Ignored.}

\item{drop, j}{ignored.}

\item{name}{The label name (partially matched).}
}
\description{
This class is used for comparing different \code{\link{featureGroups}}
objects.
}
\details{
Objects from this class are returned by \code{\link{comparison}}.
}
\section{Methods (by generic)}{
\itemize{
\item \code{names}: Obtain the labels that were given to each compared feature group.

\item \code{length}: Number of feature groups objects that were compared.

\item \code{[}: Subset on labels that were assigned to compared feature groups.

\item \code{[[}: Extract a \code{\link{featureGroups}} object by its label.

\item \code{$}: Extract a compound table for a feature group.
}}

\section{Slots}{

\describe{
\item{\code{fGroupsList}}{A \code{list} of \code{\link{featureGroups}} object that
were compared}

\item{\code{comparedFGroups}}{A \emph{pseudo} \code{featureGroups} object containing
grouped feature groups.}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TP-library.R
\name{generateTPsLibrary}
\alias{generateTPsLibrary}
\title{Obtain transformation products (TPs) from a library}
\usage{
generateTPsLibrary(
  parents = NULL,
  TPLibrary = NULL,
  skipInvalid = TRUE,
  matchParentsBy = "InChIKey"
)
}
\arguments{
\item{parents}{The parents for which transformation products should be obtained. This can be (1) a suspect list (see
\link[=suspect-screening]{suspect screening} for more information), (2) the resulting output
\code{\link{screenSuspects}} or (3) a \code{\link{compounds}} annotation object. In the former two cases, the
suspect (hits) are used as parents, whereas in the latter case all candidates are used as parents.
If \code{NULL} then TPs for all parents in the library are obtained.}

\item{TPLibrary}{If \code{NULL}, a default \href{https://doi.org/10.5281/zenodo.5644560}{PubChem} based library is
used. Otherwise, \code{TPLibrary} should be a \code{data.frame}. See the details below.}

\item{skipInvalid}{If set to \code{TRUE} then the parents will be skipped (with a warning) for which insufficient
information (\emph{e.g.} SMILES) is available.}

\item{matchParentsBy}{A \code{character} that specifies how the input parents are matched with the data from the TP
library. Valid options are: \code{"InChIKey"}, \code{"InChIKey1"}, \code{"InChI"}, \code{"SMILES"}.}
}
\value{
The TPs are stored in an object from the \code{\link{transformationProductsLibrary}} class.
}
\description{
Automatically obtains transformation products from a library.
}
\details{
This function uses a library to obtain transformation products. This function is called when calling \code{generateTPs} with
  \code{algorithm="library"}.

By default, a library is used that is based on data from
  \href{https://doi.org/10.5281/zenodo.5644560}{PubChem}. However, it also possible to use your own library.

An important advantage of this algorithm is that it provides structural information for generated TPs.
  However, this also means that if the input is from a parent suspect list or screening then either \acronym{SMILES}
  or \acronym{InChI} information must be available for the parents.
}
\note{
When the \code{parents} argument is a \code{\link{compounds}} object, the candidate library \code{identifier}
  is used in case the candidate has no defined \code{compoundName}.
}
\section{Custom TP libraries}{
 The \code{TPLibrary} argument is used to specify a custom TP library. This should be a
  \code{data.frame} where each row specifies a TP for a parent, with the following columns: \itemize{

  \item \code{parent_name} and \code{TP_name}: The name of the parent/TP.

  \item \code{parent_SMILES} and \code{TP_SMILES} The \acronym{SMILES} of the parent structure.

  \item \code{parent_LogP} and \code{TP_LogP} The \code{log P} values for the parent/TP. (\strong{optional})

  \item \code{LogPDiff} The difference between parent and TP \code{Log P} values. Ignored if \emph{both}
  \code{parent_LogP} and \code{TP_LogP} are specified. (\strong{optional})

  }

  Other columns are allowed, and will be included in the final object. Multiple TPs for a single parent are specified
  by repeating the value within \code{parent_} columns.
}

\seealso{
\code{\link{generateTPs}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups.R, R/feature_groups-set.R
\name{groupFeatures}
\alias{groupFeatures}
\alias{groupFeatures,features-method}
\alias{groupFeatures,data.frame-method}
\alias{groupFeatures,featuresSet-method}
\title{Grouping of features}
\usage{
\S4method{groupFeatures}{features}(obj, algorithm, ..., verbose = TRUE)

\S4method{groupFeatures}{data.frame}(obj, algorithm, ..., verbose = TRUE)

\S4method{groupFeatures}{featuresSet}(obj, algorithm, ..., verbose = TRUE)
}
\arguments{
\item{obj}{Either a \code{\link{features}} object to be grouped, or a \code{data.frame} with
\link[=analysis-information]{analysis info} to be passed to \code{groupFeaturesSIRIUS}}

\item{algorithm}{A \code{character} that specifies the algorithm to be used: either \code{"openms"}, \code{"xcms"},
\code{"xcms3"} or \code{"kpic2"} (\code{features method}), or \code{"sirius"} (\code{data.frame} method).}

\item{\dots}{Further parameters passed to the selected grouping algorithm.}

\item{verbose}{if \code{FALSE} then no text output will be shown.}
}
\value{
An object of a class which is derived from \code{\link{featureGroups}}.

The \code{featuresSet} method (for \link[=sets-workflow]{sets workflows}) returns a
  \code{\link{featureGroupsSet}} object.
}
\description{
Group equal features across analyses.
}
\details{
After \link[=findFeatures]{features have been found}, the next step is to align and group them across analyses. This
process is necessary to allow comparison of features between multiple analyses, which otherwise would be difficult
due to small deviations in retention and mass data. Thus, algorithms of 'feature groupers' are used to collect
features with similar retention and mass data. In addition, advanced retention time alignment algorithms exist to
enhance grouping of features even with relative large retention time deviations (\emph{e.g.} possibly observed from
analyses collected over a long period). Like \link{findFeatures}, various algorithms are supported which may have
many parameters that can be fine-tuned. This fine-tuning is likely to be necessary, since optimal settings often
depend on applied methodology and instrumentation.

\code{groupFeatures} is a generic function that will groupFeatures by one of the supported algorithms. The actual
  functionality is provided by algorithm specific functions such as \code{groupFeaturesOpenMS} and \code{groupFeaturesXCMS3}. While these
  functions may be called directly, \code{groupFeatures} provides a generic interface and is therefore usually preferred.

The \code{data.frame} method for \code{groupFeatures} is a special case that currently only supports the
  \code{"sirius"} algorithm.
}
\seealso{
The \code{\link{featureGroups}} output class and its methods and the algorithm specific functions:
  \code{\link{groupFeaturesOpenMS}}, \code{\link{groupFeaturesXCMS}}, \code{\link{groupFeaturesXCMS3}}, \code{\link{groupFeaturesKPIC2}}, \code{\link{groupFeaturesSIRIUS}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-adduct.R
\name{adduct-utils}
\alias{adduct-utils}
\alias{GenFormAdducts}
\alias{MetFragAdducts}
\alias{as.adduct}
\alias{calculateIonFormula}
\alias{calculateNeutralFormula}
\title{Adduct utilities}
\usage{
GenFormAdducts()

MetFragAdducts()

as.adduct(x, format = "generic", isPositive = NULL, charge = NULL)

calculateIonFormula(formula, adduct)

calculateNeutralFormula(formula, adduct)
}
\arguments{
\item{x}{The object that should be converted. Should be a \code{character}
string, a \code{numeric} \command{MetFrag} adduct identifier
(\code{adduct_mode} column obtained with \code{MetFragAdducts}) or an
\code{\link{adduct}} object (in which case no conversion occurs).}

\item{format}{A \code{character} that specifies the source format.

  \code{"generic"} is an internally used generic format that supports full
  textual conversion. Examples: \code{"[M+H]+"}, \code{"[2M+H]+"},
  \code{"[M+3H]3+"}.

  \code{"sirius"} Is the format used by \command{SIRIUS}. It is similar to
  \code{generic} but does not allow multiple charges/molecules. See the
  SIRIUS manual for more details.

  \code{"genform"} and \code{"metfrag"} support fixed types of adducts
  which can be obtained with the \code{GenFormAdducts} and
  \code{MetFragAdducts} functions, respectively.
  
  \code{"openms"} is the format used by the \command{MetaboliteAdductDecharger} tool.
  
  \code{"cliquems"} is the format used by \pkg{cliqueMS}.}

\item{isPositive}{A logical that specifies whether the adduct should be
positive. Should only be set when \code{format="metfrag"} and \code{x} is a
\code{numeric} identifier.}

\item{charge}{The final charge. Only needs to be set when \code{format="openms"}.}

\item{formula}{A \code{character} vector with formulae to convert.}

\item{adduct}{An \code{\link{adduct}} object (or something that can be converted to it with \code{\link{as.adduct}}).
Examples: \code{"[M-H]-"}, \code{"[M+Na]+"}.}
}
\description{
Several utility functions to work with adducts.
}
\details{
\code{GenFormAdducts} returns a table with information on adducts
  supported by \command{GenForm}.

\code{MetFragAdducts} returns a table with information on adducts
  supported by \command{MetFrag}.

\code{as.adduct} Converts an object in to an \code{\link{adduct}}
  object.

\code{calculateIonFormula} Converts one or more neutral formulae to
  adduct ions.

\code{calculateNeutralFormula} Converts one or more adduct ions to
  neutral formulae.
}
\examples{
as.adduct("[M+H]+")
as.adduct("[M+H2]2+")
as.adduct("[2M+H]+")
as.adduct("[M-H]-")
as.adduct("+H", format = "genform")
as.adduct(1, isPositive = TRUE, format = "metfrag") # MetFrag adduct ID 1 --> returns [M+H]+

calculateIonFormula("C2H4O", "[M+H]+") # C2H5O
calculateNeutralFormula("C2H5O", "[M+H]+") # C2H4O

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/components-nontarget.R,
%   R/components-nontarget-set.R
\docType{class}
\name{componentsNT-class}
\alias{componentsNT-class}
\alias{componentsNT}
\alias{plotGraph,componentsNT-method}
\alias{componentsNTSet-class}
\alias{componentsNTSet}
\alias{plotGraph,componentsNTSet-method}
\alias{componentsNTUnset-class}
\alias{componentsNTUnset}
\alias{unset,componentsNTSet-method}
\title{Components class for homologous series.}
\usage{
\S4method{plotGraph}{componentsNT}(obj, onlyLinked)

\S4method{plotGraph}{componentsNTSet}(obj, onlyLinked, set)

\S4method{unset}{componentsNTSet}(obj, set)
}
\arguments{
\item{obj}{The \code{componentsNT} object to plot.}

\item{onlyLinked}{If \code{TRUE} then only components with links are plotted.}

\item{set}{\setsWF The name of the set.}
}
\value{
\code{plotGraph} returns the result of \code{\link{visNetwork}}.
}
\description{
This class is derived from \code{\link{components}} and is used to store
results from unsupervised homolog detection with the \pkg{\link{nontarget}}
package.
}
\details{
Objects from this class are generated by
\code{\link{generateComponentsNontarget}}
}
\section{Methods (by generic)}{
\itemize{
\item \code{plotGraph}: Plots an interactive network graph for linked
homologous series (\emph{i.e.} series with (partial) overlap which could
not be merged). The resulting graph can be browsed interactively and allows
quick inspection of series which may be related. The graph is constructed
with the \pkg{\link{igraph}} package and rendered with
\pkg{\link{visNetwork}}.
}}

\section{Slots}{

\describe{
\item{\code{homol}}{A \code{list} with \code{homol} objects for each replicate group
as returned by \code{\link{homol.search}}}
}}

\section{Sets workflows}{
 \setsWFClass{componentsNTSet}{componentsNT}

  \setsWFNewMethodsSO{componentsNTUnset}{Only the components in the specified set are kept. Furthermore, the
  component names are restored to non-set specific names (see \code{\link{generateComponents}} for more details).}

  \setsWFChangedMethods{

  \item \code{plotGraph} Currently can only create graph networks from one set (specified by the \code{set}
  argument).

  }

  Note that the \code{componentsNTSet} class does not have a \code{homol} slot. Instead, the \code{\link{setObjects}}
  method can be used to access this data for a specific set.
}

\references{
\addCitations{nontarget}{1} \cr\cr \addCitations{enviPat}{1}

\addCitations{igraph}{1} \cr \cr \addCitations{visNetwork}{1}
}
\seealso{
\code{\link{components}} and \code{\link{generateComponents}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features-envipick.R
\name{importFeaturesEnviMass}
\alias{importFeaturesEnviMass}
\title{Imports features from enviMass}
\usage{
importFeaturesEnviMass(analysisInfo, enviProjPath)
}
\arguments{
\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}

\item{enviProjPath}{The path of the enviMass project.}
}
\value{
An object of a class which is derived from \code{\link{features}}.
}
\description{
Imports features from a project generated by the \pkg{enviMass} package.
}
\details{
This function imports data from enviMass. This function is called when calling \code{importFeatures} with
  \code{type="envimass"}.
}
\note{
This functionality has only been tested with older versions of \pkg{enviMass}.
}
\seealso{
\code{\link{importFeatures}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compounds-metfrag.R
\name{generateCompoundsMetFrag}
\alias{generateCompoundsMetFrag}
\alias{generateCompoundsMetFrag,featureGroups-method}
\alias{generateCompoundsMetFrag,featureGroupsSet-method}
\title{Compound annotation with MetFrag}
\usage{
\S4method{generateCompoundsMetFrag}{featureGroups}(
  fGroups,
  MSPeakLists,
  method = "CL",
  timeout = 300,
  timeoutRetries = 2,
  errorRetries = 2,
  topMost = 100,
  dbRelMzDev = 5,
  fragRelMzDev = 5,
  fragAbsMzDev = 0.002,
  adduct = NULL,
  database = "pubchem",
  extendedPubChem = "auto",
  chemSpiderToken = "",
  scoreTypes = compoundScorings("metfrag", database, onlyDefault = TRUE)$name,
  scoreWeights = 1,
  preProcessingFilters = c("UnconnectedCompoundFilter", "IsotopeFilter"),
  postProcessingFilters = c("InChIKeyFilter"),
  maxCandidatesToStop = 2500,
  identifiers = NULL,
  extraOpts = NULL
)

\S4method{generateCompoundsMetFrag}{featureGroupsSet}(
  fGroups,
  MSPeakLists,
  method = "CL",
  timeout = 300,
  timeoutRetries = 2,
  errorRetries = 2,
  topMost = 100,
  dbRelMzDev = 5,
  fragRelMzDev = 5,
  fragAbsMzDev = 0.002,
  adduct = NULL,
  ...,
  setThreshold = 0,
  setThresholdAnn = 0
)
}
\arguments{
\item{fGroups}{\code{\link{featureGroups}} object which should be annotated. This should be the same or a subset of
the object that was used to create the specified \code{MSPeakLists}. In the case of a subset only the remaining
feature groups in the subset are considered.}

\item{MSPeakLists}{A \code{\link{MSPeakLists}} object that was generated for the supplied \code{fGroups}.}

\item{method}{Which method should be used for MetFrag execution: \code{"CL"} for \command{MetFragCL} and \code{"R"}
for \command{MetFragR}. The former is usually much faster and recommended.}

\item{timeout}{Maximum time (in seconds) before a metFrag query for a feature group is stopped. Also see
\code{timeoutRetries} argument.}

\item{timeoutRetries}{Maximum number of retries after reaching a timeout before completely skipping the metFrag query
for a feature group. Also see \code{timeout} argument.}

\item{errorRetries}{Maximum number of retries after an error occurred. This may be useful to handle e.g. connection
errors.}

\item{topMost}{Only keep this number of candidates (per feature group) with highest score. Set to \code{NULL} to
always keep all candidates, however, please note that this may result in significant usage of CPU/RAM resources for
large numbers of candidates.}

\item{dbRelMzDev}{Relative mass deviation (in ppm) for database search. Sets the
\option{DatabaseSearchRelativeMassDeviation} option.}

\item{fragRelMzDev}{Relative mass deviation (in ppm) for fragment matching. Sets the
\option{FragmentPeakMatchRelativeMassDeviation} option.}

\item{fragAbsMzDev}{Absolute mass deviation (in Da) for fragment matching. Sets the
\option{FragmentPeakMatchAbsoluteMassDeviation} option.}

\item{adduct}{An \code{\link{adduct}} object (or something that can be converted to it with \code{\link{as.adduct}}).
  Examples: \code{"[M-H]-"}, \code{"[M+Na]+"}. If the \code{featureGroups} object has
  adduct annotations then these are used if \code{adducts=NULL}.

  \setsWF The \code{adduct} argument is not supported for sets workflows, since the
  adduct annotations will then always be used.}

\item{database}{Compound database to use. Valid values are: \code{"pubchem"}, \code{"chemspider"},
\code{"for-ident"}, \code{"comptox"}, \code{"pubchemlite"}, \code{"kegg"}, \code{"sdf"}, \code{"psv"} and
\code{"csv"}. See section below for more information. Sets the \code{MetFragDatabaseType} option.}

\item{extendedPubChem}{If \code{database="pubchem"}: whether to use the \emph{extended} database that includes
information for compound scoring (\emph{i.e.} number of patents/PubMed references). Note that downloading
candidates from this database might take extra time. Valid values are: \code{FALSE} (never use it), \code{TRUE}
(always use it) or \code{"auto"} (default, use if specified scorings demand it).}

\item{chemSpiderToken}{A character string with the \href{http://www.chemspider.com/AboutServices.aspx}{ChemSpider
security token} that should be set when the ChemSpider database is used. Sets the \option{ChemSpiderToken} option.}

\item{scoreTypes}{A character vector defining the scoring types. See the \verb{Scorings} section below for more
information. Note that both generic and \command{MetFrag} specific names are accepted (\emph{i.e.} \code{name} and
\code{metfrag} columns returned by \code{\link{compoundScorings}}). When a local database is used, the name should
match what is given there (\code{e.g} column names when \code{database=csv}). Note that MetFrag may still report
other scoring data, however, these are not used for ranking. Sets the \option{MetFragScoreTypes} option.}

\item{scoreWeights}{Numeric vector containing weights of the used scoring types. Order is the same as set in
\code{scoreTypes}. Values are recycled if necessary. Sets the \option{MetFragScoreWeights} option.}

\item{preProcessingFilters, postProcessingFilters}{A character vector defining pre/post filters applied before/after
fragmentation and scoring (\emph{e.g.} \code{"UnconnectedCompoundFilter"}, \code{"IsotopeFilter"},
\code{"ElementExclusionFilter"}). Some methods require further options to be set. For all filters and more
information refer to the \verb{Candidate Filters} section on the
\href{http://ipb-halle.github.io/MetFrag/projects/metfragr/}{MetFragR homepage}. Sets the
\option{MetFragPreProcessingCandidateFilter} and \code{MetFragPostProcessingCandidateFilter} options.}

\item{maxCandidatesToStop}{If more than this number of candidate structures are found then processing will be aborted
and no results this feature group will be reported. Low values increase the chance of missing data, whereas too
high values will use too much computer resources and signficantly slowdown the process. Sets the
\option{MaxCandidateLimitToStop} option.}

\item{identifiers}{A \code{list} containing for each feature group a character vector with database identifiers that
should be used to find candidates for a feature group (the list should be named by feature group names). If
\code{NULL} all relevant candidates will be retrieved from the specified database. An example usage scenario is to
obtain the list of candidate identifiers from a \code{\link{compounds}} object obtained with
\code{\link{generateCompoundsSIRIUS}} using the \code{\link{identifiers}} method. This way, only those candidates
will be searched by MetFrag that were generated by SIRIUS+CSI:FingerID. Sets the \option{PrecursorCompoundIDs}
option.}

\item{extraOpts}{A named \code{list} containing further settings \command{MetFrag}. See the
\href{http://ipb-halle.github.io/MetFrag/projects/metfragr/}{MetFragR} and
\href{http://ipb-halle.github.io/MetFrag/projects/metfragcl/}{MetFrag CL} homepages for all available options. Set
to \code{NULL} to ignore.}

\item{\dots}{\setsWF Further arguments passed to the non-sets workflow method.}

\item{setThreshold}{\setsWF Minimum abundance for a candidate among all sets (\samp{0-1}). For instance, a value of
\samp{1} means that the candidate needs to be present in all the set data.}

\item{setThresholdAnn}{\setsWF As \code{setThreshold}, but only taking into account the set data that contain
annotations for the feature group of the candidate.}
}
\value{
\code{generateCompoundsMetFrag} returns a \code{\link{compoundsMF}} object.
}
\description{
Uses the \pkg{metfRag} package or \command{MetFrag CL} for compound identification (see
\url{http://ipb-halle.github.io/MetFrag/}).
}
\details{
This function uses MetFrag to generate compound candidates. This function is called when calling \code{generateCompounds} with
  \code{algorithm="metfrag"}.

Several online compound databases such as \href{https://pubchem.ncbi.nlm.nih.gov/}{PubChem} and
  \href{http://www.chemspider.com/}{ChemSpider} may be chosen for retrieval of candidate structures. This method
  requires the availability of MS/MS data, and feature groups without it will be ignored. Many options exist to score
  and filter resulting data, and it is highly suggested to optimize these to improve results. The \command{MetFrag}
  options \code{PeakList}, \code{IonizedPrecursorMass} and \code{ExperimentalRetentionTimeValue} (in minutes) fields
  are automatically set from feature data.
}
\section{Scorings}{
 \command{MetFrag} supports \emph{many} different scorings to rank candidates. The
  \code{\link{compoundScorings}} function can be used to get an overview: (some columns are omitted)

 \tabular{lll}{
  \strong{name}       \tab \strong{metfrag}                                \tab \strong{database}\cr
  score                \tab Score                                            \tab                   \cr
  fragScore            \tab FragmenterScore                                  \tab                   \cr
  metFusionScore       \tab OfflineMetFusionScore                            \tab                   \cr
  individualMoNAScore  \tab OfflineIndividualMoNAScore                       \tab                   \cr
  numberPatents        \tab PubChemNumberPatents                             \tab pubchem           \cr
  numberPatents        \tab Patent_Count                                     \tab pubchemlite       \cr
  pubMedReferences     \tab PubChemNumberPubMedReferences                    \tab pubchem           \cr
  pubMedReferences     \tab ChemSpiderNumberPubMedReferences                 \tab chemspider        \cr
  pubMedReferences     \tab NUMBER_OF_PUBMED_ARTICLES                        \tab comptox           \cr
  pubMedReferences     \tab PubMed_Count                                     \tab pubchemlite       \cr
  extReferenceCount    \tab ChemSpiderNumberExternalReferences               \tab chemspider        \cr
  dataSourceCount      \tab ChemSpiderDataSourceCount                        \tab chemspider        \cr
  referenceCount       \tab ChemSpiderReferenceCount                         \tab chemspider        \cr
  RSCCount             \tab ChemSpiderRSCCount                               \tab chemspider        \cr
  smartsInclusionScore \tab SmartsSubstructureInclusionScore                 \tab                   \cr
  smartsExclusionScore \tab SmartsSubstructureExclusionScore                 \tab                   \cr
  suspectListScore     \tab SuspectListScore                                 \tab                   \cr
  retentionTimeScore   \tab RetentionTimeScore                               \tab                   \cr
  CPDATCount           \tab CPDAT_COUNT                                      \tab comptox           \cr
  TOXCASTActive        \tab TOXCAST_PERCENT_ACTIVE                           \tab comptox           \cr
  dataSources          \tab DATA_SOURCES                                     \tab comptox           \cr
  pubChemDataSources   \tab PUBCHEM_DATA_SOURCES                             \tab comptox           \cr
  EXPOCASTPredExpo     \tab EXPOCAST_MEDIAN_EXPOSURE_PREDICTION_MG/KG-BW/DAY \tab comptox           \cr
  ECOTOX               \tab ECOTOX                                           \tab comptox           \cr
  NORMANSUSDAT         \tab NORMANSUSDAT                                     \tab comptox           \cr
  MASSBANKEU           \tab MASSBANKEU                                       \tab comptox           \cr
  TOX21SL              \tab TOX21SL                                          \tab comptox           \cr
  TOXCAST              \tab TOXCAST                                          \tab comptox           \cr
  KEMIMARKET           \tab KEMIMARKET                                       \tab comptox           \cr
  MZCLOUD              \tab MZCLOUD                                          \tab comptox           \cr
  pubMedNeuro          \tab PubMedNeuro                                      \tab comptox           \cr
  CIGARETTES           \tab CIGARETTES                                       \tab comptox           \cr
  INDOORCT16           \tab INDOORCT16                                       \tab comptox           \cr
  SRM2585DUST          \tab SRM2585DUST                                      \tab comptox           \cr
  SLTCHEMDB            \tab SLTCHEMDB                                        \tab comptox           \cr
  THSMOKE              \tab THSMOKE                                          \tab comptox           \cr
  ITNANTIBIOTIC        \tab ITNANTIBIOTIC                                    \tab comptox           \cr
  STOFFIDENT           \tab STOFFIDENT                                       \tab comptox           \cr
  KEMIMARKET_EXPO      \tab KEMIMARKET_EXPO                                  \tab comptox           \cr
  KEMIMARKET_HAZ       \tab KEMIMARKET_HAZ                                   \tab comptox           \cr
  REACH2017            \tab REACH2017                                        \tab comptox           \cr
  KEMIWW_WDUIndex      \tab KEMIWW_WDUIndex                                  \tab comptox           \cr
  KEMIWW_StpSE         \tab KEMIWW_StpSE                                     \tab comptox           \cr
  KEMIWW_SEHitsOverDL  \tab KEMIWW_SEHitsOverDL                              \tab comptox           \cr
  ZINC15PHARMA         \tab ZINC15PHARMA                                     \tab comptox           \cr
  PFASMASTER           \tab PFASMASTER                                       \tab comptox           \cr
  peakFingerprintScore \tab AutomatedPeakFingerprintAnnotationScore          \tab                   \cr
  lossFingerprintScore \tab AutomatedLossFingerprintAnnotationScore          \tab                   \cr
  agroChemInfo         \tab AgroChemInfo                                     \tab pubchemlite       \cr
  bioPathway           \tab BioPathway                                       \tab pubchemlite       \cr
  drugMedicInfo        \tab DrugMedicInfo                                    \tab pubchemlite       \cr
  foodRelated          \tab FoodRelated                                      \tab pubchemlite       \cr
  pharmacoInfo         \tab PharmacoInfo                                     \tab pubchemlite       \cr
  safetyInfo           \tab SafetyInfo                                       \tab pubchemlite       \cr
  toxicityInfo         \tab ToxicityInfo                                     \tab pubchemlite       \cr
  knownUse             \tab KnownUse                                         \tab pubchemlite       \cr
  disorderDisease      \tab DisorderDisease                                  \tab pubchemlite       \cr
  identification       \tab Identification                                   \tab pubchemlite       \cr
  annoTypeCount        \tab FPSum                                            \tab pubchemlite       \cr
  annoTypeCount        \tab AnnoTypeCount                                    \tab pubchemlite       
}

 In addition, the \code{\link{compoundScorings}} function is also useful to programmatically
  generate a set of scorings to be used for ranking with \command{MetFrag}. For instance, the following can be given
  to the \code{scoreTypes} argument to use all default scorings for PubChem: \code{compoundScorings("metfrag",
  "pubchem", onlyDefault=TRUE)$name}.

  For all \command{MetFrag} scoring types refer to the \verb{Candidate Scores} section on the
  \href{http://ipb-halle.github.io/MetFrag/projects/metfragr/}{MetFragR homepage}.
}

\section{Usage of MetFrag databases}{
 When \code{database="chemspider"} setting the \code{chemSpiderToken} argument is
  mandatory.

  When a local database is set (\emph{i.e.} \code{sdf}, \code{psv}, \code{csv}, \code{comptox}, \code{pubchemlite})
  the file location of the database should be set in the \code{LocalDatabasePath} value via the \code{extraOpts}
  argument or using the \code{patRoon.path.MetFragCompTox}/\code{patRoon.path.MetFragPubChemLite} option (only when
  \code{database="comptox"} or \code{database="pubchemlite"}).

  Examples: \verb{options(patRoon.path.MetFragCompTox = "C:/CompTox_17March2019_SelectMetaData.csv")}

  \verb{extraOpts = list(LocalDatabasePath = "C:/myDB.csv")}.

  For \code{database="comptox"} the files can be obtained from
  \href{ftp://newftp.epa.gov/COMPTOX/Sustainable_Chemistry_Data/Chemistry_Dashboard/MetFrag_metadata_files}{here}.
  Furthermore, the files with additions for \href{smoking}{https://zenodo.org/record/3364464#.XnjM-XLvKUk} and
  \href{wastewater}{https://zenodo.org/record/3472781#.XnjMAHLvKUk} metadata are also supported. For
  \code{database="pubchemlite"} the \file{.csv} database can be downloaded from
  \href{https://zenodo.org/record/4432124/files/PubChemLite_01Jan2021_exposomics.csv}{here}. Note that only recent
  \command{MetFrag} versions (>= \samp{2.4.5}) support these libraries.
}

\section{Parallelization}{
 generateCompoundsMetFrag uses multiprocessing to parallelize
  computations. Please see the parallelization section in the handbook for
  more details and \link[=patRoon-package]{patRoon options} for configuration
  options.

 When local database files are used with \code{generateCompoundsMetFrag} (\emph{e.g.} when
  \code{database} is set to \code{"pubchemlite"}, \code{"csv"} etc.) and \option{patRoon.MP.method="future"}, then
  the database file must be present on all the nodes. When \code{pubchemlite} or \code{comptox} is used, the location
  for these databases can be configured on the host with the respective package options
  (\option{patRoon.path.MetFragPubChemLite} and \option{patRoon.path.MetFragCompTox}). Note that these files must
  \emph{also} be present on the local host computer, even if it is not participating in computations.
}

\references{
\insertRef{Ruttkies2016}{patRoon}
}
\seealso{
\code{\link{generateCompounds}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adduct.R
\docType{class}
\name{adduct-class}
\alias{adduct-class}
\alias{adduct}
\alias{show,adduct-method}
\alias{as.character,adduct-method}
\title{Generic adduct class}
\usage{
adduct(...)

\S4method{show}{adduct}(object)

\S4method{as.character}{adduct}(x, format = "generic", err = TRUE)
}
\arguments{
\item{x, object}{An \code{adduct} object.}

\item{format}{A \code{character} that specifies the source format.

  \code{"generic"} is an internally used generic format that supports full
  textual conversion. Examples: \code{"[M+H]+"}, \code{"[2M+H]+"},
  \code{"[M+3H]3+"}.

  \code{"sirius"} Is the format used by \command{SIRIUS}. It is similar to
  \code{generic} but does not allow multiple charges/molecules. See the
  SIRIUS manual for more details.

  \code{"genform"} and \code{"metfrag"} support fixed types of adducts
  which can be obtained with the \code{GenFormAdducts} and
  \code{MetFragAdducts} functions, respectively.
  
  \code{"openms"} is the format used by the \command{MetaboliteAdductDecharger} tool.
  
  \code{"cliquems"} is the format used by \pkg{cliqueMS}.}

\item{err}{If \code{TRUE} then an error will be thrown if conversion fails,
otherwise returns \code{NA}.}

\item{\dots}{Any of \code{add}, \code{sub}, \code{molMult} and/or \code{charge}. See \verb{Slots}.}
}
\description{
Objects from this class are used to specify adduct information in an
algorithm independent way.
}
\section{Methods (by generic)}{
\itemize{
\item \code{show}: Shows summary information for this object.

\item \code{as.character}: Converts an \code{adduct} object to a specified
\code{character} format.
}}

\section{Slots}{

\describe{
\item{\code{add,sub}}{A \code{character} with one or more formulas to add/subtract.}

\item{\code{molMult}}{How many times the original molecule is present in this
molecule (\emph{e.g.} for a dimer this would be \samp{2}). Default is \samp{1}.}

\item{\code{charge}}{The final charge of the adduct (default \samp{1}).}
}}

\examples{
adduct("H") # [M+H]+
adduct(sub = "H", charge = -1) # [M-H]-
adduct(add = "K", sub = "H2", charge = -1) # [M+K-H2]+
adduct(add = "H3", charge = 3) # [M+H3]3+
adduct(add = "H", molMult = 2) # [2M+H]+

as.character(adduct("H")) # returns "[M+H]+"

}
\seealso{
\code{\link{as.adduct}} for easy creation of \code{adduct} objects
  and \link[=adduct-utils]{adduct utilities} for other adduct functionality.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R, R/features-set.R,
%   R/feature_groups-set.R
\name{makeSet}
\alias{makeSet}
\alias{makeSet,features-method}
\alias{makeSet,featuresSet-method}
\alias{makeSet,featureGroups-method}
\alias{makeSet,featureGroupsSet-method}
\title{Initiate sets workflows}
\usage{
\S4method{makeSet}{features}(obj, ..., adducts, labels = NULL)

\S4method{makeSet}{featuresSet}(obj, ...)

\S4method{makeSet}{featureGroups}(
  obj,
  ...,
  groupAlgo,
  groupArgs = NULL,
  verbose = TRUE,
  adducts = NULL,
  labels = NULL
)

\S4method{makeSet}{featureGroupsSet}(obj, ...)
}
\arguments{
\item{obj, \dots}{\code{\link{features}} or \code{\link{featureGroups}} objects that should be used for the
\link[=sets-workflow]{sets workflow}.}

\item{adducts}{The adduct assignments to each set. Should either be a \code{list} with \code{\link{adduct}} objects
  or a \code{character} vector (\emph{e.g.} \code{"[M+H]+"}). The order should follow that of the objects given to
  the \code{obj} and \code{\dots} arguments.

  For the \code{\link{featureGroups}} method: if \code{NULL} then adduct annotations are used.}

\item{labels}{The labels, or \emph{set names}, for each set to be created. The order should follow that of the
objects given to the \code{obj} and \code{\dots} arguments. If \code{NULL}, then labels are automatically generated
from the polarity of the specified \code{adducts} argument (\emph{e.g.} \code{"positive"}, \code{"negative"}).}

\item{groupAlgo}{groupAlgo The name of the feature grouping algorithm. See the \code{algorithm} argument of
\code{\link{groupFeatures}} for details.}

\item{groupArgs}{A \code{list} with arguments directly passed to \code{groupFeatures} (can be named). Example:
\code{groupArgs=list(maxAlignMZ=0.002)}.}

\item{verbose}{If set to \code{FALSE} then no text output is shown.}
}
\value{
Either a \code{\link{featuresSet}} object (\code{features} method) or \code{\link{featureGroupsSet}} object
  (\code{featureGroups} method).
}
\description{
Initiate sets workflows from specified feature data.
}
\details{
The \code{makeSet} method function is used to initiate a \link[=sets-workflow]{sets workflow}. The features from
input objects are combined and then \emph{neutralized} by replacing their \emph{m/z} values by neutral monoisotopic
masses. After neutralization features measured with \emph{e.g.} different ionization polarities can be grouped since
their neutral mass will be the same.

The \link[=analysis-information]{analysis information} for this object is updated with all analyses, and a \code{set}
column is added to designate the set of each analysis. Note that currently, all analyses names \strong{must} be
unique across different sets.

\code{makeSet} supports two types of input: \enumerate{

\item \code{\link{features}} objects: \code{makeSet} combines the input objects into a \code{featuresSet} object,
which is then grouped in the 'usual way' with \code{\link{groupFeatures}}.

\item \code{\link{featureGroups}} objects: In this case the features from the input objects are first neutralized and
all features are then (re-)grouped with \code{groupFeatures}.

}

The advantage of the \code{featureGroups} method is that it preserves any adduct annotations already present
(\emph{e.g.} as set by \code{selectIons} or \code{adducts<-}). Furthermore, this approach allows more advanced
workflows where the input \code{featureGroups} are first pre-treated with \emph{e.g.} filter before the sets object
is made. On the other hand, the \code{features} method is easier, as it doesn't require intermediate feature grouping
steps and is often sufficient since adduct annotations can be made afterwards with \code{selectIons}/\code{adducts<-}
and most \code{filter} operations do not need to be done per individual set.

The adduct information used for feature neutralization is specified through the \code{adducts} argument.
Alternatively, when the \code{featureGroups} method of \code{makeSet} is used, then the adduct annotations already
present in the input objects can also by used by setting \code{adducts=NULL}. The adduct information is also used to
add adduct annotations to the output of \code{makeSet}.
}
\note{
Initiating a sets workflow recursively, \emph{i.e.} with \code{featuresSet} or \code{featureGroupsSet} objects
  as input, is currently not supported.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/components.R, R/components-camera.R,
%   R/components-features.R, R/components-cliquems.R, R/components-set.R,
%   R/components-openms.R, R/components-ramclustr.R
\docType{class}
\name{components-class}
\alias{components-class}
\alias{components}
\alias{componentTable,components-method}
\alias{componentTable}
\alias{componentInfo,components-method}
\alias{componentInfo}
\alias{groupNames,components-method}
\alias{length,components-method}
\alias{names,components-method}
\alias{show,components-method}
\alias{[,components,ANY,ANY,missing-method}
\alias{[[,components,ANY,ANY-method}
\alias{$,components-method}
\alias{delete,components-method}
\alias{as.data.table,components-method}
\alias{filter,components-method}
\alias{findFGroup,components-method}
\alias{findFGroup}
\alias{plotSpectrum,components-method}
\alias{plotChroms,components-method}
\alias{consensus,components-method}
\alias{componentsCamera-class}
\alias{componentsCamera}
\alias{componentsFeatures-class}
\alias{componentsFeatures}
\alias{show,componentsFeatures-method}
\alias{componentsCliqueMS-class}
\alias{componentsCliqueMS}
\alias{componentsSet-class}
\alias{componentsSet}
\alias{show,componentsSet-method}
\alias{[,componentsSet,ANY,ANY,missing-method}
\alias{filter,componentsSet-method}
\alias{delete,componentsSet-method}
\alias{consensus,componentsSet-method}
\alias{componentsUnset-class}
\alias{componentsUnset}
\alias{unset,componentsSet-method}
\alias{componentsOpenMS-class}
\alias{componentsOpenMS}
\alias{componentsRC-class}
\alias{componentsRC}
\title{Component class}
\usage{
\S4method{componentTable}{components}(obj)

\S4method{componentInfo}{components}(obj)

\S4method{groupNames}{components}(obj)

\S4method{length}{components}(x)

\S4method{names}{components}(x)

\S4method{show}{components}(object)

\S4method{[}{components,ANY,ANY,missing}(x, i, j, ..., drop = TRUE)

\S4method{[[}{components,ANY,ANY}(x, i, j)

\S4method{$}{components}(x, name)

\S4method{delete}{components}(obj, i = NULL, j = NULL, ...)

\S4method{as.data.table}{components}(x)

\S4method{filter}{components}(
  obj,
  size = NULL,
  adducts = NULL,
  isotopes = NULL,
  rtIncrement = NULL,
  mzIncrement = NULL,
  checkComponentsSession = NULL,
  negate = FALSE,
  verbose = TRUE
)

\S4method{findFGroup}{components}(obj, fGroup)

\S4method{plotSpectrum}{components}(obj, index, markFGroup = NULL, xlim = NULL, ylim = NULL, ...)

\S4method{plotChroms}{components}(obj, index, fGroups, rtWindow = 5, ...)

\S4method{consensus}{components}(obj, ...)

\S4method{show}{componentsFeatures}(object)

\S4method{show}{componentsSet}(object)

\S4method{[}{componentsSet,ANY,ANY,missing}(x, i, j, ..., sets = NULL, drop = TRUE)

\S4method{filter}{componentsSet}(obj, ..., negate = FALSE, sets = NULL)

\S4method{delete}{componentsSet}(obj, i = NULL, j = NULL, ...)

\S4method{consensus}{componentsSet}(obj, ...)

\S4method{unset}{componentsSet}(obj, set)
}
\arguments{
\item{obj, object, x}{The \code{component} object.}

\item{i, j}{For \code{[}/\code{[[}: A numeric or character value which is used to select components/feature groups by
their index or name, respectively (for the order/names see \code{names()/groupNames()}).\cr\cr For \code{[}: Can also be logical to perform logical selection
(similar to regular vectors). If missing all components/feature groups are selected.\cr\cr For \code{[[}: should be a scalar value. \code{j} is optional.\cr\cr For \code{delete}: The data to remove from. \code{i} are the
components as numeric index, logical or character, \code{j} the feature groups as numeric index/logical (relative to component) or character. If either is
\code{NULL} then data for all is removed. \code{j} may also be a function: it will be called for each 
component, with the component (a \code{data.table}) as first argument, the component name as second argument, and any other arguments passed as
\code{\dots} to \code{delete}. The return value of this function specifies the feature groups to be removed (same format as \code{j}).}

\item{\dots}{For \code{delete}: passed to the function specified as \code{j}.

  For \code{plotChroms}: Further (optional) arguments passed to the
  \code{plotChroms} method for the \code{\link{featureGroups}} class. Note that
  the \code{colourBy}, \code{showPeakArea}, \code{showFGroupRect} and
  \code{topMost} arguments cannot be set as these are set by this method.

  For \code{plotSpectrum}: Further arguments passed to
  \code{\link[graphics]{plot}}.

  For \code{consensus}: \code{components} objects that should be used to
  generate the consensus.
  
  \setsPassedArgs1{components}}

\item{drop}{ignored.}

\item{name}{The component name (partially matched).}

\item{size}{Should be a two sized vector with the minimum/maximum size of a component. Set to \code{NULL} to ignore.}

\item{adducts}{Remove any feature groups within components that do not match given adduct rules. If \code{adducts} is
a logical then only results are kept when an adduct is assigned (\code{adducts=TRUE}) or not assigned
(\code{adducts=FALSE}). Otherwise, if \code{adducts} contains one or more \code{\link{adduct}} objects (or
something that can be converted to it with \code{\link{as.adduct}}) then only results are kept that match the given
adducts. Set to \code{NULL} to ignore this filter.}

\item{isotopes}{Only keep results that match a given isotope rule. If \code{isotopes} is a logical then only results
are kept with (\code{isotopes=TRUE}) or without (\code{isotopes=FALSE}) isotope assignment. Otherwise
\code{isotopes} should be a numeric vector with isotope identifiers to keep (\emph{e.g.} \samp{0} for monoisotopic
results, \samp{1} for \samp{M+1} results etc.). Set to \code{NULL} to ignore this filter.}

\item{rtIncrement, mzIncrement}{Should be a two sized vector with the minimum/maximum retention or mz increment of a
homologous series. Set to \code{NULL} to ignore.}

\item{checkComponentsSession}{If set then components and/or feature groups are removed that were selected for removal
(see \link{check-GUI} and the \code{\link{checkComponents}} function). The value of \code{checkComponentsSession}
should either by a path to the session file or \code{TRUE}, in which case the default session file name is used. If
\code{negate=TRUE} then all non-selected data is removed instead.}

\item{negate}{If \code{TRUE} then filters are applied in opposite manner.}

\item{verbose}{If set to \code{FALSE} then no text output is shown.}

\item{fGroup}{The name (thus a character) of the feature group that should be
searched for.}

\item{index}{The index of the component. Can be a numeric index or a
character with its name.}

\item{markFGroup}{If specified (\emph{i.e.} not \code{NULL}) this argument
can be used to mark a feature group in the plotted spectrum. The value
should be a character with the name of the feature group. Setting this to
\code{NULL} will not mark any peak.}

\item{xlim, ylim}{Sets the plot size limits used by
\code{\link[graphics]{plot}}. Set to \code{NULL} for automatic plot sizing.}

\item{fGroups}{The \code{\link{featureGroups}} object that was used to
generate the components.}

\item{rtWindow}{Retention window: see the \code{plotChroms} method for the
\code{\link{featureGroups}} class.}

\item{sets}{\setsWF A \code{character} with name(s) of the sets to keep (or remove if \code{negate=TRUE}).}

\item{set}{\setsWF The name of the set.}
}
\value{
\code{delete} returns the object for which the specified data was removed.

\code{consensus} returns a \code{components} object that is produced
  by merging multiple specified \code{components} objects.
}
\description{
Contains data for feature groups that are related in some way. These
\emph{components} commonly include adducts, isotopes and homologues.
}
\details{
\code{components} objects are obtained from \code{\link{generateComponents}}.
}
\section{Methods (by generic)}{
\itemize{
\item \code{componentTable}: Accessor method for the \code{components} slot of a
\code{components} class. Each component is stored as a
\code{\link{data.table}}.

\item \code{componentInfo}: Accessor method for the \code{componentInfo} slot of a
\code{components} class.

\item \code{groupNames}: returns a \code{character} vector with the names of the
feature groups for which data is present in this object.

\item \code{length}: Obtain total number of components.

\item \code{names}: Obtain the names of all components.

\item \code{show}: Show summary information for this object.

\item \code{[}: Subset on components/feature groups.

\item \code{[[}: Extracts a component table, optionally filtered by a feature group.

\item \code{$}: Extracts a component table by component name.

\item \code{delete}: Completely deletes specified (parts of) components.

\item \code{as.data.table}: Returns all component data in a table.

\item \code{filter}: Provides rule based filtering for components.

\item \code{findFGroup}: Returns the component id(s) to which a feature group
belongs.

\item \code{plotSpectrum}: Plot a \emph{pseudo} mass spectrum for a single
component.

\item \code{plotChroms}: Plot an extracted ion chromatogram (EIC) for all
feature groups within a single component.

\item \code{consensus}: Generates a consensus from multiple \code{components}
objects. At this point results are simply combined and no attempt is made to
merge similar components.
}}

\section{Slots}{

\describe{
\item{\code{components}}{List of all components in this object. Use the
\code{componentTable} method for access.}

\item{\code{componentInfo}}{A \code{\link{data.table}} containing general information
for each component. Use the \code{componentInfo} method for access.}
}}

\note{
\code{filter} Applies only those filters for which a component has data available. For instance, filtering by
  adduct will only filter any results within a component if that component contains adduct information.
}
\section{S4 class hierarchy}{
 \itemize{   \item{\code{\link{workflowStep}}}   \itemize{     \item{\strong{\code{\link{components}}}}     \itemize{       \item{\code{\link{componentsCamera}}}       \item{\code{\link{componentsFeatures}}}       \itemize{         \item{\code{\link{componentsCliqueMS}}}         \item{\code{\link{componentsOpenMS}}}       }       \item{\code{\link{componentsClust}}}       \itemize{         \item{\code{\link{componentsIntClust}}}         \item{\code{\link{componentsSpecClust}}}       }       \item{\code{\link{componentsSet}}}       \itemize{         \item{\code{\link{componentsNTSet}}}       }       \item{\code{\link{componentsUnset}}}       \item{\code{\link{componentsNT}}}       \itemize{         \item{\code{\link{componentsNTUnset}}}       }       \item{\code{\link{componentsRC}}}       \item{\code{\link{componentsTPs}}}     }   } }
}

\section{Sets workflows}{
 \setsWFClass{componentsSet}{components}

  \setsWFNewMethodsSO{componentsUnset}{Only the components in the specified set are kept.}

  \setsWFChangedMethods{

  \item \code{filter} and the subset operator (\code{[}) Can be used to select components that are only present for
  selected sets.

  }
}

\seealso{
\code{\link{generateComponents}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{generics}
\alias{generics}
\alias{adducts}
\alias{adducts<-}
\alias{algorithm}
\alias{analysisInfo}
\alias{analyses}
\alias{annotatedPeakList}
\alias{annotations}
\alias{calculatePeakQualities}
\alias{clusterProperties}
\alias{clusters}
\alias{consensus}
\alias{convertToMFDB}
\alias{convertToSuspects}
\alias{cutClusters}
\alias{defaultExclNormScores}
\alias{export}
\alias{featureTable}
\alias{filter}
\alias{getFeatures}
\alias{getMCS}
\alias{groupNames}
\alias{plotChord}
\alias{plotChroms}
\alias{plotGraph}
\alias{plotInt}
\alias{plotScores}
\alias{plotSilhouettes}
\alias{plotSpectrum}
\alias{plotStructure}
\alias{plotVenn}
\alias{plotUpSet}
\alias{delete}
\alias{plotVolcano}
\alias{replicateGroups}
\alias{setObjects}
\alias{sets}
\alias{treeCut}
\alias{treeCutDynamic}
\alias{unset}
\alias{[}
\alias{[[}
\alias{$}
\alias{as.data.table}
\alias{as.data.frame}
\alias{length}
\alias{lengths}
\alias{names}
\alias{plot}
\alias{show}
\title{Miscellaneous generics}
\usage{
adducts(obj, ...)

adducts(obj, ...) <- value

algorithm(obj)

analysisInfo(obj)

analyses(obj)

annotatedPeakList(obj, ...)

annotations(obj, ...)

calculatePeakQualities(obj, weights = NULL, flatnessFactor = 0.05, ...)

clusterProperties(obj)

clusters(obj)

consensus(obj, ...)

convertToMFDB(TPs, out, ...)

convertToSuspects(TPs, ...)

cutClusters(obj)

defaultExclNormScores(obj)

export(obj, type, out, ...)

featureTable(obj, ...)

filter(obj, ...)

getFeatures(obj)

getMCS(obj, ...)

groupNames(obj)

plotChord(obj, addSelfLinks = FALSE, addRetMzPlots = TRUE, ...)

plotChroms(obj, ...)

plotGraph(obj, onlyLinked = TRUE, ...)

plotInt(obj, ...)

plotScores(obj, ...)

plotSilhouettes(obj, kSeq, ...)

plotSpectrum(obj, ...)

plotStructure(obj, ...)

plotVenn(obj, ...)

plotUpSet(obj, ...)

delete(obj, ...)

plotVolcano(obj, ...)

replicateGroups(obj)

setObjects(obj)

sets(obj)

treeCut(obj, k = NULL, h = NULL, ...)

treeCutDynamic(
  obj,
  maxTreeHeight = 1,
  deepSplit = TRUE,
  minModuleSize = 1,
  ...
)

unset(obj, set)
}
\arguments{
\item{obj}{The object the generic should be applied to.}

\item{\dots}{Any further method specific arguments. See method documentation
for details.}

\item{value}{The replacement value.}

\item{weights, flatnessFactor}{See method documentation.}

\item{TPs}{The \code{\link{transformationProducts}} derived object.}

\item{out}{Output file.}

\item{type}{The export type.}

\item{addSelfLinks}{If \code{TRUE} then 'self-links' are added which
represent non-shared data.}

\item{addRetMzPlots}{Set to \code{TRUE} to enable \emph{m/z} \emph{vs}
retention time scatter plots.}

\item{onlyLinked}{Only plots linked objects if \code{TRUE}.}

\item{kSeq}{An integer vector containing the sequence that should be used for
average silhouette width calculation.}

\item{k, h}{Desired numbers of clusters. See \code{\link{cutree}}.}

\item{maxTreeHeight, deepSplit, minModuleSize}{Arguments used by
\code{\link{cutreeDynamicTree}}.}

\item{set}{The name of the set.}
}
\description{
Various (S4) generic functions providing a common interface for common tasks
such as plotting and filtering data. The actual functionality and function
arguments are often specific for the implemented methods, for this reason,
please refer to the linked method documentation for each generic.
}
\details{
\code{adducts} returns assigned adducts of the object.

\itemize{
    \item Methods are defined for: \code{\link[=adducts,featureGroups-method]{featureGroups}}; \code{\link[=adducts,featureGroupsSet-method]{featureGroupsSet}}.
}

\code{adducts<-} sets adducts of the object.

\itemize{
    \item Methods are defined for: \code{\link[=adducts<-,featureGroups-method]{featureGroups}}; \code{\link[=adducts<-,featureGroupsSet-method]{featureGroupsSet}}.
}

\code{algorithm} returns the algorithm that was used to generate the object.

\itemize{
    \item Methods are defined for: \code{\link[=algorithm,optimizationResult-method]{optimizationResult}}; \code{\link[=algorithm,workflowStep-method]{workflowStep}}.
}

\code{analysisInfo} returns the \link[=analysis-information]{analysis information} from an object.

\itemize{
    \item Methods are defined for: \code{\link[=analysisInfo,featureGroups-method]{featureGroups}}; \code{\link[=analysisInfo,features-method]{features}}; \code{\link[=analysisInfo,MSPeakListsSet-method]{MSPeakListsSet}}.
}

\code{analyses} returns a \code{character} vector with the analyses for which data is present in this object.

\itemize{
    \item Methods are defined for: \code{\link[=analyses,featureGroups-method]{featureGroups}}; \code{\link[=analyses,features-method]{features}}; \code{\link[=analyses,formulas-method]{formulas}}; \code{\link[=analyses,MSPeakLists-method]{MSPeakLists}}.
}

\code{annotatedPeakList} returns an annotated MS peak list.

\itemize{
    \item Methods are defined for: \code{\link[=annotatedPeakList,compounds-method]{compounds}}; \code{\link[=annotatedPeakList,compoundsSet-method]{compoundsSet}}; \code{\link[=annotatedPeakList,formulas-method]{formulas}}; \code{\link[=annotatedPeakList,formulasSet-method]{formulasSet}}.
}

\code{annotations} returns annotations.

\itemize{
    \item Methods are defined for: \code{\link[=annotations,featureAnnotations-method]{featureAnnotations}}; \code{\link[=annotations,featureGroups-method]{featureGroups}}; \code{\link[=annotations,formulas-method]{formulas}}.
}

\code{calculatePeakQualities} calculates chromatographic peak qualities and scores.

\itemize{
    \item Methods are defined for: \code{\link[=calculatePeakQualities,featureGroups-method]{featureGroups}}; \code{\link[=calculatePeakQualities,features-method]{features}}.
}

\code{clusterProperties} Obtain a list with properties of the generated cluster(s).

\itemize{
    \item Methods are defined for: \code{\link[=clusterProperties,componentsClust-method]{componentsClust}}; \code{\link[=clusterProperties,compoundsCluster-method]{compoundsCluster}}.
}

\code{clusters} Obtain clustering object(s).

\itemize{
    \item Methods are defined for: \code{\link[=clusters,componentsClust-method]{componentsClust}}; \code{\link[=clusters,compoundsCluster-method]{compoundsCluster}}.
}

\code{consensus} combines and merges data from various algorithms to generate a consensus.

\itemize{
    \item Methods are defined for: \code{\link[=consensus,components-method]{components}}; \code{\link[=consensus,componentsSet-method]{componentsSet}}; \code{\link[=consensus,compounds-method]{compounds}}; \code{\link[=consensus,compoundsSet-method]{compoundsSet}}; \code{\link[=consensus,featureGroupsComparison-method]{featureGroupsComparison}}; \code{\link[=consensus,featureGroupsComparisonSet-method]{featureGroupsComparisonSet}}; \code{\link[=consensus,formulas-method]{formulas}}; \code{\link[=consensus,formulasSet-method]{formulasSet}}.
}

\code{convertToMFDB} Exports the object to a local database that can be used with \command{MetFrag}.

\itemize{
    \item Methods are defined for: .
}

\code{convertToSuspects} Converts an object to a suspect list.

\itemize{
    \item Methods are defined for: .
}

\code{cutClusters} Returns assigned cluster indices of a cut cluster.

\itemize{
    \item Methods are defined for: \code{\link[=cutClusters,componentsClust-method]{componentsClust}}; \code{\link[=cutClusters,compoundsCluster-method]{compoundsCluster}}.
}

\code{defaultExclNormScores} Returns default scorings that are excluded from normalization.

\itemize{
    \item Methods are defined for: \code{\link[=defaultExclNormScores,compounds-method]{compounds}}; \code{\link[=defaultExclNormScores,formulas-method]{formulas}}.
}

\code{export} exports workflow data to a given format.

\itemize{
    \item Methods are defined for: \code{\link[=export,featureGroups-method]{featureGroups}}; \code{\link[=export,featureGroupsSet-method]{featureGroupsSet}}.
}

\code{featureTable} returns feature information.

\itemize{
    \item Methods are defined for: \code{\link[=featureTable,featureGroups-method]{featureGroups}}; \code{\link[=featureTable,featureGroupsSet-method]{featureGroupsSet}}; \code{\link[=featureTable,features-method]{features}}.
}

\code{filter} provides various functionality to do post-filtering of data.

\itemize{
    \item Methods are defined for: \code{\link[=filter,components-method]{components}}; \code{\link[=filter,componentsSet-method]{componentsSet}}; \code{\link[=filter,componentsTPs-method]{componentsTPs}}; \code{\link[=filter,compounds-method]{compounds}}; \code{\link[=filter,compoundsSet-method]{compoundsSet}}; \code{\link[=filter,featureAnnotations-method]{featureAnnotations}}; \code{\link[=filter,featureGroups-method]{featureGroups}}; \code{\link[=filter,featureGroupsScreening-method]{featureGroupsScreening}}; \code{\link[=filter,featureGroupsScreeningSet-method]{featureGroupsScreeningSet}}; \code{\link[=filter,featureGroupsSet-method]{featureGroupsSet}}; \code{\link[=filter,features-method]{features}}; \code{\link[=filter,featuresSet-method]{featuresSet}}; \code{\link[=filter,formulasSet-method]{formulasSet}}; \code{\link[=filter,MSPeakLists-method]{MSPeakLists}}; \code{\link[=filter,MSPeakListsSet-method]{MSPeakListsSet}}; \code{\link[=filter,transformationProductsBT-method]{transformationProductsBT}}.
}

\code{getFeatures} returns the object's \code{\link{features}} object.

\itemize{
    \item Methods are defined for: \code{\link[=getFeatures,featureGroups-method]{featureGroups}}.
}

\code{getMCS} Calculcates the maximum common substructure.

\itemize{
    \item Methods are defined for: \code{\link[=getMCS,compounds-method]{compounds}}; \code{\link[=getMCS,compoundsCluster-method]{compoundsCluster}}.
}

\code{groupNames} returns a \code{character} vector with the names of the feature groups for which data is present in this object.

\itemize{
    \item Methods are defined for: \code{\link[=groupNames,components-method]{components}}; \code{\link[=groupNames,compoundsCluster-method]{compoundsCluster}}; \code{\link[=groupNames,featureAnnotations-method]{featureAnnotations}}; \code{\link[=groupNames,featureGroups-method]{featureGroups}}; \code{\link[=groupNames,MSPeakLists-method]{MSPeakLists}}.
}

\code{plotChord} plots a Chord diagram to assess overlapping data.

\itemize{
    \item Methods are defined for: \code{\link[=plotChord,featureGroups-method]{featureGroups}}; \code{\link[=plotChord,featureGroupsComparison-method]{featureGroupsComparison}}.
}

\code{plotChroms} plots extracted ion chromatogram(s).

\itemize{
    \item Methods are defined for: \code{\link[=plotChroms,components-method]{components}}; \code{\link[=plotChroms,featureGroups-method]{featureGroups}}.
}

\code{plotGraph} Plots an interactive network graph.

\itemize{
    \item Methods are defined for: \code{\link[=plotGraph,componentsNT-method]{componentsNT}}; \code{\link[=plotGraph,componentsNTSet-method]{componentsNTSet}}; \code{\link[=plotGraph,componentsTPs-method]{componentsTPs}}.
}

\code{plotInt} plots the intensity of all contained features.

\itemize{
    \item Methods are defined for: \code{\link[=plotInt,componentsIntClust-method]{componentsIntClust}}; \code{\link[=plotInt,featureGroups-method]{featureGroups}}; \code{\link[=plotInt,featureGroupsSet-method]{featureGroupsSet}}.
}

\code{plotScores} plots candidate scorings.

\itemize{
    \item Methods are defined for: \code{\link[=plotScores,compounds-method]{compounds}}; \code{\link[=plotScores,formulas-method]{formulas}}.
}

\code{plotSilhouettes} plots silhouette widths to evaluate the desired cluster size.

\itemize{
    \item Methods are defined for: \code{\link[=plotSilhouettes,componentsClust-method]{componentsClust}}; \code{\link[=plotSilhouettes,compoundsCluster-method]{compoundsCluster}}.
}

\code{plotSpectrum} plots a (annotated) spectrum.

\itemize{
    \item Methods are defined for: \code{\link[=plotSpectrum,components-method]{components}}; \code{\link[=plotSpectrum,compounds-method]{compounds}}; \code{\link[=plotSpectrum,compoundsSet-method]{compoundsSet}}; \code{\link[=plotSpectrum,formulas-method]{formulas}}; \code{\link[=plotSpectrum,formulasSet-method]{formulasSet}}; \code{\link[=plotSpectrum,MSPeakLists-method]{MSPeakLists}}; \code{\link[=plotSpectrum,MSPeakListsSet-method]{MSPeakListsSet}}.
}

\code{plotStructure} plots a chemical structure.

\itemize{
    \item Methods are defined for: \code{\link[=plotStructure,compounds-method]{compounds}}; \code{\link[=plotStructure,compoundsCluster-method]{compoundsCluster}}.
}

\code{plotVenn} plots a Venn diagram to assess unique and overlapping data.

\itemize{
    \item Methods are defined for: \code{\link[=plotVenn,featureAnnotations-method]{featureAnnotations}}; \code{\link[=plotVenn,featureGroups-method]{featureGroups}}; \code{\link[=plotVenn,featureGroupsComparison-method]{featureGroupsComparison}}; \code{\link[=plotVenn,featureGroupsSet-method]{featureGroupsSet}}.
}

\code{plotUpSet} plots an UpSet diagram to assess unique and overlapping data.

\itemize{
    \item Methods are defined for: \code{\link[=plotUpSet,featureAnnotations-method]{featureAnnotations}}; \code{\link[=plotUpSet,featureGroups-method]{featureGroups}}; \code{\link[=plotUpSet,featureGroupsComparison-method]{featureGroupsComparison}}.
}

\code{delete} Deletes results.

\itemize{
    \item Methods are defined for: \code{\link[=delete,components-method]{components}}; \code{\link[=delete,componentsClust-method]{componentsClust}}; \code{\link[=delete,componentsSet-method]{componentsSet}}; \code{\link[=delete,compoundsSet-method]{compoundsSet}}; \code{\link[=delete,featureAnnotations-method]{featureAnnotations}}; \code{\link[=delete,featureGroups-method]{featureGroups}}; \code{\link[=delete,featureGroupsKPIC2-method]{featureGroupsKPIC2}}; \code{\link[=delete,featureGroupsScreening-method]{featureGroupsScreening}}; \code{\link[=delete,featureGroupsScreeningSet-method]{featureGroupsScreeningSet}}; \code{\link[=delete,featureGroupsSet-method]{featureGroupsSet}}; \code{\link[=delete,featureGroupsXCMS-method]{featureGroupsXCMS}}; \code{\link[=delete,featureGroupsXCMS3-method]{featureGroupsXCMS3}}; \code{\link[=delete,features-method]{features}}; \code{\link[=delete,featuresKPIC2-method]{featuresKPIC2}}; \code{\link[=delete,featuresXCMS-method]{featuresXCMS}}; \code{\link[=delete,featuresXCMS3-method]{featuresXCMS3}}; \code{\link[=delete,formulas-method]{formulas}}; \code{\link[=delete,formulasSet-method]{formulasSet}}.
}

\code{plotVolcano} plots a volcano plot.

\itemize{
    \item Methods are defined for: \code{\link[=plotVolcano,featureGroups-method]{featureGroups}}.
}

\code{replicateGroups} returns a \code{character} vector with the analyses for which data is present in this object.

\itemize{
    \item Methods are defined for: \code{\link[=replicateGroups,featureGroups-method]{featureGroups}}; \code{\link[=replicateGroups,features-method]{features}}.
}

\code{setObjects} returns the \emph{set objects} of this object. See the documentation of \code{\link{workflowStepSet}}.

\itemize{
    \item Methods are defined for: \code{\link[=setObjects,workflowStepSet-method]{workflowStepSet}}.
}

\code{sets} returns the names of the sets inside this object. See the documentation for \link[=sets-workflow]{sets workflows}.

\itemize{
    \item Methods are defined for: \code{\link[=sets,featureGroupsSet-method]{featureGroupsSet}}; \code{\link[=sets,featuresSet-method]{featuresSet}}; \code{\link[=sets,workflowStepSet-method]{workflowStepSet}}.
}

\code{treeCut} Manually cut a cluster.

\itemize{
    \item Methods are defined for: \code{\link[=treeCut,componentsClust-method]{componentsClust}}; \code{\link[=treeCut,compoundsCluster-method]{compoundsCluster}}.
}

\code{treeCutDynamic} Automatically cut a cluster.

\itemize{
    \item Methods are defined for: \code{\link[=treeCutDynamic,componentsClust-method]{componentsClust}}; \code{\link[=treeCutDynamic,compoundsCluster-method]{compoundsCluster}}.
}

\code{unset} Converts this object to a regular non-set object. See the documentation for \link[=sets-workflow]{sets workflows}.

\itemize{
    \item Methods are defined for: \code{\link[=unset,componentsNTSet-method]{componentsNTSet}}; \code{\link[=unset,componentsSet-method]{componentsSet}}; \code{\link[=unset,compoundsConsensusSet-method]{compoundsConsensusSet}}; \code{\link[=unset,compoundsSet-method]{compoundsSet}}; \code{\link[=unset,featureGroupsScreeningSet-method]{featureGroupsScreeningSet}}; \code{\link[=unset,featureGroupsSet-method]{featureGroupsSet}}; \code{\link[=unset,featuresSet-method]{featuresSet}}; \code{\link[=unset,formulasConsensusSet-method]{formulasConsensusSet}}; \code{\link[=unset,formulasSet-method]{formulasSet}}; \code{\link[=unset,MSPeakListsSet-method]{MSPeakListsSet}}.
}
}
\section{Other generics}{
 Below are methods that are defined for existing
  generics (\emph{e.g.} defined in \code{base}). Please see method specific
  documentation for more details.

 \code{[} Subsets data within an object.

\itemize{
    \item Methods are defined for: \code{\link[=[,components,ANY,ANY,missing-method]{components,ANY,ANY,missing}}; \code{\link[=[,componentsSet,ANY,ANY,missing-method]{componentsSet,ANY,ANY,missing}}; \code{\link[=[,compoundsCluster,ANY,missing,missing-method]{compoundsCluster,ANY,missing,missing}}; \code{\link[=[,compoundsSet,ANY,missing,missing-method]{compoundsSet,ANY,missing,missing}}; \code{\link[=[,featureAnnotations,ANY,missing,missing-method]{featureAnnotations,ANY,missing,missing}}; \code{\link[=[,featureGroups,ANY,ANY,missing-method]{featureGroups,ANY,ANY,missing}}; \code{\link[=[,featureGroupsComparison,ANY,missing,missing-method]{featureGroupsComparison,ANY,missing,missing}}; \code{\link[=[,featureGroupsScreening,ANY,ANY,missing-method]{featureGroupsScreening,ANY,ANY,missing}}; \code{\link[=[,featureGroupsScreeningSet,ANY,ANY,missing-method]{featureGroupsScreeningSet,ANY,ANY,missing}}; \code{\link[=[,featureGroupsSet,ANY,ANY,missing-method]{featureGroupsSet,ANY,ANY,missing}}; \code{\link[=[,features,ANY,missing,missing-method]{features,ANY,missing,missing}}; \code{\link[=[,featuresSet,ANY,missing,missing-method]{featuresSet,ANY,missing,missing}}; \code{\link[=[,formulasSet,ANY,missing,missing-method]{formulasSet,ANY,missing,missing}}; \code{\link[=[,MSPeakLists,ANY,ANY,missing-method]{MSPeakLists,ANY,ANY,missing}}; \code{\link[=[,MSPeakListsSet,ANY,ANY,missing-method]{MSPeakListsSet,ANY,ANY,missing}}; \code{\link[=[,transformationProducts,ANY,missing,missing-method]{transformationProducts,ANY,missing,missing}}.
}

 \code{[[} Extract data from an object.

\itemize{
    \item Methods are defined for: \code{\link[=[[,components,ANY,ANY-method]{components,ANY,ANY}}; \code{\link[=[[,featureAnnotations,ANY,missing-method]{featureAnnotations,ANY,missing}}; \code{\link[=[[,featureGroups,ANY,ANY-method]{featureGroups,ANY,ANY}}; \code{\link[=[[,featureGroupsComparison,ANY,missing-method]{featureGroupsComparison,ANY,missing}}; \code{\link[=[[,features,ANY,missing-method]{features,ANY,missing}}; \code{\link[=[[,formulas,ANY,ANY-method]{formulas,ANY,ANY}}; \code{\link[=[[,MSPeakLists,ANY,ANY-method]{MSPeakLists,ANY,ANY}}; \code{\link[=[[,transformationProducts,ANY,missing-method]{transformationProducts,ANY,missing}}.
}

 \code{$} Extract data from an object.

\itemize{
    \item Methods are defined for: \code{\link[=$,components-method]{components}}; \code{\link[=$,featureAnnotations-method]{featureAnnotations}}; \code{\link[=$,featureGroups-method]{featureGroups}}; \code{\link[=$,featureGroupsComparison-method]{featureGroupsComparison}}; \code{\link[=$,features-method]{features}}; \code{\link[=$,MSPeakLists-method]{MSPeakLists}}; \code{\link[=$,transformationProducts-method]{transformationProducts}}.
}

 \code{as.data.table} Converts an object to a table (\code{\link{data.table}}).

\itemize{
    \item Methods are defined for: \code{\link[=as.data.table,components-method]{components}}; \code{\link[=as.data.table,componentsTPs-method]{componentsTPs}}; \code{\link[=as.data.table,featureAnnotations-method]{featureAnnotations}}; \code{\link[=as.data.table,featureGroups-method]{featureGroups}}; \code{\link[=as.data.table,featureGroupsScreening-method]{featureGroupsScreening}}; \code{\link[=as.data.table,featureGroupsScreeningSet-method]{featureGroupsScreeningSet}}; \code{\link[=as.data.table,featureGroupsSet-method]{featureGroupsSet}}; \code{\link[=as.data.table,features-method]{features}}; \code{\link[=as.data.table,featuresSet-method]{featuresSet}}; \code{\link[=as.data.table,formulas-method]{formulas}}; \code{\link[=as.data.table,MSPeakLists-method]{MSPeakLists}}; \code{\link[=as.data.table,MSPeakListsSet-method]{MSPeakListsSet}}; \code{\link[=as.data.table,transformationProducts-method]{transformationProducts}}; \code{\link[=as.data.table,workflowStep-method]{workflowStep}}.
}

 \code{as.data.frame} Converts an object to a table (\code{data.frame}).

\itemize{
    \item Methods are defined for: \code{\link[=as.data.frame,workflowStep-method]{workflowStep}}.
}

 \code{length} Returns the length of an object.

\itemize{
    \item Methods are defined for: \code{\link[=length,components-method]{components}}; \code{\link[=length,compoundsCluster-method]{compoundsCluster}}; \code{\link[=length,featureAnnotations-method]{featureAnnotations}}; \code{\link[=length,featureGroups-method]{featureGroups}}; \code{\link[=length,featureGroupsComparison-method]{featureGroupsComparison}}; \code{\link[=length,features-method]{features}}; \code{\link[=length,MSPeakLists-method]{MSPeakLists}}; \code{\link[=length,optimizationResult-method]{optimizationResult}}; \code{\link[=length,transformationProducts-method]{transformationProducts}}.
}

 \code{lengths} Returns the lengths of elements within this object.

\itemize{
    \item Methods are defined for: \code{\link[=lengths,compoundsCluster-method]{compoundsCluster}}; \code{\link[=lengths,optimizationResult-method]{optimizationResult}}.
}

 \code{names} Return names for this object.

\itemize{
    \item Methods are defined for: \code{\link[=names,components-method]{components}}; \code{\link[=names,featureGroups-method]{featureGroups}}; \code{\link[=names,featureGroupsComparison-method]{featureGroupsComparison}}; \code{\link[=names,transformationProducts-method]{transformationProducts}}.
}

 \code{plot} Generates a plot for an object.

\itemize{
    \item Methods are defined for: \code{\link[=plot,componentsClust,missing-method]{componentsClust,missing}}; \code{\link[=plot,compoundsCluster,missing-method]{compoundsCluster,missing}}; \code{\link[=plot,featureGroups,missing-method]{featureGroups,missing}}; \code{\link[=plot,featureGroupsComparison,missing-method]{featureGroupsComparison,missing}}; \code{\link[=plot,optimizationResult,missing-method]{optimizationResult,missing}}.
}

 \code{show} Prints information about this object.

\itemize{
    \item Methods are defined for: \code{\link[=show,adduct-method]{adduct}}; \code{\link[=show,components-method]{components}}; \code{\link[=show,componentsFeatures-method]{componentsFeatures}}; \code{\link[=show,componentsSet-method]{componentsSet}}; \code{\link[=show,compounds-method]{compounds}}; \code{\link[=show,compoundsCluster-method]{compoundsCluster}}; \code{\link[=show,compoundsSet-method]{compoundsSet}}; \code{\link[=show,featureGroups-method]{featureGroups}}; \code{\link[=show,featureGroupsScreening-method]{featureGroupsScreening}}; \code{\link[=show,featureGroupsScreeningSet-method]{featureGroupsScreeningSet}}; \code{\link[=show,featureGroupsSet-method]{featureGroupsSet}}; \code{\link[=show,features-method]{features}}; \code{\link[=show,featuresSet-method]{featuresSet}}; \code{\link[=show,formulas-method]{formulas}}; \code{\link[=show,formulasSet-method]{formulasSet}}; \code{\link[=show,MSPeakLists-method]{MSPeakLists}}; \code{\link[=show,MSPeakListsSet-method]{MSPeakListsSet}}; \code{\link[=show,optimizationResult-method]{optimizationResult}}; \code{\link[=show,transformationProducts-method]{transformationProducts}}; \code{\link[=show,workflowStep-method]{workflowStep}}; \code{\link[=show,workflowStepSet-method]{workflowStepSet}}.
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/components-camera.R
\name{generateComponentsCAMERA}
\alias{generateComponentsCAMERA}
\alias{generateComponentsCAMERA,featureGroups-method}
\alias{generateComponentsCAMERA,featureGroupsSet-method}
\title{Componentization of adducts, isotopes etc. with CAMERA}
\usage{
\S4method{generateComponentsCAMERA}{featureGroups}(
  fGroups,
  ionization = NULL,
  onlyIsotopes = FALSE,
  minSize = 2,
  relMinReplicates = 0.5,
  extraOpts = NULL
)

\S4method{generateComponentsCAMERA}{featureGroupsSet}(fGroups, ionization = NULL, ...)
}
\arguments{
\item{fGroups}{\code{\link{featureGroups}} object for which components should be generated.}

\item{ionization}{Which ionization polarity was used to generate the data: should be \code{"positive"}
  or \code{"negative"}. If the \code{featureGroups} object has adduct annotations, and \code{ionization=NULL}, the
  ionization will be detected automatically.

  \setsWF This parameter is not supported for sets workflows, as the ionization will always be detected
  automatically.}

\item{onlyIsotopes}{Logical value. If \code{TRUE} only isotopes are considered when generating components (faster).
Corresponds to \code{quick} argument of \code{\link[CAMERA:annotate-methods]{CAMERA::annotate}}.}

\item{minSize}{The minimum size of a component. Smaller components than this size will be removed. See note below.}

\item{relMinReplicates}{Feature groups within a component are only kept when they contain data for at least this
(relative) amount of replicate analyses. For instance, \samp{0.5} means that at least half of the replicates should
contain data for a particular feature group in a component. In this calculation replicates that are fully absent
within a component are not taken in to account. See note below.}

\item{extraOpts}{Named character vector with extra arguments directly passed to
\code{\link[CAMERA:annotate-methods]{CAMERA::annotate}}. Set to \code{NULL} to ignore.}

\item{\dots}{\setsWF Further arguments passed to the non-sets workflow method.}
}
\value{
A \code{\link{components}} (derived) object containing all generated components.
}
\description{
Interfaces with \href{https://bioconductor.org/packages/release/bioc/html/CAMERA.html}{CAMERA} to generate components
from known adducts, isotopes and in-source fragments.
}
\details{
This function uses CAMERA to generate components. This function is called when calling \code{generateComponents} with
  \code{algorithm="camera"}.

The specified \code{featureGroups} object is automatically converted to an \code{\link{xcmsSet}} object
  using \code{\link{getXCMSSet}}.
}
\note{
The default value for \code{minSize} and \code{relMinReplicates} results in
  extra filtering, hence, the final results may be different than what the algorithm normally would return.
}
\section{Sets workflows}{
 In a \link[=sets-workflow]{sets workflow} the componentization is first performed for each
  set independently. The resulting components are then all combined in a \code{\link{componentsSet}} object. Note that
  the components themselves are never merged. The components are renamed to include the set name from which they were
  generated (\emph{e.g.} \code{"CMP1"} becomes \code{"CMP1-positive"}).
}

\references{
\addCitations{CAMERA}{1}
}
\seealso{
\code{\link{generateComponents}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TP-library.R
\docType{class}
\name{transformationProductsLibrary-class}
\alias{transformationProductsLibrary-class}
\alias{transformationProductsLibrary}
\alias{convertToMFDB,transformationProductsLibrary-method}
\title{Class to store transformation products (TPs) obtained from a library}
\usage{
\S4method{convertToMFDB}{transformationProductsLibrary}(TPs, out, includeParents = FALSE)
}
\arguments{
\item{TPs}{\code{transformationProductsLibrary} object to be accessed}

\item{out}{The file name of the the output \command{MetFrag} database.}

\item{includeParents}{Set to \code{TRUE} to include the parents in the database.}
}
\description{
This class is used to store prediction results that are available in a TP library.
}
\details{
Objects from this class are generate with \code{\link{generateTPsLibrary}}. This class is derived from the
\code{\link{transformationProducts}} base class, please see its documentation for more details.
}
\section{Methods (by generic)}{
\itemize{
\item \code{convertToMFDB}: Exports this object as a \file{.csv} file that can be used as a \command{MetFrag} local
database.
}}

\section{S4 class hierarchy}{
 \itemize{   \item{\code{\link{transformationProducts}}}   \itemize{     \item{\strong{\code{\link{transformationProductsLibrary}}}}   } }
}

\seealso{
The base class \code{\link{transformationProducts}} for more relevant methods and
  \code{\link{generateTPs}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{clearCache}
\alias{clearCache}
\title{Clearing of cached data}
\usage{
clearCache(what = NULL, file = NULL)
}
\arguments{
\item{what}{This argument describes what should be done. When \code{what =
NULL} this function will list which tables are present along with an
indication of their size (database rows). If \code{what = "all"} then the
complete file will be removed. Otherwise, \code{what} should be a character
string (a regular expression) that is used to match the table names that
should be removed.}

\item{file}{The cache file. If \code{NULL} then the value of the
\code{patRoon.cache.fileName} option is used.}
}
\description{
Remove (part of) the cache database used to store (intermediate) processing
results.
}
\details{
This function will either remove one or more tables within the cache
\code{sqlite} database or simply wipe the whole cache file. Removing tables
will \code{VACUUM} the database, which may take some time for large cache
files.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups-tasq.R
\name{importFeatureGroupsBrukerTASQ}
\alias{importFeatureGroupsBrukerTASQ}
\title{Imports feature groups from Bruker TASQ}
\usage{
importFeatureGroupsBrukerTASQ(path, analysisInfo, clusterRTWindow = 12)
}
\arguments{
\item{path}{The file path to an Excel export of the Global results table from TASQ, converted to \file{.csv} format.}

\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}

\item{clusterRTWindow}{This retention time window (in seconds) is used to group hits across analyses together. See
also the details section.}
}
\value{
A new \code{featureGroups} object containing converted screening results from Bruker TASQ.
}
\description{
Imports screening results from Bruker TASQ as feature groups.
}
\details{
This function imports data from Bruker TASQ. This function is called when calling \code{importFeatureGroups} with
  \code{type="brukertasq"}.

The feature groups across analyses are formed based on the name of suspects and their closeness in retention
  time. The latter is necessary because TASQ does not necessarily perform checks on retention times and may therefore
  assign a suspect to peaks with different retention times across analyses (or within a single analysis). Hence,
  suspects with equal names are hierarchically clustered on their retention times (using \pkg{\link{fastcluster}}) to
  form the feature groups. The cut-off value for this is specified by the \code{clusterRTWindow} argument. The input
  for this function is obtained by generating an Excel export of the 'global' results and subsequently converting the
  file to \file{.csv} format.
}
\note{
This function uses estimated min/max values for retention times and dummy min/max \emph{m/z} values for
  conversion to features, since this information is not (readily) available. Hence, when plotting, for instance,
  extracted ion chromatograms (with \code{\link{plotChroms}}) the integrated chromatographic peak range shown is
  incorrect.

  This function may use suspect names to base file names used for reporting, logging etc. Therefore, it is important
  that these are file-compatible names.
}
\references{
\addCitations{fastcluster}{1}
}
\seealso{
\code{\link{importFeatureGroups}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mspeaklists.R, R/mspeaklists-set.R,
%   R/utils-mspeaklists.R
\docType{class}
\name{MSPeakLists-class}
\alias{MSPeakLists-class}
\alias{MSPeakLists}
\alias{peakLists,MSPeakLists-method}
\alias{peakLists}
\alias{averagedPeakLists,MSPeakLists-method}
\alias{averagedPeakLists}
\alias{analyses,MSPeakLists-method}
\alias{groupNames,MSPeakLists-method}
\alias{length,MSPeakLists-method}
\alias{show,MSPeakLists-method}
\alias{[,MSPeakLists,ANY,ANY,missing-method}
\alias{[[,MSPeakLists,ANY,ANY-method}
\alias{$,MSPeakLists-method}
\alias{as.data.table,MSPeakLists-method}
\alias{filter,MSPeakLists-method}
\alias{plotSpectrum,MSPeakLists-method}
\alias{spectrumSimilarity,MSPeakLists-method}
\alias{spectrumSimilarity}
\alias{MSPeakListsSet-class}
\alias{MSPeakListsSet}
\alias{analysisInfo,MSPeakListsSet-method}
\alias{show,MSPeakListsSet-method}
\alias{[,MSPeakListsSet,ANY,ANY,missing-method}
\alias{as.data.table,MSPeakListsSet-method}
\alias{filter,MSPeakListsSet-method}
\alias{plotSpectrum,MSPeakListsSet-method}
\alias{spectrumSimilarity,MSPeakListsSet-method}
\alias{MSPeakListsUnset-class}
\alias{MSPeakListsUnset}
\alias{unset,MSPeakListsSet-method}
\alias{getDefIsolatePrecParams}
\title{Class containing MS Peak Lists}
\usage{
\S4method{peakLists}{MSPeakLists}(obj)

\S4method{averagedPeakLists}{MSPeakLists}(obj)

\S4method{analyses}{MSPeakLists}(obj)

\S4method{groupNames}{MSPeakLists}(obj)

\S4method{length}{MSPeakLists}(x)

\S4method{show}{MSPeakLists}(object)

\S4method{[}{MSPeakLists,ANY,ANY,missing}(x, i, j, ..., reAverage = FALSE, drop = TRUE)

\S4method{[[}{MSPeakLists,ANY,ANY}(x, i, j)

\S4method{$}{MSPeakLists}(x, name)

\S4method{as.data.table}{MSPeakLists}(x, fGroups = NULL, averaged = TRUE)

\S4method{filter}{MSPeakLists}(
  obj,
  absMSIntThr = NULL,
  absMSMSIntThr = NULL,
  relMSIntThr = NULL,
  relMSMSIntThr = NULL,
  topMSPeaks = NULL,
  topMSMSPeaks = NULL,
  minMSMSPeaks = NULL,
  isolatePrec = NULL,
  deIsotopeMS = FALSE,
  deIsotopeMSMS = FALSE,
  withMSMS = FALSE,
  annotatedBy = NULL,
  retainPrecursorMSMS = TRUE,
  reAverage = FALSE,
  negate = FALSE
)

\S4method{plotSpectrum}{MSPeakLists}(
  obj,
  groupName,
  analysis = NULL,
  MSLevel = 1,
  title = NULL,
  specSimParams = getDefSpecSimParams(),
  xlim = NULL,
  ylim = NULL,
  ...
)

\S4method{spectrumSimilarity}{MSPeakLists}(
  obj,
  groupName1,
  groupName2 = NULL,
  analysis1 = NULL,
  analysis2 = NULL,
  MSLevel = 1,
  specSimParams = getDefSpecSimParams(),
  NAToZero = FALSE,
  drop = TRUE
)

\S4method{analysisInfo}{MSPeakListsSet}(obj)

\S4method{show}{MSPeakListsSet}(object)

\S4method{[}{MSPeakListsSet,ANY,ANY,missing}(x, i, j, ..., reAverage = FALSE, sets = NULL, drop = TRUE)

\S4method{as.data.table}{MSPeakListsSet}(x, fGroups = NULL, averaged = TRUE)

\S4method{filter}{MSPeakListsSet}(
  obj,
  ...,
  annotatedBy = NULL,
  retainPrecursorMSMS = TRUE,
  reAverage = FALSE,
  negate = FALSE,
  sets = NULL
)

\S4method{plotSpectrum}{MSPeakListsSet}(
  obj,
  groupName,
  analysis = NULL,
  MSLevel = 1,
  title = NULL,
  specSimParams = getDefSpecSimParams(),
  xlim = NULL,
  ylim = NULL,
  perSet = TRUE,
  mirror = TRUE,
  ...
)

\S4method{spectrumSimilarity}{MSPeakListsSet}(
  obj,
  groupName1,
  groupName2 = NULL,
  analysis1 = NULL,
  analysis2 = NULL,
  MSLevel = 1,
  specSimParams = getDefSpecSimParams(),
  NAToZero = FALSE,
  drop = TRUE
)

\S4method{unset}{MSPeakListsSet}(obj, set)

getDefIsolatePrecParams(...)
}
\arguments{
\item{obj, x, object}{The \code{\link{MSPeakLists}} object to access.}

\item{i, j}{For \code{[}/\code{[[}: A numeric or character value which is used to select analyses/feature groups by
their index or name, respectively (for the order/names see \code{analyses()/groupNames()}).\cr\cr For \code{[}: Can also be logical to perform logical selection
(similar to regular vectors). If missing all analyses/feature groups are selected.\cr\cr For \code{[[}: should be a scalar value. If \code{j} is not specified, \code{i} selects by feature groups instead.}

\item{\dots}{Further arguments passed to \code{\link[graphics]{plot}}.

  \setsPassedArgs1{MSPeakLists}}

\item{reAverage}{Set to \code{TRUE} to regenerate group averaged MS peak lists. \strong{NOTE} it is very important
that any annotation data relying on MS peak lists (formulae/compounds) are regenerated afterwards! Otherwise it is
likely that \emph{e.g.} plotting methods will use wrong MS/MS data.}

\item{drop}{If set to \code{TRUE} and if the comparison is made between two spectra then \code{\link{drop}} is used
to reduce the \code{matrix} return value to a \code{numeric} vector.}

\item{name}{The feature group name (partially matched).}

\item{fGroups}{The \code{\link{featureGroups}} object that was used to
generate this object. If not \code{NULL} it is used to add feature group
information (retention and \emph{m/z} values).}

\item{averaged}{If \code{TRUE} then feature group averaged peak list data is
used.}

\item{absMSIntThr, absMSMSIntThr, relMSIntThr, relMSMSIntThr}{Absolute/relative intensity threshold for MS or MS/MS peak
lists. \code{NULL} for none.}

\item{topMSPeaks, topMSMSPeaks}{Only consider this amount of MS or MS/MS peaks with highest intensity. \code{NULL} to
consider all.}

\item{minMSMSPeaks}{If the number of peaks in an MS/MS peak list (\strong{excluding} the precursor peak) is lower
than this it will be completely removed. Set to \code{NULL} to ignore.}

\item{isolatePrec}{If not \code{NULL} then value should be a \code{list} with parameters used for isolating the
precursor and its isotopes in MS peak lists (see \verb{Isolating precursor data}). Alternatively, \code{TRUE} to
apply the filter with default settings (as given with \code{getDefIsolatePrecParams}).}

\item{deIsotopeMS, deIsotopeMSMS}{Remove any isotopic peaks in MS or MS/MS peak lists. This may improve data
processing steps which do not assume the presence of isotopic peaks (e.g. MetFrag for MS/MS). Note that
\code{getMzRPeakLists} does not (yet) support flagging of isotopes.}

\item{withMSMS}{If set to \code{TRUE} then only results will be retained for which MS/MS data is available. if
\code{negate=TRUE} then only results \emph{without} MS/MS data will be retained.}

\item{annotatedBy}{Either a \code{\link{formulas}} or \code{\link{compounds}} object, or a \code{list} with both. Any
MS/MS peaks that are \emph{not} annotated by any of the candidates in the specified objects are removed.}

\item{retainPrecursorMSMS}{If \code{TRUE} then precursor peaks will never be filtered out from MS/MS peak lists (note
that precursors are never removed from MS peak lists). The \code{negate} argument does not affect this setting.}

\item{negate}{If \code{TRUE} then filters are applied in opposite manner.}

\item{groupName}{The name of the feature group for which a plot should be made. To compare spectra, two group names
can be specified.}

\item{analysis}{The name of the analysis for which a plot should be made. If \code{NULL} then data from the feature
group averaged peak list is used. When comparing spectra, either \code{NULL} or the analyses for both spectra
should be specified.}

\item{MSLevel}{The MS level: \samp{1} for regular MS, \samp{2} for MSMS.}

\item{title}{The title of the plot. If \code{NULL} a title will be automatically made.}

\item{specSimParams}{A named \code{list} with parameters that influence the calculation of MS spectra similarities.
See the \link[=specSimParams]{spectral similarity parameters} documentation for more details.}

\item{xlim, ylim}{Sets the plot size limits used by
\code{\link[graphics]{plot}}. Set to \code{NULL} for automatic plot sizing.}

\item{groupName1, groupName2}{The names of the feature groups for which the comparison should be made. If both
arguments are specified then a comparison is made with the spectra specified by \code{groupName1} \emph{vs} those
specified by \code{groupName2}. The length of either can be \samp{>1} to generate a comparison matrix.
Alternatively, if \code{groupName2} is \code{NULL} then all the spectra specified in \code{groupName1} will be
compared with eachother, \emph{i.e.} resulting in a square similarity matrix.}

\item{analysis1, analysis2}{The name of the analysis (analyses) for the comparison. If \code{NULL} then data from the
feature group averaged peak list is used. Otherwise, should be the same length as
\code{groupName1}/\code{groupName2}.}

\item{NAToZero}{Set to \code{TRUE} to convert \code{NA} similarities (\emph{i.e.} when no similarity could be
calculated) to zero values.}

\item{sets}{\setsWF A \code{character} with name(s) of the sets to keep (or remove if \code{negate=TRUE}).}

\item{perSet, mirror}{\setsWF If \code{perSet=TRUE} then the set specific mass peaks are annotated separately.
Furthermore, if \code{mirror=TRUE} (and there are two sets in the object) then a mirror plot is generated.}

\item{set}{\setsWF The name of the set.}
}
\value{
\code{peakLists} returns a nested list containing MS (and MS/MS where
  available) peak lists per feature group and per analysis. The format is:
  \code{[[analysis]][[featureGroupName]][[MSType]][[PeakLists]]} where
  \code{MSType} is either \code{"MS"} or \code{"MSMS"} and \code{PeakLists} a
  \code{\link{data.table}} containing all \emph{m/z} values (\code{mz}
  column) and their intensities (\code{intensity} column). In addition, the
  peak list tables may contain a \code{cmp} column which contains an unique
  alphabetical identifier to which isotopic cluster (or "compound") a mass
  belongs (only supported by MS peak lists generated by Bruker tools at the
  moment).

\code{averagedPeakLists} returns a nested list of feature group
  averaged peak lists in a similar format as \code{peakLists}.
}
\description{
Contains all MS (and MS/MS where available) peak lists for a \code{\link{featureGroups}} object.
}
\details{
Objects for this class are returned by \code{\link{generateMSPeakLists}}.

The \code{getDefIsolatePrecParams} is used to create a parameter
  list for isolating the precursor and its isotopes (see \verb{Isolating precursor data}).
}
\section{Methods (by generic)}{
\itemize{
\item \code{peakLists}: Accessor method to obtain the MS peak lists.

\item \code{averagedPeakLists}: Accessor method to obtain the feature group averaged
MS peak lists.

\item \code{analyses}: returns a \code{character} vector with the names of the
analyses for which data is present in this object.

\item \code{groupNames}: returns a \code{character} vector with the names of the
feature groups for which data is present in this object.

\item \code{length}: Obtain total number of \emph{m/z} values.

\item \code{show}: Shows summary information for this object.

\item \code{[}: Subset on analyses/feature groups.

\item \code{[[}: Extract a list with MS and MS/MS (if available) peak
lists. If the second argument (\code{j}) is not specified the averaged peak
lists for the group specified by the first argument (\code{i}) will be
returned.

\item \code{$}: Extract group averaged MS peaklists for a feature group.

\item \code{as.data.table}: Returns all MS peak list data in a table.

\item \code{filter}: provides post filtering of generated MS peak lists, which may further enhance quality of
subsequent workflow steps (\emph{e.g.} formulae calculation and compounds identification) and/or speed up these
processes. The filters are applied to all peak lists for each analysis. These peak lists are subsequently averaged
to update group averaged peak lists. However, since version \samp{1.1}, the resulting feature group lists are
\emph{not} filtered afterwards.

\item \code{plotSpectrum}: Plots a spectrum using MS or MS/MS peak lists for a given feature group. Two spectra can be
compared when two feature groups are specified.

\item \code{spectrumSimilarity}: Calculates the spectral similarity between two or more spectra.
}}

\section{Slots}{

\describe{
\item{\code{peakLists}}{Contains a list of all MS (and MS/MS) peak lists. Use the \code{peakLists} method for access.}

\item{\code{metadata}}{Metadata for all spectra used to generate peak lists. Follows the format of the \code{peakLists} slot.}

\item{\code{averagedPeakLists}}{A \code{list} with averaged MS (and MS/MS) peak lists for each feature group.}

\item{\code{avgPeakListArgs}}{A \code{list} with arguments used to generate feature group averaged MS(/MS) peak lists.}

\item{\code{origFGNames}}{A \code{character} with the original input feature group names.}

\item{\code{analysisInfo}}{\setsWF  \link[=analysis-information]{Analysis information}. Use the \code{analysisInfo} method
for access.}
}}

\section{Isolating precursor data}{
 Formula calculation typically relies on evaluating the measured isotopic pattern
  from the precursor to score candidates. Some algorithms (currently only \command{GenForm}) penalize candidates if
  mass peaks are present in MS1 spectra that do not contribute to the isotopic pattern. Since these spectra are
  typically very 'noisy' due to background and co-eluting ions, an additional filtering step may be recommended prior
  to formula calculation. During this precursor isolation step all mass peaks are removed that are (1) not the
  precursor and (2) not likely to be an isotopologue of the precursor. To determine potential isotopic peaks the
  following parameters are used:

  \itemize{

  \item \code{maxIsotopes} The maximum number of isotopes to consider. For instance, a value of \samp{5} means that
  \code{M+0} (\emph{i.e.} the monoisotopic peak) till \code{M+5} is considered. All mass peaks outside this range are
  removed.

  \item \code{mzDefectRange} A two-sized \code{vector} specifying the minimum (can be negative) and maximum
  \emph{m/z} defect deviation compared to the precursor \emph{m/z} defect. When chlorinated, brominated or other
  compounds with strong \emph{m/z} defect in their isotopologues are to be considered a higher range may be desired.
  On the other hand, for natural compounds this range may be tightened. Note that the search range is propegated with
  increasing distance from the precursor, \emph{e.g.} the search range is doubled for \code{M+2}, tripled for
  \code{M+3} etc.

  \item \code{intRange} A two-sized \code{vector} specifying the minimum and maximum relative intensity range
  compared to the precursor. For instance, \code{c(0.001, 2)} removes all peaks that have an intensity below 0.1\% or
  above 200\% of that of the precursor.

  \item \code{z} The \code{z} value (\emph{i.e.} absolute charge) to be considerd. For instance, a value of \code{2}
  would look for \code{M+0.5}, \code{M+1} etc. Note that the \code{mzDefectRange} is adjusted accordingly
  (\emph{e.g.} halved if \code{z=2}).

  \item \code{maxGap} The maximum number of missing adjacent isotopic peaks ('gaps'). If the (rounded) \emph{m/z}
  difference to the previous peak exceeds this value then this and all next peaks will be removed. Similar to
  \code{z}, the maximum gap is automatically adjusted for \code{charge}.

  }

  These parameters should be in a \code{list} that is passed to the \code{isolatePrec} argument to \code{filter}. The
  default values can be obtained with the \code{getDefIsolatePrecParams} function:

\code{maxIsotopes=5}; \code{mzDefectRange=c(-0.01, 0.01)}; \code{intRange=c(0.001, 2)}; \code{z=1}; \code{maxGap=2}
}

\section{S4 class hierarchy}{
 \itemize{   \item{\code{\link{workflowStep}}}   \itemize{     \item{\strong{\code{\link{MSPeakLists}}}}     \itemize{       \item{\code{\link{MSPeakListsSet}}}       \item{\code{\link{MSPeakListsUnset}}}     }   } }
}

\section{Source}{
 \code{spectrumSimilarity}: The principles of spectral binning and cosine similarity calculations
  were loosely was based on the code from \code{SpectrumSimilarity()} function of \pkg{OrgMassSpecR}.
}

\section{Sets workflows}{
 \setsWFClass{MSPeakListsSet}{MSPeakLists}

  \setsWFNewMethodsSOExtra{MSPeakListsUnset}{Only the MS peaks that are present in the specified set are kept.}{

  \item \code{analysisInfo} Returns the  \link[=analysis-information]{analysis info} for this object.

  }

  \setsWFChangedMethods{

  \item \code{filter} and the subset operator (\code{[}) Can be used to select data that is only present for selected
  sets. The \code{filter} method is applied for each set individually, and afterwards the results are combined again
  (see \code{\link{generateMSPeakLists}}). Note that this has important implications for \emph{e.g.} relative
  intensity filters (\code{relMSIntThr}/\code{relMSMSIntThr}), \code{topMSPeaks}/\code{topMSMSPeaks} and
  \code{minMSMSPeaks}. Similarly, when the \code{annotatedBy} filter is applied, each set specific MS peak list is
  filtered by the annotation results from only that set.

  \item \code{plotSpectrum} Is able to highlight set specific mass peaks (\code{perSet} and \code{mirror} arguments).

  \item \code{spectrumSimilarity} First calculates similarities for each spectral pair per set (\emph{e.g.} all
  positive mode spectra are compared and then all negative mode spectra are compared). This data is then combined
  into an overall similarity value. How this combination is performed depends on the \code{setCombineMethod} field of
  the \code{\link{specSimParams}} argument.

  }
}

\author{
For \code{spectrumSimilarity}: major contributions by Bas van de Velde for spectral binning and similarity
  calculation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R, R/feature_groups-screening.R,
%   R/feature_groups-screening-set.R, R/feature_groups-tasq.R, R/features-tasq.R
\docType{class}
\name{suspect-screening}
\alias{suspect-screening}
\alias{screenSuspects,featureGroups-method}
\alias{screenSuspects}
\alias{screenSuspects,featureGroupsScreening-method}
\alias{numericIDLevel}
\alias{genIDLevelRulesFile}
\alias{screenSuspects,featureGroupsSet-method}
\alias{screenSuspects,featureGroupsScreeningSet-method}
\alias{featureGroupsBrukerTASQ-class}
\alias{featureGroupsBrukerTASQ}
\alias{featuresBrukerTASQ-class}
\alias{featuresBrukerTASQ}
\title{Target and suspect screening}
\usage{
\S4method{screenSuspects}{featureGroups}(
  fGroups,
  suspects,
  rtWindow,
  mzWindow,
  adduct,
  skipInvalid,
  onlyHits
)

\S4method{screenSuspects}{featureGroupsScreening}(
  fGroups,
  suspects,
  rtWindow,
  mzWindow,
  adduct,
  skipInvalid,
  onlyHits,
  amend = FALSE
)

numericIDLevel(level)

genIDLevelRulesFile(out, inLevels = NULL, exLevels = NULL)

\S4method{screenSuspects}{featureGroupsSet}(
  fGroups,
  suspects,
  rtWindow,
  mzWindow,
  adduct,
  skipInvalid,
  onlyHits
)

\S4method{screenSuspects}{featureGroupsScreeningSet}(
  fGroups,
  suspects,
  rtWindow,
  mzWindow,
  adduct,
  skipInvalid,
  onlyHits,
  amend = FALSE
)
}
\arguments{
\item{fGroups}{The \code{\link{featureGroups}} object that should be screened.}

\item{suspects}{A \code{data.frame} with suspect information. See the \verb{Suspect list format} section below.

  \setsWF Can also be a \code{list} with suspect lists to be used for each set (otherwise the same suspect lists is
  used for all sets). The \code{list} can be named with the sets names to mark which suspect list is to be used with
  which set (\emph{e.g.} \code{suspects=list(positive=suspsPos, negative=suspsNeg)}).}

\item{rtWindow, mzWindow}{The retention time window (in seconds) and \emph{m/z} window that will be used for matching
a suspect (+/- feature data).}

\item{adduct}{An \code{\link{adduct}} object (or something that can be converted to it with \code{\link{as.adduct}}).
Examples: \code{"[M-H]-"}, \code{"[M+Na]+"}. May be \code{NULL}, see \verb{Suspect list format} and \verb{Matching
of suspect masses} sections below.}

\item{skipInvalid}{If set to \code{TRUE} then suspects with invalid data (\emph{e.g.} missing names or other missing
data) will be ignored with a warning. Similarly, any suspects for which mass calculation failed (when no \code{mz}
column is present in the suspect list), for instance, due to invalid \code{SMILES}, will be ignored with a warning.}

\item{onlyHits}{If \code{TRUE} then all feature groups not matched by any of the suspects will be removed.}

\item{amend}{If \code{TRUE} then screening results will be \emph{amended} to the original object.}

\item{level}{The identification level to be converted.}

\item{out}{The file path to the target file.}

\item{inLevels, exLevels}{A \link[=regex]{regular expression} for the
identification levels to include or exclude, respectively. For instance,
\code{exLevels="4|5"} would exclude level 4 and 5 from the output file. Set
to \code{NULL} to ignore.}
}
\value{
\code{screenSuspects} returns a \code{\link{featureGroupsScreening}} object, which is a copy of the input
  \code{fGroups} object amended with additional screening information.
}
\description{
Utilities to screen for analytes with known or suspected identity.
}
\details{
Besides 'full non-target analysis', where compounds may be identified with little to no prior knowledge, a common
strategy is to screen for compounds with known or suspected identity. This may be a generally favorable approach if
possible, as it can significantly reduce the load on data interpretation.

\code{screenSuspects} is used to perform suspect screening. The input \code{\link{featureGroups}} object
  will be screened for suspects by \emph{m/z} values and optionally retention times. Afterwards, any feature groups
  not matched may be kept or removed, depending whether a full non-target analysis is desired.

\code{numericIDLevel} Extracts the numeric part of a given
  identification level (\emph{e.g.} \code{"3a"} becomes \samp{3}).

\code{genIDLevelRulesFile} Generates a template YAML file that is
  used to configure the rules for automatic estimation of identification
  levels. This file can then be used as input for
  \code{\link{annotateSuspects}}.
}
\note{
Both \code{screenSuspects} may use the suspect names to base file names used for reporting, logging etc.
  Therefore, it is important that these are file-compatible names. For this purpose, \code{screenSuspects} will
  automatically try to convert long, non-unique and/or otherwise incompatible suspect names.

For \code{screenSuspects} in some cases you may need to install
  \href{http://openbabel.org/wiki/Main_Page}{OpenBabel} (\emph{e.g.} when only InChI data is available for mass
  calculation).
}
\section{Sets workflows}{
 In a \link[=sets-workflow]{sets workflow}, \code{screenSuspects} performs suspect screening
  for each set separately, and the screening results are combined afterwards. The \code{sets} column in the
  \code{screenInfo} data marks in which sets the suspect hit was found.
}

\section{Suspect list format}{
 the \code{suspects} argument for \code{screenSuspects} should be a \code{data.frame}
  with the following mandatory and optional columns:

  \itemize{

  \item \code{name} The suspect name. Must be file-compatible. (\strong{mandatory})

  \item \code{rt} The retention time (in seconds) for the suspect. If specified the suspect will only be matched if
  its retention matches the experimental value (tolerance defined by the \code{rtWindow} argument).
  (\strong{optional})

  \item \code{neutralMass},\code{formula},\code{SMILES},\code{InChI} The neutral monoisotopic mass, chemical formula,
  SMILES or InChI for the suspect. (data from one of these columns are \strong{mandatory} in case no value from the
  \code{mz} column is available for a suspect)

  \item \code{mz} The ionized \emph{m/z} of the suspect. (\strong{mandatory} unless it can be calculated from one of
  the aforementioned columns)

  \item \code{adduct} A \code{character} that can be converted with \code{\link{as.adduct}}. Can be used to
  automatically calculate values for the \code{mz} column. (\strong{mandatory} unless data from the \code{mz} column
  is available, the \code{adduct} argument is set or \code{fGroups} has adduct annotations)

  \item \code{fragments_mz},\code{fragments_formula} One or more MS/MS fragments (specified as \emph{m/z} or
  formulae, respectively). Multiple values can be specified by separating them with a semicolon (\verb{;}). This data
  is used by \code{\link{annotateSuspects}} to report detected MS/MS fragments and calculate identification levels.
  (\strong{optional})

  }
}

\section{Matching of suspect masses}{
 How the mass of a suspect is matched with the mass of a feature depends on the
  available data: \itemize{

  \item If the suspect has data from the \code{mz} column of the suspect list, then this data is matched with the
  detected feature \emph{m/z}.

  \item Otherwise, if the suspect has data in the \code{adduct} column of the suspect list, this data is used to
  calculate its \code{mz} value, which is then used like above.

  \item In the last case, the neutral mass of the suspect is matched with the neutral mass of the feature. Hence,
  either the \code{adduct} argument needs to be specified, or the \code{featureGroups} input object must have adduct
  annotations.

  }
}

\references{
\insertRef{OBoyle2011}{patRoon}
}
\seealso{
\code{featureGroupsScreening}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{patRoon-defunct}
\alias{patRoon-defunct}
\alias{compoundViewer,featureGroups,MSPeakLists,compounds-method}
\alias{compoundViewer}
\alias{groupFeaturesScreening}
\alias{checkChromatograms}
\title{Defunct functions.}
\usage{
\S4method{compoundViewer}{featureGroups,MSPeakLists,compounds}(fGroups, MSPeakLists, compounds)

groupFeaturesScreening(...)

checkChromatograms(...)
}
\arguments{
\item{MSPeakLists}{A \code{\link{MSPeakLists}} object.}

\item{compounds}{A \code{\link{compounds}} object.}

\item{...}{Ignore.}
}
\description{
These functions do not work anymore and may be updated in the future.
future.
}
\details{
The \code{compoundViewer} method is used to view compound
  identification results. It will display available candidate information
  such as scorings and identifiers, MS/MS spectra with explained peaks and
  chemical structures.

\code{groupFeaturesScreening} was replaced by
  \code{\link{screenSuspects}}. Please refer to this function and the
  handbook for the updated interface.

\code{checkChromatograms} was replaced by
  \code{\link{checkFeatures}}. Please refer to this function and the
  handbook for the updated interface.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/project-tool.R
\name{newProject}
\alias{newProject}
\title{Easily create new \pkg{patRoon} projects}
\usage{
newProject(destPath = NULL)
}
\arguments{
\item{destPath}{Set destination path value to this value (useful for debugging). Set to \code{NULL} for a default
value.}
}
\description{
The \code{newProject} function is used to quickly generate a processing R script. This tool allows the user to
quickly select the targeted analyses, workflow steps and configuring some of their common parameters. This function
requires to be run within a \href{https://www.rstudio.com/}{RStudio} session. The resulting script is either added to
the current open file or to a new file. The \link[=analysis-information]{analysis information} will be written to a
\file{.csv} file so that it can easily be modified afterwards.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R, R/check_components.R,
%   R/check_features.R
\name{check-GUI}
\alias{check-GUI}
\alias{checkComponents,components-method}
\alias{checkComponents}
\alias{importCheckFeaturesSession}
\alias{checkFeatures,featureGroups-method}
\alias{checkFeatures}
\alias{getMCTrainData}
\alias{predictCheckFeaturesSession}
\title{Interactive GUI utilities to check workflow data}
\usage{
\S4method{checkComponents}{components}(
  components,
  fGroups,
  session = "checked-components.yml",
  rtWindow = 30,
  clearSession = FALSE
)

importCheckFeaturesSession(
  sessionIn,
  sessionOut,
  fGroups,
  rtWindow = 6,
  mzWindow = 0.002,
  overWrite = FALSE
)

\S4method{checkFeatures}{featureGroups}(
  fGroups,
  session = "checked-features.yml",
  rtWindow = 30,
  clearSession = FALSE
)

getMCTrainData(fGroups, session)

predictCheckFeaturesSession(fGroups, session, model = NULL, overWrite = FALSE)
}
\arguments{
\item{components}{The \code{\link{components}} to be checked.}

\item{fGroups}{A \code{\link{featureGroups}} object.

  This should be the 'new' object for \code{importCheckFeaturesSession} for which the session needs to be imported.}

\item{session}{The session file name.}

\item{rtWindow}{For \code{checkFeatures} and \code{checkComponents}: the retention time (in seconds) that will be
  subtracted/added to respectively the minimum and maximum retention time of the plotted feature groups. Thus,
  setting this value to a positive value will 'zoom out' on the retention time axis.

  For \code{importCheckFeaturesSession}: the retention time window (seconds) used to relate 'old' with 'new' feature
  groups.}

\item{clearSession}{If \code{TRUE} the session will be completely cleared before starting the GUI. This effectively
removes all selections for data removal.}

\item{sessionIn, sessionOut}{The file names for the input and output sessions.}

\item{mzWindow}{The \emph{m/z} window (in Da) used to relate 'old' with 'new' feature groups.}

\item{overWrite}{Set to \code{TRUE} to overwrite the output session file if it already exists. If \code{FALSE}, the
function will stop with an error message.}

\item{model}{The model that was created with \pkg{MetaClean} and that should be used to predict pass/fail data. If
\code{NULL}, the example model of the \pkg{MetaCleanData} package is used.}
}
\description{
These functions provide interactive utilities to explore and review workflow data using a \pkg{\link{shiny}}
graphical user interface (GUI). In addition, unsatisfactory data (\emph{e.g.} noise identified as a feature and
unrelated feature groups in a component) can easily be selected for removal.
}
\details{
The data selected for removal is stored in \emph{sessions}. These are \file{YAML} files to allow easy external
manipulation. The sessions can be used to restore the selections that were made for data removal when the GUI tool is
executed again. Furthermore, functionality is provided to import and export sessions. To actually remove the data the
\code{\link{filter}} method should be used with the session file as input.

\code{checkComponents} is used to review components and their feature groups contained within. A typical use
  case is to verify that peaks from features that were annotated as related adducts and/or isotopes are correctly
  aligned.

\code{importCheckFeaturesSession} is used to import a session file that was generated from a different
  \code{\link{featureGroups}} object. This is useful to avoid re-doing manual interpretation of chromatographic peaks
  when, for instance, feature group data is re-created with different parameters.

\code{checkFeatures} is used to review chromatographic information for feature groups. Its main purpose is
  to assist in reviewing the quality of detected feature (groups) and easily select unwanted data such as features
  with poor peak shapes or noise.

\code{getMCTrainData} converts a session created by \code{checkFeatures} to a \code{data.frame} that can be
  used by the \pkg{MetaClean} to train a new model. The output format is comparable to that from
  \code{\link[MetaClean]{getPeakQualityMetrics}}.

\code{predictCheckFeaturesSession} Uses ML data from \pkg{MetaClean} to predict the quality (Pass/Fail) of
  feature group data, and converts this to a session which can be reviewed with \code{checkFeatures} and used to
  remove unwanted feature groups by \code{\link[=filter,featureGroups-method]{filter}}.
}
\note{
\code{checkComponents}: Some componentization algorithms (\emph{e.g.} \code{\link{generateComponentsNontarget}}
  and \code{\link{generateComponentsTPs}}) may output components where the same feature group in a component is
  present multiple times, for instance, when multiple TPs are matched to the same feature group. If such a feature
  group is selected for removal, then \emph{all} of its result in the component will be marked for removal.

\code{getMCTrainData} only uses session data for selected feature groups. Selected features for removal are
  ignored, as this is not supported by \pkg{MetaClean}.
}
\references{
\insertRef{Chetnik2020}{patRoon}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups-xcms3.R
\name{groupFeaturesXCMS3}
\alias{groupFeaturesXCMS3}
\alias{groupFeaturesXCMS3,features-method}
\alias{groupFeaturesXCMS3,featuresSet-method}
\title{Group features using XCMS (new interface)}
\usage{
\S4method{groupFeaturesXCMS3}{features}(
  feat,
  rtalign = TRUE,
  loadRawData = TRUE,
  groupParam = xcms::PeakDensityParam(sampleGroups = analysisInfo(feat)$group),
  preGroupParam = groupParam,
  retAlignParam = xcms::ObiwarpParam(),
  verbose = TRUE
)

\S4method{groupFeaturesXCMS3}{featuresSet}(
  feat,
  groupParam = xcms::PeakDensityParam(sampleGroups = analysisInfo(feat)$group),
  verbose = TRUE
)
}
\arguments{
\item{feat}{The \code{\link{features}} object with the features to be grouped.}

\item{rtalign}{Set to \code{TRUE} to enable retention time alignment.}

\item{loadRawData}{Set to \code{TRUE} if analyses are available as \code{mzXML} or \code{mzML} files. Otherwise MS
data is not loaded, and some dummy data (\emph{e.g.} file paths) is used in the returned object.}

\item{groupParam, retAlignParam}{parameter object that is directly passed to
\code{\link[xcms:groupChromPeaks]{xcms::groupChromPeaks}} and \code{\link[xcms:adjustRtime]{xcms::adjustRtime}},
respectively.}

\item{preGroupParam}{grouping parameters applied when features are grouped \emph{prior} to alignment (only with peak
groups alignment).}

\item{verbose}{if \code{FALSE} then no text output will be shown.}
}
\value{
An object of a class which is derived from \code{\link{featureGroups}}.

The \code{featuresSet} method (for \link[=sets-workflow]{sets workflows}) returns a
  \code{\link{featureGroupsSet}} object.
}
\description{
Uses the new \code{xcms3} interface from the \pkg{xcms} package to find features.
}
\details{
This function uses XCMS3 to group features. This function is called when calling \code{groupFeatures} with
  \code{algorithm="xcms3"}.

Grouping of features and alignment of their retention times are performed with the
  \code{\link[xcms:groupChromPeaks]{xcms::groupChromPeaks}} and \code{\link[xcms:adjustRtime]{xcms::adjustRtime}}
  functions, respectively. Both of these functions support an extensive amount of parameters that modify their
  behavior and may therefore require optimization.
}
\section{Sets workflows}{
 \code{loadRawData} and arguments related to retention time alignment are currently not
  supported for \link[=sets-workflow]{sets workflows}.
}

\references{
\addCitations{xcms}{1} \cr\cr \addCitations{xcms}{2} \cr\cr \addCitations{xcms}{3}
}
\seealso{
\code{\link{groupFeatures}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compounds.R, R/compounds-set.R
\docType{class}
\name{compounds-class}
\alias{compounds-class}
\alias{compounds}
\alias{compoundsConsensus-class}
\alias{compoundsConsensus}
\alias{defaultExclNormScores,compounds-method}
\alias{show,compounds-method}
\alias{identifiers,compounds-method}
\alias{identifiers}
\alias{filter,compounds-method}
\alias{addFormulaScoring,compounds-method}
\alias{addFormulaScoring}
\alias{getMCS,compounds-method}
\alias{plotStructure,compounds-method}
\alias{plotScores,compounds-method}
\alias{annotatedPeakList,compounds-method}
\alias{plotSpectrum,compounds-method}
\alias{consensus,compounds-method}
\alias{compoundsSet-class}
\alias{compoundsSet}
\alias{compoundsConsensusSet-class}
\alias{compoundsConsensusSet}
\alias{show,compoundsSet-method}
\alias{delete,compoundsSet-method}
\alias{[,compoundsSet,ANY,missing,missing-method}
\alias{filter,compoundsSet-method}
\alias{plotSpectrum,compoundsSet-method}
\alias{addFormulaScoring,compoundsSet-method}
\alias{annotatedPeakList,compoundsSet-method}
\alias{consensus,compoundsSet-method}
\alias{compoundsUnset-class}
\alias{compoundsUnset}
\alias{unset,compoundsSet-method}
\alias{unset,compoundsConsensusSet-method}
\title{Compound annotations class}
\usage{
\S4method{defaultExclNormScores}{compounds}(obj)

\S4method{show}{compounds}(object)

\S4method{identifiers}{compounds}(compounds)

\S4method{filter}{compounds}(
  obj,
  minExplainedPeaks = NULL,
  minScore = NULL,
  minFragScore = NULL,
  minFormulaScore = NULL,
  scoreLimits = NULL,
  ...
)

\S4method{addFormulaScoring}{compounds}(
  compounds,
  formulas,
  updateScore = FALSE,
  formulaScoreWeight = 1
)

\S4method{getMCS}{compounds}(obj, index, groupName)

\S4method{plotStructure}{compounds}(obj, index, groupName, width = 500, height = 500)

\S4method{plotScores}{compounds}(
  obj,
  index,
  groupName,
  normalizeScores = "max",
  excludeNormScores = defaultExclNormScores(obj),
  onlyUsed = TRUE
)

\S4method{annotatedPeakList}{compounds}(
  obj,
  index,
  groupName,
  MSPeakLists,
  formulas = NULL,
  onlyAnnotated = FALSE
)

\S4method{plotSpectrum}{compounds}(
  obj,
  index,
  groupName,
  MSPeakLists,
  formulas = NULL,
  plotStruct = FALSE,
  title = NULL,
  specSimParams = getDefSpecSimParams(),
  mincex = 0.9,
  xlim = NULL,
  ylim = NULL,
  maxMolSize = c(0.2, 0.4),
  molRes = c(100, 100),
  ...
)

\S4method{consensus}{compounds}(
  obj,
  ...,
  absMinAbundance = NULL,
  relMinAbundance = NULL,
  uniqueFrom = NULL,
  uniqueOuter = FALSE,
  rankWeights = 1,
  labels = NULL
)

\S4method{show}{compoundsSet}(object)

\S4method{delete}{compoundsSet}(obj, i, j, ...)

\S4method{[}{compoundsSet,ANY,missing,missing}(x, i, j, ..., sets = NULL, updateConsensus = FALSE, drop = TRUE)

\S4method{filter}{compoundsSet}(obj, ..., sets = NULL, updateConsensus = FALSE, negate = FALSE)

\S4method{plotSpectrum}{compoundsSet}(
  obj,
  index,
  groupName,
  MSPeakLists,
  formulas = NULL,
  plotStruct = FALSE,
  title = NULL,
  specSimParams = getDefSpecSimParams(),
  mincex = 0.9,
  xlim = NULL,
  ylim = NULL,
  maxMolSize = c(0.2, 0.4),
  molRes = c(100, 100),
  perSet = TRUE,
  mirror = TRUE,
  ...
)

\S4method{addFormulaScoring}{compoundsSet}(
  compounds,
  formulas,
  updateScore = FALSE,
  formulaScoreWeight = 1
)

\S4method{annotatedPeakList}{compoundsSet}(obj, index, groupName, MSPeakLists, formulas = NULL, ...)

\S4method{consensus}{compoundsSet}(
  obj,
  ...,
  absMinAbundance = NULL,
  relMinAbundance = NULL,
  uniqueFrom = NULL,
  uniqueOuter = FALSE,
  rankWeights = 1,
  labels = NULL,
  filterSets = FALSE,
  setThreshold = 0,
  setThresholdAnn = 0
)

\S4method{unset}{compoundsSet}(obj, set)

\S4method{unset}{compoundsConsensusSet}(obj, set)
}
\arguments{
\item{obj, object, compounds, x}{The \code{compound} object.}

\item{minExplainedPeaks, scoreLimits}{Passed to the
\code{\link[=filter,featureAnnotations-method]{featureAnnotations}} method.}

\item{minScore, minFragScore, minFormulaScore}{Minimum overall score, in-silico fragmentation score and formula score,
respectively. Set to \code{NULL} to ignore. The \code{scoreLimits} argument allows for more advanced score
filtering.}

\item{\dots}{For \code{plotSpectrum}: Further arguments passed to \code{\link[graphics]{plot}}.

  For \code{delete}: passed to the function specified as \code{j}.

  for \code{filter}: passed to the \code{\link[=filter,featureAnnotations-method]{featureAnnotations}} method.

  For \code{consensus}: any further (and unique) \code{compounds} objects.

  \setsPassedArgs1{compounds}}

\item{formulas}{The \code{\link{formulas}} object that should be used for scoring/annotation. For \code{plotSpectrum}
and \code{annotatedPeakList}: set to \code{NULL} to ignore.}

\item{updateScore}{If set to \code{TRUE} then the \code{score} column is
updated by adding the normalized \option{formulaScore} (weighted by
\option{formulaScoreWeight}). Currently, this \strong{only} makes sense for
\command{MetFrag} results!}

\item{formulaScoreWeight}{Weight used to update scoring (see
\code{updateScore} parameter).}

\item{index}{The numeric index of the candidate structure.

  For \code{plotStructure} and \code{getMCS}: multiple indices (\emph{i.e.} vector with length >=2) should be
  specified to plot/calculate the most common substructure (MCS). Alternatively, \samp{-1} may be specified to select
  all candidates.

  For \code{plotSpectrum}: two indices can be specified to compare spectra. In this case \code{groupName} should
  specify values for the spectra to compare.}

\item{groupName}{The name of the feature group (or feature groups when comparing spectra) to which the candidate
belongs.}

\item{width, height}{The dimensions (in pixels) of the raster image that
should be plotted.}

\item{normalizeScores}{A \code{character} that specifies how normalization of
annotation scorings occurs. Either
\code{"none"} (no normalization),
\code{"max"} (normalize to max value) or \code{"minmax"} (perform min-max
normalization). Note that normalization of negative scores (e.g. output by
\command{SIRIUS}) is always performed as min-max. Furthermore, currently
normalization for \code{compounds} takes the original min/max scoring
values into account when candidates were generated. Thus, for
\code{compounds} scoring, normalization is not affected when candidate
results were removed after they were generated (\emph{e.g.} by use of
\code{filter}).}

\item{excludeNormScores}{A
  \code{character} vector specifying any compound scoring names that
  should \emph{not} be normalized. Set to \code{NULL} to normalize all
  scorings. Note that whether any normalization occurs is set by the
  \code{excludeNormScores} argument.

  For \code{compounds}: By default \code{score} and
  \code{individualMoNAScore} are set to mimic the behavior of the
  \command{MetFrag} web interface.}

\item{onlyUsed}{If \code{TRUE} then only scorings are plotted that actually
have been used to rank data (see the \code{scoreTypes} argument to
\code{\link{generateCompoundsMetFrag}} for more details).}

\item{MSPeakLists}{The \code{\link{MSPeakLists}} object that was used to generate the candidate}

\item{onlyAnnotated}{Set to \code{TRUE} to filter out any peaks that could
not be annotated.}

\item{plotStruct}{If \code{TRUE} then the candidate structure is drawn in the spectrum. Currently not supported when
comparing spectra.}

\item{title}{The title of the plot. If \code{NULL} a title will be automatically made.}

\item{specSimParams}{A named \code{list} with parameters that influence the calculation of MS spectra similarities.
See the \link[=specSimParams]{spectral similarity parameters} documentation for more details.}

\item{mincex}{The formula annotation labels are automatically scaled. The \code{mincex} argument forces a minimum
\code{cex} value for readability.}

\item{xlim, ylim}{Sets the plot size limits used by
\code{\link[graphics]{plot}}. Set to \code{NULL} for automatic plot sizing.}

\item{maxMolSize}{Numeric vector of size two with the maximum width/height of the candidate structure (relative to
the plot size).}

\item{molRes}{Numeric vector of size two with the resolution of the candidate structure (in pixels).}

\item{absMinAbundance, relMinAbundance}{Minimum absolute or relative
(\samp{0-1}) abundance across objects for a result to be kept. For
instance, \code{relMinAbundance=0.5} means that a result should be present
in at least half of the number of compared objects. Set to \samp{NULL} to
ignore and keep all results. Limits cannot be set when \code{uniqueFrom} is
not \code{NULL}.}

\item{uniqueFrom}{Set this argument to only retain compounds that are unique
within one or more of the objects for which the consensus is made.
Selection is done by setting the value of \code{uniqueFrom} to a
\code{logical} (values are recycled), \code{numeric} (select by index) or a
\code{character} (as obtained with \code{algorithm(obj)}). For
\code{logical} and \code{numeric} values the order corresponds to the order
of the objects given for the consensus. Set to \code{NULL} to ignore.}

\item{uniqueOuter}{If \code{uniqueFrom} is not \code{NULL} and if
\code{uniqueOuter=TRUE}: only retain data that are also unique between
objects specified in \code{uniqueFrom}.}

\item{rankWeights}{A numeric vector with weights of to calculate the mean
ranking score for each candidate. The value will be re-cycled if necessary,
hence, the default value of \samp{1} means equal weights for all considered
objects.}

\item{labels}{A \code{character} with names to use for labelling. If \code{NULL} labels are automatically generated.}

\item{i, j, drop}{Passed to the \code{\link[=[,featureAnnotations,ANY,missing,missing-method]{featureAnnotations}}
method.}

\item{sets}{\setsWF A \code{character} with name(s) of the sets to keep (or remove if \code{negate=TRUE}). Note: if
\code{updateConsensus=FALSE} then the \code{setCoverage} column of the annotation results is not updated.}

\item{updateConsensus}{\setsWF If \code{TRUE} then the annonation consensus among set results is updated. See the
\verb{Sets workflows} section for more details.}

\item{negate}{Passed to the \code{\link[=filter,featureAnnotations-method]{featureAnnotations}} method.}

\item{perSet, mirror}{\setsWF If \code{perSet=TRUE} then the set specific mass peaks are annotated separately.
Furthermore, if \code{mirror=TRUE} (and there are two sets in the object) then a mirror plot is generated.}

\item{filterSets}{\setsWF Controls how algorithms concensus abundance filters are applied. See the \verb{Sets
workflows} section below.}

\item{setThreshold, setThresholdAnn}{\setsWF Thresholds used to create the annotation set consensus. See
\code{\link{generateCompounds}}.}

\item{set}{\setsWF The name of the set.}
}
\value{
\code{addFormulaScoring} returns a \code{compounds} object updated
  with formula scoring.

\code{getMCS} returns an \CRANpkg{rcdk} molecule object
  (\code{IAtomContainer}).

\code{consensus} returns a \code{compounds} object that is produced by merging multiple specified
  \code{compounds} objects.
}
\description{
Contains data for compound annotations for feature groups.
}
\details{
\code{compounds} objects are obtained from \link[=generateCompounds]{compound generators}. This class is derived from
the \code{\link{featureAnnotations}} class, please see its documentation for more methods and other details.
}
\section{Methods (by generic)}{
\itemize{
\item \code{defaultExclNormScores}: Returns default scorings that are excluded from normalization.

\item \code{show}: Show summary information for this object.

\item \code{identifiers}: Returns a list containing for each feature group a
character vector with database identifiers for all candidate compounds. The
list is named by feature group names, and is typically used with the
\code{identifiers} option of \code{\link{generateCompoundsMetFrag}}.

\item \code{filter}: Provides rule based filtering for generated compounds. Useful to eliminate unlikely candidates
and speed up further processing. Also see the \code{\link[=filter,featureAnnotations-method]{featureAnnotations}}
method.

\item \code{addFormulaScoring}: Adds formula ranking data from a \code{\link{formulas}}
object as an extra compound candidate scoring (\code{formulaScore} column).
The formula score for each compound candidate is between \samp{0-1}, where
\emph{zero} means no match with any formula candidates, and \emph{one}
means that the compound candidate's formula is the highest ranked.

\item \code{getMCS}: Calculates the maximum common substructure (MCS)
for two or more candidate structures for a feature group. This method uses
the \code{\link{get.mcs}} function from \CRANpkg{rcdk}.

\item \code{plotStructure}: Plots a structure of a candidate compound using the
\CRANpkg{rcdk} package. If multiple candidates are specified (\emph{i.e.}
by specifying a \code{vector} for \code{index}) then the maximum common
substructure (MCS) of the selected candidates is drawn.

\item \code{plotScores}: Plots a barplot with scoring of a candidate compound.

\item \code{annotatedPeakList}: Returns an MS/MS peak list annotated with data from a
given candidate compound for a feature group.

\item \code{plotSpectrum}: Plots an annotated spectrum for a given candidate compound for a feature group. Two spectra can
be compared by specifying a two-sized vector for the \code{index} and \code{groupName} arguments.

\item \code{consensus}: Generates a consensus of results from multiple
objects. In order to rank the consensus candidates, first
each of the candidates are scored based on their original ranking
(the scores are normalized and the highest ranked candidate gets value
\samp{1}). The (weighted) mean is then calculated for all scorings of each
candidate to derive the final ranking (if an object lacks the candidate its
score will be \samp{0}). The original rankings for each object is stored in
the \code{rank} columns.
}}

\section{Slots}{

\describe{
\item{\code{setThreshold,setThresholdAnn}}{\setsWF A copy of the equally named arguments that were passed when this object
was created by \code{\link{generateCompounds}}.}

\item{\code{origFGNames}}{\setsWF The original (order of) names of the \code{\link{featureGroups}} object that was used to
create this object.}
}}

\note{
The values ranges in the \code{scoreLimits} slot, which are used for normalization of scores, are based on the
  \emph{original} scorings when the compounds were generated (\emph{prior} to employing the \code{topMost} filter to
  \code{\link{generateCompounds}}).
}
\section{S4 class hierarchy}{
 \itemize{   \item{\code{\link{featureAnnotations}}}   \itemize{     \item{\strong{\code{\link{compounds}}}}     \itemize{       \item{\code{\link{compoundsConsensus}}}       \item{\code{\link{compoundsMF}}}       \item{\code{\link{compoundsSet}}}       \itemize{         \item{\code{\link{compoundsConsensusSet}}}       }       \item{\code{\link{compoundsUnset}}}     }   } }
}

\section{Source}{
 Subscripting of formulae for plots generated by
  \code{plotSpectrum} is based on the \code{chemistry2expression} function
  from the \href{https://github.com/schymane/ReSOLUTION}{ReSOLUTION} package.
}

\section{Sets workflows}{
 \setsWFClass{compoundsSet}{compounds}

  \setsWFNewMethodsSO{compoundsUnset}{Only the annotation results that are present in the specified set are kept
  (based on the set consensus, see below for implications).}

  \setsWFChangedMethods{

  \item \code{filter} and the subset operator (\code{[}) Can be used to select data that is only present for selected
  sets. Depending on the \code{updateConsenus}, both either operate on set consensus or original data (see below for
  implications).

  \item \code{annotatedPeakList} Returns a combined annotation table with all sets.

  \item \code{plotSpectrum} Is able to highlight set specific mass peaks (\code{perSet} and \code{mirror} arguments).

  \item \code{consensus} Creates the algorithm consensus based on the original annotation data (see below for
  implications). Then, like the sets workflow method for \code{\link{generateCompounds}}, a consensus is made for all
  sets, which can be controlled with the \code{setThreshold} and \code{setThresholdAnn} arguments. The candidate
  coverage among the different algorithms is calculated for each set (\emph{e.g.} \code{coverage-positive} column)
  and for all sets (\code{coverage} column), which is based on the presence of a candidate in all the algorithms from
  all sets data. The \code{consensus} method for sets workflow data supports the \code{filterSets} argument. This
  controls how the algorithm consensus abundance filters (\code{absMinAbundance}/\code{relMinAbundance}) are applied:
  if \code{filterSets=TRUE} then the minimum of all \code{coverage} set specific columns is used to obtain the
  algorithm abundance. Otherwise the overall \code{coverage} column is used. For instance, consider a consensus
  object to be generated from two objects generated by different algorithms (\emph{e.g.} \command{SIRIUS} and
  \command{MetFrag}), which both have a positive and negative set. Then, if a candidate occurs with both
  algorithms for the positive mode set, but only with the first algorithm in the negative mode set,
  \code{relMinAbundance=1} will remove the candidate if \code{filterSets=TRUE} (because the minimum relative
  algorithm abundance is \samp{0.5}), while \code{filterSets=FALSE} will not remove the candidate (because based on
  all sets data the candidate occurs in both algorithms).

  \item \code{addFormulaScoring} Adds the formula scorings to the original data and re-creates the annotation set consensus (see below for implications).

  }

  Two types of annotation data are stored in a \code{compoundsSet} object: \enumerate{

  \item Annotations that are produced from a consensus between set results (see \code{generateCompounds}).

  \item The 'original' annotation data per set, prior to when the set consensus was made. This includes candidates
  that were filtered out because of the thresholds set by \code{setThreshold} and \code{setThresholdAnn}. However,
  when \code{filter} or subsetting (\code{[}) operations are performed, the original data is also updated.

  }

  In most cases the first data is used. However, in a few cases the original annotation data is used (as indicated
  above), for instance, to re-create the set consensus. It is important to realize that the original annotation data
  may have \emph{additional} candidates, and a newly created set consensus may therefore have 'new' candidates. For
  instance, when the object consists of the sets \code{"positive"} and \code{"negative"} and \code{setThreshold=1}
  was used to create it, then \code{compounds[, sets = "positive", updateConsensus = TRUE]} may now have additional
  candidates, \emph{i.e.} those that were not present in the \code{"negative"} set and were previously removed due to
  the consensus threshold filter.
}

\references{
\addCitations{rcdk}{1}
}
\seealso{
The \code{\link{featureAnnotations}} base class for more relevant methods and
  \code{\link{generateCompounds}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features-bruker.R
\name{findFeaturesBruker}
\alias{findFeaturesBruker}
\title{Find features using Bruker DataAnalysis}
\usage{
findFeaturesBruker(
  analysisInfo,
  doFMF = "auto",
  startRange = 0,
  endRange = 0,
  save = TRUE,
  close = save,
  verbose = TRUE
)
}
\arguments{
\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}

\item{doFMF}{Run the 'Find Molecular Features' algorithm before loading compounds. Valid options are: \code{"auto"}
(run FMF automatically if current results indicate it is necessary) and \code{"force"} (run FMF \emph{always}, even
if cached results exist). Note that checks done if \code{doFMF="auto"} are fairly simplistic, hence set
\code{doFMF="force"} if feature data needs to be updated.}

\item{startRange, endRange}{Start/End retention range (seconds) from which to collect features. A 0 (zero) for
\code{endRange} marks the end of the analysis.}

\item{close, save}{If \code{TRUE} then Bruker files are closed and saved after
processing with DataAnalysis, respectively. Setting \code{close=TRUE}
prevents that many analyses might be opened simultaneously in DataAnalysis,
which otherwise may use excessive memory or become slow. By default
\code{save} is \code{TRUE} when \code{close} is \code{TRUE}, which is
likely what you want as otherwise any processed data is lost.}

\item{verbose}{If set to \code{FALSE} then no text output is shown.}
}
\value{
An object of a class which is derived from \code{\link{features}}.
}
\description{
Uses the 'Find Molecular Features' (FMF) algorithm of Bruker DataAnalysis vendor software to find features.
}
\details{
This function uses Bruker to automatically find features. This function is called when calling \code{findFeatures} with
  \code{algorithm="bruker"}.

The resulting 'compounds' are transferred from DataAnalysis and stored as features.

  This algorithm only works with Bruker data files (\code{.d} extension) and requires Bruker DataAnalysis
  and the \pkg{RDCOMClient} package to be installed. Furthermore, DataAnalysis combines multiple related masses in a
  feature (\emph{e.g.} isotopes, adducts) but does not report the actual (monoisotopic) mass of the feature.
  Therefore, it is simply assumed that the feature mass equals that of the highest intensity mass peak.
}
\note{
If any errors related to \command{DCOM} appear it might be necessary to
  terminate DataAnalysis (note that DataAnalysis might still be running as a
  background process). The \command{ProcessCleaner} application installed
  with DataAnalayis can be used for this.
}
\seealso{
\code{\link{findFeatures}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-exported.R
\name{createCOMReference}
\alias{createCOMReference}
\title{Internal fix for \pkg{RDCOMClient}, ignore.}
\usage{
createCOMReference(ref, className)
}
\arguments{
\item{ref, className}{ignore}
}
\description{
Internal fix for \pkg{RDCOMClient}, ignore.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formulas.R
\name{generateFormulas}
\alias{generateFormulas}
\alias{generateFormulas,featureGroups-method}
\title{Automatic chemical formula generation}
\usage{
\S4method{generateFormulas}{featureGroups}(fGroups, MSPeakLists, algorithm, ...)
}
\arguments{
\item{fGroups}{\code{\link{featureGroups}} object for which formulae should be generated. This should be the same or
a subset of the object that was used to create the specified \code{MSPeakLists}. In the case of a subset only the
remaining feature groups in the subset are considered.}

\item{MSPeakLists}{An \code{\link{MSPeakLists}} object that was generated for the supplied \code{fGroups}.}

\item{algorithm}{A character string describing the algorithm that should be
used: \code{"bruker"}, \code{"genform"}, \code{"sirius"}}

\item{\dots}{Any parameters to be passed to the selected formula generation algorithm.}
}
\value{
A \code{\link{formulas}} object containing all generated formulae.
}
\description{
Automatically calculate chemical formulae for all feature groups.
}
\details{
Several algorithms are provided to automatically generate formulae for given feature groups. All algorithms use the
accurate mass of a feature to back-calculate candidate formulae. Depending on the algorithm and data availability,
other data such as isotopic pattern and MS/MS fragments may be used to further improve formula assignment and
ranking.

\code{generateFormulas} is a generic function that will generateFormulas by one of the supported algorithms. The actual
  functionality is provided by algorithm specific functions such as \code{generateFormulasDA} and \code{generateFormulasGenForm}. While these
  functions may be called directly, \code{generateFormulas} provides a generic interface and is therefore usually preferred.
}
\section{Candidate assignment}{
 Formula candidate assignment occurs in one of the following ways: \itemize{

  \item Candidates are first generated for each feature and then pooled to form consensus candidates for the feature
  group.

  \item Candidates are directly generated for each feature group by group averaged MS peak list data.

  }

  With approach (1), scorings and mass errors are averaged and outliers are removed (controlled by
  \code{featThreshold} and \code{featThresholdAnn} arguments). Other candidate properties that cannot be averaged are
  from the feature from the analysis as specified in the \code{"analysis"} column of the results. The second approach only generates candidate formulae once for every feature group, and is therefore generally much
  faster. However, this inherently prevents removal of outliers.

  Note that with either approach subsequent workflow steps that use formula data (\emph{e.g.}
  \code{\link{addFormulaScoring}} and \link{reporting} functions) only use formula data that was eventually assigned
  to feature groups.
}

\section{Scorings}{
 Each algorithm implements their own scoring system. Their names have been harmonized where
  possible. An overview is obtained with the \code{\link{formulaScorings}} function:
  \Sexpr[results=rd,echo=FALSE,stage=build]{patRoon:::tabularRD(patRoon::formulaScorings())}
}

\section{Sets workflows}{
 With a \link[=sets-workflow]{sets workflow}, annotation is first performed for each set.
  This is important, since the annotation algorithms typically cannot work with data from mixed ionization modes. The
  annotation results are then combined to generate a \emph{sets consensus}: \itemize{

  \item The annotation tables for each feature group from the set specific data are combined. Rows with overlapping
  candidates (determined by the neutral formula) are merged.

  \item Set specific data (\emph{e.g.} the ionic formula) is retained by renaming their columns with set specific
  names.

  \item The MS/MS fragment annotations (\code{fragInfo} column) from each set are combined.

  \item The scorings for each set are averaged to calculate overall scores.

  \item The candidates are re-ranked based on their average ranking among the set data (if a candidate is absent in a
  set it is assigned the poorest rank in that set).

  \item The coverage of each candidate among sets is calculated. Depending on the \code{setThreshold} and
  \code{setThresholdAnn} arguments, candidates with low abundance are removed.

  }
}

\seealso{
The \code{\link{formulas}} output class and its methods and the algorithm specific functions:
  \code{\link{generateFormulasDA}}, \code{\link{generateFormulasGenForm}}, \code{\link{generateFormulasSIRIUS}}

The \href{https://www.researchgate.net/publication/307964728_MOLGEN-MSMS_Software_User_Manual}{GenForm
  manual} (also known as MOLGEN-MSMS).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups.R, R/feature_groups-set.R,
%   R/feature_groups-bruker.R, R/feature_groups-envimass.R,
%   R/feature_groups-kpic2.R, R/feature_groups-openms.R,
%   R/feature_groups-sirius.R, R/feature_groups-xcms.R, R/feature_groups-xcms3.R
\docType{class}
\name{featureGroups-class}
\alias{featureGroups-class}
\alias{featureGroups}
\alias{names,featureGroups-method}
\alias{analyses,featureGroups-method}
\alias{replicateGroups,featureGroups-method}
\alias{groupNames,featureGroups-method}
\alias{length,featureGroups-method}
\alias{show,featureGroups-method}
\alias{groupTable,featureGroups-method}
\alias{groupTable}
\alias{analysisInfo,featureGroups-method}
\alias{groupInfo,featureGroups-method}
\alias{groupInfo}
\alias{featureTable,featureGroups-method}
\alias{getFeatures,featureGroups-method}
\alias{groupFeatIndex,featureGroups-method}
\alias{groupFeatIndex}
\alias{groupQualities,featureGroups-method}
\alias{groupQualities}
\alias{groupScores,featureGroups-method}
\alias{groupScores}
\alias{annotations,featureGroups-method}
\alias{adducts,featureGroups-method}
\alias{adducts<-,featureGroups-method}
\alias{[,featureGroups,ANY,ANY,missing-method}
\alias{[[,featureGroups,ANY,ANY-method}
\alias{$,featureGroups-method}
\alias{delete,featureGroups-method}
\alias{export,featureGroups-method}
\alias{as.data.table,featureGroups-method}
\alias{unique,featureGroups-method}
\alias{overlap,featureGroups-method}
\alias{overlap}
\alias{calculatePeakQualities,featureGroups-method}
\alias{selectIons,featureGroups-method}
\alias{selectIons}
\alias{featureGroupsSet-class}
\alias{featureGroupsSet}
\alias{sets,featureGroupsSet-method}
\alias{adducts,featureGroupsSet-method}
\alias{adducts<-,featureGroupsSet-method}
\alias{delete,featureGroupsSet-method}
\alias{show,featureGroupsSet-method}
\alias{featureTable,featureGroupsSet-method}
\alias{[,featureGroupsSet,ANY,ANY,missing-method}
\alias{export,featureGroupsSet-method}
\alias{as.data.table,featureGroupsSet-method}
\alias{unique,featureGroupsSet-method}
\alias{overlap,featureGroupsSet-method}
\alias{selectIons,featureGroupsSet-method}
\alias{featureGroupsUnset-class}
\alias{featureGroupsUnset}
\alias{unset,featureGroupsSet-method}
\alias{featureGroupsBruker-class}
\alias{featureGroupsBruker}
\alias{featureGroupsEnviMass-class}
\alias{featureGroupsEnviMass}
\alias{featureGroupsKPIC2-class}
\alias{featureGroupsKPIC2}
\alias{delete,featureGroupsKPIC2-method}
\alias{featureGroupsOpenMS-class}
\alias{featureGroupsOpenMS}
\alias{featureGroupsSIRIUS-class}
\alias{featureGroupsSIRIUS}
\alias{featureGroupsXCMS-class}
\alias{featureGroupsXCMS}
\alias{delete,featureGroupsXCMS-method}
\alias{featureGroupsXCMS3-class}
\alias{featureGroupsXCMS3}
\alias{delete,featureGroupsXCMS3-method}
\title{Base class for grouped features.}
\usage{
\S4method{names}{featureGroups}(x)

\S4method{analyses}{featureGroups}(obj)

\S4method{replicateGroups}{featureGroups}(obj)

\S4method{groupNames}{featureGroups}(obj)

\S4method{length}{featureGroups}(x)

\S4method{show}{featureGroups}(object)

\S4method{groupTable}{featureGroups}(object, areas = FALSE)

\S4method{analysisInfo}{featureGroups}(obj)

\S4method{groupInfo}{featureGroups}(fGroups)

\S4method{featureTable}{featureGroups}(obj)

\S4method{getFeatures}{featureGroups}(obj)

\S4method{groupFeatIndex}{featureGroups}(fGroups)

\S4method{groupQualities}{featureGroups}(fGroups)

\S4method{groupScores}{featureGroups}(fGroups)

\S4method{annotations}{featureGroups}(obj)

\S4method{adducts}{featureGroups}(obj)

\S4method{adducts}{featureGroups}(obj) <- value

\S4method{[}{featureGroups,ANY,ANY,missing}(x, i, j, ..., rGroups, results, drop = TRUE)

\S4method{[[}{featureGroups,ANY,ANY}(x, i, j)

\S4method{$}{featureGroups}(x, name)

\S4method{delete}{featureGroups}(obj, i = NULL, j = NULL, ...)

\S4method{export}{featureGroups}(obj, type, out)

\S4method{as.data.table}{featureGroups}(
  x,
  average = FALSE,
  areas = FALSE,
  features = FALSE,
  qualities = FALSE,
  regression = FALSE,
  averageFunc = mean,
  normFunc = NULL,
  FCParams = NULL
)

\S4method{unique}{featureGroups}(x, which, relativeTo = NULL, outer = FALSE)

\S4method{overlap}{featureGroups}(fGroups, which, exclusive)

\S4method{calculatePeakQualities}{featureGroups}(
  obj,
  weights,
  flatnessFactor,
  avgFunc = mean,
  parallel = TRUE
)

\S4method{selectIons}{featureGroups}(
  fGroups,
  components,
  prefAdduct,
  onlyMonoIso = TRUE,
  chargeMismatch = "adduct"
)

\S4method{sets}{featureGroupsSet}(obj)

\S4method{adducts}{featureGroupsSet}(obj, set, ...)

\S4method{adducts}{featureGroupsSet}(obj, set, reGroup = TRUE) <- value

\S4method{delete}{featureGroupsSet}(obj, i = NULL, j = NULL, ...)

\S4method{show}{featureGroupsSet}(object)

\S4method{featureTable}{featureGroupsSet}(obj)

\S4method{[}{featureGroupsSet,ANY,ANY,missing}(x, i, j, ..., rGroups, sets = NULL, drop = TRUE)

\S4method{export}{featureGroupsSet}(obj, type, out, set)

\S4method{as.data.table}{featureGroupsSet}(
  x,
  average = FALSE,
  areas = FALSE,
  features = FALSE,
  qualities = FALSE,
  regression = FALSE,
  averageFunc = mean,
  normFunc = NULL,
  FCParams = NULL
)

\S4method{unique}{featureGroupsSet}(x, which, ..., sets = FALSE)

\S4method{overlap}{featureGroupsSet}(fGroups, which, exclusive, sets = FALSE)

\S4method{selectIons}{featureGroupsSet}(fGroups, components, prefAdduct, ...)

\S4method{unset}{featureGroupsSet}(obj, set)

\S4method{delete}{featureGroupsKPIC2}(obj, ...)

\S4method{delete}{featureGroupsXCMS}(obj, ...)

\S4method{delete}{featureGroupsXCMS3}(obj, ...)
}
\arguments{
\item{areas}{If set to \code{TRUE} then areas are considered instead of peak intensities.

  For \code{as.data.table}: ignored if \code{features=TRUE}, as areas of features are always reported.}

\item{fGroups, obj, x, object}{\code{featureGroups} object to be accessed.}

\item{value}{For \code{adducts<-}: A \code{character} with adduct annotations assigned to each feature group. The
length should equal the number of feature groups. Can be named with feature group names to customize the assignment
order.}

\item{i, j}{For \code{[}/\code{[[}: A numeric or character value which is used to select analyses/feature groups by
their index or name, respectively (for the order/names see \code{analyses()/names()}).\cr\cr For \code{[}: Can also be logical to perform logical selection
(similar to regular vectors). If missing all analyses/feature groups are selected.\cr\cr For \code{[[}: should be a scalar value. If \code{j} is not specified, \code{i} selects by feature groups instead.\cr\cr For \code{delete}: The data to remove from. \code{i} are the
analyses as numeric index, logical or character, \code{j} the feature groups as numeric index, logical or character. If either is
\code{NULL} then data for all is removed. \code{j} may also be a function: it will be called for each 
feature group, with a vector of the group intensities as first argument, the group name as second argument, and any other arguments passed as
\code{\dots} to \code{delete}. The return value of this function specifies the analyses of the features in the group to be removed (same format as \code{i}).}

\item{\dots}{For the \code{"["} operator: ignored.

  For \code{delete}: passed to the function specified as \code{j}.

  \setsPassedArgs1{featureGroups}}

\item{rGroups}{For \code{[}: An optional \code{character} vector: if specified only keep results for the given
replicate groups (equivalent to the \code{rGroups} argument to \code{\link[=filter,featureGroups-method]{filter}}).}

\item{results}{Optional argument. If specified only feature groups with results in the specified object are kept. The
class of \code{results} should be \code{\link{featureAnnotations}} or \code{\link{components}}. Multiple objects
can be specified in a \code{list}: in this case a feature group is kept if it has a result in \emph{any} of the
objects (equivalent to the \code{results} argument to \code{\link[=filter,featureGroups-method]{filter}}).}

\item{drop}{ignored.}

\item{name}{The feature group name (partially matched).}

\item{type}{The export type: \code{"brukerpa"} (Bruker ProfileAnalysis), \code{"brukertasq"} (Bruker TASQ) or
\code{"mzmine"} (MZmine).}

\item{out}{The destination file for the exported data.}

\item{average}{If \code{TRUE} then data within replicate groups are averaged.

  For \code{as.data.table}: if \code{features=TRUE} other feature properties are also averaged.}

\item{features}{If \code{TRUE} then feature specific data will be added. If \code{average=TRUE} this data will be
averaged for each feature group.}

\item{qualities}{Adds feature (group) qualities (\code{qualities="quality"}), scores (\code{qualities="score"}) or
both (\code{qualities="both"}), if this data is available (\emph{i.e.} from \code{calculatePeakQualities}). If
\code{qualities=FALSE} then nothing is reported.}

\item{regression}{Set to \code{TRUE} to add regression data for each feature group. For this a linear model is
created (intensity/area [depending on \code{areas} argument] \emph{vs} concentration). The model concentrations
(e.g. of a set of standards) is derived from the \code{conc} column of the \link[=analysis-information]{analysis
information}. From this model the intercept, slope and R2 is added to the output. In addition, when
\code{features=TRUE}, concentrations for each feature are added. Note that no regression information is added when
no \code{conc} column is present in the analysis information or when less than two concentrations are specified
(\emph{i.e.} the minimum amount).}

\item{averageFunc}{Function used for averaging. Only used when \code{average=TRUE} or \code{FCParams != NULL}.}

\item{normFunc}{Function that should be used for normalization of data. The function is called for all
intensities/areas of a feature group and these quantities are divided by the result of the function call. For
example, when \code{\link{max}} is used normalized intensities will be between zero and one. If all quantities are
zero then the function will not be called. Set to \code{NULL} to perform no normalization.}

\item{FCParams}{A parameter list to calculate Fold change data. See \code{getFCParams} for more details. Set to
\code{NULL} to not perform FC calculations.}

\item{which}{A character vector with replicate groups used for comparison.}

\item{relativeTo}{A character vector with replicate groups that should be
used for unique comparison. If \code{NULL} then all replicate groups are
used for comparison. Replicate groups specified in \code{which} are
ignored.}

\item{outer}{If \code{TRUE} then only feature groups are kept which do not
overlap between the specified replicate groups for the \code{which}
parameter.}

\item{exclusive}{If \code{TRUE} then all feature groups are removed that are
not unique to the given replicate groups.}

\item{weights}{A named \code{numeric} vector that defines the weight for each score to calculate the
\verb{totalScore}. The names of the vector follow the score names. Unspecified weights are defaulted to \samp{1}.
Example: \code{weights=c(ApexBoundaryRatioScore=0.5, GaussianSimilarityScore=2)}.}

\item{flatnessFactor}{Passed to \pkg{MetaClean} as the \code{flatness.factor} argument to
\code{\link[MetaClean]{calculateJaggedness}} and \code{\link[MetaClean]{calculateModality}}.}

\item{avgFunc}{The function used to average the peak qualities and scores for each feature group.}

\item{parallel}{If set to \code{TRUE} then code is executed in parallel through the \CRANpkg{futures} package. Please
see the parallelization section in the handbook for more details.}

\item{components}{The \code{components} object that was generated for the given \code{featureGroups} object.
Obviously, the components must be created with algorithms that support adduct/isotope annotations, such as those
from \pkg{RAMClustR} and \pkg{cliqueMS}.}

\item{prefAdduct}{The 'preferred adduct' (see method description). This is often \code{"[M+H]+"} or \code{"[M-H]-"}.}

\item{onlyMonoIso}{Set to \code{TRUE} to only keep feature groups that were annotated as monoisotopic. Feature groups
are never removed by this setting if no isotope annotations are available.}

\item{chargeMismatch}{Specifies how to deal with a mismatch in charge between adduct and isotope annotations. Valid
values are: \code{"adduct"} (ignore isotope annotation), \code{"isotope"} (ignore adduct annotation), \code{"none"}
(ignore both annotations) and \code{"ignore"} (don't check for charge mismatches). \emph{Important}: when
\command{OpenMS} is used to find features, it already removes any detected non-monoisotopic features by default.
Hence, in such case setting \code{chargeMismatch="adduct"} is more appropriate.}

\item{set}{\setsWF The name of the set.}

\item{reGroup}{\setsWF Set to \code{TRUE} to re-group the features after the adduct annotations are changed. See the
\verb{Sets workflow} section for more details.}

\item{sets}{\setsWF For \code{[}: a \code{character} with name(s) of the sets to keep.

  For \code{overlap} and \code{unique}: If \code{TRUE} then the \code{which} argument changes its meaning and is used
  to specify the names of the sets to be compared.}
}
\value{
\code{delete} returns the object for which the specified data was removed.

\code{calculatePeakQualities} returns a modified object amended with peak qualities and scores.

\code{selectIons} returns a \code{featureGroups} object with only the selected feature groups and amended
  with adduct annotations.
}
\description{
This class holds all the information for grouped features.
}
\details{
The \code{featureGroup} class is the workhorse of \pkg{patRoon}: almost all functionality operate on its instantiated
objects. The class holds all information from grouped features (obtained from \code{\link{features}}). This class
itself is \code{virtual}, hence, objects are not created directly from it. Instead, 'feature groupers' such as
\code{\link{groupFeaturesXCMS}} return a \code{featureGroups} derived object after performing the actual grouping of
features across analyses.
}
\section{Methods (by generic)}{
\itemize{
\item \code{names}: Obtain feature group names.

\item \code{analyses}: returns a \code{character} vector with the names of the
analyses for which data is present in this object.

\item \code{replicateGroups}: returns a \code{character} vector with the names of the
replicate groups for which data is present in this object.

\item \code{groupNames}: Same as \code{names}. Provided for consistency to other classes.

\item \code{length}: Obtain number of feature groups.

\item \code{show}: Shows summary information for this object.

\item \code{groupTable}: Accessor for \code{groups} slot.

\item \code{analysisInfo}: Obtain analysisInfo (see analysisInfo slot in \code{\link{features}}).

\item \code{groupInfo}: Accessor for \code{groupInfo} slot.

\item \code{featureTable}: Obtain feature information (see \code{\link{features}}).

\item \code{getFeatures}: Accessor for \code{features} slot.

\item \code{groupFeatIndex}: Accessor for \code{ftindex} slot.

\item \code{groupQualities}: Accessor for \code{groupQualities} slot.

\item \code{groupScores}: Accessor for \code{groupScores} slot.

\item \code{annotations}: Accessor for \code{annotations} slot.

\item \code{adducts}: Returns a named \code{character} with adduct annotations assigned to each feature group (if
available).

\item \code{adducts<-}: Sets adduct annotations for feature groups.

\item \code{[}: Subset on analyses/feature groups.

\item \code{[[}: Extract intensity values.

\item \code{$}: Extract intensity values for a feature group.

\item \code{delete}: Completely deletes specified feature groups.

\item \code{export}: Exports feature groups to a \file{.csv} file that is readable to Bruker ProfileAnalysis (a
'bucket table'), Bruker TASQ (an analyte database) or that is suitable as input for the \verb{Targeted peak
detection} functionality of \href{http://mzmine.github.io/}{MZmine}.

\item \code{as.data.table}: Obtain a summary table (a \code{\link{data.table}}) with retention, \emph{m/z}, intensity
and optionally other feature data.

\item \code{unique}: Obtain a subset with unique feature groups
present in one or more specified replicate group(s).

\item \code{overlap}: Obtain a subset with feature groups that overlap
between a set of specified replicate group(s).

\item \code{calculatePeakQualities}: Calculates peak and group qualities for all features and feature groups. The peak qualities
(and scores) are calculated with the \link[=calculatePeakQualities,features-method]{features method of this
function}, and subsequently averaged per feature group. Then, \pkg{MetaClean} is used to calculate the
\verb{Elution Shift} and \verb{Retention Time Consistency} group quality metrics (see the \pkg{MetaClean}
publication cited below for more details). Similarly to the \code{\link{features}} method, these metrics are scored
by normalizing qualities among all groups and scaling them from \samp{0} (worst) to \samp{1} (best). The
\verb{totalScore} for each group is then calculated as the weighted sum from all feature (group) scores. The
\code{\link{getMCTrainData}} and \code{\link{predictCheckFeaturesSession}} functions can be used to train and apply
Pass/Fail ML models from \pkg{MetaClean}.

\item \code{selectIons}: uses \link[=generateComponents]{componentization} results to select feature groups with
preferred adduct ion and/or isotope annotation. Typically, this means that only feature groups are kept if they are
(de-)protonated adducts and are monoisotopic. The adduct annotation assignments for the selected feature groups are
copied from the components to the \code{annotations} slot. If the adduct for a feature group is unknown, its
annotation is defaulted to the 'preferred' adduct, and hence, the feature group will never be removed. Furthermore,
if a component does not contain an annotation with the preferred adduct, the most intense feature group is selected
instead. Similarly, if no isotope annotation is available, the feature group is assumed to be monoisotopic and thus
not removed. An important advantage of \code{selectIons} is that it may considerably simplify your dataset.
Furthermore, the adduct assignments allow formula/compound annotation steps later in the workflow to improve their
annotation accuracy. On the other hand, it is important the componentization results are reliable. Hence, it is
highly recommended that, prior to calling \code{selectIons}, the settings to \code{\link{generateComponents}} are
optimized and its results are reviewed with \code{\link{checkComponents}}. Finally, the \code{adducts<-} method can
be used to manually correct adduct assignments afterwards if necessary.
}}

\section{Slots}{

\describe{
\item{\code{groups}}{Matrix (\code{\link{data.table}}) with intensities for each feature group (columns) per analysis (rows).
Access with \code{groups} method.}

\item{\code{analysisInfo,features}}{\link[=analysis-information]{Analysis info} and \code{\link{features}} class associated
with this object. Access with \code{analysisInfo} and \code{featureTable} methods, respectively.}

\item{\code{groupInfo}}{\code{data.frame} with retention time (\code{rts} column, in seconds) and \emph{m/z} (\code{mzs}
column) for each feature group. Access with \code{groupInfo} method.}

\item{\code{ftindex}}{Matrix (\code{\link{data.table}}) with feature indices for each feature group (columns) per analysis
(rows). Each index corresponds to the row within the feature table of the analysis (see
\code{\link{featureTable}}).}

\item{\code{groupQualities,groupScores}}{A \code{\link{data.table}} with qualities/scores for each feature group (see the
\code{calculatePeakQualities} method).}

\item{\code{annotations}}{A \code{\link{data.table}} with adduct annotations for each group (see the \code{selectIons}
method).}

\item{\code{groupAlgo,groupArgs,groupVerbose}}{\setsWF Grouping parameters that were used when this object was created. Used
by \code{adducts<-} and \code{selectIons} when these methods perform a re-grouping of features.}

\item{\code{annotations}}{\setsWF As the \code{featureGroups} slot, but contains the annotation data per set.}
}}

\section{S4 class hierarchy}{
 \itemize{   \item{\code{\link{workflowStep}}}   \itemize{     \item{\strong{\code{\link{featureGroups}}}}     \itemize{       \item{\code{\link{featureGroupsSet}}}       \itemize{         \item{\code{\link{featureGroupsScreeningSet}}}       }       \item{\code{\link{featureGroupsUnset}}}       \item{\code{\link{featureGroupsScreening}}}       \itemize{         \item{\code{\link{featureGroupsSetScreeningUnset}}}       }       \item{\code{\link{featureGroupsBruker}}}       \item{\code{\link{featureGroupsConsensus}}}       \item{\code{\link{featureGroupsEnviMass}}}       \item{\code{\link{featureGroupsKPIC2}}}       \item{\code{\link{featureGroupsOpenMS}}}       \item{\code{\link{featureGroupsSIRIUS}}}       \item{\code{\link{featureGroupsBrukerTASQ}}}       \item{\code{\link{featureGroupsXCMS}}}       \item{\code{\link{featureGroupsXCMS3}}}     }   } }
}

\section{Sets workflows}{
 \setsWFClass{featureGroupsSet}{featureGroups}

  \setsWFNewMethodsFeat{featureGroupsUnset}{The adduct annotations for the selected set are used to convert all
  feature (group) masses to ionic \emph{m/z} values. The annotations persist in the converted object. }

  \setsWFChangedMethods{

  \item \code{adducts}, \code{adducts<-} require the \code{set} argument. The order of the data that is
  returned/changed follows that of the \code{annotations} slot. Furthermore, \code{adducts<-} will perform a
  re-grouping of features when its \code{reGroup} parameter is set to \code{TRUE}. The implications for this are
  discussed below.

  \item the subset operator (\code{[}) has specific arguments to choose (feature presence in) sets. See the argument
  descriptions.

  \item \code{as.data.table}: normalization of intensities is performed per set.

  \item \code{export} Only allows to export data from one set. The \code{unset} method is used prior to exporting the
  data.

  \item \code{overlap} and \code{unique} allow to handle data per set. See the \code{sets} argument description.

  \item \code{selectIons} Will perform a re-grouping of features. The implications of this are discussed below.

  }

  A re-grouping of features occurs if \code{selectIons} is called or \code{adducts<-} is used with
  \code{reGroup=TRUE}. Afterwards, it is very likely that feature group names are changed. Since data generated later
  in the workflow (\emph{e.g.} annotation steps) rely on feature group names, these objects are \strong{not valid}
  anymore, and \strong{must} be re-generated.
}

\references{
\insertRef{Chetnik2020}{patRoon}
}
\seealso{
\code{\link{groupFeatures}} for generating feature groups, \link{feature-filtering} and
  \link{feature-plotting} for more advanced \code{featureGroups} methods.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups-plot.R
\name{feature-plotting}
\alias{feature-plotting}
\alias{plot,featureGroups,missing-method}
\alias{plotInt,featureGroups-method}
\alias{plotInt,featureGroupsSet-method}
\alias{plotChord,featureGroups-method}
\alias{plotChroms,featureGroups-method}
\alias{plotVenn,featureGroups-method}
\alias{plotVenn,featureGroupsSet-method}
\alias{plotUpSet,featureGroups-method}
\alias{plotVolcano,featureGroups-method}
\title{Plotting of grouped features}
\usage{
\S4method{plot}{featureGroups,missing}(
  x,
  colourBy = c("none", "rGroups", "fGroups"),
  onlyUnique = FALSE,
  retMin = FALSE,
  showLegend = TRUE,
  col = NULL,
  pch = NULL,
  ...
)

\S4method{plotInt}{featureGroups}(
  obj,
  average = FALSE,
  normFunc = NULL,
  xnames = TRUE,
  showLegend = FALSE,
  pch = 20,
  type = "b",
  lty = 3,
  col = NULL,
  ...
)

\S4method{plotInt}{featureGroupsSet}(
  obj,
  average = FALSE,
  normFunc = NULL,
  xnames = !sets,
  showLegend = sets,
  pch = 20,
  type = "b",
  lty = 3,
  col = NULL,
  ...,
  sets = FALSE
)

\S4method{plotChord}{featureGroups}(
  obj,
  addSelfLinks = FALSE,
  addRetMzPlots = TRUE,
  average = FALSE,
  outerGroups = NULL,
  addIntraOuterGroupLinks = FALSE,
  ...
)

\S4method{plotChroms}{featureGroups}(
  obj,
  rtWindow = 30,
  mzExpWindow = 0.001,
  retMin = FALSE,
  topMost = NULL,
  topMostByRGroup = FALSE,
  EICs = NULL,
  showPeakArea = FALSE,
  showFGroupRect = TRUE,
  title = NULL,
  colourBy = c("none", "rGroups", "fGroups"),
  showLegend = TRUE,
  onlyPresent = TRUE,
  annotate = c("none", "ret", "mz"),
  showProgress = FALSE,
  xlim = NULL,
  ylim = NULL,
  ...
)

\S4method{plotVenn}{featureGroups}(obj, which = NULL, ...)

\S4method{plotVenn}{featureGroupsSet}(obj, which = NULL, ..., sets = FALSE)

\S4method{plotUpSet}{featureGroups}(obj, which = NULL, nsets = length(which), nintersects = NA, ...)

\S4method{plotVolcano}{featureGroups}(
  obj,
  FCParams,
  showLegend = TRUE,
  averageFunc = mean,
  col = NULL,
  pch = 19,
  ...
)
}
\arguments{
\item{x}{\code{featureGroups} object to be accessed.}

\item{colourBy}{Sets the automatic colour selection: \code{"none"} for a single colour or
\code{"rGroups"}/\code{"fGroups"} for a distinct colour per replicate/feature group.}

\item{onlyUnique}{If \code{TRUE} and \code{colourBy="rGroups"} then only
feature groups that are unique to a replicate group are plotted.}

\item{retMin}{Plot retention time in minutes (instead of seconds).}

\item{showLegend}{Plot a legend if \code{TRUE}.}

\item{col}{Colour(s) used. If \code{col=NULL} then colours are automatically generated.}

\item{pch, type, lty}{Common plotting parameters passed to \emph{e.g.} \code{\link[graphics]{plot}}. For \code{plot}:
if \code{pch=NULL} then values are automatically assigned.}

\item{\dots}{passed to \code{\link[base]{plot}} (\code{plot} and \code{plotChroms}), \code{\link[graphics]{lines}}
(\code{plotInt}), \pkg{\link{VennDiagram}} plotting functions (\code{plotVenn}), \code{\link{chordDiagram}}
(\code{plotChord}) or \code{\link[UpSetR]{upset}} (\code{plotUpSet}).}

\item{obj}{\code{featureGroups} object to be accessed.}

\item{average}{If \code{TRUE} then data within replicate groups are averaged.

  For \code{as.data.table}: if \code{features=TRUE} other feature properties are also averaged.}

\item{normFunc}{Function that should be used for normalization of data. The function is called for all
intensities/areas of a feature group and these quantities are divided by the result of the function call. For
example, when \code{\link{max}} is used normalized intensities will be between zero and one. If all quantities are
zero then the function will not be called. Set to \code{NULL} to perform no normalization.}

\item{xnames}{Plot analysis (or replicate group if \code{average=TRUE}) names on the x axis.}

\item{sets}{\setsWF For \code{plotInt}: if \code{TRUE} then feature intensities are plot per set (order follows the
  \link[=analysis-information]{analysis information}).

  For \code{plotVenn}: If \code{TRUE} then the \code{which} argument changes its meaning and is used to specify the
  names of the sets to be compared.}

\item{addSelfLinks}{If \code{TRUE} then 'self-links' are added which
represent non-shared data.}

\item{addRetMzPlots}{Set to \code{TRUE} to enable \emph{m/z} \emph{vs}
retention time scatter plots.}

\item{outerGroups}{Character vector of names to be used as outer groups. The
values in the specified vector should be named by analysis names
(\code{average} set to \code{FALSE}) or replicate group names
(\code{average} set to \code{TRUE}), for instance: \code{c(analysis1 =
"group1", analysis2 = "group1", analysis3 = "group2")}. Set to \code{NULL}
to disable outer groups.}

\item{addIntraOuterGroupLinks}{If \code{TRUE} then links will be added within
outer groups.}

\item{rtWindow}{Retention time (in seconds) that will be subtracted/added to respectively the minimum and maximum
retention time of the plotted feature groups. Thus, setting this value to a positive value will 'zoom out' on the
retention time axis.}

\item{mzExpWindow}{In case the \emph{m/z} window to plot an EIC for a particular analysis is not known (\emph{i.e.}
no feature was detected of the feature group to be plot and \code{onlyPresent=FALSE}) then the EIC \emph{m/z} range
is estimated from the range for the complete feature group and expanded by the offset defined by
\code{mzExpWindow}.}

\item{topMost}{Only plot EICs from features within this number of top most intense analyses. If \code{NULL} then all
analyses are used for plotted.}

\item{topMostByRGroup}{If set to \code{TRUE} and \code{topMost} is set: only plot EICs for the top most features in
each replicate group. For instance, when \code{topMost=1} and \code{topMostByRGroup=TRUE}, then EICs will be
plotted for the most intense feature of each replicate group.}

\item{EICs}{Internal parameter for now and should be kept at \code{NULL} (default).}

\item{showPeakArea}{Set to \code{TRUE} to display integrated chromatographic peak ranges by filling (shading) their
areas.}

\item{showFGroupRect}{Set to \code{TRUE} to mark the full retention/intensity range of all features within a feature
group by drawing a rectangle around it.}

\item{title}{Character string used for title of the plot. If \code{NULL} a title will be automatically generated.}

\item{onlyPresent}{If \code{TRUE} then EICs will only be generated for analyses in which a particular feature group
was detected. Disabling this option might be useful to see if any features were 'missed'.}

\item{annotate}{If set to \code{"ret"} and/or \code{"mz"} then retention and/or \emph{m/z} values will be drawn for
each plotted feature group.}

\item{showProgress}{if set to \code{TRUE} then a text progressbar will be displayed when all EICs are being plot. Set
to \code{"none"} to disable any annotation.}

\item{xlim, ylim}{Sets the plot size limits used by
\code{\link[graphics]{plot}}. Set to \code{NULL} for automatic plot sizing.}

\item{which}{A character vector with replicate groups used for comparison. Set to \code{NULL} to ignore.

  For \code{plotVenn}: alternatively a named \code{list} containing elements of \code{character} vectors with
  replicate groups to compare. For instance, \code{which=list(infl = c("influent-A", "influent-B"), effl =
  c("effluent-A", "effluent-B"))}, will compare the features in replicate groups \samp{"influent-A/B"} against those
  in \samp{"effluent-A/B"}. The names of the list are used for labelling in the plot, and will be made automatically
  if not specified.}

\item{nsets, nintersects}{See \code{\link[UpSetR]{upset}}.}

\item{FCParams}{A parameter list to calculate Fold change data. See \code{getFCParams} for more details.}

\item{averageFunc}{Function used for averaging.}
}
\value{
\code{plotVenn} (invisibly) returns a list with the following fields: \itemize{
\item \code{gList} the \code{gList} object that was returned by
  the utilized \pkg{\link{VennDiagram}} plotting function.
\item \code{areas} The total area for each plotted group.
\item \code{intersectionCounts} The number of intersections between groups.
}

The order for the \code{areas} and \code{intersectionCounts} fields is the same as the parameter order
from the used plotting function (see \emph{e.g.} \code{\link{draw.pairwise.venn}} and
\code{\link{draw.triple.venn}}).
}
\description{
Various plotting functions for feature group data.
}
\details{
\code{plot} Generates an \emph{m/z} \emph{vs} retention time
  plot for all featue groups. Optionally highlights unique/overlapping
  presence amongst replicate groups.

\code{plotInt} Generates a line plot for the (averaged) intensity
  of feature groups within all analyses

\code{plotChord} Generates a chord diagram which can be used to
  visualize shared presence of feature groups between analyses or replicate
  groups. In addition, analyses/replicates sharing similar properties
  (\emph{e.g.} location, age, type) may be grouped to enhance visualization
  between these 'outer groups'.

\code{plotChroms} Plots extracted ion chromatograms (EICs) of feature groups.

\code{plotVenn} plots a Venn diagram (using \pkg{\link{VennDiagram}}) outlining unique and shared feature
  groups between up to five replicate groups.

\code{plotUpSet} plots an UpSet diagram (using the \code{\link[UpSetR]{upset}} function) outlining unique
  and shared feature groups between given replicate groups.

\code{plotVolcano} Plots Fold change data in a 'Volcano plot'.
}
\section{Sets workflows}{
 \setsWFChangedMethods{

  \item \code{plotVenn} and \code{plotInt} allow to handle data per set. See the \code{sets} argument description.

  }
}

\references{
\addCitations{circlize}{1}

\insertRef{Conway2017}{patRoon} \cr\cr \insertRef{Lex2014}{patRoon}
}
\seealso{
\code{\link{featureGroups-class}}, \code{\link{groupFeatures}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups-kpic2.R
\name{groupFeaturesKPIC2}
\alias{groupFeaturesKPIC2}
\alias{groupFeaturesKPIC2,features-method}
\alias{groupFeaturesKPIC2,featuresSet-method}
\title{Group features using KPIC2}
\usage{
\S4method{groupFeaturesKPIC2}{features}(
  feat,
  rtalign = TRUE,
  loadRawData = TRUE,
  groupArgs = list(tolerance = c(0.005, 12)),
  alignArgs = list(),
  verbose = TRUE
)

\S4method{groupFeaturesKPIC2}{featuresSet}(
  feat,
  groupArgs = list(tolerance = c(0.005, 12)),
  verbose = TRUE
)
}
\arguments{
\item{feat}{The \code{\link{features}} object with the features to be grouped.}

\item{rtalign}{Set to \code{TRUE} to enable retention time alignment.}

\item{loadRawData}{Set to \code{TRUE} if analyses are available as \code{mzXML} or \code{mzML} files. Otherwise MS
data is not loaded, and some dummy data (\emph{e.g.} file paths) is used in the returned object.}

\item{groupArgs, alignArgs}{Named \code{character} vector that may contain extra parameters to be used by
\code{\link[KPIC:PICset.group]{KPIC::PICset.group}} and \code{\link[KPIC:PICset.align]{KPIC::PICset.align}},
respectively.}

\item{verbose}{if \code{FALSE} then no text output will be shown.}
}
\value{
An object of a class which is derived from \code{\link{featureGroups}}.

The \code{featuresSet} method (for \link[=sets-workflow]{sets workflows}) returns a
  \code{\link{featureGroupsSet}} object.
}
\description{
Uses the the \href{https://github.com/hcji/KPIC2}{KPIC2} \R package for grouping of features.
}
\details{
This function uses KPIC2 to group features. This function is called when calling \code{groupFeatures} with
  \code{algorithm="kpic2"}.

Grouping of features and alignment of their retention times are performed with the
  \code{\link[KPIC:PICset.group]{KPIC::PICset.group}} and \code{\link[KPIC:PICset.align]{KPIC::PICset.align}}
  functions, respectively.
}
\section{Sets workflows}{
 \code{loadRawData} and arguments related to retention time alignment are currently not
  supported for \link[=sets-workflow]{sets workflows}.
}

\references{
\insertRef{Ji2017}{patRoon}
}
\seealso{
\code{\link{groupFeatures}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups-xcms.R
\name{importFeatureGroupsXCMS}
\alias{importFeatureGroupsXCMS}
\title{Imports feature groups from XCMS (old interface)}
\usage{
importFeatureGroupsXCMS(xs, analysisInfo)
}
\arguments{
\item{xs}{An \code{\link{xcmsSet}} object.}

\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}
}
\value{
An object of a class which is derived from \code{\link{featureGroups}}.

The \code{featuresSet} method (for \link[=sets-workflow]{sets workflows}) returns a
  \code{\link{featureGroupsSet}} object.
}
\description{
Imports grouped features from a legacy \code{\link{xcmsSet}} object from the \pkg{xcms} package.
}
\references{
\addCitations{xcms}{1} \cr\cr \addCitations{xcms}{2} \cr\cr \addCitations{xcms}{3}
}
\seealso{
\code{\link{importFeaturesXCMS3}} and \code{\link{groupFeatures}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-exported.R
\name{verifyDependencies}
\alias{verifyDependencies}
\title{Verifies if all dependencies are installed properly and instructs the user if
this is not the case.}
\usage{
verifyDependencies()
}
\description{
Verifies if all dependencies are installed properly and instructs the user if
this is not the case.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features.R
\name{findFeatures}
\alias{findFeatures}
\title{Finding features}
\usage{
findFeatures(analysisInfo, algorithm, ..., verbose = TRUE)
}
\arguments{
\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}

\item{algorithm}{A character string describing the algorithm that should be
used: \code{"bruker"}, \code{"openms"}, \code{"xcms"}, \code{"xcms3"}, \code{"envipick"}, \code{"sirius"}, \code{"kpic2"}, \code{"safd"}}

\item{\dots}{Further parameters passed to the selected feature finding algorithms.}

\item{verbose}{If set to \code{FALSE} then no text output is shown.}
}
\value{
An object of a class which is derived from \code{\link{features}}.
}
\description{
Automatically find features.
}
\details{
Several functions exist to collect features (\emph{i.e.} retention and MS information that represent potential
compounds) from a set of analyses. All 'feature finders' return an object derived from the \code{\link{features}}
base class. The next step in a general workflow is to group and align these features across analyses with
\code{\link{groupFeatures}}. Note that some feature finders have a plethora of options which sometimes may have a
large effect on the quality of results. Fine-tuning parameters is therefore important, and the optimum is largely
dependent upon applied analysis methodology and instrumentation.

\code{findFeatures} is a generic function that will find features by one of the supported algorithms. The actual
  functionality is provided by algorithm specific functions such as \code{findFeaturesOpenMS} and \code{findFeaturesXCMS}. While these
  functions may be called directly, \code{findFeatures} provides a generic interface and is therefore usually preferred.
}
\note{
In most cases it will be necessary to centroid your MS input files. The only exception is \command{Bruker},
  however, you will still need centroided \file{mzXML}/\file{mzML} files for \emph{e.g.} plotting chromatograms. In
  this case the centroided MS files should be stored in the same directory as the raw \command{Bruker} \file{.d}
  files. The \code{\link{convertMSFiles}} function can be used to centroid data.
}
\seealso{
The \code{\link{features}} output class and its methods and the algorithm specific functions:
  \code{\link{findFeaturesBruker}}, \code{\link{findFeaturesOpenMS}}, \code{\link{findFeaturesXCMS}}, \code{\link{findFeaturesXCMS3}}, \code{\link{findFeaturesEnviPick}}, \code{\link{findFeaturesSIRIUS}}, \code{\link{findFeaturesKPIC2}}, \code{\link{findFeaturesSAFD}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/components-tps.R
\name{generateComponentsTPs}
\alias{generateComponentsTPs}
\alias{generateComponentsTPs,featureGroups-method}
\alias{generateComponentsTPs,featureGroupsSet-method}
\title{Generate components of transformation products}
\usage{
\S4method{generateComponentsTPs}{featureGroups}(
  fGroups,
  fGroupsTPs = fGroups,
  ignoreParents = FALSE,
  TPs = NULL,
  MSPeakLists = NULL,
  formulas = NULL,
  compounds = NULL,
  minRTDiff = 20,
  specSimParams = getDefSpecSimParams()
)

\S4method{generateComponentsTPs}{featureGroupsSet}(
  fGroups,
  fGroupsTPs = fGroups,
  ignoreParents = FALSE,
  TPs = NULL,
  MSPeakLists = NULL,
  formulas = NULL,
  compounds = NULL,
  minRTDiff = 20,
  specSimParams = getDefSpecSimParams()
)
}
\arguments{
\item{fGroups}{The input \code{\link{featureGroups}} for componentization. See \code{fGroupsTPs}.}

\item{fGroupsTPs}{A \code{\link{featureGroups}} object containing the feature groups that are expected to be
transformation products. If a distinction between parents and TPs is not yet known, \code{fGroupsTPs} should equal
the \code{fGroups} argument. Otherwise, \code{fGroups} should only contain the parent feature groups, and both
\code{fGroups} and \code{fGroupsTPs} \emph{must} be a subset of the same \code{\link{featureGroups}} object.}

\item{ignoreParents}{If \code{TRUE} then feature groups present in both \code{fGroups} and \code{fGroupsTPs} are not
considered as TPs.}

\item{TPs}{A \code{\link{transformationProducts}} object. Set to \code{NULL} to perform linking without this data.}

\item{MSPeakLists, formulas, compounds}{A \code{\link{MSPeakLists}}/\code{\link{formulas}}/\code{\link{compounds}}
object to calculate MS/MS or annotation similarities between parents and TPs. If \code{NULL} then this data is not
calculated. For more details see the \verb{Linking parents and transformation products} section below.}

\item{minRTDiff}{Minimum retention time (in seconds) difference between the parent and a TP to determine whether a TP
elutes prior/after the parent (to calculate \code{retDir} values, see Details in \link{componentsTPs}))}

\item{specSimParams}{A named \code{list} with parameters that influence the calculation of MS spectra similarities.
See the \link[=specSimParams]{spectral similarity parameters} documentation for more details.}
}
\value{
The components are stored in objects derived from \code{\link{componentsTPs}}.
}
\description{
Generates components by linking feature groups of transformation products and their parents.
}
\details{
This function uses transformation product screening to generate components. This function is called when calling \code{generateComponents} with
  \code{algorithm="tp"}.

This method typically employs data from \link[=generateTPs]{generated transformation products} to find
  parents and their TPs. However, this data is not necessary, and components can also be made based on MS/MS
  similarity and/or other annotation similarities between the parent and its TPs. For more details see the
  \verb{Linking parents and transformation products} section below.
}
\note{
The \code{shift} parameter of \code{specSimParams} is ignored by \code{generateComponentsTPs}, since it always
  calculates similarities with all supported options.
}
\section{Linking parents and transformation products}{
 Each component consists of feature groups that are considered
  to be transformation products for one parent (the parent that 'belongs' to the component can be retrieved with the
  \code{\link{componentInfo}} method). The parent feature groups are taken from the \code{fGroups} parameter, while
  the feature groups for TPs are taken from \code{fGroupsTPs}. If a feature group occurs in both variables, it may
  therefore be considered as both a parent or TP.

  If transformation product data is given, \emph{i.e.} the \code{TPs} argument is set, then a suspect screening of
  the TPs must be performed in advance (see \code{\link{screenSuspects}} and \code{\link{convertToSuspects}} to
  create the suspect list). Furthermore, if TPs were generated with \code{\link{generateTPsBioTransformer}} or
  \code{\link{generateTPsLibrary}} then the suspect screening must also include the parents (\emph{e.g.} by setting
  \code{includeParents=TRUE} when calling \code{convertToSuspects} or by amending results by setting
  \code{amend=TRUE} to \code{screenSuspects}). The suspect screening is necessary for the componentization algorithm
  to map the feature groups of the parent or TP. If the the suspect screening yields multiple TP hits, all will be
  reported. Similarly, if the suspect screening contains multiple hits for a parent, a component is made for each of
  the parent hits.

  In case no transformation product data is provided (\code{TPs=NULL}), the componentization algorithm simply assumes
  that each feature group from \code{fGroupsTPs} is a potential TP for every parent feature group in \code{fGroups}.
  For this reason, it is highly recommended to specify which feature groups are parents/TPs (see the
  \code{fGroupsTPs} argument description above) and \emph{crucial} that the data is post-processed, for instance by
  only retaining TPs that have high annotation similarity with their parents (see the
  \code{\link[=filter,componentsTPs-method]{filter}} method for \code{\link{componentsTPs}}).

  A typical way to distinguish which feature groups are parents or TPs from two different (groups of) samples is by
  calculating Fold Changes (see the \code{\link[=as.data.table,featureGroups-method]{as.data.table}} method for
  feature groups and \code{\link{plotVolcano}}). Of course, other statistical techniques from \R are also suitable.

  During componentization, several characteristics are calculated which may be useful for post-processing: \itemize{

  \item \code{specSimilarity}: the MS/MS spectral similarity between the feature groups of the TP and its parent
  (\samp{0-1}).

  \item \code{specSimilarityPrec},\code{specSimilarityBoth}: as \code{specSimilarity}, but calculated with binned
  data using the \code{"precursor"} and \code{"both"} method, respectively (see \link[=specSimParams]{MS spectral
  similarity parameters} for more details).

  \item \code{fragmentMatches} The number of MS/MS fragment formula annotations that overlap between the TP and
  parent. If both the \code{formulas} and \code{compounds} arguments are specified then the annotation data is pooled
  prior to calculation. Note that only unique matches are counted. Furthermore, note that annotations from \emph{all}
  candidates are considered, even if the formula/structure of the parent/TP is known. Hence, \code{fragmentMatches}
  is mainly useful when little or no chemical information is known on the parents/TPs, \emph{i.e.}, when
  \code{TPs=NULL} or originates from \code{\link{generateTPsLogic}}. Since annotations for all candidates are used,
  it is highly recommended that the annotation objects are first processed with the \code{\link{filter}} method, for
  instance, to select only the top ranked candidates.

  \item \code{neutralLossMatches} As \code{fragmentMatches}, but counting overlapping neutral loss formulae.

  \item \code{retDir} The retention time direction of the TP relative to its parent. See Details in
  \link{componentsTPs}. If TP data was specified, the expected direction is stored in \code{TP_retDir}.

  \item \code{retDiff},\code{mzDiff},\code{formulaDiff} The retention time, \emph{m/z} and formula difference between
  the parent and TP (latter only available if data TP formula is available).

  }
}

\section{Sets workflows}{
 In a \link[=sets-workflow]{sets workflow} the component tables are amended with extra
  information such as overall/specific set spectrum similarities. As sets data is mixed, transformation products are
  able to be linked with a parent, even if they were not measured in the same set.
}

\seealso{
\code{\link{generateComponents}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TP.R
\name{generateTPs}
\alias{generateTPs}
\title{Generation of transformation products (TPs)}
\usage{
generateTPs(algorithm, ...)
}
\arguments{
\item{algorithm}{A character string describing the algorithm that should be
used: \code{"biotransformer"}, \code{"logic"}, \code{"library"}}

\item{\dots}{Any parameters to be passed to the selected TP generation algorithm.}
}
\value{
A \code{\link{transformationProducts}} (derived) object containing all generated TPs.
}
\description{
Functionality to automatically obtain transformation products for a given set of parent compounds.
}
\details{
\code{generateTPs} is a generic function that will generate transformation products by one of the supported algorithms. The actual
  functionality is provided by algorithm specific functions such as \code{generateTPsBioTransformer} and \code{generateTPsLogic}. While these
  functions may be called directly, \code{generateTPs} provides a generic interface and is therefore usually preferred.
}
\seealso{
The \code{\link{transformationProducts}} output class and its methods and the algorithm specific functions:
  \code{\link{generateTPsBioTransformer}}, \code{\link{generateTPsLogic}}, \code{\link{generateTPsLibrary}}

In addition, the derived classes \code{\link{transformationProductsBT}} and
  \code{\link{transformationProductsLibrary}} for algorithm specific methods to post-process TP data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features-xcms.R
\name{importFeaturesXCMS}
\alias{importFeaturesXCMS}
\title{Imports features from XCMS (old interface)}
\usage{
importFeaturesXCMS(xs, analysisInfo)
}
\arguments{
\item{xs}{An \code{\link{xcmsSet}} object.}

\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}
}
\value{
An object of a class which is derived from \code{\link{features}}.
}
\description{
Imports feature data generated with the legacy \code{\link[xcms]{xcmsSet}} function from the \pkg{xcms} package.
}
\details{
This function imports data from XCMS. This function is called when calling \code{importFeatures} with
  \code{type="xcms"}.
}
\references{
\addCitations{xcms}{1} \cr\cr \addCitations{xcms}{2} \cr\cr \addCitations{xcms}{3}
}
\seealso{
\code{\link{importFeatures}} for more details and other algorithms.

\code{\link{importFeaturesXCMS3}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doe-optimizer.R
\docType{class}
\name{optimizationResult-class}
\alias{optimizationResult-class}
\alias{optimizationResult}
\alias{algorithm,optimizationResult-method}
\alias{length,optimizationResult-method}
\alias{lengths,optimizationResult-method}
\alias{show,optimizationResult-method}
\alias{plot,optimizationResult,missing-method}
\alias{optimizedParameters,optimizationResult-method}
\alias{optimizedParameters}
\alias{optimizedObject,optimizationResult-method}
\alias{optimizedObject}
\alias{scores,optimizationResult-method}
\alias{scores}
\alias{experimentInfo,optimizationResult-method}
\alias{experimentInfo}
\title{Class containing optimization results.}
\usage{
\S4method{algorithm}{optimizationResult}(obj)

\S4method{length}{optimizationResult}(x)

\S4method{lengths}{optimizationResult}(x, use.names = FALSE)

\S4method{show}{optimizationResult}(object)

\S4method{plot}{optimizationResult,missing}(
  x,
  paramSet,
  DoEIteration,
  paramsToPlot = NULL,
  maxCols = NULL,
  type = "contour",
  image = TRUE,
  contours = "colors",
  ...
)

\S4method{optimizedParameters}{optimizationResult}(object, paramSet = NULL, DoEIteration = NULL)

\S4method{optimizedObject}{optimizationResult}(object, paramSet = NULL)

\S4method{scores}{optimizationResult}(object, paramSet = NULL, DoEIteration = NULL)

\S4method{experimentInfo}{optimizationResult}(object, paramSet, DoEIteration)
}
\arguments{
\item{obj, x, object}{An \code{optimizationResult} object.}

\item{use.names}{Ignored.}

\item{paramSet}{Numeric index of the parameter set (\emph{i.e.} the first
parameter set gets index \samp{1}). For some methods optional: if
\code{NULL} the best will be selected.}

\item{DoEIteration}{Numeric index specifying the DoE iteration within the
specified \code{paramSet}. For some methods optional: if \code{NULL} the
best will be selected.}

\item{paramsToPlot}{Which parameters relations should be plot. If \code{NULL}
all will be plot. Alternatively, a \code{list} containing one or more
\code{character} vectors specifying each two parameters that should be
plotted. Finally, if only one pair should be plotted, can be a
\code{character} vector specifying both parameters.}

\item{maxCols}{Multiple parameter pairs are plotted in a grid. The maximum
number of columns can be set with this argument. Set to \code{NULL} for no
limit.}

\item{type}{The type of plots to be generated: \code{"contour"},
\code{"image"} or \code{"persp"}. The equally named functions will be
called for plotting.}

\item{image}{Passed to \code{\link{contour}} (if \code{type="contour"}).}

\item{contours}{Passed to \code{\link{persp}} (if \code{type="persp"}).}

\item{\dots}{Further arguments passed to \code{\link{contour}},
\code{\link{image}} or \code{\link{persp}} (depending on \code{type}).}
}
\description{
Objects from this class contain optimization results resulting from design of
experiment (DoE).
}
\details{
Objects from this class are returned by \code{\link{optimizeFeatureFinding}} and
\code{\link{optimizeFeatureGrouping}}.
}
\section{Methods (by generic)}{
\itemize{
\item \code{algorithm}: Returns the algorithm that was used for finding features.

\item \code{length}: Obtain total number of experimental design iterations performed.

\item \code{lengths}: Obtain number of experimental design iterations performed for each parameter set.

\item \code{show}: Shows summary information for this object.

\item \code{plot}: Generates response plots for all or a selected
set of parameters.

\item \code{optimizedParameters}: Returns parameter set yielding optimal
results. The \code{paramSet} and \code{DoEIteration} arguments can be
\code{NULL}.

\item \code{optimizedObject}: Returns the object (\emph{i.e.} a
\code{\link{features}} or \code{\link{featureGroups}} object) that was
generated with optimized parameters. The \code{paramSet} argument can be
\code{NULL}.

\item \code{scores}: Returns optimization scores. The
\code{paramSet} and \code{DoEIteration} arguments can be \code{NULL}.

\item \code{experimentInfo}: Returns a \code{list} with optimization
information from an DoE iteration.
}}

\section{Slots}{

\describe{
\item{\code{algorithm}}{A character specifying the algorithm that was optimized.}

\item{\code{paramSets}}{A \code{list} with detailed results from each parameter set
that was tested.}

\item{\code{bestParamSet}}{Numeric index of the parameter set yielding the best
response.}
}}

\examples{
\dontrun{
# ftOpt is an optimization object.

# plot contour of all parameter pairs from the first parameter set/iteration.
plot(ftOpt, paramSet = 1, DoEIteration = 1)

# as above, but only plot two parameter pairs
plot(ftOpt, paramSet = 1, DoEIteration = 1,
     paramsToPlot = list(c("mzPPM", "chromFWHM"), c("chromFWHM", "chromSNR")))

# plot 3d perspective plots
plot(ftOpt, paramSet = 1, DoEIteration = 1, type = "persp")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compounds.R
\name{generateCompounds}
\alias{generateCompounds}
\alias{generateCompounds,featureGroups-method}
\title{Automatic compound annotation}
\usage{
\S4method{generateCompounds}{featureGroups}(fGroups, MSPeakLists, algorithm, ...)
}
\arguments{
\item{fGroups}{\code{\link{featureGroups}} object which should be annotated. This should be the same or a subset of
the object that was used to create the specified \code{MSPeakLists}. In the case of a subset only the remaining
feature groups in the subset are considered.}

\item{MSPeakLists}{A \code{\link{MSPeakLists}} object that was generated for the supplied \code{fGroups}.}

\item{algorithm}{A character string describing the algorithm that should be
used: \code{"metfrag"}, \code{"sirius"}}

\item{\dots}{Any parameters to be passed to the selected compound generation algorithm.}
}
\value{
A \code{\link{compounds}} derived object containing all compound annotations.
}
\description{
Automatically perform chemical compound annotation for feature groups.
}
\details{
Several algorithms are provided to automatically perform compound annotation for feature groups. To this end,
measured masses for all feature groups are searched within online database(s) (\emph{e.g.}
\href{https://pubchem.ncbi.nlm.nih.gov/}{PubChem}) to retrieve a list of potential candidate chemical compounds.
Depending on the algorithm and its parameters, further scoring of candidates is then performed using, for instance,
matching of measured and theoretical isotopic patterns, presence within other data sources such as patent databases
and similarity of measured and in-silico predicted MS/MS fragments. Note that this process is often quite time
consuming, especially for large feature group sets. Therefore, this is often one of the last steps within the
workflow and not performed before feature groups have been prioritized.

\code{generateCompounds} is a generic function that will generateCompounds by one of the supported algorithms. The actual
  functionality is provided by algorithm specific functions such as \code{generateCompoundsMetFrag} and \code{generateCompoundsSIRIUS}. While these
  functions may be called directly, \code{generateCompounds} provides a generic interface and is therefore usually preferred.
}
\section{Scorings}{
 Each algorithm implements their own scoring system. Their names have been simplified and
  harmonized where possible. The \code{\link{compoundScorings}} function can be used to get an overview of both the
  algorithm specific and generic scoring names.
}

\section{Sets workflows}{
 With a \link[=sets-workflow]{sets workflow}, annotation is first performed for each set.
  This is important, since the annotation algorithms typically cannot work with data from mixed ionization modes. The
  annotation results are then combined to generate a \emph{sets consensus}: \itemize{

  \item The annotation tables for each feature group from the set specific data are combined. Rows with overlapping
  candidates (determined by the first-block \acronym{InChIKey}) are merged.

  \item Set specific data (\emph{e.g.} the ionic formula) is retained by renaming their columns with set specific
  names.

  \item The MS/MS fragment annotations (\code{fragInfo} column) from each set are combined.

  \item The scorings for each set are averaged to calculate overall scores.

  \item The candidates are re-ranked based on their average ranking among the set data (if a candidate is absent in a
  set it is assigned the poorest rank in that set).

  \item The coverage of each candidate among sets is calculated. Depending on the \code{setThreshold} and
  \code{setThresholdAnn} arguments, candidates with low abundance are removed.

  }
}

\seealso{
The \code{\link{compounds}} output class and its methods and the algorithm specific functions:
  \code{\link{generateCompoundsMetFrag}}, \code{\link{generateCompoundsSIRIUS}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{patRoon-deprecated}
\alias{patRoon-deprecated}
\alias{reportMD}
\alias{exportDAFiles}
\alias{plotEIC}
\alias{groups}
\alias{plotSpec}
\alias{formulaTable}
\alias{compoundTable}
\title{Deprecated and renamed functions.}
\usage{
reportMD(...)

exportDAFiles(
  anaInfo,
  format = "mzML",
  exportLine = TRUE,
  outPath = anaInfo$path,
  overWrite = FALSE
)

plotEIC(obj, ...)

groups(object, ...)

plotSpec(obj, ...)

formulaTable(...)

compoundTable(...)
}
\arguments{
\item{\dots}{Passed to successor function.}

\item{format}{The output format of exported files. Should be either
\code{"mzXML"}, \code{"mzML"} or \code{"mzData"}.}

\item{exportLine}{Export line spectra (\code{TRUE}) or profile spectra
(\code{FALSE}). Usually line spectra are preferred, since profile spectra
use signficantly more disk space and increase required memory during
processing.}

\item{outPath}{Character vector of output paths for exported analyses. Will
be recycled if necessary.}

\item{overWrite}{If \code{TRUE} existing files will be overwritten.}
}
\description{
Please do not use these functions anymore since they may be removed in the
future.
}
\details{
\code{reportMD} performs HTML reporting, please use
  \code{\link{reportHTML}} instead.

\code{exportDAFiles} will export a set of analyses either in
  \file{.mzXML} or \file{.mzML} formats.

Please use \code{\link{plotChroms}} instead.

Please use \code{\link{groupTable}} instead.

Please use \code{\link{plotSpectrum}} instead.

Please use \code{\link{annotations}} instead.

Please use \code{\link{annotations}} instead.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups-bruker.R
\name{importFeatureGroupsBrukerPA}
\alias{importFeatureGroupsBrukerPA}
\title{Imports feature groups from Bruker ProfileAnalysis}
\usage{
importFeatureGroupsBrukerPA(
  path,
  feat,
  rtWindow = 12,
  mzWindow = 0.005,
  intWindow = 5,
  warn = TRUE
)
}
\arguments{
\item{path}{The file path to a exported 'bucket table' \file{.txt} file from PA.}

\item{feat}{The \code{\link{features}} object obtained with \code{\link{findFeaturesBruker}}.}

\item{rtWindow, mzWindow, intWindow}{Search window values for retention time (seconds), \emph{m/z} (Da) and intensity
used to find back features within feature groups from PA (+/- the retention/mass/intensity value of a feature).}

\item{warn}{Warn about missing or duplicate features when relating them back from grouped features.}
}
\value{
An object of a class which is derived from \code{\link{featureGroups}}.

The \code{featuresSet} method (for \link[=sets-workflow]{sets workflows}) returns a
  \code{\link{featureGroupsSet}} object.
}
\description{
Imports a 'bucket table' produced by Bruker ProfileAnalysis (PA)
}
\details{
This function imports data from Bruker ProfileAnalysis. This function is called when calling \code{importFeatureGroups} with
  \code{type="brukerpa"}.

The 'bucket table' should be exported as \file{.txt} file. Please note that this function only supports
  features generated by \code{\link{findFeaturesBruker}} and it is \strong{crucial} that DataAnalysis files remain
  unchanged when features are collected and the bucket table is generated. Furthermore, please note that PA does not
  retain information about originating features for generated buckets. For this reason, this function tries to find
  back the original features and care must be taken to correctly specify search parameters (\code{rtWindow},
  \code{mzWindow}, \code{intWindow}).
}
\seealso{
\code{\link{importFeatureGroups}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features-envipick.R
\name{findFeaturesEnviPick}
\alias{findFeaturesEnviPick}
\title{Find features using enviPick}
\usage{
findFeaturesEnviPick(analysisInfo, ..., parallel = TRUE, verbose = TRUE)
}
\arguments{
\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}

\item{\dots}{Further parameters passed to \code{\link[enviPick]{enviPickwrap}}.}

\item{parallel}{If set to \code{TRUE} then code is executed in parallel through the \CRANpkg{futures} package. Please
see the parallelization section in the handbook for more details.}

\item{verbose}{If set to \code{FALSE} then no text output is shown.}
}
\value{
An object of a class which is derived from \code{\link{features}}.
}
\description{
Uses the \code{\link[enviPick]{enviPickwrap}} function from the \pkg{enviPick} R package to extract features.
}
\details{
This function uses enviPick to automatically find features. This function is called when calling \code{findFeatures} with
  \code{algorithm="envipick"}.

The input MS data files need to be centroided. The \code{\link{convertMSFiles}} function can be used to
  centroid data.
}
\note{
The analysis files must be in the \code{mzXML} format.
}
\seealso{
\code{\link{findFeatures}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/components.R
\name{generateComponents}
\alias{generateComponents}
\alias{generateComponents,featureGroups-method}
\title{Grouping feature groups in components}
\usage{
\S4method{generateComponents}{featureGroups}(fGroups, algorithm, ...)
}
\arguments{
\item{fGroups}{\code{\link{featureGroups}} object for which components should be generated.}

\item{algorithm}{A character string describing the algorithm that should be
used: \code{"ramclustr"}, \code{"camera"}, \code{"nontarget"}, \code{"intclust"}, \code{"openms"}, \code{"cliquems"}, \code{"specclust"}, \code{"tp"}}

\item{\dots}{Any parameters to be passed to the selected component generation algorithm.}
}
\value{
A \code{\link{components}} (derived) object containing all generated components.
}
\description{
Functionality to automatically group related feature groups (\emph{e.g.} isotopes, adducts and homologues) to assist
and simplify annotation.
}
\details{
Several algorithms are provided to group feature groups that are related in some (chemical) way to each other. How
feature groups are related depends on the algorithm: examples include adducts, statistics and parents/transformation
products. The linking of this data is generally useful for annotation purposes and reducing data complexity.

\code{generateComponents} is a generic function that will generateComponents by one of the supported algorithms. The actual
  functionality is provided by algorithm specific functions such as \code{generateComponentsRAMClustR} and \code{generateComponentsNontarget}. While these
  functions may be called directly, \code{generateComponents} provides a generic interface and is therefore usually preferred.
}
\section{Sets workflows}{
 In a \link[=sets-workflow]{sets workflow} the componentization data is generated differently
  depending on the used algorithm. Please see the details in the algorithm specific functions linked in the \verb{See Also} section.
}

\seealso{
The \code{\link{components}} output class and its methods and the algorithm specific functions:
  \code{\link{generateComponentsRAMClustR}}, \code{\link{generateComponentsCAMERA}}, \code{\link{generateComponentsNontarget}}, \code{\link{generateComponentsIntClust}}, \code{\link{generateComponentsOpenMS}}, \code{\link{generateComponentsCliqueMS}}, \code{\link{generateComponentsSpecClust}}, \code{\link{generateComponentsTPs}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formulas-genform.R
\name{generateFormulasGenForm}
\alias{generateFormulasGenForm}
\alias{generateFormulasGenForm,featureGroups-method}
\alias{generateFormulasGenForm,featureGroupsSet-method}
\title{Generate formula with GenForm}
\usage{
\S4method{generateFormulasGenForm}{featureGroups}(
  fGroups,
  MSPeakLists,
  relMzDev = 5,
  adduct = NULL,
  elements = "CHNOP",
  hetero = TRUE,
  oc = FALSE,
  extraOpts = NULL,
  calculateFeatures = TRUE,
  featThreshold = 0,
  featThresholdAnn = 0.75,
  absAlignMzDev = 0.002,
  MSMode = "both",
  isolatePrec = TRUE,
  timeout = 120,
  topMost = 50,
  batchSize = 8
)

\S4method{generateFormulasGenForm}{featureGroupsSet}(
  fGroups,
  MSPeakLists,
  relMzDev = 5,
  adduct = NULL,
  ...,
  setThreshold = 0,
  setThresholdAnn = 0
)
}
\arguments{
\item{fGroups}{\code{\link{featureGroups}} object for which formulae should be generated. This should be the same or
a subset of the object that was used to create the specified \code{MSPeakLists}. In the case of a subset only the
remaining feature groups in the subset are considered.}

\item{MSPeakLists}{An \code{\link{MSPeakLists}} object that was generated for the supplied \code{fGroups}.}

\item{relMzDev}{Maximum relative deviation between the measured and candidate formula \emph{m/z} values (in ppm).
Sets the \option{ppm} command line option.}

\item{adduct}{An \code{\link{adduct}} object (or something that can be converted to it with \code{\link{as.adduct}}).
  Examples: \code{"[M-H]-"}, \code{"[M+Na]+"}. If the \code{featureGroups} object has
  adduct annotations then these are used if \code{adducts=NULL}.

  \setsWF The \code{adduct} argument is not supported for sets workflows, since the
  adduct annotations will then always be used.}

\item{elements}{Elements to be considered for formulae calculation. This will heavily affects the number of
candidates! Always try to work with a minimal set by excluding elements you don't expect. Sets the \option{el}
command line option.}

\item{hetero}{Only consider formulae with at least one hetero atom. Sets the \option{het} commandline option.}

\item{oc}{Only consider organic formulae (\emph{i.e.} with at least one carbon atom). Sets the \option{oc}
commandline option.}

\item{extraOpts}{An optional character vector with any other command line options that will be passed to
\command{GenForm}. See the \verb{GenForm options} section for all available command line options.}

\item{calculateFeatures}{If \code{TRUE} fomulae are first calculated for all features
prior to feature group assignment (see \verb{Candidate assignment} in \code{\link{generateFormulas}}).}

\item{featThreshold}{If \code{calculateFeatures=TRUE}: minimum presence (\samp{0-1}) of a formula in all features
before it is considered as a candidate for a feature group. For instance, \code{featThreshold=0.75} dictates that a
formula should be present in at least 75\% of the features inside a feature group.}

\item{featThresholdAnn}{As \code{featThreshold}, but only considers features with annotations. For instance,
\code{featThresholdAnn=0.75} dictates that a formula should be present in at least 75\% of the features with
annotations inside a feature group. @param topMost Only keep this number of candidates
(per feature group) with highest score.}

\item{absAlignMzDev}{When the group formula annotation consensus is made from feature annotations, the \emph{m/z}
values of annotated MS/MS fragments may slightly deviate from those of the corresponding group MS/MS peak list. The
\code{absAlignMzDev} argument specifies the maximum \emph{m/z} window used to re-align the mass peaks.}

\item{MSMode}{Whether formulae should be generated only from MS data (\code{"ms"}), MS/MS data (\code{"msms"}) or
both (\code{"both"}). Selecting \code{"both"} will fall back to formula calculation with only MS data in case no
MS/MS data is available.}

\item{isolatePrec}{Settings used for isolation of precursor mass peaks and their isotopes. This isolation is highly
important for accurate isotope scoring of candidates, as non-relevant mass peaks will dramatically decrease the
score. The value of \code{isolatePrec} should either be a \code{list} with parameters (see the
\code{\link[=filter,MSPeakLists-method]{filter method}} for \code{MSPeakLists} for more details), \code{TRUE} for
default parameters or \code{FALSE} for no isolation (\emph{e.g.} when you already performed isolation with the
\code{filter} method). The \code{z} parameter (charge) is automatically deduced from the adduct used for annotation
(unless \code{isolatePrec=FALSE}), hence any custom \code{z} setting is ignored.}

\item{timeout}{Maximum time (in seconds) that a \command{GenForm} command is allowed to execute. If this time is
exceeded a warning is emitted and the command is terminated. See the notes section for more information on the need
of timeouts.}

\item{topMost}{Only keep this number of candidates (per feature group) with highest
score.}

\item{batchSize}{Maximum number of \command{GenForm} commands that should be run sequentially in each parallel
process. Combining commands with short runtimes (such as \command{GenForm}) can significantly increase parallel
performance. For more information see \code{\link{executeMultiProcess}}. Note that this is ignored if
\option{patRoon.MP.method="future"}.}

\item{\dots}{\setsWF Further arguments passed to the non-sets workflow method.}

\item{setThreshold}{\setsWF Minimum abundance for a candidate among all sets (\samp{0-1}). For instance, a value of
\samp{1} means that the candidate needs to be present in all the set data.}

\item{setThresholdAnn}{\setsWF As \code{setThreshold}, but only taking into account the set data that contain
annotations for the feature group of the candidate.}
}
\value{
A \code{\link{formulas}} object containing all generated formulae.
}
\description{
Uses \href{https://sourceforge.net/projects/genform/}{GenForm} to generate chemical formula candidates.
}
\details{
This function uses genform to generate formula candidates. This function is called when calling \code{generateFormulas} with
  \code{algorithm="genform"}.

When MS/MS data is available it will be used to score candidate formulae by presence of 'fitting' fragments.
}
\note{
This function always sets the \option{exist} and \option{oei} \command{GenForm} command line options.

  Formula calculation with \command{GenForm} may produce an excessive number of candidates for high \emph{m/z} values
  (\emph{e.g.} above 600) and/or many elemental combinations (set by \code{elements}). In this scenario formula
  calculation may need a very long time. Timeouts are used to avoid excessive computational times by terminating long
  running commands (set by the \code{timeout} argument).
}
\section{GenForm options}{
 Below is a list of options (generated by running \command{GenForm} without commandline
  options) which can be set by the \code{extraOpts} parameter.

 \preformatted{Formula calculation from MS and MS/MS data as described in
Meringer et al (2011) MATCH Commun Math Comput Chem 65: 259-290
Usage: GenForm ms=<filename> [msms=<filename>] [out=<filename>]
        [exist[=mv]] [m=<number>] [ion=-e|+e|-H|+H|+Na] [cha=<number>]
        [ppm=<number>] [msmv=ndp|nsse|nsae] [acc=<number>] [rej=<number>]
        [thms=<number>] [thmsms=<number>] [thcomb=<number>]
        [sort[=ppm|msmv|msmsmv|combmv]] [el=<elements> [oc]] [ff=<fuzzy formula>]
        [vsp[=<even|odd>]] [vsm2mv[=<value>]] [vsm2ap2[=<value>]] [hcf]
        [wm[=lin|sqrt|log]] [wi[=lin|sqrt|log]] [exp=<number>] [oei]
        [dbeexc=<number>] [ivsm2mv=<number>] [vsm2ap2=<number>]
        [oms[=<filename>]] [omsms[=<filename>]] [oclean[=<filename>]]
        [analyze [loss] [intens]] [dbe] [cm] [pc] [sc]
Explanation:
        ms      : filename of MS data (*.txt)
        msms    : filename of MS/MS data (*.txt)
        out     : output generated formulas
        exist   : allow only molecular formulas for that at least one
                  structural formula exists;overrides vsp, vsm2mv, vsm2ap2;
                  argument mv enables multiple valencies for P and S
        m       : experimental molecular mass (default: mass of MS basepeak)
        ion     : type of ion measured (default: M+H)
        ppm     : accuracy of measurement in parts per million (default: 5)
        msmv    : MS match value based on normalized dot product, normalized
                  sum of squared or absolute errors (default: nsae)
        acc     : allowed deviation for full acceptance of MS/MS peak in ppm
                  (default: 2)
        rej     : allowed deviation for total rejection of MS/MS peak in ppm
                  (default: 4)
        thms    : threshold for the MS match value
        thmsms  : threshold for the MS/MS match value
        thcomb  : threshold for the combined match value
        sort    : sort generated formulas according to mass deviation in ppm,
                  MS match value, MS/MS match value or combined match value
        el      : used chemical elements (default: CHBrClFINOPSSi)
        oc      : only organic compounds, i.e. with at least one C atom
        ff      : overwrites el and oc and uses fuzzy formula for limits of
                  element multiplicities
        het     : formulas must have at least one hetero atom
        vsp     : valency sum parity (even for graphical formulas)
        vsm2mv  : lower bound for valency sum - 2 * maximum valency
                  (>=0 for graphical formulas)
        vsm2ap2 : lower bound for valency sum - 2 * number of atoms + 2
                  (>=0 for graphical connected formulas)
        hcf     : apply Heuerding-Clerc filter
        wm      : m/z weighting for MS/MS match value
        wi      : intensity weighting for MS/MS match value
        exp     : exponent used, when wi is set to log
        oei     : allow odd electron ions for explaining MS/MS peaks
        dbeexc  : excess of double bond equivalent for ions
        ivsm2mv : lower bound for valency sum - 2 * maximum valency
                  for fragment ions
        ivsm2ap2: lower bound for valency sum - 2 * number of atoms + 2
                  for fragment ions
        oms     : write scaled MS peaks to output
        omsms   : write weighted MS/MS peaks to output
        oclean  : write explained MS/MS peaks to output
        analyze : write explanations for MS/MS peaks to output
        loss    : for analyzing MS/MS peaks write losses instead of fragments
        intens  : write intensities of MS/MS peaks to output
        dbe     : write double bond equivalents to output
        cm      : write calculated ion masses to output
        pc      : output match values in percent
        sc      : strip calculated isotope distributions
        noref   : hide the reference information
}
}

\section{Parallelization}{
 generateFormulasGenForm uses multiprocessing to parallelize
  computations. Please see the parallelization section in the handbook for
  more details and \link[=patRoon-package]{patRoon options} for configuration
  options.

 When \code{futures} are used for parallel processing (\code{patRoon.MP.method="future"}),
  calculations with \command{GenForm} are done with batch mode disabled (see \code{batchSize} argument), which
  generally limit overall performance.
}

\references{
\insertRef{Meringer2011}{patRoon}
}
\seealso{
\code{\link{generateFormulas}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features-sirius.R
\name{findFeaturesSIRIUS}
\alias{findFeaturesSIRIUS}
\title{Find features using SIRIUS}
\usage{
findFeaturesSIRIUS(analysisInfo, verbose = TRUE)
}
\arguments{
\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}

\item{verbose}{If set to \code{FALSE} then no text output is shown.}
}
\value{
An object of a class which is derived from \code{\link{features}}.
}
\description{
Uses \href{https://bio.informatik.uni-jena.de/software/sirius/}{SIRIUS} to find features.
}
\details{
This function uses SIRIUS to automatically find features. This function is called when calling \code{findFeatures} with
  \code{algorithm="sirius"}.

The features are collected by running the \command{lcms-align} \command{SIRIUS} command for every analysis.

  The MS files should be in the \file{mzML} or \file{mzXML} format. Furthermore, this algorithms requires the
  presence of (data-dependent) MS/MS data.

The input MS data files need to be centroided. The \code{\link{convertMSFiles}} function can be used to
  centroid data.
}
\section{Parallelization}{
 \code{findFeaturesSIRIUS} uses multiprocessing to parallelize
  computations. Please see the parallelization section in the handbook for
  more details and \link[=patRoon-package]{patRoon options} for configuration
  options.

 Note that for caching purposes, the analyses files must always exist on the local host
  computer, even if it is not participating in computations.
}

\references{
\insertRef{Dhrkop2019}{patRoon}
}
\seealso{
\code{\link{findFeatures}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compounds-cluster.R
\docType{class}
\name{compoundsCluster-class}
\alias{compoundsCluster-class}
\alias{compoundsCluster}
\alias{clusters,compoundsCluster-method}
\alias{cutClusters,compoundsCluster-method}
\alias{clusterProperties,compoundsCluster-method}
\alias{groupNames,compoundsCluster-method}
\alias{length,compoundsCluster-method}
\alias{lengths,compoundsCluster-method}
\alias{show,compoundsCluster-method}
\alias{[,compoundsCluster,ANY,missing,missing-method}
\alias{treeCut,compoundsCluster-method}
\alias{treeCutDynamic,compoundsCluster-method}
\alias{plot,compoundsCluster,missing-method}
\alias{getMCS,compoundsCluster-method}
\alias{plotStructure,compoundsCluster-method}
\alias{plotSilhouettes,compoundsCluster-method}
\title{Compounds cluster class}
\usage{
\S4method{clusters}{compoundsCluster}(obj)

\S4method{cutClusters}{compoundsCluster}(obj)

\S4method{clusterProperties}{compoundsCluster}(obj)

\S4method{groupNames}{compoundsCluster}(obj)

\S4method{length}{compoundsCluster}(x)

\S4method{lengths}{compoundsCluster}(x, use.names = TRUE)

\S4method{show}{compoundsCluster}(object)

\S4method{[}{compoundsCluster,ANY,missing,missing}(x, i, j, ..., drop = TRUE)

\S4method{treeCut}{compoundsCluster}(obj, k = NULL, h = NULL, groupName)

\S4method{treeCutDynamic}{compoundsCluster}(obj, maxTreeHeight, deepSplit, minModuleSize, groupName)

\S4method{plot}{compoundsCluster,missing}(
  x,
  ...,
  groupName,
  pal = "Paired",
  colourBranches = lengths(x)[groupName] < 50,
  showLegend = lengths(x)[groupName] < 20
)

\S4method{getMCS}{compoundsCluster}(obj, groupName, cluster)

\S4method{plotStructure}{compoundsCluster}(
  obj,
  groupName,
  cluster,
  width = 500,
  height = 500,
  withTitle = TRUE
)

\S4method{plotSilhouettes}{compoundsCluster}(obj, kSeq, groupName, pch = 16, type = "b", ...)
}
\arguments{
\item{obj, x, object}{A \code{compoundsCluster} object.}

\item{use.names}{A logical value specifying whether the returned vector
should be named with the feature group names.}

\item{i}{For \code{[}: A numeric or character value which is used to select feature groups by
their index or name, respectively (for the order/names see \code{groupNames()}). Can also be logical to perform logical selection
(similar to regular vectors). If missing all feature groups are selected.}

\item{\dots}{Further arguments passed directly to the plotting function
(\code{plot} or \code{\link{plot.dendrogram}}).}

\item{drop, j}{ignored.}

\item{k, h}{Desired number of clusters or tree height to be used for cutting
the dendrogram, respecitively. One or the other must be specified.
Analogous to \code{\link{cutree}}.}

\item{groupName}{A character specifying the feature group name.}

\item{maxTreeHeight, deepSplit, minModuleSize}{Arguments used by
\code{\link{cutreeDynamicTree}}.}

\item{pal}{Colour palette to be used from \pkg{\link{RColorBrewer}}.}

\item{colourBranches}{Whether branches from cut clusters (and their labels)
should be coloured. Might be slow with large numbers of clusters, hence,
the default is only \code{TRUE} when this is not the case.}

\item{showLegend}{If \code{TRUE} and \code{colourBranches} is also
\code{TRUE} then a legend will be shown which outlines cluster numbers and
their colours. By default \code{TRUE} for small amount of clusters to avoid
overflowing the plot.}

\item{cluster}{A numeric value specifying the cluster.}

\item{width, height}{The dimensions (in pixels) of the raster image that
should be plotted.}

\item{withTitle}{A logical value specifying whether a title should be added.}

\item{kSeq}{An integer vector containing the sequence that should be used for
average silhouette width calculation.}

\item{pch, type}{Passed to \code{\link[graphics]{plot}}.}
}
\value{
\code{cutTree} and \code{cutTreeDynamic} return the modified
  \code{compoundsCluster} object.

\code{getMCS} returns an \CRANpkg{rcdk} molecule object
  (\code{IAtomContainer}).
}
\description{
Objects from this class are used to store hierarchical clustering data of
candidate structures within \code{\link{compounds}} objects.
}
\details{
Objects from this type are returned by the \code{compounds} method for
\code{\link[=makeHCluster,compounds-method]{makeHCluster}}.
}
\section{Methods (by generic)}{
\itemize{
\item \code{clusters}: Accessor method to the \code{clusters} slot.
Returns a list that contains for each feature group an object as returned
by \code{\link{hclust}}.

\item \code{cutClusters}: Accessor method to the \code{cutClusters} slot.
Returns a list that contains for each feature group a vector with cluster
membership for each candidate (format as \code{\link{cutree}}).

\item \code{clusterProperties}: Returns a list with properties on how the
clustering was performed.

\item \code{groupNames}: returns a \code{character} vector with the names of the
feature groups for which data is present in this object.

\item \code{length}: Returns the total number of clusters.

\item \code{lengths}: Returns a \code{vector} with the number of
clusters per feature group.

\item \code{show}: Show summary information for this object.

\item \code{[}: Subset on feature groups.

\item \code{treeCut}: Manually (re-)cut a dendrogram that was
generated for a feature group.

\item \code{treeCutDynamic}: Automatically (re-)cut a dendrogram that was
generated for a feature group using the \code{\link{cutreeDynamicTree}}
function from \pkg{\link{dynamicTreeCut}}.

\item \code{plot}: Plot the dendrogram for clustered compounds of a
feature group. Clusters are highlighted using \CRANpkg{dendextend}.

\item \code{getMCS}: Calculates the maximum common substructure (MCS)
for all candidate structures within a specified cluster. This method uses
the \code{\link{get.mcs}} function from \CRANpkg{rcdk}.

\item \code{plotStructure}: Plots the maximum common substructure (MCS) for
all candidate structures within a specified cluster.

\item \code{plotSilhouettes}: Plots the average silhouette width when the
clusters are cut by a sequence of k numbers. The k value with the highest
value (marked in the plot) may be considered as the optimal number of
clusters.
}}

\section{Slots}{

\describe{
\item{\code{clusters}}{A \code{list} with \code{\link{hclust}} objects for each
feature group.}

\item{\code{dists}}{A \code{list} with distance matrices for each feature group.}

\item{\code{SMILES}}{A \code{list} containing a vector with \code{SMILES} for all
candidate structures per feature group.}

\item{\code{cutClusters}}{A \code{list} with assigned clusters for all candidates per
feature group (same format as what \code{\link{cutree}} returns).}

\item{\code{properties}}{A list containing general properties and parameters used for
clustering.}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups-kpic2.R
\name{importFeatureGroupsKPIC2}
\alias{importFeatureGroupsKPIC2}
\title{Imports feature groups from KPIC2}
\usage{
importFeatureGroupsKPIC2(picsSetGrouped, analysisInfo)
}
\arguments{
\item{picsSetGrouped}{A grouped \code{PIC set} object (\emph{e.g.} as returned by
\code{\link[KPIC:PICset.group]{KPIC::PICset.group}}).}

\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}
}
\value{
An object of a class which is derived from \code{\link{featureGroups}}.

The \code{featuresSet} method (for \link[=sets-workflow]{sets workflows}) returns a
  \code{\link{featureGroupsSet}} object.
}
\description{
Imports grouped features from an \pkg{KPIC} object.
}
\references{
\insertRef{Ji2017}{patRoon}
}
\seealso{
\code{\link{groupFeatures}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mspeaklists-mzr.R
\name{generateMSPeakListsMzR}
\alias{generateMSPeakListsMzR}
\alias{generateMSPeakListsMzR,featureGroups-method}
\alias{generateMSPeakListsMzR,featureGroupsSet-method}
\title{Generate peak lists with mzR}
\usage{
\S4method{generateMSPeakListsMzR}{featureGroups}(
  fGroups,
  maxMSRtWindow = 5,
  precursorMzWindow = 4,
  topMost = NULL,
  avgFeatParams = getDefAvgPListParams(),
  avgFGroupParams = getDefAvgPListParams()
)

\S4method{generateMSPeakListsMzR}{featureGroupsSet}(fGroups, ...)
}
\arguments{
\item{fGroups}{The \code{\link{featureGroups}} object from which MS peak lists should be extracted.}

\item{maxMSRtWindow}{Maximum chromatographic peak window used for spectrum averaging (in seconds, +/- retention
time). If \code{NULL} all spectra from a feature will be taken into account. Lower to decrease processing time.}

\item{precursorMzWindow}{The \emph{m/z} window (in Da) to find MS/MS spectra of a precursor. This is typically used
for Data-Dependent like MS/MS data and should correspond to the isolation \emph{m/z} window (\emph{i.e.} +/- the
precursor \emph{m/z}) that was used to collect the data. For Data-Independent MS/MS experiments, where precursor
ions are not isolated prior to fragmentation (\emph{e.g.} bbCID, MSe, all-ion, ...) the value should be
\code{NULL}.}

\item{topMost}{Only extract MS peak lists from a maximum of \code{topMost} analyses with highest intensity. If
\code{NULL} all analyses will be used.}

\item{avgFeatParams}{Parameters used for averaging MS peak lists of individual features. Analogous to
\code{avgFGroupParams}.}

\item{avgFGroupParams}{A \code{list} with parameters used for averaging of peak lists for feature groups. See
\code{\link{getDefAvgPListParams}} for more details.}

\item{\dots}{\setsWF Further arguments passed to the non-sets workflow method.}
}
\value{
A \code{\link{MSPeakLists}} object.
}
\description{
Uses the \pkg{mzR} package to read the MS data needed for MS peak lists.
}
\details{
This function uses mzR to generate MS peak lists. This function is called when calling \code{generateMSPeakLists} with
  \code{algorithm="mzr"}.

The MS data files should be either in \file{.mzXML} or \file{.mzML} format.

The input MS data files need to be centroided. The \code{\link{convertMSFiles}} function can be used to
  centroid data.
}
\references{
\addCitations{mzR}{1} \cr\cr

  \addCitations{mzR}{2} \cr\cr

  \addCitations{mzR}{3} \cr\cr

  \addCitations{mzR}{4} \cr\cr

  \addCitations{mzR}{5} \cr\cr
}
\seealso{
\code{\link{generateMSPeakLists}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_groups-sirius.R
\name{groupFeaturesSIRIUS}
\alias{groupFeaturesSIRIUS}
\title{Group features using SIRIUS}
\usage{
groupFeaturesSIRIUS(analysisInfo, verbose = TRUE)
}
\arguments{
\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}

\item{verbose}{if \code{FALSE} then no text output will be shown.}
}
\value{
An object of a class which is derived from \code{\link{featureGroups}}.

The \code{featuresSet} method (for \link[=sets-workflow]{sets workflows}) returns a
  \code{\link{featureGroupsSet}} object.
}
\description{
Uses \href{https://bio.informatik.uni-jena.de/software/sirius/}{SIRIUS} to find \emph{and} group features.
}
\details{
This function uses SIRIUS to group features. This function is called when calling \code{groupFeatures} with
  \code{algorithm="sirius"}.

Finding and grouping features is done by running the \command{lcms-align} command on every analyses at once.
  For this reason, grouping feature data from other algorithms than \command{SIRIUS} is not supported.

  The MS files should be in the \file{mzML} or \file{mzXML} format. Furthermore, this algorithms requires the
  presence of (data-dependent) MS/MS data.

The input MS data files need to be centroided. The \code{\link{convertMSFiles}} function can be used to
  centroid data.
}
\references{
\insertRef{Dhrkop2019}{patRoon}
}
\seealso{
\code{\link{groupFeatures}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-exported.R
\name{withOpt}
\alias{withOpt}
\title{Temporarily changes package options}
\usage{
withOpt(code, ..., prefix = "patRoon.")
}
\arguments{
\item{code}{The code to be executed.}

\item{\dots}{Named arguments with options to change.}

\item{prefix}{A \code{character} that will be used to prefix given option
names.}
}
\description{
This function is inspired by
\code{\link[withr:with_options]{withr::with_options}}: it can be used to
execute some code where package options are temporarily changed. This
function uses a shortened syntax, especially when changing options for
\code{patRoon}.
}
\examples{
\dontrun{
# Set max parallel processes to five while performing formula calculations
withOpt(MP.maxProcs = 5, {
    formulas <- generateFormulas(fGroups, "genform", ...)
})
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/components-nontarget.R
\name{generateComponentsNontarget}
\alias{generateComponentsNontarget}
\alias{generateComponentsNontarget,featureGroups-method}
\alias{generateComponentsNontarget,featureGroupsSet-method}
\title{Componentization of homologous series with nontarget}
\usage{
\S4method{generateComponentsNontarget}{featureGroups}(
  fGroups,
  ionization = NULL,
  rtRange = c(-120, 120),
  mzRange = c(5, 120),
  elements = c("C", "H", "O"),
  rtDev = 30,
  absMzDev = 0.002,
  absMzDevLink = absMzDev * 2,
  traceHack = all(R.Version()[c("major", "minor")] >= c(3, 4)),
  ...
)

\S4method{generateComponentsNontarget}{featureGroupsSet}(fGroups, ionization = NULL, ...)
}
\arguments{
\item{fGroups}{\code{\link{featureGroups}} object for which components should be generated.}

\item{ionization}{Which ionization polarity was used to generate the data: should be \code{"positive"}
  or \code{"negative"}. If the \code{featureGroups} object has adduct annotations, and \code{ionization=NULL}, the
  ionization will be detected automatically.

  \setsWF This parameter is not supported for sets workflows, as the ionization will always be detected
  automatically.}

\item{rtRange}{A numeric vector containing the minimum and maximum retention time (in seconds) between homologues.
Series are always considered from low to high \emph{m/z}, thus, a negative minimum retention time allows detection
of homologous series with increasing \emph{m/z} and decreasing retention times. These values set the \code{minrt}
and \code{maxrt} arguments of \code{\link{homol.search}}.}

\item{mzRange}{A numeric vector specifying the minimum and maximum \emph{m/z} increment of a homologous series. Sets
the \code{minmz} and \code{maxmz} arguments of \code{\link{homol.search}}.}

\item{elements}{A character vector with elements to be considered for detection of repeating units. Sets the
\code{elements} argument of \code{\link{homol.search}} function.}

\item{rtDev}{Maximum retention time deviation. Sets the \code{rttol} to \code{\link{homol.search}}.}

\item{absMzDev}{Maximum absolute \emph{m/z} deviation. Sets the \code{mztol} argument to \code{\link{homol.search}}}

\item{absMzDevLink}{Maximum absolute \emph{m/z} deviation when linking series. This should usually be a bit higher
than \code{absMzDev} to ensure proper linkage.}

\item{traceHack}{Currently \code{\link{homol.search}} does not work with \R \samp{>3.3.3}. This flag, which is
enabled by default on these R versions, implements a (messy) workaround
(\href{https://github.com/blosloos/nontarget/issues/6}{more details here}).}

\item{\dots}{Any further arguments passed to \code{\link{homol.search}}.\cr\cr \setsWF Further arguments passed to the non-sets workflow method.}
}
\value{
The generated comnponents are returned as an object from the \code{\link{componentsNT}} class.
}
\description{
Uses \href{https://cran.r-project.org/web/packages/nontarget/index.html}{the nontarget R package} to generate
components by unsupervised detection of homologous series.
}
\details{
This function uses nontarget to generate components. This function is called when calling \code{generateComponents} with
  \code{algorithm="nontarget"}.

In the first step the \code{\link{homol.search}} function is used to detect all homologous series within
  each replicate group (analyses within each replicate group are averaged prior to detection). Then, homologous
  series across replicate groups are merged in case of full overlap or when merging of partial overlapping series
  causes no conflicts.
}
\section{Sets workflows}{
 In a \link[=sets-workflow]{sets workflow} the componentization is first performed for each
  set independently. The resulting components are then all combined in a \code{\link{componentsNTSet}} object. Note that
  the components themselves are never merged. The components are renamed to include the set name from which they were
  generated (\emph{e.g.} \code{"CMP1"} becomes \code{"CMP1-positive"}).

 The output class supports additional methods such as \code{plotGraph}.
}

\references{
\addCitations{nontarget}{1} \cr\cr \addCitations{enviPat}{1}
}
\seealso{
\code{\link{generateComponents}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features-xcms.R
\name{findFeaturesXCMS}
\alias{findFeaturesXCMS}
\title{Find features using XCMS (old interface)}
\usage{
findFeaturesXCMS(analysisInfo, method = "centWave", ..., verbose = TRUE)
}
\arguments{
\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}

\item{method}{The method setting used by XCMS peak finding, see \code{\link[xcms:findPeaks-methods]{xcms::findPeaks}}}

\item{\dots}{Further parameters passed to \code{\link[xcms]{xcmsSet}}.}

\item{verbose}{If set to \code{FALSE} then no text output is shown.}
}
\value{
An object of a class which is derived from \code{\link{features}}.
}
\description{
Uses the legacy \code{\link[xcms]{xcmsSet}} function from the \pkg{xcms} package to find features.
}
\details{
This function uses XCMS to automatically find features. This function is called when calling \code{findFeatures} with
  \code{algorithm="xcms"}.

This function uses the legacy interface of \pkg{xcms}. It is recommended to use
  \code{\link{findFeaturesXCMS3}} instead.

  The file format of analyses must be \code{mzML} or \code{mzXML}.

The input MS data files need to be centroided. The \code{\link{convertMSFiles}} function can be used to
  centroid data.
}
\references{
\addCitations{xcms}{1} \cr\cr \addCitations{xcms}{2} \cr\cr \addCitations{xcms}{3}
}
\seealso{
\code{\link{findFeatures}} for more details and other algorithms.

\code{\link{findFeaturesXCMS3}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/components-openms.R
\name{generateComponentsOpenMS}
\alias{generateComponentsOpenMS}
\alias{generateComponentsOpenMS,featureGroups-method}
\alias{generateComponentsOpenMS,featureGroupsSet-method}
\title{Componentization of adducts, isotopes etc. with OpenMS}
\usage{
\S4method{generateComponentsOpenMS}{featureGroups}(
  fGroups,
  ionization = NULL,
  chargeMin = 1,
  chargeMax = 1,
  chargeSpan = 3,
  qTry = "heuristic",
  potentialAdducts = NULL,
  minRTOverlap = 0.66,
  retWindow = 1,
  absMzDev = 0.005,
  minSize = 2,
  relMinAdductAbundance = 0.75,
  adductConflictsUsePref = TRUE,
  NMConflicts = c("preferential", "mostAbundant", "mostIntense"),
  prefAdducts = c("[M+H]+", "[M-H]-"),
  extraOpts = NULL
)

\S4method{generateComponentsOpenMS}{featureGroupsSet}(
  fGroups,
  ionization = NULL,
  chargeMin = 1,
  chargeMax = 1,
  chargeSpan = 3,
  qTry = "heuristic",
  potentialAdducts = NULL,
  ...
)
}
\arguments{
\item{fGroups}{\code{\link{featureGroups}} object for which components should be generated.}

\item{ionization}{Which ionization polarity was used to generate the data: should be \code{"positive"}
  or \code{"negative"}. If the \code{featureGroups} object has adduct annotations, and \code{ionization=NULL}, the
  ionization will be detected automatically.

  \setsWF This parameter is not supported for sets workflows, as the ionization will always be detected
  automatically.}

\item{chargeMin, chargeMax}{The minimum/maximum charge to consider. Corresponds to the
\command{algorithm:MetaboliteFeatureDeconvolution:charge_min}/\command{algorithm:MetaboliteFeatureDeconvolution:charge_min}
 options.}

\item{chargeSpan}{The maximum charge span for a single analyte. Corresponds to
\command{algorithm:MetaboliteFeatureDeconvolution:charge_span_max}.}

\item{qTry}{Sets how charges are determined. Corresponds to \command{algorithm:MetaboliteFeatureDeconvolution:q_try}.
Valid options are \code{"heuristic"} and \code{"all"} (the \code{"feature"} option from \command{OpenMS} is
currently not supported).}

\item{potentialAdducts}{The adducts to consider. Should be a \code{numeric} vector with probabilities for each
  adduct, \emph{e.g.} \code{potentialAdducts=c("[M+H]+" = 0.8, "[M+Na]+" = 0.2)}. Note that the sum of probabilities
  should always be \samp{1}. Furthermore, note that additions of multiple adducts should be controlled by the
  \code{chargeMin}/\code{chargeMax} arguments (and \emph{not} with \code{potentialAdducts}), \emph{e.g.} if
  \code{chargeMax=2} then both \code{[M+H]+} and \code{[2M+H]2+} may be considered. Please see the
  \command{algorithm:MetaboliteFeatureDeconvolution:potential_adducts} option of
  \href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/UTILS_MetaboliteAdductDecharger.html}{MetaboliteAdductDecharger}
   for more details. If \code{NULL} then the a default is chosen with \code{\link{defaultOpenMSAdducts}} (which is
  \emph{not} the same as \command{OpenMS}).

  \setsWF Should be a \code{list} where each entry specifies the potential adducts for a set. Should either be named
  with the sets names or follow the same order as \code{sets(fGroups)}. Example:
  \code{potentialAdducts=list(positive=c("[M+H]+" = 0.8, "[M+Na]+" = 0.2), negative=c("[M-H]-" = 0.8, "[M-H2O-H]-" =
  0.2))}}

\item{minRTOverlap, retWindow}{Sets feature retention tolerances when grouping features. Sets the
\command{"algorithm:MetaboliteFeatureDeconvolution:retention_max_diff"} and
\command{algorithm:MetaboliteFeatureDeconvolution:min_rt_overlap} options.}

\item{absMzDev}{Maximum absolute \emph{m/z} deviation. Sets the \command{algorithm:MetaboliteFeatureDeconvolution:mass_max_diff} option}

\item{minSize}{The minimum size of a component. Smaller components than this size will be removed. See note below.}

\item{relMinAdductAbundance}{The minimum relative abundance (\samp{0-1}) that an adduct should be assigned to
features within the same feature group. See the \verb{Feature components} section for more details.}

\item{adductConflictsUsePref}{If set to \code{TRUE}, and not all adduct assigments to the features within a feature
group are equal and at least one of those adducts is a preferential adduct (\code{prefAdducts} argument), then only
the features with (the lowest ranked) preferential adduct are considered. In all other cases or when
\code{adductConflictsUsePref=FALSE} only features with the most frequently assigned adduct is considered. See the
\verb{Feature components} section for more details.}

\item{NMConflicts}{The strategies to employ when not all neutral masses within a component are equal. Valid options
are: \code{"preferential"}, \code{"mostAbundant"} and \code{"mostIntense"}. Multiple strategies are possible, and
will be executed in the given order until one succeeds. See the \verb{Feature components} section for more details.}

\item{prefAdducts}{A \code{character} vector with one or more \emph{preferential adducts}. See the \verb{Feature
components} section for more details.}

\item{extraOpts}{Named character vector with extra command line parameters directly passed to
\command{MetaboliteAdductDecharger}. Set to \code{NULL} to ignore.}

\item{\dots}{\setsWF Further arguments passed to the non-sets workflow method.}
}
\value{
A \code{\link{componentsFeatures}} derived object.
}
\description{
Uses the
\href{https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/UTILS_MetaboliteAdductDecharger.html}{MetaboliteAdductDecharger}
utility (see \url{http://www.openms.de}) to generate components.
}
\details{
This function uses OpenMS to generate components. This function is called when calling \code{generateComponents} with
  \code{algorithm="openms"}.

Features that show highly similar chromatographic elution profiles are grouped, and subsequently annotated
  with their adducts.
}
\section{Feature components}{
 The returned components are based on so called \emph{feature components}. Unlike other
  algorithms, components are first made on a feature level (per analysis), instead of for complete feature groups. In
  the final step the feature components are converted to 'regular' components by employing a consensus approach with
  the following steps:

  \enumerate{

  \item If an adduct assigned to a feature only occurs as a minority compared to other adduct assigments within the
  same feature group, it is considered as an outlier and removed accordingly (controlled by the
  \code{relMinAdductAbundance} argument).

  \item For features within a feature group, only keep their adduct assignment if it occurs as the most frequent or
  is preferential (controlled by \code{adductConflictsUsePref} and \code{prefAdducts} arguments).

  \item Components are made by combining the feature groups for which at least one of their features are jointly
  present in the same feature component.

  \item Conflicts of neutral mass assignments within a component (\emph{i.e.} not all are the same) are dealt with.
  Firstly, all feature groups with an unknown neutral mass are split in another component. Then, if conflicts still
  occur, the feature groups with similar neutral mass (determined by \code{absMzDev} argument) are grouped. Depending
  on the \code{NMConflicts} argument, the group with one or more preferential adduct(s) or that is the largest or
  most intense is selected, whereas others are removed from the component. In case multiple groups contain
  preferential adducts, and \samp{>1} preferential adducts are available, the group with the adduct that matches
  first in \code{prefAdducts} 'wins'. In case of ties, one of the next strategies in \code{NMConflicts} is tried.

  \item If a feature group occurs in multiple components it will be removed completely.

  \item the \code{minSize} filter is applied.

  }
}

\section{Sets workflows}{
 In a \link[=sets-workflow]{sets workflow} the componentization is first performed for each
  set independently. The resulting components are then all combined in a \code{\link{componentsSet}} object. Note that
  the components themselves are never merged. The components are renamed to include the set name from which they were
  generated (\emph{e.g.} \code{"CMP1"} becomes \code{"CMP1-positive"}).
}

\section{Parallelization}{
 generateComponentsOpenMS uses multiprocessing to parallelize
  computations. Please see the parallelization section in the handbook for
  more details and \link[=patRoon-package]{patRoon options} for configuration
  options.
}

\references{
\insertRef{Bielow2010}{patRoon}
}
\seealso{
\code{\link{generateComponents}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/components-intclust.R
\docType{class}
\name{componentsIntClust-class}
\alias{componentsIntClust-class}
\alias{componentsIntClust}
\alias{plotHeatMap,componentsIntClust-method}
\alias{plotHeatMap}
\alias{plotInt,componentsIntClust-method}
\title{Components based on clustered intensity profiles.}
\usage{
\S4method{plotHeatMap}{componentsIntClust}(
  obj,
  interactive = FALSE,
  col = NULL,
  margins = c(6, 2),
  cexCol = 1,
  ...
)

\S4method{plotInt}{componentsIntClust}(obj, index, pch = 20, type = "b", lty = 3, col = NULL, ...)
}
\arguments{
\item{obj}{A \code{componentsIntClust} object.}

\item{interactive}{If \code{TRUE} an interactive heatmap will be drawn (with
\code{\link{heatmaply}}).}

\item{col}{The colour used for plotting. Set to \code{NULL} for automatic colours.}

\item{margins, cexCol}{Passed to \code{\link{heatmap.2}}}

\item{\dots}{Further options passed to \code{\link{heatmap.2}} / \code{\link{heatmaply}} (\code{plotHeatMap}),
\code{\link[graphics]{plot}} (\code{plotInt}).}

\item{index}{Numeric component/cluster index.}

\item{pch, type, lty}{Passed to \code{\link{lines}}.}
}
\value{
\code{plotHeatMap} returns the same as \code{\link{heatmap.2}} or
  \code{\link{heatmaply}}.
}
\description{
This class is derived from \code{\link{componentsClust}} and is used to store hierarchical clustering information
from intensity profiles of feature groups.
}
\details{
Objects from this class are generated by \code{\link{generateComponentsIntClust}}
}
\section{Methods (by generic)}{
\itemize{
\item \code{plotHeatMap}: draws a heatmap using the
\code{\link{heatmap.2}} or \code{\link{heatmaply}} function.

\item \code{plotInt}: makes a plot for all (normalized) intensity
profiles of the feature groups within a given cluster.
}}

\section{Slots}{

\describe{
\item{\code{clusterm}}{Numeric matrix with normalized feature group intensities that was used for clustering.}
}}

\note{
When the object is altered (\emph{e.g.} by filtering or subsetting it), methods that need the original
  clustered data such as plotting methods do not work anymore and stop with an error.
}
\section{S4 class hierarchy}{
 \itemize{   \item{\code{\link{componentsClust}}}   \itemize{     \item{\strong{\code{\link{componentsIntClust}}}}   } }
}

\seealso{
\code{\link{componentsClust}} for other relevant methods and \code{\link{generateComponents}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/components-cliquems.R
\name{generateComponentsCliqueMS}
\alias{generateComponentsCliqueMS}
\alias{generateComponentsCliqueMS,featureGroups-method}
\alias{generateComponentsCliqueMS,featureGroupsSet-method}
\title{Componentization of adducts, isotopes etc. with cliqueMS}
\usage{
\S4method{generateComponentsCliqueMS}{featureGroups}(
  fGroups,
  ionization = NULL,
  maxCharge = 1,
  maxGrade = 2,
  ppm = 10,
  adductInfo = NULL,
  absMzDev = 0.005,
  minSize = 2,
  relMinAdductAbundance = 0.75,
  adductConflictsUsePref = TRUE,
  NMConflicts = c("preferential", "mostAbundant", "mostIntense"),
  prefAdducts = c("[M+H]+", "[M-H]-"),
  extraOptsCli = NULL,
  extraOptsIso = NULL,
  extraOptsAnn = NULL,
  parallel = TRUE
)

\S4method{generateComponentsCliqueMS}{featureGroupsSet}(fGroups, ionization = NULL, ...)
}
\arguments{
\item{fGroups}{\code{\link{featureGroups}} object for which components should be generated.}

\item{ionization}{Which ionization polarity was used to generate the data: should be \code{"positive"}
  or \code{"negative"}. If the \code{featureGroups} object has adduct annotations, and \code{ionization=NULL}, the
  ionization will be detected automatically.

  \setsWF This parameter is not supported for sets workflows, as the ionization will always be detected
  automatically.}

\item{maxCharge, maxGrade, ppm}{Arguments passed to \code{\link[cliqueMS:getIsotopes]{cliqueMS::getIsotopes}} and/or
\code{\link[cliqueMS:getAnnotation]{cliqueMS::getAnnotation}}.}

\item{adductInfo}{Sets the \code{adinfo} argument to \code{\link[cliqueMS:getAnnotation]{cliqueMS::getAnnotation}}.
If \code{NULL} then the default adduct information from \pkg{cliqueMS} is used (\emph{i.e.} the
\code{positive.adinfo}/\code{negative.adinfo} package datasets).}

\item{absMzDev}{Maximum absolute \\emph{m/z} deviation.}

\item{minSize}{The minimum size of a component. Smaller components than this size will be removed. See note below.}

\item{relMinAdductAbundance}{The minimum relative abundance (\samp{0-1}) that an adduct should be assigned to
features within the same feature group. See the \verb{Feature components} section for more details.}

\item{adductConflictsUsePref}{If set to \code{TRUE}, and not all adduct assigments to the features within a feature
group are equal and at least one of those adducts is a preferential adduct (\code{prefAdducts} argument), then only
the features with (the lowest ranked) preferential adduct are considered. In all other cases or when
\code{adductConflictsUsePref=FALSE} only features with the most frequently assigned adduct is considered. See the
\verb{Feature components} section for more details.}

\item{NMConflicts}{The strategies to employ when not all neutral masses within a component are equal. Valid options
are: \code{"preferential"}, \code{"mostAbundant"} and \code{"mostIntense"}. Multiple strategies are possible, and
will be executed in the given order until one succeeds. See the \verb{Feature components} section for more details.}

\item{prefAdducts}{A \code{character} vector with one or more \emph{preferential adducts}. See the \verb{Feature
components} section for more details.}

\item{extraOptsCli, extraOptsIso, extraOptsAnn}{Named \code{list} with further arguments to be passed to
\code{\link[cliqueMS:getCliques]{cliqueMS::getCliques}}, \code{\link[cliqueMS:getIsotopes]{cliqueMS::getIsotopes}}
and \code{\link[cliqueMS:getAnnotation]{cliqueMS::getAnnotation}}, respectively. Set to \code{NULL} to ignore.}

\item{parallel}{If set to \code{TRUE} then code is executed in parallel through the \CRANpkg{futures} package. Please
see the parallelization section in the handbook for more details.}

\item{\dots}{\setsWF Further arguments passed to the non-sets workflow method.}
}
\value{
A \code{\link{componentsFeatures}} derived object.
}
\description{
Uses \href{https://github.com/osenan/cliqueMS}{cliqueMS} to generate components using the
\code{\link[cliqueMS:getCliques]{cliqueMS::getCliques}} function.
}
\details{
This function uses cliqueMS to generate components. This function is called when calling \code{generateComponents} with
  \code{algorithm="cliquems"}.

The grouping of features in each component ('clique') is based on high similarity of chromatographic elution
  profiles. All features in each component are then annotated with the
  \code{\link[cliqueMS:getIsotopes]{cliqueMS::getIsotopes}} and
  \code{\link[cliqueMS:getAnnotation]{cliqueMS::getAnnotation}} functions.
}
\section{Feature components}{
 The returned components are based on so called \emph{feature components}. Unlike other
  algorithms, components are first made on a feature level (per analysis), instead of for complete feature groups. In
  the final step the feature components are converted to 'regular' components by employing a consensus approach with
  the following steps:

  \enumerate{

  \item If an adduct assigned to a feature only occurs as a minority compared to other adduct assigments within the
  same feature group, it is considered as an outlier and removed accordingly (controlled by the
  \code{relMinAdductAbundance} argument).

  \item For features within a feature group, only keep their adduct assignment if it occurs as the most frequent or
  is preferential (controlled by \code{adductConflictsUsePref} and \code{prefAdducts} arguments).

  \item Components are made by combining the feature groups for which at least one of their features are jointly
  present in the same feature component.

  \item Conflicts of neutral mass assignments within a component (\emph{i.e.} not all are the same) are dealt with.
  Firstly, all feature groups with an unknown neutral mass are split in another component. Then, if conflicts still
  occur, the feature groups with similar neutral mass (determined by \code{absMzDev} argument) are grouped. Depending
  on the \code{NMConflicts} argument, the group with one or more preferential adduct(s) or that is the largest or
  most intense is selected, whereas others are removed from the component. In case multiple groups contain
  preferential adducts, and \samp{>1} preferential adducts are available, the group with the adduct that matches
  first in \code{prefAdducts} 'wins'. In case of ties, one of the next strategies in \code{NMConflicts} is tried.

  \item If a feature group occurs in multiple components it will be removed completely.

  \item the \code{minSize} filter is applied.

  }
}

\section{Sets workflows}{
 In a \link[=sets-workflow]{sets workflow} the componentization is first performed for each
  set independently. The resulting components are then all combined in a \code{\link{componentsSet}} object. Note that
  the components themselves are never merged. The components are renamed to include the set name from which they were
  generated (\emph{e.g.} \code{"CMP1"} becomes \code{"CMP1-positive"}).
}

\references{
\insertRef{Senan2019}{patRoon}
}
\seealso{
\code{\link{generateComponents}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-formulas.R
\name{formulaScorings}
\alias{formulaScorings}
\title{Scorings terms for formula candidates}
\usage{
formulaScorings()
}
\description{
Returns a \code{data.frame} with information on which scoring terms are used and what their algorithm specific name
is.
}
\seealso{
generateFormulas
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features-xcms3.R
\name{findFeaturesXCMS3}
\alias{findFeaturesXCMS3}
\title{Find features using XCMS (new interface)}
\usage{
findFeaturesXCMS3(
  analysisInfo,
  param = xcms::CentWaveParam(),
  ...,
  verbose = TRUE
)
}
\arguments{
\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}

\item{param}{The method parameters used by XCMS peak finding, see
\code{\link[xcms:findChromPeaks]{xcms::findChromPeaks}}}

\item{\dots}{Further parameters passed to \code{\link[xcms:findChromPeaks]{xcms::findChromPeaks}}.}

\item{verbose}{If set to \code{FALSE} then no text output is shown.}
}
\value{
An object of a class which is derived from \code{\link{features}}.
}
\description{
Uses the new \code{xcms3} interface from the \pkg{xcms} package to find features.
}
\details{
This function uses XCMS3 to automatically find features. This function is called when calling \code{findFeatures} with
  \code{algorithm="xcms3"}.

The file format of analyses must be \code{mzML} or \code{mzXML}.

The input MS data files need to be centroided. The \code{\link{convertMSFiles}} function can be used to
  centroid data.
}
\references{
\addCitations{xcms}{1} \cr\cr \addCitations{xcms}{2} \cr\cr \addCitations{xcms}{3}
}
\seealso{
\code{\link{findFeatures}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TP-biotransformer.R
\name{generateTPsBioTransformer}
\alias{generateTPsBioTransformer}
\title{Obtain transformation products (TPs) with BioTransformer}
\usage{
generateTPsBioTransformer(
  parents,
  type = "env",
  steps = 2,
  extraOpts = NULL,
  skipInvalid = TRUE,
  fpType = "extended",
  fpSimMethod = "tanimoto",
  MP = FALSE
)
}
\arguments{
\item{parents}{The parents for which transformation products should be obtained. This can be (1) a suspect list (see
\link[=suspect-screening]{suspect screening} for more information), (2) the resulting output
\code{\link{screenSuspects}} or (3) a \code{\link{compounds}} annotation object. In the former two cases, the
suspect (hits) are used as parents, whereas in the latter case all candidates are used as parents.}

\item{type}{The type of prediction. Valid values are: \code{"env"}, \code{"ecbased"}, \code{"cyp450"},
\code{"phaseII"}, \code{"hgut"}, \code{"superbio"}, \code{"allHuman"}. Sets the \command{-b} command line option.}

\item{steps}{The number of steps for the predictions. Sets the \command{-s} command line option.}

\item{extraOpts}{A \code{character} with extra command line options passed to the \command{biotransformer.jar} tool.}

\item{skipInvalid}{If set to \code{TRUE} then the parents will be skipped (with a warning) for which insufficient
information (\emph{e.g.} SMILES) is available.}

\item{fpType}{The type of structural fingerprint that should be calculated. See the \code{type} argument of the
\code{\link{get.fingerprint}} function of \CRANpkg{rcdk}.}

\item{fpSimMethod}{The method for calculating similarities (i.e. not dissimilarity!). See the \code{method} argument
of the \code{\link{fp.sim.matrix}} function of the \CRANpkg{fingerprint} package.}

\item{MP}{If \code{TRUE} then multiprocessing is enabled. Since \command{BioTransformer} supports native
parallelization, additional multiprocessing generally doesn't lead to significant reduction in computational times.
Furthermore, enabling multiprocessing can lead to very high CPU/RAM usage.}
}
\value{
The TPs are stored in an object from the \code{\link{transformationProductsBT}} class.
}
\description{
Uses \href{http://biotransformer.ca/}{BioTransformer} to predict TPs
}
\details{
This function uses BioTransformer to obtain transformation products. This function is called when calling \code{generateTPs} with
  \code{algorithm="biotransformer"}.

Structural similarities between the parent and its TPs are calculated, which can be used to
  \link[=filter,transformationProductsBT-method]{filter} the results.

  In order to use this function the \file{.jar} command line utility should be installed and specified in the
  \code{\link[=patRoon-package]{patRoon.path.BioTransformer}} option. The \file{.jar} file can be obtained via
  \url{https://bitbucket.org/djoumbou/biotransformer/src/master}.

An important advantage of this algorithm is that it provides structural information for generated TPs.
  However, this also means that if the input is from a parent suspect list or screening then either \acronym{SMILES}
  or \acronym{InChI} information must be available for the parents.
}
\note{
When the \code{parents} argument is a \code{\link{compounds}} object, the candidate library \code{identifier}
  is used in case the candidate has no defined \code{compoundName}.
}
\section{Parallelization}{
 \code{generateTPsBioTransformer} uses multiprocessing to parallelize
  computations. Please see the parallelization section in the handbook for
  more details and \link[=patRoon-package]{patRoon options} for configuration
  options.
}

\references{
\insertRef{DjoumbouFeunang2019}{patRoon} \cr\cr \insertRef{Wicker2015}{patRoon} \cr\cr
  \addCitations{rcdk}{1}
}
\seealso{
\code{\link{generateTPs}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mspeaklists-bruker.R
\name{generateMSPeakListsDA}
\alias{generateMSPeakListsDA}
\alias{generateMSPeakListsDA,featureGroups-method}
\alias{generateMSPeakListsDA,featureGroupsSet-method}
\title{Generate peak lists with Bruker DataAnalysis}
\usage{
\S4method{generateMSPeakListsDA}{featureGroups}(
  fGroups,
  bgsubtr = TRUE,
  maxMSRtWindow = 5,
  minMSIntensity = 500,
  minMSMSIntensity = 500,
  clear = TRUE,
  close = TRUE,
  save = close,
  MSMSType = "MSMS",
  avgFGroupParams = getDefAvgPListParams()
)

\S4method{generateMSPeakListsDA}{featureGroupsSet}(fGroups, ...)
}
\arguments{
\item{fGroups}{The \code{\link{featureGroups}} object from which MS peak lists should be extracted.}

\item{bgsubtr}{If \code{TRUE} background will be subtracted using the 'spectral' algorithm.}

\item{maxMSRtWindow}{Maximum chromatographic peak window used for spectrum averaging (in seconds, +/- retention
time). If \code{NULL} all spectra from a feature will be taken into account. Lower to decrease processing time.}

\item{minMSIntensity, minMSMSIntensity}{Minimum intensity for peak lists obtained with DataAnalysis. Highly
recommended to set \samp{>0} as DA tends to report many very low intensity peaks.}

\item{clear}{Remove any existing chromatogram traces/mass spectra prior to making new ones.}

\item{close, save}{If \code{TRUE} then Bruker files are closed and saved after
processing with DataAnalysis, respectively. Setting \code{close=TRUE}
prevents that many analyses might be opened simultaneously in DataAnalysis,
which otherwise may use excessive memory or become slow. By default
\code{save} is \code{TRUE} when \code{close} is \code{TRUE}, which is
likely what you want as otherwise any processed data is lost.}

\item{MSMSType}{The type of MS/MS experiment performed: \code{"MSMS"} for MRM/AutoMSMS or \code{"BBCID"} for
broadband CID.}

\item{avgFGroupParams}{A \code{list} with parameters used for averaging of peak lists for feature groups. See
\code{\link{getDefAvgPListParams}} for more details.}

\item{\dots}{\setsWF Further arguments passed to the non-sets workflow method.}
}
\value{
A \code{\link{MSPeakLists}} object.
}
\description{
Uses Bruker DataAnalysis to read the data needed to generate MS peak lists.
}
\details{
This function uses Bruker DataAnalysis to generate MS peak lists. This function is called when calling \code{generateMSPeakLists} with
  \code{algorithm="bruker"}.

The MS data should be in the Bruker data format (\file{.d}). This function leverages DataAnalysis
  functionality to support averaging of spectra, background subtraction and identification of isotopes. In order to
  obtain mass spectra TICs will be added in DataAnalysis of the MS and relevant MS/MS signals.
}
\note{
The \option{Component} column should be active (Method-->Parameters-->Layouts-->Mass List Layout) in order to
  add isotopologue information.

If any errors related to \command{DCOM} appear it might be necessary to
  terminate DataAnalysis (note that DataAnalysis might still be running as a
  background process). The \command{ProcessCleaner} application installed
  with DataAnalayis can be used for this.
}
\seealso{
\code{\link{generateMSPeakLists}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/components-specclust.R
\name{generateComponentsSpecClust}
\alias{generateComponentsSpecClust}
\alias{generateComponentsSpecClust,featureGroups-method}
\title{Generate components based on MS/MS similarity}
\usage{
\S4method{generateComponentsSpecClust}{featureGroups}(
  fGroups,
  MSPeakLists,
  method = "complete",
  specSimParams = getDefSpecSimParams(),
  maxTreeHeight = 1,
  deepSplit = TRUE,
  minModuleSize = 1
)
}
\arguments{
\item{fGroups}{\code{\link{featureGroups}} object for which components should be generated.}

\item{MSPeakLists}{The \code{\link{MSPeakLists}} object for the given feature groups that should be used for MS
spectral similarity calculations.}

\item{method}{Clustering method that should be applied (passed to
\code{\link[fastcluster:hclust]{fastcluster::hclust}}).}

\item{specSimParams}{A named \code{list} with parameters that influence the calculation of MS spectra similarities.
See the \link[=specSimParams]{spectral similarity parameters} documentation for more details.}

\item{maxTreeHeight, deepSplit, minModuleSize}{Arguments used by
\code{\link{cutreeDynamicTree}}.}
}
\value{
The components are stored in objects derived from \code{\link{componentsSpecClust}}.
}
\description{
Generates components based on MS/MS similarity between feature groups.
}
\details{
This function uses hierarchical clustering of MS/MS spectra to generate components. This function is called when calling \code{generateComponents} with
  \code{algorithm="specclust"}.

The similarities are converted to a distance matrix and used as input for hierarchical clustering, and the
  resulting dendrogram is automatically cut with \code{\link{cutreeDynamicTree}}. The clustering is performed with
  \code{\link[fastcluster:hclust]{fastcluster::hclust}}.
}
\section{Sets workflows}{
 In a \link[=sets-workflow]{sets workflow} the spectral similarities for each set are
  combined as is described for the \code{\link[=spectrumSimilarity,MSPeakListsSet-method]{spectrumSimilarity}} method
  for sets workflows.
}

\references{
\addCitations{fastcluster}{1}
}
\seealso{
\code{\link{generateComponents}} for more details and other algorithms.
}
\author{
Rick Helmus <\email{r.helmus@uva.nl}> and Bas van de Velde (major contributions to spectral binning and
  similarity calculation).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formulas.R, R/formulas-set.R
\docType{class}
\name{formulas-class}
\alias{formulas-class}
\alias{formulas}
\alias{formulasConsensus-class}
\alias{formulasConsensus}
\alias{annotations,formulas-method}
\alias{analyses,formulas-method}
\alias{defaultExclNormScores,formulas-method}
\alias{show,formulas-method}
\alias{[[,formulas,ANY,ANY-method}
\alias{delete,formulas-method}
\alias{as.data.table,formulas-method}
\alias{annotatedPeakList,formulas-method}
\alias{plotSpectrum,formulas-method}
\alias{plotScores,formulas-method}
\alias{consensus,formulas-method}
\alias{formulasSet-class}
\alias{formulasSet}
\alias{formulasConsensusSet-class}
\alias{formulasConsensusSet}
\alias{show,formulasSet-method}
\alias{delete,formulasSet-method}
\alias{[,formulasSet,ANY,missing,missing-method}
\alias{filter,formulasSet-method}
\alias{plotSpectrum,formulasSet-method}
\alias{annotatedPeakList,formulasSet-method}
\alias{consensus,formulasSet-method}
\alias{formulasUnset-class}
\alias{formulasUnset}
\alias{unset,formulasSet-method}
\alias{unset,formulasConsensusSet-method}
\title{Formula annotations class}
\usage{
\S4method{annotations}{formulas}(obj, features = FALSE)

\S4method{analyses}{formulas}(obj)

\S4method{defaultExclNormScores}{formulas}(obj)

\S4method{show}{formulas}(object)

\S4method{[[}{formulas,ANY,ANY}(x, i, j)

\S4method{delete}{formulas}(obj, i = NULL, j = NULL, ...)

\S4method{as.data.table}{formulas}(
  x,
  fGroups = NULL,
  fragments = FALSE,
  countElements = NULL,
  countFragElements = NULL,
  OM = FALSE,
  normalizeScores = "none",
  excludeNormScores = defaultExclNormScores(x),
  average = FALSE
)

\S4method{annotatedPeakList}{formulas}(
  obj,
  index,
  groupName,
  analysis = NULL,
  MSPeakLists,
  onlyAnnotated = FALSE
)

\S4method{plotSpectrum}{formulas}(
  obj,
  index,
  groupName,
  analysis = NULL,
  MSPeakLists,
  title = NULL,
  specSimParams = getDefSpecSimParams(),
  mincex = 0.9,
  xlim = NULL,
  ylim = NULL,
  ...
)

\S4method{plotScores}{formulas}(
  obj,
  index,
  groupName,
  analysis = NULL,
  normalizeScores = "max",
  excludeNormScores = defaultExclNormScores(obj)
)

\S4method{consensus}{formulas}(
  obj,
  ...,
  absMinAbundance = NULL,
  relMinAbundance = NULL,
  uniqueFrom = NULL,
  uniqueOuter = FALSE,
  rankWeights = 1,
  labels = NULL
)

\S4method{show}{formulasSet}(object)

\S4method{delete}{formulasSet}(obj, i, j, ...)

\S4method{[}{formulasSet,ANY,missing,missing}(x, i, j, ..., sets = NULL, updateConsensus = FALSE, drop = TRUE)

\S4method{filter}{formulasSet}(obj, ..., sets = NULL, updateConsensus = FALSE, negate = FALSE)

\S4method{plotSpectrum}{formulasSet}(
  obj,
  index,
  groupName,
  analysis = NULL,
  MSPeakLists,
  title = NULL,
  specSimParams = getDefSpecSimParams(),
  mincex = 0.9,
  xlim = NULL,
  ylim = NULL,
  perSet = TRUE,
  mirror = TRUE,
  ...
)

\S4method{annotatedPeakList}{formulasSet}(obj, index, groupName, analysis = NULL, MSPeakLists, ...)

\S4method{consensus}{formulasSet}(
  obj,
  ...,
  absMinAbundance = NULL,
  relMinAbundance = NULL,
  uniqueFrom = NULL,
  uniqueOuter = FALSE,
  rankWeights = 1,
  labels = NULL,
  filterSets = FALSE,
  setThreshold = 0,
  setThresholdAnn = 0
)

\S4method{unset}{formulasSet}(obj, set)

\S4method{unset}{formulasConsensusSet}(obj, set)
}
\arguments{
\item{obj, x, object}{The \code{formulas} object.}

\item{features}{If \code{TRUE} returns formula data for features, otherwise
for feature groups.}

\item{i, j}{For \code{[[}: If both \code{i} and \code{j} are specified then \code{i} specifies the analysis and
  \code{j} the feature group of the feature for which annotations should be returned. Otherwise \code{i} specifies
  the feature group for which group annotations should be returned. \code{i}/\code{j} can be specified as
  \code{integer} index or as a \code{character} name.

  Otherwise passed to the \code{\link[=filter,featureAnnotations-method]{featureAnnotations}} method.}

\item{\dots}{For \code{plotSpectrum}: Further arguments passed to \code{\link[graphics]{plot}}.

  For \code{delete}: passed to the function specified as \code{j}.

  For \code{consensus}: Any further (and unique) \code{formulas} objects.

  \setsPassedArgs1{formulas}}

\item{fGroups, fragments, countElements, countFragElements, OM}{Passed to the
\code{\link[=as.data.table,featureAnnotations-method]{featureAnnotations}} method.}

\item{normalizeScores}{A \code{character} that specifies how normalization of
annotation scorings occurs. Either
\code{"none"} (no normalization),
\code{"max"} (normalize to max value) or \code{"minmax"} (perform min-max
normalization). Note that normalization of negative scores (e.g. output by
\command{SIRIUS}) is always performed as min-max. Furthermore, currently
normalization for \code{compounds} takes the original min/max scoring
values into account when candidates were generated. Thus, for
\code{compounds} scoring, normalization is not affected when candidate
results were removed after they were generated (\emph{e.g.} by use of
\code{filter}).}

\item{excludeNormScores}{A
  \code{character} vector specifying any compound scoring names that
  should \emph{not} be normalized. Set to \code{NULL} to normalize all
  scorings. Note that whether any normalization occurs is set by the
  \code{excludeNormScores} argument.

  For \code{compounds}: By default \code{score} and
  \code{individualMoNAScore} are set to mimic the behavior of the
  \command{MetFrag} web interface.}

\item{average}{If set to \code{TRUE} an 'average formula' is generated for each feature group by combining all
elements from all candidates and averaging their amounts. This obviously leads to non-existing formulae, however,
this data may be useful to deal with multiple candidate formulae per feature group when performing elemental
characterization. Setting this to \code{TRUE} disables reporting of most other data.}

\item{index}{The candidate index (row). For \code{plotSpectrum} two indices can be specified to compare spectra. In
this case \code{groupName} and \code{analysis} (if not \code{NULL}) should specify values for the spectra to
compare.}

\item{groupName}{The name of the feature group (or feature groups when comparing spectra) to which the candidate
belongs.}

\item{analysis}{A \code{character} specifying the analysis (or analyses when comparing spectra) for which the
annotated spectrum should be plotted. If \code{NULL} then annotation results for the complete feature group will be
plotted.}

\item{MSPeakLists}{The \code{\link{MSPeakLists}} object that was used to generate the candidate}

\item{onlyAnnotated}{Set to \code{TRUE} to filter out any peaks that could
not be annotated.}

\item{title}{The title of the plot. Set to \code{NULL} for an automatically generated title.}

\item{specSimParams}{A named \code{list} with parameters that influence the calculation of MS spectra similarities.
See the \link[=specSimParams]{spectral similarity parameters} documentation for more details.}

\item{mincex}{The formula annotation labels are automatically scaled. The \code{mincex} argument forces a minimum
\code{cex} value for readability.}

\item{xlim, ylim}{Sets the plot size limits used by
\code{\link[graphics]{plot}}. Set to \code{NULL} for automatic plot sizing.}

\item{absMinAbundance, relMinAbundance}{Minimum absolute or relative
(\samp{0-1}) abundance across objects for a result to be kept. For
instance, \code{relMinAbundance=0.5} means that a result should be present
in at least half of the number of compared objects. Set to \samp{NULL} to
ignore and keep all results. Limits cannot be set when \code{uniqueFrom} is
not \code{NULL}.}

\item{uniqueFrom}{Set this argument to only retain formulas that are unique
within one or more of the objects for which the consensus is made.
Selection is done by setting the value of \code{uniqueFrom} to a
\code{logical} (values are recycled), \code{numeric} (select by index) or a
\code{character} (as obtained with \code{algorithm(obj)}). For
\code{logical} and \code{numeric} values the order corresponds to the order
of the objects given for the consensus. Set to \code{NULL} to ignore.}

\item{uniqueOuter}{If \code{uniqueFrom} is not \code{NULL} and if
\code{uniqueOuter=TRUE}: only retain data that are also unique between
objects specified in \code{uniqueFrom}.}

\item{rankWeights}{A numeric vector with weights of to calculate the mean
ranking score for each candidate. The value will be re-cycled if necessary,
hence, the default value of \samp{1} means equal weights for all considered
objects.}

\item{labels}{A \code{character} with names to use for labelling. If \code{NULL} labels are automatically generated.}

\item{sets}{\setsWF A \code{character} with name(s) of the sets to keep (or remove if \code{negate=TRUE}). Note: if
\code{updateConsensus=FALSE} then the \code{setCoverage} column of the annotation results is not updated.}

\item{updateConsensus}{\setsWF If \code{TRUE} then the annonation consensus among set results is updated. See the
\verb{Sets workflows} section for more details.}

\item{drop}{Passed to the \code{\link[=filter,featureAnnotations-method]{featureAnnotations}} method.}

\item{negate}{Passed to the \code{\link[=filter,featureAnnotations-method]{featureAnnotations}} method.}

\item{perSet, mirror}{\setsWF If \code{perSet=TRUE} then the set specific mass peaks are annotated separately.
Furthermore, if \code{mirror=TRUE} (and there are two sets in the object) then a mirror plot is generated.}

\item{filterSets}{\setsWF Controls how algorithms concensus abundance filters are applied. See the \verb{Sets
workflows} section below.}

\item{setThreshold, setThresholdAnn}{\setsWF Thresholds used to create the annotation set consensus. See
\code{\link{generateFormulas}}.}

\item{set}{\setsWF The name of the set.}
}
\value{
\code{annotations} returns a \code{list} containing for each feature
  group (or feature if \code{features=TRUE}) a \code{\link{data.table}}
  with an overview of all generated formulae and other data such as candidate
  scoring and MS/MS fragments.

\code{consensus} returns a \code{formulas} object that is produced by
  merging results from multiple \code{formulas} objects.
}
\description{
Contains data of generated chemical formulae for given feature groups.
}
\details{
\code{formulas} objects are obtained with \code{\link{generateFormulas}}. This class is derived from the
\code{\link{featureAnnotations}} class, please see its documentation for more methods and other details.
}
\section{Methods (by generic)}{
\itemize{
\item \code{annotations}: Accessor method to obtain generated formulae.

\item \code{analyses}: returns a \code{character} vector with the names of the
analyses for which data is present in this object.

\item \code{defaultExclNormScores}: Returns default scorings that are excluded from normalization.

\item \code{show}: Show summary information for this object.

\item \code{[[}: Extracts a formula table, either for a feature group or for features in an analysis.

\item \code{as.data.table}: Generates a table with all candidate formulae for each feature group and other information such
as element counts.

\item \code{annotatedPeakList}: Returns an MS/MS peak list annotated with data from a
given candidate formula.

\item \code{plotSpectrum}: Plots an annotated spectrum for a given candidate formula of a feature or feature group. Two
spectra can be compared by specifying a two-sized vector for the \code{index}, \code{groupName} and (if desired)
\code{analysis} arguments.

\item \code{plotScores}: Plots a barplot with scoring of a candidate formula.

\item \code{consensus}: Generates a consensus of results from multiple
objects. In order to rank the consensus candidates, first
each of the candidates are scored based on their original ranking
(the scores are normalized and the highest ranked candidate gets value
\samp{1}). The (weighted) mean is then calculated for all scorings of each
candidate to derive the final ranking (if an object lacks the candidate its
score will be \samp{0}). The original rankings for each object is stored in
the \code{rank} columns.
}}

\section{Slots}{

\describe{
\item{\code{featureFormulas}}{A \code{list} with all generated formulae for each analysis/feature group. Use the
\code{annotations} method for access.}

\item{\code{setThreshold,setThresholdAnn}}{\setsWF A copy of the equally named arguments that were passed when this object
was created by \code{\link{generateFormulas}}.}

\item{\code{origFGNames}}{\setsWF The original (order of) names of the \code{\link{featureGroups}} object that was used to
create this object.}
}}

\section{S4 class hierarchy}{
 \itemize{   \item{\code{\link{featureAnnotations}}}   \itemize{     \item{\strong{\code{\link{formulas}}}}     \itemize{       \item{\code{\link{formulasConsensus}}}       \item{\code{\link{formulasSet}}}       \itemize{         \item{\code{\link{formulasConsensusSet}}}       }       \item{\code{\link{formulasUnset}}}     }   } }
}

\section{Source}{
 Subscripting of formulae for plots generated by
  \code{plotSpectrum} is based on the \code{chemistry2expression} function
  from the \href{https://github.com/schymane/ReSOLUTION}{ReSOLUTION} package.
}

\section{Sets workflows}{
 \setsWFClass{formulasSet}{formulas}

  \setsWFNewMethodsSO{formulasUnset}{Only the annotation results that are present in the specified set are kept
  (based on the set consensus, see below for implications).}

  \setsWFChangedMethods{

  \item \code{filter} and the subset operator (\code{[}) Can be used to select data that is only present for selected
  sets. Depending on the \code{updateConsenus}, both either operate on set consensus or original data (see below for
  implications).

  \item \code{annotatedPeakList} Returns a combined annotation table with all sets.

  \item \code{plotSpectrum} Is able to highlight set specific mass peaks (\code{perSet} and \code{mirror} arguments).

  \item \code{consensus} Creates the algorithm consensus based on the original annotation data (see below for
  implications). Then, like the sets workflow method for \code{\link{generateFormulas}}, a consensus is made for all
  sets, which can be controlled with the \code{setThreshold} and \code{setThresholdAnn} arguments. The candidate
  coverage among the different algorithms is calculated for each set (\emph{e.g.} \code{coverage-positive} column)
  and for all sets (\code{coverage} column), which is based on the presence of a candidate in all the algorithms from
  all sets data. The \code{consensus} method for sets workflow data supports the \code{filterSets} argument. This
  controls how the algorithm consensus abundance filters (\code{absMinAbundance}/\code{relMinAbundance}) are applied:
  if \code{filterSets=TRUE} then the minimum of all \code{coverage} set specific columns is used to obtain the
  algorithm abundance. Otherwise the overall \code{coverage} column is used. For instance, consider a consensus
  object to be generated from two objects generated by different algorithms (\emph{e.g.} \command{SIRIUS} and
  \command{GenForm}), which both have a positive and negative set. Then, if a candidate occurs with both
  algorithms for the positive mode set, but only with the first algorithm in the negative mode set,
  \code{relMinAbundance=1} will remove the candidate if \code{filterSets=TRUE} (because the minimum relative
  algorithm abundance is \samp{0.5}), while \code{filterSets=FALSE} will not remove the candidate (because based on
  all sets data the candidate occurs in both algorithms).

  

  }

  Two types of annotation data are stored in a \code{formulasSet} object: \enumerate{

  \item Annotations that are produced from a consensus between set results (see \code{generateFormulas}).

  \item The 'original' annotation data per set, prior to when the set consensus was made. This includes candidates
  that were filtered out because of the thresholds set by \code{setThreshold} and \code{setThresholdAnn}. However,
  when \code{filter} or subsetting (\code{[}) operations are performed, the original data is also updated.

  }

  In most cases the first data is used. However, in a few cases the original annotation data is used (as indicated
  above), for instance, to re-create the set consensus. It is important to realize that the original annotation data
  may have \emph{additional} candidates, and a newly created set consensus may therefore have 'new' candidates. For
  instance, when the object consists of the sets \code{"positive"} and \code{"negative"} and \code{setThreshold=1}
  was used to create it, then \code{formulas[, sets = "positive", updateConsensus = TRUE]} may now have additional
  candidates, \emph{i.e.} those that were not present in the \code{"negative"} set and were previously removed due to
  the consensus threshold filter.
}

\seealso{
The \code{\link{featureAnnotations}} base class for more relevant methods and
  \code{\link{generateFormulas}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compounds-sirius.R
\name{generateCompoundsSIRIUS}
\alias{generateCompoundsSIRIUS}
\alias{generateCompoundsSIRIUS,featureGroups-method}
\alias{generateCompoundsSIRIUS,featureGroupsSet-method}
\title{Compound annotation with SIRIUS}
\usage{
\S4method{generateCompoundsSIRIUS}{featureGroups}(
  fGroups,
  MSPeakLists,
  relMzDev = 5,
  adduct = NULL,
  elements = "CHNOP",
  profile = "qtof",
  formulaDatabase = NULL,
  fingerIDDatabase = "pubchem",
  noise = NULL,
  errorRetries = 2,
  cores = NULL,
  topMost = 100,
  topMostFormulas = 5,
  extraOptsGeneral = NULL,
  extraOptsFormula = NULL,
  verbose = TRUE,
  splitBatches = FALSE
)

\S4method{generateCompoundsSIRIUS}{featureGroupsSet}(
  fGroups,
  MSPeakLists,
  relMzDev = 5,
  adduct = NULL,
  ...,
  setThreshold = 0,
  setThresholdAnn = 0
)
}
\arguments{
\item{fGroups}{\code{\link{featureGroups}} object which should be annotated. This should be the same or a subset of
the object that was used to create the specified \code{MSPeakLists}. In the case of a subset only the remaining
feature groups in the subset are considered.}

\item{MSPeakLists}{A \code{\link{MSPeakLists}} object that was generated for the supplied \code{fGroups}.}

\item{relMzDev}{Maximum relative deviation between the measured and candidate formula \emph{m/z} values (in ppm).
Sets the \option{--ppm-max} command line option.}

\item{adduct}{An \code{\link{adduct}} object (or something that can be converted to it with \code{\link{as.adduct}}).
  Examples: \code{"[M-H]-"}, \code{"[M+Na]+"}. If the \code{featureGroups} object has
  adduct annotations then these are used if \code{adducts=NULL}.

  \setsWF The \code{adduct} argument is not supported for sets workflows, since the
  adduct annotations will then always be used.}

\item{elements}{Elements to be considered for formulae calculation. This will heavily affects the number of
candidates! Always try to work with a minimal set by excluding elements you don't expect. The minimum/maximum
number of elements can also be specified, for example: a value of \code{"C[5]H[10-15]O"} will only consider
formulae with up to five carbon atoms, between ten and fifteen hydrogen atoms and any amount of oxygen atoms. Sets
the \option{--elements} command line  option.}

\item{profile}{Name of the configuration profile, for example: \option{"qtof"}, \option{"orbitrap"},
\option{"fticr"}. Sets the \option{--profile} commandline option.}

\item{formulaDatabase}{If not \code{NULL}, use a database for retrieval of formula
candidates. Possible values are: \option{"pubchem"}, \option{"bio"}, \option{"kegg"}, \option{"hmdb"}. Sets the
\option{--database} commandline option.}

\item{fingerIDDatabase}{Database specifically used for \command{CSI:FingerID}. If \code{NULL}, the value of the
\code{formulaDatabase} parameter will be used or \code{"pubchem"} when that is also \code{NULL}. Sets the
\option{--fingerid-db} option.}

\item{noise}{Median intensity of the noise (\code{NULL} ignores this parameter). Sets the \option{--noise}
commandline option.}

\item{errorRetries}{Maximum number of retries after an error occurred. This may be useful to handle e.g. connection
errors.}

\item{cores}{The number of cores \command{SIRIUS} will use. If \code{NULL} then the default of all cores will be
used.}

\item{topMost}{Only keep this number of candidates (per feature group) with highest score. Set to \code{NULL} to
always keep all candidates, however, please note that this may result in significant usage of CPU/RAM resources for
large numbers of candidates.}

\item{topMostFormulas}{Do not return more than this number of candidate formulae. Note that only compounds for these
formulae will be searched. Sets the \option{--candidates} commandline option.}

\item{extraOptsGeneral, extraOptsFormula}{a \code{character} vector with any extra commandline parameters for
\command{SIRIUS}. For \command{SIRIUS} versions \code{<4.4} there is no distinction between general and formula
options. Otherwise commandline options specified in \code{extraOptsGeneral} are added prior to the \code{formula}
command, while options specified in \code{extraOptsFormula} are added in afterwards. See the \command{SIRIUS}
manual for more details. Set to \code{NULL} to ignore.}

\item{verbose}{If \code{TRUE} then more output is shown in the terminal.}

\item{splitBatches}{If \code{TRUE} then the calculations done by \command{SIRIUS} will be evenly split over multiple
\command{SIRIUS} calls (which may be run in parallel depending on the \link[=patRoon-package]{set package
options}). If \code{splitBatches=FALSE} then all feature calculations are performed from a single \command{SIRIUS}
execution, which is often the fastest if calculations are performed on a single computer.}

\item{\dots}{\setsWF Further arguments passed to the non-sets workflow method.}

\item{setThreshold}{\setsWF Minimum abundance for a candidate among all sets (\samp{0-1}). For instance, a value of
\samp{1} means that the candidate needs to be present in all the set data.}

\item{setThresholdAnn}{\setsWF As \code{setThreshold}, but only taking into account the set data that contain
annotations for the feature group of the candidate.}
}
\value{
A \code{\link{compounds}} derived object containing all compound annotations.
}
\description{
Uses \href{https://bio.informatik.uni-jena.de/software/sirius/}{SIRIUS} in combination with
\href{https://www.csi-fingerid.uni-jena.de/}{CSI:FingerID} for compound annotation.
}
\details{
This function uses SIRIUS to generate compound candidates. This function is called when calling \code{generateCompounds} with
  \code{algorithm="sirius"}.

Similar to \code{\link{generateFormulasSIRIUS}}, candidate formulae are generated with SIRIUS. These results
  are then feed to CSI:FingerID to acquire candidate structures. This method requires the availability of MS/MS data,
  and feature groups without it will be ignored.
}
\note{
For annotations performed with \command{SIRIUS} it is often the fastest to keep the default
  \code{splitBatches=FALSE}. In this case, all \command{SIRIUS} output will be printed to the terminal (unless
  \code{verbose=FALSE} or \option{patRoon.MP.method="future"}). Furthermore, please note that only annotations to be
  performed for the same adduct are grouped in a single batch execution.
}
\section{Parallelization}{
 generateCompoundsSIRIUS uses multiprocessing to parallelize
  computations. Please see the parallelization section in the handbook for
  more details and \link[=patRoon-package]{patRoon options} for configuration
  options.
}

\references{
\insertRef{Dhrkop2019}{patRoon} \cr\cr \insertRef{Duhrkop2015}{patRoon} \cr\cr
  \insertRef{Duhrkop2015-2}{patRoon} \cr\cr \insertRef{Bcker2008}{patRoon}
}
\seealso{
\code{\link{generateCompounds}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features-kpic2.R
\name{importFeaturesKPIC2}
\alias{importFeaturesKPIC2}
\title{Imports features from KPIC2}
\usage{
importFeaturesKPIC2(picsList, analysisInfo)
}
\arguments{
\item{picsList}{A \code{list} with a \code{pics} objects obtained with \code{\link[KPIC]{getPIC}} or
\code{\link[KPIC]{getPIC.kmeans}} for each analysis.}

\item{analysisInfo}{A \code{data.frame} with \link[=analysis-information]{Analysis information}.}
}
\value{
An object of a class which is derived from \code{\link{features}}.
}
\description{
Imports feature data generated by the \pkg{KPIC2} package.
}
\details{
This function imports data from KPIC2. This function is called when calling \code{importFeatures} with
  \code{type="kpic2"}.
}
\references{
\insertRef{Ji2017}{patRoon}
}
\seealso{
\code{\link{importFeatures}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compounds-cluster.R
\name{compounds-cluster}
\alias{compounds-cluster}
\alias{makeHCluster,compounds-method}
\alias{makeHCluster}
\title{Hierarchical clustering of compounds}
\usage{
\S4method{makeHCluster}{compounds}(
  obj,
  method,
  fpType = "extended",
  fpSimMethod = "tanimoto",
  maxTreeHeight = 1,
  deepSplit = TRUE,
  minModuleSize = 1
)
}
\arguments{
\item{obj}{The \code{\link{compounds}} object to be clustered.}

\item{method}{The clustering method passed to \code{\link{hclust}}.}

\item{fpType}{The type of structural fingerprint that should be calculated. See the \code{type} argument of the
\code{\link{get.fingerprint}} function of \CRANpkg{rcdk}.}

\item{fpSimMethod}{The method for calculating similarities (i.e. not dissimilarity!). See the \code{method} argument
of the \code{\link{fp.sim.matrix}} function of the \CRANpkg{fingerprint} package.}

\item{maxTreeHeight, deepSplit, minModuleSize}{Arguments used by
\code{\link{cutreeDynamicTree}}.}
}
\value{
\code{makeHCluster} returns an \code{\link{compoundsCluster}} object.
}
\description{
Perform hierarchical clustering of structure candidates based on chemical
similarity and obtain overall structural information based on the maximum
common structure (MCS).
}
\details{
Often many possible chemical structure candidates are found for each feature
group when performing \link[=generateCompounds]{compound annotation}.
Therefore, it may be useful to obtain an overview of their general structural
properties. One strategy is to perform hierarchical clustering based on their
chemical (dis)similarity, for instance, using the Tanimoto score. The
resulting clusters can then be characterized by evaluating their
\emph{maximum common substructure} (MCS).

\code{makeHCluster} performs hierarchical clustering of all
  structure candidates for each feature group within a
  \code{\link{compounds}} object. The resulting dendrograms are automatically
  cut using the \code{\link{cutreeDynamicTree}} function from the
  \pkg{\link{dynamicTreeCut}} package. The returned
  \code{\link{compoundsCluster}} object can then be used, for instance, for
  plotting dendrograms and MCS structures and manually re-cutting specific
  clusters.
}
\section{Source}{
 The methodology applied here has been largely derived from
  \file{chemclust.R} from the \pkg{metfRag} package and the package vignette
  of \CRANpkg{rcdk}.
}

\references{
\addCitations{rcdk}{1}
}
\seealso{
compoundsCluster
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formulas-bruker.R
\name{generateFormulasDA}
\alias{generateFormulasDA}
\alias{generateFormulasDA,featureGroups-method}
\alias{generateFormulasDA,featureGroupsSet-method}
\title{Generate formula with Bruker DataAnalysis}
\usage{
\S4method{generateFormulasDA}{featureGroups}(
  fGroups,
  MSPeakLists,
  precursorMzSearchWindow = 0.002,
  MSMode = "both",
  adduct = NULL,
  featThreshold = 0,
  featThresholdAnn = 0.75,
  absAlignMzDev = 0.002,
  save = TRUE,
  close = save
)

\S4method{generateFormulasDA}{featureGroupsSet}(
  fGroups,
  MSPeakLists,
  precursorMzSearchWindow = 0.002,
  MSMode = "both",
  adduct = NULL,
  ...,
  setThreshold = 0,
  setThresholdAnn = 0
)
}
\arguments{
\item{fGroups}{\code{\link{featureGroups}} object for which formulae should be generated. This should be the same or
a subset of the object that was used to create the specified \code{MSPeakLists}. In the case of a subset only the
remaining feature groups in the subset are considered.}

\item{MSPeakLists}{An \code{\link{MSPeakLists}} object that was generated for the supplied \code{fGroups}.}

\item{precursorMzSearchWindow}{Search window for \emph{m/z} values (+/- the feature \emph{m/z}) used to find back
feature data of precursor/parent ions from MS/MS spectra (this data is not readily available from
\command{SmartFormula3D} results).}

\item{MSMode}{Whether formulae should be generated only from MS data (\code{"ms"}), MS/MS data (\code{"msms"}) or
both (\code{"both"}). Selecting "both" will calculate formulae from MS data and MS/MS data and combines the results
(duplicated formulae are removed). This is useful when poor MS/MS data would exclude proper candidates.}

\item{adduct}{An \code{\link{adduct}} object (or something that can be converted to it with \code{\link{as.adduct}}).
  Examples: \code{"[M-H]-"}, \code{"[M+Na]+"}. If the \code{featureGroups} object has
  adduct annotations then these are used if \code{adducts=NULL}.

  \setsWF The \code{adduct} argument is not supported for sets workflows, since the
  adduct annotations will then always be used.}

\item{featThreshold}{If \code{calculateFeatures=TRUE}: minimum presence (\samp{0-1}) of a formula in all features
before it is considered as a candidate for a feature group. For instance, \code{featThreshold=0.75} dictates that a
formula should be present in at least 75\% of the features inside a feature group.}

\item{featThresholdAnn}{As \code{featThreshold}, but only considers features with annotations. For instance,
\code{featThresholdAnn=0.75} dictates that a formula should be present in at least 75\% of the features with
annotations inside a feature group.}

\item{absAlignMzDev}{When the group formula annotation consensus is made from feature annotations, the \emph{m/z}
values of annotated MS/MS fragments may slightly deviate from those of the corresponding group MS/MS peak list. The
\code{absAlignMzDev} argument specifies the maximum \emph{m/z} window used to re-align the mass peaks.}

\item{close, save}{If \code{TRUE} then Bruker files are closed and saved after
processing with DataAnalysis, respectively. Setting \code{close=TRUE}
prevents that many analyses might be opened simultaneously in DataAnalysis,
which otherwise may use excessive memory or become slow. By default
\code{save} is \code{TRUE} when \code{close} is \code{TRUE}, which is
likely what you want as otherwise any processed data is lost.}

\item{\dots}{\setsWF Further arguments passed to the non-sets workflow method.}

\item{setThreshold}{\setsWF Minimum abundance for a candidate among all sets (\samp{0-1}). For instance, a value of
\samp{1} means that the candidate needs to be present in all the set data.}

\item{setThresholdAnn}{\setsWF As \code{setThreshold}, but only taking into account the set data that contain
annotations for the feature group of the candidate.}
}
\value{
A \code{\link{formulas}} object containing all generated formulae.
}
\description{
Uses Bruker DataAnalysis to generate chemical formulae.
}
\details{
This function uses bruker to generate formula candidates. This function is called when calling \code{generateFormulas} with
  \code{algorithm="bruker"}.

This method supports scoring based on overlap between measured and theoretical isotopic patterns (both MS
  and MS/MS data) and the presence of 'fitting' MS/MS fragments. The method will iterate through all features (or
  "Compounds" in DataAnalysis terms) and call \command{SmartFormula} (and \command{SmartFormula3D} if MS/MS data is
  available) to generate all formulae. Parameters affecting formula calculation have to be set in advance within the
  DataAnalysis method for each analysis (\emph{e.g.} by \code{\link{setDAMethod}}).

  This method requires that features were obtained with \code{\link{findFeaturesBruker}}. It is recommended, but not
  mandatory, that the \code{\link{MSPeakLists}} are also generated by DataAnalysis.

  Calculation of formulae with DataAnalysis always occurs with the 'feature approach' (see \verb{Candidate
  assignment} in \code{\link{generateFormulas}}).
}
\note{
If any errors related to \command{DCOM} appear it might be necessary to
  terminate DataAnalysis (note that DataAnalysis might still be running as a
  background process). The \command{ProcessCleaner} application installed
  with DataAnalayis can be used for this.
}
\seealso{
\code{\link{generateFormulas}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formulas-sirius.R
\name{generateFormulasSIRIUS}
\alias{generateFormulasSIRIUS}
\alias{generateFormulasSIRIUS,featureGroups-method}
\alias{generateFormulasSIRIUS,featureGroupsSet-method}
\title{Generate formula with SIRIUS}
\usage{
\S4method{generateFormulasSIRIUS}{featureGroups}(
  fGroups,
  MSPeakLists,
  relMzDev = 5,
  adduct = NULL,
  elements = "CHNOP",
  profile = "qtof",
  database = NULL,
  noise = NULL,
  cores = NULL,
  topMost = 100,
  extraOptsGeneral = NULL,
  extraOptsFormula = NULL,
  calculateFeatures = TRUE,
  featThreshold = 0,
  featThresholdAnn = 0.75,
  absAlignMzDev = 0.002,
  verbose = TRUE,
  splitBatches = FALSE
)

\S4method{generateFormulasSIRIUS}{featureGroupsSet}(
  fGroups,
  MSPeakLists,
  relMzDev = 5,
  adduct = NULL,
  ...,
  setThreshold = 0,
  setThresholdAnn = 0
)
}
\arguments{
\item{fGroups}{\code{\link{featureGroups}} object for which formulae should be generated. This should be the same or
a subset of the object that was used to create the specified \code{MSPeakLists}. In the case of a subset only the
remaining feature groups in the subset are considered.}

\item{MSPeakLists}{An \code{\link{MSPeakLists}} object that was generated for the supplied \code{fGroups}.}

\item{relMzDev}{Maximum relative deviation between the measured and candidate formula \emph{m/z} values (in ppm).
Sets the \option{--ppm-max} command line option.}

\item{adduct}{An \code{\link{adduct}} object (or something that can be converted to it with \code{\link{as.adduct}}).
  Examples: \code{"[M-H]-"}, \code{"[M+Na]+"}. If the \code{featureGroups} object has
  adduct annotations then these are used if \code{adducts=NULL}.

  \setsWF The \code{adduct} argument is not supported for sets workflows, since the
  adduct annotations will then always be used.}

\item{elements}{Elements to be considered for formulae calculation. This will heavily affects the number of
candidates! Always try to work with a minimal set by excluding elements you don't expect. The minimum/maximum
number of elements can also be specified, for example: a value of \code{"C[5]H[10-15]O"} will only consider
formulae with up to five carbon atoms, between ten and fifteen hydrogen atoms and any amount of oxygen atoms. Sets
the \option{--elements} command line  option.}

\item{profile}{Name of the configuration profile, for example: \option{"qtof"}, \option{"orbitrap"},
\option{"fticr"}. Sets the \option{--profile} commandline option.}

\item{database}{If not \code{NULL}, use a database for retrieval of formula
candidates. Possible values are: \option{"pubchem"}, \option{"bio"}, \option{"kegg"}, \option{"hmdb"}. Sets the
\option{--database} commandline option.}

\item{noise}{Median intensity of the noise (\code{NULL} ignores this parameter). Sets the \option{--noise}
commandline option.}

\item{cores}{The number of cores \command{SIRIUS} will use. If \code{NULL} then the default of all cores will be
used.}

\item{topMost}{Only keep this number of candidates (per feature group) with highest
score. Sets the \option{--candidates} command line option.}

\item{extraOptsGeneral, extraOptsFormula}{a \code{character} vector with any extra commandline parameters for
\command{SIRIUS}. For \command{SIRIUS} versions \code{<4.4} there is no distinction between general and formula
options. Otherwise commandline options specified in \code{extraOptsGeneral} are added prior to the \code{formula}
command, while options specified in \code{extraOptsFormula} are added in afterwards. See the \command{SIRIUS}
manual for more details. Set to \code{NULL} to ignore.}

\item{calculateFeatures}{If \code{TRUE} fomulae are first calculated for all features
prior to feature group assignment (see \verb{Candidate assignment} in \code{\link{generateFormulas}}).}

\item{featThreshold}{If \code{calculateFeatures=TRUE}: minimum presence (\samp{0-1}) of a formula in all features
before it is considered as a candidate for a feature group. For instance, \code{featThreshold=0.75} dictates that a
formula should be present in at least 75\% of the features inside a feature group.}

\item{featThresholdAnn}{As \code{featThreshold}, but only considers features with annotations. For instance,
\code{featThresholdAnn=0.75} dictates that a formula should be present in at least 75\% of the features with
annotations inside a feature group. @param topMost Only keep this number of candidates
(per feature group) with highest score. Sets the \option{--candidates} command line
option.}

\item{absAlignMzDev}{When the group formula annotation consensus is made from feature annotations, the \emph{m/z}
values of annotated MS/MS fragments may slightly deviate from those of the corresponding group MS/MS peak list. The
\code{absAlignMzDev} argument specifies the maximum \emph{m/z} window used to re-align the mass peaks.}

\item{verbose}{If \code{TRUE} then more output is shown in the terminal.}

\item{splitBatches}{If \code{TRUE} then the calculations done by \command{SIRIUS} will be evenly split over multiple
\command{SIRIUS} calls (which may be run in parallel depending on the \link[=patRoon-package]{set package
options}). If \code{splitBatches=FALSE} then all feature calculations are performed from a single \command{SIRIUS}
execution, which is often the fastest if calculations are performed on a single computer.}

\item{\dots}{\setsWF Further arguments passed to the non-sets workflow method.}

\item{setThreshold}{\setsWF Minimum abundance for a candidate among all sets (\samp{0-1}). For instance, a value of
\samp{1} means that the candidate needs to be present in all the set data.}

\item{setThresholdAnn}{\setsWF As \code{setThreshold}, but only taking into account the set data that contain
annotations for the feature group of the candidate.}
}
\value{
A \code{\link{formulas}} object containing all generated formulae.
}
\description{
Uses \href{https://bio.informatik.uni-jena.de/software/sirius/}{SIRIUS} to generate chemical formulae candidates.
}
\details{
This function uses sirius to generate formula candidates. This function is called when calling \code{generateFormulas} with
  \code{algorithm="sirius"}.

Similarity of measured and theoretical isotopic patterns will be used for scoring candidates. Note that
  \command{SIRIUS} requires availability of MS/MS data.
}
\note{
For annotations performed with \command{SIRIUS} it is often the fastest to keep the default
  \code{splitBatches=FALSE}. In this case, all \command{SIRIUS} output will be printed to the terminal (unless
  \code{verbose=FALSE} or \option{patRoon.MP.method="future"}). Furthermore, please note that only annotations to be
  performed for the same adduct are grouped in a single batch execution.
}
\section{Parallelization}{
 generateFormulasSIRIUS uses multiprocessing to parallelize
  computations. Please see the parallelization section in the handbook for
  more details and \link[=patRoon-package]{patRoon options} for configuration
  options.
}

\references{
\insertRef{Dhrkop2019}{patRoon} \cr\cr \insertRef{Duhrkop2015}{patRoon} \cr\cr
  \insertRef{Duhrkop2015-2}{patRoon} \cr\cr \insertRef{Bcker2008}{patRoon}
}
\seealso{
\code{\link{generateFormulas}} for more details and other algorithms.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/components-openms.R
\name{defaultOpenMSAdducts}
\alias{defaultOpenMSAdducts}
\title{Default adducts for OpenMS componentization}
\usage{
defaultOpenMSAdducts(ionization)
}
\arguments{
\item{ionization}{The ionization polarity: either \code{"positive"} or \code{"negative"}.}
}
\description{
Returns the default adducts and their probabilities when the OpenMS algorithm is used for componentization.
}
\details{
See the \code{potentialAdducts} argument of \code{\link{generateComponentsOpenMS}} for more details.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TP-biotransformer.R
\docType{class}
\name{transformationProductsBT-class}
\alias{transformationProductsBT-class}
\alias{transformationProductsBT}
\alias{convertToMFDB,transformationProductsBT-method}
\alias{filter,transformationProductsBT-method}
\title{Class to store transformation products (TPs) predicted by BioTransformer}
\usage{
\S4method{convertToMFDB}{transformationProductsBT}(TPs, out, includeParents = FALSE)

\S4method{filter}{transformationProductsBT}(
  obj,
  removeParentIsomers = FALSE,
  removeTPIsomers = FALSE,
  removeDuplicates = FALSE,
  minSimilarity = NULL,
  negate = FALSE
)
}
\arguments{
\item{out}{The file name of the the output \command{MetFrag} database.}

\item{includeParents}{Set to \code{TRUE} to include the parents in the database.}

\item{obj, TPs}{\code{transformationProductsBTs} object to be accessed}

\item{removeParentIsomers}{If \code{TRUE} then TPs with an equal formula as their parent (isomers) are removed.}

\item{removeTPIsomers}{If \code{TRUE} then all TPs with equal formula as any sibling TPs (isomers) are removed.
Unlike \code{removeDuplicates}, \emph{all} TP candidates are removed (including the first match). This filter
automatically sets \code{removeDuplicates=TRUE} to avoid complete removal of TPs with equal structure.}

\item{removeDuplicates}{If \code{TRUE} then the TPs of a parent with duplicate structures (\acronym{SMILES}) are
removed. Such duplicates may occur when different transformation pathways yield the same TPs. The first TP
candidate with duplicate structure will be kept.}

\item{minSimilarity}{Minimum structure similarity (\samp{0-1}) that a TP should have relative to its parent. For
details on how these similarities are calculated, see the \code{\link{generateTPsBioTransformer}} function. May be
useful under the assumption that parents and TPs who have a high structural similarity, also likely have a high
MS/MS spectral similarity (which can be evaluated after componentization with \code{\link{generateComponentsTPs}}.}

\item{negate}{If \code{TRUE} then filters are performed in opposite manner.}
}
\value{
\code{filter} returns a filtered \code{transformationProductsBT} object.
}
\description{
This class is used to store prediction results that are generated with
\href{http://biotransformer.ca/}{BioTransformer}.
}
\details{
Objects from this class are generate with \code{\link{generateTPsBioTransformer}}. This class is derived from the
\code{\link{transformationProducts}} base class, please see its documentation for more details.
}
\section{Methods (by generic)}{
\itemize{
\item \code{convertToMFDB}: Exports this object as a \file{.csv} file that can be used as a \command{MetFrag} local
database.

\item \code{filter}: Performs rule-based filtering of the \command{BioTransformer} predictions.
Useful to simplify and clean-up the data.
}}

\section{S4 class hierarchy}{
 \itemize{   \item{\code{\link{transformationProducts}}}   \itemize{     \item{\strong{\code{\link{transformationProductsBT}}}}   } }
}

\references{
\insertRef{DjoumbouFeunang2019}{patRoon} \cr\cr \insertRef{Wicker2015}{patRoon}
}
\seealso{
The base class \code{\link{transformationProducts}} for more relevant methods and
  \code{\link{generateTPs}}
}
\newcommand{\addCitations}{\Sexpr[results=text,echo=FALSE,strip.white=FALSE]{ret <- format(citation("#1")[[#2]], "textVersion"); if (requireNamespace("pkgdown", quietly = TRUE) && pkgdown::in_pkgdown()) ret <- gsub("<", "&lt;", ret, fixed = TRUE); return(ret)}}

\newcommand{\setsWF}{(\strong{sets workflow})}
\newcommand{\setsWFClass}{The \code{#1} class is applicable for \link[=sets-workflow]{sets workflows}. This class is derived from \code{#2} and therefore largely follows the same user interface.}
\newcommand{\setsWFUnset}{\item \code{unset} Converts the object data for a specified set into a 'non-set' object (\code{#1}), which allows it to be used in 'regular' workflows. #2}
\newcommand{\setsWFNewMethodsFeat}{The following methods are specifically defined for sets workflows: \itemize{ \item \code{sets} Returns the set names for this object. \setsWFUnset{#1}{#2} }}
\newcommand{\setsWFNewMethodsSO}{The following methods are specifically defined for sets workflows: \itemize{\item All the methods from base class \code{\link{workflowStepSet}}. \setsWFUnset{#1}{#2} }}
\newcommand{\setsWFNewMethodsSOExtra}{The following methods are specifically defined for sets workflows: \itemize{\item All the methods from base class \code{\link{workflowStepSet}}. \setsWFUnset{#1}{#2} #3 }}
\newcommand{\setsWFChangedMethods}{The following methods are changed or with new functionality: \itemize{#1}}
\newcommand{\setsPassedArgs}{For \link[=sets-workflow]{sets workflow} methods: further arguments passed to the base #1.}
\newcommand{\setsPassedArgs1}{\setsPassedArgs{\code{\link{#1}} method}}

