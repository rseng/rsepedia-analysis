[![DOI](https://joss.theoj.org/papers/10.21105/joss.03416/status.svg)](https://doi.org/10.21105/joss.03416)
[![arXiv:2109.04193](https://img.shields.io/badge/arXiv-2109.04193-b31b1b.svg)](https://arxiv.org/abs/2109.04193)
[![License: MIT](https://img.shields.io/github/license/bshoshany/thread-pool)](https://github.com/bshoshany/OGRe/blob/master/LICENSE.txt)
![Language: Mathematica 12](https://img.shields.io/badge/Language-Mathematica%2012-yellow)
![File size in bytes](https://img.shields.io/github/size/bshoshany/OGRe/OGRe.m)
![GitHub last commit](https://img.shields.io/github/last-commit/bshoshany/OGRe)
[![GitHub repo stars](https://img.shields.io/github/stars/bshoshany/OGRe?style=social)](https://github.com/bshoshany/OGRe)
[![Twitter @BarakShoshany](https://img.shields.io/twitter/follow/BarakShoshany?style=social)](https://twitter.com/BarakShoshany)
[![Open in Visual Studio Code](https://open.vscode.dev/badges/open-in-vscode.svg)](https://open.vscode.dev/bshoshany/OGRe)

# OGRe: An Object-Oriented General Relativity Package for Mathematica

**Barak Shoshany**\
Department of Physics, Brock University,\
1812 Sir Isaac Brock Way, St. Catharines, Ontario, L2S 3A1, Canada\
[bshoshany@brocku.ca](mailto:bshoshany@brocku.ca) | [https://baraksh.com/](https://baraksh.com/)\
DOI: [10.21105/joss.03416](https://doi.org/10.21105/joss.03416)

* [Version history](#version-history)
    * [v1.7.0 (2021-09-17)](#v170-2021-09-17)
    * [v1.6.1 (2021-09-01)](#v161-2021-09-01)
    * [v1.6 (2021-08-07)](#v16-2021-08-07)
    * [v1.5 (2021-06-07)](#v15-2021-06-07)
    * [v1.4 (2021-05-09)](#v14-2021-05-09)
    * [v1.3 (2021-05-06)](#v13-2021-05-06)
    * [v1.2 (2021-04-28)](#v12-2021-04-28)
    * [v1.1 (2021-04-15)](#v11-2021-04-15)
    * [v1.0 (2021-02-10)](#v10-2021-02-10)

## Version history

### v1.7.0 (2021-09-17)

* New modules:
    * `TCalcGeodesicWithTimeParameter[]` calculates the geodesic equations with respect to the time parameter (which is assumed to be the first coordinate) instead of an affine curve parameter.
    * `TCalcNormSquared[]` calculates the norm-squared of a tensor with respect to its metric, that is, the tensor contracted with itself in all indices. For example, for a vector v<sup>&mu;</sup> the norm-squared will be v<sup>&mu;</sup>v<sub>&mu;</sub> and for a rank-2 tensor T<sup>&mu;&nu;</sup> the result will be T<sup>&mu;&nu;</sup> T<sub>&mu;&nu;</sub>.
    * `TCite[]` displays information on how to cite this package in published research. Thank you for citing my work! :)
* Changes to existing modules:
    * `TCalc*` modules no longer display a progress bar for the calculation. The progress bar wasn't too useful, since almost all tensor calculations are very fast even on an average laptop. Simplifications are the only operations which take any considerable amount of time, and thus progress bars are now displayed only for simplifications.
    * `TCovariantD` and `TPartialD`:
        * If `TCovariantD` is used on a scalar, it is now replaced automatically with `TPartialD` to improve performance, since the covariant derivative of a scalar is just a partial derivative anyway.
        * Fixed bug where acting on scalars incorrectly returned a list of lists instead of a vector.
    * `TGetComponents` now applies a function given in the last argument to the components, and then simplifies them (in parallel, if parallelization is enabled), before returning the components. Note that `TShow` and `TList` have already had this option for some time, but now `TGetComponents` has it too.
    * `TLineElement` now simplifies the line element before returning it.
    * `TSimplify` now allows the user to simplify any expression, not just tensor objects. `TSimplify[expression]` simplifies `expression` based on the user-defined simplification assumptions set using `TSetAssumptions`. If `expression` is a `List`, the components will be simplified in parallel. The user can thus make use of OGRe's optimized simplification process to simplify any Mathematica expression.
* Other changes:
    * Fixed a bug where simplification assumptions were not applied correctly if parallelization was enabled.
    * Added arXiv badge to `README.md` and `CHANGELOG.md`.
    * Added a `CITATION.cff` file (in YAML format) to the GitHub repository. This should add [an option to get a citation in different formats](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/creating-a-repository-on-github/about-citation-files) directly from GitHub repository by clicking on "cite this repository" on the sidebar to the right.
    * Added a `CITATION.bib` file (in BibTeX format) to the GitHub repository. You can use it to easily cite this package in your papers.

### v1.6.1 (2021-09-01)

* This package is now [published in the Journal of Open Source Software](https://joss.theoj.org/papers/10.21105/joss.03416)! I added links to the paper and updated the citation information. The old Zenodo citation should not be used anymore.
* If you use this package in published research, please cite the JOSS paper as detailed in the documentation. A BibTeX entry is also provided.
* The source code itself has not been changed.

### v1.6 (2021-08-07)

* New modules:
    * `TCalcGeodesicFromChristoffel`:
        * Creates a new rank-1 tensor object containing the geodesic equations obtained for each of the coordinates using the Christoffel symbols of the given metric: <i>&#7821;</i><sup>&sigma;</sup> + &Gamma;<sup>&sigma;</sup><sub>&mu;&nu;</sub><i>&#7819;</i><sup>&mu;</sup><i>&#7819;</i><sup>&nu;</sup> = 0.
        * The Christoffel symbols will be calculated automatically using `TCalcChristoffel` if they have not already been calculated.
    * `TCalcGeodesicFromLagrangian`:
        * Creates a new rank-1 tensor object containing the geodesic equations obtained for each of the coordinates by applying the Euler-Lagrange equations to the curve Lagrangian.
        * The Lagrangian will be calculated automatically using `TCalcLagrangian` (see below) if it has not already been calculated.
        * This module leaves the derivatives with respect to the curve parameter in the Euler-Lagrange equation unevaluated (using `Inactive`), which can sometimes help solve the geodesic equations by inspection. Use `Activate` to evaluate the derivatives. (Recall that `TList` and `TShow` can apply a function to the tensor's components before displaying them, so you can write e.g. `TList["ID", Activate]`.)
        * Often the equations obtained in this way will be different from the ones obtained using `TCalcGeodesicFromChristoffel`, but they will always have the same solutions. Usually, one of `TCalcGeodesicFromChristoffel` or `TCalcGeodesicFromLagrangian` will generate simpler equations for a given metric and/or coordinate system.
    * `TCalcLagrangian`: Calculates the curve Lagrangian of a metric, defined as the norm-squared of the tangent to the curve: <i>g</i><sub>&mu;&nu;</sub><i>&#7819;</i><sup>&mu;</sup><i>&#7819;</i><sup>&nu;</sup>. Taking the square root of (the absolute value of) the Lagrangian yields the integrand of the curve length functional. Varying the Lagrangian using the Euler-Lagrange equations yields the geodesic equations (see `TCalcGeodesicFromLagrangian` above).
    * `TMessage`:
        * Not really a module, just a placeholder symbol to which messages not associated with any specific OGRe module are attached.
        * In particular, when a private module (called only internally within the package) invokes `Message`, the message will now be displayed as `TMessage::<message_name>` instead of the awkward ``OGRe`Private`<module_name>::<message_name>``.
        * Not all modules use `TMessage` yet; the transition will be performed gradually in the upcoming releases.
    * `TSetAllowOverwrite`:
        * Allows or disallows overwriting tensors. The default value is `False`, which means you cannot create a new tensor with the same ID as an existing tensor. Calling `TSetAllowOverwrite[True]` will allow overwriting tensors, which is more convenient, but can result in loss of data.
        * You will be warned whenever a tensor is being overwritten, but this warning can be turned off (like any other `Message`) using `Off[TMessage::WarningOverwrite]`.
        * This setting is persistent between sessions.
    * `TSetCurveParameter`:
        * Sets the curve parameter used by `TCalcGeodesicFromChristoffel`, `TCalcGeodesicFromLagrangian`, and `TCalcLagrangian`. These modules will produce results in terms of the coordinate symbols as functions of the curve parameter and their derivatives with respect to this parameter. The default value is &lambda;.
        * If the Lagrangian or geodesic equation vector is displayed using `TList` or `TShow`, the arguments of the coordinate functions are omitted (e.g. <i>x</i> instead of <i>x</i>[&lambda;]) and derivatives with respect to the curve parameter are displayed in Newton (dot) notation (e.g. <i>&#7819;</i> instead of <i>x</i>'[&lambda;]) for improved readability. However, extracting the components using `TGetComponents` will produce the full expressions (e.g. to be used with `DSolve`).
        * When the curve parameter is changed, the parameter of the coordinate functions in all of the tensors calculated so far will be changed accordingly.
    * `TSetReservedSymbols`:
        * Works similar to `TInitializeSymbols`, which has now been removed. However, `TSetReservedSymbols` also saves the reserved symbols so they can be exported and then imported in a later session.
        * If the reserved symbol is a function of the coordinates, `TList` and `TShow` will not show the arguments of the function when displaying the components of a tensor, for improved readability.
    * `TVolumeElementSquared`: Calculates the determinant of a given metric. The square root of the determinant (or its negative, for a pseudo-Riemannian metric) is the volume element.
* Changes to existing modules:
    * All `TCalc*` modules now check if the metric exists first.
    * `TGetComponents`: This module now gets the components of the tensor in the default index configuration and/or coordinate system if either or both are not specified. However, if the default value is used, a message will let you know which representation the components are given in, to avoid confusion.
    * `TInitializeSymbols` has been removed and replaced with `TSetReservedSymbols` (see above).
    * `TList` and `TShow`:
        * Partial derivatives are now displayed in compact notation for improved readability.
        * `TList` will no longer list the same element twice if it is non-zero but equal to minus itself (e.g. `ComplexInfinity`).
        * See `TSetCurveParameter` and `TSetReservedSymbols` above for other changes.
    * `TNewMetric`: If the new metric overrides a previous metric with the same ID, all of the curvature tensors calculated from the metric being overwritten will be automatically deleted, for consistency.
    * `TSetParallelization`:
        * Now uses `$MaxLicenseSubprocesses` instead of the deprecated (as of Mathematica 12.3) `$ConfiguredKernels` to determine how many kernels to launch when enabling parallelization.
        * Disabling parallelization now also closes the kernels.
        * Tensor simplifications will no longer invoke parallelization if the tensor only has one component, to avoid unnecessary overhead.
* Other changes:
    * A button to open the GitHub repository directly in Visual Studio Code has been added to the badges in `README.md`.
* This release is dedicated to my grandfather Yona Shoshany, who taught me BASIC, my first programming language, in my early childhood. He passed away a day before this release was published.

### v1.5 (2021-06-07)

* New modules:
    * `TLineElement`: Displays the line element of a given metric in a coordinate system of your choice.
    * `TSetAutoUpdates`: Enables or disables automatic checks for updates at startup. Note that this setting is persistent between sessions.
* Changes to existing modules:
    * `TInfo`:
        * Whenever this module prints out the name of a related tensor (e.g. the metric or the default coordinates), clicking on the tensor's name will execute `TInfo` for that tensor.
        * This module now provides clickable links for printing out the components of the tensor using `TList` or `TShow`.
        * Calling `TInfo[]` with no parameters lists all the tensors created so far in this session: coordinate systems, metrics, and the tensors associated with each metric.
    * `TSetParallelization`: Calling `TSetParallelization[]` with no parameters now returns the current setting, `True` or `False`.
* Other changes:
    * The startup message is now more compact. Instead of saying "To do X, type Y or click here", it now just says "To do X, type Y", where Y is a clickable link which executes the relevant command.
    * The previous release of the package introduced an incompatibility with Mathematica 12.0 due to the use of the `Splice` function, which only works in v12.1 and later. This has now been fixed, and compatibility with v12.0 has been restored. Although I always develop the package using the latest version of Mathematica, I have now installed Mathematica 12.0 on my computer alongside the latest version, and will use it to ensure compatibility with v12.0 in all future releases of the package.
    * Dynamic content was removed from `OGRe_Documentation.nb` to prevent the message *"This file contains potentially unsafe dynamic content"* from appearing when opening it.
    * The version history has become too long to be included in `README.md`, so I moved it to a separate file, `CHANGELOG.md`.

### v1.4 (2021-05-09)

* Changes to existing modules:
    * `TAddCoordTransformation` now has the clearer syntax `TAddCoordTransformation[sourceID -> targetID, rules]`. The old syntax (with `,` instead of `->`) can still be used.
    * `TChangeID` now has the clearer syntax `TChangeID[oldID -> newID]`. The old syntax (with `,` instead of `->`) can still be used.
    * `TChristoffel`, `TEinsteinTensor`, `TRicciScalar`, `TRicciTensor`, and `TRiemannTensor` have been renamed to `TCalcChristoffel`, `TCalcEinsteinTensor`, `TCalcRicciScalar`, `TCalcRicciTensor`, and `TCalcRiemannTensor` respectively, to group them all together and clarify that they all calculate specific tensors using `TCalc`.
    * `TIndexLetters`, `TParallelize`, and `TSimplifyAssumptions` have been renamed to `TSetIndexLetters`, `TSetParallelization`, and `TSetAssumptions` respectively, to group them all together and clarify that they all change various settings.
    * `TCalc`, `TCalcChristoffel`, `TCalcEinsteinTensor`, `TCalcRicciScalar`, `TCalcRicciTensor`, `TCalcRiemannTensor`, `TChangeID`, `TImport`, `TNewCoordinates`, `TNewMetric`, and `TNewTensor` no longer overwrite existing tensors, to prevent data loss. The user will be instructed to rename the existing tensor using `TChangeID` or delete it using `TDelete` first.
* Other changes:
    * Usage messages:
        * Many usage messages have been simplified, improved, or clarified.
        * Fixed a bug where loading the package directly from GitHub led to usage messages not being displayed properly due to doubling of line breaks.
        * The way that usage messages are defined internally has also been simplified.
    * Automatic checks for update:
        * The package now checks for new versions automatically at startup. No information (personal or otherwise) is sent or collected; the package simply checks the GitHub repository for updates. This check is done asynchronously, to prevent delaying the loading of the package.
        * An upcoming release will introduce persistent storage of user settings, including a setting which will allow disabling the checks for update at startup. Until then, users interested in turning off these automatic checks can change the line `SessionSubmit[StartupCheckForUpdates[]];` in `OGRe.m` to `UpdateMessage = "Checking for updates is disabled.";`.
    * Fixed several lines in the code where a semicolon was missing at the end of the line, causing the code to not run correctly since Mathematica interprets newlines as multiplication by default.
    * This package now has a DOI for citation purposes. Information on how to cite it in publications has been added to the source code and to `README.md`.
    * Added GitHub badges to `README.md`.

### v1.3 (2021-05-06)

* Changes to existing modules:
    * `TExport`, `TExportAll`, `TImport`, and `TImportAll`:
        * Added a clarification in the usage messages for `TExportAll` and `TImportAll` that the directory where the file will be saved or read from can be changed using `SetDirectory`.
        * Importing tensors from a different OGRe version (whether older or newer) will now issue a warning that compatibility issues may occur.
        * `Options` will now be populated correctly with the specific keys relevant to this version, even if you import from other versions which may use different `Options` keys.
* Other changes:
    * The previous update of the package made use of the new Mathematica syntax for anonymous functions, `x |-> body`. However, since this is a new feature introduced in Mathematica 12.2, this accidentally made the package incompatible with earlier Mathematica versions. I changed it to the old syntax `Function[x, body]` in order to maintain compatibility.

### v1.2 (2021-04-28)

* New modules:
    * `TDocs`: Opens the Mathematica notebook `OGRe_Documentation.nb` directly from the GitHub repository. This allows users to instantly access the latest documentation from any Mathematica session, whether the package itself was loaded locally or from GitHub.
    * `TParallelize`:
        * Enables or disable the parallelization of tensor simplifications (see below). It is disabled by default.
        * As a rule of thumb, if simplifications are taking less than a few seconds, then you should leave parallelization off, as it has a small overhead and may actually impede performance in that case. However, if simplifications are taking more than a few seconds, then it is highly recommended to enable parallelization for a significant performance boost.
        * This setting will be exported when using `TExportAll`, so it will be imported automatically with `TImportAll` in future sessions.
        * Aspects of the package other than simplification are not parallelized, since they are typically short calculations with few elements, and would not benefit from parallelization.
* Changes to existing modules:
    * `TCalc` now displays a progress bar (see below).
    * `TCheckForUpdates`:
        * Fixed a bug where this module falsely reported that a new version is available.
        * Fixed a bug where this module issued an error when giving the option to download `OGRe.m` to the current notebook directory if the notebook is an Untitled notebook, meaning it is not an actual file in the file system and thus does not have a notebook directory.
        * This module now gives the user the option to reload the new version directly from GitHub without saving it locally.
    * `TList` and `TShow`:
        * The output of these modules is now uneditable. I made this change after I noticed that some users were editing the output, thinking that they are editing the tensor itself. (An option to actually edit the elements of a tensor in `MatrixForm` may be added in a future version.)
        * These modules now take a fourth optional argument: a function to map to each of the tensor's elements before they are displayed. Typically this would be `ReplaceAll[rules]` to apply the `rules` to the elements, but any function can be used. Note that applying `ReplaceAll` to the output of `TShow` or `TList` directly won't work, since the output is formatted.
    * `TNewMetric` now verifies that the metric components are an invertible, symmetric, and square matrix.
* Other changes:
    * The performance of simplifications has been improved considerably!
        * Simplifications are now performed automatically less often, to prevent simplifying the same expression multiple times.
        * In addition, simplifications can now be parallelized. This is enabled using the new module `TParallelize` (see above). The parallelization is implemented using a thread pool model (with `ParallelSubmit` and `WaitAll`), which in my benchmarks was found to provide better performance than the simpler `Parallelize` or `ParallelMap`.
        * Parallelization can potentially provide a speedup proportional to the number of parallel kernels. Note that this number is determined by your Mathematica license, and it may be less than the number of cores in your CPU.
    * Calculations with `TCalc` and simplifications now both display progress bars, so you can keep track of long calculations. This was a popular feature request, especially for simplifications, since they are usually the most time-consuming operations.
    * Improved the appearance of the welcome message.
    * Messages from OGRe are now displayed with a special cell label.
    * Information on how to load the package is now displayed in the `README.md` file.

### v1.1 (2021-04-15)

* New modules:
    * `TCheckForUpdates`: Automatically checks the GitHub repository for updates. If a new version of the package is available, the module will offer an option to download or install the new version and reload the package.
    * `TGetComponents`: Extracts the components of a tensor object in a specific representation as a `List`.
    * `TInfo`: Displays information about a tensor object, and any other objects associated with it, in human-readable form.
* Changes to existing modules:
    * `TAddCoordTransformation` now also calculates the Jacobian of the transformation when it is executed, and stores the result for future use within the tensor object of the source coordinates. This improves performance by using the pre-calculated Jacobian whenever a coordinate transformation is performed, instead of calculating it from scratch every time, as was the case in v1.0.
    * `TChristoffel` now explicitly marks the resulting tensor object as having special transformation properties. The Levi-Civita connection, whose components are the Christoffel symbols, does not transform like a tensor under a coordinate transformation, and OGRe now automatically knows to use the correct transformation rule. **Please note that tensors created with `TChristoffel` in v1.0 will not transform correctly, so they should be recalculated after updating.**
    * `TExport` now adds a key named `"OGReVersion"` which records the version of the package used to create the tensor being exported. This will be used in future versions to ensure backwards compatibility.
    * `TExportAll` now exports, in addition to all the tensors defined so far, a special key, `Options`, containing information about the current session. Currently, this key stores the version of the package, the index letters to use, and the simplification assumptions set by the user. When importing the data in another session using `TImportAll`, the version number will be used to ensure backwards compatibility, and the other options will be used to restore any user-defined index letters and simplification assumptions made during the session.
    * `TNewTensor` now allows defining the components of the new tensor in any coordinate system. In v1.0, the components had to be defined in the default coordinate system of the associated metric. This is still the default behavior if a coordinate system is not specified, for compatibility with v1.0, but it is recommended to always specify the coordinate system explicitly, to avoid accidentally defining the tensor with the wrong components.
    * `TSimplifyAssumptions` now appends new simplification assumptions to the list of previously added  assumptions, instead of replacing it. Also, OGRe now automatically assumes that all variables are real, which helps simplify certain expressions. If you are using more exotic variables, use `TSimplifyAssumptions[!Reals]` to disable this assumption.
* Other changes:
    * The `"Role"` key of each tensor object now indicates how that object was created. In v1.0, tensors created with `TNewCoordinates` had the role `"Coordinates"` and tensors created with `TNewMetric` had the role `"Metric"`, but all other tensors had the role `"General"`. Now tensors will have the roles `"Tensor"`, `"Calculated"`, `"Christoffel"`, `"Riemann"`, `"Ricci Tensor"`, `"Ricci Scalar"`, or `"Einstein"` if they were created using `TNewTensor`, `TCalc`, `TChristoffel`, `TRiemannTensor`, `TRicciTensor`, `TRicciScalar`, or `TEinstein` respectively. This is currently just for bookkeeping, but may have other uses in future versions. Note that tensors imported from v1.0 will still have the role `"General"`.
    * Improved the formatting of the usage messages for all OGRe modules. They no longer break in the middle of words.
    * Debug mode has been removed. If the package is loaded more than once in a single session, it will redefine all symbols but keep any previously defined tensors intact. This is useful both for debugging and for reloading the package after an update.
* Bug fixes:
    * Fixed a bug where the partial and covariant derivatives of a scalar were not calculated correctly.

### v1.0 (2021-02-10)

* Initial release.
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03416/status.svg)](https://doi.org/10.21105/joss.03416)
[![arXiv:2109.04193](https://img.shields.io/badge/arXiv-2109.04193-b31b1b.svg)](https://arxiv.org/abs/2109.04193)
[![License: MIT](https://img.shields.io/github/license/bshoshany/thread-pool)](https://github.com/bshoshany/OGRe/blob/master/LICENSE.txt)
![Language: Mathematica 12](https://img.shields.io/badge/Language-Mathematica%2012-yellow)
![File size in bytes](https://img.shields.io/github/size/bshoshany/OGRe/OGRe.m)
![GitHub last commit](https://img.shields.io/github/last-commit/bshoshany/OGRe)
[![GitHub repo stars](https://img.shields.io/github/stars/bshoshany/OGRe?style=social)](https://github.com/bshoshany/OGRe)
[![Twitter @BarakShoshany](https://img.shields.io/twitter/follow/BarakShoshany?style=social)](https://twitter.com/BarakShoshany)
[![Open in Visual Studio Code](https://open.vscode.dev/badges/open-in-vscode.svg)](https://open.vscode.dev/bshoshany/OGRe)

# OGRe: An Object-Oriented General Relativity Package for Mathematica

**Barak Shoshany**\
Department of Physics, Brock University,\
1812 Sir Isaac Brock Way, St. Catharines, Ontario, L2S 3A1, Canada\
[bshoshany@brocku.ca](mailto:bshoshany@brocku.ca) | [https://baraksh.com/](https://baraksh.com/)\
DOI: [10.21105/joss.03416](https://doi.org/10.21105/joss.03416)

* [Summary](#summary)
* [Features](#features)
* [Installing and loading the package](#installing-and-loading-the-package)
* [Documentation](#documentation)
* [Issue and pull request policy](#issue-and-pull-request-policy)
* [Copyright and citing](#copyright-and-citing)

## Summary

OGRe is a modern Mathematica package for differential geometry and tensor calculus, designed to be both powerful and user-friendly. It can be used in a variety of contexts where tensor calculations are needed, in both mathematics and physics, but it is especially suitable for general relativity.

Tensors are abstract objects, which can be represented as multi-dimensional arrays once a choice of index configuration and coordinate system is made. OGRe stays true to this definition, but takes away the complexities that come with combining tensors in different representations. This is done using an object-oriented programming approach, taking advantage of principles such as encapsulation and class invariants to eliminate the possibility of user error.

The user initially defines each tensor in OGRe using its explicit components in any single representation. Operations on this tensor are then done abstractly, without needing to specify which representation to use. Possible operations include addition of tensors, multiplication of tensor by scalar, trace, contraction, and partial and covariant derivatives.

OGRe will automatically choose which representation to use for each tensor based on how the tensors are combined. For example, if two tensors are added, then OGRe will automatically use the same index configuration for both. Similarly, if two tensors are contracted, then OGRe will automatically ensure that the contracted indices are one upper (contravariant) and one lower (covariant). OGRe will also automatically transform all tensors being operated on to the same coordinate system.

Transformations between representations are done behind the scenes; all the user has to do is specify which metric to use for raising and lowering indices, and how to transform between the coordinate systems being used. This also means that there is no room for user error. The user cannot mistakenly perform "illegal" operations such as 2A<sup>&mu;&nu;</sup>+B<sub>&mu;&lambda;</sub>C<sub>&lambda;&nu;</sub>. Instead, the user simply inputs the names of the tensors, the order (but **not** the configuration) of indices for each, and the operations to perform - and the correct combination 2A<sup>&mu;&nu;</sup>+B<sup>&mu;</sup><sub>&lambda;</sub>C<sup>&lambda;&nu;</sup> will be automatically deduced.

I initially created OGRe for use in my own research, so I made it as flexible and powerful as possible. I also wanted my students to be able to use it easily and efficiently, even if they only have minimal experience with Mathematica and/or general relativity, so I made it simple to learn and easy to use. As a result, this package is equally suitable for both experienced and novice researchers.

## Features

* Define coordinate systems and the transformation rules between them. Tensor components are then transformed automatically between coordinates behind the scenes as needed.
* Each tensor is associated with a specific metric. Tensor components are then transformed automatically between different index configurations, raising and lowering indices behind the scenes as needed.
* Display any tensor in any index configuration and coordinate system, either in vector/matrix form or as a list of all unique non-zero elements. Metrics can also be displayed as a line element.
* Automatically simplify tensor components, optionally with user-defined simplification assumptions. Simplifications can be parallelized for a significant performance boost. The user may utilize the package's simplification algorithm for any Mathematica expression, not just tensors.
* Export tensors to a Mathematica notebook or to a file, so they can later be imported into another Mathematica session without having to redefine them from scratch.
* Highly customizable. User settings are exported and imported along with the tensors. Some settings are persistent between sessions.
* Easily calculate arbitrary tensor formulas using any combination of addition, multiplication by scalar, trace, contraction, partial derivative, and covariant derivative.
* Built-in modules for calculating the Christoffel symbols (Levi-Civita connection), Riemann tensor, Ricci tensor and scalar, Einstein tensor, curve Lagrangian, and volume element from a metric.
* Calculate the norm-squared of tensors of any rank.
* Calculate the geodesic equations in terms of a user-defined affine curve parameter, in two different ways: from the Christoffel symbols or from the curve Lagrangian. For spacetime metrics, the geodesic equations can be calculated in terms of the time coordinate.
* Designed with speed and performance in mind, using optimized algorithms designed specifically for this package.
* Fully portable. Can be imported directly from the web into any Mathematica notebook, without downloading or installing anything. Integrates seamlessly with the Wolfram Cloud.
* Clear and detailed documentation, with many examples, in both [Mathematica notebook](OGRe_Documentation.nb) and [PDF](OGRe_Documentation.pdf) format. Detailed usage messages are also provided.
* Open source. The code is extensively documented; please feel free to fork and modify it as you see fit.
* Under continuous and active development. Bug reports and feature requests are welcome, and should be made via [GitHub issues](https://github.com/bshoshany/OGRe/issues).

## Installing and loading the package

This package is compatible with Mathematica 12.0 or newer. It consists of only one file, `OGRe.m`. There are several different ways to load the package:

* **Run from local file with installation:** This is the recommended option, as it allows you to permanently use the package offline from any Mathematica notebook. Download the file `OGRe.m` from [the GitHub repository](https://github.com/bshoshany/OGRe) and copy it to the directory given by `FileNameJoin[{$UserBaseDirectory, "Applications"}]`. The package may now be loaded from any notebook by writing ``Needs["OGRe`"]`` (note the backtick <code>&#96;</code> following the word OGRe).

* **Run from local file without installation:** This option allows you to use the package in a portable fashion, without installing it in the `Applications` directory. Download the file `OGRe.m` from [the GitHub repository](https://github.com/bshoshany/OGRe), place it in the same directory as the notebook you would like to use, and use the command `Get["OGRe.m", Path->NotebookDirectory[]]` to load the package.

* **Run from web with installation:** This option allows you to automatically download and install the package on any computer. Simply run the command `URLDownload["https://raw.githubusercontent.com/bshoshany/OGRe/master/OGRe.m", FileNameJoin[{$UserBaseDirectory, "Applications", "OGRe.m"}]]` from any Mathematica notebook to permanently install the package. Then use ``Needs["OGRe`"]`` from any notebook to load it.

* **Run from web without installation:** This option allows you to use the package from any Mathematica notebook on any computer, without manually downloading or installing it, as long as you have a working Internet connection. It also ensures that you always use the latest version of the package, but be aware that updates may sometimes not be fully backwards compatible. Simply write `Get["https://raw.githubusercontent.com/bshoshany/OGRe/master/OGRe.m"]` in any Mathematica notebook to load the package.

To uninstall the package, just delete the file from the `Applications` directory, which can be done from within Mathematica using the command `DeleteFile[FileNameJoin[{$UserBaseDirectory, "Applications", "OGRe.m"}]]`.

## Documentation

The full and detailed documentation for this package may be found in the following repository files:

* `OGRe_Documentation.nb`: An interactive Mathematica notebook. Requires Mathematica to open.
* `OGRe_Documentation.pdf`: A PDF version of the notebook. Can be viewed with any PDF reader.

Once the package is loaded, the documentation can be easily accessed by executing the command `TDocs[]`, which automatically downloads the file `OGRe_Documentation.nb` from GitHub and opens it in Mathematica.

## Issue and pull request policy

This package is under continuous and active development. If you encounter any bugs, or if you would like to request any additional features, please feel free to [open a new issue on GitHub](https://github.com/bshoshany/OGRe/issues) and I will look into it as soon as I can.

Contributions are always welcome. However, I release my projects in cumulative updates after editing and testing them locally on my system, so my policy is not to accept any pull requests. If you open a pull request, and I decide to incorporate your suggestion into the project, I will first modify your code to comply with the project's coding conventions (formatting, syntax, naming, comments, programming practices, etc.), and perform some tests to ensure that the change doesn't break anything. I will then merge it into the next release of the project, possibly together with some other changes. The new release will also include a note in `CHANGELOG.md` with a link to your pull request, and modifications to the documentation in `README.md` as needed.

## Copyright and citing

Copyright (c) 2021 [Barak Shoshany](http://baraksh.com). Licensed under the [MIT license](LICENSE.txt).

If you use this package in published research, please cite it as follows:

* Shoshany, B., (2021). OGRe: An Object-Oriented General Relativity Package for Mathematica. Journal of Open Source Software, 6(65), 3416, <https://doi.org/10.21105/joss.03416>

You can also use the following BibTeX entry:

```none
@article{Shoshany2021_OGRe,
    author    = {Barak Shoshany},
    doi       = {10.21105/joss.03416},
    journal   = {Journal of Open Source Software},
    number    = {65},
    pages     = {3416},
    publisher = {The Open Journal},
    title     = {OGRe: An Object-Oriented General Relativity Package for Mathematica},
    url       = {https://doi.org/10.21105/joss.03416},
    volume    = {6},
    year      = {2021}
}
```
