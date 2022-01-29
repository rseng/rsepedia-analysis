# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Changes since the previous release can be found in the [news](./news) directory.
A raw git diff can be seen [here][unreleased].

<!--next-->

## [4.0.0] 2021-08-09

## Feature

-   [#849](https://github.com/SMI/SmiServices/pull/849) by jas88. CTPAnonymiser
    refactoring
    -   Reduce memory footprint (issue #837)
    -   Simplify RabbitMQ message handling
    -   Stop creating temporary copy of input file - no longer needed in EPCC
        environment without Lustre FS (per issue #836)
    -   Add checks input file is readable not just extant, hopefully fixing
        issue #533
-   [#861](https://github.com/SMI/SmiServices/pull/861) by rkm. Add Equ to
    automatically implement equality members for classes.
-   [#878](https://github.com/SMI/SmiServices/pull/878) by rkm. Update RDMP
    packages with replacement of System.Data.SqlClient with
    Microsoft.Data.SqlClient. Replace usages of same in codebase

## Bugfix

-   [#764](https://github.com/SMI/SmiServices/pull/764) by howff. Clean up the
    Python code lint after running pylint3
-   [#841](https://github.com/SMI/SmiServices/pull/841) by tznind. Fixed bug
    when disposing `CsvDestination` instances that have not begun writing any
    output
-   [#880](https://github.com/SMI/SmiServices/pull/880) by tznind. Fixed edge
    case in IdentifierMapper when a dicom tag has illegal multiplicity in
    PatientID field

## Meta

-   [#843](https://github.com/SMI/SmiServices/pull/843) by rkm. Add pre-commit
    and codespell. Fix all current spelling mistakes
-   [#844](https://github.com/SMI/SmiServices/pull/844) by rkm. Fixup regex in
    codespell config
-   [#855](https://github.com/SMI/SmiServices/pull/855) by rkm. Specify dotnet
    SDK version in global.json
-   [#859](https://github.com/SMI/SmiServices/pull/859) by rkm. Bump LangVersion
    to 9.0
-   [#876](https://github.com/SMI/SmiServices/pull/876) by rkm. Add check and
    error message for missing coveralls token
-   [#877](https://github.com/SMI/SmiServices/pull/877) by rkm. Fix setting
    replication for MongoDB in Windows CI pipelines

## Removal

-   [#848](https://github.com/SMI/SmiServices/pull/848) by rkm. Removed
    NationalPACSAccessionNumber from all metadata

## [3.2.1] 2021-07-07

## Bugfix

-   [#834](https://github.com/SMI/SmiServices/pull/834) by tznind. Improved
    logging and fixed yaml options not being respected in IsIdentifiableReviewer

## [3.2.0] 2021-07-05

## Feature

-   [#795](https://github.com/SMI/SmiServices/pull/795) by tznind. Added support
    for specifying IsIdentifiable CLI options in the yaml config files instead
    of command line (command line will always take precedence if both are
    specified)
-   [#797](https://github.com/SMI/SmiServices/pull/797) by tznind. Added custom
    themes to IsIdentifiableReviewer. Use flag `--theme mytheme.yaml`
-   [#801](https://github.com/SMI/SmiServices/pull/801) by tznind. Added help
    and cancel buttons to custom pattern dialog in reviewer
-   [#802](https://github.com/SMI/SmiServices/pull/802) by tznind. Added
    IsIdentifiableReviewer settings into main yaml config
-   [#818](https://github.com/SMI/SmiServices/pull/818) by tznind. Added Rules
    Manager to IsIdentifiable Reviewer

## Bugfix

-   [#787](https://github.com/SMI/SmiServices/pull/787) by rkm. Fix the call to
    the release changelog script
-   [#788](https://github.com/SMI/SmiServices/pull/788) by rkm. Require all CI
    tests to pass before packaging runs
-   [#789](https://github.com/SMI/SmiServices/pull/789) by rkm. Don't upload
    coverage for tagged builds
-   [#796](https://github.com/SMI/SmiServices/pull/796) by tznind. Fixed bug
    opening corrupted reports in IsIdentifiableReviewer crashing the application
-   [#806](https://github.com/SMI/SmiServices/pull/806) by tznind. Replace use
    of GlobalColorScheme with the proper static members in Terminal.Gui that
    will propagate correctly for all windows without having to set them
    manually.
-   [#811](https://github.com/SMI/SmiServices/pull/811) by rkm. Fix coverage
    task always running even if a previous task failed

## Doc

-   [#823](https://github.com/SMI/SmiServices/pull/823) by tznind. Refresh
    documentation for IsIdentifiableReviewer

## [3.1.0] 2021-06-11

## Feature

-   [#759](https://github.com/SMI/SmiServices/pull/759) by tznind. Added
    parallelisation to load process in IsIdentifiableReviewer rules view

## Bugfix

-   [#756](https://github.com/SMI/SmiServices/pull/756) by howff. Open CSV file
    read-only
-   [#771](https://github.com/SMI/SmiServices/pull/771) by tznind.
    IsIdentifiableReviewer:
    -   Added --usc (UseSystemConsole) for alternative display driver based on
        System.Console
    -   Removed modal dialog that could cause errors opening a previously
        completed report
    -   Added label with currently opened file and fixed ignore/update labels
    -   Added spinner indicator for when loading the Next report in sequential
        mode takes a while
    -   Fixed bug where Ctrl+Q in Ignore/Update with custom patterns in
        RulesView results in hard crash
-   [#785](https://github.com/SMI/SmiServices/pull/785) by tznind. Fixed bug
    with multiple enumeration during loading very large failure reports in
    IsIdentifiableReviewer

## Meta

-   [#773](https://github.com/SMI/SmiServices/pull/773) by rkm. Add code
    coverage
-   [#775](https://github.com/SMI/SmiServices/pull/775) by rkm. Move useful
    scripts from .azure-pipelines/scripts to utils. Update utils/README.md.
-   [#781](https://github.com/SMI/SmiServices/pull/781) by rkm. Fixup coverage
    variables between pushes/PRs

## [3.0.2] 2021-05-14

## Bugfix

-   [#745](https://github.com/SMI/SmiServices/pull/745) by tznind. Fixed
    reviewer tree view refresh/async code

## [3.0.1] 2021-05-06

## Feature

-   [#738](https://github.com/SMI/SmiServices/pull/738) by rkm. Improvements to
    CLI user experience
    -   Immediately verify that the config file (GlobalOptions) we've loaded is
        somewhat valid
    -   Improve visibility of exception messages on CLI exit

## Bugfix

-   [#737](https://github.com/SMI/SmiServices/pull/737) by rkm. Switch ordering
    of annotations in ExtractImagesCliOption to fix a runtime exception.

## Meta

-   [#736](https://github.com/SMI/SmiServices/pull/736) by rkm. Remove .NET Core
    2.2 runtime from the Azure Pipelines builds

## [3.0.0] 2021-05-06

## Feature

-   [#702](https://github.com/SMI/SmiServices/pull/702) by rkm. Replace Java
    ExtractorCLI with C# ExtractImages service
    -   _Breaking_ Existing scripts and documentation
    -   _Breaking_ Change ExtractorClOptions -> ExtractImagesOptions in YAML
        configs
-   [#713](https://github.com/SMI/SmiServices/pull/713) by rkm. Upgrade all
    projects and related CI scripts to `net5`.
-   [#734](https://github.com/SMI/SmiServices/pull/734) by tznind. Updated to
    RDMP 5.0.0 (and Dicom Plugin 3.0.0)

## Bugfix

-   [#708](https://github.com/SMI/SmiServices/pull/708) by rkm. Add nuget.config
    file to fix flaky issues with Azure CI Runners. Ref:
    -   https://github.com/NuGet/Home/issues/10586
    -   https://github.com/actions/virtual-environments/issues/3038
-   [#714](https://github.com/SMI/SmiServices/pull/714) by rkm. Fixup current
    LGTM alerts.
-   [#715](https://github.com/SMI/SmiServices/pull/715) by rkm. Revert to the
    hack-y method of fixing the Nuget cache for Azure Windows agents.
-   [#722](https://github.com/SMI/SmiServices/pull/722) by rkm. Remove
    workaround for Windows agents in Azure CI since fixed upstream. Ref:
    -   https://github.com/actions/virtual-environments/issues/3038

## [2.1.1] 2021-04-07

## Bugfix

-   [#697](https://github.com/SMI/SmiServices/pull/697) by rkm. Fixes #695.
    Removes the checks preventing modality being specified with other extraction
    keys
-   [#698](https://github.com/SMI/SmiServices/pull/698) by rkm. CohortExtractor
    database queries that crash during execution are now logged
-   [#701](https://github.com/SMI/SmiServices/pull/701) by howff. Several
    improvements to Python code for handling unusually-formatted SR documents.
-   [#704](https://github.com/SMI/SmiServices/pull/704) by rkm. Fix
    ReportNewLine being incorrectly set to a pre-escaped string. Fixes #703

## Doc

-   [#672](https://github.com/SMI/SmiServices/pull/672) by howff.
    IsIdentifiableReviewer document updated

## [2.1.0] 2021-03-30

## Feature

-   [#676](https://github.com/SMI/SmiServices/pull/676) by tznind. Improvements
    to IsIdentifiable Reviewer
    -   Tab based navigation
    -   Better pattern generation for overlapping failure parts
    -   Fixed sequential mode always showing failures covered by existing UPDATE
        rules
    -   Fixed tree view not adapting as new rules are added e.g. when used to
        interactively process many failures
-   [#690](https://github.com/SMI/SmiServices/pull/690) by tznind. Added
    ModalitySpecificRejectorOptions

## Bugfix

-   [#678](https://github.com/SMI/SmiServices/pull/678) by rkm. Fix log
    directory naming for single-entrypoint app. Fixes #677
-   [#681](https://github.com/SMI/SmiServices/pull/681) by howff. Fix Python
    packaging when run from repo root dir
-   [#684](https://github.com/SMI/SmiServices/pull/684) by howff.
    -   Minor fixes to python for compatibility (eg. semehr using python2)
    -   Use a unique temporary directory for input and output files
    -   Enable extraction from MongoDB by StudyDate (the only field other than
        FilePath which has an index)
-   [#688](https://github.com/SMI/SmiServices/pull/688) by howff. Reword the
    CohortExtractor log/info messages as it's not just patients which can be
    rejected

## Meta

-   [#668](https://github.com/SMI/SmiServices/pull/668) by rkm. Add note to the
    releasing documentation to alert the Mattermost channel

## [2.0.0] 2021-03-13

## Feature

-   [#618](https://github.com/SMI/SmiServices/pull/618) by tznind.
    IsIdentifiableReviewer
    -   Added progress when loading large files (with cancellation support)
    -   Now groups outstanding failures by column
    -   Fixed rules being flagged as 'Identical' when classifying different
        input columns
-   [#634](https://github.com/SMI/SmiServices/pull/634) by rkm. Convert to
    single entry-point app
    -   Breaking: Existing scripts and processes which reference the old
        applications
-   [#636](https://github.com/SMI/SmiServices/pull/636) by howff. Improvements
    to Python scripts, tests, documentation
-   [#647](https://github.com/SMI/SmiServices/pull/647) by howff. Improved
    testing of SR anonymiser with a standalone stub
-   [#662](https://github.com/SMI/SmiServices/pull/662) by tznind. Added ability
    to ignore multiple failures at once in Is Identifiable Reviewer

## Bugfix

-   [#597](https://github.com/SMI/SmiServices/pull/597) by tznind. Fixed
    ConsensusRules not being run
-   [#619](https://github.com/SMI/SmiServices/pull/619) by jas88. Reduce memory
    consumption in nerd
-   [#632](https://github.com/SMI/SmiServices/pull/632) by rkm. Normalise the
    IsIdentifiableReviewer namespace to match the other Applications
-   [#646](https://github.com/SMI/SmiServices/pull/646) by howff. Call the
    garbage collector (and report progress) every 1000 records processed from a
    database source in IsIdentifiable.
-   [#650](https://github.com/SMI/SmiServices/pull/650) by howff. Create python
    package called SmiServices so a wheel can be created.
    -   Rename all the imports for the new package name
    -   Remove references to PYTHONPATH Replace all paths and relative imports
        to be independent of current directory Fixes to python tests
-   [#656](https://github.com/SMI/SmiServices/pull/656) by howff. SRAnonTool
    updated to handle the output from the latest SemEHR anonymiser and ignore
    None-type annotations
-   [#661](https://github.com/SMI/SmiServices/pull/661) by tznind. Fixed layout
    of main window so it no longer obscures classification/type
-   [#665](https://github.com/SMI/SmiServices/pull/665) by tznind. Fixed tree
    view losing selected index when updating/ignoring a failure (in tree view)
-   [#666](https://github.com/SMI/SmiServices/pull/666) by howff. Silence
    deprecation warning from newer Python as noted by the azure pipeline test

## Meta

-   [#588](https://github.com/SMI/SmiServices/pull/588) by rkm. Prevent
    additional language packs being included in published packages. Reduces
    overall package size a bit.
-   [#592](https://github.com/SMI/SmiServices/pull/592) by rkm. Manually update
    csvhelper to 25.0.0, fo-dicom to 4.0.7
-   [#616](https://github.com/SMI/SmiServices/pull/616) by rkm. Check for
    clobbered files during package build
-   [#620](https://github.com/SMI/SmiServices/pull/620) by rkm. Replace the
    legacy SecurityCodeScan with SecurityCodeScan.VS2019
-   [#637](https://github.com/SMI/SmiServices/pull/637) by rkm. Change to
    tracking PR changes in news fragment files. Add a script to auto-update the
    CHANGELOG from these files.
-   [#648](https://github.com/SMI/SmiServices/pull/648) by rkm. Remove the
    temporary reference to BadMedicine added in
    [#592](https://github.com/SMI/SmiServices/pull/592)
-   [#654](https://github.com/SMI/SmiServices/pull/654) by howff. Add Azure
    Pipelines test and packaging for the Python scripts

## [1.15.1] 2021-02-17

### Fixed

-   [#610](https://github.com/SMI/SmiServices/pull/610) by `howff`. Fixed CTP
    logging

## [1.15.0] 2021-02-17

### Changed

-   [#575](https://github.com/SMI/SmiServices/pull/575) by `rkm`. Standardised
    logging setup and Program entries across whole solution
    -   Breaking: YAML config change required
    -   Removes the `SMI_LOGS_ROOT` variable - now in YAML config
    -   Removes the `--trace-logging` CLI option - now in YAML config
    -   All invocations of IsIdentifiable now require a YAML config to ensure
        logging is properly configured
-   [#577](https://github.com/SMI/SmiServices/pull/577) by `rkm`. Simplify
    branch workflow by dropping develop

### Fixed

-   [#581](https://github.com/SMI/SmiServices/pull/581) by `rkm`. Fixed a bug
    where newlines would never be correctly parsed from the config option in
    CohortPackager
-   [#597](https://github.com/SMI/SmiServices/pull/597) by `tznind`. Fixed
    ConsensusRules not being run

### Dependencies

-   Bump CsvHelper from 22.1.1 to 22.1.2
-   Bump HIC.RDMP.Plugin from 4.2.3 to 4.2.4
-   Bump HIC.RDMP.Plugin.Test from 4.2.3 to 4.2.4
-   Bump Magick.NET-Q16-AnyCPU from 7.23.1 to 7.23.2
-   Bump SecurityCodeScan from 3.5.3 to 3.5.4
-   Bump System.Drawing.Common from 5.0.0 to 5.0.1
-   Bump System.IO.Abstractions from 13.2.9 to 13.2.11
-   Bump System.IO.Abstractions.TestingHelpers from 13.2.9 to 13.2.11
-   Bump jansi from 2.2.0 to 2.3.1
-   Bump junit from 4.13.1 to 4.13.2

## [1.14.1] - 2021-02-04

### Fixed

-   [#576](https://github.com/SMI/SmiServices/pull/576) Fixup Windows package
    build

## [1.14.0] - 2021-02-04

### Added

-   Added total job duration to extraction report header
-   IsIdentifiableReviewer rule review screen
-   Added CSV input support for IsIdentifiable, use verb `file` from command
    line
-   Updater microservice now audits performance of queries (cumulative affected
    rows, queries executed etc)
-   Added `-f` option to DicomTagReader to read a single zip/dicom file
-   Added some Python library code

### Changed

-   Clarified the CLI help text for `--format` in CohortPackager
-   CTP calls an external program to anonymise Structured Reports
    -   Requires an addition to `default.yaml`:
        `CTPAnonymiserOptions.SRAnonTool` (adding this does not break existing
        programs).
-   Consolidate System.IO.Abstractions.TestingHelpers package references into
    the Smi.Common.Tests package
-   Tidy common csproj options into `Directory.Build.props` files for all, src,
    and test projects
-   Replace TravisCI and AppVeyor builds with Azure Pipelines

### Fixed

-   CohortPackager: Don't try and create the jobId file when recreating an
    existing report
-   CohortPackager.Tests: Fix a flaky test caused by NUnit setup/teardown code
    when running tests in parallel
-   CohortPackager.Tests: Fix a flaky test caused by using the same MongoDB
    database name when running tests in parallel
-   Fixed the existing CTPAnonymiser tests which had not been updated for the
    SRAnonTool changes
-   Fixed executable name on UpdateValues microservice

### Dependencies

-   Bump CsvHelper from 17.0.1 to 22.1.1
-   Bump HIC.RDMP.Dicom from 2.1.11 to 2.2.2
-   Bump HIC.RDMP.Plugin from 4.1.9 to 4.2.3
-   Bump HIC.RDMP.Plugin.Test from 4.1.9 to 4.2.3
-   Bump Magick.NET-Q16-AnyCPU 7.22.2.2 to 7.23.1
-   Bump Microsoft.NET.Test.Sdk 16.8.0 to 16.8.3
-   Bump Moq from 4.15.2 to 4.16.0
-   Bump NUnit from 3.12.0 to 3.13.1
-   Bump NunitXml.TestLogger from 2.1.80 to 3.0.91
-   Bump System.IO.Abstractions from 13.2.2 to 13.2.9
-   Bump System.IO.Abstractions.TestingHelpers from 13.2.2 to 13.2.9
-   Bump YamlDotNet from 9.1.0 to 9.1.4

## [1.13.0] - 2020-12-03

### Added

-   Added new command line application TriggerUpdates for detecting and issuing
    UpdateValuesMessages (e.g. ECHI mapping changes)
-   Added new service UpdateValues which propagates changes (e.g. ECHI mapping
    changes) throughout the deployed database tables.
-   ConsensusRule for combining 2+ other rules e.g. SocketRules (See
    IsIdentifiable Readme.md for more details)
-   Added runtime and total failures count to IsIdentifiable logs
-   Added NoSuffixProjectPathResolver which generates anonymous image path names
    that do not contain "-an" (which is the default behaviour).
    -   To use, set `CohortExtractorOptions.ProjectPathResolverType` to
        `Microservices.CohortExtractor.Execution.ProjectPathResolvers.NoSuffixProjectPathResolver`
    -   For identifiable extractions, the NoSuffixProjectPathResolver is now
        used
-   Validation reports can now be created as either "Combined" (single report as
    before" or "Split" (a
    [pack](src/microservices/Microservices.CohortPackager/README.md) of reports
    including CSVs suitable for post-processing). This is configurable in the
    YAML config and can also be specified on the CLI when recreating reports for
    an extraction
-   Added JobCompletedAt to the validation reports
-   IsIdentifiable: Add support for ignoring OCR output less than `n` characters
    in length
-   IsIdentifiable: Add a test case for burned-in image text

### Changed

-   Update docs and make more keywords links to the relevant docs (#440)
-   Reduce memory usage on long-running microservices even when .Net assumes RAM
    is plentiful
-   Validation reports are now written to the project reports directory, instead
    of to a central reports directory

### Fixed

-   Fix mismatch in Java/C# messages for ExtractionModality
-   ExtractionFileCopier: Copy files relative to the extraction root not the
    global filesystem root
-   Fix implementation of minimum OCR length (before being reported) #471

### Dependencies

-   Bump CsvHelper from 17.0.0 to 17.0.1
-   Bump System.IO.Abstractions from 13.2.1 to 13.2.2
-   Bump Moq from 4.15.1 to 4.15.2
-   Bump System.IO.Abstractions.TestingHelpers from 13.2.1 to 13.2.2
-   Bump CsvHelper from 16.2.0 to 17.0.0
-   Bump JetBrains.Annotations from 2020.1.0 to 2020.3.0
-   Bump jackson-dataformat-yaml from 2.11.3 to 2.12.0
-   Bump jackson-databind from 2.11.3 to 2.12.0

## [1.12.2] - 2020-09-18

-   Fix missing JSON fields from CTP output

## [1.12.1] - 2020-09-15

-   Remove reference to MongoDB.Driver in Smi.Common.MongoDb.csproj since it
    caused a version conflict in the output packages

## [1.12.0] - 2020-09-14

### Added

-   [breaking] Add identifiable extraction support
    -   New service "FileCopier" which sits in place of CTP for identifiable
        extractions and copies source files to their output dirs
    -   Changes to MongoDB extraction schema, but backwards compatibility has
        been tested
    -   RabbitMQ extraction config has been refactored. Queues and service
        config files need to be updated
-   Add [SecurityCodeScan](https://security-code-scan.github.io/) tool to build
    chain for .NET code
-   Add "no filters" extraction support. If specified when running ExtractorCLI,
    no file rejection filters will be applied by CohortExtractor. True by
    default for identifiable extractions
-   Added caching of values looked up in NLP/rulesbase for IsIdentifiable tool
-   Added new rejector that throws out values (e.g. patient IDs) whose IDs are
    stored in a database table. Set `RejectColumnInfos` option in yaml to enable
    this
-   Added a check to QueryToExecuteResult for RejectReason being null when
    Reject is true.

### Changed

-   [breaking] Environment variables are no longer required. Previous settings
    now appear in configuration file
    -   Environment variable `SMI_LOGS_ROOT` is now `GlobalOptions.LogsRoot`
    -   Environment variable `MONGO_SERVICE_PASSWORD` is now
        `MongoDbOptions.Password`
    -   Removed `ISIDENTIFIABLE_NUMTHREADS` as it didn't work correctly anyway
-   Extraction report: Group PixelData separately and sort by length
-   IsIdentifiable Reviewer 'Symbols' rule factory now supports digits only or
    characters only mode (e.g. use `\d` for digits but leave characters
    verbatim)
-   IsIdentifiable Reviewer 'symbols' option when building Regex now builds
    capture groups and matches only the failing parts of the input string not
    the full ProblemValue. For example `MR Head 12-11-20` would return
    `(\d\d-\d\d-\d\d)$`

### Fixed

-   Fix the extraction output directory to be
    `<projId>/extractions/<extractname>`

### Dependencies

-   Bump fo-dicom.Drawing from 4.0.5 to 4.0.6
-   Bump fo-dicom.NetCore from 4.0.5 to 4.0.6
-   Bump HIC.BadMedicine.Dicom from 0.0.6 to 0.0.7
-   Bump HIC.DicomTypeTranslation from 2.3.0 to 2.3.1
-   Bump HIC.FAnsiSql from 1.0.2 to 1.0.5
-   Bump HIC.RDMP.Dicom from 2.1.6 to 2.1.10
-   Bump HIC.RDMP.Plugin from 4.1.6 to 4.1.8
-   Bump HIC.RDMP.Plugin.Test from 4.1.6 to 4.1.8
-   Bump Microsoft.CodeAnalysis.CSharp.Scripting from 3.6.0 to 3.7.0
-   Bump Microsoft.Extensions.Caching.Memory from 3.1.6 to 3.1.8
-   Bump Microsoft.NET.Test.Sdk from 16.6.1 to 16.7.1
-   Bump MongoDB.Driver from 2.11.0 to 2.11.2
-   Bump System.IO.Abstractions from 12.1.1 to 12.1.9
-   Bump System.IO.Abstractions.TestingHelpers from 12.1.1 to 12.1.9
-   Bump Terminal.Gui from 0.81.0 to 0.89.4

## [1.11.1] - 2020-08-12

-   Set PublishTrimmed to false to fix bug with missing assemblies in prod.

## [1.11.0] - 2020-08-06

### Added

-   DicomDirectoryProcessor and TagReader support for zip archives
    -   Expressed in notation `/mydrive/myfolder/myzip.zip!somesubdir/my.dcm`
    -   Requires command line `-f zips`

### Changed

-   Improved the extraction report by summarising verification failures
-   Start MongoDB in replication mode in the Travis builds
-   Switch to self-contained .Net binaries to avoid dependency on host runtime
    package
-   NationalPACSAccessionNumber is now allowed to be null in all messages

### Dependencies

-   Bump HIC.RDMP.Plugin from 4.1.5 to 4.1.6
-   Bump MongoDB.Driver from 2.10.4 to 2.11.0
-   Bump System.IO.Abstractions from 12.0.10 to 12.1.1
-   Bump System.IO.Abstractions.TestingHelpers from 12.0.10 to 12.1.1
-   Bump jackson-dataformat-yaml from 2.11.1 to 2.11.2

## [1.10.0] - 2020-07-31

### Changed

-   Updated the extraction report to be more human-readable #320, #328
-   Add CLI option to CohortPackager to allow an existing report to be recreated
    #321
-   Added a runsettings file for NUnit to allow configuration of test output.
    Fixes an issue with TravisCI and NUnit3TestAdapter v3.17.0, which caused the
    test output to spill to over 20k lines.

### Dependencies

-   Bump HIC.FAnsiSql from 0.11.1 to 1.0.2
-   Bump HIC.RDMP.Dicom from 2.1.5 to 2.1.6
-   Bump HIC.RDMP.Plugin from 4.1.3 to 4.1.5
-   Bump Magick.NET-Q16-AnyCPU from 7.20.0 to 7.21.1
-   Bump Microsoft.Extensions.Caching.Memory from 3.1.5 to 3.1.6
-   Bump System.IO.Abstractions from 12.0.1 to 12.0.10
-   Bump System.IO.Abstractions from 12.0.1 to 12.0.2
-   Bump System.IO.Abstractions.TestingHelpers from 12.0.1 to 12.0.2
-   Bump com.fasterxml.jackson.dataformat.jackson-dataformat-yaml from 2.11.0 to
    2.11.1
-   Bump org.mockito.mockito-core from 3.3.3 to 3.4.6

## [1.9.0] - 2020-06-22

### Added

-   Added image extraction blacklist rejector.
    -   Configure with `Blacklists` option (specify a list of Catalogue IDs)
    -   Catalogues listed must include one or more column(s) StudyInstanceUID,
        SeriesInstanceUID, SOPInstanceUID.
    -   Records in the referenced table will blacklist where any UID is found
        (StudyInstanceUID, SeriesInstanceUID or SOPInstanceUID). This allows
        blacklisting an entire study or only specific images.
    -   [breaking] Config on live system may need updated
-   Change the extraction directory generation to be
    `<projname>/image-requests/<extractname>`. Fixes
    [MVP Service #159](https://dev.azure.com/smiops/MVP%20Service/_workitems/edit/159/)

### Fixed

-   Fixed IsIdentifiable rule order being the order the files are detected in
    rules directory (Now goes IgnoreRules=>ReportRules=>SocketRules)
-   Adjust log handling in CTP anonymiser to use SMIlogging setup
-   IsIdentifiable case-sensitive rules now implemented with property
-   Bufix for fo-dicom image handling race condition in Release mode builds
    (issue #238)

### Changed

-   Refactored `WhiteListRule` to inherit from `IsIdentifiableRule` (affects
    serialization).
    -   Parent property `As` replaces `IfClassification`
    -   `CaseSensitive` replaces `IfPatternCaseSensitive` and
        `IfPartPatternCaseSensitive` (Also fixes serialization bug)
-   Bump CommandLineParser from 2.7.82 to 2.8.0
-   Bump CsvHelper from 15.0.4 to 15.0.5
-   Bump HIC.BadMedicine.Dicom from 0.0.5 to 0.0.6
-   Bump HIC.DicomTypeTranslation from 2.2.2 to 2.3.0
-   Bump HIC.RDMP.Dicom from 2.0.9 to 2.1.5
-   Bump HIC.RDMP.Plugin from 4.0.2 to 4.1.3
-   Bump Magick.NET-Q16-AnyCPU from 7.16.0 to 7.20.0
-   Bump Microsoft.CodeAnalysis.CSharp.Scripting from 3.5.0 to 3.6.0
-   Bump Microsoft.Extensions.Caching.Memory from 3.1.3 to 3.1.5
-   Bump MongoDB.Driver from 2.10.3 to 2.10.4
-   Bump StackExchange.Redis from 2.1.30 to 2.1.58
-   Bump System.IO.Abstractions from 10.0.8 to 12.0.1
-   Bump YamlDotNet from 8.1.0 to 8.1.2
-   Bump fo-dicom.Drawing from 4.0.4 to 4.0.5
-   Pinned fo-dicom.NetCore to 4.0.5

## [1.8.1] - 2020-04-17

### Fixed

-   Fix null check bug in CohortPackager when no files match the extraction
    filter

## [1.8.0] - 2020-04-16

### Added

-   Added Terminal.Gui at version 0.81.0
-   Added data/IsIdentifiableRules

### Changed

-   \[Breaking\] Promote the PT modality to its own collection in MongoDB
-   \[Breaking\] Renamed `RedisHost` to `RedisConnectionString` in the config
    options for clarity
-   Update to .Net Core 3.1 (supported until Dec 2022) since 2.2 support ended
    last year
-   Switch CohortExtractor to use batched message producers
-   Simplify the Travis build script
-   Fail any integration tests in CI if a required service is not available
    (instead of skipping)
-   Specified LangVersion 8.0 in all project files
-   Upgraded CommandLineParser from 2.5.0 to 2.7.82
-   Upgraded CsvHelper from 12.1.2 to 15.0.4
-   Upgraded HIC.Rdmp.Dicom from 2.0.8 to 2.0.9
-   Upgraded JetBrains.Annotations from 2019.1.3 to 2020.1.0
-   Upgraded Magick.NET-Q16-AnyCPU from 7.15.1 to 7.16.0
-   Upgraded Microsoft.CodeAnalysis.CSharp.Scripting from 3.5.0-beta2-final to
    3.5.0
-   Upgraded MongoDB.Driver from 2.9.3 to 2.10.3
-   Upgraded StackExchange.Redis from 2.0.601 to 2.1.30
-   Upgraded System.Drawing.Common from 4.6.0 to 4.7.0
-   Upgraded System.IO.Abstractions from 7.0.7 to 10.0.8
-   Upgraded YamlDotNet from 6.0.0 to 8.1.0

### Fixed

-   Fixed logging to directories in the Java services

## [1.7.0] - 2020-03-30

### Added

-   Added undo feature to IsIdentifiableReviewer
-   Java microservices now log to SMI_LOGS_ROOT

### Changed

-   Upgraded HIC.DicomTypeTranslation from `2.1.2` to `2.2.0`
    -   This includes an upgrade to fo-dicom from `4.0.1` to `4.0.4`
-   Upgraded fo-dicom.Drawing from `4.0.1` to `4.0.4`
-   Upgraded HIC.RdmpDicom from `2.0.7` to `2.0.8`

## [1.6.0] - 2020-03-17

### Changed

-   Update CohortPackager for new extraction design

    -   Consume messages from CTP (failed anonymisation) and IsIdentifiable
        (verification)
    -   Add support for extraction by modality
    -   Remove the final check for the anonymised file. IsIdentifiable handles
        this already
    -   Refactor tests

-   Start to refactor core RabbitMqAdapter code to allow unit testing

## [1.5.2] - 2020-03-12

### Added

-   IsIdentifiableReviewer considers rule capture groups when performing
    redactions (e.g. can now handle custom rules like `^(Ninewells)`)
-   IsIdentifiableReviewer adds comment with time/user to rules file e.g.
    `#TZNind - 3/10/2020 1:17:17 PM`
-   IsIdentifiableReviewer checks custom patterns match the original Failure
-   IsIdentifiable microservice was started with --service but can now be
    started with the service verb allowing it to take additional options. It
    should now be started with `service -y file.yaml`
-   IsIdentifiable no longer reads Rules.yaml from the current directory. It now
    has a command line option --RulesDirectory, to go with the already existing
    --RulesFile. That will read all \*.yaml files in the given directory.
    However when run as a microservice the yaml file specifies a DataDirectory;
    the RulesDirectory will implicitly be a subdirectory called
    IsIdentifiableRules from which all \*.yaml files will be read.

### Changed

-   IsIdentifiableReviewer now tries to isolate 'Problem Words' when generating
    it's suggested Updater Regex rules (e.g. now suggests `^Ninewells` instead
    of `^Ninewells\ Spike\ CT$`.)

## [1.5.1] - 2020-03-06

-   Improved usability of IsIdentifiableReviewer

## [1.5.0] - 2020-03-05

-   \[Breaking\] Updated RabbitMQ extraction config to match extraction plan v2
-   Refactor Java exception handling and use of threads
-   `TessDirectory` option in [IsIdentifiable] now expects tesseract models file
    to exist (no longer downloads it on demand)
-   Added support for outsourcing classification (e.g. NLP) to other processes
    via TCP (entered in [SocketRules] in `Rules.yaml`)
-   IsIdentifiable NLP text classification now outsourced via TCP to any
    services configured in
    -   [StanfordNER implementation written in java](./src/microservices/uk.ac.dundee.hic.nerd/README.md)
-   New CohortExtractor yaml config option `ProjectPathResolverType` which
    determines the folder structure for extracted images
-   Added [script](./utils/rabbitmq-config-tester/rabbitmq-config-tester.py) to
    verify RabbitMQ config files
-   Added `DynamicRejector` which takes its cohort extraction rules from a
    script file (of CSharp code)
-   Added new application for reviewing IsIdentifiable output files

### Fixed

-   Corrected the GetHashCode implementation in the MessageHeader class

## [1.4.5] - 2020-02-26

-   Add clean shutdown hook for IdentifierMapper to clean up the worker threads

## [1.4.4] - 2020-02-25

-   Update Travis config and Java library install shell script to resolve some
    Travis stability issues
-   Adjust batching so workers queue replies/acks while a worker thread commits
    those asynchronously, allowing elastic batch sizes (qosprefetch setting now
    controls maximum batch size, parallelism capped at 50)

## [1.4.3] - 2020-02-21

### Changed

-   Batch up RabbitMQ messages/acks in IdentifierMapper to avoid contention with
    the message publishing persistence

## [1.4.2] - 2020-02-18

### Added

-   Added unit test for AccessionDirectoryLister as part of
    DicomDirectoryProcessor tests

### Changed

-   Make performance counters in RedisSwapper atomic for thread-safety
-   Clean up threads when using threaded mode in RabbitMQAdapter
-   Use explicit threads rather than Task queueing in IdentifierMapper

## [1.4.1] - 2020-02-17

### Added

-   Added randomisation in the retry delay on DicomRelationalMapper (and set
    minimum wait duration to 10s)

### Fixed

-   Fixed DLE Payload state being wrong when retrying batches (when it is half /
    completely consumed)
-   Added lock on producer sending messages in IdentifierMapper

## [1.4.0] - 2020-02-14

### Added

-   Added in memory caching of the last 1024 values when using Redis wrapper for
    an IdentifierSwapper
-   Added some parallelism and marshalling of backend queries to improve
    throughput in IdentifierSwapper
-   Added temporary flag for RabbitMQAdapter parallelism for the above. Only
    enabled for the IdentifierMapper for now
-   Added new mode to DicomDirectoryProcessor which allows reading in a list of
    accession directories

## [1.3.1] - 2020-02-13

### Changed

-   Pinned fo-dicom to v4.0.1

## [1.3.0] - 2020-02-06

### Added

-   Added (optional) DicomFileSize property to ETL pipeline. Add to template(s)
    with:

```yaml
- ColumnName: DicomFileSize
  AllowNulls: true
  Type:
      CSharpType: System.Int64
```

-   Added new microservice IsIdentifiable which scans for personally
    identifiable information (in databases and dicom files)
-   Added support for custom rules in IsIdentifiable (entered in `Rules.yaml`)
    -   Rules are applied in the order they appear in this file
    -   Rules are applied before any other classifiers (i.e. to allow
        whitelisting rules)
-   Added `RedisSwapper` which caches answers from any other swapper. Set
    `RedisHost` option in yaml to use.

### Changed

-   Updated RDMP and Dicom plugins
-   Refactor Java exception handling and use of threads

## [1.2.3] - 2020-01-09

### Changed

-   RabbitMQAdapter: Improve handling of timeouts on connection startup

### Added

-   Improved logging in IdentifierSwappers

### Changed

-   Guid swapper no longer limits input identifiers to a maximum of 10
    characters

### Fixed

-   Fixed DicomRelationalMapper not cleaning up STAGING table remnants from
    previously failed loads (leading to crash)

## [1.2.2] - 2020-01-08

### Fixed

-   RAW to STAGING migration now lists columns explicitly (previously used
    `SELECT *` which could cause problems if RAW and STAGING column orders
    somehow differed)

## [1.2.1] - 2020-01-06

### Added

-   Added the `set-sleep-time-ms` control message to DicomReprocessor

### Changed

-   Updated Rdmp.Dicom nuget package to 2.0.6

## [1.2.0] - 2019-12-12

### Added

-   Improved travis deployment
-   (Re-)added Smi.NLog.config in builds
-   Added better CLI argument descriptions for DicomReprocessor
-   Added error logging for RabbitMQ bad Ack responses
    -   Previously: `BasicReturn for TEST.IdentifiableImageExchange`
    -   Now :
        `BasicReturn for Exchange 'TEST.IdentifiableImageExchange' Routing Key 'reprocessed' ReplyCode '312' (NO_ROUTE)`
-   Added new swapper `TableLookupWithGuidFallbackSwapper` which performs lookup
    substitutions but allocates guids for lookup misses
-   Added Travis CI build & deploy for all services

### Changed

-   Make exceptions on startup clearer
-   Updated to latest RDMP API (4.0.1)
-   `TableLookupSwapper` now throws consistent error if the provided table does
    not exist during `Setup` (previously it would error with DBMS specific error
    message at lookup time)

### Fixed

-   Fixed freeze condition when exchanges are not mapped to queues
-   IdentifierMapper now loads all FAnsi database implementations up front on
    startup

## [1.1.0] - 2019-11-22

### Added

-   Improvements to unit and integration tests
-   Documentation fixes
-   Config file for Dependabot
-   Test for DicomFile SkipLargeTags option. Closes
    [#19](https://dev.azure.com/SmiOps/MVP%20Service/_workitems/edit/19)

### Changed

### C\# dependencies

-   Bumped HIC.DicomTypeTranslation from 1.0.0.3 to 2.1.2
-   Bumped HIC.RDMP.Plugin from 3.1.1 to 4.0.1-rc2
-   Bumped Newtonsoft.Json from 12.0.2 to 12.0.3
-   Bumped RabbitMQ.Client from 5.1.0 to 5.1.2
-   Bumped System.IO.Abstractions from 4.2.17 to 7.0.7
-   Bumped MongoDB.Driver from 2.8.0 to 2.9.3

### Java dependencies

-   Bumped jackson-databind from 2.9.6 to 2.9.10.0

## [1.0.0] - 2019-11-18

First stable release after importing the repository from the private
[SMIPlugin](https://github.com/HicServices/SMIPlugin) repo.

### Added

-   ForGuidIdentifierSwapper automatically creates it's mapping database if it
    does not exist on the server referenced (previously only table was
    automatically created)

### Changed

-   Updated to
    [Rdmp.Dicom 2.0.2](https://github.com/HicServices/RdmpDicom/blob/master/CHANGELOG.md#202-2019-11-13)
-   Updated to
    [Rdmp.Core 3.2.1](https://github.com/HicServices/RDMP/blob/develop/CHANGELOG.md#321---2019-10-30)

### Removed

-   Anonymous `MappingTableName` must now be fully specified to pass validation
    (e.g. `mydb.mytbl`). Previously skipping database portion was supported.

[unreleased]: https://github.com/SMI/SmiServices/compare/v4.0.0...master
[4.0.0]: https://github.com/SMI/SmiServices/compare/v3.2.1...v4.0.0
[3.2.1]: https://github.com/SMI/SmiServices/compare/v3.2.0...v3.2.1
[3.2.0]: https://github.com/SMI/SmiServices/compare/v3.1.0...v3.2.0
[3.1.0]: https://github.com/SMI/SmiServices/compare/v3.0.2...v3.1.0
[3.0.2]: https://github.com/SMI/SmiServices/compare/v3.0.1...v3.0.2
[3.0.1]: https://github.com/SMI/SmiServices/compare/v3.0.0...v3.0.1
[3.0.0]: https://github.com/SMI/SmiServices/compare/v2.1.1...v3.0.0
[2.1.1]: https://github.com/SMI/SmiServices/compare/v2.1.0...v2.1.1
[2.1.0]: https://github.com/SMI/SmiServices/compare/v2.0.0...v2.1.0
[2.0.0]: https://github.com/SMI/SmiServices/compare/v1.15.1...v2.0.0
[1.15.1]: https://github.com/SMI/SmiServices/compare/v1.15.0...v1.15.1
[1.15.0]: https://github.com/SMI/SmiServices/compare/v1.14.1...v1.15.0
[1.14.1]: https://github.com/SMI/SmiServices/compare/v1.14.0...v1.14.1
[1.14.0]: https://github.com/SMI/SmiServices/compare/v1.13.0...v1.14.0
[1.13.0]: https://github.com/SMI/SmiServices/compare/v1.12.2...v1.13.0
[1.12.2]: https://github.com/SMI/SmiServices/compare/v1.12.1...v1.12.2
[1.12.1]: https://github.com/SMI/SmiServices/compare/v1.12.0...v1.12.1
[1.12.0]: https://github.com/SMI/SmiServices/compare/v1.11.1...v1.12.0
[1.11.1]: https://github.com/SMI/SmiServices/compare/v1.11.0...v1.11.1
[1.11.0]: https://github.com/SMI/SmiServices/compare/v1.10.0...v1.11.0
[1.10.0]: https://github.com/SMI/SmiServices/compare/v1.9.0...v1.10.0
[1.9.0]: https://github.com/SMI/SmiServices/compare/v1.8.1...v1.9.0
[1.8.1]: https://github.com/SMI/SmiServices/compare/v1.8.0...v1.8.1
[1.8.0]: https://github.com/SMI/SmiServices/compare/v1.7.0...v1.8.0
[1.7.0]: https://github.com/SMI/SmiServices/compare/v1.6.0...v1.7.0
[1.6.0]: https://github.com/SMI/SmiServices/compare/v1.5.2...v1.6.0
[1.5.2]: https://github.com/SMI/SmiServices/compare/v1.5.1...v1.5.2
[1.5.1]: https://github.com/SMI/SmiServices/compare/v1.5.0...v1.5.1
[1.5.0]: https://github.com/SMI/SmiServices/compare/v1.4.5...v1.5.0
[1.4.5]: https://github.com/SMI/SmiServices/compare/v1.4.4...v1.4.5
[1.4.4]: https://github.com/SMI/SmiServices/compare/v1.4.3...v1.4.4
[1.4.3]: https://github.com/SMI/SmiServices/compare/v1.4.2...v1.4.3
[1.4.2]: https://github.com/SMI/SmiServices/compare/v1.4.1...v1.4.2
[1.4.1]: https://github.com/SMI/SmiServices/compare/v1.4.0...v1.4.1
[1.4.0]: https://github.com/SMI/SmiServices/compare/v1.3.1...v1.4.0
[1.3.1]: https://github.com/SMI/SmiServices/compare/v1.3.0...v1.3.1
[1.3.0]: https://github.com/SMI/SmiServices/compare/v1.2.3...v1.3.0
[1.2.3]: https://github.com/SMI/SmiServices/compare/v1.2.2...v1.2.3
[1.2.2]: https://github.com/SMI/SmiServices/compare/v1.2.1...v1.2.2
[1.2.1]: https://github.com/SMI/SmiServices/compare/1.2.0...v1.2.1
[1.2.0]: https://github.com/SMI/SmiServices/compare/1.1.0...1.2.0
[1.1.0]: https://github.com/SMI/SmiServices/compare/1.0.0...1.1.0
[1.0.0]: https://github.com/SMI/SmiServices/releases/tag/1.0.0
[isidentifiable]: ./src/microservices/Microservices.IsIdentifiable/README.md
[socketrules]:
    ./src/microservices/Microservices.IsIdentifiable/README.md#socket-rules

<img src="https://avatars2.githubusercontent.com/u/56437605?s=200"/>

# Contributors to SMI Services

## [EPCC](https://github.com/EPCCed)

- [Andrew Brooks](https://github.com/howff)
- [Ally Hume](https://github.com/allyhume)
- [Ruairidh MacLeod](https://github.com/rkm)
- [Bianca Prodan](https://github.com/2bPro)

### Former

- [Paul Graham](https://github.com/pjgraham)

## [HIC](https://github.com/HicServices)

- [Douglas Hardy](https://github.com/dhardy-hic)
- [Thomas Nind](https://github.com/tznind)
- [Leandro Tramma](https://github.com/Tallmaris)
- [James Sutherland](https://github.com/jas88)

### Former

- [Trevor Carpenter](https://github.com/trevor-carpenter)
- [?](https://github.com/maabdelatif)

## Independent


# Packages Used

### Risk Assessment common to all:
1. Packages on NuGet are virus scanned by the NuGet site.
2. This package is widely used and is actively maintained.
3. It is open source.

| Package | Source Code |  License | Purpose
| ------- | ------------| ------- | ------- | 
| CommandLineParser | [GitHub](https://github.com/commandlineparser/commandline) | [MIT](https://opensource.org/licenses/MIT)| Command line argument parsing |
| CsvHelper | [GitHub](https://github.com/JoshClose/CsvHelper) | [MS-PL and Apache 2.0](https://github.com/JoshClose/CsvHelper/blob/master/LICENSE.txt)| Writing reports out to CSV reports |
| Equ | [GitHub](https://github.com/thedmi/Equ) | [2.3.0](https://www.nuget.org/packages/Equ/2.3.0) | [MIT](https://opensource.org/licenses/MIT) | Automatic equality functions |
| fo-dicom.NetCore | [GitHub](https://github.com/fo-dicom/fo-dicom) | [MS-PL](https://opensource.org/licenses/MS-PL) | |
| HIC.DicomTypeTranslation | [GitHub](https://github.com/HicServices/DicomTypeTranslation) | [GPL 3.0](https://www.gnu.org/licenses/gpl-3.0.html) | Translate dicom types into C# / database types |
| HIC.FAnsiSql | [GitHub](https://github.com/HicServices/FansiSql) | [GPL 3.0](https://www.gnu.org/licenses/gpl-3.0.html) | Database abstraction layer |
| HIC.RDMP.Dicom | [GitHub](https://github.com/HicServices/RdmpDicom) | [GPL 3.0](https://www.gnu.org/licenses/gpl-3.0.html) | RDMP Plugin containing data load / pipeline components for imaging, reading dicom files etc |
| HIC.RDMP.Plugin | [GitHub](https://github.com/HicServices/RDMP) | [GPL 3.0](https://www.gnu.org/licenses/gpl-3.0.html) | Interact with RDMP objects, base classes for plugin components etc |
| JetBrains.Annotations | |[MIT](https://opensource.org/licenses/MIT) | Static analysis tool |
| Magick.NET-Q16-AnyCPU | [GitHub](https://github.com/dlemstra/Magick.NET) | [Apache License v2](https://github.com/dlemstra/Magick.NET/blob/master/License.txt) | The .NET library for [ImageMagick](https://imagemagick.org/index.php) |
| Microsoft.CodeAnalysis.CSharp.Scripting | [GitHub](https://github.com/dotnet/roslyn) | [MIT](https://opensource.org/licenses/MIT)  | Supports dynamic rules for cohort extraction logic |
| Microsoft.Extensions.Caching.Memory | [GitHub](https://github.com/dotnet/extensions) | [Apache 2.0](https://www.nuget.org/packages/Microsoft.Extensions.Caching.Memory/3.1.7/License) | Caching ID mappings retrieved from Redis/MySQL
| NLog | [GitHub](https://github.com/NLog/NLog) | [BSD 3-Clause](https://github.com/NLog/NLog/blob/dev/LICENSE.txt) | Flexible user configurable logging |
| Newtonsoft.Json | [GitHub](https://github.com/JamesNK/Newtonsoft.Json) | [MIT](https://opensource.org/licenses/MIT) | Serialization of objects for sharing/transmission
| RabbitMQ.Client | [GitHub](https://github.com/rabbitmq/rabbitmq-dotnet-client) | [Apache License v2 / MPL 1.1](https://github.com/rabbitmq/rabbitmq-dotnet-client/blob/master/LICENSE) | Handles messaging between microservices |
| SecurityCodeScan.VS2019 | [GitHub](https://security-code-scan.github.io/) | [LGPL 3.0](https://opensource.org/licenses/lgpl-3.0.html) | Scans code for security issues during build |
| StackExchange.Redis | [GitHub](https://github.com/StackExchange/StackExchange.Redis) |[MIT](https://opensource.org/licenses/MIT) | Required for RedisSwapper |
| Stanford.NLP.CoreNLP | [GitHub Pages](https://sergey-tihon.github.io/Stanford.NLP.NET/) | [GNU v2](https://github.com/sergey-tihon/Stanford.NLP.NET/blob/master/LICENSE.txt)| Name / Organisation detection in text |
| System.Drawing.Common | [GitHub](https://github.com/dotnet/corefx) | [MIT](https://opensource.org/licenses/MIT)  | Supports reading pixel data |
| System.IO.Abstractions | [GitHub](https://github.com/System-IO-Abstractions/System.IO.Abstractions) | [MIT](https://opensource.org/licenses/MIT) | Makes file system injectable in tests |
| System.IO.FileSystem | [GitHub](https://github.com/dotnet/corefx) |[MIT](https://opensource.org/licenses/MIT)  | File I/O |
| Terminal.Gui | [GitHub](https://github.com/migueldeicaza/gui.cs/) |[MIT](https://opensource.org/licenses/MIT) | Console GUI library |
| Tesseract | [GitHub](https://github.com/charlesw/tesseract/) |[Apache License v2](https://github.com/charlesw/tesseract/blob/master/LICENSE.txt)  | Optical Character Recognition in Dicom Pixel data|
| YamlDotNet | [GitHub](https://github.com/aaubry/YamlDotNet)  | [MIT](https://opensource.org/licenses/MIT) |Loading configuration files
| fo-dicom.Drawing | [GitHub](https://github.com/fo-dicom/fo-dicom) | [MS-PL](https://opensource.org/licenses/MS-PL)| Support library for reading DICOM pixel data |
| coveralls.io | [GitHub](https://github.com/coveralls-net/coveralls.net) | [GNU](https://github.com/coveralls-net/coveralls.net#license)| Uploader for dot net coverage reports to Coveralls.io |
| OpenCover | [GitHub](https://github.com/OpenCover/opencover) |[MIT Compatible](https://github.com/OpenCover/opencover/blob/master/LICENSE)  | Calculates code coverage for tests|
## CI Status

[![Coverage Status](https://coveralls.io/repos/github/SMI/SmiServices/badge.svg)](https://coveralls.io/github/SMI/SmiServices) ![GitHub](https://img.shields.io/github/license/SMI/SmiServices) [![Total alerts](https://img.shields.io/lgtm/alerts/g/SMI/SmiServices.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/SMI/SmiServices/alerts/)

| OS      | CI Pipeline |
| :---        |    :----:   |
| Windows      | [![Build Status](https://dev.azure.com/SmiOps/Public/_apis/build/status/SmiServices%20Windows?branchName=master)](https://dev.azure.com/SmiOps/Public/_build/latest?definitionId=4&branchName=master)       |
| Linux   | [![Build Status](https://dev.azure.com/SmiOps/Public/_apis/build/status/SmiServices%20Linux?branchName=master)](https://dev.azure.com/SmiOps/Public/_build/latest?definitionId=3&branchName=master) |


Version: `4.0.0`

# SMI Services

## Contents

1. [Introduction](#10-introduction)
   1. [Overview](#11-overview)
   2. [Glossary or Terminology](#12-glossary-or-terminology)
   3. [Background and Context](#13-background-and-context)
   4. [Goals and Technical Requirements](#14-goals-and-technical-requirements)
   5. [Out of Scope](#15-out-of-scope)
   6. [Assumptions](#16-assumptions)
2. [Solutions](#20-solutions)
   1. [Deliverable Solution / Design](#21-deliverable-solution--design)
   1. [Test Plan](#22-test-plan)
   1. [Monitoring and Alerting Plan](#23-monitoring-and-alerting-plan)
   1. [Release / Roll-out and Deployment Plan](#24-release--roll-out-and-deployment-plan)
   1. [Rollback Plan](#25-rollback-plan)
   1. [Associated Documentation](#26-associated-documentation)
1. [Further Considerations](#30-further-considerations)
   1. [Human Resource Requirements](#31-human-resource-requirements)
   2. [Security Considerations](#32-security-considerations)
   1. [Legal and Ethical considerations](#33-legal-and-ethical-considerations)
   1. [Accessibility considerations](#34-accessibility-considerations)
   1. [Support considerations](#35-support-considerations)
1. [Future Work](#40-future-work)
   1. [Residual work estimates and timelines](#41-residual-work-estimates-and-timelines)
   1. [Open Questions](#42-open-questions)
1. [End Matter](#50-end-matter)
   1. [Related Work](#51-related-work)
   1. [Acknowledgments](#52-acknowledgments)
1. [Microservices](#microservices)
   1. [Data Load Microservices](#data-load-microservices)
   1. [Image Extraction Microservices](#image-extraction-microservices)
1. [.sln Overivew](#sln-overview)
1. [Building](#building)
1. [Running](#running)
1. [Dependencies](#dependencies)
1. [Scalability](#scalability)


## 1.0 Introduction

### 1.1 Overview

SMI Services is a suite of tools designed to deliver scalable dicom image indexing for cohort building and extraction in anonymised sub sets (e.g. for research).  It is an Extract, Transform and Load tool (ETL) for imaging data.

![loaddiagram](./docs/Images/SmiFlow.svg)

The problem addressed is how to enable linking of dicom metadata with other clinical data (e.g. electronic health records - EHR).  The context in which it was developed is the loading and anonymisation of metadata for 10 years of Scottish national clinical imaging data (2 petabytes).

The following processes are solved by the suite:

- Robust parallel loading of dicom metadata into a relational database where it can be linked to EHR data
- Anonymisation of the dicom metadata within the relational database
- Cohort building by linking the relational database tables with other EHR data (e.g. biochemistry results, prescriptions etc)
- Producing anonymous dicom image files for a subset of the repository

Stakeholders likely to interact with this suite include Research Coordinators (building research extracts) and Data Analysts (loading data, verifying anonymisation etc)

### 1.2 Glossary or Terminology

For RDMP terms see the [RDMP Glossary](https://github.com/HicServices/RDMP/blob/develop/Documentation/CodeTutorials/Glossary.md)

For DICOM specific terms see the [DICOM tag browser](https://dicom.innolitics.com/ciods) or [DICOM specification]

Dicom tags are key/value pairs that assist in understanding the provenance of a dicom image.  They contain information such as PatientID, PatientAge, StudyDescription etc.  The [DICOM specification] describes what tags are required and which are optional.  Each tag has a specific datatype (date, text, float etc).  The specification supports tree structure tags called Sequences.  A Sequence can contain subtags which are also Sequences resulting in a tree data structure.  Tags can also support multiple elements (i.e. array datatype), this is called multiplicity.

### 1.3 Background and Context

Historically dicom images are held in a clinical PACS or in an imaging informatics platform such as [XNAT].  SmiServices can function standalone or side by side with such tools.  The unique features of SmiServices are it's ability to present large imaging datasets as indexed flat tables in a relational database (MySql, Sql Server, Postgres or Oracle) where they can be linked with cohorts/EHR datasets. This is a worthy addition since it allows for cohort building and extraction using existing tools that data analysts are familiar with (R, Sql scripts, [RDMP Cohort Builder](https://github.com/HicServices/RDMP) etc).

### 1.4 Goals and Technical Requirements

The goals are to load dicom metadata, build cohorts and extract anonymous image subsets.

This requires dotnet, RabbitMQ, RDMP, MongoDb and a Relational Database.  For more info on setting up SmiServices see [Deployment](#deployment).

SmiServices benefits from:

- Running on a cluster (many VMs running many copies of each service)
- A Parallel File System (e.g. [BeeGFS](https://en.wikipedia.org/wiki/BeeGFS) or [Lustre](https://www.lustre.org/))

### 1.5 Out of Scope

SmiServices does not support imaging workflows (e.g. running image algorithms).

It also does not have an API for external communication (e.g. Dicom Query Retrieve or FHIR/HL7 etc).  The imaging metadata produced by SmiServices can be queried using MongoDb queries or SQL.

Once the image metadata is in the relational database then cohorts can be created using standard cohort building tools (e.g. RDMP Cohort Builder, R, SQL).  The specifics of cohort building are not covered in this document as that is covered elsewhere (in the documentation of each tool).

### 1.6 Assumptions

SmiServices assumes that database servers are optimised and properly resourced to store the volume of image metadata anticipated.  ETL is robust and can deal with outages and database stability issues (e.g. lock collisions) but these can errode system performance.

The solution is designed for large image collections (i.e. billions of images).  It supports flexible schema definitions such that only the tags required for cohort building (and the image file paths) are loaded.  Therefore successful usage of the tool requires a basic understanding of dicom tag significance and an appropriately large body of images to justify its use.

## 2.0 Solutions

### 2.1 Deliverable Solution / Design

Microservices is a design in which each component is self contained and has a single responsibility.  Decoupling processes allows for easier maintenance and testability.  Any number of copies of each service can be run at the same time allowing scalability through parallel processing.

Communication between services is through RabbitMQ.   RabbitMQ is one of the most popular open source message brokers, supporting both high-scale and high-availability requirements. 

Data is promoted sequentially between places (file system, databases etc).  Each promotion is carried out by one or more service types in a chain.  For example promoting the dicom tag data to the MongoDb document store involves the [ProcessDirectory], [DicomTagReader] and [MongoDBPopulator] services.  Since only tag data is promoted, there are no excessive storage overheads associated with data existing in multiple places at once within the system (duplication).  Each place supports specific user processes (see below) and the technology chosen is based on it's suitability for those processes:

![loaddiagram](./docs/Images/processes.svg)


| User Process      | Description | Application(s) |
| ----------- | ----------- |----------- |
| ETL Scheduling | Users uses command line tool to pick specific directories of images for loading to MongoDb |  [ProcessDirectory] |
| Data Exploration   | Users explore full dicom metadata to identify anonymisation requirements, useful tags for cohort building and trigger promotion of subsets of images to the relational database (e.g. Modality CT for 20010) | mongo command line, MongoDB Compass etc and [DicomReprocessor] |
| Anonymisation Verification | Users review reports on the effectiveness of the MongoDb->Relational Database anonymisation pipeline and create new anonymisation/redaction rules where the system has failed to correctly redact the tag data | [IsIdentifiableReviewer] |
| Cohort Building | Users view link EHR data with the relational database data to produce extractable image subsets of the archive that fit specific research projects (e.g. CT head scans between 2015-2018 where the patient had 3+ prescriptions for drug X) | R, SQL, [RDMP], [ExtractImages] |

MongoDb was chosen for the 'Data Exploration' process because it is a 'document store' - which means it is able to store the full tree structure of dicom metadata.  MongoDb improves access speed to the metadata tags and allows for aggregation and search activities. MongoDb is designed to scale well and supports sharding and replication for when the number of dicom files held grows beyond single db instance capabilities.

The slowest part of ETL is reading dicom files from disk (even with a parallel file system).  The use of MongoDb allows this to be done only once per image regardless of how it is subsequently processed downstream.  Pixel data tags are not loaded to MongoDb since these would inflate storage and processing requirements and are less useful to for the task of Data Exploration.

Relational Databases were chosen for the Anonymous Tag Store since this is the format most commonly used with EHR datasets and cohort building tools.  The ability to directly link anonymous dicom metadata to research study lists and other EHR datasets held in relational databases is a core design requirement of SmiServices.  

The data load service wraps the [RDMP] data load engine and so supports MySql, Sql Server, Oracle and Postgres.  In addition it allows tailoring how corrupt/duplicate data is loaded.  Using an [RDMP] load also allows for tailoring the load process after deployment and for expanding upon the cohort building table schema(s) over time as new tags are identified as useful for cohort building.

Implicit in the [DICOM Specification] is the hierarchical layout of tags (Patient, Study, Series, Image).  Each dicom image file contains the complete tag list for its study/series.  This means that tag data is replicated in each file and only the Image level tags are uniqiue.  The use of [RDMP] and it's imaging plugin [RDMP Dicom] allows for automatic aggregation of these Patient/Study/Series level tags.  This results in a far smaller table for cohort building which improves query performance where there is not a need to query image level tags.

The exact tags required for the Cohort Building process will vary over time and may not be known at the outset.  For this reason the ETL microservices support any schema so long as the table column names match a known DICOM tag.  To assist in schema building, a standalone application [Dicom Template Builder] was created.

In order to protect patient privacy, all tag data in the relational database should be anonymous.  This is supported by several services and design choices:

- The Cohort building schema you create should contain only tags that have low volumes of identifiable data (e.g. StudyDescription but not PatientName )
- Identifiable data can be detected using the [IsIdentifiable] tool
- Identifiable data can be summarised and redacted using the [IsIdentifiableReviewer] tool
- PatientID can be substituted for an anonymous identifier with the ETL service [IdentifierMapper].

The user process of ensuring these steps have been undertaken correctly is called 'Anonymisation Verification'.

Error recovery is handled through RabbitMQ.  When a service fails to acknowledge the successful processing of a message it is automatically requeued and sent to a different service.  If a message cannot be processed after several attempts it is sent to a 'dead letter queue' where it can be evaluated later.  Services subscribe to a 'Control queue' which allows for safe shutdown of the service e.g. for system maintenance.

The most error prone section of ETL is entry to the relational database which is where primary key collisions and corrupt data must be reconciled (e.g. 2 images in the same study containing conflicting definitions of StudyDescription).  [RDMP Dicom] contains several modules designed to mitigate these issues, more information about these can be found in the [RDMP Dicom data load documentation](https://github.com/HicServices/RdmpDicom/blob/develop/Documentation/DataLoad.md).

The use of microservices not only ensures scalability and error recovery but also provides a degree of future proofing.  If a requirement emerges for a new step in ETL that cannot be handled by RDMP then a new microservice can be slotted into the load chain.  Additionally if a step is not needed for a given deployment it can be cut out (e.g. removing the [IdentifierMapper] step of ETL) by editing the SmiServices configuration files.

Configuration of the services comes from three places:  

- The command arguments given to the service on startup
- A [YAML configuration file](./data/microserviceConfigs/default.yaml)
- [RDMP]

### 2.2 Test Plan

SmiServices and [RDMP] contain both automated unit and integration tests.  These tests are automatically run on each code commit to the repository.  New features are written in a 'pull request' which is independently tested and approved.  Pull requests can be from developers working on the project (branches) or from external collaborators (forks).

User testing of the services can be done using the [Docker Image](https://github.com/jas88/smideploy).

A tool ([BadMedicine.Dicom]) has been created that generates synthetic test DICOM images.  [BadMedicine.Dicom] can be used to generate images for testing the service.  It includes support for generating placeholder pixel data so that file size can be modelled.  This helps with non functional testing of hardware when architecting an SmiServices deployment.  

Cohort building can be tested by generating synthetic EHR data with the sibling tool [BadMedicine].  When used with the same seed as [BadMedicine.Dicom] relational database tables or CSV files of synthetic medical can be generated (e.g. biochemistry, prescribing, demography).

Code coverage metrics are collected and hosted on [Coveralls](https://coveralls.io/github/SMI/SmiServices).  This shows what proportion of lines of code in the codebase are covered by automated testing.  It also allows visualisation of which parts of the codebase have less coverage.  This ensures that no complex areas of suite are untested.  Automated alerts are generated on pull requests that substantially lower the code coverage (add a lot of code without tests).

Static analysis of the codebase is performed with [LGTM](https://lgtm.com/projects/g/SMI/SmiServices/alerts/).  This identifies common coding errors such as missing null checks.  Automated alerts are generated on pull requests in which such errors are detected.

Manually testing and debugging integration/unit tests requires having the relevant tool dependencies installed.  Tests are decorated with an attribute that indicates which dependencies(if any) are required.  These include:

- RequiresRelationalDb (Microsoft Sql Server / MySql)
- RequiresMongoDb (MongoDb)
- RequiresRabbit (RabbitMQ Server)

Tests with the respective attributes will only run when these services exist in the test/development environment.  Connection strings/ports for these services can be found in:

- TestDatabases.txt (Relational Databases)
- default.yaml (RabbitMQ / MongoDb)
- Mongo.yaml
- Rabbit.yaml
- RelationalDatabases.yaml

For setting up the RDMP platform databases see https://github.com/HicServices/RDMP/blob/master/Documentation/CodeTutorials/Tests.md

### 2.3 Monitoring and Alerting Plan

Both SmiServices and [RDMP] use [NLog] for logging.  A number of additional systems are incorporated into SmiServices for logging/audit:

- The health of running services can be monitored using the RabbitMQ admin console
- Each message sent by a service has a GUID associated with it.  Records loaded into the relational database include the GUIDs of all services that acted on the tag data during ETL (see [Logging through the IMessageHeader](./src/common/Smi.Common/README.md#logging-through-the-imessageheader)).
- Anonymisation effectiveness is verified using [IsIdentifiable] and can be visualised using the [IsIdentifiableReviewer]
- All ETL activities are audited in the hierarchical RDMP logging database.  This includes run duration, record counts, error messages during load etc.

### 2.4 Release / Roll-out and Deployment Plan

The latest binaries can be downloaded from the [GitHub releases page](https://github.com/SMI/SmiServices/releases/latest).  

Each release has a [Changelog] describing all changes made to the codebase (additions, bugfixes etc).  Changelog entries include a link to the 'git diff' which shows code changes and more technical descriptions of changes.

The data load microservices wrap the [RDMP] DLE.  It is therefore important to ensure that the RDMP gui/CLI client and platform databases are maintained at a compatible version with the binary shipping with SmiServices.  RDMP updates are backwards compatible and natively support running old clients against new versions of the platform database.  Dependencies are managed with [dependabot](https://github.blog/2020-06-01-keep-all-your-packages-up-to-date-with-dependabot/) which automatically generates a 'pull request' when updates are available to [RDMP Dicom].  If there are any compatibility issues this will surface in the automated integration testing of the pull request.

The easiest way to consume SmiServices is through the [Docker Image](https://github.com/jas88/smideploy).

For a more adaptable/scalable setup or to use existing infrastructure (databases etc), you will need to:

 - Install RabbitMq
 - Install MongoDb
 - Install a relational database (MySql, SqlServer, Postgres or Oracle)
 - [Install RDMP](https://github.com/HicServices/RDMP) and setup platform databases on a Sql Server.

After all services are in place:
 - Import the [queue definitions](./data/rabbitmqConfigs) into RabbitMQ via the admin web interface (or from command line).  Make any changes for vhost etc if desired
 - Download the software assets for the latest [SmiServices Release](https://github.com/SMI/SmiServices/releases) and unzip the appropriate service e.g. `smi-services-v4.0.0-linux-x64.tgz`

Configure 'default.yaml'
 - Update the credentials to match your RabbitMQ and MongoDb instance (you can ignore RDMP, Redis etc for now).
 - Remove the `TEST.` prefix on queue names

### 2.5 Rollback Plan

SmiServices releases are backwards compatible and released as a self contained executable package.  Rolling back to an earlier version of the software involves only deleting the new binary and restoring the old one.  

Breaking changes to SmiServices can be expected only when the major version number is incremented.  This is consistent with [semantic versioning](https://semver.org/).  In such cases the changelog should describe the changes and how to rollback.

### 2.6 Associated Documentation

Each SmiServices dependency contains its own documentation.  Follow the link from the package name for help.

| Package | Status |
| :---        |    :----:   |
| [HIC.TypeGuesser](https://github.com/HicServices/TypeGuesser) |   [![Build, test and package](https://github.com/HicServices/TypeGuesser/actions/workflows/dotnet.yml/badge.svg)](https://github.com/HicServices/TypeGuesser/actions/workflows/dotnet.yml)  [![Coverage Status](https://coveralls.io/repos/github/HicServices/TypeGuesser/badge.svg?branch=master)](https://coveralls.io/github/HicServices/TypeGuesser?branch=master)  [![Total alerts](https://img.shields.io/lgtm/alerts/g/HicServices/TypeGuesser.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/HicServices/TypeGuesser/alerts/)  [![NuGet Badge](https://buildstats.info/nuget/HIC.TypeGuesser)](https://buildstats.info/nuget/HIC.TypeGuesser)
| [FAnsiSql](https://github.com/HicServices/FAnsiSql) | [![.NET Core](https://github.com/HicServices/FAnsiSql/actions/workflows/dotnet-core.yml/badge.svg)](https://github.com/HicServices/FAnsiSql/actions/workflows/dotnet-core.yml)  [![Total alerts](https://img.shields.io/lgtm/alerts/g/HicServices/FAnsiSql.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/HicServices/FAnsiSql/alerts/) [![NuGet Badge](https://buildstats.info/nuget/HIC.FAnsiSql)](https://www.nuget.org/packages/HIC.FansiSql/) |
| [DicomTypeTranslation](https://github.com/HicServices/DicomTypeTranslation) |  [![.NET Core](https://github.com/HicServices/DicomTypeTranslation/actions/workflows/dotnet-core.yml/badge.svg)](https://github.com/HicServices/DicomTypeTranslation/actions/workflows/dotnet-core.yml) [![Total alerts](https://img.shields.io/lgtm/alerts/g/HicServices/DicomTypeTranslation.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/HicServices/DicomTypeTranslation/alerts/) [![NuGet Badge](https://buildstats.info/nuget/HIC.DicomTypeTranslation)](https://buildstats.info/nuget/HIC.DicomTypeTranslation) |
| [RDMP](https://github.com/HicServices/RDMP) |   [![Build status](https://github.com/HicServices/RDMP/workflows/Build/badge.svg)](https://github.com/HicServices/RDMP/actions?query=workflow%3ABuild) [![Coverage Status](https://coveralls.io/repos/github/HicServices/RDMP/badge.svg?branch=develop)](https://coveralls.io/github/HicServices/RDMP?branch=develop) [![Total alerts](https://img.shields.io/lgtm/alerts/g/HicServices/RDMP.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/HicServices/RDMP/alerts/) [![NuGet Badge](https://buildstats.info/nuget/HIC.RDMP.Plugin)](https://buildstats.info/nuget/HIC.RDMP.Plugin) |
| [RDMP.Dicom](https://github.com/HicServices/RdmpDicom) | [![Build status](https://github.com/HicServices/RdmpDicom/actions/workflows/dotnet-core.yml/badge.svg)](https://github.com/HicServices/RdmpDicom/actions/workflows/dotnet-core.yml) [![Total alerts](https://img.shields.io/lgtm/alerts/g/HicServices/RdmpDicom.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/HicServices/RdmpDicom/alerts/) [![NuGet Badge](https://buildstats.info/nuget/HIC.RDMP.Dicom)](https://buildstats.info/nuget/HIC.RDMP.Dicom) |
| [Dicom Template Builder] | [![Build Status](https://travis-ci.org/HicServices/DicomTemplateBuilder.svg?branch=master)](https://travis-ci.org/HicServices/DicomTemplateBuilder)  [![Total alerts](https://img.shields.io/lgtm/alerts/g/HicServices/DicomTemplateBuilder.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/HicServices/DicomTemplateBuilder/alerts/)  |

## 3.0 Further Considerations

### 3.1 Human Resource Requirements

Setting up an SmiServices deployment on a large volume of imaging data requires a number of staff with skills in the following areas:

- Parallel file system setup and optimisation
- Database administration (MongoDb and chosen Relational Database)
- Deployment and configuration skills ([YAML](https://yaml.org/), )
- Infrastructure management (setting up VMs, server instances etc)

Once running SmiServices is designed to be used by staff with a background in cohort building and data analysis.  Familiarity with SQL is beneficial but not required by all analysts if the SQL is abstracted into library methods or gui applications (e.g. [RDMP] Cohort Builder)

The task of loading identifiable data from disk and understanding it in the MongoDb requires:

- Knowledge of the MongoDb query language
- Appropriate governance training for working with identifiable data (e.g. NHS contract)

### 3.2 Security considerations

SmiServices is designed to run in a secure offline environment in which network firewalls prevent interaction between any VMs and Servers with the outside world.  There are no web frontends or management consoles designed for external access.  All user interaction with tools in the suite should be done over secure VPN via SSH or secure Remote Desktop (for graphical applications).

SmiServices relies on the security of the environment into which it is deployed. Database access permissions and file permissions should be configured appropriately to the users operating the system.  Where possible connection strings and credentials should use 'integrated security' (i.e. the user account credentials not username/password).  SmiServices also supports supplying passwords in environment variables for when that is an acceptable solution.  If supplying passwords in the configuration file directly then this file should be protected from all non root users (i.e. with appropriate file access permissions).

### 3.3 Legal and Ethical considerations

SmiServices users must have appropriate governance permissions for accessing data.

Users handling the ETL and anonymisation processes must have approval to view identifiable information.

Cohort Building takes place in the anonymised relational database and therefore may require a lower level of governance approval.

### 3.4 Accessibility considerations

SmiServices command line tools are written in dotnet and work with any windows or linux terminals (bash, powershell etc).  Accessibility of terminals is typically good with inbuilt or compatible addons for controlling text size, converting text to speech etc.

All commands run in the RDMP windows gui client can [also be run directly on the command line](https://github.com/HicServices/RDMP/blob/develop/Documentation/CodeTutorials/RdmpCommandLine.md) which allows for scripting and alternate access routes.  RDMP also features a Terminal User Interface (TUI) which mirrors the design of the gui client. 

![rdmp terminal user interface](./docs/Images/rdmp-tui.png)

Both [IsIdentifiableReviewer] and the RDMP TUI support specifying alternative colour schemes for text where contrast or colors used are an accessibility issue for users.

### 3.5 Support considerations

Support for users of SmiServices is via [GitHub Issues](https://github.com/SMI/SmiServices/issues).  Tickets can also be raised in [RDMP] if the issue does not seem to be related specifically to the microservices component of the suite.

## 4.0 Future Work

### 4.1 Residual work estimates and timelines

SmiServices is feature complete for the loading, linking, cohort building and extraction of anonymous images.  Future work will focus on improvements such as:

- Better pixel anonymisation 
- Support for natural language processing (NLP) of structured reports in cohort building
- Project specific UID anonymisation

Not all of the above features will be directly implemented in the SmiServices repository and may be consumed as external pluggable resources (e.g. through RDMP cohort builder plugins).

### 4.2 Open Questions

Outstanding questions and feature debates can be seen in the [GitHub Issues](https://github.com/SMI/SmiServices/issues) and [GitHub Discussions](https://github.com/SMI/SmiServices/discussions) pages of SmiServices repository.

## 5.0 End Matter

### 5.1 Related Work 

Many image processing tools have relational database backends and most support export in interchange formats such as CSV / XLS or querying APIs.  This means that tools such as XNAT, [ORTHANC](https://www.orthanc-server.com/) and PACS servers can already be part of a cohort building process.  However this can be slow and cumbersome especially for large linkage operations involving multiple other datasets.  For example answering questions such as:

> How many patients had a head CT within 6 months of being prescribed Drug X and then were subsequently hospitalised again within 3 months of the CT StudyDate

Answering such a question using DICOM QR or querying a backend schema would be complex and slow for large cohorts.  It would also likely involve multiple steps e.g. running C-FIND queries to retrieve study dates and loading to the EHR relational database.

SmiServices improves upon this approach by loading an imaging archive to a single anonymous indexed table per modality with aggregate
tables containing Study and Series level information.  This format is the fastest queryable representation for large scale linkage and extraction.

### 5.2 Acknowledgments

Contributors to the SmiServices repository can be seen in the [GitHub Contributors page](https://github.com/SMI/SmiServices/graphs/contributors).

Contributors to the [upstream dependencies](./PACKAGES.md) can be seen in the corresponding section of those repositorys.

This work has been funded under the PICTURES programme.  A 5-year, £3.8M programme of work funded by the Medical Research Council (MRC) with additional support from the Engineering and Physical Sciences Research Council (EPSRC) as part of Health Data Research UK (HDR UK).

The authors acknowledge the support from the Farr Institute of Health Informatics Research and Dundee University Medical School. This work was supported by the Medical Research Council (MRC) grant No. MR/M501633/1 (PI: Andrew Morris) and the Wellcome Trust grant No. WT086113 through the Scottish Health Informatics Programme (SHIP) (PI: Andrew Morris). SHIP is a collaboration between the Universities of Aberdeen, Dundee, Edinburgh, Glasgow, and St Andrews, and the Information Services Division of NHS Scotland. This project has also been supported by MRC and EPSRC (grant No. MR/S010351/1) and by the Chief Scientist Office of the Scottish Government Health and Social Care Directorates via a leverage grant made to Farr Scotland. The project was also supported by the Scottish Government through the “Imaging AI” grant award.

This work was supported by Health Data Research UK, which receives its funding from HDR UK Ltd (HDR-5012) funded by the UK MRC, Engineering and Physical Sciences Research Council, Economic and Social Research Council, Department of Health and Social Care (England), Chief Scientist Office of the Scottish Government Health and Social Care Directorates, Health and Social Care Research and Development Division (Welsh Government), Public Health Agency (Northern Ireland), British Heart Foundation (BHF), and the Wellcome Trust

## Microservices

All microservices [follow the same design pattern](./src/common/Smi.Common/README.md).

The following  microservices have been written.  Microservices are loosely coupled, usually reading and writing only a single kind of message.  Each Queue and Exchange as implemented supports only one Type of `Smi.Common.Messages.IMessage`.

Microservices can be configured through [the configuration file](./data/microserviceConfigs/default.yaml).

A control queue is provided for controlling Microservices during runtime.  It supports a [limited number of commands](./docs/control-queues.md).

### Data Load Microservices

![loaddiagram](./docs/Images/LoadMicroservices.png)

| Microservice / Console App| Description |
| ------------- | ------------- |
| [ProcessDirectory]  | Command line application that finds dicom files on disk and [queues them for execution in RabbitMQ](./src/common/Smi.Common/Messages/AccessionDirectoryMessage.cs).|
| [DicomTagReader] | Opens queued dicom files on disk and [converts them to JSON](./src/common/Smi.Common/Messages/DicomFileMessage.cs).  Also creates a [summary record of the whole series](./src/common/Smi.Common/Messages/SeriesMessage.cs).|
| [IdentifierMapper] (Optional)  | Replaces the `PatientID` dicom Tag in a [DicomFileMessage] using a specified mapping table.|
| [MongoDBPopulator]  | Persists the dicom Tag data in [DicomFileMessage] and/or [SeriesMessage] into a MongoDB database document store. |
| [DicomRelationalMapper] | Persists the dicom Tag data (and file paths) in [DicomFileMessage] into a [relational database](https://github.com/HicServices/RDMP/blob/develop/Documentation/CodeTutorials/FAQ.md#databases).  ETL pipeline is controlled by an [RDMP] data load configuration.|
| [DicomReprocessor] | Runs a MongoDB query on the database populated by [MongoDBPopulator] and converts the results back into [DicomFileMessage] for (re)loading by [DicomRelationalMapper].|

### Image Extraction Microservices

![extractiondiagram](./docs/Images/ExtractionMicroservices.png)

| Microservice / Console App| Description |
| ------------- | ------------- |
| [IsIdentifiable]  | Evaluates data being prepared for extraction for personally identifiable data (PII).  See also [IsIdentifiableReviewer]|
| [ExtractImages] | Reads UIDs from a CSV file and generates [ExtractionRequestMessage] and audit message [ExtractionRequestInfoMessage].|
| [CohortExtractor] | Looks up SeriesInstanceUIDs in [ExtractionRequestMessage] and does relational database lookup(s) to resolve into physical image file location.  Generates  [ExtractFileMessage] and audit message [ExtractFileCollectionInfoMessage].|
| [CTPAnonymiser]  | Microservice wrapper for [CTP](https://github.com/johnperry/CTP).  Anonymises images specified in  [ExtractFileMessage] and copies to specified output directory.  Generates audit message [ExtractedFileStatusMessage].|
| [CohortPackager] | Records all audit messages and determines when jobs are complete.|

### Audit and Logging Systems

| Audit System | Description|
| ------------- | ------------- |
| [NLog](http://nlog-project.org/) | All Microservices log all activity to NLog, the manifestation of these logs can be to file/console/server etc as configured in the app.config file.|
| [Message Audit](./src/common/Smi.Common/README.md#logging) | Every message sent by a microservice has a unique Guid associated with it.  When a message is issued in response to an input message (all but the first message in a chain) the list of legacy message Guids is maintained.  This list is output as part of NLog logging.|
| [Data Load Audit](./src/microservices/Microservices.DicomRelationalMapper/Readme.md#7-audit)|The final Message Guid of every file identified for loading is recorded in the relational database image table.  In addition a valid from / data load ID field is recorded and any UPDATEs that take place (e.g. due to reprocessing a file) results in a persistence record being created in a shadow archive table.|
| [Extraction Audit (MongoDB)](./src/microservices/Microservices.CohortPackager/README.md) | CohortPackager is responsible for auditing extraction Info messages from all extraction services, recording which images have been requested and when image anonymisation has been completed.  This is currently implemented through `IExtractJobStore`.|
| CohortExtractor Audit | Obsolete interface `IAuditExtractions` previously existed to record the linkage results and patient release identifiers.|
| Fatal Error Logging | All Microservices that crash or log a fatal error are shut down and log a message to the Fatal Error Logging Exchange.  TODO: Nobody listens to this currently.|
| Quarantine | TODO: Doesn't exist yet.|

## .sln Overview

Apart from the Microservices (documented above) the following library classes are also included in the solution:

| Project Name | Path | Description|
| ------------- | ----- | ------------- |
| Dicom File Tester |/Applications| Application for testing DICOM files compatibility with Dicom<->JSON and Dicom to various database type conversions and back. It basically takes a file and pushes it through the various converters to see what breaks |
| Dicom Repopulator |/Applications| [See Microservices](#image-extraction-microservices) |
| [Dicom Template Builder] | /Applications| GUI tool for building modality database schema templates.  Supports viewing and exploring dicom tags in files|
| Smi.MongoDB.Common | /Reusable | Library containing methods for interacting with MongoDb |


## Building

### Building the C# Projects

Building requires a [.NET Core SDK](https://dotnet.microsoft.com/download/dotnet-core), at a compatible version as specified in the [`global.json`](/global.json) file.

To build the entire solution from the project root, run:

```bash
$ dotnet build [-r RID]
```

_The RID argument is optional. Use this if you want to build for a different platform e.g. `-r linux-x64` to build for Linux from a Windows machine. See [here](https://docs.microsoft.com/en-us/dotnet/core/rid-catalog) for more info on runtime identifiers._

To build an individual sub-project:

```bash
$ cd src/microservices/Microservices.DicomTagReader/
$ dotnet build
```

This will automatically rebuild any dependent projects which have changes as well.

### Building the Java Projects

Building the Java projects requires Java JDK `>= 1.7` (OpenJDK recommended 🙂), and Maven.

The CTP dependency first needs to be manually installed:

- Linux

```bash
$ cd lib/java/
$ ./installDat.sh
```

- Windows

```bash
$ cd lib\java\
$ .\installDat.bat
```

The projects can then be built and tested by returning to the top level directory and running:

```bash
$ mvn -f src/common/com.smi.microservices.parent/pom.xml clean test
```

This will compile and run the tests for the projects. The full test suite requires a local RabbitMQ server, however these can be skipped by passing `-PunitTests`. The entire test suite can be skipped by instead running `compile`, or by passing `-DskipTests`.

To build a single project and its dependencies, you can do:

```bash
$ mvn -f src/common/com.smi.microservices.parent/pom.xml test -pl com.smi.microservices:ctpanonymiser -am
```

Note: If you have Maven `>=3.6.1` then you can pass `-ntp` to each of the above commands in order to hide the large volume of messages related to the downloading of dependencies.

## Running

All applications and services are runnable through the `smi` program. This is available either in the binary distribution, or in the `src/Applications/Applications.SmiRunner/bin` directory if developing locally. See the SmiRunner [README](/src/applications/Applications.SmiRunner/README.md) for more information.

## Developing

### C# Projects

Development requires Visual Studio 2017 or later. Simply open the SmiServices.sln file.

To run the tests for IsIdentifiable, the Stanford NER classifier is required. This can be downloaded with the included script:

```bash
$ cd data/stanford-ner
$ ./download.sh
```

### Java Projects

Development requires Java JDK `>= 1.7`, and Maven.

### pre-commit

This repo uses [pre-commit] to manage and automatically run a series of linters
and code formatters. After cloning the repo and changing into the directory, run
this once to setup pre-commit.

```console
$ pip install pre-commit
$ pre-commit install
```

This will then run the checks before every commit. It can also be run manually
at any time:

```console
$ pre-commit run [<hook>] (--all-files | --files <file list>)
```

Running pre-commit locally is optional, since it is also run during any PR. To remove
pre-commit from your repo clone, simply run:

```console
$ pre-commit uninstall
```

## Note On Versioning

The C# projects share the same release version, which is controlled by the [SharedAssemblyInfo.cs](src/SharedAssemblyInfo.cs) file. The Java projects are versioned independently, set in their pom files, however in practice they follow the release version of the repo overall.

## Scalability

The services in this repository have been successfully used to load all medical imaging data captured in Scotland's National PACS archive.

Scalability is handled through parallel process execution (using [RabbitMQ]).  This allows slow processes (e.g. reading dicom tags from files on disk) to have more running instances while faster processes have less.  Scalability of large operations (e.g. linkage / cohort identification) is done within the [DBMS] layer.

[RabbitMQ]: https://www.rabbitmq.com/
[DBMS]: https://github.com/HicServices/RDMP/blob/develop/Documentation/CodeTutorials/Glossary.md#DBMS
[Dicom]: ./Glossary.md#dicom
[Dicom tags]: ./Glossary.md#dicom-tags
[IsIdentifiable]: ./src/microservices/Microservices.IsIdentifiable/README.md
[IsIdentifiableReviewer]: ./src/applications/Applications.IsIdentifiableReviewer/README.md
[DicomFileMessage]: ./src/common/Smi.Common/Messages/DicomFileMessage.cs
[SeriesMessage]: ./src/common/Smi.Common/Messages/SeriesMessage.cs
[ExtractionRequestMessage]: ./src/common/Smi.Common/Messages/Extraction/ExtractionRequestMessage.cs
[ExtractionRequestInfoMessage]: ./src/common/Smi.Common/Messages/Extraction/ExtractionRequestInfoMessage.cs
[ExtractFileMessage]: ./src/common/Smi.Common/Messages/Extraction/ExtractFileMessage.cs
[ExtractFileCollectionInfoMessage]: ./src/common/Smi.Common/Messages/Extraction/ExtractFileCollectionInfoMessage.cs
[ExtractedFileStatusMessage]: ./src/common/Smi.Common/Messages/Extraction/ExtractedFileStatusMessage.cs
[RDMP]: https://github.com/HicServices/RDMP
[ProcessDirectory]: ./src/applications/Applications.DicomDirectoryProcessor/README.md
[DicomTagReader]: ./src/microservices/Microservices.DicomTagReader/README.md
[IdentifierMapper]: ./src/microservices/Microservices.IdentifierMapper/Readme.md
[MongoDBPopulator]: ./src/microservices/Microservices.MongoDbPopulator/Readme.md
[DicomRelationalMapper]: ./src/microservices/Microservices.DicomRelationalMapper/Readme.md
[DicomReprocessor]: ./src/microservices/Microservices.DicomReprocessor/README.md
[ExtractImages]: ./src/applications/Applications.ExtractImages/README.md
[CohortExtractor]: ./src/microservices/Microservices.CohortExtractor/README.md
[CTPAnonymiser]: ./src/microservices/com.smi.microservices.ctpanonymiser/README.md
[CohortPackager]: ./src/microservices/Microservices.CohortPackager/README.md
[pre-commit]: https://pre-commit.com
[ExtractImages]: ./src/applications/Applications.ExtractImages/README.md
[DICOM specification]: https://www.dicomstandard.org/
[RDMP Dicom]: https://github.com/HicServices/RdmpDicom
[Dicom Template Builder]: https://github.com/HicServices/DicomTemplateBuilder
[BadMedicine.Dicom]: https://github.com/HicServices/BadMedicine.Dicom
[BadMedicine]: https://github.com/HicServices/BadMedicine
[NLog]: https://nlog-project.org/
[Changelog]: ./CHANGELOG.md
[XNAT]: https://www.xnat.org/
# SMI Services Contributing Guidelines

- TODO 😊
# Glossary

## Dicom

DICOM (Digital Imaging and Communications in Medicine) is the international standard format for storing medical imaging.  In the standard DICOM files contain both pixel data (the image) and [tag data].

## Dicom Tags

Images stored in the [Dicom] format include a metadata information (tags).  Tags include information about the patient as well as information about the device used, scan settings, modality etc.  Some tags are required for the file to be considered valid while others are optional.  The format includes support for both tree data structures (Sequences) and arrays (Multiplicity).

Tags can be either part of the [standard set](https://dicom.innolitics.com/ciods) or invented by device manufacturers (private tags).  Adherence to the standard varies widely (both by manufacturer, healthboard and over time).

## Loading

[Dicom] images can be very large but most of this space is taken up by the pixel data.  Since it is the [tag data] that is most useful for generating cohorts it is these that are loaded.

'Loading' in the context of the SMI repository means extracting the dicom [tag data] and persisting it along with the relative file path of the image into a database.  

Loading may also include typical data load operations e.g. summarization, error detect, anonymisation etc.

[Dicom]: #dicom
[Tag Data]: #dicom-tags## Tests

SMI Microservices use RabbitMQ, MongoDb, Sql Server and MySql.  These tests are run automatically by [Travis CI](https://travis-ci.org/SMI/SmiServices) (see [travis.yml](.travis.yml)).

## Connection Strings

Once you have set up the above dependencies you will need to set the connection strings in the following files (which should not be committed to the repository):

- [TestDatabases.txt](./tests/common/Smi.Common.Tests/TestDatabases.txt)
- [RelationalDatabases.yaml](./tests/common/Smi.Common.Tests/RelationalDatabases.yaml)
- [Mongo.yaml](./tests/common/Smi.Common.Tests/Mongo.yaml)
- [Rabbit.yaml](./tests/common/Smi.Common.Tests/Rabbit.yaml)

Tests involving [RDMP](https://github.com/HicServices/RDMP/) require the [RDMP databases to be set up](https://github.com/HicServices/RDMP/blob/develop/Documentation/CodeTutorials/Tests.md#database-tests) on the Sql Server.  For example:

If running on linux / travis the following command will do that (use latest version in url):

```
$ wget https://github.com/HicServices/RDMP/releases/download/v3.2.1/rdmp-cli-linux-x64.zip
$ unzip -d rdmp-cli rdmp-cli-linux-x64.zip || true # Ignore exit code since unzip returns 1 for a warning we don't care about
$ cd rdmp-cli
$ chmod +x rdmp
$ ./rdmp install localhost TEST_ -u sa -p 'YourStrongPassw0rd'
```
Adds the basis for a new "DicomAnonymiser" microservice, which can be used in place of the existing Java CTP service. It supports pluggable anonymisers through the `IDicomAnonymiser` interface. No implementations are provided in this PR.
Treat all build warnings as errors, and fix or disable existing ones. Also remove unused System.Security.AccessControl package.
Allow newer minor SDK versions to build the sln. global.json specifies the minimum version which will be used in the CI.
Fix #891; jobId not created until extraction is completed.
Fix deprecated python collections import for py310+
Update CI script to pull tessdata file from the new branch name (main not master)
Ensure control exchange exists for ExtractImagesHostTests
# News files

This directory contains files describing changes since the previous release of SmiServices.

When a release is built, these files are automatically combined into the main [CHANGELOG](/CHANGELOG.md), and are then deleted.

## File naming

News file names should be of the form

```txt
<pr#>-<type>.md
```

e.g.

```txt
1234-feature.md
```

Where `type` is one of

-   `feature`
-   `change`
-   `bugfix`
-   `doc`
-   `removal`
-   `meta`

*Note* Ensure that the file is named with the _PR_ number, rather than any associated _issue_ number.

Quick tip: You can get the most recent issue or PR number with the following one-liner. Then add one to determine the new one for your PR (so long as you're quick!)

```console
$ curl -s "https://api.github.com/repos/smi/smiservices/issues?sort=created&direction=desc&per_page=1&page=1" | jq .[].number
702
```

The file should contain a short description of the patch as one or more lines of markdown, either as a top-level sentence

```md
Fixed a foobar
```

or, if more detail is required, multiple lines formatted as a sub-list

```md
Fixed a foobar
-   Requires users to change xyz
```
Fix #913; error if file extraction command is re-run after being cancelled
Fixes the current CI issues by restricting the creation of the ControlExchange to the RequiresRabbit test decorator.
Remove reference to System.Drawing.Common
Fix issue #921 - erroneous stripping of root path to empty string in GlobalOptions
Ensure the python version used in CI runs is exactly what we specify
Always convert PHI to XML regardless of annotation_mode setting.
Add missing reference to NLog.
# SMI Python library

Some useful functions in Python

## Requirements

```
deepmerge
pika
pydicom
pymongo
PyYAML
xml.etree (comes with python)
```

## Installation

Run `python3 setup.py bdist_wheel` to create `Smi_Services_Python-0.0.0-py3-none-any.whl`

Run `python3 setup.py install` to install (including dependencies) into your python site-packages
(whether that be global or inside a current virtualenv).

Note that the version number is read from AssemblyInfo.cs in a parent directory.

## Testing

Test all modules:

```
pytest SmiServices/*.py
```

Test each module individually, for example:
```
python3 -m pytest SmiService/Dicom.py
python3 -m pytest SmiService/DicomText.py
python3 -m pytest SmiService/StructuredReport.py
```

## Usage

For example:

```
if 'SMI_ROOT' in os.environ:     # $SMI_ROOT/lib/python3
    sys.path.append(os.path.join(os.environ['SMI_ROOT'], 'lib', 'python3'))
from SmiServices import Mongo
from SmiServices import Rabbit
from SmiServices import Dicom
from SmiServices import DicomText
from SmiServices import StructuredReport as SR
```

## Dicom.py

Mostly low-level functions for reading DICOM files that are used by the DicomText module.

## DicomText.py

Provides a DicomText class which assists in parsing a DICOM Structured Report.
Also has functions for redacting the text given a set of annotations.
Uses the pydicom library internally.

Typical usage:

```
dicomtext = Dicom.DicomText(dcmname) # Reads the raw DICOM file
dicomtext.parse()                    # Analyses the text inside the ContentSequence
xmldictlist = Knowtator.annotation_xml_to_dict(xml.etree.ElementTree.parse(xmlfilename).getroot())
dicomtext.redact(xmldictlist)        # Redacts the parsed text using the annotations
dicomtext.write(redacted_dcmname)    # Writes out the redacted DICOM file
OR
write_redacted_text_into_dicom_file  # to rewrite a second file with redacted text
```

## Knowtator.py

Provides a function for parsing the XML files containing annotations
as output by the SemEHR anonymiser and input to eHOST. The files are typically
named `.knowtator.xml` and have the format:

```
 <annotation>
  <mention id="anon.txt-1"/>
  <annotator id="semehr">semehr</annotator>
  <span end="44" start="34"/>
  <spannedText>16 year old</spannedText>
  <creationDate>Wed November 11 13:04:51 2020</creationDate>
 </annotation>
 <classMention id="anon.txt-1">
  <mentionClass id="semehr_sensitive_info">16 year old</mentionClass>
 </classMention>
```

The function `annotation_xml_to_dict` parses the XML and returns a suitable Python dict.

Also contains a function to write such XML files, useful when testing, or when converting from a Phi file.

## Mongo.py

Very simple wrapper around pymongo specifically for SMI.
Provides a `SmiPyMongoCollection` class. Typical usage:

```
mongodb = Mongo.SmiPyMongoCollection(mongo_host)
mongodb.setImageCollection('SR')
mongojson = mongodb.DicomFilePathToJSON('/path/file2')
print('MONGO: %s' % mongojson)
```

## Rabbit.py

Python interface to the SMI RabbitMQ messaging system.
Provides a class `smiMessage` which is inherited by task-specific classes
`CTP_Start_Message` and `IsIdentifiable_Start_Message`.
Provides classes `RabbitProducer` and `RabbitConsumer`.
Has sample functions `send_CTP_Start_Message` and `get_CTP_Output_Message` for testing.

One known problem with using RabbitMQ from Python is the lack of data types,
in particular no concept of a difference between 16-bit and 32-bit integers.
The `pika` library tries to be efficient by constructing a message using a 16-bit
integer if its value will fit, but the C# and Java program have been written to
explicitly expect a 32-bit integer, and they crash if given a 16-bit one.
The pika library currently has no way to request a 32-bit integer as the data
type is determined dynamically so we have to omit the Timestamp field from the
messages.

## StructuredReport.py

Provides a function `SR_parse` which can parse a Python dict containing a DICOM
Structured Report and return the content as a usable string.  The dict can be
read from a DICOM file using pydicom or can be obtained from the MongoDB database
which represents the data in a similar but different format (i.e. no VR tag).

Some utility functions are used by the DicomText module.


# Parsing Structured Reports

```
# Read the DICOM file manually and extract to JSON
  dicom_raw = pydicom.dcmread(filename)
  dicom_raw_json = dicom_raw.to_json_dict(None, 0)
# OR read the JSON from MongoDB
  mongodb = Mongo.SmiPyMongoCollection(mongo_host)
  mongodb.setImageCollection('SR')
  mongojson = mongodb.DicomFilePathToJSON(args.input)
```


## Method 1 - use the StructuredReport module

```
  SR.SR_parse(dicom_raw_json, document_name, output_fd)
  SR.SR_parse(mongojson, document_name, output_fd)    
```

## Method 2 - use the pydicom walk method

```
def decode(filename):
    dicom_raw = pydicom.dcmread(filename)

    def dataset_callback(dataset, data_element):
    	if data_element.VR == 'SQ':
    		False
    	elif data_element.VR in ['SH', 'CS']:
    		False
    	elif data_element.VR == 'LO':
    		print('[[%s]]' % str(data_element.value))
    	else:
    		print('%s' % (str(data_element.value)))

    # Recurse only the values inside the ContentSequence
    for content_sequence_item in dicom_raw.ContentSequence:
    	content_sequence_item.walk(dataset_callback)
```

## Method 3 - use the pydicom recurse_tree method

```
def decode(filename):

    def recurse_tree(tree, dataset, parent):
    	for data_element in dataset:
    		# the node_id could be used as a unique reference
    		node_id = parent + "." + hex(id(data_element))
    		if isinstance(data_element.value, str):
    			# Useless data types are SH (eg. RE.05 or 99_OFFIS), CS (eg. CONTAINS)
    			if data_element.VR in ['SH', 'CS']:
    				False
    			elif data_element.VR == 'LO':
    				# LO is like a heading
    				print('[[%s]]' % (str(data_element.value)))
    			else:
    				# UT is text, DA is date, PN is name
    				print('%s' % (str(data_element.value)))
    		else:
    			# Non-string values are useless, sequences are handled below anyway
    			False #print('%s val = %s' % (node_id, str(data_element.value)))
    		if data_element.VR == "SQ":  # a sequence
    			for i, dataset in enumerate(data_element.value):
    				item_id = node_id + "." + str(i + 1)
    				sq_item_description = data_element.name.replace(" Sequence", "")  # XXX not i18n
    				item_text = "{0:s} {1:d}".format(sq_item_description, i + 1)
    				#print('%s seq = %s' % (item_id, item_text))
    				recurse_tree(tree, dataset, item_id)

    dicom_raw = pydicom.dcmread(filename)
    dicom_raw.decode()   # XXX should we decode to UTF?

    # Recurse only the values inside the ContentSequence
    # (to recurse the whole DICOM pass dicom_raw as second param).
    for content_sequence_item in dicom_raw.ContentSequence:
    	recurse_tree(None, content_sequence_item, '')
```

## Method 4 - use the DicomText module

To read an original:

```
    dicomtext = DicomText.DicomText(input)
    dicomtext.parse()
```

To redact, given an input_xml file:

```
    xmlroot = xml.etree.ElementTree.parse(args.input_xml).getroot()
    xmldictlist = Knowtator.annotation_xml_to_dict(xmlroot)
    dicomtext.redact(xmldictlist)
    dicomtext.write_redacted_text_into_dicom_file(args.output_dcm)
```

# First Install

For the moment you have to manually install DAT and its dependencies to your local Maven repo for the build to find it. To do this, run installDat.bat (installDat.sh on linux) once before you try and build for the first time.

# Maven commands

You can do all this from inside Eclipse but will have to set up the Run Configurations yourself. Otherwise, from the SMIPlugin/java/Microservices directory, the commands are:

**To build all projects**:
`mvn clean package`

Can append `-DskipTests` to that if you want to build without running tests, or change the `skipTests` option in `SMIPlugin/Java/Microservices/pom.xml` file to `true`.

**To build a single project**

`mvn clean package -pl Microservices.ExtractorCL -am`

To build only the ExtractorCL program as well as the required `Microservices.Common` project, for example. Maven flags used are:

`-pl,--projects <arg> // Comma-delimited list of specified reactor projects to build instead of all projects. A project can be specified by [groupId]:artifactId or by its relative path`

`-am, --also-make // If project list is specified, also build projects required by the list`

**To run**:

First `cd` into one of the `Microservice.*` directories, then:

`java -jar target\***-portable-***.jar ...`

where `...` are any command line options you wish to supply the microservice. 

**To package deployable archives**:
Back in `SMIPlugin/Java/Microservices`, build the projects as above, then run:

`mvn assembly:single@create-deployable`

this will create zip archives in each of the target directories containing all the resources and jars needed to deploy the microservices on another machine.

# Configs

The Java microservices share the same yaml config files as the C# microservices. Documentation can be found [here](../Microservices/Microservices.Common/Options/RabbitMqConfigOptions.md).

# Logging notes

We use Logback in Java to log to both file and console. Logback seems to be the best successor to log4j, see this excellent resource to learn more: https://stackify.com/logging-logback/

To properly set up the logging across the microservices, we use the following method:

1. The file `SmiLogbackConfig.xml` is copied from `SMIPlugin\Java\Microservices\res` to each of the output target directories.
2. At the beginning of the main method of each of the microservices, we call `SmiLogging.Setup();`, which is a helper class to configure everything.
3. In every class of the program, you can then create a logger object with `private static final Logger _logger = LoggerFactory.getLogger(MyClass.class);`.
4. Then just implement your log messages appropriately.

Everything is logged to file, `INFO`, `WARN`, and `ERROR` are also logged to console.

# Microservice Hosts

## Contents
- [Implementing a Host](#implementing-a-host)
- [Implementing a Consumer](#implementing-a-consumer)
- [Logging](#logging)
- [Rules of Microservice Club](#rules-of-microservice-club)
  - [The First Rule](#the-first-rule)
  - [The Second Rule](#the-second-rule)
- [Class Diagram](#class-diagram)

## Implementing a Host
First load an instance of `GlobalOptions` in your `Program.cs`. 
This is the base class for specifying all your configuration options.

```csharp
public class Program
{
    public static int Main(string[] args)
    {
        var options = new GlobalOptionsFactory().Load();
        
        // ...
    }
}
```

Next you need to decide which queues you want to read from and which you want to write to. 
For this example let's assume you want to consume only from a single queue. 
Your Consumer Options will be available from the GlobaOptions. 
The `FatalLoggingProducerOptions` instance is created by the base host class.

```csharp
public class Program
{
    public static int Main(string[] args)
    {
        var options = new GlobalOptionsFactory().Load(); // will use the 'default.yaml' file
        
        var consumerOptions = options.MyHostOptions; // you don't really need this here...

        // ...
    }
}
```

For this to work you will need to update [default.yaml](../../../data/microserviceConfigs/default.yaml)

```yaml
# ... other stuff above

MyHostOptions: #you can also put this the following into a subclass to avoid cramming too many things at the root level
    QueueName: 'MyQueueName'
    ConsumerTag: 'MyQueueTag'
    QoSPrefetchCount: 1
    AutoAck: false
    # other options you may need

# ... other stuff below
```

If this is a brand new Host, also add the relevant bit into the `GlobalOptions`:

```csharp
public class GlobalOptions
{
    // SNIP LOTS OF CODE
        
    #region AllOptions

    // ... other stuff above
    public MyHostOptions MyHostOptions { get; set; }

    #endregion
}    
    
// new class for the new options
public class MyHostOptions : ConsumerOptions
{
    // other options go here. ConsumerOPtions are inherited.
}
```

Next create a derived class of `MicroserviceHost` this class should take all options needed to do it's job

```csharp
public class MyHost : MicroserviceHost
{
    private readonly ConsumerOptions _consumerOptions;

    public MyHost(GlobalOptions options, bool loadSmiLogConfig = true)
        : base(options, loadSmiLogConfig)
    {
        // Load all the options you need, do all the checks you need.
        _consumerOptions = options.MyHostOptions;
        //consumer = new MyConsumer(consumerOptions); //see next section for how to implement this
    }

    public override void Start()
    {
        RabbitMQAdapter.StartConsumer(_consumerOptions, consumer);
    }
}
```

Finally create an instance of `MicroserviceHostBootstrapper`.
This class will handle construction / connection issues during construction/startup of your `MicroserviceHost`.

```csharp
public class Program
{
    public static int Main(string[] args)
    {
        var options = new GlobalOptionsFactory().Load();
            
        var bootstrapper = new MicroserviceHostBootstrapper(
            () => new CohortPackagerHost(options));
        return bootstrapper.Main();
    }
}
```

### Expected Results
At this stage running the program is meaningful, it should give you sensible logs complaining about missing exchanges on your RabbitMQ Server. 
You can now explore how to change (in the yaml file) / create these yourself.

When you have resolved the exchanges/queues you should get an error relating to `consumer` being null (we commented it out remember). 
Proceed to the next section to see how to implement an `IConsumer`

## Implementing a Consumer

A consumer is a class which listens to a RabbitMQ queue and does something based on the messages that appear. 
Create a new class derived from the `Consumer` abstract class.

```csharp
public class MyConsumer : Consumer
{
    // Add any parameters required
    public MyConsumer()
    {
        // Your constructor setup
    }

    protected override void ProcessMessageImpl(IMessageHeader header, IModel model, BasicDeliverEventArgs basicDeliverEventArgs)
    {
        // Deserialize the message from the delivery arguments

        MyMessage message;
        if (!SafeDeserializeToMessage(header, deliverEventArgs, out message))
            return;
    
        // Do your work here, Ack or Nack depending on result

        if(success)
            Ack(header, deliverEventArgs);
        else
            ErrorAndNack(header, deliverEventArgs, message, exception)
    }
}
```

The `ProcessMessageImpl` method is where you will do your processing and must either `Ack` or `ErrorAndNack`.

The `IMessageHeader` contains provenance information about the message being dequeued.  You can use it for logging (see below) 
and must also supply it when producing new messages (this ensures the message audit chain is kept in tact).

## Logging

All microservices should be derived from the `MicroserviceHost` class, which ensures that a standard logging config is applied. The logging configuration is loaded from the NLog configuration file
specified in the `FileSystemOptions.LogConfigFile` config variable if it exists and is not an empty string, or otherwise from the file
`Smi.NLog.config` in the current directory. This configuration details how file based / console based logging happens and what levels are
processed / ignored etc.

This means that you can get a logger at any time by calling:
```csharp
var logger = NLog.LogManager.GetCurrentClassLogger();
```

The meaning of log levels (Trace, Debug, Info etc) are exactly as defined in the NLog standard:

https://github.com/NLog/NLog/wiki/Configuration-file#log-levels

The `Trace` logging level should be reserved for fine grained timing/performance code only - the lowest level for logging operational
messages should therefore be `Debug`. Trace logging will be disabled unless the CLI option `--trace-logging` is provided.

### Logging through the IMessageHeader
In addition to using the `Log` methods to log routine events, you can log message specific events via `IMessageHeader`:

```csharp
protected override void ProcessMessageImpl(IMessageHeader header, IModel model, BasicDeliverEventArgs basicDeliverEventArgs)
{
    var logger = LogManager.GetCurrentClassLogger();
    header.Log(logger,LogLevel.Warn, "Message was all caps, had to call .Lower on it");
}
```

Logging through a header means that the Guid of the message (and the Guid all previous messages in the chain) will appear in the log e.g.
![Class Diagram](Images/MessageGuidLogs.png)

Logging through the header is recommended whenever the audited fact relates specifically to the content of the message (e.g. couldn't open a file referenced in a `DicomFileMessage`).  Logging through the header automatically happens when sending and acknowledging messages, this results in a view of every message the system sent and the relationship tree of knock on messages (see image above).


## Rules of Microservice Club

### The First Rule

The first rule of Microservice Club is that `LogLevel.Fatal` means game over. Do not log to this level, instead you should call the `Fatal` method:

```csharp
protected override void ProcessMessageImpl(IMessageHeader header, IModel model, BasicDeliverEventArgs basicDeliverEventArgs)
{
    var logger = LogManager.GetCurrentClassLogger();
    header.Log(logger,LogLevel.Warn, "Message was all caps, had to call .Lower on it");

    Fatal("Ran out of memory", new InsufficientMemoryException());

    //ErrorAndNack(header,model,basicDeliverEventArgs,"Something went wrong", new Exception("What went wrong"));
    Ack(header, model, basicDeliverEventArgs);
}
```

The Fatal method is in both `Consumer` and `MicroserviceHost` and causes the current Microservice to shutdown cleanly 
and log a `FatalErrorMessage` to the RabbitMQ fatal message exchange (See `FatalLoggingProducerOptions`).

If your `ProcessMessageImpl` throws an unhandled Exception then a Fatal shutdown will automatically occur.

### The Second Rule
The second rule of Microservice Club is you don't nack messages without giving a reason. This is facilitated through the `Consumer.ErrorAndNack` method.  
This will log an error to NLog and Nack the message for you.

```csharp
protected override void ProcessMessageImpl(IMessageHeader header, IModel model, BasicDeliverEventArgs basicDeliverEventArgs)
{
    ErrorAndNack(header,model,basicDeliverEventArgs,"Something went wrong", new Exception("What went wrong"));
}
```


## Class Diagram

![Class Diagram](Images/ClassDiagram.png)
﻿
FIXME: This is probably obsolete

## RabbitMQ Config Options

All the configuration settings are stored in a .yaml file. The default set of options is stored in the default.yaml file.

```yaml
    RabbitOptions:
        RabbitMqHostName: 'localhost'
        RabbitMqHostPort: 5672
        RabbitMqVirtualHost: '/'
        RabbitMqUserName: 'guest'
        RabbitMqPassword: 'guest'
        RabbitMqControlExchangeName: 'ControlExchange'

        FatalLoggingExchange: 'FatalLoggingExchange'
        
    FileSystemOptions:
        FileSystemRoot: 'C:\temp'
        ExtractRoot: 'C:\temp'
            
    RDMPOptions:
        Server: 'localhost\sqlexpress'
        CatalogueDb: 'RDMP_Catalogue'
        ExportDb: 'RDMP_DataExport'
    
    MongoDatabases:
        DicomStoreOptions:
            HostName: 'localhost'
            Port: 27017
            DatabaseName: 'dicom'
            SeriesCollection: 'series'
            ImageCollection: 'image'
        ExtractionStoreOptions:
            HostName: 'localhost'
            Port: 27017
            DatabaseName: 'extraction'

    DicomRelationalMapperOptions:
        Guid: '6ff062af-5538-473f-801c-ed2b751c7897'
        QueueName: 'AnonymousImageQueue'
        ConsumerTag: 'DicomRelationalMapper'
        QoSPrefetchCount: 10000
        AutoAck: false
        LoadMetadataId: 1
        DatabaseNamerType: 'GuidDatabaseNamer'

    CohortExtractorOptions:
        QueueName: 'RequestQueue'
        ConsumerTag: 'CohortExtractor'
        QoSPrefetchCount: 10000
        AutoAck: false
        AllCatalogues: true
        OnlyCatalogues: [1,2,3]
        # List of IDs of Catalogues to extract from (in ascending order).
        # Ignored if "AllCatalogues == true"
        #    - 2
        #    - 4
        #    - 5
        # also doable on a single line with [2,4,5] :)
        AuditorType: 'Smi.CohortExtractor.Audit.NullAuditExtractions'
        RequestFulfillerType: 'Smi.CohortExtractor.Execution.RequestFulfillers.FromCataloguesExtractionRequestFulfiller'
        # Writes (Producer) to this exchange
        ExtractFilesProducerOptions: 
            Label: 'ExtractFiles'
            ExchangeName: 'ExtractFileExchange'
        # And audits this too
        ExtractFilesInfoProducerOptions: 
            Label: 'ExtractFilesInfos'
            ExchangeName: 'FileCollectionInfoExchange'

    CohortPackagerOptions:
        JobWatcherTickrate: 30
        ExtractRequestInfoOptions: 
            QueueName: 'RequestInfoQueue'
            ConsumerTag: 'CohortPackager'
            QoSPrefetchCount: 1
            AutoAck: false
        ExtractFilesInfoOptions:
            QueueName: 'FileCollectionInfoQueue'
            ConsumerTag: 'CohortPackager'
            QoSPrefetchCount: 1
            AutoAck: false
        AnonImageStatusOptions:
            QueueName: 'FileStatusQueue'
            ConsumerTag: 'CohortPackager'
            QoSPrefetchCount: 1
            AutoAck: false

    DicomReprocessorOptions:
        ProcessingMode: 'ImageReprocessing'
        QueryScriptPath: '.\Scripts\abcd.txt'
        IdentImageProducerOptions: 
            Label: 'Image'
            ExchangeName: 'IdentifiableImageExchange'
        TagPromotionProducerOptions: 
            Label: 'TagPromotion'
            ExchangeName: 'TagPromotionExchange'

    DicomTagReaderOptions:
        QueueName: 'AccessionDirectoryQueue'
        ConsumerTag: 'TagReader'
        QoSPrefetchCount: 1
        AutoAck: false
        NackIfAnyFileErrors: true
        ImageProducerOptions: 
            Label: 'ImageProducer'
            ExchangeName: 'ImageExchange'
        SeriesProducerOptions: 
            Label: 'SeriesProducer'
            ExchangeName: 'SeriesExchange'

    IdentifierMapperOptions:
        QueueName: 'IdentifiableImageQueue'
        ConsumerTag: 'IdentifierMapper'
        QoSPrefetchCount: 1000
        AutoAck: false
        AnonImagesProducerOptions: 
            Label: 'AnonImages'
            ExchangeName: 'AnonymousImageExchange'
        MappingConnectionString: 'Server=localhost\sqlexpress;Integrated Security=true;Initial Catalog=MappingDatabase;'
        MappingDatabaseType: 'MicrosoftSQLServer'
        MappingTableName: 'MappingTable'
        SwapColumnName: 'CHI'
        ReplacementColumnName: 'ECHI'
        SwapperType: 'Smi.IdentifierMapper.Execution.ForGuidIdentifierSwapper'

    MongoDbPopulatorOptions:
        SeriesQueueConsumerOptions:
            QueueName: 'MongoSeriesQueue'
            ConsumerTag: 'MongoDBPopulator'
            QoSPrefetchCount: 1000
            AutoAck: false
        ImageQueueConsumerOptions:
            QueueName: 'MongoImageQueue'
            ConsumerTag: 'MongoDBPopulator'
            QoSPrefetchCount: 10000
            AutoAck: false
        MongoDbFlushTime: 1000
        FailedWriteLimit: 5

    ProcessDirectoryOptions:
        AccessionDirectoryProducerOptions:
            Label: 'AnonImages'
            ExchangeName: 'AccessionDirectoryExchange'
```

## How to use it

### C#

Simply use the static `Load()` method on the `GlobalOptions` class and you will receive a strongly typed object with all
of the above settings.

By default this will load the `default.yaml` file, but the `Load()` method has an overload which will accept a filename as a string (with or without the `.yaml` extension)
and will load it for you if it exists). THis is useful if you want to override the settings from the command line arguments.

You can then pass the GlobalOptions instance to your specific host in the bootstrapper sequence:

```csharp
    private static int Main(string[] args)
    {
        var options = new GlobalOptionsFactory().Load();
        Parser.Default.ParseArguments(args, cli);

        var options = new GlobalOptionsFactory().Load(cli.YamlFile);
            
        var bootstrapper = new MicroserviceHostBootstrapper(
            () => new DicomTagReaderHost(options));
        return bootstrapper.Main();
    }
```

### Java

To be completed...

### Development, test and production

There are multiple yaml files into the folder. During development work and for debugging it is better to simply edit the `default.yaml` to match your development 
environment.

During CI builds, the build script will replace the `default.yaml` contents with the `test.yaml` file, so it is important to keep them in sync if you
add or remove properties.

Production ready builds will replace the `default.yaml` with the `production.yaml` contents, although this file can also be edited after deployment and 
you can use command line arguments to load other files. Again, it is important to keep them in sync, mainly in relation to added properties to avoid 
exception while parsing.

## How to edit the yaml file

YAML syntax reference is abundant online. YAML is basically a sequence of key-value pairs separated by a colon `:`. Indentation is meaningful,
indented lines will be properties of the previous outdented line. Dashes indicate array elements, although you can use square brackets as well (see above).

If you want to add one option to an existing section, you need to create the corresponding property in the `GlobalOptions` class, as
everything is strongly typed and it will error if a property in the YAML file is not matched by the class (the viceversa instead is allowed).

If you create a new microservice, put its specific settings into a new top-level element in the YAML file.
﻿# IsIdentifiable Reviewer

Primary Author: [Thomas](https://github.com/tznind)

## Contents
 1. [Overview](#1-overview)
 2. [Setup / Installation](#2-setup--installation)
 3. [Usage](#3-usage)
    1. [Reviewing the output of IsIdentifiable]
    2. [Redacting the database]
    3. [Managing the rulebase]

## 1. Overview

Is Identifiable Reviewer is a cross platform text based UI for managing the anonymisation processes in which [PII] is detected and removed.  It serves as a management console for the
rulesbase of [IsIdentifiable] and as a downstream process for validating the results/redacting the database.

![Screenshot](./images/Role.png)
_The review process of potentially PII_

There are 3 activities that can be undertaken using the reviewer:

- [Reviewing the output of IsIdentifiable]
- [Redacting the database]
- [Managing the rulebase]

## 2. Setup / Installation

The application runs as a sub verb of `smi` (See [SmiRunner]). You can see the application help by running:

```
.\smi is-identifiable-reviewer --help
```

The following parts of the global yaml config file interact with this tool:

| YAML Section  | Purpose |
| ------------- | ------------- |
| IsIdentifiableOptions | What files/directories to find the current rulesbase for the [IsIdentifiable] process.  This used when [Managing the rulebase] |
| IsIdentifiableReviewerGlobalOptions |  |

| Command Line Options | Purpose |
| ------------- | ------------- |
|IsIdentifiableReviewerOptions | Allows overriding of which yaml file is loaded and runtime arguments e.g. which table to redact / report to load|

## 3. Usage

### Reviewing the output of IsIdentifiable

The IsIdentifiable tool applies NLP and the rules base to identify [PII] data in the database.  A sample output file is included: [ExampleReport](./ExampleReport.csv) is included.

Open the report using the `-f somefile.csv` command line option or `File->Open Report`.

Once loaded you can iterate the reports sequentially using the 'Sequential' tab or get an overview of all the issues encountered (aggregated by frequency) in the 'Tree View' tab.

![Screenshot](./images/Screenshot1.png)
_Sequential View_

The 'Sequential' view operates on one failure at a time. It shows the full string at the top, with the failures highlighted in green. At the bottom left is the classification of the failure: Person, Organisation, Date, etc. At the bottom right is the column (or DICOM tag) where the failure was found. It is important to check this column because, for example, you should Ignore a hospital name if the column is InstitutionName, but Update it if the column is StudyDescription.

The `Next` and `Prev` buttons move sequentially through the failures, i.e. `Next` does not skip over failures that are matched by existing rules.

![ScreenshotOfTree](./images/tree.png)
_'Tree View' showing PII detected aggregated by unique failing value and column where PII was found_

The 'Tree View' sorts all of the failures by number of occurrences. This tree view shows all the categories of rules and then all the categories of failures. It also shows the list of Conflicting rules which is where a failure matches both an Ignore and an Update rule.

Each instance of potential [PII] found by IsIdentifiable is termed a 'failure' (the existing anonymisation process has failed to strip this PII).  A 'failure' can be either a false positive or a genuine case of [PII].  Make a decision for each failure whether to ignore it or 'report' it.

#### Report/Ignore

Review the reports and mark either `Ignore` (this is a false positive) or `Update` (this is PII and needs to be redacted).  This will result in a new rule being added to either `NewRules.yaml` (Ignore) or `RedList.yaml` (Update).  Once  a rule is written it will be applied automatically to future reports loaded eliminating the lead to make duplicate decisions. After using `Ignore` or `Update` the display moves onto the next failure, skipping over those which are matched by existing rules.

Conceptually these rules are slightly different from the IsIdentifiable rules. IsIdentifiable first uses rules to spot known PII. Then it uses a NLP(NER) tool which attempts to find more PII. Finally it uses whitelist rules to ignore known false positives. Ideally these rules should be fine-tuned to reduce the work of the reviewer so, for example, if the reviewer shows 90% of failures are due to `Manufacturer=AGFA` it would be wise to manually edit IsIdentifiable rules. The Reviewer rules are different in that they are used filter the IsIdentifiable output and either ignore or redact its failure reports. The syntax of the rules files looks similar but is used differently, and has no effect on future runs of IsIdentifiable, only on future Reviews.

The menu `Options | Custom Patterns` menu, when ticked, will provide the opportunity to edit the Ignore/Update rule before it is saved. This allows you to make fine adjustments to the exact pattern which will be redacted. Note that all bracketed patterns are redacted so you can add (or remove) any as necessary. For example, if the full string is `John Smith Hospital^MRI Head^(20/11/2020)` but only the date has been detected you could still redact the hospital name as well by editing the pattern to be `(John Smith Hospital)^.*^\((\d\d/\d\d/\d\d\d\d)\)$` (i.e. adding the name in brackets).

The Custom Patterns window provides several options to edit the pattern:

* `x` - clears currently typed pattern
* `F` - creates a regex pattern that matches the full input value
* `G` - creates a regex pattern that matches only the failing part(s)
* `\d` - replaces all digits with regex wildcards
* `\c` - replaces all characters with regex wildcards
* `\d\c` - replaces all digits and characters with regex wildcards

### Redacting the database

Once all 'failures' in a report have been processed and either ignored or a 'report' rule generated you can redact the database.  This is done by running the application using the `-u` and `-t` flags.

Since you may have several servers / databases that are processed using this tool, it is necessary to indicate where UPDATE commands should be run.  This is done by putting the connection string in a 'targets' file:

```yaml
- Name: My Server
  ConnectionString: Server=localhost;Username=root;Password=zombie
  DatabaseType: MySql
```
_Example targets file_

The following flags should be combined to successfully redact the database:


| Flag | Example | Purpose |
| ------------- | ------------- |------------- |
| -f | -f ./ExampleReport.csv | Indicates which [IsIdentifiable] output report to redact.  You must have completed the [review process] for this report|
| -u | -u ./misses.csv  | Indicates that you want to update the database.  The file value must be included and is where reports that are not covered by rules generated in the [review process] are output.  If you have completed the [review process] correctly this file should be empty after execution completes |
| -t | -t z:\temp\targets.yaml | Path to a file containing the connection string (and DMBS type) of the relational database server that has the table requiring redaction|

```bash
smi.exe is-identifiable-reviewer -f ./ExampleReport.csv -u ./misses.csv -t z:\temp\targets.yaml
```
_Example redaction command_

### Managing the rulebase

Over time the number of rules in [IsIdentifiable] and the reviewer will increase.  It can be beneficial to move ignore rules upstream from the reviewer to the [IsIdentifiable] rulebase especially for commonly encountered reports.  This will reduce the number of false positives and the size of report files.

The 'Rules Manager' tab provides visualisation and control over the rules used by the [IsIdentifiable] tool ('Analyser Rules') and the Reviewer ('Reviewer Rules').  You should periodically review the rules base to ensure there are no mistakes and to identify candidates for pushing upstream into the analyser.

![RulesManagerScreenshot](./images/rulesmanager.png)
_Rules Manager View_
 
| Key | Function |
| ------------- | ------------- |
| `<Delete>` | Removes a rule from the rulesbase |
| `<Enter>` | Opens menu (if any) for interacting with rule(s) highlighted |

[IsIdentifiable]: ../../microservices/Microservices.IsIdentifiable/README.md
[PII]: https://en.wikipedia.org/wiki/Personal_data
[SmiRunner]: ../Applications.SmiRunner/
[Managing the rulebase]: #managing-the-rulebase
[review process]: #reviewing-the-output-of-IsIdentifiable
[Reviewing the output of IsIdentifiable]: #reviewing-the-output-of-isidentifiable
[Redacting the database]: #redacting-the-database
Queues a set of images for extraction based on the dicom UIDs in a file# SmiRunner

This CLI application serves as the single entry point for all other C# applications and servcices in the SmiServices repo.

## Usage

Running `--help` will show a list of all the services, and a link to the main README for that service

```console
$ ./smi --help
smi 1.15.1
Copyright  SMI Project 2018-2020

  trigger-updates              See here at your release version:
                               https://github.com/SMI/SmiServices/tree/master/src/applications/Applications.TriggerUpdates

  ...
```

Each service can then be run the same as previously by passing its own specific set of parameters or verbs

```console
$ ./smi dicom-tag-reader -y foo.yaml
...
```

```console
$ ./smi is-identifiable db -y default.yaml ...
...
```

## Supporting a new app or service

In order to add a new app or service to the runner, first create the csproj as normal in the appropriate "Applications" or "Microservices" directory/namespace.

It might be useful to first create this as a Console project, with a main entrypoint etc. for ease of initial testing. Once it's stable, change it to a library project by specifying `<OutputType>Library</OutputType>` in the csproj file.

To add the project to SmiRunner, the process is then to

-   Add reference to the project in the SmiRunner csproj
-   Add a new [ServiceVerb](./ServiceVerbs.cs)
-   Add the new verb to [Program](./Program.cs)
    -   Add to one of the static arrays `AllServices` or `AllApplications`
    -   Add a case statement which points to the program entry point
# SRAnonTool

The SRAnonTool is a set of programs to assist with anonymising DICOM Structured Reports.

## Requirements

Python package requirements:
* pydicom
* pymongo
* deepmerge

Python library requirements:
* The SmiServices python library, see `src/common/Smi_Common_Python/`

External tool requirements:
* The SmiServices CTP anonymiser
* SemEHR/CogStack anonymiser (or the test stub)
* dcm2json (for testing; optional; from the dcmtk package)
* jq (for testing; optional)
* diff (for testing)

## Installation

The scripts require dependencies which need to be found via the `PATH` or `PYTHONPATH`.

Copy the scripts to `$SMI_ROOT/scripts`

Install the python library `SmiServices` to `$SMI_ROOT/lib/python3/` or a virtualenv

Ensure the python package dependencies are installed system-wide or in a virtualenv on the host machine.

Modify the `default.yaml` file: in the section `CTPAnonymiserOptions` add `SRAnonTool: /path/to/SRAnonTool.sh`

Ensure the `default.yaml` file contains the necessary `FileSystemOptions`, `LoggingOptions>LogsRoot`, `MongoDatabases`, `RabbitOptions`, etc.

Install the SemEHR/CogStack anonymiser, which currently uses the following directories (which may be symbolic links):

* `/opt/semehr/CogStack` - contains the scripts
* `/opt/semehr/data/input_docs` - raw text copied here will be anonymised
* `/opt/semehr/data/anonymised` - output anonymous text and xml files
* `/opt/gcp/bio-yodie-1-2-1` - the UMLS dictionary
* `/opt/gcp/gcp-2.5-18658` - java libraries

The old SemEHR anonymiser requires Python2; all the other scripts require Python3, including the new SemEHR anonymiser.
If using the test stub then only the data directories are required and Python2 is not required.

## Usage as part of CTP

Configure CTP to call the script SRAnonTool.sh when it detects a DICOM file with `SR` in the `Modality` tag, by editing `default.yaml` as above. CTP will call the script with two options:
* `-i input.dcm` - the raw DICOM file before anonymisation
* `-o output.dcm` - the DICOM file which CTP has already anonymised

The script will extract the text from the `input.dcm` file, anonymise it, and write the redacted text into the `output.dcm` file, which must already exist.

## Standalone usage

The script `SRAnonTool.sh` calls three components:

* `CTP_DicomToText.py` - extracts the text from the raw DICOM file into a format suitable for SemEHR-CogStack.
* `clinical_doc_wrapper.py` - this is the component within SemEHR-CogStack which anonymises the text.
* `CTP_XMLToDicom.py` - redacts the text from the raw DICOM file and write the redacted text into the output DICOM file.

Usage: `[-e virtualenv] [-s semehr_dir]  -i read_from.dcm  -o write_into.dcm`

The `-e` option can be used to activate a virtual environment by sourcing the `bin/activate` script.
The default is to set `PYTHONPATH` to `$SMI_ROOT/lib/python3`

The SemEHR directory (`/opt/semehr`) can be changed with the `-s` option for testing (it's not set when called by CTP).

### `CTP_DicomToText.py`

This program can be used as part of the SRAnonTool pipeline or it can be used standalone to extract documents in bulk for later SemEHR processing.

Usage: `-y default.yaml -i input.dcm -o outfile [--semehr-unique]`

`-y default.yaml` - may be specified more than once if the configuration parameters are spread across multiple yaml files.

`-i input.dcm` - full path to the input DICOM file, or a partial path to be extracted from MongoDB, or a StudyDate to extract all records that day from MongoDB.

`-o output` - full path to the output text file, or directory for multiple files.

`--semehr-unique` - if extracting a StudyDate from MongoDB then ignore any documents which have a SOPInstanceUID that is already in the SemEHR MongoDB database. This is intended to allow reprocessing of any documents that previously failed without having to reprocess the whole day.

The MongoDB configuration read from the yaml files needs to be in `MongoDatabases | DicomStoreOptions` and `SemEHRStoreOptions`. The former is to read DICOM documents from the `dicom.image_SR` database.collection; the latter is to check if the SOPInstanceUID is already in the `semehr.semehr_results` database.collection.

Examples:

```
* CTP_DicomToText.py -i /path/to/file.dcm -o output.txt
* CTP_DicomToText.py -i 2015/01/01/AccNum/file.dcm -o output.txt -y smi_dataLoad.yaml
* CTP_DicomToText.py -i 20150101 -o output_dir -y smi_dataLoad.yaml
```

### `clinical_doc_wrapper.py`

This script performs the anonymisation.

Usage: `[-s semehr_dir] [-i input_docs] [-o anonymised]` in the stub version

It must be called with the current directory being the location of the script.

It reads all the files in the `/data/input_docs` directory. For each input file it write a slightly modified file with the same name into the `/data/anonymised` directory, basically the text with some header metadata removed, plus it writes a file of the same name plus `.knowtator.xml` appended containing annotations in XML format.

The SemEHR version requires Python2; the test stub requires Python3. Note that this script is no longer used in the new SemEHR anonymiser.

The test stub of this program has no requirement on current directory. It is best suited when tested with the given test DICOM file as it only fakes the anonymisation of the word `Baker`. The `-s` option can be used to specify the SemEHR directory instead of `/opt/semehr` which is useful when testing; this option is not present in the original.

## `anonymiser.py`

This script performs the anonymisation.

Usage: `./anonymiser.py /path/to/anonymisation_task.json`

It must be called with the current directory being the location of the script.

The template configuration file is typically in the `anonymisation/conf` directory
but should have the following elements modified for each run of the anonymiser:

```
.text_data_path=${semehr_input_dir}
.anonymisation_output=${semehr_output_dir}
.extracted_phi=${semehr_output_dir}/phi
.grouped_phi_output=${semehr_output_dir}/phi_grouped
.logging_file=${semehr_output_dir}/log
.annotation_mode=false # temporary false until the knowtator XML output is fixed
```

### `CTP_XMLToDicom.py`

Usage: `-y default.yaml -i input.dcm -x input.xml -o output.dcm`

`-y default.yaml` - may be specified more than once if the configuration parameters are spread across multiple yaml files

`-i input.dcm` - full path to the raw DICOM file

`-x input.xml` - full path to the XML file containing annotations

`-o output.dcm` - full path to the anonymised DICOM file, which must already exist, where the redacted text is written


## Testing

In the test subdirectory, run

```
./CTP_SRAnonTool_test.py [-s semehr_dir] [-d file.dcm] [-p pattern_to_redact] [-y default.yaml]
```

That will read `report10html.dcm` and run the above scripts, checking that the output matches what is expected.
It will print `SUCCESS` and exit 0 if successful, exit 1 if the output is not as expected.

The defaults are:

`-s semehr_dir` - `/opt/semehr`

`-d file.dcm` - `report10html.dcm` (a public sample file manually edited to include HTML fragments)

`-p pattern_to_redact` - `Baker` (to suit the example DICOM file)

`-y default.yaml` - `../../../../data/microserviceConfigs/default.yaml`
# TriggerUpdates

Primary Author: [Thomas Nind](https://github.com/tznind)

## Contents
 1. [Overview](#1-overview)
 2. [Setup / Installation](#2-setup--installation)
 3. [Exchange and Queue Settings](#3-exchange-and-queue-settings)
 4. [Config](#4-config)
 5. [Expectations](#5-expectations)
     1. [Aliases](#51-aliases)
 6. [Class Diagram](#6-class-diagram)

### 1. Overview
The TriggerUpdates app is a console app that runs and terminates rather than a microservice that runs forever. The app issues update messages designed for consumption by the [UpdateValues microservice] .

This application uses the verb system e.g. `./TriggerUpdates mapper` or `./TriggerUpdates mongo`.

### 2. Setup / Installation
 - Clone the project and build. Any NuGet dependencies should be automatically downloaded
 - Edit the yaml.default with the configuration for your environment
 - Run `TriggerUpdates.exe mapper -d <somedate> -f PatientID` from a commandline to issue messages for updating the live database(s) with changed mappings.

### 3. Exchange and Queue Settings
| Read/Write | Type | Config setting |
| ------------- | ------------- |------------- |
| Write | `UpdateValuesMessage` | `TriggerUpdatesOptions` |

### 4. Config
| YAML Section  | Purpose |
| ------------- | ------------- |
| RabbitOptions | Describes the location of the rabbit server for sending messages to |
| IdentifierMapperOptions | Describes the location of the mapping table when using the `mapper` verb with this command. |
| TriggerUpdatesOptions | The exchange name that `UpdateValuesMessage` should be sent to when detecting updates |

Arguments vary by verb, use `./TriggerUpdates mapper --help` to see required arguments.

### 5. Expectations
Errors are [logged as normal for a MicroserviceHost](../../common/Smi.Common/README.md#logging)

#### 5.1 Aliases

When using mapping updates it is possible for certain corner case sequences to result in crossed mappings, especially when aliases are permitted and where those aliases change over time

Initial Lookup Table

| Private | Release|
|---|---|
| A | 111|
| B | 222|

_Any live database values for A or B will have only the Release identifiers (111 and 222)_

Update lookup table with fact A=B

| Private | Release|
|---|---|
| A | 111|
| B | 111|

_Triggering an update at this point will merge 111 and 222 in the live database_

Once an alias has been established the lookup cannot successfully be updated to reverse the alias e.g. reverting it back to the initial state.


### 6. Class Diagram
![Class Diagram](./Images/ClassDiagram.png)

[UpdateValues microservice]: ../../microservices/Updating/Microservices.UpdateValues/README.md
# ProcessDirectory

Primary Author: [Ally Hume](https://github.com/allyhume)

## Contents
 1. [Overview](#1-overview)
 2. [Setup / Installation](#2-setup--installation)
 3. [Exchange and Queue Settings](#3-exchange-and-queue-settings)
 4. [Config](#4-config)
 5. [Expectations](#5-expectations)
 6. [Class Diagram](#6-class-diagram)
 7. [Directory Scan Modes](#7-directory-scan-modes)

### 1. Overview
The ProcessDirectory app is a console app that runs and terminates rather than a microservice that runs forever. The app searches in and below a specified directory to find all directories that contain DICOM files. Each directory found which contains at least 1 dicom file will result in an `AccessionDirectoryMessage` being created. The behaviour of the scan can be changed with the `-f` option (see [#7](#7-directory-scan-modes)). When all the directories below the specified directory have been processed the program terminates.

### 2. Setup / Installation
 - Clone the project and build. Any NuGet dependencies should be automatically downloaded
 - Edit the yaml.default with the configuration for your environment
 - Run `ProcessDirctory.exe` from a commandline with the top level directory as the only argument.

### 3. Exchange and Queue Settings
| Read/Write | Type | Config setting |
| ------------- | ------------- |------------- |
| Write | AccessionDirectoryMessage | `ProcessDirectoryOptions.AccessionDirectoryProducerOptions` |

### 4. Config
| YAML Section  | Purpose |
| ------------- | ------------- |
| RabbitOptions | Describes the location of the rabbit server for sending messages to |
| FileSystemOptions | Describes the root location of all images this program will ever load, all directories provided in command line arguments must be subdirectories of this root in order that file names can be expressed as relative paths for downstream microservices. |
| ProcessDirectoryOptions | The exchange name that `AccessionDirectoryMessage` should be sent to when finding subdirectories with dicom files in them |

| Argument | Command Line Options | Purpose |
| ------------- | ------------- | ------------- |
|-d| Directory root | (Required) The directory of the image archive in which to begin folder discovery |
|-f| Directory search format |(Optional) The directory scan mode to use. See [] |

### 5. Expectations
Errors are [logged as normal for a MicroserviceHost](../../common/Smi.Common/README.md#logging)

### 6. Class Diagram
![Class Diagram](./Images/ClassDiagram.png)

### 7. Directory Scan Modes

Specified by the `-f` CLI argument. Options are:

- [Default](Execution/DirectoryFinders/BasicDicomDirectoryFinder.cs). Performs a general recursive scan for files. The program does not look inside any subdirectories of directories which contain at least 1 file. This scan mode searches for DICOM files with the extension specified in the `FileSystemOptions.DicomSearchPattern` config option.

- [PACS](Execution/DirectoryFinders/PacsDirectoryFinder.cs). Performs a scan which assumes files are located inside a particular directory structure. The `PACS` directory structure is of the form `<any root>/YYYY/MM/DD/ACC/<dicom>`. `ACC` represents accession directories. For each directory found which matches this pattern, an `AccessionDirectoryMessage` is produced. Note that this scan mode does not actually assert that there are any files inside the accession directories.

- [List](Execution/DirectoryFinders/AccessionDirectoryLister.cs). Receives a file containing a list of accession directory paths. The accession directory path structure is expected to be of the form `<any root>/YYYY/MM/DD/ACC/`. For each path that matches this pattern, the existence of the directory and whether it contains DICOM files is checked. If the path meets all of the requirements, an `AccessionDirectoryMessage` is produced. Note that the input file path is passed in the same way as directory paths are passed to other operational modes.
# Cohort Packager

Primary Author: [Ruairidh MacLeod](https://github.com/rkm)

## Contents

 1. [Overview](#1-overview)
 2. [Setup / Installation](#2-setup--installation)
 3. [Queue Settings](#3-queue-settings)
 4. [Config](#4-config)
 5. [Expectations](#5-expectations)
 6. [Reports](#6-reports)

### 1. Overview

Collects all information regarding an extraction job, and monitors the filesystem for the anonymised files. Persists all information to a MongoDB collection.

Produces validation reports for each extraction suitable for review by research coordinators before the extraction files are released. See [reports section](#6-reports). Reports are created automatically when an extraction is detected as being complete, and can also be manually recreated on the CLI by passing the `-r` or `--recreate-reports` flag with the corresponding extraction GUID.

### 2. Setup / Installation

- Clone the project and build. Any NuGet dependencies should be automatically downloaded
- Setup a yaml file with the configuration for your environment
- Run `CohortPackager.exe` with your yaml config

### 3. Exchange and Queue Settings

| Read/Write | Type | Config setting |
| ------------- | ------------- |------------- |
|Read|ExtractRequestMessage|`DicomReprocessorOptions.ExtractRequestInfoOptions`|
|Read|ExtractRequestInfoMessage|`DicomReprocessorOptions.ExtractFilesInfoOptions`|
|Read|ExtractFileStatusMessage | `DicomReprocessorOptions.AnonImageStatusOptions` |

### 4. Config

| YAML Section  | Purpose |
| ------------- | ------------- |
|JobWatcherTickrate|How often the filesystem is checked for anonymised files (in seconds)|

### 5. Expectations

Errors are [logged as normal for a MicroserviceHost](../../common/Smi.Common/README.md#logging)

### 6. Reports

When an extraction is completed, report(s) are created detailing any errors or validation failures relating to the set of files that have been produced. There are currently 2 formats:
-   `Combined`. A single file useful for smaller extractions which can be manually reviewed. Format:
    ```markdown
    # SMI extraction validation report for testProj1/extract1

    Job info:
    -   Job submitted at:             2020-11-19T17:48:59
    -   Job completed at:             2020-11-19T17:49:07
    -   Job extraction id:            62c1b363-8181-435c-9eca-11e1f095043f
    -   Extraction tag:               SeriesInstanceUID
    -   Extraction modality:          Unspecified
    -   Requested identifier count:   2
    -   Identifiable extraction:      No
    -   Filtered extraction:          Yes

    Report contents:

    -   Verification failures
        -   Summary
        -   Full Details
    -   Blocked files
    -   Anonymisation failures

    ## Verification failures

    ### Summary

    -   Tag: ScanOptions (2 total occurrence(s))
        -   Value: 'FOO' (2 occurrence(s))


    ### Full details

    -   Tag: ScanOptions (1 total occurrence(s))
        -   Value: 'FOO' (1 occurrence(s))
            -   series-2-anon-2.dcm
            -   series-2-anon-3.dcm


    ## Blocked files

    -   ID: series-1
        -   1x 'rejected - blah'

    ## Anonymisation failures

    -   file 'series-2-anon-1.dcm': 'Couldn't anonymise'

    --- end of report ---
    ```

-   `Split`. A pack (directory) of files which can be used for larger extractions to script the review process by reading them into external analysis tools. 6 files are currently produced for each extraction:
    -   `README.md`. Summary information for the extract
        ```markdown
        # SMI extraction validation report for testProj1/extract1

        Job info:
        -   Job submitted at:             2020-11-19T17:50:41
        -   Job completed at:             2020-11-19T17:50:48
        -   Job extraction id:            08f22a41-2a3d-4fa3-b9d2-380aabc5a59b
        -   Extraction tag:               SeriesInstanceUID
        -   Extraction modality:          Unspecified
        -   Requested identifier count:   2
        -   Identifiable extraction:      No
        -   Filtered extraction:          Yes

        Files included:
        -   README.md (this file)
        -   pixel_data_summary.csv
        -   pixel_data_full.csv
        -   pixel_data_word_length_frequencies.csv
        -   tag_data_summary.csv
        -   tag_data_full.csv

        This file contents:
        -   Blocked files
        -   Anonymisation failures

        ## Blocked files

        -   ID: series-1
            -   1x 'rejected - blah'

        ## Anonymisation failures

        -   file 'series-2-anon-1.dcm': 'Couldn't anonymise'

        --- end of report ---
        ```

    -   `pixel_data_full.csv`. Full details for all potential PII contained in the DICOM pixel data from the post-anonymisation OCR scan
        ```csv
        TagName,FailureValue, FilePath
        PixelData,Mr. Foobar,path/to/file.dcm
        ```

    -   `pixel_data_summary.csv`. Summary details for all potential PII contained in the DICOM pixel data from the post-anonymisation OCR scan
        ```csv
        TagName,FailureValue,Occurrences,RelativeFrequencyInTag,RelativeFrequencyInReport
        PixelData,Mr. Foobar,1,1,1
        ```

    -   `pixel_data_word_frequencies.csv`. Frequency analysis of the word lengths contained in the DICOM pixel data from the post-anonymisation OCR scan
        ```csv
        WordLength,Count,RelativeFrequencyInReport
        10,1,1
        ```

    -   `tag_data_full.csv`. Full details for all potential PII contained in the DICOM tag data from the post-anonymisation NER scan
        ```csv
        TagName,FailureValue,FilePath
        ScanOptions,FOO,series-2-anon-2.dcm
        ```

    -   `tag_data_summary.csv`. Summary details for all potential PII contained in the DICOM tag data from the post-anonymisation NER scan
        ```csv
        TagName,FailureValue,Occurrences,RelativeFrequencyInTag,RelativeFrequencyInReport
        ScanOptions,FOO,1,1,1
        ```

Note that for identifiable extractions (where no anonymisation is applied), only the combined report is supported. This will have the format:

```markdown
# SMI extraction validation report for testProj1/extract1

Job info:
-   Job submitted at:             2020-11-19T17:57:12
-   Job completed at:             2020-11-19T17:57:20
-   Job extraction id:            3f400e06-19cb-45ae-9545-3a5310f426f3
-   Extraction tag:               StudyInstanceUID
-   Extraction modality:          MR
-   Requested identifier count:   1
-   Identifiable extraction:      Yes
-   Filtered extraction:          Yes

Report contents:

-   Missing file list (files which were selected from an input ID but could not be found)

## Missing file list

-   study-1-orig-2.dcm

--- end of report ---
```# MongoDB Populator

Primary Author: [Ruairidh](https://github.com/rkm)

## Contents
 1. [Overview](#1-overview)
 2. [Setup / Installation](#2-setup--installation)
 3. [Exchange and Queue Settings](#3-exchange-and-queue-settings)
 4. [Config](#4-config)
 5. [Expectations](#5-expectations)
 6. [Class Diagram](#6-class-diagram)
 
### 1. Overview
Stores the Dicom Tag data in `DicomFileMessage` and/or `SeriesMessage` into a MongoDB database document store.

### 2. Setup / Installation

 - Install [MongoDB](https://docs.mongodb.com/manual/installation/), steps will be specific to your environment
 - Optional: Install a GUI such as [Compass](https://www.mongodb.com/products/compass) to easily work with MongoDB
 - Clone the project and build. Any NuGet dependencies should be automatically downloaded
 - Edit the yaml.default with the configuration for your environment
 - Ensure MongoDB is running and run MongoDBPopulator.exe from a commandline

### 3. Exchange and Queue Settings
| Read/Write | Type | Config setting |
| ------------- | ------------- |------------- |
| Read | DicomFileMessage| MongoDbPopulatorOptions.ImageQueueConsumerOptions.QueueName |
| Read | SeriesMessage| MongoDbPopulatorOptions.SeriesQueueConsumerOptions.QueueName |

### 4. Config

| YAML Section  | Purpose |
| ------------- | ------------- |
| RabbitOptions | Describes the location of the rabbit server for sending messages to |
| MongoDatabases | Contains the connection strings and database names to write tag data to |
| MongoDbPopulatorOptions | Queue names to read messages from, error threshold and how often to push to mongo |

| Command Line Options | Purpose |
| ------------- | ------------- |
|CliOptions | Allows overriding of which yaml file is loaded. |

### 5. Expectations
Errors are [logged as normal for a MicroserviceHost](../../common/Smi.Common/README.md#logging)

#### Data Failure States
- Operation on receiving corrupt message.
	 - Program will send Nack with `requeue` flag set to `false`. This will indicate to RabbitMQ to send the message to the dead letter exchange for analysis.

#### Environmental Failure States
  - Operation on loss of MongoDB connection:
	 - Program will attempt to reconnect and send a warning message every \<x\> minutes until the connection is recovered. If no connection found after \<y\> minutes, it will Nack all messages it has received and enter a paused state, sending a `fatal` error or similar.
 

### 6. Class Diagram
![Class Diagram](./Images/ClassDiagram.png)
# DicomRelationalMapper

Primary Author: [Thomas](https://github.com/tznind)

## Contents
 1. [Overview](#1-overview)
 2. [Setup / Installation](#2-setup--installation)
 3. [Queue Settings](#3-queue-settings)
 4. [Config](#4-config)
 5. [Expectations](#5-expectations)
 6. [Adding Tags](#6-adding-tags)
 7. [Audit](#7-audit)
 8. [Class Diagram](#8-class-diagram)
 
 
### 1. Overview
Runs an RDMP data load configuration (`LoadMetadata`) with a batch of `DicomFileMessage` to load Dicom Tag data into a relational database (MySql or Microsoft Sql Server).  It is designed to work on many images at once for performance.  The load configuration is configurable through the main RDMP client.

### 2. Setup / Installation
#### 2.1 RDMP / Imaging dataset load setup
- [Install and setup RMDP](https://github.com/HicServices/RDMP#install)
- [Install Rdmp.Dicom plugin and setup imaging tables](https://github.com/HicServices/RdmpDicom/tree/develop#using-plugin) (Ensure Json sources is selected for the pipeline)
- Clone this repository and build. Any NuGet dependencies should be automatically downloaded
- Edit the yaml.default with the configuration for your environment
- Run DicomRelationalMapper

#### 2.2 Tests setup
- In order to run database tests (and after each RDMP platform API schema change) you will need to [create TEST_ platform databases](https://github.com/HicServices/RDMP/blob/master/Documentation/CodeTutorials/Tests.md).  DatabaseCreation.exe will be in the bin directory of any Test project referencing HIC.RDMP.Plugin.Test.

### 3. Queue Settings
| Read/Write | Type | Config setting |
| ------------- | ------------- |------------- |
| Read | DicomFileMessage| DicomRelationalMapperOptions.QueueName |

### 4. Config
| YAML Section  | Purpose |
| ------------- | ------------- |
| RabbitOptions | Describes the location of the rabbit server for sending messages to |
| RDMPOptions | Describes the location of the Microsoft Sql Server RDMP platform databases which keep track of load configurations, available datasets (tables) etc |
| DicomRelationalMapperOptions | The queue from which to read DicomFileMessage and the ID of the `LoadMetadata` load configuration.  A load configuration is a sequence of steps to modify/clean data such that it is loadable into the final live tables.  The LoadMetadata is designed to be modified through the RMDP user interface and is persisted in the LoadMetadata table (and other related tables) of the RDMP platform database |

| Command Line Options | Purpose |
| ------------- | ------------- |
|CliOptions | Allows overriding of which yaml file is loaded. |

### 5. Expectations
DicomRelationalMapper is expected to be robust and able to load all Dicom files it is given as input (The input is anonymised images from another service so we know the images are processable).  The load logic relies on the LoadMetadata configuration which is customisable through the RDMP application.

Errors are [logged as normal for a MicroserviceHost](../../common/Smi.Common/README.md#logging)

#### Data Failure States
- Dicom file that is not a valid .dcm file (`TestLoadingOneImage_SingleFileMessage`)
	- Ignores the image file and processes rest of directory
- Empty Directory
	- All invalid messages are automatically Nacked with reque flag false (`TestDodgyDirectories`)
- Dicom data corrupt (missing primary keys, tag conflicts etc)
	- Data batch is failed, exe stops accepting new messages until someone has resolved the problem*
	
*This is the current system behaviour because we are assuming the dicoms are all relatively processable and that we can build in a system for resolving collisions that makes sense in the load configuration relatively easily.  Therefore any problems that arise during load that cause it to crash out should be rare and require fixes to be applied across the board.

#### Environmental Failure States
 - Operation on loss of RabbitMQ connection during a load:
	 - In progress data loads will complete and process will crash on BasicAck (leaving messages unacknowledged).  Repeat processing of stale messages (by other instances) will not introduce duplication.
 - Operation on loss of Microsoft SQL Server connection during a load:
	 - Data load batch will fail and all messages will be Nacked.  Live data integrity will be maintained since final load step is applied in a transaction (Staging=>Live merge).
 - Operation on data load error (e.g. corrupt files / file system down).
	 - Data load batch will fail and all messages will be Nacked.  Raw/Staging will be left available for diagnostics/debugging
	 
### 6. Adding Tags
Relational database tables have an initial schema out of the box based on a [image template](https://github.com/HicServices/RdmpDicom/blob/develop/Documentation/DataLoad.md#Image-Tables).

Once the system has gone live, data analysts will still be able to add new tags to existing database tables through the RDMP user interface using [tag promotion](https://github.com/HicServices/RdmpDicom/blob/develop/Documentation/DataLoad.md#Image-Tables).

### 7. Audit
In addition to logging to NLog like other microservices, the data load itself will be audited in the RDMP relational logging database.  This includes facts such as how many records were loaded (UPDATES / INSERTS) and any problems encountered.

![Load Metadata Logs](Images/LoadMetadataLogst.png)

In order to improve traceability the 'image' table of every set of imaging database tables (Study + Series + Image) has a field `messageguid` which contains the Guid of the input message to DicomRelationalMapper that resulted in the record being part of the batch.  This guid can be traced back through the file [logs](../Microservices.Common/README.md#logging) to see all microservices that acted to result in that record being in the final live database.

Finally all tables (Study, Series and Image) have a validFrom and a dataLoadRunId field which record which load batch they were last part of and when it was executed.

In the event that a load fails (e.g. due to primary key collisions) the RAW / STAGING databases that contain data being worked on in the load are left intact for debugging.  In such a case all messages in the batch are Nacked.

![Extra Fields](Images/LoadBubbles.png)

### 8. Class Diagram
![Class Diagram](./Images/ClassDiagram.png)
# CTPAnonymiser

Primary Author: [Paul Graham] (https://github.com/pjgraham)

## Contents
 1. [Overview](#1-overview)
 2. [Setup / Installation](#2-setup-installation)
 3. [Queue Settings](#3-exchange-and-queue-settings)
 4. [Config](#4-config)
 5. [Expectations](#5-expectations)

### 1. Overview

The CTP Anonymiser app is a RabbitMQ-based tool for anonymising DICOM image data. Launched from the command line, it consumes from a RabbitMQ queue, receiving directory names of DICOM files, anonymises them using the CTP DICOM [Anonymizer](https://mircwiki.rsna.org/index.php?title=The_CTP_DICOM_Anonymizer), and then produces the output anonymised filenames to another RabbitMQ queue. The choice of fields etc to anonymise is determined by a provided anonymisation script file.

### 2. Setup / Installation

The anonymiser is installed via Maven as per the other Java apps, so clone the project from git and run the maven install (see [here](https://github.com/SMI/SmiServices/blob/master/src/common/com.smi.microservices.parent/README.md) for details).

### 3. Exchange and Queue Settings

| Read/Write | Type | Config setting |
| ------------- | ------------- |------------- |
| Read| ExtractFileMessage | `CTPAnonymiserOptions.AnonFileConsumerOptions` |
| Write| ExtractFileStatusMessage|`CTPAnonymiserOptions.ExtractFileStatusProducerOptions`|

The ExtractFileMessage contents include the name of a directory of DICOM files to be anonymised. As the files are anonymised, ExtractFileStatusMessage messages are produced indicating success or otherwise.

### 4. Config

As for the other Java and C# microservices, the anonymiser uses the yaml config files documented [here](https://github.com/SMI/SmiServices/blob/master/src/common/Smi.Common/Options/RabbitMqConfigOptions.md).

In particular, the anonymiser uses:

* RabbitOptions - location of the RabbitMQ host etc
* FileSystemOptions - Root directories where files will be discovered/anonymised to
* CTPAnonymiserOptions - the config for the anonymiser consumer and producer

| CLI Options| Switch | Required | Purpose |
| :----- | :-------: |:----: |:---|
|Yaml config|-y| No| Allows overriding of which yaml file to use|
|Anon script|-a| Yes | The anonymisation script file to use |

The anonymisation script is based on the example provided by CTP, but it has been further restricted, for example, to exclude all fields which *may* include user-defined text.

The current anonymisation script can be viewed [here](https://github.com/SMI/SmiServices/blob/master/data/ctp/ctp-whitelist.script), file `dicom-whitelist.script`

Note that the CTP also supports pixel anonymisation, but at the moment this is not being exploited.

### 5. Expectations

As each DICOM file is processed, a RabbitMQ status [message](https://github.com/SMI/SmiServices/blob/master/src/microservices/com.smi.microservices.ctpanonymiser/src/main/java/org/smi/ctpanonymiser/messages/ExtractedFileStatusMessage.java) is produced, indicating success or [otherwise](https://github.com/SMI/SmiServices/blob/master/src/microservices/com.smi.microservices.ctpanonymiser/src/main/java/org/smi/ctpanonymiser/util/ExtractedFileStatus.java) of the anonymisation. Unsuccessful anonymisation attempts are tagged to be either retried, or  not retried.

The status message also includes details of the path to the anonymised file, project number, job id etc. 

### 6. Structured Report anonymisation

The procedure for anonymising Structured Reports (the SR modality) is to have
CTP proceed as normal, removing almost all tags which may contain PII, then
extract the text from the original DICOM file and pass it through an external
anonymisation process, and then insert the redacted text into CTP's output file,
thus reconstructing the text content part of the DICOM with anonymous text.

To do this requires an external tool to be configured as `SRAnonTool` in the yaml,
which should be the full path to the tool, unless it is in the current directory.
The tool is called with `-i input.dcm -o output.dcm`, where input.dcm is the raw
original DICOM file and output.dcm is the output from CTP. The tool is only called
for DICOM files which have the `Modality` tag equal to `SR`.

Any failures in the execution of the tool will result in a complete failure of
the CTP anonymisation process for this file. Even though CTP may have successfully
removed all PII from the file, it will not proceed if it cannot insert redacted text.
# IdentifierMapper

Primary Author: [Thomas](https://github.com/tznind)

## Contents
 1. [Overview](#1-overview)
 2. [Setup / Installation](#2-setup--installation)
 3. [Exchange and Queue Settings](#3-exchange-and-queue-settings)
 4. [Config](#4-config)
 5. [Expectations](#5-expectations)
 6. [Class Diagram](#6-class-diagram)

### 1. Overview
This service takes serialized Dicom file as `DicomFileMessage` messages and uses an `ISwapIdentifiers` to replace the top level `DicomTag.PatientID` tag for an anonymous representation.  If there is no PatientID then the message is nacked. If the `ISwapIdentifiers` returns null then the message is nacked with the reason provided by `string GetSubstitutionFor(string toSwap, out string reason)`.

### 2. Setup / Installation
- Clone the project and build. Any NuGet dependencies should be automatically downloaded
- Edit the yaml.default with the configuration for your environment
	- Pick an implementation of `ISwapIdentifiers` e.g. `Microservices.IdentifierMapper.Execution.IdentifierSwapper` and enter the full Type name into default.yaml `IdentifierMapperOptions.SwapperType`.
	- Specify the mapping table database details*
		- MappingConnectionString, the connection string to use to connect to the mapping server 
		- MappingDatabaseType, either MicrosoftSQLServer or MYSQLServer
		- MappingTableName, the table on the mapping server that contains the identifier mapping table* (identifiable=>anonymous)
		- SwapColumnName, the column in the `MappingTableName` that will contain the expected (identifiable) input values to replace.
		- ReplacementColumnName, the column in the `MappingTableName` that contains the replacement values.

- Decide if you want to [use a Redis](#Redis) cache.

*The table/connection string details are at the disposal of the `ISwapIdentifiers` chosen.  Some might ignore them completely or might manually create the mapping table themselves (e.g. `ForGuidIdentifierSwapper`)

### 3. Exchange and Queue Settings
| Read/Write | Type | Config setting |
| ------------- | ------------- |------------- |
| Read | DicomFileMessage | `IdentifierMapperOptions.QueueName` |
| Write | DicomFileMessage | `IdentifierMapperOptions.AnonImagesProducerOptions` |

Expects to receive control messages informing it to 'refresh' it's mapping table.

### 4. Config
| YAML Section  | Purpose |
| ------------- | ------------- |
| RabbitOptions | Describes the location of the rabbit server for sending messages to. |
| IdentifierMapperOptions | Describes location of the mapping table (e.g. CHI=>EUPI), the Type responsible for mapping and all queue/exchange settings.|

| Command Line Options | Purpose |
| ------------- | ------------- |
|CliOptions | Allows overriding of which yaml file is loaded. |

### 5. Redis

If you are using an `ISwapper` implementation that consults a large mapping database e.g. 10 million then you may benefit from using a Redis caching database.  Install Redis and set the `RedisHost` option in the config file e.g.

```yaml
IdentifierMapperOptions:
    SwapperType: 'Microservices.IdentifierMapper.Execution.Swappers.TableLookupSwapper'
    RedisHost: localhost
```

All lookup results will be cached in the Redis server (both successful lookups and misses).

If you update your lookup table you will have to manually flush the Redis server (if desired).

### 6. Expectations

All identifier allocation is handled by the chosen `ISwapIdentifiers`.  The only field considered is `DicomTag.PatientID` which should be the only patient identifier field loaded by any downstream processes.

You should ensure the integrity of your mapping table yourself through an appropriate schema (e.g. not null, primary keys) to prevent 1->Many mappings and blank mappings etc at data load time.  It is not the job of the `ISwapIdentifiers` to handle this.

#### Data Failure States
- Dicom file doesn't contain a `DicomTag.PatientID`
	- Nacked `Test_BlankPatientIdentifier`
- PatientID not found by `ISwapIdentifiers`
	- Nacked `Test_NoMatchingIdentifierFound`

#### Environmental Failure States
 - Operation on loss of RabbitMQ connection:
	- No special logic
- Operation on loss of access to mapping table:
	- Any Exception thrown by the `ISwapIdentifiers` will not be caught triggering a Fatal on the `Consumer`.

### 6. Class Diagram
![Class Diagram](./Images/ClassDiagram.png)

# DicomReprocessor

Primary Author: [Ruairidh MacLeod](https://github.com/rkm)


## Contents

 1. [Overview](#1-overview)
 2. [Setup / Installation](#2-setup--installation)
 3. [Config](#4-config)
 4. [Queue Settings](#3-queue-settings)
 5. [Expectations](#5-expectations)
 6. [Class Diagram](#6-class-diagram)


### 1. Overview

Pulls documents from a MongoDB collection and republishes them to RabbitMQ as DicomFileMessages.

In future this will have to option to run in 'TagPromotion' mode, where only selected tags are republished.


### 2. Setup / Installation

- Clone the project and build. Any NuGet dependencies should be automatically downloaded
- Setup a yaml file with the configuration for your environment
- Run the normal load pipeline so that your MongoDb has some records in it
- Run `dotnet DicomReprocessor.dll -y <config> -c <collection> [options]`


### 3. Config

Requires the following fields from the default YAML config:

- `RabbitOptions`
- `MongoDatabases.DicomStoreOptions`
- `DicomReprocessorOptions`


### 4. Exchanges and Queues

Consumes messages from:

- `None`

Writes messages to:

- `DicomReprocessorOptions.ReprocessingProducerOptions`


### 5. Expectations

Errors are [logged as normal for a MicroserviceHost](../../common/Smi.Common/README.md#logging)

### 6. Class Diagram

![Class Diagram](./Images/ClassDiagram.png)
﻿# IsIdentifiable

Primary Author: [Thomas](https://github.com/tznind)

## Contents
 1. [Overview](#overview)
 1. [Setup](#setup)
 1. [Invocation](#invocation)
 1. [Rules](#rules) 
    1. [Basic Rules](#basic-rules) 
    2. [Socket Rules](#socket-rules) 
    3. [Consensus Rules](#consensus-rules) 
    4. [White List Rules](#white-list-rules) 
 1. [Exchange and Queue Settings](#exchange-and-queue-settings)
 1. [Expectations](#expectations)
 1. [Class Diagram](#class-diagram)

## Overview
This service evaluates 'data' for personally identifiable values (e.g. names).  It can source data from a veriety of places (e.g. databases, file system).

## Setup

To run IsIdentifiable you must first build the microservice then download the required data models for NER and OCR.
Rules must be placed into a suitable directory. Data and rules are supplied in the SmiServices data directory.

### Downloads

The following downloads are required to run the software:

| File     | Destination |  Windows Script |  Linux Script  |
|----------|-------------|-------- |------|
| [Tesseract Data files (pixel OCR models)](https://github.com/tesseract-ocr/tessdata/raw/master/eng.traineddata) | `./data/tessdata` |  [download.ps1](../../../data/tessdata/download.ps1)|  [download.sh](../../../data/tessdata/download.sh)|
| [Stanford NER Classifiers](http://nlp.stanford.edu/software/stanford-ner-2016-10-31.zip)*    |  `./data/stanford-ner`     | [download.ps1](../../../data/stanford-ner/download.ps1)  | [download.sh](../../../data/stanford-ner/download.sh) |

_*Required for NERDaemon_
 
## Invocation

IsIdentifiable can be run in one of several modes:

 * As a microservice host to process DICOM files named in RabbitMQ messages
 * Interactively to process a DICOM file or a directory of DICOM files
 * Interactively to process a every row of every column in a database table

To run as a service use `dotnet IsIdentifiable.dll service -y default.yaml [options]`

To run on files use `dotnet IsIdentifiable.dll dir -d /path/dir [--pattern *.dcm] [options]`

To run on a database use `dotnet IsIdentifiable.dll db -d database -t table -p dbtype [options]`

### Examples

```bash
dotnet publish

cd ./src/microservices/Microservices.IsIdentifiable/bin/AnyCPU/Debug/netcoreapp2.2/

#Generic help (lists modes)
dotnet IsIdentifiable.dll --help

#Specific help (for a given mode e.g. 'db')
dotnet ./IsIdentifiable.dll db --help
```

An example command (evaluate all the images in `C:\MassiveImageArchive`) would be as follows:

```
dotnet IsIdentifiable.dll dir -d C:/MassiveImageArchive --storereport
```

The outputs of this (on some anonymised data):

```
Resource,ResourcePrimaryKey,ProblemField,ProblemValue,PartWords,PartClassifications,PartOffsets
C:\MassiveImageArchive\DOI\000001.dcm,1.3.6.1.4.1.9590.100.1.2.64408251011211630124074907290278463475,"(0008,0005)",ISO_IR 100,ISO,Organization,0
C:\MassiveImageArchive\DOI\000001.dcm,1.3.6.1.4.1.9590.100.1.2.64408251011211630124074907290278463475,"(0018,1016)",MathWorks,MathWorks,Organization,0
C:\MassiveImageArchive\DOI\000001.dcm,1.3.6.1.4.1.9590.100.1.2.64408251011211630124074907290278463475,"(0018,1018)",MATLAB,MATLAB,Person,0
C:\MassiveImageArchive\DOI\000001.dcm,1.3.6.1.4.1.9590.100.1.2.64408251011211630124074907290278463475,"(0020,0010)",DDSM,DDSM,Person,0
C:\MassiveImageArchive\DOI\000002.dcm,1.3.6.1.4.1.9590.100.1.2.423893162212842428532864042250901777433,"(0008,0005)",ISO_IR 100,ISO,Organization,0
C:\MassiveImageArchive\DOI\000002.dcm,1.3.6.1.4.1.9590.100.1.2.423893162212842428532864042250901777433,"(0018,1018)",MATLAB,MATLAB,Person,0
C:\MassiveImageArchive\DOI\000002.dcm,1.3.6.1.4.1.9590.100.1.2.423893162212842428532864042250901777433,"(0020,0010)",DDSM,DDSM,Person,0
C:\MassiveImageArchive\DOI\000003.dcm,1.3.6.1.4.1.9590.100.1.2.84709658512632788123980174250729731712,"(0008,0005)",ISO_IR 100,ISO,Organization,0
C:\MassiveImageArchive\DOI\000003.dcm,1.3.6.1.4.1.9590.100.1.2.84709658512632788123980174250729731712,"(0018,1016)",MathWorks,MathWorks,Organization,0
C:\MassiveImageArchive\DOI\000003.dcm,1.3.6.1.4.1.9590.100.1.2.84709658512632788123980174250729731712,"(0018,1018)",MATLAB,MATLAB,Person,0
C:\MassiveImageArchive\DOI\000003.dcm,1.3.6.1.4.1.9590.100.1.2.84709658512632788123980174250729731712,"(0020,0010)",DDSM,DDSM,Person,0
C:\MassiveImageArchive\DOI\Calc-Test_P_00038_LEFT_CC\1.3.6.1.4.1.9590.100.1.2.85935434310203356712688695661986996009\1.3.6.1.4.1.9590.100.1.2.374115997511889073021386151921807063992\000000.dcm,1.3.6.1.4.1.9590.100.1.2.289923739312470966435676008311959891294,"(0008,0005)",ISO_IR 100,ISO,Organization,0
C:\MassiveImageArchive\DOI\Calc-Test_P_00038_LEFT_CC\1.3.6.1.4.1.9590.100.1.2.85935434310203356712688695661986996009\1.3.6.1.4.1.9590.100.1.2.374115997511889073021386151921807063992\000000.dcm,1.3.6.1.4.1.9590.100.1.2.289923739312470966435676008311959891294,"(0018,1016)",MathWorks,MathWorks,Organization,0
C:\MassiveImageArchive\DOI\Calc-Test_P_00038_LEFT_CC\1.3.6.1.4.1.9590.100.1.2.85935434310203356712688695661986996009\1.3.6.1.4.1.9590.100.1.2.374115997511889073021386151921807063992\000000.dcm,1.3.6.1.4.1.9590.100.1.2.289923739312470966435676008311959891294,"(0018,1018)",MATLAB,MATLAB,Person,0
C:\MassiveImageArchive\DOI\Calc-Test_P_00038_LEFT_CC\1.3.6.1.4.1.9590.100.1.2.85935434310203356712688695661986996009\1.3.6.1.4.1.9590.100.1.2.374115997511889073021386151921807063992\000000.dcm,1.3.6.1.4.1.9590.100.1.2.289923739312470966435676008311959891294,"(0020,0010)",DDSM,DDSM,Person,0
[...]
```

You can run pixel data (OCR) by passing the `--tessdirectory` flag:

```
dotnet IsIdentifiable.dll dir -d C:\MassiveImageArchive --storereport --tessdirectory E:/SmiServices/data/tessdata/
```

The directory must be named `tessdata` and contain a file named `eng.traineddata`

## Rules

Rules can be used to customise the way failures are handled.
A failure is a fragment of text (or image) which contains identifiable data.
It can either be ignored (because it is a false positive) or reported.

Some rules come out of the box (e.g. CHI/Postcode/Date) but for the rest you must configure rules in a rules.yaml file.
There are three classes of rule: BasicRules, SocketRules and WhiteListRules. See below for more details of each.
They are applied in that order, so if a value is Ignored in a Basic rule it will not be passed to any further rules.
If a value fails in a SocketRule (eg. the NER daemon labels it as a Person), then a subsequent WhiteList rule can Ignore it.
Not all Ignore rules go into the WhiteListRules section; this is intended only for white-listing fragments which NERd has incorrectly reported as failures.

Rules can be read from one or more yaml files. Each file can have zero or one set of BasicRules, plus zero or one set of WhiteListRules.
All of the BasicRules from all of the files will be merged to form a single set of BasicRules; similarly for WhiteListRules.

When running in service mode (as a microservice host) the rules are read from all `*.yaml` files found in the `IsIdentifiableRules` directory inside the data directory. The data directory path is defined in the program yaml file (which was specified with -y) using the key `IsIdentifiableOptions|DataDirectory`.

When running in file or database mode the yaml file option -y is not used so the command line option `--rulesdirectory` needs to be specified. This should be the path to the IsIdentifiableRules directory.

### Basic Rules

These can either result in a value being Reported or Ignored (i.e. not passed to any downstream classifiers).  Rules can apply to all columns (e.g. Ignore the Modality column) or only those values that match a Regex. The regex is specified using `IfPattern`.  Regex case sensitivity is determined by the `CaseSensitive` property (defaults to false).

```yaml
BasicRules: 
  # Report any values which contain 2 digits as a PrivateIdentifier
  - IfPattern: "[0-9][0-9]"
    Action: Report
    As: PrivateIdentifier

  # Do not run any classifiers on the Modality column
  - Action: Ignore
    IfColumn: Modality
```

### Socket Rules

You can outsource the classification to separate application(s) (e.g. NERDaemon) by adding `Socket Rules`

```yaml
SocketRules:   
  - Host: 127.0.123.123
    Port: 1234
```

The TCP protocol starts with IsIdentifiable sending the word for classification i.e.

```
Sender: word or sentence\0
```

The service is expected to respond with 0 or more classifications of bits in the word that are problematic.  These take the format:

```
Responder: Classification\0Offset\0Offending Word(s)\0
```

Once the responder has decided there are no more offending sections (or there were none to begin with) it sends a double null terminator.  This indicates that the original word or sentence has been fully processed and the Sender can send the next value requiring validation.

```
Responder: \0\0
```

### Consensus Rules

If you have two or more rules that you want to cooperate when determining whether data is identifiable or not e.g. 2 NLP Name Entity Recognizers you can use a ConsensusRule.  These rules require all subrules to agree on whether to Report or Ignore a given cell value.  If there is disagreement then the rule is not applied (`RuleAction.None` i.e. take no action).

You can configure a consensus rule using the following yaml:
```yaml
ConsensusRules:
    - Rules:
      - !SocketRule
          Host: 127.0.123.123
          Port: 1234
      - !SocketRule
          Host: 127.0.123.123
          Port: 567
```

Consensus rules are specifically designed for intersecting two or more over matching rules e.g. NLP classifications.  If only one rule flags something as identifiable it will be ignored (both must agree).

### White List Rules

White list rules are a last chance filter on the final output of all other rules.  They allow discarding rules based on the whole string or the specific failing part.

The Action for a White List rule must be Ignore because it is intended to allow values previously reported to be ignored as false positives.

All of the constraints must match in order for the rule to Ignore the value.

As soon as a value matches a white list rule no further white list rules are needed.
Unlike a BasicRule whose Pattern matches the full value of a field (column or DICOM tag) the whitelist rule has two Patterns, IfPattern which has the same behaviour and IfPartPattern which matches only the substring that failed. This feature allows context to be specified, see the second example below.
A whitelist rule can also match the failure classification (`PrivateIdentifier`, `Location`, `Person`, `Organization`, `Money`, `Percent`, `Date`, `Time`, `PixelText`, `Postcode`).
For example, if SIEMENS has been reported as a Person found in the the Manufacturer column,

```yaml
WhiteListRules:
 - Action: Ignore
   As: Person
   IfColumn: Manufacturer
   IfPartPattern: ^SIEMENS$
```
For example, what seems like a name Brian can be ignored if it occurs in the exact phrase "MR Brian And Skull" using:

```yaml
IfPartPattern: ^Brian$
IfPattern: MR Brian And Skull
```

Note that there is no need to specify ^ and $ in IfPattern as other text before or after it will not change the meaning.

## Exchange and Queue Settings

In order to run as a microservice you should call it with the `service` option

| Read/Write | Type | Config setting |
| ------------- | ------------- |------------- |
| Read | ExtractFileMessage | IsIdentifiableOptions.QueueName |
| Write | IsIdentifiableMessage | IsIdentifiableOptions.IsIdentifiableProducerOptions.ExchangeName |

## Config

| YAML Section  | Purpose |
| ------------- | ------------- |
| RabbitOptions | Describes the location of the rabbit server for sending messages to |
| IsIdentifiableOptions | Describes what `IClassifier` to run and where the classifier models are stored. The key `DataDirectory` specifies the path to the data directory. The key `ClassifierType` specifies which classifier to run, typically `Microservices.IsIdentifiable.Service.TesseractStanfordDicomFileClassifier` |

## Expectations

> TODO: 

### Data Failure States

> TODO: 

### Environmental Failure States
 
> TODO: 

## Class Diagram
![Class Diagram](./IsIdentifiable.png)
# Name Entity Recognition Daemon

Primary Author: [James A Sutherland](https://github.com/jas88)
Python author: Andrew Brooks

## Contents

1. [Background](#background)
1. [Setup](#setup)
1. [Running](#running)


## Background

This standalone process is designed to classify text strings sent by the [IsIdentifiable](../Microservices.IsIdentifiable/README.md#socket-rules) application. It accepts TCP connections on localhost port 1881, returning classification results as expected by the IsIdentifiable microservice.

There are two implementations, one in Java which uses the [Stanford CoreNLP](https://stanfordnlp.github.io/CoreNLP/) library and one in Python which uses the [SpaCy](https://spacy.io/) library.

The [Stanford CoreNLP NER algorithm](https://stanfordnlp.github.io/CoreNLP/ner.html#description) recognizes named (PERSON, LOCATION, ORGANIZATION, MISC), numerical (MONEY, NUMBER, ORDINAL, PERCENT), and temporal (DATE, TIME, DURATION, SET) entities (12 classes). Adding the regexner annotator and using the supplied RegexNER pattern files adds support for the fine-grained and additional entity classes but the additional annotator is not currently used here. 

The Python version recognises a different set of entities, and with different language models can recognise drugs, diseases, cells, etc. although these are not required for named entities. It tries to map its entity type names to match those returned by CoreNLP so that both can be used in a ConsensusRule in IsIdentifiable.

## Setup

No setup is required for the Java version, just run the jar file as documented below. The "<&- &" will cause it to disconnect from the terminal once initialised and run as a daemon. (For development use, you can also skip that and terminate it with ctrl-C on the console when finished.)

The Python program now requires a YAML configuration file to determine the log file location. It also requires that the SpaCy and optionally the SciSpaCy packages have been installed, if not globally then into a virtual environment. The same environment must also have the required SpaCy language model installed. Note that SpaCy version 2 (eg. 2.2.1) and SciSpacy version 0.2.4 must be used (as of Feb 2021) because SpaCy v3 uses a new architecture and SciSpaCy has not caught up yet (at least, not for NER).

## Running

`java -jar nerd.jar <&- &`

To shutdown again:

`fuser -k -TERM -n tcp 1881`

The Python version is:

`ner_daemon_spacy.py --yaml default.yaml --model en_core_web_md --port 1882`

The yaml must have `LoggingOptions | LogsRoot`.

The port should be different from NERd if both are to run in parallel.

Models are: `en_core_web_md` (the default), `en_core_sci_md` (scispacy, but warning: does not return entity types so no use for finding PII), `en_ner_bionlp13cg_md` and `en_ner_bc5cdr_md` (scispacy models for drugs and diseases). Additional SpaCy models can be used which are much larger, eg. `en_core_web_lg`.
# DicomTagReader

Primary Author: [Ruairidh MacLeod](https://github.com/rkm)

## Contents
 1. [Overview](#1-overview)
 2. [Setup / Installation](#2-setup--installation)
 3. [Exchange and Queue Settings](#3-exchange-and-queue-settings)
 4. [Config](#4-config)
 5. [Expectations](#5-expectations)
 6. [Class Diagram](#6-class-diagram)

### 1. Overview
Opens dicom files found in `AccessionDirectoryMessage` directories and converts to JSON as a `DicomFileMessage`. Also creates a summary record of the whole series as a `SeriesMessage`.

### 2. Setup / Installation
 - Clone the project and build. Any NuGet dependencies should be automatically downloaded
 - Edit the default.yaml with the configuration for your environment

### 3. Exchange and Queue Settings
| Read/Write | Type | Config setting |
| ------------- | ------------- |------------- |
| Read | AccessionDirectoryMessage | `DicomTagReaderOptions.QueueName` |
| Write | DicomFileMessage | `DicomTagReaderOptions.ImageProducerOptions` |
| Write | SeriesMessage | `DicomTagReaderOptions.SeriesProducerOptions` |

### 4. Config
| YAML Section  | Purpose |
| ------------- | ------------- |
| RabbitOptions | Describes the location of the rabbit server for sending messages to. |
| FileSystemOptions | Describes the root location of all images this program will load, all `AccessionDirectoryMessage` processed by this microservice should be assumed to refer to images relative to this root. |
| DicomTagReaderOptions | Contains names of the series and image exchanges that serialized image tag data will be written to.|

| Command Line Options | Switch |  Purpose |
| ------------- | ------------- | ------------- |
|CliOptions | -y, --yaml-file| Allows overriding of which yaml file is loaded. |

### 5. Expectations
Errors are [logged as normal for a MicroserviceHost](../../common/Smi.Common/README.md#logging)

### 6. Class Diagram
![Class Diagram](./Images/ClassDiagram.png)

# Cohort Extractor

Primary Author: [Thomas](https://github.com/tznind)

## Contents
 1. [Overview](#1-overview)
 2. [Setup / Installation](#2-setup--installation)
 3. [Exchange and Queue Settings](#3-exchange-and-queue-settings)
 4. [Config](#4-config)
    - [Fulfiller]
    - [Rejector]
 5. [Expectations](#5-expectations)
 6. [Class Diagram](#6-class-diagram)

### 1. Overview
This service services `ExtractionRequestMessage` which is a request to extract a given set images identified by a key tag (e.g. `SeriesInstanceUID`) collection (e.g. 5000 SeriesInstanceUID values).  It is the job of the Cohort Extractor to identify the images which correspond to the specified key values requested (e.g. the `SeriesInstanceUID`) and generate output messages to downstream processes responsible for anonymising the images.

There can be multiple datasets in which matching images should be sourced e.g. MR / CT which could even reside on different servers.  Datasets are identified and distinguished from one another through RDMP `ICatalogue` which exists already as part of the data load process (See `DicomRelationalMapper`). 

### 2. Setup / Installation
- Clone the project and build. Any NuGet dependencies should be automatically downloaded
- Edit the yaml.default with the configuration for your environment
- Pick an implementation of `IAuditExtractions` e.g. `Microservices.CohortExtractor.Audit.NullAuditExtractions` and enter the full Type name into default.yaml `AuditorType`
- Pick an implementation of `IExtractionRequestFulfiller` e.g. `Microservices.CohortExtractor.Execution.RequestFulfillers.FromCataloguesExtractionRequestFulfiller` and enter the full Type  name into default.yaml `RequestFulfillerType`.
- Specify the mapping RDMP catalogue database
- Optionally specify a list of Catalogue IDs in CataloguesToExtractFrom (or set it to * to use any).  Depending on your `IExtractionRequestFulfiller` this value might be ignored.

### 3. Exchange and Queue Settings

| Read/Write | Type | Config setting |
| ------------- | ------------- |------------- |
| Read | ExtractionRequestMessage | CohortExtractorOptions.QueueName |
| Write | ExtractFileMessage| CohortExtractorOptions.ExtractFilesProducerOptions |
| Write | ExtractFileCollectionInfoMessage| CohortExtractorOptions.ExtractFilesInfoProducerOptions |


| Command Line Options | Purpose |
| ------------- | ------------- |
|CliOptions | Allows overriding of which yaml file is loaded. |

### 4. Config
| YAML Section  | Purpose |
| ------------- | ------------- |
| RabbitOptions | Describes the location of the rabbit server for sending messages to |
| RDMPOptions | Describes the location of the Microsoft Sql Server RDMP platform databases which keep track of load configurations, available datasets to extract images from (tables) etc |
| CohortExtractorOptions | Which Catalogues to extract, which classes to instantiate to do the extraction |

| Command Line Options | Purpose |
| ------------- | ------------- |
|CliOptions | Allows overriding of which yaml file is loaded. |

#### Fulfiller

The set of images that __could__ be extracted is controlled by the `IExtractionRequestFulfiller`.  

The current recommended implementation is [FromCataloguesExtractionRequestFulfiller].  This fulfiller will look up one or more tables or multi table joins (Catalogues) and search for the provided extraction key (e.g. SeriesInstanceUID = x)

The matched records are what will be reported on e.g. "for x UIDs we found y available images".  From this result set a subset will be rejected (because you have made row level decisions not to extract particular images).  This is handled by the [Rejector]

Configure the fulfiller in your options yaml:

```yaml
OnlyCatalogues: 1,2,3
RequestFulfillerType: Microservices.CohortExtractor.Execution.RequestFulfillers.FromCataloguesExtractionRequestFulfiller 
```

#### Rejector

Records matched by the [Fulfiller] are passed to the `IRejector` (if any is configured).  This class can make last minute decisions on a row by row level to either extract or forbid (with a specific provided reason) the extraction of an image.

The currently recommended implementation is [DynamicRejector]. To use the dynamic rejector edit your options yaml as follows:

```yaml
RejectorType: Microservices.CohortExtractor.Execution.RequestFulfillers.Dynamic.DynamicRejector
```

Using the [DynamicRejector] also requires you to configure a file [DynamicRules.txt] in the execution directory of your binary.  An example is provided (see [DynamicRules.txt]).

Rules are written in C# and can only index fields that appear in the records returned by the [Fulfiller].

### 5. Expectations

All matching of request criteria is handled by `IExtractionRequestFulfiller`.

All audit is handled by `IAuditExtractions`.

The extraction destination is handled by `IProjectPathResolver`

#### Data Failure States

- No files matching a given tag X
	- ???
- No value in patient id substitution Y
	- ???
- Others? TODO


#### Environmental Failure States

 - Operation on loss of RabbitMQ connection:
	- No special logic
- Operation on loss of access to catalogues:
	- Any Exception thrown by the `ISwapIdentifiers` will not be caught triggering a Fatal on the `Consumer`.

	
### 6. Class Diagram
![Class Diagram](./Images/ClassDiagram.png)

[Rejector]: #rejector
[Fulfiller]: #fulfiller
[DynamicRules.txt]: ./DynamicRules.txt
[DynamicRejector]: ./Execution/RequestFulfillers/Dynamic/DynamicRejector.cs
[FromCataloguesExtractionRequestFulfiller]: ./Execution/RequestFulfillers/FromCataloguesExtractionRequestFulfiller.cs
﻿
# Update Values

Primary Author: [Thomas](https://github.com/tznind)

## Contents
 1. [Overview](#1-overview)
 2. [Setup / Installation](#2-setup--installation)
 3. [Exchange and Queue Settings](#3-exchange-and-queue-settings)
 4. [Config](#4-config)
 5. [Expectations](#5-expectations)
 6. [Class Diagram](#6-class-diagram)

### 1. Overview
This service services `UpdateValuesMessage` which is a request to update a concept e.g. PatientID in one or more tables.  Each message describes a single query (although multiple values can be updated at once).

### 2. Setup / Installation
- Clone the project and build. Any NuGet dependencies should be automatically downloaded
- Edit the yaml.default with the configuration for your environment
- Optionally specify a list of **TableInfo** IDs in TableInfosToUpdate (otherwise all currently configured TableInfo will be used)

### 3. Exchange and Queue Settings

| Read/Write | Type | Config setting |
| ------------- | ------------- |------------- |
| Read | UpdateValuesMessage | UpdateValuesOptions.QueueName |

| Command Line Options | Purpose |
| ------------- | ------------- |
|CliOptions | Allows overriding of which yaml file is loaded. |

### 4. Config
| YAML Section  | Purpose |
| ------------- | ------------- |
| RabbitOptions | Describes the location of the rabbit server to pulling messages from |
| RDMPOptions | Describes the location of the Microsoft Sql Server RDMP platform databases which keep track of load configurations, available datasets to extract images from (tables) etc |
| UpdateValuesOptions | Which tables to update, timeouts etc |

| Command Line Options | Purpose |
| ------------- | ------------- |
|CliOptions | Allows overriding of which yaml file is loaded. |

### 5. Expectations

Each message processed will result in one or more UPDATE SQL statements being run.  These may run on different servers and different DBMS (Oracle, MySql, Sql Server, Postgres).  This ensures a relatively consistent 'whole system' update of a given fact.

#### Data Failure States

- No tables exist that contain the updated field(s)
- The TableInfo in RDMP maps to a non existent table (e.g. if the database server has been shutdown or the table renamed/deleted)

#### Environmental Failure States

 - Operation on loss of RabbitMQ connection:
	- No special logic

	
### 6. Class Diagram
![Class Diagram](./Images/ClassDiagram.png)
# Azure Pipelines Build Notes

This directory contains the Azure Pipelines definitions and all other files used to test SmiServices.

## Pipelines

We currently have the following pipelines:
-   [Linux test & package](https://dev.azure.com/smiops/Public/_build?definitionId=3). Defined [here](/.azure-pipelines/linux.yml). Runs the C# and Java tests on Linux. If the build is for a tag, also packages both sets of services and uploads them to a GitHub release
-   [Windows test & package](https://dev.azure.com/smiops/Public/_build?definitionId=4). Defined [here](/.azure-pipelines/windows.yml). Same as above, but on Windows. Runs a reduced set of tests, since some services are not available (see below).
-   [Create GitHub Release](https://dev.azure.com/smiops/Public/_build?definitionId=5). Defined [here](/.azure-pipelines/create-gh-release.yml). Runs automatically for tags on master to create a GitHub release which the other pipelines will publish packages to

## Directory Contents

-   `docker-compose/` - docker-compose files and lockfiles (see below)
-   `jobs/` - job definitions
-   `scripts/` - scripts used by these pipelines
-   `steps/` - individual task step definitions
-   `*.yml` - The top-level pipeline definitions
-   `vars.json` - Variables for the pipelines

Note that yaml files which are used as templates have the extension `.tmpl.yml`.

## Variables

Variables for the pipelines are loaded from the `vars.json` file at the start of each run. These are often combined with [pre-defined variables](https://docs.microsoft.com/en-us/azure/devops/pipelines/build/variables).

## Services

Service | Linux Provider (`ubuntu-18.04`) | Windows Provider (`windows-2019`)
 ------ | -------------- | ----------------
RabbitMQ | Docker | -
MongoDB | Docker | pre-installed
MsSQL | Docker | [SqlLocalDB](https://docs.microsoft.com/en-us/sql/database-engine/configure-windows/sql-server-express-localdb?view=sql-server-ver15)
MariaDB | Docker | -
Redis | Docker | Unavailable

## Docker

We use Docker containers for running the external services wherever possible. This currently means "when on Linux", since the Azure Windows OS currently doesn't support WSL2 / running Linux containers. Other options could be investigated for Windows, e.g. see a similar discussion [here](https://github.com/opensafely/job-runner/issues/76).

The docker-compose files reference the `latest` tag for each image. During the pipeline run however, these are replaced with specific image digest versions using the [docker-lock](https://github.com/safe-waters/docker-lock) tool. This is so the docker-compose files can be used as cache keys to enable repeatable builds.

## Caches

The [Cache](https://docs.microsoft.com/en-us/azure/devops/pipelines/release/caching?view=azure-devops) task is used in multiple cases to speed-up the build by re-using previously built/downloaded assets. These are restored based on the `key` value, which can contain references to files which are hashed to generate the final key.

## Notes

-   `set -x` in bash tasks may interfere with the `##vso` syntax, as AP will interpret any stdout containing that string as a command. This can cause problems with variables being set twice and having the wrong quote escaping. See:
    -   https://github.com/Microsoft/azure-pipelines-tasks/issues/10165
    -   https://github.com/microsoft/azure-pipelines-tasks/issues/10331
    -   https://developercommunity.visualstudio.com/content/problem/375679/pipeline-variable-incorrectly-inserts-single-quote.html

-   Dealing with variables containing path separators across multiple platforms (i.e. Linux and Windows) can be tricky. Variables pre-defined by AP will have the correct separators for the current OS. When used in `Bash` tasks, it's safe enough for either OS to just use Linux (`/`) separators when combining paths, so long as they are always quoted. However, if the variable is also used as part of a cache `path`, it may need to be manually re-created with Windows-style (`\`) separators before use.

-   The `##vso[task.setvariable variable=foo]bar` syntax requires that the variable name be lowercase. This is not documented explicitly but results in `command not found` errors if uppercase names are used.
# Utils

## CI build scripts

-   buildArtefacts.py. Creates packages for all apps in this repo
-   createReleaseChangelog.py. Creates the changelog text for upload to GitHub's Releases page
-   dotnet-build.bash. Builds the dotnet services in Release mode
-   install-ctp.bash. Installs the required libs for CTP
-   run-java-tests.bash. Runs the Java tests
-   runDotnetTests.py. Runs the dotnet tests in Release mode and generates a coverage report
-   updateChangelog.py. Updates the CHANGELOG from the news files

## Other utils

-   RabbitMqTidyQueues. Automatically deletes any Control queues which have no consumers (those that have not been deleted by their host due to a crash).
-   RabbitMQ Dump Queue. Allows dumping of all messages in a RabbitMQ queue to JSON files. It can be downloaded from [here](https://github.com/dubek/rabbitmq-dump-queue).

# RabbitMQ Tidy Queues Utility

Language: `Go`

Deletes any control queues that have not been properly cleaned by services which have crashed.

## Building

Download and install [Go](https://golang.org/dl/) and its associated build tools.

### Build for your local system

```bash
> go build
```

### Build on Windows for Linux

Run the [WindowsBuildForLinux.bat](./WindowsBuildForLinux.bat) script. This builds an ELF binary for `amd64-linux`.

In general, you need to appropriately set `GOARCH` and `GOOS` to build for another system, and then unset them afterwards to ensure they get returned back to their default values for your system.

## Usage

Run the build output on the same node as the RabbitMQ server.

# TODO Move these from the old RDMP tests project (if they exist)
# TODO Separate the integration tests to here
## Proposed Changes

Summarise your proposed changes here, including any notes for reviewers.

## Types of changes

What types of changes does your code introduce? Tick all that apply.

- [ ] Bugfix (non-breaking change which fixes an issue)
- [ ] New Feature (non-breaking change which adds functionality)
- [ ] Breaking Change (fix or feature that would cause existing functionality to not work as expected)
- [ ] Documentation-Only Update (if none of the other choices apply)
  - In this case, ensure that the message of the head commit from the source branch is prefixed with `[skip ci]`

## Checklist

By opening this PR, I confirm that I have:

- [ ] Reviewed the [contributing](https://github.com/SMI/SmiServices/blob/master/CONTRIBUTING.md) guidelines for this repository
- [ ] Ensured that the PR branch is in sync with the target branch (i.e. it is automatically merge-able)
- [ ] Updated any relevant API documentation
- [ ] Created or updated any tests if relevant
- [ ] Created a [news](https://github.com/SMI/SmiServices/blob/master/news/README.md) file
    -   NOTE: This ***must*** include any changes to any of the following files: default.yaml, any of the RabbitMQ server configurations, GlobalOptions.cs
- [ ] Listed myself in the [CONTRIBUTORS](https://github.com/SMI/SmiServices/blob/master/CONTRIBUTORS.md) file 🚀
- [ ] Requested a review by one of the repository maintainers

## Issues

If relevant, tag any issues that are *expected* to be resolved with this PR. E.g.:

- Closes #\<issue-number>
- ...

# TODO Move all dicom files and other common resources here
# SMI Services Configuration

TODO - Add docs
# RabbitMQ Configuration

## Reset the server

The following will **completely** reset the state of your RabbitMQ instance, including removing all users, vhosts, exchanges, queues, vhosts, and messages.

Find your installation folder (default appears to be `C:\Program Files\RabbitMQ Server\rabbitmq_server-3.7.3\sbin` on Windows). From there, run the following:

```
rabbitmqctl.bat stop_app

rabbitmqctl.bat reset

rabbitmqctl.bat start_app
```

Can also just run the `ResetRabbitMQ.ps1` script included in this directory.

For Linux, remove the `.bat` extension from the commands above, and run as sudo.

## Uploading broker defnintions

From Linux, the exchange and queue definitions can be uploaded using curl:

```bash
$ curl \
-u guest:guest \
-XPOST -H"content-type:application/json" \
-d"@/full/path/to/defaultExtractConfig.json" \
http://0.0.0.0:15672/api/definitions
```

## Deleting a vhost

```bash
> curl \
    -u guest:guest \
    -XDELETE \
    http://0.0.0.0:15672/api/vhosts/<name>
```

## Filter the default exchanges

RabbitMQ has predefined exchanges which can't be removed from the management UI. To filter these out, tick the `regex` checkbox next to the search box on the `Exchanges` tab, then use this regex:

```regex
^(?!amq\.)(.+)$
```

## Filter the control queues

The microservice control queues can be filtered out from the default display with the following regex:

```regex
^(?!Control\.)(.+)$
```
# Notes for developing SmiServices

## Generating local HTML coverage reports

Following the process described [here](https://docs.microsoft.com/en-us/dotnet/core/testing/unit-testing-code-coverage?tabs=linux#generate-reports)

```console
$ dotnet tool update -g dotnet-reportgenerator-globaltool
...
$ utils/runDotnetTests.py
...

$ reportgenerator \
    -reporttypes:html \
    -targetdir:htmlcov \
    -reports:coverage/coverage.cobertura.xml
...
2021-06-01T17:27:04: Writing report file 'htmlcov/index.html'
```

This can then be hosted on a local webpage using Python:

```console
$ python3 -m http.server --directory htmlcov/
Serving HTTP on 0.0.0.0 port 8000 (http://0.0.0.0:8000/) ...
```
$Name

$Description

| Name | Description |
|------|-------------|
$foreach CatalogueItem
| $Name | $Description |
$end
# Microservice Control Queues

This describes how the services can be controlled via RabbitMQ messages.

### Contents

- [Commands](#commands)
- [Sending a message](#sending-a-message)
- [Implementing a new control command handler](#implementing-a-new-control-command-handler)
- [Control Queues and Cleanup](#control-queues-and-cleanup)

## Commands

Commands are sent by publishing a message to the ControlExchange (specified in your config by `RabbitOptions.RabbitMqControlExchangeName`) with a specific routing key. This allows you to easily send them from the RabbitMQ web management page, or via a CLI.

RabbitMQ message routing keys are used to control which services receive the message. The current format for routing keys is `smi.control.<who>.<what>`. Where `<who>` is the name of the service, and `<what>` is some defined action. Note that all keys must be specified in lowercase. The currently defined actions are:

### General - any service

- `stop` - Stops the service
- `ping` - Logs a `pong` message. Useful for debugging

### DicomReprocessor

- `set-sleep-time-ms` - Sets the sleep time between batches. This also requires the new value to be set in the message body

### IdentifierMapper

- `refresh` - Refreshes any caches in use

### CohortPackager

- `processjobs` - Checks if any in progress jobs are complete

## Sending a message

Messages can be sent either via the Web UI or via a CLI (see below for details). In either case, the following applies:

- The `<who>` field must exactly match the name of the microservice process (e.g. `identifiermapper`)
- All routing keys should be lowercase
- `all` can be used as the `<who>` keyword to control all services
- A specific service can be messaged by including its `PID` at the end of the routing key. This is currently the only way to control a specific service instance rather than all services of a certain type

Examples of some routing keys:

```text
smi.control.all.stop # Stop all services
smi.control.dicomtagreader.stop # Stop all DicomTagReader services
smi.control.identifiermapper.refresh1234 # Refresh the IdentifierMapper service with PID `1234`
```

Note that some services may take some time to finish their current operation and exit after recieveing a `shutdown` command.


### Via the Web UI

On your RabbitMQ Management interface (`http://<rabbit host>:15672`), click `Exchanges` then `Control Exchange`. Expand the `Publish message` box then enter the message info. Any required content should be entered into the `Payload` box in plain text. Example:

![test](Images/control-queue-publish.PNG)

### Via the CLI

`TODO`

## Implementing A New Control Command Handler

Implement a class which contains a method with the following signature:

```c#
void MyControlHandler(string action, string message)
```

Then, instantiate your class and register its event in your host (must be a subclass of `MicroserviceHost`):

```c#
var controlClass = new MyControlClass(...);
AddControlHandler(controlClass.MyControlHandler);
```

That's it! Now you will be passed the full routing key for any control message addressed to your specific microservice type (i.e. where the `<who>` part of the routing key matches your microservice name), and any message content.

## Control Queues and Cleanup

The actual implementation of the control queues works as follows:

- When each service starts up, it creates a new queue named with its service name and process ID
- It then binds this queue to the global `ControlExchange`. Two bindings are created:
  - `smi.control.all.*`: Matches any "send to all" routing keys
  - `smi.control.<process_name>.*`: Matches "all services of my type" routing keys
- On shutdown (when the RMQ connection is closed), the control queue should be automatically deleted by the server

The creation of the control queue is performed during a single ad-hoc connection, and is not part of the standard Consumer process (for _reasons_). One consequence of this is that if a microservice crashes _after_ the control queue is created, but _before_ the actual subscription to the queue is started (i.e. at some point during startup before RabbitMQAdapter.StartConsumer is called), then the control queue may not be automatically deleted. This isn't really an issue other than causing visual clutter on the RabbitMQ management interface. These dangling queues can be manually deleted with the [TidyQueues](../utils/RabbitMqTidyQueues) utility tool.

# SMI Services Documentation

This is the documentation for the SMIServices platform. It should (hopefully) contain enough information to run your own instance of the service.

The platform is currently deployed in the National Safe Haven, so some documentation may specifically refer to that environment. The software should be deployable in any environment though, so please open an [issue](https://github.com/SMI/SmiServices/issues) if anything isn't clear.


### Contents

- [Controlling the services](#controlling-the-services)
- [TODO](#todo)


## Controlling the services

The services can be controlled by sending messages to the RabbitMQ control exchange with specific routing keys. See the [main doc](control-queues.md)


## TODO

- Figure out what documentation can be (safely) imported from the old private repo
# Extraction

Describes the build-install-test procedure, not the deployment into production.

Main docs:
https://github.com/SMI/SmiServices/tree/master/src/microservices/com.smi.microservices.ctpanonymiser

https://github.com/SMI/SmiServices/tree/release/1.2.0#image-extraction-microservices

See also: the extraction-refactoring branch
https://github.com/SMI/SmiServices/tree/feature/extraction-refactoring/docs/extraction

Other docs:
https://uoe.sharepoint.com/sites/SMI/Shared%20Documents/Forms/AllItems.aspx
https://git.ecdf.ed.ac.uk/SMI/SmiServiceOps/blob/master/Planning/ExtractionFlags
https://github.com/HicServices/SMIPlugin/blob/master/Documentation/Images/ExtractionMicroservices.png

# Building

See elsewhere the documents for building the Java programs.

# Prerequisites for testing

A RabbitMQ instance is required - you can run a test version inside a Docker container:

```
sudo docker run -d --hostname my-rabbit --name some-rabbit-mgt -p 5671:5671 -p 5672:5672 -p 5673:5673 -p 15671:15671 -p 15672:15672 -p 25672:25672 rabbitmq:3-management
```

# ExtractImages

```console
# Example CSV file
cat > extractme.csv << _EOF
SeriesInstanceUID,foo
1.2.826.0.1.3680043.2.1125.1.78969117856457473538394301521877227,1
_EOF

Edit default.yaml (RabbitOptions and FileSystemOptions)
You could do this programmatically with
`yq_linux_amd64 write --inplace d FileSystemOptions.FileSystemRoot /tmp`
although the current version of yq loses comments and unnecessary quotes.

Login to rabbit (localhost:15672) and create exchanges:
TEST.RequestExchange
TEST.RequestInfoExchange
Add bindings from those exchanges to any queue (TEST.xxx)

ProjectNum=001
rmdir /tmp/${ProjectNum}/tmp  # program gives error if dir already exists

# Interactive - answer y to create messages, or run with --non-interactive
$ ./smi extract-images -y default.yaml -p "$ProjectNum" extractme.csv
```

Two messages are created:

```
{"KeyTag":"SeriesInstanceUID","ExtractionIdentifiers":["1.2.826.0.1.3680043.2.1125.1.78969117856457473538394301521877227"],"ExtractionJobIdentifier":"bb1cbed5-a666-4307-a781-5b83926eaa81","ProjectNumber":"001","ExtractionDirectory":"001/tmp","JobSubmittedAt":"2019-12-19T10:49Z"}
```
and
```
{"KeyTag":"SeriesInstanceUID","KeyValueCount":1,"ExtractionJobIdentifier":"bb1cbed5-a666-4307-a781-5b83926eaa81","ProjectNumber":"001","ExtractionDirectory":"001/tmp","JobSubmittedAt":"2019-12-19T10:49Z"}
```

# CohortExtractor

Requires MySQL instance? so not described (yet).

Creates messages containing fields:
DicomFilePath: Path to the original file
ExtractionDirectory: Extraction directory relative to the extract root
OutputPath: Output path for the anonymised file, relative to the extraction directory
See: ~/src/SmiServices/src/common/Smi.Common/Messages/Extraction/ExtractFileMessage.cs
Inherits ExtractMessage so:
Guid ExtractionJobIdentifier
string ProjectNumber
string ExtractionDirectory
DateTime JobSubmittedAt

# CTPanonymiser

`cd ~/src/SmiServices/src/microservices/com.smi.microservices.ctpanonymiser/target`

A whitelist is required, available from the old repo as:
https://raw.githubusercontent.com/HicServices/SMIPlugin/develop/Documentation/Anon/dicom-whitelist.script
https://raw.githubusercontent.com/HicServices/SMIPlugin/develop/Documentation/Anon/dicom-whitelist.script.new
(possibly identical content, apart from whitespace/newlines??)
or from the new repo in the directories:
```
SmiServices/src/applications/com.smi.applications.extractorcli/anonScript.txt
SmiServices/src/microservices/com.smi.microservices.ctpanonymiser/src/test/resources/dicom-anonymizer.script
```
Haven't yet determined which one is correct.

Edit `default.yaml` (RabbitOptions and FileSystemOptions)

Login to rabbit (http://localhost:15672/) and create exchanges:
TEST.ControlExchange
TEST.FatalLoggingExchange
TEST.FileStatusExchange
and queue: TEST.ExtractFileQueue
Check: do we need to add bindings from those exchanges to the queue?

Copy an input file into the directory relative to the root in default.yaml:
```
cp src/SmiServices/src/microservices/com.smi.microservices.ctpanonymiser/src/test/resources/image-000001.dcm /tmp
```

Create a fake message and send to TEST.ControlExchange:
```
python3 -m pip install pika
#!/usr/bin/env python3
msg_json = '{ "DicomFilePath": "image-000001.dcm", "ExtractionDirectory": "001/tmp/extractiondir/", "OutputPath": "output.dcm", "ExtractionJobIdentifier":"bb1cbed5-a666-4307-a781-5b83926eaa81", 
"ProjectNumber":"001", "ExtractionDirectory":"001/tmp", "JobSubmittedAt":"2019-12-19T10:49Z" }'
hdr={'MessageGuid':'', 'OriginalPublishTimestamp':'', 'ProducerExacutableName':'test.py', 'ProducerProcessID': '0'}
import pika
connection = pika.BlockingConnection(pika.ConnectionParameters('localhost'))
channel = connection.channel()
# exchange='TEST.ControlExchange', '' to make binding straight to routing_key queue
channel.basic_publish(exchange='', routing_key='TEST.ExtractFileQueue', body=msg_json, properties=pika.BasicProperties(content_type='application/json', headers=hdr) )
```

Run:
```
java -jar CTPAnonymiser-portable-1.0.0.jar -a dicom-whitelist.script.new -y default.yaml
```

The output is written to `/tmp/001/tmp/output.dcm` in this example and the log file is in `logs/YYYY-MM-DD-hhmmss.log`

A 'success' message is published to TEST.FileStatusExchange containing:
```
{"DicomFilePath":"image-000001.dcm","AnonymisedFileName":"output.dcm","Status":0,"ExtractionJobIdentifier":"bb1cbed5-a666-4307-a781-5b83926eaa81","ProjectNumber":"001","ExtractionDirectory":"001/tmp","JobSubmittedAt":"2019-12-19T10:49Z"}
```

# IsIdentifiable

See the netcoreapp2.2 branch of IsIdentifiable here:
https://github.com/HicServices/IsIdentifiable/tree/netcoreapp2.2
with the changes required to build and run on dotnet core 2.2 Linux
(until such time as it's merged into master).
# Data Flow

## Background

This document describes the flow of DICOM tag and pixel data through the SmiServices system.  This is an example deployment scenario, the actual implementation can be tailored according to needs.

![Where tags flow through various zones](./Images/dataflow.svg "Flow of tag and pixel data in the SMI codebase")

Key 

| Entity        | Purpose       |
| ------------- |:-------------:|
| Disk      | All original DICOM files are kept unchanged on disk |
| Mongo Database | All DICOM tags (except pixel data) are stored in JSON format |
| Pixel Algorithms | Validated production ready algorithms run on identifiable dicom pixel data and output results useful for cohort building into the relational database |
| NLP Algorithms | Validated production ready algorithms run on identifiable free text data (e.g. Dose Reports, Structured Reports) and output results useful for cohort building into the relational database |
| Relational Database | Only tags useful for cohort building that are easily (and reliably) anonymised (e.g. 5% of all tags) |
| Cohort building | Only tags useful for cohort building and only study/series level (e.g. 4% of all tags) |
| Researcher Zone (from CTP) | DICOM images with tags anonymised by CTP.  These include technical tags (some of which are not loaded to Relational / used in cohort building), date tags, anonymised patient ID e.g. 12% of original tags|
 | Researcher Zone (from NLP) | DICOM files containing full clinical reports (e.g. Dose Reports, Structured Reports). These reports would need to be be anonymised with a dedicated NLP tool as free text report redacting is not something CTP is set up to do.  A high proportion of the clinical report content would need to remain in these files for most free text research activities e.g. 98%|
# SmiServices Release Process

The steps to cut a new release of SmiServices are as follows.

All development is done via a simple branching workflow, which are merged into `master` via a reviewed PR. `master` therefore contains all the latest reviewed changes since the previous release, and the CI checks should always be passing. It is not possible to push to `master` directly.

The release workflow is to checkout a new `release/` branch from master, update the `CHANGELOG` etc. as per below, then open a release PR with just those updates. Once this is merged, a tag is pushed to `master`. This triggers a pipeline in Azure DevOps which creates a GitHub release. The other pipelines will then push artefacts to this release when they pass.

## Creating A Normal Release

-   Review all open PRs and check if any have been approved and can be merged to be included in the release.

-   Check that a [news file][news_files] is present for each merged PR since the previous release. To do this, checkout the latest `master` commit and list all the merged PRs since the last release, e.g.:
    ```console
    $ git checkout master && git pull
    $ git log  --oneline <previous_tag>.. | grep -vP "dependabot|Bump|pre-commit-ci" | grep -P '#\d+'
    ec182696 Merge pull request #430 from SMI/feature/extraction-fixes
    051a134e Merge pull request #444 from SMI/feature/trigger-updates
    65fcfe41 Merge pull request #440 from SMI/feature/value-updater
    8515f059 Merge pull request #438 from SMI/feature/consensus-tidyup
    38e25c2a Merge pull request #434 from SMI/feature/consensus-rule
    9d81b942 Merge branch 'master' into develop
    ee39d850 Merge pull request #408 from SMI/feature/isidentifiable-more-logs
    c709f7ed Merge pull request #404 from SMI/feature/no-an-suffix
    d7d90f4a Merge pull request #402 from SMI/feature/update-docs
    830bac67 Merge branch 'master' into develop
    ```
    Go through these PRs and check each has an accurate [news file][news_files] entry. Create any missing files if needed.

-   Identify the next release version. This can be determined by looking at the previous release and deciding if the new code to be released is a major, minor, or patch change as per [semver](https://semver.org). E.g. if the previous release was `v1.2.3` and only new non-breaking features are in the news files directory, then the next release should be`v1.3.0`. The definition of "breaking" can often be subjective though, so ask other members of the project if you're unsure.

-   Ensure you are on the latest commit on the `master` branch , and create a new release branch:

    ```console
    $ git fetch
    $ git status
    On branch master
    Your branch is up to date with 'origin/master'.

    nothing to commit, working tree clean

    $ git checkout -b release/v1.2.3
    Switched to a new branch 'release/v1.2.3'
    ```

-   Update the [CHANGELOG](/CHANGELOG.md) for the new release. This involves running the `utils/updateChangelog.py` script. Review the diff and check for any obvious errors.

-   Update any other files referencing the version. To see an example, check the previous release PR. At time of writing, these are:
    -   `README.md`: Bump the version in the header
    -   `src/SharedAssemblyInfo.cs`: Bump the versions in each property

-   Commit these changes and push the new branch
-   Open a PR for this branch with the title `Release <version>`. Request a review from `@tznind` and `@rkm`
-   If there are any further changes which need to be included in the release PR, then these can be merged into the release branch from `master`
-   Wait for the PR to be reviewed and merged
-   Checkout `master` and pull the merge commit
-   Tag the release, e.g.:
    ```console
    $ git tag v1.2.3
    $ git push origin v1.2.3
    ```
-   Delete the release branch
-   Wait for Azure Pipelines to build the release
-   Check that the built binaries are added to the [releases](https://github.com/SMI/SmiServices/releases) page.
-   (Internal) Ping the Mattermost ~developers channel to let everyone know there is a release available, and to not start any long-running tasks

## Creating A Hotfix Release

Hotfixes are small patches which are created in response to some show-stopper bug in the previous release.

The process is similar to above, except:

-   The branch name should be `hotfix/v...`
-   The PR should be titled `Hotfix <version>`

<!-- Links -->

[news_files]: /news/README.md
# Data Loading

__All data in this demo is synthetic (generated by [BadDicom])__

## Contents

- [Background](#background)
  - [MongoDb]
  - [RelationalDb]
- [Preparation](#preparation)
  - [Publish Binaries](#publish-binaries)
- [MongoDb Loading Microservices](#mongodb-loading-microservices)
  - [DicomDirectoryProcessor](#dicomdirectoryprocessor)
  - [DicomTagReader](#dicomtagreader)
    - [Dead Letter Exchange](#dead-letter-exchange)
  - [DicomTagReader Continued](#dicomtagreader-continued)
  - [MongoDbPopulator](#mongodbpopulator)
- [RelationalDb Loading Microservices](#relationaldb-loading-microservices)
  - [DicomReprocessor](#dicomreprocessor)
  - [IdentifierMapper](#identifiermapper)
  - [DicomRelationalMapper](#dicomrelationalmapper)
    - [Installing RDMP](#installing-rdmp)
    - [Picking the Schema](#picking-the-schema)
    - [Building the load](#building-the-load)
  - [DicomRelationalMapper Continued](#dicomrelationalmapper-continued)

## Background

This document describes all the steps required to setup data load microservices and use them to load a collection of Dicom images.

Microservices are designed to execute in parallel and scale to support hundreds of millions of dicom image files.

The data load process populates two databases:

 - Mongo Db (identifiable)
 - Relational Db (anonymous)

### MongoDb

The Mongo Db database stores all (non pixel) dicom tags and file paths for all dicom images.  It can be used for understanding what data you have (modalities, date ranges, image types etc) and serves as the input source for the subsequent relational database.

MongoDb is used because it is designed to store wide (many columns) and tree structures (e.g. dicom tags with the value representation SQ - sequence).

### RelationalDb

The Relational Db (e.g. Sql Server, MySql, Oracle or Postgres) stores only the tags required for cohort creation and image extraction (e.g. StudyDate, [PatientID], Modality, StudyDescription etc).  The RelationalDb should only be loaded with images that are fit for release (can be anonymised) and only anonymised tags should be loaded (this includes performing identifier substitution e.g. for PatientID).

A relational database is used because it allows easier linking with other traditional EHR data (e.g. prescribing, biochemistry etc) held by a safehaven.


## Preparation

Download [BadDicom] and use it to generate some test images on disk:

```
BadDicom.exe c:\temp\testdicoms
```

![Test files in file explorer (windows)](./Images/DataLoading/testfiles.png)

Ensure Mongo Db is running e.g.:

```
C:\Program Files\MongoDB\Server\3.6\bin> ./mongod
```

Ensure RabbitMQ is running e.g.:

```
C:\Program Files\RabbitMQ Server\rabbitmq_server-3.7.3\sbin> .\rabbitmq-server.bat start
```

Ensure the target DBMS is running e.g.:

```
E:\mysql-5.7.19-winx64\bin> ./mysqld
```

Delete all RabbitMQ exchanges and queues:

```
http://127.0.0.1:15672/#/queues
```

Follow instructions listed in https://stackoverflow.com/a/52002145/4824531

### Publish Binaries

For each microservice run `dotnet publish -r win-x64` e.g.

```
E:\SmiServices\src\applications\Applications.DicomDirectoryProcessor> dotnet publish -r win-x64
```


## MongoDb Loading Microservices

The following process are responsible for loading the [MongoDb] with

### DicomDirectoryProcessor

Run `DicomDirectoryProcessor` with the directory you created test dicom files in e.g.:

```
E:\SmiServices\src\applications\Applications.DicomDirectoryProcessor\bin\AnyCPU\Debug\netcoreapp2.2\win-x64> .\DicomDirectoryProcessor.exe -d C:\temp\testdicoms
```

This may cause the following error:

```
Failed to construct host:
System.IO.FileNotFoundException: Could not find the logging configuration in the current directory (Smi.NLog.config),
```

Copy and modify (if needed) [Smi.NLog.config] to the binary directory


Run the application again, this time you should see:

```
Failed to construct host:
System.ApplicationException: The given control exchange was not found on the server: "TEST.ControlExchange"
```

Create the exchange:

---

![Create Exchange](./Images/DataLoading/TEST.ControlExchange.png)

---

This is the exchange by which you can send runtime messages (e.g. shutdown) to the service

Now when it is run you will see an error relating to another missing exchange (probably `TEST.AccessionDirectoryExchange`)

Create the following exchanges:

- TEST.AccessionDirectoryExchange
- TEST.FatalLoggingExchange

Now when running you should see an error:

```
2019-12-02 13:18:50.6045|FATAL|DicomDirectoryProcessorHost|Could not confirm message published after timeout|System.ApplicationException: Could not confirm message published after timeout
```

This is because there is no queue associated with the output exchange.  Create a queue `TEST.AccessionDirectoryQueue`

---

![Create Exchange](./Images/DataLoading/TEST.AccessionDirectoryQueue.png)

---
Bind the `TEST.AccessionDirectoryExchange` exchange with the queue `TEST.AccessionDirectoryQueue`:

---

![Bind Exchange To Queue](./Images/DataLoading/BindExchange.png)

---

Once you have done this you should see output from the program like:

```
PS E:\SmiServices\src\applications\Applications.DicomDirectoryProcessor\bin\AnyCPU\Debug\netcoreapp2.2\win-x64> .\DicomDirectoryProcessor.exe -d C:\temp\testdicoms
Bootstrapper -> Main called, constructing host
2019-12-02 13:25:45.6886| INFO|DicomDirectoryProcessorHost|Host logger created with SMI logging config|||
2019-12-02 13:25:45.7365| INFO|DicomDirectoryProcessorHost|Started DicomDirectoryProcessor:5468|||
2019-12-02 13:25:45.8932| INFO|DicomDirectoryProcessorHost|Creating basic directory finder|||
Bootstrapper -> Host constructed, starting aux connections
Bootstrapper -> Host aux connections started, calling Start()
2019-12-02 13:25:45.9435| INFO|BasicDicomDirectoryFinder|Starting directory scan of: C:\temp\testdicoms|||
2019-12-02 13:25:46.1277| INFO|BasicDicomDirectoryFinder|Directory scan finished|||
2019-12-02 13:25:46.1277| INFO|BasicDicomDirectoryFinder|Total messages sent: 10|||
2019-12-02 13:25:46.1277| INFO|BasicDicomDirectoryFinder|Largest stack size was: 10|||
2019-12-02 13:25:46.1277| INFO|BasicDicomDirectoryFinder|Averages:
NewDirInfo:     0ms
EnumFiles:      0ms
FirstOrDef:     0ms
FoundNewDir:    17ms
EnumDirs:       0ms
PushDirs:       0ms
|||
2019-12-02 13:25:46.1277| INFO|DicomDirectoryProcessorHost|Host Stop called: Directory scan completed|||
2019-12-02 13:25:46.4812| INFO|DicomDirectoryProcessorHost|Host stop completed|||
Bootstrapper -> Host started
Bootstrapper -> Exiting main
```

There should be 1 message per folder in the your test dicoms directory:

---

![10 messages queued](./Images/DataLoading/AfterAccessionDirectory.png)

---

If you use GetMessages in the rabbit MQ interface you can see what was the messages contain:

---

![Example message from output queue](./Images/DataLoading/PeekAccessionDirectory.png)

---

That's right, all this work was just to get a __directory listing__ into RabbitMQ! But now that you have the basics of creating exchanges / queues down it should be much easier to get the rest of the services running (see below).

To change the exchange/queue names you should edit `default.yaml` (ensuring your RabbitMQ server has the correct entries)

### DicomTagReader

Publish and run DicomTagReader (copy across [Smi.NLog.config] if needed) e.g.:

```
E:\SmiServices\src\microservices\Microservices.DicomTagReader\bin\AnyCPU\Debug\netcoreapp2.2\win-x64> ./DicomTagReader.exe
```

This should result in an error about `TEST.IdentifiableSeriesExchange`.  Create the following exchanges:

- TEST.IdentifiableSeriesExchange
- TEST.IdentifiableImageExchange

This should cause our old friend:

```
Could not confirm message published after timeout
```

Notice also that a queue message has still been consumed and we have 1 less message in the `TEST.AccessionDirectoryQueue`

---

![One message has been consumed](./Images/DataLoading/LostMessages.png)

_RabbitMQ queue graph are 1 less message available for processing_

---

Messages that cannot be processed are 'nacked' and not returned to the processing queue.  This prevents 'bad' messages getting served up repeatedly to consumers and degrading system performance.  To prevent message loss we can set up a dead letter exchange.


#### Dead Letter Exchange

Create an internal exchange and queue for the dead letters:

- DeadLetterExchange
  - DeadLetterQueue

Now create a policy for all queues to send nacked messages to this exchange:

---

![Configuring a dead letter exchange policy](./Images/DataLoading/DeadLetterExchange.png)

---

Run DicomTagReader again to force another failure (because we still have no bound output queue for our successfully processed messages)

This should result in the 'lost' message being sent to the dead letter queue:

---

![Dead letter queue now shows 1 message waiting (the failing message) and the input queue has the remaining 3 messages](./Images/DataLoading/DeadLetterExchangeAfterError.png)

_The unprocessed message now resides in the dead letter queue_

---

### DicomTagReader Continued

Create the output queues for the tag reader exchanges (make sure to bind them to the correct exchanges):

- TEST.IdentifiableSeriesExchange
  - TEST.IdentifiableSeriesQueue
- TEST.IdentifiableImageExchange
  - TEST.IdentifiableImageQueue

This should produce the following output:

```
Bootstrapper -> Main called, constructing host
2019-12-03 09:09:42.1992| INFO|DicomTagReaderHost|Host logger created with SMI logging config|||
2019-12-03 09:09:42.2524| INFO|DicomTagReaderHost|Started DicomTagReader:17464|||
Bootstrapper -> Host constructed, starting aux connections
2019-12-03 09:09:42.5133| INFO|SerialTagReader|Stopwatch implementation - IsHighResolution: True. Frequency: 10000000 ticks/s|||
Bootstrapper -> Host aux connections started, calling Start()
Bootstrapper -> Host started
Bootstrapper -> Exiting main
2019-12-03 09:09:43.2117| INFO|SerialTagReader|Sending 8 DicomFileMessage(s)|||
2019-12-03 09:09:43.2553| INFO|SerialTagReader|Sending 2 SeriesMessage(s)|||
2019-12-03 09:09:43.3095| INFO|SerialTagReader|Sending 8 DicomFileMessage(s)|||
2019-12-03 09:09:43.3353| INFO|SerialTagReader|Sending 2 SeriesMessage(s)|||
2019-12-03 09:09:43.4184| INFO|SerialTagReader|Sending 8 DicomFileMessage(s)|||
2019-12-03 09:09:43.4330| INFO|SerialTagReader|Sending 2 SeriesMessage(s)|||
```

The binary will not exit by default (it will wait for more messages).  Use Ctrl+C to trigger shutdown of the binary.

```
2019-12-03 09:10:59.0453| INFO|SerialTagReader|Lock released, no more messages will be processed|||
2019-12-03 09:10:59.0453| INFO|SerialTagReader|Average rates - enumerate dir (per acc. message): 0.001034s, file process: 0.002731s, send messages: 0.003534s, overall: 0.211906s|||
2019-12-03 09:10:59.0453| INFO|DicomTagReaderHost|Host Stop called: Ctrl+C pressed|||
2019-12-03 09:10:59.5284| INFO|DicomTagReaderHost|Host stop completed|||
```

After execution the queues should look like:

---

![Output queues with 1 message per image + 1 message per series](./Images/DataLoading/AfterDicomTagReader.png)

_Output queues from a successful run of DicomTagReader_

---

If you peek at the messages in the `TEST.IdentifiableImageExchange`.  You should see the JSON representation of a dicom image (tags only - no pixel data):

```
Exchange 	TEST.IdentifiableImageExchange
timestamp:	1575364183
delivery_mode:	2
headers:	
MessageGuid:	a5f2ad28-6f87-49b6-b417-e87c425d74c1
OriginalPublishTimestamp:	1575293146
Parents:	90429bb7-5650-4ea0-922e-082cfc7befcc
ProducerExecutableName:	DicomTagReader
ProducerProcessID:	17464
content_encoding:	UTF-8
content_type:	application/json
Payload
2802 bytes
Encoding: string
	
{"DicomFilePath":"testdicoms\\1987\\12\\6\\2.25.176347174691273338913144606255096043339.dcm","StudyInstanceUID":"2.25.124865355268738178415667314856224778478","SeriesInstanceUID":"2.25.268908360241165259234396963267293474168","SOPInstanceUID":"2.25.176347174691273338913144606255096043339","DicomDataset":"{\"00080008\":{\"vr\":\"CS\",\"val\":\"ORIGINAL\\\\PRIMARY\\\\AXIAL\"},\"00080016\":{\"vr\":\"UI\",\"val\":\"1.2.840.10008.5.1.4.1.1.7\"},\"00080018\":{\"vr\":\"UI\",\"val\":\"2.25.176347174691273338913144606255096043339\"},\"00080020\":{\"vr\":\"DA\",\"val\":\"19871206\"},\"00080021\":{\"vr\":\"DA\",\"val\":\"19871206\"},\"00080022\":{\"vr\":\"DA\",\"val\":\"19871206\"},\"00080030\":{\"vr\":\"TM\",\"val\":\"180631\"},\"00080031\":{\"vr\":\"TM\",\"val\":\"180631\"},\"00080032\":{\"vr\":\"TM\",\"val\":\"180631\"},\"00080060\":{\"vr\":\"CS\",\"val\":\"CT\"},\"00080061\":{\"vr\":\"CS\",\"val\":\"CT\"},\"00081030\":{\"vr\":\"LO\",\"val\":\"CT Thorax & abdo & pel\"},\"00100010\":{\"vr\":\"PN\",\"val\":\"LUCA Price\"},\"00100020\":{\"vr\":\"LO\",\"val\":\"3003863640\"},\"00100030\":{\"vr\":\"DA\",\"val\":\"19860330\"},\"00101010\":{\"vr\":\"AS\",\"val\":\"001Y\"},\"00101040\":{\"vr\":\"LO\",\"val\":\"76 Foggyley Place Brechin and Edzell Angus DD9 6ES\"},\"00180050\":{\"vr\":\"DS\"},\"00180060\":{\"vr\":\"DS\",\"val\":\"0\"},\"00180088\":{\"vr\":\"DS\"},\"00181149\":{\"vr\":\"IS\",\"val\":\"0\"},\"00181150\":{\"vr\":\"IS\",\"val\":\"0\"},\"00181151\":{\"vr\":\"IS\",\"val\":\"0\"},\"00181152\":{\"vr\":\"IS\",\"val\":\"0\"},\"00189311\":{\"vr\":\"FD\",\"val\":[0.0]},\"00189461\":{\"vr\":\"FL\",\"val\":[0.0]},\"0020000D\":{\"vr\":\"UI\",\"val\":\"2.25.124865355268738178415667314856224778478\"},\"0020000E\":{\"vr\":\"UI\",\"val\":\"2.25.268908360241165259234396963267293474168\"},\"00200011\":{\"vr\":\"IS\",\"val\":\"0\"},\"00200012\":{\"vr\":\"IS\",\"val\":\"0\"},\"00200032\":{\"vr\":\"DS\",\"val\":\"0\\\\0\\\\0\"},\"00201041\":{\"vr\":\"DS\"},\"00201208\":{\"vr\":\"IS\",\"val\":\"2\"},\"00201209\":{\"vr\":\"IS\",\"val\":\"4\"},\"00280002\":{\"vr\":\"US\",\"val\":[3]},\"00280004\":{\"vr\":\"CS\"},\"00280006\":{\"vr\":\"US\",\"val\":[0]},\"00280008\":{\"vr\":\"IS\",\"val\":\"1\"},\"00280010\":{\"vr\":\"US\",\"val\":[500]},\"00280011\":{\"vr\":\"US\",\"val\":[500]},\"00280030\":{\"vr\":\"DS\",\"val\":\"0.3\\\\0.25\"},\"00280100\":{\"vr\":\"US\",\"val\":[8]},\"00280101\":{\"vr\":\"US\",\"val\":[8]},\"00280102\":{\"vr\":\"US\",\"val\":[7]},\"00280103\":{\"vr\":\"US\",\"val\":[0]},\"00280301\":{\"vr\":\"CS\",\"val\":\"NO\"},\"00282110\":{\"vr\":\"CS\",\"val\":\"00\"},\"00282112\":{\"vr\":\"DS\",\"val\":\"1\"},\"00282114\":{\"vr\":\"CS\",\"val\":\"ISO_10918_1\"},\"00400253\":{\"vr\":\"SH\",\"val\":\"0\"},\"7FE00010\":{\"vr\":\"OB\"}}"}
```

_A JSON serialized dicom dataset in RabbitMQ (this is __synthetic test data__ made up by the [BadDicom] tool)_

### MongoDbPopulator

The next microservice is responsible for persisting the dicom tag data into a MongoDb database.

Install and launch MongoDb Compas e.g.:

```
C:\Users\tznind\AppData\Local\MongoDBCompassCommunity\MongoDBCompassCommunity.exe
```

Your MongoDb instance should be blank (contain no imaging datasets at least):


![Mongo Db Compass showing no user databases](./Images/DataLoading/MongoDbCompassAtStart.png)

Publish and run `MongoDbPopulator` (making sure to copy across [Smi.NLog.config] if required)

```
E:\SmiServices\src\microservices\Microservices.MongoDbPopulator\bin\AnyCPU\Debug\netcoreapp2.2\win-x64> .\MongoDbPopulator.exe
```

This should give the following error `Expected queue "TEST.MongoSeriesQueue" to exist`.

MongoDbPopulator is designed to read the outputs from `DicomTagReader` but these outputs can also be forked to the relational database loading services to load both databases simultaneously.  For now lets stick with loading MongoDb only (we can always start RelationalDb loading from mongo collections anyway).

Open `default.yaml` in the exe directory of MongoDbPopulator and set the `QueueName` entries under `MongoDbPopulatorOptions` to `TEST.IdentifiableSeriesQueue` and `TEST.IdentifiableImageQueue`

```yaml
MongoDbPopulatorOptions:
    SeriesQueueConsumerOptions:
        QueueName: 'TEST.IdentifiableSeriesQueue'
        QoSPrefetchCount: 1000
        AutoAck: false
    ImageQueueConsumerOptions:
        QueueName: 'TEST.IdentifiableImageQueue'
        QoSPrefetchCount: 10000
        AutoAck: false
    MongoDbFlushTime: 30 # Seconds
    FailedWriteLimit: 5
```

_The MongoDbPopulator section of default.yaml should now look like this_

Run `MongoDbPopulator` again.  It should result in the following:

```
Bootstrapper -> Main called, constructing host
2019-12-03 10:32:40.4703| INFO|MongoDbPopulatorHost|Host logger created with SMI logging config|||
2019-12-03 10:32:40.5180| INFO|MongoDbPopulatorHost|Started MongoDbPopulator:7616|||
Bootstrapper -> Host constructed, starting aux connections
Bootstrapper -> Host aux connections started, calling Start()
2019-12-03 10:32:41.1959| INFO|MongoDbPopulatorHost|Starting consumers|||
Bootstrapper -> Host started
Bootstrapper -> Exiting main
2019-12-03 10:32:41.2108| INFO|MongoDbPopulatorHost|Consumers successfully started|||
2019-12-03 10:33:11.1579| INFO|Microservices.MongoDBPopulator.Execution.MongoDbAdapter|Attempting bulk write of 6 documents to dicom.series|||
2019-12-03 10:33:11.1579| INFO|ImageMessageProcessor|Queue contains 24 message to write|||
2019-12-03 10:33:11.1838| INFO|Microservices.MongoDBPopulator.Execution.MongoDbAdapter|Attempting bulk write of 24 documents to dicom.image_CT|||
```

Your MongoDb instance should now have 2 new collections `image_CT` and `series`.  The queues should also be fully drained of messages.

![Mongo Db Compass showing imaging databases](./Images/DataLoading/MongoDbCompassAtEnd.png)

_Mongo Db after MongoDbPopulator has run_

## RelationalDb Loading Microservices

The following microservices are responsible for loading the [RelationalDb] with anonymised tag data (and file paths) for downstream cohort creation, linkage and extraction processes.

### DicomReprocessor

This application is responsible for fetching records from [MongoDb] collections and queuing them for processing in RabbitMQ.

Publish and run `DicomReprocessor` (making sure to copy across [Smi.NLog.config] if required)

```
E:\SmiServices\src\microservices\Microservices.DicomReprocessor\bin\AnyCPU\Debug\netcoreapp2.2\win-x64> .\DicomReprocessor.exe
```

This should display the following helpful prompt:

```
ERROR(S):
  Required option 'c, collection-name' is missing.
USAGE:
Normal Scenario:
  DicomReprocessor --collection-name image_CT
```

Add the missing parameter, this should be the name of the collection loaded by [MongoDbPopulator] e.g. `image_CT`

```
DicomReprocessor.exe -c image_CT --auto-run
```

This should give us the error:

```
BasicReturn for Exchange 'TEST.IdentifiableImageExchange' Routing Key 'reprocessed' ReplyCode '312' (NO_ROUTE)
````

This is because `DicomReprocessor` is designed to feed images identified in MongoDb to any number of downstream processes of which loading the [RelationalDb] is only one.  To this end it requires a RabbitMQ exchange that can handle routing keys (i.e. to different destination queues).

Since we only want to send the messages on to one queue we can create a `fanout` exchange as the destination.  Create a new `fanout` exchange with a destination queue:

- TEST.DicomReprocessorExchange `(fanout)`
  - TEST.DicomReprocessorQueue

Edit `default.yaml` and set the `DicomReprocessorOptions` to use the new exchange.

```yaml
DicomReprocessorOptions:
    ProcessingMode: 'ImageReprocessing'
    ReprocessingProducerOptions: 
        ExchangeName: 'TEST.DicomReprocessorExchange'
        MaxConfirmAttempts: 1
```

This should give the following successful output (see below) and there should be image messages in the `TEST.DicomReprocessorQueue`

```
Bootstrapper -> Main called, constructing host
2019-12-04 09:09:47.6406| INFO|DicomReprocessorHost|Host logger created with SMI logging config|||
2019-12-04 09:09:47.6793| INFO|DicomReprocessorHost|Started DicomReprocessor:34748|||
2019-12-04 09:09:47.8743| INFO|DicomReprocessorHost|Documents will be reprocessed to TEST.DicomReprocessorExchange on vhost / with routing key "reprocessed"|||
Bootstrapper -> Host constructed, starting aux connections
Bootstrapper -> Host aux connections started, calling Start()
2019-12-04 09:09:48.1199| INFO|Smi.Common.MongoDB.MongoQueryParser|No query specified, fetching all records in collection|||
2019-12-04 09:09:48.3837| INFO|MongoDbReader|Using MaxDegreeOfParallelism: 4|||
2019-12-04 09:09:48.3837| INFO|MongoDbReader|Batch size is: unspecified|||
2019-12-04 09:09:48.3837| INFO|MongoDbReader|Sleeping for 0ms between batches|||
2019-12-04 09:09:48.3837| INFO|MongoDbReader|Starting reprocess operation|||
2019-12-04 09:09:48.9120| INFO|MongoDbReader|Reprocessing finished or cancelled, time elapsed: 0:00:00.529213|||
2019-12-04 09:09:48.9120| INFO|DicomReprocessorHost|Total messages sent: 24|||
2019-12-04 09:09:48.9120| INFO|DicomReprocessorHost|Total failed to reprocess : 0|||
2019-12-04 09:09:48.9120| INFO|DicomReprocessorHost|Average documents processed per second: 45|||
2019-12-04 09:09:48.9120| INFO|MongoDbReader|Cancelling the running query|||
2019-12-04 09:09:48.9120| INFO|DicomReprocessorHost|Host Stop called: Reprocessing completed|||
2019-12-04 09:09:49.1352| INFO|DicomReprocessorHost|Host stop completed|||
Bootstrapper -> Host started
Bootstrapper -> Exiting main
```

## IdentifierMapper

The next component in the load is the `IdentifierMapper`.  It's job is to anonymise the [PatientID] tag in the JSON extracted by [DicomReprocessor].  This change only occurs in the messages in the rabbit (as they are written to the output queue).  This prepares them for loading into the [RelationalDb].  

At no point are the original Dicom files on disk opened or changed.


Publish and run `IdentifierMapper` (making sure to copy across [Smi.NLog.config] if required)

```
E:\SmiServices\src\microservices\Microservices.IdentifierMapper\bin\AnyCPU\Debug\netcoreapp2.2\win-x64> .\IdentifierMapper.exe
```

Amongst the error messages should be the following interesting bits:

```
2019-12-04 09:27:51.3405| INFO|Smi.Common.Helpers.MicroserviceObjectFactory|Successfully constructed Type 'Microservices.IdentifierMapper.Execution.Swappers.ForGuidIdentifierSwapper'|||
2019-12-04 09:27:51.3405| INFO|IdentifierMapperHost|Calling Setup on swapper|||
Failed to construct host:
System.ArgumentException: MappingTableName did not contain the database/user section:'MappingTable'
```

`IdentifierMapper` uses a [strategy pattern] to determine how identifiers are substituted.  The following implementations are provided out of the box:

  - [ForGuidIdentifierSwapper]
  - [TableLookupSwapper]

We will use the [ForGuidIdentifierSwapper] because it doesn't require us to create an identiifer mapping up front.  Open `default.yaml` and edit the `IdentifierMapperOptions` settings e.g.:

```yaml
IdentifierMapperOptions:
    QueueName: 'TEST.DicomReprocessorQueue'
    QoSPrefetchCount: 1000
    AutoAck: false
    AnonImagesProducerOptions: 
        ExchangeName: 'TEST.AnonymousImageExchange'
        MaxConfirmAttempts: 1
    MappingConnectionString: 'Server=localhost\sqlexpress;Integrated Security=true;Initial Catalog=MappingDatabase;'
    MappingDatabaseType: 'MicrosoftSQLServer'
    MappingTableName: 'MappingDatabase.MappingTable'
    TimeoutInSeconds: 600
    SwapColumnName: 'PatientID'
    ReplacementColumnName: 'GuidPatientID'
    SwapperType: 'Microservices.IdentifierMapper.Execution.Swappers.ForGuidIdentifierSwapper'
    AllowRegexMatching: false
```

Make sure the `QueueName` is set to the output queue of [DicomReprocessor] (e.g. `TEST.DicomReprocessorQueue`) and that the `MappingConnectionString` is correct for your Sql Server instance.

Create the output exchange and queue

- TEST.AnonymousImageExchange
  - TEST.AnonymousImageQueue

Run `IdentifierMapper` again with the new yaml settings.  It should complete and have written all messages to the output queue `TEST.AnonymousImageQueue`:

```
Bootstrapper -> Main called, constructing host
2019-12-04 09:49:34.1412| INFO|IdentifierMapperHost|Host logger created with SMI logging config|||
2019-12-04 09:49:34.1909| INFO|IdentifierMapperHost|Started IdentifierMapper:26312|||
2019-12-04 09:49:34.4561| INFO|IdentifierMapperHost|Not passed a swapper, creating one of type Microservices.IdentifierMapper.Execution.Swappers.ForGuidIdentifierSwapper|||
2019-12-04 09:49:34.4561| INFO|Smi.Common.Helpers.MicroserviceObjectFactory|Successfully constructed Type 'Microservices.IdentifierMapper.Execution.Swappers.ForGuidIdentifierSwapper'|||
2019-12-04 09:49:34.4561| INFO|IdentifierMapperHost|Calling Setup on swapper|||
2019-12-04 09:49:35.5203| INFO|Microservices.IdentifierMapper.Execution.Swappers.ForGuidIdentifierSwapper|Guid mapping table does not exist, creating it now|||
2019-12-04 09:49:35.5574| INFO|Microservices.IdentifierMapper.Execution.Swappers.ForGuidIdentifierSwapper|Guid mapping table exist (MappingTable)|||
2019-12-04 09:49:35.5574| INFO|Microservices.IdentifierMapper.Execution.Swappers.ForGuidIdentifierSwapper|Checking for column PatientID|||
2019-12-04 09:49:35.6333| INFO|Microservices.IdentifierMapper.Execution.Swappers.ForGuidIdentifierSwapper|Checking for column GuidPatientID|||
Bootstrapper -> Host constructed, starting aux connections
Bootstrapper -> Host aux connections started, calling Start()
Bootstrapper -> Host started
Bootstrapper -> Exiting main
```

If you look in your Sql Server database you should see a persistent record of the anonmised mapping.  This ensures that patients have consistent identifiers over time and no aliases are generated.

![Sql server table containing mapped identifiers](./Images/DataLoading/SqlServerIdentifierMapperMappingTable.png)

Notice that the PatientID is a primary key column to prevent aliases ever forming.  Multiple swappers can execute in parallel without risking aliases (e.g. due to race conditions) because lookup is performed in a single atomic transaction (SELECT if not exists INSERT).

If you peek at the messages in the `TEST.AnonymousImageQueue` you can see the new PatientID tag value:

```
Exchange 	TEST.AnonymousImageExchange
timestamp:	1575452976
delivery_mode:	2
headers:	
MessageGuid:	c1b441fe-eb20-4060-a1db-3d7a6b9cb5ec
OriginalPublishTimestamp:	1575293146
Parents:	91775bf2-196d-499f-b07a-26fa098d142a->68a6435d-33ac-4895-88e9-2e18d03671f9->13a78be9-8cbf-40a2-9e28-d7745561159f->204e1e39-c408-44c9-b484-edddf4e0b2a6
ProducerExecutableName:	IdentifierMapper
ProducerProcessID:	26312
content_encoding:	UTF-8
content_type:	application/json
Payload
2813 bytes
Encoding: string
	
{"DicomFilePath":"testdicoms\\1958\\6\\9\\2.25.24425987081552946032369646390860332875.dcm","StudyInstanceUID":"2.25.39490047024894726575859340458560242419","SeriesInstanceUID":"2.25.200404735591937354156655516808405291743","SOPInstanceUID":"2.25.24425987081552946032369646390860332875","DicomDataset":"{\"00080008\":{\"vr\":\"CS\",\"val\":\"ORIGINAL\\\\PRIMARY\\\\AXIAL\"},\"00080016\":{\"vr\":\"UI\",\"val\":\"1.2.840.10008.5.1.4.1.1.7\"},\"00080018\":{\"vr\":\"UI\",\"val\":\"2.25.24425987081552946032369646390860332875\"},\"00080020\":{\"vr\":\"DA\",\"val\":\"19580609\"},\"00080021\":{\"vr\":\"DA\",\"val\":\"19580609\"},\"00080022\":{\"vr\":\"DA\",\"val\":\"19580609\"},\"00080030\":{\"vr\":\"TM\",\"val\":\"085353\"},\"00080031\":{\"vr\":\"TM\",\"val\":\"085353\"},\"00080032\":{\"vr\":\"TM\",\"val\":\"085353\"},\"00080060\":{\"vr\":\"CS\",\"val\":\"CT\"},\"00080061\":{\"vr\":\"CS\",\"val\":\"CT\"},\"00081030\":{\"vr\":\"LO\",\"val\":\"CT Head\"},\"00100010\":{\"vr\":\"PN\",\"val\":\"TOMMY Moore\"},\"00100020\":{\"vr\":\"LO\",\"val\":\"8f774005-6fe6-4a31-96d7-c454fbf7e323\"},\"00100030\":{\"vr\":\"DA\",\"val\":\"19560904\"},\"00101010\":{\"vr\":\"AS\",\"val\":\"001Y\"},\"00101040\":{\"vr\":\"LO\",\"val\":\"0 Tulloch Court Arbroath East and Lunan Angus  DD11 4RU\"},\"00180050\":{\"vr\":\"DS\"},\"00180060\":{\"vr\":\"DS\",\"val\":\"0\"},\"00180088\":{\"vr\":\"DS\"},\"00181149\":{\"vr\":\"IS\",\"val\":\"0\"},\"00181150\":{\"vr\":\"IS\",\"val\":\"0\"},\"00181151\":{\"vr\":\"IS\",\"val\":\"0\"},\"00181152\":{\"vr\":\"IS\",\"val\":\"0\"},\"00189311\":{\"vr\":\"FD\",\"val\":[0.0]},\"00189461\":{\"vr\":\"FL\",\"val\":[0.0]},\"0020000D\":{\"vr\":\"UI\",\"val\":\"2.25.39490047024894726575859340458560242419\"},\"0020000E\":{\"vr\":\"UI\",\"val\":\"2.25.200404735591937354156655516808405291743\"},\"00200011\":{\"vr\":\"IS\",\"val\":\"0\"},\"00200012\":{\"vr\":\"IS\",\"val\":\"0\"},\"00200032\":{\"vr\":\"DS\",\"val\":\"0\\\\0\\\\0\"},\"00201041\":{\"vr\":\"DS\"},\"00201208\":{\"vr\":\"IS\",\"val\":\"2\"},\"00201209\":{\"vr\":\"IS\",\"val\":\"4\"},\"00280002\":{\"vr\":\"US\",\"val\":[3]},\"00280004\":{\"vr\":\"CS\"},\"00280006\":{\"vr\":\"US\",\"val\":[0]},\"00280008\":{\"vr\":\"IS\",\"val\":\"1\"},\"00280010\":{\"vr\":\"US\",\"val\":[500]},\"00280011\":{\"vr\":\"US\",\"val\":[500]},\"00280030\":{\"vr\":\"DS\",\"val\":\"0.3\\\\0.25\"},\"00280100\":{\"vr\":\"US\",\"val\":[8]},\"00280101\":{\"vr\":\"US\",\"val\":[8]},\"00280102\":{\"vr\":\"US\",\"val\":[7]},\"00280103\":{\"vr\":\"US\",\"val\":[0]},\"00280301\":{\"vr\":\"CS\",\"val\":\"NO\"},\"00282110\":{\"vr\":\"CS\",\"val\":\"00\"},\"00282112\":{\"vr\":\"DS\",\"val\":\"1\"},\"00282114\":{\"vr\":\"CS\",\"val\":\"ISO_10918_1\"},\"00400253\":{\"vr\":\"SH\",\"val\":\"0\"},\"7FE00010\":{\"vr\":\"OB\"}}"}
```

The critical section in this JSON is `\"00100020\":{\"vr\":\"LO\",\"val\":\"8f774005-6fe6-4a31-96d7-c454fbf7e323\"}`.  The [PatientID] tag in dicom is `00100020` and we can see the that the value is the guid assigned by the swapper.

### DicomRelationalMapper

`DicomRelationalMapper` is responsible for loading the [RelationalDb] with the JSON serialized dicom images in it's RabbitMQ input queue.  

Loading research ready relational databases without introducing duplication is complicated and highly dependent on user requirements (e.g. desired table schema, any aggregate/computed columns etc).

To ensure maximum flexibility `DicomRelationalMapper` wraps the [RDMP] Data Load Engine (using the [Rdmp.Dicom] plugin)

#### Installing RDMP
    
In order to set up the data load (but not run it) you will need to install [RDMP].  This can be done through the windows client by following the instructions in the [RDMP User Manual].

Alternatively you can use the [command line client](https://github.com/HicServices/RDMP/releases).  Download and unzip the latest CLI package for your operating system (e.g. `rdmp-cli-win-x64.zip`)

Run the install command (if you need to use sql authentication use the -u and -p flags too)

```
./rdmp.exe install localhost\sqlexpress TEST_
```

This will create all the databases required for [RDMP] to run (including running data loads)

![RDMP platform databases installed on localhost sql server](./Images/DataLoading/RdmpPlatformDatabases.png)

Edit `Databases.yaml` in the RDMP CLI directory so that the connection strings are correct for your server e.g.

```yaml
CatalogueConnectionString: Server=localhost\sqlexpress;Database=TEST_Catalogue;Trusted_Connection=True;
DataExportConnectionString: Server=localhost\sqlexpress;Database=TEST_DataExport;Trusted_Connection=True;
```

Next download the latest [Rdmp.Dicom] plugin version compatible with your RDMP binary (See the release notes of the plugin for version compatibility).  For example [Rdmp.Dicom.2.0.4.nupkg](https://github.com/HicServices/RdmpDicom/releases/tag/v2.0.4)

Add the plugin to RDMP with the CLI pack command (or through the [windows client](https://github.com/HicServices/RDMP/blob/develop/Documentation/CodeTutorials/Images/AddPluginContextMenu.png))

```
./rdmp pack --file C:\Users\tznind\Downloads\Rdmp.Dicom.2.0.4.nupkg
``` 

You should now find the command `CreateNewImagingDatasetSuite` is listed when you run

```
./rdmp cmd ListSupportedCommands
```

This confirms that both [RDMP] and the [Rdmp.Dicom] plugin are correctly installed.  If you don't see the command, try running with `--logstartup` to see any errors.

#### Picking the Schema

The [RelationalDb] can have any schema you want and of any [DBMS] supported by [RDMP].  `DicomRelationalMapper` will load all tables configured and populate any columns which are named after dicom tags (e.g. [PatientID], StudyDate, SeriesDate etc).  You can define a single large table or tables at different levels of granularity (e.g. Study / Series / Image).  Only unique records will be loaded (i.e. a study table would contain 1 entry per study not 1 per image with rampant duplication).

Schemas are defined in yaml.  There are a number of [existing templates](https://github.com/HicServices/DicomTypeTranslation/tree/master/Templates) you can pick from.  Alternatively you can build your own using the [Dicom Template Builder].

Download the [CT.it](https://raw.githubusercontent.com/HicServices/DicomTypeTranslation/master/Templates/CT.it) template.  Copy the contents into a [yaml validator](http://www.yamllint.com/) to ensure there were no problems with leading/trailing whitespace etc when you downloaded it.

#### Building the load
  
Once you have picked a schema you will need to run the `CreateNewImagingDatasetSuite` command in RDMP.  This can be [done through the windows client](https://github.com/HicServices/RdmpDicom/blob/develop/Documentation/DataLoad.md).

Alternatively you can run it from the rdmp console gui:

```
./rdmp gui
```

Press `F9` to access the menu and select `R` (Run).  Locate the `CreateNewImagingDatasetSuite` command and run it:

![Create suite command listed in rdmp gui](./Images/DataLoading/RdmpGuiCreateSuite.png)

Enter the server/database name you want the tables created into.  __Make sure to select Create Database__ (if the database does not already exist).

![Create the target database in rdmp gui](./Images/DataLoading/RdmpGuiCreateDatabase.png)

- Select a 'project directory' this should be a directory on disk (that must exist) in which load scripts can be added later (all [RDMP] loads require a load folder regardless of whether they actually use it)
- Select DicomDatasetCollectionSource for the `dicomSourceType`
- Enter `CT_` when prompted for table prefix
- Select the `CT.it` template file when prompted

![Pick the template file in rdmp gui](./Images/DataLoading/RdmpGuiTemplateFile.png) 

- Enter Yes for `persistentRAW` (this allows parallel loading)
- Enter Yes for `createLoad`

Run Refresh (`F9` Refresh `f`)

Now open the tree (`F9` Tree `t`) and search for "loadmetadata"

![The load created in rdmp gui](./Images/DataLoading/RdmpGuiOpenTree.png) 

Select "SMI Image Loading CT", this will show all the operations created in the load

![The load created in rdmp gui detail](./Images/DataLoading/RdmpGuiLoadCreated.png) 

Take note of the ID of the load (in this case `1`).  You can see this in the search or by opening the load.

You can check the load by exiting (`F9` Quit `Q`) and running the command line checks

```
 ./rdmp dle -l 1 --command check
```

In Sql Server you should see tables (and archive tables) that match your template

![The imaging tables and archive tables](./Images/DataLoading/RdmpAllDatabases.png) 

### DicomRelationalMapper Continued

Publish (but do not run) `DicomRelationalMapper` (making sure to copy across [Smi.NLog.config] if required)

Now that we have a valid load configured in [RDMP] we can reference it in `default.yaml`.  We need to set the connection strings to the RDMP databases

```yaml
RDMPOptions:
    CatalogueConnectionString: 'server=localhost\sqlexpress;integrated security=true;database=TEST_Catalogue'
    DataExportConnectionString: 'server=localhost\sqlexpress;integrated security=true;database=TEST_DataExport'
```

Next set the ID of the load created above (e.g. 1).  Also make sure to set `MinimumBatchSize` to 1

```yaml
DicomRelationalMapperOptions:
    Guid: '6ff062af-5538-473f-801c-ed2b751c7897'
    QueueName: 'TEST.AnonymousImageQueue'
    QoSPrefetchCount: 10000
    AutoAck: false
    LoadMetadataId: 1
    DatabaseNamerType: 'GuidDatabaseNamer'
    MinimumBatchSize: 1
    UseInsertIntoForRAWMigration: true
    RetryOnFailureCount: 1
    RetryDelayInSeconds: 60
    RunChecks: true
```

Now run `DicomRelationalMapper` e.g.

```
E:\SmiServices\src\microservices\Microservices.DicomRelationalMapper\bin\AnyCPU\Debug\netcoreapp2.2\win-x64> .\DicomRelationalMapper.exe
```

After it finishes executing (which will produce copious logs) you should see an entry along the lines of:

```
2019-12-04 12:31:46.1759| INFO|Microservices.DicomRelationalMapper.Execution.ParallelDLEHost|Migrate table [ImagingDb].dbo.[CT_ImageTable] from STAGING to ImagingDb: 23 inserts, 0 updates|||
```
_Successful loading of images to the relational database (23 images because it ran two batches in sequence - 1 of 1 and one of 23)._

Your live tables should have the following (or rough equivellents):

![Final state of the study level tags](./Images/DataLoading/FinalStudyTable.png) 

_Final live study table containing aggregate tag data at the study level_

![The final state of the image level](./Images/DataLoading/FinalImageTable.png) 

_Final live image table containing an entry for each image_

Notice that there is a column `RelativeFileArchiveURI`.  This contains the image paths that we originally loaded:

```
testdicoms/1963/11/18/2.25.106207137390777977751707424522924063465.dcm
testdicoms/1958/6/9/2.25.110656297039029576033942135898101468288.dcm
testdicoms/1963/11/18/2.25.166202323220770548415572214738761810149.dcm
testdicoms/1958/6/9/2.25.171673898064475963802645879189030191609.dcm
testdicoms/1987/12/6/2.25.176347174691273338913144606255096043339.dcm
[...]
```

[Smi.NLog.config]: ../data/logging/Smi.NLog.config
[BadDicom]: https://github.com/HicServices/BadMedicine.Dicom/releases
[MongoDb]: #mongodb
[RelationalDb]: #relationaldb
[MongoDbPopulator]: #MongoDbPopulator
[DicomReprocessor]: #DicomReprocessor
[PatientID]: https://dicom.innolitics.com/ciods/rt-plan/patient/00100020
[strategy pattern]: https://en.wikipedia.org/wiki/Strategy_pattern
[ForGuidIdentifierSwapper]: ../src/microservices/Microservices.IdentifierMapper/Execution/Swappers/ForGuidIdentifierSwapper.cs 
[TableLookupSwapper]: ../src/microservices/Microservices.IdentifierMapper/Execution/Swappers/TableLookupSwapper.cs 
[RDMP]: https://github.com/HicServices/RDMP
[Rdmp.Dicom]: https://github.com/HicServices/RdmpDicom
[RDMP User Manual]: https://github.com/HicServices/RDMP#research-data-management-platform
[DBMS]: https://github.com/HicServices/RDMP/blob/develop/Documentation/CodeTutorials/Glossary.md#dbms
[Dicom Template Builder]: https://github.com/HicServices/DicomTemplateBuilder

Note: This was copied from a rough spec. document which was last updated on 2018-11-29.

# Image Extraction Flags & Refactoring

__Version 2.0 - 2018-11-29__

In the initial dataset we published to, only CT images with ImageType like `ORIGINAL\PRIMARY` were copied from our internal catalogue. This was our initial method for controlling what images were valid for extraction. It is important however that researchers can be informed of the data that will be included (and excluded) from any cohort that is generated, and also that we can disable certain images from being extracted due to 'dodgy' or corrupt data.

This document describes the proposed changes to our catalogue schemas and extraction process to allow for it. It also includes changes and bug-fixes which were discovered during testing.

## Metadata catalogue changes

### Image file table

We have discussed this when dealing with managing separate tables for each modality, but the general idea is to have a single table containing all the dicom file paths and any other information that is related to the specific file rather than any piece of metadata. We can use this here to store a flag for specific images we have deemed "not extractable". The schema might look something like:

```
SOPInstanceUID 			- VARCHAR(64) NOT NULL PRIMARY KEY
RelativeFileArchiveURI	- VARCHAR(512) NOT NULL
ExtractableFlag 		- BIT DEFAULT TRUE
ExtractableReason 		- TEXT DEFAULT NULL
```

In the case where an image is loaded through the extraction pipeline again after previously being marked not extractable, the flag should remain unset and not be reset to the default. This won't be needed if we also mark the image as not extractable in MongoDB in some way.

To aid the RC team, we could also add an Extractable flag at the Series/Study level which would indicate that every image in series / series in study has been 'disabled'. This could be automatically generated as part of a stored procedure, or when image(s) are manually disabled.

### Extractable white-list rules table

This is the other part of our definition what is extractable. This would consist of a very simple table containing rules that can be used to determine if an image is allowed to be extracted. An image would be considered extractable if the result of applying all the rules in a row to its dicom tags is true, for any row. For our current rules, we would have something like:

ImageType | Modality
---|---
LIKE '%ORIGINAL\\PRIMARY%' | EQ 'CT'
LIKE '%ORIGINAL\\PRIMARY%' | EQ 'MR'

Exact syntax of this to be decided. For example, any given image would be extractable if matched the filter `ImageType LIKE '%ORIGINAL\\PRIMARY%' AND Modality = 'CT'` or `ImageType LIKE '%ORIGINAL\\PRIMARY%' AND Modality = 'MR'`.

- ❓ How would this evolve over time? (performance for many rules, ordering of rules etc.)
- ❓ This should be fairly easy to deal with in the application layer, but how easy will it be calculate if an image is extractable from a database query?
- ❓ How to manage adding new columns to this?

## Extraction changes

Majority of this is software-related, however we should really create a short extraction guide for the RCs to use. It should include:

- Valid identifier formats permitted in the ID request file (i.e. "SOPInstanceUID" not "imageId")
- Notes on use of extractable flags

### CohortExtractor

At extraction, the CohortExtractor will start and load the table of white-list rules. It will then query the catalogues as normal and return a set of file paths matching the identifier(s) given in the message. It will then apply the extractable rules to the returned data; firstly checking if the ExtractableFlag is set, then applying the rules from the loaded white-list. It will then emit messages for any extractable images as normal, but also emit audit messages for anything not extractable, which will be recorded by the CohortPackager. This will require some messages to be updated.

As an aside, we would also like to refactor our output directory format to match `<ExtractRoot>/<EUPI>/<StudyInstanceUID>/<SeriesInstanceUID>/*.dcm` as standard. This means we need 2 dicom UIDs no matter which one we are using as the key for extraction.

This would mean we need to join across 3 tables to get all the information required for extracting an image (Study/Series/SOPInstanceUID, FilePath, Extractable info).

### CohortPackager

The CohortPackager now needs to handle the messages where we don't expect an anonymised image to be produced. We should be able to provide a summary of number of images requested/extractable/extracted and generate some sort of report of not extracted / reason / counts. For auditing, we lso need to be able to determine which extractions contain a specific image if requested.
