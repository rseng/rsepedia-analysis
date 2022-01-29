Contact [Stephanie Singer](http://symmetrysinger.com/index.php?id=contact). [![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/4078/badge)](https://bestpractices.coreinfrastructure.org/projects/4078)

# Overview
This repository provides tools for consolidation and analysis of raw election results from the most reliable sources -- the election agencies themselves. 
 * Consolidation: take as input election results files from a wide variety of sources and load the data into a relational database
 * Export: create consistent-format export files of results sets rolled up to any desired intermediate geography
   * tabular (tab-separated text)
   * xml (following NIST Election Results Reporting Common Data Format V2)
   * json (following NIST Election Results Reporting Common Data Format V2)
 * Analysis: 
   * Curates one-county outliers of interest
   * Calculates difference-in-difference for results available by vote type
 * Visualization: 
   * Scatter plots
   * Bar charts

# Target Audience
This system is intended to be of use to news media, campaigns, election officials, students of politics and elections, and anyone else who is interested in assembling and understanding election results. If you have ideas for using this system or if you would like to stay updated on the progress of this project, [we'd like to hear from you](CONTACT_US.md). 

# How to use the app
See [documentation directory](docs), which includes
 * for users
   * [Installation instructions](docs/Installation.md)
   * Instructions for a [sample dataloading session](docs/Sample_Session.md)
   * Detailed [User Guide](docs/User_Guide.md)
 
# How to Contribute Code
See [CONTRIBUTING.MD](CONTRIBUTING.md).

# Contributors
 * [Stephanie Singer](http://campaignscientific.com/), Hatfield School of Government (Portland State University), former Chair, Philadelphia County Board of Elections
 * Janaki Raghuram Srungavarapu, Hatfield School of Government (Portland State University)
 * Eric Tsai, Hatfield School of Government (Portland State University)
 * Bryan Loy
 * Jon Wolgamott
 * Elliot Meyerson

# Funding
Funding provided October 2019 - November 2021 by the National Science Foundation
 * Award #1936809, "EAGER: Data Science for Election Verification" 
 * Award #2027089, "RAPID: Election Result Anomaly Detection for 2020"
Data collection and consolidation for the 2020 US General Election funded in part by the Verified Voting Foundation.

# License
See [LICENSE.md](./LICENSE.md)


# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who contribute through reporting issues, posting feature requests, updating documentation, submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for everyone, regardless of level of experience, gender, gender identity and expression, sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or imagery, derogatory comments or personal attacks, trolling, public or private harassment, insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the [Contributor Covenant](http:contributor-covenant.org), version 1.0.0, available at https://www.contributor-covenant.org/version/1/0/0/code-of-conduct.html
Please note that electiondata is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). 
By contributing to this project, 
you agree to abide by its terms.

# Contributing to electiondata development

## Report a Problem

To report a problem, [file an issue on GitHub](https://github.com/ElectionDataAnalysis/electiondata/issues). When filing an issue, the most important thing is to include a minimal 
reproducible example so that we can quickly verify the problem, and then figure 
out how to fix it. Please include

1.  Link to the exact version of electiondata you are using.
  
1.  A copy of your working directory, following the model in the [Sample Session](docs/Sample_Session.md). Include:
    * a `run_time.ini` file (Feel free to redact the login information for your local postgres instance.) 
    * all subdirectories referenced in the `run_time.ini` file. 
    * If the issue involves loading a particular results data file, be sure to include that file. If possible, avoid submitting large results files -- if your file is precinct-based, for example, see if you can demonstrate the issue with a truncated results file.
    * all error and warnings files placed by the system into the `reports_and_plots_dir` specified in `run_time.ini`.
  
1.  A transcript of the python session, where python is called from the working directory.

You can check you have actually made a reproducible example by:
1. creating a virtual environment
1. installing the indicated version of electiondata
1. if files or folders were moved by the system to the archive directory, move them back to the input directory, removing any timestamps from directory names.
1. navigating to the working directory, calling python, and producing the behavior in question.

## Revise the Code

To contribute a change to `electiondata` follow these steps:

1. Create a branch in git and make your changes.
1. Run unit tests to check for broken functionality. See [pytest instructions](docs/Testing_Code_with_pytest.md).
1. Run python `black` to format your code to our standard.
1. Push branch to github and issue pull request (PR).
1. Discuss the pull request.
1. Iterate until either we accept the PR or decide that it's not
   a good fit for `electiondata`.

Each of these steps are described in more detail below. This might feel 
overwhelming the first time you get set up, but it gets easier with practice. 
If you get stuck at any point, feel free to [contact us](CONTACT_US.md) for help.

If you're not familiar with git or github, please read a tutorial such as [https://realpython.com/python-git-github-intro/](https://realpython.com/python-git-github-intro/).


Pull requests will be evaluated against this checklist:

1.  __Motivation__. Your pull request should clearly and concisely motivate the
    need for change.

1.  __Only related changes__. Before you submit your pull request, please
    check to make sure that you haven't accidentally included any unrelated
    changes. These make it harder to see exactly what's changed, and to
    evaluate any unexpected side effects.

    Each PR corresponds to a git branch, so if you expect to submit
    multiple changes make sure to create multiple branches. If you have
    multiple changes that depend on each other, start with the first one
    and don't submit any others until the first one has been processed.
    
1.  __Documentation__ Any new parameters or a new functions must be documented both in the code and in the [User Guide](docs/User_Guide.md), [Sample Session](docs/Sample_Session.md) and any other relevant documents. If you're adding a new graphical or analytical feature, please add a short example to [Sample Session](docs/Sample_Session.md).

1.  __Tests__ If fixing a bug or adding a new feature to a non-graphical function,
    please add a [pytest](https://docs.pytest.org) unit test. Document the new test in [pytest instructions](docs/Testing_Code_with_pytest.md).

This seems like a lot of work but don't worry if your pull request isn't perfect.
Unless you've submitted a few in the
past it's unlikely that your pull request will be accepted as is. All PRs require
review and approval from at least one member of the `electiondata` development team 
before merge.

# Environment
You will need:
 * `python3.9`. 
 * `postgresql`. You do *not* need to create a database; you just need to make sure a `postgresql` is installed, and you need to have the hostname, port and login credentials with default permissions. 
 * all required python packages: run `python3.9 -m pip install -r requirements.txt` from the root folder of the repository.
   
For more detail, see the [Installation Instructions](./Installation.md).

To use other varieties of SQL, you would need to modify the routines in the `database` module.

# Installation
The one-line version: From the root folder of your repository run `python3 setup.py install`. For more detail, see the [Installation Instructions](./Installation.md).

## Main parameter file
You will need a main parameter file to specify paths and database connection information specific to your local computing environment. This file is necessary for the three main classes:
 * `JurisdictionPrepper` for preparing jurisdiction files
 * `DataLoader` for loading data
 *  `Analyzer` for exporting and analyzing results
See the [template file](../src/parameter_file_templates/run_time.ini.template) for required parameters. Avoid percent signs and line breaks in the parameter values.
   
## Other recommended files
To avoid the overhead of deriving the major subdivision type for each jurisdiction from the database, make sure that your repository has a [000_major_subjurisdiction_types.txt](../src/jurisdictions/000_for_all_jurisdictions/major_subjurisdiction_types.txt) in the [jurisdictions directory](../src/jurisdictions/). This file allows the user to specify other major subdivisions. For example, it may make sense to consider towns as the major subdivisions in Connecticut rather than counties. Or a user may wish to use congressional districts as the major subdivision -- though such a user should not assume that the nesting relationships (say, of precincts within congressional districts) have been coded in the [`ReportingUnit.txt` file](../src/jurisdictions/Connecticut/ReportingUnit.txt) or the database.

# Data Load Preparation
In order to process a raw results data file, the system needs specified:
   * a "munger" with file format details -- e.g., where to find the candidate name associated with a particular vote count 
   * a "dictionary.txt" file to translate from naming conventions of the file into the naming conventions for the database
   * other "jurisdiction files" listing legitimate contests, candidates, etc., with information about each (e.g. election district or party)

## Mungers
Election result data comes in a variety of file formats. Even when the basic format is the same, file columns may have different interpretations. The code is built to ease -- as much as possible -- the chore of processing and interpreting each format. Following the [Jargon File](http://catb.org/jargon/html/M/munge.html), which gives one meaning of "munge" as "modify data in some way the speaker doesn't need to go into right now or cannot describe succinctly," we call each set of basic information about interpreting an election result file a "munger". 

If the munger for the format of your results file doesn't already exist:
 * pick a name for your munger 
 * create a file with that name and extension `.munger` in the `mungers` directory (e.g., `me_excel.munger`) with sections and parameters described below. You may find it helpful to work with the template from `src/mungers/000_template.munger`. 
 
 The file with munger parameters has one or more sections, each with a header:
  * (required) `[format]` for the main parameters
  * (required) `[munge formulas]`
  * (required if formulas use foreign keys) one section for each foreign key appearing in munge formulas. Foreign keys are preceded by " from  ", e.g., in the formula
  * (optional) `[ignore]`
 
### \[format\]
 There are two required format parameters: `file_type` and `count_location`. 
 The `file_type` parameter controls which function from the python `pandas` module reads the file contents. Related optional and required parameters must be given under the `[format]` header. Acceptable values are 'flat_text', 'excel', 'xml', 'json-nested'. The `count_location` parameter indicates where the vote counts are to be found. For 'flat_text' or 'excel' file types, either `count_location=by_name:<list of names of columns containing vote counts>` or `count_location=by_number:<list of positions of columns containing vote counts`. 
  * 'flat_text': Any tab-, comma-, or other-separated table in a plain tabular text file.
    * (required) a field delimiter `flat_text_delimiter` to be specified (usually `flat_text_delimiter=,` for csv or `flat_text_delimiter=tab` for .txt)
      
  * 'excel'
    * (optional) a list `sheets_to_read_names` (and/or `sheets_to_read_numbers`) of spreadsheets to read, 
    * (optional) a list `sheets_to_skip_names` of names of spreadsheets to skip
    * Default is to read all sheets
    * (optional) `merged_cells` If there are merged cells in the meaningful header rows, set `merged_cells=yes`. 
  * for both 'flat_text' and 'excel':
    * (required if `count_location=by_name`) specify location of field names for count columns. with integer `count_field_name_row` (NB: top row not skipped is 0, next row is 1, etc.)
    * (required):
        * Either `all_rows=data` or designate row containing column names for the candidate, reporting unit, etc. with the `noncount_header_row` parameter. (NB: top row not skipped is 0, next row is 1, etc.)
    
  *  'xml'

  * 'json-nested'

   Available if appropriate for any file type, under the `[format]` header:
   * (required if any munging information needs to be read from the `<results>.ini` file) `constant_over_file`, a comma-separated list of elements to be read, e.g., `constant_over_file=CandidateContest,CountItemType`.
   * (optional) `thousands_separator`. In the US the separator is almost always ',' if it is used. Watch out for Excel files which may show a comma when you look at them in Excel -- there may not be a comma in the underlying data.
   * (optional) `encoding` (If not specified or recognized, a default encoding will be used. Recognized encodings are limited [python's list of recognized encodings and aliases](https://docs.python.org/3/library/codecs.html#standard-encodings).)

   Available for flat_text and Excel file types:
   * (optional - use only for multi-block) `rows_to_skip` An integer giving the number of rows to skip at the top before processing blocks.
   * (optional) `all_rows` If the file has no column headers but only data rows with counts, set this parameter to 'data'
   * (optional) `multi_block` if there are multiple blocks of data per page, each with its own headers, set this parameter to 'yes'. For multi-block sheets, munge parameters refer to the blocks (and must be the same for all blocks).
   * (optional) `max_blocks` if `multi_block=yes`, `max_blocks` is an integer telling the system how many blocks at most to read off of each sheet.
 
### \[munge formulas\]
Put each formula for parsing information from the results file into the `[munge formulas]` section. Constant items can be give either:
 * as comma separated list in constant_over_sheet parameter in .munger file, with values in the .ini file
 * as a constant formula in the `[munge formulas]` section, in which case a corresponding entry must be made in the jurisdiction's `dictionary.txt`.

For many results files, it is enough to create concatenation formulas, referencing field names from your file by putting them in angle brackets (<>. The available fields are:
  * `<count_header_0>`, `<count_header_1>`, etc. to denote information from the headers of the columns containing the counts.
  * `<row_0>`, `<row_1>`, etc., to read information constant over a sheet or block from one of the header rows. The system recognizes the leftmost non-blank cell as the content to be read.
  * `<sheetname>` to denote the name of an Excel spreadsheet within a (possibly multi-sheet) workbook.
  * any field name from the file itself, e.g., `<SELECTION>`
  * for tree-structured files like json and xml, see below for field labeling conventions.
  
Some characters are reserved to indicate parsing and interpretation in the munge formulas. Angle brackets (`<>`), braces (`{}`), commas (`,`) and the word `from` should not be used in any other way.

 Consider this snippet from a comma-separated flat-text Philadelphia, Pennsylvania voting results file:
```
WARD,DIVISION,VOTE TYPE,CATEGORY,SELECTION,PARTY,VOTE COUNT
01,01,A,JUDGE OF THE SUPERIOR COURT,AMANDA GREEN-HAWKINS,DEMOCRATIC,2
01,01,M,JUDGE OF THE SUPERIOR COURT,AMANDA GREEN-HAWKINS,DEMOCRATIC,146
01,01,P,JUDGE OF THE SUPERIOR COURT,AMANDA GREEN-HAWKINS,DEMOCRATIC,0
```
The formula `ReportingUnit=Ward <WARD>;Division <DIVISION>` would yield 'Ward 01;Division 01'.

### the `count_location` parameter for xml and json files
In tree-structured file types like json and xml files, the `count_location` parameter is a path to the count.

Consider this snippet from a Georgia voting results xml file:
 ```
 <ElectionResult>
    <Contest>
        <Choice key="25" text="Raphael Warnock (Dem)" totalVotes="2279559">
            <VoteType name="Election Day Votes" votes="485573">
                <County name="Appling" votes="391" />
                <County name="Atkinson" votes="262" />
                <County name="Bacon" votes="162" />
```

Because the number of votes is an attribute, we set
```count_location=ElectionResult/Contest/Choice/VoteType/County.votes```
The munge formula `<VoteType.name>` would yield 'Election Day Votes'. It is important that the first element in any munge formula (in this example, `VoteType`) be in the `count_location` path. 


Consider this snippet from a `json` file from Virginia (with line breaks added for clarity):
 ```
{
"ElectionName":"2020 November General",
"ElectionDate":"2020-11-03T00:00:00",
"CreateDate":"2021-01-21T15:03:04.979827-05:00",
"RaceName":"Member House of Representatives (03)",
"NumberOfSeats":1,
"Localities":[
    {
    "Locality":{"LocalityName":"CHESAPEAKE CITY","LocalityCode":"550"},
    "PrecinctsReporting":31,"PrecinctsParticipating":31,
    "LastModified":"2020-11-09T18:12:13.92",
    "Candidates":[
        {
        "BallotName":"Robert C. \"Bobby\" Scott","BallotOrder":1,
        "Votes":35042,"Percentage":"63.70%",
        "PoliticalParty":"Democratic"
        },
 ...
```
Each pair of braces ({}) and its contents is a json "object", while each pair of brackets ([]) and its contents is a json "array". The `count_locations` path indicates which arrays one must step into to find the count, as well as the label for the count itself. For the Virginia example above: 
`count_location=Localities/Candidate/Votes`. Munge formulas referencing the top level are simply the label for the information (e.g., `CandidateContest=<RaceName>`). Munge formulas referencing information within an array must start with that array name, and provide the path within that array to the desired value (e.g., `ReportingUnit=<Localities.Locality.LocalityName>`).

Note that `count_location` is a `/`-separated path, possibly with a period (`.`) at the end for xml files where the count is in an attribute. Munge formulas use only periods.


### flat text 
Consider this snippet from a tab-separated North Carolina voting results file:
```
County	Election Date	Precinct	Contest Group ID	Contest Type	Contest Name	Choice	Choice Party	Vote For	Election Day	One Stop	Absentee by Mail	Provisional	Total Votes	Real Precinct
ALAMANCE	11/06/2018	064	1228	S	NC COURT OF APPEALS JUDGE SEAT 3	Michael Monaco, Sr.	LIB	1	59	65	2	1	127	Y
ALAMANCE	11/06/2018	03N	1228	S	NC COURT OF APPEALS JUDGE SEAT 3	Michael Monaco, Sr.	LIB	1	59	38	1	0	98	Y
ALAMANCE	11/06/2018	03S	1228	S	NC COURT OF APPEALS JUDGE SEAT 3	Michael Monaco, Sr.	LIB	1	106	108	0	3	217	Y
```
Here the CountItemType value ('Election Day','One Stop' a.k.a. early voting, 'Absentee by Mail','Provisional') must be read from the column headers, i.e., the information in row 0 of the file. For the first data row, the formula `CountItemType=<count_header_0>` would yield CountItemType 'Election Day' for the VoteCount of 59, 'One Stop' for the vote count of 65, etc.

### Multiple Sheets or Blocks of Data
To be read automatically, information that is constant over a sheet (or block) must be read either from the sheet name (using `<sheet_name>`) or from the left-most, non-blank entry in a row of the sheet using `<row_j>`, where `j` is the row number. Row numbers start with `0` after skipping the number of rows given in `rows_to_skip`.

NB: the system assumes that blocks are separated by blank lines. This means that files with blank lines internal to the blocks must be revised before processing.


### Regular Expressions
Sometimes it is necessary to use regular expressions to extract information from fields in the results file.  For example, in a primary both the Party and the Candidate may be in the same string (e.g., "Robin Perez (DEM)"). Braces ({}) indicate that regex analysis is needed. Inside the curly brackets there are two parts, separated by a comma. The first part is the field name to pull a string from the file. The second is a python regex formula whose first group (enclosed by parentheses) marks the desired substring.
```Candidate	{<count_header_1>,^(.*) \([a-zA-Z]{3}\)$}	row```

The system will report (in the `.warnings` files) any strings that did not match the regex. 

For an interactive sandbox for learning python regex, try `regex101.com`.

### \[\<foreign key\> lookup\]

If any of the munge formulas depend on information from other files, munger must specify lookup information. For each foreign key, there must be a separate section with corresponding header (foreign key, plus " lookup", e.g. `[CANDIDATE NAME lookup]`) if the results file has a `CANDIDATE NAME` field. This section needs:
  * `lookup_id` is the single field name holding the key to the lookup table. (If there are no headers in the lookup source file, use, e.g., `column_0`)
  * `source_file` the path to the source file, relative to the results directory given in `run_time.ini`
  * all the usual format parameters except `count_location`

For example, here is a munger for Texas 2020 General election results using lookups.
```
[format]
file_type=excel
count_location=by_name:TOTAL VOTES PER OFFICE PER COUNTY

encoding=utf_8
thousands_separator=,
count_column_numbers=4
noncount_header_row=0
count_field_name_row=0

[munge formulas]
ReportingUnit=<COUNTY NAME>
Party=<PARTY from CANDIDATE NAME>
CandidateContest=<OFFICE NAME>
Candidate=<CANDIDATE NAME>
CountItemType=total

[CANDIDATE NAME lookup]
source_file=Texas/Party_by_Candidate_20g.xlsx
lookup_id=CANDIDATE NAME
file_type=excel
noncount_header_row=0
```
NB: if there are multiple rows in the lookup file with the same values for the lookup id columns, the system will arbitrarily use the first row and ignore the others.

### \[ignore\]
 Unrecognized Contests, Candidates and Parties are collected as "none or unknown". Some states (e.g., Wisconsin 2018 General) report total votes over a contest next to individual candidates' votes. The system may read, e.g., "Total Votes Cast" as an unrecognized party name. In this case include the lines:
  ```
[ignore]
Party=Total Votes Cast
```
and similarly, if necessary, for any Contest or Selection. If there is more than one Candidate (e.g.) to be ignored, use a comma-separated list: `Candidate=Total Votes Cast,Registered Voters`

 You may find it helpful to follow the example of the mungers in the repository.

## Jurisdiction Files and `dictionary.txt`
Because each original raw results file comes from a particular election agency, and each election agency has a fixed jurisdiction, we organize information by jurisdiction. The  [`000_for_all_jurisdictions` folder](../src/jurisdictions/000_for_all_jurisdictions) holds information pertinent to all jurisdictions: the list of elections in [`Election.txt](../src/jurisdictions/000_for_all_jurisdictions/Election.txt)

It's easiest to use the JurisdictionPrepper() object to create or update jurisdiction files.

 (0) Create a `jurisdiction_prep.ini` file, following the example in `src/parameter_file_templates/jurisdiction_prep.ini.template`. You will need to specify the number of congressional, state house and state senate districts.

 (1) From the directory containing `jurisdiction_prep.ini`, open a python interpreter. Import the package and initialize a JurisdictionPrepper(), e.g.:
```
>>> import electiondata as ed
>>> jp = ed.JurisdictionPrepper()
```
 (2) Call new_juris_files(), which will create the necessary files in the jurisdiction directory, as well as a starter dictionary file (`XX_starter_dictionary.txt`) in the current directory.
```
>>> jp.new_juris_files()
```
The routine `new_juris_files` creates the necessary files in a folder (at the location `jurisdiction_path` specified in `jurisdiction_prep.ini`). Several of these files are seeded with information that can be deduced from the other information in `jurisdiction_prep.ini`.
 
In addition, `new_juris_files` creates a starter dictionary `XX_starter_dictionary.txt` in your current directory. Eventually the `dictionary.txt` file in your jurisdiction directory will need to contain all the mappings necessary for the system to match the data read from the results file ("raw_identifiers") with the internal database names specified in the other `.txt` files in the jurisdiction directory. The starter dictionary maps the internal database names to themselves, which is usually not helpful. In the steps below, you will correct (or add) lines to `dictionary.txt` following the conventions in the file. The system does not try to guess how internal database names are related to names in the files. 

NB: it is perfectly OK to have more than one raw_identifier_value for a single element. This can be necessary if, say, different counties use different names for a single contest. What can cause problems are lines with the same cdf_element and same raw_identifier_value, but different cdf_internal_names.

 (3) Add all counties to the `ReportingUnit.txt` file and `XX_starter_dictionary.txt`. 
 
The `ReportingUnit.txt` file determines the names used internally in your database. (We have found that including "County" in the name helps with disambiguation.) You must obey the semicolon convention so that the system will know that the counties are subunits of the jurisdiction. For example:
```
Name	ReportingUnitType
Texas;Angelina County	county
Texas;Gregg County	county
Texas;Harrison County	county
```

NB: in some jurisdictions, the major subdivision type is not 'county'. For instance, Louisiana's major subdivisions are called 'parish'. In the `elections.analyze` module, several routines roll up results to the major subdivision -- usually counties. By default, the ReportingUnitType of the major subdivision is read from the file [major_subjurisdiction_types.txt](../src/jurisdictions/000_for_all_jurisdictions/major_subjurisdiction_types.txt) if possible; if that file is missing, or does not provide a subdivision type for the particular jurisdiction in question, the system will try to deduce the major subdivision type from the database. A different file of subdivision types can be specified with the optional `major_subdivision_file` parameter in `Analyzer()` or `DataLoader()`

The system assumes that internal database names of ReportingUnits carry information about the nesting of the basic ReportingUnits (e.g., counties, towns, wards, etc., but not congressional districts) via semicolons. For example: `
 * `Pennsylvania;Philadelphia;Ward 8;Division 6` is a precinct in 
 * `Pennsylvania;Philadelphia;Ward 8`, which is a ward in
 * `Pennsylvania;Philadelphia`, which is a county in
 * `Pennsylvania`, which is a state.
 
Other nesting relationships (e.g., `Pennsylvania;Philadelphia;Ward 8;Division 6` is in `Pennsylvania;PA Senate District 1`) are not yet recorded in the system (as of 4/2/2021).

To add the necessary lines to the starter dictionary `XX_dictionary.txt` you will need to map the names used in the raw results file(s) to the internal names used in `ReportingUnit.txt`. We call the names used in the raw results file 'raw identifiers'. Look in your results files to see how counties are written. For example, if your results file looks like this (example from Texas):
```
ELECTION DATE-NAME	OFFICE NAME	CANDIDATE NAME	COUNTY NAME	TOTAL VOTES PER OFFICE PER COUNTY
03/03/2020 - 2020 MARCH 3RD REPUBLICAN PRIMARY	U. S. REPRESENTATIVE DISTRICT 1	JOHNATHAN KYLE DAVIDSON	ANGELINA	1,660
03/03/2020 - 2020 MARCH 3RD REPUBLICAN PRIMARY	U. S. REPRESENTATIVE DISTRICT 1	LOUIE GOHMERT	ANGELINA	10,968
03/03/2020 - 2020 MARCH 3RD REPUBLICAN PRIMARY	U. S. REPRESENTATIVE DISTRICT 1	JOHNATHAN KYLE DAVIDSON	GREGG	914
03/03/2020 - 2020 MARCH 3RD REPUBLICAN PRIMARY	U. S. REPRESENTATIVE DISTRICT 1	LOUIE GOHMERT	GREGG	9,944
03/03/2020 - 2020 MARCH 3RD REPUBLICAN PRIMARY	U. S. REPRESENTATIVE DISTRICT 1	JOHNATHAN KYLE DAVIDSON	HARRISON	774
03/03/2020 - 2020 MARCH 3RD REPUBLICAN PRIMARY	U. S. REPRESENTATIVE DISTRICT 1	LOUIE GOHMERT	HARRISON	7,449
```
you would want lines in your dictionary file like this:
```
cdf_element	cdf_internal_name	raw_identifier_value
ReportingUnit	Texas;Angelina County	ANGELINA
ReportingUnit	Texas;Gregg County	GREGG
ReportingUnit	Texas;Harrison County	HARRISON
```
Note that the entries in the `cdf_internal_name` column exactly match the entries in the `Name` column in `ReportingUnit.txt`.

 (4) As necessary, revise `CandidateContest.txt` (along with `Office.txt` and `XX_starter_dictionary.txt`). 
 * The offices and candidate-contests added by `new_juris_files()` are quite generic. For instance, your jurisdiction may have a 'Chief Financial Officer' rather than an 'Treasurer'. Use the jurisdiction's official titles, from an official government source. Add any missing offices. Our convention is to preface state, district or territory offices with the two-letter postal abbreviation. For example (in `Office.txt`):
```
Name	ElectionDistrict
US President (FL)	Florida
FL Governor	Florida
US Senate FL	Florida
FL Attorney General	Florida
FL Chief Financial Officer	Florida
FL Commissioner of Agriculture	Florida
```
If you are interested in local contests offices (such as County Commissioner), you will need to add them. If the ElectionDistrict for any added contest is not already in `ReportingUnit.txt`, you will need to add it. For judicial retentions and other Note that judicial retention elections are yes/no, so they should be handled as BallotMeasureContests, not CandidateContests. NB: If you want to add Offices in bulk from a results file, you can wait and do it more easily following instructions below.

For each new or revised Office, add or revise entries in `CandidateContest.txt`. Leave the PrimaryParty column empty. Do not add primaries at this point -- they can be added in bulk below.  For example (in `CandidateContest.txt`):
 ```
US President (FL)	1	US President (FL)	
FL Governor	1	FL Governor		
US Senate FL	1	US Senate FL	
FL Attorney General	1	FL Attorney General	
FL Chief Financial Officer	1	FL Chief Financial Officer	
FL Commissioner of Agriculture	1	FL Commissioner of Agriculture	
```

Finally, look in your results files to see what naming conventions are used for candidate contests. Add lines to the starter dictionary. For example, using data from official Florida election results files:
```
cdf_element	cdf_internal_name	raw_identifier_value
CandidateContest	US President (FL)	President of the United States
CandidateContest	US House FL District 1	Representative in Congress District 1
CandidateContest	US House FL District 2	Representative in Congress District 2
```

 (5) Make any necessary additions or changes to the more straightforward elements (e.g., `Candidate.txt`). It's often easier to add these in bulk later directly from the results files (see 'Error reporting' below) -- unless you want to use internal names that differ from the names in the results file.
  * `Party.txt`. You may be able to find a list of officially recognized parties on the Board of Election's website.
  * `BallotMeasure.txt`. A BallotMeasure is any yes/no question on the ballot, including judicial retention. Each BallotMeasure must have an ElectionDistrict and an Election matching an entry in the `ReportingUnit.txt` or `Election.txt` file.

 (6) Revise `XX_starter_dictionary.txt` so that it has entries for any of the items created in the steps above (except that there is no need to add Elections to the dictionary, as they are never munged from the contents of the results file). The 'cdf_internal_name' column should match the names in the jurisdiction files. The 'raw_identifier_value' column should hold the corresponding names that will be created from the results file via the munger. 
    * It is helpful to edit the starter dictionary in an application where you can use formulas, or to manipulate the file with regular expression replacement. If you are not fluent in manipulating text some other way, you may want to use Excel and its various text manipulation formulas (such as =CONCAT()). However, beware of Excel's tendency to revise formats on the sly. You may want to check `.txt` and `.csv` files manipulated by Excel in a plain text editor if you run into problems. (If you've been curious to learn regex replacement, now's a good time!)
 
 (7) Revise the 'CountItemType' entries in the starter dictionary to match any words or phrases used in the results files. E.g., for North Carolina
```
cdf_element	cdf_internal_name	raw_identifier_value
CountItemType	election-day	Election Day
CountItemType	early	One Stop
CountItemType	absentee-mail	Absentee by Mail
CountItemType	provisional	Provisional
CountItemType	total	Total Votes
CountItemType	total	total
```
Note that, unlike for Candidates, Parties, etc., the internal database names of the CountItemTypes (a.k.a. vote types) are specified here, in the dictionary file. The starter dictionary has internal names from the NIST Common Data Format. Using other internal names won't break anything, but will throw a warning.

 (8) Add any existing content from `dictionary.txt` to the starter dictionary. If the jurisdiction is brand new there won't be any existing content. 

 (9) Move `XX_starter_dictionary.txt` from the current directory and to the jurisdiction's directory, and rename it to `dictionary.txt`. 

 (10) If your results file is precinct based instead of county based, and you would like to use the results file to add precincts to your `ReportingUnit.txt` file, run `add_sub_county_rus`, e.g.: 
```
>>> jp.add_sub_county_rus('my_results.ini')
```
These will be added as precincts, unless another reporting unit type is specified with the optional argument `sub_ru_type`, e.g.:
```
>>> jp.add_sub_county_rus_from_multi_results_file('my_results.ini',sub_ru_type='town')
```
If the jurisdiction's major subdivision is not county but something else (e.g., state house district, as in Alaska, or parish as in Louisiana), use the optional argument `county_type`:
```
>>> jp.add_sub_county_rus_from_multi_results_file('my_results.ini',county_type='state-house')
```
Finally, Look at the newly added items in `ReportingUnit.txt` and `dictionary.txt`, and remove or revise as appropriate.

### Miscellaneous Notes
  * Candidate
      * Look for possible variant names (e.g., 'Fred S. Martin' and 'Fred Martin' for the same candidate in two different counties. If you find variations, pick an internal database name and put a line for each raw_identifier_value variation into `dictionary.txt`.
      * Look for non-candidate items treated like candidates in the raw file. For example, a file may have "Undervotes" in the candidate column. To exclude such lines from your results, exclude them from the `Candidate.txt` file. You can also exclude themby adding a line, e.g., `Candidate=Undervotes`,  to the `[ignore]` section of the munger. 
      * Our convention for internal names for candidates with quoted nicknames is to use single quotes. Make sure there are no double-quotes in the Name column in `Candidate.txt` and none in the cdf_internal_name column of `dictionary.txt`. E.g., use `Rosa Maria 'Rosy' Palomino`, not `Rosa Maria "Rosy" Palomino`. Note that if your results file has `Rosa Maria "Rosy" Palomino`, the system will read the double-quotes as single-quotes, so the dictionary.txt line `Candidate	Rosa Maria 'Rosy' Palomino	Rosa Maria "Rosy" Palomino` will munge both `Rosa Maria 'Rosy' Palomino` and `Rosa Maria "Rosy" Palomino` correctly.
     * Our convention for internal names for multiple-candidate tickets (e.g., 'Biden/Harris' is to use the full name of the top candidate, e.g., 'Joseph R. Biden'). There should be a line in `dictionary.txt` for each variation used in the results files, e.g.:
```
cdf_element	cdf_internal_name	raw_identifier_value
Candidate	Joseph R. Biden	Biden / Harris
Candidate	Joseph R. Biden	Joseph R. Biden
```

  * CandidateContest: For any new CandidateContest you do want to keep you will need to add the corresponding line to `Office.txt` (and the ElectionDistrict to `ReportingUnit.txt` if it is not already there). 
  * Primary Elections: if you will be munging primary elections, we recommend making each primary a separate election (e.g., "2021 Democratic Primary", "2021 Republican Primary")

#### Required Conventions
For ReportingUnits, the naming convention is to list as much of the composing information as possible in the name of the element, using `;` as a separator. E.g., 
 * `North Carolina` -- the state of NC
 * `North Carolina;Alamance County` -- Alamance County, which is contained in North Carolina
 * `North Carolina;Alamance County;Precinct 12W` -- Precinct 12W in Alamance County
The semicolons are used by the code to roll up results from smaller Reporting Units into larger Reporting Units.
   
Replace any double-quotes in Candidate.txt and dictionary.txt with single quotes. I.e., `Rosa Maria 'Rosy' Palomino`, not `Rosa Maria "Rosy" Palomino`.


#### Optional Conventions
The jurisdiction files in this repository follow certain conventions. Many of these are optional; using different conventions in another copy of the system will not break anything. Internal database names names are standardized as much as possible, regardless of state, following these models:

    * Office
        * `US Senate CO`
        * `US House FL District 5`
        * `PA State Senate District 1`
        * `NC State House District 22`
        * `PA Berks County Commissioner`
        * `DC City Council District 2`
        * `OR Portland City Commissioner Seat 3`
       
    * Party
        * `Constitution Party` (regardless of state or other jurisdiction)
        
#### Beware of:
 - hyphens in formal names of jurisdictions or elections -- this may break the testing.
 - Different names for same contest in different counties (if munging from a batch of county-level files)
 - Different names for candidates, especially candidates with name suffixes or middle/maiden names
 - Different "party" names for candidates without a party affiliation 
 - Any item with an internal comma (e.g., 'John Sawyer, III')
 - A county that uses all caps (e.g., Seminole County FL)
 - % signs in parameter files, particularly as web addresses for results_source (e.g.,https://elections.wi.gov/sites/elections.wi.gov/files/2020-11/UNOFFICIAL%20WI%20Election%20Results%202020%20by%20County%2011-5-2020.xlsx) -- system cannot read ini files with % signs
 - Line breaks in parameter files may interfere with software parsing the following lines
 - Files that total over all candidates (e.g., Nebraska 2020 General). Make sure not to include the totals in the counts as a nominal "candidate".
 - Excel files that show a thousands separator when you view them (2,931) but don't have a thousands separator under the hood (2931). If all your count are zero, try adding or removing the 'thousands-separator' parameter in `format.config`.
 - the parser for multi-block flat files assumes that any line without counts separates blocks. Beware of stray count-less lines (e.g., can be created by utilities pulling tables from pdfs).


If your sheets or files have a variable number of count columns (e.g., if columns are labeled by candidates), err on the side of including extra columns in count_column_numbers. Columns without data will be ignored. Be careful, however, not to include in your count columns any columns containing strings needed for munging.

If your excel file has merged cells across lines, it may not be clear which line holds the information. Save a sheet as tab-separated text to see which line holds which information from merged cells.

If not all rows are data, and some string fields to be munged have blank headers (e.g., often the counties are in the first column without a cell above reading "County"). In this case use '<column_i>' in the munge formulas to denote the i-th column (numbering from 0, as usual). For example, if counties are in the leftmost column and the header is blank, use '<column_0>' for the name of the county. See for example `wy_gen.munger`.

If there are hidden columns in an Excel file, you may need to omit the hidden columns from various counts.

# Loading Data
Each results file to be loaded must be designated in a `*.ini` file inside its jurisdiction's corresponding subfolder of `ini_files_for_results` in the repository. The `*.ini` files currrently in this repository correspond to [official raw data files for the US 2020 General Election](https://doi.org/10.7910/DVN/0GKBGT). These should load directly with the munger and jurisdiction files from the `electiondata` repository. (Note, however, that due to Excel corruption issues, Vermont and Wisconsin files may fail to load; Connecticut, Maryland and Pennsylvania will load but may fail some of the tests because of inconsistencies within their official agencies' materials.)
 
The DataLoader class allows batch uploading of all data in the directory indicated by the `results_dir` parameter in the main parameter file. The subdirectories of this file should be named for the jurisdictions (with hyphens replacing spaces, as in 'US-Virgin-Islands'. The `DataLoader.load_all()` method will upload every result file that appears, as long as its path (relative to the `results_dir`) is the `results_file` parameter for some `*.ini` file in `ini_files_for_results`. 

```
import electiondata as ed
ed.load_or_reload_all()
```

Some results files may need to be munged with multiple mungers, e.g., if they have combined absentee results by county with election-day results by precinct. If the `.ini` file for that results file has `munger_list` set to a comma-separated list of mungers, then all those mungers will be run on that one file.

## Error reporting
It is difficult to get the munger, jurisdiction and ini files all correct on the first try. The package has robust error- and warning-reporting, designed to help the user make necessary corrections. Often the most efficient way to proceed is to examine all the .error and .warning files and use their contents to modify the munger, ini or jurisdiction files.

All errors and warnings will be reported to the a subdirectory named by the database and timestamp within the directory specified by the `reports_and_plots_dir` parameter in the main parameter file.

Even when the upload has worked, there may be warnings about lines not loaded. The system will ignore lines that cannot be munged. For example, the only contests whose results are uploaded will be those in the `CandidateContest.txt` or `BallotMeasureContest.txt` files that are correctly described in `dictionary.txt`.

If there are no errors, the results files will be moved to a subdirectory of the directory specified by the `archive_dir` parameter in the main parameter file.

## Unloading and Reloading with `reload_juris_election()`
To unload existing data for a given jurisdiction and a given election you can use the routine 
```ed.reload_juris_election(jurisdiction, election, report_dir)```
This function takes optional arguments:
 * rollup (defaults to  false), if true, rolls up results within the to the major subdivision
 * dbname, name of database if given; otherwise database name taken from parameter file
 * param_file, path to parameter file for dataloading if given; otherwise parameter file
            path assumed to be 'run_time.ini'
 * move_files (defaults to True), if true, move all files to archive directory if loading (& testing, if done)
            are successful
 * run_tests (default to True), if true, run tests on data to be loaded, and load only if the tests are passed.

Results of the test will be reported in a the directory specified by the `reports_and_plots_dir` parameter in the parameter file.

## Loading Results from a Single File with Multiple Elections or Jurisdictions
Sometimes it is useful to load a single file with results from several elections or jurisdictions. For example, a secondary source may have combined results information into one file. The method `DataLoader.load_multielection_from_ini()` method allows this kind of upload. This method requires an initialization file, with all the usual required parameters except `election` and `jurisdiction`, and with the additional parameter `secondary_source`. The value of `secondary-source` should be the name of a subfolder of `src/secondary_sources` in the repository, containing files listing the elections and jurisdictions. E.g., 
```
[election_results]
results_file=MEDSL/county_2018.csv
munger_list=medsl_2018
secondary_source=MIT-Election-Data-Science-Lab
results_short_name=medsl_2018_county
results_download_date=2021-07-10
results_source=https://github.com/MEDSL/2018-elections-official/blob/master/county_2018.csv
results_note=county-level
```
The method `DataLoader.load_multielection_from_ini()` takes two optional parameters: 
 * `overwrite_existing` (default `False`): if True, will delete from database any existing results for each election-jurisdiction pair represented in the results file
 * `load_jurisdictions` (default `False`): if True, will load or update database with the jurisdiction information (from `src/jurisdictions`) for each jurisdiction represented in the results file


# Pulling Information from the Database
## The Anaylzer Class
The Analyzer class is used to export data and plots from the database of election results.
The Analyzer class takes two optional parameters:
 * `param_file`: path to the parameter file for defining the Analyzer directories, database connection, etc. if not specified, the default is the `run_time.ini` file in the directory from which the call to Analyzer() is made
 * `dbname`: name of database to use. If not specified, the default is the database specified in the `param_file`

To get an instance of an analyzer, you can call the Analyzer class directly:
```
import electiondata as ed
analyzer = ed.Analyzer(param_file=param_file, dbname=dbname)
```
or, since every instance of the DataLoader class creates its own analyzer:
```
import electiondata as ed
dl = ed.DataLoader(param_file=param_file, dbname=dbname)
analyzer = dl.analyzer
```

## Tabular Export
The Analyzer class has a number of functions that allow you to aggregate the data for analysis purposes. For example, to export all 2020 General results in your database to a tab-separated file `tabular_results.tsv`:
```
analyzer.export_election_to_tsv("tabular_results.tsv", "2020 General")
```
To export results for a single jurisdiction, use, e.g.:
```
analyzer.export_election_to_tsv("tabular_results.tsv", "2020 General", "South Carolina")
```

This code will produce all South Carolina data from the 2018 general election, grouped by contest, county, and vote type (total, early, absentee, etc).

## NIST Common Data Format Export
This package provides functionality to export the data to xml or json according to the [NIST election results reporting schema (Version 2)](https://github.com/usnistgov/ElectionResultsReporting/raw/version2/NIST_V2_election_results_reporting.xsd). 

This is as simple as identifying an election and jurisdiction of interest. For xml:
```
import electiondata as ed
analyzer = ed.Analyzer()
election_report = analyzer.export_nist_xml_as_string("2020 General", "Georgia")
```
The output is a string, the contents of the xml file.

And for json:
```
analyzer = ed.Analyzer()
analyzer.export_nist_json_as_string("2020 General","Georgia")
```
The output is a string, the contents of the json file.
The subdivision type for the roll-up is determined by the [`000_major_subjurisdiction_type.txt file](../src/jurisdictions/000_for_all_jurisdictions/major_subjurisdiction_types.txt).


## Difference-in-Difference Calculations
The system provides a way to calculate difference-in-difference statistics. For any particular election, `Analyzer.diff_in_diff_dem_vs_rep` produces a dataframe of values for any county with results by vote type, with Democratic or Republican candidates, and any comparable pair of contests both on some ballots in the county. Contests are considered "comparable" if their districts are of the same geographical district type -- e.g., both statewide, or both state-house, etc. The method also returns a list of jurisdictions for which vote counts were zero or missing.
```
dbname = "test_0314_1836"
election = "2020 General"
an = eda.Analyzer(dbname=dbname)
diff_in_diff_dem_vs_rep, missing = an.diff_in_diff_dem_vs_rep(election)
```

Specifically, for a fixed county and party, for a fixed pair of vote types and for a fixed pair of contests, we calculate the difference-in-difference value to be
```abs(abs(pct[0][0] - pct[1][0]) - abs(pct[0][1] - pct[1][1]))```
where `pct[i][j]` denotes the percentage of the total vote share earned by the party's candidate in contest `i` on ballots in the county of vote type `j`. The vote share is of the votes for all candidates, not just Democratic or Republican. However, we omit contests that don't have both Republican and Democratic candidates

For more information and context about difference-in-difference calculations for election results, see Michael C. Herron's article [Mail-In Absentee Ballot Anomalies in North Carolina's 9th Congressional District](http://doi.org/10.1089/elj.2019.0544). Note that he uses signed difference-in-difference, while we take the absolute value.

## Plots
See [Sample_Session](./Sample_Session.md).# Results files for dataloading tests
The dataloading tests rely on having some raw results data to load. And the results data should be various enough to test the various components of the data-loading code. In other words, effective testing requires a reasonable variety of input files. The repository does not contain sufficient results data for testing. A test set is available in a separate repository, [TestingData](https://github.com/ElectionDataAnalysis/TestingData). If [test_dataloading_by_ej.py](../tests/dataloading_tests/test_dataloading_by_ej.py) does not find results data, it will default to downloading the files from that repository.

# Sample Testing Session

## Directory and File Structure
 Call the tests from a working directory with the following structure and files:
```
.
+-- input_results
|   +-- Alabama
|   |   + <Alabama results file>
|   |   + <maybe another Alabama results file>
|   +-- Alaska
|   |   + <Alaska results file>
|   +-- American-Samoa
|   |   + <American Samoa results file>
|   +-- <etc>
|
+-- reports_and_plots
+-- run_time.ini
```
The file `run_time.ini` can be the same as in the [Sample Dataloading Session](Sample_Session.md). 

## Note on dataloading tests
The tests in [test_dataloading_by_ej.py](../tests/dataloading_tests/test_dataloading_by_ej.py) will attempt to load all raw results files in `input_results` that are specified by some file in the [`ini_file_for_results` directory](../src/ini_files_for_results). You can check which jurisdictions had files loaded:
 * if the test is successful, look at the `compare_*` directories in the `reports_and_plots` directory.
 * if the test fails, look in the output from the test.

## Running the tests
You will need pytest to be installed on your system (see [pytest installation instructions](https://docs.pytest.org/en/6.2.x/getting-started.html) if necessary). Commands are run from the shell, referencing the local path to the repository
 * dataloading routines: `pytest path/to/repo/tests/dataloading_tests`
 * jurisdiction prep routines: `pytest path/to/repo/tests/jurisdiction_prepper_tests/`
 * analysis routines: `pytest path/to/repo/tests/analyzer_tests/  `# Sample Dataloading Session

This document walks the reader through a simple example, from setting up project directories, through loading data and performing analyses. We assume that the package has been installed in an environment with all the necessary components (as described in [Installation.md](Installation.md)). As an example, we will load the xml results file from Georgia in the repository at [tests/000_data_for_pytest/2020-General/Georgia/GA_detail_20201120_1237.xml](../tests/000_data_for_pytest/2020-General/Georgia/GA_detail_20201120_1237.xml).

## Directory and File Structure
The package offers a fair amount of flexibility in the directory structures used. For this sample session, we assume the user will call the program from a working directory with the following structure and files:
```
.
+-- input_directory
|   +-- Georgia
|   |   +-- GA_detail_20201120_1237.xml
+-- archive_directory
+-- reports_and_plots_directory
+-- run_time.ini
```
    
Note that during processing the package uses information from the repository. In other words, the repository contains not only the code necessary to compile the package, but also files called by the package as it functions -- files with information about jurisdictions, mungers, results and result files. So the user will need to know the absolute path to the repository content root `src`. Below we will call this path `<path/to/src>`.

### Contents of `run_time.ini`
```
[electiondata]
results_dir=input_results
archive_dir=archive_directory
repository_content_root=<path/to/src>
reports_and_plots_dir=reports_and_plots

[postgresql]
host=localhost
port=5432
dbname=ga_test
user=postgres
password=
```
You may wish to check that these postgresql credentials will work on your system via the command `psql -h localhost -p 5432 -U postgres postgres`. If this command fails, or if it prompts you for a password, you will need to find the correct connection parameters specific to your postgresql instance.  (Note that the `dbname` parameter is arbitrary, and determines only the name of the postgresql database created by the package.)

### Contents of `GA_detail_20201120_1237.xml`
Copy the file of the same name in the repository: [000_template.mungerGA_detail_20201120_1237.xml](../tests/000_data_for_pytest/2020-General/Georgia/GA_detail_20201120_1237.xml)

## Load Data
```
>>> import electiondata as ed
>>> ed.load_or_reload_all()
```
After this command executes successfully, you will see a database in your postgres instance whose name matches the value of `dbname` given in `run_time.ini`. The file structure will now be:
```
.
+-- archive_directory
|   +-- Georgia_2020-11-10
|   |   +-- GA_detail_20201120_1237.xml
+-- input_directory
+-- reports_and_plots_directory
|   +-- compare_to_Georgia_xxxx_xxxx
|   |   +-- parameters.txt
|   |   +-- not_found_in_db.tsv
|   |   +-- ok.tsv
|   |   +-- wrong.tsv
|   +-- load_or_reload_all_xxxx_xxxx
|   |   +-- Georgia_jurisdiction_dictionary.txt.warnings
+-- run_time.ini
```
 Note that the `Georgia` results folder has been moved to the archive directory, with a date stamp from the date of the results file. The `reports_and_plots_directory` contains:
 * a subdirectory `compare_to_Georgia_xxxx_xxxx` with the results of the comparison to the [reference results](../src/reference_results/Georgia.tsv)
 * a subdirectory `load_or_reload_all_xxxx_xxxx` with warnings from the data uploading. You may wish to look at the file `Georgia_jurisdiction_dictionary.txt.warnings`, which lists the contests present in the xml results file that were not recognized during processing. If these contests and their candidates had been added to the Georgia-specific information in the repository (see "Creating Jurisdiction files" below, or in the [User Guide](User_Guide.md)), these warnings would not appear.

## Export and analyze data
To pull data out, you will need to use the Analyzer class:
```
>>> an = ed.Analyzer()
```

### Export
You can export results in tabular form:
```
>>> an.export_election_to_tsv("GA_results.tsv", "2020 General")
```
The  file `GA_results.tsv` containing the tab-separated results will be created in the working directory:
```
.
+-- archive_directory
|   +-- Georgia_2020-11-10
|   |   +-- GA_detail_20201120_1237.xml
+-- GA_results.tsv
+-- input_directory
+-- reports_and_plots_directory
|   +-- compare_to_Georgia_xxxx_xxxx
|   |   +-- parameters.txt
|   |   +-- not_found_in_db.tsv
|   |   +-- ok.tsv
|   |   +-- wrong.tsv
|   +-- load_or_reload_all_xxxx_xxxx
|   |   +-- Georgia_jurisdiction_dictionary.txt.warnings
+-- run_time.ini
```

The program (v.2.0.1 and higher) can also produce a string of data in the NIST Common Data Format Version 2.0, in either json or xml format:
```
>>> an.export_nist_xml_as_string("2020 General", "Georgia")
>>> an.export_nist_json_as_string("2020 General", "Georgia")
```

### Plots
To draw pictures automatically, you will need [`orca` installed on your system](https://github.com/plotly/orca). Plots will be exported to the `reports_and_plots_directory` and may also appear in a browser window.

If `orca` is not installed, the system will not automatically create pictures. Note that in any case the routines return text information sufficient to create plots in any system you may wish to use. 

#### Scatter Plots
You can create scatter plots of results by county. For example, create a jpeg comparing Biden's vote totals to Trump's vote totals with:
```
>>> biden_v_trump = an.scatter("Georgia","2020 General","Candidate total","Joseph R. Biden","2020 General","Candidate total","Donald J. Trump",fig_type="jpeg")
```
![Biden vs. Trump scatter plot of total votes by county for Georgia, 2020 General Election](images/scatter_Joseph-R-Biden-2020-General-US-President-GA-_Donald-J-Trump-2020-General-US-President-GA.jpeg)

Or compare Biden's votes on election day with votes on absentee mail ballots:
```
>>> biden_eday_v_abs = an.scatter("Georgia","2020 General","Candidate election-day","Joseph R. Biden","2020 General","Candidate absentee-mail","Joseph R. Biden",fig_type="jpeg")
```
In each case the value returned is a string with all the information necessary to draw the plot in a string (following the python dictionary format). To return such a string without calling `orca`, simply omit the `fig_type`, e.g., 
```
>>> biden_eday_v_abs = an.scatter("Georgia","2020 General","Candidate election-day","Joseph R. Biden","2020 General","Candidate absentee-mail","Joseph R. Biden")
```

The arguments (in order) are:
* the jurisdiction ("Georgia")
* three items to define the horizontal-axis count:
  * the election ("2020 General")
  * the count category, one of:
    ```Candidate absentee-mail
        Candidate early
        Candidate election-day
        Candidate provisional
        Candidate total
        Contest absentee-mail
        Contest early
        Contest election-day
        Contest provisional
        Contest total
        Party absentee-mail
        Party early
        Party election-day
        Party provisional
        Party total
        ```
  * the specific count within the category
* three items to define the vertical-axis count (same as for horizontal) 

Note that once the category has been chosen, e.g., "Party absentee-mail", the list of possibilities for the specific count can be obtained from this line of code:
```
>>> [entry["name"] for entry in an.display_options("count",["2020 General","Georgia","Party absentee-mail"])]
```
Use any category name in place of "Party absentee-mail" to see counts available for that category.

Categories starting with "Contest" give number of votes tallied in that contest in each county, lumping all candidates together. Categories starting with "Party" give number of votes tallied for members of that party in a particular contest type (e.g., "Libertarian congressional"). 

#### Curated One-County Outlier Bar Charts

The system attempts to find interesting one-county outliers within the election results. The specific algorithm is described in [an article by Singer, Srungavarapu & Tsai in _MAA Focus_, Feb/March 2021, pp. 10-13](http://digitaleditions.walsworthprintgroup.com/publication/?m=7656&i=694516&p=10&ver=html5)

For example:
```
>>> outliers = an.bar("2020 General","Georgia",contest_type="congressional",fig_type="png")
```
This will export up to three bar charts (to `reports_and_plots_directory`) showing a pair of candidates for which the vote shares in one county, for some type of ballot, differs significantly from the vote shares in other counties in the same district. The output of the command above includes: a chart showing how Clarke County differs from other Georgia counties in the 9th US House District. 

![Chart showing Pandy outperforming Clyde in Clarke County GA on early ballots, while results in all other counties favor Clyde](images/Andrew-Clyde-R-_Devin-Pandy-D-_early_US-House-GA-District-9.png)

The output variable `outliers` contains even more information, including the contest margin and an estimate of the votes at stake if the outlier were brought in line with the other counties.

Options for `contest_type` for this data set are: `congressional`, `state` (for statewide contests), `state-house` and `state-senate`. 


#### Difference-in-Difference Analysis
The program offers difference-in-difference analysis where results are available by vote type, following [Herron's analysis of congressional contests](https://www.liebertpub.com/doi/full/10.1089/elj.2019.0544). The following code will create a tab-separated file `GA_diffs.tsv` in the working directory:
```
>>> (did_frame, missing) = an.diff_in_diff_dem_vs_rep("2020 General")
>>> did_frame.to_csv("GA_diff_in_diff.tsv","\t")
```
(Note: the variable `missing` is a list of diff-in-diff comparisons that failed.)

## Optional steps
The sample session above uses the information already in the repository about Georgia and the particular results file. If you wish to create these files yourself from scratch, follow these optional steps.

### Specify contests totals for data quality check
The data loading process includes checking contest totals against reference totals in [src/reference_results/Georgia.tsv](../src/reference_results/Georgia.tsv). You may wish to add some reference totals to this file.

### Create a munger file
If you wish to practice creating  munger file, follow the steps below. 

Your munger file should live in the folder [src/mungers](../src/mungers), and its name needs to have the extension `.munger`. 

1. (Optional) Delete the munger file [ga_xml.munger](../src/mungers/ga_xml.munger) from your local copy of the repository. This step is not strictly necessary, but will help ensure that you don't accidentally use the existing munger.
1. Make a copy of [src/mungers/000_template.munger](../src/mungers/000_template.munger) named `my_Georgia_test.munger`, inside the same folder.  
2. Fill in required parameter values to specify the organization of the results file `000_template.munger`
  * `file_type=xml` indicates that the file is in xml format.
  * `count_location=ElectionResult/Contest/Choice/VoteType/County.votes` indicates where the vote counts are to be found within the xml nesting structure.
  *  Specify the location of the info defining each vote count. Each location must start with one of the nodes in `count_location`
    * Because the file contains a variety of contests, candidates, parties, vote types (a.k.a. CountItemTypes) and geographies (a.k.a. ReportingUnits), there is no need to use the `constant_over_file` parameter.
    * In the `[munge formulas]` section, specify where the other information is found. While the ReportingUnit, CandidateContest and CountItemType can be read simply from quoted strings in the file, (`ReportingUnit` is in `County.name`, `CandidateContest` is in `Contest.text` and `CountItemType` is in `VoteType.name`), the `Candidate` and `Party` must both be read out of `Choice.text` with python's regular expression ("regex") syntax. In the package syntax, the location of the string in the file is enclosed in angle brackets `<>`, and if a regular expression is needed, the location and the regular expression are given as a pair within braces `{}`.
     
```
[format]
file_type=xml
count_location=ElectionResult/Contest/Choice/VoteType/County.votes

[munge formulas]
ReportingUnit=<County.name>
Party={<Choice.text>,^.* \((.*)\)$}
CandidateContest=<Contest.text>
Candidate={<Choice.text>,^(.*) \(.*\)$}
CountItemType=<VoteType.name>
```

### Create an initialization file for the results file (optional)
If you wish to practice creating an initialization file for results, follow the steps below. Otherwise the system will use the initialization file [ga20g_20201120_1237.ini](../src/ini_files_for_results/Georgia/ga20g_20201120_1237.ini).

Because your `input_directory` folder has a subfolder `Georgia`, the dataloading routines will look to the folder [src/ini_files_for_results/Georgia](../src/ini_files_for_results/Georgia) for information about any results files in `input_directory/Georgia` in your working directory. 

1. Delete [ga20g_20201120_1237.ini](../src/ini_files_for_results/Georgia/ga20g_20201120_1237.ini) from the repository.
2. Copy [src/ini_files_for_results/single_election_jurisdiction_template.ini](../src/ini_files_for_results/single_election_jurisdiction_template.ini) to a file with extension `.ini` in the folder [src/ini_files_for_results/Georgia](../src/ini_files_for_results/Georgia).
3. Define the parameters in the `[election_results]` section. 
  * If you did not create a new munger file, use `munger_list=ga_xml` instead of `munger_list=my_Georgia_test`
  * Use the actual download date instead of `2021-11-09`
```
[election_results]
results_file=Georgia/GA_detail_20201120_1237.xml
munger_list=my_Georgia_test
jurisdiction=Georgia
election=2020 General
results_short_name=any_alphanumeric_string
results_download_date=2021-11-09
results_source=electiondata repository
results_note=
is_preliminary=False
```

### Creating jurisdiction files from scratch
If you choose to create your Georgia jurisdiction files from scratch (rather than using the ones provided in the repository) you will need to specify various information about Georgia in a file called `jurisdiction_prep.ini` in your working directory. The working directory will have this structure:
```
.
+-- jurisdiction_prep.ini
+-- run_time.ini
+-- input_directory
|   +-- Georgia
|   |   +-- GA_detail_20201120_1237.xml
+-- archive_directory
+-- reports_and_plots_directory
```
and the contents of `jurisdiction_prep.ini` will be:
```
[electiondata]
name=Georgia
reporting_unit_type=state
abbreviated_name=GA
count_of_state_house_districts=180
count_of_state_senate_districts=56
count_of_us_house_districts=14
```
Because the `input_directory` folder has a subfolder `Georgia`, the dataloading routines will look for information specific to Georgia in the repository folder [src/jurisdictions/Georgia](../src/jurisdictions/Georgia). 

1. Delete the folder [src/jurisdictions/Georgia](../src/jurisdictions/Georgia).
2. Navigate to your working directory.
5. Follow the instructions in the "Create or Improve a Jurisdiction" section of the [User_Guide](User_Guide.md). As you work through them, the following examples may be helpful.
  * Sample county lines in `Reporting.txt`:
```
Georgia;Appling County	county
Georgia;Atkinson County	county
Georgia;Bacon County	county
```
  * Sample county lines in `GA_starter_dictionary.txt` (or `dictionary.txt`):
```
ReportingUnit	Georgia;Appling County	Appling
ReportingUnit	Georgia;Atkinson County	Atkinson
ReportingUnit	Georgia;Bacon County	Bacon
```

4. Try loading your results! 
```
>>> import electiondata as ed
>>> ed.load_or_reload_all(move_files=False)
```
It is useful to prevent the system from archiving files while you are dealing with warnings with the `move_files=False` option.
Errors and warnings are par for the course, so don't be discouraged if the result is something like this:
```
Jurisdiction Georgia did not load. See .error file.
Jurisdiction errors written to reports_and_plots_directory/load_or_reload_all_1113_0757/_jurisdiction_Georgia.errors
>>> 
```
Take a look at the file(s) indicated, which should give you a good idea of what needs to be fixed. (If it doesn't please report the unhelpful message(s) and your question as an [issue on GitHub](https://github.com/ElectionDataAnalysis/electiondata/issues)). Repeat this step -- loading and looking at errors -- as often as necessary! (In rare cases, you may want to start over with a new database, either by erasing the existing `ga_test` from your postgres instance, or changing the value of `dbname` in the `run_time.ini` file.)

### Sample errors and warnings while building jurisdiction files
#### Contests
If no contests were recognized, or no candidates were recognized, the system reports an error:
```
Jurisdiction errors (Georgia):
Could not add totals for 2020 General because no data was found
No contest-selection pairs recognized via munger my_Georgia_test.munger
```
Focus one contest, and one candidate in that contest. Look in the `.errors` and `.warnings` files. If the name of the contest or the candidate appears, the file will tell you what went wrong. If the name of the contest or the candidate does not appear in the `.errors` or `.munger` file, then there is an issue with the munger named in the results initialization file.

Contests that were parsed from the file but not recognized in the dictionary will be listed in a `.warnings` file, e.g.:
```
CandidateContests (found with munger my_Georgia_test) not found in dictionary.txt :
US Senate (Perdue)
US Senate (Loeffler) - Special
Statewide Referendum A
State Senate District 39 - Special Democratic Primary
State Senate Dist 9/Senador Estatal Dist 9
```
You may very well choose to omit certain contests (such as Statewide Referendum A and the Special Democratic Primary). For the other contests, you will need to make an entry in `dictionary.txt` -- and maybe `CandidateContest.txt` and `Office.txt` as well.
    * 'State Senate Dist 9/Senador Estatal Dist 9': if you see this after following the instructions in this sample session, the corresponding Office and CandidateContest 'GA Senate District 9' should already exist, so simply add another line to the `dictionary.txt` file: `CandidateContest	GA Senate District 9	State Senate Dist 9/Senador Estatal Dist 9`. (It may save time to add this version of all the GA Senate districts at this point.)
    * 'US Senate (Perdue)' and 'US Senate (Loeffler) - Special': Be sure to disambiguate these, either by creating two separate Offices and corresponding CandidateContests for the two US Senate positions, or by creating two separate CandidateContests corresponding to the single office 'US Senate GA'. We choose the latter, adding a row to  `CandidateContest.txt` (`US Senate GA (partial term)	1	US Senate GA	`) and two rows to `dictionary.txt`:
```
CandidateContest	US Senate GA (partial term)	US Senate (Loeffler) - Special
CandidateContest	GA Attorney General	GA Attorney General
```

#### Candidates
Candidates not found in the dictionary will be listed in a `.warnings` file. Copy and paste these names into the `Candidate.txt` file as a single BallotName column:
```
BallotName
Zulma Lopez
Zachary Perry
Yasmin Neal
```
and add corresponding rows to `dictionary.txt`:
```
Candidate	Zulma Lopez	Zulma Lopez
Candidate	Zachary Perry	Zachary Perry
Candidate	Yasmin Neal	Yasmin Neal
```
You are free to choose another convention for the candidate names in the database, as long as you specify the mapping in the dictionary. E.g., you could instead have
```
BallotName
Lopez, Zulma
Perry, Zachary
Neal, Yasmin
```
if your dictionary has
```
Candidate	Lopez, Zulma	Zulma Lopez
Candidate	Perry, Zachary	Zachary Perry
Candidate	Neal, Yasmin	Yasmin Neal
```

#### Parties
Parties not recognized will also be listed in the `.warnings` file, e.g.:
```Partys (found with munger my_Georgia_test.munger) not found in dictionary.txt :
Rep
Lib
Ind
Grn
Dem
```
All these should be mapped in `dictionary.txt`
```
Party	Republican Party	Rep
Party	Libertarian Party	Lib
Party	Independent Party	Ind
Party	Green Party	Grn
Party	Democratic Party	Dem
```
to your chosen internal database names which should be in the name column of `Party.txt`:
```
Name
Democratic Party
Republican Party
Libertarian Party
Green Party
Independent Party
```
Note: Some of the routines in the `analyze` submodule assume that every party name ends with ' Party'.


# About the Code 

## Documentation in Progress! Proceed with caution!

## Code components
### About the `CDF_schema_def_info` directory:
The information in this directory determines the structure of the database created by the system to store election results information. Subdirectories and their contents are:
 * `elements` subdirectory contains a subdirectory for each main tables in the database. Most of these correspond to classes in the Common Data Format; other tables (e.g., `_datafile`) start with an underscore. 
 * `enumerations` subdirectory contains a file for each relevant enumerated list from by the Common Data Format. We treat `BallotMeasureSelection` as an enumerated list.
 * `joins` subdirectory contains a subdirectory for each join table in the database.
 
### Some hard-coded items
 Some lists are hard-coded in one place, so could be changed.
  * Recognized jurisdictions `states_and_such` in `database/__init__.py`. Anything not on this list will not appear in the output of `database.display_jurisdictions()`, even if there is corresponding data in the database.
  * Recognized contest types `contest_types_model` in `database/__init__.py`.
  * Recognized encodings `recognized_encodings` in `userinterface/__init__.py`
  * Ballot measure selections `bmselections` in `database/create_cdf_db/__init__.py`
  
### Conventions
 * Some jurisdiction names (e.g., "District of Columbia") contain spaces, while it is inconvenient to use spaces in directory and file names. So we distinguish between a jurisdiction's "true name" and its "system name", which replaces all spaces by hyphens (e.g., "District-of-Columbia").


### Testing
There are `pytest` routines to test the dataloading, analyzer and jurisdiction-prepping functions. The data required to test the latter two is in the [`000_data_for_pytest`](../tests/000_data_for_pytest) folder. These depend on a file `tests/run_time.ini`. Note that the run_time.ini file is *not* part of the repository, as it contains database connection information. You will need to create it yourself. Or you can point to a different parameter file with the custom option `--param_file` for  pytest.

The [dataloader test](../tests/dataloading_tests/test_dataloading.py) depends on information outside the [tests folder] (../tests):
 - election results files in the results directory specified in by the `results_dir` parameter in `dataloading_tests/run_time.ini`. 
    - if the results directory does not exist, the test will create it and pull files from [`https://github.com/ElectionDataAnalysis/TestingData.git`](https://github.com/ElectionDataAnalysis/TestingData.git). You can specify a different url with the custom pytest option `--test_data_url`
 - reference files in the [`reference_results` folder](../src/reference_results). By convention, these are tabs-separated and named for the jurisdiction, e.g., `Virgina.tsv` or `American-Samoa.tsv`. Note the hyphens. If there are no reference results for a given election-jurisdiction pair, the test will fail. The reference files must have columns `Jurisdiction,Election,Contest,ReportingUnit,VoteType,Count`. 

Note that the `analyzer_tests` and `dataloader_tests` directories each have a `conftest.py` file. This may cause a problem if you try to run them simultaneously via `pytest` from the `test` directory. Running them separately works:
```
tests % pytest analyzer_tests
tests % pytest dataloader_tests
```




  
  

# Installation Instructions

## Environment
You will need:
 * `python3.9` 
 * `postgresql`
 * all python packages given in [requirements.txt](../requirements.txt), and any packages those packages may require

Not absolutely required, but recommended, is a package manager, such as `homebrew` for macOS.

### python3.9   
If your environment has a `python` command without a specified version, that `python` may or may not point to `python3.9`. In linux-flavored shells, you can check the version with the command `python --version`. Similarly, your system may have a `python3` command, whose version can be checked with `python3 --version`.

### postgresql
If postgresql is present, the command `postgres --version` will yield the version number; otherwise the command will fail. The `electiondata` package has been tested with postgres version 13, but probably any reasonably recent version will do. If `postgresql` is not present, you should be able to install it with any reasonable package manager. (On macOS with package manager `homebrew`, use the command `brew install postgresql`) The default values you will need to connect to your `postgresql` instance are:
 * host: `localhost`
 * port: `5432`
 * user: `postgres`
 * password: (leave the password blank)

### python packages
To install the required packages, run `python3.9 -m pip install -r requirements.txt` from the [root folder of the repository](../).  Because some of the required packages have requirements of their own, which may or may not be installed already, your system prompt you to install some other packages. If so, install the suggested packages and try `python3.9 -m pip install -r requirements.txt` again.

## Installation
From the [root folder of the repository](../) run `python3.9 setup.py install`. (You may be able to use `python setup.py install` or `python3 setup.py install` instead, if those point to `python3.9`, as described above.)No reports or responses yet.
