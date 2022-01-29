[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1051094.svg)](https://doi.org/10.5281/zenodo.1051094)

ECE2CMOR3 Python code to CMORize and post-process EC-Earth output data.

## Required python packages:

* cmor-3.5.0 (see cmor [dependencies](https://anaconda.org/conda-forge/cmor/files))
* eccodes/gribapi (for filtering IFS output GRIB files)
* dreq (the CMIP6 data request tool drq)
* netCDF4
* cdo (version 1.9.6; only for atmosphere post-processing)
* nose, testfixtures (only for testing)
* pip (for installing python packages)
* f90nml (only for fortran namelist I/O)
* xlrd (for reading *.xlsx excel sheets)
* XlsxWriter (for writing *.xlsx excel sheets)

## Installation:

More extensive installation description can be found [here](https://dev.ec-earth.org/projects/cmip6/wiki/Installation_of_ece2cmor3) at the EC-Earth portal, including the link to an [example of running ece2cmor](https://dev.ec-earth.org/projects/cmip6/wiki/Step-by-step_guide_for_making_CMIP6_experiments#Cmorisation-with-ece2cmor-v120). The basic ece2cmor3 installation description follows below.

#### Installation & running with miniconda (strongly recommended):
The Miniconda python distribution should be installed. With miniconda all the packages can be installed within one go by the package manager `conda`. This applies also to systems where one is not allowed to install complementary python packages to the default python distribution.

##### If Miniconda is not yet installed:

Download [miniconda](https://repo.continuum.io/miniconda/) (e.g. take the latest miniconda version for python 2.7) by using `wget` and install with `bash`:
 ```shell
 wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
 bash Miniconda2-latest-Linux-x86_64.sh -b -u -p /$HOME/miniconda2
 ```
One could consider to add the following aliases in the `.bashrc` file:
 ```shell
 mincondapath=${HOME}/miniconda2/
 alias activateminiconda='source ${mincondapath}/etc/profile.d/conda.sh; export PATH="${mincondapath}/bin:$PATH"'
 alias activateece2cmor3='activateminiconda; conda activate ece2cmor3;'
 ```


##### Download ece3cmor3 by a git checkout

For example we create the directoy ${HOME}/cmorize/ for the ece2cmor tool:

```shell
git clone https://github.com/EC-Earth/ece2cmor3.git
cd ece2cmor3
git submodule update --init --recursive
./download-b2share-dataset.sh ${HOME}/cmorize/ece2cmor3/ece2cmor3/resources/b2share-data
```
Note that Github depricates the `https` clone method, therefore see how to [migrate from https to ssh](https://github.com/EC-Earth/ece2cmor3/wiki/instruction-how-to-change-from-https-to-ssh).

##### Creating a virtual conda environment and installing ece3cmor3 therein:
In the ece2cmor3 git checkout directory, type
```shell
activateminiconda                         # The alias as defined above
conda update -n base -c defaults conda    # for updating conda itself
conda env create -f environment.yml       # for linux & mac os
conda activate ece2cmor3
python setup.py install
```

##### Running ece2cmor3 inside the conda environment:

```shell
 conda activate ece2cmor3
 ece2cmor -h
 checkvars -h
 etc.
 conda deactivate
```

#### Note that the nested CMOR tables require an update once in a while: 

The CMOR tables are maintained via a nested git repository inside the ece2cmor3 git repository. 
Once in a while one of the ece2cmor3 developers will update the nested repository of the CMOR tables. 
This will be visible from the ece2cmor3 repository by a git status call, it will tell that there are "new updates" in these tables. 
In that case one has to repeat the following inside the git checkout directory:
```shell
git submodule update --init --recursive
```

#### Note for developers: 

To avoid many installation calls during development, you can symlink the installed modules to the source directory by executing
```shell
python setup.py develop;
```

#### Updating the nested CMOR table repository by maintainers:
Navigate to your git checkout directory and execute
```shell
cd ${HOME}/cmorize/ece2cmor3/ece2cmor3/resources/tables/
git pull origin master
cd ../; git add cmip6-cmor-tables
git commit cmip6-cmor-tables -m 'Update the nested CMOR tables for their updates'
git push
```

## Design:

The package consists for 2 main modules, ifs2cmor and nemo2cmor. The main api module ece2cmorlib calls initialization and processing functions in these ocean and atmosphere specific codes. The full workload is divided into tasks, which consist of a source (an IFS grib code or NEMO parameter id) and a target (a cmor3 CMIP6 table entry). The tasks are constructed by the Fortran namelist legacy loader (namloader.py) or by the new json-loader (default). The working is similar to the previous ece2cmor tool: the loader reads parameter tables and creates tasks as it receives a dictionary of desired targets from the caller script.

At execution, the nemo2cmor module searches for the sources in the NEMO output files and streams the data to cmor to rewrite it according to CMIP6 conventions. For the IFS component, the module first performs the necessary post-processing steps, creating a list of intermediate netcdf files that contain time-averaged selections of the data. Special treatment such as unit conversions and post-processing formulas are attached to the tasks as attributes, as well as the file path in which the source data resides.
# Contributing guidelines

We welcome any kind of contribution to our software, from simple comment or question to a full fledged 
[pull request](https://help.github.com/articles/about-pull-requests/). 
Please read and follow our [Code of Conduct](CODE_OF_CONDUCT.md).

A contribution can be one of the following cases:

1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update 
documentation).

The sections below outline the steps in each case.

## You have a question

1. use the search functionality [here](https://github.com/EC-Earth/ece2cmor3/issues) to see if someone already filed the 
same issue;
1. if your issue search did not yield any relevant results, make a new issue;
1. apply the "Question" label; apply other labels when relevant.

## You think you may have found a bug

1. use the search functionality [here](https://github.com/EC-Earth/ece2cmor3/issues) to see if someone already filed the 
same issue;
1. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information 
to the rest of the community to understand the cause and context of the problem. 
Depending on the issue, you may want to include:
    - the [SHA hashcode](https://help.github.com/articles/autolinked-references-and-urls/#commit-shas) of the commit 
    that is causing your problem;
    - information about the nature of the EC-Earth output you are working with and possibly the metadata that you are 
    using; 
    - information about the HPC platform or OS you are using for your job;
1. apply relevant labels to the newly created issue.

## You want to make some kind of change to the code base

1. (**important**) announce your plan to the rest of the community _before you start working_. 
This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
1. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest 
master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in 
changes, possibly from the 'upstream' repository (follow the instructions 
[here](https://help.github.com/articles/configuring-a-remote-for-a-fork/) and 
[here](https://help.github.com/articles/syncing-a-fork/));
1. make sure the existing tests still work by running ``nosetests``;
1. add your own tests (if necessary);
1. update or expand the documentation;
1. [push](http://rogerdudler.github.io/git-guide/) your feature branch to the master branch in your fork.
1. create the pull request, e.g. following the instructions [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how 
to generate the documentation: don't let this discourage you from making the pull request, just go 
ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull 
request.

### Adding a new component to ece2cmor3
This are some initial guidelines for people that want to add support for output of other EC-Earth components, such as 
TM5, LPJ-GUESS or PISM. In the ece2cmor3/ece2cmor3 directory, type
```bash
 grep -iHn NEWCOMPONENT *.py
```
This gives an idea where some action is needed, as we labeled all points by ``NEWCOMPONENT``. You will need to add an
item to the dictionary ``models`` in the ``components.py``, e.g.
```python
models = {"ifs": {realms: ["atmos", "atmosChem", "land", "landIce"],
                  table_file: os.path.join(os.path.dirname(__file__), "resources", "ifspar.json"),
                  script_flags: ("atm", 'a')},
          "nemo": {realms: ["ocean", "ocnBgchem", "seaIce"],
                   table_file: os.path.join(os.path.dirname(__file__), "resources", "nemopar.json"),
                   script_flags: ("oce", 'o')},
          "your_model" : {realms: ["your_cmip_realm"],
                   table_file: os.path.join(os.path.dirname(__file__), "resources", "your_table_file.json")}
          }
```
This registers ``your_model`` as a new component during the loading of the tasks, and it will give preference to your 
model if you claim to process data for ``your_cmip_realm`` for such variables. It will look for these variables in the 
json-file ``your_table_file.json``, so you will have to create that, please look at e.g. ``nemopar.json`` for the
schema that we expect; you will need to specify a ``source`` attribute for each variable, e.g. 
```json
...
{
  "source" : "your_var_name",
  "target" : "co2"
},
...
```
This entry will create a task with a ``cmor_source`` containing ``your_var_name`` as variable name and ``co2`` as a 
target. Upon execution ``ece2cmorlib`` contains a list of tasks, some of which belong to your model component. 
From that point, you have to implement the processing of the tasks yourself, although you may want to get inspiration 
from e.g. ``nemo2cmor.py``. Moreover you can use some utilities in the code such as the ``cdo`` python wrapper. 
List of tests
=============

These following tests must be passed by any EC-Earth data set before it is
published on the ESGF data nodes.


- Run `ece2cmor`, check logs for errors  
  (some `ece2cmor` errors may be accepted, see related discussions on
  the EC-Earth Portal and ece2cmor Github page)

- Scripted checks on file-level:

  + empty files or directories:  
    `find . -empty`

  + auxiliary files:  
    `find . -type f ! -name \*.nc`
    (e.g. `*.nc.copy` files)  
    `find . -name '*_r1i1p1f1_g[rn][a-zA-Z0-9]*'` or something like   
    `find . -type f -regextype sed ! -regex '.*/.*_EC-Earth3-Veg_[0-9a-zA-Z]*_r1i1p1f1_g[nr][_0-9-]*\.nc'`
    (temporary files from CMOR library)

  + versions check:  
    `versions.sh -l .` (check present versions)  
    `versions.sh -v vYYYYMMDD .` (set version, dry-run)  
    `versions.sh -v vYYYYMMDD -m .` (set version)

  + Check number of files per year:  
    `files-per-year <FIRST_YEAR> <LAST_YEAR> .`

  + Check number of files per variable:  
    `files-per-variable <VERSION> <VAR_DIRS>`

- Manual check against do-no-publish list  
  (for now, found in [this google document](https://docs.google.com/spreadsheets/d/1b69NCgHSjWNGqalOVWBTdTvXl7R_-EwCRfJqlwYB37Y))

- Fix ocean fraction in `Ofx/sftof` (check out `sftof.py` under `scripts`)
- Fix land fraction in `fx/sftlf` (check under `recipes` for help)

- Suggestion: set all files/directories *read-only* from here on

- Run `nctime`, checking time axes and ranges:  
  `nctcck -p cmip6 -i <ESGINI_DIR> -l <LOG_DIR> <DATA_DIR>`  
  `nctxck -p cmip6 -i <ESGINI_DIR> -l <LOG_DIR> <DATA_DIR>`

- Run `qa-dkrz`, for comprehensive QA checks:  
  checks enabled:
  + DRS, DRS_F, DRS_P
  + CF
  + CNSTY
  + TIME
  + DATA