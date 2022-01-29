[![PyPI](https://img.shields.io/pypi/v/spottool)](https://pypi.python.org/pypi/spottool)
[![DOI](https://zenodo.org/badge/212620019.svg)](https://zenodo.org/badge/latestdoi/212620019)
[![Build Status](https://travis-ci.org/ali4006/spot.svg?branch=develop)](https://travis-ci.org/big-data-lab-team/spot)
[![Coverage Status](https://coveralls.io/repos/github/big-data-lab-team/spot/badge.svg?branch=develop)](https://coveralls.io/github/big-data-lab-team/spot?branch=develop)

# Spot

Spot identifies the processes in a pipeline that produce different results in different
execution conditions.

<!-- TABLE OF CONTENTS -->
## Table of Contents

* [Installation](#installation)
* [Pre-requisites](#pre-requisites)
* [First example](#first-example)
* [HCP example](#hcp-example)
* [How to Contribute](#how-to-contribute)
* [License](#license)


## Installation

Simply install the package with `pip`

    $ pip install spottool

## Pre-requisites

* Install and start [Docker](http://www.docker.com)
* Build Docker images for the pipelines in different conditions (see [Dockerfile](spot/pfs-example/centos7/Dockerfile) as an example of PreFreeSurfer pipeline in CentOS7)
* Create [Boutiques](https://boutiques.github.io) descriptors for the pipeline, in each condition (see [descriptor.json](spot/pfs-example/descriptor_centos7.json) and [invocation.json](spot/pfs-example/invocation_centos7.json) as an example of PreFreeSurfer pipeline)
* Get provenance information using [ReproZip](http://docs.reprozip.org/en/1.0.x/packing.html) tool in one condition
by running: ```reprozip trace <CMD>```

The `auto_spot` command finds processes that create differences in results obtained in different conditions and reports them in a JSON file.

## First example

In this example, we run a bash script that calls the grep command multiple times, creating different output files when run on different OSes. We use spot to compare the outputs obtained in CentOS 7 and Debian 10.

The example can be run in this Git repository as follows:
```
git clone https://github.com/big-data-lab-team/spot.git
cd spot
pip install .

docker build . -f spot/example/centos7/Dockerfile -t spot_centos_latest
docker build . -f spot/example/debian/Dockerfile -t spot_debian_latest

cd spot/example 

auto_spot -d descriptor_centos7.json -i invocation_centos7.json -d2 descriptor_debian10.json -i2 invocation_debian10.json -s trace_test.sqlite3 -c conditions.txt -e exclude_items.txt -o commands.json .
```
In this command:
* `descriptor_<distro>.json` is the Boutiques descriptor of the application executed in OS `<distro>`.
* `invocation_<distro>.json` is the Boutiques invocation of the application executed in OS `<distro>`, containing the input files.
* `trace.sqlite3` is a ReproZip trace of the application, acquired in CentOS 7.
* `condition.txt` contains the result folder for each condition.
* `exclude_items.txt` contains the list of items to be ignored while parsing the files and directories.

The command produces the following outputs:
*  `commands_captured_c.json` contains the list of processes with temporary files and files written by multiple processes. 
*  `commands.json` contains the list of processes that create differences in two conditions. Attribute `total_commands_multi` contains processes that write files written by multiple processes and `total_commands` contains the other processes.

## HCP example

In this example, we run a short PreFreeSurfer pipeline that includes only the ACPC-Alignment step 
to process only the T1w-image of one subject. The results will show the `FLIRT` tool as the non-reproducible process in the pipeline when running on different versions of CentOS. We use `spot` to compare the outputs obtained in CentOS 7 and CentOS 6.

This example takes ~12 mins running and needs ~500 MB space in total.
Before running the example, make sure `git-lfs` is installed on your operating system
(See the [link](https://github.com/git-lfs/git-lfs/wiki/Installation) ).

The example can be run in this Git repository as follows:
```
git lfs install
git clone https://github.com/big-data-lab-team/spot.git
cd spot
pip install .

docker build . -f spot/pfs-example/centos7/Dockerfile -t short-pfs-spot-centos7
docker build . -f spot/pfs-example/centos6/Dockerfile -t short-pfs-spot-centos6

cd spot/pfs-example 

auto_spot -d descriptor_centos7.json -i invocation_centos7.json -d2 descriptor_centos6.json -i2 invocation_centos6.json -s trace.sqlite3 -c conditions.txt -e exclude_items.txt -o commands.json .
```

Furthermore, we can reorder the executions and then merge the identified processes in two different orders by running:
```
auto_spot -d2 descriptor_centos7.json -i2 invocation_centos7.json -d descriptor_centos6.json -i invocation_centos6.json -s trace.sqlite3 -c conditions2.txt -e exclude_items.txt -o commands2.json .

python ../merge_jsons.py commands.json commands2.json merged.json
```

The command produces the following output:
*  `merged.json` contains the list of processes that create differences in each order of executions.

### Expected output
The `merged.json` file should be similar to the [merged_reference.json](spot/pfs-example/merged_reference.json). 
In this file, the `flirt` process is identified under the attribute `total_commands` as a process that creates different result files, `roi2std.mat` and `acpc_final.nii.gz`.

## How to Contribute

1. Clone repo and create a new branch: `$ git checkout https://github.com/big-data-lab-team/spot -b name_for_new_branch`.
2. Make changes and test
3. Submit Pull Request with comprehensive description of changes


## License

[MIT](LICENSE) © /bin Lab
# Metrics

A set of scripts that take two files as argument and return a number on their standard output.

* `dice.sh`: returns the Dice coefficient between images. Requires `fslmaths`, `fslstats` and Freesurfer's `mri_convert` (environment variable FREESURFER_HOME should be set to a directory containing a valid Freesurfer license in `license.txt`).
* `nrmse.sh`: returns the root mean square error, normalized by the amplitude of the image (max intensity - min intensity).
