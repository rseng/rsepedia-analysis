**Dev**

**2.2.1**
- Update app title

**2.2.0**
- Remove base URL

**2.1.1**
- Update Dockerfile with Biocontainers template

**2.1.0**
- Fix style in help/documentation
- Update documentation
- Set default max run time to 48 hours
- Update conda version in Dockerfile
- Update conda dependencies

**2.0.0**
- Remove e-mail support

**1.2.0**
- Fix multiple contents in input log file.
- Results are public
- Use local css and js files for Bootstrap and jQuery
- Update ubuntu and conda version in Docker image

**1.1.0**
- Minor updates
- Documentation of default configuration 
  
**1.0.0**
- Use pathlib in `export_results.py`
- Use conda env instead of pipenv
- Move config file to the config/ directory
- Use logging for messages
- Results are destroyed automatically after some time

**0.1.3**
- Fix bug for discrete date type

**0.1.2**
- Fix name on autoclass failure
- Append multiple log files upon failure
- Increase upload file size up to 100 Mb

**0.1.1**
- Format running time as HH:MM:SS
- Do not send results if job failed
- Mount Docker volumes in the current directory
- Display version in help
- Increase to 20 Gb max size for file upload
# AutoClassWeb

[![License: BSD](https://img.shields.io/badge/License-BSD-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5215902.svg)](https://doi.org/10.5281/zenodo.5215902)
[![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/pierrepo/autoclassweb/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/pierrepo/autoclassweb)

AutoClassWeb is a web interface to [AutoClass C](https://ti.arc.nasa.gov/tech/rse/synthesis-projects-applications/autoclass/autoclass-c/), an unsupervised Bayesian classification system developed by the NASA.

It utilizes [AutoClassWrapper](https://github.com/pierrepo/autoclasswrapper), a Python wrapper for AutoClass C. 

## Installation for use on a local machine

See step-by-step instructions: <https://github.com/pierrepo/autoclassweb-app>


## Installation for use on a web server

See step-by-step instructions: <https://github.com/pierrepo/autoclassweb-server>


## Installation for development

Clone the project:
```bash
git clone https://github.com/pierrepo/autoclassweb.git
cd autoclassweb
```

Create and activate a conda environment:
```bash
conda env create -f environment.yml
conda activate autoclassweb
```

Install [AutoClass C](https://ti.arc.nasa.gov/tech/rse/synthesis-projects-applications/autoclass/autoclass-c/):

```bash
wget https://ti.arc.nasa.gov/m/project/autoclass/autoclass-c-3-3-6.tar.gz
tar zxvf autoclass-c-3-3-6.tar.gz
rm -f autoclass-c-3-3-6.tar.gz
export PATH=$PATH:$(pwd)/autoclass-c
```
If you use a 64-bit operating system, install the standard 32-bit C libraries:
```bash
sudo apt install -y libc6-i386
```

Copy config template and update config file `config/autoclassweb.cfg` accordingly:
```bash
cp config/autoclassweb-template.cfg config/autoclassweb.cfg
```

Run AutoClassWeb alone:
```bash
make run
```

or with gunicorn:
```bash
make run-gunicorn
```

AutoClassWeb is then available at <http://127.0.0.1:5000>

## Docker 

Install Docker with the following [instructions](https://docs.docker.com/install/linux/docker-ce/ubuntu/).

Build image:
```bash
make docker-build
```

Run container:
```bash
make docker-run
```

AutoClassWeb is then available at <http://127.0.0.1:5000>

Clean unused images:
```bash
make docker-clean
```

Official Docker images of AutoClassWeb are available in the [Biocontainers](https://biocontainers.pro/) docker repository: 

<https://hub.docker.com/r/biocontainers/autoclassweb>


# How to release

## Setup

Install required packages:
```bash
$ conda env create -f environment.yml
```

Or if needed, update your conda environment:
```bash
$ conda env update -f environment.yml
```

For Zenodo integration, see [Making Your Code Citable](https://guides.github.com/activities/citable-code/).

## Lock dependencies (if needed)

```
$ conda env export -n autoclassweb --no-builds  | grep -v "^prefix:" > environment-lock.yml
$ git commit -a -m "Update conda dependencies"
```

## Update version number

We use `bump2version` to update and synchronize the version number across different files.

For patch update (x.y.z → x.y.**z+1**):
```
$ bump2version --verbose --config-file devtools/bumpversion.cfg patch
```

For minor update (x.y.z → x.**y+1**.0):
```
$ bump2version --verbose --config-file devtools/bumpversion.cfg minor
```

For major update (x.y.z → **x+1**.0.0):
```
$ bump2version --verbose --config-file devtools/bumpversion.cfg major
```

Remark:

1. For a dry run with `bump2version`, use option `-n`.
2. `bump2version` will fail if the git working directory is not clean, i.e. all changes are not commited.

Once version number is updated, push everything to GitHub:
```
$ git push origin
$ git push origin --tags
```


## Add new release on GitHub

On [GitHub release page](https://github.com/pierrepo/autoclassweb/releases) :

- Click the *Draft a new release* button.
- Select the latest version as *tag version*.
- Add release version as *Release title* (e.g.: v1.3.7).
- Copy and paste the content of the `CHANGELOG.md` in the *Describe this release* field.
- Hit the *Publish release* button :rocket:.


## Zenodo integration

After the creation of the new release in GitHub, check the archive has been creating on [Zenodo](https://doi.org/10.5281/zenodo.5215902).


## Publish docker image 

Docker images of AutoClassWeb are hosted in the [Biocontainers](https://biocontainers.pro/) docker repository:

<https://hub.docker.com/r/biocontainers/autoclassweb>

To push new images:

-  Clone <https://github.com/BioContainers/containers> to <https://github.com/pierrepo/containers>
-  Add the Dockerfile for the new release in the `autoclassweb` subdirectory
-  Make a [pull request](https://github.com/BioContainers/containers/compare/master...pierrepo:master) for merging.

