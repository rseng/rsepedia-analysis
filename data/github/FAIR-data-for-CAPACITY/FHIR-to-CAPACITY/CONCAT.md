## Badges

| fair-software.eu recommendations | |
| :-- | :--  |
| code repository              | [![github repo badge](https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue)](https://github.com/FAIR-data-for-CAPACITY/FHIR-to-CAPACITY) |
| license                      | [![github license badge](https://img.shields.io/github/license/FAIR-data-for-CAPACITY/FHIR-to-CAPACITY)](https://github.com/FAIR-data-for-CAPACITY/FHIR-to-CAPACITY) |
| community registry           | [![RSD](https://img.shields.io/badge/rsd-FHIR_TO_CAPACITY-00a3e3.svg)](https://www.research-software.nl/software/fhir-to-capacity)  |
| citation                     | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4818370.svg)](https://doi.org/10.5281/zenodo.4818370) |
| howfairis                          | [![fair-software badge](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu) |
| **Other best practices**           | &nbsp; |
| Coverage                           | [![Coverage Status](https://coveralls.io/repos/github/FAIR-data-for-CAPACITY/FHIR-to-CAPACITY/badge.svg?branch=master)](https://coveralls.io/github/FAIR-data-for-CAPACITY/FHIR-to-CAPACITY?branch=master)|
| **GitHub Actions**                 | &nbsp; |
| Build                              | [![Build](https://github.com/FAIR-data-for-CAPACITY/FHIR-to-CAPACITY/actions/workflows/build.yml/badge.svg)](https://github.com/FAIR-data-for-CAPACITY/FHIR-to-CAPACITY/actions/workflows/build.yml)|

## Installation

To install FHIR-to-CAPACITY from GitHub repository, do:

```console
git clone https://github.com/FAIR-data-for-CAPACITY/FHIR-to-CAPACITY.git
cd FHIR-to-CAPACITY
python3 -m pip install .
```



## Usage

To upload records from a FHIR server to a CAPACITY registry, execute the following command:
```console
fhir-to-capacity FHIR_BASE_URL CAPACITY_API_URL CAPACITY_TOKEN
```

Using the `--help` option will give you more information:
```console
fhir-to-capacity --help

Usage: scripts/fhir_to_capacity fhir-base-url capacity-url capacity-token

Queries a FHIR server for patients, converts this to records matching the CAPACITY codebook, and uploads the result to a REDCAP server.

Arguments:
  fhir-base-url
  capacity-url
  capacity-token

Other actions:
  -h, --help       Show the help

```

## Development
There is an additional command available for filling a FHIR server with synthetic data. It can 
be run as follows:
```console
fill-server FHIR_BASE_URL
```

The `--help` option will show you more information:
```console
Usage: scripts/fill_server fhir-base [n]

Fills the server with n ramdom patients.

Arguments:
  fhir-base    the url to the FHIR api
  n            number of patients to create (type: INT, default: 10)

Other actions:
  -h, --help   Show the help

```
## Contributing

If you want to contribute to the development of FHIR-to-CAPACITY,
have a look at the [contribution guidelines](CONTRIBUTING.rst).

## Credits

This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter) and the [NLeSC/python-template](https://github.com/NLeSC/python-template).
