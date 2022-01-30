# cffconvert GitHub Action

GitHub action to validate CITATION.cff files, and convert to other citation formats using dockerized version of [cffconvert](https://pypi.org/project/cffconvert/).

[![github repo badge](https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue)](https://github.com/citation-file-format/cffconvert-github-action)
[![github license badge](https://img.shields.io/github/license/citation-file-format/cffconvert-github-action)](https://github.com/citation-file-format/cffconvert-github-action)
[![github marketplace badge](https://img.shields.io/badge/github-marketplace-000.svg?logo=github&labelColor=gray&color=blue)](https://github.com/marketplace/actions/cffconvert)
[![Research Software Directory](https://img.shields.io/badge/rsd-cffconvert--github--action-00a3e3.svg)](https://www.research-software.nl/software/cffconvert-github-action)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3993241.svg)](https://doi.org/10.5281/zenodo.3993241)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu)
[![testing](https://github.com/citation-file-format/cffconvert-github-action/workflows/selftest/badge.svg)](https://github.com/citation-file-format/cffconvert-github-action/actions?query=workflow%3A%22selftest%22)
[![citation metadata](https://github.com/citation-file-format/cffconvert-github-action/workflows/cffconvert/badge.svg)](https://github.com/citation-file-format/cffconvert-github-action/actions?query=workflow%3A%22cffconvert%22)
[![links](https://github.com/citation-file-format/cffconvert-github-action/actions/workflows/link-check.yml/badge.svg)](https://github.com/citation-file-format/cffconvert-github-action/actions/workflows/link-check.yml)
[![fair software badge](https://github.com/citation-file-format/cffconvert-github-action/actions/workflows/fair-software.yml/badge.svg)](https://github.com/citation-file-format/cffconvert-github-action/actions/workflows/fair-software.yml)

## Usage

1. Save the snippet below as ``.github/workflows/cffconvert.yml`` to validate your CITATION.cff on each push.
1. ``git add``, ``commit`` and ``push`` to your GitHub repository
1. Check the _Actions_ tab on your repository's page to check the action's output

```yaml
name: cffconvert

on: push

jobs:
  validate:
    name: "validate"
    runs-on: ubuntu-latest
    steps:
      - name: Check out a copy of the repository
        uses: actions/checkout@v2

      - name: Check whether the citation metadata from CITATION.cff is valid
        uses: citation-file-format/cffconvert-github-action@2.0.0
        with:
          args: "--validate"

```

**You can also look into [advanced examples](README.advanced.md).**
# Advanced examples

You can pass any arguments in the `args:` key that you would to the `cffconvert` CLI tool.
This allows for some interesting usage.

The list of options is right below, and you can see some examples in the next sections.
```
Options:
  -i, --infile PATH               Path to the CITATION.cff input file. If this
                                  option is omitted, './CITATION.cff' is used.
  -o, --outfile PATH              Path to the output file.
  -f, --format [apalike|bibtex|cff|codemeta|endnote|ris|schema.org|zenodo]
                                  Output format.
  -u, --url TEXT                  URL to the CITATION.cff input file.
  -h, --help                      Show help and exit.
  --show-trace                    Show error trace.
  --validate                      Validate the CITATION.cff file and exit.
  --version                       Print version and exit.
```

## Validating from a subdirectory

```yaml
name: cffconvert

on: push

jobs:
  validate:
    name: "validate"
    runs-on: ubuntu-latest
    steps:
      - name: Check out a copy of the repository
        uses: actions/checkout@v2

      - name: Validate a CITATION.cff from a subdirectory
        uses: citation-file-format/cffconvert-github-action@2.0.0
        with:
          args: "--infile ./tests/subdirectory/CITATION.cff --validate"

```

## Converting CITATION.cff to Zenodo metadata format and pushing to the repo

```yaml
name: cffconvert

on: push

jobs:
  convert:
    name: "convert"
    runs-on: ubuntu-latest
    steps:
      - name: Check out a copy of the repository
        uses: actions/checkout@v2

      - name: Convert CITATION.cff to Zenodo metadata format
        uses: citation-file-format/cffconvert-github-action@2.0.0
        with:
          args: "--infile ./CITATION.cff --format zenodo --outfile .zenodo.json"

      - name: Commit and push Zenodo metadata
        run: |
          git config --global user.name 'cffconvert GitHub Action'
          git config --global user.email 'cffconvert@users.noreply.github.com'
          git add .zenodo.json
          git commit -m "Automated update of Zenodo metadata"
          git push

```

## Converting CITATION.cff to Zenodo metadata format and opening a pull request

```yaml
name: cffconvert

on: push

jobs:
  convert:
    name: "convert"
    runs-on: ubuntu-latest
    steps:
      - name: Check out a copy of the repository
        uses: actions/checkout@v2

      - name: Convert CITATION.cff to Zenodo metadata format
        uses: citation-file-format/cffconvert-github-action@2.0.0
        with:
          args: "--infile ./CITATION.cff --format zenodo --outfile .zenodo.json"

      - name: Commit and create a Pull Request with the Zenodo file
        uses: peter-evans/create-pull-request@v3
        with:
          commit-message: ":robot: Update .zenodo.json"
          committer: "cffconvert GitHub Action <cffconvert@users.noreply.github.com>"
          title: "[auto] Update .zenodo.json"

```
There's neither a CITATION.cff nor a .zenodo.json in this directory.**Description**

<!-- Description of PR -->

**Related issues**:
- ...

**Instructions to review the pull request**

<!--
Clone and verify
```
cd $(mktemp -d --tmpdir cffaction-XXXXXX)
git clone https://github.com/citation-file-format/cffconvert-github-action .
git checkout <this-branch>
```
-->

<!--
Review online.
-->
