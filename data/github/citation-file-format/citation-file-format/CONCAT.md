# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.2.0] - 2021-05-31

### Added

- introduced `preferred-citation`
- root document now has `type:` `software` or `dataset`
- `identifiers` have gained an optional key `description`
- added regex validation for identifiers of type `swh`
- `description`s and `examples` added to schema
- schema has more checks for empty strings or empty arrays (e.g. `authors`, `abstract`, `keywords`)

### Changed

- switched from YAML schema to JSON schema
- `issue`, `number`, `version`, `post-code`, `section` are more lenient now with type union `str|number`
- `loc-start`, `loc-end`, `start`, `end`, `number-volumes`, `volume`, `pages`, `year`, `year-original` are more lenient now with type union `int|str`
- `month` is more lenient now with type union `int[1-12]|enum["1"-"12"]`
- `version`, `date-released` no longer required elements of a minimal CFF
- list of valid license SPDX codes updated to version 2021-05-14
- `language` more lenient (for performance)
- regex for `url` simplified for easier maintenance
- url regex now also allows for `sftp`
- regex for `isbn` simplified and changed for easier maintenance
- dates are now (only) strings of `format` `date` & `pattern` YYYY-MM-DD

### Fixed

- `references.term` was added

## [1.1.0] - 2021-05-31

### Added

- `identifiers` field that accepts `doi`, `swh` (short for "Software Heritage"), `url` and `other` identifiers
- `alias` for *person* objects to accept aliases such as user/screen names, etc.

### Changed

- Make `given-names` and `family-names` non-required in *person* objects

### Fixed

- Note in schema comment that DOI can include slashes

## [1.0.3-4] - 2019-11-05

### Added

- Link to concept DOI rather than single version DOI (https://github.com/citation-file-format/citation-file-format.github.io/commit/4a55d34fc7fb760a0930c2f934f97d4cbe1eb52f)
- Note arbitrary order of keys (https://github.com/citation-file-format/citation-file-format.github.io/commit/aa63008257011c86f120a6498abaa9c08e04ce0d)
- Add Neil Chue Hong to authors
- Add Jurriaan Spaaks to authors

### Fixed

- `orcid` accepts URLs (https://github.com/citation-file-format/citation-file-format/issues/43)
- Fix some typos in the specifications
- Improve rendering of documentation

## [1.0.2] - 2017-12-20

### Fixed

- Fix examples in versions 1.0.0 and 1.0.1 (https://github.com/citation-file-format/citation-file-format/issues/41)
- Indicate required keys for software metadata in the specifications document (https://github.com/citation-file-format/citation-file-format/issues/42)
- Add known bugs to specifications versions 1.0.0 and 1.0.1 (https://github.com/citation-file-format/citation-file-format.github.io/commit/de801adaf86496d4d1948f90e61a699de208508f)

## [1.0.1] - 2017-12-18

### Fixed

- `patent-states` takes a sequence of strings, not a single string (https://github.com/citation-file-format/citation-file-format/issues/39)
- `senders` also accepts entity objects (https://github.com/citation-file-format/citation-file-format/issues/40)

## [1.0.0] - 2017-12-12

### Changed

- Initial full release of the Citation File Format specifications

## [0.9-RC1] - 2017-10-06

### Added

- Initial pre-release of the specifications for the Citation File Format

[unreleased]: https://github.com/citation-file-format/citation-file-format/compare/1.2.0...HEAD
[1.2.0]: https://doi.org/10.5281/zenodo.5171937
[1.1.0]: https://doi.org/10.5281/zenodo.4813122
[1.0.3-4]: https://doi.org/10.5281/zenodo.3515946
[1.0.2]: https://doi.org/10.5281/zenodo.1120256
[1.0.1]: https://doi.org/10.5281/zenodo.1117789
[1.0.0]: https://doi.org/10.5281/zenodo.1108269
[0.9-RC1]: https://doi.org/10.5281/zenodo.1003150
# Guide to Citation File Format schema version 1.2.0

## General structure of a `CITATION.cff` file

Valid Citation File Format files

1. must be named `CITATION.cff` (note the capitalization);
1. are valid according to the Citation File Format schema version 1.2.0 outlined in [schema.json](schema.json);
1. are valid YAML 1.2 ([specification](http://yaml.org/spec/1.2/spec.html), [validator](http://www.yamllint.com/)).

<a name="yaml-strings"></a>**String quoting:** Note that in YAML you generally don't need to quote strings.
But you should use `"` quotes when a string value
contains whitespace,
contains special characters (e.g., any of `:{}[],&*#?|-<>=!%`, or any of `` ` `` and `@` at the beginning),
consists only of numbers (e.g., is the string `"42"`, not the number `42`),
or is `"true"`, `"false"`, `"yes"` or `"no"`.

In short: When a string value doesn't behave as expected, try putting it in `"` quotes.

### Minimal example

A minimal example of a valid `CITATION.cff` file, that contains only the required keys, could look like this:

```yaml
authors:
  - family-names: Druskat
    given-names: Stephan
cff-version: 1.2.0
message: "If you use this software, please cite it using these metadata."
title: "My Research Software"
```

### Typical example

For most software however, it is relatively easy to expand the minimal case with important information like the version, the date when it was last published, some keywords, etc.:

```yaml
abstract: "This is my awesome research software. It does many things."
authors:
  - family-names: Druskat
    given-names: Stephan
    orcid: "https://orcid.org/0000-0003-4925-7248"
cff-version: 1.2.0
date-released: "2021-07-18"
identifiers:
  - description: "This is the collection of archived snapshots of all versions of My Research Software"
    type: doi
    value: 10.5281/zenodo.123456
  - description: "This is the archived snapshot of version 0.11.2 of My Research Software"
    type: doi
    value: 10.5281/zenodo.123457
keywords:
  - "awesome software"
  - research
license: Apache-2.0
message: "If you use this software, please cite it using these metadata."
repository-code: "https://github.com/citation-file-format/my-research-software"
title: "My Research Software"
version: 0.11.2
```

### Referencing other work

When your software or data builds on what others have already done, it is good practice to add a `references` section to your
`CITATION.cff` file. This way, you give credit to the authors of works that your own work builds on.

The `references` section works like a reference list in a paper,
but can also contain references to other software
(such as the dependencies of your software), other datasets,
and other (research) outputs.

```yaml
authors:
  - family-names: Druskat
    given-names: Stephan
cff-version: 1.2.0
message: "If you use this software, please cite it using these metadata."
references:
  - authors:
      - family-names: Spaaks
        given-names: "Jurriaan H."
    title: "The foundation of Research Software"
    type: software
  - authors:
      - family-names: Haines
        given-names: Robert
    title: "Ruby CFF Library"
    type: software
    version: 1.0
title: "My Research Software"
```

### Credit redirection

Sometimes you want to redirect any credit your work may receive towards a second work (typically one of your own). A
common example is, that when you write software and then write a paper about it, you may want to be credited for the paper
instead of for the software itself. For this case, your `CITATION.cff` should contain metadata about the software at
the root of the `CITATION.cff` file, but additionally, you can add a `preferred-citation` key with the metadata of
the paper (or other work) you want people to cite. Usually, the `message` also reflects the authors' wishes on how they want to be credited.

```yaml
authors:
  - family-names: Druskat
    given-names: Stephan
cff-version: 1.2.0
message: "If you use this software, please cite both the article from preferred-citation and the software itself."
preferred-citation:
  authors:
    - family-names: Druskat
      given-names: Stephan
  title: "Software paper about My Research Software"
  type: article
title: "My Research Software"
```

The next sections explain each key in more detail.

## Valid keys

This section describes the valid keys in a `CITATION.cff` file.

### Index

- [`abstract`](#abstract)
- [`authors`](#authors) (array of objects)
- [`cff-version`](#cff-version)
- [`commit`](#commit)
- [`contact`](#contact) (object)
- [`date-released`](#date-released)
- [`doi`](#doi)
- [`identifiers`](#identifiers) (array of objects)
- [`keywords`](#keywords)
- [`license`](#license)
- [`license-url`](#license-url)
- [`message`](#message)
- [`preferred-citation`](#preferred-citation) (object)
- [`references`](#references) (array of objects)
- [`repository`](#repository)
- [`repository-artifact`](#repository-artifact)
- [`repository-code`](#repository-code)
- [`title`](#title)
- [`type`](#type)
- [`url`](#url)
- [`version`](#version)

### `abstract`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: A description of the software or dataset.
- **usage**:<br><br>
    ```yaml
    abstract: This software implements methods to do things.
    ```

### `authors`

- **type**: Array of [`definitions.person`](#definitionsperson) and/or [`definitions.entity`](#definitionsentity) objects.
- **required**: `true`
- **description**: The authors of a software or dataset.
(See also [How to deal with unknown individual authors?](#how-to-deal-with-unknown-individual-authors))
- **usage**:<br><br>
    ```yaml
    authors:
      - given-names: Stephan
        family-names: Druskat
    ```
    ```yaml
    authors:
      - name: "The Research Software project"
    ```
    ```yaml
    authors:
      - given-names: Stephan
        family-names: Druskat
      - name: "The Research Software project"
    ```

### `cff-version`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `true`
- **description**: The Citation File Format schema version that the `CITATION.cff` file adheres to for providing the citation metadata.
- **usage**:<br><br>
    ```yaml
    cff-version: 1.2.0
    ```
    ```yaml
    cff-version: "1.2.0"
    ```

### `commit`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The commit hash or revision number of the software version.
- **usage**:<br><br>
    ```yaml
    commit: 1ff847d81f29c45a3a1a5ce73d38e45c2f319bba
    ```
    ```yaml
    commit: "Revision: 8612"
    ```

### `contact`

- **type**: Array of [`definitions.person`](#definitionsperson) and/or [`definitions.entity`](#definitionsentity) objects.
- **required**: `false`
- **description**: The contact person, group, company, etc. for the software or dataset.
- **usage**:<br><br>
    ```yaml
    contact:
      - affiliation: "German Aerospace Center (DLR)"
        email: mail@research-project.org
        family-names: Druskat
        given-names: Stephan
    ```
    ```yaml
    contact:
      - email: mail@research-project.org
        name: "The Research Software project"
    ```
    ```yaml
    contact:
      - email: mail@research-project.org
        family-names: Druskat
        given-names: Stephan
      - email: mail@research-project.org
        name: "The Research Software project"
    ```

### `date-released`

- **type**: [`definitions.date`](#definitionsdate)
- **required**: `false`
- **description**: The date the software or data set has been released. Format is 4-digit year, 2-digit month, 2-digit day of month, separated by dashes.
- **usage**:<br><br>
    ```yaml
    date-released: "2020-01-31"
    ```

### `doi`

- **type**: [`definitions.doi`](#definitionsdoi)
- **required**: `false`
- **description**: The [DOI](https://en.wikipedia.org/wiki/Digital_object_identifier) of the software or dataset. This notation is most useful when there is just one DOI you want to include. In
that case, `doi` can be used as shorthand for something like:<br><br>
    ```yaml
    identifiers:
      - description: "The concept DOI of the work."
        type: doi
        value: 10.5281/zenodo.1003149
   ```
    or
    ```yaml
    identifiers:
      - description: "The versioned DOI of the work."
        type: doi
        value: 10.5281/zenodo.4813122
   ```
- **usage**:<br><br>
    ```yaml
    doi: 10.5281/zenodo.1003149
    ```
    ```yaml
    doi: "10.5281/zenodo.4813122"
    ```

### `identifiers`

- **type**: Array of [`definitions.identifier`](#definitionsidentifier) objects.
- **required**: `false`
- **description**: The identifiers of the software or dataset.
- **usage**:<br><br>
    ```yaml
    identifiers:
      - description: "The concept DOI of the work."
        type: doi
        value: 10.5281/zenodo.1003149
    ```

    ```yaml
    identifiers:
      - description: "The versioned DOI for version 1.1.0 of the work."
        type: doi
        value: 10.5281/zenodo.4813122
    ```

    ```yaml
    identifiers:
      - description: "The concept DOI of the work."
        type: doi
        value: 10.5281/zenodo.1003149
      - description: "The versioned DOI for version 1.1.0 of the work."
        type: doi
        value: 10.5281/zenodo.4813122
    ```

    ```yaml
    identifiers:
      - description: "The concept DOI of the work."
        type: doi
        value: 10.5281/zenodo.1003149
      - description: "The versioned DOI for version 1.1.0 of the work."
        type: doi
        value: 10.5281/zenodo.4813122
      - description: "The Software Heritage identifier for version 1.1.0 of the work."
        type: swh
        value: "swh:1:dir:bc286860f423ea7ced246ba7458eef4b4541cf2d"
      - description: "The GitHub release URL of tag 1.1.0."
        type: url
        value: "https://github.com/citation-file-format/citation-file-format/releases/tag/1.1.0"
      - description: "The GitHub release URL of the commit tagged with 1.1.0."
        type: url
        value: "https://github.com/citation-file-format/citation-file-format/tree/16192bf05e99bcb35d5c3e085047807b5720fafc"
    ```

### `keywords`

- **type**: Array of [nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: Keywords that describe the work.
- **usage**:<br><br>
    ```yaml
    keywords:
     - thefirstkeyword
     - thesecondkeyword
     - "a third keyword"
    ```

### `license`

- **type**: (Array of) [`definitions.license-enum`](#definitionslicense-enum).
- **required**: `false`
- **description**: The [SPDX license identifier(s)](https://spdx.dev/ids/) for the license(s) under which the work is made available. When there are multiple
licenses, it is assumed their relationship is OR, not AND.
- **usage**:<br><br>
    ```yaml
    license: Apache-2.0
    ```
    ```yaml
    license:
     - Apache-2.0
     - MIT
    ```
    ```yaml
    license:
     - GPL-3.0
     - GPL-3.0-or-later
    ```

### `license-url`

- **type**: [`definitions.url`](#definitionsurl)
- **required**: `false`
- **description**: The URL of the license text under which the software or dataset is licensed (only for non-standard licenses not included in the [SPDX License List](#definitionslicense-enum)).
- **usage**:<br><br>
    ```yaml
    license-url: "https://obscure-licenses.com?id=1234"
    ```

### `message`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `true`
- **default**: `If you use this software, please cite it using the metadata from this file.`
- **description**: A message to the human reader of the `CITATION.cff` file to let them know what to do with the citation metadata.
- **usage**:<br><br>
    ```yaml
    message: "If you use this software, please cite it using the metadata from this file."
    ```
    ```yaml
    message: "Please cite this software using these metadata."
    ```
    ```yaml
    message: "Please cite this software using the metadata from 'preferred-citation'."
    ```
    ```yaml
    message: "If you use this dataset, please cite it using the metadata from this file."
    ```
    ```yaml
    message: "Please cite this dataset using these metadata."
    ```
    ```yaml
    message: "Please cite this dataset using the metadata from 'preferred-citation'."
    ```

### `preferred-citation`

- **type**: A [`definitions.reference`](#definitionsreference) object.
- **required**: `false`
- **description**: A reference to another work that should be cited instead of the software or dataset itself.
Note that the principles of [software citation](https://doi.org/10.7717/peerj-cs.86) and [data citation](https://doi.org/10.25490/a97f-egyk) require that
software should be cited on the same basis as any other research product such as a paper or a book.
Adding a different preferred citation may result in a violation of the respective
primary principle, "Importance", when others cite this work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      authors:
        - family-names: Famnames
          given-names: "Given Nam E."
      title: "Title of the work."
      type: generic
      year: 2021
    ```

### `references`

- **type**: Array of [`definitions.reference`](#definitionsreference) objects.
- **required**: `false`
- **description**: Reference(s) to other creative works. Similar to a list of references in a paper, references of the software or dataset may include other software (dependencies), or other research products that the software or dataset builds on, but not work describing the software or dataset.
- **usage**:<br><br>
    ```yaml
    references:
      - authors:
          - name: "The Dependency Project"
        date-released: "2021-07-26"
        doi: 10.5281/zenodo.x1234567
        repository-code: "https://github.com/dependency-project/dependency"
        title: Dependency
        type: software
        version: 0.13.4
      - authors:
          - family-names: Bielefeld
            given-names: Arthur
            name-particle: von
        doi: 10.9999/hardscifi-lang.42132
        issue: 13
        journal: "Journal of Hard Science Fiction"
        scope: "Cite this paper if you want to reference the general concepts of the software."
        title: "Towards a 100% accuracy syntax parser for all languages"
        type: article
        volume: 42
        year: 2099
    ```

### `repository`

- **type**: [`definitions.url`](#definitionsurl)
- **required**: `false`
- **description**: The URL of the software or dataset in a repository/archive (when the repository is neither a source code repository nor a build artifact repository).
- **usage**:<br><br>
    ```yaml
    repository: "https://ascl.net/2105.013"
    ```

### `repository-artifact`

- **type**: [`definitions.url`](#definitionsurl)
- **required**: `false`
- **description**: The URL of the work in a build artifact/binary repository (when the work is software).
- **usage**:<br><br>
    ```yaml
    repository-artifact: "https://search.maven.org/artifact/org.corpus-tools/cff-maven-plugin/0.4.0/maven-plugin"
    ```

### `repository-code`

- **type**: [`definitions.url`](#definitionsurl)
- **required**: `false`
- **description**: The URL of the work in a source code repository.
- **usage**:<br><br>
    ```yaml
    repository-code: "https://github.com/citation-file-format/cff-converter-python"
    ```

### `title`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `true`
- **description**: The name of the software or dataset.
- **usage**:<br><br>
    ```yaml
    title: "cffconvert"
    ```

### `type`

- **type**: enum (`software` or `dataset`)
- **default**: `software`
- **required**: `false`
- **description**: The type of the work that is being described by this `CITATION.cff` file.
- **usage**:<br><br>
    ```yaml
    type: software
    ```
    ```yaml
    type: dataset
    ```

### `url`

- **type**: [`definitions.url`](#definitionsurl)
- **required**: `false`
- **description**: The URL of a landing page/website for the software or dataset.
- **usage**:<br><br>
    ```yaml
    url: "https://citation-file-format.github.io/"
    ```

### `version`

- **type**: [`definitions.version`](#definitionsversion)
- **required**: `false`
- **description**: The version of the software or dataset.
- **usage**:<br><br>
    ```yaml
    version: "1.2.0"
    ```
    ```yaml
    version: 1.2
    ```
    ```yaml
    version: "21.10 (Impish Indri)"
    ```

## Definitions

Some values in CFF files are valid in different keys.
For example, `repository-code`, `url` and `license-url` all take URLs as values.

The schema therefore has [*definitions*](https://json-schema.org/understanding-json-schema/structuring.html#definitions)
of smaller subschemas, that can be reused in the schema from the respective key,
instead of having to duplicate them.
For example, there is one definition for a valid URL value ([`definitions.url`](#definitionsurl)),
that is being referenced from `repository-code`, `url` and `license-url`.

**Note:** `definitions` and its subkeys like `definitions.alias` or `definitions.entity.alias` should not be used as keys in `CITATION.cff` files:
```yaml
# incorrect
authors:
  - definitions.alias: sdruskat
```
```yaml
# incorrect
authors:
  - definitions:
      alias: sdruskat
```
```yaml
# correct
authors:
  - alias: sdruskat
```

### Index

- [`definitions.address`](#definitionsaddress)
- [`definitions.alias`](#definitionsalias)
- [`definitions.city`](#definitionscity)
- [`definitions.commit`](#definitionscommit)
- [`definitions.country`](#definitionscountry)
- [`definitions.date`](#definitionsdate)
- [`definitions.doi`](#definitionsdoi)
- [`definitions.email`](#definitionsemail)
- [`definitions.entity`](#definitionsentity) (object)
- [`definitions.fax`](#definitionsfax)
- [`definitions.identifier`](#definitionsidentifier) (object)
- [`definitions.identifier-description`](#definitionsidentifier-description)
- [`definitions.license`](#definitionslicense)
- [`definitions.license-enum`](#definitionslicense-enum)
- [`definitions.orcid`](#definitionsorcid)
- [`definitions.person`](#definitionsperson) (object)
- [`definitions.post-code`](#definitionspost-code)
- [`definitions.reference`](#definitionsreference) (object)
- [`definitions.region`](#definitionsregion)
- [`definitions.swh-identifier`](#definitionsswh-identifier)
- [`definitions.tel`](#definitionstel)
- [`definitions.url`](#definitionsurl)
- [`definitions.version`](#definitionsversion)

### `definitions.address`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: N/A
- **description**: An address.
- **usage**:<br><br>
    ```yaml
    authors:
      - name: "The Research Software Project"
        address: "742 Evergreen Terrace"
    ```

### `definitions.alias`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: N/A
- **description**: An alias.
- **usage**:<br><br>
    ```yaml
    authors:
      - name: "The Research Software Project"
        alias: "RSP"
    ```

### `definitions.city`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: N/A
- **description**: A city.
- **usage**:<br><br>
    ```yaml
    authors:
      - name: "The Research Software Project"
        city: "Berlin"
    ```

### `definitions.commit`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: N/A
- **description**: The (e.g., Git) commit hash or (e.g., Subversion) revision number of the work.
- **usage**:<br><br>
    ```yaml
    commit: "1ff847d81f29c45a3a1a5ce73d38e45c2f319bba"
    ```
    ```yaml
    commit: "Revision: 8612"
    ```

### `definitions.country`

- **type**: `enum` with values
    - `AD`
    - `AE`
    - `AF`
    - `AG`
    - `AI`
    - `AL`
    - `AM`
    - `AO`
    - `AQ`
    - `AR`
    - `AS`
    - `AT`
    - `AU`
    - `AW`
    - `AX`
    - `AZ`
    - `BA`
    - `BB`
    - `BD`
    - `BE`
    - `BF`
    - `BG`
    - `BH`
    - `BI`
    - `BJ`
    - `BL`
    - `BM`
    - `BN`
    - `BO`
    - `BQ`
    - `BR`
    - `BS`
    - `BT`
    - `BV`
    - `BW`
    - `BY`
    - `BZ`
    - `CA`
    - `CC`
    - `CD`
    - `CF`
    - `CG`
    - `CH`
    - `CI`
    - `CK`
    - `CL`
    - `CM`
    - `CN`
    - `CO`
    - `CR`
    - `CU`
    - `CV`
    - `CW`
    - `CX`
    - `CY`
    - `CZ`
    - `DE`
    - `DJ`
    - `DK`
    - `DM`
    - `DO`
    - `DZ`
    - `EC`
    - `EE`
    - `EG`
    - `EH`
    - `ER`
    - `ES`
    - `ET`
    - `FI`
    - `FJ`
    - `FK`
    - `FM`
    - `FO`
    - `FR`
    - `GA`
    - `GB`
    - `GD`
    - `GE`
    - `GF`
    - `GG`
    - `GH`
    - `GI`
    - `GL`
    - `GM`
    - `GN`
    - `GP`
    - `GQ`
    - `GR`
    - `GS`
    - `GT`
    - `GU`
    - `GW`
    - `GY`
    - `HK`
    - `HM`
    - `HN`
    - `HR`
    - `HT`
    - `HU`
    - `ID`
    - `IE`
    - `IL`
    - `IM`
    - `IN`
    - `IO`
    - `IQ`
    - `IR`
    - `IS`
    - `IT`
    - `JE`
    - `JM`
    - `JO`
    - `JP`
    - `KE`
    - `KG`
    - `KH`
    - `KI`
    - `KM`
    - `KN`
    - `KP`
    - `KR`
    - `KW`
    - `KY`
    - `KZ`
    - `LA`
    - `LB`
    - `LC`
    - `LI`
    - `LK`
    - `LR`
    - `LS`
    - `LT`
    - `LU`
    - `LV`
    - `LY`
    - `MA`
    - `MC`
    - `MD`
    - `ME`
    - `MF`
    - `MG`
    - `MH`
    - `MK`
    - `ML`
    - `MM`
    - `MN`
    - `MO`
    - `MP`
    - `MQ`
    - `MR`
    - `MS`
    - `MT`
    - `MU`
    - `MV`
    - `MW`
    - `MX`
    - `MY`
    - `MZ`
    - `NA`
    - `NC`
    - `NE`
    - `NF`
    - `NG`
    - `NI`
    - `NL`
    - `NO`
    - `NP`
    - `NR`
    - `NU`
    - `NZ`
    - `OM`
    - `PA`
    - `PE`
    - `PF`
    - `PG`
    - `PH`
    - `PK`
    - `PL`
    - `PM`
    - `PN`
    - `PR`
    - `PS`
    - `PT`
    - `PW`
    - `PY`
    - `QA`
    - `RE`
    - `RO`
    - `RS`
    - `RU`
    - `RW`
    - `SA`
    - `SB`
    - `SC`
    - `SD`
    - `SE`
    - `SG`
    - `SH`
    - `SI`
    - `SJ`
    - `SK`
    - `SL`
    - `SM`
    - `SN`
    - `SO`
    - `SR`
    - `SS`
    - `ST`
    - `SV`
    - `SX`
    - `SY`
    - `SZ`
    - `TC`
    - `TD`
    - `TF`
    - `TG`
    - `TH`
    - `TJ`
    - `TK`
    - `TL`
    - `TM`
    - `TN`
    - `TO`
    - `TR`
    - `TT`
    - `TV`
    - `TW`
    - `TZ`
    - `UA`
    - `UG`
    - `UM`
    - `US`
    - `UY`
    - `UZ`
    - `VA`
    - `VC`
    - `VE`
    - `VG`
    - `VI`
    - `VN`
    - `VU`
    - `WF`
    - `WS`
    - `YE`
    - `YT`
    - `ZA`
    - `ZM`
    - `ZW`
- **required**: N/A
- **description**: The [ISO 3166-1 alpha-2 country code](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2) for a country.
- **usage**:<br><br>
    ```yaml
    authors:
      - country: NL
        name: "The Authors Team"
    ```
    ```yaml
    references:
      - conference:
          country: DE
        type: conference
    ```

### `definitions.date`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: N/A
- **description**: A date. Format is 4-digit year, followed by 2-digit month, followed by 2-digit day of month, and separated by dashes. Note to tool implementers: it is necessary to cast YAML `date` objects to `string` objects when validating against the schema.
- **usage**:<br><br>
    ```yaml
    date-released: "2020-01-31"
    ```
    ```yaml
    references:
      - date-accessed: "2020-01-31"
        type: generic
    ```
    ```yaml
    references:
      - date-downloaded: "2020-01-31"
        type: generic
    ```
    ```yaml
    references:
      - date-published: "2020-01-31"
        type: generic
    ```
    ```yaml
    references:
      - date-released: "2020-01-31"
        type: generic
    ```
    ```yaml
    references:
      - conference:
            date-end: "2020-02-02"
            date-start: "2020-01-31"
        type: conference
    ```

### `definitions.doi`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: N/A
- **description**: The [DOI](https://en.wikipedia.org/wiki/Digital_object_identifier) of the work (i.e., `10.5281/zenodo.1003150`, not the resolver URL `http://doi.org/10.5281/zenodo.1003150`).
- **usage**:<br><br>
    ```yaml
    doi: 10.5281/zenodo.1003150
    ```

### `definitions.email`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: N/A
- **description**: An email address.
- **usage**:<br><br>
    ```yaml
    authors:
      - email: "mail@research-project.org"
        name: "The Research Software project"
    ```

### `definitions.entity`

- **type**: `object` with the following keys:
    - [`address`](#definitionsentityaddress)
    - [`alias`](#definitionsentityalias)
    - [`city`](#definitionsentitycity)
    - [`country`](#definitionsentitycountry)
    - [`date-end`](#definitionsentitydate-end)
    - [`date-start`](#definitionsentitydate-start)
    - [`email`](#definitionsentityemail)
    - [`fax`](#definitionsentityfax)
    - [`location`](#definitionsentitylocation)
    - [`name`](#definitionsentityname)
    - [`orcid`](#definitionsentityorcid)
    - [`post-code`](#definitionsentitypost-code)
    - [`region`](#definitionsentityregion)
    - [`tel`](#definitionsentitytel)
    - [`website`](#definitionsentitywebsite)
- **required**: N/A
- **description**: An entity. Entities are used in keys that can also take [`definitions.person`](#definitionsperson) objects. An entity can represent different types of entities, such as a team, an institution, a company, a conference, etc.
- **usage**:<br><br>
    ```yaml
    authors:
      - address: "742 Evergreen Terrace"
        alias: RSP
        city: Berlin
        country: DE
        date-end: "2018-07-27"
        date-start: "2021-07-27"
        email: team@research-project.org
        fax: +12-3456-7890
        location: "Lovelace Building, room 0.42"
        name: "The Research Software Project"
        orcid: "https://orcid.org/0000-0003-4925-7248"
        post-code: 90210
        region: Renfrewshire
        tel: +12-345-6789098
        website: "https://research-software-project.org"
    ```
    ```yaml
    contact:
      - name: "The Research Software Project team"
    ```
    ```yaml
    references:
      - type: generic
        title: "A reference showing different keys that take entity objects"
        authors:
          - name: "The Research Software Project team"
        conference:
          name: "RC21 - Research Conference 2021"
        contact:
          - name: "Customer Support"
        database-provider:
          name: "Database Provider"
        editors:
          - name: "The Publication Editing Team"
        editors-series:
          - name: "The Publication Series Editing Team"
        institution:
          name: "Department of Research, Random University"
        location:
          name: "Museum of Postmodern Art"
        publisher:
          name: "Open Access Publishing House"
        recipients:
          - name: "The recipient institution of a personal communication"
        senders:
          - name: "The team sending a personal communication"
        translators:
          - name: "Research Translators, Ltd."
    ```

### `definitions.entity.address`

- **type**: [`definitions.address`](#definitionsaddress).
- **required**: `false`
- **description**: The entity's address.
- **usage**:<br><br>
    ```yaml
    authors:
      - address: "742 Evergreen Terrace"
        name: "The Research Software Project"
    ```

### `definitions.entity.alias`

- **type**: [`definitions.alias`](#definitionsalias).
- **required**: `false`
- **description**: The entity's alias.
- **usage**:<br><br>
    ```yaml
    authors:
      - alias: NASA
        name: "National Aeronautics and Space Administration"
    ```

### `definitions.entity.city`

- **type**: [`definitions.city`](#definitionscity).
- **required**: `false`
- **description**: The entity's city.
- **usage**:<br><br>
    ```yaml
    authors:
      - city: Berlin
        name: "The Research Software Project"
    ```

### `definitions.entity.country`

- **type**: [`definitions.country`](#definitionscountry).
- **required**: `false`
- **description**: The entity's country.
- **usage**:<br><br>
    ```yaml
    authors:
      - name: "The Research Software Project"
        country: DE
    ```

### `definitions.entity.date-end`

- **type**: [`definitions.date`](#definitionsdate).
- **required**: `false`
- **description**: The entity's ending date, e.g. when the entity is a conference.
- **usage**:<br><br>
    ```yaml
    references:
      - authors:
          - name: "The Research Software Project"
        conference:
          date-end: "2021-07-27"
          name: "Research Conference 2021"
        title: "Conference Paper"
        type: conference-paper
    ```

### `definitions.entity.date-start`

- **type**: [`definitions.date`](#definitionsdate).
- **required**: `false`
- **description**: The entity's starting date, e.g. when the entity is a conference.
- **usage**:<br><br>
    ```yaml
    references:
      - type: conference-paper
        title: "Conference Paper"
        authors:
          - name: "The Research Software Project"
        conference:
          name: "Research Conference 2021"
          date-start: "2021-07-27"
    ```

### `definitions.entity.email`

- **type**: [`definitions.email`](#definitionsemail).
- **required**: `false`
- **description**: The entity's email address.
- **usage**:<br><br>
    ```yaml
  authors:
    - email: team@research-project.org
      name: "The Research Software Project"
    ```

### `definitions.entity.fax`

- **type**: [`definitions.fax`](#definitionsfax).
- **required**: `false`
- **description**: The entity's fax number.
- **usage**:<br><br>
    ```yaml
    authors:
      - fax: +12-3456-7890
        name: "The Research Software Project"
    ```

### `definitions.entity.location`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The entity's location.
- **usage**:<br><br>
    ```yaml
    authors:
      - location: "Lovelace Building, room 0.42"
        name: "The Research Software Project"
    ```

### `definitions.entity.name`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `true`
- **description**: The entity's name.
- **usage**:<br><br>
    ```yaml
    authors:
      - name: "The Research Software Project"
    ```

### `definitions.entity.orcid`

- **type**: [`definitions.orcid`](#definitionsorcid).
- **required**: `false`
- **description**: The entity's [ORCID](https://orcid.org) identifier.
- **usage**:<br><br>
    ```yaml
    authors:
      - name: "The Research Software Project"
        orcid: "https://orcid.org/0000-0003-4925-7248"
    ```

### `definitions.entity.post-code`

- **type**: [`definitions.post-code`](#definitionspost-code).
- **required**: `false`
- **description**: The entity's post code.
- **usage**:<br><br>
    ```yaml
    authors:
      - name: "The Research Software Project"
        post-code: 90210
    ```

    ```yaml
    authors:
      - name: "The Research Software Project"
        post-code: "90210"
    ```

### `definitions.entity.region`

- **type**: [`definitions.region`](#definitionsregion).
- **required**: `false`
- **description**: The entity's region.
- **usage**:<br><br>
    ```yaml
    authors:
      - name: "The Research Software Project"
        region: Renfrewshire
    ```

### `definitions.entity.tel`

- **type**: [`definitions.tel`](#definitionstel).
- **required**: `false`
- **description**: The entity's telephone number.
- **usage**:<br><br>
    ```yaml
    authors:
      - name: "The Research Software Project"
        tel: +12-345-6789098
    ```

### `definitions.entity.website`

- **type**: [`definitions.url`](#definitionsurl).
- **required**: `false`
- **description**: The entity's website.
- **usage**:<br><br>
    ```yaml
    authors:
      - name: "The Research Software Project"
        website: "https://research-software-project.org"
    ```

### `definitions.fax`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: N/A
- **description**: A fax number.
- **usage**:<br><br>
    ```yaml
    authors:
      - name: "The Research Software Project"
        fax: +12-3456-7890
    ```
    ```yaml
    authors:
      - family-names: Druskat
        fax: +12-3456-7890
        given-names: Stephan
    ```

### `definitions.identifier`

- **type**: One of the following `object` types (click to expand/collapse):<br><br>
    1. <details>
          <summary>DOI</summary>
          <br>
          <ul>
            <li>
              <code>type</code>:
              <ul>
                <li><strong>type</strong>: <code>enum</code> with singular value <code>doi</code></li>
                <li><strong>required</strong>: <code>true</code></li>
                <li><strong>description</strong>: The type of identifier.</li>
              </ul>
            </li>
            <li>
              <code>value</code>:
              <ul>
                <li><strong>type</strong>: <a href="#definitionsdoi"><code>definitions.doi</code></a></li>
                <li><strong>required</strong>: <code>true</code></li>
                <li><strong>description</strong>: The value of the DOI, e.g. <code>10.5281/zenodo.1003149</code></li>
              </ul>
            </li>
            <li>
              <code>description</code>:
              <ul>
                <li><strong>type</strong>: <a href="#definitionsidentifier-description"><code>definitions.identifier-description</code></a></li>
                <li><strong>required</strong>: <code>false</code></li>
                <li><strong>description</strong>: The description of the DOI, e.g. <code>This is the DOI for version 0.11.4.</code></li>
              </ul>
            </li>
          </ul>
        </details>
    1. <details>
          <summary>URL</summary>
          <br>
          <ul>
            <li>
              <code>type</code>:
              <ul>
                <li><strong>type</strong>: <code>enum</code> with singular value <code>url</code></li>
                <li><strong>required</strong>: <code>true</code></li>
                <li><strong>description</strong>: The type of identifier.</li>
              </ul>
            </li>
            <li>
              <code>value</code>:
              <ul>
                <li><strong>type</strong>: <a href="#definitionsurl"><code>definitions.url</code></a></li>
                <li><strong>required</strong>: <code>true</code></li>
                <li><strong>description</strong>: The value of the URL, e.g. <code>https://github.com/citation-file-format/citation-file-format</code>.</li>
              </ul>
            </li>
            <li>
              <code>description</code>:
              <ul>
                <li><strong>type</strong>: <a href="#definitionsidentifier-description"><code>definitions.identifier-description</code></a></li>
                <li><strong>required</strong>: <code>false</code></li>
                <li><strong>description</strong>: The description of the URL, e.g. <code>The homepage for the project</code>.</li>
              </ul>
            </li>
          </ul>
        </details>
    1. <details>
          <summary>Software Heritage identifier</summary>
          <br>
          <ul>
            <li>
              <code>type</code>:
              <ul>
                <li><strong>type</strong>: <code>enum</code> with singular value <code>swh</code></li>
                <li><strong>required</strong>: <code>true</code></li>
                <li><strong>description</strong>: The type of identifier.</li>
              </ul>
            </li>
            <li>
              <code>value</code>:
              <ul>
                <li><strong>type</strong>: <a href="#definitionsswh-identifier"><code>definitions.swh-identifier</code></a></li>
                <li><strong>required</strong>: <code>true</code></li>
                <li><strong>description</strong>: The value of the Software Heritage identifier, e.g. <code>swh:1:dir:bc286860f423ea7ced246ba7458eef4b4541cf2d</code>.</li>
              </ul>
            </li>
            <li>
              <code>description</code>:
              <ul>
                <li><strong>type</strong>: <a href="#definitionsidentifier-description"><code>definitions.identifier-description</code></a></li>
                <li><strong>required</strong>: <code>false</code></li>
                <li><strong>description</strong>: The description of the Software Heritage identifier, e.g. <code>The directory object of the repository as stored on Software Heritage.</code>.</li>
              </ul>
            </li>
          </ul>
        </details>
    1. <details>
          <summary>Other</summary>
          <br>
          <ul>
            <li>
              <code>type</code>:
              <ul>
                <li><strong>type</strong>: <code>enum</code> with singular value <code>other</code></li>
                <li><strong>required</strong>: <code>true</code></li>
                <li><strong>description</strong>: The type of identifier.</li>
              </ul>
            </li>
            <li>
              <code>value</code>:
              <ul>
                <li><strong>type</strong>: <a href="#yaml-strings">Nonempty <code>string</code></a>.</li>
                <li><strong>required</strong>: <code>true</code></li>
                <li><strong>description</strong>: The value of the identifier, e.g. <code>arXiv:2103.06681</code>.</li>
              </ul>
            </li>
            <li>
              <code>description</code>:
              <ul>
                <li><strong>type</strong>: <a href="#definitionsidentifier-description"><code>definitions.identifier-description</code></a></li>
                <li><strong>required</strong>: <code>false</code></li>
                <li><strong>description</strong>: The description of the identifier, e.g. <code>The ArXiv preprint of the paper.</code>.</li>
              </ul>
            </li>
          </ul>
        </details>
- **required**: N/A
- **description**: An identifier.
- **usage**:<br><br>
    ```yaml
    identifiers:
      - type: doi
        value: 10.5281/zenodo.1003149
        description: "The concept DOI of the work."
    ```
    ```yaml
    identifiers:
      - type: doi
        value: 10.5281/zenodo.4813122
        description: "The versioned DOI for version 1.1.0 of the work."
    ```
    ```yaml
    identifiers:
      - type: doi
        value: 10.5281/zenodo.1003149
        description: "The concept DOI of the work."
      - type: doi
        value: 10.5281/zenodo.4813122
        description: "The versioned DOI for version 1.1.0 of the work."
    ```
    ```yaml
    identifiers:
      - type: doi
        value: 10.5281/zenodo.1003149
        description: "The concept DOI of the work."
      - type: doi
        value: 10.5281/zenodo.4813122
        description: "The versioned DOI for version 1.1.0 of the work."
      - type: swh
        value: "swh:1:dir:bc286860f423ea7ced246ba7458eef4b4541cf2d"
        description: "The Software Heritage identifier for version 1.1.0 of the work."
      - type: url
        value: "https://github.com/citation-file-format/citation-file-format/releases/tag/1.1.0"
        description: "The GitHub release URL of tag 1.1.0."
      - type: url
        value: "https://github.com/citation-file-format/citation-file-format/tree/16192bf05e99bcb35d5c3e085047807b5720fafc"
        description: "The GitHub release URL of the commit tagged with 1.1.0."
    ```
    ```yaml
    preferred-citation:
      identifiers:
        - type: other
          value: "arXiv:2103.06681"
          description: The ArXiv preprint of the paper
    ```

### `definitions.identifier-description`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: N/A
- **description**: A description for a specific identifier value.
- **usage**:<br><br>
    ```yaml
    doi: 10.5281/zenodo.4813121
    identifiers:
      - type: doi
        value: 10.5281/zenodo.4813122
        description: "The version DOI for this version, which has a relation childOf with the concept DOI specified in the doi key in the root of this file."
      - type: other
        value: "ar:1234/5678.ABCD"
        description: "The identifier provided by Archival Repository, which points to this version of the software."
    ```

### `definitions.license`

- **type**: (Array of) [`definitions.license-enum`](#definitions.license-enum) objects.
- **required**: N/A
- **description**: The [SPDX license identifier(s)](https://spdx.dev/ids/) for the license(s) under which a work is made available. When there are multiple licenses, it is assumed their relationship is OR, not AND.
- **usage**:<br><br>
    ```yaml
    license: Apache-2.0
    ```
    ```yaml
    license:
      - Apache-2.0
      - MIT
    ```

### `definitions.license-enum`

- **type**: `enum` with values:
    - `0BSD`
    - `AAL`
    - `Abstyles`
    - `Adobe-2006`
    - `Adobe-Glyph`
    - `ADSL`
    - `AFL-1.1`
    - `AFL-1.2`
    - `AFL-2.0`
    - `AFL-2.1`
    - `AFL-3.0`
    - `Afmparse`
    - `AGPL-1.0`
    - `AGPL-1.0-only`
    - `AGPL-1.0-or-later`
    - `AGPL-3.0`
    - `AGPL-3.0-only`
    - `AGPL-3.0-or-later`
    - `Aladdin`
    - `AMDPLPA`
    - `AML`
    - `AMPAS`
    - `ANTLR-PD`
    - `ANTLR-PD-fallback`
    - `Apache-1.0`
    - `Apache-1.1`
    - `Apache-2.0`
    - `APAFML`
    - `APL-1.0`
    - `APSL-1.0`
    - `APSL-1.1`
    - `APSL-1.2`
    - `APSL-2.0`
    - `Artistic-1.0`
    - `Artistic-1.0-cl8`
    - `Artistic-1.0-Perl`
    - `Artistic-2.0`
    - `Bahyph`
    - `Barr`
    - `Beerware`
    - `BitTorrent-1.0`
    - `BitTorrent-1.1`
    - `blessing`
    - `BlueOak-1.0.0`
    - `Borceux`
    - `BSD-1-Clause`
    - `BSD-2-Clause`
    - `BSD-2-Clause-FreeBSD`
    - `BSD-2-Clause-NetBSD`
    - `BSD-2-Clause-Patent`
    - `BSD-2-Clause-Views`
    - `BSD-3-Clause`
    - `BSD-3-Clause-Attribution`
    - `BSD-3-Clause-Clear`
    - `BSD-3-Clause-LBNL`
    - `BSD-3-Clause-Modification`
    - `BSD-3-Clause-No-Nuclear-License`
    - `BSD-3-Clause-No-Nuclear-License-2014`
    - `BSD-3-Clause-No-Nuclear-Warranty`
    - `BSD-3-Clause-Open-MPI`
    - `BSD-4-Clause`
    - `BSD-4-Clause-Shortened`
    - `BSD-4-Clause-UC`
    - `BSD-Protection`
    - `BSD-Source-Code`
    - `BSL-1.0`
    - `BUSL-1.1`
    - `bzip2-1.0.5`
    - `bzip2-1.0.6`
    - `C-UDA-1.0`
    - `CAL-1.0`
    - `CAL-1.0-Combined-Work-Exception`
    - `Caldera`
    - `CATOSL-1.1`
    - `CC-BY-1.0`
    - `CC-BY-2.0`
    - `CC-BY-2.5`
    - `CC-BY-3.0`
    - `CC-BY-3.0-AT`
    - `CC-BY-3.0-US`
    - `CC-BY-4.0`
    - `CC-BY-NC-1.0`
    - `CC-BY-NC-2.0`
    - `CC-BY-NC-2.5`
    - `CC-BY-NC-3.0`
    - `CC-BY-NC-4.0`
    - `CC-BY-NC-ND-1.0`
    - `CC-BY-NC-ND-2.0`
    - `CC-BY-NC-ND-2.5`
    - `CC-BY-NC-ND-3.0`
    - `CC-BY-NC-ND-3.0-IGO`
    - `CC-BY-NC-ND-4.0`
    - `CC-BY-NC-SA-1.0`
    - `CC-BY-NC-SA-2.0`
    - `CC-BY-NC-SA-2.5`
    - `CC-BY-NC-SA-3.0`
    - `CC-BY-NC-SA-4.0`
    - `CC-BY-ND-1.0`
    - `CC-BY-ND-2.0`
    - `CC-BY-ND-2.5`
    - `CC-BY-ND-3.0`
    - `CC-BY-ND-4.0`
    - `CC-BY-SA-1.0`
    - `CC-BY-SA-2.0`
    - `CC-BY-SA-2.0-UK`
    - `CC-BY-SA-2.1-JP`
    - `CC-BY-SA-2.5`
    - `CC-BY-SA-3.0`
    - `CC-BY-SA-3.0-AT`
    - `CC-BY-SA-4.0`
    - `CC-PDDC`
    - `CC0-1.0`
    - `CDDL-1.0`
    - `CDDL-1.1`
    - `CDL-1.0`
    - `CDLA-Permissive-1.0`
    - `CDLA-Sharing-1.0`
    - `CECILL-1.0`
    - `CECILL-1.1`
    - `CECILL-2.0`
    - `CECILL-2.1`
    - `CECILL-B`
    - `CECILL-C`
    - `CERN-OHL-1.1`
    - `CERN-OHL-1.2`
    - `CERN-OHL-P-2.0`
    - `CERN-OHL-S-2.0`
    - `CERN-OHL-W-2.0`
    - `ClArtistic`
    - `CNRI-Jython`
    - `CNRI-Python`
    - `CNRI-Python-GPL-Compatible`
    - `Condor-1.1`
    - `copyleft-next-0.3.0`
    - `copyleft-next-0.3.1`
    - `CPAL-1.0`
    - `CPL-1.0`
    - `CPOL-1.02`
    - `Crossword`
    - `CrystalStacker`
    - `CUA-OPL-1.0`
    - `Cube`
    - `curl`
    - `D-FSL-1.0`
    - `diffmark`
    - `DOC`
    - `Dotseqn`
    - `DRL-1.0`
    - `DSDP`
    - `dvipdfm`
    - `ECL-1.0`
    - `ECL-2.0`
    - `eCos-2.0`
    - `EFL-1.0`
    - `EFL-2.0`
    - `eGenix`
    - `Entessa`
    - `EPICS`
    - `EPL-1.0`
    - `EPL-2.0`
    - `ErlPL-1.1`
    - `etalab-2.0`
    - `EUDatagrid`
    - `EUPL-1.0`
    - `EUPL-1.1`
    - `EUPL-1.2`
    - `Eurosym`
    - `Fair`
    - `Frameworx-1.0`
    - `FreeBSD-DOC`
    - `FreeImage`
    - `FSFAP`
    - `FSFUL`
    - `FSFULLR`
    - `FTL`
    - `GD`
    - `GFDL-1.1`
    - `GFDL-1.1-invariants-only`
    - `GFDL-1.1-invariants-or-later`
    - `GFDL-1.1-no-invariants-only`
    - `GFDL-1.1-no-invariants-or-later`
    - `GFDL-1.1-only`
    - `GFDL-1.1-or-later`
    - `GFDL-1.2`
    - `GFDL-1.2-invariants-only`
    - `GFDL-1.2-invariants-or-later`
    - `GFDL-1.2-no-invariants-only`
    - `GFDL-1.2-no-invariants-or-later`
    - `GFDL-1.2-only`
    - `GFDL-1.2-or-later`
    - `GFDL-1.3`
    - `GFDL-1.3-invariants-only`
    - `GFDL-1.3-invariants-or-later`
    - `GFDL-1.3-no-invariants-only`
    - `GFDL-1.3-no-invariants-or-later`
    - `GFDL-1.3-only`
    - `GFDL-1.3-or-later`
    - `Giftware`
    - `GL2PS`
    - `Glide`
    - `Glulxe`
    - `GLWTPL`
    - `gnuplot`
    - `GPL-1.0`
    - `GPL-1.0-only`
    - `GPL-1.0-or-later`
    - `GPL-1.0+`
    - `GPL-2.0`
    - `GPL-2.0-only`
    - `GPL-2.0-or-later`
    - `GPL-2.0-with-autoconf-exception`
    - `GPL-2.0-with-bison-exception`
    - `GPL-2.0-with-classpath-exception`
    - `GPL-2.0-with-font-exception`
    - `GPL-2.0-with-GCC-exception`
    - `GPL-2.0+`
    - `GPL-3.0`
    - `GPL-3.0-only`
    - `GPL-3.0-or-later`
    - `GPL-3.0-with-autoconf-exception`
    - `GPL-3.0-with-GCC-exception`
    - `GPL-3.0+`
    - `gSOAP-1.3b`
    - `HaskellReport`
    - `Hippocratic-2.1`
    - `HPND`
    - `HPND-sell-variant`
    - `HTMLTIDY`
    - `IBM-pibs`
    - `ICU`
    - `IJG`
    - `ImageMagick`
    - `iMatix`
    - `Imlib2`
    - `Info-ZIP`
    - `Intel`
    - `Intel-ACPI`
    - `Interbase-1.0`
    - `IPA`
    - `IPL-1.0`
    - `ISC`
    - `JasPer-2.0`
    - `JPNIC`
    - `JSON`
    - `LAL-1.2`
    - `LAL-1.3`
    - `Latex2e`
    - `Leptonica`
    - `LGPL-2.0`
    - `LGPL-2.0-only`
    - `LGPL-2.0-or-later`
    - `LGPL-2.0+`
    - `LGPL-2.1`
    - `LGPL-2.1-only`
    - `LGPL-2.1-or-later`
    - `LGPL-2.1+`
    - `LGPL-3.0`
    - `LGPL-3.0-only`
    - `LGPL-3.0-or-later`
    - `LGPL-3.0+`
    - `LGPLLR`
    - `Libpng`
    - `libpng-2.0`
    - `libselinux-1.0`
    - `libtiff`
    - `LiLiQ-P-1.1`
    - `LiLiQ-R-1.1`
    - `LiLiQ-Rplus-1.1`
    - `Linux-OpenIB`
    - `LPL-1.0`
    - `LPL-1.02`
    - `LPPL-1.0`
    - `LPPL-1.1`
    - `LPPL-1.2`
    - `LPPL-1.3a`
    - `LPPL-1.3c`
    - `MakeIndex`
    - `MirOS`
    - `MIT`
    - `MIT-0`
    - `MIT-advertising`
    - `MIT-CMU`
    - `MIT-enna`
    - `MIT-feh`
    - `MIT-Modern-Variant`
    - `MIT-open-group`
    - `MITNFA`
    - `Motosoto`
    - `mpich2`
    - `MPL-1.0`
    - `MPL-1.1`
    - `MPL-2.0`
    - `MPL-2.0-no-copyleft-exception`
    - `MS-PL`
    - `MS-RL`
    - `MTLL`
    - `MulanPSL-1.0`
    - `MulanPSL-2.0`
    - `Multics`
    - `Mup`
    - `NAIST-2003`
    - `NASA-1.3`
    - `Naumen`
    - `NBPL-1.0`
    - `NCGL-UK-2.0`
    - `NCSA`
    - `Net-SNMP`
    - `NetCDF`
    - `Newsletr`
    - `NGPL`
    - `NIST-PD`
    - `NIST-PD-fallback`
    - `NLOD-1.0`
    - `NLPL`
    - `Nokia`
    - `NOSL`
    - `Noweb`
    - `NPL-1.0`
    - `NPL-1.1`
    - `NPOSL-3.0`
    - `NRL`
    - `NTP`
    - `NTP-0`
    - `Nunit`
    - `O-UDA-1.0`
    - `OCCT-PL`
    - `OCLC-2.0`
    - `ODbL-1.0`
    - `ODC-By-1.0`
    - `OFL-1.0`
    - `OFL-1.0-no-RFN`
    - `OFL-1.0-RFN`
    - `OFL-1.1`
    - `OFL-1.1-no-RFN`
    - `OFL-1.1-RFN`
    - `OGC-1.0`
    - `OGDL-Taiwan-1.0`
    - `OGL-Canada-2.0`
    - `OGL-UK-1.0`
    - `OGL-UK-2.0`
    - `OGL-UK-3.0`
    - `OGTSL`
    - `OLDAP-1.1`
    - `OLDAP-1.2`
    - `OLDAP-1.3`
    - `OLDAP-1.4`
    - `OLDAP-2.0`
    - `OLDAP-2.0.1`
    - `OLDAP-2.1`
    - `OLDAP-2.2`
    - `OLDAP-2.2.1`
    - `OLDAP-2.2.2`
    - `OLDAP-2.3`
    - `OLDAP-2.4`
    - `OLDAP-2.5`
    - `OLDAP-2.6`
    - `OLDAP-2.7`
    - `OLDAP-2.8`
    - `OML`
    - `OpenSSL`
    - `OPL-1.0`
    - `OSET-PL-2.1`
    - `OSL-1.0`
    - `OSL-1.1`
    - `OSL-2.0`
    - `OSL-2.1`
    - `OSL-3.0`
    - `Parity-6.0.0`
    - `Parity-7.0.0`
    - `PDDL-1.0`
    - `PHP-3.0`
    - `PHP-3.01`
    - `Plexus`
    - `PolyForm-Noncommercial-1.0.0`
    - `PolyForm-Small-Business-1.0.0`
    - `PostgreSQL`
    - `PSF-2.0`
    - `psfrag`
    - `psutils`
    - `Python-2.0`
    - `Qhull`
    - `QPL-1.0`
    - `Rdisc`
    - `RHeCos-1.1`
    - `RPL-1.1`
    - `RPL-1.5`
    - `RPSL-1.0`
    - `RSA-MD`
    - `RSCPL`
    - `Ruby`
    - `SAX-PD`
    - `Saxpath`
    - `SCEA`
    - `Sendmail`
    - `Sendmail-8.23`
    - `SGI-B-1.0`
    - `SGI-B-1.1`
    - `SGI-B-2.0`
    - `SHL-0.5`
    - `SHL-0.51`
    - `SimPL-2.0`
    - `SISSL`
    - `SISSL-1.2`
    - `Sleepycat`
    - `SMLNJ`
    - `SMPPL`
    - `SNIA`
    - `Spencer-86`
    - `Spencer-94`
    - `Spencer-99`
    - `SPL-1.0`
    - `SSH-OpenSSH`
    - `SSH-short`
    - `SSPL-1.0`
    - `StandardML-NJ`
    - `SugarCRM-1.1.3`
    - `SWL`
    - `TAPR-OHL-1.0`
    - `TCL`
    - `TCP-wrappers`
    - `TMate`
    - `TORQUE-1.1`
    - `TOSL`
    - `TU-Berlin-1.0`
    - `TU-Berlin-2.0`
    - `UCL-1.0`
    - `Unicode-DFS-2015`
    - `Unicode-DFS-2016`
    - `Unicode-TOU`
    - `Unlicense`
    - `UPL-1.0`
    - `Vim`
    - `VOSTROM`
    - `VSL-1.0`
    - `W3C`
    - `W3C-19980720`
    - `W3C-20150513`
    - `Watcom-1.0`
    - `Wsuipa`
    - `WTFPL`
    - `wxWindows`
    - `X11`
    - `Xerox`
    - `XFree86-1.1`
    - `xinetd`
    - `Xnet`
    - `xpp`
    - `XSkat`
    - `YPL-1.0`
    - `YPL-1.1`
    - `Zed`
    - `Zend-2.0`
    - `Zimbra-1.3`
    - `Zimbra-1.4`
    - `Zlib`
    - `zlib-acknowledgement`
    - `ZPL-1.1`
    - `ZPL-2.0`
    - `ZPL-2.1`
- **required**: N/A
- **description**: [SPDX identifier](https://spdx.dev/ids/) for the license under which a work is made available. The list of identifiers originates from https://github.com/spdx/license-list-data/blob/bd8e963a41b13524b2ccb67f9335d2dd397c378e/json/licenses.json.
- **usage**:<br><br>
    ```yaml
    license: Apache-2.0
    ```
    ```yaml
    license:
      - Apache-2.0
      - MIT
    ```

### `definitions.orcid`

- **type**: `uri` with pattern [`https://orcid\.org/[0-9]{4}-[0-9]{4}-[0-9]{4}-[0-9]{3}[0-9X]{1}`](https://regex101.com/library/wvvVYE)
- **required**: N/A
- **description**: An [ORCID](https://orcid.org) identifier.
- **usage**:<br><br>
    ```yaml
    authors:
      - family-names: Druskat
        given-names: Stephan
        orcid: "https://orcid.org/0000-0003-4925-7248"
    ```

### `definitions.person`

- **type**: `object` with the following keys:
    - [`address`](#definitionspersonaddress)
    - [`affiliation`](#definitionspersonaffiliation)
    - [`alias`](#definitionspersonalias)
    - [`city`](#definitionspersoncity)
    - [`country`](#definitionspersoncountry)
    - [`email`](#definitionspersonemail)
    - [`family-names`](#definitionspersonfamily-names)
    - [`fax`](#definitionspersonfax)
    - [`given-names`](#definitionspersongiven-names)
    - [`name-particle`](#definitionspersonname-particle)
    - [`name-suffix`](#definitionspersonname-suffix)
    - [`orcid`](#definitionspersonorcid)
    - [`post-code`](#definitionspersonpost-code)
    - [`region`](#definitionspersonregion)
    - [`tel`](#definitionspersontel)
    - [`website`](#definitionspersonwebsite)
- **required**: N/A
- **description**: A person.
The Citation File Format aims to implement a culturally neutral model for personal names, according to the [suggestions on splitting personal names by the W3C](https://www.w3.org/International/questions/qa-personal-names) and the implementation of personal name splitting in BibTeX ([Hufflen, 2006](https://www.tug.org/TUGboat/tb27-2/tb87hufflen.pdf)). To this end, the Citation File Format provides four generic keys to specify personal names:

    1. Values for `family-names` specify family names, including combinations of given and patronymic forms, such as *Guðmundsdóttir* or *bin Osman*; double names with or without hyphen, such as *Leutheusser-Schnarrenberger* or *Sánchez Vicario*. It can potentially also specify names that include prepositions or (nobiliary) particles, especially if they occur in between family names such as in Spanish- or Portuguese-origin names, such as *Fernández de Córdoba*.
    2. Values for `given-names` specify given and any other names.
    3. Values for `name-particle` specify [nobiliary particles](https://en.wikipedia.org/wiki/Nobiliary_particle) and prepositions, such as in Ludwig *van* Beethoven or Rafael *van der* Vaart.
    4. Values for `name-suffix` specify suffixes such as *Jr.* or *III* (as in Frank Edwin Wright *III*).

Note that these keys may still not be optimal for, e.g., Icelandic names which do not have the concept of family names, or Chinese generation names, but represent a best effort.
- **usage**:<br><br>
    ```yaml
    authors:
      - address: "742 Evergreen Terrace"
        affiliation: "German Aerospace Center (DLR)"
        alias: sdruskat
        city: Berlin
        country: DE
        email: sdruskat@research-project.org
        family-names: Druskat
        fax: +12-3456-7890
        given-names: Stephan
        name-particle: von
        name-suffix: III
        orcid: "https://orcid.org/0000-0003-4925-7248"
        post-code: 90210
        region: Renfrewshire
        tel: +12-345-6789098
        website: "https://research-project.org"
    ```

### `definitions.person.address`

- **type**: [`definitions.address`](#definitionsaddress)
- **required**: `false`
- **description**: The person's address.
- **usage**:<br><br>
    ```yaml
    authors:
      - family-names: Druskat
        given-names: Stephan
        address: "742 Evergreen Terrace"
    ```

### `definitions.person.affiliation`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The person's affiliation.
- **usage**:<br><br>
    ```yaml
    authors:
      - family-names: Druskat
        given-names: Stephan
        affiliation: "German Aerospace Center (DLR)"
    ```

### `definitions.person.alias`

- **type**: [`definitions.alias`](#definitionsalias)
- **required**: `false`
- **description**: The person's alias.
- **usage**:<br><br>
    ```yaml
    authors:
      - family-names: Druskat
        given-names: Stephan
        alias: sdruskat
    ```

### `definitions.person.city`

- **type**: [`definitions.city`](#definitionscity)
- **required**: `false`
- **description**: The person's city.
- **usage**:<br><br>
    ```yaml
    authors:
      - family-names: Druskat
        given-names: Stephan
        city: Berlin
    ```

### `definitions.person.country`

- **type**: [`definitions.country`](#definitioncountry)
- **required**: `false`
- **description**: The person's country.
- **usage**:<br><br>
    ```yaml
    authors:
      - family-names: Druskat
        given-names: Stephan
        country: DE
    ```

### `definitions.person.email`

- **type**: [`definitions.email`](#definitionsemail)
- **required**: `false`
- **description**: The person's email address.
- **usage**:<br><br>
    ```yaml
    authors:
      - family-names: Druskat
        given-names: Stephan
        email: mail@research-project.org
    ```

### `definitions.person.family-names`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The person's family names.
- **usage**:<br><br>
    ```yaml
    authors:
      - family-names: Druskat
        given-names: Stephan
    ```

### `definitions.person.fax`

- **type**: [`definitions.fax`](#definitionsfax)
- **required**: `false`
- **description**: The person's fax number.
- **usage**:<br><br>
    ```yaml
    authors:
      - family-names: Druskat
        given-names: Stephan
        fax: +12-3456-7890
    ```

### `definitions.person.given-names`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The person's given names.
- **usage**:<br><br>
    ```yaml
    authors:
      - family-names: Druskat
        given-names: Stephan
    ```

### `definitions.person.name-particle`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The person's name particle, e.g., a [nobiliary particle](https://en.wikipedia.org/wiki/Nobiliary_particle) or a [preposition] meaning 'of' or 'from' (for example 'von' in 'Alexander von Humboldt').
- **usage**:<br><br>
    ```yaml
    authors:
      - family-names: Humboldt
        given-names: Alexander
        name-particle: von
    ```

### `definitions.person.name-suffix`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The person's [name suffix](https://en.wikipedia.org/wiki/Suffix_(name)), e.g. 'Jr.' for Sammy Davis Jr. or 'III' for Frank Edwin Wright III.
- **usage**:<br><br>
    ```yaml
    authors:
      - family-names: Davis
        given-names: Sammy
        name-suffix: Jr.
    ```

### `definitions.person.orcid`

- **type**: [`definitions.orcid`](#definitionsorcid)
- **required**: `false`
- **description**: The person's [ORCID](https://orcid.org) identifier.
- **usage**:<br><br>
    ```yaml
    authors:
      - family-names: Druskat
        given-names: Stephan
        orcid: "https://orcid.org/0000-0003-4925-7248"
    ```

### `definitions.person.post-code`

- **type**: [`definitions.post-code`](#definitionspost-code)
- **required**: `false`
- **description**: The person's post code.
- **usage**:<br><br>
    ```yaml
    authors:
      - family-names: Druskat
        given-names: Stephan
        post-code: 90210
    ```
    ```yaml
    authors:
      - family-names: Druskat
        given-names: Stephan
        post-code: "90210"
    ```

### `definitions.person.region`

- **type**: [`definitions.region`](#definitionsregion)
- **required**: `false`
- **description**: The person's region.
- **usage**:<br><br>
    ```yaml
    authors:
      - family-names: Druskat
        given-names: Stephan
        region: Renfrewshire
    ```

### `definitions.person.tel`

- **type**: [`definitions.tel`](#definitionstel)
- **required**: `false`
- **description**: The person's telephone number.
- **usage**:<br><br>
    ```yaml
    authors:
      - family-names: Druskat
        given-names: Stephan
        tel: +12-345-6789098
    ```

### `definitions.person.website`

- **type**: [`definitions.url`](#definitionsurl)
- **required**: `false`
- **description**: The person's website.
- **usage**:<br><br>
    ```yaml
    authors:
      - family-names: Druskat
        given-names: Stephan
        website: "https://research-project.org"
    ```

### `definitions.post-code`

- **type**: `string` or `number`
- **required**: N/A
- **description**: A post code.
- **usage**:<br><br>
    ```yaml
    authors:
      - name: "The Research Software Project"
        post-code: "G42 2PN"
    ```
    ```yaml
    authors:
      - name: "The Research Software Project"
        post-code: 12053
    ```

### `definitions.reference`

- **type**: `object` with the following keys:
    - [`abbreviation`](#definitionsreferenceabbreviation)
    - [`abstract`](#definitionsreferenceabstract)
    - [`authors`](#definitionsreferenceauthors)
    - [`collection-doi`](#definitionsreferencecollection-doi)
    - [`collection-title`](#definitionsreferencecollection-title)
    - [`collection-type`](#definitionsreferencecollection-type)
    - [`commit`](#definitionsreferencecommit)
    - [`conference`](#definitionsreferenceconference)
    - [`contact`](#definitionsreferencecontact)
    - [`copyright`](#definitionsreferencecopyright)
    - [`data-type`](#definitionsreferencedata-type)
    - [`database-provider`](#definitionsreferencedatabase-provider)
    - [`database`](#definitionsreferencedatabase)
    - [`date-accessed`](#definitionsreferencedate-accessed)
    - [`date-downloaded`](#definitionsreferencedate-downloaded)
    - [`date-published`](#definitionsreferencedate-published)
    - [`date-released`](#definitionsreferencedate-released)
    - [`department`](#definitionsreferencedepartment)
    - [`doi`](#definitionsreferencedoi)
    - [`edition`](#definitionsreferenceedition)
    - [`editors`](#definitionsreferenceeditors)
    - [`editors-series`](#definitionsreferenceeditors-series)
    - [`end`](#definitionsreferenceend)
    - [`entry`](#definitionsreferenceentry)
    - [`filename`](#definitionsreferencefilename)
    - [`format`](#definitionsreferenceformat)
    - [`identifiers`](#definitionsreferenceidentifiers)
    - [`institution`](#definitionsreferenceinstitution)
    - [`isbn`](#definitionsreferenceisbn)
    - [`issn`](#definitionsreferenceissn)
    - [`issue`](#definitionsreferenceissue)
    - [`issue-date`](#definitionsreferenceissue-date)
    - [`issue-title`](#definitionsreferenceissue-title)
    - [`journal`](#definitionsreferencejournal)
    - [`keywords`](#definitionsreferencekeywords)
    - [`languages`](#definitionsreferencelanguages)
    - [`license`](#definitionsreferencelicense)
    - [`license-url`](#definitionsreferencelicense-url)
    - [`loc-end`](#definitionsreferenceloc-end)
    - [`loc-start`](#definitionsreferenceloc-start)
    - [`location`](#definitionsreferencelocation)
    - [`medium`](#definitionsreferencemedium)
    - [`month`](#definitionsreferencemonth)
    - [`nihmsid`](#definitionsreferencenihmsid)
    - [`notes`](#definitionsreferencenotes)
    - [`number`](#definitionsreferencenumber)
    - [`number-volumes`](#definitionsreferencenumber-volumes)
    - [`pages`](#definitionsreferencepages)
    - [`patent-states`](#definitionsreferencepatent-states)
    - [`pmcid`](#definitionsreferencepmcid)
    - [`publisher`](#definitionsreferencepublisher)
    - [`recipients`](#definitionsreferencerecipients)
    - [`repository`](#definitionsreferencerepository)
    - [`repository-artifact`](#definitionsreferencerepository-artifact)
    - [`repository-code`](#definitionsreferencerepository-code)
    - [`scope`](#definitionsreferencescope)
    - [`section`](#definitionsreferencesection)
    - [`senders`](#definitionsreferencesenders)
    - [`start`](#definitionsreferencestart)
    - [`status`](#definitionsreferencestatus)
    - [`term`](#definitionsreferenceterm)
    - [`thesis-type`](#definitionsreferencethesis-type)
    - [`title`](#definitionsreferencetitle)
    - [`translators`](#definitionsreferencetranslators)
    - [`type`](#definitionsreferencetype)
    - [`url`](#definitionsreferenceurl)
    - [`version`](#definitionsreferenceversion)
    - [`volume`](#definitionsreferencevolume)
    - [`volume-title`](#definitionsreferencevolume-title)
    - [`year`](#definitionsreferenceyear)
    - [`year-original`](#definitionsreferenceyear-original)
- **required**: N/A
- **description**: A reference.
- **usage**:<br><br>
    ```yaml
    references:
      - authors:
          - family-names: Smith
            given-names: "A. M."
          - family-names: Katz
            given-names: "D. S."
          - family-names: Niemeyer
            given-names: "K. E."
          - name: "FORCE11 Software Citation Working Group"
        doi: 10.7717/peerj-cs.86
        journal: "PeerJ Computer Science"
        month: 9
        start: e86
        title: "Software citation principles"
        type: article
        volume: 2
        year: 2016
    ```
    ```yaml
    references:
      - abbreviation: NP
        abstract: "This is a non-sensical example reference to show how the many different keys in a reference object can be used."
        authors:
          - name: "The Research Software Project"
        collection-doi: 10.5281/zenodo.1003149
        collection-title: "Proceedings of the Research Conference 2021"
        collection-type: proceedings
        commit: 16192bf05e99bcb35d5c3e085047807b5720fafc
        conference:
          name: "Research Conference 2021"
        contact:
          - name: "The RC21 Organizing Committee"
        copyright: "© 2021 The Research Software Project team"
        data-type: YAML
        database: "Research Database"
        database-provider:
          name: "Research Databases Ltd."
        date-accessed: "2021-07-27"
        date-downloaded: "2021-07-27"
        date-published: "2021-07-26"
        date-released: "2021-07-26"
        department: "Department of Hard Science Fiction"
        doi: 10.5281/zenodo.4813122
        edition: "2nd abridged edition"
        editors:
          - family-names: Inchief
            given-names: Editor
          - name: "The RCProc Editorial Team"
        editors-series:
          - family-names: Editor
            given-names: Series
          - name: "The Series Editors"
        end: 42
        entry: "Citation <n. 1>"
        filename: CITATION.cff
        format: "Citation File Format"
        identifiers:
          - description: "The concept DOI of the work."
            type: doi
            value: 10.5281/zenodo.1003149
          - description: "The versioned DOI for version 1.1.0 of the work."
            type: doi
            value: 10.5281/zenodo.4813122
          - description: "The Software Heritage identifier for version 1.1.0 of the work."
            type: swh
            value: "swh:1:dir:bc286860f423ea7ced246ba7458eef4b4541cf2d"
          - description: "The GitHub release URL of tag 1.1.0."
            type: url
            value: "https://github.com/citation-file-format/citation-file-format/releases/tag/1.1.0"
        institution:
          name: "University of Arcadia"
        isbn: "9781603095075"
        issn: 2475-9066
        issue: 42
        issue-date: "November/December 2021"
        issue-title: "Special Issue: Software Citation"
        journal: "Journal of Open Source Software"
        keywords:
          - "software citation"
          - "citation file format"
          - research
        languages:
          - en
          - haw
        license: Apache-2.0
        license-url: "https://obscure-licenses.com?id=1234"
        loc-end: 42
        loc-start: 21
        location:
          name: "Library of the Unseen University"
        medium: "5¼-inch floppy disk"
        month: 7
        nihmsid: NIHMS236863
        notes: "Excellent reference! TODO Read for thesis."
        number: 12053
        number-volumes: 7
        pages: 78
        patent-states:
          - Canada
        pmcid: PMC3134971
        publisher:
          name: "Open Access Publishing House"
        recipients:
          - name: "Recipient entity of personal communication"
          - family-names: Recipient
            given-names: Communication
        repository: "https://ascl.net/2105.013"
        repository-artifact: "https://search.maven.org/artifact/org.corpus-tools/cff-maven-plugin/0.4.0/maven-plugin"
        repository-code: "https://github.com/citation-file-format/my-research-software"
        scope: "Supplement 2: Additional material"
        section: 7
        senders:
          - name: "Sender entity of personal communication"
          - family-names: Sender
            given-names: Communication
        start: 17
        status: submitted
        term: Citation
        thesis-type: "PhD thesis"
        title: "Towards better software citation"
        translators:
          - name: "Research Translators Ltd."
          - family-names: Lator
            given-names: Trans
        type: conference-paper
        url: "https://citation-file-format.github.io/"
        version: 0.3.12
        volume: 2
        volume-title: "Volume II: How it went on"
        year: 2021
        year-original: 1978
    ```

### `definitions.reference.abbreviation`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The abbreviation of a work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      abbreviation: ABC
      type: generic
    ```
    ```yaml
    references:
      - abbreviation: DEF
        type: generic
    ```

### `definitions.reference.abstract`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The abstract of the work.
    - If the work is a journal paper or other academic work: The abstract of the work.
    - If the work is a film, broadcast or similar: The synopsis of the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      abstract: "This work describes the software or dataset that should be actually cited. etc."
      type: generic
    ```
    ```yaml
    references:
      - abstract: "This work implements an algorithm that we use in our software. etc."
        type: generic
    ```

### `definitions.reference.authors`

- **type**: Array of [`definitions.person`](#definitionsperson) and/or [`definitions.entity`](#definitionsentity) objects.
- **required**: `true`
- **description**: The authors of the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      authors:
        - name: "The Research Software Project team"
        - family-names: Druskat
          given-names: Stephan
      type: generic
    ```
    ```yaml
    references:
      - authors:
          - name: "The Research Software Project team"
          - family-names: Druskat
            given-names: Stephan
        type: generic
    ```

#### How to deal with unknown individual authors?

To enable credit for the individuals that have created a work,
it is good practice to cite the respective individuals as authors.
Sometimes you may not be able to determine the person names of the relevant individuals to create
[`definitions.person`](#definitionsperson) objects for them,
for example when a software you cite does not provide a `CITATION.cff` file.
Then, the next best thing is to refer to those that you could not determine person names for
collectively as a "team" or "project" using the title of the work
in a [`definitions.entity`](#definitionsentity) object:

```yaml
authors:
  - name: "The Research Software project"
```
```yaml
authors:
  - family-names: Spaaks
    given-names: "Jurriaan H."
  - family-names: Druskat
    given-names: Stephan
  - name: "The Research Software team"
```

This still represents the maximum knowledge you have about the authors.
It's also a good starting point to enable users of the citation metadata to determine the correct names themselves in the future,
especially if you also provide a source of information such as a `repository-code` or a `url`.

If the authors of a work are truly anonymous,
you can represent this in the same way:

```yaml
authors:
  - name: anonymous
```

### `definitions.reference.collection-doi`

- **type**: [`definitions.doi`](#definitionsdoi)
- **required**: `false`
- **description**: The [DOI](https://en.wikipedia.org/wiki/Digital_object_identifier) of a collection containing the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      collection-doi: 10.5281/zenodo.1003149
      type: generic
    ```
    ```yaml
    references:
      - collection-doi: 10.5281/zenodo.1003149
        type: generic
    ```

### `definitions.reference.collection-title`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The title of a collection or proceedings.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      collection-title: "Proceedings of the Research Conference 2021"
      type: generic
    ```
    ```yaml
    references:
      - collection-title: "Proceedings of the Research Conference 2021"
        type: generic
    ```

### `definitions.reference.collection-type`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The type of a collection.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      collection-type: proceedings
      type: generic
    ```
    ```yaml
    references:
      - collection-type: proceedings
        type: generic
    ```

### `definitions.reference.commit`

- **type**: [`definitions.commit`](#definitionscommit)
- **required**: `false`
- **description**: The (e.g., Git) commit hash or (e.g., Subversion) revision number of the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      commit: "16192bf05e99bcb35d5c3e085047807b5720fafc"
      type: generic
    ```
    ```yaml
    references:
      - commit: "16192bf05e99bcb35d5c3e085047807b5720fafc"
        type: generic
    ```

### `definitions.reference.conference`

- **type**: [`definitions.entity`](#definitionsentity)
- **required**: `false`
- **description**: The conference where the work was presented.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      conference:
        name: "Research Conference 2021"
      type: generic
    ```
    ```yaml
    references:
      - conference:
          name: "Research Conference 2021"
        type: generic
    ```

### `definitions.reference.contact`

- **type**: Array of [`definitions.person`](#definitionsperson) and/or [`definitions.entity`](#definitionsentity) objects.
- **required**: `false`
- **description**: The contact person, group, company, etc. for a work.
- **usage**:<br><br>
    ```yaml
    references:
      - contact:
        - name: "The RC21 Organizing Committee"
        - family-names: Druskat
          given-names: Stephan
        type: generic
    ```
    ```yaml
    preferred-citation:
      contact:
        - name: "The RC21 Organizing Committee"
        - family-names: Druskat
          given-names: Stephan
      type: generic
    ```

### `definitions.reference.copyright`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The copyright information pertaining to the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      copyright: "© 2021 The Research Software Project team"
      type: generic
    ```
    ```yaml
    references:
      - copyright: "© 2021 The Research Software Project team"
        type: generic
    ```

### `definitions.reference.data-type`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The data type of a data set.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      data-type: YAML
      type: generic
    ```
    ```yaml
    references:
      - data-type: YAML
        type: generic
    ```

### `definitions.reference.database-provider`

- **type**: [`definitions.entity`](#definitionsentity)
- **required**: `false`
- **description**: The provider of the database where a work was accessed/is stored.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      database-provider:
        name: "Research Databases Ltd."
      type: generic
    ```
    ```yaml
    references:
      - database-provider:
          name: "Research Databases Ltd."
        type: generic
    ```

### `definitions.reference.database`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The name of the database where a work was accessed/is stored.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      database: "Research Database"
      type: generic
    ```
    ```yaml
    references:
      - database: "Research Database"
        type: generic
    ```

### `definitions.reference.date-accessed`

- **type**: [`definitions.date`](#definitionsdate)
- **required**: `false`
- **description**: The date the work was accessed.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      date-accessed: "2021-07-27"
      type: generic
    ```
    ```yaml
    references:
      - date-accessed: "2021-07-27"
        type: generic
    ```

### `definitions.reference.date-downloaded`

- **type**: [`definitions.date`](#definitionsdate)
- **required**: `false`
- **description**: The date the work has been downloaded.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      date-downloaded: "2021-07-27"
      type: generic
    ```
    ```yaml
    references:
      - date-downloaded: "2021-07-27"
        type: generic
    ```

### `definitions.reference.date-published`

- **type**: [`definitions.date`](#definitionsdate)
- **required**: `false`
- **description**: The date the work has been published.
- **usage**:<br><br>
    ```yaml
    references:
      - date-published: "2021-07-27"
        type: generic
    ```
    ```yaml
    preferred-citation:
      date-published: "2021-07-27"
      type: generic
    ```

### `definitions.reference.date-released`

- **type**: [`definitions.date`](#definitionsdate)
- **required**: `false`
- **description**: The date the work has been released.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      date-released: "2021-07-27"
      type: generic
    ```
    ```yaml
    references:
      - date-released: "2021-07-27"
        type: generic
    ```

### `definitions.reference.department`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The department where a work has been produced.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      department: "Department of Hard Science Fiction"
      type: generic
    ```
    ```yaml
    references:
      - department: "Department of Hard Science Fiction"
        type: generic
    ```

### `definitions.reference.doi`

- **type**: [`definitions.doi`](#definitionsdoi)
- **required**: `false`
- **description**: The [DOI](https://en.wikipedia.org/wiki/Digital_object_identifier) of the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      doi: 10.5281/zenodo.4813122
      type: generic
    ```
    ```yaml
    references:
      - doi: 10.5281/zenodo.4813122
        type: generic
    ```

### `definitions.reference.edition`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The edition of the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      edition: "2nd abridged edition"
      type: generic
    ```
    ```yaml
    references:
      - edition: "2nd abridged edition"
        type: generic
    ```

### `definitions.reference.editors`

- **type**: Array of [`definitions.person`](#definitionsperson) and/or [`definitions.entity`](#definitionsentity) objects.
- **required**: `false`
- **description**: The editor(s) of a work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      editors:
        - family-names: Inchief
          given-names: Editor
        - name: "The RCProc Editorial Team"
      type: generic
    ```
    ```yaml
    references:
      - editors:
        - family-names: Inchief
          given-names: Editor
        - name: "The RCProc Editorial Team"
        type: generic
    ```

### `definitions.reference.editors-series`

- **type**: Array of [`definitions.person`](#definitionsperson) and/or [`definitions.entity`](#definitionsentity) objects.
- **required**: `false`
- **description**: The editor(s) of a series in which a work has been published.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      editors-series:
        - family-names: Editor
          given-names: Series
        - name: "The Series Editors"
      type: generic
    ```
    ```yaml
    references:
      - editors-series:
        - family-names: Editor
          given-names: Series
        - name: "The Series Editors"
        type: generic
    ```

### `definitions.reference.end`

- **type**: [Nonempty `string`](#yaml-strings) or `integer`
- **required**: `false`
- **description**: The end page of the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      end: "42"
      type: generic
    ```
    ```yaml
    preferred-citation:
      end: 42
      type: generic
    ```
    ```yaml
    references:
      - end: "42"
        type: generic
    ```
    ```yaml
    references:
      - end: 42
        type: generic
    ```

### `definitions.reference.entry`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: An entry in the collection that constitutes the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      entry: "Citation <n. 1>"
      type: generic
    ```
    ```yaml
    references:
      - entry: "Citation <n. 1>"
        type: generic
    ```

### `definitions.reference.filename`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The name of the electronic file containing the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      filename: CITATION.cff
      type: generic
    ```
    ```yaml
    references:
      - filename: CITATION.cff
        type: generic
    ```

### `definitions.reference.format`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The format in which a work is represented.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      format: "Citation File Format"
      type: generic
    ```
    ```yaml
    references:
      - format: "Citation File Format"
        type: generic
    ```

### `definitions.reference.identifiers`

- **type**: Array of [`definitions.identifier`](#definitionsidentifier) objects.
- **required**: `false`
- **description**: The identifier(s) of the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      identifiers:
        - description: "The concept DOI of the work."
          type: doi
          value: 10.5281/zenodo.1003149
        - description: "The versioned DOI for version 1.1.0 of the work."
          type: doi
          value: 10.5281/zenodo.4813122
        - description: "The Software Heritage identifier for version 1.1.0 of the work."
          type: swh
          value: "swh:1:dir:bc286860f423ea7ced246ba7458eef4b4541cf2d"
        - description: "The GitHub release URL of tag 1.1.0."
          type: url
          value: "https://github.com/citation-file-format/citation-file-format/releases/tag/1.1.0"
      type: generic
    ```
    ```yaml
    references:
      - identifiers:
          - description: "The concept DOI of the work."
            type: doi
            value: 10.5281/zenodo.1003149
          - description: "The versioned DOI for version 1.1.0 of the work."
            type: doi
            value: 10.5281/zenodo.4813122
          - description: "The Software Heritage identifier for version 1.1.0 of the work."
            type: swh
            value: "swh:1:dir:bc286860f423ea7ced246ba7458eef4b4541cf2d"
          - description: "The GitHub release URL of tag 1.1.0."
            type: url
            value: "https://github.com/citation-file-format/citation-file-format/releases/tag/1.1.0"
        type: generic
    ```

### `definitions.reference.institution`

- **type**: [`definitions.entity`](#definitionsentity)
- **required**: `false`
- **description**: The institution where a work has been produced or published.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      institution:
        name: "University of Arcadia"
      type: generic
    ```
    ```yaml
    references:
      - institution:
          name: "University of Arcadia"
        type: generic
    ```

### `definitions.reference.isbn`

- **type**: `string` with pattern [`^[0-9\- ]{10,17}X?$`](https://regex101.com/library/6oS1PA)
- **required**: `false`
- **description**: The [ISBN](https://en.wikipedia.org/wiki/International_Standard_Book_Number) of the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      isbn: "9781603095075"
      type: generic
    ```
    ```yaml
    references:
      - isbn: "9781603095075"
        type: generic
    ```

### `definitions.reference.issn`

- **type**: `string` with pattern [`^\d{4}-\d{3}[\dxX]$`](https://regex101.com/library/jqobq9)
- **required**: `false`
- **description**: The [ISSN](https://en.wikipedia.org/wiki/International_Standard_Serial_Number) of the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      issn: "2475-9066"
      type: generic
    ```
    ```yaml
    references:
      - issn: "2475-9066"
        type: generic
    ```

### `definitions.reference.issue`

- **type**: [Nonempty `string`](#yaml-strings) or `number`
- **required**: `false`
- **description**: The issue of a periodical in which a work appeared.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      issue: "42"
      type: generic
    ```
    ```yaml
    references:
      - issue: "42"
        type: generic
    ```
    ```yaml
    preferred-citation:
      issue: 42
      type: generic
    ```
    ```yaml
    references:
      - issue: 42
        type: generic
    ```

### `definitions.reference.issue-date`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The publication date of the issue of a periodical in which a work appeared.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      issue-date: "November/December 2021"
      type: generic
    ```
    ```yaml
    references:
      - issue-date: "November/December 2021"
        type: generic
    ```
    ```yaml
    preferred-citation:
      issue-date: "2021-12-31"
      type: generic
    ```
    ```yaml
    references:
      - issue-date: "2021-12-31"
        type: generic
    ```

### `definitions.reference.issue-title`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The name of the issue of a periodical in which the work appeared.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      issue-title: "Special Issue: Software Citation"
      type: generic
    ```
    ```yaml
    references:
      - issue-title: "Special Issue: Software Citation"
        type: generic
    ```

### `definitions.reference.journal`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The name of the journal/magazine/newspaper/periodical where the work was published.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      journal: "Journal of Open Source Software"
      type: generic
    ```
    ```yaml
    references:
      - journal: "Journal of Open Source Software"
        type: generic
    ```

### `definitions.reference.keywords`

- **type**: Array of [nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: Keywords pertaining to the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      keywords:
        - "software citation"
        - "citation file format"
        - research
      type: generic
    ```
    ```yaml
    references:
      - keywords:
          - "software citation"
          - "citation file format"
          - research
        type: generic
    ```

### `definitions.reference.languages`

- **type**: Array of [ISO 639](https://en.wikipedia.org/wiki/ISO_639) `string` with 2 or 3 characters and pattern [`^[a-z]{2,3}$`](https://regex101.com/library/aMqWLH)
- **required**: `false`
- **description**: The language identifier(s) of the work according to ISO 639 language strings.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      languages:
        - en
        - haw
      type: generic
    ```
    ```yaml
    references:
      - languages:
          - en
          - haw
        type: generic
    ```

### `definitions.reference.license`

- **type**: (Array of) [`definitions.license-enum`](#definitionslicense-enum).
- **required**: `false`
- **description**: The [SPDX license identifier(s)](https://spdx.dev/ids/) for the license(s) under which the work is made available.
When there are multiple licenses, it is assumed their relationship is OR, not AND.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      license: Apache-2.0
      type: generic
    ```
    ```yaml
    references:
      - license: Apache-2.0
        type: generic
    ```
    ```yaml
    preferred-citation:
      license:
        - Apache-2.0
        - MIT
      type: generic
    ```
    ```yaml
    references:
      - license:
          - Apache-2.0
          - MIT
        type: generic
    ```

### `definitions.reference.license-url`

- **type**: [`definitions.url`](#definitionsurl)
- **required**: `false`
- **description**: The URL of the license text under which the work is licensed (only for non-standard licenses not included in the [SPDX License List](#definitionslicense-enum)).
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      license-url: "https://obscure-licenses.com?id=1234"
      type: generic
    ```
    ```yaml
    references:
      - license-url: "https://obscure-licenses.com?id=1234"
        type: generic
    ```

### `definitions.reference.loc-end`

- **type**: [Nonempty `string`](#yaml-strings) or `integer`
- **required**: `false`
- **description**: The line of code in the file where the work ends.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      loc-end: "42"
      type: generic
    ```
    ```yaml
    preferred-citation:
      loc-end: 42
      type: generic
    ```
    ```yaml
    references:
      - loc-end: "42"
        type: generic
    ```
    ```yaml
    references:
      - loc-end: 42
        type: generic
    ```

### `definitions.reference.loc-start`

- **type**: [Nonempty `string`](#yaml-strings) or `integer`
- **required**: `false`
- **description**: The line of code in the file where the work starts.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      loc-start: "21"
      type: generic
    ```
    ```yaml
    preferred-citation:
      loc-start: 21
      type: generic
    ```
    ```yaml
    references:
      - loc-start: "21"
        type: generic
    ```
    ```yaml
    references:
      - loc-start: 21
        type: generic
    ```

### `definitions.reference.location`

- **type**: [`definitions.entity`](#definitionsentity)
- **required**: `false`
- **description**: The location of the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      location:
        name: "Library of the Unseen University"
      type: generic
    ```
    ```yaml
    references:
      - location:
          name: "Library of the Unseen University"
        type: generic
    ```

### `definitions.reference.medium`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The medium of the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      medium: "5¼-inch floppy disk"
      type: generic
    ```
    ```yaml
    references:
      - medium: "5¼-inch floppy disk"
        type: generic
    ```

### `definitions.reference.month`

- **type**: `integer` in range `1-12` or `enum` with values:
    - `1`
    - `2`
    - `3`
    - `4`
    - `5`
    - `6`
    - `7`
    - `8`
    - `9`
    - `10`
    - `11`
    - `12`
- **required**: `false`
- **description**: The month in which a work has been published.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      month: 7
      type: generic
    ```
    ```yaml
    preferred-citation:
      month: "7"
      type: generic
    ```
    ```yaml
    references:
      - month: 7
        type: generic
    ```
    ```yaml
    references:
      - month: "7"
        type: generic
    ```

### `definitions.reference.nihmsid`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The [NIHMSID](https://web.archive.org/web/20210802210057/https://www.ncbi.nlm.nih.gov/pmc/about/public-access-info/) of a work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      nihmsid: NIHMS236863
      type: generic
    ```
    ```yaml
    references:
      - nihmsid: NIHMS236863
        type: generic
    ```

### `definitions.reference.notes`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: Notes pertaining to the work. Note that this key should contain notes that may be picked up by some downstream tooling (e.g., reference managers), but not others (e.g., a software index).
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      notes: "Excellent reference! TODO Read for thesis."
      type: generic
    ```
    ```yaml
    references:
      - notes: "Excellent reference! TODO Read for thesis."
        type: generic
    ```

### `definitions.reference.number`

- **type**: `string` or `number`
- **required**: `false`
- **description**: The (library) [accession number](https://en.wikipedia.org/wiki/Accession_number) for a work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      number: "005.1-ABC"
      type: generic
    ```
    ```yaml
    references:
      - number: "005.1-ABC"
        type: generic
    ```
    ```yaml
    preferred-citation:
      number: 1234567
      type: generic
    ```
    ```yaml
    references:
      - number: 1234567
        type: generic
    ```
    ```yaml
    preferred-citation:
      number: 1.4
      type: generic
    ```
    ```yaml
    references:
      - number: 1.4
        type: generic
    ```

### `definitions.reference.number-volumes`

- **type**: [Nonempty `string`](#yaml-strings) or `integer`
- **required**: `false`
- **description**: The number of volumes making up the collection in which the work has been published.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      number-volumes: "4"
      type: generic
    ```
    ```yaml
    preferred-citation:
      number-volumes: 4
      type: generic
    ```
    ```yaml
    references:
      - number-volumes: "4"
        type: generic
    ```
    ```yaml
    references:
      - number-volumes: 4
        type: generic
    ```

### `definitions.reference.pages`

- **type**: [Nonempty `string`](#yaml-strings) or `integer`
- **required**: `false`
- **description**: The number of pages of the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      pages: "123"
      type: generic
    ```
    ```yaml
    preferred-citation:
      pages: 123
      type: generic
    ```
    ```yaml
    references:
      - pages: "123"
        type: generic
    ```
    ```yaml
    references:
      - pages: 123
        type: generic
    ```

### `definitions.reference.patent-states`

- **type**: Array of nonemtpy `string`
- **required**: `false`
- **description**: The states for which a patent is granted.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      patent-states:
        - Canada
        - IL
      type: generic
    ```
    ```yaml
    references:
      - patent-states:
          - Canada
          - IL
        type: generic
    ```

### `definitions.reference.pmcid`

- **type**: `string` with pattern [`^PMC[0-9]{7}$`](https://regex101.com/library/EsU1QH)
- **required**: `false`
- **description**: The [PMCID](https://web.archive.org/web/20210802210057/https://www.ncbi.nlm.nih.gov/pmc/about/public-access-info/) of a work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      pmcid: PMC3134971
      type: generic
    ```
    ```yaml
    references:
      - pmcid: PMC3134971
        type: generic
    ```

### `definitions.reference.publisher`

- **type**: [`definitions.entity`](#definitionsentity)
- **required**: `false`
- **description**: The publisher who has published the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      publisher:
        name: "Open Access Publishing House"
      type: generic
    ```
    ```yaml
    references:
      - publisher:
          name: "Open Access Publishing House"
        type: generic
    ```

### `definitions.reference.recipients`

- **type**: Array of [`definitions.person`](#definitionsperson) and/or [`definitions.entity`](#definitionsentity) objects.
- **required**: `false`
- **description**: The recipient(s) of a personal communication.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      recipients:
        - name: "Recipient entity of personal communication"
        - family-names: Recipient
          given-names: Communication
      type: generic
    ```
    ```yaml
    references:
      - recipients:
          - name: "Recipient entity of personal communication"
          - family-names: Recipient
            given-names: Communication
        type: generic
    ```

### `definitions.reference.repository`

- **type**: [`definitions.url`](#definitionsurl)
- **required**: `false`
- **description**: The URL of the work in a repository/archive (when the repository is neither a source code repository nor a build artifact repository).
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      repository: "https://ascl.net/2105.013"
      type: generic
    ```
    ```yaml
    references:
      - repository: "https://ascl.net/2105.013"
        type: generic
    ```

### `definitions.reference.repository-artifact`

- **type**: [`definitions.url`](#definitionsurl)
- **required**: `false`
- **description**: The URL of the work in a build artifact/binary repository.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      repository-artifact: "https://search.maven.org/artifact/org.corpus-tools/cff-maven-plugin/0.4.0/maven-plugin"
      type: generic
    ```
    ```yaml
    references:
      - repository-artifact: "https://search.maven.org/artifact/org.corpus-tools/cff-maven-plugin/0.4.0/maven-plugin"
        type: generic
    ```

### `definitions.reference.repository-code`

- **type**: [`definitions.url`](#definitionsurl)
- **required**: `false`
- **description**: The URL of the work in a source code repository.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      repository-code: "https://github.com/citation-file-format/my-research-software"
      type: generic
    ```
    ```yaml
    references:
      - repository-code: "https://github.com/citation-file-format/my-research-software"
        type: generic
    ```

### `definitions.reference.scope`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The scope of the reference, e.g., the section of the work it adheres to.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      scope: "Supplement 2: Additional material"
      type: generic
    ```
    ```yaml
    references:
      - scope: "Supplement 2: Additional material"
        type: generic
    ```

### `definitions.reference.section`

- **type**: `string` or `number`
- **required**: `false`
- **description**: The section of a work that is referenced.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      section: "Section 7: Software citation"
      type: generic
    ```
    ```yaml
    references:
      - section: "Section 7: Software citation"
        type: generic
    ```
    ```yaml
    preferred-citation:
      section: 7
      type: generic
    ```
    ```yaml
    references:
      - section: 7
        type: generic
    ```
    ```yaml
    preferred-citation:
      section: 7.1
      type: generic
    ```
    ```yaml
    references:
      - section: 7.1
        type: generic
    ```

### `definitions.reference.senders`

- **type**: Array of [`definitions.person`](#definitionsperson) and/or [`definitions.entity`](#definitionsentity) objects.
- **required**: `false`
- **description**: The sender(s) of a personal communication.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      senders:
        - name: "Sender entity of personal communication"
        - family-names: Sender
          given-names: Communication
      type: generic
    ```
    ```yaml
    references:
      - senders:
          - name: "Sender entity of personal communication"
          - family-names: Sender
            given-names: Communication
        type: generic
    ```

### `definitions.reference.start`

- **type**: [Nonempty `string`](#yaml-strings) or `integer`
- **required**: `false`
- **description**: The start page of the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      start: "42"
      type: generic
    ```
    ```yaml
    preferred-citation:
      start: 42
      type: generic
    ```
    ```yaml
    references:
      - start: "42"
        type: generic
    ```
    ```yaml
    references:
      - start: 42
        type: generic
    ```

### `definitions.reference.status`

- **type**: `enum` with values:
    - `abstract`
    - `advance-online`
    - `in-preparation`
    - `in-press`
    - `preprint`
    - `submitted`
- **required**: `false`
- **description**: The publication status of the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      status: submitted
      type: generic
    ```
    ```yaml
    references:
      - status: submitted
        type: generic
    ```

### `definitions.reference.term`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The term being referenced if the work is a dictionary or encyclopedia.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      term: Citation
      type: generic
    ```
    ```yaml
    references:
      - term: Citation
        type: generic
    ```

### `definitions.reference.thesis-type`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The type of the thesis that is the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      thesis-type: "PhD thesis"
      type: generic
    ```
    ```yaml
    references:
      - thesis-type: "PhD thesis"
        type: generic
    ```

### `definitions.reference.title`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `true`
- **description**: The title of the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      title: "Towards better software citation"
      type: generic
    ```
    ```yaml
    references:
      - title: "Towards better software citation"
        type: generic
    ```

### `definitions.reference.translators`

- **type**: Array of [`definitions.person`](#definitionsperson) and/or [`definitions.entity`](#definitionsentity) objects.
- **required**: `false`
- **description**: The translator(s) of a work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      translators:
        - name: "Research Translators Ltd."
        - family-names: Lator
          given-names: Trans
      type: generic
    ```
    ```yaml
    references:
      - translators:
          - name: "Research Translators Ltd."
          - family-names: Lator
            given-names: Trans
        type: generic
    ```

### `definitions.reference.type`

- **type**: `enum` with values:
    - `art`
    - `article`
    - `audiovisual`
    - `bill`
    - `blog`
    - `book`
    - `catalogue`
    - `conference-paper`
    - `conference`
    - `data`
    - `database`
    - `dictionary`
    - `edited-work`
    - `encyclopedia`
    - `film-broadcast`
    - `generic`
    - `government-document`
    - `grant`
    - `hearing`
    - `historical-work`
    - `legal-case`
    - `legal-rule`
    - `magazine-article`
    - `manual`
    - `map`
    - `multimedia`
    - `music`
    - `newspaper-article`
    - `pamphlet`
    - `patent`
    - `personal-communication`
    - `proceedings`
    - `report`
    - `serial`
    - `slides`
    - `software-code`
    - `software-container`
    - `software-executable`
    - `software-virtual-machine`
    - `software`
    - `sound-recording`
    - `standard`
    - `statute`
    - `thesis`
    - `unpublished`
    - `video`
    - `website`
- **required**: `true`
- **description**: The type of the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
        type: generic
    ```
    ```yaml
    references:
      - type: generic
    ```

### `definitions.reference.url`

- **type**: [`definitions.url`](#definitionsurl)
- **required**: `false`
- **description**: The URL of the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      type: generic
      url: "https://citation-file-format.github.io/"
    ```
    ```yaml
    references:
      - type: generic
        url: "https://citation-file-format.github.io/"
    ```

### `definitions.reference.version`

- **type**: [`definitions.version`](#definitionsversion)
- **required**: `false`
- **description**: The version of the work.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      type: generic
      version: 0.3.12
    ```
    ```yaml
    references:
      - type: generic
        version: 0.3.12
    ```

### `definitions.reference.volume`

- **type**: [Nonempty `string`](#yaml-strings) or `integer`
- **required**: `false`
- **description**: The volume of the periodical in which a work appeared.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      type: generic
      volume: "5"
    ```
    ```yaml
    preferred-citation:
      type: generic
      volume: 5
    ```
    ```yaml
    references:
      - type: generic
        volume: "5"
    ```
    ```yaml
    references:
      - type: generic
        volume: 5
    ```

### `definitions.reference.volume-title`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: `false`
- **description**: The title of the volume in which the work appeared.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      type: generic
      volume-title: "Volume II: How it went on"
    ```
    ```yaml
    references:
      - type: generic
        volume-title: "Volume II: How it went on"
    ```

### `definitions.reference.year`

- **type**: [Nonempty `string`](#yaml-strings) or `integer`
- **required**: `false`
- **description**: The year in which a work has been published.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      type: generic
      year: "2021"
    ```
    ```yaml
    preferred-citation:
      type: generic
      year: 2021
    ```
    ```yaml
    references:
      - type: generic
        year: "2021"
    ```
    ```yaml
    references:
      - type: generic
        year: 2021
    ```

### `definitions.reference.year-original`

- **type**: [Nonempty `string`](#yaml-strings) or `integer`
- **required**: `false`
- **description**: The year of the original publication.
- **usage**:<br><br>
    ```yaml
    preferred-citation:
      type: generic
      year-original: "1978"
    ```
    ```yaml
    preferred-citation:
      type: generic
      year-original: 1978
    ```
    ```yaml
    references:
      - type: generic
        year-original: "1978"
    ```
    ```yaml
    references:
      - type: generic
        year-original: 1978
    ```

### `definitions.region`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: N/A
- **description**: A region.
- **usage**:<br><br>
    ```yaml
    authors:
      - name: "The Research Software Project"
        region: Renfrewshire
    ```

### `definitions.swh-identifier`

- **type**: `string` with pattern [`^swh:1:(snp|rel|rev|dir|cnt):[0-9a-fA-F]{40}$`](https://regex101.com/library/o399MX)
- **required**: N/A
- **description**: The [Software Heritage](https://www.softwareheritage.org/) identifier (without further qualifiers such as origin, visit, anchor, path). Note: Software Heritage identifiers are documented here: https://docs.softwareheritage.org/devel/swh-model/persistent-identifiers.html.
- **usage**:<br><br>
    ```yaml
    identifiers:
      - type: swh
        value: "swh:1:rev:309cf2674ee7a0749978cf8265ab91a60aea0f7d"
    ```

### `definitions.tel`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: N/A
- **description**: A telephone number.
- **usage**:<br><br>
    ```yaml
    authors:
      - family-names: Druskat
        given-names: Stephan
        tel: +12-345-6789098
    ```
    ```yaml
    authors:
      - name: "The Research Software Project"
        tel: +12-345-6789098
    ```

### `definitions.url`

- **type**: [Nonempty `string`](#yaml-strings)
- **required**: N/A
- **description**: A URL. Supported URLs start with one of:
    - `https://`
    - `http://`
    - `ftp://`
    - `sftp://`
- **usage**
    ```yaml
    url: "https://citation-file-format.github.io/"
    ```
    ```yaml
    authors:
      - name: "The Research Software Project"
        website: "https://research-software-project.org"
    ```
    ```yaml
    references:
      - type: generic
        url: "sftp://files.research-software-project.org"
    ```

### `definitions.version`

- **type**: [Nonempty `string`](#yaml-strings) or `number`
- **required**: N/A
- **description**: The version of a work.
- **usage**:<br><br>
    ```yaml
    version: "1.2.0"
    ```
    ```yaml
    version: 1.2
    ```
    ```yaml
    version: "21.10 (Impish Indri)"
    ```
# Citation File Format

[![Build Status](https://github.com/citation-file-format/citation-file-format/workflows/testing/badge.svg?branch=main)](https://github.com/citation-file-format/citation-file-format/actions/workflows/testing.yml?query=branch%3Amain)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1003149.svg)](https://doi.org/10.5281/zenodo.1003149)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Project homepage](https://img.shields.io/badge/Project%20homepage-citation--file--format.github.io-ff0080)](https://citation-file-format.github.io)

The Citation File Format lets you provide citation metadata for software or datasets 
in plaintext files that are easy to read by both humans and machines.

## Structure

You can specify citation metadata for your software (or dataset) in a file named `CITATION.cff`. 
This is what a typical `CITATION.cff` file may look like for research software:

```yaml
cff-version: 1.2.0
message: If you use this software, please cite it using these metadata.
title: My Research Software
abstract: This is my awesome research software. It does many things.
authors:
  - family-names: Druskat
    given-names: Stephan
    orcid: "https://orcid.org/0000-0003-4925-7248"
  - name: "The Research Software project"
version: 0.11.2
date-released: "2021-07-18"
identifiers:
  - description: This is the collection of archived snapshots of all versions of My Research Software
    type: doi
    value: "10.5281/zenodo.123456"
  - description: This is the archived snapshot of version 0.11.2 of My Research Software
    type: doi
    value: "10.5281/zenodo.123457"
license: Apache-2.0
repository-code: "https://github.com/citation-file-format/my-research-software"
```

In addition, the Citation File Format allows you to

- provide references to works that your software or dataset builds on ([see here for more info](schema-guide.md#referencing-other-work));
- ask people to cite a different, related work instead of the software or dataset itself ([see here for more info](schema-guide.md#credit-redirection)).

## Format specifications :books:

**You can find the complete format specifications in the [Guide to the Citation File Format schema](schema-guide.md).**

## Why should I add a `CITATION.cff` file to my repository? :bulb:

When you do this, great things may happen:

1. Users of your software can easily cite it using the metadata from `CITATION.cff`!
2. If your repository is hosted on GitHub, they will [show the citation information in the sidebar](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/creating-a-repository-on-github/about-citation-files), which makes it easy for visitors to cite your software or dataset correctly.
3. When you publish your software on Zenodo via the [GitHub-Zenodo integration](https://guides.github.com/activities/citable-code/), they will use the metadata from your `CITATION.cff` file.
4. People can import the correct reference to your software into the [Zotero](https://www.zotero.org) reference manager via a [browser plugin](https://www.zotero.org/download/).

## Creation :heavy_plus_sign:

To create a `CITATION.cff` file, you can 

- use the [**cffinit** website](https://citation-file-format.github.io/cff-initializer-javascript/#/),
- copy and paste the [example snippet](#structure), and adapt it to your needs, or
- create a new file called `CITATION.cff` using the *Add file* button on GitHub, and use the template they provide.

## Validation :heavy_check_mark:

You can validate your `CITATION.cff` file on the command line with the [`cffconvert` Python package](https://pypi.org/project/cffconvert/):

```shell
# Install cffconvert with pip in user space
python3 -m pip install --user cffconvert

# Validate your CFF file
cffconvert --validate
```

If you get a Traceback with error messages, look for the relevant validation error and fix it.
If the output is very long, it may help if you search it for lines starting with `jsonschema.exceptions.ValidationError`.

If you prefer to use Docker, you can use the [`cffconvert` Docker image](https://hub.docker.com/r/citationcff/cffconvert):

```bash
cd <directory-containing-your-CITATION.cff>
docker run --rm -v ${PWD}:/app citationcff/cffconvert --validate
```

<!-- Later, this should link to tutorials -->

## Tools to work with `CITATION.cff` files :wrench:

There is tooling available to work with `CITATION.cff` files to do different things:
create new files, edit existing files, validate existing files, convert files from the Citation File Format into another format.
The following table gives an overview of the tools that we know about. If there is a tool missing from this table, please [open a new issue](https://github.com/citation-file-format/citation-file-format/issues/new/choose) and let us know.

|                | Creation                                                                        | Editing/Updating                                                    | Validation                                                                      | Conversion                                                                                                                                                       |
| -------------- | ------------------------------------------------------------------------------- | ------------------------------------------------------------------- | ------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Command line   |                                                                                 |                                                                     | • [cffconvert](#validation-heavy_check_mark)                                    | • [cffconvert](https://pypi.org/project/cffconvert/)<br> • [bibtex-to-cff](https://github.com/monperrus/bibtexbrowser/)                                          |
| GitHub Actions |                                                                                 |                                                                     | [cff-validator](https://github.com/marketplace/actions/cff-validator)           | • [cffconvert](https://github.com/marketplace/actions/cffconvert)<br>• [codemeta2cff](https://github.com/caltechlibrary/codemeta2cff)                            |
| GitHub Bot     |                                                                                 |                                                                     | [#238](https://github.com/citation-file-format/citation-file-format/issues/238) |                                                                                                                                                                  |
| Docker         |                                                                                 |                                                                     | [cffconvert Docker image](#validation-heavy_check_mark)                         | [cffconvert Docker image](https://hub.docker.com/r/citationcff/cffconvert)                                                                                       |
| Go             |                                                                                 |                                                                     |                                                                                 | • [datatools/codemeta2cff](https://github.com/caltechlibrary/datatools/)                                                                                         |
| Java           | • [CFF Maven plugin](https://github.com/hexatomic/cff-maven-plugin)             | • [CFF Maven plugin](https://github.com/hexatomic/cff-maven-plugin) |                                                                                 | • [CFF Maven plugin](https://github.com/hexatomic/cff-maven-plugin)                                                                                              |
| JavaScript     |                                                                                 |                                                                     |                                                                                 | • [Citation.js](https://citation.js.org/) [plugin](https://www.npmjs.com/package/@citation-js/plugin-software-formats)                                           |
| Julia          |                                                                                 |                                                                     | • [Bibliography.jl](https://github.com/Humans-of-Julia/Bibliography.jl)                    | • [Bibliography.jl](https://github.com/Humans-of-Julia/Bibliography.jl)                                                                             |
| PHP            |                                                                                 |                                                                     |                                                                                 | • [bibtex-to-cff](https://github.com/monperrus/bibtexbrowser/)                                                                                                   |
| Python         |                                                                                 | • [doi2cff](https://github.com/citation-file-format/doi2cff)        | • [cffconvert](#validation-heavy_check_mark)                                    | • [cffconvert](https://github.com/citation-file-format/cff-converter-python)<br>• [doi2cff](https://github.com/citation-file-format/doi2cff)                     |
| R              |                                                                                 |                                                                     |                                                                                 | • [citation](https://cran.r-project.org/web/packages/citation/)<br>• [r2cff](https://github.com/ocbe-uio/RCFF)<br>• [handlr](https://github.com/ropensci/handlr)<br>• [cffr](https://CRAN.R-project.org/package=cffr) |
| Ruby           | • [ruby-cff](https://github.com/citation-file-format/ruby-cff)                  | • [ruby-cff](https://github.com/citation-file-format/ruby-cff)      | • [ruby-cff](https://github.com/citation-file-format/ruby-cff)                  | • [ruby-cff](https://github.com/citation-file-format/ruby-cff)                                                                                                   |
| TypeScript     |                                                                                 |                                                                     |                                                                                 | [#28](https://github.com/citation-file-format/citation-file-format/issues/28#issuecomment-892105342)                                                             |
| Website        | • [cffinit](https://citation-file-format.github.io/cff-initializer-javascript/) |                                                                     |                                                                                 |                                                                                                                                                                  |

## Maintainers :nerd_face:

The Citation File Format schema is maintained by

- Stephan Druskat ([@sdruskat](https://github.com/sdruskat/))
- Jurriaan H. Spaaks ([@jspaaks](https://github.com/jspaaks/))

## Contributing :handshake:

The Citation File Format is a collaborative project and we welcome suggestions and contributions. We hope one of the invitations below works for you, but if not, please let us know!

:running: **I'm busy, I only have 1 minute**
- Tell a friend about the Citation File Format, or tweet about it!
- Give the project a star :star:!

:hourglass_flowing_sand: **I've got 10 minutes - tell me what I should do**
- Create a `CITATION.cff` file for your repository.
- Suggest ideas for how you would like to use the Citation File Format, or for an improvement to the format or its tooling.
- If you know how to validate `CITATION.cff` files, help someone with a validation problem and look at the [issues labeled ![GitHub labels](https://img.shields.io/github/labels/citation-file-format/citation-file-format/validation)](https://github.com/citation-file-format/citation-file-format/issues?q=is%3Aopen+is%3Aissue+label%3A%22help+wanted%22+label%3Avalidation)

:computer: **I've got a few hours to work on this**
- Help create tooling for the community by looking at the [issues labeled ![GitHub labels](https://img.shields.io/github/labels/citation-file-format/citation-file-format/tooling)](https://github.com/citation-file-format/citation-file-format/issues?q=is%3Aopen+is%3Aissue+label%3A%22help+wanted%22+label%3Atooling)

:tada: **I want to help grow the community**
- Write a blog post or news item for your own community.
- Organise a hack event or workshop to help others use or improve the Citation File Format.

Please read the more detailed [contributing guidelines](CONTRIBUTING.md) and [open a GitHub issue](https://github.com/citation-file-format/citation-file-format/issues) to suggest a new idea or let us know about bugs. Please put up pull requests for changes to the format and schema against the `develop` branch!

## License :balance_scale:

Copyright © 2016ff. The Citation File Format Contributors

This work is licensed under a [Creative Commons Attribution 4.0 International (CC-BY-4.0)](https://creativecommons.org/licenses/by/4.0/legalcode) license.

## Acknowledgments :pray:

**We'd like to thank everyone who has contributed to the Citation File Format!**  
They are listed in the [`CITATION.cff`](CITATION.cff) file for this repository. Please open an issue if you find that you are missing from the file.

We gratefully acknowledge support from:

- The [Institute for Software Technology](https://www.dlr.de/sc/en/desktopdefault.aspx/) of the [German Aerospace Center (DLR)](https://www.dlr.de/en/)
- The [Netherlands eScience Center](https://www.esciencecenter.nl/)
- The [Software Sustainability Institute](https://software.ac.uk/)

### Research notice
Please note that this repository is participating in a study into sustainability
of open source projects. Data will be gathered about this repository for
approximately the next 12 months, starting from June 2021.

Data collected will include number of contributors, number of PRs, time taken to
close/merge these PRs, and issues closed.

For more information, please visit
[our informational page](https://sustainable-open-science-and-software.github.io/) or download the [participant information sheet](https://sustainable-open-science-and-software.github.io/assets/PIS_sustainable_software.pdf).
# Code of Conduct

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of gender,
gender identity and expression, age, sexual identity and orientation,
disability, physical appearance, body size, race, ethnicity, religion (or lack
thereof), technology choices, level of experience, or other group status.

## Standards

Examples of behavior that contributes to creating a positive environment
include:

- Using welcoming and inclusive language
- Being respectful of differing viewpoints and experiences
- Gracefully accepting constructive criticism
- Focusing on what is best for the community
- Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

- The use of sexualized language or imagery and unwelcome sexual attention or
  advances
- Trolling, insulting/derogatory comments, and personal or political attacks
- Public or private harassment
- Publishing others' private information, such as a physical or electronic
  address, without explicit permission
- Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at *cff-coc /at\ sdruskat \dot/ net*. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's maintenance team.

## Attribution

This Code of Conduct is adapted from the 
[Contributor Covenant](https://www.contributor-covenant.org/), 
version 1.4, available at 
https://www.contributor-covenant.org/version/1/4/code-of-conduct.html, 
and the [WSSSPE5.1](http://wssspe.researchcomputing.org.uk/wssspe5-1/) 
Code of Conduct, available at 
http://wssspe.researchcomputing.org.uk/wssspe5-1/wssspe5-1-code-of-conduct/.# Introduction

**Thank you** for considering to contribute to the **Citation File Format (CFF) and its schema**!  
If you intend to contribute to another part of the Citation File Format project, 
for example a specific software package for working with CFF files,
please contribute to the respective repository ([list of repositories in the `citation-file-format` GitHub organization](https://github.com/orgs/citation-file-format/repositories)).

**Please follow these guidelines.** Their purpose is to make both contributing and accepting contributions easier for all parties involved.

There are many **ways to contribute**, e.g.:

- Tell a friend or colleague about the Citation File Format, or tweet about it
- Write blog posts, tutorials, etc. about the Citation File Format
- Review the format and its schema and documentation
- Improve wording in any prose output, including the specifications
- Create a new, better version of the schema and specifications
- Improve automated tests, continuous integration, documentation, etc.

# Ground Rules

Your contribution to CFF is valued, and it should be an enjoyable experience.
To ensure this, there is the CFF 
[Code of Conduct](https://github.com/citation-file-format/citation-file-format/blob/master/CODE_OF_CONDUCT.md), which you are required to follow.

Please always start any contribution that will change the contents of this repository by [creating a new issue](https://github.com/citation-file-format/citation-file-format/issues/new/choose). 
This way, 

- you can make sure that you don't invest your valuable time in something that may not be merged;
- we can make sure that your contribution is something that will improve CFF, 
is in scope, 
and aligns with the roadmap for the Citation File Format.

# Your First Contribution

If you are unsure where to begin with your contribution to CFF, have a look at the
[open issues in this repository](https://github.com/citation-file-format/citation-file-format/issues), 
and see if you can identify one that you would like to work on.

If you have never contributed to an open source project, you may find this tutorial helpful:
[How to Contribute to an Open Source Project on GitHub](https://app.egghead.io/playlists/how-to-contribute-to-an-open-source-project-on-github).

# Getting started

This is the workflow for contributions to this repository:

1. [Create a new issue](https://github.com/citation-file-format/citation-file-format/issues/new/choose) 
to discuss the changes you want to make with the maintainers and community
2. Fork the repository
3. Create a branch in your fork of the repository:
  - If you want to make changes to the format and its schema: checkout `develop` and create your branch from it.
  - If you want to make changes that don't affect the format and its schema, like fixing typos or formatting: checkout `main` and create your branch from it.
  - When in doubt which branch to fork from, ask in the issue you have created.
4. Make changes in the new branch in your fork
5. Take note of the [code of conduct](https://github.com/citation-file-format/citation-file-format/blob/main/CODE_OF_CONDUCT.md)
6. Create a pull request
7. Address any comments that have come up during review
8. If and when your pull request has been merged, you can delete your branch (or the whole forked repository)

This workflow is loosely based on GitHub flow, and you can find more information in the [GitHub flow documentation](https://docs.github.com/en/get-started/quickstart/github-flow).

## Working with examples and tests

There is a collection of test `CITATION.cff` files in the `examples` directory.
It may be a good idea to add a test if you modify the specification,
or if you discover a part of the specification that is not covered
by tests. If you modify anything in the `tests` or `examples` directory, it is
advised to run these tests locally on your computer prior to submitting
a pull request. However, if that's not possible, you still can submit
the pull request and later check the status of the tests for your
pull requests on GitHub. Please see [`examples/README.md`](examples/README.md) file for further
details about tests and instructions how to run them locally.

<!--
TODO Include a link to README.dev.md here once it exists! 
See https://github.com/citation-file-format/citation-file-format/issues/301
-->

# How to submit an issue

If you find a security vulnerability anywhere, do NOT open an issue. Email *cff-
security /at\ sdruskat \dot/ net* instead.

This repository uses issue templates, so that when
you create a new issue, you can simply fill it in and everything that is
needed to work on your issue will be there.

If no issue template exists for your specific case, please use the template for "Any other issue".

# How to suggest a feature or enhancement

See [How to submit an issue](#how-to-submit-an-issue).

# FAQ

- **These guidelines do not address aspect XYZ! What should I do now?**

Please [submit an issue](https://github.com/citation-file-format/citation-file-format/issues/new/choose), 
asking for clarification and addition to the guidelines.
**Related issues**

Refs: #ISSUE_NUMBER

(For autoclosure of issues when PR is merged use `Fixes #<issue-number>` syntax)

**Describe the changes made in this pull request**

**Review checklist**

- [ ] Please check if the pull request is against the correct branch  
(format/schema/semantic documentation changes: `develop`; typos, meta files, etc.: `main`)
- [ ] Please check if all changes are recorded in `CHANGELOG.md` and adapt if necessary.
- [ ] Please run tests locally.
<!-- 
CONTRIBUTORS: Please replace <do other things> in the snippet below 
with something that reviewers should do to test and review your contribution!
-->
```bash
cd $(mktemp -d --tmpdir cff.XXXXXX)
git clone https://github.com/citation-file-format/citation-file-format .
git checkout <branch>
python3 -m venv env
source env/bin/activate
pip install --upgrade pip wheel setuptools
pip install -r requirements.txt
pytest
<do other things>
```
<!-- 
CONTRIBUTORS: Please replace `<do other things>` in the checklist item below 
with something that reviewers should do additionally
to test and review your contribution!
-->
- [ ] `<do other things>`
---
name: Any other issue
about: Any other issue
title: ''
labels: ''
assignees: ''

---
---
name: Idea for a schema/format enhancement
about: I have an idea for how to improve the Citation File Format schema
title: ''
labels: 'enhancement'
assignees: ''

---
---
name: Bug report
about: I think I may have found a bug!
title: ''
labels: 'bug'
assignees: ''

---
---
name: Validation issue
about: My CITATION.cff file is not valid and I need help
title: ''
labels:
  - 'help wanted'
  - validation
assignees: ''

---

<!--
Please have a look at the section on validating your CITATION.cff file before submitting this issue: https://github.com/citation-file-format/citation-file-format#validation-heavy_check_mark.
-->

<!-- Don't delete the following note! -->
Please note that issues with the validity of single CITATION.cff files may take some time to be picked up by the Citation File Format maintainers themselves. Therefore, if you are reading this issue and know how to validate `CITATION.cff` files, please help out if you can!

## Invalid `CITATION.cff` file

```yaml
Paste your invalid CITATION.cff file here!
```

## Context

- Please describe what you want to do (create a valid `CITATION.cff` file, show the citation information in the sidebar of your GitHub repository, ...)
---
name: Tool idea
about: I have an idea for, or need, specific tooling for CITATION.cff files
title: ''
labels: 
  - 'help wanted'
  - tooling
assignees: ''

---

<!--
Please have a look at the section on tooling for CITATION.cff files before submitting this issue: https://github.com/citation-file-format/citation-file-format#tools-to-work-with-citationcff-files-wrench.
-->

### Who needs the new tool?

<!-- 
Describe the type of user that would need the new tool. 
Examples:
- Researchers
- A research software registry
-->

### What should the new tool enable users to do?

<!-- 
Describe a goal that users could achieve with the new tool.
You can start with "I want to ..." or similar. 
Example:
I want to convert the metadata in CITATION.cff files to BibTeX.
-->

### What benefit would the new tool add?

<!-- 
Describe the benefit the tool will give users. 
Example:
Users will be able to use the metadata from CITATION.cff files to cite software in papers they write in LaTeX.
-->

### Implementation suggestions

<!-- 
Do you have a specific technology or programming language in mind for the implementation of the new tool? 
If so, please describe why.
Examples:
- This could be a website hosted by the software registries that want to use this functionality.
- This could be a Python package with a command-line interface.
-->

### Can you help?

<!--
Indicate if you will be able to work on the implementation yourself.
-->
- [ ] Yes
- [ ] No
# Testing 

## Overview

We maintain a collection of `.cff` files for testing. In these files we try to cover
various possible scenarious described in the specification. During the test, these
files are validated using `cffconvert` (a command line utility to read a `.cff` file
and convert it to BibTex and other formats). We check that all files are validated or
rejected as expected.

## Dependencies

- Python 3.6 or higher
- [cffconvert](https://pypi.org/project/cffconvert/) version 1.0.0 or higher
- [pytest](https://pypi.org/project/pytest/)

You can install `cffconvert` and `pytest` using the `requirements.txt` file from
this directory by entering the following command:

    python3 -m pip install -r requirements.txt

in the root directory of the clone of this repository.

## Files and directories

The tests are organised into subdirectories named after minimal versions of the CFF
specification supporting the tested features. 

Each test is placed in a directory with some informative name. We suggest that a name
of the directory with a test for failing validation should start with `fail-`, otherwise
its name may be arbitrary, but we suggest to follow the same style as for existing tests.

Inside the directory for a specific test, there should be two files: a `CITATION.cff`
file to be validated during the test, and a Python script to run the test. The name of
the latter must start with `test_`, followed by some string with underscores, based on
the name of the directory for this specific test (see existing tests for some examples).
To create such a Python script, you can copy and rename an existing test script from some
other test directory, only paying attention whether it checks that validation should
pass or fail - there are only two kinds of test scripts at the moment.

## Running tests

To run tests, enter

    pytest

You should expect to see output similar to this:

```
tests/validate.py::test[./tests/1.1.0/reference-article/CITATION.cff] PASSED                         [  2%]
tests/validate.py::test[./tests/1.1.0/reference-art/CITATION.cff] PASSED                             [  5%]
tests/validate.py::test[./tests/1.1.0/reference-book/CITATION.cff] PASSED                            [  7%]
tests/validate.py::test[./tests/1.1.0/fail-bad-identifier-type-in-root/CITATION.cff] PASSED          [ 10%]

...

tests/validate.py::test[./tests/1.0.3/reference-blog/CITATION.cff] PASSED                            [ 94%]
tests/validate.py::test[./tests/1.0.3/software-container/CITATION.cff] PASSED                        [ 97%]
tests/validate.py::test[./tests/1.0.3/software-executable/CITATION.cff] PASSED                       [100%]

```

surrounded with some additional information messages.
# Software with a DOI

Note that [Smith et al., 2016](https://doi.org/10.7717/peerj-cs.86), p. 12 recommend

> [...] the use of DOIs as the unique identifier due to their common usage and
acceptance, particularly as they are the standard for other digital products
such as publications.

Furthermore, DOIs should point to a "unique, specific software version"
([Smith et al., 2016](https://doi.org/10.7717/peerj-cs.86), p. 12). Also it is recommended ([Smith et al., 2016](https://doi.org/10.7717/peerj-cs.86), p. 13) that:

> the [DOI] should resolve to a persistent landing page that contains metadata
and a link to the software itself, rather than directly to the source code files,
repository, or executable.

Therefore, a minimal `CITATION.cff` file in such a case would look similar to [CITATION.cff](CITATION.cff).
# Software with a doi (expanded)

Same as [../software-with-a-doi/README.md](../software-with-a-doi/README.md), but with additional properties.
# software with reference of type edited-work

Note that the editors of the edited work must be specified under the `authors`
key. Specific citation styles may or may not attach a suffix to the authors,
such as ", eds." or similar.
For software without a DOI, it is recommended that "the metadata should still
provide information on how to access the specific software, but this may be a
company’s product number or a link to a website that allows the software be
purchased." [Smith et al., 2016](https://doi.org/10.7717/peerj-cs.86), p. 13. Furthermore, "if the version number and
release date are not available, the download date can be used. Similarly, the
contact name/email is an alternative to the location/repository."
([Smith et al., 2016](https://doi.org/10.7717/peerj-cs.86), p. 7).

Hence, for closed-source software without a DOI for which the version number
and release date cannot be determined, a `CITATION.cff` file could look like
[./CITATION.cff](./CITATION.cff).
# Software with a further reference

Where authors wish to encourage citation of an outline paper with citation
of their software, we recommend the use of [reference keys](#references-optional)
to highlight the existence of further references.
We recognize that there are certain situations where it may not be possible to
follow the recommended best-practice. For example, if (1) the software authors
did not register a DOI and/or release a specific version, or (2) the version of
the software used does not match what is available to cite. In those cases,
falling back on a combination of the repository URL and version number/commit
hash would be an appropriate way to cite the software used ([Smith et al., 2016](https://doi.org/10.7717/peerj-cs.86), p. 12).
Borrowed from [`ruby-cff`](https://github.com/citation-file-format/ruby-cff)'s fixtures https://github.com/citation-file-format/ruby-cff/tree/ec3aeff98855690d56ac2cd35f7181e5d9ee26f3/test/files.Borrowed from [`ruby-cff`](https://github.com/citation-file-format/ruby-cff)'s fixtures https://github.com/citation-file-format/ruby-cff/tree/ec3aeff98855690d56ac2cd35f7181e5d9ee26f3/test/files.

- Updated `cff-version` to `"1.2.0"`

Fails because
- `date-released` is not in YYYY-MM-DD format
- `pages` should be an integer describing the number of pages
Borrowed from [`ruby-cff`](https://github.com/citation-file-format/ruby-cff)'s fixtures https://github.com/citation-file-format/ruby-cff/tree/ec3aeff98855690d56ac2cd35f7181e5d9ee26f3/test/files.Borrowed from [`ruby-cff`](https://github.com/citation-file-format/ruby-cff)'s fixtures https://github.com/citation-file-format/ruby-cff/tree/ec3aeff98855690d56ac2cd35f7181e5d9ee26f3/test/files.# Software with a DOI

Note that [Smith et al., 2016](https://doi.org/10.7717/peerj-cs.86), p. 12 recommend

> [...] the use of DOIs as the unique identifier due to their common usage and
acceptance, particularly as they are the standard for other digital products
such as publications.

Furthermore, DOIs should point to a "unique, specific software version"
([Smith et al., 2016](https://doi.org/10.7717/peerj-cs.86), p. 12). Also it is recommended ([Smith et al., 2016](https://doi.org/10.7717/peerj-cs.86), p. 13) that:

> the [DOI] should resolve to a persistent landing page that contains metadata
and a link to the software itself, rather than directly to the source code files,
repository, or executable.

Therefore, a minimal `CITATION.cff` file in such a case would look similar to [CITATION.cff](CITATION.cff).
Borrowed from [`ruby-cff`](https://github.com/citation-file-format/ruby-cff)'s fixtures https://github.com/citation-file-format/ruby-cff/tree/ec3aeff98855690d56ac2cd35f7181e5d9ee26f3/test/files.

- Updated `cff-version` to `"1.2.0"`
- Updated `date-released` to comply with `YYYY-MM-DD` format
- Updated `issue` to be a string
- Replaced `pages` with `start` and `end`
Borrowed from [`ruby-cff`](https://github.com/citation-file-format/ruby-cff)'s fixtures https://github.com/citation-file-format/ruby-cff/tree/ec3aeff98855690d56ac2cd35f7181e5d9ee26f3/test/files, with minor changes:

1. updated `cff-version` to "1.2.0"
1. made `version` a string

# Software with a doi (expanded)

Same as [../software-with-a-doi/README.md](../software-with-a-doi/README.md), but with additional properties.
# software with reference of type edited-work

Note that the editors of the edited work must be specified under the `authors`
key. Specific citation styles may or may not attach a suffix to the authors,
such as ", eds." or similar.
Borrowed from [`ruby-cff`](https://github.com/citation-file-format/ruby-cff)'s fixtures https://github.com/citation-file-format/ruby-cff/tree/ec3aeff98855690d56ac2cd35f7181e5d9ee26f3/test/files.For software without a DOI, it is recommended that "the metadata should still
provide information on how to access the specific software, but this may be a
company’s product number or a link to a website that allows the software be
purchased." [Smith et al., 2016](https://doi.org/10.7717/peerj-cs.86), p. 13. Furthermore, "if the version number and
release date are not available, the download date can be used. Similarly, the
contact name/email is an alternative to the location/repository."
([Smith et al., 2016](https://doi.org/10.7717/peerj-cs.86), p. 7).

Hence, for closed-source software without a DOI for which the version number
and release date cannot be determined, a `CITATION.cff` file could look like
[./CITATION.cff](./CITATION.cff).
# Software with a further reference

Where authors wish to encourage citation of an outline paper with citation
of their software, we recommend the use of [reference keys](#references-optional)
to highlight the existence of further references.
Borrowed from [`ruby-cff`](https://github.com/citation-file-format/ruby-cff)'s fixtures https://github.com/citation-file-format/ruby-cff/tree/ec3aeff98855690d56ac2cd35f7181e5d9ee26f3/test/files.We recognize that there are certain situations where it may not be possible to
follow the recommended best-practice. For example, if (1) the software authors
did not register a DOI and/or release a specific version, or (2) the version of
the software used does not match what is available to cite. In those cases,
falling back on a combination of the repository URL and version number/commit
hash would be an appropriate way to cite the software used ([Smith et al., 2016](https://doi.org/10.7717/peerj-cs.86), p. 12).
Borrowed from [`ruby-cff`](https://github.com/citation-file-format/ruby-cff)'s fixtures https://github.com/citation-file-format/ruby-cff/tree/ec3aeff98855690d56ac2cd35f7181e5d9ee26f3/test/files, with some changes.
Borrowed from [`ruby-cff`](https://github.com/citation-file-format/ruby-cff)'s fixtures https://github.com/citation-file-format/ruby-cff/tree/ec3aeff98855690d56ac2cd35f7181e5d9ee26f3/test/files.Borrowed from [`ruby-cff`](https://github.com/citation-file-format/ruby-cff)'s fixtures https://github.com/citation-file-format/ruby-cff/tree/ec3aeff98855690d56ac2cd35f7181e5d9ee26f3/test/files.

- Updated `cff-version` to `"1.2.0"`
- Updated `version` to be a string

Borrowed from [`ruby-cff`](https://github.com/citation-file-format/ruby-cff)'s fixtures https://github.com/citation-file-format/ruby-cff/tree/ec3aeff98855690d56ac2cd35f7181e5d9ee26f3/test/files.