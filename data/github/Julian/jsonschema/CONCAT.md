# Security Policy

## Supported Versions

In general, only the latest released ``jsonschema`` version is supported
and will receive updates.

## Reporting a Vulnerability

To report a security vulnerability, please send an email to
``Julian+Security@GrayVines.com`` with subject line ``SECURITY
(jsonschema)``.

I will do my best to respond within 48 hours to acknowledge the message
and discuss further steps.

If the vulnerability is accepted, an advisory will be sent out via
GitHub's security advisory functionality.

For non-sensitive discussion related to this policy itself, feel free to
open an issue on the issue tracker.
# JSON Schema Test Suite 
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](https://github.com/json-schema-org/.github/blob/main/CODE_OF_CONDUCT.md)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Financial Contributors on Open Collective](https://opencollective.com/json-schema/all/badge.svg?label=financial+contributors)](https://opencollective.com/json-schema)

[![Build Status](https://github.com/json-schema-org/JSON-Schema-Test-Suite/workflows/Test%20Suite%20Sanity%20Checking/badge.svg)](https://github.com/json-schema-org/JSON-Schema-Test-Suite/actions?query=workflow%3A%22Test+Suite+Sanity+Checking%22)

This repository contains a set of JSON objects that implementors of JSON Schema
validation libraries can use to test their validators.

It is meant to be language agnostic and should require only a JSON parser.

The conversion of the JSON objects into tests within your test framework of
choice is still the job of the validator implementor.

## Structure of a Test

The tests in this suite are contained in the `tests` directory at the root of
this repository. Inside that directory is a subdirectory for each draft or
version of the specification.

Inside each draft directory, there are a number of `.json` files and one or more
special subdirectories. The subdirectories contain `.json` files meant for a
specific testing purpose, and each `.json` file logically groups a set of test
cases together. Often the grouping is by property under test, but not always.

The subdirectories are described in the next section.

Inside each `.json` file is a single array containing objects. It's easiest to
illustrate the structure of these with an example:

```json
{
    "description": "The description of the test case",
    "schema": {
        "description": "The schema against which the data in each test is validated",
        "type": "string"
    },
    "tests": [
        {
            "description": "Test for a valid instance",
            "data": "the instance to validate",
            "valid": true
        },
        {
            "description": "Test for an invalid instance",
            "data": 15,
            "valid": false
        }
    ]
}
```

In short: a description, a schema under test, and some tests, where each test
in the `tests` array is an objects with a description of the case itself, the
instance under test, and a boolean indicating whether it should be valid
or invalid.

## Test Subdirectories

There is currently only one subdirectory that may exist within each draft
directory. This is:

1. `optional/`: Contains tests that are considered optional.

## Coverage

All JSON Schema specification releases should be well covered by this suite,
including drafts 2020-12, 2019-09, 07, 06, 04 and 03. Additional coverage is
always welcome, particularly for bugs encountered in real-world
implementations.

Drafts 04 and 03 are considered "frozen" in that less effort is put in to
backport new tests to these versions.

Contributions are very welcome, especially from implementers as they add support
to their own implementations.

If you see anything missing from the current supported drafts, or incorrect on
any draft still accepting bug fixes, please
[file an issue](https://github.com/json-schema-org/JSON-Schema-Test-Suite/issues)
or [submit a PR](https://github.com/json-schema-org/JSON-Schema-Test-Suite).

## Who Uses the Test Suite

This suite is being used by:

### Clojure

* [jinx](https://github.com/juxt/jinx)
* [json-schema](https://github.com/tatut/json-schema)

### Coffeescript

* [jsck](https://github.com/pandastrike/jsck)

### Common Lisp

* [json-schema](https://github.com/fisxoj/json-schema)

### C++

* [Modern C++ JSON schema validator](https://github.com/pboettch/json-schema-validator)

### Dart

* [json\_schema](https://github.com/patefacio/json_schema)

### Elixir

* [ex\_json\_schema](https://github.com/jonasschmidt/ex_json_schema)

### Erlang

* [jesse](https://github.com/for-GET/jesse)

### Go

* [gojsonschema](https://github.com/sigu-399/gojsonschema)
* [validate-json](https://github.com/cesanta/validate-json)

### Haskell

* [aeson-schema](https://github.com/timjb/aeson-schema)
* [hjsonschema](https://github.com/seagreen/hjsonschema)

### Java

* [json-schema-validator](https://github.com/daveclayton/json-schema-validator)
* [everit-org/json-schema](https://github.com/everit-org/json-schema)
* [networknt/json-schema-validator](https://github.com/networknt/json-schema-validator)
* [Justify](https://github.com/leadpony/justify)
* [Snow](https://github.com/ssilverman/snowy-json)
* [jsonschemafriend](https://github.com/jimblackler/jsonschemafriend)

### JavaScript

* [json-schema-benchmark](https://github.com/Muscula/json-schema-benchmark)
* [direct-schema](https://github.com/IreneKnapp/direct-schema)
* [is-my-json-valid](https://github.com/mafintosh/is-my-json-valid)
* [jassi](https://github.com/iclanzan/jassi)
* [JaySchema](https://github.com/natesilva/jayschema)
* [json-schema-valid](https://github.com/ericgj/json-schema-valid)
* [Jsonary](https://github.com/jsonary-js/jsonary)
* [jsonschema](https://github.com/tdegrunt/jsonschema)
* [request-validator](https://github.com/bugventure/request-validator)
* [skeemas](https://github.com/Prestaul/skeemas)
* [tv4](https://github.com/geraintluff/tv4)
* [z-schema](https://github.com/zaggino/z-schema)
* [jsen](https://github.com/bugventure/jsen)
* [ajv](https://github.com/epoberezkin/ajv)
* [djv](https://github.com/korzio/djv)

### Node.js

For node.js developers, the suite is also available as an
[npm](https://www.npmjs.com/package/@json-schema-org/tests) package.

Node-specific support is maintained in a [separate
repository](https://github.com/json-schema-org/json-schema-test-suite-npm)
which also welcomes your contributions!

### .NET

* [Newtonsoft.Json.Schema](https://github.com/JamesNK/Newtonsoft.Json.Schema)
* [Manatee.Json](https://github.com/gregsdennis/Manatee.Json)

### Perl

* [JSON::Schema::Draft201909](https://github.com/karenetheridge/JSON-Schema-Draft201909)
* [JSON::Schema::Tiny](https://github.com/karenetheridge/JSON-Schema-Tiny)
* [Test::JSON::Schema::Acceptance](https://github.com/karenetheridge/Test-JSON-Schema-Acceptance)

### PHP

* [opis/json-schema](https://github.com/opis/json-schema)
* [json-schema](https://github.com/justinrainbow/json-schema)
* [json-guard](https://github.com/thephpleague/json-guard)

### PostgreSQL

* [postgres-json-schema](https://github.com/gavinwahl/postgres-json-schema)
* [is\_jsonb\_valid](https://github.com/furstenheim/is_jsonb_valid)

### Python

* [jsonschema](https://github.com/Julian/jsonschema)
* [fastjsonschema](https://github.com/seznam/python-fastjsonschema)
* [hypothesis-jsonschema](https://github.com/Zac-HD/hypothesis-jsonschema)
* [jschon](https://github.com/marksparkza/jschon)

### Ruby

* [json-schema](https://github.com/hoxworth/json-schema)
* [json\_schemer](https://github.com/davishmcclurg/json_schemer)

### Rust

* [jsonschema](https://github.com/Stranger6667/jsonschema-rs)
* [valico](https://github.com/rustless/valico)

### Swift

* [JSONSchema](https://github.com/kylef/JSONSchema.swift)

If you use it as well, please fork and send a pull request adding yourself to
the list :).

## Contributing

If you see something missing or incorrect, a pull request is most welcome!

There are some sanity checks in place for testing the test suite. You can run
them with `bin/jsonschema_suite check` or `tox`. They will be run automatically
by [GitHub Actions](https://github.com/json-schema-org/JSON-Schema-Test-Suite/actions?query=workflow%3A%22Test+Suite+Sanity+Checking%22)
as well.
v4.4.0
------

* Add ``mypy`` support (#892)
* Add support for Python 3.11

v4.3.3
------

* Properly report deprecation warnings at the right stack level (#899)

v4.3.2
------

* Additional performance improvements for resolving refs (#896)

v4.3.1
------

* Resolving refs has had performance improvements (#893)

v4.3.0
------

* Fix undesired fallback to brute force container uniqueness check on
  certain input types (#893)
* Implement a PEP544 Protocol for validator classes (#890)

v4.2.1
------

* Pin ``importlib.resources`` from below (#877)

v4.2.0
------

* Use ``importlib.resources`` to load schemas (#873)
* Ensure all elements of arrays are verified for uniqueness by ``uniqueItems``
  (#866)

v4.1.2
------

* Fix ``dependentSchemas`` to properly consider non-object instances to be
  valid (#850)

v4.1.1
------

* Fix ``prefixItems`` not indicating which item was invalid within the instance
  path (#862)

v4.1.0
------

* Add Python 3.10 to the list of supported Python versions

v4.0.1
------

* Fix the declaration of minimum supported Python version (#846)

v4.0.0
------

* Partial support for Draft 2020-12 (as well as 2019-09).
  Thanks to Thomas Schmidt and Harald Nezbeda.
* ``False`` and ``0`` are now properly considered non-equal even
  recursively within a container (#686). As part of this change,
  ``uniqueItems`` validation may be *slower* in some cases. Please feel
  free to report any significant performance regressions, though in
  some cases they may be difficult to address given the specification
  requirement.
* The CLI has been improved, and in particular now supports a ``--output``
  option (with ``plain`` (default) or ``pretty`` arguments) to control the
  output format. Future work may add additional machine-parsable output
  formats.
* Code surrounding ``DEFAULT_TYPES`` and the legacy mechanism for
  specifying types to validators have been removed, as per the deprecation
  policy. Validators should use the ``TypeChecker`` object to customize
  the set of Python types corresponding to JSON Schema types.
* Validation errors now have a ``json_path`` attribute, describing their
  location in JSON path format
* Support for the IP address and domain name formats has been improved
* Support for Python 2 and 3.6 has been dropped, with ``python_requires``
  properly set.
* ``multipleOf`` could overflow when given sufficiently large numbers. Now,
  when an overflow occurs, ``jsonschema`` will fall back to using fraction
  division (#746).
* ``jsonschema.__version__``, ``jsonschema.validators.validators``,
  ``jsonschema.validators.meta_schemas`` and
  ``jsonschema.RefResolver.in_scope`` have been deprecated, as has
  passing a second-argument schema to ``Validator.iter_errors`` and
  ``Validator.is_valid``.

v3.2.0
------

* Added a ``format_nongpl`` setuptools extra, which installs only ``format``
  dependencies that are non-GPL (#619).

v3.1.1
------

* Temporarily revert the switch to ``js-regex`` until #611 and #612 are
  resolved.

v3.1.0
------

* Regular expressions throughout schemas now respect the ECMA 262 dialect, as
  recommended by the specification (#609).

v3.0.2
------

* Fixed a bug where ``0`` and ``False`` were considered equal by
  ``const`` and ``enum`` (#575).

v3.0.1
------

* Fixed a bug where extending validators did not preserve their notion
  of which validator property contains ``$id`` information.

v3.0.0
------

* Support for Draft 6 and Draft 7
* Draft 7 is now the default
* New ``TypeChecker`` object for more complex type definitions (and overrides)
* Falling back to isodate for the date-time format checker is no longer
  attempted, in accordance with the specification

v2.6.0
------

* Support for Python 2.6 has been dropped.
* Improve a few error messages for ``uniqueItems`` (#224) and
  ``additionalProperties`` (#317)
* Fixed an issue with ``ErrorTree``'s handling of multiple errors (#288)

v2.5.0
------

* Improved performance on CPython by adding caching around ref resolution
  (#203)

v2.4.0
------

* Added a CLI (#134)
* Added absolute path and absolute schema path to errors (#120)
* Added ``relevance``
* Meta-schemas are now loaded via ``pkgutil``

v2.3.0
------

* Added ``by_relevance`` and ``best_match`` (#91)
* Fixed ``format`` to allow adding formats for non-strings (#125)
* Fixed the ``uri`` format to reject URI references (#131)

v2.2.0
------

* Compile the host name regex (#127)
* Allow arbitrary objects to be types (#129)

v2.1.0
------

* Support RFC 3339 datetimes in conformance with the spec
* Fixed error paths for additionalItems + items (#122)
* Fixed wording for min / maxProperties (#117)


v2.0.0
------

* Added ``create`` and ``extend`` to ``jsonschema.validators``
* Removed ``ValidatorMixin``
* Fixed array indices ref resolution (#95)
* Fixed unknown scheme defragmenting and handling (#102)


v1.3.0
------

* Better error tracebacks (#83)
* Raise exceptions in ``ErrorTree``\s for keys not in the instance (#92)
* __cause__ (#93)


v1.2.0
------

* More attributes for ValidationError (#86)
* Added ``ValidatorMixin.descend``
* Fixed bad ``RefResolutionError`` message (#82)


v1.1.0
------

* Canonicalize URIs (#70)
* Allow attaching exceptions to ``format`` errors (#77)


v1.0.0
------

* Support for Draft 4
* Support for format
* Longs are ints too!
* Fixed a number of issues with ``$ref`` support (#66)
* Draft4Validator is now the default
* ``ValidationError.path`` is now in sequential order
* Added ``ValidatorMixin``


v0.8.0
------

* Full support for JSON References
* ``validates`` for registering new validators
* Documentation
* Bugfixes

    * uniqueItems not so unique (#34)
    * Improper any (#47)


v0.7
----

* Partial support for (JSON Pointer) ``$ref``
* Deprecations

  * ``Validator`` is replaced by ``Draft3Validator`` with a slightly different
    interface
  * ``validator(meta_validate=False)``


v0.6
----

* Bugfixes

  * Issue #30 - Wrong behavior for the dependencies property validation
  * Fixed a miswritten test


v0.5
----

* Bugfixes

  * Issue #17 - require path for error objects
  * Issue #18 - multiple type validation for non-objects


v0.4
----

* Preliminary support for programmatic access to error details (Issue #5).
  There are certainly some corner cases that don't do the right thing yet, but
  this works mostly.

    In order to make this happen (and also to clean things up a bit), a number
    of deprecations are necessary:

        * ``stop_on_error`` is deprecated in ``Validator.__init__``. Use
          ``Validator.iter_errors()`` instead.
        * ``number_types`` and ``string_types`` are deprecated there as well.
          Use ``types={"number" : ..., "string" : ...}`` instead.
        * ``meta_validate`` is also deprecated, and instead is now accepted as
          an argument to ``validate``, ``iter_errors`` and ``is_valid``.

* A bugfix or two


v0.3
----

* Default for unknown types and properties is now to *not* error (consistent
  with the schema).
* Python 3 support
* Removed dependency on SecureTypes now that the hash bug has been resolved.
* "Numerous bug fixes" -- most notably, a divisibleBy error for floats and a
  bunch of missing typechecks for irrelevant properties.
==========
jsonschema
==========

|PyPI| |Pythons| |CI| |ReadTheDocs| |Precommit| |Zenodo|

.. |PyPI| image:: https://img.shields.io/pypi/v/jsonschema.svg
   :alt: PyPI version
   :target: https://pypi.org/project/jsonschema/

.. |Pythons| image:: https://img.shields.io/pypi/pyversions/jsonschema.svg
   :alt: Supported Python versions
   :target: https://pypi.org/project/jsonschema/

.. |CI| image:: https://github.com/Julian/jsonschema/workflows/CI/badge.svg
  :alt: Build status
  :target: https://github.com/Julian/jsonschema/actions?query=workflow%3ACI

.. |ReadTheDocs| image:: https://readthedocs.org/projects/python-jsonschema/badge/?version=stable&style=flat
   :alt: ReadTheDocs status
   :target: https://python-jsonschema.readthedocs.io/en/stable/

.. |Precommit| image:: https://results.pre-commit.ci/badge/github/Julian/jsonschema/main.svg
   :alt: pre-commit.ci status
   :target: https://results.pre-commit.ci/latest/github/Julian/jsonschema/main

.. |Zenodo| image:: https://zenodo.org/badge/3072629.svg
   :target: https://zenodo.org/badge/latestdoi/3072629


``jsonschema`` is an implementation of the `JSON Schema
<https://json-schema.org>`_ specification for Python.

.. code-block:: python

    >>> from jsonschema import validate

    >>> # A sample schema, like what we'd get from json.load()
    >>> schema = {
    ...     "type" : "object",
    ...     "properties" : {
    ...         "price" : {"type" : "number"},
    ...         "name" : {"type" : "string"},
    ...     },
    ... }

    >>> # If no exception is raised by validate(), the instance is valid.
    >>> validate(instance={"name" : "Eggs", "price" : 34.99}, schema=schema)

    >>> validate(
    ...     instance={"name" : "Eggs", "price" : "Invalid"}, schema=schema,
    ... )                                   # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
        ...
    ValidationError: 'Invalid' is not of type 'number'

It can also be used from console:

.. code-block:: bash

    $ jsonschema --instance sample.json sample.schema

Features
--------

* Partial support for
  `Draft 2020-12 <https://python-jsonschema.readthedocs.io/en/latest/validate/#jsonschema.Draft202012Validator>`_ and
  `Draft 2019-09 <https://python-jsonschema.readthedocs.io/en/latest/validate/#jsonschema.Draft201909Validator>`_,
  except for ``dynamicRef`` / ``recursiveRef`` and ``$vocabulary`` (in-progress).
  Full support for
  `Draft 7 <https://python-jsonschema.readthedocs.io/en/latest/validate/#jsonschema.Draft7Validator>`_,
  `Draft 6 <https://python-jsonschema.readthedocs.io/en/latest/validate/#jsonschema.Draft6Validator>`_,
  `Draft 4 <https://python-jsonschema.readthedocs.io/en/latest/validate/#jsonschema.Draft4Validator>`_
  and
  `Draft 3 <https://python-jsonschema.readthedocs.io/en/latest/validate/#jsonschema.Draft3Validator>`_

* `Lazy validation <https://python-jsonschema.readthedocs.io/en/latest/validate/#jsonschema.protocols.Validator.iter_errors>`_
  that can iteratively report *all* validation errors.

* `Programmatic querying <https://python-jsonschema.readthedocs.io/en/latest/errors/>`_
  of which properties or items failed validation.


Installation
------------

``jsonschema`` is available on `PyPI <https://pypi.org/project/jsonschema/>`_. You can install using `pip <https://pip.pypa.io/en/stable/>`_:

.. code-block:: bash

    $ pip install jsonschema


Running the Test Suite
----------------------

If you have ``tox`` installed (perhaps via ``pip install tox`` or your
package manager), running ``tox`` in the directory of your source
checkout will run ``jsonschema``'s test suite on all of the versions
of Python ``jsonschema`` supports. If you don't have all of the
versions that ``jsonschema`` is tested under, you'll likely want to run
using ``tox``'s ``--skip-missing-interpreters`` option.

Of course you're also free to just run the tests on a single version with your
favorite test runner. The tests live in the ``jsonschema.tests`` package.


Benchmarks
----------

``jsonschema``'s benchmarks make use of `pyperf
<https://pyperf.readthedocs.io>`_. Running them can be done via::

      $ tox -e perf


Community
---------

The JSON Schema specification has `a Slack
<https://json-schema.slack.com>`_, with an `invite link on its home page
<https://json-schema.org/>`_. Many folks knowledgeable on authoring
schemas can be found there.

Otherwise, asking questions on Stack Overflow is another means of
getting help if you're stuck.

Contributing
------------

I'm Julian Berman.

``jsonschema`` is on `GitHub <https://github.com/Julian/jsonschema>`_.

Get in touch, via GitHub or otherwise, if you've got something to contribute,
it'd be most welcome!

You can also generally find me on Libera (nick: ``Julian``) in various
channels, including ``#python``.

If you feel overwhelmingly grateful, you can also `sponsor me
<https://github.com/sponsors/Julian/>`_.

And for companies who appreciate ``jsonschema`` and its continued support
and growth, ``jsonschema`` is also now supportable via `TideLift
<https://tidelift.com/subscription/pkg/pypi-jsonschema?utm_source=pypi-j
sonschema&utm_medium=referral&utm_campaign=readme>`_.
=================
Schema Validation
=================


.. currentmodule:: jsonschema


The Basics
----------

The simplest way to validate an instance under a given schema is to use the
:func:`validate` function.

.. autofunction:: validate

.. [#] For information on creating JSON schemas to validate
    your data, there is a good introduction to JSON Schema
    fundamentals underway at `Understanding JSON Schema
    <https://json-schema.org/understanding-json-schema/>`_

.. _validator-protocol:

The Validator Protocol
-----------------------

`jsonschema` defines a protocol that all validator
classes should adhere to.

.. autoclass:: jsonschema.protocols.Validator
    :members:

All of the `versioned validators <versioned-validators>` that are included with
`jsonschema` adhere to the protocol, and implementers of validator classes
that extend or complement the ones included should adhere to it as well. For
more information see `creating-validators`.

Type Checking
-------------

To handle JSON Schema's :validator:`type` property, a `Validator` uses
an associated `TypeChecker`. The type checker provides an immutable
mapping between names of types and functions that can test if an instance is
of that type. The defaults are suitable for most users - each of the
`versioned validators <versioned-validators>` that are included with
`jsonschema` have a `TypeChecker` that can correctly handle their respective
versions.

.. seealso:: `validating-types`

    For an example of providing a custom type check.

.. autoclass:: TypeChecker
    :members:

.. autoexception:: jsonschema.exceptions.UndefinedTypeCheck

    Raised when trying to remove a type check that is not known to this
    TypeChecker, or when calling `jsonschema.TypeChecker.is_type`
    directly.

.. _validating-types:

Validating With Additional Types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Occasionally it can be useful to provide additional or alternate types when
validating the JSON Schema's :validator:`type` property.

`jsonschema` tries to strike a balance between performance in the common
case and generality. For instance, JSON Schema defines a ``number`` type, which
can be validated with a schema such as ``{"type" : "number"}``. By default,
this will accept instances of Python `numbers.Number`. This includes in
particular `int`\s and `float`\s, along with
`decimal.Decimal` objects, `complex` numbers etc. For
``integer`` and ``object``, however, rather than checking for
`numbers.Integral` and `collections.abc.Mapping`,
`jsonschema` simply checks for `int` and `dict`, since the
more general instance checks can introduce significant slowdown, especially
given how common validating these types are.

If you *do* want the generality, or just want to add a few specific additional
types as being acceptable for a validator object, then you should update an
existing `TypeChecker` or create a new one. You may then create a new
`Validator` via `jsonschema.validators.extend`.

.. code-block:: python

    class MyInteger(object):
        pass

    def is_my_int(checker, instance):
        return (
            Draft3Validator.TYPE_CHECKER.is_type(instance, "number") or
            isinstance(instance, MyInteger)
        )

    type_checker = Draft3Validator.TYPE_CHECKER.redefine("number", is_my_int)

    CustomValidator = extend(Draft3Validator, type_checker=type_checker)
    validator = CustomValidator(schema={"type" : "number"})


.. autoexception:: jsonschema.exceptions.UnknownType

.. _versioned-validators:

Versioned Validators
--------------------

`jsonschema` ships with validator classes for various versions of
the JSON Schema specification. For details on the methods and attributes
that each validator class provides see the `Validator` protocol,
which each included validator class implements.

.. autoclass:: Draft202012Validator

.. autoclass:: Draft201909Validator

.. autoclass:: Draft7Validator

.. autoclass:: Draft6Validator

.. autoclass:: Draft4Validator

.. autoclass:: Draft3Validator


For example, if you wanted to validate a schema you created against the
Draft 7 meta-schema, you could use:

.. code-block:: python

    from jsonschema import Draft7Validator

    schema = {
        "$schema": "http://json-schema.org/draft-07/schema#",

        "type": "object",
        "properties": {
            "name": {"type": "string"},
            "email": {"type": "string"},
        },
        "required": ["email"]
    }
    Draft7Validator.check_schema(schema)


.. _validating formats:

Validating Formats
------------------

JSON Schema defines the :validator:`format` property which can be used to check
if primitive types (``string``\s, ``number``\s, ``boolean``\s) conform to
well-defined formats. By default, no validation is enforced, but optionally,
validation can be enabled by hooking in a format-checking object into an
`Validator`.

.. doctest::

    >>> validate("127.0.0.1", {"format" : "ipv4"})
    >>> validate(
    ...     instance="-12",
    ...     schema={"format" : "ipv4"},
    ...     format_checker=draft7_format_checker,
    ... )
    Traceback (most recent call last):
        ...
    ValidationError: "-12" is not a "ipv4"

.. autoclass:: FormatChecker
    :members:
    :exclude-members: cls_checks

    .. attribute:: checkers

        A mapping of currently known formats to tuple of functions that
        validate them and errors that should be caught. New checkers can be
        added and removed either per-instance or globally for all checkers
        using the `FormatChecker.checks` or `FormatChecker.cls_checks`
        decorators respectively.

    .. classmethod:: cls_checks(format, raises=())

        Register a decorated function as *globally* validating a new format.

        Any instance created after this function is called will pick up the
        supplied checker.

        :argument str format: the format that the decorated function will check
        :argument Exception raises: the exception(s) raised
            by the decorated function when an invalid instance is
            found. The exception object will be accessible as the
            `jsonschema.exceptions.ValidationError.cause` attribute
            of the resulting validation error.


.. autoexception:: FormatError
    :members:


There are a number of default checkers that `FormatChecker`\s know how
to validate. Their names can be viewed by inspecting the
`FormatChecker.checkers` attribute. Certain checkers will only be
available if an appropriate package is available for use. The easiest way to
ensure you have what is needed is to install ``jsonschema`` using the
``format`` or ``format_nongpl`` setuptools extra -- i.e.

.. code-block:: sh

   $ pip install jsonschema[format]

which will install all of the below dependencies for all formats.

Or if you want to install MIT-license compatible dependencies only:

.. code-block:: sh

   $ pip install jsonschema[format_nongpl]

The non-GPL extra is intended to not install any direct dependencies
that are GPL (but that of course end-users should do their own verification).
At the moment, it supports all the available checkers except for ``iri`` and
``iri-reference``.

The more specific list of available checkers, along with their requirement
(if any,) are listed below.

.. note::

    If the following packages are not installed when using a checker
    that requires it, validation will succeed without throwing an error,
    as specified by the JSON Schema specification.

=========================  ====================
Checker                    Notes
=========================  ====================
``color``                  requires webcolors_
``date``
``date-time``              requires rfc3339-validator_
``duration``               requires isoduration_
``email``
``hostname``               requires fqdn_
``idn-hostname``           requires idna_
``ipv4``
``ipv6``                   OS must have `socket.inet_pton` function
``iri``                    requires rfc3987_
``iri-reference``          requires rfc3987_
``json-pointer``           requires jsonpointer_
``regex``
``relative-json-pointer``  requires jsonpointer_
``time``                   requires rfc3339-validator_
``uri``                    requires rfc3987_ or rfc3986-validator_
``uri-reference``          requires rfc3987_ or rfc3986-validator_
``uri-template``           requires uri-template_
=========================  ====================


.. _fqdn: https://pypi.org/pypi/fqdn/
.. _idna: https://pypi.org/pypi/idna/
.. _isoduration: https://pypi.org/pypi/isoduration/
.. _jsonpointer: https://pypi.org/pypi/jsonpointer/
.. _rfc3339-validator: https://pypi.org/project/rfc3339-validator/
.. _rfc3986-validator: https://pypi.org/project/rfc3986-validator/
.. _rfc3987: https://pypi.org/pypi/rfc3987/
.. _rfc5322: https://tools.ietf.org/html/rfc5322#section-3.4.1
.. _uri-template: https://pypi.org/pypi/uri-template/
.. _webcolors: https://pypi.org/pypi/webcolors/

.. note::

    Since in most cases "validating" an email address is an attempt
    instead to confirm that mail sent to it will deliver to a recipient,
    and that that recipient is the correct one the email is intended
    for, and since many valid email addresses are in many places
    incorrectly rejected, and many invalid email addresses are in many
    places incorrectly accepted, the ``email`` format validator only
    provides a sanity check, not full rfc5322_ validation.

    The same applies to the ``idn-email`` format.
.. currentmodule:: jsonschema.validators

.. _creating-validators:

=======================================
Creating or Extending Validator Classes
=======================================

.. autofunction:: create

.. autofunction:: extend

.. autofunction:: validator_for

.. autofunction:: validates


Creating Validation Errors
--------------------------

Any validating function that validates against a subschema should call
``descend``, rather than ``iter_errors``. If it recurses into the
instance, or schema, it should pass one or both of the ``path`` or
``schema_path`` arguments to ``descend`` in order to properly maintain
where in the instance or schema respectively the error occurred.

The Validator Protocol
----------------------

``jsonschema`` defines a `protocol <typing.Protocol>`,
`jsonschema.protocols.Validator` which can be used in type annotations to
describe the type of a validator object.

For full details, see `validator-protocol`.
==========================
Frequently Asked Questions
==========================


My schema specifies format validation. Why do invalid instances seem valid?
---------------------------------------------------------------------------

The :validator:`format` validator can be a bit of a stumbling block for new
users working with JSON Schema.

In a schema such as:

.. code-block:: json

    {"type": "string", "format": "date"}

JSON Schema specifications have historically differentiated between the
:validator:`format` validator and other validators. In particular, the
:validator:`format` validator was specified to be *informational* as much
as it may be used for validation.

In other words, for many use cases, schema authors may wish to use
values for the :validator:`format` validator but have no expectation
they be validated alongside other required assertions in a schema.

Of course this does not represent all or even most use cases -- many
schema authors *do* wish to assert that instances conform fully, even to
the specific format mentioned.

In drafts prior to ``draft2019-09``, the decision on whether to
automatically enable :validator:`format` validation was left up to
validation implementations such as this one.

This library made the choice to leave it off by default, for two reasons:

    * for forward compatibility and implementation complexity reasons
      -- if :validator:`format` validation were on by default, and a
      future draft of JSON Schema introduced a hard-to-implement format,
      either the implementation of that format would block releases of
      this library until it were implemented, or the behavior surrounding
      :validator:`format` would need to be even more complex than simply
      defaulting to be on. It therefore was safer to start with it off,
      and defend against the expectation that a given format would always
      automatically work.

    * given that a common use of JSON Schema is for portability across
      languages (and therefore implementations of JSON Schema), so that
      users be aware of this point itself regarding :validator:`format`
      validation, and therefore remember to check any *other*
      implementations they were using to ensure they too were explicitly
      enabled for :validator:`format` validation.

As of ``draft2019-09`` however, the opt-out by default behavior
mentioned here is now *required* for all validators.

Difficult as this may sound for new users, at this point it at least
means they should expect the same behavior that has always been
implemented here, across any other implementation they encounter.

.. seealso::

    `Draft 2019-09's release notes on format <https://json-schema.org/draft/2019-09/release-notes.html#format-vocabulary>`_

        for upstream details on the behavior of format and how it has changed
        in ``draft2019-09``

    `validating formats`

        for details on how to enable format validation

    `jsonschema.FormatChecker`

        the object which implements format validation


Why doesn't my schema's default property set the default on my instance?
------------------------------------------------------------------------

The basic answer is that the specification does not require that
:validator:`default` actually do anything.

For an inkling as to *why* it doesn't actually do anything, consider
that none of the other validators modify the instance either. More
importantly, having :validator:`default` modify the instance can produce
quite peculiar things. It's perfectly valid (and perhaps even useful)
to have a default that is not valid under the schema it lives in! So an
instance modified by the default would pass validation the first time,
but fail the second!

Still, filling in defaults is a thing that is useful. `jsonschema`
allows you to `define your own validator classes and callables
<creating>`, so you can easily create an `jsonschema.protocols.Validator`
that does do default setting. Here's some code to get you started. (In
this code, we add the default properties to each object *before* the
properties are validated, so the default values themselves will need to
be valid under the schema.)

    .. code-block:: python

        from jsonschema import Draft7Validator, validators


        def extend_with_default(validator_class):
            validate_properties = validator_class.VALIDATORS["properties"]

            def set_defaults(validator, properties, instance, schema):
                for property, subschema in properties.items():
                    if "default" in subschema:
                        instance.setdefault(property, subschema["default"])

                for error in validate_properties(
                    validator, properties, instance, schema,
                ):
                    yield error

            return validators.extend(
                validator_class, {"properties" : set_defaults},
            )


        DefaultValidatingDraft7Validator = extend_with_default(Draft7Validator)


        # Example usage:
        obj = {}
        schema = {'properties': {'foo': {'default': 'bar'}}}
        # Note jsonschem.validate(obj, schema, cls=DefaultValidatingDraft7Validator)
        # will not work because the metaschema contains `default` directives.
        DefaultValidatingDraft7Validator(schema).validate(obj)
        assert obj == {'foo': 'bar'}


See the above-linked document for more info on how this works, but
basically, it just extends the :validator:`properties` validator on
a `jsonschema.Draft7Validator` to then go ahead and update all the
defaults.

.. note::

    If you're interested in a more interesting solution to a larger
    class of these types of transformations, keep an eye on `Seep
    <https://github.com/Julian/Seep>`_, which is an experimental
    data transformation and extraction library written on top of
    `jsonschema`.


.. hint::

    The above code can provide default values for an entire object and
    all of its properties, but only if your schema provides a default
    value for the object itself, like so:

    .. code-block:: python

        schema = {
            "type": "object",
            "properties": {
                "outer-object": {
                    "type": "object",
                    "properties" : {
                        "inner-object": {
                            "type": "string",
                            "default": "INNER-DEFAULT"
                        }
                    },
                    "default": {} # <-- MUST PROVIDE DEFAULT OBJECT
                }
            }
        }

        obj = {}
        DefaultValidatingDraft7Validator(schema).validate(obj)
        assert obj == {'outer-object': {'inner-object': 'INNER-DEFAULT'}}

    ...but if you don't provide a default value for your object, then
    it won't be instantiated at all, much less populated with default
    properties.

    .. code-block:: python

        del schema["properties"]["outer-object"]["default"]
        obj2 = {}
        DefaultValidatingDraft7Validator(schema).validate(obj2)
        assert obj2 == {} # whoops


How do jsonschema version numbers work?
---------------------------------------

``jsonschema`` tries to follow the `Semantic Versioning
<https://semver.org/>`_ specification.

This means broadly that no backwards-incompatible changes should be made
in minor releases (and certainly not in dot releases).

The full picture requires defining what constitutes a
backwards-incompatible change.

The following are simple examples of things considered public API,
and therefore should *not* be changed without bumping a major version
number:

    * module names and contents, when not marked private by Python
      convention (a single leading underscore)

    * function and object signature (parameter order and name)

The following are *not* considered public API and may change without
notice:

    * the exact wording and contents of error messages; typical reasons
      to rely on this seem to involve downstream tests in packages using
      `jsonschema`. These use cases are encouraged to use the extensive
      introspection provided in `jsonschema.exceptions.ValidationError`\s
      instead to make meaningful assertions about what failed rather than
      relying on *how* what failed is explained to a human.

    * the order in which validation errors are returned or raised

    * the contents of the ``jsonschema.tests`` package

    * the contents of the ``jsonschema.benchmarks`` package

    * the specific non-zero error codes presented by the command line
      interface

    * the exact representation of errors presented by the command line
      interface, other than that errors represented by the plain outputter
      will be reported one per line

    * anything marked private

With the exception of the last two of those, flippant changes are
avoided, but changes can and will be made if there is improvement to be
had. Feel free to open an issue ticket if there is a specific issue or
question worth raising.
.. module:: jsonschema
.. include:: ../README.rst


Contents
--------

.. toctree::
    :maxdepth: 2

    validate
    errors
    references
    creating
    faq


Indices and tables
==================

* `genindex`
* `search`
=========================
Resolving JSON References
=========================


.. currentmodule:: jsonschema

.. autoclass:: RefResolver
    :members:

.. autoexception:: RefResolutionError

    A JSON reference failed to resolve.
==========================
Handling Validation Errors
==========================

.. currentmodule:: jsonschema.exceptions

When an invalid instance is encountered, a `ValidationError` will be
raised or returned, depending on which method or function is used.

.. autoexception:: ValidationError

    The information carried by an error roughly breaks down into:

    ===============  =================  ========================
     What Happened   Why Did It Happen  What Was Being Validated
    ===============  =================  ========================
    `message`        `context`          `instance`

                     `cause`            `json_path`

                                        `path`

                                        `schema`

                                        `schema_path`

                                        `validator`

                                        `validator_value`
    ===============  =================  ========================


    .. attribute:: message

        A human readable message explaining the error.

    .. attribute:: validator

        The name of the failed `validator
        <https://json-schema.org/draft-04/json-schema-validation.html#rfc.section.5>`_.

    .. attribute:: validator_value

        The value for the failed validator in the schema.

    .. attribute:: schema

        The full schema that this error came from. This is potentially a
        subschema from within the schema that was passed in originally,
        or even an entirely different schema if a :validator:`$ref` was
        followed.

    .. attribute:: relative_schema_path

        A `collections.deque` containing the path to the failed
        validator within the schema.

    .. attribute:: absolute_schema_path

        A `collections.deque` containing the path to the failed
        validator within the schema, but always relative to the
        *original* schema as opposed to any subschema (i.e. the one
        originally passed into a validator class, *not* `schema`\).

    .. attribute:: schema_path

        Same as `relative_schema_path`.

    .. attribute:: relative_path

        A `collections.deque` containing the path to the
        offending element within the instance. The deque can be empty if
        the error happened at the root of the instance.

    .. attribute:: absolute_path

        A `collections.deque` containing the path to the
        offending element within the instance. The absolute path
        is always relative to the *original* instance that was
        validated (i.e. the one passed into a validation method, *not*
        `instance`\). The deque can be empty if the error happened
        at the root of the instance.

    .. attribute:: json_path

        A `JSON path <https://goessner.net/articles/JsonPath/index.html>`_
        to the offending element within the instance.

    .. attribute:: path

        Same as `relative_path`.

    .. attribute:: instance

        The instance that was being validated. This will differ from
        the instance originally passed into ``validate`` if the
        validator object was in the process of validating a (possibly
        nested) element within the top-level instance. The path within
        the top-level instance (i.e. `ValidationError.path`) could
        be used to find this object, but it is provided for convenience.

    .. attribute:: context

        If the error was caused by errors in subschemas, the list of errors
        from the subschemas will be available on this property. The
        `schema_path` and `path` of these errors will be relative
        to the parent error.

    .. attribute:: cause

        If the error was caused by a *non*-validation error, the
        exception object will be here. Currently this is only used
        for the exception raised by a failed format checker in
        `jsonschema.FormatChecker.check`.

    .. attribute:: parent

        A validation error which this error is the `context` of.
        ``None`` if there wasn't one.


In case an invalid schema itself is encountered, a `SchemaError` is
raised.

.. autoexception:: SchemaError

    The same attributes are present as for `ValidationError`\s.


These attributes can be clarified with a short example:

.. testcode::

    schema = {
        "items": {
            "anyOf": [
                {"type": "string", "maxLength": 2},
                {"type": "integer", "minimum": 5}
            ]
        }
    }
    instance = [{}, 3, "foo"]
    v = Draft7Validator(schema)
    errors = sorted(v.iter_errors(instance), key=lambda e: e.path)

The error messages in this situation are not very helpful on their own.

.. testcode::

    for error in errors:
        print(error.message)

outputs:

.. testoutput::

    {} is not valid under any of the given schemas
    3 is not valid under any of the given schemas
    'foo' is not valid under any of the given schemas

If we look at `ValidationError.path` on each of the errors, we can find
out which elements in the instance correspond to each of the errors. In
this example, `ValidationError.path` will have only one element, which
will be the index in our list.

.. testcode::

    for error in errors:
        print(list(error.path))

.. testoutput::

    [0]
    [1]
    [2]

Since our schema contained nested subschemas, it can be helpful to look at
the specific part of the instance and subschema that caused each of the errors.
This can be seen with the `ValidationError.instance` and
`ValidationError.schema` attributes.

With validators like :validator:`anyOf`, the `ValidationError.context`
attribute can be used to see the sub-errors which caused the failure. Since
these errors actually came from two separate subschemas, it can be helpful to
look at the `ValidationError.schema_path` attribute as well to see where
exactly in the schema each of these errors come from. In the case of sub-errors
from the `ValidationError.context` attribute, this path will be relative
to the `ValidationError.schema_path` of the parent error.

.. testcode::

    for error in errors:
        for suberror in sorted(error.context, key=lambda e: e.schema_path):
            print(list(suberror.schema_path), suberror.message, sep=", ")

.. testoutput::

    [0, 'type'], {} is not of type 'string'
    [1, 'type'], {} is not of type 'integer'
    [0, 'type'], 3 is not of type 'string'
    [1, 'minimum'], 3 is less than the minimum of 5
    [0, 'maxLength'], 'foo' is too long
    [1, 'type'], 'foo' is not of type 'integer'

The string representation of an error combines some of these attributes for
easier debugging.

.. testcode::

    print(errors[1])

.. testoutput::

    3 is not valid under any of the given schemas

    Failed validating 'anyOf' in schema['items']:
        {'anyOf': [{'maxLength': 2, 'type': 'string'},
                   {'minimum': 5, 'type': 'integer'}]}

    On instance[1]:
        3


ErrorTrees
----------

If you want to programmatically be able to query which properties or validators
failed when validating a given instance, you probably will want to do so using
`jsonschema.exceptions.ErrorTree` objects.

.. autoclass:: jsonschema.exceptions.ErrorTree
    :members:
    :special-members:
    :exclude-members: __dict__,__weakref__

    .. attribute:: errors

        The mapping of validator names to the error objects (usually
        `jsonschema.exceptions.ValidationError`\s) at this level
        of the tree.

Consider the following example:

.. testcode::

    schema = {
        "type" : "array",
        "items" : {"type" : "number", "enum" : [1, 2, 3]},
        "minItems" : 3,
    }
    instance = ["spam", 2]

For clarity's sake, the given instance has three errors under this schema:

.. testcode::

    v = Draft3Validator(schema)
    for error in sorted(v.iter_errors(["spam", 2]), key=str):
        print(error.message)

.. testoutput::

    'spam' is not of type 'number'
    'spam' is not one of [1, 2, 3]
    ['spam', 2] is too short

Let's construct an `jsonschema.exceptions.ErrorTree` so that we
can query the errors a bit more easily than by just iterating over the
error objects.

.. testcode::

    tree = ErrorTree(v.iter_errors(instance))

As you can see, `jsonschema.exceptions.ErrorTree` takes an
iterable of `ValidationError`\s when constructing a tree so
you can directly pass it the return value of a validator object's
`jsonschema.protocols.Validator.iter_errors` method.

`ErrorTree`\s support a number of useful operations. The first one we
might want to perform is to check whether a given element in our instance
failed validation. We do so using the :keyword:`in` operator:

.. doctest::

    >>> 0 in tree
    True

    >>> 1 in tree
    False

The interpretation here is that the 0th index into the instance (``"spam"``)
did have an error (in fact it had 2), while the 1th index (``2``) did not (i.e.
it was valid).

If we want to see which errors a child had, we index into the tree and look at
the `ErrorTree.errors` attribute.

.. doctest::

    >>> sorted(tree[0].errors)
    ['enum', 'type']

Here we see that the :validator:`enum` and :validator:`type` validators failed
for index ``0``. In fact `ErrorTree.errors` is a dict, whose values are
the `ValidationError`\s, so we can get at those directly if we want
them.

.. doctest::

    >>> print(tree[0].errors["type"].message)
    'spam' is not of type 'number'

Of course this means that if we want to know if a given named
validator failed for a given index, we check for its presence in
`ErrorTree.errors`:

.. doctest::

    >>> "enum" in tree[0].errors
    True

    >>> "minimum" in tree[0].errors
    False

Finally, if you were paying close enough attention, you'll notice that we
haven't seen our :validator:`minItems` error appear anywhere yet. This is
because :validator:`minItems` is an error that applies globally to the instance
itself. So it appears in the root node of the tree.

.. doctest::

    >>> "minItems" in tree.errors
    True

That's all you need to know to use error trees.

To summarize, each tree contains child trees that can be accessed by
indexing the tree to get the corresponding child tree for a given index
into the instance. Each tree and child has a `ErrorTree.errors`
attribute, a dict, that maps the failed validator name to the
corresponding validation error.


best_match and relevance
------------------------

The `best_match` function is a simple but useful function for attempting
to guess the most relevant error in a given bunch.

.. doctest::

        >>> from jsonschema import Draft7Validator
        >>> from jsonschema.exceptions import best_match

        >>> schema = {
        ...     "type": "array",
        ...     "minItems": 3,
        ... }
        >>> print(best_match(Draft7Validator(schema).iter_errors(11)).message)
        11 is not of type 'array'


.. autofunction:: best_match


.. function:: relevance(validation_error)

    A key function that sorts errors based on heuristic relevance.

    If you want to sort a bunch of errors entirely, you can use
    this function to do so. Using this function as a key to e.g.
    `sorted` or `max` will cause more relevant errors to be
    considered greater than less relevant ones.

    Within the different validators that can fail, this function
    considers :validator:`anyOf` and :validator:`oneOf` to be *weak*
    validation errors, and will sort them lower than other validators at
    the same level in the instance.

    If you want to change the set of weak [or strong] validators you can create
    a custom version of this function with `by_relevance` and provide a
    different set of each.

.. doctest::

    >>> schema = {
    ...     "properties": {
    ...         "name": {"type": "string"},
    ...         "phones": {
    ...             "properties": {
    ...                 "home": {"type": "string"}
    ...             },
    ...         },
    ...     },
    ... }
    >>> instance = {"name": 123, "phones": {"home": [123]}}
    >>> errors = Draft7Validator(schema).iter_errors(instance)
    >>> [
    ...     e.path[-1]
    ...     for e in sorted(errors, key=exceptions.relevance)
    ... ]
    ['home', 'name']


.. autofunction:: by_relevance
