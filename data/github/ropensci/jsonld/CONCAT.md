# jsonld

> JSON for Linking Data

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/ropensci/jsonld.svg?branch=master)](https://travis-ci.org/ropensci/jsonld)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/jsonld?branch=master&svg=true)](https://ci.appveyor.com/project/jeroen/jsonld)
[![Coverage Status](https://codecov.io/github/ropensci/jsonld/coverage.svg?branch=master)](https://codecov.io/github/ropensci/jsonld?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/jsonld)](https://cran.r-project.org/package=jsonld)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/jsonld)](https://cran.r-project.org/package=jsonld)
[![Github Stars](https://img.shields.io/github/stars/ropensci/jsonld.svg?style=social&label=Github)](https://github.com/ropensci/jsonld)


JSON-LD is a light-weight syntax for expressing linked data. It is primarily
intended for web-based programming environments, interoperable web services and for 
storing linked data in JSON-based databases. This package provides bindings to the 
JavaScript library for converting, expanding and compacting JSON-LD documents.

## Hello World



Example from https://github.com/digitalbazaar/jsonld.js#quick-examples. Example data:


```r
doc <- '{
  "http://schema.org/name": "Manu Sporny",
  "http://schema.org/url": {"@id": "http://manu.sporny.org/"},
  "http://schema.org/image": {"@id": "http://manu.sporny.org/images/manu.png"}
}'

context <- '{
  "name": "http://schema.org/name",
  "homepage": {"@id": "http://schema.org/url", "@type": "@id"},
  "image": {"@id": "http://schema.org/image", "@type": "@id"}
}'
```

### Compact and expand:


```r
(out <- jsonld_compact(doc, context))
```

```
{
  "@context": {
    "name": "http://schema.org/name",
    "homepage": {
      "@id": "http://schema.org/url",
      "@type": "@id"
    },
    "image": {
      "@id": "http://schema.org/image",
      "@type": "@id"
    }
  },
  "image": "http://manu.sporny.org/images/manu.png",
  "name": "Manu Sporny",
  "homepage": "http://manu.sporny.org/"
} 
```

```r
(expanded <- jsonld_expand(out))
```

```
[
  {
    "http://schema.org/url": [
      {
        "@id": "http://manu.sporny.org/"
      }
    ],
    "http://schema.org/image": [
      {
        "@id": "http://manu.sporny.org/images/manu.png"
      }
    ],
    "http://schema.org/name": [
      {
        "@value": "Manu Sporny"
      }
    ]
  }
] 
```

### Convert between JSON and RDF:


```r
cat(nquads <- jsonld_to_rdf(doc))
```

```
_:b0 <http://schema.org/image> <http://manu.sporny.org/images/manu.png> .
_:b0 <http://schema.org/name> "Manu Sporny" .
_:b0 <http://schema.org/url> <http://manu.sporny.org/> .
```

```r
jsonld_from_rdf(nquads)
```

```
[
  {
    "@id": "_:b0",
    "http://schema.org/image": [
      {
        "@id": "http://manu.sporny.org/images/manu.png"
      }
    ],
    "http://schema.org/name": [
      {
        "@value": "Manu Sporny"
      }
    ],
    "http://schema.org/url": [
      {
        "@id": "http://manu.sporny.org/"
      }
    ]
  }
] 
```

### Other utilities:


```r
jsonld_flatten(doc)
```

```
[
  {
    "@id": "_:b0",
    "http://schema.org/image": [
      {
        "@id": "http://manu.sporny.org/images/manu.png"
      }
    ],
    "http://schema.org/name": [
      {
        "@value": "Manu Sporny"
      }
    ],
    "http://schema.org/url": [
      {
        "@id": "http://manu.sporny.org/"
      }
    ]
  }
] 
```

```r
cat(jsonld_normalize(doc, algorithm = 'URDNA2015', format = 'application/nquads'))
```

```
_:c14n0 <http://schema.org/image> <http://manu.sporny.org/images/manu.png> .
_:c14n0 <http://schema.org/name> "Manu Sporny" .
_:c14n0 <http://schema.org/url> <http://manu.sporny.org/> .
```

# R bindings for jsonld.js

Tests at: 2019-02-01 12:30:49 

## Failures for jsonld.compact


### Test: `compact-0038-in.jsonld`

Expected:
```json
{
    "@context": {
        "site": "http://example.com/",
        "site-cd": "site:site-schema/content-deployment/",
        "title": {
            "@id": "site-cd:node/article/title",
            "@container": "@index"
        },
        "body": {
            "@id": "site-cd:node/article/body",
            "@container": "@index"
        },
        "field_tags": {
            "@id": "site-cd:node/article/field_tags",
            "@container": "@index"
        }
    },
    "@id": "site:node/1",
    "@type": "site-cd:node/article",
    "title": {
        "en": {
            "@type": "site-cd:field-types/title_field",
            "title:/value": "This is the English title"
        },
        "es": {
            "@type": "site-cd:field-types/title_field",
            "title:/value": "Este es el t’tulo espa–ol"
        }
    },
    "body": {
        "en": {
            "@type": "site-cd:field-types/text_with_summary",
            "body:/value": "This is the English body. There is no Spanish body, so this will be displayed for both the English and Spanish versions.",
            "body:/summary": "This is the teaser for the body.",
            "body:/format": "full_html"
        }
    },
    "field_tags": {
        "en": {
            "@type": "site-cd:taxonomy/term",
            "@id": "site:taxonomy/term/1",
            "site-cd:taxonomy/term/uuid": "e34b982c-98ac-4862-9b00-fa771a388010"
        },
        "es": [
            {
                "@type": "site-cd:taxonomy/term",
                "@id": "site:taxonomy/term/1",
                "site-cd:taxonomy/term/uuid": "e34b982c-98ac-4862-9b00-fa771a388010"
            },
            {
                "@type": "site-cd:taxonomy/term",
                "@id": "site:taxonomy/term/2",
                "site-cd:taxonomy/term/uuid": "a55b982c-58ac-4862-9b00-aa221a388010"
            }
        ]
    }
}

```

Found output:
```json
{
    "@context": {
        "site": "http://example.com/",
        "site-cd": "site:site-schema/content-deployment/",
        "title": {
            "@id": "site-cd:node/article/title",
            "@container": "@index"
        },
        "body": {
            "@id": "site-cd:node/article/body",
            "@container": "@index"
        },
        "field_tags": {
            "@id": "site-cd:node/article/field_tags",
            "@container": "@index"
        }
    },
    "@id": "site:node/1",
    "@type": "site-cd:node/article",
    "body": {
        "en": {
            "@type": "site-cd:field-types/text_with_summary",
            "site-cd:node/article/body/format": "full_html",
            "site-cd:node/article/body/summary": "This is the teaser for the body.",
            "site-cd:node/article/body/value": "This is the English body. There is no Spanish body, so this will be displayed for both the English and Spanish versions."
        }
    },
    "field_tags": {
        "en": {
            "@id": "site:taxonomy/term/1",
            "@type": "site-cd:taxonomy/term",
            "site-cd:taxonomy/term/uuid": "e34b982c-98ac-4862-9b00-fa771a388010"
        },
        "es": [
            {
                "@id": "site:taxonomy/term/1",
                "@type": "site-cd:taxonomy/term",
                "site-cd:taxonomy/term/uuid": "e34b982c-98ac-4862-9b00-fa771a388010"
            },
            {
                "@id": "site:taxonomy/term/2",
                "@type": "site-cd:taxonomy/term",
                "site-cd:taxonomy/term/uuid": "a55b982c-58ac-4862-9b00-aa221a388010"
            }
        ]
    },
    "title": {
        "en": {
            "@type": "site-cd:field-types/title_field",
            "site-cd:node/article/title/value": "This is the English title"
        },
        "es": {
            "@type": "site-cd:field-types/title_field",
            "site-cd:node/article/title/value": "Este es el t’tulo espa–ol"
        }
    }
}

```



### Test: `compact-0045-in.jsonld`

Expected:
```json
{
    "@context": {
        "term": "http://example.com/terms-are-not-considered-in-id",
        "compact-iris": "http://example.com/compact-iris-",
        "property": "http://example.com/property",
        "@vocab": "http://example.org/vocab-is-not-considered-for-id"
    },
    "@id": "term",
    "property": [
        {
            "@id": "compact-iris:are-considered",
            "property": "@id supports the following values: relative, absolute, and compact IRIs"
        },
        {
            "@id": "../parent-node",
            "property": "relative IRIs get resolved against the document's base IRI"
        }
    ]
}

```

Found output:
```json
{
    "@context": {
        "term": "http://example.com/terms-are-not-considered-in-id",
        "compact-iris": "http://example.com/compact-iris-",
        "property": "http://example.com/property",
        "@vocab": "http://example.org/vocab-is-not-considered-for-id"
    },
    "@id": "term",
    "property": [
        {
            "@id": "http://example.com/compact-iris-are-considered",
            "property": "@id supports the following values: relative, absolute, and compact IRIs"
        },
        {
            "@id": "../parent-node",
            "property": "relative IRIs get resolved against the document's base IRI"
        }
    ]
}

```


## Failures for error messages


### Test: `error-0042-in.jsonld`

Failed to raise error!


### Test: `error-0043-in.jsonld`

Failed to raise error!

## Failures for jsonld.expand

## Failures for jsonld.flatten

## Failures for jsonld.frame


### Test: `frame-0010-in.jsonld`

**RUNTIME ERROR!!**:  x$success isn't true. 


### Test: `frame-0020-in.jsonld`

**RUNTIME ERROR!!**:  x$success isn't true. 


### Test: `frame-0046-in.jsonld`

**RUNTIME ERROR!!**:  x$success isn't true. 


### Test: `frame-0047-in.jsonld`

**RUNTIME ERROR!!**:  x$success isn't true. 


### Test: `frame-0048-in.jsonld`

**RUNTIME ERROR!!**:  x$success isn't true. 


### Test: `frame-0049-in.jsonld`

**RUNTIME ERROR!!**:  x$success isn't true. 


### Test: `frame-0050-in.jsonld`

**RUNTIME ERROR!!**:  x$success isn't true. 

Official JSON-LD Test Suite
===========================

The JSON-LD Test Suite is a set of tests that can be used to verify
[specification conformance of JSON-LD Processors](http://json-ld.org/test-suite/reports/).
The test suite provides an easy and comprehensive testing solution for
developers implementing JSON-LD Processors.

More information about how to run the test suite or contribute new test
cases can be found at http://json-ld.org/test-suite/
