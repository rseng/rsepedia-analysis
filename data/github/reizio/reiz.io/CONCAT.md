<p align="center"><img src="https://github.com/reizio/reiz.io/blob/master/docs/logo.png"></p>

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03296/status.svg)](https://doi.org/10.21105/joss.03296)

# reiz.io

reiz.io is a structural source code search engine for Python. Compared to the
popular alternatives (e.g Github Code Search) it executes queries over the
syntax trees (instead of raw source code) and tries to retrive structural
knowledge (no semantics applied). For more information, please see the
[docs](https://reizio.readthedocs.io/en/latest/).

## A gentle introduction

Reiz is the code search framework that reiz.io is built a top on. Due to it's
nature, it solely works with the ASTs and intentionally avoids doing any
semantical work.

```{note}
Some ASTs attach a bit of contextual knowledge (e.g `Name(ctx=...)`
on python) which can be queried through simple matcher queries but
reiz.io doesn't include them when comparing references (see
matchers#reference-matcher for details).
```

Here is a simple ReizQL query that searches for a function that ends with a try
statement where we return a call to a function that has the same name as the
function we are within.

```python
FunctionDef(~func, body=[*..., Try(body=[Return(Call(Name(~func)))])])
```

which would match the following;

```py
def foo(spam):
    eggs = bar()
    try:
        return foo(spam + eggs)
    except ValueError:
        return None
```

In the very basic sense, it is generating the AST of the code above and checks
whether it fits the *pattern* (ReizQL query) or not.;

```py
FunctionDef(
    name='foo',
    args=arguments(
        posonlyargs=[],
        args=[arg(arg='spam', annotation=None, type_comment=None)],
        vararg=None,
        kwonlyargs=[],
        kw_defaults=[],
        kwarg=None,
        defaults=[],
    ),
    body=[
        Assign(
            targets=[Name(id='eggs', ctx=Store())],
            value=Call(
                func=Name(id='bar', ctx=Load()),
                args=[],
                keywords=[],
            ),
            type_comment=None,
        ),
        Try(
            body=[
                Return(
                    value=Call(
                        func=Name(id='foo', ctx=Load()),
                        args=[
                            BinOp(
                                left=Name(id='spam', ctx=Load()),
                                op=Add(),
                                right=Name(id='eggs', ctx=Load()),
                            ),
                        ],
                        keywords=[],
                    ),
                ),
            ],
            handlers=[
                ExceptHandler(
                    type=Name(id='ValueError', ctx=Load()),
                    name=None,
                    body=[
                        Return(
                            value=Constant(value=None, kind=None),
                        ),
                    ],
                ),
            ],
            orelse=[],
            finalbody=[],
        ),
    ],
    decorator_list=[],
    returns=None,
    type_comment=None,
)
```
# reiz.io

reiz.io is a structural source code search engine for Python. Compared to the
popular alternatives (e.g Github Code Search) it executes queries over the
syntax trees (instead of raw source code) and tries to retrive structural
knowledge (no semantics applied).

```{toctree}
:hidden:
:maxdepth: 1

reizql
internals
performance
contributing
installation-to-aws
related-projects
```

## Installation via Docker

A local instance of Reiz can be installed through `docker` and `docker-compose`
without any need to setup anythin else. It will run on a very small dataset (~75
files from 10 different projects) and will come with a bundled web interface.
Steps;

Get a fresh reiz clone

```
$ git clone https://github.com/reizio/reiz.io
```

Enter the directory and run `docker-compose up`

```
$ cd reiz.io
$ docker-compose up --build --remove-orphans
```

It will take about six to seven minutes for Reiz to initially build necessary
packages, install requirements, sample some packages for the dataset, prepare
the database and apply the schema, serialize those packages and finally run the
API.

```
reiz_1    | ... is inserted
reiz_1    | + python -m reiz.web.api
reiz_1    | [2021-04-24 21:35:38 +0000] [157] [INFO] Goin' Fast @ http://0.0.0.0:8000
reiz_1    | [2021-04-24 21:35:38,597] _helper         --- Goin' Fast @ http://0.0.0.0:8000
reiz_1    | [2021-04-24 21:35:38 +0000] [157] [INFO] Starting worker [157]
reiz_1    | [2021-04-24 21:35:38,929] serve           --- Starting worker [157]
```

After seeing the `Goin' Fast @ ...` message, you can open up your browser and
visit `http://localhost:8000/` and be greeted by the web UI.

## A gentle introduction

Reiz is the code search framework that reiz.io is built a top on. Due to it's
nature, it solely works with the ASTs and intentionally avoids doing any
semantical work.

```{note}
Some ASTs attach a bit of contextual knowledge (e.g `Name(ctx=...)`
on python) which can be queried through simple matcher queries but
reiz.io doesn't include them when comparing references (see
matchers#reference-matcher for details).
```

Here is a simple ReizQL query that searches for a function that ends with a try
statement where we return a call to a function that has the same name as the
function we are within.

```python
FunctionDef(~func, body=[*..., Try(body=[Return(Call(Name(~func)))])])
```

which would match the following;

```py
def foo(spam):
    eggs = bar()
    try:
        return foo(spam + eggs)
    except ValueError:
        return None
```

In the very basic sense, it is generating the AST of the code above and checks
whether it fits the *pattern* (ReizQL query) or not.;

```py
FunctionDef(
    name='foo',
    args=arguments(
        posonlyargs=[],
        args=[arg(arg='spam', annotation=None, type_comment=None)],
        vararg=None,
        kwonlyargs=[],
        kw_defaults=[],
        kwarg=None,
        defaults=[],
    ),
    body=[
        Assign(
            targets=[Name(id='eggs', ctx=Store())],
            value=Call(
                func=Name(id='bar', ctx=Load()),
                args=[],
                keywords=[],
            ),
            type_comment=None,
        ),
        Try(
            body=[
                Return(
                    value=Call(
                        func=Name(id='foo', ctx=Load()),
                        args=[
                            BinOp(
                                left=Name(id='spam', ctx=Load()),
                                op=Add(),
                                right=Name(id='eggs', ctx=Load()),
                            ),
                        ],
                        keywords=[],
                    ),
                ),
            ],
            handlers=[
                ExceptHandler(
                    type=Name(id='ValueError', ctx=Load()),
                    name=None,
                    body=[
                        Return(
                            value=Constant(value=None, kind=None),
                        ),
                    ],
                ),
            ],
            orelse=[],
            finalbody=[],
        ),
    ],
    decorator_list=[],
    returns=None,
    type_comment=None,
)
```

## QA

> What sort of source warehouses Reiz can index? Can I use my project on GitHub?

Short answer, currently it comes embedded with an indexer for the most common
python packages. And no, there is no formal support for individual package
deployment.

The source code sampling process for Reiz happens in 3 stages, which the initial
one is actually about indexing the target locations. Currently the embedded
tooling within reiz.io can index PyPI, through `reiz.sampling.get_packages`
module. The module's output is a JSON file that describes certain metadata
regarding the indexed repositories;

```json
app@frankfurt:~/reiz.io$ cat data/dataset.json | head
[
    {
        "name": "six",
        "downloads": 885510678,
        "git_source": "https://github.com/benjaminp/six",
        "git_revision": "3974f0c4f6700a5821b451abddff8b3ba6b2a04f",
        "license_type": "MIT"
    },
    ...
]
```

This list can be modified through hand, or you could simply write a script that
can generate the same sample for packages/repositories from other platforms. It
is also possible for you to just create a dataset of a single repository, that
being the project that you want to index. Though stating again, reiz.io as a
search engine for Python does not support you to add custom repositories in a
formal way.
# Contributing

Reiz is an open source project and would welcome any contributions under the
project code of conduct.

## Steps

- Create a fork of the `reizio/reiz.io` repository
- Clone your fork, add a new remote called `upstream` that points to the original
  one
- Sync your fork with the latest changes (`git fetch upstream && git rebase --hard upstream/master`)
- Find an issue to work on ([issues tagged with 'good first issue'](https://github.com/reizio/reiz.io/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22%5D)
- Create a new git branch (`git checkout -b <branch-name>`, the name can be either
  the `issue-<issue id>` or something that describes the change in a couple words)
- Write your PR, add new tests under `tests/` by following the pre-existing
  convention (one source file to `tests/dataset` and the correspondant query to
  the `tests/queries`)
- Run the tests to confirm everything workers (`python tests/runner.py`)
- Push your changes and create a PR

## Contributor Covenant Code of Conduct

### Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to make participation in our project and our
community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and
expression, level of experience, education, socio-economic status, nationality,
personal appearance, race, religion, or sexual identity and orientation.

### Our Standards

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

### Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, or to ban temporarily or permanently any
contributor for other behaviors that they deem inappropriate, threatening,
offensive, or harmful.

### Scope

This Code of Conduct applies within all project spaces, and it also applies when
an individual is representing the project or its community in public spaces.
Examples of representing a project or community include using an official
project e-mail address, posting via an official social media account, or acting
as an appointed representative at an online or offline event. Representation of
a project may be further defined and clarified by project maintainers.

### Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at
`isidentical [plus] reiz [at] gmail [dot] com`. All complaints will be reviewed
and investigated and will result in a response that is deemed necessary and
appropriate to the circumstances. The project team is obligated to maintain
confidentiality with regard to the reporter of an incident. Further details of
specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

### Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 1.4, available at
https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq

[homepage]: https://www.contributor-covenant.org
# ReizQL Language

ReizQL is a declarative query language for building AST matchers that work on
the Reiz platform.

:::{hint} Here is an example ReizQL query that searches for an if statement
where the body consists from a single assignment statement that assigns the
result of `requests.get(...)` call's result into a variable named `response`

```
if cache.invalidated:
    response = requests.get('https://api.reiz.io/refresh')
```

:::

```py
If(
    body = [
        Assign(
            targets = [
                Name('response')
            ],
            value = Call(
                Attribute(
                    Name('requests'),
                    'get'
                )
            )
        )
    ]
)
```

## Full Grammar

```bnf
start                   ::= match_pattern

pattern                 ::= negate_pattern
                         | or_pattern
                         | and_pattern
                         | match_pattern
                         | sequential_pattern
                         | reference_pattern
                         | match_string_pattern
                         | atom_pattern

negate_pattern          ::= "not" pattern
or_pattern              ::= pattern "|" pattern
and_pattern             ::= pattern "&" pattern
match_pattern           ::= NAME "(" ",".argument+ ")"
sequential_pattern      ::= "[" ",".(pattern | "*" IGNORE)+ "]"
reference_pattern       ::= "~" NAME
atom_pattern            ::= NONE
                         | STRING
                         | NUMBER
                         | IGNORE
                         | "f" STRING

argument                ::= pattern
                         | NAME "=" pattern

NONE                    ::= "None"
IGNORE                  ::= "..."
NAME                    ::= "a".."Z"
NUMBER                  ::= INTEGER | FLOAT
```

## Match Patterns

```bnf
match_pattern           ::= NAME "(" ",".argument+ ")"
```

Match patterns are the most fundamental part of the query expression. They
consist from an identifier (matcher name) which corresponds to an AST node type,
additionally they take any number of fields to be matched (values, optionally
attached with the corresponding field names).

All node types and fields are described in the
[Abstract Grammar](https://docs.python.org/3.8/library/ast.html#abstract-grammar)
of Python. Here are some entries from the ASDL;

```
module Python
{
    ...
    stmt = FunctionDef(identifier name, arguments args,
                       stmt* body, expr* decorator_list, expr? returns,
                       string? type_comment)
          | While(expr test, stmt* body, stmt* orelse)
          | If(expr test, stmt* body, stmt* orelse)
          | With(withitem* items, stmt* body, string? type_comment)

    expr = BoolOp(boolop op, expr* values)
         | NamedExpr(expr target, expr value)
         | BinOp(expr left, operator op, expr right)
         | UnaryOp(unaryop op, expr operand)
         | Lambda(arguments args, expr body)
         | IfExp(expr test, expr body, expr orelse)
         | Dict(expr* keys, expr* values)
```

The left hand side is the name of the base type, `stmt` would be a matcher that
could match all of the types in its right hand side (e.g `stmt()` would match
`FunctionDef()` / `While()` / `If()` / `With()`). Each element on the right hand
side are concrete matchers for that element in syntax. For example a `BinOp()`
represents a binary operation (2 operands), like `2 + 2` or `a % b()`.

Each element on the right hand side have different fields with types attached to
them. So the `BinOp()` node has 3 fields: `left`, `op`, `right` (respectively
they mean left hand side, operator, right hand side of an arithmetic operation).
`left` and the `right` must be another matcher from the `expr` base type (`BoolOp`
/ `NamedExpr`, ...). The star (`*`) at the end of type declaration implies that
it requires a [sequential pattern](#list-patterns) where the member types
inherit from that base type (e.g `stmt*` might be something like
`[If(), If(), While()]`). The question mark (`?`) indicates the value is
optional and can be `None`.

If the values are not named (e.g `BinOp(Constant())`) then the name will be
positionally given (`BinOp(Constant(), Add())` will be transformed to
`BinOp(left=Constant(), op=Add()`).

### Example Queries

- Match the `1994` literal

```py
Constant(1994)
```

- Match a binary operation where both sides are literals

```py
BinOp(left=Constant(), right=Constant())
```

- Match an (ternary) if expression that checks `a.b`'s truthness

```py
IfExp(
    test = Attribute(
        Name('a'),
        attr = 'b'
    )
)
```

## Sequential Patterns

```bnf
sequential_pattern      ::= "[" ",".(pattern | "*" IGNORE)+ "]"
```

Sequential patterns represent a list of subpatterns that are combined together
to match a list on the host AST. If we want to search a function definition
where there are 2 statements, the first one being an if statement and the second
one is a return of an identifier named `status` then we simply describe this
query like this;

```py
FunctionDef(
    body = [
        If(),
        Return(
            Name('status')
        )
    ]
)
```

Sequential patterns are ordered, and matched one-to-one unless a
[ignore star](#ignore-star) is seen.

### Ignore Star

If any of the elements on the sequence pattern is a star (`*`) followed by an
[ignore](#ignore-atom) then the matchers before the ignore-star are relative
from the beginning and the matchers after the ignore-star are relative to the
end of the sequence. This implies that there is no maximum limit of items (in
contrast to normal sequential patterns, where the number of elements is always
fixed to amount of patterns seen) and the minimum being the total amount of
matchers (excluding the ignore star).

Let's say we want to find a function that starts with an if statement, and then
ends with a call to `fetch` function.

```py
FunctionDef(
    body = [
        If(),
        *...,
        Return(
            Call(
                Name(
                    'fetch'
                )
            )
        )
    ]
)
```

There might be any number of elements between the if statement and the return,
and it simply won't care.

:::{note} If you need a filler value (for example you want the minimum number of
statements to be 3 instead of 2 in the case above) you can use
[ignore atom](#ignore-atom).

```py
FunctionDef(
    body = [
        If(),
        ...,
        *...,
        Return(
            Call(
                Name(
                    'fetch'
                )
            )
        )
    ]
)
```

:::

### Example Queries

- Match all functions that have 2 statements and the last being a return

```py
FunctionDef(
    body = [
        ...,
        Return()
    ]
)
```

- Match all try/except's where the last except handler is a bare `except: ...`

```py
Try(
    handlers = [
        *...,
        ExceptHandler(
            type = None
        )
    ]
)
```

## Logical Patterns

Logical patterns are different patterns connected together in the sense of some
logical operation (either `AND` or `OR`)

### AND patterns

```bnf
and_pattern             ::= pattern "&" pattern
```

`AND` patterns chains 2 different pattern together and matches the host value if
it can be matched by both of the connected patterns.

### OR patterns

```bnf
or_pattern              ::= pattern "|" pattern
```

`OR` patterns chains 2 different pattern together and matches the host value if
it can be matched by either of the connected patterns.

### Example Queries

- Match a return statement that either returns a list literal or a tuple literal

```py
Return(List() | Tuple())
```

- Match an if statement where the first statement is an assign and the total
  number of statements lower/equal than 5

```py
If(
    body = [
        Assign(),
        *...
    ] & LEN(max=5)
)
```

## NOT Patterns

```bnf
negate_pattern          ::= "not" pattern
```

For checking whether a certain pattern *does not* match on the host AST, the
negation operator can be used.

:::{hint} If a value is described as an optional (`?`) on the ASDL, then the
existence of value can be denoted via `not None` pattern.

### Example Queries

- Match a return statement that doesn't return a call

```py
Return(not Call())
```

- A list that doesn't start with tuples or sets

```py
List(
    elts = [
        (not Tuple()) & (not List()),
        *...
    ]
)
```

## Reference Patterns

```bnf
reference_pattern       ::= "~" NAME
```

Reference patterns are query-bound variables that can be referred elsewhere and
the truthness determined by checking whether all the references point to the
same expression (structurally, not semantically) or not.

### Example Queries

- Match a function definition where the last statement calls another function with
  the same name

```py
FunctionDef(
    ~name,
    body = [
        *...,
        Expr(Call(Name(~name)))
    ]
)
```

- Match an if statement where the test is a compare operation with the same
  lhs/rhs (`a == a` / `b() is b()`)

```py
If(
    test = CompareOp(
        left=~comp_expr,
        comparators = [
            ~comp_expr
        ]
    )
)
```

## Atom Patterns

```bnf
atom_pattern            ::= NONE
                         | STRING
                         | NUMBER
                         | IGNORE
                         | "f" STRING
```

Atoms represents basic values (like integers, strings) and also some
ReizQL-flavored constructs.

### Ignore

```
IGNORE                  ::= "..."
```

Ignore is a construct that just omits matching that field/element (in contrary
to None, where it means that value does not exist).

### None

```
NONE                    ::= "None"
```

None represents the absence of the value

### Match String

```
MATCH_STRING            ::= "f" STRING
```

| Pattern | Interpretation |
|----------------|------------------------------------| | `%` | matches zero or
more characters | | `_` | matches exactly one character | | `\%`/`\_` | matches
the literal `%`/`_` |

Match strings can match alike strings via denoting some basic structures (like
starts/ends with some static text).

### Example Queries

- Match a string that starts with `http://` or `https://`

```py
Constant(f'http://%' & f'https://%') 
```

- Match an arg that doesn't have any type annotations

```py
arg(annotation = None)
```

## Builtin Matchers

There are a couple of builtin matchers (builtin functions) that can match
against certain conditions.

### ALL/ANY

**Signature**: `ALL($0: pattern)` / `ANY($0: pattern)`

Apply the given matcher (`$0`) to a sequence. `ALL` would check whether all
elements can be matched through the given argument, and any would check if any
of the elements would be matched.

### LEN

**Signature**: `LEN($0: Opt[INTEGER], $1: Opt[INTEGER])`

Checks whether the length of the sequence fits into `$0 <= <host AST> <= $1`.
The `$0`/`$1` are optional values but at least one of them should be specified.

### META

**Signature**: `META(**metadata_providers)`

Checks for various metadata information (like file names, project names, parent
types, etc).

### I

**Signature**: `I($0: atom_pattern[MATCH_STRING])`

Supports case insensitive match through match strings.

### Example Queries

- Match a tuple where all members are literals

```py
Tuple(ALL(Constant()))
```

- Match a function where one of the top level statements is an if statement

```py
FunctionDef(
    body = ANY(
        If()
    )
)
```

- Match a function call where there are minimum 3 positional arguments and 5
  maximum keyword arguments

```py
Call(
    args = LEN(min=3),
    keywords = LEN(max=5)
)
```

- Match a string in an case insensitive way

```py
Constant(I(f"foo"))
```
# Demo mode installation to an AWS machine

## Requirements

- AWS account
- SSH client
- Browser (firefox/chrome)

## Instance

- Image: `Ubuntu Server 18.04 LTS (HVM), SSD Volume Type`
- Type: `t2.large` (2 vCPU, 8GB ram)
- Storage: 30 GB

## Install Docker / Docker-Compose

Find the instance IP address and then login the instance via SSH (with
forwarding the `8000` port);

```
$ ssh -L 8000:localhost:8000 ubuntu@<ip addr> -i <identity file>
Welcome to Ubuntu 18.04.5 LTS (GNU/Linux 5.4.0-1045-aws x86_64)

  System information as of Tue Jun  8 17:14:49 UTC 2021

  System load:  0.23              Processes:           117
  Usage of /:   3.9% of 29.02GB   Users logged in:     0
  Memory usage: 2%                IP address for eth0: 172.31.9.52
  Swap usage:   0%

ubuntu@ip-172-31-9-52:~$ 
```

After the login, go and install docker.

```
$ sudo apt update
$ sudo apt install apt-transport-https ca-certificates curl software-properties-common
$ curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
$ sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu bionic stable"
$ sudo apt update
$ sudo apt install docker-ce
```

Add your user to the `docker` group to not to use `sudo` for every command

```
$ sudo groupadd docker
$ sudo usermod -aG docker $USER
$ newgrp docker
$ docker --version
Docker version 20.10.7, build f0df350
```

Ensure that you use either the same or a newer version. Then install
docker-compose

```
$ sudo curl -L https://github.com/docker/compose/releases/download/1.29.2/docker-compose-`uname -s`-`uname -m` -o /usr/local/bin/docker-compose
$ sudo chmod +x /usr/local/bin/docker-compose
$ docker-compose --version
docker-compose version 1.29.2, build 5becea4c
```

## Install Reiz

For installing reiz, you only need `git` and `docker-compose`.

Clone and ensure that you have the latest revision;

```
$ git clone https://github.com/reizio/reiz.io
$ cd reiz.io/
~/reiz.io$ git fetch origin
~/reiz.io$ git reset --hard origin/master
```

And finally start the reiz

```
$ docker-compose up
```

In case of you are interested, here are the
[full logs](https://gist.github.com/isidentical/bf6b4e2dbdd60407a4b51d6fbfc8e28a).
You should pretty much get the similiar stuff. After getting these lines;

```
reiz_1    | [2021-06-08 17:29:34,186] insert_file     --- 'pip/tests/lib/certs.py' has been inserted successfully
reiz_1    | [2021-06-08 17:29:37,175] insert_file     --- 'pip/tests/lib/local_repos.py' has been inserted successfully
reiz_1    | [2021-06-08 17:29:40,128] insert_file     --- 'pip/tests/lib/configuration_helpers.py' has been inserted successfully
reiz_1    | + python -m reiz.web.api
reiz_1    | [2021-06-08 17:29:40 +0000] [151] [INFO] Goin' Fast @ http://0.0.0.0:8000
reiz_1    | [2021-06-08 17:29:40,944] _helper         --- Goin' Fast @ http://0.0.0.0:8000
reiz_1    | [2021-06-08 17:29:41 +0000] [151] [INFO] Starting worker [151]
reiz_1    | [2021-06-08 17:29:41,226] serve           --- Starting worker [151]
```

You can go to the `localhost:8000` on your browser and be greeted by the web
page;

![image](https://user-images.githubusercontent.com/47358913/121231015-911b4600-c898-11eb-9c99-5d46efb4d356.png)

After that you could either click on one of pre-selected queries or write your
own;

![image](https://user-images.githubusercontent.com/47358913/121231095-a98b6080-c898-11eb-9fcd-d250ae44ad5f.png)

If you want to see the size of dataset, you could go to the
`localhost:8000/stats`. Mine for example indexes 10k nodes;

```
{
    "status": "success",
    "results": {
        "Module": 91,
        "AST": 10785,
        "stmt": 2204,
        "expr": 8426
    },
    "exception": null
}
```

Every time you do `docker-compose up`, it will index more files from the 10
projects it downloaded (~75, +/- 20).
# Reiz Internals

## Flow

![flowgraph](flowgraph.png)

## Components

```
Warehouse Preparation (reiz.schema):
    -> reiz.schema.builders.$DB
        Build a schema from the given ASDL file
    -> reiz.schema.$DB
        Store schema metadata for the specific database

Data Collection (reiz.sampling):
    -> reiz.sampling.get_dataset
        Get a list of possible projects with cross references to their SCM pages
    -> reiz.sampling.fetch_dataset
        Download the given list of projects
    -> reiz.sampling.sanitize_dataset
        Remove everything beside valid Python 3 source code files

Data Serialization (reiz.serialization):
    -> reiz.serialization.transformers
        Transform and annotate the raw language AST for querying
    -> reiz.serialization.serializer
        Serialize all source files in a single project to the database
    -> reiz.serialization.serialize
        Serialize all downloaded projects to the database

Data Querying [ReizQL] (reiz.reizql):
    -> reiz.reizql.parser
        -> reiz.reizql.parser.grammar
            Represent Reiz AST
        -> reiz.reizql.parser.parse
            Generate Reiz AST from ReizQL
    -> reiz.reizql.compiler
        Generate IR from Reiz AST
```
# Performance Evaluation

| query                                            | timing (s) |
| ------------------------------------------------ | ---------- |
| `Call(Name("len"))`                              | 0.025985   |
| `BinOp(op=Add() \| Sub())`                       | 0.030508   |
| `Try(handlers=LEN(min=3, max=5))`                | 0.033486   |
| `BinOp(left=Constant(), right=Constant())`       | 0.146516   |
| `FunctionDef(f"run_%", returns = not None)`      | 0.0216     |
| `ClassDef(body=[Assign(), *..., FunctionDef()])` | 0.28737    |

## Analysis

There are 2 major points that cost nearly %95 of the whole query operation. The first,
and the obvious point is the actually running the query in the database. There are a couple
points that Reizc an do to optimize this step, including trying to generate the best
possible query while being in a linear motion (for supporting constructs like reference
variables). The code generator (`reiz.reizql.compiler`) went through a couple major
refactors for performance reasons (e.g [#12](https://github.com/reizio/reiz.io/pull/12)).
Also there is a simple/naive [AST optimization pass](https://github.com/reizio/reiz.io/blob/cff3cc6eaad532ac1a956c1f7c7a58d97ea00e4b/reiz/ir/backends/edgeql.py#L461-L513) on
the IR (EdgeQL) itself.

The second part is the actually retrieving the code snippets from the disk itself. We
already store a lot of metadata (like start/end positions, github project etc.) but
the actual 'source' is still on the disk. So after retrieving the filenames from the
query, we simply go and read those files and get the related segments. This is an area
that is open to more optimizations (we could statically determine the byte-range and
only fetch it, we could parallelize this for multiple matches \[the default resultset
come with 10 matches\], ...), though these won't have the same effects as in getting
a better speed in the DB.

Of course alongside these, there have been tons of ways to optimize postgresql itself
for different workloads, though it is outside of the Reiz project.

## Setup

Machine;

|              |                        |
| ------------ | ---------------------- |
| provider     | digital ocean          |
| service type | droplet (basic plan)   |
| cpu          | (shared) 2vCPU         |
| ram          | 2GB                    |
| disk         | regular SSD (not NVME) |

IndexDB;

|                 |            |
| --------------- | ---------- |
| total files     | 53k        |
| total AST nodes | 17 521 894 |

Benchmark script is present at the source checkout (`scripts/benchmark_doc.py`).
# Related Technologies / Projects

Reiz is preceded by many different projects, and there are also some others
which might fit to some scopes better than Reiz depending on your use case. Here
are some of the related projects;

- [irun](https://github.com/reizio/irun)
- [Sourcer](https://www.sciencedirect.com/science/article/pii/S016764231200072X)
- [GitHub Code Search](github.com/search)
- [Debian Code Search](http://codesearch.debian.net/research/bsc-thesis.pdf)
- [grep.app](https://grep.app)
- [krugle](krugle.com)
- [astpath](https://github.com/hchasestevens/astpath)
- [astsearch](https://github.com/takluyver/astsearch)
